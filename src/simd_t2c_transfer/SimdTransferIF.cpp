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


#include "SimdTransferIF.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_hrr_if(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t id, const size_t kd, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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
    const auto *id_109 = buffer.data(id + 109);
    const auto *id_110 = buffer.data(id + 110);
    const auto *id_111 = buffer.data(id + 111);
    const auto *id_112 = buffer.data(id + 112);
    const auto *id_113 = buffer.data(id + 113);
    const auto *id_114 = buffer.data(id + 114);
    const auto *id_115 = buffer.data(id + 115);
    const auto *id_116 = buffer.data(id + 116);
    const auto *id_117 = buffer.data(id + 117);
    const auto *id_118 = buffer.data(id + 118);
    const auto *id_119 = buffer.data(id + 119);
    const auto *id_120 = buffer.data(id + 120);
    const auto *id_121 = buffer.data(id + 121);
    const auto *id_122 = buffer.data(id + 122);
    const auto *id_123 = buffer.data(id + 123);
    const auto *id_124 = buffer.data(id + 124);
    const auto *id_125 = buffer.data(id + 125);
    const auto *id_126 = buffer.data(id + 126);
    const auto *id_127 = buffer.data(id + 127);
    const auto *id_128 = buffer.data(id + 128);
    const auto *id_129 = buffer.data(id + 129);
    const auto *id_130 = buffer.data(id + 130);
    const auto *id_131 = buffer.data(id + 131);
    const auto *id_132 = buffer.data(id + 132);
    const auto *id_133 = buffer.data(id + 133);
    const auto *id_134 = buffer.data(id + 134);
    const auto *id_135 = buffer.data(id + 135);
    const auto *id_136 = buffer.data(id + 136);
    const auto *id_137 = buffer.data(id + 137);
    const auto *id_138 = buffer.data(id + 138);
    const auto *id_139 = buffer.data(id + 139);
    const auto *id_140 = buffer.data(id + 140);
    const auto *id_141 = buffer.data(id + 141);
    const auto *id_142 = buffer.data(id + 142);
    const auto *id_143 = buffer.data(id + 143);
    const auto *id_144 = buffer.data(id + 144);
    const auto *id_145 = buffer.data(id + 145);
    const auto *id_146 = buffer.data(id + 146);
    const auto *id_147 = buffer.data(id + 147);
    const auto *id_148 = buffer.data(id + 148);
    const auto *id_149 = buffer.data(id + 149);
    const auto *id_150 = buffer.data(id + 150);
    const auto *id_151 = buffer.data(id + 151);
    const auto *id_152 = buffer.data(id + 152);
    const auto *id_153 = buffer.data(id + 153);
    const auto *id_154 = buffer.data(id + 154);
    const auto *id_155 = buffer.data(id + 155);
    const auto *id_156 = buffer.data(id + 156);
    const auto *id_157 = buffer.data(id + 157);
    const auto *id_158 = buffer.data(id + 158);
    const auto *id_159 = buffer.data(id + 159);
    const auto *id_160 = buffer.data(id + 160);
    const auto *id_161 = buffer.data(id + 161);
    const auto *id_162 = buffer.data(id + 162);
    const auto *id_163 = buffer.data(id + 163);
    const auto *id_164 = buffer.data(id + 164);
    const auto *id_165 = buffer.data(id + 165);
    const auto *id_166 = buffer.data(id + 166);
    const auto *id_167 = buffer.data(id + 167);

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
    const auto *kd_93 = buffer.data(kd + 93);
    const auto *kd_94 = buffer.data(kd + 94);
    const auto *kd_95 = buffer.data(kd + 95);
    const auto *kd_96 = buffer.data(kd + 96);
    const auto *kd_97 = buffer.data(kd + 97);
    const auto *kd_98 = buffer.data(kd + 98);
    const auto *kd_99 = buffer.data(kd + 99);
    const auto *kd_100 = buffer.data(kd + 100);
    const auto *kd_101 = buffer.data(kd + 101);
    const auto *kd_102 = buffer.data(kd + 102);
    const auto *kd_103 = buffer.data(kd + 103);
    const auto *kd_104 = buffer.data(kd + 104);
    const auto *kd_105 = buffer.data(kd + 105);
    const auto *kd_106 = buffer.data(kd + 106);
    const auto *kd_107 = buffer.data(kd + 107);
    const auto *kd_108 = buffer.data(kd + 108);
    const auto *kd_109 = buffer.data(kd + 109);
    const auto *kd_110 = buffer.data(kd + 110);
    const auto *kd_111 = buffer.data(kd + 111);
    const auto *kd_112 = buffer.data(kd + 112);
    const auto *kd_113 = buffer.data(kd + 113);
    const auto *kd_114 = buffer.data(kd + 114);
    const auto *kd_115 = buffer.data(kd + 115);
    const auto *kd_116 = buffer.data(kd + 116);
    const auto *kd_117 = buffer.data(kd + 117);
    const auto *kd_118 = buffer.data(kd + 118);
    const auto *kd_119 = buffer.data(kd + 119);
    const auto *kd_120 = buffer.data(kd + 120);
    const auto *kd_121 = buffer.data(kd + 121);
    const auto *kd_122 = buffer.data(kd + 122);
    const auto *kd_123 = buffer.data(kd + 123);
    const auto *kd_124 = buffer.data(kd + 124);
    const auto *kd_125 = buffer.data(kd + 125);
    const auto *kd_126 = buffer.data(kd + 126);
    const auto *kd_127 = buffer.data(kd + 127);
    const auto *kd_128 = buffer.data(kd + 128);
    const auto *kd_129 = buffer.data(kd + 129);
    const auto *kd_130 = buffer.data(kd + 130);
    const auto *kd_131 = buffer.data(kd + 131);
    const auto *kd_132 = buffer.data(kd + 132);
    const auto *kd_133 = buffer.data(kd + 133);
    const auto *kd_134 = buffer.data(kd + 134);
    const auto *kd_135 = buffer.data(kd + 135);
    const auto *kd_136 = buffer.data(kd + 136);
    const auto *kd_137 = buffer.data(kd + 137);
    const auto *kd_138 = buffer.data(kd + 138);
    const auto *kd_139 = buffer.data(kd + 139);
    const auto *kd_140 = buffer.data(kd + 140);
    const auto *kd_141 = buffer.data(kd + 141);
    const auto *kd_142 = buffer.data(kd + 142);
    const auto *kd_143 = buffer.data(kd + 143);
    const auto *kd_144 = buffer.data(kd + 144);
    const auto *kd_145 = buffer.data(kd + 145);
    const auto *kd_146 = buffer.data(kd + 146);
    const auto *kd_147 = buffer.data(kd + 147);
    const auto *kd_148 = buffer.data(kd + 148);
    const auto *kd_149 = buffer.data(kd + 149);
    const auto *kd_150 = buffer.data(kd + 150);
    const auto *kd_151 = buffer.data(kd + 151);
    const auto *kd_152 = buffer.data(kd + 152);
    const auto *kd_153 = buffer.data(kd + 153);
    const auto *kd_154 = buffer.data(kd + 154);
    const auto *kd_155 = buffer.data(kd + 155);
    const auto *kd_156 = buffer.data(kd + 156);
    const auto *kd_157 = buffer.data(kd + 157);
    const auto *kd_158 = buffer.data(kd + 158);
    const auto *kd_159 = buffer.data(kd + 159);
    const auto *kd_160 = buffer.data(kd + 160);
    const auto *kd_161 = buffer.data(kd + 161);
    const auto *kd_162 = buffer.data(kd + 162);
    const auto *kd_163 = buffer.data(kd + 163);
    const auto *kd_164 = buffer.data(kd + 164);
    const auto *kd_165 = buffer.data(kd + 165);
    const auto *kd_166 = buffer.data(kd + 166);
    const auto *kd_167 = buffer.data(kd + 167);
    const auto *kd_171 = buffer.data(kd + 171);
    const auto *kd_172 = buffer.data(kd + 172);
    const auto *kd_173 = buffer.data(kd + 173);
    const auto *kd_177 = buffer.data(kd + 177);
    const auto *kd_178 = buffer.data(kd + 178);
    const auto *kd_179 = buffer.data(kd + 179);
    const auto *kd_183 = buffer.data(kd + 183);
    const auto *kd_184 = buffer.data(kd + 184);
    const auto *kd_185 = buffer.data(kd + 185);
    const auto *kd_189 = buffer.data(kd + 189);
    const auto *kd_190 = buffer.data(kd + 190);
    const auto *kd_191 = buffer.data(kd + 191);
    const auto *kd_195 = buffer.data(kd + 195);
    const auto *kd_196 = buffer.data(kd + 196);
    const auto *kd_197 = buffer.data(kd + 197);
    const auto *kd_201 = buffer.data(kd + 201);
    const auto *kd_202 = buffer.data(kd + 202);
    const auto *kd_203 = buffer.data(kd + 203);
    const auto *kd_207 = buffer.data(kd + 207);
    const auto *kd_208 = buffer.data(kd + 208);
    const auto *kd_209 = buffer.data(kd + 209);
    const auto *kd_215 = buffer.data(kd + 215);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, id_0, id_1, id_2, id_3, id_4, kd_0, \
                         kd_1, kd_2, kd_3, kd_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * id_0[k]
                 + kd_0[k];

        t_1[k] = ab_x[k] * id_1[k]
                 + kd_1[k];

        t_2[k] = ab_x[k] * id_2[k]
                 + kd_2[k];

        t_3[k] = ab_x[k] * id_3[k]
                 + kd_3[k];

        t_4[k] = ab_x[k] * id_4[k]
                 + kd_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ab_y, ab_z, id_3, id_4, id_5, kd_5, \
                         kd_9, kd_10, kd_11, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = ab_x[k] * id_5[k]
                 + kd_5[k];

        t_6[k] = ab_y[k] * id_3[k]
                 + kd_9[k];

        t_7[k] = ab_y[k] * id_4[k]
                 + kd_10[k];

        t_8[k] = ab_y[k] * id_5[k]
                 + kd_11[k];

        t_9[k] = ab_z[k] * id_5[k]
                 + kd_17[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, id_6, id_7, id_8, id_9, id_10, \
                         kd_6, kd_7, kd_8, kd_9, kd_10 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = ab_x[k] * id_6[k]
                  + kd_6[k];

        t_11[k] = ab_x[k] * id_7[k]
                  + kd_7[k];

        t_12[k] = ab_x[k] * id_8[k]
                  + kd_8[k];

        t_13[k] = ab_x[k] * id_9[k]
                  + kd_9[k];

        t_14[k] = ab_x[k] * id_10[k]
                  + kd_10[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, ab_y, ab_z, id_9, id_10, id_11, \
                         kd_11, kd_21, kd_22, kd_23, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = ab_x[k] * id_11[k]
                  + kd_11[k];

        t_16[k] = ab_y[k] * id_9[k]
                  + kd_21[k];

        t_17[k] = ab_y[k] * id_10[k]
                  + kd_22[k];

        t_18[k] = ab_y[k] * id_11[k]
                  + kd_23[k];

        t_19[k] = ab_z[k] * id_11[k]
                  + kd_29[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, id_12, id_13, id_14, id_15, \
                         id_16, kd_12, kd_13, kd_14, kd_15, kd_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = ab_x[k] * id_12[k]
                  + kd_12[k];

        t_21[k] = ab_x[k] * id_13[k]
                  + kd_13[k];

        t_22[k] = ab_x[k] * id_14[k]
                  + kd_14[k];

        t_23[k] = ab_x[k] * id_15[k]
                  + kd_15[k];

        t_24[k] = ab_x[k] * id_16[k]
                  + kd_16[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, ab_y, ab_z, id_15, id_16, id_17, \
                         kd_17, kd_27, kd_28, kd_29, kd_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = ab_x[k] * id_17[k]
                  + kd_17[k];

        t_26[k] = ab_y[k] * id_15[k]
                  + kd_27[k];

        t_27[k] = ab_y[k] * id_16[k]
                  + kd_28[k];

        t_28[k] = ab_y[k] * id_17[k]
                  + kd_29[k];

        t_29[k] = ab_z[k] * id_17[k]
                  + kd_35[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, id_18, id_19, id_20, id_21, \
                         id_22, kd_18, kd_19, kd_20, kd_21, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = ab_x[k] * id_18[k]
                  + kd_18[k];

        t_31[k] = ab_x[k] * id_19[k]
                  + kd_19[k];

        t_32[k] = ab_x[k] * id_20[k]
                  + kd_20[k];

        t_33[k] = ab_x[k] * id_21[k]
                  + kd_21[k];

        t_34[k] = ab_x[k] * id_22[k]
                  + kd_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, ab_y, ab_z, id_21, id_22, id_23, \
                         kd_23, kd_39, kd_40, kd_41, kd_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = ab_x[k] * id_23[k]
                  + kd_23[k];

        t_36[k] = ab_y[k] * id_21[k]
                  + kd_39[k];

        t_37[k] = ab_y[k] * id_22[k]
                  + kd_40[k];

        t_38[k] = ab_y[k] * id_23[k]
                  + kd_41[k];

        t_39[k] = ab_z[k] * id_23[k]
                  + kd_47[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, id_24, id_25, id_26, id_27, \
                         id_28, kd_24, kd_25, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = ab_x[k] * id_24[k]
                  + kd_24[k];

        t_41[k] = ab_x[k] * id_25[k]
                  + kd_25[k];

        t_42[k] = ab_x[k] * id_26[k]
                  + kd_26[k];

        t_43[k] = ab_x[k] * id_27[k]
                  + kd_27[k];

        t_44[k] = ab_x[k] * id_28[k]
                  + kd_28[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, ab_y, ab_z, id_27, id_28, id_29, \
                         kd_29, kd_45, kd_46, kd_47, kd_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = ab_x[k] * id_29[k]
                  + kd_29[k];

        t_46[k] = ab_y[k] * id_27[k]
                  + kd_45[k];

        t_47[k] = ab_y[k] * id_28[k]
                  + kd_46[k];

        t_48[k] = ab_y[k] * id_29[k]
                  + kd_47[k];

        t_49[k] = ab_z[k] * id_29[k]
                  + kd_53[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, id_30, id_31, id_32, id_33, \
                         id_34, kd_30, kd_31, kd_32, kd_33, kd_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = ab_x[k] * id_30[k]
                  + kd_30[k];

        t_51[k] = ab_x[k] * id_31[k]
                  + kd_31[k];

        t_52[k] = ab_x[k] * id_32[k]
                  + kd_32[k];

        t_53[k] = ab_x[k] * id_33[k]
                  + kd_33[k];

        t_54[k] = ab_x[k] * id_34[k]
                  + kd_34[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, id_33, id_34, id_35, \
                         kd_35, kd_51, kd_52, kd_53, kd_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = ab_x[k] * id_35[k]
                  + kd_35[k];

        t_56[k] = ab_y[k] * id_33[k]
                  + kd_51[k];

        t_57[k] = ab_y[k] * id_34[k]
                  + kd_52[k];

        t_58[k] = ab_y[k] * id_35[k]
                  + kd_53[k];

        t_59[k] = ab_z[k] * id_35[k]
                  + kd_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, id_36, id_37, id_38, id_39, \
                         id_40, kd_36, kd_37, kd_38, kd_39, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = ab_x[k] * id_36[k]
                  + kd_36[k];

        t_61[k] = ab_x[k] * id_37[k]
                  + kd_37[k];

        t_62[k] = ab_x[k] * id_38[k]
                  + kd_38[k];

        t_63[k] = ab_x[k] * id_39[k]
                  + kd_39[k];

        t_64[k] = ab_x[k] * id_40[k]
                  + kd_40[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, ab_y, ab_z, id_39, id_40, id_41, \
                         kd_41, kd_63, kd_64, kd_65, kd_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = ab_x[k] * id_41[k]
                  + kd_41[k];

        t_66[k] = ab_y[k] * id_39[k]
                  + kd_63[k];

        t_67[k] = ab_y[k] * id_40[k]
                  + kd_64[k];

        t_68[k] = ab_y[k] * id_41[k]
                  + kd_65[k];

        t_69[k] = ab_z[k] * id_41[k]
                  + kd_71[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, id_42, id_43, id_44, id_45, \
                         id_46, kd_42, kd_43, kd_44, kd_45, kd_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = ab_x[k] * id_42[k]
                  + kd_42[k];

        t_71[k] = ab_x[k] * id_43[k]
                  + kd_43[k];

        t_72[k] = ab_x[k] * id_44[k]
                  + kd_44[k];

        t_73[k] = ab_x[k] * id_45[k]
                  + kd_45[k];

        t_74[k] = ab_x[k] * id_46[k]
                  + kd_46[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, ab_y, ab_z, id_45, id_46, id_47, \
                         kd_47, kd_69, kd_70, kd_71, kd_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = ab_x[k] * id_47[k]
                  + kd_47[k];

        t_76[k] = ab_y[k] * id_45[k]
                  + kd_69[k];

        t_77[k] = ab_y[k] * id_46[k]
                  + kd_70[k];

        t_78[k] = ab_y[k] * id_47[k]
                  + kd_71[k];

        t_79[k] = ab_z[k] * id_47[k]
                  + kd_77[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, id_48, id_49, id_50, id_51, \
                         id_52, kd_48, kd_49, kd_50, kd_51, kd_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = ab_x[k] * id_48[k]
                  + kd_48[k];

        t_81[k] = ab_x[k] * id_49[k]
                  + kd_49[k];

        t_82[k] = ab_x[k] * id_50[k]
                  + kd_50[k];

        t_83[k] = ab_x[k] * id_51[k]
                  + kd_51[k];

        t_84[k] = ab_x[k] * id_52[k]
                  + kd_52[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, id_51, id_52, id_53, \
                         kd_53, kd_75, kd_76, kd_77, kd_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = ab_x[k] * id_53[k]
                  + kd_53[k];

        t_86[k] = ab_y[k] * id_51[k]
                  + kd_75[k];

        t_87[k] = ab_y[k] * id_52[k]
                  + kd_76[k];

        t_88[k] = ab_y[k] * id_53[k]
                  + kd_77[k];

        t_89[k] = ab_z[k] * id_53[k]
                  + kd_83[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, id_54, id_55, id_56, id_57, \
                         id_58, kd_54, kd_55, kd_56, kd_57, kd_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = ab_x[k] * id_54[k]
                  + kd_54[k];

        t_91[k] = ab_x[k] * id_55[k]
                  + kd_55[k];

        t_92[k] = ab_x[k] * id_56[k]
                  + kd_56[k];

        t_93[k] = ab_x[k] * id_57[k]
                  + kd_57[k];

        t_94[k] = ab_x[k] * id_58[k]
                  + kd_58[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, ab_y, ab_z, id_57, id_58, id_59, \
                         kd_59, kd_81, kd_82, kd_83, kd_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = ab_x[k] * id_59[k]
                  + kd_59[k];

        t_96[k] = ab_y[k] * id_57[k]
                  + kd_81[k];

        t_97[k] = ab_y[k] * id_58[k]
                  + kd_82[k];

        t_98[k] = ab_y[k] * id_59[k]
                  + kd_83[k];

        t_99[k] = ab_z[k] * id_59[k]
                  + kd_89[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, id_60, id_61, id_62, id_63, \
                         id_64, kd_60, kd_61, kd_62, kd_63, kd_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_100[k] = ab_x[k] * id_60[k]
                   + kd_60[k];

        t_101[k] = ab_x[k] * id_61[k]
                   + kd_61[k];

        t_102[k] = ab_x[k] * id_62[k]
                   + kd_62[k];

        t_103[k] = ab_x[k] * id_63[k]
                   + kd_63[k];

        t_104[k] = ab_x[k] * id_64[k]
                   + kd_64[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, ab_y, ab_z, id_63, id_64, \
                         id_65, kd_65, kd_93, kd_94, kd_95, kd_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_105[k] = ab_x[k] * id_65[k]
                   + kd_65[k];

        t_106[k] = ab_y[k] * id_63[k]
                   + kd_93[k];

        t_107[k] = ab_y[k] * id_64[k]
                   + kd_94[k];

        t_108[k] = ab_y[k] * id_65[k]
                   + kd_95[k];

        t_109[k] = ab_z[k] * id_65[k]
                   + kd_101[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, id_66, id_67, id_68, id_69, \
                         id_70, kd_66, kd_67, kd_68, kd_69, kd_70 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_110[k] = ab_x[k] * id_66[k]
                   + kd_66[k];

        t_111[k] = ab_x[k] * id_67[k]
                   + kd_67[k];

        t_112[k] = ab_x[k] * id_68[k]
                   + kd_68[k];

        t_113[k] = ab_x[k] * id_69[k]
                   + kd_69[k];

        t_114[k] = ab_x[k] * id_70[k]
                   + kd_70[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, ab_y, ab_z, id_69, id_70, \
                         id_71, kd_71, kd_99, kd_100, kd_101, kd_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_115[k] = ab_x[k] * id_71[k]
                   + kd_71[k];

        t_116[k] = ab_y[k] * id_69[k]
                   + kd_99[k];

        t_117[k] = ab_y[k] * id_70[k]
                   + kd_100[k];

        t_118[k] = ab_y[k] * id_71[k]
                   + kd_101[k];

        t_119[k] = ab_z[k] * id_71[k]
                   + kd_107[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, id_72, id_73, id_74, id_75, \
                         id_76, kd_72, kd_73, kd_74, kd_75, kd_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_120[k] = ab_x[k] * id_72[k]
                   + kd_72[k];

        t_121[k] = ab_x[k] * id_73[k]
                   + kd_73[k];

        t_122[k] = ab_x[k] * id_74[k]
                   + kd_74[k];

        t_123[k] = ab_x[k] * id_75[k]
                   + kd_75[k];

        t_124[k] = ab_x[k] * id_76[k]
                   + kd_76[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, ab_y, ab_z, id_75, id_76, \
                         id_77, kd_77, kd_105, kd_106, kd_107, kd_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_125[k] = ab_x[k] * id_77[k]
                   + kd_77[k];

        t_126[k] = ab_y[k] * id_75[k]
                   + kd_105[k];

        t_127[k] = ab_y[k] * id_76[k]
                   + kd_106[k];

        t_128[k] = ab_y[k] * id_77[k]
                   + kd_107[k];

        t_129[k] = ab_z[k] * id_77[k]
                   + kd_113[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, id_78, id_79, id_80, id_81, \
                         id_82, kd_78, kd_79, kd_80, kd_81, kd_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_130[k] = ab_x[k] * id_78[k]
                   + kd_78[k];

        t_131[k] = ab_x[k] * id_79[k]
                   + kd_79[k];

        t_132[k] = ab_x[k] * id_80[k]
                   + kd_80[k];

        t_133[k] = ab_x[k] * id_81[k]
                   + kd_81[k];

        t_134[k] = ab_x[k] * id_82[k]
                   + kd_82[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, ab_y, ab_z, id_81, id_82, \
                         id_83, kd_83, kd_111, kd_112, kd_113, kd_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = ab_x[k] * id_83[k]
                   + kd_83[k];

        t_136[k] = ab_y[k] * id_81[k]
                   + kd_111[k];

        t_137[k] = ab_y[k] * id_82[k]
                   + kd_112[k];

        t_138[k] = ab_y[k] * id_83[k]
                   + kd_113[k];

        t_139[k] = ab_z[k] * id_83[k]
                   + kd_119[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, id_84, id_85, id_86, id_87, \
                         id_88, kd_84, kd_85, kd_86, kd_87, kd_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_140[k] = ab_x[k] * id_84[k]
                   + kd_84[k];

        t_141[k] = ab_x[k] * id_85[k]
                   + kd_85[k];

        t_142[k] = ab_x[k] * id_86[k]
                   + kd_86[k];

        t_143[k] = ab_x[k] * id_87[k]
                   + kd_87[k];

        t_144[k] = ab_x[k] * id_88[k]
                   + kd_88[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, ab_y, ab_z, id_87, id_88, \
                         id_89, kd_89, kd_117, kd_118, kd_119, kd_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_145[k] = ab_x[k] * id_89[k]
                   + kd_89[k];

        t_146[k] = ab_y[k] * id_87[k]
                   + kd_117[k];

        t_147[k] = ab_y[k] * id_88[k]
                   + kd_118[k];

        t_148[k] = ab_y[k] * id_89[k]
                   + kd_119[k];

        t_149[k] = ab_z[k] * id_89[k]
                   + kd_125[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, id_90, id_91, id_92, id_93, \
                         id_94, kd_90, kd_91, kd_92, kd_93, kd_94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_150[k] = ab_x[k] * id_90[k]
                   + kd_90[k];

        t_151[k] = ab_x[k] * id_91[k]
                   + kd_91[k];

        t_152[k] = ab_x[k] * id_92[k]
                   + kd_92[k];

        t_153[k] = ab_x[k] * id_93[k]
                   + kd_93[k];

        t_154[k] = ab_x[k] * id_94[k]
                   + kd_94[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, ab_y, ab_z, id_93, id_94, \
                         id_95, kd_95, kd_129, kd_130, kd_131, kd_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_155[k] = ab_x[k] * id_95[k]
                   + kd_95[k];

        t_156[k] = ab_y[k] * id_93[k]
                   + kd_129[k];

        t_157[k] = ab_y[k] * id_94[k]
                   + kd_130[k];

        t_158[k] = ab_y[k] * id_95[k]
                   + kd_131[k];

        t_159[k] = ab_z[k] * id_95[k]
                   + kd_137[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, id_96, id_97, id_98, id_99, \
                         id_100, kd_96, kd_97, kd_98, kd_99, kd_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_160[k] = ab_x[k] * id_96[k]
                   + kd_96[k];

        t_161[k] = ab_x[k] * id_97[k]
                   + kd_97[k];

        t_162[k] = ab_x[k] * id_98[k]
                   + kd_98[k];

        t_163[k] = ab_x[k] * id_99[k]
                   + kd_99[k];

        t_164[k] = ab_x[k] * id_100[k]
                   + kd_100[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, ab_y, ab_z, id_99, id_100, \
                         id_101, kd_101, kd_135, kd_136, kd_137, \
                         kd_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_165[k] = ab_x[k] * id_101[k]
                   + kd_101[k];

        t_166[k] = ab_y[k] * id_99[k]
                   + kd_135[k];

        t_167[k] = ab_y[k] * id_100[k]
                   + kd_136[k];

        t_168[k] = ab_y[k] * id_101[k]
                   + kd_137[k];

        t_169[k] = ab_z[k] * id_101[k]
                   + kd_143[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, id_102, id_103, id_104, \
                         id_105, id_106, kd_102, kd_103, kd_104, kd_105, \
                         kd_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_170[k] = ab_x[k] * id_102[k]
                   + kd_102[k];

        t_171[k] = ab_x[k] * id_103[k]
                   + kd_103[k];

        t_172[k] = ab_x[k] * id_104[k]
                   + kd_104[k];

        t_173[k] = ab_x[k] * id_105[k]
                   + kd_105[k];

        t_174[k] = ab_x[k] * id_106[k]
                   + kd_106[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, ab_y, ab_z, id_105, id_106, \
                         id_107, kd_107, kd_141, kd_142, kd_143, \
                         kd_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_175[k] = ab_x[k] * id_107[k]
                   + kd_107[k];

        t_176[k] = ab_y[k] * id_105[k]
                   + kd_141[k];

        t_177[k] = ab_y[k] * id_106[k]
                   + kd_142[k];

        t_178[k] = ab_y[k] * id_107[k]
                   + kd_143[k];

        t_179[k] = ab_z[k] * id_107[k]
                   + kd_149[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, id_108, id_109, id_110, \
                         id_111, id_112, kd_108, kd_109, kd_110, kd_111, \
                         kd_112 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_180[k] = ab_x[k] * id_108[k]
                   + kd_108[k];

        t_181[k] = ab_x[k] * id_109[k]
                   + kd_109[k];

        t_182[k] = ab_x[k] * id_110[k]
                   + kd_110[k];

        t_183[k] = ab_x[k] * id_111[k]
                   + kd_111[k];

        t_184[k] = ab_x[k] * id_112[k]
                   + kd_112[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, ab_y, ab_z, id_111, id_112, \
                         id_113, kd_113, kd_147, kd_148, kd_149, \
                         kd_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_185[k] = ab_x[k] * id_113[k]
                   + kd_113[k];

        t_186[k] = ab_y[k] * id_111[k]
                   + kd_147[k];

        t_187[k] = ab_y[k] * id_112[k]
                   + kd_148[k];

        t_188[k] = ab_y[k] * id_113[k]
                   + kd_149[k];

        t_189[k] = ab_z[k] * id_113[k]
                   + kd_155[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, id_114, id_115, id_116, \
                         id_117, id_118, kd_114, kd_115, kd_116, kd_117, \
                         kd_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_190[k] = ab_x[k] * id_114[k]
                   + kd_114[k];

        t_191[k] = ab_x[k] * id_115[k]
                   + kd_115[k];

        t_192[k] = ab_x[k] * id_116[k]
                   + kd_116[k];

        t_193[k] = ab_x[k] * id_117[k]
                   + kd_117[k];

        t_194[k] = ab_x[k] * id_118[k]
                   + kd_118[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, ab_y, ab_z, id_117, id_118, \
                         id_119, kd_119, kd_153, kd_154, kd_155, \
                         kd_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_195[k] = ab_x[k] * id_119[k]
                   + kd_119[k];

        t_196[k] = ab_y[k] * id_117[k]
                   + kd_153[k];

        t_197[k] = ab_y[k] * id_118[k]
                   + kd_154[k];

        t_198[k] = ab_y[k] * id_119[k]
                   + kd_155[k];

        t_199[k] = ab_z[k] * id_119[k]
                   + kd_161[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, id_120, id_121, id_122, \
                         id_123, id_124, kd_120, kd_121, kd_122, kd_123, \
                         kd_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_200[k] = ab_x[k] * id_120[k]
                   + kd_120[k];

        t_201[k] = ab_x[k] * id_121[k]
                   + kd_121[k];

        t_202[k] = ab_x[k] * id_122[k]
                   + kd_122[k];

        t_203[k] = ab_x[k] * id_123[k]
                   + kd_123[k];

        t_204[k] = ab_x[k] * id_124[k]
                   + kd_124[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, ab_y, ab_z, id_123, id_124, \
                         id_125, kd_125, kd_159, kd_160, kd_161, \
                         kd_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_205[k] = ab_x[k] * id_125[k]
                   + kd_125[k];

        t_206[k] = ab_y[k] * id_123[k]
                   + kd_159[k];

        t_207[k] = ab_y[k] * id_124[k]
                   + kd_160[k];

        t_208[k] = ab_y[k] * id_125[k]
                   + kd_161[k];

        t_209[k] = ab_z[k] * id_125[k]
                   + kd_167[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, id_126, id_127, id_128, \
                         id_129, id_130, kd_126, kd_127, kd_128, kd_129, \
                         kd_130 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_210[k] = ab_x[k] * id_126[k]
                   + kd_126[k];

        t_211[k] = ab_x[k] * id_127[k]
                   + kd_127[k];

        t_212[k] = ab_x[k] * id_128[k]
                   + kd_128[k];

        t_213[k] = ab_x[k] * id_129[k]
                   + kd_129[k];

        t_214[k] = ab_x[k] * id_130[k]
                   + kd_130[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, ab_y, ab_z, id_129, id_130, \
                         id_131, kd_131, kd_171, kd_172, kd_173, \
                         kd_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_215[k] = ab_x[k] * id_131[k]
                   + kd_131[k];

        t_216[k] = ab_y[k] * id_129[k]
                   + kd_171[k];

        t_217[k] = ab_y[k] * id_130[k]
                   + kd_172[k];

        t_218[k] = ab_y[k] * id_131[k]
                   + kd_173[k];

        t_219[k] = ab_z[k] * id_131[k]
                   + kd_179[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, id_132, id_133, id_134, \
                         id_135, id_136, kd_132, kd_133, kd_134, kd_135, \
                         kd_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_220[k] = ab_x[k] * id_132[k]
                   + kd_132[k];

        t_221[k] = ab_x[k] * id_133[k]
                   + kd_133[k];

        t_222[k] = ab_x[k] * id_134[k]
                   + kd_134[k];

        t_223[k] = ab_x[k] * id_135[k]
                   + kd_135[k];

        t_224[k] = ab_x[k] * id_136[k]
                   + kd_136[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, ab_y, ab_z, id_135, id_136, \
                         id_137, kd_137, kd_177, kd_178, kd_179, \
                         kd_185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_225[k] = ab_x[k] * id_137[k]
                   + kd_137[k];

        t_226[k] = ab_y[k] * id_135[k]
                   + kd_177[k];

        t_227[k] = ab_y[k] * id_136[k]
                   + kd_178[k];

        t_228[k] = ab_y[k] * id_137[k]
                   + kd_179[k];

        t_229[k] = ab_z[k] * id_137[k]
                   + kd_185[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, id_138, id_139, id_140, \
                         id_141, id_142, kd_138, kd_139, kd_140, kd_141, \
                         kd_142 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_230[k] = ab_x[k] * id_138[k]
                   + kd_138[k];

        t_231[k] = ab_x[k] * id_139[k]
                   + kd_139[k];

        t_232[k] = ab_x[k] * id_140[k]
                   + kd_140[k];

        t_233[k] = ab_x[k] * id_141[k]
                   + kd_141[k];

        t_234[k] = ab_x[k] * id_142[k]
                   + kd_142[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_x, ab_y, ab_z, id_141, id_142, \
                         id_143, kd_143, kd_183, kd_184, kd_185, \
                         kd_191 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_235[k] = ab_x[k] * id_143[k]
                   + kd_143[k];

        t_236[k] = ab_y[k] * id_141[k]
                   + kd_183[k];

        t_237[k] = ab_y[k] * id_142[k]
                   + kd_184[k];

        t_238[k] = ab_y[k] * id_143[k]
                   + kd_185[k];

        t_239[k] = ab_z[k] * id_143[k]
                   + kd_191[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, id_144, id_145, id_146, \
                         id_147, id_148, kd_144, kd_145, kd_146, kd_147, \
                         kd_148 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_240[k] = ab_x[k] * id_144[k]
                   + kd_144[k];

        t_241[k] = ab_x[k] * id_145[k]
                   + kd_145[k];

        t_242[k] = ab_x[k] * id_146[k]
                   + kd_146[k];

        t_243[k] = ab_x[k] * id_147[k]
                   + kd_147[k];

        t_244[k] = ab_x[k] * id_148[k]
                   + kd_148[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, ab_y, ab_z, id_147, id_148, \
                         id_149, kd_149, kd_189, kd_190, kd_191, \
                         kd_197 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_245[k] = ab_x[k] * id_149[k]
                   + kd_149[k];

        t_246[k] = ab_y[k] * id_147[k]
                   + kd_189[k];

        t_247[k] = ab_y[k] * id_148[k]
                   + kd_190[k];

        t_248[k] = ab_y[k] * id_149[k]
                   + kd_191[k];

        t_249[k] = ab_z[k] * id_149[k]
                   + kd_197[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_x, id_150, id_151, id_152, \
                         id_153, id_154, kd_150, kd_151, kd_152, kd_153, \
                         kd_154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_250[k] = ab_x[k] * id_150[k]
                   + kd_150[k];

        t_251[k] = ab_x[k] * id_151[k]
                   + kd_151[k];

        t_252[k] = ab_x[k] * id_152[k]
                   + kd_152[k];

        t_253[k] = ab_x[k] * id_153[k]
                   + kd_153[k];

        t_254[k] = ab_x[k] * id_154[k]
                   + kd_154[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, ab_y, ab_z, id_153, id_154, \
                         id_155, kd_155, kd_195, kd_196, kd_197, \
                         kd_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_255[k] = ab_x[k] * id_155[k]
                   + kd_155[k];

        t_256[k] = ab_y[k] * id_153[k]
                   + kd_195[k];

        t_257[k] = ab_y[k] * id_154[k]
                   + kd_196[k];

        t_258[k] = ab_y[k] * id_155[k]
                   + kd_197[k];

        t_259[k] = ab_z[k] * id_155[k]
                   + kd_203[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, id_156, id_157, id_158, \
                         id_159, id_160, kd_156, kd_157, kd_158, kd_159, \
                         kd_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_260[k] = ab_x[k] * id_156[k]
                   + kd_156[k];

        t_261[k] = ab_x[k] * id_157[k]
                   + kd_157[k];

        t_262[k] = ab_x[k] * id_158[k]
                   + kd_158[k];

        t_263[k] = ab_x[k] * id_159[k]
                   + kd_159[k];

        t_264[k] = ab_x[k] * id_160[k]
                   + kd_160[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, ab_y, ab_z, id_159, id_160, \
                         id_161, kd_161, kd_201, kd_202, kd_203, \
                         kd_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_265[k] = ab_x[k] * id_161[k]
                   + kd_161[k];

        t_266[k] = ab_y[k] * id_159[k]
                   + kd_201[k];

        t_267[k] = ab_y[k] * id_160[k]
                   + kd_202[k];

        t_268[k] = ab_y[k] * id_161[k]
                   + kd_203[k];

        t_269[k] = ab_z[k] * id_161[k]
                   + kd_209[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, id_162, id_163, id_164, \
                         id_165, id_166, kd_162, kd_163, kd_164, kd_165, \
                         kd_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_270[k] = ab_x[k] * id_162[k]
                   + kd_162[k];

        t_271[k] = ab_x[k] * id_163[k]
                   + kd_163[k];

        t_272[k] = ab_x[k] * id_164[k]
                   + kd_164[k];

        t_273[k] = ab_x[k] * id_165[k]
                   + kd_165[k];

        t_274[k] = ab_x[k] * id_166[k]
                   + kd_166[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, ab_y, ab_z, id_165, id_166, \
                         id_167, kd_167, kd_207, kd_208, kd_209, \
                         kd_215 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_275[k] = ab_x[k] * id_167[k]
                   + kd_167[k];

        t_276[k] = ab_y[k] * id_165[k]
                   + kd_207[k];

        t_277[k] = ab_y[k] * id_166[k]
                   + kd_208[k];

        t_278[k] = ab_y[k] * id_167[k]
                   + kd_209[k];

        t_279[k] = ab_z[k] * id_167[k]
                   + kd_215[k];
    }
}

auto
compute_hrr_if_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t id, const size_t kd,
                   const size_t nmax) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.28125 * std::sqrt(1155.0);
    const auto f_1 = 0.09375 * std::sqrt(1155.0);
    const auto f_2 = 0.9375 * std::sqrt(1155.0);
    const auto f_3 = 0.3125 * std::sqrt(1155.0);
    const auto f_4 = 0.5625 * std::sqrt(770.0);
    const auto f_5 = 1.875 * std::sqrt(770.0);
    const auto f_6 = 0.28125 * std::sqrt(77.0);
    const auto f_7 = 1.125 * std::sqrt(77.0);
    const auto f_8 = 0.9375 * std::sqrt(77.0);
    const auto f_9 = 3.75 * std::sqrt(77.0);
    const auto f_10 = 0.28125 * std::sqrt(462.0);
    const auto f_11 = 0.1875 * std::sqrt(462.0);
    const auto f_12 = 0.9375 * std::sqrt(462.0);
    const auto f_13 = 0.625 * std::sqrt(462.0);
    const auto f_14 = 0.28125 * std::sqrt(770.0);
    const auto f_15 = 0.9375 * std::sqrt(770.0);
    const auto f_16 = 1.40625 * std::sqrt(385.0);
    const auto f_17 = 0.46875 * std::sqrt(385.0);
    const auto f_18 = 2.8125 * std::sqrt(385.0);
    const auto f_19 = 0.9375 * std::sqrt(385.0);
    const auto f_20 = 0.28125 * std::sqrt(385.0);
    const auto f_21 = 0.09375 * std::sqrt(385.0);
    const auto f_22 = 0.9375 * std::sqrt(2310.0);
    const auto f_23 = 1.875 * std::sqrt(2310.0);
    const auto f_24 = 0.1875 * std::sqrt(2310.0);
    const auto f_25 = 0.46875 * std::sqrt(231.0);
    const auto f_26 = 1.875 * std::sqrt(231.0);
    const auto f_27 = 0.9375 * std::sqrt(231.0);
    const auto f_28 = 3.75 * std::sqrt(231.0);
    const auto f_29 = 0.09375 * std::sqrt(231.0);
    const auto f_30 = 0.375 * std::sqrt(231.0);
    const auto f_31 = 1.40625 * std::sqrt(154.0);
    const auto f_32 = 0.9375 * std::sqrt(154.0);
    const auto f_33 = 2.8125 * std::sqrt(154.0);
    const auto f_34 = 1.875 * std::sqrt(154.0);
    const auto f_35 = 0.28125 * std::sqrt(154.0);
    const auto f_36 = 0.1875 * std::sqrt(154.0);
    const auto f_37 = 0.46875 * std::sqrt(2310.0);
    const auto f_38 = 0.09375 * std::sqrt(2310.0);
    const auto f_39 = 0.5625 * std::sqrt(70.0);
    const auto f_40 = 0.1875 * std::sqrt(70.0);
    const auto f_41 = 5.625 * std::sqrt(70.0);
    const auto f_42 = 1.875 * std::sqrt(70.0);
    const auto f_43 = 0.75 * std::sqrt(105.0);
    const auto f_44 = 7.5 * std::sqrt(105.0);
    const auto f_45 = 0.1875 * std::sqrt(42.0);
    const auto f_46 = 0.75 * std::sqrt(42.0);
    const auto f_47 = 1.875 * std::sqrt(42.0);
    const auto f_48 = 7.5 * std::sqrt(42.0);
    const auto f_49 = 1.125 * std::sqrt(7.0);
    const auto f_50 = 0.75 * std::sqrt(7.0);
    const auto f_51 = 11.25 * std::sqrt(7.0);
    const auto f_52 = 7.5 * std::sqrt(7.0);
    const auto f_53 = 0.375 * std::sqrt(105.0);
    const auto f_54 = 3.75 * std::sqrt(105.0);
    const auto f_55 = 4.21875 * std::sqrt(21.0);
    const auto f_56 = 1.40625 * std::sqrt(21.0);
    const auto f_57 = 2.8125 * std::sqrt(21.0);
    const auto f_58 = 0.9375 * std::sqrt(21.0);
    const auto f_59 = 11.25 * std::sqrt(21.0);
    const auto f_60 = 3.75 * std::sqrt(21.0);
    const auto f_61 = 0.46875 * std::sqrt(21.0);
    const auto f_62 = 1.25 * std::sqrt(21.0);
    const auto f_63 = 8.4375 * std::sqrt(14.0);
    const auto f_64 = 5.625 * std::sqrt(14.0);
    const auto f_65 = 22.5 * std::sqrt(14.0);
    const auto f_66 = 2.8125 * std::sqrt(14.0);
    const auto f_67 = 7.5 * std::sqrt(14.0);
    const auto f_68 = 0.84375 * std::sqrt(35.0);
    const auto f_69 = 3.375 * std::sqrt(35.0);
    const auto f_70 = 0.5625 * std::sqrt(35.0);
    const auto f_71 = 2.25 * std::sqrt(35.0);
    const auto f_72 = 9.0 * std::sqrt(35.0);
    const auto f_73 = 0.28125 * std::sqrt(35.0);
    const auto f_74 = 1.125 * std::sqrt(35.0);
    const auto f_75 = 0.75 * std::sqrt(35.0);
    const auto f_76 = 3.0 * std::sqrt(35.0);
    const auto f_77 = 0.84375 * std::sqrt(210.0);
    const auto f_78 = 0.5625 * std::sqrt(210.0);
    const auto f_79 = 0.375 * std::sqrt(210.0);
    const auto f_80 = 2.25 * std::sqrt(210.0);
    const auto f_81 = 1.5 * std::sqrt(210.0);
    const auto f_82 = 0.28125 * std::sqrt(210.0);
    const auto f_83 = 0.1875 * std::sqrt(210.0);
    const auto f_84 = 0.75 * std::sqrt(210.0);
    const auto f_85 = 0.5 * std::sqrt(210.0);
    const auto f_86 = 4.21875 * std::sqrt(14.0);
    const auto f_87 = 11.25 * std::sqrt(14.0);
    const auto f_88 = 1.40625 * std::sqrt(14.0);
    const auto f_89 = 3.75 * std::sqrt(14.0);
    const auto f_90 = 0.15625 * std::sqrt(21.0);
    const auto f_91 = 0.3125 * std::sqrt(21.0);
    const auto f_92 = 7.5 * std::sqrt(21.0);
    const auto f_93 = 2.5 * std::sqrt(21.0);
    const auto f_94 = 0.9375 * std::sqrt(14.0);
    const auto f_95 = 1.875 * std::sqrt(14.0);
    const auto f_96 = 15.0 * std::sqrt(14.0);
    const auto f_97 = 0.09375 * std::sqrt(35.0);
    const auto f_98 = 0.375 * std::sqrt(35.0);
    const auto f_99 = 0.1875 * std::sqrt(35.0);
    const auto f_100 = 1.5 * std::sqrt(35.0);
    const auto f_101 = 6.0 * std::sqrt(35.0);
    const auto f_102 = 0.09375 * std::sqrt(210.0);
    const auto f_103 = 0.0625 * std::sqrt(210.0);
    const auto f_104 = 0.125 * std::sqrt(210.0);
    const auto f_105 = std::sqrt(210.0);
    const auto f_106 = 0.46875 * std::sqrt(14.0);
    const auto f_107 = 0.46875 * std::sqrt(210.0);
    const auto f_108 = 0.15625 * std::sqrt(210.0);
    const auto f_109 = 0.9375 * std::sqrt(210.0);
    const auto f_110 = 0.3125 * std::sqrt(210.0);
    const auto f_111 = 1.875 * std::sqrt(210.0);
    const auto f_112 = 0.625 * std::sqrt(210.0);
    const auto f_113 = 0.25 * std::sqrt(210.0);
    const auto f_114 = 1.875 * std::sqrt(35.0);
    const auto f_115 = 3.75 * std::sqrt(35.0);
    const auto f_116 = 7.5 * std::sqrt(35.0);
    const auto f_117 = 0.75 * std::sqrt(14.0);
    const auto f_118 = 3.0 * std::sqrt(14.0);
    const auto f_119 = 0.625 * std::sqrt(21.0);
    const auto f_120 = 1.875 * std::sqrt(21.0);
    const auto f_121 = 1.5 * std::sqrt(21.0);
    const auto f_122 = std::sqrt(21.0);
    const auto f_123 = 0.9375 * std::sqrt(35.0);
    const auto f_124 = 0.234375 * std::sqrt(10.0);
    const auto f_125 = 0.078125 * std::sqrt(10.0);
    const auto f_126 = 0.703125 * std::sqrt(10.0);
    const auto f_127 = 4.21875 * std::sqrt(10.0);
    const auto f_128 = 1.40625 * std::sqrt(10.0);
    const auto f_129 = 8.4375 * std::sqrt(10.0);
    const auto f_130 = 2.8125 * std::sqrt(10.0);
    const auto f_131 = 5.625 * std::sqrt(10.0);
    const auto f_132 = 1.875 * std::sqrt(10.0);
    const auto f_133 = 0.75 * std::sqrt(10.0);
    const auto f_134 = 0.25 * std::sqrt(10.0);
    const auto f_135 = 0.3125 * std::sqrt(15.0);
    const auto f_136 = 0.9375 * std::sqrt(15.0);
    const auto f_137 = 5.625 * std::sqrt(15.0);
    const auto f_138 = 11.25 * std::sqrt(15.0);
    const auto f_139 = 7.5 * std::sqrt(15.0);
    const auto f_140 = std::sqrt(15.0);
    const auto f_141 = 0.078125 * std::sqrt(6.0);
    const auto f_142 = 0.3125 * std::sqrt(6.0);
    const auto f_143 = 0.234375 * std::sqrt(6.0);
    const auto f_144 = 0.9375 * std::sqrt(6.0);
    const auto f_145 = 1.40625 * std::sqrt(6.0);
    const auto f_146 = 5.625 * std::sqrt(6.0);
    const auto f_147 = 2.8125 * std::sqrt(6.0);
    const auto f_148 = 11.25 * std::sqrt(6.0);
    const auto f_149 = 1.875 * std::sqrt(6.0);
    const auto f_150 = 7.5 * std::sqrt(6.0);
    const auto f_151 = 0.25 * std::sqrt(6.0);
    const auto f_152 = std::sqrt(6.0);
    const auto f_153 = 0.15625 * std::sqrt(15.0);
    const auto f_154 = 0.46875 * std::sqrt(15.0);
    const auto f_155 = 2.8125 * std::sqrt(15.0);
    const auto f_156 = 3.75 * std::sqrt(15.0);
    const auto f_157 = 0.5 * std::sqrt(15.0);
    const auto f_158 = 0.234375 * std::sqrt(21.0);
    const auto f_159 = 0.078125 * std::sqrt(21.0);
    const auto f_160 = 0.046875 * std::sqrt(35.0);
    const auto f_161 = 0.046875 * std::sqrt(210.0);
    const auto f_162 = 0.03125 * std::sqrt(210.0);
    const auto f_163 = 0.234375 * std::sqrt(14.0);
    const auto f_164 = 0.140625 * std::sqrt(70.0);
    const auto f_165 = 0.046875 * std::sqrt(70.0);
    const auto f_166 = 0.703125 * std::sqrt(70.0);
    const auto f_167 = 0.234375 * std::sqrt(70.0);
    const auto f_168 = 1.40625 * std::sqrt(70.0);
    const auto f_169 = 0.46875 * std::sqrt(70.0);
    const auto f_170 = 8.4375 * std::sqrt(70.0);
    const auto f_171 = 2.8125 * std::sqrt(70.0);
    const auto f_172 = 0.1875 * std::sqrt(105.0);
    const auto f_173 = 0.9375 * std::sqrt(105.0);
    const auto f_174 = 1.875 * std::sqrt(105.0);
    const auto f_175 = 11.25 * std::sqrt(105.0);
    const auto f_176 = 0.046875 * std::sqrt(42.0);
    const auto f_177 = 0.234375 * std::sqrt(42.0);
    const auto f_178 = 0.9375 * std::sqrt(42.0);
    const auto f_179 = 0.46875 * std::sqrt(42.0);
    const auto f_180 = 2.8125 * std::sqrt(42.0);
    const auto f_181 = 11.25 * std::sqrt(42.0);
    const auto f_182 = 0.28125 * std::sqrt(7.0);
    const auto f_183 = 0.1875 * std::sqrt(7.0);
    const auto f_184 = 1.40625 * std::sqrt(7.0);
    const auto f_185 = 0.9375 * std::sqrt(7.0);
    const auto f_186 = 2.8125 * std::sqrt(7.0);
    const auto f_187 = 1.875 * std::sqrt(7.0);
    const auto f_188 = 16.875 * std::sqrt(7.0);
    const auto f_189 = 0.09375 * std::sqrt(105.0);
    const auto f_190 = 0.46875 * std::sqrt(105.0);
    const auto f_191 = 5.625 * std::sqrt(105.0);
    const auto f_192 = 0.046875 * std::sqrt(1155.0);
    const auto f_193 = 0.015625 * std::sqrt(1155.0);
    const auto f_194 = 0.703125 * std::sqrt(1155.0);
    const auto f_195 = 0.234375 * std::sqrt(1155.0);
    const auto f_196 = 0.09375 * std::sqrt(770.0);
    const auto f_197 = 1.40625 * std::sqrt(770.0);
    const auto f_198 = 0.046875 * std::sqrt(77.0);
    const auto f_199 = 0.1875 * std::sqrt(77.0);
    const auto f_200 = 0.703125 * std::sqrt(77.0);
    const auto f_201 = 2.8125 * std::sqrt(77.0);
    const auto f_202 = 0.046875 * std::sqrt(462.0);
    const auto f_203 = 0.03125 * std::sqrt(462.0);
    const auto f_204 = 0.703125 * std::sqrt(462.0);
    const auto f_205 = 0.46875 * std::sqrt(462.0);
    const auto f_206 = 0.046875 * std::sqrt(770.0);
    const auto f_207 = 0.703125 * std::sqrt(770.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;
    auto *g_7 = values + 7 * nvalues;
    auto *g_8 = values + 8 * nvalues;
    auto *g_9 = values + 9 * nvalues;
    auto *g_10 = values + 10 * nvalues;
    auto *g_11 = values + 11 * nvalues;
    auto *g_12 = values + 12 * nvalues;
    auto *g_13 = values + 13 * nvalues;
    auto *g_14 = values + 14 * nvalues;
    auto *g_15 = values + 15 * nvalues;
    auto *g_16 = values + 16 * nvalues;
    auto *g_17 = values + 17 * nvalues;
    auto *g_18 = values + 18 * nvalues;
    auto *g_19 = values + 19 * nvalues;
    auto *g_20 = values + 20 * nvalues;
    auto *g_21 = values + 21 * nvalues;
    auto *g_22 = values + 22 * nvalues;
    auto *g_23 = values + 23 * nvalues;
    auto *g_24 = values + 24 * nvalues;
    auto *g_25 = values + 25 * nvalues;
    auto *g_26 = values + 26 * nvalues;
    auto *g_27 = values + 27 * nvalues;
    auto *g_28 = values + 28 * nvalues;
    auto *g_29 = values + 29 * nvalues;
    auto *g_30 = values + 30 * nvalues;
    auto *g_31 = values + 31 * nvalues;
    auto *g_32 = values + 32 * nvalues;
    auto *g_33 = values + 33 * nvalues;
    auto *g_34 = values + 34 * nvalues;
    auto *g_35 = values + 35 * nvalues;
    auto *g_36 = values + 36 * nvalues;
    auto *g_37 = values + 37 * nvalues;
    auto *g_38 = values + 38 * nvalues;
    auto *g_39 = values + 39 * nvalues;
    auto *g_40 = values + 40 * nvalues;
    auto *g_41 = values + 41 * nvalues;
    auto *g_42 = values + 42 * nvalues;
    auto *g_43 = values + 43 * nvalues;
    auto *g_44 = values + 44 * nvalues;
    auto *g_45 = values + 45 * nvalues;
    auto *g_46 = values + 46 * nvalues;
    auto *g_47 = values + 47 * nvalues;
    auto *g_48 = values + 48 * nvalues;
    auto *g_49 = values + 49 * nvalues;
    auto *g_50 = values + 50 * nvalues;
    auto *g_51 = values + 51 * nvalues;
    auto *g_52 = values + 52 * nvalues;
    auto *g_53 = values + 53 * nvalues;
    auto *g_54 = values + 54 * nvalues;
    auto *g_55 = values + 55 * nvalues;
    auto *g_56 = values + 56 * nvalues;
    auto *g_57 = values + 57 * nvalues;
    auto *g_58 = values + 58 * nvalues;
    auto *g_59 = values + 59 * nvalues;
    auto *g_60 = values + 60 * nvalues;
    auto *g_61 = values + 61 * nvalues;
    auto *g_62 = values + 62 * nvalues;
    auto *g_63 = values + 63 * nvalues;
    auto *g_64 = values + 64 * nvalues;
    auto *g_65 = values + 65 * nvalues;
    auto *g_66 = values + 66 * nvalues;
    auto *g_67 = values + 67 * nvalues;
    auto *g_68 = values + 68 * nvalues;
    auto *g_69 = values + 69 * nvalues;
    auto *g_70 = values + 70 * nvalues;
    auto *g_71 = values + 71 * nvalues;
    auto *g_72 = values + 72 * nvalues;
    auto *g_73 = values + 73 * nvalues;
    auto *g_74 = values + 74 * nvalues;
    auto *g_75 = values + 75 * nvalues;
    auto *g_76 = values + 76 * nvalues;
    auto *g_77 = values + 77 * nvalues;
    auto *g_78 = values + 78 * nvalues;
    auto *g_79 = values + 79 * nvalues;
    auto *g_80 = values + 80 * nvalues;
    auto *g_81 = values + 81 * nvalues;
    auto *g_82 = values + 82 * nvalues;
    auto *g_83 = values + 83 * nvalues;
    auto *g_84 = values + 84 * nvalues;
    auto *g_85 = values + 85 * nvalues;
    auto *g_86 = values + 86 * nvalues;
    auto *g_87 = values + 87 * nvalues;
    auto *g_88 = values + 88 * nvalues;
    auto *g_89 = values + 89 * nvalues;
    auto *g_90 = values + 90 * nvalues;

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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
    const auto *id_109 = buffer.data(id + 109);
    const auto *id_110 = buffer.data(id + 110);
    const auto *id_111 = buffer.data(id + 111);
    const auto *id_112 = buffer.data(id + 112);
    const auto *id_113 = buffer.data(id + 113);
    const auto *id_114 = buffer.data(id + 114);
    const auto *id_115 = buffer.data(id + 115);
    const auto *id_116 = buffer.data(id + 116);
    const auto *id_117 = buffer.data(id + 117);
    const auto *id_118 = buffer.data(id + 118);
    const auto *id_119 = buffer.data(id + 119);
    const auto *id_120 = buffer.data(id + 120);
    const auto *id_121 = buffer.data(id + 121);
    const auto *id_122 = buffer.data(id + 122);
    const auto *id_123 = buffer.data(id + 123);
    const auto *id_124 = buffer.data(id + 124);
    const auto *id_125 = buffer.data(id + 125);
    const auto *id_126 = buffer.data(id + 126);
    const auto *id_127 = buffer.data(id + 127);
    const auto *id_128 = buffer.data(id + 128);
    const auto *id_129 = buffer.data(id + 129);
    const auto *id_130 = buffer.data(id + 130);
    const auto *id_131 = buffer.data(id + 131);
    const auto *id_132 = buffer.data(id + 132);
    const auto *id_133 = buffer.data(id + 133);
    const auto *id_134 = buffer.data(id + 134);
    const auto *id_135 = buffer.data(id + 135);
    const auto *id_136 = buffer.data(id + 136);
    const auto *id_137 = buffer.data(id + 137);
    const auto *id_138 = buffer.data(id + 138);
    const auto *id_139 = buffer.data(id + 139);
    const auto *id_140 = buffer.data(id + 140);
    const auto *id_141 = buffer.data(id + 141);
    const auto *id_142 = buffer.data(id + 142);
    const auto *id_143 = buffer.data(id + 143);
    const auto *id_144 = buffer.data(id + 144);
    const auto *id_145 = buffer.data(id + 145);
    const auto *id_146 = buffer.data(id + 146);
    const auto *id_147 = buffer.data(id + 147);
    const auto *id_148 = buffer.data(id + 148);
    const auto *id_149 = buffer.data(id + 149);
    const auto *id_150 = buffer.data(id + 150);
    const auto *id_151 = buffer.data(id + 151);
    const auto *id_152 = buffer.data(id + 152);
    const auto *id_153 = buffer.data(id + 153);
    const auto *id_154 = buffer.data(id + 154);
    const auto *id_155 = buffer.data(id + 155);
    const auto *id_156 = buffer.data(id + 156);
    const auto *id_157 = buffer.data(id + 157);
    const auto *id_158 = buffer.data(id + 158);
    const auto *id_159 = buffer.data(id + 159);
    const auto *id_160 = buffer.data(id + 160);
    const auto *id_161 = buffer.data(id + 161);
    const auto *id_162 = buffer.data(id + 162);
    const auto *id_163 = buffer.data(id + 163);
    const auto *id_164 = buffer.data(id + 164);
    const auto *id_165 = buffer.data(id + 165);
    const auto *id_166 = buffer.data(id + 166);
    const auto *id_167 = buffer.data(id + 167);

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
    const auto *kd_93 = buffer.data(kd + 93);
    const auto *kd_94 = buffer.data(kd + 94);
    const auto *kd_95 = buffer.data(kd + 95);
    const auto *kd_96 = buffer.data(kd + 96);
    const auto *kd_97 = buffer.data(kd + 97);
    const auto *kd_98 = buffer.data(kd + 98);
    const auto *kd_99 = buffer.data(kd + 99);
    const auto *kd_100 = buffer.data(kd + 100);
    const auto *kd_101 = buffer.data(kd + 101);
    const auto *kd_102 = buffer.data(kd + 102);
    const auto *kd_103 = buffer.data(kd + 103);
    const auto *kd_104 = buffer.data(kd + 104);
    const auto *kd_105 = buffer.data(kd + 105);
    const auto *kd_106 = buffer.data(kd + 106);
    const auto *kd_107 = buffer.data(kd + 107);
    const auto *kd_108 = buffer.data(kd + 108);
    const auto *kd_109 = buffer.data(kd + 109);
    const auto *kd_110 = buffer.data(kd + 110);
    const auto *kd_111 = buffer.data(kd + 111);
    const auto *kd_112 = buffer.data(kd + 112);
    const auto *kd_113 = buffer.data(kd + 113);
    const auto *kd_114 = buffer.data(kd + 114);
    const auto *kd_115 = buffer.data(kd + 115);
    const auto *kd_116 = buffer.data(kd + 116);
    const auto *kd_117 = buffer.data(kd + 117);
    const auto *kd_118 = buffer.data(kd + 118);
    const auto *kd_119 = buffer.data(kd + 119);
    const auto *kd_120 = buffer.data(kd + 120);
    const auto *kd_121 = buffer.data(kd + 121);
    const auto *kd_122 = buffer.data(kd + 122);
    const auto *kd_123 = buffer.data(kd + 123);
    const auto *kd_124 = buffer.data(kd + 124);
    const auto *kd_125 = buffer.data(kd + 125);
    const auto *kd_126 = buffer.data(kd + 126);
    const auto *kd_127 = buffer.data(kd + 127);
    const auto *kd_128 = buffer.data(kd + 128);
    const auto *kd_129 = buffer.data(kd + 129);
    const auto *kd_130 = buffer.data(kd + 130);
    const auto *kd_131 = buffer.data(kd + 131);
    const auto *kd_132 = buffer.data(kd + 132);
    const auto *kd_133 = buffer.data(kd + 133);
    const auto *kd_134 = buffer.data(kd + 134);
    const auto *kd_135 = buffer.data(kd + 135);
    const auto *kd_136 = buffer.data(kd + 136);
    const auto *kd_137 = buffer.data(kd + 137);
    const auto *kd_138 = buffer.data(kd + 138);
    const auto *kd_139 = buffer.data(kd + 139);
    const auto *kd_140 = buffer.data(kd + 140);
    const auto *kd_141 = buffer.data(kd + 141);
    const auto *kd_142 = buffer.data(kd + 142);
    const auto *kd_143 = buffer.data(kd + 143);
    const auto *kd_144 = buffer.data(kd + 144);
    const auto *kd_145 = buffer.data(kd + 145);
    const auto *kd_146 = buffer.data(kd + 146);
    const auto *kd_147 = buffer.data(kd + 147);
    const auto *kd_148 = buffer.data(kd + 148);
    const auto *kd_149 = buffer.data(kd + 149);
    const auto *kd_150 = buffer.data(kd + 150);
    const auto *kd_151 = buffer.data(kd + 151);
    const auto *kd_152 = buffer.data(kd + 152);
    const auto *kd_153 = buffer.data(kd + 153);
    const auto *kd_154 = buffer.data(kd + 154);
    const auto *kd_155 = buffer.data(kd + 155);
    const auto *kd_156 = buffer.data(kd + 156);
    const auto *kd_157 = buffer.data(kd + 157);
    const auto *kd_158 = buffer.data(kd + 158);
    const auto *kd_159 = buffer.data(kd + 159);
    const auto *kd_160 = buffer.data(kd + 160);
    const auto *kd_161 = buffer.data(kd + 161);
    const auto *kd_162 = buffer.data(kd + 162);
    const auto *kd_163 = buffer.data(kd + 163);
    const auto *kd_164 = buffer.data(kd + 164);
    const auto *kd_165 = buffer.data(kd + 165);
    const auto *kd_166 = buffer.data(kd + 166);
    const auto *kd_167 = buffer.data(kd + 167);
    const auto *kd_171 = buffer.data(kd + 171);
    const auto *kd_172 = buffer.data(kd + 172);
    const auto *kd_173 = buffer.data(kd + 173);
    const auto *kd_177 = buffer.data(kd + 177);
    const auto *kd_178 = buffer.data(kd + 178);
    const auto *kd_179 = buffer.data(kd + 179);
    const auto *kd_183 = buffer.data(kd + 183);
    const auto *kd_184 = buffer.data(kd + 184);
    const auto *kd_185 = buffer.data(kd + 185);
    const auto *kd_189 = buffer.data(kd + 189);
    const auto *kd_190 = buffer.data(kd + 190);
    const auto *kd_191 = buffer.data(kd + 191);
    const auto *kd_195 = buffer.data(kd + 195);
    const auto *kd_196 = buffer.data(kd + 196);
    const auto *kd_197 = buffer.data(kd + 197);
    const auto *kd_201 = buffer.data(kd + 201);
    const auto *kd_202 = buffer.data(kd + 202);
    const auto *kd_203 = buffer.data(kd + 203);
    const auto *kd_207 = buffer.data(kd + 207);
    const auto *kd_208 = buffer.data(kd + 208);
    const auto *kd_209 = buffer.data(kd + 209);
    const auto *kd_215 = buffer.data(kd + 215);

#pragma omp simd aligned(ab_x, ab_y, id_7, id_9, id_37, id_39, id_91, id_93, kd_7, kd_21, \
                         kd_37, kd_63, kd_91, kd_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * ab_x[k] * id_7[k]
                 - f_1 * ab_y[k] * id_9[k]
                 - f_2 * ab_x[k] * id_37[k]
                 + f_3 * ab_y[k] * id_39[k]
                 + f_0 * ab_x[k] * id_91[k]
                 - f_1 * ab_y[k] * id_93[k]
                 + f_0 * kd_7[k]
                 - f_1 * kd_21[k]
                 - f_2 * kd_37[k]
                 + f_3 * kd_63[k]
                 + f_0 * kd_91[k]
                 - f_1 * kd_129[k];
    }

#pragma omp simd aligned(ab_x, id_10, id_40, id_94, kd_10, kd_40, \
                         kd_94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = f_4 * ab_x[k] * id_10[k]
                 - f_5 * ab_x[k] * id_40[k]
                 + f_4 * ab_x[k] * id_94[k]
                 + f_4 * kd_10[k]
                 - f_5 * kd_40[k]
                 + f_4 * kd_94[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_7, id_9, id_11, id_37, id_39, id_41, id_91, id_93, \
                         id_95, kd_7, kd_21, kd_23, kd_37, kd_63, kd_65, kd_91, kd_129, \
                         kd_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_6 * ab_x[k] * id_7[k]
                 - f_6 * ab_y[k] * id_9[k]
                 + f_7 * ab_y[k] * id_11[k]
                 + f_8 * ab_x[k] * id_37[k]
                 + f_8 * ab_y[k] * id_39[k]
                 - f_9 * ab_y[k] * id_41[k]
                 - f_6 * ab_x[k] * id_91[k]
                 - f_6 * ab_y[k] * id_93[k]
                 + f_7 * ab_y[k] * id_95[k]
                 - f_6 * kd_7[k]
                 - f_6 * kd_21[k]
                 + f_7 * kd_23[k]
                 + f_8 * kd_37[k]
                 + f_8 * kd_63[k]
                 - f_9 * kd_65[k]
                 - f_6 * kd_91[k]
                 - f_6 * kd_129[k]
                 + f_7 * kd_131[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, id_8, id_10, id_11, id_38, id_40, id_41, id_92, \
                         id_94, id_95, kd_8, kd_22, kd_29, kd_38, kd_64, kd_71, kd_92, kd_130, \
                         kd_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_10 * ab_x[k] * id_8[k]
                 - f_10 * ab_y[k] * id_10[k]
                 + f_11 * ab_z[k] * id_11[k]
                 + f_12 * ab_x[k] * id_38[k]
                 + f_12 * ab_y[k] * id_40[k]
                 - f_13 * ab_z[k] * id_41[k]
                 - f_10 * ab_x[k] * id_92[k]
                 - f_10 * ab_y[k] * id_94[k]
                 + f_11 * ab_z[k] * id_95[k]
                 - f_10 * kd_8[k]
                 - f_10 * kd_22[k]
                 + f_11 * kd_29[k]
                 + f_12 * kd_38[k]
                 + f_12 * kd_64[k]
                 - f_13 * kd_71[k]
                 - f_10 * kd_92[k]
                 - f_10 * kd_130[k]
                 + f_11 * kd_137[k];
    }

#pragma omp simd aligned(ab_x, id_6, id_9, id_11, id_36, id_39, id_41, id_90, id_93, id_95, \
                         kd_6, kd_9, kd_11, kd_36, kd_39, kd_41, kd_90, kd_93, \
                         kd_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_6 * ab_x[k] * id_6[k]
                 - f_6 * ab_x[k] * id_9[k]
                 + f_7 * ab_x[k] * id_11[k]
                 + f_8 * ab_x[k] * id_36[k]
                 + f_8 * ab_x[k] * id_39[k]
                 - f_9 * ab_x[k] * id_41[k]
                 - f_6 * ab_x[k] * id_90[k]
                 - f_6 * ab_x[k] * id_93[k]
                 + f_7 * ab_x[k] * id_95[k]
                 - f_6 * kd_6[k]
                 - f_6 * kd_9[k]
                 + f_7 * kd_11[k]
                 + f_8 * kd_36[k]
                 + f_8 * kd_39[k]
                 - f_9 * kd_41[k]
                 - f_6 * kd_90[k]
                 - f_6 * kd_93[k]
                 + f_7 * kd_95[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_8, id_10, id_38, id_40, id_92, id_94, kd_8, kd_22, \
                         kd_38, kd_64, kd_92, kd_130 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_14 * ab_x[k] * id_8[k]
                 - f_14 * ab_y[k] * id_10[k]
                 - f_15 * ab_x[k] * id_38[k]
                 + f_15 * ab_y[k] * id_40[k]
                 + f_14 * ab_x[k] * id_92[k]
                 - f_14 * ab_y[k] * id_94[k]
                 + f_14 * kd_8[k]
                 - f_14 * kd_22[k]
                 - f_15 * kd_38[k]
                 + f_15 * kd_64[k]
                 + f_14 * kd_92[k]
                 - f_14 * kd_130[k];
    }

#pragma omp simd aligned(ab_x, id_6, id_9, id_36, id_39, id_90, id_93, kd_6, kd_9, kd_36, \
                         kd_39, kd_90, kd_93 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_1 * ab_x[k] * id_6[k]
                 - f_0 * ab_x[k] * id_9[k]
                 - f_3 * ab_x[k] * id_36[k]
                 + f_2 * ab_x[k] * id_39[k]
                 + f_1 * ab_x[k] * id_90[k]
                 - f_0 * ab_x[k] * id_93[k]
                 + f_1 * kd_6[k]
                 - f_0 * kd_9[k]
                 - f_3 * kd_36[k]
                 + f_2 * kd_39[k]
                 + f_1 * kd_90[k]
                 - f_0 * kd_93[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_25, id_27, id_67, id_69, id_133, id_135, kd_25, kd_45, \
                         kd_67, kd_99, kd_133, kd_177 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_16 * ab_x[k] * id_25[k]
                 - f_17 * ab_y[k] * id_27[k]
                 - f_18 * ab_x[k] * id_67[k]
                 + f_19 * ab_y[k] * id_69[k]
                 + f_20 * ab_x[k] * id_133[k]
                 - f_21 * ab_y[k] * id_135[k]
                 + f_16 * kd_25[k]
                 - f_17 * kd_45[k]
                 - f_18 * kd_67[k]
                 + f_19 * kd_99[k]
                 + f_20 * kd_133[k]
                 - f_21 * kd_177[k];
    }

#pragma omp simd aligned(ab_x, id_28, id_70, id_136, kd_28, kd_70, \
                         kd_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_22 * ab_x[k] * id_28[k]
                 - f_23 * ab_x[k] * id_70[k]
                 + f_24 * ab_x[k] * id_136[k]
                 + f_22 * kd_28[k]
                 - f_23 * kd_70[k]
                 + f_24 * kd_136[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_25, id_27, id_29, id_67, id_69, id_71, id_133, id_135, \
                         id_137, kd_25, kd_45, kd_47, kd_67, kd_99, kd_101, kd_133, kd_177, \
                         kd_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_25 * ab_x[k] * id_25[k]
                 - f_25 * ab_y[k] * id_27[k]
                 + f_26 * ab_y[k] * id_29[k]
                 + f_27 * ab_x[k] * id_67[k]
                 + f_27 * ab_y[k] * id_69[k]
                 - f_28 * ab_y[k] * id_71[k]
                 - f_29 * ab_x[k] * id_133[k]
                 - f_29 * ab_y[k] * id_135[k]
                 + f_30 * ab_y[k] * id_137[k]
                 - f_25 * kd_25[k]
                 - f_25 * kd_45[k]
                 + f_26 * kd_47[k]
                 + f_27 * kd_67[k]
                 + f_27 * kd_99[k]
                 - f_28 * kd_101[k]
                 - f_29 * kd_133[k]
                 - f_29 * kd_177[k]
                 + f_30 * kd_179[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, id_26, id_28, id_29, id_68, id_70, id_71, id_134, \
                         id_136, id_137, kd_26, kd_46, kd_53, kd_68, kd_100, kd_107, kd_134, \
                         kd_178, kd_185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_31 * ab_x[k] * id_26[k]
                  - f_31 * ab_y[k] * id_28[k]
                  + f_32 * ab_z[k] * id_29[k]
                  + f_33 * ab_x[k] * id_68[k]
                  + f_33 * ab_y[k] * id_70[k]
                  - f_34 * ab_z[k] * id_71[k]
                  - f_35 * ab_x[k] * id_134[k]
                  - f_35 * ab_y[k] * id_136[k]
                  + f_36 * ab_z[k] * id_137[k]
                  - f_31 * kd_26[k]
                  - f_31 * kd_46[k]
                  + f_32 * kd_53[k]
                  + f_33 * kd_68[k]
                  + f_33 * kd_100[k]
                  - f_34 * kd_107[k]
                  - f_35 * kd_134[k]
                  - f_35 * kd_178[k]
                  + f_36 * kd_185[k];
    }

#pragma omp simd aligned(ab_x, id_24, id_27, id_29, id_66, id_69, id_71, id_132, id_135, \
                         id_137, kd_24, kd_27, kd_29, kd_66, kd_69, kd_71, kd_132, kd_135, \
                         kd_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_25 * ab_x[k] * id_24[k]
                  - f_25 * ab_x[k] * id_27[k]
                  + f_26 * ab_x[k] * id_29[k]
                  + f_27 * ab_x[k] * id_66[k]
                  + f_27 * ab_x[k] * id_69[k]
                  - f_28 * ab_x[k] * id_71[k]
                  - f_29 * ab_x[k] * id_132[k]
                  - f_29 * ab_x[k] * id_135[k]
                  + f_30 * ab_x[k] * id_137[k]
                  - f_25 * kd_24[k]
                  - f_25 * kd_27[k]
                  + f_26 * kd_29[k]
                  + f_27 * kd_66[k]
                  + f_27 * kd_69[k]
                  - f_28 * kd_71[k]
                  - f_29 * kd_132[k]
                  - f_29 * kd_135[k]
                  + f_30 * kd_137[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_26, id_28, id_68, id_70, id_134, id_136, kd_26, kd_46, \
                         kd_68, kd_100, kd_134, kd_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_37 * ab_x[k] * id_26[k]
                  - f_37 * ab_y[k] * id_28[k]
                  - f_22 * ab_x[k] * id_68[k]
                  + f_22 * ab_y[k] * id_70[k]
                  + f_38 * ab_x[k] * id_134[k]
                  - f_38 * ab_y[k] * id_136[k]
                  + f_37 * kd_26[k]
                  - f_37 * kd_46[k]
                  - f_22 * kd_68[k]
                  + f_22 * kd_100[k]
                  + f_38 * kd_134[k]
                  - f_38 * kd_178[k];
    }

#pragma omp simd aligned(ab_x, id_24, id_27, id_66, id_69, id_132, id_135, kd_24, kd_27, \
                         kd_66, kd_69, kd_132, kd_135 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_17 * ab_x[k] * id_24[k]
                  - f_16 * ab_x[k] * id_27[k]
                  - f_19 * ab_x[k] * id_66[k]
                  + f_18 * ab_x[k] * id_69[k]
                  + f_21 * ab_x[k] * id_132[k]
                  - f_20 * ab_x[k] * id_135[k]
                  + f_17 * kd_24[k]
                  - f_16 * kd_27[k]
                  - f_19 * kd_66[k]
                  + f_18 * kd_69[k]
                  + f_21 * kd_132[k]
                  - f_20 * kd_135[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_7, id_9, id_49, id_51, id_91, id_93, id_103, id_105, \
                         kd_7, kd_21, kd_49, kd_75, kd_91, kd_103, kd_129, \
                         kd_141 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_39 * ab_x[k] * id_7[k]
                  + f_40 * ab_y[k] * id_9[k]
                  + f_41 * ab_x[k] * id_49[k]
                  - f_42 * ab_y[k] * id_51[k]
                  + f_39 * ab_x[k] * id_91[k]
                  - f_40 * ab_y[k] * id_93[k]
                  - f_41 * ab_x[k] * id_103[k]
                  + f_42 * ab_y[k] * id_105[k]
                  - f_39 * kd_7[k]
                  + f_40 * kd_21[k]
                  + f_41 * kd_49[k]
                  - f_42 * kd_75[k]
                  + f_39 * kd_91[k]
                  - f_41 * kd_103[k]
                  - f_40 * kd_129[k]
                  + f_42 * kd_141[k];
    }

#pragma omp simd aligned(ab_x, id_10, id_52, id_94, id_106, kd_10, kd_52, kd_94, \
                         kd_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_43 * ab_x[k] * id_10[k]
                  + f_44 * ab_x[k] * id_52[k]
                  + f_43 * ab_x[k] * id_94[k]
                  - f_44 * ab_x[k] * id_106[k]
                  - f_43 * kd_10[k]
                  + f_44 * kd_52[k]
                  + f_43 * kd_94[k]
                  - f_44 * kd_106[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_7, id_9, id_11, id_49, id_51, id_53, id_91, id_93, \
                         id_95, id_103, id_105, id_107, kd_7, kd_21, kd_23, kd_49, kd_75, \
                         kd_77, kd_91, kd_103, kd_129, kd_131, kd_141, \
                         kd_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_45 * ab_x[k] * id_7[k]
                  + f_45 * ab_y[k] * id_9[k]
                  - f_46 * ab_y[k] * id_11[k]
                  - f_47 * ab_x[k] * id_49[k]
                  - f_47 * ab_y[k] * id_51[k]
                  + f_48 * ab_y[k] * id_53[k]
                  - f_45 * ab_x[k] * id_91[k]
                  - f_45 * ab_y[k] * id_93[k]
                  + f_46 * ab_y[k] * id_95[k]
                  + f_47 * ab_x[k] * id_103[k]
                  + f_47 * ab_y[k] * id_105[k]
                  - f_48 * ab_y[k] * id_107[k]
                  + f_45 * kd_7[k]
                  + f_45 * kd_21[k]
                  - f_46 * kd_23[k]
                  - f_47 * kd_49[k]
                  - f_47 * kd_75[k]
                  + f_48 * kd_77[k]
                  - f_45 * kd_91[k]
                  + f_47 * kd_103[k]
                  - f_45 * kd_129[k]
                  + f_46 * kd_131[k]
                  + f_47 * kd_141[k]
                  - f_48 * kd_143[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, id_8, id_10, id_11, id_50, id_52, id_53, id_92, \
                         id_94, id_95, id_104, id_106, id_107, kd_8, kd_22, kd_29, kd_50, \
                         kd_76, kd_83, kd_92, kd_104, kd_130, kd_137, kd_142, \
                         kd_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_49 * ab_x[k] * id_8[k]
                  + f_49 * ab_y[k] * id_10[k]
                  - f_50 * ab_z[k] * id_11[k]
                  - f_51 * ab_x[k] * id_50[k]
                  - f_51 * ab_y[k] * id_52[k]
                  + f_52 * ab_z[k] * id_53[k]
                  - f_49 * ab_x[k] * id_92[k]
                  - f_49 * ab_y[k] * id_94[k]
                  + f_50 * ab_z[k] * id_95[k]
                  + f_51 * ab_x[k] * id_104[k]
                  + f_51 * ab_y[k] * id_106[k]
                  - f_52 * ab_z[k] * id_107[k]
                  + f_49 * kd_8[k]
                  + f_49 * kd_22[k]
                  - f_50 * kd_29[k]
                  - f_51 * kd_50[k]
                  - f_51 * kd_76[k]
                  + f_52 * kd_83[k]
                  - f_49 * kd_92[k]
                  + f_51 * kd_104[k]
                  - f_49 * kd_130[k]
                  + f_50 * kd_137[k]
                  + f_51 * kd_142[k]
                  - f_52 * kd_149[k];
    }

#pragma omp simd aligned(ab_x, id_6, id_9, id_11, id_48, id_51, id_53, id_90, id_93, id_95, \
                         id_102, id_105, id_107, kd_6, kd_9, kd_11, kd_48, kd_51, kd_53, \
                         kd_90, kd_93, kd_95, kd_102, kd_105, kd_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_45 * ab_x[k] * id_6[k]
                  + f_45 * ab_x[k] * id_9[k]
                  - f_46 * ab_x[k] * id_11[k]
                  - f_47 * ab_x[k] * id_48[k]
                  - f_47 * ab_x[k] * id_51[k]
                  + f_48 * ab_x[k] * id_53[k]
                  - f_45 * ab_x[k] * id_90[k]
                  - f_45 * ab_x[k] * id_93[k]
                  + f_46 * ab_x[k] * id_95[k]
                  + f_47 * ab_x[k] * id_102[k]
                  + f_47 * ab_x[k] * id_105[k]
                  - f_48 * ab_x[k] * id_107[k]
                  + f_45 * kd_6[k]
                  + f_45 * kd_9[k]
                  - f_46 * kd_11[k]
                  - f_47 * kd_48[k]
                  - f_47 * kd_51[k]
                  + f_48 * kd_53[k]
                  - f_45 * kd_90[k]
                  - f_45 * kd_93[k]
                  + f_46 * kd_95[k]
                  + f_47 * kd_102[k]
                  + f_47 * kd_105[k]
                  - f_48 * kd_107[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_8, id_10, id_50, id_52, id_92, id_94, id_104, id_106, \
                         kd_8, kd_22, kd_50, kd_76, kd_92, kd_104, kd_130, \
                         kd_142 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_53 * ab_x[k] * id_8[k]
                  + f_53 * ab_y[k] * id_10[k]
                  + f_54 * ab_x[k] * id_50[k]
                  - f_54 * ab_y[k] * id_52[k]
                  + f_53 * ab_x[k] * id_92[k]
                  - f_53 * ab_y[k] * id_94[k]
                  - f_54 * ab_x[k] * id_104[k]
                  + f_54 * ab_y[k] * id_106[k]
                  - f_53 * kd_8[k]
                  + f_53 * kd_22[k]
                  + f_54 * kd_50[k]
                  - f_54 * kd_76[k]
                  + f_53 * kd_92[k]
                  - f_54 * kd_104[k]
                  - f_53 * kd_130[k]
                  + f_54 * kd_142[k];
    }

#pragma omp simd aligned(ab_x, id_6, id_9, id_48, id_51, id_90, id_93, id_102, id_105, kd_6, \
                         kd_9, kd_48, kd_51, kd_90, kd_93, kd_102, \
                         kd_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_40 * ab_x[k] * id_6[k]
                  + f_39 * ab_x[k] * id_9[k]
                  + f_42 * ab_x[k] * id_48[k]
                  - f_41 * ab_x[k] * id_51[k]
                  + f_40 * ab_x[k] * id_90[k]
                  - f_39 * ab_x[k] * id_93[k]
                  - f_42 * ab_x[k] * id_102[k]
                  + f_41 * ab_x[k] * id_105[k]
                  - f_40 * kd_6[k]
                  + f_39 * kd_9[k]
                  + f_42 * kd_48[k]
                  - f_41 * kd_51[k]
                  + f_40 * kd_90[k]
                  - f_39 * kd_93[k]
                  - f_42 * kd_102[k]
                  + f_41 * kd_105[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_25, id_27, id_67, id_69, id_79, id_81, id_133, id_135, \
                         id_145, id_147, kd_25, kd_45, kd_67, kd_79, kd_99, kd_111, kd_133, \
                         kd_145, kd_177, kd_189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_55 * ab_x[k] * id_25[k]
                  + f_56 * ab_y[k] * id_27[k]
                  - f_57 * ab_x[k] * id_67[k]
                  + f_58 * ab_y[k] * id_69[k]
                  + f_59 * ab_x[k] * id_79[k]
                  - f_60 * ab_y[k] * id_81[k]
                  + f_56 * ab_x[k] * id_133[k]
                  - f_61 * ab_y[k] * id_135[k]
                  - f_60 * ab_x[k] * id_145[k]
                  + f_62 * ab_y[k] * id_147[k]
                  - f_55 * kd_25[k]
                  + f_56 * kd_45[k]
                  - f_57 * kd_67[k]
                  + f_59 * kd_79[k]
                  + f_58 * kd_99[k]
                  - f_60 * kd_111[k]
                  + f_56 * kd_133[k]
                  - f_60 * kd_145[k]
                  - f_61 * kd_177[k]
                  + f_62 * kd_189[k];
    }

#pragma omp simd aligned(ab_x, id_28, id_70, id_82, id_136, id_148, kd_28, kd_70, kd_82, \
                         kd_136, kd_148 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_63 * ab_x[k] * id_28[k]
                  - f_64 * ab_x[k] * id_70[k]
                  + f_65 * ab_x[k] * id_82[k]
                  + f_66 * ab_x[k] * id_136[k]
                  - f_67 * ab_x[k] * id_148[k]
                  - f_63 * kd_28[k]
                  - f_64 * kd_70[k]
                  + f_65 * kd_82[k]
                  + f_66 * kd_136[k]
                  - f_67 * kd_148[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_25, id_27, id_29, id_67, id_69, id_71, id_79, id_81, \
                         id_83, id_133, id_135, id_137, id_145, id_147, id_149, kd_25, kd_45, \
                         kd_47, kd_67, kd_79, kd_99, kd_101, kd_111, kd_113, kd_133, kd_145, \
                         kd_177, kd_179, kd_189, kd_191 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_68 * ab_x[k] * id_25[k]
                  + f_68 * ab_y[k] * id_27[k]
                  - f_69 * ab_y[k] * id_29[k]
                  + f_70 * ab_x[k] * id_67[k]
                  + f_70 * ab_y[k] * id_69[k]
                  - f_71 * ab_y[k] * id_71[k]
                  - f_71 * ab_x[k] * id_79[k]
                  - f_71 * ab_y[k] * id_81[k]
                  + f_72 * ab_y[k] * id_83[k]
                  - f_73 * ab_x[k] * id_133[k]
                  - f_73 * ab_y[k] * id_135[k]
                  + f_74 * ab_y[k] * id_137[k]
                  + f_75 * ab_x[k] * id_145[k]
                  + f_75 * ab_y[k] * id_147[k]
                  - f_76 * ab_y[k] * id_149[k]
                  + f_68 * kd_25[k]
                  + f_68 * kd_45[k]
                  - f_69 * kd_47[k]
                  + f_70 * kd_67[k]
                  - f_71 * kd_79[k]
                  + f_70 * kd_99[k]
                  - f_71 * kd_101[k]
                  - f_71 * kd_111[k]
                  + f_72 * kd_113[k]
                  - f_73 * kd_133[k]
                  + f_75 * kd_145[k]
                  - f_73 * kd_177[k]
                  + f_74 * kd_179[k]
                  + f_75 * kd_189[k]
                  - f_76 * kd_191[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, id_26, id_28, id_29, id_68, id_70, id_71, id_80, \
                         id_82, id_83, id_134, id_136, id_137, id_146, id_148, id_149, kd_26, \
                         kd_46, kd_53, kd_68, kd_80, kd_100, kd_107, kd_112, kd_119, kd_134, \
                         kd_146, kd_178, kd_185, kd_190, kd_197 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_77 * ab_x[k] * id_26[k]
                  + f_77 * ab_y[k] * id_28[k]
                  - f_78 * ab_z[k] * id_29[k]
                  + f_78 * ab_x[k] * id_68[k]
                  + f_78 * ab_y[k] * id_70[k]
                  - f_79 * ab_z[k] * id_71[k]
                  - f_80 * ab_x[k] * id_80[k]
                  - f_80 * ab_y[k] * id_82[k]
                  + f_81 * ab_z[k] * id_83[k]
                  - f_82 * ab_x[k] * id_134[k]
                  - f_82 * ab_y[k] * id_136[k]
                  + f_83 * ab_z[k] * id_137[k]
                  + f_84 * ab_x[k] * id_146[k]
                  + f_84 * ab_y[k] * id_148[k]
                  - f_85 * ab_z[k] * id_149[k]
                  + f_77 * kd_26[k]
                  + f_77 * kd_46[k]
                  - f_78 * kd_53[k]
                  + f_78 * kd_68[k]
                  - f_80 * kd_80[k]
                  + f_78 * kd_100[k]
                  - f_79 * kd_107[k]
                  - f_80 * kd_112[k]
                  + f_81 * kd_119[k]
                  - f_82 * kd_134[k]
                  + f_84 * kd_146[k]
                  - f_82 * kd_178[k]
                  + f_83 * kd_185[k]
                  + f_84 * kd_190[k]
                  - f_85 * kd_197[k];
    }

#pragma omp simd aligned(ab_x, id_24, id_27, id_29, id_66, id_69, id_71, id_78, id_81, id_83, \
                         id_132, id_135, id_137, id_144, id_147, id_149, kd_24, kd_27, kd_29, \
                         kd_66, kd_69, kd_71, kd_78, kd_81, kd_83, kd_132, kd_135, kd_137, \
                         kd_144, kd_147, kd_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_68 * ab_x[k] * id_24[k]
                  + f_68 * ab_x[k] * id_27[k]
                  - f_69 * ab_x[k] * id_29[k]
                  + f_70 * ab_x[k] * id_66[k]
                  + f_70 * ab_x[k] * id_69[k]
                  - f_71 * ab_x[k] * id_71[k]
                  - f_71 * ab_x[k] * id_78[k]
                  - f_71 * ab_x[k] * id_81[k]
                  + f_72 * ab_x[k] * id_83[k]
                  - f_73 * ab_x[k] * id_132[k]
                  - f_73 * ab_x[k] * id_135[k]
                  + f_74 * ab_x[k] * id_137[k]
                  + f_75 * ab_x[k] * id_144[k]
                  + f_75 * ab_x[k] * id_147[k]
                  - f_76 * ab_x[k] * id_149[k]
                  + f_68 * kd_24[k]
                  + f_68 * kd_27[k]
                  - f_69 * kd_29[k]
                  + f_70 * kd_66[k]
                  + f_70 * kd_69[k]
                  - f_71 * kd_71[k]
                  - f_71 * kd_78[k]
                  - f_71 * kd_81[k]
                  + f_72 * kd_83[k]
                  - f_73 * kd_132[k]
                  - f_73 * kd_135[k]
                  + f_74 * kd_137[k]
                  + f_75 * kd_144[k]
                  + f_75 * kd_147[k]
                  - f_76 * kd_149[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_26, id_28, id_68, id_70, id_80, id_82, id_134, id_136, \
                         id_146, id_148, kd_26, kd_46, kd_68, kd_80, kd_100, kd_112, kd_134, \
                         kd_146, kd_178, kd_190 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_86 * ab_x[k] * id_26[k]
                  + f_86 * ab_y[k] * id_28[k]
                  - f_66 * ab_x[k] * id_68[k]
                  + f_66 * ab_y[k] * id_70[k]
                  + f_87 * ab_x[k] * id_80[k]
                  - f_87 * ab_y[k] * id_82[k]
                  + f_88 * ab_x[k] * id_134[k]
                  - f_88 * ab_y[k] * id_136[k]
                  - f_89 * ab_x[k] * id_146[k]
                  + f_89 * ab_y[k] * id_148[k]
                  - f_86 * kd_26[k]
                  + f_86 * kd_46[k]
                  - f_66 * kd_68[k]
                  + f_87 * kd_80[k]
                  + f_66 * kd_100[k]
                  - f_87 * kd_112[k]
                  + f_88 * kd_134[k]
                  - f_89 * kd_146[k]
                  - f_88 * kd_178[k]
                  + f_89 * kd_190[k];
    }

#pragma omp simd aligned(ab_x, id_24, id_27, id_66, id_69, id_78, id_81, id_132, id_135, \
                         id_144, id_147, kd_24, kd_27, kd_66, kd_69, kd_78, kd_81, kd_132, \
                         kd_135, kd_144, kd_147 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_56 * ab_x[k] * id_24[k]
                  + f_55 * ab_x[k] * id_27[k]
                  - f_58 * ab_x[k] * id_66[k]
                  + f_57 * ab_x[k] * id_69[k]
                  + f_60 * ab_x[k] * id_78[k]
                  - f_59 * ab_x[k] * id_81[k]
                  + f_61 * ab_x[k] * id_132[k]
                  - f_56 * ab_x[k] * id_135[k]
                  - f_62 * ab_x[k] * id_144[k]
                  + f_60 * ab_x[k] * id_147[k]
                  - f_56 * kd_24[k]
                  + f_55 * kd_27[k]
                  - f_58 * kd_66[k]
                  + f_57 * kd_69[k]
                  + f_60 * kd_78[k]
                  - f_59 * kd_81[k]
                  + f_61 * kd_132[k]
                  - f_56 * kd_135[k]
                  - f_62 * kd_144[k]
                  + f_60 * kd_147[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_7, id_9, id_37, id_39, id_49, id_51, id_91, id_93, \
                         id_103, id_105, id_115, id_117, kd_7, kd_21, kd_37, kd_49, kd_63, \
                         kd_75, kd_91, kd_103, kd_115, kd_129, kd_141, \
                         kd_153 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_61 * ab_x[k] * id_7[k]
                  - f_90 * ab_y[k] * id_9[k]
                  + f_58 * ab_x[k] * id_37[k]
                  - f_91 * ab_y[k] * id_39[k]
                  - f_92 * ab_x[k] * id_49[k]
                  + f_93 * ab_y[k] * id_51[k]
                  + f_61 * ab_x[k] * id_91[k]
                  - f_90 * ab_y[k] * id_93[k]
                  - f_92 * ab_x[k] * id_103[k]
                  + f_93 * ab_y[k] * id_105[k]
                  + f_92 * ab_x[k] * id_115[k]
                  - f_93 * ab_y[k] * id_117[k]
                  + f_61 * kd_7[k]
                  - f_90 * kd_21[k]
                  + f_58 * kd_37[k]
                  - f_92 * kd_49[k]
                  - f_91 * kd_63[k]
                  + f_93 * kd_75[k]
                  + f_61 * kd_91[k]
                  - f_92 * kd_103[k]
                  + f_92 * kd_115[k]
                  - f_90 * kd_129[k]
                  + f_93 * kd_141[k]
                  - f_93 * kd_153[k];
    }

#pragma omp simd aligned(ab_x, id_10, id_40, id_52, id_94, id_106, id_118, kd_10, kd_40, \
                         kd_52, kd_94, kd_106, kd_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_94 * ab_x[k] * id_10[k]
                  + f_95 * ab_x[k] * id_40[k]
                  - f_96 * ab_x[k] * id_52[k]
                  + f_94 * ab_x[k] * id_94[k]
                  - f_96 * ab_x[k] * id_106[k]
                  + f_96 * ab_x[k] * id_118[k]
                  + f_94 * kd_10[k]
                  + f_95 * kd_40[k]
                  - f_96 * kd_52[k]
                  + f_94 * kd_94[k]
                  - f_96 * kd_106[k]
                  + f_96 * kd_118[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_7, id_9, id_11, id_37, id_39, id_41, id_49, id_51, \
                         id_53, id_91, id_93, id_95, id_103, id_105, id_107, id_115, id_117, \
                         id_119, kd_7, kd_21, kd_23, kd_37, kd_49, kd_63, kd_65, kd_75, kd_77, \
                         kd_91, kd_103, kd_115, kd_129, kd_131, kd_141, kd_143, kd_153, \
                         kd_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_97 * ab_x[k] * id_7[k]
                  - f_97 * ab_y[k] * id_9[k]
                  + f_98 * ab_y[k] * id_11[k]
                  - f_99 * ab_x[k] * id_37[k]
                  - f_99 * ab_y[k] * id_39[k]
                  + f_75 * ab_y[k] * id_41[k]
                  + f_100 * ab_x[k] * id_49[k]
                  + f_100 * ab_y[k] * id_51[k]
                  - f_101 * ab_y[k] * id_53[k]
                  - f_97 * ab_x[k] * id_91[k]
                  - f_97 * ab_y[k] * id_93[k]
                  + f_98 * ab_y[k] * id_95[k]
                  + f_100 * ab_x[k] * id_103[k]
                  + f_100 * ab_y[k] * id_105[k]
                  - f_101 * ab_y[k] * id_107[k]
                  - f_100 * ab_x[k] * id_115[k]
                  - f_100 * ab_y[k] * id_117[k]
                  + f_101 * ab_y[k] * id_119[k]
                  - f_97 * kd_7[k]
                  - f_97 * kd_21[k]
                  + f_98 * kd_23[k]
                  - f_99 * kd_37[k]
                  + f_100 * kd_49[k]
                  - f_99 * kd_63[k]
                  + f_75 * kd_65[k]
                  + f_100 * kd_75[k]
                  - f_101 * kd_77[k]
                  - f_97 * kd_91[k]
                  + f_100 * kd_103[k]
                  - f_100 * kd_115[k]
                  - f_97 * kd_129[k]
                  + f_98 * kd_131[k]
                  + f_100 * kd_141[k]
                  - f_101 * kd_143[k]
                  - f_100 * kd_153[k]
                  + f_101 * kd_155[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, id_8, id_10, id_11, id_38, id_40, id_41, id_50, \
                         id_52, id_53, id_92, id_94, id_95, id_104, id_106, id_107, id_116, \
                         id_118, id_119, kd_8, kd_22, kd_29, kd_38, kd_50, kd_64, kd_71, \
                         kd_76, kd_83, kd_92, kd_104, kd_116, kd_130, kd_137, kd_142, kd_149, \
                         kd_154, kd_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_102 * ab_x[k] * id_8[k]
                  - f_102 * ab_y[k] * id_10[k]
                  + f_103 * ab_z[k] * id_11[k]
                  - f_83 * ab_x[k] * id_38[k]
                  - f_83 * ab_y[k] * id_40[k]
                  + f_104 * ab_z[k] * id_41[k]
                  + f_81 * ab_x[k] * id_50[k]
                  + f_81 * ab_y[k] * id_52[k]
                  - f_105 * ab_z[k] * id_53[k]
                  - f_102 * ab_x[k] * id_92[k]
                  - f_102 * ab_y[k] * id_94[k]
                  + f_103 * ab_z[k] * id_95[k]
                  + f_81 * ab_x[k] * id_104[k]
                  + f_81 * ab_y[k] * id_106[k]
                  - f_105 * ab_z[k] * id_107[k]
                  - f_81 * ab_x[k] * id_116[k]
                  - f_81 * ab_y[k] * id_118[k]
                  + f_105 * ab_z[k] * id_119[k]
                  - f_102 * kd_8[k]
                  - f_102 * kd_22[k]
                  + f_103 * kd_29[k]
                  - f_83 * kd_38[k]
                  + f_81 * kd_50[k]
                  - f_83 * kd_64[k]
                  + f_104 * kd_71[k]
                  + f_81 * kd_76[k]
                  - f_105 * kd_83[k]
                  - f_102 * kd_92[k]
                  + f_81 * kd_104[k]
                  - f_81 * kd_116[k]
                  - f_102 * kd_130[k]
                  + f_103 * kd_137[k]
                  + f_81 * kd_142[k]
                  - f_105 * kd_149[k]
                  - f_81 * kd_154[k]
                  + f_105 * kd_161[k];
    }

#pragma omp simd aligned(ab_x, id_6, id_9, id_11, id_36, id_39, id_41, id_48, id_51, id_53, \
                         id_90, id_93, id_95, id_102, id_105, id_107, id_114, id_117, id_119, \
                         kd_6, kd_9, kd_11, kd_36, kd_39, kd_41, kd_48, kd_51, kd_53, kd_90, \
                         kd_93, kd_95, kd_102, kd_105, kd_107, kd_114, kd_117, \
                         kd_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_97 * ab_x[k] * id_6[k]
                  - f_97 * ab_x[k] * id_9[k]
                  + f_98 * ab_x[k] * id_11[k]
                  - f_99 * ab_x[k] * id_36[k]
                  - f_99 * ab_x[k] * id_39[k]
                  + f_75 * ab_x[k] * id_41[k]
                  + f_100 * ab_x[k] * id_48[k]
                  + f_100 * ab_x[k] * id_51[k]
                  - f_101 * ab_x[k] * id_53[k]
                  - f_97 * ab_x[k] * id_90[k]
                  - f_97 * ab_x[k] * id_93[k]
                  + f_98 * ab_x[k] * id_95[k]
                  + f_100 * ab_x[k] * id_102[k]
                  + f_100 * ab_x[k] * id_105[k]
                  - f_101 * ab_x[k] * id_107[k]
                  - f_100 * ab_x[k] * id_114[k]
                  - f_100 * ab_x[k] * id_117[k]
                  + f_101 * ab_x[k] * id_119[k]
                  - f_97 * kd_6[k]
                  - f_97 * kd_9[k]
                  + f_98 * kd_11[k]
                  - f_99 * kd_36[k]
                  - f_99 * kd_39[k]
                  + f_75 * kd_41[k]
                  + f_100 * kd_48[k]
                  + f_100 * kd_51[k]
                  - f_101 * kd_53[k]
                  - f_97 * kd_90[k]
                  - f_97 * kd_93[k]
                  + f_98 * kd_95[k]
                  + f_100 * kd_102[k]
                  + f_100 * kd_105[k]
                  - f_101 * kd_107[k]
                  - f_100 * kd_114[k]
                  - f_100 * kd_117[k]
                  + f_101 * kd_119[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_8, id_10, id_38, id_40, id_50, id_52, id_92, id_94, \
                         id_104, id_106, id_116, id_118, kd_8, kd_22, kd_38, kd_50, kd_64, \
                         kd_76, kd_92, kd_104, kd_116, kd_130, kd_142, \
                         kd_154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_106 * ab_x[k] * id_8[k]
                  - f_106 * ab_y[k] * id_10[k]
                  + f_94 * ab_x[k] * id_38[k]
                  - f_94 * ab_y[k] * id_40[k]
                  - f_67 * ab_x[k] * id_50[k]
                  + f_67 * ab_y[k] * id_52[k]
                  + f_106 * ab_x[k] * id_92[k]
                  - f_106 * ab_y[k] * id_94[k]
                  - f_67 * ab_x[k] * id_104[k]
                  + f_67 * ab_y[k] * id_106[k]
                  + f_67 * ab_x[k] * id_116[k]
                  - f_67 * ab_y[k] * id_118[k]
                  + f_106 * kd_8[k]
                  - f_106 * kd_22[k]
                  + f_94 * kd_38[k]
                  - f_67 * kd_50[k]
                  - f_94 * kd_64[k]
                  + f_67 * kd_76[k]
                  + f_106 * kd_92[k]
                  - f_67 * kd_104[k]
                  + f_67 * kd_116[k]
                  - f_106 * kd_130[k]
                  + f_67 * kd_142[k]
                  - f_67 * kd_154[k];
    }

#pragma omp simd aligned(ab_x, id_6, id_9, id_36, id_39, id_48, id_51, id_90, id_93, id_102, \
                         id_105, id_114, id_117, kd_6, kd_9, kd_36, kd_39, kd_48, kd_51, \
                         kd_90, kd_93, kd_102, kd_105, kd_114, kd_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_90 * ab_x[k] * id_6[k]
                  - f_61 * ab_x[k] * id_9[k]
                  + f_91 * ab_x[k] * id_36[k]
                  - f_58 * ab_x[k] * id_39[k]
                  - f_93 * ab_x[k] * id_48[k]
                  + f_92 * ab_x[k] * id_51[k]
                  + f_90 * ab_x[k] * id_90[k]
                  - f_61 * ab_x[k] * id_93[k]
                  - f_93 * ab_x[k] * id_102[k]
                  + f_92 * ab_x[k] * id_105[k]
                  + f_93 * ab_x[k] * id_114[k]
                  - f_92 * ab_x[k] * id_117[k]
                  + f_90 * kd_6[k]
                  - f_61 * kd_9[k]
                  + f_91 * kd_36[k]
                  - f_58 * kd_39[k]
                  - f_93 * kd_48[k]
                  + f_92 * kd_51[k]
                  + f_90 * kd_90[k]
                  - f_61 * kd_93[k]
                  - f_93 * kd_102[k]
                  + f_92 * kd_105[k]
                  + f_93 * kd_114[k]
                  - f_92 * kd_117[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_25, id_27, id_67, id_69, id_79, id_81, id_133, id_135, \
                         id_145, id_147, id_157, id_159, kd_25, kd_45, kd_67, kd_79, kd_99, \
                         kd_111, kd_133, kd_145, kd_157, kd_177, kd_189, \
                         kd_201 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_107 * ab_x[k] * id_25[k]
                  - f_108 * ab_y[k] * id_27[k]
                  + f_109 * ab_x[k] * id_67[k]
                  - f_110 * ab_y[k] * id_69[k]
                  - f_111 * ab_x[k] * id_79[k]
                  + f_112 * ab_y[k] * id_81[k]
                  + f_107 * ab_x[k] * id_133[k]
                  - f_108 * ab_y[k] * id_135[k]
                  - f_111 * ab_x[k] * id_145[k]
                  + f_112 * ab_y[k] * id_147[k]
                  + f_84 * ab_x[k] * id_157[k]
                  - f_113 * ab_y[k] * id_159[k]
                  + f_107 * kd_25[k]
                  - f_108 * kd_45[k]
                  + f_109 * kd_67[k]
                  - f_111 * kd_79[k]
                  - f_110 * kd_99[k]
                  + f_112 * kd_111[k]
                  + f_107 * kd_133[k]
                  - f_111 * kd_145[k]
                  + f_84 * kd_157[k]
                  - f_108 * kd_177[k]
                  + f_112 * kd_189[k]
                  - f_113 * kd_201[k];
    }

#pragma omp simd aligned(ab_x, id_28, id_70, id_82, id_136, id_148, id_160, kd_28, kd_70, \
                         kd_82, kd_136, kd_148, kd_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_114 * ab_x[k] * id_28[k]
                  + f_115 * ab_x[k] * id_70[k]
                  - f_116 * ab_x[k] * id_82[k]
                  + f_114 * ab_x[k] * id_136[k]
                  - f_116 * ab_x[k] * id_148[k]
                  + f_76 * ab_x[k] * id_160[k]
                  + f_114 * kd_28[k]
                  + f_115 * kd_70[k]
                  - f_116 * kd_82[k]
                  + f_114 * kd_136[k]
                  - f_116 * kd_148[k]
                  + f_76 * kd_160[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_25, id_27, id_29, id_67, id_69, id_71, id_79, id_81, \
                         id_83, id_133, id_135, id_137, id_145, id_147, id_149, id_157, \
                         id_159, id_161, kd_25, kd_45, kd_47, kd_67, kd_79, kd_99, kd_101, \
                         kd_111, kd_113, kd_133, kd_145, kd_157, kd_177, kd_179, kd_189, \
                         kd_191, kd_201, kd_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_106 * ab_x[k] * id_25[k]
                  - f_106 * ab_y[k] * id_27[k]
                  + f_95 * ab_y[k] * id_29[k]
                  - f_94 * ab_x[k] * id_67[k]
                  - f_94 * ab_y[k] * id_69[k]
                  + f_89 * ab_y[k] * id_71[k]
                  + f_95 * ab_x[k] * id_79[k]
                  + f_95 * ab_y[k] * id_81[k]
                  - f_67 * ab_y[k] * id_83[k]
                  - f_106 * ab_x[k] * id_133[k]
                  - f_106 * ab_y[k] * id_135[k]
                  + f_95 * ab_y[k] * id_137[k]
                  + f_95 * ab_x[k] * id_145[k]
                  + f_95 * ab_y[k] * id_147[k]
                  - f_67 * ab_y[k] * id_149[k]
                  - f_117 * ab_x[k] * id_157[k]
                  - f_117 * ab_y[k] * id_159[k]
                  + f_118 * ab_y[k] * id_161[k]
                  - f_106 * kd_25[k]
                  - f_106 * kd_45[k]
                  + f_95 * kd_47[k]
                  - f_94 * kd_67[k]
                  + f_95 * kd_79[k]
                  - f_94 * kd_99[k]
                  + f_89 * kd_101[k]
                  + f_95 * kd_111[k]
                  - f_67 * kd_113[k]
                  - f_106 * kd_133[k]
                  + f_95 * kd_145[k]
                  - f_117 * kd_157[k]
                  - f_106 * kd_177[k]
                  + f_95 * kd_179[k]
                  + f_95 * kd_189[k]
                  - f_67 * kd_191[k]
                  - f_117 * kd_201[k]
                  + f_118 * kd_203[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, id_26, id_28, id_29, id_68, id_70, id_71, id_80, \
                         id_82, id_83, id_134, id_136, id_137, id_146, id_148, id_149, id_158, \
                         id_160, id_161, kd_26, kd_46, kd_53, kd_68, kd_80, kd_100, kd_107, \
                         kd_112, kd_119, kd_134, kd_146, kd_158, kd_178, kd_185, kd_190, \
                         kd_197, kd_202, kd_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_58 * ab_x[k] * id_26[k]
                  - f_58 * ab_y[k] * id_28[k]
                  + f_119 * ab_z[k] * id_29[k]
                  - f_120 * ab_x[k] * id_68[k]
                  - f_120 * ab_y[k] * id_70[k]
                  + f_62 * ab_z[k] * id_71[k]
                  + f_60 * ab_x[k] * id_80[k]
                  + f_60 * ab_y[k] * id_82[k]
                  - f_93 * ab_z[k] * id_83[k]
                  - f_58 * ab_x[k] * id_134[k]
                  - f_58 * ab_y[k] * id_136[k]
                  + f_119 * ab_z[k] * id_137[k]
                  + f_60 * ab_x[k] * id_146[k]
                  + f_60 * ab_y[k] * id_148[k]
                  - f_93 * ab_z[k] * id_149[k]
                  - f_121 * ab_x[k] * id_158[k]
                  - f_121 * ab_y[k] * id_160[k]
                  + f_122 * ab_z[k] * id_161[k]
                  - f_58 * kd_26[k]
                  - f_58 * kd_46[k]
                  + f_119 * kd_53[k]
                  - f_120 * kd_68[k]
                  + f_60 * kd_80[k]
                  - f_120 * kd_100[k]
                  + f_62 * kd_107[k]
                  + f_60 * kd_112[k]
                  - f_93 * kd_119[k]
                  - f_58 * kd_134[k]
                  + f_60 * kd_146[k]
                  - f_121 * kd_158[k]
                  - f_58 * kd_178[k]
                  + f_119 * kd_185[k]
                  + f_60 * kd_190[k]
                  - f_93 * kd_197[k]
                  - f_121 * kd_202[k]
                  + f_122 * kd_209[k];
    }

#pragma omp simd aligned(ab_x, id_24, id_27, id_29, id_66, id_69, id_71, id_78, id_81, id_83, \
                         id_132, id_135, id_137, id_144, id_147, id_149, id_156, id_159, \
                         id_161, kd_24, kd_27, kd_29, kd_66, kd_69, kd_71, kd_78, kd_81, \
                         kd_83, kd_132, kd_135, kd_137, kd_144, kd_147, kd_149, kd_156, \
                         kd_159, kd_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_106 * ab_x[k] * id_24[k]
                  - f_106 * ab_x[k] * id_27[k]
                  + f_95 * ab_x[k] * id_29[k]
                  - f_94 * ab_x[k] * id_66[k]
                  - f_94 * ab_x[k] * id_69[k]
                  + f_89 * ab_x[k] * id_71[k]
                  + f_95 * ab_x[k] * id_78[k]
                  + f_95 * ab_x[k] * id_81[k]
                  - f_67 * ab_x[k] * id_83[k]
                  - f_106 * ab_x[k] * id_132[k]
                  - f_106 * ab_x[k] * id_135[k]
                  + f_95 * ab_x[k] * id_137[k]
                  + f_95 * ab_x[k] * id_144[k]
                  + f_95 * ab_x[k] * id_147[k]
                  - f_67 * ab_x[k] * id_149[k]
                  - f_117 * ab_x[k] * id_156[k]
                  - f_117 * ab_x[k] * id_159[k]
                  + f_118 * ab_x[k] * id_161[k]
                  - f_106 * kd_24[k]
                  - f_106 * kd_27[k]
                  + f_95 * kd_29[k]
                  - f_94 * kd_66[k]
                  - f_94 * kd_69[k]
                  + f_89 * kd_71[k]
                  + f_95 * kd_78[k]
                  + f_95 * kd_81[k]
                  - f_67 * kd_83[k]
                  - f_106 * kd_132[k]
                  - f_106 * kd_135[k]
                  + f_95 * kd_137[k]
                  + f_95 * kd_144[k]
                  + f_95 * kd_147[k]
                  - f_67 * kd_149[k]
                  - f_117 * kd_156[k]
                  - f_117 * kd_159[k]
                  + f_118 * kd_161[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_26, id_28, id_68, id_70, id_80, id_82, id_134, id_136, \
                         id_146, id_148, id_158, id_160, kd_26, kd_46, kd_68, kd_80, kd_100, \
                         kd_112, kd_134, kd_146, kd_158, kd_178, kd_190, \
                         kd_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_123 * ab_x[k] * id_26[k]
                  - f_123 * ab_y[k] * id_28[k]
                  + f_114 * ab_x[k] * id_68[k]
                  - f_114 * ab_y[k] * id_70[k]
                  - f_115 * ab_x[k] * id_80[k]
                  + f_115 * ab_y[k] * id_82[k]
                  + f_123 * ab_x[k] * id_134[k]
                  - f_123 * ab_y[k] * id_136[k]
                  - f_115 * ab_x[k] * id_146[k]
                  + f_115 * ab_y[k] * id_148[k]
                  + f_100 * ab_x[k] * id_158[k]
                  - f_100 * ab_y[k] * id_160[k]
                  + f_123 * kd_26[k]
                  - f_123 * kd_46[k]
                  + f_114 * kd_68[k]
                  - f_115 * kd_80[k]
                  - f_114 * kd_100[k]
                  + f_115 * kd_112[k]
                  + f_123 * kd_134[k]
                  - f_115 * kd_146[k]
                  + f_100 * kd_158[k]
                  - f_123 * kd_178[k]
                  + f_115 * kd_190[k]
                  - f_100 * kd_202[k];
    }

#pragma omp simd aligned(ab_x, id_24, id_27, id_66, id_69, id_78, id_81, id_132, id_135, \
                         id_144, id_147, id_156, id_159, kd_24, kd_27, kd_66, kd_69, kd_78, \
                         kd_81, kd_132, kd_135, kd_144, kd_147, kd_156, \
                         kd_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_108 * ab_x[k] * id_24[k]
                  - f_107 * ab_x[k] * id_27[k]
                  + f_110 * ab_x[k] * id_66[k]
                  - f_109 * ab_x[k] * id_69[k]
                  - f_112 * ab_x[k] * id_78[k]
                  + f_111 * ab_x[k] * id_81[k]
                  + f_108 * ab_x[k] * id_132[k]
                  - f_107 * ab_x[k] * id_135[k]
                  - f_112 * ab_x[k] * id_144[k]
                  + f_111 * ab_x[k] * id_147[k]
                  + f_113 * ab_x[k] * id_156[k]
                  - f_84 * ab_x[k] * id_159[k]
                  + f_108 * kd_24[k]
                  - f_107 * kd_27[k]
                  + f_110 * kd_66[k]
                  - f_109 * kd_69[k]
                  - f_112 * kd_78[k]
                  + f_111 * kd_81[k]
                  + f_108 * kd_132[k]
                  - f_107 * kd_135[k]
                  - f_112 * kd_144[k]
                  + f_111 * kd_147[k]
                  + f_113 * kd_156[k]
                  - f_84 * kd_159[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_1, id_3, id_19, id_21, id_31, id_33, id_61, id_63, \
                         id_73, id_75, id_85, id_87, id_127, id_129, id_139, id_141, id_151, \
                         id_153, id_163, id_165, kd_1, kd_9, kd_19, kd_31, kd_39, kd_51, \
                         kd_61, kd_73, kd_85, kd_93, kd_105, kd_117, kd_127, kd_139, kd_151, \
                         kd_163, kd_171, kd_183, kd_195, kd_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_124 * ab_x[k] * id_1[k]
                  + f_125 * ab_y[k] * id_3[k]
                  - f_126 * ab_x[k] * id_19[k]
                  + f_124 * ab_y[k] * id_21[k]
                  + f_127 * ab_x[k] * id_31[k]
                  - f_128 * ab_y[k] * id_33[k]
                  - f_126 * ab_x[k] * id_61[k]
                  + f_124 * ab_y[k] * id_63[k]
                  + f_129 * ab_x[k] * id_73[k]
                  - f_130 * ab_y[k] * id_75[k]
                  - f_131 * ab_x[k] * id_85[k]
                  + f_132 * ab_y[k] * id_87[k]
                  - f_124 * ab_x[k] * id_127[k]
                  + f_125 * ab_y[k] * id_129[k]
                  + f_127 * ab_x[k] * id_139[k]
                  - f_128 * ab_y[k] * id_141[k]
                  - f_131 * ab_x[k] * id_151[k]
                  + f_132 * ab_y[k] * id_153[k]
                  + f_133 * ab_x[k] * id_163[k]
                  - f_134 * ab_y[k] * id_165[k]
                  - f_124 * kd_1[k]
                  + f_125 * kd_9[k]
                  - f_126 * kd_19[k]
                  + f_127 * kd_31[k]
                  + f_124 * kd_39[k]
                  - f_128 * kd_51[k]
                  - f_126 * kd_61[k]
                  + f_129 * kd_73[k]
                  - f_131 * kd_85[k]
                  + f_124 * kd_93[k]
                  - f_130 * kd_105[k]
                  + f_132 * kd_117[k]
                  - f_124 * kd_127[k]
                  + f_127 * kd_139[k]
                  - f_131 * kd_151[k]
                  + f_133 * kd_163[k]
                  + f_125 * kd_171[k]
                  - f_128 * kd_183[k]
                  + f_132 * kd_195[k]
                  - f_134 * kd_207[k];
    }

#pragma omp simd aligned(ab_x, id_4, id_22, id_34, id_64, id_76, id_88, id_130, id_142, \
                         id_154, id_166, kd_4, kd_22, kd_34, kd_64, kd_76, kd_88, kd_130, \
                         kd_142, kd_154, kd_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_135 * ab_x[k] * id_4[k]
                  - f_136 * ab_x[k] * id_22[k]
                  + f_137 * ab_x[k] * id_34[k]
                  - f_136 * ab_x[k] * id_64[k]
                  + f_138 * ab_x[k] * id_76[k]
                  - f_139 * ab_x[k] * id_88[k]
                  - f_135 * ab_x[k] * id_130[k]
                  + f_137 * ab_x[k] * id_142[k]
                  - f_139 * ab_x[k] * id_154[k]
                  + f_140 * ab_x[k] * id_166[k]
                  - f_135 * kd_4[k]
                  - f_136 * kd_22[k]
                  + f_137 * kd_34[k]
                  - f_136 * kd_64[k]
                  + f_138 * kd_76[k]
                  - f_139 * kd_88[k]
                  - f_135 * kd_130[k]
                  + f_137 * kd_142[k]
                  - f_139 * kd_154[k]
                  + f_140 * kd_166[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_1, id_3, id_5, id_19, id_21, id_23, id_31, id_33, \
                         id_35, id_61, id_63, id_65, id_73, id_75, id_77, id_85, id_87, id_89, \
                         id_127, id_129, id_131, id_139, id_141, id_143, id_151, id_153, \
                         id_155, id_163, id_165, id_167, kd_1, kd_9, kd_11, kd_19, kd_31, \
                         kd_39, kd_41, kd_51, kd_53, kd_61, kd_73, kd_85, kd_93, kd_95, \
                         kd_105, kd_107, kd_117, kd_119, kd_127, kd_139, kd_151, kd_163, \
                         kd_171, kd_173, kd_183, kd_185, kd_195, kd_197, kd_207, \
                         kd_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_141 * ab_x[k] * id_1[k]
                  + f_141 * ab_y[k] * id_3[k]
                  - f_142 * ab_y[k] * id_5[k]
                  + f_143 * ab_x[k] * id_19[k]
                  + f_143 * ab_y[k] * id_21[k]
                  - f_144 * ab_y[k] * id_23[k]
                  - f_145 * ab_x[k] * id_31[k]
                  - f_145 * ab_y[k] * id_33[k]
                  + f_146 * ab_y[k] * id_35[k]
                  + f_143 * ab_x[k] * id_61[k]
                  + f_143 * ab_y[k] * id_63[k]
                  - f_144 * ab_y[k] * id_65[k]
                  - f_147 * ab_x[k] * id_73[k]
                  - f_147 * ab_y[k] * id_75[k]
                  + f_148 * ab_y[k] * id_77[k]
                  + f_149 * ab_x[k] * id_85[k]
                  + f_149 * ab_y[k] * id_87[k]
                  - f_150 * ab_y[k] * id_89[k]
                  + f_141 * ab_x[k] * id_127[k]
                  + f_141 * ab_y[k] * id_129[k]
                  - f_142 * ab_y[k] * id_131[k]
                  - f_145 * ab_x[k] * id_139[k]
                  - f_145 * ab_y[k] * id_141[k]
                  + f_146 * ab_y[k] * id_143[k]
                  + f_149 * ab_x[k] * id_151[k]
                  + f_149 * ab_y[k] * id_153[k]
                  - f_150 * ab_y[k] * id_155[k]
                  - f_151 * ab_x[k] * id_163[k]
                  - f_151 * ab_y[k] * id_165[k]
                  + f_152 * ab_y[k] * id_167[k]
                  + f_141 * kd_1[k]
                  + f_141 * kd_9[k]
                  - f_142 * kd_11[k]
                  + f_143 * kd_19[k]
                  - f_145 * kd_31[k]
                  + f_143 * kd_39[k]
                  - f_144 * kd_41[k]
                  - f_145 * kd_51[k]
                  + f_146 * kd_53[k]
                  + f_143 * kd_61[k]
                  - f_147 * kd_73[k]
                  + f_149 * kd_85[k]
                  + f_143 * kd_93[k]
                  - f_144 * kd_95[k]
                  - f_147 * kd_105[k]
                  + f_148 * kd_107[k]
                  + f_149 * kd_117[k]
                  - f_150 * kd_119[k]
                  + f_141 * kd_127[k]
                  - f_145 * kd_139[k]
                  + f_149 * kd_151[k]
                  - f_151 * kd_163[k]
                  + f_141 * kd_171[k]
                  - f_142 * kd_173[k]
                  - f_145 * kd_183[k]
                  + f_146 * kd_185[k]
                  + f_149 * kd_195[k]
                  - f_150 * kd_197[k]
                  - f_151 * kd_207[k]
                  + f_152 * kd_209[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, id_2, id_4, id_5, id_20, id_22, id_23, id_32, \
                         id_34, id_35, id_62, id_64, id_65, id_74, id_76, id_77, id_86, id_88, \
                         id_89, id_128, id_130, id_131, id_140, id_142, id_143, id_152, \
                         id_154, id_155, id_164, id_166, id_167, kd_2, kd_10, kd_17, kd_20, \
                         kd_32, kd_40, kd_47, kd_52, kd_59, kd_62, kd_74, kd_86, kd_94, \
                         kd_101, kd_106, kd_113, kd_118, kd_125, kd_128, kd_140, kd_152, \
                         kd_164, kd_172, kd_179, kd_184, kd_191, kd_196, kd_203, kd_208, \
                         kd_215 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = 0.46875 * ab_x[k] * id_2[k]
                  + 0.46875 * ab_y[k] * id_4[k]
                  - 0.3125 * ab_z[k] * id_5[k]
                  + 1.40625 * ab_x[k] * id_20[k]
                  + 1.40625 * ab_y[k] * id_22[k]
                  - 0.9375 * ab_z[k] * id_23[k]
                  - 8.4375 * ab_x[k] * id_32[k]
                  - 8.4375 * ab_y[k] * id_34[k]
                  + 5.625 * ab_z[k] * id_35[k]
                  + 1.40625 * ab_x[k] * id_62[k]
                  + 1.40625 * ab_y[k] * id_64[k]
                  - 0.9375 * ab_z[k] * id_65[k]
                  - 16.875 * ab_x[k] * id_74[k]
                  - 16.875 * ab_y[k] * id_76[k]
                  + 11.25 * ab_z[k] * id_77[k]
                  + 11.25 * ab_x[k] * id_86[k]
                  + 11.25 * ab_y[k] * id_88[k]
                  - 7.5 * ab_z[k] * id_89[k]
                  + 0.46875 * ab_x[k] * id_128[k]
                  + 0.46875 * ab_y[k] * id_130[k]
                  - 0.3125 * ab_z[k] * id_131[k]
                  - 8.4375 * ab_x[k] * id_140[k]
                  - 8.4375 * ab_y[k] * id_142[k]
                  + 5.625 * ab_z[k] * id_143[k]
                  + 11.25 * ab_x[k] * id_152[k]
                  + 11.25 * ab_y[k] * id_154[k]
                  - 7.5 * ab_z[k] * id_155[k]
                  - 1.5 * ab_x[k] * id_164[k]
                  - 1.5 * ab_y[k] * id_166[k]
                  + ab_z[k] * id_167[k]
                  + 0.46875 * kd_2[k]
                  + 0.46875 * kd_10[k]
                  - 0.3125 * kd_17[k]
                  + 1.40625 * kd_20[k]
                  - 8.4375 * kd_32[k]
                  + 1.40625 * kd_40[k]
                  - 0.9375 * kd_47[k]
                  - 8.4375 * kd_52[k]
                  + 5.625 * kd_59[k]
                  + 1.40625 * kd_62[k]
                  - 16.875 * kd_74[k]
                  + 11.25 * kd_86[k]
                  + 1.40625 * kd_94[k]
                  - 0.9375 * kd_101[k]
                  - 16.875 * kd_106[k]
                  + 11.25 * kd_113[k]
                  + 11.25 * kd_118[k]
                  - 7.5 * kd_125[k]
                  + 0.46875 * kd_128[k]
                  - 8.4375 * kd_140[k]
                  + 11.25 * kd_152[k]
                  - 1.5 * kd_164[k]
                  + 0.46875 * kd_172[k]
                  - 0.3125 * kd_179[k]
                  - 8.4375 * kd_184[k]
                  + 5.625 * kd_191[k]
                  + 11.25 * kd_196[k]
                  - 7.5 * kd_203[k]
                  - 1.5 * kd_208[k]
                  + kd_215[k];
    }

#pragma omp simd aligned(ab_x, id_0, id_3, id_5, id_18, id_21, id_23, id_30, id_33, id_35, \
                         id_60, id_63, id_65, id_72, id_75, id_77, id_84, id_87, id_89, \
                         id_126, id_129, id_131, id_138, id_141, id_143, id_150, id_153, \
                         id_155, id_162, id_165, id_167, kd_0, kd_3, kd_5, kd_18, kd_21, \
                         kd_23, kd_30, kd_33, kd_35, kd_60, kd_63, kd_65, kd_72, kd_75, kd_77, \
                         kd_84, kd_87, kd_89, kd_126, kd_129, kd_131, kd_138, kd_141, kd_143, \
                         kd_150, kd_153, kd_155, kd_162, kd_165, \
                         kd_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_141 * ab_x[k] * id_0[k]
                  + f_141 * ab_x[k] * id_3[k]
                  - f_142 * ab_x[k] * id_5[k]
                  + f_143 * ab_x[k] * id_18[k]
                  + f_143 * ab_x[k] * id_21[k]
                  - f_144 * ab_x[k] * id_23[k]
                  - f_145 * ab_x[k] * id_30[k]
                  - f_145 * ab_x[k] * id_33[k]
                  + f_146 * ab_x[k] * id_35[k]
                  + f_143 * ab_x[k] * id_60[k]
                  + f_143 * ab_x[k] * id_63[k]
                  - f_144 * ab_x[k] * id_65[k]
                  - f_147 * ab_x[k] * id_72[k]
                  - f_147 * ab_x[k] * id_75[k]
                  + f_148 * ab_x[k] * id_77[k]
                  + f_149 * ab_x[k] * id_84[k]
                  + f_149 * ab_x[k] * id_87[k]
                  - f_150 * ab_x[k] * id_89[k]
                  + f_141 * ab_x[k] * id_126[k]
                  + f_141 * ab_x[k] * id_129[k]
                  - f_142 * ab_x[k] * id_131[k]
                  - f_145 * ab_x[k] * id_138[k]
                  - f_145 * ab_x[k] * id_141[k]
                  + f_146 * ab_x[k] * id_143[k]
                  + f_149 * ab_x[k] * id_150[k]
                  + f_149 * ab_x[k] * id_153[k]
                  - f_150 * ab_x[k] * id_155[k]
                  - f_151 * ab_x[k] * id_162[k]
                  - f_151 * ab_x[k] * id_165[k]
                  + f_152 * ab_x[k] * id_167[k]
                  + f_141 * kd_0[k]
                  + f_141 * kd_3[k]
                  - f_142 * kd_5[k]
                  + f_143 * kd_18[k]
                  + f_143 * kd_21[k]
                  - f_144 * kd_23[k]
                  - f_145 * kd_30[k]
                  - f_145 * kd_33[k]
                  + f_146 * kd_35[k]
                  + f_143 * kd_60[k]
                  + f_143 * kd_63[k]
                  - f_144 * kd_65[k]
                  - f_147 * kd_72[k]
                  - f_147 * kd_75[k]
                  + f_148 * kd_77[k]
                  + f_149 * kd_84[k]
                  + f_149 * kd_87[k]
                  - f_150 * kd_89[k]
                  + f_141 * kd_126[k]
                  + f_141 * kd_129[k]
                  - f_142 * kd_131[k]
                  - f_145 * kd_138[k]
                  - f_145 * kd_141[k]
                  + f_146 * kd_143[k]
                  + f_149 * kd_150[k]
                  + f_149 * kd_153[k]
                  - f_150 * kd_155[k]
                  - f_151 * kd_162[k]
                  - f_151 * kd_165[k]
                  + f_152 * kd_167[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_2, id_4, id_20, id_22, id_32, id_34, id_62, id_64, \
                         id_74, id_76, id_86, id_88, id_128, id_130, id_140, id_142, id_152, \
                         id_154, id_164, id_166, kd_2, kd_10, kd_20, kd_32, kd_40, kd_52, \
                         kd_62, kd_74, kd_86, kd_94, kd_106, kd_118, kd_128, kd_140, kd_152, \
                         kd_164, kd_172, kd_184, kd_196, kd_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_153 * ab_x[k] * id_2[k]
                  + f_153 * ab_y[k] * id_4[k]
                  - f_154 * ab_x[k] * id_20[k]
                  + f_154 * ab_y[k] * id_22[k]
                  + f_155 * ab_x[k] * id_32[k]
                  - f_155 * ab_y[k] * id_34[k]
                  - f_154 * ab_x[k] * id_62[k]
                  + f_154 * ab_y[k] * id_64[k]
                  + f_137 * ab_x[k] * id_74[k]
                  - f_137 * ab_y[k] * id_76[k]
                  - f_156 * ab_x[k] * id_86[k]
                  + f_156 * ab_y[k] * id_88[k]
                  - f_153 * ab_x[k] * id_128[k]
                  + f_153 * ab_y[k] * id_130[k]
                  + f_155 * ab_x[k] * id_140[k]
                  - f_155 * ab_y[k] * id_142[k]
                  - f_156 * ab_x[k] * id_152[k]
                  + f_156 * ab_y[k] * id_154[k]
                  + f_157 * ab_x[k] * id_164[k]
                  - f_157 * ab_y[k] * id_166[k]
                  - f_153 * kd_2[k]
                  + f_153 * kd_10[k]
                  - f_154 * kd_20[k]
                  + f_155 * kd_32[k]
                  + f_154 * kd_40[k]
                  - f_155 * kd_52[k]
                  - f_154 * kd_62[k]
                  + f_137 * kd_74[k]
                  - f_156 * kd_86[k]
                  + f_154 * kd_94[k]
                  - f_137 * kd_106[k]
                  + f_156 * kd_118[k]
                  - f_153 * kd_128[k]
                  + f_155 * kd_140[k]
                  - f_156 * kd_152[k]
                  + f_157 * kd_164[k]
                  + f_153 * kd_172[k]
                  - f_155 * kd_184[k]
                  + f_156 * kd_196[k]
                  - f_157 * kd_208[k];
    }

#pragma omp simd aligned(ab_x, id_0, id_3, id_18, id_21, id_30, id_33, id_60, id_63, id_72, \
                         id_75, id_84, id_87, id_126, id_129, id_138, id_141, id_150, id_153, \
                         id_162, id_165, kd_0, kd_3, kd_18, kd_21, kd_30, kd_33, kd_60, kd_63, \
                         kd_72, kd_75, kd_84, kd_87, kd_126, kd_129, kd_138, kd_141, kd_150, \
                         kd_153, kd_162, kd_165 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_125 * ab_x[k] * id_0[k]
                  + f_124 * ab_x[k] * id_3[k]
                  - f_124 * ab_x[k] * id_18[k]
                  + f_126 * ab_x[k] * id_21[k]
                  + f_128 * ab_x[k] * id_30[k]
                  - f_127 * ab_x[k] * id_33[k]
                  - f_124 * ab_x[k] * id_60[k]
                  + f_126 * ab_x[k] * id_63[k]
                  + f_130 * ab_x[k] * id_72[k]
                  - f_129 * ab_x[k] * id_75[k]
                  - f_132 * ab_x[k] * id_84[k]
                  + f_131 * ab_x[k] * id_87[k]
                  - f_125 * ab_x[k] * id_126[k]
                  + f_124 * ab_x[k] * id_129[k]
                  + f_128 * ab_x[k] * id_138[k]
                  - f_127 * ab_x[k] * id_141[k]
                  - f_132 * ab_x[k] * id_150[k]
                  + f_131 * ab_x[k] * id_153[k]
                  + f_134 * ab_x[k] * id_162[k]
                  - f_133 * ab_x[k] * id_165[k]
                  - f_125 * kd_0[k]
                  + f_124 * kd_3[k]
                  - f_124 * kd_18[k]
                  + f_126 * kd_21[k]
                  + f_128 * kd_30[k]
                  - f_127 * kd_33[k]
                  - f_124 * kd_60[k]
                  + f_126 * kd_63[k]
                  + f_130 * kd_72[k]
                  - f_129 * kd_75[k]
                  - f_132 * kd_84[k]
                  + f_131 * kd_87[k]
                  - f_125 * kd_126[k]
                  + f_124 * kd_129[k]
                  + f_128 * kd_138[k]
                  - f_127 * kd_141[k]
                  - f_132 * kd_150[k]
                  + f_131 * kd_153[k]
                  + f_134 * kd_162[k]
                  - f_133 * kd_165[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_13, id_15, id_43, id_45, id_55, id_57, id_97, id_99, \
                         id_109, id_111, id_121, id_123, kd_13, kd_27, kd_43, kd_55, kd_69, \
                         kd_81, kd_97, kd_109, kd_121, kd_135, kd_147, \
                         kd_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_107 * ab_x[k] * id_13[k]
                  - f_108 * ab_y[k] * id_15[k]
                  + f_109 * ab_x[k] * id_43[k]
                  - f_110 * ab_y[k] * id_45[k]
                  - f_111 * ab_x[k] * id_55[k]
                  + f_112 * ab_y[k] * id_57[k]
                  + f_107 * ab_x[k] * id_97[k]
                  - f_108 * ab_y[k] * id_99[k]
                  - f_111 * ab_x[k] * id_109[k]
                  + f_112 * ab_y[k] * id_111[k]
                  + f_84 * ab_x[k] * id_121[k]
                  - f_113 * ab_y[k] * id_123[k]
                  + f_107 * kd_13[k]
                  - f_108 * kd_27[k]
                  + f_109 * kd_43[k]
                  - f_111 * kd_55[k]
                  - f_110 * kd_69[k]
                  + f_112 * kd_81[k]
                  + f_107 * kd_97[k]
                  - f_111 * kd_109[k]
                  + f_84 * kd_121[k]
                  - f_108 * kd_135[k]
                  + f_112 * kd_147[k]
                  - f_113 * kd_159[k];
    }

#pragma omp simd aligned(ab_x, id_16, id_46, id_58, id_100, id_112, id_124, kd_16, kd_46, \
                         kd_58, kd_100, kd_112, kd_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = f_114 * ab_x[k] * id_16[k]
                  + f_115 * ab_x[k] * id_46[k]
                  - f_116 * ab_x[k] * id_58[k]
                  + f_114 * ab_x[k] * id_100[k]
                  - f_116 * ab_x[k] * id_112[k]
                  + f_76 * ab_x[k] * id_124[k]
                  + f_114 * kd_16[k]
                  + f_115 * kd_46[k]
                  - f_116 * kd_58[k]
                  + f_114 * kd_100[k]
                  - f_116 * kd_112[k]
                  + f_76 * kd_124[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_13, id_15, id_17, id_43, id_45, id_47, id_55, id_57, \
                         id_59, id_97, id_99, id_101, id_109, id_111, id_113, id_121, id_123, \
                         id_125, kd_13, kd_27, kd_29, kd_43, kd_55, kd_69, kd_71, kd_81, \
                         kd_83, kd_97, kd_109, kd_121, kd_135, kd_137, kd_147, kd_149, kd_159, \
                         kd_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_106 * ab_x[k] * id_13[k]
                  - f_106 * ab_y[k] * id_15[k]
                  + f_95 * ab_y[k] * id_17[k]
                  - f_94 * ab_x[k] * id_43[k]
                  - f_94 * ab_y[k] * id_45[k]
                  + f_89 * ab_y[k] * id_47[k]
                  + f_95 * ab_x[k] * id_55[k]
                  + f_95 * ab_y[k] * id_57[k]
                  - f_67 * ab_y[k] * id_59[k]
                  - f_106 * ab_x[k] * id_97[k]
                  - f_106 * ab_y[k] * id_99[k]
                  + f_95 * ab_y[k] * id_101[k]
                  + f_95 * ab_x[k] * id_109[k]
                  + f_95 * ab_y[k] * id_111[k]
                  - f_67 * ab_y[k] * id_113[k]
                  - f_117 * ab_x[k] * id_121[k]
                  - f_117 * ab_y[k] * id_123[k]
                  + f_118 * ab_y[k] * id_125[k]
                  - f_106 * kd_13[k]
                  - f_106 * kd_27[k]
                  + f_95 * kd_29[k]
                  - f_94 * kd_43[k]
                  + f_95 * kd_55[k]
                  - f_94 * kd_69[k]
                  + f_89 * kd_71[k]
                  + f_95 * kd_81[k]
                  - f_67 * kd_83[k]
                  - f_106 * kd_97[k]
                  + f_95 * kd_109[k]
                  - f_117 * kd_121[k]
                  - f_106 * kd_135[k]
                  + f_95 * kd_137[k]
                  + f_95 * kd_147[k]
                  - f_67 * kd_149[k]
                  - f_117 * kd_159[k]
                  + f_118 * kd_161[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, id_14, id_16, id_17, id_44, id_46, id_47, id_56, \
                         id_58, id_59, id_98, id_100, id_101, id_110, id_112, id_113, id_122, \
                         id_124, id_125, kd_14, kd_28, kd_35, kd_44, kd_56, kd_70, kd_77, \
                         kd_82, kd_89, kd_98, kd_110, kd_122, kd_136, kd_143, kd_148, kd_155, \
                         kd_160, kd_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_58 * ab_x[k] * id_14[k]
                  - f_58 * ab_y[k] * id_16[k]
                  + f_119 * ab_z[k] * id_17[k]
                  - f_120 * ab_x[k] * id_44[k]
                  - f_120 * ab_y[k] * id_46[k]
                  + f_62 * ab_z[k] * id_47[k]
                  + f_60 * ab_x[k] * id_56[k]
                  + f_60 * ab_y[k] * id_58[k]
                  - f_93 * ab_z[k] * id_59[k]
                  - f_58 * ab_x[k] * id_98[k]
                  - f_58 * ab_y[k] * id_100[k]
                  + f_119 * ab_z[k] * id_101[k]
                  + f_60 * ab_x[k] * id_110[k]
                  + f_60 * ab_y[k] * id_112[k]
                  - f_93 * ab_z[k] * id_113[k]
                  - f_121 * ab_x[k] * id_122[k]
                  - f_121 * ab_y[k] * id_124[k]
                  + f_122 * ab_z[k] * id_125[k]
                  - f_58 * kd_14[k]
                  - f_58 * kd_28[k]
                  + f_119 * kd_35[k]
                  - f_120 * kd_44[k]
                  + f_60 * kd_56[k]
                  - f_120 * kd_70[k]
                  + f_62 * kd_77[k]
                  + f_60 * kd_82[k]
                  - f_93 * kd_89[k]
                  - f_58 * kd_98[k]
                  + f_60 * kd_110[k]
                  - f_121 * kd_122[k]
                  - f_58 * kd_136[k]
                  + f_119 * kd_143[k]
                  + f_60 * kd_148[k]
                  - f_93 * kd_155[k]
                  - f_121 * kd_160[k]
                  + f_122 * kd_167[k];
    }

#pragma omp simd aligned(ab_x, id_12, id_15, id_17, id_42, id_45, id_47, id_54, id_57, id_59, \
                         id_96, id_99, id_101, id_108, id_111, id_113, id_120, id_123, id_125, \
                         kd_12, kd_15, kd_17, kd_42, kd_45, kd_47, kd_54, kd_57, kd_59, kd_96, \
                         kd_99, kd_101, kd_108, kd_111, kd_113, kd_120, kd_123, \
                         kd_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = -f_106 * ab_x[k] * id_12[k]
                  - f_106 * ab_x[k] * id_15[k]
                  + f_95 * ab_x[k] * id_17[k]
                  - f_94 * ab_x[k] * id_42[k]
                  - f_94 * ab_x[k] * id_45[k]
                  + f_89 * ab_x[k] * id_47[k]
                  + f_95 * ab_x[k] * id_54[k]
                  + f_95 * ab_x[k] * id_57[k]
                  - f_67 * ab_x[k] * id_59[k]
                  - f_106 * ab_x[k] * id_96[k]
                  - f_106 * ab_x[k] * id_99[k]
                  + f_95 * ab_x[k] * id_101[k]
                  + f_95 * ab_x[k] * id_108[k]
                  + f_95 * ab_x[k] * id_111[k]
                  - f_67 * ab_x[k] * id_113[k]
                  - f_117 * ab_x[k] * id_120[k]
                  - f_117 * ab_x[k] * id_123[k]
                  + f_118 * ab_x[k] * id_125[k]
                  - f_106 * kd_12[k]
                  - f_106 * kd_15[k]
                  + f_95 * kd_17[k]
                  - f_94 * kd_42[k]
                  - f_94 * kd_45[k]
                  + f_89 * kd_47[k]
                  + f_95 * kd_54[k]
                  + f_95 * kd_57[k]
                  - f_67 * kd_59[k]
                  - f_106 * kd_96[k]
                  - f_106 * kd_99[k]
                  + f_95 * kd_101[k]
                  + f_95 * kd_108[k]
                  + f_95 * kd_111[k]
                  - f_67 * kd_113[k]
                  - f_117 * kd_120[k]
                  - f_117 * kd_123[k]
                  + f_118 * kd_125[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_14, id_16, id_44, id_46, id_56, id_58, id_98, id_100, \
                         id_110, id_112, id_122, id_124, kd_14, kd_28, kd_44, kd_56, kd_70, \
                         kd_82, kd_98, kd_110, kd_122, kd_136, kd_148, \
                         kd_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_123 * ab_x[k] * id_14[k]
                  - f_123 * ab_y[k] * id_16[k]
                  + f_114 * ab_x[k] * id_44[k]
                  - f_114 * ab_y[k] * id_46[k]
                  - f_115 * ab_x[k] * id_56[k]
                  + f_115 * ab_y[k] * id_58[k]
                  + f_123 * ab_x[k] * id_98[k]
                  - f_123 * ab_y[k] * id_100[k]
                  - f_115 * ab_x[k] * id_110[k]
                  + f_115 * ab_y[k] * id_112[k]
                  + f_100 * ab_x[k] * id_122[k]
                  - f_100 * ab_y[k] * id_124[k]
                  + f_123 * kd_14[k]
                  - f_123 * kd_28[k]
                  + f_114 * kd_44[k]
                  - f_115 * kd_56[k]
                  - f_114 * kd_70[k]
                  + f_115 * kd_82[k]
                  + f_123 * kd_98[k]
                  - f_115 * kd_110[k]
                  + f_100 * kd_122[k]
                  - f_123 * kd_136[k]
                  + f_115 * kd_148[k]
                  - f_100 * kd_160[k];
    }

#pragma omp simd aligned(ab_x, id_12, id_15, id_42, id_45, id_54, id_57, id_96, id_99, id_108, \
                         id_111, id_120, id_123, kd_12, kd_15, kd_42, kd_45, kd_54, kd_57, \
                         kd_96, kd_99, kd_108, kd_111, kd_120, kd_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_108 * ab_x[k] * id_12[k]
                  - f_107 * ab_x[k] * id_15[k]
                  + f_110 * ab_x[k] * id_42[k]
                  - f_109 * ab_x[k] * id_45[k]
                  - f_112 * ab_x[k] * id_54[k]
                  + f_111 * ab_x[k] * id_57[k]
                  + f_108 * ab_x[k] * id_96[k]
                  - f_107 * ab_x[k] * id_99[k]
                  - f_112 * ab_x[k] * id_108[k]
                  + f_111 * ab_x[k] * id_111[k]
                  + f_113 * ab_x[k] * id_120[k]
                  - f_84 * ab_x[k] * id_123[k]
                  + f_108 * kd_12[k]
                  - f_107 * kd_15[k]
                  + f_110 * kd_42[k]
                  - f_109 * kd_45[k]
                  - f_112 * kd_54[k]
                  + f_111 * kd_57[k]
                  + f_108 * kd_96[k]
                  - f_107 * kd_99[k]
                  - f_112 * kd_108[k]
                  + f_111 * kd_111[k]
                  + f_113 * kd_120[k]
                  - f_84 * kd_123[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_1, id_3, id_19, id_21, id_31, id_33, id_61, id_63, \
                         id_85, id_87, id_127, id_129, id_139, id_141, id_151, id_153, kd_1, \
                         kd_9, kd_19, kd_31, kd_39, kd_51, kd_61, kd_85, kd_93, kd_117, \
                         kd_127, kd_139, kd_151, kd_171, kd_183, \
                         kd_195 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_158 * ab_x[k] * id_1[k]
                  - f_159 * ab_y[k] * id_3[k]
                  + f_158 * ab_x[k] * id_19[k]
                  - f_159 * ab_y[k] * id_21[k]
                  - f_60 * ab_x[k] * id_31[k]
                  + f_62 * ab_y[k] * id_33[k]
                  - f_158 * ab_x[k] * id_61[k]
                  + f_159 * ab_y[k] * id_63[k]
                  + f_60 * ab_x[k] * id_85[k]
                  - f_62 * ab_y[k] * id_87[k]
                  - f_158 * ab_x[k] * id_127[k]
                  + f_159 * ab_y[k] * id_129[k]
                  + f_60 * ab_x[k] * id_139[k]
                  - f_62 * ab_y[k] * id_141[k]
                  - f_60 * ab_x[k] * id_151[k]
                  + f_62 * ab_y[k] * id_153[k]
                  + f_158 * kd_1[k]
                  - f_159 * kd_9[k]
                  + f_158 * kd_19[k]
                  - f_60 * kd_31[k]
                  - f_159 * kd_39[k]
                  + f_62 * kd_51[k]
                  - f_158 * kd_61[k]
                  + f_60 * kd_85[k]
                  + f_159 * kd_93[k]
                  - f_62 * kd_117[k]
                  - f_158 * kd_127[k]
                  + f_60 * kd_139[k]
                  - f_60 * kd_151[k]
                  + f_159 * kd_171[k]
                  - f_62 * kd_183[k]
                  + f_62 * kd_195[k];
    }

#pragma omp simd aligned(ab_x, id_4, id_22, id_34, id_64, id_88, id_130, id_142, id_154, kd_4, \
                         kd_22, kd_34, kd_64, kd_88, kd_130, kd_142, \
                         kd_154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_106 * ab_x[k] * id_4[k]
                  + f_106 * ab_x[k] * id_22[k]
                  - f_67 * ab_x[k] * id_34[k]
                  - f_106 * ab_x[k] * id_64[k]
                  + f_67 * ab_x[k] * id_88[k]
                  - f_106 * ab_x[k] * id_130[k]
                  + f_67 * ab_x[k] * id_142[k]
                  - f_67 * ab_x[k] * id_154[k]
                  + f_106 * kd_4[k]
                  + f_106 * kd_22[k]
                  - f_67 * kd_34[k]
                  - f_106 * kd_64[k]
                  + f_67 * kd_88[k]
                  - f_106 * kd_130[k]
                  + f_67 * kd_142[k]
                  - f_67 * kd_154[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_1, id_3, id_5, id_19, id_21, id_23, id_31, id_33, \
                         id_35, id_61, id_63, id_65, id_85, id_87, id_89, id_127, id_129, \
                         id_131, id_139, id_141, id_143, id_151, id_153, id_155, kd_1, kd_9, \
                         kd_11, kd_19, kd_31, kd_39, kd_41, kd_51, kd_53, kd_61, kd_85, kd_93, \
                         kd_95, kd_117, kd_119, kd_127, kd_139, kd_151, kd_171, kd_173, \
                         kd_183, kd_185, kd_195, kd_197 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_160 * ab_x[k] * id_1[k]
                  - f_160 * ab_y[k] * id_3[k]
                  + f_99 * ab_y[k] * id_5[k]
                  - f_160 * ab_x[k] * id_19[k]
                  - f_160 * ab_y[k] * id_21[k]
                  + f_99 * ab_y[k] * id_23[k]
                  + f_75 * ab_x[k] * id_31[k]
                  + f_75 * ab_y[k] * id_33[k]
                  - f_76 * ab_y[k] * id_35[k]
                  + f_160 * ab_x[k] * id_61[k]
                  + f_160 * ab_y[k] * id_63[k]
                  - f_99 * ab_y[k] * id_65[k]
                  - f_75 * ab_x[k] * id_85[k]
                  - f_75 * ab_y[k] * id_87[k]
                  + f_76 * ab_y[k] * id_89[k]
                  + f_160 * ab_x[k] * id_127[k]
                  + f_160 * ab_y[k] * id_129[k]
                  - f_99 * ab_y[k] * id_131[k]
                  - f_75 * ab_x[k] * id_139[k]
                  - f_75 * ab_y[k] * id_141[k]
                  + f_76 * ab_y[k] * id_143[k]
                  + f_75 * ab_x[k] * id_151[k]
                  + f_75 * ab_y[k] * id_153[k]
                  - f_76 * ab_y[k] * id_155[k]
                  - f_160 * kd_1[k]
                  - f_160 * kd_9[k]
                  + f_99 * kd_11[k]
                  - f_160 * kd_19[k]
                  + f_75 * kd_31[k]
                  - f_160 * kd_39[k]
                  + f_99 * kd_41[k]
                  + f_75 * kd_51[k]
                  - f_76 * kd_53[k]
                  + f_160 * kd_61[k]
                  - f_75 * kd_85[k]
                  + f_160 * kd_93[k]
                  - f_99 * kd_95[k]
                  - f_75 * kd_117[k]
                  + f_76 * kd_119[k]
                  + f_160 * kd_127[k]
                  - f_75 * kd_139[k]
                  + f_75 * kd_151[k]
                  + f_160 * kd_171[k]
                  - f_99 * kd_173[k]
                  - f_75 * kd_183[k]
                  + f_76 * kd_185[k]
                  + f_75 * kd_195[k]
                  - f_76 * kd_197[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, id_2, id_4, id_5, id_20, id_22, id_23, id_32, \
                         id_34, id_35, id_62, id_64, id_65, id_86, id_88, id_89, id_128, \
                         id_130, id_131, id_140, id_142, id_143, id_152, id_154, id_155, kd_2, \
                         kd_10, kd_17, kd_20, kd_32, kd_40, kd_47, kd_52, kd_59, kd_62, kd_86, \
                         kd_94, kd_101, kd_118, kd_125, kd_128, kd_140, kd_152, kd_172, \
                         kd_179, kd_184, kd_191, kd_196, kd_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_161 * ab_x[k] * id_2[k]
                  - f_161 * ab_y[k] * id_4[k]
                  + f_162 * ab_z[k] * id_5[k]
                  - f_161 * ab_x[k] * id_20[k]
                  - f_161 * ab_y[k] * id_22[k]
                  + f_162 * ab_z[k] * id_23[k]
                  + f_84 * ab_x[k] * id_32[k]
                  + f_84 * ab_y[k] * id_34[k]
                  - f_85 * ab_z[k] * id_35[k]
                  + f_161 * ab_x[k] * id_62[k]
                  + f_161 * ab_y[k] * id_64[k]
                  - f_162 * ab_z[k] * id_65[k]
                  - f_84 * ab_x[k] * id_86[k]
                  - f_84 * ab_y[k] * id_88[k]
                  + f_85 * ab_z[k] * id_89[k]
                  + f_161 * ab_x[k] * id_128[k]
                  + f_161 * ab_y[k] * id_130[k]
                  - f_162 * ab_z[k] * id_131[k]
                  - f_84 * ab_x[k] * id_140[k]
                  - f_84 * ab_y[k] * id_142[k]
                  + f_85 * ab_z[k] * id_143[k]
                  + f_84 * ab_x[k] * id_152[k]
                  + f_84 * ab_y[k] * id_154[k]
                  - f_85 * ab_z[k] * id_155[k]
                  - f_161 * kd_2[k]
                  - f_161 * kd_10[k]
                  + f_162 * kd_17[k]
                  - f_161 * kd_20[k]
                  + f_84 * kd_32[k]
                  - f_161 * kd_40[k]
                  + f_162 * kd_47[k]
                  + f_84 * kd_52[k]
                  - f_85 * kd_59[k]
                  + f_161 * kd_62[k]
                  - f_84 * kd_86[k]
                  + f_161 * kd_94[k]
                  - f_162 * kd_101[k]
                  - f_84 * kd_118[k]
                  + f_85 * kd_125[k]
                  + f_161 * kd_128[k]
                  - f_84 * kd_140[k]
                  + f_84 * kd_152[k]
                  + f_161 * kd_172[k]
                  - f_162 * kd_179[k]
                  - f_84 * kd_184[k]
                  + f_85 * kd_191[k]
                  + f_84 * kd_196[k]
                  - f_85 * kd_203[k];
    }

#pragma omp simd aligned(ab_x, id_0, id_3, id_5, id_18, id_21, id_23, id_30, id_33, id_35, \
                         id_60, id_63, id_65, id_84, id_87, id_89, id_126, id_129, id_131, \
                         id_138, id_141, id_143, id_150, id_153, id_155, kd_0, kd_3, kd_5, \
                         kd_18, kd_21, kd_23, kd_30, kd_33, kd_35, kd_60, kd_63, kd_65, kd_84, \
                         kd_87, kd_89, kd_126, kd_129, kd_131, kd_138, kd_141, kd_143, kd_150, \
                         kd_153, kd_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -f_160 * ab_x[k] * id_0[k]
                  - f_160 * ab_x[k] * id_3[k]
                  + f_99 * ab_x[k] * id_5[k]
                  - f_160 * ab_x[k] * id_18[k]
                  - f_160 * ab_x[k] * id_21[k]
                  + f_99 * ab_x[k] * id_23[k]
                  + f_75 * ab_x[k] * id_30[k]
                  + f_75 * ab_x[k] * id_33[k]
                  - f_76 * ab_x[k] * id_35[k]
                  + f_160 * ab_x[k] * id_60[k]
                  + f_160 * ab_x[k] * id_63[k]
                  - f_99 * ab_x[k] * id_65[k]
                  - f_75 * ab_x[k] * id_84[k]
                  - f_75 * ab_x[k] * id_87[k]
                  + f_76 * ab_x[k] * id_89[k]
                  + f_160 * ab_x[k] * id_126[k]
                  + f_160 * ab_x[k] * id_129[k]
                  - f_99 * ab_x[k] * id_131[k]
                  - f_75 * ab_x[k] * id_138[k]
                  - f_75 * ab_x[k] * id_141[k]
                  + f_76 * ab_x[k] * id_143[k]
                  + f_75 * ab_x[k] * id_150[k]
                  + f_75 * ab_x[k] * id_153[k]
                  - f_76 * ab_x[k] * id_155[k]
                  - f_160 * kd_0[k]
                  - f_160 * kd_3[k]
                  + f_99 * kd_5[k]
                  - f_160 * kd_18[k]
                  - f_160 * kd_21[k]
                  + f_99 * kd_23[k]
                  + f_75 * kd_30[k]
                  + f_75 * kd_33[k]
                  - f_76 * kd_35[k]
                  + f_160 * kd_60[k]
                  + f_160 * kd_63[k]
                  - f_99 * kd_65[k]
                  - f_75 * kd_84[k]
                  - f_75 * kd_87[k]
                  + f_76 * kd_89[k]
                  + f_160 * kd_126[k]
                  + f_160 * kd_129[k]
                  - f_99 * kd_131[k]
                  - f_75 * kd_138[k]
                  - f_75 * kd_141[k]
                  + f_76 * kd_143[k]
                  + f_75 * kd_150[k]
                  + f_75 * kd_153[k]
                  - f_76 * kd_155[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_2, id_4, id_20, id_22, id_32, id_34, id_62, id_64, \
                         id_86, id_88, id_128, id_130, id_140, id_142, id_152, id_154, kd_2, \
                         kd_10, kd_20, kd_32, kd_40, kd_52, kd_62, kd_86, kd_94, kd_118, \
                         kd_128, kd_140, kd_152, kd_172, kd_184, \
                         kd_196 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_163 * ab_x[k] * id_2[k]
                  - f_163 * ab_y[k] * id_4[k]
                  + f_163 * ab_x[k] * id_20[k]
                  - f_163 * ab_y[k] * id_22[k]
                  - f_89 * ab_x[k] * id_32[k]
                  + f_89 * ab_y[k] * id_34[k]
                  - f_163 * ab_x[k] * id_62[k]
                  + f_163 * ab_y[k] * id_64[k]
                  + f_89 * ab_x[k] * id_86[k]
                  - f_89 * ab_y[k] * id_88[k]
                  - f_163 * ab_x[k] * id_128[k]
                  + f_163 * ab_y[k] * id_130[k]
                  + f_89 * ab_x[k] * id_140[k]
                  - f_89 * ab_y[k] * id_142[k]
                  - f_89 * ab_x[k] * id_152[k]
                  + f_89 * ab_y[k] * id_154[k]
                  + f_163 * kd_2[k]
                  - f_163 * kd_10[k]
                  + f_163 * kd_20[k]
                  - f_89 * kd_32[k]
                  - f_163 * kd_40[k]
                  + f_89 * kd_52[k]
                  - f_163 * kd_62[k]
                  + f_89 * kd_86[k]
                  + f_163 * kd_94[k]
                  - f_89 * kd_118[k]
                  - f_163 * kd_128[k]
                  + f_89 * kd_140[k]
                  - f_89 * kd_152[k]
                  + f_163 * kd_172[k]
                  - f_89 * kd_184[k]
                  + f_89 * kd_196[k];
    }

#pragma omp simd aligned(ab_x, id_0, id_3, id_18, id_21, id_30, id_33, id_60, id_63, id_84, \
                         id_87, id_126, id_129, id_138, id_141, id_150, id_153, kd_0, kd_3, \
                         kd_18, kd_21, kd_30, kd_33, kd_60, kd_63, kd_84, kd_87, kd_126, \
                         kd_129, kd_138, kd_141, kd_150, kd_153 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_159 * ab_x[k] * id_0[k]
                  - f_158 * ab_x[k] * id_3[k]
                  + f_159 * ab_x[k] * id_18[k]
                  - f_158 * ab_x[k] * id_21[k]
                  - f_62 * ab_x[k] * id_30[k]
                  + f_60 * ab_x[k] * id_33[k]
                  - f_159 * ab_x[k] * id_60[k]
                  + f_158 * ab_x[k] * id_63[k]
                  + f_62 * ab_x[k] * id_84[k]
                  - f_60 * ab_x[k] * id_87[k]
                  - f_159 * ab_x[k] * id_126[k]
                  + f_158 * ab_x[k] * id_129[k]
                  + f_62 * ab_x[k] * id_138[k]
                  - f_60 * ab_x[k] * id_141[k]
                  - f_62 * ab_x[k] * id_150[k]
                  + f_60 * ab_x[k] * id_153[k]
                  + f_159 * kd_0[k]
                  - f_158 * kd_3[k]
                  + f_159 * kd_18[k]
                  - f_158 * kd_21[k]
                  - f_62 * kd_30[k]
                  + f_60 * kd_33[k]
                  - f_159 * kd_60[k]
                  + f_158 * kd_63[k]
                  + f_62 * kd_84[k]
                  - f_60 * kd_87[k]
                  - f_159 * kd_126[k]
                  + f_158 * kd_129[k]
                  + f_62 * kd_138[k]
                  - f_60 * kd_141[k]
                  - f_62 * kd_150[k]
                  + f_60 * kd_153[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_13, id_15, id_43, id_45, id_55, id_57, id_97, id_99, \
                         id_109, id_111, kd_13, kd_27, kd_43, kd_55, kd_69, kd_81, kd_97, \
                         kd_109, kd_135, kd_147 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_56 * ab_x[k] * id_13[k]
                  + f_61 * ab_y[k] * id_15[k]
                  + f_57 * ab_x[k] * id_43[k]
                  - f_58 * ab_y[k] * id_45[k]
                  + f_60 * ab_x[k] * id_55[k]
                  - f_62 * ab_y[k] * id_57[k]
                  + f_55 * ab_x[k] * id_97[k]
                  - f_56 * ab_y[k] * id_99[k]
                  - f_59 * ab_x[k] * id_109[k]
                  + f_60 * ab_y[k] * id_111[k]
                  - f_56 * kd_13[k]
                  + f_61 * kd_27[k]
                  + f_57 * kd_43[k]
                  + f_60 * kd_55[k]
                  - f_58 * kd_69[k]
                  - f_62 * kd_81[k]
                  + f_55 * kd_97[k]
                  - f_59 * kd_109[k]
                  - f_56 * kd_135[k]
                  + f_60 * kd_147[k];
    }

#pragma omp simd aligned(ab_x, id_16, id_46, id_58, id_100, id_112, kd_16, kd_46, kd_58, \
                         kd_100, kd_112 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_66 * ab_x[k] * id_16[k]
                  + f_64 * ab_x[k] * id_46[k]
                  + f_67 * ab_x[k] * id_58[k]
                  + f_63 * ab_x[k] * id_100[k]
                  - f_65 * ab_x[k] * id_112[k]
                  - f_66 * kd_16[k]
                  + f_64 * kd_46[k]
                  + f_67 * kd_58[k]
                  + f_63 * kd_100[k]
                  - f_65 * kd_112[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_13, id_15, id_17, id_43, id_45, id_47, id_55, id_57, \
                         id_59, id_97, id_99, id_101, id_109, id_111, id_113, kd_13, kd_27, \
                         kd_29, kd_43, kd_55, kd_69, kd_71, kd_81, kd_83, kd_97, kd_109, \
                         kd_135, kd_137, kd_147, kd_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_73 * ab_x[k] * id_13[k]
                  + f_73 * ab_y[k] * id_15[k]
                  - f_74 * ab_y[k] * id_17[k]
                  - f_70 * ab_x[k] * id_43[k]
                  - f_70 * ab_y[k] * id_45[k]
                  + f_71 * ab_y[k] * id_47[k]
                  - f_75 * ab_x[k] * id_55[k]
                  - f_75 * ab_y[k] * id_57[k]
                  + f_76 * ab_y[k] * id_59[k]
                  - f_68 * ab_x[k] * id_97[k]
                  - f_68 * ab_y[k] * id_99[k]
                  + f_69 * ab_y[k] * id_101[k]
                  + f_71 * ab_x[k] * id_109[k]
                  + f_71 * ab_y[k] * id_111[k]
                  - f_72 * ab_y[k] * id_113[k]
                  + f_73 * kd_13[k]
                  + f_73 * kd_27[k]
                  - f_74 * kd_29[k]
                  - f_70 * kd_43[k]
                  - f_75 * kd_55[k]
                  - f_70 * kd_69[k]
                  + f_71 * kd_71[k]
                  - f_75 * kd_81[k]
                  + f_76 * kd_83[k]
                  - f_68 * kd_97[k]
                  + f_71 * kd_109[k]
                  - f_68 * kd_135[k]
                  + f_69 * kd_137[k]
                  + f_71 * kd_147[k]
                  - f_72 * kd_149[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, id_14, id_16, id_17, id_44, id_46, id_47, id_56, \
                         id_58, id_59, id_98, id_100, id_101, id_110, id_112, id_113, kd_14, \
                         kd_28, kd_35, kd_44, kd_56, kd_70, kd_77, kd_82, kd_89, kd_98, \
                         kd_110, kd_136, kd_143, kd_148, kd_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = f_82 * ab_x[k] * id_14[k]
                  + f_82 * ab_y[k] * id_16[k]
                  - f_83 * ab_z[k] * id_17[k]
                  - f_78 * ab_x[k] * id_44[k]
                  - f_78 * ab_y[k] * id_46[k]
                  + f_79 * ab_z[k] * id_47[k]
                  - f_84 * ab_x[k] * id_56[k]
                  - f_84 * ab_y[k] * id_58[k]
                  + f_85 * ab_z[k] * id_59[k]
                  - f_77 * ab_x[k] * id_98[k]
                  - f_77 * ab_y[k] * id_100[k]
                  + f_78 * ab_z[k] * id_101[k]
                  + f_80 * ab_x[k] * id_110[k]
                  + f_80 * ab_y[k] * id_112[k]
                  - f_81 * ab_z[k] * id_113[k]
                  + f_82 * kd_14[k]
                  + f_82 * kd_28[k]
                  - f_83 * kd_35[k]
                  - f_78 * kd_44[k]
                  - f_84 * kd_56[k]
                  - f_78 * kd_70[k]
                  + f_79 * kd_77[k]
                  - f_84 * kd_82[k]
                  + f_85 * kd_89[k]
                  - f_77 * kd_98[k]
                  + f_80 * kd_110[k]
                  - f_77 * kd_136[k]
                  + f_78 * kd_143[k]
                  + f_80 * kd_148[k]
                  - f_81 * kd_155[k];
    }

#pragma omp simd aligned(ab_x, id_12, id_15, id_17, id_42, id_45, id_47, id_54, id_57, id_59, \
                         id_96, id_99, id_101, id_108, id_111, id_113, kd_12, kd_15, kd_17, \
                         kd_42, kd_45, kd_47, kd_54, kd_57, kd_59, kd_96, kd_99, kd_101, \
                         kd_108, kd_111, kd_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = f_73 * ab_x[k] * id_12[k]
                  + f_73 * ab_x[k] * id_15[k]
                  - f_74 * ab_x[k] * id_17[k]
                  - f_70 * ab_x[k] * id_42[k]
                  - f_70 * ab_x[k] * id_45[k]
                  + f_71 * ab_x[k] * id_47[k]
                  - f_75 * ab_x[k] * id_54[k]
                  - f_75 * ab_x[k] * id_57[k]
                  + f_76 * ab_x[k] * id_59[k]
                  - f_68 * ab_x[k] * id_96[k]
                  - f_68 * ab_x[k] * id_99[k]
                  + f_69 * ab_x[k] * id_101[k]
                  + f_71 * ab_x[k] * id_108[k]
                  + f_71 * ab_x[k] * id_111[k]
                  - f_72 * ab_x[k] * id_113[k]
                  + f_73 * kd_12[k]
                  + f_73 * kd_15[k]
                  - f_74 * kd_17[k]
                  - f_70 * kd_42[k]
                  - f_70 * kd_45[k]
                  + f_71 * kd_47[k]
                  - f_75 * kd_54[k]
                  - f_75 * kd_57[k]
                  + f_76 * kd_59[k]
                  - f_68 * kd_96[k]
                  - f_68 * kd_99[k]
                  + f_69 * kd_101[k]
                  + f_71 * kd_108[k]
                  + f_71 * kd_111[k]
                  - f_72 * kd_113[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_14, id_16, id_44, id_46, id_56, id_58, id_98, id_100, \
                         id_110, id_112, kd_14, kd_28, kd_44, kd_56, kd_70, kd_82, kd_98, \
                         kd_110, kd_136, kd_148 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = -f_88 * ab_x[k] * id_14[k]
                  + f_88 * ab_y[k] * id_16[k]
                  + f_66 * ab_x[k] * id_44[k]
                  - f_66 * ab_y[k] * id_46[k]
                  + f_89 * ab_x[k] * id_56[k]
                  - f_89 * ab_y[k] * id_58[k]
                  + f_86 * ab_x[k] * id_98[k]
                  - f_86 * ab_y[k] * id_100[k]
                  - f_87 * ab_x[k] * id_110[k]
                  + f_87 * ab_y[k] * id_112[k]
                  - f_88 * kd_14[k]
                  + f_88 * kd_28[k]
                  + f_66 * kd_44[k]
                  + f_89 * kd_56[k]
                  - f_66 * kd_70[k]
                  - f_89 * kd_82[k]
                  + f_86 * kd_98[k]
                  - f_87 * kd_110[k]
                  - f_86 * kd_136[k]
                  + f_87 * kd_148[k];
    }

#pragma omp simd aligned(ab_x, id_12, id_15, id_42, id_45, id_54, id_57, id_96, id_99, id_108, \
                         id_111, kd_12, kd_15, kd_42, kd_45, kd_54, kd_57, kd_96, kd_99, \
                         kd_108, kd_111 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = -f_61 * ab_x[k] * id_12[k]
                  + f_56 * ab_x[k] * id_15[k]
                  + f_58 * ab_x[k] * id_42[k]
                  - f_57 * ab_x[k] * id_45[k]
                  + f_62 * ab_x[k] * id_54[k]
                  - f_60 * ab_x[k] * id_57[k]
                  + f_56 * ab_x[k] * id_96[k]
                  - f_55 * ab_x[k] * id_99[k]
                  - f_60 * ab_x[k] * id_108[k]
                  + f_59 * ab_x[k] * id_111[k]
                  - f_61 * kd_12[k]
                  + f_56 * kd_15[k]
                  + f_58 * kd_42[k]
                  - f_57 * kd_45[k]
                  + f_62 * kd_54[k]
                  - f_60 * kd_57[k]
                  + f_56 * kd_96[k]
                  - f_55 * kd_99[k]
                  - f_60 * kd_108[k]
                  + f_59 * kd_111[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_1, id_3, id_19, id_21, id_31, id_33, id_61, id_63, \
                         id_73, id_75, id_127, id_129, id_139, id_141, kd_1, kd_9, kd_19, \
                         kd_31, kd_39, kd_51, kd_61, kd_73, kd_93, kd_105, kd_127, kd_139, \
                         kd_171, kd_183 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_164 * ab_x[k] * id_1[k]
                  + f_165 * ab_y[k] * id_3[k]
                  + f_166 * ab_x[k] * id_19[k]
                  - f_167 * ab_y[k] * id_21[k]
                  + f_168 * ab_x[k] * id_31[k]
                  - f_169 * ab_y[k] * id_33[k]
                  + f_166 * ab_x[k] * id_61[k]
                  - f_167 * ab_y[k] * id_63[k]
                  - f_170 * ab_x[k] * id_73[k]
                  + f_171 * ab_y[k] * id_75[k]
                  - f_164 * ab_x[k] * id_127[k]
                  + f_165 * ab_y[k] * id_129[k]
                  + f_168 * ab_x[k] * id_139[k]
                  - f_169 * ab_y[k] * id_141[k]
                  - f_164 * kd_1[k]
                  + f_165 * kd_9[k]
                  + f_166 * kd_19[k]
                  + f_168 * kd_31[k]
                  - f_167 * kd_39[k]
                  - f_169 * kd_51[k]
                  + f_166 * kd_61[k]
                  - f_170 * kd_73[k]
                  - f_167 * kd_93[k]
                  + f_171 * kd_105[k]
                  - f_164 * kd_127[k]
                  + f_168 * kd_139[k]
                  + f_165 * kd_171[k]
                  - f_169 * kd_183[k];
    }

#pragma omp simd aligned(ab_x, id_4, id_22, id_34, id_64, id_76, id_130, id_142, kd_4, kd_22, \
                         kd_34, kd_64, kd_76, kd_130, kd_142 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_172 * ab_x[k] * id_4[k]
                  + f_173 * ab_x[k] * id_22[k]
                  + f_174 * ab_x[k] * id_34[k]
                  + f_173 * ab_x[k] * id_64[k]
                  - f_175 * ab_x[k] * id_76[k]
                  - f_172 * ab_x[k] * id_130[k]
                  + f_174 * ab_x[k] * id_142[k]
                  - f_172 * kd_4[k]
                  + f_173 * kd_22[k]
                  + f_174 * kd_34[k]
                  + f_173 * kd_64[k]
                  - f_175 * kd_76[k]
                  - f_172 * kd_130[k]
                  + f_174 * kd_142[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_1, id_3, id_5, id_19, id_21, id_23, id_31, id_33, \
                         id_35, id_61, id_63, id_65, id_73, id_75, id_77, id_127, id_129, \
                         id_131, id_139, id_141, id_143, kd_1, kd_9, kd_11, kd_19, kd_31, \
                         kd_39, kd_41, kd_51, kd_53, kd_61, kd_73, kd_93, kd_95, kd_105, \
                         kd_107, kd_127, kd_139, kd_171, kd_173, kd_183, \
                         kd_185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = f_176 * ab_x[k] * id_1[k]
                  + f_176 * ab_y[k] * id_3[k]
                  - f_45 * ab_y[k] * id_5[k]
                  - f_177 * ab_x[k] * id_19[k]
                  - f_177 * ab_y[k] * id_21[k]
                  + f_178 * ab_y[k] * id_23[k]
                  - f_179 * ab_x[k] * id_31[k]
                  - f_179 * ab_y[k] * id_33[k]
                  + f_47 * ab_y[k] * id_35[k]
                  - f_177 * ab_x[k] * id_61[k]
                  - f_177 * ab_y[k] * id_63[k]
                  + f_178 * ab_y[k] * id_65[k]
                  + f_180 * ab_x[k] * id_73[k]
                  + f_180 * ab_y[k] * id_75[k]
                  - f_181 * ab_y[k] * id_77[k]
                  + f_176 * ab_x[k] * id_127[k]
                  + f_176 * ab_y[k] * id_129[k]
                  - f_45 * ab_y[k] * id_131[k]
                  - f_179 * ab_x[k] * id_139[k]
                  - f_179 * ab_y[k] * id_141[k]
                  + f_47 * ab_y[k] * id_143[k]
                  + f_176 * kd_1[k]
                  + f_176 * kd_9[k]
                  - f_45 * kd_11[k]
                  - f_177 * kd_19[k]
                  - f_179 * kd_31[k]
                  - f_177 * kd_39[k]
                  + f_178 * kd_41[k]
                  - f_179 * kd_51[k]
                  + f_47 * kd_53[k]
                  - f_177 * kd_61[k]
                  + f_180 * kd_73[k]
                  - f_177 * kd_93[k]
                  + f_178 * kd_95[k]
                  + f_180 * kd_105[k]
                  - f_181 * kd_107[k]
                  + f_176 * kd_127[k]
                  - f_179 * kd_139[k]
                  + f_176 * kd_171[k]
                  - f_45 * kd_173[k]
                  - f_179 * kd_183[k]
                  + f_47 * kd_185[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, id_2, id_4, id_5, id_20, id_22, id_23, id_32, \
                         id_34, id_35, id_62, id_64, id_65, id_74, id_76, id_77, id_128, \
                         id_130, id_131, id_140, id_142, id_143, kd_2, kd_10, kd_17, kd_20, \
                         kd_32, kd_40, kd_47, kd_52, kd_59, kd_62, kd_74, kd_94, kd_101, \
                         kd_106, kd_113, kd_128, kd_140, kd_172, kd_179, kd_184, \
                         kd_191 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_182 * ab_x[k] * id_2[k]
                  + f_182 * ab_y[k] * id_4[k]
                  - f_183 * ab_z[k] * id_5[k]
                  - f_184 * ab_x[k] * id_20[k]
                  - f_184 * ab_y[k] * id_22[k]
                  + f_185 * ab_z[k] * id_23[k]
                  - f_186 * ab_x[k] * id_32[k]
                  - f_186 * ab_y[k] * id_34[k]
                  + f_187 * ab_z[k] * id_35[k]
                  - f_184 * ab_x[k] * id_62[k]
                  - f_184 * ab_y[k] * id_64[k]
                  + f_185 * ab_z[k] * id_65[k]
                  + f_188 * ab_x[k] * id_74[k]
                  + f_188 * ab_y[k] * id_76[k]
                  - f_51 * ab_z[k] * id_77[k]
                  + f_182 * ab_x[k] * id_128[k]
                  + f_182 * ab_y[k] * id_130[k]
                  - f_183 * ab_z[k] * id_131[k]
                  - f_186 * ab_x[k] * id_140[k]
                  - f_186 * ab_y[k] * id_142[k]
                  + f_187 * ab_z[k] * id_143[k]
                  + f_182 * kd_2[k]
                  + f_182 * kd_10[k]
                  - f_183 * kd_17[k]
                  - f_184 * kd_20[k]
                  - f_186 * kd_32[k]
                  - f_184 * kd_40[k]
                  + f_185 * kd_47[k]
                  - f_186 * kd_52[k]
                  + f_187 * kd_59[k]
                  - f_184 * kd_62[k]
                  + f_188 * kd_74[k]
                  - f_184 * kd_94[k]
                  + f_185 * kd_101[k]
                  + f_188 * kd_106[k]
                  - f_51 * kd_113[k]
                  + f_182 * kd_128[k]
                  - f_186 * kd_140[k]
                  + f_182 * kd_172[k]
                  - f_183 * kd_179[k]
                  - f_186 * kd_184[k]
                  + f_187 * kd_191[k];
    }

#pragma omp simd aligned(ab_x, id_0, id_3, id_5, id_18, id_21, id_23, id_30, id_33, id_35, \
                         id_60, id_63, id_65, id_72, id_75, id_77, id_126, id_129, id_131, \
                         id_138, id_141, id_143, kd_0, kd_3, kd_5, kd_18, kd_21, kd_23, kd_30, \
                         kd_33, kd_35, kd_60, kd_63, kd_65, kd_72, kd_75, kd_77, kd_126, \
                         kd_129, kd_131, kd_138, kd_141, kd_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_176 * ab_x[k] * id_0[k]
                  + f_176 * ab_x[k] * id_3[k]
                  - f_45 * ab_x[k] * id_5[k]
                  - f_177 * ab_x[k] * id_18[k]
                  - f_177 * ab_x[k] * id_21[k]
                  + f_178 * ab_x[k] * id_23[k]
                  - f_179 * ab_x[k] * id_30[k]
                  - f_179 * ab_x[k] * id_33[k]
                  + f_47 * ab_x[k] * id_35[k]
                  - f_177 * ab_x[k] * id_60[k]
                  - f_177 * ab_x[k] * id_63[k]
                  + f_178 * ab_x[k] * id_65[k]
                  + f_180 * ab_x[k] * id_72[k]
                  + f_180 * ab_x[k] * id_75[k]
                  - f_181 * ab_x[k] * id_77[k]
                  + f_176 * ab_x[k] * id_126[k]
                  + f_176 * ab_x[k] * id_129[k]
                  - f_45 * ab_x[k] * id_131[k]
                  - f_179 * ab_x[k] * id_138[k]
                  - f_179 * ab_x[k] * id_141[k]
                  + f_47 * ab_x[k] * id_143[k]
                  + f_176 * kd_0[k]
                  + f_176 * kd_3[k]
                  - f_45 * kd_5[k]
                  - f_177 * kd_18[k]
                  - f_177 * kd_21[k]
                  + f_178 * kd_23[k]
                  - f_179 * kd_30[k]
                  - f_179 * kd_33[k]
                  + f_47 * kd_35[k]
                  - f_177 * kd_60[k]
                  - f_177 * kd_63[k]
                  + f_178 * kd_65[k]
                  + f_180 * kd_72[k]
                  + f_180 * kd_75[k]
                  - f_181 * kd_77[k]
                  + f_176 * kd_126[k]
                  + f_176 * kd_129[k]
                  - f_45 * kd_131[k]
                  - f_179 * kd_138[k]
                  - f_179 * kd_141[k]
                  + f_47 * kd_143[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_2, id_4, id_20, id_22, id_32, id_34, id_62, id_64, \
                         id_74, id_76, id_128, id_130, id_140, id_142, kd_2, kd_10, kd_20, \
                         kd_32, kd_40, kd_52, kd_62, kd_74, kd_94, kd_106, kd_128, kd_140, \
                         kd_172, kd_184 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_189 * ab_x[k] * id_2[k]
                  + f_189 * ab_y[k] * id_4[k]
                  + f_190 * ab_x[k] * id_20[k]
                  - f_190 * ab_y[k] * id_22[k]
                  + f_173 * ab_x[k] * id_32[k]
                  - f_173 * ab_y[k] * id_34[k]
                  + f_190 * ab_x[k] * id_62[k]
                  - f_190 * ab_y[k] * id_64[k]
                  - f_191 * ab_x[k] * id_74[k]
                  + f_191 * ab_y[k] * id_76[k]
                  - f_189 * ab_x[k] * id_128[k]
                  + f_189 * ab_y[k] * id_130[k]
                  + f_173 * ab_x[k] * id_140[k]
                  - f_173 * ab_y[k] * id_142[k]
                  - f_189 * kd_2[k]
                  + f_189 * kd_10[k]
                  + f_190 * kd_20[k]
                  + f_173 * kd_32[k]
                  - f_190 * kd_40[k]
                  - f_173 * kd_52[k]
                  + f_190 * kd_62[k]
                  - f_191 * kd_74[k]
                  - f_190 * kd_94[k]
                  + f_191 * kd_106[k]
                  - f_189 * kd_128[k]
                  + f_173 * kd_140[k]
                  + f_189 * kd_172[k]
                  - f_173 * kd_184[k];
    }

#pragma omp simd aligned(ab_x, id_0, id_3, id_18, id_21, id_30, id_33, id_60, id_63, id_72, \
                         id_75, id_126, id_129, id_138, id_141, kd_0, kd_3, kd_18, kd_21, \
                         kd_30, kd_33, kd_60, kd_63, kd_72, kd_75, kd_126, kd_129, kd_138, \
                         kd_141 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_165 * ab_x[k] * id_0[k]
                  + f_164 * ab_x[k] * id_3[k]
                  + f_167 * ab_x[k] * id_18[k]
                  - f_166 * ab_x[k] * id_21[k]
                  + f_169 * ab_x[k] * id_30[k]
                  - f_168 * ab_x[k] * id_33[k]
                  + f_167 * ab_x[k] * id_60[k]
                  - f_166 * ab_x[k] * id_63[k]
                  - f_171 * ab_x[k] * id_72[k]
                  + f_170 * ab_x[k] * id_75[k]
                  - f_165 * ab_x[k] * id_126[k]
                  + f_164 * ab_x[k] * id_129[k]
                  + f_169 * ab_x[k] * id_138[k]
                  - f_168 * ab_x[k] * id_141[k]
                  - f_165 * kd_0[k]
                  + f_164 * kd_3[k]
                  + f_167 * kd_18[k]
                  - f_166 * kd_21[k]
                  + f_169 * kd_30[k]
                  - f_168 * kd_33[k]
                  + f_167 * kd_60[k]
                  - f_166 * kd_63[k]
                  - f_171 * kd_72[k]
                  + f_170 * kd_75[k]
                  - f_165 * kd_126[k]
                  + f_164 * kd_129[k]
                  + f_169 * kd_138[k]
                  - f_168 * kd_141[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_13, id_15, id_43, id_45, id_97, id_99, kd_13, kd_27, \
                         kd_43, kd_69, kd_97, kd_135 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = f_20 * ab_x[k] * id_13[k]
                  - f_21 * ab_y[k] * id_15[k]
                  - f_18 * ab_x[k] * id_43[k]
                  + f_19 * ab_y[k] * id_45[k]
                  + f_16 * ab_x[k] * id_97[k]
                  - f_17 * ab_y[k] * id_99[k]
                  + f_20 * kd_13[k]
                  - f_21 * kd_27[k]
                  - f_18 * kd_43[k]
                  + f_19 * kd_69[k]
                  + f_16 * kd_97[k]
                  - f_17 * kd_135[k];
    }

#pragma omp simd aligned(ab_x, id_16, id_46, id_100, kd_16, kd_46, \
                         kd_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = f_24 * ab_x[k] * id_16[k]
                  - f_23 * ab_x[k] * id_46[k]
                  + f_22 * ab_x[k] * id_100[k]
                  + f_24 * kd_16[k]
                  - f_23 * kd_46[k]
                  + f_22 * kd_100[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_13, id_15, id_17, id_43, id_45, id_47, id_97, id_99, \
                         id_101, kd_13, kd_27, kd_29, kd_43, kd_69, kd_71, kd_97, kd_135, \
                         kd_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = -f_29 * ab_x[k] * id_13[k]
                  - f_29 * ab_y[k] * id_15[k]
                  + f_30 * ab_y[k] * id_17[k]
                  + f_27 * ab_x[k] * id_43[k]
                  + f_27 * ab_y[k] * id_45[k]
                  - f_28 * ab_y[k] * id_47[k]
                  - f_25 * ab_x[k] * id_97[k]
                  - f_25 * ab_y[k] * id_99[k]
                  + f_26 * ab_y[k] * id_101[k]
                  - f_29 * kd_13[k]
                  - f_29 * kd_27[k]
                  + f_30 * kd_29[k]
                  + f_27 * kd_43[k]
                  + f_27 * kd_69[k]
                  - f_28 * kd_71[k]
                  - f_25 * kd_97[k]
                  - f_25 * kd_135[k]
                  + f_26 * kd_137[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, id_14, id_16, id_17, id_44, id_46, id_47, id_98, \
                         id_100, id_101, kd_14, kd_28, kd_35, kd_44, kd_70, kd_77, kd_98, \
                         kd_136, kd_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = -f_35 * ab_x[k] * id_14[k]
                  - f_35 * ab_y[k] * id_16[k]
                  + f_36 * ab_z[k] * id_17[k]
                  + f_33 * ab_x[k] * id_44[k]
                  + f_33 * ab_y[k] * id_46[k]
                  - f_34 * ab_z[k] * id_47[k]
                  - f_31 * ab_x[k] * id_98[k]
                  - f_31 * ab_y[k] * id_100[k]
                  + f_32 * ab_z[k] * id_101[k]
                  - f_35 * kd_14[k]
                  - f_35 * kd_28[k]
                  + f_36 * kd_35[k]
                  + f_33 * kd_44[k]
                  + f_33 * kd_70[k]
                  - f_34 * kd_77[k]
                  - f_31 * kd_98[k]
                  - f_31 * kd_136[k]
                  + f_32 * kd_143[k];
    }

#pragma omp simd aligned(ab_x, id_12, id_15, id_17, id_42, id_45, id_47, id_96, id_99, id_101, \
                         kd_12, kd_15, kd_17, kd_42, kd_45, kd_47, kd_96, kd_99, \
                         kd_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = -f_29 * ab_x[k] * id_12[k]
                  - f_29 * ab_x[k] * id_15[k]
                  + f_30 * ab_x[k] * id_17[k]
                  + f_27 * ab_x[k] * id_42[k]
                  + f_27 * ab_x[k] * id_45[k]
                  - f_28 * ab_x[k] * id_47[k]
                  - f_25 * ab_x[k] * id_96[k]
                  - f_25 * ab_x[k] * id_99[k]
                  + f_26 * ab_x[k] * id_101[k]
                  - f_29 * kd_12[k]
                  - f_29 * kd_15[k]
                  + f_30 * kd_17[k]
                  + f_27 * kd_42[k]
                  + f_27 * kd_45[k]
                  - f_28 * kd_47[k]
                  - f_25 * kd_96[k]
                  - f_25 * kd_99[k]
                  + f_26 * kd_101[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_14, id_16, id_44, id_46, id_98, id_100, kd_14, kd_28, \
                         kd_44, kd_70, kd_98, kd_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = f_38 * ab_x[k] * id_14[k]
                  - f_38 * ab_y[k] * id_16[k]
                  - f_22 * ab_x[k] * id_44[k]
                  + f_22 * ab_y[k] * id_46[k]
                  + f_37 * ab_x[k] * id_98[k]
                  - f_37 * ab_y[k] * id_100[k]
                  + f_38 * kd_14[k]
                  - f_38 * kd_28[k]
                  - f_22 * kd_44[k]
                  + f_22 * kd_70[k]
                  + f_37 * kd_98[k]
                  - f_37 * kd_136[k];
    }

#pragma omp simd aligned(ab_x, id_12, id_15, id_42, id_45, id_96, id_99, kd_12, kd_15, kd_42, \
                         kd_45, kd_96, kd_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = f_21 * ab_x[k] * id_12[k]
                  - f_20 * ab_x[k] * id_15[k]
                  - f_19 * ab_x[k] * id_42[k]
                  + f_18 * ab_x[k] * id_45[k]
                  + f_17 * ab_x[k] * id_96[k]
                  - f_16 * ab_x[k] * id_99[k]
                  + f_21 * kd_12[k]
                  - f_20 * kd_15[k]
                  - f_19 * kd_42[k]
                  + f_18 * kd_45[k]
                  + f_17 * kd_96[k]
                  - f_16 * kd_99[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_1, id_3, id_19, id_21, id_61, id_63, id_127, id_129, \
                         kd_1, kd_9, kd_19, kd_39, kd_61, kd_93, kd_127, \
                         kd_171 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = f_192 * ab_x[k] * id_1[k]
                  - f_193 * ab_y[k] * id_3[k]
                  - f_194 * ab_x[k] * id_19[k]
                  + f_195 * ab_y[k] * id_21[k]
                  + f_194 * ab_x[k] * id_61[k]
                  - f_195 * ab_y[k] * id_63[k]
                  - f_192 * ab_x[k] * id_127[k]
                  + f_193 * ab_y[k] * id_129[k]
                  + f_192 * kd_1[k]
                  - f_193 * kd_9[k]
                  - f_194 * kd_19[k]
                  + f_195 * kd_39[k]
                  + f_194 * kd_61[k]
                  - f_195 * kd_93[k]
                  - f_192 * kd_127[k]
                  + f_193 * kd_171[k];
    }

#pragma omp simd aligned(ab_x, id_4, id_22, id_64, id_130, kd_4, kd_22, kd_64, \
                         kd_130 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = f_196 * ab_x[k] * id_4[k]
                  - f_197 * ab_x[k] * id_22[k]
                  + f_197 * ab_x[k] * id_64[k]
                  - f_196 * ab_x[k] * id_130[k]
                  + f_196 * kd_4[k]
                  - f_197 * kd_22[k]
                  + f_197 * kd_64[k]
                  - f_196 * kd_130[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_1, id_3, id_5, id_19, id_21, id_23, id_61, id_63, \
                         id_65, id_127, id_129, id_131, kd_1, kd_9, kd_11, kd_19, kd_39, \
                         kd_41, kd_61, kd_93, kd_95, kd_127, kd_171, \
                         kd_173 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_198 * ab_x[k] * id_1[k]
                  - f_198 * ab_y[k] * id_3[k]
                  + f_199 * ab_y[k] * id_5[k]
                  + f_200 * ab_x[k] * id_19[k]
                  + f_200 * ab_y[k] * id_21[k]
                  - f_201 * ab_y[k] * id_23[k]
                  - f_200 * ab_x[k] * id_61[k]
                  - f_200 * ab_y[k] * id_63[k]
                  + f_201 * ab_y[k] * id_65[k]
                  + f_198 * ab_x[k] * id_127[k]
                  + f_198 * ab_y[k] * id_129[k]
                  - f_199 * ab_y[k] * id_131[k]
                  - f_198 * kd_1[k]
                  - f_198 * kd_9[k]
                  + f_199 * kd_11[k]
                  + f_200 * kd_19[k]
                  + f_200 * kd_39[k]
                  - f_201 * kd_41[k]
                  - f_200 * kd_61[k]
                  - f_200 * kd_93[k]
                  + f_201 * kd_95[k]
                  + f_198 * kd_127[k]
                  + f_198 * kd_171[k]
                  - f_199 * kd_173[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, id_2, id_4, id_5, id_20, id_22, id_23, id_62, \
                         id_64, id_65, id_128, id_130, id_131, kd_2, kd_10, kd_17, kd_20, \
                         kd_40, kd_47, kd_62, kd_94, kd_101, kd_128, kd_172, \
                         kd_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = -f_202 * ab_x[k] * id_2[k]
                  - f_202 * ab_y[k] * id_4[k]
                  + f_203 * ab_z[k] * id_5[k]
                  + f_204 * ab_x[k] * id_20[k]
                  + f_204 * ab_y[k] * id_22[k]
                  - f_205 * ab_z[k] * id_23[k]
                  - f_204 * ab_x[k] * id_62[k]
                  - f_204 * ab_y[k] * id_64[k]
                  + f_205 * ab_z[k] * id_65[k]
                  + f_202 * ab_x[k] * id_128[k]
                  + f_202 * ab_y[k] * id_130[k]
                  - f_203 * ab_z[k] * id_131[k]
                  - f_202 * kd_2[k]
                  - f_202 * kd_10[k]
                  + f_203 * kd_17[k]
                  + f_204 * kd_20[k]
                  + f_204 * kd_40[k]
                  - f_205 * kd_47[k]
                  - f_204 * kd_62[k]
                  - f_204 * kd_94[k]
                  + f_205 * kd_101[k]
                  + f_202 * kd_128[k]
                  + f_202 * kd_172[k]
                  - f_203 * kd_179[k];
    }

#pragma omp simd aligned(ab_x, id_0, id_3, id_5, id_18, id_21, id_23, id_60, id_63, id_65, \
                         id_126, id_129, id_131, kd_0, kd_3, kd_5, kd_18, kd_21, kd_23, kd_60, \
                         kd_63, kd_65, kd_126, kd_129, kd_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = -f_198 * ab_x[k] * id_0[k]
                  - f_198 * ab_x[k] * id_3[k]
                  + f_199 * ab_x[k] * id_5[k]
                  + f_200 * ab_x[k] * id_18[k]
                  + f_200 * ab_x[k] * id_21[k]
                  - f_201 * ab_x[k] * id_23[k]
                  - f_200 * ab_x[k] * id_60[k]
                  - f_200 * ab_x[k] * id_63[k]
                  + f_201 * ab_x[k] * id_65[k]
                  + f_198 * ab_x[k] * id_126[k]
                  + f_198 * ab_x[k] * id_129[k]
                  - f_199 * ab_x[k] * id_131[k]
                  - f_198 * kd_0[k]
                  - f_198 * kd_3[k]
                  + f_199 * kd_5[k]
                  + f_200 * kd_18[k]
                  + f_200 * kd_21[k]
                  - f_201 * kd_23[k]
                  - f_200 * kd_60[k]
                  - f_200 * kd_63[k]
                  + f_201 * kd_65[k]
                  + f_198 * kd_126[k]
                  + f_198 * kd_129[k]
                  - f_199 * kd_131[k];
    }

#pragma omp simd aligned(ab_x, ab_y, id_2, id_4, id_20, id_22, id_62, id_64, id_128, id_130, \
                         kd_2, kd_10, kd_20, kd_40, kd_62, kd_94, kd_128, \
                         kd_172 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = f_206 * ab_x[k] * id_2[k]
                  - f_206 * ab_y[k] * id_4[k]
                  - f_207 * ab_x[k] * id_20[k]
                  + f_207 * ab_y[k] * id_22[k]
                  + f_207 * ab_x[k] * id_62[k]
                  - f_207 * ab_y[k] * id_64[k]
                  - f_206 * ab_x[k] * id_128[k]
                  + f_206 * ab_y[k] * id_130[k]
                  + f_206 * kd_2[k]
                  - f_206 * kd_10[k]
                  - f_207 * kd_20[k]
                  + f_207 * kd_40[k]
                  + f_207 * kd_62[k]
                  - f_207 * kd_94[k]
                  - f_206 * kd_128[k]
                  + f_206 * kd_172[k];
    }

#pragma omp simd aligned(ab_x, id_0, id_3, id_18, id_21, id_60, id_63, id_126, id_129, kd_0, \
                         kd_3, kd_18, kd_21, kd_60, kd_63, kd_126, \
                         kd_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = f_193 * ab_x[k] * id_0[k]
                  - f_192 * ab_x[k] * id_3[k]
                  - f_195 * ab_x[k] * id_18[k]
                  + f_194 * ab_x[k] * id_21[k]
                  + f_195 * ab_x[k] * id_60[k]
                  - f_194 * ab_x[k] * id_63[k]
                  - f_193 * ab_x[k] * id_126[k]
                  + f_192 * ab_x[k] * id_129[k]
                  + f_193 * kd_0[k]
                  - f_192 * kd_3[k]
                  - f_195 * kd_18[k]
                  + f_194 * kd_21[k]
                  + f_195 * kd_60[k]
                  - f_194 * kd_63[k]
                  - f_193 * kd_126[k]
                  + f_192 * kd_129[k];
    }
}

}  // namespace simdovl
