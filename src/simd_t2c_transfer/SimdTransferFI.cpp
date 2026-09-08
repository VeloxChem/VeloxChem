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


#include "SimdTransferFI.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_fi(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t di, const size_t dk, const size_t nmax) -> void
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
    const auto *di_30 = buffer.data(di + 30);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_35 = buffer.data(di + 35);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_39 = buffer.data(di + 39);
    const auto *di_40 = buffer.data(di + 40);
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
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_57 = buffer.data(di + 57);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_59 = buffer.data(di + 59);
    const auto *di_60 = buffer.data(di + 60);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_62 = buffer.data(di + 62);
    const auto *di_63 = buffer.data(di + 63);
    const auto *di_64 = buffer.data(di + 64);
    const auto *di_65 = buffer.data(di + 65);
    const auto *di_66 = buffer.data(di + 66);
    const auto *di_67 = buffer.data(di + 67);
    const auto *di_68 = buffer.data(di + 68);
    const auto *di_69 = buffer.data(di + 69);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_71 = buffer.data(di + 71);
    const auto *di_72 = buffer.data(di + 72);
    const auto *di_73 = buffer.data(di + 73);
    const auto *di_74 = buffer.data(di + 74);
    const auto *di_75 = buffer.data(di + 75);
    const auto *di_76 = buffer.data(di + 76);
    const auto *di_77 = buffer.data(di + 77);
    const auto *di_78 = buffer.data(di + 78);
    const auto *di_79 = buffer.data(di + 79);
    const auto *di_80 = buffer.data(di + 80);
    const auto *di_81 = buffer.data(di + 81);
    const auto *di_82 = buffer.data(di + 82);
    const auto *di_83 = buffer.data(di + 83);
    const auto *di_84 = buffer.data(di + 84);
    const auto *di_85 = buffer.data(di + 85);
    const auto *di_86 = buffer.data(di + 86);
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
    const auto *di_99 = buffer.data(di + 99);
    const auto *di_100 = buffer.data(di + 100);
    const auto *di_101 = buffer.data(di + 101);
    const auto *di_102 = buffer.data(di + 102);
    const auto *di_103 = buffer.data(di + 103);
    const auto *di_104 = buffer.data(di + 104);
    const auto *di_105 = buffer.data(di + 105);
    const auto *di_106 = buffer.data(di + 106);
    const auto *di_107 = buffer.data(di + 107);
    const auto *di_108 = buffer.data(di + 108);
    const auto *di_109 = buffer.data(di + 109);
    const auto *di_110 = buffer.data(di + 110);
    const auto *di_111 = buffer.data(di + 111);
    const auto *di_112 = buffer.data(di + 112);
    const auto *di_113 = buffer.data(di + 113);
    const auto *di_114 = buffer.data(di + 114);
    const auto *di_115 = buffer.data(di + 115);
    const auto *di_116 = buffer.data(di + 116);
    const auto *di_117 = buffer.data(di + 117);
    const auto *di_118 = buffer.data(di + 118);
    const auto *di_119 = buffer.data(di + 119);
    const auto *di_120 = buffer.data(di + 120);
    const auto *di_121 = buffer.data(di + 121);
    const auto *di_122 = buffer.data(di + 122);
    const auto *di_123 = buffer.data(di + 123);
    const auto *di_124 = buffer.data(di + 124);
    const auto *di_125 = buffer.data(di + 125);
    const auto *di_126 = buffer.data(di + 126);
    const auto *di_127 = buffer.data(di + 127);
    const auto *di_128 = buffer.data(di + 128);
    const auto *di_129 = buffer.data(di + 129);
    const auto *di_130 = buffer.data(di + 130);
    const auto *di_131 = buffer.data(di + 131);
    const auto *di_132 = buffer.data(di + 132);
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
    const auto *di_155 = buffer.data(di + 155);
    const auto *di_156 = buffer.data(di + 156);
    const auto *di_157 = buffer.data(di + 157);
    const auto *di_158 = buffer.data(di + 158);
    const auto *di_159 = buffer.data(di + 159);
    const auto *di_160 = buffer.data(di + 160);
    const auto *di_161 = buffer.data(di + 161);
    const auto *di_162 = buffer.data(di + 162);
    const auto *di_163 = buffer.data(di + 163);
    const auto *di_164 = buffer.data(di + 164);
    const auto *di_165 = buffer.data(di + 165);
    const auto *di_166 = buffer.data(di + 166);
    const auto *di_167 = buffer.data(di + 167);

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
    const auto *dk_129 = buffer.data(dk + 129);
    const auto *dk_130 = buffer.data(dk + 130);
    const auto *dk_131 = buffer.data(dk + 131);
    const auto *dk_132 = buffer.data(dk + 132);
    const auto *dk_133 = buffer.data(dk + 133);
    const auto *dk_134 = buffer.data(dk + 134);
    const auto *dk_135 = buffer.data(dk + 135);
    const auto *dk_136 = buffer.data(dk + 136);
    const auto *dk_137 = buffer.data(dk + 137);
    const auto *dk_138 = buffer.data(dk + 138);
    const auto *dk_139 = buffer.data(dk + 139);
    const auto *dk_140 = buffer.data(dk + 140);
    const auto *dk_141 = buffer.data(dk + 141);
    const auto *dk_142 = buffer.data(dk + 142);
    const auto *dk_144 = buffer.data(dk + 144);
    const auto *dk_145 = buffer.data(dk + 145);
    const auto *dk_146 = buffer.data(dk + 146);
    const auto *dk_147 = buffer.data(dk + 147);
    const auto *dk_148 = buffer.data(dk + 148);
    const auto *dk_149 = buffer.data(dk + 149);
    const auto *dk_150 = buffer.data(dk + 150);
    const auto *dk_151 = buffer.data(dk + 151);
    const auto *dk_152 = buffer.data(dk + 152);
    const auto *dk_153 = buffer.data(dk + 153);
    const auto *dk_154 = buffer.data(dk + 154);
    const auto *dk_155 = buffer.data(dk + 155);
    const auto *dk_156 = buffer.data(dk + 156);
    const auto *dk_157 = buffer.data(dk + 157);
    const auto *dk_158 = buffer.data(dk + 158);
    const auto *dk_159 = buffer.data(dk + 159);
    const auto *dk_160 = buffer.data(dk + 160);
    const auto *dk_161 = buffer.data(dk + 161);
    const auto *dk_162 = buffer.data(dk + 162);
    const auto *dk_163 = buffer.data(dk + 163);
    const auto *dk_164 = buffer.data(dk + 164);
    const auto *dk_165 = buffer.data(dk + 165);
    const auto *dk_166 = buffer.data(dk + 166);
    const auto *dk_167 = buffer.data(dk + 167);
    const auto *dk_168 = buffer.data(dk + 168);
    const auto *dk_169 = buffer.data(dk + 169);
    const auto *dk_170 = buffer.data(dk + 170);
    const auto *dk_171 = buffer.data(dk + 171);
    const auto *dk_172 = buffer.data(dk + 172);
    const auto *dk_173 = buffer.data(dk + 173);
    const auto *dk_174 = buffer.data(dk + 174);
    const auto *dk_175 = buffer.data(dk + 175);
    const auto *dk_176 = buffer.data(dk + 176);
    const auto *dk_177 = buffer.data(dk + 177);
    const auto *dk_178 = buffer.data(dk + 178);
    const auto *dk_180 = buffer.data(dk + 180);
    const auto *dk_181 = buffer.data(dk + 181);
    const auto *dk_182 = buffer.data(dk + 182);
    const auto *dk_183 = buffer.data(dk + 183);
    const auto *dk_184 = buffer.data(dk + 184);
    const auto *dk_185 = buffer.data(dk + 185);
    const auto *dk_186 = buffer.data(dk + 186);
    const auto *dk_187 = buffer.data(dk + 187);
    const auto *dk_188 = buffer.data(dk + 188);
    const auto *dk_189 = buffer.data(dk + 189);
    const auto *dk_190 = buffer.data(dk + 190);
    const auto *dk_191 = buffer.data(dk + 191);
    const auto *dk_192 = buffer.data(dk + 192);
    const auto *dk_193 = buffer.data(dk + 193);
    const auto *dk_194 = buffer.data(dk + 194);
    const auto *dk_195 = buffer.data(dk + 195);
    const auto *dk_196 = buffer.data(dk + 196);
    const auto *dk_197 = buffer.data(dk + 197);
    const auto *dk_198 = buffer.data(dk + 198);
    const auto *dk_199 = buffer.data(dk + 199);
    const auto *dk_200 = buffer.data(dk + 200);
    const auto *dk_201 = buffer.data(dk + 201);
    const auto *dk_202 = buffer.data(dk + 202);
    const auto *dk_203 = buffer.data(dk + 203);
    const auto *dk_204 = buffer.data(dk + 204);
    const auto *dk_205 = buffer.data(dk + 205);
    const auto *dk_206 = buffer.data(dk + 206);
    const auto *dk_207 = buffer.data(dk + 207);
    const auto *dk_208 = buffer.data(dk + 208);
    const auto *dk_209 = buffer.data(dk + 209);
    const auto *dk_210 = buffer.data(dk + 210);
    const auto *dk_211 = buffer.data(dk + 211);
    const auto *dk_212 = buffer.data(dk + 212);
    const auto *dk_213 = buffer.data(dk + 213);
    const auto *dk_214 = buffer.data(dk + 214);
    const auto *dk_215 = buffer.data(dk + 215);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, di_0, di_1, di_2, di_3, di_4, dk_0, \
                         dk_1, dk_2, dk_3, dk_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * di_0[k]
                 + dk_0[k];

        t_1[k] = -ab_x[k] * di_1[k]
                 + dk_1[k];

        t_2[k] = -ab_x[k] * di_2[k]
                 + dk_2[k];

        t_3[k] = -ab_x[k] * di_3[k]
                 + dk_3[k];

        t_4[k] = -ab_x[k] * di_4[k]
                 + dk_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, di_5, di_6, di_7, di_8, di_9, dk_5, \
                         dk_6, dk_7, dk_8, dk_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * di_5[k]
                 + dk_5[k];

        t_6[k] = -ab_x[k] * di_6[k]
                 + dk_6[k];

        t_7[k] = -ab_x[k] * di_7[k]
                 + dk_7[k];

        t_8[k] = -ab_x[k] * di_8[k]
                 + dk_8[k];

        t_9[k] = -ab_x[k] * di_9[k]
                 + dk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, di_10, di_11, di_12, di_13, \
                         di_14, dk_10, dk_11, dk_12, dk_13, dk_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * di_10[k]
                  + dk_10[k];

        t_11[k] = -ab_x[k] * di_11[k]
                  + dk_11[k];

        t_12[k] = -ab_x[k] * di_12[k]
                  + dk_12[k];

        t_13[k] = -ab_x[k] * di_13[k]
                  + dk_13[k];

        t_14[k] = -ab_x[k] * di_14[k]
                  + dk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, di_15, di_16, di_17, di_18, \
                         di_19, dk_15, dk_16, dk_17, dk_18, dk_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * di_15[k]
                  + dk_15[k];

        t_16[k] = -ab_x[k] * di_16[k]
                  + dk_16[k];

        t_17[k] = -ab_x[k] * di_17[k]
                  + dk_17[k];

        t_18[k] = -ab_x[k] * di_18[k]
                  + dk_18[k];

        t_19[k] = -ab_x[k] * di_19[k]
                  + dk_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, di_20, di_21, di_22, di_23, \
                         di_24, dk_20, dk_21, dk_22, dk_23, dk_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * di_20[k]
                  + dk_20[k];

        t_21[k] = -ab_x[k] * di_21[k]
                  + dk_21[k];

        t_22[k] = -ab_x[k] * di_22[k]
                  + dk_22[k];

        t_23[k] = -ab_x[k] * di_23[k]
                  + dk_23[k];

        t_24[k] = -ab_x[k] * di_24[k]
                  + dk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, di_25, di_26, di_27, di_28, \
                         di_29, dk_25, dk_26, dk_27, dk_36, dk_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * di_25[k]
                  + dk_25[k];

        t_26[k] = -ab_x[k] * di_26[k]
                  + dk_26[k];

        t_27[k] = -ab_x[k] * di_27[k]
                  + dk_27[k];

        t_28[k] = -ab_x[k] * di_28[k]
                  + dk_36[k];

        t_29[k] = -ab_x[k] * di_29[k]
                  + dk_37[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, di_30, di_31, di_32, di_33, \
                         di_34, dk_38, dk_39, dk_40, dk_41, dk_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_x[k] * di_30[k]
                  + dk_38[k];

        t_31[k] = -ab_x[k] * di_31[k]
                  + dk_39[k];

        t_32[k] = -ab_x[k] * di_32[k]
                  + dk_40[k];

        t_33[k] = -ab_x[k] * di_33[k]
                  + dk_41[k];

        t_34[k] = -ab_x[k] * di_34[k]
                  + dk_42[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, di_35, di_36, di_37, di_38, \
                         di_39, dk_43, dk_44, dk_45, dk_46, dk_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_x[k] * di_35[k]
                  + dk_43[k];

        t_36[k] = -ab_x[k] * di_36[k]
                  + dk_44[k];

        t_37[k] = -ab_x[k] * di_37[k]
                  + dk_45[k];

        t_38[k] = -ab_x[k] * di_38[k]
                  + dk_46[k];

        t_39[k] = -ab_x[k] * di_39[k]
                  + dk_47[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, di_40, di_41, di_42, di_43, \
                         di_44, dk_48, dk_49, dk_50, dk_51, dk_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_x[k] * di_40[k]
                  + dk_48[k];

        t_41[k] = -ab_x[k] * di_41[k]
                  + dk_49[k];

        t_42[k] = -ab_x[k] * di_42[k]
                  + dk_50[k];

        t_43[k] = -ab_x[k] * di_43[k]
                  + dk_51[k];

        t_44[k] = -ab_x[k] * di_44[k]
                  + dk_52[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, di_45, di_46, di_47, di_48, \
                         di_49, dk_53, dk_54, dk_55, dk_56, dk_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = -ab_x[k] * di_45[k]
                  + dk_53[k];

        t_46[k] = -ab_x[k] * di_46[k]
                  + dk_54[k];

        t_47[k] = -ab_x[k] * di_47[k]
                  + dk_55[k];

        t_48[k] = -ab_x[k] * di_48[k]
                  + dk_56[k];

        t_49[k] = -ab_x[k] * di_49[k]
                  + dk_57[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, di_50, di_51, di_52, di_53, \
                         di_54, dk_58, dk_59, dk_60, dk_61, dk_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = -ab_x[k] * di_50[k]
                  + dk_58[k];

        t_51[k] = -ab_x[k] * di_51[k]
                  + dk_59[k];

        t_52[k] = -ab_x[k] * di_52[k]
                  + dk_60[k];

        t_53[k] = -ab_x[k] * di_53[k]
                  + dk_61[k];

        t_54[k] = -ab_x[k] * di_54[k]
                  + dk_62[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, di_55, di_56, di_57, di_58, \
                         di_59, dk_63, dk_72, dk_73, dk_74, dk_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = -ab_x[k] * di_55[k]
                  + dk_63[k];

        t_56[k] = -ab_x[k] * di_56[k]
                  + dk_72[k];

        t_57[k] = -ab_x[k] * di_57[k]
                  + dk_73[k];

        t_58[k] = -ab_x[k] * di_58[k]
                  + dk_74[k];

        t_59[k] = -ab_x[k] * di_59[k]
                  + dk_75[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, di_60, di_61, di_62, di_63, \
                         di_64, dk_76, dk_77, dk_78, dk_79, dk_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = -ab_x[k] * di_60[k]
                  + dk_76[k];

        t_61[k] = -ab_x[k] * di_61[k]
                  + dk_77[k];

        t_62[k] = -ab_x[k] * di_62[k]
                  + dk_78[k];

        t_63[k] = -ab_x[k] * di_63[k]
                  + dk_79[k];

        t_64[k] = -ab_x[k] * di_64[k]
                  + dk_80[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, di_65, di_66, di_67, di_68, \
                         di_69, dk_81, dk_82, dk_83, dk_84, dk_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = -ab_x[k] * di_65[k]
                  + dk_81[k];

        t_66[k] = -ab_x[k] * di_66[k]
                  + dk_82[k];

        t_67[k] = -ab_x[k] * di_67[k]
                  + dk_83[k];

        t_68[k] = -ab_x[k] * di_68[k]
                  + dk_84[k];

        t_69[k] = -ab_x[k] * di_69[k]
                  + dk_85[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, di_70, di_71, di_72, di_73, \
                         di_74, dk_86, dk_87, dk_88, dk_89, dk_90 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = -ab_x[k] * di_70[k]
                  + dk_86[k];

        t_71[k] = -ab_x[k] * di_71[k]
                  + dk_87[k];

        t_72[k] = -ab_x[k] * di_72[k]
                  + dk_88[k];

        t_73[k] = -ab_x[k] * di_73[k]
                  + dk_89[k];

        t_74[k] = -ab_x[k] * di_74[k]
                  + dk_90[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, di_75, di_76, di_77, di_78, \
                         di_79, dk_91, dk_92, dk_93, dk_94, dk_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = -ab_x[k] * di_75[k]
                  + dk_91[k];

        t_76[k] = -ab_x[k] * di_76[k]
                  + dk_92[k];

        t_77[k] = -ab_x[k] * di_77[k]
                  + dk_93[k];

        t_78[k] = -ab_x[k] * di_78[k]
                  + dk_94[k];

        t_79[k] = -ab_x[k] * di_79[k]
                  + dk_95[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, di_80, di_81, di_82, di_83, \
                         di_84, dk_96, dk_97, dk_98, dk_99, dk_108 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = -ab_x[k] * di_80[k]
                  + dk_96[k];

        t_81[k] = -ab_x[k] * di_81[k]
                  + dk_97[k];

        t_82[k] = -ab_x[k] * di_82[k]
                  + dk_98[k];

        t_83[k] = -ab_x[k] * di_83[k]
                  + dk_99[k];

        t_84[k] = -ab_x[k] * di_84[k]
                  + dk_108[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, di_85, di_86, di_87, di_88, \
                         di_89, dk_109, dk_110, dk_111, dk_112, \
                         dk_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = -ab_x[k] * di_85[k]
                  + dk_109[k];

        t_86[k] = -ab_x[k] * di_86[k]
                  + dk_110[k];

        t_87[k] = -ab_x[k] * di_87[k]
                  + dk_111[k];

        t_88[k] = -ab_x[k] * di_88[k]
                  + dk_112[k];

        t_89[k] = -ab_x[k] * di_89[k]
                  + dk_113[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, di_90, di_91, di_92, di_93, \
                         di_94, dk_114, dk_115, dk_116, dk_117, \
                         dk_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = -ab_x[k] * di_90[k]
                  + dk_114[k];

        t_91[k] = -ab_x[k] * di_91[k]
                  + dk_115[k];

        t_92[k] = -ab_x[k] * di_92[k]
                  + dk_116[k];

        t_93[k] = -ab_x[k] * di_93[k]
                  + dk_117[k];

        t_94[k] = -ab_x[k] * di_94[k]
                  + dk_118[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, di_95, di_96, di_97, di_98, \
                         di_99, dk_119, dk_120, dk_121, dk_122, \
                         dk_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = -ab_x[k] * di_95[k]
                  + dk_119[k];

        t_96[k] = -ab_x[k] * di_96[k]
                  + dk_120[k];

        t_97[k] = -ab_x[k] * di_97[k]
                  + dk_121[k];

        t_98[k] = -ab_x[k] * di_98[k]
                  + dk_122[k];

        t_99[k] = -ab_x[k] * di_99[k]
                  + dk_123[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, di_100, di_101, di_102, \
                         di_103, di_104, dk_124, dk_125, dk_126, dk_127, \
                         dk_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_100[k] = -ab_x[k] * di_100[k]
                   + dk_124[k];

        t_101[k] = -ab_x[k] * di_101[k]
                   + dk_125[k];

        t_102[k] = -ab_x[k] * di_102[k]
                   + dk_126[k];

        t_103[k] = -ab_x[k] * di_103[k]
                   + dk_127[k];

        t_104[k] = -ab_x[k] * di_104[k]
                   + dk_128[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, di_105, di_106, di_107, \
                         di_108, di_109, dk_129, dk_130, dk_131, dk_132, \
                         dk_133 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_105[k] = -ab_x[k] * di_105[k]
                   + dk_129[k];

        t_106[k] = -ab_x[k] * di_106[k]
                   + dk_130[k];

        t_107[k] = -ab_x[k] * di_107[k]
                   + dk_131[k];

        t_108[k] = -ab_x[k] * di_108[k]
                   + dk_132[k];

        t_109[k] = -ab_x[k] * di_109[k]
                   + dk_133[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, di_110, di_111, di_112, \
                         di_113, di_114, dk_134, dk_135, dk_144, dk_145, \
                         dk_146 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_110[k] = -ab_x[k] * di_110[k]
                   + dk_134[k];

        t_111[k] = -ab_x[k] * di_111[k]
                   + dk_135[k];

        t_112[k] = -ab_x[k] * di_112[k]
                   + dk_144[k];

        t_113[k] = -ab_x[k] * di_113[k]
                   + dk_145[k];

        t_114[k] = -ab_x[k] * di_114[k]
                   + dk_146[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, di_115, di_116, di_117, \
                         di_118, di_119, dk_147, dk_148, dk_149, dk_150, \
                         dk_151 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_115[k] = -ab_x[k] * di_115[k]
                   + dk_147[k];

        t_116[k] = -ab_x[k] * di_116[k]
                   + dk_148[k];

        t_117[k] = -ab_x[k] * di_117[k]
                   + dk_149[k];

        t_118[k] = -ab_x[k] * di_118[k]
                   + dk_150[k];

        t_119[k] = -ab_x[k] * di_119[k]
                   + dk_151[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, di_120, di_121, di_122, \
                         di_123, di_124, dk_152, dk_153, dk_154, dk_155, \
                         dk_156 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_120[k] = -ab_x[k] * di_120[k]
                   + dk_152[k];

        t_121[k] = -ab_x[k] * di_121[k]
                   + dk_153[k];

        t_122[k] = -ab_x[k] * di_122[k]
                   + dk_154[k];

        t_123[k] = -ab_x[k] * di_123[k]
                   + dk_155[k];

        t_124[k] = -ab_x[k] * di_124[k]
                   + dk_156[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, di_125, di_126, di_127, \
                         di_128, di_129, dk_157, dk_158, dk_159, dk_160, \
                         dk_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_125[k] = -ab_x[k] * di_125[k]
                   + dk_157[k];

        t_126[k] = -ab_x[k] * di_126[k]
                   + dk_158[k];

        t_127[k] = -ab_x[k] * di_127[k]
                   + dk_159[k];

        t_128[k] = -ab_x[k] * di_128[k]
                   + dk_160[k];

        t_129[k] = -ab_x[k] * di_129[k]
                   + dk_161[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, di_130, di_131, di_132, \
                         di_133, di_134, dk_162, dk_163, dk_164, dk_165, \
                         dk_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_130[k] = -ab_x[k] * di_130[k]
                   + dk_162[k];

        t_131[k] = -ab_x[k] * di_131[k]
                   + dk_163[k];

        t_132[k] = -ab_x[k] * di_132[k]
                   + dk_164[k];

        t_133[k] = -ab_x[k] * di_133[k]
                   + dk_165[k];

        t_134[k] = -ab_x[k] * di_134[k]
                   + dk_166[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, di_135, di_136, di_137, \
                         di_138, di_139, dk_167, dk_168, dk_169, dk_170, \
                         dk_171 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = -ab_x[k] * di_135[k]
                   + dk_167[k];

        t_136[k] = -ab_x[k] * di_136[k]
                   + dk_168[k];

        t_137[k] = -ab_x[k] * di_137[k]
                   + dk_169[k];

        t_138[k] = -ab_x[k] * di_138[k]
                   + dk_170[k];

        t_139[k] = -ab_x[k] * di_139[k]
                   + dk_171[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, di_140, di_141, di_142, \
                         di_143, di_144, dk_180, dk_181, dk_182, dk_183, \
                         dk_184 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_140[k] = -ab_x[k] * di_140[k]
                   + dk_180[k];

        t_141[k] = -ab_x[k] * di_141[k]
                   + dk_181[k];

        t_142[k] = -ab_x[k] * di_142[k]
                   + dk_182[k];

        t_143[k] = -ab_x[k] * di_143[k]
                   + dk_183[k];

        t_144[k] = -ab_x[k] * di_144[k]
                   + dk_184[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, di_145, di_146, di_147, \
                         di_148, di_149, dk_185, dk_186, dk_187, dk_188, \
                         dk_189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_145[k] = -ab_x[k] * di_145[k]
                   + dk_185[k];

        t_146[k] = -ab_x[k] * di_146[k]
                   + dk_186[k];

        t_147[k] = -ab_x[k] * di_147[k]
                   + dk_187[k];

        t_148[k] = -ab_x[k] * di_148[k]
                   + dk_188[k];

        t_149[k] = -ab_x[k] * di_149[k]
                   + dk_189[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, di_150, di_151, di_152, \
                         di_153, di_154, dk_190, dk_191, dk_192, dk_193, \
                         dk_194 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_150[k] = -ab_x[k] * di_150[k]
                   + dk_190[k];

        t_151[k] = -ab_x[k] * di_151[k]
                   + dk_191[k];

        t_152[k] = -ab_x[k] * di_152[k]
                   + dk_192[k];

        t_153[k] = -ab_x[k] * di_153[k]
                   + dk_193[k];

        t_154[k] = -ab_x[k] * di_154[k]
                   + dk_194[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, di_155, di_156, di_157, \
                         di_158, di_159, dk_195, dk_196, dk_197, dk_198, \
                         dk_199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_155[k] = -ab_x[k] * di_155[k]
                   + dk_195[k];

        t_156[k] = -ab_x[k] * di_156[k]
                   + dk_196[k];

        t_157[k] = -ab_x[k] * di_157[k]
                   + dk_197[k];

        t_158[k] = -ab_x[k] * di_158[k]
                   + dk_198[k];

        t_159[k] = -ab_x[k] * di_159[k]
                   + dk_199[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, di_160, di_161, di_162, \
                         di_163, di_164, dk_200, dk_201, dk_202, dk_203, \
                         dk_204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_160[k] = -ab_x[k] * di_160[k]
                   + dk_200[k];

        t_161[k] = -ab_x[k] * di_161[k]
                   + dk_201[k];

        t_162[k] = -ab_x[k] * di_162[k]
                   + dk_202[k];

        t_163[k] = -ab_x[k] * di_163[k]
                   + dk_203[k];

        t_164[k] = -ab_x[k] * di_164[k]
                   + dk_204[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, ab_x, ab_y, di_84, di_165, di_166, \
                         di_167, dk_109, dk_205, dk_206, dk_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_165[k] = -ab_x[k] * di_165[k]
                   + dk_205[k];

        t_166[k] = -ab_x[k] * di_166[k]
                   + dk_206[k];

        t_167[k] = -ab_x[k] * di_167[k]
                   + dk_207[k];

        t_168[k] = -ab_y[k] * di_84[k]
                   + dk_109[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, ab_y, di_85, di_86, di_87, di_88, \
                         di_89, dk_111, dk_112, dk_114, dk_115, \
                         dk_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_169[k] = -ab_y[k] * di_85[k]
                   + dk_111[k];

        t_170[k] = -ab_y[k] * di_86[k]
                   + dk_112[k];

        t_171[k] = -ab_y[k] * di_87[k]
                   + dk_114[k];

        t_172[k] = -ab_y[k] * di_88[k]
                   + dk_115[k];

        t_173[k] = -ab_y[k] * di_89[k]
                   + dk_116[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, ab_y, di_90, di_91, di_92, di_93, \
                         di_94, dk_118, dk_119, dk_120, dk_121, \
                         dk_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_174[k] = -ab_y[k] * di_90[k]
                   + dk_118[k];

        t_175[k] = -ab_y[k] * di_91[k]
                   + dk_119[k];

        t_176[k] = -ab_y[k] * di_92[k]
                   + dk_120[k];

        t_177[k] = -ab_y[k] * di_93[k]
                   + dk_121[k];

        t_178[k] = -ab_y[k] * di_94[k]
                   + dk_123[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, ab_y, di_95, di_96, di_97, di_98, \
                         di_99, dk_124, dk_125, dk_126, dk_127, \
                         dk_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_179[k] = -ab_y[k] * di_95[k]
                   + dk_124[k];

        t_180[k] = -ab_y[k] * di_96[k]
                   + dk_125[k];

        t_181[k] = -ab_y[k] * di_97[k]
                   + dk_126[k];

        t_182[k] = -ab_y[k] * di_98[k]
                   + dk_127[k];

        t_183[k] = -ab_y[k] * di_99[k]
                   + dk_129[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, ab_y, di_100, di_101, di_102, \
                         di_103, di_104, dk_130, dk_131, dk_132, dk_133, \
                         dk_134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_184[k] = -ab_y[k] * di_100[k]
                   + dk_130[k];

        t_185[k] = -ab_y[k] * di_101[k]
                   + dk_131[k];

        t_186[k] = -ab_y[k] * di_102[k]
                   + dk_132[k];

        t_187[k] = -ab_y[k] * di_103[k]
                   + dk_133[k];

        t_188[k] = -ab_y[k] * di_104[k]
                   + dk_134[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, ab_y, di_105, di_106, di_107, \
                         di_108, di_109, dk_136, dk_137, dk_138, dk_139, \
                         dk_140 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_189[k] = -ab_y[k] * di_105[k]
                   + dk_136[k];

        t_190[k] = -ab_y[k] * di_106[k]
                   + dk_137[k];

        t_191[k] = -ab_y[k] * di_107[k]
                   + dk_138[k];

        t_192[k] = -ab_y[k] * di_108[k]
                   + dk_139[k];

        t_193[k] = -ab_y[k] * di_109[k]
                   + dk_140[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, ab_y, di_110, di_111, di_112, \
                         di_113, di_114, dk_141, dk_142, dk_145, dk_147, \
                         dk_148 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_194[k] = -ab_y[k] * di_110[k]
                   + dk_141[k];

        t_195[k] = -ab_y[k] * di_111[k]
                   + dk_142[k];

        t_196[k] = -ab_y[k] * di_112[k]
                   + dk_145[k];

        t_197[k] = -ab_y[k] * di_113[k]
                   + dk_147[k];

        t_198[k] = -ab_y[k] * di_114[k]
                   + dk_148[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, ab_y, di_115, di_116, di_117, \
                         di_118, di_119, dk_150, dk_151, dk_152, dk_154, \
                         dk_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_199[k] = -ab_y[k] * di_115[k]
                   + dk_150[k];

        t_200[k] = -ab_y[k] * di_116[k]
                   + dk_151[k];

        t_201[k] = -ab_y[k] * di_117[k]
                   + dk_152[k];

        t_202[k] = -ab_y[k] * di_118[k]
                   + dk_154[k];

        t_203[k] = -ab_y[k] * di_119[k]
                   + dk_155[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, ab_y, di_120, di_121, di_122, \
                         di_123, di_124, dk_156, dk_157, dk_159, dk_160, \
                         dk_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_204[k] = -ab_y[k] * di_120[k]
                   + dk_156[k];

        t_205[k] = -ab_y[k] * di_121[k]
                   + dk_157[k];

        t_206[k] = -ab_y[k] * di_122[k]
                   + dk_159[k];

        t_207[k] = -ab_y[k] * di_123[k]
                   + dk_160[k];

        t_208[k] = -ab_y[k] * di_124[k]
                   + dk_161[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, ab_y, di_125, di_126, di_127, \
                         di_128, di_129, dk_162, dk_163, dk_165, dk_166, \
                         dk_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_209[k] = -ab_y[k] * di_125[k]
                   + dk_162[k];

        t_210[k] = -ab_y[k] * di_126[k]
                   + dk_163[k];

        t_211[k] = -ab_y[k] * di_127[k]
                   + dk_165[k];

        t_212[k] = -ab_y[k] * di_128[k]
                   + dk_166[k];

        t_213[k] = -ab_y[k] * di_129[k]
                   + dk_167[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, ab_y, di_130, di_131, di_132, \
                         di_133, di_134, dk_168, dk_169, dk_170, dk_172, \
                         dk_173 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_214[k] = -ab_y[k] * di_130[k]
                   + dk_168[k];

        t_215[k] = -ab_y[k] * di_131[k]
                   + dk_169[k];

        t_216[k] = -ab_y[k] * di_132[k]
                   + dk_170[k];

        t_217[k] = -ab_y[k] * di_133[k]
                   + dk_172[k];

        t_218[k] = -ab_y[k] * di_134[k]
                   + dk_173[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, ab_y, di_135, di_136, di_137, \
                         di_138, di_139, dk_174, dk_175, dk_176, dk_177, \
                         dk_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_219[k] = -ab_y[k] * di_135[k]
                   + dk_174[k];

        t_220[k] = -ab_y[k] * di_136[k]
                   + dk_175[k];

        t_221[k] = -ab_y[k] * di_137[k]
                   + dk_176[k];

        t_222[k] = -ab_y[k] * di_138[k]
                   + dk_177[k];

        t_223[k] = -ab_y[k] * di_139[k]
                   + dk_178[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, ab_y, di_140, di_141, di_142, \
                         di_143, di_144, dk_181, dk_183, dk_184, dk_186, \
                         dk_187 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_224[k] = -ab_y[k] * di_140[k]
                   + dk_181[k];

        t_225[k] = -ab_y[k] * di_141[k]
                   + dk_183[k];

        t_226[k] = -ab_y[k] * di_142[k]
                   + dk_184[k];

        t_227[k] = -ab_y[k] * di_143[k]
                   + dk_186[k];

        t_228[k] = -ab_y[k] * di_144[k]
                   + dk_187[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, ab_y, di_145, di_146, di_147, \
                         di_148, di_149, dk_188, dk_190, dk_191, dk_192, \
                         dk_193 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_229[k] = -ab_y[k] * di_145[k]
                   + dk_188[k];

        t_230[k] = -ab_y[k] * di_146[k]
                   + dk_190[k];

        t_231[k] = -ab_y[k] * di_147[k]
                   + dk_191[k];

        t_232[k] = -ab_y[k] * di_148[k]
                   + dk_192[k];

        t_233[k] = -ab_y[k] * di_149[k]
                   + dk_193[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, ab_y, di_150, di_151, di_152, \
                         di_153, di_154, dk_195, dk_196, dk_197, dk_198, \
                         dk_199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_234[k] = -ab_y[k] * di_150[k]
                   + dk_195[k];

        t_235[k] = -ab_y[k] * di_151[k]
                   + dk_196[k];

        t_236[k] = -ab_y[k] * di_152[k]
                   + dk_197[k];

        t_237[k] = -ab_y[k] * di_153[k]
                   + dk_198[k];

        t_238[k] = -ab_y[k] * di_154[k]
                   + dk_199[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, ab_y, di_155, di_156, di_157, \
                         di_158, di_159, dk_201, dk_202, dk_203, dk_204, \
                         dk_205 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_239[k] = -ab_y[k] * di_155[k]
                   + dk_201[k];

        t_240[k] = -ab_y[k] * di_156[k]
                   + dk_202[k];

        t_241[k] = -ab_y[k] * di_157[k]
                   + dk_203[k];

        t_242[k] = -ab_y[k] * di_158[k]
                   + dk_204[k];

        t_243[k] = -ab_y[k] * di_159[k]
                   + dk_205[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, t_248, ab_y, di_160, di_161, di_162, \
                         di_163, di_164, dk_206, dk_208, dk_209, dk_210, \
                         dk_211 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_244[k] = -ab_y[k] * di_160[k]
                   + dk_206[k];

        t_245[k] = -ab_y[k] * di_161[k]
                   + dk_208[k];

        t_246[k] = -ab_y[k] * di_162[k]
                   + dk_209[k];

        t_247[k] = -ab_y[k] * di_163[k]
                   + dk_210[k];

        t_248[k] = -ab_y[k] * di_164[k]
                   + dk_211[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, ab_y, ab_z, di_140, di_165, di_166, \
                         di_167, dk_182, dk_212, dk_213, dk_214 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_249[k] = -ab_y[k] * di_165[k]
                   + dk_212[k];

        t_250[k] = -ab_y[k] * di_166[k]
                   + dk_213[k];

        t_251[k] = -ab_y[k] * di_167[k]
                   + dk_214[k];

        t_252[k] = -ab_z[k] * di_140[k]
                   + dk_182[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, t_257, ab_z, di_141, di_142, di_143, \
                         di_144, di_145, dk_184, dk_185, dk_187, dk_188, \
                         dk_189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_253[k] = -ab_z[k] * di_141[k]
                   + dk_184[k];

        t_254[k] = -ab_z[k] * di_142[k]
                   + dk_185[k];

        t_255[k] = -ab_z[k] * di_143[k]
                   + dk_187[k];

        t_256[k] = -ab_z[k] * di_144[k]
                   + dk_188[k];

        t_257[k] = -ab_z[k] * di_145[k]
                   + dk_189[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, ab_z, di_146, di_147, di_148, \
                         di_149, di_150, dk_191, dk_192, dk_193, dk_194, \
                         dk_196 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_258[k] = -ab_z[k] * di_146[k]
                   + dk_191[k];

        t_259[k] = -ab_z[k] * di_147[k]
                   + dk_192[k];

        t_260[k] = -ab_z[k] * di_148[k]
                   + dk_193[k];

        t_261[k] = -ab_z[k] * di_149[k]
                   + dk_194[k];

        t_262[k] = -ab_z[k] * di_150[k]
                   + dk_196[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, ab_z, di_151, di_152, di_153, \
                         di_154, di_155, dk_197, dk_198, dk_199, dk_200, \
                         dk_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_263[k] = -ab_z[k] * di_151[k]
                   + dk_197[k];

        t_264[k] = -ab_z[k] * di_152[k]
                   + dk_198[k];

        t_265[k] = -ab_z[k] * di_153[k]
                   + dk_199[k];

        t_266[k] = -ab_z[k] * di_154[k]
                   + dk_200[k];

        t_267[k] = -ab_z[k] * di_155[k]
                   + dk_202[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, ab_z, di_156, di_157, di_158, \
                         di_159, di_160, dk_203, dk_204, dk_205, dk_206, \
                         dk_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_268[k] = -ab_z[k] * di_156[k]
                   + dk_203[k];

        t_269[k] = -ab_z[k] * di_157[k]
                   + dk_204[k];

        t_270[k] = -ab_z[k] * di_158[k]
                   + dk_205[k];

        t_271[k] = -ab_z[k] * di_159[k]
                   + dk_206[k];

        t_272[k] = -ab_z[k] * di_160[k]
                   + dk_207[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, ab_z, di_161, di_162, di_163, \
                         di_164, di_165, dk_209, dk_210, dk_211, dk_212, \
                         dk_213 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_273[k] = -ab_z[k] * di_161[k]
                   + dk_209[k];

        t_274[k] = -ab_z[k] * di_162[k]
                   + dk_210[k];

        t_275[k] = -ab_z[k] * di_163[k]
                   + dk_211[k];

        t_276[k] = -ab_z[k] * di_164[k]
                   + dk_212[k];

        t_277[k] = -ab_z[k] * di_165[k]
                   + dk_213[k];
    }

#pragma omp simd aligned(t_278, t_279, ab_z, di_166, di_167, dk_214, \
                         dk_215 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_278[k] = -ab_z[k] * di_166[k]
                   + dk_214[k];

        t_279[k] = -ab_z[k] * di_167[k]
                   + dk_215[k];
    }
}

}  // namespace simdtrf
