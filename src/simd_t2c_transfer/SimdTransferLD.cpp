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


#include "SimdTransferLD.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_ld(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t lp, const size_t mp, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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
    const auto *lp_101 = buffer.data(lp + 101);
    const auto *lp_102 = buffer.data(lp + 102);
    const auto *lp_103 = buffer.data(lp + 103);
    const auto *lp_104 = buffer.data(lp + 104);
    const auto *lp_105 = buffer.data(lp + 105);
    const auto *lp_106 = buffer.data(lp + 106);
    const auto *lp_107 = buffer.data(lp + 107);
    const auto *lp_108 = buffer.data(lp + 108);
    const auto *lp_109 = buffer.data(lp + 109);
    const auto *lp_110 = buffer.data(lp + 110);
    const auto *lp_111 = buffer.data(lp + 111);
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
    const auto *lp_129 = buffer.data(lp + 129);
    const auto *lp_130 = buffer.data(lp + 130);
    const auto *lp_131 = buffer.data(lp + 131);
    const auto *lp_132 = buffer.data(lp + 132);
    const auto *lp_133 = buffer.data(lp + 133);
    const auto *lp_134 = buffer.data(lp + 134);

    const auto *mp_0 = buffer.data(mp + 0);
    const auto *mp_1 = buffer.data(mp + 1);
    const auto *mp_2 = buffer.data(mp + 2);
    const auto *mp_3 = buffer.data(mp + 3);
    const auto *mp_4 = buffer.data(mp + 4);
    const auto *mp_5 = buffer.data(mp + 5);
    const auto *mp_6 = buffer.data(mp + 6);
    const auto *mp_7 = buffer.data(mp + 7);
    const auto *mp_8 = buffer.data(mp + 8);
    const auto *mp_9 = buffer.data(mp + 9);
    const auto *mp_10 = buffer.data(mp + 10);
    const auto *mp_11 = buffer.data(mp + 11);
    const auto *mp_12 = buffer.data(mp + 12);
    const auto *mp_13 = buffer.data(mp + 13);
    const auto *mp_14 = buffer.data(mp + 14);
    const auto *mp_15 = buffer.data(mp + 15);
    const auto *mp_16 = buffer.data(mp + 16);
    const auto *mp_17 = buffer.data(mp + 17);
    const auto *mp_18 = buffer.data(mp + 18);
    const auto *mp_19 = buffer.data(mp + 19);
    const auto *mp_20 = buffer.data(mp + 20);
    const auto *mp_21 = buffer.data(mp + 21);
    const auto *mp_22 = buffer.data(mp + 22);
    const auto *mp_23 = buffer.data(mp + 23);
    const auto *mp_24 = buffer.data(mp + 24);
    const auto *mp_25 = buffer.data(mp + 25);
    const auto *mp_26 = buffer.data(mp + 26);
    const auto *mp_27 = buffer.data(mp + 27);
    const auto *mp_28 = buffer.data(mp + 28);
    const auto *mp_29 = buffer.data(mp + 29);
    const auto *mp_30 = buffer.data(mp + 30);
    const auto *mp_31 = buffer.data(mp + 31);
    const auto *mp_32 = buffer.data(mp + 32);
    const auto *mp_33 = buffer.data(mp + 33);
    const auto *mp_34 = buffer.data(mp + 34);
    const auto *mp_35 = buffer.data(mp + 35);
    const auto *mp_36 = buffer.data(mp + 36);
    const auto *mp_37 = buffer.data(mp + 37);
    const auto *mp_38 = buffer.data(mp + 38);
    const auto *mp_39 = buffer.data(mp + 39);
    const auto *mp_40 = buffer.data(mp + 40);
    const auto *mp_41 = buffer.data(mp + 41);
    const auto *mp_42 = buffer.data(mp + 42);
    const auto *mp_43 = buffer.data(mp + 43);
    const auto *mp_44 = buffer.data(mp + 44);
    const auto *mp_45 = buffer.data(mp + 45);
    const auto *mp_46 = buffer.data(mp + 46);
    const auto *mp_47 = buffer.data(mp + 47);
    const auto *mp_48 = buffer.data(mp + 48);
    const auto *mp_49 = buffer.data(mp + 49);
    const auto *mp_50 = buffer.data(mp + 50);
    const auto *mp_51 = buffer.data(mp + 51);
    const auto *mp_52 = buffer.data(mp + 52);
    const auto *mp_53 = buffer.data(mp + 53);
    const auto *mp_54 = buffer.data(mp + 54);
    const auto *mp_55 = buffer.data(mp + 55);
    const auto *mp_56 = buffer.data(mp + 56);
    const auto *mp_57 = buffer.data(mp + 57);
    const auto *mp_58 = buffer.data(mp + 58);
    const auto *mp_59 = buffer.data(mp + 59);
    const auto *mp_60 = buffer.data(mp + 60);
    const auto *mp_61 = buffer.data(mp + 61);
    const auto *mp_62 = buffer.data(mp + 62);
    const auto *mp_63 = buffer.data(mp + 63);
    const auto *mp_64 = buffer.data(mp + 64);
    const auto *mp_65 = buffer.data(mp + 65);
    const auto *mp_66 = buffer.data(mp + 66);
    const auto *mp_67 = buffer.data(mp + 67);
    const auto *mp_68 = buffer.data(mp + 68);
    const auto *mp_69 = buffer.data(mp + 69);
    const auto *mp_70 = buffer.data(mp + 70);
    const auto *mp_71 = buffer.data(mp + 71);
    const auto *mp_72 = buffer.data(mp + 72);
    const auto *mp_73 = buffer.data(mp + 73);
    const auto *mp_74 = buffer.data(mp + 74);
    const auto *mp_75 = buffer.data(mp + 75);
    const auto *mp_76 = buffer.data(mp + 76);
    const auto *mp_77 = buffer.data(mp + 77);
    const auto *mp_78 = buffer.data(mp + 78);
    const auto *mp_79 = buffer.data(mp + 79);
    const auto *mp_80 = buffer.data(mp + 80);
    const auto *mp_81 = buffer.data(mp + 81);
    const auto *mp_82 = buffer.data(mp + 82);
    const auto *mp_83 = buffer.data(mp + 83);
    const auto *mp_84 = buffer.data(mp + 84);
    const auto *mp_85 = buffer.data(mp + 85);
    const auto *mp_86 = buffer.data(mp + 86);
    const auto *mp_87 = buffer.data(mp + 87);
    const auto *mp_88 = buffer.data(mp + 88);
    const auto *mp_89 = buffer.data(mp + 89);
    const auto *mp_90 = buffer.data(mp + 90);
    const auto *mp_91 = buffer.data(mp + 91);
    const auto *mp_92 = buffer.data(mp + 92);
    const auto *mp_93 = buffer.data(mp + 93);
    const auto *mp_94 = buffer.data(mp + 94);
    const auto *mp_95 = buffer.data(mp + 95);
    const auto *mp_96 = buffer.data(mp + 96);
    const auto *mp_97 = buffer.data(mp + 97);
    const auto *mp_98 = buffer.data(mp + 98);
    const auto *mp_99 = buffer.data(mp + 99);
    const auto *mp_100 = buffer.data(mp + 100);
    const auto *mp_101 = buffer.data(mp + 101);
    const auto *mp_102 = buffer.data(mp + 102);
    const auto *mp_103 = buffer.data(mp + 103);
    const auto *mp_104 = buffer.data(mp + 104);
    const auto *mp_105 = buffer.data(mp + 105);
    const auto *mp_106 = buffer.data(mp + 106);
    const auto *mp_107 = buffer.data(mp + 107);
    const auto *mp_108 = buffer.data(mp + 108);
    const auto *mp_109 = buffer.data(mp + 109);
    const auto *mp_110 = buffer.data(mp + 110);
    const auto *mp_111 = buffer.data(mp + 111);
    const auto *mp_112 = buffer.data(mp + 112);
    const auto *mp_113 = buffer.data(mp + 113);
    const auto *mp_114 = buffer.data(mp + 114);
    const auto *mp_115 = buffer.data(mp + 115);
    const auto *mp_116 = buffer.data(mp + 116);
    const auto *mp_117 = buffer.data(mp + 117);
    const auto *mp_118 = buffer.data(mp + 118);
    const auto *mp_119 = buffer.data(mp + 119);
    const auto *mp_120 = buffer.data(mp + 120);
    const auto *mp_121 = buffer.data(mp + 121);
    const auto *mp_122 = buffer.data(mp + 122);
    const auto *mp_123 = buffer.data(mp + 123);
    const auto *mp_124 = buffer.data(mp + 124);
    const auto *mp_125 = buffer.data(mp + 125);
    const auto *mp_126 = buffer.data(mp + 126);
    const auto *mp_127 = buffer.data(mp + 127);
    const auto *mp_128 = buffer.data(mp + 128);
    const auto *mp_129 = buffer.data(mp + 129);
    const auto *mp_130 = buffer.data(mp + 130);
    const auto *mp_131 = buffer.data(mp + 131);
    const auto *mp_132 = buffer.data(mp + 132);
    const auto *mp_133 = buffer.data(mp + 133);
    const auto *mp_134 = buffer.data(mp + 134);
    const auto *mp_136 = buffer.data(mp + 136);
    const auto *mp_137 = buffer.data(mp + 137);
    const auto *mp_139 = buffer.data(mp + 139);
    const auto *mp_140 = buffer.data(mp + 140);
    const auto *mp_142 = buffer.data(mp + 142);
    const auto *mp_143 = buffer.data(mp + 143);
    const auto *mp_145 = buffer.data(mp + 145);
    const auto *mp_146 = buffer.data(mp + 146);
    const auto *mp_148 = buffer.data(mp + 148);
    const auto *mp_149 = buffer.data(mp + 149);
    const auto *mp_151 = buffer.data(mp + 151);
    const auto *mp_152 = buffer.data(mp + 152);
    const auto *mp_154 = buffer.data(mp + 154);
    const auto *mp_155 = buffer.data(mp + 155);
    const auto *mp_157 = buffer.data(mp + 157);
    const auto *mp_158 = buffer.data(mp + 158);
    const auto *mp_160 = buffer.data(mp + 160);
    const auto *mp_161 = buffer.data(mp + 161);
    const auto *mp_164 = buffer.data(mp + 164);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, lp_0, lp_1, lp_2, mp_0, mp_1, \
                         mp_2, mp_4, mp_5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * lp_0[k]
                 + mp_0[k];

        t_1[k] = ab_x[k] * lp_1[k]
                 + mp_1[k];

        t_2[k] = ab_x[k] * lp_2[k]
                 + mp_2[k];

        t_3[k] = ab_y[k] * lp_1[k]
                 + mp_4[k];

        t_4[k] = ab_y[k] * lp_2[k]
                 + mp_5[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_z, lp_2, lp_3, lp_4, lp_5, mp_3, mp_4, \
                         mp_5, mp_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = ab_z[k] * lp_2[k]
                 + mp_8[k];

        t_6[k] = ab_x[k] * lp_3[k]
                 + mp_3[k];

        t_7[k] = ab_x[k] * lp_4[k]
                 + mp_4[k];

        t_8[k] = ab_x[k] * lp_5[k]
                 + mp_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_y, ab_z, lp_4, lp_5, lp_6, mp_6, \
                         mp_10, mp_11, mp_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_9[k] = ab_y[k] * lp_4[k]
                 + mp_10[k];

        t_10[k] = ab_y[k] * lp_5[k]
                  + mp_11[k];

        t_11[k] = ab_z[k] * lp_5[k]
                  + mp_14[k];

        t_12[k] = ab_x[k] * lp_6[k]
                  + mp_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, lp_7, lp_8, mp_7, \
                         mp_8, mp_13, mp_14, mp_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_13[k] = ab_x[k] * lp_7[k]
                  + mp_7[k];

        t_14[k] = ab_x[k] * lp_8[k]
                  + mp_8[k];

        t_15[k] = ab_y[k] * lp_7[k]
                  + mp_13[k];

        t_16[k] = ab_y[k] * lp_8[k]
                  + mp_14[k];

        t_17[k] = ab_z[k] * lp_8[k]
                  + mp_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, lp_9, lp_10, lp_11, mp_9, \
                         mp_10, mp_11, mp_19, mp_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_18[k] = ab_x[k] * lp_9[k]
                  + mp_9[k];

        t_19[k] = ab_x[k] * lp_10[k]
                  + mp_10[k];

        t_20[k] = ab_x[k] * lp_11[k]
                  + mp_11[k];

        t_21[k] = ab_y[k] * lp_10[k]
                  + mp_19[k];

        t_22[k] = ab_y[k] * lp_11[k]
                  + mp_20[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_z, lp_11, lp_12, lp_13, lp_14, \
                         mp_12, mp_13, mp_14, mp_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_23[k] = ab_z[k] * lp_11[k]
                  + mp_23[k];

        t_24[k] = ab_x[k] * lp_12[k]
                  + mp_12[k];

        t_25[k] = ab_x[k] * lp_13[k]
                  + mp_13[k];

        t_26[k] = ab_x[k] * lp_14[k]
                  + mp_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, lp_13, lp_14, lp_15, mp_15, \
                         mp_22, mp_23, mp_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_27[k] = ab_y[k] * lp_13[k]
                  + mp_22[k];

        t_28[k] = ab_y[k] * lp_14[k]
                  + mp_23[k];

        t_29[k] = ab_z[k] * lp_14[k]
                  + mp_26[k];

        t_30[k] = ab_x[k] * lp_15[k]
                  + mp_15[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, lp_16, lp_17, mp_16, \
                         mp_17, mp_25, mp_26, mp_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_31[k] = ab_x[k] * lp_16[k]
                  + mp_16[k];

        t_32[k] = ab_x[k] * lp_17[k]
                  + mp_17[k];

        t_33[k] = ab_y[k] * lp_16[k]
                  + mp_25[k];

        t_34[k] = ab_y[k] * lp_17[k]
                  + mp_26[k];

        t_35[k] = ab_z[k] * lp_17[k]
                  + mp_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_x, ab_y, lp_18, lp_19, lp_20, mp_18, \
                         mp_19, mp_20, mp_31, mp_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_36[k] = ab_x[k] * lp_18[k]
                  + mp_18[k];

        t_37[k] = ab_x[k] * lp_19[k]
                  + mp_19[k];

        t_38[k] = ab_x[k] * lp_20[k]
                  + mp_20[k];

        t_39[k] = ab_y[k] * lp_19[k]
                  + mp_31[k];

        t_40[k] = ab_y[k] * lp_20[k]
                  + mp_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_x, ab_z, lp_20, lp_21, lp_22, lp_23, \
                         mp_21, mp_22, mp_23, mp_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_41[k] = ab_z[k] * lp_20[k]
                  + mp_35[k];

        t_42[k] = ab_x[k] * lp_21[k]
                  + mp_21[k];

        t_43[k] = ab_x[k] * lp_22[k]
                  + mp_22[k];

        t_44[k] = ab_x[k] * lp_23[k]
                  + mp_23[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, lp_22, lp_23, lp_24, mp_24, \
                         mp_34, mp_35, mp_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = ab_y[k] * lp_22[k]
                  + mp_34[k];

        t_46[k] = ab_y[k] * lp_23[k]
                  + mp_35[k];

        t_47[k] = ab_z[k] * lp_23[k]
                  + mp_38[k];

        t_48[k] = ab_x[k] * lp_24[k]
                  + mp_24[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, lp_25, lp_26, mp_25, \
                         mp_26, mp_37, mp_38, mp_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_49[k] = ab_x[k] * lp_25[k]
                  + mp_25[k];

        t_50[k] = ab_x[k] * lp_26[k]
                  + mp_26[k];

        t_51[k] = ab_y[k] * lp_25[k]
                  + mp_37[k];

        t_52[k] = ab_y[k] * lp_26[k]
                  + mp_38[k];

        t_53[k] = ab_z[k] * lp_26[k]
                  + mp_41[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_x, ab_y, lp_27, lp_28, lp_29, mp_27, \
                         mp_28, mp_29, mp_40, mp_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_54[k] = ab_x[k] * lp_27[k]
                  + mp_27[k];

        t_55[k] = ab_x[k] * lp_28[k]
                  + mp_28[k];

        t_56[k] = ab_x[k] * lp_29[k]
                  + mp_29[k];

        t_57[k] = ab_y[k] * lp_28[k]
                  + mp_40[k];

        t_58[k] = ab_y[k] * lp_29[k]
                  + mp_41[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, ab_x, ab_z, lp_29, lp_30, lp_31, lp_32, \
                         mp_30, mp_31, mp_32, mp_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_59[k] = ab_z[k] * lp_29[k]
                  + mp_44[k];

        t_60[k] = ab_x[k] * lp_30[k]
                  + mp_30[k];

        t_61[k] = ab_x[k] * lp_31[k]
                  + mp_31[k];

        t_62[k] = ab_x[k] * lp_32[k]
                  + mp_32[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, ab_x, ab_y, ab_z, lp_31, lp_32, lp_33, mp_33, \
                         mp_46, mp_47, mp_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_63[k] = ab_y[k] * lp_31[k]
                  + mp_46[k];

        t_64[k] = ab_y[k] * lp_32[k]
                  + mp_47[k];

        t_65[k] = ab_z[k] * lp_32[k]
                  + mp_50[k];

        t_66[k] = ab_x[k] * lp_33[k]
                  + mp_33[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, ab_x, ab_y, ab_z, lp_34, lp_35, mp_34, \
                         mp_35, mp_49, mp_50, mp_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_67[k] = ab_x[k] * lp_34[k]
                  + mp_34[k];

        t_68[k] = ab_x[k] * lp_35[k]
                  + mp_35[k];

        t_69[k] = ab_y[k] * lp_34[k]
                  + mp_49[k];

        t_70[k] = ab_y[k] * lp_35[k]
                  + mp_50[k];

        t_71[k] = ab_z[k] * lp_35[k]
                  + mp_53[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, lp_36, lp_37, lp_38, mp_36, \
                         mp_37, mp_38, mp_52, mp_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_72[k] = ab_x[k] * lp_36[k]
                  + mp_36[k];

        t_73[k] = ab_x[k] * lp_37[k]
                  + mp_37[k];

        t_74[k] = ab_x[k] * lp_38[k]
                  + mp_38[k];

        t_75[k] = ab_y[k] * lp_37[k]
                  + mp_52[k];

        t_76[k] = ab_y[k] * lp_38[k]
                  + mp_53[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, ab_x, ab_z, lp_38, lp_39, lp_40, lp_41, \
                         mp_39, mp_40, mp_41, mp_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_77[k] = ab_z[k] * lp_38[k]
                  + mp_56[k];

        t_78[k] = ab_x[k] * lp_39[k]
                  + mp_39[k];

        t_79[k] = ab_x[k] * lp_40[k]
                  + mp_40[k];

        t_80[k] = ab_x[k] * lp_41[k]
                  + mp_41[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, ab_x, ab_y, ab_z, lp_40, lp_41, lp_42, mp_42, \
                         mp_55, mp_56, mp_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_81[k] = ab_y[k] * lp_40[k]
                  + mp_55[k];

        t_82[k] = ab_y[k] * lp_41[k]
                  + mp_56[k];

        t_83[k] = ab_z[k] * lp_41[k]
                  + mp_59[k];

        t_84[k] = ab_x[k] * lp_42[k]
                  + mp_42[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, lp_43, lp_44, mp_43, \
                         mp_44, mp_58, mp_59, mp_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = ab_x[k] * lp_43[k]
                  + mp_43[k];

        t_86[k] = ab_x[k] * lp_44[k]
                  + mp_44[k];

        t_87[k] = ab_y[k] * lp_43[k]
                  + mp_58[k];

        t_88[k] = ab_y[k] * lp_44[k]
                  + mp_59[k];

        t_89[k] = ab_z[k] * lp_44[k]
                  + mp_62[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, ab_y, lp_45, lp_46, lp_47, mp_45, \
                         mp_46, mp_47, mp_64, mp_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = ab_x[k] * lp_45[k]
                  + mp_45[k];

        t_91[k] = ab_x[k] * lp_46[k]
                  + mp_46[k];

        t_92[k] = ab_x[k] * lp_47[k]
                  + mp_47[k];

        t_93[k] = ab_y[k] * lp_46[k]
                  + mp_64[k];

        t_94[k] = ab_y[k] * lp_47[k]
                  + mp_65[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, ab_x, ab_z, lp_47, lp_48, lp_49, lp_50, \
                         mp_48, mp_49, mp_50, mp_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = ab_z[k] * lp_47[k]
                  + mp_68[k];

        t_96[k] = ab_x[k] * lp_48[k]
                  + mp_48[k];

        t_97[k] = ab_x[k] * lp_49[k]
                  + mp_49[k];

        t_98[k] = ab_x[k] * lp_50[k]
                  + mp_50[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, ab_x, ab_y, ab_z, lp_49, lp_50, lp_51, \
                         mp_51, mp_67, mp_68, mp_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_99[k] = ab_y[k] * lp_49[k]
                  + mp_67[k];

        t_100[k] = ab_y[k] * lp_50[k]
                   + mp_68[k];

        t_101[k] = ab_z[k] * lp_50[k]
                   + mp_71[k];

        t_102[k] = ab_x[k] * lp_51[k]
                   + mp_51[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, ab_x, ab_y, ab_z, lp_52, lp_53, \
                         mp_52, mp_53, mp_70, mp_71, mp_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_103[k] = ab_x[k] * lp_52[k]
                   + mp_52[k];

        t_104[k] = ab_x[k] * lp_53[k]
                   + mp_53[k];

        t_105[k] = ab_y[k] * lp_52[k]
                   + mp_70[k];

        t_106[k] = ab_y[k] * lp_53[k]
                   + mp_71[k];

        t_107[k] = ab_z[k] * lp_53[k]
                   + mp_74[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, ab_x, ab_y, lp_54, lp_55, lp_56, \
                         mp_54, mp_55, mp_56, mp_73, mp_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_108[k] = ab_x[k] * lp_54[k]
                   + mp_54[k];

        t_109[k] = ab_x[k] * lp_55[k]
                   + mp_55[k];

        t_110[k] = ab_x[k] * lp_56[k]
                   + mp_56[k];

        t_111[k] = ab_y[k] * lp_55[k]
                   + mp_73[k];

        t_112[k] = ab_y[k] * lp_56[k]
                   + mp_74[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, ab_x, ab_z, lp_56, lp_57, lp_58, lp_59, \
                         mp_57, mp_58, mp_59, mp_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_113[k] = ab_z[k] * lp_56[k]
                   + mp_77[k];

        t_114[k] = ab_x[k] * lp_57[k]
                   + mp_57[k];

        t_115[k] = ab_x[k] * lp_58[k]
                   + mp_58[k];

        t_116[k] = ab_x[k] * lp_59[k]
                   + mp_59[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, ab_x, ab_y, ab_z, lp_58, lp_59, lp_60, \
                         mp_60, mp_76, mp_77, mp_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_117[k] = ab_y[k] * lp_58[k]
                   + mp_76[k];

        t_118[k] = ab_y[k] * lp_59[k]
                   + mp_77[k];

        t_119[k] = ab_z[k] * lp_59[k]
                   + mp_80[k];

        t_120[k] = ab_x[k] * lp_60[k]
                   + mp_60[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, ab_x, ab_y, ab_z, lp_61, lp_62, \
                         mp_61, mp_62, mp_79, mp_80, mp_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_121[k] = ab_x[k] * lp_61[k]
                   + mp_61[k];

        t_122[k] = ab_x[k] * lp_62[k]
                   + mp_62[k];

        t_123[k] = ab_y[k] * lp_61[k]
                   + mp_79[k];

        t_124[k] = ab_y[k] * lp_62[k]
                   + mp_80[k];

        t_125[k] = ab_z[k] * lp_62[k]
                   + mp_83[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, ab_x, ab_y, lp_63, lp_64, lp_65, \
                         mp_63, mp_64, mp_65, mp_85, mp_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_126[k] = ab_x[k] * lp_63[k]
                   + mp_63[k];

        t_127[k] = ab_x[k] * lp_64[k]
                   + mp_64[k];

        t_128[k] = ab_x[k] * lp_65[k]
                   + mp_65[k];

        t_129[k] = ab_y[k] * lp_64[k]
                   + mp_85[k];

        t_130[k] = ab_y[k] * lp_65[k]
                   + mp_86[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, ab_x, ab_z, lp_65, lp_66, lp_67, lp_68, \
                         mp_66, mp_67, mp_68, mp_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_131[k] = ab_z[k] * lp_65[k]
                   + mp_89[k];

        t_132[k] = ab_x[k] * lp_66[k]
                   + mp_66[k];

        t_133[k] = ab_x[k] * lp_67[k]
                   + mp_67[k];

        t_134[k] = ab_x[k] * lp_68[k]
                   + mp_68[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, ab_x, ab_y, ab_z, lp_67, lp_68, lp_69, \
                         mp_69, mp_88, mp_89, mp_92 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = ab_y[k] * lp_67[k]
                   + mp_88[k];

        t_136[k] = ab_y[k] * lp_68[k]
                   + mp_89[k];

        t_137[k] = ab_z[k] * lp_68[k]
                   + mp_92[k];

        t_138[k] = ab_x[k] * lp_69[k]
                   + mp_69[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, ab_x, ab_y, ab_z, lp_70, lp_71, \
                         mp_70, mp_71, mp_91, mp_92, mp_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_139[k] = ab_x[k] * lp_70[k]
                   + mp_70[k];

        t_140[k] = ab_x[k] * lp_71[k]
                   + mp_71[k];

        t_141[k] = ab_y[k] * lp_70[k]
                   + mp_91[k];

        t_142[k] = ab_y[k] * lp_71[k]
                   + mp_92[k];

        t_143[k] = ab_z[k] * lp_71[k]
                   + mp_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, ab_x, ab_y, lp_72, lp_73, lp_74, \
                         mp_72, mp_73, mp_74, mp_94, mp_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_144[k] = ab_x[k] * lp_72[k]
                   + mp_72[k];

        t_145[k] = ab_x[k] * lp_73[k]
                   + mp_73[k];

        t_146[k] = ab_x[k] * lp_74[k]
                   + mp_74[k];

        t_147[k] = ab_y[k] * lp_73[k]
                   + mp_94[k];

        t_148[k] = ab_y[k] * lp_74[k]
                   + mp_95[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, ab_x, ab_z, lp_74, lp_75, lp_76, lp_77, \
                         mp_75, mp_76, mp_77, mp_98 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_149[k] = ab_z[k] * lp_74[k]
                   + mp_98[k];

        t_150[k] = ab_x[k] * lp_75[k]
                   + mp_75[k];

        t_151[k] = ab_x[k] * lp_76[k]
                   + mp_76[k];

        t_152[k] = ab_x[k] * lp_77[k]
                   + mp_77[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, ab_x, ab_y, ab_z, lp_76, lp_77, lp_78, \
                         mp_78, mp_97, mp_98, mp_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_153[k] = ab_y[k] * lp_76[k]
                   + mp_97[k];

        t_154[k] = ab_y[k] * lp_77[k]
                   + mp_98[k];

        t_155[k] = ab_z[k] * lp_77[k]
                   + mp_101[k];

        t_156[k] = ab_x[k] * lp_78[k]
                   + mp_78[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, ab_x, ab_y, ab_z, lp_79, lp_80, \
                         mp_79, mp_80, mp_100, mp_101, mp_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_157[k] = ab_x[k] * lp_79[k]
                   + mp_79[k];

        t_158[k] = ab_x[k] * lp_80[k]
                   + mp_80[k];

        t_159[k] = ab_y[k] * lp_79[k]
                   + mp_100[k];

        t_160[k] = ab_y[k] * lp_80[k]
                   + mp_101[k];

        t_161[k] = ab_z[k] * lp_80[k]
                   + mp_104[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, ab_x, ab_y, lp_81, lp_82, lp_83, \
                         mp_81, mp_82, mp_83, mp_103, mp_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_162[k] = ab_x[k] * lp_81[k]
                   + mp_81[k];

        t_163[k] = ab_x[k] * lp_82[k]
                   + mp_82[k];

        t_164[k] = ab_x[k] * lp_83[k]
                   + mp_83[k];

        t_165[k] = ab_y[k] * lp_82[k]
                   + mp_103[k];

        t_166[k] = ab_y[k] * lp_83[k]
                   + mp_104[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, ab_x, ab_z, lp_83, lp_84, lp_85, lp_86, \
                         mp_84, mp_85, mp_86, mp_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_167[k] = ab_z[k] * lp_83[k]
                   + mp_107[k];

        t_168[k] = ab_x[k] * lp_84[k]
                   + mp_84[k];

        t_169[k] = ab_x[k] * lp_85[k]
                   + mp_85[k];

        t_170[k] = ab_x[k] * lp_86[k]
                   + mp_86[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, ab_x, ab_y, ab_z, lp_85, lp_86, lp_87, \
                         mp_87, mp_109, mp_110, mp_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_171[k] = ab_y[k] * lp_85[k]
                   + mp_109[k];

        t_172[k] = ab_y[k] * lp_86[k]
                   + mp_110[k];

        t_173[k] = ab_z[k] * lp_86[k]
                   + mp_113[k];

        t_174[k] = ab_x[k] * lp_87[k]
                   + mp_87[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, ab_y, ab_z, lp_88, lp_89, \
                         mp_88, mp_89, mp_112, mp_113, mp_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_175[k] = ab_x[k] * lp_88[k]
                   + mp_88[k];

        t_176[k] = ab_x[k] * lp_89[k]
                   + mp_89[k];

        t_177[k] = ab_y[k] * lp_88[k]
                   + mp_112[k];

        t_178[k] = ab_y[k] * lp_89[k]
                   + mp_113[k];

        t_179[k] = ab_z[k] * lp_89[k]
                   + mp_116[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, ab_y, lp_90, lp_91, lp_92, \
                         mp_90, mp_91, mp_92, mp_115, mp_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_180[k] = ab_x[k] * lp_90[k]
                   + mp_90[k];

        t_181[k] = ab_x[k] * lp_91[k]
                   + mp_91[k];

        t_182[k] = ab_x[k] * lp_92[k]
                   + mp_92[k];

        t_183[k] = ab_y[k] * lp_91[k]
                   + mp_115[k];

        t_184[k] = ab_y[k] * lp_92[k]
                   + mp_116[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, ab_x, ab_z, lp_92, lp_93, lp_94, lp_95, \
                         mp_93, mp_94, mp_95, mp_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_185[k] = ab_z[k] * lp_92[k]
                   + mp_119[k];

        t_186[k] = ab_x[k] * lp_93[k]
                   + mp_93[k];

        t_187[k] = ab_x[k] * lp_94[k]
                   + mp_94[k];

        t_188[k] = ab_x[k] * lp_95[k]
                   + mp_95[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, ab_x, ab_y, ab_z, lp_94, lp_95, lp_96, \
                         mp_96, mp_118, mp_119, mp_122 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_189[k] = ab_y[k] * lp_94[k]
                   + mp_118[k];

        t_190[k] = ab_y[k] * lp_95[k]
                   + mp_119[k];

        t_191[k] = ab_z[k] * lp_95[k]
                   + mp_122[k];

        t_192[k] = ab_x[k] * lp_96[k]
                   + mp_96[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, ab_x, ab_y, ab_z, lp_97, lp_98, \
                         mp_97, mp_98, mp_121, mp_122, mp_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_193[k] = ab_x[k] * lp_97[k]
                   + mp_97[k];

        t_194[k] = ab_x[k] * lp_98[k]
                   + mp_98[k];

        t_195[k] = ab_y[k] * lp_97[k]
                   + mp_121[k];

        t_196[k] = ab_y[k] * lp_98[k]
                   + mp_122[k];

        t_197[k] = ab_z[k] * lp_98[k]
                   + mp_125[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, ab_x, ab_y, lp_99, lp_100, lp_101, \
                         mp_99, mp_100, mp_101, mp_124, mp_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_198[k] = ab_x[k] * lp_99[k]
                   + mp_99[k];

        t_199[k] = ab_x[k] * lp_100[k]
                   + mp_100[k];

        t_200[k] = ab_x[k] * lp_101[k]
                   + mp_101[k];

        t_201[k] = ab_y[k] * lp_100[k]
                   + mp_124[k];

        t_202[k] = ab_y[k] * lp_101[k]
                   + mp_125[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, ab_x, ab_z, lp_101, lp_102, lp_103, \
                         lp_104, mp_102, mp_103, mp_104, mp_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_203[k] = ab_z[k] * lp_101[k]
                   + mp_128[k];

        t_204[k] = ab_x[k] * lp_102[k]
                   + mp_102[k];

        t_205[k] = ab_x[k] * lp_103[k]
                   + mp_103[k];

        t_206[k] = ab_x[k] * lp_104[k]
                   + mp_104[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, ab_x, ab_y, ab_z, lp_103, lp_104, lp_105, \
                         mp_105, mp_127, mp_128, mp_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_207[k] = ab_y[k] * lp_103[k]
                   + mp_127[k];

        t_208[k] = ab_y[k] * lp_104[k]
                   + mp_128[k];

        t_209[k] = ab_z[k] * lp_104[k]
                   + mp_131[k];

        t_210[k] = ab_x[k] * lp_105[k]
                   + mp_105[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, ab_x, ab_y, ab_z, lp_106, lp_107, \
                         mp_106, mp_107, mp_130, mp_131, mp_134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_211[k] = ab_x[k] * lp_106[k]
                   + mp_106[k];

        t_212[k] = ab_x[k] * lp_107[k]
                   + mp_107[k];

        t_213[k] = ab_y[k] * lp_106[k]
                   + mp_130[k];

        t_214[k] = ab_y[k] * lp_107[k]
                   + mp_131[k];

        t_215[k] = ab_z[k] * lp_107[k]
                   + mp_134[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, ab_x, ab_y, lp_108, lp_109, \
                         lp_110, mp_108, mp_109, mp_110, mp_136, \
                         mp_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_216[k] = ab_x[k] * lp_108[k]
                   + mp_108[k];

        t_217[k] = ab_x[k] * lp_109[k]
                   + mp_109[k];

        t_218[k] = ab_x[k] * lp_110[k]
                   + mp_110[k];

        t_219[k] = ab_y[k] * lp_109[k]
                   + mp_136[k];

        t_220[k] = ab_y[k] * lp_110[k]
                   + mp_137[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, ab_x, ab_z, lp_110, lp_111, lp_112, \
                         lp_113, mp_111, mp_112, mp_113, mp_140 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_221[k] = ab_z[k] * lp_110[k]
                   + mp_140[k];

        t_222[k] = ab_x[k] * lp_111[k]
                   + mp_111[k];

        t_223[k] = ab_x[k] * lp_112[k]
                   + mp_112[k];

        t_224[k] = ab_x[k] * lp_113[k]
                   + mp_113[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, ab_x, ab_y, ab_z, lp_112, lp_113, lp_114, \
                         mp_114, mp_139, mp_140, mp_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_225[k] = ab_y[k] * lp_112[k]
                   + mp_139[k];

        t_226[k] = ab_y[k] * lp_113[k]
                   + mp_140[k];

        t_227[k] = ab_z[k] * lp_113[k]
                   + mp_143[k];

        t_228[k] = ab_x[k] * lp_114[k]
                   + mp_114[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, ab_x, ab_y, ab_z, lp_115, lp_116, \
                         mp_115, mp_116, mp_142, mp_143, mp_146 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_229[k] = ab_x[k] * lp_115[k]
                   + mp_115[k];

        t_230[k] = ab_x[k] * lp_116[k]
                   + mp_116[k];

        t_231[k] = ab_y[k] * lp_115[k]
                   + mp_142[k];

        t_232[k] = ab_y[k] * lp_116[k]
                   + mp_143[k];

        t_233[k] = ab_z[k] * lp_116[k]
                   + mp_146[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, ab_x, ab_y, lp_117, lp_118, \
                         lp_119, mp_117, mp_118, mp_119, mp_145, \
                         mp_146 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_234[k] = ab_x[k] * lp_117[k]
                   + mp_117[k];

        t_235[k] = ab_x[k] * lp_118[k]
                   + mp_118[k];

        t_236[k] = ab_x[k] * lp_119[k]
                   + mp_119[k];

        t_237[k] = ab_y[k] * lp_118[k]
                   + mp_145[k];

        t_238[k] = ab_y[k] * lp_119[k]
                   + mp_146[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, ab_x, ab_z, lp_119, lp_120, lp_121, \
                         lp_122, mp_120, mp_121, mp_122, mp_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_239[k] = ab_z[k] * lp_119[k]
                   + mp_149[k];

        t_240[k] = ab_x[k] * lp_120[k]
                   + mp_120[k];

        t_241[k] = ab_x[k] * lp_121[k]
                   + mp_121[k];

        t_242[k] = ab_x[k] * lp_122[k]
                   + mp_122[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, ab_x, ab_y, ab_z, lp_121, lp_122, lp_123, \
                         mp_123, mp_148, mp_149, mp_152 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_243[k] = ab_y[k] * lp_121[k]
                   + mp_148[k];

        t_244[k] = ab_y[k] * lp_122[k]
                   + mp_149[k];

        t_245[k] = ab_z[k] * lp_122[k]
                   + mp_152[k];

        t_246[k] = ab_x[k] * lp_123[k]
                   + mp_123[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, ab_x, ab_y, ab_z, lp_124, lp_125, \
                         mp_124, mp_125, mp_151, mp_152, mp_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_247[k] = ab_x[k] * lp_124[k]
                   + mp_124[k];

        t_248[k] = ab_x[k] * lp_125[k]
                   + mp_125[k];

        t_249[k] = ab_y[k] * lp_124[k]
                   + mp_151[k];

        t_250[k] = ab_y[k] * lp_125[k]
                   + mp_152[k];

        t_251[k] = ab_z[k] * lp_125[k]
                   + mp_155[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, ab_x, ab_y, lp_126, lp_127, \
                         lp_128, mp_126, mp_127, mp_128, mp_154, \
                         mp_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_252[k] = ab_x[k] * lp_126[k]
                   + mp_126[k];

        t_253[k] = ab_x[k] * lp_127[k]
                   + mp_127[k];

        t_254[k] = ab_x[k] * lp_128[k]
                   + mp_128[k];

        t_255[k] = ab_y[k] * lp_127[k]
                   + mp_154[k];

        t_256[k] = ab_y[k] * lp_128[k]
                   + mp_155[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, ab_x, ab_z, lp_128, lp_129, lp_130, \
                         lp_131, mp_129, mp_130, mp_131, mp_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_257[k] = ab_z[k] * lp_128[k]
                   + mp_158[k];

        t_258[k] = ab_x[k] * lp_129[k]
                   + mp_129[k];

        t_259[k] = ab_x[k] * lp_130[k]
                   + mp_130[k];

        t_260[k] = ab_x[k] * lp_131[k]
                   + mp_131[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, ab_x, ab_y, ab_z, lp_130, lp_131, lp_132, \
                         mp_132, mp_157, mp_158, mp_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_261[k] = ab_y[k] * lp_130[k]
                   + mp_157[k];

        t_262[k] = ab_y[k] * lp_131[k]
                   + mp_158[k];

        t_263[k] = ab_z[k] * lp_131[k]
                   + mp_161[k];

        t_264[k] = ab_x[k] * lp_132[k]
                   + mp_132[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, ab_y, ab_z, lp_133, lp_134, \
                         mp_133, mp_134, mp_160, mp_161, mp_164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_265[k] = ab_x[k] * lp_133[k]
                   + mp_133[k];

        t_266[k] = ab_x[k] * lp_134[k]
                   + mp_134[k];

        t_267[k] = ab_y[k] * lp_133[k]
                   + mp_160[k];

        t_268[k] = ab_y[k] * lp_134[k]
                   + mp_161[k];

        t_269[k] = ab_z[k] * lp_134[k]
                   + mp_164[k];
    }
}

}  // namespace simdtrf
