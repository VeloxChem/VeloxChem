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


#include "SimdElectronRepulsionVrrRecGG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_gg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg0, const size_t dg1,
                                     const size_t ff, const size_t fg, const size_t gd0,
                                     const size_t gd1, const size_t gf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_1 = buffer.data(dg1 + 1);
    const auto *dg1_2 = buffer.data(dg1 + 2);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_41 = buffer.data(ff + 41);
    const auto *ff_42 = buffer.data(ff + 42);
    const auto *ff_43 = buffer.data(ff + 43);
    const auto *ff_44 = buffer.data(ff + 44);
    const auto *ff_45 = buffer.data(ff + 45);
    const auto *ff_46 = buffer.data(ff + 46);
    const auto *ff_47 = buffer.data(ff + 47);
    const auto *ff_48 = buffer.data(ff + 48);
    const auto *ff_49 = buffer.data(ff + 49);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_51 = buffer.data(ff + 51);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_53 = buffer.data(ff + 53);
    const auto *ff_54 = buffer.data(ff + 54);
    const auto *ff_55 = buffer.data(ff + 55);
    const auto *ff_56 = buffer.data(ff + 56);
    const auto *ff_57 = buffer.data(ff + 57);
    const auto *ff_58 = buffer.data(ff + 58);
    const auto *ff_59 = buffer.data(ff + 59);
    const auto *ff_60 = buffer.data(ff + 60);
    const auto *ff_61 = buffer.data(ff + 61);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
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
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
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
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_17 = buffer.data(gd0 + 17);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_17 = buffer.data(gd1 + 17);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_62 = buffer.data(gf + 62);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_64 = buffer.data(gf + 64);
    const auto *gf_65 = buffer.data(gf + 65);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_67 = buffer.data(gf + 67);
    const auto *gf_68 = buffer.data(gf + 68);
    const auto *gf_69 = buffer.data(gf + 69);
    const auto *gf_70 = buffer.data(gf + 70);
    const auto *gf_71 = buffer.data(gf + 71);
    const auto *gf_72 = buffer.data(gf + 72);
    const auto *gf_73 = buffer.data(gf + 73);
    const auto *gf_74 = buffer.data(gf + 74);
    const auto *gf_75 = buffer.data(gf + 75);
    const auto *gf_76 = buffer.data(gf + 76);
    const auto *gf_77 = buffer.data(gf + 77);
    const auto *gf_78 = buffer.data(gf + 78);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_80 = buffer.data(gf + 80);
    const auto *gf_81 = buffer.data(gf + 81);
    const auto *gf_82 = buffer.data(gf + 82);
    const auto *gf_83 = buffer.data(gf + 83);
    const auto *gf_84 = buffer.data(gf + 84);
    const auto *gf_85 = buffer.data(gf + 85);
    const auto *gf_86 = buffer.data(gf + 86);
    const auto *gf_87 = buffer.data(gf + 87);
    const auto *gf_88 = buffer.data(gf + 88);
    const auto *gf_89 = buffer.data(gf + 89);
    const auto *gf_90 = buffer.data(gf + 90);
    const auto *gf_91 = buffer.data(gf + 91);
    const auto *gf_92 = buffer.data(gf + 92);
    const auto *gf_93 = buffer.data(gf + 93);
    const auto *gf_94 = buffer.data(gf + 94);
    const auto *gf_95 = buffer.data(gf + 95);
    const auto *gf_96 = buffer.data(gf + 96);
    const auto *gf_97 = buffer.data(gf + 97);
    const auto *gf_98 = buffer.data(gf + 98);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, ff_0, gd0_0, gd1_0, \
                         gf_0, gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = pb_y[k] * gf_2[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, ff_3, ff_6, gd0_1, gd1_1, \
                         gf_3, gf_4, gf_5, gf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * ff_3[k]
                 + pb_x[k] * gf_5[k];

        t_7[k] = pb_z[k] * gf_3[k];

        t_8[k] = pb_y[k] * gf_4[k];

        t_9[k] = f_0 * ff_6[k]
                 + pb_x[k] * gf_7[k];

        t_10[k] = f_1 * gd0_1[k]
                  - f_2 * gd1_1[k]
                  + pb_y[k] * gf_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, fg_0, gd0_2, gd1_2, \
                         gf_5, gf_6, gf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * gf_5[k];

        t_12[k] = f_3 * gd0_2[k]
                  - f_4 * gd1_2[k]
                  + pb_y[k] * gf_6[k];

        t_13[k] = pb_y[k] * gf_7[k];

        t_14[k] = f_1 * gd0_2[k]
                  - f_2 * gd1_2[k]
                  + pb_z[k] * gf_7[k];

        t_15[k] = pa_y[k] * fg_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pb_y, pb_z, ff_0, ff_1, fg_1, \
                         fg_2, gf_8, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * ff_0[k]
                  + pb_y[k] * gf_8[k];

        t_17[k] = pb_z[k] * gf_8[k];

        t_18[k] = f_6 * ff_1[k]
                  + pa_y[k] * fg_1[k];

        t_19[k] = pb_z[k] * gf_9[k];

        t_20[k] = pa_y[k] * fg_2[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pb_x, pb_z, ff_3, ff_8, ff_9, \
                         fg_4, fg_5, gf_10, gf_11, gf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * ff_8[k]
                  + pb_x[k] * gf_11[k];

        t_22[k] = pb_z[k] * gf_10[k];

        t_23[k] = f_7 * ff_9[k]
                  + pb_x[k] * gf_12[k];

        t_24[k] = pa_y[k] * fg_4[k];

        t_25[k] = f_0 * ff_3[k]
                  + pa_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, pb_z, ff_5, ff_6, \
                         fg_0, fg_6, fg_7, gf_11, gf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * gf_11[k];

        t_27[k] = f_6 * ff_5[k]
                  + pa_y[k] * fg_6[k];

        t_28[k] = f_5 * ff_6[k]
                  + pb_y[k] * gf_13[k];

        t_29[k] = pa_y[k] * fg_7[k];

        t_30[k] = pa_z[k] * fg_0[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, ff_0, ff_2, \
                         fg_1, fg_2, fg_3, gf_14, gf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * gf_14[k];

        t_32[k] = f_5 * ff_0[k]
                  + pb_z[k] * gf_14[k];

        t_33[k] = pa_z[k] * fg_1[k];

        t_34[k] = pb_y[k] * gf_15[k];

        t_35[k] = f_6 * ff_2[k]
                  + pa_z[k] * fg_2[k];

        t_36[k] = pa_z[k] * fg_3[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, pb_y, ff_14, ff_16, fg_5, gf_16, \
                         gf_18, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_7 * ff_14[k]
                  + pb_x[k] * gf_18[k];

        t_38[k] = pb_y[k] * gf_16[k];

        t_39[k] = f_7 * ff_16[k]
                  + pb_x[k] * gf_19[k];

        t_40[k] = pa_z[k] * fg_5[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_z, pb_y, pb_z, ff_3, ff_4, ff_6, fg_6, \
                         fg_7, gf_17, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * ff_3[k]
                  + pb_z[k] * gf_17[k];

        t_42[k] = f_6 * ff_4[k]
                  + pa_z[k] * fg_6[k];

        t_43[k] = pb_y[k] * gf_19[k];

        t_44[k] = f_0 * ff_6[k]
                  + pa_z[k] * fg_7[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_y, pb_z, dg0_0, dg1_0, ff_7, fg_8, \
                         gf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_8 * dg0_0[k]
                  - f_9 * dg1_0[k]
                  + pa_y[k] * fg_8[k];

        t_46[k] = f_6 * ff_7[k]
                  + pb_y[k] * gf_20[k];

        t_47[k] = pb_z[k] * gf_20[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_x, pb_z, ff_18, ff_19, gd0_3, gd0_4, \
                         gd1_3, gd1_4, gf_21, gf_22, gf_23, gf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_6 * ff_18[k]
                  + f_3 * gd0_4[k]
                  - f_4 * gd1_4[k]
                  + pb_x[k] * gf_23[k];

        t_49[k] = pb_z[k] * gf_21[k];

        t_50[k] = f_3 * gd0_3[k]
                  - f_4 * gd1_3[k]
                  + pb_z[k] * gf_22[k];

        t_51[k] = f_6 * ff_19[k]
                  + pb_x[k] * gf_24[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_z, dg0_1, dg1_1, ff_20, ff_21, \
                         fg_23, gf_23, gf_26, gf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_z[k] * gf_23[k];

        t_53[k] = f_6 * ff_20[k]
                  + pb_x[k] * gf_26[k];

        t_54[k] = f_6 * ff_21[k]
                  + pb_x[k] * gf_27[k];

        t_55[k] = f_8 * dg0_1[k]
                  - f_9 * dg1_1[k]
                  + pa_x[k] * fg_23[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pb_z, ff_10, gd0_4, gd0_5, gd1_4, \
                         gd1_5, gf_24, gf_25, gf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * gf_24[k];

        t_57[k] = f_3 * gd0_4[k]
                  - f_4 * gd1_4[k]
                  + pb_z[k] * gf_25[k];

        t_58[k] = f_6 * ff_10[k]
                  + pb_y[k] * gf_27[k];

        t_59[k] = f_1 * gd0_5[k]
                  - f_2 * gd1_5[k]
                  + pb_z[k] * gf_27[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, ff_12, fg_9, \
                         fg_10, fg_13, fg_14, fg_15, gf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * fg_13[k];

        t_61[k] = pa_z[k] * fg_9[k];

        t_62[k] = pa_y[k] * fg_14[k];

        t_63[k] = pa_z[k] * fg_10[k];

        t_64[k] = f_5 * ff_12[k]
                  + pb_y[k] * gf_28[k];

        t_65[k] = pa_y[k] * fg_15[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, ff_23, ff_24, fg_11, \
                         fg_12, fg_16, gf_30, gf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * fg_11[k];

        t_67[k] = f_6 * ff_23[k]
                  + pb_x[k] * gf_30[k];

        t_68[k] = f_6 * ff_24[k]
                  + pb_x[k] * gf_31[k];

        t_69[k] = pa_y[k] * fg_16[k];

        t_70[k] = pa_z[k] * fg_12[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, pb_z, ff_8, ff_15, ff_16, fg_17, \
                         fg_18, gf_29, gf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * ff_8[k]
                  + pb_z[k] * gf_29[k];

        t_72[k] = f_6 * ff_15[k]
                  + pa_y[k] * fg_17[k];

        t_73[k] = f_5 * ff_16[k]
                  + pb_y[k] * gf_32[k];

        t_74[k] = pa_y[k] * fg_18[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, dg0_0, dg1_0, ff_11, fg_13, \
                         gd0_6, gd1_6, gf_33, gf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_8 * dg0_0[k]
                  - f_9 * dg1_0[k]
                  + pa_z[k] * fg_13[k];

        t_76[k] = pb_y[k] * gf_33[k];

        t_77[k] = f_6 * ff_11[k]
                  + pb_z[k] * gf_33[k];

        t_78[k] = f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_y[k] * gf_34[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pb_x, pb_y, ff_27, ff_28, ff_29, gd0_8, \
                         gd1_8, gf_35, gf_36, gf_37, gf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_y[k] * gf_35[k];

        t_80[k] = f_6 * ff_27[k]
                  + f_3 * gd0_8[k]
                  - f_4 * gd1_8[k]
                  + pb_x[k] * gf_36[k];

        t_81[k] = f_6 * ff_28[k]
                  + pb_x[k] * gf_37[k];

        t_82[k] = f_6 * ff_29[k]
                  + pb_x[k] * gf_38[k];

        t_83[k] = pb_y[k] * gf_36[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_x, pb_y, pb_z, ff_13, ff_30, gd0_7, gd0_8, \
                         gd1_7, gd1_8, gf_37, gf_39, gf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_6 * ff_30[k]
                  + pb_x[k] * gf_40[k];

        t_85[k] = f_1 * gd0_7[k]
                  - f_2 * gd1_7[k]
                  + pb_y[k] * gf_37[k];

        t_86[k] = f_6 * ff_13[k]
                  + pb_z[k] * gf_37[k];

        t_87[k] = f_3 * gd0_8[k]
                  - f_4 * gd1_8[k]
                  + pb_y[k] * gf_39[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_x, pb_y, pb_z, dg0_2, dg1_2, ff_17, \
                         ff_31, fg_28, fg_29, gf_40, gf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pb_y[k] * gf_40[k];

        t_89[k] = f_8 * dg0_2[k]
                  - f_9 * dg1_2[k]
                  + pa_x[k] * fg_28[k];

        t_90[k] = f_0 * ff_31[k]
                  + pa_x[k] * fg_29[k];

        t_91[k] = f_7 * ff_17[k]
                  + pb_y[k] * gf_41[k];

        t_92[k] = pb_z[k] * gf_41[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pa_x, pb_x, pb_z, ff_33, ff_34, ff_35, \
                         fg_31, fg_32, gf_42, gf_43, gf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_6 * ff_33[k]
                  + pa_x[k] * fg_31[k];

        t_94[k] = pb_z[k] * gf_42[k];

        t_95[k] = f_6 * ff_34[k]
                  + pa_x[k] * fg_32[k];

        t_96[k] = f_5 * ff_35[k]
                  + pb_x[k] * gf_44[k];

        t_97[k] = pb_z[k] * gf_43[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, pa_x, pb_x, pb_z, ff_37, ff_38, \
                         fg_33, fg_34, gf_44, gf_45, gf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_5 * ff_37[k]
                  + pb_x[k] * gf_45[k];

        t_99[k] = f_5 * ff_38[k]
                  + pb_x[k] * gf_46[k];

        t_100[k] = pa_x[k] * fg_33[k];

        t_101[k] = pb_z[k] * gf_44[k];

        t_102[k] = pa_x[k] * fg_34[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, t_108, pa_x, pa_z, pb_z, ff_17, \
                         fg_19, fg_20, fg_21, fg_35, fg_36, gf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = pa_x[k] * fg_35[k];

        t_104[k] = pa_x[k] * fg_36[k];

        t_105[k] = pa_z[k] * fg_19[k];

        t_106[k] = pa_z[k] * fg_20[k];

        t_107[k] = f_5 * ff_17[k]
                   + pb_z[k] * gf_47[k];

        t_108[k] = pa_z[k] * fg_21[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_x, pa_z, pb_x, pb_y, ff_22, ff_41, \
                         ff_43, fg_22, fg_37, gf_48, gf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_6 * ff_22[k]
                   + pb_y[k] * gf_48[k];

        t_110[k] = f_6 * ff_41[k]
                   + pa_x[k] * fg_37[k];

        t_111[k] = pa_z[k] * fg_22[k];

        t_112[k] = f_5 * ff_43[k]
                   + pb_x[k] * gf_49[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, t_118, pa_x, pb_x, ff_44, ff_45, \
                         fg_38, fg_39, fg_40, fg_41, gf_50, gf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_5 * ff_44[k]
                   + pb_x[k] * gf_50[k];

        t_114[k] = f_5 * ff_45[k]
                   + pb_x[k] * gf_51[k];

        t_115[k] = pa_x[k] * fg_38[k];

        t_116[k] = pa_x[k] * fg_39[k];

        t_117[k] = pa_x[k] * fg_40[k];

        t_118[k] = pa_x[k] * fg_41[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, pa_x, pa_y, pb_y, ff_25, ff_48, \
                         fg_24, fg_25, fg_42, fg_43, gf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pa_x[k] * fg_42[k];

        t_120[k] = pa_y[k] * fg_24[k];

        t_121[k] = f_5 * ff_25[k]
                   + pb_y[k] * gf_52[k];

        t_122[k] = pa_y[k] * fg_25[k];

        t_123[k] = f_6 * ff_48[k]
                   + pa_x[k] * fg_43[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pb_x, pb_y, ff_26, ff_49, ff_50, \
                         fg_26, gf_53, gf_54, gf_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_5 * ff_26[k]
                   + pb_y[k] * gf_53[k];

        t_125[k] = pa_y[k] * fg_26[k];

        t_126[k] = f_5 * ff_49[k]
                   + pb_x[k] * gf_54[k];

        t_127[k] = f_5 * ff_50[k]
                   + pb_x[k] * gf_55[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, pa_x, pa_y, pb_x, ff_51, \
                         fg_27, fg_44, fg_45, fg_46, fg_47, gf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_5 * ff_51[k]
                   + pb_x[k] * gf_56[k];

        t_129[k] = pa_y[k] * fg_27[k];

        t_130[k] = pa_x[k] * fg_44[k];

        t_131[k] = pa_x[k] * fg_45[k];

        t_132[k] = pa_x[k] * fg_46[k];

        t_133[k] = pa_x[k] * fg_47[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, pa_x, pb_y, pb_z, ff_25, ff_53, \
                         ff_56, fg_48, fg_49, fg_51, gf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = pa_x[k] * fg_48[k];

        t_135[k] = f_0 * ff_53[k]
                   + pa_x[k] * fg_49[k];

        t_136[k] = pb_y[k] * gf_57[k];

        t_137[k] = f_7 * ff_25[k]
                   + pb_z[k] * gf_57[k];

        t_138[k] = f_6 * ff_56[k]
                   + pa_x[k] * fg_51[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pa_x, pb_x, pb_y, ff_57, ff_58, \
                         ff_59, fg_52, gf_58, gf_59, gf_60, gf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * gf_58[k];

        t_140[k] = f_6 * ff_57[k]
                   + pa_x[k] * fg_52[k];

        t_141[k] = f_5 * ff_58[k]
                   + pb_x[k] * gf_60[k];

        t_142[k] = f_5 * ff_59[k]
                   + pb_x[k] * gf_61[k];

        t_143[k] = pb_y[k] * gf_59[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, t_149, pa_x, pb_x, pb_y, ff_61, \
                         fg_53, fg_54, fg_55, fg_56, gf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_5 * ff_61[k]
                   + pb_x[k] * gf_62[k];

        t_145[k] = pa_x[k] * fg_53[k];

        t_146[k] = pa_x[k] * fg_54[k];

        t_147[k] = pa_x[k] * fg_55[k];

        t_148[k] = pb_y[k] * gf_62[k];

        t_149[k] = pa_x[k] * fg_56[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, pb_x, pb_y, pb_z, ff_31, gd0_9, \
                         gd0_10, gd1_9, gd1_10, gf_63, gf_64, gf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_1 * gd0_9[k]
                   - f_2 * gd1_9[k]
                   + pb_x[k] * gf_63[k];

        t_151[k] = f_0 * ff_31[k]
                   + pb_y[k] * gf_63[k];

        t_152[k] = pb_z[k] * gf_63[k];

        t_153[k] = f_3 * gd0_10[k]
                   - f_4 * gd1_10[k]
                   + pb_x[k] * gf_65[k];

        t_154[k] = pb_z[k] * gf_64[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pb_x, gd0_11, gd1_11, gf_66, \
                         gf_67, gf_68, gf_69, gf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_3 * gd0_11[k]
                   - f_4 * gd1_11[k]
                   + pb_x[k] * gf_66[k];

        t_156[k] = pb_x[k] * gf_67[k];

        t_157[k] = pb_x[k] * gf_68[k];

        t_158[k] = pb_x[k] * gf_69[k];

        t_159[k] = pb_x[k] * gf_70[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, pb_y, pb_z, ff_35, ff_38, gd0_10, \
                         gd0_11, gd1_10, gd1_11, gf_67, gf_68, gf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_0 * ff_35[k]
                   + f_1 * gd0_10[k]
                   - f_2 * gd1_10[k]
                   + pb_y[k] * gf_67[k];

        t_161[k] = pb_z[k] * gf_67[k];

        t_162[k] = f_3 * gd0_10[k]
                   - f_4 * gd1_10[k]
                   + pb_z[k] * gf_68[k];

        t_163[k] = f_0 * ff_38[k]
                   + pb_y[k] * gf_70[k];

        t_164[k] = f_1 * gd0_11[k]
                   - f_2 * gd1_11[k]
                   + pb_z[k] * gf_70[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pa_z, pb_y, pb_z, ff_31, ff_40, \
                         fg_29, fg_30, fg_31, gf_71, gf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = pa_z[k] * fg_29[k];

        t_166[k] = pa_z[k] * fg_30[k];

        t_167[k] = f_5 * ff_31[k]
                   + pb_z[k] * gf_71[k];

        t_168[k] = pa_z[k] * fg_31[k];

        t_169[k] = f_7 * ff_40[k]
                   + pb_y[k] * gf_72[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, t_175, pa_z, pb_x, ff_32, fg_32, \
                         fg_33, gf_73, gf_74, gf_75, gf_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_6 * ff_32[k]
                   + pa_z[k] * fg_32[k];

        t_171[k] = pb_x[k] * gf_73[k];

        t_172[k] = pb_x[k] * gf_74[k];

        t_173[k] = pb_x[k] * gf_75[k];

        t_174[k] = pb_x[k] * gf_76[k];

        t_175[k] = pa_z[k] * fg_33[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_z, pb_y, pb_z, ff_35, ff_36, ff_38, \
                         ff_45, fg_34, fg_36, gf_73, gf_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_5 * ff_35[k]
                   + pb_z[k] * gf_73[k];

        t_177[k] = f_6 * ff_36[k]
                   + pa_z[k] * fg_34[k];

        t_178[k] = f_7 * ff_45[k]
                   + pb_y[k] * gf_76[k];

        t_179[k] = f_0 * ff_38[k]
                   + pa_z[k] * fg_36[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pb_x, pb_y, pb_z, ff_39, ff_46, gd0_12, \
                         gd0_13, gd1_12, gd1_13, gf_77, gf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_1 * gd0_12[k]
                   - f_2 * gd1_12[k]
                   + pb_x[k] * gf_77[k];

        t_181[k] = f_6 * ff_46[k]
                   + pb_y[k] * gf_77[k];

        t_182[k] = f_6 * ff_39[k]
                   + pb_z[k] * gf_77[k];

        t_183[k] = f_3 * gd0_13[k]
                   - f_4 * gd1_13[k]
                   + pb_x[k] * gf_79[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pb_x, pb_y, ff_47, gd0_14, gd1_14, \
                         gf_78, gf_80, gf_81, gf_82, gf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_6 * ff_47[k]
                   + pb_y[k] * gf_78[k];

        t_185[k] = f_3 * gd0_14[k]
                   - f_4 * gd1_14[k]
                   + pb_x[k] * gf_80[k];

        t_186[k] = pb_x[k] * gf_81[k];

        t_187[k] = pb_x[k] * gf_82[k];

        t_188[k] = pb_x[k] * gf_83[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pa_z, pb_x, pb_z, dg0_1, dg1_1, ff_42, fg_38, \
                         gf_81, gf_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pb_x[k] * gf_84[k];

        t_190[k] = f_8 * dg0_1[k]
                   - f_9 * dg1_1[k]
                   + pa_z[k] * fg_38[k];

        t_191[k] = f_6 * ff_42[k]
                   + pb_z[k] * gf_81[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_y, pb_y, dg0_2, dg1_2, ff_51, ff_52, \
                         fg_48, fg_49, gd0_14, gd1_14, gf_83, gf_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_6 * ff_51[k]
                   + f_3 * gd0_14[k]
                   - f_4 * gd1_14[k]
                   + pb_y[k] * gf_83[k];

        t_193[k] = f_6 * ff_52[k]
                   + pb_y[k] * gf_84[k];

        t_194[k] = f_8 * dg0_2[k]
                   - f_9 * dg1_2[k]
                   + pa_y[k] * fg_48[k];

        t_195[k] = pa_y[k] * fg_49[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, pa_y, pb_y, ff_53, ff_54, ff_55, \
                         fg_50, fg_51, fg_52, gf_85, gf_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_5 * ff_53[k]
                   + pb_y[k] * gf_85[k];

        t_197[k] = pa_y[k] * fg_50[k];

        t_198[k] = f_6 * ff_54[k]
                   + pa_y[k] * fg_51[k];

        t_199[k] = f_5 * ff_55[k]
                   + pb_y[k] * gf_86[k];

        t_200[k] = pa_y[k] * fg_52[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, t_206, pa_y, pb_x, pb_z, ff_49, \
                         ff_58, fg_53, gf_87, gf_88, gf_89, gf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pb_x[k] * gf_87[k];

        t_202[k] = pb_x[k] * gf_88[k];

        t_203[k] = pb_x[k] * gf_89[k];

        t_204[k] = pb_x[k] * gf_90[k];

        t_205[k] = f_0 * ff_58[k]
                   + pa_y[k] * fg_53[k];

        t_206[k] = f_7 * ff_49[k]
                   + pb_z[k] * gf_87[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, pa_y, pb_x, pb_y, ff_60, ff_61, \
                         fg_55, fg_56, gd0_15, gd1_15, gf_90, gf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_6 * ff_60[k]
                   + pa_y[k] * fg_55[k];

        t_208[k] = f_5 * ff_61[k]
                   + pb_y[k] * gf_90[k];

        t_209[k] = pa_y[k] * fg_56[k];

        t_210[k] = f_1 * gd0_15[k]
                   - f_2 * gd1_15[k]
                   + pb_x[k] * gf_91[k];

        t_211[k] = pb_y[k] * gf_91[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pb_x, pb_y, pb_z, ff_53, gd0_16, gd0_17, \
                         gd1_16, gd1_17, gf_91, gf_92, gf_93, gf_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_0 * ff_53[k]
                   + pb_z[k] * gf_91[k];

        t_213[k] = f_3 * gd0_16[k]
                   - f_4 * gd1_16[k]
                   + pb_x[k] * gf_93[k];

        t_214[k] = pb_y[k] * gf_92[k];

        t_215[k] = f_3 * gd0_17[k]
                   - f_4 * gd1_17[k]
                   + pb_x[k] * gf_94[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, t_221, pb_x, pb_y, pb_z, ff_58, \
                         gd0_16, gd1_16, gf_95, gf_96, gf_97, gf_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = pb_x[k] * gf_95[k];

        t_217[k] = pb_x[k] * gf_96[k];

        t_218[k] = pb_x[k] * gf_97[k];

        t_219[k] = pb_x[k] * gf_98[k];

        t_220[k] = f_1 * gd0_16[k]
                   - f_2 * gd1_16[k]
                   + pb_y[k] * gf_95[k];

        t_221[k] = f_0 * ff_58[k]
                   + pb_z[k] * gf_95[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pb_y, pb_z, ff_61, gd0_17, gd1_17, gf_97, \
                         gf_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_3 * gd0_17[k]
                   - f_4 * gd1_17[k]
                   + pb_y[k] * gf_97[k];

        t_223[k] = pb_y[k] * gf_98[k];

        t_224[k] = f_0 * ff_61[k]
                   + f_1 * gd0_17[k]
                   - f_2 * gd1_17[k]
                   + pb_z[k] * gf_98[k];
    }
}

auto
compute_prim_gg_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg0, const size_t dg1,
                                     const size_t ff, const size_t fg, const size_t gd0,
                                     const size_t gd1, const size_t gf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_21 = buffer.data(dg1 + 21);
    const auto *dg1_41 = buffer.data(dg1 + 41);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_41 = buffer.data(ff + 41);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
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
    const auto *fg_29 = buffer.data(fg + 29);
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
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_56 = buffer.data(fg + 56);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_17 = buffer.data(gd0 + 17);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_17 = buffer.data(gd1 + 17);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ff_0, gd0_0, gd1_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, ff_3, ff_6, gd0_1, gd0_2, gd1_1, \
                         gd1_2, gf_3, gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * ff_3[k]
                 + pb_x[k] * gf_3[k];

        t_6[k] = f_0 * ff_6[k]
                 + pb_x[k] * gf_5[k];

        t_7[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];

        t_8[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_4[k];

        t_9[k] = pb_y[k] * gf_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, ff_0, ff_1, fg_0, fg_3, \
                         gd0_2, gd1_2, gf_5, gf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * gd0_2[k]
                  - f_2 * gd1_2[k]
                  + pb_z[k] * gf_5[k];

        t_11[k] = pa_y[k] * fg_0[k];

        t_12[k] = f_5 * ff_0[k]
                  + pb_y[k] * gf_6[k];

        t_13[k] = f_6 * ff_1[k]
                  + pa_y[k] * fg_3[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pb_x, ff_3, ff_5, ff_8, fg_4, fg_5, \
                         fg_6, gf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_y[k] * fg_4[k];

        t_15[k] = f_7 * ff_8[k]
                  + pb_x[k] * gf_7[k];

        t_16[k] = f_0 * ff_3[k]
                  + pa_y[k] * fg_5[k];

        t_17[k] = f_6 * ff_5[k]
                  + pa_y[k] * fg_6[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_y, pa_z, pb_y, pb_z, ff_0, ff_6, \
                         fg_0, fg_3, fg_8, gf_8, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * ff_6[k]
                  + pb_y[k] * gf_8[k];

        t_19[k] = pa_y[k] * fg_8[k];

        t_20[k] = pa_z[k] * fg_0[k];

        t_21[k] = f_5 * ff_0[k]
                  + pb_z[k] * gf_9[k];

        t_22[k] = pa_z[k] * fg_3[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pb_x, pb_z, ff_2, ff_3, ff_13, fg_4, \
                         fg_5, gf_10, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * ff_2[k]
                  + pa_z[k] * fg_4[k];

        t_24[k] = f_7 * ff_13[k]
                  + pb_x[k] * gf_11[k];

        t_25[k] = pa_z[k] * fg_5[k];

        t_26[k] = f_5 * ff_3[k]
                  + pb_z[k] * gf_10[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, dg0_0, dg1_0, ff_4, ff_6, \
                         ff_7, fg_6, fg_8, fg_9, gf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_6 * ff_4[k]
                  + pa_z[k] * fg_6[k];

        t_28[k] = f_0 * ff_6[k]
                  + pa_z[k] * fg_8[k];

        t_29[k] = f_8 * dg0_0[k]
                  - f_9 * dg1_0[k]
                  + pa_y[k] * fg_9[k];

        t_30[k] = f_6 * ff_7[k]
                  + pb_y[k] * gf_12[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pb_x, pb_z, ff_15, ff_16, gd0_3, gd0_4, \
                         gd1_3, gd1_4, gf_12, gf_13, gf_14, gf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_z[k] * gf_12[k];

        t_32[k] = f_6 * ff_15[k]
                  + f_3 * gd0_4[k]
                  - f_4 * gd1_4[k]
                  + pb_x[k] * gf_14[k];

        t_33[k] = f_3 * gd0_3[k]
                  - f_4 * gd1_3[k]
                  + pb_z[k] * gf_13[k];

        t_34[k] = f_6 * ff_16[k]
                  + pb_x[k] * gf_15[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_y, pb_z, dg0_1, dg1_21, ff_9, fg_19, \
                         gd0_4, gd1_4, gf_15, gf_16, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * dg0_1[k]
                  - f_9 * dg1_21[k]
                  + pa_x[k] * fg_19[k];

        t_36[k] = pb_z[k] * gf_15[k];

        t_37[k] = f_3 * gd0_4[k]
                  - f_4 * gd1_4[k]
                  + pb_z[k] * gf_16[k];

        t_38[k] = f_6 * ff_9[k]
                  + pb_y[k] * gf_17[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_y, pa_z, pb_z, fg_10, fg_11, fg_13, \
                         fg_14, gd0_5, gd1_5, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_1 * gd0_5[k]
                  - f_2 * gd1_5[k]
                  + pb_z[k] * gf_17[k];

        t_40[k] = pa_y[k] * fg_13[k];

        t_41[k] = pa_z[k] * fg_10[k];

        t_42[k] = pa_y[k] * fg_14[k];

        t_43[k] = pa_z[k] * fg_11[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_y, pb_y, pb_z, ff_8, ff_12, ff_13, fg_15, \
                         fg_16, gf_18, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_5 * ff_8[k]
                  + pb_z[k] * gf_18[k];

        t_45[k] = f_6 * ff_12[k]
                  + pa_y[k] * fg_15[k];

        t_46[k] = f_5 * ff_13[k]
                  + pb_y[k] * gf_19[k];

        t_47[k] = pa_y[k] * fg_16[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_y, pb_z, dg0_0, dg1_0, ff_10, fg_12, \
                         gd0_6, gd1_6, gf_20, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_8 * dg0_0[k]
                  - f_9 * dg1_0[k]
                  + pa_z[k] * fg_12[k];

        t_49[k] = pb_y[k] * gf_20[k];

        t_50[k] = f_6 * ff_10[k]
                  + pb_z[k] * gf_20[k];

        t_51[k] = f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_y[k] * gf_21[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_x, pb_y, ff_18, ff_19, gd0_7, gd0_8, gd1_7, \
                         gd1_8, gf_22, gf_23, gf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_6 * ff_18[k]
                  + f_3 * gd0_8[k]
                  - f_4 * gd1_8[k]
                  + pb_x[k] * gf_22[k];

        t_53[k] = f_6 * ff_19[k]
                  + pb_x[k] * gf_25[k];

        t_54[k] = f_1 * gd0_7[k]
                  - f_2 * gd1_7[k]
                  + pb_y[k] * gf_23[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_x, pb_y, pb_z, dg0_2, dg1_41, ff_11, \
                         fg_23, gd0_8, gd1_8, gf_23, gf_24, gf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_6 * ff_11[k]
                  + pb_z[k] * gf_23[k];

        t_56[k] = f_3 * gd0_8[k]
                  - f_4 * gd1_8[k]
                  + pb_y[k] * gf_24[k];

        t_57[k] = pb_y[k] * gf_25[k];

        t_58[k] = f_8 * dg0_2[k]
                  - f_9 * dg1_41[k]
                  + pa_x[k] * fg_23[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_x, pb_y, ff_14, ff_20, ff_22, ff_23, \
                         fg_24, fg_25, fg_26, gf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_0 * ff_20[k]
                  + pa_x[k] * fg_24[k];

        t_60[k] = f_7 * ff_14[k]
                  + pb_y[k] * gf_26[k];

        t_61[k] = f_6 * ff_22[k]
                  + pa_x[k] * fg_25[k];

        t_62[k] = f_6 * ff_23[k]
                  + pa_x[k] * fg_26[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, t_68, pa_x, pa_z, pb_x, ff_24, fg_17, \
                         fg_29, fg_31, fg_32, fg_33, gf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_5 * ff_24[k]
                  + pb_x[k] * gf_27[k];

        t_64[k] = pa_x[k] * fg_29[k];

        t_65[k] = pa_x[k] * fg_31[k];

        t_66[k] = pa_x[k] * fg_32[k];

        t_67[k] = pa_x[k] * fg_33[k];

        t_68[k] = pa_z[k] * fg_17[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pa_x, pa_z, pb_z, ff_14, ff_28, fg_18, \
                         fg_34, fg_36, fg_37, gf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_5 * ff_14[k]
                  + pb_z[k] * gf_28[k];

        t_70[k] = pa_z[k] * fg_18[k];

        t_71[k] = f_6 * ff_28[k]
                  + pa_x[k] * fg_34[k];

        t_72[k] = pa_x[k] * fg_36[k];

        t_73[k] = pa_x[k] * fg_37[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, pa_x, pa_y, ff_31, fg_20, fg_21, \
                         fg_22, fg_38, fg_39, fg_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pa_x[k] * fg_38[k];

        t_75[k] = pa_x[k] * fg_39[k];

        t_76[k] = pa_y[k] * fg_20[k];

        t_77[k] = pa_y[k] * fg_21[k];

        t_78[k] = f_6 * ff_31[k]
                  + pa_x[k] * fg_40[k];

        t_79[k] = pa_y[k] * fg_22[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, t_85, pa_x, pb_z, ff_17, ff_35, fg_41, \
                         fg_42, fg_43, fg_44, fg_46, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pa_x[k] * fg_41[k];

        t_81[k] = pa_x[k] * fg_42[k];

        t_82[k] = pa_x[k] * fg_43[k];

        t_83[k] = pa_x[k] * fg_44[k];

        t_84[k] = f_0 * ff_35[k]
                  + pa_x[k] * fg_46[k];

        t_85[k] = f_7 * ff_17[k]
                  + pb_z[k] * gf_29[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, pa_x, pb_x, ff_37, ff_38, ff_41, fg_48, \
                         fg_49, fg_52, fg_53, gf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_6 * ff_37[k]
                  + pa_x[k] * fg_48[k];

        t_87[k] = f_6 * ff_38[k]
                  + pa_x[k] * fg_49[k];

        t_88[k] = f_5 * ff_41[k]
                  + pb_x[k] * gf_30[k];

        t_89[k] = pa_x[k] * fg_52[k];

        t_90[k] = pa_x[k] * fg_53[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_x, pb_x, pb_y, ff_20, fg_54, fg_56, gd0_9, \
                         gd1_9, gf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = pa_x[k] * fg_54[k];

        t_92[k] = pa_x[k] * fg_56[k];

        t_93[k] = f_1 * gd0_9[k]
                  - f_2 * gd1_9[k]
                  + pb_x[k] * gf_31[k];

        t_94[k] = f_0 * ff_20[k]
                  + pb_y[k] * gf_31[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pb_x, pb_y, ff_24, gd0_10, gd0_11, \
                         gd1_10, gd1_11, gf_32, gf_33, gf_34, gf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_x[k] * gf_32[k];

        t_96[k] = f_3 * gd0_11[k]
                  - f_4 * gd1_11[k]
                  + pb_x[k] * gf_33[k];

        t_97[k] = pb_x[k] * gf_34[k];

        t_98[k] = pb_x[k] * gf_36[k];

        t_99[k] = f_0 * ff_24[k]
                  + f_1 * gd0_10[k]
                  - f_2 * gd1_10[k]
                  + pb_y[k] * gf_34[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pb_y, pb_z, ff_26, gd0_10, gd0_11, \
                         gd1_10, gd1_11, gf_34, gf_35, gf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pb_z[k] * gf_34[k];

        t_101[k] = f_3 * gd0_10[k]
                   - f_4 * gd1_10[k]
                   + pb_z[k] * gf_35[k];

        t_102[k] = f_0 * ff_26[k]
                   + pb_y[k] * gf_36[k];

        t_103[k] = f_1 * gd0_11[k]
                   - f_2 * gd1_11[k]
                   + pb_z[k] * gf_36[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pa_z, pb_z, ff_20, ff_21, fg_24, \
                         fg_25, fg_26, fg_29, gf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_z[k] * fg_24[k];

        t_105[k] = f_5 * ff_20[k]
                   + pb_z[k] * gf_37[k];

        t_106[k] = pa_z[k] * fg_25[k];

        t_107[k] = f_6 * ff_21[k]
                   + pa_z[k] * fg_26[k];

        t_108[k] = pa_z[k] * fg_29[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_z, pb_y, pb_z, ff_24, ff_25, ff_26, \
                         ff_30, fg_31, fg_33, gf_38, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_5 * ff_24[k]
                   + pb_z[k] * gf_38[k];

        t_110[k] = f_6 * ff_25[k]
                   + pa_z[k] * fg_31[k];

        t_111[k] = f_7 * ff_30[k]
                   + pb_y[k] * gf_39[k];

        t_112[k] = f_0 * ff_26[k]
                   + pa_z[k] * fg_33[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pb_x, pb_z, ff_27, gd0_12, gd0_13, \
                         gd0_14, gd1_12, gd1_13, gd1_14, gf_40, gf_41, \
                         gf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_1 * gd0_12[k]
                   - f_2 * gd1_12[k]
                   + pb_x[k] * gf_40[k];

        t_114[k] = f_6 * ff_27[k]
                   + pb_z[k] * gf_40[k];

        t_115[k] = f_3 * gd0_13[k]
                   - f_4 * gd1_13[k]
                   + pb_x[k] * gf_41[k];

        t_116[k] = f_3 * gd0_14[k]
                   - f_4 * gd1_14[k]
                   + pb_x[k] * gf_42[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_z, pb_x, pb_z, dg0_1, dg1_21, ff_29, \
                         fg_35, gf_43, gf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pb_x[k] * gf_43[k];

        t_118[k] = pb_x[k] * gf_45[k];

        t_119[k] = f_8 * dg0_1[k]
                   - f_9 * dg1_21[k]
                   + pa_z[k] * fg_35[k];

        t_120[k] = f_6 * ff_29[k]
                   + pb_z[k] * gf_43[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pa_y, pb_y, dg0_2, dg1_41, ff_33, ff_34, \
                         fg_45, fg_46, gd0_14, gd1_14, gf_44, gf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_6 * ff_33[k]
                   + f_3 * gd0_14[k]
                   - f_4 * gd1_14[k]
                   + pb_y[k] * gf_44[k];

        t_122[k] = f_6 * ff_34[k]
                   + pb_y[k] * gf_45[k];

        t_123[k] = f_8 * dg0_2[k]
                   - f_9 * dg1_41[k]
                   + pa_y[k] * fg_45[k];

        t_124[k] = pa_y[k] * fg_46[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, pa_y, pb_z, ff_32, ff_36, ff_39, \
                         fg_47, fg_48, fg_49, fg_52, gf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pa_y[k] * fg_47[k];

        t_126[k] = f_6 * ff_36[k]
                   + pa_y[k] * fg_48[k];

        t_127[k] = pa_y[k] * fg_49[k];

        t_128[k] = f_0 * ff_39[k]
                   + pa_y[k] * fg_52[k];

        t_129[k] = f_7 * ff_32[k]
                   + pb_z[k] * gf_46[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pa_y, pb_x, pb_y, ff_40, ff_41, fg_54, \
                         fg_56, gd0_15, gd1_15, gf_47, gf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_6 * ff_40[k]
                   + pa_y[k] * fg_54[k];

        t_131[k] = f_5 * ff_41[k]
                   + pb_y[k] * gf_47[k];

        t_132[k] = pa_y[k] * fg_56[k];

        t_133[k] = f_1 * gd0_15[k]
                   - f_2 * gd1_15[k]
                   + pb_x[k] * gf_48[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pb_x, pb_z, ff_35, gd0_16, gd0_17, \
                         gd1_16, gd1_17, gf_48, gf_49, gf_50, gf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_0 * ff_35[k]
                   + pb_z[k] * gf_48[k];

        t_135[k] = f_3 * gd0_16[k]
                   - f_4 * gd1_16[k]
                   + pb_x[k] * gf_49[k];

        t_136[k] = f_3 * gd0_17[k]
                   - f_4 * gd1_17[k]
                   + pb_x[k] * gf_50[k];

        t_137[k] = pb_x[k] * gf_51[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, pb_x, pb_y, pb_z, ff_39, gd0_16, \
                         gd0_17, gd1_16, gd1_17, gf_51, gf_52, gf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = pb_x[k] * gf_53[k];

        t_139[k] = f_1 * gd0_16[k]
                   - f_2 * gd1_16[k]
                   + pb_y[k] * gf_51[k];

        t_140[k] = f_0 * ff_39[k]
                   + pb_z[k] * gf_51[k];

        t_141[k] = f_3 * gd0_17[k]
                   - f_4 * gd1_17[k]
                   + pb_y[k] * gf_52[k];

        t_142[k] = pb_y[k] * gf_53[k];
    }

#pragma omp simd aligned(t_143, pb_z, ff_41, gd0_17, gd1_17, gf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_0 * ff_41[k]
                   + f_1 * gd0_17[k]
                   - f_2 * gd1_17[k]
                   + pb_z[k] * gf_53[k];
    }
}

auto
compute_prim_gg_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg0, const size_t dg1,
                                     const size_t ff, const size_t fg, const size_t gd0,
                                     const size_t gd1, const size_t gf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / p;

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

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_1 = buffer.data(dg1 + 1);
    const auto *dg1_2 = buffer.data(dg1 + 2);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_23 = buffer.data(ff + 23);

    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_17 = buffer.data(gd0 + 17);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_17 = buffer.data(gd1 + 17);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ff_0, gd0_0, gd1_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, gd0_1, gd0_2, gd1_1, gd1_2, gf_3, \
                         gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];

        t_6[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_4[k];

        t_7[k] = pb_y[k] * gf_5[k];

        t_8[k] = f_1 * gd0_2[k]
                 - f_2 * gd1_2[k]
                 + pb_z[k] * gf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pb_z, dg0_0, dg1_0, ff_6, fg_9, gd0_4, \
                         gd1_4, gf_6, gf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * dg0_0[k]
                 - f_6 * dg1_0[k]
                 + pa_y[k] * fg_9[k];

        t_10[k] = pb_z[k] * gf_6[k];

        t_11[k] = f_7 * ff_6[k]
                  + f_3 * gd0_4[k]
                  - f_4 * gd1_4[k]
                  + pb_x[k] * gf_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, pb_z, dg0_1, dg1_1, ff_7, fg_11, \
                         gd0_3, gd1_3, gf_7, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * gd0_3[k]
                  - f_4 * gd1_3[k]
                  + pb_z[k] * gf_7[k];

        t_13[k] = f_7 * ff_7[k]
                  + pb_x[k] * gf_9[k];

        t_14[k] = f_5 * dg0_1[k]
                  - f_6 * dg1_1[k]
                  + pa_x[k] * fg_11[k];

        t_15[k] = pb_z[k] * gf_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_z, pb_z, dg0_0, dg1_0, fg_10, gd0_4, gd0_5, \
                         gd1_4, gd1_5, gf_10, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * gd0_4[k]
                  - f_4 * gd1_4[k]
                  + pb_z[k] * gf_10[k];

        t_17[k] = f_1 * gd0_5[k]
                  - f_2 * gd1_5[k]
                  + pb_z[k] * gf_11[k];

        t_18[k] = f_5 * dg0_0[k]
                  - f_6 * dg1_0[k]
                  + pa_z[k] * fg_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_x, pb_y, ff_8, ff_9, gd0_6, gd0_8, gd1_6, \
                         gd1_8, gf_12, gf_13, gf_14, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * gf_12[k];

        t_20[k] = f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_y[k] * gf_13[k];

        t_21[k] = f_7 * ff_8[k]
                  + f_3 * gd0_8[k]
                  - f_4 * gd1_8[k]
                  + pb_x[k] * gf_14[k];

        t_22[k] = f_7 * ff_9[k]
                  + pb_x[k] * gf_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_y, dg0_2, dg1_2, fg_12, gd0_7, \
                         gd0_8, gd1_7, gd1_8, gf_15, gf_16, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * gd0_7[k]
                  - f_2 * gd1_7[k]
                  + pb_y[k] * gf_15[k];

        t_24[k] = f_3 * gd0_8[k]
                  - f_4 * gd1_8[k]
                  + pb_y[k] * gf_16[k];

        t_25[k] = pb_y[k] * gf_17[k];

        t_26[k] = f_5 * dg0_2[k]
                  - f_6 * dg1_2[k]
                  + pa_x[k] * fg_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_x, gd0_9, gd0_10, gd0_11, gd1_9, gd1_10, \
                         gd1_11, gf_18, gf_19, gf_20, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * gd0_9[k]
                  - f_2 * gd1_9[k]
                  + pb_x[k] * gf_18[k];

        t_28[k] = f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_x[k] * gf_19[k];

        t_29[k] = f_3 * gd0_11[k]
                  - f_4 * gd1_11[k]
                  + pb_x[k] * gf_20[k];

        t_30[k] = pb_x[k] * gf_21[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pb_x, pb_y, pb_z, ff_13, gd0_10, \
                         gd0_11, gd1_10, gd1_11, gf_21, gf_22, gf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_x[k] * gf_23[k];

        t_32[k] = f_0 * ff_13[k]
                  + f_1 * gd0_10[k]
                  - f_2 * gd1_10[k]
                  + pb_y[k] * gf_21[k];

        t_33[k] = pb_z[k] * gf_21[k];

        t_34[k] = f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_z[k] * gf_22[k];

        t_35[k] = f_1 * gd0_11[k]
                  - f_2 * gd1_11[k]
                  + pb_z[k] * gf_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_x, gd0_12, gd0_13, gd0_14, gd1_12, gd1_13, \
                         gd1_14, gf_24, gf_25, gf_26, gf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * gd0_12[k]
                  - f_2 * gd1_12[k]
                  + pb_x[k] * gf_24[k];

        t_37[k] = f_3 * gd0_13[k]
                  - f_4 * gd1_13[k]
                  + pb_x[k] * gf_25[k];

        t_38[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_x[k] * gf_26[k];

        t_39[k] = pb_x[k] * gf_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pb_x, pb_y, dg0_1, dg1_1, ff_16, ff_17, \
                         fg_22, gd0_14, gd1_14, gf_28, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_x[k] * gf_29[k];

        t_41[k] = f_5 * dg0_1[k]
                  - f_6 * dg1_1[k]
                  + pa_z[k] * fg_22[k];

        t_42[k] = f_7 * ff_16[k]
                  + f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_y[k] * gf_28[k];

        t_43[k] = f_7 * ff_17[k]
                  + pb_y[k] * gf_29[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, pb_x, dg0_2, dg1_2, fg_23, gd0_15, gd0_16, \
                         gd1_15, gd1_16, gf_30, gf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_5 * dg0_2[k]
                  - f_6 * dg1_2[k]
                  + pa_y[k] * fg_23[k];

        t_45[k] = f_1 * gd0_15[k]
                  - f_2 * gd1_15[k]
                  + pb_x[k] * gf_30[k];

        t_46[k] = f_3 * gd0_16[k]
                  - f_4 * gd1_16[k]
                  + pb_x[k] * gf_31[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, pb_x, pb_y, gd0_16, gd0_17, \
                         gd1_16, gd1_17, gf_32, gf_33, gf_34, gf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * gd0_17[k]
                  - f_4 * gd1_17[k]
                  + pb_x[k] * gf_32[k];

        t_48[k] = pb_x[k] * gf_33[k];

        t_49[k] = pb_x[k] * gf_35[k];

        t_50[k] = f_1 * gd0_16[k]
                  - f_2 * gd1_16[k]
                  + pb_y[k] * gf_33[k];

        t_51[k] = f_3 * gd0_17[k]
                  - f_4 * gd1_17[k]
                  + pb_y[k] * gf_34[k];

        t_52[k] = pb_y[k] * gf_35[k];
    }

#pragma omp simd aligned(t_53, pb_z, ff_23, gd0_17, gd1_17, gf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * ff_23[k]
                  + f_1 * gd0_17[k]
                  - f_2 * gd1_17[k]
                  + pb_z[k] * gf_35[k];
    }
}

auto
compute_prim_gg_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg0, const size_t dg1,
                                     const size_t ff, const size_t fg, const size_t gd0,
                                     const size_t gd1, const size_t gf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / p;

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

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_17 = buffer.data(dg1 + 17);
    const auto *dg1_32 = buffer.data(dg1 + 32);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_26 = buffer.data(ff + 26);

    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_48 = buffer.data(fg + 48);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_17 = buffer.data(gd0 + 17);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_17 = buffer.data(gd1 + 17);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ff_0, gd0_0, gd1_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, gd0_1, gd0_2, gd1_1, gd1_2, gf_3, \
                         gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];

        t_6[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_4[k];

        t_7[k] = pb_y[k] * gf_5[k];

        t_8[k] = f_1 * gd0_2[k]
                 - f_2 * gd1_2[k]
                 + pb_z[k] * gf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pb_z, dg0_0, dg1_0, ff_8, fg_10, gd0_4, \
                         gd1_4, gf_6, gf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * dg0_0[k]
                 - f_6 * dg1_0[k]
                 + pa_y[k] * fg_10[k];

        t_10[k] = pb_z[k] * gf_6[k];

        t_11[k] = f_7 * ff_8[k]
                  + f_3 * gd0_4[k]
                  - f_4 * gd1_4[k]
                  + pb_x[k] * gf_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, pb_z, dg0_1, dg1_17, ff_9, fg_22, \
                         gd0_3, gd1_3, gf_7, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * gd0_3[k]
                  - f_4 * gd1_3[k]
                  + pb_z[k] * gf_7[k];

        t_13[k] = f_7 * ff_9[k]
                  + pb_x[k] * gf_9[k];

        t_14[k] = f_5 * dg0_1[k]
                  - f_6 * dg1_17[k]
                  + pa_x[k] * fg_22[k];

        t_15[k] = pb_z[k] * gf_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_z, pb_z, dg0_0, dg1_0, fg_13, gd0_4, gd0_5, \
                         gd1_4, gd1_5, gf_10, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * gd0_4[k]
                  - f_4 * gd1_4[k]
                  + pb_z[k] * gf_10[k];

        t_17[k] = f_1 * gd0_5[k]
                  - f_2 * gd1_5[k]
                  + pb_z[k] * gf_11[k];

        t_18[k] = f_5 * dg0_0[k]
                  - f_6 * dg1_0[k]
                  + pa_z[k] * fg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_x, pb_y, ff_10, ff_11, gd0_6, gd0_8, \
                         gd1_6, gd1_8, gf_12, gf_13, gf_14, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * gf_12[k];

        t_20[k] = f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_y[k] * gf_13[k];

        t_21[k] = f_7 * ff_10[k]
                  + f_3 * gd0_8[k]
                  - f_4 * gd1_8[k]
                  + pb_x[k] * gf_14[k];

        t_22[k] = f_7 * ff_11[k]
                  + pb_x[k] * gf_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_y, dg0_2, dg1_32, fg_26, gd0_7, \
                         gd0_8, gd1_7, gd1_8, gf_15, gf_16, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * gd0_7[k]
                  - f_2 * gd1_7[k]
                  + pb_y[k] * gf_15[k];

        t_24[k] = f_3 * gd0_8[k]
                  - f_4 * gd1_8[k]
                  + pb_y[k] * gf_16[k];

        t_25[k] = pb_y[k] * gf_17[k];

        t_26[k] = f_5 * dg0_2[k]
                  - f_6 * dg1_32[k]
                  + pa_x[k] * fg_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_x, gd0_9, gd0_10, gd0_11, gd1_9, gd1_10, \
                         gd1_11, gf_18, gf_19, gf_20, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * gd0_9[k]
                  - f_2 * gd1_9[k]
                  + pb_x[k] * gf_18[k];

        t_28[k] = f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_x[k] * gf_19[k];

        t_29[k] = f_3 * gd0_11[k]
                  - f_4 * gd1_11[k]
                  + pb_x[k] * gf_20[k];

        t_30[k] = pb_x[k] * gf_21[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pb_x, pb_y, pb_z, ff_15, gd0_10, \
                         gd0_11, gd1_10, gd1_11, gf_21, gf_22, gf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_x[k] * gf_23[k];

        t_32[k] = f_0 * ff_15[k]
                  + f_1 * gd0_10[k]
                  - f_2 * gd1_10[k]
                  + pb_y[k] * gf_21[k];

        t_33[k] = pb_z[k] * gf_21[k];

        t_34[k] = f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_z[k] * gf_22[k];

        t_35[k] = f_1 * gd0_11[k]
                  - f_2 * gd1_11[k]
                  + pb_z[k] * gf_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_x, gd0_12, gd0_13, gd0_14, gd1_12, gd1_13, \
                         gd1_14, gf_24, gf_25, gf_26, gf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * gd0_12[k]
                  - f_2 * gd1_12[k]
                  + pb_x[k] * gf_24[k];

        t_37[k] = f_3 * gd0_13[k]
                  - f_4 * gd1_13[k]
                  + pb_x[k] * gf_25[k];

        t_38[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_x[k] * gf_26[k];

        t_39[k] = pb_x[k] * gf_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pb_x, pb_y, dg0_1, dg1_17, ff_19, \
                         ff_20, fg_40, gd0_14, gd1_14, gf_28, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_x[k] * gf_29[k];

        t_41[k] = f_5 * dg0_1[k]
                  - f_6 * dg1_17[k]
                  + pa_z[k] * fg_40[k];

        t_42[k] = f_7 * ff_19[k]
                  + f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_y[k] * gf_28[k];

        t_43[k] = f_7 * ff_20[k]
                  + pb_y[k] * gf_29[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, pb_x, dg0_2, dg1_32, fg_48, gd0_15, gd0_16, \
                         gd1_15, gd1_16, gf_30, gf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_5 * dg0_2[k]
                  - f_6 * dg1_32[k]
                  + pa_y[k] * fg_48[k];

        t_45[k] = f_1 * gd0_15[k]
                  - f_2 * gd1_15[k]
                  + pb_x[k] * gf_30[k];

        t_46[k] = f_3 * gd0_16[k]
                  - f_4 * gd1_16[k]
                  + pb_x[k] * gf_31[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, pb_x, pb_y, gd0_16, gd0_17, \
                         gd1_16, gd1_17, gf_32, gf_33, gf_34, gf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * gd0_17[k]
                  - f_4 * gd1_17[k]
                  + pb_x[k] * gf_32[k];

        t_48[k] = pb_x[k] * gf_33[k];

        t_49[k] = pb_x[k] * gf_35[k];

        t_50[k] = f_1 * gd0_16[k]
                  - f_2 * gd1_16[k]
                  + pb_y[k] * gf_33[k];

        t_51[k] = f_3 * gd0_17[k]
                  - f_4 * gd1_17[k]
                  + pb_y[k] * gf_34[k];

        t_52[k] = pb_y[k] * gf_35[k];
    }

#pragma omp simd aligned(t_53, pb_z, ff_26, gd0_17, gd1_17, gf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * ff_26[k]
                  + f_1 * gd0_17[k]
                  - f_2 * gd1_17[k]
                  + pb_z[k] * gf_35[k];
    }
}

auto
compute_prim_gg_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg0, const size_t dg1,
                                     const size_t ff, const size_t fg, const size_t gd0,
                                     const size_t gd1, const size_t gf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / p;

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

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_17 = buffer.data(dg0 + 17);
    const auto *dg0_32 = buffer.data(dg0 + 32);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_17 = buffer.data(dg1 + 17);
    const auto *dg1_32 = buffer.data(dg1 + 32);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_29 = buffer.data(ff + 29);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_38 = buffer.data(fg + 38);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_17 = buffer.data(gd0 + 17);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_17 = buffer.data(gd1 + 17);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ff_0, gd0_0, gd1_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, fg_0, gd0_1, gd0_2, gd1_1, \
                         gd1_2, gf_3, gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];

        t_6[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_4[k];

        t_7[k] = pb_y[k] * gf_5[k];

        t_8[k] = f_1 * gd0_2[k]
                 - f_2 * gd1_2[k]
                 + pb_z[k] * gf_5[k];

        t_9[k] = pa_y[k] * fg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pa_y, pa_z, dg0_0, dg1_0, ff_3, \
                         ff_5, fg_0, fg_5, fg_8, fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * ff_3[k]
                  + pa_y[k] * fg_5[k];

        t_11[k] = pa_y[k] * fg_8[k];

        t_12[k] = pa_z[k] * fg_0[k];

        t_13[k] = pa_z[k] * fg_5[k];

        t_14[k] = f_0 * ff_5[k]
                  + pa_z[k] * fg_8[k];

        t_15[k] = f_5 * dg0_0[k]
                  - f_6 * dg1_0[k]
                  + pa_y[k] * fg_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_x, pb_z, ff_9, ff_10, gd0_3, gd0_4, gd1_3, \
                         gd1_4, gf_6, gf_7, gf_8, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_z[k] * gf_6[k];

        t_17[k] = f_7 * ff_9[k]
                  + f_3 * gd0_4[k]
                  - f_4 * gd1_4[k]
                  + pb_x[k] * gf_8[k];

        t_18[k] = f_3 * gd0_3[k]
                  - f_4 * gd1_3[k]
                  + pb_z[k] * gf_7[k];

        t_19[k] = f_7 * ff_10[k]
                  + pb_x[k] * gf_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_z, dg0_17, dg1_17, fg_14, gd0_4, \
                         gd0_5, gd1_4, gd1_5, gf_9, gf_10, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * dg0_17[k]
                  - f_6 * dg1_17[k]
                  + pa_x[k] * fg_14[k];

        t_21[k] = pb_z[k] * gf_9[k];

        t_22[k] = f_3 * gd0_4[k]
                  - f_4 * gd1_4[k]
                  + pb_z[k] * gf_10[k];

        t_23[k] = f_1 * gd0_5[k]
                  - f_2 * gd1_5[k]
                  + pb_z[k] * gf_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pa_z, pb_y, dg0_0, dg1_0, fg_10, fg_11, \
                         fg_12, gf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_z[k] * fg_10[k];

        t_25[k] = pa_y[k] * fg_12[k];

        t_26[k] = f_5 * dg0_0[k]
                  - f_6 * dg1_0[k]
                  + pa_z[k] * fg_11[k];

        t_27[k] = pb_y[k] * gf_12[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pb_x, pb_y, ff_11, ff_12, gd0_6, gd0_8, gd1_6, \
                         gd1_8, gf_13, gf_14, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_y[k] * gf_13[k];

        t_29[k] = f_7 * ff_11[k]
                  + f_3 * gd0_8[k]
                  - f_4 * gd1_8[k]
                  + pb_x[k] * gf_14[k];

        t_30[k] = f_7 * ff_12[k]
                  + pb_x[k] * gf_17[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pb_y, dg0_32, dg1_32, fg_16, gd0_7, \
                         gd0_8, gd1_7, gd1_8, gf_15, gf_16, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_1 * gd0_7[k]
                  - f_2 * gd1_7[k]
                  + pb_y[k] * gf_15[k];

        t_32[k] = f_3 * gd0_8[k]
                  - f_4 * gd1_8[k]
                  + pb_y[k] * gf_16[k];

        t_33[k] = pb_y[k] * gf_17[k];

        t_34[k] = f_5 * dg0_32[k]
                  - f_6 * dg1_32[k]
                  + pa_x[k] * fg_16[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_x, pa_z, ff_13, ff_24, fg_13, fg_17, \
                         fg_22, fg_30, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * ff_13[k]
                  + pa_x[k] * fg_17[k];

        t_36[k] = pa_x[k] * fg_22[k];

        t_37[k] = pa_z[k] * fg_13[k];

        t_38[k] = f_0 * ff_24[k]
                  + pa_x[k] * fg_30[k];

        t_39[k] = pa_x[k] * fg_38[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_x, gd0_9, gd0_10, gd0_11, gd1_9, gd1_10, \
                         gd1_11, gf_18, gf_19, gf_20, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_1 * gd0_9[k]
                  - f_2 * gd1_9[k]
                  + pb_x[k] * gf_18[k];

        t_41[k] = f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_x[k] * gf_19[k];

        t_42[k] = f_3 * gd0_11[k]
                  - f_4 * gd1_11[k]
                  + pb_x[k] * gf_20[k];

        t_43[k] = pb_x[k] * gf_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pb_x, pb_y, pb_z, ff_16, gd0_10, \
                         gd0_11, gd1_10, gd1_11, gf_21, gf_22, gf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * gf_23[k];

        t_45[k] = f_0 * ff_16[k]
                  + f_1 * gd0_10[k]
                  - f_2 * gd1_10[k]
                  + pb_y[k] * gf_21[k];

        t_46[k] = pb_z[k] * gf_21[k];

        t_47[k] = f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_z[k] * gf_22[k];

        t_48[k] = f_1 * gd0_11[k]
                  - f_2 * gd1_11[k]
                  + pb_z[k] * gf_23[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_z, pb_x, ff_18, fg_17, fg_22, fg_25, \
                         gd0_12, gd1_12, gf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pa_z[k] * fg_17[k];

        t_50[k] = pa_z[k] * fg_22[k];

        t_51[k] = f_0 * ff_18[k]
                  + pa_z[k] * fg_25[k];

        t_52[k] = f_1 * gd0_12[k]
                  - f_2 * gd1_12[k]
                  + pb_x[k] * gf_24[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_x, gd0_13, gd0_14, gd1_13, gd1_14, gf_25, \
                         gf_26, gf_27, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_3 * gd0_13[k]
                  - f_4 * gd1_13[k]
                  + pb_x[k] * gf_25[k];

        t_54[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_x[k] * gf_26[k];

        t_55[k] = pb_x[k] * gf_27[k];

        t_56[k] = pb_x[k] * gf_29[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_z, pb_y, dg0_17, dg1_17, ff_22, ff_23, fg_26, \
                         gd0_14, gd1_14, gf_28, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_5 * dg0_17[k]
                  - f_6 * dg1_17[k]
                  + pa_z[k] * fg_26[k];

        t_58[k] = f_7 * ff_22[k]
                  + f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_y[k] * gf_28[k];

        t_59[k] = f_7 * ff_23[k]
                  + pb_y[k] * gf_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pb_x, dg0_32, dg1_32, ff_27, fg_29, \
                         fg_35, fg_38, gd0_15, gd1_15, gf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_5 * dg0_32[k]
                  - f_6 * dg1_32[k]
                  + pa_y[k] * fg_29[k];

        t_61[k] = f_0 * ff_27[k]
                  + pa_y[k] * fg_35[k];

        t_62[k] = pa_y[k] * fg_38[k];

        t_63[k] = f_1 * gd0_15[k]
                  - f_2 * gd1_15[k]
                  + pb_x[k] * gf_30[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, pb_x, pb_y, gd0_16, gd0_17, gd1_16, \
                         gd1_17, gf_31, gf_32, gf_33, gf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * gd0_16[k]
                  - f_4 * gd1_16[k]
                  + pb_x[k] * gf_31[k];

        t_65[k] = f_3 * gd0_17[k]
                  - f_4 * gd1_17[k]
                  + pb_x[k] * gf_32[k];

        t_66[k] = pb_x[k] * gf_33[k];

        t_67[k] = pb_x[k] * gf_35[k];

        t_68[k] = f_1 * gd0_16[k]
                  - f_2 * gd1_16[k]
                  + pb_y[k] * gf_33[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_y, pb_z, ff_29, gd0_17, gd1_17, gf_34, \
                         gf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * gd0_17[k]
                  - f_4 * gd1_17[k]
                  + pb_y[k] * gf_34[k];

        t_70[k] = pb_y[k] * gf_35[k];

        t_71[k] = f_0 * ff_29[k]
                  + f_1 * gd0_17[k]
                  - f_2 * gd1_17[k]
                  + pb_z[k] * gf_35[k];
    }
}

auto
compute_prim_gg_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg0, const size_t dg1,
                                     const size_t ff, const size_t fg, const size_t gd0,
                                     const size_t gd1, const size_t gf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_1 = buffer.data(dg1 + 1);
    const auto *dg1_2 = buffer.data(dg1 + 2);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_40 = buffer.data(ff + 40);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_17 = buffer.data(gd0 + 17);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_19 = buffer.data(gd0 + 19);
    const auto *gd0_21 = buffer.data(gd0 + 21);
    const auto *gd0_22 = buffer.data(gd0 + 22);
    const auto *gd0_23 = buffer.data(gd0 + 23);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_22 = buffer.data(gd1 + 22);
    const auto *gd1_23 = buffer.data(gd1 + 23);
    const auto *gd1_24 = buffer.data(gd1 + 24);
    const auto *gd1_28 = buffer.data(gd1 + 28);
    const auto *gd1_29 = buffer.data(gd1 + 29);
    const auto *gd1_30 = buffer.data(gd1 + 30);
    const auto *gd1_33 = buffer.data(gd1 + 33);
    const auto *gd1_34 = buffer.data(gd1 + 34);
    const auto *gd1_35 = buffer.data(gd1 + 35);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, ff_0, ff_3, gd0_0, gd1_0, gf_0, \
                         gf_1, gf_2, gf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_2[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];

        t_3[k] = f_0 * ff_3[k]
                 + pb_x[k] * gf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, ff_5, gd0_1, gd0_2, gd1_1, \
                         gd1_2, gf_3, gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_0 * ff_5[k]
                 + pb_x[k] * gf_5[k];

        t_5[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];

        t_6[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_4[k];

        t_7[k] = f_1 * gd0_2[k]
                 - f_2 * gd1_2[k]
                 + pb_z[k] * gf_5[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_y, pb_x, pb_y, ff_0, ff_1, ff_7, fg_0, fg_1, \
                         gf_6, gf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = pa_y[k] * fg_0[k];

        t_9[k] = f_5 * ff_0[k]
                 + pb_y[k] * gf_6[k];

        t_10[k] = f_6 * ff_1[k]
                  + pa_y[k] * fg_1[k];

        t_11[k] = f_7 * ff_7[k]
                  + pb_x[k] * gf_7[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_y, pa_z, pb_z, ff_0, ff_2, ff_3, fg_0, \
                         fg_2, fg_3, gf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * ff_3[k]
                  + pa_y[k] * fg_3[k];

        t_13[k] = pa_z[k] * fg_0[k];

        t_14[k] = f_5 * ff_0[k]
                  + pb_z[k] * gf_8[k];

        t_15[k] = f_6 * ff_2[k]
                  + pa_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_y, pa_z, pb_x, dg0_0, dg1_0, ff_4, ff_5, \
                         ff_10, fg_4, fg_5, fg_6, gf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * ff_10[k]
                  + pb_x[k] * gf_10[k];

        t_17[k] = f_6 * ff_4[k]
                  + pa_z[k] * fg_4[k];

        t_18[k] = f_0 * ff_5[k]
                  + pa_z[k] * fg_5[k];

        t_19[k] = f_8 * dg0_0[k]
                  - f_9 * dg1_0[k]
                  + pa_y[k] * fg_6[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_x, pb_y, pb_z, ff_6, ff_12, gd0_5, gd0_6, gd1_9, \
                         gd1_10, gf_11, gf_12, gf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_6 * ff_6[k]
                  + pb_y[k] * gf_11[k];

        t_21[k] = f_6 * ff_12[k]
                  + f_3 * gd0_6[k]
                  - f_4 * gd1_10[k]
                  + pb_x[k] * gf_13[k];

        t_22[k] = f_3 * gd0_5[k]
                  - f_4 * gd1_9[k]
                  + pb_z[k] * gf_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pb_x, pb_z, dg0_1, dg1_1, ff_13, fg_8, gd0_6, \
                         gd1_10, gf_14, gf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * ff_13[k]
                  + pb_x[k] * gf_14[k];

        t_24[k] = f_8 * dg0_1[k]
                  - f_9 * dg1_1[k]
                  + pa_x[k] * fg_8[k];

        t_25[k] = f_3 * gd0_6[k]
                  - f_4 * gd1_10[k]
                  + pb_z[k] * gf_15[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_z, dg0_0, dg1_0, ff_8, fg_7, gd0_7, \
                         gd1_11, gf_16, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * gd0_7[k]
                  - f_2 * gd1_11[k]
                  + pb_z[k] * gf_16[k];

        t_27[k] = f_8 * dg0_0[k]
                  - f_9 * dg1_0[k]
                  + pa_z[k] * fg_7[k];

        t_28[k] = f_6 * ff_8[k]
                  + pb_z[k] * gf_17[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, pb_y, ff_16, ff_17, gd0_8, gd0_10, gd1_14, \
                         gd1_16, gf_18, gf_20, gf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_3 * gd0_8[k]
                  - f_4 * gd1_14[k]
                  + pb_y[k] * gf_18[k];

        t_30[k] = f_6 * ff_16[k]
                  + f_3 * gd0_10[k]
                  - f_4 * gd1_16[k]
                  + pb_x[k] * gf_20[k];

        t_31[k] = f_6 * ff_17[k]
                  + pb_x[k] * gf_23[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_x, pb_y, dg0_2, dg1_2, fg_9, gd0_9, gd0_10, \
                         gd1_15, gd1_16, gf_21, gf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * gd0_9[k]
                  - f_2 * gd1_15[k]
                  + pb_y[k] * gf_21[k];

        t_33[k] = f_3 * gd0_10[k]
                  - f_4 * gd1_16[k]
                  + pb_y[k] * gf_22[k];

        t_34[k] = f_8 * dg0_2[k]
                  - f_9 * dg1_2[k]
                  + pa_x[k] * fg_9[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_x, pb_y, ff_11, ff_18, ff_20, ff_21, \
                         fg_10, fg_11, gf_24, gf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * ff_18[k]
                  + pa_x[k] * fg_10[k];

        t_36[k] = f_7 * ff_11[k]
                  + pb_y[k] * gf_24[k];

        t_37[k] = f_6 * ff_20[k]
                  + pa_x[k] * fg_11[k];

        t_38[k] = f_5 * ff_21[k]
                  + pb_x[k] * gf_25[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pb_z, ff_14, ff_33, ff_36, fg_13, \
                         fg_18, fg_20, gf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pa_x[k] * fg_13[k];

        t_40[k] = f_0 * ff_33[k]
                  + pa_x[k] * fg_18[k];

        t_41[k] = f_7 * ff_14[k]
                  + pb_z[k] * gf_26[k];

        t_42[k] = f_6 * ff_36[k]
                  + pa_x[k] * fg_20[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pb_x, pb_y, ff_18, ff_40, fg_23, \
                         gd0_13, gd1_22, gf_28, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_5 * ff_40[k]
                  + pb_x[k] * gf_28[k];

        t_44[k] = pa_x[k] * fg_23[k];

        t_45[k] = f_1 * gd0_13[k]
                  - f_2 * gd1_22[k]
                  + pb_x[k] * gf_29[k];

        t_46[k] = f_0 * ff_18[k]
                  + pb_y[k] * gf_29[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_x, pb_y, pb_z, ff_21, gd0_14, gd0_15, \
                         gd1_23, gd1_24, gf_30, gf_31, gf_32, gf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_23[k]
                  + pb_x[k] * gf_30[k];

        t_48[k] = f_3 * gd0_15[k]
                  - f_4 * gd1_24[k]
                  + pb_x[k] * gf_31[k];

        t_49[k] = f_0 * ff_21[k]
                  + f_1 * gd0_14[k]
                  - f_2 * gd1_23[k]
                  + pb_y[k] * gf_32[k];

        t_50[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_23[k]
                  + pb_z[k] * gf_33[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_y, pb_z, ff_19, ff_24, fg_12, fg_13, \
                         gd0_15, gd1_24, gf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * ff_24[k]
                  + pb_y[k] * gf_35[k];

        t_52[k] = f_1 * gd0_15[k]
                  - f_2 * gd1_24[k]
                  + pb_z[k] * gf_35[k];

        t_53[k] = f_6 * ff_19[k]
                  + pa_z[k] * fg_12[k];

        t_54[k] = pa_z[k] * fg_13[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_z, pb_y, pb_z, ff_21, ff_22, ff_24, ff_28, \
                         fg_14, fg_15, gf_36, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_5 * ff_21[k]
                  + pb_z[k] * gf_36[k];

        t_56[k] = f_6 * ff_22[k]
                  + pa_z[k] * fg_14[k];

        t_57[k] = f_7 * ff_28[k]
                  + pb_y[k] * gf_39[k];

        t_58[k] = f_0 * ff_24[k]
                  + pa_z[k] * fg_15[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pb_x, gd0_17, gd0_18, gd0_19, gd1_28, gd1_29, \
                         gd1_30, gf_40, gf_41, gf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * gd0_17[k]
                  - f_2 * gd1_28[k]
                  + pb_x[k] * gf_40[k];

        t_60[k] = f_3 * gd0_18[k]
                  - f_4 * gd1_29[k]
                  + pb_x[k] * gf_41[k];

        t_61[k] = f_3 * gd0_19[k]
                  - f_4 * gd1_30[k]
                  + pb_x[k] * gf_42[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pa_z, pb_y, pb_z, dg0_1, dg1_1, ff_25, ff_31, \
                         fg_16, gd0_19, gd1_30, gf_43, gf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_8 * dg0_1[k]
                  - f_9 * dg1_1[k]
                  + pa_z[k] * fg_16[k];

        t_63[k] = f_6 * ff_25[k]
                  + pb_z[k] * gf_43[k];

        t_64[k] = f_6 * ff_31[k]
                  + f_3 * gd0_19[k]
                  - f_4 * gd1_30[k]
                  + pb_y[k] * gf_45[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_y, pb_y, dg0_2, dg1_2, ff_32, ff_34, \
                         ff_37, fg_17, fg_19, fg_21, gf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_6 * ff_32[k]
                  + pb_y[k] * gf_46[k];

        t_66[k] = f_8 * dg0_2[k]
                  - f_9 * dg1_2[k]
                  + pa_y[k] * fg_17[k];

        t_67[k] = f_6 * ff_34[k]
                  + pa_y[k] * fg_19[k];

        t_68[k] = f_0 * ff_37[k]
                  + pa_y[k] * fg_21[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_y, pb_y, pb_z, ff_29, ff_39, ff_40, fg_22, \
                         fg_23, gf_47, gf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_7 * ff_29[k]
                  + pb_z[k] * gf_47[k];

        t_70[k] = f_6 * ff_39[k]
                  + pa_y[k] * fg_22[k];

        t_71[k] = f_5 * ff_40[k]
                  + pb_y[k] * gf_50[k];

        t_72[k] = pa_y[k] * fg_23[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pb_x, pb_z, ff_33, gd0_21, gd0_22, gd0_23, \
                         gd1_33, gd1_34, gd1_35, gf_51, gf_53, gf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_1 * gd0_21[k]
                  - f_2 * gd1_33[k]
                  + pb_x[k] * gf_51[k];

        t_74[k] = f_0 * ff_33[k]
                  + pb_z[k] * gf_51[k];

        t_75[k] = f_3 * gd0_22[k]
                  - f_4 * gd1_34[k]
                  + pb_x[k] * gf_53[k];

        t_76[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_35[k]
                  + pb_x[k] * gf_54[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pb_y, pb_z, ff_37, ff_40, gd0_22, gd0_23, \
                         gd1_34, gd1_35, gf_55, gf_57, gf_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_1 * gd0_22[k]
                  - f_2 * gd1_34[k]
                  + pb_y[k] * gf_55[k];

        t_78[k] = f_0 * ff_37[k]
                  + pb_z[k] * gf_55[k];

        t_79[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_35[k]
                  + pb_y[k] * gf_57[k];

        t_80[k] = f_0 * ff_40[k]
                  + f_1 * gd0_23[k]
                  - f_2 * gd1_35[k]
                  + pb_z[k] * gf_58[k];
    }
}

auto
compute_prim_gg_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg0, const size_t dg1,
                                     const size_t ff, const size_t fg, const size_t gd0,
                                     const size_t gd1, const size_t gf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / p;
    const auto f_7 = 0.5 / alpha;
    const auto f_8 = 0.5 * beta / (alpha * p);
    const auto f_9 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_1 = buffer.data(dg1 + 1);
    const auto *dg1_2 = buffer.data(dg1 + 2);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_41 = buffer.data(ff + 41);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_31 = buffer.data(fg + 31);
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
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_58 = buffer.data(fg + 58);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_19 = buffer.data(gd0 + 19);
    const auto *gd0_20 = buffer.data(gd0 + 20);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_21 = buffer.data(gd1 + 21);
    const auto *gd1_22 = buffer.data(gd1 + 22);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_62 = buffer.data(gf + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ff_0, gd0_0, gd1_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, gd0_1, gd0_2, gd1_1, gd1_2, \
                         gf_3, gf_5, gf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];

        t_6[k] = pb_z[k] * gf_3[k];

        t_7[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_5[k];

        t_8[k] = pb_y[k] * gf_6[k];

        t_9[k] = f_1 * gd0_2[k]
                 - f_2 * gd1_2[k]
                 + pb_z[k] * gf_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, ff_1, ff_3, ff_5, fg_0, fg_3, \
                         fg_4, fg_5, fg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * fg_0[k];

        t_11[k] = f_5 * ff_1[k]
                  + pa_y[k] * fg_3[k];

        t_12[k] = pa_y[k] * fg_4[k];

        t_13[k] = f_0 * ff_3[k]
                  + pa_y[k] * fg_5[k];

        t_14[k] = f_5 * ff_5[k]
                  + pa_y[k] * fg_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pa_z, pb_y, pb_z, ff_0, ff_6, \
                         fg_0, fg_3, fg_8, gf_9, gf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * ff_6[k]
                  + pb_y[k] * gf_9[k];

        t_16[k] = pa_y[k] * fg_8[k];

        t_17[k] = pa_z[k] * fg_0[k];

        t_18[k] = f_6 * ff_0[k]
                  + pb_z[k] * gf_10[k];

        t_19[k] = pa_z[k] * fg_3[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, ff_2, ff_3, ff_4, \
                         fg_4, fg_5, fg_7, gf_11, gf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * ff_2[k]
                  + pa_z[k] * fg_4[k];

        t_21[k] = pa_z[k] * fg_5[k];

        t_22[k] = f_6 * ff_3[k]
                  + pb_z[k] * gf_11[k];

        t_23[k] = f_5 * ff_4[k]
                  + pa_z[k] * fg_7[k];

        t_24[k] = pb_y[k] * gf_13[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pb_z, dg0_0, dg1_0, ff_6, fg_8, fg_9, \
                         gf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * ff_6[k]
                  + pa_z[k] * fg_8[k];

        t_26[k] = f_7 * dg0_0[k]
                  - f_8 * dg1_0[k]
                  + pa_y[k] * fg_9[k];

        t_27[k] = pb_z[k] * gf_14[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pb_x, pb_z, ff_15, ff_16, gd0_3, gd0_4, gd1_5, \
                         gd1_6, gf_15, gf_16, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_5 * ff_15[k]
                  + f_3 * gd0_4[k]
                  - f_4 * gd1_6[k]
                  + pb_x[k] * gf_16[k];

        t_29[k] = f_3 * gd0_3[k]
                  - f_4 * gd1_5[k]
                  + pb_z[k] * gf_15[k];

        t_30[k] = f_5 * ff_16[k]
                  + pb_x[k] * gf_17[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pb_y, pb_z, dg0_1, dg1_1, ff_9, fg_20, \
                         gd0_4, gd1_6, gf_17, gf_18, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_7 * dg0_1[k]
                  - f_8 * dg1_1[k]
                  + pa_x[k] * fg_20[k];

        t_32[k] = pb_z[k] * gf_17[k];

        t_33[k] = f_3 * gd0_4[k]
                  - f_4 * gd1_6[k]
                  + pb_z[k] * gf_18[k];

        t_34[k] = f_5 * ff_9[k]
                  + pb_y[k] * gf_19[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_y, pa_z, pb_z, fg_10, fg_11, fg_13, \
                         fg_14, gd0_5, gd1_7, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * gd0_5[k]
                  - f_2 * gd1_7[k]
                  + pb_z[k] * gf_19[k];

        t_36[k] = pa_y[k] * fg_13[k];

        t_37[k] = pa_z[k] * fg_10[k];

        t_38[k] = pa_y[k] * fg_14[k];

        t_39[k] = pa_z[k] * fg_11[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pb_y, pb_z, ff_8, ff_12, ff_13, fg_15, \
                         fg_16, gf_20, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_6 * ff_8[k]
                  + pb_z[k] * gf_20[k];

        t_41[k] = f_5 * ff_12[k]
                  + pa_y[k] * fg_15[k];

        t_42[k] = f_6 * ff_13[k]
                  + pb_y[k] * gf_21[k];

        t_43[k] = pa_y[k] * fg_16[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pb_y, pb_z, dg0_0, dg1_0, ff_10, fg_12, \
                         gd0_6, gd1_8, gf_22, gf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_7 * dg0_0[k]
                  - f_8 * dg1_0[k]
                  + pa_z[k] * fg_12[k];

        t_45[k] = pb_y[k] * gf_22[k];

        t_46[k] = f_5 * ff_10[k]
                  + pb_z[k] * gf_22[k];

        t_47[k] = f_3 * gd0_6[k]
                  - f_4 * gd1_8[k]
                  + pb_y[k] * gf_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_x, pb_y, ff_18, ff_19, gd0_7, gd0_8, gd1_9, \
                         gd1_10, gf_24, gf_25, gf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * ff_18[k]
                  + f_3 * gd0_8[k]
                  - f_4 * gd1_10[k]
                  + pb_x[k] * gf_24[k];

        t_49[k] = f_5 * ff_19[k]
                  + pb_x[k] * gf_27[k];

        t_50[k] = f_1 * gd0_7[k]
                  - f_2 * gd1_9[k]
                  + pb_y[k] * gf_25[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_x, pb_y, pb_z, dg0_2, dg1_2, ff_11, fg_25, \
                         gd0_8, gd1_10, gf_25, gf_26, gf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_5 * ff_11[k]
                  + pb_z[k] * gf_25[k];

        t_52[k] = f_3 * gd0_8[k]
                  - f_4 * gd1_10[k]
                  + pb_y[k] * gf_26[k];

        t_53[k] = pb_y[k] * gf_27[k];

        t_54[k] = f_7 * dg0_2[k]
                  - f_8 * dg1_2[k]
                  + pa_x[k] * fg_25[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_x, pb_x, ff_20, ff_22, ff_23, ff_24, \
                         fg_26, fg_28, fg_29, fg_31, gf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_0 * ff_20[k]
                  + pa_x[k] * fg_26[k];

        t_56[k] = f_5 * ff_22[k]
                  + pa_x[k] * fg_28[k];

        t_57[k] = f_5 * ff_23[k]
                  + pa_x[k] * fg_29[k];

        t_58[k] = f_6 * ff_24[k]
                  + pb_x[k] * gf_30[k];

        t_59[k] = pa_x[k] * fg_31[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_x, pa_z, pb_z, ff_14, fg_17, \
                         fg_18, fg_33, fg_34, fg_35, gf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_x[k] * fg_33[k];

        t_61[k] = pa_x[k] * fg_34[k];

        t_62[k] = pa_x[k] * fg_35[k];

        t_63[k] = pa_z[k] * fg_17[k];

        t_64[k] = f_6 * ff_14[k]
                  + pb_z[k] * gf_31[k];

        t_65[k] = pa_z[k] * fg_18[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, pa_x, pa_y, ff_28, fg_21, fg_36, \
                         fg_38, fg_39, fg_40, fg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_5 * ff_28[k]
                  + pa_x[k] * fg_36[k];

        t_67[k] = pa_x[k] * fg_38[k];

        t_68[k] = pa_x[k] * fg_39[k];

        t_69[k] = pa_x[k] * fg_40[k];

        t_70[k] = pa_x[k] * fg_41[k];

        t_71[k] = pa_y[k] * fg_21[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, pa_x, pa_y, ff_31, fg_22, fg_23, \
                         fg_42, fg_43, fg_44, fg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_y[k] * fg_22[k];

        t_73[k] = f_5 * ff_31[k]
                  + pa_x[k] * fg_42[k];

        t_74[k] = pa_y[k] * fg_23[k];

        t_75[k] = pa_x[k] * fg_43[k];

        t_76[k] = pa_x[k] * fg_44[k];

        t_77[k] = pa_x[k] * fg_45[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_x, pb_z, ff_17, ff_35, ff_37, ff_38, \
                         fg_46, fg_48, fg_51, fg_52, gf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_x[k] * fg_46[k];

        t_79[k] = f_0 * ff_35[k]
                  + pa_x[k] * fg_48[k];

        t_80[k] = f_9 * ff_17[k]
                  + pb_z[k] * gf_32[k];

        t_81[k] = f_5 * ff_37[k]
                  + pa_x[k] * fg_51[k];

        t_82[k] = f_5 * ff_38[k]
                  + pa_x[k] * fg_52[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_x, pb_x, ff_41, fg_54, fg_55, fg_56, \
                         fg_58, gf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_6 * ff_41[k]
                  + pb_x[k] * gf_34[k];

        t_84[k] = pa_x[k] * fg_54[k];

        t_85[k] = pa_x[k] * fg_55[k];

        t_86[k] = pa_x[k] * fg_56[k];

        t_87[k] = pa_x[k] * fg_58[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pb_x, pb_z, gd0_11, gd0_12, gd0_13, gd1_13, \
                         gd1_14, gd1_15, gf_35, gf_37, gf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_1 * gd0_11[k]
                  - f_2 * gd1_13[k]
                  + pb_x[k] * gf_35[k];

        t_89[k] = pb_z[k] * gf_35[k];

        t_90[k] = f_3 * gd0_12[k]
                  - f_4 * gd1_14[k]
                  + pb_x[k] * gf_37[k];

        t_91[k] = f_3 * gd0_13[k]
                  - f_4 * gd1_15[k]
                  + pb_x[k] * gf_38[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, t_97, pb_x, pb_y, pb_z, ff_24, ff_26, \
                         gd0_12, gd1_14, gf_39, gf_40, gf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_x[k] * gf_39[k];

        t_93[k] = pb_x[k] * gf_41[k];

        t_94[k] = f_0 * ff_24[k]
                  + f_1 * gd0_12[k]
                  - f_2 * gd1_14[k]
                  + pb_y[k] * gf_39[k];

        t_95[k] = pb_z[k] * gf_39[k];

        t_96[k] = f_3 * gd0_12[k]
                  - f_4 * gd1_14[k]
                  + pb_z[k] * gf_40[k];

        t_97[k] = f_0 * ff_26[k]
                  + pb_y[k] * gf_41[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, pa_z, pb_z, ff_20, ff_21, fg_26, \
                         fg_28, fg_29, gd0_13, gd1_15, gf_41, gf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_1 * gd0_13[k]
                  - f_2 * gd1_15[k]
                  + pb_z[k] * gf_41[k];

        t_99[k] = pa_z[k] * fg_26[k];

        t_100[k] = f_6 * ff_20[k]
                   + pb_z[k] * gf_42[k];

        t_101[k] = pa_z[k] * fg_28[k];

        t_102[k] = f_5 * ff_21[k]
                   + pa_z[k] * fg_29[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_x, pb_y, pb_z, ff_24, \
                         ff_25, ff_30, fg_31, fg_33, gf_44, gf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = pb_x[k] * gf_45[k];

        t_104[k] = pa_z[k] * fg_31[k];

        t_105[k] = f_6 * ff_24[k]
                   + pb_z[k] * gf_44[k];

        t_106[k] = f_5 * ff_25[k]
                   + pa_z[k] * fg_33[k];

        t_107[k] = f_9 * ff_30[k]
                   + pb_y[k] * gf_45[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pa_z, pb_x, pb_z, ff_26, ff_27, fg_35, \
                         gd0_14, gd0_15, gd1_17, gd1_18, gf_46, gf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_0 * ff_26[k]
                   + pa_z[k] * fg_35[k];

        t_109[k] = f_1 * gd0_14[k]
                   - f_2 * gd1_17[k]
                   + pb_x[k] * gf_46[k];

        t_110[k] = f_5 * ff_27[k]
                   + pb_z[k] * gf_46[k];

        t_111[k] = f_3 * gd0_15[k]
                   - f_4 * gd1_18[k]
                   + pb_x[k] * gf_47[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pa_z, pb_x, dg0_1, dg1_1, fg_37, gd0_16, \
                         gd1_19, gf_48, gf_49, gf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_3 * gd0_16[k]
                   - f_4 * gd1_19[k]
                   + pb_x[k] * gf_48[k];

        t_113[k] = pb_x[k] * gf_49[k];

        t_114[k] = pb_x[k] * gf_51[k];

        t_115[k] = f_7 * dg0_1[k]
                   - f_8 * dg1_1[k]
                   + pa_z[k] * fg_37[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pb_y, pb_z, ff_29, ff_33, ff_34, gd0_16, gd1_19, \
                         gf_49, gf_50, gf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_5 * ff_29[k]
                   + pb_z[k] * gf_49[k];

        t_117[k] = f_5 * ff_33[k]
                   + f_3 * gd0_16[k]
                   - f_4 * gd1_19[k]
                   + pb_y[k] * gf_50[k];

        t_118[k] = f_5 * ff_34[k]
                   + pb_y[k] * gf_51[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, pa_y, dg0_2, dg1_2, ff_36, fg_47, \
                         fg_48, fg_50, fg_51, fg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_7 * dg0_2[k]
                   - f_8 * dg1_2[k]
                   + pa_y[k] * fg_47[k];

        t_120[k] = pa_y[k] * fg_48[k];

        t_121[k] = pa_y[k] * fg_50[k];

        t_122[k] = f_5 * ff_36[k]
                   + pa_y[k] * fg_51[k];

        t_123[k] = pa_y[k] * fg_52[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pb_x, pb_z, ff_32, ff_39, ff_40, \
                         fg_54, fg_56, gf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_x[k] * gf_53[k];

        t_125[k] = f_0 * ff_39[k]
                   + pa_y[k] * fg_54[k];

        t_126[k] = f_9 * ff_32[k]
                   + pb_z[k] * gf_53[k];

        t_127[k] = f_5 * ff_40[k]
                   + pa_y[k] * fg_56[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pa_y, pb_x, pb_y, pb_z, ff_35, \
                         ff_41, fg_58, gd0_18, gd1_21, gf_55, gf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_6 * ff_41[k]
                   + pb_y[k] * gf_55[k];

        t_129[k] = pa_y[k] * fg_58[k];

        t_130[k] = f_1 * gd0_18[k]
                   - f_2 * gd1_21[k]
                   + pb_x[k] * gf_56[k];

        t_131[k] = pb_y[k] * gf_56[k];

        t_132[k] = f_0 * ff_35[k]
                   + pb_z[k] * gf_56[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, pb_x, pb_y, gd0_19, gd0_20, \
                         gd1_22, gd1_23, gf_58, gf_59, gf_60, gf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_3 * gd0_19[k]
                   - f_4 * gd1_22[k]
                   + pb_x[k] * gf_58[k];

        t_134[k] = f_3 * gd0_20[k]
                   - f_4 * gd1_23[k]
                   + pb_x[k] * gf_59[k];

        t_135[k] = pb_x[k] * gf_60[k];

        t_136[k] = pb_x[k] * gf_62[k];

        t_137[k] = f_1 * gd0_19[k]
                   - f_2 * gd1_22[k]
                   + pb_y[k] * gf_60[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pb_y, pb_z, ff_39, ff_41, gd0_20, gd1_23, \
                         gf_60, gf_61, gf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_0 * ff_39[k]
                   + pb_z[k] * gf_60[k];

        t_139[k] = f_3 * gd0_20[k]
                   - f_4 * gd1_23[k]
                   + pb_y[k] * gf_61[k];

        t_140[k] = pb_y[k] * gf_62[k];

        t_141[k] = f_0 * ff_41[k]
                   + f_1 * gd0_20[k]
                   - f_2 * gd1_23[k]
                   + pb_z[k] * gf_62[k];
    }
}

auto
compute_prim_gg_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg0, const size_t dg1,
                                     const size_t ff, const size_t fg, const size_t gd0,
                                     const size_t gd1, const size_t gf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / p;
    const auto f_7 = 0.5 / alpha;
    const auto f_8 = 0.5 * beta / (alpha * p);
    const auto f_9 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_9 = buffer.data(dg1 + 9);
    const auto *dg1_17 = buffer.data(dg1 + 17);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_32 = buffer.data(fg + 32);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_17 = buffer.data(gd0 + 17);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_19 = buffer.data(gd0 + 19);
    const auto *gd0_21 = buffer.data(gd0 + 21);
    const auto *gd0_22 = buffer.data(gd0 + 22);
    const auto *gd0_23 = buffer.data(gd0 + 23);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_20 = buffer.data(gd1 + 20);
    const auto *gd1_21 = buffer.data(gd1 + 21);
    const auto *gd1_24 = buffer.data(gd1 + 24);
    const auto *gd1_25 = buffer.data(gd1 + 25);
    const auto *gd1_26 = buffer.data(gd1 + 26);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ff_0, gd0_0, gd1_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, fg_0, gd0_1, gd0_2, gd1_1, \
                         gd1_2, gf_3, gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];

        t_6[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_4[k];

        t_7[k] = pb_y[k] * gf_5[k];

        t_8[k] = f_1 * gd0_2[k]
                 - f_2 * gd1_2[k]
                 + pb_z[k] * gf_5[k];

        t_9[k] = pa_y[k] * fg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_z, ff_0, ff_1, ff_3, fg_0, \
                         fg_3, fg_5, gf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * ff_1[k]
                  + pa_y[k] * fg_3[k];

        t_11[k] = f_0 * ff_3[k]
                  + pa_y[k] * fg_5[k];

        t_12[k] = pa_z[k] * fg_0[k];

        t_13[k] = f_6 * ff_0[k]
                  + pb_z[k] * gf_8[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pa_z, dg0_0, dg1_0, ff_2, ff_4, ff_6, \
                         fg_4, fg_6, fg_8, fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * ff_2[k]
                  + pa_z[k] * fg_4[k];

        t_15[k] = f_5 * ff_4[k]
                  + pa_z[k] * fg_6[k];

        t_16[k] = f_0 * ff_6[k]
                  + pa_z[k] * fg_8[k];

        t_17[k] = f_7 * dg0_0[k]
                  - f_8 * dg1_0[k]
                  + pa_y[k] * fg_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_x, pb_z, ff_12, ff_13, gd0_5, gd0_6, \
                         gd1_6, gd1_7, gf_10, gf_11, gf_12, gf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_z[k] * gf_10[k];

        t_19[k] = f_5 * ff_12[k]
                  + f_3 * gd0_6[k]
                  - f_4 * gd1_7[k]
                  + pb_x[k] * gf_12[k];

        t_20[k] = f_3 * gd0_5[k]
                  - f_4 * gd1_6[k]
                  + pb_z[k] * gf_11[k];

        t_21[k] = f_5 * ff_13[k]
                  + pb_x[k] * gf_13[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_z, dg0_1, dg1_9, fg_11, gd0_6, \
                         gd0_7, gd1_7, gd1_8, gf_13, gf_14, gf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_7 * dg0_1[k]
                  - f_8 * dg1_9[k]
                  + pa_x[k] * fg_11[k];

        t_23[k] = pb_z[k] * gf_13[k];

        t_24[k] = f_3 * gd0_6[k]
                  - f_4 * gd1_7[k]
                  + pb_z[k] * gf_14[k];

        t_25[k] = f_1 * gd0_7[k]
                  - f_2 * gd1_8[k]
                  + pb_z[k] * gf_15[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_y, pb_z, dg0_0, dg1_0, ff_9, fg_10, \
                         gd0_8, gd1_9, gf_16, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * dg0_0[k]
                  - f_8 * dg1_0[k]
                  + pa_z[k] * fg_10[k];

        t_27[k] = pb_y[k] * gf_16[k];

        t_28[k] = f_5 * ff_9[k]
                  + pb_z[k] * gf_16[k];

        t_29[k] = f_3 * gd0_8[k]
                  - f_4 * gd1_9[k]
                  + pb_y[k] * gf_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_x, pb_y, ff_15, ff_16, gd0_9, gd0_10, \
                         gd1_10, gd1_11, gf_18, gf_19, gf_20, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * ff_15[k]
                  + f_3 * gd0_10[k]
                  - f_4 * gd1_11[k]
                  + pb_x[k] * gf_18[k];

        t_31[k] = f_5 * ff_16[k]
                  + pb_x[k] * gf_21[k];

        t_32[k] = f_1 * gd0_9[k]
                  - f_2 * gd1_10[k]
                  + pb_y[k] * gf_19[k];

        t_33[k] = f_3 * gd0_10[k]
                  - f_4 * gd1_11[k]
                  + pb_y[k] * gf_20[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_x, pb_y, dg0_2, dg1_17, ff_17, \
                         ff_19, fg_12, fg_13, fg_14, fg_18, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * gf_21[k];

        t_35[k] = f_7 * dg0_2[k]
                  - f_8 * dg1_17[k]
                  + pa_x[k] * fg_12[k];

        t_36[k] = f_0 * ff_17[k]
                  + pa_x[k] * fg_13[k];

        t_37[k] = f_5 * ff_19[k]
                  + pa_x[k] * fg_14[k];

        t_38[k] = pa_x[k] * fg_18[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pb_z, ff_14, ff_29, ff_32, fg_24, \
                         fg_26, fg_32, gf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * ff_29[k]
                  + pa_x[k] * fg_24[k];

        t_40[k] = f_9 * ff_14[k]
                  + pb_z[k] * gf_24[k];

        t_41[k] = f_5 * ff_32[k]
                  + pa_x[k] * fg_26[k];

        t_42[k] = pa_x[k] * fg_32[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pb_x, gd0_13, gd0_14, gd0_15, gd1_14, gd1_15, \
                         gd1_16, gf_26, gf_27, gf_28, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * gd0_13[k]
                  - f_2 * gd1_14[k]
                  + pb_x[k] * gf_26[k];

        t_44[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_15[k]
                  + pb_x[k] * gf_27[k];

        t_45[k] = f_3 * gd0_15[k]
                  - f_4 * gd1_16[k]
                  + pb_x[k] * gf_28[k];

        t_46[k] = pb_x[k] * gf_29[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pb_x, pb_y, pb_z, ff_21, ff_23, gd0_14, \
                         gd1_15, gf_29, gf_30, gf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_x[k] * gf_31[k];

        t_48[k] = f_0 * ff_21[k]
                  + f_1 * gd0_14[k]
                  - f_2 * gd1_15[k]
                  + pb_y[k] * gf_29[k];

        t_49[k] = pb_z[k] * gf_29[k];

        t_50[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_15[k]
                  + pb_z[k] * gf_30[k];

        t_51[k] = f_0 * ff_23[k]
                  + pb_y[k] * gf_31[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_z, pb_z, ff_18, ff_21, fg_15, fg_18, \
                         gd0_15, gd1_16, gf_31, gf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * gd0_15[k]
                  - f_2 * gd1_16[k]
                  + pb_z[k] * gf_31[k];

        t_53[k] = f_5 * ff_18[k]
                  + pa_z[k] * fg_15[k];

        t_54[k] = pa_z[k] * fg_18[k];

        t_55[k] = f_6 * ff_21[k]
                  + pb_z[k] * gf_32[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_z, pb_x, pb_y, ff_22, ff_23, ff_25, fg_20, \
                         fg_21, gd0_17, gd1_19, gf_33, gf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_5 * ff_22[k]
                  + pa_z[k] * fg_20[k];

        t_57[k] = f_9 * ff_25[k]
                  + pb_y[k] * gf_33[k];

        t_58[k] = f_0 * ff_23[k]
                  + pa_z[k] * fg_21[k];

        t_59[k] = f_1 * gd0_17[k]
                  - f_2 * gd1_19[k]
                  + pb_x[k] * gf_34[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pb_x, gd0_18, gd0_19, gd1_20, gd1_21, gf_35, \
                         gf_36, gf_37, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_3 * gd0_18[k]
                  - f_4 * gd1_20[k]
                  + pb_x[k] * gf_35[k];

        t_61[k] = f_3 * gd0_19[k]
                  - f_4 * gd1_21[k]
                  + pb_x[k] * gf_36[k];

        t_62[k] = pb_x[k] * gf_37[k];

        t_63[k] = pb_x[k] * gf_39[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_z, pb_y, pb_z, dg0_1, dg1_9, ff_24, ff_27, \
                         fg_22, gd0_19, gd1_21, gf_37, gf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_7 * dg0_1[k]
                  - f_8 * dg1_9[k]
                  + pa_z[k] * fg_22[k];

        t_65[k] = f_5 * ff_24[k]
                  + pb_z[k] * gf_37[k];

        t_66[k] = f_5 * ff_27[k]
                  + f_3 * gd0_19[k]
                  - f_4 * gd1_21[k]
                  + pb_y[k] * gf_38[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_y, pb_y, dg0_2, dg1_17, ff_28, ff_30, \
                         ff_33, fg_23, fg_25, fg_29, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_5 * ff_28[k]
                  + pb_y[k] * gf_39[k];

        t_68[k] = f_7 * dg0_2[k]
                  - f_8 * dg1_17[k]
                  + pa_y[k] * fg_23[k];

        t_69[k] = f_5 * ff_30[k]
                  + pa_y[k] * fg_25[k];

        t_70[k] = f_0 * ff_33[k]
                  + pa_y[k] * fg_29[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, pb_z, ff_26, ff_34, ff_35, fg_30, \
                         fg_32, gf_40, gf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_9 * ff_26[k]
                  + pb_z[k] * gf_40[k];

        t_72[k] = f_5 * ff_34[k]
                  + pa_y[k] * fg_30[k];

        t_73[k] = f_6 * ff_35[k]
                  + pb_y[k] * gf_41[k];

        t_74[k] = pa_y[k] * fg_32[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pb_x, pb_z, ff_29, gd0_21, gd0_22, gd0_23, \
                         gd1_24, gd1_25, gd1_26, gf_42, gf_43, gf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * gd0_21[k]
                  - f_2 * gd1_24[k]
                  + pb_x[k] * gf_42[k];

        t_76[k] = f_0 * ff_29[k]
                  + pb_z[k] * gf_42[k];

        t_77[k] = f_3 * gd0_22[k]
                  - f_4 * gd1_25[k]
                  + pb_x[k] * gf_43[k];

        t_78[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_26[k]
                  + pb_x[k] * gf_44[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pb_x, pb_y, pb_z, ff_33, gd0_22, \
                         gd0_23, gd1_25, gd1_26, gf_45, gf_46, gf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_x[k] * gf_45[k];

        t_80[k] = pb_x[k] * gf_47[k];

        t_81[k] = f_1 * gd0_22[k]
                  - f_2 * gd1_25[k]
                  + pb_y[k] * gf_45[k];

        t_82[k] = f_0 * ff_33[k]
                  + pb_z[k] * gf_45[k];

        t_83[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_26[k]
                  + pb_y[k] * gf_46[k];
    }

#pragma omp simd aligned(t_84, t_85, pb_y, pb_z, ff_35, gd0_23, gd1_26, \
                         gf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = pb_y[k] * gf_47[k];

        t_85[k] = f_0 * ff_35[k]
                  + f_1 * gd0_23[k]
                  - f_2 * gd1_26[k]
                  + pb_z[k] * gf_47[k];
    }
}

auto
compute_prim_gg_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg0, const size_t dg1,
                                     const size_t ff, const size_t fg, const size_t gd0,
                                     const size_t gd1, const size_t gf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / p;

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

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_1 = buffer.data(dg1 + 1);
    const auto *dg1_2 = buffer.data(dg1 + 2);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_26 = buffer.data(ff + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_56 = buffer.data(fg + 56);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_17 = buffer.data(gd0 + 17);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_19 = buffer.data(gd0 + 19);
    const auto *gd0_21 = buffer.data(gd0 + 21);
    const auto *gd0_22 = buffer.data(gd0 + 22);
    const auto *gd0_23 = buffer.data(gd0 + 23);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_21 = buffer.data(gd1 + 21);
    const auto *gd1_22 = buffer.data(gd1 + 22);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ff_0, gd0_0, gd1_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, fg_0, gd0_1, gd0_2, gd1_1, \
                         gd1_2, gf_3, gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];

        t_6[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_4[k];

        t_7[k] = pb_y[k] * gf_5[k];

        t_8[k] = f_1 * gd0_2[k]
                 - f_2 * gd1_2[k]
                 + pb_z[k] * gf_5[k];

        t_9[k] = pa_y[k] * fg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pa_z, pb_z, dg0_0, dg1_0, fg_0, fg_10, \
                         gf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * fg_0[k];

        t_11[k] = f_5 * dg0_0[k]
                  - f_6 * dg1_0[k]
                  + pa_y[k] * fg_10[k];

        t_12[k] = pb_z[k] * gf_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, pb_z, ff_8, ff_9, gd0_5, gd0_6, gd1_5, gd1_6, \
                         gf_9, gf_10, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_7 * ff_8[k]
                  + f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_x[k] * gf_10[k];

        t_14[k] = f_3 * gd0_5[k]
                  - f_4 * gd1_5[k]
                  + pb_z[k] * gf_9[k];

        t_15[k] = f_7 * ff_9[k]
                  + pb_x[k] * gf_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_z, dg0_1, dg1_1, fg_21, gd0_6, \
                         gd0_7, gd1_6, gd1_7, gf_11, gf_12, gf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * dg0_1[k]
                  - f_6 * dg1_1[k]
                  + pa_x[k] * fg_21[k];

        t_17[k] = pb_z[k] * gf_11[k];

        t_18[k] = f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_z[k] * gf_12[k];

        t_19[k] = f_1 * gd0_7[k]
                  - f_2 * gd1_7[k]
                  + pb_z[k] * gf_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, dg0_0, dg1_0, fg_13, gd0_8, gd1_8, \
                         gf_14, gf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * dg0_0[k]
                  - f_6 * dg1_0[k]
                  + pa_z[k] * fg_13[k];

        t_21[k] = pb_y[k] * gf_14[k];

        t_22[k] = f_3 * gd0_8[k]
                  - f_4 * gd1_8[k]
                  + pb_y[k] * gf_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_x, pb_y, ff_10, ff_11, gd0_9, gd0_10, \
                         gd1_9, gd1_10, gf_16, gf_17, gf_18, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_7 * ff_10[k]
                  + f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_x[k] * gf_16[k];

        t_24[k] = f_7 * ff_11[k]
                  + pb_x[k] * gf_19[k];

        t_25[k] = f_1 * gd0_9[k]
                  - f_2 * gd1_9[k]
                  + pb_y[k] * gf_17[k];

        t_26[k] = f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_y[k] * gf_18[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pb_y, dg0_2, dg1_2, fg_25, fg_32, \
                         fg_56, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_y[k] * gf_19[k];

        t_28[k] = f_5 * dg0_2[k]
                  - f_6 * dg1_2[k]
                  + pa_x[k] * fg_25[k];

        t_29[k] = pa_x[k] * fg_32[k];

        t_30[k] = pa_x[k] * fg_56[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pb_x, gd0_13, gd0_14, gd0_15, gd1_13, gd1_14, \
                         gd1_15, gf_24, gf_25, gf_26, gf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_1 * gd0_13[k]
                  - f_2 * gd1_13[k]
                  + pb_x[k] * gf_24[k];

        t_32[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_x[k] * gf_25[k];

        t_33[k] = f_3 * gd0_15[k]
                  - f_4 * gd1_15[k]
                  + pb_x[k] * gf_26[k];

        t_34[k] = pb_x[k] * gf_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pb_x, pb_y, pb_z, ff_15, gd0_14, \
                         gd0_15, gd1_14, gd1_15, gf_27, gf_28, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pb_x[k] * gf_29[k];

        t_36[k] = f_0 * ff_15[k]
                  + f_1 * gd0_14[k]
                  - f_2 * gd1_14[k]
                  + pb_y[k] * gf_27[k];

        t_37[k] = pb_z[k] * gf_27[k];

        t_38[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_z[k] * gf_28[k];

        t_39[k] = f_1 * gd0_15[k]
                  - f_2 * gd1_15[k]
                  + pb_z[k] * gf_29[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pb_x, fg_32, gd0_17, gd0_18, gd0_19, \
                         gd1_17, gd1_18, gd1_19, gf_31, gf_32, gf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_z[k] * fg_32[k];

        t_41[k] = f_1 * gd0_17[k]
                  - f_2 * gd1_17[k]
                  + pb_x[k] * gf_31[k];

        t_42[k] = f_3 * gd0_18[k]
                  - f_4 * gd1_18[k]
                  + pb_x[k] * gf_32[k];

        t_43[k] = f_3 * gd0_19[k]
                  - f_4 * gd1_19[k]
                  + pb_x[k] * gf_33[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pb_x, pb_y, dg0_1, dg1_1, ff_19, fg_39, \
                         gd0_19, gd1_19, gf_34, gf_35, gf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * gf_34[k];

        t_45[k] = pb_x[k] * gf_36[k];

        t_46[k] = f_5 * dg0_1[k]
                  - f_6 * dg1_1[k]
                  + pa_z[k] * fg_39[k];

        t_47[k] = f_7 * ff_19[k]
                  + f_3 * gd0_19[k]
                  - f_4 * gd1_19[k]
                  + pb_y[k] * gf_35[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pb_x, pb_y, dg0_2, dg1_2, ff_20, fg_46, \
                         fg_56, gd0_21, gd1_21, gf_36, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_7 * ff_20[k]
                  + pb_y[k] * gf_36[k];

        t_49[k] = f_5 * dg0_2[k]
                  - f_6 * dg1_2[k]
                  + pa_y[k] * fg_46[k];

        t_50[k] = pa_y[k] * fg_56[k];

        t_51[k] = f_1 * gd0_21[k]
                  - f_2 * gd1_21[k]
                  + pb_x[k] * gf_39[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pb_x, pb_y, gd0_22, gd0_23, gd1_22, \
                         gd1_23, gf_40, gf_41, gf_42, gf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * gd0_22[k]
                  - f_4 * gd1_22[k]
                  + pb_x[k] * gf_40[k];

        t_53[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pb_x[k] * gf_41[k];

        t_54[k] = pb_x[k] * gf_42[k];

        t_55[k] = pb_x[k] * gf_44[k];

        t_56[k] = f_1 * gd0_22[k]
                  - f_2 * gd1_22[k]
                  + pb_y[k] * gf_42[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_y, pb_z, ff_26, gd0_23, gd1_23, gf_43, \
                         gf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pb_y[k] * gf_43[k];

        t_58[k] = pb_y[k] * gf_44[k];

        t_59[k] = f_0 * ff_26[k]
                  + f_1 * gd0_23[k]
                  - f_2 * gd1_23[k]
                  + pb_z[k] * gf_44[k];
    }
}

auto
compute_prim_gg_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg0, const size_t dg1,
                                     const size_t ff, const size_t fg, const size_t gd0,
                                     const size_t gd1, const size_t gf, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_15 = buffer.data(dg1 + 15);
    const auto *dg1_29 = buffer.data(dg1 + 29);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_47 = buffer.data(fg + 47);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_17 = buffer.data(gd0 + 17);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_19 = buffer.data(gd0 + 19);
    const auto *gd0_21 = buffer.data(gd0 + 21);
    const auto *gd0_22 = buffer.data(gd0 + 22);
    const auto *gd0_23 = buffer.data(gd0 + 23);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_21 = buffer.data(gd1 + 21);
    const auto *gd1_22 = buffer.data(gd1 + 22);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ff_0, gd0_0, gd1_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, gd0_1, gd0_2, gd1_1, gd1_2, \
                         gf_3, gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];

        t_6[k] = pb_z[k] * gf_3[k];

        t_7[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_4[k];

        t_8[k] = pb_y[k] * gf_5[k];

        t_9[k] = f_1 * gd0_2[k]
                 - f_2 * gd1_2[k]
                 + pb_z[k] * gf_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pa_y, pa_z, ff_2, ff_3, fg_0, \
                         fg_4, fg_5, fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * fg_0[k];

        t_11[k] = f_0 * ff_3[k]
                  + pa_y[k] * fg_5[k];

        t_12[k] = pa_y[k] * fg_9[k];

        t_13[k] = pa_z[k] * fg_0[k];

        t_14[k] = f_5 * ff_2[k]
                  + pa_z[k] * fg_4[k];

        t_15[k] = pa_z[k] * fg_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_y, pa_z, pb_y, dg0_0, dg1_0, ff_4, ff_6, \
                         fg_7, fg_9, fg_10, gf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * ff_4[k]
                  + pa_z[k] * fg_7[k];

        t_17[k] = pb_y[k] * gf_8[k];

        t_18[k] = f_0 * ff_6[k]
                  + pa_z[k] * fg_9[k];

        t_19[k] = f_6 * dg0_0[k]
                  - f_7 * dg1_0[k]
                  + pa_y[k] * fg_10[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pb_x, pb_z, ff_10, ff_11, gd0_5, gd0_6, \
                         gd1_5, gd1_6, gf_9, gf_10, gf_11, gf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_z[k] * gf_9[k];

        t_21[k] = f_5 * ff_10[k]
                  + f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_x[k] * gf_11[k];

        t_22[k] = f_3 * gd0_5[k]
                  - f_4 * gd1_5[k]
                  + pb_z[k] * gf_10[k];

        t_23[k] = f_5 * ff_11[k]
                  + pb_x[k] * gf_12[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_x, pb_z, dg0_1, dg1_15, fg_17, gd0_6, \
                         gd0_7, gd1_6, gd1_7, gf_12, gf_13, gf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_6 * dg0_1[k]
                  - f_7 * dg1_15[k]
                  + pa_x[k] * fg_17[k];

        t_25[k] = pb_z[k] * gf_12[k];

        t_26[k] = f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_z[k] * gf_13[k];

        t_27[k] = f_1 * gd0_7[k]
                  - f_2 * gd1_7[k]
                  + pb_z[k] * gf_14[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pa_z, pb_y, dg0_0, dg1_0, fg_11, fg_12, \
                         fg_13, gf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_z[k] * fg_11[k];

        t_29[k] = pa_y[k] * fg_13[k];

        t_30[k] = f_6 * dg0_0[k]
                  - f_7 * dg1_0[k]
                  + pa_z[k] * fg_12[k];

        t_31[k] = pb_y[k] * gf_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, ff_12, ff_13, gd0_8, gd0_10, gd1_8, \
                         gd1_10, gf_16, gf_17, gf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * gd0_8[k]
                  - f_4 * gd1_8[k]
                  + pb_y[k] * gf_16[k];

        t_33[k] = f_5 * ff_12[k]
                  + f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_x[k] * gf_17[k];

        t_34[k] = f_5 * ff_13[k]
                  + pb_x[k] * gf_20[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_y, dg0_2, dg1_29, fg_21, gd0_9, \
                         gd0_10, gd1_9, gd1_10, gf_18, gf_19, gf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * gd0_9[k]
                  - f_2 * gd1_9[k]
                  + pb_y[k] * gf_18[k];

        t_36[k] = f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_y[k] * gf_19[k];

        t_37[k] = pb_y[k] * gf_20[k];

        t_38[k] = f_6 * dg0_2[k]
                  - f_7 * dg1_29[k]
                  + pa_x[k] * fg_21[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_x, pa_z, pb_x, ff_14, ff_16, ff_18, \
                         fg_14, fg_22, fg_24, fg_28, gf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * ff_14[k]
                  + pa_x[k] * fg_22[k];

        t_40[k] = f_5 * ff_16[k]
                  + pa_x[k] * fg_24[k];

        t_41[k] = f_8 * ff_18[k]
                  + pb_x[k] * gf_22[k];

        t_42[k] = pa_x[k] * fg_28[k];

        t_43[k] = pa_z[k] * fg_14[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, ff_26, ff_29, ff_32, fg_38, \
                         fg_41, fg_47, gf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * ff_26[k]
                  + pa_x[k] * fg_38[k];

        t_45[k] = f_5 * ff_29[k]
                  + pa_x[k] * fg_41[k];

        t_46[k] = f_8 * ff_32[k]
                  + pb_x[k] * gf_24[k];

        t_47[k] = pa_x[k] * fg_47[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_x, pb_z, gd0_13, gd0_14, gd0_15, gd1_13, \
                         gd1_14, gd1_15, gf_25, gf_26, gf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * gd0_13[k]
                  - f_2 * gd1_13[k]
                  + pb_x[k] * gf_25[k];

        t_49[k] = pb_z[k] * gf_25[k];

        t_50[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_x[k] * gf_26[k];

        t_51[k] = f_3 * gd0_15[k]
                  - f_4 * gd1_15[k]
                  + pb_x[k] * gf_27[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pb_x, pb_y, pb_z, ff_18, gd0_14, \
                         gd1_14, gf_28, gf_29, gf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_x[k] * gf_28[k];

        t_53[k] = pb_x[k] * gf_30[k];

        t_54[k] = f_0 * ff_18[k]
                  + f_1 * gd0_14[k]
                  - f_2 * gd1_14[k]
                  + pb_y[k] * gf_28[k];

        t_55[k] = pb_z[k] * gf_28[k];

        t_56[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_z[k] * gf_29[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pa_z, pb_x, pb_z, ff_15, fg_22, fg_25, \
                         fg_28, gd0_15, gd1_15, gf_30, gf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_1 * gd0_15[k]
                  - f_2 * gd1_15[k]
                  + pb_z[k] * gf_30[k];

        t_58[k] = pa_z[k] * fg_22[k];

        t_59[k] = f_5 * ff_15[k]
                  + pa_z[k] * fg_25[k];

        t_60[k] = pb_x[k] * gf_32[k];

        t_61[k] = pa_z[k] * fg_28[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_z, pb_x, ff_19, ff_20, fg_30, fg_31, \
                         gd0_17, gd0_18, gd1_17, gd1_18, gf_33, gf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_5 * ff_19[k]
                  + pa_z[k] * fg_30[k];

        t_63[k] = f_0 * ff_20[k]
                  + pa_z[k] * fg_31[k];

        t_64[k] = f_1 * gd0_17[k]
                  - f_2 * gd1_17[k]
                  + pb_x[k] * gf_33[k];

        t_65[k] = f_3 * gd0_18[k]
                  - f_4 * gd1_18[k]
                  + pb_x[k] * gf_34[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_z, pb_x, dg0_1, dg1_15, fg_32, gd0_19, \
                         gd1_19, gf_35, gf_36, gf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_3 * gd0_19[k]
                  - f_4 * gd1_19[k]
                  + pb_x[k] * gf_35[k];

        t_67[k] = pb_x[k] * gf_36[k];

        t_68[k] = pb_x[k] * gf_38[k];

        t_69[k] = f_6 * dg0_1[k]
                  - f_7 * dg1_15[k]
                  + pa_z[k] * fg_32[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_y, pb_y, dg0_2, dg1_29, ff_24, ff_25, fg_37, \
                         gd0_19, gd1_19, gf_37, gf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_5 * ff_24[k]
                  + f_3 * gd0_19[k]
                  - f_4 * gd1_19[k]
                  + pb_y[k] * gf_37[k];

        t_71[k] = f_5 * ff_25[k]
                  + pb_y[k] * gf_38[k];

        t_72[k] = f_6 * dg0_2[k]
                  - f_7 * dg1_29[k]
                  + pa_y[k] * fg_37[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_y, pb_x, ff_27, ff_30, ff_31, fg_40, \
                         fg_44, fg_45, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_5 * ff_27[k]
                  + pa_y[k] * fg_40[k];

        t_74[k] = pb_x[k] * gf_39[k];

        t_75[k] = f_0 * ff_30[k]
                  + pa_y[k] * fg_44[k];

        t_76[k] = f_5 * ff_31[k]
                  + pa_y[k] * fg_45[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pb_x, pb_y, ff_32, fg_47, gd0_21, \
                         gd1_21, gf_41, gf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_8 * ff_32[k]
                  + pb_y[k] * gf_41[k];

        t_78[k] = pa_y[k] * fg_47[k];

        t_79[k] = f_1 * gd0_21[k]
                  - f_2 * gd1_21[k]
                  + pb_x[k] * gf_42[k];

        t_80[k] = pb_y[k] * gf_42[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, pb_x, pb_y, gd0_22, gd0_23, gd1_22, \
                         gd1_23, gf_43, gf_44, gf_45, gf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_3 * gd0_22[k]
                  - f_4 * gd1_22[k]
                  + pb_x[k] * gf_43[k];

        t_82[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pb_x[k] * gf_44[k];

        t_83[k] = pb_x[k] * gf_45[k];

        t_84[k] = pb_x[k] * gf_47[k];

        t_85[k] = f_1 * gd0_22[k]
                  - f_2 * gd1_22[k]
                  + pb_y[k] * gf_45[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pb_y, pb_z, ff_32, gd0_23, gd1_23, gf_46, \
                         gf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pb_y[k] * gf_46[k];

        t_87[k] = pb_y[k] * gf_47[k];

        t_88[k] = f_0 * ff_32[k]
                  + f_1 * gd0_23[k]
                  - f_2 * gd1_23[k]
                  + pb_z[k] * gf_47[k];
    }
}

auto
compute_prim_gg_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg0, const size_t dg1,
                                      const size_t ff, const size_t fg, const size_t gd0,
                                      const size_t gd1, const size_t gf, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_15 = buffer.data(dg0 + 15);
    const auto *dg0_29 = buffer.data(dg0 + 29);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_14 = buffer.data(dg1 + 14);
    const auto *dg1_26 = buffer.data(dg1 + 26);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_26 = buffer.data(ff + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_32 = buffer.data(fg + 32);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_17 = buffer.data(gd0 + 17);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_19 = buffer.data(gd0 + 19);
    const auto *gd0_21 = buffer.data(gd0 + 21);
    const auto *gd0_22 = buffer.data(gd0 + 22);
    const auto *gd0_23 = buffer.data(gd0 + 23);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_21 = buffer.data(gd1 + 21);
    const auto *gd1_22 = buffer.data(gd1 + 22);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ff_0, gd0_0, gd1_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, fg_0, gd0_1, gd0_2, gd1_1, \
                         gd1_2, gf_3, gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];

        t_6[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_4[k];

        t_7[k] = pb_y[k] * gf_5[k];

        t_8[k] = f_1 * gd0_2[k]
                 - f_2 * gd1_2[k]
                 + pb_z[k] * gf_5[k];

        t_9[k] = pa_y[k] * fg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, dg0_0, dg1_0, ff_3, ff_5, fg_0, \
                         fg_5, fg_8, fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * ff_3[k]
                  + pa_y[k] * fg_5[k];

        t_11[k] = pa_z[k] * fg_0[k];

        t_12[k] = f_0 * ff_5[k]
                  + pa_z[k] * fg_8[k];

        t_13[k] = f_5 * dg0_0[k]
                  - f_6 * dg1_0[k]
                  + pa_y[k] * fg_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_x, pb_z, ff_8, ff_9, gd0_5, gd0_6, gd1_5, \
                         gd1_6, gf_8, gf_9, gf_10, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pb_z[k] * gf_8[k];

        t_15[k] = f_7 * ff_8[k]
                  + f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_x[k] * gf_10[k];

        t_16[k] = f_3 * gd0_5[k]
                  - f_4 * gd1_5[k]
                  + pb_z[k] * gf_9[k];

        t_17[k] = f_7 * ff_9[k]
                  + pb_x[k] * gf_11[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pb_z, dg0_15, dg1_14, fg_11, gd0_6, \
                         gd0_7, gd1_6, gd1_7, gf_11, gf_12, gf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * dg0_15[k]
                  - f_6 * dg1_14[k]
                  + pa_x[k] * fg_11[k];

        t_19[k] = pb_z[k] * gf_11[k];

        t_20[k] = f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_z[k] * gf_12[k];

        t_21[k] = f_1 * gd0_7[k]
                  - f_2 * gd1_7[k]
                  + pb_z[k] * gf_13[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_z, pb_y, dg0_0, dg1_0, fg_10, gd0_8, gd1_8, \
                         gf_14, gf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * dg0_0[k]
                  - f_6 * dg1_0[k]
                  + pa_z[k] * fg_10[k];

        t_23[k] = pb_y[k] * gf_14[k];

        t_24[k] = f_3 * gd0_8[k]
                  - f_4 * gd1_8[k]
                  + pb_y[k] * gf_15[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pb_x, pb_y, ff_10, ff_11, gd0_9, gd0_10, \
                         gd1_9, gd1_10, gf_16, gf_17, gf_18, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * ff_10[k]
                  + f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_x[k] * gf_16[k];

        t_26[k] = f_7 * ff_11[k]
                  + pb_x[k] * gf_19[k];

        t_27[k] = f_1 * gd0_9[k]
                  - f_2 * gd1_9[k]
                  + pb_y[k] * gf_17[k];

        t_28[k] = f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_y[k] * gf_18[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pa_x, pb_y, dg0_29, dg1_26, ff_12, \
                         ff_21, fg_12, fg_13, fg_18, fg_24, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_y[k] * gf_19[k];

        t_30[k] = f_5 * dg0_29[k]
                  - f_6 * dg1_26[k]
                  + pa_x[k] * fg_12[k];

        t_31[k] = f_0 * ff_12[k]
                  + pa_x[k] * fg_13[k];

        t_32[k] = pa_x[k] * fg_18[k];

        t_33[k] = f_0 * ff_21[k]
                  + pa_x[k] * fg_24[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, fg_32, gd0_13, gd0_14, gd0_15, \
                         gd1_13, gd1_14, gd1_15, gf_22, gf_23, gf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_x[k] * fg_32[k];

        t_35[k] = f_1 * gd0_13[k]
                  - f_2 * gd1_13[k]
                  + pb_x[k] * gf_22[k];

        t_36[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_x[k] * gf_23[k];

        t_37[k] = f_3 * gd0_15[k]
                  - f_4 * gd1_15[k]
                  + pb_x[k] * gf_24[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pb_x, pb_y, pb_z, ff_15, gd0_14, \
                         gd1_14, gf_25, gf_26, gf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_x[k] * gf_25[k];

        t_39[k] = pb_x[k] * gf_27[k];

        t_40[k] = f_0 * ff_15[k]
                  + f_1 * gd0_14[k]
                  - f_2 * gd1_14[k]
                  + pb_y[k] * gf_25[k];

        t_41[k] = pb_z[k] * gf_25[k];

        t_42[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_z[k] * gf_26[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_z, pb_x, pb_z, ff_17, fg_18, fg_21, \
                         gd0_15, gd0_17, gd1_15, gd1_17, gf_27, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * gd0_15[k]
                  - f_2 * gd1_15[k]
                  + pb_z[k] * gf_27[k];

        t_44[k] = pa_z[k] * fg_18[k];

        t_45[k] = f_0 * ff_17[k]
                  + pa_z[k] * fg_21[k];

        t_46[k] = f_1 * gd0_17[k]
                  - f_2 * gd1_17[k]
                  + pb_x[k] * gf_29[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_x, gd0_18, gd0_19, gd1_18, gd1_19, gf_30, \
                         gf_31, gf_32, gf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * gd0_18[k]
                  - f_4 * gd1_18[k]
                  + pb_x[k] * gf_30[k];

        t_48[k] = f_3 * gd0_19[k]
                  - f_4 * gd1_19[k]
                  + pb_x[k] * gf_31[k];

        t_49[k] = pb_x[k] * gf_32[k];

        t_50[k] = pb_x[k] * gf_34[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pa_z, pb_y, dg0_15, dg1_14, ff_19, ff_20, fg_22, \
                         gd0_19, gd1_19, gf_33, gf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_5 * dg0_15[k]
                  - f_6 * dg1_14[k]
                  + pa_z[k] * fg_22[k];

        t_52[k] = f_7 * ff_19[k]
                  + f_3 * gd0_19[k]
                  - f_4 * gd1_19[k]
                  + pb_y[k] * gf_33[k];

        t_53[k] = f_7 * ff_20[k]
                  + pb_y[k] * gf_34[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pb_x, dg0_29, dg1_26, ff_24, fg_23, \
                         fg_29, fg_32, gd0_21, gd1_21, gf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_5 * dg0_29[k]
                  - f_6 * dg1_26[k]
                  + pa_y[k] * fg_23[k];

        t_55[k] = f_0 * ff_24[k]
                  + pa_y[k] * fg_29[k];

        t_56[k] = pa_y[k] * fg_32[k];

        t_57[k] = f_1 * gd0_21[k]
                  - f_2 * gd1_21[k]
                  + pb_x[k] * gf_36[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_x, pb_y, gd0_22, gd0_23, gd1_22, \
                         gd1_23, gf_37, gf_38, gf_39, gf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_3 * gd0_22[k]
                  - f_4 * gd1_22[k]
                  + pb_x[k] * gf_37[k];

        t_59[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pb_x[k] * gf_38[k];

        t_60[k] = pb_x[k] * gf_39[k];

        t_61[k] = pb_x[k] * gf_41[k];

        t_62[k] = f_1 * gd0_22[k]
                  - f_2 * gd1_22[k]
                  + pb_y[k] * gf_39[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_y, pb_z, ff_26, gd0_23, gd1_23, gf_40, \
                         gf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pb_y[k] * gf_40[k];

        t_64[k] = pb_y[k] * gf_41[k];

        t_65[k] = f_0 * ff_26[k]
                  + f_1 * gd0_23[k]
                  - f_2 * gd1_23[k]
                  + pb_z[k] * gf_41[k];
    }
}

auto
compute_prim_gg_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg0, const size_t dg1,
                                      const size_t ff, const size_t fg, const size_t gd0,
                                      const size_t gd1, const size_t gf, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_1 = buffer.data(dg1 + 1);
    const auto *dg1_2 = buffer.data(dg1 + 2);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_11 = buffer.data(ff + 11);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_17 = buffer.data(gf + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dg0_0, dg1_0, ff_0, fg_0, fg_1, \
                         gd0_0, gd1_0, gf_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pa_y[k] * fg_0[k];

        t_2[k] = pa_z[k] * fg_0[k];

        t_3[k] = f_3 * dg0_0[k]
                 - f_4 * dg1_0[k]
                 + pa_y[k] * fg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_z, pb_x, dg0_0, dg0_1, dg1_0, dg1_1, ff_3, \
                         fg_2, fg_3, gd0_4, gd1_4, gf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * ff_3[k]
                 + f_6 * gd0_4[k]
                 - f_7 * gd1_4[k]
                 + pb_x[k] * gf_4[k];

        t_5[k] = f_3 * dg0_1[k]
                 - f_4 * dg1_1[k]
                 + pa_x[k] * fg_3[k];

        t_6[k] = f_3 * dg0_0[k]
                 - f_4 * dg1_0[k]
                 + pa_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pa_x, pb_x, dg0_2, dg1_2, ff_5, fg_4, fg_5, \
                         fg_8, gd0_6, gd1_6, gf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * ff_5[k]
                 + f_6 * gd0_6[k]
                 - f_7 * gd1_6[k]
                 + pb_x[k] * gf_7[k];

        t_8[k] = f_3 * dg0_2[k]
                 - f_4 * dg1_2[k]
                 + pa_x[k] * fg_4[k];

        t_9[k] = pa_x[k] * fg_5[k];

        t_10[k] = pa_x[k] * fg_8[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_z, pb_y, dg0_1, dg1_1, ff_7, fg_5, fg_6, gd0_9, \
                         gd1_9, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * ff_7[k]
                  + f_1 * gd0_9[k]
                  - f_2 * gd1_9[k]
                  + pb_y[k] * gf_11[k];

        t_12[k] = pa_z[k] * fg_5[k];

        t_13[k] = f_3 * dg0_1[k]
                  - f_4 * dg1_1[k]
                  + pa_z[k] * fg_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_y, pb_y, dg0_2, dg1_2, ff_9, fg_7, fg_8, gd0_12, \
                         gd1_12, gf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * ff_9[k]
                  + f_6 * gd0_12[k]
                  - f_7 * gd1_12[k]
                  + pb_y[k] * gf_14[k];

        t_15[k] = f_3 * dg0_2[k]
                  - f_4 * dg1_2[k]
                  + pa_y[k] * fg_7[k];

        t_16[k] = pa_y[k] * fg_8[k];
    }

#pragma omp simd aligned(t_17, pb_z, ff_11, gd0_14, gd1_14, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * ff_11[k]
                  + f_1 * gd0_14[k]
                  - f_2 * gd1_14[k]
                  + pb_z[k] * gf_17[k];
    }
}

auto
compute_prim_gg_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg0, const size_t dg1,
                                      const size_t ff, const size_t fg, const size_t gd0,
                                      const size_t gd1, const size_t gf, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_1 = buffer.data(dg1 + 1);
    const auto *dg1_2 = buffer.data(dg1 + 2);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_11 = buffer.data(ff + 11);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_23 = buffer.data(gd1 + 23);
    const auto *gd1_30 = buffer.data(gd1 + 30);
    const auto *gd1_35 = buffer.data(gd1 + 35);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_54 = buffer.data(gf + 54);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dg0_0, dg1_0, ff_0, fg_0, fg_1, \
                         gd0_0, gd1_0, gf_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pa_y[k] * fg_0[k];

        t_2[k] = pa_z[k] * fg_0[k];

        t_3[k] = f_3 * dg0_0[k]
                 - f_4 * dg1_0[k]
                 + pa_y[k] * fg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_z, pb_x, dg0_0, dg0_1, dg1_0, dg1_1, ff_3, \
                         fg_2, fg_3, gd0_4, gd1_10, gf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * ff_3[k]
                 + f_6 * gd0_4[k]
                 - f_7 * gd1_10[k]
                 + pb_x[k] * gf_12[k];

        t_5[k] = f_3 * dg0_1[k]
                 - f_4 * dg1_1[k]
                 + pa_x[k] * fg_3[k];

        t_6[k] = f_3 * dg0_0[k]
                 - f_4 * dg1_0[k]
                 + pa_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pa_x, pb_x, dg0_2, dg1_2, ff_5, fg_4, fg_5, \
                         fg_8, gd0_6, gd1_16, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_5 * ff_5[k]
                 + f_6 * gd0_6[k]
                 - f_7 * gd1_16[k]
                 + pb_x[k] * gf_17[k];

        t_8[k] = f_3 * dg0_2[k]
                 - f_4 * dg1_2[k]
                 + pa_x[k] * fg_4[k];

        t_9[k] = pa_x[k] * fg_5[k];

        t_10[k] = pa_x[k] * fg_8[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_z, pb_y, dg0_1, dg1_1, ff_7, fg_5, fg_6, gd0_9, \
                         gd1_23, gf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * ff_7[k]
                  + f_1 * gd0_9[k]
                  - f_2 * gd1_23[k]
                  + pb_y[k] * gf_30[k];

        t_12[k] = pa_z[k] * fg_5[k];

        t_13[k] = f_3 * dg0_1[k]
                  - f_4 * dg1_1[k]
                  + pa_z[k] * fg_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_y, pb_y, dg0_2, dg1_2, ff_9, fg_7, fg_8, gd0_12, \
                         gd1_30, gf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * ff_9[k]
                  + f_6 * gd0_12[k]
                  - f_7 * gd1_30[k]
                  + pb_y[k] * gf_41[k];

        t_15[k] = f_3 * dg0_2[k]
                  - f_4 * dg1_2[k]
                  + pa_y[k] * fg_7[k];

        t_16[k] = pa_y[k] * fg_8[k];
    }

#pragma omp simd aligned(t_17, pb_z, ff_11, gd0_14, gd1_35, gf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * ff_11[k]
                  + f_1 * gd0_14[k]
                  - f_2 * gd1_35[k]
                  + pb_z[k] * gf_54[k];
    }
}

auto
compute_prim_gg_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg0, const size_t dg1,
                                      const size_t ff, const size_t fg, const size_t gd0,
                                      const size_t gd1, const size_t gf, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_1 = buffer.data(dg1 + 1);
    const auto *dg1_2 = buffer.data(dg1 + 2);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_17 = buffer.data(ff + 17);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_22 = buffer.data(gd0 + 22);
    const auto *gd0_23 = buffer.data(gd0 + 23);
    const auto *gd0_24 = buffer.data(gd0 + 24);
    const auto *gd0_30 = buffer.data(gd0 + 30);
    const auto *gd0_33 = buffer.data(gd0 + 33);
    const auto *gd0_34 = buffer.data(gd0 + 34);
    const auto *gd0_35 = buffer.data(gd0 + 35);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_24 = buffer.data(gd1 + 24);
    const auto *gd1_27 = buffer.data(gd1 + 27);
    const auto *gd1_28 = buffer.data(gd1 + 28);
    const auto *gd1_29 = buffer.data(gd1 + 29);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, ff_0, gd0_0, gd0_1, gd1_0, \
                         gd1_1, gf_0, gf_1, gf_2, gf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_2[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];

        t_3[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, pb_y, pb_z, fg_0, gd0_2, gd1_2, gf_4, \
                         gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_4[k];

        t_5[k] = f_1 * gd0_2[k]
                 - f_2 * gd1_2[k]
                 + pb_z[k] * gf_5[k];

        t_6[k] = pa_y[k] * fg_0[k];

        t_7[k] = pa_z[k] * fg_0[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_x, pa_y, pb_x, dg0_0, dg0_1, dg1_0, dg1_1, ff_5, \
                         fg_1, fg_3, gd0_10, gd1_8, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * dg0_0[k]
                 - f_6 * dg1_0[k]
                 + pa_y[k] * fg_1[k];

        t_9[k] = f_7 * ff_5[k]
                 + f_3 * gd0_10[k]
                 - f_4 * gd1_8[k]
                 + pb_x[k] * gf_9[k];

        t_10[k] = f_5 * dg0_1[k]
                  - f_6 * dg1_1[k]
                  + pa_x[k] * fg_3[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_z, pb_x, dg0_0, dg0_2, dg1_0, dg1_2, ff_7, \
                         fg_2, fg_4, gd0_16, gd1_12, gf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_5 * dg0_0[k]
                  - f_6 * dg1_0[k]
                  + pa_z[k] * fg_2[k];

        t_12[k] = f_7 * ff_7[k]
                  + f_3 * gd0_16[k]
                  - f_4 * gd1_12[k]
                  + pb_x[k] * gf_12[k];

        t_13[k] = f_5 * dg0_2[k]
                  - f_6 * dg1_2[k]
                  + pa_x[k] * fg_4[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pb_x, fg_5, fg_8, gd0_22, gd0_23, \
                         gd1_17, gd1_18, gf_16, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_x[k] * fg_5[k];

        t_15[k] = pa_x[k] * fg_8[k];

        t_16[k] = f_1 * gd0_22[k]
                  - f_2 * gd1_17[k]
                  + pb_x[k] * gf_16[k];

        t_17[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_18[k]
                  + pb_x[k] * gf_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_x, pb_y, pb_z, ff_10, gd0_23, gd0_24, \
                         gd1_18, gd1_19, gf_18, gf_19, gf_20, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * gd0_24[k]
                  - f_4 * gd1_19[k]
                  + pb_x[k] * gf_18[k];

        t_19[k] = f_0 * ff_10[k]
                  + f_1 * gd0_23[k]
                  - f_2 * gd1_18[k]
                  + pb_y[k] * gf_19[k];

        t_20[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_18[k]
                  + pb_z[k] * gf_20[k];

        t_21[k] = f_1 * gd0_24[k]
                  - f_2 * gd1_19[k]
                  + pb_z[k] * gf_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_z, pb_y, dg0_1, dg1_1, ff_13, fg_5, fg_6, \
                         gd0_30, gd1_24, gf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pa_z[k] * fg_5[k];

        t_23[k] = f_5 * dg0_1[k]
                  - f_6 * dg1_1[k]
                  + pa_z[k] * fg_6[k];

        t_24[k] = f_7 * ff_13[k]
                  + f_3 * gd0_30[k]
                  - f_4 * gd1_24[k]
                  + pb_y[k] * gf_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, pb_x, dg0_2, dg1_2, fg_7, fg_8, gd0_33, \
                         gd0_34, gd1_27, gd1_28, gf_27, gf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_5 * dg0_2[k]
                  - f_6 * dg1_2[k]
                  + pa_y[k] * fg_7[k];

        t_26[k] = pa_y[k] * fg_8[k];

        t_27[k] = f_1 * gd0_33[k]
                  - f_2 * gd1_27[k]
                  + pb_x[k] * gf_27[k];

        t_28[k] = f_3 * gd0_34[k]
                  - f_4 * gd1_28[k]
                  + pb_x[k] * gf_28[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pb_x, pb_y, pb_z, ff_17, gd0_34, gd0_35, \
                         gd1_28, gd1_29, gf_29, gf_30, gf_31, gf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_3 * gd0_35[k]
                  - f_4 * gd1_29[k]
                  + pb_x[k] * gf_29[k];

        t_30[k] = f_1 * gd0_34[k]
                  - f_2 * gd1_28[k]
                  + pb_y[k] * gf_30[k];

        t_31[k] = f_3 * gd0_35[k]
                  - f_4 * gd1_29[k]
                  + pb_y[k] * gf_31[k];

        t_32[k] = f_0 * ff_17[k]
                  + f_1 * gd0_35[k]
                  - f_2 * gd1_29[k]
                  + pb_z[k] * gf_32[k];
    }
}

auto
compute_prim_gg_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg0, const size_t dg1,
                                      const size_t ff, const size_t fg, const size_t gd0,
                                      const size_t gd1, const size_t gf, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_2 = buffer.data(dg1 + 2);
    const auto *dg1_5 = buffer.data(dg1 + 5);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_62 = buffer.data(gf + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dg0_0, dg1_0, ff_0, fg_0, fg_1, \
                         gd0_0, gd1_0, gf_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pa_y[k] * fg_0[k];

        t_2[k] = pa_z[k] * fg_0[k];

        t_3[k] = f_3 * dg0_0[k]
                 - f_4 * dg1_0[k]
                 + pa_y[k] * fg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pb_x, dg0_1, dg1_2, ff_3, ff_4, fg_5, gd0_4, \
                         gd1_6, gf_16, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * ff_3[k]
                 + f_6 * gd0_4[k]
                 - f_7 * gd1_6[k]
                 + pb_x[k] * gf_16[k];

        t_5[k] = f_5 * ff_4[k]
                 + pb_x[k] * gf_17[k];

        t_6[k] = f_3 * dg0_1[k]
                 - f_4 * dg1_2[k]
                 + pa_x[k] * fg_5[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_x, dg0_0, dg1_0, ff_5, ff_6, fg_2, gd0_6, \
                         gd1_10, gf_24, gf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * dg0_0[k]
                 - f_4 * dg1_0[k]
                 + pa_z[k] * fg_2[k];

        t_8[k] = f_5 * ff_5[k]
                 + f_6 * gd0_6[k]
                 - f_7 * gd1_10[k]
                 + pb_x[k] * gf_24[k];

        t_9[k] = f_5 * ff_6[k]
                 + pb_x[k] * gf_27[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_y, dg0_2, dg1_5, ff_7, fg_8, fg_9, \
                         fg_14, gd0_9, gd1_14, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dg0_2[k]
                  - f_4 * dg1_5[k]
                  + pa_x[k] * fg_8[k];

        t_11[k] = pa_x[k] * fg_9[k];

        t_12[k] = pa_x[k] * fg_14[k];

        t_13[k] = f_0 * ff_7[k]
                  + f_1 * gd0_9[k]
                  - f_2 * gd1_14[k]
                  + pb_y[k] * gf_39[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_y, dg0_1, dg1_2, ff_9, ff_10, fg_9, \
                         fg_10, gd0_12, gd1_19, gf_50, gf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * fg_9[k];

        t_15[k] = f_3 * dg0_1[k]
                  - f_4 * dg1_2[k]
                  + pa_z[k] * fg_10[k];

        t_16[k] = f_5 * ff_9[k]
                  + f_6 * gd0_12[k]
                  - f_7 * gd1_19[k]
                  + pb_y[k] * gf_50[k];

        t_17[k] = f_5 * ff_10[k]
                  + pb_y[k] * gf_51[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pb_z, dg0_2, dg1_5, ff_11, fg_13, fg_14, \
                         gd0_14, gd1_23, gf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * dg0_2[k]
                  - f_4 * dg1_5[k]
                  + pa_y[k] * fg_13[k];

        t_19[k] = pa_y[k] * fg_14[k];

        t_20[k] = f_0 * ff_11[k]
                  + f_1 * gd0_14[k]
                  - f_2 * gd1_23[k]
                  + pb_z[k] * gf_62[k];
    }
}

auto
compute_prim_gg_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg0, const size_t dg1,
                                      const size_t ff, const size_t fg, const size_t gd0,
                                      const size_t gd1, const size_t gf, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / p;
    const auto f_7 = 0.5 / alpha;
    const auto f_8 = 0.5 * beta / (alpha * p);
    const auto f_9 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_2 = buffer.data(dg0 + 2);
    const auto *dg0_5 = buffer.data(dg0 + 5);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_2 = buffer.data(dg1 + 2);
    const auto *dg1_5 = buffer.data(dg1 + 5);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_17 = buffer.data(gd0 + 17);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_19 = buffer.data(gd0 + 19);
    const auto *gd0_21 = buffer.data(gd0 + 21);
    const auto *gd0_22 = buffer.data(gd0 + 22);
    const auto *gd0_23 = buffer.data(gd0 + 23);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_20 = buffer.data(gd1 + 20);
    const auto *gd1_21 = buffer.data(gd1 + 21);
    const auto *gd1_24 = buffer.data(gd1 + 24);
    const auto *gd1_25 = buffer.data(gd1 + 25);
    const auto *gd1_26 = buffer.data(gd1 + 26);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ff_0, gd0_0, gd1_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, fg_0, gd0_1, gd0_2, gd1_1, \
                         gd1_2, gf_3, gf_5, gf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];

        t_6[k] = pb_z[k] * gf_3[k];

        t_7[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_5[k];

        t_8[k] = f_1 * gd0_2[k]
                 - f_2 * gd1_2[k]
                 + pb_z[k] * gf_6[k];

        t_9[k] = pa_y[k] * fg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_z, ff_0, ff_1, ff_3, fg_0, \
                         fg_1, fg_3, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * ff_1[k]
                  + pa_y[k] * fg_1[k];

        t_11[k] = f_0 * ff_3[k]
                  + pa_y[k] * fg_3[k];

        t_12[k] = pa_z[k] * fg_0[k];

        t_13[k] = f_6 * ff_0[k]
                  + pb_z[k] * gf_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pa_z, dg0_0, dg1_0, ff_2, ff_4, ff_5, \
                         fg_2, fg_4, fg_5, fg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * ff_2[k]
                  + pa_z[k] * fg_2[k];

        t_15[k] = f_5 * ff_4[k]
                  + pa_z[k] * fg_4[k];

        t_16[k] = f_0 * ff_5[k]
                  + pa_z[k] * fg_5[k];

        t_17[k] = f_7 * dg0_0[k]
                  - f_8 * dg1_0[k]
                  + pa_y[k] * fg_6[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_x, pb_z, ff_11, ff_12, gd0_5, gd0_6, gd1_6, \
                         gd1_7, gf_12, gf_13, gf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * ff_11[k]
                  + f_3 * gd0_6[k]
                  - f_4 * gd1_7[k]
                  + pb_x[k] * gf_13[k];

        t_19[k] = f_3 * gd0_5[k]
                  - f_4 * gd1_6[k]
                  + pb_z[k] * gf_12[k];

        t_20[k] = f_5 * ff_12[k]
                  + pb_x[k] * gf_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pb_z, dg0_2, dg1_2, fg_10, gd0_6, gd0_7, \
                         gd1_7, gd1_8, gf_15, gf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * dg0_2[k]
                  - f_8 * dg1_2[k]
                  + pa_x[k] * fg_10[k];

        t_22[k] = f_3 * gd0_6[k]
                  - f_4 * gd1_7[k]
                  + pb_z[k] * gf_15[k];

        t_23[k] = f_1 * gd0_7[k]
                  - f_2 * gd1_8[k]
                  + pb_z[k] * gf_16[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_z, pb_y, pb_z, dg0_0, dg1_0, ff_8, fg_7, gd0_8, \
                         gd1_9, gf_17, gf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_7 * dg0_0[k]
                  - f_8 * dg1_0[k]
                  + pa_z[k] * fg_7[k];

        t_25[k] = f_5 * ff_8[k]
                  + pb_z[k] * gf_17[k];

        t_26[k] = f_3 * gd0_8[k]
                  - f_4 * gd1_9[k]
                  + pb_y[k] * gf_18[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_x, pb_y, ff_14, ff_15, gd0_9, gd0_10, \
                         gd1_10, gd1_11, gf_19, gf_20, gf_21, gf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * ff_14[k]
                  + f_3 * gd0_10[k]
                  - f_4 * gd1_11[k]
                  + pb_x[k] * gf_19[k];

        t_28[k] = f_5 * ff_15[k]
                  + pb_x[k] * gf_22[k];

        t_29[k] = f_1 * gd0_9[k]
                  - f_2 * gd1_10[k]
                  + pb_y[k] * gf_20[k];

        t_30[k] = f_3 * gd0_10[k]
                  - f_4 * gd1_11[k]
                  + pb_y[k] * gf_21[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pb_x, dg0_5, dg1_5, ff_16, ff_18, \
                         ff_19, fg_13, fg_14, fg_15, gf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_7 * dg0_5[k]
                  - f_8 * dg1_5[k]
                  + pa_x[k] * fg_13[k];

        t_32[k] = f_0 * ff_16[k]
                  + pa_x[k] * fg_14[k];

        t_33[k] = f_5 * ff_18[k]
                  + pa_x[k] * fg_15[k];

        t_34[k] = f_6 * ff_19[k]
                  + pb_x[k] * gf_25[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_z, ff_13, ff_27, ff_29, fg_17, \
                         fg_24, fg_26, gf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_x[k] * fg_17[k];

        t_36[k] = f_0 * ff_27[k]
                  + pa_x[k] * fg_24[k];

        t_37[k] = f_9 * ff_13[k]
                  + pb_z[k] * gf_26[k];

        t_38[k] = f_5 * ff_29[k]
                  + pa_x[k] * fg_26[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pb_x, pb_z, ff_32, fg_29, gd0_13, \
                         gd1_14, gf_28, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_6 * ff_32[k]
                  + pb_x[k] * gf_28[k];

        t_40[k] = pa_x[k] * fg_29[k];

        t_41[k] = f_1 * gd0_13[k]
                  - f_2 * gd1_14[k]
                  + pb_x[k] * gf_29[k];

        t_42[k] = pb_z[k] * gf_29[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pb_x, pb_y, pb_z, ff_19, gd0_14, \
                         gd0_15, gd1_15, gd1_16, gf_31, gf_32, gf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_15[k]
                  + pb_x[k] * gf_31[k];

        t_44[k] = f_3 * gd0_15[k]
                  - f_4 * gd1_16[k]
                  + pb_x[k] * gf_32[k];

        t_45[k] = pb_x[k] * gf_33[k];

        t_46[k] = f_0 * ff_19[k]
                  + f_1 * gd0_14[k]
                  - f_2 * gd1_15[k]
                  + pb_y[k] * gf_33[k];

        t_47[k] = pb_z[k] * gf_33[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_y, pb_z, ff_17, ff_21, fg_16, \
                         gd0_14, gd0_15, gd1_15, gd1_16, gf_34, gf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_15[k]
                  + pb_z[k] * gf_34[k];

        t_49[k] = f_0 * ff_21[k]
                  + pb_y[k] * gf_35[k];

        t_50[k] = f_1 * gd0_15[k]
                  - f_2 * gd1_16[k]
                  + pb_z[k] * gf_35[k];

        t_51[k] = f_5 * ff_17[k]
                  + pa_z[k] * fg_16[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_z, pb_y, pb_z, ff_19, ff_20, ff_23, fg_17, \
                         fg_18, gf_36, gf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pa_z[k] * fg_17[k];

        t_53[k] = f_6 * ff_19[k]
                  + pb_z[k] * gf_36[k];

        t_54[k] = f_5 * ff_20[k]
                  + pa_z[k] * fg_18[k];

        t_55[k] = f_9 * ff_23[k]
                  + pb_y[k] * gf_37[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pa_z, pb_x, ff_21, fg_19, gd0_17, gd0_18, gd1_19, \
                         gd1_20, gf_38, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * ff_21[k]
                  + pa_z[k] * fg_19[k];

        t_57[k] = f_1 * gd0_17[k]
                  - f_2 * gd1_19[k]
                  + pb_x[k] * gf_38[k];

        t_58[k] = f_3 * gd0_18[k]
                  - f_4 * gd1_20[k]
                  + pb_x[k] * gf_39[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pa_z, pb_x, pb_z, dg0_2, dg1_2, ff_22, fg_20, \
                         gd0_19, gd1_21, gf_40, gf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_3 * gd0_19[k]
                  - f_4 * gd1_21[k]
                  + pb_x[k] * gf_40[k];

        t_60[k] = f_7 * dg0_2[k]
                  - f_8 * dg1_2[k]
                  + pa_z[k] * fg_20[k];

        t_61[k] = f_5 * ff_22[k]
                  + pb_z[k] * gf_41[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pa_y, pb_y, dg0_5, dg1_5, ff_25, ff_26, fg_23, \
                         gd0_19, gd1_21, gf_42, gf_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_5 * ff_25[k]
                  + f_3 * gd0_19[k]
                  - f_4 * gd1_21[k]
                  + pb_y[k] * gf_42[k];

        t_63[k] = f_5 * ff_26[k]
                  + pb_y[k] * gf_43[k];

        t_64[k] = f_7 * dg0_5[k]
                  - f_8 * dg1_5[k]
                  + pa_y[k] * fg_23[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_y, pb_z, ff_24, ff_28, ff_30, ff_31, \
                         fg_25, fg_27, fg_28, gf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_5 * ff_28[k]
                  + pa_y[k] * fg_25[k];

        t_66[k] = f_0 * ff_30[k]
                  + pa_y[k] * fg_27[k];

        t_67[k] = f_9 * ff_24[k]
                  + pb_z[k] * gf_44[k];

        t_68[k] = f_5 * ff_31[k]
                  + pa_y[k] * fg_28[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pa_y, pb_x, pb_y, pb_z, ff_27, ff_32, \
                         fg_29, gd0_21, gd1_24, gf_46, gf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_6 * ff_32[k]
                  + pb_y[k] * gf_46[k];

        t_70[k] = pa_y[k] * fg_29[k];

        t_71[k] = f_1 * gd0_21[k]
                  - f_2 * gd1_24[k]
                  + pb_x[k] * gf_47[k];

        t_72[k] = pb_y[k] * gf_47[k];

        t_73[k] = f_0 * ff_27[k]
                  + pb_z[k] * gf_47[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_x, pb_y, gd0_22, gd0_23, gd1_25, gd1_26, \
                         gf_49, gf_50, gf_51, gf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_3 * gd0_22[k]
                  - f_4 * gd1_25[k]
                  + pb_x[k] * gf_49[k];

        t_75[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_26[k]
                  + pb_x[k] * gf_50[k];

        t_76[k] = pb_x[k] * gf_53[k];

        t_77[k] = f_1 * gd0_22[k]
                  - f_2 * gd1_25[k]
                  + pb_y[k] * gf_51[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_y, pb_z, ff_30, ff_32, gd0_23, gd1_26, \
                         gf_51, gf_52, gf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * ff_30[k]
                  + pb_z[k] * gf_51[k];

        t_79[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_26[k]
                  + pb_y[k] * gf_52[k];

        t_80[k] = pb_y[k] * gf_53[k];

        t_81[k] = f_0 * ff_32[k]
                  + f_1 * gd0_23[k]
                  - f_2 * gd1_26[k]
                  + pb_z[k] * gf_53[k];
    }
}

auto
compute_prim_gg_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg0, const size_t dg1,
                                      const size_t ff, const size_t fg, const size_t gd0,
                                      const size_t gd1, const size_t gf, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / p;

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

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_2 = buffer.data(dg0 + 2);
    const auto *dg0_5 = buffer.data(dg0 + 5);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_2 = buffer.data(dg1 + 2);
    const auto *dg1_5 = buffer.data(dg1 + 5);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_26 = buffer.data(ff + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_21 = buffer.data(gd0 + 21);
    const auto *gd0_24 = buffer.data(gd0 + 24);
    const auto *gd0_25 = buffer.data(gd0 + 25);
    const auto *gd0_26 = buffer.data(gd0 + 26);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_21 = buffer.data(gd1 + 21);
    const auto *gd1_22 = buffer.data(gd1 + 22);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ff_0, gd0_0, gd1_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, fg_0, gd0_1, gd0_2, gd1_1, \
                         gd1_2, gf_3, gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];

        t_6[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_4[k];

        t_7[k] = pb_y[k] * gf_5[k];

        t_8[k] = f_1 * gd0_2[k]
                 - f_2 * gd1_2[k]
                 + pb_z[k] * gf_5[k];

        t_9[k] = pa_y[k] * fg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pa_z, pb_x, dg0_0, dg1_0, ff_8, fg_0, fg_1, \
                         gd0_7, gd1_6, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * fg_0[k];

        t_11[k] = f_5 * dg0_0[k]
                  - f_6 * dg1_0[k]
                  + pa_y[k] * fg_1[k];

        t_12[k] = f_7 * ff_8[k]
                  + f_3 * gd0_7[k]
                  - f_4 * gd1_6[k]
                  + pb_x[k] * gf_9[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pa_z, pb_x, dg0_0, dg0_2, dg1_0, dg1_2, ff_9, \
                         fg_2, fg_3, gf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_7 * ff_9[k]
                  + pb_x[k] * gf_10[k];

        t_14[k] = f_5 * dg0_2[k]
                  - f_6 * dg1_2[k]
                  + pa_x[k] * fg_3[k];

        t_15[k] = f_5 * dg0_0[k]
                  - f_6 * dg1_0[k]
                  + pa_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_x, dg0_5, dg1_5, ff_10, ff_11, fg_4, \
                         fg_5, gd0_11, gd1_10, gf_12, gf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * ff_10[k]
                  + f_3 * gd0_11[k]
                  - f_4 * gd1_10[k]
                  + pb_x[k] * gf_12[k];

        t_17[k] = f_7 * ff_11[k]
                  + pb_x[k] * gf_13[k];

        t_18[k] = f_5 * dg0_5[k]
                  - f_6 * dg1_5[k]
                  + pa_x[k] * fg_4[k];

        t_19[k] = pa_x[k] * fg_5[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, fg_8, gd0_14, gd0_15, gd0_16, \
                         gd1_13, gd1_14, gd1_15, gf_16, gf_17, gf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_x[k] * fg_8[k];

        t_21[k] = f_1 * gd0_14[k]
                  - f_2 * gd1_13[k]
                  + pb_x[k] * gf_16[k];

        t_22[k] = f_3 * gd0_15[k]
                  - f_4 * gd1_14[k]
                  + pb_x[k] * gf_17[k];

        t_23[k] = f_3 * gd0_16[k]
                  - f_4 * gd1_15[k]
                  + pb_x[k] * gf_18[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pb_x, pb_y, pb_z, ff_15, gd0_15, \
                         gd1_14, gf_19, gf_20, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_x[k] * gf_19[k];

        t_25[k] = pb_x[k] * gf_21[k];

        t_26[k] = f_0 * ff_15[k]
                  + f_1 * gd0_15[k]
                  - f_2 * gd1_14[k]
                  + pb_y[k] * gf_19[k];

        t_27[k] = pb_z[k] * gf_19[k];

        t_28[k] = f_3 * gd0_15[k]
                  - f_4 * gd1_14[k]
                  + pb_z[k] * gf_20[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_z, pb_z, dg0_2, dg1_2, fg_5, fg_6, gd0_16, \
                         gd1_15, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * gd0_16[k]
                  - f_2 * gd1_15[k]
                  + pb_z[k] * gf_21[k];

        t_30[k] = pa_z[k] * fg_5[k];

        t_31[k] = f_5 * dg0_2[k]
                  - f_6 * dg1_2[k]
                  + pa_z[k] * fg_6[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_y, pb_y, dg0_5, dg1_5, ff_19, ff_20, fg_7, \
                         fg_8, gd0_21, gd1_19, gf_24, gf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_7 * ff_19[k]
                  + f_3 * gd0_21[k]
                  - f_4 * gd1_19[k]
                  + pb_y[k] * gf_24[k];

        t_33[k] = f_7 * ff_20[k]
                  + pb_y[k] * gf_25[k];

        t_34[k] = f_5 * dg0_5[k]
                  - f_6 * dg1_5[k]
                  + pa_y[k] * fg_7[k];

        t_35[k] = pa_y[k] * fg_8[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_x, gd0_24, gd0_25, gd0_26, gd1_21, gd1_22, \
                         gd1_23, gf_27, gf_28, gf_29, gf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * gd0_24[k]
                  - f_2 * gd1_21[k]
                  + pb_x[k] * gf_27[k];

        t_37[k] = f_3 * gd0_25[k]
                  - f_4 * gd1_22[k]
                  + pb_x[k] * gf_28[k];

        t_38[k] = f_3 * gd0_26[k]
                  - f_4 * gd1_23[k]
                  + pb_x[k] * gf_29[k];

        t_39[k] = pb_x[k] * gf_30[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pb_x, pb_y, pb_z, ff_26, gd0_25, \
                         gd0_26, gd1_22, gd1_23, gf_30, gf_31, gf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_x[k] * gf_32[k];

        t_41[k] = f_1 * gd0_25[k]
                  - f_2 * gd1_22[k]
                  + pb_y[k] * gf_30[k];

        t_42[k] = f_3 * gd0_26[k]
                  - f_4 * gd1_23[k]
                  + pb_y[k] * gf_31[k];

        t_43[k] = pb_y[k] * gf_32[k];

        t_44[k] = f_0 * ff_26[k]
                  + f_1 * gd0_26[k]
                  - f_2 * gd1_23[k]
                  + pb_z[k] * gf_32[k];
    }
}

auto
compute_prim_gg_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg0, const size_t dg1,
                                      const size_t ff, const size_t fg, const size_t gd0,
                                      const size_t gd1, const size_t gf, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_1 = buffer.data(dg1 + 1);
    const auto *dg1_2 = buffer.data(dg1 + 2);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_17 = buffer.data(gf + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dg0_0, dg1_0, ff_0, fg_0, fg_1, \
                         gd0_0, gd1_0, gf_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pa_y[k] * fg_0[k];

        t_2[k] = pa_z[k] * fg_0[k];

        t_3[k] = f_3 * dg0_0[k]
                 - f_4 * dg1_0[k]
                 + pa_y[k] * fg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pb_x, dg0_1, dg1_1, ff_3, ff_4, fg_3, gd0_4, \
                         gd1_4, gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * ff_3[k]
                 + f_6 * gd0_4[k]
                 - f_7 * gd1_4[k]
                 + pb_x[k] * gf_4[k];

        t_5[k] = f_5 * ff_4[k]
                 + pb_x[k] * gf_5[k];

        t_6[k] = f_3 * dg0_1[k]
                 - f_4 * dg1_1[k]
                 + pa_x[k] * fg_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_x, dg0_0, dg1_0, ff_5, ff_6, fg_2, gd0_6, \
                         gd1_6, gf_7, gf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * dg0_0[k]
                 - f_4 * dg1_0[k]
                 + pa_z[k] * fg_2[k];

        t_8[k] = f_5 * ff_5[k]
                 + f_6 * gd0_6[k]
                 - f_7 * gd1_6[k]
                 + pb_x[k] * gf_7[k];

        t_9[k] = f_5 * ff_6[k]
                 + pb_x[k] * gf_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_y, dg0_2, dg1_2, ff_7, fg_4, fg_5, \
                         fg_8, gd0_9, gd1_9, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dg0_2[k]
                  - f_4 * dg1_2[k]
                  + pa_x[k] * fg_4[k];

        t_11[k] = pa_x[k] * fg_5[k];

        t_12[k] = pa_x[k] * fg_8[k];

        t_13[k] = f_0 * ff_7[k]
                  + f_1 * gd0_9[k]
                  - f_2 * gd1_9[k]
                  + pb_y[k] * gf_11[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_y, dg0_1, dg1_1, ff_9, ff_10, fg_5, \
                         fg_6, gd0_12, gd1_12, gf_14, gf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * fg_5[k];

        t_15[k] = f_3 * dg0_1[k]
                  - f_4 * dg1_1[k]
                  + pa_z[k] * fg_6[k];

        t_16[k] = f_5 * ff_9[k]
                  + f_6 * gd0_12[k]
                  - f_7 * gd1_12[k]
                  + pb_y[k] * gf_14[k];

        t_17[k] = f_5 * ff_10[k]
                  + pb_y[k] * gf_15[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pb_z, dg0_2, dg1_2, ff_11, fg_7, fg_8, \
                         gd0_14, gd1_14, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * dg0_2[k]
                  - f_4 * dg1_2[k]
                  + pa_y[k] * fg_7[k];

        t_19[k] = pa_y[k] * fg_8[k];

        t_20[k] = f_0 * ff_11[k]
                  + f_1 * gd0_14[k]
                  - f_2 * gd1_14[k]
                  + pb_z[k] * gf_17[k];
    }
}

auto
compute_prim_gg_electron_repulsion_18(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg0, const size_t dg1,
                                      const size_t ff, const size_t fg, const size_t gd0,
                                      const size_t gd1, const size_t gf, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / beta;
    const auto f_7 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_3 = buffer.data(dg1 + 3);
    const auto *dg1_8 = buffer.data(dg1 + 8);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_44 = buffer.data(gf + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, dg0_0, dg1_0, ff_0, fg_0, fg_1, \
                         gd0_0, gd1_0, gf_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pa_y[k] * fg_0[k];

        t_2[k] = pa_z[k] * fg_0[k];

        t_3[k] = f_3 * dg0_0[k]
                 - f_4 * dg1_0[k]
                 + pa_y[k] * fg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pb_x, dg0_1, dg1_3, ff_3, ff_4, fg_5, gd0_4, \
                         gd1_6, gf_10, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * ff_3[k]
                 + f_6 * gd0_4[k]
                 - f_7 * gd1_6[k]
                 + pb_x[k] * gf_10[k];

        t_5[k] = f_5 * ff_4[k]
                 + pb_x[k] * gf_11[k];

        t_6[k] = f_3 * dg0_1[k]
                 - f_4 * dg1_3[k]
                 + pa_x[k] * fg_5[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_x, dg0_0, dg1_0, ff_5, ff_6, fg_2, gd0_6, \
                         gd1_10, gf_16, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * dg0_0[k]
                 - f_4 * dg1_0[k]
                 + pa_z[k] * fg_2[k];

        t_8[k] = f_5 * ff_5[k]
                 + f_6 * gd0_6[k]
                 - f_7 * gd1_10[k]
                 + pb_x[k] * gf_16[k];

        t_9[k] = f_5 * ff_6[k]
                 + pb_x[k] * gf_19[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_y, dg0_2, dg1_8, ff_7, fg_8, fg_9, \
                         fg_14, gd0_9, gd1_14, gf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dg0_2[k]
                  - f_4 * dg1_8[k]
                  + pa_x[k] * fg_8[k];

        t_11[k] = pa_x[k] * fg_9[k];

        t_12[k] = pa_x[k] * fg_14[k];

        t_13[k] = f_0 * ff_7[k]
                  + f_1 * gd0_9[k]
                  - f_2 * gd1_14[k]
                  + pb_y[k] * gf_27[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_y, dg0_1, dg1_3, ff_9, ff_10, fg_9, \
                         fg_10, gd0_12, gd1_19, gf_35, gf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * fg_9[k];

        t_15[k] = f_3 * dg0_1[k]
                  - f_4 * dg1_3[k]
                  + pa_z[k] * fg_10[k];

        t_16[k] = f_5 * ff_9[k]
                  + f_6 * gd0_12[k]
                  - f_7 * gd1_19[k]
                  + pb_y[k] * gf_35[k];

        t_17[k] = f_5 * ff_10[k]
                  + pb_y[k] * gf_36[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pb_z, dg0_2, dg1_8, ff_11, fg_13, fg_14, \
                         gd0_14, gd1_23, gf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * dg0_2[k]
                  - f_4 * dg1_8[k]
                  + pa_y[k] * fg_13[k];

        t_19[k] = pa_y[k] * fg_14[k];

        t_20[k] = f_0 * ff_11[k]
                  + f_1 * gd0_14[k]
                  - f_2 * gd1_23[k]
                  + pb_z[k] * gf_44[k];
    }
}

auto
compute_prim_gg_electron_repulsion_19(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg0, const size_t dg1,
                                      const size_t ff, const size_t fg, const size_t gd0,
                                      const size_t gd1, const size_t gf, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 0.5 / p;

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

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_3 = buffer.data(dg0 + 3);
    const auto *dg0_8 = buffer.data(dg0 + 8);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_3 = buffer.data(dg1 + 3);
    const auto *dg1_8 = buffer.data(dg1 + 8);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_44 = buffer.data(fg + 44);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_17 = buffer.data(gd0 + 17);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_19 = buffer.data(gd0 + 19);
    const auto *gd0_21 = buffer.data(gd0 + 21);
    const auto *gd0_22 = buffer.data(gd0 + 22);
    const auto *gd0_23 = buffer.data(gd0 + 23);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_21 = buffer.data(gd1 + 21);
    const auto *gd1_22 = buffer.data(gd1 + 22);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ff_0, gd0_0, gd1_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, gd0_1, gd0_2, gd1_1, gd1_2, \
                         gf_3, gf_5, gf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];

        t_6[k] = pb_z[k] * gf_3[k];

        t_7[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_5[k];

        t_8[k] = pb_y[k] * gf_6[k];

        t_9[k] = f_1 * gd0_2[k]
                 - f_2 * gd1_2[k]
                 + pb_z[k] * gf_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pa_y, pa_z, ff_3, ff_4, fg_0, \
                         fg_5, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * fg_0[k];

        t_11[k] = f_0 * ff_3[k]
                  + pa_y[k] * fg_5[k];

        t_12[k] = pa_y[k] * fg_8[k];

        t_13[k] = pa_z[k] * fg_0[k];

        t_14[k] = pa_z[k] * fg_5[k];

        t_15[k] = f_5 * ff_4[k]
                  + pa_z[k] * fg_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_y, pa_z, pb_y, pb_z, dg0_0, dg1_0, ff_6, \
                         fg_8, fg_9, gf_9, gf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * gf_9[k];

        t_17[k] = f_0 * ff_6[k]
                  + pa_z[k] * fg_8[k];

        t_18[k] = f_6 * dg0_0[k]
                  - f_7 * dg1_0[k]
                  + pa_y[k] * fg_9[k];

        t_19[k] = pb_z[k] * gf_10[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_x, pb_z, ff_9, ff_10, gd0_5, gd0_6, gd1_5, \
                         gd1_6, gf_11, gf_12, gf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * ff_9[k]
                  + f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_x[k] * gf_12[k];

        t_21[k] = f_3 * gd0_5[k]
                  - f_4 * gd1_5[k]
                  + pb_z[k] * gf_11[k];

        t_22[k] = f_5 * ff_10[k]
                  + pb_x[k] * gf_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_z, dg0_3, dg1_3, fg_16, gd0_6, \
                         gd0_7, gd1_6, gd1_7, gf_13, gf_14, gf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * dg0_3[k]
                  - f_7 * dg1_3[k]
                  + pa_x[k] * fg_16[k];

        t_24[k] = pb_z[k] * gf_13[k];

        t_25[k] = f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_z[k] * gf_14[k];

        t_26[k] = f_1 * gd0_7[k]
                  - f_2 * gd1_7[k]
                  + pb_z[k] * gf_15[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, dg0_0, dg1_0, fg_10, fg_11, \
                         fg_12, gf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_z[k] * fg_10[k];

        t_28[k] = pa_y[k] * fg_12[k];

        t_29[k] = f_6 * dg0_0[k]
                  - f_7 * dg1_0[k]
                  + pa_z[k] * fg_11[k];

        t_30[k] = pb_y[k] * gf_16[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pb_x, pb_y, ff_11, ff_12, gd0_8, gd0_10, gd1_8, \
                         gd1_10, gf_17, gf_18, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_3 * gd0_8[k]
                  - f_4 * gd1_8[k]
                  + pb_y[k] * gf_17[k];

        t_32[k] = f_5 * ff_11[k]
                  + f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_x[k] * gf_18[k];

        t_33[k] = f_5 * ff_12[k]
                  + pb_x[k] * gf_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_y, dg0_8, dg1_8, fg_20, gd0_9, \
                         gd0_10, gd1_9, gd1_10, gf_19, gf_20, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_1 * gd0_9[k]
                  - f_2 * gd1_9[k]
                  + pb_y[k] * gf_19[k];

        t_35[k] = f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_y[k] * gf_20[k];

        t_36[k] = pb_y[k] * gf_21[k];

        t_37[k] = f_6 * dg0_8[k]
                  - f_7 * dg1_8[k]
                  + pa_x[k] * fg_20[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_x, pa_z, pb_x, ff_13, ff_15, ff_17, \
                         fg_13, fg_21, fg_23, fg_26, gf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * ff_13[k]
                  + pa_x[k] * fg_21[k];

        t_39[k] = f_5 * ff_15[k]
                  + pa_x[k] * fg_23[k];

        t_40[k] = f_8 * ff_17[k]
                  + pb_x[k] * gf_23[k];

        t_41[k] = pa_x[k] * fg_26[k];

        t_42[k] = pa_z[k] * fg_13[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pb_x, ff_23, ff_26, ff_29, fg_36, \
                         fg_39, fg_44, gf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_0 * ff_23[k]
                  + pa_x[k] * fg_36[k];

        t_44[k] = f_5 * ff_26[k]
                  + pa_x[k] * fg_39[k];

        t_45[k] = f_8 * ff_29[k]
                  + pb_x[k] * gf_25[k];

        t_46[k] = pa_x[k] * fg_44[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_x, pb_z, gd0_13, gd0_14, gd0_15, gd1_13, \
                         gd1_14, gd1_15, gf_26, gf_28, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_1 * gd0_13[k]
                  - f_2 * gd1_13[k]
                  + pb_x[k] * gf_26[k];

        t_48[k] = pb_z[k] * gf_26[k];

        t_49[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_x[k] * gf_28[k];

        t_50[k] = f_3 * gd0_15[k]
                  - f_4 * gd1_15[k]
                  + pb_x[k] * gf_29[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pb_x, pb_y, pb_z, ff_17, gd0_14, \
                         gd1_14, gf_30, gf_31, gf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_x[k] * gf_30[k];

        t_52[k] = pb_x[k] * gf_32[k];

        t_53[k] = f_0 * ff_17[k]
                  + f_1 * gd0_14[k]
                  - f_2 * gd1_14[k]
                  + pb_y[k] * gf_30[k];

        t_54[k] = pb_z[k] * gf_30[k];

        t_55[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_z[k] * gf_31[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_z, pb_x, pb_z, ff_14, fg_21, fg_24, \
                         fg_26, gd0_15, gd1_15, gf_32, gf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * gd0_15[k]
                  - f_2 * gd1_15[k]
                  + pb_z[k] * gf_32[k];

        t_57[k] = pa_z[k] * fg_21[k];

        t_58[k] = f_5 * ff_14[k]
                  + pa_z[k] * fg_24[k];

        t_59[k] = pb_x[k] * gf_34[k];

        t_60[k] = pa_z[k] * fg_26[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pa_z, pb_x, ff_19, fg_29, gd0_17, gd0_18, gd1_17, \
                         gd1_18, gf_35, gf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_0 * ff_19[k]
                  + pa_z[k] * fg_29[k];

        t_62[k] = f_1 * gd0_17[k]
                  - f_2 * gd1_17[k]
                  + pb_x[k] * gf_35[k];

        t_63[k] = f_3 * gd0_18[k]
                  - f_4 * gd1_18[k]
                  + pb_x[k] * gf_36[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_z, pb_x, dg0_3, dg1_3, fg_30, gd0_19, \
                         gd1_19, gf_37, gf_38, gf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * gd0_19[k]
                  - f_4 * gd1_19[k]
                  + pb_x[k] * gf_37[k];

        t_65[k] = pb_x[k] * gf_38[k];

        t_66[k] = pb_x[k] * gf_40[k];

        t_67[k] = f_6 * dg0_3[k]
                  - f_7 * dg1_3[k]
                  + pa_z[k] * fg_30[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pa_y, pb_y, dg0_8, dg1_8, ff_21, ff_22, fg_35, \
                         gd0_19, gd1_19, gf_39, gf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_5 * ff_21[k]
                  + f_3 * gd0_19[k]
                  - f_4 * gd1_19[k]
                  + pb_y[k] * gf_39[k];

        t_69[k] = f_5 * ff_22[k]
                  + pb_y[k] * gf_40[k];

        t_70[k] = f_6 * dg0_8[k]
                  - f_7 * dg1_8[k]
                  + pa_y[k] * fg_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_x, ff_24, ff_27, ff_28, fg_38, \
                         fg_41, fg_42, gf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * ff_24[k]
                  + pa_y[k] * fg_38[k];

        t_72[k] = pb_x[k] * gf_41[k];

        t_73[k] = f_0 * ff_27[k]
                  + pa_y[k] * fg_41[k];

        t_74[k] = f_5 * ff_28[k]
                  + pa_y[k] * fg_42[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_y, pb_x, pb_y, ff_29, fg_44, gd0_21, \
                         gd1_21, gf_43, gf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_8 * ff_29[k]
                  + pb_y[k] * gf_43[k];

        t_76[k] = pa_y[k] * fg_44[k];

        t_77[k] = f_1 * gd0_21[k]
                  - f_2 * gd1_21[k]
                  + pb_x[k] * gf_44[k];

        t_78[k] = pb_y[k] * gf_44[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pb_x, pb_y, gd0_22, gd0_23, gd1_22, \
                         gd1_23, gf_46, gf_47, gf_48, gf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * gd0_22[k]
                  - f_4 * gd1_22[k]
                  + pb_x[k] * gf_46[k];

        t_80[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pb_x[k] * gf_47[k];

        t_81[k] = pb_x[k] * gf_48[k];

        t_82[k] = pb_x[k] * gf_50[k];

        t_83[k] = f_1 * gd0_22[k]
                  - f_2 * gd1_22[k]
                  + pb_y[k] * gf_48[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pb_y, pb_z, ff_29, gd0_23, gd1_23, gf_49, \
                         gf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pb_y[k] * gf_49[k];

        t_85[k] = pb_y[k] * gf_50[k];

        t_86[k] = f_0 * ff_29[k]
                  + f_1 * gd0_23[k]
                  - f_2 * gd1_23[k]
                  + pb_z[k] * gf_50[k];
    }
}

auto
compute_prim_gg_electron_repulsion_20(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg0, const size_t dg1,
                                      const size_t ff, const size_t fg, const size_t gd0,
                                      const size_t gd1, const size_t gf, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / p;
    const auto f_8 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_3 = buffer.data(dg0 + 3);
    const auto *dg0_8 = buffer.data(dg0 + 8);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_6 = buffer.data(dg1 + 6);
    const auto *dg1_14 = buffer.data(dg1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_38 = buffer.data(fg + 38);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_17 = buffer.data(gd0 + 17);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_19 = buffer.data(gd0 + 19);
    const auto *gd0_21 = buffer.data(gd0 + 21);
    const auto *gd0_22 = buffer.data(gd0 + 22);
    const auto *gd0_23 = buffer.data(gd0 + 23);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_21 = buffer.data(gd1 + 21);
    const auto *gd1_22 = buffer.data(gd1 + 22);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ff_0, gd0_0, gd1_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, gd0_1, gd0_2, gd1_1, gd1_2, \
                         gf_3, gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];

        t_6[k] = pb_z[k] * gf_3[k];

        t_7[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_4[k];

        t_8[k] = pb_y[k] * gf_5[k];

        t_9[k] = f_1 * gd0_2[k]
                 - f_2 * gd1_2[k]
                 + pb_z[k] * gf_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pa_z, dg0_0, dg1_0, ff_3, ff_5, \
                         fg_0, fg_5, fg_8, fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * fg_0[k];

        t_11[k] = f_0 * ff_3[k]
                  + pa_y[k] * fg_5[k];

        t_12[k] = pa_z[k] * fg_0[k];

        t_13[k] = f_0 * ff_5[k]
                  + pa_z[k] * fg_8[k];

        t_14[k] = f_5 * dg0_0[k]
                  - f_6 * dg1_0[k]
                  + pa_y[k] * fg_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_x, pb_z, ff_8, ff_9, gd0_5, gd0_6, gd1_5, \
                         gd1_6, gf_8, gf_9, gf_10, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_z[k] * gf_8[k];

        t_16[k] = f_7 * ff_8[k]
                  + f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_x[k] * gf_10[k];

        t_17[k] = f_3 * gd0_5[k]
                  - f_4 * gd1_5[k]
                  + pb_z[k] * gf_9[k];

        t_18[k] = f_7 * ff_9[k]
                  + pb_x[k] * gf_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_z, dg0_3, dg1_6, fg_13, gd0_6, \
                         gd0_7, gd1_6, gd1_7, gf_11, gf_12, gf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_5 * dg0_3[k]
                  - f_6 * dg1_6[k]
                  + pa_x[k] * fg_13[k];

        t_20[k] = pb_z[k] * gf_11[k];

        t_21[k] = f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_z[k] * gf_12[k];

        t_22[k] = f_1 * gd0_7[k]
                  - f_2 * gd1_7[k]
                  + pb_z[k] * gf_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_z, pb_y, dg0_0, dg1_0, fg_10, gd0_8, gd1_8, \
                         gf_14, gf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * dg0_0[k]
                  - f_6 * dg1_0[k]
                  + pa_z[k] * fg_10[k];

        t_24[k] = pb_y[k] * gf_14[k];

        t_25[k] = f_3 * gd0_8[k]
                  - f_4 * gd1_8[k]
                  + pb_y[k] * gf_15[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_x, pb_y, ff_10, ff_11, gd0_9, gd0_10, \
                         gd1_9, gd1_10, gf_16, gf_17, gf_18, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * ff_10[k]
                  + f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_x[k] * gf_16[k];

        t_27[k] = f_7 * ff_11[k]
                  + pb_x[k] * gf_19[k];

        t_28[k] = f_1 * gd0_9[k]
                  - f_2 * gd1_9[k]
                  + pb_y[k] * gf_17[k];

        t_29[k] = f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_y[k] * gf_18[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_y, dg0_8, dg1_14, ff_12, ff_13, \
                         fg_16, fg_17, fg_18, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * gf_19[k];

        t_31[k] = f_5 * dg0_8[k]
                  - f_6 * dg1_14[k]
                  + pa_x[k] * fg_16[k];

        t_32[k] = f_0 * ff_12[k]
                  + pa_x[k] * fg_17[k];

        t_33[k] = f_7 * ff_13[k]
                  + pa_x[k] * fg_18[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_x, pb_x, ff_15, ff_21, ff_23, ff_26, \
                         fg_22, fg_30, fg_32, gf_21, gf_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_8 * ff_15[k]
                  + pb_x[k] * gf_21[k];

        t_35[k] = pa_x[k] * fg_22[k];

        t_36[k] = f_0 * ff_21[k]
                  + pa_x[k] * fg_30[k];

        t_37[k] = f_7 * ff_23[k]
                  + pa_x[k] * fg_32[k];

        t_38[k] = f_8 * ff_26[k]
                  + pb_x[k] * gf_23[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pb_x, pb_z, fg_38, gd0_13, gd0_14, \
                         gd1_13, gd1_14, gf_24, gf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pa_x[k] * fg_38[k];

        t_40[k] = f_1 * gd0_13[k]
                  - f_2 * gd1_13[k]
                  + pb_x[k] * gf_24[k];

        t_41[k] = pb_z[k] * gf_24[k];

        t_42[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_x[k] * gf_25[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pb_x, pb_y, pb_z, ff_15, gd0_14, \
                         gd0_15, gd1_14, gd1_15, gf_26, gf_27, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * gd0_15[k]
                  - f_4 * gd1_15[k]
                  + pb_x[k] * gf_26[k];

        t_44[k] = pb_x[k] * gf_27[k];

        t_45[k] = pb_x[k] * gf_29[k];

        t_46[k] = f_0 * ff_15[k]
                  + f_1 * gd0_14[k]
                  - f_2 * gd1_14[k]
                  + pb_y[k] * gf_27[k];

        t_47[k] = pb_z[k] * gf_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_z, ff_17, fg_22, fg_25, gd0_14, \
                         gd0_15, gd1_14, gd1_15, gf_28, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_z[k] * gf_28[k];

        t_49[k] = f_1 * gd0_15[k]
                  - f_2 * gd1_15[k]
                  + pb_z[k] * gf_29[k];

        t_50[k] = pa_z[k] * fg_22[k];

        t_51[k] = f_0 * ff_17[k]
                  + pa_z[k] * fg_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_x, gd0_17, gd0_18, gd0_19, gd1_17, gd1_18, \
                         gd1_19, gf_31, gf_32, gf_33, gf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * gd0_17[k]
                  - f_2 * gd1_17[k]
                  + pb_x[k] * gf_31[k];

        t_53[k] = f_3 * gd0_18[k]
                  - f_4 * gd1_18[k]
                  + pb_x[k] * gf_32[k];

        t_54[k] = f_3 * gd0_19[k]
                  - f_4 * gd1_19[k]
                  + pb_x[k] * gf_33[k];

        t_55[k] = pb_x[k] * gf_34[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_z, pb_x, pb_y, dg0_3, dg1_6, ff_19, ff_20, \
                         fg_26, gd0_19, gd1_19, gf_35, gf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * gf_36[k];

        t_57[k] = f_5 * dg0_3[k]
                  - f_6 * dg1_6[k]
                  + pa_z[k] * fg_26[k];

        t_58[k] = f_7 * ff_19[k]
                  + f_3 * gd0_19[k]
                  - f_4 * gd1_19[k]
                  + pb_y[k] * gf_35[k];

        t_59[k] = f_7 * ff_20[k]
                  + pb_y[k] * gf_36[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pb_y, dg0_8, dg1_14, ff_24, ff_25, \
                         ff_26, fg_29, fg_35, fg_36, gf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_5 * dg0_8[k]
                  - f_6 * dg1_14[k]
                  + pa_y[k] * fg_29[k];

        t_61[k] = f_0 * ff_24[k]
                  + pa_y[k] * fg_35[k];

        t_62[k] = f_7 * ff_25[k]
                  + pa_y[k] * fg_36[k];

        t_63[k] = f_8 * ff_26[k]
                  + pb_y[k] * gf_38[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, fg_38, gd0_21, gd0_22, \
                         gd1_21, gd1_22, gf_39, gf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pa_y[k] * fg_38[k];

        t_65[k] = f_1 * gd0_21[k]
                  - f_2 * gd1_21[k]
                  + pb_x[k] * gf_39[k];

        t_66[k] = pb_y[k] * gf_39[k];

        t_67[k] = f_3 * gd0_22[k]
                  - f_4 * gd1_22[k]
                  + pb_x[k] * gf_40[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, pb_x, pb_y, gd0_22, gd0_23, \
                         gd1_22, gd1_23, gf_41, gf_42, gf_43, gf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pb_x[k] * gf_41[k];

        t_69[k] = pb_x[k] * gf_42[k];

        t_70[k] = pb_x[k] * gf_44[k];

        t_71[k] = f_1 * gd0_22[k]
                  - f_2 * gd1_22[k]
                  + pb_y[k] * gf_42[k];

        t_72[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pb_y[k] * gf_43[k];

        t_73[k] = pb_y[k] * gf_44[k];
    }

#pragma omp simd aligned(t_74, pb_z, ff_26, gd0_23, gd1_23, gf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * ff_26[k]
                  + f_1 * gd0_23[k]
                  - f_2 * gd1_23[k]
                  + pb_z[k] * gf_44[k];
    }
}

auto
compute_prim_gg_electron_repulsion_21(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg0, const size_t dg1,
                                      const size_t ff, const size_t fg, const size_t gd0,
                                      const size_t gd1, const size_t gf, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / p;

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

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_6 = buffer.data(dg0 + 6);
    const auto *dg0_14 = buffer.data(dg0 + 14);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_3 = buffer.data(dg1 + 3);
    const auto *dg1_8 = buffer.data(dg1 + 8);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_26 = buffer.data(ff + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_19 = buffer.data(gd0 + 19);
    const auto *gd0_21 = buffer.data(gd0 + 21);
    const auto *gd0_22 = buffer.data(gd0 + 22);
    const auto *gd0_23 = buffer.data(gd0 + 23);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_21 = buffer.data(gd1 + 21);
    const auto *gd1_22 = buffer.data(gd1 + 22);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, ff_0, gd0_0, gd1_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd0_0[k]
                 - f_2 * gd1_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, fg_0, gd0_1, gd0_2, gd1_1, \
                         gd1_2, gf_3, gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * gd0_1[k]
                 - f_2 * gd1_1[k]
                 + pb_y[k] * gf_3[k];

        t_6[k] = f_3 * gd0_2[k]
                 - f_4 * gd1_2[k]
                 + pb_y[k] * gf_4[k];

        t_7[k] = pb_y[k] * gf_5[k];

        t_8[k] = f_1 * gd0_2[k]
                 - f_2 * gd1_2[k]
                 + pb_z[k] * gf_5[k];

        t_9[k] = pa_y[k] * fg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pa_z, pb_x, dg0_0, dg1_0, ff_8, fg_0, fg_1, \
                         gd0_6, gd1_6, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * fg_0[k];

        t_11[k] = f_5 * dg0_0[k]
                  - f_6 * dg1_0[k]
                  + pa_y[k] * fg_1[k];

        t_12[k] = f_7 * ff_8[k]
                  + f_3 * gd0_6[k]
                  - f_4 * gd1_6[k]
                  + pb_x[k] * gf_9[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pa_z, pb_x, dg0_0, dg0_6, dg1_0, dg1_3, ff_9, \
                         fg_2, fg_3, gf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_7 * ff_9[k]
                  + pb_x[k] * gf_10[k];

        t_14[k] = f_5 * dg0_6[k]
                  - f_6 * dg1_3[k]
                  + pa_x[k] * fg_3[k];

        t_15[k] = f_5 * dg0_0[k]
                  - f_6 * dg1_0[k]
                  + pa_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_x, dg0_14, dg1_8, ff_10, ff_11, \
                         fg_4, fg_5, gd0_10, gd1_10, gf_12, gf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * ff_10[k]
                  + f_3 * gd0_10[k]
                  - f_4 * gd1_10[k]
                  + pb_x[k] * gf_12[k];

        t_17[k] = f_7 * ff_11[k]
                  + pb_x[k] * gf_13[k];

        t_18[k] = f_5 * dg0_14[k]
                  - f_6 * dg1_8[k]
                  + pa_x[k] * fg_4[k];

        t_19[k] = pa_x[k] * fg_5[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, fg_8, gd0_13, gd0_14, gd0_15, \
                         gd1_13, gd1_14, gd1_15, gf_16, gf_17, gf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_x[k] * fg_8[k];

        t_21[k] = f_1 * gd0_13[k]
                  - f_2 * gd1_13[k]
                  + pb_x[k] * gf_16[k];

        t_22[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_x[k] * gf_17[k];

        t_23[k] = f_3 * gd0_15[k]
                  - f_4 * gd1_15[k]
                  + pb_x[k] * gf_18[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pb_x, pb_y, pb_z, ff_15, gd0_14, \
                         gd1_14, gf_19, gf_20, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_x[k] * gf_19[k];

        t_25[k] = pb_x[k] * gf_21[k];

        t_26[k] = f_0 * ff_15[k]
                  + f_1 * gd0_14[k]
                  - f_2 * gd1_14[k]
                  + pb_y[k] * gf_19[k];

        t_27[k] = pb_z[k] * gf_19[k];

        t_28[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pb_z[k] * gf_20[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_z, pb_z, dg0_6, dg1_3, fg_5, fg_6, gd0_15, \
                         gd1_15, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * gd0_15[k]
                  - f_2 * gd1_15[k]
                  + pb_z[k] * gf_21[k];

        t_30[k] = pa_z[k] * fg_5[k];

        t_31[k] = f_5 * dg0_6[k]
                  - f_6 * dg1_3[k]
                  + pa_z[k] * fg_6[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_y, pb_y, dg0_14, dg1_8, ff_19, ff_20, \
                         fg_7, fg_8, gd0_19, gd1_19, gf_24, gf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_7 * ff_19[k]
                  + f_3 * gd0_19[k]
                  - f_4 * gd1_19[k]
                  + pb_y[k] * gf_24[k];

        t_33[k] = f_7 * ff_20[k]
                  + pb_y[k] * gf_25[k];

        t_34[k] = f_5 * dg0_14[k]
                  - f_6 * dg1_8[k]
                  + pa_y[k] * fg_7[k];

        t_35[k] = pa_y[k] * fg_8[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_x, gd0_21, gd0_22, gd0_23, gd1_21, gd1_22, \
                         gd1_23, gf_27, gf_28, gf_29, gf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * gd0_21[k]
                  - f_2 * gd1_21[k]
                  + pb_x[k] * gf_27[k];

        t_37[k] = f_3 * gd0_22[k]
                  - f_4 * gd1_22[k]
                  + pb_x[k] * gf_28[k];

        t_38[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pb_x[k] * gf_29[k];

        t_39[k] = pb_x[k] * gf_30[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pb_x, pb_y, pb_z, ff_26, gd0_22, \
                         gd0_23, gd1_22, gd1_23, gf_30, gf_31, gf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_x[k] * gf_32[k];

        t_41[k] = f_1 * gd0_22[k]
                  - f_2 * gd1_22[k]
                  + pb_y[k] * gf_30[k];

        t_42[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pb_y[k] * gf_31[k];

        t_43[k] = pb_y[k] * gf_32[k];

        t_44[k] = f_0 * ff_26[k]
                  + f_1 * gd0_23[k]
                  - f_2 * gd1_23[k]
                  + pb_z[k] * gf_32[k];
    }
}

}  // namespace simdt2ceri
