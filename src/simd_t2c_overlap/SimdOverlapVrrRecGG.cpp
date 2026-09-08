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


#include "SimdOverlapVrrRecGG.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_gg_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dg, const size_t ff, const size_t fg,
                          const size_t gd, const size_t gf, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_16 = buffer.data(dg + 16);

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
    const auto *ff_29 = buffer.data(ff + 29);
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
    const auto *ff_49 = buffer.data(ff + 49);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_51 = buffer.data(ff + 51);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_53 = buffer.data(ff + 53);

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
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_41 = buffer.data(gd + 41);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_43 = buffer.data(gd + 43);
    const auto *gd_44 = buffer.data(gd + 44);
    const auto *gd_45 = buffer.data(gd + 45);

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
    const auto *gf_99 = buffer.data(gf + 99);
    const auto *gf_100 = buffer.data(gf + 100);
    const auto *gf_101 = buffer.data(gf + 101);
    const auto *gf_102 = buffer.data(gf + 102);
    const auto *gf_103 = buffer.data(gf + 103);
    const auto *gf_104 = buffer.data(gf + 104);
    const auto *gf_105 = buffer.data(gf + 105);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, ff_0, gd_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_2 * gd_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = pb_y[k] * gf_2[k];

        t_5[k] = f_2 * gd_0[k]
                 + pb_z[k] * gf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pb_x, pb_y, pb_z, ff_3, ff_4, gd_1, \
                         gf_3, gf_4, gf_5, gf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * ff_3[k]
                 + pb_x[k] * gf_5[k];

        t_7[k] = pb_z[k] * gf_3[k];

        t_8[k] = pb_y[k] * gf_4[k];

        t_9[k] = f_0 * ff_4[k]
                 + pb_x[k] * gf_7[k];

        t_10[k] = f_1 * gd_1[k]
                  + pb_y[k] * gf_5[k];

        t_11[k] = pb_z[k] * gf_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pa_y, pb_y, pb_z, ff_0, fg_0, \
                         gd_2, gf_6, gf_7, gf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_2 * gd_2[k]
                  + pb_y[k] * gf_6[k];

        t_13[k] = pb_y[k] * gf_7[k];

        t_14[k] = f_1 * gd_2[k]
                  + pb_z[k] * gf_7[k];

        t_15[k] = pa_y[k] * fg_0[k];

        t_16[k] = f_2 * ff_0[k]
                  + pb_y[k] * gf_8[k];

        t_17[k] = pb_z[k] * gf_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_y, pb_x, pb_z, ff_1, ff_6, fg_1, \
                         fg_2, gf_9, gf_10, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * ff_1[k]
                  + pa_y[k] * fg_1[k];

        t_19[k] = pb_z[k] * gf_9[k];

        t_20[k] = pa_y[k] * fg_2[k];

        t_21[k] = f_1 * ff_6[k]
                  + pb_x[k] * gf_11[k];

        t_22[k] = pb_z[k] * gf_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pa_y, pb_x, pb_z, dg_3, ff_7, fg_4, \
                         fg_11, gf_11, gf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * ff_7[k]
                  + pb_x[k] * gf_13[k];

        t_24[k] = pa_y[k] * fg_4[k];

        t_25[k] = f_3 * dg_3[k]
                  + pa_x[k] * fg_11[k];

        t_26[k] = pb_z[k] * gf_11[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_y, pa_z, pb_y, pb_z, ff_4, fg_0, \
                         fg_6, gd_4, gf_12, gf_14, gf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_2 * gd_4[k]
                  + pb_z[k] * gf_12[k];

        t_28[k] = f_2 * ff_4[k]
                  + pb_y[k] * gf_14[k];

        t_29[k] = pa_y[k] * fg_6[k];

        t_30[k] = pa_z[k] * fg_0[k];

        t_31[k] = pb_y[k] * gf_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, ff_0, ff_2, fg_1, \
                         fg_2, fg_3, gf_15, gf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_2 * ff_0[k]
                  + pb_z[k] * gf_15[k];

        t_33[k] = pa_z[k] * fg_1[k];

        t_34[k] = pb_y[k] * gf_16[k];

        t_35[k] = f_3 * ff_2[k]
                  + pa_z[k] * fg_2[k];

        t_36[k] = pa_z[k] * fg_3[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, ff_11, ff_12, fg_5, \
                         gd_7, gf_17, gf_18, gf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_1 * ff_11[k]
                  + pb_x[k] * gf_18[k];

        t_38[k] = pb_y[k] * gf_17[k];

        t_39[k] = f_1 * ff_12[k]
                  + pb_x[k] * gf_20[k];

        t_40[k] = pa_z[k] * fg_5[k];

        t_41[k] = f_3 * gd_7[k]
                  + pb_y[k] * gf_18[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pa_y, pb_y, dg_0, dg_4, fg_7, fg_16, \
                         gd_8, gf_19, gf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_2 * gd_8[k]
                  + pb_y[k] * gf_19[k];

        t_43[k] = pb_y[k] * gf_20[k];

        t_44[k] = f_3 * dg_4[k]
                  + pa_x[k] * fg_16[k];

        t_45[k] = f_2 * dg_0[k]
                  + pa_y[k] * fg_7[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pb_x, pb_y, pb_z, ff_5, ff_14, gd_9, \
                         gd_10, gf_21, gf_22, gf_23, gf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_3 * ff_5[k]
                  + pb_y[k] * gf_21[k];

        t_47[k] = pb_z[k] * gf_21[k];

        t_48[k] = f_3 * ff_14[k]
                  + f_2 * gd_10[k]
                  + pb_x[k] * gf_24[k];

        t_49[k] = pb_z[k] * gf_22[k];

        t_50[k] = f_2 * gd_9[k]
                  + pb_z[k] * gf_23[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_x, pb_z, ff_15, ff_16, ff_17, gf_24, \
                         gf_25, gf_27, gf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * ff_15[k]
                  + pb_x[k] * gf_25[k];

        t_52[k] = pb_z[k] * gf_24[k];

        t_53[k] = f_3 * ff_16[k]
                  + pb_x[k] * gf_27[k];

        t_54[k] = f_3 * ff_17[k]
                  + pb_x[k] * gf_28[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_x, pb_y, pb_z, dg_7, ff_8, fg_21, \
                         gd_10, gd_11, gf_25, gf_26, gf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_2 * dg_7[k]
                  + pa_x[k] * fg_21[k];

        t_56[k] = pb_z[k] * gf_25[k];

        t_57[k] = f_2 * gd_10[k]
                  + pb_z[k] * gf_26[k];

        t_58[k] = f_3 * ff_8[k]
                  + pb_y[k] * gf_28[k];

        t_59[k] = f_1 * gd_11[k]
                  + pb_z[k] * gf_28[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, ff_10, fg_8, \
                         fg_9, fg_12, fg_13, fg_14, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * fg_12[k];

        t_61[k] = pa_z[k] * fg_8[k];

        t_62[k] = pa_y[k] * fg_13[k];

        t_63[k] = pa_z[k] * fg_9[k];

        t_64[k] = f_2 * ff_10[k]
                  + pb_y[k] * gf_29[k];

        t_65[k] = pa_y[k] * fg_14[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, ff_19, ff_20, fg_10, \
                         fg_11, fg_15, gf_31, gf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * fg_10[k];

        t_67[k] = f_3 * ff_19[k]
                  + pb_x[k] * gf_31[k];

        t_68[k] = f_3 * ff_20[k]
                  + pb_x[k] * gf_32[k];

        t_69[k] = pa_y[k] * fg_15[k];

        t_70[k] = pa_z[k] * fg_11[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_x, pa_y, pb_y, pb_z, dg_9, ff_6, ff_12, \
                         fg_16, fg_22, gf_30, gf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_2 * ff_6[k]
                  + pb_z[k] * gf_30[k];

        t_72[k] = f_2 * dg_9[k]
                  + pa_x[k] * fg_22[k];

        t_73[k] = f_2 * ff_12[k]
                  + pb_y[k] * gf_33[k];

        t_74[k] = pa_y[k] * fg_16[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pa_z, pb_y, pb_z, dg_0, ff_9, fg_12, \
                         gd_14, gf_34, gf_35, gf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_2 * dg_0[k]
                  + pa_z[k] * fg_12[k];

        t_76[k] = pb_y[k] * gf_34[k];

        t_77[k] = f_3 * ff_9[k]
                  + pb_z[k] * gf_34[k];

        t_78[k] = f_2 * gd_14[k]
                  + pb_y[k] * gf_35[k];

        t_79[k] = pb_y[k] * gf_36[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pb_x, pb_y, ff_23, ff_24, ff_25, ff_26, \
                         gd_17, gf_37, gf_38, gf_39, gf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_3 * ff_23[k]
                  + f_2 * gd_17[k]
                  + pb_x[k] * gf_37[k];

        t_81[k] = f_3 * ff_24[k]
                  + pb_x[k] * gf_38[k];

        t_82[k] = f_3 * ff_25[k]
                  + pb_x[k] * gf_39[k];

        t_83[k] = pb_y[k] * gf_37[k];

        t_84[k] = f_3 * ff_26[k]
                  + pb_x[k] * gf_41[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pa_x, pb_y, dg_16, fg_27, gd_15, gd_16, \
                         gd_17, gf_38, gf_39, gf_40, gf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_1 * gd_15[k]
                  + pb_y[k] * gf_38[k];

        t_86[k] = f_3 * gd_16[k]
                  + pb_y[k] * gf_39[k];

        t_87[k] = f_2 * gd_17[k]
                  + pb_y[k] * gf_40[k];

        t_88[k] = pb_y[k] * gf_41[k];

        t_89[k] = f_2 * dg_16[k]
                  + pa_x[k] * fg_27[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pa_x, pb_y, pb_z, ff_13, ff_27, ff_29, \
                         fg_28, fg_30, gf_42, gf_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * ff_27[k]
                  + pa_x[k] * fg_28[k];

        t_91[k] = f_1 * ff_13[k]
                  + pb_y[k] * gf_42[k];

        t_92[k] = pb_z[k] * gf_42[k];

        t_93[k] = f_3 * ff_29[k]
                  + pa_x[k] * fg_30[k];

        t_94[k] = pb_z[k] * gf_43[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pb_x, pb_z, ff_31, ff_33, ff_34, gd_18, \
                         gf_44, gf_45, gf_46, gf_47, gf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_2 * gd_18[k]
                  + pb_z[k] * gf_44[k];

        t_96[k] = f_2 * ff_31[k]
                  + pb_x[k] * gf_46[k];

        t_97[k] = pb_z[k] * gf_45[k];

        t_98[k] = f_2 * ff_33[k]
                  + pb_x[k] * gf_47[k];

        t_99[k] = f_2 * ff_34[k]
                  + pb_x[k] * gf_48[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, t_105, pa_x, pa_z, pb_z, fg_17, \
                         fg_33, fg_34, fg_35, fg_36, gf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_x[k] * fg_33[k];

        t_101[k] = pb_z[k] * gf_46[k];

        t_102[k] = pa_x[k] * fg_34[k];

        t_103[k] = pa_x[k] * fg_35[k];

        t_104[k] = pa_x[k] * fg_36[k];

        t_105[k] = pa_z[k] * fg_17[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pa_z, pb_y, pb_z, ff_13, ff_18, fg_18, \
                         fg_19, gf_49, gf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = pa_z[k] * fg_18[k];

        t_107[k] = f_2 * ff_13[k]
                   + pb_z[k] * gf_49[k];

        t_108[k] = pa_z[k] * fg_19[k];

        t_109[k] = f_3 * ff_18[k]
                   + pb_y[k] * gf_50[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_x, pa_z, pb_x, ff_35, ff_37, ff_38, \
                         fg_20, fg_37, gf_51, gf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_3 * ff_35[k]
                   + pa_x[k] * fg_37[k];

        t_111[k] = pa_z[k] * fg_20[k];

        t_112[k] = f_2 * ff_37[k]
                   + pb_x[k] * gf_51[k];

        t_113[k] = f_2 * ff_38[k]
                   + pb_x[k] * gf_52[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, pa_x, pb_x, ff_39, fg_38, \
                         fg_39, fg_40, fg_41, fg_42, gf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_2 * ff_39[k]
                   + pb_x[k] * gf_53[k];

        t_115[k] = pa_x[k] * fg_38[k];

        t_116[k] = pa_x[k] * fg_39[k];

        t_117[k] = pa_x[k] * fg_40[k];

        t_118[k] = pa_x[k] * fg_41[k];

        t_119[k] = pa_x[k] * fg_42[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pa_x, pa_y, pb_y, ff_21, ff_22, \
                         ff_40, fg_23, fg_24, fg_43, gf_54, gf_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_y[k] * fg_23[k];

        t_121[k] = f_2 * ff_21[k]
                   + pb_y[k] * gf_54[k];

        t_122[k] = pa_y[k] * fg_24[k];

        t_123[k] = f_3 * ff_40[k]
                   + pa_x[k] * fg_43[k];

        t_124[k] = f_2 * ff_22[k]
                   + pb_y[k] * gf_55[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, pa_y, pb_x, ff_41, ff_42, ff_43, \
                         fg_25, fg_26, gf_56, gf_57, gf_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pa_y[k] * fg_25[k];

        t_126[k] = f_2 * ff_41[k]
                   + pb_x[k] * gf_56[k];

        t_127[k] = f_2 * ff_42[k]
                   + pb_x[k] * gf_57[k];

        t_128[k] = f_2 * ff_43[k]
                   + pb_x[k] * gf_58[k];

        t_129[k] = pa_y[k] * fg_26[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, t_135, pa_x, ff_45, fg_44, fg_45, \
                         fg_46, fg_47, fg_48, fg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pa_x[k] * fg_44[k];

        t_131[k] = pa_x[k] * fg_45[k];

        t_132[k] = pa_x[k] * fg_46[k];

        t_133[k] = pa_x[k] * fg_47[k];

        t_134[k] = pa_x[k] * fg_48[k];

        t_135[k] = f_0 * ff_45[k]
                   + pa_x[k] * fg_49[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, pa_x, pb_y, pb_z, ff_21, ff_49, \
                         fg_54, gd_21, gf_59, gf_60, gf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pb_y[k] * gf_59[k];

        t_137[k] = f_1 * ff_21[k]
                   + pb_z[k] * gf_59[k];

        t_138[k] = f_2 * gd_21[k]
                   + pb_y[k] * gf_60[k];

        t_139[k] = pb_y[k] * gf_61[k];

        t_140[k] = f_3 * ff_49[k]
                   + pa_x[k] * fg_54[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, pa_x, pb_x, pb_y, ff_50, ff_51, \
                         ff_53, fg_55, gf_62, gf_63, gf_64, gf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_2 * ff_50[k]
                   + pb_x[k] * gf_63[k];

        t_142[k] = f_2 * ff_51[k]
                   + pb_x[k] * gf_64[k];

        t_143[k] = pb_y[k] * gf_62[k];

        t_144[k] = f_2 * ff_53[k]
                   + pb_x[k] * gf_65[k];

        t_145[k] = pa_x[k] * fg_55[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, pa_x, pb_x, pb_y, fg_56, fg_57, \
                         fg_58, gd_23, gf_65, gf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = pa_x[k] * fg_56[k];

        t_147[k] = pa_x[k] * fg_57[k];

        t_148[k] = pb_y[k] * gf_65[k];

        t_149[k] = pa_x[k] * fg_58[k];

        t_150[k] = f_1 * gd_23[k]
                   + pb_x[k] * gf_66[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, t_156, pb_x, pb_z, gd_24, gd_25, \
                         gd_26, gf_66, gf_67, gf_68, gf_69, gf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_3 * gd_24[k]
                   + pb_x[k] * gf_67[k];

        t_152[k] = pb_z[k] * gf_66[k];

        t_153[k] = f_2 * gd_25[k]
                   + pb_x[k] * gf_68[k];

        t_154[k] = pb_z[k] * gf_67[k];

        t_155[k] = f_2 * gd_26[k]
                   + pb_x[k] * gf_69[k];

        t_156[k] = pb_x[k] * gf_70[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, t_162, pb_x, pb_y, pb_z, ff_31, \
                         gd_25, gf_70, gf_71, gf_72, gf_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = pb_x[k] * gf_71[k];

        t_158[k] = pb_x[k] * gf_72[k];

        t_159[k] = pb_x[k] * gf_73[k];

        t_160[k] = f_0 * ff_31[k]
                   + f_1 * gd_25[k]
                   + pb_y[k] * gf_70[k];

        t_161[k] = pb_z[k] * gf_70[k];

        t_162[k] = f_2 * gd_25[k]
                   + pb_z[k] * gf_71[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pa_z, pb_x, pb_y, pb_z, ff_34, \
                         fg_28, fg_29, gd_26, gd_27, gf_73, gf_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_0 * ff_34[k]
                   + pb_y[k] * gf_73[k];

        t_164[k] = f_1 * gd_26[k]
                   + pb_z[k] * gf_73[k];

        t_165[k] = pa_z[k] * fg_28[k];

        t_166[k] = pa_z[k] * fg_29[k];

        t_167[k] = f_3 * gd_27[k]
                   + pb_x[k] * gf_74[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, t_173, pa_z, pb_x, fg_30, gd_29, \
                         gd_30, gf_75, gf_76, gf_77, gf_78, gf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_z[k] * fg_30[k];

        t_169[k] = f_2 * gd_29[k]
                   + pb_x[k] * gf_75[k];

        t_170[k] = f_2 * gd_30[k]
                   + pb_x[k] * gf_76[k];

        t_171[k] = pb_x[k] * gf_77[k];

        t_172[k] = pb_x[k] * gf_78[k];

        t_173[k] = pb_x[k] * gf_79[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, pa_z, pb_x, pb_y, pb_z, ff_31, \
                         ff_32, ff_39, fg_33, fg_34, gf_77, gf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pb_x[k] * gf_80[k];

        t_175[k] = pa_z[k] * fg_33[k];

        t_176[k] = f_2 * ff_31[k]
                   + pb_z[k] * gf_77[k];

        t_177[k] = f_3 * ff_32[k]
                   + pa_z[k] * fg_34[k];

        t_178[k] = f_1 * ff_39[k]
                   + pb_y[k] * gf_80[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pa_y, pb_x, dg_10, fg_42, gd_31, gd_32, \
                         gd_33, gf_81, gf_82, gf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_3 * dg_10[k]
                   + pa_y[k] * fg_42[k];

        t_180[k] = f_1 * gd_31[k]
                   + pb_x[k] * gf_81[k];

        t_181[k] = f_3 * gd_32[k]
                   + pb_x[k] * gf_82[k];

        t_182[k] = f_3 * gd_33[k]
                   + pb_x[k] * gf_83[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, t_188, pb_x, gd_34, gd_35, gd_36, \
                         gf_84, gf_85, gf_86, gf_87, gf_88, gf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_2 * gd_34[k]
                   + pb_x[k] * gf_84[k];

        t_184[k] = f_2 * gd_35[k]
                   + pb_x[k] * gf_85[k];

        t_185[k] = f_2 * gd_36[k]
                   + pb_x[k] * gf_86[k];

        t_186[k] = pb_x[k] * gf_87[k];

        t_187[k] = pb_x[k] * gf_88[k];

        t_188[k] = pb_x[k] * gf_89[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_z, pb_x, pb_y, pb_z, dg_7, ff_36, \
                         ff_43, fg_38, gd_36, gf_87, gf_89, gf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pb_x[k] * gf_90[k];

        t_190[k] = f_2 * dg_7[k]
                   + pa_z[k] * fg_38[k];

        t_191[k] = f_3 * ff_36[k]
                   + pb_z[k] * gf_87[k];

        t_192[k] = f_3 * ff_43[k]
                   + f_2 * gd_36[k]
                   + pb_y[k] * gf_89[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, pa_y, pb_x, pb_y, dg_16, ff_44, \
                         fg_48, fg_49, fg_51, gd_37, gf_90, gf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * ff_44[k]
                   + pb_y[k] * gf_90[k];

        t_194[k] = f_2 * dg_16[k]
                   + pa_y[k] * fg_48[k];

        t_195[k] = pa_y[k] * fg_49[k];

        t_196[k] = f_3 * gd_37[k]
                   + pb_x[k] * gf_91[k];

        t_197[k] = pa_y[k] * fg_51[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, t_203, pa_y, pb_x, fg_54, gd_38, \
                         gd_39, gf_92, gf_93, gf_94, gf_95, gf_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_2 * gd_38[k]
                   + pb_x[k] * gf_92[k];

        t_199[k] = f_2 * gd_39[k]
                   + pb_x[k] * gf_93[k];

        t_200[k] = pa_y[k] * fg_54[k];

        t_201[k] = pb_x[k] * gf_94[k];

        t_202[k] = pb_x[k] * gf_95[k];

        t_203[k] = pb_x[k] * gf_96[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_y, pb_x, pb_z, ff_41, ff_50, ff_52, \
                         fg_55, fg_57, gf_94, gf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pb_x[k] * gf_97[k];

        t_205[k] = f_0 * ff_50[k]
                   + pa_y[k] * fg_55[k];

        t_206[k] = f_1 * ff_41[k]
                   + pb_z[k] * gf_94[k];

        t_207[k] = f_3 * ff_52[k]
                   + pa_y[k] * fg_57[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pa_y, pb_x, pb_y, ff_53, fg_58, \
                         gd_41, gd_42, gf_97, gf_98, gf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_2 * ff_53[k]
                   + pb_y[k] * gf_97[k];

        t_209[k] = pa_y[k] * fg_58[k];

        t_210[k] = f_1 * gd_41[k]
                   + pb_x[k] * gf_98[k];

        t_211[k] = pb_y[k] * gf_98[k];

        t_212[k] = f_3 * gd_42[k]
                   + pb_x[k] * gf_99[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, t_218, pb_x, pb_y, gd_43, gd_45, \
                         gf_99, gf_100, gf_101, gf_102, gf_103, \
                         gf_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_2 * gd_43[k]
                   + pb_x[k] * gf_100[k];

        t_214[k] = pb_y[k] * gf_99[k];

        t_215[k] = f_2 * gd_45[k]
                   + pb_x[k] * gf_101[k];

        t_216[k] = pb_x[k] * gf_102[k];

        t_217[k] = pb_x[k] * gf_103[k];

        t_218[k] = pb_x[k] * gf_104[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, pb_x, pb_y, gd_43, gd_44, gd_45, \
                         gf_102, gf_103, gf_104, gf_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = pb_x[k] * gf_105[k];

        t_220[k] = f_1 * gd_43[k]
                   + pb_y[k] * gf_102[k];

        t_221[k] = f_3 * gd_44[k]
                   + pb_y[k] * gf_103[k];

        t_222[k] = f_2 * gd_45[k]
                   + pb_y[k] * gf_104[k];

        t_223[k] = pb_y[k] * gf_105[k];
    }

#pragma omp simd aligned(t_224, pb_z, ff_53, gd_45, gf_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_0 * ff_53[k]
                   + f_1 * gd_45[k]
                   + pb_z[k] * gf_105[k];
    }
}

auto
compute_prim_gg_overlap_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dg, const size_t ff, const size_t fg,
                          const size_t gd, const size_t gf, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_38 = buffer.data(dg + 38);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_42 = buffer.data(ff + 42);
    const auto *ff_43 = buffer.data(ff + 43);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_58 = buffer.data(fg + 58);
    const auto *fg_59 = buffer.data(fg + 59);
    const auto *fg_60 = buffer.data(fg + 60);
    const auto *fg_61 = buffer.data(fg + 61);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_63 = buffer.data(fg + 63);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_66 = buffer.data(fg + 66);
    const auto *fg_67 = buffer.data(fg + 67);
    const auto *fg_68 = buffer.data(fg + 68);
    const auto *fg_69 = buffer.data(fg + 69);
    const auto *fg_70 = buffer.data(fg + 70);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_79 = buffer.data(fg + 79);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_81 = buffer.data(fg + 81);
    const auto *fg_83 = buffer.data(fg + 83);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, ff_0, ff_3, gd_0, \
                         gf_0, gf_1, gf_2, gf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_2 * gd_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_2 * gd_0[k]
                 + pb_z[k] * gf_2[k];

        t_5[k] = f_0 * ff_3[k]
                 + pb_x[k] * gf_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, ff_5, gd_1, gd_2, gf_3, \
                         gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * ff_5[k]
                 + pb_x[k] * gf_5[k];

        t_7[k] = f_1 * gd_1[k]
                 + pb_y[k] * gf_3[k];

        t_8[k] = f_2 * gd_2[k]
                 + pb_y[k] * gf_4[k];

        t_9[k] = pb_y[k] * gf_5[k];

        t_10[k] = f_1 * gd_2[k]
                  + pb_z[k] * gf_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_x, pb_y, ff_0, ff_1, ff_7, \
                         fg_0, fg_3, fg_4, gf_6, gf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_y[k] * fg_0[k];

        t_12[k] = f_2 * ff_0[k]
                  + pb_y[k] * gf_6[k];

        t_13[k] = f_3 * ff_1[k]
                  + pa_y[k] * fg_3[k];

        t_14[k] = pa_y[k] * fg_4[k];

        t_15[k] = f_1 * ff_7[k]
                  + pb_x[k] * gf_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, pb_z, dg_9, ff_5, fg_11, gd_3, \
                         gf_7, gf_8, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dg_9[k]
                  + pa_x[k] * fg_11[k];

        t_17[k] = pb_z[k] * gf_7[k];

        t_18[k] = f_2 * gd_3[k]
                  + pb_z[k] * gf_8[k];

        t_19[k] = f_2 * ff_5[k]
                  + pb_y[k] * gf_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pa_z, pb_y, pb_z, ff_0, ff_2, \
                         fg_0, fg_4, fg_7, gf_10, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_y[k] * fg_7[k];

        t_21[k] = pa_z[k] * fg_0[k];

        t_22[k] = f_2 * ff_0[k]
                  + pb_z[k] * gf_10[k];

        t_23[k] = pb_y[k] * gf_11[k];

        t_24[k] = f_3 * ff_2[k]
                  + pa_z[k] * fg_4[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_x, pb_x, pb_y, dg_13, ff_13, fg_20, \
                         gd_5, gd_6, gf_12, gf_13, gf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * ff_13[k]
                  + pb_x[k] * gf_14[k];

        t_26[k] = f_3 * gd_5[k]
                  + pb_y[k] * gf_12[k];

        t_27[k] = f_2 * gd_6[k]
                  + pb_y[k] * gf_13[k];

        t_28[k] = pb_y[k] * gf_14[k];

        t_29[k] = f_3 * dg_13[k]
                  + pa_x[k] * fg_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_x, pb_y, pb_z, dg_0, ff_6, ff_16, \
                         fg_8, gd_8, gf_15, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * dg_0[k]
                  + pa_y[k] * fg_8[k];

        t_31[k] = f_3 * ff_6[k]
                  + pb_y[k] * gf_15[k];

        t_32[k] = pb_z[k] * gf_15[k];

        t_33[k] = f_3 * ff_16[k]
                  + f_2 * gd_8[k]
                  + pb_x[k] * gf_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_x, pb_x, pb_z, dg_18, ff_17, fg_25, \
                         gd_7, gd_8, gf_16, gf_18, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_2 * gd_7[k]
                  + pb_z[k] * gf_16[k];

        t_35[k] = f_3 * ff_17[k]
                  + pb_x[k] * gf_18[k];

        t_36[k] = f_2 * dg_18[k]
                  + pa_x[k] * fg_25[k];

        t_37[k] = pb_z[k] * gf_18[k];

        t_38[k] = f_2 * gd_8[k]
                  + pb_z[k] * gf_19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_y, pa_z, pb_y, pb_z, ff_9, fg_9, \
                         fg_16, fg_17, gd_9, gf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_3 * ff_9[k]
                  + pb_y[k] * gf_20[k];

        t_40[k] = f_1 * gd_9[k]
                  + pb_z[k] * gf_20[k];

        t_41[k] = pa_y[k] * fg_16[k];

        t_42[k] = pa_z[k] * fg_9[k];

        t_43[k] = pa_y[k] * fg_17[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pa_z, pb_y, pb_z, dg_25, ff_7, ff_13, \
                         fg_11, fg_33, gf_21, gf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_z[k] * fg_11[k];

        t_45[k] = f_2 * ff_7[k]
                  + pb_z[k] * gf_21[k];

        t_46[k] = f_2 * dg_25[k]
                  + pa_x[k] * fg_33[k];

        t_47[k] = f_2 * ff_13[k]
                  + pb_y[k] * gf_22[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_y, pa_z, pb_y, pb_z, dg_0, ff_10, \
                         fg_15, fg_20, gd_10, gf_23, gf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_y[k] * fg_20[k];

        t_49[k] = f_2 * dg_0[k]
                  + pa_z[k] * fg_15[k];

        t_50[k] = pb_y[k] * gf_23[k];

        t_51[k] = f_3 * ff_10[k]
                  + pb_z[k] * gf_23[k];

        t_52[k] = f_2 * gd_10[k]
                  + pb_y[k] * gf_24[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_x, pb_y, ff_19, ff_20, gd_11, gd_13, \
                         gf_25, gf_26, gf_27, gf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_y[k] * gf_25[k];

        t_54[k] = f_3 * ff_19[k]
                  + f_2 * gd_13[k]
                  + pb_x[k] * gf_26[k];

        t_55[k] = f_3 * ff_20[k]
                  + pb_x[k] * gf_30[k];

        t_56[k] = f_1 * gd_11[k]
                  + pb_y[k] * gf_27[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pa_x, pb_y, dg_38, ff_21, fg_42, fg_43, \
                         gd_12, gd_13, gf_28, gf_29, gf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_3 * gd_12[k]
                  + pb_y[k] * gf_28[k];

        t_58[k] = f_2 * gd_13[k]
                  + pb_y[k] * gf_29[k];

        t_59[k] = pb_y[k] * gf_30[k];

        t_60[k] = f_2 * dg_38[k]
                  + pa_x[k] * fg_42[k];

        t_61[k] = f_0 * ff_21[k]
                  + pa_x[k] * fg_43[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pb_y, pb_z, ff_14, ff_23, fg_45, gd_14, \
                         gf_31, gf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_1 * ff_14[k]
                  + pb_y[k] * gf_31[k];

        t_63[k] = pb_z[k] * gf_31[k];

        t_64[k] = f_3 * ff_23[k]
                  + pa_x[k] * fg_45[k];

        t_65[k] = f_2 * gd_14[k]
                  + pb_z[k] * gf_32[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, pa_x, pa_z, pb_x, ff_25, fg_21, \
                         fg_51, fg_53, fg_54, fg_55, gf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_2 * ff_25[k]
                  + pb_x[k] * gf_33[k];

        t_67[k] = pa_x[k] * fg_51[k];

        t_68[k] = pa_x[k] * fg_53[k];

        t_69[k] = pa_x[k] * fg_54[k];

        t_70[k] = pa_x[k] * fg_55[k];

        t_71[k] = pa_z[k] * fg_21[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pa_x, pa_z, pb_z, ff_14, ff_29, fg_22, \
                         fg_56, fg_59, fg_60, gf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_2 * ff_14[k]
                  + pb_z[k] * gf_34[k];

        t_73[k] = pa_z[k] * fg_22[k];

        t_74[k] = f_3 * ff_29[k]
                  + pa_x[k] * fg_56[k];

        t_75[k] = pa_x[k] * fg_59[k];

        t_76[k] = pa_x[k] * fg_60[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, t_82, pa_x, pa_y, ff_32, fg_35, fg_36, \
                         fg_37, fg_61, fg_62, fg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_x[k] * fg_61[k];

        t_78[k] = pa_x[k] * fg_62[k];

        t_79[k] = pa_y[k] * fg_35[k];

        t_80[k] = pa_y[k] * fg_36[k];

        t_81[k] = f_3 * ff_32[k]
                  + pa_x[k] * fg_63[k];

        t_82[k] = pa_y[k] * fg_37[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, t_88, pa_x, pb_y, ff_36, fg_65, fg_66, \
                         fg_67, fg_68, fg_70, gf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pa_x[k] * fg_65[k];

        t_84[k] = pa_x[k] * fg_66[k];

        t_85[k] = pa_x[k] * fg_67[k];

        t_86[k] = pa_x[k] * fg_68[k];

        t_87[k] = f_0 * ff_36[k]
                  + pa_x[k] * fg_70[k];

        t_88[k] = pb_y[k] * gf_35[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_x, pb_y, pb_z, ff_18, ff_39, fg_75, gd_15, \
                         gf_35, gf_36, gf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_1 * ff_18[k]
                  + pb_z[k] * gf_35[k];

        t_90[k] = f_2 * gd_15[k]
                  + pb_y[k] * gf_36[k];

        t_91[k] = pb_y[k] * gf_37[k];

        t_92[k] = f_3 * ff_39[k]
                  + pa_x[k] * fg_75[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, pa_x, pb_x, ff_43, fg_79, fg_80, \
                         fg_81, fg_83, gd_16, gf_38, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_2 * ff_43[k]
                  + pb_x[k] * gf_38[k];

        t_94[k] = pa_x[k] * fg_79[k];

        t_95[k] = pa_x[k] * fg_80[k];

        t_96[k] = pa_x[k] * fg_81[k];

        t_97[k] = pa_x[k] * fg_83[k];

        t_98[k] = f_1 * gd_16[k]
                  + pb_x[k] * gf_39[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, pb_x, gd_17, gd_18, gd_19, \
                         gf_40, gf_41, gf_42, gf_43, gf_45, gf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_3 * gd_17[k]
                  + pb_x[k] * gf_40[k];

        t_100[k] = f_2 * gd_18[k]
                   + pb_x[k] * gf_41[k];

        t_101[k] = f_2 * gd_19[k]
                   + pb_x[k] * gf_42[k];

        t_102[k] = pb_x[k] * gf_43[k];

        t_103[k] = pb_x[k] * gf_45[k];

        t_104[k] = pb_x[k] * gf_46[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pb_y, pb_z, ff_25, ff_28, gd_18, \
                         gd_19, gf_43, gf_44, gf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_0 * ff_25[k]
                   + f_1 * gd_18[k]
                   + pb_y[k] * gf_43[k];

        t_106[k] = pb_z[k] * gf_43[k];

        t_107[k] = f_2 * gd_18[k]
                   + pb_z[k] * gf_44[k];

        t_108[k] = f_0 * ff_28[k]
                   + pb_y[k] * gf_46[k];

        t_109[k] = f_1 * gd_19[k]
                   + pb_z[k] * gf_46[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, t_115, pb_x, gd_20, gd_22, gd_23, \
                         gf_47, gf_48, gf_49, gf_51, gf_52, gf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_3 * gd_20[k]
                   + pb_x[k] * gf_47[k];

        t_111[k] = f_2 * gd_22[k]
                   + pb_x[k] * gf_48[k];

        t_112[k] = f_2 * gd_23[k]
                   + pb_x[k] * gf_49[k];

        t_113[k] = pb_x[k] * gf_51[k];

        t_114[k] = pb_x[k] * gf_52[k];

        t_115[k] = pb_x[k] * gf_53[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_z, pb_y, pb_z, ff_25, ff_26, ff_31, \
                         fg_51, fg_53, gf_50, gf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = pa_z[k] * fg_51[k];

        t_117[k] = f_2 * ff_25[k]
                   + pb_z[k] * gf_50[k];

        t_118[k] = f_3 * ff_26[k]
                   + pa_z[k] * fg_53[k];

        t_119[k] = f_1 * ff_31[k]
                   + pb_y[k] * gf_53[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_y, pb_x, dg_27, fg_62, gd_24, gd_25, \
                         gd_26, gf_54, gf_55, gf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_3 * dg_27[k]
                   + pa_y[k] * fg_62[k];

        t_121[k] = f_1 * gd_24[k]
                   + pb_x[k] * gf_54[k];

        t_122[k] = f_3 * gd_25[k]
                   + pb_x[k] * gf_55[k];

        t_123[k] = f_3 * gd_26[k]
                   + pb_x[k] * gf_56[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, t_129, pb_x, gd_27, gd_28, gd_29, \
                         gf_57, gf_58, gf_59, gf_60, gf_61, gf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_2 * gd_27[k]
                   + pb_x[k] * gf_57[k];

        t_125[k] = f_2 * gd_28[k]
                   + pb_x[k] * gf_58[k];

        t_126[k] = f_2 * gd_29[k]
                   + pb_x[k] * gf_59[k];

        t_127[k] = pb_x[k] * gf_60[k];

        t_128[k] = pb_x[k] * gf_61[k];

        t_129[k] = pb_x[k] * gf_62[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pa_z, pb_x, pb_y, pb_z, dg_18, ff_30, \
                         ff_34, fg_58, gd_29, gf_60, gf_62, gf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pb_x[k] * gf_63[k];

        t_131[k] = f_2 * dg_18[k]
                   + pa_z[k] * fg_58[k];

        t_132[k] = f_3 * ff_30[k]
                   + pb_z[k] * gf_60[k];

        t_133[k] = f_3 * ff_34[k]
                   + f_2 * gd_29[k]
                   + pb_y[k] * gf_62[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_y, pb_x, pb_y, dg_38, ff_35, fg_69, \
                         gd_30, gd_31, gf_63, gf_64, gf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_3 * ff_35[k]
                   + pb_y[k] * gf_63[k];

        t_135[k] = f_2 * dg_38[k]
                   + pa_y[k] * fg_69[k];

        t_136[k] = f_3 * gd_30[k]
                   + pb_x[k] * gf_64[k];

        t_137[k] = f_2 * gd_31[k]
                   + pb_x[k] * gf_65[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, pa_y, pb_x, ff_40, fg_79, gd_32, \
                         gf_66, gf_67, gf_68, gf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_2 * gd_32[k]
                   + pb_x[k] * gf_66[k];

        t_139[k] = pb_x[k] * gf_67[k];

        t_140[k] = pb_x[k] * gf_68[k];

        t_141[k] = pb_x[k] * gf_69[k];

        t_142[k] = f_0 * ff_40[k]
                   + pa_y[k] * fg_79[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pa_y, pb_y, pb_z, ff_33, ff_42, ff_43, \
                         fg_81, fg_83, gf_67, gf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_1 * ff_33[k]
                   + pb_z[k] * gf_67[k];

        t_144[k] = f_3 * ff_42[k]
                   + pa_y[k] * fg_81[k];

        t_145[k] = f_2 * ff_43[k]
                   + pb_y[k] * gf_70[k];

        t_146[k] = pa_y[k] * fg_83[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, pb_x, gd_34, gd_35, gd_36, gd_38, \
                         gf_71, gf_72, gf_73, gf_74, gf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_1 * gd_34[k]
                   + pb_x[k] * gf_71[k];

        t_148[k] = f_3 * gd_35[k]
                   + pb_x[k] * gf_72[k];

        t_149[k] = f_2 * gd_36[k]
                   + pb_x[k] * gf_73[k];

        t_150[k] = f_2 * gd_38[k]
                   + pb_x[k] * gf_74[k];

        t_151[k] = pb_x[k] * gf_75[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, t_157, pb_x, pb_y, gd_36, gd_37, \
                         gd_38, gf_75, gf_76, gf_77, gf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = pb_x[k] * gf_76[k];

        t_153[k] = pb_x[k] * gf_78[k];

        t_154[k] = f_1 * gd_36[k]
                   + pb_y[k] * gf_75[k];

        t_155[k] = f_3 * gd_37[k]
                   + pb_y[k] * gf_76[k];

        t_156[k] = f_2 * gd_38[k]
                   + pb_y[k] * gf_77[k];

        t_157[k] = pb_y[k] * gf_78[k];
    }

#pragma omp simd aligned(t_158, pb_z, ff_43, gd_38, gf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_0 * ff_43[k]
                   + f_1 * gd_38[k]
                   + pb_z[k] * gf_78[k];
    }
}

auto
compute_prim_gg_overlap_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dg, const size_t ff, const size_t fg,
                          const size_t gd, const size_t gf, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_31 = buffer.data(dg + 31);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_41 = buffer.data(ff + 41);
    const auto *ff_42 = buffer.data(ff + 42);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_60 = buffer.data(fg + 60);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_64 = buffer.data(fg + 64);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, ff_0, gd_0, gd_1, \
                         gf_0, gf_1, gf_2, gf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_2 * gd_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_2 * gd_0[k]
                 + pb_z[k] * gf_2[k];

        t_5[k] = f_1 * gd_1[k]
                 + pb_y[k] * gf_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_y, pb_y, pb_z, ff_1, fg_0, fg_3, gd_2, \
                         gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * gd_2[k]
                 + pb_y[k] * gf_4[k];

        t_7[k] = pb_y[k] * gf_5[k];

        t_8[k] = f_1 * gd_2[k]
                 + pb_z[k] * gf_5[k];

        t_9[k] = pa_y[k] * fg_0[k];

        t_10[k] = f_3 * ff_1[k]
                  + pa_y[k] * fg_3[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_x, pa_y, pa_z, pb_z, dg_7, fg_0, \
                         fg_8, fg_10, gd_3, gf_6, gf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * dg_7[k]
                  + pa_x[k] * fg_10[k];

        t_12[k] = pb_z[k] * gf_6[k];

        t_13[k] = f_2 * gd_3[k]
                  + pb_z[k] * gf_7[k];

        t_14[k] = pa_y[k] * fg_8[k];

        t_15[k] = pa_z[k] * fg_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_z, pb_y, pb_z, ff_0, ff_2, fg_4, gd_5, \
                         gf_8, gf_9, gf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * ff_0[k]
                  + pb_z[k] * gf_8[k];

        t_17[k] = pb_y[k] * gf_9[k];

        t_18[k] = f_3 * ff_2[k]
                  + pa_z[k] * fg_4[k];

        t_19[k] = f_3 * gd_5[k]
                  + pb_y[k] * gf_10[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pa_y, pb_y, dg_0, dg_9, fg_9, fg_19, \
                         gd_6, gf_11, gf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_2 * gd_6[k]
                  + pb_y[k] * gf_11[k];

        t_21[k] = pb_y[k] * gf_12[k];

        t_22[k] = f_3 * dg_9[k]
                  + pa_x[k] * fg_19[k];

        t_23[k] = f_2 * dg_0[k]
                  + pa_y[k] * fg_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_x, pb_z, ff_15, ff_16, gd_7, gd_8, gf_13, \
                         gf_14, gf_15, gf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_z[k] * gf_13[k];

        t_25[k] = f_3 * ff_15[k]
                  + f_2 * gd_8[k]
                  + pb_x[k] * gf_15[k];

        t_26[k] = f_2 * gd_7[k]
                  + pb_z[k] * gf_14[k];

        t_27[k] = f_3 * ff_16[k]
                  + pb_x[k] * gf_16[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_x, pa_z, pb_z, dg_15, fg_10, fg_25, \
                         gd_8, gd_9, gf_16, gf_17, gf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_2 * dg_15[k]
                  + pa_x[k] * fg_25[k];

        t_29[k] = pb_z[k] * gf_16[k];

        t_30[k] = f_2 * gd_8[k]
                  + pb_z[k] * gf_17[k];

        t_31[k] = f_1 * gd_9[k]
                  + pb_z[k] * gf_18[k];

        t_32[k] = pa_z[k] * fg_10[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pa_z, pb_y, pb_z, dg_0, ff_9, \
                         fg_14, fg_19, gd_10, gf_19, gf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_y[k] * fg_19[k];

        t_34[k] = f_2 * dg_0[k]
                  + pa_z[k] * fg_14[k];

        t_35[k] = pb_y[k] * gf_19[k];

        t_36[k] = f_3 * ff_9[k]
                  + pb_z[k] * gf_19[k];

        t_37[k] = f_2 * gd_10[k]
                  + pb_y[k] * gf_20[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_y, ff_18, ff_19, gd_11, gd_13, \
                         gf_21, gf_22, gf_23, gf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * gf_21[k];

        t_39[k] = f_3 * ff_18[k]
                  + f_2 * gd_13[k]
                  + pb_x[k] * gf_22[k];

        t_40[k] = f_3 * ff_19[k]
                  + pb_x[k] * gf_26[k];

        t_41[k] = f_1 * gd_11[k]
                  + pb_y[k] * gf_23[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_x, pb_y, dg_31, ff_20, fg_29, fg_30, \
                         gd_12, gd_13, gf_24, gf_25, gf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * gd_12[k]
                  + pb_y[k] * gf_24[k];

        t_43[k] = f_2 * gd_13[k]
                  + pb_y[k] * gf_25[k];

        t_44[k] = pb_y[k] * gf_26[k];

        t_45[k] = f_2 * dg_31[k]
                  + pa_x[k] * fg_29[k];

        t_46[k] = f_0 * ff_20[k]
                  + pa_x[k] * fg_30[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pa_x, pa_z, pb_z, ff_22, fg_20, fg_32, \
                         fg_37, gd_14, gf_27, gf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_z[k] * gf_27[k];

        t_48[k] = f_3 * ff_22[k]
                  + pa_x[k] * fg_32[k];

        t_49[k] = f_2 * gd_14[k]
                  + pb_z[k] * gf_28[k];

        t_50[k] = pa_x[k] * fg_37[k];

        t_51[k] = pa_z[k] * fg_20[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pa_x, pb_y, pb_z, ff_17, ff_35, fg_53, \
                         gd_15, gf_29, gf_30, gf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_0 * ff_35[k]
                  + pa_x[k] * fg_53[k];

        t_53[k] = pb_y[k] * gf_29[k];

        t_54[k] = f_1 * ff_17[k]
                  + pb_z[k] * gf_29[k];

        t_55[k] = f_2 * gd_15[k]
                  + pb_y[k] * gf_30[k];

        t_56[k] = pb_y[k] * gf_31[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pa_x, pb_x, ff_38, fg_56, fg_64, gd_16, \
                         gd_17, gd_18, gf_32, gf_33, gf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_3 * ff_38[k]
                  + pa_x[k] * fg_56[k];

        t_58[k] = pa_x[k] * fg_64[k];

        t_59[k] = f_1 * gd_16[k]
                  + pb_x[k] * gf_32[k];

        t_60[k] = f_3 * gd_17[k]
                  + pb_x[k] * gf_33[k];

        t_61[k] = f_2 * gd_18[k]
                  + pb_x[k] * gf_34[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, t_67, pb_x, pb_y, pb_z, ff_24, gd_18, \
                         gd_19, gf_35, gf_36, gf_38, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_2 * gd_19[k]
                  + pb_x[k] * gf_35[k];

        t_63[k] = pb_x[k] * gf_36[k];

        t_64[k] = pb_x[k] * gf_38[k];

        t_65[k] = pb_x[k] * gf_39[k];

        t_66[k] = f_0 * ff_24[k]
                  + f_1 * gd_18[k]
                  + pb_y[k] * gf_36[k];

        t_67[k] = pb_z[k] * gf_36[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, ff_27, gd_18, gd_19, gd_20, \
                         gf_37, gf_39, gf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_2 * gd_18[k]
                  + pb_z[k] * gf_37[k];

        t_69[k] = f_0 * ff_27[k]
                  + pb_y[k] * gf_39[k];

        t_70[k] = f_1 * gd_19[k]
                  + pb_z[k] * gf_39[k];

        t_71[k] = f_3 * gd_20[k]
                  + pb_x[k] * gf_40[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, pa_z, pb_x, fg_37, gd_22, gd_23, \
                         gf_41, gf_42, gf_44, gf_45, gf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_2 * gd_22[k]
                  + pb_x[k] * gf_41[k];

        t_73[k] = f_2 * gd_23[k]
                  + pb_x[k] * gf_42[k];

        t_74[k] = pb_x[k] * gf_44[k];

        t_75[k] = pb_x[k] * gf_45[k];

        t_76[k] = pb_x[k] * gf_46[k];

        t_77[k] = pa_z[k] * fg_37[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_y, pa_z, pb_y, pb_z, dg_20, ff_24, ff_25, \
                         ff_30, fg_39, fg_46, gf_43, gf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_2 * ff_24[k]
                  + pb_z[k] * gf_43[k];

        t_79[k] = f_3 * ff_25[k]
                  + pa_z[k] * fg_39[k];

        t_80[k] = f_1 * ff_30[k]
                  + pb_y[k] * gf_46[k];

        t_81[k] = f_3 * dg_20[k]
                  + pa_y[k] * fg_46[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pb_x, gd_24, gd_25, gd_26, gd_27, \
                         gd_28, gf_47, gf_48, gf_49, gf_50, gf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_1 * gd_24[k]
                  + pb_x[k] * gf_47[k];

        t_83[k] = f_3 * gd_25[k]
                  + pb_x[k] * gf_48[k];

        t_84[k] = f_3 * gd_26[k]
                  + pb_x[k] * gf_49[k];

        t_85[k] = f_2 * gd_27[k]
                  + pb_x[k] * gf_50[k];

        t_86[k] = f_2 * gd_28[k]
                  + pb_x[k] * gf_51[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, pa_z, pb_x, dg_15, fg_44, gd_29, \
                         gf_52, gf_53, gf_54, gf_55, gf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_2 * gd_29[k]
                  + pb_x[k] * gf_52[k];

        t_88[k] = pb_x[k] * gf_53[k];

        t_89[k] = pb_x[k] * gf_54[k];

        t_90[k] = pb_x[k] * gf_55[k];

        t_91[k] = pb_x[k] * gf_56[k];

        t_92[k] = f_2 * dg_15[k]
                  + pa_z[k] * fg_44[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pa_y, pb_y, pb_z, dg_31, ff_29, ff_33, ff_34, \
                         fg_52, gd_29, gf_53, gf_55, gf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_3 * ff_29[k]
                  + pb_z[k] * gf_53[k];

        t_94[k] = f_3 * ff_33[k]
                  + f_2 * gd_29[k]
                  + pb_y[k] * gf_55[k];

        t_95[k] = f_3 * ff_34[k]
                  + pb_y[k] * gf_56[k];

        t_96[k] = f_2 * dg_31[k]
                  + pa_y[k] * fg_52[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, t_102, pb_x, gd_30, gd_31, gd_32, \
                         gf_57, gf_58, gf_59, gf_60, gf_61, gf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_3 * gd_30[k]
                  + pb_x[k] * gf_57[k];

        t_98[k] = f_2 * gd_31[k]
                  + pb_x[k] * gf_58[k];

        t_99[k] = f_2 * gd_32[k]
                  + pb_x[k] * gf_59[k];

        t_100[k] = pb_x[k] * gf_60[k];

        t_101[k] = pb_x[k] * gf_61[k];

        t_102[k] = pb_x[k] * gf_62[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pa_y, pb_y, pb_z, ff_32, ff_39, ff_41, \
                         ff_42, fg_60, fg_62, gf_60, gf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_0 * ff_39[k]
                   + pa_y[k] * fg_60[k];

        t_104[k] = f_1 * ff_32[k]
                   + pb_z[k] * gf_60[k];

        t_105[k] = f_3 * ff_41[k]
                   + pa_y[k] * fg_62[k];

        t_106[k] = f_2 * ff_42[k]
                   + pb_y[k] * gf_63[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, pa_y, pb_x, fg_64, gd_34, gd_35, \
                         gd_36, gd_38, gf_64, gf_65, gf_66, gf_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pa_y[k] * fg_64[k];

        t_108[k] = f_1 * gd_34[k]
                   + pb_x[k] * gf_64[k];

        t_109[k] = f_3 * gd_35[k]
                   + pb_x[k] * gf_65[k];

        t_110[k] = f_2 * gd_36[k]
                   + pb_x[k] * gf_66[k];

        t_111[k] = f_2 * gd_38[k]
                   + pb_x[k] * gf_67[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, t_117, t_118, pb_x, pb_y, gd_36, \
                         gd_37, gd_38, gf_68, gf_69, gf_70, gf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = pb_x[k] * gf_68[k];

        t_113[k] = pb_x[k] * gf_69[k];

        t_114[k] = pb_x[k] * gf_71[k];

        t_115[k] = f_1 * gd_36[k]
                   + pb_y[k] * gf_68[k];

        t_116[k] = f_3 * gd_37[k]
                   + pb_y[k] * gf_69[k];

        t_117[k] = f_2 * gd_38[k]
                   + pb_y[k] * gf_70[k];

        t_118[k] = pb_y[k] * gf_71[k];
    }

#pragma omp simd aligned(t_119, pb_z, ff_42, gd_38, gf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_0 * ff_42[k]
                   + f_1 * gd_38[k]
                   + pb_z[k] * gf_71[k];
    }
}

auto
compute_prim_gg_overlap_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dg, const size_t ff, const size_t fg,
                          const size_t gd, const size_t gf, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_10 = buffer.data(dg + 10);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);

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
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_40 = buffer.data(gd + 40);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_64 = buffer.data(gf + 64);
    const auto *gf_65 = buffer.data(gf + 65);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_67 = buffer.data(gf + 67);
    const auto *gf_69 = buffer.data(gf + 69);
    const auto *gf_70 = buffer.data(gf + 70);
    const auto *gf_71 = buffer.data(gf + 71);
    const auto *gf_72 = buffer.data(gf + 72);
    const auto *gf_75 = buffer.data(gf + 75);
    const auto *gf_76 = buffer.data(gf + 76);
    const auto *gf_78 = buffer.data(gf + 78);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_80 = buffer.data(gf + 80);
    const auto *gf_81 = buffer.data(gf + 81);
    const auto *gf_82 = buffer.data(gf + 82);
    const auto *gf_83 = buffer.data(gf + 83);
    const auto *gf_84 = buffer.data(gf + 84);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, ff_0, ff_3, gd_0, gf_0, gf_1, \
                         gf_2, gf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = f_2 * gd_0[k]
                 + pb_y[k] * gf_1[k];

        t_2[k] = f_2 * gd_0[k]
                 + pb_z[k] * gf_2[k];

        t_3[k] = f_0 * ff_3[k]
                 + pb_x[k] * gf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pb_x, pb_y, pb_z, ff_4, fg_0, gd_1, gd_2, \
                         gf_3, gf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_0 * ff_4[k]
                 + pb_x[k] * gf_4[k];

        t_5[k] = f_1 * gd_1[k]
                 + pb_y[k] * gf_3[k];

        t_6[k] = f_1 * gd_2[k]
                 + pb_z[k] * gf_4[k];

        t_7[k] = pa_y[k] * fg_0[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_y, pb_x, pb_y, dg_2, ff_0, ff_1, ff_6, \
                         fg_1, fg_4, gf_5, gf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * ff_0[k]
                 + pb_y[k] * gf_5[k];

        t_9[k] = f_3 * ff_1[k]
                 + pa_y[k] * fg_1[k];

        t_10[k] = f_1 * ff_6[k]
                  + pb_x[k] * gf_6[k];

        t_11[k] = f_3 * dg_2[k]
                  + pa_x[k] * fg_4[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_z, ff_0, ff_2, ff_9, fg_0, \
                         fg_2, gf_9, gf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_z[k] * fg_0[k];

        t_13[k] = f_2 * ff_0[k]
                  + pb_z[k] * gf_9[k];

        t_14[k] = f_3 * ff_2[k]
                  + pa_z[k] * fg_2[k];

        t_15[k] = f_1 * ff_9[k]
                  + pb_x[k] * gf_12[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, pb_y, dg_0, dg_3, ff_5, fg_3, fg_7, \
                         gf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * dg_3[k]
                  + pa_x[k] * fg_7[k];

        t_17[k] = f_2 * dg_0[k]
                  + pa_y[k] * fg_3[k];

        t_18[k] = f_3 * ff_5[k]
                  + pb_y[k] * gf_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pa_y, pb_x, dg_4, ff_11, ff_12, fg_6, \
                         fg_8, gd_10, gf_14, gf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * ff_11[k]
                  + f_2 * gd_10[k]
                  + pb_x[k] * gf_14[k];

        t_20[k] = f_3 * ff_12[k]
                  + pb_x[k] * gf_15[k];

        t_21[k] = f_2 * dg_4[k]
                  + pa_x[k] * fg_8[k];

        t_22[k] = pa_y[k] * fg_6[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pa_z, pb_y, pb_z, dg_0, dg_6, ff_7, \
                         fg_5, fg_9, gd_14, gf_23, gf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_2 * dg_6[k]
                  + pa_x[k] * fg_9[k];

        t_24[k] = f_2 * dg_0[k]
                  + pa_z[k] * fg_5[k];

        t_25[k] = f_3 * ff_7[k]
                  + pb_z[k] * gf_23[k];

        t_26[k] = f_2 * gd_14[k]
                  + pb_y[k] * gf_24[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pb_x, dg_10, ff_15, ff_16, ff_17, \
                         fg_10, fg_11, gd_17, gf_26, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_3 * ff_15[k]
                  + f_2 * gd_17[k]
                  + pb_x[k] * gf_26[k];

        t_28[k] = f_3 * ff_16[k]
                  + pb_x[k] * gf_29[k];

        t_29[k] = f_2 * dg_10[k]
                  + pa_x[k] * fg_10[k];

        t_30[k] = f_0 * ff_17[k]
                  + pa_x[k] * fg_11[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_x, pb_x, pb_y, ff_10, ff_18, ff_19, \
                         fg_12, fg_13, fg_16, gf_30, gf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_1 * ff_10[k]
                  + pb_y[k] * gf_30[k];

        t_32[k] = f_3 * ff_18[k]
                  + pa_x[k] * fg_12[k];

        t_33[k] = f_2 * ff_19[k]
                  + pb_x[k] * gf_32[k];

        t_34[k] = pa_x[k] * fg_13[k];

        t_35[k] = pa_x[k] * fg_16[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pa_x, pb_z, ff_13, ff_31, fg_17, fg_18, \
                         fg_19, fg_21, gf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_x[k] * fg_17[k];

        t_37[k] = pa_x[k] * fg_18[k];

        t_38[k] = pa_x[k] * fg_19[k];

        t_39[k] = f_0 * ff_31[k]
                  + pa_x[k] * fg_21[k];

        t_40[k] = f_1 * ff_13[k]
                  + pb_z[k] * gf_45[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_x, pb_x, ff_34, ff_38, fg_23, fg_26, \
                         gd_25, gd_26, gf_50, gf_51, gf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_3 * ff_34[k]
                  + pa_x[k] * fg_23[k];

        t_42[k] = f_2 * ff_38[k]
                  + pb_x[k] * gf_50[k];

        t_43[k] = pa_x[k] * fg_26[k];

        t_44[k] = f_1 * gd_25[k]
                  + pb_x[k] * gf_51[k];

        t_45[k] = f_3 * gd_26[k]
                  + pb_x[k] * gf_52[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pb_x, pb_y, pb_z, ff_19, gd_27, gd_28, \
                         gf_52, gf_53, gf_54, gf_55, gf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_2 * gd_27[k]
                  + pb_x[k] * gf_53[k];

        t_47[k] = pb_z[k] * gf_52[k];

        t_48[k] = f_2 * gd_28[k]
                  + pb_x[k] * gf_54[k];

        t_49[k] = f_0 * ff_19[k]
                  + f_1 * gd_27[k]
                  + pb_y[k] * gf_55[k];

        t_50[k] = f_2 * gd_27[k]
                  + pb_z[k] * gf_56[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_x, pb_y, pb_z, ff_22, fg_13, gd_28, \
                         gd_30, gf_58, gf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * ff_22[k]
                  + pb_y[k] * gf_58[k];

        t_52[k] = f_1 * gd_28[k]
                  + pb_z[k] * gf_58[k];

        t_53[k] = f_2 * gd_30[k]
                  + pb_x[k] * gf_59[k];

        t_54[k] = pa_z[k] * fg_13[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pa_z, pb_y, pb_z, dg_7, ff_19, ff_20, \
                         ff_26, fg_14, fg_17, gf_60, gf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_2 * ff_19[k]
                  + pb_z[k] * gf_60[k];

        t_56[k] = f_3 * ff_20[k]
                  + pa_z[k] * fg_14[k];

        t_57[k] = f_1 * ff_26[k]
                  + pb_y[k] * gf_63[k];

        t_58[k] = f_3 * dg_7[k]
                  + pa_y[k] * fg_17[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_z, pb_x, dg_4, fg_15, gd_31, gd_32, gd_33, \
                         gf_64, gf_65, gf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * gd_31[k]
                  + pb_x[k] * gf_64[k];

        t_60[k] = f_2 * gd_32[k]
                  + pb_x[k] * gf_65[k];

        t_61[k] = f_2 * gd_33[k]
                  + pb_x[k] * gf_66[k];

        t_62[k] = f_2 * dg_4[k]
                  + pa_z[k] * fg_15[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_y, pb_y, pb_z, dg_10, ff_23, ff_29, ff_30, \
                         fg_20, gd_33, gf_67, gf_69, gf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_3 * ff_23[k]
                  + pb_z[k] * gf_67[k];

        t_64[k] = f_3 * ff_29[k]
                  + f_2 * gd_33[k]
                  + pb_y[k] * gf_69[k];

        t_65[k] = f_3 * ff_30[k]
                  + pb_y[k] * gf_70[k];

        t_66[k] = f_2 * dg_10[k]
                  + pa_y[k] * fg_20[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_y, pb_x, pb_z, ff_27, ff_35, ff_37, fg_24, \
                         fg_25, gd_34, gf_71, gf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_2 * gd_34[k]
                  + pb_x[k] * gf_71[k];

        t_68[k] = f_0 * ff_35[k]
                  + pa_y[k] * fg_24[k];

        t_69[k] = f_1 * ff_27[k]
                  + pb_z[k] * gf_72[k];

        t_70[k] = f_3 * ff_37[k]
                  + pa_y[k] * fg_25[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, pa_y, pb_x, pb_y, ff_38, fg_26, gd_36, \
                         gd_37, gf_75, gf_76, gf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_2 * ff_38[k]
                  + pb_y[k] * gf_75[k];

        t_72[k] = pa_y[k] * fg_26[k];

        t_73[k] = f_1 * gd_36[k]
                  + pb_x[k] * gf_76[k];

        t_74[k] = pb_y[k] * gf_76[k];

        t_75[k] = f_3 * gd_37[k]
                  + pb_x[k] * gf_78[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pb_x, pb_y, gd_38, gd_39, gd_40, gf_78, \
                         gf_79, gf_80, gf_81, gf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_2 * gd_38[k]
                  + pb_x[k] * gf_79[k];

        t_77[k] = pb_y[k] * gf_78[k];

        t_78[k] = f_2 * gd_40[k]
                  + pb_x[k] * gf_80[k];

        t_79[k] = f_1 * gd_38[k]
                  + pb_y[k] * gf_81[k];

        t_80[k] = f_3 * gd_39[k]
                  + pb_y[k] * gf_82[k];
    }

#pragma omp simd aligned(t_81, t_82, pb_y, pb_z, ff_38, gd_40, gf_83, \
                         gf_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_2 * gd_40[k]
                  + pb_y[k] * gf_83[k];

        t_82[k] = f_0 * ff_38[k]
                  + f_1 * gd_40[k]
                  + pb_z[k] * gf_84[k];
    }
}

auto
compute_prim_gg_overlap_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dg, const size_t ff, const size_t fg,
                          const size_t gd, const size_t gf, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_19 = buffer.data(dg + 19);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_4 = buffer.data(ff + 4);
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
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
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
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_28 = buffer.data(fg + 28);
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
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_55 = buffer.data(fg + 55);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);

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
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_36 = buffer.data(gf + 36);
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
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_62 = buffer.data(gf + 62);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_64 = buffer.data(gf + 64);
    const auto *gf_65 = buffer.data(gf + 65);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_67 = buffer.data(gf + 67);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, ff_0, gd_0, gd_1, \
                         gf_0, gf_1, gf_2, gf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_2 * gd_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_2 * gd_0[k]
                 + pb_z[k] * gf_2[k];

        t_5[k] = f_1 * gd_1[k]
                 + pb_y[k] * gf_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_y, pb_y, pb_z, ff_1, fg_0, fg_3, fg_4, \
                         gd_2, gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * gd_2[k]
                 + pb_y[k] * gf_4[k];

        t_7[k] = f_1 * gd_2[k]
                 + pb_z[k] * gf_5[k];

        t_8[k] = pa_y[k] * fg_0[k];

        t_9[k] = f_3 * ff_1[k]
                 + pa_y[k] * fg_3[k];

        t_10[k] = pa_y[k] * fg_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pa_x, pa_y, pb_y, pb_z, dg_4, ff_4, fg_6, \
                         fg_9, gd_4, gf_8, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * dg_4[k]
                  + pa_x[k] * fg_9[k];

        t_12[k] = f_2 * gd_4[k]
                  + pb_z[k] * gf_8[k];

        t_13[k] = f_2 * ff_4[k]
                  + pb_y[k] * gf_9[k];

        t_14[k] = pa_y[k] * fg_6[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_z, pb_y, pb_z, ff_0, ff_2, fg_0, fg_4, \
                         gd_6, gf_10, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_z[k] * fg_0[k];

        t_16[k] = f_2 * ff_0[k]
                  + pb_z[k] * gf_10[k];

        t_17[k] = f_3 * ff_2[k]
                  + pa_z[k] * fg_4[k];

        t_18[k] = f_3 * gd_6[k]
                  + pb_y[k] * gf_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pa_y, pb_y, dg_0, dg_5, fg_7, fg_13, gd_7, \
                         gf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * gd_7[k]
                  + pb_y[k] * gf_12[k];

        t_20[k] = f_3 * dg_5[k]
                  + pa_x[k] * fg_13[k];

        t_21[k] = f_2 * dg_0[k]
                  + pa_y[k] * fg_7[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, dg_8, ff_11, ff_12, fg_17, \
                         gd_8, gd_9, gf_15, gf_16, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * ff_11[k]
                  + f_2 * gd_9[k]
                  + pb_x[k] * gf_16[k];

        t_23[k] = f_2 * gd_8[k]
                  + pb_z[k] * gf_15[k];

        t_24[k] = f_3 * ff_12[k]
                  + pb_x[k] * gf_17[k];

        t_25[k] = f_2 * dg_8[k]
                  + pa_x[k] * fg_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, pb_z, ff_7, fg_8, \
                         fg_11, gd_9, gd_10, gf_18, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * gd_9[k]
                  + pb_z[k] * gf_18[k];

        t_27[k] = f_3 * ff_7[k]
                  + pb_y[k] * gf_19[k];

        t_28[k] = f_1 * gd_10[k]
                  + pb_z[k] * gf_19[k];

        t_29[k] = pa_y[k] * fg_11[k];

        t_30[k] = pa_z[k] * fg_8[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pa_y, pa_z, pb_z, dg_11, ff_6, fg_9, \
                         fg_12, fg_18, gf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pa_y[k] * fg_12[k];

        t_32[k] = pa_z[k] * fg_9[k];

        t_33[k] = f_2 * ff_6[k]
                  + pb_z[k] * gf_20[k];

        t_34[k] = f_2 * dg_11[k]
                  + pa_x[k] * fg_18[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_y, pa_z, pb_y, pb_z, dg_0, ff_8, \
                         ff_9, fg_10, fg_13, gf_21, gf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_2 * ff_9[k]
                  + pb_y[k] * gf_21[k];

        t_36[k] = pa_y[k] * fg_13[k];

        t_37[k] = f_2 * dg_0[k]
                  + pa_z[k] * fg_10[k];

        t_38[k] = pb_y[k] * gf_22[k];

        t_39[k] = f_3 * ff_8[k]
                  + pb_z[k] * gf_22[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_x, pb_y, ff_14, ff_15, gd_11, gd_12, \
                         gd_14, gf_23, gf_24, gf_25, gf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_2 * gd_11[k]
                  + pb_y[k] * gf_23[k];

        t_41[k] = f_3 * ff_14[k]
                  + f_2 * gd_14[k]
                  + pb_x[k] * gf_24[k];

        t_42[k] = f_3 * ff_15[k]
                  + pb_x[k] * gf_28[k];

        t_43[k] = f_1 * gd_12[k]
                  + pb_y[k] * gf_25[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_y, dg_19, ff_16, fg_23, fg_24, \
                         gd_13, gd_14, gf_26, gf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_3 * gd_13[k]
                  + pb_y[k] * gf_26[k];

        t_45[k] = f_2 * gd_14[k]
                  + pb_y[k] * gf_27[k];

        t_46[k] = f_2 * dg_19[k]
                  + pa_x[k] * fg_23[k];

        t_47[k] = f_0 * ff_16[k]
                  + pa_x[k] * fg_24[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_x, pb_x, pb_z, ff_17, ff_19, fg_25, \
                         fg_28, fg_30, gd_15, gf_30, gf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * ff_17[k]
                  + pa_x[k] * fg_25[k];

        t_49[k] = f_2 * gd_15[k]
                  + pb_z[k] * gf_30[k];

        t_50[k] = f_2 * ff_19[k]
                  + pb_x[k] * gf_32[k];

        t_51[k] = pa_x[k] * fg_28[k];

        t_52[k] = pa_x[k] * fg_30[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_x, pa_z, pb_z, ff_10, fg_14, fg_15, \
                         fg_31, fg_32, gf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pa_x[k] * fg_31[k];

        t_54[k] = pa_x[k] * fg_32[k];

        t_55[k] = pa_z[k] * fg_14[k];

        t_56[k] = f_2 * ff_10[k]
                  + pb_z[k] * gf_33[k];

        t_57[k] = pa_z[k] * fg_15[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, t_63, pa_x, pa_y, ff_22, fg_19, fg_33, \
                         fg_35, fg_36, fg_37, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_3 * ff_22[k]
                  + pa_x[k] * fg_33[k];

        t_59[k] = pa_x[k] * fg_35[k];

        t_60[k] = pa_x[k] * fg_36[k];

        t_61[k] = pa_x[k] * fg_37[k];

        t_62[k] = pa_x[k] * fg_38[k];

        t_63[k] = pa_y[k] * fg_19[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, t_69, pa_x, pa_y, ff_25, fg_20, fg_21, \
                         fg_39, fg_40, fg_41, fg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pa_y[k] * fg_20[k];

        t_65[k] = f_3 * ff_25[k]
                  + pa_x[k] * fg_39[k];

        t_66[k] = pa_y[k] * fg_21[k];

        t_67[k] = pa_x[k] * fg_40[k];

        t_68[k] = pa_x[k] * fg_41[k];

        t_69[k] = pa_x[k] * fg_42[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_x, pb_z, ff_13, ff_29, ff_32, fg_43, \
                         fg_45, fg_48, gf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_x[k] * fg_43[k];

        t_71[k] = f_0 * ff_29[k]
                  + pa_x[k] * fg_45[k];

        t_72[k] = f_1 * ff_13[k]
                  + pb_z[k] * gf_36[k];

        t_73[k] = f_3 * ff_32[k]
                  + pa_x[k] * fg_48[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, pa_x, pb_x, ff_36, fg_51, fg_52, \
                         fg_53, fg_55, gd_19, gf_38, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_2 * ff_36[k]
                  + pb_x[k] * gf_38[k];

        t_75[k] = pa_x[k] * fg_51[k];

        t_76[k] = pa_x[k] * fg_52[k];

        t_77[k] = pa_x[k] * fg_53[k];

        t_78[k] = pa_x[k] * fg_55[k];

        t_79[k] = f_1 * gd_19[k]
                  + pb_x[k] * gf_39[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, t_85, pb_x, pb_z, gd_20, gd_21, gd_22, \
                         gf_40, gf_41, gf_42, gf_43, gf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_3 * gd_20[k]
                  + pb_x[k] * gf_40[k];

        t_81[k] = f_2 * gd_21[k]
                  + pb_x[k] * gf_41[k];

        t_82[k] = pb_z[k] * gf_40[k];

        t_83[k] = f_2 * gd_22[k]
                  + pb_x[k] * gf_42[k];

        t_84[k] = pb_x[k] * gf_43[k];

        t_85[k] = pb_x[k] * gf_45[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, t_91, pb_x, pb_y, pb_z, ff_19, ff_21, \
                         gd_21, gd_22, gf_43, gf_44, gf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_x[k] * gf_46[k];

        t_87[k] = f_0 * ff_19[k]
                  + f_1 * gd_21[k]
                  + pb_y[k] * gf_43[k];

        t_88[k] = pb_z[k] * gf_43[k];

        t_89[k] = f_2 * gd_21[k]
                  + pb_z[k] * gf_44[k];

        t_90[k] = f_0 * ff_21[k]
                  + pb_y[k] * gf_46[k];

        t_91[k] = f_1 * gd_22[k]
                  + pb_z[k] * gf_46[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, pa_z, pb_x, pb_z, ff_19, ff_20, fg_28, \
                         fg_30, gd_24, gf_47, gf_48, gf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_2 * gd_24[k]
                  + pb_x[k] * gf_47[k];

        t_93[k] = pb_x[k] * gf_49[k];

        t_94[k] = pa_z[k] * fg_28[k];

        t_95[k] = f_2 * ff_19[k]
                  + pb_z[k] * gf_48[k];

        t_96[k] = f_3 * ff_20[k]
                  + pa_z[k] * fg_30[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_y, pb_x, pb_y, dg_12, ff_24, fg_38, \
                         gd_25, gd_26, gf_49, gf_50, gf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_1 * ff_24[k]
                  + pb_y[k] * gf_49[k];

        t_98[k] = f_3 * dg_12[k]
                  + pa_y[k] * fg_38[k];

        t_99[k] = f_1 * gd_25[k]
                  + pb_x[k] * gf_50[k];

        t_100[k] = f_2 * gd_26[k]
                   + pb_x[k] * gf_51[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, pa_z, pb_x, pb_z, dg_8, ff_23, \
                         fg_34, gd_27, gf_52, gf_53, gf_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_2 * gd_27[k]
                   + pb_x[k] * gf_52[k];

        t_102[k] = pb_x[k] * gf_53[k];

        t_103[k] = pb_x[k] * gf_55[k];

        t_104[k] = f_2 * dg_8[k]
                   + pa_z[k] * fg_34[k];

        t_105[k] = f_3 * ff_23[k]
                   + pb_z[k] * gf_53[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pa_y, pb_x, pb_y, dg_19, ff_27, ff_28, \
                         fg_44, gd_27, gd_28, gf_54, gf_55, gf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_3 * ff_27[k]
                   + f_2 * gd_27[k]
                   + pb_y[k] * gf_54[k];

        t_107[k] = f_3 * ff_28[k]
                   + pb_y[k] * gf_55[k];

        t_108[k] = f_2 * dg_19[k]
                   + pa_y[k] * fg_44[k];

        t_109[k] = f_2 * gd_28[k]
                   + pb_x[k] * gf_56[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_y, pb_x, pb_z, ff_26, ff_33, ff_35, \
                         fg_51, fg_53, gf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pb_x[k] * gf_57[k];

        t_111[k] = f_0 * ff_33[k]
                   + pa_y[k] * fg_51[k];

        t_112[k] = f_1 * ff_26[k]
                   + pb_z[k] * gf_57[k];

        t_113[k] = f_3 * ff_35[k]
                   + pa_y[k] * fg_53[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, pa_y, pb_x, pb_y, ff_36, fg_55, \
                         gd_30, gd_31, gf_59, gf_60, gf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_2 * ff_36[k]
                   + pb_y[k] * gf_59[k];

        t_115[k] = pa_y[k] * fg_55[k];

        t_116[k] = f_1 * gd_30[k]
                   + pb_x[k] * gf_60[k];

        t_117[k] = pb_y[k] * gf_60[k];

        t_118[k] = f_3 * gd_31[k]
                   + pb_x[k] * gf_61[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, t_124, pb_x, pb_y, gd_32, gd_34, \
                         gf_61, gf_62, gf_63, gf_64, gf_65, gf_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_2 * gd_32[k]
                   + pb_x[k] * gf_62[k];

        t_120[k] = pb_y[k] * gf_61[k];

        t_121[k] = f_2 * gd_34[k]
                   + pb_x[k] * gf_63[k];

        t_122[k] = pb_x[k] * gf_64[k];

        t_123[k] = pb_x[k] * gf_65[k];

        t_124[k] = pb_x[k] * gf_67[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, pb_y, pb_z, ff_36, gd_32, gd_33, \
                         gd_34, gf_64, gf_65, gf_66, gf_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_1 * gd_32[k]
                   + pb_y[k] * gf_64[k];

        t_126[k] = f_3 * gd_33[k]
                   + pb_y[k] * gf_65[k];

        t_127[k] = f_2 * gd_34[k]
                   + pb_y[k] * gf_66[k];

        t_128[k] = pb_y[k] * gf_67[k];

        t_129[k] = f_0 * ff_36[k]
                   + f_1 * gd_34[k]
                   + pb_z[k] * gf_67[k];
    }
}

auto
compute_prim_gg_overlap_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dg, const size_t ff, const size_t fg,
                          const size_t gd, const size_t gf, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_27 = buffer.data(dg + 27);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_49 = buffer.data(fg + 49);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
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
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_29 = buffer.data(gf + 29);
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
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, ff_0, gd_0, gd_1, \
                         gf_0, gf_1, gf_2, gf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_2 * gd_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_2 * gd_0[k]
                 + pb_z[k] * gf_2[k];

        t_5[k] = f_1 * gd_1[k]
                 + pb_y[k] * gf_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pb_y, pb_z, dg_6, fg_0, fg_8, \
                         gd_2, gf_4, gf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * gd_2[k]
                 + pb_y[k] * gf_4[k];

        t_7[k] = pb_y[k] * gf_5[k];

        t_8[k] = f_1 * gd_2[k]
                 + pb_z[k] * gf_5[k];

        t_9[k] = pa_y[k] * fg_0[k];

        t_10[k] = f_3 * dg_6[k]
                  + pa_x[k] * fg_8[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pa_z, pb_z, ff_2, fg_0, fg_4, \
                         fg_6, gd_4, gf_7, gf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * gf_7[k];

        t_12[k] = f_2 * gd_4[k]
                  + pb_z[k] * gf_8[k];

        t_13[k] = pa_y[k] * fg_6[k];

        t_14[k] = pa_z[k] * fg_0[k];

        t_15[k] = f_3 * ff_2[k]
                  + pa_z[k] * fg_4[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, dg_8, fg_11, gd_6, gd_7, gf_10, \
                         gf_11, gf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * gd_6[k]
                  + pb_y[k] * gf_10[k];

        t_17[k] = f_2 * gd_7[k]
                  + pb_y[k] * gf_11[k];

        t_18[k] = pb_y[k] * gf_12[k];

        t_19[k] = f_3 * dg_8[k]
                  + pa_x[k] * fg_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_y, pb_x, pb_z, dg_0, ff_9, fg_7, gd_8, \
                         gd_9, gf_13, gf_14, gf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_2 * dg_0[k]
                  + pa_y[k] * fg_7[k];

        t_21[k] = pb_z[k] * gf_13[k];

        t_22[k] = f_3 * ff_9[k]
                  + f_2 * gd_9[k]
                  + pb_x[k] * gf_15[k];

        t_23[k] = f_2 * gd_8[k]
                  + pb_z[k] * gf_14[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_x, pb_x, pb_z, dg_12, ff_10, fg_15, \
                         gd_9, gd_10, gf_16, gf_17, gf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_3 * ff_10[k]
                  + pb_x[k] * gf_16[k];

        t_25[k] = f_2 * dg_12[k]
                  + pa_x[k] * fg_15[k];

        t_26[k] = pb_z[k] * gf_16[k];

        t_27[k] = f_2 * gd_9[k]
                  + pb_z[k] * gf_17[k];

        t_28[k] = f_1 * gd_10[k]
                  + pb_z[k] * gf_18[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pa_y, pa_z, pb_y, dg_0, fg_8, fg_10, \
                         fg_11, gd_11, gf_19, gf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_z[k] * fg_8[k];

        t_30[k] = pa_y[k] * fg_11[k];

        t_31[k] = f_2 * dg_0[k]
                  + pa_z[k] * fg_10[k];

        t_32[k] = pb_y[k] * gf_19[k];

        t_33[k] = f_2 * gd_11[k]
                  + pb_y[k] * gf_20[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_x, pb_y, ff_11, ff_12, gd_12, gd_13, \
                         gd_14, gf_21, gf_22, gf_23, gf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * ff_11[k]
                  + f_2 * gd_14[k]
                  + pb_x[k] * gf_21[k];

        t_35[k] = f_3 * ff_12[k]
                  + pb_x[k] * gf_25[k];

        t_36[k] = f_1 * gd_12[k]
                  + pb_y[k] * gf_22[k];

        t_37[k] = f_3 * gd_13[k]
                  + pb_y[k] * gf_23[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_x, pb_y, pb_z, dg_27, ff_13, fg_19, \
                         fg_20, gd_14, gf_24, gf_25, gf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_2 * gd_14[k]
                  + pb_y[k] * gf_24[k];

        t_39[k] = pb_y[k] * gf_25[k];

        t_40[k] = f_2 * dg_27[k]
                  + pa_x[k] * fg_19[k];

        t_41[k] = f_0 * ff_13[k]
                  + pa_x[k] * fg_20[k];

        t_42[k] = pb_z[k] * gf_26[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pb_x, pb_z, ff_14, ff_16, fg_21, fg_25, \
                         gd_15, gf_27, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * ff_14[k]
                  + pa_x[k] * fg_21[k];

        t_44[k] = f_2 * gd_15[k]
                  + pb_z[k] * gf_27[k];

        t_45[k] = f_2 * ff_16[k]
                  + pb_x[k] * gf_29[k];

        t_46[k] = pa_x[k] * fg_25[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pa_x, pa_z, ff_26, ff_29, fg_12, fg_32, \
                         fg_35, fg_39, fg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pa_z[k] * fg_12[k];

        t_48[k] = pa_x[k] * fg_32[k];

        t_49[k] = pa_x[k] * fg_35[k];

        t_50[k] = f_0 * ff_26[k]
                  + pa_x[k] * fg_39[k];

        t_51[k] = f_3 * ff_29[k]
                  + pa_x[k] * fg_42[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pa_x, pb_x, ff_33, fg_49, gd_19, gd_20, \
                         gd_21, gf_32, gf_33, gf_34, gf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_2 * ff_33[k]
                  + pb_x[k] * gf_32[k];

        t_53[k] = pa_x[k] * fg_49[k];

        t_54[k] = f_1 * gd_19[k]
                  + pb_x[k] * gf_33[k];

        t_55[k] = f_3 * gd_20[k]
                  + pb_x[k] * gf_34[k];

        t_56[k] = f_2 * gd_21[k]
                  + pb_x[k] * gf_35[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, t_62, pb_x, pb_y, pb_z, ff_16, gd_21, \
                         gd_22, gf_36, gf_37, gf_39, gf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_2 * gd_22[k]
                  + pb_x[k] * gf_36[k];

        t_58[k] = pb_x[k] * gf_37[k];

        t_59[k] = pb_x[k] * gf_39[k];

        t_60[k] = pb_x[k] * gf_40[k];

        t_61[k] = f_0 * ff_16[k]
                  + f_1 * gd_21[k]
                  + pb_y[k] * gf_37[k];

        t_62[k] = pb_z[k] * gf_37[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pb_x, pb_y, pb_z, ff_18, gd_21, gd_22, \
                         gd_24, gf_38, gf_40, gf_41, gf_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_2 * gd_21[k]
                  + pb_z[k] * gf_38[k];

        t_64[k] = f_0 * ff_18[k]
                  + pb_y[k] * gf_40[k];

        t_65[k] = f_1 * gd_22[k]
                  + pb_z[k] * gf_40[k];

        t_66[k] = f_2 * gd_24[k]
                  + pb_x[k] * gf_41[k];

        t_67[k] = pb_x[k] * gf_43[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_y, pa_z, pb_x, dg_17, ff_17, fg_25, fg_27, \
                         fg_32, gd_25, gf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pa_z[k] * fg_25[k];

        t_69[k] = f_3 * ff_17[k]
                  + pa_z[k] * fg_27[k];

        t_70[k] = f_3 * dg_17[k]
                  + pa_y[k] * fg_32[k];

        t_71[k] = f_1 * gd_25[k]
                  + pb_x[k] * gf_44[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pa_z, pb_x, dg_12, fg_31, gd_26, gd_27, \
                         gf_45, gf_46, gf_47, gf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_2 * gd_26[k]
                  + pb_x[k] * gf_45[k];

        t_73[k] = f_2 * gd_27[k]
                  + pb_x[k] * gf_46[k];

        t_74[k] = pb_x[k] * gf_47[k];

        t_75[k] = pb_x[k] * gf_49[k];

        t_76[k] = f_2 * dg_12[k]
                  + pa_z[k] * fg_31[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pb_x, pb_y, dg_27, ff_24, ff_25, fg_38, \
                         gd_27, gd_28, gf_48, gf_49, gf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_3 * ff_24[k]
                  + f_2 * gd_27[k]
                  + pb_y[k] * gf_48[k];

        t_78[k] = f_3 * ff_25[k]
                  + pb_y[k] * gf_49[k];

        t_79[k] = f_2 * dg_27[k]
                  + pa_y[k] * fg_38[k];

        t_80[k] = f_2 * gd_28[k]
                  + pb_x[k] * gf_50[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, pa_y, pb_x, pb_y, ff_30, ff_32, ff_33, \
                         fg_45, fg_47, fg_49, gf_51, gf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pb_x[k] * gf_51[k];

        t_82[k] = f_0 * ff_30[k]
                  + pa_y[k] * fg_45[k];

        t_83[k] = f_3 * ff_32[k]
                  + pa_y[k] * fg_47[k];

        t_84[k] = f_2 * ff_33[k]
                  + pb_y[k] * gf_53[k];

        t_85[k] = pa_y[k] * fg_49[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, pb_x, gd_30, gd_31, gd_32, gd_34, \
                         gf_54, gf_55, gf_56, gf_57, gf_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_1 * gd_30[k]
                  + pb_x[k] * gf_54[k];

        t_87[k] = f_3 * gd_31[k]
                  + pb_x[k] * gf_55[k];

        t_88[k] = f_2 * gd_32[k]
                  + pb_x[k] * gf_56[k];

        t_89[k] = f_2 * gd_34[k]
                  + pb_x[k] * gf_57[k];

        t_90[k] = pb_x[k] * gf_58[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, t_96, pb_x, pb_y, gd_32, gd_33, gd_34, \
                         gf_58, gf_59, gf_60, gf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = pb_x[k] * gf_59[k];

        t_92[k] = pb_x[k] * gf_61[k];

        t_93[k] = f_1 * gd_32[k]
                  + pb_y[k] * gf_58[k];

        t_94[k] = f_3 * gd_33[k]
                  + pb_y[k] * gf_59[k];

        t_95[k] = f_2 * gd_34[k]
                  + pb_y[k] * gf_60[k];

        t_96[k] = pb_y[k] * gf_61[k];
    }

#pragma omp simd aligned(t_97, pb_z, ff_33, gd_34, gf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_0 * ff_33[k]
                  + f_1 * gd_34[k]
                  + pb_z[k] * gf_61[k];
    }
}

auto
compute_prim_gg_overlap_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dg, const size_t ff, const size_t fg,
                          const size_t gd, const size_t gf, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_10 = buffer.data(dg + 10);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);

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

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_35 = buffer.data(gd + 35);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, ff_0, fg_0, gd_0, gf_0, \
                         gf_1, gf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = f_2 * gd_0[k]
                 + pb_y[k] * gf_1[k];

        t_2[k] = f_2 * gd_0[k]
                 + pb_z[k] * gf_2[k];

        t_3[k] = pa_y[k] * fg_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_z, dg_2, dg_3, ff_1, fg_0, fg_1, fg_3, \
                         fg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dg_2[k]
                 + pa_x[k] * fg_3[k];

        t_5[k] = pa_z[k] * fg_0[k];

        t_6[k] = f_3 * ff_1[k]
                 + pa_z[k] * fg_1[k];

        t_7[k] = f_3 * dg_3[k]
                 + pa_x[k] * fg_6[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_y, pb_x, dg_0, dg_4, ff_6, fg_2, fg_5, \
                         fg_7, gd_9, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * dg_0[k]
                 + pa_y[k] * fg_2[k];

        t_9[k] = f_3 * ff_6[k]
                 + f_2 * gd_9[k]
                 + pb_x[k] * gf_11[k];

        t_10[k] = f_2 * dg_4[k]
                  + pa_x[k] * fg_7[k];

        t_11[k] = pa_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pa_z, pb_x, dg_0, dg_6, dg_10, ff_8, \
                         fg_4, fg_8, fg_9, gd_14, gf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_2 * dg_6[k]
                  + pa_x[k] * fg_8[k];

        t_13[k] = f_2 * dg_0[k]
                  + pa_z[k] * fg_4[k];

        t_14[k] = f_3 * ff_8[k]
                  + f_2 * gd_14[k]
                  + pb_x[k] * gf_15[k];

        t_15[k] = f_2 * dg_10[k]
                  + pa_x[k] * fg_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, pa_x, fg_10, fg_13, fg_14, fg_15, \
                         fg_16, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_x[k] * fg_10[k];

        t_17[k] = pa_x[k] * fg_13[k];

        t_18[k] = pa_x[k] * fg_14[k];

        t_19[k] = pa_x[k] * fg_15[k];

        t_20[k] = pa_x[k] * fg_16[k];

        t_21[k] = pa_x[k] * fg_20[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pb_x, pb_y, pb_z, ff_11, gd_22, gd_23, gf_26, \
                         gf_27, gf_28, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_1 * gd_22[k]
                  + pb_x[k] * gf_26[k];

        t_23[k] = f_2 * gd_23[k]
                  + pb_x[k] * gf_27[k];

        t_24[k] = f_0 * ff_11[k]
                  + f_1 * gd_23[k]
                  + pb_y[k] * gf_28[k];

        t_25[k] = f_2 * gd_23[k]
                  + pb_z[k] * gf_29[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, dg_4, dg_7, ff_12, fg_10, fg_11, \
                         fg_12, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_z[k] * fg_10[k];

        t_27[k] = f_3 * ff_12[k]
                  + pa_z[k] * fg_11[k];

        t_28[k] = f_3 * dg_7[k]
                  + pa_y[k] * fg_14[k];

        t_29[k] = f_2 * dg_4[k]
                  + pa_z[k] * fg_12[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_y, dg_10, ff_16, ff_19, ff_20, \
                         fg_17, fg_18, fg_19, gd_29, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * ff_16[k]
                  + f_2 * gd_29[k]
                  + pb_y[k] * gf_39[k];

        t_31[k] = f_2 * dg_10[k]
                  + pa_y[k] * fg_17[k];

        t_32[k] = f_0 * ff_19[k]
                  + pa_y[k] * fg_18[k];

        t_33[k] = f_3 * ff_20[k]
                  + pa_y[k] * fg_19[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_y, pb_x, pb_y, fg_20, gd_32, gd_33, \
                         gd_35, gf_45, gf_47, gf_48, gf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_y[k] * fg_20[k];

        t_35[k] = f_1 * gd_32[k]
                  + pb_x[k] * gf_45[k];

        t_36[k] = f_2 * gd_33[k]
                  + pb_x[k] * gf_47[k];

        t_37[k] = f_2 * gd_35[k]
                  + pb_x[k] * gf_48[k];

        t_38[k] = f_1 * gd_33[k]
                  + pb_y[k] * gf_49[k];
    }

#pragma omp simd aligned(t_39, t_40, pb_y, pb_z, ff_21, gd_35, gf_51, \
                         gf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_2 * gd_35[k]
                  + pb_y[k] * gf_51[k];

        t_40[k] = f_0 * ff_21[k]
                  + f_1 * gd_35[k]
                  + pb_z[k] * gf_52[k];
    }
}

auto
compute_prim_gg_overlap_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dg, const size_t ff, const size_t fg,
                          const size_t gd, const size_t gf, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_15 = buffer.data(dg + 15);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_6 = buffer.data(ff + 6);
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
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_34 = buffer.data(fg + 34);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_26 = buffer.data(gf + 26);
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
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, ff_0, gd_0, gd_1, \
                         gf_0, gf_1, gf_2, gf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_2 * gd_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_2 * gd_0[k]
                 + pb_z[k] * gf_2[k];

        t_5[k] = f_1 * gd_1[k]
                 + pb_y[k] * gf_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, pb_z, dg_3, ff_1, fg_0, \
                         fg_2, fg_5, gd_2, gf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * gd_2[k]
                 + pb_z[k] * gf_4[k];

        t_7[k] = pa_y[k] * fg_0[k];

        t_8[k] = f_3 * ff_1[k]
                 + pa_y[k] * fg_2[k];

        t_9[k] = f_3 * dg_3[k]
                 + pa_x[k] * fg_5[k];

        t_10[k] = pa_z[k] * fg_0[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pa_x, pa_y, pa_z, pb_z, dg_0, dg_4, ff_0, \
                         ff_2, fg_3, fg_4, fg_8, gf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * ff_0[k]
                  + pb_z[k] * gf_8[k];

        t_12[k] = f_3 * ff_2[k]
                  + pa_z[k] * fg_3[k];

        t_13[k] = f_3 * dg_4[k]
                  + pa_x[k] * fg_8[k];

        t_14[k] = f_2 * dg_0[k]
                  + pa_y[k] * fg_4[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pa_y, pb_x, dg_6, ff_9, ff_10, fg_7, \
                         fg_11, gd_8, gf_11, gf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_3 * ff_9[k]
                  + f_2 * gd_8[k]
                  + pb_x[k] * gf_11[k];

        t_16[k] = f_3 * ff_10[k]
                  + pb_x[k] * gf_12[k];

        t_17[k] = f_2 * dg_6[k]
                  + pa_x[k] * fg_11[k];

        t_18[k] = pa_y[k] * fg_7[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pa_z, pb_y, pb_z, dg_0, dg_9, ff_6, \
                         fg_6, fg_12, gd_10, gf_16, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * dg_9[k]
                  + pa_x[k] * fg_12[k];

        t_20[k] = f_2 * dg_0[k]
                  + pa_z[k] * fg_6[k];

        t_21[k] = f_3 * ff_6[k]
                  + pb_z[k] * gf_16[k];

        t_22[k] = f_2 * gd_10[k]
                  + pb_y[k] * gf_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_x, dg_15, ff_12, ff_13, ff_14, \
                         fg_15, fg_16, gd_11, gf_18, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * ff_12[k]
                  + f_2 * gd_11[k]
                  + pb_x[k] * gf_18[k];

        t_24[k] = f_3 * ff_13[k]
                  + pb_x[k] * gf_19[k];

        t_25[k] = f_2 * dg_15[k]
                  + pa_x[k] * fg_15[k];

        t_26[k] = f_0 * ff_14[k]
                  + pa_x[k] * fg_16[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, t_32, pa_x, pb_x, ff_15, ff_16, fg_17, \
                         fg_18, fg_22, fg_23, fg_24, gf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_3 * ff_15[k]
                  + pa_x[k] * fg_17[k];

        t_28[k] = f_2 * ff_16[k]
                  + pb_x[k] * gf_22[k];

        t_29[k] = pa_x[k] * fg_18[k];

        t_30[k] = pa_x[k] * fg_22[k];

        t_31[k] = pa_x[k] * fg_23[k];

        t_32[k] = pa_x[k] * fg_24[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_x, pb_z, ff_11, ff_24, ff_27, fg_25, \
                         fg_28, fg_30, gf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_x[k] * fg_25[k];

        t_34[k] = f_0 * ff_24[k]
                  + pa_x[k] * fg_28[k];

        t_35[k] = f_1 * ff_11[k]
                  + pb_z[k] * gf_26[k];

        t_36[k] = f_3 * ff_27[k]
                  + pa_x[k] * fg_30[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pa_x, pb_x, ff_30, fg_34, gd_14, gd_15, \
                         gd_16, gf_28, gf_29, gf_30, gf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_2 * ff_30[k]
                  + pb_x[k] * gf_28[k];

        t_38[k] = pa_x[k] * fg_34[k];

        t_39[k] = f_1 * gd_14[k]
                  + pb_x[k] * gf_29[k];

        t_40[k] = f_2 * gd_15[k]
                  + pb_x[k] * gf_30[k];

        t_41[k] = f_2 * gd_16[k]
                  + pb_x[k] * gf_31[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pb_x, pb_y, pb_z, ff_16, ff_18, \
                         gd_15, gd_16, gf_32, gf_33, gf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_x[k] * gf_32[k];

        t_43[k] = f_0 * ff_16[k]
                  + f_1 * gd_15[k]
                  + pb_y[k] * gf_32[k];

        t_44[k] = pb_z[k] * gf_32[k];

        t_45[k] = f_2 * gd_15[k]
                  + pb_z[k] * gf_33[k];

        t_46[k] = f_0 * ff_18[k]
                  + pb_y[k] * gf_34[k];

        t_47[k] = f_1 * gd_16[k]
                  + pb_z[k] * gf_34[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_z, ff_16, ff_17, fg_18, fg_20, \
                         gd_18, gf_35, gf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_2 * gd_18[k]
                  + pb_x[k] * gf_35[k];

        t_49[k] = pa_z[k] * fg_18[k];

        t_50[k] = f_2 * ff_16[k]
                  + pb_z[k] * gf_36[k];

        t_51[k] = f_3 * ff_17[k]
                  + pa_z[k] * fg_20[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_y, pb_x, pb_y, dg_10, ff_20, fg_23, gd_19, \
                         gd_20, gf_37, gf_38, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * ff_20[k]
                  + pb_y[k] * gf_37[k];

        t_53[k] = f_3 * dg_10[k]
                  + pa_y[k] * fg_23[k];

        t_54[k] = f_1 * gd_19[k]
                  + pb_x[k] * gf_38[k];

        t_55[k] = f_2 * gd_20[k]
                  + pb_x[k] * gf_39[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_z, pb_x, pb_y, pb_z, dg_6, ff_19, ff_22, \
                         fg_21, gd_21, gf_40, gf_41, gf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_2 * gd_21[k]
                  + pb_x[k] * gf_40[k];

        t_57[k] = f_2 * dg_6[k]
                  + pa_z[k] * fg_21[k];

        t_58[k] = f_3 * ff_19[k]
                  + pb_z[k] * gf_41[k];

        t_59[k] = f_3 * ff_22[k]
                  + f_2 * gd_21[k]
                  + pb_y[k] * gf_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pb_x, pb_y, dg_15, ff_23, ff_28, fg_27, \
                         fg_31, gd_22, gf_43, gf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_3 * ff_23[k]
                  + pb_y[k] * gf_43[k];

        t_61[k] = f_2 * dg_15[k]
                  + pa_y[k] * fg_27[k];

        t_62[k] = f_2 * gd_22[k]
                  + pb_x[k] * gf_44[k];

        t_63[k] = f_0 * ff_28[k]
                  + pa_y[k] * fg_31[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_y, pb_z, ff_21, ff_29, ff_30, fg_32, \
                         fg_34, gf_45, gf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_1 * ff_21[k]
                  + pb_z[k] * gf_45[k];

        t_65[k] = f_3 * ff_29[k]
                  + pa_y[k] * fg_32[k];

        t_66[k] = f_2 * ff_30[k]
                  + pb_y[k] * gf_47[k];

        t_67[k] = pa_y[k] * fg_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, pb_x, pb_y, gd_24, gd_25, gd_27, \
                         gf_48, gf_50, gf_51, gf_52, gf_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_1 * gd_24[k]
                  + pb_x[k] * gf_48[k];

        t_69[k] = pb_y[k] * gf_48[k];

        t_70[k] = f_2 * gd_25[k]
                  + pb_x[k] * gf_50[k];

        t_71[k] = f_2 * gd_27[k]
                  + pb_x[k] * gf_51[k];

        t_72[k] = pb_x[k] * gf_52[k];

        t_73[k] = pb_x[k] * gf_55[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, pb_y, pb_z, ff_30, gd_25, gd_26, gd_27, \
                         gf_52, gf_53, gf_54, gf_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_1 * gd_25[k]
                  + pb_y[k] * gf_52[k];

        t_75[k] = f_3 * gd_26[k]
                  + pb_y[k] * gf_53[k];

        t_76[k] = f_2 * gd_27[k]
                  + pb_y[k] * gf_54[k];

        t_77[k] = pb_y[k] * gf_55[k];

        t_78[k] = f_0 * ff_30[k]
                  + f_1 * gd_27[k]
                  + pb_z[k] * gf_55[k];
    }
}

auto
compute_prim_gg_overlap_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dg, const size_t ff, const size_t fg,
                          const size_t gd, const size_t gf, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_18 = buffer.data(dg + 18);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_40 = buffer.data(fg + 40);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
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
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, ff_0, gd_0, gd_1, \
                         gf_0, gf_1, gf_2, gf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 + f_1 * gd_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = pb_y[k] * gf_0[k];

        t_2[k] = pb_z[k] * gf_0[k];

        t_3[k] = f_2 * gd_0[k]
                 + pb_y[k] * gf_1[k];

        t_4[k] = f_2 * gd_0[k]
                 + pb_z[k] * gf_2[k];

        t_5[k] = f_1 * gd_1[k]
                 + pb_y[k] * gf_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, pb_z, dg_3, fg_0, fg_5, \
                         fg_7, gd_2, gf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * gd_2[k]
                 + pb_z[k] * gf_4[k];

        t_7[k] = pa_y[k] * fg_0[k];

        t_8[k] = f_3 * dg_3[k]
                 + pa_x[k] * fg_7[k];

        t_9[k] = pa_y[k] * fg_5[k];

        t_10[k] = pa_z[k] * fg_0[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pa_x, pa_y, pb_x, dg_0, dg_4, ff_8, ff_9, \
                         fg_6, fg_9, gd_8, gf_10, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * dg_4[k]
                  + pa_x[k] * fg_9[k];

        t_12[k] = f_2 * dg_0[k]
                  + pa_y[k] * fg_6[k];

        t_13[k] = f_3 * ff_8[k]
                  + f_2 * gd_8[k]
                  + pb_x[k] * gf_10[k];

        t_14[k] = f_3 * ff_9[k]
                  + pb_x[k] * gf_11[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pa_y, pa_z, pb_z, dg_8, fg_7, fg_9, \
                         fg_13, gd_9, gf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_2 * dg_8[k]
                  + pa_x[k] * fg_13[k];

        t_16[k] = f_1 * gd_9[k]
                  + pb_z[k] * gf_12[k];

        t_17[k] = pa_z[k] * fg_7[k];

        t_18[k] = pa_y[k] * fg_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_z, pb_x, pb_y, dg_0, ff_10, fg_8, gd_10, \
                         gd_11, gf_13, gf_14, gf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * dg_0[k]
                  + pa_z[k] * fg_8[k];

        t_20[k] = pb_y[k] * gf_13[k];

        t_21[k] = f_2 * gd_10[k]
                  + pb_y[k] * gf_14[k];

        t_22[k] = f_3 * ff_10[k]
                  + f_2 * gd_11[k]
                  + pb_x[k] * gf_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_x, dg_18, ff_11, ff_12, ff_13, \
                         fg_17, fg_18, fg_19, gf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * ff_11[k]
                  + pb_x[k] * gf_16[k];

        t_24[k] = f_2 * dg_18[k]
                  + pa_x[k] * fg_17[k];

        t_25[k] = f_0 * ff_12[k]
                  + pa_x[k] * fg_18[k];

        t_26[k] = f_3 * ff_13[k]
                  + pa_x[k] * fg_19[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_x, pa_z, pb_x, ff_14, fg_10, fg_21, \
                         fg_26, fg_27, gf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_2 * ff_14[k]
                  + pb_x[k] * gf_18[k];

        t_28[k] = pa_x[k] * fg_21[k];

        t_29[k] = pa_z[k] * fg_10[k];

        t_30[k] = pa_x[k] * fg_26[k];

        t_31[k] = pa_x[k] * fg_27[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, pa_x, pb_x, ff_21, ff_24, ff_27, fg_31, \
                         fg_34, fg_40, gd_14, gf_20, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * ff_21[k]
                  + pa_x[k] * fg_31[k];

        t_33[k] = f_3 * ff_24[k]
                  + pa_x[k] * fg_34[k];

        t_34[k] = f_2 * ff_27[k]
                  + pb_x[k] * gf_20[k];

        t_35[k] = pa_x[k] * fg_40[k];

        t_36[k] = f_1 * gd_14[k]
                  + pb_x[k] * gf_21[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, t_42, pb_x, pb_y, pb_z, ff_14, gd_15, \
                         gd_16, gf_22, gf_23, gf_24, gf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_2 * gd_15[k]
                  + pb_x[k] * gf_22[k];

        t_38[k] = f_2 * gd_16[k]
                  + pb_x[k] * gf_23[k];

        t_39[k] = pb_x[k] * gf_24[k];

        t_40[k] = pb_x[k] * gf_26[k];

        t_41[k] = f_0 * ff_14[k]
                  + f_1 * gd_15[k]
                  + pb_y[k] * gf_24[k];

        t_42[k] = pb_z[k] * gf_24[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pa_z, pb_x, pb_z, fg_21, gd_15, gd_16, \
                         gd_18, gf_25, gf_26, gf_27, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_2 * gd_15[k]
                  + pb_z[k] * gf_25[k];

        t_44[k] = f_1 * gd_16[k]
                  + pb_z[k] * gf_26[k];

        t_45[k] = f_2 * gd_18[k]
                  + pb_x[k] * gf_27[k];

        t_46[k] = pb_x[k] * gf_29[k];

        t_47[k] = pa_z[k] * fg_21[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_y, pb_x, dg_11, fg_26, gd_19, gd_20, \
                         gd_21, gf_30, gf_31, gf_32, gf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * dg_11[k]
                  + pa_y[k] * fg_26[k];

        t_49[k] = f_1 * gd_19[k]
                  + pb_x[k] * gf_30[k];

        t_50[k] = f_2 * gd_20[k]
                  + pb_x[k] * gf_31[k];

        t_51[k] = f_2 * gd_21[k]
                  + pb_x[k] * gf_32[k];

        t_52[k] = pb_x[k] * gf_33[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_z, pb_x, pb_y, dg_8, ff_19, ff_20, fg_25, \
                         gd_21, gf_34, gf_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_x[k] * gf_35[k];

        t_54[k] = f_2 * dg_8[k]
                  + pa_z[k] * fg_25[k];

        t_55[k] = f_3 * ff_19[k]
                  + f_2 * gd_21[k]
                  + pb_y[k] * gf_34[k];

        t_56[k] = f_3 * ff_20[k]
                  + pb_y[k] * gf_35[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pa_y, pb_x, dg_18, ff_25, ff_26, fg_30, \
                         fg_37, fg_38, gd_22, gf_36, gf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_2 * dg_18[k]
                  + pa_y[k] * fg_30[k];

        t_58[k] = f_2 * gd_22[k]
                  + pb_x[k] * gf_36[k];

        t_59[k] = pb_x[k] * gf_37[k];

        t_60[k] = f_0 * ff_25[k]
                  + pa_y[k] * fg_37[k];

        t_61[k] = f_3 * ff_26[k]
                  + pa_y[k] * fg_38[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, pa_y, pb_x, pb_y, ff_27, fg_40, gd_24, \
                         gd_25, gf_39, gf_40, gf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_2 * ff_27[k]
                  + pb_y[k] * gf_39[k];

        t_63[k] = pa_y[k] * fg_40[k];

        t_64[k] = f_1 * gd_24[k]
                  + pb_x[k] * gf_40[k];

        t_65[k] = pb_y[k] * gf_40[k];

        t_66[k] = f_2 * gd_25[k]
                  + pb_x[k] * gf_42[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, t_72, pb_x, pb_y, gd_25, gd_26, gd_27, \
                         gf_43, gf_44, gf_45, gf_46, gf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_2 * gd_27[k]
                  + pb_x[k] * gf_43[k];

        t_68[k] = pb_x[k] * gf_44[k];

        t_69[k] = pb_x[k] * gf_47[k];

        t_70[k] = f_1 * gd_25[k]
                  + pb_y[k] * gf_44[k];

        t_71[k] = f_3 * gd_26[k]
                  + pb_y[k] * gf_45[k];

        t_72[k] = f_2 * gd_27[k]
                  + pb_y[k] * gf_46[k];
    }

#pragma omp simd aligned(t_73, t_74, pb_y, pb_z, ff_27, gd_27, gf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pb_y[k] * gf_47[k];

        t_74[k] = f_0 * ff_27[k]
                  + f_1 * gd_27[k]
                  + pb_z[k] * gf_47[k];
    }
}

}  // namespace simdovl
