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
    const auto *dg0_55 = buffer.data(dg0 + 55);
    const auto *dg0_89 = buffer.data(dg0 + 89);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_55 = buffer.data(dg1 + 55);
    const auto *dg1_89 = buffer.data(dg1 + 89);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_42 = buffer.data(ff + 42);
    const auto *ff_47 = buffer.data(ff + 47);
    const auto *ff_48 = buffer.data(ff + 48);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_55 = buffer.data(ff + 55);
    const auto *ff_56 = buffer.data(ff + 56);
    const auto *ff_57 = buffer.data(ff + 57);
    const auto *ff_59 = buffer.data(ff + 59);
    const auto *ff_60 = buffer.data(ff + 60);
    const auto *ff_62 = buffer.data(ff + 62);
    const auto *ff_63 = buffer.data(ff + 63);
    const auto *ff_65 = buffer.data(ff + 65);
    const auto *ff_66 = buffer.data(ff + 66);
    const auto *ff_67 = buffer.data(ff + 67);
    const auto *ff_68 = buffer.data(ff + 68);
    const auto *ff_69 = buffer.data(ff + 69);
    const auto *ff_70 = buffer.data(ff + 70);
    const auto *ff_72 = buffer.data(ff + 72);
    const auto *ff_75 = buffer.data(ff + 75);
    const auto *ff_76 = buffer.data(ff + 76);
    const auto *ff_77 = buffer.data(ff + 77);
    const auto *ff_78 = buffer.data(ff + 78);
    const auto *ff_79 = buffer.data(ff + 79);
    const auto *ff_80 = buffer.data(ff + 80);
    const auto *ff_82 = buffer.data(ff + 82);
    const auto *ff_83 = buffer.data(ff + 83);
    const auto *ff_86 = buffer.data(ff + 86);
    const auto *ff_87 = buffer.data(ff + 87);
    const auto *ff_88 = buffer.data(ff + 88);
    const auto *ff_89 = buffer.data(ff + 89);
    const auto *ff_90 = buffer.data(ff + 90);
    const auto *ff_91 = buffer.data(ff + 91);
    const auto *ff_92 = buffer.data(ff + 92);
    const auto *ff_93 = buffer.data(ff + 93);
    const auto *ff_95 = buffer.data(ff + 95);
    const auto *ff_96 = buffer.data(ff + 96);
    const auto *ff_97 = buffer.data(ff + 97);
    const auto *ff_98 = buffer.data(ff + 98);
    const auto *ff_99 = buffer.data(ff + 99);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_84 = buffer.data(fg + 84);
    const auto *fg_89 = buffer.data(fg + 89);
    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_91 = buffer.data(fg + 91);
    const auto *fg_93 = buffer.data(fg + 93);
    const auto *fg_95 = buffer.data(fg + 95);
    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_102 = buffer.data(fg + 102);
    const auto *fg_103 = buffer.data(fg + 103);
    const auto *fg_104 = buffer.data(fg + 104);
    const auto *fg_110 = buffer.data(fg + 110);
    const auto *fg_115 = buffer.data(fg + 115);
    const auto *fg_116 = buffer.data(fg + 116);
    const auto *fg_117 = buffer.data(fg + 117);
    const auto *fg_118 = buffer.data(fg + 118);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_123 = buffer.data(fg + 123);
    const auto *fg_130 = buffer.data(fg + 130);
    const auto *fg_131 = buffer.data(fg + 131);
    const auto *fg_132 = buffer.data(fg + 132);
    const auto *fg_133 = buffer.data(fg + 133);
    const auto *fg_134 = buffer.data(fg + 134);
    const auto *fg_135 = buffer.data(fg + 135);
    const auto *fg_137 = buffer.data(fg + 137);
    const auto *fg_138 = buffer.data(fg + 138);
    const auto *fg_140 = buffer.data(fg + 140);
    const auto *fg_145 = buffer.data(fg + 145);
    const auto *fg_146 = buffer.data(fg + 146);
    const auto *fg_147 = buffer.data(fg + 147);
    const auto *fg_149 = buffer.data(fg + 149);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_21 = buffer.data(gd0 + 21);
    const auto *gd0_23 = buffer.data(gd0 + 23);
    const auto *gd0_30 = buffer.data(gd0 + 30);
    const auto *gd0_33 = buffer.data(gd0 + 33);
    const auto *gd0_35 = buffer.data(gd0 + 35);
    const auto *gd0_60 = buffer.data(gd0 + 60);
    const auto *gd0_63 = buffer.data(gd0 + 63);
    const auto *gd0_65 = buffer.data(gd0 + 65);
    const auto *gd0_72 = buffer.data(gd0 + 72);
    const auto *gd0_75 = buffer.data(gd0 + 75);
    const auto *gd0_77 = buffer.data(gd0 + 77);
    const auto *gd0_84 = buffer.data(gd0 + 84);
    const auto *gd0_87 = buffer.data(gd0 + 87);
    const auto *gd0_89 = buffer.data(gd0 + 89);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_21 = buffer.data(gd1 + 21);
    const auto *gd1_23 = buffer.data(gd1 + 23);
    const auto *gd1_30 = buffer.data(gd1 + 30);
    const auto *gd1_33 = buffer.data(gd1 + 33);
    const auto *gd1_35 = buffer.data(gd1 + 35);
    const auto *gd1_60 = buffer.data(gd1 + 60);
    const auto *gd1_63 = buffer.data(gd1 + 63);
    const auto *gd1_65 = buffer.data(gd1 + 65);
    const auto *gd1_72 = buffer.data(gd1 + 72);
    const auto *gd1_75 = buffer.data(gd1 + 75);
    const auto *gd1_77 = buffer.data(gd1 + 77);
    const auto *gd1_84 = buffer.data(gd1 + 84);
    const auto *gd1_87 = buffer.data(gd1 + 87);
    const auto *gd1_89 = buffer.data(gd1 + 89);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_68 = buffer.data(gf + 68);
    const auto *gf_69 = buffer.data(gf + 69);
    const auto *gf_70 = buffer.data(gf + 70);
    const auto *gf_72 = buffer.data(gf + 72);
    const auto *gf_77 = buffer.data(gf + 77);
    const auto *gf_78 = buffer.data(gf + 78);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_80 = buffer.data(gf + 80);
    const auto *gf_82 = buffer.data(gf + 82);
    const auto *gf_86 = buffer.data(gf + 86);
    const auto *gf_87 = buffer.data(gf + 87);
    const auto *gf_88 = buffer.data(gf + 88);
    const auto *gf_90 = buffer.data(gf + 90);
    const auto *gf_92 = buffer.data(gf + 92);
    const auto *gf_95 = buffer.data(gf + 95);
    const auto *gf_96 = buffer.data(gf + 96);
    const auto *gf_97 = buffer.data(gf + 97);
    const auto *gf_99 = buffer.data(gf + 99);
    const auto *gf_100 = buffer.data(gf + 100);
    const auto *gf_101 = buffer.data(gf + 101);
    const auto *gf_103 = buffer.data(gf + 103);
    const auto *gf_105 = buffer.data(gf + 105);
    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_107 = buffer.data(gf + 107);
    const auto *gf_108 = buffer.data(gf + 108);
    const auto *gf_109 = buffer.data(gf + 109);
    const auto *gf_110 = buffer.data(gf + 110);
    const auto *gf_112 = buffer.data(gf + 112);
    const auto *gf_116 = buffer.data(gf + 116);
    const auto *gf_117 = buffer.data(gf + 117);
    const auto *gf_118 = buffer.data(gf + 118);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_120 = buffer.data(gf + 120);
    const auto *gf_122 = buffer.data(gf + 122);
    const auto *gf_123 = buffer.data(gf + 123);
    const auto *gf_125 = buffer.data(gf + 125);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_127 = buffer.data(gf + 127);
    const auto *gf_128 = buffer.data(gf + 128);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_130 = buffer.data(gf + 130);
    const auto *gf_132 = buffer.data(gf + 132);
    const auto *gf_136 = buffer.data(gf + 136);
    const auto *gf_137 = buffer.data(gf + 137);
    const auto *gf_138 = buffer.data(gf + 138);
    const auto *gf_139 = buffer.data(gf + 139);
    const auto *gf_140 = buffer.data(gf + 140);
    const auto *gf_142 = buffer.data(gf + 142);
    const auto *gf_143 = buffer.data(gf + 143);
    const auto *gf_145 = buffer.data(gf + 145);
    const auto *gf_146 = buffer.data(gf + 146);
    const auto *gf_147 = buffer.data(gf + 147);
    const auto *gf_148 = buffer.data(gf + 148);
    const auto *gf_149 = buffer.data(gf + 149);

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

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, ff_6, ff_9, gd0_3, gd1_3, \
                         gf_3, gf_5, gf_6, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * ff_6[k]
                 + pb_x[k] * gf_6[k];

        t_7[k] = pb_z[k] * gf_3[k];

        t_8[k] = pb_y[k] * gf_5[k];

        t_9[k] = f_0 * ff_9[k]
                 + pb_x[k] * gf_9[k];

        t_10[k] = f_1 * gd0_3[k]
                  - f_2 * gd1_3[k]
                  + pb_y[k] * gf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, fg_0, gd0_5, gd1_5, \
                         gf_6, gf_8, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * gf_6[k];

        t_12[k] = f_3 * gd0_5[k]
                  - f_4 * gd1_5[k]
                  + pb_y[k] * gf_8[k];

        t_13[k] = pb_y[k] * gf_9[k];

        t_14[k] = f_1 * gd0_5[k]
                  - f_2 * gd1_5[k]
                  + pb_z[k] * gf_9[k];

        t_15[k] = pa_y[k] * fg_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pb_y, pb_z, ff_0, ff_1, fg_3, \
                         fg_5, gf_10, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * ff_0[k]
                  + pb_y[k] * gf_10[k];

        t_17[k] = pb_z[k] * gf_10[k];

        t_18[k] = f_6 * ff_1[k]
                  + pa_y[k] * fg_3[k];

        t_19[k] = pb_z[k] * gf_11[k];

        t_20[k] = pa_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pb_x, pb_z, ff_6, ff_16, ff_18, \
                         fg_9, fg_10, gf_13, gf_16, gf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * ff_16[k]
                  + pb_x[k] * gf_16[k];

        t_22[k] = pb_z[k] * gf_13[k];

        t_23[k] = f_7 * ff_18[k]
                  + pb_x[k] * gf_18[k];

        t_24[k] = pa_y[k] * fg_9[k];

        t_25[k] = f_0 * ff_6[k]
                  + pa_y[k] * fg_10[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, pb_z, ff_8, ff_9, \
                         fg_0, fg_12, fg_14, gf_16, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * gf_16[k];

        t_27[k] = f_6 * ff_8[k]
                  + pa_y[k] * fg_12[k];

        t_28[k] = f_5 * ff_9[k]
                  + pb_y[k] * gf_19[k];

        t_29[k] = pa_y[k] * fg_14[k];

        t_30[k] = pa_z[k] * fg_0[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, ff_0, ff_2, \
                         fg_3, fg_5, fg_6, gf_20, gf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * gf_20[k];

        t_32[k] = f_5 * ff_0[k]
                  + pb_z[k] * gf_20[k];

        t_33[k] = pa_z[k] * fg_3[k];

        t_34[k] = pb_y[k] * gf_22[k];

        t_35[k] = f_6 * ff_2[k]
                  + pa_z[k] * fg_5[k];

        t_36[k] = pa_z[k] * fg_6[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, pb_y, ff_27, ff_29, fg_10, gf_25, \
                         gf_27, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_7 * ff_27[k]
                  + pb_x[k] * gf_27[k];

        t_38[k] = pb_y[k] * gf_25[k];

        t_39[k] = f_7 * ff_29[k]
                  + pb_x[k] * gf_29[k];

        t_40[k] = pa_z[k] * fg_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_z, pb_y, pb_z, ff_6, ff_7, ff_9, fg_12, \
                         fg_14, gf_26, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * ff_6[k]
                  + pb_z[k] * gf_26[k];

        t_42[k] = f_6 * ff_7[k]
                  + pa_z[k] * fg_12[k];

        t_43[k] = pb_y[k] * gf_29[k];

        t_44[k] = f_0 * ff_9[k]
                  + pa_z[k] * fg_14[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_y, pb_z, dg0_0, dg1_0, ff_10, fg_15, \
                         gf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_8 * dg0_0[k]
                  - f_9 * dg1_0[k]
                  + pa_y[k] * fg_15[k];

        t_46[k] = f_6 * ff_10[k]
                  + pb_y[k] * gf_30[k];

        t_47[k] = pb_z[k] * gf_30[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_x, pb_z, ff_33, ff_36, gd0_18, gd0_21, \
                         gd1_18, gd1_21, gf_31, gf_32, gf_33, gf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_6 * ff_33[k]
                  + f_3 * gd0_21[k]
                  - f_4 * gd1_21[k]
                  + pb_x[k] * gf_33[k];

        t_49[k] = pb_z[k] * gf_31[k];

        t_50[k] = f_3 * gd0_18[k]
                  - f_4 * gd1_18[k]
                  + pb_z[k] * gf_32[k];

        t_51[k] = f_6 * ff_36[k]
                  + pb_x[k] * gf_36[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_z, dg0_55, dg1_55, ff_38, \
                         ff_39, fg_55, gf_33, gf_38, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_z[k] * gf_33[k];

        t_53[k] = f_6 * ff_38[k]
                  + pb_x[k] * gf_38[k];

        t_54[k] = f_6 * ff_39[k]
                  + pb_x[k] * gf_39[k];

        t_55[k] = f_8 * dg0_55[k]
                  - f_9 * dg1_55[k]
                  + pa_x[k] * fg_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pb_z, ff_19, gd0_21, gd0_23, gd1_21, \
                         gd1_23, gf_36, gf_37, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * gf_36[k];

        t_57[k] = f_3 * gd0_21[k]
                  - f_4 * gd1_21[k]
                  + pb_z[k] * gf_37[k];

        t_58[k] = f_6 * ff_19[k]
                  + pb_y[k] * gf_39[k];

        t_59[k] = f_1 * gd0_23[k]
                  - f_2 * gd1_23[k]
                  + pb_z[k] * gf_39[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, ff_22, fg_16, \
                         fg_18, fg_30, fg_32, fg_35, gf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * fg_30[k];

        t_61[k] = pa_z[k] * fg_16[k];

        t_62[k] = pa_y[k] * fg_32[k];

        t_63[k] = pa_z[k] * fg_18[k];

        t_64[k] = f_5 * ff_22[k]
                  + pb_y[k] * gf_42[k];

        t_65[k] = pa_y[k] * fg_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, ff_47, ff_48, fg_21, \
                         fg_25, fg_39, gf_47, gf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * fg_21[k];

        t_67[k] = f_6 * ff_47[k]
                  + pb_x[k] * gf_47[k];

        t_68[k] = f_6 * ff_48[k]
                  + pb_x[k] * gf_48[k];

        t_69[k] = pa_y[k] * fg_39[k];

        t_70[k] = pa_z[k] * fg_25[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, pb_z, ff_16, ff_28, ff_29, fg_42, \
                         fg_44, gf_46, gf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * ff_16[k]
                  + pb_z[k] * gf_46[k];

        t_72[k] = f_6 * ff_28[k]
                  + pa_y[k] * fg_42[k];

        t_73[k] = f_5 * ff_29[k]
                  + pb_y[k] * gf_49[k];

        t_74[k] = pa_y[k] * fg_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, dg0_0, dg1_0, ff_20, fg_30, \
                         gd0_30, gd1_30, gf_50, gf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_8 * dg0_0[k]
                  - f_9 * dg1_0[k]
                  + pa_z[k] * fg_30[k];

        t_76[k] = pb_y[k] * gf_50[k];

        t_77[k] = f_6 * ff_20[k]
                  + pb_z[k] * gf_50[k];

        t_78[k] = f_3 * gd0_30[k]
                  - f_4 * gd1_30[k]
                  + pb_y[k] * gf_51[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pb_x, pb_y, ff_55, ff_56, ff_57, \
                         gd0_35, gd1_35, gf_52, gf_55, gf_56, gf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_y[k] * gf_52[k];

        t_80[k] = f_6 * ff_55[k]
                  + f_3 * gd0_35[k]
                  - f_4 * gd1_35[k]
                  + pb_x[k] * gf_55[k];

        t_81[k] = f_6 * ff_56[k]
                  + pb_x[k] * gf_56[k];

        t_82[k] = f_6 * ff_57[k]
                  + pb_x[k] * gf_57[k];

        t_83[k] = pb_y[k] * gf_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_x, pb_y, pb_z, ff_26, ff_59, gd0_33, \
                         gd0_35, gd1_33, gd1_35, gf_56, gf_58, gf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_6 * ff_59[k]
                  + pb_x[k] * gf_59[k];

        t_85[k] = f_1 * gd0_33[k]
                  - f_2 * gd1_33[k]
                  + pb_y[k] * gf_56[k];

        t_86[k] = f_6 * ff_26[k]
                  + pb_z[k] * gf_56[k];

        t_87[k] = f_3 * gd0_35[k]
                  - f_4 * gd1_35[k]
                  + pb_y[k] * gf_58[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_x, pb_y, pb_z, dg0_89, dg1_89, \
                         ff_30, ff_60, fg_89, fg_90, gf_59, gf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pb_y[k] * gf_59[k];

        t_89[k] = f_8 * dg0_89[k]
                  - f_9 * dg1_89[k]
                  + pa_x[k] * fg_89[k];

        t_90[k] = f_0 * ff_60[k]
                  + pa_x[k] * fg_90[k];

        t_91[k] = f_7 * ff_30[k]
                  + pb_y[k] * gf_60[k];

        t_92[k] = pb_z[k] * gf_60[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pa_x, pb_x, pb_z, ff_63, ff_65, ff_66, \
                         fg_93, fg_95, gf_61, gf_63, gf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_6 * ff_63[k]
                  + pa_x[k] * fg_93[k];

        t_94[k] = pb_z[k] * gf_61[k];

        t_95[k] = f_6 * ff_65[k]
                  + pa_x[k] * fg_95[k];

        t_96[k] = f_5 * ff_66[k]
                  + pb_x[k] * gf_66[k];

        t_97[k] = pb_z[k] * gf_63[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, pa_x, pb_x, pb_z, ff_68, ff_69, \
                         fg_100, fg_102, gf_66, gf_68, gf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_5 * ff_68[k]
                  + pb_x[k] * gf_68[k];

        t_99[k] = f_5 * ff_69[k]
                  + pb_x[k] * gf_69[k];

        t_100[k] = pa_x[k] * fg_100[k];

        t_101[k] = pb_z[k] * gf_66[k];

        t_102[k] = pa_x[k] * fg_102[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, t_108, pa_x, pa_z, pb_z, ff_30, \
                         fg_45, fg_46, fg_48, fg_103, fg_104, gf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = pa_x[k] * fg_103[k];

        t_104[k] = pa_x[k] * fg_104[k];

        t_105[k] = pa_z[k] * fg_45[k];

        t_106[k] = pa_z[k] * fg_46[k];

        t_107[k] = f_5 * ff_30[k]
                   + pb_z[k] * gf_70[k];

        t_108[k] = pa_z[k] * fg_48[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_x, pa_z, pb_x, pb_y, ff_42, ff_75, \
                         ff_77, fg_51, fg_110, gf_72, gf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_6 * ff_42[k]
                   + pb_y[k] * gf_72[k];

        t_110[k] = f_6 * ff_75[k]
                   + pa_x[k] * fg_110[k];

        t_111[k] = pa_z[k] * fg_51[k];

        t_112[k] = f_5 * ff_77[k]
                   + pb_x[k] * gf_77[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, t_118, pa_x, pb_x, ff_78, ff_79, \
                         fg_115, fg_116, fg_117, fg_118, gf_78, gf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_5 * ff_78[k]
                   + pb_x[k] * gf_78[k];

        t_114[k] = f_5 * ff_79[k]
                   + pb_x[k] * gf_79[k];

        t_115[k] = pa_x[k] * fg_115[k];

        t_116[k] = pa_x[k] * fg_116[k];

        t_117[k] = pa_x[k] * fg_117[k];

        t_118[k] = pa_x[k] * fg_118[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, pa_x, pa_y, pb_y, ff_50, ff_83, \
                         fg_75, fg_77, fg_119, fg_123, gf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pa_x[k] * fg_119[k];

        t_120[k] = pa_y[k] * fg_75[k];

        t_121[k] = f_5 * ff_50[k]
                   + pb_y[k] * gf_80[k];

        t_122[k] = pa_y[k] * fg_77[k];

        t_123[k] = f_6 * ff_83[k]
                   + pa_x[k] * fg_123[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pb_x, pb_y, ff_52, ff_86, ff_87, \
                         fg_80, gf_82, gf_86, gf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_5 * ff_52[k]
                   + pb_y[k] * gf_82[k];

        t_125[k] = pa_y[k] * fg_80[k];

        t_126[k] = f_5 * ff_86[k]
                   + pb_x[k] * gf_86[k];

        t_127[k] = f_5 * ff_87[k]
                   + pb_x[k] * gf_87[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, pa_x, pa_y, pb_x, ff_88, \
                         fg_84, fg_130, fg_131, fg_132, fg_133, gf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_5 * ff_88[k]
                   + pb_x[k] * gf_88[k];

        t_129[k] = pa_y[k] * fg_84[k];

        t_130[k] = pa_x[k] * fg_130[k];

        t_131[k] = pa_x[k] * fg_131[k];

        t_132[k] = pa_x[k] * fg_132[k];

        t_133[k] = pa_x[k] * fg_133[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, pa_x, pb_y, pb_z, ff_50, ff_90, \
                         ff_93, fg_134, fg_135, fg_138, gf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = pa_x[k] * fg_134[k];

        t_135[k] = f_0 * ff_90[k]
                   + pa_x[k] * fg_135[k];

        t_136[k] = pb_y[k] * gf_90[k];

        t_137[k] = f_7 * ff_50[k]
                   + pb_z[k] * gf_90[k];

        t_138[k] = f_6 * ff_93[k]
                   + pa_x[k] * fg_138[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pa_x, pb_x, pb_y, ff_95, ff_96, \
                         ff_97, fg_140, gf_92, gf_95, gf_96, gf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * gf_92[k];

        t_140[k] = f_6 * ff_95[k]
                   + pa_x[k] * fg_140[k];

        t_141[k] = f_5 * ff_96[k]
                   + pb_x[k] * gf_96[k];

        t_142[k] = f_5 * ff_97[k]
                   + pb_x[k] * gf_97[k];

        t_143[k] = pb_y[k] * gf_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, t_149, pa_x, pb_x, pb_y, ff_99, \
                         fg_145, fg_146, fg_147, fg_149, gf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_5 * ff_99[k]
                   + pb_x[k] * gf_99[k];

        t_145[k] = pa_x[k] * fg_145[k];

        t_146[k] = pa_x[k] * fg_146[k];

        t_147[k] = pa_x[k] * fg_147[k];

        t_148[k] = pb_y[k] * gf_99[k];

        t_149[k] = pa_x[k] * fg_149[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, pb_x, pb_y, pb_z, ff_60, gd0_60, \
                         gd0_63, gd1_60, gd1_63, gf_100, gf_101, \
                         gf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_1 * gd0_60[k]
                   - f_2 * gd1_60[k]
                   + pb_x[k] * gf_100[k];

        t_151[k] = f_0 * ff_60[k]
                   + pb_y[k] * gf_100[k];

        t_152[k] = pb_z[k] * gf_100[k];

        t_153[k] = f_3 * gd0_63[k]
                   - f_4 * gd1_63[k]
                   + pb_x[k] * gf_103[k];

        t_154[k] = pb_z[k] * gf_101[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pb_x, gd0_65, gd1_65, gf_105, \
                         gf_106, gf_107, gf_108, gf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_3 * gd0_65[k]
                   - f_4 * gd1_65[k]
                   + pb_x[k] * gf_105[k];

        t_156[k] = pb_x[k] * gf_106[k];

        t_157[k] = pb_x[k] * gf_107[k];

        t_158[k] = pb_x[k] * gf_108[k];

        t_159[k] = pb_x[k] * gf_109[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, pb_y, pb_z, ff_66, ff_69, gd0_63, \
                         gd0_65, gd1_63, gd1_65, gf_106, gf_107, \
                         gf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_0 * ff_66[k]
                   + f_1 * gd0_63[k]
                   - f_2 * gd1_63[k]
                   + pb_y[k] * gf_106[k];

        t_161[k] = pb_z[k] * gf_106[k];

        t_162[k] = f_3 * gd0_63[k]
                   - f_4 * gd1_63[k]
                   + pb_z[k] * gf_107[k];

        t_163[k] = f_0 * ff_69[k]
                   + pb_y[k] * gf_109[k];

        t_164[k] = f_1 * gd0_65[k]
                   - f_2 * gd1_65[k]
                   + pb_z[k] * gf_109[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pa_z, pb_y, pb_z, ff_60, ff_72, \
                         fg_90, fg_91, fg_93, gf_110, gf_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = pa_z[k] * fg_90[k];

        t_166[k] = pa_z[k] * fg_91[k];

        t_167[k] = f_5 * ff_60[k]
                   + pb_z[k] * gf_110[k];

        t_168[k] = pa_z[k] * fg_93[k];

        t_169[k] = f_7 * ff_72[k]
                   + pb_y[k] * gf_112[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, t_175, pa_z, pb_x, ff_62, fg_95, \
                         fg_100, gf_116, gf_117, gf_118, gf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_6 * ff_62[k]
                   + pa_z[k] * fg_95[k];

        t_171[k] = pb_x[k] * gf_116[k];

        t_172[k] = pb_x[k] * gf_117[k];

        t_173[k] = pb_x[k] * gf_118[k];

        t_174[k] = pb_x[k] * gf_119[k];

        t_175[k] = pa_z[k] * fg_100[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_z, pb_y, pb_z, ff_66, ff_67, ff_69, \
                         ff_79, fg_102, fg_104, gf_116, gf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_5 * ff_66[k]
                   + pb_z[k] * gf_116[k];

        t_177[k] = f_6 * ff_67[k]
                   + pa_z[k] * fg_102[k];

        t_178[k] = f_7 * ff_79[k]
                   + pb_y[k] * gf_119[k];

        t_179[k] = f_0 * ff_69[k]
                   + pa_z[k] * fg_104[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pb_x, pb_y, pb_z, ff_70, ff_80, gd0_72, \
                         gd0_75, gd1_72, gd1_75, gf_120, gf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_1 * gd0_72[k]
                   - f_2 * gd1_72[k]
                   + pb_x[k] * gf_120[k];

        t_181[k] = f_6 * ff_80[k]
                   + pb_y[k] * gf_120[k];

        t_182[k] = f_6 * ff_70[k]
                   + pb_z[k] * gf_120[k];

        t_183[k] = f_3 * gd0_75[k]
                   - f_4 * gd1_75[k]
                   + pb_x[k] * gf_123[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pb_x, pb_y, ff_82, gd0_77, gd1_77, \
                         gf_122, gf_125, gf_126, gf_127, gf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_6 * ff_82[k]
                   + pb_y[k] * gf_122[k];

        t_185[k] = f_3 * gd0_77[k]
                   - f_4 * gd1_77[k]
                   + pb_x[k] * gf_125[k];

        t_186[k] = pb_x[k] * gf_126[k];

        t_187[k] = pb_x[k] * gf_127[k];

        t_188[k] = pb_x[k] * gf_128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pa_z, pb_x, pb_z, dg0_55, dg1_55, ff_76, fg_115, \
                         gf_126, gf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pb_x[k] * gf_129[k];

        t_190[k] = f_8 * dg0_55[k]
                   - f_9 * dg1_55[k]
                   + pa_z[k] * fg_115[k];

        t_191[k] = f_6 * ff_76[k]
                   + pb_z[k] * gf_126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_y, pb_y, dg0_89, dg1_89, ff_88, ff_89, \
                         fg_134, fg_135, gd0_77, gd1_77, gf_128, \
                         gf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_6 * ff_88[k]
                   + f_3 * gd0_77[k]
                   - f_4 * gd1_77[k]
                   + pb_y[k] * gf_128[k];

        t_193[k] = f_6 * ff_89[k]
                   + pb_y[k] * gf_129[k];

        t_194[k] = f_8 * dg0_89[k]
                   - f_9 * dg1_89[k]
                   + pa_y[k] * fg_134[k];

        t_195[k] = pa_y[k] * fg_135[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, pa_y, pb_y, ff_90, ff_91, ff_92, \
                         fg_137, fg_138, fg_140, gf_130, gf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_5 * ff_90[k]
                   + pb_y[k] * gf_130[k];

        t_197[k] = pa_y[k] * fg_137[k];

        t_198[k] = f_6 * ff_91[k]
                   + pa_y[k] * fg_138[k];

        t_199[k] = f_5 * ff_92[k]
                   + pb_y[k] * gf_132[k];

        t_200[k] = pa_y[k] * fg_140[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, t_206, pa_y, pb_x, pb_z, ff_86, \
                         ff_96, fg_145, gf_136, gf_137, gf_138, \
                         gf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pb_x[k] * gf_136[k];

        t_202[k] = pb_x[k] * gf_137[k];

        t_203[k] = pb_x[k] * gf_138[k];

        t_204[k] = pb_x[k] * gf_139[k];

        t_205[k] = f_0 * ff_96[k]
                   + pa_y[k] * fg_145[k];

        t_206[k] = f_7 * ff_86[k]
                   + pb_z[k] * gf_136[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, pa_y, pb_x, pb_y, ff_98, ff_99, \
                         fg_147, fg_149, gd0_84, gd1_84, gf_139, \
                         gf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_6 * ff_98[k]
                   + pa_y[k] * fg_147[k];

        t_208[k] = f_5 * ff_99[k]
                   + pb_y[k] * gf_139[k];

        t_209[k] = pa_y[k] * fg_149[k];

        t_210[k] = f_1 * gd0_84[k]
                   - f_2 * gd1_84[k]
                   + pb_x[k] * gf_140[k];

        t_211[k] = pb_y[k] * gf_140[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pb_x, pb_y, pb_z, ff_90, gd0_87, gd0_89, \
                         gd1_87, gd1_89, gf_140, gf_142, gf_143, \
                         gf_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_0 * ff_90[k]
                   + pb_z[k] * gf_140[k];

        t_213[k] = f_3 * gd0_87[k]
                   - f_4 * gd1_87[k]
                   + pb_x[k] * gf_143[k];

        t_214[k] = pb_y[k] * gf_142[k];

        t_215[k] = f_3 * gd0_89[k]
                   - f_4 * gd1_89[k]
                   + pb_x[k] * gf_145[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, t_221, pb_x, pb_y, pb_z, ff_96, \
                         gd0_87, gd1_87, gf_146, gf_147, gf_148, \
                         gf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = pb_x[k] * gf_146[k];

        t_217[k] = pb_x[k] * gf_147[k];

        t_218[k] = pb_x[k] * gf_148[k];

        t_219[k] = pb_x[k] * gf_149[k];

        t_220[k] = f_1 * gd0_87[k]
                   - f_2 * gd1_87[k]
                   + pb_y[k] * gf_146[k];

        t_221[k] = f_0 * ff_96[k]
                   + pb_z[k] * gf_146[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pb_y, pb_z, ff_99, gd0_89, gd1_89, gf_148, \
                         gf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_3 * gd0_89[k]
                   - f_4 * gd1_89[k]
                   + pb_y[k] * gf_148[k];

        t_223[k] = pb_y[k] * gf_149[k];

        t_224[k] = f_0 * ff_99[k]
                   + f_1 * gd0_89[k]
                   - f_2 * gd1_89[k]
                   + pb_z[k] * gf_149[k];
    }
}

}  // namespace simdt2ceri
