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


#include "SimdElectronRepulsionVrrRecFH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_fh_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg, const size_t dh,
                                     const size_t ff0, const size_t ff1, const size_t fg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_40 = buffer.data(dg + 40);
    const auto *dg_41 = buffer.data(dg + 41);
    const auto *dg_42 = buffer.data(dg + 42);
    const auto *dg_43 = buffer.data(dg + 43);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_46 = buffer.data(dg + 46);
    const auto *dg_47 = buffer.data(dg + 47);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_49 = buffer.data(dg + 49);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_51 = buffer.data(dg + 51);
    const auto *dg_52 = buffer.data(dg + 52);
    const auto *dg_53 = buffer.data(dg + 53);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_37 = buffer.data(dh + 37);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_41 = buffer.data(dh + 41);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_43 = buffer.data(dh + 43);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_45 = buffer.data(dh + 45);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_47 = buffer.data(dh + 47);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_9 = buffer.data(ff0 + 9);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_9 = buffer.data(ff1 + 9);
    const auto *ff1_10 = buffer.data(ff1 + 10);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);

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
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);
    const auto *fg_59 = buffer.data(fg + 59);
    const auto *fg_60 = buffer.data(fg + 60);
    const auto *fg_61 = buffer.data(fg + 61);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_63 = buffer.data(fg + 63);
    const auto *fg_64 = buffer.data(fg + 64);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_66 = buffer.data(fg + 66);
    const auto *fg_67 = buffer.data(fg + 67);
    const auto *fg_68 = buffer.data(fg + 68);
    const auto *fg_69 = buffer.data(fg + 69);
    const auto *fg_70 = buffer.data(fg + 70);
    const auto *fg_71 = buffer.data(fg + 71);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_73 = buffer.data(fg + 73);
    const auto *fg_74 = buffer.data(fg + 74);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_76 = buffer.data(fg + 76);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_78 = buffer.data(fg + 78);
    const auto *fg_79 = buffer.data(fg + 79);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_81 = buffer.data(fg + 81);
    const auto *fg_82 = buffer.data(fg + 82);
    const auto *fg_83 = buffer.data(fg + 83);
    const auto *fg_84 = buffer.data(fg + 84);
    const auto *fg_85 = buffer.data(fg + 85);
    const auto *fg_86 = buffer.data(fg + 86);
    const auto *fg_87 = buffer.data(fg + 87);
    const auto *fg_88 = buffer.data(fg + 88);
    const auto *fg_89 = buffer.data(fg + 89);
    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_91 = buffer.data(fg + 91);
    const auto *fg_92 = buffer.data(fg + 92);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, \
                         fg_0, fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = pb_y[k] * fg_2[k];

        t_5[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, dg_5, ff0_1, ff0_2, \
                         ff1_1, ff1_2, fg_3, fg_4, fg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_7[k] = pb_z[k] * fg_3[k];

        t_8[k] = pb_y[k] * fg_4[k];

        t_9[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_10[k] = f_0 * dg_5[k]
                  + pb_x[k] * fg_7[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pb_x, pb_y, pb_z, dg_7, dg_9, fg_5, fg_6, \
                         fg_8, fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * fg_5[k];

        t_12[k] = f_0 * dg_7[k]
                  + pb_x[k] * fg_8[k];

        t_13[k] = pb_y[k] * fg_6[k];

        t_14[k] = f_0 * dg_9[k]
                  + pb_x[k] * fg_10[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pb_z, ff0_3, ff0_4, ff0_5, ff1_3, \
                         ff1_4, ff1_5, fg_7, fg_8, fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * ff0_3[k]
                  - f_2 * ff1_3[k]
                  + pb_y[k] * fg_7[k];

        t_16[k] = pb_z[k] * fg_7[k];

        t_17[k] = f_5 * ff0_4[k]
                  - f_6 * ff1_4[k]
                  + pb_y[k] * fg_8[k];

        t_18[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pb_y, pb_z, dg_0, dh_0, ff0_5, \
                         ff1_5, fg_10, fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * fg_10[k];

        t_20[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_10[k];

        t_21[k] = pa_y[k] * dh_0[k];

        t_22[k] = f_7 * dg_0[k]
                  + pb_y[k] * fg_11[k];

        t_23[k] = pb_z[k] * fg_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_y, pb_z, dg_1, dg_3, dh_1, dh_2, \
                         dh_3, fg_12, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * dg_1[k]
                  + pa_y[k] * dh_1[k];

        t_25[k] = pb_z[k] * fg_12[k];

        t_26[k] = pa_y[k] * dh_2[k];

        t_27[k] = f_0 * dg_3[k]
                  + pa_y[k] * dh_3[k];

        t_28[k] = pb_z[k] * fg_13[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, dg_4, dg_13, dh_4, \
                         fg_14, fg_15, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_7 * dg_4[k]
                  + pb_y[k] * fg_14[k];

        t_30[k] = pa_y[k] * dh_4[k];

        t_31[k] = f_8 * dg_13[k]
                  + pb_x[k] * fg_16[k];

        t_32[k] = pb_z[k] * fg_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pb_x, pb_z, dg_5, dg_14, dg_15, \
                         dh_6, dh_7, fg_16, fg_17, fg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_8 * dg_14[k]
                  + pb_x[k] * fg_17[k];

        t_34[k] = f_8 * dg_15[k]
                  + pb_x[k] * fg_18[k];

        t_35[k] = pa_y[k] * dh_6[k];

        t_36[k] = f_9 * dg_5[k]
                  + pa_y[k] * dh_7[k];

        t_37[k] = pb_z[k] * fg_16[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, dg_7, dg_8, dg_9, \
                         dh_0, dh_8, dh_9, dh_10, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * dg_7[k]
                  + pa_y[k] * dh_8[k];

        t_39[k] = f_8 * dg_8[k]
                  + pa_y[k] * dh_9[k];

        t_40[k] = f_7 * dg_9[k]
                  + pb_y[k] * fg_19[k];

        t_41[k] = pa_y[k] * dh_10[k];

        t_42[k] = pa_z[k] * dh_0[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, t_48, pa_z, pb_y, pb_z, dg_0, dg_2, \
                         dh_1, dh_2, dh_3, fg_20, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * fg_20[k];

        t_44[k] = f_7 * dg_0[k]
                  + pb_z[k] * fg_20[k];

        t_45[k] = pa_z[k] * dh_1[k];

        t_46[k] = pb_y[k] * fg_21[k];

        t_47[k] = f_8 * dg_2[k]
                  + pa_z[k] * dh_2[k];

        t_48[k] = pa_z[k] * dh_3[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_z, pb_y, pb_z, dg_3, dg_4, dh_4, dh_5, \
                         fg_22, fg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_7 * dg_3[k]
                  + pb_z[k] * fg_22[k];

        t_50[k] = pb_y[k] * fg_23[k];

        t_51[k] = f_0 * dg_4[k]
                  + pa_z[k] * dh_4[k];

        t_52[k] = pa_z[k] * dh_5[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_z, pb_x, pb_y, dg_20, dg_21, dg_22, \
                         dh_7, fg_24, fg_26, fg_27, fg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_8 * dg_20[k]
                  + pb_x[k] * fg_26[k];

        t_54[k] = f_8 * dg_21[k]
                  + pb_x[k] * fg_27[k];

        t_55[k] = pb_y[k] * fg_24[k];

        t_56[k] = f_8 * dg_22[k]
                  + pb_x[k] * fg_28[k];

        t_57[k] = pa_z[k] * dh_7[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pb_y, pb_z, dg_5, dg_6, dg_7, dh_8, \
                         dh_9, fg_25, fg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_7 * dg_5[k]
                  + pb_z[k] * fg_25[k];

        t_59[k] = f_8 * dg_6[k]
                  + pa_z[k] * dh_8[k];

        t_60[k] = f_0 * dg_7[k]
                  + pa_z[k] * dh_9[k];

        t_61[k] = pb_y[k] * fg_28[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pa_z, pb_y, pb_z, dg_9, dg_10, dg_23, \
                         dh_10, dh_20, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_9 * dg_9[k]
                  + pa_z[k] * dh_10[k];

        t_63[k] = f_9 * dg_23[k]
                  + pa_x[k] * dh_20[k];

        t_64[k] = f_8 * dg_10[k]
                  + pb_y[k] * fg_29[k];

        t_65[k] = pb_z[k] * fg_29[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_x, pb_z, dg_25, dg_26, dg_27, dh_22, \
                         dh_23, dh_24, fg_30, fg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_0 * dg_25[k]
                  + pa_x[k] * dh_22[k];

        t_67[k] = pb_z[k] * fg_30[k];

        t_68[k] = f_0 * dg_26[k]
                  + pa_x[k] * dh_23[k];

        t_69[k] = f_8 * dg_27[k]
                  + pa_x[k] * dh_24[k];

        t_70[k] = pb_z[k] * fg_31[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_x, pb_x, pb_y, pb_z, dg_12, dg_28, dg_29, \
                         dh_25, fg_32, fg_33, fg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_8 * dg_12[k]
                  + pb_y[k] * fg_32[k];

        t_72[k] = f_8 * dg_28[k]
                  + pa_x[k] * dh_25[k];

        t_73[k] = f_7 * dg_29[k]
                  + pb_x[k] * fg_34[k];

        t_74[k] = pb_z[k] * fg_33[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pa_x, pb_x, pb_z, dg_31, dg_32, dg_33, \
                         dh_26, fg_34, fg_35, fg_36, fg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_7 * dg_31[k]
                  + pb_x[k] * fg_35[k];

        t_76[k] = f_7 * dg_32[k]
                  + pb_x[k] * fg_36[k];

        t_77[k] = f_7 * dg_33[k]
                  + pb_x[k] * fg_37[k];

        t_78[k] = pa_x[k] * dh_26[k];

        t_79[k] = pb_z[k] * fg_34[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, t_85, pa_x, pa_y, pa_z, dh_11, dh_15, \
                         dh_27, dh_28, dh_29, dh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pa_x[k] * dh_27[k];

        t_81[k] = pa_x[k] * dh_28[k];

        t_82[k] = pa_x[k] * dh_29[k];

        t_83[k] = pa_x[k] * dh_30[k];

        t_84[k] = pa_y[k] * dh_15[k];

        t_85[k] = pa_z[k] * dh_11[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, pa_y, pa_z, pb_y, dg_17, dh_12, dh_13, \
                         dh_16, dh_17, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pa_y[k] * dh_16[k];

        t_87[k] = pa_z[k] * dh_12[k];

        t_88[k] = f_7 * dg_17[k]
                  + pb_y[k] * fg_38[k];

        t_89[k] = pa_y[k] * dh_17[k];

        t_90[k] = pa_z[k] * dh_13[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pa_z, pb_y, pb_z, dg_11, dg_19, dh_14, \
                         dh_18, fg_39, fg_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_7 * dg_11[k]
                  + pb_z[k] * fg_39[k];

        t_92[k] = f_7 * dg_19[k]
                  + pb_y[k] * fg_40[k];

        t_93[k] = pa_y[k] * dh_18[k];

        t_94[k] = pa_z[k] * dh_14[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pa_x, pa_y, pb_x, dg_38, dg_39, dg_40, \
                         dh_19, dh_31, fg_41, fg_42, fg_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_7 * dg_38[k]
                  + pb_x[k] * fg_41[k];

        t_96[k] = f_7 * dg_39[k]
                  + pb_x[k] * fg_42[k];

        t_97[k] = f_7 * dg_40[k]
                  + pb_x[k] * fg_43[k];

        t_98[k] = pa_y[k] * dh_19[k];

        t_99[k] = pa_x[k] * dh_31[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, t_105, pa_x, dg_42, dh_32, dh_33, \
                         dh_34, dh_35, dh_36, dh_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_x[k] * dh_32[k];

        t_101[k] = pa_x[k] * dh_33[k];

        t_102[k] = pa_x[k] * dh_34[k];

        t_103[k] = pa_x[k] * dh_35[k];

        t_104[k] = pa_x[k] * dh_36[k];

        t_105[k] = f_9 * dg_42[k]
                   + pa_x[k] * dh_37[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, pa_x, pb_y, pb_z, dg_16, dg_45, \
                         dg_46, dh_39, dh_40, fg_44, fg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = pb_y[k] * fg_44[k];

        t_107[k] = f_8 * dg_16[k]
                   + pb_z[k] * fg_44[k];

        t_108[k] = f_0 * dg_45[k]
                   + pa_x[k] * dh_39[k];

        t_109[k] = pb_y[k] * fg_45[k];

        t_110[k] = f_0 * dg_46[k]
                   + pa_x[k] * dh_40[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_x, pb_y, pb_z, dg_18, dg_47, dg_48, \
                         dh_41, dh_42, fg_46, fg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_8 * dg_47[k]
                   + pa_x[k] * dh_41[k];

        t_112[k] = f_8 * dg_18[k]
                   + pb_z[k] * fg_46[k];

        t_113[k] = pb_y[k] * fg_47[k];

        t_114[k] = f_8 * dg_48[k]
                   + pa_x[k] * dh_42[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pb_x, pb_y, dg_49, dg_50, dg_51, \
                         dg_53, fg_48, fg_49, fg_50, fg_51, fg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_7 * dg_49[k]
                   + pb_x[k] * fg_49[k];

        t_116[k] = f_7 * dg_50[k]
                   + pb_x[k] * fg_50[k];

        t_117[k] = f_7 * dg_51[k]
                   + pb_x[k] * fg_51[k];

        t_118[k] = pb_y[k] * fg_48[k];

        t_119[k] = f_7 * dg_53[k]
                   + pb_x[k] * fg_52[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, t_125, pa_x, pb_y, dh_43, dh_44, \
                         dh_45, dh_46, dh_47, fg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_x[k] * dh_43[k];

        t_121[k] = pa_x[k] * dh_44[k];

        t_122[k] = pa_x[k] * dh_45[k];

        t_123[k] = pa_x[k] * dh_46[k];

        t_124[k] = pb_y[k] * fg_52[k];

        t_125[k] = pa_x[k] * dh_47[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, pb_x, pb_y, pb_z, dg_23, ff0_6, \
                         ff0_7, ff1_6, ff1_7, fg_53, fg_54, fg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_1 * ff0_6[k]
                   - f_2 * ff1_6[k]
                   + pb_x[k] * fg_53[k];

        t_127[k] = f_0 * dg_23[k]
                   + pb_y[k] * fg_53[k];

        t_128[k] = pb_z[k] * fg_53[k];

        t_129[k] = f_5 * ff0_7[k]
                   - f_6 * ff1_7[k]
                   + pb_x[k] * fg_55[k];

        t_130[k] = pb_z[k] * fg_54[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pb_x, pb_y, pb_z, dg_26, ff0_8, ff0_9, \
                         ff1_8, ff1_9, fg_55, fg_56, fg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_5 * ff0_8[k]
                   - f_6 * ff1_8[k]
                   + pb_x[k] * fg_56[k];

        t_132[k] = f_3 * ff0_9[k]
                   - f_4 * ff1_9[k]
                   + pb_x[k] * fg_57[k];

        t_133[k] = pb_z[k] * fg_55[k];

        t_134[k] = f_0 * dg_26[k]
                   + pb_y[k] * fg_56[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, t_140, pb_x, ff0_11, ff1_11, \
                         fg_58, fg_59, fg_60, fg_61, fg_62, fg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_3 * ff0_11[k]
                   - f_4 * ff1_11[k]
                   + pb_x[k] * fg_58[k];

        t_136[k] = pb_x[k] * fg_59[k];

        t_137[k] = pb_x[k] * fg_60[k];

        t_138[k] = pb_x[k] * fg_61[k];

        t_139[k] = pb_x[k] * fg_62[k];

        t_140[k] = pb_x[k] * fg_63[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pb_y, pb_z, dg_29, ff0_9, ff0_10, ff1_9, \
                         ff1_10, fg_59, fg_60, fg_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_0 * dg_29[k]
                   + f_1 * ff0_9[k]
                   - f_2 * ff1_9[k]
                   + pb_y[k] * fg_59[k];

        t_142[k] = pb_z[k] * fg_59[k];

        t_143[k] = f_3 * ff0_9[k]
                   - f_4 * ff1_9[k]
                   + pb_z[k] * fg_60[k];

        t_144[k] = f_5 * ff0_10[k]
                   - f_6 * ff1_10[k]
                   + pb_z[k] * fg_61[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pa_z, pb_y, pb_z, dg_23, dg_33, \
                         dh_20, dh_21, ff0_11, ff1_11, fg_63, fg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_0 * dg_33[k]
                   + pb_y[k] * fg_63[k];

        t_146[k] = f_1 * ff0_11[k]
                   - f_2 * ff1_11[k]
                   + pb_z[k] * fg_63[k];

        t_147[k] = pa_z[k] * dh_20[k];

        t_148[k] = pa_z[k] * dh_21[k];

        t_149[k] = f_7 * dg_23[k]
                   + pb_z[k] * fg_64[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, pa_z, pb_y, pb_z, dg_24, dg_25, \
                         dg_34, dh_22, dh_23, dh_24, fg_65, fg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * dh_22[k];

        t_151[k] = f_8 * dg_34[k]
                   + pb_y[k] * fg_65[k];

        t_152[k] = f_8 * dg_24[k]
                   + pa_z[k] * dh_23[k];

        t_153[k] = pa_z[k] * dh_24[k];

        t_154[k] = f_7 * dg_25[k]
                   + pb_z[k] * fg_66[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pa_z, pb_x, pb_y, dg_26, dg_36, \
                         dh_25, fg_67, fg_68, fg_69, fg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_8 * dg_36[k]
                   + pb_y[k] * fg_67[k];

        t_156[k] = f_0 * dg_26[k]
                   + pa_z[k] * dh_25[k];

        t_157[k] = pb_x[k] * fg_68[k];

        t_158[k] = pb_x[k] * fg_69[k];

        t_159[k] = pb_x[k] * fg_70[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, pa_z, pb_x, pb_z, dg_29, dg_30, \
                         dh_26, dh_27, fg_68, fg_71, fg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = pb_x[k] * fg_71[k];

        t_161[k] = pb_x[k] * fg_72[k];

        t_162[k] = pa_z[k] * dh_26[k];

        t_163[k] = f_7 * dg_29[k]
                   + pb_z[k] * fg_68[k];

        t_164[k] = f_8 * dg_30[k]
                   + pa_z[k] * dh_27[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, pa_y, pa_z, pb_y, dg_31, dg_33, dg_41, \
                         dh_28, dh_30, dh_37, fg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_0 * dg_31[k]
                   + pa_z[k] * dh_28[k];

        t_166[k] = f_8 * dg_41[k]
                   + pb_y[k] * fg_72[k];

        t_167[k] = f_9 * dg_33[k]
                   + pa_z[k] * dh_30[k];

        t_168[k] = pa_y[k] * dh_37[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, pa_y, pb_y, dg_42, dg_43, dg_44, \
                         dh_38, dh_39, dh_40, fg_73, fg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_7 * dg_42[k]
                   + pb_y[k] * fg_73[k];

        t_170[k] = pa_y[k] * dh_38[k];

        t_171[k] = f_8 * dg_43[k]
                   + pa_y[k] * dh_39[k];

        t_172[k] = f_7 * dg_44[k]
                   + pb_y[k] * fg_74[k];

        t_173[k] = pa_y[k] * dh_40[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pa_y, pb_y, pb_z, dg_35, dg_45, dg_46, \
                         dh_41, dh_42, fg_75, fg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_0 * dg_45[k]
                   + pa_y[k] * dh_41[k];

        t_175[k] = f_8 * dg_35[k]
                   + pb_z[k] * fg_75[k];

        t_176[k] = f_7 * dg_46[k]
                   + pb_y[k] * fg_76[k];

        t_177[k] = pa_y[k] * dh_42[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, t_183, pa_y, pb_x, dg_49, dh_43, \
                         fg_77, fg_78, fg_79, fg_80, fg_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = pb_x[k] * fg_77[k];

        t_179[k] = pb_x[k] * fg_78[k];

        t_180[k] = pb_x[k] * fg_79[k];

        t_181[k] = pb_x[k] * fg_80[k];

        t_182[k] = pb_x[k] * fg_81[k];

        t_183[k] = f_9 * dg_49[k]
                   + pa_y[k] * dh_43[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pb_y, pb_z, dg_37, dg_51, dg_52, \
                         dg_53, dh_45, dh_46, fg_77, fg_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_8 * dg_37[k]
                   + pb_z[k] * fg_77[k];

        t_185[k] = f_0 * dg_51[k]
                   + pa_y[k] * dh_45[k];

        t_186[k] = f_8 * dg_52[k]
                   + pa_y[k] * dh_46[k];

        t_187[k] = f_7 * dg_53[k]
                   + pb_y[k] * fg_81[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_y, pb_x, pb_y, pb_z, dg_42, dh_47, \
                         ff0_12, ff1_12, fg_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pa_y[k] * dh_47[k];

        t_189[k] = f_1 * ff0_12[k]
                   - f_2 * ff1_12[k]
                   + pb_x[k] * fg_82[k];

        t_190[k] = pb_y[k] * fg_82[k];

        t_191[k] = f_0 * dg_42[k]
                   + pb_z[k] * fg_82[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pb_x, pb_y, ff0_13, ff0_14, ff0_15, \
                         ff1_13, ff1_14, ff1_15, fg_83, fg_84, fg_85, \
                         fg_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_5 * ff0_13[k]
                   - f_6 * ff1_13[k]
                   + pb_x[k] * fg_84[k];

        t_193[k] = pb_y[k] * fg_83[k];

        t_194[k] = f_5 * ff0_14[k]
                   - f_6 * ff1_14[k]
                   + pb_x[k] * fg_85[k];

        t_195[k] = f_3 * ff0_15[k]
                   - f_4 * ff1_15[k]
                   + pb_x[k] * fg_86[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, pb_x, pb_y, pb_z, dg_45, ff0_17, \
                         ff1_17, fg_84, fg_85, fg_87, fg_88, fg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_0 * dg_45[k]
                   + pb_z[k] * fg_84[k];

        t_197[k] = pb_y[k] * fg_85[k];

        t_198[k] = f_3 * ff0_17[k]
                   - f_4 * ff1_17[k]
                   + pb_x[k] * fg_87[k];

        t_199[k] = pb_x[k] * fg_88[k];

        t_200[k] = pb_x[k] * fg_89[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, pb_x, pb_y, pb_z, dg_49, ff0_15, \
                         ff1_15, fg_88, fg_90, fg_91, fg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pb_x[k] * fg_90[k];

        t_202[k] = pb_x[k] * fg_91[k];

        t_203[k] = pb_x[k] * fg_92[k];

        t_204[k] = f_1 * ff0_15[k]
                   - f_2 * ff1_15[k]
                   + pb_y[k] * fg_88[k];

        t_205[k] = f_0 * dg_49[k]
                   + pb_z[k] * fg_88[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pb_y, pb_z, dg_53, ff0_16, ff0_17, \
                         ff1_16, ff1_17, fg_90, fg_91, fg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_5 * ff0_16[k]
                   - f_6 * ff1_16[k]
                   + pb_y[k] * fg_90[k];

        t_207[k] = f_3 * ff0_17[k]
                   - f_4 * ff1_17[k]
                   + pb_y[k] * fg_91[k];

        t_208[k] = pb_y[k] * fg_92[k];

        t_209[k] = f_0 * dg_53[k]
                   + f_1 * ff0_17[k]
                   - f_2 * ff1_17[k]
                   + pb_z[k] * fg_92[k];
    }
}

auto
compute_prim_fh_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg, const size_t dh,
                                     const size_t ff0, const size_t ff1, const size_t fg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_37 = buffer.data(dh + 37);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_43 = buffer.data(dh + 43);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_45 = buffer.data(dh + 45);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_54 = buffer.data(dh + 54);
    const auto *dh_56 = buffer.data(dh + 56);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_9 = buffer.data(ff0 + 9);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_9 = buffer.data(ff1 + 9);
    const auto *ff1_10 = buffer.data(ff1 + 10);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, dg_5, ff0_1, ff0_2, ff1_1, \
                         ff1_2, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_0 * dg_5[k]
                 + pb_x[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, pb_y, dg_9, ff0_3, ff0_4, ff1_3, ff1_4, fg_5, \
                         fg_6, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * dg_9[k]
                 + pb_x[k] * fg_8[k];

        t_10[k] = f_1 * ff0_3[k]
                  - f_2 * ff1_3[k]
                  + pb_y[k] * fg_5[k];

        t_11[k] = f_5 * ff0_4[k]
                  - f_6 * ff1_4[k]
                  + pb_y[k] * fg_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pa_y, pb_y, pb_z, dg_0, dh_0, ff0_5, \
                         ff1_5, fg_7, fg_8, fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_13[k] = pb_y[k] * fg_8[k];

        t_14[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_15[k] = pa_y[k] * dh_0[k];

        t_16[k] = f_7 * dg_0[k]
                  + pb_y[k] * fg_9[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_y, pb_x, dg_1, dg_3, dg_11, dh_3, \
                         dh_4, dh_5, dh_7, fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * dg_1[k]
                  + pa_y[k] * dh_3[k];

        t_18[k] = pa_y[k] * dh_4[k];

        t_19[k] = f_0 * dg_3[k]
                  + pa_y[k] * dh_5[k];

        t_20[k] = pa_y[k] * dh_7[k];

        t_21[k] = f_8 * dg_11[k]
                  + pb_x[k] * fg_10[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pa_y, pb_y, dg_5, dg_7, dg_8, dg_9, \
                         dh_8, dh_9, dh_10, dh_12, fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_9 * dg_5[k]
                  + pa_y[k] * dh_8[k];

        t_23[k] = f_0 * dg_7[k]
                  + pa_y[k] * dh_9[k];

        t_24[k] = f_8 * dg_8[k]
                  + pa_y[k] * dh_10[k];

        t_25[k] = f_7 * dg_9[k]
                  + pb_y[k] * fg_11[k];

        t_26[k] = pa_y[k] * dh_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_z, pb_z, dg_0, dg_2, dh_0, dh_3, \
                         dh_4, dh_5, fg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_z[k] * dh_0[k];

        t_28[k] = f_7 * dg_0[k]
                  + pb_z[k] * fg_12[k];

        t_29[k] = pa_z[k] * dh_3[k];

        t_30[k] = f_8 * dg_2[k]
                  + pa_z[k] * dh_4[k];

        t_31[k] = pa_z[k] * dh_5[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_z, pb_x, pb_z, dg_4, dg_5, dg_13, dh_7, \
                         dh_8, fg_13, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * dg_4[k]
                  + pa_z[k] * dh_7[k];

        t_33[k] = f_8 * dg_13[k]
                  + pb_x[k] * fg_14[k];

        t_34[k] = pa_z[k] * dh_8[k];

        t_35[k] = f_7 * dg_5[k]
                  + pb_z[k] * fg_13[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pa_z, dg_6, dg_7, dg_9, dg_14, dh_9, \
                         dh_10, dh_12, dh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_8 * dg_6[k]
                  + pa_z[k] * dh_9[k];

        t_37[k] = f_0 * dg_7[k]
                  + pa_z[k] * dh_10[k];

        t_38[k] = f_9 * dg_9[k]
                  + pa_z[k] * dh_12[k];

        t_39[k] = f_9 * dg_14[k]
                  + pa_x[k] * dh_22[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_y, dg_10, dg_16, dg_17, dg_18, \
                         dh_23, dh_24, dh_25, fg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_8 * dg_10[k]
                  + pb_y[k] * fg_15[k];

        t_41[k] = f_0 * dg_16[k]
                  + pa_x[k] * dh_23[k];

        t_42[k] = f_0 * dg_17[k]
                  + pa_x[k] * dh_24[k];

        t_43[k] = f_8 * dg_18[k]
                  + pa_x[k] * dh_25[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, pa_x, pb_x, dg_19, dg_20, dh_26, \
                         dh_30, dh_32, dh_33, dh_34, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_8 * dg_19[k]
                  + pa_x[k] * dh_26[k];

        t_45[k] = f_7 * dg_20[k]
                  + pb_x[k] * fg_16[k];

        t_46[k] = pa_x[k] * dh_30[k];

        t_47[k] = pa_x[k] * dh_32[k];

        t_48[k] = pa_x[k] * dh_33[k];

        t_49[k] = pa_x[k] * dh_34[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, t_55, pa_x, pa_y, pa_z, dh_14, dh_15, \
                         dh_18, dh_19, dh_20, dh_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_x[k] * dh_35[k];

        t_51[k] = pa_y[k] * dh_18[k];

        t_52[k] = pa_z[k] * dh_14[k];

        t_53[k] = pa_y[k] * dh_19[k];

        t_54[k] = pa_z[k] * dh_15[k];

        t_55[k] = pa_y[k] * dh_20[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, pa_x, pb_z, dg_12, dg_27, dh_37, \
                         dh_38, dh_39, dh_40, dh_42, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_x[k] * dh_37[k];

        t_57[k] = pa_x[k] * dh_38[k];

        t_58[k] = pa_x[k] * dh_39[k];

        t_59[k] = pa_x[k] * dh_40[k];

        t_60[k] = f_9 * dg_27[k]
                  + pa_x[k] * dh_42[k];

        t_61[k] = f_8 * dg_12[k]
                  + pb_z[k] * fg_17[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, dg_29, dg_30, dg_31, dg_32, dh_44, \
                         dh_45, dh_46, dh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_0 * dg_29[k]
                  + pa_x[k] * dh_44[k];

        t_63[k] = f_0 * dg_30[k]
                  + pa_x[k] * dh_45[k];

        t_64[k] = f_8 * dg_31[k]
                  + pa_x[k] * dh_46[k];

        t_65[k] = f_8 * dg_32[k]
                  + pa_x[k] * dh_47[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, pa_x, pb_x, dg_36, dh_51, dh_52, \
                         dh_53, dh_54, dh_56, fg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_7 * dg_36[k]
                  + pb_x[k] * fg_18[k];

        t_67[k] = pa_x[k] * dh_51[k];

        t_68[k] = pa_x[k] * dh_52[k];

        t_69[k] = pa_x[k] * dh_53[k];

        t_70[k] = pa_x[k] * dh_54[k];

        t_71[k] = pa_x[k] * dh_56[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pb_x, pb_y, dg_14, ff0_6, ff0_7, ff0_8, \
                         ff1_6, ff1_7, ff1_8, fg_19, fg_20, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * ff0_6[k]
                  - f_2 * ff1_6[k]
                  + pb_x[k] * fg_19[k];

        t_73[k] = f_0 * dg_14[k]
                  + pb_y[k] * fg_19[k];

        t_74[k] = f_5 * ff0_7[k]
                  - f_6 * ff1_7[k]
                  + pb_x[k] * fg_20[k];

        t_75[k] = f_5 * ff0_8[k]
                  - f_6 * ff1_8[k]
                  + pb_x[k] * fg_21[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pb_x, ff0_9, ff0_11, ff1_9, ff1_11, \
                         fg_22, fg_23, fg_24, fg_26, fg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_3 * ff0_9[k]
                  - f_4 * ff1_9[k]
                  + pb_x[k] * fg_22[k];

        t_77[k] = f_3 * ff0_11[k]
                  - f_4 * ff1_11[k]
                  + pb_x[k] * fg_23[k];

        t_78[k] = pb_x[k] * fg_24[k];

        t_79[k] = pb_x[k] * fg_26[k];

        t_80[k] = pb_x[k] * fg_27[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pb_y, pb_z, dg_20, ff0_9, ff0_10, ff1_9, \
                         ff1_10, fg_24, fg_25, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_0 * dg_20[k]
                  + f_1 * ff0_9[k]
                  - f_2 * ff1_9[k]
                  + pb_y[k] * fg_24[k];

        t_82[k] = pb_z[k] * fg_24[k];

        t_83[k] = f_3 * ff0_9[k]
                  - f_4 * ff1_9[k]
                  + pb_z[k] * fg_25[k];

        t_84[k] = f_5 * ff0_10[k]
                  - f_6 * ff1_10[k]
                  + pb_z[k] * fg_26[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pa_z, pb_y, pb_z, dg_14, dg_23, dh_22, \
                         dh_23, ff0_11, ff1_11, fg_27, fg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_0 * dg_23[k]
                  + pb_y[k] * fg_27[k];

        t_86[k] = f_1 * ff0_11[k]
                  - f_2 * ff1_11[k]
                  + pb_z[k] * fg_27[k];

        t_87[k] = pa_z[k] * dh_22[k];

        t_88[k] = f_7 * dg_14[k]
                  + pb_z[k] * fg_28[k];

        t_89[k] = pa_z[k] * dh_23[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pa_z, pb_z, dg_15, dg_17, dg_20, dh_24, \
                         dh_25, dh_26, dh_30, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_8 * dg_15[k]
                  + pa_z[k] * dh_24[k];

        t_91[k] = pa_z[k] * dh_25[k];

        t_92[k] = f_0 * dg_17[k]
                  + pa_z[k] * dh_26[k];

        t_93[k] = pa_z[k] * dh_30[k];

        t_94[k] = f_7 * dg_20[k]
                  + pb_z[k] * fg_29[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_z, pb_y, dg_21, dg_22, dg_23, dg_26, \
                         dh_32, dh_33, dh_35, fg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_8 * dg_21[k]
                  + pa_z[k] * dh_32[k];

        t_96[k] = f_0 * dg_22[k]
                  + pa_z[k] * dh_33[k];

        t_97[k] = f_8 * dg_26[k]
                  + pb_y[k] * fg_30[k];

        t_98[k] = f_9 * dg_23[k]
                  + pa_z[k] * dh_35[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, pa_y, dg_28, dg_29, dh_42, \
                         dh_43, dh_44, dh_45, dh_46, dh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_y[k] * dh_42[k];

        t_100[k] = pa_y[k] * dh_43[k];

        t_101[k] = f_8 * dg_28[k]
                   + pa_y[k] * dh_44[k];

        t_102[k] = pa_y[k] * dh_45[k];

        t_103[k] = f_0 * dg_29[k]
                   + pa_y[k] * dh_46[k];

        t_104[k] = pa_y[k] * dh_47[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_y, pb_z, dg_24, dg_33, dg_34, dg_35, \
                         dh_51, dh_53, dh_54, fg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_9 * dg_33[k]
                   + pa_y[k] * dh_51[k];

        t_106[k] = f_8 * dg_24[k]
                   + pb_z[k] * fg_31[k];

        t_107[k] = f_0 * dg_34[k]
                   + pa_y[k] * dh_53[k];

        t_108[k] = f_8 * dg_35[k]
                   + pa_y[k] * dh_54[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_y, pb_x, pb_y, pb_z, dg_27, dg_36, \
                         dh_56, ff0_12, ff1_12, fg_32, fg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_7 * dg_36[k]
                   + pb_y[k] * fg_32[k];

        t_110[k] = pa_y[k] * dh_56[k];

        t_111[k] = f_1 * ff0_12[k]
                   - f_2 * ff1_12[k]
                   + pb_x[k] * fg_33[k];

        t_112[k] = f_0 * dg_27[k]
                   + pb_z[k] * fg_33[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pb_x, ff0_13, ff0_14, ff0_15, ff1_13, ff1_14, \
                         ff1_15, fg_34, fg_35, fg_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_5 * ff0_13[k]
                   - f_6 * ff1_13[k]
                   + pb_x[k] * fg_34[k];

        t_114[k] = f_5 * ff0_14[k]
                   - f_6 * ff1_14[k]
                   + pb_x[k] * fg_35[k];

        t_115[k] = f_3 * ff0_15[k]
                   - f_4 * ff1_15[k]
                   + pb_x[k] * fg_36[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, pb_x, pb_y, ff0_15, ff0_17, \
                         ff1_15, ff1_17, fg_37, fg_38, fg_39, fg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_3 * ff0_17[k]
                   - f_4 * ff1_17[k]
                   + pb_x[k] * fg_37[k];

        t_117[k] = pb_x[k] * fg_38[k];

        t_118[k] = pb_x[k] * fg_39[k];

        t_119[k] = pb_x[k] * fg_41[k];

        t_120[k] = f_1 * ff0_15[k]
                   - f_2 * ff1_15[k]
                   + pb_y[k] * fg_38[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_y, pb_z, dg_33, ff0_16, ff0_17, \
                         ff1_16, ff1_17, fg_38, fg_39, fg_40, fg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_0 * dg_33[k]
                   + pb_z[k] * fg_38[k];

        t_122[k] = f_5 * ff0_16[k]
                   - f_6 * ff1_16[k]
                   + pb_y[k] * fg_39[k];

        t_123[k] = f_3 * ff0_17[k]
                   - f_4 * ff1_17[k]
                   + pb_y[k] * fg_40[k];

        t_124[k] = pb_y[k] * fg_41[k];
    }

#pragma omp simd aligned(t_125, pb_z, dg_36, ff0_17, ff1_17, fg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_0 * dg_36[k]
                   + f_1 * ff0_17[k]
                   - f_2 * ff1_17[k]
                   + pb_z[k] * fg_41[k];
    }
}

auto
compute_prim_fh_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t dg, const size_t ff0, const size_t ff1,
                                     const size_t fg, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_9 = buffer.data(ff0 + 9);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_9 = buffer.data(ff1 + 9);
    const auto *ff1_10 = buffer.data(ff1 + 10);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, \
                         ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, ff0_4, ff0_5, ff1_4, ff1_5, fg_6, \
                         fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_10[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = pb_y[k] * fg_8[k];

        t_12[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, ff0_6, ff0_7, ff0_8, ff1_6, ff1_7, ff1_8, \
                         fg_9, fg_10, fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * ff0_6[k]
                  - f_2 * ff1_6[k]
                  + pb_x[k] * fg_9[k];

        t_14[k] = f_5 * ff0_7[k]
                  - f_6 * ff1_7[k]
                  + pb_x[k] * fg_10[k];

        t_15[k] = f_5 * ff0_8[k]
                  - f_6 * ff1_8[k]
                  + pb_x[k] * fg_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pb_x, ff0_9, ff0_11, ff1_9, ff1_11, \
                         fg_12, fg_13, fg_14, fg_16, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * ff0_9[k]
                  - f_4 * ff1_9[k]
                  + pb_x[k] * fg_12[k];

        t_17[k] = f_3 * ff0_11[k]
                  - f_4 * ff1_11[k]
                  + pb_x[k] * fg_13[k];

        t_18[k] = pb_x[k] * fg_14[k];

        t_19[k] = pb_x[k] * fg_16[k];

        t_20[k] = pb_x[k] * fg_17[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, pb_z, dg_14, ff0_9, ff0_10, ff1_9, \
                         ff1_10, fg_14, fg_15, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * dg_14[k]
                  + f_1 * ff0_9[k]
                  - f_2 * ff1_9[k]
                  + pb_y[k] * fg_14[k];

        t_22[k] = pb_z[k] * fg_14[k];

        t_23[k] = f_3 * ff0_9[k]
                  - f_4 * ff1_9[k]
                  + pb_z[k] * fg_15[k];

        t_24[k] = f_5 * ff0_10[k]
                  - f_6 * ff1_10[k]
                  + pb_z[k] * fg_16[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pb_x, pb_z, ff0_11, ff0_12, ff0_13, ff1_11, ff1_12, \
                         ff1_13, fg_17, fg_18, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * ff0_11[k]
                  - f_2 * ff1_11[k]
                  + pb_z[k] * fg_17[k];

        t_26[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_18[k];

        t_27[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_19[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_x, ff0_14, ff0_15, ff0_17, ff1_14, ff1_15, \
                         ff1_17, fg_20, fg_21, fg_22, fg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_20[k];

        t_29[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_21[k];

        t_30[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_22[k];

        t_31[k] = pb_x[k] * fg_23[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_x, pb_y, ff0_15, ff0_16, ff1_15, ff1_16, \
                         fg_23, fg_24, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_x[k] * fg_24[k];

        t_33[k] = pb_x[k] * fg_26[k];

        t_34[k] = f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_23[k];

        t_35[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_y[k] * fg_24[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_y, pb_z, dg_26, ff0_17, ff1_17, fg_25, \
                         fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_y[k] * fg_25[k];

        t_37[k] = pb_y[k] * fg_26[k];

        t_38[k] = f_0 * dg_26[k]
                  + f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_26[k];
    }
}

auto
compute_prim_fh_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg, const size_t dh,
                                     const size_t ff0, const size_t ff1, const size_t fg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 2.5 / p;

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_37 = buffer.data(dh + 37);
    const auto *dh_41 = buffer.data(dh + 41);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_9 = buffer.data(ff0 + 9);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_9 = buffer.data(ff1 + 9);
    const auto *ff1_10 = buffer.data(ff1 + 10);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, \
                         ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_10[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = pb_y[k] * fg_8[k];

        t_12[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_13[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, t_19, pa_x, pa_y, pa_z, dg_5, dg_8, \
                         dg_9, dh_0, dh_8, dh_12, dh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_7 * dg_5[k]
                  + pa_y[k] * dh_8[k];

        t_15[k] = pa_y[k] * dh_12[k];

        t_16[k] = pa_z[k] * dh_0[k];

        t_17[k] = pa_z[k] * dh_8[k];

        t_18[k] = f_7 * dg_8[k]
                  + pa_z[k] * dh_12[k];

        t_19[k] = f_7 * dg_9[k]
                  + pa_x[k] * dh_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, dg_18, dh_22, dh_29, dh_41, \
                         ff0_6, ff1_6, fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_x[k] * dh_22[k];

        t_21[k] = f_7 * dg_18[k]
                  + pa_x[k] * dh_29[k];

        t_22[k] = pa_x[k] * dh_41[k];

        t_23[k] = f_1 * ff0_6[k]
                  - f_2 * ff1_6[k]
                  + pb_x[k] * fg_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_x, ff0_7, ff0_8, ff0_9, ff1_7, ff1_8, ff1_9, \
                         fg_10, fg_11, fg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * ff0_7[k]
                  - f_6 * ff1_7[k]
                  + pb_x[k] * fg_10[k];

        t_25[k] = f_5 * ff0_8[k]
                  - f_6 * ff1_8[k]
                  + pb_x[k] * fg_11[k];

        t_26[k] = f_3 * ff0_9[k]
                  - f_4 * ff1_9[k]
                  + pb_x[k] * fg_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pb_x, pb_y, dg_14, ff0_9, ff0_11, \
                         ff1_9, ff1_11, fg_13, fg_14, fg_16, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_3 * ff0_11[k]
                  - f_4 * ff1_11[k]
                  + pb_x[k] * fg_13[k];

        t_28[k] = pb_x[k] * fg_14[k];

        t_29[k] = pb_x[k] * fg_16[k];

        t_30[k] = pb_x[k] * fg_17[k];

        t_31[k] = f_0 * dg_14[k]
                  + f_1 * ff0_9[k]
                  - f_2 * ff1_9[k]
                  + pb_y[k] * fg_14[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_z, ff0_9, ff0_10, ff0_11, ff1_9, ff1_10, \
                         ff1_11, fg_14, fg_15, fg_16, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_z[k] * fg_14[k];

        t_33[k] = f_3 * ff0_9[k]
                  - f_4 * ff1_9[k]
                  + pb_z[k] * fg_15[k];

        t_34[k] = f_5 * ff0_10[k]
                  - f_6 * ff1_10[k]
                  + pb_z[k] * fg_16[k];

        t_35[k] = f_1 * ff0_11[k]
                  - f_2 * ff1_11[k]
                  + pb_z[k] * fg_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pa_y, pa_z, dg_17, dg_23, dh_14, dh_22, \
                         dh_26, dh_37, dh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_z[k] * dh_14[k];

        t_37[k] = pa_z[k] * dh_22[k];

        t_38[k] = f_7 * dg_17[k]
                  + pa_z[k] * dh_26[k];

        t_39[k] = f_7 * dg_23[k]
                  + pa_y[k] * dh_37[k];

        t_40[k] = pa_y[k] * dh_41[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pb_x, ff0_12, ff0_13, ff0_14, ff1_12, ff1_13, \
                         ff1_14, fg_18, fg_19, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_18[k];

        t_42[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_19[k];

        t_43[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_20[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pb_x, ff0_15, ff0_17, ff1_15, ff1_17, \
                         fg_21, fg_22, fg_23, fg_24, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_21[k];

        t_45[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_22[k];

        t_46[k] = pb_x[k] * fg_23[k];

        t_47[k] = pb_x[k] * fg_24[k];

        t_48[k] = pb_x[k] * fg_26[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pb_y, ff0_15, ff0_16, ff0_17, ff1_15, ff1_16, \
                         ff1_17, fg_23, fg_24, fg_25, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_23[k];

        t_50[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_y[k] * fg_24[k];

        t_51[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_y[k] * fg_25[k];

        t_52[k] = pb_y[k] * fg_26[k];
    }

#pragma omp simd aligned(t_53, pb_z, dg_26, ff0_17, ff1_17, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * dg_26[k]
                  + f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_26[k];
    }
}

auto
compute_prim_fh_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg, const size_t dh,
                                     const size_t ff0, const size_t ff1, const size_t fg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 2.5 / p;

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
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_40 = buffer.data(dg + 40);
    const auto *dg_41 = buffer.data(dg + 41);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_18 = buffer.data(ff0 + 18);
    const auto *ff0_19 = buffer.data(ff0 + 19);
    const auto *ff0_20 = buffer.data(ff0 + 20);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_20 = buffer.data(ff1 + 20);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);
    const auto *ff1_35 = buffer.data(ff1 + 35);
    const auto *ff1_37 = buffer.data(ff1 + 37);
    const auto *ff1_38 = buffer.data(ff1 + 38);
    const auto *ff1_39 = buffer.data(ff1 + 39);
    const auto *ff1_40 = buffer.data(ff1 + 40);
    const auto *ff1_41 = buffer.data(ff1 + 41);

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
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, dg_0, ff0_0, ff0_1, ff1_0, \
                         ff1_1, fg_0, fg_1, fg_2, fg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_2[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];

        t_3[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, dg_5, dg_8, ff0_2, ff0_3, \
                         ff1_2, ff1_3, fg_4, fg_5, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_5[k] = f_0 * dg_5[k]
                 + pb_x[k] * fg_5[k];

        t_6[k] = f_0 * dg_8[k]
                 + pb_x[k] * fg_8[k];

        t_7[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, ff1_5, \
                         ff1_6, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_5[k]
                 + pb_y[k] * fg_6[k];

        t_9[k] = f_3 * ff0_5[k]
                 - f_4 * ff1_6[k]
                 + pb_y[k] * fg_7[k];

        t_10[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_6[k]
                  + pb_z[k] * fg_8[k];

        t_11[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_y, pb_x, pb_y, dg_0, dg_1, dg_3, dg_11, \
                         dh_1, dh_3, fg_9, fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * dg_0[k]
                  + pb_y[k] * fg_9[k];

        t_13[k] = f_8 * dg_1[k]
                  + pa_y[k] * dh_1[k];

        t_14[k] = f_0 * dg_3[k]
                  + pa_y[k] * dh_3[k];

        t_15[k] = f_8 * dg_11[k]
                  + pb_x[k] * fg_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_y, pa_z, pb_z, dg_0, dg_2, dg_5, dh_0, \
                         dh_2, dh_5, fg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_9 * dg_5[k]
                  + pa_y[k] * dh_5[k];

        t_17[k] = pa_z[k] * dh_0[k];

        t_18[k] = f_7 * dg_0[k]
                  + pb_z[k] * fg_12[k];

        t_19[k] = f_8 * dg_2[k]
                  + pa_z[k] * dh_2[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_z, pb_x, dg_4, dg_6, dg_7, dg_15, dh_4, \
                         dh_6, dh_7, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * dg_4[k]
                  + pa_z[k] * dh_4[k];

        t_21[k] = f_8 * dg_15[k]
                  + pb_x[k] * fg_16[k];

        t_22[k] = f_8 * dg_6[k]
                  + pa_z[k] * dh_6[k];

        t_23[k] = f_0 * dg_7[k]
                  + pa_z[k] * dh_7[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_x, pa_z, pb_y, dg_8, dg_9, dg_16, dg_18, \
                         dh_8, dh_9, dh_10, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_9 * dg_8[k]
                  + pa_z[k] * dh_8[k];

        t_25[k] = f_9 * dg_16[k]
                  + pa_x[k] * dh_9[k];

        t_26[k] = f_8 * dg_9[k]
                  + pb_y[k] * fg_17[k];

        t_27[k] = f_0 * dg_18[k]
                  + pa_x[k] * dh_10[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_x, dg_20, dg_21, dg_31, dh_12, \
                         dh_14, dh_18, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_8 * dg_20[k]
                  + pa_x[k] * dh_12[k];

        t_29[k] = f_7 * dg_21[k]
                  + pb_x[k] * fg_19[k];

        t_30[k] = pa_x[k] * dh_14[k];

        t_31[k] = f_9 * dg_31[k]
                  + pa_x[k] * dh_18[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, dg_12, dg_35, dg_36, dg_41, \
                         dh_20, dh_22, fg_20, fg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_8 * dg_12[k]
                  + pb_z[k] * fg_20[k];

        t_33[k] = f_0 * dg_35[k]
                  + pa_x[k] * dh_20[k];

        t_34[k] = f_8 * dg_36[k]
                  + pa_x[k] * dh_22[k];

        t_35[k] = f_7 * dg_41[k]
                  + pb_x[k] * fg_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pb_x, pb_y, dg_16, dh_26, ff0_10, \
                         ff0_11, ff1_20, ff1_22, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_x[k] * dh_26[k];

        t_37[k] = f_1 * ff0_10[k]
                  - f_2 * ff1_20[k]
                  + pb_x[k] * fg_24[k];

        t_38[k] = f_0 * dg_16[k]
                  + pb_y[k] * fg_24[k];

        t_39[k] = f_5 * ff0_11[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_25[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_x, ff0_12, ff0_13, ff0_15, ff1_23, ff1_24, \
                         ff1_26, fg_26, fg_27, fg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * ff0_12[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_26[k];

        t_41[k] = f_3 * ff0_13[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_27[k];

        t_42[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_28[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pb_y, pb_z, dg_21, dg_25, ff0_13, ff0_14, \
                         ff1_24, ff1_25, fg_29, fg_30, fg_31, fg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_0 * dg_21[k]
                  + f_1 * ff0_13[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_29[k];

        t_44[k] = f_3 * ff0_13[k]
                  - f_4 * ff1_24[k]
                  + pb_z[k] * fg_30[k];

        t_45[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_25[k]
                  + pb_z[k] * fg_31[k];

        t_46[k] = f_0 * dg_25[k]
                  + pb_y[k] * fg_33[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_z, pb_z, dg_17, dg_19, dh_11, dh_13, \
                         dh_14, ff0_15, ff1_26, fg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_1 * ff0_15[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_33[k];

        t_48[k] = f_8 * dg_17[k]
                  + pa_z[k] * dh_11[k];

        t_49[k] = f_0 * dg_19[k]
                  + pa_z[k] * dh_13[k];

        t_50[k] = pa_z[k] * dh_14[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_y, pb_z, dg_21, dg_22, dg_23, dg_30, \
                         dh_15, dh_16, fg_35, fg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_7 * dg_21[k]
                  + pb_z[k] * fg_35[k];

        t_52[k] = f_8 * dg_22[k]
                  + pa_z[k] * dh_15[k];

        t_53[k] = f_0 * dg_23[k]
                  + pa_z[k] * dh_16[k];

        t_54[k] = f_8 * dg_30[k]
                  + pb_y[k] * fg_39[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pa_z, dg_25, dg_32, dg_34, dg_37, \
                         dh_17, dh_19, dh_21, dh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_9 * dg_25[k]
                  + pa_z[k] * dh_17[k];

        t_56[k] = f_8 * dg_32[k]
                  + pa_y[k] * dh_19[k];

        t_57[k] = f_0 * dg_34[k]
                  + pa_y[k] * dh_21[k];

        t_58[k] = f_9 * dg_37[k]
                  + pa_y[k] * dh_23[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pb_y, pb_z, dg_26, dg_39, dg_40, dg_41, \
                         dh_24, dh_25, fg_41, fg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_8 * dg_26[k]
                  + pb_z[k] * fg_41[k];

        t_60[k] = f_0 * dg_39[k]
                  + pa_y[k] * dh_24[k];

        t_61[k] = f_8 * dg_40[k]
                  + pa_y[k] * dh_25[k];

        t_62[k] = f_7 * dg_41[k]
                  + pb_y[k] * fg_45[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_y, pb_x, pb_z, dg_31, dh_26, ff0_18, \
                         ff0_19, ff1_35, ff1_37, fg_46, fg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_y[k] * dh_26[k];

        t_64[k] = f_1 * ff0_18[k]
                  - f_2 * ff1_35[k]
                  + pb_x[k] * fg_46[k];

        t_65[k] = f_0 * dg_31[k]
                  + pb_z[k] * fg_46[k];

        t_66[k] = f_5 * ff0_19[k]
                  - f_6 * ff1_37[k]
                  + pb_x[k] * fg_48[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pb_x, pb_y, ff0_20, ff0_21, ff0_23, ff1_38, \
                         ff1_39, ff1_41, fg_49, fg_50, fg_51, fg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_5 * ff0_20[k]
                  - f_6 * ff1_38[k]
                  + pb_x[k] * fg_49[k];

        t_68[k] = f_3 * ff0_21[k]
                  - f_4 * ff1_39[k]
                  + pb_x[k] * fg_50[k];

        t_69[k] = f_3 * ff0_23[k]
                  - f_4 * ff1_41[k]
                  + pb_x[k] * fg_51[k];

        t_70[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_39[k]
                  + pb_y[k] * fg_52[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pb_y, pb_z, dg_37, dg_41, ff0_22, ff0_23, \
                         ff1_40, ff1_41, fg_52, fg_54, fg_55, fg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_0 * dg_37[k]
                  + pb_z[k] * fg_52[k];

        t_72[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_40[k]
                  + pb_y[k] * fg_54[k];

        t_73[k] = f_3 * ff0_23[k]
                  - f_4 * ff1_41[k]
                  + pb_y[k] * fg_55[k];

        t_74[k] = f_0 * dg_41[k]
                  + f_1 * ff0_23[k]
                  - f_2 * ff1_41[k]
                  + pb_z[k] * fg_56[k];
    }
}

auto
compute_prim_fh_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg, const size_t dh,
                                     const size_t ff0, const size_t ff1, const size_t fg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.0 / p;
    const auto f_8 = 2.5 / p;
    const auto f_9 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_37 = buffer.data(dh + 37);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_45 = buffer.data(dh + 45);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_49 = buffer.data(dh + 49);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_54 = buffer.data(dh + 54);
    const auto *dh_56 = buffer.data(dh + 56);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_18 = buffer.data(ff0 + 18);
    const auto *ff0_19 = buffer.data(ff0 + 19);
    const auto *ff0_20 = buffer.data(ff0 + 20);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, \
                         ff1_2, ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_z[k] * fg_3[k];

        t_7[k] = pb_y[k] * fg_4[k];

        t_8[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_9[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pb_y, pb_z, ff0_4, ff0_5, ff1_4, ff1_5, \
                         fg_5, fg_7, fg_8, fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_z[k] * fg_5[k];

        t_11[k] = f_5 * ff0_4[k]
                  - f_6 * ff1_4[k]
                  + pb_y[k] * fg_7[k];

        t_12[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_8[k];

        t_13[k] = pb_y[k] * fg_9[k];

        t_14[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, t_20, pa_y, dg_1, dg_3, dg_5, dh_0, \
                         dh_3, dh_4, dh_5, dh_7, dh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * dh_0[k];

        t_16[k] = f_7 * dg_1[k]
                  + pa_y[k] * dh_3[k];

        t_17[k] = pa_y[k] * dh_4[k];

        t_18[k] = f_0 * dg_3[k]
                  + pa_y[k] * dh_5[k];

        t_19[k] = pa_y[k] * dh_7[k];

        t_20[k] = f_8 * dg_5[k]
                  + pa_y[k] * dh_8[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pa_z, pb_y, dg_7, dg_8, dg_9, \
                         dh_0, dh_10, dh_11, dh_12, fg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * dg_7[k]
                  + pa_y[k] * dh_10[k];

        t_22[k] = f_7 * dg_8[k]
                  + pa_y[k] * dh_11[k];

        t_23[k] = f_9 * dg_9[k]
                  + pb_y[k] * fg_12[k];

        t_24[k] = pa_y[k] * dh_12[k];

        t_25[k] = pa_z[k] * dh_0[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_z, pb_y, pb_z, dg_0, dg_2, dh_3, \
                         dh_4, dh_5, fg_13, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_9 * dg_0[k]
                  + pb_z[k] * fg_13[k];

        t_27[k] = pa_z[k] * dh_3[k];

        t_28[k] = f_7 * dg_2[k]
                  + pa_z[k] * dh_4[k];

        t_29[k] = pa_z[k] * dh_5[k];

        t_30[k] = pb_y[k] * fg_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_z, pb_z, dg_4, dg_5, dg_6, dg_7, \
                         dh_7, dh_8, dh_10, dh_11, fg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * dg_4[k]
                  + pa_z[k] * dh_7[k];

        t_32[k] = pa_z[k] * dh_8[k];

        t_33[k] = f_9 * dg_5[k]
                  + pb_z[k] * fg_15[k];

        t_34[k] = f_7 * dg_6[k]
                  + pa_z[k] * dh_10[k];

        t_35[k] = f_0 * dg_7[k]
                  + pa_z[k] * dh_11[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pa_z, pb_y, dg_9, dg_11, dg_13, dh_12, \
                         dh_22, dh_24, fg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_y[k] * fg_18[k];

        t_37[k] = f_8 * dg_9[k]
                  + pa_z[k] * dh_12[k];

        t_38[k] = f_8 * dg_11[k]
                  + pa_x[k] * dh_22[k];

        t_39[k] = f_0 * dg_13[k]
                  + pa_x[k] * dh_24[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pa_x, pb_x, dg_14, dg_15, dg_16, dg_17, \
                         dh_25, dh_26, dh_28, dh_30, fg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * dg_14[k]
                  + pa_x[k] * dh_25[k];

        t_41[k] = f_7 * dg_15[k]
                  + pa_x[k] * dh_26[k];

        t_42[k] = f_7 * dg_16[k]
                  + pa_x[k] * dh_28[k];

        t_43[k] = f_9 * dg_17[k]
                  + pb_x[k] * fg_22[k];

        t_44[k] = pa_x[k] * dh_30[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, pa_x, pa_y, pa_z, dh_14, dh_18, \
                         dh_32, dh_33, dh_34, dh_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_x[k] * dh_32[k];

        t_46[k] = pa_x[k] * dh_33[k];

        t_47[k] = pa_x[k] * dh_34[k];

        t_48[k] = pa_x[k] * dh_35[k];

        t_49[k] = pa_y[k] * dh_18[k];

        t_50[k] = pa_z[k] * dh_14[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, pa_x, pa_y, pa_z, dh_15, dh_19, \
                         dh_20, dh_37, dh_38, dh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_y[k] * dh_19[k];

        t_52[k] = pa_z[k] * dh_15[k];

        t_53[k] = pa_y[k] * dh_20[k];

        t_54[k] = pa_x[k] * dh_37[k];

        t_55[k] = pa_x[k] * dh_38[k];

        t_56[k] = pa_x[k] * dh_39[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pa_x, pb_z, dg_10, dg_23, dg_25, dg_26, \
                         dh_40, dh_42, dh_45, dh_46, fg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_x[k] * dh_40[k];

        t_58[k] = f_8 * dg_23[k]
                  + pa_x[k] * dh_42[k];

        t_59[k] = f_7 * dg_10[k]
                  + pb_z[k] * fg_23[k];

        t_60[k] = f_0 * dg_25[k]
                  + pa_x[k] * dh_45[k];

        t_61[k] = f_0 * dg_26[k]
                  + pa_x[k] * dh_46[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, pa_x, pb_x, dg_27, dg_28, dg_32, dh_47, \
                         dh_49, dh_51, dh_52, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_7 * dg_27[k]
                  + pa_x[k] * dh_47[k];

        t_63[k] = f_7 * dg_28[k]
                  + pa_x[k] * dh_49[k];

        t_64[k] = f_9 * dg_32[k]
                  + pb_x[k] * fg_26[k];

        t_65[k] = pa_x[k] * dh_51[k];

        t_66[k] = pa_x[k] * dh_52[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, pa_x, pb_x, pb_z, dh_53, dh_54, dh_56, \
                         ff0_10, ff1_12, fg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pa_x[k] * dh_53[k];

        t_68[k] = pa_x[k] * dh_54[k];

        t_69[k] = pa_x[k] * dh_56[k];

        t_70[k] = f_1 * ff0_10[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_27[k];

        t_71[k] = pb_z[k] * fg_27[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pb_x, pb_z, ff0_11, ff0_12, ff0_13, ff1_13, \
                         ff1_14, ff1_15, fg_29, fg_30, fg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_5 * ff0_11[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_29[k];

        t_73[k] = f_5 * ff0_12[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_30[k];

        t_74[k] = f_3 * ff0_13[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_31[k];

        t_75[k] = pb_z[k] * fg_29[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pb_x, pb_y, dg_17, ff0_13, ff0_15, \
                         ff1_15, ff1_17, fg_32, fg_33, fg_35, fg_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_32[k];

        t_77[k] = pb_x[k] * fg_33[k];

        t_78[k] = pb_x[k] * fg_35[k];

        t_79[k] = pb_x[k] * fg_36[k];

        t_80[k] = f_0 * dg_17[k]
                  + f_1 * ff0_13[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_33[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pb_y, pb_z, dg_20, ff0_13, ff0_14, ff1_15, \
                         ff1_16, fg_33, fg_34, fg_35, fg_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pb_z[k] * fg_33[k];

        t_82[k] = f_3 * ff0_13[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_34[k];

        t_83[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_35[k];

        t_84[k] = f_0 * dg_20[k]
                  + pb_y[k] * fg_36[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pa_z, pb_z, dg_11, dg_12, dh_22, dh_24, \
                         dh_25, ff0_15, ff1_17, fg_36, fg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_1 * ff0_15[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_36[k];

        t_86[k] = pa_z[k] * dh_22[k];

        t_87[k] = f_9 * dg_11[k]
                  + pb_z[k] * fg_37[k];

        t_88[k] = pa_z[k] * dh_24[k];

        t_89[k] = f_7 * dg_12[k]
                  + pa_z[k] * dh_25[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pa_z, pb_x, dg_14, dh_26, dh_28, dh_30, \
                         fg_41, fg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pa_z[k] * dh_26[k];

        t_91[k] = f_0 * dg_14[k]
                  + pa_z[k] * dh_28[k];

        t_92[k] = pb_x[k] * fg_41[k];

        t_93[k] = pb_x[k] * fg_42[k];

        t_94[k] = pa_z[k] * dh_30[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_z, pb_y, pb_z, dg_17, dg_18, dg_19, dg_22, \
                         dh_32, dh_33, fg_40, fg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_9 * dg_17[k]
                  + pb_z[k] * fg_40[k];

        t_96[k] = f_7 * dg_18[k]
                  + pa_z[k] * dh_32[k];

        t_97[k] = f_0 * dg_19[k]
                  + pa_z[k] * dh_33[k];

        t_98[k] = f_7 * dg_22[k]
                  + pb_y[k] * fg_42[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pa_y, pa_z, dg_20, dg_24, dh_35, \
                         dh_42, dh_44, dh_45, dh_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_8 * dg_20[k]
                  + pa_z[k] * dh_35[k];

        t_100[k] = pa_y[k] * dh_42[k];

        t_101[k] = pa_y[k] * dh_44[k];

        t_102[k] = f_7 * dg_24[k]
                   + pa_y[k] * dh_45[k];

        t_103[k] = pa_y[k] * dh_46[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pa_y, pb_x, dg_25, dg_29, dh_47, \
                         dh_49, dh_51, fg_45, fg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_0 * dg_25[k]
                   + pa_y[k] * dh_47[k];

        t_105[k] = pa_y[k] * dh_49[k];

        t_106[k] = pb_x[k] * fg_45[k];

        t_107[k] = pb_x[k] * fg_46[k];

        t_108[k] = f_8 * dg_29[k]
                   + pa_y[k] * dh_51[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_y, pb_y, pb_z, dg_21, dg_30, dg_31, \
                         dg_32, dh_53, dh_54, fg_45, fg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_7 * dg_21[k]
                   + pb_z[k] * fg_45[k];

        t_110[k] = f_0 * dg_30[k]
                   + pa_y[k] * dh_53[k];

        t_111[k] = f_7 * dg_31[k]
                   + pa_y[k] * dh_54[k];

        t_112[k] = f_9 * dg_32[k]
                   + pb_y[k] * fg_48[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_y, pb_x, pb_y, pb_z, dg_23, dh_56, \
                         ff0_18, ff1_21, fg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_y[k] * dh_56[k];

        t_114[k] = f_1 * ff0_18[k]
                   - f_2 * ff1_21[k]
                   + pb_x[k] * fg_49[k];

        t_115[k] = pb_y[k] * fg_49[k];

        t_116[k] = f_0 * dg_23[k]
                   + pb_z[k] * fg_49[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_x, pb_y, ff0_19, ff0_20, ff0_21, \
                         ff1_22, ff1_23, ff1_24, fg_51, fg_52, fg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_5 * ff0_19[k]
                   - f_6 * ff1_22[k]
                   + pb_x[k] * fg_51[k];

        t_118[k] = f_5 * ff0_20[k]
                   - f_6 * ff1_23[k]
                   + pb_x[k] * fg_52[k];

        t_119[k] = f_3 * ff0_21[k]
                   - f_4 * ff1_24[k]
                   + pb_x[k] * fg_53[k];

        t_120[k] = pb_y[k] * fg_52[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pb_x, pb_y, ff0_21, ff0_23, \
                         ff1_24, ff1_26, fg_54, fg_55, fg_56, fg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_3 * ff0_23[k]
                   - f_4 * ff1_26[k]
                   + pb_x[k] * fg_54[k];

        t_122[k] = pb_x[k] * fg_55[k];

        t_123[k] = pb_x[k] * fg_56[k];

        t_124[k] = pb_x[k] * fg_58[k];

        t_125[k] = f_1 * ff0_21[k]
                   - f_2 * ff1_24[k]
                   + pb_y[k] * fg_55[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pb_y, pb_z, dg_29, ff0_22, ff0_23, \
                         ff1_25, ff1_26, fg_55, fg_56, fg_57, fg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_0 * dg_29[k]
                   + pb_z[k] * fg_55[k];

        t_127[k] = f_5 * ff0_22[k]
                   - f_6 * ff1_25[k]
                   + pb_y[k] * fg_56[k];

        t_128[k] = f_3 * ff0_23[k]
                   - f_4 * ff1_26[k]
                   + pb_y[k] * fg_57[k];

        t_129[k] = pb_y[k] * fg_58[k];
    }

#pragma omp simd aligned(t_130, pb_z, dg_32, ff0_23, ff1_26, fg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_0 * dg_32[k]
                   + f_1 * ff0_23[k]
                   - f_2 * ff1_26[k]
                   + pb_z[k] * fg_58[k];
    }
}

auto
compute_prim_fh_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg, const size_t dh,
                                     const size_t ff0, const size_t ff1, const size_t fg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.0 / p;
    const auto f_8 = 2.5 / p;
    const auto f_9 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_38 = buffer.data(dh + 38);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_18 = buffer.data(ff1 + 18);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);
    const auto *ff1_27 = buffer.data(ff1 + 27);
    const auto *ff1_28 = buffer.data(ff1 + 28);
    const auto *ff1_29 = buffer.data(ff1 + 29);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_15 = buffer.data(fg + 15);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, \
                         ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_10[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = pb_y[k] * fg_8[k];

        t_12[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_13[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pa_z, dg_1, dg_3, dg_5, dh_0, dh_3, \
                         dh_5, dh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_7 * dg_1[k]
                  + pa_y[k] * dh_3[k];

        t_15[k] = f_0 * dg_3[k]
                  + pa_y[k] * dh_5[k];

        t_16[k] = f_8 * dg_5[k]
                  + pa_y[k] * dh_8[k];

        t_17[k] = pa_z[k] * dh_0[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_z, pb_z, dg_0, dg_2, dg_4, dg_6, dh_4, \
                         dh_7, dh_9, fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_9 * dg_0[k]
                  + pb_z[k] * fg_11[k];

        t_19[k] = f_7 * dg_2[k]
                  + pa_z[k] * dh_4[k];

        t_20[k] = f_0 * dg_4[k]
                  + pa_z[k] * dh_7[k];

        t_21[k] = f_7 * dg_6[k]
                  + pa_z[k] * dh_9[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pa_z, dg_7, dg_9, dg_11, dg_13, dh_10, \
                         dh_12, dh_13, dh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_0 * dg_7[k]
                  + pa_z[k] * dh_10[k];

        t_23[k] = f_8 * dg_9[k]
                  + pa_z[k] * dh_12[k];

        t_24[k] = f_8 * dg_11[k]
                  + pa_x[k] * dh_13[k];

        t_25[k] = f_0 * dg_13[k]
                  + pa_x[k] * dh_14[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_x, pb_z, dg_10, dg_15, dg_23, dg_26, \
                         dh_16, dh_21, dh_26, dh_28, fg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * dg_15[k]
                  + pa_x[k] * dh_16[k];

        t_27[k] = pa_x[k] * dh_21[k];

        t_28[k] = f_8 * dg_23[k]
                  + pa_x[k] * dh_26[k];

        t_29[k] = f_7 * dg_10[k]
                  + pb_z[k] * fg_15[k];

        t_30[k] = f_0 * dg_26[k]
                  + pa_x[k] * dh_28[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pb_x, dg_28, dh_30, dh_38, ff0_12, \
                         ff0_13, ff1_13, ff1_14, fg_17, fg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_7 * dg_28[k]
                  + pa_x[k] * dh_30[k];

        t_32[k] = pa_x[k] * dh_38[k];

        t_33[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_13[k]
                  + pb_x[k] * fg_17[k];

        t_34[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_18[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, ff0_14, ff0_15, ff0_17, ff1_15, ff1_16, \
                         ff1_18, fg_19, fg_20, fg_21, fg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_15[k]
                  + pb_x[k] * fg_19[k];

        t_36[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_16[k]
                  + pb_x[k] * fg_20[k];

        t_37[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_18[k]
                  + pb_x[k] * fg_21[k];

        t_38[k] = pb_x[k] * fg_22[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pb_x, pb_y, pb_z, dg_17, ff0_15, \
                         ff1_16, fg_22, fg_23, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pb_x[k] * fg_24[k];

        t_40[k] = pb_x[k] * fg_25[k];

        t_41[k] = f_0 * dg_17[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_16[k]
                  + pb_y[k] * fg_22[k];

        t_42[k] = pb_z[k] * fg_22[k];

        t_43[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_16[k]
                  + pb_z[k] * fg_23[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pb_y, pb_z, dg_12, dg_20, dh_15, \
                         ff0_16, ff0_17, ff1_17, ff1_18, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_17[k]
                  + pb_z[k] * fg_24[k];

        t_45[k] = f_0 * dg_20[k]
                  + pb_y[k] * fg_25[k];

        t_46[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_18[k]
                  + pb_z[k] * fg_25[k];

        t_47[k] = f_7 * dg_12[k]
                  + pa_z[k] * dh_15[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_z, pb_z, dg_14, dg_17, dg_18, dg_19, \
                         dh_17, dh_21, dh_23, dh_24, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * dg_14[k]
                  + pa_z[k] * dh_17[k];

        t_49[k] = pa_z[k] * dh_21[k];

        t_50[k] = f_9 * dg_17[k]
                  + pb_z[k] * fg_26[k];

        t_51[k] = f_7 * dg_18[k]
                  + pa_z[k] * dh_23[k];

        t_52[k] = f_0 * dg_19[k]
                  + pa_z[k] * dh_24[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_y, pa_z, pb_y, dg_20, dg_22, dg_24, dg_25, \
                         dh_25, dh_27, dh_29, fg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_7 * dg_22[k]
                  + pb_y[k] * fg_27[k];

        t_54[k] = f_8 * dg_20[k]
                  + pa_z[k] * dh_25[k];

        t_55[k] = f_7 * dg_24[k]
                  + pa_y[k] * dh_27[k];

        t_56[k] = f_0 * dg_25[k]
                  + pa_y[k] * dh_29[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_y, pb_z, dg_21, dg_29, dg_30, dg_31, \
                         dh_34, dh_35, dh_36, fg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_8 * dg_29[k]
                  + pa_y[k] * dh_34[k];

        t_58[k] = f_7 * dg_21[k]
                  + pb_z[k] * fg_28[k];

        t_59[k] = f_0 * dg_30[k]
                  + pa_y[k] * dh_35[k];

        t_60[k] = f_7 * dg_31[k]
                  + pa_y[k] * dh_36[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_y, pb_x, pb_y, pb_z, dg_23, dg_32, dh_38, \
                         ff0_21, ff1_24, fg_29, fg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_9 * dg_32[k]
                  + pb_y[k] * fg_29[k];

        t_62[k] = pa_y[k] * dh_38[k];

        t_63[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_24[k]
                  + pb_x[k] * fg_30[k];

        t_64[k] = f_0 * dg_23[k]
                  + pb_z[k] * fg_30[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pb_x, ff0_22, ff0_23, ff0_24, ff1_25, ff1_26, \
                         ff1_27, fg_31, fg_32, fg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_25[k]
                  + pb_x[k] * fg_31[k];

        t_66[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_26[k]
                  + pb_x[k] * fg_32[k];

        t_67[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_27[k]
                  + pb_x[k] * fg_33[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pb_x, pb_y, ff0_24, ff0_26, ff1_27, \
                         ff1_29, fg_34, fg_35, fg_36, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_29[k]
                  + pb_x[k] * fg_34[k];

        t_69[k] = pb_x[k] * fg_35[k];

        t_70[k] = pb_x[k] * fg_36[k];

        t_71[k] = pb_x[k] * fg_38[k];

        t_72[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_27[k]
                  + pb_y[k] * fg_35[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pb_y, pb_z, dg_29, ff0_25, ff0_26, ff1_28, \
                         ff1_29, fg_35, fg_36, fg_37, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_0 * dg_29[k]
                  + pb_z[k] * fg_35[k];

        t_74[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_28[k]
                  + pb_y[k] * fg_36[k];

        t_75[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_29[k]
                  + pb_y[k] * fg_37[k];

        t_76[k] = pb_y[k] * fg_38[k];
    }

#pragma omp simd aligned(t_77, pb_z, dg_32, ff0_26, ff1_29, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_0 * dg_32[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_29[k]
                  + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg, const size_t dh,
                                     const size_t ff0, const size_t ff1, const size_t fg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_47 = buffer.data(dh + 47);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_18 = buffer.data(ff0 + 18);
    const auto *ff0_19 = buffer.data(ff0 + 19);
    const auto *ff0_20 = buffer.data(ff0 + 20);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, \
                         ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_10[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = pb_y[k] * fg_8[k];

        t_12[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_13[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_z, pb_x, dh_0, dh_26, dh_47, ff0_10, \
                         ff1_12, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * dh_0[k];

        t_15[k] = pa_x[k] * dh_26[k];

        t_16[k] = pa_x[k] * dh_47[k];

        t_17[k] = f_1 * ff0_10[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_x, ff0_11, ff0_12, ff0_13, ff1_13, ff1_14, \
                         ff1_15, fg_18, fg_19, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * ff0_11[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_18[k];

        t_19[k] = f_5 * ff0_12[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_19[k];

        t_20[k] = f_3 * ff0_13[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_y, dg_14, ff0_13, ff0_15, \
                         ff1_15, ff1_17, fg_21, fg_22, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_21[k];

        t_22[k] = pb_x[k] * fg_22[k];

        t_23[k] = pb_x[k] * fg_24[k];

        t_24[k] = pb_x[k] * fg_25[k];

        t_25[k] = f_0 * dg_14[k]
                  + f_1 * ff0_13[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_22[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, ff0_13, ff0_14, ff0_15, ff1_15, ff1_16, \
                         ff1_17, fg_22, fg_23, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * fg_22[k];

        t_27[k] = f_3 * ff0_13[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_23[k];

        t_28[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_24[k];

        t_29[k] = f_1 * ff0_15[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_25[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pa_z, pb_x, dh_26, dh_47, ff0_18, \
                         ff0_19, ff1_21, ff1_22, fg_30, fg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_z[k] * dh_26[k];

        t_31[k] = pa_y[k] * dh_47[k];

        t_32[k] = f_1 * ff0_18[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_30[k];

        t_33[k] = f_5 * ff0_19[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_31[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_x, ff0_20, ff0_21, ff0_23, ff1_23, ff1_24, \
                         ff1_26, fg_32, fg_33, fg_34, fg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_5 * ff0_20[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_32[k];

        t_35[k] = f_3 * ff0_21[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_33[k];

        t_36[k] = f_3 * ff0_23[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_34[k];

        t_37[k] = pb_x[k] * fg_35[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_y, ff0_21, ff0_22, ff1_24, ff1_25, \
                         fg_35, fg_36, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_x[k] * fg_36[k];

        t_39[k] = pb_x[k] * fg_38[k];

        t_40[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_35[k];

        t_41[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_36[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_y, pb_z, dg_26, ff0_23, ff1_26, fg_37, \
                         fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * ff0_23[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_37[k];

        t_43[k] = pb_y[k] * fg_38[k];

        t_44[k] = f_0 * dg_26[k]
                  + f_1 * ff0_23[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg, const size_t dh,
                                     const size_t ff0, const size_t ff1, const size_t fg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 2.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 0.5 / p;

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
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_41 = buffer.data(dh + 41);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_44 = buffer.data(dh + 44);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_14 = buffer.data(fg + 14);
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
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, \
                         ff1_2, ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];

        t_9[k] = pb_z[k] * fg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * ff0_4[k]
                  - f_6 * ff1_4[k]
                  + pb_y[k] * fg_6[k];

        t_11[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_12[k] = pb_y[k] * fg_8[k];

        t_13[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_14[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, t_20, pa_y, pa_z, dg_2, dg_4, dg_5, \
                         dh_0, dh_4, dh_7, dh_8, dh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_7 * dg_5[k]
                  + pa_y[k] * dh_8[k];

        t_16[k] = pa_y[k] * dh_13[k];

        t_17[k] = pa_z[k] * dh_0[k];

        t_18[k] = f_8 * dg_2[k]
                  + pa_z[k] * dh_4[k];

        t_19[k] = f_0 * dg_4[k]
                  + pa_z[k] * dh_7[k];

        t_20[k] = pa_z[k] * dh_8[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_z, pb_y, dg_6, dg_7, dg_9, dh_10, dh_11, \
                         dh_13, fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * dg_6[k]
                  + pa_z[k] * dh_10[k];

        t_22[k] = f_0 * dg_7[k]
                  + pa_z[k] * dh_11[k];

        t_23[k] = pb_y[k] * fg_11[k];

        t_24[k] = f_7 * dg_9[k]
                  + pa_z[k] * dh_13[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_x, pb_x, dg_10, dg_12, dg_14, dg_16, \
                         dh_15, dh_17, dh_19, dh_24, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * dg_10[k]
                  + pa_x[k] * dh_15[k];

        t_26[k] = f_0 * dg_12[k]
                  + pa_x[k] * dh_17[k];

        t_27[k] = f_8 * dg_14[k]
                  + pa_x[k] * dh_19[k];

        t_28[k] = f_9 * dg_16[k]
                  + pb_x[k] * fg_14[k];

        t_29[k] = pa_x[k] * dh_24[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pa_x, pb_x, dg_20, dg_23, dg_25, dg_29, \
                         dh_31, dh_34, dh_36, dh_44, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * dg_20[k]
                  + pa_x[k] * dh_31[k];

        t_31[k] = f_0 * dg_23[k]
                  + pa_x[k] * dh_34[k];

        t_32[k] = f_8 * dg_25[k]
                  + pa_x[k] * dh_36[k];

        t_33[k] = f_9 * dg_29[k]
                  + pb_x[k] * fg_17[k];

        t_34[k] = pa_x[k] * dh_44[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, pb_z, ff0_12, ff0_13, ff0_14, ff1_12, \
                         ff1_13, ff1_14, fg_18, fg_19, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_18[k];

        t_36[k] = pb_z[k] * fg_18[k];

        t_37[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_19[k];

        t_38[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_20[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pb_x, ff0_15, ff0_17, ff1_15, ff1_17, \
                         fg_21, fg_22, fg_23, fg_25, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_21[k];

        t_40[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_22[k];

        t_41[k] = pb_x[k] * fg_23[k];

        t_42[k] = pb_x[k] * fg_25[k];

        t_43[k] = pb_x[k] * fg_26[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_y, pb_z, dg_16, ff0_15, ff0_16, ff1_15, \
                         ff1_16, fg_23, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * dg_16[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_23[k];

        t_45[k] = pb_z[k] * fg_23[k];

        t_46[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_24[k];

        t_47[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_25[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_z, dg_11, dg_13, dh_15, dh_18, \
                         dh_20, ff0_17, ff1_17, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_26[k];

        t_49[k] = pa_z[k] * dh_15[k];

        t_50[k] = f_8 * dg_11[k]
                  + pa_z[k] * dh_18[k];

        t_51[k] = f_0 * dg_13[k]
                  + pa_z[k] * dh_20[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pa_z, pb_x, dg_17, dg_18, dg_19, dh_24, \
                         dh_26, dh_27, dh_28, fg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_x[k] * fg_28[k];

        t_53[k] = pa_z[k] * dh_24[k];

        t_54[k] = f_8 * dg_17[k]
                  + pa_z[k] * dh_26[k];

        t_55[k] = f_0 * dg_18[k]
                  + pa_z[k] * dh_27[k];

        t_56[k] = f_7 * dg_19[k]
                  + pa_z[k] * dh_28[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pa_y, pb_x, dg_21, dg_22, dg_26, dg_27, \
                         dh_33, dh_35, dh_40, dh_41, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_8 * dg_21[k]
                  + pa_y[k] * dh_33[k];

        t_58[k] = f_0 * dg_22[k]
                  + pa_y[k] * dh_35[k];

        t_59[k] = pb_x[k] * fg_29[k];

        t_60[k] = f_7 * dg_26[k]
                  + pa_y[k] * dh_40[k];

        t_61[k] = f_0 * dg_27[k]
                  + pa_y[k] * dh_41[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, pa_y, pb_x, pb_y, dg_28, dg_29, dh_42, \
                         dh_44, ff0_21, ff1_21, fg_32, fg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_8 * dg_28[k]
                  + pa_y[k] * dh_42[k];

        t_63[k] = f_9 * dg_29[k]
                  + pb_y[k] * fg_32[k];

        t_64[k] = pa_y[k] * dh_44[k];

        t_65[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_33[k];

        t_66[k] = pb_y[k] * fg_33[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pb_x, ff0_22, ff0_23, ff0_24, ff1_22, ff1_23, \
                         ff1_24, fg_34, fg_35, fg_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_34[k];

        t_68[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_35[k];

        t_69[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_36[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pb_x, pb_y, ff0_24, ff0_26, ff1_24, \
                         ff1_26, fg_37, fg_38, fg_39, fg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_37[k];

        t_71[k] = pb_x[k] * fg_38[k];

        t_72[k] = pb_x[k] * fg_39[k];

        t_73[k] = pb_x[k] * fg_41[k];

        t_74[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_38[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pb_y, pb_z, dg_29, ff0_25, ff0_26, ff1_25, \
                         ff1_26, fg_39, fg_40, fg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_39[k];

        t_76[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_40[k];

        t_77[k] = pb_y[k] * fg_41[k];

        t_78[k] = f_0 * dg_29[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_41[k];
    }
}

auto
compute_prim_fh_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dg, const size_t dh,
                                     const size_t ff0, const size_t ff1, const size_t fg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_38 = buffer.data(dh + 38);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, \
                         ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_10[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = pb_y[k] * fg_8[k];

        t_12[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_13[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_x, pa_y, pa_z, dg_5, dg_8, dg_9, \
                         dh_0, dh_8, dh_12, dh_13, dh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_7 * dg_5[k]
                  + pa_y[k] * dh_8[k];

        t_15[k] = pa_z[k] * dh_0[k];

        t_16[k] = f_7 * dg_8[k]
                  + pa_z[k] * dh_12[k];

        t_17[k] = f_7 * dg_9[k]
                  + pa_x[k] * dh_13[k];

        t_18[k] = pa_x[k] * dh_21[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, dg_18, dh_26, dh_38, ff0_12, \
                         ff0_13, ff1_12, ff1_13, fg_13, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * dg_18[k]
                  + pa_x[k] * dh_26[k];

        t_20[k] = pa_x[k] * dh_38[k];

        t_21[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_13[k];

        t_22[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_14[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_x, ff0_14, ff0_15, ff0_17, ff1_14, ff1_15, \
                         ff1_17, fg_15, fg_16, fg_17, fg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_15[k];

        t_24[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_16[k];

        t_25[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_17[k];

        t_26[k] = pb_x[k] * fg_18[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pb_x, pb_y, pb_z, dg_14, ff0_15, \
                         ff1_15, fg_18, fg_19, fg_20, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_x[k] * fg_20[k];

        t_28[k] = pb_x[k] * fg_21[k];

        t_29[k] = f_0 * dg_14[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_18[k];

        t_30[k] = pb_z[k] * fg_18[k];

        t_31[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_19[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_z, pb_z, dg_17, dh_21, dh_25, ff0_16, \
                         ff0_17, ff1_16, ff1_17, fg_20, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_20[k];

        t_33[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_21[k];

        t_34[k] = pa_z[k] * dh_21[k];

        t_35[k] = f_7 * dg_17[k]
                  + pa_z[k] * dh_25[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_y, pb_x, dg_23, dh_34, dh_38, ff0_21, \
                         ff0_22, ff1_21, ff1_22, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_7 * dg_23[k]
                  + pa_y[k] * dh_34[k];

        t_37[k] = pa_y[k] * dh_38[k];

        t_38[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_24[k];

        t_39[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_25[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_x, ff0_23, ff0_24, ff0_26, ff1_23, ff1_24, \
                         ff1_26, fg_26, fg_27, fg_28, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_26[k];

        t_41[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_27[k];

        t_42[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_28[k];

        t_43[k] = pb_x[k] * fg_29[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_y, ff0_24, ff0_25, ff1_24, ff1_25, \
                         fg_29, fg_30, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * fg_30[k];

        t_45[k] = pb_x[k] * fg_32[k];

        t_46[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_29[k];

        t_47[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_30[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_y, pb_z, dg_26, ff0_26, ff1_26, fg_31, \
                         fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_31[k];

        t_49[k] = pb_y[k] * fg_32[k];

        t_50[k] = f_0 * dg_26[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_32[k];
    }
}

auto
compute_prim_fh_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_38 = buffer.data(dh + 38);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_18 = buffer.data(ff0 + 18);
    const auto *ff0_19 = buffer.data(ff0 + 19);
    const auto *ff0_20 = buffer.data(ff0 + 20);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, \
                         ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, ff0_4, ff0_5, ff1_4, ff1_5, fg_6, \
                         fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_10[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = pb_y[k] * fg_8[k];

        t_12[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_x, pb_x, dg_10, dg_12, dg_14, dg_20, \
                         dh_14, dh_16, dh_21, dh_28, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dg_10[k]
                  + pa_x[k] * dh_14[k];

        t_14[k] = f_7 * dg_12[k]
                  + pa_x[k] * dh_16[k];

        t_15[k] = f_8 * dg_14[k]
                  + pb_x[k] * fg_13[k];

        t_16[k] = pa_x[k] * dh_21[k];

        t_17[k] = f_0 * dg_20[k]
                  + pa_x[k] * dh_28[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pb_x, dg_22, dg_26, dh_30, dh_38, \
                         ff0_10, ff1_12, fg_16, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_7 * dg_22[k]
                  + pa_x[k] * dh_30[k];

        t_19[k] = f_8 * dg_26[k]
                  + pb_x[k] * fg_16[k];

        t_20[k] = pa_x[k] * dh_38[k];

        t_21[k] = f_1 * ff0_10[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_17[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pb_x, ff0_11, ff0_12, ff0_13, ff1_13, ff1_14, \
                         ff1_15, fg_18, fg_19, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * ff0_11[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_18[k];

        t_23[k] = f_5 * ff0_12[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_19[k];

        t_24[k] = f_3 * ff0_13[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_20[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_x, pb_y, dg_14, ff0_13, ff0_15, \
                         ff1_15, ff1_17, fg_21, fg_22, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_21[k];

        t_26[k] = pb_x[k] * fg_22[k];

        t_27[k] = pb_x[k] * fg_24[k];

        t_28[k] = pb_x[k] * fg_25[k];

        t_29[k] = f_0 * dg_14[k]
                  + f_1 * ff0_13[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_22[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_z, ff0_13, ff0_14, ff0_15, ff1_15, ff1_16, \
                         ff1_17, fg_22, fg_23, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_z[k] * fg_22[k];

        t_31[k] = f_3 * ff0_13[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_23[k];

        t_32[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_24[k];

        t_33[k] = f_1 * ff0_15[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_25[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pb_y, dg_24, dg_25, dg_26, dh_35, \
                         dh_36, dh_38, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_0 * dg_24[k]
                  + pa_y[k] * dh_35[k];

        t_35[k] = f_7 * dg_25[k]
                  + pa_y[k] * dh_36[k];

        t_36[k] = f_8 * dg_26[k]
                  + pb_y[k] * fg_29[k];

        t_37[k] = pa_y[k] * dh_38[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pb_x, ff0_18, ff0_19, ff0_20, ff1_21, ff1_22, \
                         ff1_23, fg_30, fg_31, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_1 * ff0_18[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_30[k];

        t_39[k] = f_5 * ff0_19[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_31[k];

        t_40[k] = f_5 * ff0_20[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pb_x, ff0_21, ff0_23, ff1_24, ff1_26, \
                         fg_33, fg_34, fg_35, fg_36, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_3 * ff0_21[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_33[k];

        t_42[k] = f_3 * ff0_23[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_34[k];

        t_43[k] = pb_x[k] * fg_35[k];

        t_44[k] = pb_x[k] * fg_36[k];

        t_45[k] = pb_x[k] * fg_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_y, ff0_21, ff0_22, ff0_23, ff1_24, ff1_25, \
                         ff1_26, fg_35, fg_36, fg_37, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_35[k];

        t_47[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_36[k];

        t_48[k] = f_3 * ff0_23[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_37[k];

        t_49[k] = pb_y[k] * fg_38[k];
    }

#pragma omp simd aligned(t_50, pb_z, dg_26, ff0_23, ff1_26, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_0 * dg_26[k]
                  + f_1 * ff0_23[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_41 = buffer.data(dh + 41);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, \
                         ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_10[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = pb_y[k] * fg_8[k];

        t_12[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_13[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_x, pa_z, pb_x, dg_10, dg_12, dg_14, \
                         dh_0, dh_16, dh_18, dh_23, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * dh_0[k];

        t_15[k] = f_0 * dg_10[k]
                  + pa_x[k] * dh_16[k];

        t_16[k] = f_7 * dg_12[k]
                  + pa_x[k] * dh_18[k];

        t_17[k] = f_8 * dg_14[k]
                  + pb_x[k] * fg_13[k];

        t_18[k] = pa_x[k] * dh_23[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, dg_20, dg_22, dg_26, dh_31, \
                         dh_33, dh_41, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * dg_20[k]
                  + pa_x[k] * dh_31[k];

        t_20[k] = f_7 * dg_22[k]
                  + pa_x[k] * dh_33[k];

        t_21[k] = f_8 * dg_26[k]
                  + pb_x[k] * fg_16[k];

        t_22[k] = pa_x[k] * dh_41[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_x, ff0_12, ff0_13, ff0_14, ff1_12, ff1_13, \
                         ff1_14, fg_17, fg_18, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_17[k];

        t_24[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_18[k];

        t_25[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_19[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pb_x, ff0_15, ff0_17, ff1_15, ff1_17, \
                         fg_20, fg_21, fg_22, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_20[k];

        t_27[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_21[k];

        t_28[k] = pb_x[k] * fg_22[k];

        t_29[k] = pb_x[k] * fg_24[k];

        t_30[k] = pb_x[k] * fg_25[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pb_y, pb_z, dg_14, ff0_15, ff0_16, ff1_15, \
                         ff1_16, fg_22, fg_23, fg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * dg_14[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_22[k];

        t_32[k] = pb_z[k] * fg_22[k];

        t_33[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_23[k];

        t_34[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_24[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_y, pa_z, pb_z, dg_24, dg_25, dh_23, dh_38, \
                         dh_39, ff0_17, ff1_17, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_25[k];

        t_36[k] = pa_z[k] * dh_23[k];

        t_37[k] = f_0 * dg_24[k]
                  + pa_y[k] * dh_38[k];

        t_38[k] = f_7 * dg_25[k]
                  + pa_y[k] * dh_39[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_y, pb_x, pb_y, dg_26, dh_41, ff0_21, \
                         ff0_22, ff1_21, ff1_22, fg_29, fg_30, fg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_8 * dg_26[k]
                  + pb_y[k] * fg_29[k];

        t_40[k] = pa_y[k] * dh_41[k];

        t_41[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_30[k];

        t_42[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_31[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pb_x, ff0_23, ff0_24, ff0_26, ff1_23, ff1_24, \
                         ff1_26, fg_32, fg_33, fg_34, fg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_32[k];

        t_44[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_33[k];

        t_45[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_34[k];

        t_46[k] = pb_x[k] * fg_35[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_x, pb_y, ff0_24, ff0_25, ff1_24, ff1_25, \
                         fg_35, fg_36, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_x[k] * fg_36[k];

        t_48[k] = pb_x[k] * fg_38[k];

        t_49[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_35[k];

        t_50[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_36[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_y, pb_z, dg_26, ff0_26, ff1_26, fg_37, \
                         fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_37[k];

        t_52[k] = pb_y[k] * fg_38[k];

        t_53[k] = f_0 * dg_26[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 2.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_38 = buffer.data(dh + 38);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, \
                         ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_10[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = pb_y[k] * fg_8[k];

        t_12[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_13[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_z, dg_8, dg_10, dg_12, dh_0, dh_12, \
                         dh_14, dh_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * dh_0[k];

        t_15[k] = f_7 * dg_8[k]
                  + pa_z[k] * dh_12[k];

        t_16[k] = f_0 * dg_10[k]
                  + pa_x[k] * dh_14[k];

        t_17[k] = f_8 * dg_12[k]
                  + pa_x[k] * dh_16[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_x, pb_x, dg_14, dg_20, dg_22, dg_26, \
                         dh_21, dh_28, dh_30, fg_13, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_9 * dg_14[k]
                  + pb_x[k] * fg_13[k];

        t_19[k] = pa_x[k] * dh_21[k];

        t_20[k] = f_0 * dg_20[k]
                  + pa_x[k] * dh_28[k];

        t_21[k] = f_8 * dg_22[k]
                  + pa_x[k] * dh_30[k];

        t_22[k] = f_9 * dg_26[k]
                  + pb_x[k] * fg_16[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_x, dh_38, ff0_12, ff0_13, ff0_14, \
                         ff1_12, ff1_13, ff1_14, fg_17, fg_18, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_x[k] * dh_38[k];

        t_24[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_17[k];

        t_25[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_18[k];

        t_26[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_19[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pb_x, ff0_15, ff0_17, ff1_15, ff1_17, \
                         fg_20, fg_21, fg_22, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_20[k];

        t_28[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_21[k];

        t_29[k] = pb_x[k] * fg_22[k];

        t_30[k] = pb_x[k] * fg_24[k];

        t_31[k] = pb_x[k] * fg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_y, pb_z, dg_14, ff0_15, ff0_16, ff1_15, \
                         ff1_16, fg_22, fg_23, fg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * dg_14[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_22[k];

        t_33[k] = pb_z[k] * fg_22[k];

        t_34[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_23[k];

        t_35[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_24[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_y, pa_z, pb_z, dg_17, dg_23, dh_21, dh_25, \
                         dh_34, ff0_17, ff1_17, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_25[k];

        t_37[k] = pa_z[k] * dh_21[k];

        t_38[k] = f_7 * dg_17[k]
                  + pa_z[k] * dh_25[k];

        t_39[k] = f_7 * dg_23[k]
                  + pa_y[k] * dh_34[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pb_y, dg_24, dg_25, dg_26, dh_35, \
                         dh_36, dh_38, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * dg_24[k]
                  + pa_y[k] * dh_35[k];

        t_41[k] = f_8 * dg_25[k]
                  + pa_y[k] * dh_36[k];

        t_42[k] = f_9 * dg_26[k]
                  + pb_y[k] * fg_29[k];

        t_43[k] = pa_y[k] * dh_38[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pb_x, ff0_21, ff0_22, ff0_23, ff1_21, ff1_22, \
                         ff1_23, fg_30, fg_31, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_30[k];

        t_45[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_31[k];

        t_46[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pb_x, ff0_24, ff0_26, ff1_24, ff1_26, \
                         fg_33, fg_34, fg_35, fg_36, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_33[k];

        t_48[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_34[k];

        t_49[k] = pb_x[k] * fg_35[k];

        t_50[k] = pb_x[k] * fg_36[k];

        t_51[k] = pb_x[k] * fg_38[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_y, ff0_24, ff0_25, ff0_26, ff1_24, ff1_25, \
                         ff1_26, fg_35, fg_36, fg_37, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_35[k];

        t_53[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_36[k];

        t_54[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_37[k];

        t_55[k] = pb_y[k] * fg_38[k];
    }

#pragma omp simd aligned(t_56, pb_z, dg_26, ff0_26, ff1_26, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * dg_26[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_38 = buffer.data(dh + 38);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, \
                         ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_10[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = pb_y[k] * fg_8[k];

        t_12[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_13[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_z, pb_x, dh_0, dh_21, dh_38, ff0_12, \
                         ff1_12, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * dh_0[k];

        t_15[k] = pa_x[k] * dh_21[k];

        t_16[k] = pa_x[k] * dh_38[k];

        t_17[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_13[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_x, ff0_13, ff0_14, ff0_15, ff1_13, ff1_14, \
                         ff1_15, fg_14, fg_15, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_14[k];

        t_19[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_15[k];

        t_20[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_16[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_y, dg_14, ff0_15, ff0_17, \
                         ff1_15, ff1_17, fg_17, fg_18, fg_20, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_17[k];

        t_22[k] = pb_x[k] * fg_18[k];

        t_23[k] = pb_x[k] * fg_20[k];

        t_24[k] = pb_x[k] * fg_21[k];

        t_25[k] = f_0 * dg_14[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_18[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, ff0_15, ff0_16, ff0_17, ff1_15, ff1_16, \
                         ff1_17, fg_18, fg_19, fg_20, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * fg_18[k];

        t_27[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_19[k];

        t_28[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_20[k];

        t_29[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pa_z, pb_x, dh_21, dh_38, ff0_21, \
                         ff0_22, ff1_21, ff1_22, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_z[k] * dh_21[k];

        t_31[k] = pa_y[k] * dh_38[k];

        t_32[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_24[k];

        t_33[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_25[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_x, ff0_23, ff0_24, ff0_26, ff1_23, ff1_24, \
                         ff1_26, fg_26, fg_27, fg_28, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_26[k];

        t_35[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_27[k];

        t_36[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_28[k];

        t_37[k] = pb_x[k] * fg_29[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_y, ff0_24, ff0_25, ff1_24, ff1_25, \
                         fg_29, fg_30, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_x[k] * fg_30[k];

        t_39[k] = pb_x[k] * fg_32[k];

        t_40[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_29[k];

        t_41[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_y, pb_z, dg_26, ff0_26, ff1_26, fg_31, \
                         fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_31[k];

        t_43[k] = pb_y[k] * fg_32[k];

        t_44[k] = f_0 * dg_26[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_32[k];
    }
}

auto
compute_prim_fh_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_8 = buffer.data(ff0 + 8);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_8 = buffer.data(ff1 + 8);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_1, dh_1, dh_2, ff0_5, ff1_5, \
                         fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_1[k]
                 + f_1 * ff0_5[k]
                 - f_2 * ff1_5[k]
                 + pb_y[k] * fg_5[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_2, ff0_8, ff1_8, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_2[k]
                 + f_1 * ff0_8[k]
                 - f_2 * ff1_8[k]
                 + pb_z[k] * fg_8[k];
    }
}

auto
compute_prim_fh_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_8 = buffer.data(ff0 + 8);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_41 = buffer.data(ff1 + 41);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_58 = buffer.data(fg + 58);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_1, dh_1, dh_2, ff0_5, ff1_24, \
                         fg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_1[k]
                 + f_1 * ff0_5[k]
                 - f_2 * ff1_24[k]
                 + pb_y[k] * fg_31[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_2, ff0_8, ff1_41, fg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_2[k]
                 + f_1 * ff0_8[k]
                 - f_2 * ff1_41[k]
                 + pb_z[k] * fg_58[k];
    }
}

auto
compute_prim_fh_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_17 = buffer.data(dg + 17);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_20 = buffer.data(ff0 + 20);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);
    const auto *ff0_35 = buffer.data(ff0 + 35);
    const auto *ff0_37 = buffer.data(ff0 + 37);
    const auto *ff0_38 = buffer.data(ff0 + 38);
    const auto *ff0_39 = buffer.data(ff0 + 39);
    const auto *ff0_40 = buffer.data(ff0 + 40);
    const auto *ff0_41 = buffer.data(ff0 + 41);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_19 = buffer.data(ff1 + 19);
    const auto *ff1_20 = buffer.data(ff1 + 20);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_29 = buffer.data(ff1 + 29);
    const auto *ff1_31 = buffer.data(ff1 + 31);
    const auto *ff1_32 = buffer.data(ff1 + 32);
    const auto *ff1_33 = buffer.data(ff1 + 33);
    const auto *ff1_34 = buffer.data(ff1 + 34);
    const auto *ff1_35 = buffer.data(ff1 + 35);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, dg_0, ff0_0, ff0_1, ff1_0, \
                         ff1_1, fg_0, fg_1, fg_2, fg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_2[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];

        t_3[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, ff0_2, ff0_3, ff0_5, ff1_2, ff1_3, ff1_5, \
                         fg_4, fg_5, fg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_5[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];

        t_6[k] = f_5 * ff0_5[k]
                 - f_6 * ff1_5[k]
                 + pb_y[k] * fg_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, pa_x, pa_y, pa_z, pb_y, pb_z, dh_0, dh_1, \
                         ff0_6, ff1_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * ff0_6[k]
                 - f_4 * ff1_6[k]
                 + pb_y[k] * fg_7[k];

        t_8[k] = f_1 * ff0_6[k]
                 - f_2 * ff1_6[k]
                 + pb_z[k] * fg_8[k];

        t_9[k] = pa_y[k] * dh_0[k];

        t_10[k] = pa_z[k] * dh_0[k];

        t_11[k] = pa_x[k] * dh_1[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, dh_2, ff0_20, ff0_22, ff0_23, \
                         ff1_17, ff1_19, ff1_20, fg_13, fg_14, fg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_x[k] * dh_2[k];

        t_13[k] = f_1 * ff0_20[k]
                  - f_2 * ff1_17[k]
                  + pb_x[k] * fg_13[k];

        t_14[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_19[k]
                  + pb_x[k] * fg_14[k];

        t_15[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_20[k]
                  + pb_x[k] * fg_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_x, pb_y, pb_z, dg_9, ff0_24, ff0_26, \
                         ff1_21, ff1_23, fg_16, fg_17, fg_18, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_21[k]
                  + pb_x[k] * fg_16[k];

        t_17[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_23[k]
                  + pb_x[k] * fg_17[k];

        t_18[k] = f_0 * dg_9[k]
                  + f_1 * ff0_24[k]
                  - f_2 * ff1_21[k]
                  + pb_y[k] * fg_18[k];

        t_19[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_21[k]
                  + pb_z[k] * fg_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_y, pa_z, pb_z, dh_1, dh_2, ff0_25, ff0_26, \
                         ff1_22, ff1_23, fg_20, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_22[k]
                  + pb_z[k] * fg_20[k];

        t_21[k] = f_1 * ff0_26[k]
                  - f_2 * ff1_23[k]
                  + pb_z[k] * fg_21[k];

        t_22[k] = pa_z[k] * dh_1[k];

        t_23[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_x, ff0_35, ff0_37, ff0_38, ff1_29, ff1_31, \
                         ff1_32, fg_24, fg_25, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_1 * ff0_35[k]
                  - f_2 * ff1_29[k]
                  + pb_x[k] * fg_24[k];

        t_25[k] = f_5 * ff0_37[k]
                  - f_6 * ff1_31[k]
                  + pb_x[k] * fg_25[k];

        t_26[k] = f_5 * ff0_38[k]
                  - f_6 * ff1_32[k]
                  + pb_x[k] * fg_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_x, pb_y, ff0_39, ff0_40, ff0_41, ff1_33, \
                         ff1_34, ff1_35, fg_27, fg_28, fg_29, fg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_3 * ff0_39[k]
                  - f_4 * ff1_33[k]
                  + pb_x[k] * fg_27[k];

        t_28[k] = f_3 * ff0_41[k]
                  - f_4 * ff1_35[k]
                  + pb_x[k] * fg_28[k];

        t_29[k] = f_1 * ff0_39[k]
                  - f_2 * ff1_33[k]
                  + pb_y[k] * fg_29[k];

        t_30[k] = f_5 * ff0_40[k]
                  - f_6 * ff1_34[k]
                  + pb_y[k] * fg_30[k];
    }

#pragma omp simd aligned(t_31, t_32, pb_y, pb_z, dg_17, ff0_41, ff1_35, fg_31, \
                         fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_3 * ff0_41[k]
                  - f_4 * ff1_35[k]
                  + pb_y[k] * fg_31[k];

        t_32[k] = f_0 * dg_17[k]
                  + f_1 * ff0_41[k]
                  - f_2 * ff1_35[k]
                  + pb_z[k] * fg_32[k];
    }
}

auto
compute_prim_fh_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_8 = buffer.data(ff0 + 8);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_56 = buffer.data(fg + 56);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_1, dh_1, dh_2, ff0_5, ff1_15, \
                         fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_1[k]
                 + f_1 * ff0_5[k]
                 - f_2 * ff1_15[k]
                 + pb_y[k] * fg_32[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_2, ff0_8, ff1_26, fg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_2[k]
                 + f_1 * ff0_8[k]
                 - f_2 * ff1_26[k]
                 + pb_z[k] * fg_56[k];
    }
}

auto
compute_prim_fh_electron_repulsion_18(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.0 / p;
    const auto f_8 = 2.5 / p;
    const auto f_9 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_18 = buffer.data(ff1 + 18);
    const auto *ff1_19 = buffer.data(ff1 + 19);
    const auto *ff1_20 = buffer.data(ff1 + 20);
    const auto *ff1_26 = buffer.data(ff1 + 26);
    const auto *ff1_28 = buffer.data(ff1 + 28);
    const auto *ff1_29 = buffer.data(ff1 + 29);
    const auto *ff1_30 = buffer.data(ff1 + 30);
    const auto *ff1_31 = buffer.data(ff1 + 31);
    const auto *ff1_32 = buffer.data(ff1 + 32);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
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
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, \
                         ff1_2, ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_z[k] * fg_3[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];

        t_9[k] = pb_z[k] * fg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, ff1_5, \
                         ff1_6, fg_7, fg_8, fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * ff0_4[k]
                  - f_6 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_6[k]
                  + pb_y[k] * fg_8[k];

        t_12[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_6[k]
                  + pb_z[k] * fg_9[k];

        t_13[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pa_z, dg_1, dg_3, dg_5, dh_0, dh_1, \
                         dh_3, dh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_7 * dg_1[k]
                  + pa_y[k] * dh_1[k];

        t_15[k] = f_0 * dg_3[k]
                  + pa_y[k] * dh_3[k];

        t_16[k] = f_8 * dg_5[k]
                  + pa_y[k] * dh_5[k];

        t_17[k] = pa_z[k] * dh_0[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_z, pb_z, dg_0, dg_2, dg_4, dg_6, dh_2, \
                         dh_4, dh_6, fg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_9 * dg_0[k]
                  + pb_z[k] * fg_12[k];

        t_19[k] = f_7 * dg_2[k]
                  + pa_z[k] * dh_2[k];

        t_20[k] = f_0 * dg_4[k]
                  + pa_z[k] * dh_4[k];

        t_21[k] = f_7 * dg_6[k]
                  + pa_z[k] * dh_6[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pa_z, dg_7, dg_8, dg_10, dg_12, dh_7, \
                         dh_8, dh_9, dh_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_0 * dg_7[k]
                  + pa_z[k] * dh_7[k];

        t_23[k] = f_8 * dg_8[k]
                  + pa_z[k] * dh_8[k];

        t_24[k] = f_8 * dg_10[k]
                  + pa_x[k] * dh_9[k];

        t_25[k] = f_0 * dg_12[k]
                  + pa_x[k] * dh_10[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_x, pb_x, dg_14, dg_15, dg_21, dh_12, \
                         dh_14, dh_18, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * dg_14[k]
                  + pa_x[k] * dh_12[k];

        t_27[k] = f_9 * dg_15[k]
                  + pb_x[k] * fg_17[k];

        t_28[k] = pa_x[k] * dh_14[k];

        t_29[k] = f_8 * dg_21[k]
                  + pa_x[k] * dh_18[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, pb_z, dg_9, dg_24, dg_25, dg_29, \
                         dh_20, dh_22, fg_18, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * dg_9[k]
                  + pb_z[k] * fg_18[k];

        t_31[k] = f_0 * dg_24[k]
                  + pa_x[k] * dh_20[k];

        t_32[k] = f_7 * dg_25[k]
                  + pa_x[k] * dh_22[k];

        t_33[k] = f_9 * dg_29[k]
                  + pb_x[k] * fg_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, dh_26, ff0_12, ff0_13, \
                         ff1_14, ff1_16, fg_22, fg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_x[k] * dh_26[k];

        t_35[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_14[k]
                  + pb_x[k] * fg_22[k];

        t_36[k] = pb_z[k] * fg_22[k];

        t_37[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_16[k]
                  + pb_x[k] * fg_24[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_z, ff0_14, ff0_15, ff0_17, ff1_17, \
                         ff1_18, ff1_20, fg_24, fg_25, fg_26, fg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_17[k]
                  + pb_x[k] * fg_25[k];

        t_39[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_18[k]
                  + pb_x[k] * fg_26[k];

        t_40[k] = pb_z[k] * fg_24[k];

        t_41[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_20[k]
                  + pb_x[k] * fg_27[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_x, pb_y, pb_z, dg_15, ff0_15, \
                         ff0_16, ff1_18, ff1_19, fg_28, fg_29, fg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_x[k] * fg_28[k];

        t_43[k] = f_0 * dg_15[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_18[k]
                  + pb_y[k] * fg_28[k];

        t_44[k] = pb_z[k] * fg_28[k];

        t_45[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_18[k]
                  + pb_z[k] * fg_29[k];

        t_46[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_19[k]
                  + pb_z[k] * fg_30[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_z, pb_y, pb_z, dg_11, dg_13, dg_18, dh_11, \
                         dh_13, ff0_17, ff1_20, fg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * dg_18[k]
                  + pb_y[k] * fg_31[k];

        t_48[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_20[k]
                  + pb_z[k] * fg_31[k];

        t_49[k] = f_7 * dg_11[k]
                  + pa_z[k] * dh_11[k];

        t_50[k] = f_0 * dg_13[k]
                  + pa_z[k] * dh_13[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_z, dg_15, dg_16, dg_17, dh_14, \
                         dh_15, dh_16, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_z[k] * dh_14[k];

        t_52[k] = f_9 * dg_15[k]
                  + pb_z[k] * fg_32[k];

        t_53[k] = f_7 * dg_16[k]
                  + pa_z[k] * dh_15[k];

        t_54[k] = f_0 * dg_17[k]
                  + pa_z[k] * dh_16[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pa_z, pb_y, dg_18, dg_20, dg_22, dg_23, \
                         dh_17, dh_19, dh_21, fg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_7 * dg_20[k]
                  + pb_y[k] * fg_33[k];

        t_56[k] = f_8 * dg_18[k]
                  + pa_z[k] * dh_17[k];

        t_57[k] = f_7 * dg_22[k]
                  + pa_y[k] * dh_19[k];

        t_58[k] = f_0 * dg_23[k]
                  + pa_y[k] * dh_21[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pb_z, dg_19, dg_26, dg_27, dg_28, \
                         dh_23, dh_24, dh_25, fg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_8 * dg_26[k]
                  + pa_y[k] * dh_23[k];

        t_60[k] = f_7 * dg_19[k]
                  + pb_z[k] * fg_34[k];

        t_61[k] = f_0 * dg_27[k]
                  + pa_y[k] * dh_24[k];

        t_62[k] = f_7 * dg_28[k]
                  + pa_y[k] * dh_25[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, pb_z, dg_21, dg_29, \
                         dh_26, ff0_21, ff1_26, fg_37, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_9 * dg_29[k]
                  + pb_y[k] * fg_37[k];

        t_64[k] = pa_y[k] * dh_26[k];

        t_65[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_26[k]
                  + pb_x[k] * fg_38[k];

        t_66[k] = pb_y[k] * fg_38[k];

        t_67[k] = f_0 * dg_21[k]
                  + pb_z[k] * fg_38[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, ff0_22, ff0_23, ff0_24, ff1_28, \
                         ff1_29, ff1_30, fg_40, fg_41, fg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_28[k]
                  + pb_x[k] * fg_40[k];

        t_69[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_29[k]
                  + pb_x[k] * fg_41[k];

        t_70[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_30[k]
                  + pb_x[k] * fg_42[k];

        t_71[k] = pb_y[k] * fg_41[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pb_x, pb_y, pb_z, dg_26, ff0_24, ff0_26, \
                         ff1_30, ff1_32, fg_43, fg_44, fg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_32[k]
                  + pb_x[k] * fg_43[k];

        t_73[k] = pb_x[k] * fg_47[k];

        t_74[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_30[k]
                  + pb_y[k] * fg_44[k];

        t_75[k] = f_0 * dg_26[k]
                  + pb_z[k] * fg_44[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pb_y, pb_z, dg_29, ff0_25, ff0_26, ff1_31, \
                         ff1_32, fg_45, fg_46, fg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_31[k]
                  + pb_y[k] * fg_45[k];

        t_77[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_32[k]
                  + pb_y[k] * fg_46[k];

        t_78[k] = pb_y[k] * fg_47[k];

        t_79[k] = f_0 * dg_29[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_32[k]
                  + pb_z[k] * fg_47[k];
    }
}

auto
compute_prim_fh_electron_repulsion_19(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_18 = buffer.data(ff0 + 18);
    const auto *ff0_19 = buffer.data(ff0 + 19);
    const auto *ff0_20 = buffer.data(ff0 + 20);
    const auto *ff0_26 = buffer.data(ff0 + 26);
    const auto *ff0_28 = buffer.data(ff0 + 28);
    const auto *ff0_29 = buffer.data(ff0 + 29);
    const auto *ff0_30 = buffer.data(ff0 + 30);
    const auto *ff0_31 = buffer.data(ff0 + 31);
    const auto *ff0_32 = buffer.data(ff0 + 32);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, \
                         ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, dh_0, ff0_5, ff0_6, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ff0_5[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_10[k] = f_3 * ff0_6[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = pb_y[k] * fg_8[k];

        t_12[k] = f_1 * ff0_6[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_13[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_z, pb_x, dh_0, dh_1, dh_2, ff0_14, \
                         ff1_12, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * dh_0[k];

        t_15[k] = pa_x[k] * dh_1[k];

        t_16[k] = pa_x[k] * dh_2[k];

        t_17[k] = f_1 * ff0_14[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_13[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_x, ff0_16, ff0_17, ff0_18, ff1_13, ff1_14, \
                         ff1_15, fg_14, fg_15, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_14[k];

        t_19[k] = f_5 * ff0_17[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_15[k];

        t_20[k] = f_3 * ff0_18[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_16[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_y, dg_14, ff0_18, ff0_20, \
                         ff1_15, ff1_17, fg_17, fg_18, fg_20, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_3 * ff0_20[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_17[k];

        t_22[k] = pb_x[k] * fg_18[k];

        t_23[k] = pb_x[k] * fg_20[k];

        t_24[k] = pb_x[k] * fg_21[k];

        t_25[k] = f_0 * dg_14[k]
                  + f_1 * ff0_18[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_18[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, ff0_18, ff0_19, ff0_20, ff1_15, ff1_16, \
                         ff1_17, fg_18, fg_19, fg_20, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * fg_18[k];

        t_27[k] = f_3 * ff0_18[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_19[k];

        t_28[k] = f_5 * ff0_19[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_20[k];

        t_29[k] = f_1 * ff0_20[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pa_z, pb_x, dh_1, dh_2, ff0_26, ff0_28, \
                         ff1_21, ff1_22, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_z[k] * dh_1[k];

        t_31[k] = pa_y[k] * dh_2[k];

        t_32[k] = f_1 * ff0_26[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_24[k];

        t_33[k] = f_5 * ff0_28[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_25[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_x, ff0_29, ff0_30, ff0_32, ff1_23, ff1_24, \
                         ff1_26, fg_26, fg_27, fg_28, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_5 * ff0_29[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_26[k];

        t_35[k] = f_3 * ff0_30[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_27[k];

        t_36[k] = f_3 * ff0_32[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_28[k];

        t_37[k] = pb_x[k] * fg_29[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_y, ff0_30, ff0_31, ff1_24, ff1_25, \
                         fg_29, fg_30, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_x[k] * fg_30[k];

        t_39[k] = pb_x[k] * fg_32[k];

        t_40[k] = f_1 * ff0_30[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_29[k];

        t_41[k] = f_5 * ff0_31[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_y, pb_z, dg_26, ff0_32, ff1_26, fg_31, \
                         fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * ff0_32[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_31[k];

        t_43[k] = pb_y[k] * fg_32[k];

        t_44[k] = f_0 * dg_26[k]
                  + f_1 * ff0_32[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_32[k];
    }
}

auto
compute_prim_fh_electron_repulsion_20(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_8 = buffer.data(ff0 + 8);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_38 = buffer.data(fg + 38);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_1, dh_1, dh_2, ff0_5, ff1_15, \
                         fg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_1[k]
                 + f_1 * ff0_5[k]
                 - f_2 * ff1_15[k]
                 + pb_y[k] * fg_22[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_2, ff0_8, ff1_26, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_2[k]
                 + f_1 * ff0_8[k]
                 - f_2 * ff1_26[k]
                 + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_21(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 2.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_38 = buffer.data(dh + 38);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, \
                         ff1_2, ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_z[k] * fg_3[k];

        t_7[k] = pb_y[k] * fg_4[k];

        t_8[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_9[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pb_y, pb_z, ff0_4, ff0_5, ff1_4, ff1_5, \
                         fg_5, fg_7, fg_8, fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_z[k] * fg_5[k];

        t_11[k] = f_5 * ff0_4[k]
                  - f_6 * ff1_4[k]
                  + pb_y[k] * fg_7[k];

        t_12[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_8[k];

        t_13[k] = pb_y[k] * fg_9[k];

        t_14[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, t_20, pa_y, pa_z, dg_5, dg_6, dh_0, \
                         dh_7, dh_9, dh_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * dh_0[k];

        t_16[k] = f_7 * dg_5[k]
                  + pa_y[k] * dh_7[k];

        t_17[k] = pa_y[k] * dh_11[k];

        t_18[k] = pa_z[k] * dh_0[k];

        t_19[k] = pa_z[k] * dh_7[k];

        t_20[k] = f_8 * dg_6[k]
                  + pa_z[k] * dh_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_x, pa_z, pb_y, dg_7, dg_9, dg_10, dh_10, \
                         dh_11, dh_13, fg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * dg_7[k]
                  + pa_z[k] * dh_10[k];

        t_22[k] = pb_y[k] * fg_12[k];

        t_23[k] = f_7 * dg_9[k]
                  + pa_z[k] * dh_11[k];

        t_24[k] = f_7 * dg_10[k]
                  + pa_x[k] * dh_13[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_x, pb_x, dg_12, dg_14, dg_16, dg_20, \
                         dh_15, dh_17, dh_20, dh_27, fg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * dg_12[k]
                  + pa_x[k] * dh_15[k];

        t_26[k] = f_8 * dg_14[k]
                  + pa_x[k] * dh_17[k];

        t_27[k] = f_9 * dg_16[k]
                  + pb_x[k] * fg_15[k];

        t_28[k] = pa_x[k] * dh_20[k];

        t_29[k] = f_7 * dg_20[k]
                  + pa_x[k] * dh_27[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, dg_23, dg_25, dg_29, dh_30, \
                         dh_32, dh_38, fg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * dg_23[k]
                  + pa_x[k] * dh_30[k];

        t_31[k] = f_8 * dg_25[k]
                  + pa_x[k] * dh_32[k];

        t_32[k] = f_9 * dg_29[k]
                  + pb_x[k] * fg_18[k];

        t_33[k] = pa_x[k] * dh_38[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_x, pb_z, ff0_12, ff0_13, ff0_14, ff1_12, \
                         ff1_13, ff1_14, fg_19, fg_21, fg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_19[k];

        t_35[k] = pb_z[k] * fg_19[k];

        t_36[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_21[k];

        t_37[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_22[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pb_x, pb_z, ff0_15, ff0_17, ff1_15, \
                         ff1_17, fg_21, fg_23, fg_24, fg_25, fg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_23[k];

        t_39[k] = pb_z[k] * fg_21[k];

        t_40[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_24[k];

        t_41[k] = pb_x[k] * fg_25[k];

        t_42[k] = pb_x[k] * fg_27[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pb_x, pb_y, pb_z, dg_16, ff0_15, ff1_15, \
                         fg_25, fg_26, fg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_x[k] * fg_28[k];

        t_44[k] = f_0 * dg_16[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_25[k];

        t_45[k] = pb_z[k] * fg_25[k];

        t_46[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_26[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_z, pb_z, dg_11, dh_13, dh_16, ff0_16, \
                         ff0_17, ff1_16, ff1_17, fg_27, fg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_27[k];

        t_48[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_28[k];

        t_49[k] = pa_z[k] * dh_13[k];

        t_50[k] = f_8 * dg_11[k]
                  + pa_z[k] * dh_16[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pa_y, pa_z, pb_x, dg_13, dg_19, dg_21, \
                         dh_18, dh_20, dh_24, dh_29, fg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * dg_13[k]
                  + pa_z[k] * dh_18[k];

        t_52[k] = pb_x[k] * fg_30[k];

        t_53[k] = pa_z[k] * dh_20[k];

        t_54[k] = f_7 * dg_19[k]
                  + pa_z[k] * dh_24[k];

        t_55[k] = f_8 * dg_21[k]
                  + pa_y[k] * dh_29[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_y, pb_x, dg_22, dg_26, dg_27, dg_28, \
                         dh_31, dh_34, dh_35, dh_36, fg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * dg_22[k]
                  + pa_y[k] * dh_31[k];

        t_57[k] = pb_x[k] * fg_31[k];

        t_58[k] = f_7 * dg_26[k]
                  + pa_y[k] * dh_34[k];

        t_59[k] = f_0 * dg_27[k]
                  + pa_y[k] * dh_35[k];

        t_60[k] = f_8 * dg_28[k]
                  + pa_y[k] * dh_36[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_y, pb_x, pb_y, dg_29, dh_38, ff0_21, \
                         ff1_21, fg_34, fg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_9 * dg_29[k]
                  + pb_y[k] * fg_34[k];

        t_62[k] = pa_y[k] * dh_38[k];

        t_63[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_35[k];

        t_64[k] = pb_y[k] * fg_35[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_x, pb_y, ff0_22, ff0_23, ff0_24, ff1_22, \
                         ff1_23, ff1_24, fg_37, fg_38, fg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_37[k];

        t_66[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_38[k];

        t_67[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_39[k];

        t_68[k] = pb_y[k] * fg_38[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pb_x, pb_y, ff0_24, ff0_26, ff1_24, \
                         ff1_26, fg_40, fg_41, fg_42, fg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_40[k];

        t_70[k] = pb_x[k] * fg_41[k];

        t_71[k] = pb_x[k] * fg_42[k];

        t_72[k] = pb_x[k] * fg_44[k];

        t_73[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_41[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_y, pb_z, dg_29, ff0_25, ff0_26, ff1_25, \
                         ff1_26, fg_42, fg_43, fg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_42[k];

        t_75[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_43[k];

        t_76[k] = pb_y[k] * fg_44[k];

        t_77[k] = f_0 * dg_29[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_44[k];
    }
}

auto
compute_prim_fh_electron_repulsion_22(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 2.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_35 = buffer.data(dh + 35);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, \
                         ff1_2, ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];

        t_9[k] = pb_z[k] * fg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * ff0_4[k]
                  - f_6 * ff1_4[k]
                  + pb_y[k] * fg_6[k];

        t_11[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_12[k] = pb_y[k] * fg_8[k];

        t_13[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_14[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pa_y, pa_z, dg_5, dg_8, dg_9, dh_0, \
                         dh_7, dh_11, dh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_7 * dg_5[k]
                  + pa_y[k] * dh_7[k];

        t_16[k] = pa_z[k] * dh_0[k];

        t_17[k] = f_7 * dg_8[k]
                  + pa_z[k] * dh_11[k];

        t_18[k] = f_7 * dg_9[k]
                  + pa_x[k] * dh_12[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_x, pb_x, dg_10, dg_12, dg_14, dg_18, \
                         dh_13, dh_15, dh_19, dh_24, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * dg_10[k]
                  + pa_x[k] * dh_13[k];

        t_20[k] = f_8 * dg_12[k]
                  + pa_x[k] * dh_15[k];

        t_21[k] = f_9 * dg_14[k]
                  + pb_x[k] * fg_13[k];

        t_22[k] = pa_x[k] * dh_19[k];

        t_23[k] = f_7 * dg_18[k]
                  + pa_x[k] * dh_24[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_x, pb_x, dg_20, dg_22, dg_26, dh_26, \
                         dh_28, dh_35, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * dg_20[k]
                  + pa_x[k] * dh_26[k];

        t_25[k] = f_8 * dg_22[k]
                  + pa_x[k] * dh_28[k];

        t_26[k] = f_9 * dg_26[k]
                  + pb_x[k] * fg_16[k];

        t_27[k] = pa_x[k] * dh_35[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_x, pb_z, ff0_12, ff0_13, ff0_14, ff1_12, \
                         ff1_13, ff1_14, fg_17, fg_18, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_17[k];

        t_29[k] = pb_z[k] * fg_17[k];

        t_30[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_18[k];

        t_31[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_19[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, pb_x, ff0_15, ff0_17, ff1_15, ff1_17, \
                         fg_20, fg_21, fg_22, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_20[k];

        t_33[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_21[k];

        t_34[k] = pb_x[k] * fg_22[k];

        t_35[k] = pb_x[k] * fg_24[k];

        t_36[k] = pb_x[k] * fg_25[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_y, pb_z, dg_14, ff0_15, ff0_16, ff1_15, \
                         ff1_16, fg_22, fg_23, fg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_0 * dg_14[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_22[k];

        t_38[k] = pb_z[k] * fg_22[k];

        t_39[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_23[k];

        t_40[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_24[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_y, pa_z, pb_z, dg_17, dg_23, dh_19, dh_23, \
                         dh_31, ff0_17, ff1_17, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_25[k];

        t_42[k] = pa_z[k] * dh_19[k];

        t_43[k] = f_7 * dg_17[k]
                  + pa_z[k] * dh_23[k];

        t_44[k] = f_7 * dg_23[k]
                  + pa_y[k] * dh_31[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_y, pb_y, dg_24, dg_25, dg_26, dh_32, \
                         dh_33, dh_35, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_0 * dg_24[k]
                  + pa_y[k] * dh_32[k];

        t_46[k] = f_8 * dg_25[k]
                  + pa_y[k] * dh_33[k];

        t_47[k] = f_9 * dg_26[k]
                  + pb_y[k] * fg_29[k];

        t_48[k] = pa_y[k] * dh_35[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pb_x, pb_y, ff0_21, ff0_22, ff0_23, ff1_21, \
                         ff1_22, ff1_23, fg_30, fg_31, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_30[k];

        t_50[k] = pb_y[k] * fg_30[k];

        t_51[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_31[k];

        t_52[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_32[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pb_x, ff0_24, ff0_26, ff1_24, ff1_26, \
                         fg_33, fg_34, fg_35, fg_36, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_33[k];

        t_54[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_34[k];

        t_55[k] = pb_x[k] * fg_35[k];

        t_56[k] = pb_x[k] * fg_36[k];

        t_57[k] = pb_x[k] * fg_38[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_y, ff0_24, ff0_25, ff0_26, ff1_24, ff1_25, \
                         ff1_26, fg_35, fg_36, fg_37, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_35[k];

        t_59[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_36[k];

        t_60[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_37[k];

        t_61[k] = pb_y[k] * fg_38[k];
    }

#pragma omp simd aligned(t_62, pb_z, dg_26, ff0_26, ff1_26, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_0 * dg_26[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_23(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, \
                         ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_10[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = pb_y[k] * fg_8[k];

        t_12[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_13[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_z, pb_x, dh_0, dh_1, dh_2, ff0_12, \
                         ff1_12, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * dh_0[k];

        t_15[k] = pa_x[k] * dh_1[k];

        t_16[k] = pa_x[k] * dh_2[k];

        t_17[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_13[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_x, ff0_13, ff0_14, ff0_15, ff1_13, ff1_14, \
                         ff1_15, fg_14, fg_15, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_14[k];

        t_19[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_15[k];

        t_20[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_16[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_y, dg_14, ff0_15, ff0_17, \
                         ff1_15, ff1_17, fg_17, fg_18, fg_20, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_17[k];

        t_22[k] = pb_x[k] * fg_18[k];

        t_23[k] = pb_x[k] * fg_20[k];

        t_24[k] = pb_x[k] * fg_21[k];

        t_25[k] = f_0 * dg_14[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_18[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, ff0_15, ff0_16, ff0_17, ff1_15, ff1_16, \
                         ff1_17, fg_18, fg_19, fg_20, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * fg_18[k];

        t_27[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_19[k];

        t_28[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_20[k];

        t_29[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pa_z, pb_x, dh_1, dh_2, ff0_21, ff0_22, \
                         ff1_21, ff1_22, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_z[k] * dh_1[k];

        t_31[k] = pa_y[k] * dh_2[k];

        t_32[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_24[k];

        t_33[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_25[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_x, ff0_23, ff0_24, ff0_26, ff1_23, ff1_24, \
                         ff1_26, fg_26, fg_27, fg_28, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_26[k];

        t_35[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_27[k];

        t_36[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_28[k];

        t_37[k] = pb_x[k] * fg_29[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_y, ff0_24, ff0_25, ff1_24, ff1_25, \
                         fg_29, fg_30, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_x[k] * fg_30[k];

        t_39[k] = pb_x[k] * fg_32[k];

        t_40[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_29[k];

        t_41[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_y, pb_z, dg_26, ff0_26, ff1_26, fg_31, \
                         fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_31[k];

        t_43[k] = pb_y[k] * fg_32[k];

        t_44[k] = f_0 * dg_26[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_32[k];
    }
}

auto
compute_prim_fh_electron_repulsion_24(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_38 = buffer.data(dh + 38);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, \
                         ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_10[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = pb_y[k] * fg_8[k];

        t_12[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_13[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_x, pa_z, pb_x, dg_10, dg_12, dg_14, \
                         dh_0, dh_15, dh_17, dh_21, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * dh_0[k];

        t_15[k] = f_0 * dg_10[k]
                  + pa_x[k] * dh_15[k];

        t_16[k] = f_7 * dg_12[k]
                  + pa_x[k] * dh_17[k];

        t_17[k] = f_8 * dg_14[k]
                  + pb_x[k] * fg_13[k];

        t_18[k] = pa_x[k] * dh_21[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, dg_20, dg_22, dg_26, dh_29, \
                         dh_31, dh_38, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * dg_20[k]
                  + pa_x[k] * dh_29[k];

        t_20[k] = f_7 * dg_22[k]
                  + pa_x[k] * dh_31[k];

        t_21[k] = f_8 * dg_26[k]
                  + pb_x[k] * fg_16[k];

        t_22[k] = pa_x[k] * dh_38[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_x, ff0_12, ff0_13, ff0_14, ff1_12, ff1_13, \
                         ff1_14, fg_17, fg_18, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_17[k];

        t_24[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_18[k];

        t_25[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_19[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pb_x, ff0_15, ff0_17, ff1_15, ff1_17, \
                         fg_20, fg_21, fg_22, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_20[k];

        t_27[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_21[k];

        t_28[k] = pb_x[k] * fg_22[k];

        t_29[k] = pb_x[k] * fg_24[k];

        t_30[k] = pb_x[k] * fg_25[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pb_y, pb_z, dg_14, ff0_15, ff0_16, ff1_15, \
                         ff1_16, fg_22, fg_23, fg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * dg_14[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_22[k];

        t_32[k] = pb_z[k] * fg_22[k];

        t_33[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_23[k];

        t_34[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_24[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_y, pa_z, pb_z, dg_24, dg_25, dh_21, dh_35, \
                         dh_36, ff0_17, ff1_17, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_25[k];

        t_36[k] = pa_z[k] * dh_21[k];

        t_37[k] = f_0 * dg_24[k]
                  + pa_y[k] * dh_35[k];

        t_38[k] = f_7 * dg_25[k]
                  + pa_y[k] * dh_36[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_y, pb_x, pb_y, dg_26, dh_38, ff0_21, \
                         ff0_22, ff1_21, ff1_22, fg_29, fg_30, fg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_8 * dg_26[k]
                  + pb_y[k] * fg_29[k];

        t_40[k] = pa_y[k] * dh_38[k];

        t_41[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_30[k];

        t_42[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_31[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pb_x, ff0_23, ff0_24, ff0_26, ff1_23, ff1_24, \
                         ff1_26, fg_32, fg_33, fg_34, fg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_32[k];

        t_44[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_33[k];

        t_45[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_34[k];

        t_46[k] = pb_x[k] * fg_35[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_x, pb_y, ff0_24, ff0_25, ff1_24, ff1_25, \
                         fg_35, fg_36, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_x[k] * fg_36[k];

        t_48[k] = pb_x[k] * fg_38[k];

        t_49[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_35[k];

        t_50[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_36[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_y, pb_z, dg_26, ff0_26, ff1_26, fg_37, \
                         fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_37[k];

        t_52[k] = pb_y[k] * fg_38[k];

        t_53[k] = f_0 * dg_26[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_25(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 2.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 0.5 / p;

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_35 = buffer.data(dh + 35);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, \
                         ff1_2, ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];

        t_9[k] = pb_z[k] * fg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * ff0_4[k]
                  - f_6 * ff1_4[k]
                  + pb_y[k] * fg_6[k];

        t_11[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_12[k] = pb_y[k] * fg_8[k];

        t_13[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_14[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pa_z, dg_8, dg_10, dg_12, dh_0, dh_11, \
                         dh_13, dh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_z[k] * dh_0[k];

        t_16[k] = f_7 * dg_8[k]
                  + pa_z[k] * dh_11[k];

        t_17[k] = f_0 * dg_10[k]
                  + pa_x[k] * dh_13[k];

        t_18[k] = f_8 * dg_12[k]
                  + pa_x[k] * dh_15[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_x, pb_x, dg_14, dg_20, dg_22, dg_26, \
                         dh_19, dh_26, dh_28, fg_13, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_9 * dg_14[k]
                  + pb_x[k] * fg_13[k];

        t_20[k] = pa_x[k] * dh_19[k];

        t_21[k] = f_0 * dg_20[k]
                  + pa_x[k] * dh_26[k];

        t_22[k] = f_8 * dg_22[k]
                  + pa_x[k] * dh_28[k];

        t_23[k] = f_9 * dg_26[k]
                  + pb_x[k] * fg_16[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_x, pb_x, pb_z, dh_35, ff0_12, ff0_13, \
                         ff1_12, ff1_13, fg_17, fg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_x[k] * dh_35[k];

        t_25[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_17[k];

        t_26[k] = pb_z[k] * fg_17[k];

        t_27[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_18[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_x, ff0_14, ff0_15, ff0_17, ff1_14, ff1_15, \
                         ff1_17, fg_19, fg_20, fg_21, fg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_19[k];

        t_29[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_20[k];

        t_30[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_21[k];

        t_31[k] = pb_x[k] * fg_22[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, pb_x, pb_y, pb_z, dg_14, ff0_15, \
                         ff1_15, fg_22, fg_23, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_x[k] * fg_24[k];

        t_33[k] = pb_x[k] * fg_25[k];

        t_34[k] = f_0 * dg_14[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_22[k];

        t_35[k] = pb_z[k] * fg_22[k];

        t_36[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_23[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_z, dg_17, dh_19, dh_23, ff0_16, \
                         ff0_17, ff1_16, ff1_17, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_24[k];

        t_38[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_25[k];

        t_39[k] = pa_z[k] * dh_19[k];

        t_40[k] = f_7 * dg_17[k]
                  + pa_z[k] * dh_23[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_y, pb_y, dg_23, dg_24, dg_25, dg_26, \
                         dh_31, dh_32, dh_33, dh_35, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * dg_23[k]
                  + pa_y[k] * dh_31[k];

        t_42[k] = f_0 * dg_24[k]
                  + pa_y[k] * dh_32[k];

        t_43[k] = f_8 * dg_25[k]
                  + pa_y[k] * dh_33[k];

        t_44[k] = f_9 * dg_26[k]
                  + pb_y[k] * fg_29[k];

        t_45[k] = pa_y[k] * dh_35[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_x, pb_y, ff0_21, ff0_22, ff0_23, ff1_21, \
                         ff1_22, ff1_23, fg_30, fg_31, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_30[k];

        t_47[k] = pb_y[k] * fg_30[k];

        t_48[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_31[k];

        t_49[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_32[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pb_x, ff0_24, ff0_26, ff1_24, ff1_26, \
                         fg_33, fg_34, fg_35, fg_36, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_33[k];

        t_51[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_34[k];

        t_52[k] = pb_x[k] * fg_35[k];

        t_53[k] = pb_x[k] * fg_36[k];

        t_54[k] = pb_x[k] * fg_38[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_y, ff0_24, ff0_25, ff0_26, ff1_24, ff1_25, \
                         ff1_26, fg_35, fg_36, fg_37, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_35[k];

        t_56[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_36[k];

        t_57[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_37[k];

        t_58[k] = pb_y[k] * fg_38[k];
    }

#pragma omp simd aligned(t_59, pb_z, dg_26, ff0_26, ff1_26, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_0 * dg_26[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_26(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_35 = buffer.data(dh + 35);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, \
                         ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = pb_y[k] * fg_4[k];

        t_7[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_10[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = pb_y[k] * fg_8[k];

        t_12[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_13[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_x, pa_z, pb_x, dg_10, dg_12, dg_14, \
                         dh_0, dh_13, dh_15, dh_19, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * dh_0[k];

        t_15[k] = f_0 * dg_10[k]
                  + pa_x[k] * dh_13[k];

        t_16[k] = f_7 * dg_12[k]
                  + pa_x[k] * dh_15[k];

        t_17[k] = f_8 * dg_14[k]
                  + pb_x[k] * fg_13[k];

        t_18[k] = pa_x[k] * dh_19[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, dg_20, dg_22, dg_26, dh_26, \
                         dh_28, dh_35, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * dg_20[k]
                  + pa_x[k] * dh_26[k];

        t_20[k] = f_7 * dg_22[k]
                  + pa_x[k] * dh_28[k];

        t_21[k] = f_8 * dg_26[k]
                  + pb_x[k] * fg_16[k];

        t_22[k] = pa_x[k] * dh_35[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_x, ff0_12, ff0_13, ff0_14, ff1_12, ff1_13, \
                         ff1_14, fg_17, fg_18, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_17[k];

        t_24[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_18[k];

        t_25[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_19[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pb_x, ff0_15, ff0_17, ff1_15, ff1_17, \
                         fg_20, fg_21, fg_22, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_20[k];

        t_27[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_21[k];

        t_28[k] = pb_x[k] * fg_22[k];

        t_29[k] = pb_x[k] * fg_24[k];

        t_30[k] = pb_x[k] * fg_25[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pb_y, pb_z, dg_14, ff0_15, ff0_16, ff1_15, \
                         ff1_16, fg_22, fg_23, fg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * dg_14[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_22[k];

        t_32[k] = pb_z[k] * fg_22[k];

        t_33[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_23[k];

        t_34[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_24[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_y, pa_z, pb_z, dg_24, dg_25, dh_19, dh_32, \
                         dh_33, ff0_17, ff1_17, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_25[k];

        t_36[k] = pa_z[k] * dh_19[k];

        t_37[k] = f_0 * dg_24[k]
                  + pa_y[k] * dh_32[k];

        t_38[k] = f_7 * dg_25[k]
                  + pa_y[k] * dh_33[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_y, pb_x, pb_y, dg_26, dh_35, ff0_21, \
                         ff0_22, ff1_21, ff1_22, fg_29, fg_30, fg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_8 * dg_26[k]
                  + pb_y[k] * fg_29[k];

        t_40[k] = pa_y[k] * dh_35[k];

        t_41[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_30[k];

        t_42[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_31[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pb_x, ff0_23, ff0_24, ff0_26, ff1_23, ff1_24, \
                         ff1_26, fg_32, fg_33, fg_34, fg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_32[k];

        t_44[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_33[k];

        t_45[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_34[k];

        t_46[k] = pb_x[k] * fg_35[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_x, pb_y, ff0_24, ff0_25, ff1_24, ff1_25, \
                         fg_35, fg_36, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_x[k] * fg_36[k];

        t_48[k] = pb_x[k] * fg_38[k];

        t_49[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_35[k];

        t_50[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_36[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_y, pb_z, dg_26, ff0_26, ff1_26, fg_37, \
                         fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_37[k];

        t_52[k] = pb_y[k] * fg_38[k];

        t_53[k] = f_0 * dg_26[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_27(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_5 = buffer.data(dg + 5);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_8 = buffer.data(ff0 + 8);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_11 = buffer.data(ff1 + 11);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_14 = buffer.data(fg + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_2, dh_1, dh_2, ff0_5, ff1_7, \
                         fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_2[k]
                 + f_1 * ff0_5[k]
                 - f_2 * ff1_7[k]
                 + pb_y[k] * fg_9[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_5, ff0_8, ff1_11, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_5[k]
                 + f_1 * ff0_8[k]
                 - f_2 * ff1_11[k]
                 + pb_z[k] * fg_14[k];
    }
}

auto
compute_prim_fh_electron_repulsion_28(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_5 = buffer.data(dg + 5);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_11 = buffer.data(ff0 + 11);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_19 = buffer.data(ff1 + 19);
    const auto *ff1_32 = buffer.data(ff1 + 32);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_29 = buffer.data(fg + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_2, dh_1, dh_2, ff0_7, ff1_19, \
                         fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_2[k]
                 + f_1 * ff0_7[k]
                 - f_2 * ff1_19[k]
                 + pb_y[k] * fg_17[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_5, ff0_11, ff1_32, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_5[k]
                 + f_1 * ff0_11[k]
                 - f_2 * ff1_32[k]
                 + pb_z[k] * fg_29[k];
    }
}

auto
compute_prim_fh_electron_repulsion_29(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_5 = buffer.data(dg + 5);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_19 = buffer.data(ff0 + 19);
    const auto *ff0_32 = buffer.data(ff0 + 32);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_2, dh_1, dh_2, ff0_19, ff1_15, \
                         fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_2[k]
                 + f_1 * ff0_19[k]
                 - f_2 * ff1_15[k]
                 + pb_y[k] * fg_5[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_5, ff0_32, ff1_26, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_5[k]
                 + f_1 * ff0_32[k]
                 - f_2 * ff1_26[k]
                 + pb_z[k] * fg_8[k];
    }
}

auto
compute_prim_fh_electron_repulsion_30(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_11 = buffer.data(ff0 + 11);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_11 = buffer.data(ff1 + 11);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, dg_0, dg_1, dh_0, dh_1, \
                         ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = f_0 * dg_1[k]
                 + pa_x[k] * dh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, dg_2, dg_3, dg_4, dg_5, dh_2, \
                         dh_3, dh_4, dh_5, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dg_2[k]
                 + pa_x[k] * dh_2[k];

        t_5[k] = f_4 * dg_3[k]
                 + pb_x[k] * fg_5[k];

        t_6[k] = pa_x[k] * dh_3[k];

        t_7[k] = f_0 * dg_4[k]
                 + pa_x[k] * dh_4[k];

        t_8[k] = f_3 * dg_5[k]
                 + pa_x[k] * dh_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pa_z, pb_x, pb_y, dg_3, dg_8, dh_3, \
                         dh_8, ff0_7, ff1_7, fg_8, fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dg_8[k]
                 + pb_x[k] * fg_8[k];

        t_10[k] = pa_x[k] * dh_8[k];

        t_11[k] = f_0 * dg_3[k]
                  + f_1 * ff0_7[k]
                  - f_2 * ff1_7[k]
                  + pb_y[k] * fg_9[k];

        t_12[k] = pa_z[k] * dh_3[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pb_y, dg_6, dg_7, dg_8, dh_6, dh_7, \
                         dh_8, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dg_6[k]
                  + pa_y[k] * dh_6[k];

        t_14[k] = f_3 * dg_7[k]
                  + pa_y[k] * dh_7[k];

        t_15[k] = f_4 * dg_8[k]
                  + pb_y[k] * fg_13[k];

        t_16[k] = pa_y[k] * dh_8[k];
    }

#pragma omp simd aligned(t_17, pb_z, dg_8, ff0_11, ff1_11, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * dg_8[k]
                  + f_1 * ff0_11[k]
                  - f_2 * ff1_11[k]
                  + pb_z[k] * fg_14[k];
    }
}

auto
compute_prim_fh_electron_repulsion_31(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_11 = buffer.data(ff0 + 11);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_29 = buffer.data(ff1 + 29);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_44 = buffer.data(fg + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, dg_0, dg_1, dh_0, dh_1, \
                         ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = f_0 * dg_1[k]
                 + pa_x[k] * dh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, dg_2, dg_3, dg_4, dg_5, dh_2, \
                         dh_3, dh_4, dh_5, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dg_2[k]
                 + pa_x[k] * dh_2[k];

        t_5[k] = f_4 * dg_3[k]
                 + pb_x[k] * fg_16[k];

        t_6[k] = pa_x[k] * dh_3[k];

        t_7[k] = f_0 * dg_4[k]
                 + pa_x[k] * dh_4[k];

        t_8[k] = f_3 * dg_5[k]
                 + pa_x[k] * dh_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pa_z, pb_x, pb_y, dg_3, dg_8, dh_3, \
                         dh_8, ff0_7, ff1_17, fg_20, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dg_8[k]
                 + pb_x[k] * fg_20[k];

        t_10[k] = pa_x[k] * dh_8[k];

        t_11[k] = f_0 * dg_3[k]
                  + f_1 * ff0_7[k]
                  - f_2 * ff1_17[k]
                  + pb_y[k] * fg_26[k];

        t_12[k] = pa_z[k] * dh_3[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pb_y, dg_6, dg_7, dg_8, dh_6, dh_7, \
                         dh_8, fg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dg_6[k]
                  + pa_y[k] * dh_6[k];

        t_14[k] = f_3 * dg_7[k]
                  + pa_y[k] * dh_7[k];

        t_15[k] = f_4 * dg_8[k]
                  + pb_y[k] * fg_35[k];

        t_16[k] = pa_y[k] * dh_8[k];
    }

#pragma omp simd aligned(t_17, pb_z, dg_8, ff0_11, ff1_29, fg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * dg_8[k]
                  + f_1 * ff0_11[k]
                  - f_2 * ff1_29[k]
                  + pb_z[k] * fg_44[k];
    }
}

auto
compute_prim_fh_electron_repulsion_32(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_18 = buffer.data(ff0 + 18);
    const auto *ff0_19 = buffer.data(ff0 + 19);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);
    const auto *ff0_27 = buffer.data(ff0 + 27);
    const auto *ff0_28 = buffer.data(ff0 + 28);
    const auto *ff0_29 = buffer.data(ff0 + 29);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, dg_0, ff0_0, ff0_1, ff1_0, \
                         ff1_1, fg_0, fg_1, fg_2, fg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_2[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];

        t_3[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, ff0_2, ff0_3, ff0_5, ff1_2, ff1_3, ff1_4, \
                         fg_4, fg_5, fg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_5[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];

        t_6[k] = f_5 * ff0_5[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pa_y, pa_z, pb_y, pb_z, dh_0, ff0_6, ff1_5, \
                         fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * ff0_6[k]
                 - f_4 * ff1_5[k]
                 + pb_y[k] * fg_7[k];

        t_8[k] = f_1 * ff0_6[k]
                 - f_2 * ff1_5[k]
                 + pb_z[k] * fg_8[k];

        t_9[k] = pa_y[k] * dh_0[k];

        t_10[k] = pa_z[k] * dh_0[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_x, pb_x, dg_4, dg_5, dg_6, dg_9, \
                         dh_1, dh_2, dh_3, dh_4, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * dg_4[k]
                  + pa_x[k] * dh_1[k];

        t_12[k] = f_7 * dg_5[k]
                  + pa_x[k] * dh_2[k];

        t_13[k] = f_8 * dg_6[k]
                  + pb_x[k] * fg_13[k];

        t_14[k] = pa_x[k] * dh_3[k];

        t_15[k] = f_0 * dg_9[k]
                  + pa_x[k] * dh_4[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_x, dg_10, dg_14, dh_5, dh_8, ff0_13, \
                         ff1_12, fg_16, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * dg_10[k]
                  + pa_x[k] * dh_5[k];

        t_17[k] = f_8 * dg_14[k]
                  + pb_x[k] * fg_16[k];

        t_18[k] = pa_x[k] * dh_8[k];

        t_19[k] = f_1 * ff0_13[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_17[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_x, ff0_15, ff0_16, ff0_17, ff1_13, ff1_14, \
                         ff1_15, fg_18, fg_19, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * ff0_15[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_18[k];

        t_21[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_19[k];

        t_22[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_20[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_x, pb_y, pb_z, dg_6, ff0_17, ff0_19, ff1_15, \
                         ff1_17, fg_21, fg_22, fg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * ff0_19[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_21[k];

        t_24[k] = f_0 * dg_6[k]
                  + f_1 * ff0_17[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_22[k];

        t_25[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_23[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, pb_z, dg_12, dh_3, dh_6, ff0_18, \
                         ff0_19, ff1_16, ff1_17, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * ff0_18[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_24[k];

        t_27[k] = f_1 * ff0_19[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_25[k];

        t_28[k] = pa_z[k] * dh_3[k];

        t_29[k] = f_0 * dg_12[k]
                  + pa_y[k] * dh_6[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_x, pb_y, dg_13, dg_14, dh_7, dh_8, \
                         ff0_23, ff1_21, fg_29, fg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * dg_13[k]
                  + pa_y[k] * dh_7[k];

        t_31[k] = f_8 * dg_14[k]
                  + pb_y[k] * fg_29[k];

        t_32[k] = pa_y[k] * dh_8[k];

        t_33[k] = f_1 * ff0_23[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_30[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pb_x, ff0_25, ff0_26, ff0_27, ff1_22, ff1_23, \
                         ff1_24, fg_31, fg_32, fg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_31[k];

        t_35[k] = f_5 * ff0_26[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_32[k];

        t_36[k] = f_3 * ff0_27[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_33[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_x, pb_y, ff0_27, ff0_28, ff0_29, ff1_24, \
                         ff1_25, ff1_26, fg_34, fg_35, fg_36, fg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * ff0_29[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_34[k];

        t_38[k] = f_1 * ff0_27[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_35[k];

        t_39[k] = f_5 * ff0_28[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_36[k];

        t_40[k] = f_3 * ff0_29[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_37[k];
    }

#pragma omp simd aligned(t_41, pb_z, dg_14, ff0_29, ff1_26, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_0 * dg_14[k]
                  + f_1 * ff0_29[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_33(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_3, dh_1, dh_2, ff0_15, ff1_15, \
                         fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_3[k]
                 + f_1 * ff0_15[k]
                 - f_2 * ff1_15[k]
                 + pb_y[k] * fg_5[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_8, ff0_26, ff1_26, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_8[k]
                 + f_1 * ff0_26[k]
                 - f_2 * ff1_26[k]
                 + pb_z[k] * fg_8[k];
    }
}

auto
compute_prim_fh_electron_repulsion_34(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_11 = buffer.data(dh + 11);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_8 = buffer.data(ff0 + 8);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_11 = buffer.data(ff1 + 11);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_14 = buffer.data(fg + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_4, \
                         dh_11, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_4[k];

        t_4[k] = pa_x[k] * dh_11[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_3, dh_4, dh_11, ff0_5, ff1_7, \
                         fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_3[k]
                 + f_1 * ff0_5[k]
                 - f_2 * ff1_7[k]
                 + pb_y[k] * fg_9[k];

        t_6[k] = pa_z[k] * dh_4[k];

        t_7[k] = pa_y[k] * dh_11[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_8, ff0_8, ff1_11, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_8[k]
                 + f_1 * ff0_8[k]
                 - f_2 * ff1_11[k]
                 + pb_z[k] * fg_14[k];
    }
}

auto
compute_prim_fh_electron_repulsion_35(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_11 = buffer.data(dh + 11);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_11 = buffer.data(ff0 + 11);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_41 = buffer.data(fg + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, dg_0, dg_1, dh_0, dh_1, \
                         ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = f_0 * dg_1[k]
                 + pa_x[k] * dh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, dg_2, dg_3, dg_4, dg_5, dh_2, \
                         dh_4, dh_5, dh_6, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dg_2[k]
                 + pa_x[k] * dh_2[k];

        t_5[k] = f_4 * dg_3[k]
                 + pb_x[k] * fg_14[k];

        t_6[k] = pa_x[k] * dh_4[k];

        t_7[k] = f_0 * dg_4[k]
                 + pa_x[k] * dh_5[k];

        t_8[k] = f_3 * dg_5[k]
                 + pa_x[k] * dh_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pa_z, pb_x, pb_y, dg_3, dg_8, dh_4, \
                         dh_11, ff0_7, ff1_15, fg_17, fg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dg_8[k]
                 + pb_x[k] * fg_17[k];

        t_10[k] = pa_x[k] * dh_11[k];

        t_11[k] = f_0 * dg_3[k]
                  + f_1 * ff0_7[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_24[k];

        t_12[k] = pa_z[k] * dh_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pb_y, dg_6, dg_7, dg_8, dh_8, dh_9, \
                         dh_11, fg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dg_6[k]
                  + pa_y[k] * dh_8[k];

        t_14[k] = f_3 * dg_7[k]
                  + pa_y[k] * dh_9[k];

        t_15[k] = f_4 * dg_8[k]
                  + pb_y[k] * fg_31[k];

        t_16[k] = pa_y[k] * dh_11[k];
    }

#pragma omp simd aligned(t_17, pb_z, dg_8, ff0_11, ff1_26, fg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * dg_8[k]
                  + f_1 * ff0_11[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_41[k];
    }
}

auto
compute_prim_fh_electron_repulsion_36(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 2.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_17 = buffer.data(dh + 17);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, \
                         ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_7[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];

        t_8[k] = pb_z[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, ff1_4, \
                         ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_10[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_12[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pa_y, pa_z, dg_1, dg_2, dg_3, dh_0, \
                         dh_1, dh_2, dh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_7 * dg_1[k]
                  + pa_y[k] * dh_1[k];

        t_14[k] = pa_z[k] * dh_0[k];

        t_15[k] = f_7 * dg_2[k]
                  + pa_z[k] * dh_2[k];

        t_16[k] = f_7 * dg_3[k]
                  + pa_x[k] * dh_3[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_x, pb_x, dg_4, dg_5, dg_6, dg_8, \
                         dh_4, dh_5, dh_7, dh_9, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * dg_4[k]
                  + pa_x[k] * dh_4[k];

        t_18[k] = f_8 * dg_5[k]
                  + pa_x[k] * dh_5[k];

        t_19[k] = f_9 * dg_6[k]
                  + pb_x[k] * fg_13[k];

        t_20[k] = pa_x[k] * dh_7[k];

        t_21[k] = f_7 * dg_8[k]
                  + pa_x[k] * dh_9[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, dg_9, dg_10, dg_14, dh_10, dh_11, \
                         dh_17, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_0 * dg_9[k]
                  + pa_x[k] * dh_10[k];

        t_23[k] = f_8 * dg_10[k]
                  + pa_x[k] * dh_11[k];

        t_24[k] = f_9 * dg_14[k]
                  + pb_x[k] * fg_16[k];

        t_25[k] = pa_x[k] * dh_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_x, pb_z, ff0_12, ff0_13, ff0_14, ff1_12, \
                         ff1_13, ff1_14, fg_17, fg_18, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_17[k];

        t_27[k] = pb_z[k] * fg_17[k];

        t_28[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_18[k];

        t_29[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_19[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_x, pb_y, pb_z, dg_6, ff0_15, ff0_17, \
                         ff1_15, ff1_17, fg_20, fg_21, fg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_20[k];

        t_31[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_21[k];

        t_32[k] = pb_x[k] * fg_22[k];

        t_33[k] = f_0 * dg_6[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_22[k];

        t_34[k] = pb_z[k] * fg_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_z, pb_z, dh_7, ff0_15, ff0_16, ff0_17, \
                         ff1_15, ff1_16, ff1_17, fg_23, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_23[k];

        t_36[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_24[k];

        t_37[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_25[k];

        t_38[k] = pa_z[k] * dh_7[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_y, pa_z, dg_7, dg_11, dg_12, dg_13, dh_8, \
                         dh_13, dh_14, dh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_7 * dg_7[k]
                  + pa_z[k] * dh_8[k];

        t_40[k] = f_7 * dg_11[k]
                  + pa_y[k] * dh_13[k];

        t_41[k] = f_0 * dg_12[k]
                  + pa_y[k] * dh_14[k];

        t_42[k] = f_8 * dg_13[k]
                  + pa_y[k] * dh_15[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_y, pb_x, pb_y, dg_14, dh_17, ff0_21, \
                         ff1_21, fg_29, fg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_9 * dg_14[k]
                  + pb_y[k] * fg_29[k];

        t_44[k] = pa_y[k] * dh_17[k];

        t_45[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_30[k];

        t_46[k] = pb_y[k] * fg_30[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pb_x, ff0_22, ff0_23, ff0_24, ff1_22, ff1_23, \
                         ff1_24, fg_31, fg_32, fg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_31[k];

        t_48[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_32[k];

        t_49[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_33[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pb_x, pb_y, ff0_24, ff0_25, ff0_26, ff1_24, \
                         ff1_25, ff1_26, fg_34, fg_35, fg_36, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_34[k];

        t_51[k] = pb_x[k] * fg_38[k];

        t_52[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_35[k];

        t_53[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_36[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_y, pb_z, dg_14, ff0_26, ff1_26, fg_37, \
                         fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_37[k];

        t_55[k] = pb_y[k] * fg_38[k];

        t_56[k] = f_0 * dg_14[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_37(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_11 = buffer.data(dh + 11);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, ff1_3, \
                         fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_7[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_9[k] = f_3 * ff0_5[k]
                 - f_4 * ff1_5[k]
                 + pb_y[k] * fg_7[k];

        t_10[k] = pb_y[k] * fg_8[k];

        t_11[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_12[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_x, pa_z, pb_x, dg_4, dg_5, dg_6, \
                         dh_0, dh_1, dh_2, dh_4, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_z[k] * dh_0[k];

        t_14[k] = f_0 * dg_4[k]
                  + pa_x[k] * dh_1[k];

        t_15[k] = f_7 * dg_5[k]
                  + pa_x[k] * dh_2[k];

        t_16[k] = f_8 * dg_6[k]
                  + pb_x[k] * fg_13[k];

        t_17[k] = pa_x[k] * dh_4[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pb_x, dg_9, dg_10, dg_14, dh_5, dh_6, \
                         dh_11, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * dg_9[k]
                  + pa_x[k] * dh_5[k];

        t_19[k] = f_7 * dg_10[k]
                  + pa_x[k] * dh_6[k];

        t_20[k] = f_8 * dg_14[k]
                  + pb_x[k] * fg_16[k];

        t_21[k] = pa_x[k] * dh_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pb_x, ff0_12, ff0_13, ff0_14, ff1_12, ff1_13, \
                         ff1_14, fg_17, fg_18, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_17[k];

        t_23[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_18[k];

        t_24[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_x, pb_y, dg_6, ff0_15, ff0_17, \
                         ff1_15, ff1_17, fg_20, fg_21, fg_22, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_20[k];

        t_26[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_21[k];

        t_27[k] = pb_x[k] * fg_22[k];

        t_28[k] = pb_x[k] * fg_25[k];

        t_29[k] = f_0 * dg_6[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_22[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_z, ff0_15, ff0_16, ff0_17, ff1_15, ff1_16, \
                         ff1_17, fg_22, fg_23, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_z[k] * fg_22[k];

        t_31[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_23[k];

        t_32[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_24[k];

        t_33[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_25[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_y, pa_z, pb_y, dg_12, dg_13, dg_14, \
                         dh_4, dh_8, dh_9, dh_11, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_z[k] * dh_4[k];

        t_35[k] = f_0 * dg_12[k]
                  + pa_y[k] * dh_8[k];

        t_36[k] = f_7 * dg_13[k]
                  + pa_y[k] * dh_9[k];

        t_37[k] = f_8 * dg_14[k]
                  + pb_y[k] * fg_29[k];

        t_38[k] = pa_y[k] * dh_11[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_x, ff0_21, ff0_22, ff0_23, ff1_21, ff1_22, \
                         ff1_23, fg_30, fg_31, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_30[k];

        t_40[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_31[k];

        t_41[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_32[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_x, pb_y, ff0_24, ff0_26, ff1_24, \
                         ff1_26, fg_33, fg_34, fg_35, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_33[k];

        t_43[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_34[k];

        t_44[k] = pb_x[k] * fg_35[k];

        t_45[k] = pb_x[k] * fg_38[k];

        t_46[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_35[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_y, pb_z, dg_14, ff0_25, ff0_26, ff1_25, \
                         ff1_26, fg_36, fg_37, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_36[k];

        t_48[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_37[k];

        t_49[k] = pb_y[k] * fg_38[k];

        t_50[k] = f_0 * dg_14[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_38(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_11 = buffer.data(dh + 11);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_11 = buffer.data(ff0 + 11);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_11 = buffer.data(ff1 + 11);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, dg_0, dg_1, dh_0, dh_1, \
                         ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = f_0 * dg_1[k]
                 + pa_x[k] * dh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, dg_2, dg_3, dg_4, dg_5, dh_2, \
                         dh_4, dh_5, dh_6, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dg_2[k]
                 + pa_x[k] * dh_2[k];

        t_5[k] = f_4 * dg_3[k]
                 + pb_x[k] * fg_5[k];

        t_6[k] = pa_x[k] * dh_4[k];

        t_7[k] = f_0 * dg_4[k]
                 + pa_x[k] * dh_5[k];

        t_8[k] = f_3 * dg_5[k]
                 + pa_x[k] * dh_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pa_z, pb_x, pb_y, dg_3, dg_8, dh_4, \
                         dh_11, ff0_7, ff1_7, fg_8, fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dg_8[k]
                 + pb_x[k] * fg_8[k];

        t_10[k] = pa_x[k] * dh_11[k];

        t_11[k] = f_0 * dg_3[k]
                  + f_1 * ff0_7[k]
                  - f_2 * ff1_7[k]
                  + pb_y[k] * fg_9[k];

        t_12[k] = pa_z[k] * dh_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pb_y, dg_6, dg_7, dg_8, dh_8, dh_9, \
                         dh_11, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dg_6[k]
                  + pa_y[k] * dh_8[k];

        t_14[k] = f_3 * dg_7[k]
                  + pa_y[k] * dh_9[k];

        t_15[k] = f_4 * dg_8[k]
                  + pb_y[k] * fg_13[k];

        t_16[k] = pa_y[k] * dh_11[k];
    }

#pragma omp simd aligned(t_17, pb_z, dg_8, ff0_11, ff1_11, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * dg_8[k]
                  + f_1 * ff0_11[k]
                  - f_2 * ff1_11[k]
                  + pb_z[k] * fg_14[k];
    }
}

auto
compute_prim_fh_electron_repulsion_39(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_11 = buffer.data(dh + 11);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_11 = buffer.data(ff0 + 11);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_38 = buffer.data(fg + 38);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, dg_0, dg_1, dh_0, dh_1, \
                         ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = f_0 * dg_1[k]
                 + pa_x[k] * dh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, dg_2, dg_3, dg_4, dg_5, dh_2, \
                         dh_4, dh_5, dh_6, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dg_2[k]
                 + pa_x[k] * dh_2[k];

        t_5[k] = f_4 * dg_3[k]
                 + pb_x[k] * fg_13[k];

        t_6[k] = pa_x[k] * dh_4[k];

        t_7[k] = f_0 * dg_4[k]
                 + pa_x[k] * dh_5[k];

        t_8[k] = f_3 * dg_5[k]
                 + pa_x[k] * dh_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pa_z, pb_x, pb_y, dg_3, dg_8, dh_4, \
                         dh_11, ff0_7, ff1_15, fg_16, fg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dg_8[k]
                 + pb_x[k] * fg_16[k];

        t_10[k] = pa_x[k] * dh_11[k];

        t_11[k] = f_0 * dg_3[k]
                  + f_1 * ff0_7[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_22[k];

        t_12[k] = pa_z[k] * dh_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pb_y, dg_6, dg_7, dg_8, dh_8, dh_9, \
                         dh_11, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dg_6[k]
                  + pa_y[k] * dh_8[k];

        t_14[k] = f_3 * dg_7[k]
                  + pa_y[k] * dh_9[k];

        t_15[k] = f_4 * dg_8[k]
                  + pb_y[k] * fg_29[k];

        t_16[k] = pa_y[k] * dh_11[k];
    }

#pragma omp simd aligned(t_17, pb_z, dg_8, ff0_11, ff1_26, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * dg_8[k]
                  + f_1 * ff0_11[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_40(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_17 = buffer.data(dh + 17);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, \
                         ff1_3, fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_7[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];

        t_8[k] = pb_z[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_10[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pb_y[k] * fg_7[k];

        t_11[k] = pb_y[k] * fg_8[k];

        t_12[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_13[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_x, pa_z, pb_x, dg_4, dg_5, dg_6, \
                         dh_0, dh_4, dh_5, dh_7, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * dh_0[k];

        t_15[k] = f_0 * dg_4[k]
                  + pa_x[k] * dh_4[k];

        t_16[k] = f_7 * dg_5[k]
                  + pa_x[k] * dh_5[k];

        t_17[k] = f_8 * dg_6[k]
                  + pb_x[k] * fg_13[k];

        t_18[k] = pa_x[k] * dh_7[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, dg_9, dg_10, dg_14, dh_10, dh_11, \
                         dh_17, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * dg_9[k]
                  + pa_x[k] * dh_10[k];

        t_20[k] = f_7 * dg_10[k]
                  + pa_x[k] * dh_11[k];

        t_21[k] = f_8 * dg_14[k]
                  + pb_x[k] * fg_16[k];

        t_22[k] = pa_x[k] * dh_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_x, pb_z, ff0_12, ff0_13, ff0_14, ff1_12, \
                         ff1_13, ff1_14, fg_17, fg_18, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_17[k];

        t_24[k] = pb_z[k] * fg_17[k];

        t_25[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_18[k];

        t_26[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_19[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pb_x, pb_y, dg_6, ff0_15, ff0_17, \
                         ff1_15, ff1_17, fg_20, fg_21, fg_22, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_20[k];

        t_28[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_21[k];

        t_29[k] = pb_x[k] * fg_22[k];

        t_30[k] = pb_x[k] * fg_25[k];

        t_31[k] = f_0 * dg_6[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_22[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_z, ff0_15, ff0_16, ff0_17, ff1_15, ff1_16, \
                         ff1_17, fg_22, fg_23, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_z[k] * fg_22[k];

        t_33[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_23[k];

        t_34[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_24[k];

        t_35[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_25[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pa_y, pa_z, pb_y, dg_12, dg_13, dg_14, \
                         dh_7, dh_14, dh_15, dh_17, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_z[k] * dh_7[k];

        t_37[k] = f_0 * dg_12[k]
                  + pa_y[k] * dh_14[k];

        t_38[k] = f_7 * dg_13[k]
                  + pa_y[k] * dh_15[k];

        t_39[k] = f_8 * dg_14[k]
                  + pb_y[k] * fg_29[k];

        t_40[k] = pa_y[k] * dh_17[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pb_x, pb_y, ff0_21, ff0_22, ff0_23, ff1_21, \
                         ff1_22, ff1_23, fg_30, fg_31, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_30[k];

        t_42[k] = pb_y[k] * fg_30[k];

        t_43[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_31[k];

        t_44[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_32[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pb_x, pb_y, ff0_24, ff0_26, ff1_24, \
                         ff1_26, fg_33, fg_34, fg_35, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_33[k];

        t_46[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_34[k];

        t_47[k] = pb_x[k] * fg_35[k];

        t_48[k] = pb_x[k] * fg_38[k];

        t_49[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_35[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pb_y, pb_z, dg_14, ff0_25, ff0_26, ff1_25, \
                         ff1_26, fg_36, fg_37, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_36[k];

        t_51[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_37[k];

        t_52[k] = pb_y[k] * fg_38[k];

        t_53[k] = f_0 * dg_14[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_41(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_17 = buffer.data(dh + 17);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_23 = buffer.data(ff0 + 23);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_23 = buffer.data(ff1 + 23);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_13 = buffer.data(fg + 13);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dg_0, ff0_0, ff1_0, fg_0, \
                         fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pb_y[k] * fg_0[k];

        t_2[k] = pb_z[k] * fg_0[k];

        t_3[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_y[k] * fg_1[k];

        t_4[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pb_y, pb_z, ff0_1, ff0_2, ff0_3, ff1_1, ff1_2, ff1_3, \
                         fg_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * ff0_1[k]
                 - f_6 * ff1_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = f_5 * ff0_2[k]
                 - f_6 * ff1_2[k]
                 + pb_z[k] * fg_4[k];

        t_7[k] = f_1 * ff0_3[k]
                 - f_2 * ff1_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, pa_y, pb_y, pb_z, dh_0, ff0_4, ff0_5, \
                         ff1_4, ff1_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * ff0_4[k]
                 - f_6 * ff1_4[k]
                 + pb_y[k] * fg_6[k];

        t_9[k] = f_3 * ff0_5[k]
                 - f_4 * ff1_5[k]
                 + pb_y[k] * fg_7[k];

        t_10[k] = pb_y[k] * fg_8[k];

        t_11[k] = f_1 * ff0_5[k]
                  - f_2 * ff1_5[k]
                  + pb_z[k] * fg_8[k];

        t_12[k] = pa_y[k] * dh_0[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_x, pa_z, pb_x, dg_4, dg_5, dg_6, \
                         dh_0, dh_4, dh_5, dh_7, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_z[k] * dh_0[k];

        t_14[k] = f_0 * dg_4[k]
                  + pa_x[k] * dh_4[k];

        t_15[k] = f_7 * dg_5[k]
                  + pa_x[k] * dh_5[k];

        t_16[k] = f_8 * dg_6[k]
                  + pb_x[k] * fg_13[k];

        t_17[k] = pa_x[k] * dh_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pb_x, dg_9, dg_10, dg_14, dh_10, dh_11, \
                         dh_17, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * dg_9[k]
                  + pa_x[k] * dh_10[k];

        t_19[k] = f_7 * dg_10[k]
                  + pa_x[k] * dh_11[k];

        t_20[k] = f_8 * dg_14[k]
                  + pb_x[k] * fg_16[k];

        t_21[k] = pa_x[k] * dh_17[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pb_x, ff0_12, ff0_13, ff0_14, ff1_12, ff1_13, \
                         ff1_14, fg_17, fg_18, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_x[k] * fg_17[k];

        t_23[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_x[k] * fg_18[k];

        t_24[k] = f_5 * ff0_14[k]
                  - f_6 * ff1_14[k]
                  + pb_x[k] * fg_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_x, pb_y, dg_6, ff0_15, ff0_17, \
                         ff1_15, ff1_17, fg_20, fg_21, fg_22, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_x[k] * fg_20[k];

        t_26[k] = f_3 * ff0_17[k]
                  - f_4 * ff1_17[k]
                  + pb_x[k] * fg_21[k];

        t_27[k] = pb_x[k] * fg_22[k];

        t_28[k] = pb_x[k] * fg_25[k];

        t_29[k] = f_0 * dg_6[k]
                  + f_1 * ff0_15[k]
                  - f_2 * ff1_15[k]
                  + pb_y[k] * fg_22[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_z, ff0_15, ff0_16, ff0_17, ff1_15, ff1_16, \
                         ff1_17, fg_22, fg_23, fg_24, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_z[k] * fg_22[k];

        t_31[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pb_z[k] * fg_23[k];

        t_32[k] = f_5 * ff0_16[k]
                  - f_6 * ff1_16[k]
                  + pb_z[k] * fg_24[k];

        t_33[k] = f_1 * ff0_17[k]
                  - f_2 * ff1_17[k]
                  + pb_z[k] * fg_25[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_y, pa_z, pb_y, dg_12, dg_13, dg_14, \
                         dh_7, dh_14, dh_15, dh_17, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_z[k] * dh_7[k];

        t_35[k] = f_0 * dg_12[k]
                  + pa_y[k] * dh_14[k];

        t_36[k] = f_7 * dg_13[k]
                  + pa_y[k] * dh_15[k];

        t_37[k] = f_8 * dg_14[k]
                  + pb_y[k] * fg_29[k];

        t_38[k] = pa_y[k] * dh_17[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_x, ff0_21, ff0_22, ff0_23, ff1_21, ff1_22, \
                         ff1_23, fg_30, fg_31, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_1 * ff0_21[k]
                  - f_2 * ff1_21[k]
                  + pb_x[k] * fg_30[k];

        t_40[k] = f_5 * ff0_22[k]
                  - f_6 * ff1_22[k]
                  + pb_x[k] * fg_31[k];

        t_41[k] = f_5 * ff0_23[k]
                  - f_6 * ff1_23[k]
                  + pb_x[k] * fg_32[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_x, pb_y, ff0_24, ff0_26, ff1_24, \
                         ff1_26, fg_33, fg_34, fg_35, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * ff0_24[k]
                  - f_4 * ff1_24[k]
                  + pb_x[k] * fg_33[k];

        t_43[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_x[k] * fg_34[k];

        t_44[k] = pb_x[k] * fg_35[k];

        t_45[k] = pb_x[k] * fg_38[k];

        t_46[k] = f_1 * ff0_24[k]
                  - f_2 * ff1_24[k]
                  + pb_y[k] * fg_35[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_y, pb_z, dg_14, ff0_25, ff0_26, ff1_25, \
                         ff1_26, fg_36, fg_37, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_5 * ff0_25[k]
                  - f_6 * ff1_25[k]
                  + pb_y[k] * fg_36[k];

        t_48[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pb_y[k] * fg_37[k];

        t_49[k] = pb_y[k] * fg_38[k];

        t_50[k] = f_0 * dg_14[k]
                  + f_1 * ff0_26[k]
                  - f_2 * ff1_26[k]
                  + pb_z[k] * fg_38[k];
    }
}

auto
compute_prim_fh_electron_repulsion_42(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_8 = buffer.data(ff0 + 8);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_14 = buffer.data(ff1 + 14);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_17 = buffer.data(fg + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_1, dh_1, dh_2, ff0_5, ff1_8, \
                         fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_1[k]
                 + f_1 * ff0_5[k]
                 - f_2 * ff1_8[k]
                 + pb_y[k] * fg_10[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_2, ff0_8, ff1_14, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_2[k]
                 + f_1 * ff0_8[k]
                 - f_2 * ff1_14[k]
                 + pb_z[k] * fg_17[k];
    }
}

auto
compute_prim_fh_electron_repulsion_43(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_5 = buffer.data(dg + 5);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_14 = buffer.data(ff0 + 14);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_14 = buffer.data(ff1 + 14);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_17 = buffer.data(fg + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_2, dh_1, dh_2, ff0_8, ff1_8, \
                         fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_2[k]
                 + f_1 * ff0_8[k]
                 - f_2 * ff1_8[k]
                 + pb_y[k] * fg_10[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_5, ff0_14, ff1_14, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_5[k]
                 + f_1 * ff0_14[k]
                 - f_2 * ff1_14[k]
                 + pb_z[k] * fg_17[k];
    }
}

auto
compute_prim_fh_electron_repulsion_44(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_5 = buffer.data(dg + 5);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_14 = buffer.data(ff0 + 14);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_20 = buffer.data(ff1 + 20);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_17 = buffer.data(fg + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_2, dh_1, dh_2, ff0_8, ff1_11, \
                         fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_2[k]
                 + f_1 * ff0_8[k]
                 - f_2 * ff1_11[k]
                 + pb_y[k] * fg_10[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_5, ff0_14, ff1_20, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_5[k]
                 + f_1 * ff0_14[k]
                 - f_2 * ff1_20[k]
                 + pb_z[k] * fg_17[k];
    }
}

auto
compute_prim_fh_electron_repulsion_45(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_5 = buffer.data(dg + 5);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_20 = buffer.data(ff0 + 20);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_14 = buffer.data(ff1 + 14);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_2, dh_1, dh_2, ff0_11, ff1_8, \
                         fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_2[k]
                 + f_1 * ff0_11[k]
                 - f_2 * ff1_8[k]
                 + pb_y[k] * fg_5[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_5, ff0_20, ff1_14, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_5[k]
                 + f_1 * ff0_20[k]
                 - f_2 * ff1_14[k]
                 + pb_z[k] * fg_8[k];
    }
}

auto
compute_prim_fh_electron_repulsion_46(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_8 = buffer.data(ff0 + 8);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_14 = buffer.data(ff1 + 14);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_20 = buffer.data(fg + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_1, dh_1, dh_2, ff0_5, ff1_8, \
                         fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_1[k]
                 + f_1 * ff0_5[k]
                 - f_2 * ff1_8[k]
                 + pb_y[k] * fg_11[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_2, ff0_8, ff1_14, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_2[k]
                 + f_1 * ff0_8[k]
                 - f_2 * ff1_14[k]
                 + pb_z[k] * fg_20[k];
    }
}

auto
compute_prim_fh_electron_repulsion_47(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / beta;
    const auto f_8 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, dg_0, dg_1, dh_0, dh_1, \
                         ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = f_0 * dg_1[k]
                 + pa_x[k] * dh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, dg_2, dg_3, dg_4, dg_5, dh_2, \
                         dh_3, dh_4, dh_5, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dg_2[k]
                 + pa_x[k] * dh_2[k];

        t_5[k] = f_4 * dg_3[k]
                 + pb_x[k] * fg_5[k];

        t_6[k] = pa_x[k] * dh_3[k];

        t_7[k] = f_0 * dg_4[k]
                 + pa_x[k] * dh_4[k];

        t_8[k] = f_3 * dg_5[k]
                 + pa_x[k] * dh_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pb_x, dg_8, dh_8, ff0_7, ff0_8, ff1_7, \
                         ff1_8, fg_8, fg_9, fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dg_8[k]
                 + pb_x[k] * fg_8[k];

        t_10[k] = pa_x[k] * dh_8[k];

        t_11[k] = f_5 * ff0_7[k]
                  - f_6 * ff1_7[k]
                  + pb_x[k] * fg_9[k];

        t_12[k] = f_7 * ff0_8[k]
                  - f_8 * ff1_8[k]
                  + pb_x[k] * fg_10[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pa_z, pb_y, dg_3, dg_6, dg_7, dh_3, \
                         dh_6, dh_7, ff0_8, ff1_8, fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dg_3[k]
                  + f_1 * ff0_8[k]
                  - f_2 * ff1_8[k]
                  + pb_y[k] * fg_11[k];

        t_14[k] = pa_z[k] * dh_3[k];

        t_15[k] = f_0 * dg_6[k]
                  + pa_y[k] * dh_6[k];

        t_16[k] = f_3 * dg_7[k]
                  + pa_y[k] * dh_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pb_x, pb_y, dg_8, dh_8, ff0_12, ff0_14, \
                         ff1_12, ff1_14, fg_15, fg_16, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_4 * dg_8[k]
                  + pb_y[k] * fg_15[k];

        t_18[k] = pa_y[k] * dh_8[k];

        t_19[k] = f_5 * ff0_12[k]
                  - f_6 * ff1_12[k]
                  + pb_x[k] * fg_16[k];

        t_20[k] = f_7 * ff0_14[k]
                  - f_8 * ff1_14[k]
                  + pb_x[k] * fg_17[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_y, pb_z, dg_8, ff0_13, ff0_14, ff1_13, ff1_14, \
                         fg_18, fg_19, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_y[k] * fg_18[k];

        t_22[k] = f_7 * ff0_14[k]
                  - f_8 * ff1_14[k]
                  + pb_y[k] * fg_19[k];

        t_23[k] = f_0 * dg_8[k]
                  + f_1 * ff0_14[k]
                  - f_2 * ff1_14[k]
                  + pb_z[k] * fg_20[k];
    }
}

auto
compute_prim_fh_electron_repulsion_48(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / beta;
    const auto f_8 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_10 = buffer.data(ff1 + 10);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_19 = buffer.data(ff1 + 19);
    const auto *ff1_20 = buffer.data(ff1 + 20);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, dg_0, dg_1, dh_0, dh_1, \
                         ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = f_0 * dg_1[k]
                 + pa_x[k] * dh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, dg_2, dg_3, dg_4, dg_5, dh_2, \
                         dh_3, dh_4, dh_5, fg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dg_2[k]
                 + pa_x[k] * dh_2[k];

        t_5[k] = f_4 * dg_3[k]
                 + pb_x[k] * fg_7[k];

        t_6[k] = pa_x[k] * dh_3[k];

        t_7[k] = f_0 * dg_4[k]
                 + pa_x[k] * dh_4[k];

        t_8[k] = f_3 * dg_5[k]
                 + pa_x[k] * dh_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pb_x, dg_8, dh_8, ff0_7, ff0_8, ff1_10, \
                         ff1_11, fg_10, fg_12, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dg_8[k]
                 + pb_x[k] * fg_10[k];

        t_10[k] = pa_x[k] * dh_8[k];

        t_11[k] = f_5 * ff0_7[k]
                  - f_6 * ff1_10[k]
                  + pb_x[k] * fg_12[k];

        t_12[k] = f_7 * ff0_8[k]
                  - f_8 * ff1_11[k]
                  + pb_x[k] * fg_13[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pa_z, pb_y, dg_3, dg_6, dg_7, dh_3, \
                         dh_6, dh_7, ff0_8, ff1_11, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dg_3[k]
                  + f_1 * ff0_8[k]
                  - f_2 * ff1_11[k]
                  + pb_y[k] * fg_14[k];

        t_14[k] = pa_z[k] * dh_3[k];

        t_15[k] = f_0 * dg_6[k]
                  + pa_y[k] * dh_6[k];

        t_16[k] = f_3 * dg_7[k]
                  + pa_y[k] * dh_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pb_x, pb_y, dg_8, dh_8, ff0_12, ff0_14, \
                         ff1_17, ff1_20, fg_19, fg_21, fg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_4 * dg_8[k]
                  + pb_y[k] * fg_19[k];

        t_18[k] = pa_y[k] * dh_8[k];

        t_19[k] = f_5 * ff0_12[k]
                  - f_6 * ff1_17[k]
                  + pb_x[k] * fg_21[k];

        t_20[k] = f_7 * ff0_14[k]
                  - f_8 * ff1_20[k]
                  + pb_x[k] * fg_22[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_y, pb_z, dg_8, ff0_13, ff0_14, ff1_19, ff1_20, \
                         fg_24, fg_25, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_19[k]
                  + pb_y[k] * fg_24[k];

        t_22[k] = f_7 * ff0_14[k]
                  - f_8 * ff1_20[k]
                  + pb_y[k] * fg_25[k];

        t_23[k] = f_0 * dg_8[k]
                  + f_1 * ff0_14[k]
                  - f_2 * ff1_20[k]
                  + pb_z[k] * fg_26[k];
    }
}

auto
compute_prim_fh_electron_repulsion_49(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / beta;
    const auto f_8 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_19 = buffer.data(ff0 + 19);
    const auto *ff0_20 = buffer.data(ff0 + 20);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_10 = buffer.data(ff1 + 10);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_19 = buffer.data(ff1 + 19);
    const auto *ff1_20 = buffer.data(ff1 + 20);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, dg_0, dg_1, dh_0, dh_1, \
                         ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = f_0 * dg_1[k]
                 + pa_x[k] * dh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, dg_2, dg_3, dg_4, dg_5, dh_2, \
                         dh_3, dh_4, dh_5, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dg_2[k]
                 + pa_x[k] * dh_2[k];

        t_5[k] = f_4 * dg_3[k]
                 + pb_x[k] * fg_5[k];

        t_6[k] = pa_x[k] * dh_3[k];

        t_7[k] = f_0 * dg_4[k]
                 + pa_x[k] * dh_4[k];

        t_8[k] = f_3 * dg_5[k]
                 + pa_x[k] * dh_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pb_x, dg_8, dh_8, ff0_10, ff0_11, \
                         ff1_10, ff1_11, fg_8, fg_9, fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dg_8[k]
                 + pb_x[k] * fg_8[k];

        t_10[k] = pa_x[k] * dh_8[k];

        t_11[k] = f_5 * ff0_10[k]
                  - f_6 * ff1_10[k]
                  + pb_x[k] * fg_9[k];

        t_12[k] = f_7 * ff0_11[k]
                  - f_8 * ff1_11[k]
                  + pb_x[k] * fg_10[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pa_z, pb_y, dg_3, dg_6, dg_7, dh_3, \
                         dh_6, dh_7, ff0_11, ff1_11, fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dg_3[k]
                  + f_1 * ff0_11[k]
                  - f_2 * ff1_11[k]
                  + pb_y[k] * fg_11[k];

        t_14[k] = pa_z[k] * dh_3[k];

        t_15[k] = f_0 * dg_6[k]
                  + pa_y[k] * dh_6[k];

        t_16[k] = f_3 * dg_7[k]
                  + pa_y[k] * dh_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pb_x, pb_y, dg_8, dh_8, ff0_17, ff0_20, \
                         ff1_17, ff1_20, fg_15, fg_16, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_4 * dg_8[k]
                  + pb_y[k] * fg_15[k];

        t_18[k] = pa_y[k] * dh_8[k];

        t_19[k] = f_5 * ff0_17[k]
                  - f_6 * ff1_17[k]
                  + pb_x[k] * fg_16[k];

        t_20[k] = f_7 * ff0_20[k]
                  - f_8 * ff1_20[k]
                  + pb_x[k] * fg_17[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_y, pb_z, dg_8, ff0_19, ff0_20, ff1_19, ff1_20, \
                         fg_18, fg_19, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * ff0_19[k]
                  - f_6 * ff1_19[k]
                  + pb_y[k] * fg_18[k];

        t_22[k] = f_7 * ff0_20[k]
                  - f_8 * ff1_20[k]
                  + pb_y[k] * fg_19[k];

        t_23[k] = f_0 * dg_8[k]
                  + f_1 * ff0_20[k]
                  - f_2 * ff1_20[k]
                  + pb_z[k] * fg_20[k];
    }
}

auto
compute_prim_fh_electron_repulsion_50(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_20 = buffer.data(ff0 + 20);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_14 = buffer.data(ff1 + 14);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_3, dh_1, dh_2, ff0_11, ff1_8, \
                         fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_3[k]
                 + f_1 * ff0_11[k]
                 - f_2 * ff1_8[k]
                 + pb_y[k] * fg_5[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_8, ff0_20, ff1_14, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_8[k]
                 + f_1 * ff0_20[k]
                 - f_2 * ff1_14[k]
                 + pb_z[k] * fg_8[k];
    }
}

auto
compute_prim_fh_electron_repulsion_51(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / beta;
    const auto f_8 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_11 = buffer.data(dh + 11);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_14 = buffer.data(ff1 + 14);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, dg_0, dg_1, dh_0, dh_1, \
                         ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = f_0 * dg_1[k]
                 + pa_x[k] * dh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, dg_2, dg_3, dg_4, dg_5, dh_2, \
                         dh_4, dh_5, dh_6, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dg_2[k]
                 + pa_x[k] * dh_2[k];

        t_5[k] = f_4 * dg_3[k]
                 + pb_x[k] * fg_5[k];

        t_6[k] = pa_x[k] * dh_4[k];

        t_7[k] = f_0 * dg_4[k]
                 + pa_x[k] * dh_5[k];

        t_8[k] = f_3 * dg_5[k]
                 + pa_x[k] * dh_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pb_x, dg_8, dh_11, ff0_7, ff0_8, ff1_7, \
                         ff1_8, fg_8, fg_9, fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dg_8[k]
                 + pb_x[k] * fg_8[k];

        t_10[k] = pa_x[k] * dh_11[k];

        t_11[k] = f_5 * ff0_7[k]
                  - f_6 * ff1_7[k]
                  + pb_x[k] * fg_9[k];

        t_12[k] = f_7 * ff0_8[k]
                  - f_8 * ff1_8[k]
                  + pb_x[k] * fg_10[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pa_z, pb_x, pb_y, dg_3, dg_6, dh_4, \
                         dh_8, ff0_8, ff1_8, fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_x[k] * fg_11[k];

        t_14[k] = f_0 * dg_3[k]
                  + f_1 * ff0_8[k]
                  - f_2 * ff1_8[k]
                  + pb_y[k] * fg_11[k];

        t_15[k] = pa_z[k] * dh_4[k];

        t_16[k] = f_0 * dg_6[k]
                  + pa_y[k] * dh_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pb_x, pb_y, dg_7, dg_8, dh_9, dh_11, \
                         ff0_12, ff1_12, fg_15, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * dg_7[k]
                  + pa_y[k] * dh_9[k];

        t_18[k] = f_4 * dg_8[k]
                  + pb_y[k] * fg_15[k];

        t_19[k] = pa_y[k] * dh_11[k];

        t_20[k] = f_5 * ff0_12[k]
                  - f_6 * ff1_12[k]
                  + pb_x[k] * fg_16[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_y, ff0_13, ff0_14, ff1_13, \
                         ff1_14, fg_17, fg_18, fg_19, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * ff0_14[k]
                  - f_8 * ff1_14[k]
                  + pb_x[k] * fg_17[k];

        t_22[k] = pb_x[k] * fg_20[k];

        t_23[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_13[k]
                  + pb_y[k] * fg_18[k];

        t_24[k] = f_7 * ff0_14[k]
                  - f_8 * ff1_14[k]
                  + pb_y[k] * fg_19[k];

        t_25[k] = pb_y[k] * fg_20[k];
    }

#pragma omp simd aligned(t_26, pb_z, dg_8, ff0_14, ff1_14, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * dg_8[k]
                  + f_1 * ff0_14[k]
                  - f_2 * ff1_14[k]
                  + pb_z[k] * fg_20[k];
    }
}

auto
compute_prim_fh_electron_repulsion_52(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / beta;
    const auto f_8 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_11 = buffer.data(dh + 11);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_14 = buffer.data(ff0 + 14);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_10 = buffer.data(ff1 + 10);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_19 = buffer.data(ff1 + 19);
    const auto *ff1_20 = buffer.data(ff1 + 20);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, dg_0, dg_1, dh_0, dh_1, \
                         ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = f_0 * dg_1[k]
                 + pa_x[k] * dh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, dg_2, dg_3, dg_4, dg_5, dh_2, \
                         dh_4, dh_5, dh_6, fg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dg_2[k]
                 + pa_x[k] * dh_2[k];

        t_5[k] = f_4 * dg_3[k]
                 + pb_x[k] * fg_7[k];

        t_6[k] = pa_x[k] * dh_4[k];

        t_7[k] = f_0 * dg_4[k]
                 + pa_x[k] * dh_5[k];

        t_8[k] = f_3 * dg_5[k]
                 + pa_x[k] * dh_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pb_x, dg_8, dh_11, ff0_7, ff0_8, ff1_10, \
                         ff1_11, fg_10, fg_12, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dg_8[k]
                 + pb_x[k] * fg_10[k];

        t_10[k] = pa_x[k] * dh_11[k];

        t_11[k] = f_5 * ff0_7[k]
                  - f_6 * ff1_10[k]
                  + pb_x[k] * fg_12[k];

        t_12[k] = f_7 * ff0_8[k]
                  - f_8 * ff1_11[k]
                  + pb_x[k] * fg_13[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pa_z, pb_x, pb_y, dg_3, dg_6, dh_4, \
                         dh_8, ff0_8, ff1_11, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_x[k] * fg_14[k];

        t_14[k] = f_0 * dg_3[k]
                  + f_1 * ff0_8[k]
                  - f_2 * ff1_11[k]
                  + pb_y[k] * fg_14[k];

        t_15[k] = pa_z[k] * dh_4[k];

        t_16[k] = f_0 * dg_6[k]
                  + pa_y[k] * dh_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pb_x, pb_y, dg_7, dg_8, dh_9, dh_11, \
                         ff0_12, ff1_17, fg_19, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * dg_7[k]
                  + pa_y[k] * dh_9[k];

        t_18[k] = f_4 * dg_8[k]
                  + pb_y[k] * fg_19[k];

        t_19[k] = pa_y[k] * dh_11[k];

        t_20[k] = f_5 * ff0_12[k]
                  - f_6 * ff1_17[k]
                  + pb_x[k] * fg_21[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_y, ff0_13, ff0_14, ff1_19, \
                         ff1_20, fg_22, fg_24, fg_25, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * ff0_14[k]
                  - f_8 * ff1_20[k]
                  + pb_x[k] * fg_22[k];

        t_22[k] = pb_x[k] * fg_26[k];

        t_23[k] = f_5 * ff0_13[k]
                  - f_6 * ff1_19[k]
                  + pb_y[k] * fg_24[k];

        t_24[k] = f_7 * ff0_14[k]
                  - f_8 * ff1_20[k]
                  + pb_y[k] * fg_25[k];

        t_25[k] = pb_y[k] * fg_26[k];
    }

#pragma omp simd aligned(t_26, pb_z, dg_8, ff0_14, ff1_20, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * dg_8[k]
                  + f_1 * ff0_14[k]
                  - f_2 * ff1_20[k]
                  + pb_z[k] * fg_26[k];
    }
}

auto
compute_prim_fh_electron_repulsion_53(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / beta;
    const auto f_8 = 0.5 * alpha / (beta * p);

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_11 = buffer.data(dh + 11);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_9 = buffer.data(ff0 + 9);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_18 = buffer.data(ff0 + 18);
    const auto *ff0_19 = buffer.data(ff0 + 19);
    const auto *ff0_20 = buffer.data(ff0 + 20);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_9 = buffer.data(ff1 + 9);
    const auto *ff1_10 = buffer.data(ff1 + 10);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_18 = buffer.data(ff1 + 18);
    const auto *ff1_19 = buffer.data(ff1 + 19);
    const auto *ff1_20 = buffer.data(ff1 + 20);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dg_0, ff0_0, ff0_1, ff0_2, ff1_0, \
                         ff1_1, ff1_2, fg_0, fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = f_1 * ff0_1[k]
                 - f_2 * ff1_1[k]
                 + pb_y[k] * fg_1[k];

        t_2[k] = f_1 * ff0_2[k]
                 - f_2 * ff1_2[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, pb_x, dg_1, dg_2, dg_3, \
                         dh_0, dh_1, dh_2, fg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * dh_0[k];

        t_4[k] = pa_z[k] * dh_0[k];

        t_5[k] = f_0 * dg_1[k]
                 + pa_x[k] * dh_1[k];

        t_6[k] = f_3 * dg_2[k]
                 + pa_x[k] * dh_2[k];

        t_7[k] = f_4 * dg_3[k]
                 + pb_x[k] * fg_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, pa_x, pb_x, dg_4, dg_5, dg_8, dh_4, dh_5, \
                         dh_6, dh_11, fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = pa_x[k] * dh_4[k];

        t_9[k] = f_0 * dg_4[k]
                 + pa_x[k] * dh_5[k];

        t_10[k] = f_3 * dg_5[k]
                  + pa_x[k] * dh_6[k];

        t_11[k] = f_4 * dg_8[k]
                  + pb_x[k] * fg_10[k];

        t_12[k] = pa_x[k] * dh_11[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pb_x, ff0_9, ff0_10, ff0_11, ff1_9, ff1_10, \
                         ff1_11, fg_11, fg_12, fg_13, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * ff0_9[k]
                  - f_2 * ff1_9[k]
                  + pb_x[k] * fg_11[k];

        t_14[k] = f_5 * ff0_10[k]
                  - f_6 * ff1_10[k]
                  + pb_x[k] * fg_12[k];

        t_15[k] = f_7 * ff0_11[k]
                  - f_8 * ff1_11[k]
                  + pb_x[k] * fg_13[k];

        t_16[k] = pb_x[k] * fg_14[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_z, pb_y, pb_z, dg_3, dh_4, ff0_11, ff0_12, \
                         ff1_11, ff1_12, fg_14, fg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * dg_3[k]
                  + f_1 * ff0_11[k]
                  - f_2 * ff1_11[k]
                  + pb_y[k] * fg_14[k];

        t_18[k] = f_1 * ff0_12[k]
                  - f_2 * ff1_12[k]
                  + pb_z[k] * fg_15[k];

        t_19[k] = pa_z[k] * dh_4[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_y, pb_y, dg_6, dg_7, dg_8, dh_8, dh_9, \
                         dh_11, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * dg_6[k]
                  + pa_y[k] * dh_8[k];

        t_21[k] = f_3 * dg_7[k]
                  + pa_y[k] * dh_9[k];

        t_22[k] = f_4 * dg_8[k]
                  + pb_y[k] * fg_19[k];

        t_23[k] = pa_y[k] * dh_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_x, ff0_16, ff0_17, ff0_20, ff1_16, ff1_17, \
                         ff1_20, fg_20, fg_21, fg_22, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_1 * ff0_16[k]
                  - f_2 * ff1_16[k]
                  + pb_x[k] * fg_20[k];

        t_25[k] = f_5 * ff0_17[k]
                  - f_6 * ff1_17[k]
                  + pb_x[k] * fg_21[k];

        t_26[k] = f_7 * ff0_20[k]
                  - f_8 * ff1_20[k]
                  + pb_x[k] * fg_22[k];

        t_27[k] = pb_x[k] * fg_26[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, ff0_18, ff0_19, ff0_20, ff1_18, ff1_19, \
                         ff1_20, fg_23, fg_24, fg_25, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * ff0_18[k]
                  - f_2 * ff1_18[k]
                  + pb_y[k] * fg_23[k];

        t_29[k] = f_5 * ff0_19[k]
                  - f_6 * ff1_19[k]
                  + pb_y[k] * fg_24[k];

        t_30[k] = f_7 * ff0_20[k]
                  - f_8 * ff1_20[k]
                  + pb_y[k] * fg_25[k];

        t_31[k] = pb_y[k] * fg_26[k];
    }

#pragma omp simd aligned(t_32, pb_z, dg_8, ff0_20, ff1_20, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * dg_8[k]
                  + f_1 * ff0_20[k]
                  - f_2 * ff1_20[k]
                  + pb_z[k] * fg_26[k];
    }
}

auto
compute_prim_fh_electron_repulsion_54(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / beta;
    const auto f_8 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_11 = buffer.data(dh + 11);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_19 = buffer.data(ff0 + 19);
    const auto *ff0_20 = buffer.data(ff0 + 20);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_10 = buffer.data(ff1 + 10);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_19 = buffer.data(ff1 + 19);
    const auto *ff1_20 = buffer.data(ff1 + 20);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, dg_0, dg_1, dh_0, dh_1, \
                         ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = f_0 * dg_1[k]
                 + pa_x[k] * dh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, dg_2, dg_3, dg_4, dg_5, dh_2, \
                         dh_4, dh_5, dh_6, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dg_2[k]
                 + pa_x[k] * dh_2[k];

        t_5[k] = f_4 * dg_3[k]
                 + pb_x[k] * fg_5[k];

        t_6[k] = pa_x[k] * dh_4[k];

        t_7[k] = f_0 * dg_4[k]
                 + pa_x[k] * dh_5[k];

        t_8[k] = f_3 * dg_5[k]
                 + pa_x[k] * dh_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pb_x, dg_8, dh_11, ff0_10, ff0_11, \
                         ff1_10, ff1_11, fg_8, fg_9, fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dg_8[k]
                 + pb_x[k] * fg_8[k];

        t_10[k] = pa_x[k] * dh_11[k];

        t_11[k] = f_5 * ff0_10[k]
                  - f_6 * ff1_10[k]
                  + pb_x[k] * fg_9[k];

        t_12[k] = f_7 * ff0_11[k]
                  - f_8 * ff1_11[k]
                  + pb_x[k] * fg_10[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pa_z, pb_x, pb_y, dg_3, dg_6, dh_4, \
                         dh_8, ff0_11, ff1_11, fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_x[k] * fg_11[k];

        t_14[k] = f_0 * dg_3[k]
                  + f_1 * ff0_11[k]
                  - f_2 * ff1_11[k]
                  + pb_y[k] * fg_11[k];

        t_15[k] = pa_z[k] * dh_4[k];

        t_16[k] = f_0 * dg_6[k]
                  + pa_y[k] * dh_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pb_x, pb_y, dg_7, dg_8, dh_9, dh_11, \
                         ff0_17, ff1_17, fg_15, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * dg_7[k]
                  + pa_y[k] * dh_9[k];

        t_18[k] = f_4 * dg_8[k]
                  + pb_y[k] * fg_15[k];

        t_19[k] = pa_y[k] * dh_11[k];

        t_20[k] = f_5 * ff0_17[k]
                  - f_6 * ff1_17[k]
                  + pb_x[k] * fg_16[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_y, ff0_19, ff0_20, ff1_19, \
                         ff1_20, fg_17, fg_18, fg_19, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * ff0_20[k]
                  - f_8 * ff1_20[k]
                  + pb_x[k] * fg_17[k];

        t_22[k] = pb_x[k] * fg_20[k];

        t_23[k] = f_5 * ff0_19[k]
                  - f_6 * ff1_19[k]
                  + pb_y[k] * fg_18[k];

        t_24[k] = f_7 * ff0_20[k]
                  - f_8 * ff1_20[k]
                  + pb_y[k] * fg_19[k];

        t_25[k] = pb_y[k] * fg_20[k];
    }

#pragma omp simd aligned(t_26, pb_z, dg_8, ff0_20, ff1_20, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * dg_8[k]
                  + f_1 * ff0_20[k]
                  - f_2 * ff1_20[k]
                  + pb_z[k] * fg_20[k];
    }
}

auto
compute_prim_fh_electron_repulsion_55(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_5 = buffer.data(dg + 5);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_11 = buffer.data(ff0 + 11);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_14 = buffer.data(ff1 + 14);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_17 = buffer.data(fg + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_2, dh_1, dh_2, ff0_7, ff1_8, \
                         fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_2[k]
                 + f_1 * ff0_7[k]
                 - f_2 * ff1_8[k]
                 + pb_y[k] * fg_10[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_5, ff0_11, ff1_14, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_5[k]
                 + f_1 * ff0_11[k]
                 - f_2 * ff1_14[k]
                 + pb_z[k] * fg_17[k];
    }
}

auto
compute_prim_fh_electron_repulsion_56(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_5 = buffer.data(dg + 5);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_14 = buffer.data(ff0 + 14);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_14 = buffer.data(ff1 + 14);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_2, dh_1, dh_2, ff0_8, ff1_8, \
                         fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_2[k]
                 + f_1 * ff0_8[k]
                 - f_2 * ff1_8[k]
                 + pb_y[k] * fg_5[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_5, ff0_14, ff1_14, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_5[k]
                 + f_1 * ff0_14[k]
                 - f_2 * ff1_14[k]
                 + pb_z[k] * fg_8[k];
    }
}

auto
compute_prim_fh_electron_repulsion_57(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_11 = buffer.data(ff0 + 11);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_14 = buffer.data(ff1 + 14);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_20 = buffer.data(fg + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, dg_0, dg_1, dh_0, dh_1, \
                         ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = f_0 * dg_1[k]
                 + pa_x[k] * dh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, dg_2, dg_3, dg_4, dg_5, dh_2, \
                         dh_3, dh_4, dh_5, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dg_2[k]
                 + pa_x[k] * dh_2[k];

        t_5[k] = f_4 * dg_3[k]
                 + pb_x[k] * fg_5[k];

        t_6[k] = pa_x[k] * dh_3[k];

        t_7[k] = f_0 * dg_4[k]
                 + pa_x[k] * dh_4[k];

        t_8[k] = f_3 * dg_5[k]
                 + pa_x[k] * dh_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pa_z, pb_x, pb_y, dg_3, dg_8, dh_3, \
                         dh_8, ff0_7, ff1_8, fg_8, fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dg_8[k]
                 + pb_x[k] * fg_8[k];

        t_10[k] = pa_x[k] * dh_8[k];

        t_11[k] = f_0 * dg_3[k]
                  + f_1 * ff0_7[k]
                  - f_2 * ff1_8[k]
                  + pb_y[k] * fg_11[k];

        t_12[k] = pa_z[k] * dh_3[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pb_y, dg_6, dg_7, dg_8, dh_6, dh_7, \
                         dh_8, fg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dg_6[k]
                  + pa_y[k] * dh_6[k];

        t_14[k] = f_3 * dg_7[k]
                  + pa_y[k] * dh_7[k];

        t_15[k] = f_4 * dg_8[k]
                  + pb_y[k] * fg_15[k];

        t_16[k] = pa_y[k] * dh_8[k];
    }

#pragma omp simd aligned(t_17, pb_z, dg_8, ff0_11, ff1_14, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * dg_8[k]
                  + f_1 * ff0_11[k]
                  - f_2 * ff1_14[k]
                  + pb_z[k] * fg_20[k];
    }
}

auto
compute_prim_fh_electron_repulsion_58(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_14 = buffer.data(ff0 + 14);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_14 = buffer.data(ff1 + 14);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, dg_0, dh_0, dh_1, \
                         dh_2, ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = pa_x[k] * dh_1[k];

        t_4[k] = pa_x[k] * dh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, dg_3, dh_1, dh_2, ff0_8, ff1_8, \
                         fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dg_3[k]
                 + f_1 * ff0_8[k]
                 - f_2 * ff1_8[k]
                 + pb_y[k] * fg_5[k];

        t_6[k] = pa_z[k] * dh_1[k];

        t_7[k] = pa_y[k] * dh_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, dg_8, ff0_14, ff1_14, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_8[k]
                 + f_1 * ff0_14[k]
                 - f_2 * ff1_14[k]
                 + pb_z[k] * fg_8[k];
    }
}

auto
compute_prim_fh_electron_repulsion_59(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t dg, const size_t dh,
                                      const size_t ff0, const size_t ff1, const size_t fg,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;

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

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_11 = buffer.data(dh + 11);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_11 = buffer.data(ff0 + 11);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_14 = buffer.data(ff1 + 14);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_20 = buffer.data(fg + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, dg_0, dg_1, dh_0, dh_1, \
                         ff0_0, ff1_0, fg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 + f_1 * ff0_0[k]
                 - f_2 * ff1_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = pa_y[k] * dh_0[k];

        t_2[k] = pa_z[k] * dh_0[k];

        t_3[k] = f_0 * dg_1[k]
                 + pa_x[k] * dh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, dg_2, dg_3, dg_4, dg_5, dh_2, \
                         dh_4, dh_5, dh_6, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dg_2[k]
                 + pa_x[k] * dh_2[k];

        t_5[k] = f_4 * dg_3[k]
                 + pb_x[k] * fg_5[k];

        t_6[k] = pa_x[k] * dh_4[k];

        t_7[k] = f_0 * dg_4[k]
                 + pa_x[k] * dh_5[k];

        t_8[k] = f_3 * dg_5[k]
                 + pa_x[k] * dh_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pa_z, pb_x, pb_y, dg_3, dg_8, dh_4, \
                         dh_11, ff0_7, ff1_8, fg_8, fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dg_8[k]
                 + pb_x[k] * fg_8[k];

        t_10[k] = pa_x[k] * dh_11[k];

        t_11[k] = f_0 * dg_3[k]
                  + f_1 * ff0_7[k]
                  - f_2 * ff1_8[k]
                  + pb_y[k] * fg_11[k];

        t_12[k] = pa_z[k] * dh_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pb_y, dg_6, dg_7, dg_8, dh_8, dh_9, \
                         dh_11, fg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dg_6[k]
                  + pa_y[k] * dh_8[k];

        t_14[k] = f_3 * dg_7[k]
                  + pa_y[k] * dh_9[k];

        t_15[k] = f_4 * dg_8[k]
                  + pb_y[k] * fg_15[k];

        t_16[k] = pa_y[k] * dh_11[k];
    }

#pragma omp simd aligned(t_17, pb_z, dg_8, ff0_11, ff1_14, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * dg_8[k]
                  + f_1 * ff0_11[k]
                  - f_2 * ff1_14[k]
                  + pb_z[k] * fg_20[k];
    }
}

}  // namespace simdt2ceri
