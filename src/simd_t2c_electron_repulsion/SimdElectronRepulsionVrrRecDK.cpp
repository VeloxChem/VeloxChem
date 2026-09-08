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


#include "SimdElectronRepulsionVrrRecDK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_dk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pi, const size_t pk,
                                     const size_t dh0, const size_t dh1, const size_t di,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 2.5 / p;
    const auto f_13 = 2.0 / p;
    const auto f_14 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_11 = buffer.data(pi + 11);
    const auto *pi_12 = buffer.data(pi + 12);
    const auto *pi_13 = buffer.data(pi + 13);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_16 = buffer.data(pi + 16);
    const auto *pi_17 = buffer.data(pi + 17);
    const auto *pi_18 = buffer.data(pi + 18);
    const auto *pi_19 = buffer.data(pi + 19);
    const auto *pi_20 = buffer.data(pi + 20);
    const auto *pi_21 = buffer.data(pi + 21);
    const auto *pi_22 = buffer.data(pi + 22);
    const auto *pi_23 = buffer.data(pi + 23);
    const auto *pi_24 = buffer.data(pi + 24);
    const auto *pi_25 = buffer.data(pi + 25);
    const auto *pi_26 = buffer.data(pi + 26);
    const auto *pi_27 = buffer.data(pi + 27);
    const auto *pi_28 = buffer.data(pi + 28);
    const auto *pi_29 = buffer.data(pi + 29);
    const auto *pi_30 = buffer.data(pi + 30);
    const auto *pi_31 = buffer.data(pi + 31);
    const auto *pi_32 = buffer.data(pi + 32);
    const auto *pi_33 = buffer.data(pi + 33);
    const auto *pi_34 = buffer.data(pi + 34);
    const auto *pi_35 = buffer.data(pi + 35);
    const auto *pi_36 = buffer.data(pi + 36);
    const auto *pi_37 = buffer.data(pi + 37);
    const auto *pi_38 = buffer.data(pi + 38);
    const auto *pi_39 = buffer.data(pi + 39);
    const auto *pi_40 = buffer.data(pi + 40);
    const auto *pi_41 = buffer.data(pi + 41);
    const auto *pi_42 = buffer.data(pi + 42);
    const auto *pi_43 = buffer.data(pi + 43);
    const auto *pi_44 = buffer.data(pi + 44);
    const auto *pi_45 = buffer.data(pi + 45);
    const auto *pi_46 = buffer.data(pi + 46);
    const auto *pi_47 = buffer.data(pi + 47);
    const auto *pi_48 = buffer.data(pi + 48);
    const auto *pi_49 = buffer.data(pi + 49);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);
    const auto *pk_3 = buffer.data(pk + 3);
    const auto *pk_4 = buffer.data(pk + 4);
    const auto *pk_5 = buffer.data(pk + 5);
    const auto *pk_6 = buffer.data(pk + 6);
    const auto *pk_7 = buffer.data(pk + 7);
    const auto *pk_8 = buffer.data(pk + 8);
    const auto *pk_9 = buffer.data(pk + 9);
    const auto *pk_10 = buffer.data(pk + 10);
    const auto *pk_11 = buffer.data(pk + 11);
    const auto *pk_12 = buffer.data(pk + 12);
    const auto *pk_13 = buffer.data(pk + 13);
    const auto *pk_14 = buffer.data(pk + 14);
    const auto *pk_15 = buffer.data(pk + 15);
    const auto *pk_16 = buffer.data(pk + 16);
    const auto *pk_17 = buffer.data(pk + 17);
    const auto *pk_18 = buffer.data(pk + 18);
    const auto *pk_19 = buffer.data(pk + 19);
    const auto *pk_20 = buffer.data(pk + 20);
    const auto *pk_21 = buffer.data(pk + 21);
    const auto *pk_22 = buffer.data(pk + 22);
    const auto *pk_23 = buffer.data(pk + 23);
    const auto *pk_24 = buffer.data(pk + 24);
    const auto *pk_25 = buffer.data(pk + 25);
    const auto *pk_26 = buffer.data(pk + 26);
    const auto *pk_27 = buffer.data(pk + 27);
    const auto *pk_28 = buffer.data(pk + 28);
    const auto *pk_29 = buffer.data(pk + 29);
    const auto *pk_30 = buffer.data(pk + 30);
    const auto *pk_31 = buffer.data(pk + 31);
    const auto *pk_32 = buffer.data(pk + 32);
    const auto *pk_33 = buffer.data(pk + 33);
    const auto *pk_34 = buffer.data(pk + 34);
    const auto *pk_35 = buffer.data(pk + 35);
    const auto *pk_36 = buffer.data(pk + 36);
    const auto *pk_37 = buffer.data(pk + 37);
    const auto *pk_38 = buffer.data(pk + 38);
    const auto *pk_39 = buffer.data(pk + 39);
    const auto *pk_40 = buffer.data(pk + 40);
    const auto *pk_41 = buffer.data(pk + 41);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);
    const auto *dh0_3 = buffer.data(dh0 + 3);
    const auto *dh0_4 = buffer.data(dh0 + 4);
    const auto *dh0_5 = buffer.data(dh0 + 5);
    const auto *dh0_6 = buffer.data(dh0 + 6);
    const auto *dh0_7 = buffer.data(dh0 + 7);
    const auto *dh0_8 = buffer.data(dh0 + 8);
    const auto *dh0_9 = buffer.data(dh0 + 9);
    const auto *dh0_10 = buffer.data(dh0 + 10);
    const auto *dh0_11 = buffer.data(dh0 + 11);
    const auto *dh0_12 = buffer.data(dh0 + 12);
    const auto *dh0_13 = buffer.data(dh0 + 13);
    const auto *dh0_14 = buffer.data(dh0 + 14);
    const auto *dh0_15 = buffer.data(dh0 + 15);
    const auto *dh0_16 = buffer.data(dh0 + 16);
    const auto *dh0_17 = buffer.data(dh0 + 17);
    const auto *dh0_18 = buffer.data(dh0 + 18);
    const auto *dh0_19 = buffer.data(dh0 + 19);
    const auto *dh0_20 = buffer.data(dh0 + 20);
    const auto *dh0_21 = buffer.data(dh0 + 21);
    const auto *dh0_22 = buffer.data(dh0 + 22);
    const auto *dh0_23 = buffer.data(dh0 + 23);
    const auto *dh0_24 = buffer.data(dh0 + 24);
    const auto *dh0_25 = buffer.data(dh0 + 25);
    const auto *dh0_26 = buffer.data(dh0 + 26);
    const auto *dh0_27 = buffer.data(dh0 + 27);
    const auto *dh0_28 = buffer.data(dh0 + 28);
    const auto *dh0_29 = buffer.data(dh0 + 29);
    const auto *dh0_30 = buffer.data(dh0 + 30);
    const auto *dh0_31 = buffer.data(dh0 + 31);
    const auto *dh0_32 = buffer.data(dh0 + 32);
    const auto *dh0_33 = buffer.data(dh0 + 33);
    const auto *dh0_34 = buffer.data(dh0 + 34);
    const auto *dh0_35 = buffer.data(dh0 + 35);
    const auto *dh0_36 = buffer.data(dh0 + 36);
    const auto *dh0_37 = buffer.data(dh0 + 37);
    const auto *dh0_38 = buffer.data(dh0 + 38);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);
    const auto *dh1_3 = buffer.data(dh1 + 3);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_5 = buffer.data(dh1 + 5);
    const auto *dh1_6 = buffer.data(dh1 + 6);
    const auto *dh1_7 = buffer.data(dh1 + 7);
    const auto *dh1_8 = buffer.data(dh1 + 8);
    const auto *dh1_9 = buffer.data(dh1 + 9);
    const auto *dh1_10 = buffer.data(dh1 + 10);
    const auto *dh1_11 = buffer.data(dh1 + 11);
    const auto *dh1_12 = buffer.data(dh1 + 12);
    const auto *dh1_13 = buffer.data(dh1 + 13);
    const auto *dh1_14 = buffer.data(dh1 + 14);
    const auto *dh1_15 = buffer.data(dh1 + 15);
    const auto *dh1_16 = buffer.data(dh1 + 16);
    const auto *dh1_17 = buffer.data(dh1 + 17);
    const auto *dh1_18 = buffer.data(dh1 + 18);
    const auto *dh1_19 = buffer.data(dh1 + 19);
    const auto *dh1_20 = buffer.data(dh1 + 20);
    const auto *dh1_21 = buffer.data(dh1 + 21);
    const auto *dh1_22 = buffer.data(dh1 + 22);
    const auto *dh1_23 = buffer.data(dh1 + 23);
    const auto *dh1_24 = buffer.data(dh1 + 24);
    const auto *dh1_25 = buffer.data(dh1 + 25);
    const auto *dh1_26 = buffer.data(dh1 + 26);
    const auto *dh1_27 = buffer.data(dh1 + 27);
    const auto *dh1_28 = buffer.data(dh1 + 28);
    const auto *dh1_29 = buffer.data(dh1 + 29);
    const auto *dh1_30 = buffer.data(dh1 + 30);
    const auto *dh1_31 = buffer.data(dh1 + 31);
    const auto *dh1_32 = buffer.data(dh1 + 32);
    const auto *dh1_33 = buffer.data(dh1 + 33);
    const auto *dh1_34 = buffer.data(dh1 + 34);
    const auto *dh1_35 = buffer.data(dh1 + 35);
    const auto *dh1_36 = buffer.data(dh1 + 36);
    const auto *dh1_37 = buffer.data(dh1 + 37);
    const auto *dh1_38 = buffer.data(dh1 + 38);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, pi_0, dh0_0, dh1_0, \
                         di_0, di_1, di_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = pb_y[k] * di_0[k];

        t_2[k] = pb_z[k] * di_0[k];

        t_3[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_y[k] * di_1[k];

        t_4[k] = pb_y[k] * di_2[k];

        t_5[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_z[k] * di_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, dh0_1, dh0_2, dh0_3, dh1_1, \
                         dh1_2, dh1_3, di_3, di_4, di_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * dh0_1[k]
                 - f_6 * dh1_1[k]
                 + pb_y[k] * di_3[k];

        t_7[k] = pb_z[k] * di_3[k];

        t_8[k] = pb_y[k] * di_4[k];

        t_9[k] = f_5 * dh0_2[k]
                 - f_6 * dh1_2[k]
                 + pb_z[k] * di_4[k];

        t_10[k] = f_7 * dh0_3[k]
                  - f_8 * dh1_3[k]
                  + pb_y[k] * di_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, dh0_4, dh0_5, dh1_4, \
                         dh1_5, di_5, di_6, di_7, di_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * di_5[k];

        t_12[k] = f_3 * dh0_4[k]
                  - f_4 * dh1_4[k]
                  + pb_y[k] * di_6[k];

        t_13[k] = pb_y[k] * di_7[k];

        t_14[k] = f_7 * dh0_4[k]
                  - f_8 * dh1_4[k]
                  + pb_z[k] * di_7[k];

        t_15[k] = f_9 * dh0_5[k]
                  - f_10 * dh1_5[k]
                  + pb_y[k] * di_8[k];

        t_16[k] = pb_z[k] * di_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, dh0_6, dh0_7, dh1_6, dh1_7, di_9, \
                         di_10, di_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * dh0_6[k]
                  - f_6 * dh1_6[k]
                  + pb_y[k] * di_9[k];

        t_18[k] = f_3 * dh0_7[k]
                  - f_4 * dh1_7[k]
                  + pb_y[k] * di_10[k];

        t_19[k] = pb_y[k] * di_11[k];

        t_20[k] = f_9 * dh0_7[k]
                  - f_10 * dh1_7[k]
                  + pb_z[k] * di_11[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_z, pi_7, pi_8, pi_9, pi_10, \
                         di_12, di_14, di_15, di_16, di_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * pi_7[k]
                  + pb_x[k] * di_14[k];

        t_22[k] = pb_z[k] * di_12[k];

        t_23[k] = f_0 * pi_8[k]
                  + pb_x[k] * di_15[k];

        t_24[k] = f_0 * pi_9[k]
                  + pb_x[k] * di_16[k];

        t_25[k] = f_0 * pi_10[k]
                  + pb_x[k] * di_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, pi_11, dh0_8, dh1_8, di_13, \
                         di_14, di_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_y[k] * di_13[k];

        t_27[k] = f_0 * pi_11[k]
                  + pb_x[k] * di_19[k];

        t_28[k] = f_1 * dh0_8[k]
                  - f_2 * dh1_8[k]
                  + pb_y[k] * di_14[k];

        t_29[k] = pb_z[k] * di_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_y, dh0_9, dh0_10, dh0_11, dh1_9, dh1_10, dh1_11, \
                         di_15, di_16, di_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_9 * dh0_9[k]
                  - f_10 * dh1_9[k]
                  + pb_y[k] * di_15[k];

        t_31[k] = f_7 * dh0_10[k]
                  - f_8 * dh1_10[k]
                  + pb_y[k] * di_16[k];

        t_32[k] = f_5 * dh0_11[k]
                  - f_6 * dh1_11[k]
                  + pb_y[k] * di_17[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pa_y, pb_y, pb_z, pi_0, pk_0, \
                         dh0_12, dh1_12, di_18, di_19, di_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * dh0_12[k]
                  - f_4 * dh1_12[k]
                  + pb_y[k] * di_18[k];

        t_34[k] = pb_y[k] * di_19[k];

        t_35[k] = f_1 * dh0_12[k]
                  - f_2 * dh1_12[k]
                  + pb_z[k] * di_19[k];

        t_36[k] = pa_y[k] * pk_0[k];

        t_37[k] = f_11 * pi_0[k]
                  + pb_y[k] * di_20[k];

        t_38[k] = pb_z[k] * di_20[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_x, pa_y, pb_z, pi_13, pi_15, pk_2, \
                         pk_12, pk_13, di_21, di_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * pi_13[k]
                  + pa_x[k] * pk_12[k];

        t_40[k] = pb_z[k] * di_21[k];

        t_41[k] = pa_y[k] * pk_2[k];

        t_42[k] = f_13 * pi_15[k]
                  + pa_x[k] * pk_13[k];

        t_43[k] = pb_z[k] * di_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pa_y, pb_y, pb_z, pi_2, pi_17, pk_4, \
                         pk_14, di_23, di_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_11 * pi_2[k]
                  + pb_y[k] * di_23[k];

        t_45[k] = pa_y[k] * pk_4[k];

        t_46[k] = f_14 * pi_17[k]
                  + pa_x[k] * pk_14[k];

        t_47[k] = pb_z[k] * di_24[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pa_y, pb_y, pi_4, pi_18, pi_20, pk_6, \
                         pk_15, pk_16, di_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_14 * pi_18[k]
                  + pa_x[k] * pk_15[k];

        t_49[k] = f_11 * pi_4[k]
                  + pb_y[k] * di_25[k];

        t_50[k] = pa_y[k] * pk_6[k];

        t_51[k] = f_0 * pi_20[k]
                  + pa_x[k] * pk_16[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_y, pb_z, pi_6, pi_21, pi_22, pk_17, \
                         pk_18, di_26, di_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_z[k] * di_26[k];

        t_53[k] = f_0 * pi_21[k]
                  + pa_x[k] * pk_17[k];

        t_54[k] = f_0 * pi_22[k]
                  + pa_x[k] * pk_18[k];

        t_55[k] = f_11 * pi_6[k]
                  + pb_y[k] * di_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_y, pb_x, pb_z, pi_23, pi_24, pi_25, \
                         pk_8, di_28, di_29, di_30, di_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_y[k] * pk_8[k];

        t_57[k] = f_11 * pi_23[k]
                  + pb_x[k] * di_29[k];

        t_58[k] = pb_z[k] * di_28[k];

        t_59[k] = f_11 * pi_24[k]
                  + pb_x[k] * di_30[k];

        t_60[k] = f_11 * pi_25[k]
                  + pb_x[k] * di_31[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pa_x, pa_y, pb_x, pb_z, pi_26, pi_27, \
                         pk_10, pk_19, di_29, di_32, di_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_11 * pi_26[k]
                  + pb_x[k] * di_32[k];

        t_62[k] = f_11 * pi_27[k]
                  + pb_x[k] * di_33[k];

        t_63[k] = pa_y[k] * pk_10[k];

        t_64[k] = pa_x[k] * pk_19[k];

        t_65[k] = pb_z[k] * di_29[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, t_72, pa_x, pa_z, pk_0, pk_20, \
                         pk_21, pk_22, pk_23, pk_24, pk_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_x[k] * pk_20[k];

        t_67[k] = pa_x[k] * pk_21[k];

        t_68[k] = pa_x[k] * pk_22[k];

        t_69[k] = pa_x[k] * pk_23[k];

        t_70[k] = pa_x[k] * pk_24[k];

        t_71[k] = pa_x[k] * pk_25[k];

        t_72[k] = pa_z[k] * pk_0[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_x, pa_z, pb_y, pb_z, pi_0, pi_32, \
                         pk_1, pk_28, di_34, di_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pb_y[k] * di_34[k];

        t_74[k] = f_11 * pi_0[k]
                  + pb_z[k] * di_34[k];

        t_75[k] = pa_z[k] * pk_1[k];

        t_76[k] = pb_y[k] * di_35[k];

        t_77[k] = f_12 * pi_32[k]
                  + pa_x[k] * pk_28[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_x, pa_z, pb_y, pb_z, pi_1, pi_35, \
                         pk_3, pk_5, pk_29, di_36, di_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_z[k] * pk_3[k];

        t_79[k] = f_11 * pi_1[k]
                  + pb_z[k] * di_36[k];

        t_80[k] = pb_y[k] * di_37[k];

        t_81[k] = f_13 * pi_35[k]
                  + pa_x[k] * pk_29[k];

        t_82[k] = pa_z[k] * pk_5[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_x, pb_y, pb_z, pi_3, pi_37, pi_39, pk_30, \
                         pk_31, di_38, di_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_11 * pi_3[k]
                  + pb_z[k] * di_38[k];

        t_84[k] = f_14 * pi_37[k]
                  + pa_x[k] * pk_30[k];

        t_85[k] = pb_y[k] * di_39[k];

        t_86[k] = f_14 * pi_39[k]
                  + pa_x[k] * pk_31[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_x, pa_z, pb_z, pi_5, pi_40, pi_41, pk_7, \
                         pk_32, pk_33, di_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * pk_7[k];

        t_88[k] = f_11 * pi_5[k]
                  + pb_z[k] * di_40[k];

        t_89[k] = f_0 * pi_40[k]
                  + pa_x[k] * pk_32[k];

        t_90[k] = f_0 * pi_41[k]
                  + pa_x[k] * pk_33[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_x, pa_z, pb_x, pb_y, pi_42, pi_44, pk_9, \
                         pk_34, di_41, di_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = pb_y[k] * di_41[k];

        t_92[k] = f_0 * pi_42[k]
                  + pa_x[k] * pk_34[k];

        t_93[k] = pa_z[k] * pk_9[k];

        t_94[k] = f_11 * pi_44[k]
                  + pb_x[k] * di_43[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pb_x, pb_y, pi_45, pi_46, pi_47, pi_49, \
                         di_42, di_44, di_45, di_46, di_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_11 * pi_45[k]
                  + pb_x[k] * di_44[k];

        t_96[k] = f_11 * pi_46[k]
                  + pb_x[k] * di_45[k];

        t_97[k] = f_11 * pi_47[k]
                  + pb_x[k] * di_46[k];

        t_98[k] = pb_y[k] * di_42[k];

        t_99[k] = f_11 * pi_49[k]
                  + pb_x[k] * di_47[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, t_105, t_106, pa_x, pb_y, pk_35, \
                         pk_36, pk_37, pk_38, pk_39, pk_40, di_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_x[k] * pk_35[k];

        t_101[k] = pa_x[k] * pk_36[k];

        t_102[k] = pa_x[k] * pk_37[k];

        t_103[k] = pa_x[k] * pk_38[k];

        t_104[k] = pa_x[k] * pk_39[k];

        t_105[k] = pa_x[k] * pk_40[k];

        t_106[k] = pb_y[k] * di_47[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_x, pb_x, pb_y, pb_z, pi_12, pk_41, \
                         dh0_13, dh1_13, di_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pa_x[k] * pk_41[k];

        t_108[k] = f_1 * dh0_13[k]
                   - f_2 * dh1_13[k]
                   + pb_x[k] * di_48[k];

        t_109[k] = f_0 * pi_12[k]
                   + pb_y[k] * di_48[k];

        t_110[k] = pb_z[k] * di_48[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pb_x, pb_z, dh0_14, dh0_15, dh0_16, \
                         dh1_14, dh1_15, dh1_16, di_49, di_50, di_51, \
                         di_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_9 * dh0_14[k]
                   - f_10 * dh1_14[k]
                   + pb_x[k] * di_50[k];

        t_112[k] = pb_z[k] * di_49[k];

        t_113[k] = f_9 * dh0_15[k]
                   - f_10 * dh1_15[k]
                   + pb_x[k] * di_51[k];

        t_114[k] = f_7 * dh0_16[k]
                   - f_8 * dh1_16[k]
                   + pb_x[k] * di_52[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pb_x, pb_y, pb_z, pi_14, dh0_17, dh0_18, \
                         dh1_17, dh1_18, di_50, di_51, di_53, di_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pb_z[k] * di_50[k];

        t_116[k] = f_0 * pi_14[k]
                   + pb_y[k] * di_51[k];

        t_117[k] = f_7 * dh0_17[k]
                   - f_8 * dh1_17[k]
                   + pb_x[k] * di_53[k];

        t_118[k] = f_5 * dh0_18[k]
                   - f_6 * dh1_18[k]
                   + pb_x[k] * di_54[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pb_x, pb_y, pb_z, pi_16, dh0_19, dh0_20, \
                         dh1_19, dh1_20, di_52, di_53, di_55, di_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pb_z[k] * di_52[k];

        t_120[k] = f_5 * dh0_19[k]
                   - f_6 * dh1_19[k]
                   + pb_x[k] * di_55[k];

        t_121[k] = f_0 * pi_16[k]
                   + pb_y[k] * di_53[k];

        t_122[k] = f_5 * dh0_20[k]
                   - f_6 * dh1_20[k]
                   + pb_x[k] * di_56[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pb_x, pb_z, dh0_21, dh0_23, dh0_24, \
                         dh1_21, dh1_23, dh1_24, di_54, di_57, di_58, \
                         di_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_3 * dh0_21[k]
                   - f_4 * dh1_21[k]
                   + pb_x[k] * di_57[k];

        t_124[k] = pb_z[k] * di_54[k];

        t_125[k] = f_3 * dh0_23[k]
                   - f_4 * dh1_23[k]
                   + pb_x[k] * di_58[k];

        t_126[k] = f_3 * dh0_24[k]
                   - f_4 * dh1_24[k]
                   + pb_x[k] * di_59[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pb_x, pb_y, pi_19, dh0_25, dh1_25, \
                         di_56, di_60, di_61, di_62, di_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_0 * pi_19[k]
                   + pb_y[k] * di_56[k];

        t_128[k] = f_3 * dh0_25[k]
                   - f_4 * dh1_25[k]
                   + pb_x[k] * di_60[k];

        t_129[k] = pb_x[k] * di_61[k];

        t_130[k] = pb_x[k] * di_62[k];

        t_131[k] = pb_x[k] * di_63[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, pb_x, pb_y, pi_23, dh0_21, dh1_21, \
                         di_61, di_64, di_65, di_66, di_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pb_x[k] * di_64[k];

        t_133[k] = pb_x[k] * di_65[k];

        t_134[k] = pb_x[k] * di_66[k];

        t_135[k] = pb_x[k] * di_67[k];

        t_136[k] = f_0 * pi_23[k]
                   + f_1 * dh0_21[k]
                   - f_2 * dh1_21[k]
                   + pb_y[k] * di_61[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pb_z, dh0_21, dh0_22, dh0_23, dh1_21, \
                         dh1_22, dh1_23, di_61, di_62, di_63, di_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = pb_z[k] * di_61[k];

        t_138[k] = f_3 * dh0_21[k]
                   - f_4 * dh1_21[k]
                   + pb_z[k] * di_62[k];

        t_139[k] = f_5 * dh0_22[k]
                   - f_6 * dh1_22[k]
                   + pb_z[k] * di_63[k];

        t_140[k] = f_7 * dh0_23[k]
                   - f_8 * dh1_23[k]
                   + pb_z[k] * di_64[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_y, pb_y, pb_z, pi_28, pk_26, dh0_24, \
                         dh0_25, dh1_24, dh1_25, di_65, di_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_9 * dh0_24[k]
                   - f_10 * dh1_24[k]
                   + pb_z[k] * di_65[k];

        t_142[k] = f_0 * pi_28[k]
                   + pb_y[k] * di_67[k];

        t_143[k] = f_1 * dh0_25[k]
                   - f_2 * dh1_25[k]
                   + pb_z[k] * di_67[k];

        t_144[k] = pa_y[k] * pk_26[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, t_150, pa_y, pa_z, pb_y, pi_30, \
                         pk_11, pk_12, pk_13, pk_27, pk_28, di_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = pa_z[k] * pk_11[k];

        t_146[k] = pa_y[k] * pk_27[k];

        t_147[k] = pa_z[k] * pk_12[k];

        t_148[k] = f_11 * pi_30[k]
                   + pb_y[k] * di_68[k];

        t_149[k] = pa_y[k] * pk_28[k];

        t_150[k] = pa_z[k] * pk_13[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pa_y, pa_z, pb_y, pb_z, pi_13, pi_32, \
                         pk_14, pk_29, di_69, di_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_11 * pi_13[k]
                   + pb_z[k] * di_69[k];

        t_152[k] = f_11 * pi_32[k]
                   + pb_y[k] * di_70[k];

        t_153[k] = pa_y[k] * pk_29[k];

        t_154[k] = pa_z[k] * pk_14[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_y, pb_y, pb_z, pi_15, pi_34, pi_35, \
                         pk_30, pk_31, di_71, di_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_11 * pi_15[k]
                   + pb_z[k] * di_71[k];

        t_156[k] = f_0 * pi_34[k]
                   + pa_y[k] * pk_30[k];

        t_157[k] = f_11 * pi_35[k]
                   + pb_y[k] * di_72[k];

        t_158[k] = pa_y[k] * pk_31[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_y, pa_z, pb_z, pi_17, pi_37, pi_38, \
                         pk_16, pk_32, pk_33, di_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * pk_16[k];

        t_160[k] = f_11 * pi_17[k]
                   + pb_z[k] * di_73[k];

        t_161[k] = f_14 * pi_37[k]
                   + pa_y[k] * pk_32[k];

        t_162[k] = f_0 * pi_38[k]
                   + pa_y[k] * pk_33[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, t_168, pa_y, pb_x, pb_y, pi_39, \
                         pk_34, di_74, di_75, di_76, di_77, di_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_11 * pi_39[k]
                   + pb_y[k] * di_74[k];

        t_164[k] = pa_y[k] * pk_34[k];

        t_165[k] = pb_x[k] * di_75[k];

        t_166[k] = pb_x[k] * di_76[k];

        t_167[k] = pb_x[k] * di_77[k];

        t_168[k] = pb_x[k] * di_78[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, pa_z, pb_x, pb_z, pi_23, pk_19, \
                         di_75, di_79, di_80, di_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = pb_x[k] * di_79[k];

        t_170[k] = pb_x[k] * di_80[k];

        t_171[k] = pb_x[k] * di_81[k];

        t_172[k] = pa_z[k] * pk_19[k];

        t_173[k] = f_11 * pi_23[k]
                   + pb_z[k] * di_75[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pa_y, pi_45, pi_46, pi_47, pi_48, pk_37, \
                         pk_38, pk_39, pk_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_12 * pi_45[k]
                   + pa_y[k] * pk_37[k];

        t_175[k] = f_13 * pi_46[k]
                   + pa_y[k] * pk_38[k];

        t_176[k] = f_14 * pi_47[k]
                   + pa_y[k] * pk_39[k];

        t_177[k] = f_0 * pi_48[k]
                   + pa_y[k] * pk_40[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pa_y, pb_x, pb_y, pb_z, pi_29, \
                         pi_49, pk_41, dh0_26, dh1_26, di_81, di_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_11 * pi_49[k]
                   + pb_y[k] * di_81[k];

        t_179[k] = pa_y[k] * pk_41[k];

        t_180[k] = f_1 * dh0_26[k]
                   - f_2 * dh1_26[k]
                   + pb_x[k] * di_82[k];

        t_181[k] = pb_y[k] * di_82[k];

        t_182[k] = f_0 * pi_29[k]
                   + pb_z[k] * di_82[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, pb_x, pb_y, dh0_27, dh0_28, dh0_29, \
                         dh1_27, dh1_28, dh1_29, di_83, di_84, di_85, \
                         di_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * dh0_27[k]
                   - f_10 * dh1_27[k]
                   + pb_x[k] * di_84[k];

        t_184[k] = pb_y[k] * di_83[k];

        t_185[k] = f_9 * dh0_28[k]
                   - f_10 * dh1_28[k]
                   + pb_x[k] * di_85[k];

        t_186[k] = f_7 * dh0_29[k]
                   - f_8 * dh1_29[k]
                   + pb_x[k] * di_86[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, pb_x, pb_y, pb_z, pi_31, dh0_30, dh0_31, \
                         dh1_30, dh1_31, di_84, di_85, di_87, di_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_0 * pi_31[k]
                   + pb_z[k] * di_84[k];

        t_188[k] = pb_y[k] * di_85[k];

        t_189[k] = f_7 * dh0_30[k]
                   - f_8 * dh1_30[k]
                   + pb_x[k] * di_87[k];

        t_190[k] = f_5 * dh0_31[k]
                   - f_6 * dh1_31[k]
                   + pb_x[k] * di_88[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pb_x, pb_y, pb_z, pi_33, dh0_32, dh0_33, \
                         dh1_32, dh1_33, di_86, di_87, di_89, di_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_0 * pi_33[k]
                   + pb_z[k] * di_86[k];

        t_192[k] = f_5 * dh0_32[k]
                   - f_6 * dh1_32[k]
                   + pb_x[k] * di_89[k];

        t_193[k] = pb_y[k] * di_87[k];

        t_194[k] = f_5 * dh0_33[k]
                   - f_6 * dh1_33[k]
                   + pb_x[k] * di_90[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pb_x, pb_z, pi_36, dh0_34, dh0_35, dh1_34, \
                         dh1_35, di_88, di_91, di_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_3 * dh0_34[k]
                   - f_4 * dh1_34[k]
                   + pb_x[k] * di_91[k];

        t_196[k] = f_0 * pi_36[k]
                   + pb_z[k] * di_88[k];

        t_197[k] = f_3 * dh0_35[k]
                   - f_4 * dh1_35[k]
                   + pb_x[k] * di_92[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, pb_x, pb_y, dh0_36, dh0_38, \
                         dh1_36, dh1_38, di_90, di_93, di_94, di_95, \
                         di_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_3 * dh0_36[k]
                   - f_4 * dh1_36[k]
                   + pb_x[k] * di_93[k];

        t_199[k] = pb_y[k] * di_90[k];

        t_200[k] = f_3 * dh0_38[k]
                   - f_4 * dh1_38[k]
                   + pb_x[k] * di_94[k];

        t_201[k] = pb_x[k] * di_95[k];

        t_202[k] = pb_x[k] * di_96[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, t_208, pb_x, pb_y, dh0_34, dh1_34, \
                         di_95, di_97, di_98, di_99, di_100, di_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = pb_x[k] * di_97[k];

        t_204[k] = pb_x[k] * di_98[k];

        t_205[k] = pb_x[k] * di_99[k];

        t_206[k] = pb_x[k] * di_100[k];

        t_207[k] = pb_x[k] * di_101[k];

        t_208[k] = f_1 * dh0_34[k]
                   - f_2 * dh1_34[k]
                   + pb_y[k] * di_95[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, pb_y, pb_z, pi_43, dh0_35, dh0_36, dh1_35, \
                         dh1_36, di_95, di_97, di_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_0 * pi_43[k]
                   + pb_z[k] * di_95[k];

        t_210[k] = f_9 * dh0_35[k]
                   - f_10 * dh1_35[k]
                   + pb_y[k] * di_97[k];

        t_211[k] = f_7 * dh0_36[k]
                   - f_8 * dh1_36[k]
                   + pb_y[k] * di_98[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pb_y, pb_z, pi_49, dh0_37, dh0_38, \
                         dh1_37, dh1_38, di_99, di_100, di_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_5 * dh0_37[k]
                   - f_6 * dh1_37[k]
                   + pb_y[k] * di_99[k];

        t_213[k] = f_3 * dh0_38[k]
                   - f_4 * dh1_38[k]
                   + pb_y[k] * di_100[k];

        t_214[k] = pb_y[k] * di_101[k];

        t_215[k] = f_0 * pi_49[k]
                   + f_1 * dh0_38[k]
                   - f_2 * dh1_38[k]
                   + pb_z[k] * di_101[k];
    }
}

auto
compute_prim_dk_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pi, const size_t pk,
                                     const size_t dh0, const size_t dh1, const size_t di,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 2.5 / p;
    const auto f_13 = 2.0 / p;
    const auto f_14 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_12 = buffer.data(pi + 12);
    const auto *pi_13 = buffer.data(pi + 13);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_16 = buffer.data(pi + 16);
    const auto *pi_17 = buffer.data(pi + 17);
    const auto *pi_18 = buffer.data(pi + 18);
    const auto *pi_19 = buffer.data(pi + 19);
    const auto *pi_20 = buffer.data(pi + 20);
    const auto *pi_21 = buffer.data(pi + 21);
    const auto *pi_22 = buffer.data(pi + 22);
    const auto *pi_23 = buffer.data(pi + 23);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);
    const auto *pk_3 = buffer.data(pk + 3);
    const auto *pk_4 = buffer.data(pk + 4);
    const auto *pk_5 = buffer.data(pk + 5);
    const auto *pk_6 = buffer.data(pk + 6);
    const auto *pk_7 = buffer.data(pk + 7);
    const auto *pk_8 = buffer.data(pk + 8);
    const auto *pk_12 = buffer.data(pk + 12);
    const auto *pk_13 = buffer.data(pk + 13);
    const auto *pk_14 = buffer.data(pk + 14);
    const auto *pk_15 = buffer.data(pk + 15);
    const auto *pk_16 = buffer.data(pk + 16);
    const auto *pk_17 = buffer.data(pk + 17);
    const auto *pk_18 = buffer.data(pk + 18);
    const auto *pk_19 = buffer.data(pk + 19);
    const auto *pk_20 = buffer.data(pk + 20);
    const auto *pk_21 = buffer.data(pk + 21);
    const auto *pk_24 = buffer.data(pk + 24);
    const auto *pk_25 = buffer.data(pk + 25);
    const auto *pk_26 = buffer.data(pk + 26);
    const auto *pk_27 = buffer.data(pk + 27);
    const auto *pk_28 = buffer.data(pk + 28);
    const auto *pk_30 = buffer.data(pk + 30);
    const auto *pk_31 = buffer.data(pk + 31);
    const auto *pk_32 = buffer.data(pk + 32);
    const auto *pk_33 = buffer.data(pk + 33);
    const auto *pk_34 = buffer.data(pk + 34);
    const auto *pk_35 = buffer.data(pk + 35);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);
    const auto *dh0_3 = buffer.data(dh0 + 3);
    const auto *dh0_4 = buffer.data(dh0 + 4);
    const auto *dh0_5 = buffer.data(dh0 + 5);
    const auto *dh0_6 = buffer.data(dh0 + 6);
    const auto *dh0_7 = buffer.data(dh0 + 7);
    const auto *dh0_8 = buffer.data(dh0 + 8);
    const auto *dh0_9 = buffer.data(dh0 + 9);
    const auto *dh0_10 = buffer.data(dh0 + 10);
    const auto *dh0_11 = buffer.data(dh0 + 11);
    const auto *dh0_12 = buffer.data(dh0 + 12);
    const auto *dh0_13 = buffer.data(dh0 + 13);
    const auto *dh0_14 = buffer.data(dh0 + 14);
    const auto *dh0_15 = buffer.data(dh0 + 15);
    const auto *dh0_16 = buffer.data(dh0 + 16);
    const auto *dh0_17 = buffer.data(dh0 + 17);
    const auto *dh0_18 = buffer.data(dh0 + 18);
    const auto *dh0_19 = buffer.data(dh0 + 19);
    const auto *dh0_20 = buffer.data(dh0 + 20);
    const auto *dh0_21 = buffer.data(dh0 + 21);
    const auto *dh0_22 = buffer.data(dh0 + 22);
    const auto *dh0_23 = buffer.data(dh0 + 23);
    const auto *dh0_24 = buffer.data(dh0 + 24);
    const auto *dh0_25 = buffer.data(dh0 + 25);
    const auto *dh0_26 = buffer.data(dh0 + 26);
    const auto *dh0_27 = buffer.data(dh0 + 27);
    const auto *dh0_28 = buffer.data(dh0 + 28);
    const auto *dh0_29 = buffer.data(dh0 + 29);
    const auto *dh0_30 = buffer.data(dh0 + 30);
    const auto *dh0_31 = buffer.data(dh0 + 31);
    const auto *dh0_32 = buffer.data(dh0 + 32);
    const auto *dh0_33 = buffer.data(dh0 + 33);
    const auto *dh0_34 = buffer.data(dh0 + 34);
    const auto *dh0_35 = buffer.data(dh0 + 35);
    const auto *dh0_36 = buffer.data(dh0 + 36);
    const auto *dh0_37 = buffer.data(dh0 + 37);
    const auto *dh0_38 = buffer.data(dh0 + 38);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);
    const auto *dh1_3 = buffer.data(dh1 + 3);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_5 = buffer.data(dh1 + 5);
    const auto *dh1_6 = buffer.data(dh1 + 6);
    const auto *dh1_7 = buffer.data(dh1 + 7);
    const auto *dh1_8 = buffer.data(dh1 + 8);
    const auto *dh1_9 = buffer.data(dh1 + 9);
    const auto *dh1_10 = buffer.data(dh1 + 10);
    const auto *dh1_11 = buffer.data(dh1 + 11);
    const auto *dh1_12 = buffer.data(dh1 + 12);
    const auto *dh1_13 = buffer.data(dh1 + 13);
    const auto *dh1_14 = buffer.data(dh1 + 14);
    const auto *dh1_15 = buffer.data(dh1 + 15);
    const auto *dh1_16 = buffer.data(dh1 + 16);
    const auto *dh1_17 = buffer.data(dh1 + 17);
    const auto *dh1_18 = buffer.data(dh1 + 18);
    const auto *dh1_19 = buffer.data(dh1 + 19);
    const auto *dh1_20 = buffer.data(dh1 + 20);
    const auto *dh1_21 = buffer.data(dh1 + 21);
    const auto *dh1_22 = buffer.data(dh1 + 22);
    const auto *dh1_23 = buffer.data(dh1 + 23);
    const auto *dh1_24 = buffer.data(dh1 + 24);
    const auto *dh1_25 = buffer.data(dh1 + 25);
    const auto *dh1_26 = buffer.data(dh1 + 26);
    const auto *dh1_27 = buffer.data(dh1 + 27);
    const auto *dh1_28 = buffer.data(dh1 + 28);
    const auto *dh1_29 = buffer.data(dh1 + 29);
    const auto *dh1_30 = buffer.data(dh1 + 30);
    const auto *dh1_31 = buffer.data(dh1 + 31);
    const auto *dh1_32 = buffer.data(dh1 + 32);
    const auto *dh1_33 = buffer.data(dh1 + 33);
    const auto *dh1_34 = buffer.data(dh1 + 34);
    const auto *dh1_35 = buffer.data(dh1 + 35);
    const auto *dh1_36 = buffer.data(dh1 + 36);
    const auto *dh1_37 = buffer.data(dh1 + 37);
    const auto *dh1_38 = buffer.data(dh1 + 38);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pi_0, dh0_0, dh1_0, di_0, \
                         di_1, di_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = pb_y[k] * di_0[k];

        t_2[k] = pb_z[k] * di_0[k];

        t_3[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_y[k] * di_1[k];

        t_4[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_z[k] * di_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, dh0_1, dh0_2, dh0_3, dh1_1, dh1_2, \
                         dh1_3, di_3, di_4, di_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * dh0_1[k]
                 - f_6 * dh1_1[k]
                 + pb_y[k] * di_3[k];

        t_6[k] = pb_y[k] * di_4[k];

        t_7[k] = f_5 * dh0_2[k]
                 - f_6 * dh1_2[k]
                 + pb_z[k] * di_4[k];

        t_8[k] = f_7 * dh0_3[k]
                 - f_8 * dh1_3[k]
                 + pb_y[k] * di_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, dh0_4, dh0_5, dh1_4, dh1_5, di_6, \
                         di_7, di_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * dh0_4[k]
                 - f_4 * dh1_4[k]
                 + pb_y[k] * di_6[k];

        t_10[k] = pb_y[k] * di_7[k];

        t_11[k] = f_7 * dh0_4[k]
                  - f_8 * dh1_4[k]
                  + pb_z[k] * di_7[k];

        t_12[k] = f_9 * dh0_5[k]
                  - f_10 * dh1_5[k]
                  + pb_y[k] * di_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pb_y, pb_z, dh0_6, dh0_7, dh1_6, dh1_7, di_9, \
                         di_10, di_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * dh0_6[k]
                  - f_6 * dh1_6[k]
                  + pb_y[k] * di_9[k];

        t_14[k] = f_3 * dh0_7[k]
                  - f_4 * dh1_7[k]
                  + pb_y[k] * di_10[k];

        t_15[k] = pb_y[k] * di_11[k];

        t_16[k] = f_9 * dh0_7[k]
                  - f_10 * dh1_7[k]
                  + pb_z[k] * di_11[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_x, pb_y, pi_1, pi_2, dh0_8, dh0_9, dh1_8, \
                         dh1_9, di_12, di_13, di_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * pi_1[k]
                  + pb_x[k] * di_12[k];

        t_18[k] = f_0 * pi_2[k]
                  + pb_x[k] * di_17[k];

        t_19[k] = f_1 * dh0_8[k]
                  - f_2 * dh1_8[k]
                  + pb_y[k] * di_12[k];

        t_20[k] = f_9 * dh0_9[k]
                  - f_10 * dh1_9[k]
                  + pb_y[k] * di_13[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, dh0_10, dh0_11, dh0_12, dh1_10, dh1_11, \
                         dh1_12, di_14, di_15, di_16, di_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * dh0_10[k]
                  - f_8 * dh1_10[k]
                  + pb_y[k] * di_14[k];

        t_22[k] = f_5 * dh0_11[k]
                  - f_6 * dh1_11[k]
                  + pb_y[k] * di_15[k];

        t_23[k] = f_3 * dh0_12[k]
                  - f_4 * dh1_12[k]
                  + pb_y[k] * di_16[k];

        t_24[k] = pb_y[k] * di_17[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pa_y, pb_y, pb_z, pi_0, pi_4, pk_0, \
                         pk_12, dh0_12, dh1_12, di_17, di_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * dh0_12[k]
                  - f_2 * dh1_12[k]
                  + pb_z[k] * di_17[k];

        t_26[k] = pa_y[k] * pk_0[k];

        t_27[k] = f_11 * pi_0[k]
                  + pb_y[k] * di_18[k];

        t_28[k] = f_12 * pi_4[k]
                  + pa_x[k] * pk_12[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pa_x, pa_y, pi_5, pi_6, pk_2, pk_4, \
                         pk_6, pk_13, pk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * pk_2[k];

        t_30[k] = f_13 * pi_5[k]
                  + pa_x[k] * pk_13[k];

        t_31[k] = pa_y[k] * pk_4[k];

        t_32[k] = f_14 * pi_6[k]
                  + pa_x[k] * pk_14[k];

        t_33[k] = pa_y[k] * pk_6[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_x, pa_y, pb_x, pi_7, pi_8, pk_8, \
                         pk_15, pk_16, pk_17, di_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_0 * pi_7[k]
                  + pa_x[k] * pk_15[k];

        t_35[k] = pa_y[k] * pk_8[k];

        t_36[k] = f_11 * pi_8[k]
                  + pb_x[k] * di_19[k];

        t_37[k] = pa_x[k] * pk_16[k];

        t_38[k] = pa_x[k] * pk_17[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, t_44, pa_x, pa_z, pb_z, pi_0, pk_0, \
                         pk_18, pk_19, pk_20, pk_21, di_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pa_x[k] * pk_18[k];

        t_40[k] = pa_x[k] * pk_19[k];

        t_41[k] = pa_x[k] * pk_20[k];

        t_42[k] = pa_x[k] * pk_21[k];

        t_43[k] = pa_z[k] * pk_0[k];

        t_44[k] = f_11 * pi_0[k]
                  + pb_z[k] * di_20[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_x, pa_z, pi_14, pi_15, pk_1, pk_3, \
                         pk_5, pk_25, pk_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_z[k] * pk_1[k];

        t_46[k] = f_12 * pi_14[k]
                  + pa_x[k] * pk_25[k];

        t_47[k] = pa_z[k] * pk_3[k];

        t_48[k] = f_13 * pi_15[k]
                  + pa_x[k] * pk_26[k];

        t_49[k] = pa_z[k] * pk_5[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_x, pa_z, pb_x, pi_16, pi_17, pi_23, \
                         pk_7, pk_27, pk_28, pk_30, di_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_14 * pi_16[k]
                  + pa_x[k] * pk_27[k];

        t_51[k] = pa_z[k] * pk_7[k];

        t_52[k] = f_0 * pi_17[k]
                  + pa_x[k] * pk_28[k];

        t_53[k] = f_11 * pi_23[k]
                  + pb_x[k] * di_21[k];

        t_54[k] = pa_x[k] * pk_30[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, pa_x, pb_x, pk_31, pk_32, pk_33, \
                         pk_34, pk_35, dh0_13, dh1_13, di_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_x[k] * pk_31[k];

        t_56[k] = pa_x[k] * pk_32[k];

        t_57[k] = pa_x[k] * pk_33[k];

        t_58[k] = pa_x[k] * pk_34[k];

        t_59[k] = pa_x[k] * pk_35[k];

        t_60[k] = f_1 * dh0_13[k]
                  - f_2 * dh1_13[k]
                  + pb_x[k] * di_22[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pb_x, pb_y, pi_3, dh0_14, dh0_15, dh1_14, dh1_15, \
                         di_22, di_23, di_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_0 * pi_3[k]
                  + pb_y[k] * di_22[k];

        t_62[k] = f_9 * dh0_14[k]
                  - f_10 * dh1_14[k]
                  + pb_x[k] * di_23[k];

        t_63[k] = f_9 * dh0_15[k]
                  - f_10 * dh1_15[k]
                  + pb_x[k] * di_24[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_x, dh0_16, dh0_17, dh0_18, dh1_16, dh1_17, \
                         dh1_18, di_25, di_26, di_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_7 * dh0_16[k]
                  - f_8 * dh1_16[k]
                  + pb_x[k] * di_25[k];

        t_65[k] = f_7 * dh0_17[k]
                  - f_8 * dh1_17[k]
                  + pb_x[k] * di_26[k];

        t_66[k] = f_5 * dh0_18[k]
                  - f_6 * dh1_18[k]
                  + pb_x[k] * di_27[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pb_x, dh0_19, dh0_20, dh0_21, dh1_19, dh1_20, \
                         dh1_21, di_28, di_29, di_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_5 * dh0_19[k]
                  - f_6 * dh1_19[k]
                  + pb_x[k] * di_28[k];

        t_68[k] = f_5 * dh0_20[k]
                  - f_6 * dh1_20[k]
                  + pb_x[k] * di_29[k];

        t_69[k] = f_3 * dh0_21[k]
                  - f_4 * dh1_21[k]
                  + pb_x[k] * di_30[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_x, dh0_23, dh0_24, dh0_25, dh1_23, dh1_24, \
                         dh1_25, di_31, di_32, di_33, di_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * dh0_23[k]
                  - f_4 * dh1_23[k]
                  + pb_x[k] * di_31[k];

        t_71[k] = f_3 * dh0_24[k]
                  - f_4 * dh1_24[k]
                  + pb_x[k] * di_32[k];

        t_72[k] = f_3 * dh0_25[k]
                  - f_4 * dh1_25[k]
                  + pb_x[k] * di_33[k];

        t_73[k] = pb_x[k] * di_34[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, pb_x, pb_y, pi_8, dh0_21, dh1_21, \
                         di_34, di_36, di_37, di_38, di_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pb_x[k] * di_36[k];

        t_75[k] = pb_x[k] * di_37[k];

        t_76[k] = pb_x[k] * di_38[k];

        t_77[k] = pb_x[k] * di_39[k];

        t_78[k] = f_0 * pi_8[k]
                  + f_1 * dh0_21[k]
                  - f_2 * dh1_21[k]
                  + pb_y[k] * di_34[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_z, dh0_21, dh0_22, dh0_23, dh1_21, dh1_22, \
                         dh1_23, di_34, di_35, di_36, di_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_z[k] * di_34[k];

        t_80[k] = f_3 * dh0_21[k]
                  - f_4 * dh1_21[k]
                  + pb_z[k] * di_35[k];

        t_81[k] = f_5 * dh0_22[k]
                  - f_6 * dh1_22[k]
                  + pb_z[k] * di_36[k];

        t_82[k] = f_7 * dh0_23[k]
                  - f_8 * dh1_23[k]
                  + pb_z[k] * di_37[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pb_y, pb_z, pi_12, pk_24, dh0_24, \
                         dh0_25, dh1_24, dh1_25, di_38, di_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_9 * dh0_24[k]
                  - f_10 * dh1_24[k]
                  + pb_z[k] * di_38[k];

        t_84[k] = f_0 * pi_12[k]
                  + pb_y[k] * di_39[k];

        t_85[k] = f_1 * dh0_25[k]
                  - f_2 * dh1_25[k]
                  + pb_z[k] * di_39[k];

        t_86[k] = pa_y[k] * pk_24[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, t_93, pa_y, pa_z, pk_12, pk_13, \
                         pk_14, pk_15, pk_25, pk_26, pk_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * pk_12[k];

        t_88[k] = pa_y[k] * pk_25[k];

        t_89[k] = pa_z[k] * pk_13[k];

        t_90[k] = pa_y[k] * pk_26[k];

        t_91[k] = pa_z[k] * pk_14[k];

        t_92[k] = pa_y[k] * pk_27[k];

        t_93[k] = pa_z[k] * pk_15[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pa_y, pa_z, pb_z, pi_8, pi_19, pi_20, \
                         pk_16, pk_28, pk_31, pk_32, di_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = pa_y[k] * pk_28[k];

        t_95[k] = pa_z[k] * pk_16[k];

        t_96[k] = f_11 * pi_8[k]
                  + pb_z[k] * di_40[k];

        t_97[k] = f_12 * pi_19[k]
                  + pa_y[k] * pk_31[k];

        t_98[k] = f_13 * pi_20[k]
                  + pa_y[k] * pk_32[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_y, pb_y, pi_21, pi_22, pi_23, pk_33, \
                         pk_34, pk_35, di_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_14 * pi_21[k]
                  + pa_y[k] * pk_33[k];

        t_100[k] = f_0 * pi_22[k]
                   + pa_y[k] * pk_34[k];

        t_101[k] = f_11 * pi_23[k]
                   + pb_y[k] * di_41[k];

        t_102[k] = pa_y[k] * pk_35[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pb_x, pb_z, pi_13, dh0_26, dh0_27, \
                         dh0_28, dh1_26, dh1_27, dh1_28, di_42, di_43, \
                         di_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_1 * dh0_26[k]
                   - f_2 * dh1_26[k]
                   + pb_x[k] * di_42[k];

        t_104[k] = f_0 * pi_13[k]
                   + pb_z[k] * di_42[k];

        t_105[k] = f_9 * dh0_27[k]
                   - f_10 * dh1_27[k]
                   + pb_x[k] * di_43[k];

        t_106[k] = f_9 * dh0_28[k]
                   - f_10 * dh1_28[k]
                   + pb_x[k] * di_44[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pb_x, dh0_29, dh0_30, dh0_31, dh1_29, dh1_30, \
                         dh1_31, di_45, di_46, di_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_7 * dh0_29[k]
                   - f_8 * dh1_29[k]
                   + pb_x[k] * di_45[k];

        t_108[k] = f_7 * dh0_30[k]
                   - f_8 * dh1_30[k]
                   + pb_x[k] * di_46[k];

        t_109[k] = f_5 * dh0_31[k]
                   - f_6 * dh1_31[k]
                   + pb_x[k] * di_47[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pb_x, dh0_32, dh0_33, dh0_34, dh1_32, dh1_33, \
                         dh1_34, di_48, di_49, di_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_5 * dh0_32[k]
                   - f_6 * dh1_32[k]
                   + pb_x[k] * di_48[k];

        t_111[k] = f_5 * dh0_33[k]
                   - f_6 * dh1_33[k]
                   + pb_x[k] * di_49[k];

        t_112[k] = f_3 * dh0_34[k]
                   - f_4 * dh1_34[k]
                   + pb_x[k] * di_50[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pb_x, dh0_35, dh0_36, dh0_38, dh1_35, \
                         dh1_36, dh1_38, di_51, di_52, di_53, di_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_3 * dh0_35[k]
                   - f_4 * dh1_35[k]
                   + pb_x[k] * di_51[k];

        t_114[k] = f_3 * dh0_36[k]
                   - f_4 * dh1_36[k]
                   + pb_x[k] * di_52[k];

        t_115[k] = f_3 * dh0_38[k]
                   - f_4 * dh1_38[k]
                   + pb_x[k] * di_53[k];

        t_116[k] = pb_x[k] * di_54[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, pb_x, pb_y, dh0_34, dh1_34, di_54, \
                         di_55, di_56, di_57, di_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pb_x[k] * di_55[k];

        t_118[k] = pb_x[k] * di_56[k];

        t_119[k] = pb_x[k] * di_57[k];

        t_120[k] = pb_x[k] * di_59[k];

        t_121[k] = f_1 * dh0_34[k]
                   - f_2 * dh1_34[k]
                   + pb_y[k] * di_54[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pb_y, pb_z, pi_18, dh0_35, dh0_36, dh1_35, \
                         dh1_36, di_54, di_55, di_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_0 * pi_18[k]
                   + pb_z[k] * di_54[k];

        t_123[k] = f_9 * dh0_35[k]
                   - f_10 * dh1_35[k]
                   + pb_y[k] * di_55[k];

        t_124[k] = f_7 * dh0_36[k]
                   - f_8 * dh1_36[k]
                   + pb_y[k] * di_56[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pb_z, pi_23, dh0_37, dh0_38, \
                         dh1_37, dh1_38, di_57, di_58, di_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_5 * dh0_37[k]
                   - f_6 * dh1_37[k]
                   + pb_y[k] * di_57[k];

        t_126[k] = f_3 * dh0_38[k]
                   - f_4 * dh1_38[k]
                   + pb_y[k] * di_58[k];

        t_127[k] = pb_y[k] * di_59[k];

        t_128[k] = f_0 * pi_23[k]
                   + f_1 * dh0_38[k]
                   - f_2 * dh1_38[k]
                   + pb_z[k] * di_59[k];
    }
}

auto
compute_prim_dk_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pi, const size_t pk,
                                     const size_t dh0, const size_t dh1, const size_t di,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 2.5 / p;
    const auto f_13 = 2.0 / p;
    const auto f_14 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_11 = buffer.data(pi + 11);
    const auto *pi_12 = buffer.data(pi + 12);
    const auto *pi_13 = buffer.data(pi + 13);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_19 = buffer.data(pi + 19);
    const auto *pi_20 = buffer.data(pi + 20);
    const auto *pi_22 = buffer.data(pi + 22);
    const auto *pi_23 = buffer.data(pi + 23);
    const auto *pi_24 = buffer.data(pi + 24);
    const auto *pi_25 = buffer.data(pi + 25);
    const auto *pi_26 = buffer.data(pi + 26);
    const auto *pi_27 = buffer.data(pi + 27);
    const auto *pi_28 = buffer.data(pi + 28);
    const auto *pi_29 = buffer.data(pi + 29);
    const auto *pi_31 = buffer.data(pi + 31);
    const auto *pi_32 = buffer.data(pi + 32);
    const auto *pi_33 = buffer.data(pi + 33);
    const auto *pi_34 = buffer.data(pi + 34);
    const auto *pi_35 = buffer.data(pi + 35);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);
    const auto *pk_3 = buffer.data(pk + 3);
    const auto *pk_4 = buffer.data(pk + 4);
    const auto *pk_5 = buffer.data(pk + 5);
    const auto *pk_6 = buffer.data(pk + 6);
    const auto *pk_7 = buffer.data(pk + 7);
    const auto *pk_8 = buffer.data(pk + 8);
    const auto *pk_9 = buffer.data(pk + 9);
    const auto *pk_10 = buffer.data(pk + 10);
    const auto *pk_11 = buffer.data(pk + 11);
    const auto *pk_12 = buffer.data(pk + 12);
    const auto *pk_13 = buffer.data(pk + 13);
    const auto *pk_14 = buffer.data(pk + 14);
    const auto *pk_15 = buffer.data(pk + 15);
    const auto *pk_16 = buffer.data(pk + 16);
    const auto *pk_17 = buffer.data(pk + 17);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);
    const auto *dh0_3 = buffer.data(dh0 + 3);
    const auto *dh0_4 = buffer.data(dh0 + 4);
    const auto *dh0_5 = buffer.data(dh0 + 5);
    const auto *dh0_6 = buffer.data(dh0 + 6);
    const auto *dh0_7 = buffer.data(dh0 + 7);
    const auto *dh0_8 = buffer.data(dh0 + 8);
    const auto *dh0_9 = buffer.data(dh0 + 9);
    const auto *dh0_10 = buffer.data(dh0 + 10);
    const auto *dh0_11 = buffer.data(dh0 + 11);
    const auto *dh0_12 = buffer.data(dh0 + 12);
    const auto *dh0_13 = buffer.data(dh0 + 13);
    const auto *dh0_14 = buffer.data(dh0 + 14);
    const auto *dh0_15 = buffer.data(dh0 + 15);
    const auto *dh0_16 = buffer.data(dh0 + 16);
    const auto *dh0_17 = buffer.data(dh0 + 17);
    const auto *dh0_18 = buffer.data(dh0 + 18);
    const auto *dh0_19 = buffer.data(dh0 + 19);
    const auto *dh0_20 = buffer.data(dh0 + 20);
    const auto *dh0_21 = buffer.data(dh0 + 21);
    const auto *dh0_22 = buffer.data(dh0 + 22);
    const auto *dh0_23 = buffer.data(dh0 + 23);
    const auto *dh0_24 = buffer.data(dh0 + 24);
    const auto *dh0_25 = buffer.data(dh0 + 25);
    const auto *dh0_26 = buffer.data(dh0 + 26);
    const auto *dh0_27 = buffer.data(dh0 + 27);
    const auto *dh0_28 = buffer.data(dh0 + 28);
    const auto *dh0_29 = buffer.data(dh0 + 29);
    const auto *dh0_30 = buffer.data(dh0 + 30);
    const auto *dh0_31 = buffer.data(dh0 + 31);
    const auto *dh0_32 = buffer.data(dh0 + 32);
    const auto *dh0_33 = buffer.data(dh0 + 33);
    const auto *dh0_34 = buffer.data(dh0 + 34);
    const auto *dh0_35 = buffer.data(dh0 + 35);
    const auto *dh0_36 = buffer.data(dh0 + 36);
    const auto *dh0_37 = buffer.data(dh0 + 37);
    const auto *dh0_38 = buffer.data(dh0 + 38);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);
    const auto *dh1_3 = buffer.data(dh1 + 3);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_5 = buffer.data(dh1 + 5);
    const auto *dh1_6 = buffer.data(dh1 + 6);
    const auto *dh1_7 = buffer.data(dh1 + 7);
    const auto *dh1_8 = buffer.data(dh1 + 8);
    const auto *dh1_10 = buffer.data(dh1 + 10);
    const auto *dh1_11 = buffer.data(dh1 + 11);
    const auto *dh1_12 = buffer.data(dh1 + 12);
    const auto *dh1_13 = buffer.data(dh1 + 13);
    const auto *dh1_18 = buffer.data(dh1 + 18);
    const auto *dh1_20 = buffer.data(dh1 + 20);
    const auto *dh1_21 = buffer.data(dh1 + 21);
    const auto *dh1_22 = buffer.data(dh1 + 22);
    const auto *dh1_23 = buffer.data(dh1 + 23);
    const auto *dh1_24 = buffer.data(dh1 + 24);
    const auto *dh1_25 = buffer.data(dh1 + 25);
    const auto *dh1_26 = buffer.data(dh1 + 26);
    const auto *dh1_27 = buffer.data(dh1 + 27);
    const auto *dh1_28 = buffer.data(dh1 + 28);
    const auto *dh1_29 = buffer.data(dh1 + 29);
    const auto *dh1_30 = buffer.data(dh1 + 30);
    const auto *dh1_31 = buffer.data(dh1 + 31);
    const auto *dh1_36 = buffer.data(dh1 + 36);
    const auto *dh1_38 = buffer.data(dh1 + 38);
    const auto *dh1_39 = buffer.data(dh1 + 39);
    const auto *dh1_40 = buffer.data(dh1 + 40);
    const auto *dh1_41 = buffer.data(dh1 + 41);
    const auto *dh1_42 = buffer.data(dh1 + 42);
    const auto *dh1_43 = buffer.data(dh1 + 43);
    const auto *dh1_44 = buffer.data(dh1 + 44);
    const auto *dh1_45 = buffer.data(dh1 + 45);
    const auto *dh1_46 = buffer.data(dh1 + 46);
    const auto *dh1_47 = buffer.data(dh1 + 47);
    const auto *dh1_48 = buffer.data(dh1 + 48);
    const auto *dh1_49 = buffer.data(dh1 + 49);

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
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_23 = buffer.data(di + 23);
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
    const auto *di_47 = buffer.data(di + 47);
    const auto *di_48 = buffer.data(di + 48);
    const auto *di_54 = buffer.data(di + 54);
    const auto *di_55 = buffer.data(di + 55);
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
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_71 = buffer.data(di + 71);
    const auto *di_72 = buffer.data(di + 72);
    const auto *di_73 = buffer.data(di + 73);
    const auto *di_74 = buffer.data(di + 74);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pi_0, dh0_0, dh0_1, dh1_0, \
                         dh1_1, di_0, di_1, di_2, di_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_y[k] * di_1[k];

        t_2[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_z[k] * di_2[k];

        t_3[k] = f_5 * dh0_1[k]
                 - f_6 * dh1_1[k]
                 + pb_y[k] * di_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, dh0_2, dh0_3, dh0_4, dh1_2, dh1_3, \
                         dh1_4, di_4, di_5, di_6, di_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * dh0_2[k]
                 - f_6 * dh1_2[k]
                 + pb_z[k] * di_4[k];

        t_5[k] = f_7 * dh0_3[k]
                 - f_8 * dh1_3[k]
                 + pb_y[k] * di_5[k];

        t_6[k] = f_3 * dh0_4[k]
                 - f_4 * dh1_4[k]
                 + pb_y[k] * di_6[k];

        t_7[k] = f_7 * dh0_4[k]
                 - f_8 * dh1_4[k]
                 + pb_z[k] * di_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, dh0_5, dh0_6, dh0_7, dh1_5, dh1_6, \
                         dh1_7, di_8, di_9, di_10, di_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_9 * dh0_5[k]
                 - f_10 * dh1_5[k]
                 + pb_y[k] * di_8[k];

        t_9[k] = f_5 * dh0_6[k]
                 - f_6 * dh1_6[k]
                 + pb_y[k] * di_9[k];

        t_10[k] = f_3 * dh0_7[k]
                  - f_4 * dh1_7[k]
                  + pb_y[k] * di_10[k];

        t_11[k] = f_9 * dh0_7[k]
                  - f_10 * dh1_7[k]
                  + pb_z[k] * di_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_x, pb_y, pi_7, pi_8, dh0_8, dh0_9, dh1_8, \
                         dh1_10, di_12, di_13, di_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * pi_7[k]
                  + pb_x[k] * di_12[k];

        t_13[k] = f_0 * pi_8[k]
                  + pb_x[k] * di_17[k];

        t_14[k] = f_1 * dh0_8[k]
                  - f_2 * dh1_8[k]
                  + pb_y[k] * di_12[k];

        t_15[k] = f_9 * dh0_9[k]
                  - f_10 * dh1_10[k]
                  + pb_y[k] * di_13[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_y, pb_z, dh0_10, dh0_11, dh0_12, dh1_11, \
                         dh1_12, dh1_13, di_14, di_15, di_16, di_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * dh0_10[k]
                  - f_8 * dh1_11[k]
                  + pb_y[k] * di_14[k];

        t_17[k] = f_5 * dh0_11[k]
                  - f_6 * dh1_12[k]
                  + pb_y[k] * di_15[k];

        t_18[k] = f_3 * dh0_12[k]
                  - f_4 * dh1_13[k]
                  + pb_y[k] * di_16[k];

        t_19[k] = f_1 * dh0_12[k]
                  - f_2 * dh1_13[k]
                  + pb_z[k] * di_17[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_y, pi_0, pi_10, pi_11, pi_12, pk_1, \
                         pk_2, pk_3, di_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_11 * pi_0[k]
                  + pb_y[k] * di_18[k];

        t_21[k] = f_12 * pi_10[k]
                  + pa_x[k] * pk_1[k];

        t_22[k] = f_13 * pi_11[k]
                  + pa_x[k] * pk_2[k];

        t_23[k] = f_14 * pi_12[k]
                  + pa_x[k] * pk_3[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_x, pa_z, pb_x, pb_z, pi_0, pi_13, pi_14, \
                         pk_0, pk_4, di_22, di_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * pi_13[k]
                  + pa_x[k] * pk_4[k];

        t_25[k] = f_11 * pi_14[k]
                  + pb_x[k] * di_22[k];

        t_26[k] = pa_z[k] * pk_0[k];

        t_27[k] = f_11 * pi_0[k]
                  + pb_z[k] * di_23[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pi_22, pi_24, pi_27, pi_28, pk_6, pk_7, \
                         pk_9, pk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_12 * pi_22[k]
                  + pa_x[k] * pk_6[k];

        t_29[k] = f_13 * pi_24[k]
                  + pa_x[k] * pk_7[k];

        t_30[k] = f_14 * pi_27[k]
                  + pa_x[k] * pk_9[k];

        t_31[k] = f_0 * pi_28[k]
                  + pa_x[k] * pk_12[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_x, pb_y, pi_9, pi_35, dh0_13, dh0_14, \
                         dh1_18, dh1_20, di_28, di_29, di_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_11 * pi_35[k]
                  + pb_x[k] * di_28[k];

        t_33[k] = f_1 * dh0_13[k]
                  - f_2 * dh1_18[k]
                  + pb_x[k] * di_29[k];

        t_34[k] = f_0 * pi_9[k]
                  + pb_y[k] * di_29[k];

        t_35[k] = f_9 * dh0_14[k]
                  - f_10 * dh1_20[k]
                  + pb_x[k] * di_30[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_x, dh0_15, dh0_16, dh0_17, dh1_21, dh1_22, \
                         dh1_23, di_31, di_32, di_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * dh0_15[k]
                  - f_10 * dh1_21[k]
                  + pb_x[k] * di_31[k];

        t_37[k] = f_7 * dh0_16[k]
                  - f_8 * dh1_22[k]
                  + pb_x[k] * di_32[k];

        t_38[k] = f_7 * dh0_17[k]
                  - f_8 * dh1_23[k]
                  + pb_x[k] * di_33[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_x, dh0_18, dh0_19, dh0_20, dh1_24, dh1_25, \
                         dh1_26, di_34, di_35, di_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_5 * dh0_18[k]
                  - f_6 * dh1_24[k]
                  + pb_x[k] * di_34[k];

        t_40[k] = f_5 * dh0_19[k]
                  - f_6 * dh1_25[k]
                  + pb_x[k] * di_35[k];

        t_41[k] = f_5 * dh0_20[k]
                  - f_6 * dh1_26[k]
                  + pb_x[k] * di_36[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, dh0_21, dh0_23, dh0_24, dh1_27, dh1_29, \
                         dh1_30, di_37, di_38, di_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * dh0_21[k]
                  - f_4 * dh1_27[k]
                  + pb_x[k] * di_37[k];

        t_43[k] = f_3 * dh0_23[k]
                  - f_4 * dh1_29[k]
                  + pb_x[k] * di_38[k];

        t_44[k] = f_3 * dh0_24[k]
                  - f_4 * dh1_30[k]
                  + pb_x[k] * di_39[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pb_x, pb_y, pb_z, pi_14, dh0_21, dh0_25, dh1_27, \
                         dh1_31, di_40, di_41, di_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_3 * dh0_25[k]
                  - f_4 * dh1_31[k]
                  + pb_x[k] * di_40[k];

        t_46[k] = f_0 * pi_14[k]
                  + f_1 * dh0_21[k]
                  - f_2 * dh1_27[k]
                  + pb_y[k] * di_41[k];

        t_47[k] = f_3 * dh0_21[k]
                  - f_4 * dh1_27[k]
                  + pb_z[k] * di_42[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_z, dh0_22, dh0_23, dh0_24, dh1_28, dh1_29, \
                         dh1_30, di_43, di_44, di_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * dh0_22[k]
                  - f_6 * dh1_28[k]
                  + pb_z[k] * di_43[k];

        t_49[k] = f_7 * dh0_23[k]
                  - f_8 * dh1_29[k]
                  + pb_z[k] * di_44[k];

        t_50[k] = f_9 * dh0_24[k]
                  - f_10 * dh1_30[k]
                  + pb_z[k] * di_45[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pb_y, pb_z, pi_19, pi_23, pi_25, pk_8, \
                         pk_10, dh0_25, dh1_31, di_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * pi_19[k]
                  + pb_y[k] * di_47[k];

        t_52[k] = f_1 * dh0_25[k]
                  - f_2 * dh1_31[k]
                  + pb_z[k] * di_47[k];

        t_53[k] = f_0 * pi_23[k]
                  + pa_y[k] * pk_8[k];

        t_54[k] = f_14 * pi_25[k]
                  + pa_y[k] * pk_10[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pa_z, pb_z, pi_14, pi_26, pi_31, pk_5, \
                         pk_11, pk_13, di_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_0 * pi_26[k]
                  + pa_y[k] * pk_11[k];

        t_56[k] = pa_z[k] * pk_5[k];

        t_57[k] = f_11 * pi_14[k]
                  + pb_z[k] * di_48[k];

        t_58[k] = f_12 * pi_31[k]
                  + pa_y[k] * pk_13[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, t_63, pa_y, pb_y, pi_32, pi_33, pi_34, pi_35, \
                         pk_14, pk_15, pk_16, pk_17, di_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_13 * pi_32[k]
                  + pa_y[k] * pk_14[k];

        t_60[k] = f_14 * pi_33[k]
                  + pa_y[k] * pk_15[k];

        t_61[k] = f_0 * pi_34[k]
                  + pa_y[k] * pk_16[k];

        t_62[k] = f_11 * pi_35[k]
                  + pb_y[k] * di_54[k];

        t_63[k] = pa_y[k] * pk_17[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pb_x, pb_z, pi_20, dh0_26, dh0_27, dh0_28, \
                         dh1_36, dh1_38, dh1_39, di_55, di_57, di_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_1 * dh0_26[k]
                  - f_2 * dh1_36[k]
                  + pb_x[k] * di_55[k];

        t_65[k] = f_0 * pi_20[k]
                  + pb_z[k] * di_55[k];

        t_66[k] = f_9 * dh0_27[k]
                  - f_10 * dh1_38[k]
                  + pb_x[k] * di_57[k];

        t_67[k] = f_9 * dh0_28[k]
                  - f_10 * dh1_39[k]
                  + pb_x[k] * di_58[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pb_x, dh0_29, dh0_30, dh0_31, dh1_40, dh1_41, \
                         dh1_42, di_59, di_60, di_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_7 * dh0_29[k]
                  - f_8 * dh1_40[k]
                  + pb_x[k] * di_59[k];

        t_69[k] = f_7 * dh0_30[k]
                  - f_8 * dh1_41[k]
                  + pb_x[k] * di_60[k];

        t_70[k] = f_5 * dh0_31[k]
                  - f_6 * dh1_42[k]
                  + pb_x[k] * di_61[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pb_x, dh0_32, dh0_33, dh0_34, dh1_43, dh1_44, \
                         dh1_45, di_62, di_63, di_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * dh0_32[k]
                  - f_6 * dh1_43[k]
                  + pb_x[k] * di_62[k];

        t_72[k] = f_5 * dh0_33[k]
                  - f_6 * dh1_44[k]
                  + pb_x[k] * di_63[k];

        t_73[k] = f_3 * dh0_34[k]
                  - f_4 * dh1_45[k]
                  + pb_x[k] * di_64[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pb_x, dh0_35, dh0_36, dh0_38, dh1_46, dh1_47, \
                         dh1_49, di_65, di_66, di_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_3 * dh0_35[k]
                  - f_4 * dh1_46[k]
                  + pb_x[k] * di_65[k];

        t_75[k] = f_3 * dh0_36[k]
                  - f_4 * dh1_47[k]
                  + pb_x[k] * di_66[k];

        t_76[k] = f_3 * dh0_38[k]
                  - f_4 * dh1_49[k]
                  + pb_x[k] * di_67[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pb_y, pb_z, pi_29, dh0_34, dh0_35, dh0_36, \
                         dh1_45, dh1_46, dh1_47, di_68, di_70, di_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_1 * dh0_34[k]
                  - f_2 * dh1_45[k]
                  + pb_y[k] * di_68[k];

        t_78[k] = f_0 * pi_29[k]
                  + pb_z[k] * di_68[k];

        t_79[k] = f_9 * dh0_35[k]
                  - f_10 * dh1_46[k]
                  + pb_y[k] * di_70[k];

        t_80[k] = f_7 * dh0_36[k]
                  - f_8 * dh1_47[k]
                  + pb_y[k] * di_71[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_y, pb_z, pi_35, dh0_37, dh0_38, dh1_48, dh1_49, \
                         di_72, di_73, di_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * dh0_37[k]
                  - f_6 * dh1_48[k]
                  + pb_y[k] * di_72[k];

        t_82[k] = f_3 * dh0_38[k]
                  - f_4 * dh1_49[k]
                  + pb_y[k] * di_73[k];

        t_83[k] = f_0 * pi_35[k]
                  + f_1 * dh0_38[k]
                  - f_2 * dh1_49[k]
                  + pb_z[k] * di_74[k];
    }
}

auto
compute_prim_dk_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pi, const size_t pk,
                                     const size_t dh0, const size_t dh1, const size_t di,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_11 = buffer.data(pi + 11);
    const auto *pi_12 = buffer.data(pi + 12);
    const auto *pi_13 = buffer.data(pi + 13);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_16 = buffer.data(pi + 16);
    const auto *pi_17 = buffer.data(pi + 17);
    const auto *pi_18 = buffer.data(pi + 18);
    const auto *pi_19 = buffer.data(pi + 19);
    const auto *pi_20 = buffer.data(pi + 20);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);
    const auto *pk_3 = buffer.data(pk + 3);
    const auto *pk_4 = buffer.data(pk + 4);
    const auto *pk_5 = buffer.data(pk + 5);
    const auto *pk_6 = buffer.data(pk + 6);
    const auto *pk_7 = buffer.data(pk + 7);
    const auto *pk_8 = buffer.data(pk + 8);
    const auto *pk_12 = buffer.data(pk + 12);
    const auto *pk_13 = buffer.data(pk + 13);
    const auto *pk_14 = buffer.data(pk + 14);
    const auto *pk_15 = buffer.data(pk + 15);
    const auto *pk_17 = buffer.data(pk + 17);
    const auto *pk_18 = buffer.data(pk + 18);
    const auto *pk_19 = buffer.data(pk + 19);
    const auto *pk_20 = buffer.data(pk + 20);
    const auto *pk_21 = buffer.data(pk + 21);
    const auto *pk_22 = buffer.data(pk + 22);
    const auto *pk_25 = buffer.data(pk + 25);
    const auto *pk_26 = buffer.data(pk + 26);
    const auto *pk_28 = buffer.data(pk + 28);
    const auto *pk_29 = buffer.data(pk + 29);
    const auto *pk_31 = buffer.data(pk + 31);
    const auto *pk_32 = buffer.data(pk + 32);
    const auto *pk_33 = buffer.data(pk + 33);
    const auto *pk_35 = buffer.data(pk + 35);
    const auto *pk_38 = buffer.data(pk + 38);
    const auto *pk_39 = buffer.data(pk + 39);
    const auto *pk_40 = buffer.data(pk + 40);
    const auto *pk_41 = buffer.data(pk + 41);
    const auto *pk_42 = buffer.data(pk + 42);
    const auto *pk_44 = buffer.data(pk + 44);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);
    const auto *dh0_3 = buffer.data(dh0 + 3);
    const auto *dh0_4 = buffer.data(dh0 + 4);
    const auto *dh0_5 = buffer.data(dh0 + 5);
    const auto *dh0_6 = buffer.data(dh0 + 6);
    const auto *dh0_7 = buffer.data(dh0 + 7);
    const auto *dh0_8 = buffer.data(dh0 + 8);
    const auto *dh0_9 = buffer.data(dh0 + 9);
    const auto *dh0_10 = buffer.data(dh0 + 10);
    const auto *dh0_11 = buffer.data(dh0 + 11);
    const auto *dh0_12 = buffer.data(dh0 + 12);
    const auto *dh0_13 = buffer.data(dh0 + 13);
    const auto *dh0_14 = buffer.data(dh0 + 14);
    const auto *dh0_15 = buffer.data(dh0 + 15);
    const auto *dh0_16 = buffer.data(dh0 + 16);
    const auto *dh0_17 = buffer.data(dh0 + 17);
    const auto *dh0_18 = buffer.data(dh0 + 18);
    const auto *dh0_19 = buffer.data(dh0 + 19);
    const auto *dh0_20 = buffer.data(dh0 + 20);
    const auto *dh0_21 = buffer.data(dh0 + 21);
    const auto *dh0_22 = buffer.data(dh0 + 22);
    const auto *dh0_23 = buffer.data(dh0 + 23);
    const auto *dh0_24 = buffer.data(dh0 + 24);
    const auto *dh0_25 = buffer.data(dh0 + 25);
    const auto *dh0_26 = buffer.data(dh0 + 26);
    const auto *dh0_27 = buffer.data(dh0 + 27);
    const auto *dh0_28 = buffer.data(dh0 + 28);
    const auto *dh0_29 = buffer.data(dh0 + 29);
    const auto *dh0_30 = buffer.data(dh0 + 30);
    const auto *dh0_31 = buffer.data(dh0 + 31);
    const auto *dh0_32 = buffer.data(dh0 + 32);
    const auto *dh0_33 = buffer.data(dh0 + 33);
    const auto *dh0_34 = buffer.data(dh0 + 34);
    const auto *dh0_35 = buffer.data(dh0 + 35);
    const auto *dh0_36 = buffer.data(dh0 + 36);
    const auto *dh0_37 = buffer.data(dh0 + 37);
    const auto *dh0_38 = buffer.data(dh0 + 38);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);
    const auto *dh1_3 = buffer.data(dh1 + 3);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_5 = buffer.data(dh1 + 5);
    const auto *dh1_6 = buffer.data(dh1 + 6);
    const auto *dh1_7 = buffer.data(dh1 + 7);
    const auto *dh1_8 = buffer.data(dh1 + 8);
    const auto *dh1_9 = buffer.data(dh1 + 9);
    const auto *dh1_10 = buffer.data(dh1 + 10);
    const auto *dh1_11 = buffer.data(dh1 + 11);
    const auto *dh1_12 = buffer.data(dh1 + 12);
    const auto *dh1_13 = buffer.data(dh1 + 13);
    const auto *dh1_14 = buffer.data(dh1 + 14);
    const auto *dh1_15 = buffer.data(dh1 + 15);
    const auto *dh1_16 = buffer.data(dh1 + 16);
    const auto *dh1_17 = buffer.data(dh1 + 17);
    const auto *dh1_18 = buffer.data(dh1 + 18);
    const auto *dh1_19 = buffer.data(dh1 + 19);
    const auto *dh1_20 = buffer.data(dh1 + 20);
    const auto *dh1_21 = buffer.data(dh1 + 21);
    const auto *dh1_22 = buffer.data(dh1 + 22);
    const auto *dh1_23 = buffer.data(dh1 + 23);
    const auto *dh1_24 = buffer.data(dh1 + 24);
    const auto *dh1_25 = buffer.data(dh1 + 25);
    const auto *dh1_26 = buffer.data(dh1 + 26);
    const auto *dh1_27 = buffer.data(dh1 + 27);
    const auto *dh1_28 = buffer.data(dh1 + 28);
    const auto *dh1_29 = buffer.data(dh1 + 29);
    const auto *dh1_30 = buffer.data(dh1 + 30);
    const auto *dh1_31 = buffer.data(dh1 + 31);
    const auto *dh1_32 = buffer.data(dh1 + 32);
    const auto *dh1_33 = buffer.data(dh1 + 33);
    const auto *dh1_34 = buffer.data(dh1 + 34);
    const auto *dh1_35 = buffer.data(dh1 + 35);
    const auto *dh1_36 = buffer.data(dh1 + 36);
    const auto *dh1_37 = buffer.data(dh1 + 37);
    const auto *dh1_38 = buffer.data(dh1 + 38);

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
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_23 = buffer.data(di + 23);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pi_0, dh0_0, dh1_0, di_0, \
                         di_1, di_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = pb_y[k] * di_0[k];

        t_2[k] = pb_z[k] * di_0[k];

        t_3[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_y[k] * di_1[k];

        t_4[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_z[k] * di_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, dh0_1, dh0_2, dh0_3, dh1_1, \
                         dh1_2, dh1_3, di_3, di_4, di_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * dh0_1[k]
                 - f_6 * dh1_1[k]
                 + pb_y[k] * di_3[k];

        t_6[k] = pb_z[k] * di_3[k];

        t_7[k] = pb_y[k] * di_4[k];

        t_8[k] = f_5 * dh0_2[k]
                 - f_6 * dh1_2[k]
                 + pb_z[k] * di_4[k];

        t_9[k] = f_7 * dh0_3[k]
                 - f_8 * dh1_3[k]
                 + pb_y[k] * di_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pb_y, pb_z, dh0_4, dh0_5, dh1_4, \
                         dh1_5, di_5, di_6, di_7, di_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_z[k] * di_5[k];

        t_11[k] = f_3 * dh0_4[k]
                  - f_4 * dh1_4[k]
                  + pb_y[k] * di_6[k];

        t_12[k] = pb_y[k] * di_7[k];

        t_13[k] = f_7 * dh0_4[k]
                  - f_8 * dh1_4[k]
                  + pb_z[k] * di_7[k];

        t_14[k] = f_9 * dh0_5[k]
                  - f_10 * dh1_5[k]
                  + pb_y[k] * di_8[k];

        t_15[k] = pb_z[k] * di_8[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_y, pb_z, dh0_6, dh0_7, dh1_6, dh1_7, di_9, \
                         di_10, di_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * dh0_6[k]
                  - f_6 * dh1_6[k]
                  + pb_y[k] * di_9[k];

        t_17[k] = f_3 * dh0_7[k]
                  - f_4 * dh1_7[k]
                  + pb_y[k] * di_10[k];

        t_18[k] = pb_y[k] * di_11[k];

        t_19[k] = f_9 * dh0_7[k]
                  - f_10 * dh1_7[k]
                  + pb_z[k] * di_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pb_y, pb_z, dh0_8, dh0_9, dh0_10, dh1_8, \
                         dh1_9, dh1_10, di_12, di_14, di_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * dh0_8[k]
                  - f_2 * dh1_8[k]
                  + pb_y[k] * di_12[k];

        t_21[k] = pb_z[k] * di_12[k];

        t_22[k] = f_9 * dh0_9[k]
                  - f_10 * dh1_9[k]
                  + pb_y[k] * di_14[k];

        t_23[k] = f_7 * dh0_10[k]
                  - f_8 * dh1_10[k]
                  + pb_y[k] * di_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_y, pb_y, pb_z, pk_0, dh0_11, dh0_12, \
                         dh1_11, dh1_12, di_16, di_17, di_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * dh0_11[k]
                  - f_6 * dh1_11[k]
                  + pb_y[k] * di_16[k];

        t_25[k] = f_3 * dh0_12[k]
                  - f_4 * dh1_12[k]
                  + pb_y[k] * di_17[k];

        t_26[k] = pb_y[k] * di_18[k];

        t_27[k] = f_1 * dh0_12[k]
                  - f_2 * dh1_12[k]
                  + pb_z[k] * di_18[k];

        t_28[k] = pa_y[k] * pk_0[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pa_x, pa_y, pi_1, pi_2, pi_3, pk_2, \
                         pk_4, pk_12, pk_13, pk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_11 * pi_1[k]
                  + pa_x[k] * pk_12[k];

        t_30[k] = pa_y[k] * pk_2[k];

        t_31[k] = f_12 * pi_2[k]
                  + pa_x[k] * pk_13[k];

        t_32[k] = pa_y[k] * pk_4[k];

        t_33[k] = f_13 * pi_3[k]
                  + pa_x[k] * pk_14[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, t_39, pa_x, pa_y, pi_4, pk_6, pk_8, \
                         pk_15, pk_17, pk_18, pk_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_y[k] * pk_6[k];

        t_35[k] = f_0 * pi_4[k]
                  + pa_x[k] * pk_15[k];

        t_36[k] = pa_y[k] * pk_8[k];

        t_37[k] = pa_x[k] * pk_17[k];

        t_38[k] = pa_x[k] * pk_18[k];

        t_39[k] = pa_x[k] * pk_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, pa_x, pa_z, pb_z, pi_0, pk_0, \
                         pk_1, pk_20, pk_21, pk_22, di_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_x[k] * pk_20[k];

        t_41[k] = pa_x[k] * pk_21[k];

        t_42[k] = pa_x[k] * pk_22[k];

        t_43[k] = pa_z[k] * pk_0[k];

        t_44[k] = f_14 * pi_0[k]
                  + pb_z[k] * di_21[k];

        t_45[k] = pa_z[k] * pk_1[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_x, pa_z, pi_8, pi_10, pi_13, pk_3, \
                         pk_5, pk_26, pk_28, pk_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_11 * pi_8[k]
                  + pa_x[k] * pk_26[k];

        t_47[k] = pa_z[k] * pk_3[k];

        t_48[k] = f_12 * pi_10[k]
                  + pa_x[k] * pk_28[k];

        t_49[k] = pa_z[k] * pk_5[k];

        t_50[k] = f_13 * pi_13[k]
                  + pa_x[k] * pk_31[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, pa_x, pa_z, pi_14, pk_7, pk_35, \
                         pk_38, pk_39, pk_40, pk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_z[k] * pk_7[k];

        t_52[k] = f_0 * pi_14[k]
                  + pa_x[k] * pk_35[k];

        t_53[k] = pa_x[k] * pk_38[k];

        t_54[k] = pa_x[k] * pk_39[k];

        t_55[k] = pa_x[k] * pk_40[k];

        t_56[k] = pa_x[k] * pk_41[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pa_x, pb_x, pb_z, pk_42, pk_44, dh0_13, \
                         dh0_14, dh1_13, dh1_14, di_23, di_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_x[k] * pk_42[k];

        t_58[k] = pa_x[k] * pk_44[k];

        t_59[k] = f_1 * dh0_13[k]
                  - f_2 * dh1_13[k]
                  + pb_x[k] * di_23[k];

        t_60[k] = pb_z[k] * di_23[k];

        t_61[k] = f_9 * dh0_14[k]
                  - f_10 * dh1_14[k]
                  + pb_x[k] * di_25[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pb_x, pb_z, dh0_15, dh0_16, dh0_17, dh1_15, \
                         dh1_16, dh1_17, di_25, di_26, di_27, di_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_9 * dh0_15[k]
                  - f_10 * dh1_15[k]
                  + pb_x[k] * di_26[k];

        t_63[k] = f_7 * dh0_16[k]
                  - f_8 * dh1_16[k]
                  + pb_x[k] * di_27[k];

        t_64[k] = pb_z[k] * di_25[k];

        t_65[k] = f_7 * dh0_17[k]
                  - f_8 * dh1_17[k]
                  + pb_x[k] * di_28[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pb_x, pb_z, dh0_18, dh0_19, dh0_20, dh1_18, \
                         dh1_19, dh1_20, di_27, di_29, di_30, di_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_5 * dh0_18[k]
                  - f_6 * dh1_18[k]
                  + pb_x[k] * di_29[k];

        t_67[k] = pb_z[k] * di_27[k];

        t_68[k] = f_5 * dh0_19[k]
                  - f_6 * dh1_19[k]
                  + pb_x[k] * di_30[k];

        t_69[k] = f_5 * dh0_20[k]
                  - f_6 * dh1_20[k]
                  + pb_x[k] * di_31[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_x, pb_z, dh0_21, dh0_23, dh0_24, dh1_21, \
                         dh1_23, dh1_24, di_29, di_32, di_33, di_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * dh0_21[k]
                  - f_4 * dh1_21[k]
                  + pb_x[k] * di_32[k];

        t_71[k] = pb_z[k] * di_29[k];

        t_72[k] = f_3 * dh0_23[k]
                  - f_4 * dh1_23[k]
                  + pb_x[k] * di_33[k];

        t_73[k] = f_3 * dh0_24[k]
                  - f_4 * dh1_24[k]
                  + pb_x[k] * di_34[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, pb_x, dh0_25, dh1_25, di_35, \
                         di_36, di_38, di_39, di_40, di_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_3 * dh0_25[k]
                  - f_4 * dh1_25[k]
                  + pb_x[k] * di_35[k];

        t_75[k] = pb_x[k] * di_36[k];

        t_76[k] = pb_x[k] * di_38[k];

        t_77[k] = pb_x[k] * di_39[k];

        t_78[k] = pb_x[k] * di_40[k];

        t_79[k] = pb_x[k] * di_41[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pb_y, pb_z, pi_5, dh0_21, dh0_22, dh1_21, \
                         dh1_22, di_36, di_37, di_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * pi_5[k]
                  + f_1 * dh0_21[k]
                  - f_2 * dh1_21[k]
                  + pb_y[k] * di_36[k];

        t_81[k] = pb_z[k] * di_36[k];

        t_82[k] = f_3 * dh0_21[k]
                  - f_4 * dh1_21[k]
                  + pb_z[k] * di_37[k];

        t_83[k] = f_5 * dh0_22[k]
                  - f_6 * dh1_22[k]
                  + pb_z[k] * di_38[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_y, pb_z, pi_6, dh0_23, dh0_24, dh0_25, \
                         dh1_23, dh1_24, dh1_25, di_39, di_40, di_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_7 * dh0_23[k]
                  - f_8 * dh1_23[k]
                  + pb_z[k] * di_39[k];

        t_85[k] = f_9 * dh0_24[k]
                  - f_10 * dh1_24[k]
                  + pb_z[k] * di_40[k];

        t_86[k] = f_0 * pi_6[k]
                  + pb_y[k] * di_41[k];

        t_87[k] = f_1 * dh0_25[k]
                  - f_2 * dh1_25[k]
                  + pb_z[k] * di_41[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, t_93, pa_y, pa_z, pk_12, pk_13, pk_14, \
                         pk_25, pk_26, pk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pa_y[k] * pk_25[k];

        t_89[k] = pa_z[k] * pk_12[k];

        t_90[k] = pa_y[k] * pk_26[k];

        t_91[k] = pa_z[k] * pk_13[k];

        t_92[k] = pa_y[k] * pk_28[k];

        t_93[k] = pa_z[k] * pk_14[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pa_y, pa_z, pi_9, pi_11, pi_12, pk_15, \
                         pk_29, pk_31, pk_32, pk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_0 * pi_9[k]
                  + pa_y[k] * pk_29[k];

        t_95[k] = pa_y[k] * pk_31[k];

        t_96[k] = pa_z[k] * pk_15[k];

        t_97[k] = f_13 * pi_11[k]
                  + pa_y[k] * pk_32[k];

        t_98[k] = f_0 * pi_12[k]
                  + pa_y[k] * pk_33[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pa_y, pa_z, pb_x, pk_17, pk_35, \
                         di_43, di_44, di_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_y[k] * pk_35[k];

        t_100[k] = pb_x[k] * di_43[k];

        t_101[k] = pb_x[k] * di_44[k];

        t_102[k] = pb_x[k] * di_45[k];

        t_103[k] = pa_z[k] * pk_17[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_y, pb_z, pi_5, pi_16, pi_17, pi_18, \
                         pk_39, pk_40, pk_41, di_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_14 * pi_5[k]
                   + pb_z[k] * di_42[k];

        t_105[k] = f_11 * pi_16[k]
                   + pa_y[k] * pk_39[k];

        t_106[k] = f_12 * pi_17[k]
                   + pa_y[k] * pk_40[k];

        t_107[k] = f_13 * pi_18[k]
                   + pa_y[k] * pk_41[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_y, pb_x, pb_y, pi_19, pi_20, \
                         pk_42, pk_44, dh0_26, dh1_26, di_46, di_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_0 * pi_19[k]
                   + pa_y[k] * pk_42[k];

        t_109[k] = f_14 * pi_20[k]
                   + pb_y[k] * di_46[k];

        t_110[k] = pa_y[k] * pk_44[k];

        t_111[k] = f_1 * dh0_26[k]
                   - f_2 * dh1_26[k]
                   + pb_x[k] * di_47[k];

        t_112[k] = pb_y[k] * di_47[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pb_x, pb_z, pi_7, dh0_27, dh0_28, dh1_27, \
                         dh1_28, di_47, di_49, di_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_0 * pi_7[k]
                   + pb_z[k] * di_47[k];

        t_114[k] = f_9 * dh0_27[k]
                   - f_10 * dh1_27[k]
                   + pb_x[k] * di_49[k];

        t_115[k] = f_9 * dh0_28[k]
                   - f_10 * dh1_28[k]
                   + pb_x[k] * di_50[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_x, pb_y, dh0_29, dh0_30, dh0_31, \
                         dh1_29, dh1_30, dh1_31, di_50, di_51, di_52, \
                         di_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_7 * dh0_29[k]
                   - f_8 * dh1_29[k]
                   + pb_x[k] * di_51[k];

        t_117[k] = pb_y[k] * di_50[k];

        t_118[k] = f_7 * dh0_30[k]
                   - f_8 * dh1_30[k]
                   + pb_x[k] * di_52[k];

        t_119[k] = f_5 * dh0_31[k]
                   - f_6 * dh1_31[k]
                   + pb_x[k] * di_53[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_x, pb_y, dh0_32, dh0_33, dh0_34, \
                         dh1_32, dh1_33, dh1_34, di_52, di_54, di_55, \
                         di_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_5 * dh0_32[k]
                   - f_6 * dh1_32[k]
                   + pb_x[k] * di_54[k];

        t_121[k] = pb_y[k] * di_52[k];

        t_122[k] = f_5 * dh0_33[k]
                   - f_6 * dh1_33[k]
                   + pb_x[k] * di_55[k];

        t_123[k] = f_3 * dh0_34[k]
                   - f_4 * dh1_34[k]
                   + pb_x[k] * di_56[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pb_x, pb_y, dh0_35, dh0_36, dh0_38, \
                         dh1_35, dh1_36, dh1_38, di_55, di_57, di_58, \
                         di_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_3 * dh0_35[k]
                   - f_4 * dh1_35[k]
                   + pb_x[k] * di_57[k];

        t_125[k] = f_3 * dh0_36[k]
                   - f_4 * dh1_36[k]
                   + pb_x[k] * di_58[k];

        t_126[k] = pb_y[k] * di_55[k];

        t_127[k] = f_3 * dh0_38[k]
                   - f_4 * dh1_38[k]
                   + pb_x[k] * di_59[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, pb_x, pb_y, dh0_34, dh1_34, \
                         di_60, di_61, di_62, di_63, di_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = pb_x[k] * di_60[k];

        t_129[k] = pb_x[k] * di_61[k];

        t_130[k] = pb_x[k] * di_62[k];

        t_131[k] = pb_x[k] * di_63[k];

        t_132[k] = pb_x[k] * di_65[k];

        t_133[k] = f_1 * dh0_34[k]
                   - f_2 * dh1_34[k]
                   + pb_y[k] * di_60[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, pb_y, pb_z, pi_15, dh0_35, dh0_36, dh1_35, \
                         dh1_36, di_60, di_61, di_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_0 * pi_15[k]
                   + pb_z[k] * di_60[k];

        t_135[k] = f_9 * dh0_35[k]
                   - f_10 * dh1_35[k]
                   + pb_y[k] * di_61[k];

        t_136[k] = f_7 * dh0_36[k]
                   - f_8 * dh1_36[k]
                   + pb_y[k] * di_62[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pb_y, pb_z, pi_20, dh0_37, dh0_38, \
                         dh1_37, dh1_38, di_63, di_64, di_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_5 * dh0_37[k]
                   - f_6 * dh1_37[k]
                   + pb_y[k] * di_63[k];

        t_138[k] = f_3 * dh0_38[k]
                   - f_4 * dh1_38[k]
                   + pb_y[k] * di_64[k];

        t_139[k] = pb_y[k] * di_65[k];

        t_140[k] = f_0 * pi_20[k]
                   + f_1 * dh0_38[k]
                   - f_2 * dh1_38[k]
                   + pb_z[k] * di_65[k];
    }
}

auto
compute_prim_dk_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pi, const size_t pk,
                                     const size_t dh0, const size_t dh1, const size_t di,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_11 = buffer.data(pi + 11);
    const auto *pi_12 = buffer.data(pi + 12);
    const auto *pi_13 = buffer.data(pi + 13);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_16 = buffer.data(pi + 16);
    const auto *pi_17 = buffer.data(pi + 17);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);
    const auto *pk_3 = buffer.data(pk + 3);
    const auto *pk_4 = buffer.data(pk + 4);
    const auto *pk_5 = buffer.data(pk + 5);
    const auto *pk_6 = buffer.data(pk + 6);
    const auto *pk_7 = buffer.data(pk + 7);
    const auto *pk_8 = buffer.data(pk + 8);
    const auto *pk_9 = buffer.data(pk + 9);
    const auto *pk_10 = buffer.data(pk + 10);
    const auto *pk_11 = buffer.data(pk + 11);
    const auto *pk_12 = buffer.data(pk + 12);
    const auto *pk_13 = buffer.data(pk + 13);
    const auto *pk_14 = buffer.data(pk + 14);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);
    const auto *dh0_3 = buffer.data(dh0 + 3);
    const auto *dh0_4 = buffer.data(dh0 + 4);
    const auto *dh0_5 = buffer.data(dh0 + 5);
    const auto *dh0_6 = buffer.data(dh0 + 6);
    const auto *dh0_7 = buffer.data(dh0 + 7);
    const auto *dh0_8 = buffer.data(dh0 + 8);
    const auto *dh0_9 = buffer.data(dh0 + 9);
    const auto *dh0_10 = buffer.data(dh0 + 10);
    const auto *dh0_11 = buffer.data(dh0 + 11);
    const auto *dh0_12 = buffer.data(dh0 + 12);
    const auto *dh0_13 = buffer.data(dh0 + 13);
    const auto *dh0_14 = buffer.data(dh0 + 14);
    const auto *dh0_15 = buffer.data(dh0 + 15);
    const auto *dh0_16 = buffer.data(dh0 + 16);
    const auto *dh0_17 = buffer.data(dh0 + 17);
    const auto *dh0_18 = buffer.data(dh0 + 18);
    const auto *dh0_19 = buffer.data(dh0 + 19);
    const auto *dh0_20 = buffer.data(dh0 + 20);
    const auto *dh0_21 = buffer.data(dh0 + 21);
    const auto *dh0_22 = buffer.data(dh0 + 22);
    const auto *dh0_23 = buffer.data(dh0 + 23);
    const auto *dh0_24 = buffer.data(dh0 + 24);
    const auto *dh0_25 = buffer.data(dh0 + 25);
    const auto *dh0_26 = buffer.data(dh0 + 26);
    const auto *dh0_27 = buffer.data(dh0 + 27);
    const auto *dh0_28 = buffer.data(dh0 + 28);
    const auto *dh0_29 = buffer.data(dh0 + 29);
    const auto *dh0_30 = buffer.data(dh0 + 30);
    const auto *dh0_31 = buffer.data(dh0 + 31);
    const auto *dh0_32 = buffer.data(dh0 + 32);
    const auto *dh0_33 = buffer.data(dh0 + 33);
    const auto *dh0_34 = buffer.data(dh0 + 34);
    const auto *dh0_35 = buffer.data(dh0 + 35);
    const auto *dh0_36 = buffer.data(dh0 + 36);
    const auto *dh0_37 = buffer.data(dh0 + 37);
    const auto *dh0_38 = buffer.data(dh0 + 38);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);
    const auto *dh1_3 = buffer.data(dh1 + 3);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_5 = buffer.data(dh1 + 5);
    const auto *dh1_6 = buffer.data(dh1 + 6);
    const auto *dh1_7 = buffer.data(dh1 + 7);
    const auto *dh1_8 = buffer.data(dh1 + 8);
    const auto *dh1_9 = buffer.data(dh1 + 9);
    const auto *dh1_10 = buffer.data(dh1 + 10);
    const auto *dh1_11 = buffer.data(dh1 + 11);
    const auto *dh1_12 = buffer.data(dh1 + 12);
    const auto *dh1_13 = buffer.data(dh1 + 13);
    const auto *dh1_14 = buffer.data(dh1 + 14);
    const auto *dh1_15 = buffer.data(dh1 + 15);
    const auto *dh1_16 = buffer.data(dh1 + 16);
    const auto *dh1_17 = buffer.data(dh1 + 17);
    const auto *dh1_18 = buffer.data(dh1 + 18);
    const auto *dh1_19 = buffer.data(dh1 + 19);
    const auto *dh1_20 = buffer.data(dh1 + 20);
    const auto *dh1_21 = buffer.data(dh1 + 21);
    const auto *dh1_22 = buffer.data(dh1 + 22);
    const auto *dh1_23 = buffer.data(dh1 + 23);
    const auto *dh1_24 = buffer.data(dh1 + 24);
    const auto *dh1_25 = buffer.data(dh1 + 25);
    const auto *dh1_26 = buffer.data(dh1 + 26);
    const auto *dh1_27 = buffer.data(dh1 + 27);
    const auto *dh1_28 = buffer.data(dh1 + 28);
    const auto *dh1_29 = buffer.data(dh1 + 29);
    const auto *dh1_30 = buffer.data(dh1 + 30);
    const auto *dh1_31 = buffer.data(dh1 + 31);
    const auto *dh1_32 = buffer.data(dh1 + 32);
    const auto *dh1_33 = buffer.data(dh1 + 33);
    const auto *dh1_34 = buffer.data(dh1 + 34);
    const auto *dh1_35 = buffer.data(dh1 + 35);
    const auto *dh1_36 = buffer.data(dh1 + 36);
    const auto *dh1_37 = buffer.data(dh1 + 37);
    const auto *dh1_38 = buffer.data(dh1 + 38);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pi_0, dh0_0, dh1_0, di_0, \
                         di_1, di_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = pb_y[k] * di_0[k];

        t_2[k] = pb_z[k] * di_0[k];

        t_3[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_y[k] * di_1[k];

        t_4[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_z[k] * di_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, dh0_1, dh0_2, dh0_3, dh1_1, dh1_2, \
                         dh1_3, di_3, di_4, di_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * dh0_1[k]
                 - f_6 * dh1_1[k]
                 + pb_y[k] * di_3[k];

        t_6[k] = pb_y[k] * di_4[k];

        t_7[k] = f_5 * dh0_2[k]
                 - f_6 * dh1_2[k]
                 + pb_z[k] * di_4[k];

        t_8[k] = f_7 * dh0_3[k]
                 - f_8 * dh1_3[k]
                 + pb_y[k] * di_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, dh0_4, dh0_5, dh1_4, dh1_5, di_6, \
                         di_7, di_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * dh0_4[k]
                 - f_4 * dh1_4[k]
                 + pb_y[k] * di_6[k];

        t_10[k] = pb_y[k] * di_7[k];

        t_11[k] = f_7 * dh0_4[k]
                  - f_8 * dh1_4[k]
                  + pb_z[k] * di_7[k];

        t_12[k] = f_9 * dh0_5[k]
                  - f_10 * dh1_5[k]
                  + pb_y[k] * di_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pb_y, pb_z, dh0_6, dh0_7, dh1_6, dh1_7, di_9, \
                         di_10, di_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * dh0_6[k]
                  - f_6 * dh1_6[k]
                  + pb_y[k] * di_9[k];

        t_14[k] = f_3 * dh0_7[k]
                  - f_4 * dh1_7[k]
                  + pb_y[k] * di_10[k];

        t_15[k] = pb_y[k] * di_11[k];

        t_16[k] = f_9 * dh0_7[k]
                  - f_10 * dh1_7[k]
                  + pb_z[k] * di_11[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_y, dh0_8, dh0_9, dh0_10, dh1_8, dh1_9, dh1_10, \
                         di_12, di_13, di_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * dh0_8[k]
                  - f_2 * dh1_8[k]
                  + pb_y[k] * di_12[k];

        t_18[k] = f_9 * dh0_9[k]
                  - f_10 * dh1_9[k]
                  + pb_y[k] * di_13[k];

        t_19[k] = f_7 * dh0_10[k]
                  - f_8 * dh1_10[k]
                  + pb_y[k] * di_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pb_y, pb_z, pk_0, dh0_11, dh0_12, \
                         dh1_11, dh1_12, di_15, di_16, di_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * dh0_11[k]
                  - f_6 * dh1_11[k]
                  + pb_y[k] * di_15[k];

        t_21[k] = f_3 * dh0_12[k]
                  - f_4 * dh1_12[k]
                  + pb_y[k] * di_16[k];

        t_22[k] = pb_y[k] * di_17[k];

        t_23[k] = f_1 * dh0_12[k]
                  - f_2 * dh1_12[k]
                  + pb_z[k] * di_17[k];

        t_24[k] = pa_y[k] * pk_0[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_x, pi_1, pi_2, pi_3, pi_4, pk_1, \
                         pk_2, pk_3, pk_4, pk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_11 * pi_1[k]
                  + pa_x[k] * pk_1[k];

        t_26[k] = f_12 * pi_2[k]
                  + pa_x[k] * pk_2[k];

        t_27[k] = f_13 * pi_3[k]
                  + pa_x[k] * pk_3[k];

        t_28[k] = f_0 * pi_4[k]
                  + pa_x[k] * pk_4[k];

        t_29[k] = pa_x[k] * pk_5[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_z, pb_z, pi_0, pi_8, pi_9, pk_0, \
                         pk_6, pk_7, di_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_z[k] * pk_0[k];

        t_31[k] = f_14 * pi_0[k]
                  + pb_z[k] * di_18[k];

        t_32[k] = f_11 * pi_8[k]
                  + pa_x[k] * pk_6[k];

        t_33[k] = f_12 * pi_9[k]
                  + pa_x[k] * pk_7[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pi_10, pi_11, pk_8, pk_9, pk_14, \
                         dh0_13, dh1_13, di_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_13 * pi_10[k]
                  + pa_x[k] * pk_8[k];

        t_35[k] = f_0 * pi_11[k]
                  + pa_x[k] * pk_9[k];

        t_36[k] = pa_x[k] * pk_14[k];

        t_37[k] = f_1 * dh0_13[k]
                  - f_2 * dh1_13[k]
                  + pb_x[k] * di_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pb_x, dh0_14, dh0_15, dh0_16, dh1_14, dh1_15, \
                         dh1_16, di_20, di_21, di_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_9 * dh0_14[k]
                  - f_10 * dh1_14[k]
                  + pb_x[k] * di_20[k];

        t_39[k] = f_9 * dh0_15[k]
                  - f_10 * dh1_15[k]
                  + pb_x[k] * di_21[k];

        t_40[k] = f_7 * dh0_16[k]
                  - f_8 * dh1_16[k]
                  + pb_x[k] * di_22[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pb_x, dh0_17, dh0_18, dh0_19, dh1_17, dh1_18, \
                         dh1_19, di_23, di_24, di_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * dh0_17[k]
                  - f_8 * dh1_17[k]
                  + pb_x[k] * di_23[k];

        t_42[k] = f_5 * dh0_18[k]
                  - f_6 * dh1_18[k]
                  + pb_x[k] * di_24[k];

        t_43[k] = f_5 * dh0_19[k]
                  - f_6 * dh1_19[k]
                  + pb_x[k] * di_25[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pb_x, dh0_20, dh0_21, dh0_23, dh1_20, dh1_21, \
                         dh1_23, di_26, di_27, di_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_5 * dh0_20[k]
                  - f_6 * dh1_20[k]
                  + pb_x[k] * di_26[k];

        t_45[k] = f_3 * dh0_21[k]
                  - f_4 * dh1_21[k]
                  + pb_x[k] * di_27[k];

        t_46[k] = f_3 * dh0_23[k]
                  - f_4 * dh1_23[k]
                  + pb_x[k] * di_28[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pb_x, dh0_24, dh0_25, dh1_24, dh1_25, \
                         di_29, di_30, di_31, di_33, di_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * dh0_24[k]
                  - f_4 * dh1_24[k]
                  + pb_x[k] * di_29[k];

        t_48[k] = f_3 * dh0_25[k]
                  - f_4 * dh1_25[k]
                  + pb_x[k] * di_30[k];

        t_49[k] = pb_x[k] * di_31[k];

        t_50[k] = pb_x[k] * di_33[k];

        t_51[k] = pb_x[k] * di_34[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pb_x, pb_y, pb_z, pi_5, dh0_21, dh1_21, \
                         di_31, di_32, di_35, di_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_x[k] * di_35[k];

        t_53[k] = pb_x[k] * di_36[k];

        t_54[k] = f_0 * pi_5[k]
                  + f_1 * dh0_21[k]
                  - f_2 * dh1_21[k]
                  + pb_y[k] * di_31[k];

        t_55[k] = pb_z[k] * di_31[k];

        t_56[k] = f_3 * dh0_21[k]
                  - f_4 * dh1_21[k]
                  + pb_z[k] * di_32[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_z, dh0_22, dh0_23, dh0_24, dh1_22, dh1_23, \
                         dh1_24, di_33, di_34, di_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_5 * dh0_22[k]
                  - f_6 * dh1_22[k]
                  + pb_z[k] * di_33[k];

        t_58[k] = f_7 * dh0_23[k]
                  - f_8 * dh1_23[k]
                  + pb_z[k] * di_34[k];

        t_59[k] = f_9 * dh0_24[k]
                  - f_10 * dh1_24[k]
                  + pb_z[k] * di_35[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_y, pb_z, pi_5, pi_6, pk_5, dh0_25, \
                         dh1_25, di_36, di_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * pi_6[k]
                  + pb_y[k] * di_36[k];

        t_61[k] = f_1 * dh0_25[k]
                  - f_2 * dh1_25[k]
                  + pb_z[k] * di_36[k];

        t_62[k] = pa_z[k] * pk_5[k];

        t_63[k] = f_14 * pi_5[k]
                  + pb_z[k] * di_37[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pi_13, pi_14, pi_15, pi_16, pk_10, \
                         pk_11, pk_12, pk_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_11 * pi_13[k]
                  + pa_y[k] * pk_10[k];

        t_65[k] = f_12 * pi_14[k]
                  + pa_y[k] * pk_11[k];

        t_66[k] = f_13 * pi_15[k]
                  + pa_y[k] * pk_12[k];

        t_67[k] = f_0 * pi_16[k]
                  + pa_y[k] * pk_13[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_y, pb_x, pb_y, pb_z, pi_7, pi_17, pk_14, \
                         dh0_26, dh1_26, di_38, di_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_14 * pi_17[k]
                  + pb_y[k] * di_38[k];

        t_69[k] = pa_y[k] * pk_14[k];

        t_70[k] = f_1 * dh0_26[k]
                  - f_2 * dh1_26[k]
                  + pb_x[k] * di_39[k];

        t_71[k] = f_0 * pi_7[k]
                  + pb_z[k] * di_39[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_x, dh0_27, dh0_28, dh0_29, dh1_27, dh1_28, \
                         dh1_29, di_40, di_41, di_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_9 * dh0_27[k]
                  - f_10 * dh1_27[k]
                  + pb_x[k] * di_40[k];

        t_73[k] = f_9 * dh0_28[k]
                  - f_10 * dh1_28[k]
                  + pb_x[k] * di_41[k];

        t_74[k] = f_7 * dh0_29[k]
                  - f_8 * dh1_29[k]
                  + pb_x[k] * di_42[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_x, dh0_30, dh0_31, dh0_32, dh1_30, dh1_31, \
                         dh1_32, di_43, di_44, di_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_7 * dh0_30[k]
                  - f_8 * dh1_30[k]
                  + pb_x[k] * di_43[k];

        t_76[k] = f_5 * dh0_31[k]
                  - f_6 * dh1_31[k]
                  + pb_x[k] * di_44[k];

        t_77[k] = f_5 * dh0_32[k]
                  - f_6 * dh1_32[k]
                  + pb_x[k] * di_45[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pb_x, dh0_33, dh0_34, dh0_35, dh1_33, dh1_34, \
                         dh1_35, di_46, di_47, di_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_5 * dh0_33[k]
                  - f_6 * dh1_33[k]
                  + pb_x[k] * di_46[k];

        t_79[k] = f_3 * dh0_34[k]
                  - f_4 * dh1_34[k]
                  + pb_x[k] * di_47[k];

        t_80[k] = f_3 * dh0_35[k]
                  - f_4 * dh1_35[k]
                  + pb_x[k] * di_48[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, pb_x, dh0_36, dh0_38, dh1_36, dh1_38, \
                         di_49, di_50, di_51, di_52, di_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_3 * dh0_36[k]
                  - f_4 * dh1_36[k]
                  + pb_x[k] * di_49[k];

        t_82[k] = f_3 * dh0_38[k]
                  - f_4 * dh1_38[k]
                  + pb_x[k] * di_50[k];

        t_83[k] = pb_x[k] * di_51[k];

        t_84[k] = pb_x[k] * di_52[k];

        t_85[k] = pb_x[k] * di_53[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_x, pb_y, pb_z, pi_12, dh0_34, dh1_34, \
                         di_51, di_54, di_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_x[k] * di_54[k];

        t_87[k] = pb_x[k] * di_56[k];

        t_88[k] = f_1 * dh0_34[k]
                  - f_2 * dh1_34[k]
                  + pb_y[k] * di_51[k];

        t_89[k] = f_0 * pi_12[k]
                  + pb_z[k] * di_51[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_y, dh0_35, dh0_36, dh0_37, dh1_35, dh1_36, \
                         dh1_37, di_52, di_53, di_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_9 * dh0_35[k]
                  - f_10 * dh1_35[k]
                  + pb_y[k] * di_52[k];

        t_91[k] = f_7 * dh0_36[k]
                  - f_8 * dh1_36[k]
                  + pb_y[k] * di_53[k];

        t_92[k] = f_5 * dh0_37[k]
                  - f_6 * dh1_37[k]
                  + pb_y[k] * di_54[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pb_y, pb_z, pi_17, dh0_38, dh1_38, di_55, \
                         di_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_3 * dh0_38[k]
                  - f_4 * dh1_38[k]
                  + pb_y[k] * di_55[k];

        t_94[k] = pb_y[k] * di_56[k];

        t_95[k] = f_0 * pi_17[k]
                  + f_1 * dh0_38[k]
                  - f_2 * dh1_38[k]
                  + pb_z[k] * di_56[k];
    }
}

auto
compute_prim_dk_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pi, const size_t dh0, const size_t dh1,
                                     const size_t di, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pi_0, pi_1, dh0_0, dh0_1, dh1_0, dh1_1, di_0, \
                         di_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_0 * pi_1[k]
                 + f_1 * dh0_1[k]
                 - f_2 * dh1_1[k]
                 + pb_y[k] * di_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pi_2, dh0_2, dh1_2, di_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pi_2[k]
                 + f_1 * dh0_2[k]
                 - f_2 * dh1_2[k]
                 + pb_z[k] * di_2[k];
    }
}

auto
compute_prim_dk_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pi, const size_t dh0, const size_t dh1,
                                     const size_t di, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_26 = buffer.data(dh1 + 26);
    const auto *dh1_47 = buffer.data(dh1 + 47);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_74 = buffer.data(di + 74);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pi_0, pi_1, dh0_0, dh0_1, dh1_0, dh1_26, di_0, \
                         di_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_0 * pi_1[k]
                 + f_1 * dh0_1[k]
                 - f_2 * dh1_26[k]
                 + pb_y[k] * di_41[k];
    }

#pragma omp simd aligned(t_2, pb_z, pi_2, dh0_2, dh1_47, di_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pi_2[k]
                 + f_1 * dh0_2[k]
                 - f_2 * dh1_47[k]
                 + pb_z[k] * di_74[k];
    }
}

auto
compute_prim_dk_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pi, const size_t dh0, const size_t dh1,
                                     const size_t di, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_11 = buffer.data(pi + 11);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);
    const auto *dh0_3 = buffer.data(dh0 + 3);
    const auto *dh0_4 = buffer.data(dh0 + 4);
    const auto *dh0_5 = buffer.data(dh0 + 5);
    const auto *dh0_7 = buffer.data(dh0 + 7);
    const auto *dh0_8 = buffer.data(dh0 + 8);
    const auto *dh0_9 = buffer.data(dh0 + 9);
    const auto *dh0_11 = buffer.data(dh0 + 11);
    const auto *dh0_12 = buffer.data(dh0 + 12);
    const auto *dh0_13 = buffer.data(dh0 + 13);
    const auto *dh0_14 = buffer.data(dh0 + 14);
    const auto *dh0_16 = buffer.data(dh0 + 16);
    const auto *dh0_18 = buffer.data(dh0 + 18);
    const auto *dh0_19 = buffer.data(dh0 + 19);
    const auto *dh0_20 = buffer.data(dh0 + 20);
    const auto *dh0_22 = buffer.data(dh0 + 22);
    const auto *dh0_23 = buffer.data(dh0 + 23);
    const auto *dh0_24 = buffer.data(dh0 + 24);
    const auto *dh0_25 = buffer.data(dh0 + 25);
    const auto *dh0_26 = buffer.data(dh0 + 26);
    const auto *dh0_27 = buffer.data(dh0 + 27);
    const auto *dh0_28 = buffer.data(dh0 + 28);
    const auto *dh0_29 = buffer.data(dh0 + 29);
    const auto *dh0_30 = buffer.data(dh0 + 30);
    const auto *dh0_33 = buffer.data(dh0 + 33);
    const auto *dh0_35 = buffer.data(dh0 + 35);
    const auto *dh0_36 = buffer.data(dh0 + 36);
    const auto *dh0_37 = buffer.data(dh0 + 37);
    const auto *dh0_39 = buffer.data(dh0 + 39);
    const auto *dh0_40 = buffer.data(dh0 + 40);
    const auto *dh0_41 = buffer.data(dh0 + 41);
    const auto *dh0_42 = buffer.data(dh0 + 42);
    const auto *dh0_43 = buffer.data(dh0 + 43);
    const auto *dh0_44 = buffer.data(dh0 + 44);
    const auto *dh0_45 = buffer.data(dh0 + 45);
    const auto *dh0_46 = buffer.data(dh0 + 46);
    const auto *dh0_47 = buffer.data(dh0 + 47);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);
    const auto *dh1_3 = buffer.data(dh1 + 3);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_5 = buffer.data(dh1 + 5);
    const auto *dh1_6 = buffer.data(dh1 + 6);
    const auto *dh1_7 = buffer.data(dh1 + 7);
    const auto *dh1_8 = buffer.data(dh1 + 8);
    const auto *dh1_10 = buffer.data(dh1 + 10);
    const auto *dh1_11 = buffer.data(dh1 + 11);
    const auto *dh1_12 = buffer.data(dh1 + 12);
    const auto *dh1_13 = buffer.data(dh1 + 13);
    const auto *dh1_15 = buffer.data(dh1 + 15);
    const auto *dh1_17 = buffer.data(dh1 + 17);
    const auto *dh1_18 = buffer.data(dh1 + 18);
    const auto *dh1_19 = buffer.data(dh1 + 19);
    const auto *dh1_20 = buffer.data(dh1 + 20);
    const auto *dh1_21 = buffer.data(dh1 + 21);
    const auto *dh1_22 = buffer.data(dh1 + 22);
    const auto *dh1_23 = buffer.data(dh1 + 23);
    const auto *dh1_24 = buffer.data(dh1 + 24);
    const auto *dh1_25 = buffer.data(dh1 + 25);
    const auto *dh1_26 = buffer.data(dh1 + 26);
    const auto *dh1_27 = buffer.data(dh1 + 27);
    const auto *dh1_28 = buffer.data(dh1 + 28);
    const auto *dh1_31 = buffer.data(dh1 + 31);
    const auto *dh1_33 = buffer.data(dh1 + 33);
    const auto *dh1_34 = buffer.data(dh1 + 34);
    const auto *dh1_35 = buffer.data(dh1 + 35);
    const auto *dh1_36 = buffer.data(dh1 + 36);
    const auto *dh1_37 = buffer.data(dh1 + 37);
    const auto *dh1_38 = buffer.data(dh1 + 38);
    const auto *dh1_39 = buffer.data(dh1 + 39);
    const auto *dh1_40 = buffer.data(dh1 + 40);
    const auto *dh1_41 = buffer.data(dh1 + 41);
    const auto *dh1_42 = buffer.data(dh1 + 42);
    const auto *dh1_43 = buffer.data(dh1 + 43);
    const auto *dh1_44 = buffer.data(dh1 + 44);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pi_0, dh0_0, dh0_1, dh1_0, \
                         dh1_1, di_0, di_1, di_2, di_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_y[k] * di_1[k];

        t_2[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_z[k] * di_2[k];

        t_3[k] = f_5 * dh0_1[k]
                 - f_6 * dh1_1[k]
                 + pb_y[k] * di_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, dh0_2, dh0_3, dh0_4, dh1_2, dh1_3, \
                         dh1_4, di_4, di_5, di_6, di_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * dh0_2[k]
                 - f_6 * dh1_2[k]
                 + pb_z[k] * di_4[k];

        t_5[k] = f_7 * dh0_3[k]
                 - f_8 * dh1_3[k]
                 + pb_y[k] * di_5[k];

        t_6[k] = f_3 * dh0_4[k]
                 - f_4 * dh1_4[k]
                 + pb_y[k] * di_6[k];

        t_7[k] = f_7 * dh0_4[k]
                 - f_8 * dh1_4[k]
                 + pb_z[k] * di_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, dh0_5, dh0_7, dh0_8, dh1_5, dh1_6, \
                         dh1_7, di_8, di_9, di_10, di_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_9 * dh0_5[k]
                 - f_10 * dh1_5[k]
                 + pb_y[k] * di_8[k];

        t_9[k] = f_5 * dh0_7[k]
                 - f_6 * dh1_6[k]
                 + pb_y[k] * di_9[k];

        t_10[k] = f_3 * dh0_8[k]
                  - f_4 * dh1_7[k]
                  + pb_y[k] * di_10[k];

        t_11[k] = f_9 * dh0_8[k]
                  - f_10 * dh1_7[k]
                  + pb_z[k] * di_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, dh0_9, dh0_11, dh0_12, dh1_8, dh1_10, dh1_11, \
                         di_12, di_13, di_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_1 * dh0_9[k]
                  - f_2 * dh1_8[k]
                  + pb_y[k] * di_12[k];

        t_13[k] = f_9 * dh0_11[k]
                  - f_10 * dh1_10[k]
                  + pb_y[k] * di_13[k];

        t_14[k] = f_7 * dh0_12[k]
                  - f_8 * dh1_11[k]
                  + pb_y[k] * di_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_y, pb_z, dh0_13, dh0_14, dh1_12, dh1_13, di_15, \
                         di_16, di_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * dh0_13[k]
                  - f_6 * dh1_12[k]
                  + pb_y[k] * di_15[k];

        t_16[k] = f_3 * dh0_14[k]
                  - f_4 * dh1_13[k]
                  + pb_y[k] * di_16[k];

        t_17[k] = f_1 * dh0_14[k]
                  - f_2 * dh1_13[k]
                  + pb_z[k] * di_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_x, dh0_16, dh0_18, dh0_19, dh1_15, dh1_17, \
                         dh1_18, di_18, di_19, di_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_1 * dh0_16[k]
                  - f_2 * dh1_15[k]
                  + pb_x[k] * di_18[k];

        t_19[k] = f_9 * dh0_18[k]
                  - f_10 * dh1_17[k]
                  + pb_x[k] * di_19[k];

        t_20[k] = f_9 * dh0_19[k]
                  - f_10 * dh1_18[k]
                  + pb_x[k] * di_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_x, dh0_20, dh0_22, dh0_23, dh1_19, dh1_20, \
                         dh1_21, di_21, di_22, di_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * dh0_20[k]
                  - f_8 * dh1_19[k]
                  + pb_x[k] * di_21[k];

        t_22[k] = f_7 * dh0_22[k]
                  - f_8 * dh1_20[k]
                  + pb_x[k] * di_22[k];

        t_23[k] = f_5 * dh0_23[k]
                  - f_6 * dh1_21[k]
                  + pb_x[k] * di_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_x, dh0_24, dh0_25, dh0_26, dh1_22, dh1_23, \
                         dh1_24, di_24, di_25, di_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * dh0_24[k]
                  - f_6 * dh1_22[k]
                  + pb_x[k] * di_24[k];

        t_25[k] = f_5 * dh0_25[k]
                  - f_6 * dh1_23[k]
                  + pb_x[k] * di_25[k];

        t_26[k] = f_3 * dh0_26[k]
                  - f_4 * dh1_24[k]
                  + pb_x[k] * di_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_x, dh0_28, dh0_29, dh0_30, dh1_26, dh1_27, \
                         dh1_28, di_27, di_28, di_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_3 * dh0_28[k]
                  - f_4 * dh1_26[k]
                  + pb_x[k] * di_27[k];

        t_28[k] = f_3 * dh0_29[k]
                  - f_4 * dh1_27[k]
                  + pb_x[k] * di_28[k];

        t_29[k] = f_3 * dh0_30[k]
                  - f_4 * dh1_28[k]
                  + pb_x[k] * di_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_y, pb_z, pi_4, dh0_26, dh0_27, dh1_24, dh1_25, \
                         di_30, di_31, di_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * pi_4[k]
                  + f_1 * dh0_26[k]
                  - f_2 * dh1_24[k]
                  + pb_y[k] * di_30[k];

        t_31[k] = f_3 * dh0_26[k]
                  - f_4 * dh1_24[k]
                  + pb_z[k] * di_31[k];

        t_32[k] = f_5 * dh0_27[k]
                  - f_6 * dh1_25[k]
                  + pb_z[k] * di_32[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_z, dh0_28, dh0_29, dh0_30, dh1_26, dh1_27, \
                         dh1_28, di_33, di_34, di_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_7 * dh0_28[k]
                  - f_8 * dh1_26[k]
                  + pb_z[k] * di_33[k];

        t_34[k] = f_9 * dh0_29[k]
                  - f_10 * dh1_27[k]
                  + pb_z[k] * di_34[k];

        t_35[k] = f_1 * dh0_30[k]
                  - f_2 * dh1_28[k]
                  + pb_z[k] * di_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_x, dh0_33, dh0_35, dh0_36, dh1_31, dh1_33, \
                         dh1_34, di_36, di_37, di_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * dh0_33[k]
                  - f_2 * dh1_31[k]
                  + pb_x[k] * di_36[k];

        t_37[k] = f_9 * dh0_35[k]
                  - f_10 * dh1_33[k]
                  + pb_x[k] * di_37[k];

        t_38[k] = f_9 * dh0_36[k]
                  - f_10 * dh1_34[k]
                  + pb_x[k] * di_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_x, dh0_37, dh0_39, dh0_40, dh1_35, dh1_36, \
                         dh1_37, di_39, di_40, di_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_7 * dh0_37[k]
                  - f_8 * dh1_35[k]
                  + pb_x[k] * di_39[k];

        t_40[k] = f_7 * dh0_39[k]
                  - f_8 * dh1_36[k]
                  + pb_x[k] * di_40[k];

        t_41[k] = f_5 * dh0_40[k]
                  - f_6 * dh1_37[k]
                  + pb_x[k] * di_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, dh0_41, dh0_42, dh0_43, dh1_38, dh1_39, \
                         dh1_40, di_42, di_43, di_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * dh0_41[k]
                  - f_6 * dh1_38[k]
                  + pb_x[k] * di_42[k];

        t_43[k] = f_5 * dh0_42[k]
                  - f_6 * dh1_39[k]
                  + pb_x[k] * di_43[k];

        t_44[k] = f_3 * dh0_43[k]
                  - f_4 * dh1_40[k]
                  + pb_x[k] * di_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pb_x, dh0_44, dh0_45, dh0_47, dh1_41, dh1_42, \
                         dh1_44, di_45, di_46, di_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_3 * dh0_44[k]
                  - f_4 * dh1_41[k]
                  + pb_x[k] * di_45[k];

        t_46[k] = f_3 * dh0_45[k]
                  - f_4 * dh1_42[k]
                  + pb_x[k] * di_46[k];

        t_47[k] = f_3 * dh0_47[k]
                  - f_4 * dh1_44[k]
                  + pb_x[k] * di_47[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_y, dh0_43, dh0_44, dh0_45, dh1_40, dh1_41, \
                         dh1_42, di_48, di_49, di_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * dh0_43[k]
                  - f_2 * dh1_40[k]
                  + pb_y[k] * di_48[k];

        t_49[k] = f_9 * dh0_44[k]
                  - f_10 * dh1_41[k]
                  + pb_y[k] * di_49[k];

        t_50[k] = f_7 * dh0_45[k]
                  - f_8 * dh1_42[k]
                  + pb_y[k] * di_50[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_y, pb_z, pi_11, dh0_46, dh0_47, dh1_43, dh1_44, \
                         di_51, di_52, di_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_5 * dh0_46[k]
                  - f_6 * dh1_43[k]
                  + pb_y[k] * di_51[k];

        t_52[k] = f_3 * dh0_47[k]
                  - f_4 * dh1_44[k]
                  + pb_y[k] * di_52[k];

        t_53[k] = f_0 * pi_11[k]
                  + f_1 * dh0_47[k]
                  - f_2 * dh1_44[k]
                  + pb_z[k] * di_53[k];
    }
}

auto
compute_prim_dk_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pi, const size_t dh0, const size_t dh1,
                                     const size_t di, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_21 = buffer.data(dh1 + 21);
    const auto *dh1_38 = buffer.data(dh1 + 38);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_65 = buffer.data(di + 65);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pi_0, pi_1, dh0_0, dh0_1, dh1_0, dh1_21, di_0, \
                         di_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_0 * pi_1[k]
                 + f_1 * dh0_1[k]
                 - f_2 * dh1_21[k]
                 + pb_y[k] * di_37[k];
    }

#pragma omp simd aligned(t_2, pb_z, pi_2, dh0_2, dh1_38, di_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pi_2[k]
                 + f_1 * dh0_2[k]
                 - f_2 * dh1_38[k]
                 + pb_z[k] * di_65[k];
    }
}

auto
compute_prim_dk_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pi, const size_t pk,
                                     const size_t dh0, const size_t dh1, const size_t di,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_11 = buffer.data(pi + 11);
    const auto *pi_12 = buffer.data(pi + 12);
    const auto *pi_13 = buffer.data(pi + 13);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_16 = buffer.data(pi + 16);
    const auto *pi_17 = buffer.data(pi + 17);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);
    const auto *pk_3 = buffer.data(pk + 3);
    const auto *pk_4 = buffer.data(pk + 4);
    const auto *pk_5 = buffer.data(pk + 5);
    const auto *pk_6 = buffer.data(pk + 6);
    const auto *pk_7 = buffer.data(pk + 7);
    const auto *pk_8 = buffer.data(pk + 8);
    const auto *pk_9 = buffer.data(pk + 9);
    const auto *pk_10 = buffer.data(pk + 10);
    const auto *pk_11 = buffer.data(pk + 11);
    const auto *pk_12 = buffer.data(pk + 12);
    const auto *pk_13 = buffer.data(pk + 13);
    const auto *pk_14 = buffer.data(pk + 14);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);
    const auto *dh0_3 = buffer.data(dh0 + 3);
    const auto *dh0_4 = buffer.data(dh0 + 4);
    const auto *dh0_5 = buffer.data(dh0 + 5);
    const auto *dh0_6 = buffer.data(dh0 + 6);
    const auto *dh0_7 = buffer.data(dh0 + 7);
    const auto *dh0_8 = buffer.data(dh0 + 8);
    const auto *dh0_9 = buffer.data(dh0 + 9);
    const auto *dh0_10 = buffer.data(dh0 + 10);
    const auto *dh0_11 = buffer.data(dh0 + 11);
    const auto *dh0_12 = buffer.data(dh0 + 12);
    const auto *dh0_13 = buffer.data(dh0 + 13);
    const auto *dh0_14 = buffer.data(dh0 + 14);
    const auto *dh0_15 = buffer.data(dh0 + 15);
    const auto *dh0_16 = buffer.data(dh0 + 16);
    const auto *dh0_17 = buffer.data(dh0 + 17);
    const auto *dh0_18 = buffer.data(dh0 + 18);
    const auto *dh0_19 = buffer.data(dh0 + 19);
    const auto *dh0_20 = buffer.data(dh0 + 20);
    const auto *dh0_21 = buffer.data(dh0 + 21);
    const auto *dh0_22 = buffer.data(dh0 + 22);
    const auto *dh0_23 = buffer.data(dh0 + 23);
    const auto *dh0_24 = buffer.data(dh0 + 24);
    const auto *dh0_25 = buffer.data(dh0 + 25);
    const auto *dh0_26 = buffer.data(dh0 + 26);
    const auto *dh0_27 = buffer.data(dh0 + 27);
    const auto *dh0_28 = buffer.data(dh0 + 28);
    const auto *dh0_29 = buffer.data(dh0 + 29);
    const auto *dh0_30 = buffer.data(dh0 + 30);
    const auto *dh0_31 = buffer.data(dh0 + 31);
    const auto *dh0_32 = buffer.data(dh0 + 32);
    const auto *dh0_33 = buffer.data(dh0 + 33);
    const auto *dh0_34 = buffer.data(dh0 + 34);
    const auto *dh0_35 = buffer.data(dh0 + 35);
    const auto *dh0_36 = buffer.data(dh0 + 36);
    const auto *dh0_37 = buffer.data(dh0 + 37);
    const auto *dh0_38 = buffer.data(dh0 + 38);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);
    const auto *dh1_3 = buffer.data(dh1 + 3);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_5 = buffer.data(dh1 + 5);
    const auto *dh1_6 = buffer.data(dh1 + 6);
    const auto *dh1_7 = buffer.data(dh1 + 7);
    const auto *dh1_8 = buffer.data(dh1 + 8);
    const auto *dh1_10 = buffer.data(dh1 + 10);
    const auto *dh1_11 = buffer.data(dh1 + 11);
    const auto *dh1_12 = buffer.data(dh1 + 12);
    const auto *dh1_13 = buffer.data(dh1 + 13);
    const auto *dh1_14 = buffer.data(dh1 + 14);
    const auto *dh1_16 = buffer.data(dh1 + 16);
    const auto *dh1_17 = buffer.data(dh1 + 17);
    const auto *dh1_18 = buffer.data(dh1 + 18);
    const auto *dh1_19 = buffer.data(dh1 + 19);
    const auto *dh1_20 = buffer.data(dh1 + 20);
    const auto *dh1_21 = buffer.data(dh1 + 21);
    const auto *dh1_22 = buffer.data(dh1 + 22);
    const auto *dh1_23 = buffer.data(dh1 + 23);
    const auto *dh1_24 = buffer.data(dh1 + 24);
    const auto *dh1_25 = buffer.data(dh1 + 25);
    const auto *dh1_26 = buffer.data(dh1 + 26);
    const auto *dh1_27 = buffer.data(dh1 + 27);
    const auto *dh1_28 = buffer.data(dh1 + 28);
    const auto *dh1_30 = buffer.data(dh1 + 30);
    const auto *dh1_31 = buffer.data(dh1 + 31);
    const auto *dh1_32 = buffer.data(dh1 + 32);
    const auto *dh1_33 = buffer.data(dh1 + 33);
    const auto *dh1_34 = buffer.data(dh1 + 34);
    const auto *dh1_35 = buffer.data(dh1 + 35);
    const auto *dh1_36 = buffer.data(dh1 + 36);
    const auto *dh1_37 = buffer.data(dh1 + 37);
    const auto *dh1_38 = buffer.data(dh1 + 38);
    const auto *dh1_39 = buffer.data(dh1 + 39);
    const auto *dh1_40 = buffer.data(dh1 + 40);
    const auto *dh1_41 = buffer.data(dh1 + 41);

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
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pi_0, dh0_0, dh1_0, di_0, \
                         di_1, di_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = pb_y[k] * di_0[k];

        t_2[k] = pb_z[k] * di_0[k];

        t_3[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_y[k] * di_1[k];

        t_4[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_z[k] * di_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, dh0_1, dh0_2, dh0_3, dh1_1, \
                         dh1_2, dh1_3, di_3, di_4, di_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * dh0_1[k]
                 - f_6 * dh1_1[k]
                 + pb_y[k] * di_3[k];

        t_6[k] = pb_z[k] * di_3[k];

        t_7[k] = f_5 * dh0_2[k]
                 - f_6 * dh1_2[k]
                 + pb_z[k] * di_4[k];

        t_8[k] = f_7 * dh0_3[k]
                 - f_8 * dh1_3[k]
                 + pb_y[k] * di_5[k];

        t_9[k] = pb_z[k] * di_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pb_z, dh0_4, dh0_5, dh1_4, dh1_5, di_6, \
                         di_7, di_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dh0_4[k]
                  - f_4 * dh1_4[k]
                  + pb_y[k] * di_6[k];

        t_11[k] = f_7 * dh0_4[k]
                  - f_8 * dh1_4[k]
                  + pb_z[k] * di_7[k];

        t_12[k] = f_9 * dh0_5[k]
                  - f_10 * dh1_5[k]
                  + pb_y[k] * di_8[k];

        t_13[k] = pb_z[k] * di_8[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_y, pb_z, dh0_6, dh0_7, dh0_8, dh1_6, \
                         dh1_7, dh1_8, di_9, di_10, di_11, di_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * dh0_6[k]
                  - f_6 * dh1_6[k]
                  + pb_y[k] * di_9[k];

        t_15[k] = f_3 * dh0_7[k]
                  - f_4 * dh1_7[k]
                  + pb_y[k] * di_10[k];

        t_16[k] = f_9 * dh0_7[k]
                  - f_10 * dh1_7[k]
                  + pb_z[k] * di_11[k];

        t_17[k] = f_1 * dh0_8[k]
                  - f_2 * dh1_8[k]
                  + pb_y[k] * di_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_y, pb_z, dh0_9, dh0_10, dh0_11, dh1_10, \
                         dh1_11, dh1_12, di_12, di_14, di_15, di_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_z[k] * di_12[k];

        t_19[k] = f_9 * dh0_9[k]
                  - f_10 * dh1_10[k]
                  + pb_y[k] * di_14[k];

        t_20[k] = f_7 * dh0_10[k]
                  - f_8 * dh1_11[k]
                  + pb_y[k] * di_15[k];

        t_21[k] = f_5 * dh0_11[k]
                  - f_6 * dh1_12[k]
                  + pb_y[k] * di_16[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pa_y, pb_y, pb_z, pi_1, pk_0, pk_1, \
                         dh0_12, dh1_13, di_17, di_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * dh0_12[k]
                  - f_4 * dh1_13[k]
                  + pb_y[k] * di_17[k];

        t_23[k] = f_1 * dh0_12[k]
                  - f_2 * dh1_13[k]
                  + pb_z[k] * di_18[k];

        t_24[k] = pa_y[k] * pk_0[k];

        t_25[k] = f_11 * pi_1[k]
                  + pa_x[k] * pk_1[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_x, pa_z, pi_2, pi_3, pi_4, pk_0, \
                         pk_2, pk_3, pk_4, pk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_12 * pi_2[k]
                  + pa_x[k] * pk_2[k];

        t_27[k] = f_13 * pi_3[k]
                  + pa_x[k] * pk_3[k];

        t_28[k] = f_0 * pi_4[k]
                  + pa_x[k] * pk_4[k];

        t_29[k] = pa_x[k] * pk_5[k];

        t_30[k] = pa_z[k] * pk_0[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pb_z, pi_0, pi_8, pi_9, pi_10, pk_6, \
                         pk_7, pk_8, di_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_14 * pi_0[k]
                  + pb_z[k] * di_19[k];

        t_32[k] = f_11 * pi_8[k]
                  + pa_x[k] * pk_6[k];

        t_33[k] = f_12 * pi_9[k]
                  + pa_x[k] * pk_7[k];

        t_34[k] = f_13 * pi_10[k]
                  + pa_x[k] * pk_8[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_x, pb_z, pi_11, pk_9, pk_14, dh0_13, \
                         dh1_14, di_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * pi_11[k]
                  + pa_x[k] * pk_9[k];

        t_36[k] = pa_x[k] * pk_14[k];

        t_37[k] = f_1 * dh0_13[k]
                  - f_2 * dh1_14[k]
                  + pb_x[k] * di_20[k];

        t_38[k] = pb_z[k] * di_20[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pb_x, pb_z, dh0_14, dh0_15, dh0_16, dh1_16, \
                         dh1_17, dh1_18, di_22, di_23, di_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * dh0_14[k]
                  - f_10 * dh1_16[k]
                  + pb_x[k] * di_22[k];

        t_40[k] = f_9 * dh0_15[k]
                  - f_10 * dh1_17[k]
                  + pb_x[k] * di_23[k];

        t_41[k] = f_7 * dh0_16[k]
                  - f_8 * dh1_18[k]
                  + pb_x[k] * di_24[k];

        t_42[k] = pb_z[k] * di_22[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pb_x, pb_z, dh0_17, dh0_18, dh0_19, dh1_19, \
                         dh1_20, dh1_21, di_24, di_25, di_26, di_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_7 * dh0_17[k]
                  - f_8 * dh1_19[k]
                  + pb_x[k] * di_25[k];

        t_44[k] = f_5 * dh0_18[k]
                  - f_6 * dh1_20[k]
                  + pb_x[k] * di_26[k];

        t_45[k] = pb_z[k] * di_24[k];

        t_46[k] = f_5 * dh0_19[k]
                  - f_6 * dh1_21[k]
                  + pb_x[k] * di_27[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_x, pb_z, dh0_20, dh0_21, dh0_23, dh1_22, \
                         dh1_23, dh1_25, di_26, di_28, di_29, di_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_5 * dh0_20[k]
                  - f_6 * dh1_22[k]
                  + pb_x[k] * di_28[k];

        t_48[k] = f_3 * dh0_21[k]
                  - f_4 * dh1_23[k]
                  + pb_x[k] * di_29[k];

        t_49[k] = pb_z[k] * di_26[k];

        t_50[k] = f_3 * dh0_23[k]
                  - f_4 * dh1_25[k]
                  + pb_x[k] * di_30[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_x, pb_y, pi_5, dh0_21, dh0_24, dh0_25, \
                         dh1_23, dh1_26, dh1_27, di_31, di_32, di_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * dh0_24[k]
                  - f_4 * dh1_26[k]
                  + pb_x[k] * di_31[k];

        t_52[k] = f_3 * dh0_25[k]
                  - f_4 * dh1_27[k]
                  + pb_x[k] * di_32[k];

        t_53[k] = pb_x[k] * di_33[k];

        t_54[k] = f_0 * pi_5[k]
                  + f_1 * dh0_21[k]
                  - f_2 * dh1_23[k]
                  + pb_y[k] * di_33[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_z, dh0_21, dh0_22, dh0_23, dh1_23, dh1_24, \
                         dh1_25, di_33, di_34, di_35, di_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pb_z[k] * di_33[k];

        t_56[k] = f_3 * dh0_21[k]
                  - f_4 * dh1_23[k]
                  + pb_z[k] * di_34[k];

        t_57[k] = f_5 * dh0_22[k]
                  - f_6 * dh1_24[k]
                  + pb_z[k] * di_35[k];

        t_58[k] = f_7 * dh0_23[k]
                  - f_8 * dh1_25[k]
                  + pb_z[k] * di_36[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_z, pb_y, pb_z, pi_6, pk_5, dh0_24, dh0_25, \
                         dh1_26, dh1_27, di_37, di_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_9 * dh0_24[k]
                  - f_10 * dh1_26[k]
                  + pb_z[k] * di_37[k];

        t_60[k] = f_0 * pi_6[k]
                  + pb_y[k] * di_38[k];

        t_61[k] = f_1 * dh0_25[k]
                  - f_2 * dh1_27[k]
                  + pb_z[k] * di_38[k];

        t_62[k] = pa_z[k] * pk_5[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_y, pb_z, pi_5, pi_13, pi_14, pi_15, pk_10, \
                         pk_11, pk_12, di_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_14 * pi_5[k]
                  + pb_z[k] * di_39[k];

        t_64[k] = f_11 * pi_13[k]
                  + pa_y[k] * pk_10[k];

        t_65[k] = f_12 * pi_14[k]
                  + pa_y[k] * pk_11[k];

        t_66[k] = f_13 * pi_15[k]
                  + pa_y[k] * pk_12[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, pa_y, pb_x, pb_y, pi_16, pi_17, pk_13, \
                         pk_14, dh0_26, dh1_28, di_40, di_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_0 * pi_16[k]
                  + pa_y[k] * pk_13[k];

        t_68[k] = f_14 * pi_17[k]
                  + pb_y[k] * di_40[k];

        t_69[k] = pa_y[k] * pk_14[k];

        t_70[k] = f_1 * dh0_26[k]
                  - f_2 * dh1_28[k]
                  + pb_x[k] * di_41[k];

        t_71[k] = pb_y[k] * di_41[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_x, pb_z, pi_7, dh0_27, dh0_28, dh1_30, dh1_31, \
                         di_41, di_43, di_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_0 * pi_7[k]
                  + pb_z[k] * di_41[k];

        t_73[k] = f_9 * dh0_27[k]
                  - f_10 * dh1_30[k]
                  + pb_x[k] * di_43[k];

        t_74[k] = f_9 * dh0_28[k]
                  - f_10 * dh1_31[k]
                  + pb_x[k] * di_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pb_x, pb_y, dh0_29, dh0_30, dh0_31, dh1_32, \
                         dh1_33, dh1_34, di_44, di_45, di_46, di_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_7 * dh0_29[k]
                  - f_8 * dh1_32[k]
                  + pb_x[k] * di_45[k];

        t_76[k] = pb_y[k] * di_44[k];

        t_77[k] = f_7 * dh0_30[k]
                  - f_8 * dh1_33[k]
                  + pb_x[k] * di_46[k];

        t_78[k] = f_5 * dh0_31[k]
                  - f_6 * dh1_34[k]
                  + pb_x[k] * di_47[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_x, pb_y, dh0_32, dh0_33, dh0_34, dh1_35, \
                         dh1_36, dh1_37, di_46, di_48, di_49, di_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_5 * dh0_32[k]
                  - f_6 * dh1_35[k]
                  + pb_x[k] * di_48[k];

        t_80[k] = pb_y[k] * di_46[k];

        t_81[k] = f_5 * dh0_33[k]
                  - f_6 * dh1_36[k]
                  + pb_x[k] * di_49[k];

        t_82[k] = f_3 * dh0_34[k]
                  - f_4 * dh1_37[k]
                  + pb_x[k] * di_50[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pb_x, pb_y, dh0_35, dh0_36, dh0_38, dh1_38, \
                         dh1_39, dh1_41, di_49, di_51, di_52, di_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_3 * dh0_35[k]
                  - f_4 * dh1_38[k]
                  + pb_x[k] * di_51[k];

        t_84[k] = f_3 * dh0_36[k]
                  - f_4 * dh1_39[k]
                  + pb_x[k] * di_52[k];

        t_85[k] = pb_y[k] * di_49[k];

        t_86[k] = f_3 * dh0_38[k]
                  - f_4 * dh1_41[k]
                  + pb_x[k] * di_53[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pb_x, pb_y, pb_z, pi_12, dh0_34, dh0_35, \
                         dh1_37, dh1_38, di_54, di_55, di_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pb_x[k] * di_59[k];

        t_88[k] = f_1 * dh0_34[k]
                  - f_2 * dh1_37[k]
                  + pb_y[k] * di_54[k];

        t_89[k] = f_0 * pi_12[k]
                  + pb_z[k] * di_54[k];

        t_90[k] = f_9 * dh0_35[k]
                  - f_10 * dh1_38[k]
                  + pb_y[k] * di_55[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pb_y, dh0_36, dh0_37, dh0_38, dh1_39, dh1_40, \
                         dh1_41, di_56, di_57, di_58, di_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_7 * dh0_36[k]
                  - f_8 * dh1_39[k]
                  + pb_y[k] * di_56[k];

        t_92[k] = f_5 * dh0_37[k]
                  - f_6 * dh1_40[k]
                  + pb_y[k] * di_57[k];

        t_93[k] = f_3 * dh0_38[k]
                  - f_4 * dh1_41[k]
                  + pb_y[k] * di_58[k];

        t_94[k] = pb_y[k] * di_59[k];
    }

#pragma omp simd aligned(t_95, pb_z, pi_17, dh0_38, dh1_41, di_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_0 * pi_17[k]
                  + f_1 * dh0_38[k]
                  - f_2 * dh1_41[k]
                  + pb_z[k] * di_59[k];
    }
}

auto
compute_prim_dk_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pi, const size_t dh0, const size_t dh1,
                                      const size_t di, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);
    const auto *dh0_3 = buffer.data(dh0 + 3);
    const auto *dh0_4 = buffer.data(dh0 + 4);
    const auto *dh0_5 = buffer.data(dh0 + 5);
    const auto *dh0_6 = buffer.data(dh0 + 6);
    const auto *dh0_7 = buffer.data(dh0 + 7);
    const auto *dh0_8 = buffer.data(dh0 + 8);
    const auto *dh0_10 = buffer.data(dh0 + 10);
    const auto *dh0_11 = buffer.data(dh0 + 11);
    const auto *dh0_12 = buffer.data(dh0 + 12);
    const auto *dh0_13 = buffer.data(dh0 + 13);
    const auto *dh0_14 = buffer.data(dh0 + 14);
    const auto *dh0_16 = buffer.data(dh0 + 16);
    const auto *dh0_17 = buffer.data(dh0 + 17);
    const auto *dh0_18 = buffer.data(dh0 + 18);
    const auto *dh0_19 = buffer.data(dh0 + 19);
    const auto *dh0_20 = buffer.data(dh0 + 20);
    const auto *dh0_21 = buffer.data(dh0 + 21);
    const auto *dh0_22 = buffer.data(dh0 + 22);
    const auto *dh0_23 = buffer.data(dh0 + 23);
    const auto *dh0_24 = buffer.data(dh0 + 24);
    const auto *dh0_25 = buffer.data(dh0 + 25);
    const auto *dh0_26 = buffer.data(dh0 + 26);
    const auto *dh0_27 = buffer.data(dh0 + 27);
    const auto *dh0_28 = buffer.data(dh0 + 28);
    const auto *dh0_30 = buffer.data(dh0 + 30);
    const auto *dh0_31 = buffer.data(dh0 + 31);
    const auto *dh0_32 = buffer.data(dh0 + 32);
    const auto *dh0_33 = buffer.data(dh0 + 33);
    const auto *dh0_34 = buffer.data(dh0 + 34);
    const auto *dh0_35 = buffer.data(dh0 + 35);
    const auto *dh0_36 = buffer.data(dh0 + 36);
    const auto *dh0_37 = buffer.data(dh0 + 37);
    const auto *dh0_38 = buffer.data(dh0 + 38);
    const auto *dh0_39 = buffer.data(dh0 + 39);
    const auto *dh0_40 = buffer.data(dh0 + 40);
    const auto *dh0_41 = buffer.data(dh0 + 41);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);
    const auto *dh1_3 = buffer.data(dh1 + 3);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_5 = buffer.data(dh1 + 5);
    const auto *dh1_6 = buffer.data(dh1 + 6);
    const auto *dh1_7 = buffer.data(dh1 + 7);
    const auto *dh1_8 = buffer.data(dh1 + 8);
    const auto *dh1_9 = buffer.data(dh1 + 9);
    const auto *dh1_10 = buffer.data(dh1 + 10);
    const auto *dh1_11 = buffer.data(dh1 + 11);
    const auto *dh1_12 = buffer.data(dh1 + 12);
    const auto *dh1_13 = buffer.data(dh1 + 13);
    const auto *dh1_14 = buffer.data(dh1 + 14);
    const auto *dh1_15 = buffer.data(dh1 + 15);
    const auto *dh1_16 = buffer.data(dh1 + 16);
    const auto *dh1_17 = buffer.data(dh1 + 17);
    const auto *dh1_18 = buffer.data(dh1 + 18);
    const auto *dh1_19 = buffer.data(dh1 + 19);
    const auto *dh1_20 = buffer.data(dh1 + 20);
    const auto *dh1_21 = buffer.data(dh1 + 21);
    const auto *dh1_22 = buffer.data(dh1 + 22);
    const auto *dh1_23 = buffer.data(dh1 + 23);
    const auto *dh1_24 = buffer.data(dh1 + 24);
    const auto *dh1_25 = buffer.data(dh1 + 25);
    const auto *dh1_26 = buffer.data(dh1 + 26);
    const auto *dh1_27 = buffer.data(dh1 + 27);
    const auto *dh1_28 = buffer.data(dh1 + 28);
    const auto *dh1_29 = buffer.data(dh1 + 29);
    const auto *dh1_30 = buffer.data(dh1 + 30);
    const auto *dh1_31 = buffer.data(dh1 + 31);
    const auto *dh1_32 = buffer.data(dh1 + 32);
    const auto *dh1_33 = buffer.data(dh1 + 33);
    const auto *dh1_34 = buffer.data(dh1 + 34);
    const auto *dh1_35 = buffer.data(dh1 + 35);
    const auto *dh1_36 = buffer.data(dh1 + 36);
    const auto *dh1_37 = buffer.data(dh1 + 37);
    const auto *dh1_38 = buffer.data(dh1 + 38);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pi_0, dh0_0, dh1_0, di_0, \
                         di_1, di_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = pb_y[k] * di_0[k];

        t_2[k] = pb_z[k] * di_0[k];

        t_3[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_y[k] * di_1[k];

        t_4[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_z[k] * di_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, dh0_1, dh0_2, dh0_3, dh1_1, dh1_2, \
                         dh1_3, di_3, di_4, di_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * dh0_1[k]
                 - f_6 * dh1_1[k]
                 + pb_y[k] * di_3[k];

        t_6[k] = pb_y[k] * di_4[k];

        t_7[k] = f_5 * dh0_2[k]
                 - f_6 * dh1_2[k]
                 + pb_z[k] * di_4[k];

        t_8[k] = f_7 * dh0_3[k]
                 - f_8 * dh1_3[k]
                 + pb_y[k] * di_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, dh0_4, dh0_5, dh1_4, dh1_5, di_6, \
                         di_7, di_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * dh0_4[k]
                 - f_4 * dh1_4[k]
                 + pb_y[k] * di_6[k];

        t_10[k] = pb_y[k] * di_7[k];

        t_11[k] = f_7 * dh0_4[k]
                  - f_8 * dh1_4[k]
                  + pb_z[k] * di_7[k];

        t_12[k] = f_9 * dh0_5[k]
                  - f_10 * dh1_5[k]
                  + pb_y[k] * di_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pb_y, pb_z, dh0_6, dh0_7, dh1_6, dh1_7, di_9, \
                         di_10, di_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * dh0_6[k]
                  - f_6 * dh1_6[k]
                  + pb_y[k] * di_9[k];

        t_14[k] = f_3 * dh0_7[k]
                  - f_4 * dh1_7[k]
                  + pb_y[k] * di_10[k];

        t_15[k] = pb_y[k] * di_11[k];

        t_16[k] = f_9 * dh0_7[k]
                  - f_10 * dh1_7[k]
                  + pb_z[k] * di_11[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_y, dh0_8, dh0_10, dh0_11, dh1_8, dh1_9, dh1_10, \
                         di_12, di_13, di_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * dh0_8[k]
                  - f_2 * dh1_8[k]
                  + pb_y[k] * di_12[k];

        t_18[k] = f_9 * dh0_10[k]
                  - f_10 * dh1_9[k]
                  + pb_y[k] * di_13[k];

        t_19[k] = f_7 * dh0_11[k]
                  - f_8 * dh1_10[k]
                  + pb_y[k] * di_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pb_y, pb_z, dh0_12, dh0_13, dh1_11, dh1_12, \
                         di_15, di_16, di_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * dh0_12[k]
                  - f_6 * dh1_11[k]
                  + pb_y[k] * di_15[k];

        t_21[k] = f_3 * dh0_13[k]
                  - f_4 * dh1_12[k]
                  + pb_y[k] * di_16[k];

        t_22[k] = pb_y[k] * di_17[k];

        t_23[k] = f_1 * dh0_13[k]
                  - f_2 * dh1_12[k]
                  + pb_z[k] * di_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_x, dh0_14, dh0_16, dh0_17, dh1_13, dh1_14, \
                         dh1_15, di_18, di_19, di_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_1 * dh0_14[k]
                  - f_2 * dh1_13[k]
                  + pb_x[k] * di_18[k];

        t_25[k] = f_9 * dh0_16[k]
                  - f_10 * dh1_14[k]
                  + pb_x[k] * di_19[k];

        t_26[k] = f_9 * dh0_17[k]
                  - f_10 * dh1_15[k]
                  + pb_x[k] * di_20[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_x, dh0_18, dh0_19, dh0_20, dh1_16, dh1_17, \
                         dh1_18, di_21, di_22, di_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * dh0_18[k]
                  - f_8 * dh1_16[k]
                  + pb_x[k] * di_21[k];

        t_28[k] = f_7 * dh0_19[k]
                  - f_8 * dh1_17[k]
                  + pb_x[k] * di_22[k];

        t_29[k] = f_5 * dh0_20[k]
                  - f_6 * dh1_18[k]
                  + pb_x[k] * di_23[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_x, dh0_21, dh0_22, dh0_23, dh1_19, dh1_20, \
                         dh1_21, di_24, di_25, di_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * dh0_21[k]
                  - f_6 * dh1_19[k]
                  + pb_x[k] * di_24[k];

        t_31[k] = f_5 * dh0_22[k]
                  - f_6 * dh1_20[k]
                  + pb_x[k] * di_25[k];

        t_32[k] = f_3 * dh0_23[k]
                  - f_4 * dh1_21[k]
                  + pb_x[k] * di_26[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_x, dh0_25, dh0_26, dh0_27, dh1_23, dh1_24, \
                         dh1_25, di_27, di_28, di_29, di_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * dh0_25[k]
                  - f_4 * dh1_23[k]
                  + pb_x[k] * di_27[k];

        t_34[k] = f_3 * dh0_26[k]
                  - f_4 * dh1_24[k]
                  + pb_x[k] * di_28[k];

        t_35[k] = f_3 * dh0_27[k]
                  - f_4 * dh1_25[k]
                  + pb_x[k] * di_29[k];

        t_36[k] = pb_x[k] * di_30[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pb_x, pb_y, pi_1, dh0_23, dh1_21, \
                         di_30, di_32, di_33, di_34, di_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_x[k] * di_32[k];

        t_38[k] = pb_x[k] * di_33[k];

        t_39[k] = pb_x[k] * di_34[k];

        t_40[k] = pb_x[k] * di_35[k];

        t_41[k] = f_0 * pi_1[k]
                  + f_1 * dh0_23[k]
                  - f_2 * dh1_21[k]
                  + pb_y[k] * di_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_z, dh0_23, dh0_24, dh0_25, dh1_21, dh1_22, \
                         dh1_23, di_30, di_31, di_32, di_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * di_30[k];

        t_43[k] = f_3 * dh0_23[k]
                  - f_4 * dh1_21[k]
                  + pb_z[k] * di_31[k];

        t_44[k] = f_5 * dh0_24[k]
                  - f_6 * dh1_22[k]
                  + pb_z[k] * di_32[k];

        t_45[k] = f_7 * dh0_25[k]
                  - f_8 * dh1_23[k]
                  + pb_z[k] * di_33[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_x, pb_z, dh0_26, dh0_27, dh0_28, dh1_24, dh1_25, \
                         dh1_26, di_34, di_35, di_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_9 * dh0_26[k]
                  - f_10 * dh1_24[k]
                  + pb_z[k] * di_34[k];

        t_47[k] = f_1 * dh0_27[k]
                  - f_2 * dh1_25[k]
                  + pb_z[k] * di_35[k];

        t_48[k] = f_1 * dh0_28[k]
                  - f_2 * dh1_26[k]
                  + pb_x[k] * di_36[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_x, dh0_30, dh0_31, dh0_32, dh1_27, dh1_28, \
                         dh1_29, di_37, di_38, di_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_9 * dh0_30[k]
                  - f_10 * dh1_27[k]
                  + pb_x[k] * di_37[k];

        t_50[k] = f_9 * dh0_31[k]
                  - f_10 * dh1_28[k]
                  + pb_x[k] * di_38[k];

        t_51[k] = f_7 * dh0_32[k]
                  - f_8 * dh1_29[k]
                  + pb_x[k] * di_39[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_x, dh0_33, dh0_34, dh0_35, dh1_30, dh1_31, \
                         dh1_32, di_40, di_41, di_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_7 * dh0_33[k]
                  - f_8 * dh1_30[k]
                  + pb_x[k] * di_40[k];

        t_53[k] = f_5 * dh0_34[k]
                  - f_6 * dh1_31[k]
                  + pb_x[k] * di_41[k];

        t_54[k] = f_5 * dh0_35[k]
                  - f_6 * dh1_32[k]
                  + pb_x[k] * di_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pb_x, dh0_36, dh0_37, dh0_38, dh1_33, dh1_34, \
                         dh1_35, di_43, di_44, di_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_5 * dh0_36[k]
                  - f_6 * dh1_33[k]
                  + pb_x[k] * di_43[k];

        t_56[k] = f_3 * dh0_37[k]
                  - f_4 * dh1_34[k]
                  + pb_x[k] * di_44[k];

        t_57[k] = f_3 * dh0_38[k]
                  - f_4 * dh1_35[k]
                  + pb_x[k] * di_45[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_x, dh0_39, dh0_41, dh1_36, dh1_38, \
                         di_46, di_47, di_48, di_49, di_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_3 * dh0_39[k]
                  - f_4 * dh1_36[k]
                  + pb_x[k] * di_46[k];

        t_59[k] = f_3 * dh0_41[k]
                  - f_4 * dh1_38[k]
                  + pb_x[k] * di_47[k];

        t_60[k] = pb_x[k] * di_48[k];

        t_61[k] = pb_x[k] * di_49[k];

        t_62[k] = pb_x[k] * di_50[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pb_x, pb_y, dh0_37, dh0_38, dh1_34, dh1_35, \
                         di_48, di_49, di_51, di_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_x[k] * di_51[k];

        t_64[k] = pb_x[k] * di_53[k];

        t_65[k] = f_1 * dh0_37[k]
                  - f_2 * dh1_34[k]
                  + pb_y[k] * di_48[k];

        t_66[k] = f_9 * dh0_38[k]
                  - f_10 * dh1_35[k]
                  + pb_y[k] * di_49[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pb_y, dh0_39, dh0_40, dh0_41, dh1_36, dh1_37, \
                         dh1_38, di_50, di_51, di_52, di_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_7 * dh0_39[k]
                  - f_8 * dh1_36[k]
                  + pb_y[k] * di_50[k];

        t_68[k] = f_5 * dh0_40[k]
                  - f_6 * dh1_37[k]
                  + pb_y[k] * di_51[k];

        t_69[k] = f_3 * dh0_41[k]
                  - f_4 * dh1_38[k]
                  + pb_y[k] * di_52[k];

        t_70[k] = pb_y[k] * di_53[k];
    }

#pragma omp simd aligned(t_71, pb_z, pi_2, dh0_41, dh1_38, di_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_0 * pi_2[k]
                  + f_1 * dh0_41[k]
                  - f_2 * dh1_38[k]
                  + pb_z[k] * di_53[k];
    }
}

auto
compute_prim_dk_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pi, const size_t dh0, const size_t dh1,
                                      const size_t di, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_20 = buffer.data(dh1 + 20);
    const auto *dh1_38 = buffer.data(dh1 + 38);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_35 = buffer.data(di + 35);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pi_0, pi_1, dh0_0, dh0_1, dh1_0, dh1_20, di_0, \
                         di_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_0 * pi_1[k]
                 + f_1 * dh0_1[k]
                 - f_2 * dh1_20[k]
                 + pb_y[k] * di_19[k];
    }

#pragma omp simd aligned(t_2, pb_z, pi_2, dh0_2, dh1_38, di_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pi_2[k]
                 + f_1 * dh0_2[k]
                 - f_2 * dh1_38[k]
                 + pb_z[k] * di_35[k];
    }
}

auto
compute_prim_dk_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pi, const size_t dh0, const size_t dh1,
                                      const size_t di, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_20 = buffer.data(dh0 + 20);
    const auto *dh0_38 = buffer.data(dh0 + 38);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_19 = buffer.data(dh1 + 19);
    const auto *dh1_35 = buffer.data(dh1 + 35);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pi_0, pi_1, dh0_0, dh0_20, dh1_0, dh1_19, di_0, \
                         di_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_0 * pi_1[k]
                 + f_1 * dh0_20[k]
                 - f_2 * dh1_19[k]
                 + pb_y[k] * di_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pi_2, dh0_38, dh1_35, di_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pi_2[k]
                 + f_1 * dh0_38[k]
                 - f_2 * dh1_35[k]
                 + pb_z[k] * di_2[k];
    }
}

auto
compute_prim_dk_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pi, const size_t dh0, const size_t dh1,
                                      const size_t di, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_21 = buffer.data(dh1 + 21);
    const auto *dh1_38 = buffer.data(dh1 + 38);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_25 = buffer.data(di + 25);
    const auto *di_47 = buffer.data(di + 47);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pi_0, pi_1, dh0_0, dh0_1, dh1_0, dh1_21, di_0, \
                         di_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_0 * pi_1[k]
                 + f_1 * dh0_1[k]
                 - f_2 * dh1_21[k]
                 + pb_y[k] * di_25[k];
    }

#pragma omp simd aligned(t_2, pb_z, pi_2, dh0_2, dh1_38, di_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pi_2[k]
                 + f_1 * dh0_2[k]
                 - f_2 * dh1_38[k]
                 + pb_z[k] * di_47[k];
    }
}

auto
compute_prim_dk_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pi, const size_t dh0, const size_t dh1,
                                      const size_t di, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);
    const auto *dh0_3 = buffer.data(dh0 + 3);
    const auto *dh0_4 = buffer.data(dh0 + 4);
    const auto *dh0_5 = buffer.data(dh0 + 5);
    const auto *dh0_6 = buffer.data(dh0 + 6);
    const auto *dh0_7 = buffer.data(dh0 + 7);
    const auto *dh0_9 = buffer.data(dh0 + 9);
    const auto *dh0_10 = buffer.data(dh0 + 10);
    const auto *dh0_11 = buffer.data(dh0 + 11);
    const auto *dh0_12 = buffer.data(dh0 + 12);
    const auto *dh0_13 = buffer.data(dh0 + 13);
    const auto *dh0_15 = buffer.data(dh0 + 15);
    const auto *dh0_16 = buffer.data(dh0 + 16);
    const auto *dh0_17 = buffer.data(dh0 + 17);
    const auto *dh0_18 = buffer.data(dh0 + 18);
    const auto *dh0_19 = buffer.data(dh0 + 19);
    const auto *dh0_20 = buffer.data(dh0 + 20);
    const auto *dh0_21 = buffer.data(dh0 + 21);
    const auto *dh0_22 = buffer.data(dh0 + 22);
    const auto *dh0_23 = buffer.data(dh0 + 23);
    const auto *dh0_24 = buffer.data(dh0 + 24);
    const auto *dh0_25 = buffer.data(dh0 + 25);
    const auto *dh0_26 = buffer.data(dh0 + 26);
    const auto *dh0_28 = buffer.data(dh0 + 28);
    const auto *dh0_29 = buffer.data(dh0 + 29);
    const auto *dh0_30 = buffer.data(dh0 + 30);
    const auto *dh0_31 = buffer.data(dh0 + 31);
    const auto *dh0_32 = buffer.data(dh0 + 32);
    const auto *dh0_33 = buffer.data(dh0 + 33);
    const auto *dh0_34 = buffer.data(dh0 + 34);
    const auto *dh0_35 = buffer.data(dh0 + 35);
    const auto *dh0_36 = buffer.data(dh0 + 36);
    const auto *dh0_37 = buffer.data(dh0 + 37);
    const auto *dh0_38 = buffer.data(dh0 + 38);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);
    const auto *dh1_3 = buffer.data(dh1 + 3);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_5 = buffer.data(dh1 + 5);
    const auto *dh1_6 = buffer.data(dh1 + 6);
    const auto *dh1_7 = buffer.data(dh1 + 7);
    const auto *dh1_8 = buffer.data(dh1 + 8);
    const auto *dh1_9 = buffer.data(dh1 + 9);
    const auto *dh1_10 = buffer.data(dh1 + 10);
    const auto *dh1_11 = buffer.data(dh1 + 11);
    const auto *dh1_12 = buffer.data(dh1 + 12);
    const auto *dh1_13 = buffer.data(dh1 + 13);
    const auto *dh1_14 = buffer.data(dh1 + 14);
    const auto *dh1_15 = buffer.data(dh1 + 15);
    const auto *dh1_16 = buffer.data(dh1 + 16);
    const auto *dh1_17 = buffer.data(dh1 + 17);
    const auto *dh1_18 = buffer.data(dh1 + 18);
    const auto *dh1_19 = buffer.data(dh1 + 19);
    const auto *dh1_20 = buffer.data(dh1 + 20);
    const auto *dh1_21 = buffer.data(dh1 + 21);
    const auto *dh1_22 = buffer.data(dh1 + 22);
    const auto *dh1_23 = buffer.data(dh1 + 23);
    const auto *dh1_24 = buffer.data(dh1 + 24);
    const auto *dh1_25 = buffer.data(dh1 + 25);
    const auto *dh1_26 = buffer.data(dh1 + 26);
    const auto *dh1_27 = buffer.data(dh1 + 27);
    const auto *dh1_28 = buffer.data(dh1 + 28);
    const auto *dh1_29 = buffer.data(dh1 + 29);
    const auto *dh1_30 = buffer.data(dh1 + 30);
    const auto *dh1_31 = buffer.data(dh1 + 31);
    const auto *dh1_32 = buffer.data(dh1 + 32);
    const auto *dh1_33 = buffer.data(dh1 + 33);
    const auto *dh1_34 = buffer.data(dh1 + 34);
    const auto *dh1_35 = buffer.data(dh1 + 35);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pi_0, dh0_0, dh0_1, dh1_0, \
                         dh1_1, di_0, di_1, di_2, di_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_y[k] * di_1[k];

        t_2[k] = f_3 * dh0_0[k]
                 - f_4 * dh1_0[k]
                 + pb_z[k] * di_2[k];

        t_3[k] = f_5 * dh0_1[k]
                 - f_6 * dh1_1[k]
                 + pb_y[k] * di_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, dh0_2, dh0_3, dh0_4, dh1_2, dh1_3, dh1_4, \
                         di_4, di_5, di_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * dh0_2[k]
                 - f_6 * dh1_2[k]
                 + pb_z[k] * di_4[k];

        t_5[k] = f_7 * dh0_3[k]
                 - f_8 * dh1_3[k]
                 + pb_y[k] * di_5[k];

        t_6[k] = f_7 * dh0_4[k]
                 - f_8 * dh1_4[k]
                 + pb_z[k] * di_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_y, pb_z, dh0_5, dh0_6, dh0_7, dh1_5, dh1_6, dh1_7, \
                         di_7, di_8, di_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_9 * dh0_5[k]
                 - f_10 * dh1_5[k]
                 + pb_y[k] * di_7[k];

        t_8[k] = f_9 * dh0_6[k]
                 - f_10 * dh1_6[k]
                 + pb_z[k] * di_8[k];

        t_9[k] = f_1 * dh0_7[k]
                 - f_2 * dh1_7[k]
                 + pb_y[k] * di_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, dh0_9, dh0_10, dh0_11, dh1_8, dh1_9, dh1_10, \
                         di_10, di_11, di_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_9 * dh0_9[k]
                  - f_10 * dh1_8[k]
                  + pb_y[k] * di_10[k];

        t_11[k] = f_7 * dh0_10[k]
                  - f_8 * dh1_9[k]
                  + pb_y[k] * di_11[k];

        t_12[k] = f_5 * dh0_11[k]
                  - f_6 * dh1_10[k]
                  + pb_y[k] * di_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, pb_y, pb_z, dh0_12, dh0_13, dh1_11, dh1_12, \
                         di_13, di_14, di_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * dh0_12[k]
                  - f_4 * dh1_11[k]
                  + pb_y[k] * di_13[k];

        t_14[k] = f_1 * dh0_12[k]
                  - f_2 * dh1_11[k]
                  + pb_z[k] * di_14[k];

        t_15[k] = f_1 * dh0_13[k]
                  - f_2 * dh1_12[k]
                  + pb_x[k] * di_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, dh0_15, dh0_16, dh0_17, dh1_13, dh1_14, \
                         dh1_15, di_16, di_17, di_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_9 * dh0_15[k]
                  - f_10 * dh1_13[k]
                  + pb_x[k] * di_16[k];

        t_17[k] = f_9 * dh0_16[k]
                  - f_10 * dh1_14[k]
                  + pb_x[k] * di_17[k];

        t_18[k] = f_7 * dh0_17[k]
                  - f_8 * dh1_15[k]
                  + pb_x[k] * di_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_x, dh0_18, dh0_19, dh0_20, dh1_16, dh1_17, \
                         dh1_18, di_19, di_20, di_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * dh0_18[k]
                  - f_8 * dh1_16[k]
                  + pb_x[k] * di_19[k];

        t_20[k] = f_5 * dh0_19[k]
                  - f_6 * dh1_17[k]
                  + pb_x[k] * di_20[k];

        t_21[k] = f_5 * dh0_20[k]
                  - f_6 * dh1_18[k]
                  + pb_x[k] * di_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pb_x, pb_y, pb_z, pi_1, dh0_21, dh0_25, \
                         dh1_19, dh1_23, di_22, di_23, di_24, di_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * dh0_21[k]
                  - f_4 * dh1_19[k]
                  + pb_x[k] * di_22[k];

        t_23[k] = f_3 * dh0_25[k]
                  - f_4 * dh1_23[k]
                  + pb_x[k] * di_23[k];

        t_24[k] = f_0 * pi_1[k]
                  + f_1 * dh0_21[k]
                  - f_2 * dh1_19[k]
                  + pb_y[k] * di_24[k];

        t_25[k] = f_3 * dh0_21[k]
                  - f_4 * dh1_19[k]
                  + pb_z[k] * di_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pb_z, dh0_22, dh0_23, dh0_24, dh1_20, dh1_21, \
                         dh1_22, di_26, di_27, di_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * dh0_22[k]
                  - f_6 * dh1_20[k]
                  + pb_z[k] * di_26[k];

        t_27[k] = f_7 * dh0_23[k]
                  - f_8 * dh1_21[k]
                  + pb_z[k] * di_27[k];

        t_28[k] = f_9 * dh0_24[k]
                  - f_10 * dh1_22[k]
                  + pb_z[k] * di_28[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, pb_z, dh0_25, dh0_26, dh0_28, dh1_23, dh1_24, \
                         dh1_25, di_29, di_30, di_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * dh0_25[k]
                  - f_2 * dh1_23[k]
                  + pb_z[k] * di_29[k];

        t_30[k] = f_1 * dh0_26[k]
                  - f_2 * dh1_24[k]
                  + pb_x[k] * di_30[k];

        t_31[k] = f_9 * dh0_28[k]
                  - f_10 * dh1_25[k]
                  + pb_x[k] * di_31[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, dh0_29, dh0_30, dh0_31, dh1_26, dh1_27, \
                         dh1_28, di_32, di_33, di_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_9 * dh0_29[k]
                  - f_10 * dh1_26[k]
                  + pb_x[k] * di_32[k];

        t_33[k] = f_7 * dh0_30[k]
                  - f_8 * dh1_27[k]
                  + pb_x[k] * di_33[k];

        t_34[k] = f_7 * dh0_31[k]
                  - f_8 * dh1_28[k]
                  + pb_x[k] * di_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_x, dh0_32, dh0_33, dh0_34, dh1_29, dh1_30, \
                         dh1_31, di_35, di_36, di_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_5 * dh0_32[k]
                  - f_6 * dh1_29[k]
                  + pb_x[k] * di_35[k];

        t_36[k] = f_5 * dh0_33[k]
                  - f_6 * dh1_30[k]
                  + pb_x[k] * di_36[k];

        t_37[k] = f_3 * dh0_34[k]
                  - f_4 * dh1_31[k]
                  + pb_x[k] * di_37[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pb_x, pb_y, dh0_34, dh0_35, dh0_38, dh1_31, dh1_32, \
                         dh1_35, di_38, di_39, di_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * dh0_38[k]
                  - f_4 * dh1_35[k]
                  + pb_x[k] * di_38[k];

        t_39[k] = f_1 * dh0_34[k]
                  - f_2 * dh1_31[k]
                  + pb_y[k] * di_39[k];

        t_40[k] = f_9 * dh0_35[k]
                  - f_10 * dh1_32[k]
                  + pb_y[k] * di_40[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pb_y, dh0_36, dh0_37, dh0_38, dh1_33, dh1_34, \
                         dh1_35, di_41, di_42, di_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * dh0_36[k]
                  - f_8 * dh1_33[k]
                  + pb_y[k] * di_41[k];

        t_42[k] = f_5 * dh0_37[k]
                  - f_6 * dh1_34[k]
                  + pb_y[k] * di_42[k];

        t_43[k] = f_3 * dh0_38[k]
                  - f_4 * dh1_35[k]
                  + pb_y[k] * di_43[k];
    }

#pragma omp simd aligned(t_44, pb_z, pi_2, dh0_38, dh1_35, di_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * pi_2[k]
                  + f_1 * dh0_38[k]
                  - f_2 * dh1_35[k]
                  + pb_z[k] * di_44[k];
    }
}

auto
compute_prim_dk_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pi, const size_t dh0, const size_t dh1,
                                      const size_t di, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_19 = buffer.data(dh0 + 19);
    const auto *dh0_35 = buffer.data(dh0 + 35);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_19 = buffer.data(dh1 + 19);
    const auto *dh1_35 = buffer.data(dh1 + 35);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pi_0, pi_1, dh0_0, dh0_19, dh1_0, dh1_19, di_0, \
                         di_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_0 * pi_1[k]
                 + f_1 * dh0_19[k]
                 - f_2 * dh1_19[k]
                 + pb_y[k] * di_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pi_2, dh0_35, dh1_35, di_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pi_2[k]
                 + f_1 * dh0_35[k]
                 - f_2 * dh1_35[k]
                 + pb_z[k] * di_2[k];
    }
}

auto
compute_prim_dk_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pi, const size_t dh0, const size_t dh1,
                                      const size_t di, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_11 = buffer.data(dh1 + 11);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_11 = buffer.data(di + 11);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pi_0, pi_1, dh0_0, dh0_1, dh1_0, dh1_4, di_0, \
                         di_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_0 * pi_1[k]
                 + f_1 * dh0_1[k]
                 - f_2 * dh1_4[k]
                 + pb_y[k] * di_4[k];
    }

#pragma omp simd aligned(t_2, pb_z, pi_2, dh0_2, dh1_11, di_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pi_2[k]
                 + f_1 * dh0_2[k]
                 - f_2 * dh1_11[k]
                 + pb_z[k] * di_11[k];
    }
}

auto
compute_prim_dk_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pi, const size_t dh0, const size_t dh1,
                                      const size_t di, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_4 = buffer.data(dh0 + 4);
    const auto *dh0_11 = buffer.data(dh0 + 11);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_11 = buffer.data(dh1 + 11);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_11 = buffer.data(di + 11);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pi_0, pi_1, dh0_0, dh0_4, dh1_0, dh1_4, di_0, \
                         di_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_0 * pi_1[k]
                 + f_1 * dh0_4[k]
                 - f_2 * dh1_4[k]
                 + pb_y[k] * di_4[k];
    }

#pragma omp simd aligned(t_2, pb_z, pi_2, dh0_11, dh1_11, di_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pi_2[k]
                 + f_1 * dh0_11[k]
                 - f_2 * dh1_11[k]
                 + pb_z[k] * di_11[k];
    }
}

auto
compute_prim_dk_electron_repulsion_18(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pi, const size_t dh0, const size_t dh1,
                                      const size_t di, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_4 = buffer.data(dh0 + 4);
    const auto *dh0_11 = buffer.data(dh0 + 11);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_7 = buffer.data(dh1 + 7);
    const auto *dh1_17 = buffer.data(dh1 + 17);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_11 = buffer.data(di + 11);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pi_0, pi_1, dh0_0, dh0_4, dh1_0, dh1_7, di_0, \
                         di_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_0 * pi_1[k]
                 + f_1 * dh0_4[k]
                 - f_2 * dh1_7[k]
                 + pb_y[k] * di_4[k];
    }

#pragma omp simd aligned(t_2, pb_z, pi_2, dh0_11, dh1_17, di_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pi_2[k]
                 + f_1 * dh0_11[k]
                 - f_2 * dh1_17[k]
                 + pb_z[k] * di_11[k];
    }
}

auto
compute_prim_dk_electron_repulsion_19(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pi, const size_t dh0, const size_t dh1,
                                      const size_t di, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_7 = buffer.data(dh0 + 7);
    const auto *dh0_17 = buffer.data(dh0 + 17);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_11 = buffer.data(dh1 + 11);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pi_0, pi_1, dh0_0, dh0_7, dh1_0, dh1_4, di_0, \
                         di_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_0 * pi_1[k]
                 + f_1 * dh0_7[k]
                 - f_2 * dh1_4[k]
                 + pb_y[k] * di_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pi_2, dh0_17, dh1_11, di_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pi_2[k]
                 + f_1 * dh0_17[k]
                 - f_2 * dh1_11[k]
                 + pb_z[k] * di_2[k];
    }
}

auto
compute_prim_dk_electron_repulsion_20(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pi, const size_t dh0, const size_t dh1,
                                      const size_t di, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_11 = buffer.data(dh1 + 11);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_14 = buffer.data(di + 14);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pi_0, pi_1, dh0_0, dh0_1, dh1_0, dh1_4, di_0, \
                         di_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_0 * pi_1[k]
                 + f_1 * dh0_1[k]
                 - f_2 * dh1_4[k]
                 + pb_y[k] * di_5[k];
    }

#pragma omp simd aligned(t_2, pb_z, pi_2, dh0_2, dh1_11, di_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pi_2[k]
                 + f_1 * dh0_2[k]
                 - f_2 * dh1_11[k]
                 + pb_z[k] * di_14[k];
    }
}

auto
compute_prim_dk_electron_repulsion_21(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pi, const size_t dh0, const size_t dh1,
                                      const size_t di, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 2.0 / beta;
    const auto f_4 = 2.0 * alpha / (beta * p);
    const auto f_5 = 1.5 / beta;
    const auto f_6 = 1.5 * alpha / (beta * p);
    const auto f_7 = 1.0 / beta;
    const auto f_8 = alpha / (beta * p);
    const auto f_9 = 0.5 / beta;
    const auto f_10 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);
    const auto *dh0_3 = buffer.data(dh0 + 3);
    const auto *dh0_4 = buffer.data(dh0 + 4);
    const auto *dh0_5 = buffer.data(dh0 + 5);
    const auto *dh0_6 = buffer.data(dh0 + 6);
    const auto *dh0_7 = buffer.data(dh0 + 7);
    const auto *dh0_8 = buffer.data(dh0 + 8);
    const auto *dh0_9 = buffer.data(dh0 + 9);
    const auto *dh0_10 = buffer.data(dh0 + 10);
    const auto *dh0_11 = buffer.data(dh0 + 11);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);
    const auto *dh1_3 = buffer.data(dh1 + 3);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_5 = buffer.data(dh1 + 5);
    const auto *dh1_6 = buffer.data(dh1 + 6);
    const auto *dh1_7 = buffer.data(dh1 + 7);
    const auto *dh1_8 = buffer.data(dh1 + 8);
    const auto *dh1_9 = buffer.data(dh1 + 9);
    const auto *dh1_10 = buffer.data(dh1 + 10);
    const auto *dh1_11 = buffer.data(dh1 + 11);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pi_0, dh0_0, dh0_1, dh0_2, dh1_0, dh1_1, dh1_2, \
                         di_0, di_1, di_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_3 * dh0_1[k]
                 - f_4 * dh1_1[k]
                 + pb_x[k] * di_1[k];

        t_2[k] = f_5 * dh0_2[k]
                 - f_6 * dh1_2[k]
                 + pb_x[k] * di_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, pi_1, dh0_3, dh0_4, dh1_3, dh1_4, di_3, \
                         di_4, di_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_7 * dh0_3[k]
                 - f_8 * dh1_3[k]
                 + pb_x[k] * di_3[k];

        t_4[k] = f_9 * dh0_4[k]
                 - f_10 * dh1_4[k]
                 + pb_x[k] * di_4[k];

        t_5[k] = f_0 * pi_1[k]
                 + f_1 * dh0_4[k]
                 - f_2 * dh1_4[k]
                 + pb_y[k] * di_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, dh0_5, dh0_6, dh0_7, dh1_5, dh1_6, dh1_7, di_6, \
                         di_7, di_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * dh0_5[k]
                 - f_4 * dh1_5[k]
                 + pb_x[k] * di_6[k];

        t_7[k] = f_5 * dh0_6[k]
                 - f_6 * dh1_6[k]
                 + pb_x[k] * di_7[k];

        t_8[k] = f_7 * dh0_7[k]
                 - f_8 * dh1_7[k]
                 + pb_x[k] * di_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, pb_y, dh0_8, dh0_9, dh0_11, dh1_8, dh1_9, \
                         dh1_11, di_9, di_10, di_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_9 * dh0_11[k]
                 - f_10 * dh1_11[k]
                 + pb_x[k] * di_9[k];

        t_10[k] = f_3 * dh0_8[k]
                  - f_4 * dh1_8[k]
                  + pb_y[k] * di_10[k];

        t_11[k] = f_5 * dh0_9[k]
                  - f_6 * dh1_9[k]
                  + pb_y[k] * di_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, pb_z, pi_2, dh0_10, dh0_11, dh1_10, dh1_11, \
                         di_12, di_13, di_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * dh0_10[k]
                  - f_8 * dh1_10[k]
                  + pb_y[k] * di_12[k];

        t_13[k] = f_9 * dh0_11[k]
                  - f_10 * dh1_11[k]
                  + pb_y[k] * di_13[k];

        t_14[k] = f_0 * pi_2[k]
                  + f_1 * dh0_11[k]
                  - f_2 * dh1_11[k]
                  + pb_z[k] * di_14[k];
    }
}

auto
compute_prim_dk_electron_repulsion_22(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pi, const size_t dh0, const size_t dh1,
                                      const size_t di, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 2.0 / beta;
    const auto f_4 = 2.0 * alpha / (beta * p);
    const auto f_5 = 1.5 / beta;
    const auto f_6 = 1.5 * alpha / (beta * p);
    const auto f_7 = 1.0 / beta;
    const auto f_8 = alpha / (beta * p);
    const auto f_9 = 0.5 / beta;
    const auto f_10 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);
    const auto *dh0_3 = buffer.data(dh0 + 3);
    const auto *dh0_4 = buffer.data(dh0 + 4);
    const auto *dh0_5 = buffer.data(dh0 + 5);
    const auto *dh0_6 = buffer.data(dh0 + 6);
    const auto *dh0_7 = buffer.data(dh0 + 7);
    const auto *dh0_8 = buffer.data(dh0 + 8);
    const auto *dh0_9 = buffer.data(dh0 + 9);
    const auto *dh0_10 = buffer.data(dh0 + 10);
    const auto *dh0_11 = buffer.data(dh0 + 11);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_5 = buffer.data(dh1 + 5);
    const auto *dh1_6 = buffer.data(dh1 + 6);
    const auto *dh1_7 = buffer.data(dh1 + 7);
    const auto *dh1_10 = buffer.data(dh1 + 10);
    const auto *dh1_11 = buffer.data(dh1 + 11);
    const auto *dh1_12 = buffer.data(dh1 + 12);
    const auto *dh1_14 = buffer.data(dh1 + 14);
    const auto *dh1_15 = buffer.data(dh1 + 15);
    const auto *dh1_16 = buffer.data(dh1 + 16);
    const auto *dh1_17 = buffer.data(dh1 + 17);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pi_0, dh0_0, dh0_1, dh0_2, dh1_0, dh1_4, dh1_5, \
                         di_0, di_4, di_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_3 * dh0_1[k]
                 - f_4 * dh1_4[k]
                 + pb_x[k] * di_4[k];

        t_2[k] = f_5 * dh0_2[k]
                 - f_6 * dh1_5[k]
                 + pb_x[k] * di_5[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, pi_1, dh0_3, dh0_4, dh1_6, dh1_7, di_6, \
                         di_7, di_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_7 * dh0_3[k]
                 - f_8 * dh1_6[k]
                 + pb_x[k] * di_6[k];

        t_4[k] = f_9 * dh0_4[k]
                 - f_10 * dh1_7[k]
                 + pb_x[k] * di_7[k];

        t_5[k] = f_0 * pi_1[k]
                 + f_1 * dh0_4[k]
                 - f_2 * dh1_7[k]
                 + pb_y[k] * di_8[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, dh0_5, dh0_6, dh0_7, dh1_10, dh1_11, dh1_12, \
                         di_11, di_12, di_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * dh0_5[k]
                 - f_4 * dh1_10[k]
                 + pb_x[k] * di_11[k];

        t_7[k] = f_5 * dh0_6[k]
                 - f_6 * dh1_11[k]
                 + pb_x[k] * di_12[k];

        t_8[k] = f_7 * dh0_7[k]
                 - f_8 * dh1_12[k]
                 + pb_x[k] * di_13[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, pb_y, dh0_8, dh0_9, dh0_11, dh1_14, dh1_15, \
                         dh1_17, di_14, di_16, di_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_9 * dh0_11[k]
                 - f_10 * dh1_17[k]
                 + pb_x[k] * di_14[k];

        t_10[k] = f_3 * dh0_8[k]
                  - f_4 * dh1_14[k]
                  + pb_y[k] * di_16[k];

        t_11[k] = f_5 * dh0_9[k]
                  - f_6 * dh1_15[k]
                  + pb_y[k] * di_17[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, pb_z, pi_2, dh0_10, dh0_11, dh1_16, dh1_17, \
                         di_18, di_19, di_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * dh0_10[k]
                  - f_8 * dh1_16[k]
                  + pb_y[k] * di_18[k];

        t_13[k] = f_9 * dh0_11[k]
                  - f_10 * dh1_17[k]
                  + pb_y[k] * di_19[k];

        t_14[k] = f_0 * pi_2[k]
                  + f_1 * dh0_11[k]
                  - f_2 * dh1_17[k]
                  + pb_z[k] * di_20[k];
    }
}

auto
compute_prim_dk_electron_repulsion_23(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pi, const size_t dh0, const size_t dh1,
                                      const size_t di, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 2.0 / beta;
    const auto f_4 = 2.0 * alpha / (beta * p);
    const auto f_5 = 1.5 / beta;
    const auto f_6 = 1.5 * alpha / (beta * p);
    const auto f_7 = 1.0 / beta;
    const auto f_8 = alpha / (beta * p);
    const auto f_9 = 0.5 / beta;
    const auto f_10 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_4 = buffer.data(dh0 + 4);
    const auto *dh0_5 = buffer.data(dh0 + 5);
    const auto *dh0_6 = buffer.data(dh0 + 6);
    const auto *dh0_7 = buffer.data(dh0 + 7);
    const auto *dh0_10 = buffer.data(dh0 + 10);
    const auto *dh0_11 = buffer.data(dh0 + 11);
    const auto *dh0_12 = buffer.data(dh0 + 12);
    const auto *dh0_14 = buffer.data(dh0 + 14);
    const auto *dh0_15 = buffer.data(dh0 + 15);
    const auto *dh0_16 = buffer.data(dh0 + 16);
    const auto *dh0_17 = buffer.data(dh0 + 17);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_5 = buffer.data(dh1 + 5);
    const auto *dh1_6 = buffer.data(dh1 + 6);
    const auto *dh1_7 = buffer.data(dh1 + 7);
    const auto *dh1_10 = buffer.data(dh1 + 10);
    const auto *dh1_11 = buffer.data(dh1 + 11);
    const auto *dh1_12 = buffer.data(dh1 + 12);
    const auto *dh1_14 = buffer.data(dh1 + 14);
    const auto *dh1_15 = buffer.data(dh1 + 15);
    const auto *dh1_16 = buffer.data(dh1 + 16);
    const auto *dh1_17 = buffer.data(dh1 + 17);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pi_0, dh0_0, dh0_4, dh0_5, dh1_0, dh1_4, dh1_5, \
                         di_0, di_1, di_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_3 * dh0_4[k]
                 - f_4 * dh1_4[k]
                 + pb_x[k] * di_1[k];

        t_2[k] = f_5 * dh0_5[k]
                 - f_6 * dh1_5[k]
                 + pb_x[k] * di_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, pi_1, dh0_6, dh0_7, dh1_6, dh1_7, di_3, \
                         di_4, di_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_7 * dh0_6[k]
                 - f_8 * dh1_6[k]
                 + pb_x[k] * di_3[k];

        t_4[k] = f_9 * dh0_7[k]
                 - f_10 * dh1_7[k]
                 + pb_x[k] * di_4[k];

        t_5[k] = f_0 * pi_1[k]
                 + f_1 * dh0_7[k]
                 - f_2 * dh1_7[k]
                 + pb_y[k] * di_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, dh0_10, dh0_11, dh0_12, dh1_10, dh1_11, dh1_12, \
                         di_6, di_7, di_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * dh0_10[k]
                 - f_4 * dh1_10[k]
                 + pb_x[k] * di_6[k];

        t_7[k] = f_5 * dh0_11[k]
                 - f_6 * dh1_11[k]
                 + pb_x[k] * di_7[k];

        t_8[k] = f_7 * dh0_12[k]
                 - f_8 * dh1_12[k]
                 + pb_x[k] * di_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, pb_y, dh0_14, dh0_15, dh0_17, dh1_14, dh1_15, \
                         dh1_17, di_9, di_10, di_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_9 * dh0_17[k]
                 - f_10 * dh1_17[k]
                 + pb_x[k] * di_9[k];

        t_10[k] = f_3 * dh0_14[k]
                  - f_4 * dh1_14[k]
                  + pb_y[k] * di_10[k];

        t_11[k] = f_5 * dh0_15[k]
                  - f_6 * dh1_15[k]
                  + pb_y[k] * di_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, pb_z, pi_2, dh0_16, dh0_17, dh1_16, dh1_17, \
                         di_12, di_13, di_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * dh0_16[k]
                  - f_8 * dh1_16[k]
                  + pb_y[k] * di_12[k];

        t_13[k] = f_9 * dh0_17[k]
                  - f_10 * dh1_17[k]
                  + pb_y[k] * di_13[k];

        t_14[k] = f_0 * pi_2[k]
                  + f_1 * dh0_17[k]
                  - f_2 * dh1_17[k]
                  + pb_z[k] * di_14[k];
    }
}

auto
compute_prim_dk_electron_repulsion_24(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pi, const size_t dh0, const size_t dh1,
                                      const size_t di, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_4 = buffer.data(dh0 + 4);
    const auto *dh0_11 = buffer.data(dh0 + 11);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_4 = buffer.data(dh1 + 4);
    const auto *dh1_11 = buffer.data(dh1 + 11);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pi_0, pi_1, dh0_0, dh0_4, dh1_0, dh1_4, di_0, \
                         di_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pi_0[k]
                 + f_1 * dh0_0[k]
                 - f_2 * dh1_0[k]
                 + pb_x[k] * di_0[k];

        t_1[k] = f_0 * pi_1[k]
                 + f_1 * dh0_4[k]
                 - f_2 * dh1_4[k]
                 + pb_y[k] * di_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pi_2, dh0_11, dh1_11, di_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pi_2[k]
                 + f_1 * dh0_11[k]
                 - f_2 * dh1_11[k]
                 + pb_z[k] * di_2[k];
    }
}

}  // namespace simdt2ceri
