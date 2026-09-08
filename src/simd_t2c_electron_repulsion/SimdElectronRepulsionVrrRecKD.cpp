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


#include "SimdElectronRepulsionVrrRecKD.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_kd_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hd0, const size_t hd1,
                                     const size_t ip, const size_t id, const size_t ks0,
                                     const size_t ks1, const size_t kp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 3.0 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 2.5 / p;
    const auto f_8 = 2.0 / alpha;
    const auto f_9 = 2.0 * beta / (alpha * p);
    const auto f_10 = 0.5 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);
    const auto f_13 = 2.0 / p;
    const auto f_14 = 1.5 / alpha;
    const auto f_15 = 1.5 * beta / (alpha * p);
    const auto f_16 = 1.5 / p;

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_23 = buffer.data(hd0 + 23);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_23 = buffer.data(hd1 + 23);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_21 = buffer.data(ip + 21);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_27 = buffer.data(ip + 27);
    const auto *ip_28 = buffer.data(ip + 28);
    const auto *ip_29 = buffer.data(ip + 29);
    const auto *ip_30 = buffer.data(ip + 30);
    const auto *ip_31 = buffer.data(ip + 31);
    const auto *ip_32 = buffer.data(ip + 32);
    const auto *ip_33 = buffer.data(ip + 33);
    const auto *ip_34 = buffer.data(ip + 34);
    const auto *ip_35 = buffer.data(ip + 35);
    const auto *ip_36 = buffer.data(ip + 36);
    const auto *ip_37 = buffer.data(ip + 37);
    const auto *ip_38 = buffer.data(ip + 38);
    const auto *ip_39 = buffer.data(ip + 39);
    const auto *ip_40 = buffer.data(ip + 40);
    const auto *ip_41 = buffer.data(ip + 41);
    const auto *ip_42 = buffer.data(ip + 42);
    const auto *ip_43 = buffer.data(ip + 43);
    const auto *ip_44 = buffer.data(ip + 44);
    const auto *ip_45 = buffer.data(ip + 45);
    const auto *ip_46 = buffer.data(ip + 46);
    const auto *ip_47 = buffer.data(ip + 47);
    const auto *ip_48 = buffer.data(ip + 48);
    const auto *ip_49 = buffer.data(ip + 49);
    const auto *ip_50 = buffer.data(ip + 50);
    const auto *ip_51 = buffer.data(ip + 51);

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

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_9 = buffer.data(ks0 + 9);
    const auto *ks0_10 = buffer.data(ks0 + 10);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_9 = buffer.data(ks1 + 9);
    const auto *ks1_10 = buffer.data(ks1 + 10);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_21 = buffer.data(kp + 21);
    const auto *kp_22 = buffer.data(kp + 22);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_24 = buffer.data(kp + 24);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_28 = buffer.data(kp + 28);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_30 = buffer.data(kp + 30);
    const auto *kp_31 = buffer.data(kp + 31);
    const auto *kp_32 = buffer.data(kp + 32);
    const auto *kp_33 = buffer.data(kp + 33);
    const auto *kp_34 = buffer.data(kp + 34);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_36 = buffer.data(kp + 36);
    const auto *kp_37 = buffer.data(kp + 37);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_39 = buffer.data(kp + 39);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_42 = buffer.data(kp + 42);
    const auto *kp_43 = buffer.data(kp + 43);
    const auto *kp_44 = buffer.data(kp + 44);
    const auto *kp_45 = buffer.data(kp + 45);
    const auto *kp_46 = buffer.data(kp + 46);
    const auto *kp_47 = buffer.data(kp + 47);
    const auto *kp_48 = buffer.data(kp + 48);
    const auto *kp_49 = buffer.data(kp + 49);
    const auto *kp_50 = buffer.data(kp + 50);
    const auto *kp_51 = buffer.data(kp + 51);
    const auto *kp_52 = buffer.data(kp + 52);
    const auto *kp_53 = buffer.data(kp + 53);
    const auto *kp_54 = buffer.data(kp + 54);
    const auto *kp_55 = buffer.data(kp + 55);
    const auto *kp_56 = buffer.data(kp + 56);
    const auto *kp_57 = buffer.data(kp + 57);
    const auto *kp_58 = buffer.data(kp + 58);
    const auto *kp_59 = buffer.data(kp + 59);
    const auto *kp_60 = buffer.data(kp + 60);
    const auto *kp_61 = buffer.data(kp + 61);
    const auto *kp_62 = buffer.data(kp + 62);
    const auto *kp_63 = buffer.data(kp + 63);
    const auto *kp_64 = buffer.data(kp + 64);
    const auto *kp_65 = buffer.data(kp + 65);
    const auto *kp_66 = buffer.data(kp + 66);
    const auto *kp_67 = buffer.data(kp + 67);
    const auto *kp_68 = buffer.data(kp + 68);
    const auto *kp_69 = buffer.data(kp + 69);
    const auto *kp_70 = buffer.data(kp + 70);
    const auto *kp_71 = buffer.data(kp + 71);
    const auto *kp_72 = buffer.data(kp + 72);
    const auto *kp_73 = buffer.data(kp + 73);
    const auto *kp_74 = buffer.data(kp + 74);
    const auto *kp_75 = buffer.data(kp + 75);
    const auto *kp_76 = buffer.data(kp + 76);
    const auto *kp_77 = buffer.data(kp + 77);
    const auto *kp_78 = buffer.data(kp + 78);
    const auto *kp_79 = buffer.data(kp + 79);
    const auto *kp_80 = buffer.data(kp + 80);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, ip_0, ks0_0, ks1_0, \
                         kp_0, kp_1, kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = pb_y[k] * kp_0[k];

        t_2[k] = pb_z[k] * kp_0[k];

        t_3[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_4[k] = pb_y[k] * kp_2[k];

        t_5[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pa_y, pb_x, pb_z, ip_1, ip_3, id_0, \
                         id_1, id_2, kp_3, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * id_0[k];

        t_7[k] = f_3 * ip_3[k]
                 + pb_x[k] * kp_4[k];

        t_8[k] = pb_z[k] * kp_3[k];

        t_9[k] = f_4 * ip_1[k]
                 + pa_y[k] * id_1[k];

        t_10[k] = pb_z[k] * kp_4[k];

        t_11[k] = pa_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, ip_2, ip_4, \
                         id_0, id_1, id_2, kp_5, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_z[k] * id_0[k];

        t_13[k] = pb_y[k] * kp_5[k];

        t_14[k] = f_3 * ip_4[k]
                  + pb_x[k] * kp_6[k];

        t_15[k] = pa_z[k] * id_1[k];

        t_16[k] = pb_y[k] * kp_6[k];

        t_17[k] = f_4 * ip_2[k]
                  + pa_z[k] * id_2[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pb_x, pb_z, hd0_0, hd1_0, ip_5, id_3, kp_7, \
                         kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * hd0_0[k]
                  - f_6 * hd1_0[k]
                  + pa_y[k] * id_3[k];

        t_19[k] = f_7 * ip_5[k]
                  + pb_x[k] * kp_8[k];

        t_20[k] = pb_z[k] * kp_7[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_x, pa_y, pb_z, hd0_4, hd1_4, id_6, id_11, \
                         ks0_1, ks1_1, kp_8, kp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * hd0_4[k]
                  - f_9 * hd1_4[k]
                  + pa_x[k] * id_11[k];

        t_22[k] = pb_z[k] * kp_8[k];

        t_23[k] = f_1 * ks0_1[k]
                  - f_2 * ks1_1[k]
                  + pb_z[k] * kp_9[k];

        t_24[k] = pa_y[k] * id_6[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pa_z, pb_y, ip_4, id_4, id_5, \
                         id_7, id_8, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_z[k] * id_4[k];

        t_26[k] = pa_y[k] * id_7[k];

        t_27[k] = pa_z[k] * id_5[k];

        t_28[k] = f_10 * ip_4[k]
                  + pb_y[k] * kp_10[k];

        t_29[k] = pa_y[k] * id_8[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_z, pb_x, pb_y, hd0_0, hd1_0, ip_9, id_6, \
                         ks0_2, ks1_2, kp_11, kp_12, kp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * hd0_0[k]
                  - f_6 * hd1_0[k]
                  + pa_z[k] * id_6[k];

        t_31[k] = pb_y[k] * kp_11[k];

        t_32[k] = f_7 * ip_9[k]
                  + pb_x[k] * kp_13[k];

        t_33[k] = f_1 * ks0_2[k]
                  - f_2 * ks1_2[k]
                  + pb_y[k] * kp_12[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_x, pa_y, pb_y, hd0_1, hd0_6, hd1_1, hd1_6, id_9, \
                         id_16, kp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * kp_13[k];

        t_35[k] = f_8 * hd0_6[k]
                  - f_9 * hd1_6[k]
                  + pa_x[k] * id_16[k];

        t_36[k] = f_11 * hd0_1[k]
                  - f_12 * hd1_1[k]
                  + pa_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_x, pb_x, pb_z, hd0_8, hd1_8, ip_10, id_19, \
                         kp_14, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_13 * ip_10[k]
                  + pb_x[k] * kp_15[k];

        t_38[k] = pb_z[k] * kp_14[k];

        t_39[k] = f_14 * hd0_8[k]
                  - f_15 * hd1_8[k]
                  + pa_x[k] * id_19[k];

        t_40[k] = pb_z[k] * kp_15[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_z, pb_x, pb_z, ip_12, id_9, id_10, \
                         id_11, ks0_3, ks1_3, kp_16, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_16[k];

        t_42[k] = pa_z[k] * id_9[k];

        t_43[k] = pa_z[k] * id_10[k];

        t_44[k] = f_13 * ip_12[k]
                  + pb_x[k] * kp_17[k];

        t_45[k] = pa_z[k] * id_11[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pa_z, pb_x, pb_y, ip_6, ip_7, ip_13, \
                         id_12, id_13, kp_17, kp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_4 * ip_7[k]
                  + pb_y[k] * kp_17[k];

        t_47[k] = f_4 * ip_6[k]
                  + pa_z[k] * id_12[k];

        t_48[k] = pa_y[k] * id_13[k];

        t_49[k] = f_13 * ip_13[k]
                  + pb_x[k] * kp_18[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pb_y, ip_8, ip_9, id_14, id_15, id_16, \
                         kp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * id_14[k];

        t_51[k] = f_4 * ip_8[k]
                  + pa_y[k] * id_15[k];

        t_52[k] = f_10 * ip_9[k]
                  + pb_y[k] * kp_19[k];

        t_53[k] = pa_y[k] * id_16[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_z, pb_x, pb_y, hd0_2, hd1_2, ip_16, id_13, \
                         ks0_4, ks1_4, kp_20, kp_21, kp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * hd0_2[k]
                  - f_12 * hd1_2[k]
                  + pa_z[k] * id_13[k];

        t_55[k] = pb_y[k] * kp_20[k];

        t_56[k] = f_13 * ip_16[k]
                  + pb_x[k] * kp_22[k];

        t_57[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_21[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_x, pa_y, pb_y, hd0_3, hd0_11, hd1_3, hd1_11, \
                         id_17, id_25, kp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_y[k] * kp_22[k];

        t_59[k] = f_14 * hd0_11[k]
                  - f_15 * hd1_11[k]
                  + pa_x[k] * id_25[k];

        t_60[k] = f_14 * hd0_3[k]
                  - f_15 * hd1_3[k]
                  + pa_y[k] * id_17[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pb_x, pb_z, hd0_12, hd1_12, ip_17, \
                         id_28, kp_23, kp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_16 * ip_17[k]
                  + pb_x[k] * kp_24[k];

        t_62[k] = pb_z[k] * kp_23[k];

        t_63[k] = f_11 * hd0_12[k]
                  - f_12 * hd1_12[k]
                  + pa_x[k] * id_28[k];

        t_64[k] = pb_z[k] * kp_24[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pa_z, pb_x, pb_z, ip_19, id_17, id_18, \
                         id_19, ks0_5, ks1_5, kp_25, kp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_25[k];

        t_66[k] = pa_z[k] * id_17[k];

        t_67[k] = pa_z[k] * id_18[k];

        t_68[k] = f_16 * ip_19[k]
                  + pb_x[k] * kp_26[k];

        t_69[k] = pa_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_y, pa_z, pb_y, hd0_5, hd1_5, ip_11, ip_12, \
                         id_20, id_21, kp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_16 * ip_12[k]
                  + pb_y[k] * kp_26[k];

        t_71[k] = f_4 * ip_11[k]
                  + pa_z[k] * id_20[k];

        t_72[k] = f_5 * hd0_5[k]
                  - f_6 * hd1_5[k]
                  + pa_y[k] * id_21[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pb_x, pb_y, hd0_13, hd1_13, ip_14, \
                         ip_20, ip_21, id_31, kp_27, kp_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_16 * ip_20[k]
                  + pb_x[k] * kp_27[k];

        t_74[k] = f_16 * ip_21[k]
                  + pb_x[k] * kp_28[k];

        t_75[k] = f_11 * hd0_13[k]
                  - f_12 * hd1_13[k]
                  + pa_x[k] * id_31[k];

        t_76[k] = f_4 * ip_14[k]
                  + pb_y[k] * kp_28[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_x, pa_y, pb_x, hd0_14, hd1_14, ip_22, \
                         id_22, id_23, id_32, kp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_11 * hd0_14[k]
                  - f_12 * hd1_14[k]
                  + pa_x[k] * id_32[k];

        t_78[k] = pa_y[k] * id_22[k];

        t_79[k] = f_16 * ip_22[k]
                  + pb_x[k] * kp_29[k];

        t_80[k] = pa_y[k] * id_23[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_y, pa_z, pb_y, hd0_5, hd1_5, ip_15, ip_16, \
                         id_22, id_24, id_25, kp_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * ip_15[k]
                  + pa_y[k] * id_24[k];

        t_82[k] = f_10 * ip_16[k]
                  + pb_y[k] * kp_30[k];

        t_83[k] = pa_y[k] * id_25[k];

        t_84[k] = f_14 * hd0_5[k]
                  - f_15 * hd1_5[k]
                  + pa_z[k] * id_22[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_x, pb_y, ip_25, ks0_6, ks1_6, kp_31, \
                         kp_32, kp_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pb_y[k] * kp_31[k];

        t_86[k] = f_16 * ip_25[k]
                  + pb_x[k] * kp_33[k];

        t_87[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_32[k];

        t_88[k] = pb_y[k] * kp_33[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_x, pa_y, pb_x, hd0_7, hd0_15, hd1_7, hd1_15, \
                         ip_26, id_26, id_37, kp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_11 * hd0_15[k]
                  - f_12 * hd1_15[k]
                  + pa_x[k] * id_37[k];

        t_90[k] = f_8 * hd0_7[k]
                  - f_9 * hd1_7[k]
                  + pa_y[k] * id_26[k];

        t_91[k] = f_4 * ip_26[k]
                  + pb_x[k] * kp_35[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_x, pb_z, hd0_16, hd1_16, id_40, ks0_7, \
                         ks1_7, kp_34, kp_35, kp_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_z[k] * kp_34[k];

        t_93[k] = f_5 * hd0_16[k]
                  - f_6 * hd1_16[k]
                  + pa_x[k] * id_40[k];

        t_94[k] = pb_z[k] * kp_35[k];

        t_95[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_36[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pa_z, pb_x, pb_y, ip_19, ip_27, id_26, \
                         id_27, id_28, kp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * id_26[k];

        t_97[k] = pa_z[k] * id_27[k];

        t_98[k] = f_4 * ip_27[k]
                  + pb_x[k] * kp_37[k];

        t_99[k] = pa_z[k] * id_28[k];

        t_100[k] = f_13 * ip_19[k]
                   + pb_y[k] * kp_37[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_y, pa_z, pb_x, hd0_9, hd1_9, ip_18, \
                         ip_28, ip_29, id_29, id_30, kp_38, kp_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_4 * ip_18[k]
                   + pa_z[k] * id_29[k];

        t_102[k] = f_11 * hd0_9[k]
                   - f_12 * hd1_9[k]
                   + pa_y[k] * id_30[k];

        t_103[k] = f_4 * ip_28[k]
                   + pb_x[k] * kp_38[k];

        t_104[k] = f_4 * ip_29[k]
                   + pb_x[k] * kp_39[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pa_x, pb_y, hd0_18, hd0_19, hd1_18, hd1_19, \
                         ip_21, id_41, id_42, kp_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_5 * hd0_18[k]
                   - f_6 * hd1_18[k]
                   + pa_x[k] * id_41[k];

        t_106[k] = f_16 * ip_21[k]
                   + pb_y[k] * kp_39[k];

        t_107[k] = f_5 * hd0_19[k]
                   - f_6 * hd1_19[k]
                   + pa_x[k] * id_42[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_y, pb_x, hd0_10, hd1_10, ip_30, ip_31, id_33, \
                         kp_40, kp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_5 * hd0_10[k]
                   - f_6 * hd1_10[k]
                   + pa_y[k] * id_33[k];

        t_109[k] = f_4 * ip_30[k]
                   + pb_x[k] * kp_40[k];

        t_110[k] = f_4 * ip_31[k]
                   + pb_x[k] * kp_41[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_x, pa_y, pb_y, hd0_20, hd0_21, hd1_20, \
                         hd1_21, ip_23, id_34, id_43, id_44, kp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_5 * hd0_20[k]
                   - f_6 * hd1_20[k]
                   + pa_x[k] * id_43[k];

        t_112[k] = f_4 * ip_23[k]
                   + pb_y[k] * kp_41[k];

        t_113[k] = f_5 * hd0_21[k]
                   - f_6 * hd1_21[k]
                   + pa_x[k] * id_44[k];

        t_114[k] = pa_y[k] * id_34[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pa_y, pb_x, pb_y, ip_24, ip_25, \
                         ip_32, id_35, id_36, id_37, kp_42, kp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_4 * ip_32[k]
                   + pb_x[k] * kp_42[k];

        t_116[k] = pa_y[k] * id_35[k];

        t_117[k] = f_4 * ip_24[k]
                   + pa_y[k] * id_36[k];

        t_118[k] = f_10 * ip_25[k]
                   + pb_y[k] * kp_43[k];

        t_119[k] = pa_y[k] * id_37[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_z, pb_x, pb_y, hd0_10, hd1_10, ip_33, \
                         id_34, ks0_8, ks1_8, kp_44, kp_45, kp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_8 * hd0_10[k]
                   - f_9 * hd1_10[k]
                   + pa_z[k] * id_34[k];

        t_121[k] = pb_y[k] * kp_44[k];

        t_122[k] = f_4 * ip_33[k]
                   + pb_x[k] * kp_46[k];

        t_123[k] = f_1 * ks0_8[k]
                   - f_2 * ks1_8[k]
                   + pb_y[k] * kp_45[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_x, pb_x, pb_y, hd0_23, hd1_23, ip_34, \
                         ip_35, id_47, id_48, kp_46, kp_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * kp_46[k];

        t_125[k] = f_5 * hd0_23[k]
                   - f_6 * hd1_23[k]
                   + pa_x[k] * id_47[k];

        t_126[k] = f_4 * ip_34[k]
                   + pa_x[k] * id_48[k];

        t_127[k] = f_10 * ip_35[k]
                   + pb_x[k] * kp_48[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, pa_x, pa_z, pb_z, id_38, \
                         id_39, id_49, id_50, kp_47, kp_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = pb_z[k] * kp_47[k];

        t_129[k] = pa_x[k] * id_49[k];

        t_130[k] = pb_z[k] * kp_48[k];

        t_131[k] = pa_x[k] * id_50[k];

        t_132[k] = pa_z[k] * id_38[k];

        t_133[k] = pa_z[k] * id_39[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, pa_x, pb_x, ip_37, ip_38, id_51, \
                         id_52, id_53, id_54, kp_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_10 * ip_37[k]
                   + pb_x[k] * kp_49[k];

        t_135[k] = pa_x[k] * id_51[k];

        t_136[k] = pa_x[k] * id_52[k];

        t_137[k] = pa_x[k] * id_53[k];

        t_138[k] = f_4 * ip_38[k]
                   + pa_x[k] * id_54[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pa_x, pb_x, ip_39, ip_40, id_55, \
                         id_56, id_57, kp_50, kp_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_10 * ip_39[k]
                   + pb_x[k] * kp_50[k];

        t_140[k] = f_10 * ip_40[k]
                   + pb_x[k] * kp_51[k];

        t_141[k] = pa_x[k] * id_55[k];

        t_142[k] = pa_x[k] * id_56[k];

        t_143[k] = pa_x[k] * id_57[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, pa_x, pb_x, ip_41, ip_42, ip_43, \
                         id_58, id_59, id_60, kp_52, kp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_4 * ip_41[k]
                   + pa_x[k] * id_58[k];

        t_145[k] = f_10 * ip_42[k]
                   + pb_x[k] * kp_52[k];

        t_146[k] = f_10 * ip_43[k]
                   + pb_x[k] * kp_53[k];

        t_147[k] = pa_x[k] * id_59[k];

        t_148[k] = pa_x[k] * id_60[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, pa_x, pb_x, ip_44, ip_45, ip_46, \
                         id_61, id_62, id_63, kp_54, kp_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pa_x[k] * id_61[k];

        t_150[k] = f_4 * ip_44[k]
                   + pa_x[k] * id_62[k];

        t_151[k] = f_10 * ip_45[k]
                   + pb_x[k] * kp_54[k];

        t_152[k] = f_10 * ip_46[k]
                   + pb_x[k] * kp_55[k];

        t_153[k] = pa_x[k] * id_63[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, t_159, pa_x, pa_y, pb_x, ip_47, \
                         id_45, id_46, id_64, id_65, id_66, kp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pa_x[k] * id_64[k];

        t_155[k] = pa_x[k] * id_65[k];

        t_156[k] = pa_y[k] * id_45[k];

        t_157[k] = f_10 * ip_47[k]
                   + pb_x[k] * kp_56[k];

        t_158[k] = pa_y[k] * id_46[k];

        t_159[k] = pa_x[k] * id_66[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, pa_x, pb_x, pb_y, ip_49, ip_51, \
                         id_67, id_68, id_69, kp_57, kp_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = pa_x[k] * id_67[k];

        t_161[k] = pa_x[k] * id_68[k];

        t_162[k] = f_4 * ip_49[k]
                   + pa_x[k] * id_69[k];

        t_163[k] = pb_y[k] * kp_57[k];

        t_164[k] = f_10 * ip_51[k]
                   + pb_x[k] * kp_58[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pa_x, pb_x, pb_y, id_70, id_71, \
                         ks0_9, ks1_9, kp_58, kp_59, kp_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = pa_x[k] * id_70[k];

        t_166[k] = pb_y[k] * kp_58[k];

        t_167[k] = pa_x[k] * id_71[k];

        t_168[k] = f_1 * ks0_9[k]
                   - f_2 * ks1_9[k]
                   + pb_x[k] * kp_59[k];

        t_169[k] = pb_x[k] * kp_60[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, pa_z, pb_x, pb_y, pb_z, ip_35, \
                         id_48, ks0_9, ks1_9, kp_60, kp_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = pb_x[k] * kp_61[k];

        t_171[k] = f_0 * ip_35[k]
                   + f_1 * ks0_9[k]
                   - f_2 * ks1_9[k]
                   + pb_y[k] * kp_60[k];

        t_172[k] = pb_z[k] * kp_60[k];

        t_173[k] = f_1 * ks0_9[k]
                   - f_2 * ks1_9[k]
                   + pb_z[k] * kp_61[k];

        t_174[k] = pa_z[k] * id_48[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, pa_z, pb_x, pb_y, ip_36, ip_37, \
                         id_49, id_50, kp_62, kp_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pb_x[k] * kp_62[k];

        t_176[k] = pb_x[k] * kp_63[k];

        t_177[k] = pa_z[k] * id_49[k];

        t_178[k] = f_3 * ip_37[k]
                   + pb_y[k] * kp_63[k];

        t_179[k] = f_4 * ip_36[k]
                   + pa_z[k] * id_50[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_z, pb_x, hd0_16, hd1_16, id_51, \
                         ks0_10, ks1_10, kp_64, kp_65, kp_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_1 * ks0_10[k]
                   - f_2 * ks1_10[k]
                   + pb_x[k] * kp_64[k];

        t_181[k] = pb_x[k] * kp_65[k];

        t_182[k] = pb_x[k] * kp_66[k];

        t_183[k] = f_5 * hd0_16[k]
                   - f_6 * hd1_16[k]
                   + pa_z[k] * id_51[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pb_x, pb_y, hd0_19, hd1_19, ip_40, \
                         id_57, ks0_11, ks1_11, kp_66, kp_67, kp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_7 * ip_40[k]
                   + pb_y[k] * kp_66[k];

        t_185[k] = f_8 * hd0_19[k]
                   - f_9 * hd1_19[k]
                   + pa_y[k] * id_57[k];

        t_186[k] = f_1 * ks0_11[k]
                   - f_2 * ks1_11[k]
                   + pb_x[k] * kp_67[k];

        t_187[k] = pb_x[k] * kp_68[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_y, pa_z, pb_x, pb_y, hd0_17, hd0_21, \
                         hd1_17, hd1_21, ip_43, id_55, id_61, kp_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pb_x[k] * kp_69[k];

        t_189[k] = f_11 * hd0_17[k]
                   - f_12 * hd1_17[k]
                   + pa_z[k] * id_55[k];

        t_190[k] = f_13 * ip_43[k]
                   + pb_y[k] * kp_69[k];

        t_191[k] = f_14 * hd0_21[k]
                   - f_15 * hd1_21[k]
                   + pa_y[k] * id_61[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_z, pb_x, hd0_18, hd1_18, id_59, \
                         ks0_12, ks1_12, kp_70, kp_71, kp_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_1 * ks0_12[k]
                   - f_2 * ks1_12[k]
                   + pb_x[k] * kp_70[k];

        t_193[k] = pb_x[k] * kp_71[k];

        t_194[k] = pb_x[k] * kp_72[k];

        t_195[k] = f_14 * hd0_18[k]
                   - f_15 * hd1_18[k]
                   + pa_z[k] * id_59[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_y, pb_x, pb_y, hd0_22, hd1_22, ip_46, \
                         id_65, ks0_13, ks1_13, kp_72, kp_73, kp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_16 * ip_46[k]
                   + pb_y[k] * kp_72[k];

        t_197[k] = f_11 * hd0_22[k]
                   - f_12 * hd1_22[k]
                   + pa_y[k] * id_65[k];

        t_198[k] = f_1 * ks0_13[k]
                   - f_2 * ks1_13[k]
                   + pb_x[k] * kp_73[k];

        t_199[k] = pb_x[k] * kp_74[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_y, pa_z, pb_x, pb_y, hd0_20, hd0_23, \
                         hd1_20, hd1_23, ip_48, id_63, id_68, kp_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pb_x[k] * kp_75[k];

        t_201[k] = f_8 * hd0_20[k]
                   - f_9 * hd1_20[k]
                   + pa_z[k] * id_63[k];

        t_202[k] = f_4 * ip_48[k]
                   + pb_y[k] * kp_75[k];

        t_203[k] = f_5 * hd0_23[k]
                   - f_6 * hd1_23[k]
                   + pa_y[k] * id_68[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, t_209, pa_y, pb_x, pb_y, ip_50, \
                         ip_51, id_69, id_70, id_71, kp_76, kp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_y[k] * id_69[k];

        t_205[k] = pb_x[k] * kp_76[k];

        t_206[k] = pb_x[k] * kp_77[k];

        t_207[k] = f_4 * ip_50[k]
                   + pa_y[k] * id_70[k];

        t_208[k] = f_10 * ip_51[k]
                   + pb_y[k] * kp_77[k];

        t_209[k] = pa_y[k] * id_71[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, pb_x, pb_y, pb_z, ip_51, \
                         ks0_14, ks1_14, kp_78, kp_79, kp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_1 * ks0_14[k]
                   - f_2 * ks1_14[k]
                   + pb_x[k] * kp_78[k];

        t_211[k] = pb_x[k] * kp_79[k];

        t_212[k] = pb_x[k] * kp_80[k];

        t_213[k] = f_1 * ks0_14[k]
                   - f_2 * ks1_14[k]
                   + pb_y[k] * kp_79[k];

        t_214[k] = pb_y[k] * kp_80[k];

        t_215[k] = f_0 * ip_51[k]
                   + f_1 * ks0_14[k]
                   - f_2 * ks1_14[k]
                   + pb_z[k] * kp_80[k];
    }
}

auto
compute_prim_kd_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hd0, const size_t hd1,
                                     const size_t ip, const size_t id, const size_t ks0,
                                     const size_t ks1, const size_t kp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 3.0 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 2.5 / p;
    const auto f_8 = 2.0 / alpha;
    const auto f_9 = 2.0 * beta / (alpha * p);
    const auto f_10 = 0.5 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);
    const auto f_13 = 2.0 / p;
    const auto f_14 = 1.5 / alpha;
    const auto f_15 = 1.5 * beta / (alpha * p);
    const auto f_16 = 1.5 / p;

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_23 = buffer.data(hd0 + 23);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_29 = buffer.data(hd1 + 29);
    const auto *hd1_33 = buffer.data(hd1 + 33);
    const auto *hd1_35 = buffer.data(hd1 + 35);
    const auto *hd1_37 = buffer.data(hd1 + 37);
    const auto *hd1_39 = buffer.data(hd1 + 39);
    const auto *hd1_42 = buffer.data(hd1 + 42);
    const auto *hd1_45 = buffer.data(hd1 + 45);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_21 = buffer.data(ip + 21);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_27 = buffer.data(ip + 27);
    const auto *ip_28 = buffer.data(ip + 28);
    const auto *ip_29 = buffer.data(ip + 29);
    const auto *ip_30 = buffer.data(ip + 30);
    const auto *ip_31 = buffer.data(ip + 31);
    const auto *ip_32 = buffer.data(ip + 32);
    const auto *ip_33 = buffer.data(ip + 33);
    const auto *ip_34 = buffer.data(ip + 34);
    const auto *ip_35 = buffer.data(ip + 35);
    const auto *ip_36 = buffer.data(ip + 36);
    const auto *ip_37 = buffer.data(ip + 37);
    const auto *ip_38 = buffer.data(ip + 38);

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

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_9 = buffer.data(ks0 + 9);
    const auto *ks0_10 = buffer.data(ks0 + 10);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_9 = buffer.data(ks1 + 9);
    const auto *ks1_10 = buffer.data(ks1 + 10);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_21 = buffer.data(kp + 21);
    const auto *kp_22 = buffer.data(kp + 22);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_24 = buffer.data(kp + 24);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_28 = buffer.data(kp + 28);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_30 = buffer.data(kp + 30);
    const auto *kp_31 = buffer.data(kp + 31);
    const auto *kp_32 = buffer.data(kp + 32);
    const auto *kp_33 = buffer.data(kp + 33);
    const auto *kp_34 = buffer.data(kp + 34);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_36 = buffer.data(kp + 36);
    const auto *kp_37 = buffer.data(kp + 37);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_39 = buffer.data(kp + 39);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_42 = buffer.data(kp + 42);
    const auto *kp_43 = buffer.data(kp + 43);
    const auto *kp_44 = buffer.data(kp + 44);
    const auto *kp_45 = buffer.data(kp + 45);
    const auto *kp_46 = buffer.data(kp + 46);
    const auto *kp_47 = buffer.data(kp + 47);
    const auto *kp_48 = buffer.data(kp + 48);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_y, pb_x, pb_y, pb_z, ip_0, id_0, ks0_0, \
                         ks1_0, kp_0, kp_1, kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = pb_z[k] * kp_0[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_3[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];

        t_4[k] = pa_y[k] * id_0[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_x, ip_1, ip_3, ip_4, id_0, \
                         id_1, id_2, kp_3, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * ip_3[k]
                 + pb_x[k] * kp_3[k];

        t_6[k] = f_4 * ip_1[k]
                 + pa_y[k] * id_1[k];

        t_7[k] = pa_y[k] * id_2[k];

        t_8[k] = pa_z[k] * id_0[k];

        t_9[k] = f_3 * ip_4[k]
                 + pb_x[k] * kp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_x, hd0_0, hd1_0, ip_2, ip_5, \
                         id_1, id_2, id_3, kp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * id_1[k];

        t_11[k] = f_4 * ip_2[k]
                  + pa_z[k] * id_2[k];

        t_12[k] = f_5 * hd0_0[k]
                  - f_6 * hd1_0[k]
                  + pa_y[k] * id_3[k];

        t_13[k] = f_7 * ip_5[k]
                  + pb_x[k] * kp_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_z, hd0_4, hd1_8, id_4, id_8, ks0_1, \
                         ks1_1, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_8 * hd0_4[k]
                  - f_9 * hd1_8[k]
                  + pa_x[k] * id_8[k];

        t_15[k] = f_1 * ks0_1[k]
                  - f_2 * ks1_1[k]
                  + pb_z[k] * kp_6[k];

        t_16[k] = pa_z[k] * id_4[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pa_z, pb_x, pb_y, hd0_0, hd1_0, ip_4, \
                         ip_9, id_5, id_6, kp_7, kp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_10 * ip_4[k]
                  + pb_y[k] * kp_7[k];

        t_18[k] = pa_y[k] * id_6[k];

        t_19[k] = f_5 * hd0_0[k]
                  - f_6 * hd1_0[k]
                  + pa_z[k] * id_5[k];

        t_20[k] = f_7 * ip_9[k]
                  + pb_x[k] * kp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_y, hd0_1, hd0_6, hd1_3, hd1_12, \
                         id_7, id_12, ks0_2, ks1_2, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * ks0_2[k]
                  - f_2 * ks1_2[k]
                  + pb_y[k] * kp_8[k];

        t_22[k] = f_8 * hd0_6[k]
                  - f_9 * hd1_12[k]
                  + pa_x[k] * id_12[k];

        t_23[k] = f_11 * hd0_1[k]
                  - f_12 * hd1_3[k]
                  + pa_y[k] * id_7[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pb_x, pb_z, hd0_8, hd1_14, ip_10, id_14, \
                         ks0_3, ks1_3, kp_10, kp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_13 * ip_10[k]
                  + pb_x[k] * kp_10[k];

        t_25[k] = f_14 * hd0_8[k]
                  - f_15 * hd1_14[k]
                  + pa_x[k] * id_14[k];

        t_26[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_11[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_y, pa_z, pb_y, ip_6, ip_7, id_7, \
                         id_8, id_9, id_10, kp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_z[k] * id_7[k];

        t_28[k] = pa_z[k] * id_8[k];

        t_29[k] = f_4 * ip_7[k]
                  + pb_y[k] * kp_12[k];

        t_30[k] = f_4 * ip_6[k]
                  + pa_z[k] * id_9[k];

        t_31[k] = pa_y[k] * id_10[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_y, pa_z, pb_y, hd0_2, hd1_5, ip_8, ip_9, \
                         id_10, id_11, id_12, kp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_4 * ip_8[k]
                  + pa_y[k] * id_11[k];

        t_33[k] = f_10 * ip_9[k]
                  + pb_y[k] * kp_13[k];

        t_34[k] = pa_y[k] * id_12[k];

        t_35[k] = f_11 * hd0_2[k]
                  - f_12 * hd1_5[k]
                  + pa_z[k] * id_10[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_x, pb_x, pb_y, hd0_11, hd1_19, ip_15, id_19, \
                         ks0_4, ks1_4, kp_14, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_13 * ip_15[k]
                  + pb_x[k] * kp_15[k];

        t_37[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_14[k];

        t_38[k] = f_14 * hd0_11[k]
                  - f_15 * hd1_19[k]
                  + pa_x[k] * id_19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_x, pa_y, pb_x, hd0_3, hd0_12, hd1_7, hd1_21, \
                         ip_16, id_13, id_21, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_14 * hd0_3[k]
                  - f_15 * hd1_7[k]
                  + pa_y[k] * id_13[k];

        t_40[k] = f_16 * ip_16[k]
                  + pb_x[k] * kp_16[k];

        t_41[k] = f_11 * hd0_12[k]
                  - f_12 * hd1_21[k]
                  + pa_x[k] * id_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_z, pb_y, pb_z, ip_12, id_13, id_14, ks0_5, \
                         ks1_5, kp_17, kp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_17[k];

        t_43[k] = pa_z[k] * id_13[k];

        t_44[k] = pa_z[k] * id_14[k];

        t_45[k] = f_16 * ip_12[k]
                  + pb_y[k] * kp_18[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_x, pa_y, pa_z, hd0_5, hd0_13, hd1_10, hd1_22, \
                         ip_11, id_15, id_16, id_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_4 * ip_11[k]
                  + pa_z[k] * id_15[k];

        t_47[k] = f_5 * hd0_5[k]
                  - f_6 * hd1_10[k]
                  + pa_y[k] * id_16[k];

        t_48[k] = f_11 * hd0_13[k]
                  - f_12 * hd1_22[k]
                  + pa_x[k] * id_24[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_x, pa_y, pb_y, hd0_14, hd1_23, ip_13, \
                         ip_14, id_17, id_18, id_25, kp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_4 * ip_13[k]
                  + pb_y[k] * kp_19[k];

        t_50[k] = f_11 * hd0_14[k]
                  - f_12 * hd1_23[k]
                  + pa_x[k] * id_25[k];

        t_51[k] = pa_y[k] * id_17[k];

        t_52[k] = f_4 * ip_14[k]
                  + pa_y[k] * id_18[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_y, pa_z, pb_x, pb_y, hd0_5, hd1_10, ip_15, \
                         ip_22, id_17, id_19, kp_20, kp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_10 * ip_15[k]
                  + pb_y[k] * kp_20[k];

        t_54[k] = pa_y[k] * id_19[k];

        t_55[k] = f_14 * hd0_5[k]
                  - f_15 * hd1_10[k]
                  + pa_z[k] * id_17[k];

        t_56[k] = f_16 * ip_22[k]
                  + pb_x[k] * kp_22[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_x, pa_y, pb_y, hd0_7, hd0_15, hd1_13, hd1_25, \
                         id_20, id_29, ks0_6, ks1_6, kp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_21[k];

        t_58[k] = f_11 * hd0_15[k]
                  - f_12 * hd1_25[k]
                  + pa_x[k] * id_29[k];

        t_59[k] = f_8 * hd0_7[k]
                  - f_9 * hd1_13[k]
                  + pa_y[k] * id_20[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_x, pb_x, pb_z, hd0_16, hd1_27, ip_23, id_31, \
                         ks0_7, ks1_7, kp_23, kp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_4 * ip_23[k]
                  + pb_x[k] * kp_23[k];

        t_61[k] = f_5 * hd0_16[k]
                  - f_6 * hd1_27[k]
                  + pa_x[k] * id_31[k];

        t_62[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_24[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pb_y, ip_17, ip_18, id_20, id_21, \
                         id_22, kp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_z[k] * id_20[k];

        t_64[k] = pa_z[k] * id_21[k];

        t_65[k] = f_13 * ip_18[k]
                  + pb_y[k] * kp_25[k];

        t_66[k] = f_4 * ip_17[k]
                  + pa_z[k] * id_22[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pa_x, pa_y, pb_y, hd0_9, hd0_18, hd1_16, hd1_33, \
                         ip_19, id_23, id_32, kp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_11 * hd0_9[k]
                  - f_12 * hd1_16[k]
                  + pa_y[k] * id_23[k];

        t_68[k] = f_5 * hd0_18[k]
                  - f_6 * hd1_33[k]
                  + pa_x[k] * id_32[k];

        t_69[k] = f_16 * ip_19[k]
                  + pb_y[k] * kp_26[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_x, pa_y, hd0_10, hd0_19, hd0_20, hd1_17, hd1_35, \
                         hd1_37, id_26, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_5 * hd0_19[k]
                  - f_6 * hd1_35[k]
                  + pa_x[k] * id_33[k];

        t_71[k] = f_5 * hd0_10[k]
                  - f_6 * hd1_17[k]
                  + pa_y[k] * id_26[k];

        t_72[k] = f_5 * hd0_20[k]
                  - f_6 * hd1_37[k]
                  + pa_x[k] * id_34[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pa_y, pb_y, hd0_21, hd1_39, ip_20, \
                         ip_21, id_27, id_28, id_35, kp_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_4 * ip_20[k]
                  + pb_y[k] * kp_27[k];

        t_74[k] = f_5 * hd0_21[k]
                  - f_6 * hd1_39[k]
                  + pa_x[k] * id_35[k];

        t_75[k] = pa_y[k] * id_27[k];

        t_76[k] = f_4 * ip_21[k]
                  + pa_y[k] * id_28[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pa_z, pb_x, pb_y, hd0_10, hd1_17, \
                         ip_22, ip_24, id_27, id_29, kp_28, kp_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_10 * ip_22[k]
                  + pb_y[k] * kp_28[k];

        t_78[k] = pa_y[k] * id_29[k];

        t_79[k] = f_8 * hd0_10[k]
                  - f_9 * hd1_17[k]
                  + pa_z[k] * id_27[k];

        t_80[k] = f_4 * ip_24[k]
                  + pb_x[k] * kp_30[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_x, pb_y, hd0_23, hd1_45, ip_25, id_37, id_38, \
                         ks0_8, ks1_8, kp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_29[k];

        t_82[k] = f_5 * hd0_23[k]
                  - f_6 * hd1_45[k]
                  + pa_x[k] * id_37[k];

        t_83[k] = f_4 * ip_25[k]
                  + pa_x[k] * id_38[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, pa_x, pa_z, pb_x, ip_26, id_30, \
                         id_39, id_40, id_42, id_43, kp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_10 * ip_26[k]
                  + pb_x[k] * kp_31[k];

        t_85[k] = pa_x[k] * id_39[k];

        t_86[k] = pa_x[k] * id_40[k];

        t_87[k] = pa_z[k] * id_30[k];

        t_88[k] = pa_x[k] * id_42[k];

        t_89[k] = pa_x[k] * id_43[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, pa_x, ip_29, ip_31, id_44, id_45, \
                         id_46, id_47, id_48, id_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_4 * ip_29[k]
                  + pa_x[k] * id_44[k];

        t_91[k] = pa_x[k] * id_45[k];

        t_92[k] = pa_x[k] * id_46[k];

        t_93[k] = pa_x[k] * id_47[k];

        t_94[k] = f_4 * ip_31[k]
                  + pa_x[k] * id_48[k];

        t_95[k] = pa_x[k] * id_49[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, t_101, pa_x, ip_33, id_50, id_51, \
                         id_52, id_53, id_54, id_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_x[k] * id_50[k];

        t_97[k] = pa_x[k] * id_51[k];

        t_98[k] = f_4 * ip_33[k]
                  + pa_x[k] * id_52[k];

        t_99[k] = pa_x[k] * id_53[k];

        t_100[k] = pa_x[k] * id_54[k];

        t_101[k] = pa_x[k] * id_55[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, pa_x, pa_y, pb_x, ip_36, ip_38, \
                         id_36, id_56, id_57, id_59, kp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = pa_y[k] * id_36[k];

        t_103[k] = pa_x[k] * id_56[k];

        t_104[k] = pa_x[k] * id_57[k];

        t_105[k] = f_4 * ip_36[k]
                   + pa_x[k] * id_59[k];

        t_106[k] = f_10 * ip_38[k]
                   + pb_x[k] * kp_32[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, pa_x, pb_x, pb_y, pb_z, ip_26, \
                         id_60, id_61, ks0_9, ks1_9, kp_33, kp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pa_x[k] * id_60[k];

        t_108[k] = pa_x[k] * id_61[k];

        t_109[k] = f_1 * ks0_9[k]
                   - f_2 * ks1_9[k]
                   + pb_x[k] * kp_33[k];

        t_110[k] = f_0 * ip_26[k]
                   + f_1 * ks0_9[k]
                   - f_2 * ks1_9[k]
                   + pb_y[k] * kp_34[k];

        t_111[k] = pb_z[k] * kp_34[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pa_z, pb_y, pb_z, ip_28, id_38, id_39, \
                         ks0_9, ks1_9, kp_35, kp_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_1 * ks0_9[k]
                   - f_2 * ks1_9[k]
                   + pb_z[k] * kp_35[k];

        t_113[k] = pa_z[k] * id_38[k];

        t_114[k] = pa_z[k] * id_39[k];

        t_115[k] = f_3 * ip_28[k]
                   + pb_y[k] * kp_36[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pa_z, pb_x, hd0_16, hd1_27, ip_27, id_40, id_41, \
                         ks0_10, ks1_10, kp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_4 * ip_27[k]
                   + pa_z[k] * id_40[k];

        t_117[k] = f_1 * ks0_10[k]
                   - f_2 * ks1_10[k]
                   + pb_x[k] * kp_37[k];

        t_118[k] = f_5 * hd0_16[k]
                   - f_6 * hd1_27[k]
                   + pa_z[k] * id_41[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pa_y, pb_x, pb_y, hd0_19, hd1_35, ip_30, id_47, \
                         ks0_11, ks1_11, kp_38, kp_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_7 * ip_30[k]
                   + pb_y[k] * kp_38[k];

        t_120[k] = f_8 * hd0_19[k]
                   - f_9 * hd1_35[k]
                   + pa_y[k] * id_47[k];

        t_121[k] = f_1 * ks0_11[k]
                   - f_2 * ks1_11[k]
                   + pb_x[k] * kp_39[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pa_y, pa_z, pb_y, hd0_17, hd0_21, hd1_29, \
                         hd1_39, ip_32, id_45, id_51, kp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_11 * hd0_17[k]
                   - f_12 * hd1_29[k]
                   + pa_z[k] * id_45[k];

        t_123[k] = f_13 * ip_32[k]
                   + pb_y[k] * kp_40[k];

        t_124[k] = f_14 * hd0_21[k]
                   - f_15 * hd1_39[k]
                   + pa_y[k] * id_51[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pa_z, pb_x, pb_y, hd0_18, hd1_33, ip_34, id_49, \
                         ks0_12, ks1_12, kp_41, kp_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_1 * ks0_12[k]
                   - f_2 * ks1_12[k]
                   + pb_x[k] * kp_41[k];

        t_126[k] = f_14 * hd0_18[k]
                   - f_15 * hd1_33[k]
                   + pa_z[k] * id_49[k];

        t_127[k] = f_16 * ip_34[k]
                   + pb_y[k] * kp_42[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, pa_y, pa_z, pb_x, hd0_20, hd0_22, hd1_37, \
                         hd1_42, id_53, id_55, ks0_13, ks1_13, kp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_11 * hd0_22[k]
                   - f_12 * hd1_42[k]
                   + pa_y[k] * id_55[k];

        t_129[k] = f_1 * ks0_13[k]
                   - f_2 * ks1_13[k]
                   + pb_x[k] * kp_43[k];

        t_130[k] = f_8 * hd0_20[k]
                   - f_9 * hd1_37[k]
                   + pa_z[k] * id_53[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_y, pb_y, hd0_23, hd1_45, ip_35, ip_37, \
                         id_58, id_59, id_60, kp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_4 * ip_35[k]
                   + pb_y[k] * kp_44[k];

        t_132[k] = f_5 * hd0_23[k]
                   - f_6 * hd1_45[k]
                   + pa_y[k] * id_58[k];

        t_133[k] = pa_y[k] * id_59[k];

        t_134[k] = f_4 * ip_37[k]
                   + pa_y[k] * id_60[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, pa_y, pb_x, pb_y, ip_38, id_61, \
                         ks0_14, ks1_14, kp_45, kp_46, kp_47, kp_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_10 * ip_38[k]
                   + pb_y[k] * kp_45[k];

        t_136[k] = pa_y[k] * id_61[k];

        t_137[k] = f_1 * ks0_14[k]
                   - f_2 * ks1_14[k]
                   + pb_x[k] * kp_46[k];

        t_138[k] = f_1 * ks0_14[k]
                   - f_2 * ks1_14[k]
                   + pb_y[k] * kp_47[k];

        t_139[k] = pb_y[k] * kp_48[k];
    }

#pragma omp simd aligned(t_140, pb_z, ip_38, ks0_14, ks1_14, kp_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * ip_38[k]
                   + f_1 * ks0_14[k]
                   - f_2 * ks1_14[k]
                   + pb_z[k] * kp_48[k];
    }
}

auto
compute_prim_kd_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hd0, const size_t hd1,
                                     const size_t ip, const size_t id, const size_t ks0,
                                     const size_t ks1, const size_t kp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / alpha;
    const auto f_6 = 2.0 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 1.5 / alpha;
    const auto f_10 = 1.5 * beta / (alpha * p);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_17 = buffer.data(ip + 17);

    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_9 = buffer.data(ks0 + 9);
    const auto *ks0_10 = buffer.data(ks0 + 10);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_9 = buffer.data(ks1 + 9);
    const auto *ks1_10 = buffer.data(ks1 + 10);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, ip_0, ks0_0, ks1_0, kp_0, kp_1, \
                         kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, hd0_0, hd0_4, hd1_0, hd1_4, id_3, \
                         id_6, ks0_1, ks1_1, kp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_y[k] * id_3[k];

        t_4[k] = f_5 * hd0_4[k]
                 - f_6 * hd1_4[k]
                 + pa_x[k] * id_6[k];

        t_5[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_1[k]
                 + pb_z[k] * kp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, hd0_0, hd0_6, hd1_0, hd1_6, id_4, \
                         id_10, ks0_2, ks1_2, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_z[k] * id_4[k];

        t_7[k] = f_1 * ks0_2[k]
                 - f_2 * ks1_2[k]
                 + pb_y[k] * kp_4[k];

        t_8[k] = f_5 * hd0_6[k]
                 - f_6 * hd1_6[k]
                 + pa_x[k] * id_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, hd0_1, hd0_8, hd1_1, hd1_8, id_5, \
                         id_12, ks0_3, ks1_3, kp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * hd0_1[k]
                 - f_8 * hd1_1[k]
                 + pa_y[k] * id_5[k];

        t_10[k] = f_9 * hd0_8[k]
                  - f_10 * hd1_8[k]
                  + pa_x[k] * id_12[k];

        t_11[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, hd0_2, hd0_10, hd1_2, hd1_10, \
                         id_8, id_16, ks0_4, ks1_4, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * hd0_2[k]
                  - f_8 * hd1_2[k]
                  + pa_z[k] * id_8[k];

        t_13[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];

        t_14[k] = f_9 * hd0_10[k]
                  - f_10 * hd1_10[k]
                  + pa_x[k] * id_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, hd0_3, hd0_11, hd1_3, hd1_11, \
                         id_11, id_18, ks0_5, ks1_5, kp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_9 * hd0_3[k]
                  - f_10 * hd1_3[k]
                  + pa_y[k] * id_11[k];

        t_16[k] = f_7 * hd0_11[k]
                  - f_8 * hd1_11[k]
                  + pa_x[k] * id_18[k];

        t_17[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, hd0_5, hd0_12, hd1_5, hd1_12, \
                         id_14, id_22, ks0_6, ks1_6, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_9 * hd0_5[k]
                  - f_10 * hd1_5[k]
                  + pa_z[k] * id_14[k];

        t_19[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];

        t_20[k] = f_7 * hd0_12[k]
                  - f_8 * hd1_12[k]
                  + pa_x[k] * id_22[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_z, hd0_7, hd0_13, hd1_7, hd1_13, \
                         id_17, id_23, ks0_7, ks1_7, kp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * hd0_7[k]
                  - f_6 * hd1_7[k]
                  + pa_y[k] * id_17[k];

        t_22[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_13[k]
                  + pa_x[k] * id_23[k];

        t_23[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_z, pb_y, hd0_9, hd0_20, hd1_9, hd1_20, \
                         id_20, id_24, ks0_8, ks1_8, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * hd0_9[k]
                  - f_6 * hd1_9[k]
                  + pa_z[k] * id_20[k];

        t_25[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];

        t_26[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pa_x[k] * id_24[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_x, pb_y, pb_z, ip_10, ks0_9, ks0_10, \
                         ks1_9, ks1_10, kp_11, kp_12, kp_13, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_x[k] * kp_11[k];

        t_28[k] = f_0 * ip_10[k]
                  + f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_y[k] * kp_12[k];

        t_29[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_z[k] * kp_13[k];

        t_30[k] = f_1 * ks0_10[k]
                  - f_2 * ks1_10[k]
                  + pb_x[k] * kp_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pa_z, pb_x, hd0_13, hd0_16, hd1_13, hd1_16, \
                         id_28, id_31, ks0_11, ks1_11, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_13[k]
                  + pa_z[k] * id_28[k];

        t_32[k] = f_5 * hd0_16[k]
                  - f_6 * hd1_16[k]
                  + pa_y[k] * id_31[k];

        t_33[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_15[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pa_z, pb_x, hd0_14, hd0_18, hd1_14, hd1_18, \
                         id_30, id_34, ks0_12, ks1_12, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_7 * hd0_14[k]
                  - f_8 * hd1_14[k]
                  + pa_z[k] * id_30[k];

        t_35[k] = f_9 * hd0_18[k]
                  - f_10 * hd1_18[k]
                  + pa_y[k] * id_34[k];

        t_36[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_16[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, hd0_15, hd0_19, hd1_15, hd1_19, \
                         id_33, id_37, ks0_13, ks1_13, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * hd0_15[k]
                  - f_10 * hd1_15[k]
                  + pa_z[k] * id_33[k];

        t_38[k] = f_7 * hd0_19[k]
                  - f_8 * hd1_19[k]
                  + pa_y[k] * id_37[k];

        t_39[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_17[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_x, hd0_17, hd0_20, hd1_17, hd1_20, \
                         id_36, id_38, ks0_14, ks1_14, kp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * hd0_17[k]
                  - f_6 * hd1_17[k]
                  + pa_z[k] * id_36[k];

        t_41[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pa_y[k] * id_38[k];

        t_42[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_18[k];
    }

#pragma omp simd aligned(t_43, t_44, pb_y, pb_z, ip_17, ks0_14, ks1_14, kp_19, \
                         kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_y[k] * kp_19[k];

        t_44[k] = f_0 * ip_17[k]
                  + f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hd0, const size_t hd1,
                                     const size_t ip, const size_t id, const size_t ks0,
                                     const size_t ks1, const size_t kp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / alpha;
    const auto f_6 = 2.0 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 1.5 / alpha;
    const auto f_10 = 1.5 * beta / (alpha * p);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_26 = buffer.data(hd1 + 26);
    const auto *hd1_28 = buffer.data(hd1 + 28);
    const auto *hd1_30 = buffer.data(hd1 + 30);
    const auto *hd1_36 = buffer.data(hd1 + 36);
    const auto *hd1_38 = buffer.data(hd1 + 38);
    const auto *hd1_41 = buffer.data(hd1 + 41);
    const auto *hd1_44 = buffer.data(hd1 + 44);
    const auto *hd1_45 = buffer.data(hd1 + 45);
    const auto *hd1_47 = buffer.data(hd1 + 47);
    const auto *hd1_48 = buffer.data(hd1 + 48);
    const auto *hd1_50 = buffer.data(hd1 + 50);
    const auto *hd1_53 = buffer.data(hd1 + 53);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_17 = buffer.data(ip + 17);

    const auto *id_3 = buffer.data(id + 3);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_62 = buffer.data(id + 62);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_65 = buffer.data(id + 65);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_9 = buffer.data(ks0 + 9);
    const auto *ks0_10 = buffer.data(ks0 + 10);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_9 = buffer.data(ks1 + 9);
    const auto *ks1_10 = buffer.data(ks1 + 10);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, ip_0, ks0_0, ks1_0, kp_0, kp_1, \
                         kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, hd0_0, hd0_4, hd1_0, hd1_10, id_3, \
                         id_10, ks0_1, ks1_1, kp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_y[k] * id_3[k];

        t_4[k] = f_5 * hd0_4[k]
                 - f_6 * hd1_10[k]
                 + pa_x[k] * id_10[k];

        t_5[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_1[k]
                 + pb_z[k] * kp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, hd0_0, hd0_6, hd1_0, hd1_16, id_6, \
                         id_16, ks0_2, ks1_2, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_z[k] * id_6[k];

        t_7[k] = f_1 * ks0_2[k]
                 - f_2 * ks1_2[k]
                 + pb_y[k] * kp_4[k];

        t_8[k] = f_5 * hd0_6[k]
                 - f_6 * hd1_16[k]
                 + pa_x[k] * id_16[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, hd0_1, hd0_8, hd1_3, hd1_18, id_9, \
                         id_18, ks0_3, ks1_3, kp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * hd0_1[k]
                 - f_8 * hd1_3[k]
                 + pa_y[k] * id_9[k];

        t_10[k] = f_9 * hd0_8[k]
                  - f_10 * hd1_18[k]
                  + pa_x[k] * id_18[k];

        t_11[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, hd0_2, hd0_10, hd1_6, hd1_28, \
                         id_14, id_27, ks0_4, ks1_4, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * hd0_2[k]
                  - f_8 * hd1_6[k]
                  + pa_z[k] * id_14[k];

        t_13[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];

        t_14[k] = f_9 * hd0_10[k]
                  - f_10 * hd1_28[k]
                  + pa_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, hd0_3, hd0_11, hd1_9, hd1_30, \
                         id_17, id_29, ks0_5, ks1_5, kp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_9 * hd0_3[k]
                  - f_10 * hd1_9[k]
                  + pa_y[k] * id_17[k];

        t_16[k] = f_7 * hd0_11[k]
                  - f_8 * hd1_30[k]
                  + pa_x[k] * id_29[k];

        t_17[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, hd0_5, hd0_12, hd1_14, hd1_36, \
                         id_25, id_41, ks0_6, ks1_6, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_9 * hd0_5[k]
                  - f_10 * hd1_14[k]
                  + pa_z[k] * id_25[k];

        t_19[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];

        t_20[k] = f_7 * hd0_12[k]
                  - f_8 * hd1_36[k]
                  + pa_x[k] * id_41[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_z, hd0_7, hd0_13, hd1_17, hd1_38, \
                         id_28, id_43, ks0_7, ks1_7, kp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * hd0_7[k]
                  - f_6 * hd1_17[k]
                  + pa_y[k] * id_28[k];

        t_22[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_38[k]
                  + pa_x[k] * id_43[k];

        t_23[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_z, pb_y, hd0_9, hd0_20, hd1_26, hd1_53, \
                         id_39, id_48, ks0_8, ks1_8, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * hd0_9[k]
                  - f_6 * hd1_26[k]
                  + pa_z[k] * id_39[k];

        t_25[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];

        t_26[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_53[k]
                  + pa_x[k] * id_48[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_x, pb_y, pb_z, ip_10, ks0_9, ks0_10, \
                         ks1_9, ks1_10, kp_11, kp_12, kp_13, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_x[k] * kp_11[k];

        t_28[k] = f_0 * ip_10[k]
                  + f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_y[k] * kp_12[k];

        t_29[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_z[k] * kp_13[k];

        t_30[k] = f_1 * ks0_10[k]
                  - f_2 * ks1_10[k]
                  + pb_x[k] * kp_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pa_z, pb_x, hd0_13, hd0_16, hd1_38, hd1_45, \
                         id_53, id_57, ks0_11, ks1_11, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_38[k]
                  + pa_z[k] * id_53[k];

        t_32[k] = f_5 * hd0_16[k]
                  - f_6 * hd1_45[k]
                  + pa_y[k] * id_57[k];

        t_33[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_15[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pa_z, pb_x, hd0_14, hd0_18, hd1_41, hd1_48, \
                         id_56, id_60, ks0_12, ks1_12, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_7 * hd0_14[k]
                  - f_8 * hd1_41[k]
                  + pa_z[k] * id_56[k];

        t_35[k] = f_9 * hd0_18[k]
                  - f_10 * hd1_48[k]
                  + pa_y[k] * id_60[k];

        t_36[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_16[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, hd0_15, hd0_19, hd1_44, hd1_50, \
                         id_59, id_63, ks0_13, ks1_13, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * hd0_15[k]
                  - f_10 * hd1_44[k]
                  + pa_z[k] * id_59[k];

        t_38[k] = f_7 * hd0_19[k]
                  - f_8 * hd1_50[k]
                  + pa_y[k] * id_63[k];

        t_39[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_17[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_x, hd0_17, hd0_20, hd1_47, hd1_53, \
                         id_62, id_65, ks0_14, ks1_14, kp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * hd0_17[k]
                  - f_6 * hd1_47[k]
                  + pa_z[k] * id_62[k];

        t_41[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_53[k]
                  + pa_y[k] * id_65[k];

        t_42[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_18[k];
    }

#pragma omp simd aligned(t_43, t_44, pb_y, pb_z, ip_17, ks0_14, ks1_14, kp_19, \
                         kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_y[k] * kp_19[k];

        t_44[k] = f_0 * ip_17[k]
                  + f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hd0, const size_t hd1,
                                     const size_t ip, const size_t id, const size_t ks0,
                                     const size_t ks1, const size_t kp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / alpha;
    const auto f_11 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_23 = buffer.data(hd0 + 23);
    const auto *hd0_26 = buffer.data(hd0 + 26);
    const auto *hd0_28 = buffer.data(hd0 + 28);
    const auto *hd0_30 = buffer.data(hd0 + 30);
    const auto *hd0_33 = buffer.data(hd0 + 33);
    const auto *hd0_34 = buffer.data(hd0 + 34);
    const auto *hd0_36 = buffer.data(hd0 + 36);
    const auto *hd0_38 = buffer.data(hd0 + 38);
    const auto *hd0_41 = buffer.data(hd0 + 41);
    const auto *hd0_44 = buffer.data(hd0 + 44);
    const auto *hd0_45 = buffer.data(hd0 + 45);
    const auto *hd0_47 = buffer.data(hd0 + 47);
    const auto *hd0_48 = buffer.data(hd0 + 48);
    const auto *hd0_50 = buffer.data(hd0 + 50);
    const auto *hd0_53 = buffer.data(hd0 + 53);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_29 = buffer.data(hd1 + 29);
    const auto *hd1_32 = buffer.data(hd1 + 32);
    const auto *hd1_33 = buffer.data(hd1 + 33);
    const auto *hd1_35 = buffer.data(hd1 + 35);
    const auto *hd1_36 = buffer.data(hd1 + 36);
    const auto *hd1_38 = buffer.data(hd1 + 38);
    const auto *hd1_41 = buffer.data(hd1 + 41);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);

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
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_9 = buffer.data(ks0 + 9);
    const auto *ks0_10 = buffer.data(ks0 + 10);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_9 = buffer.data(ks1 + 9);
    const auto *ks1_10 = buffer.data(ks1 + 10);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, ip_0, id_0, ks0_0, ks1_0, \
                         kp_0, kp_1, kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];

        t_3[k] = pa_y[k] * id_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, hd0_0, hd1_0, ip_1, ip_2, \
                         id_0, id_1, id_2, id_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * ip_1[k]
                 + pa_y[k] * id_1[k];

        t_5[k] = pa_y[k] * id_2[k];

        t_6[k] = pa_z[k] * id_0[k];

        t_7[k] = pa_z[k] * id_1[k];

        t_8[k] = f_3 * ip_2[k]
                 + pa_z[k] * id_2[k];

        t_9[k] = f_4 * hd0_0[k]
                 - f_5 * hd1_0[k]
                 + pa_y[k] * id_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, pa_z, pb_z, hd0_10, hd1_8, id_4, \
                         id_6, id_8, ks0_1, ks1_1, kp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_6 * hd0_10[k]
                  - f_7 * hd1_8[k]
                  + pa_x[k] * id_8[k];

        t_11[k] = f_1 * ks0_1[k]
                  - f_2 * ks1_1[k]
                  + pb_z[k] * kp_3[k];

        t_12[k] = pa_z[k] * id_4[k];

        t_13[k] = pa_y[k] * id_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_y, hd0_0, hd0_16, hd1_0, hd1_12, \
                         id_5, id_12, ks0_2, ks1_2, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_4 * hd0_0[k]
                  - f_5 * hd1_0[k]
                  + pa_z[k] * id_5[k];

        t_15[k] = f_1 * ks0_2[k]
                  - f_2 * ks1_2[k]
                  + pb_y[k] * kp_4[k];

        t_16[k] = f_6 * hd0_16[k]
                  - f_7 * hd1_12[k]
                  + pa_x[k] * id_12[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, pb_z, hd0_3, hd0_18, hd1_3, hd1_14, \
                         id_7, id_14, ks0_3, ks1_3, kp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * hd0_3[k]
                  - f_9 * hd1_3[k]
                  + pa_y[k] * id_7[k];

        t_18[k] = f_10 * hd0_18[k]
                  - f_11 * hd1_14[k]
                  + pa_x[k] * id_14[k];

        t_19[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pa_z, ip_3, ip_4, id_7, id_8, \
                         id_9, id_11, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * id_7[k];

        t_21[k] = pa_z[k] * id_8[k];

        t_22[k] = f_3 * ip_3[k]
                  + pa_z[k] * id_9[k];

        t_23[k] = f_3 * ip_4[k]
                  + pa_y[k] * id_11[k];

        t_24[k] = pa_y[k] * id_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_x, pa_z, pb_y, hd0_6, hd0_28, hd1_5, hd1_19, \
                         id_10, id_19, ks0_4, ks1_4, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * hd0_6[k]
                  - f_9 * hd1_5[k]
                  + pa_z[k] * id_10[k];

        t_26[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];

        t_27[k] = f_10 * hd0_28[k]
                  - f_11 * hd1_19[k]
                  + pa_x[k] * id_19[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_x, pa_y, pb_z, hd0_9, hd0_30, hd1_7, hd1_21, \
                         id_13, id_21, ks0_5, ks1_5, kp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_10 * hd0_9[k]
                  - f_11 * hd1_7[k]
                  + pa_y[k] * id_13[k];

        t_29[k] = f_8 * hd0_30[k]
                  - f_9 * hd1_21[k]
                  + pa_x[k] * id_21[k];

        t_30[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pa_z, hd0_14, hd1_10, ip_5, id_13, \
                         id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pa_z[k] * id_13[k];

        t_32[k] = pa_z[k] * id_14[k];

        t_33[k] = f_3 * ip_5[k]
                  + pa_z[k] * id_15[k];

        t_34[k] = f_4 * hd0_14[k]
                  - f_5 * hd1_10[k]
                  + pa_y[k] * id_16[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pa_y, hd0_33, hd0_34, hd1_22, hd1_23, \
                         ip_6, id_18, id_19, id_24, id_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * hd0_33[k]
                  - f_9 * hd1_22[k]
                  + pa_x[k] * id_24[k];

        t_36[k] = f_8 * hd0_34[k]
                  - f_9 * hd1_23[k]
                  + pa_x[k] * id_25[k];

        t_37[k] = f_3 * ip_6[k]
                  + pa_y[k] * id_18[k];

        t_38[k] = pa_y[k] * id_19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_x, pa_z, pb_y, hd0_14, hd0_36, hd1_10, hd1_25, \
                         id_17, id_29, ks0_6, ks1_6, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * hd0_14[k]
                  - f_11 * hd1_10[k]
                  + pa_z[k] * id_17[k];

        t_40[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];

        t_41[k] = f_8 * hd0_36[k]
                  - f_9 * hd1_25[k]
                  + pa_x[k] * id_29[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_x, pa_y, pb_z, hd0_17, hd0_38, hd1_13, hd1_27, \
                         id_20, id_31, ks0_7, ks1_7, kp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_6 * hd0_17[k]
                  - f_7 * hd1_13[k]
                  + pa_y[k] * id_20[k];

        t_43[k] = f_4 * hd0_38[k]
                  - f_5 * hd1_27[k]
                  + pa_x[k] * id_31[k];

        t_44[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_y, pa_z, hd0_23, hd1_16, ip_7, id_20, \
                         id_21, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_z[k] * id_20[k];

        t_46[k] = pa_z[k] * id_21[k];

        t_47[k] = f_3 * ip_7[k]
                  + pa_z[k] * id_22[k];

        t_48[k] = f_8 * hd0_23[k]
                  - f_9 * hd1_16[k]
                  + pa_y[k] * id_23[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_x, pa_y, hd0_26, hd0_44, hd0_45, hd1_17, hd1_32, \
                         hd1_33, id_26, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_4 * hd0_44[k]
                  - f_5 * hd1_32[k]
                  + pa_x[k] * id_32[k];

        t_50[k] = f_4 * hd0_45[k]
                  - f_5 * hd1_33[k]
                  + pa_x[k] * id_33[k];

        t_51[k] = f_4 * hd0_26[k]
                  - f_5 * hd1_17[k]
                  + pa_y[k] * id_26[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pa_y, hd0_47, hd0_48, hd1_35, hd1_36, \
                         ip_8, id_28, id_29, id_34, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_4 * hd0_47[k]
                  - f_5 * hd1_35[k]
                  + pa_x[k] * id_34[k];

        t_53[k] = f_4 * hd0_48[k]
                  - f_5 * hd1_36[k]
                  + pa_x[k] * id_35[k];

        t_54[k] = f_3 * ip_8[k]
                  + pa_y[k] * id_28[k];

        t_55[k] = pa_y[k] * id_29[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pa_x, pa_z, pb_y, hd0_26, hd0_53, hd1_17, hd1_41, \
                         id_27, id_37, ks0_8, ks1_8, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_6 * hd0_26[k]
                  - f_7 * hd1_17[k]
                  + pa_z[k] * id_27[k];

        t_57[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];

        t_58[k] = f_4 * hd0_53[k]
                  - f_5 * hd1_41[k]
                  + pa_x[k] * id_37[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, t_63, pa_x, pa_z, ip_9, ip_12, ip_13, id_30, \
                         id_38, id_39, id_43, id_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_3 * ip_9[k]
                  + pa_x[k] * id_38[k];

        t_60[k] = pa_x[k] * id_39[k];

        t_61[k] = pa_z[k] * id_30[k];

        t_62[k] = f_3 * ip_12[k]
                  + pa_x[k] * id_43[k];

        t_63[k] = f_3 * ip_13[k]
                  + pa_x[k] * id_46[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, ip_14, ip_15, id_49, id_54, \
                         id_56, ks0_9, ks1_9, kp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * ip_14[k]
                  + pa_x[k] * id_49[k];

        t_65[k] = f_3 * ip_15[k]
                  + pa_x[k] * id_54[k];

        t_66[k] = pa_x[k] * id_56[k];

        t_67[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_x[k] * kp_11[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_z, pb_y, pb_z, ip_10, id_38, id_39, ks0_9, \
                         ks1_9, kp_12, kp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * ip_10[k]
                  + f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_y[k] * kp_12[k];

        t_69[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_z[k] * kp_13[k];

        t_70[k] = pa_z[k] * id_38[k];

        t_71[k] = pa_z[k] * id_39[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_z, pb_x, hd0_38, hd1_27, ip_11, id_40, id_41, \
                         ks0_10, ks1_10, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_3 * ip_11[k]
                  + pa_z[k] * id_40[k];

        t_73[k] = f_1 * ks0_10[k]
                  - f_2 * ks1_10[k]
                  + pb_x[k] * kp_14[k];

        t_74[k] = f_4 * hd0_38[k]
                  - f_5 * hd1_27[k]
                  + pa_z[k] * id_41[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_y, pa_z, pb_x, hd0_41, hd0_45, hd1_29, hd1_33, \
                         id_44, id_45, ks0_11, ks1_11, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_6 * hd0_45[k]
                  - f_7 * hd1_33[k]
                  + pa_y[k] * id_45[k];

        t_76[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_15[k];

        t_77[k] = f_8 * hd0_41[k]
                  - f_9 * hd1_29[k]
                  + pa_z[k] * id_44[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pa_y, pa_z, pb_x, hd0_44, hd0_48, hd1_32, hd1_36, \
                         id_47, id_48, ks0_12, ks1_12, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * hd0_48[k]
                  - f_11 * hd1_36[k]
                  + pa_y[k] * id_48[k];

        t_79[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_16[k];

        t_80[k] = f_10 * hd0_44[k]
                  - f_11 * hd1_32[k]
                  + pa_z[k] * id_47[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_y, pa_z, pb_x, hd0_47, hd0_50, hd1_35, hd1_38, \
                         id_50, id_51, ks0_13, ks1_13, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_8 * hd0_50[k]
                  - f_9 * hd1_38[k]
                  + pa_y[k] * id_51[k];

        t_82[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_17[k];

        t_83[k] = f_6 * hd0_47[k]
                  - f_7 * hd1_35[k]
                  + pa_z[k] * id_50[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_y, pb_x, hd0_53, hd1_41, ip_16, id_53, \
                         id_55, id_56, ks0_14, ks1_14, kp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * hd0_53[k]
                  - f_5 * hd1_41[k]
                  + pa_y[k] * id_53[k];

        t_85[k] = f_3 * ip_16[k]
                  + pa_y[k] * id_55[k];

        t_86[k] = pa_y[k] * id_56[k];

        t_87[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_18[k];
    }

#pragma omp simd aligned(t_88, t_89, pb_y, pb_z, ip_17, ks0_14, ks1_14, kp_19, \
                         kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_y[k] * kp_19[k];

        t_89[k] = f_0 * ip_17[k]
                  + f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hd0, const size_t hd1,
                                     const size_t ip, const size_t id, const size_t ks0,
                                     const size_t ks1, const size_t kp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / alpha;
    const auto f_6 = 2.0 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 1.5 / alpha;
    const auto f_10 = 1.5 * beta / (alpha * p);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_28 = buffer.data(hd1 + 28);
    const auto *hd1_29 = buffer.data(hd1 + 29);
    const auto *hd1_32 = buffer.data(hd1 + 32);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_17 = buffer.data(ip + 17);

    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_9 = buffer.data(ks0 + 9);
    const auto *ks0_10 = buffer.data(ks0 + 10);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_9 = buffer.data(ks1 + 9);
    const auto *ks1_10 = buffer.data(ks1 + 10);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, ip_0, ks0_0, ks1_0, kp_0, kp_1, \
                         kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, hd0_0, hd0_4, hd1_0, hd1_6, id_3, \
                         id_6, ks0_1, ks1_1, kp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_y[k] * id_3[k];

        t_4[k] = f_5 * hd0_4[k]
                 - f_6 * hd1_6[k]
                 + pa_x[k] * id_6[k];

        t_5[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_1[k]
                 + pb_z[k] * kp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, hd0_0, hd0_6, hd1_0, hd1_10, id_4, \
                         id_10, ks0_2, ks1_2, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_z[k] * id_4[k];

        t_7[k] = f_1 * ks0_2[k]
                 - f_2 * ks1_2[k]
                 + pb_y[k] * kp_4[k];

        t_8[k] = f_5 * hd0_6[k]
                 - f_6 * hd1_10[k]
                 + pa_x[k] * id_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, hd0_1, hd0_8, hd1_3, hd1_12, id_5, \
                         id_12, ks0_3, ks1_3, kp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * hd0_1[k]
                 - f_8 * hd1_3[k]
                 + pa_y[k] * id_5[k];

        t_10[k] = f_9 * hd0_8[k]
                  - f_10 * hd1_12[k]
                  + pa_x[k] * id_12[k];

        t_11[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, hd0_2, hd0_10, hd1_4, hd1_16, \
                         id_8, id_16, ks0_4, ks1_4, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * hd0_2[k]
                  - f_8 * hd1_4[k]
                  + pa_z[k] * id_8[k];

        t_13[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];

        t_14[k] = f_9 * hd0_10[k]
                  - f_10 * hd1_16[k]
                  + pa_x[k] * id_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, hd0_3, hd0_11, hd1_5, hd1_17, \
                         id_11, id_18, ks0_5, ks1_5, kp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_9 * hd0_3[k]
                  - f_10 * hd1_5[k]
                  + pa_y[k] * id_11[k];

        t_16[k] = f_7 * hd0_11[k]
                  - f_8 * hd1_17[k]
                  + pa_x[k] * id_18[k];

        t_17[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, hd0_5, hd0_12, hd1_8, hd1_18, \
                         id_14, id_22, ks0_6, ks1_6, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_9 * hd0_5[k]
                  - f_10 * hd1_8[k]
                  + pa_z[k] * id_14[k];

        t_19[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];

        t_20[k] = f_7 * hd0_12[k]
                  - f_8 * hd1_18[k]
                  + pa_x[k] * id_22[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_z, hd0_7, hd0_13, hd1_11, hd1_20, \
                         id_17, id_23, ks0_7, ks1_7, kp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * hd0_7[k]
                  - f_6 * hd1_11[k]
                  + pa_y[k] * id_17[k];

        t_22[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_20[k]
                  + pa_x[k] * id_23[k];

        t_23[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_z, pb_y, hd0_9, hd0_20, hd1_14, hd1_32, \
                         id_20, id_24, ks0_8, ks1_8, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * hd0_9[k]
                  - f_6 * hd1_14[k]
                  + pa_z[k] * id_20[k];

        t_25[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];

        t_26[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_32[k]
                  + pa_x[k] * id_24[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_x, pb_y, pb_z, ip_10, ks0_9, ks0_10, \
                         ks1_9, ks1_10, kp_11, kp_12, kp_13, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_x[k] * kp_11[k];

        t_28[k] = f_0 * ip_10[k]
                  + f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_y[k] * kp_12[k];

        t_29[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_z[k] * kp_13[k];

        t_30[k] = f_1 * ks0_10[k]
                  - f_2 * ks1_10[k]
                  + pb_x[k] * kp_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pa_z, pb_x, hd0_13, hd0_16, hd1_20, hd1_25, \
                         id_28, id_31, ks0_11, ks1_11, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_20[k]
                  + pa_z[k] * id_28[k];

        t_32[k] = f_5 * hd0_16[k]
                  - f_6 * hd1_25[k]
                  + pa_y[k] * id_31[k];

        t_33[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_15[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pa_z, pb_x, hd0_14, hd0_18, hd1_22, hd1_28, \
                         id_30, id_34, ks0_12, ks1_12, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_7 * hd0_14[k]
                  - f_8 * hd1_22[k]
                  + pa_z[k] * id_30[k];

        t_35[k] = f_9 * hd0_18[k]
                  - f_10 * hd1_28[k]
                  + pa_y[k] * id_34[k];

        t_36[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_16[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, hd0_15, hd0_19, hd1_24, hd1_29, \
                         id_33, id_37, ks0_13, ks1_13, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * hd0_15[k]
                  - f_10 * hd1_24[k]
                  + pa_z[k] * id_33[k];

        t_38[k] = f_7 * hd0_19[k]
                  - f_8 * hd1_29[k]
                  + pa_y[k] * id_37[k];

        t_39[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_17[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_x, hd0_17, hd0_20, hd1_27, hd1_32, \
                         id_36, id_38, ks0_14, ks1_14, kp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * hd0_17[k]
                  - f_6 * hd1_27[k]
                  + pa_z[k] * id_36[k];

        t_41[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_32[k]
                  + pa_y[k] * id_38[k];

        t_42[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_18[k];
    }

#pragma omp simd aligned(t_43, t_44, pb_y, pb_z, ip_17, ks0_14, ks1_14, kp_19, \
                         kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_y[k] * kp_19[k];

        t_44[k] = f_0 * ip_17[k]
                  + f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hd0, const size_t hd1,
                                     const size_t ip, const size_t id, const size_t ks0,
                                     const size_t ks1, const size_t kp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / alpha;
    const auto f_6 = 2.0 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 1.5 / alpha;
    const auto f_10 = 1.5 * beta / (alpha * p);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_27 = buffer.data(hd0 + 27);
    const auto *hd0_28 = buffer.data(hd0 + 28);
    const auto *hd0_29 = buffer.data(hd0 + 29);
    const auto *hd0_32 = buffer.data(hd0 + 32);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_26 = buffer.data(hd1 + 26);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_29 = buffer.data(hd1 + 29);
    const auto *hd1_30 = buffer.data(hd1 + 30);
    const auto *hd1_32 = buffer.data(hd1 + 32);
    const auto *hd1_35 = buffer.data(hd1 + 35);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_17 = buffer.data(ip + 17);

    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_41 = buffer.data(id + 41);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_9 = buffer.data(ks0 + 9);
    const auto *ks0_10 = buffer.data(ks0 + 10);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_9 = buffer.data(ks1 + 9);
    const auto *ks1_10 = buffer.data(ks1 + 10);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, ip_0, ks0_0, ks1_0, kp_0, kp_1, \
                         kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, hd0_0, hd0_6, hd1_0, hd1_7, id_3, \
                         id_7, ks0_1, ks1_1, kp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_y[k] * id_3[k];

        t_4[k] = f_5 * hd0_6[k]
                 - f_6 * hd1_7[k]
                 + pa_x[k] * id_7[k];

        t_5[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_1[k]
                 + pb_z[k] * kp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, hd0_0, hd0_10, hd1_0, hd1_11, id_4, \
                         id_11, ks0_2, ks1_2, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_z[k] * id_4[k];

        t_7[k] = f_1 * ks0_2[k]
                 - f_2 * ks1_2[k]
                 + pb_y[k] * kp_4[k];

        t_8[k] = f_5 * hd0_10[k]
                 - f_6 * hd1_11[k]
                 + pa_x[k] * id_11[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, hd0_3, hd0_12, hd1_3, hd1_13, \
                         id_6, id_13, ks0_3, ks1_3, kp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * hd0_3[k]
                 - f_8 * hd1_3[k]
                 + pa_y[k] * id_6[k];

        t_10[k] = f_9 * hd0_12[k]
                  - f_10 * hd1_13[k]
                  + pa_x[k] * id_13[k];

        t_11[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, hd0_4, hd0_16, hd1_4, hd1_17, \
                         id_9, id_17, ks0_4, ks1_4, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * hd0_4[k]
                  - f_8 * hd1_4[k]
                  + pa_z[k] * id_9[k];

        t_13[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];

        t_14[k] = f_9 * hd0_16[k]
                  - f_10 * hd1_17[k]
                  + pa_x[k] * id_17[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, hd0_5, hd0_17, hd1_6, hd1_18, \
                         id_12, id_19, ks0_5, ks1_5, kp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_9 * hd0_5[k]
                  - f_10 * hd1_6[k]
                  + pa_y[k] * id_12[k];

        t_16[k] = f_7 * hd0_17[k]
                  - f_8 * hd1_18[k]
                  + pa_x[k] * id_19[k];

        t_17[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, hd0_8, hd0_18, hd1_9, hd1_19, \
                         id_15, id_23, ks0_6, ks1_6, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_9 * hd0_8[k]
                  - f_10 * hd1_9[k]
                  + pa_z[k] * id_15[k];

        t_19[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];

        t_20[k] = f_7 * hd0_18[k]
                  - f_8 * hd1_19[k]
                  + pa_x[k] * id_23[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_z, hd0_11, hd0_20, hd1_12, hd1_21, \
                         id_18, id_24, ks0_7, ks1_7, kp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * hd0_11[k]
                  - f_6 * hd1_12[k]
                  + pa_y[k] * id_18[k];

        t_22[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_21[k]
                  + pa_x[k] * id_24[k];

        t_23[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_z, pb_y, hd0_14, hd0_32, hd1_15, hd1_35, \
                         id_21, id_25, ks0_8, ks1_8, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * hd0_14[k]
                  - f_6 * hd1_15[k]
                  + pa_z[k] * id_21[k];

        t_25[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];

        t_26[k] = f_3 * hd0_32[k]
                  - f_4 * hd1_35[k]
                  + pa_x[k] * id_25[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_x, pb_y, pb_z, ip_10, ks0_9, ks0_10, \
                         ks1_9, ks1_10, kp_11, kp_12, kp_13, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_x[k] * kp_11[k];

        t_28[k] = f_0 * ip_10[k]
                  + f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_y[k] * kp_12[k];

        t_29[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_z[k] * kp_13[k];

        t_30[k] = f_1 * ks0_10[k]
                  - f_2 * ks1_10[k]
                  + pb_x[k] * kp_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pa_z, pb_x, hd0_20, hd0_25, hd1_21, hd1_27, \
                         id_29, id_33, ks0_11, ks1_11, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_21[k]
                  + pa_z[k] * id_29[k];

        t_32[k] = f_5 * hd0_25[k]
                  - f_6 * hd1_27[k]
                  + pa_y[k] * id_33[k];

        t_33[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_15[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pa_z, pb_x, hd0_22, hd0_28, hd1_23, hd1_30, \
                         id_32, id_36, ks0_12, ks1_12, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_7 * hd0_22[k]
                  - f_8 * hd1_23[k]
                  + pa_z[k] * id_32[k];

        t_35[k] = f_9 * hd0_28[k]
                  - f_10 * hd1_30[k]
                  + pa_y[k] * id_36[k];

        t_36[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_16[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, hd0_24, hd0_29, hd1_26, hd1_32, \
                         id_35, id_39, ks0_13, ks1_13, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * hd0_24[k]
                  - f_10 * hd1_26[k]
                  + pa_z[k] * id_35[k];

        t_38[k] = f_7 * hd0_29[k]
                  - f_8 * hd1_32[k]
                  + pa_y[k] * id_39[k];

        t_39[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_17[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_x, hd0_27, hd0_32, hd1_29, hd1_35, \
                         id_38, id_41, ks0_14, ks1_14, kp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * hd0_27[k]
                  - f_6 * hd1_29[k]
                  + pa_z[k] * id_38[k];

        t_41[k] = f_3 * hd0_32[k]
                  - f_4 * hd1_35[k]
                  + pa_y[k] * id_41[k];

        t_42[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_18[k];
    }

#pragma omp simd aligned(t_43, t_44, pb_y, pb_z, ip_17, ks0_14, ks1_14, kp_19, \
                         kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_y[k] * kp_19[k];

        t_44[k] = f_0 * ip_17[k]
                  + f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hd0, const size_t hd1,
                                     const size_t ip, const size_t id, const size_t ks0,
                                     const size_t ks1, const size_t kp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / alpha;
    const auto f_6 = 2.0 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 1.5 / alpha;
    const auto f_10 = 1.5 * beta / (alpha * p);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_23 = buffer.data(hd0 + 23);
    const auto *hd0_26 = buffer.data(hd0 + 26);
    const auto *hd0_27 = buffer.data(hd0 + 27);
    const auto *hd0_29 = buffer.data(hd0 + 29);
    const auto *hd0_30 = buffer.data(hd0 + 30);
    const auto *hd0_32 = buffer.data(hd0 + 32);
    const auto *hd0_35 = buffer.data(hd0 + 35);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_28 = buffer.data(hd1 + 28);
    const auto *hd1_29 = buffer.data(hd1 + 29);
    const auto *hd1_32 = buffer.data(hd1 + 32);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_17 = buffer.data(ip + 17);

    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_9 = buffer.data(ks0 + 9);
    const auto *ks0_10 = buffer.data(ks0 + 10);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_9 = buffer.data(ks1 + 9);
    const auto *ks1_10 = buffer.data(ks1 + 10);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, ip_0, ks0_0, ks1_0, kp_0, kp_1, \
                         kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, hd0_0, hd0_7, hd1_0, hd1_6, id_3, \
                         id_6, ks0_1, ks1_1, kp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_y[k] * id_3[k];

        t_4[k] = f_5 * hd0_7[k]
                 - f_6 * hd1_6[k]
                 + pa_x[k] * id_6[k];

        t_5[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_1[k]
                 + pb_z[k] * kp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, hd0_0, hd0_11, hd1_0, hd1_10, id_4, \
                         id_10, ks0_2, ks1_2, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_z[k] * id_4[k];

        t_7[k] = f_1 * ks0_2[k]
                 - f_2 * ks1_2[k]
                 + pb_y[k] * kp_4[k];

        t_8[k] = f_5 * hd0_11[k]
                 - f_6 * hd1_10[k]
                 + pa_x[k] * id_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, hd0_3, hd0_13, hd1_3, hd1_12, \
                         id_5, id_12, ks0_3, ks1_3, kp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * hd0_3[k]
                 - f_8 * hd1_3[k]
                 + pa_y[k] * id_5[k];

        t_10[k] = f_9 * hd0_13[k]
                  - f_10 * hd1_12[k]
                  + pa_x[k] * id_12[k];

        t_11[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, hd0_4, hd0_17, hd1_4, hd1_16, \
                         id_8, id_16, ks0_4, ks1_4, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * hd0_4[k]
                  - f_8 * hd1_4[k]
                  + pa_z[k] * id_8[k];

        t_13[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];

        t_14[k] = f_9 * hd0_17[k]
                  - f_10 * hd1_16[k]
                  + pa_x[k] * id_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, hd0_6, hd0_18, hd1_5, hd1_17, \
                         id_11, id_18, ks0_5, ks1_5, kp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_9 * hd0_6[k]
                  - f_10 * hd1_5[k]
                  + pa_y[k] * id_11[k];

        t_16[k] = f_7 * hd0_18[k]
                  - f_8 * hd1_17[k]
                  + pa_x[k] * id_18[k];

        t_17[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, hd0_9, hd0_19, hd1_8, hd1_18, \
                         id_14, id_22, ks0_6, ks1_6, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_9 * hd0_9[k]
                  - f_10 * hd1_8[k]
                  + pa_z[k] * id_14[k];

        t_19[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];

        t_20[k] = f_7 * hd0_19[k]
                  - f_8 * hd1_18[k]
                  + pa_x[k] * id_22[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_z, hd0_12, hd0_21, hd1_11, hd1_20, \
                         id_17, id_23, ks0_7, ks1_7, kp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * hd0_12[k]
                  - f_6 * hd1_11[k]
                  + pa_y[k] * id_17[k];

        t_22[k] = f_3 * hd0_21[k]
                  - f_4 * hd1_20[k]
                  + pa_x[k] * id_23[k];

        t_23[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_z, pb_y, hd0_15, hd0_35, hd1_14, hd1_32, \
                         id_20, id_24, ks0_8, ks1_8, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * hd0_15[k]
                  - f_6 * hd1_14[k]
                  + pa_z[k] * id_20[k];

        t_25[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];

        t_26[k] = f_3 * hd0_35[k]
                  - f_4 * hd1_32[k]
                  + pa_x[k] * id_24[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_x, pb_y, pb_z, ip_10, ks0_9, ks0_10, \
                         ks1_9, ks1_10, kp_11, kp_12, kp_13, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_x[k] * kp_11[k];

        t_28[k] = f_0 * ip_10[k]
                  + f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_y[k] * kp_12[k];

        t_29[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_z[k] * kp_13[k];

        t_30[k] = f_1 * ks0_10[k]
                  - f_2 * ks1_10[k]
                  + pb_x[k] * kp_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pa_z, pb_x, hd0_21, hd0_27, hd1_20, hd1_25, \
                         id_28, id_31, ks0_11, ks1_11, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_3 * hd0_21[k]
                  - f_4 * hd1_20[k]
                  + pa_z[k] * id_28[k];

        t_32[k] = f_5 * hd0_27[k]
                  - f_6 * hd1_25[k]
                  + pa_y[k] * id_31[k];

        t_33[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_15[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pa_z, pb_x, hd0_23, hd0_30, hd1_22, hd1_28, \
                         id_30, id_34, ks0_12, ks1_12, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_7 * hd0_23[k]
                  - f_8 * hd1_22[k]
                  + pa_z[k] * id_30[k];

        t_35[k] = f_9 * hd0_30[k]
                  - f_10 * hd1_28[k]
                  + pa_y[k] * id_34[k];

        t_36[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_16[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, hd0_26, hd0_32, hd1_24, hd1_29, \
                         id_33, id_37, ks0_13, ks1_13, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * hd0_26[k]
                  - f_10 * hd1_24[k]
                  + pa_z[k] * id_33[k];

        t_38[k] = f_7 * hd0_32[k]
                  - f_8 * hd1_29[k]
                  + pa_y[k] * id_37[k];

        t_39[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_17[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_x, hd0_29, hd0_35, hd1_27, hd1_32, \
                         id_36, id_38, ks0_14, ks1_14, kp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * hd0_29[k]
                  - f_6 * hd1_27[k]
                  + pa_z[k] * id_36[k];

        t_41[k] = f_3 * hd0_35[k]
                  - f_4 * hd1_32[k]
                  + pa_y[k] * id_38[k];

        t_42[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_18[k];
    }

#pragma omp simd aligned(t_43, t_44, pb_y, pb_z, ip_17, ks0_14, ks1_14, kp_19, \
                         kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_y[k] * kp_19[k];

        t_44[k] = f_0 * ip_17[k]
                  + f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hd0, const size_t hd1,
                                     const size_t ip, const size_t id, const size_t ks0,
                                     const size_t ks1, const size_t kp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / alpha;
    const auto f_6 = 2.0 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 1.5 / alpha;
    const auto f_10 = 1.5 * beta / (alpha * p);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_27 = buffer.data(hd0 + 27);
    const auto *hd0_28 = buffer.data(hd0 + 28);
    const auto *hd0_29 = buffer.data(hd0 + 29);
    const auto *hd0_32 = buffer.data(hd0 + 32);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_28 = buffer.data(hd1 + 28);
    const auto *hd1_29 = buffer.data(hd1 + 29);
    const auto *hd1_32 = buffer.data(hd1 + 32);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_17 = buffer.data(ip + 17);

    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_9 = buffer.data(ks0 + 9);
    const auto *ks0_10 = buffer.data(ks0 + 10);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_9 = buffer.data(ks1 + 9);
    const auto *ks1_10 = buffer.data(ks1 + 10);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, ip_0, ks0_0, ks1_0, kp_0, kp_1, \
                         kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, hd0_0, hd0_6, hd1_0, hd1_6, id_3, \
                         id_6, ks0_1, ks1_1, kp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_y[k] * id_3[k];

        t_4[k] = f_5 * hd0_6[k]
                 - f_6 * hd1_6[k]
                 + pa_x[k] * id_6[k];

        t_5[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_1[k]
                 + pb_z[k] * kp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, hd0_0, hd0_10, hd1_0, hd1_10, id_4, \
                         id_10, ks0_2, ks1_2, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_z[k] * id_4[k];

        t_7[k] = f_1 * ks0_2[k]
                 - f_2 * ks1_2[k]
                 + pb_y[k] * kp_4[k];

        t_8[k] = f_5 * hd0_10[k]
                 - f_6 * hd1_10[k]
                 + pa_x[k] * id_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, hd0_3, hd0_12, hd1_3, hd1_12, \
                         id_5, id_12, ks0_3, ks1_3, kp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * hd0_3[k]
                 - f_8 * hd1_3[k]
                 + pa_y[k] * id_5[k];

        t_10[k] = f_9 * hd0_12[k]
                  - f_10 * hd1_12[k]
                  + pa_x[k] * id_12[k];

        t_11[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, hd0_4, hd0_16, hd1_4, hd1_16, \
                         id_8, id_16, ks0_4, ks1_4, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * hd0_4[k]
                  - f_8 * hd1_4[k]
                  + pa_z[k] * id_8[k];

        t_13[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];

        t_14[k] = f_9 * hd0_16[k]
                  - f_10 * hd1_16[k]
                  + pa_x[k] * id_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, hd0_5, hd0_17, hd1_5, hd1_17, \
                         id_11, id_18, ks0_5, ks1_5, kp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_9 * hd0_5[k]
                  - f_10 * hd1_5[k]
                  + pa_y[k] * id_11[k];

        t_16[k] = f_7 * hd0_17[k]
                  - f_8 * hd1_17[k]
                  + pa_x[k] * id_18[k];

        t_17[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, hd0_8, hd0_18, hd1_8, hd1_18, \
                         id_14, id_22, ks0_6, ks1_6, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_9 * hd0_8[k]
                  - f_10 * hd1_8[k]
                  + pa_z[k] * id_14[k];

        t_19[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];

        t_20[k] = f_7 * hd0_18[k]
                  - f_8 * hd1_18[k]
                  + pa_x[k] * id_22[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_z, hd0_11, hd0_20, hd1_11, hd1_20, \
                         id_17, id_23, ks0_7, ks1_7, kp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * hd0_11[k]
                  - f_6 * hd1_11[k]
                  + pa_y[k] * id_17[k];

        t_22[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pa_x[k] * id_23[k];

        t_23[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_z, pb_y, hd0_14, hd0_32, hd1_14, hd1_32, \
                         id_20, id_24, ks0_8, ks1_8, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * hd0_14[k]
                  - f_6 * hd1_14[k]
                  + pa_z[k] * id_20[k];

        t_25[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];

        t_26[k] = f_3 * hd0_32[k]
                  - f_4 * hd1_32[k]
                  + pa_x[k] * id_24[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_x, pb_y, pb_z, ip_10, ks0_9, ks0_10, \
                         ks1_9, ks1_10, kp_11, kp_12, kp_13, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_x[k] * kp_11[k];

        t_28[k] = f_0 * ip_10[k]
                  + f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_y[k] * kp_12[k];

        t_29[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_z[k] * kp_13[k];

        t_30[k] = f_1 * ks0_10[k]
                  - f_2 * ks1_10[k]
                  + pb_x[k] * kp_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pa_z, pb_x, hd0_20, hd0_25, hd1_20, hd1_25, \
                         id_28, id_31, ks0_11, ks1_11, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pa_z[k] * id_28[k];

        t_32[k] = f_5 * hd0_25[k]
                  - f_6 * hd1_25[k]
                  + pa_y[k] * id_31[k];

        t_33[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_15[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pa_z, pb_x, hd0_22, hd0_28, hd1_22, hd1_28, \
                         id_30, id_34, ks0_12, ks1_12, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_7 * hd0_22[k]
                  - f_8 * hd1_22[k]
                  + pa_z[k] * id_30[k];

        t_35[k] = f_9 * hd0_28[k]
                  - f_10 * hd1_28[k]
                  + pa_y[k] * id_34[k];

        t_36[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_16[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, hd0_24, hd0_29, hd1_24, hd1_29, \
                         id_33, id_37, ks0_13, ks1_13, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * hd0_24[k]
                  - f_10 * hd1_24[k]
                  + pa_z[k] * id_33[k];

        t_38[k] = f_7 * hd0_29[k]
                  - f_8 * hd1_29[k]
                  + pa_y[k] * id_37[k];

        t_39[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_17[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_x, hd0_27, hd0_32, hd1_27, hd1_32, \
                         id_36, id_38, ks0_14, ks1_14, kp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * hd0_27[k]
                  - f_6 * hd1_27[k]
                  + pa_z[k] * id_36[k];

        t_41[k] = f_3 * hd0_32[k]
                  - f_4 * hd1_32[k]
                  + pa_y[k] * id_38[k];

        t_42[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_18[k];
    }

#pragma omp simd aligned(t_43, t_44, pb_y, pb_z, ip_17, ks0_14, ks1_14, kp_19, \
                         kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_y[k] * kp_19[k];

        t_44[k] = f_0 * ip_17[k]
                  + f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hd0, const size_t hd1,
                                     const size_t ip, const size_t id, const size_t ks0,
                                     const size_t ks1, const size_t kp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 3.0 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 2.5 / p;
    const auto f_8 = 2.0 / alpha;
    const auto f_9 = 2.0 * beta / (alpha * p);
    const auto f_10 = 1.0 / alpha;
    const auto f_11 = beta / (alpha * p);
    const auto f_12 = 2.0 / p;
    const auto f_13 = 1.5 / alpha;
    const auto f_14 = 1.5 * beta / (alpha * p);
    const auto f_15 = 1.5 / p;
    const auto f_16 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_23 = buffer.data(hd0 + 23);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_23 = buffer.data(hd1 + 23);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_27 = buffer.data(ip + 27);

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

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);
    const auto *ks0_15 = buffer.data(ks0 + 15);
    const auto *ks0_17 = buffer.data(ks0 + 17);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_9 = buffer.data(ks1 + 9);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);
    const auto *ks1_19 = buffer.data(ks1 + 19);
    const auto *ks1_22 = buffer.data(ks1 + 22);
    const auto *ks1_24 = buffer.data(ks1 + 24);
    const auto *ks1_25 = buffer.data(ks1 + 25);
    const auto *ks1_26 = buffer.data(ks1 + 26);
    const auto *ks1_27 = buffer.data(ks1 + 27);
    const auto *ks1_29 = buffer.data(ks1 + 29);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_21 = buffer.data(kp + 21);
    const auto *kp_22 = buffer.data(kp + 22);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_24 = buffer.data(kp + 24);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_30 = buffer.data(kp + 30);
    const auto *kp_32 = buffer.data(kp + 32);
    const auto *kp_33 = buffer.data(kp + 33);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_36 = buffer.data(kp + 36);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_42 = buffer.data(kp + 42);
    const auto *kp_43 = buffer.data(kp + 43);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, ip_0, id_0, ks0_0, ks1_0, \
                         kp_0, kp_1, kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];

        t_3[k] = pa_y[k] * id_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, pb_x, ip_1, ip_3, ip_4, id_0, id_1, \
                         kp_3, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * ip_3[k]
                 + pb_x[k] * kp_3[k];

        t_5[k] = f_4 * ip_1[k]
                 + pa_y[k] * id_1[k];

        t_6[k] = pa_z[k] * id_0[k];

        t_7[k] = f_3 * ip_4[k]
                 + pb_x[k] * kp_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_y, pa_z, pb_x, hd0_0, hd1_0, ip_2, ip_5, id_2, \
                         id_3, kp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * ip_2[k]
                 + pa_z[k] * id_2[k];

        t_9[k] = f_5 * hd0_0[k]
                 - f_6 * hd1_0[k]
                 + pa_y[k] * id_3[k];

        t_10[k] = f_7 * ip_5[k]
                  + pb_x[k] * kp_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_z, pb_z, hd0_0, hd0_4, hd1_0, hd1_4, id_4, \
                         id_6, ks0_1, ks1_3, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_8 * hd0_4[k]
                  - f_9 * hd1_4[k]
                  + pa_x[k] * id_6[k];

        t_12[k] = f_1 * ks0_1[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_6[k];

        t_13[k] = f_5 * hd0_0[k]
                  - f_6 * hd1_0[k]
                  + pa_z[k] * id_4[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pb_x, pb_y, hd0_6, hd1_6, ip_6, id_8, ks0_2, \
                         ks1_4, kp_7, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_7 * ip_6[k]
                  + pb_x[k] * kp_8[k];

        t_15[k] = f_1 * ks0_2[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_7[k];

        t_16[k] = f_8 * hd0_6[k]
                  - f_9 * hd1_6[k]
                  + pa_x[k] * id_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, pb_x, hd0_1, hd0_8, hd1_1, hd1_8, ip_7, \
                         id_5, id_10, kp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_10 * hd0_1[k]
                  - f_11 * hd1_1[k]
                  + pa_y[k] * id_5[k];

        t_18[k] = f_12 * ip_7[k]
                  + pb_x[k] * kp_9[k];

        t_19[k] = f_13 * hd0_8[k]
                  - f_14 * hd1_8[k]
                  + pa_x[k] * id_10[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_y, pa_z, pb_x, pb_z, hd0_2, hd1_2, ip_8, \
                         id_7, ks0_3, ks1_5, kp_10, kp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_10[k];

        t_21[k] = pa_y[k] * id_7[k];

        t_22[k] = f_10 * hd0_2[k]
                  - f_11 * hd1_2[k]
                  + pa_z[k] * id_7[k];

        t_23[k] = f_12 * ip_8[k]
                  + pb_x[k] * kp_12[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_y, pb_y, hd0_3, hd0_11, hd1_3, hd1_11, \
                         id_9, id_13, ks0_4, ks1_8, kp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_11[k];

        t_25[k] = f_13 * hd0_11[k]
                  - f_14 * hd1_11[k]
                  + pa_x[k] * id_13[k];

        t_26[k] = f_13 * hd0_3[k]
                  - f_14 * hd1_3[k]
                  + pa_y[k] * id_9[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pb_x, pb_z, hd0_12, hd1_12, ip_9, id_15, \
                         ks0_5, ks1_9, kp_13, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_15 * ip_9[k]
                  + pb_x[k] * kp_13[k];

        t_28[k] = f_10 * hd0_12[k]
                  - f_11 * hd1_12[k]
                  + pa_x[k] * id_15[k];

        t_29[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_9[k]
                  + pb_z[k] * kp_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, hd0_5, hd0_13, hd0_14, hd1_5, \
                         hd1_13, hd1_14, id_11, id_12, id_17, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * hd0_5[k]
                  - f_6 * hd1_5[k]
                  + pa_y[k] * id_11[k];

        t_31[k] = f_10 * hd0_13[k]
                  - f_11 * hd1_13[k]
                  + pa_x[k] * id_17[k];

        t_32[k] = f_10 * hd0_14[k]
                  - f_11 * hd1_14[k]
                  + pa_x[k] * id_18[k];

        t_33[k] = pa_y[k] * id_12[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_x, pb_y, hd0_5, hd1_5, ip_10, id_12, \
                         ks0_6, ks1_13, kp_15, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_13 * hd0_5[k]
                  - f_14 * hd1_5[k]
                  + pa_z[k] * id_12[k];

        t_35[k] = f_15 * ip_10[k]
                  + pb_x[k] * kp_16[k];

        t_36[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_13[k]
                  + pb_y[k] * kp_15[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_x, pa_y, pb_x, hd0_7, hd0_15, hd1_7, hd1_15, \
                         ip_11, id_14, id_21, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_10 * hd0_15[k]
                  - f_11 * hd1_15[k]
                  + pa_x[k] * id_21[k];

        t_38[k] = f_8 * hd0_7[k]
                  - f_9 * hd1_7[k]
                  + pa_y[k] * id_14[k];

        t_39[k] = f_4 * ip_11[k]
                  + pb_x[k] * kp_17[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_x, pa_y, pb_z, hd0_9, hd0_16, hd1_9, hd1_16, \
                         id_16, id_22, ks0_7, ks1_14, kp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * hd0_16[k]
                  - f_6 * hd1_16[k]
                  + pa_x[k] * id_22[k];

        t_41[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_14[k]
                  + pb_z[k] * kp_18[k];

        t_42[k] = f_10 * hd0_9[k]
                  - f_11 * hd1_9[k]
                  + pa_y[k] * id_16[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_x, pa_y, hd0_10, hd0_18, hd0_19, hd1_10, hd1_18, \
                         hd1_19, id_19, id_23, id_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_5 * hd0_18[k]
                  - f_6 * hd1_18[k]
                  + pa_x[k] * id_23[k];

        t_44[k] = f_5 * hd0_19[k]
                  - f_6 * hd1_19[k]
                  + pa_x[k] * id_24[k];

        t_45[k] = f_5 * hd0_10[k]
                  - f_6 * hd1_10[k]
                  + pa_y[k] * id_19[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pa_y, pa_z, hd0_10, hd0_20, hd0_21, \
                         hd1_10, hd1_20, hd1_21, id_20, id_25, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_5 * hd0_20[k]
                  - f_6 * hd1_20[k]
                  + pa_x[k] * id_25[k];

        t_47[k] = f_5 * hd0_21[k]
                  - f_6 * hd1_21[k]
                  + pa_x[k] * id_26[k];

        t_48[k] = pa_y[k] * id_20[k];

        t_49[k] = f_8 * hd0_10[k]
                  - f_9 * hd1_10[k]
                  + pa_z[k] * id_20[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_x, pb_x, pb_y, hd0_23, hd1_23, ip_12, id_27, \
                         ks0_8, ks1_19, kp_19, kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_4 * ip_12[k]
                  + pb_x[k] * kp_20[k];

        t_51[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_19[k]
                  + pb_y[k] * kp_19[k];

        t_52[k] = f_5 * hd0_23[k]
                  - f_6 * hd1_23[k]
                  + pa_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, t_58, pa_x, pb_x, ip_13, ip_14, id_28, \
                         id_29, id_32, id_33, id_34, kp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_4 * ip_13[k]
                  + pa_x[k] * id_28[k];

        t_54[k] = f_16 * ip_14[k]
                  + pb_x[k] * kp_21[k];

        t_55[k] = pa_x[k] * id_29[k];

        t_56[k] = pa_x[k] * id_32[k];

        t_57[k] = pa_x[k] * id_33[k];

        t_58[k] = pa_x[k] * id_34[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, t_63, t_64, pa_x, pb_x, ip_25, ip_27, id_35, \
                         id_36, id_37, id_39, id_41, kp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_x[k] * id_35[k];

        t_60[k] = pa_x[k] * id_36[k];

        t_61[k] = pa_x[k] * id_37[k];

        t_62[k] = f_4 * ip_25[k]
                  + pa_x[k] * id_39[k];

        t_63[k] = f_16 * ip_27[k]
                  + pb_x[k] * kp_22[k];

        t_64[k] = pa_x[k] * id_41[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_z, pb_x, pb_y, pb_z, ip_14, id_29, ks0_11, \
                         ks1_22, kp_23, kp_24, kp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_22[k]
                  + pb_x[k] * kp_23[k];

        t_66[k] = f_0 * ip_14[k]
                  + f_1 * ks0_11[k]
                  - f_2 * ks1_22[k]
                  + pb_y[k] * kp_24[k];

        t_67[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_22[k]
                  + pb_z[k] * kp_25[k];

        t_68[k] = pa_z[k] * id_29[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_z, pb_x, pb_y, ip_15, ip_16, id_30, ks0_12, \
                         ks1_24, kp_26, kp_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * ip_16[k]
                  + pb_y[k] * kp_26[k];

        t_70[k] = f_4 * ip_15[k]
                  + pa_z[k] * id_30[k];

        t_71[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_24[k]
                  + pb_x[k] * kp_27[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_y, pa_z, pb_y, hd0_16, hd0_19, hd1_16, hd1_19, \
                         ip_18, id_31, id_33, kp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_5 * hd0_16[k]
                  - f_6 * hd1_16[k]
                  + pa_z[k] * id_31[k];

        t_73[k] = f_7 * ip_18[k]
                  + pb_y[k] * kp_29[k];

        t_74[k] = f_8 * hd0_19[k]
                  - f_9 * hd1_19[k]
                  + pa_y[k] * id_33[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_z, pb_x, pb_y, hd0_17, hd1_17, ip_20, id_32, \
                         ks0_13, ks1_25, kp_30, kp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_25[k]
                  + pb_x[k] * kp_30[k];

        t_76[k] = f_10 * hd0_17[k]
                  - f_11 * hd1_17[k]
                  + pa_z[k] * id_32[k];

        t_77[k] = f_12 * ip_20[k]
                  + pb_y[k] * kp_32[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pa_y, pa_z, pb_x, hd0_18, hd0_21, hd1_18, hd1_21, \
                         id_34, id_35, ks0_14, ks1_26, kp_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_13 * hd0_21[k]
                  - f_14 * hd1_21[k]
                  + pa_y[k] * id_35[k];

        t_79[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_26[k]
                  + pb_x[k] * kp_33[k];

        t_80[k] = f_13 * hd0_18[k]
                  - f_14 * hd1_18[k]
                  + pa_z[k] * id_34[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_y, pb_x, pb_y, hd0_22, hd1_22, ip_22, id_37, \
                         ks0_15, ks1_27, kp_35, kp_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_15 * ip_22[k]
                  + pb_y[k] * kp_35[k];

        t_82[k] = f_10 * hd0_22[k]
                  - f_11 * hd1_22[k]
                  + pa_y[k] * id_37[k];

        t_83[k] = f_1 * ks0_15[k]
                  - f_2 * ks1_27[k]
                  + pb_x[k] * kp_36[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_y, pa_z, pb_y, hd0_20, hd0_23, hd1_20, hd1_23, \
                         ip_24, id_36, id_38, kp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_8 * hd0_20[k]
                  - f_9 * hd1_20[k]
                  + pa_z[k] * id_36[k];

        t_85[k] = f_4 * ip_24[k]
                  + pb_y[k] * kp_38[k];

        t_86[k] = f_5 * hd0_23[k]
                  - f_6 * hd1_23[k]
                  + pa_y[k] * id_38[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pb_x, pb_y, ip_26, ip_27, id_40, id_41, \
                         ks0_17, ks1_29, kp_40, kp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * ip_26[k]
                  + pa_y[k] * id_40[k];

        t_88[k] = f_16 * ip_27[k]
                  + pb_y[k] * kp_40[k];

        t_89[k] = pa_y[k] * id_41[k];

        t_90[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_29[k]
                  + pb_x[k] * kp_41[k];
    }

#pragma omp simd aligned(t_91, t_92, pb_y, pb_z, ip_27, ks0_17, ks1_29, kp_42, \
                         kp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_29[k]
                  + pb_y[k] * kp_42[k];

        t_92[k] = f_0 * ip_27[k]
                  + f_1 * ks0_17[k]
                  - f_2 * ks1_29[k]
                  + pb_z[k] * kp_43[k];
    }
}

auto
compute_prim_kd_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hd0, const size_t hd1,
                                      const size_t ip, const size_t id, const size_t ks0,
                                      const size_t ks1, const size_t kp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 0.5 / p;
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 3.0 / p;
    const auto f_16 = 2.5 / p;

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_23 = buffer.data(hd0 + 23);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_23 = buffer.data(hd1 + 23);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_21 = buffer.data(ip + 21);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_27 = buffer.data(ip + 27);
    const auto *ip_28 = buffer.data(ip + 28);
    const auto *ip_29 = buffer.data(ip + 29);
    const auto *ip_30 = buffer.data(ip + 30);
    const auto *ip_31 = buffer.data(ip + 31);
    const auto *ip_32 = buffer.data(ip + 32);

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

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_9 = buffer.data(ks0 + 9);
    const auto *ks0_10 = buffer.data(ks0 + 10);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);
    const auto *ks1_15 = buffer.data(ks1 + 15);
    const auto *ks1_17 = buffer.data(ks1 + 17);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_21 = buffer.data(kp + 21);
    const auto *kp_24 = buffer.data(kp + 24);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_28 = buffer.data(kp + 28);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_33 = buffer.data(kp + 33);
    const auto *kp_34 = buffer.data(kp + 34);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_36 = buffer.data(kp + 36);
    const auto *kp_37 = buffer.data(kp + 37);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_39 = buffer.data(kp + 39);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_42 = buffer.data(kp + 42);
    const auto *kp_43 = buffer.data(kp + 43);
    const auto *kp_44 = buffer.data(kp + 44);
    const auto *kp_45 = buffer.data(kp + 45);
    const auto *kp_46 = buffer.data(kp + 46);
    const auto *kp_47 = buffer.data(kp + 47);
    const auto *kp_48 = buffer.data(kp + 48);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, ip_0, id_0, ks0_0, ks1_0, \
                         kp_0, kp_1, kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];

        t_3[k] = pa_y[k] * id_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, hd0_0, hd1_0, ip_1, ip_2, \
                         id_0, id_1, id_2, id_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * ip_1[k]
                 + pa_y[k] * id_1[k];

        t_5[k] = pa_y[k] * id_2[k];

        t_6[k] = pa_z[k] * id_0[k];

        t_7[k] = pa_z[k] * id_1[k];

        t_8[k] = f_3 * ip_2[k]
                 + pa_z[k] * id_2[k];

        t_9[k] = f_4 * hd0_0[k]
                 - f_5 * hd1_0[k]
                 + pa_y[k] * id_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_z, pb_z, hd0_4, hd1_4, id_4, id_8, ks0_1, \
                         ks1_1, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_6 * hd0_4[k]
                  - f_7 * hd1_4[k]
                  + pa_x[k] * id_8[k];

        t_11[k] = f_1 * ks0_1[k]
                  - f_2 * ks1_1[k]
                  + pb_z[k] * kp_6[k];

        t_12[k] = pa_z[k] * id_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pa_z, pb_y, hd0_0, hd1_0, ip_3, id_5, \
                         id_6, ks0_2, ks1_2, kp_7, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_8 * ip_3[k]
                  + pb_y[k] * kp_7[k];

        t_14[k] = pa_y[k] * id_6[k];

        t_15[k] = f_4 * hd0_0[k]
                  - f_5 * hd1_0[k]
                  + pa_z[k] * id_5[k];

        t_16[k] = f_1 * ks0_2[k]
                  - f_2 * ks1_2[k]
                  + pb_y[k] * kp_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, hd0_1, hd0_6, hd0_8, hd1_1, hd1_6, \
                         hd1_8, id_7, id_12, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * hd0_6[k]
                  - f_7 * hd1_6[k]
                  + pa_x[k] * id_12[k];

        t_18[k] = f_9 * hd0_1[k]
                  - f_10 * hd1_1[k]
                  + pa_y[k] * id_7[k];

        t_19[k] = f_11 * hd0_8[k]
                  - f_12 * hd1_8[k]
                  + pa_x[k] * id_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_z, pb_y, pb_z, ip_5, id_7, id_8, ks0_3, \
                         ks1_3, kp_11, kp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_11[k];

        t_21[k] = pa_z[k] * id_7[k];

        t_22[k] = pa_z[k] * id_8[k];

        t_23[k] = f_3 * ip_5[k]
                  + pb_y[k] * kp_12[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_y, pa_z, pb_y, ip_4, ip_6, ip_7, \
                         id_9, id_10, id_11, id_12, kp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_3 * ip_4[k]
                  + pa_z[k] * id_9[k];

        t_25[k] = pa_y[k] * id_10[k];

        t_26[k] = f_3 * ip_6[k]
                  + pa_y[k] * id_11[k];

        t_27[k] = f_8 * ip_7[k]
                  + pb_y[k] * kp_13[k];

        t_28[k] = pa_y[k] * id_12[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_x, pa_z, pb_y, hd0_2, hd0_11, hd1_2, hd1_11, \
                         id_10, id_19, ks0_4, ks1_4, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_9 * hd0_2[k]
                  - f_10 * hd1_2[k]
                  + pa_z[k] * id_10[k];

        t_30[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_14[k];

        t_31[k] = f_11 * hd0_11[k]
                  - f_12 * hd1_11[k]
                  + pa_x[k] * id_19[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_x, pa_y, pb_z, hd0_3, hd0_12, hd1_3, hd1_12, \
                         id_13, id_21, ks0_5, ks1_5, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_11 * hd0_3[k]
                  - f_12 * hd1_3[k]
                  + pa_y[k] * id_13[k];

        t_33[k] = f_9 * hd0_12[k]
                  - f_10 * hd1_12[k]
                  + pa_x[k] * id_21[k];

        t_34[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_z, pb_y, ip_8, ip_9, id_13, id_14, id_15, \
                         kp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_z[k] * id_13[k];

        t_36[k] = pa_z[k] * id_14[k];

        t_37[k] = f_13 * ip_9[k]
                  + pb_y[k] * kp_18[k];

        t_38[k] = f_3 * ip_8[k]
                  + pa_z[k] * id_15[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_x, pa_y, pb_y, hd0_5, hd0_13, hd1_5, hd1_13, \
                         ip_10, id_16, id_24, kp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_4 * hd0_5[k]
                  - f_5 * hd1_5[k]
                  + pa_y[k] * id_16[k];

        t_40[k] = f_9 * hd0_13[k]
                  - f_10 * hd1_13[k]
                  + pa_x[k] * id_24[k];

        t_41[k] = f_3 * ip_10[k]
                  + pb_y[k] * kp_19[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pa_y, pb_y, hd0_14, hd1_14, ip_11, \
                         ip_12, id_17, id_18, id_25, kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_9 * hd0_14[k]
                  - f_10 * hd1_14[k]
                  + pa_x[k] * id_25[k];

        t_43[k] = pa_y[k] * id_17[k];

        t_44[k] = f_3 * ip_11[k]
                  + pa_y[k] * id_18[k];

        t_45[k] = f_8 * ip_12[k]
                  + pb_y[k] * kp_20[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_y, pa_z, pb_y, hd0_5, hd1_5, id_17, id_19, \
                         ks0_6, ks1_6, kp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_y[k] * id_19[k];

        t_47[k] = f_11 * hd0_5[k]
                  - f_12 * hd1_5[k]
                  + pa_z[k] * id_17[k];

        t_48[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_21[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_x, pa_y, hd0_7, hd0_15, hd0_16, hd1_7, hd1_15, \
                         hd1_16, id_20, id_29, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_9 * hd0_15[k]
                  - f_10 * hd1_15[k]
                  + pa_x[k] * id_29[k];

        t_50[k] = f_6 * hd0_7[k]
                  - f_7 * hd1_7[k]
                  + pa_y[k] * id_20[k];

        t_51[k] = f_4 * hd0_16[k]
                  - f_5 * hd1_16[k]
                  + pa_x[k] * id_31[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_z, pb_y, pb_z, ip_14, id_20, id_21, ks0_7, \
                         ks1_7, kp_24, kp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_24[k];

        t_53[k] = pa_z[k] * id_20[k];

        t_54[k] = pa_z[k] * id_21[k];

        t_55[k] = f_14 * ip_14[k]
                  + pb_y[k] * kp_25[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pa_x, pa_y, pa_z, hd0_9, hd0_18, hd1_9, hd1_18, \
                         ip_13, id_22, id_23, id_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_3 * ip_13[k]
                  + pa_z[k] * id_22[k];

        t_57[k] = f_9 * hd0_9[k]
                  - f_10 * hd1_9[k]
                  + pa_y[k] * id_23[k];

        t_58[k] = f_4 * hd0_18[k]
                  - f_5 * hd1_18[k]
                  + pa_x[k] * id_32[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pa_x, pa_y, pb_y, hd0_10, hd0_19, hd1_10, hd1_19, \
                         ip_15, id_26, id_33, kp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_13 * ip_15[k]
                  + pb_y[k] * kp_26[k];

        t_60[k] = f_4 * hd0_19[k]
                  - f_5 * hd1_19[k]
                  + pa_x[k] * id_33[k];

        t_61[k] = f_4 * hd0_10[k]
                  - f_5 * hd1_10[k]
                  + pa_y[k] * id_26[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pa_y, pb_y, hd0_20, hd0_21, hd1_20, \
                         hd1_21, ip_16, id_27, id_34, id_35, kp_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_4 * hd0_20[k]
                  - f_5 * hd1_20[k]
                  + pa_x[k] * id_34[k];

        t_63[k] = f_3 * ip_16[k]
                  + pb_y[k] * kp_27[k];

        t_64[k] = f_4 * hd0_21[k]
                  - f_5 * hd1_21[k]
                  + pa_x[k] * id_35[k];

        t_65[k] = pa_y[k] * id_27[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_y, pa_z, pb_y, hd0_10, hd1_10, ip_17, \
                         ip_18, id_27, id_28, id_29, kp_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_3 * ip_17[k]
                  + pa_y[k] * id_28[k];

        t_67[k] = f_8 * ip_18[k]
                  + pb_y[k] * kp_28[k];

        t_68[k] = pa_y[k] * id_29[k];

        t_69[k] = f_6 * hd0_10[k]
                  - f_7 * hd1_10[k]
                  + pa_z[k] * id_27[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_x, pb_y, hd0_23, hd1_23, ip_19, id_37, \
                         id_38, id_39, ks0_8, ks1_8, kp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_29[k];

        t_71[k] = f_4 * hd0_23[k]
                  - f_5 * hd1_23[k]
                  + pa_x[k] * id_37[k];

        t_72[k] = f_3 * ip_19[k]
                  + pa_x[k] * id_38[k];

        t_73[k] = pa_x[k] * id_39[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, pa_x, pa_z, ip_23, id_30, id_40, \
                         id_42, id_43, id_44, id_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pa_x[k] * id_40[k];

        t_75[k] = pa_z[k] * id_30[k];

        t_76[k] = pa_x[k] * id_42[k];

        t_77[k] = pa_x[k] * id_43[k];

        t_78[k] = f_3 * ip_23[k]
                  + pa_x[k] * id_44[k];

        t_79[k] = pa_x[k] * id_45[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, t_85, pa_x, ip_25, id_46, id_47, id_48, \
                         id_49, id_50, id_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pa_x[k] * id_46[k];

        t_81[k] = pa_x[k] * id_47[k];

        t_82[k] = f_3 * ip_25[k]
                  + pa_x[k] * id_48[k];

        t_83[k] = pa_x[k] * id_49[k];

        t_84[k] = pa_x[k] * id_50[k];

        t_85[k] = pa_x[k] * id_51[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, t_91, pa_x, pa_y, ip_27, id_36, id_52, \
                         id_53, id_54, id_55, id_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_3 * ip_27[k]
                  + pa_x[k] * id_52[k];

        t_87[k] = pa_x[k] * id_53[k];

        t_88[k] = pa_x[k] * id_54[k];

        t_89[k] = pa_x[k] * id_55[k];

        t_90[k] = pa_y[k] * id_36[k];

        t_91[k] = pa_x[k] * id_56[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, pa_x, pb_x, ip_30, id_57, id_59, id_60, \
                         id_61, ks0_9, ks1_11, kp_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_x[k] * id_57[k];

        t_93[k] = f_3 * ip_30[k]
                  + pa_x[k] * id_59[k];

        t_94[k] = pa_x[k] * id_60[k];

        t_95[k] = pa_x[k] * id_61[k];

        t_96[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_33[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_z, pb_y, pb_z, ip_20, id_38, id_39, \
                         ks0_9, ks1_11, kp_34, kp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_0 * ip_20[k]
                  + f_1 * ks0_9[k]
                  - f_2 * ks1_11[k]
                  + pb_y[k] * kp_34[k];

        t_98[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_11[k]
                  + pb_z[k] * kp_35[k];

        t_99[k] = pa_z[k] * id_38[k];

        t_100[k] = pa_z[k] * id_39[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pa_z, pb_x, pb_y, ip_21, ip_22, id_40, ks0_10, \
                         ks1_12, kp_36, kp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_15 * ip_22[k]
                   + pb_y[k] * kp_36[k];

        t_102[k] = f_3 * ip_21[k]
                   + pa_z[k] * id_40[k];

        t_103[k] = f_1 * ks0_10[k]
                   - f_2 * ks1_12[k]
                   + pb_x[k] * kp_37[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_y, pa_z, pb_y, hd0_16, hd0_19, hd1_16, \
                         hd1_19, ip_24, id_41, id_47, kp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_4 * hd0_16[k]
                   - f_5 * hd1_16[k]
                   + pa_z[k] * id_41[k];

        t_105[k] = f_16 * ip_24[k]
                   + pb_y[k] * kp_38[k];

        t_106[k] = f_6 * hd0_19[k]
                   - f_7 * hd1_19[k]
                   + pa_y[k] * id_47[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pa_z, pb_x, pb_y, hd0_17, hd1_17, ip_26, id_45, \
                         ks0_11, ks1_13, kp_39, kp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_1 * ks0_11[k]
                   - f_2 * ks1_13[k]
                   + pb_x[k] * kp_39[k];

        t_108[k] = f_9 * hd0_17[k]
                   - f_10 * hd1_17[k]
                   + pa_z[k] * id_45[k];

        t_109[k] = f_14 * ip_26[k]
                   + pb_y[k] * kp_40[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pa_y, pa_z, pb_x, hd0_18, hd0_21, hd1_18, \
                         hd1_21, id_49, id_51, ks0_12, ks1_14, kp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_11 * hd0_21[k]
                   - f_12 * hd1_21[k]
                   + pa_y[k] * id_51[k];

        t_111[k] = f_1 * ks0_12[k]
                   - f_2 * ks1_14[k]
                   + pb_x[k] * kp_41[k];

        t_112[k] = f_11 * hd0_18[k]
                   - f_12 * hd1_18[k]
                   + pa_z[k] * id_49[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pa_y, pb_x, pb_y, hd0_22, hd1_22, ip_28, id_55, \
                         ks0_13, ks1_15, kp_42, kp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_13 * ip_28[k]
                   + pb_y[k] * kp_42[k];

        t_114[k] = f_9 * hd0_22[k]
                   - f_10 * hd1_22[k]
                   + pa_y[k] * id_55[k];

        t_115[k] = f_1 * ks0_13[k]
                   - f_2 * ks1_15[k]
                   + pb_x[k] * kp_43[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_y, pa_z, pb_y, hd0_20, hd0_23, hd1_20, \
                         hd1_23, ip_29, id_53, id_58, id_59, kp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_6 * hd0_20[k]
                   - f_7 * hd1_20[k]
                   + pa_z[k] * id_53[k];

        t_117[k] = f_3 * ip_29[k]
                   + pb_y[k] * kp_44[k];

        t_118[k] = f_4 * hd0_23[k]
                   - f_5 * hd1_23[k]
                   + pa_y[k] * id_58[k];

        t_119[k] = pa_y[k] * id_59[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_y, pb_x, pb_y, ip_31, ip_32, id_60, \
                         id_61, ks0_14, ks1_17, kp_45, kp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_3 * ip_31[k]
                   + pa_y[k] * id_60[k];

        t_121[k] = f_8 * ip_32[k]
                   + pb_y[k] * kp_45[k];

        t_122[k] = pa_y[k] * id_61[k];

        t_123[k] = f_1 * ks0_14[k]
                   - f_2 * ks1_17[k]
                   + pb_x[k] * kp_46[k];
    }

#pragma omp simd aligned(t_124, t_125, pb_y, pb_z, ip_32, ks0_14, ks1_17, kp_47, \
                         kp_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_1 * ks0_14[k]
                   - f_2 * ks1_17[k]
                   + pb_y[k] * kp_47[k];

        t_125[k] = f_0 * ip_32[k]
                   + f_1 * ks0_14[k]
                   - f_2 * ks1_17[k]
                   + pb_z[k] * kp_48[k];
    }
}

auto
compute_prim_kd_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hd0, const size_t hd1,
                                      const size_t ip, const size_t id, const size_t ks0,
                                      const size_t ks1, const size_t kp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / alpha;
    const auto f_11 = 1.5 * beta / (alpha * p);
    const auto f_12 = 3.0 / p;
    const auto f_13 = 2.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 0.5 / p;

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_23 = buffer.data(hd0 + 23);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_26 = buffer.data(hd1 + 26);
    const auto *hd1_29 = buffer.data(hd1 + 29);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);

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

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);
    const auto *ks0_15 = buffer.data(ks0 + 15);
    const auto *ks0_17 = buffer.data(ks0 + 17);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_16 = buffer.data(ks1 + 16);
    const auto *ks1_19 = buffer.data(ks1 + 19);
    const auto *ks1_21 = buffer.data(ks1 + 21);
    const auto *ks1_22 = buffer.data(ks1 + 22);
    const auto *ks1_23 = buffer.data(ks1 + 23);
    const auto *ks1_24 = buffer.data(ks1 + 24);
    const auto *ks1_26 = buffer.data(ks1 + 26);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_21 = buffer.data(kp + 21);
    const auto *kp_22 = buffer.data(kp + 22);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_24 = buffer.data(kp + 24);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_28 = buffer.data(kp + 28);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_30 = buffer.data(kp + 30);
    const auto *kp_31 = buffer.data(kp + 31);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, ip_0, id_0, ks0_0, ks1_0, \
                         kp_0, kp_1, kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];

        t_3[k] = pa_y[k] * id_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, hd0_0, hd1_0, ip_1, ip_2, id_0, id_1, \
                         id_2, id_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * ip_1[k]
                 + pa_y[k] * id_1[k];

        t_5[k] = pa_z[k] * id_0[k];

        t_6[k] = f_3 * ip_2[k]
                 + pa_z[k] * id_2[k];

        t_7[k] = f_4 * hd0_0[k]
                 - f_5 * hd1_0[k]
                 + pa_y[k] * id_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_x, pa_z, pb_z, hd0_0, hd0_4, hd1_0, hd1_6, id_4, \
                         id_6, ks0_1, ks1_3, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_6 * hd0_4[k]
                 - f_7 * hd1_6[k]
                 + pa_x[k] * id_6[k];

        t_9[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_3[k]
                 + pb_z[k] * kp_4[k];

        t_10[k] = f_4 * hd0_0[k]
                  - f_5 * hd1_0[k]
                  + pa_z[k] * id_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, pb_y, hd0_1, hd0_6, hd1_3, hd1_8, id_5, \
                         id_8, ks0_2, ks1_4, kp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_1 * ks0_2[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_5[k];

        t_12[k] = f_6 * hd0_6[k]
                  - f_7 * hd1_8[k]
                  + pa_x[k] * id_8[k];

        t_13[k] = f_8 * hd0_1[k]
                  - f_9 * hd1_3[k]
                  + pa_y[k] * id_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_y, pb_z, hd0_8, hd1_10, id_7, id_10, \
                         ks0_3, ks1_5, kp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_10 * hd0_8[k]
                  - f_11 * hd1_10[k]
                  + pa_x[k] * id_10[k];

        t_15[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];

        t_16[k] = pa_y[k] * id_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_z, pb_y, hd0_2, hd0_11, hd1_4, hd1_13, \
                         id_7, id_13, ks0_4, ks1_7, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * hd0_2[k]
                  - f_9 * hd1_4[k]
                  + pa_z[k] * id_7[k];

        t_18[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_7[k]
                  + pb_y[k] * kp_8[k];

        t_19[k] = f_10 * hd0_11[k]
                  - f_11 * hd1_13[k]
                  + pa_x[k] * id_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_x, pa_y, pb_z, hd0_3, hd0_12, hd1_5, hd1_14, \
                         id_9, id_15, ks0_5, ks1_8, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * hd0_3[k]
                  - f_11 * hd1_5[k]
                  + pa_y[k] * id_9[k];

        t_21[k] = f_8 * hd0_12[k]
                  - f_9 * hd1_14[k]
                  + pa_x[k] * id_15[k];

        t_22[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_8[k]
                  + pb_z[k] * kp_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pa_y, hd0_5, hd0_13, hd0_14, hd1_7, \
                         hd1_15, hd1_16, id_11, id_12, id_17, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * hd0_5[k]
                  - f_5 * hd1_7[k]
                  + pa_y[k] * id_11[k];

        t_24[k] = f_8 * hd0_13[k]
                  - f_9 * hd1_15[k]
                  + pa_x[k] * id_17[k];

        t_25[k] = f_8 * hd0_14[k]
                  - f_9 * hd1_16[k]
                  + pa_x[k] * id_18[k];

        t_26[k] = pa_y[k] * id_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_z, pb_y, hd0_5, hd0_15, hd1_7, hd1_17, \
                         id_12, id_21, ks0_6, ks1_11, kp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_10 * hd0_5[k]
                  - f_11 * hd1_7[k]
                  + pa_z[k] * id_12[k];

        t_28[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_11[k]
                  + pb_y[k] * kp_11[k];

        t_29[k] = f_8 * hd0_15[k]
                  - f_9 * hd1_17[k]
                  + pa_x[k] * id_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_x, pa_y, pb_z, hd0_7, hd0_16, hd1_9, hd1_19, \
                         id_14, id_22, ks0_7, ks1_12, kp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_6 * hd0_7[k]
                  - f_7 * hd1_9[k]
                  + pa_y[k] * id_14[k];

        t_31[k] = f_4 * hd0_16[k]
                  - f_5 * hd1_19[k]
                  + pa_x[k] * id_22[k];

        t_32[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_12[k]
                  + pb_z[k] * kp_13[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_x, pa_y, hd0_9, hd0_18, hd0_19, hd1_11, hd1_22, \
                         hd1_23, id_16, id_23, id_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_8 * hd0_9[k]
                  - f_9 * hd1_11[k]
                  + pa_y[k] * id_16[k];

        t_34[k] = f_4 * hd0_18[k]
                  - f_5 * hd1_22[k]
                  + pa_x[k] * id_23[k];

        t_35[k] = f_4 * hd0_19[k]
                  - f_5 * hd1_23[k]
                  + pa_x[k] * id_24[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pa_y, hd0_10, hd0_20, hd0_21, hd1_12, \
                         hd1_24, hd1_25, id_19, id_20, id_25, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_4 * hd0_10[k]
                  - f_5 * hd1_12[k]
                  + pa_y[k] * id_19[k];

        t_37[k] = f_4 * hd0_20[k]
                  - f_5 * hd1_24[k]
                  + pa_x[k] * id_25[k];

        t_38[k] = f_4 * hd0_21[k]
                  - f_5 * hd1_25[k]
                  + pa_x[k] * id_26[k];

        t_39[k] = pa_y[k] * id_20[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_x, pa_z, pb_y, hd0_10, hd0_23, hd1_12, hd1_29, \
                         id_20, id_27, ks0_8, ks1_16, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_6 * hd0_10[k]
                  - f_7 * hd1_12[k]
                  + pa_z[k] * id_20[k];

        t_41[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_16[k]
                  + pb_y[k] * kp_14[k];

        t_42[k] = f_4 * hd0_23[k]
                  - f_5 * hd1_29[k]
                  + pa_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, t_48, t_49, pa_x, ip_13, id_28, id_29, \
                         id_32, id_33, id_34, id_35, id_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * ip_13[k]
                  + pa_x[k] * id_28[k];

        t_44[k] = pa_x[k] * id_29[k];

        t_45[k] = pa_x[k] * id_32[k];

        t_46[k] = pa_x[k] * id_33[k];

        t_47[k] = pa_x[k] * id_34[k];

        t_48[k] = pa_x[k] * id_35[k];

        t_49[k] = pa_x[k] * id_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_x, pb_x, ip_24, id_37, id_39, id_41, \
                         ks0_11, ks1_19, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_x[k] * id_37[k];

        t_51[k] = f_3 * ip_24[k]
                  + pa_x[k] * id_39[k];

        t_52[k] = pa_x[k] * id_41[k];

        t_53[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_19[k]
                  + pb_x[k] * kp_16[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_z, pb_y, pb_z, ip_14, ip_16, id_29, \
                         ks0_11, ks1_19, kp_17, kp_18, kp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_0 * ip_14[k]
                  + f_1 * ks0_11[k]
                  - f_2 * ks1_19[k]
                  + pb_y[k] * kp_17[k];

        t_55[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_19[k]
                  + pb_z[k] * kp_18[k];

        t_56[k] = pa_z[k] * id_29[k];

        t_57[k] = f_12 * ip_16[k]
                  + pb_y[k] * kp_19[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_z, pb_x, hd0_16, hd1_19, ip_15, id_30, id_31, \
                         ks0_12, ks1_21, kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_3 * ip_15[k]
                  + pa_z[k] * id_30[k];

        t_59[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_21[k]
                  + pb_x[k] * kp_20[k];

        t_60[k] = f_4 * hd0_16[k]
                  - f_5 * hd1_19[k]
                  + pa_z[k] * id_31[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pa_y, pb_x, pb_y, hd0_19, hd1_23, ip_18, id_33, \
                         ks0_13, ks1_22, kp_21, kp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_13 * ip_18[k]
                  + pb_y[k] * kp_21[k];

        t_62[k] = f_6 * hd0_19[k]
                  - f_7 * hd1_23[k]
                  + pa_y[k] * id_33[k];

        t_63[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_22[k]
                  + pb_x[k] * kp_22[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_y, pa_z, pb_y, hd0_17, hd0_21, hd1_21, hd1_25, \
                         ip_20, id_32, id_35, kp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_8 * hd0_17[k]
                  - f_9 * hd1_21[k]
                  + pa_z[k] * id_32[k];

        t_65[k] = f_14 * ip_20[k]
                  + pb_y[k] * kp_23[k];

        t_66[k] = f_10 * hd0_21[k]
                  - f_11 * hd1_25[k]
                  + pa_y[k] * id_35[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pa_z, pb_x, pb_y, hd0_18, hd1_22, ip_22, id_34, \
                         ks0_14, ks1_23, kp_24, kp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_23[k]
                  + pb_x[k] * kp_24[k];

        t_68[k] = f_10 * hd0_18[k]
                  - f_11 * hd1_22[k]
                  + pa_z[k] * id_34[k];

        t_69[k] = f_15 * ip_22[k]
                  + pb_y[k] * kp_25[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_y, pa_z, pb_x, hd0_20, hd0_22, hd1_24, hd1_26, \
                         id_36, id_37, ks0_15, ks1_24, kp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_8 * hd0_22[k]
                  - f_9 * hd1_26[k]
                  + pa_y[k] * id_37[k];

        t_71[k] = f_1 * ks0_15[k]
                  - f_2 * ks1_24[k]
                  + pb_x[k] * kp_26[k];

        t_72[k] = f_6 * hd0_20[k]
                  - f_7 * hd1_24[k]
                  + pa_z[k] * id_36[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_y, pb_y, hd0_23, hd1_29, ip_23, ip_25, \
                         ip_26, id_38, id_40, kp_27, kp_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_3 * ip_23[k]
                  + pb_y[k] * kp_27[k];

        t_74[k] = f_4 * hd0_23[k]
                  - f_5 * hd1_29[k]
                  + pa_y[k] * id_38[k];

        t_75[k] = f_3 * ip_25[k]
                  + pa_y[k] * id_40[k];

        t_76[k] = f_16 * ip_26[k]
                  + pb_y[k] * kp_28[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pb_x, pb_y, pb_z, ip_26, id_41, ks0_17, \
                         ks1_26, kp_29, kp_30, kp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_y[k] * id_41[k];

        t_78[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_26[k]
                  + pb_x[k] * kp_29[k];

        t_79[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_26[k]
                  + pb_y[k] * kp_30[k];

        t_80[k] = f_0 * ip_26[k]
                  + f_1 * ks0_17[k]
                  - f_2 * ks1_26[k]
                  + pb_z[k] * kp_31[k];
    }
}

auto
compute_prim_kd_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hd0, const size_t hd1,
                                      const size_t ip, const size_t id, const size_t ks0,
                                      const size_t ks1, const size_t kp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / alpha;
    const auto f_6 = 2.0 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 1.5 / alpha;
    const auto f_10 = 1.5 * beta / (alpha * p);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_23 = buffer.data(hd1 + 23);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_17 = buffer.data(ip + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_62 = buffer.data(id + 62);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_68 = buffer.data(id + 68);
    const auto *id_69 = buffer.data(id + 69);
    const auto *id_71 = buffer.data(id + 71);
    const auto *id_74 = buffer.data(id + 74);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);
    const auto *ks0_15 = buffer.data(ks0 + 15);
    const auto *ks0_17 = buffer.data(ks0 + 17);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);
    const auto *ks1_15 = buffer.data(ks1 + 15);
    const auto *ks1_17 = buffer.data(ks1 + 17);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, ip_0, id_0, ks0_0, ks1_0, \
                         kp_0, kp_1, kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];

        t_3[k] = pa_y[k] * id_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, hd0_0, hd0_4, hd1_0, hd1_4, id_0, \
                         id_3, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * id_0[k];

        t_5[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_y[k] * id_3[k];

        t_6[k] = f_5 * hd0_4[k]
                 - f_6 * hd1_4[k]
                 + pa_x[k] * id_10[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, hd0_0, hd1_0, id_6, ks0_1, ks0_2, \
                         ks1_1, ks1_2, kp_3, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_1[k]
                 + pb_z[k] * kp_3[k];

        t_8[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_z[k] * id_6[k];

        t_9[k] = f_1 * ks0_2[k]
                 - f_2 * ks1_2[k]
                 + pb_y[k] * kp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, hd0_1, hd0_6, hd0_8, hd1_1, hd1_6, \
                         hd1_8, id_9, id_16, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd0_6[k]
                  - f_6 * hd1_6[k]
                  + pa_x[k] * id_16[k];

        t_11[k] = f_7 * hd0_1[k]
                  - f_8 * hd1_1[k]
                  + pa_y[k] * id_9[k];

        t_12[k] = f_9 * hd0_8[k]
                  - f_10 * hd1_8[k]
                  + pa_x[k] * id_18[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, hd0_2, hd1_2, id_14, ks0_3, \
                         ks0_4, ks1_3, ks1_4, kp_5, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];

        t_14[k] = f_7 * hd0_2[k]
                  - f_8 * hd1_2[k]
                  + pa_z[k] * id_14[k];

        t_15[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, hd0_3, hd0_10, hd0_11, hd1_3, hd1_11, \
                         hd1_12, id_17, id_28, id_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_9 * hd0_10[k]
                  - f_10 * hd1_11[k]
                  + pa_x[k] * id_28[k];

        t_17[k] = f_9 * hd0_3[k]
                  - f_10 * hd1_3[k]
                  + pa_y[k] * id_17[k];

        t_18[k] = f_7 * hd0_11[k]
                  - f_8 * hd1_12[k]
                  + pa_x[k] * id_30[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_y, pb_z, hd0_5, hd1_5, id_26, ks0_5, \
                         ks0_6, ks1_5, ks1_6, kp_7, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];

        t_20[k] = f_9 * hd0_5[k]
                  - f_10 * hd1_5[k]
                  + pa_z[k] * id_26[k];

        t_21[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pa_y, hd0_7, hd0_12, hd0_13, hd1_7, hd1_15, \
                         hd1_16, id_29, id_43, id_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_7 * hd0_12[k]
                  - f_8 * hd1_15[k]
                  + pa_x[k] * id_43[k];

        t_23[k] = f_5 * hd0_7[k]
                  - f_6 * hd1_7[k]
                  + pa_y[k] * id_29[k];

        t_24[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_16[k]
                  + pa_x[k] * id_45[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_z, pb_y, pb_z, hd0_9, hd1_10, id_41, ks0_7, \
                         ks0_8, ks1_7, ks1_8, kp_9, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];

        t_26[k] = f_5 * hd0_9[k]
                  - f_6 * hd1_10[k]
                  + pa_z[k] * id_41[k];

        t_27[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_x, hd0_20, hd1_23, id_54, id_56, \
                         id_74, ks0_11, ks1_11, kp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_23[k]
                  + pa_x[k] * id_54[k];

        t_29[k] = pa_x[k] * id_56[k];

        t_30[k] = pa_x[k] * id_74[k];

        t_31[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_11[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_z, pb_y, pb_z, ip_10, id_56, ks0_11, ks1_11, \
                         kp_12, kp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * ip_10[k]
                  + f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_y[k] * kp_12[k];

        t_33[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_z[k] * kp_13[k];

        t_34[k] = pa_z[k] * id_56[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pa_z, pb_x, hd0_13, hd0_16, hd1_16, hd1_19, \
                         id_59, id_63, ks0_12, ks1_12, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_14[k];

        t_36[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_16[k]
                  + pa_z[k] * id_59[k];

        t_37[k] = f_5 * hd0_16[k]
                  - f_6 * hd1_19[k]
                  + pa_y[k] * id_63[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pa_z, pb_x, hd0_14, hd0_18, hd1_17, hd1_21, \
                         id_62, id_66, ks0_13, ks1_13, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_15[k];

        t_39[k] = f_7 * hd0_14[k]
                  - f_8 * hd1_17[k]
                  + pa_z[k] * id_62[k];

        t_40[k] = f_9 * hd0_18[k]
                  - f_10 * hd1_21[k]
                  + pa_y[k] * id_66[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_y, pa_z, pb_x, hd0_15, hd0_19, hd1_18, hd1_22, \
                         id_65, id_69, ks0_14, ks1_14, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_16[k];

        t_42[k] = f_9 * hd0_15[k]
                  - f_10 * hd1_18[k]
                  + pa_z[k] * id_65[k];

        t_43[k] = f_7 * hd0_19[k]
                  - f_8 * hd1_22[k]
                  + pa_y[k] * id_69[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, pa_z, pb_x, hd0_17, hd0_20, hd1_20, hd1_23, \
                         id_68, id_71, ks0_15, ks1_15, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_1 * ks0_15[k]
                  - f_2 * ks1_15[k]
                  + pb_x[k] * kp_17[k];

        t_45[k] = f_5 * hd0_17[k]
                  - f_6 * hd1_20[k]
                  + pa_z[k] * id_68[k];

        t_46[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_23[k]
                  + pa_y[k] * id_71[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pb_x, pb_y, pb_z, ip_17, id_74, ks0_17, \
                         ks1_17, kp_18, kp_19, kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pa_y[k] * id_74[k];

        t_48[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_x[k] * kp_18[k];

        t_49[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_y[k] * kp_19[k];

        t_50[k] = f_0 * ip_17[k]
                  + f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hd0, const size_t hd1,
                                      const size_t ip, const size_t id, const size_t ks0,
                                      const size_t ks1, const size_t kp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / alpha;
    const auto f_11 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_23 = buffer.data(hd0 + 23);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_29 = buffer.data(hd1 + 29);
    const auto *hd1_32 = buffer.data(hd1 + 32);
    const auto *hd1_33 = buffer.data(hd1 + 33);
    const auto *hd1_35 = buffer.data(hd1 + 35);
    const auto *hd1_36 = buffer.data(hd1 + 36);
    const auto *hd1_38 = buffer.data(hd1 + 38);
    const auto *hd1_41 = buffer.data(hd1 + 41);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);

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
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);
    const auto *ks0_15 = buffer.data(ks0 + 15);
    const auto *ks0_17 = buffer.data(ks0 + 17);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);
    const auto *ks1_15 = buffer.data(ks1 + 15);
    const auto *ks1_17 = buffer.data(ks1 + 17);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, ip_0, id_0, ks0_0, ks1_0, \
                         kp_0, kp_1, kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];

        t_3[k] = pa_y[k] * id_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, hd0_0, hd1_0, ip_1, ip_2, \
                         id_0, id_1, id_2, id_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * ip_1[k]
                 + pa_y[k] * id_1[k];

        t_5[k] = pa_y[k] * id_2[k];

        t_6[k] = pa_z[k] * id_0[k];

        t_7[k] = pa_z[k] * id_1[k];

        t_8[k] = f_3 * ip_2[k]
                 + pa_z[k] * id_2[k];

        t_9[k] = f_4 * hd0_0[k]
                 - f_5 * hd1_0[k]
                 + pa_y[k] * id_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, pa_z, pb_z, hd0_4, hd1_8, id_4, \
                         id_6, id_8, ks0_1, ks1_1, kp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_6 * hd0_4[k]
                  - f_7 * hd1_8[k]
                  + pa_x[k] * id_8[k];

        t_11[k] = f_1 * ks0_1[k]
                  - f_2 * ks1_1[k]
                  + pb_z[k] * kp_3[k];

        t_12[k] = pa_z[k] * id_4[k];

        t_13[k] = pa_y[k] * id_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_y, hd0_0, hd0_6, hd1_0, hd1_12, \
                         id_5, id_12, ks0_2, ks1_2, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_4 * hd0_0[k]
                  - f_5 * hd1_0[k]
                  + pa_z[k] * id_5[k];

        t_15[k] = f_1 * ks0_2[k]
                  - f_2 * ks1_2[k]
                  + pb_y[k] * kp_4[k];

        t_16[k] = f_6 * hd0_6[k]
                  - f_7 * hd1_12[k]
                  + pa_x[k] * id_12[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, pb_z, hd0_1, hd0_8, hd1_3, hd1_14, \
                         id_7, id_14, ks0_3, ks1_3, kp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * hd0_1[k]
                  - f_9 * hd1_3[k]
                  + pa_y[k] * id_7[k];

        t_18[k] = f_10 * hd0_8[k]
                  - f_11 * hd1_14[k]
                  + pa_x[k] * id_14[k];

        t_19[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pa_z, ip_3, ip_4, id_7, id_8, \
                         id_9, id_11, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * id_7[k];

        t_21[k] = pa_z[k] * id_8[k];

        t_22[k] = f_3 * ip_3[k]
                  + pa_z[k] * id_9[k];

        t_23[k] = f_3 * ip_4[k]
                  + pa_y[k] * id_11[k];

        t_24[k] = pa_y[k] * id_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_x, pa_z, pb_y, hd0_2, hd0_11, hd1_5, hd1_19, \
                         id_10, id_19, ks0_4, ks1_4, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * hd0_2[k]
                  - f_9 * hd1_5[k]
                  + pa_z[k] * id_10[k];

        t_26[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];

        t_27[k] = f_10 * hd0_11[k]
                  - f_11 * hd1_19[k]
                  + pa_x[k] * id_19[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_x, pa_y, pb_z, hd0_3, hd0_12, hd1_7, hd1_21, \
                         id_13, id_21, ks0_5, ks1_5, kp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_10 * hd0_3[k]
                  - f_11 * hd1_7[k]
                  + pa_y[k] * id_13[k];

        t_29[k] = f_8 * hd0_12[k]
                  - f_9 * hd1_21[k]
                  + pa_x[k] * id_21[k];

        t_30[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pa_z, hd0_5, hd1_10, ip_5, id_13, \
                         id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pa_z[k] * id_13[k];

        t_32[k] = pa_z[k] * id_14[k];

        t_33[k] = f_3 * ip_5[k]
                  + pa_z[k] * id_15[k];

        t_34[k] = f_4 * hd0_5[k]
                  - f_5 * hd1_10[k]
                  + pa_y[k] * id_16[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pa_y, hd0_13, hd0_14, hd1_22, hd1_23, \
                         ip_6, id_18, id_19, id_24, id_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * hd0_13[k]
                  - f_9 * hd1_22[k]
                  + pa_x[k] * id_24[k];

        t_36[k] = f_8 * hd0_14[k]
                  - f_9 * hd1_23[k]
                  + pa_x[k] * id_25[k];

        t_37[k] = f_3 * ip_6[k]
                  + pa_y[k] * id_18[k];

        t_38[k] = pa_y[k] * id_19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_x, pa_z, pb_y, hd0_5, hd0_15, hd1_10, hd1_25, \
                         id_17, id_29, ks0_6, ks1_6, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * hd0_5[k]
                  - f_11 * hd1_10[k]
                  + pa_z[k] * id_17[k];

        t_40[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];

        t_41[k] = f_8 * hd0_15[k]
                  - f_9 * hd1_25[k]
                  + pa_x[k] * id_29[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_x, pa_y, pb_z, hd0_7, hd0_16, hd1_13, hd1_27, \
                         id_20, id_31, ks0_7, ks1_7, kp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_6 * hd0_7[k]
                  - f_7 * hd1_13[k]
                  + pa_y[k] * id_20[k];

        t_43[k] = f_4 * hd0_16[k]
                  - f_5 * hd1_27[k]
                  + pa_x[k] * id_31[k];

        t_44[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_y, pa_z, hd0_9, hd1_16, ip_7, id_20, \
                         id_21, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_z[k] * id_20[k];

        t_46[k] = pa_z[k] * id_21[k];

        t_47[k] = f_3 * ip_7[k]
                  + pa_z[k] * id_22[k];

        t_48[k] = f_8 * hd0_9[k]
                  - f_9 * hd1_16[k]
                  + pa_y[k] * id_23[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_x, pa_y, hd0_10, hd0_18, hd0_19, hd1_17, hd1_32, \
                         hd1_33, id_26, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_4 * hd0_18[k]
                  - f_5 * hd1_32[k]
                  + pa_x[k] * id_32[k];

        t_50[k] = f_4 * hd0_19[k]
                  - f_5 * hd1_33[k]
                  + pa_x[k] * id_33[k];

        t_51[k] = f_4 * hd0_10[k]
                  - f_5 * hd1_17[k]
                  + pa_y[k] * id_26[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pa_y, hd0_20, hd0_21, hd1_35, hd1_36, \
                         ip_8, id_28, id_29, id_34, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_4 * hd0_20[k]
                  - f_5 * hd1_35[k]
                  + pa_x[k] * id_34[k];

        t_53[k] = f_4 * hd0_21[k]
                  - f_5 * hd1_36[k]
                  + pa_x[k] * id_35[k];

        t_54[k] = f_3 * ip_8[k]
                  + pa_y[k] * id_28[k];

        t_55[k] = pa_y[k] * id_29[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pa_x, pa_z, pb_y, hd0_10, hd0_23, hd1_17, hd1_41, \
                         id_27, id_37, ks0_8, ks1_8, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_6 * hd0_10[k]
                  - f_7 * hd1_17[k]
                  + pa_z[k] * id_27[k];

        t_57[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];

        t_58[k] = f_4 * hd0_23[k]
                  - f_5 * hd1_41[k]
                  + pa_x[k] * id_37[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, t_63, pa_x, pa_z, ip_9, ip_12, ip_13, id_30, \
                         id_38, id_39, id_43, id_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_3 * ip_9[k]
                  + pa_x[k] * id_38[k];

        t_60[k] = pa_x[k] * id_39[k];

        t_61[k] = pa_z[k] * id_30[k];

        t_62[k] = f_3 * ip_12[k]
                  + pa_x[k] * id_43[k];

        t_63[k] = f_3 * ip_13[k]
                  + pa_x[k] * id_46[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, ip_14, ip_15, id_49, id_54, \
                         id_56, ks0_11, ks1_11, kp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * ip_14[k]
                  + pa_x[k] * id_49[k];

        t_65[k] = f_3 * ip_15[k]
                  + pa_x[k] * id_54[k];

        t_66[k] = pa_x[k] * id_56[k];

        t_67[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_11[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_z, pb_y, pb_z, ip_10, id_38, id_39, \
                         ks0_11, ks1_11, kp_12, kp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * ip_10[k]
                  + f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_y[k] * kp_12[k];

        t_69[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_z[k] * kp_13[k];

        t_70[k] = pa_z[k] * id_38[k];

        t_71[k] = pa_z[k] * id_39[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_z, pb_x, hd0_16, hd1_27, ip_11, id_40, id_41, \
                         ks0_12, ks1_12, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_3 * ip_11[k]
                  + pa_z[k] * id_40[k];

        t_73[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_14[k];

        t_74[k] = f_4 * hd0_16[k]
                  - f_5 * hd1_27[k]
                  + pa_z[k] * id_41[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_y, pa_z, pb_x, hd0_17, hd0_19, hd1_29, hd1_33, \
                         id_44, id_45, ks0_13, ks1_13, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_6 * hd0_19[k]
                  - f_7 * hd1_33[k]
                  + pa_y[k] * id_45[k];

        t_76[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_15[k];

        t_77[k] = f_8 * hd0_17[k]
                  - f_9 * hd1_29[k]
                  + pa_z[k] * id_44[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pa_y, pa_z, pb_x, hd0_18, hd0_21, hd1_32, hd1_36, \
                         id_47, id_48, ks0_14, ks1_14, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * hd0_21[k]
                  - f_11 * hd1_36[k]
                  + pa_y[k] * id_48[k];

        t_79[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_16[k];

        t_80[k] = f_10 * hd0_18[k]
                  - f_11 * hd1_32[k]
                  + pa_z[k] * id_47[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_y, pa_z, pb_x, hd0_20, hd0_22, hd1_35, hd1_38, \
                         id_50, id_51, ks0_15, ks1_15, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_8 * hd0_22[k]
                  - f_9 * hd1_38[k]
                  + pa_y[k] * id_51[k];

        t_82[k] = f_1 * ks0_15[k]
                  - f_2 * ks1_15[k]
                  + pb_x[k] * kp_17[k];

        t_83[k] = f_6 * hd0_20[k]
                  - f_7 * hd1_35[k]
                  + pa_z[k] * id_50[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_y, pb_x, hd0_23, hd1_41, ip_16, id_53, \
                         id_55, id_56, ks0_17, ks1_17, kp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * hd0_23[k]
                  - f_5 * hd1_41[k]
                  + pa_y[k] * id_53[k];

        t_85[k] = f_3 * ip_16[k]
                  + pa_y[k] * id_55[k];

        t_86[k] = pa_y[k] * id_56[k];

        t_87[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_x[k] * kp_18[k];
    }

#pragma omp simd aligned(t_88, t_89, pb_y, pb_z, ip_17, ks0_17, ks1_17, kp_19, \
                         kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_y[k] * kp_19[k];

        t_89[k] = f_0 * ip_17[k]
                  + f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hd0, const size_t hd1,
                                      const size_t ip, const size_t id, const size_t ks0,
                                      const size_t ks1, const size_t kp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / alpha;
    const auto f_11 = 1.5 * beta / (alpha * p);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_23 = buffer.data(hd0 + 23);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_27 = buffer.data(hd0 + 27);
    const auto *hd0_29 = buffer.data(hd0 + 29);
    const auto *hd0_32 = buffer.data(hd0 + 32);
    const auto *hd0_33 = buffer.data(hd0 + 33);
    const auto *hd0_35 = buffer.data(hd0 + 35);
    const auto *hd0_36 = buffer.data(hd0 + 36);
    const auto *hd0_38 = buffer.data(hd0 + 38);
    const auto *hd0_41 = buffer.data(hd0 + 41);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_26 = buffer.data(hd1 + 26);
    const auto *hd1_29 = buffer.data(hd1 + 29);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);

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

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);
    const auto *ks0_15 = buffer.data(ks0 + 15);
    const auto *ks0_17 = buffer.data(ks0 + 17);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);
    const auto *ks1_15 = buffer.data(ks1 + 15);
    const auto *ks1_17 = buffer.data(ks1 + 17);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, ip_0, id_0, ks0_0, ks1_0, \
                         kp_0, kp_1, kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];

        t_3[k] = pa_y[k] * id_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, hd0_0, hd1_0, ip_1, ip_2, id_0, id_1, \
                         id_2, id_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * ip_1[k]
                 + pa_y[k] * id_1[k];

        t_5[k] = pa_z[k] * id_0[k];

        t_6[k] = f_3 * ip_2[k]
                 + pa_z[k] * id_2[k];

        t_7[k] = f_4 * hd0_0[k]
                 - f_5 * hd1_0[k]
                 + pa_y[k] * id_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_x, pa_z, pb_z, hd0_0, hd0_8, hd1_0, hd1_6, id_4, \
                         id_6, ks0_1, ks1_1, kp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_6 * hd0_8[k]
                 - f_7 * hd1_6[k]
                 + pa_x[k] * id_6[k];

        t_9[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_1[k]
                 + pb_z[k] * kp_3[k];

        t_10[k] = f_4 * hd0_0[k]
                  - f_5 * hd1_0[k]
                  + pa_z[k] * id_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, pb_y, hd0_3, hd0_12, hd1_3, hd1_8, \
                         id_5, id_8, ks0_2, ks1_2, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_1 * ks0_2[k]
                  - f_2 * ks1_2[k]
                  + pb_y[k] * kp_4[k];

        t_12[k] = f_6 * hd0_12[k]
                  - f_7 * hd1_8[k]
                  + pa_x[k] * id_8[k];

        t_13[k] = f_8 * hd0_3[k]
                  - f_9 * hd1_3[k]
                  + pa_y[k] * id_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_y, pb_z, hd0_14, hd1_10, id_7, id_10, \
                         ks0_3, ks1_3, kp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_10 * hd0_14[k]
                  - f_11 * hd1_10[k]
                  + pa_x[k] * id_10[k];

        t_15[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];

        t_16[k] = pa_y[k] * id_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_z, pb_y, hd0_5, hd0_19, hd1_4, hd1_13, \
                         id_7, id_13, ks0_4, ks1_4, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * hd0_5[k]
                  - f_9 * hd1_4[k]
                  + pa_z[k] * id_7[k];

        t_18[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];

        t_19[k] = f_10 * hd0_19[k]
                  - f_11 * hd1_13[k]
                  + pa_x[k] * id_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_x, pa_y, pb_z, hd0_7, hd0_21, hd1_5, hd1_14, \
                         id_9, id_15, ks0_5, ks1_5, kp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * hd0_7[k]
                  - f_11 * hd1_5[k]
                  + pa_y[k] * id_9[k];

        t_21[k] = f_8 * hd0_21[k]
                  - f_9 * hd1_14[k]
                  + pa_x[k] * id_15[k];

        t_22[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pa_y, hd0_10, hd0_22, hd0_23, hd1_7, \
                         hd1_15, hd1_16, id_11, id_12, id_17, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * hd0_10[k]
                  - f_5 * hd1_7[k]
                  + pa_y[k] * id_11[k];

        t_24[k] = f_8 * hd0_22[k]
                  - f_9 * hd1_15[k]
                  + pa_x[k] * id_17[k];

        t_25[k] = f_8 * hd0_23[k]
                  - f_9 * hd1_16[k]
                  + pa_x[k] * id_18[k];

        t_26[k] = pa_y[k] * id_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_z, pb_y, hd0_10, hd0_25, hd1_7, hd1_17, \
                         id_12, id_21, ks0_6, ks1_6, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_10 * hd0_10[k]
                  - f_11 * hd1_7[k]
                  + pa_z[k] * id_12[k];

        t_28[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];

        t_29[k] = f_8 * hd0_25[k]
                  - f_9 * hd1_17[k]
                  + pa_x[k] * id_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_x, pa_y, pb_z, hd0_13, hd0_27, hd1_9, hd1_19, \
                         id_14, id_22, ks0_7, ks1_7, kp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_6 * hd0_13[k]
                  - f_7 * hd1_9[k]
                  + pa_y[k] * id_14[k];

        t_31[k] = f_4 * hd0_27[k]
                  - f_5 * hd1_19[k]
                  + pa_x[k] * id_22[k];

        t_32[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_x, pa_y, hd0_16, hd0_32, hd0_33, hd1_11, hd1_22, \
                         hd1_23, id_16, id_23, id_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_8 * hd0_16[k]
                  - f_9 * hd1_11[k]
                  + pa_y[k] * id_16[k];

        t_34[k] = f_4 * hd0_32[k]
                  - f_5 * hd1_22[k]
                  + pa_x[k] * id_23[k];

        t_35[k] = f_4 * hd0_33[k]
                  - f_5 * hd1_23[k]
                  + pa_x[k] * id_24[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pa_y, hd0_17, hd0_35, hd0_36, hd1_12, \
                         hd1_24, hd1_25, id_19, id_20, id_25, id_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_4 * hd0_17[k]
                  - f_5 * hd1_12[k]
                  + pa_y[k] * id_19[k];

        t_37[k] = f_4 * hd0_35[k]
                  - f_5 * hd1_24[k]
                  + pa_x[k] * id_25[k];

        t_38[k] = f_4 * hd0_36[k]
                  - f_5 * hd1_25[k]
                  + pa_x[k] * id_26[k];

        t_39[k] = pa_y[k] * id_20[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_x, pa_z, pb_y, hd0_17, hd0_41, hd1_12, hd1_29, \
                         id_20, id_27, ks0_8, ks1_8, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_6 * hd0_17[k]
                  - f_7 * hd1_12[k]
                  + pa_z[k] * id_20[k];

        t_41[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];

        t_42[k] = f_4 * hd0_41[k]
                  - f_5 * hd1_29[k]
                  + pa_x[k] * id_27[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, t_48, t_49, pa_x, ip_9, id_28, id_29, \
                         id_32, id_33, id_34, id_35, id_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * ip_9[k]
                  + pa_x[k] * id_28[k];

        t_44[k] = pa_x[k] * id_29[k];

        t_45[k] = pa_x[k] * id_32[k];

        t_46[k] = pa_x[k] * id_33[k];

        t_47[k] = pa_x[k] * id_34[k];

        t_48[k] = pa_x[k] * id_35[k];

        t_49[k] = pa_x[k] * id_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_x, pb_x, ip_15, id_37, id_39, id_41, \
                         ks0_11, ks1_11, kp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_x[k] * id_37[k];

        t_51[k] = f_3 * ip_15[k]
                  + pa_x[k] * id_39[k];

        t_52[k] = pa_x[k] * id_41[k];

        t_53[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_11[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_z, pb_y, pb_z, ip_10, ip_11, id_29, id_30, \
                         ks0_11, ks1_11, kp_12, kp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_0 * ip_10[k]
                  + f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_y[k] * kp_12[k];

        t_55[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_z[k] * kp_13[k];

        t_56[k] = pa_z[k] * id_29[k];

        t_57[k] = f_3 * ip_11[k]
                  + pa_z[k] * id_30[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_y, pa_z, pb_x, hd0_27, hd0_33, hd1_19, hd1_23, \
                         id_31, id_33, ks0_12, ks1_12, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_14[k];

        t_59[k] = f_4 * hd0_27[k]
                  - f_5 * hd1_19[k]
                  + pa_z[k] * id_31[k];

        t_60[k] = f_6 * hd0_33[k]
                  - f_7 * hd1_23[k]
                  + pa_y[k] * id_33[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pa_y, pa_z, pb_x, hd0_29, hd0_36, hd1_21, hd1_25, \
                         id_32, id_35, ks0_13, ks1_13, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_15[k];

        t_62[k] = f_8 * hd0_29[k]
                  - f_9 * hd1_21[k]
                  + pa_z[k] * id_32[k];

        t_63[k] = f_10 * hd0_36[k]
                  - f_11 * hd1_25[k]
                  + pa_y[k] * id_35[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_y, pa_z, pb_x, hd0_32, hd0_38, hd1_22, hd1_26, \
                         id_34, id_37, ks0_14, ks1_14, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_16[k];

        t_65[k] = f_10 * hd0_32[k]
                  - f_11 * hd1_22[k]
                  + pa_z[k] * id_34[k];

        t_66[k] = f_8 * hd0_38[k]
                  - f_9 * hd1_26[k]
                  + pa_y[k] * id_37[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pa_y, pa_z, pb_x, hd0_35, hd0_41, hd1_24, hd1_29, \
                         id_36, id_38, ks0_15, ks1_15, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_1 * ks0_15[k]
                  - f_2 * ks1_15[k]
                  + pb_x[k] * kp_17[k];

        t_68[k] = f_6 * hd0_35[k]
                  - f_7 * hd1_24[k]
                  + pa_z[k] * id_36[k];

        t_69[k] = f_4 * hd0_41[k]
                  - f_5 * hd1_29[k]
                  + pa_y[k] * id_38[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_y, pb_x, pb_y, ip_16, id_40, id_41, \
                         ks0_17, ks1_17, kp_18, kp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * ip_16[k]
                  + pa_y[k] * id_40[k];

        t_71[k] = pa_y[k] * id_41[k];

        t_72[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_x[k] * kp_18[k];

        t_73[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_y[k] * kp_19[k];
    }

#pragma omp simd aligned(t_74, pb_z, ip_17, ks0_17, ks1_17, kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * ip_17[k]
                  + f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hd0, const size_t hd1,
                                      const size_t ip, const size_t id, const size_t ks0,
                                      const size_t ks1, const size_t kp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / alpha;
    const auto f_6 = 2.0 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 1.5 / alpha;
    const auto f_10 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_17 = buffer.data(ip + 17);

    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_9 = buffer.data(ks0 + 9);
    const auto *ks0_10 = buffer.data(ks0 + 10);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);
    const auto *ks1_15 = buffer.data(ks1 + 15);
    const auto *ks1_17 = buffer.data(ks1 + 17);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, ip_0, ks0_0, ks1_0, kp_0, kp_1, \
                         kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, hd0_0, hd0_4, hd1_0, hd1_4, id_3, \
                         id_6, ks0_1, ks1_1, kp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_y[k] * id_3[k];

        t_4[k] = f_5 * hd0_4[k]
                 - f_6 * hd1_4[k]
                 + pa_x[k] * id_6[k];

        t_5[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_1[k]
                 + pb_z[k] * kp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, hd0_0, hd0_6, hd1_0, hd1_6, id_4, \
                         id_10, ks0_2, ks1_2, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_z[k] * id_4[k];

        t_7[k] = f_1 * ks0_2[k]
                 - f_2 * ks1_2[k]
                 + pb_y[k] * kp_4[k];

        t_8[k] = f_5 * hd0_6[k]
                 - f_6 * hd1_6[k]
                 + pa_x[k] * id_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, hd0_1, hd0_8, hd1_1, hd1_8, id_5, \
                         id_12, ks0_3, ks1_3, kp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * hd0_1[k]
                 - f_8 * hd1_1[k]
                 + pa_y[k] * id_5[k];

        t_10[k] = f_9 * hd0_8[k]
                  - f_10 * hd1_8[k]
                  + pa_x[k] * id_12[k];

        t_11[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, hd0_2, hd0_10, hd1_2, hd1_10, \
                         id_8, id_16, ks0_4, ks1_4, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * hd0_2[k]
                  - f_8 * hd1_2[k]
                  + pa_z[k] * id_8[k];

        t_13[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];

        t_14[k] = f_9 * hd0_10[k]
                  - f_10 * hd1_10[k]
                  + pa_x[k] * id_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, hd0_3, hd0_11, hd1_3, hd1_11, \
                         id_11, id_18, ks0_5, ks1_5, kp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_9 * hd0_3[k]
                  - f_10 * hd1_3[k]
                  + pa_y[k] * id_11[k];

        t_16[k] = f_7 * hd0_11[k]
                  - f_8 * hd1_11[k]
                  + pa_x[k] * id_18[k];

        t_17[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, hd0_5, hd0_12, hd1_5, hd1_12, \
                         id_14, id_22, ks0_6, ks1_6, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_9 * hd0_5[k]
                  - f_10 * hd1_5[k]
                  + pa_z[k] * id_14[k];

        t_19[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];

        t_20[k] = f_7 * hd0_12[k]
                  - f_8 * hd1_12[k]
                  + pa_x[k] * id_22[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_z, hd0_7, hd0_13, hd1_7, hd1_13, \
                         id_17, id_23, ks0_7, ks1_7, kp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * hd0_7[k]
                  - f_6 * hd1_7[k]
                  + pa_y[k] * id_17[k];

        t_22[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_13[k]
                  + pa_x[k] * id_23[k];

        t_23[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_z, pb_y, hd0_9, hd0_20, hd1_9, hd1_20, \
                         id_20, id_24, ks0_8, ks1_8, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * hd0_9[k]
                  - f_6 * hd1_9[k]
                  + pa_z[k] * id_20[k];

        t_25[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];

        t_26[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pa_x[k] * id_24[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pb_x, pb_y, ip_10, id_26, id_41, ks0_9, \
                         ks1_11, kp_11, kp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_x[k] * id_26[k];

        t_28[k] = pa_x[k] * id_41[k];

        t_29[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_11[k];

        t_30[k] = f_0 * ip_10[k]
                  + f_1 * ks0_9[k]
                  - f_2 * ks1_11[k]
                  + pb_y[k] * kp_12[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_z, pb_x, pb_z, hd0_13, hd1_13, id_28, ks0_9, \
                         ks0_10, ks1_11, ks1_12, kp_13, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_11[k]
                  + pb_z[k] * kp_13[k];

        t_32[k] = f_1 * ks0_10[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_14[k];

        t_33[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_13[k]
                  + pa_z[k] * id_28[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pa_z, pb_x, hd0_14, hd0_16, hd1_14, hd1_16, \
                         id_30, id_31, ks0_11, ks1_13, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_5 * hd0_16[k]
                  - f_6 * hd1_16[k]
                  + pa_y[k] * id_31[k];

        t_35[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_15[k];

        t_36[k] = f_7 * hd0_14[k]
                  - f_8 * hd1_14[k]
                  + pa_z[k] * id_30[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, hd0_15, hd0_18, hd1_15, hd1_18, \
                         id_33, id_34, ks0_12, ks1_14, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * hd0_18[k]
                  - f_10 * hd1_18[k]
                  + pa_y[k] * id_34[k];

        t_38[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_16[k];

        t_39[k] = f_9 * hd0_15[k]
                  - f_10 * hd1_15[k]
                  + pa_z[k] * id_33[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_x, hd0_17, hd0_19, hd1_17, hd1_19, \
                         id_36, id_37, ks0_13, ks1_15, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_7 * hd0_19[k]
                  - f_8 * hd1_19[k]
                  + pa_y[k] * id_37[k];

        t_41[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_15[k]
                  + pb_x[k] * kp_17[k];

        t_42[k] = f_5 * hd0_17[k]
                  - f_6 * hd1_17[k]
                  + pa_z[k] * id_36[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_y, pb_x, pb_y, hd0_20, hd1_20, id_38, \
                         id_41, ks0_14, ks1_17, kp_18, kp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pa_y[k] * id_38[k];

        t_44[k] = pa_y[k] * id_41[k];

        t_45[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_17[k]
                  + pb_x[k] * kp_18[k];

        t_46[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_17[k]
                  + pb_y[k] * kp_19[k];
    }

#pragma omp simd aligned(t_47, pb_z, ip_17, ks0_14, ks1_17, kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * ip_17[k]
                  + f_1 * ks0_14[k]
                  - f_2 * ks1_17[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hd0, const size_t hd1,
                                      const size_t ip, const size_t id, const size_t ks0,
                                      const size_t ks1, const size_t kp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / alpha;
    const auto f_6 = 2.0 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 1.5 / alpha;
    const auto f_10 = 1.5 * beta / (alpha * p);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_28 = buffer.data(hd1 + 28);
    const auto *hd1_29 = buffer.data(hd1 + 29);
    const auto *hd1_32 = buffer.data(hd1 + 32);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_17 = buffer.data(ip + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_44 = buffer.data(id + 44);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);
    const auto *ks0_15 = buffer.data(ks0 + 15);
    const auto *ks0_17 = buffer.data(ks0 + 17);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);
    const auto *ks1_15 = buffer.data(ks1 + 15);
    const auto *ks1_17 = buffer.data(ks1 + 17);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, ip_0, id_0, ks0_0, ks1_0, \
                         kp_0, kp_1, kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];

        t_3[k] = pa_y[k] * id_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, hd0_0, hd0_4, hd1_0, hd1_6, id_0, \
                         id_3, id_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * id_0[k];

        t_5[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_y[k] * id_3[k];

        t_6[k] = f_5 * hd0_4[k]
                 - f_6 * hd1_6[k]
                 + pa_x[k] * id_7[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, hd0_0, hd1_0, id_4, ks0_1, ks0_2, \
                         ks1_1, ks1_2, kp_3, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_1[k]
                 + pb_z[k] * kp_3[k];

        t_8[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_z[k] * id_4[k];

        t_9[k] = f_1 * ks0_2[k]
                 - f_2 * ks1_2[k]
                 + pb_y[k] * kp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, hd0_1, hd0_6, hd0_8, hd1_3, hd1_10, \
                         hd1_12, id_6, id_11, id_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd0_6[k]
                  - f_6 * hd1_10[k]
                  + pa_x[k] * id_11[k];

        t_11[k] = f_7 * hd0_1[k]
                  - f_8 * hd1_3[k]
                  + pa_y[k] * id_6[k];

        t_12[k] = f_9 * hd0_8[k]
                  - f_10 * hd1_12[k]
                  + pa_x[k] * id_13[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, hd0_2, hd1_4, id_9, ks0_3, ks0_4, \
                         ks1_3, ks1_4, kp_5, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];

        t_14[k] = f_7 * hd0_2[k]
                  - f_8 * hd1_4[k]
                  + pa_z[k] * id_9[k];

        t_15[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, hd0_3, hd0_10, hd0_11, hd1_5, hd1_16, \
                         hd1_17, id_12, id_17, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_9 * hd0_10[k]
                  - f_10 * hd1_16[k]
                  + pa_x[k] * id_17[k];

        t_17[k] = f_9 * hd0_3[k]
                  - f_10 * hd1_5[k]
                  + pa_y[k] * id_12[k];

        t_18[k] = f_7 * hd0_11[k]
                  - f_8 * hd1_17[k]
                  + pa_x[k] * id_19[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_y, pb_z, hd0_5, hd1_8, id_15, ks0_5, \
                         ks0_6, ks1_5, ks1_6, kp_7, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];

        t_20[k] = f_9 * hd0_5[k]
                  - f_10 * hd1_8[k]
                  + pa_z[k] * id_15[k];

        t_21[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pa_y, hd0_7, hd0_12, hd0_13, hd1_11, hd1_18, \
                         hd1_20, id_18, id_23, id_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_7 * hd0_12[k]
                  - f_8 * hd1_18[k]
                  + pa_x[k] * id_23[k];

        t_23[k] = f_5 * hd0_7[k]
                  - f_6 * hd1_11[k]
                  + pa_y[k] * id_18[k];

        t_24[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_20[k]
                  + pa_x[k] * id_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_z, pb_y, pb_z, hd0_9, hd1_14, id_21, ks0_7, \
                         ks0_8, ks1_7, ks1_8, kp_9, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];

        t_26[k] = f_5 * hd0_9[k]
                  - f_6 * hd1_14[k]
                  + pa_z[k] * id_21[k];

        t_27[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_x, hd0_20, hd1_32, id_25, id_27, \
                         id_44, ks0_11, ks1_11, kp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_32[k]
                  + pa_x[k] * id_25[k];

        t_29[k] = pa_x[k] * id_27[k];

        t_30[k] = pa_x[k] * id_44[k];

        t_31[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_11[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_z, pb_y, pb_z, ip_10, id_27, ks0_11, ks1_11, \
                         kp_12, kp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * ip_10[k]
                  + f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_y[k] * kp_12[k];

        t_33[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_z[k] * kp_13[k];

        t_34[k] = pa_z[k] * id_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pa_z, pb_x, hd0_13, hd0_16, hd1_20, hd1_25, \
                         id_29, id_33, ks0_12, ks1_12, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_14[k];

        t_36[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_20[k]
                  + pa_z[k] * id_29[k];

        t_37[k] = f_5 * hd0_16[k]
                  - f_6 * hd1_25[k]
                  + pa_y[k] * id_33[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pa_z, pb_x, hd0_14, hd0_18, hd1_22, hd1_28, \
                         id_32, id_36, ks0_13, ks1_13, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_15[k];

        t_39[k] = f_7 * hd0_14[k]
                  - f_8 * hd1_22[k]
                  + pa_z[k] * id_32[k];

        t_40[k] = f_9 * hd0_18[k]
                  - f_10 * hd1_28[k]
                  + pa_y[k] * id_36[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_y, pa_z, pb_x, hd0_15, hd0_19, hd1_24, hd1_29, \
                         id_35, id_39, ks0_14, ks1_14, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_16[k];

        t_42[k] = f_9 * hd0_15[k]
                  - f_10 * hd1_24[k]
                  + pa_z[k] * id_35[k];

        t_43[k] = f_7 * hd0_19[k]
                  - f_8 * hd1_29[k]
                  + pa_y[k] * id_39[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, pa_z, pb_x, hd0_17, hd0_20, hd1_27, hd1_32, \
                         id_38, id_41, ks0_15, ks1_15, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_1 * ks0_15[k]
                  - f_2 * ks1_15[k]
                  + pb_x[k] * kp_17[k];

        t_45[k] = f_5 * hd0_17[k]
                  - f_6 * hd1_27[k]
                  + pa_z[k] * id_38[k];

        t_46[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_32[k]
                  + pa_y[k] * id_41[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pb_x, pb_y, pb_z, ip_17, id_44, ks0_17, \
                         ks1_17, kp_18, kp_19, kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pa_y[k] * id_44[k];

        t_48[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_x[k] * kp_18[k];

        t_49[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_y[k] * kp_19[k];

        t_50[k] = f_0 * ip_17[k]
                  + f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hd0, const size_t hd1,
                                      const size_t ip, const size_t id, const size_t ks0,
                                      const size_t ks1, const size_t kp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 1.5 / alpha;
    const auto f_11 = 1.5 * beta / (alpha * p);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_27 = buffer.data(hd0 + 27);
    const auto *hd0_28 = buffer.data(hd0 + 28);
    const auto *hd0_29 = buffer.data(hd0 + 29);
    const auto *hd0_32 = buffer.data(hd0 + 32);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_28 = buffer.data(hd1 + 28);
    const auto *hd1_29 = buffer.data(hd1 + 29);
    const auto *hd1_32 = buffer.data(hd1 + 32);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);
    const auto *ks0_15 = buffer.data(ks0 + 15);
    const auto *ks0_17 = buffer.data(ks0 + 17);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);
    const auto *ks1_15 = buffer.data(ks1 + 15);
    const auto *ks1_17 = buffer.data(ks1 + 17);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, ip_0, id_0, ks0_0, ks1_0, \
                         kp_0, kp_1, kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];

        t_3[k] = pa_y[k] * id_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, hd0_0, hd0_6, hd1_0, hd1_6, \
                         ip_2, id_0, id_2, id_3, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * id_0[k];

        t_5[k] = f_3 * ip_2[k]
                 + pa_z[k] * id_2[k];

        t_6[k] = f_4 * hd0_0[k]
                 - f_5 * hd1_0[k]
                 + pa_y[k] * id_3[k];

        t_7[k] = f_6 * hd0_6[k]
                 - f_7 * hd1_6[k]
                 + pa_x[k] * id_6[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_z, pb_y, pb_z, hd0_0, hd1_0, id_4, ks0_1, ks0_2, \
                         ks1_1, ks1_2, kp_3, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_1[k]
                 + pb_z[k] * kp_3[k];

        t_9[k] = f_4 * hd0_0[k]
                 - f_5 * hd1_0[k]
                 + pa_z[k] * id_4[k];

        t_10[k] = f_1 * ks0_2[k]
                  - f_2 * ks1_2[k]
                  + pb_y[k] * kp_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, hd0_3, hd0_10, hd0_12, hd1_3, hd1_10, \
                         hd1_12, id_5, id_10, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_6 * hd0_10[k]
                  - f_7 * hd1_10[k]
                  + pa_x[k] * id_10[k];

        t_12[k] = f_8 * hd0_3[k]
                  - f_9 * hd1_3[k]
                  + pa_y[k] * id_5[k];

        t_13[k] = f_10 * hd0_12[k]
                  - f_11 * hd1_12[k]
                  + pa_x[k] * id_12[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_z, pb_y, pb_z, hd0_4, hd1_4, id_8, ks0_3, ks0_4, \
                         ks1_3, ks1_4, kp_5, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];

        t_15[k] = f_8 * hd0_4[k]
                  - f_9 * hd1_4[k]
                  + pa_z[k] * id_8[k];

        t_16[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, hd0_5, hd0_16, hd0_17, hd1_5, hd1_16, \
                         hd1_17, id_11, id_16, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_10 * hd0_16[k]
                  - f_11 * hd1_16[k]
                  + pa_x[k] * id_16[k];

        t_18[k] = f_10 * hd0_5[k]
                  - f_11 * hd1_5[k]
                  + pa_y[k] * id_11[k];

        t_19[k] = f_8 * hd0_17[k]
                  - f_9 * hd1_17[k]
                  + pa_x[k] * id_18[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, pb_z, hd0_8, hd1_8, id_14, ks0_5, \
                         ks0_6, ks1_5, ks1_6, kp_7, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];

        t_21[k] = f_10 * hd0_8[k]
                  - f_11 * hd1_8[k]
                  + pa_z[k] * id_14[k];

        t_22[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pa_y, hd0_11, hd0_18, hd0_20, hd1_11, hd1_18, \
                         hd1_20, id_17, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_8 * hd0_18[k]
                  - f_9 * hd1_18[k]
                  + pa_x[k] * id_22[k];

        t_24[k] = f_6 * hd0_11[k]
                  - f_7 * hd1_11[k]
                  + pa_y[k] * id_17[k];

        t_25[k] = f_4 * hd0_20[k]
                  - f_5 * hd1_20[k]
                  + pa_x[k] * id_23[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, pb_z, hd0_14, hd1_14, id_20, ks0_7, \
                         ks0_8, ks1_7, ks1_8, kp_9, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];

        t_27[k] = f_6 * hd0_14[k]
                  - f_7 * hd1_14[k]
                  + pa_z[k] * id_20[k];

        t_28[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_x, pb_x, hd0_32, hd1_32, id_24, id_26, \
                         id_41, ks0_11, ks1_11, kp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_4 * hd0_32[k]
                  - f_5 * hd1_32[k]
                  + pa_x[k] * id_24[k];

        t_30[k] = pa_x[k] * id_26[k];

        t_31[k] = pa_x[k] * id_41[k];

        t_32[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_11[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, ip_10, ip_11, id_26, id_27, \
                         ks0_11, ks1_11, kp_12, kp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * ip_10[k]
                  + f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_y[k] * kp_12[k];

        t_34[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_z[k] * kp_13[k];

        t_35[k] = pa_z[k] * id_26[k];

        t_36[k] = f_3 * ip_11[k]
                  + pa_z[k] * id_27[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, hd0_20, hd0_25, hd1_20, hd1_25, \
                         id_28, id_31, ks0_12, ks1_12, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_14[k];

        t_38[k] = f_4 * hd0_20[k]
                  - f_5 * hd1_20[k]
                  + pa_z[k] * id_28[k];

        t_39[k] = f_6 * hd0_25[k]
                  - f_7 * hd1_25[k]
                  + pa_y[k] * id_31[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_x, hd0_22, hd0_28, hd1_22, hd1_28, \
                         id_30, id_34, ks0_13, ks1_13, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_15[k];

        t_41[k] = f_8 * hd0_22[k]
                  - f_9 * hd1_22[k]
                  + pa_z[k] * id_30[k];

        t_42[k] = f_10 * hd0_28[k]
                  - f_11 * hd1_28[k]
                  + pa_y[k] * id_34[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_y, pa_z, pb_x, hd0_24, hd0_29, hd1_24, hd1_29, \
                         id_33, id_37, ks0_14, ks1_14, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_16[k];

        t_44[k] = f_10 * hd0_24[k]
                  - f_11 * hd1_24[k]
                  + pa_z[k] * id_33[k];

        t_45[k] = f_8 * hd0_29[k]
                  - f_9 * hd1_29[k]
                  + pa_y[k] * id_37[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_y, pa_z, pb_x, hd0_27, hd0_32, hd1_27, hd1_32, \
                         id_36, id_38, ks0_15, ks1_15, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_1 * ks0_15[k]
                  - f_2 * ks1_15[k]
                  + pb_x[k] * kp_17[k];

        t_47[k] = f_6 * hd0_27[k]
                  - f_7 * hd1_27[k]
                  + pa_z[k] * id_36[k];

        t_48[k] = f_4 * hd0_32[k]
                  - f_5 * hd1_32[k]
                  + pa_y[k] * id_38[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_y, pb_x, pb_y, ip_16, id_40, id_41, \
                         ks0_17, ks1_17, kp_18, kp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * ip_16[k]
                  + pa_y[k] * id_40[k];

        t_50[k] = pa_y[k] * id_41[k];

        t_51[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_x[k] * kp_18[k];

        t_52[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_y[k] * kp_19[k];
    }

#pragma omp simd aligned(t_53, pb_z, ip_17, ks0_17, ks1_17, kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * ip_17[k]
                  + f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_18(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hd0, const size_t hd1,
                                      const size_t ip, const size_t id, const size_t ks0,
                                      const size_t ks1, const size_t kp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / alpha;
    const auto f_6 = 2.0 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 1.5 / alpha;
    const auto f_10 = 1.5 * beta / (alpha * p);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_27 = buffer.data(hd0 + 27);
    const auto *hd0_28 = buffer.data(hd0 + 28);
    const auto *hd0_29 = buffer.data(hd0 + 29);
    const auto *hd0_32 = buffer.data(hd0 + 32);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_23 = buffer.data(hd1 + 23);
    const auto *hd1_26 = buffer.data(hd1 + 26);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_17 = buffer.data(ip + 17);

    const auto *id_0 = buffer.data(id + 0);
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
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_32 = buffer.data(id + 32);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);
    const auto *ks0_15 = buffer.data(ks0 + 15);
    const auto *ks0_17 = buffer.data(ks0 + 17);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);
    const auto *ks1_15 = buffer.data(ks1 + 15);
    const auto *ks1_17 = buffer.data(ks1 + 17);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, ip_0, id_0, ks0_0, ks1_0, \
                         kp_0, kp_1, kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];

        t_3[k] = pa_y[k] * id_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, hd0_0, hd0_6, hd1_0, hd1_6, id_0, \
                         id_3, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * id_0[k];

        t_5[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_y[k] * id_3[k];

        t_6[k] = f_5 * hd0_6[k]
                 - f_6 * hd1_6[k]
                 + pa_x[k] * id_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, hd0_0, hd1_0, id_4, ks0_1, ks0_2, \
                         ks1_1, ks1_2, kp_3, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_1[k]
                 + pb_z[k] * kp_3[k];

        t_8[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_z[k] * id_4[k];

        t_9[k] = f_1 * ks0_2[k]
                 - f_2 * ks1_2[k]
                 + pb_y[k] * kp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, hd0_3, hd0_10, hd0_12, hd1_3, hd1_8, \
                         hd1_10, id_5, id_8, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd0_10[k]
                  - f_6 * hd1_8[k]
                  + pa_x[k] * id_8[k];

        t_11[k] = f_7 * hd0_3[k]
                  - f_8 * hd1_3[k]
                  + pa_y[k] * id_5[k];

        t_12[k] = f_9 * hd0_12[k]
                  - f_10 * hd1_10[k]
                  + pa_x[k] * id_10[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, hd0_4, hd1_4, id_7, ks0_3, ks0_4, \
                         ks1_3, ks1_4, kp_5, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];

        t_14[k] = f_7 * hd0_4[k]
                  - f_8 * hd1_4[k]
                  + pa_z[k] * id_7[k];

        t_15[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, hd0_5, hd0_16, hd0_17, hd1_5, hd1_12, \
                         hd1_13, id_9, id_12, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_9 * hd0_16[k]
                  - f_10 * hd1_12[k]
                  + pa_x[k] * id_12[k];

        t_17[k] = f_9 * hd0_5[k]
                  - f_10 * hd1_5[k]
                  + pa_y[k] * id_9[k];

        t_18[k] = f_7 * hd0_17[k]
                  - f_8 * hd1_13[k]
                  + pa_x[k] * id_14[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_y, pb_z, hd0_8, hd1_7, id_11, ks0_5, \
                         ks0_6, ks1_5, ks1_6, kp_7, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];

        t_20[k] = f_9 * hd0_8[k]
                  - f_10 * hd1_7[k]
                  + pa_z[k] * id_11[k];

        t_21[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pa_y, hd0_11, hd0_18, hd0_20, hd1_9, hd1_14, \
                         hd1_16, id_13, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_7 * hd0_18[k]
                  - f_8 * hd1_14[k]
                  + pa_x[k] * id_16[k];

        t_23[k] = f_5 * hd0_11[k]
                  - f_6 * hd1_9[k]
                  + pa_y[k] * id_13[k];

        t_24[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_16[k]
                  + pa_x[k] * id_17[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_z, pb_y, pb_z, hd0_14, hd1_11, id_15, ks0_7, \
                         ks0_8, ks1_7, ks1_8, kp_9, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];

        t_26[k] = f_5 * hd0_14[k]
                  - f_6 * hd1_11[k]
                  + pa_z[k] * id_15[k];

        t_27[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_x, hd0_32, hd1_26, id_18, id_20, \
                         id_32, ks0_11, ks1_11, kp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * hd0_32[k]
                  - f_4 * hd1_26[k]
                  + pa_x[k] * id_18[k];

        t_29[k] = pa_x[k] * id_20[k];

        t_30[k] = pa_x[k] * id_32[k];

        t_31[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_11[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_z, pb_y, pb_z, ip_10, id_20, ks0_11, ks1_11, \
                         kp_12, kp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * ip_10[k]
                  + f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_y[k] * kp_12[k];

        t_33[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_z[k] * kp_13[k];

        t_34[k] = pa_z[k] * id_20[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pa_z, pb_x, hd0_20, hd0_25, hd1_16, hd1_20, \
                         id_22, id_24, ks0_12, ks1_12, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_14[k];

        t_36[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_16[k]
                  + pa_z[k] * id_22[k];

        t_37[k] = f_5 * hd0_25[k]
                  - f_6 * hd1_20[k]
                  + pa_y[k] * id_24[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pa_z, pb_x, hd0_22, hd0_28, hd1_18, hd1_22, \
                         id_23, id_26, ks0_13, ks1_13, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_15[k];

        t_39[k] = f_7 * hd0_22[k]
                  - f_8 * hd1_18[k]
                  + pa_z[k] * id_23[k];

        t_40[k] = f_9 * hd0_28[k]
                  - f_10 * hd1_22[k]
                  + pa_y[k] * id_26[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_y, pa_z, pb_x, hd0_24, hd0_29, hd1_19, hd1_23, \
                         id_25, id_28, ks0_14, ks1_14, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_16[k];

        t_42[k] = f_9 * hd0_24[k]
                  - f_10 * hd1_19[k]
                  + pa_z[k] * id_25[k];

        t_43[k] = f_7 * hd0_29[k]
                  - f_8 * hd1_23[k]
                  + pa_y[k] * id_28[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, pa_z, pb_x, hd0_27, hd0_32, hd1_21, hd1_26, \
                         id_27, id_29, ks0_15, ks1_15, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_1 * ks0_15[k]
                  - f_2 * ks1_15[k]
                  + pb_x[k] * kp_17[k];

        t_45[k] = f_5 * hd0_27[k]
                  - f_6 * hd1_21[k]
                  + pa_z[k] * id_27[k];

        t_46[k] = f_3 * hd0_32[k]
                  - f_4 * hd1_26[k]
                  + pa_y[k] * id_29[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pb_x, pb_y, pb_z, ip_17, id_32, ks0_17, \
                         ks1_17, kp_18, kp_19, kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pa_y[k] * id_32[k];

        t_48[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_x[k] * kp_18[k];

        t_49[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_y[k] * kp_19[k];

        t_50[k] = f_0 * ip_17[k]
                  + f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_19(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hd0, const size_t hd1,
                                      const size_t ip, const size_t id, const size_t ks0,
                                      const size_t ks1, const size_t kp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / alpha;
    const auto f_6 = 2.0 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 1.5 / alpha;
    const auto f_10 = 1.5 * beta / (alpha * p);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_1 = buffer.data(hd1 + 1);
    const auto *hd1_2 = buffer.data(hd1 + 2);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_7 = buffer.data(hd1 + 7);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_9 = buffer.data(hd1 + 9);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_13 = buffer.data(hd1 + 13);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_15 = buffer.data(hd1 + 15);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_19 = buffer.data(hd1 + 19);
    const auto *hd1_20 = buffer.data(hd1 + 20);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_17 = buffer.data(ip + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);
    const auto *ks0_15 = buffer.data(ks0 + 15);
    const auto *ks0_17 = buffer.data(ks0 + 17);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);
    const auto *ks1_15 = buffer.data(ks1 + 15);
    const auto *ks1_17 = buffer.data(ks1 + 17);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, ip_0, id_0, ks0_0, ks1_0, \
                         kp_0, kp_1, kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];

        t_3[k] = pa_y[k] * id_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, hd0_0, hd0_4, hd1_0, hd1_4, id_0, \
                         id_3, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * id_0[k];

        t_5[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_y[k] * id_3[k];

        t_6[k] = f_5 * hd0_4[k]
                 - f_6 * hd1_4[k]
                 + pa_x[k] * id_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, hd0_0, hd1_0, id_4, ks0_1, ks0_2, \
                         ks1_1, ks1_2, kp_3, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_1[k]
                 + pb_z[k] * kp_3[k];

        t_8[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_z[k] * id_4[k];

        t_9[k] = f_1 * ks0_2[k]
                 - f_2 * ks1_2[k]
                 + pb_y[k] * kp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, hd0_1, hd0_6, hd0_8, hd1_1, hd1_6, \
                         hd1_8, id_5, id_10, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd0_6[k]
                  - f_6 * hd1_6[k]
                  + pa_x[k] * id_10[k];

        t_11[k] = f_7 * hd0_1[k]
                  - f_8 * hd1_1[k]
                  + pa_y[k] * id_5[k];

        t_12[k] = f_9 * hd0_8[k]
                  - f_10 * hd1_8[k]
                  + pa_x[k] * id_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, hd0_2, hd1_2, id_8, ks0_3, ks0_4, \
                         ks1_3, ks1_4, kp_5, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];

        t_14[k] = f_7 * hd0_2[k]
                  - f_8 * hd1_2[k]
                  + pa_z[k] * id_8[k];

        t_15[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, hd0_3, hd0_10, hd0_11, hd1_3, hd1_10, \
                         hd1_11, id_11, id_16, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_9 * hd0_10[k]
                  - f_10 * hd1_10[k]
                  + pa_x[k] * id_16[k];

        t_17[k] = f_9 * hd0_3[k]
                  - f_10 * hd1_3[k]
                  + pa_y[k] * id_11[k];

        t_18[k] = f_7 * hd0_11[k]
                  - f_8 * hd1_11[k]
                  + pa_x[k] * id_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_y, pb_z, hd0_5, hd1_5, id_14, ks0_5, \
                         ks0_6, ks1_5, ks1_6, kp_7, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];

        t_20[k] = f_9 * hd0_5[k]
                  - f_10 * hd1_5[k]
                  + pa_z[k] * id_14[k];

        t_21[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pa_y, hd0_7, hd0_12, hd0_13, hd1_7, hd1_12, \
                         hd1_13, id_17, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_7 * hd0_12[k]
                  - f_8 * hd1_12[k]
                  + pa_x[k] * id_22[k];

        t_23[k] = f_5 * hd0_7[k]
                  - f_6 * hd1_7[k]
                  + pa_y[k] * id_17[k];

        t_24[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_13[k]
                  + pa_x[k] * id_23[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_z, pb_y, pb_z, hd0_9, hd1_9, id_20, ks0_7, \
                         ks0_8, ks1_7, ks1_8, kp_9, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];

        t_26[k] = f_5 * hd0_9[k]
                  - f_6 * hd1_9[k]
                  + pa_z[k] * id_20[k];

        t_27[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_x, hd0_20, hd1_20, id_24, id_26, \
                         id_41, ks0_11, ks1_11, kp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pa_x[k] * id_24[k];

        t_29[k] = pa_x[k] * id_26[k];

        t_30[k] = pa_x[k] * id_41[k];

        t_31[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_11[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_z, pb_y, pb_z, ip_10, id_26, ks0_11, ks1_11, \
                         kp_12, kp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * ip_10[k]
                  + f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_y[k] * kp_12[k];

        t_33[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_z[k] * kp_13[k];

        t_34[k] = pa_z[k] * id_26[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pa_z, pb_x, hd0_13, hd0_16, hd1_13, hd1_16, \
                         id_28, id_31, ks0_12, ks1_12, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_14[k];

        t_36[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_13[k]
                  + pa_z[k] * id_28[k];

        t_37[k] = f_5 * hd0_16[k]
                  - f_6 * hd1_16[k]
                  + pa_y[k] * id_31[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pa_z, pb_x, hd0_14, hd0_18, hd1_14, hd1_18, \
                         id_30, id_34, ks0_13, ks1_13, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_15[k];

        t_39[k] = f_7 * hd0_14[k]
                  - f_8 * hd1_14[k]
                  + pa_z[k] * id_30[k];

        t_40[k] = f_9 * hd0_18[k]
                  - f_10 * hd1_18[k]
                  + pa_y[k] * id_34[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_y, pa_z, pb_x, hd0_15, hd0_19, hd1_15, hd1_19, \
                         id_33, id_37, ks0_14, ks1_14, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_16[k];

        t_42[k] = f_9 * hd0_15[k]
                  - f_10 * hd1_15[k]
                  + pa_z[k] * id_33[k];

        t_43[k] = f_7 * hd0_19[k]
                  - f_8 * hd1_19[k]
                  + pa_y[k] * id_37[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, pa_z, pb_x, hd0_17, hd0_20, hd1_17, hd1_20, \
                         id_36, id_38, ks0_15, ks1_15, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_1 * ks0_15[k]
                  - f_2 * ks1_15[k]
                  + pb_x[k] * kp_17[k];

        t_45[k] = f_5 * hd0_17[k]
                  - f_6 * hd1_17[k]
                  + pa_z[k] * id_36[k];

        t_46[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pa_y[k] * id_38[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pb_x, pb_y, pb_z, ip_17, id_41, ks0_17, \
                         ks1_17, kp_18, kp_19, kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pa_y[k] * id_41[k];

        t_48[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_x[k] * kp_18[k];

        t_49[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_y[k] * kp_19[k];

        t_50[k] = f_0 * ip_17[k]
                  + f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_20(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hd0, const size_t hd1,
                                      const size_t ip, const size_t id, const size_t ks0,
                                      const size_t ks1, const size_t kp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / alpha;
    const auto f_6 = 2.0 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 1.5 / alpha;
    const auto f_10 = 1.5 * beta / (alpha * p);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_1 = buffer.data(hd0 + 1);
    const auto *hd0_2 = buffer.data(hd0 + 2);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_7 = buffer.data(hd0 + 7);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_9 = buffer.data(hd0 + 9);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_13 = buffer.data(hd0 + 13);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_15 = buffer.data(hd0 + 15);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_19 = buffer.data(hd0 + 19);
    const auto *hd0_20 = buffer.data(hd0 + 20);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_28 = buffer.data(hd1 + 28);
    const auto *hd1_29 = buffer.data(hd1 + 29);
    const auto *hd1_32 = buffer.data(hd1 + 32);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_17 = buffer.data(ip + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);
    const auto *ks0_15 = buffer.data(ks0 + 15);
    const auto *ks0_17 = buffer.data(ks0 + 17);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);
    const auto *ks1_15 = buffer.data(ks1 + 15);
    const auto *ks1_17 = buffer.data(ks1 + 17);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, ip_0, id_0, ks0_0, ks1_0, \
                         kp_0, kp_1, kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];

        t_3[k] = pa_y[k] * id_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, hd0_0, hd0_4, hd1_0, hd1_6, id_0, \
                         id_3, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * id_0[k];

        t_5[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_y[k] * id_3[k];

        t_6[k] = f_5 * hd0_4[k]
                 - f_6 * hd1_6[k]
                 + pa_x[k] * id_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, hd0_0, hd1_0, id_4, ks0_1, ks0_2, \
                         ks1_1, ks1_2, kp_3, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_1[k]
                 + pb_z[k] * kp_3[k];

        t_8[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_z[k] * id_4[k];

        t_9[k] = f_1 * ks0_2[k]
                 - f_2 * ks1_2[k]
                 + pb_y[k] * kp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, hd0_1, hd0_6, hd0_8, hd1_3, hd1_10, \
                         hd1_12, id_5, id_10, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd0_6[k]
                  - f_6 * hd1_10[k]
                  + pa_x[k] * id_10[k];

        t_11[k] = f_7 * hd0_1[k]
                  - f_8 * hd1_3[k]
                  + pa_y[k] * id_5[k];

        t_12[k] = f_9 * hd0_8[k]
                  - f_10 * hd1_12[k]
                  + pa_x[k] * id_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, hd0_2, hd1_4, id_8, ks0_3, ks0_4, \
                         ks1_3, ks1_4, kp_5, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];

        t_14[k] = f_7 * hd0_2[k]
                  - f_8 * hd1_4[k]
                  + pa_z[k] * id_8[k];

        t_15[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, hd0_3, hd0_10, hd0_11, hd1_5, hd1_16, \
                         hd1_17, id_11, id_16, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_9 * hd0_10[k]
                  - f_10 * hd1_16[k]
                  + pa_x[k] * id_16[k];

        t_17[k] = f_9 * hd0_3[k]
                  - f_10 * hd1_5[k]
                  + pa_y[k] * id_11[k];

        t_18[k] = f_7 * hd0_11[k]
                  - f_8 * hd1_17[k]
                  + pa_x[k] * id_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_y, pb_z, hd0_5, hd1_8, id_14, ks0_5, \
                         ks0_6, ks1_5, ks1_6, kp_7, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];

        t_20[k] = f_9 * hd0_5[k]
                  - f_10 * hd1_8[k]
                  + pa_z[k] * id_14[k];

        t_21[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pa_y, hd0_7, hd0_12, hd0_13, hd1_11, hd1_18, \
                         hd1_20, id_17, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_7 * hd0_12[k]
                  - f_8 * hd1_18[k]
                  + pa_x[k] * id_22[k];

        t_23[k] = f_5 * hd0_7[k]
                  - f_6 * hd1_11[k]
                  + pa_y[k] * id_17[k];

        t_24[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_20[k]
                  + pa_x[k] * id_23[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_z, pb_y, pb_z, hd0_9, hd1_14, id_20, ks0_7, \
                         ks0_8, ks1_7, ks1_8, kp_9, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];

        t_26[k] = f_5 * hd0_9[k]
                  - f_6 * hd1_14[k]
                  + pa_z[k] * id_20[k];

        t_27[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_x, hd0_20, hd1_32, id_24, id_26, \
                         id_41, ks0_11, ks1_11, kp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_32[k]
                  + pa_x[k] * id_24[k];

        t_29[k] = pa_x[k] * id_26[k];

        t_30[k] = pa_x[k] * id_41[k];

        t_31[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_11[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_z, pb_y, pb_z, ip_10, id_26, ks0_11, ks1_11, \
                         kp_12, kp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * ip_10[k]
                  + f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_y[k] * kp_12[k];

        t_33[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_z[k] * kp_13[k];

        t_34[k] = pa_z[k] * id_26[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pa_z, pb_x, hd0_13, hd0_16, hd1_20, hd1_25, \
                         id_28, id_31, ks0_12, ks1_12, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_14[k];

        t_36[k] = f_3 * hd0_13[k]
                  - f_4 * hd1_20[k]
                  + pa_z[k] * id_28[k];

        t_37[k] = f_5 * hd0_16[k]
                  - f_6 * hd1_25[k]
                  + pa_y[k] * id_31[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pa_z, pb_x, hd0_14, hd0_18, hd1_22, hd1_28, \
                         id_30, id_34, ks0_13, ks1_13, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_15[k];

        t_39[k] = f_7 * hd0_14[k]
                  - f_8 * hd1_22[k]
                  + pa_z[k] * id_30[k];

        t_40[k] = f_9 * hd0_18[k]
                  - f_10 * hd1_28[k]
                  + pa_y[k] * id_34[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_y, pa_z, pb_x, hd0_15, hd0_19, hd1_24, hd1_29, \
                         id_33, id_37, ks0_14, ks1_14, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_16[k];

        t_42[k] = f_9 * hd0_15[k]
                  - f_10 * hd1_24[k]
                  + pa_z[k] * id_33[k];

        t_43[k] = f_7 * hd0_19[k]
                  - f_8 * hd1_29[k]
                  + pa_y[k] * id_37[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, pa_z, pb_x, hd0_17, hd0_20, hd1_27, hd1_32, \
                         id_36, id_38, ks0_15, ks1_15, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_1 * ks0_15[k]
                  - f_2 * ks1_15[k]
                  + pb_x[k] * kp_17[k];

        t_45[k] = f_5 * hd0_17[k]
                  - f_6 * hd1_27[k]
                  + pa_z[k] * id_36[k];

        t_46[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_32[k]
                  + pa_y[k] * id_38[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pb_x, pb_y, pb_z, ip_17, id_41, ks0_17, \
                         ks1_17, kp_18, kp_19, kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pa_y[k] * id_41[k];

        t_48[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_x[k] * kp_18[k];

        t_49[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_y[k] * kp_19[k];

        t_50[k] = f_0 * ip_17[k]
                  + f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_z[k] * kp_20[k];
    }
}

auto
compute_prim_kd_electron_repulsion_21(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hd0, const size_t hd1,
                                      const size_t ip, const size_t id, const size_t ks0,
                                      const size_t ks1, const size_t kp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.0 / alpha;
    const auto f_6 = 2.0 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 1.5 / alpha;
    const auto f_10 = 1.5 * beta / (alpha * p);

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

    const auto *hd0_0 = buffer.data(hd0 + 0);
    const auto *hd0_3 = buffer.data(hd0 + 3);
    const auto *hd0_4 = buffer.data(hd0 + 4);
    const auto *hd0_5 = buffer.data(hd0 + 5);
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_8 = buffer.data(hd0 + 8);
    const auto *hd0_10 = buffer.data(hd0 + 10);
    const auto *hd0_11 = buffer.data(hd0 + 11);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_14 = buffer.data(hd0 + 14);
    const auto *hd0_16 = buffer.data(hd0 + 16);
    const auto *hd0_17 = buffer.data(hd0 + 17);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_20 = buffer.data(hd0 + 20);
    const auto *hd0_22 = buffer.data(hd0 + 22);
    const auto *hd0_24 = buffer.data(hd0 + 24);
    const auto *hd0_25 = buffer.data(hd0 + 25);
    const auto *hd0_27 = buffer.data(hd0 + 27);
    const auto *hd0_28 = buffer.data(hd0 + 28);
    const auto *hd0_29 = buffer.data(hd0 + 29);
    const auto *hd0_32 = buffer.data(hd0 + 32);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_3 = buffer.data(hd1 + 3);
    const auto *hd1_4 = buffer.data(hd1 + 4);
    const auto *hd1_5 = buffer.data(hd1 + 5);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_8 = buffer.data(hd1 + 8);
    const auto *hd1_10 = buffer.data(hd1 + 10);
    const auto *hd1_11 = buffer.data(hd1 + 11);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_14 = buffer.data(hd1 + 14);
    const auto *hd1_16 = buffer.data(hd1 + 16);
    const auto *hd1_17 = buffer.data(hd1 + 17);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_20 = buffer.data(hd1 + 20);
    const auto *hd1_22 = buffer.data(hd1 + 22);
    const auto *hd1_24 = buffer.data(hd1 + 24);
    const auto *hd1_25 = buffer.data(hd1 + 25);
    const auto *hd1_27 = buffer.data(hd1 + 27);
    const auto *hd1_28 = buffer.data(hd1 + 28);
    const auto *hd1_29 = buffer.data(hd1 + 29);
    const auto *hd1_32 = buffer.data(hd1 + 32);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_17 = buffer.data(ip + 17);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_1 = buffer.data(ks0 + 1);
    const auto *ks0_2 = buffer.data(ks0 + 2);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_4 = buffer.data(ks0 + 4);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_7 = buffer.data(ks0 + 7);
    const auto *ks0_8 = buffer.data(ks0 + 8);
    const auto *ks0_11 = buffer.data(ks0 + 11);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_13 = buffer.data(ks0 + 13);
    const auto *ks0_14 = buffer.data(ks0 + 14);
    const auto *ks0_15 = buffer.data(ks0 + 15);
    const auto *ks0_17 = buffer.data(ks0 + 17);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_1 = buffer.data(ks1 + 1);
    const auto *ks1_2 = buffer.data(ks1 + 2);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_4 = buffer.data(ks1 + 4);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_7 = buffer.data(ks1 + 7);
    const auto *ks1_8 = buffer.data(ks1 + 8);
    const auto *ks1_11 = buffer.data(ks1 + 11);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_13 = buffer.data(ks1 + 13);
    const auto *ks1_14 = buffer.data(ks1 + 14);
    const auto *ks1_15 = buffer.data(ks1 + 15);
    const auto *ks1_17 = buffer.data(ks1 + 17);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_5 = buffer.data(kp + 5);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_7 = buffer.data(kp + 7);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_12 = buffer.data(kp + 12);
    const auto *kp_13 = buffer.data(kp + 13);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, ip_0, id_0, ks0_0, ks1_0, \
                         kp_0, kp_1, kp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_0[k]
                 + f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_x[k] * kp_0[k];

        t_1[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_y[k] * kp_1[k];

        t_2[k] = f_1 * ks0_0[k]
                 - f_2 * ks1_0[k]
                 + pb_z[k] * kp_2[k];

        t_3[k] = pa_y[k] * id_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, hd0_0, hd0_6, hd1_0, hd1_6, id_0, \
                         id_3, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * id_0[k];

        t_5[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_y[k] * id_3[k];

        t_6[k] = f_5 * hd0_6[k]
                 - f_6 * hd1_6[k]
                 + pa_x[k] * id_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, hd0_0, hd1_0, id_4, ks0_1, ks0_2, \
                         ks1_1, ks1_2, kp_3, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * ks0_1[k]
                 - f_2 * ks1_1[k]
                 + pb_z[k] * kp_3[k];

        t_8[k] = f_3 * hd0_0[k]
                 - f_4 * hd1_0[k]
                 + pa_z[k] * id_4[k];

        t_9[k] = f_1 * ks0_2[k]
                 - f_2 * ks1_2[k]
                 + pb_y[k] * kp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, hd0_3, hd0_10, hd0_12, hd1_3, hd1_10, \
                         hd1_12, id_5, id_10, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * hd0_10[k]
                  - f_6 * hd1_10[k]
                  + pa_x[k] * id_10[k];

        t_11[k] = f_7 * hd0_3[k]
                  - f_8 * hd1_3[k]
                  + pa_y[k] * id_5[k];

        t_12[k] = f_9 * hd0_12[k]
                  - f_10 * hd1_12[k]
                  + pa_x[k] * id_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, hd0_4, hd1_4, id_8, ks0_3, ks0_4, \
                         ks1_3, ks1_4, kp_5, kp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_5[k];

        t_14[k] = f_7 * hd0_4[k]
                  - f_8 * hd1_4[k]
                  + pa_z[k] * id_8[k];

        t_15[k] = f_1 * ks0_4[k]
                  - f_2 * ks1_4[k]
                  + pb_y[k] * kp_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, hd0_5, hd0_16, hd0_17, hd1_5, hd1_16, \
                         hd1_17, id_11, id_16, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_9 * hd0_16[k]
                  - f_10 * hd1_16[k]
                  + pa_x[k] * id_16[k];

        t_17[k] = f_9 * hd0_5[k]
                  - f_10 * hd1_5[k]
                  + pa_y[k] * id_11[k];

        t_18[k] = f_7 * hd0_17[k]
                  - f_8 * hd1_17[k]
                  + pa_x[k] * id_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_y, pb_z, hd0_8, hd1_8, id_14, ks0_5, \
                         ks0_6, ks1_5, ks1_6, kp_7, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_z[k] * kp_7[k];

        t_20[k] = f_9 * hd0_8[k]
                  - f_10 * hd1_8[k]
                  + pa_z[k] * id_14[k];

        t_21[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_y[k] * kp_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pa_y, hd0_11, hd0_18, hd0_20, hd1_11, hd1_18, \
                         hd1_20, id_17, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_7 * hd0_18[k]
                  - f_8 * hd1_18[k]
                  + pa_x[k] * id_22[k];

        t_23[k] = f_5 * hd0_11[k]
                  - f_6 * hd1_11[k]
                  + pa_y[k] * id_17[k];

        t_24[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pa_x[k] * id_23[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_z, pb_y, pb_z, hd0_14, hd1_14, id_20, ks0_7, \
                         ks0_8, ks1_7, ks1_8, kp_9, kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * ks0_7[k]
                  - f_2 * ks1_7[k]
                  + pb_z[k] * kp_9[k];

        t_26[k] = f_5 * hd0_14[k]
                  - f_6 * hd1_14[k]
                  + pa_z[k] * id_20[k];

        t_27[k] = f_1 * ks0_8[k]
                  - f_2 * ks1_8[k]
                  + pb_y[k] * kp_10[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_x, hd0_32, hd1_32, id_24, id_26, \
                         id_41, ks0_11, ks1_11, kp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * hd0_32[k]
                  - f_4 * hd1_32[k]
                  + pa_x[k] * id_24[k];

        t_29[k] = pa_x[k] * id_26[k];

        t_30[k] = pa_x[k] * id_41[k];

        t_31[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_x[k] * kp_11[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_z, pb_y, pb_z, ip_10, id_26, ks0_11, ks1_11, \
                         kp_12, kp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * ip_10[k]
                  + f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_y[k] * kp_12[k];

        t_33[k] = f_1 * ks0_11[k]
                  - f_2 * ks1_11[k]
                  + pb_z[k] * kp_13[k];

        t_34[k] = pa_z[k] * id_26[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pa_z, pb_x, hd0_20, hd0_25, hd1_20, hd1_25, \
                         id_28, id_31, ks0_12, ks1_12, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * ks0_12[k]
                  - f_2 * ks1_12[k]
                  + pb_x[k] * kp_14[k];

        t_36[k] = f_3 * hd0_20[k]
                  - f_4 * hd1_20[k]
                  + pa_z[k] * id_28[k];

        t_37[k] = f_5 * hd0_25[k]
                  - f_6 * hd1_25[k]
                  + pa_y[k] * id_31[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pa_z, pb_x, hd0_22, hd0_28, hd1_22, hd1_28, \
                         id_30, id_34, ks0_13, ks1_13, kp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_1 * ks0_13[k]
                  - f_2 * ks1_13[k]
                  + pb_x[k] * kp_15[k];

        t_39[k] = f_7 * hd0_22[k]
                  - f_8 * hd1_22[k]
                  + pa_z[k] * id_30[k];

        t_40[k] = f_9 * hd0_28[k]
                  - f_10 * hd1_28[k]
                  + pa_y[k] * id_34[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_y, pa_z, pb_x, hd0_24, hd0_29, hd1_24, hd1_29, \
                         id_33, id_37, ks0_14, ks1_14, kp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_x[k] * kp_16[k];

        t_42[k] = f_9 * hd0_24[k]
                  - f_10 * hd1_24[k]
                  + pa_z[k] * id_33[k];

        t_43[k] = f_7 * hd0_29[k]
                  - f_8 * hd1_29[k]
                  + pa_y[k] * id_37[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, pa_z, pb_x, hd0_27, hd0_32, hd1_27, hd1_32, \
                         id_36, id_38, ks0_15, ks1_15, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_1 * ks0_15[k]
                  - f_2 * ks1_15[k]
                  + pb_x[k] * kp_17[k];

        t_45[k] = f_5 * hd0_27[k]
                  - f_6 * hd1_27[k]
                  + pa_z[k] * id_36[k];

        t_46[k] = f_3 * hd0_32[k]
                  - f_4 * hd1_32[k]
                  + pa_y[k] * id_38[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pb_x, pb_y, pb_z, ip_17, id_41, ks0_17, \
                         ks1_17, kp_18, kp_19, kp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pa_y[k] * id_41[k];

        t_48[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_x[k] * kp_18[k];

        t_49[k] = f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_y[k] * kp_19[k];

        t_50[k] = f_0 * ip_17[k]
                  + f_1 * ks0_17[k]
                  - f_2 * ks1_17[k]
                  + pb_z[k] * kp_20[k];
    }
}

}  // namespace simdt2ceri
