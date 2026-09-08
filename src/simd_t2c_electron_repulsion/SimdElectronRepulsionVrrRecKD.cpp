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
    const auto *hd0_6 = buffer.data(hd0 + 6);
    const auto *hd0_12 = buffer.data(hd0 + 12);
    const auto *hd0_18 = buffer.data(hd0 + 18);
    const auto *hd0_21 = buffer.data(hd0 + 21);
    const auto *hd0_30 = buffer.data(hd0 + 30);
    const auto *hd0_35 = buffer.data(hd0 + 35);
    const auto *hd0_36 = buffer.data(hd0 + 36);
    const auto *hd0_39 = buffer.data(hd0 + 39);
    const auto *hd0_48 = buffer.data(hd0 + 48);
    const auto *hd0_54 = buffer.data(hd0 + 54);
    const auto *hd0_59 = buffer.data(hd0 + 59);
    const auto *hd0_63 = buffer.data(hd0 + 63);
    const auto *hd0_75 = buffer.data(hd0 + 75);
    const auto *hd0_77 = buffer.data(hd0 + 77);
    const auto *hd0_89 = buffer.data(hd0 + 89);
    const auto *hd0_93 = buffer.data(hd0 + 93);
    const auto *hd0_99 = buffer.data(hd0 + 99);
    const auto *hd0_105 = buffer.data(hd0 + 105);
    const auto *hd0_107 = buffer.data(hd0 + 107);
    const auto *hd0_111 = buffer.data(hd0 + 111);
    const auto *hd0_113 = buffer.data(hd0 + 113);
    const auto *hd0_119 = buffer.data(hd0 + 119);
    const auto *hd0_125 = buffer.data(hd0 + 125);

    const auto *hd1_0 = buffer.data(hd1 + 0);
    const auto *hd1_6 = buffer.data(hd1 + 6);
    const auto *hd1_12 = buffer.data(hd1 + 12);
    const auto *hd1_18 = buffer.data(hd1 + 18);
    const auto *hd1_21 = buffer.data(hd1 + 21);
    const auto *hd1_30 = buffer.data(hd1 + 30);
    const auto *hd1_35 = buffer.data(hd1 + 35);
    const auto *hd1_36 = buffer.data(hd1 + 36);
    const auto *hd1_39 = buffer.data(hd1 + 39);
    const auto *hd1_48 = buffer.data(hd1 + 48);
    const auto *hd1_54 = buffer.data(hd1 + 54);
    const auto *hd1_59 = buffer.data(hd1 + 59);
    const auto *hd1_63 = buffer.data(hd1 + 63);
    const auto *hd1_75 = buffer.data(hd1 + 75);
    const auto *hd1_77 = buffer.data(hd1 + 77);
    const auto *hd1_89 = buffer.data(hd1 + 89);
    const auto *hd1_93 = buffer.data(hd1 + 93);
    const auto *hd1_99 = buffer.data(hd1 + 99);
    const auto *hd1_105 = buffer.data(hd1 + 105);
    const auto *hd1_107 = buffer.data(hd1 + 107);
    const auto *hd1_111 = buffer.data(hd1 + 111);
    const auto *hd1_113 = buffer.data(hd1 + 113);
    const auto *hd1_119 = buffer.data(hd1 + 119);
    const auto *hd1_125 = buffer.data(hd1 + 125);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_28 = buffer.data(ip + 28);
    const auto *ip_29 = buffer.data(ip + 29);
    const auto *ip_31 = buffer.data(ip + 31);
    const auto *ip_32 = buffer.data(ip + 32);
    const auto *ip_35 = buffer.data(ip + 35);
    const auto *ip_37 = buffer.data(ip + 37);
    const auto *ip_38 = buffer.data(ip + 38);
    const auto *ip_40 = buffer.data(ip + 40);
    const auto *ip_41 = buffer.data(ip + 41);
    const auto *ip_43 = buffer.data(ip + 43);
    const auto *ip_44 = buffer.data(ip + 44);
    const auto *ip_46 = buffer.data(ip + 46);
    const auto *ip_50 = buffer.data(ip + 50);
    const auto *ip_52 = buffer.data(ip + 52);
    const auto *ip_53 = buffer.data(ip + 53);
    const auto *ip_55 = buffer.data(ip + 55);
    const auto *ip_56 = buffer.data(ip + 56);
    const auto *ip_58 = buffer.data(ip + 58);
    const auto *ip_62 = buffer.data(ip + 62);
    const auto *ip_63 = buffer.data(ip + 63);
    const auto *ip_64 = buffer.data(ip + 64);
    const auto *ip_65 = buffer.data(ip + 65);
    const auto *ip_68 = buffer.data(ip + 68);
    const auto *ip_69 = buffer.data(ip + 69);
    const auto *ip_70 = buffer.data(ip + 70);
    const auto *ip_71 = buffer.data(ip + 71);
    const auto *ip_72 = buffer.data(ip + 72);
    const auto *ip_73 = buffer.data(ip + 73);
    const auto *ip_74 = buffer.data(ip + 74);
    const auto *ip_75 = buffer.data(ip + 75);
    const auto *ip_76 = buffer.data(ip + 76);
    const auto *ip_77 = buffer.data(ip + 77);
    const auto *ip_79 = buffer.data(ip + 79);
    const auto *ip_80 = buffer.data(ip + 80);
    const auto *ip_81 = buffer.data(ip + 81);
    const auto *ip_82 = buffer.data(ip + 82);
    const auto *ip_83 = buffer.data(ip + 83);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_72 = buffer.data(id + 72);
    const auto *id_75 = buffer.data(id + 75);
    const auto *id_77 = buffer.data(id + 77);
    const auto *id_78 = buffer.data(id + 78);
    const auto *id_84 = buffer.data(id + 84);
    const auto *id_86 = buffer.data(id + 86);
    const auto *id_87 = buffer.data(id + 87);
    const auto *id_89 = buffer.data(id + 89);
    const auto *id_90 = buffer.data(id + 90);
    const auto *id_91 = buffer.data(id + 91);
    const auto *id_93 = buffer.data(id + 93);
    const auto *id_105 = buffer.data(id + 105);
    const auto *id_107 = buffer.data(id + 107);
    const auto *id_111 = buffer.data(id + 111);
    const auto *id_113 = buffer.data(id + 113);
    const auto *id_120 = buffer.data(id + 120);
    const auto *id_122 = buffer.data(id + 122);
    const auto *id_125 = buffer.data(id + 125);
    const auto *id_126 = buffer.data(id + 126);
    const auto *id_129 = buffer.data(id + 129);
    const auto *id_131 = buffer.data(id + 131);
    const auto *id_135 = buffer.data(id + 135);
    const auto *id_136 = buffer.data(id + 136);
    const auto *id_137 = buffer.data(id + 137);
    const auto *id_138 = buffer.data(id + 138);
    const auto *id_141 = buffer.data(id + 141);
    const auto *id_142 = buffer.data(id + 142);
    const auto *id_143 = buffer.data(id + 143);
    const auto *id_144 = buffer.data(id + 144);
    const auto *id_147 = buffer.data(id + 147);
    const auto *id_148 = buffer.data(id + 148);
    const auto *id_149 = buffer.data(id + 149);
    const auto *id_150 = buffer.data(id + 150);
    const auto *id_153 = buffer.data(id + 153);
    const auto *id_154 = buffer.data(id + 154);
    const auto *id_155 = buffer.data(id + 155);
    const auto *id_159 = buffer.data(id + 159);
    const auto *id_160 = buffer.data(id + 160);
    const auto *id_161 = buffer.data(id + 161);
    const auto *id_162 = buffer.data(id + 162);
    const auto *id_165 = buffer.data(id + 165);
    const auto *id_167 = buffer.data(id + 167);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_9 = buffer.data(ks0 + 9);
    const auto *ks0_10 = buffer.data(ks0 + 10);
    const auto *ks0_14 = buffer.data(ks0 + 14);
    const auto *ks0_15 = buffer.data(ks0 + 15);
    const auto *ks0_20 = buffer.data(ks0 + 20);
    const auto *ks0_28 = buffer.data(ks0 + 28);
    const auto *ks0_30 = buffer.data(ks0 + 30);
    const auto *ks0_31 = buffer.data(ks0 + 31);
    const auto *ks0_32 = buffer.data(ks0 + 32);
    const auto *ks0_33 = buffer.data(ks0 + 33);
    const auto *ks0_35 = buffer.data(ks0 + 35);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_9 = buffer.data(ks1 + 9);
    const auto *ks1_10 = buffer.data(ks1 + 10);
    const auto *ks1_14 = buffer.data(ks1 + 14);
    const auto *ks1_15 = buffer.data(ks1 + 15);
    const auto *ks1_20 = buffer.data(ks1 + 20);
    const auto *ks1_28 = buffer.data(ks1 + 28);
    const auto *ks1_30 = buffer.data(ks1 + 30);
    const auto *ks1_31 = buffer.data(ks1 + 31);
    const auto *ks1_32 = buffer.data(ks1 + 32);
    const auto *ks1_33 = buffer.data(ks1 + 33);
    const auto *ks1_35 = buffer.data(ks1 + 35);

    const auto *kp_0 = buffer.data(kp + 0);
    const auto *kp_1 = buffer.data(kp + 1);
    const auto *kp_2 = buffer.data(kp + 2);
    const auto *kp_3 = buffer.data(kp + 3);
    const auto *kp_4 = buffer.data(kp + 4);
    const auto *kp_6 = buffer.data(kp + 6);
    const auto *kp_8 = buffer.data(kp + 8);
    const auto *kp_9 = buffer.data(kp + 9);
    const auto *kp_10 = buffer.data(kp + 10);
    const auto *kp_11 = buffer.data(kp + 11);
    const auto *kp_14 = buffer.data(kp + 14);
    const auto *kp_15 = buffer.data(kp + 15);
    const auto *kp_16 = buffer.data(kp + 16);
    const auto *kp_17 = buffer.data(kp + 17);
    const auto *kp_18 = buffer.data(kp + 18);
    const auto *kp_19 = buffer.data(kp + 19);
    const auto *kp_20 = buffer.data(kp + 20);
    const auto *kp_23 = buffer.data(kp + 23);
    const auto *kp_25 = buffer.data(kp + 25);
    const auto *kp_26 = buffer.data(kp + 26);
    const auto *kp_27 = buffer.data(kp + 27);
    const auto *kp_28 = buffer.data(kp + 28);
    const auto *kp_29 = buffer.data(kp + 29);
    const auto *kp_30 = buffer.data(kp + 30);
    const auto *kp_31 = buffer.data(kp + 31);
    const auto *kp_32 = buffer.data(kp + 32);
    const auto *kp_35 = buffer.data(kp + 35);
    const auto *kp_37 = buffer.data(kp + 37);
    const auto *kp_38 = buffer.data(kp + 38);
    const auto *kp_40 = buffer.data(kp + 40);
    const auto *kp_41 = buffer.data(kp + 41);
    const auto *kp_42 = buffer.data(kp + 42);
    const auto *kp_43 = buffer.data(kp + 43);
    const auto *kp_44 = buffer.data(kp + 44);
    const auto *kp_45 = buffer.data(kp + 45);
    const auto *kp_46 = buffer.data(kp + 46);
    const auto *kp_47 = buffer.data(kp + 47);
    const auto *kp_50 = buffer.data(kp + 50);
    const auto *kp_52 = buffer.data(kp + 52);
    const auto *kp_53 = buffer.data(kp + 53);
    const auto *kp_55 = buffer.data(kp + 55);
    const auto *kp_56 = buffer.data(kp + 56);
    const auto *kp_58 = buffer.data(kp + 58);
    const auto *kp_59 = buffer.data(kp + 59);
    const auto *kp_60 = buffer.data(kp + 60);
    const auto *kp_61 = buffer.data(kp + 61);
    const auto *kp_62 = buffer.data(kp + 62);
    const auto *kp_63 = buffer.data(kp + 63);
    const auto *kp_64 = buffer.data(kp + 64);
    const auto *kp_68 = buffer.data(kp + 68);
    const auto *kp_70 = buffer.data(kp + 70);
    const auto *kp_71 = buffer.data(kp + 71);
    const auto *kp_73 = buffer.data(kp + 73);
    const auto *kp_74 = buffer.data(kp + 74);
    const auto *kp_76 = buffer.data(kp + 76);
    const auto *kp_77 = buffer.data(kp + 77);
    const auto *kp_79 = buffer.data(kp + 79);
    const auto *kp_81 = buffer.data(kp + 81);
    const auto *kp_83 = buffer.data(kp + 83);
    const auto *kp_84 = buffer.data(kp + 84);
    const auto *kp_85 = buffer.data(kp + 85);
    const auto *kp_86 = buffer.data(kp + 86);
    const auto *kp_88 = buffer.data(kp + 88);
    const auto *kp_89 = buffer.data(kp + 89);
    const auto *kp_90 = buffer.data(kp + 90);
    const auto *kp_91 = buffer.data(kp + 91);
    const auto *kp_92 = buffer.data(kp + 92);
    const auto *kp_93 = buffer.data(kp + 93);
    const auto *kp_94 = buffer.data(kp + 94);
    const auto *kp_95 = buffer.data(kp + 95);
    const auto *kp_96 = buffer.data(kp + 96);
    const auto *kp_97 = buffer.data(kp + 97);
    const auto *kp_98 = buffer.data(kp + 98);
    const auto *kp_99 = buffer.data(kp + 99);
    const auto *kp_100 = buffer.data(kp + 100);
    const auto *kp_101 = buffer.data(kp + 101);
    const auto *kp_103 = buffer.data(kp + 103);
    const auto *kp_104 = buffer.data(kp + 104);
    const auto *kp_105 = buffer.data(kp + 105);
    const auto *kp_106 = buffer.data(kp + 106);
    const auto *kp_107 = buffer.data(kp + 107);

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

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pa_y, pb_x, pb_z, ip_1, ip_4, id_0, \
                         id_3, id_5, kp_3, kp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * id_0[k];

        t_7[k] = f_3 * ip_4[k]
                 + pb_x[k] * kp_4[k];

        t_8[k] = pb_z[k] * kp_3[k];

        t_9[k] = f_4 * ip_1[k]
                 + pa_y[k] * id_3[k];

        t_10[k] = pb_z[k] * kp_4[k];

        t_11[k] = pa_y[k] * id_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, ip_2, ip_8, \
                         id_0, id_3, id_5, kp_6, kp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_z[k] * id_0[k];

        t_13[k] = pb_y[k] * kp_6[k];

        t_14[k] = f_3 * ip_8[k]
                  + pb_x[k] * kp_8[k];

        t_15[k] = pa_z[k] * id_3[k];

        t_16[k] = pb_y[k] * kp_8[k];

        t_17[k] = f_4 * ip_2[k]
                  + pa_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pb_x, pb_z, hd0_0, hd1_0, ip_10, id_6, kp_9, \
                         kp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * hd0_0[k]
                  - f_6 * hd1_0[k]
                  + pa_y[k] * id_6[k];

        t_19[k] = f_7 * ip_10[k]
                  + pb_x[k] * kp_10[k];

        t_20[k] = pb_z[k] * kp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_x, pa_y, pb_z, hd0_21, hd1_21, id_12, \
                         id_21, ks0_3, ks1_3, kp_10, kp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * hd0_21[k]
                  - f_9 * hd1_21[k]
                  + pa_x[k] * id_21[k];

        t_22[k] = pb_z[k] * kp_10[k];

        t_23[k] = f_1 * ks0_3[k]
                  - f_2 * ks1_3[k]
                  + pb_z[k] * kp_11[k];

        t_24[k] = pa_y[k] * id_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pa_z, pb_y, ip_8, id_7, id_9, \
                         id_14, id_17, kp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_z[k] * id_7[k];

        t_26[k] = pa_y[k] * id_14[k];

        t_27[k] = pa_z[k] * id_9[k];

        t_28[k] = f_10 * ip_8[k]
                  + pb_y[k] * kp_14[k];

        t_29[k] = pa_y[k] * id_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_z, pb_x, pb_y, hd0_0, hd1_0, ip_17, id_12, \
                         ks0_5, ks1_5, kp_15, kp_16, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * hd0_0[k]
                  - f_6 * hd1_0[k]
                  + pa_z[k] * id_12[k];

        t_31[k] = pb_y[k] * kp_15[k];

        t_32[k] = f_7 * ip_17[k]
                  + pb_x[k] * kp_17[k];

        t_33[k] = f_1 * ks0_5[k]
                  - f_2 * ks1_5[k]
                  + pb_y[k] * kp_16[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_x, pa_y, pb_y, hd0_6, hd0_35, hd1_6, hd1_35, \
                         id_18, id_35, kp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * kp_17[k];

        t_35[k] = f_8 * hd0_35[k]
                  - f_9 * hd1_35[k]
                  + pa_x[k] * id_35[k];

        t_36[k] = f_11 * hd0_6[k]
                  - f_12 * hd1_6[k]
                  + pa_y[k] * id_18[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_x, pb_x, pb_z, hd0_39, hd1_39, ip_19, \
                         id_39, kp_18, kp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_13 * ip_19[k]
                  + pb_x[k] * kp_19[k];

        t_38[k] = pb_z[k] * kp_18[k];

        t_39[k] = f_14 * hd0_39[k]
                  - f_15 * hd1_39[k]
                  + pa_x[k] * id_39[k];

        t_40[k] = pb_z[k] * kp_19[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_z, pb_x, pb_z, ip_23, id_18, id_19, \
                         id_21, ks0_6, ks1_6, kp_20, kp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ks0_6[k]
                  - f_2 * ks1_6[k]
                  + pb_z[k] * kp_20[k];

        t_42[k] = pa_z[k] * id_18[k];

        t_43[k] = pa_z[k] * id_19[k];

        t_44[k] = f_13 * ip_23[k]
                  + pb_x[k] * kp_23[k];

        t_45[k] = pa_z[k] * id_21[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pa_z, pb_x, pb_y, ip_11, ip_14, ip_25, \
                         id_23, id_30, kp_23, kp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_4 * ip_14[k]
                  + pb_y[k] * kp_23[k];

        t_47[k] = f_4 * ip_11[k]
                  + pa_z[k] * id_23[k];

        t_48[k] = pa_y[k] * id_30[k];

        t_49[k] = f_13 * ip_25[k]
                  + pb_x[k] * kp_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pb_y, ip_16, ip_17, id_32, id_33, \
                         id_35, kp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * id_32[k];

        t_51[k] = f_4 * ip_16[k]
                  + pa_y[k] * id_33[k];

        t_52[k] = f_10 * ip_17[k]
                  + pb_y[k] * kp_26[k];

        t_53[k] = pa_y[k] * id_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_z, pb_x, pb_y, hd0_12, hd1_12, ip_29, \
                         id_30, ks0_9, ks1_9, kp_27, kp_28, kp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * hd0_12[k]
                  - f_12 * hd1_12[k]
                  + pa_z[k] * id_30[k];

        t_55[k] = pb_y[k] * kp_27[k];

        t_56[k] = f_13 * ip_29[k]
                  + pb_x[k] * kp_29[k];

        t_57[k] = f_1 * ks0_9[k]
                  - f_2 * ks1_9[k]
                  + pb_y[k] * kp_28[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_x, pa_y, pb_y, hd0_18, hd0_59, hd1_18, hd1_59, \
                         id_36, id_59, kp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_y[k] * kp_29[k];

        t_59[k] = f_14 * hd0_59[k]
                  - f_15 * hd1_59[k]
                  + pa_x[k] * id_59[k];

        t_60[k] = f_14 * hd0_18[k]
                  - f_15 * hd1_18[k]
                  + pa_y[k] * id_36[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pb_x, pb_z, hd0_63, hd1_63, ip_31, \
                         id_63, kp_30, kp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_16 * ip_31[k]
                  + pb_x[k] * kp_31[k];

        t_62[k] = pb_z[k] * kp_30[k];

        t_63[k] = f_11 * hd0_63[k]
                  - f_12 * hd1_63[k]
                  + pa_x[k] * id_63[k];

        t_64[k] = pb_z[k] * kp_31[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pa_z, pb_x, pb_z, ip_35, id_36, id_37, \
                         id_39, ks0_10, ks1_10, kp_32, kp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_1 * ks0_10[k]
                  - f_2 * ks1_10[k]
                  + pb_z[k] * kp_32[k];

        t_66[k] = pa_z[k] * id_36[k];

        t_67[k] = pa_z[k] * id_37[k];

        t_68[k] = f_16 * ip_35[k]
                  + pb_x[k] * kp_35[k];

        t_69[k] = pa_z[k] * id_39[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_y, pa_z, pb_y, hd0_30, hd1_30, ip_20, ip_23, \
                         id_41, id_48, kp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_16 * ip_23[k]
                  + pb_y[k] * kp_35[k];

        t_71[k] = f_4 * ip_20[k]
                  + pa_z[k] * id_41[k];

        t_72[k] = f_5 * hd0_30[k]
                  - f_6 * hd1_30[k]
                  + pa_y[k] * id_48[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pb_x, pb_y, hd0_75, hd1_75, ip_26, \
                         ip_37, ip_38, id_75, kp_37, kp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_16 * ip_37[k]
                  + pb_x[k] * kp_37[k];

        t_74[k] = f_16 * ip_38[k]
                  + pb_x[k] * kp_38[k];

        t_75[k] = f_11 * hd0_75[k]
                  - f_12 * hd1_75[k]
                  + pa_x[k] * id_75[k];

        t_76[k] = f_4 * ip_26[k]
                  + pb_y[k] * kp_38[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_x, pa_y, pb_x, hd0_77, hd1_77, ip_40, \
                         id_54, id_56, id_77, kp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_11 * hd0_77[k]
                  - f_12 * hd1_77[k]
                  + pa_x[k] * id_77[k];

        t_78[k] = pa_y[k] * id_54[k];

        t_79[k] = f_16 * ip_40[k]
                  + pb_x[k] * kp_40[k];

        t_80[k] = pa_y[k] * id_56[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_y, pa_z, pb_y, hd0_30, hd1_30, ip_28, \
                         ip_29, id_54, id_57, id_59, kp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * ip_28[k]
                  + pa_y[k] * id_57[k];

        t_82[k] = f_10 * ip_29[k]
                  + pb_y[k] * kp_41[k];

        t_83[k] = pa_y[k] * id_59[k];

        t_84[k] = f_14 * hd0_30[k]
                  - f_15 * hd1_30[k]
                  + pa_z[k] * id_54[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_x, pb_y, ip_44, ks0_14, ks1_14, kp_42, \
                         kp_43, kp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pb_y[k] * kp_42[k];

        t_86[k] = f_16 * ip_44[k]
                  + pb_x[k] * kp_44[k];

        t_87[k] = f_1 * ks0_14[k]
                  - f_2 * ks1_14[k]
                  + pb_y[k] * kp_43[k];

        t_88[k] = pb_y[k] * kp_44[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_x, pa_y, pb_x, hd0_36, hd0_89, hd1_36, hd1_89, \
                         ip_46, id_60, id_89, kp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_11 * hd0_89[k]
                  - f_12 * hd1_89[k]
                  + pa_x[k] * id_89[k];

        t_90[k] = f_8 * hd0_36[k]
                  - f_9 * hd1_36[k]
                  + pa_y[k] * id_60[k];

        t_91[k] = f_4 * ip_46[k]
                  + pb_x[k] * kp_46[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_x, pb_z, hd0_93, hd1_93, id_93, ks0_15, \
                         ks1_15, kp_45, kp_46, kp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_z[k] * kp_45[k];

        t_93[k] = f_5 * hd0_93[k]
                  - f_6 * hd1_93[k]
                  + pa_x[k] * id_93[k];

        t_94[k] = pb_z[k] * kp_46[k];

        t_95[k] = f_1 * ks0_15[k]
                  - f_2 * ks1_15[k]
                  + pb_z[k] * kp_47[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pa_z, pb_x, pb_y, ip_35, ip_50, id_60, \
                         id_61, id_63, kp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * id_60[k];

        t_97[k] = pa_z[k] * id_61[k];

        t_98[k] = f_4 * ip_50[k]
                  + pb_x[k] * kp_50[k];

        t_99[k] = pa_z[k] * id_63[k];

        t_100[k] = f_13 * ip_35[k]
                   + pb_y[k] * kp_50[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_y, pa_z, pb_x, hd0_48, hd1_48, ip_32, \
                         ip_52, ip_53, id_65, id_72, kp_52, kp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_4 * ip_32[k]
                   + pa_z[k] * id_65[k];

        t_102[k] = f_11 * hd0_48[k]
                   - f_12 * hd1_48[k]
                   + pa_y[k] * id_72[k];

        t_103[k] = f_4 * ip_52[k]
                   + pb_x[k] * kp_52[k];

        t_104[k] = f_4 * ip_53[k]
                   + pb_x[k] * kp_53[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pa_x, pb_y, hd0_105, hd0_107, hd1_105, hd1_107, \
                         ip_38, id_105, id_107, kp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_5 * hd0_105[k]
                   - f_6 * hd1_105[k]
                   + pa_x[k] * id_105[k];

        t_106[k] = f_16 * ip_38[k]
                   + pb_y[k] * kp_53[k];

        t_107[k] = f_5 * hd0_107[k]
                   - f_6 * hd1_107[k]
                   + pa_x[k] * id_107[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_y, pb_x, hd0_54, hd1_54, ip_55, ip_56, id_78, \
                         kp_55, kp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_5 * hd0_54[k]
                   - f_6 * hd1_54[k]
                   + pa_y[k] * id_78[k];

        t_109[k] = f_4 * ip_55[k]
                   + pb_x[k] * kp_55[k];

        t_110[k] = f_4 * ip_56[k]
                   + pb_x[k] * kp_56[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_x, pa_y, pb_y, hd0_111, hd0_113, \
                         hd1_111, hd1_113, ip_41, id_84, id_111, id_113, \
                         kp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_5 * hd0_111[k]
                   - f_6 * hd1_111[k]
                   + pa_x[k] * id_111[k];

        t_112[k] = f_4 * ip_41[k]
                   + pb_y[k] * kp_56[k];

        t_113[k] = f_5 * hd0_113[k]
                   - f_6 * hd1_113[k]
                   + pa_x[k] * id_113[k];

        t_114[k] = pa_y[k] * id_84[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pa_y, pb_x, pb_y, ip_43, ip_44, \
                         ip_58, id_86, id_87, id_89, kp_58, kp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_4 * ip_58[k]
                   + pb_x[k] * kp_58[k];

        t_116[k] = pa_y[k] * id_86[k];

        t_117[k] = f_4 * ip_43[k]
                   + pa_y[k] * id_87[k];

        t_118[k] = f_10 * ip_44[k]
                   + pb_y[k] * kp_59[k];

        t_119[k] = pa_y[k] * id_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_z, pb_x, pb_y, hd0_54, hd1_54, ip_62, \
                         id_84, ks0_20, ks1_20, kp_60, kp_61, kp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_8 * hd0_54[k]
                   - f_9 * hd1_54[k]
                   + pa_z[k] * id_84[k];

        t_121[k] = pb_y[k] * kp_60[k];

        t_122[k] = f_4 * ip_62[k]
                   + pb_x[k] * kp_62[k];

        t_123[k] = f_1 * ks0_20[k]
                   - f_2 * ks1_20[k]
                   + pb_y[k] * kp_61[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_x, pb_x, pb_y, hd0_125, hd1_125, \
                         ip_63, ip_64, id_125, id_126, kp_62, kp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * kp_62[k];

        t_125[k] = f_5 * hd0_125[k]
                   - f_6 * hd1_125[k]
                   + pa_x[k] * id_125[k];

        t_126[k] = f_4 * ip_63[k]
                   + pa_x[k] * id_126[k];

        t_127[k] = f_10 * ip_64[k]
                   + pb_x[k] * kp_64[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, pa_x, pa_z, pb_z, id_90, \
                         id_91, id_129, id_131, kp_63, kp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = pb_z[k] * kp_63[k];

        t_129[k] = pa_x[k] * id_129[k];

        t_130[k] = pb_z[k] * kp_64[k];

        t_131[k] = pa_x[k] * id_131[k];

        t_132[k] = pa_z[k] * id_90[k];

        t_133[k] = pa_z[k] * id_91[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, pa_x, pb_x, ip_68, ip_69, id_135, \
                         id_136, id_137, id_138, kp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_10 * ip_68[k]
                   + pb_x[k] * kp_68[k];

        t_135[k] = pa_x[k] * id_135[k];

        t_136[k] = pa_x[k] * id_136[k];

        t_137[k] = pa_x[k] * id_137[k];

        t_138[k] = f_4 * ip_69[k]
                   + pa_x[k] * id_138[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pa_x, pb_x, ip_70, ip_71, id_141, \
                         id_142, id_143, kp_70, kp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_10 * ip_70[k]
                   + pb_x[k] * kp_70[k];

        t_140[k] = f_10 * ip_71[k]
                   + pb_x[k] * kp_71[k];

        t_141[k] = pa_x[k] * id_141[k];

        t_142[k] = pa_x[k] * id_142[k];

        t_143[k] = pa_x[k] * id_143[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, pa_x, pb_x, ip_72, ip_73, ip_74, \
                         id_144, id_147, id_148, kp_73, kp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_4 * ip_72[k]
                   + pa_x[k] * id_144[k];

        t_145[k] = f_10 * ip_73[k]
                   + pb_x[k] * kp_73[k];

        t_146[k] = f_10 * ip_74[k]
                   + pb_x[k] * kp_74[k];

        t_147[k] = pa_x[k] * id_147[k];

        t_148[k] = pa_x[k] * id_148[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, pa_x, pb_x, ip_75, ip_76, ip_77, \
                         id_149, id_150, id_153, kp_76, kp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pa_x[k] * id_149[k];

        t_150[k] = f_4 * ip_75[k]
                   + pa_x[k] * id_150[k];

        t_151[k] = f_10 * ip_76[k]
                   + pb_x[k] * kp_76[k];

        t_152[k] = f_10 * ip_77[k]
                   + pb_x[k] * kp_77[k];

        t_153[k] = pa_x[k] * id_153[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, t_159, pa_x, pa_y, pb_x, ip_79, \
                         id_120, id_122, id_154, id_155, id_159, \
                         kp_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pa_x[k] * id_154[k];

        t_155[k] = pa_x[k] * id_155[k];

        t_156[k] = pa_y[k] * id_120[k];

        t_157[k] = f_10 * ip_79[k]
                   + pb_x[k] * kp_79[k];

        t_158[k] = pa_y[k] * id_122[k];

        t_159[k] = pa_x[k] * id_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, pa_x, pb_x, pb_y, ip_81, ip_83, \
                         id_160, id_161, id_162, kp_81, kp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = pa_x[k] * id_160[k];

        t_161[k] = pa_x[k] * id_161[k];

        t_162[k] = f_4 * ip_81[k]
                   + pa_x[k] * id_162[k];

        t_163[k] = pb_y[k] * kp_81[k];

        t_164[k] = f_10 * ip_83[k]
                   + pb_x[k] * kp_83[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pa_x, pb_x, pb_y, id_165, id_167, \
                         ks0_28, ks1_28, kp_83, kp_84, kp_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = pa_x[k] * id_165[k];

        t_166[k] = pb_y[k] * kp_83[k];

        t_167[k] = pa_x[k] * id_167[k];

        t_168[k] = f_1 * ks0_28[k]
                   - f_2 * ks1_28[k]
                   + pb_x[k] * kp_84[k];

        t_169[k] = pb_x[k] * kp_85[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, pa_z, pb_x, pb_y, pb_z, ip_64, \
                         id_126, ks0_28, ks1_28, kp_85, kp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = pb_x[k] * kp_86[k];

        t_171[k] = f_0 * ip_64[k]
                   + f_1 * ks0_28[k]
                   - f_2 * ks1_28[k]
                   + pb_y[k] * kp_85[k];

        t_172[k] = pb_z[k] * kp_85[k];

        t_173[k] = f_1 * ks0_28[k]
                   - f_2 * ks1_28[k]
                   + pb_z[k] * kp_86[k];

        t_174[k] = pa_z[k] * id_126[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, pa_z, pb_x, pb_y, ip_65, ip_68, \
                         id_129, id_131, kp_88, kp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pb_x[k] * kp_88[k];

        t_176[k] = pb_x[k] * kp_89[k];

        t_177[k] = pa_z[k] * id_129[k];

        t_178[k] = f_3 * ip_68[k]
                   + pb_y[k] * kp_89[k];

        t_179[k] = f_4 * ip_65[k]
                   + pa_z[k] * id_131[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_z, pb_x, hd0_93, hd1_93, id_135, \
                         ks0_30, ks1_30, kp_90, kp_91, kp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_1 * ks0_30[k]
                   - f_2 * ks1_30[k]
                   + pb_x[k] * kp_90[k];

        t_181[k] = pb_x[k] * kp_91[k];

        t_182[k] = pb_x[k] * kp_92[k];

        t_183[k] = f_5 * hd0_93[k]
                   - f_6 * hd1_93[k]
                   + pa_z[k] * id_135[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pb_x, pb_y, hd0_107, hd1_107, \
                         ip_71, id_143, ks0_31, ks1_31, kp_92, kp_93, \
                         kp_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_7 * ip_71[k]
                   + pb_y[k] * kp_92[k];

        t_185[k] = f_8 * hd0_107[k]
                   - f_9 * hd1_107[k]
                   + pa_y[k] * id_143[k];

        t_186[k] = f_1 * ks0_31[k]
                   - f_2 * ks1_31[k]
                   + pb_x[k] * kp_93[k];

        t_187[k] = pb_x[k] * kp_94[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_y, pa_z, pb_x, pb_y, hd0_99, hd0_113, \
                         hd1_99, hd1_113, ip_74, id_141, id_149, \
                         kp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pb_x[k] * kp_95[k];

        t_189[k] = f_11 * hd0_99[k]
                   - f_12 * hd1_99[k]
                   + pa_z[k] * id_141[k];

        t_190[k] = f_13 * ip_74[k]
                   + pb_y[k] * kp_95[k];

        t_191[k] = f_14 * hd0_113[k]
                   - f_15 * hd1_113[k]
                   + pa_y[k] * id_149[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_z, pb_x, hd0_105, hd1_105, id_147, \
                         ks0_32, ks1_32, kp_96, kp_97, kp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_1 * ks0_32[k]
                   - f_2 * ks1_32[k]
                   + pb_x[k] * kp_96[k];

        t_193[k] = pb_x[k] * kp_97[k];

        t_194[k] = pb_x[k] * kp_98[k];

        t_195[k] = f_14 * hd0_105[k]
                   - f_15 * hd1_105[k]
                   + pa_z[k] * id_147[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_y, pb_x, pb_y, hd0_119, hd1_119, \
                         ip_77, id_155, ks0_33, ks1_33, kp_98, kp_99, \
                         kp_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_16 * ip_77[k]
                   + pb_y[k] * kp_98[k];

        t_197[k] = f_11 * hd0_119[k]
                   - f_12 * hd1_119[k]
                   + pa_y[k] * id_155[k];

        t_198[k] = f_1 * ks0_33[k]
                   - f_2 * ks1_33[k]
                   + pb_x[k] * kp_99[k];

        t_199[k] = pb_x[k] * kp_100[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_y, pa_z, pb_x, pb_y, hd0_111, hd0_125, \
                         hd1_111, hd1_125, ip_80, id_153, id_161, \
                         kp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pb_x[k] * kp_101[k];

        t_201[k] = f_8 * hd0_111[k]
                   - f_9 * hd1_111[k]
                   + pa_z[k] * id_153[k];

        t_202[k] = f_4 * ip_80[k]
                   + pb_y[k] * kp_101[k];

        t_203[k] = f_5 * hd0_125[k]
                   - f_6 * hd1_125[k]
                   + pa_y[k] * id_161[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, t_209, pa_y, pb_x, pb_y, ip_82, \
                         ip_83, id_162, id_165, id_167, kp_103, \
                         kp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_y[k] * id_162[k];

        t_205[k] = pb_x[k] * kp_103[k];

        t_206[k] = pb_x[k] * kp_104[k];

        t_207[k] = f_4 * ip_82[k]
                   + pa_y[k] * id_165[k];

        t_208[k] = f_10 * ip_83[k]
                   + pb_y[k] * kp_104[k];

        t_209[k] = pa_y[k] * id_167[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, pb_x, pb_y, pb_z, ip_83, \
                         ks0_35, ks1_35, kp_105, kp_106, kp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_1 * ks0_35[k]
                   - f_2 * ks1_35[k]
                   + pb_x[k] * kp_105[k];

        t_211[k] = pb_x[k] * kp_106[k];

        t_212[k] = pb_x[k] * kp_107[k];

        t_213[k] = f_1 * ks0_35[k]
                   - f_2 * ks1_35[k]
                   + pb_y[k] * kp_106[k];

        t_214[k] = pb_y[k] * kp_107[k];

        t_215[k] = f_0 * ip_83[k]
                   + f_1 * ks0_35[k]
                   - f_2 * ks1_35[k]
                   + pb_z[k] * kp_107[k];
    }
}

}  // namespace simdt2ceri
