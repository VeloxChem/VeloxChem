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
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_21 = buffer.data(pi + 21);
    const auto *pi_23 = buffer.data(pi + 23);
    const auto *pi_24 = buffer.data(pi + 24);
    const auto *pi_25 = buffer.data(pi + 25);
    const auto *pi_27 = buffer.data(pi + 27);
    const auto *pi_28 = buffer.data(pi + 28);
    const auto *pi_31 = buffer.data(pi + 31);
    const auto *pi_33 = buffer.data(pi + 33);
    const auto *pi_34 = buffer.data(pi + 34);
    const auto *pi_37 = buffer.data(pi + 37);
    const auto *pi_38 = buffer.data(pi + 38);
    const auto *pi_40 = buffer.data(pi + 40);
    const auto *pi_42 = buffer.data(pi + 42);
    const auto *pi_43 = buffer.data(pi + 43);
    const auto *pi_45 = buffer.data(pi + 45);
    const auto *pi_46 = buffer.data(pi + 46);
    const auto *pi_49 = buffer.data(pi + 49);
    const auto *pi_51 = buffer.data(pi + 51);
    const auto *pi_52 = buffer.data(pi + 52);
    const auto *pi_53 = buffer.data(pi + 53);
    const auto *pi_54 = buffer.data(pi + 54);
    const auto *pi_55 = buffer.data(pi + 55);
    const auto *pi_56 = buffer.data(pi + 56);
    const auto *pi_58 = buffer.data(pi + 58);
    const auto *pi_59 = buffer.data(pi + 59);
    const auto *pi_61 = buffer.data(pi + 61);
    const auto *pi_62 = buffer.data(pi + 62);
    const auto *pi_64 = buffer.data(pi + 64);
    const auto *pi_65 = buffer.data(pi + 65);
    const auto *pi_66 = buffer.data(pi + 66);
    const auto *pi_68 = buffer.data(pi + 68);
    const auto *pi_69 = buffer.data(pi + 69);
    const auto *pi_70 = buffer.data(pi + 70);
    const auto *pi_73 = buffer.data(pi + 73);
    const auto *pi_74 = buffer.data(pi + 74);
    const auto *pi_76 = buffer.data(pi + 76);
    const auto *pi_77 = buffer.data(pi + 77);
    const auto *pi_78 = buffer.data(pi + 78);
    const auto *pi_79 = buffer.data(pi + 79);
    const auto *pi_80 = buffer.data(pi + 80);
    const auto *pi_81 = buffer.data(pi + 81);
    const auto *pi_82 = buffer.data(pi + 82);
    const auto *pi_83 = buffer.data(pi + 83);

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_3 = buffer.data(pk + 3);
    const auto *pk_5 = buffer.data(pk + 5);
    const auto *pk_6 = buffer.data(pk + 6);
    const auto *pk_9 = buffer.data(pk + 9);
    const auto *pk_10 = buffer.data(pk + 10);
    const auto *pk_14 = buffer.data(pk + 14);
    const auto *pk_15 = buffer.data(pk + 15);
    const auto *pk_20 = buffer.data(pk + 20);
    const auto *pk_21 = buffer.data(pk + 21);
    const auto *pk_27 = buffer.data(pk + 27);
    const auto *pk_37 = buffer.data(pk + 37);
    const auto *pk_39 = buffer.data(pk + 39);
    const auto *pk_42 = buffer.data(pk + 42);
    const auto *pk_46 = buffer.data(pk + 46);
    const auto *pk_48 = buffer.data(pk + 48);
    const auto *pk_51 = buffer.data(pk + 51);
    const auto *pk_53 = buffer.data(pk + 53);
    const auto *pk_54 = buffer.data(pk + 54);
    const auto *pk_64 = buffer.data(pk + 64);
    const auto *pk_66 = buffer.data(pk + 66);
    const auto *pk_67 = buffer.data(pk + 67);
    const auto *pk_68 = buffer.data(pk + 68);
    const auto *pk_69 = buffer.data(pk + 69);
    const auto *pk_70 = buffer.data(pk + 70);
    const auto *pk_71 = buffer.data(pk + 71);
    const auto *pk_72 = buffer.data(pk + 72);
    const auto *pk_74 = buffer.data(pk + 74);
    const auto *pk_77 = buffer.data(pk + 77);
    const auto *pk_81 = buffer.data(pk + 81);
    const auto *pk_84 = buffer.data(pk + 84);
    const auto *pk_86 = buffer.data(pk + 86);
    const auto *pk_89 = buffer.data(pk + 89);
    const auto *pk_90 = buffer.data(pk + 90);
    const auto *pk_92 = buffer.data(pk + 92);
    const auto *pk_100 = buffer.data(pk + 100);
    const auto *pk_101 = buffer.data(pk + 101);
    const auto *pk_102 = buffer.data(pk + 102);
    const auto *pk_103 = buffer.data(pk + 103);
    const auto *pk_104 = buffer.data(pk + 104);
    const auto *pk_105 = buffer.data(pk + 105);
    const auto *pk_107 = buffer.data(pk + 107);

    const auto *dh0_0 = buffer.data(dh0 + 0);
    const auto *dh0_1 = buffer.data(dh0 + 1);
    const auto *dh0_2 = buffer.data(dh0 + 2);
    const auto *dh0_3 = buffer.data(dh0 + 3);
    const auto *dh0_5 = buffer.data(dh0 + 5);
    const auto *dh0_6 = buffer.data(dh0 + 6);
    const auto *dh0_8 = buffer.data(dh0 + 8);
    const auto *dh0_9 = buffer.data(dh0 + 9);
    const auto *dh0_15 = buffer.data(dh0 + 15);
    const auto *dh0_17 = buffer.data(dh0 + 17);
    const auto *dh0_18 = buffer.data(dh0 + 18);
    const auto *dh0_19 = buffer.data(dh0 + 19);
    const auto *dh0_20 = buffer.data(dh0 + 20);
    const auto *dh0_63 = buffer.data(dh0 + 63);
    const auto *dh0_66 = buffer.data(dh0 + 66);
    const auto *dh0_68 = buffer.data(dh0 + 68);
    const auto *dh0_69 = buffer.data(dh0 + 69);
    const auto *dh0_72 = buffer.data(dh0 + 72);
    const auto *dh0_73 = buffer.data(dh0 + 73);
    const auto *dh0_75 = buffer.data(dh0 + 75);
    const auto *dh0_77 = buffer.data(dh0 + 77);
    const auto *dh0_78 = buffer.data(dh0 + 78);
    const auto *dh0_79 = buffer.data(dh0 + 79);
    const auto *dh0_80 = buffer.data(dh0 + 80);
    const auto *dh0_81 = buffer.data(dh0 + 81);
    const auto *dh0_83 = buffer.data(dh0 + 83);
    const auto *dh0_105 = buffer.data(dh0 + 105);
    const auto *dh0_108 = buffer.data(dh0 + 108);
    const auto *dh0_110 = buffer.data(dh0 + 110);
    const auto *dh0_111 = buffer.data(dh0 + 111);
    const auto *dh0_114 = buffer.data(dh0 + 114);
    const auto *dh0_115 = buffer.data(dh0 + 115);
    const auto *dh0_117 = buffer.data(dh0 + 117);
    const auto *dh0_119 = buffer.data(dh0 + 119);
    const auto *dh0_120 = buffer.data(dh0 + 120);
    const auto *dh0_122 = buffer.data(dh0 + 122);
    const auto *dh0_123 = buffer.data(dh0 + 123);
    const auto *dh0_124 = buffer.data(dh0 + 124);
    const auto *dh0_125 = buffer.data(dh0 + 125);

    const auto *dh1_0 = buffer.data(dh1 + 0);
    const auto *dh1_1 = buffer.data(dh1 + 1);
    const auto *dh1_2 = buffer.data(dh1 + 2);
    const auto *dh1_3 = buffer.data(dh1 + 3);
    const auto *dh1_5 = buffer.data(dh1 + 5);
    const auto *dh1_6 = buffer.data(dh1 + 6);
    const auto *dh1_8 = buffer.data(dh1 + 8);
    const auto *dh1_9 = buffer.data(dh1 + 9);
    const auto *dh1_15 = buffer.data(dh1 + 15);
    const auto *dh1_17 = buffer.data(dh1 + 17);
    const auto *dh1_18 = buffer.data(dh1 + 18);
    const auto *dh1_19 = buffer.data(dh1 + 19);
    const auto *dh1_20 = buffer.data(dh1 + 20);
    const auto *dh1_63 = buffer.data(dh1 + 63);
    const auto *dh1_66 = buffer.data(dh1 + 66);
    const auto *dh1_68 = buffer.data(dh1 + 68);
    const auto *dh1_69 = buffer.data(dh1 + 69);
    const auto *dh1_72 = buffer.data(dh1 + 72);
    const auto *dh1_73 = buffer.data(dh1 + 73);
    const auto *dh1_75 = buffer.data(dh1 + 75);
    const auto *dh1_77 = buffer.data(dh1 + 77);
    const auto *dh1_78 = buffer.data(dh1 + 78);
    const auto *dh1_79 = buffer.data(dh1 + 79);
    const auto *dh1_80 = buffer.data(dh1 + 80);
    const auto *dh1_81 = buffer.data(dh1 + 81);
    const auto *dh1_83 = buffer.data(dh1 + 83);
    const auto *dh1_105 = buffer.data(dh1 + 105);
    const auto *dh1_108 = buffer.data(dh1 + 108);
    const auto *dh1_110 = buffer.data(dh1 + 110);
    const auto *dh1_111 = buffer.data(dh1 + 111);
    const auto *dh1_114 = buffer.data(dh1 + 114);
    const auto *dh1_115 = buffer.data(dh1 + 115);
    const auto *dh1_117 = buffer.data(dh1 + 117);
    const auto *dh1_119 = buffer.data(dh1 + 119);
    const auto *dh1_120 = buffer.data(dh1 + 120);
    const auto *dh1_122 = buffer.data(dh1 + 122);
    const auto *dh1_123 = buffer.data(dh1 + 123);
    const auto *dh1_124 = buffer.data(dh1 + 124);
    const auto *dh1_125 = buffer.data(dh1 + 125);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_23 = buffer.data(di + 23);
    const auto *di_24 = buffer.data(di + 24);
    const auto *di_25 = buffer.data(di + 25);
    const auto *di_26 = buffer.data(di + 26);
    const auto *di_27 = buffer.data(di + 27);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_51 = buffer.data(di + 51);
    const auto *di_52 = buffer.data(di + 52);
    const auto *di_53 = buffer.data(di + 53);
    const auto *di_54 = buffer.data(di + 54);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_59 = buffer.data(di + 59);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_62 = buffer.data(di + 62);
    const auto *di_65 = buffer.data(di + 65);
    const auto *di_66 = buffer.data(di + 66);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_76 = buffer.data(di + 76);
    const auto *di_78 = buffer.data(di + 78);
    const auto *di_79 = buffer.data(di + 79);
    const auto *di_80 = buffer.data(di + 80);
    const auto *di_81 = buffer.data(di + 81);
    const auto *di_83 = buffer.data(di + 83);
    const auto *di_84 = buffer.data(di + 84);
    const auto *di_85 = buffer.data(di + 85);
    const auto *di_87 = buffer.data(di + 87);
    const auto *di_89 = buffer.data(di + 89);
    const auto *di_90 = buffer.data(di + 90);
    const auto *di_93 = buffer.data(di + 93);
    const auto *di_94 = buffer.data(di + 94);
    const auto *di_96 = buffer.data(di + 96);
    const auto *di_98 = buffer.data(di + 98);
    const auto *di_99 = buffer.data(di + 99);
    const auto *di_101 = buffer.data(di + 101);
    const auto *di_102 = buffer.data(di + 102);
    const auto *di_104 = buffer.data(di + 104);
    const auto *di_105 = buffer.data(di + 105);
    const auto *di_106 = buffer.data(di + 106);
    const auto *di_107 = buffer.data(di + 107);
    const auto *di_108 = buffer.data(di + 108);
    const auto *di_109 = buffer.data(di + 109);
    const auto *di_110 = buffer.data(di + 110);
    const auto *di_111 = buffer.data(di + 111);
    const auto *di_114 = buffer.data(di + 114);
    const auto *di_115 = buffer.data(di + 115);
    const auto *di_117 = buffer.data(di + 117);
    const auto *di_118 = buffer.data(di + 118);
    const auto *di_121 = buffer.data(di + 121);
    const auto *di_122 = buffer.data(di + 122);
    const auto *di_126 = buffer.data(di + 126);
    const auto *di_133 = buffer.data(di + 133);
    const auto *di_134 = buffer.data(di + 134);
    const auto *di_135 = buffer.data(di + 135);
    const auto *di_136 = buffer.data(di + 136);
    const auto *di_137 = buffer.data(di + 137);
    const auto *di_138 = buffer.data(di + 138);
    const auto *di_139 = buffer.data(di + 139);
    const auto *di_140 = buffer.data(di + 140);
    const auto *di_142 = buffer.data(di + 142);
    const auto *di_143 = buffer.data(di + 143);
    const auto *di_145 = buffer.data(di + 145);
    const auto *di_146 = buffer.data(di + 146);
    const auto *di_149 = buffer.data(di + 149);
    const auto *di_150 = buffer.data(di + 150);
    const auto *di_152 = buffer.data(di + 152);
    const auto *di_154 = buffer.data(di + 154);
    const auto *di_155 = buffer.data(di + 155);
    const auto *di_157 = buffer.data(di + 157);
    const auto *di_158 = buffer.data(di + 158);
    const auto *di_160 = buffer.data(di + 160);
    const auto *di_161 = buffer.data(di + 161);
    const auto *di_162 = buffer.data(di + 162);
    const auto *di_163 = buffer.data(di + 163);
    const auto *di_164 = buffer.data(di + 164);
    const auto *di_165 = buffer.data(di + 165);
    const auto *di_166 = buffer.data(di + 166);
    const auto *di_167 = buffer.data(di + 167);

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
                         dh1_2, dh1_3, di_3, di_5, di_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * dh0_1[k]
                 - f_6 * dh1_1[k]
                 + pb_y[k] * di_3[k];

        t_7[k] = pb_z[k] * di_3[k];

        t_8[k] = pb_y[k] * di_5[k];

        t_9[k] = f_5 * dh0_2[k]
                 - f_6 * dh1_2[k]
                 + pb_z[k] * di_5[k];

        t_10[k] = f_7 * dh0_3[k]
                  - f_8 * dh1_3[k]
                  + pb_y[k] * di_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, dh0_5, dh0_6, dh1_5, \
                         dh1_6, di_6, di_8, di_9, di_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * di_6[k];

        t_12[k] = f_3 * dh0_5[k]
                  - f_4 * dh1_5[k]
                  + pb_y[k] * di_8[k];

        t_13[k] = pb_y[k] * di_9[k];

        t_14[k] = f_7 * dh0_5[k]
                  - f_8 * dh1_5[k]
                  + pb_z[k] * di_9[k];

        t_15[k] = f_9 * dh0_6[k]
                  - f_10 * dh1_6[k]
                  + pb_y[k] * di_10[k];

        t_16[k] = pb_z[k] * di_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, dh0_8, dh0_9, dh1_8, dh1_9, \
                         di_12, di_13, di_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * dh0_8[k]
                  - f_6 * dh1_8[k]
                  + pb_y[k] * di_12[k];

        t_18[k] = f_3 * dh0_9[k]
                  - f_4 * dh1_9[k]
                  + pb_y[k] * di_13[k];

        t_19[k] = pb_y[k] * di_14[k];

        t_20[k] = f_9 * dh0_9[k]
                  - f_10 * dh1_9[k]
                  + pb_z[k] * di_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_z, pi_21, pi_23, pi_24, pi_25, \
                         di_15, di_21, di_23, di_24, di_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * pi_21[k]
                  + pb_x[k] * di_21[k];

        t_22[k] = pb_z[k] * di_15[k];

        t_23[k] = f_0 * pi_23[k]
                  + pb_x[k] * di_23[k];

        t_24[k] = f_0 * pi_24[k]
                  + pb_x[k] * di_24[k];

        t_25[k] = f_0 * pi_25[k]
                  + pb_x[k] * di_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, pi_27, dh0_15, dh1_15, \
                         di_20, di_21, di_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_y[k] * di_20[k];

        t_27[k] = f_0 * pi_27[k]
                  + pb_x[k] * di_27[k];

        t_28[k] = f_1 * dh0_15[k]
                  - f_2 * dh1_15[k]
                  + pb_y[k] * di_21[k];

        t_29[k] = pb_z[k] * di_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_y, dh0_17, dh0_18, dh0_19, dh1_17, dh1_18, \
                         dh1_19, di_23, di_24, di_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_9 * dh0_17[k]
                  - f_10 * dh1_17[k]
                  + pb_y[k] * di_23[k];

        t_31[k] = f_7 * dh0_18[k]
                  - f_8 * dh1_18[k]
                  + pb_y[k] * di_24[k];

        t_32[k] = f_5 * dh0_19[k]
                  - f_6 * dh1_19[k]
                  + pb_y[k] * di_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pa_y, pb_y, pb_z, pi_0, pk_0, \
                         dh0_20, dh1_20, di_26, di_27, di_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * dh0_20[k]
                  - f_4 * dh1_20[k]
                  + pb_y[k] * di_26[k];

        t_34[k] = pb_y[k] * di_27[k];

        t_35[k] = f_1 * dh0_20[k]
                  - f_2 * dh1_20[k]
                  + pb_z[k] * di_27[k];

        t_36[k] = pa_y[k] * pk_0[k];

        t_37[k] = f_11 * pi_0[k]
                  + pb_y[k] * di_28[k];

        t_38[k] = pb_z[k] * di_28[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_x, pa_y, pb_z, pi_31, pi_34, pk_5, \
                         pk_39, pk_42, di_29, di_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * pi_31[k]
                  + pa_x[k] * pk_39[k];

        t_40[k] = pb_z[k] * di_29[k];

        t_41[k] = pa_y[k] * pk_5[k];

        t_42[k] = f_13 * pi_34[k]
                  + pa_x[k] * pk_42[k];

        t_43[k] = pb_z[k] * di_31[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pa_y, pb_y, pb_z, pi_5, pi_38, pk_9, \
                         pk_46, di_33, di_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_11 * pi_5[k]
                  + pb_y[k] * di_33[k];

        t_45[k] = pa_y[k] * pk_9[k];

        t_46[k] = f_14 * pi_38[k]
                  + pa_x[k] * pk_46[k];

        t_47[k] = pb_z[k] * di_34[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pa_y, pb_y, pi_9, pi_40, pi_43, pk_14, \
                         pk_48, pk_51, di_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_14 * pi_40[k]
                  + pa_x[k] * pk_48[k];

        t_49[k] = f_11 * pi_9[k]
                  + pb_y[k] * di_37[k];

        t_50[k] = pa_y[k] * pk_14[k];

        t_51[k] = f_0 * pi_43[k]
                  + pa_x[k] * pk_51[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_y, pb_z, pi_14, pi_45, pi_46, pk_53, \
                         pk_54, di_38, di_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_z[k] * di_38[k];

        t_53[k] = f_0 * pi_45[k]
                  + pa_x[k] * pk_53[k];

        t_54[k] = f_0 * pi_46[k]
                  + pa_x[k] * pk_54[k];

        t_55[k] = f_11 * pi_14[k]
                  + pb_y[k] * di_42[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_y, pb_x, pb_z, pi_49, pi_51, pi_52, \
                         pk_20, di_43, di_49, di_51, di_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_y[k] * pk_20[k];

        t_57[k] = f_11 * pi_49[k]
                  + pb_x[k] * di_49[k];

        t_58[k] = pb_z[k] * di_43[k];

        t_59[k] = f_11 * pi_51[k]
                  + pb_x[k] * di_51[k];

        t_60[k] = f_11 * pi_52[k]
                  + pb_x[k] * di_52[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pa_x, pa_y, pb_x, pb_z, pi_53, pi_54, \
                         pk_27, pk_64, di_49, di_53, di_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_11 * pi_53[k]
                  + pb_x[k] * di_53[k];

        t_62[k] = f_11 * pi_54[k]
                  + pb_x[k] * di_54[k];

        t_63[k] = pa_y[k] * pk_27[k];

        t_64[k] = pa_x[k] * pk_64[k];

        t_65[k] = pb_z[k] * di_49[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, t_72, pa_x, pa_z, pk_0, pk_66, \
                         pk_67, pk_68, pk_69, pk_70, pk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_x[k] * pk_66[k];

        t_67[k] = pa_x[k] * pk_67[k];

        t_68[k] = pa_x[k] * pk_68[k];

        t_69[k] = pa_x[k] * pk_69[k];

        t_70[k] = pa_x[k] * pk_70[k];

        t_71[k] = pa_x[k] * pk_71[k];

        t_72[k] = pa_z[k] * pk_0[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_x, pa_z, pb_y, pb_z, pi_0, pi_61, \
                         pk_3, pk_77, di_56, di_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pb_y[k] * di_56[k];

        t_74[k] = f_11 * pi_0[k]
                  + pb_z[k] * di_56[k];

        t_75[k] = pa_z[k] * pk_3[k];

        t_76[k] = pb_y[k] * di_58[k];

        t_77[k] = f_12 * pi_61[k]
                  + pa_x[k] * pk_77[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_x, pa_z, pb_y, pb_z, pi_3, pi_65, \
                         pk_6, pk_10, pk_81, di_59, di_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_z[k] * pk_6[k];

        t_79[k] = f_11 * pi_3[k]
                  + pb_z[k] * di_59[k];

        t_80[k] = pb_y[k] * di_61[k];

        t_81[k] = f_13 * pi_65[k]
                  + pa_x[k] * pk_81[k];

        t_82[k] = pa_z[k] * pk_10[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_x, pb_y, pb_z, pi_6, pi_68, pi_70, pk_84, \
                         pk_86, di_62, di_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_11 * pi_6[k]
                  + pb_z[k] * di_62[k];

        t_84[k] = f_14 * pi_68[k]
                  + pa_x[k] * pk_84[k];

        t_85[k] = pb_y[k] * di_65[k];

        t_86[k] = f_14 * pi_70[k]
                  + pa_x[k] * pk_86[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_x, pa_z, pb_z, pi_10, pi_73, pi_74, pk_15, \
                         pk_89, pk_90, di_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * pk_15[k];

        t_88[k] = f_11 * pi_10[k]
                  + pb_z[k] * di_66[k];

        t_89[k] = f_0 * pi_73[k]
                  + pa_x[k] * pk_89[k];

        t_90[k] = f_0 * pi_74[k]
                  + pa_x[k] * pk_90[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_x, pa_z, pb_x, pb_y, pi_76, pi_78, pk_21, \
                         pk_92, di_70, di_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = pb_y[k] * di_70[k];

        t_92[k] = f_0 * pi_76[k]
                  + pa_x[k] * pk_92[k];

        t_93[k] = pa_z[k] * pk_21[k];

        t_94[k] = f_11 * pi_78[k]
                  + pb_x[k] * di_78[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pb_x, pb_y, pi_79, pi_80, pi_81, pi_83, \
                         di_76, di_79, di_80, di_81, di_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_11 * pi_79[k]
                  + pb_x[k] * di_79[k];

        t_96[k] = f_11 * pi_80[k]
                  + pb_x[k] * di_80[k];

        t_97[k] = f_11 * pi_81[k]
                  + pb_x[k] * di_81[k];

        t_98[k] = pb_y[k] * di_76[k];

        t_99[k] = f_11 * pi_83[k]
                  + pb_x[k] * di_83[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, t_105, t_106, pa_x, pb_y, pk_100, \
                         pk_101, pk_102, pk_103, pk_104, pk_105, \
                         di_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_x[k] * pk_100[k];

        t_101[k] = pa_x[k] * pk_101[k];

        t_102[k] = pa_x[k] * pk_102[k];

        t_103[k] = pa_x[k] * pk_103[k];

        t_104[k] = pa_x[k] * pk_104[k];

        t_105[k] = pa_x[k] * pk_105[k];

        t_106[k] = pb_y[k] * di_83[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_x, pb_x, pb_y, pb_z, pi_28, pk_107, \
                         dh0_63, dh1_63, di_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pa_x[k] * pk_107[k];

        t_108[k] = f_1 * dh0_63[k]
                   - f_2 * dh1_63[k]
                   + pb_x[k] * di_84[k];

        t_109[k] = f_0 * pi_28[k]
                   + pb_y[k] * di_84[k];

        t_110[k] = pb_z[k] * di_84[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pb_x, pb_z, dh0_66, dh0_68, dh0_69, \
                         dh1_66, dh1_68, dh1_69, di_85, di_87, di_89, \
                         di_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_9 * dh0_66[k]
                   - f_10 * dh1_66[k]
                   + pb_x[k] * di_87[k];

        t_112[k] = pb_z[k] * di_85[k];

        t_113[k] = f_9 * dh0_68[k]
                   - f_10 * dh1_68[k]
                   + pb_x[k] * di_89[k];

        t_114[k] = f_7 * dh0_69[k]
                   - f_8 * dh1_69[k]
                   + pb_x[k] * di_90[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pb_x, pb_y, pb_z, pi_33, dh0_72, dh0_73, \
                         dh1_72, dh1_73, di_87, di_89, di_93, di_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pb_z[k] * di_87[k];

        t_116[k] = f_0 * pi_33[k]
                   + pb_y[k] * di_89[k];

        t_117[k] = f_7 * dh0_72[k]
                   - f_8 * dh1_72[k]
                   + pb_x[k] * di_93[k];

        t_118[k] = f_5 * dh0_73[k]
                   - f_6 * dh1_73[k]
                   + pb_x[k] * di_94[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pb_x, pb_y, pb_z, pi_37, dh0_75, dh0_77, \
                         dh1_75, dh1_77, di_90, di_93, di_96, di_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pb_z[k] * di_90[k];

        t_120[k] = f_5 * dh0_75[k]
                   - f_6 * dh1_75[k]
                   + pb_x[k] * di_96[k];

        t_121[k] = f_0 * pi_37[k]
                   + pb_y[k] * di_93[k];

        t_122[k] = f_5 * dh0_77[k]
                   - f_6 * dh1_77[k]
                   + pb_x[k] * di_98[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pb_x, pb_z, dh0_78, dh0_80, dh0_81, \
                         dh1_78, dh1_80, dh1_81, di_94, di_99, di_101, \
                         di_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_3 * dh0_78[k]
                   - f_4 * dh1_78[k]
                   + pb_x[k] * di_99[k];

        t_124[k] = pb_z[k] * di_94[k];

        t_125[k] = f_3 * dh0_80[k]
                   - f_4 * dh1_80[k]
                   + pb_x[k] * di_101[k];

        t_126[k] = f_3 * dh0_81[k]
                   - f_4 * dh1_81[k]
                   + pb_x[k] * di_102[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pb_x, pb_y, pi_42, dh0_83, dh1_83, \
                         di_98, di_104, di_105, di_106, di_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_0 * pi_42[k]
                   + pb_y[k] * di_98[k];

        t_128[k] = f_3 * dh0_83[k]
                   - f_4 * dh1_83[k]
                   + pb_x[k] * di_104[k];

        t_129[k] = pb_x[k] * di_105[k];

        t_130[k] = pb_x[k] * di_106[k];

        t_131[k] = pb_x[k] * di_107[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, pb_x, pb_y, pi_49, dh0_78, dh1_78, \
                         di_105, di_108, di_109, di_110, di_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pb_x[k] * di_108[k];

        t_133[k] = pb_x[k] * di_109[k];

        t_134[k] = pb_x[k] * di_110[k];

        t_135[k] = pb_x[k] * di_111[k];

        t_136[k] = f_0 * pi_49[k]
                   + f_1 * dh0_78[k]
                   - f_2 * dh1_78[k]
                   + pb_y[k] * di_105[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pb_z, dh0_78, dh0_79, dh0_80, dh1_78, \
                         dh1_79, dh1_80, di_105, di_106, di_107, \
                         di_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = pb_z[k] * di_105[k];

        t_138[k] = f_3 * dh0_78[k]
                   - f_4 * dh1_78[k]
                   + pb_z[k] * di_106[k];

        t_139[k] = f_5 * dh0_79[k]
                   - f_6 * dh1_79[k]
                   + pb_z[k] * di_107[k];

        t_140[k] = f_7 * dh0_80[k]
                   - f_8 * dh1_80[k]
                   + pb_z[k] * di_108[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_y, pb_y, pb_z, pi_55, pk_72, dh0_81, \
                         dh0_83, dh1_81, dh1_83, di_109, di_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_9 * dh0_81[k]
                   - f_10 * dh1_81[k]
                   + pb_z[k] * di_109[k];

        t_142[k] = f_0 * pi_55[k]
                   + pb_y[k] * di_111[k];

        t_143[k] = f_1 * dh0_83[k]
                   - f_2 * dh1_83[k]
                   + pb_z[k] * di_111[k];

        t_144[k] = pa_y[k] * pk_72[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, t_150, pa_y, pa_z, pb_y, pi_58, \
                         pk_37, pk_39, pk_42, pk_74, pk_77, di_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = pa_z[k] * pk_37[k];

        t_146[k] = pa_y[k] * pk_74[k];

        t_147[k] = pa_z[k] * pk_39[k];

        t_148[k] = f_11 * pi_58[k]
                   + pb_y[k] * di_114[k];

        t_149[k] = pa_y[k] * pk_77[k];

        t_150[k] = pa_z[k] * pk_42[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pa_y, pa_z, pb_y, pb_z, pi_31, pi_61, \
                         pk_46, pk_81, di_115, di_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_11 * pi_31[k]
                   + pb_z[k] * di_115[k];

        t_152[k] = f_11 * pi_61[k]
                   + pb_y[k] * di_117[k];

        t_153[k] = pa_y[k] * pk_81[k];

        t_154[k] = pa_z[k] * pk_46[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_y, pb_y, pb_z, pi_34, pi_64, pi_65, \
                         pk_84, pk_86, di_118, di_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_11 * pi_34[k]
                   + pb_z[k] * di_118[k];

        t_156[k] = f_0 * pi_64[k]
                   + pa_y[k] * pk_84[k];

        t_157[k] = f_11 * pi_65[k]
                   + pb_y[k] * di_121[k];

        t_158[k] = pa_y[k] * pk_86[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_y, pa_z, pb_z, pi_38, pi_68, pi_69, \
                         pk_51, pk_89, pk_90, di_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * pk_51[k];

        t_160[k] = f_11 * pi_38[k]
                   + pb_z[k] * di_122[k];

        t_161[k] = f_14 * pi_68[k]
                   + pa_y[k] * pk_89[k];

        t_162[k] = f_0 * pi_69[k]
                   + pa_y[k] * pk_90[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, t_168, pa_y, pb_x, pb_y, pi_70, \
                         pk_92, di_126, di_133, di_134, di_135, \
                         di_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_11 * pi_70[k]
                   + pb_y[k] * di_126[k];

        t_164[k] = pa_y[k] * pk_92[k];

        t_165[k] = pb_x[k] * di_133[k];

        t_166[k] = pb_x[k] * di_134[k];

        t_167[k] = pb_x[k] * di_135[k];

        t_168[k] = pb_x[k] * di_136[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, pa_z, pb_x, pb_z, pi_49, pk_64, \
                         di_133, di_137, di_138, di_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = pb_x[k] * di_137[k];

        t_170[k] = pb_x[k] * di_138[k];

        t_171[k] = pb_x[k] * di_139[k];

        t_172[k] = pa_z[k] * pk_64[k];

        t_173[k] = f_11 * pi_49[k]
                   + pb_z[k] * di_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pa_y, pi_79, pi_80, pi_81, pi_82, pk_102, \
                         pk_103, pk_104, pk_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_12 * pi_79[k]
                   + pa_y[k] * pk_102[k];

        t_175[k] = f_13 * pi_80[k]
                   + pa_y[k] * pk_103[k];

        t_176[k] = f_14 * pi_81[k]
                   + pa_y[k] * pk_104[k];

        t_177[k] = f_0 * pi_82[k]
                   + pa_y[k] * pk_105[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pa_y, pb_x, pb_y, pb_z, pi_56, \
                         pi_83, pk_107, dh0_105, dh1_105, di_139, \
                         di_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_11 * pi_83[k]
                   + pb_y[k] * di_139[k];

        t_179[k] = pa_y[k] * pk_107[k];

        t_180[k] = f_1 * dh0_105[k]
                   - f_2 * dh1_105[k]
                   + pb_x[k] * di_140[k];

        t_181[k] = pb_y[k] * di_140[k];

        t_182[k] = f_0 * pi_56[k]
                   + pb_z[k] * di_140[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, pb_x, pb_y, dh0_108, dh0_110, dh0_111, \
                         dh1_108, dh1_110, dh1_111, di_142, di_143, di_145, \
                         di_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * dh0_108[k]
                   - f_10 * dh1_108[k]
                   + pb_x[k] * di_143[k];

        t_184[k] = pb_y[k] * di_142[k];

        t_185[k] = f_9 * dh0_110[k]
                   - f_10 * dh1_110[k]
                   + pb_x[k] * di_145[k];

        t_186[k] = f_7 * dh0_111[k]
                   - f_8 * dh1_111[k]
                   + pb_x[k] * di_146[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, pb_x, pb_y, pb_z, pi_59, dh0_114, \
                         dh0_115, dh1_114, dh1_115, di_143, di_145, di_149, \
                         di_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_0 * pi_59[k]
                   + pb_z[k] * di_143[k];

        t_188[k] = pb_y[k] * di_145[k];

        t_189[k] = f_7 * dh0_114[k]
                   - f_8 * dh1_114[k]
                   + pb_x[k] * di_149[k];

        t_190[k] = f_5 * dh0_115[k]
                   - f_6 * dh1_115[k]
                   + pb_x[k] * di_150[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pb_x, pb_y, pb_z, pi_62, dh0_117, \
                         dh0_119, dh1_117, dh1_119, di_146, di_149, di_152, \
                         di_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_0 * pi_62[k]
                   + pb_z[k] * di_146[k];

        t_192[k] = f_5 * dh0_117[k]
                   - f_6 * dh1_117[k]
                   + pb_x[k] * di_152[k];

        t_193[k] = pb_y[k] * di_149[k];

        t_194[k] = f_5 * dh0_119[k]
                   - f_6 * dh1_119[k]
                   + pb_x[k] * di_154[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pb_x, pb_z, pi_66, dh0_120, dh0_122, dh1_120, \
                         dh1_122, di_150, di_155, di_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_3 * dh0_120[k]
                   - f_4 * dh1_120[k]
                   + pb_x[k] * di_155[k];

        t_196[k] = f_0 * pi_66[k]
                   + pb_z[k] * di_150[k];

        t_197[k] = f_3 * dh0_122[k]
                   - f_4 * dh1_122[k]
                   + pb_x[k] * di_157[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, pb_x, pb_y, dh0_123, dh0_125, \
                         dh1_123, dh1_125, di_154, di_158, di_160, di_161, \
                         di_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_3 * dh0_123[k]
                   - f_4 * dh1_123[k]
                   + pb_x[k] * di_158[k];

        t_199[k] = pb_y[k] * di_154[k];

        t_200[k] = f_3 * dh0_125[k]
                   - f_4 * dh1_125[k]
                   + pb_x[k] * di_160[k];

        t_201[k] = pb_x[k] * di_161[k];

        t_202[k] = pb_x[k] * di_162[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, t_208, pb_x, pb_y, dh0_120, \
                         dh1_120, di_161, di_163, di_164, di_165, di_166, \
                         di_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = pb_x[k] * di_163[k];

        t_204[k] = pb_x[k] * di_164[k];

        t_205[k] = pb_x[k] * di_165[k];

        t_206[k] = pb_x[k] * di_166[k];

        t_207[k] = pb_x[k] * di_167[k];

        t_208[k] = f_1 * dh0_120[k]
                   - f_2 * dh1_120[k]
                   + pb_y[k] * di_161[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, pb_y, pb_z, pi_77, dh0_122, dh0_123, dh1_122, \
                         dh1_123, di_161, di_163, di_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_0 * pi_77[k]
                   + pb_z[k] * di_161[k];

        t_210[k] = f_9 * dh0_122[k]
                   - f_10 * dh1_122[k]
                   + pb_y[k] * di_163[k];

        t_211[k] = f_7 * dh0_123[k]
                   - f_8 * dh1_123[k]
                   + pb_y[k] * di_164[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pb_y, pb_z, pi_83, dh0_124, dh0_125, \
                         dh1_124, dh1_125, di_165, di_166, di_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_5 * dh0_124[k]
                   - f_6 * dh1_124[k]
                   + pb_y[k] * di_165[k];

        t_213[k] = f_3 * dh0_125[k]
                   - f_4 * dh1_125[k]
                   + pb_y[k] * di_166[k];

        t_214[k] = pb_y[k] * di_167[k];

        t_215[k] = f_0 * pi_83[k]
                   + f_1 * dh0_125[k]
                   - f_2 * dh1_125[k]
                   + pb_z[k] * di_167[k];
    }
}

}  // namespace simdt2ceri
