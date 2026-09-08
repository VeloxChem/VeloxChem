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


#include "SimdElectronRepulsionVrrRecHF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_hf_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ff0, const size_t ff1,
                                     const size_t gd, const size_t gf, const size_t hp0,
                                     const size_t hp1, const size_t hd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);

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

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_8 = buffer.data(ff1 + 8);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
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
    const auto *gd_21 = buffer.data(gd + 21);
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
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_40 = buffer.data(gd + 40);
    const auto *gd_41 = buffer.data(gd + 41);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_43 = buffer.data(gd + 43);
    const auto *gd_44 = buffer.data(gd + 44);
    const auto *gd_45 = buffer.data(gd + 45);
    const auto *gd_46 = buffer.data(gd + 46);
    const auto *gd_47 = buffer.data(gd + 47);

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

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_55 = buffer.data(hd + 55);
    const auto *hd_56 = buffer.data(hd + 56);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_58 = buffer.data(hd + 58);
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_60 = buffer.data(hd + 60);
    const auto *hd_61 = buffer.data(hd + 61);
    const auto *hd_62 = buffer.data(hd + 62);
    const auto *hd_63 = buffer.data(hd + 63);
    const auto *hd_64 = buffer.data(hd + 64);
    const auto *hd_65 = buffer.data(hd + 65);
    const auto *hd_66 = buffer.data(hd + 66);
    const auto *hd_67 = buffer.data(hd + 67);
    const auto *hd_68 = buffer.data(hd + 68);
    const auto *hd_69 = buffer.data(hd + 69);
    const auto *hd_70 = buffer.data(hd + 70);
    const auto *hd_71 = buffer.data(hd + 71);
    const auto *hd_72 = buffer.data(hd + 72);
    const auto *hd_73 = buffer.data(hd + 73);
    const auto *hd_74 = buffer.data(hd + 74);
    const auto *hd_75 = buffer.data(hd + 75);
    const auto *hd_76 = buffer.data(hd + 76);
    const auto *hd_77 = buffer.data(hd + 77);
    const auto *hd_78 = buffer.data(hd + 78);
    const auto *hd_79 = buffer.data(hd + 79);
    const auto *hd_80 = buffer.data(hd + 80);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, gd_1, hp0_0, hp1_0, \
                         hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_0 * gd_1[k]
                 + pb_x[k] * hd_2[k];

        t_4[k] = pb_y[k] * hd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, gd_2, hp0_1, hp0_2, hp1_1, \
                         hp1_2, hd_2, hd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * gd_2[k]
                 + pb_x[k] * hd_3[k];

        t_6[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_2[k];

        t_7[k] = pb_z[k] * hd_2[k];

        t_8[k] = pb_y[k] * hd_3[k];

        t_9[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, pb_y, pb_z, gd_0, gd_4, \
                         gf_0, hd_4, hd_5, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * gf_0[k];

        t_11[k] = f_3 * gd_0[k]
                  + pb_y[k] * hd_4[k];

        t_12[k] = pb_z[k] * hd_4[k];

        t_13[k] = f_4 * gd_4[k]
                  + pb_x[k] * hd_6[k];

        t_14[k] = pb_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, gd_1, gd_2, gf_2, \
                         gf_3, gf_4, hd_6, hd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * gf_2[k];

        t_16[k] = f_5 * gd_1[k]
                  + pa_y[k] * gf_3[k];

        t_17[k] = pb_z[k] * hd_6[k];

        t_18[k] = f_3 * gd_2[k]
                  + pb_y[k] * hd_7[k];

        t_19[k] = pa_y[k] * gf_4[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, gd_0, gf_0, gf_1, \
                         hd_8, hd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * gf_0[k];

        t_21[k] = pb_y[k] * hd_8[k];

        t_22[k] = f_3 * gd_0[k]
                  + pb_z[k] * hd_8[k];

        t_23[k] = pa_z[k] * gf_1[k];

        t_24[k] = pb_y[k] * hd_9[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, pb_z, gd_1, gd_2, \
                         gd_8, gf_3, gf_4, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * gd_8[k]
                  + pb_x[k] * hd_11[k];

        t_26[k] = pa_z[k] * gf_3[k];

        t_27[k] = f_3 * gd_1[k]
                  + pb_z[k] * hd_10[k];

        t_28[k] = pb_y[k] * hd_11[k];

        t_29[k] = f_5 * gd_2[k]
                  + pa_z[k] * gf_4[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_x, pb_y, pb_z, ff0_0, ff1_0, gd_3, \
                         gd_10, gf_5, hd_12, hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_6 * ff0_0[k]
                  - f_7 * ff1_0[k]
                  + pa_y[k] * gf_5[k];

        t_31[k] = f_8 * gd_3[k]
                  + pb_y[k] * hd_12[k];

        t_32[k] = pb_z[k] * hd_12[k];

        t_33[k] = f_5 * gd_10[k]
                  + pb_x[k] * hd_14[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, ff0_3, ff1_3, gd_11, gf_16, \
                         hd_13, hd_14, hd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_z[k] * hd_13[k];

        t_35[k] = f_5 * gd_11[k]
                  + pb_x[k] * hd_15[k];

        t_36[k] = f_9 * ff0_3[k]
                  - f_10 * ff1_3[k]
                  + pa_x[k] * gf_16[k];

        t_37[k] = pb_z[k] * hd_14[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, pb_z, gd_5, gf_6, \
                         gf_9, gf_10, hp0_3, hp1_3, hd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * gd_5[k]
                  + pb_y[k] * hd_15[k];

        t_39[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_15[k];

        t_40[k] = pa_y[k] * gf_9[k];

        t_41[k] = pa_z[k] * gf_6[k];

        t_42[k] = pa_y[k] * gf_10[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pa_y, pa_z, pb_x, pb_z, gd_4, gd_13, \
                         gf_7, gf_8, gf_11, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pa_z[k] * gf_7[k];

        t_44[k] = f_5 * gd_13[k]
                  + pb_x[k] * hd_17[k];

        t_45[k] = pa_y[k] * gf_11[k];

        t_46[k] = pa_z[k] * gf_8[k];

        t_47[k] = f_3 * gd_4[k]
                  + pb_z[k] * hd_16[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pa_z, pb_y, ff0_0, ff1_0, gd_8, gf_9, \
                         gf_12, hd_18, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * gd_8[k]
                  + pb_y[k] * hd_18[k];

        t_49[k] = pa_y[k] * gf_12[k];

        t_50[k] = f_6 * ff0_0[k]
                  - f_7 * ff1_0[k]
                  + pa_z[k] * gf_9[k];

        t_51[k] = pb_y[k] * hd_19[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_x, pb_y, pb_z, gd_6, gd_16, gd_17, hd_19, \
                         hd_20, hd_21, hd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_8 * gd_6[k]
                  + pb_z[k] * hd_19[k];

        t_53[k] = f_5 * gd_16[k]
                  + pb_x[k] * hd_21[k];

        t_54[k] = pb_y[k] * hd_20[k];

        t_55[k] = f_5 * gd_17[k]
                  + pb_x[k] * hd_22[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_y, pb_z, ff0_4, ff1_4, gd_7, gf_22, \
                         hp0_4, hp1_4, hd_21, hd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_21[k];

        t_57[k] = f_8 * gd_7[k]
                  + pb_z[k] * hd_21[k];

        t_58[k] = pb_y[k] * hd_22[k];

        t_59[k] = f_9 * ff0_4[k]
                  - f_10 * ff1_4[k]
                  + pa_x[k] * gf_22[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pb_x, pb_y, pb_z, ff0_1, ff1_1, gd_9, \
                         gd_19, gf_13, hd_23, hd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_9 * ff0_1[k]
                  - f_10 * ff1_1[k]
                  + pa_y[k] * gf_13[k];

        t_61[k] = f_5 * gd_9[k]
                  + pb_y[k] * hd_23[k];

        t_62[k] = pb_z[k] * hd_23[k];

        t_63[k] = f_8 * gd_19[k]
                  + pb_x[k] * hd_25[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, pb_z, ff0_5, ff1_5, gd_20, gf_26, \
                         hd_24, hd_25, hd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_z[k] * hd_24[k];

        t_65[k] = f_8 * gd_20[k]
                  + pb_x[k] * hd_26[k];

        t_66[k] = f_6 * ff0_5[k]
                  - f_7 * ff1_5[k]
                  + pa_x[k] * gf_26[k];

        t_67[k] = pb_z[k] * hd_25[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_z, pb_y, pb_z, gd_9, gd_11, gf_13, \
                         gf_14, hp0_5, hp1_5, hd_26, hd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_5 * gd_11[k]
                  + pb_y[k] * hd_26[k];

        t_69[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_26[k];

        t_70[k] = pa_z[k] * gf_13[k];

        t_71[k] = pa_z[k] * gf_14[k];

        t_72[k] = f_3 * gd_9[k]
                  + pb_z[k] * hd_27[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_z, pb_x, pb_z, gd_10, gd_22, gd_23, \
                         gf_15, gf_16, hd_28, hd_29, hd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pa_z[k] * gf_15[k];

        t_74[k] = f_8 * gd_22[k]
                  + pb_x[k] * hd_29[k];

        t_75[k] = f_8 * gd_23[k]
                  + pb_x[k] * hd_30[k];

        t_76[k] = pa_z[k] * gf_16[k];

        t_77[k] = f_3 * gd_10[k]
                  + pb_z[k] * hd_28[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pa_z, pb_y, gd_11, gd_14, gd_15, \
                         gf_17, gf_18, gf_19, hd_30, hd_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_8 * gd_14[k]
                  + pb_y[k] * hd_30[k];

        t_79[k] = f_5 * gd_11[k]
                  + pa_z[k] * gf_17[k];

        t_80[k] = pa_y[k] * gf_18[k];

        t_81[k] = f_3 * gd_15[k]
                  + pb_y[k] * hd_31[k];

        t_82[k] = pa_y[k] * gf_19[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, pb_x, pb_z, gd_12, gd_16, gd_25, \
                         gd_26, gf_20, gf_21, hd_32, hd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_8 * gd_25[k]
                  + pb_x[k] * hd_32[k];

        t_84[k] = f_8 * gd_26[k]
                  + pb_x[k] * hd_33[k];

        t_85[k] = pa_y[k] * gf_20[k];

        t_86[k] = f_5 * gd_16[k]
                  + pa_y[k] * gf_21[k];

        t_87[k] = f_8 * gd_12[k]
                  + pb_z[k] * hd_32[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pa_z, pb_y, ff0_2, ff1_2, gd_17, gf_18, \
                         gf_22, hd_34, hd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_3 * gd_17[k]
                  + pb_y[k] * hd_34[k];

        t_89[k] = pa_y[k] * gf_22[k];

        t_90[k] = f_9 * ff0_2[k]
                  - f_10 * ff1_2[k]
                  + pa_z[k] * gf_18[k];

        t_91[k] = pb_y[k] * hd_35[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, gd_15, gd_28, gd_29, hd_35, \
                         hd_36, hd_37, hd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_5 * gd_15[k]
                  + pb_z[k] * hd_35[k];

        t_93[k] = f_8 * gd_28[k]
                  + pb_x[k] * hd_37[k];

        t_94[k] = pb_y[k] * hd_36[k];

        t_95[k] = f_8 * gd_29[k]
                  + pb_x[k] * hd_38[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_x, pb_y, pb_z, ff0_8, ff1_8, gd_16, gf_30, \
                         hp0_6, hp1_6, hd_37, hd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_37[k];

        t_97[k] = f_5 * gd_16[k]
                  + pb_z[k] * hd_37[k];

        t_98[k] = pb_y[k] * hd_38[k];

        t_99[k] = f_6 * ff0_8[k]
                  - f_7 * ff1_8[k]
                  + pa_x[k] * gf_30[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pa_x, pb_x, pb_y, pb_z, gd_18, \
                         gd_30, gd_31, gf_31, hd_39, hd_40, hd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_5 * gd_30[k]
                   + pa_x[k] * gf_31[k];

        t_101[k] = f_4 * gd_18[k]
                   + pb_y[k] * hd_39[k];

        t_102[k] = pb_z[k] * hd_39[k];

        t_103[k] = f_3 * gd_31[k]
                   + pb_x[k] * hd_41[k];

        t_104[k] = pb_z[k] * hd_40[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pa_x, pb_x, pb_z, gd_32, gf_33, \
                         gf_34, gf_35, hd_41, hd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_3 * gd_32[k]
                   + pb_x[k] * hd_42[k];

        t_106[k] = pa_x[k] * gf_33[k];

        t_107[k] = pb_z[k] * hd_41[k];

        t_108[k] = pa_x[k] * gf_34[k];

        t_109[k] = pa_x[k] * gf_35[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, pa_z, pb_x, pb_z, gd_18, gd_35, \
                         gf_23, gf_24, gf_25, hd_43, hd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * gf_23[k];

        t_111[k] = pa_z[k] * gf_24[k];

        t_112[k] = f_3 * gd_18[k]
                   + pb_z[k] * hd_43[k];

        t_113[k] = pa_z[k] * gf_25[k];

        t_114[k] = f_3 * gd_35[k]
                   + pb_x[k] * hd_44[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, t_120, pa_x, pb_x, gd_36, gd_37, \
                         gf_36, gf_37, gf_38, gf_39, gf_40, hd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_3 * gd_36[k]
                   + pb_x[k] * hd_45[k];

        t_116[k] = pa_x[k] * gf_36[k];

        t_117[k] = pa_x[k] * gf_37[k];

        t_118[k] = pa_x[k] * gf_38[k];

        t_119[k] = pa_x[k] * gf_39[k];

        t_120[k] = f_5 * gd_37[k]
                   + pa_x[k] * gf_40[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_x, pb_y, pb_z, gd_21, gd_24, gd_38, \
                         gd_39, hd_46, hd_47, hd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_8 * gd_24[k]
                   + pb_y[k] * hd_46[k];

        t_122[k] = f_8 * gd_21[k]
                   + pb_z[k] * hd_46[k];

        t_123[k] = f_3 * gd_38[k]
                   + pb_x[k] * hd_47[k];

        t_124[k] = f_3 * gd_39[k]
                   + pb_x[k] * hd_48[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, t_130, pa_x, pa_y, pb_x, gd_40, \
                         gf_27, gf_41, gf_42, gf_43, gf_44, hd_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_3 * gd_40[k]
                   + pb_x[k] * hd_49[k];

        t_126[k] = pa_x[k] * gf_41[k];

        t_127[k] = pa_x[k] * gf_42[k];

        t_128[k] = pa_x[k] * gf_43[k];

        t_129[k] = pa_x[k] * gf_44[k];

        t_130[k] = pa_y[k] * gf_27[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, pa_y, pb_x, pb_y, gd_27, gd_42, \
                         gd_43, gf_28, gf_29, hd_50, hd_51, hd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_3 * gd_27[k]
                   + pb_y[k] * hd_50[k];

        t_132[k] = pa_y[k] * gf_28[k];

        t_133[k] = f_3 * gd_42[k]
                   + pb_x[k] * hd_51[k];

        t_134[k] = f_3 * gd_43[k]
                   + pb_x[k] * hd_52[k];

        t_135[k] = pa_y[k] * gf_29[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, t_141, pa_x, pb_y, gd_45, gf_45, \
                         gf_46, gf_47, gf_48, gf_49, hd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pa_x[k] * gf_45[k];

        t_137[k] = pa_x[k] * gf_46[k];

        t_138[k] = pa_x[k] * gf_47[k];

        t_139[k] = pa_x[k] * gf_48[k];

        t_140[k] = f_5 * gd_45[k]
                   + pa_x[k] * gf_49[k];

        t_141[k] = pb_y[k] * hd_53[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pb_x, pb_y, pb_z, gd_27, gd_46, gd_47, \
                         hd_53, hd_54, hd_55, hd_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_4 * gd_27[k]
                   + pb_z[k] * hd_53[k];

        t_143[k] = f_3 * gd_46[k]
                   + pb_x[k] * hd_55[k];

        t_144[k] = pb_y[k] * hd_54[k];

        t_145[k] = f_3 * gd_47[k]
                   + pb_x[k] * hd_56[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, pa_x, pb_x, pb_y, gf_51, gf_52, \
                         gf_53, hp0_7, hp1_7, hd_56, hd_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = pa_x[k] * gf_51[k];

        t_147[k] = pa_x[k] * gf_52[k];

        t_148[k] = pb_y[k] * hd_56[k];

        t_149[k] = pa_x[k] * gf_53[k];

        t_150[k] = f_1 * hp0_7[k]
                   - f_2 * hp1_7[k]
                   + pb_x[k] * hd_57[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, pb_x, pb_y, pb_z, gd_30, hd_57, \
                         hd_58, hd_59, hd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_0 * gd_30[k]
                   + pb_y[k] * hd_57[k];

        t_152[k] = pb_z[k] * hd_57[k];

        t_153[k] = pb_x[k] * hd_58[k];

        t_154[k] = pb_x[k] * hd_59[k];

        t_155[k] = pb_x[k] * hd_60[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pb_y, pb_z, gd_31, gd_32, hp0_8, hp0_9, \
                         hp1_8, hp1_9, hd_58, hd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_0 * gd_31[k]
                   + f_1 * hp0_8[k]
                   - f_2 * hp1_8[k]
                   + pb_y[k] * hd_58[k];

        t_157[k] = pb_z[k] * hd_58[k];

        t_158[k] = f_0 * gd_32[k]
                   + pb_y[k] * hd_60[k];

        t_159[k] = f_1 * hp0_9[k]
                   - f_2 * hp1_9[k]
                   + pb_z[k] * hd_60[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, t_165, pa_z, pb_x, pb_z, gd_30, \
                         gf_31, gf_32, hd_61, hd_62, hd_63, hd_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = pa_z[k] * gf_31[k];

        t_161[k] = pa_z[k] * gf_32[k];

        t_162[k] = f_3 * gd_30[k]
                   + pb_z[k] * hd_61[k];

        t_163[k] = pb_x[k] * hd_62[k];

        t_164[k] = pb_x[k] * hd_63[k];

        t_165[k] = pb_x[k] * hd_64[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_z, pb_y, pb_z, gd_31, gd_32, gd_36, \
                         gf_33, gf_35, hd_62, hd_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pa_z[k] * gf_33[k];

        t_167[k] = f_3 * gd_31[k]
                   + pb_z[k] * hd_62[k];

        t_168[k] = f_4 * gd_36[k]
                   + pb_y[k] * hd_64[k];

        t_169[k] = f_5 * gd_32[k]
                   + pa_z[k] * gf_35[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, pb_x, pb_y, pb_z, gd_33, gd_37, \
                         hp0_10, hp1_10, hd_65, hd_66, hd_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_1 * hp0_10[k]
                   - f_2 * hp1_10[k]
                   + pb_x[k] * hd_65[k];

        t_171[k] = f_5 * gd_37[k]
                   + pb_y[k] * hd_65[k];

        t_172[k] = f_8 * gd_33[k]
                   + pb_z[k] * hd_65[k];

        t_173[k] = pb_x[k] * hd_66[k];

        t_174[k] = pb_x[k] * hd_67[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_z, pb_x, pb_y, pb_z, ff0_5, ff1_5, \
                         gd_34, gd_40, gf_36, hd_66, hd_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pb_x[k] * hd_68[k];

        t_176[k] = f_6 * ff0_5[k]
                   - f_7 * ff1_5[k]
                   + pa_z[k] * gf_36[k];

        t_177[k] = f_8 * gd_34[k]
                   + pb_z[k] * hd_66[k];

        t_178[k] = f_5 * gd_40[k]
                   + pb_y[k] * hd_68[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pa_y, pb_x, pb_y, pb_z, ff0_7, ff1_7, \
                         gd_37, gd_41, gf_44, hp0_11, hp1_11, hd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_9 * ff0_7[k]
                   - f_10 * ff1_7[k]
                   + pa_y[k] * gf_44[k];

        t_180[k] = f_1 * hp0_11[k]
                   - f_2 * hp1_11[k]
                   + pb_x[k] * hd_69[k];

        t_181[k] = f_8 * gd_41[k]
                   + pb_y[k] * hd_69[k];

        t_182[k] = f_5 * gd_37[k]
                   + pb_z[k] * hd_69[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, pa_z, pb_x, pb_z, ff0_6, ff1_6, \
                         gd_38, gf_41, hd_70, hd_71, hd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = pb_x[k] * hd_70[k];

        t_184[k] = pb_x[k] * hd_71[k];

        t_185[k] = pb_x[k] * hd_72[k];

        t_186[k] = f_9 * ff0_6[k]
                   - f_10 * ff1_6[k]
                   + pa_z[k] * gf_41[k];

        t_187[k] = f_5 * gd_38[k]
                   + pb_z[k] * hd_70[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, pa_y, pb_y, ff0_8, ff1_8, gd_44, \
                         gd_45, gf_48, gf_49, gf_50, hd_72, hd_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_8 * gd_44[k]
                   + pb_y[k] * hd_72[k];

        t_189[k] = f_6 * ff0_8[k]
                   - f_7 * ff1_8[k]
                   + pa_y[k] * gf_48[k];

        t_190[k] = pa_y[k] * gf_49[k];

        t_191[k] = f_3 * gd_45[k]
                   + pb_y[k] * hd_73[k];

        t_192[k] = pa_y[k] * gf_50[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, pa_y, pb_x, pb_z, gd_42, gd_46, \
                         gf_51, hd_74, hd_75, hd_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pb_x[k] * hd_74[k];

        t_194[k] = pb_x[k] * hd_75[k];

        t_195[k] = pb_x[k] * hd_76[k];

        t_196[k] = f_5 * gd_46[k]
                   + pa_y[k] * gf_51[k];

        t_197[k] = f_4 * gd_42[k]
                   + pb_z[k] * hd_74[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, pa_y, pb_x, pb_y, pb_z, gd_45, \
                         gd_47, gf_53, hp0_12, hp1_12, hd_76, hd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_3 * gd_47[k]
                   + pb_y[k] * hd_76[k];

        t_199[k] = pa_y[k] * gf_53[k];

        t_200[k] = f_1 * hp0_12[k]
                   - f_2 * hp1_12[k]
                   + pb_x[k] * hd_77[k];

        t_201[k] = pb_y[k] * hd_77[k];

        t_202[k] = f_0 * gd_45[k]
                   + pb_z[k] * hd_77[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, t_208, pb_x, pb_y, pb_z, gd_46, \
                         hp0_13, hp1_13, hd_78, hd_79, hd_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = pb_x[k] * hd_78[k];

        t_204[k] = pb_x[k] * hd_79[k];

        t_205[k] = pb_x[k] * hd_80[k];

        t_206[k] = f_1 * hp0_13[k]
                   - f_2 * hp1_13[k]
                   + pb_y[k] * hd_78[k];

        t_207[k] = f_0 * gd_46[k]
                   + pb_z[k] * hd_78[k];

        t_208[k] = pb_y[k] * hd_80[k];
    }

#pragma omp simd aligned(t_209, pb_z, gd_47, hp0_14, hp1_14, hd_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_0 * gd_47[k]
                   + f_1 * hp0_14[k]
                   - f_2 * hp1_14[k]
                   + pb_z[k] * hd_80[k];
    }
}

auto
compute_prim_hf_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ff0, const size_t ff1,
                                     const size_t gd, const size_t gf, const size_t hp0,
                                     const size_t hp1, const size_t hd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);

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

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_32 = buffer.data(ff1 + 32);
    const auto *ff1_40 = buffer.data(ff1 + 40);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
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
    const auto *gd_21 = buffer.data(gd + 21);
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
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
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
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_58 = buffer.data(gf + 58);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, gd_1, gd_2, hp0_0, \
                         hp1_0, hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_0 * gd_1[k]
                 + pb_x[k] * hd_1[k];

        t_4[k] = f_0 * gd_2[k]
                 + pb_x[k] * hd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_y, pb_z, gf_0, hp0_1, hp0_2, hp1_1, \
                         hp1_2, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_6[k] = pb_y[k] * hd_2[k];

        t_7[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_8[k] = pa_y[k] * gf_0[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_y, pb_x, pb_y, gd_0, gd_1, gd_2, gd_4, \
                         gf_3, hd_3, hd_4, hd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * gd_0[k]
                 + pb_y[k] * hd_3[k];

        t_10[k] = f_4 * gd_4[k]
                  + pb_x[k] * hd_4[k];

        t_11[k] = f_5 * gd_1[k]
                  + pa_y[k] * gf_3[k];

        t_12[k] = f_3 * gd_2[k]
                  + pb_y[k] * hd_5[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_y, pa_z, pb_x, pb_z, gd_0, gd_8, \
                         gf_0, gf_3, gf_5, hd_6, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * gf_5[k];

        t_14[k] = pa_z[k] * gf_0[k];

        t_15[k] = f_3 * gd_0[k]
                  + pb_z[k] * hd_6[k];

        t_16[k] = f_4 * gd_8[k]
                  + pb_x[k] * hd_8[k];

        t_17[k] = pa_z[k] * gf_3[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pa_z, pb_z, ff0_0, ff1_0, gd_1, gd_2, gf_5, \
                         gf_6, hd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * gd_1[k]
                  + pb_z[k] * hd_7[k];

        t_19[k] = f_5 * gd_2[k]
                  + pa_z[k] * gf_5[k];

        t_20[k] = f_6 * ff0_0[k]
                  - f_7 * ff1_0[k]
                  + pa_y[k] * gf_6[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_x, pb_x, pb_y, pb_z, ff0_3, ff1_13, \
                         gd_3, gd_10, gf_14, hd_9, hd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * gd_3[k]
                  + pb_y[k] * hd_9[k];

        t_22[k] = pb_z[k] * hd_9[k];

        t_23[k] = f_5 * gd_10[k]
                  + pb_x[k] * hd_10[k];

        t_24[k] = f_9 * ff0_3[k]
                  - f_10 * ff1_13[k]
                  + pa_x[k] * gf_14[k];

        t_25[k] = pb_z[k] * hd_10[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, pb_y, pb_z, gd_5, gf_7, gf_9, \
                         hp0_3, hp1_3, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * gd_5[k]
                  + pb_y[k] * hd_11[k];

        t_27[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_11[k];

        t_28[k] = pa_y[k] * gf_9[k];

        t_29[k] = pa_z[k] * gf_7[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pa_z, pb_y, pb_z, ff0_0, ff1_0, gd_4, \
                         gd_8, gf_8, gf_10, hd_12, hd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * gd_4[k]
                  + pb_z[k] * hd_12[k];

        t_31[k] = f_3 * gd_8[k]
                  + pb_y[k] * hd_13[k];

        t_32[k] = pa_y[k] * gf_10[k];

        t_33[k] = f_6 * ff0_0[k]
                  - f_7 * ff1_0[k]
                  + pa_z[k] * gf_8[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, gd_6, gd_7, gd_16, \
                         hp0_4, hp1_4, hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * hd_14[k];

        t_35[k] = f_8 * gd_6[k]
                  + pb_z[k] * hd_14[k];

        t_36[k] = f_5 * gd_16[k]
                  + pb_x[k] * hd_16[k];

        t_37[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_15[k];

        t_38[k] = f_8 * gd_7[k]
                  + pb_z[k] * hd_15[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pa_y, pb_y, ff0_1, ff0_4, ff1_6, \
                         ff1_17, gd_9, gf_11, gf_23, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pb_y[k] * hd_16[k];

        t_40[k] = f_9 * ff0_4[k]
                  - f_10 * ff1_17[k]
                  + pa_x[k] * gf_23[k];

        t_41[k] = f_9 * ff0_1[k]
                  - f_10 * ff1_6[k]
                  + pa_y[k] * gf_11[k];

        t_42[k] = f_5 * gd_9[k]
                  + pb_y[k] * hd_17[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pb_x, pb_z, ff0_5, ff1_21, gd_18, \
                         gf_25, hd_17, hd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_z[k] * hd_17[k];

        t_44[k] = f_8 * gd_18[k]
                  + pb_x[k] * hd_18[k];

        t_45[k] = f_6 * ff0_5[k]
                  - f_7 * ff1_21[k]
                  + pa_x[k] * gf_25[k];

        t_46[k] = pb_z[k] * hd_18[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pa_z, pb_y, pb_z, gd_9, gd_11, gf_11, \
                         gf_14, hp0_5, hp1_5, hd_19, hd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_5 * gd_11[k]
                  + pb_y[k] * hd_19[k];

        t_48[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_19[k];

        t_49[k] = pa_z[k] * gf_11[k];

        t_50[k] = f_3 * gd_9[k]
                  + pb_z[k] * hd_20[k];

        t_51[k] = pa_z[k] * gf_14[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_y, pa_z, pb_y, pb_z, gd_10, gd_11, gd_13, \
                         gf_16, gf_17, hd_21, hd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * gd_10[k]
                  + pb_z[k] * hd_21[k];

        t_53[k] = f_8 * gd_13[k]
                  + pb_y[k] * hd_22[k];

        t_54[k] = f_5 * gd_11[k]
                  + pa_z[k] * gf_16[k];

        t_55[k] = pa_y[k] * gf_17[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_y, pb_y, pb_z, gd_12, gd_15, gd_16, \
                         gf_19, gf_21, gf_23, hd_23, hd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_y[k] * gf_19[k];

        t_57[k] = f_5 * gd_15[k]
                  + pa_y[k] * gf_21[k];

        t_58[k] = f_8 * gd_12[k]
                  + pb_z[k] * hd_23[k];

        t_59[k] = f_3 * gd_16[k]
                  + pb_y[k] * hd_24[k];

        t_60[k] = pa_y[k] * gf_23[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_z, pb_x, pb_y, pb_z, ff0_2, ff1_8, gd_14, \
                         gd_21, gf_17, hd_25, hd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_9 * ff0_2[k]
                  - f_10 * ff1_8[k]
                  + pa_z[k] * gf_17[k];

        t_62[k] = pb_y[k] * hd_25[k];

        t_63[k] = f_5 * gd_14[k]
                  + pb_z[k] * hd_25[k];

        t_64[k] = f_8 * gd_21[k]
                  + pb_x[k] * hd_27[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pb_y, pb_z, ff0_8, ff1_40, gd_15, \
                         gf_28, hp0_6, hp1_6, hd_26, hd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_26[k];

        t_66[k] = f_5 * gd_15[k]
                  + pb_z[k] * hd_26[k];

        t_67[k] = pb_y[k] * hd_27[k];

        t_68[k] = f_6 * ff0_8[k]
                  - f_7 * ff1_40[k]
                  + pa_x[k] * gf_28[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pa_x, pb_x, pb_y, gd_17, gd_22, gd_23, \
                         gf_29, gf_32, gf_34, hd_28, hd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_5 * gd_22[k]
                  + pa_x[k] * gf_29[k];

        t_70[k] = f_4 * gd_17[k]
                  + pb_y[k] * hd_28[k];

        t_71[k] = f_3 * gd_23[k]
                  + pb_x[k] * hd_29[k];

        t_72[k] = pa_x[k] * gf_32[k];

        t_73[k] = pa_x[k] * gf_34[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, pa_x, pa_z, pb_z, gd_17, gf_24, \
                         gf_35, gf_37, gf_38, gf_39, hd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pa_x[k] * gf_35[k];

        t_75[k] = pa_z[k] * gf_24[k];

        t_76[k] = f_3 * gd_17[k]
                  + pb_z[k] * hd_30[k];

        t_77[k] = pa_x[k] * gf_37[k];

        t_78[k] = pa_x[k] * gf_38[k];

        t_79[k] = pa_x[k] * gf_39[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, t_85, pa_x, pb_z, gd_19, gd_28, gf_40, \
                         gf_43, gf_44, gf_45, gf_46, hd_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_5 * gd_28[k]
                  + pa_x[k] * gf_40[k];

        t_81[k] = f_8 * gd_19[k]
                  + pb_z[k] * hd_31[k];

        t_82[k] = pa_x[k] * gf_43[k];

        t_83[k] = pa_x[k] * gf_44[k];

        t_84[k] = pa_x[k] * gf_45[k];

        t_85[k] = pa_x[k] * gf_46[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, t_91, pa_x, pa_y, gd_33, gf_26, gf_27, \
                         gf_47, gf_48, gf_49, gf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pa_y[k] * gf_26[k];

        t_87[k] = pa_y[k] * gf_27[k];

        t_88[k] = pa_x[k] * gf_47[k];

        t_89[k] = pa_x[k] * gf_48[k];

        t_90[k] = pa_x[k] * gf_49[k];

        t_91[k] = f_5 * gd_33[k]
                  + pa_x[k] * gf_51[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, pa_x, pb_x, pb_z, gd_20, gd_35, gf_55, \
                         gf_56, gf_58, hd_32, hd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_4 * gd_20[k]
                  + pb_z[k] * hd_32[k];

        t_93[k] = f_3 * gd_35[k]
                  + pb_x[k] * hd_33[k];

        t_94[k] = pa_x[k] * gf_55[k];

        t_95[k] = pa_x[k] * gf_56[k];

        t_96[k] = pa_x[k] * gf_58[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, pb_x, pb_y, gd_22, gd_23, hp0_7, \
                         hp0_8, hp1_7, hp1_8, hd_34, hd_35, hd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_34[k];

        t_98[k] = f_0 * gd_22[k]
                  + pb_y[k] * hd_34[k];

        t_99[k] = pb_x[k] * hd_35[k];

        t_100[k] = pb_x[k] * hd_36[k];

        t_101[k] = f_0 * gd_23[k]
                   + f_1 * hp0_8[k]
                   - f_2 * hp1_8[k]
                   + pb_y[k] * hd_35[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, pa_z, pb_y, pb_z, gd_22, gd_24, \
                         gf_29, hp0_9, hp1_9, hd_35, hd_36, hd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = pb_z[k] * hd_35[k];

        t_103[k] = f_0 * gd_24[k]
                   + pb_y[k] * hd_36[k];

        t_104[k] = f_1 * hp0_9[k]
                   - f_2 * hp1_9[k]
                   + pb_z[k] * hd_36[k];

        t_105[k] = pa_z[k] * gf_29[k];

        t_106[k] = f_3 * gd_22[k]
                   + pb_z[k] * hd_37[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_z, pb_y, pb_z, gd_23, gd_24, gd_27, \
                         gf_32, gf_35, hd_38, hd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pa_z[k] * gf_32[k];

        t_108[k] = f_3 * gd_23[k]
                   + pb_z[k] * hd_38[k];

        t_109[k] = f_4 * gd_27[k]
                   + pb_y[k] * hd_39[k];

        t_110[k] = f_5 * gd_24[k]
                   + pa_z[k] * gf_35[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pb_x, pb_z, gd_25, hp0_10, hp1_10, hd_40, \
                         hd_41, hd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_1 * hp0_10[k]
                   - f_2 * hp1_10[k]
                   + pb_x[k] * hd_40[k];

        t_112[k] = f_8 * gd_25[k]
                   + pb_z[k] * hd_40[k];

        t_113[k] = pb_x[k] * hd_41[k];

        t_114[k] = pb_x[k] * hd_42[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pa_z, pb_y, pb_z, ff0_5, ff1_21, gd_26, gd_30, \
                         gf_36, hd_41, hd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_6 * ff0_5[k]
                   - f_7 * ff1_21[k]
                   + pa_z[k] * gf_36[k];

        t_116[k] = f_8 * gd_26[k]
                   + pb_z[k] * hd_41[k];

        t_117[k] = f_5 * gd_30[k]
                   + pb_y[k] * hd_42[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pb_x, pb_z, ff0_7, ff1_32, gd_28, \
                         gf_46, hp0_11, hp1_11, hd_43, hd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_9 * ff0_7[k]
                   - f_10 * ff1_32[k]
                   + pa_y[k] * gf_46[k];

        t_119[k] = f_1 * hp0_11[k]
                   - f_2 * hp1_11[k]
                   + pb_x[k] * hd_43[k];

        t_120[k] = f_5 * gd_28[k]
                   + pb_z[k] * hd_43[k];

        t_121[k] = pb_x[k] * hd_44[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pa_z, pb_x, pb_y, pb_z, ff0_6, ff1_25, \
                         gd_29, gd_32, gf_43, hd_44, hd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = pb_x[k] * hd_45[k];

        t_123[k] = f_9 * ff0_6[k]
                   - f_10 * ff1_25[k]
                   + pa_z[k] * gf_43[k];

        t_124[k] = f_5 * gd_29[k]
                   + pb_z[k] * hd_44[k];

        t_125[k] = f_8 * gd_32[k]
                   + pb_y[k] * hd_45[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, pa_y, pb_z, ff0_8, ff1_40, gd_31, \
                         gd_34, gf_50, gf_51, gf_52, gf_55, hd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_6 * ff0_8[k]
                   - f_7 * ff1_40[k]
                   + pa_y[k] * gf_50[k];

        t_127[k] = pa_y[k] * gf_51[k];

        t_128[k] = pa_y[k] * gf_52[k];

        t_129[k] = f_5 * gd_34[k]
                   + pa_y[k] * gf_55[k];

        t_130[k] = f_4 * gd_31[k]
                   + pb_z[k] * hd_46[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_y, pb_x, pb_y, pb_z, gd_33, gd_35, \
                         gf_58, hp0_12, hp1_12, hd_47, hd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_3 * gd_35[k]
                   + pb_y[k] * hd_47[k];

        t_132[k] = pa_y[k] * gf_58[k];

        t_133[k] = f_1 * hp0_12[k]
                   - f_2 * hp1_12[k]
                   + pb_x[k] * hd_48[k];

        t_134[k] = f_0 * gd_33[k]
                   + pb_z[k] * hd_48[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, pb_x, pb_y, pb_z, gd_34, hp0_13, \
                         hp1_13, hd_49, hd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = pb_x[k] * hd_49[k];

        t_136[k] = pb_x[k] * hd_50[k];

        t_137[k] = f_1 * hp0_13[k]
                   - f_2 * hp1_13[k]
                   + pb_y[k] * hd_49[k];

        t_138[k] = f_0 * gd_34[k]
                   + pb_z[k] * hd_49[k];

        t_139[k] = pb_y[k] * hd_50[k];
    }

#pragma omp simd aligned(t_140, pb_z, gd_35, hp0_14, hp1_14, hd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * gd_35[k]
                   + f_1 * hp0_14[k]
                   - f_2 * hp1_14[k]
                   + pb_z[k] * hd_50[k];
    }
}

auto
compute_prim_hf_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ff0, const size_t ff1,
                                     const size_t gd, const size_t gf, const size_t hp0,
                                     const size_t hp1, const size_t hd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / p;
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 1.0 / p;

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

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_8 = buffer.data(ff1 + 8);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_20 = buffer.data(gd + 20);

    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, hp0_0, hp0_1, hp1_0, \
                         hp1_1, hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_4[k] = pb_y[k] * hd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, ff0_0, ff1_0, gd_4, gf_6, \
                         hp0_2, hp1_2, hd_2, hd_3, hd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_6[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pa_y[k] * gf_6[k];

        t_7[k] = pb_z[k] * hd_3[k];

        t_8[k] = f_5 * gd_4[k]
                 + pb_x[k] * hd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, ff0_3, ff1_3, gf_11, hp0_3, hp1_3, hd_4, \
                         hd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * ff0_3[k]
                 - f_7 * ff1_3[k]
                 + pa_x[k] * gf_11[k];

        t_10[k] = pb_z[k] * hd_4[k];

        t_11[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, ff0_0, ff1_0, gd_8, gf_7, \
                         hp0_4, hp1_4, hd_6, hd_7, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * ff0_0[k]
                  - f_4 * ff1_0[k]
                  + pa_z[k] * gf_7[k];

        t_13[k] = pb_y[k] * hd_6[k];

        t_14[k] = f_5 * gd_8[k]
                  + pb_x[k] * hd_8[k];

        t_15[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, ff0_1, ff0_4, ff1_1, \
                         ff1_4, gf_8, gf_19, hd_8, hd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * hd_8[k];

        t_17[k] = f_6 * ff0_4[k]
                  - f_7 * ff1_4[k]
                  + pa_x[k] * gf_19[k];

        t_18[k] = f_6 * ff0_1[k]
                  - f_7 * ff1_1[k]
                  + pa_y[k] * gf_8[k];

        t_19[k] = pb_z[k] * hd_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, ff0_5, ff1_5, gd_9, gf_20, \
                         hp0_5, hp1_5, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_8 * gd_9[k]
                  + pb_x[k] * hd_10[k];

        t_21[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pa_x[k] * gf_20[k];

        t_22[k] = pb_z[k] * hd_10[k];

        t_23[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, ff0_2, ff1_2, gd_10, gf_14, \
                         hp0_6, hp1_6, hd_12, hd_13, hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_6 * ff0_2[k]
                  - f_7 * ff1_2[k]
                  + pa_z[k] * gf_14[k];

        t_25[k] = pb_y[k] * hd_12[k];

        t_26[k] = f_8 * gd_10[k]
                  + pb_x[k] * hd_14[k];

        t_27[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_x, pb_y, ff0_8, ff1_8, gf_21, hp0_7, \
                         hp1_7, hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * hd_14[k];

        t_29[k] = f_3 * ff0_8[k]
                  - f_4 * ff1_8[k]
                  + pa_x[k] * gf_21[k];

        t_30[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_15[k];

        t_31[k] = pb_x[k] * hd_16[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_x, pb_y, pb_z, gd_12, hp0_8, hp0_9, hp1_8, \
                         hp1_9, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_x[k] * hd_17[k];

        t_33[k] = f_0 * gd_12[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_16[k];

        t_34[k] = pb_z[k] * hd_16[k];

        t_35[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_9[k]
                  + pb_z[k] * hd_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, ff0_5, ff1_5, gf_28, hp0_10, \
                         hp1_10, hd_18, hd_19, hd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * hp0_10[k]
                  - f_2 * hp1_10[k]
                  + pb_x[k] * hd_18[k];

        t_37[k] = pb_x[k] * hd_19[k];

        t_38[k] = pb_x[k] * hd_20[k];

        t_39[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_5[k]
                  + pa_z[k] * gf_28[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pb_x, pb_y, ff0_7, ff1_7, gd_16, gf_34, \
                         hp0_11, hp1_11, hd_20, hd_21, hd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * gd_16[k]
                  + pb_y[k] * hd_20[k];

        t_41[k] = f_6 * ff0_7[k]
                  - f_7 * ff1_7[k]
                  + pa_y[k] * gf_34[k];

        t_42[k] = f_1 * hp0_11[k]
                  - f_2 * hp1_11[k]
                  + pb_x[k] * hd_21[k];

        t_43[k] = pb_x[k] * hd_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_y, pa_z, pb_x, pb_y, ff0_6, ff0_8, ff1_6, \
                         ff1_8, gd_17, gf_32, gf_35, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * hd_23[k];

        t_45[k] = f_6 * ff0_6[k]
                  - f_7 * ff1_6[k]
                  + pa_z[k] * gf_32[k];

        t_46[k] = f_8 * gd_17[k]
                  + pb_y[k] * hd_23[k];

        t_47[k] = f_3 * ff0_8[k]
                  - f_4 * ff1_8[k]
                  + pa_y[k] * gf_35[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pb_x, pb_y, hp0_12, hp0_13, hp1_12, \
                         hp1_13, hd_24, hd_25, hd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * hp0_12[k]
                  - f_2 * hp1_12[k]
                  + pb_x[k] * hd_24[k];

        t_49[k] = pb_x[k] * hd_25[k];

        t_50[k] = pb_x[k] * hd_26[k];

        t_51[k] = f_1 * hp0_13[k]
                  - f_2 * hp1_13[k]
                  + pb_y[k] * hd_25[k];

        t_52[k] = pb_y[k] * hd_26[k];
    }

#pragma omp simd aligned(t_53, pb_z, gd_20, hp0_14, hp1_14, hd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * gd_20[k]
                  + f_1 * hp0_14[k]
                  - f_2 * hp1_14[k]
                  + pb_z[k] * hd_26[k];
    }
}

auto
compute_prim_hf_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ff0, const size_t ff1,
                                     const size_t gd, const size_t gf, const size_t hp0,
                                     const size_t hp1, const size_t hd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / p;
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 1.0 / p;

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

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_10 = buffer.data(ff1 + 10);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_19 = buffer.data(ff1 + 19);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_29 = buffer.data(ff1 + 29);
    const auto *ff1_34 = buffer.data(ff1 + 34);
    const auto *ff1_41 = buffer.data(ff1 + 41);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_55 = buffer.data(gf + 55);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, hp0_0, hp0_1, hp1_0, \
                         hp1_1, hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_4[k] = pb_y[k] * hd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, ff0_0, ff1_0, gd_6, gf_7, \
                         hp0_2, hp1_2, hd_2, hd_3, hd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_6[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pa_y[k] * gf_7[k];

        t_7[k] = pb_z[k] * hd_3[k];

        t_8[k] = f_5 * gd_6[k]
                 + pb_x[k] * hd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, ff0_3, ff1_16, gf_17, hp0_3, hp1_3, \
                         hd_4, hd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * ff0_3[k]
                 - f_7 * ff1_16[k]
                 + pa_x[k] * gf_17[k];

        t_10[k] = pb_z[k] * hd_4[k];

        t_11[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, ff0_0, ff1_0, gd_10, gf_10, \
                         hp0_4, hp1_4, hd_6, hd_7, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * ff0_0[k]
                  - f_4 * ff1_0[k]
                  + pa_z[k] * gf_10[k];

        t_13[k] = pb_y[k] * hd_6[k];

        t_14[k] = f_5 * gd_10[k]
                  + pb_x[k] * hd_8[k];

        t_15[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, ff0_1, ff0_4, ff1_7, \
                         ff1_19, gf_14, gf_27, hd_8, hd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * hd_8[k];

        t_17[k] = f_6 * ff0_4[k]
                  - f_7 * ff1_19[k]
                  + pa_x[k] * gf_27[k];

        t_18[k] = f_6 * ff0_1[k]
                  - f_7 * ff1_7[k]
                  + pa_y[k] * gf_14[k];

        t_19[k] = pb_z[k] * hd_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, ff0_5, ff1_24, gd_11, \
                         gf_30, hp0_5, hp1_5, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_8 * gd_11[k]
                  + pb_x[k] * hd_10[k];

        t_21[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_24[k]
                  + pa_x[k] * gf_30[k];

        t_22[k] = pb_z[k] * hd_10[k];

        t_23[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, ff0_2, ff1_10, gd_12, \
                         gf_22, hp0_6, hp1_6, hd_12, hd_13, hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_6 * ff0_2[k]
                  - f_7 * ff1_10[k]
                  + pa_z[k] * gf_22[k];

        t_25[k] = pb_y[k] * hd_12[k];

        t_26[k] = f_8 * gd_12[k]
                  + pb_x[k] * hd_14[k];

        t_27[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_x, pb_y, ff0_8, ff1_41, gf_34, \
                         hp0_7, hp1_7, hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * hd_14[k];

        t_29[k] = f_3 * ff0_8[k]
                  - f_4 * ff1_41[k]
                  + pa_x[k] * gf_34[k];

        t_30[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_15[k];

        t_31[k] = pb_x[k] * hd_16[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_x, pb_y, pb_z, gd_14, hp0_8, hp0_9, hp1_8, \
                         hp1_9, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_x[k] * hd_17[k];

        t_33[k] = f_0 * gd_14[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_16[k];

        t_34[k] = pb_z[k] * hd_16[k];

        t_35[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_9[k]
                  + pb_z[k] * hd_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, ff0_5, ff1_24, gf_44, hp0_10, \
                         hp1_10, hd_18, hd_19, hd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * hp0_10[k]
                  - f_2 * hp1_10[k]
                  + pb_x[k] * hd_18[k];

        t_37[k] = pb_x[k] * hd_19[k];

        t_38[k] = pb_x[k] * hd_20[k];

        t_39[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_24[k]
                  + pa_z[k] * gf_44[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pb_x, pb_y, ff0_7, ff1_34, gd_19, \
                         gf_51, hp0_11, hp1_11, hd_20, hd_21, hd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * gd_19[k]
                  + pb_y[k] * hd_20[k];

        t_41[k] = f_6 * ff0_7[k]
                  - f_7 * ff1_34[k]
                  + pa_y[k] * gf_51[k];

        t_42[k] = f_1 * hp0_11[k]
                  - f_2 * hp1_11[k]
                  + pb_x[k] * hd_21[k];

        t_43[k] = pb_x[k] * hd_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_y, pa_z, pb_x, pb_y, ff0_6, ff0_8, ff1_29, \
                         ff1_41, gd_20, gf_49, gf_55, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * hd_23[k];

        t_45[k] = f_6 * ff0_6[k]
                  - f_7 * ff1_29[k]
                  + pa_z[k] * gf_49[k];

        t_46[k] = f_8 * gd_20[k]
                  + pb_y[k] * hd_23[k];

        t_47[k] = f_3 * ff0_8[k]
                  - f_4 * ff1_41[k]
                  + pa_y[k] * gf_55[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pb_x, pb_y, hp0_12, hp0_13, hp1_12, \
                         hp1_13, hd_24, hd_25, hd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * hp0_12[k]
                  - f_2 * hp1_12[k]
                  + pb_x[k] * hd_24[k];

        t_49[k] = pb_x[k] * hd_25[k];

        t_50[k] = pb_x[k] * hd_26[k];

        t_51[k] = f_1 * hp0_13[k]
                  - f_2 * hp1_13[k]
                  + pb_y[k] * hd_25[k];

        t_52[k] = pb_y[k] * hd_26[k];
    }

#pragma omp simd aligned(t_53, pb_z, gd_23, hp0_14, hp1_14, hd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * gd_23[k]
                  + f_1 * hp0_14[k]
                  - f_2 * hp1_14[k]
                  + pb_z[k] * hd_26[k];
    }
}

auto
compute_prim_hf_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ff0, const size_t ff1,
                                     const size_t gd, const size_t gf, const size_t hp0,
                                     const size_t hp1, const size_t hd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 1.0 / p;

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

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_16 = buffer.data(ff0 + 16);
    const auto *ff0_19 = buffer.data(ff0 + 19);
    const auto *ff0_24 = buffer.data(ff0 + 24);
    const auto *ff0_29 = buffer.data(ff0 + 29);
    const auto *ff0_34 = buffer.data(ff0 + 34);
    const auto *ff0_41 = buffer.data(ff0 + 41);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_9 = buffer.data(ff1 + 9);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_16 = buffer.data(ff1 + 16);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_24 = buffer.data(ff1 + 24);
    const auto *ff1_28 = buffer.data(ff1 + 28);
    const auto *ff1_35 = buffer.data(ff1 + 35);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_47 = buffer.data(gf + 47);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, hp0_0, hp0_1, hp1_0, \
                         hp1_1, hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_4[k] = pb_y[k] * hd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, t_10, pa_y, pa_z, pb_z, gd_1, gf_0, gf_3, \
                         gf_5, hp0_2, hp1_2, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_6[k] = pa_y[k] * gf_0[k];

        t_7[k] = f_3 * gd_1[k]
                 + pa_y[k] * gf_3[k];

        t_8[k] = pa_y[k] * gf_5[k];

        t_9[k] = pa_z[k] * gf_0[k];

        t_10[k] = pa_z[k] * gf_3[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pa_y, pa_z, pb_x, pb_z, ff0_0, ff1_0, gd_2, \
                         gd_7, gf_5, gf_6, hd_3, hd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * gd_2[k]
                  + pa_z[k] * gf_5[k];

        t_12[k] = f_4 * ff0_0[k]
                  - f_5 * ff1_0[k]
                  + pa_y[k] * gf_6[k];

        t_13[k] = pb_z[k] * hd_3[k];

        t_14[k] = f_3 * gd_7[k]
                  + pb_x[k] * hd_4[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pa_z, pb_z, ff0_16, ff1_13, gf_7, \
                         gf_13, hp0_3, hp1_3, hd_4, hd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * ff0_16[k]
                  - f_7 * ff1_13[k]
                  + pa_x[k] * gf_13[k];

        t_16[k] = pb_z[k] * hd_4[k];

        t_17[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_5[k];

        t_18[k] = pa_z[k] * gf_7[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_x, pb_y, ff0_0, ff1_0, gd_11, \
                         gf_8, gf_9, hd_6, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * gf_9[k];

        t_20[k] = f_4 * ff0_0[k]
                  - f_5 * ff1_0[k]
                  + pa_z[k] * gf_8[k];

        t_21[k] = pb_y[k] * hd_6[k];

        t_22[k] = f_3 * gd_11[k]
                  + pb_x[k] * hd_8[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pb_y, ff0_19, ff1_16, gf_21, hp0_4, hp1_4, \
                         hd_7, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_7[k];

        t_24[k] = pb_y[k] * hd_8[k];

        t_25[k] = f_6 * ff0_19[k]
                  - f_7 * ff1_16[k]
                  + pa_x[k] * gf_21[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_x, pb_z, ff0_7, ff1_7, gd_12, gf_10, hd_9, \
                         hd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_6 * ff0_7[k]
                  - f_7 * ff1_7[k]
                  + pa_y[k] * gf_10[k];

        t_27[k] = pb_z[k] * hd_9[k];

        t_28[k] = f_8 * gd_12[k]
                  + pb_x[k] * hd_10[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_x, pa_z, pb_z, ff0_24, ff1_21, gf_10, \
                         gf_23, hp0_5, hp1_5, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_4 * ff0_24[k]
                  - f_5 * ff1_21[k]
                  + pa_x[k] * gf_23[k];

        t_30[k] = pb_z[k] * hd_10[k];

        t_31[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_11[k];

        t_32[k] = pa_z[k] * gf_10[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pa_z, ff0_10, ff1_9, gd_8, gd_10, \
                         gf_13, gf_15, gf_16, gf_19, gf_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * gf_13[k];

        t_34[k] = f_3 * gd_8[k]
                  + pa_z[k] * gf_15[k];

        t_35[k] = f_3 * gd_10[k]
                  + pa_y[k] * gf_19[k];

        t_36[k] = pa_y[k] * gf_21[k];

        t_37[k] = f_6 * ff0_10[k]
                  - f_7 * ff1_9[k]
                  + pa_z[k] * gf_16[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_y, gd_13, hp0_6, hp1_6, hd_12, \
                         hd_13, hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * hd_12[k];

        t_39[k] = f_8 * gd_13[k]
                  + pb_x[k] * hd_14[k];

        t_40[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_13[k];

        t_41[k] = pb_y[k] * hd_14[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_x, pa_z, ff0_41, ff1_35, gd_14, \
                         gd_19, gf_22, gf_25, gf_26, gf_29, gf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_4 * ff0_41[k]
                  - f_5 * ff1_35[k]
                  + pa_x[k] * gf_25[k];

        t_43[k] = f_3 * gd_14[k]
                  + pa_x[k] * gf_26[k];

        t_44[k] = pa_x[k] * gf_29[k];

        t_45[k] = pa_z[k] * gf_22[k];

        t_46[k] = f_3 * gd_19[k]
                  + pa_x[k] * gf_34[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pa_x, pb_x, gd_24, gf_42, gf_47, hp0_7, \
                         hp1_7, hd_15, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * gd_24[k]
                  + pa_x[k] * gf_42[k];

        t_48[k] = pa_x[k] * gf_47[k];

        t_49[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_15[k];

        t_50[k] = pb_x[k] * hd_16[k];

        t_51[k] = pb_x[k] * hd_17[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_z, pb_y, pb_z, gd_15, gf_26, hp0_8, hp0_9, \
                         hp1_8, hp1_9, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_0 * gd_15[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_16[k];

        t_53[k] = pb_z[k] * hd_16[k];

        t_54[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_9[k]
                  + pb_z[k] * hd_17[k];

        t_55[k] = pa_z[k] * gf_26[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_z, pb_x, gd_16, gf_29, gf_31, \
                         hp0_10, hp1_10, hd_18, hd_19, hd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_z[k] * gf_29[k];

        t_57[k] = f_3 * gd_16[k]
                  + pa_z[k] * gf_31[k];

        t_58[k] = f_1 * hp0_10[k]
                  - f_2 * hp1_10[k]
                  + pb_x[k] * hd_18[k];

        t_59[k] = pb_x[k] * hd_19[k];

        t_60[k] = pb_x[k] * hd_20[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pa_y, pa_z, pb_y, ff0_24, ff0_34, ff1_21, ff1_28, \
                         gd_21, gf_32, gf_39, hd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_4 * ff0_24[k]
                  - f_5 * ff1_21[k]
                  + pa_z[k] * gf_32[k];

        t_62[k] = f_3 * gd_21[k]
                  + pb_y[k] * hd_20[k];

        t_63[k] = f_6 * ff0_34[k]
                  - f_7 * ff1_28[k]
                  + pa_y[k] * gf_39[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_z, pb_x, ff0_29, ff1_24, gf_37, hp0_11, \
                         hp1_11, hd_21, hd_22, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_1 * hp0_11[k]
                  - f_2 * hp1_11[k]
                  + pb_x[k] * hd_21[k];

        t_65[k] = pb_x[k] * hd_22[k];

        t_66[k] = pb_x[k] * hd_23[k];

        t_67[k] = f_6 * ff0_29[k]
                  - f_7 * ff1_24[k]
                  + pa_z[k] * gf_37[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_y, pb_y, ff0_41, ff1_35, gd_23, gd_25, \
                         gf_41, gf_45, gf_47, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_8 * gd_23[k]
                  + pb_y[k] * hd_23[k];

        t_69[k] = f_4 * ff0_41[k]
                  - f_5 * ff1_35[k]
                  + pa_y[k] * gf_41[k];

        t_70[k] = f_3 * gd_25[k]
                  + pa_y[k] * gf_45[k];

        t_71[k] = pa_y[k] * gf_47[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pb_x, pb_y, hp0_12, hp0_13, hp1_12, \
                         hp1_13, hd_24, hd_25, hd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * hp0_12[k]
                  - f_2 * hp1_12[k]
                  + pb_x[k] * hd_24[k];

        t_73[k] = pb_x[k] * hd_25[k];

        t_74[k] = pb_x[k] * hd_26[k];

        t_75[k] = f_1 * hp0_13[k]
                  - f_2 * hp1_13[k]
                  + pb_y[k] * hd_25[k];

        t_76[k] = pb_y[k] * hd_26[k];
    }

#pragma omp simd aligned(t_77, pb_z, gd_26, hp0_14, hp1_14, hd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_0 * gd_26[k]
                  + f_1 * hp0_14[k]
                  - f_2 * hp1_14[k]
                  + pb_z[k] * hd_26[k];
    }
}

auto
compute_prim_hf_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ff0, const size_t ff1,
                                     const size_t gd, const size_t gf, const size_t hp0,
                                     const size_t hp1, const size_t hd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / p;
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 1.0 / p;

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

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_9 = buffer.data(ff1 + 9);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_18 = buffer.data(ff1 + 18);
    const auto *ff1_20 = buffer.data(ff1 + 20);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_38 = buffer.data(gf + 38);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, hp0_0, hp0_1, hp1_0, \
                         hp1_1, hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_4[k] = pb_y[k] * hd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, ff0_0, ff1_0, gd_6, gf_6, \
                         hp0_2, hp1_2, hd_2, hd_3, hd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_6[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pa_y[k] * gf_6[k];

        t_7[k] = pb_z[k] * hd_3[k];

        t_8[k] = f_5 * gd_6[k]
                 + pb_x[k] * hd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, ff0_3, ff1_9, gf_11, hp0_3, hp1_3, hd_4, \
                         hd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * ff0_3[k]
                 - f_7 * ff1_9[k]
                 + pa_x[k] * gf_11[k];

        t_10[k] = pb_z[k] * hd_4[k];

        t_11[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, ff0_0, ff1_0, gd_10, gf_7, \
                         hp0_4, hp1_4, hd_6, hd_7, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * ff0_0[k]
                  - f_4 * ff1_0[k]
                  + pa_z[k] * gf_7[k];

        t_13[k] = pb_y[k] * hd_6[k];

        t_14[k] = f_5 * gd_10[k]
                  + pb_x[k] * hd_8[k];

        t_15[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, ff0_1, ff0_4, ff1_6, \
                         ff1_11, gf_8, gf_19, hd_8, hd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * hd_8[k];

        t_17[k] = f_6 * ff0_4[k]
                  - f_7 * ff1_11[k]
                  + pa_x[k] * gf_19[k];

        t_18[k] = f_6 * ff0_1[k]
                  - f_7 * ff1_6[k]
                  + pa_y[k] * gf_8[k];

        t_19[k] = pb_z[k] * hd_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, ff0_5, ff1_15, gd_11, \
                         gf_21, hp0_5, hp1_5, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_8 * gd_11[k]
                  + pb_x[k] * hd_10[k];

        t_21[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_15[k]
                  + pa_x[k] * gf_21[k];

        t_22[k] = pb_z[k] * hd_10[k];

        t_23[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, ff0_2, ff1_7, gd_12, gf_14, \
                         hp0_6, hp1_6, hd_12, hd_13, hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_6 * ff0_2[k]
                  - f_7 * ff1_7[k]
                  + pa_z[k] * gf_14[k];

        t_25[k] = pb_y[k] * hd_12[k];

        t_26[k] = f_8 * gd_12[k]
                  + pb_x[k] * hd_14[k];

        t_27[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_x, pb_y, ff0_8, ff1_26, gf_23, \
                         hp0_7, hp1_7, hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * hd_14[k];

        t_29[k] = f_3 * ff0_8[k]
                  - f_4 * ff1_26[k]
                  + pa_x[k] * gf_23[k];

        t_30[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_15[k];

        t_31[k] = pb_x[k] * hd_16[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_x, pb_y, pb_z, gd_14, hp0_8, hp0_9, hp1_8, \
                         hp1_9, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_x[k] * hd_17[k];

        t_33[k] = f_0 * gd_14[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_16[k];

        t_34[k] = pb_z[k] * hd_16[k];

        t_35[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_9[k]
                  + pb_z[k] * hd_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, ff0_5, ff1_15, gf_30, hp0_10, \
                         hp1_10, hd_18, hd_19, hd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * hp0_10[k]
                  - f_2 * hp1_10[k]
                  + pb_x[k] * hd_18[k];

        t_37[k] = pb_x[k] * hd_19[k];

        t_38[k] = pb_x[k] * hd_20[k];

        t_39[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_15[k]
                  + pa_z[k] * gf_30[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pb_x, pb_y, ff0_7, ff1_20, gd_19, \
                         gf_36, hp0_11, hp1_11, hd_20, hd_21, hd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * gd_19[k]
                  + pb_y[k] * hd_20[k];

        t_41[k] = f_6 * ff0_7[k]
                  - f_7 * ff1_20[k]
                  + pa_y[k] * gf_36[k];

        t_42[k] = f_1 * hp0_11[k]
                  - f_2 * hp1_11[k]
                  + pb_x[k] * hd_21[k];

        t_43[k] = pb_x[k] * hd_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_y, pa_z, pb_x, pb_y, ff0_6, ff0_8, ff1_18, \
                         ff1_26, gd_20, gf_34, gf_38, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * hd_23[k];

        t_45[k] = f_6 * ff0_6[k]
                  - f_7 * ff1_18[k]
                  + pa_z[k] * gf_34[k];

        t_46[k] = f_8 * gd_20[k]
                  + pb_y[k] * hd_23[k];

        t_47[k] = f_3 * ff0_8[k]
                  - f_4 * ff1_26[k]
                  + pa_y[k] * gf_38[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pb_x, pb_y, hp0_12, hp0_13, hp1_12, \
                         hp1_13, hd_24, hd_25, hd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * hp0_12[k]
                  - f_2 * hp1_12[k]
                  + pb_x[k] * hd_24[k];

        t_49[k] = pb_x[k] * hd_25[k];

        t_50[k] = pb_x[k] * hd_26[k];

        t_51[k] = f_1 * hp0_13[k]
                  - f_2 * hp1_13[k]
                  + pb_y[k] * hd_25[k];

        t_52[k] = pb_y[k] * hd_26[k];
    }

#pragma omp simd aligned(t_53, pb_z, gd_23, hp0_14, hp1_14, hd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * gd_23[k]
                  + f_1 * hp0_14[k]
                  - f_2 * hp1_14[k]
                  + pb_z[k] * hd_26[k];
    }
}

auto
compute_prim_hf_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ff0, const size_t ff1,
                                     const size_t gd, const size_t gf, const size_t hp0,
                                     const size_t hp1, const size_t hd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / p;
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 1.0 / p;

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

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_9 = buffer.data(ff0 + 9);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_18 = buffer.data(ff0 + 18);
    const auto *ff0_20 = buffer.data(ff0 + 20);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_13 = buffer.data(ff1 + 13);
    const auto *ff1_18 = buffer.data(ff1 + 18);
    const auto *ff1_21 = buffer.data(ff1 + 21);
    const auto *ff1_25 = buffer.data(ff1 + 25);
    const auto *ff1_32 = buffer.data(ff1 + 32);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_41 = buffer.data(gf + 41);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, hp0_0, hp0_1, hp1_0, \
                         hp1_1, hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_4[k] = pb_y[k] * hd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, ff0_0, ff1_0, gd_6, gf_6, \
                         hp0_2, hp1_2, hd_2, hd_3, hd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_6[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pa_y[k] * gf_6[k];

        t_7[k] = pb_z[k] * hd_3[k];

        t_8[k] = f_5 * gd_6[k]
                 + pb_x[k] * hd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, ff0_9, ff1_11, gf_12, hp0_3, hp1_3, \
                         hd_4, hd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * ff0_9[k]
                 - f_7 * ff1_11[k]
                 + pa_x[k] * gf_12[k];

        t_10[k] = pb_z[k] * hd_4[k];

        t_11[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, ff0_0, ff1_0, gd_10, gf_7, \
                         hp0_4, hp1_4, hd_6, hd_7, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * ff0_0[k]
                  - f_4 * ff1_0[k]
                  + pa_z[k] * gf_7[k];

        t_13[k] = pb_y[k] * hd_6[k];

        t_14[k] = f_5 * gd_10[k]
                  + pb_x[k] * hd_8[k];

        t_15[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, ff0_6, ff0_11, ff1_7, \
                         ff1_13, gf_9, gf_20, hd_8, hd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * hd_8[k];

        t_17[k] = f_6 * ff0_11[k]
                  - f_7 * ff1_13[k]
                  + pa_x[k] * gf_20[k];

        t_18[k] = f_6 * ff0_6[k]
                  - f_7 * ff1_7[k]
                  + pa_y[k] * gf_9[k];

        t_19[k] = pb_z[k] * hd_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, ff0_15, ff1_18, gd_11, \
                         gf_22, hp0_5, hp1_5, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_8 * gd_11[k]
                  + pb_x[k] * hd_10[k];

        t_21[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_18[k]
                  + pa_x[k] * gf_22[k];

        t_22[k] = pb_z[k] * hd_10[k];

        t_23[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, ff0_7, ff1_8, gd_12, gf_15, \
                         hp0_6, hp1_6, hd_12, hd_13, hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_6 * ff0_7[k]
                  - f_7 * ff1_8[k]
                  + pa_z[k] * gf_15[k];

        t_25[k] = pb_y[k] * hd_12[k];

        t_26[k] = f_8 * gd_12[k]
                  + pb_x[k] * hd_14[k];

        t_27[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_x, pb_y, ff0_26, ff1_32, gf_24, \
                         hp0_7, hp1_7, hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * hd_14[k];

        t_29[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_32[k]
                  + pa_x[k] * gf_24[k];

        t_30[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_15[k];

        t_31[k] = pb_x[k] * hd_16[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_x, pb_y, pb_z, gd_14, hp0_8, hp0_9, hp1_8, \
                         hp1_9, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_x[k] * hd_17[k];

        t_33[k] = f_0 * gd_14[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_16[k];

        t_34[k] = pb_z[k] * hd_16[k];

        t_35[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_9[k]
                  + pb_z[k] * hd_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, ff0_15, ff1_18, gf_31, hp0_10, \
                         hp1_10, hd_18, hd_19, hd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * hp0_10[k]
                  - f_2 * hp1_10[k]
                  + pb_x[k] * hd_18[k];

        t_37[k] = pb_x[k] * hd_19[k];

        t_38[k] = pb_x[k] * hd_20[k];

        t_39[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_18[k]
                  + pa_z[k] * gf_31[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pb_x, pb_y, ff0_20, ff1_25, gd_19, \
                         gf_38, hp0_11, hp1_11, hd_20, hd_21, hd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * gd_19[k]
                  + pb_y[k] * hd_20[k];

        t_41[k] = f_6 * ff0_20[k]
                  - f_7 * ff1_25[k]
                  + pa_y[k] * gf_38[k];

        t_42[k] = f_1 * hp0_11[k]
                  - f_2 * hp1_11[k]
                  + pb_x[k] * hd_21[k];

        t_43[k] = pb_x[k] * hd_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_y, pa_z, pb_x, pb_y, ff0_18, ff0_26, \
                         ff1_21, ff1_32, gd_20, gf_36, gf_41, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * hd_23[k];

        t_45[k] = f_6 * ff0_18[k]
                  - f_7 * ff1_21[k]
                  + pa_z[k] * gf_36[k];

        t_46[k] = f_8 * gd_20[k]
                  + pb_y[k] * hd_23[k];

        t_47[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_32[k]
                  + pa_y[k] * gf_41[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pb_x, pb_y, hp0_12, hp0_13, hp1_12, \
                         hp1_13, hd_24, hd_25, hd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * hp0_12[k]
                  - f_2 * hp1_12[k]
                  + pb_x[k] * hd_24[k];

        t_49[k] = pb_x[k] * hd_25[k];

        t_50[k] = pb_x[k] * hd_26[k];

        t_51[k] = f_1 * hp0_13[k]
                  - f_2 * hp1_13[k]
                  + pb_y[k] * hd_25[k];

        t_52[k] = pb_y[k] * hd_26[k];
    }

#pragma omp simd aligned(t_53, pb_z, gd_23, hp0_14, hp1_14, hd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * gd_23[k]
                  + f_1 * hp0_14[k]
                  - f_2 * hp1_14[k]
                  + pb_z[k] * hd_26[k];
    }
}

auto
compute_prim_hf_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ff0, const size_t ff1,
                                     const size_t gd, const size_t gf, const size_t hp0,
                                     const size_t hp1, const size_t hd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / p;
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 1.0 / p;

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

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_13 = buffer.data(ff0 + 13);
    const auto *ff0_18 = buffer.data(ff0 + 18);
    const auto *ff0_21 = buffer.data(ff0 + 21);
    const auto *ff0_25 = buffer.data(ff0 + 25);
    const auto *ff0_32 = buffer.data(ff0 + 32);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_9 = buffer.data(ff1 + 9);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_18 = buffer.data(ff1 + 18);
    const auto *ff1_20 = buffer.data(ff1 + 20);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, hp0_0, hp0_1, hp1_0, \
                         hp1_1, hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_4[k] = pb_y[k] * hd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, ff0_0, ff1_0, gd_6, gf_6, \
                         hp0_2, hp1_2, hd_2, hd_3, hd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_6[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pa_y[k] * gf_6[k];

        t_7[k] = pb_z[k] * hd_3[k];

        t_8[k] = f_5 * gd_6[k]
                 + pb_x[k] * hd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, ff0_11, ff1_9, gf_11, hp0_3, hp1_3, \
                         hd_4, hd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * ff0_11[k]
                 - f_7 * ff1_9[k]
                 + pa_x[k] * gf_11[k];

        t_10[k] = pb_z[k] * hd_4[k];

        t_11[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, ff0_0, ff1_0, gd_10, gf_7, \
                         hp0_4, hp1_4, hd_6, hd_7, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * ff0_0[k]
                  - f_4 * ff1_0[k]
                  + pa_z[k] * gf_7[k];

        t_13[k] = pb_y[k] * hd_6[k];

        t_14[k] = f_5 * gd_10[k]
                  + pb_x[k] * hd_8[k];

        t_15[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, ff0_7, ff0_13, ff1_6, \
                         ff1_11, gf_8, gf_19, hd_8, hd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * hd_8[k];

        t_17[k] = f_6 * ff0_13[k]
                  - f_7 * ff1_11[k]
                  + pa_x[k] * gf_19[k];

        t_18[k] = f_6 * ff0_7[k]
                  - f_7 * ff1_6[k]
                  + pa_y[k] * gf_8[k];

        t_19[k] = pb_z[k] * hd_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, ff0_18, ff1_15, gd_11, \
                         gf_20, hp0_5, hp1_5, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_8 * gd_11[k]
                  + pb_x[k] * hd_10[k];

        t_21[k] = f_3 * ff0_18[k]
                  - f_4 * ff1_15[k]
                  + pa_x[k] * gf_20[k];

        t_22[k] = pb_z[k] * hd_10[k];

        t_23[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, ff0_8, ff1_7, gd_12, gf_14, \
                         hp0_6, hp1_6, hd_12, hd_13, hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_6 * ff0_8[k]
                  - f_7 * ff1_7[k]
                  + pa_z[k] * gf_14[k];

        t_25[k] = pb_y[k] * hd_12[k];

        t_26[k] = f_8 * gd_12[k]
                  + pb_x[k] * hd_14[k];

        t_27[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_x, pb_y, ff0_32, ff1_26, gf_21, \
                         hp0_7, hp1_7, hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * hd_14[k];

        t_29[k] = f_3 * ff0_32[k]
                  - f_4 * ff1_26[k]
                  + pa_x[k] * gf_21[k];

        t_30[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_15[k];

        t_31[k] = pb_x[k] * hd_16[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_x, pb_y, pb_z, gd_14, hp0_8, hp0_9, hp1_8, \
                         hp1_9, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_x[k] * hd_17[k];

        t_33[k] = f_0 * gd_14[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_16[k];

        t_34[k] = pb_z[k] * hd_16[k];

        t_35[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_9[k]
                  + pb_z[k] * hd_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, ff0_18, ff1_15, gf_28, hp0_10, \
                         hp1_10, hd_18, hd_19, hd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * hp0_10[k]
                  - f_2 * hp1_10[k]
                  + pb_x[k] * hd_18[k];

        t_37[k] = pb_x[k] * hd_19[k];

        t_38[k] = pb_x[k] * hd_20[k];

        t_39[k] = f_3 * ff0_18[k]
                  - f_4 * ff1_15[k]
                  + pa_z[k] * gf_28[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pb_x, pb_y, ff0_25, ff1_20, gd_19, \
                         gf_34, hp0_11, hp1_11, hd_20, hd_21, hd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * gd_19[k]
                  + pb_y[k] * hd_20[k];

        t_41[k] = f_6 * ff0_25[k]
                  - f_7 * ff1_20[k]
                  + pa_y[k] * gf_34[k];

        t_42[k] = f_1 * hp0_11[k]
                  - f_2 * hp1_11[k]
                  + pb_x[k] * hd_21[k];

        t_43[k] = pb_x[k] * hd_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_y, pa_z, pb_x, pb_y, ff0_21, ff0_32, \
                         ff1_18, ff1_26, gd_20, gf_32, gf_35, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * hd_23[k];

        t_45[k] = f_6 * ff0_21[k]
                  - f_7 * ff1_18[k]
                  + pa_z[k] * gf_32[k];

        t_46[k] = f_8 * gd_20[k]
                  + pb_y[k] * hd_23[k];

        t_47[k] = f_3 * ff0_32[k]
                  - f_4 * ff1_26[k]
                  + pa_y[k] * gf_35[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pb_x, pb_y, hp0_12, hp0_13, hp1_12, \
                         hp1_13, hd_24, hd_25, hd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * hp0_12[k]
                  - f_2 * hp1_12[k]
                  + pb_x[k] * hd_24[k];

        t_49[k] = pb_x[k] * hd_25[k];

        t_50[k] = pb_x[k] * hd_26[k];

        t_51[k] = f_1 * hp0_13[k]
                  - f_2 * hp1_13[k]
                  + pb_y[k] * hd_25[k];

        t_52[k] = pb_y[k] * hd_26[k];
    }

#pragma omp simd aligned(t_53, pb_z, gd_23, hp0_14, hp1_14, hd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * gd_23[k]
                  + f_1 * hp0_14[k]
                  - f_2 * hp1_14[k]
                  + pb_z[k] * hd_26[k];
    }
}

auto
compute_prim_hf_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ff0, const size_t ff1,
                                     const size_t gd, const size_t gf, const size_t hp0,
                                     const size_t hp1, const size_t hd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);

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

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_5 = buffer.data(ff1 + 5);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_8 = buffer.data(ff1 + 8);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
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
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);

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

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_14 = buffer.data(hp1 + 14);
    const auto *hp1_18 = buffer.data(hp1 + 18);
    const auto *hp1_19 = buffer.data(hp1 + 19);
    const auto *hp1_20 = buffer.data(hp1 + 20);
    const auto *hp1_22 = buffer.data(hp1 + 22);
    const auto *hp1_24 = buffer.data(hp1 + 24);
    const auto *hp1_27 = buffer.data(hp1 + 27);
    const auto *hp1_28 = buffer.data(hp1 + 28);
    const auto *hp1_29 = buffer.data(hp1 + 29);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, gd_0, gd_1, gd_2, hp0_0, hp0_1, \
                         hp1_0, hp1_1, hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = f_0 * gd_1[k]
                 + pb_x[k] * hd_1[k];

        t_2[k] = f_0 * gd_2[k]
                 + pb_x[k] * hd_2[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pb_x, pb_y, pb_z, gd_0, gd_4, gf_0, hp0_2, \
                         hp1_2, hd_2, hd_3, hd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_5[k] = pa_y[k] * gf_0[k];

        t_6[k] = f_3 * gd_0[k]
                 + pb_y[k] * hd_3[k];

        t_7[k] = f_4 * gd_4[k]
                 + pb_x[k] * hd_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_y, pa_z, pb_x, pb_z, gd_0, gd_1, gd_6, gf_0, \
                         gf_1, hd_5, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * gd_1[k]
                 + pa_y[k] * gf_1[k];

        t_9[k] = pa_z[k] * gf_0[k];

        t_10[k] = f_3 * gd_0[k]
                  + pb_z[k] * hd_5[k];

        t_11[k] = f_4 * gd_6[k]
                  + pb_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pa_z, pb_y, ff0_0, ff1_0, gd_2, gd_3, gf_2, \
                         gf_3, hd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * gd_2[k]
                  + pa_z[k] * gf_2[k];

        t_13[k] = f_6 * ff0_0[k]
                  - f_7 * ff1_0[k]
                  + pa_y[k] * gf_3[k];

        t_14[k] = f_8 * gd_3[k]
                  + pb_y[k] * hd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pb_x, pb_z, ff0_3, ff1_3, gd_8, gf_6, hp0_3, \
                         hp1_6, hd_8, hd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * gd_8[k]
                  + pb_x[k] * hd_8[k];

        t_16[k] = f_9 * ff0_3[k]
                  - f_10 * ff1_3[k]
                  + pa_x[k] * gf_6[k];

        t_17[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_6[k]
                  + pb_z[k] * hd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_z, pb_x, pb_z, ff0_0, ff1_0, gd_5, gd_12, gf_4, \
                         hd_10, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_6 * ff0_0[k]
                  - f_7 * ff1_0[k]
                  + pa_z[k] * gf_4[k];

        t_19[k] = f_8 * gd_5[k]
                  + pb_z[k] * hd_10[k];

        t_20[k] = f_5 * gd_12[k]
                  + pb_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_y, ff0_1, ff0_4, ff1_1, ff1_4, gf_5, \
                         gf_8, hp0_4, hp1_8, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_11[k];

        t_22[k] = f_9 * ff0_4[k]
                  - f_10 * ff1_4[k]
                  + pa_x[k] * gf_8[k];

        t_23[k] = f_9 * ff0_1[k]
                  - f_10 * ff1_1[k]
                  + pa_y[k] * gf_5[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pb_x, pb_y, ff0_5, ff1_5, gd_7, gd_14, gf_9, \
                         hd_13, hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * gd_7[k]
                  + pb_y[k] * hd_13[k];

        t_25[k] = f_8 * gd_14[k]
                  + pb_x[k] * hd_14[k];

        t_26[k] = f_6 * ff0_5[k]
                  - f_7 * ff1_5[k]
                  + pa_x[k] * gf_9[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_y, pa_z, pb_z, ff0_2, ff1_2, gd_10, gf_7, \
                         hp0_5, hp1_11, hd_15, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_11[k]
                  + pb_z[k] * hd_15[k];

        t_28[k] = pa_y[k] * gf_7[k];

        t_29[k] = f_9 * ff0_2[k]
                  - f_10 * ff1_2[k]
                  + pa_z[k] * gf_7[k];

        t_30[k] = f_5 * gd_10[k]
                  + pb_z[k] * hd_17[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_x, pb_x, pb_y, ff0_8, ff1_8, gd_16, gf_10, \
                         hp0_6, hp1_14, hd_18, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_8 * gd_16[k]
                  + pb_x[k] * hd_19[k];

        t_32[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_14[k]
                  + pb_y[k] * hd_18[k];

        t_33[k] = f_6 * ff0_8[k]
                  - f_7 * ff1_8[k]
                  + pa_x[k] * gf_10[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_x, pb_x, pb_y, gd_13, gd_17, gd_18, \
                         gf_11, gf_12, gf_15, hd_20, hd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_5 * gd_17[k]
                  + pa_x[k] * gf_11[k];

        t_35[k] = f_4 * gd_13[k]
                  + pb_y[k] * hd_20[k];

        t_36[k] = f_3 * gd_18[k]
                  + pb_x[k] * hd_21[k];

        t_37[k] = pa_x[k] * gf_12[k];

        t_38[k] = pa_x[k] * gf_15[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_x, pb_x, pb_z, gd_15, gd_30, gd_32, \
                         gf_16, gf_18, gf_20, hd_24, hd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pa_x[k] * gf_16[k];

        t_40[k] = f_5 * gd_30[k]
                  + pa_x[k] * gf_18[k];

        t_41[k] = f_4 * gd_15[k]
                  + pb_z[k] * hd_24[k];

        t_42[k] = f_3 * gd_32[k]
                  + pb_x[k] * hd_25[k];

        t_43[k] = pa_x[k] * gf_20[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_y, gd_17, gd_18, gd_19, hp0_7, \
                         hp0_8, hp1_18, hp1_19, hd_26, hd_27, hd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_18[k]
                  + pb_x[k] * hd_26[k];

        t_45[k] = f_0 * gd_17[k]
                  + pb_y[k] * hd_26[k];

        t_46[k] = f_0 * gd_18[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_19[k]
                  + pb_y[k] * hd_27[k];

        t_47[k] = f_0 * gd_19[k]
                  + pb_y[k] * hd_28[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_y, pb_z, gd_18, gd_22, gf_12, hp0_9, \
                         hp1_20, hd_28, hd_29, hd_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_20[k]
                  + pb_z[k] * hd_28[k];

        t_49[k] = pa_z[k] * gf_12[k];

        t_50[k] = f_3 * gd_18[k]
                  + pb_z[k] * hd_29[k];

        t_51[k] = f_4 * gd_22[k]
                  + pb_y[k] * hd_31[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_z, pb_x, ff0_5, ff1_5, gd_19, gf_13, gf_14, \
                         hp0_10, hp1_22, hd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_5 * gd_19[k]
                  + pa_z[k] * gf_13[k];

        t_53[k] = f_1 * hp0_10[k]
                  - f_2 * hp1_22[k]
                  + pb_x[k] * hd_32[k];

        t_54[k] = f_6 * ff0_5[k]
                  - f_7 * ff1_5[k]
                  + pa_z[k] * gf_14[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_y, pb_y, pb_z, ff0_7, ff1_7, gd_20, gd_26, \
                         gf_16, hd_33, hd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_8 * gd_20[k]
                  + pb_z[k] * hd_33[k];

        t_56[k] = f_5 * gd_26[k]
                  + pb_y[k] * hd_35[k];

        t_57[k] = f_9 * ff0_7[k]
                  - f_10 * ff1_7[k]
                  + pa_y[k] * gf_16[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_z, pb_x, pb_z, ff0_6, ff1_6, gd_24, gf_15, \
                         hp0_11, hp1_24, hd_36, hd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * hp0_11[k]
                  - f_2 * hp1_24[k]
                  + pb_x[k] * hd_36[k];

        t_59[k] = f_9 * ff0_6[k]
                  - f_10 * ff1_6[k]
                  + pa_z[k] * gf_15[k];

        t_60[k] = f_5 * gd_24[k]
                  + pb_z[k] * hd_37[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_y, pb_y, pb_z, ff0_8, ff1_8, gd_27, gd_29, \
                         gd_31, gf_17, gf_19, hd_39, hd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_8 * gd_29[k]
                  + pb_y[k] * hd_39[k];

        t_62[k] = f_6 * ff0_8[k]
                  - f_7 * ff1_8[k]
                  + pa_y[k] * gf_17[k];

        t_63[k] = f_5 * gd_31[k]
                  + pa_y[k] * gf_19[k];

        t_64[k] = f_4 * gd_27[k]
                  + pb_z[k] * hd_40[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_y, pb_x, pb_y, pb_z, gd_30, gd_32, gf_20, \
                         hp0_12, hp1_27, hd_42, hd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_3 * gd_32[k]
                  + pb_y[k] * hd_42[k];

        t_66[k] = pa_y[k] * gf_20[k];

        t_67[k] = f_1 * hp0_12[k]
                  - f_2 * hp1_27[k]
                  + pb_x[k] * hd_43[k];

        t_68[k] = f_0 * gd_30[k]
                  + pb_z[k] * hd_43[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_y, pb_z, gd_31, gd_32, hp0_13, hp0_14, hp1_28, \
                         hp1_29, hd_44, hd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * hp0_13[k]
                  - f_2 * hp1_28[k]
                  + pb_y[k] * hd_44[k];

        t_70[k] = f_0 * gd_31[k]
                  + pb_z[k] * hd_44[k];

        t_71[k] = f_0 * gd_32[k]
                  + f_1 * hp0_14[k]
                  - f_2 * hp1_29[k]
                  + pb_z[k] * hd_45[k];
    }
}

auto
compute_prim_hf_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ff0, const size_t ff1,
                                     const size_t gd, const size_t gf, const size_t hp0,
                                     const size_t hp1, const size_t hd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 1.0 / p;
    const auto f_10 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_10 = buffer.data(ff1 + 10);
    const auto *ff1_11 = buffer.data(ff1 + 11);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_4 = buffer.data(gd + 4);
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
    const auto *gd_21 = buffer.data(gd + 21);
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
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_30 = buffer.data(gf + 30);
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
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_54 = buffer.data(gf + 54);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, hp0_0, hp0_1, hp1_0, \
                         hp1_1, hd_0, hd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_4[k] = pb_z[k] * hd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, gd_1, gd_2, gf_0, gf_3, \
                         hp0_2, hp1_2, hd_2, hd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_y[k] * hd_2[k];

        t_6[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_7[k] = pa_y[k] * gf_0[k];

        t_8[k] = f_3 * gd_1[k]
                 + pa_y[k] * gf_3[k];

        t_9[k] = f_4 * gd_2[k]
                 + pb_y[k] * hd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pa_z, pb_z, gd_0, gd_1, gf_0, \
                         gf_3, gf_5, hd_6, hd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * gf_5[k];

        t_11[k] = pa_z[k] * gf_0[k];

        t_12[k] = f_4 * gd_0[k]
                  + pb_z[k] * hd_6[k];

        t_13[k] = pa_z[k] * gf_3[k];

        t_14[k] = f_4 * gd_1[k]
                  + pb_z[k] * hd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_y, pa_z, pb_y, pb_z, ff0_0, ff1_0, gd_2, \
                         gf_5, gf_6, hd_8, hd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_y[k] * hd_8[k];

        t_16[k] = f_3 * gd_2[k]
                  + pa_z[k] * gf_5[k];

        t_17[k] = f_5 * ff0_0[k]
                  - f_6 * ff1_0[k]
                  + pa_y[k] * gf_6[k];

        t_18[k] = pb_z[k] * hd_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, pb_y, pb_z, ff0_4, ff1_4, gd_5, \
                         gd_10, gf_13, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * gd_10[k]
                  + pb_x[k] * hd_10[k];

        t_20[k] = f_7 * ff0_4[k]
                  - f_8 * ff1_4[k]
                  + pa_x[k] * gf_13[k];

        t_21[k] = pb_z[k] * hd_10[k];

        t_22[k] = f_9 * gd_5[k]
                  + pb_y[k] * hd_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pa_z, pb_z, gd_4, gf_7, gf_9, hp0_3, \
                         hp1_3, hd_11, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_11[k];

        t_24[k] = pa_y[k] * gf_9[k];

        t_25[k] = pa_z[k] * gf_7[k];

        t_26[k] = f_4 * gd_4[k]
                  + pb_z[k] * hd_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, ff0_0, ff1_0, gd_8, gf_8, \
                         gf_10, hd_13, hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_4 * gd_8[k]
                  + pb_y[k] * hd_13[k];

        t_28[k] = pa_y[k] * gf_10[k];

        t_29[k] = f_5 * ff0_0[k]
                  - f_6 * ff1_0[k]
                  + pa_z[k] * gf_8[k];

        t_30[k] = pb_y[k] * hd_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pb_x, pb_y, pb_z, gd_6, gd_7, gd_16, \
                         hp0_4, hp1_4, hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_9 * gd_6[k]
                  + pb_z[k] * hd_14[k];

        t_32[k] = f_3 * gd_16[k]
                  + pb_x[k] * hd_16[k];

        t_33[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_15[k];

        t_34[k] = f_9 * gd_7[k]
                  + pb_z[k] * hd_15[k];

        t_35[k] = pb_y[k] * hd_16[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_x, pa_y, pb_z, ff0_1, ff0_6, ff1_1, ff1_6, \
                         gf_11, gf_19, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_7 * ff0_6[k]
                  - f_8 * ff1_6[k]
                  + pa_x[k] * gf_19[k];

        t_37[k] = f_7 * ff0_1[k]
                  - f_8 * ff1_1[k]
                  + pa_y[k] * gf_11[k];

        t_38[k] = pb_z[k] * hd_17[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pb_x, pb_y, pb_z, ff0_7, ff1_7, gd_11, \
                         gd_18, gf_22, hd_18, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * gd_18[k]
                  + pb_x[k] * hd_18[k];

        t_40[k] = f_5 * ff0_7[k]
                  - f_6 * ff1_7[k]
                  + pa_x[k] * gf_22[k];

        t_41[k] = pb_z[k] * hd_18[k];

        t_42[k] = f_3 * gd_11[k]
                  + pb_y[k] * hd_19[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pa_z, pb_z, gd_9, gd_10, gf_11, gf_13, \
                         hp0_5, hp1_5, hd_19, hd_20, hd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_19[k];

        t_44[k] = pa_z[k] * gf_11[k];

        t_45[k] = f_4 * gd_9[k]
                  + pb_z[k] * hd_20[k];

        t_46[k] = pa_z[k] * gf_13[k];

        t_47[k] = f_4 * gd_10[k]
                  + pb_z[k] * hd_21[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_y, pa_z, pb_y, gd_11, gd_13, gd_15, \
                         gf_14, gf_15, gf_16, gf_18, hd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_9 * gd_13[k]
                  + pb_y[k] * hd_22[k];

        t_49[k] = f_3 * gd_11[k]
                  + pa_z[k] * gf_14[k];

        t_50[k] = pa_y[k] * gf_15[k];

        t_51[k] = pa_y[k] * gf_16[k];

        t_52[k] = f_3 * gd_15[k]
                  + pa_y[k] * gf_18[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_y, pa_z, pb_y, pb_z, ff0_2, ff1_2, gd_12, \
                         gd_16, gf_15, gf_19, hd_23, hd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_9 * gd_12[k]
                  + pb_z[k] * hd_23[k];

        t_54[k] = f_4 * gd_16[k]
                  + pb_y[k] * hd_24[k];

        t_55[k] = pa_y[k] * gf_19[k];

        t_56[k] = f_7 * ff0_2[k]
                  - f_8 * ff1_2[k]
                  + pa_z[k] * gf_15[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pb_x, pb_y, pb_z, gd_14, gd_15, gd_21, \
                         hp0_6, hp1_6, hd_25, hd_26, hd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_y[k] * hd_25[k];

        t_58[k] = f_3 * gd_14[k]
                  + pb_z[k] * hd_25[k];

        t_59[k] = f_9 * gd_21[k]
                  + pb_x[k] * hd_27[k];

        t_60[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_26[k];

        t_61[k] = f_3 * gd_15[k]
                  + pb_z[k] * hd_26[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pb_x, pb_y, ff0_11, ff1_11, gd_22, \
                         gd_23, gf_26, gf_27, hd_27, hd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_y[k] * hd_27[k];

        t_63[k] = f_5 * ff0_11[k]
                  - f_6 * ff1_11[k]
                  + pa_x[k] * gf_26[k];

        t_64[k] = f_3 * gd_22[k]
                  + pa_x[k] * gf_27[k];

        t_65[k] = f_4 * gd_23[k]
                  + pb_x[k] * hd_29[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, pa_x, pa_z, pb_z, gd_17, gf_20, \
                         gf_30, gf_32, gf_33, gf_35, hd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_x[k] * gf_30[k];

        t_67[k] = pa_x[k] * gf_32[k];

        t_68[k] = pa_x[k] * gf_33[k];

        t_69[k] = pa_z[k] * gf_20[k];

        t_70[k] = f_4 * gd_17[k]
                  + pb_z[k] * hd_30[k];

        t_71[k] = pa_x[k] * gf_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, pa_x, pb_z, gd_19, gd_28, gf_36, \
                         gf_37, gf_38, gf_39, gf_40, hd_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_x[k] * gf_36[k];

        t_73[k] = pa_x[k] * gf_37[k];

        t_74[k] = f_3 * gd_28[k]
                  + pa_x[k] * gf_38[k];

        t_75[k] = f_9 * gd_19[k]
                  + pb_z[k] * hd_31[k];

        t_76[k] = pa_x[k] * gf_39[k];

        t_77[k] = pa_x[k] * gf_40[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, t_84, pa_x, pa_y, gf_23, gf_24, \
                         gf_41, gf_42, gf_43, gf_44, gf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_x[k] * gf_41[k];

        t_79[k] = pa_x[k] * gf_42[k];

        t_80[k] = pa_y[k] * gf_23[k];

        t_81[k] = pa_y[k] * gf_24[k];

        t_82[k] = pa_x[k] * gf_43[k];

        t_83[k] = pa_x[k] * gf_44[k];

        t_84[k] = pa_x[k] * gf_45[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pa_x, pb_x, pb_z, gd_20, gd_33, gd_35, \
                         gf_47, gf_51, gf_52, hd_32, hd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_3 * gd_33[k]
                  + pa_x[k] * gf_47[k];

        t_86[k] = f_10 * gd_20[k]
                  + pb_z[k] * hd_32[k];

        t_87[k] = f_4 * gd_35[k]
                  + pb_x[k] * hd_33[k];

        t_88[k] = pa_x[k] * gf_51[k];

        t_89[k] = pa_x[k] * gf_52[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pa_x, pb_x, pb_z, gf_54, hp0_7, hp1_7, \
                         hd_34, hd_35, hd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pa_x[k] * gf_54[k];

        t_91[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_34[k];

        t_92[k] = pb_z[k] * hd_34[k];

        t_93[k] = pb_x[k] * hd_35[k];

        t_94[k] = pb_x[k] * hd_36[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pb_y, pb_z, gd_23, gd_24, hp0_8, hp0_9, \
                         hp1_8, hp1_9, hd_35, hd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_0 * gd_23[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_35[k];

        t_96[k] = pb_z[k] * hd_35[k];

        t_97[k] = f_0 * gd_24[k]
                  + pb_y[k] * hd_36[k];

        t_98[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_9[k]
                  + pb_z[k] * hd_36[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pa_z, pb_x, pb_z, gd_22, gd_23, \
                         gf_27, gf_30, hd_37, hd_38, hd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * gf_27[k];

        t_100[k] = f_4 * gd_22[k]
                   + pb_z[k] * hd_37[k];

        t_101[k] = pb_x[k] * hd_39[k];

        t_102[k] = pa_z[k] * gf_30[k];

        t_103[k] = f_4 * gd_23[k]
                   + pb_z[k] * hd_38[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_z, pb_x, pb_y, pb_z, gd_24, gd_25, \
                         gd_27, gf_33, hp0_10, hp1_10, hd_39, hd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_10 * gd_27[k]
                   + pb_y[k] * hd_39[k];

        t_105[k] = f_3 * gd_24[k]
                   + pa_z[k] * gf_33[k];

        t_106[k] = f_1 * hp0_10[k]
                   - f_2 * hp1_10[k]
                   + pb_x[k] * hd_40[k];

        t_107[k] = f_9 * gd_25[k]
                   + pb_z[k] * hd_40[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_z, pb_x, pb_y, pb_z, ff0_7, \
                         ff1_7, gd_26, gd_30, gf_34, hd_41, hd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pb_x[k] * hd_41[k];

        t_109[k] = pb_x[k] * hd_42[k];

        t_110[k] = f_5 * ff0_7[k]
                   - f_6 * ff1_7[k]
                   + pa_z[k] * gf_34[k];

        t_111[k] = f_9 * gd_26[k]
                   + pb_z[k] * hd_41[k];

        t_112[k] = f_3 * gd_30[k]
                   + pb_y[k] * hd_42[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_y, pb_x, pb_z, ff0_10, ff1_10, gd_28, \
                         gf_42, hp0_11, hp1_11, hd_43, hd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_7 * ff0_10[k]
                   - f_8 * ff1_10[k]
                   + pa_y[k] * gf_42[k];

        t_114[k] = f_1 * hp0_11[k]
                   - f_2 * hp1_11[k]
                   + pb_x[k] * hd_43[k];

        t_115[k] = f_3 * gd_28[k]
                   + pb_z[k] * hd_43[k];

        t_116[k] = pb_x[k] * hd_44[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_z, pb_x, pb_y, pb_z, ff0_8, ff1_8, \
                         gd_29, gd_32, gf_39, hd_44, hd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pb_x[k] * hd_45[k];

        t_118[k] = f_7 * ff0_8[k]
                   - f_8 * ff1_8[k]
                   + pa_z[k] * gf_39[k];

        t_119[k] = f_3 * gd_29[k]
                   + pb_z[k] * hd_44[k];

        t_120[k] = f_9 * gd_32[k]
                   + pb_y[k] * hd_45[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pa_y, pb_x, ff0_11, ff1_11, gd_34, \
                         gf_46, gf_47, gf_49, gf_51, hd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_5 * ff0_11[k]
                   - f_6 * ff1_11[k]
                   + pa_y[k] * gf_46[k];

        t_122[k] = pa_y[k] * gf_47[k];

        t_123[k] = pa_y[k] * gf_49[k];

        t_124[k] = pb_x[k] * hd_46[k];

        t_125[k] = f_3 * gd_34[k]
                   + pa_y[k] * gf_51[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pa_y, pb_x, pb_y, pb_z, gd_31, gd_35, \
                         gf_54, hp0_12, hp1_12, hd_46, hd_47, hd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_10 * gd_31[k]
                   + pb_z[k] * hd_46[k];

        t_127[k] = f_4 * gd_35[k]
                   + pb_y[k] * hd_47[k];

        t_128[k] = pa_y[k] * gf_54[k];

        t_129[k] = f_1 * hp0_12[k]
                   - f_2 * hp1_12[k]
                   + pb_x[k] * hd_48[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, t_135, pb_x, pb_y, pb_z, gd_33, \
                         gd_34, hp0_13, hp1_13, hd_48, hd_49, hd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pb_y[k] * hd_48[k];

        t_131[k] = f_0 * gd_33[k]
                   + pb_z[k] * hd_48[k];

        t_132[k] = pb_x[k] * hd_49[k];

        t_133[k] = pb_x[k] * hd_50[k];

        t_134[k] = f_1 * hp0_13[k]
                   - f_2 * hp1_13[k]
                   + pb_y[k] * hd_49[k];

        t_135[k] = f_0 * gd_34[k]
                   + pb_z[k] * hd_49[k];
    }

#pragma omp simd aligned(t_136, t_137, pb_y, pb_z, gd_35, hp0_14, hp1_14, \
                         hd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pb_y[k] * hd_50[k];

        t_137[k] = f_0 * gd_35[k]
                   + f_1 * hp0_14[k]
                   - f_2 * hp1_14[k]
                   + pb_z[k] * hd_50[k];
    }
}

auto
compute_prim_hf_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t ff0, const size_t ff1,
                                      const size_t gd, const size_t gf, const size_t hp0,
                                      const size_t hp1, const size_t hd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);
    const auto f_9 = 1.0 / p;
    const auto f_10 = 2.0 / p;

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

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_3 = buffer.data(ff1 + 3);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_10 = buffer.data(ff1 + 10);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_14 = buffer.data(ff1 + 14);
    const auto *ff1_17 = buffer.data(ff1 + 17);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_32 = buffer.data(gf + 32);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, hp0_0, hp0_1, hp1_0, \
                         hp1_1, hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_4[k] = pb_y[k] * hd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, gd_0, gd_1, gf_0, gf_3, \
                         hp0_2, hp1_2, hd_2, hd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_6[k] = pa_y[k] * gf_0[k];

        t_7[k] = f_3 * gd_1[k]
                 + pa_y[k] * gf_3[k];

        t_8[k] = pa_z[k] * gf_0[k];

        t_9[k] = f_4 * gd_0[k]
                 + pb_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_x, pb_z, ff0_0, ff1_0, gd_2, \
                         gd_8, gf_5, gf_6, hd_7, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * gd_2[k]
                  + pa_z[k] * gf_5[k];

        t_11[k] = f_5 * ff0_0[k]
                  - f_6 * ff1_0[k]
                  + pa_y[k] * gf_6[k];

        t_12[k] = pb_z[k] * hd_7[k];

        t_13[k] = f_3 * gd_8[k]
                  + pb_x[k] * hd_8[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pb_z, ff0_4, ff1_6, gf_10, hp0_3, hp1_3, \
                         hd_8, hd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_7 * ff0_4[k]
                  - f_8 * ff1_6[k]
                  + pa_x[k] * gf_10[k];

        t_15[k] = pb_z[k] * hd_8[k];

        t_16[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_9[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_z, pb_x, pb_y, pb_z, ff0_0, ff1_0, gd_5, \
                         gd_12, gf_7, hd_10, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * ff0_0[k]
                  - f_6 * ff1_0[k]
                  + pa_z[k] * gf_7[k];

        t_18[k] = pb_y[k] * hd_10[k];

        t_19[k] = f_9 * gd_5[k]
                  + pb_z[k] * hd_10[k];

        t_20[k] = f_3 * gd_12[k]
                  + pb_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pb_y, ff0_6, ff1_8, gf_13, hp0_4, hp1_4, \
                         hd_11, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_11[k];

        t_22[k] = pb_y[k] * hd_12[k];

        t_23[k] = f_7 * ff0_6[k]
                  - f_8 * ff1_8[k]
                  + pa_x[k] * gf_13[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_y, pb_x, pb_z, ff0_1, ff1_3, gd_14, gf_8, hd_13, \
                         hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_7 * ff0_1[k]
                  - f_8 * ff1_3[k]
                  + pa_y[k] * gf_8[k];

        t_25[k] = pb_z[k] * hd_13[k];

        t_26[k] = f_9 * gd_14[k]
                  + pb_x[k] * hd_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pa_y, pb_z, ff0_7, ff1_10, gf_11, \
                         gf_14, hp0_5, hp1_5, hd_14, hd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * ff0_7[k]
                  - f_6 * ff1_10[k]
                  + pa_x[k] * gf_14[k];

        t_28[k] = pb_z[k] * hd_14[k];

        t_29[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_15[k];

        t_30[k] = pa_y[k] * gf_11[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_z, pb_x, pb_y, pb_z, ff0_2, ff1_4, gd_10, \
                         gd_16, gf_11, hd_17, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_7 * ff0_2[k]
                  - f_8 * ff1_4[k]
                  + pa_z[k] * gf_11[k];

        t_32[k] = pb_y[k] * hd_17[k];

        t_33[k] = f_3 * gd_10[k]
                  + pb_z[k] * hd_17[k];

        t_34[k] = f_9 * gd_16[k]
                  + pb_x[k] * hd_19[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_y, ff0_11, ff1_17, gd_17, gf_15, \
                         gf_16, hp0_6, hp1_6, hd_18, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_18[k];

        t_36[k] = pb_y[k] * hd_19[k];

        t_37[k] = f_5 * ff0_11[k]
                  - f_6 * ff1_17[k]
                  + pa_x[k] * gf_15[k];

        t_38[k] = f_3 * gd_17[k]
                  + pa_x[k] * gf_16[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, t_44, pa_x, pb_z, gd_15, gd_27, gf_19, \
                         gf_23, gf_25, gf_27, gf_32, hd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pa_x[k] * gf_19[k];

        t_40[k] = pa_x[k] * gf_23[k];

        t_41[k] = pa_x[k] * gf_25[k];

        t_42[k] = f_3 * gd_27[k]
                  + pa_x[k] * gf_27[k];

        t_43[k] = f_10 * gd_15[k]
                  + pb_z[k] * hd_24[k];

        t_44[k] = pa_x[k] * gf_32[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pb_x, pb_y, pb_z, gd_18, hp0_7, hp0_8, \
                         hp1_7, hp1_8, hd_26, hd_27, hd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_26[k];

        t_46[k] = pb_x[k] * hd_27[k];

        t_47[k] = pb_x[k] * hd_28[k];

        t_48[k] = f_0 * gd_18[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_27[k];

        t_49[k] = pb_z[k] * hd_27[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_z, pb_y, pb_z, gd_18, gd_19, gf_19, hp0_9, \
                         hp1_9, hd_28, hd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_0 * gd_19[k]
                  + pb_y[k] * hd_28[k];

        t_51[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_9[k]
                  + pb_z[k] * hd_28[k];

        t_52[k] = pa_z[k] * gf_19[k];

        t_53[k] = f_4 * gd_18[k]
                  + pb_z[k] * hd_29[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_z, pb_x, pb_y, gd_19, gd_21, gf_21, \
                         hp0_10, hp1_10, hd_30, hd_31, hd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_10 * gd_21[k]
                  + pb_y[k] * hd_30[k];

        t_55[k] = f_3 * gd_19[k]
                  + pa_z[k] * gf_21[k];

        t_56[k] = f_1 * hp0_10[k]
                  - f_2 * hp1_10[k]
                  + pb_x[k] * hd_31[k];

        t_57[k] = pb_x[k] * hd_32[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pb_x, pb_y, pb_z, ff0_7, ff1_10, gd_20, \
                         gd_24, gf_22, hd_32, hd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_x[k] * hd_33[k];

        t_59[k] = f_5 * ff0_7[k]
                  - f_6 * ff1_10[k]
                  + pa_z[k] * gf_22[k];

        t_60[k] = f_9 * gd_20[k]
                  + pb_z[k] * hd_32[k];

        t_61[k] = f_3 * gd_24[k]
                  + pb_y[k] * hd_33[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_y, pb_x, ff0_10, ff1_14, gf_25, hp0_11, \
                         hp1_11, hd_34, hd_35, hd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_7 * ff0_10[k]
                  - f_8 * ff1_14[k]
                  + pa_y[k] * gf_25[k];

        t_63[k] = f_1 * hp0_11[k]
                  - f_2 * hp1_11[k]
                  + pb_x[k] * hd_34[k];

        t_64[k] = pb_x[k] * hd_35[k];

        t_65[k] = pb_x[k] * hd_36[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_z, pb_y, pb_z, ff0_8, ff1_12, gd_23, gd_26, \
                         gf_23, hd_35, hd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_7 * ff0_8[k]
                  - f_8 * ff1_12[k]
                  + pa_z[k] * gf_23[k];

        t_67[k] = f_3 * gd_23[k]
                  + pb_z[k] * hd_35[k];

        t_68[k] = f_9 * gd_26[k]
                  + pb_y[k] * hd_36[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_y, pb_y, pb_z, ff0_11, ff1_17, gd_25, \
                         gd_28, gd_29, gf_26, gf_30, hd_37, hd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_5 * ff0_11[k]
                  - f_6 * ff1_17[k]
                  + pa_y[k] * gf_26[k];

        t_70[k] = f_3 * gd_28[k]
                  + pa_y[k] * gf_30[k];

        t_71[k] = f_10 * gd_25[k]
                  + pb_z[k] * hd_37[k];

        t_72[k] = f_4 * gd_29[k]
                  + pb_y[k] * hd_38[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_y, pb_x, pb_z, gd_27, gf_32, hp0_12, \
                         hp1_12, hd_39, hd_40, hd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pa_y[k] * gf_32[k];

        t_74[k] = f_1 * hp0_12[k]
                  - f_2 * hp1_12[k]
                  + pb_x[k] * hd_39[k];

        t_75[k] = f_0 * gd_27[k]
                  + pb_z[k] * hd_39[k];

        t_76[k] = pb_x[k] * hd_40[k];

        t_77[k] = pb_x[k] * hd_41[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_y, pb_z, gd_28, gd_29, hp0_13, hp0_14, \
                         hp1_13, hp1_14, hd_40, hd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_1 * hp0_13[k]
                  - f_2 * hp1_13[k]
                  + pb_y[k] * hd_40[k];

        t_79[k] = f_0 * gd_28[k]
                  + pb_z[k] * hd_40[k];

        t_80[k] = pb_y[k] * hd_41[k];

        t_81[k] = f_0 * gd_29[k]
                  + f_1 * hp0_14[k]
                  - f_2 * hp1_14[k]
                  + pb_z[k] * hd_41[k];
    }
}

auto
compute_prim_hf_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t ff0, const size_t ff1,
                                      const size_t gd, const size_t gf, const size_t hp0,
                                      const size_t hp1, const size_t hd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / p;
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 1.0 / p;

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

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_3 = buffer.data(ff0 + 3);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_5 = buffer.data(ff0 + 5);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_10 = buffer.data(ff1 + 10);
    const auto *ff1_11 = buffer.data(ff1 + 11);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_62 = buffer.data(gf + 62);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, hp0_0, hp0_1, hp1_0, \
                         hp1_1, hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_4[k] = pb_y[k] * hd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, ff0_0, ff1_0, gf_0, gf_7, \
                         hp0_2, hp1_2, hd_2, hd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_6[k] = pa_y[k] * gf_0[k];

        t_7[k] = pa_z[k] * gf_0[k];

        t_8[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pa_y[k] * gf_7[k];

        t_9[k] = pb_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, ff0_3, ff1_4, gd_6, gf_17, \
                         hp0_3, hp1_3, hd_6, hd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gd_6[k]
                  + pb_x[k] * hd_6[k];

        t_11[k] = f_6 * ff0_3[k]
                  - f_7 * ff1_4[k]
                  + pa_x[k] * gf_17[k];

        t_12[k] = pb_z[k] * hd_6[k];

        t_13[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, ff0_0, ff1_0, gd_10, gf_10, \
                         hp0_4, hp1_4, hd_8, hd_9, hd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * ff0_0[k]
                  - f_4 * ff1_0[k]
                  + pa_z[k] * gf_10[k];

        t_15[k] = pb_y[k] * hd_8[k];

        t_16[k] = f_5 * gd_10[k]
                  + pb_x[k] * hd_10[k];

        t_17[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, ff0_1, ff0_4, ff1_1, \
                         ff1_6, gf_14, gf_27, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * hd_10[k];

        t_19[k] = f_6 * ff0_4[k]
                  - f_7 * ff1_6[k]
                  + pa_x[k] * gf_27[k];

        t_20[k] = f_6 * ff0_1[k]
                  - f_7 * ff1_1[k]
                  + pa_y[k] * gf_14[k];

        t_21[k] = pb_z[k] * hd_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, ff0_5, ff1_7, gd_11, gf_30, \
                         hp0_5, hp1_5, hd_12, hd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_8 * gd_11[k]
                  + pb_x[k] * hd_12[k];

        t_23[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_7[k]
                  + pa_x[k] * gf_30[k];

        t_24[k] = pb_z[k] * hd_12[k];

        t_25[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, ff0_2, ff1_2, gd_12, gf_22, \
                         hp0_6, hp1_6, hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_6 * ff0_2[k]
                  - f_7 * ff1_2[k]
                  + pa_z[k] * gf_22[k];

        t_27[k] = pb_y[k] * hd_14[k];

        t_28[k] = f_8 * gd_12[k]
                  + pb_x[k] * hd_16[k];

        t_29[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_y, ff0_8, ff1_11, gf_34, gf_39, \
                         gf_62, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * hd_16[k];

        t_31[k] = f_3 * ff0_8[k]
                  - f_4 * ff1_11[k]
                  + pa_x[k] * gf_34[k];

        t_32[k] = pa_x[k] * gf_39[k];

        t_33[k] = pa_x[k] * gf_62[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, gd_14, hp0_7, hp0_8, \
                         hp1_7, hp1_8, hd_19, hd_20, hd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_19[k];

        t_35[k] = pb_x[k] * hd_20[k];

        t_36[k] = pb_x[k] * hd_21[k];

        t_37[k] = f_0 * gd_14[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_20[k];

        t_38[k] = pb_z[k] * hd_20[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pb_x, pb_z, gf_39, hp0_9, hp0_10, \
                         hp1_9, hp1_10, hd_21, hd_23, hd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_9[k]
                  + pb_z[k] * hd_21[k];

        t_40[k] = pa_z[k] * gf_39[k];

        t_41[k] = f_1 * hp0_10[k]
                  - f_2 * hp1_10[k]
                  + pb_x[k] * hd_23[k];

        t_42[k] = pb_x[k] * hd_24[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_y, pa_z, pb_x, pb_y, ff0_5, ff0_7, ff1_7, \
                         ff1_10, gd_19, gf_44, gf_51, hd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_x[k] * hd_25[k];

        t_44[k] = f_3 * ff0_5[k]
                  - f_4 * ff1_7[k]
                  + pa_z[k] * gf_44[k];

        t_45[k] = f_5 * gd_19[k]
                  + pb_y[k] * hd_25[k];

        t_46[k] = f_6 * ff0_7[k]
                  - f_7 * ff1_10[k]
                  + pa_y[k] * gf_51[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_z, pb_x, ff0_6, ff1_8, gf_49, hp0_11, \
                         hp1_11, hd_26, hd_27, hd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_1 * hp0_11[k]
                  - f_2 * hp1_11[k]
                  + pb_x[k] * hd_26[k];

        t_48[k] = pb_x[k] * hd_27[k];

        t_49[k] = pb_x[k] * hd_28[k];

        t_50[k] = f_6 * ff0_6[k]
                  - f_7 * ff1_8[k]
                  + pa_z[k] * gf_49[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pb_x, pb_y, ff0_8, ff1_11, gd_20, \
                         gf_55, gf_62, hp0_12, hp1_12, hd_28, hd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_8 * gd_20[k]
                  + pb_y[k] * hd_28[k];

        t_52[k] = f_3 * ff0_8[k]
                  - f_4 * ff1_11[k]
                  + pa_y[k] * gf_55[k];

        t_53[k] = pa_y[k] * gf_62[k];

        t_54[k] = f_1 * hp0_12[k]
                  - f_2 * hp1_12[k]
                  + pb_x[k] * hd_30[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pb_x, pb_y, pb_z, gd_23, hp0_13, \
                         hp0_14, hp1_13, hp1_14, hd_31, hd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pb_x[k] * hd_31[k];

        t_56[k] = pb_x[k] * hd_32[k];

        t_57[k] = f_1 * hp0_13[k]
                  - f_2 * hp1_13[k]
                  + pb_y[k] * hd_31[k];

        t_58[k] = pb_y[k] * hd_32[k];

        t_59[k] = f_0 * gd_23[k]
                  + f_1 * hp0_14[k]
                  - f_2 * hp1_14[k]
                  + pb_z[k] * hd_32[k];
    }
}

auto
compute_prim_hf_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t ff0, const size_t ff1,
                                      const size_t gd, const size_t gf, const size_t hp0,
                                      const size_t hp1, const size_t hd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
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

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_19 = buffer.data(ff1 + 19);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_26 = buffer.data(ff1 + 26);
    const auto *ff1_32 = buffer.data(ff1 + 32);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_53 = buffer.data(gf + 53);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, hp0_0, hp0_1, hp1_0, \
                         hp1_1, hd_0, hd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_4[k] = pb_z[k] * hd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, gd_1, gf_0, gf_3, gf_6, \
                         hp0_2, hp1_2, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_y[k] * hd_2[k];

        t_6[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_7[k] = pa_y[k] * gf_0[k];

        t_8[k] = f_3 * gd_1[k]
                 + pa_y[k] * gf_3[k];

        t_9[k] = pa_y[k] * gf_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pa_z, pb_y, ff0_0, ff1_0, gd_2, \
                         gf_0, gf_3, gf_6, gf_7, hd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * gf_0[k];

        t_11[k] = pa_z[k] * gf_3[k];

        t_12[k] = pb_y[k] * hd_5[k];

        t_13[k] = f_3 * gd_2[k]
                  + pa_z[k] * gf_6[k];

        t_14[k] = f_4 * ff0_0[k]
                  - f_5 * ff1_0[k]
                  + pa_y[k] * gf_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pb_x, pb_z, ff0_4, ff1_12, gd_7, gf_14, \
                         hd_6, hd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_z[k] * hd_6[k];

        t_16[k] = f_3 * gd_7[k]
                  + pb_x[k] * hd_7[k];

        t_17[k] = f_6 * ff0_4[k]
                  - f_7 * ff1_12[k]
                  + pa_x[k] * gf_14[k];

        t_18[k] = pb_z[k] * hd_7[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_z, ff0_0, ff1_0, gf_8, gf_9, \
                         gf_10, hp0_3, hp1_3, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_8[k];

        t_20[k] = pa_z[k] * gf_8[k];

        t_21[k] = pa_y[k] * gf_10[k];

        t_22[k] = f_4 * ff0_0[k]
                  - f_5 * ff1_0[k]
                  + pa_z[k] * gf_9[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_x, pb_y, gd_11, hp0_4, hp1_4, hd_9, hd_10, \
                         hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pb_y[k] * hd_9[k];

        t_24[k] = f_3 * gd_11[k]
                  + pb_x[k] * hd_11[k];

        t_25[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_10[k];

        t_26[k] = pb_y[k] * hd_11[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_y, pb_z, ff0_1, ff0_6, ff1_6, ff1_15, \
                         gf_11, gf_22, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_6 * ff0_6[k]
                  - f_7 * ff1_15[k]
                  + pa_x[k] * gf_22[k];

        t_28[k] = f_6 * ff0_1[k]
                  - f_7 * ff1_6[k]
                  + pa_y[k] * gf_11[k];

        t_29[k] = pb_z[k] * hd_12[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, pb_z, ff0_7, ff1_19, gd_12, \
                         gf_25, hp0_5, hp1_5, hd_13, hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_8 * gd_12[k]
                  + pb_x[k] * hd_13[k];

        t_31[k] = f_4 * ff0_7[k]
                  - f_5 * ff1_19[k]
                  + pa_x[k] * gf_25[k];

        t_32[k] = pb_z[k] * hd_13[k];

        t_33[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_14[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_y, pa_z, gd_8, gd_10, gf_11, gf_14, \
                         gf_16, gf_20, gf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_z[k] * gf_11[k];

        t_35[k] = pa_z[k] * gf_14[k];

        t_36[k] = f_3 * gd_8[k]
                  + pa_z[k] * gf_16[k];

        t_37[k] = f_3 * gd_10[k]
                  + pa_y[k] * gf_20[k];

        t_38[k] = pa_y[k] * gf_22[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pb_x, pb_y, ff0_2, ff1_8, gd_13, gf_17, \
                         hp0_6, hp1_6, hd_15, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_6 * ff0_2[k]
                  - f_7 * ff1_8[k]
                  + pa_z[k] * gf_17[k];

        t_40[k] = pb_y[k] * hd_15[k];

        t_41[k] = f_8 * gd_13[k]
                  + pb_x[k] * hd_17[k];

        t_42[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_16[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pb_x, pb_y, ff0_11, ff1_32, gd_14, \
                         gd_15, gf_28, gf_29, hd_17, hd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * hd_17[k];

        t_44[k] = f_4 * ff0_11[k]
                  - f_5 * ff1_32[k]
                  + pa_x[k] * gf_28[k];

        t_45[k] = f_3 * gd_14[k]
                  + pa_x[k] * gf_29[k];

        t_46[k] = f_9 * gd_15[k]
                  + pb_x[k] * hd_18[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pa_x, pa_z, pb_x, gd_19, gd_24, gd_26, \
                         gf_23, gf_33, gf_38, gf_47, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pa_x[k] * gf_33[k];

        t_48[k] = pa_z[k] * gf_23[k];

        t_49[k] = f_3 * gd_19[k]
                  + pa_x[k] * gf_38[k];

        t_50[k] = f_3 * gd_24[k]
                  + pa_x[k] * gf_47[k];

        t_51[k] = f_9 * gd_26[k]
                  + pb_x[k] * hd_19[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pa_x, pb_x, pb_z, gf_53, hp0_7, hp1_7, \
                         hd_20, hd_21, hd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pa_x[k] * gf_53[k];

        t_53[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_20[k];

        t_54[k] = pb_z[k] * hd_20[k];

        t_55[k] = pb_x[k] * hd_21[k];

        t_56[k] = pb_x[k] * hd_22[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_z, pb_y, pb_z, gd_15, gf_29, hp0_8, hp0_9, \
                         hp1_8, hp1_9, hd_21, hd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_0 * gd_15[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_21[k];

        t_58[k] = pb_z[k] * hd_21[k];

        t_59[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_9[k]
                  + pb_z[k] * hd_22[k];

        t_60[k] = pa_z[k] * gf_29[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pa_z, pb_x, gd_16, gf_33, gf_35, \
                         hp0_10, hp1_10, hd_24, hd_25, hd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_x[k] * hd_24[k];

        t_62[k] = pa_z[k] * gf_33[k];

        t_63[k] = f_3 * gd_16[k]
                  + pa_z[k] * gf_35[k];

        t_64[k] = f_1 * hp0_10[k]
                  - f_2 * hp1_10[k]
                  + pb_x[k] * hd_25[k];

        t_65[k] = pb_x[k] * hd_26[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_y, pa_z, pb_x, pb_y, ff0_7, ff0_10, \
                         ff1_19, ff1_26, gd_21, gf_36, gf_43, hd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pb_x[k] * hd_27[k];

        t_67[k] = f_4 * ff0_7[k]
                  - f_5 * ff1_19[k]
                  + pa_z[k] * gf_36[k];

        t_68[k] = f_3 * gd_21[k]
                  + pb_y[k] * hd_27[k];

        t_69[k] = f_6 * ff0_10[k]
                  - f_7 * ff1_26[k]
                  + pa_y[k] * gf_43[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_z, pb_x, ff0_8, ff1_22, gf_41, hp0_11, \
                         hp1_11, hd_28, hd_29, hd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_1 * hp0_11[k]
                  - f_2 * hp1_11[k]
                  + pb_x[k] * hd_28[k];

        t_71[k] = pb_x[k] * hd_29[k];

        t_72[k] = pb_x[k] * hd_30[k];

        t_73[k] = f_6 * ff0_8[k]
                  - f_7 * ff1_22[k]
                  + pa_z[k] * gf_41[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_y, pb_x, pb_y, ff0_11, ff1_32, gd_23, \
                         gd_25, gf_46, gf_51, hd_30, hd_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_8 * gd_23[k]
                  + pb_y[k] * hd_30[k];

        t_75[k] = f_4 * ff0_11[k]
                  - f_5 * ff1_32[k]
                  + pa_y[k] * gf_46[k];

        t_76[k] = pb_x[k] * hd_31[k];

        t_77[k] = f_3 * gd_25[k]
                  + pa_y[k] * gf_51[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pb_x, pb_y, gd_26, gf_53, hp0_12, \
                         hp1_12, hd_32, hd_33, hd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_9 * gd_26[k]
                  + pb_y[k] * hd_32[k];

        t_79[k] = pa_y[k] * gf_53[k];

        t_80[k] = f_1 * hp0_12[k]
                  - f_2 * hp1_12[k]
                  + pb_x[k] * hd_33[k];

        t_81[k] = pb_y[k] * hd_33[k];

        t_82[k] = pb_x[k] * hd_34[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pb_x, pb_y, pb_z, gd_26, hp0_13, hp0_14, \
                         hp1_13, hp1_14, hd_34, hd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pb_x[k] * hd_35[k];

        t_84[k] = f_1 * hp0_13[k]
                  - f_2 * hp1_13[k]
                  + pb_y[k] * hd_34[k];

        t_85[k] = pb_y[k] * hd_35[k];

        t_86[k] = f_0 * gd_26[k]
                  + f_1 * hp0_14[k]
                  - f_2 * hp1_14[k]
                  + pb_z[k] * hd_35[k];
    }
}

auto
compute_prim_hf_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t ff0, const size_t ff1,
                                      const size_t gd, const size_t gf, const size_t hp0,
                                      const size_t hp1, const size_t hd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_19 = buffer.data(ff0 + 19);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_26 = buffer.data(ff0 + 26);
    const auto *ff0_32 = buffer.data(ff0 + 32);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_9 = buffer.data(ff1 + 9);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_18 = buffer.data(ff1 + 18);
    const auto *ff1_20 = buffer.data(ff1 + 20);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_6 = buffer.data(gd + 6);
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
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_32 = buffer.data(gf + 32);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, hp0_0, hp0_1, hp1_0, \
                         hp1_1, hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_4[k] = pb_y[k] * hd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, gd_1, gd_2, gf_0, gf_3, \
                         gf_5, hp0_2, hp1_2, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_6[k] = pa_y[k] * gf_0[k];

        t_7[k] = f_3 * gd_1[k]
                 + pa_y[k] * gf_3[k];

        t_8[k] = pa_z[k] * gf_0[k];

        t_9[k] = f_3 * gd_2[k]
                 + pa_z[k] * gf_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pb_x, pb_z, ff0_0, ff1_0, gd_6, gf_6, hd_5, \
                         hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_4 * ff0_0[k]
                  - f_5 * ff1_0[k]
                  + pa_y[k] * gf_6[k];

        t_11[k] = pb_z[k] * hd_5[k];

        t_12[k] = f_3 * gd_6[k]
                  + pb_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pb_z, ff0_12, ff1_9, gf_10, hp0_3, hp1_3, \
                         hd_6, hd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_6 * ff0_12[k]
                  - f_7 * ff1_9[k]
                  + pa_x[k] * gf_10[k];

        t_14[k] = pb_z[k] * hd_6[k];

        t_15[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_z, pb_x, pb_y, ff0_0, ff1_0, gd_10, gf_7, \
                         hp0_4, hp1_4, hd_8, hd_9, hd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * ff0_0[k]
                  - f_5 * ff1_0[k]
                  + pa_z[k] * gf_7[k];

        t_17[k] = pb_y[k] * hd_8[k];

        t_18[k] = f_3 * gd_10[k]
                  + pb_x[k] * hd_10[k];

        t_19[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pa_y, pb_y, pb_z, ff0_6, ff0_15, ff1_6, \
                         ff1_11, gf_8, gf_13, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_y[k] * hd_10[k];

        t_21[k] = f_6 * ff0_15[k]
                  - f_7 * ff1_11[k]
                  + pa_x[k] * gf_13[k];

        t_22[k] = f_6 * ff0_6[k]
                  - f_7 * ff1_6[k]
                  + pa_y[k] * gf_8[k];

        t_23[k] = pb_z[k] * hd_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_x, pb_x, pb_z, ff0_19, ff1_15, gd_11, \
                         gf_14, hp0_5, hp1_5, hd_12, hd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * gd_11[k]
                  + pb_x[k] * hd_12[k];

        t_25[k] = f_4 * ff0_19[k]
                  - f_5 * ff1_15[k]
                  + pa_x[k] * gf_14[k];

        t_26[k] = pb_z[k] * hd_12[k];

        t_27[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pa_z, pb_x, pb_y, ff0_8, ff1_7, gd_12, \
                         gf_11, hd_14, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_y[k] * gf_11[k];

        t_29[k] = f_6 * ff0_8[k]
                  - f_7 * ff1_7[k]
                  + pa_z[k] * gf_11[k];

        t_30[k] = pb_y[k] * hd_14[k];

        t_31[k] = f_8 * gd_12[k]
                  + pb_x[k] * hd_16[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_y, ff0_32, ff1_26, gd_13, gf_15, \
                         gf_16, hp0_6, hp1_6, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_15[k];

        t_33[k] = pb_y[k] * hd_16[k];

        t_34[k] = f_4 * ff0_32[k]
                  - f_5 * ff1_26[k]
                  + pa_x[k] * gf_15[k];

        t_35[k] = f_3 * gd_13[k]
                  + pa_x[k] * gf_16[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pa_x, gd_21, gf_19, gf_23, gf_25, \
                         gf_27, gf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_x[k] * gf_19[k];

        t_37[k] = pa_x[k] * gf_23[k];

        t_38[k] = pa_x[k] * gf_25[k];

        t_39[k] = f_3 * gd_21[k]
                  + pa_x[k] * gf_27[k];

        t_40[k] = pa_x[k] * gf_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pb_x, pb_y, pb_z, gd_14, hp0_7, hp0_8, \
                         hp1_7, hp1_8, hd_19, hd_20, hd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_19[k];

        t_42[k] = pb_x[k] * hd_20[k];

        t_43[k] = pb_x[k] * hd_21[k];

        t_44[k] = f_0 * gd_14[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_20[k];

        t_45[k] = pb_z[k] * hd_20[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_z, pb_x, pb_z, gd_15, gf_19, gf_21, hp0_9, \
                         hp0_10, hp1_9, hp1_10, hd_21, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_9[k]
                  + pb_z[k] * hd_21[k];

        t_47[k] = pa_z[k] * gf_19[k];

        t_48[k] = f_3 * gd_15[k]
                  + pa_z[k] * gf_21[k];

        t_49[k] = f_1 * hp0_10[k]
                  - f_2 * hp1_10[k]
                  + pb_x[k] * hd_23[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_z, pb_x, pb_y, ff0_19, ff1_15, gd_19, \
                         gf_22, hd_24, hd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pb_x[k] * hd_24[k];

        t_51[k] = pb_x[k] * hd_25[k];

        t_52[k] = f_4 * ff0_19[k]
                  - f_5 * ff1_15[k]
                  + pa_z[k] * gf_22[k];

        t_53[k] = f_3 * gd_19[k]
                  + pb_y[k] * hd_25[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pb_x, ff0_26, ff1_20, gf_25, hp0_11, \
                         hp1_11, hd_26, hd_27, hd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_6 * ff0_26[k]
                  - f_7 * ff1_20[k]
                  + pa_y[k] * gf_25[k];

        t_55[k] = f_1 * hp0_11[k]
                  - f_2 * hp1_11[k]
                  + pb_x[k] * hd_26[k];

        t_56[k] = pb_x[k] * hd_27[k];

        t_57[k] = pb_x[k] * hd_28[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_y, pa_z, pb_y, ff0_22, ff0_32, ff1_18, ff1_26, \
                         gd_20, gf_23, gf_26, hd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_6 * ff0_22[k]
                  - f_7 * ff1_18[k]
                  + pa_z[k] * gf_23[k];

        t_59[k] = f_8 * gd_20[k]
                  + pb_y[k] * hd_28[k];

        t_60[k] = f_4 * ff0_32[k]
                  - f_5 * ff1_26[k]
                  + pa_y[k] * gf_26[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pa_y, pb_x, gd_22, gf_30, gf_32, \
                         hp0_12, hp1_12, hd_30, hd_31, hd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * gd_22[k]
                  + pa_y[k] * gf_30[k];

        t_62[k] = pa_y[k] * gf_32[k];

        t_63[k] = f_1 * hp0_12[k]
                  - f_2 * hp1_12[k]
                  + pb_x[k] * hd_30[k];

        t_64[k] = pb_x[k] * hd_31[k];

        t_65[k] = pb_x[k] * hd_32[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_y, pb_z, gd_23, hp0_13, hp0_14, hp1_13, hp1_14, \
                         hd_31, hd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_1 * hp0_13[k]
                  - f_2 * hp1_13[k]
                  + pb_y[k] * hd_31[k];

        t_67[k] = pb_y[k] * hd_32[k];

        t_68[k] = f_0 * gd_23[k]
                  + f_1 * hp0_14[k]
                  - f_2 * hp1_14[k]
                  + pb_z[k] * hd_32[k];
    }
}

auto
compute_prim_hf_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t ff0, const size_t ff1,
                                      const size_t gd, const size_t gf, const size_t hp0,
                                      const size_t hp1, const size_t hd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / p;
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
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

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_1 = buffer.data(ff1 + 1);
    const auto *ff1_2 = buffer.data(ff1 + 2);
    const auto *ff1_4 = buffer.data(ff1 + 4);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_10 = buffer.data(ff1 + 10);
    const auto *ff1_11 = buffer.data(ff1 + 11);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_44 = buffer.data(gf + 44);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, hp0_0, hp0_1, hp1_0, \
                         hp1_1, hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_4[k] = pb_y[k] * hd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, ff0_0, ff1_0, gd_6, gf_6, \
                         hp0_2, hp1_2, hd_2, hd_5, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_6[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pa_y[k] * gf_6[k];

        t_7[k] = pb_z[k] * hd_5[k];

        t_8[k] = f_5 * gd_6[k]
                 + pb_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, ff0_4, ff1_4, gf_11, hp0_3, hp1_3, hd_6, \
                         hd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * ff0_4[k]
                 - f_7 * ff1_4[k]
                 + pa_x[k] * gf_11[k];

        t_10[k] = pb_z[k] * hd_6[k];

        t_11[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_7[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, ff0_0, ff1_0, gd_10, gf_7, \
                         hp0_4, hp1_4, hd_8, hd_9, hd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * ff0_0[k]
                  - f_4 * ff1_0[k]
                  + pa_z[k] * gf_7[k];

        t_13[k] = pb_y[k] * hd_8[k];

        t_14[k] = f_5 * gd_10[k]
                  + pb_x[k] * hd_10[k];

        t_15[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, ff0_1, ff0_6, ff1_1, \
                         ff1_6, gf_8, gf_19, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * hd_10[k];

        t_17[k] = f_6 * ff0_6[k]
                  - f_7 * ff1_6[k]
                  + pa_x[k] * gf_19[k];

        t_18[k] = f_6 * ff0_1[k]
                  - f_7 * ff1_1[k]
                  + pa_y[k] * gf_8[k];

        t_19[k] = pb_z[k] * hd_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, ff0_7, ff1_7, gd_11, gf_21, \
                         hp0_5, hp1_5, hd_12, hd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_8 * gd_11[k]
                  + pb_x[k] * hd_12[k];

        t_21[k] = f_3 * ff0_7[k]
                  - f_4 * ff1_7[k]
                  + pa_x[k] * gf_21[k];

        t_22[k] = pb_z[k] * hd_12[k];

        t_23[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_13[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, ff0_2, ff1_2, gd_12, gf_14, \
                         hp0_6, hp1_6, hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_6 * ff0_2[k]
                  - f_7 * ff1_2[k]
                  + pa_z[k] * gf_14[k];

        t_25[k] = pb_y[k] * hd_14[k];

        t_26[k] = f_8 * gd_12[k]
                  + pb_x[k] * hd_16[k];

        t_27[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_15[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_x, pb_y, ff0_11, ff1_11, gd_14, \
                         gf_23, gf_27, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * hd_16[k];

        t_29[k] = f_3 * ff0_11[k]
                  - f_4 * ff1_11[k]
                  + pa_x[k] * gf_23[k];

        t_30[k] = f_9 * gd_14[k]
                  + pb_x[k] * hd_17[k];

        t_31[k] = pa_x[k] * gf_27[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, pa_x, pb_x, gd_23, gf_44, hp0_7, hp1_7, \
                         hd_18, hd_19, hd_20, hd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_9 * gd_23[k]
                  + pb_x[k] * hd_18[k];

        t_33[k] = pa_x[k] * gf_44[k];

        t_34[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_19[k];

        t_35[k] = pb_x[k] * hd_20[k];

        t_36[k] = pb_x[k] * hd_21[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_y, pb_z, gd_14, hp0_8, hp0_9, hp1_8, hp1_9, \
                         hd_20, hd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_0 * gd_14[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_20[k];

        t_38[k] = pb_z[k] * hd_20[k];

        t_39[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_9[k]
                  + pb_z[k] * hd_21[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pb_x, ff0_7, ff1_7, gf_30, hp0_10, \
                         hp1_10, hd_23, hd_24, hd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_1 * hp0_10[k]
                  - f_2 * hp1_10[k]
                  + pb_x[k] * hd_23[k];

        t_41[k] = pb_x[k] * hd_24[k];

        t_42[k] = pb_x[k] * hd_25[k];

        t_43[k] = f_3 * ff0_7[k]
                  - f_4 * ff1_7[k]
                  + pa_z[k] * gf_30[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_y, pb_x, pb_y, ff0_10, ff1_10, gd_19, \
                         gf_36, hp0_11, hp1_11, hd_25, hd_26, hd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_5 * gd_19[k]
                  + pb_y[k] * hd_25[k];

        t_45[k] = f_6 * ff0_10[k]
                  - f_7 * ff1_10[k]
                  + pa_y[k] * gf_36[k];

        t_46[k] = f_1 * hp0_11[k]
                  - f_2 * hp1_11[k]
                  + pb_x[k] * hd_26[k];

        t_47[k] = pb_x[k] * hd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pa_z, pb_x, pb_y, ff0_8, ff0_11, ff1_8, \
                         ff1_11, gd_20, gf_34, gf_38, hd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pb_x[k] * hd_28[k];

        t_49[k] = f_6 * ff0_8[k]
                  - f_7 * ff1_8[k]
                  + pa_z[k] * gf_34[k];

        t_50[k] = f_8 * gd_20[k]
                  + pb_y[k] * hd_28[k];

        t_51[k] = f_3 * ff0_11[k]
                  - f_4 * ff1_11[k]
                  + pa_y[k] * gf_38[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pa_y, pb_x, pb_y, gd_23, gf_44, hp0_12, \
                         hp1_12, hd_29, hd_30, hd_31, hd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_9 * gd_23[k]
                  + pb_y[k] * hd_29[k];

        t_53[k] = pa_y[k] * gf_44[k];

        t_54[k] = f_1 * hp0_12[k]
                  - f_2 * hp1_12[k]
                  + pb_x[k] * hd_30[k];

        t_55[k] = pb_x[k] * hd_31[k];

        t_56[k] = pb_x[k] * hd_32[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_y, pb_z, gd_23, hp0_13, hp0_14, hp1_13, hp1_14, \
                         hd_31, hd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_1 * hp0_13[k]
                  - f_2 * hp1_13[k]
                  + pb_y[k] * hd_31[k];

        t_58[k] = pb_y[k] * hd_32[k];

        t_59[k] = f_0 * gd_23[k]
                  + f_1 * hp0_14[k]
                  - f_2 * hp1_14[k]
                  + pb_z[k] * hd_32[k];
    }
}

auto
compute_prim_hf_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t ff0, const size_t ff1,
                                      const size_t gd, const size_t gf, const size_t hp0,
                                      const size_t hp1, const size_t hd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / p;
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
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

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_1 = buffer.data(ff0 + 1);
    const auto *ff0_2 = buffer.data(ff0 + 2);
    const auto *ff0_4 = buffer.data(ff0 + 4);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_11 = buffer.data(ff0 + 11);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_8 = buffer.data(ff1 + 8);
    const auto *ff1_10 = buffer.data(ff1 + 10);
    const auto *ff1_12 = buffer.data(ff1 + 12);
    const auto *ff1_17 = buffer.data(ff1 + 17);
    const auto *ff1_20 = buffer.data(ff1 + 20);
    const auto *ff1_22 = buffer.data(ff1 + 22);
    const auto *ff1_29 = buffer.data(ff1 + 29);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_50 = buffer.data(gf + 50);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, hp0_0, hp0_1, hp1_0, \
                         hp1_1, hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_4[k] = pb_y[k] * hd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, ff0_0, ff1_0, gf_0, gf_7, \
                         hp0_2, hp1_2, hd_2, hd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_6[k] = pa_y[k] * gf_0[k];

        t_7[k] = pa_z[k] * gf_0[k];

        t_8[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pa_y[k] * gf_7[k];

        t_9[k] = pb_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, ff0_4, ff1_10, gd_6, gf_13, \
                         hp0_3, hp1_3, hd_6, hd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gd_6[k]
                  + pb_x[k] * hd_6[k];

        t_11[k] = f_6 * ff0_4[k]
                  - f_7 * ff1_10[k]
                  + pa_x[k] * gf_13[k];

        t_12[k] = pb_z[k] * hd_6[k];

        t_13[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, ff0_0, ff1_0, gd_10, gf_8, \
                         hp0_4, hp1_4, hd_8, hd_9, hd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * ff0_0[k]
                  - f_4 * ff1_0[k]
                  + pa_z[k] * gf_8[k];

        t_15[k] = pb_y[k] * hd_8[k];

        t_16[k] = f_5 * gd_10[k]
                  + pb_x[k] * hd_10[k];

        t_17[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, ff0_1, ff0_6, ff1_7, \
                         ff1_12, gf_10, gf_21, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * hd_10[k];

        t_19[k] = f_6 * ff0_6[k]
                  - f_7 * ff1_12[k]
                  + pa_x[k] * gf_21[k];

        t_20[k] = f_6 * ff0_1[k]
                  - f_7 * ff1_7[k]
                  + pa_y[k] * gf_10[k];

        t_21[k] = pb_z[k] * hd_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, ff0_7, ff1_17, gd_11, \
                         gf_23, hp0_5, hp1_5, hd_12, hd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_8 * gd_11[k]
                  + pb_x[k] * hd_12[k];

        t_23[k] = f_3 * ff0_7[k]
                  - f_4 * ff1_17[k]
                  + pa_x[k] * gf_23[k];

        t_24[k] = pb_z[k] * hd_12[k];

        t_25[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, ff0_2, ff1_8, gd_12, gf_16, \
                         hp0_6, hp1_6, hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_6 * ff0_2[k]
                  - f_7 * ff1_8[k]
                  + pa_z[k] * gf_16[k];

        t_27[k] = pb_y[k] * hd_14[k];

        t_28[k] = f_8 * gd_12[k]
                  + pb_x[k] * hd_16[k];

        t_29[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, pb_y, ff0_11, ff1_29, gd_14, \
                         gf_25, gf_30, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * hd_16[k];

        t_31[k] = f_3 * ff0_11[k]
                  - f_4 * ff1_29[k]
                  + pa_x[k] * gf_25[k];

        t_32[k] = f_9 * gd_14[k]
                  + pb_x[k] * hd_17[k];

        t_33[k] = pa_x[k] * gf_30[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_x, pb_x, gd_23, gf_50, hp0_7, hp1_7, \
                         hd_18, hd_19, hd_20, hd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_9 * gd_23[k]
                  + pb_x[k] * hd_18[k];

        t_35[k] = pa_x[k] * gf_50[k];

        t_36[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_19[k];

        t_37[k] = pb_x[k] * hd_20[k];

        t_38[k] = pb_x[k] * hd_21[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pb_y, pb_z, gd_14, gf_30, hp0_8, hp0_9, \
                         hp1_8, hp1_9, hd_20, hd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * gd_14[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_20[k];

        t_40[k] = pb_z[k] * hd_20[k];

        t_41[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_9[k]
                  + pb_z[k] * hd_21[k];

        t_42[k] = pa_z[k] * gf_30[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_z, pb_x, ff0_7, ff1_17, gf_33, hp0_10, \
                         hp1_10, hd_23, hd_24, hd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * hp0_10[k]
                  - f_2 * hp1_10[k]
                  + pb_x[k] * hd_23[k];

        t_44[k] = pb_x[k] * hd_24[k];

        t_45[k] = pb_x[k] * hd_25[k];

        t_46[k] = f_3 * ff0_7[k]
                  - f_4 * ff1_17[k]
                  + pa_z[k] * gf_33[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pb_x, pb_y, ff0_10, ff1_22, gd_19, \
                         gf_40, hp0_11, hp1_11, hd_25, hd_26, hd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_5 * gd_19[k]
                  + pb_y[k] * hd_25[k];

        t_48[k] = f_6 * ff0_10[k]
                  - f_7 * ff1_22[k]
                  + pa_y[k] * gf_40[k];

        t_49[k] = f_1 * hp0_11[k]
                  - f_2 * hp1_11[k]
                  + pb_x[k] * hd_26[k];

        t_50[k] = pb_x[k] * hd_27[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pa_z, pb_x, pb_y, ff0_8, ff0_11, \
                         ff1_20, ff1_29, gd_20, gf_38, gf_43, hd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_x[k] * hd_28[k];

        t_52[k] = f_6 * ff0_8[k]
                  - f_7 * ff1_20[k]
                  + pa_z[k] * gf_38[k];

        t_53[k] = f_8 * gd_20[k]
                  + pb_y[k] * hd_28[k];

        t_54[k] = f_3 * ff0_11[k]
                  - f_4 * ff1_29[k]
                  + pa_y[k] * gf_43[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_y, pb_x, pb_y, gd_23, gf_50, hp0_12, \
                         hp1_12, hd_29, hd_30, hd_31, hd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_9 * gd_23[k]
                  + pb_y[k] * hd_29[k];

        t_56[k] = pa_y[k] * gf_50[k];

        t_57[k] = f_1 * hp0_12[k]
                  - f_2 * hp1_12[k]
                  + pb_x[k] * hd_30[k];

        t_58[k] = pb_x[k] * hd_31[k];

        t_59[k] = pb_x[k] * hd_32[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pb_y, pb_z, gd_23, hp0_13, hp0_14, hp1_13, hp1_14, \
                         hd_31, hd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * hp0_13[k]
                  - f_2 * hp1_13[k]
                  + pb_y[k] * hd_31[k];

        t_61[k] = pb_y[k] * hd_32[k];

        t_62[k] = f_0 * gd_23[k]
                  + f_1 * hp0_14[k]
                  - f_2 * hp1_14[k]
                  + pb_z[k] * hd_32[k];
    }
}

auto
compute_prim_hf_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t ff0, const size_t ff1,
                                      const size_t gd, const size_t gf, const size_t hp0,
                                      const size_t hp1, const size_t hd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_8 = buffer.data(ff0 + 8);
    const auto *ff0_10 = buffer.data(ff0 + 10);
    const auto *ff0_12 = buffer.data(ff0 + 12);
    const auto *ff0_17 = buffer.data(ff0 + 17);
    const auto *ff0_20 = buffer.data(ff0 + 20);
    const auto *ff0_22 = buffer.data(ff0 + 22);
    const auto *ff0_29 = buffer.data(ff0 + 29);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_9 = buffer.data(ff1 + 9);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_18 = buffer.data(ff1 + 18);
    const auto *ff1_20 = buffer.data(ff1 + 20);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_44 = buffer.data(gf + 44);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, hp0_0, hp0_1, hp1_0, \
                         hp1_1, hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_4[k] = pb_y[k] * hd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pa_z, pb_z, gd_2, gf_0, gf_5, hp0_2, hp1_2, \
                         hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_6[k] = pa_y[k] * gf_0[k];

        t_7[k] = pa_z[k] * gf_0[k];

        t_8[k] = f_3 * gd_2[k]
                 + pa_z[k] * gf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pb_z, ff0_0, ff1_0, gd_6, gf_6, hd_5, \
                         hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * ff0_0[k]
                 - f_5 * ff1_0[k]
                 + pa_y[k] * gf_6[k];

        t_10[k] = pb_z[k] * hd_5[k];

        t_11[k] = f_3 * gd_6[k]
                  + pb_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_z, ff0_10, ff1_9, gf_11, hp0_3, hp1_3, \
                         hd_6, hd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_6 * ff0_10[k]
                  - f_7 * ff1_9[k]
                  + pa_x[k] * gf_11[k];

        t_13[k] = pb_z[k] * hd_6[k];

        t_14[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_z, pb_x, pb_y, ff0_0, ff1_0, gd_10, gf_7, \
                         hp0_4, hp1_4, hd_8, hd_9, hd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * ff0_0[k]
                  - f_5 * ff1_0[k]
                  + pa_z[k] * gf_7[k];

        t_16[k] = pb_y[k] * hd_8[k];

        t_17[k] = f_3 * gd_10[k]
                  + pb_x[k] * hd_10[k];

        t_18[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pa_y, pb_y, pb_z, ff0_7, ff0_12, ff1_6, \
                         ff1_11, gf_8, gf_19, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * hd_10[k];

        t_20[k] = f_6 * ff0_12[k]
                  - f_7 * ff1_11[k]
                  + pa_x[k] * gf_19[k];

        t_21[k] = f_6 * ff0_7[k]
                  - f_7 * ff1_6[k]
                  + pa_y[k] * gf_8[k];

        t_22[k] = pb_z[k] * hd_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_x, pb_z, ff0_17, ff1_15, gd_11, \
                         gf_21, hp0_5, hp1_5, hd_12, hd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_8 * gd_11[k]
                  + pb_x[k] * hd_12[k];

        t_24[k] = f_4 * ff0_17[k]
                  - f_5 * ff1_15[k]
                  + pa_x[k] * gf_21[k];

        t_25[k] = pb_z[k] * hd_12[k];

        t_26[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_13[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_z, pb_x, pb_y, ff0_8, ff1_7, gd_12, gf_14, \
                         hp0_6, hp1_6, hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_6 * ff0_8[k]
                  - f_7 * ff1_7[k]
                  + pa_z[k] * gf_14[k];

        t_28[k] = pb_y[k] * hd_14[k];

        t_29[k] = f_8 * gd_12[k]
                  + pb_x[k] * hd_16[k];

        t_30[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_15[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pb_x, pb_y, ff0_29, ff1_26, gd_14, \
                         gf_23, gf_27, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * hd_16[k];

        t_32[k] = f_4 * ff0_29[k]
                  - f_5 * ff1_26[k]
                  + pa_x[k] * gf_23[k];

        t_33[k] = f_9 * gd_14[k]
                  + pb_x[k] * hd_17[k];

        t_34[k] = pa_x[k] * gf_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_x, pb_x, gd_23, gf_44, hp0_7, hp1_7, \
                         hd_18, hd_19, hd_20, hd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_9 * gd_23[k]
                  + pb_x[k] * hd_18[k];

        t_36[k] = pa_x[k] * gf_44[k];

        t_37[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_19[k];

        t_38[k] = pb_x[k] * hd_20[k];

        t_39[k] = pb_x[k] * hd_21[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pb_y, pb_z, gd_14, gf_27, hp0_8, hp0_9, \
                         hp1_8, hp1_9, hd_20, hd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * gd_14[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_20[k];

        t_41[k] = pb_z[k] * hd_20[k];

        t_42[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_9[k]
                  + pb_z[k] * hd_21[k];

        t_43[k] = pa_z[k] * gf_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pb_x, gd_15, gf_29, hp0_10, hp1_10, \
                         hd_23, hd_24, hd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_3 * gd_15[k]
                  + pa_z[k] * gf_29[k];

        t_45[k] = f_1 * hp0_10[k]
                  - f_2 * hp1_10[k]
                  + pb_x[k] * hd_23[k];

        t_46[k] = pb_x[k] * hd_24[k];

        t_47[k] = pb_x[k] * hd_25[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pa_y, pa_z, pb_y, ff0_17, ff0_22, ff1_15, ff1_20, \
                         gd_19, gf_30, gf_36, hd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_4 * ff0_17[k]
                  - f_5 * ff1_15[k]
                  + pa_z[k] * gf_30[k];

        t_49[k] = f_3 * gd_19[k]
                  + pb_y[k] * hd_25[k];

        t_50[k] = f_6 * ff0_22[k]
                  - f_7 * ff1_20[k]
                  + pa_y[k] * gf_36[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_x, ff0_20, ff1_18, gf_34, hp0_11, \
                         hp1_11, hd_26, hd_27, hd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * hp0_11[k]
                  - f_2 * hp1_11[k]
                  + pb_x[k] * hd_26[k];

        t_52[k] = pb_x[k] * hd_27[k];

        t_53[k] = pb_x[k] * hd_28[k];

        t_54[k] = f_6 * ff0_20[k]
                  - f_7 * ff1_18[k]
                  + pa_z[k] * gf_34[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pb_y, ff0_29, ff1_26, gd_20, gd_22, \
                         gd_23, gf_38, gf_42, hd_28, hd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_8 * gd_20[k]
                  + pb_y[k] * hd_28[k];

        t_56[k] = f_4 * ff0_29[k]
                  - f_5 * ff1_26[k]
                  + pa_y[k] * gf_38[k];

        t_57[k] = f_3 * gd_22[k]
                  + pa_y[k] * gf_42[k];

        t_58[k] = f_9 * gd_23[k]
                  + pb_y[k] * hd_29[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, t_63, pa_y, pb_x, pb_y, gf_44, hp0_12, \
                         hp0_13, hp1_12, hp1_13, hd_30, hd_31, hd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_y[k] * gf_44[k];

        t_60[k] = f_1 * hp0_12[k]
                  - f_2 * hp1_12[k]
                  + pb_x[k] * hd_30[k];

        t_61[k] = pb_x[k] * hd_31[k];

        t_62[k] = pb_x[k] * hd_32[k];

        t_63[k] = f_1 * hp0_13[k]
                  - f_2 * hp1_13[k]
                  + pb_y[k] * hd_31[k];
    }

#pragma omp simd aligned(t_64, t_65, pb_y, pb_z, gd_23, hp0_14, hp1_14, \
                         hd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_y[k] * hd_32[k];

        t_65[k] = f_0 * gd_23[k]
                  + f_1 * hp0_14[k]
                  - f_2 * hp1_14[k]
                  + pb_z[k] * hd_32[k];
    }
}

auto
compute_prim_hf_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t ff0, const size_t ff1,
                                      const size_t gd, const size_t gf, const size_t hp0,
                                      const size_t hp1, const size_t hd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / p;
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 1.0 / p;

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

    const auto *ff0_0 = buffer.data(ff0 + 0);
    const auto *ff0_6 = buffer.data(ff0 + 6);
    const auto *ff0_7 = buffer.data(ff0 + 7);
    const auto *ff0_9 = buffer.data(ff0 + 9);
    const auto *ff0_11 = buffer.data(ff0 + 11);
    const auto *ff0_15 = buffer.data(ff0 + 15);
    const auto *ff0_18 = buffer.data(ff0 + 18);
    const auto *ff0_20 = buffer.data(ff0 + 20);
    const auto *ff0_26 = buffer.data(ff0 + 26);

    const auto *ff1_0 = buffer.data(ff1 + 0);
    const auto *ff1_6 = buffer.data(ff1 + 6);
    const auto *ff1_7 = buffer.data(ff1 + 7);
    const auto *ff1_9 = buffer.data(ff1 + 9);
    const auto *ff1_11 = buffer.data(ff1 + 11);
    const auto *ff1_15 = buffer.data(ff1 + 15);
    const auto *ff1_18 = buffer.data(ff1 + 18);
    const auto *ff1_20 = buffer.data(ff1 + 20);
    const auto *ff1_26 = buffer.data(ff1 + 26);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_32 = buffer.data(gf + 32);

    const auto *hp0_0 = buffer.data(hp0 + 0);
    const auto *hp0_1 = buffer.data(hp0 + 1);
    const auto *hp0_2 = buffer.data(hp0 + 2);
    const auto *hp0_3 = buffer.data(hp0 + 3);
    const auto *hp0_4 = buffer.data(hp0 + 4);
    const auto *hp0_5 = buffer.data(hp0 + 5);
    const auto *hp0_6 = buffer.data(hp0 + 6);
    const auto *hp0_7 = buffer.data(hp0 + 7);
    const auto *hp0_8 = buffer.data(hp0 + 8);
    const auto *hp0_9 = buffer.data(hp0 + 9);
    const auto *hp0_10 = buffer.data(hp0 + 10);
    const auto *hp0_11 = buffer.data(hp0 + 11);
    const auto *hp0_12 = buffer.data(hp0 + 12);
    const auto *hp0_13 = buffer.data(hp0 + 13);
    const auto *hp0_14 = buffer.data(hp0 + 14);

    const auto *hp1_0 = buffer.data(hp1 + 0);
    const auto *hp1_1 = buffer.data(hp1 + 1);
    const auto *hp1_2 = buffer.data(hp1 + 2);
    const auto *hp1_3 = buffer.data(hp1 + 3);
    const auto *hp1_4 = buffer.data(hp1 + 4);
    const auto *hp1_5 = buffer.data(hp1 + 5);
    const auto *hp1_6 = buffer.data(hp1 + 6);
    const auto *hp1_7 = buffer.data(hp1 + 7);
    const auto *hp1_8 = buffer.data(hp1 + 8);
    const auto *hp1_9 = buffer.data(hp1 + 9);
    const auto *hp1_10 = buffer.data(hp1 + 10);
    const auto *hp1_11 = buffer.data(hp1 + 11);
    const auto *hp1_12 = buffer.data(hp1 + 12);
    const auto *hp1_13 = buffer.data(hp1 + 13);
    const auto *hp1_14 = buffer.data(hp1 + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gd_0, hp0_0, hp0_1, hp1_0, \
                         hp1_1, hd_0, hd_1, hd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 + f_1 * hp0_0[k]
                 - f_2 * hp1_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = pb_y[k] * hd_0[k];

        t_2[k] = pb_z[k] * hd_0[k];

        t_3[k] = f_1 * hp0_1[k]
                 - f_2 * hp1_1[k]
                 + pb_y[k] * hd_1[k];

        t_4[k] = pb_y[k] * hd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, ff0_0, ff1_0, gf_0, gf_6, \
                         hp0_2, hp1_2, hd_2, hd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hp0_2[k]
                 - f_2 * hp1_2[k]
                 + pb_z[k] * hd_2[k];

        t_6[k] = pa_y[k] * gf_0[k];

        t_7[k] = pa_z[k] * gf_0[k];

        t_8[k] = f_3 * ff0_0[k]
                 - f_4 * ff1_0[k]
                 + pa_y[k] * gf_6[k];

        t_9[k] = pb_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, ff0_9, ff1_9, gd_6, gf_10, \
                         hp0_3, hp1_3, hd_6, hd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gd_6[k]
                  + pb_x[k] * hd_6[k];

        t_11[k] = f_6 * ff0_9[k]
                  - f_7 * ff1_9[k]
                  + pa_x[k] * gf_10[k];

        t_12[k] = pb_z[k] * hd_6[k];

        t_13[k] = f_1 * hp0_3[k]
                  - f_2 * hp1_3[k]
                  + pb_z[k] * hd_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, ff0_0, ff1_0, gd_10, gf_7, \
                         hp0_4, hp1_4, hd_8, hd_9, hd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * ff0_0[k]
                  - f_4 * ff1_0[k]
                  + pa_z[k] * gf_7[k];

        t_15[k] = pb_y[k] * hd_8[k];

        t_16[k] = f_5 * gd_10[k]
                  + pb_x[k] * hd_10[k];

        t_17[k] = f_1 * hp0_4[k]
                  - f_2 * hp1_4[k]
                  + pb_y[k] * hd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, ff0_6, ff0_11, ff1_6, \
                         ff1_11, gf_8, gf_13, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * hd_10[k];

        t_19[k] = f_6 * ff0_11[k]
                  - f_7 * ff1_11[k]
                  + pa_x[k] * gf_13[k];

        t_20[k] = f_6 * ff0_6[k]
                  - f_7 * ff1_6[k]
                  + pa_y[k] * gf_8[k];

        t_21[k] = pb_z[k] * hd_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, ff0_15, ff1_15, gd_11, \
                         gf_14, hp0_5, hp1_5, hd_12, hd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_8 * gd_11[k]
                  + pb_x[k] * hd_12[k];

        t_23[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pa_x[k] * gf_14[k];

        t_24[k] = pb_z[k] * hd_12[k];

        t_25[k] = f_1 * hp0_5[k]
                  - f_2 * hp1_5[k]
                  + pb_z[k] * hd_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, ff0_7, ff1_7, gd_12, gf_11, \
                         hp0_6, hp1_6, hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_6 * ff0_7[k]
                  - f_7 * ff1_7[k]
                  + pa_z[k] * gf_11[k];

        t_27[k] = pb_y[k] * hd_14[k];

        t_28[k] = f_8 * gd_12[k]
                  + pb_x[k] * hd_16[k];

        t_29[k] = f_1 * hp0_6[k]
                  - f_2 * hp1_6[k]
                  + pb_y[k] * hd_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_y, ff0_26, ff1_26, gf_15, gf_19, \
                         gf_32, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * hd_16[k];

        t_31[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pa_x[k] * gf_15[k];

        t_32[k] = pa_x[k] * gf_19[k];

        t_33[k] = pa_x[k] * gf_32[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, gd_14, hp0_7, hp0_8, \
                         hp1_7, hp1_8, hd_19, hd_20, hd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_1 * hp0_7[k]
                  - f_2 * hp1_7[k]
                  + pb_x[k] * hd_19[k];

        t_35[k] = pb_x[k] * hd_20[k];

        t_36[k] = pb_x[k] * hd_21[k];

        t_37[k] = f_0 * gd_14[k]
                  + f_1 * hp0_8[k]
                  - f_2 * hp1_8[k]
                  + pb_y[k] * hd_20[k];

        t_38[k] = pb_z[k] * hd_20[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pb_x, pb_z, gf_19, hp0_9, hp0_10, \
                         hp1_9, hp1_10, hd_21, hd_23, hd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_1 * hp0_9[k]
                  - f_2 * hp1_9[k]
                  + pb_z[k] * hd_21[k];

        t_40[k] = pa_z[k] * gf_19[k];

        t_41[k] = f_1 * hp0_10[k]
                  - f_2 * hp1_10[k]
                  + pb_x[k] * hd_23[k];

        t_42[k] = pb_x[k] * hd_24[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_y, pa_z, pb_x, pb_y, ff0_15, ff0_20, \
                         ff1_15, ff1_20, gd_19, gf_22, gf_25, hd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_x[k] * hd_25[k];

        t_44[k] = f_3 * ff0_15[k]
                  - f_4 * ff1_15[k]
                  + pa_z[k] * gf_22[k];

        t_45[k] = f_5 * gd_19[k]
                  + pb_y[k] * hd_25[k];

        t_46[k] = f_6 * ff0_20[k]
                  - f_7 * ff1_20[k]
                  + pa_y[k] * gf_25[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_z, pb_x, ff0_18, ff1_18, gf_23, hp0_11, \
                         hp1_11, hd_26, hd_27, hd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_1 * hp0_11[k]
                  - f_2 * hp1_11[k]
                  + pb_x[k] * hd_26[k];

        t_48[k] = pb_x[k] * hd_27[k];

        t_49[k] = pb_x[k] * hd_28[k];

        t_50[k] = f_6 * ff0_18[k]
                  - f_7 * ff1_18[k]
                  + pa_z[k] * gf_23[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pb_x, pb_y, ff0_26, ff1_26, gd_20, \
                         gf_26, gf_32, hp0_12, hp1_12, hd_28, hd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_8 * gd_20[k]
                  + pb_y[k] * hd_28[k];

        t_52[k] = f_3 * ff0_26[k]
                  - f_4 * ff1_26[k]
                  + pa_y[k] * gf_26[k];

        t_53[k] = pa_y[k] * gf_32[k];

        t_54[k] = f_1 * hp0_12[k]
                  - f_2 * hp1_12[k]
                  + pb_x[k] * hd_30[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pb_x, pb_y, pb_z, gd_23, hp0_13, \
                         hp0_14, hp1_13, hp1_14, hd_31, hd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pb_x[k] * hd_31[k];

        t_56[k] = pb_x[k] * hd_32[k];

        t_57[k] = f_1 * hp0_13[k]
                  - f_2 * hp1_13[k]
                  + pb_y[k] * hd_31[k];

        t_58[k] = pb_y[k] * hd_32[k];

        t_59[k] = f_0 * gd_23[k]
                  + f_1 * hp0_14[k]
                  - f_2 * hp1_14[k]
                  + pb_z[k] * hd_32[k];
    }
}

}  // namespace simdt2ceri
