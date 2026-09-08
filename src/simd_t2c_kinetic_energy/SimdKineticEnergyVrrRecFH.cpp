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


#include "SimdKineticEnergyVrrRecFH.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_fh_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t ph_s, const size_t ph,
                                 const size_t dg, const size_t dh, const size_t ff_s,
                                 const size_t fh_s, const size_t ff, const size_t fg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = beta / p;
    const auto f_9 = 3.0 * alpha / p;
    const auto f_10 = 2.5 / p;

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

    const auto *ph_s_3 = buffer.data(ph_s + 3);
    const auto *ph_s_8 = buffer.data(ph_s + 8);

    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_8 = buffer.data(ph + 8);

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
    const auto *dh_48 = buffer.data(dh + 48);
    const auto *dh_49 = buffer.data(dh + 49);
    const auto *dh_50 = buffer.data(dh + 50);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_54 = buffer.data(dh + 54);
    const auto *dh_55 = buffer.data(dh + 55);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_30 = buffer.data(ff_s + 30);
    const auto *ff_s_32 = buffer.data(ff_s + 32);
    const auto *ff_s_36 = buffer.data(ff_s + 36);
    const auto *ff_s_44 = buffer.data(ff_s + 44);
    const auto *ff_s_45 = buffer.data(ff_s + 45);
    const auto *ff_s_46 = buffer.data(ff_s + 46);
    const auto *ff_s_47 = buffer.data(ff_s + 47);
    const auto *ff_s_48 = buffer.data(ff_s + 48);
    const auto *ff_s_49 = buffer.data(ff_s + 49);
    const auto *ff_s_50 = buffer.data(ff_s + 50);
    const auto *ff_s_51 = buffer.data(ff_s + 51);

    const auto *fh_s_0 = buffer.data(fh_s + 0);
    const auto *fh_s_1 = buffer.data(fh_s + 1);
    const auto *fh_s_2 = buffer.data(fh_s + 2);
    const auto *fh_s_3 = buffer.data(fh_s + 3);
    const auto *fh_s_4 = buffer.data(fh_s + 4);
    const auto *fh_s_5 = buffer.data(fh_s + 5);
    const auto *fh_s_6 = buffer.data(fh_s + 6);
    const auto *fh_s_7 = buffer.data(fh_s + 7);
    const auto *fh_s_8 = buffer.data(fh_s + 8);
    const auto *fh_s_9 = buffer.data(fh_s + 9);
    const auto *fh_s_10 = buffer.data(fh_s + 10);
    const auto *fh_s_11 = buffer.data(fh_s + 11);
    const auto *fh_s_12 = buffer.data(fh_s + 12);
    const auto *fh_s_13 = buffer.data(fh_s + 13);
    const auto *fh_s_14 = buffer.data(fh_s + 14);
    const auto *fh_s_15 = buffer.data(fh_s + 15);
    const auto *fh_s_16 = buffer.data(fh_s + 16);
    const auto *fh_s_17 = buffer.data(fh_s + 17);
    const auto *fh_s_18 = buffer.data(fh_s + 18);
    const auto *fh_s_19 = buffer.data(fh_s + 19);
    const auto *fh_s_20 = buffer.data(fh_s + 20);
    const auto *fh_s_21 = buffer.data(fh_s + 21);
    const auto *fh_s_22 = buffer.data(fh_s + 22);
    const auto *fh_s_23 = buffer.data(fh_s + 23);
    const auto *fh_s_24 = buffer.data(fh_s + 24);
    const auto *fh_s_25 = buffer.data(fh_s + 25);
    const auto *fh_s_26 = buffer.data(fh_s + 26);
    const auto *fh_s_27 = buffer.data(fh_s + 27);
    const auto *fh_s_28 = buffer.data(fh_s + 28);
    const auto *fh_s_29 = buffer.data(fh_s + 29);
    const auto *fh_s_30 = buffer.data(fh_s + 30);
    const auto *fh_s_31 = buffer.data(fh_s + 31);
    const auto *fh_s_32 = buffer.data(fh_s + 32);
    const auto *fh_s_33 = buffer.data(fh_s + 33);
    const auto *fh_s_34 = buffer.data(fh_s + 34);
    const auto *fh_s_35 = buffer.data(fh_s + 35);
    const auto *fh_s_36 = buffer.data(fh_s + 36);
    const auto *fh_s_37 = buffer.data(fh_s + 37);
    const auto *fh_s_38 = buffer.data(fh_s + 38);
    const auto *fh_s_39 = buffer.data(fh_s + 39);
    const auto *fh_s_40 = buffer.data(fh_s + 40);
    const auto *fh_s_41 = buffer.data(fh_s + 41);
    const auto *fh_s_42 = buffer.data(fh_s + 42);
    const auto *fh_s_43 = buffer.data(fh_s + 43);
    const auto *fh_s_44 = buffer.data(fh_s + 44);
    const auto *fh_s_45 = buffer.data(fh_s + 45);
    const auto *fh_s_46 = buffer.data(fh_s + 46);
    const auto *fh_s_47 = buffer.data(fh_s + 47);
    const auto *fh_s_48 = buffer.data(fh_s + 48);
    const auto *fh_s_49 = buffer.data(fh_s + 49);
    const auto *fh_s_50 = buffer.data(fh_s + 50);
    const auto *fh_s_51 = buffer.data(fh_s + 51);
    const auto *fh_s_52 = buffer.data(fh_s + 52);
    const auto *fh_s_53 = buffer.data(fh_s + 53);
    const auto *fh_s_54 = buffer.data(fh_s + 54);
    const auto *fh_s_55 = buffer.data(fh_s + 55);
    const auto *fh_s_56 = buffer.data(fh_s + 56);
    const auto *fh_s_57 = buffer.data(fh_s + 57);
    const auto *fh_s_58 = buffer.data(fh_s + 58);
    const auto *fh_s_59 = buffer.data(fh_s + 59);
    const auto *fh_s_60 = buffer.data(fh_s + 60);
    const auto *fh_s_61 = buffer.data(fh_s + 61);
    const auto *fh_s_62 = buffer.data(fh_s + 62);
    const auto *fh_s_63 = buffer.data(fh_s + 63);
    const auto *fh_s_64 = buffer.data(fh_s + 64);
    const auto *fh_s_65 = buffer.data(fh_s + 65);
    const auto *fh_s_66 = buffer.data(fh_s + 66);
    const auto *fh_s_67 = buffer.data(fh_s + 67);
    const auto *fh_s_68 = buffer.data(fh_s + 68);
    const auto *fh_s_69 = buffer.data(fh_s + 69);
    const auto *fh_s_70 = buffer.data(fh_s + 70);
    const auto *fh_s_71 = buffer.data(fh_s + 71);
    const auto *fh_s_72 = buffer.data(fh_s + 72);
    const auto *fh_s_73 = buffer.data(fh_s + 73);
    const auto *fh_s_74 = buffer.data(fh_s + 74);
    const auto *fh_s_75 = buffer.data(fh_s + 75);
    const auto *fh_s_76 = buffer.data(fh_s + 76);
    const auto *fh_s_77 = buffer.data(fh_s + 77);
    const auto *fh_s_78 = buffer.data(fh_s + 78);
    const auto *fh_s_79 = buffer.data(fh_s + 79);
    const auto *fh_s_80 = buffer.data(fh_s + 80);
    const auto *fh_s_81 = buffer.data(fh_s + 81);
    const auto *fh_s_82 = buffer.data(fh_s + 82);
    const auto *fh_s_83 = buffer.data(fh_s + 83);
    const auto *fh_s_84 = buffer.data(fh_s + 84);
    const auto *fh_s_85 = buffer.data(fh_s + 85);
    const auto *fh_s_86 = buffer.data(fh_s + 86);
    const auto *fh_s_87 = buffer.data(fh_s + 87);
    const auto *fh_s_88 = buffer.data(fh_s + 88);
    const auto *fh_s_89 = buffer.data(fh_s + 89);
    const auto *fh_s_90 = buffer.data(fh_s + 90);
    const auto *fh_s_91 = buffer.data(fh_s + 91);
    const auto *fh_s_92 = buffer.data(fh_s + 92);
    const auto *fh_s_93 = buffer.data(fh_s + 93);
    const auto *fh_s_94 = buffer.data(fh_s + 94);
    const auto *fh_s_95 = buffer.data(fh_s + 95);
    const auto *fh_s_96 = buffer.data(fh_s + 96);
    const auto *fh_s_97 = buffer.data(fh_s + 97);
    const auto *fh_s_98 = buffer.data(fh_s + 98);
    const auto *fh_s_99 = buffer.data(fh_s + 99);
    const auto *fh_s_100 = buffer.data(fh_s + 100);
    const auto *fh_s_101 = buffer.data(fh_s + 101);
    const auto *fh_s_102 = buffer.data(fh_s + 102);
    const auto *fh_s_103 = buffer.data(fh_s + 103);
    const auto *fh_s_104 = buffer.data(fh_s + 104);
    const auto *fh_s_105 = buffer.data(fh_s + 105);
    const auto *fh_s_106 = buffer.data(fh_s + 106);
    const auto *fh_s_107 = buffer.data(fh_s + 107);
    const auto *fh_s_108 = buffer.data(fh_s + 108);
    const auto *fh_s_109 = buffer.data(fh_s + 109);
    const auto *fh_s_110 = buffer.data(fh_s + 110);
    const auto *fh_s_111 = buffer.data(fh_s + 111);
    const auto *fh_s_112 = buffer.data(fh_s + 112);
    const auto *fh_s_113 = buffer.data(fh_s + 113);
    const auto *fh_s_114 = buffer.data(fh_s + 114);
    const auto *fh_s_115 = buffer.data(fh_s + 115);
    const auto *fh_s_116 = buffer.data(fh_s + 116);
    const auto *fh_s_117 = buffer.data(fh_s + 117);
    const auto *fh_s_118 = buffer.data(fh_s + 118);
    const auto *fh_s_119 = buffer.data(fh_s + 119);
    const auto *fh_s_120 = buffer.data(fh_s + 120);
    const auto *fh_s_121 = buffer.data(fh_s + 121);
    const auto *fh_s_122 = buffer.data(fh_s + 122);
    const auto *fh_s_123 = buffer.data(fh_s + 123);
    const auto *fh_s_124 = buffer.data(fh_s + 124);
    const auto *fh_s_125 = buffer.data(fh_s + 125);
    const auto *fh_s_126 = buffer.data(fh_s + 126);
    const auto *fh_s_127 = buffer.data(fh_s + 127);
    const auto *fh_s_128 = buffer.data(fh_s + 128);
    const auto *fh_s_129 = buffer.data(fh_s + 129);
    const auto *fh_s_130 = buffer.data(fh_s + 130);
    const auto *fh_s_131 = buffer.data(fh_s + 131);
    const auto *fh_s_132 = buffer.data(fh_s + 132);
    const auto *fh_s_133 = buffer.data(fh_s + 133);
    const auto *fh_s_134 = buffer.data(fh_s + 134);
    const auto *fh_s_135 = buffer.data(fh_s + 135);
    const auto *fh_s_136 = buffer.data(fh_s + 136);
    const auto *fh_s_137 = buffer.data(fh_s + 137);
    const auto *fh_s_138 = buffer.data(fh_s + 138);
    const auto *fh_s_139 = buffer.data(fh_s + 139);
    const auto *fh_s_140 = buffer.data(fh_s + 140);
    const auto *fh_s_141 = buffer.data(fh_s + 141);
    const auto *fh_s_142 = buffer.data(fh_s + 142);
    const auto *fh_s_143 = buffer.data(fh_s + 143);
    const auto *fh_s_144 = buffer.data(fh_s + 144);
    const auto *fh_s_145 = buffer.data(fh_s + 145);
    const auto *fh_s_146 = buffer.data(fh_s + 146);
    const auto *fh_s_147 = buffer.data(fh_s + 147);
    const auto *fh_s_148 = buffer.data(fh_s + 148);
    const auto *fh_s_149 = buffer.data(fh_s + 149);
    const auto *fh_s_150 = buffer.data(fh_s + 150);
    const auto *fh_s_151 = buffer.data(fh_s + 151);
    const auto *fh_s_152 = buffer.data(fh_s + 152);
    const auto *fh_s_153 = buffer.data(fh_s + 153);
    const auto *fh_s_154 = buffer.data(fh_s + 154);
    const auto *fh_s_155 = buffer.data(fh_s + 155);
    const auto *fh_s_156 = buffer.data(fh_s + 156);
    const auto *fh_s_157 = buffer.data(fh_s + 157);
    const auto *fh_s_158 = buffer.data(fh_s + 158);
    const auto *fh_s_159 = buffer.data(fh_s + 159);
    const auto *fh_s_160 = buffer.data(fh_s + 160);
    const auto *fh_s_161 = buffer.data(fh_s + 161);
    const auto *fh_s_162 = buffer.data(fh_s + 162);
    const auto *fh_s_163 = buffer.data(fh_s + 163);
    const auto *fh_s_164 = buffer.data(fh_s + 164);
    const auto *fh_s_165 = buffer.data(fh_s + 165);
    const auto *fh_s_166 = buffer.data(fh_s + 166);
    const auto *fh_s_167 = buffer.data(fh_s + 167);
    const auto *fh_s_168 = buffer.data(fh_s + 168);
    const auto *fh_s_169 = buffer.data(fh_s + 169);
    const auto *fh_s_170 = buffer.data(fh_s + 170);
    const auto *fh_s_171 = buffer.data(fh_s + 171);
    const auto *fh_s_172 = buffer.data(fh_s + 172);
    const auto *fh_s_173 = buffer.data(fh_s + 173);
    const auto *fh_s_174 = buffer.data(fh_s + 174);
    const auto *fh_s_175 = buffer.data(fh_s + 175);
    const auto *fh_s_176 = buffer.data(fh_s + 176);
    const auto *fh_s_177 = buffer.data(fh_s + 177);
    const auto *fh_s_178 = buffer.data(fh_s + 178);
    const auto *fh_s_179 = buffer.data(fh_s + 179);
    const auto *fh_s_180 = buffer.data(fh_s + 180);
    const auto *fh_s_181 = buffer.data(fh_s + 181);
    const auto *fh_s_182 = buffer.data(fh_s + 182);
    const auto *fh_s_183 = buffer.data(fh_s + 183);
    const auto *fh_s_184 = buffer.data(fh_s + 184);
    const auto *fh_s_185 = buffer.data(fh_s + 185);
    const auto *fh_s_186 = buffer.data(fh_s + 186);
    const auto *fh_s_187 = buffer.data(fh_s + 187);
    const auto *fh_s_188 = buffer.data(fh_s + 188);
    const auto *fh_s_189 = buffer.data(fh_s + 189);
    const auto *fh_s_190 = buffer.data(fh_s + 190);
    const auto *fh_s_191 = buffer.data(fh_s + 191);
    const auto *fh_s_192 = buffer.data(fh_s + 192);
    const auto *fh_s_193 = buffer.data(fh_s + 193);
    const auto *fh_s_194 = buffer.data(fh_s + 194);
    const auto *fh_s_195 = buffer.data(fh_s + 195);
    const auto *fh_s_196 = buffer.data(fh_s + 196);
    const auto *fh_s_197 = buffer.data(fh_s + 197);
    const auto *fh_s_198 = buffer.data(fh_s + 198);
    const auto *fh_s_199 = buffer.data(fh_s + 199);
    const auto *fh_s_200 = buffer.data(fh_s + 200);
    const auto *fh_s_201 = buffer.data(fh_s + 201);
    const auto *fh_s_202 = buffer.data(fh_s + 202);
    const auto *fh_s_203 = buffer.data(fh_s + 203);
    const auto *fh_s_204 = buffer.data(fh_s + 204);
    const auto *fh_s_205 = buffer.data(fh_s + 205);
    const auto *fh_s_206 = buffer.data(fh_s + 206);
    const auto *fh_s_207 = buffer.data(fh_s + 207);
    const auto *fh_s_208 = buffer.data(fh_s + 208);
    const auto *fh_s_209 = buffer.data(fh_s + 209);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
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
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, dg_0, ff_s_0, fh_s_0, fh_s_1, \
                         fh_s_2, fh_s_3, ff_0, fg_0, fg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 - f_1 * ff_s_0[k]
                 + f_2 * fh_s_0[k]
                 + f_3 * ff_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = f_2 * fh_s_1[k]
                 + pb_y[k] * fg_0[k];

        t_2[k] = f_2 * fh_s_2[k]
                 + pb_z[k] * fg_0[k];

        t_3[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_3[k]
                 + f_5 * ff_0[k]
                 + pb_y[k] * fg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, ff_s_0, ff_s_1, fh_s_4, fh_s_5, \
                         fh_s_6, fh_s_7, ff_0, ff_1, fg_2, fg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * fh_s_4[k]
                 + pb_y[k] * fg_2[k];

        t_5[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_5[k]
                 + f_5 * ff_0[k]
                 + pb_z[k] * fg_2[k];

        t_6[k] = -f_6 * ff_s_1[k]
                 + f_2 * fh_s_6[k]
                 + f_7 * ff_1[k]
                 + pb_y[k] * fg_3[k];

        t_7[k] = f_2 * fh_s_7[k]
                 + pb_z[k] * fg_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pb_x, pb_y, pb_z, dg_5, ff_s_2, fh_s_8, fh_s_9, \
                         fh_s_10, ff_2, fg_4, fg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * fh_s_8[k]
                 + pb_y[k] * fg_4[k];

        t_9[k] = -f_6 * ff_s_2[k]
                 + f_2 * fh_s_9[k]
                 + f_7 * ff_2[k]
                 + pb_z[k] * fg_4[k];

        t_10[k] = f_0 * dg_5[k]
                  + f_2 * fh_s_10[k]
                  + pb_x[k] * fg_7[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pb_x, pb_y, pb_z, dg_6, fh_s_11, fh_s_12, fh_s_13, \
                         fg_5, fg_6, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * fh_s_11[k]
                  + pb_z[k] * fg_5[k];

        t_12[k] = f_0 * dg_6[k]
                  + f_2 * fh_s_12[k]
                  + pb_x[k] * fg_8[k];

        t_13[k] = f_2 * fh_s_13[k]
                  + pb_y[k] * fg_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pb_x, pb_y, pb_z, dg_7, ff_s_3, fh_s_14, fh_s_15, \
                         fh_s_16, ff_3, fg_7, fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * dg_7[k]
                  + f_2 * fh_s_14[k]
                  + pb_x[k] * fg_10[k];

        t_15[k] = -f_1 * ff_s_3[k]
                  + f_2 * fh_s_15[k]
                  + f_3 * ff_3[k]
                  + pb_y[k] * fg_7[k];

        t_16[k] = f_2 * fh_s_16[k]
                  + pb_z[k] * fg_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_y, ff_s_4, ff_s_5, fh_s_17, fh_s_18, fh_s_19, \
                         ff_4, ff_5, fg_8, fg_9, fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_6 * ff_s_4[k]
                  + f_2 * fh_s_17[k]
                  + f_7 * ff_4[k]
                  + pb_y[k] * fg_8[k];

        t_18[k] = -f_4 * ff_s_5[k]
                  + f_2 * fh_s_18[k]
                  + f_5 * ff_5[k]
                  + pb_y[k] * fg_9[k];

        t_19[k] = f_2 * fh_s_19[k]
                  + pb_y[k] * fg_10[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_y, pb_z, dg_0, dh_0, ff_s_5, fh_s_20, \
                         fh_s_21, fh_s_22, ff_5, fg_10, fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_1 * ff_s_5[k]
                  + f_2 * fh_s_20[k]
                  + f_3 * ff_5[k]
                  + pb_z[k] * fg_10[k];

        t_21[k] = pa_y[k] * dh_0[k]
                  + f_2 * fh_s_21[k];

        t_22[k] = f_5 * dg_0[k]
                  + f_2 * fh_s_22[k]
                  + pb_y[k] * fg_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pb_z, dg_1, dh_1, dh_2, fh_s_23, \
                         fh_s_24, fh_s_25, fh_s_26, fg_11, fg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_2 * fh_s_23[k]
                  + pb_z[k] * fg_11[k];

        t_24[k] = f_7 * dg_1[k]
                  + pa_y[k] * dh_1[k]
                  + f_2 * fh_s_24[k];

        t_25[k] = f_2 * fh_s_25[k]
                  + pb_z[k] * fg_12[k];

        t_26[k] = pa_y[k] * dh_2[k]
                  + f_2 * fh_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_y, pb_z, dg_3, dg_4, dh_3, fh_s_27, \
                         fh_s_28, fh_s_29, fg_13, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * dg_3[k]
                  + pa_y[k] * dh_3[k]
                  + f_2 * fh_s_27[k];

        t_28[k] = f_2 * fh_s_28[k]
                  + pb_z[k] * fg_13[k];

        t_29[k] = f_5 * dg_4[k]
                  + f_2 * fh_s_29[k]
                  + pb_y[k] * fg_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pb_x, pb_z, dg_11, dh_5, fh_s_30, fh_s_31, \
                         fh_s_32, fg_15, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * dh_5[k]
                  + f_2 * fh_s_30[k];

        t_31[k] = f_7 * dg_11[k]
                  + f_2 * fh_s_31[k]
                  + pb_x[k] * fg_16[k];

        t_32[k] = f_2 * fh_s_32[k]
                  + pb_z[k] * fg_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pb_x, dg_12, dg_13, dh_7, fh_s_33, fh_s_34, \
                         fh_s_35, fg_18, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_7 * dg_12[k]
                  + f_2 * fh_s_33[k]
                  + pb_x[k] * fg_18[k];

        t_34[k] = f_7 * dg_13[k]
                  + f_2 * fh_s_34[k]
                  + pb_x[k] * fg_19[k];

        t_35[k] = pa_y[k] * dh_7[k]
                  + f_2 * fh_s_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_x, pb_z, ph_s_3, ph_3, dh_14, ff_s_7, fh_s_36, \
                         fh_s_37, fh_s_38, ff_7, fg_16, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_8 * ph_s_3[k]
                  + f_5 * ph_3[k]
                  + pa_x[k] * dh_14[k]
                  + f_2 * fh_s_36[k];

        t_37[k] = f_2 * fh_s_37[k]
                  + pb_z[k] * fg_16[k];

        t_38[k] = -f_4 * ff_s_7[k]
                  + f_2 * fh_s_38[k]
                  + f_5 * ff_7[k]
                  + pb_z[k] * fg_17[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_y, pb_z, dg_7, dh_9, ff_s_8, fh_s_39, \
                         fh_s_40, fh_s_41, ff_8, fg_18, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -f_6 * ff_s_8[k]
                  + f_2 * fh_s_39[k]
                  + f_7 * ff_8[k]
                  + pb_z[k] * fg_18[k];

        t_40[k] = f_5 * dg_7[k]
                  + f_2 * fh_s_40[k]
                  + pb_y[k] * fg_20[k];

        t_41[k] = pa_y[k] * dh_9[k]
                  + f_2 * fh_s_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_z, pb_y, pb_z, dg_0, dh_0, dh_1, fh_s_42, \
                         fh_s_43, fh_s_44, fh_s_45, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * dh_0[k]
                  + f_2 * fh_s_42[k];

        t_43[k] = f_2 * fh_s_43[k]
                  + pb_y[k] * fg_21[k];

        t_44[k] = f_5 * dg_0[k]
                  + f_2 * fh_s_44[k]
                  + pb_z[k] * fg_21[k];

        t_45[k] = pa_z[k] * dh_1[k]
                  + f_2 * fh_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_z, pb_y, dg_2, dg_3, dh_2, dh_3, dh_4, \
                         fh_s_46, fh_s_47, fh_s_48, fh_s_49, fg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_2 * fh_s_46[k]
                  + pb_y[k] * fg_22[k];

        t_47[k] = f_7 * dg_2[k]
                  + pa_z[k] * dh_2[k]
                  + f_2 * fh_s_47[k];

        t_48[k] = pa_z[k] * dh_3[k]
                  + f_2 * fh_s_48[k];

        t_49[k] = f_5 * dg_3[k]
                  + pa_z[k] * dh_4[k]
                  + f_2 * fh_s_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_z, pb_y, dg_4, dh_5, dh_6, fh_s_50, fh_s_51, \
                         fh_s_52, fg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_2 * fh_s_50[k]
                  + pb_y[k] * fg_23[k];

        t_51[k] = f_0 * dg_4[k]
                  + pa_z[k] * dh_5[k]
                  + f_2 * fh_s_51[k];

        t_52[k] = pa_z[k] * dh_6[k]
                  + f_2 * fh_s_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_x, pb_y, dg_17, dg_18, fh_s_53, fh_s_54, \
                         fh_s_55, fg_24, fg_25, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_7 * dg_17[k]
                  + f_2 * fh_s_53[k]
                  + pb_x[k] * fg_25[k];

        t_54[k] = f_7 * dg_18[k]
                  + f_2 * fh_s_54[k]
                  + pb_x[k] * fg_26[k];

        t_55[k] = f_2 * fh_s_55[k]
                  + pb_y[k] * fg_24[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pa_z, pb_x, pb_y, dg_19, dh_8, ff_s_12, fh_s_56, \
                         fh_s_57, fh_s_58, ff_11, fg_25, fg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_7 * dg_19[k]
                  + f_2 * fh_s_56[k]
                  + pb_x[k] * fg_28[k];

        t_57[k] = pa_z[k] * dh_8[k]
                  + f_2 * fh_s_57[k];

        t_58[k] = -f_9 * ff_s_12[k]
                  + f_2 * fh_s_58[k]
                  + f_0 * ff_11[k]
                  + pb_y[k] * fg_25[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pb_y, ff_s_13, ff_s_14, fh_s_59, fh_s_60, fh_s_61, \
                         ff_12, ff_13, fg_26, fg_27, fg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = -f_6 * ff_s_13[k]
                  + f_2 * fh_s_59[k]
                  + f_7 * ff_12[k]
                  + pb_y[k] * fg_26[k];

        t_60[k] = -f_4 * ff_s_14[k]
                  + f_2 * fh_s_60[k]
                  + f_5 * ff_13[k]
                  + pb_y[k] * fg_27[k];

        t_61[k] = f_2 * fh_s_61[k]
                  + pb_y[k] * fg_28[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pa_x, pb_y, ph_s_8, ph_8, dg_8, dg_20, dh_20, \
                         dh_21, fh_s_62, fh_s_63, fh_s_64, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -f_8 * ph_s_8[k]
                  + f_5 * ph_8[k]
                  + pa_x[k] * dh_20[k]
                  + f_2 * fh_s_62[k];

        t_63[k] = f_10 * dg_20[k]
                  + pa_x[k] * dh_21[k]
                  + f_2 * fh_s_63[k];

        t_64[k] = f_7 * dg_8[k]
                  + f_2 * fh_s_64[k]
                  + pb_y[k] * fg_29[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pb_z, dg_22, dg_24, dh_23, dh_25, \
                         fh_s_65, fh_s_66, fh_s_67, fh_s_68, fg_29, \
                         fg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_2 * fh_s_65[k]
                  + pb_z[k] * fg_29[k];

        t_66[k] = f_0 * dg_22[k]
                  + pa_x[k] * dh_23[k]
                  + f_2 * fh_s_66[k];

        t_67[k] = f_2 * fh_s_67[k]
                  + pb_z[k] * fg_30[k];

        t_68[k] = f_0 * dg_24[k]
                  + pa_x[k] * dh_25[k]
                  + f_2 * fh_s_68[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_x, pb_y, pb_z, dg_10, dg_25, dh_26, fh_s_69, \
                         fh_s_70, fh_s_71, fg_31, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_7 * dg_25[k]
                  + pa_x[k] * dh_26[k]
                  + f_2 * fh_s_69[k];

        t_70[k] = f_2 * fh_s_70[k]
                  + pb_z[k] * fg_31[k];

        t_71[k] = f_7 * dg_10[k]
                  + f_2 * fh_s_71[k]
                  + pb_y[k] * fg_32[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_x, pb_x, pb_z, dg_27, dg_28, dh_29, fh_s_72, \
                         fh_s_73, fh_s_74, fg_33, fg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_7 * dg_27[k]
                  + pa_x[k] * dh_29[k]
                  + f_2 * fh_s_72[k];

        t_73[k] = f_5 * dg_28[k]
                  + f_2 * fh_s_73[k]
                  + pb_x[k] * fg_34[k];

        t_74[k] = f_2 * fh_s_74[k]
                  + pb_z[k] * fg_33[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_x, dg_30, dg_31, dg_32, fh_s_75, fh_s_76, \
                         fh_s_77, fg_35, fg_36, fg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_5 * dg_30[k]
                  + f_2 * fh_s_75[k]
                  + pb_x[k] * fg_35[k];

        t_76[k] = f_5 * dg_31[k]
                  + f_2 * fh_s_76[k]
                  + pb_x[k] * fg_36[k];

        t_77[k] = f_5 * dg_32[k]
                  + f_2 * fh_s_77[k]
                  + pb_x[k] * fg_37[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pb_z, dh_30, dh_31, dh_32, fh_s_78, \
                         fh_s_79, fh_s_80, fh_s_81, fg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_x[k] * dh_30[k]
                  + f_2 * fh_s_78[k];

        t_79[k] = f_2 * fh_s_79[k]
                  + pb_z[k] * fg_34[k];

        t_80[k] = pa_x[k] * dh_31[k]
                  + f_2 * fh_s_80[k];

        t_81[k] = pa_x[k] * dh_32[k]
                  + f_2 * fh_s_81[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, pa_y, pa_z, dh_10, dh_15, dh_33, dh_34, \
                         fh_s_82, fh_s_83, fh_s_84, fh_s_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pa_x[k] * dh_33[k]
                  + f_2 * fh_s_82[k];

        t_83[k] = pa_x[k] * dh_34[k]
                  + f_2 * fh_s_83[k];

        t_84[k] = pa_y[k] * dh_15[k]
                  + f_2 * fh_s_84[k];

        t_85[k] = pa_z[k] * dh_10[k]
                  + f_2 * fh_s_85[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_y, pa_z, pb_y, dg_15, dh_11, dh_16, dh_17, \
                         fh_s_86, fh_s_87, fh_s_88, fh_s_89, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pa_y[k] * dh_16[k]
                  + f_2 * fh_s_86[k];

        t_87[k] = pa_z[k] * dh_11[k]
                  + f_2 * fh_s_87[k];

        t_88[k] = f_5 * dg_15[k]
                  + f_2 * fh_s_88[k]
                  + pb_y[k] * fg_38[k];

        t_89[k] = pa_y[k] * dh_17[k]
                  + f_2 * fh_s_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pa_z, pb_y, pb_z, dg_9, dg_16, dh_12, fh_s_90, \
                         fh_s_91, fh_s_92, fg_39, fg_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pa_z[k] * dh_12[k]
                  + f_2 * fh_s_90[k];

        t_91[k] = f_5 * dg_9[k]
                  + f_2 * fh_s_91[k]
                  + pb_z[k] * fg_39[k];

        t_92[k] = f_5 * dg_16[k]
                  + f_2 * fh_s_92[k]
                  + pb_y[k] * fg_40[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pa_y, pa_z, pb_x, dg_34, dh_13, dh_18, fh_s_93, \
                         fh_s_94, fh_s_95, fg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_y[k] * dh_18[k]
                  + f_2 * fh_s_93[k];

        t_94[k] = pa_z[k] * dh_13[k]
                  + f_2 * fh_s_94[k];

        t_95[k] = f_5 * dg_34[k]
                  + f_2 * fh_s_95[k]
                  + pb_x[k] * fg_41[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pa_y, pb_x, dg_35, dg_36, dh_19, fh_s_96, fh_s_97, \
                         fh_s_98, fg_42, fg_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_5 * dg_35[k]
                  + f_2 * fh_s_96[k]
                  + pb_x[k] * fg_42[k];

        t_97[k] = f_5 * dg_36[k]
                  + f_2 * fh_s_97[k]
                  + pb_x[k] * fg_43[k];

        t_98[k] = pa_y[k] * dh_19[k]
                  + f_2 * fh_s_98[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pa_x, dh_35, dh_36, dh_37, dh_38, \
                         dh_39, fh_s_99, fh_s_100, fh_s_101, fh_s_102, \
                         fh_s_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_x[k] * dh_35[k]
                  + f_2 * fh_s_99[k];

        t_100[k] = pa_x[k] * dh_36[k]
                   + f_2 * fh_s_100[k];

        t_101[k] = pa_x[k] * dh_37[k]
                   + f_2 * fh_s_101[k];

        t_102[k] = pa_x[k] * dh_38[k]
                   + f_2 * fh_s_102[k];

        t_103[k] = pa_x[k] * dh_39[k]
                   + f_2 * fh_s_103[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_x, pb_y, pb_z, dg_14, dg_38, dh_40, \
                         dh_41, fh_s_104, fh_s_105, fh_s_106, fh_s_107, \
                         fg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_x[k] * dh_40[k]
                   + f_2 * fh_s_104[k];

        t_105[k] = f_10 * dg_38[k]
                   + pa_x[k] * dh_41[k]
                   + f_2 * fh_s_105[k];

        t_106[k] = f_2 * fh_s_106[k]
                   + pb_y[k] * fg_44[k];

        t_107[k] = f_7 * dg_14[k]
                   + f_2 * fh_s_107[k]
                   + pb_z[k] * fg_44[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_x, pb_y, dg_41, dg_43, dh_44, dh_46, \
                         fh_s_108, fh_s_109, fh_s_110, fg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_0 * dg_41[k]
                   + pa_x[k] * dh_44[k]
                   + f_2 * fh_s_108[k];

        t_109[k] = f_2 * fh_s_109[k]
                   + pb_y[k] * fg_45[k];

        t_110[k] = f_0 * dg_43[k]
                   + pa_x[k] * dh_46[k]
                   + f_2 * fh_s_110[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pa_x, pb_y, dg_44, dg_45, dh_47, dh_48, \
                         fh_s_111, fh_s_112, fh_s_113, fg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_7 * dg_44[k]
                   + pa_x[k] * dh_47[k]
                   + f_2 * fh_s_111[k];

        t_112[k] = f_7 * dg_45[k]
                   + pa_x[k] * dh_48[k]
                   + f_2 * fh_s_112[k];

        t_113[k] = f_2 * fh_s_113[k]
                   + pb_y[k] * fg_46[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pa_x, pb_x, dg_46, dg_47, dg_48, dh_50, \
                         fh_s_114, fh_s_115, fh_s_116, fg_48, fg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_7 * dg_46[k]
                   + pa_x[k] * dh_50[k]
                   + f_2 * fh_s_114[k];

        t_115[k] = f_5 * dg_47[k]
                   + f_2 * fh_s_115[k]
                   + pb_x[k] * fg_48[k];

        t_116[k] = f_5 * dg_48[k]
                   + f_2 * fh_s_116[k]
                   + pb_x[k] * fg_49[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pb_x, pb_y, dg_49, dg_51, fh_s_117, fh_s_118, \
                         fh_s_119, fg_47, fg_50, fg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_5 * dg_49[k]
                   + f_2 * fh_s_117[k]
                   + pb_x[k] * fg_50[k];

        t_118[k] = f_2 * fh_s_118[k]
                   + pb_y[k] * fg_47[k];

        t_119[k] = f_5 * dg_51[k]
                   + f_2 * fh_s_119[k]
                   + pb_x[k] * fg_51[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_x, dh_51, dh_52, dh_53, dh_54, \
                         fh_s_120, fh_s_121, fh_s_122, fh_s_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_x[k] * dh_51[k]
                   + f_2 * fh_s_120[k];

        t_121[k] = pa_x[k] * dh_52[k]
                   + f_2 * fh_s_121[k];

        t_122[k] = pa_x[k] * dh_53[k]
                   + f_2 * fh_s_122[k];

        t_123[k] = pa_x[k] * dh_54[k]
                   + f_2 * fh_s_123[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_x, pb_x, pb_y, dh_55, ff_s_22, fh_s_124, \
                         fh_s_125, fh_s_126, ff_18, fg_51, fg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_2 * fh_s_124[k]
                   + pb_y[k] * fg_51[k];

        t_125[k] = pa_x[k] * dh_55[k]
                   + f_2 * fh_s_125[k];

        t_126[k] = -f_1 * ff_s_22[k]
                   + f_2 * fh_s_126[k]
                   + f_3 * ff_18[k]
                   + pb_x[k] * fg_52[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pb_x, pb_z, ff_s_23, ff_s_24, fh_s_127, \
                         fh_s_128, fh_s_129, ff_19, ff_20, fg_52, fg_53, \
                         fg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -f_9 * ff_s_23[k]
                   + f_2 * fh_s_127[k]
                   + f_0 * ff_19[k]
                   + pb_x[k] * fg_53[k];

        t_128[k] = f_2 * fh_s_128[k]
                   + pb_z[k] * fg_52[k];

        t_129[k] = -f_6 * ff_s_24[k]
                   + f_2 * fh_s_129[k]
                   + f_7 * ff_20[k]
                   + pb_x[k] * fg_54[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pb_x, pb_z, ff_s_25, ff_s_26, fh_s_130, \
                         fh_s_131, fh_s_132, ff_21, ff_22, fg_53, fg_55, \
                         fg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_2 * fh_s_130[k]
                   + pb_z[k] * fg_53[k];

        t_131[k] = -f_6 * ff_s_25[k]
                   + f_2 * fh_s_131[k]
                   + f_7 * ff_21[k]
                   + pb_x[k] * fg_55[k];

        t_132[k] = -f_4 * ff_s_26[k]
                   + f_2 * fh_s_132[k]
                   + f_5 * ff_22[k]
                   + pb_x[k] * fg_56[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pb_x, pb_z, ff_s_28, ff_s_29, fh_s_133, \
                         fh_s_134, fh_s_135, ff_24, ff_25, fg_54, fg_57, \
                         fg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_2 * fh_s_133[k]
                   + pb_z[k] * fg_54[k];

        t_134[k] = -f_4 * ff_s_28[k]
                   + f_2 * fh_s_134[k]
                   + f_5 * ff_24[k]
                   + pb_x[k] * fg_57[k];

        t_135[k] = -f_4 * ff_s_29[k]
                   + f_2 * fh_s_135[k]
                   + f_5 * ff_25[k]
                   + pb_x[k] * fg_58[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, pb_x, fh_s_136, fh_s_137, \
                         fh_s_138, fh_s_139, fh_s_140, fg_59, fg_60, fg_61, fg_62, \
                         fg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_2 * fh_s_136[k]
                   + pb_x[k] * fg_59[k];

        t_137[k] = f_2 * fh_s_137[k]
                   + pb_x[k] * fg_60[k];

        t_138[k] = f_2 * fh_s_138[k]
                   + pb_x[k] * fg_61[k];

        t_139[k] = f_2 * fh_s_139[k]
                   + pb_x[k] * fg_62[k];

        t_140[k] = f_2 * fh_s_140[k]
                   + pb_x[k] * fg_63[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pb_y, pb_z, dg_28, ff_s_26, fh_s_141, fh_s_142, \
                         fh_s_143, ff_22, fg_59, fg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_0 * dg_28[k]
                   - f_1 * ff_s_26[k]
                   + f_2 * fh_s_141[k]
                   + f_3 * ff_22[k]
                   + pb_y[k] * fg_59[k];

        t_142[k] = f_2 * fh_s_142[k]
                   + pb_z[k] * fg_59[k];

        t_143[k] = -f_4 * ff_s_26[k]
                   + f_2 * fh_s_143[k]
                   + f_5 * ff_22[k]
                   + pb_z[k] * fg_60[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pb_y, pb_z, dg_32, ff_s_27, ff_s_29, fh_s_144, \
                         fh_s_145, fh_s_146, ff_23, ff_25, fg_61, \
                         fg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -f_6 * ff_s_27[k]
                   + f_2 * fh_s_144[k]
                   + f_7 * ff_23[k]
                   + pb_z[k] * fg_61[k];

        t_145[k] = f_0 * dg_32[k]
                   + f_2 * fh_s_145[k]
                   + pb_y[k] * fg_63[k];

        t_146[k] = -f_1 * ff_s_29[k]
                   + f_2 * fh_s_146[k]
                   + f_3 * ff_25[k]
                   + pb_z[k] * fg_63[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pa_z, pb_x, dh_21, dh_22, dh_23, ff_s_30, \
                         fh_s_147, fh_s_148, fh_s_149, fh_s_150, ff_26, \
                         fg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = pa_z[k] * dh_21[k]
                   + f_2 * fh_s_147[k];

        t_148[k] = pa_z[k] * dh_22[k]
                   + f_2 * fh_s_148[k];

        t_149[k] = -f_9 * ff_s_30[k]
                   + f_2 * fh_s_149[k]
                   + f_0 * ff_26[k]
                   + pb_x[k] * fg_64[k];

        t_150[k] = pa_z[k] * dh_23[k]
                   + f_2 * fh_s_150[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pa_z, pb_x, dg_21, dh_24, dh_26, ff_s_32, \
                         fh_s_151, fh_s_152, fh_s_153, ff_27, fg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_5 * dg_21[k]
                   + pa_z[k] * dh_24[k]
                   + f_2 * fh_s_151[k];

        t_152[k] = -f_6 * ff_s_32[k]
                   + f_2 * fh_s_152[k]
                   + f_7 * ff_27[k]
                   + pb_x[k] * fg_65[k];

        t_153[k] = pa_z[k] * dh_26[k]
                   + f_2 * fh_s_153[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pa_z, pb_x, dg_22, dg_23, dh_27, dh_28, ff_s_36, \
                         fh_s_154, fh_s_155, fh_s_156, ff_29, fg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_5 * dg_22[k]
                   + pa_z[k] * dh_27[k]
                   + f_2 * fh_s_154[k];

        t_155[k] = f_7 * dg_23[k]
                   + pa_z[k] * dh_28[k]
                   + f_2 * fh_s_155[k];

        t_156[k] = -f_4 * ff_s_36[k]
                   + f_2 * fh_s_156[k]
                   + f_5 * ff_29[k]
                   + pb_x[k] * fg_66[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, pb_x, fh_s_157, fh_s_158, \
                         fh_s_159, fh_s_160, fh_s_161, fg_67, fg_68, fg_69, fg_70, \
                         fg_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_2 * fh_s_157[k]
                   + pb_x[k] * fg_67[k];

        t_158[k] = f_2 * fh_s_158[k]
                   + pb_x[k] * fg_68[k];

        t_159[k] = f_2 * fh_s_159[k]
                   + pb_x[k] * fg_69[k];

        t_160[k] = f_2 * fh_s_160[k]
                   + pb_x[k] * fg_70[k];

        t_161[k] = f_2 * fh_s_161[k]
                   + pb_x[k] * fg_71[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pa_z, pb_z, dg_28, dg_29, dh_30, dh_31, \
                         fh_s_162, fh_s_163, fh_s_164, fg_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_z[k] * dh_30[k]
                   + f_2 * fh_s_162[k];

        t_163[k] = f_5 * dg_28[k]
                   + f_2 * fh_s_163[k]
                   + pb_z[k] * fg_67[k];

        t_164[k] = f_7 * dg_29[k]
                   + pa_z[k] * dh_31[k]
                   + f_2 * fh_s_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pa_y, pa_z, pb_y, ph_s_8, ph_8, dg_30, dg_37, \
                         dh_32, dh_40, fh_s_165, fh_s_166, fh_s_167, \
                         fg_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_0 * dg_30[k]
                   + pa_z[k] * dh_32[k]
                   + f_2 * fh_s_165[k];

        t_166[k] = f_7 * dg_37[k]
                   + f_2 * fh_s_166[k]
                   + pb_y[k] * fg_71[k];

        t_167[k] = -f_8 * ph_s_8[k]
                   + f_5 * ph_8[k]
                   + pa_y[k] * dh_40[k]
                   + f_2 * fh_s_167[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_y, dg_38, dg_39, dh_41, dh_42, dh_43, \
                         dh_44, fh_s_168, fh_s_169, fh_s_170, \
                         fh_s_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_y[k] * dh_41[k]
                   + f_2 * fh_s_168[k];

        t_169[k] = f_5 * dg_38[k]
                   + pa_y[k] * dh_42[k]
                   + f_2 * fh_s_169[k];

        t_170[k] = pa_y[k] * dh_43[k]
                   + f_2 * fh_s_170[k];

        t_171[k] = f_7 * dg_39[k]
                   + pa_y[k] * dh_44[k]
                   + f_2 * fh_s_171[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, dg_40, dg_41, dg_42, dh_45, dh_46, \
                         dh_47, dh_48, fh_s_172, fh_s_173, fh_s_174, \
                         fh_s_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_5 * dg_40[k]
                   + pa_y[k] * dh_45[k]
                   + f_2 * fh_s_172[k];

        t_173[k] = pa_y[k] * dh_46[k]
                   + f_2 * fh_s_173[k];

        t_174[k] = f_0 * dg_41[k]
                   + pa_y[k] * dh_47[k]
                   + f_2 * fh_s_174[k];

        t_175[k] = f_7 * dg_42[k]
                   + pa_y[k] * dh_48[k]
                   + f_2 * fh_s_175[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pb_x, dg_43, dh_49, dh_50, \
                         fh_s_176, fh_s_177, fh_s_178, fh_s_179, fg_72, \
                         fg_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_5 * dg_43[k]
                   + pa_y[k] * dh_49[k]
                   + f_2 * fh_s_176[k];

        t_177[k] = pa_y[k] * dh_50[k]
                   + f_2 * fh_s_177[k];

        t_178[k] = f_2 * fh_s_178[k]
                   + pb_x[k] * fg_72[k];

        t_179[k] = f_2 * fh_s_179[k]
                   + pb_x[k] * fg_73[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pb_x, dg_47, dh_51, fh_s_180, \
                         fh_s_181, fh_s_182, fh_s_183, fg_74, fg_75, \
                         fg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_2 * fh_s_180[k]
                   + pb_x[k] * fg_74[k];

        t_181[k] = f_2 * fh_s_181[k]
                   + pb_x[k] * fg_75[k];

        t_182[k] = f_2 * fh_s_182[k]
                   + pb_x[k] * fg_76[k];

        t_183[k] = f_10 * dg_47[k]
                   + pa_y[k] * dh_51[k]
                   + f_2 * fh_s_183[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pa_y, pb_z, dg_33, dg_49, dg_50, dh_53, dh_54, \
                         fh_s_184, fh_s_185, fh_s_186, fg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_7 * dg_33[k]
                   + f_2 * fh_s_184[k]
                   + pb_z[k] * fg_72[k];

        t_185[k] = f_0 * dg_49[k]
                   + pa_y[k] * dh_53[k]
                   + f_2 * fh_s_185[k];

        t_186[k] = f_7 * dg_50[k]
                   + pa_y[k] * dh_54[k]
                   + f_2 * fh_s_186[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pa_y, pb_x, pb_y, dg_51, dh_55, ff_s_44, \
                         fh_s_187, fh_s_188, fh_s_189, ff_32, fg_76, \
                         fg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_5 * dg_51[k]
                   + f_2 * fh_s_187[k]
                   + pb_y[k] * fg_76[k];

        t_188[k] = pa_y[k] * dh_55[k]
                   + f_2 * fh_s_188[k];

        t_189[k] = -f_1 * ff_s_44[k]
                   + f_2 * fh_s_189[k]
                   + f_3 * ff_32[k]
                   + pb_x[k] * fg_77[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pb_x, pb_y, ff_s_45, ff_s_46, fh_s_190, \
                         fh_s_191, fh_s_192, ff_33, ff_34, fg_77, fg_78, \
                         fg_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_2 * fh_s_190[k]
                   + pb_y[k] * fg_77[k];

        t_191[k] = -f_9 * ff_s_45[k]
                   + f_2 * fh_s_191[k]
                   + f_0 * ff_33[k]
                   + pb_x[k] * fg_78[k];

        t_192[k] = -f_6 * ff_s_46[k]
                   + f_2 * fh_s_192[k]
                   + f_7 * ff_34[k]
                   + pb_x[k] * fg_79[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pb_x, pb_y, ff_s_47, ff_s_48, fh_s_193, \
                         fh_s_194, fh_s_195, ff_35, ff_36, fg_78, fg_80, \
                         fg_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_2 * fh_s_193[k]
                   + pb_y[k] * fg_78[k];

        t_194[k] = -f_6 * ff_s_47[k]
                   + f_2 * fh_s_194[k]
                   + f_7 * ff_35[k]
                   + pb_x[k] * fg_80[k];

        t_195[k] = -f_4 * ff_s_48[k]
                   + f_2 * fh_s_195[k]
                   + f_5 * ff_36[k]
                   + pb_x[k] * fg_81[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, pb_x, pb_y, ff_s_49, ff_s_51, fh_s_196, \
                         fh_s_197, fh_s_198, ff_37, ff_39, fg_80, fg_82, \
                         fg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = -f_4 * ff_s_49[k]
                   + f_2 * fh_s_196[k]
                   + f_5 * ff_37[k]
                   + pb_x[k] * fg_82[k];

        t_197[k] = f_2 * fh_s_197[k]
                   + pb_y[k] * fg_80[k];

        t_198[k] = -f_4 * ff_s_51[k]
                   + f_2 * fh_s_198[k]
                   + f_5 * ff_39[k]
                   + pb_x[k] * fg_83[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, pb_x, fh_s_199, fh_s_200, \
                         fh_s_201, fh_s_202, fh_s_203, fg_84, fg_85, fg_86, fg_87, \
                         fg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_2 * fh_s_199[k]
                   + pb_x[k] * fg_84[k];

        t_200[k] = f_2 * fh_s_200[k]
                   + pb_x[k] * fg_85[k];

        t_201[k] = f_2 * fh_s_201[k]
                   + pb_x[k] * fg_86[k];

        t_202[k] = f_2 * fh_s_202[k]
                   + pb_x[k] * fg_87[k];

        t_203[k] = f_2 * fh_s_203[k]
                   + pb_x[k] * fg_88[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pb_y, ff_s_48, ff_s_49, ff_s_50, fh_s_204, \
                         fh_s_205, fh_s_206, ff_36, ff_37, ff_38, fg_84, fg_85, \
                         fg_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = -f_1 * ff_s_48[k]
                   + f_2 * fh_s_204[k]
                   + f_3 * ff_36[k]
                   + pb_y[k] * fg_84[k];

        t_205[k] = -f_9 * ff_s_49[k]
                   + f_2 * fh_s_205[k]
                   + f_0 * ff_37[k]
                   + pb_y[k] * fg_85[k];

        t_206[k] = -f_6 * ff_s_50[k]
                   + f_2 * fh_s_206[k]
                   + f_7 * ff_38[k]
                   + pb_y[k] * fg_86[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pb_y, pb_z, dg_51, ff_s_51, fh_s_207, fh_s_208, \
                         fh_s_209, ff_39, fg_87, fg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -f_4 * ff_s_51[k]
                   + f_2 * fh_s_207[k]
                   + f_5 * ff_39[k]
                   + pb_y[k] * fg_87[k];

        t_208[k] = f_2 * fh_s_208[k]
                   + pb_y[k] * fg_88[k];

        t_209[k] = f_0 * dg_51[k]
                   - f_1 * ff_s_51[k]
                   + f_2 * fh_s_209[k]
                   + f_3 * ff_39[k]
                   + pb_z[k] * fg_88[k];
    }
}

auto
compute_prim_fh_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t ph_s, const size_t ph,
                                 const size_t dg, const size_t dh, const size_t ff_s,
                                 const size_t fh_s, const size_t ff, const size_t fg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = beta / p;
    const auto f_9 = 3.0 * alpha / p;
    const auto f_10 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_s_8 = buffer.data(ph_s + 8);
    const auto *ph_s_22 = buffer.data(ph_s + 22);

    const auto *ph_8 = buffer.data(ph + 8);
    const auto *ph_21 = buffer.data(ph + 21);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_43 = buffer.data(dh + 43);
    const auto *dh_45 = buffer.data(dh + 45);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_48 = buffer.data(dh + 48);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_54 = buffer.data(dh + 54);
    const auto *dh_55 = buffer.data(dh + 55);
    const auto *dh_56 = buffer.data(dh + 56);
    const auto *dh_61 = buffer.data(dh + 61);
    const auto *dh_65 = buffer.data(dh + 65);
    const auto *dh_70 = buffer.data(dh + 70);
    const auto *dh_71 = buffer.data(dh + 71);
    const auto *dh_72 = buffer.data(dh + 72);
    const auto *dh_73 = buffer.data(dh + 73);
    const auto *dh_75 = buffer.data(dh + 75);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_32 = buffer.data(ff_s + 32);
    const auto *ff_s_40 = buffer.data(ff_s + 40);
    const auto *ff_s_41 = buffer.data(ff_s + 41);
    const auto *ff_s_42 = buffer.data(ff_s + 42);
    const auto *ff_s_43 = buffer.data(ff_s + 43);
    const auto *ff_s_44 = buffer.data(ff_s + 44);
    const auto *ff_s_45 = buffer.data(ff_s + 45);
    const auto *ff_s_46 = buffer.data(ff_s + 46);
    const auto *ff_s_47 = buffer.data(ff_s + 47);

    const auto *fh_s_0 = buffer.data(fh_s + 0);
    const auto *fh_s_1 = buffer.data(fh_s + 1);
    const auto *fh_s_2 = buffer.data(fh_s + 2);
    const auto *fh_s_3 = buffer.data(fh_s + 3);
    const auto *fh_s_4 = buffer.data(fh_s + 4);
    const auto *fh_s_5 = buffer.data(fh_s + 5);
    const auto *fh_s_6 = buffer.data(fh_s + 6);
    const auto *fh_s_7 = buffer.data(fh_s + 7);
    const auto *fh_s_8 = buffer.data(fh_s + 8);
    const auto *fh_s_9 = buffer.data(fh_s + 9);
    const auto *fh_s_10 = buffer.data(fh_s + 10);
    const auto *fh_s_11 = buffer.data(fh_s + 11);
    const auto *fh_s_12 = buffer.data(fh_s + 12);
    const auto *fh_s_13 = buffer.data(fh_s + 13);
    const auto *fh_s_14 = buffer.data(fh_s + 14);
    const auto *fh_s_15 = buffer.data(fh_s + 15);
    const auto *fh_s_16 = buffer.data(fh_s + 16);
    const auto *fh_s_17 = buffer.data(fh_s + 17);
    const auto *fh_s_18 = buffer.data(fh_s + 18);
    const auto *fh_s_19 = buffer.data(fh_s + 19);
    const auto *fh_s_21 = buffer.data(fh_s + 21);
    const auto *fh_s_22 = buffer.data(fh_s + 22);
    const auto *fh_s_23 = buffer.data(fh_s + 23);
    const auto *fh_s_24 = buffer.data(fh_s + 24);
    const auto *fh_s_25 = buffer.data(fh_s + 25);
    const auto *fh_s_26 = buffer.data(fh_s + 26);
    const auto *fh_s_27 = buffer.data(fh_s + 27);
    const auto *fh_s_28 = buffer.data(fh_s + 28);
    const auto *fh_s_29 = buffer.data(fh_s + 29);
    const auto *fh_s_30 = buffer.data(fh_s + 30);
    const auto *fh_s_32 = buffer.data(fh_s + 32);
    const auto *fh_s_35 = buffer.data(fh_s + 35);
    const auto *fh_s_36 = buffer.data(fh_s + 36);
    const auto *fh_s_37 = buffer.data(fh_s + 37);
    const auto *fh_s_38 = buffer.data(fh_s + 38);
    const auto *fh_s_39 = buffer.data(fh_s + 39);
    const auto *fh_s_40 = buffer.data(fh_s + 40);
    const auto *fh_s_41 = buffer.data(fh_s + 41);
    const auto *fh_s_42 = buffer.data(fh_s + 42);
    const auto *fh_s_43 = buffer.data(fh_s + 43);
    const auto *fh_s_45 = buffer.data(fh_s + 45);
    const auto *fh_s_46 = buffer.data(fh_s + 46);
    const auto *fh_s_47 = buffer.data(fh_s + 47);
    const auto *fh_s_49 = buffer.data(fh_s + 49);
    const auto *fh_s_50 = buffer.data(fh_s + 50);
    const auto *fh_s_51 = buffer.data(fh_s + 51);
    const auto *fh_s_52 = buffer.data(fh_s + 52);
    const auto *fh_s_53 = buffer.data(fh_s + 53);
    const auto *fh_s_54 = buffer.data(fh_s + 54);
    const auto *fh_s_55 = buffer.data(fh_s + 55);
    const auto *fh_s_56 = buffer.data(fh_s + 56);
    const auto *fh_s_57 = buffer.data(fh_s + 57);
    const auto *fh_s_58 = buffer.data(fh_s + 58);
    const auto *fh_s_59 = buffer.data(fh_s + 59);
    const auto *fh_s_60 = buffer.data(fh_s + 60);
    const auto *fh_s_61 = buffer.data(fh_s + 61);
    const auto *fh_s_62 = buffer.data(fh_s + 62);
    const auto *fh_s_63 = buffer.data(fh_s + 63);
    const auto *fh_s_64 = buffer.data(fh_s + 64);
    const auto *fh_s_65 = buffer.data(fh_s + 65);
    const auto *fh_s_67 = buffer.data(fh_s + 67);
    const auto *fh_s_70 = buffer.data(fh_s + 70);
    const auto *fh_s_74 = buffer.data(fh_s + 74);
    const auto *fh_s_75 = buffer.data(fh_s + 75);
    const auto *fh_s_76 = buffer.data(fh_s + 76);
    const auto *fh_s_77 = buffer.data(fh_s + 77);
    const auto *fh_s_78 = buffer.data(fh_s + 78);
    const auto *fh_s_79 = buffer.data(fh_s + 79);
    const auto *fh_s_80 = buffer.data(fh_s + 80);
    const auto *fh_s_81 = buffer.data(fh_s + 81);
    const auto *fh_s_82 = buffer.data(fh_s + 82);
    const auto *fh_s_83 = buffer.data(fh_s + 83);
    const auto *fh_s_84 = buffer.data(fh_s + 84);
    const auto *fh_s_85 = buffer.data(fh_s + 85);
    const auto *fh_s_86 = buffer.data(fh_s + 86);
    const auto *fh_s_87 = buffer.data(fh_s + 87);
    const auto *fh_s_88 = buffer.data(fh_s + 88);
    const auto *fh_s_89 = buffer.data(fh_s + 89);
    const auto *fh_s_90 = buffer.data(fh_s + 90);
    const auto *fh_s_91 = buffer.data(fh_s + 91);
    const auto *fh_s_92 = buffer.data(fh_s + 92);
    const auto *fh_s_93 = buffer.data(fh_s + 93);
    const auto *fh_s_94 = buffer.data(fh_s + 94);
    const auto *fh_s_95 = buffer.data(fh_s + 95);
    const auto *fh_s_96 = buffer.data(fh_s + 96);
    const auto *fh_s_97 = buffer.data(fh_s + 97);
    const auto *fh_s_98 = buffer.data(fh_s + 98);
    const auto *fh_s_100 = buffer.data(fh_s + 100);
    const auto *fh_s_103 = buffer.data(fh_s + 103);
    const auto *fh_s_107 = buffer.data(fh_s + 107);
    const auto *fh_s_108 = buffer.data(fh_s + 108);
    const auto *fh_s_109 = buffer.data(fh_s + 109);
    const auto *fh_s_110 = buffer.data(fh_s + 110);
    const auto *fh_s_111 = buffer.data(fh_s + 111);
    const auto *fh_s_112 = buffer.data(fh_s + 112);
    const auto *fh_s_113 = buffer.data(fh_s + 113);
    const auto *fh_s_124 = buffer.data(fh_s + 124);
    const auto *fh_s_125 = buffer.data(fh_s + 125);
    const auto *fh_s_126 = buffer.data(fh_s + 126);
    const auto *fh_s_127 = buffer.data(fh_s + 127);
    const auto *fh_s_128 = buffer.data(fh_s + 128);
    const auto *fh_s_129 = buffer.data(fh_s + 129);
    const auto *fh_s_130 = buffer.data(fh_s + 130);
    const auto *fh_s_131 = buffer.data(fh_s + 131);
    const auto *fh_s_132 = buffer.data(fh_s + 132);
    const auto *fh_s_133 = buffer.data(fh_s + 133);
    const auto *fh_s_134 = buffer.data(fh_s + 134);
    const auto *fh_s_135 = buffer.data(fh_s + 135);
    const auto *fh_s_136 = buffer.data(fh_s + 136);
    const auto *fh_s_137 = buffer.data(fh_s + 137);
    const auto *fh_s_138 = buffer.data(fh_s + 138);
    const auto *fh_s_139 = buffer.data(fh_s + 139);
    const auto *fh_s_140 = buffer.data(fh_s + 140);
    const auto *fh_s_141 = buffer.data(fh_s + 141);
    const auto *fh_s_142 = buffer.data(fh_s + 142);
    const auto *fh_s_143 = buffer.data(fh_s + 143);
    const auto *fh_s_144 = buffer.data(fh_s + 144);
    const auto *fh_s_145 = buffer.data(fh_s + 145);
    const auto *fh_s_146 = buffer.data(fh_s + 146);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, dg_0, ff_s_0, fh_s_0, fh_s_1, \
                         fh_s_2, fh_s_3, ff_0, fg_0, fg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 - f_1 * ff_s_0[k]
                 + f_2 * fh_s_0[k]
                 + f_3 * ff_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = f_2 * fh_s_1[k]
                 + pb_y[k] * fg_0[k];

        t_2[k] = f_2 * fh_s_2[k]
                 + pb_z[k] * fg_0[k];

        t_3[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_3[k]
                 + f_5 * ff_0[k]
                 + pb_y[k] * fg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, ff_s_0, ff_s_1, fh_s_4, fh_s_5, fh_s_6, \
                         ff_0, ff_1, fg_2, fg_3, fg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_4[k]
                 + f_5 * ff_0[k]
                 + pb_z[k] * fg_2[k];

        t_5[k] = -f_6 * ff_s_1[k]
                 + f_2 * fh_s_5[k]
                 + f_7 * ff_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = f_2 * fh_s_6[k]
                 + pb_y[k] * fg_4[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_x, pb_z, dg_5, dg_8, ff_s_2, fh_s_7, fh_s_8, \
                         fh_s_9, ff_2, fg_4, fg_5, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -f_6 * ff_s_2[k]
                 + f_2 * fh_s_7[k]
                 + f_7 * ff_2[k]
                 + pb_z[k] * fg_4[k];

        t_8[k] = f_0 * dg_5[k]
                 + f_2 * fh_s_8[k]
                 + pb_x[k] * fg_5[k];

        t_9[k] = f_0 * dg_8[k]
                 + f_2 * fh_s_9[k]
                 + pb_x[k] * fg_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, ff_s_3, ff_s_4, ff_s_5, fh_s_10, fh_s_11, \
                         fh_s_12, ff_3, ff_4, ff_5, fg_5, fg_6, fg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * ff_s_3[k]
                  + f_2 * fh_s_10[k]
                  + f_3 * ff_3[k]
                  + pb_y[k] * fg_5[k];

        t_11[k] = -f_6 * ff_s_4[k]
                  + f_2 * fh_s_11[k]
                  + f_7 * ff_4[k]
                  + pb_y[k] * fg_6[k];

        t_12[k] = -f_4 * ff_s_5[k]
                  + f_2 * fh_s_12[k]
                  + f_5 * ff_5[k]
                  + pb_y[k] * fg_7[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_y, pb_z, dh_0, ff_s_5, fh_s_13, fh_s_14, \
                         fh_s_15, ff_5, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * fh_s_13[k]
                  + pb_y[k] * fg_8[k];

        t_14[k] = -f_1 * ff_s_5[k]
                  + f_2 * fh_s_14[k]
                  + f_3 * ff_5[k]
                  + pb_z[k] * fg_8[k];

        t_15[k] = pa_y[k] * dh_0[k]
                  + f_2 * fh_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_y, pb_y, dg_0, dg_1, dh_3, dh_4, fh_s_16, \
                         fh_s_17, fh_s_18, fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * dg_0[k]
                  + f_2 * fh_s_16[k]
                  + pb_y[k] * fg_9[k];

        t_17[k] = f_7 * dg_1[k]
                  + pa_y[k] * dh_3[k]
                  + f_2 * fh_s_17[k];

        t_18[k] = pa_y[k] * dh_4[k]
                  + f_2 * fh_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pb_x, dg_3, dg_10, dh_5, dh_8, fh_s_19, \
                         fh_s_21, fh_s_22, fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * dg_3[k]
                  + pa_y[k] * dh_5[k]
                  + f_2 * fh_s_19[k];

        t_20[k] = pa_y[k] * dh_8[k]
                  + f_2 * fh_s_21[k];

        t_21[k] = f_7 * dg_10[k]
                  + f_2 * fh_s_22[k]
                  + pb_x[k] * fg_10[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pb_z, ph_s_8, ph_8, dh_18, ff_s_6, fh_s_23, \
                         fh_s_24, fh_s_25, ff_6, fg_10, fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -f_8 * ph_s_8[k]
                  + f_5 * ph_8[k]
                  + pa_x[k] * dh_18[k]
                  + f_2 * fh_s_23[k];

        t_23[k] = f_2 * fh_s_24[k]
                  + pb_z[k] * fg_10[k];

        t_24[k] = -f_4 * ff_s_6[k]
                  + f_2 * fh_s_25[k]
                  + f_5 * ff_6[k]
                  + pb_z[k] * fg_11[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pb_y, pb_z, dg_8, dh_12, ff_s_7, fh_s_26, \
                         fh_s_27, fh_s_28, ff_7, fg_12, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -f_6 * ff_s_7[k]
                  + f_2 * fh_s_26[k]
                  + f_7 * ff_7[k]
                  + pb_z[k] * fg_12[k];

        t_26[k] = f_5 * dg_8[k]
                  + f_2 * fh_s_27[k]
                  + pb_y[k] * fg_13[k];

        t_27[k] = pa_y[k] * dh_12[k]
                  + f_2 * fh_s_28[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_z, pb_z, dg_0, dg_2, dh_0, dh_4, fh_s_29, \
                         fh_s_30, fh_s_32, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_z[k] * dh_0[k]
                  + f_2 * fh_s_29[k];

        t_29[k] = f_5 * dg_0[k]
                  + f_2 * fh_s_30[k]
                  + pb_z[k] * fg_14[k];

        t_30[k] = f_7 * dg_2[k]
                  + pa_z[k] * dh_4[k]
                  + f_2 * fh_s_32[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_z, pb_x, pb_y, dg_4, dg_12, dh_8, ff_s_10, \
                         fh_s_35, fh_s_36, fh_s_37, ff_8, fg_15, \
                         fg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * dg_4[k]
                  + pa_z[k] * dh_8[k]
                  + f_2 * fh_s_35[k];

        t_32[k] = f_7 * dg_12[k]
                  + f_2 * fh_s_36[k]
                  + pb_x[k] * fg_18[k];

        t_33[k] = -f_9 * ff_s_10[k]
                  + f_2 * fh_s_37[k]
                  + f_0 * ff_8[k]
                  + pb_y[k] * fg_15[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pb_y, ff_s_11, ff_s_12, fh_s_38, fh_s_39, fh_s_40, \
                         ff_9, ff_10, fg_16, fg_17, fg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = -f_6 * ff_s_11[k]
                  + f_2 * fh_s_38[k]
                  + f_7 * ff_9[k]
                  + pb_y[k] * fg_16[k];

        t_35[k] = -f_4 * ff_s_12[k]
                  + f_2 * fh_s_39[k]
                  + f_5 * ff_10[k]
                  + pb_y[k] * fg_17[k];

        t_36[k] = f_2 * fh_s_40[k]
                  + pb_y[k] * fg_18[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_x, pb_y, ph_s_22, ph_21, dg_9, dg_13, dh_29, \
                         dh_30, fh_s_41, fh_s_42, fh_s_43, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -f_8 * ph_s_22[k]
                  + f_5 * ph_21[k]
                  + pa_x[k] * dh_29[k]
                  + f_2 * fh_s_41[k];

        t_38[k] = f_10 * dg_13[k]
                  + pa_x[k] * dh_30[k]
                  + f_2 * fh_s_42[k];

        t_39[k] = f_7 * dg_9[k]
                  + f_2 * fh_s_43[k]
                  + pb_y[k] * fg_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_x, dg_15, dg_16, dg_17, dh_32, dh_34, dh_35, \
                         fh_s_45, fh_s_46, fh_s_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * dg_15[k]
                  + pa_x[k] * dh_32[k]
                  + f_2 * fh_s_45[k];

        t_41[k] = f_0 * dg_16[k]
                  + pa_x[k] * dh_34[k]
                  + f_2 * fh_s_46[k];

        t_42[k] = f_7 * dg_17[k]
                  + pa_x[k] * dh_35[k]
                  + f_2 * fh_s_47[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pb_x, dg_19, dg_20, dh_38, dh_43, \
                         dh_45, fh_s_49, fh_s_50, fh_s_51, fh_s_52, \
                         fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_7 * dg_19[k]
                  + pa_x[k] * dh_38[k]
                  + f_2 * fh_s_49[k];

        t_44[k] = f_5 * dg_20[k]
                  + f_2 * fh_s_50[k]
                  + pb_x[k] * fg_20[k];

        t_45[k] = pa_x[k] * dh_43[k]
                  + f_2 * fh_s_51[k];

        t_46[k] = pa_x[k] * dh_45[k]
                  + f_2 * fh_s_52[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_x, pa_y, dh_23, dh_46, dh_47, dh_48, \
                         fh_s_53, fh_s_54, fh_s_55, fh_s_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pa_x[k] * dh_46[k]
                  + f_2 * fh_s_53[k];

        t_48[k] = pa_x[k] * dh_47[k]
                  + f_2 * fh_s_54[k];

        t_49[k] = pa_x[k] * dh_48[k]
                  + f_2 * fh_s_55[k];

        t_50[k] = pa_y[k] * dh_23[k]
                  + f_2 * fh_s_56[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pa_z, dh_14, dh_16, dh_24, dh_25, \
                         fh_s_57, fh_s_58, fh_s_59, fh_s_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_z[k] * dh_14[k]
                  + f_2 * fh_s_57[k];

        t_52[k] = pa_y[k] * dh_24[k]
                  + f_2 * fh_s_58[k];

        t_53[k] = pa_z[k] * dh_16[k]
                  + f_2 * fh_s_59[k];

        t_54[k] = pa_y[k] * dh_25[k]
                  + f_2 * fh_s_60[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_x, dh_51, dh_52, dh_53, dh_54, fh_s_61, \
                         fh_s_62, fh_s_63, fh_s_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_x[k] * dh_51[k]
                  + f_2 * fh_s_61[k];

        t_56[k] = pa_x[k] * dh_52[k]
                  + f_2 * fh_s_62[k];

        t_57[k] = pa_x[k] * dh_53[k]
                  + f_2 * fh_s_63[k];

        t_58[k] = pa_x[k] * dh_54[k]
                  + f_2 * fh_s_64[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pa_x, pb_z, dg_11, dg_28, dg_31, dh_56, dh_61, \
                         fh_s_65, fh_s_67, fh_s_70, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_10 * dg_28[k]
                  + pa_x[k] * dh_56[k]
                  + f_2 * fh_s_65[k];

        t_60[k] = f_7 * dg_11[k]
                  + f_2 * fh_s_67[k]
                  + pb_z[k] * fg_21[k];

        t_61[k] = f_0 * dg_31[k]
                  + pa_x[k] * dh_61[k]
                  + f_2 * fh_s_70[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pb_x, dg_34, dg_39, dh_65, dh_70, \
                         dh_71, fh_s_74, fh_s_75, fh_s_76, fh_s_77, \
                         fg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_7 * dg_34[k]
                  + pa_x[k] * dh_65[k]
                  + f_2 * fh_s_74[k];

        t_63[k] = f_5 * dg_39[k]
                  + f_2 * fh_s_75[k]
                  + pb_x[k] * fg_22[k];

        t_64[k] = pa_x[k] * dh_70[k]
                  + f_2 * fh_s_76[k];

        t_65[k] = pa_x[k] * dh_71[k]
                  + f_2 * fh_s_77[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_x, pb_x, dh_72, dh_73, dh_75, ff_s_18, \
                         fh_s_78, fh_s_79, fh_s_80, fh_s_81, ff_11, \
                         fg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_x[k] * dh_72[k]
                  + f_2 * fh_s_78[k];

        t_67[k] = pa_x[k] * dh_73[k]
                  + f_2 * fh_s_79[k];

        t_68[k] = pa_x[k] * dh_75[k]
                  + f_2 * fh_s_80[k];

        t_69[k] = -f_1 * ff_s_18[k]
                  + f_2 * fh_s_81[k]
                  + f_3 * ff_11[k]
                  + pb_x[k] * fg_23[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pb_x, ff_s_19, ff_s_20, ff_s_21, fh_s_82, fh_s_83, \
                         fh_s_84, ff_12, ff_13, ff_14, fg_24, fg_25, \
                         fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -f_9 * ff_s_19[k]
                  + f_2 * fh_s_82[k]
                  + f_0 * ff_12[k]
                  + pb_x[k] * fg_24[k];

        t_71[k] = -f_6 * ff_s_20[k]
                  + f_2 * fh_s_83[k]
                  + f_7 * ff_13[k]
                  + pb_x[k] * fg_25[k];

        t_72[k] = -f_6 * ff_s_21[k]
                  + f_2 * fh_s_84[k]
                  + f_7 * ff_14[k]
                  + pb_x[k] * fg_26[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pb_x, ff_s_22, ff_s_24, ff_s_25, fh_s_85, fh_s_86, \
                         fh_s_87, ff_15, ff_17, ff_18, fg_27, fg_28, \
                         fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = -f_4 * ff_s_22[k]
                  + f_2 * fh_s_85[k]
                  + f_5 * ff_15[k]
                  + pb_x[k] * fg_27[k];

        t_74[k] = -f_4 * ff_s_24[k]
                  + f_2 * fh_s_86[k]
                  + f_5 * ff_17[k]
                  + pb_x[k] * fg_28[k];

        t_75[k] = -f_4 * ff_s_25[k]
                  + f_2 * fh_s_87[k]
                  + f_5 * ff_18[k]
                  + pb_x[k] * fg_29[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pb_x, fh_s_88, fh_s_89, fh_s_90, fh_s_91, \
                         fg_30, fg_32, fg_33, fg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_2 * fh_s_88[k]
                  + pb_x[k] * fg_30[k];

        t_77[k] = f_2 * fh_s_89[k]
                  + pb_x[k] * fg_32[k];

        t_78[k] = f_2 * fh_s_90[k]
                  + pb_x[k] * fg_33[k];

        t_79[k] = f_2 * fh_s_91[k]
                  + pb_x[k] * fg_34[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pb_y, pb_z, dg_20, ff_s_22, fh_s_92, fh_s_93, \
                         fh_s_94, ff_15, fg_30, fg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * dg_20[k]
                  - f_1 * ff_s_22[k]
                  + f_2 * fh_s_92[k]
                  + f_3 * ff_15[k]
                  + pb_y[k] * fg_30[k];

        t_81[k] = f_2 * fh_s_93[k]
                  + pb_z[k] * fg_30[k];

        t_82[k] = -f_4 * ff_s_22[k]
                  + f_2 * fh_s_94[k]
                  + f_5 * ff_15[k]
                  + pb_z[k] * fg_31[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_y, pb_z, dg_24, ff_s_23, ff_s_25, fh_s_95, \
                         fh_s_96, fh_s_97, ff_16, ff_18, fg_32, fg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = -f_6 * ff_s_23[k]
                  + f_2 * fh_s_95[k]
                  + f_7 * ff_16[k]
                  + pb_z[k] * fg_32[k];

        t_84[k] = f_0 * dg_24[k]
                  + f_2 * fh_s_96[k]
                  + pb_y[k] * fg_34[k];

        t_85[k] = -f_1 * ff_s_25[k]
                  + f_2 * fh_s_97[k]
                  + f_3 * ff_18[k]
                  + pb_z[k] * fg_34[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pb_x, ff_s_26, ff_s_28, ff_s_32, fh_s_98, fh_s_100, \
                         fh_s_103, ff_19, ff_20, ff_21, fg_35, fg_36, \
                         fg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = -f_9 * ff_s_26[k]
                  + f_2 * fh_s_98[k]
                  + f_0 * ff_19[k]
                  + pb_x[k] * fg_35[k];

        t_87[k] = -f_6 * ff_s_28[k]
                  + f_2 * fh_s_100[k]
                  + f_7 * ff_20[k]
                  + pb_x[k] * fg_36[k];

        t_88[k] = -f_4 * ff_s_32[k]
                  + f_2 * fh_s_103[k]
                  + f_5 * ff_21[k]
                  + pb_x[k] * fg_37[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_z, pb_x, pb_z, dg_20, dh_43, fh_s_107, fh_s_108, \
                         fh_s_109, fg_38, fg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_2 * fh_s_107[k]
                  + pb_x[k] * fg_39[k];

        t_90[k] = pa_z[k] * dh_43[k]
                  + f_2 * fh_s_108[k];

        t_91[k] = f_5 * dg_20[k]
                  + f_2 * fh_s_109[k]
                  + pb_z[k] * fg_38[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, pa_z, pb_y, dg_21, dg_22, dg_27, dh_45, dh_46, \
                         fh_s_110, fh_s_111, fh_s_112, fg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_7 * dg_21[k]
                  + pa_z[k] * dh_45[k]
                  + f_2 * fh_s_110[k];

        t_93[k] = f_0 * dg_22[k]
                  + pa_z[k] * dh_46[k]
                  + f_2 * fh_s_111[k];

        t_94[k] = f_7 * dg_27[k]
                  + f_2 * fh_s_112[k]
                  + pb_y[k] * fg_39[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, pa_y, pb_z, ph_s_22, ph_21, dg_25, dg_35, dh_55, \
                         dh_70, fh_s_113, fh_s_124, fh_s_125, fg_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -f_8 * ph_s_22[k]
                  + f_5 * ph_21[k]
                  + pa_y[k] * dh_55[k]
                  + f_2 * fh_s_113[k];

        t_96[k] = f_10 * dg_35[k]
                  + pa_y[k] * dh_70[k]
                  + f_2 * fh_s_124[k];

        t_97[k] = f_7 * dg_25[k]
                  + f_2 * fh_s_125[k]
                  + pb_z[k] * fg_40[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pa_y, pb_y, dg_37, dg_38, dg_39, dh_72, dh_73, \
                         fh_s_126, fh_s_127, fh_s_128, fg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * dg_37[k]
                  + pa_y[k] * dh_72[k]
                  + f_2 * fh_s_126[k];

        t_99[k] = f_7 * dg_38[k]
                  + pa_y[k] * dh_73[k]
                  + f_2 * fh_s_127[k];

        t_100[k] = f_5 * dg_39[k]
                   + f_2 * fh_s_128[k]
                   + pb_y[k] * fg_41[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pa_y, pb_x, dh_75, ff_s_40, ff_s_41, fh_s_129, \
                         fh_s_130, fh_s_131, ff_22, ff_23, fg_42, \
                         fg_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = pa_y[k] * dh_75[k]
                   + f_2 * fh_s_129[k];

        t_102[k] = -f_1 * ff_s_40[k]
                   + f_2 * fh_s_130[k]
                   + f_3 * ff_22[k]
                   + pb_x[k] * fg_42[k];

        t_103[k] = -f_9 * ff_s_41[k]
                   + f_2 * fh_s_131[k]
                   + f_0 * ff_23[k]
                   + pb_x[k] * fg_43[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pb_x, ff_s_42, ff_s_43, ff_s_44, fh_s_132, \
                         fh_s_133, fh_s_134, ff_24, ff_25, ff_26, fg_44, fg_45, \
                         fg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -f_6 * ff_s_42[k]
                   + f_2 * fh_s_132[k]
                   + f_7 * ff_24[k]
                   + pb_x[k] * fg_44[k];

        t_105[k] = -f_6 * ff_s_43[k]
                   + f_2 * fh_s_133[k]
                   + f_7 * ff_25[k]
                   + pb_x[k] * fg_45[k];

        t_106[k] = -f_4 * ff_s_44[k]
                   + f_2 * fh_s_134[k]
                   + f_5 * ff_26[k]
                   + pb_x[k] * fg_46[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pb_x, ff_s_45, ff_s_47, fh_s_135, fh_s_136, \
                         fh_s_137, ff_27, ff_29, fg_47, fg_48, fg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -f_4 * ff_s_45[k]
                   + f_2 * fh_s_135[k]
                   + f_5 * ff_27[k]
                   + pb_x[k] * fg_47[k];

        t_108[k] = -f_4 * ff_s_47[k]
                   + f_2 * fh_s_136[k]
                   + f_5 * ff_29[k]
                   + pb_x[k] * fg_48[k];

        t_109[k] = f_2 * fh_s_137[k]
                   + pb_x[k] * fg_49[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pb_x, pb_y, ff_s_44, fh_s_138, fh_s_139, \
                         fh_s_140, fh_s_141, ff_26, fg_49, fg_50, fg_51, \
                         fg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_2 * fh_s_138[k]
                   + pb_x[k] * fg_50[k];

        t_111[k] = f_2 * fh_s_139[k]
                   + pb_x[k] * fg_51[k];

        t_112[k] = f_2 * fh_s_140[k]
                   + pb_x[k] * fg_53[k];

        t_113[k] = -f_1 * ff_s_44[k]
                   + f_2 * fh_s_141[k]
                   + f_3 * ff_26[k]
                   + pb_y[k] * fg_49[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pb_y, ff_s_45, ff_s_46, ff_s_47, fh_s_142, \
                         fh_s_143, fh_s_144, ff_27, ff_28, ff_29, fg_50, fg_51, \
                         fg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -f_9 * ff_s_45[k]
                   + f_2 * fh_s_142[k]
                   + f_0 * ff_27[k]
                   + pb_y[k] * fg_50[k];

        t_115[k] = -f_6 * ff_s_46[k]
                   + f_2 * fh_s_143[k]
                   + f_7 * ff_28[k]
                   + pb_y[k] * fg_51[k];

        t_116[k] = -f_4 * ff_s_47[k]
                   + f_2 * fh_s_144[k]
                   + f_5 * ff_29[k]
                   + pb_y[k] * fg_52[k];
    }

#pragma omp simd aligned(t_117, t_118, pb_y, pb_z, dg_39, ff_s_47, fh_s_145, fh_s_146, ff_29, \
                         fg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_2 * fh_s_145[k]
                   + pb_y[k] * fg_53[k];

        t_118[k] = f_0 * dg_39[k]
                   - f_1 * ff_s_47[k]
                   + f_2 * fh_s_146[k]
                   + f_3 * ff_29[k]
                   + pb_z[k] * fg_53[k];
    }
}

auto
compute_prim_fh_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t ph_s, const size_t ph,
                                 const size_t dg, const size_t dh, const size_t ff_s,
                                 const size_t fh_s, const size_t ff, const size_t fg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = beta / p;
    const auto f_9 = 2.5 / p;
    const auto f_10 = 3.0 * alpha / p;

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

    const auto *ph_s_0 = buffer.data(ph_s + 0);
    const auto *ph_s_3 = buffer.data(ph_s + 3);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_3 = buffer.data(ph + 3);

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
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);

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

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_31 = buffer.data(ff_s + 31);
    const auto *ff_s_36 = buffer.data(ff_s + 36);
    const auto *ff_s_37 = buffer.data(ff_s + 37);
    const auto *ff_s_38 = buffer.data(ff_s + 38);
    const auto *ff_s_39 = buffer.data(ff_s + 39);
    const auto *ff_s_40 = buffer.data(ff_s + 40);
    const auto *ff_s_41 = buffer.data(ff_s + 41);
    const auto *ff_s_42 = buffer.data(ff_s + 42);
    const auto *ff_s_43 = buffer.data(ff_s + 43);

    const auto *fh_s_0 = buffer.data(fh_s + 0);
    const auto *fh_s_1 = buffer.data(fh_s + 1);
    const auto *fh_s_2 = buffer.data(fh_s + 2);
    const auto *fh_s_3 = buffer.data(fh_s + 3);
    const auto *fh_s_4 = buffer.data(fh_s + 4);
    const auto *fh_s_5 = buffer.data(fh_s + 5);
    const auto *fh_s_6 = buffer.data(fh_s + 6);
    const auto *fh_s_7 = buffer.data(fh_s + 7);
    const auto *fh_s_8 = buffer.data(fh_s + 8);
    const auto *fh_s_9 = buffer.data(fh_s + 9);
    const auto *fh_s_10 = buffer.data(fh_s + 10);
    const auto *fh_s_11 = buffer.data(fh_s + 11);
    const auto *fh_s_12 = buffer.data(fh_s + 12);
    const auto *fh_s_13 = buffer.data(fh_s + 13);
    const auto *fh_s_14 = buffer.data(fh_s + 14);
    const auto *fh_s_15 = buffer.data(fh_s + 15);
    const auto *fh_s_16 = buffer.data(fh_s + 16);
    const auto *fh_s_17 = buffer.data(fh_s + 17);
    const auto *fh_s_18 = buffer.data(fh_s + 18);
    const auto *fh_s_19 = buffer.data(fh_s + 19);
    const auto *fh_s_20 = buffer.data(fh_s + 20);
    const auto *fh_s_21 = buffer.data(fh_s + 21);
    const auto *fh_s_22 = buffer.data(fh_s + 22);
    const auto *fh_s_23 = buffer.data(fh_s + 23);
    const auto *fh_s_24 = buffer.data(fh_s + 24);
    const auto *fh_s_25 = buffer.data(fh_s + 25);
    const auto *fh_s_26 = buffer.data(fh_s + 26);
    const auto *fh_s_27 = buffer.data(fh_s + 27);
    const auto *fh_s_28 = buffer.data(fh_s + 28);
    const auto *fh_s_29 = buffer.data(fh_s + 29);
    const auto *fh_s_30 = buffer.data(fh_s + 30);
    const auto *fh_s_31 = buffer.data(fh_s + 31);
    const auto *fh_s_32 = buffer.data(fh_s + 32);
    const auto *fh_s_33 = buffer.data(fh_s + 33);
    const auto *fh_s_34 = buffer.data(fh_s + 34);
    const auto *fh_s_35 = buffer.data(fh_s + 35);
    const auto *fh_s_36 = buffer.data(fh_s + 36);
    const auto *fh_s_37 = buffer.data(fh_s + 37);
    const auto *fh_s_38 = buffer.data(fh_s + 38);
    const auto *fh_s_39 = buffer.data(fh_s + 39);
    const auto *fh_s_40 = buffer.data(fh_s + 40);
    const auto *fh_s_41 = buffer.data(fh_s + 41);
    const auto *fh_s_42 = buffer.data(fh_s + 42);
    const auto *fh_s_43 = buffer.data(fh_s + 43);
    const auto *fh_s_44 = buffer.data(fh_s + 44);
    const auto *fh_s_45 = buffer.data(fh_s + 45);
    const auto *fh_s_46 = buffer.data(fh_s + 46);
    const auto *fh_s_47 = buffer.data(fh_s + 47);
    const auto *fh_s_48 = buffer.data(fh_s + 48);
    const auto *fh_s_49 = buffer.data(fh_s + 49);
    const auto *fh_s_50 = buffer.data(fh_s + 50);
    const auto *fh_s_51 = buffer.data(fh_s + 51);
    const auto *fh_s_52 = buffer.data(fh_s + 52);
    const auto *fh_s_53 = buffer.data(fh_s + 53);
    const auto *fh_s_54 = buffer.data(fh_s + 54);
    const auto *fh_s_55 = buffer.data(fh_s + 55);
    const auto *fh_s_56 = buffer.data(fh_s + 56);
    const auto *fh_s_57 = buffer.data(fh_s + 57);
    const auto *fh_s_58 = buffer.data(fh_s + 58);
    const auto *fh_s_59 = buffer.data(fh_s + 59);
    const auto *fh_s_60 = buffer.data(fh_s + 60);
    const auto *fh_s_61 = buffer.data(fh_s + 61);
    const auto *fh_s_62 = buffer.data(fh_s + 62);
    const auto *fh_s_63 = buffer.data(fh_s + 63);
    const auto *fh_s_64 = buffer.data(fh_s + 64);
    const auto *fh_s_65 = buffer.data(fh_s + 65);
    const auto *fh_s_66 = buffer.data(fh_s + 66);
    const auto *fh_s_67 = buffer.data(fh_s + 67);
    const auto *fh_s_68 = buffer.data(fh_s + 68);
    const auto *fh_s_69 = buffer.data(fh_s + 69);
    const auto *fh_s_70 = buffer.data(fh_s + 70);
    const auto *fh_s_71 = buffer.data(fh_s + 71);
    const auto *fh_s_72 = buffer.data(fh_s + 72);
    const auto *fh_s_73 = buffer.data(fh_s + 73);
    const auto *fh_s_74 = buffer.data(fh_s + 74);
    const auto *fh_s_75 = buffer.data(fh_s + 75);
    const auto *fh_s_76 = buffer.data(fh_s + 76);
    const auto *fh_s_77 = buffer.data(fh_s + 77);
    const auto *fh_s_78 = buffer.data(fh_s + 78);
    const auto *fh_s_79 = buffer.data(fh_s + 79);
    const auto *fh_s_80 = buffer.data(fh_s + 80);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_41 = buffer.data(ff + 41);
    const auto *ff_42 = buffer.data(ff + 42);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_69 = buffer.data(fg + 69);
    const auto *fg_70 = buffer.data(fg + 70);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_73 = buffer.data(fg + 73);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_76 = buffer.data(fg + 76);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_78 = buffer.data(fg + 78);
    const auto *fg_79 = buffer.data(fg + 79);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_81 = buffer.data(fg + 81);
    const auto *fg_82 = buffer.data(fg + 82);
    const auto *fg_83 = buffer.data(fg + 83);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dg_0, ff_s_0, fh_s_0, fh_s_1, \
                         fh_s_2, ff_0, fg_0, fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 - f_1 * ff_s_0[k]
                 + f_2 * fh_s_0[k]
                 + f_3 * ff_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_1[k]
                 + f_5 * ff_0[k]
                 + pb_y[k] * fg_1[k];

        t_2[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_2[k]
                 + f_5 * ff_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_y, pb_z, ff_s_1, ff_s_2, fh_s_3, fh_s_4, fh_s_5, \
                         ff_1, ff_2, fg_3, fg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_6 * ff_s_1[k]
                 + f_2 * fh_s_3[k]
                 + f_7 * ff_1[k]
                 + pb_y[k] * fg_3[k];

        t_4[k] = f_2 * fh_s_4[k]
                 + pb_z[k] * fg_3[k];

        t_5[k] = -f_6 * ff_s_2[k]
                 + f_2 * fh_s_5[k]
                 + f_7 * ff_2[k]
                 + pb_z[k] * fg_4[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, dg_5, dg_6, ff_s_3, fh_s_6, fh_s_7, \
                         fh_s_8, ff_3, fg_5, fg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * dg_5[k]
                 + f_2 * fh_s_6[k]
                 + pb_x[k] * fg_5[k];

        t_7[k] = f_0 * dg_6[k]
                 + f_2 * fh_s_7[k]
                 + pb_x[k] * fg_7[k];

        t_8[k] = -f_1 * ff_s_3[k]
                 + f_2 * fh_s_8[k]
                 + f_3 * ff_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_y, pb_z, dg_0, dh_0, ff_s_5, fh_s_9, \
                         fh_s_10, fh_s_11, ff_5, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_1 * ff_s_5[k]
                 + f_2 * fh_s_9[k]
                 + f_3 * ff_5[k]
                 + pb_z[k] * fg_7[k];

        t_10[k] = pa_y[k] * dh_0[k]
                  + f_2 * fh_s_10[k];

        t_11[k] = f_5 * dg_0[k]
                  + f_2 * fh_s_11[k]
                  + pb_y[k] * fg_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pb_x, dg_1, dg_3, dg_9, dh_1, dh_3, fh_s_12, \
                         fh_s_13, fh_s_14, fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * dg_1[k]
                  + pa_y[k] * dh_1[k]
                  + f_2 * fh_s_12[k];

        t_13[k] = f_0 * dg_3[k]
                  + pa_y[k] * dh_3[k]
                  + f_2 * fh_s_13[k];

        t_14[k] = f_7 * dg_9[k]
                  + f_2 * fh_s_14[k]
                  + pb_x[k] * fg_11[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_z, pb_z, ph_s_0, ph_0, dg_0, dh_0, dh_5, \
                         fh_s_15, fh_s_16, fh_s_17, fg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_8 * ph_s_0[k]
                  + f_5 * ph_0[k]
                  + pa_x[k] * dh_5[k]
                  + f_2 * fh_s_15[k];

        t_16[k] = pa_z[k] * dh_0[k]
                  + f_2 * fh_s_16[k];

        t_17[k] = f_5 * dg_0[k]
                  + f_2 * fh_s_17[k]
                  + pb_z[k] * fg_15[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_z, pb_x, dg_2, dg_4, dg_13, dh_2, dh_4, fh_s_18, \
                         fh_s_19, fh_s_20, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_7 * dg_2[k]
                  + pa_z[k] * dh_2[k]
                  + f_2 * fh_s_18[k];

        t_19[k] = f_0 * dg_4[k]
                  + pa_z[k] * dh_4[k]
                  + f_2 * fh_s_19[k];

        t_20[k] = f_7 * dg_13[k]
                  + f_2 * fh_s_20[k]
                  + pb_x[k] * fg_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pb_y, ph_s_3, ph_3, dg_7, dg_14, dh_6, dh_7, \
                         fh_s_21, fh_s_22, fh_s_23, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_8 * ph_s_3[k]
                  + f_5 * ph_3[k]
                  + pa_x[k] * dh_6[k]
                  + f_2 * fh_s_21[k];

        t_22[k] = f_9 * dg_14[k]
                  + pa_x[k] * dh_7[k]
                  + f_2 * fh_s_22[k];

        t_23[k] = f_7 * dg_7[k]
                  + f_2 * fh_s_23[k]
                  + pb_y[k] * fg_21[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pb_x, dg_15, dg_17, dg_18, dh_8, dh_9, \
                         fh_s_24, fh_s_25, fh_s_26, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * dg_15[k]
                  + pa_x[k] * dh_8[k]
                  + f_2 * fh_s_24[k];

        t_25[k] = f_7 * dg_17[k]
                  + pa_x[k] * dh_9[k]
                  + f_2 * fh_s_25[k];

        t_26[k] = f_5 * dg_18[k]
                  + f_2 * fh_s_26[k]
                  + pb_x[k] * fg_25[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, dg_28, dh_10, dh_13, dh_14, dh_16, \
                         fh_s_27, fh_s_28, fh_s_29, fh_s_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_x[k] * dh_10[k]
                  + f_2 * fh_s_27[k];

        t_28[k] = pa_x[k] * dh_13[k]
                  + f_2 * fh_s_28[k];

        t_29[k] = pa_x[k] * dh_14[k]
                  + f_2 * fh_s_29[k];

        t_30[k] = f_9 * dg_28[k]
                  + pa_x[k] * dh_16[k]
                  + f_2 * fh_s_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_x, pb_z, dg_10, dg_31, dg_33, dh_18, dh_20, \
                         fh_s_31, fh_s_32, fh_s_33, fg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_7 * dg_10[k]
                  + f_2 * fh_s_31[k]
                  + pb_z[k] * fg_35[k];

        t_32[k] = f_0 * dg_31[k]
                  + pa_x[k] * dh_18[k]
                  + f_2 * fh_s_32[k];

        t_33[k] = f_7 * dg_33[k]
                  + pa_x[k] * dh_20[k]
                  + f_2 * fh_s_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_x, pb_x, dg_38, dh_24, ff_s_21, fh_s_34, \
                         fh_s_35, fh_s_36, ff_20, fg_42, fg_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_5 * dg_38[k]
                  + f_2 * fh_s_34[k]
                  + pb_x[k] * fg_42[k];

        t_35[k] = pa_x[k] * dh_24[k]
                  + f_2 * fh_s_35[k];

        t_36[k] = -f_1 * ff_s_21[k]
                  + f_2 * fh_s_36[k]
                  + f_3 * ff_20[k]
                  + pb_x[k] * fg_43[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, pb_z, ff_s_22, ff_s_23, fh_s_37, fh_s_38, \
                         fh_s_39, ff_21, ff_22, fg_44, fg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -f_10 * ff_s_22[k]
                  + f_2 * fh_s_37[k]
                  + f_0 * ff_21[k]
                  + pb_x[k] * fg_44[k];

        t_38[k] = -f_6 * ff_s_23[k]
                  + f_2 * fh_s_38[k]
                  + f_7 * ff_22[k]
                  + pb_x[k] * fg_45[k];

        t_39[k] = f_2 * fh_s_39[k]
                  + pb_z[k] * fg_44[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_x, pb_z, ff_s_24, ff_s_25, fh_s_40, fh_s_41, \
                         fh_s_42, ff_23, ff_24, fg_45, fg_47, fg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -f_6 * ff_s_24[k]
                  + f_2 * fh_s_40[k]
                  + f_7 * ff_23[k]
                  + pb_x[k] * fg_47[k];

        t_41[k] = -f_4 * ff_s_25[k]
                  + f_2 * fh_s_41[k]
                  + f_5 * ff_24[k]
                  + pb_x[k] * fg_48[k];

        t_42[k] = f_2 * fh_s_42[k]
                  + pb_z[k] * fg_45[k];
    }

#pragma omp simd aligned(t_43, t_44, pb_x, ff_s_27, ff_s_28, fh_s_43, fh_s_44, ff_26, ff_27, \
                         fg_49, fg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -f_4 * ff_s_27[k]
                  + f_2 * fh_s_43[k]
                  + f_5 * ff_26[k]
                  + pb_x[k] * fg_49[k];

        t_44[k] = -f_4 * ff_s_28[k]
                  + f_2 * fh_s_44[k]
                  + f_5 * ff_27[k]
                  + pb_x[k] * fg_50[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pb_y, pb_z, dg_18, ff_s_25, ff_s_26, fh_s_45, \
                         fh_s_46, fh_s_47, ff_24, ff_25, fg_51, fg_52, \
                         fg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_0 * dg_18[k]
                  - f_1 * ff_s_25[k]
                  + f_2 * fh_s_45[k]
                  + f_3 * ff_24[k]
                  + pb_y[k] * fg_51[k];

        t_46[k] = -f_4 * ff_s_25[k]
                  + f_2 * fh_s_46[k]
                  + f_5 * ff_24[k]
                  + pb_z[k] * fg_52[k];

        t_47[k] = -f_6 * ff_s_26[k]
                  + f_2 * fh_s_47[k]
                  + f_7 * ff_25[k]
                  + pb_z[k] * fg_53[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_x, pb_y, pb_z, dg_22, ff_s_28, ff_s_29, fh_s_48, \
                         fh_s_49, fh_s_50, ff_27, ff_28, fg_55, fg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * dg_22[k]
                  + f_2 * fh_s_48[k]
                  + pb_y[k] * fg_55[k];

        t_49[k] = -f_1 * ff_s_28[k]
                  + f_2 * fh_s_49[k]
                  + f_3 * ff_27[k]
                  + pb_z[k] * fg_55[k];

        t_50[k] = -f_6 * ff_s_29[k]
                  + f_2 * fh_s_50[k]
                  + f_7 * ff_28[k]
                  + pb_x[k] * fg_56[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pa_z, pb_x, pb_z, dg_18, dh_10, ff_s_31, fh_s_51, \
                         fh_s_52, fh_s_53, ff_30, fg_57, fg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -f_4 * ff_s_31[k]
                  + f_2 * fh_s_51[k]
                  + f_5 * ff_30[k]
                  + pb_x[k] * fg_57[k];

        t_52[k] = pa_z[k] * dh_10[k]
                  + f_2 * fh_s_52[k];

        t_53[k] = f_5 * dg_18[k]
                  + f_2 * fh_s_53[k]
                  + pb_z[k] * fg_58[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_z, pb_y, dg_19, dg_20, dg_27, dh_11, dh_12, \
                         fh_s_54, fh_s_55, fh_s_56, fg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_7 * dg_19[k]
                  + pa_z[k] * dh_11[k]
                  + f_2 * fh_s_54[k];

        t_55[k] = f_0 * dg_20[k]
                  + pa_z[k] * dh_12[k]
                  + f_2 * fh_s_55[k];

        t_56[k] = f_7 * dg_27[k]
                  + f_2 * fh_s_56[k]
                  + pb_y[k] * fg_62[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_y, ph_s_3, ph_3, dg_29, dg_30, dh_15, dh_17, \
                         dh_19, fh_s_57, fh_s_58, fh_s_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = -f_8 * ph_s_3[k]
                  + f_5 * ph_3[k]
                  + pa_y[k] * dh_15[k]
                  + f_2 * fh_s_57[k];

        t_58[k] = f_7 * dg_29[k]
                  + pa_y[k] * dh_17[k]
                  + f_2 * fh_s_58[k];

        t_59[k] = f_0 * dg_30[k]
                  + pa_y[k] * dh_19[k]
                  + f_2 * fh_s_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_y, pb_z, dg_23, dg_34, dg_36, dh_21, dh_22, \
                         fh_s_60, fh_s_61, fh_s_62, fg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_9 * dg_34[k]
                  + pa_y[k] * dh_21[k]
                  + f_2 * fh_s_60[k];

        t_61[k] = f_7 * dg_23[k]
                  + f_2 * fh_s_61[k]
                  + pb_z[k] * fg_65[k];

        t_62[k] = f_0 * dg_36[k]
                  + pa_y[k] * dh_22[k]
                  + f_2 * fh_s_62[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pb_y, dg_37, dg_38, dh_23, dh_24, fh_s_63, \
                         fh_s_64, fh_s_65, fg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_7 * dg_37[k]
                  + pa_y[k] * dh_23[k]
                  + f_2 * fh_s_63[k];

        t_64[k] = f_5 * dg_38[k]
                  + f_2 * fh_s_64[k]
                  + pb_y[k] * fg_69[k];

        t_65[k] = pa_y[k] * dh_24[k]
                  + f_2 * fh_s_65[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_y, ff_s_36, ff_s_37, fh_s_66, fh_s_67, \
                         fh_s_68, ff_35, ff_36, fg_70, fg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -f_1 * ff_s_36[k]
                  + f_2 * fh_s_66[k]
                  + f_3 * ff_35[k]
                  + pb_x[k] * fg_70[k];

        t_67[k] = f_2 * fh_s_67[k]
                  + pb_y[k] * fg_70[k];

        t_68[k] = -f_10 * ff_s_37[k]
                  + f_2 * fh_s_68[k]
                  + f_0 * ff_36[k]
                  + pb_x[k] * fg_72[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_x, pb_y, ff_s_38, ff_s_39, fh_s_69, fh_s_70, \
                         fh_s_71, ff_37, ff_38, fg_72, fg_73, fg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = -f_6 * ff_s_38[k]
                  + f_2 * fh_s_69[k]
                  + f_7 * ff_37[k]
                  + pb_x[k] * fg_73[k];

        t_70[k] = f_2 * fh_s_70[k]
                  + pb_y[k] * fg_72[k];

        t_71[k] = -f_6 * ff_s_39[k]
                  + f_2 * fh_s_71[k]
                  + f_7 * ff_38[k]
                  + pb_x[k] * fg_75[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_x, pb_y, ff_s_40, ff_s_41, fh_s_72, fh_s_73, \
                         fh_s_74, ff_39, ff_40, fg_75, fg_76, fg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -f_4 * ff_s_40[k]
                  + f_2 * fh_s_72[k]
                  + f_5 * ff_39[k]
                  + pb_x[k] * fg_76[k];

        t_73[k] = -f_4 * ff_s_41[k]
                  + f_2 * fh_s_73[k]
                  + f_5 * ff_40[k]
                  + pb_x[k] * fg_77[k];

        t_74[k] = f_2 * fh_s_74[k]
                  + pb_y[k] * fg_75[k];
    }

#pragma omp simd aligned(t_75, t_76, pb_x, pb_y, ff_s_40, ff_s_43, fh_s_75, fh_s_76, ff_39, \
                         ff_42, fg_78, fg_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -f_4 * ff_s_43[k]
                  + f_2 * fh_s_75[k]
                  + f_5 * ff_42[k]
                  + pb_x[k] * fg_78[k];

        t_76[k] = -f_1 * ff_s_40[k]
                  + f_2 * fh_s_76[k]
                  + f_3 * ff_39[k]
                  + pb_y[k] * fg_79[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_y, ff_s_41, ff_s_42, ff_s_43, fh_s_77, fh_s_78, \
                         fh_s_79, ff_40, ff_41, ff_42, fg_80, fg_81, \
                         fg_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -f_10 * ff_s_41[k]
                  + f_2 * fh_s_77[k]
                  + f_0 * ff_40[k]
                  + pb_y[k] * fg_80[k];

        t_78[k] = -f_6 * ff_s_42[k]
                  + f_2 * fh_s_78[k]
                  + f_7 * ff_41[k]
                  + pb_y[k] * fg_81[k];

        t_79[k] = -f_4 * ff_s_43[k]
                  + f_2 * fh_s_79[k]
                  + f_5 * ff_42[k]
                  + pb_y[k] * fg_82[k];
    }

#pragma omp simd aligned(t_80, pb_z, dg_38, ff_s_43, fh_s_80, ff_42, \
                         fg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * dg_38[k]
                  - f_1 * ff_s_43[k]
                  + f_2 * fh_s_80[k]
                  + f_3 * ff_42[k]
                  + pb_z[k] * fg_83[k];
    }
}

auto
compute_prim_fh_kinetic_energy_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t ph_s, const size_t ph,
                                 const size_t dg, const size_t dh, const size_t ff_s,
                                 const size_t fh_s, const size_t ff, const size_t fg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = beta / p;
    const auto f_9 = 3.0 * alpha / p;
    const auto f_10 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_s_3 = buffer.data(ph_s + 3);
    const auto *ph_s_9 = buffer.data(ph_s + 9);

    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_9 = buffer.data(ph + 9);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_8 = buffer.data(dg + 8);
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
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_37 = buffer.data(dh + 37);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_43 = buffer.data(dh + 43);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_50 = buffer.data(dh + 50);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_55 = buffer.data(dh + 55);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_30 = buffer.data(ff_s + 30);
    const auto *ff_s_35 = buffer.data(ff_s + 35);
    const auto *ff_s_36 = buffer.data(ff_s + 36);
    const auto *ff_s_37 = buffer.data(ff_s + 37);
    const auto *ff_s_38 = buffer.data(ff_s + 38);
    const auto *ff_s_39 = buffer.data(ff_s + 39);
    const auto *ff_s_40 = buffer.data(ff_s + 40);
    const auto *ff_s_41 = buffer.data(ff_s + 41);
    const auto *ff_s_42 = buffer.data(ff_s + 42);

    const auto *fh_s_0 = buffer.data(fh_s + 0);
    const auto *fh_s_1 = buffer.data(fh_s + 1);
    const auto *fh_s_2 = buffer.data(fh_s + 2);
    const auto *fh_s_3 = buffer.data(fh_s + 3);
    const auto *fh_s_4 = buffer.data(fh_s + 4);
    const auto *fh_s_5 = buffer.data(fh_s + 5);
    const auto *fh_s_6 = buffer.data(fh_s + 6);
    const auto *fh_s_7 = buffer.data(fh_s + 7);
    const auto *fh_s_8 = buffer.data(fh_s + 8);
    const auto *fh_s_9 = buffer.data(fh_s + 9);
    const auto *fh_s_10 = buffer.data(fh_s + 10);
    const auto *fh_s_11 = buffer.data(fh_s + 11);
    const auto *fh_s_12 = buffer.data(fh_s + 12);
    const auto *fh_s_13 = buffer.data(fh_s + 13);
    const auto *fh_s_14 = buffer.data(fh_s + 14);
    const auto *fh_s_15 = buffer.data(fh_s + 15);
    const auto *fh_s_16 = buffer.data(fh_s + 16);
    const auto *fh_s_17 = buffer.data(fh_s + 17);
    const auto *fh_s_18 = buffer.data(fh_s + 18);
    const auto *fh_s_19 = buffer.data(fh_s + 19);
    const auto *fh_s_20 = buffer.data(fh_s + 20);
    const auto *fh_s_21 = buffer.data(fh_s + 21);
    const auto *fh_s_22 = buffer.data(fh_s + 22);
    const auto *fh_s_23 = buffer.data(fh_s + 23);
    const auto *fh_s_24 = buffer.data(fh_s + 24);
    const auto *fh_s_25 = buffer.data(fh_s + 25);
    const auto *fh_s_26 = buffer.data(fh_s + 26);
    const auto *fh_s_27 = buffer.data(fh_s + 27);
    const auto *fh_s_28 = buffer.data(fh_s + 28);
    const auto *fh_s_29 = buffer.data(fh_s + 29);
    const auto *fh_s_30 = buffer.data(fh_s + 30);
    const auto *fh_s_31 = buffer.data(fh_s + 31);
    const auto *fh_s_32 = buffer.data(fh_s + 32);
    const auto *fh_s_33 = buffer.data(fh_s + 33);
    const auto *fh_s_34 = buffer.data(fh_s + 34);
    const auto *fh_s_35 = buffer.data(fh_s + 35);
    const auto *fh_s_36 = buffer.data(fh_s + 36);
    const auto *fh_s_37 = buffer.data(fh_s + 37);
    const auto *fh_s_38 = buffer.data(fh_s + 38);
    const auto *fh_s_39 = buffer.data(fh_s + 39);
    const auto *fh_s_40 = buffer.data(fh_s + 40);
    const auto *fh_s_41 = buffer.data(fh_s + 41);
    const auto *fh_s_42 = buffer.data(fh_s + 42);
    const auto *fh_s_43 = buffer.data(fh_s + 43);
    const auto *fh_s_44 = buffer.data(fh_s + 44);
    const auto *fh_s_45 = buffer.data(fh_s + 45);
    const auto *fh_s_46 = buffer.data(fh_s + 46);
    const auto *fh_s_47 = buffer.data(fh_s + 47);
    const auto *fh_s_48 = buffer.data(fh_s + 48);
    const auto *fh_s_49 = buffer.data(fh_s + 49);
    const auto *fh_s_50 = buffer.data(fh_s + 50);
    const auto *fh_s_51 = buffer.data(fh_s + 51);
    const auto *fh_s_52 = buffer.data(fh_s + 52);
    const auto *fh_s_53 = buffer.data(fh_s + 53);
    const auto *fh_s_54 = buffer.data(fh_s + 54);
    const auto *fh_s_55 = buffer.data(fh_s + 55);
    const auto *fh_s_56 = buffer.data(fh_s + 56);
    const auto *fh_s_57 = buffer.data(fh_s + 57);
    const auto *fh_s_58 = buffer.data(fh_s + 58);
    const auto *fh_s_59 = buffer.data(fh_s + 59);
    const auto *fh_s_60 = buffer.data(fh_s + 60);
    const auto *fh_s_61 = buffer.data(fh_s + 61);
    const auto *fh_s_62 = buffer.data(fh_s + 62);
    const auto *fh_s_63 = buffer.data(fh_s + 63);
    const auto *fh_s_64 = buffer.data(fh_s + 64);
    const auto *fh_s_65 = buffer.data(fh_s + 65);
    const auto *fh_s_66 = buffer.data(fh_s + 66);
    const auto *fh_s_67 = buffer.data(fh_s + 67);
    const auto *fh_s_68 = buffer.data(fh_s + 68);
    const auto *fh_s_69 = buffer.data(fh_s + 69);
    const auto *fh_s_70 = buffer.data(fh_s + 70);
    const auto *fh_s_71 = buffer.data(fh_s + 71);
    const auto *fh_s_72 = buffer.data(fh_s + 72);
    const auto *fh_s_73 = buffer.data(fh_s + 73);
    const auto *fh_s_74 = buffer.data(fh_s + 74);
    const auto *fh_s_75 = buffer.data(fh_s + 75);
    const auto *fh_s_76 = buffer.data(fh_s + 76);
    const auto *fh_s_77 = buffer.data(fh_s + 77);
    const auto *fh_s_78 = buffer.data(fh_s + 78);
    const auto *fh_s_79 = buffer.data(fh_s + 79);
    const auto *fh_s_80 = buffer.data(fh_s + 80);
    const auto *fh_s_81 = buffer.data(fh_s + 81);
    const auto *fh_s_82 = buffer.data(fh_s + 82);
    const auto *fh_s_83 = buffer.data(fh_s + 83);
    const auto *fh_s_84 = buffer.data(fh_s + 84);
    const auto *fh_s_85 = buffer.data(fh_s + 85);
    const auto *fh_s_86 = buffer.data(fh_s + 86);
    const auto *fh_s_87 = buffer.data(fh_s + 87);
    const auto *fh_s_88 = buffer.data(fh_s + 88);
    const auto *fh_s_89 = buffer.data(fh_s + 89);
    const auto *fh_s_90 = buffer.data(fh_s + 90);
    const auto *fh_s_91 = buffer.data(fh_s + 91);
    const auto *fh_s_92 = buffer.data(fh_s + 92);
    const auto *fh_s_93 = buffer.data(fh_s + 93);
    const auto *fh_s_94 = buffer.data(fh_s + 94);
    const auto *fh_s_95 = buffer.data(fh_s + 95);
    const auto *fh_s_96 = buffer.data(fh_s + 96);
    const auto *fh_s_97 = buffer.data(fh_s + 97);
    const auto *fh_s_98 = buffer.data(fh_s + 98);
    const auto *fh_s_99 = buffer.data(fh_s + 99);
    const auto *fh_s_100 = buffer.data(fh_s + 100);
    const auto *fh_s_101 = buffer.data(fh_s + 101);
    const auto *fh_s_102 = buffer.data(fh_s + 102);
    const auto *fh_s_103 = buffer.data(fh_s + 103);
    const auto *fh_s_104 = buffer.data(fh_s + 104);
    const auto *fh_s_105 = buffer.data(fh_s + 105);
    const auto *fh_s_106 = buffer.data(fh_s + 106);
    const auto *fh_s_107 = buffer.data(fh_s + 107);
    const auto *fh_s_108 = buffer.data(fh_s + 108);
    const auto *fh_s_109 = buffer.data(fh_s + 109);
    const auto *fh_s_110 = buffer.data(fh_s + 110);
    const auto *fh_s_111 = buffer.data(fh_s + 111);
    const auto *fh_s_112 = buffer.data(fh_s + 112);
    const auto *fh_s_113 = buffer.data(fh_s + 113);
    const auto *fh_s_114 = buffer.data(fh_s + 114);
    const auto *fh_s_115 = buffer.data(fh_s + 115);
    const auto *fh_s_116 = buffer.data(fh_s + 116);
    const auto *fh_s_117 = buffer.data(fh_s + 117);
    const auto *fh_s_118 = buffer.data(fh_s + 118);
    const auto *fh_s_119 = buffer.data(fh_s + 119);
    const auto *fh_s_120 = buffer.data(fh_s + 120);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
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
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
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
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, dg_0, ff_s_0, fh_s_0, fh_s_1, \
                         fh_s_2, fh_s_3, ff_0, fg_0, fg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 - f_1 * ff_s_0[k]
                 + f_2 * fh_s_0[k]
                 + f_3 * ff_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = f_2 * fh_s_1[k]
                 + pb_y[k] * fg_0[k];

        t_2[k] = f_2 * fh_s_2[k]
                 + pb_z[k] * fg_0[k];

        t_3[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_3[k]
                 + f_5 * ff_0[k]
                 + pb_y[k] * fg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, ff_s_0, ff_s_1, fh_s_4, fh_s_5, fh_s_6, \
                         ff_0, ff_1, fg_2, fg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_4[k]
                 + f_5 * ff_0[k]
                 + pb_z[k] * fg_2[k];

        t_5[k] = -f_6 * ff_s_1[k]
                 + f_2 * fh_s_5[k]
                 + f_7 * ff_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = f_2 * fh_s_6[k]
                 + pb_z[k] * fg_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_y, pb_z, ff_s_2, ff_s_3, fh_s_7, fh_s_8, fh_s_9, \
                         ff_2, ff_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * fh_s_7[k]
                 + pb_y[k] * fg_4[k];

        t_8[k] = -f_6 * ff_s_2[k]
                 + f_2 * fh_s_8[k]
                 + f_7 * ff_2[k]
                 + pb_z[k] * fg_4[k];

        t_9[k] = -f_1 * ff_s_3[k]
                 + f_2 * fh_s_9[k]
                 + f_3 * ff_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, pb_z, ff_s_4, ff_s_5, fh_s_10, fh_s_11, \
                         fh_s_12, ff_4, ff_5, fg_6, fg_7, fg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_6 * ff_s_4[k]
                  + f_2 * fh_s_10[k]
                  + f_7 * ff_4[k]
                  + pb_y[k] * fg_6[k];

        t_11[k] = -f_4 * ff_s_5[k]
                  + f_2 * fh_s_11[k]
                  + f_5 * ff_5[k]
                  + pb_y[k] * fg_7[k];

        t_12[k] = -f_1 * ff_s_5[k]
                  + f_2 * fh_s_12[k]
                  + f_3 * ff_5[k]
                  + pb_z[k] * fg_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, dg_1, dg_3, dh_0, dh_3, dh_4, dh_5, \
                         fh_s_13, fh_s_14, fh_s_15, fh_s_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * dh_0[k]
                  + f_2 * fh_s_13[k];

        t_14[k] = f_7 * dg_1[k]
                  + pa_y[k] * dh_3[k]
                  + f_2 * fh_s_14[k];

        t_15[k] = pa_y[k] * dh_4[k]
                  + f_2 * fh_s_15[k];

        t_16[k] = f_0 * dg_3[k]
                  + pa_y[k] * dh_5[k]
                  + f_2 * fh_s_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, pb_z, ph_s_3, ph_3, dh_8, dh_14, \
                         ff_s_7, fh_s_17, fh_s_18, fh_s_19, ff_7, \
                         fg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = pa_y[k] * dh_8[k]
                  + f_2 * fh_s_17[k];

        t_18[k] = -f_8 * ph_s_3[k]
                  + f_5 * ph_3[k]
                  + pa_x[k] * dh_14[k]
                  + f_2 * fh_s_18[k];

        t_19[k] = -f_4 * ff_s_7[k]
                  + f_2 * fh_s_19[k]
                  + f_5 * ff_7[k]
                  + pb_z[k] * fg_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_y, pb_z, dg_6, dh_10, ff_s_8, fh_s_20, \
                         fh_s_21, fh_s_22, ff_8, fg_12, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_6 * ff_s_8[k]
                  + f_2 * fh_s_20[k]
                  + f_7 * ff_8[k]
                  + pb_z[k] * fg_12[k];

        t_21[k] = f_5 * dg_6[k]
                  + f_2 * fh_s_21[k]
                  + pb_y[k] * fg_13[k];

        t_22[k] = pa_y[k] * dh_10[k]
                  + f_2 * fh_s_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_z, pb_z, dg_0, dg_2, dh_0, dh_4, fh_s_23, \
                         fh_s_24, fh_s_25, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_z[k] * dh_0[k]
                  + f_2 * fh_s_23[k];

        t_24[k] = f_5 * dg_0[k]
                  + f_2 * fh_s_24[k]
                  + pb_z[k] * fg_14[k];

        t_25[k] = f_7 * dg_2[k]
                  + pa_z[k] * dh_4[k]
                  + f_2 * fh_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, dg_4, dh_8, ff_s_10, fh_s_26, fh_s_27, \
                         fh_s_28, ff_10, fg_15, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * fh_s_26[k]
                  + pb_y[k] * fg_15[k];

        t_27[k] = f_0 * dg_4[k]
                  + pa_z[k] * dh_8[k]
                  + f_2 * fh_s_27[k];

        t_28[k] = -f_9 * ff_s_10[k]
                  + f_2 * fh_s_28[k]
                  + f_0 * ff_10[k]
                  + pb_y[k] * fg_16[k];
    }

#pragma omp simd aligned(t_29, t_30, pb_y, ff_s_11, ff_s_12, fh_s_29, fh_s_30, ff_11, ff_12, \
                         fg_17, fg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -f_6 * ff_s_11[k]
                  + f_2 * fh_s_29[k]
                  + f_7 * ff_11[k]
                  + pb_y[k] * fg_17[k];

        t_30[k] = -f_4 * ff_s_12[k]
                  + f_2 * fh_s_30[k]
                  + f_5 * ff_12[k]
                  + pb_y[k] * fg_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_x, ph_s_9, ph_9, dg_10, dg_11, dh_19, dh_20, \
                         dh_21, fh_s_31, fh_s_32, fh_s_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -f_8 * ph_s_9[k]
                  + f_5 * ph_9[k]
                  + pa_x[k] * dh_19[k]
                  + f_2 * fh_s_31[k];

        t_32[k] = f_10 * dg_10[k]
                  + pa_x[k] * dh_20[k]
                  + f_2 * fh_s_32[k];

        t_33[k] = f_0 * dg_11[k]
                  + pa_x[k] * dh_21[k]
                  + f_2 * fh_s_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_x, dg_12, dg_13, dg_14, dh_22, dh_23, dh_25, \
                         fh_s_34, fh_s_35, fh_s_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_0 * dg_12[k]
                  + pa_x[k] * dh_22[k]
                  + f_2 * fh_s_34[k];

        t_35[k] = f_7 * dg_13[k]
                  + pa_x[k] * dh_23[k]
                  + f_2 * fh_s_35[k];

        t_36[k] = f_7 * dg_14[k]
                  + pa_x[k] * dh_25[k]
                  + f_2 * fh_s_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_x, pb_x, dg_15, dh_28, dh_30, dh_31, \
                         fh_s_37, fh_s_38, fh_s_39, fh_s_40, fg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * dg_15[k]
                  + f_2 * fh_s_37[k]
                  + pb_x[k] * fg_23[k];

        t_38[k] = pa_x[k] * dh_28[k]
                  + f_2 * fh_s_38[k];

        t_39[k] = pa_x[k] * dh_30[k]
                  + f_2 * fh_s_39[k];

        t_40[k] = pa_x[k] * dh_31[k]
                  + f_2 * fh_s_40[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pa_y, pa_z, dh_12, dh_16, dh_32, dh_33, \
                         fh_s_41, fh_s_42, fh_s_43, fh_s_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pa_x[k] * dh_32[k]
                  + f_2 * fh_s_41[k];

        t_42[k] = pa_x[k] * dh_33[k]
                  + f_2 * fh_s_42[k];

        t_43[k] = pa_y[k] * dh_16[k]
                  + f_2 * fh_s_43[k];

        t_44[k] = pa_z[k] * dh_12[k]
                  + f_2 * fh_s_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_x, pa_y, pa_z, dh_13, dh_17, dh_18, dh_35, \
                         fh_s_45, fh_s_46, fh_s_47, fh_s_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_y[k] * dh_17[k]
                  + f_2 * fh_s_45[k];

        t_46[k] = pa_z[k] * dh_13[k]
                  + f_2 * fh_s_46[k];

        t_47[k] = pa_y[k] * dh_18[k]
                  + f_2 * fh_s_47[k];

        t_48[k] = pa_x[k] * dh_35[k]
                  + f_2 * fh_s_48[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_x, dg_21, dh_36, dh_37, dh_38, dh_40, \
                         fh_s_49, fh_s_50, fh_s_51, fh_s_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pa_x[k] * dh_36[k]
                  + f_2 * fh_s_49[k];

        t_50[k] = pa_x[k] * dh_37[k]
                  + f_2 * fh_s_50[k];

        t_51[k] = pa_x[k] * dh_38[k]
                  + f_2 * fh_s_51[k];

        t_52[k] = f_10 * dg_21[k]
                  + pa_x[k] * dh_40[k]
                  + f_2 * fh_s_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_x, pb_z, dg_8, dg_24, dg_26, dh_43, dh_46, \
                         fh_s_53, fh_s_54, fh_s_55, fg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_7 * dg_8[k]
                  + f_2 * fh_s_53[k]
                  + pb_z[k] * fg_24[k];

        t_54[k] = f_0 * dg_24[k]
                  + pa_x[k] * dh_43[k]
                  + f_2 * fh_s_54[k];

        t_55[k] = f_7 * dg_26[k]
                  + pa_x[k] * dh_46[k]
                  + f_2 * fh_s_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_x, dg_31, dh_50, dh_51, dh_52, \
                         fh_s_56, fh_s_57, fh_s_58, fh_s_59, fg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_5 * dg_31[k]
                  + f_2 * fh_s_56[k]
                  + pb_x[k] * fg_27[k];

        t_57[k] = pa_x[k] * dh_50[k]
                  + f_2 * fh_s_57[k];

        t_58[k] = pa_x[k] * dh_51[k]
                  + f_2 * fh_s_58[k];

        t_59[k] = pa_x[k] * dh_52[k]
                  + f_2 * fh_s_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_x, pb_x, dh_53, dh_55, ff_s_20, fh_s_60, \
                         fh_s_61, fh_s_62, ff_17, fg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_x[k] * dh_53[k]
                  + f_2 * fh_s_60[k];

        t_61[k] = pa_x[k] * dh_55[k]
                  + f_2 * fh_s_61[k];

        t_62[k] = -f_1 * ff_s_20[k]
                  + f_2 * fh_s_62[k]
                  + f_3 * ff_17[k]
                  + pb_x[k] * fg_28[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_x, pb_z, ff_s_21, ff_s_22, fh_s_63, fh_s_64, \
                         fh_s_65, ff_18, ff_19, fg_29, fg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -f_9 * ff_s_21[k]
                  + f_2 * fh_s_63[k]
                  + f_0 * ff_18[k]
                  + pb_x[k] * fg_29[k];

        t_64[k] = -f_6 * ff_s_22[k]
                  + f_2 * fh_s_64[k]
                  + f_7 * ff_19[k]
                  + pb_x[k] * fg_30[k];

        t_65[k] = f_2 * fh_s_65[k]
                  + pb_z[k] * fg_29[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, ff_s_23, ff_s_24, fh_s_66, fh_s_67, \
                         fh_s_68, ff_20, ff_21, fg_30, fg_31, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -f_6 * ff_s_23[k]
                  + f_2 * fh_s_66[k]
                  + f_7 * ff_20[k]
                  + pb_x[k] * fg_31[k];

        t_67[k] = -f_4 * ff_s_24[k]
                  + f_2 * fh_s_67[k]
                  + f_5 * ff_21[k]
                  + pb_x[k] * fg_32[k];

        t_68[k] = f_2 * fh_s_68[k]
                  + pb_z[k] * fg_30[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_x, ff_s_26, ff_s_27, fh_s_69, fh_s_70, fh_s_71, \
                         ff_23, ff_24, fg_33, fg_34, fg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = -f_4 * ff_s_26[k]
                  + f_2 * fh_s_69[k]
                  + f_5 * ff_23[k]
                  + pb_x[k] * fg_33[k];

        t_70[k] = -f_4 * ff_s_27[k]
                  + f_2 * fh_s_70[k]
                  + f_5 * ff_24[k]
                  + pb_x[k] * fg_34[k];

        t_71[k] = f_2 * fh_s_71[k]
                  + pb_x[k] * fg_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_x, fh_s_72, fh_s_73, fh_s_74, fg_37, fg_38, \
                         fg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_2 * fh_s_72[k]
                  + pb_x[k] * fg_37[k];

        t_73[k] = f_2 * fh_s_73[k]
                  + pb_x[k] * fg_38[k];

        t_74[k] = f_2 * fh_s_74[k]
                  + pb_x[k] * fg_39[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_y, pb_z, dg_15, ff_s_24, fh_s_75, fh_s_76, \
                         fh_s_77, ff_21, fg_35, fg_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_0 * dg_15[k]
                  - f_1 * ff_s_24[k]
                  + f_2 * fh_s_75[k]
                  + f_3 * ff_21[k]
                  + pb_y[k] * fg_35[k];

        t_76[k] = f_2 * fh_s_76[k]
                  + pb_z[k] * fg_35[k];

        t_77[k] = -f_4 * ff_s_24[k]
                  + f_2 * fh_s_77[k]
                  + f_5 * ff_21[k]
                  + pb_z[k] * fg_36[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pb_y, pb_z, dg_18, ff_s_25, ff_s_27, fh_s_78, \
                         fh_s_79, fh_s_80, ff_22, ff_24, fg_37, fg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = -f_6 * ff_s_25[k]
                  + f_2 * fh_s_78[k]
                  + f_7 * ff_22[k]
                  + pb_z[k] * fg_37[k];

        t_79[k] = f_0 * dg_18[k]
                  + f_2 * fh_s_79[k]
                  + pb_y[k] * fg_39[k];

        t_80[k] = -f_1 * ff_s_27[k]
                  + f_2 * fh_s_80[k]
                  + f_3 * ff_24[k]
                  + pb_z[k] * fg_39[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_x, ff_s_28, ff_s_30, fh_s_81, fh_s_82, fh_s_83, \
                         ff_25, ff_27, fg_40, fg_41, fg_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -f_6 * ff_s_28[k]
                  + f_2 * fh_s_81[k]
                  + f_7 * ff_25[k]
                  + pb_x[k] * fg_40[k];

        t_82[k] = -f_4 * ff_s_30[k]
                  + f_2 * fh_s_82[k]
                  + f_5 * ff_27[k]
                  + pb_x[k] * fg_41[k];

        t_83[k] = f_2 * fh_s_83[k]
                  + pb_x[k] * fg_43[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_z, pb_x, pb_z, dg_15, dh_28, fh_s_84, fh_s_85, \
                         fh_s_86, fg_42, fg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_2 * fh_s_84[k]
                  + pb_x[k] * fg_44[k];

        t_85[k] = pa_z[k] * dh_28[k]
                  + f_2 * fh_s_85[k];

        t_86[k] = f_5 * dg_15[k]
                  + f_2 * fh_s_86[k]
                  + pb_z[k] * fg_42[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pa_z, pb_y, dg_16, dg_17, dg_20, dh_30, dh_31, \
                         fh_s_87, fh_s_88, fh_s_89, fg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_7 * dg_16[k]
                  + pa_z[k] * dh_30[k]
                  + f_2 * fh_s_87[k];

        t_88[k] = f_0 * dg_17[k]
                  + pa_z[k] * dh_31[k]
                  + f_2 * fh_s_88[k];

        t_89[k] = f_7 * dg_20[k]
                  + f_2 * fh_s_89[k]
                  + pb_y[k] * fg_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pa_y, ph_s_9, ph_9, dg_22, dg_23, dh_39, dh_42, \
                         dh_44, fh_s_90, fh_s_91, fh_s_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -f_8 * ph_s_9[k]
                  + f_5 * ph_9[k]
                  + pa_y[k] * dh_39[k]
                  + f_2 * fh_s_90[k];

        t_91[k] = f_7 * dg_22[k]
                  + pa_y[k] * dh_42[k]
                  + f_2 * fh_s_91[k];

        t_92[k] = f_0 * dg_23[k]
                  + pa_y[k] * dh_44[k]
                  + f_2 * fh_s_92[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pa_y, pb_x, pb_z, dg_19, dg_27, dh_50, \
                         fh_s_93, fh_s_94, fh_s_95, fh_s_96, fg_47, \
                         fg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_2 * fh_s_93[k]
                  + pb_x[k] * fg_47[k];

        t_94[k] = f_2 * fh_s_94[k]
                  + pb_x[k] * fg_48[k];

        t_95[k] = f_10 * dg_27[k]
                  + pa_y[k] * dh_50[k]
                  + f_2 * fh_s_95[k];

        t_96[k] = f_7 * dg_19[k]
                  + f_2 * fh_s_96[k]
                  + pb_z[k] * fg_47[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pa_y, pb_y, dg_29, dg_30, dg_31, dh_52, dh_53, \
                         fh_s_97, fh_s_98, fh_s_99, fg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_0 * dg_29[k]
                  + pa_y[k] * dh_52[k]
                  + f_2 * fh_s_97[k];

        t_98[k] = f_7 * dg_30[k]
                  + pa_y[k] * dh_53[k]
                  + f_2 * fh_s_98[k];

        t_99[k] = f_5 * dg_31[k]
                  + f_2 * fh_s_99[k]
                  + pb_y[k] * fg_50[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_y, pb_x, pb_y, dh_55, ff_s_35, fh_s_100, \
                         fh_s_101, fh_s_102, ff_31, fg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_y[k] * dh_55[k]
                   + f_2 * fh_s_100[k];

        t_101[k] = -f_1 * ff_s_35[k]
                   + f_2 * fh_s_101[k]
                   + f_3 * ff_31[k]
                   + pb_x[k] * fg_51[k];

        t_102[k] = f_2 * fh_s_102[k]
                   + pb_y[k] * fg_51[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, pb_x, pb_y, ff_s_36, ff_s_37, fh_s_103, \
                         fh_s_104, fh_s_105, ff_32, ff_33, fg_52, \
                         fg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = -f_9 * ff_s_36[k]
                   + f_2 * fh_s_103[k]
                   + f_0 * ff_32[k]
                   + pb_x[k] * fg_52[k];

        t_104[k] = -f_6 * ff_s_37[k]
                   + f_2 * fh_s_104[k]
                   + f_7 * ff_33[k]
                   + pb_x[k] * fg_53[k];

        t_105[k] = f_2 * fh_s_105[k]
                   + pb_y[k] * fg_52[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pb_x, ff_s_38, ff_s_39, ff_s_40, fh_s_106, \
                         fh_s_107, fh_s_108, ff_34, ff_35, ff_36, fg_54, fg_55, \
                         fg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = -f_6 * ff_s_38[k]
                   + f_2 * fh_s_106[k]
                   + f_7 * ff_34[k]
                   + pb_x[k] * fg_54[k];

        t_107[k] = -f_4 * ff_s_39[k]
                   + f_2 * fh_s_107[k]
                   + f_5 * ff_35[k]
                   + pb_x[k] * fg_55[k];

        t_108[k] = -f_4 * ff_s_40[k]
                   + f_2 * fh_s_108[k]
                   + f_5 * ff_36[k]
                   + pb_x[k] * fg_56[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pb_x, pb_y, ff_s_42, fh_s_109, fh_s_110, \
                         fh_s_111, fh_s_112, ff_38, fg_54, fg_57, fg_58, \
                         fg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_2 * fh_s_109[k]
                   + pb_y[k] * fg_54[k];

        t_110[k] = -f_4 * ff_s_42[k]
                   + f_2 * fh_s_110[k]
                   + f_5 * ff_38[k]
                   + pb_x[k] * fg_57[k];

        t_111[k] = f_2 * fh_s_111[k]
                   + pb_x[k] * fg_58[k];

        t_112[k] = f_2 * fh_s_112[k]
                   + pb_x[k] * fg_59[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pb_x, pb_y, ff_s_39, fh_s_113, fh_s_114, \
                         fh_s_115, ff_35, fg_58, fg_60, fg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_2 * fh_s_113[k]
                   + pb_x[k] * fg_60[k];

        t_114[k] = f_2 * fh_s_114[k]
                   + pb_x[k] * fg_62[k];

        t_115[k] = -f_1 * ff_s_39[k]
                   + f_2 * fh_s_115[k]
                   + f_3 * ff_35[k]
                   + pb_y[k] * fg_58[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pb_y, ff_s_40, ff_s_41, ff_s_42, fh_s_116, \
                         fh_s_117, fh_s_118, ff_36, ff_37, ff_38, fg_59, fg_60, \
                         fg_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = -f_9 * ff_s_40[k]
                   + f_2 * fh_s_116[k]
                   + f_0 * ff_36[k]
                   + pb_y[k] * fg_59[k];

        t_117[k] = -f_6 * ff_s_41[k]
                   + f_2 * fh_s_117[k]
                   + f_7 * ff_37[k]
                   + pb_y[k] * fg_60[k];

        t_118[k] = -f_4 * ff_s_42[k]
                   + f_2 * fh_s_118[k]
                   + f_5 * ff_38[k]
                   + pb_y[k] * fg_61[k];
    }

#pragma omp simd aligned(t_119, t_120, pb_y, pb_z, dg_31, ff_s_42, fh_s_119, fh_s_120, ff_38, \
                         fg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_2 * fh_s_119[k]
                   + pb_y[k] * fg_62[k];

        t_120[k] = f_0 * dg_31[k]
                   - f_1 * ff_s_42[k]
                   + f_2 * fh_s_120[k]
                   + f_3 * ff_38[k]
                   + pb_z[k] * fg_62[k];
    }
}

auto
compute_prim_fh_kinetic_energy_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t ph_s, const size_t ph,
                                 const size_t dg, const size_t dh, const size_t ff_s,
                                 const size_t fh_s, const size_t ff, const size_t fg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = beta / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_s_0 = buffer.data(ph_s + 0);
    const auto *ph_s_3 = buffer.data(ph_s + 3);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_3 = buffer.data(ph + 3);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);

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

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_31 = buffer.data(ff_s + 31);
    const auto *ff_s_32 = buffer.data(ff_s + 32);
    const auto *ff_s_33 = buffer.data(ff_s + 33);
    const auto *ff_s_35 = buffer.data(ff_s + 35);
    const auto *ff_s_36 = buffer.data(ff_s + 36);

    const auto *fh_s_0 = buffer.data(fh_s + 0);
    const auto *fh_s_1 = buffer.data(fh_s + 1);
    const auto *fh_s_2 = buffer.data(fh_s + 2);
    const auto *fh_s_3 = buffer.data(fh_s + 3);
    const auto *fh_s_4 = buffer.data(fh_s + 4);
    const auto *fh_s_5 = buffer.data(fh_s + 5);
    const auto *fh_s_6 = buffer.data(fh_s + 6);
    const auto *fh_s_7 = buffer.data(fh_s + 7);
    const auto *fh_s_8 = buffer.data(fh_s + 8);
    const auto *fh_s_9 = buffer.data(fh_s + 9);
    const auto *fh_s_10 = buffer.data(fh_s + 10);
    const auto *fh_s_11 = buffer.data(fh_s + 11);
    const auto *fh_s_12 = buffer.data(fh_s + 12);
    const auto *fh_s_13 = buffer.data(fh_s + 13);
    const auto *fh_s_14 = buffer.data(fh_s + 14);
    const auto *fh_s_15 = buffer.data(fh_s + 15);
    const auto *fh_s_16 = buffer.data(fh_s + 16);
    const auto *fh_s_17 = buffer.data(fh_s + 17);
    const auto *fh_s_18 = buffer.data(fh_s + 18);
    const auto *fh_s_19 = buffer.data(fh_s + 19);
    const auto *fh_s_20 = buffer.data(fh_s + 20);
    const auto *fh_s_21 = buffer.data(fh_s + 21);
    const auto *fh_s_22 = buffer.data(fh_s + 22);
    const auto *fh_s_23 = buffer.data(fh_s + 23);
    const auto *fh_s_24 = buffer.data(fh_s + 24);
    const auto *fh_s_25 = buffer.data(fh_s + 25);
    const auto *fh_s_26 = buffer.data(fh_s + 26);
    const auto *fh_s_27 = buffer.data(fh_s + 27);
    const auto *fh_s_28 = buffer.data(fh_s + 28);
    const auto *fh_s_29 = buffer.data(fh_s + 29);
    const auto *fh_s_30 = buffer.data(fh_s + 30);
    const auto *fh_s_31 = buffer.data(fh_s + 31);
    const auto *fh_s_32 = buffer.data(fh_s + 32);
    const auto *fh_s_33 = buffer.data(fh_s + 33);
    const auto *fh_s_34 = buffer.data(fh_s + 34);
    const auto *fh_s_35 = buffer.data(fh_s + 35);
    const auto *fh_s_36 = buffer.data(fh_s + 36);
    const auto *fh_s_37 = buffer.data(fh_s + 37);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dg_0, ff_s_0, fh_s_0, fh_s_1, \
                         fh_s_2, ff_0, fg_0, fg_1, fg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 - f_1 * ff_s_0[k]
                 + f_2 * fh_s_0[k]
                 + f_3 * ff_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_1[k]
                 + f_5 * ff_0[k]
                 + pb_y[k] * fg_1[k];

        t_2[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_2[k]
                 + f_5 * ff_0[k]
                 + pb_z[k] * fg_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_y, pb_y, pb_z, dh_0, ff_s_1, ff_s_2, fh_s_3, \
                         fh_s_4, fh_s_5, ff_1, ff_2, fg_3, fg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_6 * ff_s_1[k]
                 + f_2 * fh_s_3[k]
                 + f_7 * ff_1[k]
                 + pb_y[k] * fg_3[k];

        t_4[k] = -f_6 * ff_s_2[k]
                 + f_2 * fh_s_4[k]
                 + f_7 * ff_2[k]
                 + pb_z[k] * fg_4[k];

        t_5[k] = pa_y[k] * dh_0[k]
                 + f_2 * fh_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, ph_s_0, ph_0, dg_1, dh_0, dh_1, dh_3, \
                         fh_s_6, fh_s_7, fh_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_8 * ph_s_0[k]
                 + f_5 * ph_0[k]
                 + pa_x[k] * dh_3[k]
                 + f_2 * fh_s_6[k];

        t_7[k] = pa_z[k] * dh_0[k]
                 + f_2 * fh_s_7[k];

        t_8[k] = f_7 * dg_1[k]
                 + pa_z[k] * dh_1[k]
                 + f_2 * fh_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_z, ph_s_3, ph_3, dg_3, dh_2, dh_4, dh_5, \
                         fh_s_9, fh_s_10, fh_s_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * dg_3[k]
                 + pa_z[k] * dh_2[k]
                 + f_2 * fh_s_9[k];

        t_10[k] = -f_8 * ph_s_3[k]
                  + f_5 * ph_3[k]
                  + pa_x[k] * dh_4[k]
                  + f_2 * fh_s_10[k];

        t_11[k] = pa_x[k] * dh_5[k]
                  + f_2 * fh_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, dh_8, dh_9, dh_14, ff_s_16, \
                         fh_s_12, fh_s_13, fh_s_14, fh_s_15, ff_16, \
                         fg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_x[k] * dh_8[k]
                  + f_2 * fh_s_12[k];

        t_13[k] = pa_x[k] * dh_9[k]
                  + f_2 * fh_s_13[k];

        t_14[k] = pa_x[k] * dh_14[k]
                  + f_2 * fh_s_14[k];

        t_15[k] = -f_1 * ff_s_16[k]
                  + f_2 * fh_s_15[k]
                  + f_3 * ff_16[k]
                  + pb_x[k] * fg_24[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, pb_y, dg_8, ff_s_17, ff_s_19, fh_s_16, \
                         fh_s_17, fh_s_18, ff_17, ff_19, fg_25, fg_27, \
                         fg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_6 * ff_s_17[k]
                  + f_2 * fh_s_16[k]
                  + f_7 * ff_17[k]
                  + pb_x[k] * fg_25[k];

        t_17[k] = -f_4 * ff_s_19[k]
                  + f_2 * fh_s_17[k]
                  + f_5 * ff_19[k]
                  + pb_x[k] * fg_27[k];

        t_18[k] = f_0 * dg_8[k]
                  - f_1 * ff_s_19[k]
                  + f_2 * fh_s_18[k]
                  + f_3 * ff_19[k]
                  + pb_y[k] * fg_28[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_z, dh_5, ff_s_19, ff_s_20, fh_s_19, \
                         fh_s_20, fh_s_21, ff_19, ff_20, fg_29, fg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_4 * ff_s_19[k]
                  + f_2 * fh_s_19[k]
                  + f_5 * ff_19[k]
                  + pb_z[k] * fg_29[k];

        t_20[k] = -f_6 * ff_s_20[k]
                  + f_2 * fh_s_20[k]
                  + f_7 * ff_20[k]
                  + pb_z[k] * fg_30[k];

        t_21[k] = pa_z[k] * dh_5[k]
                  + f_2 * fh_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pa_z, ph_s_3, ph_3, dg_9, dg_10, dh_6, dh_7, \
                         dh_10, fh_s_22, fh_s_23, fh_s_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_7 * dg_9[k]
                  + pa_z[k] * dh_6[k]
                  + f_2 * fh_s_22[k];

        t_23[k] = f_0 * dg_10[k]
                  + pa_z[k] * dh_7[k]
                  + f_2 * fh_s_23[k];

        t_24[k] = -f_8 * ph_s_3[k]
                  + f_5 * ph_3[k]
                  + pa_y[k] * dh_10[k]
                  + f_2 * fh_s_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, dg_16, dg_17, dg_18, dh_11, dh_12, \
                         dh_13, dh_14, fh_s_25, fh_s_26, fh_s_27, \
                         fh_s_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_9 * dg_16[k]
                  + pa_y[k] * dh_11[k]
                  + f_2 * fh_s_25[k];

        t_26[k] = f_0 * dg_17[k]
                  + pa_y[k] * dh_12[k]
                  + f_2 * fh_s_26[k];

        t_27[k] = f_7 * dg_18[k]
                  + pa_y[k] * dh_13[k]
                  + f_2 * fh_s_27[k];

        t_28[k] = pa_y[k] * dh_14[k]
                  + f_2 * fh_s_28[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, ff_s_29, ff_s_31, ff_s_32, fh_s_29, fh_s_30, \
                         fh_s_31, ff_29, ff_31, ff_32, fg_45, fg_47, \
                         fg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -f_1 * ff_s_29[k]
                  + f_2 * fh_s_29[k]
                  + f_3 * ff_29[k]
                  + pb_x[k] * fg_45[k];

        t_30[k] = -f_6 * ff_s_31[k]
                  + f_2 * fh_s_30[k]
                  + f_7 * ff_31[k]
                  + pb_x[k] * fg_47[k];

        t_31[k] = -f_6 * ff_s_32[k]
                  + f_2 * fh_s_31[k]
                  + f_7 * ff_32[k]
                  + pb_x[k] * fg_48[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, ff_s_33, ff_s_36, fh_s_32, fh_s_33, \
                         fh_s_34, ff_33, ff_36, fg_49, fg_50, fg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -f_4 * ff_s_33[k]
                  + f_2 * fh_s_32[k]
                  + f_5 * ff_33[k]
                  + pb_x[k] * fg_49[k];

        t_33[k] = -f_4 * ff_s_36[k]
                  + f_2 * fh_s_33[k]
                  + f_5 * ff_36[k]
                  + pb_x[k] * fg_50[k];

        t_34[k] = -f_1 * ff_s_33[k]
                  + f_2 * fh_s_34[k]
                  + f_3 * ff_33[k]
                  + pb_y[k] * fg_51[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_y, pb_z, dg_19, ff_s_35, ff_s_36, fh_s_35, \
                         fh_s_36, fh_s_37, ff_35, ff_36, fg_53, fg_54, \
                         fg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -f_6 * ff_s_35[k]
                  + f_2 * fh_s_35[k]
                  + f_7 * ff_35[k]
                  + pb_y[k] * fg_53[k];

        t_36[k] = -f_4 * ff_s_36[k]
                  + f_2 * fh_s_36[k]
                  + f_5 * ff_36[k]
                  + pb_y[k] * fg_54[k];

        t_37[k] = f_0 * dg_19[k]
                  - f_1 * ff_s_36[k]
                  + f_2 * fh_s_37[k]
                  + f_3 * ff_36[k]
                  + pb_z[k] * fg_55[k];
    }
}

auto
compute_prim_fh_kinetic_energy_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t ph_s, const size_t ph,
                                 const size_t dg, const size_t dh, const size_t ff_s,
                                 const size_t fh_s, const size_t ff, const size_t fg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = beta / p;
    const auto f_9 = 2.5 / p;
    const auto f_10 = 3.0 * alpha / p;

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

    const auto *ph_s_1 = buffer.data(ph_s + 1);
    const auto *ph_s_5 = buffer.data(ph_s + 5);

    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_5 = buffer.data(ph + 5);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_7 = buffer.data(dg + 7);
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
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);

    const auto *dh_0 = buffer.data(dh + 0);
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
    const auto *dh_27 = buffer.data(dh + 27);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_15 = buffer.data(ff_s + 15);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_30 = buffer.data(ff_s + 30);
    const auto *ff_s_31 = buffer.data(ff_s + 31);
    const auto *ff_s_32 = buffer.data(ff_s + 32);
    const auto *ff_s_33 = buffer.data(ff_s + 33);

    const auto *fh_s_0 = buffer.data(fh_s + 0);
    const auto *fh_s_1 = buffer.data(fh_s + 1);
    const auto *fh_s_2 = buffer.data(fh_s + 2);
    const auto *fh_s_3 = buffer.data(fh_s + 3);
    const auto *fh_s_4 = buffer.data(fh_s + 4);
    const auto *fh_s_5 = buffer.data(fh_s + 5);
    const auto *fh_s_6 = buffer.data(fh_s + 6);
    const auto *fh_s_7 = buffer.data(fh_s + 7);
    const auto *fh_s_8 = buffer.data(fh_s + 8);
    const auto *fh_s_9 = buffer.data(fh_s + 9);
    const auto *fh_s_10 = buffer.data(fh_s + 10);
    const auto *fh_s_11 = buffer.data(fh_s + 11);
    const auto *fh_s_12 = buffer.data(fh_s + 12);
    const auto *fh_s_13 = buffer.data(fh_s + 13);
    const auto *fh_s_14 = buffer.data(fh_s + 14);
    const auto *fh_s_15 = buffer.data(fh_s + 15);
    const auto *fh_s_16 = buffer.data(fh_s + 16);
    const auto *fh_s_17 = buffer.data(fh_s + 17);
    const auto *fh_s_18 = buffer.data(fh_s + 18);
    const auto *fh_s_19 = buffer.data(fh_s + 19);
    const auto *fh_s_20 = buffer.data(fh_s + 20);
    const auto *fh_s_21 = buffer.data(fh_s + 21);
    const auto *fh_s_22 = buffer.data(fh_s + 22);
    const auto *fh_s_23 = buffer.data(fh_s + 23);
    const auto *fh_s_24 = buffer.data(fh_s + 24);
    const auto *fh_s_25 = buffer.data(fh_s + 25);
    const auto *fh_s_26 = buffer.data(fh_s + 26);
    const auto *fh_s_27 = buffer.data(fh_s + 27);
    const auto *fh_s_28 = buffer.data(fh_s + 28);
    const auto *fh_s_29 = buffer.data(fh_s + 29);
    const auto *fh_s_30 = buffer.data(fh_s + 30);
    const auto *fh_s_31 = buffer.data(fh_s + 31);
    const auto *fh_s_32 = buffer.data(fh_s + 32);
    const auto *fh_s_33 = buffer.data(fh_s + 33);
    const auto *fh_s_34 = buffer.data(fh_s + 34);
    const auto *fh_s_35 = buffer.data(fh_s + 35);
    const auto *fh_s_36 = buffer.data(fh_s + 36);
    const auto *fh_s_37 = buffer.data(fh_s + 37);
    const auto *fh_s_38 = buffer.data(fh_s + 38);
    const auto *fh_s_39 = buffer.data(fh_s + 39);
    const auto *fh_s_40 = buffer.data(fh_s + 40);
    const auto *fh_s_41 = buffer.data(fh_s + 41);
    const auto *fh_s_42 = buffer.data(fh_s + 42);
    const auto *fh_s_43 = buffer.data(fh_s + 43);
    const auto *fh_s_44 = buffer.data(fh_s + 44);
    const auto *fh_s_45 = buffer.data(fh_s + 45);
    const auto *fh_s_46 = buffer.data(fh_s + 46);
    const auto *fh_s_47 = buffer.data(fh_s + 47);
    const auto *fh_s_48 = buffer.data(fh_s + 48);
    const auto *fh_s_49 = buffer.data(fh_s + 49);
    const auto *fh_s_50 = buffer.data(fh_s + 50);
    const auto *fh_s_51 = buffer.data(fh_s + 51);
    const auto *fh_s_52 = buffer.data(fh_s + 52);
    const auto *fh_s_53 = buffer.data(fh_s + 53);
    const auto *fh_s_54 = buffer.data(fh_s + 54);
    const auto *fh_s_55 = buffer.data(fh_s + 55);
    const auto *fh_s_56 = buffer.data(fh_s + 56);
    const auto *fh_s_57 = buffer.data(fh_s + 57);
    const auto *fh_s_58 = buffer.data(fh_s + 58);
    const auto *fh_s_59 = buffer.data(fh_s + 59);
    const auto *fh_s_60 = buffer.data(fh_s + 60);
    const auto *fh_s_61 = buffer.data(fh_s + 61);
    const auto *fh_s_62 = buffer.data(fh_s + 62);
    const auto *fh_s_63 = buffer.data(fh_s + 63);
    const auto *fh_s_64 = buffer.data(fh_s + 64);
    const auto *fh_s_65 = buffer.data(fh_s + 65);
    const auto *fh_s_66 = buffer.data(fh_s + 66);
    const auto *fh_s_67 = buffer.data(fh_s + 67);
    const auto *fh_s_68 = buffer.data(fh_s + 68);
    const auto *fh_s_69 = buffer.data(fh_s + 69);
    const auto *fh_s_70 = buffer.data(fh_s + 70);
    const auto *fh_s_71 = buffer.data(fh_s + 71);
    const auto *fh_s_72 = buffer.data(fh_s + 72);
    const auto *fh_s_73 = buffer.data(fh_s + 73);
    const auto *fh_s_74 = buffer.data(fh_s + 74);
    const auto *fh_s_75 = buffer.data(fh_s + 75);
    const auto *fh_s_76 = buffer.data(fh_s + 76);
    const auto *fh_s_77 = buffer.data(fh_s + 77);
    const auto *fh_s_78 = buffer.data(fh_s + 78);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_25 = buffer.data(ff + 25);
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
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
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
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, dg_0, ff_s_0, fh_s_0, fh_s_1, \
                         fh_s_2, fh_s_3, ff_0, fg_0, fg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 - f_1 * ff_s_0[k]
                 + f_2 * fh_s_0[k]
                 + f_3 * ff_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = f_2 * fh_s_1[k]
                 + pb_y[k] * fg_0[k];

        t_2[k] = f_2 * fh_s_2[k]
                 + pb_z[k] * fg_0[k];

        t_3[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_3[k]
                 + f_5 * ff_0[k]
                 + pb_y[k] * fg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, ff_s_0, ff_s_1, fh_s_4, fh_s_5, fh_s_6, \
                         ff_0, ff_1, fg_2, fg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_4[k]
                 + f_5 * ff_0[k]
                 + pb_z[k] * fg_2[k];

        t_5[k] = -f_6 * ff_s_1[k]
                 + f_2 * fh_s_5[k]
                 + f_7 * ff_1[k]
                 + pb_y[k] * fg_3[k];

        t_6[k] = f_2 * fh_s_6[k]
                 + pb_z[k] * fg_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_y, pb_z, ff_s_2, ff_s_3, fh_s_7, fh_s_8, fh_s_9, \
                         ff_2, ff_3, fg_4, fg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * fh_s_7[k]
                 + pb_y[k] * fg_4[k];

        t_8[k] = -f_6 * ff_s_2[k]
                 + f_2 * fh_s_8[k]
                 + f_7 * ff_2[k]
                 + pb_z[k] * fg_4[k];

        t_9[k] = -f_1 * ff_s_3[k]
                 + f_2 * fh_s_9[k]
                 + f_3 * ff_3[k]
                 + pb_y[k] * fg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pb_z, dg_1, dh_0, dh_2, ff_s_4, fh_s_10, \
                         fh_s_11, fh_s_12, ff_4, fg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * ff_s_4[k]
                  + f_2 * fh_s_10[k]
                  + f_3 * ff_4[k]
                  + pb_z[k] * fg_6[k];

        t_11[k] = pa_y[k] * dh_0[k]
                  + f_2 * fh_s_11[k];

        t_12[k] = f_7 * dg_1[k]
                  + pa_y[k] * dh_2[k]
                  + f_2 * fh_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pa_y, pa_z, ph_s_1, ph_1, dg_3, dh_0, dh_4, \
                         dh_6, fh_s_13, fh_s_14, fh_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dg_3[k]
                  + pa_y[k] * dh_4[k]
                  + f_2 * fh_s_13[k];

        t_14[k] = -f_8 * ph_s_1[k]
                  + f_5 * ph_1[k]
                  + pa_x[k] * dh_6[k]
                  + f_2 * fh_s_14[k];

        t_15[k] = pa_z[k] * dh_0[k]
                  + f_2 * fh_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_z, pb_z, dg_0, dg_2, dg_4, dh_3, dh_5, fh_s_16, \
                         fh_s_17, fh_s_18, fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * dg_0[k]
                  + f_2 * fh_s_16[k]
                  + pb_z[k] * fg_10[k];

        t_17[k] = f_7 * dg_2[k]
                  + pa_z[k] * dh_3[k]
                  + f_2 * fh_s_17[k];

        t_18[k] = f_0 * dg_4[k]
                  + pa_z[k] * dh_5[k]
                  + f_2 * fh_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, ph_s_5, ph_5, dg_9, dg_10, dh_7, dh_8, dh_9, \
                         fh_s_19, fh_s_20, fh_s_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_8 * ph_s_5[k]
                  + f_5 * ph_5[k]
                  + pa_x[k] * dh_7[k]
                  + f_2 * fh_s_19[k];

        t_20[k] = f_9 * dg_9[k]
                  + pa_x[k] * dh_8[k]
                  + f_2 * fh_s_20[k];

        t_21[k] = f_0 * dg_10[k]
                  + pa_x[k] * dh_9[k]
                  + f_2 * fh_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, dg_11, dg_12, dh_10, dh_11, \
                         dh_15, fh_s_22, fh_s_23, fh_s_24, fh_s_25, \
                         fg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_7 * dg_11[k]
                  + pa_x[k] * dh_10[k]
                  + f_2 * fh_s_22[k];

        t_23[k] = f_5 * dg_12[k]
                  + f_2 * fh_s_23[k]
                  + pb_x[k] * fg_15[k];

        t_24[k] = pa_x[k] * dh_11[k]
                  + f_2 * fh_s_24[k];

        t_25[k] = pa_x[k] * dh_15[k]
                  + f_2 * fh_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_x, pb_z, dg_7, dg_18, dh_16, dh_18, fh_s_26, \
                         fh_s_27, fh_s_28, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_x[k] * dh_16[k]
                  + f_2 * fh_s_26[k];

        t_27[k] = f_9 * dg_18[k]
                  + pa_x[k] * dh_18[k]
                  + f_2 * fh_s_27[k];

        t_28[k] = f_7 * dg_7[k]
                  + f_2 * fh_s_28[k]
                  + pb_z[k] * fg_16[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_x, pb_x, dg_21, dg_23, dg_27, dh_20, dh_22, \
                         fh_s_29, fh_s_30, fh_s_31, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * dg_21[k]
                  + pa_x[k] * dh_20[k]
                  + f_2 * fh_s_29[k];

        t_30[k] = f_7 * dg_23[k]
                  + pa_x[k] * dh_22[k]
                  + f_2 * fh_s_30[k];

        t_31[k] = f_5 * dg_27[k]
                  + f_2 * fh_s_31[k]
                  + pb_x[k] * fg_19[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_x, pb_x, dh_27, ff_s_13, ff_s_14, fh_s_32, \
                         fh_s_33, fh_s_34, ff_13, ff_14, fg_20, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pa_x[k] * dh_27[k]
                  + f_2 * fh_s_32[k];

        t_33[k] = -f_1 * ff_s_13[k]
                  + f_2 * fh_s_33[k]
                  + f_3 * ff_13[k]
                  + pb_x[k] * fg_20[k];

        t_34[k] = -f_6 * ff_s_14[k]
                  + f_2 * fh_s_34[k]
                  + f_7 * ff_14[k]
                  + pb_x[k] * fg_21[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_x, pb_z, ff_s_15, ff_s_16, fh_s_35, fh_s_36, \
                         fh_s_37, ff_15, ff_16, fg_21, fg_22, fg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -f_6 * ff_s_15[k]
                  + f_2 * fh_s_35[k]
                  + f_7 * ff_15[k]
                  + pb_x[k] * fg_22[k];

        t_36[k] = -f_4 * ff_s_16[k]
                  + f_2 * fh_s_36[k]
                  + f_5 * ff_16[k]
                  + pb_x[k] * fg_23[k];

        t_37[k] = f_2 * fh_s_37[k]
                  + pb_z[k] * fg_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pb_x, ff_s_18, fh_s_38, fh_s_39, fh_s_40, ff_18, \
                         fg_24, fg_25, fg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -f_4 * ff_s_18[k]
                  + f_2 * fh_s_38[k]
                  + f_5 * ff_18[k]
                  + pb_x[k] * fg_24[k];

        t_39[k] = f_2 * fh_s_39[k]
                  + pb_x[k] * fg_25[k];

        t_40[k] = f_2 * fh_s_40[k]
                  + pb_x[k] * fg_27[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pb_y, pb_z, dg_12, ff_s_16, fh_s_41, fh_s_42, \
                         fh_s_43, ff_16, fg_25, fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_0 * dg_12[k]
                  - f_1 * ff_s_16[k]
                  + f_2 * fh_s_41[k]
                  + f_3 * ff_16[k]
                  + pb_y[k] * fg_25[k];

        t_42[k] = f_2 * fh_s_42[k]
                  + pb_z[k] * fg_25[k];

        t_43[k] = -f_4 * ff_s_16[k]
                  + f_2 * fh_s_43[k]
                  + f_5 * ff_16[k]
                  + pb_z[k] * fg_26[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pb_y, pb_z, dg_15, ff_s_17, ff_s_18, fh_s_44, \
                         fh_s_45, fh_s_46, ff_17, ff_18, fg_27, fg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = -f_6 * ff_s_17[k]
                  + f_2 * fh_s_44[k]
                  + f_7 * ff_17[k]
                  + pb_z[k] * fg_27[k];

        t_45[k] = f_0 * dg_15[k]
                  + f_2 * fh_s_45[k]
                  + pb_y[k] * fg_28[k];

        t_46[k] = -f_1 * ff_s_18[k]
                  + f_2 * fh_s_46[k]
                  + f_3 * ff_18[k]
                  + pb_z[k] * fg_28[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_z, pb_x, dh_11, ff_s_19, ff_s_21, fh_s_47, \
                         fh_s_48, fh_s_49, ff_19, ff_21, fg_29, fg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -f_6 * ff_s_19[k]
                  + f_2 * fh_s_47[k]
                  + f_7 * ff_19[k]
                  + pb_x[k] * fg_29[k];

        t_48[k] = -f_4 * ff_s_21[k]
                  + f_2 * fh_s_48[k]
                  + f_5 * ff_21[k]
                  + pb_x[k] * fg_30[k];

        t_49[k] = pa_z[k] * dh_11[k]
                  + f_2 * fh_s_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_z, pb_z, dg_12, dg_13, dg_14, dh_13, dh_14, \
                         fh_s_50, fh_s_51, fh_s_52, fg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_5 * dg_12[k]
                  + f_2 * fh_s_50[k]
                  + pb_z[k] * fg_31[k];

        t_51[k] = f_7 * dg_13[k]
                  + pa_z[k] * dh_13[k]
                  + f_2 * fh_s_51[k];

        t_52[k] = f_0 * dg_14[k]
                  + pa_z[k] * dh_14[k]
                  + f_2 * fh_s_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_y, pb_y, ph_s_5, ph_5, dg_17, dg_19, dh_17, \
                         dh_19, fh_s_53, fh_s_54, fh_s_55, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_7 * dg_17[k]
                  + f_2 * fh_s_53[k]
                  + pb_y[k] * fg_32[k];

        t_54[k] = -f_8 * ph_s_5[k]
                  + f_5 * ph_5[k]
                  + pa_y[k] * dh_17[k]
                  + f_2 * fh_s_54[k];

        t_55[k] = f_7 * dg_19[k]
                  + pa_y[k] * dh_19[k]
                  + f_2 * fh_s_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pa_y, pb_z, dg_16, dg_20, dg_24, dh_21, dh_23, \
                         fh_s_56, fh_s_57, fh_s_58, fg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * dg_20[k]
                  + pa_y[k] * dh_21[k]
                  + f_2 * fh_s_56[k];

        t_57[k] = f_9 * dg_24[k]
                  + pa_y[k] * dh_23[k]
                  + f_2 * fh_s_57[k];

        t_58[k] = f_7 * dg_16[k]
                  + f_2 * fh_s_58[k]
                  + pb_z[k] * fg_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pa_y, pb_y, dg_25, dg_26, dg_27, dh_24, dh_25, \
                         fh_s_59, fh_s_60, fh_s_61, fg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_0 * dg_25[k]
                  + pa_y[k] * dh_24[k]
                  + f_2 * fh_s_59[k];

        t_60[k] = f_7 * dg_26[k]
                  + pa_y[k] * dh_25[k]
                  + f_2 * fh_s_60[k];

        t_61[k] = f_5 * dg_27[k]
                  + f_2 * fh_s_61[k]
                  + pb_y[k] * fg_38[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pa_y, pb_x, pb_y, dh_27, ff_s_26, fh_s_62, fh_s_63, \
                         fh_s_64, ff_25, fg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_y[k] * dh_27[k]
                  + f_2 * fh_s_62[k];

        t_63[k] = -f_1 * ff_s_26[k]
                  + f_2 * fh_s_63[k]
                  + f_3 * ff_25[k]
                  + pb_x[k] * fg_39[k];

        t_64[k] = f_2 * fh_s_64[k]
                  + pb_y[k] * fg_39[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pb_x, ff_s_28, ff_s_29, ff_s_30, fh_s_65, fh_s_66, \
                         fh_s_67, ff_27, ff_28, ff_29, fg_41, fg_42, \
                         fg_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -f_6 * ff_s_28[k]
                  + f_2 * fh_s_65[k]
                  + f_7 * ff_27[k]
                  + pb_x[k] * fg_41[k];

        t_66[k] = -f_6 * ff_s_29[k]
                  + f_2 * fh_s_66[k]
                  + f_7 * ff_28[k]
                  + pb_x[k] * fg_42[k];

        t_67[k] = -f_4 * ff_s_30[k]
                  + f_2 * fh_s_67[k]
                  + f_5 * ff_29[k]
                  + pb_x[k] * fg_43[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, ff_s_33, fh_s_68, fh_s_69, \
                         fh_s_70, fh_s_71, ff_32, fg_42, fg_44, fg_45, \
                         fg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_2 * fh_s_68[k]
                  + pb_y[k] * fg_42[k];

        t_69[k] = -f_4 * ff_s_33[k]
                  + f_2 * fh_s_69[k]
                  + f_5 * ff_32[k]
                  + pb_x[k] * fg_44[k];

        t_70[k] = f_2 * fh_s_70[k]
                  + pb_x[k] * fg_45[k];

        t_71[k] = f_2 * fh_s_71[k]
                  + pb_x[k] * fg_47[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_x, pb_y, ff_s_30, ff_s_31, fh_s_72, fh_s_73, \
                         fh_s_74, ff_29, ff_30, fg_45, fg_46, fg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_2 * fh_s_72[k]
                  + pb_x[k] * fg_49[k];

        t_73[k] = -f_1 * ff_s_30[k]
                  + f_2 * fh_s_73[k]
                  + f_3 * ff_29[k]
                  + pb_y[k] * fg_45[k];

        t_74[k] = -f_10 * ff_s_31[k]
                  + f_2 * fh_s_74[k]
                  + f_0 * ff_30[k]
                  + pb_y[k] * fg_46[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_y, ff_s_32, ff_s_33, fh_s_75, fh_s_76, fh_s_77, \
                         ff_31, ff_32, fg_47, fg_48, fg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -f_6 * ff_s_32[k]
                  + f_2 * fh_s_75[k]
                  + f_7 * ff_31[k]
                  + pb_y[k] * fg_47[k];

        t_76[k] = -f_4 * ff_s_33[k]
                  + f_2 * fh_s_76[k]
                  + f_5 * ff_32[k]
                  + pb_y[k] * fg_48[k];

        t_77[k] = f_2 * fh_s_77[k]
                  + pb_y[k] * fg_49[k];
    }

#pragma omp simd aligned(t_78, pb_z, dg_27, ff_s_33, fh_s_78, ff_32, \
                         fg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * dg_27[k]
                  - f_1 * ff_s_33[k]
                  + f_2 * fh_s_78[k]
                  + f_3 * ff_32[k]
                  + pb_z[k] * fg_49[k];
    }
}

auto
compute_prim_fh_kinetic_energy_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t ph_s, const size_t ph,
                                 const size_t dg, const size_t dh, const size_t ff_s,
                                 const size_t fh_s, const size_t ff, const size_t fg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = beta / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_s_0 = buffer.data(ph_s + 0);
    const auto *ph_s_3 = buffer.data(ph_s + 3);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_3 = buffer.data(ph + 3);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);

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

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_30 = buffer.data(ff_s + 30);

    const auto *fh_s_0 = buffer.data(fh_s + 0);
    const auto *fh_s_1 = buffer.data(fh_s + 1);
    const auto *fh_s_2 = buffer.data(fh_s + 2);
    const auto *fh_s_3 = buffer.data(fh_s + 3);
    const auto *fh_s_4 = buffer.data(fh_s + 4);
    const auto *fh_s_5 = buffer.data(fh_s + 5);
    const auto *fh_s_6 = buffer.data(fh_s + 6);
    const auto *fh_s_7 = buffer.data(fh_s + 7);
    const auto *fh_s_8 = buffer.data(fh_s + 8);
    const auto *fh_s_9 = buffer.data(fh_s + 9);
    const auto *fh_s_10 = buffer.data(fh_s + 10);
    const auto *fh_s_11 = buffer.data(fh_s + 11);
    const auto *fh_s_12 = buffer.data(fh_s + 12);
    const auto *fh_s_13 = buffer.data(fh_s + 13);
    const auto *fh_s_14 = buffer.data(fh_s + 14);
    const auto *fh_s_15 = buffer.data(fh_s + 15);
    const auto *fh_s_16 = buffer.data(fh_s + 16);
    const auto *fh_s_17 = buffer.data(fh_s + 17);
    const auto *fh_s_18 = buffer.data(fh_s + 18);
    const auto *fh_s_19 = buffer.data(fh_s + 19);
    const auto *fh_s_20 = buffer.data(fh_s + 20);
    const auto *fh_s_21 = buffer.data(fh_s + 21);
    const auto *fh_s_22 = buffer.data(fh_s + 22);
    const auto *fh_s_23 = buffer.data(fh_s + 23);
    const auto *fh_s_24 = buffer.data(fh_s + 24);
    const auto *fh_s_25 = buffer.data(fh_s + 25);
    const auto *fh_s_26 = buffer.data(fh_s + 26);
    const auto *fh_s_27 = buffer.data(fh_s + 27);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_z, dg_0, ff_s_0, ff_s_2, fh_s_0, fh_s_1, \
                         fh_s_2, ff_0, ff_2, fg_0, fg_1, fg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 - f_1 * ff_s_0[k]
                 + f_2 * fh_s_0[k]
                 + f_3 * ff_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_1[k]
                 + f_5 * ff_0[k]
                 + pb_z[k] * fg_1[k];

        t_2[k] = -f_6 * ff_s_2[k]
                 + f_2 * fh_s_2[k]
                 + f_7 * ff_2[k]
                 + pb_z[k] * fg_3[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pa_z, ph_s_0, ph_0, dh_0, dh_3, fh_s_3, \
                         fh_s_4, fh_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * dh_0[k]
                 + f_2 * fh_s_3[k];

        t_4[k] = -f_8 * ph_s_0[k]
                 + f_5 * ph_0[k]
                 + pa_x[k] * dh_3[k]
                 + f_2 * fh_s_4[k];

        t_5[k] = pa_z[k] * dh_0[k]
                 + f_2 * fh_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, ph_s_3, ph_3, dg_1, dg_2, dh_1, dh_2, \
                         dh_4, fh_s_6, fh_s_7, fh_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_7 * dg_1[k]
                 + pa_z[k] * dh_1[k]
                 + f_2 * fh_s_6[k];

        t_7[k] = f_0 * dg_2[k]
                 + pa_z[k] * dh_2[k]
                 + f_2 * fh_s_7[k];

        t_8[k] = -f_8 * ph_s_3[k]
                 + f_5 * ph_3[k]
                 + pa_x[k] * dh_4[k]
                 + f_2 * fh_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, dh_5, dh_8, dh_9, dh_14, fh_s_9, \
                         fh_s_10, fh_s_11, fh_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_x[k] * dh_5[k]
                 + f_2 * fh_s_9[k];

        t_10[k] = pa_x[k] * dh_8[k]
                  + f_2 * fh_s_10[k];

        t_11[k] = pa_x[k] * dh_9[k]
                  + f_2 * fh_s_11[k];

        t_12[k] = pa_x[k] * dh_14[k]
                  + f_2 * fh_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_y, pb_z, dg_6, ff_s_16, ff_s_17, fh_s_13, \
                         fh_s_14, fh_s_15, ff_16, ff_17, fg_18, fg_19, \
                         fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dg_6[k]
                  - f_1 * ff_s_16[k]
                  + f_2 * fh_s_13[k]
                  + f_3 * ff_16[k]
                  + pb_y[k] * fg_18[k];

        t_14[k] = -f_4 * ff_s_16[k]
                  + f_2 * fh_s_14[k]
                  + f_5 * ff_16[k]
                  + pb_z[k] * fg_19[k];

        t_15[k] = -f_6 * ff_s_17[k]
                  + f_2 * fh_s_15[k]
                  + f_7 * ff_17[k]
                  + pb_z[k] * fg_20[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_z, dg_7, dg_8, dh_5, dh_6, dh_7, fh_s_16, \
                         fh_s_17, fh_s_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_z[k] * dh_5[k]
                  + f_2 * fh_s_16[k];

        t_17[k] = f_7 * dg_7[k]
                  + pa_z[k] * dh_6[k]
                  + f_2 * fh_s_17[k];

        t_18[k] = f_0 * dg_8[k]
                  + pa_z[k] * dh_7[k]
                  + f_2 * fh_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, ph_s_3, ph_3, dg_12, dg_13, dh_10, dh_11, \
                         dh_12, fh_s_19, fh_s_20, fh_s_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_8 * ph_s_3[k]
                  + f_5 * ph_3[k]
                  + pa_y[k] * dh_10[k]
                  + f_2 * fh_s_19[k];

        t_20[k] = f_9 * dg_12[k]
                  + pa_y[k] * dh_11[k]
                  + f_2 * fh_s_20[k];

        t_21[k] = f_0 * dg_13[k]
                  + pa_y[k] * dh_12[k]
                  + f_2 * fh_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pb_y, dg_14, dh_13, dh_14, ff_s_28, fh_s_22, \
                         fh_s_23, fh_s_24, ff_28, fg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_7 * dg_14[k]
                  + pa_y[k] * dh_13[k]
                  + f_2 * fh_s_22[k];

        t_23[k] = pa_y[k] * dh_14[k]
                  + f_2 * fh_s_23[k];

        t_24[k] = -f_1 * ff_s_28[k]
                  + f_2 * fh_s_24[k]
                  + f_3 * ff_28[k]
                  + pb_y[k] * fg_31[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pb_y, pb_z, dg_15, ff_s_29, ff_s_30, fh_s_25, \
                         fh_s_26, fh_s_27, ff_29, ff_30, fg_32, fg_33, \
                         fg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -f_6 * ff_s_29[k]
                  + f_2 * fh_s_25[k]
                  + f_7 * ff_29[k]
                  + pb_y[k] * fg_32[k];

        t_26[k] = -f_4 * ff_s_30[k]
                  + f_2 * fh_s_26[k]
                  + f_5 * ff_30[k]
                  + pb_y[k] * fg_33[k];

        t_27[k] = f_0 * dg_15[k]
                  - f_1 * ff_s_30[k]
                  + f_2 * fh_s_27[k]
                  + f_3 * ff_30[k]
                  + pb_z[k] * fg_34[k];
    }
}

auto
compute_prim_fh_kinetic_energy_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t ph_s, const size_t ph,
                                 const size_t dg, const size_t dh, const size_t ff_s,
                                 const size_t fh_s, const size_t ff, const size_t fg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = beta / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_s_1 = buffer.data(ph_s + 1);
    const auto *ph_s_5 = buffer.data(ph_s + 5);

    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_5 = buffer.data(ph + 5);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
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
    const auto *dh_21 = buffer.data(dh + 21);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_15 = buffer.data(ff_s + 15);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);

    const auto *fh_s_0 = buffer.data(fh_s + 0);
    const auto *fh_s_1 = buffer.data(fh_s + 1);
    const auto *fh_s_2 = buffer.data(fh_s + 2);
    const auto *fh_s_3 = buffer.data(fh_s + 3);
    const auto *fh_s_4 = buffer.data(fh_s + 4);
    const auto *fh_s_5 = buffer.data(fh_s + 5);
    const auto *fh_s_6 = buffer.data(fh_s + 6);
    const auto *fh_s_7 = buffer.data(fh_s + 7);
    const auto *fh_s_8 = buffer.data(fh_s + 8);
    const auto *fh_s_9 = buffer.data(fh_s + 9);
    const auto *fh_s_10 = buffer.data(fh_s + 10);
    const auto *fh_s_11 = buffer.data(fh_s + 11);
    const auto *fh_s_12 = buffer.data(fh_s + 12);
    const auto *fh_s_13 = buffer.data(fh_s + 13);
    const auto *fh_s_14 = buffer.data(fh_s + 14);
    const auto *fh_s_15 = buffer.data(fh_s + 15);
    const auto *fh_s_16 = buffer.data(fh_s + 16);
    const auto *fh_s_17 = buffer.data(fh_s + 17);
    const auto *fh_s_18 = buffer.data(fh_s + 18);
    const auto *fh_s_19 = buffer.data(fh_s + 19);
    const auto *fh_s_20 = buffer.data(fh_s + 20);
    const auto *fh_s_21 = buffer.data(fh_s + 21);
    const auto *fh_s_22 = buffer.data(fh_s + 22);
    const auto *fh_s_23 = buffer.data(fh_s + 23);
    const auto *fh_s_24 = buffer.data(fh_s + 24);
    const auto *fh_s_25 = buffer.data(fh_s + 25);
    const auto *fh_s_26 = buffer.data(fh_s + 26);
    const auto *fh_s_27 = buffer.data(fh_s + 27);
    const auto *fh_s_28 = buffer.data(fh_s + 28);
    const auto *fh_s_29 = buffer.data(fh_s + 29);
    const auto *fh_s_30 = buffer.data(fh_s + 30);
    const auto *fh_s_31 = buffer.data(fh_s + 31);
    const auto *fh_s_32 = buffer.data(fh_s + 32);
    const auto *fh_s_33 = buffer.data(fh_s + 33);
    const auto *fh_s_34 = buffer.data(fh_s + 34);
    const auto *fh_s_35 = buffer.data(fh_s + 35);
    const auto *fh_s_36 = buffer.data(fh_s + 36);
    const auto *fh_s_37 = buffer.data(fh_s + 37);
    const auto *fh_s_38 = buffer.data(fh_s + 38);
    const auto *fh_s_39 = buffer.data(fh_s + 39);
    const auto *fh_s_40 = buffer.data(fh_s + 40);
    const auto *fh_s_41 = buffer.data(fh_s + 41);
    const auto *fh_s_42 = buffer.data(fh_s + 42);
    const auto *fh_s_43 = buffer.data(fh_s + 43);
    const auto *fh_s_44 = buffer.data(fh_s + 44);
    const auto *fh_s_45 = buffer.data(fh_s + 45);
    const auto *fh_s_46 = buffer.data(fh_s + 46);
    const auto *fh_s_47 = buffer.data(fh_s + 47);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dg_0, ff_s_0, fh_s_0, fh_s_1, \
                         fh_s_2, ff_0, fg_0, fg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 - f_1 * ff_s_0[k]
                 + f_2 * fh_s_0[k]
                 + f_3 * ff_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = f_2 * fh_s_1[k]
                 + pb_z[k] * fg_0[k];

        t_2[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_2[k]
                 + f_5 * ff_0[k]
                 + pb_y[k] * fg_1[k];
    }

#pragma omp simd aligned(t_3, t_4, pb_y, pb_z, ff_s_0, ff_s_1, fh_s_3, fh_s_4, ff_0, ff_1, \
                         fg_2, fg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_3[k]
                 + f_5 * ff_0[k]
                 + pb_z[k] * fg_2[k];

        t_4[k] = -f_6 * ff_s_1[k]
                 + f_2 * fh_s_4[k]
                 + f_7 * ff_1[k]
                 + pb_y[k] * fg_3[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_x, pa_y, pb_z, ph_s_1, ph_1, dh_0, dh_4, ff_s_2, \
                         fh_s_5, fh_s_6, fh_s_7, ff_2, fg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -f_6 * ff_s_2[k]
                 + f_2 * fh_s_5[k]
                 + f_7 * ff_2[k]
                 + pb_z[k] * fg_4[k];

        t_6[k] = pa_y[k] * dh_0[k]
                 + f_2 * fh_s_6[k];

        t_7[k] = -f_8 * ph_s_1[k]
                 + f_5 * ph_1[k]
                 + pa_x[k] * dh_4[k]
                 + f_2 * fh_s_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_z, dg_1, dg_2, dh_0, dh_2, dh_3, fh_s_8, fh_s_9, \
                         fh_s_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = pa_z[k] * dh_0[k]
                 + f_2 * fh_s_8[k];

        t_9[k] = f_7 * dg_1[k]
                 + pa_z[k] * dh_2[k]
                 + f_2 * fh_s_9[k];

        t_10[k] = f_0 * dg_2[k]
                  + pa_z[k] * dh_3[k]
                  + f_2 * fh_s_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, ph_s_5, ph_5, dg_6, dg_7, dh_5, dh_6, dh_7, \
                         fh_s_11, fh_s_12, fh_s_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = -f_8 * ph_s_5[k]
                  + f_5 * ph_5[k]
                  + pa_x[k] * dh_5[k]
                  + f_2 * fh_s_11[k];

        t_12[k] = f_0 * dg_6[k]
                  + pa_x[k] * dh_6[k]
                  + f_2 * fh_s_12[k];

        t_13[k] = f_7 * dg_7[k]
                  + pa_x[k] * dh_7[k]
                  + f_2 * fh_s_13[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pb_x, dg_8, dh_8, dh_12, dh_13, \
                         fh_s_14, fh_s_15, fh_s_16, fh_s_17, fg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * dg_8[k]
                  + f_2 * fh_s_14[k]
                  + pb_x[k] * fg_13[k];

        t_15[k] = pa_x[k] * dh_8[k]
                  + f_2 * fh_s_15[k];

        t_16[k] = pa_x[k] * dh_12[k]
                  + f_2 * fh_s_16[k];

        t_17[k] = pa_x[k] * dh_13[k]
                  + f_2 * fh_s_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pb_x, dg_13, dg_14, dg_18, dh_15, dh_16, \
                         fh_s_18, fh_s_19, fh_s_20, fg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * dg_13[k]
                  + pa_x[k] * dh_15[k]
                  + f_2 * fh_s_18[k];

        t_19[k] = f_7 * dg_14[k]
                  + pa_x[k] * dh_16[k]
                  + f_2 * fh_s_19[k];

        t_20[k] = f_5 * dg_18[k]
                  + f_2 * fh_s_20[k]
                  + pb_x[k] * fg_17[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pb_x, dh_21, ff_s_12, ff_s_13, fh_s_21, \
                         fh_s_22, fh_s_23, ff_12, ff_13, fg_18, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_x[k] * dh_21[k]
                  + f_2 * fh_s_21[k];

        t_22[k] = -f_1 * ff_s_12[k]
                  + f_2 * fh_s_22[k]
                  + f_3 * ff_12[k]
                  + pb_x[k] * fg_18[k];

        t_23[k] = -f_6 * ff_s_13[k]
                  + f_2 * fh_s_23[k]
                  + f_7 * ff_13[k]
                  + pb_x[k] * fg_19[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_x, pb_y, pb_z, dg_8, ff_s_14, fh_s_24, fh_s_25, \
                         fh_s_26, ff_14, fg_20, fg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_4 * ff_s_14[k]
                  + f_2 * fh_s_24[k]
                  + f_5 * ff_14[k]
                  + pb_x[k] * fg_20[k];

        t_25[k] = f_0 * dg_8[k]
                  - f_1 * ff_s_14[k]
                  + f_2 * fh_s_25[k]
                  + f_3 * ff_14[k]
                  + pb_y[k] * fg_21[k];

        t_26[k] = f_2 * fh_s_26[k]
                  + pb_z[k] * fg_21[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_z, pb_z, dh_8, ff_s_14, ff_s_15, fh_s_27, \
                         fh_s_28, fh_s_29, ff_14, ff_15, fg_22, fg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_4 * ff_s_14[k]
                  + f_2 * fh_s_27[k]
                  + f_5 * ff_14[k]
                  + pb_z[k] * fg_22[k];

        t_28[k] = -f_6 * ff_s_15[k]
                  + f_2 * fh_s_28[k]
                  + f_7 * ff_15[k]
                  + pb_z[k] * fg_23[k];

        t_29[k] = pa_z[k] * dh_8[k]
                  + f_2 * fh_s_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pa_z, ph_s_5, ph_5, dg_9, dg_10, dh_10, \
                         dh_11, dh_14, fh_s_30, fh_s_31, fh_s_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * dg_9[k]
                  + pa_z[k] * dh_10[k]
                  + f_2 * fh_s_30[k];

        t_31[k] = f_0 * dg_10[k]
                  + pa_z[k] * dh_11[k]
                  + f_2 * fh_s_31[k];

        t_32[k] = -f_8 * ph_s_5[k]
                  + f_5 * ph_5[k]
                  + pa_y[k] * dh_14[k]
                  + f_2 * fh_s_32[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, dg_15, dg_16, dg_17, dh_17, dh_18, dh_19, \
                         fh_s_33, fh_s_34, fh_s_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_9 * dg_15[k]
                  + pa_y[k] * dh_17[k]
                  + f_2 * fh_s_33[k];

        t_34[k] = f_0 * dg_16[k]
                  + pa_y[k] * dh_18[k]
                  + f_2 * fh_s_34[k];

        t_35[k] = f_7 * dg_17[k]
                  + pa_y[k] * dh_19[k]
                  + f_2 * fh_s_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_y, pb_x, pb_y, dg_18, dh_21, ff_s_21, fh_s_36, \
                         fh_s_37, fh_s_38, ff_21, fg_30, fg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_5 * dg_18[k]
                  + f_2 * fh_s_36[k]
                  + pb_y[k] * fg_30[k];

        t_37[k] = pa_y[k] * dh_21[k]
                  + f_2 * fh_s_37[k];

        t_38[k] = -f_1 * ff_s_21[k]
                  + f_2 * fh_s_38[k]
                  + f_3 * ff_21[k]
                  + pb_x[k] * fg_31[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_x, ff_s_23, ff_s_24, ff_s_25, fh_s_39, fh_s_40, \
                         fh_s_41, ff_23, ff_24, ff_25, fg_33, fg_34, \
                         fg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -f_6 * ff_s_23[k]
                  + f_2 * fh_s_39[k]
                  + f_7 * ff_23[k]
                  + pb_x[k] * fg_33[k];

        t_40[k] = -f_6 * ff_s_24[k]
                  + f_2 * fh_s_40[k]
                  + f_7 * ff_24[k]
                  + pb_x[k] * fg_34[k];

        t_41[k] = -f_4 * ff_s_25[k]
                  + f_2 * fh_s_41[k]
                  + f_5 * ff_25[k]
                  + pb_x[k] * fg_35[k];
    }

#pragma omp simd aligned(t_42, t_43, pb_x, pb_y, ff_s_25, ff_s_27, fh_s_42, fh_s_43, ff_25, \
                         ff_27, fg_36, fg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = -f_4 * ff_s_27[k]
                  + f_2 * fh_s_42[k]
                  + f_5 * ff_27[k]
                  + pb_x[k] * fg_36[k];

        t_43[k] = -f_1 * ff_s_25[k]
                  + f_2 * fh_s_43[k]
                  + f_3 * ff_25[k]
                  + pb_y[k] * fg_37[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pb_y, ff_s_26, ff_s_27, fh_s_44, fh_s_45, fh_s_46, \
                         ff_26, ff_27, fg_38, fg_39, fg_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = -f_6 * ff_s_26[k]
                  + f_2 * fh_s_44[k]
                  + f_7 * ff_26[k]
                  + pb_y[k] * fg_38[k];

        t_45[k] = -f_4 * ff_s_27[k]
                  + f_2 * fh_s_45[k]
                  + f_5 * ff_27[k]
                  + pb_y[k] * fg_39[k];

        t_46[k] = f_2 * fh_s_46[k]
                  + pb_y[k] * fg_40[k];
    }

#pragma omp simd aligned(t_47, pb_z, dg_18, ff_s_27, fh_s_47, ff_27, \
                         fg_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * dg_18[k]
                  - f_1 * ff_s_27[k]
                  + f_2 * fh_s_47[k]
                  + f_3 * ff_27[k]
                  + pb_z[k] * fg_40[k];
    }
}

}  // namespace simdkin
