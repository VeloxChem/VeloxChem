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


#include "SimdTransformGG.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_gg(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t gg,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 13.125 * std::sqrt(2.0);
    const auto f_1 = 4.375 * std::sqrt(2.0);
    const auto f_2 = 1.25 * std::sqrt(7.0);
    const auto f_3 = 7.5 * std::sqrt(7.0);
    const auto f_4 = 1.875 * std::sqrt(14.0);
    const auto f_5 = 2.5 * std::sqrt(14.0);
    const auto f_6 = 0.1875 * std::sqrt(35.0);
    const auto f_7 = 0.375 * std::sqrt(35.0);
    const auto f_8 = 1.5 * std::sqrt(35.0);
    const auto f_9 = 0.5 * std::sqrt(35.0);
    const auto f_10 = 0.625 * std::sqrt(7.0);
    const auto f_11 = 3.75 * std::sqrt(7.0);
    const auto f_12 = 11.25 * std::sqrt(14.0);
    const auto f_13 = 0.625 * std::sqrt(14.0);
    const auto f_14 = 3.75 * std::sqrt(14.0);
    const auto f_15 = 5.625 * std::sqrt(7.0);
    const auto f_16 = 1.875 * std::sqrt(7.0);
    const auto f_17 = 2.5 * std::sqrt(7.0);
    const auto f_18 = 0.28125 * std::sqrt(70.0);
    const auto f_19 = 0.5625 * std::sqrt(70.0);
    const auto f_20 = 2.25 * std::sqrt(70.0);
    const auto f_21 = 0.75 * std::sqrt(70.0);
    const auto f_22 = 0.09375 * std::sqrt(70.0);
    const auto f_23 = 0.1875 * std::sqrt(70.0);
    const auto f_24 = 0.25 * std::sqrt(70.0);
    const auto f_25 = 0.9375 * std::sqrt(14.0);
    const auto f_26 = 5.625 * std::sqrt(14.0);
    const auto f_27 = 0.3125 * std::sqrt(14.0);
    const auto f_28 = 3.28125 * std::sqrt(2.0);
    const auto f_29 = 19.6875 * std::sqrt(2.0);
    const auto f_30 = 1.09375 * std::sqrt(2.0);
    const auto f_31 = 6.5625 * std::sqrt(2.0);
    const auto f_32 = 1.875 * std::sqrt(2.0);
    const auto f_33 = 2.5 * std::sqrt(2.0);
    const auto f_34 = 11.25 * std::sqrt(2.0);
    const auto f_35 = 15.0 * std::sqrt(2.0);
    const auto f_36 = 0.1875 * std::sqrt(5.0);
    const auto f_37 = 0.375 * std::sqrt(5.0);
    const auto f_38 = 1.5 * std::sqrt(5.0);
    const auto f_39 = 0.5 * std::sqrt(5.0);
    const auto f_40 = 1.125 * std::sqrt(5.0);
    const auto f_41 = 2.25 * std::sqrt(5.0);
    const auto f_42 = 9.0 * std::sqrt(5.0);
    const auto f_43 = 3.0 * std::sqrt(5.0);
    const auto f_44 = 0.3125 * std::sqrt(7.0);
    const auto f_45 = 11.25 * std::sqrt(7.0);
    const auto f_46 = 0.28125 * std::sqrt(10.0);
    const auto f_47 = 0.5625 * std::sqrt(10.0);
    const auto f_48 = 2.25 * std::sqrt(10.0);
    const auto f_49 = 0.75 * std::sqrt(10.0);
    const auto f_50 = 0.375 * std::sqrt(10.0);
    const auto f_51 = 3.0 * std::sqrt(10.0);
    const auto f_52 = std::sqrt(10.0);
    const auto f_53 = 0.9375 * std::sqrt(2.0);
    const auto f_54 = 5.625 * std::sqrt(2.0);
    const auto f_55 = 1.25 * std::sqrt(2.0);
    const auto f_56 = 7.5 * std::sqrt(2.0);
    const auto f_57 = 0.46875 * std::sqrt(14.0);
    const auto f_58 = 2.8125 * std::sqrt(14.0);
    const auto f_59 = 0.09375 * std::sqrt(5.0);
    const auto f_60 = 0.5625 * std::sqrt(5.0);
    const auto f_61 = 0.75 * std::sqrt(5.0);
    const auto f_62 = 4.5 * std::sqrt(5.0);
    const auto f_63 = 0.25 * std::sqrt(5.0);
    const auto f_64 = 0.046875 * std::sqrt(35.0);
    const auto f_65 = 0.28125 * std::sqrt(35.0);
    const auto f_66 = 0.09375 * std::sqrt(35.0);
    const auto f_67 = 0.5625 * std::sqrt(35.0);
    const auto f_68 = 2.25 * std::sqrt(35.0);
    const auto f_69 = 0.125 * std::sqrt(35.0);
    const auto f_70 = 0.75 * std::sqrt(35.0);
    const auto f_71 = 0.15625 * std::sqrt(7.0);
    const auto f_72 = 0.9375 * std::sqrt(7.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;
    auto *g_7 = values + 7 * nvalues;
    auto *g_8 = values + 8 * nvalues;
    auto *g_9 = values + 9 * nvalues;
    auto *g_10 = values + 10 * nvalues;
    auto *g_11 = values + 11 * nvalues;
    auto *g_12 = values + 12 * nvalues;
    auto *g_13 = values + 13 * nvalues;
    auto *g_14 = values + 14 * nvalues;
    auto *g_15 = values + 15 * nvalues;
    auto *g_16 = values + 16 * nvalues;
    auto *g_17 = values + 17 * nvalues;
    auto *g_18 = values + 18 * nvalues;
    auto *g_19 = values + 19 * nvalues;
    auto *g_20 = values + 20 * nvalues;
    auto *g_21 = values + 21 * nvalues;
    auto *g_22 = values + 22 * nvalues;
    auto *g_23 = values + 23 * nvalues;
    auto *g_24 = values + 24 * nvalues;
    auto *g_25 = values + 25 * nvalues;
    auto *g_26 = values + 26 * nvalues;
    auto *g_27 = values + 27 * nvalues;
    auto *g_28 = values + 28 * nvalues;
    auto *g_29 = values + 29 * nvalues;
    auto *g_30 = values + 30 * nvalues;
    auto *g_31 = values + 31 * nvalues;
    auto *g_32 = values + 32 * nvalues;
    auto *g_33 = values + 33 * nvalues;
    auto *g_34 = values + 34 * nvalues;
    auto *g_35 = values + 35 * nvalues;
    auto *g_36 = values + 36 * nvalues;
    auto *g_37 = values + 37 * nvalues;
    auto *g_38 = values + 38 * nvalues;
    auto *g_39 = values + 39 * nvalues;
    auto *g_40 = values + 40 * nvalues;
    auto *g_41 = values + 41 * nvalues;
    auto *g_42 = values + 42 * nvalues;
    auto *g_43 = values + 43 * nvalues;
    auto *g_44 = values + 44 * nvalues;
    auto *g_45 = values + 45 * nvalues;
    auto *g_46 = values + 46 * nvalues;
    auto *g_47 = values + 47 * nvalues;
    auto *g_48 = values + 48 * nvalues;
    auto *g_49 = values + 49 * nvalues;
    auto *g_50 = values + 50 * nvalues;
    auto *g_51 = values + 51 * nvalues;
    auto *g_52 = values + 52 * nvalues;
    auto *g_53 = values + 53 * nvalues;
    auto *g_54 = values + 54 * nvalues;
    auto *g_55 = values + 55 * nvalues;
    auto *g_56 = values + 56 * nvalues;
    auto *g_57 = values + 57 * nvalues;
    auto *g_58 = values + 58 * nvalues;
    auto *g_59 = values + 59 * nvalues;
    auto *g_60 = values + 60 * nvalues;
    auto *g_61 = values + 61 * nvalues;
    auto *g_62 = values + 62 * nvalues;
    auto *g_63 = values + 63 * nvalues;
    auto *g_64 = values + 64 * nvalues;
    auto *g_65 = values + 65 * nvalues;
    auto *g_66 = values + 66 * nvalues;
    auto *g_67 = values + 67 * nvalues;
    auto *g_68 = values + 68 * nvalues;
    auto *g_69 = values + 69 * nvalues;
    auto *g_70 = values + 70 * nvalues;
    auto *g_71 = values + 71 * nvalues;
    auto *g_72 = values + 72 * nvalues;
    auto *g_73 = values + 73 * nvalues;
    auto *g_74 = values + 74 * nvalues;
    auto *g_75 = values + 75 * nvalues;
    auto *g_76 = values + 76 * nvalues;
    auto *g_77 = values + 77 * nvalues;
    auto *g_78 = values + 78 * nvalues;
    auto *g_79 = values + 79 * nvalues;
    auto *g_80 = values + 80 * nvalues;

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_82 = buffer.data(gg + 82);
    const auto *gg_83 = buffer.data(gg + 83);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_88 = buffer.data(gg + 88);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_91 = buffer.data(gg + 91);
    const auto *gg_92 = buffer.data(gg + 92);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_94 = buffer.data(gg + 94);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_96 = buffer.data(gg + 96);
    const auto *gg_97 = buffer.data(gg + 97);
    const auto *gg_98 = buffer.data(gg + 98);
    const auto *gg_99 = buffer.data(gg + 99);
    const auto *gg_100 = buffer.data(gg + 100);
    const auto *gg_101 = buffer.data(gg + 101);
    const auto *gg_102 = buffer.data(gg + 102);
    const auto *gg_103 = buffer.data(gg + 103);
    const auto *gg_104 = buffer.data(gg + 104);
    const auto *gg_105 = buffer.data(gg + 105);
    const auto *gg_106 = buffer.data(gg + 106);
    const auto *gg_107 = buffer.data(gg + 107);
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_109 = buffer.data(gg + 109);
    const auto *gg_110 = buffer.data(gg + 110);
    const auto *gg_111 = buffer.data(gg + 111);
    const auto *gg_112 = buffer.data(gg + 112);
    const auto *gg_113 = buffer.data(gg + 113);
    const auto *gg_114 = buffer.data(gg + 114);
    const auto *gg_115 = buffer.data(gg + 115);
    const auto *gg_116 = buffer.data(gg + 116);
    const auto *gg_117 = buffer.data(gg + 117);
    const auto *gg_118 = buffer.data(gg + 118);
    const auto *gg_119 = buffer.data(gg + 119);
    const auto *gg_120 = buffer.data(gg + 120);
    const auto *gg_121 = buffer.data(gg + 121);
    const auto *gg_122 = buffer.data(gg + 122);
    const auto *gg_123 = buffer.data(gg + 123);
    const auto *gg_124 = buffer.data(gg + 124);
    const auto *gg_125 = buffer.data(gg + 125);
    const auto *gg_126 = buffer.data(gg + 126);
    const auto *gg_127 = buffer.data(gg + 127);
    const auto *gg_128 = buffer.data(gg + 128);
    const auto *gg_129 = buffer.data(gg + 129);
    const auto *gg_130 = buffer.data(gg + 130);
    const auto *gg_131 = buffer.data(gg + 131);
    const auto *gg_132 = buffer.data(gg + 132);
    const auto *gg_133 = buffer.data(gg + 133);
    const auto *gg_134 = buffer.data(gg + 134);
    const auto *gg_135 = buffer.data(gg + 135);
    const auto *gg_136 = buffer.data(gg + 136);
    const auto *gg_137 = buffer.data(gg + 137);
    const auto *gg_138 = buffer.data(gg + 138);
    const auto *gg_139 = buffer.data(gg + 139);
    const auto *gg_140 = buffer.data(gg + 140);
    const auto *gg_141 = buffer.data(gg + 141);
    const auto *gg_142 = buffer.data(gg + 142);
    const auto *gg_143 = buffer.data(gg + 143);
    const auto *gg_144 = buffer.data(gg + 144);
    const auto *gg_145 = buffer.data(gg + 145);
    const auto *gg_146 = buffer.data(gg + 146);
    const auto *gg_147 = buffer.data(gg + 147);
    const auto *gg_148 = buffer.data(gg + 148);
    const auto *gg_149 = buffer.data(gg + 149);
    const auto *gg_150 = buffer.data(gg + 150);
    const auto *gg_151 = buffer.data(gg + 151);
    const auto *gg_152 = buffer.data(gg + 152);
    const auto *gg_153 = buffer.data(gg + 153);
    const auto *gg_154 = buffer.data(gg + 154);
    const auto *gg_155 = buffer.data(gg + 155);
    const auto *gg_156 = buffer.data(gg + 156);
    const auto *gg_157 = buffer.data(gg + 157);
    const auto *gg_158 = buffer.data(gg + 158);
    const auto *gg_159 = buffer.data(gg + 159);
    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_161 = buffer.data(gg + 161);
    const auto *gg_162 = buffer.data(gg + 162);
    const auto *gg_163 = buffer.data(gg + 163);
    const auto *gg_164 = buffer.data(gg + 164);
    const auto *gg_165 = buffer.data(gg + 165);
    const auto *gg_166 = buffer.data(gg + 166);
    const auto *gg_167 = buffer.data(gg + 167);
    const auto *gg_168 = buffer.data(gg + 168);
    const auto *gg_169 = buffer.data(gg + 169);
    const auto *gg_170 = buffer.data(gg + 170);
    const auto *gg_171 = buffer.data(gg + 171);
    const auto *gg_172 = buffer.data(gg + 172);
    const auto *gg_173 = buffer.data(gg + 173);
    const auto *gg_174 = buffer.data(gg + 174);
    const auto *gg_175 = buffer.data(gg + 175);
    const auto *gg_176 = buffer.data(gg + 176);
    const auto *gg_177 = buffer.data(gg + 177);
    const auto *gg_178 = buffer.data(gg + 178);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_180 = buffer.data(gg + 180);
    const auto *gg_181 = buffer.data(gg + 181);
    const auto *gg_182 = buffer.data(gg + 182);
    const auto *gg_183 = buffer.data(gg + 183);
    const auto *gg_184 = buffer.data(gg + 184);
    const auto *gg_185 = buffer.data(gg + 185);
    const auto *gg_186 = buffer.data(gg + 186);
    const auto *gg_187 = buffer.data(gg + 187);
    const auto *gg_188 = buffer.data(gg + 188);
    const auto *gg_189 = buffer.data(gg + 189);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_191 = buffer.data(gg + 191);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_193 = buffer.data(gg + 193);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_195 = buffer.data(gg + 195);
    const auto *gg_196 = buffer.data(gg + 196);
    const auto *gg_197 = buffer.data(gg + 197);
    const auto *gg_198 = buffer.data(gg + 198);
    const auto *gg_199 = buffer.data(gg + 199);
    const auto *gg_200 = buffer.data(gg + 200);
    const auto *gg_201 = buffer.data(gg + 201);
    const auto *gg_202 = buffer.data(gg + 202);
    const auto *gg_203 = buffer.data(gg + 203);
    const auto *gg_204 = buffer.data(gg + 204);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_206 = buffer.data(gg + 206);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_208 = buffer.data(gg + 208);
    const auto *gg_209 = buffer.data(gg + 209);
    const auto *gg_210 = buffer.data(gg + 210);
    const auto *gg_211 = buffer.data(gg + 211);
    const auto *gg_212 = buffer.data(gg + 212);
    const auto *gg_213 = buffer.data(gg + 213);
    const auto *gg_214 = buffer.data(gg + 214);
    const auto *gg_215 = buffer.data(gg + 215);
    const auto *gg_216 = buffer.data(gg + 216);
    const auto *gg_217 = buffer.data(gg + 217);
    const auto *gg_218 = buffer.data(gg + 218);
    const auto *gg_219 = buffer.data(gg + 219);
    const auto *gg_220 = buffer.data(gg + 220);
    const auto *gg_221 = buffer.data(gg + 221);
    const auto *gg_222 = buffer.data(gg + 222);
    const auto *gg_223 = buffer.data(gg + 223);
    const auto *gg_224 = buffer.data(gg + 224);

#pragma omp simd aligned(gg_16, gg_19, gg_21, gg_23, gg_26, gg_28, gg_91, gg_94, gg_96, gg_98, \
                         gg_101, gg_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = 8.75 * gg_16[k]
                 - 8.75 * gg_21[k]
                 - 8.75 * gg_91[k]
                 + 8.75 * gg_96[k];

        g_1[k] = f_0 * gg_19[k]
                 - f_1 * gg_26[k]
                 - f_0 * gg_94[k]
                 + f_1 * gg_101[k];

        g_2[k] = -f_2 * gg_16[k]
                 - f_2 * gg_21[k]
                 + f_3 * gg_23[k]
                 + f_2 * gg_91[k]
                 + f_2 * gg_96[k]
                 - f_3 * gg_98[k];

        g_3[k] = -f_4 * gg_19[k]
                 - f_4 * gg_26[k]
                 + f_5 * gg_28[k]
                 + f_4 * gg_94[k]
                 + f_4 * gg_101[k]
                 - f_5 * gg_103[k];
    }

#pragma omp simd aligned(gg_15, gg_18, gg_20, gg_25, gg_27, gg_29, gg_90, gg_93, gg_95, \
                         gg_100, gg_102, gg_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_6 * gg_15[k]
                 + f_7 * gg_18[k]
                 - f_8 * gg_20[k]
                 + f_6 * gg_25[k]
                 - f_8 * gg_27[k]
                 + f_9 * gg_29[k]
                 - f_6 * gg_90[k]
                 - f_7 * gg_93[k]
                 + f_8 * gg_95[k]
                 - f_6 * gg_100[k]
                 + f_8 * gg_102[k]
                 - f_9 * gg_104[k];
    }

#pragma omp simd aligned(gg_15, gg_17, gg_20, gg_22, gg_24, gg_25, gg_27, gg_90, gg_92, gg_95, \
                         gg_97, gg_99, gg_100, gg_102 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -f_4 * gg_17[k]
                 - f_4 * gg_22[k]
                 + f_5 * gg_24[k]
                 + f_4 * gg_92[k]
                 + f_4 * gg_97[k]
                 - f_5 * gg_99[k];

        g_6[k] = -f_10 * gg_15[k]
                 + f_11 * gg_20[k]
                 + f_10 * gg_25[k]
                 - f_11 * gg_27[k]
                 + f_10 * gg_90[k]
                 - f_11 * gg_95[k]
                 - f_10 * gg_100[k]
                 + f_11 * gg_102[k];
    }

#pragma omp simd aligned(gg_15, gg_17, gg_18, gg_22, gg_25, gg_90, gg_92, gg_93, gg_97, \
                         gg_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_1 * gg_17[k]
                 - f_0 * gg_22[k]
                 - f_1 * gg_92[k]
                 + f_0 * gg_97[k];

        g_8[k] = 2.1875 * gg_15[k]
                 - 13.125 * gg_18[k]
                 + 2.1875 * gg_25[k]
                 - 2.1875 * gg_90[k]
                 + 13.125 * gg_93[k]
                 - 2.1875 * gg_100[k];
    }

#pragma omp simd aligned(gg_61, gg_64, gg_66, gg_68, gg_71, gg_73, gg_166, gg_169, gg_171, \
                         gg_173, gg_176, gg_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_0 * gg_61[k]
                 - f_0 * gg_66[k]
                 - f_1 * gg_166[k]
                 + f_1 * gg_171[k];

        g_10[k] = 39.375 * gg_64[k]
                  - 13.125 * gg_71[k]
                  - 13.125 * gg_169[k]
                  + 4.375 * gg_176[k];

        g_11[k] = -f_4 * gg_61[k]
                  - f_4 * gg_66[k]
                  + f_12 * gg_68[k]
                  + f_13 * gg_166[k]
                  + f_13 * gg_171[k]
                  - f_14 * gg_173[k];

        g_12[k] = -f_15 * gg_64[k]
                  - f_15 * gg_71[k]
                  + f_3 * gg_73[k]
                  + f_16 * gg_169[k]
                  + f_16 * gg_176[k]
                  - f_17 * gg_178[k];
    }

#pragma omp simd aligned(gg_60, gg_63, gg_65, gg_70, gg_72, gg_74, gg_165, gg_168, gg_170, \
                         gg_175, gg_177, gg_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_18 * gg_60[k]
                  + f_19 * gg_63[k]
                  - f_20 * gg_65[k]
                  + f_18 * gg_70[k]
                  - f_20 * gg_72[k]
                  + f_21 * gg_74[k]
                  - f_22 * gg_165[k]
                  - f_23 * gg_168[k]
                  + f_21 * gg_170[k]
                  - f_22 * gg_175[k]
                  + f_21 * gg_177[k]
                  - f_24 * gg_179[k];
    }

#pragma omp simd aligned(gg_60, gg_62, gg_65, gg_67, gg_69, gg_70, gg_72, gg_165, gg_167, \
                         gg_170, gg_172, gg_174, gg_175, gg_177 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_15 * gg_62[k]
                  - f_15 * gg_67[k]
                  + f_3 * gg_69[k]
                  + f_16 * gg_167[k]
                  + f_16 * gg_172[k]
                  - f_17 * gg_174[k];

        g_15[k] = -f_25 * gg_60[k]
                  + f_26 * gg_65[k]
                  + f_25 * gg_70[k]
                  - f_26 * gg_72[k]
                  + f_27 * gg_165[k]
                  - f_4 * gg_170[k]
                  - f_27 * gg_175[k]
                  + f_4 * gg_177[k];
    }

#pragma omp simd aligned(gg_60, gg_62, gg_63, gg_67, gg_70, gg_165, gg_167, gg_168, gg_172, \
                         gg_175 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = 13.125 * gg_62[k]
                  - 39.375 * gg_67[k]
                  - 4.375 * gg_167[k]
                  + 13.125 * gg_172[k];

        g_17[k] = f_28 * gg_60[k]
                  - f_29 * gg_63[k]
                  + f_28 * gg_70[k]
                  - f_30 * gg_165[k]
                  + f_31 * gg_168[k]
                  - f_30 * gg_175[k];
    }

#pragma omp simd aligned(gg_16, gg_19, gg_21, gg_26, gg_91, gg_94, gg_96, gg_101, gg_121, \
                         gg_124, gg_126, gg_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_2 * gg_16[k]
                  + f_2 * gg_21[k]
                  - f_2 * gg_91[k]
                  + f_2 * gg_96[k]
                  + f_3 * gg_121[k]
                  - f_3 * gg_126[k];

        g_19[k] = -f_4 * gg_19[k]
                  + f_13 * gg_26[k]
                  - f_4 * gg_94[k]
                  + f_13 * gg_101[k]
                  + f_12 * gg_124[k]
                  - f_14 * gg_131[k];
    }

#pragma omp simd aligned(gg_16, gg_21, gg_23, gg_91, gg_96, gg_98, gg_121, gg_126, \
                         gg_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = 1.25 * gg_16[k]
                  + 1.25 * gg_21[k]
                  - 7.5 * gg_23[k]
                  + 1.25 * gg_91[k]
                  + 1.25 * gg_96[k]
                  - 7.5 * gg_98[k]
                  - 7.5 * gg_121[k]
                  - 7.5 * gg_126[k]
                  + 45.0 * gg_128[k];
    }

#pragma omp simd aligned(gg_19, gg_26, gg_28, gg_94, gg_101, gg_103, gg_124, gg_131, \
                         gg_133 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_32 * gg_19[k]
                  + f_32 * gg_26[k]
                  - f_33 * gg_28[k]
                  + f_32 * gg_94[k]
                  + f_32 * gg_101[k]
                  - f_33 * gg_103[k]
                  - f_34 * gg_124[k]
                  - f_34 * gg_131[k]
                  + f_35 * gg_133[k];
    }

#pragma omp simd aligned(gg_15, gg_18, gg_20, gg_25, gg_27, gg_29, gg_90, gg_93, gg_95, \
                         gg_100, gg_102, gg_104, gg_120, gg_123, gg_125, gg_130, gg_132, \
                         gg_134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_36 * gg_15[k]
                  - f_37 * gg_18[k]
                  + f_38 * gg_20[k]
                  - f_36 * gg_25[k]
                  + f_38 * gg_27[k]
                  - f_39 * gg_29[k]
                  - f_36 * gg_90[k]
                  - f_37 * gg_93[k]
                  + f_38 * gg_95[k]
                  - f_36 * gg_100[k]
                  + f_38 * gg_102[k]
                  - f_39 * gg_104[k]
                  + f_40 * gg_120[k]
                  + f_41 * gg_123[k]
                  - f_42 * gg_125[k]
                  + f_40 * gg_130[k]
                  - f_42 * gg_132[k]
                  + f_43 * gg_134[k];
    }

#pragma omp simd aligned(gg_17, gg_22, gg_24, gg_92, gg_97, gg_99, gg_122, gg_127, \
                         gg_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_32 * gg_17[k]
                  + f_32 * gg_22[k]
                  - f_33 * gg_24[k]
                  + f_32 * gg_92[k]
                  + f_32 * gg_97[k]
                  - f_33 * gg_99[k]
                  - f_34 * gg_122[k]
                  - f_34 * gg_127[k]
                  + f_35 * gg_129[k];
    }

#pragma omp simd aligned(gg_15, gg_20, gg_25, gg_27, gg_90, gg_95, gg_100, gg_102, gg_120, \
                         gg_125, gg_130, gg_132 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = 0.625 * gg_15[k]
                  - 3.75 * gg_20[k]
                  - 0.625 * gg_25[k]
                  + 3.75 * gg_27[k]
                  + 0.625 * gg_90[k]
                  - 3.75 * gg_95[k]
                  - 0.625 * gg_100[k]
                  + 3.75 * gg_102[k]
                  - 3.75 * gg_120[k]
                  + 22.5 * gg_125[k]
                  + 3.75 * gg_130[k]
                  - 22.5 * gg_132[k];
    }

#pragma omp simd aligned(gg_17, gg_22, gg_92, gg_97, gg_122, gg_127 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_13 * gg_17[k]
                  + f_4 * gg_22[k]
                  - f_13 * gg_92[k]
                  + f_4 * gg_97[k]
                  + f_14 * gg_122[k]
                  - f_12 * gg_127[k];
    }

#pragma omp simd aligned(gg_15, gg_18, gg_25, gg_90, gg_93, gg_100, gg_120, gg_123, \
                         gg_130 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_44 * gg_15[k]
                  + f_16 * gg_18[k]
                  - f_44 * gg_25[k]
                  - f_44 * gg_90[k]
                  + f_16 * gg_93[k]
                  - f_44 * gg_100[k]
                  + f_16 * gg_120[k]
                  - f_45 * gg_123[k]
                  + f_16 * gg_130[k];
    }

#pragma omp simd aligned(gg_61, gg_64, gg_66, gg_71, gg_166, gg_169, gg_171, gg_176, gg_196, \
                         gg_199, gg_201, gg_206 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_4 * gg_61[k]
                  + f_4 * gg_66[k]
                  - f_4 * gg_166[k]
                  + f_4 * gg_171[k]
                  + f_5 * gg_196[k]
                  - f_5 * gg_201[k];

        g_28[k] = -f_15 * gg_64[k]
                  + f_16 * gg_71[k]
                  - f_15 * gg_169[k]
                  + f_16 * gg_176[k]
                  + f_3 * gg_199[k]
                  - f_17 * gg_206[k];
    }

#pragma omp simd aligned(gg_61, gg_66, gg_68, gg_166, gg_171, gg_173, gg_196, gg_201, \
                         gg_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_32 * gg_61[k]
                  + f_32 * gg_66[k]
                  - f_34 * gg_68[k]
                  + f_32 * gg_166[k]
                  + f_32 * gg_171[k]
                  - f_34 * gg_173[k]
                  - f_33 * gg_196[k]
                  - f_33 * gg_201[k]
                  + f_35 * gg_203[k];
    }

#pragma omp simd aligned(gg_64, gg_71, gg_73, gg_169, gg_176, gg_178, gg_199, gg_206, \
                         gg_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = 5.625 * gg_64[k]
                  + 5.625 * gg_71[k]
                  - 7.5 * gg_73[k]
                  + 5.625 * gg_169[k]
                  + 5.625 * gg_176[k]
                  - 7.5 * gg_178[k]
                  - 7.5 * gg_199[k]
                  - 7.5 * gg_206[k]
                  + 10.0 * gg_208[k];
    }

#pragma omp simd aligned(gg_60, gg_63, gg_65, gg_70, gg_72, gg_74, gg_165, gg_168, gg_170, \
                         gg_175, gg_177, gg_179, gg_195, gg_198, gg_200, gg_205, gg_207, \
                         gg_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_46 * gg_60[k]
                  - f_47 * gg_63[k]
                  + f_48 * gg_65[k]
                  - f_46 * gg_70[k]
                  + f_48 * gg_72[k]
                  - f_49 * gg_74[k]
                  - f_46 * gg_165[k]
                  - f_47 * gg_168[k]
                  + f_48 * gg_170[k]
                  - f_46 * gg_175[k]
                  + f_48 * gg_177[k]
                  - f_49 * gg_179[k]
                  + f_50 * gg_195[k]
                  + f_49 * gg_198[k]
                  - f_51 * gg_200[k]
                  + f_50 * gg_205[k]
                  - f_51 * gg_207[k]
                  + f_52 * gg_209[k];
    }

#pragma omp simd aligned(gg_62, gg_67, gg_69, gg_167, gg_172, gg_174, gg_197, gg_202, \
                         gg_204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = 5.625 * gg_62[k]
                  + 5.625 * gg_67[k]
                  - 7.5 * gg_69[k]
                  + 5.625 * gg_167[k]
                  + 5.625 * gg_172[k]
                  - 7.5 * gg_174[k]
                  - 7.5 * gg_197[k]
                  - 7.5 * gg_202[k]
                  + 10.0 * gg_204[k];
    }

#pragma omp simd aligned(gg_60, gg_65, gg_70, gg_72, gg_165, gg_170, gg_175, gg_177, gg_195, \
                         gg_200, gg_205, gg_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_53 * gg_60[k]
                  - f_54 * gg_65[k]
                  - f_53 * gg_70[k]
                  + f_54 * gg_72[k]
                  + f_53 * gg_165[k]
                  - f_54 * gg_170[k]
                  - f_53 * gg_175[k]
                  + f_54 * gg_177[k]
                  - f_55 * gg_195[k]
                  + f_56 * gg_200[k]
                  + f_55 * gg_205[k]
                  - f_56 * gg_207[k];
    }

#pragma omp simd aligned(gg_62, gg_67, gg_167, gg_172, gg_197, gg_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_16 * gg_62[k]
                  + f_15 * gg_67[k]
                  - f_16 * gg_167[k]
                  + f_15 * gg_172[k]
                  + f_17 * gg_197[k]
                  - f_3 * gg_202[k];
    }

#pragma omp simd aligned(gg_60, gg_63, gg_70, gg_165, gg_168, gg_175, gg_195, gg_198, \
                         gg_205 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_57 * gg_60[k]
                  + f_58 * gg_63[k]
                  - f_57 * gg_70[k]
                  - f_57 * gg_165[k]
                  + f_58 * gg_168[k]
                  - f_57 * gg_175[k]
                  + f_13 * gg_195[k]
                  - f_14 * gg_198[k]
                  + f_13 * gg_205[k];
    }

#pragma omp simd aligned(gg_1, gg_6, gg_46, gg_51, gg_76, gg_81, gg_151, gg_156, gg_181, \
                         gg_186, gg_211, gg_216 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_6 * gg_1[k]
                  - f_6 * gg_6[k]
                  + f_7 * gg_46[k]
                  - f_7 * gg_51[k]
                  - f_8 * gg_76[k]
                  + f_8 * gg_81[k]
                  + f_6 * gg_151[k]
                  - f_6 * gg_156[k]
                  - f_8 * gg_181[k]
                  + f_8 * gg_186[k]
                  + f_9 * gg_211[k]
                  - f_9 * gg_216[k];
    }

#pragma omp simd aligned(gg_4, gg_11, gg_49, gg_56, gg_79, gg_86, gg_154, gg_161, gg_184, \
                         gg_191, gg_214, gg_221 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_18 * gg_4[k]
                  - f_22 * gg_11[k]
                  + f_19 * gg_49[k]
                  - f_23 * gg_56[k]
                  - f_20 * gg_79[k]
                  + f_21 * gg_86[k]
                  + f_18 * gg_154[k]
                  - f_22 * gg_161[k]
                  - f_20 * gg_184[k]
                  + f_21 * gg_191[k]
                  + f_21 * gg_214[k]
                  - f_24 * gg_221[k];
    }

#pragma omp simd aligned(gg_1, gg_6, gg_8, gg_46, gg_51, gg_53, gg_76, gg_81, gg_83, gg_151, \
                         gg_156, gg_158, gg_181, gg_186, gg_188, gg_211, gg_216, \
                         gg_218 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_36 * gg_1[k]
                  - f_36 * gg_6[k]
                  + f_40 * gg_8[k]
                  - f_37 * gg_46[k]
                  - f_37 * gg_51[k]
                  + f_41 * gg_53[k]
                  + f_38 * gg_76[k]
                  + f_38 * gg_81[k]
                  - f_42 * gg_83[k]
                  - f_36 * gg_151[k]
                  - f_36 * gg_156[k]
                  + f_40 * gg_158[k]
                  + f_38 * gg_181[k]
                  + f_38 * gg_186[k]
                  - f_42 * gg_188[k]
                  - f_39 * gg_211[k]
                  - f_39 * gg_216[k]
                  + f_43 * gg_218[k];
    }

#pragma omp simd aligned(gg_4, gg_11, gg_13, gg_49, gg_56, gg_58, gg_79, gg_86, gg_88, gg_154, \
                         gg_161, gg_163, gg_184, gg_191, gg_193, gg_214, gg_221, \
                         gg_223 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_46 * gg_4[k]
                  - f_46 * gg_11[k]
                  + f_50 * gg_13[k]
                  - f_47 * gg_49[k]
                  - f_47 * gg_56[k]
                  + f_49 * gg_58[k]
                  + f_48 * gg_79[k]
                  + f_48 * gg_86[k]
                  - f_51 * gg_88[k]
                  - f_46 * gg_154[k]
                  - f_46 * gg_161[k]
                  + f_50 * gg_163[k]
                  + f_48 * gg_184[k]
                  + f_48 * gg_191[k]
                  - f_51 * gg_193[k]
                  - f_49 * gg_214[k]
                  - f_49 * gg_221[k]
                  + f_52 * gg_223[k];
    }

#pragma omp simd aligned(gg_0, gg_3, gg_5, gg_10, gg_12, gg_14, gg_45, gg_48, gg_50, gg_55, \
                         gg_57, gg_59, gg_75, gg_78, gg_80, gg_85, gg_87, gg_89, gg_150, \
                         gg_153, gg_155, gg_160, gg_162, gg_164, gg_180, gg_183, gg_185, \
                         gg_190, gg_192, gg_194, gg_210, gg_213, gg_215, gg_220, gg_222, \
                         gg_224 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = 0.140625 * gg_0[k]
                  + 0.28125 * gg_3[k]
                  - 1.125 * gg_5[k]
                  + 0.140625 * gg_10[k]
                  - 1.125 * gg_12[k]
                  + 0.375 * gg_14[k]
                  + 0.28125 * gg_45[k]
                  + 0.5625 * gg_48[k]
                  - 2.25 * gg_50[k]
                  + 0.28125 * gg_55[k]
                  - 2.25 * gg_57[k]
                  + 0.75 * gg_59[k]
                  - 1.125 * gg_75[k]
                  - 2.25 * gg_78[k]
                  + 9.0 * gg_80[k]
                  - 1.125 * gg_85[k]
                  + 9.0 * gg_87[k]
                  - 3.0 * gg_89[k]
                  + 0.140625 * gg_150[k]
                  + 0.28125 * gg_153[k]
                  - 1.125 * gg_155[k]
                  + 0.140625 * gg_160[k]
                  - 1.125 * gg_162[k]
                  + 0.375 * gg_164[k]
                  - 1.125 * gg_180[k]
                  - 2.25 * gg_183[k]
                  + 9.0 * gg_185[k]
                  - 1.125 * gg_190[k]
                  + 9.0 * gg_192[k]
                  - 3.0 * gg_194[k]
                  + 0.375 * gg_210[k]
                  + 0.75 * gg_213[k]
                  - 3.0 * gg_215[k]
                  + 0.375 * gg_220[k]
                  - 3.0 * gg_222[k]
                  + gg_224[k];
    }

#pragma omp simd aligned(gg_2, gg_7, gg_9, gg_47, gg_52, gg_54, gg_77, gg_82, gg_84, gg_152, \
                         gg_157, gg_159, gg_182, gg_187, gg_189, gg_212, gg_217, \
                         gg_219 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_46 * gg_2[k]
                  - f_46 * gg_7[k]
                  + f_50 * gg_9[k]
                  - f_47 * gg_47[k]
                  - f_47 * gg_52[k]
                  + f_49 * gg_54[k]
                  + f_48 * gg_77[k]
                  + f_48 * gg_82[k]
                  - f_51 * gg_84[k]
                  - f_46 * gg_152[k]
                  - f_46 * gg_157[k]
                  + f_50 * gg_159[k]
                  + f_48 * gg_182[k]
                  + f_48 * gg_187[k]
                  - f_51 * gg_189[k]
                  - f_49 * gg_212[k]
                  - f_49 * gg_217[k]
                  + f_52 * gg_219[k];
    }

#pragma omp simd aligned(gg_0, gg_5, gg_10, gg_12, gg_45, gg_50, gg_55, gg_57, gg_75, gg_80, \
                         gg_85, gg_87, gg_150, gg_155, gg_160, gg_162, gg_180, gg_185, gg_190, \
                         gg_192, gg_210, gg_215, gg_220, gg_222 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_59 * gg_0[k]
                  + f_60 * gg_5[k]
                  + f_59 * gg_10[k]
                  - f_60 * gg_12[k]
                  - f_36 * gg_45[k]
                  + f_40 * gg_50[k]
                  + f_36 * gg_55[k]
                  - f_40 * gg_57[k]
                  + f_61 * gg_75[k]
                  - f_62 * gg_80[k]
                  - f_61 * gg_85[k]
                  + f_62 * gg_87[k]
                  - f_59 * gg_150[k]
                  + f_60 * gg_155[k]
                  + f_59 * gg_160[k]
                  - f_60 * gg_162[k]
                  + f_61 * gg_180[k]
                  - f_62 * gg_185[k]
                  - f_61 * gg_190[k]
                  + f_62 * gg_192[k]
                  - f_63 * gg_210[k]
                  + f_38 * gg_215[k]
                  + f_63 * gg_220[k]
                  - f_38 * gg_222[k];
    }

#pragma omp simd aligned(gg_2, gg_7, gg_47, gg_52, gg_77, gg_82, gg_152, gg_157, gg_182, \
                         gg_187, gg_212, gg_217 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_22 * gg_2[k]
                  - f_18 * gg_7[k]
                  + f_23 * gg_47[k]
                  - f_19 * gg_52[k]
                  - f_21 * gg_77[k]
                  + f_20 * gg_82[k]
                  + f_22 * gg_152[k]
                  - f_18 * gg_157[k]
                  - f_21 * gg_182[k]
                  + f_20 * gg_187[k]
                  + f_24 * gg_212[k]
                  - f_21 * gg_217[k];
    }

#pragma omp simd aligned(gg_0, gg_3, gg_10, gg_45, gg_48, gg_55, gg_75, gg_78, gg_85, gg_150, \
                         gg_153, gg_160, gg_180, gg_183, gg_190, gg_210, gg_213, \
                         gg_220 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_64 * gg_0[k]
                  - f_65 * gg_3[k]
                  + f_64 * gg_10[k]
                  + f_66 * gg_45[k]
                  - f_67 * gg_48[k]
                  + f_66 * gg_55[k]
                  - f_7 * gg_75[k]
                  + f_68 * gg_78[k]
                  - f_7 * gg_85[k]
                  + f_64 * gg_150[k]
                  - f_65 * gg_153[k]
                  + f_64 * gg_160[k]
                  - f_7 * gg_180[k]
                  + f_68 * gg_183[k]
                  - f_7 * gg_190[k]
                  + f_69 * gg_210[k]
                  - f_70 * gg_213[k]
                  + f_69 * gg_220[k];
    }

#pragma omp simd aligned(gg_31, gg_34, gg_36, gg_41, gg_106, gg_109, gg_111, gg_116, gg_136, \
                         gg_139, gg_141, gg_146 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_4 * gg_31[k]
                  + f_4 * gg_36[k]
                  - f_4 * gg_106[k]
                  + f_4 * gg_111[k]
                  + f_5 * gg_136[k]
                  - f_5 * gg_141[k];

        g_46[k] = -f_15 * gg_34[k]
                  + f_16 * gg_41[k]
                  - f_15 * gg_109[k]
                  + f_16 * gg_116[k]
                  + f_3 * gg_139[k]
                  - f_17 * gg_146[k];
    }

#pragma omp simd aligned(gg_31, gg_36, gg_38, gg_106, gg_111, gg_113, gg_136, gg_141, \
                         gg_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = f_32 * gg_31[k]
                  + f_32 * gg_36[k]
                  - f_34 * gg_38[k]
                  + f_32 * gg_106[k]
                  + f_32 * gg_111[k]
                  - f_34 * gg_113[k]
                  - f_33 * gg_136[k]
                  - f_33 * gg_141[k]
                  + f_35 * gg_143[k];
    }

#pragma omp simd aligned(gg_34, gg_41, gg_43, gg_109, gg_116, gg_118, gg_139, gg_146, \
                         gg_148 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = 5.625 * gg_34[k]
                  + 5.625 * gg_41[k]
                  - 7.5 * gg_43[k]
                  + 5.625 * gg_109[k]
                  + 5.625 * gg_116[k]
                  - 7.5 * gg_118[k]
                  - 7.5 * gg_139[k]
                  - 7.5 * gg_146[k]
                  + 10.0 * gg_148[k];
    }

#pragma omp simd aligned(gg_30, gg_33, gg_35, gg_40, gg_42, gg_44, gg_105, gg_108, gg_110, \
                         gg_115, gg_117, gg_119, gg_135, gg_138, gg_140, gg_145, gg_147, \
                         gg_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_46 * gg_30[k]
                  - f_47 * gg_33[k]
                  + f_48 * gg_35[k]
                  - f_46 * gg_40[k]
                  + f_48 * gg_42[k]
                  - f_49 * gg_44[k]
                  - f_46 * gg_105[k]
                  - f_47 * gg_108[k]
                  + f_48 * gg_110[k]
                  - f_46 * gg_115[k]
                  + f_48 * gg_117[k]
                  - f_49 * gg_119[k]
                  + f_50 * gg_135[k]
                  + f_49 * gg_138[k]
                  - f_51 * gg_140[k]
                  + f_50 * gg_145[k]
                  - f_51 * gg_147[k]
                  + f_52 * gg_149[k];
    }

#pragma omp simd aligned(gg_32, gg_37, gg_39, gg_107, gg_112, gg_114, gg_137, gg_142, \
                         gg_144 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = 5.625 * gg_32[k]
                  + 5.625 * gg_37[k]
                  - 7.5 * gg_39[k]
                  + 5.625 * gg_107[k]
                  + 5.625 * gg_112[k]
                  - 7.5 * gg_114[k]
                  - 7.5 * gg_137[k]
                  - 7.5 * gg_142[k]
                  + 10.0 * gg_144[k];
    }

#pragma omp simd aligned(gg_30, gg_35, gg_40, gg_42, gg_105, gg_110, gg_115, gg_117, gg_135, \
                         gg_140, gg_145, gg_147 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_53 * gg_30[k]
                  - f_54 * gg_35[k]
                  - f_53 * gg_40[k]
                  + f_54 * gg_42[k]
                  + f_53 * gg_105[k]
                  - f_54 * gg_110[k]
                  - f_53 * gg_115[k]
                  + f_54 * gg_117[k]
                  - f_55 * gg_135[k]
                  + f_56 * gg_140[k]
                  + f_55 * gg_145[k]
                  - f_56 * gg_147[k];
    }

#pragma omp simd aligned(gg_32, gg_37, gg_107, gg_112, gg_137, gg_142 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_16 * gg_32[k]
                  + f_15 * gg_37[k]
                  - f_16 * gg_107[k]
                  + f_15 * gg_112[k]
                  + f_17 * gg_137[k]
                  - f_3 * gg_142[k];
    }

#pragma omp simd aligned(gg_30, gg_33, gg_40, gg_105, gg_108, gg_115, gg_135, gg_138, \
                         gg_145 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = -f_57 * gg_30[k]
                  + f_58 * gg_33[k]
                  - f_57 * gg_40[k]
                  - f_57 * gg_105[k]
                  + f_58 * gg_108[k]
                  - f_57 * gg_115[k]
                  + f_13 * gg_135[k]
                  - f_14 * gg_138[k]
                  + f_13 * gg_145[k];
    }

#pragma omp simd aligned(gg_1, gg_6, gg_76, gg_81, gg_151, gg_156, gg_181, \
                         gg_186 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_10 * gg_1[k]
                  + f_10 * gg_6[k]
                  + f_11 * gg_76[k]
                  - f_11 * gg_81[k]
                  + f_10 * gg_151[k]
                  - f_10 * gg_156[k]
                  - f_11 * gg_181[k]
                  + f_11 * gg_186[k];
    }

#pragma omp simd aligned(gg_4, gg_11, gg_79, gg_86, gg_154, gg_161, gg_184, \
                         gg_191 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_25 * gg_4[k]
                  + f_27 * gg_11[k]
                  + f_26 * gg_79[k]
                  - f_4 * gg_86[k]
                  + f_25 * gg_154[k]
                  - f_27 * gg_161[k]
                  - f_26 * gg_184[k]
                  + f_4 * gg_191[k];
    }

#pragma omp simd aligned(gg_1, gg_6, gg_8, gg_76, gg_81, gg_83, gg_151, gg_156, gg_158, \
                         gg_181, gg_186, gg_188 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = 0.625 * gg_1[k]
                  + 0.625 * gg_6[k]
                  - 3.75 * gg_8[k]
                  - 3.75 * gg_76[k]
                  - 3.75 * gg_81[k]
                  + 22.5 * gg_83[k]
                  - 0.625 * gg_151[k]
                  - 0.625 * gg_156[k]
                  + 3.75 * gg_158[k]
                  + 3.75 * gg_181[k]
                  + 3.75 * gg_186[k]
                  - 22.5 * gg_188[k];
    }

#pragma omp simd aligned(gg_4, gg_11, gg_13, gg_79, gg_86, gg_88, gg_154, gg_161, gg_163, \
                         gg_184, gg_191, gg_193 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_53 * gg_4[k]
                  + f_53 * gg_11[k]
                  - f_55 * gg_13[k]
                  - f_54 * gg_79[k]
                  - f_54 * gg_86[k]
                  + f_56 * gg_88[k]
                  - f_53 * gg_154[k]
                  - f_53 * gg_161[k]
                  + f_55 * gg_163[k]
                  + f_54 * gg_184[k]
                  + f_54 * gg_191[k]
                  - f_56 * gg_193[k];
    }

#pragma omp simd aligned(gg_0, gg_3, gg_5, gg_10, gg_12, gg_14, gg_75, gg_78, gg_80, gg_85, \
                         gg_87, gg_89, gg_150, gg_153, gg_155, gg_160, gg_162, gg_164, gg_180, \
                         gg_183, gg_185, gg_190, gg_192, gg_194 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_59 * gg_0[k]
                  - f_36 * gg_3[k]
                  + f_61 * gg_5[k]
                  - f_59 * gg_10[k]
                  + f_61 * gg_12[k]
                  - f_63 * gg_14[k]
                  + f_60 * gg_75[k]
                  + f_40 * gg_78[k]
                  - f_62 * gg_80[k]
                  + f_60 * gg_85[k]
                  - f_62 * gg_87[k]
                  + f_38 * gg_89[k]
                  + f_59 * gg_150[k]
                  + f_36 * gg_153[k]
                  - f_61 * gg_155[k]
                  + f_59 * gg_160[k]
                  - f_61 * gg_162[k]
                  + f_63 * gg_164[k]
                  - f_60 * gg_180[k]
                  - f_40 * gg_183[k]
                  + f_62 * gg_185[k]
                  - f_60 * gg_190[k]
                  + f_62 * gg_192[k]
                  - f_38 * gg_194[k];
    }

#pragma omp simd aligned(gg_2, gg_7, gg_9, gg_77, gg_82, gg_84, gg_152, gg_157, gg_159, \
                         gg_182, gg_187, gg_189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_53 * gg_2[k]
                  + f_53 * gg_7[k]
                  - f_55 * gg_9[k]
                  - f_54 * gg_77[k]
                  - f_54 * gg_82[k]
                  + f_56 * gg_84[k]
                  - f_53 * gg_152[k]
                  - f_53 * gg_157[k]
                  + f_55 * gg_159[k]
                  + f_54 * gg_182[k]
                  + f_54 * gg_187[k]
                  - f_56 * gg_189[k];
    }

#pragma omp simd aligned(gg_0, gg_5, gg_10, gg_12, gg_75, gg_80, gg_85, gg_87, gg_150, gg_155, \
                         gg_160, gg_162, gg_180, gg_185, gg_190, \
                         gg_192 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = 0.3125 * gg_0[k]
                  - 1.875 * gg_5[k]
                  - 0.3125 * gg_10[k]
                  + 1.875 * gg_12[k]
                  - 1.875 * gg_75[k]
                  + 11.25 * gg_80[k]
                  + 1.875 * gg_85[k]
                  - 11.25 * gg_87[k]
                  - 0.3125 * gg_150[k]
                  + 1.875 * gg_155[k]
                  + 0.3125 * gg_160[k]
                  - 1.875 * gg_162[k]
                  + 1.875 * gg_180[k]
                  - 11.25 * gg_185[k]
                  - 1.875 * gg_190[k]
                  + 11.25 * gg_192[k];
    }

#pragma omp simd aligned(gg_2, gg_7, gg_77, gg_82, gg_152, gg_157, gg_182, \
                         gg_187 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_27 * gg_2[k]
                  + f_25 * gg_7[k]
                  + f_4 * gg_77[k]
                  - f_26 * gg_82[k]
                  + f_27 * gg_152[k]
                  - f_25 * gg_157[k]
                  - f_4 * gg_182[k]
                  + f_26 * gg_187[k];
    }

#pragma omp simd aligned(gg_0, gg_3, gg_10, gg_75, gg_78, gg_85, gg_150, gg_153, gg_160, \
                         gg_180, gg_183, gg_190 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_71 * gg_0[k]
                  + f_72 * gg_3[k]
                  - f_71 * gg_10[k]
                  + f_72 * gg_75[k]
                  - f_15 * gg_78[k]
                  + f_72 * gg_85[k]
                  + f_71 * gg_150[k]
                  - f_72 * gg_153[k]
                  + f_71 * gg_160[k]
                  - f_72 * gg_180[k]
                  + f_15 * gg_183[k]
                  - f_72 * gg_190[k];
    }

#pragma omp simd aligned(gg_31, gg_34, gg_36, gg_38, gg_41, gg_43, gg_106, gg_109, gg_111, \
                         gg_113, gg_116, gg_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_1 * gg_31[k]
                  - f_1 * gg_36[k]
                  - f_0 * gg_106[k]
                  + f_0 * gg_111[k];

        g_64[k] = 13.125 * gg_34[k]
                  - 4.375 * gg_41[k]
                  - 39.375 * gg_109[k]
                  + 13.125 * gg_116[k];

        g_65[k] = -f_13 * gg_31[k]
                  - f_13 * gg_36[k]
                  + f_14 * gg_38[k]
                  + f_4 * gg_106[k]
                  + f_4 * gg_111[k]
                  - f_12 * gg_113[k];

        g_66[k] = -f_16 * gg_34[k]
                  - f_16 * gg_41[k]
                  + f_17 * gg_43[k]
                  + f_15 * gg_109[k]
                  + f_15 * gg_116[k]
                  - f_3 * gg_118[k];
    }

#pragma omp simd aligned(gg_30, gg_33, gg_35, gg_40, gg_42, gg_44, gg_105, gg_108, gg_110, \
                         gg_115, gg_117, gg_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = f_22 * gg_30[k]
                  + f_23 * gg_33[k]
                  - f_21 * gg_35[k]
                  + f_22 * gg_40[k]
                  - f_21 * gg_42[k]
                  + f_24 * gg_44[k]
                  - f_18 * gg_105[k]
                  - f_19 * gg_108[k]
                  + f_20 * gg_110[k]
                  - f_18 * gg_115[k]
                  + f_20 * gg_117[k]
                  - f_21 * gg_119[k];
    }

#pragma omp simd aligned(gg_30, gg_32, gg_35, gg_37, gg_39, gg_40, gg_42, gg_105, gg_107, \
                         gg_110, gg_112, gg_114, gg_115, gg_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = -f_16 * gg_32[k]
                  - f_16 * gg_37[k]
                  + f_17 * gg_39[k]
                  + f_15 * gg_107[k]
                  + f_15 * gg_112[k]
                  - f_3 * gg_114[k];

        g_69[k] = -f_27 * gg_30[k]
                  + f_4 * gg_35[k]
                  + f_27 * gg_40[k]
                  - f_4 * gg_42[k]
                  + f_25 * gg_105[k]
                  - f_26 * gg_110[k]
                  - f_25 * gg_115[k]
                  + f_26 * gg_117[k];
    }

#pragma omp simd aligned(gg_30, gg_32, gg_33, gg_37, gg_40, gg_105, gg_107, gg_108, gg_112, \
                         gg_115 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = 4.375 * gg_32[k]
                  - 13.125 * gg_37[k]
                  - 13.125 * gg_107[k]
                  + 39.375 * gg_112[k];

        g_71[k] = f_30 * gg_30[k]
                  - f_31 * gg_33[k]
                  + f_30 * gg_40[k]
                  - f_28 * gg_105[k]
                  + f_29 * gg_108[k]
                  - f_28 * gg_115[k];
    }

#pragma omp simd aligned(gg_1, gg_4, gg_6, gg_11, gg_46, gg_49, gg_51, gg_56, gg_151, gg_154, \
                         gg_156, gg_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = 2.1875 * gg_1[k]
                  - 2.1875 * gg_6[k]
                  - 13.125 * gg_46[k]
                  + 13.125 * gg_51[k]
                  + 2.1875 * gg_151[k]
                  - 2.1875 * gg_156[k];

        g_73[k] = f_28 * gg_4[k]
                  - f_30 * gg_11[k]
                  - f_29 * gg_49[k]
                  + f_31 * gg_56[k]
                  + f_28 * gg_154[k]
                  - f_30 * gg_161[k];
    }

#pragma omp simd aligned(gg_1, gg_6, gg_8, gg_46, gg_51, gg_53, gg_151, gg_156, \
                         gg_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_44 * gg_1[k]
                  - f_44 * gg_6[k]
                  + f_16 * gg_8[k]
                  + f_16 * gg_46[k]
                  + f_16 * gg_51[k]
                  - f_45 * gg_53[k]
                  - f_44 * gg_151[k]
                  - f_44 * gg_156[k]
                  + f_16 * gg_158[k];
    }

#pragma omp simd aligned(gg_4, gg_11, gg_13, gg_49, gg_56, gg_58, gg_154, gg_161, \
                         gg_163 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_57 * gg_4[k]
                  - f_57 * gg_11[k]
                  + f_13 * gg_13[k]
                  + f_58 * gg_49[k]
                  + f_58 * gg_56[k]
                  - f_14 * gg_58[k]
                  - f_57 * gg_154[k]
                  - f_57 * gg_161[k]
                  + f_13 * gg_163[k];
    }

#pragma omp simd aligned(gg_0, gg_3, gg_5, gg_10, gg_12, gg_14, gg_45, gg_48, gg_50, gg_55, \
                         gg_57, gg_59, gg_150, gg_153, gg_155, gg_160, gg_162, \
                         gg_164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = f_64 * gg_0[k]
                  + f_66 * gg_3[k]
                  - f_7 * gg_5[k]
                  + f_64 * gg_10[k]
                  - f_7 * gg_12[k]
                  + f_69 * gg_14[k]
                  - f_65 * gg_45[k]
                  - f_67 * gg_48[k]
                  + f_68 * gg_50[k]
                  - f_65 * gg_55[k]
                  + f_68 * gg_57[k]
                  - f_70 * gg_59[k]
                  + f_64 * gg_150[k]
                  + f_66 * gg_153[k]
                  - f_7 * gg_155[k]
                  + f_64 * gg_160[k]
                  - f_7 * gg_162[k]
                  + f_69 * gg_164[k];
    }

#pragma omp simd aligned(gg_2, gg_7, gg_9, gg_47, gg_52, gg_54, gg_152, gg_157, \
                         gg_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = -f_57 * gg_2[k]
                  - f_57 * gg_7[k]
                  + f_13 * gg_9[k]
                  + f_58 * gg_47[k]
                  + f_58 * gg_52[k]
                  - f_14 * gg_54[k]
                  - f_57 * gg_152[k]
                  - f_57 * gg_157[k]
                  + f_13 * gg_159[k];
    }

#pragma omp simd aligned(gg_0, gg_5, gg_10, gg_12, gg_45, gg_50, gg_55, gg_57, gg_150, gg_155, \
                         gg_160, gg_162 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = -f_71 * gg_0[k]
                  + f_72 * gg_5[k]
                  + f_71 * gg_10[k]
                  - f_72 * gg_12[k]
                  + f_72 * gg_45[k]
                  - f_15 * gg_50[k]
                  - f_72 * gg_55[k]
                  + f_15 * gg_57[k]
                  - f_71 * gg_150[k]
                  + f_72 * gg_155[k]
                  + f_71 * gg_160[k]
                  - f_72 * gg_162[k];
    }

#pragma omp simd aligned(gg_2, gg_7, gg_47, gg_52, gg_152, gg_157 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = f_30 * gg_2[k]
                  - f_28 * gg_7[k]
                  - f_31 * gg_47[k]
                  + f_29 * gg_52[k]
                  + f_30 * gg_152[k]
                  - f_28 * gg_157[k];
    }

#pragma omp simd aligned(gg_0, gg_3, gg_10, gg_45, gg_48, gg_55, gg_150, gg_153, \
                         gg_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = 0.546875 * gg_0[k]
                  - 3.28125 * gg_3[k]
                  + 0.546875 * gg_10[k]
                  - 3.28125 * gg_45[k]
                  + 19.6875 * gg_48[k]
                  - 3.28125 * gg_55[k]
                  + 0.546875 * gg_150[k]
                  - 3.28125 * gg_153[k]
                  + 0.546875 * gg_160[k];
    }
}

auto
transform_gg_tri(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t gg,
                 const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 13.125 * std::sqrt(2.0);
    const auto f_1 = 4.375 * std::sqrt(2.0);
    const auto f_2 = 1.25 * std::sqrt(7.0);
    const auto f_3 = 7.5 * std::sqrt(7.0);
    const auto f_4 = 1.875 * std::sqrt(14.0);
    const auto f_5 = 2.5 * std::sqrt(14.0);
    const auto f_6 = 0.1875 * std::sqrt(35.0);
    const auto f_7 = 0.375 * std::sqrt(35.0);
    const auto f_8 = 1.5 * std::sqrt(35.0);
    const auto f_9 = 0.5 * std::sqrt(35.0);
    const auto f_10 = 0.625 * std::sqrt(7.0);
    const auto f_11 = 3.75 * std::sqrt(7.0);
    const auto f_12 = 11.25 * std::sqrt(14.0);
    const auto f_13 = 0.625 * std::sqrt(14.0);
    const auto f_14 = 3.75 * std::sqrt(14.0);
    const auto f_15 = 5.625 * std::sqrt(7.0);
    const auto f_16 = 1.875 * std::sqrt(7.0);
    const auto f_17 = 2.5 * std::sqrt(7.0);
    const auto f_18 = 0.28125 * std::sqrt(70.0);
    const auto f_19 = 0.5625 * std::sqrt(70.0);
    const auto f_20 = 2.25 * std::sqrt(70.0);
    const auto f_21 = 0.75 * std::sqrt(70.0);
    const auto f_22 = 0.09375 * std::sqrt(70.0);
    const auto f_23 = 0.1875 * std::sqrt(70.0);
    const auto f_24 = 0.25 * std::sqrt(70.0);
    const auto f_25 = 0.9375 * std::sqrt(14.0);
    const auto f_26 = 5.625 * std::sqrt(14.0);
    const auto f_27 = 0.3125 * std::sqrt(14.0);
    const auto f_28 = 3.28125 * std::sqrt(2.0);
    const auto f_29 = 19.6875 * std::sqrt(2.0);
    const auto f_30 = 1.09375 * std::sqrt(2.0);
    const auto f_31 = 6.5625 * std::sqrt(2.0);
    const auto f_32 = 1.875 * std::sqrt(2.0);
    const auto f_33 = 2.5 * std::sqrt(2.0);
    const auto f_34 = 11.25 * std::sqrt(2.0);
    const auto f_35 = 15.0 * std::sqrt(2.0);
    const auto f_36 = 0.1875 * std::sqrt(5.0);
    const auto f_37 = 0.375 * std::sqrt(5.0);
    const auto f_38 = 1.5 * std::sqrt(5.0);
    const auto f_39 = 0.5 * std::sqrt(5.0);
    const auto f_40 = 1.125 * std::sqrt(5.0);
    const auto f_41 = 2.25 * std::sqrt(5.0);
    const auto f_42 = 9.0 * std::sqrt(5.0);
    const auto f_43 = 3.0 * std::sqrt(5.0);
    const auto f_44 = 0.3125 * std::sqrt(7.0);
    const auto f_45 = 11.25 * std::sqrt(7.0);
    const auto f_46 = 0.28125 * std::sqrt(10.0);
    const auto f_47 = 0.5625 * std::sqrt(10.0);
    const auto f_48 = 2.25 * std::sqrt(10.0);
    const auto f_49 = 0.75 * std::sqrt(10.0);
    const auto f_50 = 0.375 * std::sqrt(10.0);
    const auto f_51 = 3.0 * std::sqrt(10.0);
    const auto f_52 = std::sqrt(10.0);
    const auto f_53 = 0.9375 * std::sqrt(2.0);
    const auto f_54 = 5.625 * std::sqrt(2.0);
    const auto f_55 = 1.25 * std::sqrt(2.0);
    const auto f_56 = 7.5 * std::sqrt(2.0);
    const auto f_57 = 0.46875 * std::sqrt(14.0);
    const auto f_58 = 2.8125 * std::sqrt(14.0);
    const auto f_59 = 0.09375 * std::sqrt(5.0);
    const auto f_60 = 0.5625 * std::sqrt(5.0);
    const auto f_61 = 0.75 * std::sqrt(5.0);
    const auto f_62 = 4.5 * std::sqrt(5.0);
    const auto f_63 = 0.25 * std::sqrt(5.0);
    const auto f_64 = 0.046875 * std::sqrt(35.0);
    const auto f_65 = 0.28125 * std::sqrt(35.0);
    const auto f_66 = 0.09375 * std::sqrt(35.0);
    const auto f_67 = 0.5625 * std::sqrt(35.0);
    const auto f_68 = 2.25 * std::sqrt(35.0);
    const auto f_69 = 0.125 * std::sqrt(35.0);
    const auto f_70 = 0.75 * std::sqrt(35.0);
    const auto f_71 = 0.15625 * std::sqrt(7.0);
    const auto f_72 = 0.9375 * std::sqrt(7.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;
    auto *g_7 = values + 7 * nvalues;
    auto *g_8 = values + 8 * nvalues;
    auto *g_9 = values + 9 * nvalues;
    auto *g_10 = values + 10 * nvalues;
    auto *g_11 = values + 11 * nvalues;
    auto *g_12 = values + 12 * nvalues;
    auto *g_13 = values + 13 * nvalues;
    auto *g_14 = values + 14 * nvalues;
    auto *g_15 = values + 15 * nvalues;
    auto *g_16 = values + 16 * nvalues;
    auto *g_17 = values + 17 * nvalues;
    auto *g_18 = values + 18 * nvalues;
    auto *g_19 = values + 19 * nvalues;
    auto *g_20 = values + 20 * nvalues;
    auto *g_21 = values + 21 * nvalues;
    auto *g_22 = values + 22 * nvalues;
    auto *g_23 = values + 23 * nvalues;
    auto *g_24 = values + 24 * nvalues;
    auto *g_25 = values + 25 * nvalues;
    auto *g_26 = values + 26 * nvalues;
    auto *g_27 = values + 27 * nvalues;
    auto *g_28 = values + 28 * nvalues;
    auto *g_29 = values + 29 * nvalues;
    auto *g_30 = values + 30 * nvalues;
    auto *g_31 = values + 31 * nvalues;
    auto *g_32 = values + 32 * nvalues;
    auto *g_33 = values + 33 * nvalues;
    auto *g_34 = values + 34 * nvalues;
    auto *g_35 = values + 35 * nvalues;
    auto *g_36 = values + 36 * nvalues;
    auto *g_37 = values + 37 * nvalues;
    auto *g_38 = values + 38 * nvalues;
    auto *g_39 = values + 39 * nvalues;
    auto *g_40 = values + 40 * nvalues;
    auto *g_41 = values + 41 * nvalues;
    auto *g_42 = values + 42 * nvalues;
    auto *g_43 = values + 43 * nvalues;
    auto *g_44 = values + 44 * nvalues;
    auto *g_45 = values + 45 * nvalues;
    auto *g_46 = values + 46 * nvalues;
    auto *g_47 = values + 47 * nvalues;
    auto *g_48 = values + 48 * nvalues;
    auto *g_49 = values + 49 * nvalues;
    auto *g_50 = values + 50 * nvalues;
    auto *g_51 = values + 51 * nvalues;
    auto *g_52 = values + 52 * nvalues;
    auto *g_53 = values + 53 * nvalues;
    auto *g_54 = values + 54 * nvalues;
    auto *g_55 = values + 55 * nvalues;
    auto *g_56 = values + 56 * nvalues;
    auto *g_57 = values + 57 * nvalues;
    auto *g_58 = values + 58 * nvalues;
    auto *g_59 = values + 59 * nvalues;
    auto *g_60 = values + 60 * nvalues;
    auto *g_61 = values + 61 * nvalues;
    auto *g_62 = values + 62 * nvalues;
    auto *g_63 = values + 63 * nvalues;
    auto *g_64 = values + 64 * nvalues;
    auto *g_65 = values + 65 * nvalues;
    auto *g_66 = values + 66 * nvalues;
    auto *g_67 = values + 67 * nvalues;
    auto *g_68 = values + 68 * nvalues;
    auto *g_69 = values + 69 * nvalues;
    auto *g_70 = values + 70 * nvalues;
    auto *g_71 = values + 71 * nvalues;
    auto *g_72 = values + 72 * nvalues;
    auto *g_73 = values + 73 * nvalues;
    auto *g_74 = values + 74 * nvalues;
    auto *g_75 = values + 75 * nvalues;
    auto *g_76 = values + 76 * nvalues;
    auto *g_77 = values + 77 * nvalues;
    auto *g_78 = values + 78 * nvalues;
    auto *g_79 = values + 79 * nvalues;
    auto *g_80 = values + 80 * nvalues;

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_82 = buffer.data(gg + 82);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_91 = buffer.data(gg + 91);
    const auto *gg_92 = buffer.data(gg + 92);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_94 = buffer.data(gg + 94);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_96 = buffer.data(gg + 96);
    const auto *gg_97 = buffer.data(gg + 97);
    const auto *gg_98 = buffer.data(gg + 98);
    const auto *gg_99 = buffer.data(gg + 99);
    const auto *gg_100 = buffer.data(gg + 100);
    const auto *gg_101 = buffer.data(gg + 101);
    const auto *gg_102 = buffer.data(gg + 102);
    const auto *gg_103 = buffer.data(gg + 103);
    const auto *gg_104 = buffer.data(gg + 104);
    const auto *gg_105 = buffer.data(gg + 105);
    const auto *gg_107 = buffer.data(gg + 107);
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_110 = buffer.data(gg + 110);
    const auto *gg_112 = buffer.data(gg + 112);
    const auto *gg_114 = buffer.data(gg + 114);
    const auto *gg_115 = buffer.data(gg + 115);
    const auto *gg_117 = buffer.data(gg + 117);
    const auto *gg_120 = buffer.data(gg + 120);
    const auto *gg_121 = buffer.data(gg + 121);
    const auto *gg_122 = buffer.data(gg + 122);
    const auto *gg_123 = buffer.data(gg + 123);
    const auto *gg_124 = buffer.data(gg + 124);
    const auto *gg_125 = buffer.data(gg + 125);
    const auto *gg_126 = buffer.data(gg + 126);
    const auto *gg_127 = buffer.data(gg + 127);
    const auto *gg_128 = buffer.data(gg + 128);
    const auto *gg_129 = buffer.data(gg + 129);
    const auto *gg_130 = buffer.data(gg + 130);
    const auto *gg_131 = buffer.data(gg + 131);
    const auto *gg_132 = buffer.data(gg + 132);
    const auto *gg_133 = buffer.data(gg + 133);
    const auto *gg_134 = buffer.data(gg + 134);
    const auto *gg_135 = buffer.data(gg + 135);
    const auto *gg_137 = buffer.data(gg + 137);
    const auto *gg_138 = buffer.data(gg + 138);
    const auto *gg_140 = buffer.data(gg + 140);
    const auto *gg_142 = buffer.data(gg + 142);
    const auto *gg_144 = buffer.data(gg + 144);
    const auto *gg_145 = buffer.data(gg + 145);
    const auto *gg_147 = buffer.data(gg + 147);
    const auto *gg_150 = buffer.data(gg + 150);
    const auto *gg_152 = buffer.data(gg + 152);
    const auto *gg_153 = buffer.data(gg + 153);
    const auto *gg_155 = buffer.data(gg + 155);
    const auto *gg_157 = buffer.data(gg + 157);
    const auto *gg_159 = buffer.data(gg + 159);
    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_162 = buffer.data(gg + 162);
    const auto *gg_164 = buffer.data(gg + 164);
    const auto *gg_165 = buffer.data(gg + 165);
    const auto *gg_166 = buffer.data(gg + 166);
    const auto *gg_167 = buffer.data(gg + 167);
    const auto *gg_168 = buffer.data(gg + 168);
    const auto *gg_169 = buffer.data(gg + 169);
    const auto *gg_170 = buffer.data(gg + 170);
    const auto *gg_171 = buffer.data(gg + 171);
    const auto *gg_172 = buffer.data(gg + 172);
    const auto *gg_173 = buffer.data(gg + 173);
    const auto *gg_174 = buffer.data(gg + 174);
    const auto *gg_175 = buffer.data(gg + 175);
    const auto *gg_176 = buffer.data(gg + 176);
    const auto *gg_177 = buffer.data(gg + 177);
    const auto *gg_178 = buffer.data(gg + 178);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_180 = buffer.data(gg + 180);
    const auto *gg_182 = buffer.data(gg + 182);
    const auto *gg_183 = buffer.data(gg + 183);
    const auto *gg_185 = buffer.data(gg + 185);
    const auto *gg_187 = buffer.data(gg + 187);
    const auto *gg_189 = buffer.data(gg + 189);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_195 = buffer.data(gg + 195);
    const auto *gg_197 = buffer.data(gg + 197);
    const auto *gg_198 = buffer.data(gg + 198);
    const auto *gg_199 = buffer.data(gg + 199);
    const auto *gg_200 = buffer.data(gg + 200);
    const auto *gg_202 = buffer.data(gg + 202);
    const auto *gg_204 = buffer.data(gg + 204);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_206 = buffer.data(gg + 206);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_208 = buffer.data(gg + 208);
    const auto *gg_209 = buffer.data(gg + 209);
    const auto *gg_210 = buffer.data(gg + 210);
    const auto *gg_212 = buffer.data(gg + 212);
    const auto *gg_213 = buffer.data(gg + 213);
    const auto *gg_215 = buffer.data(gg + 215);
    const auto *gg_217 = buffer.data(gg + 217);
    const auto *gg_219 = buffer.data(gg + 219);
    const auto *gg_220 = buffer.data(gg + 220);
    const auto *gg_222 = buffer.data(gg + 222);
    const auto *gg_224 = buffer.data(gg + 224);

#pragma omp simd aligned(gg_16, gg_19, gg_21, gg_23, gg_26, gg_28, gg_91, gg_94, gg_96, gg_98, \
                         gg_101, gg_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = 8.75 * gg_16[k]
                 - 8.75 * gg_21[k]
                 - 8.75 * gg_91[k]
                 + 8.75 * gg_96[k];

        g_1[k] = f_0 * gg_19[k]
                 - f_1 * gg_26[k]
                 - f_0 * gg_94[k]
                 + f_1 * gg_101[k];
        g_9[k] = g_1[k];

        g_2[k] = -f_2 * gg_16[k]
                 - f_2 * gg_21[k]
                 + f_3 * gg_23[k]
                 + f_2 * gg_91[k]
                 + f_2 * gg_96[k]
                 - f_3 * gg_98[k];
        g_18[k] = g_2[k];

        g_3[k] = -f_4 * gg_19[k]
                 - f_4 * gg_26[k]
                 + f_5 * gg_28[k]
                 + f_4 * gg_94[k]
                 + f_4 * gg_101[k]
                 - f_5 * gg_103[k];
        g_27[k] = g_3[k];
    }

#pragma omp simd aligned(gg_15, gg_18, gg_20, gg_25, gg_27, gg_29, gg_90, gg_93, gg_95, \
                         gg_100, gg_102, gg_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_6 * gg_15[k]
                 + f_7 * gg_18[k]
                 - f_8 * gg_20[k]
                 + f_6 * gg_25[k]
                 - f_8 * gg_27[k]
                 + f_9 * gg_29[k]
                 - f_6 * gg_90[k]
                 - f_7 * gg_93[k]
                 + f_8 * gg_95[k]
                 - f_6 * gg_100[k]
                 + f_8 * gg_102[k]
                 - f_9 * gg_104[k];
        g_36[k] = g_4[k];
    }

#pragma omp simd aligned(gg_15, gg_17, gg_20, gg_22, gg_24, gg_25, gg_27, gg_90, gg_92, gg_95, \
                         gg_97, gg_99, gg_100, gg_102 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -f_4 * gg_17[k]
                 - f_4 * gg_22[k]
                 + f_5 * gg_24[k]
                 + f_4 * gg_92[k]
                 + f_4 * gg_97[k]
                 - f_5 * gg_99[k];
        g_45[k] = g_5[k];

        g_6[k] = -f_10 * gg_15[k]
                 + f_11 * gg_20[k]
                 + f_10 * gg_25[k]
                 - f_11 * gg_27[k]
                 + f_10 * gg_90[k]
                 - f_11 * gg_95[k]
                 - f_10 * gg_100[k]
                 + f_11 * gg_102[k];
        g_54[k] = g_6[k];
    }

#pragma omp simd aligned(gg_15, gg_17, gg_18, gg_22, gg_25, gg_90, gg_92, gg_93, gg_97, \
                         gg_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_1 * gg_17[k]
                 - f_0 * gg_22[k]
                 - f_1 * gg_92[k]
                 + f_0 * gg_97[k];
        g_63[k] = g_7[k];

        g_8[k] = 2.1875 * gg_15[k]
                 - 13.125 * gg_18[k]
                 + 2.1875 * gg_25[k]
                 - 2.1875 * gg_90[k]
                 + 13.125 * gg_93[k]
                 - 2.1875 * gg_100[k];
        g_72[k] = g_8[k];
    }

#pragma omp simd aligned(gg_61, gg_64, gg_66, gg_68, gg_71, gg_73, gg_166, gg_169, gg_171, \
                         gg_173, gg_176, gg_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = 39.375 * gg_64[k]
                  - 13.125 * gg_71[k]
                  - 13.125 * gg_169[k]
                  + 4.375 * gg_176[k];

        g_11[k] = -f_4 * gg_61[k]
                  - f_4 * gg_66[k]
                  + f_12 * gg_68[k]
                  + f_13 * gg_166[k]
                  + f_13 * gg_171[k]
                  - f_14 * gg_173[k];
        g_19[k] = g_11[k];

        g_12[k] = -f_15 * gg_64[k]
                  - f_15 * gg_71[k]
                  + f_3 * gg_73[k]
                  + f_16 * gg_169[k]
                  + f_16 * gg_176[k]
                  - f_17 * gg_178[k];
        g_28[k] = g_12[k];
    }

#pragma omp simd aligned(gg_60, gg_63, gg_65, gg_70, gg_72, gg_74, gg_165, gg_168, gg_170, \
                         gg_175, gg_177, gg_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_18 * gg_60[k]
                  + f_19 * gg_63[k]
                  - f_20 * gg_65[k]
                  + f_18 * gg_70[k]
                  - f_20 * gg_72[k]
                  + f_21 * gg_74[k]
                  - f_22 * gg_165[k]
                  - f_23 * gg_168[k]
                  + f_21 * gg_170[k]
                  - f_22 * gg_175[k]
                  + f_21 * gg_177[k]
                  - f_24 * gg_179[k];
        g_37[k] = g_13[k];
    }

#pragma omp simd aligned(gg_60, gg_62, gg_65, gg_67, gg_69, gg_70, gg_72, gg_165, gg_167, \
                         gg_170, gg_172, gg_174, gg_175, gg_177 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_15 * gg_62[k]
                  - f_15 * gg_67[k]
                  + f_3 * gg_69[k]
                  + f_16 * gg_167[k]
                  + f_16 * gg_172[k]
                  - f_17 * gg_174[k];
        g_46[k] = g_14[k];

        g_15[k] = -f_25 * gg_60[k]
                  + f_26 * gg_65[k]
                  + f_25 * gg_70[k]
                  - f_26 * gg_72[k]
                  + f_27 * gg_165[k]
                  - f_4 * gg_170[k]
                  - f_27 * gg_175[k]
                  + f_4 * gg_177[k];
        g_55[k] = g_15[k];
    }

#pragma omp simd aligned(gg_60, gg_62, gg_63, gg_67, gg_70, gg_165, gg_167, gg_168, gg_172, \
                         gg_175 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = 13.125 * gg_62[k]
                  - 39.375 * gg_67[k]
                  - 4.375 * gg_167[k]
                  + 13.125 * gg_172[k];
        g_64[k] = g_16[k];

        g_17[k] = f_28 * gg_60[k]
                  - f_29 * gg_63[k]
                  + f_28 * gg_70[k]
                  - f_30 * gg_165[k]
                  + f_31 * gg_168[k]
                  - f_30 * gg_175[k];
        g_73[k] = g_17[k];
    }

#pragma omp simd aligned(gg_16, gg_21, gg_23, gg_91, gg_96, gg_98, gg_121, gg_126, \
                         gg_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = 1.25 * gg_16[k]
                  + 1.25 * gg_21[k]
                  - 7.5 * gg_23[k]
                  + 1.25 * gg_91[k]
                  + 1.25 * gg_96[k]
                  - 7.5 * gg_98[k]
                  - 7.5 * gg_121[k]
                  - 7.5 * gg_126[k]
                  + 45.0 * gg_128[k];
    }

#pragma omp simd aligned(gg_19, gg_26, gg_28, gg_94, gg_101, gg_103, gg_124, gg_131, \
                         gg_133 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_32 * gg_19[k]
                  + f_32 * gg_26[k]
                  - f_33 * gg_28[k]
                  + f_32 * gg_94[k]
                  + f_32 * gg_101[k]
                  - f_33 * gg_103[k]
                  - f_34 * gg_124[k]
                  - f_34 * gg_131[k]
                  + f_35 * gg_133[k];
        g_29[k] = g_21[k];
    }

#pragma omp simd aligned(gg_15, gg_18, gg_20, gg_25, gg_27, gg_29, gg_90, gg_93, gg_95, \
                         gg_100, gg_102, gg_104, gg_120, gg_123, gg_125, gg_130, gg_132, \
                         gg_134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_36 * gg_15[k]
                  - f_37 * gg_18[k]
                  + f_38 * gg_20[k]
                  - f_36 * gg_25[k]
                  + f_38 * gg_27[k]
                  - f_39 * gg_29[k]
                  - f_36 * gg_90[k]
                  - f_37 * gg_93[k]
                  + f_38 * gg_95[k]
                  - f_36 * gg_100[k]
                  + f_38 * gg_102[k]
                  - f_39 * gg_104[k]
                  + f_40 * gg_120[k]
                  + f_41 * gg_123[k]
                  - f_42 * gg_125[k]
                  + f_40 * gg_130[k]
                  - f_42 * gg_132[k]
                  + f_43 * gg_134[k];
        g_38[k] = g_22[k];
    }

#pragma omp simd aligned(gg_17, gg_22, gg_24, gg_92, gg_97, gg_99, gg_122, gg_127, \
                         gg_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_32 * gg_17[k]
                  + f_32 * gg_22[k]
                  - f_33 * gg_24[k]
                  + f_32 * gg_92[k]
                  + f_32 * gg_97[k]
                  - f_33 * gg_99[k]
                  - f_34 * gg_122[k]
                  - f_34 * gg_127[k]
                  + f_35 * gg_129[k];
        g_47[k] = g_23[k];
    }

#pragma omp simd aligned(gg_15, gg_20, gg_25, gg_27, gg_90, gg_95, gg_100, gg_102, gg_120, \
                         gg_125, gg_130, gg_132 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = 0.625 * gg_15[k]
                  - 3.75 * gg_20[k]
                  - 0.625 * gg_25[k]
                  + 3.75 * gg_27[k]
                  + 0.625 * gg_90[k]
                  - 3.75 * gg_95[k]
                  - 0.625 * gg_100[k]
                  + 3.75 * gg_102[k]
                  - 3.75 * gg_120[k]
                  + 22.5 * gg_125[k]
                  + 3.75 * gg_130[k]
                  - 22.5 * gg_132[k];
        g_56[k] = g_24[k];
    }

#pragma omp simd aligned(gg_17, gg_22, gg_92, gg_97, gg_122, gg_127 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_13 * gg_17[k]
                  + f_4 * gg_22[k]
                  - f_13 * gg_92[k]
                  + f_4 * gg_97[k]
                  + f_14 * gg_122[k]
                  - f_12 * gg_127[k];
        g_65[k] = g_25[k];
    }

#pragma omp simd aligned(gg_15, gg_18, gg_25, gg_90, gg_93, gg_100, gg_120, gg_123, \
                         gg_130 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_44 * gg_15[k]
                  + f_16 * gg_18[k]
                  - f_44 * gg_25[k]
                  - f_44 * gg_90[k]
                  + f_16 * gg_93[k]
                  - f_44 * gg_100[k]
                  + f_16 * gg_120[k]
                  - f_45 * gg_123[k]
                  + f_16 * gg_130[k];
        g_74[k] = g_26[k];
    }

#pragma omp simd aligned(gg_64, gg_71, gg_73, gg_169, gg_176, gg_178, gg_199, gg_206, \
                         gg_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = 5.625 * gg_64[k]
                  + 5.625 * gg_71[k]
                  - 7.5 * gg_73[k]
                  + 5.625 * gg_169[k]
                  + 5.625 * gg_176[k]
                  - 7.5 * gg_178[k]
                  - 7.5 * gg_199[k]
                  - 7.5 * gg_206[k]
                  + 10.0 * gg_208[k];
    }

#pragma omp simd aligned(gg_60, gg_63, gg_65, gg_70, gg_72, gg_74, gg_165, gg_168, gg_170, \
                         gg_175, gg_177, gg_179, gg_195, gg_198, gg_200, gg_205, gg_207, \
                         gg_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_46 * gg_60[k]
                  - f_47 * gg_63[k]
                  + f_48 * gg_65[k]
                  - f_46 * gg_70[k]
                  + f_48 * gg_72[k]
                  - f_49 * gg_74[k]
                  - f_46 * gg_165[k]
                  - f_47 * gg_168[k]
                  + f_48 * gg_170[k]
                  - f_46 * gg_175[k]
                  + f_48 * gg_177[k]
                  - f_49 * gg_179[k]
                  + f_50 * gg_195[k]
                  + f_49 * gg_198[k]
                  - f_51 * gg_200[k]
                  + f_50 * gg_205[k]
                  - f_51 * gg_207[k]
                  + f_52 * gg_209[k];
        g_39[k] = g_31[k];
    }

#pragma omp simd aligned(gg_62, gg_67, gg_69, gg_167, gg_172, gg_174, gg_197, gg_202, \
                         gg_204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = 5.625 * gg_62[k]
                  + 5.625 * gg_67[k]
                  - 7.5 * gg_69[k]
                  + 5.625 * gg_167[k]
                  + 5.625 * gg_172[k]
                  - 7.5 * gg_174[k]
                  - 7.5 * gg_197[k]
                  - 7.5 * gg_202[k]
                  + 10.0 * gg_204[k];
        g_48[k] = g_32[k];
    }

#pragma omp simd aligned(gg_60, gg_65, gg_70, gg_72, gg_165, gg_170, gg_175, gg_177, gg_195, \
                         gg_200, gg_205, gg_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_53 * gg_60[k]
                  - f_54 * gg_65[k]
                  - f_53 * gg_70[k]
                  + f_54 * gg_72[k]
                  + f_53 * gg_165[k]
                  - f_54 * gg_170[k]
                  - f_53 * gg_175[k]
                  + f_54 * gg_177[k]
                  - f_55 * gg_195[k]
                  + f_56 * gg_200[k]
                  + f_55 * gg_205[k]
                  - f_56 * gg_207[k];
        g_57[k] = g_33[k];
    }

#pragma omp simd aligned(gg_62, gg_67, gg_167, gg_172, gg_197, gg_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_16 * gg_62[k]
                  + f_15 * gg_67[k]
                  - f_16 * gg_167[k]
                  + f_15 * gg_172[k]
                  + f_17 * gg_197[k]
                  - f_3 * gg_202[k];
        g_66[k] = g_34[k];
    }

#pragma omp simd aligned(gg_60, gg_63, gg_70, gg_165, gg_168, gg_175, gg_195, gg_198, \
                         gg_205 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_57 * gg_60[k]
                  + f_58 * gg_63[k]
                  - f_57 * gg_70[k]
                  - f_57 * gg_165[k]
                  + f_58 * gg_168[k]
                  - f_57 * gg_175[k]
                  + f_13 * gg_195[k]
                  - f_14 * gg_198[k]
                  + f_13 * gg_205[k];
        g_75[k] = g_35[k];
    }

#pragma omp simd aligned(gg_0, gg_3, gg_5, gg_10, gg_12, gg_14, gg_45, gg_48, gg_50, gg_55, \
                         gg_57, gg_59, gg_75, gg_78, gg_80, gg_85, gg_87, gg_89, gg_150, \
                         gg_153, gg_155, gg_160, gg_162, gg_164, gg_180, gg_183, gg_185, \
                         gg_190, gg_192, gg_194, gg_210, gg_213, gg_215, gg_220, gg_222, \
                         gg_224 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = 0.140625 * gg_0[k]
                  + 0.28125 * gg_3[k]
                  - 1.125 * gg_5[k]
                  + 0.140625 * gg_10[k]
                  - 1.125 * gg_12[k]
                  + 0.375 * gg_14[k]
                  + 0.28125 * gg_45[k]
                  + 0.5625 * gg_48[k]
                  - 2.25 * gg_50[k]
                  + 0.28125 * gg_55[k]
                  - 2.25 * gg_57[k]
                  + 0.75 * gg_59[k]
                  - 1.125 * gg_75[k]
                  - 2.25 * gg_78[k]
                  + 9.0 * gg_80[k]
                  - 1.125 * gg_85[k]
                  + 9.0 * gg_87[k]
                  - 3.0 * gg_89[k]
                  + 0.140625 * gg_150[k]
                  + 0.28125 * gg_153[k]
                  - 1.125 * gg_155[k]
                  + 0.140625 * gg_160[k]
                  - 1.125 * gg_162[k]
                  + 0.375 * gg_164[k]
                  - 1.125 * gg_180[k]
                  - 2.25 * gg_183[k]
                  + 9.0 * gg_185[k]
                  - 1.125 * gg_190[k]
                  + 9.0 * gg_192[k]
                  - 3.0 * gg_194[k]
                  + 0.375 * gg_210[k]
                  + 0.75 * gg_213[k]
                  - 3.0 * gg_215[k]
                  + 0.375 * gg_220[k]
                  - 3.0 * gg_222[k]
                  + gg_224[k];
    }

#pragma omp simd aligned(gg_2, gg_7, gg_9, gg_47, gg_52, gg_54, gg_77, gg_82, gg_84, gg_152, \
                         gg_157, gg_159, gg_182, gg_187, gg_189, gg_212, gg_217, \
                         gg_219 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_46 * gg_2[k]
                  - f_46 * gg_7[k]
                  + f_50 * gg_9[k]
                  - f_47 * gg_47[k]
                  - f_47 * gg_52[k]
                  + f_49 * gg_54[k]
                  + f_48 * gg_77[k]
                  + f_48 * gg_82[k]
                  - f_51 * gg_84[k]
                  - f_46 * gg_152[k]
                  - f_46 * gg_157[k]
                  + f_50 * gg_159[k]
                  + f_48 * gg_182[k]
                  + f_48 * gg_187[k]
                  - f_51 * gg_189[k]
                  - f_49 * gg_212[k]
                  - f_49 * gg_217[k]
                  + f_52 * gg_219[k];
        g_49[k] = g_41[k];
    }

#pragma omp simd aligned(gg_0, gg_5, gg_10, gg_12, gg_45, gg_50, gg_55, gg_57, gg_75, gg_80, \
                         gg_85, gg_87, gg_150, gg_155, gg_160, gg_162, gg_180, gg_185, gg_190, \
                         gg_192, gg_210, gg_215, gg_220, gg_222 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_59 * gg_0[k]
                  + f_60 * gg_5[k]
                  + f_59 * gg_10[k]
                  - f_60 * gg_12[k]
                  - f_36 * gg_45[k]
                  + f_40 * gg_50[k]
                  + f_36 * gg_55[k]
                  - f_40 * gg_57[k]
                  + f_61 * gg_75[k]
                  - f_62 * gg_80[k]
                  - f_61 * gg_85[k]
                  + f_62 * gg_87[k]
                  - f_59 * gg_150[k]
                  + f_60 * gg_155[k]
                  + f_59 * gg_160[k]
                  - f_60 * gg_162[k]
                  + f_61 * gg_180[k]
                  - f_62 * gg_185[k]
                  - f_61 * gg_190[k]
                  + f_62 * gg_192[k]
                  - f_63 * gg_210[k]
                  + f_38 * gg_215[k]
                  + f_63 * gg_220[k]
                  - f_38 * gg_222[k];
        g_58[k] = g_42[k];
    }

#pragma omp simd aligned(gg_2, gg_7, gg_47, gg_52, gg_77, gg_82, gg_152, gg_157, gg_182, \
                         gg_187, gg_212, gg_217 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_22 * gg_2[k]
                  - f_18 * gg_7[k]
                  + f_23 * gg_47[k]
                  - f_19 * gg_52[k]
                  - f_21 * gg_77[k]
                  + f_20 * gg_82[k]
                  + f_22 * gg_152[k]
                  - f_18 * gg_157[k]
                  - f_21 * gg_182[k]
                  + f_20 * gg_187[k]
                  + f_24 * gg_212[k]
                  - f_21 * gg_217[k];
        g_67[k] = g_43[k];
    }

#pragma omp simd aligned(gg_0, gg_3, gg_10, gg_45, gg_48, gg_55, gg_75, gg_78, gg_85, gg_150, \
                         gg_153, gg_160, gg_180, gg_183, gg_190, gg_210, gg_213, \
                         gg_220 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_64 * gg_0[k]
                  - f_65 * gg_3[k]
                  + f_64 * gg_10[k]
                  + f_66 * gg_45[k]
                  - f_67 * gg_48[k]
                  + f_66 * gg_55[k]
                  - f_7 * gg_75[k]
                  + f_68 * gg_78[k]
                  - f_7 * gg_85[k]
                  + f_64 * gg_150[k]
                  - f_65 * gg_153[k]
                  + f_64 * gg_160[k]
                  - f_7 * gg_180[k]
                  + f_68 * gg_183[k]
                  - f_7 * gg_190[k]
                  + f_69 * gg_210[k]
                  - f_70 * gg_213[k]
                  + f_69 * gg_220[k];
        g_76[k] = g_44[k];
    }

#pragma omp simd aligned(gg_32, gg_37, gg_39, gg_107, gg_112, gg_114, gg_137, gg_142, \
                         gg_144 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = 5.625 * gg_32[k]
                  + 5.625 * gg_37[k]
                  - 7.5 * gg_39[k]
                  + 5.625 * gg_107[k]
                  + 5.625 * gg_112[k]
                  - 7.5 * gg_114[k]
                  - 7.5 * gg_137[k]
                  - 7.5 * gg_142[k]
                  + 10.0 * gg_144[k];
    }

#pragma omp simd aligned(gg_30, gg_35, gg_40, gg_42, gg_105, gg_110, gg_115, gg_117, gg_135, \
                         gg_140, gg_145, gg_147 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_53 * gg_30[k]
                  - f_54 * gg_35[k]
                  - f_53 * gg_40[k]
                  + f_54 * gg_42[k]
                  + f_53 * gg_105[k]
                  - f_54 * gg_110[k]
                  - f_53 * gg_115[k]
                  + f_54 * gg_117[k]
                  - f_55 * gg_135[k]
                  + f_56 * gg_140[k]
                  + f_55 * gg_145[k]
                  - f_56 * gg_147[k];
        g_59[k] = g_51[k];
    }

#pragma omp simd aligned(gg_32, gg_37, gg_107, gg_112, gg_137, gg_142 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_16 * gg_32[k]
                  + f_15 * gg_37[k]
                  - f_16 * gg_107[k]
                  + f_15 * gg_112[k]
                  + f_17 * gg_137[k]
                  - f_3 * gg_142[k];
        g_68[k] = g_52[k];
    }

#pragma omp simd aligned(gg_30, gg_33, gg_40, gg_105, gg_108, gg_115, gg_135, gg_138, \
                         gg_145 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = -f_57 * gg_30[k]
                  + f_58 * gg_33[k]
                  - f_57 * gg_40[k]
                  - f_57 * gg_105[k]
                  + f_58 * gg_108[k]
                  - f_57 * gg_115[k]
                  + f_13 * gg_135[k]
                  - f_14 * gg_138[k]
                  + f_13 * gg_145[k];
        g_77[k] = g_53[k];
    }

#pragma omp simd aligned(gg_0, gg_5, gg_10, gg_12, gg_75, gg_80, gg_85, gg_87, gg_150, gg_155, \
                         gg_160, gg_162, gg_180, gg_185, gg_190, \
                         gg_192 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = 0.3125 * gg_0[k]
                  - 1.875 * gg_5[k]
                  - 0.3125 * gg_10[k]
                  + 1.875 * gg_12[k]
                  - 1.875 * gg_75[k]
                  + 11.25 * gg_80[k]
                  + 1.875 * gg_85[k]
                  - 11.25 * gg_87[k]
                  - 0.3125 * gg_150[k]
                  + 1.875 * gg_155[k]
                  + 0.3125 * gg_160[k]
                  - 1.875 * gg_162[k]
                  + 1.875 * gg_180[k]
                  - 11.25 * gg_185[k]
                  - 1.875 * gg_190[k]
                  + 11.25 * gg_192[k];
    }

#pragma omp simd aligned(gg_2, gg_7, gg_77, gg_82, gg_152, gg_157, gg_182, \
                         gg_187 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_27 * gg_2[k]
                  + f_25 * gg_7[k]
                  + f_4 * gg_77[k]
                  - f_26 * gg_82[k]
                  + f_27 * gg_152[k]
                  - f_25 * gg_157[k]
                  - f_4 * gg_182[k]
                  + f_26 * gg_187[k];
        g_69[k] = g_61[k];
    }

#pragma omp simd aligned(gg_0, gg_3, gg_10, gg_75, gg_78, gg_85, gg_150, gg_153, gg_160, \
                         gg_180, gg_183, gg_190 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_71 * gg_0[k]
                  + f_72 * gg_3[k]
                  - f_71 * gg_10[k]
                  + f_72 * gg_75[k]
                  - f_15 * gg_78[k]
                  + f_72 * gg_85[k]
                  + f_71 * gg_150[k]
                  - f_72 * gg_153[k]
                  + f_71 * gg_160[k]
                  - f_72 * gg_180[k]
                  + f_15 * gg_183[k]
                  - f_72 * gg_190[k];
        g_78[k] = g_62[k];
    }

#pragma omp simd aligned(gg_30, gg_32, gg_33, gg_37, gg_40, gg_105, gg_107, gg_108, gg_112, \
                         gg_115 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = 4.375 * gg_32[k]
                  - 13.125 * gg_37[k]
                  - 13.125 * gg_107[k]
                  + 39.375 * gg_112[k];

        g_71[k] = f_30 * gg_30[k]
                  - f_31 * gg_33[k]
                  + f_30 * gg_40[k]
                  - f_28 * gg_105[k]
                  + f_29 * gg_108[k]
                  - f_28 * gg_115[k];
        g_79[k] = g_71[k];
    }

#pragma omp simd aligned(gg_0, gg_3, gg_10, gg_45, gg_48, gg_55, gg_150, gg_153, \
                         gg_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = 0.546875 * gg_0[k]
                  - 3.28125 * gg_3[k]
                  + 0.546875 * gg_10[k]
                  - 3.28125 * gg_45[k]
                  + 19.6875 * gg_48[k]
                  - 3.28125 * gg_55[k]
                  + 0.546875 * gg_150[k]
                  - 3.28125 * gg_153[k]
                  + 0.546875 * gg_160[k];
    }
}

}  // namespace simdtrf
