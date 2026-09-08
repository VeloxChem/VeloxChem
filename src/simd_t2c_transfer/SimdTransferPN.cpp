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


#include "SimdTransferPN.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_pn(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t sn, const size_t so, const size_t nmax) -> void
{
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *sn_0 = buffer.data(sn + 0);
    const auto *sn_1 = buffer.data(sn + 1);
    const auto *sn_2 = buffer.data(sn + 2);
    const auto *sn_3 = buffer.data(sn + 3);
    const auto *sn_4 = buffer.data(sn + 4);
    const auto *sn_5 = buffer.data(sn + 5);
    const auto *sn_6 = buffer.data(sn + 6);
    const auto *sn_7 = buffer.data(sn + 7);
    const auto *sn_8 = buffer.data(sn + 8);
    const auto *sn_9 = buffer.data(sn + 9);
    const auto *sn_10 = buffer.data(sn + 10);
    const auto *sn_11 = buffer.data(sn + 11);
    const auto *sn_12 = buffer.data(sn + 12);
    const auto *sn_13 = buffer.data(sn + 13);
    const auto *sn_14 = buffer.data(sn + 14);
    const auto *sn_15 = buffer.data(sn + 15);
    const auto *sn_16 = buffer.data(sn + 16);
    const auto *sn_17 = buffer.data(sn + 17);
    const auto *sn_18 = buffer.data(sn + 18);
    const auto *sn_19 = buffer.data(sn + 19);
    const auto *sn_20 = buffer.data(sn + 20);
    const auto *sn_21 = buffer.data(sn + 21);
    const auto *sn_22 = buffer.data(sn + 22);
    const auto *sn_23 = buffer.data(sn + 23);
    const auto *sn_24 = buffer.data(sn + 24);
    const auto *sn_25 = buffer.data(sn + 25);
    const auto *sn_26 = buffer.data(sn + 26);
    const auto *sn_27 = buffer.data(sn + 27);
    const auto *sn_28 = buffer.data(sn + 28);
    const auto *sn_29 = buffer.data(sn + 29);
    const auto *sn_30 = buffer.data(sn + 30);
    const auto *sn_31 = buffer.data(sn + 31);
    const auto *sn_32 = buffer.data(sn + 32);
    const auto *sn_33 = buffer.data(sn + 33);
    const auto *sn_34 = buffer.data(sn + 34);
    const auto *sn_35 = buffer.data(sn + 35);
    const auto *sn_36 = buffer.data(sn + 36);
    const auto *sn_37 = buffer.data(sn + 37);
    const auto *sn_38 = buffer.data(sn + 38);
    const auto *sn_39 = buffer.data(sn + 39);
    const auto *sn_40 = buffer.data(sn + 40);
    const auto *sn_41 = buffer.data(sn + 41);
    const auto *sn_42 = buffer.data(sn + 42);
    const auto *sn_43 = buffer.data(sn + 43);
    const auto *sn_44 = buffer.data(sn + 44);
    const auto *sn_45 = buffer.data(sn + 45);
    const auto *sn_46 = buffer.data(sn + 46);
    const auto *sn_47 = buffer.data(sn + 47);
    const auto *sn_48 = buffer.data(sn + 48);
    const auto *sn_49 = buffer.data(sn + 49);
    const auto *sn_50 = buffer.data(sn + 50);
    const auto *sn_51 = buffer.data(sn + 51);
    const auto *sn_52 = buffer.data(sn + 52);
    const auto *sn_53 = buffer.data(sn + 53);
    const auto *sn_54 = buffer.data(sn + 54);
    const auto *sn_55 = buffer.data(sn + 55);
    const auto *sn_56 = buffer.data(sn + 56);
    const auto *sn_57 = buffer.data(sn + 57);
    const auto *sn_58 = buffer.data(sn + 58);
    const auto *sn_59 = buffer.data(sn + 59);
    const auto *sn_60 = buffer.data(sn + 60);
    const auto *sn_61 = buffer.data(sn + 61);
    const auto *sn_62 = buffer.data(sn + 62);
    const auto *sn_63 = buffer.data(sn + 63);
    const auto *sn_64 = buffer.data(sn + 64);
    const auto *sn_65 = buffer.data(sn + 65);

    const auto *so_0 = buffer.data(so + 0);
    const auto *so_1 = buffer.data(so + 1);
    const auto *so_2 = buffer.data(so + 2);
    const auto *so_3 = buffer.data(so + 3);
    const auto *so_4 = buffer.data(so + 4);
    const auto *so_5 = buffer.data(so + 5);
    const auto *so_6 = buffer.data(so + 6);
    const auto *so_7 = buffer.data(so + 7);
    const auto *so_8 = buffer.data(so + 8);
    const auto *so_9 = buffer.data(so + 9);
    const auto *so_10 = buffer.data(so + 10);
    const auto *so_11 = buffer.data(so + 11);
    const auto *so_12 = buffer.data(so + 12);
    const auto *so_13 = buffer.data(so + 13);
    const auto *so_14 = buffer.data(so + 14);
    const auto *so_15 = buffer.data(so + 15);
    const auto *so_16 = buffer.data(so + 16);
    const auto *so_17 = buffer.data(so + 17);
    const auto *so_18 = buffer.data(so + 18);
    const auto *so_19 = buffer.data(so + 19);
    const auto *so_20 = buffer.data(so + 20);
    const auto *so_21 = buffer.data(so + 21);
    const auto *so_22 = buffer.data(so + 22);
    const auto *so_23 = buffer.data(so + 23);
    const auto *so_24 = buffer.data(so + 24);
    const auto *so_25 = buffer.data(so + 25);
    const auto *so_26 = buffer.data(so + 26);
    const auto *so_27 = buffer.data(so + 27);
    const auto *so_28 = buffer.data(so + 28);
    const auto *so_29 = buffer.data(so + 29);
    const auto *so_30 = buffer.data(so + 30);
    const auto *so_31 = buffer.data(so + 31);
    const auto *so_32 = buffer.data(so + 32);
    const auto *so_33 = buffer.data(so + 33);
    const auto *so_34 = buffer.data(so + 34);
    const auto *so_35 = buffer.data(so + 35);
    const auto *so_36 = buffer.data(so + 36);
    const auto *so_37 = buffer.data(so + 37);
    const auto *so_38 = buffer.data(so + 38);
    const auto *so_39 = buffer.data(so + 39);
    const auto *so_40 = buffer.data(so + 40);
    const auto *so_41 = buffer.data(so + 41);
    const auto *so_42 = buffer.data(so + 42);
    const auto *so_43 = buffer.data(so + 43);
    const auto *so_44 = buffer.data(so + 44);
    const auto *so_45 = buffer.data(so + 45);
    const auto *so_46 = buffer.data(so + 46);
    const auto *so_47 = buffer.data(so + 47);
    const auto *so_48 = buffer.data(so + 48);
    const auto *so_49 = buffer.data(so + 49);
    const auto *so_50 = buffer.data(so + 50);
    const auto *so_51 = buffer.data(so + 51);
    const auto *so_52 = buffer.data(so + 52);
    const auto *so_53 = buffer.data(so + 53);
    const auto *so_54 = buffer.data(so + 54);
    const auto *so_55 = buffer.data(so + 55);
    const auto *so_56 = buffer.data(so + 56);
    const auto *so_57 = buffer.data(so + 57);
    const auto *so_58 = buffer.data(so + 58);
    const auto *so_59 = buffer.data(so + 59);
    const auto *so_60 = buffer.data(so + 60);
    const auto *so_61 = buffer.data(so + 61);
    const auto *so_62 = buffer.data(so + 62);
    const auto *so_63 = buffer.data(so + 63);
    const auto *so_64 = buffer.data(so + 64);
    const auto *so_65 = buffer.data(so + 65);
    const auto *so_66 = buffer.data(so + 66);
    const auto *so_67 = buffer.data(so + 67);
    const auto *so_68 = buffer.data(so + 68);
    const auto *so_69 = buffer.data(so + 69);
    const auto *so_70 = buffer.data(so + 70);
    const auto *so_71 = buffer.data(so + 71);
    const auto *so_72 = buffer.data(so + 72);
    const auto *so_73 = buffer.data(so + 73);
    const auto *so_74 = buffer.data(so + 74);
    const auto *so_75 = buffer.data(so + 75);
    const auto *so_76 = buffer.data(so + 76);
    const auto *so_77 = buffer.data(so + 77);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, sn_0, sn_1, sn_2, sn_3, sn_4, so_0, \
                         so_1, so_2, so_3, so_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * sn_0[k]
                 + so_0[k];

        t_1[k] = -ab_x[k] * sn_1[k]
                 + so_1[k];

        t_2[k] = -ab_x[k] * sn_2[k]
                 + so_2[k];

        t_3[k] = -ab_x[k] * sn_3[k]
                 + so_3[k];

        t_4[k] = -ab_x[k] * sn_4[k]
                 + so_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, sn_5, sn_6, sn_7, sn_8, sn_9, so_5, \
                         so_6, so_7, so_8, so_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * sn_5[k]
                 + so_5[k];

        t_6[k] = -ab_x[k] * sn_6[k]
                 + so_6[k];

        t_7[k] = -ab_x[k] * sn_7[k]
                 + so_7[k];

        t_8[k] = -ab_x[k] * sn_8[k]
                 + so_8[k];

        t_9[k] = -ab_x[k] * sn_9[k]
                 + so_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, sn_10, sn_11, sn_12, sn_13, \
                         sn_14, so_10, so_11, so_12, so_13, so_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * sn_10[k]
                  + so_10[k];

        t_11[k] = -ab_x[k] * sn_11[k]
                  + so_11[k];

        t_12[k] = -ab_x[k] * sn_12[k]
                  + so_12[k];

        t_13[k] = -ab_x[k] * sn_13[k]
                  + so_13[k];

        t_14[k] = -ab_x[k] * sn_14[k]
                  + so_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, sn_15, sn_16, sn_17, sn_18, \
                         sn_19, so_15, so_16, so_17, so_18, so_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * sn_15[k]
                  + so_15[k];

        t_16[k] = -ab_x[k] * sn_16[k]
                  + so_16[k];

        t_17[k] = -ab_x[k] * sn_17[k]
                  + so_17[k];

        t_18[k] = -ab_x[k] * sn_18[k]
                  + so_18[k];

        t_19[k] = -ab_x[k] * sn_19[k]
                  + so_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, sn_20, sn_21, sn_22, sn_23, \
                         sn_24, so_20, so_21, so_22, so_23, so_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * sn_20[k]
                  + so_20[k];

        t_21[k] = -ab_x[k] * sn_21[k]
                  + so_21[k];

        t_22[k] = -ab_x[k] * sn_22[k]
                  + so_22[k];

        t_23[k] = -ab_x[k] * sn_23[k]
                  + so_23[k];

        t_24[k] = -ab_x[k] * sn_24[k]
                  + so_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, sn_25, sn_26, sn_27, sn_28, \
                         sn_29, so_25, so_26, so_27, so_28, so_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * sn_25[k]
                  + so_25[k];

        t_26[k] = -ab_x[k] * sn_26[k]
                  + so_26[k];

        t_27[k] = -ab_x[k] * sn_27[k]
                  + so_27[k];

        t_28[k] = -ab_x[k] * sn_28[k]
                  + so_28[k];

        t_29[k] = -ab_x[k] * sn_29[k]
                  + so_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, sn_30, sn_31, sn_32, sn_33, \
                         sn_34, so_30, so_31, so_32, so_33, so_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_x[k] * sn_30[k]
                  + so_30[k];

        t_31[k] = -ab_x[k] * sn_31[k]
                  + so_31[k];

        t_32[k] = -ab_x[k] * sn_32[k]
                  + so_32[k];

        t_33[k] = -ab_x[k] * sn_33[k]
                  + so_33[k];

        t_34[k] = -ab_x[k] * sn_34[k]
                  + so_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, sn_35, sn_36, sn_37, sn_38, \
                         sn_39, so_35, so_36, so_37, so_38, so_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_x[k] * sn_35[k]
                  + so_35[k];

        t_36[k] = -ab_x[k] * sn_36[k]
                  + so_36[k];

        t_37[k] = -ab_x[k] * sn_37[k]
                  + so_37[k];

        t_38[k] = -ab_x[k] * sn_38[k]
                  + so_38[k];

        t_39[k] = -ab_x[k] * sn_39[k]
                  + so_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, sn_40, sn_41, sn_42, sn_43, \
                         sn_44, so_40, so_41, so_42, so_43, so_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_x[k] * sn_40[k]
                  + so_40[k];

        t_41[k] = -ab_x[k] * sn_41[k]
                  + so_41[k];

        t_42[k] = -ab_x[k] * sn_42[k]
                  + so_42[k];

        t_43[k] = -ab_x[k] * sn_43[k]
                  + so_43[k];

        t_44[k] = -ab_x[k] * sn_44[k]
                  + so_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, sn_45, sn_46, sn_47, sn_48, \
                         sn_49, so_45, so_46, so_47, so_48, so_49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = -ab_x[k] * sn_45[k]
                  + so_45[k];

        t_46[k] = -ab_x[k] * sn_46[k]
                  + so_46[k];

        t_47[k] = -ab_x[k] * sn_47[k]
                  + so_47[k];

        t_48[k] = -ab_x[k] * sn_48[k]
                  + so_48[k];

        t_49[k] = -ab_x[k] * sn_49[k]
                  + so_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, sn_50, sn_51, sn_52, sn_53, \
                         sn_54, so_50, so_51, so_52, so_53, so_54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = -ab_x[k] * sn_50[k]
                  + so_50[k];

        t_51[k] = -ab_x[k] * sn_51[k]
                  + so_51[k];

        t_52[k] = -ab_x[k] * sn_52[k]
                  + so_52[k];

        t_53[k] = -ab_x[k] * sn_53[k]
                  + so_53[k];

        t_54[k] = -ab_x[k] * sn_54[k]
                  + so_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, sn_55, sn_56, sn_57, sn_58, \
                         sn_59, so_55, so_56, so_57, so_58, so_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = -ab_x[k] * sn_55[k]
                  + so_55[k];

        t_56[k] = -ab_x[k] * sn_56[k]
                  + so_56[k];

        t_57[k] = -ab_x[k] * sn_57[k]
                  + so_57[k];

        t_58[k] = -ab_x[k] * sn_58[k]
                  + so_58[k];

        t_59[k] = -ab_x[k] * sn_59[k]
                  + so_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, sn_60, sn_61, sn_62, sn_63, \
                         sn_64, so_60, so_61, so_62, so_63, so_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = -ab_x[k] * sn_60[k]
                  + so_60[k];

        t_61[k] = -ab_x[k] * sn_61[k]
                  + so_61[k];

        t_62[k] = -ab_x[k] * sn_62[k]
                  + so_62[k];

        t_63[k] = -ab_x[k] * sn_63[k]
                  + so_63[k];

        t_64[k] = -ab_x[k] * sn_64[k]
                  + so_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, ab_x, ab_y, sn_0, sn_1, sn_2, sn_65, so_1, \
                         so_3, so_4, so_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = -ab_x[k] * sn_65[k]
                  + so_65[k];

        t_66[k] = -ab_y[k] * sn_0[k]
                  + so_1[k];

        t_67[k] = -ab_y[k] * sn_1[k]
                  + so_3[k];

        t_68[k] = -ab_y[k] * sn_2[k]
                  + so_4[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, ab_y, sn_3, sn_4, sn_5, sn_6, sn_7, \
                         so_6, so_7, so_8, so_10, so_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_69[k] = -ab_y[k] * sn_3[k]
                  + so_6[k];

        t_70[k] = -ab_y[k] * sn_4[k]
                  + so_7[k];

        t_71[k] = -ab_y[k] * sn_5[k]
                  + so_8[k];

        t_72[k] = -ab_y[k] * sn_6[k]
                  + so_10[k];

        t_73[k] = -ab_y[k] * sn_7[k]
                  + so_11[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, ab_y, sn_8, sn_9, sn_10, sn_11, sn_12, \
                         so_12, so_13, so_15, so_16, so_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_74[k] = -ab_y[k] * sn_8[k]
                  + so_12[k];

        t_75[k] = -ab_y[k] * sn_9[k]
                  + so_13[k];

        t_76[k] = -ab_y[k] * sn_10[k]
                  + so_15[k];

        t_77[k] = -ab_y[k] * sn_11[k]
                  + so_16[k];

        t_78[k] = -ab_y[k] * sn_12[k]
                  + so_17[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, ab_y, sn_13, sn_14, sn_15, sn_16, \
                         sn_17, so_18, so_19, so_21, so_22, so_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_79[k] = -ab_y[k] * sn_13[k]
                  + so_18[k];

        t_80[k] = -ab_y[k] * sn_14[k]
                  + so_19[k];

        t_81[k] = -ab_y[k] * sn_15[k]
                  + so_21[k];

        t_82[k] = -ab_y[k] * sn_16[k]
                  + so_22[k];

        t_83[k] = -ab_y[k] * sn_17[k]
                  + so_23[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, ab_y, sn_18, sn_19, sn_20, sn_21, \
                         sn_22, so_24, so_25, so_26, so_28, so_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_84[k] = -ab_y[k] * sn_18[k]
                  + so_24[k];

        t_85[k] = -ab_y[k] * sn_19[k]
                  + so_25[k];

        t_86[k] = -ab_y[k] * sn_20[k]
                  + so_26[k];

        t_87[k] = -ab_y[k] * sn_21[k]
                  + so_28[k];

        t_88[k] = -ab_y[k] * sn_22[k]
                  + so_29[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, ab_y, sn_23, sn_24, sn_25, sn_26, \
                         sn_27, so_30, so_31, so_32, so_33, so_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_89[k] = -ab_y[k] * sn_23[k]
                  + so_30[k];

        t_90[k] = -ab_y[k] * sn_24[k]
                  + so_31[k];

        t_91[k] = -ab_y[k] * sn_25[k]
                  + so_32[k];

        t_92[k] = -ab_y[k] * sn_26[k]
                  + so_33[k];

        t_93[k] = -ab_y[k] * sn_27[k]
                  + so_34[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, ab_y, sn_28, sn_29, sn_30, sn_31, \
                         sn_32, so_36, so_37, so_38, so_39, so_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_94[k] = -ab_y[k] * sn_28[k]
                  + so_36[k];

        t_95[k] = -ab_y[k] * sn_29[k]
                  + so_37[k];

        t_96[k] = -ab_y[k] * sn_30[k]
                  + so_38[k];

        t_97[k] = -ab_y[k] * sn_31[k]
                  + so_39[k];

        t_98[k] = -ab_y[k] * sn_32[k]
                  + so_40[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, ab_y, sn_33, sn_34, sn_35, sn_36, \
                         sn_37, so_41, so_42, so_43, so_45, so_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_99[k] = -ab_y[k] * sn_33[k]
                  + so_41[k];

        t_100[k] = -ab_y[k] * sn_34[k]
                   + so_42[k];

        t_101[k] = -ab_y[k] * sn_35[k]
                   + so_43[k];

        t_102[k] = -ab_y[k] * sn_36[k]
                   + so_45[k];

        t_103[k] = -ab_y[k] * sn_37[k]
                   + so_46[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, ab_y, sn_38, sn_39, sn_40, sn_41, \
                         sn_42, so_47, so_48, so_49, so_50, so_51 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_104[k] = -ab_y[k] * sn_38[k]
                   + so_47[k];

        t_105[k] = -ab_y[k] * sn_39[k]
                   + so_48[k];

        t_106[k] = -ab_y[k] * sn_40[k]
                   + so_49[k];

        t_107[k] = -ab_y[k] * sn_41[k]
                   + so_50[k];

        t_108[k] = -ab_y[k] * sn_42[k]
                   + so_51[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, ab_y, sn_43, sn_44, sn_45, sn_46, \
                         sn_47, so_52, so_53, so_55, so_56, so_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_109[k] = -ab_y[k] * sn_43[k]
                   + so_52[k];

        t_110[k] = -ab_y[k] * sn_44[k]
                   + so_53[k];

        t_111[k] = -ab_y[k] * sn_45[k]
                   + so_55[k];

        t_112[k] = -ab_y[k] * sn_46[k]
                   + so_56[k];

        t_113[k] = -ab_y[k] * sn_47[k]
                   + so_57[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, ab_y, sn_48, sn_49, sn_50, sn_51, \
                         sn_52, so_58, so_59, so_60, so_61, so_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_114[k] = -ab_y[k] * sn_48[k]
                   + so_58[k];

        t_115[k] = -ab_y[k] * sn_49[k]
                   + so_59[k];

        t_116[k] = -ab_y[k] * sn_50[k]
                   + so_60[k];

        t_117[k] = -ab_y[k] * sn_51[k]
                   + so_61[k];

        t_118[k] = -ab_y[k] * sn_52[k]
                   + so_62[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, ab_y, sn_53, sn_54, sn_55, sn_56, \
                         sn_57, so_63, so_64, so_66, so_67, so_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_119[k] = -ab_y[k] * sn_53[k]
                   + so_63[k];

        t_120[k] = -ab_y[k] * sn_54[k]
                   + so_64[k];

        t_121[k] = -ab_y[k] * sn_55[k]
                   + so_66[k];

        t_122[k] = -ab_y[k] * sn_56[k]
                   + so_67[k];

        t_123[k] = -ab_y[k] * sn_57[k]
                   + so_68[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, ab_y, sn_58, sn_59, sn_60, sn_61, \
                         sn_62, so_69, so_70, so_71, so_72, so_73 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_124[k] = -ab_y[k] * sn_58[k]
                   + so_69[k];

        t_125[k] = -ab_y[k] * sn_59[k]
                   + so_70[k];

        t_126[k] = -ab_y[k] * sn_60[k]
                   + so_71[k];

        t_127[k] = -ab_y[k] * sn_61[k]
                   + so_72[k];

        t_128[k] = -ab_y[k] * sn_62[k]
                   + so_73[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, ab_y, ab_z, sn_0, sn_63, sn_64, sn_65, \
                         so_2, so_74, so_75, so_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_129[k] = -ab_y[k] * sn_63[k]
                   + so_74[k];

        t_130[k] = -ab_y[k] * sn_64[k]
                   + so_75[k];

        t_131[k] = -ab_y[k] * sn_65[k]
                   + so_76[k];

        t_132[k] = -ab_z[k] * sn_0[k]
                   + so_2[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, ab_z, sn_1, sn_2, sn_3, sn_4, \
                         sn_5, so_4, so_5, so_7, so_8, so_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_133[k] = -ab_z[k] * sn_1[k]
                   + so_4[k];

        t_134[k] = -ab_z[k] * sn_2[k]
                   + so_5[k];

        t_135[k] = -ab_z[k] * sn_3[k]
                   + so_7[k];

        t_136[k] = -ab_z[k] * sn_4[k]
                   + so_8[k];

        t_137[k] = -ab_z[k] * sn_5[k]
                   + so_9[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, ab_z, sn_6, sn_7, sn_8, sn_9, \
                         sn_10, so_11, so_12, so_13, so_14, so_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_138[k] = -ab_z[k] * sn_6[k]
                   + so_11[k];

        t_139[k] = -ab_z[k] * sn_7[k]
                   + so_12[k];

        t_140[k] = -ab_z[k] * sn_8[k]
                   + so_13[k];

        t_141[k] = -ab_z[k] * sn_9[k]
                   + so_14[k];

        t_142[k] = -ab_z[k] * sn_10[k]
                   + so_16[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, t_147, ab_z, sn_11, sn_12, sn_13, sn_14, \
                         sn_15, so_17, so_18, so_19, so_20, so_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_143[k] = -ab_z[k] * sn_11[k]
                   + so_17[k];

        t_144[k] = -ab_z[k] * sn_12[k]
                   + so_18[k];

        t_145[k] = -ab_z[k] * sn_13[k]
                   + so_19[k];

        t_146[k] = -ab_z[k] * sn_14[k]
                   + so_20[k];

        t_147[k] = -ab_z[k] * sn_15[k]
                   + so_22[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, ab_z, sn_16, sn_17, sn_18, sn_19, \
                         sn_20, so_23, so_24, so_25, so_26, so_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_148[k] = -ab_z[k] * sn_16[k]
                   + so_23[k];

        t_149[k] = -ab_z[k] * sn_17[k]
                   + so_24[k];

        t_150[k] = -ab_z[k] * sn_18[k]
                   + so_25[k];

        t_151[k] = -ab_z[k] * sn_19[k]
                   + so_26[k];

        t_152[k] = -ab_z[k] * sn_20[k]
                   + so_27[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, ab_z, sn_21, sn_22, sn_23, sn_24, \
                         sn_25, so_29, so_30, so_31, so_32, so_33 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_153[k] = -ab_z[k] * sn_21[k]
                   + so_29[k];

        t_154[k] = -ab_z[k] * sn_22[k]
                   + so_30[k];

        t_155[k] = -ab_z[k] * sn_23[k]
                   + so_31[k];

        t_156[k] = -ab_z[k] * sn_24[k]
                   + so_32[k];

        t_157[k] = -ab_z[k] * sn_25[k]
                   + so_33[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, ab_z, sn_26, sn_27, sn_28, sn_29, \
                         sn_30, so_34, so_35, so_37, so_38, so_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_158[k] = -ab_z[k] * sn_26[k]
                   + so_34[k];

        t_159[k] = -ab_z[k] * sn_27[k]
                   + so_35[k];

        t_160[k] = -ab_z[k] * sn_28[k]
                   + so_37[k];

        t_161[k] = -ab_z[k] * sn_29[k]
                   + so_38[k];

        t_162[k] = -ab_z[k] * sn_30[k]
                   + so_39[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, ab_z, sn_31, sn_32, sn_33, sn_34, \
                         sn_35, so_40, so_41, so_42, so_43, so_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_163[k] = -ab_z[k] * sn_31[k]
                   + so_40[k];

        t_164[k] = -ab_z[k] * sn_32[k]
                   + so_41[k];

        t_165[k] = -ab_z[k] * sn_33[k]
                   + so_42[k];

        t_166[k] = -ab_z[k] * sn_34[k]
                   + so_43[k];

        t_167[k] = -ab_z[k] * sn_35[k]
                   + so_44[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, ab_z, sn_36, sn_37, sn_38, sn_39, \
                         sn_40, so_46, so_47, so_48, so_49, so_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_168[k] = -ab_z[k] * sn_36[k]
                   + so_46[k];

        t_169[k] = -ab_z[k] * sn_37[k]
                   + so_47[k];

        t_170[k] = -ab_z[k] * sn_38[k]
                   + so_48[k];

        t_171[k] = -ab_z[k] * sn_39[k]
                   + so_49[k];

        t_172[k] = -ab_z[k] * sn_40[k]
                   + so_50[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, ab_z, sn_41, sn_42, sn_43, sn_44, \
                         sn_45, so_51, so_52, so_53, so_54, so_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_173[k] = -ab_z[k] * sn_41[k]
                   + so_51[k];

        t_174[k] = -ab_z[k] * sn_42[k]
                   + so_52[k];

        t_175[k] = -ab_z[k] * sn_43[k]
                   + so_53[k];

        t_176[k] = -ab_z[k] * sn_44[k]
                   + so_54[k];

        t_177[k] = -ab_z[k] * sn_45[k]
                   + so_56[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, ab_z, sn_46, sn_47, sn_48, sn_49, \
                         sn_50, so_57, so_58, so_59, so_60, so_61 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_178[k] = -ab_z[k] * sn_46[k]
                   + so_57[k];

        t_179[k] = -ab_z[k] * sn_47[k]
                   + so_58[k];

        t_180[k] = -ab_z[k] * sn_48[k]
                   + so_59[k];

        t_181[k] = -ab_z[k] * sn_49[k]
                   + so_60[k];

        t_182[k] = -ab_z[k] * sn_50[k]
                   + so_61[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, ab_z, sn_51, sn_52, sn_53, sn_54, \
                         sn_55, so_62, so_63, so_64, so_65, so_67 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_183[k] = -ab_z[k] * sn_51[k]
                   + so_62[k];

        t_184[k] = -ab_z[k] * sn_52[k]
                   + so_63[k];

        t_185[k] = -ab_z[k] * sn_53[k]
                   + so_64[k];

        t_186[k] = -ab_z[k] * sn_54[k]
                   + so_65[k];

        t_187[k] = -ab_z[k] * sn_55[k]
                   + so_67[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, ab_z, sn_56, sn_57, sn_58, sn_59, \
                         sn_60, so_68, so_69, so_70, so_71, so_72 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_188[k] = -ab_z[k] * sn_56[k]
                   + so_68[k];

        t_189[k] = -ab_z[k] * sn_57[k]
                   + so_69[k];

        t_190[k] = -ab_z[k] * sn_58[k]
                   + so_70[k];

        t_191[k] = -ab_z[k] * sn_59[k]
                   + so_71[k];

        t_192[k] = -ab_z[k] * sn_60[k]
                   + so_72[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, ab_z, sn_61, sn_62, sn_63, sn_64, \
                         sn_65, so_73, so_74, so_75, so_76, so_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_193[k] = -ab_z[k] * sn_61[k]
                   + so_73[k];

        t_194[k] = -ab_z[k] * sn_62[k]
                   + so_74[k];

        t_195[k] = -ab_z[k] * sn_63[k]
                   + so_75[k];

        t_196[k] = -ab_z[k] * sn_64[k]
                   + so_76[k];

        t_197[k] = -ab_z[k] * sn_65[k]
                   + so_77[k];
    }
}

}  // namespace simdtrf
