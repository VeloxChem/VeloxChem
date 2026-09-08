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


#include "SimdTransferDI.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_di(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t pi, const size_t pk, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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
    const auto *pi_50 = buffer.data(pi + 50);
    const auto *pi_51 = buffer.data(pi + 51);
    const auto *pi_52 = buffer.data(pi + 52);
    const auto *pi_53 = buffer.data(pi + 53);
    const auto *pi_54 = buffer.data(pi + 54);
    const auto *pi_55 = buffer.data(pi + 55);
    const auto *pi_56 = buffer.data(pi + 56);
    const auto *pi_57 = buffer.data(pi + 57);
    const auto *pi_58 = buffer.data(pi + 58);
    const auto *pi_59 = buffer.data(pi + 59);
    const auto *pi_60 = buffer.data(pi + 60);
    const auto *pi_61 = buffer.data(pi + 61);
    const auto *pi_62 = buffer.data(pi + 62);
    const auto *pi_63 = buffer.data(pi + 63);
    const auto *pi_64 = buffer.data(pi + 64);
    const auto *pi_65 = buffer.data(pi + 65);
    const auto *pi_66 = buffer.data(pi + 66);
    const auto *pi_67 = buffer.data(pi + 67);
    const auto *pi_68 = buffer.data(pi + 68);
    const auto *pi_69 = buffer.data(pi + 69);
    const auto *pi_70 = buffer.data(pi + 70);
    const auto *pi_71 = buffer.data(pi + 71);
    const auto *pi_72 = buffer.data(pi + 72);
    const auto *pi_73 = buffer.data(pi + 73);
    const auto *pi_74 = buffer.data(pi + 74);
    const auto *pi_75 = buffer.data(pi + 75);
    const auto *pi_76 = buffer.data(pi + 76);
    const auto *pi_77 = buffer.data(pi + 77);
    const auto *pi_78 = buffer.data(pi + 78);
    const auto *pi_79 = buffer.data(pi + 79);
    const auto *pi_80 = buffer.data(pi + 80);
    const auto *pi_81 = buffer.data(pi + 81);
    const auto *pi_82 = buffer.data(pi + 82);
    const auto *pi_83 = buffer.data(pi + 83);

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
    const auto *pk_36 = buffer.data(pk + 36);
    const auto *pk_37 = buffer.data(pk + 37);
    const auto *pk_38 = buffer.data(pk + 38);
    const auto *pk_39 = buffer.data(pk + 39);
    const auto *pk_40 = buffer.data(pk + 40);
    const auto *pk_41 = buffer.data(pk + 41);
    const auto *pk_42 = buffer.data(pk + 42);
    const auto *pk_43 = buffer.data(pk + 43);
    const auto *pk_44 = buffer.data(pk + 44);
    const auto *pk_45 = buffer.data(pk + 45);
    const auto *pk_46 = buffer.data(pk + 46);
    const auto *pk_47 = buffer.data(pk + 47);
    const auto *pk_48 = buffer.data(pk + 48);
    const auto *pk_49 = buffer.data(pk + 49);
    const auto *pk_50 = buffer.data(pk + 50);
    const auto *pk_51 = buffer.data(pk + 51);
    const auto *pk_52 = buffer.data(pk + 52);
    const auto *pk_53 = buffer.data(pk + 53);
    const auto *pk_54 = buffer.data(pk + 54);
    const auto *pk_55 = buffer.data(pk + 55);
    const auto *pk_56 = buffer.data(pk + 56);
    const auto *pk_57 = buffer.data(pk + 57);
    const auto *pk_58 = buffer.data(pk + 58);
    const auto *pk_59 = buffer.data(pk + 59);
    const auto *pk_60 = buffer.data(pk + 60);
    const auto *pk_61 = buffer.data(pk + 61);
    const auto *pk_62 = buffer.data(pk + 62);
    const auto *pk_63 = buffer.data(pk + 63);
    const auto *pk_64 = buffer.data(pk + 64);
    const auto *pk_65 = buffer.data(pk + 65);
    const auto *pk_66 = buffer.data(pk + 66);
    const auto *pk_67 = buffer.data(pk + 67);
    const auto *pk_68 = buffer.data(pk + 68);
    const auto *pk_69 = buffer.data(pk + 69);
    const auto *pk_70 = buffer.data(pk + 70);
    const auto *pk_72 = buffer.data(pk + 72);
    const auto *pk_73 = buffer.data(pk + 73);
    const auto *pk_74 = buffer.data(pk + 74);
    const auto *pk_75 = buffer.data(pk + 75);
    const auto *pk_76 = buffer.data(pk + 76);
    const auto *pk_77 = buffer.data(pk + 77);
    const auto *pk_78 = buffer.data(pk + 78);
    const auto *pk_79 = buffer.data(pk + 79);
    const auto *pk_80 = buffer.data(pk + 80);
    const auto *pk_81 = buffer.data(pk + 81);
    const auto *pk_82 = buffer.data(pk + 82);
    const auto *pk_83 = buffer.data(pk + 83);
    const auto *pk_84 = buffer.data(pk + 84);
    const auto *pk_85 = buffer.data(pk + 85);
    const auto *pk_86 = buffer.data(pk + 86);
    const auto *pk_87 = buffer.data(pk + 87);
    const auto *pk_88 = buffer.data(pk + 88);
    const auto *pk_89 = buffer.data(pk + 89);
    const auto *pk_90 = buffer.data(pk + 90);
    const auto *pk_91 = buffer.data(pk + 91);
    const auto *pk_92 = buffer.data(pk + 92);
    const auto *pk_93 = buffer.data(pk + 93);
    const auto *pk_94 = buffer.data(pk + 94);
    const auto *pk_95 = buffer.data(pk + 95);
    const auto *pk_96 = buffer.data(pk + 96);
    const auto *pk_97 = buffer.data(pk + 97);
    const auto *pk_98 = buffer.data(pk + 98);
    const auto *pk_99 = buffer.data(pk + 99);
    const auto *pk_100 = buffer.data(pk + 100);
    const auto *pk_101 = buffer.data(pk + 101);
    const auto *pk_102 = buffer.data(pk + 102);
    const auto *pk_103 = buffer.data(pk + 103);
    const auto *pk_104 = buffer.data(pk + 104);
    const auto *pk_105 = buffer.data(pk + 105);
    const auto *pk_106 = buffer.data(pk + 106);
    const auto *pk_107 = buffer.data(pk + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, pi_0, pi_1, pi_2, pi_3, pi_4, pk_0, \
                         pk_1, pk_2, pk_3, pk_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * pi_0[k]
                 + pk_0[k];

        t_1[k] = -ab_x[k] * pi_1[k]
                 + pk_1[k];

        t_2[k] = -ab_x[k] * pi_2[k]
                 + pk_2[k];

        t_3[k] = -ab_x[k] * pi_3[k]
                 + pk_3[k];

        t_4[k] = -ab_x[k] * pi_4[k]
                 + pk_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, pi_5, pi_6, pi_7, pi_8, pi_9, pk_5, \
                         pk_6, pk_7, pk_8, pk_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * pi_5[k]
                 + pk_5[k];

        t_6[k] = -ab_x[k] * pi_6[k]
                 + pk_6[k];

        t_7[k] = -ab_x[k] * pi_7[k]
                 + pk_7[k];

        t_8[k] = -ab_x[k] * pi_8[k]
                 + pk_8[k];

        t_9[k] = -ab_x[k] * pi_9[k]
                 + pk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, pi_10, pi_11, pi_12, pi_13, \
                         pi_14, pk_10, pk_11, pk_12, pk_13, pk_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * pi_10[k]
                  + pk_10[k];

        t_11[k] = -ab_x[k] * pi_11[k]
                  + pk_11[k];

        t_12[k] = -ab_x[k] * pi_12[k]
                  + pk_12[k];

        t_13[k] = -ab_x[k] * pi_13[k]
                  + pk_13[k];

        t_14[k] = -ab_x[k] * pi_14[k]
                  + pk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, pi_15, pi_16, pi_17, pi_18, \
                         pi_19, pk_15, pk_16, pk_17, pk_18, pk_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * pi_15[k]
                  + pk_15[k];

        t_16[k] = -ab_x[k] * pi_16[k]
                  + pk_16[k];

        t_17[k] = -ab_x[k] * pi_17[k]
                  + pk_17[k];

        t_18[k] = -ab_x[k] * pi_18[k]
                  + pk_18[k];

        t_19[k] = -ab_x[k] * pi_19[k]
                  + pk_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, pi_20, pi_21, pi_22, pi_23, \
                         pi_24, pk_20, pk_21, pk_22, pk_23, pk_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * pi_20[k]
                  + pk_20[k];

        t_21[k] = -ab_x[k] * pi_21[k]
                  + pk_21[k];

        t_22[k] = -ab_x[k] * pi_22[k]
                  + pk_22[k];

        t_23[k] = -ab_x[k] * pi_23[k]
                  + pk_23[k];

        t_24[k] = -ab_x[k] * pi_24[k]
                  + pk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, pi_25, pi_26, pi_27, pi_28, \
                         pi_29, pk_25, pk_26, pk_27, pk_36, pk_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * pi_25[k]
                  + pk_25[k];

        t_26[k] = -ab_x[k] * pi_26[k]
                  + pk_26[k];

        t_27[k] = -ab_x[k] * pi_27[k]
                  + pk_27[k];

        t_28[k] = -ab_x[k] * pi_28[k]
                  + pk_36[k];

        t_29[k] = -ab_x[k] * pi_29[k]
                  + pk_37[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, pi_30, pi_31, pi_32, pi_33, \
                         pi_34, pk_38, pk_39, pk_40, pk_41, pk_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_x[k] * pi_30[k]
                  + pk_38[k];

        t_31[k] = -ab_x[k] * pi_31[k]
                  + pk_39[k];

        t_32[k] = -ab_x[k] * pi_32[k]
                  + pk_40[k];

        t_33[k] = -ab_x[k] * pi_33[k]
                  + pk_41[k];

        t_34[k] = -ab_x[k] * pi_34[k]
                  + pk_42[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, pi_35, pi_36, pi_37, pi_38, \
                         pi_39, pk_43, pk_44, pk_45, pk_46, pk_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_x[k] * pi_35[k]
                  + pk_43[k];

        t_36[k] = -ab_x[k] * pi_36[k]
                  + pk_44[k];

        t_37[k] = -ab_x[k] * pi_37[k]
                  + pk_45[k];

        t_38[k] = -ab_x[k] * pi_38[k]
                  + pk_46[k];

        t_39[k] = -ab_x[k] * pi_39[k]
                  + pk_47[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, pi_40, pi_41, pi_42, pi_43, \
                         pi_44, pk_48, pk_49, pk_50, pk_51, pk_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_x[k] * pi_40[k]
                  + pk_48[k];

        t_41[k] = -ab_x[k] * pi_41[k]
                  + pk_49[k];

        t_42[k] = -ab_x[k] * pi_42[k]
                  + pk_50[k];

        t_43[k] = -ab_x[k] * pi_43[k]
                  + pk_51[k];

        t_44[k] = -ab_x[k] * pi_44[k]
                  + pk_52[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, pi_45, pi_46, pi_47, pi_48, \
                         pi_49, pk_53, pk_54, pk_55, pk_56, pk_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = -ab_x[k] * pi_45[k]
                  + pk_53[k];

        t_46[k] = -ab_x[k] * pi_46[k]
                  + pk_54[k];

        t_47[k] = -ab_x[k] * pi_47[k]
                  + pk_55[k];

        t_48[k] = -ab_x[k] * pi_48[k]
                  + pk_56[k];

        t_49[k] = -ab_x[k] * pi_49[k]
                  + pk_57[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, pi_50, pi_51, pi_52, pi_53, \
                         pi_54, pk_58, pk_59, pk_60, pk_61, pk_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = -ab_x[k] * pi_50[k]
                  + pk_58[k];

        t_51[k] = -ab_x[k] * pi_51[k]
                  + pk_59[k];

        t_52[k] = -ab_x[k] * pi_52[k]
                  + pk_60[k];

        t_53[k] = -ab_x[k] * pi_53[k]
                  + pk_61[k];

        t_54[k] = -ab_x[k] * pi_54[k]
                  + pk_62[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, pi_55, pi_56, pi_57, pi_58, \
                         pi_59, pk_63, pk_72, pk_73, pk_74, pk_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = -ab_x[k] * pi_55[k]
                  + pk_63[k];

        t_56[k] = -ab_x[k] * pi_56[k]
                  + pk_72[k];

        t_57[k] = -ab_x[k] * pi_57[k]
                  + pk_73[k];

        t_58[k] = -ab_x[k] * pi_58[k]
                  + pk_74[k];

        t_59[k] = -ab_x[k] * pi_59[k]
                  + pk_75[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, pi_60, pi_61, pi_62, pi_63, \
                         pi_64, pk_76, pk_77, pk_78, pk_79, pk_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = -ab_x[k] * pi_60[k]
                  + pk_76[k];

        t_61[k] = -ab_x[k] * pi_61[k]
                  + pk_77[k];

        t_62[k] = -ab_x[k] * pi_62[k]
                  + pk_78[k];

        t_63[k] = -ab_x[k] * pi_63[k]
                  + pk_79[k];

        t_64[k] = -ab_x[k] * pi_64[k]
                  + pk_80[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, pi_65, pi_66, pi_67, pi_68, \
                         pi_69, pk_81, pk_82, pk_83, pk_84, pk_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = -ab_x[k] * pi_65[k]
                  + pk_81[k];

        t_66[k] = -ab_x[k] * pi_66[k]
                  + pk_82[k];

        t_67[k] = -ab_x[k] * pi_67[k]
                  + pk_83[k];

        t_68[k] = -ab_x[k] * pi_68[k]
                  + pk_84[k];

        t_69[k] = -ab_x[k] * pi_69[k]
                  + pk_85[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, pi_70, pi_71, pi_72, pi_73, \
                         pi_74, pk_86, pk_87, pk_88, pk_89, pk_90 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = -ab_x[k] * pi_70[k]
                  + pk_86[k];

        t_71[k] = -ab_x[k] * pi_71[k]
                  + pk_87[k];

        t_72[k] = -ab_x[k] * pi_72[k]
                  + pk_88[k];

        t_73[k] = -ab_x[k] * pi_73[k]
                  + pk_89[k];

        t_74[k] = -ab_x[k] * pi_74[k]
                  + pk_90[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, pi_75, pi_76, pi_77, pi_78, \
                         pi_79, pk_91, pk_92, pk_93, pk_94, pk_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = -ab_x[k] * pi_75[k]
                  + pk_91[k];

        t_76[k] = -ab_x[k] * pi_76[k]
                  + pk_92[k];

        t_77[k] = -ab_x[k] * pi_77[k]
                  + pk_93[k];

        t_78[k] = -ab_x[k] * pi_78[k]
                  + pk_94[k];

        t_79[k] = -ab_x[k] * pi_79[k]
                  + pk_95[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, ab_x, pi_80, pi_81, pi_82, pi_83, pk_96, \
                         pk_97, pk_98, pk_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = -ab_x[k] * pi_80[k]
                  + pk_96[k];

        t_81[k] = -ab_x[k] * pi_81[k]
                  + pk_97[k];

        t_82[k] = -ab_x[k] * pi_82[k]
                  + pk_98[k];

        t_83[k] = -ab_x[k] * pi_83[k]
                  + pk_99[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, ab_y, pi_28, pi_29, pi_30, pi_31, \
                         pi_32, pk_37, pk_39, pk_40, pk_42, pk_43 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_84[k] = -ab_y[k] * pi_28[k]
                  + pk_37[k];

        t_85[k] = -ab_y[k] * pi_29[k]
                  + pk_39[k];

        t_86[k] = -ab_y[k] * pi_30[k]
                  + pk_40[k];

        t_87[k] = -ab_y[k] * pi_31[k]
                  + pk_42[k];

        t_88[k] = -ab_y[k] * pi_32[k]
                  + pk_43[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, ab_y, pi_33, pi_34, pi_35, pi_36, \
                         pi_37, pk_44, pk_46, pk_47, pk_48, pk_49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_89[k] = -ab_y[k] * pi_33[k]
                  + pk_44[k];

        t_90[k] = -ab_y[k] * pi_34[k]
                  + pk_46[k];

        t_91[k] = -ab_y[k] * pi_35[k]
                  + pk_47[k];

        t_92[k] = -ab_y[k] * pi_36[k]
                  + pk_48[k];

        t_93[k] = -ab_y[k] * pi_37[k]
                  + pk_49[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, ab_y, pi_38, pi_39, pi_40, pi_41, \
                         pi_42, pk_51, pk_52, pk_53, pk_54, pk_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_94[k] = -ab_y[k] * pi_38[k]
                  + pk_51[k];

        t_95[k] = -ab_y[k] * pi_39[k]
                  + pk_52[k];

        t_96[k] = -ab_y[k] * pi_40[k]
                  + pk_53[k];

        t_97[k] = -ab_y[k] * pi_41[k]
                  + pk_54[k];

        t_98[k] = -ab_y[k] * pi_42[k]
                  + pk_55[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, ab_y, pi_43, pi_44, pi_45, pi_46, \
                         pi_47, pk_57, pk_58, pk_59, pk_60, pk_61 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_99[k] = -ab_y[k] * pi_43[k]
                  + pk_57[k];

        t_100[k] = -ab_y[k] * pi_44[k]
                   + pk_58[k];

        t_101[k] = -ab_y[k] * pi_45[k]
                   + pk_59[k];

        t_102[k] = -ab_y[k] * pi_46[k]
                   + pk_60[k];

        t_103[k] = -ab_y[k] * pi_47[k]
                   + pk_61[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, ab_y, pi_48, pi_49, pi_50, pi_51, \
                         pi_52, pk_62, pk_64, pk_65, pk_66, pk_67 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_104[k] = -ab_y[k] * pi_48[k]
                   + pk_62[k];

        t_105[k] = -ab_y[k] * pi_49[k]
                   + pk_64[k];

        t_106[k] = -ab_y[k] * pi_50[k]
                   + pk_65[k];

        t_107[k] = -ab_y[k] * pi_51[k]
                   + pk_66[k];

        t_108[k] = -ab_y[k] * pi_52[k]
                   + pk_67[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, ab_y, pi_53, pi_54, pi_55, pi_56, \
                         pi_57, pk_68, pk_69, pk_70, pk_73, pk_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_109[k] = -ab_y[k] * pi_53[k]
                   + pk_68[k];

        t_110[k] = -ab_y[k] * pi_54[k]
                   + pk_69[k];

        t_111[k] = -ab_y[k] * pi_55[k]
                   + pk_70[k];

        t_112[k] = -ab_y[k] * pi_56[k]
                   + pk_73[k];

        t_113[k] = -ab_y[k] * pi_57[k]
                   + pk_75[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, ab_y, pi_58, pi_59, pi_60, pi_61, \
                         pi_62, pk_76, pk_78, pk_79, pk_80, pk_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_114[k] = -ab_y[k] * pi_58[k]
                   + pk_76[k];

        t_115[k] = -ab_y[k] * pi_59[k]
                   + pk_78[k];

        t_116[k] = -ab_y[k] * pi_60[k]
                   + pk_79[k];

        t_117[k] = -ab_y[k] * pi_61[k]
                   + pk_80[k];

        t_118[k] = -ab_y[k] * pi_62[k]
                   + pk_82[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, ab_y, pi_63, pi_64, pi_65, pi_66, \
                         pi_67, pk_83, pk_84, pk_85, pk_87, pk_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_119[k] = -ab_y[k] * pi_63[k]
                   + pk_83[k];

        t_120[k] = -ab_y[k] * pi_64[k]
                   + pk_84[k];

        t_121[k] = -ab_y[k] * pi_65[k]
                   + pk_85[k];

        t_122[k] = -ab_y[k] * pi_66[k]
                   + pk_87[k];

        t_123[k] = -ab_y[k] * pi_67[k]
                   + pk_88[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, ab_y, pi_68, pi_69, pi_70, pi_71, \
                         pi_72, pk_89, pk_90, pk_91, pk_93, pk_94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_124[k] = -ab_y[k] * pi_68[k]
                   + pk_89[k];

        t_125[k] = -ab_y[k] * pi_69[k]
                   + pk_90[k];

        t_126[k] = -ab_y[k] * pi_70[k]
                   + pk_91[k];

        t_127[k] = -ab_y[k] * pi_71[k]
                   + pk_93[k];

        t_128[k] = -ab_y[k] * pi_72[k]
                   + pk_94[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, ab_y, pi_73, pi_74, pi_75, pi_76, \
                         pi_77, pk_95, pk_96, pk_97, pk_98, pk_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_129[k] = -ab_y[k] * pi_73[k]
                   + pk_95[k];

        t_130[k] = -ab_y[k] * pi_74[k]
                   + pk_96[k];

        t_131[k] = -ab_y[k] * pi_75[k]
                   + pk_97[k];

        t_132[k] = -ab_y[k] * pi_76[k]
                   + pk_98[k];

        t_133[k] = -ab_y[k] * pi_77[k]
                   + pk_100[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, ab_y, pi_78, pi_79, pi_80, pi_81, \
                         pi_82, pk_101, pk_102, pk_103, pk_104, \
                         pk_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_134[k] = -ab_y[k] * pi_78[k]
                   + pk_101[k];

        t_135[k] = -ab_y[k] * pi_79[k]
                   + pk_102[k];

        t_136[k] = -ab_y[k] * pi_80[k]
                   + pk_103[k];

        t_137[k] = -ab_y[k] * pi_81[k]
                   + pk_104[k];

        t_138[k] = -ab_y[k] * pi_82[k]
                   + pk_105[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, ab_y, ab_z, pi_56, pi_57, pi_58, pi_83, \
                         pk_74, pk_76, pk_77, pk_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_139[k] = -ab_y[k] * pi_83[k]
                   + pk_106[k];

        t_140[k] = -ab_z[k] * pi_56[k]
                   + pk_74[k];

        t_141[k] = -ab_z[k] * pi_57[k]
                   + pk_76[k];

        t_142[k] = -ab_z[k] * pi_58[k]
                   + pk_77[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, t_147, ab_z, pi_59, pi_60, pi_61, pi_62, \
                         pi_63, pk_79, pk_80, pk_81, pk_83, pk_84 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_143[k] = -ab_z[k] * pi_59[k]
                   + pk_79[k];

        t_144[k] = -ab_z[k] * pi_60[k]
                   + pk_80[k];

        t_145[k] = -ab_z[k] * pi_61[k]
                   + pk_81[k];

        t_146[k] = -ab_z[k] * pi_62[k]
                   + pk_83[k];

        t_147[k] = -ab_z[k] * pi_63[k]
                   + pk_84[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, ab_z, pi_64, pi_65, pi_66, pi_67, \
                         pi_68, pk_85, pk_86, pk_88, pk_89, pk_90 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_148[k] = -ab_z[k] * pi_64[k]
                   + pk_85[k];

        t_149[k] = -ab_z[k] * pi_65[k]
                   + pk_86[k];

        t_150[k] = -ab_z[k] * pi_66[k]
                   + pk_88[k];

        t_151[k] = -ab_z[k] * pi_67[k]
                   + pk_89[k];

        t_152[k] = -ab_z[k] * pi_68[k]
                   + pk_90[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, ab_z, pi_69, pi_70, pi_71, pi_72, \
                         pi_73, pk_91, pk_92, pk_94, pk_95, pk_96 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_153[k] = -ab_z[k] * pi_69[k]
                   + pk_91[k];

        t_154[k] = -ab_z[k] * pi_70[k]
                   + pk_92[k];

        t_155[k] = -ab_z[k] * pi_71[k]
                   + pk_94[k];

        t_156[k] = -ab_z[k] * pi_72[k]
                   + pk_95[k];

        t_157[k] = -ab_z[k] * pi_73[k]
                   + pk_96[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, ab_z, pi_74, pi_75, pi_76, pi_77, \
                         pi_78, pk_97, pk_98, pk_99, pk_101, pk_102 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_158[k] = -ab_z[k] * pi_74[k]
                   + pk_97[k];

        t_159[k] = -ab_z[k] * pi_75[k]
                   + pk_98[k];

        t_160[k] = -ab_z[k] * pi_76[k]
                   + pk_99[k];

        t_161[k] = -ab_z[k] * pi_77[k]
                   + pk_101[k];

        t_162[k] = -ab_z[k] * pi_78[k]
                   + pk_102[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, ab_z, pi_79, pi_80, pi_81, pi_82, \
                         pi_83, pk_103, pk_104, pk_105, pk_106, \
                         pk_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_163[k] = -ab_z[k] * pi_79[k]
                   + pk_103[k];

        t_164[k] = -ab_z[k] * pi_80[k]
                   + pk_104[k];

        t_165[k] = -ab_z[k] * pi_81[k]
                   + pk_105[k];

        t_166[k] = -ab_z[k] * pi_82[k]
                   + pk_106[k];

        t_167[k] = -ab_z[k] * pi_83[k]
                   + pk_107[k];
    }
}

}  // namespace simdtrf
