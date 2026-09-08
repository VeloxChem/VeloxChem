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


#include "SimdTransformID.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_id(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t id,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5625 * std::sqrt(154.0);
    const auto f_1 = 1.875 * std::sqrt(154.0);
    const auto f_2 = 0.09375 * std::sqrt(462.0);
    const auto f_3 = 0.1875 * std::sqrt(462.0);
    const auto f_4 = 0.3125 * std::sqrt(462.0);
    const auto f_5 = 0.625 * std::sqrt(462.0);
    const auto f_6 = 0.28125 * std::sqrt(154.0);
    const auto f_7 = 0.9375 * std::sqrt(154.0);
    const auto f_8 = 0.9375 * std::sqrt(462.0);
    const auto f_9 = 1.875 * std::sqrt(462.0);
    const auto f_10 = 0.46875 * std::sqrt(154.0);
    const auto f_11 = 0.09375 * std::sqrt(154.0);
    const auto f_12 = 0.1875 * std::sqrt(154.0);
    const auto f_13 = 0.46875 * std::sqrt(462.0);
    const auto f_14 = 0.75 * std::sqrt(21.0);
    const auto f_15 = 7.5 * std::sqrt(21.0);
    const auto f_16 = 0.375 * std::sqrt(7.0);
    const auto f_17 = 0.75 * std::sqrt(7.0);
    const auto f_18 = 3.75 * std::sqrt(7.0);
    const auto f_19 = 7.5 * std::sqrt(7.0);
    const auto f_20 = 0.375 * std::sqrt(21.0);
    const auto f_21 = 3.75 * std::sqrt(21.0);
    const auto f_22 = 1.6875 * std::sqrt(70.0);
    const auto f_23 = 1.125 * std::sqrt(70.0);
    const auto f_24 = 4.5 * std::sqrt(70.0);
    const auto f_25 = 0.5625 * std::sqrt(70.0);
    const auto f_26 = 1.5 * std::sqrt(70.0);
    const auto f_27 = 0.28125 * std::sqrt(210.0);
    const auto f_28 = 0.5625 * std::sqrt(210.0);
    const auto f_29 = 0.1875 * std::sqrt(210.0);
    const auto f_30 = 0.375 * std::sqrt(210.0);
    const auto f_31 = 0.75 * std::sqrt(210.0);
    const auto f_32 = 1.5 * std::sqrt(210.0);
    const auto f_33 = 0.09375 * std::sqrt(210.0);
    const auto f_34 = 0.25 * std::sqrt(210.0);
    const auto f_35 = 0.5 * std::sqrt(210.0);
    const auto f_36 = 0.84375 * std::sqrt(70.0);
    const auto f_37 = 2.25 * std::sqrt(70.0);
    const auto f_38 = 0.28125 * std::sqrt(70.0);
    const auto f_39 = 0.75 * std::sqrt(70.0);
    const auto f_40 = 0.1875 * std::sqrt(70.0);
    const auto f_41 = 0.375 * std::sqrt(70.0);
    const auto f_42 = 3.0 * std::sqrt(70.0);
    const auto f_43 = 0.03125 * std::sqrt(210.0);
    const auto f_44 = 0.0625 * std::sqrt(210.0);
    const auto f_45 = 0.125 * std::sqrt(210.0);
    const auto f_46 = std::sqrt(210.0);
    const auto f_47 = 0.09375 * std::sqrt(70.0);
    const auto f_48 = 1.875 * std::sqrt(7.0);
    const auto f_49 = 3.0 * std::sqrt(7.0);
    const auto f_50 = 0.3125 * std::sqrt(21.0);
    const auto f_51 = 0.625 * std::sqrt(21.0);
    const auto f_52 = 1.25 * std::sqrt(21.0);
    const auto f_53 = 2.5 * std::sqrt(21.0);
    const auto f_54 = 0.5 * std::sqrt(21.0);
    const auto f_55 = std::sqrt(21.0);
    const auto f_56 = 0.9375 * std::sqrt(7.0);
    const auto f_57 = 1.5 * std::sqrt(7.0);
    const auto f_58 = 0.3125 * std::sqrt(3.0);
    const auto f_59 = 0.9375 * std::sqrt(3.0);
    const auto f_60 = 5.625 * std::sqrt(3.0);
    const auto f_61 = 11.25 * std::sqrt(3.0);
    const auto f_62 = 7.5 * std::sqrt(3.0);
    const auto f_63 = std::sqrt(3.0);
    const auto f_64 = 0.15625 * std::sqrt(3.0);
    const auto f_65 = 0.46875 * std::sqrt(3.0);
    const auto f_66 = 2.8125 * std::sqrt(3.0);
    const auto f_67 = 3.75 * std::sqrt(3.0);
    const auto f_68 = 0.5 * std::sqrt(3.0);
    const auto f_69 = 0.015625 * std::sqrt(210.0);
    const auto f_70 = 0.046875 * std::sqrt(70.0);
    const auto f_71 = 0.1875 * std::sqrt(21.0);
    const auto f_72 = 0.9375 * std::sqrt(21.0);
    const auto f_73 = 1.875 * std::sqrt(21.0);
    const auto f_74 = 11.25 * std::sqrt(21.0);
    const auto f_75 = 0.09375 * std::sqrt(7.0);
    const auto f_76 = 0.1875 * std::sqrt(7.0);
    const auto f_77 = 0.46875 * std::sqrt(7.0);
    const auto f_78 = 5.625 * std::sqrt(7.0);
    const auto f_79 = 11.25 * std::sqrt(7.0);
    const auto f_80 = 0.09375 * std::sqrt(21.0);
    const auto f_81 = 0.46875 * std::sqrt(21.0);
    const auto f_82 = 5.625 * std::sqrt(21.0);
    const auto f_83 = 1.40625 * std::sqrt(154.0);
    const auto f_84 = 0.015625 * std::sqrt(462.0);
    const auto f_85 = 0.03125 * std::sqrt(462.0);
    const auto f_86 = 0.234375 * std::sqrt(462.0);
    const auto f_87 = 0.046875 * std::sqrt(154.0);
    const auto f_88 = 0.703125 * std::sqrt(154.0);

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
    const auto *id_72 = buffer.data(id + 72);
    const auto *id_73 = buffer.data(id + 73);
    const auto *id_74 = buffer.data(id + 74);
    const auto *id_75 = buffer.data(id + 75);
    const auto *id_76 = buffer.data(id + 76);
    const auto *id_77 = buffer.data(id + 77);
    const auto *id_78 = buffer.data(id + 78);
    const auto *id_79 = buffer.data(id + 79);
    const auto *id_80 = buffer.data(id + 80);
    const auto *id_81 = buffer.data(id + 81);
    const auto *id_82 = buffer.data(id + 82);
    const auto *id_83 = buffer.data(id + 83);
    const auto *id_84 = buffer.data(id + 84);
    const auto *id_85 = buffer.data(id + 85);
    const auto *id_86 = buffer.data(id + 86);
    const auto *id_87 = buffer.data(id + 87);
    const auto *id_88 = buffer.data(id + 88);
    const auto *id_89 = buffer.data(id + 89);
    const auto *id_90 = buffer.data(id + 90);
    const auto *id_91 = buffer.data(id + 91);
    const auto *id_92 = buffer.data(id + 92);
    const auto *id_93 = buffer.data(id + 93);
    const auto *id_94 = buffer.data(id + 94);
    const auto *id_95 = buffer.data(id + 95);
    const auto *id_96 = buffer.data(id + 96);
    const auto *id_97 = buffer.data(id + 97);
    const auto *id_98 = buffer.data(id + 98);
    const auto *id_99 = buffer.data(id + 99);
    const auto *id_100 = buffer.data(id + 100);
    const auto *id_101 = buffer.data(id + 101);
    const auto *id_102 = buffer.data(id + 102);
    const auto *id_103 = buffer.data(id + 103);
    const auto *id_104 = buffer.data(id + 104);
    const auto *id_105 = buffer.data(id + 105);
    const auto *id_106 = buffer.data(id + 106);
    const auto *id_107 = buffer.data(id + 107);
    const auto *id_108 = buffer.data(id + 108);
    const auto *id_109 = buffer.data(id + 109);
    const auto *id_110 = buffer.data(id + 110);
    const auto *id_111 = buffer.data(id + 111);
    const auto *id_112 = buffer.data(id + 112);
    const auto *id_113 = buffer.data(id + 113);
    const auto *id_114 = buffer.data(id + 114);
    const auto *id_115 = buffer.data(id + 115);
    const auto *id_116 = buffer.data(id + 116);
    const auto *id_117 = buffer.data(id + 117);
    const auto *id_118 = buffer.data(id + 118);
    const auto *id_119 = buffer.data(id + 119);
    const auto *id_120 = buffer.data(id + 120);
    const auto *id_121 = buffer.data(id + 121);
    const auto *id_122 = buffer.data(id + 122);
    const auto *id_123 = buffer.data(id + 123);
    const auto *id_124 = buffer.data(id + 124);
    const auto *id_125 = buffer.data(id + 125);
    const auto *id_126 = buffer.data(id + 126);
    const auto *id_127 = buffer.data(id + 127);
    const auto *id_128 = buffer.data(id + 128);
    const auto *id_129 = buffer.data(id + 129);
    const auto *id_130 = buffer.data(id + 130);
    const auto *id_131 = buffer.data(id + 131);
    const auto *id_132 = buffer.data(id + 132);
    const auto *id_133 = buffer.data(id + 133);
    const auto *id_134 = buffer.data(id + 134);
    const auto *id_135 = buffer.data(id + 135);
    const auto *id_136 = buffer.data(id + 136);
    const auto *id_137 = buffer.data(id + 137);
    const auto *id_138 = buffer.data(id + 138);
    const auto *id_139 = buffer.data(id + 139);
    const auto *id_140 = buffer.data(id + 140);
    const auto *id_141 = buffer.data(id + 141);
    const auto *id_142 = buffer.data(id + 142);
    const auto *id_143 = buffer.data(id + 143);
    const auto *id_144 = buffer.data(id + 144);
    const auto *id_145 = buffer.data(id + 145);
    const auto *id_146 = buffer.data(id + 146);
    const auto *id_147 = buffer.data(id + 147);
    const auto *id_148 = buffer.data(id + 148);
    const auto *id_149 = buffer.data(id + 149);
    const auto *id_150 = buffer.data(id + 150);
    const auto *id_151 = buffer.data(id + 151);
    const auto *id_152 = buffer.data(id + 152);
    const auto *id_153 = buffer.data(id + 153);
    const auto *id_154 = buffer.data(id + 154);
    const auto *id_155 = buffer.data(id + 155);
    const auto *id_156 = buffer.data(id + 156);
    const auto *id_157 = buffer.data(id + 157);
    const auto *id_158 = buffer.data(id + 158);
    const auto *id_159 = buffer.data(id + 159);
    const auto *id_160 = buffer.data(id + 160);
    const auto *id_161 = buffer.data(id + 161);
    const auto *id_162 = buffer.data(id + 162);
    const auto *id_163 = buffer.data(id + 163);
    const auto *id_164 = buffer.data(id + 164);
    const auto *id_165 = buffer.data(id + 165);
    const auto *id_166 = buffer.data(id + 166);
    const auto *id_167 = buffer.data(id + 167);

#pragma omp simd aligned(id_7, id_10, id_37, id_40, id_91, id_94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * id_7[k]
                 - f_1 * id_37[k]
                 + f_0 * id_91[k];

        g_1[k] = f_0 * id_10[k]
                 - f_1 * id_40[k]
                 + f_0 * id_94[k];
    }

#pragma omp simd aligned(id_6, id_8, id_9, id_11, id_36, id_38, id_39, id_41, id_90, id_92, \
                         id_93, id_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_2 * id_6[k]
                 - f_2 * id_9[k]
                 + f_3 * id_11[k]
                 + f_4 * id_36[k]
                 + f_4 * id_39[k]
                 - f_5 * id_41[k]
                 - f_2 * id_90[k]
                 - f_2 * id_93[k]
                 + f_3 * id_95[k];

        g_3[k] = f_0 * id_8[k]
                 - f_1 * id_38[k]
                 + f_0 * id_92[k];

        g_4[k] = f_6 * id_6[k]
                 - f_6 * id_9[k]
                 - f_7 * id_36[k]
                 + f_7 * id_39[k]
                 + f_6 * id_90[k]
                 - f_6 * id_93[k];
    }

#pragma omp simd aligned(id_25, id_28, id_67, id_70, id_133, id_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_8 * id_25[k]
                 - f_9 * id_67[k]
                 + f_3 * id_133[k];

        g_6[k] = f_8 * id_28[k]
                 - f_9 * id_70[k]
                 + f_3 * id_136[k];
    }

#pragma omp simd aligned(id_24, id_26, id_27, id_29, id_66, id_68, id_69, id_71, id_132, \
                         id_134, id_135, id_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_10 * id_24[k]
                 - f_10 * id_27[k]
                 + f_7 * id_29[k]
                 + f_7 * id_66[k]
                 + f_7 * id_69[k]
                 - f_1 * id_71[k]
                 - f_11 * id_132[k]
                 - f_11 * id_135[k]
                 + f_12 * id_137[k];

        g_8[k] = f_8 * id_26[k]
                 - f_9 * id_68[k]
                 + f_3 * id_134[k];

        g_9[k] = f_13 * id_24[k]
                 - f_13 * id_27[k]
                 - f_8 * id_66[k]
                 + f_8 * id_69[k]
                 + f_2 * id_132[k]
                 - f_2 * id_135[k];
    }

#pragma omp simd aligned(id_7, id_10, id_49, id_52, id_91, id_94, id_103, \
                         id_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_14 * id_7[k]
                  + f_15 * id_49[k]
                  + f_14 * id_91[k]
                  - f_15 * id_103[k];

        g_11[k] = -f_14 * id_10[k]
                  + f_15 * id_52[k]
                  + f_14 * id_94[k]
                  - f_15 * id_106[k];
    }

#pragma omp simd aligned(id_6, id_9, id_11, id_48, id_51, id_53, id_90, id_93, id_95, id_102, \
                         id_105, id_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_16 * id_6[k]
                  + f_16 * id_9[k]
                  - f_17 * id_11[k]
                  - f_18 * id_48[k]
                  - f_18 * id_51[k]
                  + f_19 * id_53[k]
                  - f_16 * id_90[k]
                  - f_16 * id_93[k]
                  + f_17 * id_95[k]
                  + f_18 * id_102[k]
                  + f_18 * id_105[k]
                  - f_19 * id_107[k];
    }

#pragma omp simd aligned(id_6, id_8, id_9, id_48, id_50, id_51, id_90, id_92, id_93, id_102, \
                         id_104, id_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_14 * id_8[k]
                  + f_15 * id_50[k]
                  + f_14 * id_92[k]
                  - f_15 * id_104[k];

        g_14[k] = -f_20 * id_6[k]
                  + f_20 * id_9[k]
                  + f_21 * id_48[k]
                  - f_21 * id_51[k]
                  + f_20 * id_90[k]
                  - f_20 * id_93[k]
                  - f_21 * id_102[k]
                  + f_21 * id_105[k];
    }

#pragma omp simd aligned(id_25, id_28, id_67, id_70, id_79, id_82, id_133, id_136, id_145, \
                         id_148 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_22 * id_25[k]
                  - f_23 * id_67[k]
                  + f_24 * id_79[k]
                  + f_25 * id_133[k]
                  - f_26 * id_145[k];

        g_16[k] = -f_22 * id_28[k]
                  - f_23 * id_70[k]
                  + f_24 * id_82[k]
                  + f_25 * id_136[k]
                  - f_26 * id_148[k];
    }

#pragma omp simd aligned(id_24, id_27, id_29, id_66, id_69, id_71, id_78, id_81, id_83, \
                         id_132, id_135, id_137, id_144, id_147, \
                         id_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_27 * id_24[k]
                  + f_27 * id_27[k]
                  - f_28 * id_29[k]
                  + f_29 * id_66[k]
                  + f_29 * id_69[k]
                  - f_30 * id_71[k]
                  - f_31 * id_78[k]
                  - f_31 * id_81[k]
                  + f_32 * id_83[k]
                  - f_33 * id_132[k]
                  - f_33 * id_135[k]
                  + f_29 * id_137[k]
                  + f_34 * id_144[k]
                  + f_34 * id_147[k]
                  - f_35 * id_149[k];
    }

#pragma omp simd aligned(id_26, id_68, id_80, id_134, id_146 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_22 * id_26[k]
                  - f_23 * id_68[k]
                  + f_24 * id_80[k]
                  + f_25 * id_134[k]
                  - f_26 * id_146[k];
    }

#pragma omp simd aligned(id_24, id_27, id_66, id_69, id_78, id_81, id_132, id_135, id_144, \
                         id_147 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_36 * id_24[k]
                  + f_36 * id_27[k]
                  - f_25 * id_66[k]
                  + f_25 * id_69[k]
                  + f_37 * id_78[k]
                  - f_37 * id_81[k]
                  + f_38 * id_132[k]
                  - f_38 * id_135[k]
                  - f_39 * id_144[k]
                  + f_39 * id_147[k];
    }

#pragma omp simd aligned(id_7, id_10, id_37, id_40, id_49, id_52, id_91, id_94, id_103, \
                         id_106, id_115, id_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_40 * id_7[k]
                  + f_41 * id_37[k]
                  - f_42 * id_49[k]
                  + f_40 * id_91[k]
                  - f_42 * id_103[k]
                  + f_42 * id_115[k];

        g_21[k] = f_40 * id_10[k]
                  + f_41 * id_40[k]
                  - f_42 * id_52[k]
                  + f_40 * id_94[k]
                  - f_42 * id_106[k]
                  + f_42 * id_118[k];
    }

#pragma omp simd aligned(id_6, id_9, id_11, id_36, id_39, id_41, id_48, id_51, id_53, id_90, \
                         id_93, id_95, id_102, id_105, id_107, id_114, id_117, \
                         id_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_43 * id_6[k]
                  - f_43 * id_9[k]
                  + f_44 * id_11[k]
                  - f_44 * id_36[k]
                  - f_44 * id_39[k]
                  + f_45 * id_41[k]
                  + f_35 * id_48[k]
                  + f_35 * id_51[k]
                  - f_46 * id_53[k]
                  - f_43 * id_90[k]
                  - f_43 * id_93[k]
                  + f_44 * id_95[k]
                  + f_35 * id_102[k]
                  + f_35 * id_105[k]
                  - f_46 * id_107[k]
                  - f_35 * id_114[k]
                  - f_35 * id_117[k]
                  + f_46 * id_119[k];
    }

#pragma omp simd aligned(id_8, id_38, id_50, id_92, id_104, id_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_40 * id_8[k]
                  + f_41 * id_38[k]
                  - f_42 * id_50[k]
                  + f_40 * id_92[k]
                  - f_42 * id_104[k]
                  + f_42 * id_116[k];
    }

#pragma omp simd aligned(id_6, id_9, id_36, id_39, id_48, id_51, id_90, id_93, id_102, id_105, \
                         id_114, id_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_47 * id_6[k]
                  - f_47 * id_9[k]
                  + f_40 * id_36[k]
                  - f_40 * id_39[k]
                  - f_26 * id_48[k]
                  + f_26 * id_51[k]
                  + f_47 * id_90[k]
                  - f_47 * id_93[k]
                  - f_26 * id_102[k]
                  + f_26 * id_105[k]
                  + f_26 * id_114[k]
                  - f_26 * id_117[k];
    }

#pragma omp simd aligned(id_25, id_28, id_67, id_70, id_79, id_82, id_133, id_136, id_145, \
                         id_148, id_157, id_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_48 * id_25[k]
                  + f_18 * id_67[k]
                  - f_19 * id_79[k]
                  + f_48 * id_133[k]
                  - f_19 * id_145[k]
                  + f_49 * id_157[k];

        g_26[k] = f_48 * id_28[k]
                  + f_18 * id_70[k]
                  - f_19 * id_82[k]
                  + f_48 * id_136[k]
                  - f_19 * id_148[k]
                  + f_49 * id_160[k];
    }

#pragma omp simd aligned(id_24, id_27, id_29, id_66, id_69, id_71, id_78, id_81, id_83, \
                         id_132, id_135, id_137, id_144, id_147, id_149, id_156, id_159, \
                         id_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_50 * id_24[k]
                  - f_50 * id_27[k]
                  + f_51 * id_29[k]
                  - f_51 * id_66[k]
                  - f_51 * id_69[k]
                  + f_52 * id_71[k]
                  + f_52 * id_78[k]
                  + f_52 * id_81[k]
                  - f_53 * id_83[k]
                  - f_50 * id_132[k]
                  - f_50 * id_135[k]
                  + f_51 * id_137[k]
                  + f_52 * id_144[k]
                  + f_52 * id_147[k]
                  - f_53 * id_149[k]
                  - f_54 * id_156[k]
                  - f_54 * id_159[k]
                  + f_55 * id_161[k];
    }

#pragma omp simd aligned(id_26, id_68, id_80, id_134, id_146, id_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_48 * id_26[k]
                  + f_18 * id_68[k]
                  - f_19 * id_80[k]
                  + f_48 * id_134[k]
                  - f_19 * id_146[k]
                  + f_49 * id_158[k];
    }

#pragma omp simd aligned(id_24, id_27, id_66, id_69, id_78, id_81, id_132, id_135, id_144, \
                         id_147, id_156, id_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_56 * id_24[k]
                  - f_56 * id_27[k]
                  + f_48 * id_66[k]
                  - f_48 * id_69[k]
                  - f_18 * id_78[k]
                  + f_18 * id_81[k]
                  + f_56 * id_132[k]
                  - f_56 * id_135[k]
                  - f_18 * id_144[k]
                  + f_18 * id_147[k]
                  + f_57 * id_156[k]
                  - f_57 * id_159[k];
    }

#pragma omp simd aligned(id_1, id_19, id_31, id_61, id_73, id_85, id_127, id_139, id_151, \
                         id_163 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_58 * id_1[k]
                  - f_59 * id_19[k]
                  + f_60 * id_31[k]
                  - f_59 * id_61[k]
                  + f_61 * id_73[k]
                  - f_62 * id_85[k]
                  - f_58 * id_127[k]
                  + f_60 * id_139[k]
                  - f_62 * id_151[k]
                  + f_63 * id_163[k];
    }

#pragma omp simd aligned(id_4, id_22, id_34, id_64, id_76, id_88, id_130, id_142, id_154, \
                         id_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_58 * id_4[k]
                  - f_59 * id_22[k]
                  + f_60 * id_34[k]
                  - f_59 * id_64[k]
                  + f_61 * id_76[k]
                  - f_62 * id_88[k]
                  - f_58 * id_130[k]
                  + f_60 * id_142[k]
                  - f_62 * id_154[k]
                  + f_63 * id_166[k];
    }

#pragma omp simd aligned(id_0, id_3, id_5, id_18, id_21, id_23, id_30, id_33, id_35, id_60, \
                         id_63, id_65, id_72, id_75, id_77, id_84, id_87, id_89, id_126, \
                         id_129, id_131, id_138, id_141, id_143, id_150, id_153, id_155, \
                         id_162, id_165, id_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = 0.15625 * id_0[k]
                  + 0.15625 * id_3[k]
                  - 0.3125 * id_5[k]
                  + 0.46875 * id_18[k]
                  + 0.46875 * id_21[k]
                  - 0.9375 * id_23[k]
                  - 2.8125 * id_30[k]
                  - 2.8125 * id_33[k]
                  + 5.625 * id_35[k]
                  + 0.46875 * id_60[k]
                  + 0.46875 * id_63[k]
                  - 0.9375 * id_65[k]
                  - 5.625 * id_72[k]
                  - 5.625 * id_75[k]
                  + 11.25 * id_77[k]
                  + 3.75 * id_84[k]
                  + 3.75 * id_87[k]
                  - 7.5 * id_89[k]
                  + 0.15625 * id_126[k]
                  + 0.15625 * id_129[k]
                  - 0.3125 * id_131[k]
                  - 2.8125 * id_138[k]
                  - 2.8125 * id_141[k]
                  + 5.625 * id_143[k]
                  + 3.75 * id_150[k]
                  + 3.75 * id_153[k]
                  - 7.5 * id_155[k]
                  - 0.5 * id_162[k]
                  - 0.5 * id_165[k]
                  + id_167[k];
    }

#pragma omp simd aligned(id_2, id_20, id_32, id_62, id_74, id_86, id_128, id_140, id_152, \
                         id_164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_58 * id_2[k]
                  - f_59 * id_20[k]
                  + f_60 * id_32[k]
                  - f_59 * id_62[k]
                  + f_61 * id_74[k]
                  - f_62 * id_86[k]
                  - f_58 * id_128[k]
                  + f_60 * id_140[k]
                  - f_62 * id_152[k]
                  + f_63 * id_164[k];
    }

#pragma omp simd aligned(id_0, id_3, id_18, id_21, id_30, id_33, id_60, id_63, id_72, id_75, \
                         id_84, id_87, id_126, id_129, id_138, id_141, id_150, id_153, id_162, \
                         id_165 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_64 * id_0[k]
                  + f_64 * id_3[k]
                  - f_65 * id_18[k]
                  + f_65 * id_21[k]
                  + f_66 * id_30[k]
                  - f_66 * id_33[k]
                  - f_65 * id_60[k]
                  + f_65 * id_63[k]
                  + f_60 * id_72[k]
                  - f_60 * id_75[k]
                  - f_67 * id_84[k]
                  + f_67 * id_87[k]
                  - f_64 * id_126[k]
                  + f_64 * id_129[k]
                  + f_66 * id_138[k]
                  - f_66 * id_141[k]
                  - f_67 * id_150[k]
                  + f_67 * id_153[k]
                  + f_68 * id_162[k]
                  - f_68 * id_165[k];
    }

#pragma omp simd aligned(id_13, id_16, id_43, id_46, id_55, id_58, id_97, id_100, id_109, \
                         id_112, id_121, id_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_48 * id_13[k]
                  + f_18 * id_43[k]
                  - f_19 * id_55[k]
                  + f_48 * id_97[k]
                  - f_19 * id_109[k]
                  + f_49 * id_121[k];

        g_36[k] = f_48 * id_16[k]
                  + f_18 * id_46[k]
                  - f_19 * id_58[k]
                  + f_48 * id_100[k]
                  - f_19 * id_112[k]
                  + f_49 * id_124[k];
    }

#pragma omp simd aligned(id_12, id_15, id_17, id_42, id_45, id_47, id_54, id_57, id_59, id_96, \
                         id_99, id_101, id_108, id_111, id_113, id_120, id_123, \
                         id_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_50 * id_12[k]
                  - f_50 * id_15[k]
                  + f_51 * id_17[k]
                  - f_51 * id_42[k]
                  - f_51 * id_45[k]
                  + f_52 * id_47[k]
                  + f_52 * id_54[k]
                  + f_52 * id_57[k]
                  - f_53 * id_59[k]
                  - f_50 * id_96[k]
                  - f_50 * id_99[k]
                  + f_51 * id_101[k]
                  + f_52 * id_108[k]
                  + f_52 * id_111[k]
                  - f_53 * id_113[k]
                  - f_54 * id_120[k]
                  - f_54 * id_123[k]
                  + f_55 * id_125[k];
    }

#pragma omp simd aligned(id_14, id_44, id_56, id_98, id_110, id_122 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_48 * id_14[k]
                  + f_18 * id_44[k]
                  - f_19 * id_56[k]
                  + f_48 * id_98[k]
                  - f_19 * id_110[k]
                  + f_49 * id_122[k];
    }

#pragma omp simd aligned(id_12, id_15, id_42, id_45, id_54, id_57, id_96, id_99, id_108, \
                         id_111, id_120, id_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = f_56 * id_12[k]
                  - f_56 * id_15[k]
                  + f_48 * id_42[k]
                  - f_48 * id_45[k]
                  - f_18 * id_54[k]
                  + f_18 * id_57[k]
                  + f_56 * id_96[k]
                  - f_56 * id_99[k]
                  - f_18 * id_108[k]
                  + f_18 * id_111[k]
                  + f_57 * id_120[k]
                  - f_57 * id_123[k];
    }

#pragma omp simd aligned(id_1, id_19, id_31, id_61, id_85, id_127, id_139, \
                         id_151 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_47 * id_1[k]
                  + f_47 * id_19[k]
                  - f_26 * id_31[k]
                  - f_47 * id_61[k]
                  + f_26 * id_85[k]
                  - f_47 * id_127[k]
                  + f_26 * id_139[k]
                  - f_26 * id_151[k];
    }

#pragma omp simd aligned(id_4, id_22, id_34, id_64, id_88, id_130, id_142, \
                         id_154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_47 * id_4[k]
                  + f_47 * id_22[k]
                  - f_26 * id_34[k]
                  - f_47 * id_64[k]
                  + f_26 * id_88[k]
                  - f_47 * id_130[k]
                  + f_26 * id_142[k]
                  - f_26 * id_154[k];
    }

#pragma omp simd aligned(id_0, id_3, id_5, id_18, id_21, id_23, id_30, id_33, id_35, id_60, \
                         id_63, id_65, id_84, id_87, id_89, id_126, id_129, id_131, id_138, \
                         id_141, id_143, id_150, id_153, id_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_69 * id_0[k]
                  - f_69 * id_3[k]
                  + f_43 * id_5[k]
                  - f_69 * id_18[k]
                  - f_69 * id_21[k]
                  + f_43 * id_23[k]
                  + f_34 * id_30[k]
                  + f_34 * id_33[k]
                  - f_35 * id_35[k]
                  + f_69 * id_60[k]
                  + f_69 * id_63[k]
                  - f_43 * id_65[k]
                  - f_34 * id_84[k]
                  - f_34 * id_87[k]
                  + f_35 * id_89[k]
                  + f_69 * id_126[k]
                  + f_69 * id_129[k]
                  - f_43 * id_131[k]
                  - f_34 * id_138[k]
                  - f_34 * id_141[k]
                  + f_35 * id_143[k]
                  + f_34 * id_150[k]
                  + f_34 * id_153[k]
                  - f_35 * id_155[k];
    }

#pragma omp simd aligned(id_2, id_20, id_32, id_62, id_86, id_128, id_140, \
                         id_152 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_47 * id_2[k]
                  + f_47 * id_20[k]
                  - f_26 * id_32[k]
                  - f_47 * id_62[k]
                  + f_26 * id_86[k]
                  - f_47 * id_128[k]
                  + f_26 * id_140[k]
                  - f_26 * id_152[k];
    }

#pragma omp simd aligned(id_0, id_3, id_18, id_21, id_30, id_33, id_60, id_63, id_84, id_87, \
                         id_126, id_129, id_138, id_141, id_150, \
                         id_153 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_70 * id_0[k]
                  - f_70 * id_3[k]
                  + f_70 * id_18[k]
                  - f_70 * id_21[k]
                  - f_39 * id_30[k]
                  + f_39 * id_33[k]
                  - f_70 * id_60[k]
                  + f_70 * id_63[k]
                  + f_39 * id_84[k]
                  - f_39 * id_87[k]
                  - f_70 * id_126[k]
                  + f_70 * id_129[k]
                  + f_39 * id_138[k]
                  - f_39 * id_141[k]
                  - f_39 * id_150[k]
                  + f_39 * id_153[k];
    }

#pragma omp simd aligned(id_13, id_16, id_43, id_46, id_55, id_58, id_97, id_100, id_109, \
                         id_112 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_25 * id_13[k]
                  + f_23 * id_43[k]
                  + f_26 * id_55[k]
                  + f_22 * id_97[k]
                  - f_24 * id_109[k];

        g_46[k] = -f_25 * id_16[k]
                  + f_23 * id_46[k]
                  + f_26 * id_58[k]
                  + f_22 * id_100[k]
                  - f_24 * id_112[k];
    }

#pragma omp simd aligned(id_12, id_15, id_17, id_42, id_45, id_47, id_54, id_57, id_59, id_96, \
                         id_99, id_101, id_108, id_111, id_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = f_33 * id_12[k]
                  + f_33 * id_15[k]
                  - f_29 * id_17[k]
                  - f_29 * id_42[k]
                  - f_29 * id_45[k]
                  + f_30 * id_47[k]
                  - f_34 * id_54[k]
                  - f_34 * id_57[k]
                  + f_35 * id_59[k]
                  - f_27 * id_96[k]
                  - f_27 * id_99[k]
                  + f_28 * id_101[k]
                  + f_31 * id_108[k]
                  + f_31 * id_111[k]
                  - f_32 * id_113[k];
    }

#pragma omp simd aligned(id_14, id_44, id_56, id_98, id_110 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_25 * id_14[k]
                  + f_23 * id_44[k]
                  + f_26 * id_56[k]
                  + f_22 * id_98[k]
                  - f_24 * id_110[k];
    }

#pragma omp simd aligned(id_12, id_15, id_42, id_45, id_54, id_57, id_96, id_99, id_108, \
                         id_111 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_38 * id_12[k]
                  + f_38 * id_15[k]
                  + f_25 * id_42[k]
                  - f_25 * id_45[k]
                  + f_39 * id_54[k]
                  - f_39 * id_57[k]
                  + f_36 * id_96[k]
                  - f_36 * id_99[k]
                  - f_37 * id_108[k]
                  + f_37 * id_111[k];
    }

#pragma omp simd aligned(id_1, id_4, id_19, id_22, id_31, id_34, id_61, id_64, id_73, id_76, \
                         id_127, id_130, id_139, id_142 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_71 * id_1[k]
                  + f_72 * id_19[k]
                  + f_73 * id_31[k]
                  + f_72 * id_61[k]
                  - f_74 * id_73[k]
                  - f_71 * id_127[k]
                  + f_73 * id_139[k];

        g_51[k] = -f_71 * id_4[k]
                  + f_72 * id_22[k]
                  + f_73 * id_34[k]
                  + f_72 * id_64[k]
                  - f_74 * id_76[k]
                  - f_71 * id_130[k]
                  + f_73 * id_142[k];
    }

#pragma omp simd aligned(id_0, id_3, id_5, id_18, id_21, id_23, id_30, id_33, id_35, id_60, \
                         id_63, id_65, id_72, id_75, id_77, id_126, id_129, id_131, id_138, \
                         id_141, id_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_75 * id_0[k]
                  + f_75 * id_3[k]
                  - f_76 * id_5[k]
                  - f_77 * id_18[k]
                  - f_77 * id_21[k]
                  + f_56 * id_23[k]
                  - f_56 * id_30[k]
                  - f_56 * id_33[k]
                  + f_48 * id_35[k]
                  - f_77 * id_60[k]
                  - f_77 * id_63[k]
                  + f_56 * id_65[k]
                  + f_78 * id_72[k]
                  + f_78 * id_75[k]
                  - f_79 * id_77[k]
                  + f_75 * id_126[k]
                  + f_75 * id_129[k]
                  - f_76 * id_131[k]
                  - f_56 * id_138[k]
                  - f_56 * id_141[k]
                  + f_48 * id_143[k];
    }

#pragma omp simd aligned(id_2, id_20, id_32, id_62, id_74, id_128, \
                         id_140 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = -f_71 * id_2[k]
                  + f_72 * id_20[k]
                  + f_73 * id_32[k]
                  + f_72 * id_62[k]
                  - f_74 * id_74[k]
                  - f_71 * id_128[k]
                  + f_73 * id_140[k];
    }

#pragma omp simd aligned(id_0, id_3, id_18, id_21, id_30, id_33, id_60, id_63, id_72, id_75, \
                         id_126, id_129, id_138, id_141 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_80 * id_0[k]
                  + f_80 * id_3[k]
                  + f_81 * id_18[k]
                  - f_81 * id_21[k]
                  + f_72 * id_30[k]
                  - f_72 * id_33[k]
                  + f_81 * id_60[k]
                  - f_81 * id_63[k]
                  - f_82 * id_72[k]
                  + f_82 * id_75[k]
                  - f_80 * id_126[k]
                  + f_80 * id_129[k]
                  + f_72 * id_138[k]
                  - f_72 * id_141[k];
    }

#pragma omp simd aligned(id_13, id_16, id_43, id_46, id_97, id_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_3 * id_13[k]
                  - f_9 * id_43[k]
                  + f_8 * id_97[k];

        g_56[k] = f_3 * id_16[k]
                  - f_9 * id_46[k]
                  + f_8 * id_100[k];
    }

#pragma omp simd aligned(id_12, id_14, id_15, id_17, id_42, id_44, id_45, id_47, id_96, id_98, \
                         id_99, id_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = -f_11 * id_12[k]
                  - f_11 * id_15[k]
                  + f_12 * id_17[k]
                  + f_7 * id_42[k]
                  + f_7 * id_45[k]
                  - f_1 * id_47[k]
                  - f_10 * id_96[k]
                  - f_10 * id_99[k]
                  + f_7 * id_101[k];

        g_58[k] = f_3 * id_14[k]
                  - f_9 * id_44[k]
                  + f_8 * id_98[k];

        g_59[k] = f_2 * id_12[k]
                  - f_2 * id_15[k]
                  - f_8 * id_42[k]
                  + f_8 * id_45[k]
                  + f_13 * id_96[k]
                  - f_13 * id_99[k];
    }

#pragma omp simd aligned(id_1, id_4, id_19, id_22, id_61, id_64, id_127, \
                         id_130 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_11 * id_1[k]
                  - f_83 * id_19[k]
                  + f_83 * id_61[k]
                  - f_11 * id_127[k];

        g_61[k] = f_11 * id_4[k]
                  - f_83 * id_22[k]
                  + f_83 * id_64[k]
                  - f_11 * id_130[k];
    }

#pragma omp simd aligned(id_0, id_3, id_5, id_18, id_21, id_23, id_60, id_63, id_65, id_126, \
                         id_129, id_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_84 * id_0[k]
                  - f_84 * id_3[k]
                  + f_85 * id_5[k]
                  + f_86 * id_18[k]
                  + f_86 * id_21[k]
                  - f_13 * id_23[k]
                  - f_86 * id_60[k]
                  - f_86 * id_63[k]
                  + f_13 * id_65[k]
                  + f_84 * id_126[k]
                  + f_84 * id_129[k]
                  - f_85 * id_131[k];
    }

#pragma omp simd aligned(id_0, id_2, id_3, id_18, id_20, id_21, id_60, id_62, id_63, id_126, \
                         id_128, id_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_11 * id_2[k]
                  - f_83 * id_20[k]
                  + f_83 * id_62[k]
                  - f_11 * id_128[k];

        g_64[k] = f_87 * id_0[k]
                  - f_87 * id_3[k]
                  - f_88 * id_18[k]
                  + f_88 * id_21[k]
                  + f_88 * id_60[k]
                  - f_88 * id_63[k]
                  - f_87 * id_126[k]
                  + f_87 * id_129[k];
    }
}

}  // namespace simdtrf
