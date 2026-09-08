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


#include "SimdTransformLD.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_ld(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t ld,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.1875 * std::sqrt(2145.0);
    const auto f_1 = 1.3125 * std::sqrt(2145.0);
    const auto f_2 = 0.09375 * std::sqrt(715.0);
    const auto f_3 = 0.1875 * std::sqrt(715.0);
    const auto f_4 = 0.65625 * std::sqrt(715.0);
    const auto f_5 = 1.3125 * std::sqrt(715.0);
    const auto f_6 = 0.09375 * std::sqrt(2145.0);
    const auto f_7 = 0.65625 * std::sqrt(2145.0);
    const auto f_8 = 3.28125 * std::sqrt(2145.0);
    const auto f_9 = 1.96875 * std::sqrt(2145.0);
    const auto f_10 = 0.328125 * std::sqrt(715.0);
    const auto f_11 = 1.640625 * std::sqrt(715.0);
    const auto f_12 = 3.28125 * std::sqrt(715.0);
    const auto f_13 = 0.984375 * std::sqrt(715.0);
    const auto f_14 = 1.96875 * std::sqrt(715.0);
    const auto f_15 = 0.046875 * std::sqrt(715.0);
    const auto f_16 = 0.328125 * std::sqrt(2145.0);
    const auto f_17 = 1.640625 * std::sqrt(2145.0);
    const auto f_18 = 0.984375 * std::sqrt(2145.0);
    const auto f_19 = 0.046875 * std::sqrt(2145.0);
    const auto f_20 = 0.28125 * std::sqrt(286.0);
    const auto f_21 = 0.65625 * std::sqrt(286.0);
    const auto f_22 = 3.9375 * std::sqrt(286.0);
    const auto f_23 = 13.125 * std::sqrt(286.0);
    const auto f_24 = 0.046875 * std::sqrt(858.0);
    const auto f_25 = 0.09375 * std::sqrt(858.0);
    const auto f_26 = 0.109375 * std::sqrt(858.0);
    const auto f_27 = 0.21875 * std::sqrt(858.0);
    const auto f_28 = 0.65625 * std::sqrt(858.0);
    const auto f_29 = 1.3125 * std::sqrt(858.0);
    const auto f_30 = 2.1875 * std::sqrt(858.0);
    const auto f_31 = 4.375 * std::sqrt(858.0);
    const auto f_32 = 0.140625 * std::sqrt(286.0);
    const auto f_33 = 0.328125 * std::sqrt(286.0);
    const auto f_34 = 1.96875 * std::sqrt(286.0);
    const auto f_35 = 6.5625 * std::sqrt(286.0);
    const auto f_36 = 0.46875 * std::sqrt(3003.0);
    const auto f_37 = 1.875 * std::sqrt(3003.0);
    const auto f_38 = 0.84375 * std::sqrt(3003.0);
    const auto f_39 = 3.75 * std::sqrt(3003.0);
    const auto f_40 = 0.09375 * std::sqrt(3003.0);
    const auto f_41 = 0.375 * std::sqrt(3003.0);
    const auto f_42 = 0.234375 * std::sqrt(1001.0);
    const auto f_43 = 0.46875 * std::sqrt(1001.0);
    const auto f_44 = 0.9375 * std::sqrt(1001.0);
    const auto f_45 = 1.875 * std::sqrt(1001.0);
    const auto f_46 = 0.421875 * std::sqrt(1001.0);
    const auto f_47 = 0.84375 * std::sqrt(1001.0);
    const auto f_48 = 3.75 * std::sqrt(1001.0);
    const auto f_49 = 0.046875 * std::sqrt(1001.0);
    const auto f_50 = 0.09375 * std::sqrt(1001.0);
    const auto f_51 = 0.1875 * std::sqrt(1001.0);
    const auto f_52 = 0.375 * std::sqrt(1001.0);
    const auto f_53 = 0.234375 * std::sqrt(3003.0);
    const auto f_54 = 0.9375 * std::sqrt(3003.0);
    const auto f_55 = 0.421875 * std::sqrt(3003.0);
    const auto f_56 = 0.046875 * std::sqrt(3003.0);
    const auto f_57 = 0.1875 * std::sqrt(3003.0);
    const auto f_58 = 0.1875 * std::sqrt(231.0);
    const auto f_59 = 4.5 * std::sqrt(231.0);
    const auto f_60 = 7.5 * std::sqrt(231.0);
    const auto f_61 = 0.09375 * std::sqrt(77.0);
    const auto f_62 = 0.1875 * std::sqrt(77.0);
    const auto f_63 = 2.25 * std::sqrt(77.0);
    const auto f_64 = 4.5 * std::sqrt(77.0);
    const auto f_65 = 3.75 * std::sqrt(77.0);
    const auto f_66 = 7.5 * std::sqrt(77.0);
    const auto f_67 = 0.09375 * std::sqrt(231.0);
    const auto f_68 = 2.25 * std::sqrt(231.0);
    const auto f_69 = 3.75 * std::sqrt(231.0);
    const auto f_70 = 0.84375 * std::sqrt(385.0);
    const auto f_71 = 1.40625 * std::sqrt(385.0);
    const auto f_72 = 5.625 * std::sqrt(385.0);
    const auto f_73 = 0.28125 * std::sqrt(385.0);
    const auto f_74 = 3.75 * std::sqrt(385.0);
    const auto f_75 = 4.5 * std::sqrt(385.0);
    const auto f_76 = 1.875 * std::sqrt(385.0);
    const auto f_77 = 1.5 * std::sqrt(385.0);
    const auto f_78 = 0.140625 * std::sqrt(1155.0);
    const auto f_79 = 0.28125 * std::sqrt(1155.0);
    const auto f_80 = 0.234375 * std::sqrt(1155.0);
    const auto f_81 = 0.46875 * std::sqrt(1155.0);
    const auto f_82 = 0.9375 * std::sqrt(1155.0);
    const auto f_83 = 1.875 * std::sqrt(1155.0);
    const auto f_84 = 0.046875 * std::sqrt(1155.0);
    const auto f_85 = 0.09375 * std::sqrt(1155.0);
    const auto f_86 = 0.625 * std::sqrt(1155.0);
    const auto f_87 = 1.25 * std::sqrt(1155.0);
    const auto f_88 = 0.75 * std::sqrt(1155.0);
    const auto f_89 = 1.5 * std::sqrt(1155.0);
    const auto f_90 = 0.3125 * std::sqrt(1155.0);
    const auto f_91 = 0.25 * std::sqrt(1155.0);
    const auto f_92 = 0.5 * std::sqrt(1155.0);
    const auto f_93 = 0.421875 * std::sqrt(385.0);
    const auto f_94 = 0.703125 * std::sqrt(385.0);
    const auto f_95 = 2.8125 * std::sqrt(385.0);
    const auto f_96 = 0.140625 * std::sqrt(385.0);
    const auto f_97 = 2.25 * std::sqrt(385.0);
    const auto f_98 = 0.9375 * std::sqrt(385.0);
    const auto f_99 = 0.75 * std::sqrt(385.0);
    const auto f_100 = 0.09375 * std::sqrt(210.0);
    const auto f_101 = 0.28125 * std::sqrt(210.0);
    const auto f_102 = 2.8125 * std::sqrt(210.0);
    const auto f_103 = 5.625 * std::sqrt(210.0);
    const auto f_104 = 7.5 * std::sqrt(210.0);
    const auto f_105 = 3.0 * std::sqrt(210.0);
    const auto f_106 = 0.046875 * std::sqrt(70.0);
    const auto f_107 = 0.09375 * std::sqrt(70.0);
    const auto f_108 = 0.140625 * std::sqrt(70.0);
    const auto f_109 = 0.28125 * std::sqrt(70.0);
    const auto f_110 = 1.40625 * std::sqrt(70.0);
    const auto f_111 = 2.8125 * std::sqrt(70.0);
    const auto f_112 = 5.625 * std::sqrt(70.0);
    const auto f_113 = 3.75 * std::sqrt(70.0);
    const auto f_114 = 7.5 * std::sqrt(70.0);
    const auto f_115 = 1.5 * std::sqrt(70.0);
    const auto f_116 = 3.0 * std::sqrt(70.0);
    const auto f_117 = 0.046875 * std::sqrt(210.0);
    const auto f_118 = 0.140625 * std::sqrt(210.0);
    const auto f_119 = 1.40625 * std::sqrt(210.0);
    const auto f_120 = 3.75 * std::sqrt(210.0);
    const auto f_121 = 1.5 * std::sqrt(210.0);
    const auto f_122 = 3.28125 * std::sqrt(3.0);
    const auto f_123 = 9.84375 * std::sqrt(3.0);
    const auto f_124 = 26.25 * std::sqrt(3.0);
    const auto f_125 = 52.5 * std::sqrt(3.0);
    const auto f_126 = 31.5 * std::sqrt(3.0);
    const auto f_127 = 6.0 * std::sqrt(3.0);
    const auto f_128 = 1.640625 * std::sqrt(3.0);
    const auto f_129 = 4.921875 * std::sqrt(3.0);
    const auto f_130 = 13.125 * std::sqrt(3.0);
    const auto f_131 = 15.75 * std::sqrt(3.0);
    const auto f_132 = 3.0 * std::sqrt(3.0);
    const auto f_133 = 0.2734375 * std::sqrt(3.0);
    const auto f_134 = 1.09375 * std::sqrt(3.0);
    const auto f_135 = 8.75 * std::sqrt(3.0);
    const auto f_136 = 14.0 * std::sqrt(3.0);
    const auto f_137 = std::sqrt(3.0);
    const auto f_138 = 0.13671875 * std::sqrt(3.0);
    const auto f_139 = 0.546875 * std::sqrt(3.0);
    const auto f_140 = 4.375 * std::sqrt(3.0);
    const auto f_141 = 0.8203125 * std::sqrt(3.0);
    const auto f_142 = 7.0 * std::sqrt(3.0);
    const auto f_143 = 0.5 * std::sqrt(3.0);
    const auto f_144 = 0.0234375 * std::sqrt(70.0);
    const auto f_145 = 0.703125 * std::sqrt(70.0);
    const auto f_146 = 1.875 * std::sqrt(70.0);
    const auto f_147 = 0.75 * std::sqrt(70.0);
    const auto f_148 = 0.0234375 * std::sqrt(210.0);
    const auto f_149 = 0.703125 * std::sqrt(210.0);
    const auto f_150 = 1.875 * std::sqrt(210.0);
    const auto f_151 = 0.75 * std::sqrt(210.0);
    const auto f_152 = 0.046875 * std::sqrt(231.0);
    const auto f_153 = 1.125 * std::sqrt(231.0);
    const auto f_154 = 0.46875 * std::sqrt(231.0);
    const auto f_155 = 5.625 * std::sqrt(231.0);
    const auto f_156 = 1.875 * std::sqrt(231.0);
    const auto f_157 = 11.25 * std::sqrt(231.0);
    const auto f_158 = 0.0234375 * std::sqrt(77.0);
    const auto f_159 = 0.046875 * std::sqrt(77.0);
    const auto f_160 = 0.5625 * std::sqrt(77.0);
    const auto f_161 = 1.125 * std::sqrt(77.0);
    const auto f_162 = 0.234375 * std::sqrt(77.0);
    const auto f_163 = 0.46875 * std::sqrt(77.0);
    const auto f_164 = 2.8125 * std::sqrt(77.0);
    const auto f_165 = 5.625 * std::sqrt(77.0);
    const auto f_166 = 0.9375 * std::sqrt(77.0);
    const auto f_167 = 1.875 * std::sqrt(77.0);
    const auto f_168 = 11.25 * std::sqrt(77.0);
    const auto f_169 = 0.0234375 * std::sqrt(231.0);
    const auto f_170 = 0.5625 * std::sqrt(231.0);
    const auto f_171 = 0.234375 * std::sqrt(231.0);
    const auto f_172 = 2.8125 * std::sqrt(231.0);
    const auto f_173 = 0.9375 * std::sqrt(231.0);
    const auto f_174 = 0.046875 * std::sqrt(286.0);
    const auto f_175 = 9.84375 * std::sqrt(286.0);
    const auto f_176 = 0.0078125 * std::sqrt(858.0);
    const auto f_177 = 0.015625 * std::sqrt(858.0);
    const auto f_178 = 1.640625 * std::sqrt(858.0);
    const auto f_179 = 3.28125 * std::sqrt(858.0);
    const auto f_180 = 0.0234375 * std::sqrt(286.0);
    const auto f_181 = 4.921875 * std::sqrt(286.0);
    const auto f_182 = 0.0234375 * std::sqrt(2145.0);
    const auto f_183 = 0.01171875 * std::sqrt(715.0);
    const auto f_184 = 0.0234375 * std::sqrt(715.0);
    const auto f_185 = 0.8203125 * std::sqrt(715.0);
    const auto f_186 = 0.01171875 * std::sqrt(2145.0);
    const auto f_187 = 0.8203125 * std::sqrt(2145.0);

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
    auto *g_81 = values + 81 * nvalues;
    auto *g_82 = values + 82 * nvalues;
    auto *g_83 = values + 83 * nvalues;
    auto *g_84 = values + 84 * nvalues;

    const auto *ld_0 = buffer.data(ld + 0);
    const auto *ld_1 = buffer.data(ld + 1);
    const auto *ld_2 = buffer.data(ld + 2);
    const auto *ld_3 = buffer.data(ld + 3);
    const auto *ld_4 = buffer.data(ld + 4);
    const auto *ld_5 = buffer.data(ld + 5);
    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_7 = buffer.data(ld + 7);
    const auto *ld_8 = buffer.data(ld + 8);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_10 = buffer.data(ld + 10);
    const auto *ld_11 = buffer.data(ld + 11);
    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_13 = buffer.data(ld + 13);
    const auto *ld_14 = buffer.data(ld + 14);
    const auto *ld_15 = buffer.data(ld + 15);
    const auto *ld_16 = buffer.data(ld + 16);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_19 = buffer.data(ld + 19);
    const auto *ld_20 = buffer.data(ld + 20);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_22 = buffer.data(ld + 22);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_24 = buffer.data(ld + 24);
    const auto *ld_25 = buffer.data(ld + 25);
    const auto *ld_26 = buffer.data(ld + 26);
    const auto *ld_27 = buffer.data(ld + 27);
    const auto *ld_28 = buffer.data(ld + 28);
    const auto *ld_29 = buffer.data(ld + 29);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_31 = buffer.data(ld + 31);
    const auto *ld_32 = buffer.data(ld + 32);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_34 = buffer.data(ld + 34);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_37 = buffer.data(ld + 37);
    const auto *ld_38 = buffer.data(ld + 38);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_40 = buffer.data(ld + 40);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_42 = buffer.data(ld + 42);
    const auto *ld_43 = buffer.data(ld + 43);
    const auto *ld_44 = buffer.data(ld + 44);
    const auto *ld_45 = buffer.data(ld + 45);
    const auto *ld_46 = buffer.data(ld + 46);
    const auto *ld_47 = buffer.data(ld + 47);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_49 = buffer.data(ld + 49);
    const auto *ld_50 = buffer.data(ld + 50);
    const auto *ld_51 = buffer.data(ld + 51);
    const auto *ld_52 = buffer.data(ld + 52);
    const auto *ld_53 = buffer.data(ld + 53);
    const auto *ld_54 = buffer.data(ld + 54);
    const auto *ld_55 = buffer.data(ld + 55);
    const auto *ld_56 = buffer.data(ld + 56);
    const auto *ld_57 = buffer.data(ld + 57);
    const auto *ld_58 = buffer.data(ld + 58);
    const auto *ld_59 = buffer.data(ld + 59);
    const auto *ld_60 = buffer.data(ld + 60);
    const auto *ld_61 = buffer.data(ld + 61);
    const auto *ld_62 = buffer.data(ld + 62);
    const auto *ld_63 = buffer.data(ld + 63);
    const auto *ld_64 = buffer.data(ld + 64);
    const auto *ld_65 = buffer.data(ld + 65);
    const auto *ld_66 = buffer.data(ld + 66);
    const auto *ld_67 = buffer.data(ld + 67);
    const auto *ld_68 = buffer.data(ld + 68);
    const auto *ld_69 = buffer.data(ld + 69);
    const auto *ld_70 = buffer.data(ld + 70);
    const auto *ld_71 = buffer.data(ld + 71);
    const auto *ld_72 = buffer.data(ld + 72);
    const auto *ld_73 = buffer.data(ld + 73);
    const auto *ld_74 = buffer.data(ld + 74);
    const auto *ld_75 = buffer.data(ld + 75);
    const auto *ld_76 = buffer.data(ld + 76);
    const auto *ld_77 = buffer.data(ld + 77);
    const auto *ld_78 = buffer.data(ld + 78);
    const auto *ld_79 = buffer.data(ld + 79);
    const auto *ld_80 = buffer.data(ld + 80);
    const auto *ld_81 = buffer.data(ld + 81);
    const auto *ld_82 = buffer.data(ld + 82);
    const auto *ld_83 = buffer.data(ld + 83);
    const auto *ld_84 = buffer.data(ld + 84);
    const auto *ld_85 = buffer.data(ld + 85);
    const auto *ld_86 = buffer.data(ld + 86);
    const auto *ld_87 = buffer.data(ld + 87);
    const auto *ld_88 = buffer.data(ld + 88);
    const auto *ld_89 = buffer.data(ld + 89);
    const auto *ld_90 = buffer.data(ld + 90);
    const auto *ld_91 = buffer.data(ld + 91);
    const auto *ld_92 = buffer.data(ld + 92);
    const auto *ld_93 = buffer.data(ld + 93);
    const auto *ld_94 = buffer.data(ld + 94);
    const auto *ld_95 = buffer.data(ld + 95);
    const auto *ld_96 = buffer.data(ld + 96);
    const auto *ld_97 = buffer.data(ld + 97);
    const auto *ld_98 = buffer.data(ld + 98);
    const auto *ld_99 = buffer.data(ld + 99);
    const auto *ld_100 = buffer.data(ld + 100);
    const auto *ld_101 = buffer.data(ld + 101);
    const auto *ld_102 = buffer.data(ld + 102);
    const auto *ld_103 = buffer.data(ld + 103);
    const auto *ld_104 = buffer.data(ld + 104);
    const auto *ld_105 = buffer.data(ld + 105);
    const auto *ld_106 = buffer.data(ld + 106);
    const auto *ld_107 = buffer.data(ld + 107);
    const auto *ld_108 = buffer.data(ld + 108);
    const auto *ld_109 = buffer.data(ld + 109);
    const auto *ld_110 = buffer.data(ld + 110);
    const auto *ld_111 = buffer.data(ld + 111);
    const auto *ld_112 = buffer.data(ld + 112);
    const auto *ld_113 = buffer.data(ld + 113);
    const auto *ld_114 = buffer.data(ld + 114);
    const auto *ld_115 = buffer.data(ld + 115);
    const auto *ld_116 = buffer.data(ld + 116);
    const auto *ld_117 = buffer.data(ld + 117);
    const auto *ld_118 = buffer.data(ld + 118);
    const auto *ld_119 = buffer.data(ld + 119);
    const auto *ld_120 = buffer.data(ld + 120);
    const auto *ld_121 = buffer.data(ld + 121);
    const auto *ld_122 = buffer.data(ld + 122);
    const auto *ld_123 = buffer.data(ld + 123);
    const auto *ld_124 = buffer.data(ld + 124);
    const auto *ld_125 = buffer.data(ld + 125);
    const auto *ld_126 = buffer.data(ld + 126);
    const auto *ld_127 = buffer.data(ld + 127);
    const auto *ld_128 = buffer.data(ld + 128);
    const auto *ld_129 = buffer.data(ld + 129);
    const auto *ld_130 = buffer.data(ld + 130);
    const auto *ld_131 = buffer.data(ld + 131);
    const auto *ld_132 = buffer.data(ld + 132);
    const auto *ld_133 = buffer.data(ld + 133);
    const auto *ld_134 = buffer.data(ld + 134);
    const auto *ld_135 = buffer.data(ld + 135);
    const auto *ld_136 = buffer.data(ld + 136);
    const auto *ld_137 = buffer.data(ld + 137);
    const auto *ld_138 = buffer.data(ld + 138);
    const auto *ld_139 = buffer.data(ld + 139);
    const auto *ld_140 = buffer.data(ld + 140);
    const auto *ld_141 = buffer.data(ld + 141);
    const auto *ld_142 = buffer.data(ld + 142);
    const auto *ld_143 = buffer.data(ld + 143);
    const auto *ld_144 = buffer.data(ld + 144);
    const auto *ld_145 = buffer.data(ld + 145);
    const auto *ld_146 = buffer.data(ld + 146);
    const auto *ld_147 = buffer.data(ld + 147);
    const auto *ld_148 = buffer.data(ld + 148);
    const auto *ld_149 = buffer.data(ld + 149);
    const auto *ld_150 = buffer.data(ld + 150);
    const auto *ld_151 = buffer.data(ld + 151);
    const auto *ld_152 = buffer.data(ld + 152);
    const auto *ld_153 = buffer.data(ld + 153);
    const auto *ld_154 = buffer.data(ld + 154);
    const auto *ld_155 = buffer.data(ld + 155);
    const auto *ld_156 = buffer.data(ld + 156);
    const auto *ld_157 = buffer.data(ld + 157);
    const auto *ld_158 = buffer.data(ld + 158);
    const auto *ld_159 = buffer.data(ld + 159);
    const auto *ld_160 = buffer.data(ld + 160);
    const auto *ld_161 = buffer.data(ld + 161);
    const auto *ld_162 = buffer.data(ld + 162);
    const auto *ld_163 = buffer.data(ld + 163);
    const auto *ld_164 = buffer.data(ld + 164);
    const auto *ld_165 = buffer.data(ld + 165);
    const auto *ld_166 = buffer.data(ld + 166);
    const auto *ld_167 = buffer.data(ld + 167);
    const auto *ld_168 = buffer.data(ld + 168);
    const auto *ld_169 = buffer.data(ld + 169);
    const auto *ld_170 = buffer.data(ld + 170);
    const auto *ld_171 = buffer.data(ld + 171);
    const auto *ld_172 = buffer.data(ld + 172);
    const auto *ld_173 = buffer.data(ld + 173);
    const auto *ld_174 = buffer.data(ld + 174);
    const auto *ld_175 = buffer.data(ld + 175);
    const auto *ld_176 = buffer.data(ld + 176);
    const auto *ld_177 = buffer.data(ld + 177);
    const auto *ld_178 = buffer.data(ld + 178);
    const auto *ld_179 = buffer.data(ld + 179);
    const auto *ld_180 = buffer.data(ld + 180);
    const auto *ld_181 = buffer.data(ld + 181);
    const auto *ld_182 = buffer.data(ld + 182);
    const auto *ld_183 = buffer.data(ld + 183);
    const auto *ld_184 = buffer.data(ld + 184);
    const auto *ld_185 = buffer.data(ld + 185);
    const auto *ld_186 = buffer.data(ld + 186);
    const auto *ld_187 = buffer.data(ld + 187);
    const auto *ld_188 = buffer.data(ld + 188);
    const auto *ld_189 = buffer.data(ld + 189);
    const auto *ld_190 = buffer.data(ld + 190);
    const auto *ld_191 = buffer.data(ld + 191);
    const auto *ld_192 = buffer.data(ld + 192);
    const auto *ld_193 = buffer.data(ld + 193);
    const auto *ld_194 = buffer.data(ld + 194);
    const auto *ld_195 = buffer.data(ld + 195);
    const auto *ld_196 = buffer.data(ld + 196);
    const auto *ld_197 = buffer.data(ld + 197);
    const auto *ld_198 = buffer.data(ld + 198);
    const auto *ld_199 = buffer.data(ld + 199);
    const auto *ld_200 = buffer.data(ld + 200);
    const auto *ld_201 = buffer.data(ld + 201);
    const auto *ld_202 = buffer.data(ld + 202);
    const auto *ld_203 = buffer.data(ld + 203);
    const auto *ld_204 = buffer.data(ld + 204);
    const auto *ld_205 = buffer.data(ld + 205);
    const auto *ld_206 = buffer.data(ld + 206);
    const auto *ld_207 = buffer.data(ld + 207);
    const auto *ld_208 = buffer.data(ld + 208);
    const auto *ld_209 = buffer.data(ld + 209);
    const auto *ld_210 = buffer.data(ld + 210);
    const auto *ld_211 = buffer.data(ld + 211);
    const auto *ld_212 = buffer.data(ld + 212);
    const auto *ld_213 = buffer.data(ld + 213);
    const auto *ld_214 = buffer.data(ld + 214);
    const auto *ld_215 = buffer.data(ld + 215);
    const auto *ld_216 = buffer.data(ld + 216);
    const auto *ld_217 = buffer.data(ld + 217);
    const auto *ld_218 = buffer.data(ld + 218);
    const auto *ld_219 = buffer.data(ld + 219);
    const auto *ld_220 = buffer.data(ld + 220);
    const auto *ld_221 = buffer.data(ld + 221);
    const auto *ld_222 = buffer.data(ld + 222);
    const auto *ld_223 = buffer.data(ld + 223);
    const auto *ld_224 = buffer.data(ld + 224);
    const auto *ld_225 = buffer.data(ld + 225);
    const auto *ld_226 = buffer.data(ld + 226);
    const auto *ld_227 = buffer.data(ld + 227);
    const auto *ld_228 = buffer.data(ld + 228);
    const auto *ld_229 = buffer.data(ld + 229);
    const auto *ld_230 = buffer.data(ld + 230);
    const auto *ld_231 = buffer.data(ld + 231);
    const auto *ld_232 = buffer.data(ld + 232);
    const auto *ld_233 = buffer.data(ld + 233);
    const auto *ld_234 = buffer.data(ld + 234);
    const auto *ld_235 = buffer.data(ld + 235);
    const auto *ld_236 = buffer.data(ld + 236);
    const auto *ld_237 = buffer.data(ld + 237);
    const auto *ld_238 = buffer.data(ld + 238);
    const auto *ld_239 = buffer.data(ld + 239);
    const auto *ld_240 = buffer.data(ld + 240);
    const auto *ld_241 = buffer.data(ld + 241);
    const auto *ld_242 = buffer.data(ld + 242);
    const auto *ld_243 = buffer.data(ld + 243);
    const auto *ld_244 = buffer.data(ld + 244);
    const auto *ld_245 = buffer.data(ld + 245);
    const auto *ld_246 = buffer.data(ld + 246);
    const auto *ld_247 = buffer.data(ld + 247);
    const auto *ld_248 = buffer.data(ld + 248);
    const auto *ld_249 = buffer.data(ld + 249);
    const auto *ld_250 = buffer.data(ld + 250);
    const auto *ld_251 = buffer.data(ld + 251);
    const auto *ld_252 = buffer.data(ld + 252);
    const auto *ld_253 = buffer.data(ld + 253);
    const auto *ld_254 = buffer.data(ld + 254);
    const auto *ld_255 = buffer.data(ld + 255);
    const auto *ld_256 = buffer.data(ld + 256);
    const auto *ld_257 = buffer.data(ld + 257);
    const auto *ld_258 = buffer.data(ld + 258);
    const auto *ld_259 = buffer.data(ld + 259);
    const auto *ld_260 = buffer.data(ld + 260);
    const auto *ld_261 = buffer.data(ld + 261);
    const auto *ld_262 = buffer.data(ld + 262);
    const auto *ld_263 = buffer.data(ld + 263);
    const auto *ld_264 = buffer.data(ld + 264);
    const auto *ld_265 = buffer.data(ld + 265);
    const auto *ld_266 = buffer.data(ld + 266);
    const auto *ld_267 = buffer.data(ld + 267);
    const auto *ld_268 = buffer.data(ld + 268);
    const auto *ld_269 = buffer.data(ld + 269);

#pragma omp simd aligned(ld_7, ld_10, ld_37, ld_40, ld_91, ld_94, ld_169, \
                         ld_172 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * ld_7[k]
                 - f_1 * ld_37[k]
                 + f_1 * ld_91[k]
                 - f_0 * ld_169[k];

        g_1[k] = f_0 * ld_10[k]
                 - f_1 * ld_40[k]
                 + f_1 * ld_94[k]
                 - f_0 * ld_172[k];
    }

#pragma omp simd aligned(ld_6, ld_9, ld_11, ld_36, ld_39, ld_41, ld_90, ld_93, ld_95, ld_168, \
                         ld_171, ld_173 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_2 * ld_6[k]
                 - f_2 * ld_9[k]
                 + f_3 * ld_11[k]
                 + f_4 * ld_36[k]
                 + f_4 * ld_39[k]
                 - f_5 * ld_41[k]
                 - f_4 * ld_90[k]
                 - f_4 * ld_93[k]
                 + f_5 * ld_95[k]
                 + f_2 * ld_168[k]
                 + f_2 * ld_171[k]
                 - f_3 * ld_173[k];
    }

#pragma omp simd aligned(ld_6, ld_8, ld_9, ld_36, ld_38, ld_39, ld_90, ld_92, ld_93, ld_168, \
                         ld_170, ld_171 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_0 * ld_8[k]
                 - f_1 * ld_38[k]
                 + f_1 * ld_92[k]
                 - f_0 * ld_170[k];

        g_4[k] = f_6 * ld_6[k]
                 - f_6 * ld_9[k]
                 - f_7 * ld_36[k]
                 + f_7 * ld_39[k]
                 + f_7 * ld_90[k]
                 - f_7 * ld_93[k]
                 - f_6 * ld_168[k]
                 + f_6 * ld_171[k];
    }

#pragma omp simd aligned(ld_25, ld_28, ld_67, ld_70, ld_133, ld_136, ld_223, \
                         ld_226 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_7 * ld_25[k]
                 - f_8 * ld_67[k]
                 + f_9 * ld_133[k]
                 - f_6 * ld_223[k];

        g_6[k] = f_7 * ld_28[k]
                 - f_8 * ld_70[k]
                 + f_9 * ld_136[k]
                 - f_6 * ld_226[k];
    }

#pragma omp simd aligned(ld_24, ld_27, ld_29, ld_66, ld_69, ld_71, ld_132, ld_135, ld_137, \
                         ld_222, ld_225, ld_227 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_10 * ld_24[k]
                 - f_10 * ld_27[k]
                 + f_4 * ld_29[k]
                 + f_11 * ld_66[k]
                 + f_11 * ld_69[k]
                 - f_12 * ld_71[k]
                 - f_13 * ld_132[k]
                 - f_13 * ld_135[k]
                 + f_14 * ld_137[k]
                 + f_15 * ld_222[k]
                 + f_15 * ld_225[k]
                 - f_2 * ld_227[k];
    }

#pragma omp simd aligned(ld_24, ld_26, ld_27, ld_66, ld_68, ld_69, ld_132, ld_134, ld_135, \
                         ld_222, ld_224, ld_225 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_7 * ld_26[k]
                 - f_8 * ld_68[k]
                 + f_9 * ld_134[k]
                 - f_6 * ld_224[k];

        g_9[k] = f_16 * ld_24[k]
                 - f_16 * ld_27[k]
                 - f_17 * ld_66[k]
                 + f_17 * ld_69[k]
                 + f_18 * ld_132[k]
                 - f_18 * ld_135[k]
                 - f_19 * ld_222[k]
                 + f_19 * ld_225[k];
    }

#pragma omp simd aligned(ld_7, ld_10, ld_37, ld_40, ld_49, ld_52, ld_91, ld_94, ld_103, \
                         ld_106, ld_169, ld_172, ld_181, ld_184 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_20 * ld_7[k]
                  + f_21 * ld_37[k]
                  + f_22 * ld_49[k]
                  + f_21 * ld_91[k]
                  - f_23 * ld_103[k]
                  - f_20 * ld_169[k]
                  + f_22 * ld_181[k];

        g_11[k] = -f_20 * ld_10[k]
                  + f_21 * ld_40[k]
                  + f_22 * ld_52[k]
                  + f_21 * ld_94[k]
                  - f_23 * ld_106[k]
                  - f_20 * ld_172[k]
                  + f_22 * ld_184[k];
    }

#pragma omp simd aligned(ld_6, ld_9, ld_11, ld_36, ld_39, ld_41, ld_48, ld_51, ld_53, ld_90, \
                         ld_93, ld_95, ld_102, ld_105, ld_107, ld_168, ld_171, ld_173, ld_180, \
                         ld_183, ld_185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_24 * ld_6[k]
                  + f_24 * ld_9[k]
                  - f_25 * ld_11[k]
                  - f_26 * ld_36[k]
                  - f_26 * ld_39[k]
                  + f_27 * ld_41[k]
                  - f_28 * ld_48[k]
                  - f_28 * ld_51[k]
                  + f_29 * ld_53[k]
                  - f_26 * ld_90[k]
                  - f_26 * ld_93[k]
                  + f_27 * ld_95[k]
                  + f_30 * ld_102[k]
                  + f_30 * ld_105[k]
                  - f_31 * ld_107[k]
                  + f_24 * ld_168[k]
                  + f_24 * ld_171[k]
                  - f_25 * ld_173[k]
                  - f_28 * ld_180[k]
                  - f_28 * ld_183[k]
                  + f_29 * ld_185[k];
    }

#pragma omp simd aligned(ld_8, ld_38, ld_50, ld_92, ld_104, ld_170, \
                         ld_182 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_20 * ld_8[k]
                  + f_21 * ld_38[k]
                  + f_22 * ld_50[k]
                  + f_21 * ld_92[k]
                  - f_23 * ld_104[k]
                  - f_20 * ld_170[k]
                  + f_22 * ld_182[k];
    }

#pragma omp simd aligned(ld_6, ld_9, ld_36, ld_39, ld_48, ld_51, ld_90, ld_93, ld_102, ld_105, \
                         ld_168, ld_171, ld_180, ld_183 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_32 * ld_6[k]
                  + f_32 * ld_9[k]
                  + f_33 * ld_36[k]
                  - f_33 * ld_39[k]
                  + f_34 * ld_48[k]
                  - f_34 * ld_51[k]
                  + f_33 * ld_90[k]
                  - f_33 * ld_93[k]
                  - f_35 * ld_102[k]
                  + f_35 * ld_105[k]
                  - f_32 * ld_168[k]
                  + f_32 * ld_171[k]
                  + f_34 * ld_180[k]
                  - f_34 * ld_183[k];
    }

#pragma omp simd aligned(ld_25, ld_28, ld_67, ld_70, ld_79, ld_82, ld_133, ld_136, ld_145, \
                         ld_148, ld_223, ld_226, ld_235, ld_238 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_36 * ld_25[k]
                  + f_36 * ld_67[k]
                  + f_37 * ld_79[k]
                  + f_38 * ld_133[k]
                  - f_39 * ld_145[k]
                  - f_40 * ld_223[k]
                  + f_41 * ld_235[k];

        g_16[k] = -f_36 * ld_28[k]
                  + f_36 * ld_70[k]
                  + f_37 * ld_82[k]
                  + f_38 * ld_136[k]
                  - f_39 * ld_148[k]
                  - f_40 * ld_226[k]
                  + f_41 * ld_238[k];
    }

#pragma omp simd aligned(ld_24, ld_27, ld_29, ld_66, ld_69, ld_71, ld_78, ld_81, ld_83, \
                         ld_132, ld_135, ld_137, ld_144, ld_147, ld_149, ld_222, ld_225, \
                         ld_227, ld_234, ld_237, ld_239 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_42 * ld_24[k]
                  + f_42 * ld_27[k]
                  - f_43 * ld_29[k]
                  - f_42 * ld_66[k]
                  - f_42 * ld_69[k]
                  + f_43 * ld_71[k]
                  - f_44 * ld_78[k]
                  - f_44 * ld_81[k]
                  + f_45 * ld_83[k]
                  - f_46 * ld_132[k]
                  - f_46 * ld_135[k]
                  + f_47 * ld_137[k]
                  + f_45 * ld_144[k]
                  + f_45 * ld_147[k]
                  - f_48 * ld_149[k]
                  + f_49 * ld_222[k]
                  + f_49 * ld_225[k]
                  - f_50 * ld_227[k]
                  - f_51 * ld_234[k]
                  - f_51 * ld_237[k]
                  + f_52 * ld_239[k];
    }

#pragma omp simd aligned(ld_26, ld_68, ld_80, ld_134, ld_146, ld_224, \
                         ld_236 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_36 * ld_26[k]
                  + f_36 * ld_68[k]
                  + f_37 * ld_80[k]
                  + f_38 * ld_134[k]
                  - f_39 * ld_146[k]
                  - f_40 * ld_224[k]
                  + f_41 * ld_236[k];
    }

#pragma omp simd aligned(ld_24, ld_27, ld_66, ld_69, ld_78, ld_81, ld_132, ld_135, ld_144, \
                         ld_147, ld_222, ld_225, ld_234, ld_237 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_53 * ld_24[k]
                  + f_53 * ld_27[k]
                  + f_53 * ld_66[k]
                  - f_53 * ld_69[k]
                  + f_54 * ld_78[k]
                  - f_54 * ld_81[k]
                  + f_55 * ld_132[k]
                  - f_55 * ld_135[k]
                  - f_37 * ld_144[k]
                  + f_37 * ld_147[k]
                  - f_56 * ld_222[k]
                  + f_56 * ld_225[k]
                  + f_57 * ld_234[k]
                  - f_57 * ld_237[k];
    }

#pragma omp simd aligned(ld_7, ld_37, ld_49, ld_91, ld_115, ld_169, ld_181, \
                         ld_193 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_58 * ld_7[k]
                  + f_58 * ld_37[k]
                  - f_59 * ld_49[k]
                  - f_58 * ld_91[k]
                  + f_60 * ld_115[k]
                  - f_58 * ld_169[k]
                  + f_59 * ld_181[k]
                  - f_60 * ld_193[k];
    }

#pragma omp simd aligned(ld_10, ld_40, ld_52, ld_94, ld_118, ld_172, ld_184, \
                         ld_196 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_58 * ld_10[k]
                  + f_58 * ld_40[k]
                  - f_59 * ld_52[k]
                  - f_58 * ld_94[k]
                  + f_60 * ld_118[k]
                  - f_58 * ld_172[k]
                  + f_59 * ld_184[k]
                  - f_60 * ld_196[k];
    }

#pragma omp simd aligned(ld_6, ld_9, ld_11, ld_36, ld_39, ld_41, ld_48, ld_51, ld_53, ld_90, \
                         ld_93, ld_95, ld_114, ld_117, ld_119, ld_168, ld_171, ld_173, ld_180, \
                         ld_183, ld_185, ld_192, ld_195, ld_197 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_61 * ld_6[k]
                  - f_61 * ld_9[k]
                  + f_62 * ld_11[k]
                  - f_61 * ld_36[k]
                  - f_61 * ld_39[k]
                  + f_62 * ld_41[k]
                  + f_63 * ld_48[k]
                  + f_63 * ld_51[k]
                  - f_64 * ld_53[k]
                  + f_61 * ld_90[k]
                  + f_61 * ld_93[k]
                  - f_62 * ld_95[k]
                  - f_65 * ld_114[k]
                  - f_65 * ld_117[k]
                  + f_66 * ld_119[k]
                  + f_61 * ld_168[k]
                  + f_61 * ld_171[k]
                  - f_62 * ld_173[k]
                  - f_63 * ld_180[k]
                  - f_63 * ld_183[k]
                  + f_64 * ld_185[k]
                  + f_65 * ld_192[k]
                  + f_65 * ld_195[k]
                  - f_66 * ld_197[k];
    }

#pragma omp simd aligned(ld_8, ld_38, ld_50, ld_92, ld_116, ld_170, ld_182, \
                         ld_194 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_58 * ld_8[k]
                  + f_58 * ld_38[k]
                  - f_59 * ld_50[k]
                  - f_58 * ld_92[k]
                  + f_60 * ld_116[k]
                  - f_58 * ld_170[k]
                  + f_59 * ld_182[k]
                  - f_60 * ld_194[k];
    }

#pragma omp simd aligned(ld_6, ld_9, ld_36, ld_39, ld_48, ld_51, ld_90, ld_93, ld_114, ld_117, \
                         ld_168, ld_171, ld_180, ld_183, ld_192, \
                         ld_195 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_67 * ld_6[k]
                  - f_67 * ld_9[k]
                  + f_67 * ld_36[k]
                  - f_67 * ld_39[k]
                  - f_68 * ld_48[k]
                  + f_68 * ld_51[k]
                  - f_67 * ld_90[k]
                  + f_67 * ld_93[k]
                  + f_69 * ld_114[k]
                  - f_69 * ld_117[k]
                  - f_67 * ld_168[k]
                  + f_67 * ld_171[k]
                  + f_68 * ld_180[k]
                  - f_68 * ld_183[k]
                  - f_69 * ld_192[k]
                  + f_69 * ld_195[k];
    }

#pragma omp simd aligned(ld_25, ld_67, ld_79, ld_133, ld_145, ld_157, ld_223, ld_235, \
                         ld_247 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_70 * ld_25[k]
                  + f_71 * ld_67[k]
                  - f_72 * ld_79[k]
                  + f_73 * ld_133[k]
                  - f_74 * ld_145[k]
                  + f_75 * ld_157[k]
                  - f_73 * ld_223[k]
                  + f_76 * ld_235[k]
                  - f_77 * ld_247[k];
    }

#pragma omp simd aligned(ld_28, ld_70, ld_82, ld_136, ld_148, ld_160, ld_226, ld_238, \
                         ld_250 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = f_70 * ld_28[k]
                  + f_71 * ld_70[k]
                  - f_72 * ld_82[k]
                  + f_73 * ld_136[k]
                  - f_74 * ld_148[k]
                  + f_75 * ld_160[k]
                  - f_73 * ld_226[k]
                  + f_76 * ld_238[k]
                  - f_77 * ld_250[k];
    }

#pragma omp simd aligned(ld_24, ld_27, ld_29, ld_66, ld_69, ld_71, ld_78, ld_81, ld_83, \
                         ld_132, ld_135, ld_137, ld_144, ld_147, ld_149, ld_156, ld_159, \
                         ld_161, ld_222, ld_225, ld_227, ld_234, ld_237, ld_239, ld_246, \
                         ld_249, ld_251 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_78 * ld_24[k]
                  - f_78 * ld_27[k]
                  + f_79 * ld_29[k]
                  - f_80 * ld_66[k]
                  - f_80 * ld_69[k]
                  + f_81 * ld_71[k]
                  + f_82 * ld_78[k]
                  + f_82 * ld_81[k]
                  - f_83 * ld_83[k]
                  - f_84 * ld_132[k]
                  - f_84 * ld_135[k]
                  + f_85 * ld_137[k]
                  + f_86 * ld_144[k]
                  + f_86 * ld_147[k]
                  - f_87 * ld_149[k]
                  - f_88 * ld_156[k]
                  - f_88 * ld_159[k]
                  + f_89 * ld_161[k]
                  + f_84 * ld_222[k]
                  + f_84 * ld_225[k]
                  - f_85 * ld_227[k]
                  - f_90 * ld_234[k]
                  - f_90 * ld_237[k]
                  + f_86 * ld_239[k]
                  + f_91 * ld_246[k]
                  + f_91 * ld_249[k]
                  - f_92 * ld_251[k];
    }

#pragma omp simd aligned(ld_26, ld_68, ld_80, ld_134, ld_146, ld_158, ld_224, ld_236, \
                         ld_248 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_70 * ld_26[k]
                  + f_71 * ld_68[k]
                  - f_72 * ld_80[k]
                  + f_73 * ld_134[k]
                  - f_74 * ld_146[k]
                  + f_75 * ld_158[k]
                  - f_73 * ld_224[k]
                  + f_76 * ld_236[k]
                  - f_77 * ld_248[k];
    }

#pragma omp simd aligned(ld_24, ld_27, ld_66, ld_69, ld_78, ld_81, ld_132, ld_135, ld_144, \
                         ld_147, ld_156, ld_159, ld_222, ld_225, ld_234, ld_237, ld_246, \
                         ld_249 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_93 * ld_24[k]
                  - f_93 * ld_27[k]
                  + f_94 * ld_66[k]
                  - f_94 * ld_69[k]
                  - f_95 * ld_78[k]
                  + f_95 * ld_81[k]
                  + f_96 * ld_132[k]
                  - f_96 * ld_135[k]
                  - f_76 * ld_144[k]
                  + f_76 * ld_147[k]
                  + f_97 * ld_156[k]
                  - f_97 * ld_159[k]
                  - f_96 * ld_222[k]
                  + f_96 * ld_225[k]
                  + f_98 * ld_234[k]
                  - f_98 * ld_237[k]
                  - f_99 * ld_246[k]
                  + f_99 * ld_249[k];
    }

#pragma omp simd aligned(ld_7, ld_37, ld_49, ld_91, ld_103, ld_115, ld_169, ld_181, ld_193, \
                         ld_205 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_100 * ld_7[k]
                  - f_101 * ld_37[k]
                  + f_102 * ld_49[k]
                  - f_101 * ld_91[k]
                  + f_103 * ld_103[k]
                  - f_104 * ld_115[k]
                  - f_100 * ld_169[k]
                  + f_102 * ld_181[k]
                  - f_104 * ld_193[k]
                  + f_105 * ld_205[k];
    }

#pragma omp simd aligned(ld_10, ld_40, ld_52, ld_94, ld_106, ld_118, ld_172, ld_184, ld_196, \
                         ld_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_100 * ld_10[k]
                  - f_101 * ld_40[k]
                  + f_102 * ld_52[k]
                  - f_101 * ld_94[k]
                  + f_103 * ld_106[k]
                  - f_104 * ld_118[k]
                  - f_100 * ld_172[k]
                  + f_102 * ld_184[k]
                  - f_104 * ld_196[k]
                  + f_105 * ld_208[k];
    }

#pragma omp simd aligned(ld_6, ld_9, ld_11, ld_36, ld_39, ld_41, ld_48, ld_51, ld_53, ld_90, \
                         ld_93, ld_95, ld_102, ld_105, ld_107, ld_114, ld_117, ld_119, ld_168, \
                         ld_171, ld_173, ld_180, ld_183, ld_185, ld_192, ld_195, ld_197, \
                         ld_204, ld_207, ld_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_106 * ld_6[k]
                  + f_106 * ld_9[k]
                  - f_107 * ld_11[k]
                  + f_108 * ld_36[k]
                  + f_108 * ld_39[k]
                  - f_109 * ld_41[k]
                  - f_110 * ld_48[k]
                  - f_110 * ld_51[k]
                  + f_111 * ld_53[k]
                  + f_108 * ld_90[k]
                  + f_108 * ld_93[k]
                  - f_109 * ld_95[k]
                  - f_111 * ld_102[k]
                  - f_111 * ld_105[k]
                  + f_112 * ld_107[k]
                  + f_113 * ld_114[k]
                  + f_113 * ld_117[k]
                  - f_114 * ld_119[k]
                  + f_106 * ld_168[k]
                  + f_106 * ld_171[k]
                  - f_107 * ld_173[k]
                  - f_110 * ld_180[k]
                  - f_110 * ld_183[k]
                  + f_111 * ld_185[k]
                  + f_113 * ld_192[k]
                  + f_113 * ld_195[k]
                  - f_114 * ld_197[k]
                  - f_115 * ld_204[k]
                  - f_115 * ld_207[k]
                  + f_116 * ld_209[k];
    }

#pragma omp simd aligned(ld_8, ld_38, ld_50, ld_92, ld_104, ld_116, ld_170, ld_182, ld_194, \
                         ld_206 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_100 * ld_8[k]
                  - f_101 * ld_38[k]
                  + f_102 * ld_50[k]
                  - f_101 * ld_92[k]
                  + f_103 * ld_104[k]
                  - f_104 * ld_116[k]
                  - f_100 * ld_170[k]
                  + f_102 * ld_182[k]
                  - f_104 * ld_194[k]
                  + f_105 * ld_206[k];
    }

#pragma omp simd aligned(ld_6, ld_9, ld_36, ld_39, ld_48, ld_51, ld_90, ld_93, ld_102, ld_105, \
                         ld_114, ld_117, ld_168, ld_171, ld_180, ld_183, ld_192, ld_195, \
                         ld_204, ld_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_117 * ld_6[k]
                  + f_117 * ld_9[k]
                  - f_118 * ld_36[k]
                  + f_118 * ld_39[k]
                  + f_119 * ld_48[k]
                  - f_119 * ld_51[k]
                  - f_118 * ld_90[k]
                  + f_118 * ld_93[k]
                  + f_102 * ld_102[k]
                  - f_102 * ld_105[k]
                  - f_120 * ld_114[k]
                  + f_120 * ld_117[k]
                  - f_117 * ld_168[k]
                  + f_117 * ld_171[k]
                  + f_119 * ld_180[k]
                  - f_119 * ld_183[k]
                  - f_120 * ld_192[k]
                  + f_120 * ld_195[k]
                  + f_121 * ld_204[k]
                  - f_121 * ld_207[k];
    }

#pragma omp simd aligned(ld_25, ld_67, ld_79, ld_133, ld_145, ld_157, ld_223, ld_235, ld_247, \
                         ld_259 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_122 * ld_25[k]
                  - f_123 * ld_67[k]
                  + f_124 * ld_79[k]
                  - f_123 * ld_133[k]
                  + f_125 * ld_145[k]
                  - f_126 * ld_157[k]
                  - f_122 * ld_223[k]
                  + f_124 * ld_235[k]
                  - f_126 * ld_247[k]
                  + f_127 * ld_259[k];
    }

#pragma omp simd aligned(ld_28, ld_70, ld_82, ld_136, ld_148, ld_160, ld_226, ld_238, ld_250, \
                         ld_262 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = -f_122 * ld_28[k]
                  - f_123 * ld_70[k]
                  + f_124 * ld_82[k]
                  - f_123 * ld_136[k]
                  + f_125 * ld_148[k]
                  - f_126 * ld_160[k]
                  - f_122 * ld_226[k]
                  + f_124 * ld_238[k]
                  - f_126 * ld_250[k]
                  + f_127 * ld_262[k];
    }

#pragma omp simd aligned(ld_24, ld_27, ld_29, ld_66, ld_69, ld_71, ld_78, ld_81, ld_83, \
                         ld_132, ld_135, ld_137, ld_144, ld_147, ld_149, ld_156, ld_159, \
                         ld_161, ld_222, ld_225, ld_227, ld_234, ld_237, ld_239, ld_246, \
                         ld_249, ld_251, ld_258, ld_261, ld_263 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = 1.640625 * ld_24[k]
                  + 1.640625 * ld_27[k]
                  - 3.28125 * ld_29[k]
                  + 4.921875 * ld_66[k]
                  + 4.921875 * ld_69[k]
                  - 9.84375 * ld_71[k]
                  - 13.125 * ld_78[k]
                  - 13.125 * ld_81[k]
                  + 26.25 * ld_83[k]
                  + 4.921875 * ld_132[k]
                  + 4.921875 * ld_135[k]
                  - 9.84375 * ld_137[k]
                  - 26.25 * ld_144[k]
                  - 26.25 * ld_147[k]
                  + 52.5 * ld_149[k]
                  + 15.75 * ld_156[k]
                  + 15.75 * ld_159[k]
                  - 31.5 * ld_161[k]
                  + 1.640625 * ld_222[k]
                  + 1.640625 * ld_225[k]
                  - 3.28125 * ld_227[k]
                  - 13.125 * ld_234[k]
                  - 13.125 * ld_237[k]
                  + 26.25 * ld_239[k]
                  + 15.75 * ld_246[k]
                  + 15.75 * ld_249[k]
                  - 31.5 * ld_251[k]
                  - 3.0 * ld_258[k]
                  - 3.0 * ld_261[k]
                  + 6.0 * ld_263[k];
    }

#pragma omp simd aligned(ld_26, ld_68, ld_80, ld_134, ld_146, ld_158, ld_224, ld_236, ld_248, \
                         ld_260 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_122 * ld_26[k]
                  - f_123 * ld_68[k]
                  + f_124 * ld_80[k]
                  - f_123 * ld_134[k]
                  + f_125 * ld_146[k]
                  - f_126 * ld_158[k]
                  - f_122 * ld_224[k]
                  + f_124 * ld_236[k]
                  - f_126 * ld_248[k]
                  + f_127 * ld_260[k];
    }

#pragma omp simd aligned(ld_24, ld_27, ld_66, ld_69, ld_78, ld_81, ld_132, ld_135, ld_144, \
                         ld_147, ld_156, ld_159, ld_222, ld_225, ld_234, ld_237, ld_246, \
                         ld_249, ld_258, ld_261 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_128 * ld_24[k]
                  + f_128 * ld_27[k]
                  - f_129 * ld_66[k]
                  + f_129 * ld_69[k]
                  + f_130 * ld_78[k]
                  - f_130 * ld_81[k]
                  - f_129 * ld_132[k]
                  + f_129 * ld_135[k]
                  + f_124 * ld_144[k]
                  - f_124 * ld_147[k]
                  - f_131 * ld_156[k]
                  + f_131 * ld_159[k]
                  - f_128 * ld_222[k]
                  + f_128 * ld_225[k]
                  + f_130 * ld_234[k]
                  - f_130 * ld_237[k]
                  - f_131 * ld_246[k]
                  + f_131 * ld_249[k]
                  + f_132 * ld_258[k]
                  - f_132 * ld_261[k];
    }

#pragma omp simd aligned(ld_1, ld_19, ld_31, ld_61, ld_73, ld_85, ld_127, ld_139, ld_151, \
                         ld_163, ld_217, ld_229, ld_241, ld_253, \
                         ld_265 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_133 * ld_1[k]
                  + f_134 * ld_19[k]
                  - f_135 * ld_31[k]
                  + f_128 * ld_61[k]
                  - f_124 * ld_73[k]
                  + f_124 * ld_85[k]
                  + f_134 * ld_127[k]
                  - f_124 * ld_139[k]
                  + f_125 * ld_151[k]
                  - f_136 * ld_163[k]
                  + f_133 * ld_217[k]
                  - f_135 * ld_229[k]
                  + f_124 * ld_241[k]
                  - f_136 * ld_253[k]
                  + f_137 * ld_265[k];
    }

#pragma omp simd aligned(ld_4, ld_22, ld_34, ld_64, ld_76, ld_88, ld_130, ld_142, ld_154, \
                         ld_166, ld_220, ld_232, ld_244, ld_256, \
                         ld_268 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_133 * ld_4[k]
                  + f_134 * ld_22[k]
                  - f_135 * ld_34[k]
                  + f_128 * ld_64[k]
                  - f_124 * ld_76[k]
                  + f_124 * ld_88[k]
                  + f_134 * ld_130[k]
                  - f_124 * ld_142[k]
                  + f_125 * ld_154[k]
                  - f_136 * ld_166[k]
                  + f_133 * ld_220[k]
                  - f_135 * ld_232[k]
                  + f_124 * ld_244[k]
                  - f_136 * ld_256[k]
                  + f_137 * ld_268[k];
    }

#pragma omp simd aligned(ld_0, ld_3, ld_5, ld_18, ld_21, ld_23, ld_30, ld_33, ld_35, ld_60, \
                         ld_63, ld_65, ld_72, ld_75, ld_77, ld_84, ld_87, ld_89, ld_126, \
                         ld_129, ld_131, ld_138, ld_141, ld_143, ld_150, ld_153, ld_155, \
                         ld_162, ld_165, ld_167, ld_216, ld_219, ld_221, ld_228, ld_231, \
                         ld_233, ld_240, ld_243, ld_245, ld_252, ld_255, ld_257, ld_264, \
                         ld_267, ld_269 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -0.13671875 * ld_0[k]
                  - 0.13671875 * ld_3[k]
                  + 0.2734375 * ld_5[k]
                  - 0.546875 * ld_18[k]
                  - 0.546875 * ld_21[k]
                  + 1.09375 * ld_23[k]
                  + 4.375 * ld_30[k]
                  + 4.375 * ld_33[k]
                  - 8.75 * ld_35[k]
                  - 0.8203125 * ld_60[k]
                  - 0.8203125 * ld_63[k]
                  + 1.640625 * ld_65[k]
                  + 13.125 * ld_72[k]
                  + 13.125 * ld_75[k]
                  - 26.25 * ld_77[k]
                  - 13.125 * ld_84[k]
                  - 13.125 * ld_87[k]
                  + 26.25 * ld_89[k]
                  - 0.546875 * ld_126[k]
                  - 0.546875 * ld_129[k]
                  + 1.09375 * ld_131[k]
                  + 13.125 * ld_138[k]
                  + 13.125 * ld_141[k]
                  - 26.25 * ld_143[k]
                  - 26.25 * ld_150[k]
                  - 26.25 * ld_153[k]
                  + 52.5 * ld_155[k]
                  + 7.0 * ld_162[k]
                  + 7.0 * ld_165[k]
                  - 14.0 * ld_167[k]
                  - 0.13671875 * ld_216[k]
                  - 0.13671875 * ld_219[k]
                  + 0.2734375 * ld_221[k]
                  + 4.375 * ld_228[k]
                  + 4.375 * ld_231[k]
                  - 8.75 * ld_233[k]
                  - 13.125 * ld_240[k]
                  - 13.125 * ld_243[k]
                  + 26.25 * ld_245[k]
                  + 7.0 * ld_252[k]
                  + 7.0 * ld_255[k]
                  - 14.0 * ld_257[k]
                  - 0.5 * ld_264[k]
                  - 0.5 * ld_267[k]
                  + ld_269[k];
    }

#pragma omp simd aligned(ld_2, ld_20, ld_32, ld_62, ld_74, ld_86, ld_128, ld_140, ld_152, \
                         ld_164, ld_218, ld_230, ld_242, ld_254, \
                         ld_266 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_133 * ld_2[k]
                  + f_134 * ld_20[k]
                  - f_135 * ld_32[k]
                  + f_128 * ld_62[k]
                  - f_124 * ld_74[k]
                  + f_124 * ld_86[k]
                  + f_134 * ld_128[k]
                  - f_124 * ld_140[k]
                  + f_125 * ld_152[k]
                  - f_136 * ld_164[k]
                  + f_133 * ld_218[k]
                  - f_135 * ld_230[k]
                  + f_124 * ld_242[k]
                  - f_136 * ld_254[k]
                  + f_137 * ld_266[k];
    }

#pragma omp simd aligned(ld_0, ld_3, ld_18, ld_21, ld_30, ld_33, ld_60, ld_63, ld_72, ld_75, \
                         ld_84, ld_87, ld_126, ld_129, ld_138, ld_141, ld_150, ld_153, ld_162, \
                         ld_165, ld_216, ld_219, ld_228, ld_231, ld_240, ld_243, ld_252, \
                         ld_255, ld_264, ld_267 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_138 * ld_0[k]
                  - f_138 * ld_3[k]
                  + f_139 * ld_18[k]
                  - f_139 * ld_21[k]
                  - f_140 * ld_30[k]
                  + f_140 * ld_33[k]
                  + f_141 * ld_60[k]
                  - f_141 * ld_63[k]
                  - f_130 * ld_72[k]
                  + f_130 * ld_75[k]
                  + f_130 * ld_84[k]
                  - f_130 * ld_87[k]
                  + f_139 * ld_126[k]
                  - f_139 * ld_129[k]
                  - f_130 * ld_138[k]
                  + f_130 * ld_141[k]
                  + f_124 * ld_150[k]
                  - f_124 * ld_153[k]
                  - f_142 * ld_162[k]
                  + f_142 * ld_165[k]
                  + f_138 * ld_216[k]
                  - f_138 * ld_219[k]
                  - f_140 * ld_228[k]
                  + f_140 * ld_231[k]
                  + f_130 * ld_240[k]
                  - f_130 * ld_243[k]
                  - f_142 * ld_252[k]
                  + f_142 * ld_255[k]
                  + f_143 * ld_264[k]
                  - f_143 * ld_267[k];
    }

#pragma omp simd aligned(ld_13, ld_43, ld_55, ld_97, ld_109, ld_121, ld_175, ld_187, ld_199, \
                         ld_211 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_122 * ld_13[k]
                  - f_123 * ld_43[k]
                  + f_124 * ld_55[k]
                  - f_123 * ld_97[k]
                  + f_125 * ld_109[k]
                  - f_126 * ld_121[k]
                  - f_122 * ld_175[k]
                  + f_124 * ld_187[k]
                  - f_126 * ld_199[k]
                  + f_127 * ld_211[k];
    }

#pragma omp simd aligned(ld_16, ld_46, ld_58, ld_100, ld_112, ld_124, ld_178, ld_190, ld_202, \
                         ld_214 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_122 * ld_16[k]
                  - f_123 * ld_46[k]
                  + f_124 * ld_58[k]
                  - f_123 * ld_100[k]
                  + f_125 * ld_112[k]
                  - f_126 * ld_124[k]
                  - f_122 * ld_178[k]
                  + f_124 * ld_190[k]
                  - f_126 * ld_202[k]
                  + f_127 * ld_214[k];
    }

#pragma omp simd aligned(ld_12, ld_15, ld_17, ld_42, ld_45, ld_47, ld_54, ld_57, ld_59, ld_96, \
                         ld_99, ld_101, ld_108, ld_111, ld_113, ld_120, ld_123, ld_125, \
                         ld_174, ld_177, ld_179, ld_186, ld_189, ld_191, ld_198, ld_201, \
                         ld_203, ld_210, ld_213, ld_215 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = 1.640625 * ld_12[k]
                  + 1.640625 * ld_15[k]
                  - 3.28125 * ld_17[k]
                  + 4.921875 * ld_42[k]
                  + 4.921875 * ld_45[k]
                  - 9.84375 * ld_47[k]
                  - 13.125 * ld_54[k]
                  - 13.125 * ld_57[k]
                  + 26.25 * ld_59[k]
                  + 4.921875 * ld_96[k]
                  + 4.921875 * ld_99[k]
                  - 9.84375 * ld_101[k]
                  - 26.25 * ld_108[k]
                  - 26.25 * ld_111[k]
                  + 52.5 * ld_113[k]
                  + 15.75 * ld_120[k]
                  + 15.75 * ld_123[k]
                  - 31.5 * ld_125[k]
                  + 1.640625 * ld_174[k]
                  + 1.640625 * ld_177[k]
                  - 3.28125 * ld_179[k]
                  - 13.125 * ld_186[k]
                  - 13.125 * ld_189[k]
                  + 26.25 * ld_191[k]
                  + 15.75 * ld_198[k]
                  + 15.75 * ld_201[k]
                  - 31.5 * ld_203[k]
                  - 3.0 * ld_210[k]
                  - 3.0 * ld_213[k]
                  + 6.0 * ld_215[k];
    }

#pragma omp simd aligned(ld_14, ld_44, ld_56, ld_98, ld_110, ld_122, ld_176, ld_188, ld_200, \
                         ld_212 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_122 * ld_14[k]
                  - f_123 * ld_44[k]
                  + f_124 * ld_56[k]
                  - f_123 * ld_98[k]
                  + f_125 * ld_110[k]
                  - f_126 * ld_122[k]
                  - f_122 * ld_176[k]
                  + f_124 * ld_188[k]
                  - f_126 * ld_200[k]
                  + f_127 * ld_212[k];
    }

#pragma omp simd aligned(ld_12, ld_15, ld_42, ld_45, ld_54, ld_57, ld_96, ld_99, ld_108, \
                         ld_111, ld_120, ld_123, ld_174, ld_177, ld_186, ld_189, ld_198, \
                         ld_201, ld_210, ld_213 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_128 * ld_12[k]
                  + f_128 * ld_15[k]
                  - f_129 * ld_42[k]
                  + f_129 * ld_45[k]
                  + f_130 * ld_54[k]
                  - f_130 * ld_57[k]
                  - f_129 * ld_96[k]
                  + f_129 * ld_99[k]
                  + f_124 * ld_108[k]
                  - f_124 * ld_111[k]
                  - f_131 * ld_120[k]
                  + f_131 * ld_123[k]
                  - f_128 * ld_174[k]
                  + f_128 * ld_177[k]
                  + f_130 * ld_186[k]
                  - f_130 * ld_189[k]
                  - f_131 * ld_198[k]
                  + f_131 * ld_201[k]
                  + f_132 * ld_210[k]
                  - f_132 * ld_213[k];
    }

#pragma omp simd aligned(ld_1, ld_19, ld_31, ld_73, ld_85, ld_127, ld_139, ld_163, ld_217, \
                         ld_229, ld_241, ld_253 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_117 * ld_1[k]
                  - f_100 * ld_19[k]
                  + f_119 * ld_31[k]
                  + f_119 * ld_73[k]
                  - f_120 * ld_85[k]
                  + f_100 * ld_127[k]
                  - f_119 * ld_139[k]
                  + f_121 * ld_163[k]
                  + f_117 * ld_217[k]
                  - f_119 * ld_229[k]
                  + f_120 * ld_241[k]
                  - f_121 * ld_253[k];
    }

#pragma omp simd aligned(ld_4, ld_22, ld_34, ld_76, ld_88, ld_130, ld_142, ld_166, ld_220, \
                         ld_232, ld_244, ld_256 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_117 * ld_4[k]
                  - f_100 * ld_22[k]
                  + f_119 * ld_34[k]
                  + f_119 * ld_76[k]
                  - f_120 * ld_88[k]
                  + f_100 * ld_130[k]
                  - f_119 * ld_142[k]
                  + f_121 * ld_166[k]
                  + f_117 * ld_220[k]
                  - f_119 * ld_232[k]
                  + f_120 * ld_244[k]
                  - f_121 * ld_256[k];
    }

#pragma omp simd aligned(ld_0, ld_3, ld_5, ld_18, ld_21, ld_23, ld_30, ld_33, ld_35, ld_72, \
                         ld_75, ld_77, ld_84, ld_87, ld_89, ld_126, ld_129, ld_131, ld_138, \
                         ld_141, ld_143, ld_162, ld_165, ld_167, ld_216, ld_219, ld_221, \
                         ld_228, ld_231, ld_233, ld_240, ld_243, ld_245, ld_252, ld_255, \
                         ld_257 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_144 * ld_0[k]
                  + f_144 * ld_3[k]
                  - f_106 * ld_5[k]
                  + f_106 * ld_18[k]
                  + f_106 * ld_21[k]
                  - f_107 * ld_23[k]
                  - f_145 * ld_30[k]
                  - f_145 * ld_33[k]
                  + f_110 * ld_35[k]
                  - f_145 * ld_72[k]
                  - f_145 * ld_75[k]
                  + f_110 * ld_77[k]
                  + f_146 * ld_84[k]
                  + f_146 * ld_87[k]
                  - f_113 * ld_89[k]
                  - f_106 * ld_126[k]
                  - f_106 * ld_129[k]
                  + f_107 * ld_131[k]
                  + f_145 * ld_138[k]
                  + f_145 * ld_141[k]
                  - f_110 * ld_143[k]
                  - f_147 * ld_162[k]
                  - f_147 * ld_165[k]
                  + f_115 * ld_167[k]
                  - f_144 * ld_216[k]
                  - f_144 * ld_219[k]
                  + f_106 * ld_221[k]
                  + f_145 * ld_228[k]
                  + f_145 * ld_231[k]
                  - f_110 * ld_233[k]
                  - f_146 * ld_240[k]
                  - f_146 * ld_243[k]
                  + f_113 * ld_245[k]
                  + f_147 * ld_252[k]
                  + f_147 * ld_255[k]
                  - f_115 * ld_257[k];
    }

#pragma omp simd aligned(ld_2, ld_20, ld_32, ld_74, ld_86, ld_128, ld_140, ld_164, ld_218, \
                         ld_230, ld_242, ld_254 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = -f_117 * ld_2[k]
                  - f_100 * ld_20[k]
                  + f_119 * ld_32[k]
                  + f_119 * ld_74[k]
                  - f_120 * ld_86[k]
                  + f_100 * ld_128[k]
                  - f_119 * ld_140[k]
                  + f_121 * ld_164[k]
                  + f_117 * ld_218[k]
                  - f_119 * ld_230[k]
                  + f_120 * ld_242[k]
                  - f_121 * ld_254[k];
    }

#pragma omp simd aligned(ld_0, ld_3, ld_18, ld_21, ld_30, ld_33, ld_72, ld_75, ld_84, ld_87, \
                         ld_126, ld_129, ld_138, ld_141, ld_162, ld_165, ld_216, ld_219, \
                         ld_228, ld_231, ld_240, ld_243, ld_252, \
                         ld_255 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_148 * ld_0[k]
                  + f_148 * ld_3[k]
                  - f_117 * ld_18[k]
                  + f_117 * ld_21[k]
                  + f_149 * ld_30[k]
                  - f_149 * ld_33[k]
                  + f_149 * ld_72[k]
                  - f_149 * ld_75[k]
                  - f_150 * ld_84[k]
                  + f_150 * ld_87[k]
                  + f_117 * ld_126[k]
                  - f_117 * ld_129[k]
                  - f_149 * ld_138[k]
                  + f_149 * ld_141[k]
                  + f_151 * ld_162[k]
                  - f_151 * ld_165[k]
                  + f_148 * ld_216[k]
                  - f_148 * ld_219[k]
                  - f_149 * ld_228[k]
                  + f_149 * ld_231[k]
                  + f_150 * ld_240[k]
                  - f_150 * ld_243[k]
                  - f_151 * ld_252[k]
                  + f_151 * ld_255[k];
    }

#pragma omp simd aligned(ld_13, ld_43, ld_55, ld_97, ld_109, ld_121, ld_175, ld_187, \
                         ld_199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_73 * ld_13[k]
                  - f_73 * ld_43[k]
                  - f_76 * ld_55[k]
                  - f_71 * ld_97[k]
                  + f_74 * ld_109[k]
                  + f_77 * ld_121[k]
                  - f_70 * ld_175[k]
                  + f_72 * ld_187[k]
                  - f_75 * ld_199[k];
    }

#pragma omp simd aligned(ld_16, ld_46, ld_58, ld_100, ld_112, ld_124, ld_178, ld_190, \
                         ld_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_73 * ld_16[k]
                  - f_73 * ld_46[k]
                  - f_76 * ld_58[k]
                  - f_71 * ld_100[k]
                  + f_74 * ld_112[k]
                  + f_77 * ld_124[k]
                  - f_70 * ld_178[k]
                  + f_72 * ld_190[k]
                  - f_75 * ld_202[k];
    }

#pragma omp simd aligned(ld_12, ld_15, ld_17, ld_42, ld_45, ld_47, ld_54, ld_57, ld_59, ld_96, \
                         ld_99, ld_101, ld_108, ld_111, ld_113, ld_120, ld_123, ld_125, \
                         ld_174, ld_177, ld_179, ld_186, ld_189, ld_191, ld_198, ld_201, \
                         ld_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = -f_84 * ld_12[k]
                  - f_84 * ld_15[k]
                  + f_85 * ld_17[k]
                  + f_84 * ld_42[k]
                  + f_84 * ld_45[k]
                  - f_85 * ld_47[k]
                  + f_90 * ld_54[k]
                  + f_90 * ld_57[k]
                  - f_86 * ld_59[k]
                  + f_80 * ld_96[k]
                  + f_80 * ld_99[k]
                  - f_81 * ld_101[k]
                  - f_86 * ld_108[k]
                  - f_86 * ld_111[k]
                  + f_87 * ld_113[k]
                  - f_91 * ld_120[k]
                  - f_91 * ld_123[k]
                  + f_92 * ld_125[k]
                  + f_78 * ld_174[k]
                  + f_78 * ld_177[k]
                  - f_79 * ld_179[k]
                  - f_82 * ld_186[k]
                  - f_82 * ld_189[k]
                  + f_83 * ld_191[k]
                  + f_88 * ld_198[k]
                  + f_88 * ld_201[k]
                  - f_89 * ld_203[k];
    }

#pragma omp simd aligned(ld_14, ld_44, ld_56, ld_98, ld_110, ld_122, ld_176, ld_188, \
                         ld_200 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = f_73 * ld_14[k]
                  - f_73 * ld_44[k]
                  - f_76 * ld_56[k]
                  - f_71 * ld_98[k]
                  + f_74 * ld_110[k]
                  + f_77 * ld_122[k]
                  - f_70 * ld_176[k]
                  + f_72 * ld_188[k]
                  - f_75 * ld_200[k];
    }

#pragma omp simd aligned(ld_12, ld_15, ld_42, ld_45, ld_54, ld_57, ld_96, ld_99, ld_108, \
                         ld_111, ld_120, ld_123, ld_174, ld_177, ld_186, ld_189, ld_198, \
                         ld_201 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_96 * ld_12[k]
                  - f_96 * ld_15[k]
                  - f_96 * ld_42[k]
                  + f_96 * ld_45[k]
                  - f_98 * ld_54[k]
                  + f_98 * ld_57[k]
                  - f_94 * ld_96[k]
                  + f_94 * ld_99[k]
                  + f_76 * ld_108[k]
                  - f_76 * ld_111[k]
                  + f_99 * ld_120[k]
                  - f_99 * ld_123[k]
                  - f_93 * ld_174[k]
                  + f_93 * ld_177[k]
                  + f_95 * ld_186[k]
                  - f_95 * ld_189[k]
                  - f_97 * ld_198[k]
                  + f_97 * ld_201[k];
    }

#pragma omp simd aligned(ld_1, ld_19, ld_31, ld_61, ld_73, ld_85, ld_127, ld_139, ld_151, \
                         ld_217, ld_229, ld_241 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_152 * ld_1[k]
                  - f_58 * ld_19[k]
                  - f_153 * ld_31[k]
                  - f_154 * ld_61[k]
                  + f_155 * ld_73[k]
                  + f_156 * ld_85[k]
                  - f_58 * ld_127[k]
                  + f_155 * ld_139[k]
                  - f_157 * ld_151[k]
                  + f_152 * ld_217[k]
                  - f_153 * ld_229[k]
                  + f_156 * ld_241[k];
    }

#pragma omp simd aligned(ld_4, ld_22, ld_34, ld_64, ld_76, ld_88, ld_130, ld_142, ld_154, \
                         ld_220, ld_232, ld_244 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_152 * ld_4[k]
                  - f_58 * ld_22[k]
                  - f_153 * ld_34[k]
                  - f_154 * ld_64[k]
                  + f_155 * ld_76[k]
                  + f_156 * ld_88[k]
                  - f_58 * ld_130[k]
                  + f_155 * ld_142[k]
                  - f_157 * ld_154[k]
                  + f_152 * ld_220[k]
                  - f_153 * ld_232[k]
                  + f_156 * ld_244[k];
    }

#pragma omp simd aligned(ld_0, ld_3, ld_5, ld_18, ld_21, ld_23, ld_30, ld_33, ld_35, ld_60, \
                         ld_63, ld_65, ld_72, ld_75, ld_77, ld_84, ld_87, ld_89, ld_126, \
                         ld_129, ld_131, ld_138, ld_141, ld_143, ld_150, ld_153, ld_155, \
                         ld_216, ld_219, ld_221, ld_228, ld_231, ld_233, ld_240, ld_243, \
                         ld_245 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_158 * ld_0[k]
                  - f_158 * ld_3[k]
                  + f_159 * ld_5[k]
                  + f_61 * ld_18[k]
                  + f_61 * ld_21[k]
                  - f_62 * ld_23[k]
                  + f_160 * ld_30[k]
                  + f_160 * ld_33[k]
                  - f_161 * ld_35[k]
                  + f_162 * ld_60[k]
                  + f_162 * ld_63[k]
                  - f_163 * ld_65[k]
                  - f_164 * ld_72[k]
                  - f_164 * ld_75[k]
                  + f_165 * ld_77[k]
                  - f_166 * ld_84[k]
                  - f_166 * ld_87[k]
                  + f_167 * ld_89[k]
                  + f_61 * ld_126[k]
                  + f_61 * ld_129[k]
                  - f_62 * ld_131[k]
                  - f_164 * ld_138[k]
                  - f_164 * ld_141[k]
                  + f_165 * ld_143[k]
                  + f_165 * ld_150[k]
                  + f_165 * ld_153[k]
                  - f_168 * ld_155[k]
                  - f_158 * ld_216[k]
                  - f_158 * ld_219[k]
                  + f_159 * ld_221[k]
                  + f_160 * ld_228[k]
                  + f_160 * ld_231[k]
                  - f_161 * ld_233[k]
                  - f_166 * ld_240[k]
                  - f_166 * ld_243[k]
                  + f_167 * ld_245[k];
    }

#pragma omp simd aligned(ld_2, ld_20, ld_32, ld_62, ld_74, ld_86, ld_128, ld_140, ld_152, \
                         ld_218, ld_230, ld_242 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_152 * ld_2[k]
                  - f_58 * ld_20[k]
                  - f_153 * ld_32[k]
                  - f_154 * ld_62[k]
                  + f_155 * ld_74[k]
                  + f_156 * ld_86[k]
                  - f_58 * ld_128[k]
                  + f_155 * ld_140[k]
                  - f_157 * ld_152[k]
                  + f_152 * ld_218[k]
                  - f_153 * ld_230[k]
                  + f_156 * ld_242[k];
    }

#pragma omp simd aligned(ld_0, ld_3, ld_18, ld_21, ld_30, ld_33, ld_60, ld_63, ld_72, ld_75, \
                         ld_84, ld_87, ld_126, ld_129, ld_138, ld_141, ld_150, ld_153, ld_216, \
                         ld_219, ld_228, ld_231, ld_240, ld_243 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_169 * ld_0[k]
                  - f_169 * ld_3[k]
                  - f_67 * ld_18[k]
                  + f_67 * ld_21[k]
                  - f_170 * ld_30[k]
                  + f_170 * ld_33[k]
                  - f_171 * ld_60[k]
                  + f_171 * ld_63[k]
                  + f_172 * ld_72[k]
                  - f_172 * ld_75[k]
                  + f_173 * ld_84[k]
                  - f_173 * ld_87[k]
                  - f_67 * ld_126[k]
                  + f_67 * ld_129[k]
                  + f_172 * ld_138[k]
                  - f_172 * ld_141[k]
                  - f_155 * ld_150[k]
                  + f_155 * ld_153[k]
                  + f_169 * ld_216[k]
                  - f_169 * ld_219[k]
                  - f_170 * ld_228[k]
                  + f_170 * ld_231[k]
                  + f_173 * ld_240[k]
                  - f_173 * ld_243[k];
    }

#pragma omp simd aligned(ld_13, ld_16, ld_43, ld_46, ld_55, ld_58, ld_97, ld_100, ld_109, \
                         ld_112, ld_175, ld_178, ld_187, ld_190 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = -f_40 * ld_13[k]
                  + f_38 * ld_43[k]
                  + f_41 * ld_55[k]
                  + f_36 * ld_97[k]
                  - f_39 * ld_109[k]
                  - f_36 * ld_175[k]
                  + f_37 * ld_187[k];

        g_66[k] = -f_40 * ld_16[k]
                  + f_38 * ld_46[k]
                  + f_41 * ld_58[k]
                  + f_36 * ld_100[k]
                  - f_39 * ld_112[k]
                  - f_36 * ld_178[k]
                  + f_37 * ld_190[k];
    }

#pragma omp simd aligned(ld_12, ld_15, ld_17, ld_42, ld_45, ld_47, ld_54, ld_57, ld_59, ld_96, \
                         ld_99, ld_101, ld_108, ld_111, ld_113, ld_174, ld_177, ld_179, \
                         ld_186, ld_189, ld_191 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = f_49 * ld_12[k]
                  + f_49 * ld_15[k]
                  - f_50 * ld_17[k]
                  - f_46 * ld_42[k]
                  - f_46 * ld_45[k]
                  + f_47 * ld_47[k]
                  - f_51 * ld_54[k]
                  - f_51 * ld_57[k]
                  + f_52 * ld_59[k]
                  - f_42 * ld_96[k]
                  - f_42 * ld_99[k]
                  + f_43 * ld_101[k]
                  + f_45 * ld_108[k]
                  + f_45 * ld_111[k]
                  - f_48 * ld_113[k]
                  + f_42 * ld_174[k]
                  + f_42 * ld_177[k]
                  - f_43 * ld_179[k]
                  - f_44 * ld_186[k]
                  - f_44 * ld_189[k]
                  + f_45 * ld_191[k];
    }

#pragma omp simd aligned(ld_14, ld_44, ld_56, ld_98, ld_110, ld_176, \
                         ld_188 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = -f_40 * ld_14[k]
                  + f_38 * ld_44[k]
                  + f_41 * ld_56[k]
                  + f_36 * ld_98[k]
                  - f_39 * ld_110[k]
                  - f_36 * ld_176[k]
                  + f_37 * ld_188[k];
    }

#pragma omp simd aligned(ld_12, ld_15, ld_42, ld_45, ld_54, ld_57, ld_96, ld_99, ld_108, \
                         ld_111, ld_174, ld_177, ld_186, ld_189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = -f_56 * ld_12[k]
                  + f_56 * ld_15[k]
                  + f_55 * ld_42[k]
                  - f_55 * ld_45[k]
                  + f_57 * ld_54[k]
                  - f_57 * ld_57[k]
                  + f_53 * ld_96[k]
                  - f_53 * ld_99[k]
                  - f_37 * ld_108[k]
                  + f_37 * ld_111[k]
                  - f_53 * ld_174[k]
                  + f_53 * ld_177[k]
                  + f_54 * ld_186[k]
                  - f_54 * ld_189[k];
    }

#pragma omp simd aligned(ld_1, ld_19, ld_31, ld_73, ld_127, ld_139, ld_217, \
                         ld_229 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_174 * ld_1[k]
                  + f_21 * ld_19[k]
                  + f_21 * ld_31[k]
                  - f_175 * ld_73[k]
                  - f_21 * ld_127[k]
                  + f_175 * ld_139[k]
                  + f_174 * ld_217[k]
                  - f_21 * ld_229[k];
    }

#pragma omp simd aligned(ld_4, ld_22, ld_34, ld_76, ld_130, ld_142, ld_220, \
                         ld_232 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_174 * ld_4[k]
                  + f_21 * ld_22[k]
                  + f_21 * ld_34[k]
                  - f_175 * ld_76[k]
                  - f_21 * ld_130[k]
                  + f_175 * ld_142[k]
                  + f_174 * ld_220[k]
                  - f_21 * ld_232[k];
    }

#pragma omp simd aligned(ld_0, ld_3, ld_5, ld_18, ld_21, ld_23, ld_30, ld_33, ld_35, ld_72, \
                         ld_75, ld_77, ld_126, ld_129, ld_131, ld_138, ld_141, ld_143, ld_216, \
                         ld_219, ld_221, ld_228, ld_231, ld_233 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = f_176 * ld_0[k]
                  + f_176 * ld_3[k]
                  - f_177 * ld_5[k]
                  - f_26 * ld_18[k]
                  - f_26 * ld_21[k]
                  + f_27 * ld_23[k]
                  - f_26 * ld_30[k]
                  - f_26 * ld_33[k]
                  + f_27 * ld_35[k]
                  + f_178 * ld_72[k]
                  + f_178 * ld_75[k]
                  - f_179 * ld_77[k]
                  + f_26 * ld_126[k]
                  + f_26 * ld_129[k]
                  - f_27 * ld_131[k]
                  - f_178 * ld_138[k]
                  - f_178 * ld_141[k]
                  + f_179 * ld_143[k]
                  - f_176 * ld_216[k]
                  - f_176 * ld_219[k]
                  + f_177 * ld_221[k]
                  + f_26 * ld_228[k]
                  + f_26 * ld_231[k]
                  - f_27 * ld_233[k];
    }

#pragma omp simd aligned(ld_2, ld_20, ld_32, ld_74, ld_128, ld_140, ld_218, \
                         ld_230 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = -f_174 * ld_2[k]
                  + f_21 * ld_20[k]
                  + f_21 * ld_32[k]
                  - f_175 * ld_74[k]
                  - f_21 * ld_128[k]
                  + f_175 * ld_140[k]
                  + f_174 * ld_218[k]
                  - f_21 * ld_230[k];
    }

#pragma omp simd aligned(ld_0, ld_3, ld_18, ld_21, ld_30, ld_33, ld_72, ld_75, ld_126, ld_129, \
                         ld_138, ld_141, ld_216, ld_219, ld_228, \
                         ld_231 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_180 * ld_0[k]
                  + f_180 * ld_3[k]
                  + f_33 * ld_18[k]
                  - f_33 * ld_21[k]
                  + f_33 * ld_30[k]
                  - f_33 * ld_33[k]
                  - f_181 * ld_72[k]
                  + f_181 * ld_75[k]
                  - f_33 * ld_126[k]
                  + f_33 * ld_129[k]
                  + f_181 * ld_138[k]
                  - f_181 * ld_141[k]
                  + f_180 * ld_216[k]
                  - f_180 * ld_219[k]
                  - f_33 * ld_228[k]
                  + f_33 * ld_231[k];
    }

#pragma omp simd aligned(ld_13, ld_16, ld_43, ld_46, ld_97, ld_100, ld_175, \
                         ld_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = f_6 * ld_13[k]
                  - f_9 * ld_43[k]
                  + f_8 * ld_97[k]
                  - f_7 * ld_175[k];

        g_76[k] = f_6 * ld_16[k]
                  - f_9 * ld_46[k]
                  + f_8 * ld_100[k]
                  - f_7 * ld_178[k];
    }

#pragma omp simd aligned(ld_12, ld_15, ld_17, ld_42, ld_45, ld_47, ld_96, ld_99, ld_101, \
                         ld_174, ld_177, ld_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = -f_15 * ld_12[k]
                  - f_15 * ld_15[k]
                  + f_2 * ld_17[k]
                  + f_13 * ld_42[k]
                  + f_13 * ld_45[k]
                  - f_14 * ld_47[k]
                  - f_11 * ld_96[k]
                  - f_11 * ld_99[k]
                  + f_12 * ld_101[k]
                  + f_10 * ld_174[k]
                  + f_10 * ld_177[k]
                  - f_4 * ld_179[k];
    }

#pragma omp simd aligned(ld_12, ld_14, ld_15, ld_42, ld_44, ld_45, ld_96, ld_98, ld_99, \
                         ld_174, ld_176, ld_177 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = f_6 * ld_14[k]
                  - f_9 * ld_44[k]
                  + f_8 * ld_98[k]
                  - f_7 * ld_176[k];

        g_79[k] = f_19 * ld_12[k]
                  - f_19 * ld_15[k]
                  - f_18 * ld_42[k]
                  + f_18 * ld_45[k]
                  + f_17 * ld_96[k]
                  - f_17 * ld_99[k]
                  - f_16 * ld_174[k]
                  + f_16 * ld_177[k];
    }

#pragma omp simd aligned(ld_1, ld_4, ld_19, ld_22, ld_61, ld_64, ld_127, ld_130, ld_217, \
                         ld_220 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = f_182 * ld_1[k]
                  - f_7 * ld_19[k]
                  + f_17 * ld_61[k]
                  - f_7 * ld_127[k]
                  + f_182 * ld_217[k];

        g_81[k] = f_182 * ld_4[k]
                  - f_7 * ld_22[k]
                  + f_17 * ld_64[k]
                  - f_7 * ld_130[k]
                  + f_182 * ld_220[k];
    }

#pragma omp simd aligned(ld_0, ld_3, ld_5, ld_18, ld_21, ld_23, ld_60, ld_63, ld_65, ld_126, \
                         ld_129, ld_131, ld_216, ld_219, ld_221 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = -f_183 * ld_0[k]
                  - f_183 * ld_3[k]
                  + f_184 * ld_5[k]
                  + f_10 * ld_18[k]
                  + f_10 * ld_21[k]
                  - f_4 * ld_23[k]
                  - f_185 * ld_60[k]
                  - f_185 * ld_63[k]
                  + f_11 * ld_65[k]
                  + f_10 * ld_126[k]
                  + f_10 * ld_129[k]
                  - f_4 * ld_131[k]
                  - f_183 * ld_216[k]
                  - f_183 * ld_219[k]
                  + f_184 * ld_221[k];
    }

#pragma omp simd aligned(ld_2, ld_20, ld_62, ld_128, ld_218 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = f_182 * ld_2[k]
                  - f_7 * ld_20[k]
                  + f_17 * ld_62[k]
                  - f_7 * ld_128[k]
                  + f_182 * ld_218[k];
    }

#pragma omp simd aligned(ld_0, ld_3, ld_18, ld_21, ld_60, ld_63, ld_126, ld_129, ld_216, \
                         ld_219 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = f_186 * ld_0[k]
                  - f_186 * ld_3[k]
                  - f_16 * ld_18[k]
                  + f_16 * ld_21[k]
                  + f_187 * ld_60[k]
                  - f_187 * ld_63[k]
                  - f_16 * ld_126[k]
                  + f_16 * ld_129[k]
                  + f_186 * ld_216[k]
                  - f_186 * ld_219[k];
    }
}

}  // namespace simdtrf
