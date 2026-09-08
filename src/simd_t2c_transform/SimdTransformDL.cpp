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


#include "SimdTransformDL.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_dl(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t dl,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.1875 * std::sqrt(2145.0);
    const auto f_1 = 1.3125 * std::sqrt(2145.0);
    const auto f_2 = 0.65625 * std::sqrt(2145.0);
    const auto f_3 = 3.28125 * std::sqrt(2145.0);
    const auto f_4 = 1.96875 * std::sqrt(2145.0);
    const auto f_5 = 0.09375 * std::sqrt(2145.0);
    const auto f_6 = 0.28125 * std::sqrt(286.0);
    const auto f_7 = 0.65625 * std::sqrt(286.0);
    const auto f_8 = 3.9375 * std::sqrt(286.0);
    const auto f_9 = 13.125 * std::sqrt(286.0);
    const auto f_10 = 0.46875 * std::sqrt(3003.0);
    const auto f_11 = 1.875 * std::sqrt(3003.0);
    const auto f_12 = 0.84375 * std::sqrt(3003.0);
    const auto f_13 = 3.75 * std::sqrt(3003.0);
    const auto f_14 = 0.09375 * std::sqrt(3003.0);
    const auto f_15 = 0.375 * std::sqrt(3003.0);
    const auto f_16 = 0.1875 * std::sqrt(231.0);
    const auto f_17 = 4.5 * std::sqrt(231.0);
    const auto f_18 = 7.5 * std::sqrt(231.0);
    const auto f_19 = 0.84375 * std::sqrt(385.0);
    const auto f_20 = 1.40625 * std::sqrt(385.0);
    const auto f_21 = 5.625 * std::sqrt(385.0);
    const auto f_22 = 0.28125 * std::sqrt(385.0);
    const auto f_23 = 3.75 * std::sqrt(385.0);
    const auto f_24 = 4.5 * std::sqrt(385.0);
    const auto f_25 = 1.875 * std::sqrt(385.0);
    const auto f_26 = 1.5 * std::sqrt(385.0);
    const auto f_27 = 0.09375 * std::sqrt(210.0);
    const auto f_28 = 0.28125 * std::sqrt(210.0);
    const auto f_29 = 2.8125 * std::sqrt(210.0);
    const auto f_30 = 5.625 * std::sqrt(210.0);
    const auto f_31 = 7.5 * std::sqrt(210.0);
    const auto f_32 = 3.0 * std::sqrt(210.0);
    const auto f_33 = 3.28125 * std::sqrt(3.0);
    const auto f_34 = 9.84375 * std::sqrt(3.0);
    const auto f_35 = 26.25 * std::sqrt(3.0);
    const auto f_36 = 52.5 * std::sqrt(3.0);
    const auto f_37 = 31.5 * std::sqrt(3.0);
    const auto f_38 = 6.0 * std::sqrt(3.0);
    const auto f_39 = 0.2734375 * std::sqrt(3.0);
    const auto f_40 = 1.09375 * std::sqrt(3.0);
    const auto f_41 = 8.75 * std::sqrt(3.0);
    const auto f_42 = 1.640625 * std::sqrt(3.0);
    const auto f_43 = 14.0 * std::sqrt(3.0);
    const auto f_44 = std::sqrt(3.0);
    const auto f_45 = 0.046875 * std::sqrt(210.0);
    const auto f_46 = 1.40625 * std::sqrt(210.0);
    const auto f_47 = 3.75 * std::sqrt(210.0);
    const auto f_48 = 1.5 * std::sqrt(210.0);
    const auto f_49 = 0.046875 * std::sqrt(231.0);
    const auto f_50 = 1.125 * std::sqrt(231.0);
    const auto f_51 = 0.46875 * std::sqrt(231.0);
    const auto f_52 = 5.625 * std::sqrt(231.0);
    const auto f_53 = 1.875 * std::sqrt(231.0);
    const auto f_54 = 11.25 * std::sqrt(231.0);
    const auto f_55 = 0.046875 * std::sqrt(286.0);
    const auto f_56 = 9.84375 * std::sqrt(286.0);
    const auto f_57 = 0.0234375 * std::sqrt(2145.0);
    const auto f_58 = 1.640625 * std::sqrt(2145.0);
    const auto f_59 = 0.09375 * std::sqrt(715.0);
    const auto f_60 = 0.65625 * std::sqrt(715.0);
    const auto f_61 = 0.1875 * std::sqrt(715.0);
    const auto f_62 = 1.3125 * std::sqrt(715.0);
    const auto f_63 = 0.328125 * std::sqrt(715.0);
    const auto f_64 = 1.640625 * std::sqrt(715.0);
    const auto f_65 = 0.984375 * std::sqrt(715.0);
    const auto f_66 = 0.046875 * std::sqrt(715.0);
    const auto f_67 = 3.28125 * std::sqrt(715.0);
    const auto f_68 = 1.96875 * std::sqrt(715.0);
    const auto f_69 = 0.046875 * std::sqrt(858.0);
    const auto f_70 = 0.109375 * std::sqrt(858.0);
    const auto f_71 = 0.65625 * std::sqrt(858.0);
    const auto f_72 = 2.1875 * std::sqrt(858.0);
    const auto f_73 = 0.09375 * std::sqrt(858.0);
    const auto f_74 = 0.21875 * std::sqrt(858.0);
    const auto f_75 = 1.3125 * std::sqrt(858.0);
    const auto f_76 = 4.375 * std::sqrt(858.0);
    const auto f_77 = 0.234375 * std::sqrt(1001.0);
    const auto f_78 = 0.9375 * std::sqrt(1001.0);
    const auto f_79 = 0.421875 * std::sqrt(1001.0);
    const auto f_80 = 1.875 * std::sqrt(1001.0);
    const auto f_81 = 0.046875 * std::sqrt(1001.0);
    const auto f_82 = 0.1875 * std::sqrt(1001.0);
    const auto f_83 = 0.46875 * std::sqrt(1001.0);
    const auto f_84 = 0.84375 * std::sqrt(1001.0);
    const auto f_85 = 3.75 * std::sqrt(1001.0);
    const auto f_86 = 0.09375 * std::sqrt(1001.0);
    const auto f_87 = 0.375 * std::sqrt(1001.0);
    const auto f_88 = 0.09375 * std::sqrt(77.0);
    const auto f_89 = 2.25 * std::sqrt(77.0);
    const auto f_90 = 3.75 * std::sqrt(77.0);
    const auto f_91 = 0.1875 * std::sqrt(77.0);
    const auto f_92 = 4.5 * std::sqrt(77.0);
    const auto f_93 = 7.5 * std::sqrt(77.0);
    const auto f_94 = 0.140625 * std::sqrt(1155.0);
    const auto f_95 = 0.234375 * std::sqrt(1155.0);
    const auto f_96 = 0.9375 * std::sqrt(1155.0);
    const auto f_97 = 0.046875 * std::sqrt(1155.0);
    const auto f_98 = 0.625 * std::sqrt(1155.0);
    const auto f_99 = 0.75 * std::sqrt(1155.0);
    const auto f_100 = 0.3125 * std::sqrt(1155.0);
    const auto f_101 = 0.25 * std::sqrt(1155.0);
    const auto f_102 = 0.28125 * std::sqrt(1155.0);
    const auto f_103 = 0.46875 * std::sqrt(1155.0);
    const auto f_104 = 1.875 * std::sqrt(1155.0);
    const auto f_105 = 0.09375 * std::sqrt(1155.0);
    const auto f_106 = 1.25 * std::sqrt(1155.0);
    const auto f_107 = 1.5 * std::sqrt(1155.0);
    const auto f_108 = 0.5 * std::sqrt(1155.0);
    const auto f_109 = 0.046875 * std::sqrt(70.0);
    const auto f_110 = 0.140625 * std::sqrt(70.0);
    const auto f_111 = 1.40625 * std::sqrt(70.0);
    const auto f_112 = 2.8125 * std::sqrt(70.0);
    const auto f_113 = 3.75 * std::sqrt(70.0);
    const auto f_114 = 1.5 * std::sqrt(70.0);
    const auto f_115 = 0.09375 * std::sqrt(70.0);
    const auto f_116 = 0.28125 * std::sqrt(70.0);
    const auto f_117 = 5.625 * std::sqrt(70.0);
    const auto f_118 = 7.5 * std::sqrt(70.0);
    const auto f_119 = 3.0 * std::sqrt(70.0);
    const auto f_120 = 0.0234375 * std::sqrt(70.0);
    const auto f_121 = 0.703125 * std::sqrt(70.0);
    const auto f_122 = 1.875 * std::sqrt(70.0);
    const auto f_123 = 0.75 * std::sqrt(70.0);
    const auto f_124 = 0.0234375 * std::sqrt(77.0);
    const auto f_125 = 0.5625 * std::sqrt(77.0);
    const auto f_126 = 0.234375 * std::sqrt(77.0);
    const auto f_127 = 2.8125 * std::sqrt(77.0);
    const auto f_128 = 0.9375 * std::sqrt(77.0);
    const auto f_129 = 5.625 * std::sqrt(77.0);
    const auto f_130 = 0.046875 * std::sqrt(77.0);
    const auto f_131 = 1.125 * std::sqrt(77.0);
    const auto f_132 = 0.46875 * std::sqrt(77.0);
    const auto f_133 = 1.875 * std::sqrt(77.0);
    const auto f_134 = 11.25 * std::sqrt(77.0);
    const auto f_135 = 0.0078125 * std::sqrt(858.0);
    const auto f_136 = 1.640625 * std::sqrt(858.0);
    const auto f_137 = 0.015625 * std::sqrt(858.0);
    const auto f_138 = 3.28125 * std::sqrt(858.0);
    const auto f_139 = 0.01171875 * std::sqrt(715.0);
    const auto f_140 = 0.8203125 * std::sqrt(715.0);
    const auto f_141 = 0.0234375 * std::sqrt(715.0);
    const auto f_142 = 0.328125 * std::sqrt(2145.0);
    const auto f_143 = 0.984375 * std::sqrt(2145.0);
    const auto f_144 = 0.046875 * std::sqrt(2145.0);
    const auto f_145 = 0.140625 * std::sqrt(286.0);
    const auto f_146 = 0.328125 * std::sqrt(286.0);
    const auto f_147 = 1.96875 * std::sqrt(286.0);
    const auto f_148 = 6.5625 * std::sqrt(286.0);
    const auto f_149 = 0.234375 * std::sqrt(3003.0);
    const auto f_150 = 0.9375 * std::sqrt(3003.0);
    const auto f_151 = 0.421875 * std::sqrt(3003.0);
    const auto f_152 = 0.046875 * std::sqrt(3003.0);
    const auto f_153 = 0.1875 * std::sqrt(3003.0);
    const auto f_154 = 0.09375 * std::sqrt(231.0);
    const auto f_155 = 2.25 * std::sqrt(231.0);
    const auto f_156 = 3.75 * std::sqrt(231.0);
    const auto f_157 = 0.421875 * std::sqrt(385.0);
    const auto f_158 = 0.703125 * std::sqrt(385.0);
    const auto f_159 = 2.8125 * std::sqrt(385.0);
    const auto f_160 = 0.140625 * std::sqrt(385.0);
    const auto f_161 = 2.25 * std::sqrt(385.0);
    const auto f_162 = 0.9375 * std::sqrt(385.0);
    const auto f_163 = 0.75 * std::sqrt(385.0);
    const auto f_164 = 0.140625 * std::sqrt(210.0);
    const auto f_165 = 4.921875 * std::sqrt(3.0);
    const auto f_166 = 13.125 * std::sqrt(3.0);
    const auto f_167 = 15.75 * std::sqrt(3.0);
    const auto f_168 = 3.0 * std::sqrt(3.0);
    const auto f_169 = 0.13671875 * std::sqrt(3.0);
    const auto f_170 = 0.546875 * std::sqrt(3.0);
    const auto f_171 = 4.375 * std::sqrt(3.0);
    const auto f_172 = 0.8203125 * std::sqrt(3.0);
    const auto f_173 = 7.0 * std::sqrt(3.0);
    const auto f_174 = 0.5 * std::sqrt(3.0);
    const auto f_175 = 0.0234375 * std::sqrt(210.0);
    const auto f_176 = 0.703125 * std::sqrt(210.0);
    const auto f_177 = 1.875 * std::sqrt(210.0);
    const auto f_178 = 0.75 * std::sqrt(210.0);
    const auto f_179 = 0.0234375 * std::sqrt(231.0);
    const auto f_180 = 0.5625 * std::sqrt(231.0);
    const auto f_181 = 0.234375 * std::sqrt(231.0);
    const auto f_182 = 2.8125 * std::sqrt(231.0);
    const auto f_183 = 0.9375 * std::sqrt(231.0);
    const auto f_184 = 0.0234375 * std::sqrt(286.0);
    const auto f_185 = 4.921875 * std::sqrt(286.0);
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

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);
    const auto *dl_3 = buffer.data(dl + 3);
    const auto *dl_4 = buffer.data(dl + 4);
    const auto *dl_5 = buffer.data(dl + 5);
    const auto *dl_6 = buffer.data(dl + 6);
    const auto *dl_7 = buffer.data(dl + 7);
    const auto *dl_8 = buffer.data(dl + 8);
    const auto *dl_9 = buffer.data(dl + 9);
    const auto *dl_10 = buffer.data(dl + 10);
    const auto *dl_11 = buffer.data(dl + 11);
    const auto *dl_12 = buffer.data(dl + 12);
    const auto *dl_13 = buffer.data(dl + 13);
    const auto *dl_14 = buffer.data(dl + 14);
    const auto *dl_15 = buffer.data(dl + 15);
    const auto *dl_16 = buffer.data(dl + 16);
    const auto *dl_17 = buffer.data(dl + 17);
    const auto *dl_18 = buffer.data(dl + 18);
    const auto *dl_19 = buffer.data(dl + 19);
    const auto *dl_20 = buffer.data(dl + 20);
    const auto *dl_21 = buffer.data(dl + 21);
    const auto *dl_22 = buffer.data(dl + 22);
    const auto *dl_23 = buffer.data(dl + 23);
    const auto *dl_24 = buffer.data(dl + 24);
    const auto *dl_25 = buffer.data(dl + 25);
    const auto *dl_26 = buffer.data(dl + 26);
    const auto *dl_27 = buffer.data(dl + 27);
    const auto *dl_28 = buffer.data(dl + 28);
    const auto *dl_29 = buffer.data(dl + 29);
    const auto *dl_30 = buffer.data(dl + 30);
    const auto *dl_31 = buffer.data(dl + 31);
    const auto *dl_32 = buffer.data(dl + 32);
    const auto *dl_33 = buffer.data(dl + 33);
    const auto *dl_34 = buffer.data(dl + 34);
    const auto *dl_35 = buffer.data(dl + 35);
    const auto *dl_36 = buffer.data(dl + 36);
    const auto *dl_37 = buffer.data(dl + 37);
    const auto *dl_38 = buffer.data(dl + 38);
    const auto *dl_39 = buffer.data(dl + 39);
    const auto *dl_40 = buffer.data(dl + 40);
    const auto *dl_41 = buffer.data(dl + 41);
    const auto *dl_42 = buffer.data(dl + 42);
    const auto *dl_43 = buffer.data(dl + 43);
    const auto *dl_44 = buffer.data(dl + 44);
    const auto *dl_45 = buffer.data(dl + 45);
    const auto *dl_46 = buffer.data(dl + 46);
    const auto *dl_47 = buffer.data(dl + 47);
    const auto *dl_48 = buffer.data(dl + 48);
    const auto *dl_49 = buffer.data(dl + 49);
    const auto *dl_50 = buffer.data(dl + 50);
    const auto *dl_51 = buffer.data(dl + 51);
    const auto *dl_52 = buffer.data(dl + 52);
    const auto *dl_53 = buffer.data(dl + 53);
    const auto *dl_54 = buffer.data(dl + 54);
    const auto *dl_55 = buffer.data(dl + 55);
    const auto *dl_56 = buffer.data(dl + 56);
    const auto *dl_57 = buffer.data(dl + 57);
    const auto *dl_58 = buffer.data(dl + 58);
    const auto *dl_59 = buffer.data(dl + 59);
    const auto *dl_60 = buffer.data(dl + 60);
    const auto *dl_61 = buffer.data(dl + 61);
    const auto *dl_62 = buffer.data(dl + 62);
    const auto *dl_63 = buffer.data(dl + 63);
    const auto *dl_64 = buffer.data(dl + 64);
    const auto *dl_65 = buffer.data(dl + 65);
    const auto *dl_66 = buffer.data(dl + 66);
    const auto *dl_67 = buffer.data(dl + 67);
    const auto *dl_68 = buffer.data(dl + 68);
    const auto *dl_69 = buffer.data(dl + 69);
    const auto *dl_70 = buffer.data(dl + 70);
    const auto *dl_71 = buffer.data(dl + 71);
    const auto *dl_72 = buffer.data(dl + 72);
    const auto *dl_73 = buffer.data(dl + 73);
    const auto *dl_74 = buffer.data(dl + 74);
    const auto *dl_75 = buffer.data(dl + 75);
    const auto *dl_76 = buffer.data(dl + 76);
    const auto *dl_77 = buffer.data(dl + 77);
    const auto *dl_78 = buffer.data(dl + 78);
    const auto *dl_79 = buffer.data(dl + 79);
    const auto *dl_80 = buffer.data(dl + 80);
    const auto *dl_81 = buffer.data(dl + 81);
    const auto *dl_82 = buffer.data(dl + 82);
    const auto *dl_83 = buffer.data(dl + 83);
    const auto *dl_84 = buffer.data(dl + 84);
    const auto *dl_85 = buffer.data(dl + 85);
    const auto *dl_86 = buffer.data(dl + 86);
    const auto *dl_87 = buffer.data(dl + 87);
    const auto *dl_88 = buffer.data(dl + 88);
    const auto *dl_89 = buffer.data(dl + 89);
    const auto *dl_90 = buffer.data(dl + 90);
    const auto *dl_91 = buffer.data(dl + 91);
    const auto *dl_92 = buffer.data(dl + 92);
    const auto *dl_93 = buffer.data(dl + 93);
    const auto *dl_94 = buffer.data(dl + 94);
    const auto *dl_95 = buffer.data(dl + 95);
    const auto *dl_96 = buffer.data(dl + 96);
    const auto *dl_97 = buffer.data(dl + 97);
    const auto *dl_98 = buffer.data(dl + 98);
    const auto *dl_99 = buffer.data(dl + 99);
    const auto *dl_100 = buffer.data(dl + 100);
    const auto *dl_101 = buffer.data(dl + 101);
    const auto *dl_102 = buffer.data(dl + 102);
    const auto *dl_103 = buffer.data(dl + 103);
    const auto *dl_104 = buffer.data(dl + 104);
    const auto *dl_105 = buffer.data(dl + 105);
    const auto *dl_106 = buffer.data(dl + 106);
    const auto *dl_107 = buffer.data(dl + 107);
    const auto *dl_108 = buffer.data(dl + 108);
    const auto *dl_109 = buffer.data(dl + 109);
    const auto *dl_110 = buffer.data(dl + 110);
    const auto *dl_111 = buffer.data(dl + 111);
    const auto *dl_112 = buffer.data(dl + 112);
    const auto *dl_113 = buffer.data(dl + 113);
    const auto *dl_114 = buffer.data(dl + 114);
    const auto *dl_115 = buffer.data(dl + 115);
    const auto *dl_116 = buffer.data(dl + 116);
    const auto *dl_117 = buffer.data(dl + 117);
    const auto *dl_118 = buffer.data(dl + 118);
    const auto *dl_119 = buffer.data(dl + 119);
    const auto *dl_120 = buffer.data(dl + 120);
    const auto *dl_121 = buffer.data(dl + 121);
    const auto *dl_122 = buffer.data(dl + 122);
    const auto *dl_123 = buffer.data(dl + 123);
    const auto *dl_124 = buffer.data(dl + 124);
    const auto *dl_125 = buffer.data(dl + 125);
    const auto *dl_126 = buffer.data(dl + 126);
    const auto *dl_127 = buffer.data(dl + 127);
    const auto *dl_128 = buffer.data(dl + 128);
    const auto *dl_129 = buffer.data(dl + 129);
    const auto *dl_130 = buffer.data(dl + 130);
    const auto *dl_131 = buffer.data(dl + 131);
    const auto *dl_132 = buffer.data(dl + 132);
    const auto *dl_133 = buffer.data(dl + 133);
    const auto *dl_134 = buffer.data(dl + 134);
    const auto *dl_135 = buffer.data(dl + 135);
    const auto *dl_136 = buffer.data(dl + 136);
    const auto *dl_137 = buffer.data(dl + 137);
    const auto *dl_138 = buffer.data(dl + 138);
    const auto *dl_139 = buffer.data(dl + 139);
    const auto *dl_140 = buffer.data(dl + 140);
    const auto *dl_141 = buffer.data(dl + 141);
    const auto *dl_142 = buffer.data(dl + 142);
    const auto *dl_143 = buffer.data(dl + 143);
    const auto *dl_144 = buffer.data(dl + 144);
    const auto *dl_145 = buffer.data(dl + 145);
    const auto *dl_146 = buffer.data(dl + 146);
    const auto *dl_147 = buffer.data(dl + 147);
    const auto *dl_148 = buffer.data(dl + 148);
    const auto *dl_149 = buffer.data(dl + 149);
    const auto *dl_150 = buffer.data(dl + 150);
    const auto *dl_151 = buffer.data(dl + 151);
    const auto *dl_152 = buffer.data(dl + 152);
    const auto *dl_153 = buffer.data(dl + 153);
    const auto *dl_154 = buffer.data(dl + 154);
    const auto *dl_155 = buffer.data(dl + 155);
    const auto *dl_156 = buffer.data(dl + 156);
    const auto *dl_157 = buffer.data(dl + 157);
    const auto *dl_158 = buffer.data(dl + 158);
    const auto *dl_159 = buffer.data(dl + 159);
    const auto *dl_160 = buffer.data(dl + 160);
    const auto *dl_161 = buffer.data(dl + 161);
    const auto *dl_162 = buffer.data(dl + 162);
    const auto *dl_163 = buffer.data(dl + 163);
    const auto *dl_164 = buffer.data(dl + 164);
    const auto *dl_165 = buffer.data(dl + 165);
    const auto *dl_166 = buffer.data(dl + 166);
    const auto *dl_167 = buffer.data(dl + 167);
    const auto *dl_168 = buffer.data(dl + 168);
    const auto *dl_169 = buffer.data(dl + 169);
    const auto *dl_170 = buffer.data(dl + 170);
    const auto *dl_171 = buffer.data(dl + 171);
    const auto *dl_172 = buffer.data(dl + 172);
    const auto *dl_173 = buffer.data(dl + 173);
    const auto *dl_174 = buffer.data(dl + 174);
    const auto *dl_175 = buffer.data(dl + 175);
    const auto *dl_176 = buffer.data(dl + 176);
    const auto *dl_177 = buffer.data(dl + 177);
    const auto *dl_178 = buffer.data(dl + 178);
    const auto *dl_179 = buffer.data(dl + 179);
    const auto *dl_180 = buffer.data(dl + 180);
    const auto *dl_181 = buffer.data(dl + 181);
    const auto *dl_182 = buffer.data(dl + 182);
    const auto *dl_183 = buffer.data(dl + 183);
    const auto *dl_184 = buffer.data(dl + 184);
    const auto *dl_185 = buffer.data(dl + 185);
    const auto *dl_186 = buffer.data(dl + 186);
    const auto *dl_187 = buffer.data(dl + 187);
    const auto *dl_188 = buffer.data(dl + 188);
    const auto *dl_189 = buffer.data(dl + 189);
    const auto *dl_190 = buffer.data(dl + 190);
    const auto *dl_191 = buffer.data(dl + 191);
    const auto *dl_192 = buffer.data(dl + 192);
    const auto *dl_193 = buffer.data(dl + 193);
    const auto *dl_194 = buffer.data(dl + 194);
    const auto *dl_195 = buffer.data(dl + 195);
    const auto *dl_196 = buffer.data(dl + 196);
    const auto *dl_197 = buffer.data(dl + 197);
    const auto *dl_198 = buffer.data(dl + 198);
    const auto *dl_199 = buffer.data(dl + 199);
    const auto *dl_200 = buffer.data(dl + 200);
    const auto *dl_201 = buffer.data(dl + 201);
    const auto *dl_202 = buffer.data(dl + 202);
    const auto *dl_203 = buffer.data(dl + 203);
    const auto *dl_204 = buffer.data(dl + 204);
    const auto *dl_205 = buffer.data(dl + 205);
    const auto *dl_206 = buffer.data(dl + 206);
    const auto *dl_207 = buffer.data(dl + 207);
    const auto *dl_208 = buffer.data(dl + 208);
    const auto *dl_209 = buffer.data(dl + 209);
    const auto *dl_210 = buffer.data(dl + 210);
    const auto *dl_211 = buffer.data(dl + 211);
    const auto *dl_212 = buffer.data(dl + 212);
    const auto *dl_213 = buffer.data(dl + 213);
    const auto *dl_214 = buffer.data(dl + 214);
    const auto *dl_215 = buffer.data(dl + 215);
    const auto *dl_216 = buffer.data(dl + 216);
    const auto *dl_217 = buffer.data(dl + 217);
    const auto *dl_218 = buffer.data(dl + 218);
    const auto *dl_219 = buffer.data(dl + 219);
    const auto *dl_220 = buffer.data(dl + 220);
    const auto *dl_221 = buffer.data(dl + 221);
    const auto *dl_222 = buffer.data(dl + 222);
    const auto *dl_223 = buffer.data(dl + 223);
    const auto *dl_224 = buffer.data(dl + 224);
    const auto *dl_225 = buffer.data(dl + 225);
    const auto *dl_226 = buffer.data(dl + 226);
    const auto *dl_227 = buffer.data(dl + 227);
    const auto *dl_228 = buffer.data(dl + 228);
    const auto *dl_229 = buffer.data(dl + 229);
    const auto *dl_230 = buffer.data(dl + 230);
    const auto *dl_231 = buffer.data(dl + 231);
    const auto *dl_232 = buffer.data(dl + 232);
    const auto *dl_233 = buffer.data(dl + 233);
    const auto *dl_234 = buffer.data(dl + 234);
    const auto *dl_235 = buffer.data(dl + 235);
    const auto *dl_236 = buffer.data(dl + 236);
    const auto *dl_237 = buffer.data(dl + 237);
    const auto *dl_238 = buffer.data(dl + 238);
    const auto *dl_239 = buffer.data(dl + 239);
    const auto *dl_240 = buffer.data(dl + 240);
    const auto *dl_241 = buffer.data(dl + 241);
    const auto *dl_242 = buffer.data(dl + 242);
    const auto *dl_243 = buffer.data(dl + 243);
    const auto *dl_244 = buffer.data(dl + 244);
    const auto *dl_245 = buffer.data(dl + 245);
    const auto *dl_246 = buffer.data(dl + 246);
    const auto *dl_247 = buffer.data(dl + 247);
    const auto *dl_248 = buffer.data(dl + 248);
    const auto *dl_249 = buffer.data(dl + 249);
    const auto *dl_250 = buffer.data(dl + 250);
    const auto *dl_251 = buffer.data(dl + 251);
    const auto *dl_252 = buffer.data(dl + 252);
    const auto *dl_253 = buffer.data(dl + 253);
    const auto *dl_254 = buffer.data(dl + 254);
    const auto *dl_255 = buffer.data(dl + 255);
    const auto *dl_256 = buffer.data(dl + 256);
    const auto *dl_257 = buffer.data(dl + 257);
    const auto *dl_258 = buffer.data(dl + 258);
    const auto *dl_259 = buffer.data(dl + 259);
    const auto *dl_260 = buffer.data(dl + 260);
    const auto *dl_261 = buffer.data(dl + 261);
    const auto *dl_262 = buffer.data(dl + 262);
    const auto *dl_263 = buffer.data(dl + 263);
    const auto *dl_264 = buffer.data(dl + 264);
    const auto *dl_265 = buffer.data(dl + 265);
    const auto *dl_266 = buffer.data(dl + 266);
    const auto *dl_267 = buffer.data(dl + 267);
    const auto *dl_268 = buffer.data(dl + 268);
    const auto *dl_269 = buffer.data(dl + 269);

#pragma omp simd aligned(dl_46, dl_49, dl_51, dl_53, dl_56, dl_60, dl_62, dl_67, dl_73, dl_75, \
                         dl_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * dl_46[k]
                 - f_1 * dl_51[k]
                 + f_1 * dl_60[k]
                 - f_0 * dl_73[k];

        g_1[k] = f_2 * dl_49[k]
                 - f_3 * dl_56[k]
                 + f_4 * dl_67[k]
                 - f_5 * dl_82[k];

        g_2[k] = -f_6 * dl_46[k]
                 + f_7 * dl_51[k]
                 + f_8 * dl_53[k]
                 + f_7 * dl_60[k]
                 - f_9 * dl_62[k]
                 - f_6 * dl_73[k]
                 + f_8 * dl_75[k];
    }

#pragma omp simd aligned(dl_49, dl_56, dl_58, dl_67, dl_69, dl_82, \
                         dl_84 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_10 * dl_49[k]
                 + f_10 * dl_56[k]
                 + f_11 * dl_58[k]
                 + f_12 * dl_67[k]
                 - f_13 * dl_69[k]
                 - f_14 * dl_82[k]
                 + f_15 * dl_84[k];
    }

#pragma omp simd aligned(dl_46, dl_51, dl_53, dl_60, dl_64, dl_73, dl_75, \
                         dl_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_16 * dl_46[k]
                 + f_16 * dl_51[k]
                 - f_17 * dl_53[k]
                 - f_16 * dl_60[k]
                 + f_18 * dl_64[k]
                 - f_16 * dl_73[k]
                 + f_17 * dl_75[k]
                 - f_18 * dl_77[k];
    }

#pragma omp simd aligned(dl_49, dl_56, dl_58, dl_67, dl_69, dl_71, dl_82, dl_84, \
                         dl_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_19 * dl_49[k]
                 + f_20 * dl_56[k]
                 - f_21 * dl_58[k]
                 + f_22 * dl_67[k]
                 - f_23 * dl_69[k]
                 + f_24 * dl_71[k]
                 - f_22 * dl_82[k]
                 + f_25 * dl_84[k]
                 - f_26 * dl_86[k];
    }

#pragma omp simd aligned(dl_46, dl_51, dl_53, dl_60, dl_62, dl_64, dl_73, dl_75, dl_77, \
                         dl_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_27 * dl_46[k]
                 - f_28 * dl_51[k]
                 + f_29 * dl_53[k]
                 - f_28 * dl_60[k]
                 + f_30 * dl_62[k]
                 - f_31 * dl_64[k]
                 - f_27 * dl_73[k]
                 + f_29 * dl_75[k]
                 - f_31 * dl_77[k]
                 + f_32 * dl_79[k];
    }

#pragma omp simd aligned(dl_49, dl_56, dl_58, dl_67, dl_69, dl_71, dl_82, dl_84, dl_86, \
                         dl_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_33 * dl_49[k]
                 - f_34 * dl_56[k]
                 + f_35 * dl_58[k]
                 - f_34 * dl_67[k]
                 + f_36 * dl_69[k]
                 - f_37 * dl_71[k]
                 - f_33 * dl_82[k]
                 + f_35 * dl_84[k]
                 - f_37 * dl_86[k]
                 + f_38 * dl_88[k];
    }

#pragma omp simd aligned(dl_45, dl_48, dl_50, dl_55, dl_57, dl_59, dl_66, dl_68, dl_70, dl_72, \
                         dl_81, dl_83, dl_85, dl_87, dl_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_39 * dl_45[k]
                 + f_40 * dl_48[k]
                 - f_41 * dl_50[k]
                 + f_42 * dl_55[k]
                 - f_35 * dl_57[k]
                 + f_35 * dl_59[k]
                 + f_40 * dl_66[k]
                 - f_35 * dl_68[k]
                 + f_36 * dl_70[k]
                 - f_43 * dl_72[k]
                 + f_39 * dl_81[k]
                 - f_41 * dl_83[k]
                 + f_35 * dl_85[k]
                 - f_43 * dl_87[k]
                 + f_44 * dl_89[k];
    }

#pragma omp simd aligned(dl_47, dl_52, dl_54, dl_61, dl_63, dl_65, dl_74, dl_76, dl_78, \
                         dl_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_33 * dl_47[k]
                 - f_34 * dl_52[k]
                 + f_35 * dl_54[k]
                 - f_34 * dl_61[k]
                 + f_36 * dl_63[k]
                 - f_37 * dl_65[k]
                 - f_33 * dl_74[k]
                 + f_35 * dl_76[k]
                 - f_37 * dl_78[k]
                 + f_38 * dl_80[k];
    }

#pragma omp simd aligned(dl_45, dl_48, dl_50, dl_57, dl_59, dl_66, dl_68, dl_72, dl_81, dl_83, \
                         dl_85, dl_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_45 * dl_45[k]
                  - f_27 * dl_48[k]
                  + f_46 * dl_50[k]
                  + f_46 * dl_57[k]
                  - f_47 * dl_59[k]
                  + f_27 * dl_66[k]
                  - f_46 * dl_68[k]
                  + f_48 * dl_72[k]
                  + f_45 * dl_81[k]
                  - f_46 * dl_83[k]
                  + f_47 * dl_85[k]
                  - f_48 * dl_87[k];
    }

#pragma omp simd aligned(dl_47, dl_52, dl_54, dl_61, dl_63, dl_65, dl_74, dl_76, \
                         dl_78 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_22 * dl_47[k]
                  - f_22 * dl_52[k]
                  - f_25 * dl_54[k]
                  - f_20 * dl_61[k]
                  + f_23 * dl_63[k]
                  + f_26 * dl_65[k]
                  - f_19 * dl_74[k]
                  + f_21 * dl_76[k]
                  - f_24 * dl_78[k];
    }

#pragma omp simd aligned(dl_45, dl_48, dl_50, dl_55, dl_57, dl_59, dl_66, dl_68, dl_70, dl_81, \
                         dl_83, dl_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_49 * dl_45[k]
                  - f_16 * dl_48[k]
                  - f_50 * dl_50[k]
                  - f_51 * dl_55[k]
                  + f_52 * dl_57[k]
                  + f_53 * dl_59[k]
                  - f_16 * dl_66[k]
                  + f_52 * dl_68[k]
                  - f_54 * dl_70[k]
                  + f_49 * dl_81[k]
                  - f_50 * dl_83[k]
                  + f_53 * dl_85[k];
    }

#pragma omp simd aligned(dl_47, dl_52, dl_54, dl_61, dl_63, dl_74, \
                         dl_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_14 * dl_47[k]
                  + f_12 * dl_52[k]
                  + f_15 * dl_54[k]
                  + f_10 * dl_61[k]
                  - f_13 * dl_63[k]
                  - f_10 * dl_74[k]
                  + f_11 * dl_76[k];
    }

#pragma omp simd aligned(dl_45, dl_47, dl_48, dl_50, dl_52, dl_55, dl_57, dl_61, dl_66, dl_68, \
                         dl_74, dl_81, dl_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_55 * dl_45[k]
                  + f_7 * dl_48[k]
                  + f_7 * dl_50[k]
                  - f_56 * dl_57[k]
                  - f_7 * dl_66[k]
                  + f_56 * dl_68[k]
                  + f_55 * dl_81[k]
                  - f_7 * dl_83[k];

        g_15[k] = f_5 * dl_47[k]
                  - f_4 * dl_52[k]
                  + f_3 * dl_61[k]
                  - f_2 * dl_74[k];

        g_16[k] = f_57 * dl_45[k]
                  - f_2 * dl_48[k]
                  + f_58 * dl_55[k]
                  - f_2 * dl_66[k]
                  + f_57 * dl_81[k];
    }

#pragma omp simd aligned(dl_181, dl_184, dl_186, dl_188, dl_191, dl_195, dl_197, dl_202, \
                         dl_208, dl_210, dl_217 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_0 * dl_181[k]
                  - f_1 * dl_186[k]
                  + f_1 * dl_195[k]
                  - f_0 * dl_208[k];

        g_18[k] = f_2 * dl_184[k]
                  - f_3 * dl_191[k]
                  + f_4 * dl_202[k]
                  - f_5 * dl_217[k];

        g_19[k] = -f_6 * dl_181[k]
                  + f_7 * dl_186[k]
                  + f_8 * dl_188[k]
                  + f_7 * dl_195[k]
                  - f_9 * dl_197[k]
                  - f_6 * dl_208[k]
                  + f_8 * dl_210[k];
    }

#pragma omp simd aligned(dl_184, dl_191, dl_193, dl_202, dl_204, dl_217, \
                         dl_219 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_10 * dl_184[k]
                  + f_10 * dl_191[k]
                  + f_11 * dl_193[k]
                  + f_12 * dl_202[k]
                  - f_13 * dl_204[k]
                  - f_14 * dl_217[k]
                  + f_15 * dl_219[k];
    }

#pragma omp simd aligned(dl_181, dl_186, dl_188, dl_195, dl_199, dl_208, dl_210, \
                         dl_212 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_16 * dl_181[k]
                  + f_16 * dl_186[k]
                  - f_17 * dl_188[k]
                  - f_16 * dl_195[k]
                  + f_18 * dl_199[k]
                  - f_16 * dl_208[k]
                  + f_17 * dl_210[k]
                  - f_18 * dl_212[k];
    }

#pragma omp simd aligned(dl_184, dl_191, dl_193, dl_202, dl_204, dl_206, dl_217, dl_219, \
                         dl_221 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_19 * dl_184[k]
                  + f_20 * dl_191[k]
                  - f_21 * dl_193[k]
                  + f_22 * dl_202[k]
                  - f_23 * dl_204[k]
                  + f_24 * dl_206[k]
                  - f_22 * dl_217[k]
                  + f_25 * dl_219[k]
                  - f_26 * dl_221[k];
    }

#pragma omp simd aligned(dl_181, dl_186, dl_188, dl_195, dl_197, dl_199, dl_208, dl_210, \
                         dl_212, dl_214 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_27 * dl_181[k]
                  - f_28 * dl_186[k]
                  + f_29 * dl_188[k]
                  - f_28 * dl_195[k]
                  + f_30 * dl_197[k]
                  - f_31 * dl_199[k]
                  - f_27 * dl_208[k]
                  + f_29 * dl_210[k]
                  - f_31 * dl_212[k]
                  + f_32 * dl_214[k];
    }

#pragma omp simd aligned(dl_184, dl_191, dl_193, dl_202, dl_204, dl_206, dl_217, dl_219, \
                         dl_221, dl_223 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -f_33 * dl_184[k]
                  - f_34 * dl_191[k]
                  + f_35 * dl_193[k]
                  - f_34 * dl_202[k]
                  + f_36 * dl_204[k]
                  - f_37 * dl_206[k]
                  - f_33 * dl_217[k]
                  + f_35 * dl_219[k]
                  - f_37 * dl_221[k]
                  + f_38 * dl_223[k];
    }

#pragma omp simd aligned(dl_180, dl_183, dl_185, dl_190, dl_192, dl_194, dl_201, dl_203, \
                         dl_205, dl_207, dl_216, dl_218, dl_220, dl_222, \
                         dl_224 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_39 * dl_180[k]
                  + f_40 * dl_183[k]
                  - f_41 * dl_185[k]
                  + f_42 * dl_190[k]
                  - f_35 * dl_192[k]
                  + f_35 * dl_194[k]
                  + f_40 * dl_201[k]
                  - f_35 * dl_203[k]
                  + f_36 * dl_205[k]
                  - f_43 * dl_207[k]
                  + f_39 * dl_216[k]
                  - f_41 * dl_218[k]
                  + f_35 * dl_220[k]
                  - f_43 * dl_222[k]
                  + f_44 * dl_224[k];
    }

#pragma omp simd aligned(dl_182, dl_187, dl_189, dl_196, dl_198, dl_200, dl_209, dl_211, \
                         dl_213, dl_215 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_33 * dl_182[k]
                  - f_34 * dl_187[k]
                  + f_35 * dl_189[k]
                  - f_34 * dl_196[k]
                  + f_36 * dl_198[k]
                  - f_37 * dl_200[k]
                  - f_33 * dl_209[k]
                  + f_35 * dl_211[k]
                  - f_37 * dl_213[k]
                  + f_38 * dl_215[k];
    }

#pragma omp simd aligned(dl_180, dl_183, dl_185, dl_192, dl_194, dl_201, dl_203, dl_207, \
                         dl_216, dl_218, dl_220, dl_222 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_45 * dl_180[k]
                  - f_27 * dl_183[k]
                  + f_46 * dl_185[k]
                  + f_46 * dl_192[k]
                  - f_47 * dl_194[k]
                  + f_27 * dl_201[k]
                  - f_46 * dl_203[k]
                  + f_48 * dl_207[k]
                  + f_45 * dl_216[k]
                  - f_46 * dl_218[k]
                  + f_47 * dl_220[k]
                  - f_48 * dl_222[k];
    }

#pragma omp simd aligned(dl_182, dl_187, dl_189, dl_196, dl_198, dl_200, dl_209, dl_211, \
                         dl_213 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_22 * dl_182[k]
                  - f_22 * dl_187[k]
                  - f_25 * dl_189[k]
                  - f_20 * dl_196[k]
                  + f_23 * dl_198[k]
                  + f_26 * dl_200[k]
                  - f_19 * dl_209[k]
                  + f_21 * dl_211[k]
                  - f_24 * dl_213[k];
    }

#pragma omp simd aligned(dl_180, dl_183, dl_185, dl_190, dl_192, dl_194, dl_201, dl_203, \
                         dl_205, dl_216, dl_218, dl_220 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_49 * dl_180[k]
                  - f_16 * dl_183[k]
                  - f_50 * dl_185[k]
                  - f_51 * dl_190[k]
                  + f_52 * dl_192[k]
                  + f_53 * dl_194[k]
                  - f_16 * dl_201[k]
                  + f_52 * dl_203[k]
                  - f_54 * dl_205[k]
                  + f_49 * dl_216[k]
                  - f_50 * dl_218[k]
                  + f_53 * dl_220[k];
    }

#pragma omp simd aligned(dl_182, dl_187, dl_189, dl_196, dl_198, dl_209, \
                         dl_211 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_14 * dl_182[k]
                  + f_12 * dl_187[k]
                  + f_15 * dl_189[k]
                  + f_10 * dl_196[k]
                  - f_13 * dl_198[k]
                  - f_10 * dl_209[k]
                  + f_11 * dl_211[k];
    }

#pragma omp simd aligned(dl_180, dl_182, dl_183, dl_185, dl_187, dl_190, dl_192, dl_196, \
                         dl_201, dl_203, dl_209, dl_216, dl_218 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_55 * dl_180[k]
                  + f_7 * dl_183[k]
                  + f_7 * dl_185[k]
                  - f_56 * dl_192[k]
                  - f_7 * dl_201[k]
                  + f_56 * dl_203[k]
                  + f_55 * dl_216[k]
                  - f_7 * dl_218[k];

        g_32[k] = f_5 * dl_182[k]
                  - f_4 * dl_187[k]
                  + f_3 * dl_196[k]
                  - f_2 * dl_209[k];

        g_33[k] = f_57 * dl_180[k]
                  - f_2 * dl_183[k]
                  + f_58 * dl_190[k]
                  - f_2 * dl_201[k]
                  + f_57 * dl_216[k];
    }

#pragma omp simd aligned(dl_1, dl_6, dl_15, dl_28, dl_136, dl_141, dl_150, dl_163, dl_226, \
                         dl_231, dl_240, dl_253 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_59 * dl_1[k]
                  + f_60 * dl_6[k]
                  - f_60 * dl_15[k]
                  + f_59 * dl_28[k]
                  - f_59 * dl_136[k]
                  + f_60 * dl_141[k]
                  - f_60 * dl_150[k]
                  + f_59 * dl_163[k]
                  + f_61 * dl_226[k]
                  - f_62 * dl_231[k]
                  + f_62 * dl_240[k]
                  - f_61 * dl_253[k];
    }

#pragma omp simd aligned(dl_4, dl_11, dl_22, dl_37, dl_139, dl_146, dl_157, dl_172, dl_229, \
                         dl_236, dl_247, dl_262 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_63 * dl_4[k]
                  + f_64 * dl_11[k]
                  - f_65 * dl_22[k]
                  + f_66 * dl_37[k]
                  - f_63 * dl_139[k]
                  + f_64 * dl_146[k]
                  - f_65 * dl_157[k]
                  + f_66 * dl_172[k]
                  + f_60 * dl_229[k]
                  - f_67 * dl_236[k]
                  + f_68 * dl_247[k]
                  - f_59 * dl_262[k];
    }

#pragma omp simd aligned(dl_1, dl_6, dl_8, dl_15, dl_17, dl_28, dl_30, dl_136, dl_141, dl_143, \
                         dl_150, dl_152, dl_163, dl_165, dl_226, dl_231, dl_233, dl_240, \
                         dl_242, dl_253, dl_255 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_69 * dl_1[k]
                  - f_70 * dl_6[k]
                  - f_71 * dl_8[k]
                  - f_70 * dl_15[k]
                  + f_72 * dl_17[k]
                  + f_69 * dl_28[k]
                  - f_71 * dl_30[k]
                  + f_69 * dl_136[k]
                  - f_70 * dl_141[k]
                  - f_71 * dl_143[k]
                  - f_70 * dl_150[k]
                  + f_72 * dl_152[k]
                  + f_69 * dl_163[k]
                  - f_71 * dl_165[k]
                  - f_73 * dl_226[k]
                  + f_74 * dl_231[k]
                  + f_75 * dl_233[k]
                  + f_74 * dl_240[k]
                  - f_76 * dl_242[k]
                  - f_73 * dl_253[k]
                  + f_75 * dl_255[k];
    }

#pragma omp simd aligned(dl_4, dl_11, dl_13, dl_22, dl_24, dl_37, dl_39, dl_139, dl_146, \
                         dl_148, dl_157, dl_159, dl_172, dl_174, dl_229, dl_236, dl_238, \
                         dl_247, dl_249, dl_262, dl_264 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_77 * dl_4[k]
                  - f_77 * dl_11[k]
                  - f_78 * dl_13[k]
                  - f_79 * dl_22[k]
                  + f_80 * dl_24[k]
                  + f_81 * dl_37[k]
                  - f_82 * dl_39[k]
                  + f_77 * dl_139[k]
                  - f_77 * dl_146[k]
                  - f_78 * dl_148[k]
                  - f_79 * dl_157[k]
                  + f_80 * dl_159[k]
                  + f_81 * dl_172[k]
                  - f_82 * dl_174[k]
                  - f_83 * dl_229[k]
                  + f_83 * dl_236[k]
                  + f_80 * dl_238[k]
                  + f_84 * dl_247[k]
                  - f_85 * dl_249[k]
                  - f_86 * dl_262[k]
                  + f_87 * dl_264[k];
    }

#pragma omp simd aligned(dl_1, dl_6, dl_8, dl_15, dl_19, dl_28, dl_30, dl_32, dl_136, dl_141, \
                         dl_143, dl_150, dl_154, dl_163, dl_165, dl_167, dl_226, dl_231, \
                         dl_233, dl_240, dl_244, dl_253, dl_255, \
                         dl_257 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_88 * dl_1[k]
                  - f_88 * dl_6[k]
                  + f_89 * dl_8[k]
                  + f_88 * dl_15[k]
                  - f_90 * dl_19[k]
                  + f_88 * dl_28[k]
                  - f_89 * dl_30[k]
                  + f_90 * dl_32[k]
                  - f_88 * dl_136[k]
                  - f_88 * dl_141[k]
                  + f_89 * dl_143[k]
                  + f_88 * dl_150[k]
                  - f_90 * dl_154[k]
                  + f_88 * dl_163[k]
                  - f_89 * dl_165[k]
                  + f_90 * dl_167[k]
                  + f_91 * dl_226[k]
                  + f_91 * dl_231[k]
                  - f_92 * dl_233[k]
                  - f_91 * dl_240[k]
                  + f_93 * dl_244[k]
                  - f_91 * dl_253[k]
                  + f_92 * dl_255[k]
                  - f_93 * dl_257[k];
    }

#pragma omp simd aligned(dl_4, dl_11, dl_13, dl_22, dl_24, dl_26, dl_37, dl_39, dl_41, dl_139, \
                         dl_146, dl_148, dl_157, dl_159, dl_161, dl_172, dl_174, dl_176, \
                         dl_229, dl_236, dl_238, dl_247, dl_249, dl_251, dl_262, dl_264, \
                         dl_266 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_94 * dl_4[k]
                  - f_95 * dl_11[k]
                  + f_96 * dl_13[k]
                  - f_97 * dl_22[k]
                  + f_98 * dl_24[k]
                  - f_99 * dl_26[k]
                  + f_97 * dl_37[k]
                  - f_100 * dl_39[k]
                  + f_101 * dl_41[k]
                  - f_94 * dl_139[k]
                  - f_95 * dl_146[k]
                  + f_96 * dl_148[k]
                  - f_97 * dl_157[k]
                  + f_98 * dl_159[k]
                  - f_99 * dl_161[k]
                  + f_97 * dl_172[k]
                  - f_100 * dl_174[k]
                  + f_101 * dl_176[k]
                  + f_102 * dl_229[k]
                  + f_103 * dl_236[k]
                  - f_104 * dl_238[k]
                  + f_105 * dl_247[k]
                  - f_106 * dl_249[k]
                  + f_107 * dl_251[k]
                  - f_105 * dl_262[k]
                  + f_98 * dl_264[k]
                  - f_108 * dl_266[k];
    }

#pragma omp simd aligned(dl_1, dl_6, dl_8, dl_15, dl_17, dl_19, dl_28, dl_30, dl_32, dl_34, \
                         dl_136, dl_141, dl_143, dl_150, dl_152, dl_154, dl_163, dl_165, \
                         dl_167, dl_169, dl_226, dl_231, dl_233, dl_240, dl_242, dl_244, \
                         dl_253, dl_255, dl_257, dl_259 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_109 * dl_1[k]
                  + f_110 * dl_6[k]
                  - f_111 * dl_8[k]
                  + f_110 * dl_15[k]
                  - f_112 * dl_17[k]
                  + f_113 * dl_19[k]
                  + f_109 * dl_28[k]
                  - f_111 * dl_30[k]
                  + f_113 * dl_32[k]
                  - f_114 * dl_34[k]
                  + f_109 * dl_136[k]
                  + f_110 * dl_141[k]
                  - f_111 * dl_143[k]
                  + f_110 * dl_150[k]
                  - f_112 * dl_152[k]
                  + f_113 * dl_154[k]
                  + f_109 * dl_163[k]
                  - f_111 * dl_165[k]
                  + f_113 * dl_167[k]
                  - f_114 * dl_169[k]
                  - f_115 * dl_226[k]
                  - f_116 * dl_231[k]
                  + f_112 * dl_233[k]
                  - f_116 * dl_240[k]
                  + f_117 * dl_242[k]
                  - f_118 * dl_244[k]
                  - f_115 * dl_253[k]
                  + f_112 * dl_255[k]
                  - f_118 * dl_257[k]
                  + f_119 * dl_259[k];
    }

#pragma omp simd aligned(dl_4, dl_11, dl_13, dl_22, dl_24, dl_26, dl_37, dl_39, dl_41, dl_43, \
                         dl_139, dl_146, dl_148, dl_157, dl_159, dl_161, dl_172, dl_174, \
                         dl_176, dl_178, dl_229, dl_236, dl_238, dl_247, dl_249, dl_251, \
                         dl_262, dl_264, dl_266, dl_268 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = 1.640625 * dl_4[k]
                  + 4.921875 * dl_11[k]
                  - 13.125 * dl_13[k]
                  + 4.921875 * dl_22[k]
                  - 26.25 * dl_24[k]
                  + 15.75 * dl_26[k]
                  + 1.640625 * dl_37[k]
                  - 13.125 * dl_39[k]
                  + 15.75 * dl_41[k]
                  - 3.0 * dl_43[k]
                  + 1.640625 * dl_139[k]
                  + 4.921875 * dl_146[k]
                  - 13.125 * dl_148[k]
                  + 4.921875 * dl_157[k]
                  - 26.25 * dl_159[k]
                  + 15.75 * dl_161[k]
                  + 1.640625 * dl_172[k]
                  - 13.125 * dl_174[k]
                  + 15.75 * dl_176[k]
                  - 3.0 * dl_178[k]
                  - 3.28125 * dl_229[k]
                  - 9.84375 * dl_236[k]
                  + 26.25 * dl_238[k]
                  - 9.84375 * dl_247[k]
                  + 52.5 * dl_249[k]
                  - 31.5 * dl_251[k]
                  - 3.28125 * dl_262[k]
                  + 26.25 * dl_264[k]
                  - 31.5 * dl_266[k]
                  + 6.0 * dl_268[k];
    }

#pragma omp simd aligned(dl_0, dl_3, dl_5, dl_10, dl_12, dl_14, dl_21, dl_23, dl_25, dl_27, \
                         dl_36, dl_38, dl_40, dl_42, dl_44, dl_135, dl_138, dl_140, dl_145, \
                         dl_147, dl_149, dl_156, dl_158, dl_160, dl_162, dl_171, dl_173, \
                         dl_175, dl_177, dl_179, dl_225, dl_228, dl_230, dl_235, dl_237, \
                         dl_239, dl_246, dl_248, dl_250, dl_252, dl_261, dl_263, dl_265, \
                         dl_267, dl_269 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -0.13671875 * dl_0[k]
                  - 0.546875 * dl_3[k]
                  + 4.375 * dl_5[k]
                  - 0.8203125 * dl_10[k]
                  + 13.125 * dl_12[k]
                  - 13.125 * dl_14[k]
                  - 0.546875 * dl_21[k]
                  + 13.125 * dl_23[k]
                  - 26.25 * dl_25[k]
                  + 7.0 * dl_27[k]
                  - 0.13671875 * dl_36[k]
                  + 4.375 * dl_38[k]
                  - 13.125 * dl_40[k]
                  + 7.0 * dl_42[k]
                  - 0.5 * dl_44[k]
                  - 0.13671875 * dl_135[k]
                  - 0.546875 * dl_138[k]
                  + 4.375 * dl_140[k]
                  - 0.8203125 * dl_145[k]
                  + 13.125 * dl_147[k]
                  - 13.125 * dl_149[k]
                  - 0.546875 * dl_156[k]
                  + 13.125 * dl_158[k]
                  - 26.25 * dl_160[k]
                  + 7.0 * dl_162[k]
                  - 0.13671875 * dl_171[k]
                  + 4.375 * dl_173[k]
                  - 13.125 * dl_175[k]
                  + 7.0 * dl_177[k]
                  - 0.5 * dl_179[k]
                  + 0.2734375 * dl_225[k]
                  + 1.09375 * dl_228[k]
                  - 8.75 * dl_230[k]
                  + 1.640625 * dl_235[k]
                  - 26.25 * dl_237[k]
                  + 26.25 * dl_239[k]
                  + 1.09375 * dl_246[k]
                  - 26.25 * dl_248[k]
                  + 52.5 * dl_250[k]
                  - 14.0 * dl_252[k]
                  + 0.2734375 * dl_261[k]
                  - 8.75 * dl_263[k]
                  + 26.25 * dl_265[k]
                  - 14.0 * dl_267[k]
                  + dl_269[k];
    }

#pragma omp simd aligned(dl_2, dl_7, dl_9, dl_16, dl_18, dl_20, dl_29, dl_31, dl_33, dl_35, \
                         dl_137, dl_142, dl_144, dl_151, dl_153, dl_155, dl_164, dl_166, \
                         dl_168, dl_170, dl_227, dl_232, dl_234, dl_241, dl_243, dl_245, \
                         dl_254, dl_256, dl_258, dl_260 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = 1.640625 * dl_2[k]
                  + 4.921875 * dl_7[k]
                  - 13.125 * dl_9[k]
                  + 4.921875 * dl_16[k]
                  - 26.25 * dl_18[k]
                  + 15.75 * dl_20[k]
                  + 1.640625 * dl_29[k]
                  - 13.125 * dl_31[k]
                  + 15.75 * dl_33[k]
                  - 3.0 * dl_35[k]
                  + 1.640625 * dl_137[k]
                  + 4.921875 * dl_142[k]
                  - 13.125 * dl_144[k]
                  + 4.921875 * dl_151[k]
                  - 26.25 * dl_153[k]
                  + 15.75 * dl_155[k]
                  + 1.640625 * dl_164[k]
                  - 13.125 * dl_166[k]
                  + 15.75 * dl_168[k]
                  - 3.0 * dl_170[k]
                  - 3.28125 * dl_227[k]
                  - 9.84375 * dl_232[k]
                  + 26.25 * dl_234[k]
                  - 9.84375 * dl_241[k]
                  + 52.5 * dl_243[k]
                  - 31.5 * dl_245[k]
                  - 3.28125 * dl_254[k]
                  + 26.25 * dl_256[k]
                  - 31.5 * dl_258[k]
                  + 6.0 * dl_260[k];
    }

#pragma omp simd aligned(dl_0, dl_3, dl_5, dl_12, dl_14, dl_21, dl_23, dl_27, dl_36, dl_38, \
                         dl_40, dl_42, dl_135, dl_138, dl_140, dl_147, dl_149, dl_156, dl_158, \
                         dl_162, dl_171, dl_173, dl_175, dl_177, dl_225, dl_228, dl_230, \
                         dl_237, dl_239, dl_246, dl_248, dl_252, dl_261, dl_263, dl_265, \
                         dl_267 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_120 * dl_0[k]
                  + f_109 * dl_3[k]
                  - f_121 * dl_5[k]
                  - f_121 * dl_12[k]
                  + f_122 * dl_14[k]
                  - f_109 * dl_21[k]
                  + f_121 * dl_23[k]
                  - f_123 * dl_27[k]
                  - f_120 * dl_36[k]
                  + f_121 * dl_38[k]
                  - f_122 * dl_40[k]
                  + f_123 * dl_42[k]
                  + f_120 * dl_135[k]
                  + f_109 * dl_138[k]
                  - f_121 * dl_140[k]
                  - f_121 * dl_147[k]
                  + f_122 * dl_149[k]
                  - f_109 * dl_156[k]
                  + f_121 * dl_158[k]
                  - f_123 * dl_162[k]
                  - f_120 * dl_171[k]
                  + f_121 * dl_173[k]
                  - f_122 * dl_175[k]
                  + f_123 * dl_177[k]
                  - f_109 * dl_225[k]
                  - f_115 * dl_228[k]
                  + f_111 * dl_230[k]
                  + f_111 * dl_237[k]
                  - f_113 * dl_239[k]
                  + f_115 * dl_246[k]
                  - f_111 * dl_248[k]
                  + f_114 * dl_252[k]
                  + f_109 * dl_261[k]
                  - f_111 * dl_263[k]
                  + f_113 * dl_265[k]
                  - f_114 * dl_267[k];
    }

#pragma omp simd aligned(dl_2, dl_7, dl_9, dl_16, dl_18, dl_20, dl_29, dl_31, dl_33, dl_137, \
                         dl_142, dl_144, dl_151, dl_153, dl_155, dl_164, dl_166, dl_168, \
                         dl_227, dl_232, dl_234, dl_241, dl_243, dl_245, dl_254, dl_256, \
                         dl_258 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_97 * dl_2[k]
                  + f_97 * dl_7[k]
                  + f_100 * dl_9[k]
                  + f_95 * dl_16[k]
                  - f_98 * dl_18[k]
                  - f_101 * dl_20[k]
                  + f_94 * dl_29[k]
                  - f_96 * dl_31[k]
                  + f_99 * dl_33[k]
                  - f_97 * dl_137[k]
                  + f_97 * dl_142[k]
                  + f_100 * dl_144[k]
                  + f_95 * dl_151[k]
                  - f_98 * dl_153[k]
                  - f_101 * dl_155[k]
                  + f_94 * dl_164[k]
                  - f_96 * dl_166[k]
                  + f_99 * dl_168[k]
                  + f_105 * dl_227[k]
                  - f_105 * dl_232[k]
                  - f_98 * dl_234[k]
                  - f_103 * dl_241[k]
                  + f_106 * dl_243[k]
                  + f_108 * dl_245[k]
                  - f_102 * dl_254[k]
                  + f_104 * dl_256[k]
                  - f_107 * dl_258[k];
    }

#pragma omp simd aligned(dl_0, dl_3, dl_5, dl_10, dl_12, dl_14, dl_21, dl_23, dl_25, dl_36, \
                         dl_38, dl_40, dl_135, dl_138, dl_140, dl_145, dl_147, dl_149, dl_156, \
                         dl_158, dl_160, dl_171, dl_173, dl_175, dl_225, dl_228, dl_230, \
                         dl_235, dl_237, dl_239, dl_246, dl_248, dl_250, dl_261, dl_263, \
                         dl_265 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_124 * dl_0[k]
                  + f_88 * dl_3[k]
                  + f_125 * dl_5[k]
                  + f_126 * dl_10[k]
                  - f_127 * dl_12[k]
                  - f_128 * dl_14[k]
                  + f_88 * dl_21[k]
                  - f_127 * dl_23[k]
                  + f_129 * dl_25[k]
                  - f_124 * dl_36[k]
                  + f_125 * dl_38[k]
                  - f_128 * dl_40[k]
                  - f_124 * dl_135[k]
                  + f_88 * dl_138[k]
                  + f_125 * dl_140[k]
                  + f_126 * dl_145[k]
                  - f_127 * dl_147[k]
                  - f_128 * dl_149[k]
                  + f_88 * dl_156[k]
                  - f_127 * dl_158[k]
                  + f_129 * dl_160[k]
                  - f_124 * dl_171[k]
                  + f_125 * dl_173[k]
                  - f_128 * dl_175[k]
                  + f_130 * dl_225[k]
                  - f_91 * dl_228[k]
                  - f_131 * dl_230[k]
                  - f_132 * dl_235[k]
                  + f_129 * dl_237[k]
                  + f_133 * dl_239[k]
                  - f_91 * dl_246[k]
                  + f_129 * dl_248[k]
                  - f_134 * dl_250[k]
                  + f_130 * dl_261[k]
                  - f_131 * dl_263[k]
                  + f_133 * dl_265[k];
    }

#pragma omp simd aligned(dl_2, dl_7, dl_9, dl_16, dl_18, dl_29, dl_31, dl_137, dl_142, dl_144, \
                         dl_151, dl_153, dl_164, dl_166, dl_227, dl_232, dl_234, dl_241, \
                         dl_243, dl_254, dl_256 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = f_81 * dl_2[k]
                  - f_79 * dl_7[k]
                  - f_82 * dl_9[k]
                  - f_77 * dl_16[k]
                  + f_80 * dl_18[k]
                  + f_77 * dl_29[k]
                  - f_78 * dl_31[k]
                  + f_81 * dl_137[k]
                  - f_79 * dl_142[k]
                  - f_82 * dl_144[k]
                  - f_77 * dl_151[k]
                  + f_80 * dl_153[k]
                  + f_77 * dl_164[k]
                  - f_78 * dl_166[k]
                  - f_86 * dl_227[k]
                  + f_84 * dl_232[k]
                  + f_87 * dl_234[k]
                  + f_83 * dl_241[k]
                  - f_85 * dl_243[k]
                  - f_83 * dl_254[k]
                  + f_80 * dl_256[k];
    }

#pragma omp simd aligned(dl_0, dl_3, dl_5, dl_12, dl_21, dl_23, dl_36, dl_38, dl_135, dl_138, \
                         dl_140, dl_147, dl_156, dl_158, dl_171, dl_173, dl_225, dl_228, \
                         dl_230, dl_237, dl_246, dl_248, dl_261, \
                         dl_263 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_135 * dl_0[k]
                  - f_70 * dl_3[k]
                  - f_70 * dl_5[k]
                  + f_136 * dl_12[k]
                  + f_70 * dl_21[k]
                  - f_136 * dl_23[k]
                  - f_135 * dl_36[k]
                  + f_70 * dl_38[k]
                  + f_135 * dl_135[k]
                  - f_70 * dl_138[k]
                  - f_70 * dl_140[k]
                  + f_136 * dl_147[k]
                  + f_70 * dl_156[k]
                  - f_136 * dl_158[k]
                  - f_135 * dl_171[k]
                  + f_70 * dl_173[k]
                  - f_137 * dl_225[k]
                  + f_74 * dl_228[k]
                  + f_74 * dl_230[k]
                  - f_138 * dl_237[k]
                  - f_74 * dl_246[k]
                  + f_138 * dl_248[k]
                  + f_137 * dl_261[k]
                  - f_74 * dl_263[k];
    }

#pragma omp simd aligned(dl_2, dl_7, dl_16, dl_29, dl_137, dl_142, dl_151, dl_164, dl_227, \
                         dl_232, dl_241, dl_254 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_66 * dl_2[k]
                  + f_65 * dl_7[k]
                  - f_64 * dl_16[k]
                  + f_63 * dl_29[k]
                  - f_66 * dl_137[k]
                  + f_65 * dl_142[k]
                  - f_64 * dl_151[k]
                  + f_63 * dl_164[k]
                  + f_59 * dl_227[k]
                  - f_68 * dl_232[k]
                  + f_67 * dl_241[k]
                  - f_60 * dl_254[k];
    }

#pragma omp simd aligned(dl_0, dl_3, dl_10, dl_21, dl_36, dl_135, dl_138, dl_145, dl_156, \
                         dl_171, dl_225, dl_228, dl_235, dl_246, \
                         dl_261 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_139 * dl_0[k]
                  + f_63 * dl_3[k]
                  - f_140 * dl_10[k]
                  + f_63 * dl_21[k]
                  - f_139 * dl_36[k]
                  - f_139 * dl_135[k]
                  + f_63 * dl_138[k]
                  - f_140 * dl_145[k]
                  + f_63 * dl_156[k]
                  - f_139 * dl_171[k]
                  + f_141 * dl_225[k]
                  - f_60 * dl_228[k]
                  + f_64 * dl_235[k]
                  - f_60 * dl_246[k]
                  + f_141 * dl_261[k];
    }

#pragma omp simd aligned(dl_91, dl_94, dl_96, dl_98, dl_101, dl_105, dl_107, dl_112, dl_118, \
                         dl_120, dl_127 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_0 * dl_91[k]
                  - f_1 * dl_96[k]
                  + f_1 * dl_105[k]
                  - f_0 * dl_118[k];

        g_52[k] = f_2 * dl_94[k]
                  - f_3 * dl_101[k]
                  + f_4 * dl_112[k]
                  - f_5 * dl_127[k];

        g_53[k] = -f_6 * dl_91[k]
                  + f_7 * dl_96[k]
                  + f_8 * dl_98[k]
                  + f_7 * dl_105[k]
                  - f_9 * dl_107[k]
                  - f_6 * dl_118[k]
                  + f_8 * dl_120[k];
    }

#pragma omp simd aligned(dl_94, dl_101, dl_103, dl_112, dl_114, dl_127, \
                         dl_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_10 * dl_94[k]
                  + f_10 * dl_101[k]
                  + f_11 * dl_103[k]
                  + f_12 * dl_112[k]
                  - f_13 * dl_114[k]
                  - f_14 * dl_127[k]
                  + f_15 * dl_129[k];
    }

#pragma omp simd aligned(dl_91, dl_96, dl_98, dl_105, dl_109, dl_118, dl_120, \
                         dl_122 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_16 * dl_91[k]
                  + f_16 * dl_96[k]
                  - f_17 * dl_98[k]
                  - f_16 * dl_105[k]
                  + f_18 * dl_109[k]
                  - f_16 * dl_118[k]
                  + f_17 * dl_120[k]
                  - f_18 * dl_122[k];
    }

#pragma omp simd aligned(dl_94, dl_101, dl_103, dl_112, dl_114, dl_116, dl_127, dl_129, \
                         dl_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_19 * dl_94[k]
                  + f_20 * dl_101[k]
                  - f_21 * dl_103[k]
                  + f_22 * dl_112[k]
                  - f_23 * dl_114[k]
                  + f_24 * dl_116[k]
                  - f_22 * dl_127[k]
                  + f_25 * dl_129[k]
                  - f_26 * dl_131[k];
    }

#pragma omp simd aligned(dl_91, dl_96, dl_98, dl_105, dl_107, dl_109, dl_118, dl_120, dl_122, \
                         dl_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = -f_27 * dl_91[k]
                  - f_28 * dl_96[k]
                  + f_29 * dl_98[k]
                  - f_28 * dl_105[k]
                  + f_30 * dl_107[k]
                  - f_31 * dl_109[k]
                  - f_27 * dl_118[k]
                  + f_29 * dl_120[k]
                  - f_31 * dl_122[k]
                  + f_32 * dl_124[k];
    }

#pragma omp simd aligned(dl_94, dl_101, dl_103, dl_112, dl_114, dl_116, dl_127, dl_129, \
                         dl_131, dl_133 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_33 * dl_94[k]
                  - f_34 * dl_101[k]
                  + f_35 * dl_103[k]
                  - f_34 * dl_112[k]
                  + f_36 * dl_114[k]
                  - f_37 * dl_116[k]
                  - f_33 * dl_127[k]
                  + f_35 * dl_129[k]
                  - f_37 * dl_131[k]
                  + f_38 * dl_133[k];
    }

#pragma omp simd aligned(dl_90, dl_93, dl_95, dl_100, dl_102, dl_104, dl_111, dl_113, dl_115, \
                         dl_117, dl_126, dl_128, dl_130, dl_132, \
                         dl_134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_39 * dl_90[k]
                  + f_40 * dl_93[k]
                  - f_41 * dl_95[k]
                  + f_42 * dl_100[k]
                  - f_35 * dl_102[k]
                  + f_35 * dl_104[k]
                  + f_40 * dl_111[k]
                  - f_35 * dl_113[k]
                  + f_36 * dl_115[k]
                  - f_43 * dl_117[k]
                  + f_39 * dl_126[k]
                  - f_41 * dl_128[k]
                  + f_35 * dl_130[k]
                  - f_43 * dl_132[k]
                  + f_44 * dl_134[k];
    }

#pragma omp simd aligned(dl_92, dl_97, dl_99, dl_106, dl_108, dl_110, dl_119, dl_121, dl_123, \
                         dl_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -f_33 * dl_92[k]
                  - f_34 * dl_97[k]
                  + f_35 * dl_99[k]
                  - f_34 * dl_106[k]
                  + f_36 * dl_108[k]
                  - f_37 * dl_110[k]
                  - f_33 * dl_119[k]
                  + f_35 * dl_121[k]
                  - f_37 * dl_123[k]
                  + f_38 * dl_125[k];
    }

#pragma omp simd aligned(dl_90, dl_93, dl_95, dl_102, dl_104, dl_111, dl_113, dl_117, dl_126, \
                         dl_128, dl_130, dl_132 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_45 * dl_90[k]
                  - f_27 * dl_93[k]
                  + f_46 * dl_95[k]
                  + f_46 * dl_102[k]
                  - f_47 * dl_104[k]
                  + f_27 * dl_111[k]
                  - f_46 * dl_113[k]
                  + f_48 * dl_117[k]
                  + f_45 * dl_126[k]
                  - f_46 * dl_128[k]
                  + f_47 * dl_130[k]
                  - f_48 * dl_132[k];
    }

#pragma omp simd aligned(dl_92, dl_97, dl_99, dl_106, dl_108, dl_110, dl_119, dl_121, \
                         dl_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_22 * dl_92[k]
                  - f_22 * dl_97[k]
                  - f_25 * dl_99[k]
                  - f_20 * dl_106[k]
                  + f_23 * dl_108[k]
                  + f_26 * dl_110[k]
                  - f_19 * dl_119[k]
                  + f_21 * dl_121[k]
                  - f_24 * dl_123[k];
    }

#pragma omp simd aligned(dl_90, dl_93, dl_95, dl_100, dl_102, dl_104, dl_111, dl_113, dl_115, \
                         dl_126, dl_128, dl_130 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_49 * dl_90[k]
                  - f_16 * dl_93[k]
                  - f_50 * dl_95[k]
                  - f_51 * dl_100[k]
                  + f_52 * dl_102[k]
                  + f_53 * dl_104[k]
                  - f_16 * dl_111[k]
                  + f_52 * dl_113[k]
                  - f_54 * dl_115[k]
                  + f_49 * dl_126[k]
                  - f_50 * dl_128[k]
                  + f_53 * dl_130[k];
    }

#pragma omp simd aligned(dl_92, dl_97, dl_99, dl_106, dl_108, dl_119, \
                         dl_121 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_14 * dl_92[k]
                  + f_12 * dl_97[k]
                  + f_15 * dl_99[k]
                  + f_10 * dl_106[k]
                  - f_13 * dl_108[k]
                  - f_10 * dl_119[k]
                  + f_11 * dl_121[k];
    }

#pragma omp simd aligned(dl_90, dl_92, dl_93, dl_95, dl_97, dl_100, dl_102, dl_106, dl_111, \
                         dl_113, dl_119, dl_126, dl_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = -f_55 * dl_90[k]
                  + f_7 * dl_93[k]
                  + f_7 * dl_95[k]
                  - f_56 * dl_102[k]
                  - f_7 * dl_111[k]
                  + f_56 * dl_113[k]
                  + f_55 * dl_126[k]
                  - f_7 * dl_128[k];

        g_66[k] = f_5 * dl_92[k]
                  - f_4 * dl_97[k]
                  + f_3 * dl_106[k]
                  - f_2 * dl_119[k];

        g_67[k] = f_57 * dl_90[k]
                  - f_2 * dl_93[k]
                  + f_58 * dl_100[k]
                  - f_2 * dl_111[k]
                  + f_57 * dl_126[k];
    }

#pragma omp simd aligned(dl_1, dl_6, dl_15, dl_28, dl_136, dl_141, dl_150, \
                         dl_163 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_5 * dl_1[k]
                  - f_2 * dl_6[k]
                  + f_2 * dl_15[k]
                  - f_5 * dl_28[k]
                  - f_5 * dl_136[k]
                  + f_2 * dl_141[k]
                  - f_2 * dl_150[k]
                  + f_5 * dl_163[k];
    }

#pragma omp simd aligned(dl_4, dl_11, dl_22, dl_37, dl_139, dl_146, dl_157, \
                         dl_172 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_142 * dl_4[k]
                  - f_58 * dl_11[k]
                  + f_143 * dl_22[k]
                  - f_144 * dl_37[k]
                  - f_142 * dl_139[k]
                  + f_58 * dl_146[k]
                  - f_143 * dl_157[k]
                  + f_144 * dl_172[k];
    }

#pragma omp simd aligned(dl_1, dl_6, dl_8, dl_15, dl_17, dl_28, dl_30, dl_136, dl_141, dl_143, \
                         dl_150, dl_152, dl_163, dl_165 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_145 * dl_1[k]
                  + f_146 * dl_6[k]
                  + f_147 * dl_8[k]
                  + f_146 * dl_15[k]
                  - f_148 * dl_17[k]
                  - f_145 * dl_28[k]
                  + f_147 * dl_30[k]
                  + f_145 * dl_136[k]
                  - f_146 * dl_141[k]
                  - f_147 * dl_143[k]
                  - f_146 * dl_150[k]
                  + f_148 * dl_152[k]
                  + f_145 * dl_163[k]
                  - f_147 * dl_165[k];
    }

#pragma omp simd aligned(dl_4, dl_11, dl_13, dl_22, dl_24, dl_37, dl_39, dl_139, dl_146, \
                         dl_148, dl_157, dl_159, dl_172, dl_174 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_149 * dl_4[k]
                  + f_149 * dl_11[k]
                  + f_150 * dl_13[k]
                  + f_151 * dl_22[k]
                  - f_11 * dl_24[k]
                  - f_152 * dl_37[k]
                  + f_153 * dl_39[k]
                  + f_149 * dl_139[k]
                  - f_149 * dl_146[k]
                  - f_150 * dl_148[k]
                  - f_151 * dl_157[k]
                  + f_11 * dl_159[k]
                  + f_152 * dl_172[k]
                  - f_153 * dl_174[k];
    }

#pragma omp simd aligned(dl_1, dl_6, dl_8, dl_15, dl_19, dl_28, dl_30, dl_32, dl_136, dl_141, \
                         dl_143, dl_150, dl_154, dl_163, dl_165, \
                         dl_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = f_154 * dl_1[k]
                  + f_154 * dl_6[k]
                  - f_155 * dl_8[k]
                  - f_154 * dl_15[k]
                  + f_156 * dl_19[k]
                  - f_154 * dl_28[k]
                  + f_155 * dl_30[k]
                  - f_156 * dl_32[k]
                  - f_154 * dl_136[k]
                  - f_154 * dl_141[k]
                  + f_155 * dl_143[k]
                  + f_154 * dl_150[k]
                  - f_156 * dl_154[k]
                  + f_154 * dl_163[k]
                  - f_155 * dl_165[k]
                  + f_156 * dl_167[k];
    }

#pragma omp simd aligned(dl_4, dl_11, dl_13, dl_22, dl_24, dl_26, dl_37, dl_39, dl_41, dl_139, \
                         dl_146, dl_148, dl_157, dl_159, dl_161, dl_172, dl_174, \
                         dl_176 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_157 * dl_4[k]
                  + f_158 * dl_11[k]
                  - f_159 * dl_13[k]
                  + f_160 * dl_22[k]
                  - f_25 * dl_24[k]
                  + f_161 * dl_26[k]
                  - f_160 * dl_37[k]
                  + f_162 * dl_39[k]
                  - f_163 * dl_41[k]
                  - f_157 * dl_139[k]
                  - f_158 * dl_146[k]
                  + f_159 * dl_148[k]
                  - f_160 * dl_157[k]
                  + f_25 * dl_159[k]
                  - f_161 * dl_161[k]
                  + f_160 * dl_172[k]
                  - f_162 * dl_174[k]
                  + f_163 * dl_176[k];
    }

#pragma omp simd aligned(dl_1, dl_6, dl_8, dl_15, dl_17, dl_19, dl_28, dl_30, dl_32, dl_34, \
                         dl_136, dl_141, dl_143, dl_150, dl_152, dl_154, dl_163, dl_165, \
                         dl_167, dl_169 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_45 * dl_1[k]
                  - f_164 * dl_6[k]
                  + f_46 * dl_8[k]
                  - f_164 * dl_15[k]
                  + f_29 * dl_17[k]
                  - f_47 * dl_19[k]
                  - f_45 * dl_28[k]
                  + f_46 * dl_30[k]
                  - f_47 * dl_32[k]
                  + f_48 * dl_34[k]
                  + f_45 * dl_136[k]
                  + f_164 * dl_141[k]
                  - f_46 * dl_143[k]
                  + f_164 * dl_150[k]
                  - f_29 * dl_152[k]
                  + f_47 * dl_154[k]
                  + f_45 * dl_163[k]
                  - f_46 * dl_165[k]
                  + f_47 * dl_167[k]
                  - f_48 * dl_169[k];
    }

#pragma omp simd aligned(dl_4, dl_11, dl_13, dl_22, dl_24, dl_26, dl_37, dl_39, dl_41, dl_43, \
                         dl_139, dl_146, dl_148, dl_157, dl_159, dl_161, dl_172, dl_174, \
                         dl_176, dl_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_42 * dl_4[k]
                  - f_165 * dl_11[k]
                  + f_166 * dl_13[k]
                  - f_165 * dl_22[k]
                  + f_35 * dl_24[k]
                  - f_167 * dl_26[k]
                  - f_42 * dl_37[k]
                  + f_166 * dl_39[k]
                  - f_167 * dl_41[k]
                  + f_168 * dl_43[k]
                  + f_42 * dl_139[k]
                  + f_165 * dl_146[k]
                  - f_166 * dl_148[k]
                  + f_165 * dl_157[k]
                  - f_35 * dl_159[k]
                  + f_167 * dl_161[k]
                  + f_42 * dl_172[k]
                  - f_166 * dl_174[k]
                  + f_167 * dl_176[k]
                  - f_168 * dl_178[k];
    }

#pragma omp simd aligned(dl_0, dl_3, dl_5, dl_10, dl_12, dl_14, dl_21, dl_23, dl_25, dl_27, \
                         dl_36, dl_38, dl_40, dl_42, dl_44, dl_135, dl_138, dl_140, dl_145, \
                         dl_147, dl_149, dl_156, dl_158, dl_160, dl_162, dl_171, dl_173, \
                         dl_175, dl_177, dl_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = f_169 * dl_0[k]
                  + f_170 * dl_3[k]
                  - f_171 * dl_5[k]
                  + f_172 * dl_10[k]
                  - f_166 * dl_12[k]
                  + f_166 * dl_14[k]
                  + f_170 * dl_21[k]
                  - f_166 * dl_23[k]
                  + f_35 * dl_25[k]
                  - f_173 * dl_27[k]
                  + f_169 * dl_36[k]
                  - f_171 * dl_38[k]
                  + f_166 * dl_40[k]
                  - f_173 * dl_42[k]
                  + f_174 * dl_44[k]
                  - f_169 * dl_135[k]
                  - f_170 * dl_138[k]
                  + f_171 * dl_140[k]
                  - f_172 * dl_145[k]
                  + f_166 * dl_147[k]
                  - f_166 * dl_149[k]
                  - f_170 * dl_156[k]
                  + f_166 * dl_158[k]
                  - f_35 * dl_160[k]
                  + f_173 * dl_162[k]
                  - f_169 * dl_171[k]
                  + f_171 * dl_173[k]
                  - f_166 * dl_175[k]
                  + f_173 * dl_177[k]
                  - f_174 * dl_179[k];
    }

#pragma omp simd aligned(dl_2, dl_7, dl_9, dl_16, dl_18, dl_20, dl_29, dl_31, dl_33, dl_35, \
                         dl_137, dl_142, dl_144, dl_151, dl_153, dl_155, dl_164, dl_166, \
                         dl_168, dl_170 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = -f_42 * dl_2[k]
                  - f_165 * dl_7[k]
                  + f_166 * dl_9[k]
                  - f_165 * dl_16[k]
                  + f_35 * dl_18[k]
                  - f_167 * dl_20[k]
                  - f_42 * dl_29[k]
                  + f_166 * dl_31[k]
                  - f_167 * dl_33[k]
                  + f_168 * dl_35[k]
                  + f_42 * dl_137[k]
                  + f_165 * dl_142[k]
                  - f_166 * dl_144[k]
                  + f_165 * dl_151[k]
                  - f_35 * dl_153[k]
                  + f_167 * dl_155[k]
                  + f_42 * dl_164[k]
                  - f_166 * dl_166[k]
                  + f_167 * dl_168[k]
                  - f_168 * dl_170[k];
    }

#pragma omp simd aligned(dl_0, dl_3, dl_5, dl_12, dl_14, dl_21, dl_23, dl_27, dl_36, dl_38, \
                         dl_40, dl_42, dl_135, dl_138, dl_140, dl_147, dl_149, dl_156, dl_158, \
                         dl_162, dl_171, dl_173, dl_175, dl_177 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = -f_175 * dl_0[k]
                  - f_45 * dl_3[k]
                  + f_176 * dl_5[k]
                  + f_176 * dl_12[k]
                  - f_177 * dl_14[k]
                  + f_45 * dl_21[k]
                  - f_176 * dl_23[k]
                  + f_178 * dl_27[k]
                  + f_175 * dl_36[k]
                  - f_176 * dl_38[k]
                  + f_177 * dl_40[k]
                  - f_178 * dl_42[k]
                  + f_175 * dl_135[k]
                  + f_45 * dl_138[k]
                  - f_176 * dl_140[k]
                  - f_176 * dl_147[k]
                  + f_177 * dl_149[k]
                  - f_45 * dl_156[k]
                  + f_176 * dl_158[k]
                  - f_178 * dl_162[k]
                  - f_175 * dl_171[k]
                  + f_176 * dl_173[k]
                  - f_177 * dl_175[k]
                  + f_178 * dl_177[k];
    }

#pragma omp simd aligned(dl_2, dl_7, dl_9, dl_16, dl_18, dl_20, dl_29, dl_31, dl_33, dl_137, \
                         dl_142, dl_144, dl_151, dl_153, dl_155, dl_164, dl_166, \
                         dl_168 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = f_160 * dl_2[k]
                  - f_160 * dl_7[k]
                  - f_162 * dl_9[k]
                  - f_158 * dl_16[k]
                  + f_25 * dl_18[k]
                  + f_163 * dl_20[k]
                  - f_157 * dl_29[k]
                  + f_159 * dl_31[k]
                  - f_161 * dl_33[k]
                  - f_160 * dl_137[k]
                  + f_160 * dl_142[k]
                  + f_162 * dl_144[k]
                  + f_158 * dl_151[k]
                  - f_25 * dl_153[k]
                  - f_163 * dl_155[k]
                  + f_157 * dl_164[k]
                  - f_159 * dl_166[k]
                  + f_161 * dl_168[k];
    }

#pragma omp simd aligned(dl_0, dl_3, dl_5, dl_10, dl_12, dl_14, dl_21, dl_23, dl_25, dl_36, \
                         dl_38, dl_40, dl_135, dl_138, dl_140, dl_145, dl_147, dl_149, dl_156, \
                         dl_158, dl_160, dl_171, dl_173, dl_175 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = f_179 * dl_0[k]
                  - f_154 * dl_3[k]
                  - f_180 * dl_5[k]
                  - f_181 * dl_10[k]
                  + f_182 * dl_12[k]
                  + f_183 * dl_14[k]
                  - f_154 * dl_21[k]
                  + f_182 * dl_23[k]
                  - f_52 * dl_25[k]
                  + f_179 * dl_36[k]
                  - f_180 * dl_38[k]
                  + f_183 * dl_40[k]
                  - f_179 * dl_135[k]
                  + f_154 * dl_138[k]
                  + f_180 * dl_140[k]
                  + f_181 * dl_145[k]
                  - f_182 * dl_147[k]
                  - f_183 * dl_149[k]
                  + f_154 * dl_156[k]
                  - f_182 * dl_158[k]
                  + f_52 * dl_160[k]
                  - f_179 * dl_171[k]
                  + f_180 * dl_173[k]
                  - f_183 * dl_175[k];
    }

#pragma omp simd aligned(dl_2, dl_7, dl_9, dl_16, dl_18, dl_29, dl_31, dl_137, dl_142, dl_144, \
                         dl_151, dl_153, dl_164, dl_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = -f_152 * dl_2[k]
                  + f_151 * dl_7[k]
                  + f_153 * dl_9[k]
                  + f_149 * dl_16[k]
                  - f_11 * dl_18[k]
                  - f_149 * dl_29[k]
                  + f_150 * dl_31[k]
                  + f_152 * dl_137[k]
                  - f_151 * dl_142[k]
                  - f_153 * dl_144[k]
                  - f_149 * dl_151[k]
                  + f_11 * dl_153[k]
                  + f_149 * dl_164[k]
                  - f_150 * dl_166[k];
    }

#pragma omp simd aligned(dl_0, dl_3, dl_5, dl_12, dl_21, dl_23, dl_36, dl_38, dl_135, dl_138, \
                         dl_140, dl_147, dl_156, dl_158, dl_171, \
                         dl_173 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = -f_184 * dl_0[k]
                  + f_146 * dl_3[k]
                  + f_146 * dl_5[k]
                  - f_185 * dl_12[k]
                  - f_146 * dl_21[k]
                  + f_185 * dl_23[k]
                  + f_184 * dl_36[k]
                  - f_146 * dl_38[k]
                  + f_184 * dl_135[k]
                  - f_146 * dl_138[k]
                  - f_146 * dl_140[k]
                  + f_185 * dl_147[k]
                  + f_146 * dl_156[k]
                  - f_185 * dl_158[k]
                  - f_184 * dl_171[k]
                  + f_146 * dl_173[k];
    }

#pragma omp simd aligned(dl_2, dl_7, dl_16, dl_29, dl_137, dl_142, dl_151, \
                         dl_164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = f_144 * dl_2[k]
                  - f_143 * dl_7[k]
                  + f_58 * dl_16[k]
                  - f_142 * dl_29[k]
                  - f_144 * dl_137[k]
                  + f_143 * dl_142[k]
                  - f_58 * dl_151[k]
                  + f_142 * dl_164[k];
    }

#pragma omp simd aligned(dl_0, dl_3, dl_10, dl_21, dl_36, dl_135, dl_138, dl_145, dl_156, \
                         dl_171 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = f_186 * dl_0[k]
                  - f_142 * dl_3[k]
                  + f_187 * dl_10[k]
                  - f_142 * dl_21[k]
                  + f_186 * dl_36[k]
                  - f_186 * dl_135[k]
                  + f_142 * dl_138[k]
                  - f_187 * dl_145[k]
                  + f_142 * dl_156[k]
                  - f_186 * dl_171[k];
    }
}

}  // namespace simdtrf
