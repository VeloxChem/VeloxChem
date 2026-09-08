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


#include "SimdTransformDK.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_dk(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t dk,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.65625 * std::sqrt(143.0);
    const auto f_1 = 3.28125 * std::sqrt(143.0);
    const auto f_2 = 1.96875 * std::sqrt(143.0);
    const auto f_3 = 0.09375 * std::sqrt(143.0);
    const auto f_4 = 0.5625 * std::sqrt(2002.0);
    const auto f_5 = 1.875 * std::sqrt(2002.0);
    const auto f_6 = 0.46875 * std::sqrt(77.0);
    const auto f_7 = 5.625 * std::sqrt(77.0);
    const auto f_8 = 0.84375 * std::sqrt(77.0);
    const auto f_9 = 11.25 * std::sqrt(77.0);
    const auto f_10 = 0.09375 * std::sqrt(77.0);
    const auto f_11 = 1.125 * std::sqrt(77.0);
    const auto f_12 = 2.25 * std::sqrt(77.0);
    const auto f_13 = 7.5 * std::sqrt(77.0);
    const auto f_14 = 0.84375 * std::sqrt(7.0);
    const auto f_15 = 1.40625 * std::sqrt(7.0);
    const auto f_16 = 16.875 * std::sqrt(7.0);
    const auto f_17 = 0.28125 * std::sqrt(7.0);
    const auto f_18 = 11.25 * std::sqrt(7.0);
    const auto f_19 = 22.5 * std::sqrt(7.0);
    const auto f_20 = 5.625 * std::sqrt(7.0);
    const auto f_21 = 7.5 * std::sqrt(7.0);
    const auto f_22 = 2.8125 * std::sqrt(14.0);
    const auto f_23 = 5.625 * std::sqrt(14.0);
    const auto f_24 = 15.0 * std::sqrt(14.0);
    const auto f_25 = 9.0 * std::sqrt(14.0);
    const auto f_26 = 0.15625 * std::sqrt(21.0);
    const auto f_27 = 0.46875 * std::sqrt(21.0);
    const auto f_28 = 3.75 * std::sqrt(21.0);
    const auto f_29 = 7.5 * std::sqrt(21.0);
    const auto f_30 = 2.0 * std::sqrt(21.0);
    const auto f_31 = 2.1875 * std::sqrt(3.0);
    const auto f_32 = 6.5625 * std::sqrt(3.0);
    const auto f_33 = 13.125 * std::sqrt(3.0);
    const auto f_34 = 26.25 * std::sqrt(3.0);
    const auto f_35 = 10.5 * std::sqrt(3.0);
    const auto f_36 = std::sqrt(3.0);
    const auto f_37 = 1.40625 * std::sqrt(14.0);
    const auto f_38 = 7.5 * std::sqrt(14.0);
    const auto f_39 = 4.5 * std::sqrt(14.0);
    const auto f_40 = 0.5625 * std::sqrt(77.0);
    const auto f_41 = 2.8125 * std::sqrt(77.0);
    const auto f_42 = 1.875 * std::sqrt(77.0);
    const auto f_43 = 0.09375 * std::sqrt(2002.0);
    const auto f_44 = 1.40625 * std::sqrt(2002.0);
    const auto f_45 = 0.109375 * std::sqrt(429.0);
    const auto f_46 = 0.546875 * std::sqrt(429.0);
    const auto f_47 = 0.328125 * std::sqrt(429.0);
    const auto f_48 = 0.015625 * std::sqrt(429.0);
    const auto f_49 = 0.21875 * std::sqrt(429.0);
    const auto f_50 = 1.09375 * std::sqrt(429.0);
    const auto f_51 = 0.65625 * std::sqrt(429.0);
    const auto f_52 = 0.03125 * std::sqrt(429.0);
    const auto f_53 = 0.09375 * std::sqrt(6006.0);
    const auto f_54 = 0.3125 * std::sqrt(6006.0);
    const auto f_55 = 0.1875 * std::sqrt(6006.0);
    const auto f_56 = 0.625 * std::sqrt(6006.0);
    const auto f_57 = 0.078125 * std::sqrt(231.0);
    const auto f_58 = 0.9375 * std::sqrt(231.0);
    const auto f_59 = 0.140625 * std::sqrt(231.0);
    const auto f_60 = 1.875 * std::sqrt(231.0);
    const auto f_61 = 0.015625 * std::sqrt(231.0);
    const auto f_62 = 0.1875 * std::sqrt(231.0);
    const auto f_63 = 0.15625 * std::sqrt(231.0);
    const auto f_64 = 0.28125 * std::sqrt(231.0);
    const auto f_65 = 3.75 * std::sqrt(231.0);
    const auto f_66 = 0.03125 * std::sqrt(231.0);
    const auto f_67 = 0.375 * std::sqrt(231.0);
    const auto f_68 = 1.25 * std::sqrt(231.0);
    const auto f_69 = 0.75 * std::sqrt(231.0);
    const auto f_70 = 2.5 * std::sqrt(231.0);
    const auto f_71 = 0.140625 * std::sqrt(21.0);
    const auto f_72 = 0.234375 * std::sqrt(21.0);
    const auto f_73 = 2.8125 * std::sqrt(21.0);
    const auto f_74 = 0.046875 * std::sqrt(21.0);
    const auto f_75 = 1.875 * std::sqrt(21.0);
    const auto f_76 = 0.9375 * std::sqrt(21.0);
    const auto f_77 = 1.25 * std::sqrt(21.0);
    const auto f_78 = 0.28125 * std::sqrt(21.0);
    const auto f_79 = 5.625 * std::sqrt(21.0);
    const auto f_80 = 0.09375 * std::sqrt(21.0);
    const auto f_81 = 2.5 * std::sqrt(21.0);
    const auto f_82 = 0.46875 * std::sqrt(42.0);
    const auto f_83 = 0.9375 * std::sqrt(42.0);
    const auto f_84 = 2.5 * std::sqrt(42.0);
    const auto f_85 = 1.5 * std::sqrt(42.0);
    const auto f_86 = 1.875 * std::sqrt(42.0);
    const auto f_87 = 5.0 * std::sqrt(42.0);
    const auto f_88 = 3.0 * std::sqrt(42.0);
    const auto f_89 = 0.078125 * std::sqrt(7.0);
    const auto f_90 = 0.234375 * std::sqrt(7.0);
    const auto f_91 = 1.875 * std::sqrt(7.0);
    const auto f_92 = 3.75 * std::sqrt(7.0);
    const auto f_93 = std::sqrt(7.0);
    const auto f_94 = 0.15625 * std::sqrt(7.0);
    const auto f_95 = 0.46875 * std::sqrt(7.0);
    const auto f_96 = 2.0 * std::sqrt(7.0);
    const auto f_97 = 0.234375 * std::sqrt(42.0);
    const auto f_98 = 1.25 * std::sqrt(42.0);
    const auto f_99 = 0.75 * std::sqrt(42.0);
    const auto f_100 = 0.09375 * std::sqrt(231.0);
    const auto f_101 = 0.46875 * std::sqrt(231.0);
    const auto f_102 = 0.3125 * std::sqrt(231.0);
    const auto f_103 = 0.625 * std::sqrt(231.0);
    const auto f_104 = 0.015625 * std::sqrt(6006.0);
    const auto f_105 = 0.234375 * std::sqrt(6006.0);
    const auto f_106 = 0.03125 * std::sqrt(6006.0);
    const auto f_107 = 0.46875 * std::sqrt(6006.0);
    const auto f_108 = 0.328125 * std::sqrt(143.0);
    const auto f_109 = 1.640625 * std::sqrt(143.0);
    const auto f_110 = 0.984375 * std::sqrt(143.0);
    const auto f_111 = 0.046875 * std::sqrt(143.0);
    const auto f_112 = 0.28125 * std::sqrt(2002.0);
    const auto f_113 = 0.9375 * std::sqrt(2002.0);
    const auto f_114 = 0.234375 * std::sqrt(77.0);
    const auto f_115 = 0.421875 * std::sqrt(77.0);
    const auto f_116 = 0.046875 * std::sqrt(77.0);
    const auto f_117 = 3.75 * std::sqrt(77.0);
    const auto f_118 = 0.421875 * std::sqrt(7.0);
    const auto f_119 = 0.703125 * std::sqrt(7.0);
    const auto f_120 = 8.4375 * std::sqrt(7.0);
    const auto f_121 = 0.140625 * std::sqrt(7.0);
    const auto f_122 = 2.8125 * std::sqrt(7.0);
    const auto f_123 = 0.078125 * std::sqrt(21.0);
    const auto f_124 = std::sqrt(21.0);
    const auto f_125 = 1.09375 * std::sqrt(3.0);
    const auto f_126 = 3.28125 * std::sqrt(3.0);
    const auto f_127 = 5.25 * std::sqrt(3.0);
    const auto f_128 = 0.5 * std::sqrt(3.0);
    const auto f_129 = 0.703125 * std::sqrt(14.0);
    const auto f_130 = 3.75 * std::sqrt(14.0);
    const auto f_131 = 2.25 * std::sqrt(14.0);
    const auto f_132 = 0.28125 * std::sqrt(77.0);
    const auto f_133 = 1.40625 * std::sqrt(77.0);
    const auto f_134 = 0.9375 * std::sqrt(77.0);
    const auto f_135 = 0.046875 * std::sqrt(2002.0);
    const auto f_136 = 0.703125 * std::sqrt(2002.0);

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

    const auto *dk_0 = buffer.data(dk + 0);
    const auto *dk_1 = buffer.data(dk + 1);
    const auto *dk_2 = buffer.data(dk + 2);
    const auto *dk_3 = buffer.data(dk + 3);
    const auto *dk_4 = buffer.data(dk + 4);
    const auto *dk_5 = buffer.data(dk + 5);
    const auto *dk_6 = buffer.data(dk + 6);
    const auto *dk_7 = buffer.data(dk + 7);
    const auto *dk_8 = buffer.data(dk + 8);
    const auto *dk_9 = buffer.data(dk + 9);
    const auto *dk_10 = buffer.data(dk + 10);
    const auto *dk_11 = buffer.data(dk + 11);
    const auto *dk_12 = buffer.data(dk + 12);
    const auto *dk_13 = buffer.data(dk + 13);
    const auto *dk_14 = buffer.data(dk + 14);
    const auto *dk_15 = buffer.data(dk + 15);
    const auto *dk_16 = buffer.data(dk + 16);
    const auto *dk_17 = buffer.data(dk + 17);
    const auto *dk_18 = buffer.data(dk + 18);
    const auto *dk_19 = buffer.data(dk + 19);
    const auto *dk_20 = buffer.data(dk + 20);
    const auto *dk_21 = buffer.data(dk + 21);
    const auto *dk_22 = buffer.data(dk + 22);
    const auto *dk_23 = buffer.data(dk + 23);
    const auto *dk_24 = buffer.data(dk + 24);
    const auto *dk_25 = buffer.data(dk + 25);
    const auto *dk_26 = buffer.data(dk + 26);
    const auto *dk_27 = buffer.data(dk + 27);
    const auto *dk_28 = buffer.data(dk + 28);
    const auto *dk_29 = buffer.data(dk + 29);
    const auto *dk_30 = buffer.data(dk + 30);
    const auto *dk_31 = buffer.data(dk + 31);
    const auto *dk_32 = buffer.data(dk + 32);
    const auto *dk_33 = buffer.data(dk + 33);
    const auto *dk_34 = buffer.data(dk + 34);
    const auto *dk_35 = buffer.data(dk + 35);
    const auto *dk_36 = buffer.data(dk + 36);
    const auto *dk_37 = buffer.data(dk + 37);
    const auto *dk_38 = buffer.data(dk + 38);
    const auto *dk_39 = buffer.data(dk + 39);
    const auto *dk_40 = buffer.data(dk + 40);
    const auto *dk_41 = buffer.data(dk + 41);
    const auto *dk_42 = buffer.data(dk + 42);
    const auto *dk_43 = buffer.data(dk + 43);
    const auto *dk_44 = buffer.data(dk + 44);
    const auto *dk_45 = buffer.data(dk + 45);
    const auto *dk_46 = buffer.data(dk + 46);
    const auto *dk_47 = buffer.data(dk + 47);
    const auto *dk_48 = buffer.data(dk + 48);
    const auto *dk_49 = buffer.data(dk + 49);
    const auto *dk_50 = buffer.data(dk + 50);
    const auto *dk_51 = buffer.data(dk + 51);
    const auto *dk_52 = buffer.data(dk + 52);
    const auto *dk_53 = buffer.data(dk + 53);
    const auto *dk_54 = buffer.data(dk + 54);
    const auto *dk_55 = buffer.data(dk + 55);
    const auto *dk_56 = buffer.data(dk + 56);
    const auto *dk_57 = buffer.data(dk + 57);
    const auto *dk_58 = buffer.data(dk + 58);
    const auto *dk_59 = buffer.data(dk + 59);
    const auto *dk_60 = buffer.data(dk + 60);
    const auto *dk_61 = buffer.data(dk + 61);
    const auto *dk_62 = buffer.data(dk + 62);
    const auto *dk_63 = buffer.data(dk + 63);
    const auto *dk_64 = buffer.data(dk + 64);
    const auto *dk_65 = buffer.data(dk + 65);
    const auto *dk_66 = buffer.data(dk + 66);
    const auto *dk_67 = buffer.data(dk + 67);
    const auto *dk_68 = buffer.data(dk + 68);
    const auto *dk_69 = buffer.data(dk + 69);
    const auto *dk_70 = buffer.data(dk + 70);
    const auto *dk_71 = buffer.data(dk + 71);
    const auto *dk_72 = buffer.data(dk + 72);
    const auto *dk_73 = buffer.data(dk + 73);
    const auto *dk_74 = buffer.data(dk + 74);
    const auto *dk_75 = buffer.data(dk + 75);
    const auto *dk_76 = buffer.data(dk + 76);
    const auto *dk_77 = buffer.data(dk + 77);
    const auto *dk_78 = buffer.data(dk + 78);
    const auto *dk_79 = buffer.data(dk + 79);
    const auto *dk_80 = buffer.data(dk + 80);
    const auto *dk_81 = buffer.data(dk + 81);
    const auto *dk_82 = buffer.data(dk + 82);
    const auto *dk_83 = buffer.data(dk + 83);
    const auto *dk_84 = buffer.data(dk + 84);
    const auto *dk_85 = buffer.data(dk + 85);
    const auto *dk_86 = buffer.data(dk + 86);
    const auto *dk_87 = buffer.data(dk + 87);
    const auto *dk_88 = buffer.data(dk + 88);
    const auto *dk_89 = buffer.data(dk + 89);
    const auto *dk_90 = buffer.data(dk + 90);
    const auto *dk_91 = buffer.data(dk + 91);
    const auto *dk_92 = buffer.data(dk + 92);
    const auto *dk_93 = buffer.data(dk + 93);
    const auto *dk_94 = buffer.data(dk + 94);
    const auto *dk_95 = buffer.data(dk + 95);
    const auto *dk_96 = buffer.data(dk + 96);
    const auto *dk_97 = buffer.data(dk + 97);
    const auto *dk_98 = buffer.data(dk + 98);
    const auto *dk_99 = buffer.data(dk + 99);
    const auto *dk_100 = buffer.data(dk + 100);
    const auto *dk_101 = buffer.data(dk + 101);
    const auto *dk_102 = buffer.data(dk + 102);
    const auto *dk_103 = buffer.data(dk + 103);
    const auto *dk_104 = buffer.data(dk + 104);
    const auto *dk_105 = buffer.data(dk + 105);
    const auto *dk_106 = buffer.data(dk + 106);
    const auto *dk_107 = buffer.data(dk + 107);
    const auto *dk_108 = buffer.data(dk + 108);
    const auto *dk_109 = buffer.data(dk + 109);
    const auto *dk_110 = buffer.data(dk + 110);
    const auto *dk_111 = buffer.data(dk + 111);
    const auto *dk_112 = buffer.data(dk + 112);
    const auto *dk_113 = buffer.data(dk + 113);
    const auto *dk_114 = buffer.data(dk + 114);
    const auto *dk_115 = buffer.data(dk + 115);
    const auto *dk_116 = buffer.data(dk + 116);
    const auto *dk_117 = buffer.data(dk + 117);
    const auto *dk_118 = buffer.data(dk + 118);
    const auto *dk_119 = buffer.data(dk + 119);
    const auto *dk_120 = buffer.data(dk + 120);
    const auto *dk_121 = buffer.data(dk + 121);
    const auto *dk_122 = buffer.data(dk + 122);
    const auto *dk_123 = buffer.data(dk + 123);
    const auto *dk_124 = buffer.data(dk + 124);
    const auto *dk_125 = buffer.data(dk + 125);
    const auto *dk_126 = buffer.data(dk + 126);
    const auto *dk_127 = buffer.data(dk + 127);
    const auto *dk_128 = buffer.data(dk + 128);
    const auto *dk_129 = buffer.data(dk + 129);
    const auto *dk_130 = buffer.data(dk + 130);
    const auto *dk_131 = buffer.data(dk + 131);
    const auto *dk_132 = buffer.data(dk + 132);
    const auto *dk_133 = buffer.data(dk + 133);
    const auto *dk_134 = buffer.data(dk + 134);
    const auto *dk_135 = buffer.data(dk + 135);
    const auto *dk_136 = buffer.data(dk + 136);
    const auto *dk_137 = buffer.data(dk + 137);
    const auto *dk_138 = buffer.data(dk + 138);
    const auto *dk_139 = buffer.data(dk + 139);
    const auto *dk_140 = buffer.data(dk + 140);
    const auto *dk_141 = buffer.data(dk + 141);
    const auto *dk_142 = buffer.data(dk + 142);
    const auto *dk_143 = buffer.data(dk + 143);
    const auto *dk_144 = buffer.data(dk + 144);
    const auto *dk_145 = buffer.data(dk + 145);
    const auto *dk_146 = buffer.data(dk + 146);
    const auto *dk_147 = buffer.data(dk + 147);
    const auto *dk_148 = buffer.data(dk + 148);
    const auto *dk_149 = buffer.data(dk + 149);
    const auto *dk_150 = buffer.data(dk + 150);
    const auto *dk_151 = buffer.data(dk + 151);
    const auto *dk_152 = buffer.data(dk + 152);
    const auto *dk_153 = buffer.data(dk + 153);
    const auto *dk_154 = buffer.data(dk + 154);
    const auto *dk_155 = buffer.data(dk + 155);
    const auto *dk_156 = buffer.data(dk + 156);
    const auto *dk_157 = buffer.data(dk + 157);
    const auto *dk_158 = buffer.data(dk + 158);
    const auto *dk_159 = buffer.data(dk + 159);
    const auto *dk_160 = buffer.data(dk + 160);
    const auto *dk_161 = buffer.data(dk + 161);
    const auto *dk_162 = buffer.data(dk + 162);
    const auto *dk_163 = buffer.data(dk + 163);
    const auto *dk_164 = buffer.data(dk + 164);
    const auto *dk_165 = buffer.data(dk + 165);
    const auto *dk_166 = buffer.data(dk + 166);
    const auto *dk_167 = buffer.data(dk + 167);
    const auto *dk_168 = buffer.data(dk + 168);
    const auto *dk_169 = buffer.data(dk + 169);
    const auto *dk_170 = buffer.data(dk + 170);
    const auto *dk_171 = buffer.data(dk + 171);
    const auto *dk_172 = buffer.data(dk + 172);
    const auto *dk_173 = buffer.data(dk + 173);
    const auto *dk_174 = buffer.data(dk + 174);
    const auto *dk_175 = buffer.data(dk + 175);
    const auto *dk_176 = buffer.data(dk + 176);
    const auto *dk_177 = buffer.data(dk + 177);
    const auto *dk_178 = buffer.data(dk + 178);
    const auto *dk_179 = buffer.data(dk + 179);
    const auto *dk_180 = buffer.data(dk + 180);
    const auto *dk_181 = buffer.data(dk + 181);
    const auto *dk_182 = buffer.data(dk + 182);
    const auto *dk_183 = buffer.data(dk + 183);
    const auto *dk_184 = buffer.data(dk + 184);
    const auto *dk_185 = buffer.data(dk + 185);
    const auto *dk_186 = buffer.data(dk + 186);
    const auto *dk_187 = buffer.data(dk + 187);
    const auto *dk_188 = buffer.data(dk + 188);
    const auto *dk_189 = buffer.data(dk + 189);
    const auto *dk_190 = buffer.data(dk + 190);
    const auto *dk_191 = buffer.data(dk + 191);
    const auto *dk_192 = buffer.data(dk + 192);
    const auto *dk_193 = buffer.data(dk + 193);
    const auto *dk_194 = buffer.data(dk + 194);
    const auto *dk_195 = buffer.data(dk + 195);
    const auto *dk_196 = buffer.data(dk + 196);
    const auto *dk_197 = buffer.data(dk + 197);
    const auto *dk_198 = buffer.data(dk + 198);
    const auto *dk_199 = buffer.data(dk + 199);
    const auto *dk_200 = buffer.data(dk + 200);
    const auto *dk_201 = buffer.data(dk + 201);
    const auto *dk_202 = buffer.data(dk + 202);
    const auto *dk_203 = buffer.data(dk + 203);
    const auto *dk_204 = buffer.data(dk + 204);
    const auto *dk_205 = buffer.data(dk + 205);
    const auto *dk_206 = buffer.data(dk + 206);
    const auto *dk_207 = buffer.data(dk + 207);
    const auto *dk_208 = buffer.data(dk + 208);
    const auto *dk_209 = buffer.data(dk + 209);
    const auto *dk_210 = buffer.data(dk + 210);
    const auto *dk_211 = buffer.data(dk + 211);
    const auto *dk_212 = buffer.data(dk + 212);
    const auto *dk_213 = buffer.data(dk + 213);
    const auto *dk_214 = buffer.data(dk + 214);
    const auto *dk_215 = buffer.data(dk + 215);

#pragma omp simd aligned(dk_37, dk_40, dk_42, dk_44, dk_47, dk_49, dk_51, dk_53, dk_58, dk_60, \
                         dk_64, dk_66 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * dk_37[k]
                 - f_1 * dk_42[k]
                 + f_2 * dk_51[k]
                 - f_3 * dk_64[k];

        g_1[k] = f_4 * dk_40[k]
                 - f_5 * dk_47[k]
                 + f_4 * dk_58[k];

        g_2[k] = -f_6 * dk_37[k]
                 + f_6 * dk_42[k]
                 + f_7 * dk_44[k]
                 + f_8 * dk_51[k]
                 - f_9 * dk_53[k]
                 - f_10 * dk_64[k]
                 + f_11 * dk_66[k];

        g_3[k] = -f_12 * dk_40[k]
                 + f_13 * dk_49[k]
                 + f_12 * dk_58[k]
                 - f_13 * dk_60[k];
    }

#pragma omp simd aligned(dk_37, dk_42, dk_44, dk_51, dk_53, dk_55, dk_64, dk_66, \
                         dk_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_14 * dk_37[k]
                 + f_15 * dk_42[k]
                 - f_16 * dk_44[k]
                 + f_17 * dk_51[k]
                 - f_18 * dk_53[k]
                 + f_19 * dk_55[k]
                 - f_17 * dk_64[k]
                 + f_20 * dk_66[k]
                 - f_21 * dk_68[k];
    }

#pragma omp simd aligned(dk_40, dk_47, dk_49, dk_58, dk_60, dk_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_22 * dk_40[k]
                 + f_23 * dk_47[k]
                 - f_24 * dk_49[k]
                 + f_22 * dk_58[k]
                 - f_24 * dk_60[k]
                 + f_25 * dk_62[k];
    }

#pragma omp simd aligned(dk_37, dk_42, dk_44, dk_51, dk_53, dk_55, dk_64, dk_66, dk_68, \
                         dk_70 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_26 * dk_37[k]
                 - f_27 * dk_42[k]
                 + f_28 * dk_44[k]
                 - f_27 * dk_51[k]
                 + f_29 * dk_53[k]
                 - f_29 * dk_55[k]
                 - f_26 * dk_64[k]
                 + f_28 * dk_66[k]
                 - f_29 * dk_68[k]
                 + f_30 * dk_70[k];
    }

#pragma omp simd aligned(dk_38, dk_43, dk_45, dk_52, dk_54, dk_56, dk_65, dk_67, dk_69, \
                         dk_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_31 * dk_38[k]
                 - f_32 * dk_43[k]
                 + f_33 * dk_45[k]
                 - f_32 * dk_52[k]
                 + f_34 * dk_54[k]
                 - f_35 * dk_56[k]
                 - f_31 * dk_65[k]
                 + f_33 * dk_67[k]
                 - f_35 * dk_69[k]
                 + f_36 * dk_71[k];
    }

#pragma omp simd aligned(dk_36, dk_39, dk_41, dk_46, dk_48, dk_50, dk_57, dk_59, dk_61, \
                         dk_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_26 * dk_36[k]
                 - f_27 * dk_39[k]
                 + f_28 * dk_41[k]
                 - f_27 * dk_46[k]
                 + f_29 * dk_48[k]
                 - f_29 * dk_50[k]
                 - f_26 * dk_57[k]
                 + f_28 * dk_59[k]
                 - f_29 * dk_61[k]
                 + f_30 * dk_63[k];
    }

#pragma omp simd aligned(dk_38, dk_43, dk_45, dk_52, dk_56, dk_65, dk_67, \
                         dk_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_37 * dk_38[k]
                 + f_37 * dk_43[k]
                 - f_38 * dk_45[k]
                 - f_37 * dk_52[k]
                 + f_39 * dk_56[k]
                 - f_37 * dk_65[k]
                 + f_38 * dk_67[k]
                 - f_39 * dk_69[k];
    }

#pragma omp simd aligned(dk_36, dk_39, dk_41, dk_46, dk_48, dk_50, dk_57, dk_59, \
                         dk_61 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = f_17 * dk_36[k]
                  - f_17 * dk_39[k]
                  - f_20 * dk_41[k]
                  - f_15 * dk_46[k]
                  + f_18 * dk_48[k]
                  + f_21 * dk_50[k]
                  - f_14 * dk_57[k]
                  + f_16 * dk_59[k]
                  - f_19 * dk_61[k];
    }

#pragma omp simd aligned(dk_36, dk_38, dk_39, dk_41, dk_43, dk_45, dk_46, dk_48, dk_52, dk_54, \
                         dk_57, dk_59, dk_65, dk_67 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_40 * dk_38[k]
                  + f_41 * dk_43[k]
                  + f_42 * dk_45[k]
                  + f_41 * dk_52[k]
                  - f_9 * dk_54[k]
                  - f_40 * dk_65[k]
                  + f_42 * dk_67[k];

        g_12[k] = -f_10 * dk_36[k]
                  + f_8 * dk_39[k]
                  + f_11 * dk_41[k]
                  + f_6 * dk_46[k]
                  - f_9 * dk_48[k]
                  - f_6 * dk_57[k]
                  + f_7 * dk_59[k];
    }

#pragma omp simd aligned(dk_36, dk_38, dk_39, dk_43, dk_46, dk_52, dk_57, dk_65, dk_145, \
                         dk_150, dk_159, dk_172 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_43 * dk_38[k]
                  - f_44 * dk_43[k]
                  + f_44 * dk_52[k]
                  - f_43 * dk_65[k];

        g_14[k] = f_3 * dk_36[k]
                  - f_2 * dk_39[k]
                  + f_1 * dk_46[k]
                  - f_0 * dk_57[k];

        g_15[k] = f_0 * dk_145[k]
                  - f_1 * dk_150[k]
                  + f_2 * dk_159[k]
                  - f_3 * dk_172[k];
    }

#pragma omp simd aligned(dk_145, dk_148, dk_150, dk_152, dk_155, dk_157, dk_159, dk_161, \
                         dk_166, dk_168, dk_172, dk_174 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_4 * dk_148[k]
                  - f_5 * dk_155[k]
                  + f_4 * dk_166[k];

        g_17[k] = -f_6 * dk_145[k]
                  + f_6 * dk_150[k]
                  + f_7 * dk_152[k]
                  + f_8 * dk_159[k]
                  - f_9 * dk_161[k]
                  - f_10 * dk_172[k]
                  + f_11 * dk_174[k];

        g_18[k] = -f_12 * dk_148[k]
                  + f_13 * dk_157[k]
                  + f_12 * dk_166[k]
                  - f_13 * dk_168[k];
    }

#pragma omp simd aligned(dk_145, dk_150, dk_152, dk_159, dk_161, dk_163, dk_172, dk_174, \
                         dk_176 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_14 * dk_145[k]
                  + f_15 * dk_150[k]
                  - f_16 * dk_152[k]
                  + f_17 * dk_159[k]
                  - f_18 * dk_161[k]
                  + f_19 * dk_163[k]
                  - f_17 * dk_172[k]
                  + f_20 * dk_174[k]
                  - f_21 * dk_176[k];
    }

#pragma omp simd aligned(dk_148, dk_155, dk_157, dk_166, dk_168, \
                         dk_170 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_22 * dk_148[k]
                  + f_23 * dk_155[k]
                  - f_24 * dk_157[k]
                  + f_22 * dk_166[k]
                  - f_24 * dk_168[k]
                  + f_25 * dk_170[k];
    }

#pragma omp simd aligned(dk_145, dk_150, dk_152, dk_159, dk_161, dk_163, dk_172, dk_174, \
                         dk_176, dk_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_26 * dk_145[k]
                  - f_27 * dk_150[k]
                  + f_28 * dk_152[k]
                  - f_27 * dk_159[k]
                  + f_29 * dk_161[k]
                  - f_29 * dk_163[k]
                  - f_26 * dk_172[k]
                  + f_28 * dk_174[k]
                  - f_29 * dk_176[k]
                  + f_30 * dk_178[k];
    }

#pragma omp simd aligned(dk_146, dk_151, dk_153, dk_160, dk_162, dk_164, dk_173, dk_175, \
                         dk_177, dk_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_31 * dk_146[k]
                  - f_32 * dk_151[k]
                  + f_33 * dk_153[k]
                  - f_32 * dk_160[k]
                  + f_34 * dk_162[k]
                  - f_35 * dk_164[k]
                  - f_31 * dk_173[k]
                  + f_33 * dk_175[k]
                  - f_35 * dk_177[k]
                  + f_36 * dk_179[k];
    }

#pragma omp simd aligned(dk_144, dk_147, dk_149, dk_154, dk_156, dk_158, dk_165, dk_167, \
                         dk_169, dk_171 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_26 * dk_144[k]
                  - f_27 * dk_147[k]
                  + f_28 * dk_149[k]
                  - f_27 * dk_154[k]
                  + f_29 * dk_156[k]
                  - f_29 * dk_158[k]
                  - f_26 * dk_165[k]
                  + f_28 * dk_167[k]
                  - f_29 * dk_169[k]
                  + f_30 * dk_171[k];
    }

#pragma omp simd aligned(dk_146, dk_151, dk_153, dk_160, dk_164, dk_173, dk_175, \
                         dk_177 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_37 * dk_146[k]
                  + f_37 * dk_151[k]
                  - f_38 * dk_153[k]
                  - f_37 * dk_160[k]
                  + f_39 * dk_164[k]
                  - f_37 * dk_173[k]
                  + f_38 * dk_175[k]
                  - f_39 * dk_177[k];
    }

#pragma omp simd aligned(dk_144, dk_147, dk_149, dk_154, dk_156, dk_158, dk_165, dk_167, \
                         dk_169 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_17 * dk_144[k]
                  - f_17 * dk_147[k]
                  - f_20 * dk_149[k]
                  - f_15 * dk_154[k]
                  + f_18 * dk_156[k]
                  + f_21 * dk_158[k]
                  - f_14 * dk_165[k]
                  + f_16 * dk_167[k]
                  - f_19 * dk_169[k];
    }

#pragma omp simd aligned(dk_144, dk_146, dk_147, dk_149, dk_151, dk_153, dk_154, dk_156, \
                         dk_160, dk_162, dk_165, dk_167, dk_173, \
                         dk_175 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_40 * dk_146[k]
                  + f_41 * dk_151[k]
                  + f_42 * dk_153[k]
                  + f_41 * dk_160[k]
                  - f_9 * dk_162[k]
                  - f_40 * dk_173[k]
                  + f_42 * dk_175[k];

        g_27[k] = -f_10 * dk_144[k]
                  + f_8 * dk_147[k]
                  + f_11 * dk_149[k]
                  + f_6 * dk_154[k]
                  - f_9 * dk_156[k]
                  - f_6 * dk_165[k]
                  + f_7 * dk_167[k];
    }

#pragma omp simd aligned(dk_144, dk_146, dk_147, dk_151, dk_154, dk_160, dk_165, \
                         dk_173 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_43 * dk_146[k]
                  - f_44 * dk_151[k]
                  + f_44 * dk_160[k]
                  - f_43 * dk_173[k];

        g_29[k] = f_3 * dk_144[k]
                  - f_2 * dk_147[k]
                  + f_1 * dk_154[k]
                  - f_0 * dk_165[k];
    }

#pragma omp simd aligned(dk_1, dk_6, dk_15, dk_28, dk_109, dk_114, dk_123, dk_136, dk_181, \
                         dk_186, dk_195, dk_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_45 * dk_1[k]
                  + f_46 * dk_6[k]
                  - f_47 * dk_15[k]
                  + f_48 * dk_28[k]
                  - f_45 * dk_109[k]
                  + f_46 * dk_114[k]
                  - f_47 * dk_123[k]
                  + f_48 * dk_136[k]
                  + f_49 * dk_181[k]
                  - f_50 * dk_186[k]
                  + f_51 * dk_195[k]
                  - f_52 * dk_208[k];
    }

#pragma omp simd aligned(dk_4, dk_11, dk_22, dk_112, dk_119, dk_130, dk_184, dk_191, \
                         dk_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_53 * dk_4[k]
                  + f_54 * dk_11[k]
                  - f_53 * dk_22[k]
                  - f_53 * dk_112[k]
                  + f_54 * dk_119[k]
                  - f_53 * dk_130[k]
                  + f_55 * dk_184[k]
                  - f_56 * dk_191[k]
                  + f_55 * dk_202[k];
    }

#pragma omp simd aligned(dk_1, dk_6, dk_8, dk_15, dk_17, dk_28, dk_30, dk_109, dk_114, dk_116, \
                         dk_123, dk_125, dk_136, dk_138, dk_181, dk_186, dk_188, dk_195, \
                         dk_197, dk_208, dk_210 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_57 * dk_1[k]
                  - f_57 * dk_6[k]
                  - f_58 * dk_8[k]
                  - f_59 * dk_15[k]
                  + f_60 * dk_17[k]
                  + f_61 * dk_28[k]
                  - f_62 * dk_30[k]
                  + f_57 * dk_109[k]
                  - f_57 * dk_114[k]
                  - f_58 * dk_116[k]
                  - f_59 * dk_123[k]
                  + f_60 * dk_125[k]
                  + f_61 * dk_136[k]
                  - f_62 * dk_138[k]
                  - f_63 * dk_181[k]
                  + f_63 * dk_186[k]
                  + f_60 * dk_188[k]
                  + f_64 * dk_195[k]
                  - f_65 * dk_197[k]
                  - f_66 * dk_208[k]
                  + f_67 * dk_210[k];
    }

#pragma omp simd aligned(dk_4, dk_13, dk_22, dk_24, dk_112, dk_121, dk_130, dk_132, dk_184, \
                         dk_193, dk_202, dk_204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_67 * dk_4[k]
                  - f_68 * dk_13[k]
                  - f_67 * dk_22[k]
                  + f_68 * dk_24[k]
                  + f_67 * dk_112[k]
                  - f_68 * dk_121[k]
                  - f_67 * dk_130[k]
                  + f_68 * dk_132[k]
                  - f_69 * dk_184[k]
                  + f_70 * dk_193[k]
                  + f_69 * dk_202[k]
                  - f_70 * dk_204[k];
    }

#pragma omp simd aligned(dk_1, dk_6, dk_8, dk_15, dk_17, dk_19, dk_28, dk_30, dk_32, dk_109, \
                         dk_114, dk_116, dk_123, dk_125, dk_127, dk_136, dk_138, dk_140, \
                         dk_181, dk_186, dk_188, dk_195, dk_197, dk_199, dk_208, dk_210, \
                         dk_212 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_71 * dk_1[k]
                  - f_72 * dk_6[k]
                  + f_73 * dk_8[k]
                  - f_74 * dk_15[k]
                  + f_75 * dk_17[k]
                  - f_28 * dk_19[k]
                  + f_74 * dk_28[k]
                  - f_76 * dk_30[k]
                  + f_77 * dk_32[k]
                  - f_71 * dk_109[k]
                  - f_72 * dk_114[k]
                  + f_73 * dk_116[k]
                  - f_74 * dk_123[k]
                  + f_75 * dk_125[k]
                  - f_28 * dk_127[k]
                  + f_74 * dk_136[k]
                  - f_76 * dk_138[k]
                  + f_77 * dk_140[k]
                  + f_78 * dk_181[k]
                  + f_27 * dk_186[k]
                  - f_79 * dk_188[k]
                  + f_80 * dk_195[k]
                  - f_28 * dk_197[k]
                  + f_29 * dk_199[k]
                  - f_80 * dk_208[k]
                  + f_75 * dk_210[k]
                  - f_81 * dk_212[k];
    }

#pragma omp simd aligned(dk_4, dk_11, dk_13, dk_22, dk_24, dk_26, dk_112, dk_119, dk_121, \
                         dk_130, dk_132, dk_134, dk_184, dk_191, dk_193, dk_202, dk_204, \
                         dk_206 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_82 * dk_4[k]
                  - f_83 * dk_11[k]
                  + f_84 * dk_13[k]
                  - f_82 * dk_22[k]
                  + f_84 * dk_24[k]
                  - f_85 * dk_26[k]
                  - f_82 * dk_112[k]
                  - f_83 * dk_119[k]
                  + f_84 * dk_121[k]
                  - f_82 * dk_130[k]
                  + f_84 * dk_132[k]
                  - f_85 * dk_134[k]
                  + f_83 * dk_184[k]
                  + f_86 * dk_191[k]
                  - f_87 * dk_193[k]
                  + f_83 * dk_202[k]
                  - f_87 * dk_204[k]
                  + f_88 * dk_206[k];
    }

#pragma omp simd aligned(dk_1, dk_6, dk_8, dk_15, dk_17, dk_19, dk_28, dk_30, dk_32, dk_34, \
                         dk_109, dk_114, dk_116, dk_123, dk_125, dk_127, dk_136, dk_138, \
                         dk_140, dk_142, dk_181, dk_186, dk_188, dk_195, dk_197, dk_199, \
                         dk_208, dk_210, dk_212, dk_214 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_89 * dk_1[k]
                  + f_90 * dk_6[k]
                  - f_91 * dk_8[k]
                  + f_90 * dk_15[k]
                  - f_92 * dk_17[k]
                  + f_92 * dk_19[k]
                  + f_89 * dk_28[k]
                  - f_91 * dk_30[k]
                  + f_92 * dk_32[k]
                  - f_93 * dk_34[k]
                  + f_89 * dk_109[k]
                  + f_90 * dk_114[k]
                  - f_91 * dk_116[k]
                  + f_90 * dk_123[k]
                  - f_92 * dk_125[k]
                  + f_92 * dk_127[k]
                  + f_89 * dk_136[k]
                  - f_91 * dk_138[k]
                  + f_92 * dk_140[k]
                  - f_93 * dk_142[k]
                  - f_94 * dk_181[k]
                  - f_95 * dk_186[k]
                  + f_92 * dk_188[k]
                  - f_95 * dk_195[k]
                  + f_21 * dk_197[k]
                  - f_21 * dk_199[k]
                  - f_94 * dk_208[k]
                  + f_92 * dk_210[k]
                  - f_21 * dk_212[k]
                  + f_96 * dk_214[k];
    }

#pragma omp simd aligned(dk_2, dk_7, dk_9, dk_16, dk_18, dk_20, dk_29, dk_31, dk_33, dk_35, \
                         dk_110, dk_115, dk_117, dk_124, dk_126, dk_128, dk_137, dk_139, \
                         dk_141, dk_143, dk_182, dk_187, dk_189, dk_196, dk_198, dk_200, \
                         dk_209, dk_211, dk_213, dk_215 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = 1.09375 * dk_2[k]
                  + 3.28125 * dk_7[k]
                  - 6.5625 * dk_9[k]
                  + 3.28125 * dk_16[k]
                  - 13.125 * dk_18[k]
                  + 5.25 * dk_20[k]
                  + 1.09375 * dk_29[k]
                  - 6.5625 * dk_31[k]
                  + 5.25 * dk_33[k]
                  - 0.5 * dk_35[k]
                  + 1.09375 * dk_110[k]
                  + 3.28125 * dk_115[k]
                  - 6.5625 * dk_117[k]
                  + 3.28125 * dk_124[k]
                  - 13.125 * dk_126[k]
                  + 5.25 * dk_128[k]
                  + 1.09375 * dk_137[k]
                  - 6.5625 * dk_139[k]
                  + 5.25 * dk_141[k]
                  - 0.5 * dk_143[k]
                  - 2.1875 * dk_182[k]
                  - 6.5625 * dk_187[k]
                  + 13.125 * dk_189[k]
                  - 6.5625 * dk_196[k]
                  + 26.25 * dk_198[k]
                  - 10.5 * dk_200[k]
                  - 2.1875 * dk_209[k]
                  + 13.125 * dk_211[k]
                  - 10.5 * dk_213[k]
                  + dk_215[k];
    }

#pragma omp simd aligned(dk_0, dk_3, dk_5, dk_10, dk_12, dk_14, dk_21, dk_23, dk_25, dk_27, \
                         dk_108, dk_111, dk_113, dk_118, dk_120, dk_122, dk_129, dk_131, \
                         dk_133, dk_135, dk_180, dk_183, dk_185, dk_190, dk_192, dk_194, \
                         dk_201, dk_203, dk_205, dk_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_89 * dk_0[k]
                  + f_90 * dk_3[k]
                  - f_91 * dk_5[k]
                  + f_90 * dk_10[k]
                  - f_92 * dk_12[k]
                  + f_92 * dk_14[k]
                  + f_89 * dk_21[k]
                  - f_91 * dk_23[k]
                  + f_92 * dk_25[k]
                  - f_93 * dk_27[k]
                  + f_89 * dk_108[k]
                  + f_90 * dk_111[k]
                  - f_91 * dk_113[k]
                  + f_90 * dk_118[k]
                  - f_92 * dk_120[k]
                  + f_92 * dk_122[k]
                  + f_89 * dk_129[k]
                  - f_91 * dk_131[k]
                  + f_92 * dk_133[k]
                  - f_93 * dk_135[k]
                  - f_94 * dk_180[k]
                  - f_95 * dk_183[k]
                  + f_92 * dk_185[k]
                  - f_95 * dk_190[k]
                  + f_21 * dk_192[k]
                  - f_21 * dk_194[k]
                  - f_94 * dk_201[k]
                  + f_92 * dk_203[k]
                  - f_21 * dk_205[k]
                  + f_96 * dk_207[k];
    }

#pragma omp simd aligned(dk_2, dk_7, dk_9, dk_16, dk_20, dk_29, dk_31, dk_33, dk_110, dk_115, \
                         dk_117, dk_124, dk_128, dk_137, dk_139, dk_141, dk_182, dk_187, \
                         dk_189, dk_196, dk_200, dk_209, dk_211, \
                         dk_213 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_97 * dk_2[k]
                  - f_97 * dk_7[k]
                  + f_98 * dk_9[k]
                  + f_97 * dk_16[k]
                  - f_99 * dk_20[k]
                  + f_97 * dk_29[k]
                  - f_98 * dk_31[k]
                  + f_99 * dk_33[k]
                  - f_97 * dk_110[k]
                  - f_97 * dk_115[k]
                  + f_98 * dk_117[k]
                  + f_97 * dk_124[k]
                  - f_99 * dk_128[k]
                  + f_97 * dk_137[k]
                  - f_98 * dk_139[k]
                  + f_99 * dk_141[k]
                  + f_82 * dk_182[k]
                  + f_82 * dk_187[k]
                  - f_84 * dk_189[k]
                  - f_82 * dk_196[k]
                  + f_85 * dk_200[k]
                  - f_82 * dk_209[k]
                  + f_84 * dk_211[k]
                  - f_85 * dk_213[k];
    }

#pragma omp simd aligned(dk_0, dk_3, dk_5, dk_10, dk_12, dk_14, dk_21, dk_23, dk_25, dk_108, \
                         dk_111, dk_113, dk_118, dk_120, dk_122, dk_129, dk_131, dk_133, \
                         dk_180, dk_183, dk_185, dk_190, dk_192, dk_194, dk_201, dk_203, \
                         dk_205 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_74 * dk_0[k]
                  + f_74 * dk_3[k]
                  + f_76 * dk_5[k]
                  + f_72 * dk_10[k]
                  - f_75 * dk_12[k]
                  - f_77 * dk_14[k]
                  + f_71 * dk_21[k]
                  - f_73 * dk_23[k]
                  + f_28 * dk_25[k]
                  - f_74 * dk_108[k]
                  + f_74 * dk_111[k]
                  + f_76 * dk_113[k]
                  + f_72 * dk_118[k]
                  - f_75 * dk_120[k]
                  - f_77 * dk_122[k]
                  + f_71 * dk_129[k]
                  - f_73 * dk_131[k]
                  + f_28 * dk_133[k]
                  + f_80 * dk_180[k]
                  - f_80 * dk_183[k]
                  - f_75 * dk_185[k]
                  - f_27 * dk_190[k]
                  + f_28 * dk_192[k]
                  + f_81 * dk_194[k]
                  - f_78 * dk_201[k]
                  + f_79 * dk_203[k]
                  - f_29 * dk_205[k];
    }

#pragma omp simd aligned(dk_2, dk_7, dk_9, dk_16, dk_18, dk_29, dk_31, dk_110, dk_115, dk_117, \
                         dk_124, dk_126, dk_137, dk_139, dk_182, dk_187, dk_189, dk_196, \
                         dk_198, dk_209, dk_211 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_100 * dk_2[k]
                  - f_101 * dk_7[k]
                  - f_102 * dk_9[k]
                  - f_101 * dk_16[k]
                  + f_60 * dk_18[k]
                  + f_100 * dk_29[k]
                  - f_102 * dk_31[k]
                  + f_100 * dk_110[k]
                  - f_101 * dk_115[k]
                  - f_102 * dk_117[k]
                  - f_101 * dk_124[k]
                  + f_60 * dk_126[k]
                  + f_100 * dk_137[k]
                  - f_102 * dk_139[k]
                  - f_62 * dk_182[k]
                  + f_58 * dk_187[k]
                  + f_103 * dk_189[k]
                  + f_58 * dk_196[k]
                  - f_65 * dk_198[k]
                  - f_62 * dk_209[k]
                  + f_103 * dk_211[k];
    }

#pragma omp simd aligned(dk_0, dk_3, dk_5, dk_10, dk_12, dk_21, dk_23, dk_108, dk_111, dk_113, \
                         dk_118, dk_120, dk_129, dk_131, dk_180, dk_183, dk_185, dk_190, \
                         dk_192, dk_201, dk_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_61 * dk_0[k]
                  - f_59 * dk_3[k]
                  - f_62 * dk_5[k]
                  - f_57 * dk_10[k]
                  + f_60 * dk_12[k]
                  + f_57 * dk_21[k]
                  - f_58 * dk_23[k]
                  + f_61 * dk_108[k]
                  - f_59 * dk_111[k]
                  - f_62 * dk_113[k]
                  - f_57 * dk_118[k]
                  + f_60 * dk_120[k]
                  + f_57 * dk_129[k]
                  - f_58 * dk_131[k]
                  - f_66 * dk_180[k]
                  + f_64 * dk_183[k]
                  + f_67 * dk_185[k]
                  + f_63 * dk_190[k]
                  - f_65 * dk_192[k]
                  - f_63 * dk_201[k]
                  + f_60 * dk_203[k];
    }

#pragma omp simd aligned(dk_2, dk_7, dk_16, dk_29, dk_110, dk_115, dk_124, dk_137, dk_182, \
                         dk_187, dk_196, dk_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_104 * dk_2[k]
                  + f_105 * dk_7[k]
                  - f_105 * dk_16[k]
                  + f_104 * dk_29[k]
                  - f_104 * dk_110[k]
                  + f_105 * dk_115[k]
                  - f_105 * dk_124[k]
                  + f_104 * dk_137[k]
                  + f_106 * dk_182[k]
                  - f_107 * dk_187[k]
                  + f_107 * dk_196[k]
                  - f_106 * dk_209[k];
    }

#pragma omp simd aligned(dk_0, dk_3, dk_10, dk_21, dk_108, dk_111, dk_118, dk_129, dk_180, \
                         dk_183, dk_190, dk_201 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_48 * dk_0[k]
                  + f_47 * dk_3[k]
                  - f_46 * dk_10[k]
                  + f_45 * dk_21[k]
                  - f_48 * dk_108[k]
                  + f_47 * dk_111[k]
                  - f_46 * dk_118[k]
                  + f_45 * dk_129[k]
                  + f_52 * dk_180[k]
                  - f_51 * dk_183[k]
                  + f_50 * dk_190[k]
                  - f_49 * dk_201[k];
    }

#pragma omp simd aligned(dk_73, dk_76, dk_78, dk_80, dk_83, dk_85, dk_87, dk_89, dk_94, dk_96, \
                         dk_100, dk_102 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_0 * dk_73[k]
                  - f_1 * dk_78[k]
                  + f_2 * dk_87[k]
                  - f_3 * dk_100[k];

        g_46[k] = f_4 * dk_76[k]
                  - f_5 * dk_83[k]
                  + f_4 * dk_94[k];

        g_47[k] = -f_6 * dk_73[k]
                  + f_6 * dk_78[k]
                  + f_7 * dk_80[k]
                  + f_8 * dk_87[k]
                  - f_9 * dk_89[k]
                  - f_10 * dk_100[k]
                  + f_11 * dk_102[k];

        g_48[k] = -f_12 * dk_76[k]
                  + f_13 * dk_85[k]
                  + f_12 * dk_94[k]
                  - f_13 * dk_96[k];
    }

#pragma omp simd aligned(dk_73, dk_78, dk_80, dk_87, dk_89, dk_91, dk_100, dk_102, \
                         dk_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_14 * dk_73[k]
                  + f_15 * dk_78[k]
                  - f_16 * dk_80[k]
                  + f_17 * dk_87[k]
                  - f_18 * dk_89[k]
                  + f_19 * dk_91[k]
                  - f_17 * dk_100[k]
                  + f_20 * dk_102[k]
                  - f_21 * dk_104[k];
    }

#pragma omp simd aligned(dk_76, dk_83, dk_85, dk_94, dk_96, dk_98 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = f_22 * dk_76[k]
                  + f_23 * dk_83[k]
                  - f_24 * dk_85[k]
                  + f_22 * dk_94[k]
                  - f_24 * dk_96[k]
                  + f_25 * dk_98[k];
    }

#pragma omp simd aligned(dk_73, dk_78, dk_80, dk_87, dk_89, dk_91, dk_100, dk_102, dk_104, \
                         dk_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_26 * dk_73[k]
                  - f_27 * dk_78[k]
                  + f_28 * dk_80[k]
                  - f_27 * dk_87[k]
                  + f_29 * dk_89[k]
                  - f_29 * dk_91[k]
                  - f_26 * dk_100[k]
                  + f_28 * dk_102[k]
                  - f_29 * dk_104[k]
                  + f_30 * dk_106[k];
    }

#pragma omp simd aligned(dk_74, dk_79, dk_81, dk_88, dk_90, dk_92, dk_101, dk_103, dk_105, \
                         dk_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_31 * dk_74[k]
                  - f_32 * dk_79[k]
                  + f_33 * dk_81[k]
                  - f_32 * dk_88[k]
                  + f_34 * dk_90[k]
                  - f_35 * dk_92[k]
                  - f_31 * dk_101[k]
                  + f_33 * dk_103[k]
                  - f_35 * dk_105[k]
                  + f_36 * dk_107[k];
    }

#pragma omp simd aligned(dk_72, dk_75, dk_77, dk_82, dk_84, dk_86, dk_93, dk_95, dk_97, \
                         dk_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = -f_26 * dk_72[k]
                  - f_27 * dk_75[k]
                  + f_28 * dk_77[k]
                  - f_27 * dk_82[k]
                  + f_29 * dk_84[k]
                  - f_29 * dk_86[k]
                  - f_26 * dk_93[k]
                  + f_28 * dk_95[k]
                  - f_29 * dk_97[k]
                  + f_30 * dk_99[k];
    }

#pragma omp simd aligned(dk_74, dk_79, dk_81, dk_88, dk_92, dk_101, dk_103, \
                         dk_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_37 * dk_74[k]
                  + f_37 * dk_79[k]
                  - f_38 * dk_81[k]
                  - f_37 * dk_88[k]
                  + f_39 * dk_92[k]
                  - f_37 * dk_101[k]
                  + f_38 * dk_103[k]
                  - f_39 * dk_105[k];
    }

#pragma omp simd aligned(dk_72, dk_75, dk_77, dk_82, dk_84, dk_86, dk_93, dk_95, \
                         dk_97 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_17 * dk_72[k]
                  - f_17 * dk_75[k]
                  - f_20 * dk_77[k]
                  - f_15 * dk_82[k]
                  + f_18 * dk_84[k]
                  + f_21 * dk_86[k]
                  - f_14 * dk_93[k]
                  + f_16 * dk_95[k]
                  - f_19 * dk_97[k];
    }

#pragma omp simd aligned(dk_72, dk_74, dk_75, dk_77, dk_79, dk_81, dk_82, dk_84, dk_88, dk_90, \
                         dk_93, dk_95, dk_101, dk_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = -f_40 * dk_74[k]
                  + f_41 * dk_79[k]
                  + f_42 * dk_81[k]
                  + f_41 * dk_88[k]
                  - f_9 * dk_90[k]
                  - f_40 * dk_101[k]
                  + f_42 * dk_103[k];

        g_57[k] = -f_10 * dk_72[k]
                  + f_8 * dk_75[k]
                  + f_11 * dk_77[k]
                  + f_6 * dk_82[k]
                  - f_9 * dk_84[k]
                  - f_6 * dk_93[k]
                  + f_7 * dk_95[k];
    }

#pragma omp simd aligned(dk_72, dk_74, dk_75, dk_79, dk_82, dk_88, dk_93, \
                         dk_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = f_43 * dk_74[k]
                  - f_44 * dk_79[k]
                  + f_44 * dk_88[k]
                  - f_43 * dk_101[k];

        g_59[k] = f_3 * dk_72[k]
                  - f_2 * dk_75[k]
                  + f_1 * dk_82[k]
                  - f_0 * dk_93[k];
    }

#pragma omp simd aligned(dk_1, dk_4, dk_6, dk_11, dk_15, dk_22, dk_28, dk_109, dk_112, dk_114, \
                         dk_119, dk_123, dk_130, dk_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_108 * dk_1[k]
                  - f_109 * dk_6[k]
                  + f_110 * dk_15[k]
                  - f_111 * dk_28[k]
                  - f_108 * dk_109[k]
                  + f_109 * dk_114[k]
                  - f_110 * dk_123[k]
                  + f_111 * dk_136[k];

        g_61[k] = f_112 * dk_4[k]
                  - f_113 * dk_11[k]
                  + f_112 * dk_22[k]
                  - f_112 * dk_112[k]
                  + f_113 * dk_119[k]
                  - f_112 * dk_130[k];
    }

#pragma omp simd aligned(dk_1, dk_6, dk_8, dk_15, dk_17, dk_28, dk_30, dk_109, dk_114, dk_116, \
                         dk_123, dk_125, dk_136, dk_138 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_114 * dk_1[k]
                  + f_114 * dk_6[k]
                  + f_41 * dk_8[k]
                  + f_115 * dk_15[k]
                  - f_7 * dk_17[k]
                  - f_116 * dk_28[k]
                  + f_40 * dk_30[k]
                  + f_114 * dk_109[k]
                  - f_114 * dk_114[k]
                  - f_41 * dk_116[k]
                  - f_115 * dk_123[k]
                  + f_7 * dk_125[k]
                  + f_116 * dk_136[k]
                  - f_40 * dk_138[k];
    }

#pragma omp simd aligned(dk_4, dk_13, dk_22, dk_24, dk_112, dk_121, dk_130, \
                         dk_132 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_11 * dk_4[k]
                  + f_117 * dk_13[k]
                  + f_11 * dk_22[k]
                  - f_117 * dk_24[k]
                  + f_11 * dk_112[k]
                  - f_117 * dk_121[k]
                  - f_11 * dk_130[k]
                  + f_117 * dk_132[k];
    }

#pragma omp simd aligned(dk_1, dk_6, dk_8, dk_15, dk_17, dk_19, dk_28, dk_30, dk_32, dk_109, \
                         dk_114, dk_116, dk_123, dk_125, dk_127, dk_136, dk_138, \
                         dk_140 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_118 * dk_1[k]
                  + f_119 * dk_6[k]
                  - f_120 * dk_8[k]
                  + f_121 * dk_15[k]
                  - f_20 * dk_17[k]
                  + f_18 * dk_19[k]
                  - f_121 * dk_28[k]
                  + f_122 * dk_30[k]
                  - f_92 * dk_32[k]
                  - f_118 * dk_109[k]
                  - f_119 * dk_114[k]
                  + f_120 * dk_116[k]
                  - f_121 * dk_123[k]
                  + f_20 * dk_125[k]
                  - f_18 * dk_127[k]
                  + f_121 * dk_136[k]
                  - f_122 * dk_138[k]
                  + f_92 * dk_140[k];
    }

#pragma omp simd aligned(dk_4, dk_11, dk_13, dk_22, dk_24, dk_26, dk_112, dk_119, dk_121, \
                         dk_130, dk_132, dk_134 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_37 * dk_4[k]
                  + f_22 * dk_11[k]
                  - f_38 * dk_13[k]
                  + f_37 * dk_22[k]
                  - f_38 * dk_24[k]
                  + f_39 * dk_26[k]
                  - f_37 * dk_112[k]
                  - f_22 * dk_119[k]
                  + f_38 * dk_121[k]
                  - f_37 * dk_130[k]
                  + f_38 * dk_132[k]
                  - f_39 * dk_134[k];
    }

#pragma omp simd aligned(dk_1, dk_6, dk_8, dk_15, dk_17, dk_19, dk_28, dk_30, dk_32, dk_34, \
                         dk_109, dk_114, dk_116, dk_123, dk_125, dk_127, dk_136, dk_138, \
                         dk_140, dk_142 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = -f_123 * dk_1[k]
                  - f_72 * dk_6[k]
                  + f_75 * dk_8[k]
                  - f_72 * dk_15[k]
                  + f_28 * dk_17[k]
                  - f_28 * dk_19[k]
                  - f_123 * dk_28[k]
                  + f_75 * dk_30[k]
                  - f_28 * dk_32[k]
                  + f_124 * dk_34[k]
                  + f_123 * dk_109[k]
                  + f_72 * dk_114[k]
                  - f_75 * dk_116[k]
                  + f_72 * dk_123[k]
                  - f_28 * dk_125[k]
                  + f_28 * dk_127[k]
                  + f_123 * dk_136[k]
                  - f_75 * dk_138[k]
                  + f_28 * dk_140[k]
                  - f_124 * dk_142[k];
    }

#pragma omp simd aligned(dk_2, dk_7, dk_9, dk_16, dk_18, dk_20, dk_29, dk_31, dk_33, dk_35, \
                         dk_110, dk_115, dk_117, dk_124, dk_126, dk_128, dk_137, dk_139, \
                         dk_141, dk_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_125 * dk_2[k]
                  - f_126 * dk_7[k]
                  + f_32 * dk_9[k]
                  - f_126 * dk_16[k]
                  + f_33 * dk_18[k]
                  - f_127 * dk_20[k]
                  - f_125 * dk_29[k]
                  + f_32 * dk_31[k]
                  - f_127 * dk_33[k]
                  + f_128 * dk_35[k]
                  + f_125 * dk_110[k]
                  + f_126 * dk_115[k]
                  - f_32 * dk_117[k]
                  + f_126 * dk_124[k]
                  - f_33 * dk_126[k]
                  + f_127 * dk_128[k]
                  + f_125 * dk_137[k]
                  - f_32 * dk_139[k]
                  + f_127 * dk_141[k]
                  - f_128 * dk_143[k];
    }

#pragma omp simd aligned(dk_0, dk_3, dk_5, dk_10, dk_12, dk_14, dk_21, dk_23, dk_25, dk_27, \
                         dk_108, dk_111, dk_113, dk_118, dk_120, dk_122, dk_129, dk_131, \
                         dk_133, dk_135 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = -f_123 * dk_0[k]
                  - f_72 * dk_3[k]
                  + f_75 * dk_5[k]
                  - f_72 * dk_10[k]
                  + f_28 * dk_12[k]
                  - f_28 * dk_14[k]
                  - f_123 * dk_21[k]
                  + f_75 * dk_23[k]
                  - f_28 * dk_25[k]
                  + f_124 * dk_27[k]
                  + f_123 * dk_108[k]
                  + f_72 * dk_111[k]
                  - f_75 * dk_113[k]
                  + f_72 * dk_118[k]
                  - f_28 * dk_120[k]
                  + f_28 * dk_122[k]
                  + f_123 * dk_129[k]
                  - f_75 * dk_131[k]
                  + f_28 * dk_133[k]
                  - f_124 * dk_135[k];
    }

#pragma omp simd aligned(dk_2, dk_7, dk_9, dk_16, dk_20, dk_29, dk_31, dk_33, dk_110, dk_115, \
                         dk_117, dk_124, dk_128, dk_137, dk_139, \
                         dk_141 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_129 * dk_2[k]
                  + f_129 * dk_7[k]
                  - f_130 * dk_9[k]
                  - f_129 * dk_16[k]
                  + f_131 * dk_20[k]
                  - f_129 * dk_29[k]
                  + f_130 * dk_31[k]
                  - f_131 * dk_33[k]
                  - f_129 * dk_110[k]
                  - f_129 * dk_115[k]
                  + f_130 * dk_117[k]
                  + f_129 * dk_124[k]
                  - f_131 * dk_128[k]
                  + f_129 * dk_137[k]
                  - f_130 * dk_139[k]
                  + f_131 * dk_141[k];
    }

#pragma omp simd aligned(dk_0, dk_3, dk_5, dk_10, dk_12, dk_14, dk_21, dk_23, dk_25, dk_108, \
                         dk_111, dk_113, dk_118, dk_120, dk_122, dk_129, dk_131, \
                         dk_133 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = f_121 * dk_0[k]
                  - f_121 * dk_3[k]
                  - f_122 * dk_5[k]
                  - f_119 * dk_10[k]
                  + f_20 * dk_12[k]
                  + f_92 * dk_14[k]
                  - f_118 * dk_21[k]
                  + f_120 * dk_23[k]
                  - f_18 * dk_25[k]
                  - f_121 * dk_108[k]
                  + f_121 * dk_111[k]
                  + f_122 * dk_113[k]
                  + f_119 * dk_118[k]
                  - f_20 * dk_120[k]
                  - f_92 * dk_122[k]
                  + f_118 * dk_129[k]
                  - f_120 * dk_131[k]
                  + f_18 * dk_133[k];
    }

#pragma omp simd aligned(dk_2, dk_7, dk_9, dk_16, dk_18, dk_29, dk_31, dk_110, dk_115, dk_117, \
                         dk_124, dk_126, dk_137, dk_139 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_132 * dk_2[k]
                  + f_133 * dk_7[k]
                  + f_134 * dk_9[k]
                  + f_133 * dk_16[k]
                  - f_7 * dk_18[k]
                  - f_132 * dk_29[k]
                  + f_134 * dk_31[k]
                  + f_132 * dk_110[k]
                  - f_133 * dk_115[k]
                  - f_134 * dk_117[k]
                  - f_133 * dk_124[k]
                  + f_7 * dk_126[k]
                  + f_132 * dk_137[k]
                  - f_134 * dk_139[k];
    }

#pragma omp simd aligned(dk_0, dk_3, dk_5, dk_10, dk_12, dk_21, dk_23, dk_108, dk_111, dk_113, \
                         dk_118, dk_120, dk_129, dk_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -f_116 * dk_0[k]
                  + f_115 * dk_3[k]
                  + f_40 * dk_5[k]
                  + f_114 * dk_10[k]
                  - f_7 * dk_12[k]
                  - f_114 * dk_21[k]
                  + f_41 * dk_23[k]
                  + f_116 * dk_108[k]
                  - f_115 * dk_111[k]
                  - f_40 * dk_113[k]
                  - f_114 * dk_118[k]
                  + f_7 * dk_120[k]
                  + f_114 * dk_129[k]
                  - f_41 * dk_131[k];
    }

#pragma omp simd aligned(dk_2, dk_7, dk_16, dk_29, dk_110, dk_115, dk_124, \
                         dk_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_135 * dk_2[k]
                  - f_136 * dk_7[k]
                  + f_136 * dk_16[k]
                  - f_135 * dk_29[k]
                  - f_135 * dk_110[k]
                  + f_136 * dk_115[k]
                  - f_136 * dk_124[k]
                  + f_135 * dk_137[k];
    }

#pragma omp simd aligned(dk_0, dk_3, dk_10, dk_21, dk_108, dk_111, dk_118, \
                         dk_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_111 * dk_0[k]
                  - f_110 * dk_3[k]
                  + f_109 * dk_10[k]
                  - f_108 * dk_21[k]
                  - f_111 * dk_108[k]
                  + f_110 * dk_111[k]
                  - f_109 * dk_118[k]
                  + f_108 * dk_129[k];
    }
}

}  // namespace simdtrf
