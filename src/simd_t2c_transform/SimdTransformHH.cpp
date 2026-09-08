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


#include "SimdTransformHH.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_hh(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t hh,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 9.84375 * std::sqrt(10.0);
    const auto f_1 = 19.6875 * std::sqrt(10.0);
    const auto f_2 = 1.96875 * std::sqrt(10.0);
    const auto f_3 = 2.4609375 * std::sqrt(5.0);
    const auto f_4 = 1.640625 * std::sqrt(5.0);
    const auto f_5 = 19.6875 * std::sqrt(5.0);
    const auto f_6 = 0.8203125 * std::sqrt(5.0);
    const auto f_7 = 6.5625 * std::sqrt(5.0);
    const auto f_8 = 4.921875 * std::sqrt(5.0);
    const auto f_9 = 3.28125 * std::sqrt(5.0);
    const auto f_10 = 39.375 * std::sqrt(5.0);
    const auto f_11 = 13.125 * std::sqrt(5.0);
    const auto f_12 = 0.4921875 * std::sqrt(5.0);
    const auto f_13 = 0.328125 * std::sqrt(5.0);
    const auto f_14 = 3.9375 * std::sqrt(5.0);
    const auto f_15 = 0.1640625 * std::sqrt(5.0);
    const auto f_16 = 1.3125 * std::sqrt(5.0);
    const auto f_17 = 3.28125 * std::sqrt(30.0);
    const auto f_18 = 6.5625 * std::sqrt(30.0);
    const auto f_19 = 13.125 * std::sqrt(30.0);
    const auto f_20 = 0.65625 * std::sqrt(30.0);
    const auto f_21 = 1.3125 * std::sqrt(30.0);
    const auto f_22 = 0.1171875 * std::sqrt(210.0);
    const auto f_23 = 0.234375 * std::sqrt(210.0);
    const auto f_24 = 1.40625 * std::sqrt(210.0);
    const auto f_25 = 0.9375 * std::sqrt(210.0);
    const auto f_26 = 0.46875 * std::sqrt(210.0);
    const auto f_27 = 2.8125 * std::sqrt(210.0);
    const auto f_28 = 1.875 * std::sqrt(210.0);
    const auto f_29 = 0.0234375 * std::sqrt(210.0);
    const auto f_30 = 0.046875 * std::sqrt(210.0);
    const auto f_31 = 0.28125 * std::sqrt(210.0);
    const auto f_32 = 0.1875 * std::sqrt(210.0);
    const auto f_33 = 1.7578125 * std::sqrt(14.0);
    const auto f_34 = 3.515625 * std::sqrt(14.0);
    const auto f_35 = 4.6875 * std::sqrt(14.0);
    const auto f_36 = 0.9375 * std::sqrt(14.0);
    const auto f_37 = 7.03125 * std::sqrt(14.0);
    const auto f_38 = 9.375 * std::sqrt(14.0);
    const auto f_39 = 1.875 * std::sqrt(14.0);
    const auto f_40 = 0.3515625 * std::sqrt(14.0);
    const auto f_41 = 0.703125 * std::sqrt(14.0);
    const auto f_42 = 0.1875 * std::sqrt(14.0);
    const auto f_43 = 1.640625 * std::sqrt(30.0);
    const auto f_44 = 0.328125 * std::sqrt(30.0);
    const auto f_45 = 2.4609375 * std::sqrt(10.0);
    const auto f_46 = 14.765625 * std::sqrt(10.0);
    const auto f_47 = 4.921875 * std::sqrt(10.0);
    const auto f_48 = 29.53125 * std::sqrt(10.0);
    const auto f_49 = 0.4921875 * std::sqrt(10.0);
    const auto f_50 = 2.953125 * std::sqrt(10.0);
    const auto f_51 = 9.84375 * std::sqrt(2.0);
    const auto f_52 = 6.5625 * std::sqrt(2.0);
    const auto f_53 = 78.75 * std::sqrt(2.0);
    const auto f_54 = 3.28125 * std::sqrt(2.0);
    const auto f_55 = 26.25 * std::sqrt(2.0);
    const auto f_56 = 26.25 * std::sqrt(3.0);
    const auto f_57 = 52.5 * std::sqrt(3.0);
    const auto f_58 = 0.9375 * std::sqrt(21.0);
    const auto f_59 = 1.875 * std::sqrt(21.0);
    const auto f_60 = 11.25 * std::sqrt(21.0);
    const auto f_61 = 7.5 * std::sqrt(21.0);
    const auto f_62 = 2.8125 * std::sqrt(35.0);
    const auto f_63 = 5.625 * std::sqrt(35.0);
    const auto f_64 = 7.5 * std::sqrt(35.0);
    const auto f_65 = 1.5 * std::sqrt(35.0);
    const auto f_66 = 13.125 * std::sqrt(3.0);
    const auto f_67 = 3.28125 * std::sqrt(6.0);
    const auto f_68 = 6.5625 * std::sqrt(6.0);
    const auto f_69 = 2.1875 * std::sqrt(6.0);
    const auto f_70 = 4.375 * std::sqrt(6.0);
    const auto f_71 = 26.25 * std::sqrt(6.0);
    const auto f_72 = 52.5 * std::sqrt(6.0);
    const auto f_73 = 1.09375 * std::sqrt(6.0);
    const auto f_74 = 8.75 * std::sqrt(6.0);
    const auto f_75 = 17.5 * std::sqrt(6.0);
    const auto f_76 = 0.1171875 * std::sqrt(42.0);
    const auto f_77 = 0.234375 * std::sqrt(42.0);
    const auto f_78 = 1.40625 * std::sqrt(42.0);
    const auto f_79 = 0.9375 * std::sqrt(42.0);
    const auto f_80 = 0.078125 * std::sqrt(42.0);
    const auto f_81 = 0.15625 * std::sqrt(42.0);
    const auto f_82 = 0.625 * std::sqrt(42.0);
    const auto f_83 = 1.875 * std::sqrt(42.0);
    const auto f_84 = 11.25 * std::sqrt(42.0);
    const auto f_85 = 7.5 * std::sqrt(42.0);
    const auto f_86 = 0.0390625 * std::sqrt(42.0);
    const auto f_87 = 0.46875 * std::sqrt(42.0);
    const auto f_88 = 0.3125 * std::sqrt(42.0);
    const auto f_89 = 3.75 * std::sqrt(42.0);
    const auto f_90 = 2.5 * std::sqrt(42.0);
    const auto f_91 = 0.3515625 * std::sqrt(70.0);
    const auto f_92 = 0.703125 * std::sqrt(70.0);
    const auto f_93 = 0.9375 * std::sqrt(70.0);
    const auto f_94 = 0.1875 * std::sqrt(70.0);
    const auto f_95 = 0.234375 * std::sqrt(70.0);
    const auto f_96 = 0.46875 * std::sqrt(70.0);
    const auto f_97 = 0.625 * std::sqrt(70.0);
    const auto f_98 = 0.125 * std::sqrt(70.0);
    const auto f_99 = 2.8125 * std::sqrt(70.0);
    const auto f_100 = 5.625 * std::sqrt(70.0);
    const auto f_101 = 7.5 * std::sqrt(70.0);
    const auto f_102 = 1.5 * std::sqrt(70.0);
    const auto f_103 = 0.1171875 * std::sqrt(70.0);
    const auto f_104 = 0.3125 * std::sqrt(70.0);
    const auto f_105 = 0.0625 * std::sqrt(70.0);
    const auto f_106 = 1.875 * std::sqrt(70.0);
    const auto f_107 = 2.5 * std::sqrt(70.0);
    const auto f_108 = 0.5 * std::sqrt(70.0);
    const auto f_109 = 1.640625 * std::sqrt(6.0);
    const auto f_110 = 13.125 * std::sqrt(6.0);
    const auto f_111 = 0.546875 * std::sqrt(6.0);
    const auto f_112 = 2.4609375 * std::sqrt(2.0);
    const auto f_113 = 14.765625 * std::sqrt(2.0);
    const auto f_114 = 1.640625 * std::sqrt(2.0);
    const auto f_115 = 19.6875 * std::sqrt(2.0);
    const auto f_116 = 118.125 * std::sqrt(2.0);
    const auto f_117 = 0.8203125 * std::sqrt(2.0);
    const auto f_118 = 4.921875 * std::sqrt(2.0);
    const auto f_119 = 39.375 * std::sqrt(2.0);
    const auto f_120 = 0.9375 * std::sqrt(7.0);
    const auto f_121 = 1.875 * std::sqrt(7.0);
    const auto f_122 = 11.25 * std::sqrt(7.0);
    const auto f_123 = 7.5 * std::sqrt(7.0);
    const auto f_124 = 3.75 * std::sqrt(7.0);
    const auto f_125 = 22.5 * std::sqrt(7.0);
    const auto f_126 = 15.0 * std::sqrt(7.0);
    const auto f_127 = 0.9375 * std::sqrt(105.0);
    const auto f_128 = 1.875 * std::sqrt(105.0);
    const auto f_129 = 2.5 * std::sqrt(105.0);
    const auto f_130 = 0.5 * std::sqrt(105.0);
    const auto f_131 = 3.75 * std::sqrt(105.0);
    const auto f_132 = 5.0 * std::sqrt(105.0);
    const auto f_133 = std::sqrt(105.0);
    const auto f_134 = 6.5625 * std::sqrt(3.0);
    const auto f_135 = 39.375 * std::sqrt(3.0);
    const auto f_136 = 78.75 * std::sqrt(3.0);
    const auto f_137 = 0.234375 * std::sqrt(15.0);
    const auto f_138 = 0.46875 * std::sqrt(15.0);
    const auto f_139 = 0.625 * std::sqrt(15.0);
    const auto f_140 = 0.125 * std::sqrt(15.0);
    const auto f_141 = 0.9375 * std::sqrt(15.0);
    const auto f_142 = 1.25 * std::sqrt(15.0);
    const auto f_143 = 0.25 * std::sqrt(15.0);
    const auto f_144 = 2.8125 * std::sqrt(15.0);
    const auto f_145 = 5.625 * std::sqrt(15.0);
    const auto f_146 = 7.5 * std::sqrt(15.0);
    const auto f_147 = 1.5 * std::sqrt(15.0);
    const auto f_148 = 1.875 * std::sqrt(15.0);
    const auto f_149 = 3.75 * std::sqrt(15.0);
    const auto f_150 = 5.0 * std::sqrt(15.0);
    const auto f_151 = std::sqrt(15.0);
    const auto f_152 = 0.46875 * std::sqrt(7.0);
    const auto f_153 = 5.625 * std::sqrt(7.0);
    const auto f_154 = 0.234375 * std::sqrt(21.0);
    const auto f_155 = 1.40625 * std::sqrt(21.0);
    const auto f_156 = 0.46875 * std::sqrt(21.0);
    const auto f_157 = 2.8125 * std::sqrt(21.0);
    const auto f_158 = 16.875 * std::sqrt(21.0);
    const auto f_159 = 0.46875 * std::sqrt(105.0);
    const auto f_160 = 1.25 * std::sqrt(105.0);
    const auto f_161 = 0.25 * std::sqrt(105.0);
    const auto f_162 = 0.703125 * std::sqrt(35.0);
    const auto f_163 = 4.21875 * std::sqrt(35.0);
    const auto f_164 = 1.40625 * std::sqrt(35.0);
    const auto f_165 = 8.4375 * std::sqrt(35.0);
    const auto f_166 = 1.875 * std::sqrt(35.0);
    const auto f_167 = 11.25 * std::sqrt(35.0);
    const auto f_168 = 0.375 * std::sqrt(35.0);
    const auto f_169 = 2.25 * std::sqrt(35.0);
    const auto f_170 = 3.28125 * std::sqrt(3.0);
    const auto f_171 = 19.6875 * std::sqrt(3.0);

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
    auto *g_85 = values + 85 * nvalues;
    auto *g_86 = values + 86 * nvalues;
    auto *g_87 = values + 87 * nvalues;
    auto *g_88 = values + 88 * nvalues;
    auto *g_89 = values + 89 * nvalues;
    auto *g_90 = values + 90 * nvalues;
    auto *g_91 = values + 91 * nvalues;
    auto *g_92 = values + 92 * nvalues;
    auto *g_93 = values + 93 * nvalues;
    auto *g_94 = values + 94 * nvalues;
    auto *g_95 = values + 95 * nvalues;
    auto *g_96 = values + 96 * nvalues;
    auto *g_97 = values + 97 * nvalues;
    auto *g_98 = values + 98 * nvalues;
    auto *g_99 = values + 99 * nvalues;
    auto *g_100 = values + 100 * nvalues;
    auto *g_101 = values + 101 * nvalues;
    auto *g_102 = values + 102 * nvalues;
    auto *g_103 = values + 103 * nvalues;
    auto *g_104 = values + 104 * nvalues;
    auto *g_105 = values + 105 * nvalues;
    auto *g_106 = values + 106 * nvalues;
    auto *g_107 = values + 107 * nvalues;
    auto *g_108 = values + 108 * nvalues;
    auto *g_109 = values + 109 * nvalues;
    auto *g_110 = values + 110 * nvalues;
    auto *g_111 = values + 111 * nvalues;
    auto *g_112 = values + 112 * nvalues;
    auto *g_113 = values + 113 * nvalues;
    auto *g_114 = values + 114 * nvalues;
    auto *g_115 = values + 115 * nvalues;
    auto *g_116 = values + 116 * nvalues;
    auto *g_117 = values + 117 * nvalues;
    auto *g_118 = values + 118 * nvalues;
    auto *g_119 = values + 119 * nvalues;
    auto *g_120 = values + 120 * nvalues;

    const auto *hh_0 = buffer.data(hh + 0);
    const auto *hh_1 = buffer.data(hh + 1);
    const auto *hh_2 = buffer.data(hh + 2);
    const auto *hh_3 = buffer.data(hh + 3);
    const auto *hh_4 = buffer.data(hh + 4);
    const auto *hh_5 = buffer.data(hh + 5);
    const auto *hh_6 = buffer.data(hh + 6);
    const auto *hh_7 = buffer.data(hh + 7);
    const auto *hh_8 = buffer.data(hh + 8);
    const auto *hh_9 = buffer.data(hh + 9);
    const auto *hh_10 = buffer.data(hh + 10);
    const auto *hh_11 = buffer.data(hh + 11);
    const auto *hh_12 = buffer.data(hh + 12);
    const auto *hh_13 = buffer.data(hh + 13);
    const auto *hh_14 = buffer.data(hh + 14);
    const auto *hh_15 = buffer.data(hh + 15);
    const auto *hh_16 = buffer.data(hh + 16);
    const auto *hh_17 = buffer.data(hh + 17);
    const auto *hh_18 = buffer.data(hh + 18);
    const auto *hh_19 = buffer.data(hh + 19);
    const auto *hh_20 = buffer.data(hh + 20);
    const auto *hh_21 = buffer.data(hh + 21);
    const auto *hh_22 = buffer.data(hh + 22);
    const auto *hh_23 = buffer.data(hh + 23);
    const auto *hh_24 = buffer.data(hh + 24);
    const auto *hh_25 = buffer.data(hh + 25);
    const auto *hh_26 = buffer.data(hh + 26);
    const auto *hh_27 = buffer.data(hh + 27);
    const auto *hh_28 = buffer.data(hh + 28);
    const auto *hh_29 = buffer.data(hh + 29);
    const auto *hh_30 = buffer.data(hh + 30);
    const auto *hh_31 = buffer.data(hh + 31);
    const auto *hh_32 = buffer.data(hh + 32);
    const auto *hh_33 = buffer.data(hh + 33);
    const auto *hh_34 = buffer.data(hh + 34);
    const auto *hh_35 = buffer.data(hh + 35);
    const auto *hh_36 = buffer.data(hh + 36);
    const auto *hh_37 = buffer.data(hh + 37);
    const auto *hh_38 = buffer.data(hh + 38);
    const auto *hh_39 = buffer.data(hh + 39);
    const auto *hh_40 = buffer.data(hh + 40);
    const auto *hh_41 = buffer.data(hh + 41);
    const auto *hh_42 = buffer.data(hh + 42);
    const auto *hh_43 = buffer.data(hh + 43);
    const auto *hh_44 = buffer.data(hh + 44);
    const auto *hh_45 = buffer.data(hh + 45);
    const auto *hh_46 = buffer.data(hh + 46);
    const auto *hh_47 = buffer.data(hh + 47);
    const auto *hh_48 = buffer.data(hh + 48);
    const auto *hh_49 = buffer.data(hh + 49);
    const auto *hh_50 = buffer.data(hh + 50);
    const auto *hh_51 = buffer.data(hh + 51);
    const auto *hh_52 = buffer.data(hh + 52);
    const auto *hh_53 = buffer.data(hh + 53);
    const auto *hh_54 = buffer.data(hh + 54);
    const auto *hh_55 = buffer.data(hh + 55);
    const auto *hh_56 = buffer.data(hh + 56);
    const auto *hh_57 = buffer.data(hh + 57);
    const auto *hh_58 = buffer.data(hh + 58);
    const auto *hh_59 = buffer.data(hh + 59);
    const auto *hh_60 = buffer.data(hh + 60);
    const auto *hh_61 = buffer.data(hh + 61);
    const auto *hh_62 = buffer.data(hh + 62);
    const auto *hh_63 = buffer.data(hh + 63);
    const auto *hh_64 = buffer.data(hh + 64);
    const auto *hh_65 = buffer.data(hh + 65);
    const auto *hh_66 = buffer.data(hh + 66);
    const auto *hh_67 = buffer.data(hh + 67);
    const auto *hh_68 = buffer.data(hh + 68);
    const auto *hh_69 = buffer.data(hh + 69);
    const auto *hh_70 = buffer.data(hh + 70);
    const auto *hh_71 = buffer.data(hh + 71);
    const auto *hh_72 = buffer.data(hh + 72);
    const auto *hh_73 = buffer.data(hh + 73);
    const auto *hh_74 = buffer.data(hh + 74);
    const auto *hh_75 = buffer.data(hh + 75);
    const auto *hh_76 = buffer.data(hh + 76);
    const auto *hh_77 = buffer.data(hh + 77);
    const auto *hh_78 = buffer.data(hh + 78);
    const auto *hh_79 = buffer.data(hh + 79);
    const auto *hh_80 = buffer.data(hh + 80);
    const auto *hh_81 = buffer.data(hh + 81);
    const auto *hh_82 = buffer.data(hh + 82);
    const auto *hh_83 = buffer.data(hh + 83);
    const auto *hh_84 = buffer.data(hh + 84);
    const auto *hh_85 = buffer.data(hh + 85);
    const auto *hh_86 = buffer.data(hh + 86);
    const auto *hh_87 = buffer.data(hh + 87);
    const auto *hh_88 = buffer.data(hh + 88);
    const auto *hh_89 = buffer.data(hh + 89);
    const auto *hh_90 = buffer.data(hh + 90);
    const auto *hh_91 = buffer.data(hh + 91);
    const auto *hh_92 = buffer.data(hh + 92);
    const auto *hh_93 = buffer.data(hh + 93);
    const auto *hh_94 = buffer.data(hh + 94);
    const auto *hh_95 = buffer.data(hh + 95);
    const auto *hh_96 = buffer.data(hh + 96);
    const auto *hh_97 = buffer.data(hh + 97);
    const auto *hh_98 = buffer.data(hh + 98);
    const auto *hh_99 = buffer.data(hh + 99);
    const auto *hh_100 = buffer.data(hh + 100);
    const auto *hh_101 = buffer.data(hh + 101);
    const auto *hh_102 = buffer.data(hh + 102);
    const auto *hh_103 = buffer.data(hh + 103);
    const auto *hh_104 = buffer.data(hh + 104);
    const auto *hh_105 = buffer.data(hh + 105);
    const auto *hh_106 = buffer.data(hh + 106);
    const auto *hh_107 = buffer.data(hh + 107);
    const auto *hh_108 = buffer.data(hh + 108);
    const auto *hh_109 = buffer.data(hh + 109);
    const auto *hh_110 = buffer.data(hh + 110);
    const auto *hh_111 = buffer.data(hh + 111);
    const auto *hh_112 = buffer.data(hh + 112);
    const auto *hh_113 = buffer.data(hh + 113);
    const auto *hh_114 = buffer.data(hh + 114);
    const auto *hh_115 = buffer.data(hh + 115);
    const auto *hh_116 = buffer.data(hh + 116);
    const auto *hh_117 = buffer.data(hh + 117);
    const auto *hh_118 = buffer.data(hh + 118);
    const auto *hh_119 = buffer.data(hh + 119);
    const auto *hh_120 = buffer.data(hh + 120);
    const auto *hh_121 = buffer.data(hh + 121);
    const auto *hh_122 = buffer.data(hh + 122);
    const auto *hh_123 = buffer.data(hh + 123);
    const auto *hh_124 = buffer.data(hh + 124);
    const auto *hh_125 = buffer.data(hh + 125);
    const auto *hh_126 = buffer.data(hh + 126);
    const auto *hh_127 = buffer.data(hh + 127);
    const auto *hh_128 = buffer.data(hh + 128);
    const auto *hh_129 = buffer.data(hh + 129);
    const auto *hh_130 = buffer.data(hh + 130);
    const auto *hh_131 = buffer.data(hh + 131);
    const auto *hh_132 = buffer.data(hh + 132);
    const auto *hh_133 = buffer.data(hh + 133);
    const auto *hh_134 = buffer.data(hh + 134);
    const auto *hh_135 = buffer.data(hh + 135);
    const auto *hh_136 = buffer.data(hh + 136);
    const auto *hh_137 = buffer.data(hh + 137);
    const auto *hh_138 = buffer.data(hh + 138);
    const auto *hh_139 = buffer.data(hh + 139);
    const auto *hh_140 = buffer.data(hh + 140);
    const auto *hh_141 = buffer.data(hh + 141);
    const auto *hh_142 = buffer.data(hh + 142);
    const auto *hh_143 = buffer.data(hh + 143);
    const auto *hh_144 = buffer.data(hh + 144);
    const auto *hh_145 = buffer.data(hh + 145);
    const auto *hh_146 = buffer.data(hh + 146);
    const auto *hh_147 = buffer.data(hh + 147);
    const auto *hh_148 = buffer.data(hh + 148);
    const auto *hh_149 = buffer.data(hh + 149);
    const auto *hh_150 = buffer.data(hh + 150);
    const auto *hh_151 = buffer.data(hh + 151);
    const auto *hh_152 = buffer.data(hh + 152);
    const auto *hh_153 = buffer.data(hh + 153);
    const auto *hh_154 = buffer.data(hh + 154);
    const auto *hh_155 = buffer.data(hh + 155);
    const auto *hh_156 = buffer.data(hh + 156);
    const auto *hh_157 = buffer.data(hh + 157);
    const auto *hh_158 = buffer.data(hh + 158);
    const auto *hh_159 = buffer.data(hh + 159);
    const auto *hh_160 = buffer.data(hh + 160);
    const auto *hh_161 = buffer.data(hh + 161);
    const auto *hh_162 = buffer.data(hh + 162);
    const auto *hh_163 = buffer.data(hh + 163);
    const auto *hh_164 = buffer.data(hh + 164);
    const auto *hh_165 = buffer.data(hh + 165);
    const auto *hh_166 = buffer.data(hh + 166);
    const auto *hh_167 = buffer.data(hh + 167);
    const auto *hh_168 = buffer.data(hh + 168);
    const auto *hh_169 = buffer.data(hh + 169);
    const auto *hh_170 = buffer.data(hh + 170);
    const auto *hh_171 = buffer.data(hh + 171);
    const auto *hh_172 = buffer.data(hh + 172);
    const auto *hh_173 = buffer.data(hh + 173);
    const auto *hh_174 = buffer.data(hh + 174);
    const auto *hh_175 = buffer.data(hh + 175);
    const auto *hh_176 = buffer.data(hh + 176);
    const auto *hh_177 = buffer.data(hh + 177);
    const auto *hh_178 = buffer.data(hh + 178);
    const auto *hh_179 = buffer.data(hh + 179);
    const auto *hh_180 = buffer.data(hh + 180);
    const auto *hh_181 = buffer.data(hh + 181);
    const auto *hh_182 = buffer.data(hh + 182);
    const auto *hh_183 = buffer.data(hh + 183);
    const auto *hh_184 = buffer.data(hh + 184);
    const auto *hh_185 = buffer.data(hh + 185);
    const auto *hh_186 = buffer.data(hh + 186);
    const auto *hh_187 = buffer.data(hh + 187);
    const auto *hh_188 = buffer.data(hh + 188);
    const auto *hh_189 = buffer.data(hh + 189);
    const auto *hh_190 = buffer.data(hh + 190);
    const auto *hh_191 = buffer.data(hh + 191);
    const auto *hh_192 = buffer.data(hh + 192);
    const auto *hh_193 = buffer.data(hh + 193);
    const auto *hh_194 = buffer.data(hh + 194);
    const auto *hh_195 = buffer.data(hh + 195);
    const auto *hh_196 = buffer.data(hh + 196);
    const auto *hh_197 = buffer.data(hh + 197);
    const auto *hh_198 = buffer.data(hh + 198);
    const auto *hh_199 = buffer.data(hh + 199);
    const auto *hh_200 = buffer.data(hh + 200);
    const auto *hh_201 = buffer.data(hh + 201);
    const auto *hh_202 = buffer.data(hh + 202);
    const auto *hh_203 = buffer.data(hh + 203);
    const auto *hh_204 = buffer.data(hh + 204);
    const auto *hh_205 = buffer.data(hh + 205);
    const auto *hh_206 = buffer.data(hh + 206);
    const auto *hh_207 = buffer.data(hh + 207);
    const auto *hh_208 = buffer.data(hh + 208);
    const auto *hh_209 = buffer.data(hh + 209);
    const auto *hh_210 = buffer.data(hh + 210);
    const auto *hh_211 = buffer.data(hh + 211);
    const auto *hh_212 = buffer.data(hh + 212);
    const auto *hh_213 = buffer.data(hh + 213);
    const auto *hh_214 = buffer.data(hh + 214);
    const auto *hh_215 = buffer.data(hh + 215);
    const auto *hh_216 = buffer.data(hh + 216);
    const auto *hh_217 = buffer.data(hh + 217);
    const auto *hh_218 = buffer.data(hh + 218);
    const auto *hh_219 = buffer.data(hh + 219);
    const auto *hh_220 = buffer.data(hh + 220);
    const auto *hh_221 = buffer.data(hh + 221);
    const auto *hh_222 = buffer.data(hh + 222);
    const auto *hh_223 = buffer.data(hh + 223);
    const auto *hh_224 = buffer.data(hh + 224);
    const auto *hh_225 = buffer.data(hh + 225);
    const auto *hh_226 = buffer.data(hh + 226);
    const auto *hh_227 = buffer.data(hh + 227);
    const auto *hh_228 = buffer.data(hh + 228);
    const auto *hh_229 = buffer.data(hh + 229);
    const auto *hh_230 = buffer.data(hh + 230);
    const auto *hh_231 = buffer.data(hh + 231);
    const auto *hh_232 = buffer.data(hh + 232);
    const auto *hh_233 = buffer.data(hh + 233);
    const auto *hh_234 = buffer.data(hh + 234);
    const auto *hh_235 = buffer.data(hh + 235);
    const auto *hh_236 = buffer.data(hh + 236);
    const auto *hh_237 = buffer.data(hh + 237);
    const auto *hh_238 = buffer.data(hh + 238);
    const auto *hh_239 = buffer.data(hh + 239);
    const auto *hh_240 = buffer.data(hh + 240);
    const auto *hh_241 = buffer.data(hh + 241);
    const auto *hh_242 = buffer.data(hh + 242);
    const auto *hh_243 = buffer.data(hh + 243);
    const auto *hh_244 = buffer.data(hh + 244);
    const auto *hh_245 = buffer.data(hh + 245);
    const auto *hh_246 = buffer.data(hh + 246);
    const auto *hh_247 = buffer.data(hh + 247);
    const auto *hh_248 = buffer.data(hh + 248);
    const auto *hh_249 = buffer.data(hh + 249);
    const auto *hh_250 = buffer.data(hh + 250);
    const auto *hh_251 = buffer.data(hh + 251);
    const auto *hh_252 = buffer.data(hh + 252);
    const auto *hh_253 = buffer.data(hh + 253);
    const auto *hh_254 = buffer.data(hh + 254);
    const auto *hh_255 = buffer.data(hh + 255);
    const auto *hh_256 = buffer.data(hh + 256);
    const auto *hh_257 = buffer.data(hh + 257);
    const auto *hh_258 = buffer.data(hh + 258);
    const auto *hh_259 = buffer.data(hh + 259);
    const auto *hh_260 = buffer.data(hh + 260);
    const auto *hh_261 = buffer.data(hh + 261);
    const auto *hh_262 = buffer.data(hh + 262);
    const auto *hh_263 = buffer.data(hh + 263);
    const auto *hh_264 = buffer.data(hh + 264);
    const auto *hh_265 = buffer.data(hh + 265);
    const auto *hh_266 = buffer.data(hh + 266);
    const auto *hh_267 = buffer.data(hh + 267);
    const auto *hh_268 = buffer.data(hh + 268);
    const auto *hh_269 = buffer.data(hh + 269);
    const auto *hh_270 = buffer.data(hh + 270);
    const auto *hh_271 = buffer.data(hh + 271);
    const auto *hh_272 = buffer.data(hh + 272);
    const auto *hh_273 = buffer.data(hh + 273);
    const auto *hh_274 = buffer.data(hh + 274);
    const auto *hh_275 = buffer.data(hh + 275);
    const auto *hh_276 = buffer.data(hh + 276);
    const auto *hh_277 = buffer.data(hh + 277);
    const auto *hh_278 = buffer.data(hh + 278);
    const auto *hh_279 = buffer.data(hh + 279);
    const auto *hh_280 = buffer.data(hh + 280);
    const auto *hh_281 = buffer.data(hh + 281);
    const auto *hh_282 = buffer.data(hh + 282);
    const auto *hh_283 = buffer.data(hh + 283);
    const auto *hh_284 = buffer.data(hh + 284);
    const auto *hh_285 = buffer.data(hh + 285);
    const auto *hh_286 = buffer.data(hh + 286);
    const auto *hh_287 = buffer.data(hh + 287);
    const auto *hh_288 = buffer.data(hh + 288);
    const auto *hh_289 = buffer.data(hh + 289);
    const auto *hh_290 = buffer.data(hh + 290);
    const auto *hh_291 = buffer.data(hh + 291);
    const auto *hh_292 = buffer.data(hh + 292);
    const auto *hh_293 = buffer.data(hh + 293);
    const auto *hh_294 = buffer.data(hh + 294);
    const auto *hh_295 = buffer.data(hh + 295);
    const auto *hh_296 = buffer.data(hh + 296);
    const auto *hh_297 = buffer.data(hh + 297);
    const auto *hh_298 = buffer.data(hh + 298);
    const auto *hh_299 = buffer.data(hh + 299);
    const auto *hh_300 = buffer.data(hh + 300);
    const auto *hh_301 = buffer.data(hh + 301);
    const auto *hh_302 = buffer.data(hh + 302);
    const auto *hh_303 = buffer.data(hh + 303);
    const auto *hh_304 = buffer.data(hh + 304);
    const auto *hh_305 = buffer.data(hh + 305);
    const auto *hh_306 = buffer.data(hh + 306);
    const auto *hh_307 = buffer.data(hh + 307);
    const auto *hh_308 = buffer.data(hh + 308);
    const auto *hh_309 = buffer.data(hh + 309);
    const auto *hh_310 = buffer.data(hh + 310);
    const auto *hh_311 = buffer.data(hh + 311);
    const auto *hh_312 = buffer.data(hh + 312);
    const auto *hh_313 = buffer.data(hh + 313);
    const auto *hh_314 = buffer.data(hh + 314);
    const auto *hh_315 = buffer.data(hh + 315);
    const auto *hh_316 = buffer.data(hh + 316);
    const auto *hh_317 = buffer.data(hh + 317);
    const auto *hh_318 = buffer.data(hh + 318);
    const auto *hh_319 = buffer.data(hh + 319);
    const auto *hh_320 = buffer.data(hh + 320);
    const auto *hh_321 = buffer.data(hh + 321);
    const auto *hh_322 = buffer.data(hh + 322);
    const auto *hh_323 = buffer.data(hh + 323);
    const auto *hh_324 = buffer.data(hh + 324);
    const auto *hh_325 = buffer.data(hh + 325);
    const auto *hh_326 = buffer.data(hh + 326);
    const auto *hh_327 = buffer.data(hh + 327);
    const auto *hh_328 = buffer.data(hh + 328);
    const auto *hh_329 = buffer.data(hh + 329);
    const auto *hh_330 = buffer.data(hh + 330);
    const auto *hh_331 = buffer.data(hh + 331);
    const auto *hh_332 = buffer.data(hh + 332);
    const auto *hh_333 = buffer.data(hh + 333);
    const auto *hh_334 = buffer.data(hh + 334);
    const auto *hh_335 = buffer.data(hh + 335);
    const auto *hh_336 = buffer.data(hh + 336);
    const auto *hh_337 = buffer.data(hh + 337);
    const auto *hh_338 = buffer.data(hh + 338);
    const auto *hh_339 = buffer.data(hh + 339);
    const auto *hh_340 = buffer.data(hh + 340);
    const auto *hh_341 = buffer.data(hh + 341);
    const auto *hh_342 = buffer.data(hh + 342);
    const auto *hh_343 = buffer.data(hh + 343);
    const auto *hh_344 = buffer.data(hh + 344);
    const auto *hh_345 = buffer.data(hh + 345);
    const auto *hh_346 = buffer.data(hh + 346);
    const auto *hh_347 = buffer.data(hh + 347);
    const auto *hh_348 = buffer.data(hh + 348);
    const auto *hh_349 = buffer.data(hh + 349);
    const auto *hh_350 = buffer.data(hh + 350);
    const auto *hh_351 = buffer.data(hh + 351);
    const auto *hh_352 = buffer.data(hh + 352);
    const auto *hh_353 = buffer.data(hh + 353);
    const auto *hh_354 = buffer.data(hh + 354);
    const auto *hh_355 = buffer.data(hh + 355);
    const auto *hh_356 = buffer.data(hh + 356);
    const auto *hh_357 = buffer.data(hh + 357);
    const auto *hh_358 = buffer.data(hh + 358);
    const auto *hh_359 = buffer.data(hh + 359);
    const auto *hh_360 = buffer.data(hh + 360);
    const auto *hh_361 = buffer.data(hh + 361);
    const auto *hh_362 = buffer.data(hh + 362);
    const auto *hh_363 = buffer.data(hh + 363);
    const auto *hh_364 = buffer.data(hh + 364);
    const auto *hh_365 = buffer.data(hh + 365);
    const auto *hh_366 = buffer.data(hh + 366);
    const auto *hh_367 = buffer.data(hh + 367);
    const auto *hh_368 = buffer.data(hh + 368);
    const auto *hh_369 = buffer.data(hh + 369);
    const auto *hh_370 = buffer.data(hh + 370);
    const auto *hh_371 = buffer.data(hh + 371);
    const auto *hh_372 = buffer.data(hh + 372);
    const auto *hh_373 = buffer.data(hh + 373);
    const auto *hh_374 = buffer.data(hh + 374);
    const auto *hh_375 = buffer.data(hh + 375);
    const auto *hh_376 = buffer.data(hh + 376);
    const auto *hh_377 = buffer.data(hh + 377);
    const auto *hh_378 = buffer.data(hh + 378);
    const auto *hh_379 = buffer.data(hh + 379);
    const auto *hh_380 = buffer.data(hh + 380);
    const auto *hh_381 = buffer.data(hh + 381);
    const auto *hh_382 = buffer.data(hh + 382);
    const auto *hh_383 = buffer.data(hh + 383);
    const auto *hh_384 = buffer.data(hh + 384);
    const auto *hh_385 = buffer.data(hh + 385);
    const auto *hh_386 = buffer.data(hh + 386);
    const auto *hh_387 = buffer.data(hh + 387);
    const auto *hh_388 = buffer.data(hh + 388);
    const auto *hh_389 = buffer.data(hh + 389);
    const auto *hh_390 = buffer.data(hh + 390);
    const auto *hh_391 = buffer.data(hh + 391);
    const auto *hh_392 = buffer.data(hh + 392);
    const auto *hh_393 = buffer.data(hh + 393);
    const auto *hh_394 = buffer.data(hh + 394);
    const auto *hh_395 = buffer.data(hh + 395);
    const auto *hh_396 = buffer.data(hh + 396);
    const auto *hh_397 = buffer.data(hh + 397);
    const auto *hh_398 = buffer.data(hh + 398);
    const auto *hh_399 = buffer.data(hh + 399);
    const auto *hh_400 = buffer.data(hh + 400);
    const auto *hh_401 = buffer.data(hh + 401);
    const auto *hh_402 = buffer.data(hh + 402);
    const auto *hh_403 = buffer.data(hh + 403);
    const auto *hh_404 = buffer.data(hh + 404);
    const auto *hh_405 = buffer.data(hh + 405);
    const auto *hh_406 = buffer.data(hh + 406);
    const auto *hh_407 = buffer.data(hh + 407);
    const auto *hh_408 = buffer.data(hh + 408);
    const auto *hh_409 = buffer.data(hh + 409);
    const auto *hh_410 = buffer.data(hh + 410);
    const auto *hh_411 = buffer.data(hh + 411);
    const auto *hh_412 = buffer.data(hh + 412);
    const auto *hh_413 = buffer.data(hh + 413);
    const auto *hh_414 = buffer.data(hh + 414);
    const auto *hh_415 = buffer.data(hh + 415);
    const auto *hh_416 = buffer.data(hh + 416);
    const auto *hh_417 = buffer.data(hh + 417);
    const auto *hh_418 = buffer.data(hh + 418);
    const auto *hh_419 = buffer.data(hh + 419);
    const auto *hh_420 = buffer.data(hh + 420);
    const auto *hh_421 = buffer.data(hh + 421);
    const auto *hh_422 = buffer.data(hh + 422);
    const auto *hh_423 = buffer.data(hh + 423);
    const auto *hh_424 = buffer.data(hh + 424);
    const auto *hh_425 = buffer.data(hh + 425);
    const auto *hh_426 = buffer.data(hh + 426);
    const auto *hh_427 = buffer.data(hh + 427);
    const auto *hh_428 = buffer.data(hh + 428);
    const auto *hh_429 = buffer.data(hh + 429);
    const auto *hh_430 = buffer.data(hh + 430);
    const auto *hh_431 = buffer.data(hh + 431);
    const auto *hh_432 = buffer.data(hh + 432);
    const auto *hh_433 = buffer.data(hh + 433);
    const auto *hh_434 = buffer.data(hh + 434);
    const auto *hh_435 = buffer.data(hh + 435);
    const auto *hh_436 = buffer.data(hh + 436);
    const auto *hh_437 = buffer.data(hh + 437);
    const auto *hh_438 = buffer.data(hh + 438);
    const auto *hh_439 = buffer.data(hh + 439);
    const auto *hh_440 = buffer.data(hh + 440);

#pragma omp simd aligned(hh_22, hh_27, hh_36, hh_127, hh_132, hh_141, hh_316, hh_321, \
                         hh_330 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = 12.3046875 * hh_22[k]
                 - 24.609375 * hh_27[k]
                 + 2.4609375 * hh_36[k]
                 - 24.609375 * hh_127[k]
                 + 49.21875 * hh_132[k]
                 - 4.921875 * hh_141[k]
                 + 2.4609375 * hh_316[k]
                 - 4.921875 * hh_321[k]
                 + 0.4921875 * hh_330[k];
    }

#pragma omp simd aligned(hh_25, hh_32, hh_130, hh_137, hh_319, hh_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = f_0 * hh_25[k]
                 - f_0 * hh_32[k]
                 - f_1 * hh_130[k]
                 + f_1 * hh_137[k]
                 + f_2 * hh_319[k]
                 - f_2 * hh_326[k];
    }

#pragma omp simd aligned(hh_22, hh_27, hh_29, hh_36, hh_38, hh_127, hh_132, hh_134, hh_141, \
                         hh_143, hh_316, hh_321, hh_323, hh_330, \
                         hh_332 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_3 * hh_22[k]
                 - f_4 * hh_27[k]
                 + f_5 * hh_29[k]
                 + f_6 * hh_36[k]
                 - f_7 * hh_38[k]
                 + f_8 * hh_127[k]
                 + f_9 * hh_132[k]
                 - f_10 * hh_134[k]
                 - f_4 * hh_141[k]
                 + f_11 * hh_143[k]
                 - f_12 * hh_316[k]
                 - f_13 * hh_321[k]
                 + f_14 * hh_323[k]
                 + f_15 * hh_330[k]
                 - f_16 * hh_332[k];
    }

#pragma omp simd aligned(hh_25, hh_32, hh_34, hh_130, hh_137, hh_139, hh_319, hh_326, \
                         hh_328 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_17 * hh_25[k]
                 - f_17 * hh_32[k]
                 + f_18 * hh_34[k]
                 + f_18 * hh_130[k]
                 + f_18 * hh_137[k]
                 - f_19 * hh_139[k]
                 - f_20 * hh_319[k]
                 - f_20 * hh_326[k]
                 + f_21 * hh_328[k];
    }

#pragma omp simd aligned(hh_22, hh_27, hh_29, hh_36, hh_38, hh_40, hh_127, hh_132, hh_134, \
                         hh_141, hh_143, hh_145, hh_316, hh_321, hh_323, hh_330, hh_332, \
                         hh_334 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_22 * hh_22[k]
                 + f_23 * hh_27[k]
                 - f_24 * hh_29[k]
                 + f_22 * hh_36[k]
                 - f_24 * hh_38[k]
                 + f_25 * hh_40[k]
                 - f_23 * hh_127[k]
                 - f_26 * hh_132[k]
                 + f_27 * hh_134[k]
                 - f_23 * hh_141[k]
                 + f_27 * hh_143[k]
                 - f_28 * hh_145[k]
                 + f_29 * hh_316[k]
                 + f_30 * hh_321[k]
                 - f_31 * hh_323[k]
                 + f_29 * hh_330[k]
                 - f_31 * hh_332[k]
                 + f_32 * hh_334[k];
    }

#pragma omp simd aligned(hh_23, hh_28, hh_30, hh_37, hh_39, hh_41, hh_128, hh_133, hh_135, \
                         hh_142, hh_144, hh_146, hh_317, hh_322, hh_324, hh_331, hh_333, \
                         hh_335 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_33 * hh_23[k]
                 + f_34 * hh_28[k]
                 - f_35 * hh_30[k]
                 + f_33 * hh_37[k]
                 - f_35 * hh_39[k]
                 + f_36 * hh_41[k]
                 - f_34 * hh_128[k]
                 - f_37 * hh_133[k]
                 + f_38 * hh_135[k]
                 - f_34 * hh_142[k]
                 + f_38 * hh_144[k]
                 - f_39 * hh_146[k]
                 + f_40 * hh_317[k]
                 + f_41 * hh_322[k]
                 - f_36 * hh_324[k]
                 + f_40 * hh_331[k]
                 - f_36 * hh_333[k]
                 + f_42 * hh_335[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_26, hh_31, hh_33, hh_35, hh_126, hh_129, hh_131, \
                         hh_136, hh_138, hh_140, hh_315, hh_318, hh_320, hh_325, hh_327, \
                         hh_329 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_22 * hh_21[k]
                 + f_23 * hh_24[k]
                 - f_24 * hh_26[k]
                 + f_22 * hh_31[k]
                 - f_24 * hh_33[k]
                 + f_25 * hh_35[k]
                 - f_23 * hh_126[k]
                 - f_26 * hh_129[k]
                 + f_27 * hh_131[k]
                 - f_23 * hh_136[k]
                 + f_27 * hh_138[k]
                 - f_28 * hh_140[k]
                 + f_29 * hh_315[k]
                 + f_30 * hh_318[k]
                 - f_31 * hh_320[k]
                 + f_29 * hh_325[k]
                 - f_31 * hh_327[k]
                 + f_32 * hh_329[k];
    }

#pragma omp simd aligned(hh_23, hh_30, hh_37, hh_39, hh_128, hh_135, hh_142, hh_144, hh_317, \
                         hh_324, hh_331, hh_333 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_43 * hh_23[k]
                 + f_17 * hh_30[k]
                 + f_43 * hh_37[k]
                 - f_17 * hh_39[k]
                 + f_17 * hh_128[k]
                 - f_18 * hh_135[k]
                 - f_17 * hh_142[k]
                 + f_18 * hh_144[k]
                 - f_44 * hh_317[k]
                 + f_20 * hh_324[k]
                 + f_44 * hh_331[k]
                 - f_20 * hh_333[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_26, hh_31, hh_33, hh_126, hh_129, hh_131, hh_136, \
                         hh_138, hh_315, hh_318, hh_320, hh_325, \
                         hh_327 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_6 * hh_21[k]
                 + f_4 * hh_24[k]
                 + f_7 * hh_26[k]
                 + f_3 * hh_31[k]
                 - f_5 * hh_33[k]
                 + f_4 * hh_126[k]
                 - f_9 * hh_129[k]
                 - f_11 * hh_131[k]
                 - f_8 * hh_136[k]
                 + f_10 * hh_138[k]
                 - f_15 * hh_315[k]
                 + f_13 * hh_318[k]
                 + f_16 * hh_320[k]
                 + f_12 * hh_325[k]
                 - f_14 * hh_327[k];
    }

#pragma omp simd aligned(hh_23, hh_28, hh_37, hh_128, hh_133, hh_142, hh_317, hh_322, \
                         hh_331 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_45 * hh_23[k]
                 - f_46 * hh_28[k]
                 + f_45 * hh_37[k]
                 - f_47 * hh_128[k]
                 + f_48 * hh_133[k]
                 - f_47 * hh_142[k]
                 + f_49 * hh_317[k]
                 - f_50 * hh_322[k]
                 + f_49 * hh_331[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_31, hh_126, hh_129, hh_136, hh_315, hh_318, \
                         hh_325 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = 2.4609375 * hh_21[k]
                  - 24.609375 * hh_24[k]
                  + 12.3046875 * hh_31[k]
                  - 4.921875 * hh_126[k]
                  + 49.21875 * hh_129[k]
                  - 24.609375 * hh_136[k]
                  + 0.4921875 * hh_315[k]
                  - 4.921875 * hh_318[k]
                  + 2.4609375 * hh_325[k];
    }

#pragma omp simd aligned(hh_85, hh_88, hh_90, hh_95, hh_99, hh_232, hh_235, hh_237, hh_242, \
                         hh_246 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_0 * hh_85[k]
                  - f_1 * hh_90[k]
                  + f_2 * hh_99[k]
                  - f_0 * hh_232[k]
                  + f_1 * hh_237[k]
                  - f_2 * hh_246[k];

        g_12[k] = 78.75 * hh_88[k]
                  - 78.75 * hh_95[k]
                  - 78.75 * hh_235[k]
                  + 78.75 * hh_242[k];
    }

#pragma omp simd aligned(hh_85, hh_90, hh_92, hh_99, hh_101, hh_232, hh_237, hh_239, hh_246, \
                         hh_248 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_51 * hh_85[k]
                  - f_52 * hh_90[k]
                  + f_53 * hh_92[k]
                  + f_54 * hh_99[k]
                  - f_55 * hh_101[k]
                  + f_51 * hh_232[k]
                  + f_52 * hh_237[k]
                  - f_53 * hh_239[k]
                  - f_54 * hh_246[k]
                  + f_55 * hh_248[k];
    }

#pragma omp simd aligned(hh_88, hh_95, hh_97, hh_235, hh_242, hh_244 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_56 * hh_88[k]
                  - f_56 * hh_95[k]
                  + f_57 * hh_97[k]
                  + f_56 * hh_235[k]
                  + f_56 * hh_242[k]
                  - f_57 * hh_244[k];
    }

#pragma omp simd aligned(hh_85, hh_90, hh_92, hh_99, hh_101, hh_103, hh_232, hh_237, hh_239, \
                         hh_246, hh_248, hh_250 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_58 * hh_85[k]
                  + f_59 * hh_90[k]
                  - f_60 * hh_92[k]
                  + f_58 * hh_99[k]
                  - f_60 * hh_101[k]
                  + f_61 * hh_103[k]
                  - f_58 * hh_232[k]
                  - f_59 * hh_237[k]
                  + f_60 * hh_239[k]
                  - f_58 * hh_246[k]
                  + f_60 * hh_248[k]
                  - f_61 * hh_250[k];
    }

#pragma omp simd aligned(hh_86, hh_91, hh_93, hh_100, hh_102, hh_104, hh_233, hh_238, hh_240, \
                         hh_247, hh_249, hh_251 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_62 * hh_86[k]
                  + f_63 * hh_91[k]
                  - f_64 * hh_93[k]
                  + f_62 * hh_100[k]
                  - f_64 * hh_102[k]
                  + f_65 * hh_104[k]
                  - f_62 * hh_233[k]
                  - f_63 * hh_238[k]
                  + f_64 * hh_240[k]
                  - f_62 * hh_247[k]
                  + f_64 * hh_249[k]
                  - f_65 * hh_251[k];
    }

#pragma omp simd aligned(hh_84, hh_87, hh_89, hh_94, hh_96, hh_98, hh_231, hh_234, hh_236, \
                         hh_241, hh_243, hh_245 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_58 * hh_84[k]
                  + f_59 * hh_87[k]
                  - f_60 * hh_89[k]
                  + f_58 * hh_94[k]
                  - f_60 * hh_96[k]
                  + f_61 * hh_98[k]
                  - f_58 * hh_231[k]
                  - f_59 * hh_234[k]
                  + f_60 * hh_236[k]
                  - f_58 * hh_241[k]
                  + f_60 * hh_243[k]
                  - f_61 * hh_245[k];
    }

#pragma omp simd aligned(hh_86, hh_93, hh_100, hh_102, hh_233, hh_240, hh_247, \
                         hh_249 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_66 * hh_86[k]
                  + f_56 * hh_93[k]
                  + f_66 * hh_100[k]
                  - f_56 * hh_102[k]
                  + f_66 * hh_233[k]
                  - f_56 * hh_240[k]
                  - f_66 * hh_247[k]
                  + f_56 * hh_249[k];
    }

#pragma omp simd aligned(hh_84, hh_87, hh_89, hh_94, hh_96, hh_231, hh_234, hh_236, hh_241, \
                         hh_243 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_54 * hh_84[k]
                  + f_52 * hh_87[k]
                  + f_55 * hh_89[k]
                  + f_51 * hh_94[k]
                  - f_53 * hh_96[k]
                  + f_54 * hh_231[k]
                  - f_52 * hh_234[k]
                  - f_55 * hh_236[k]
                  - f_51 * hh_241[k]
                  + f_53 * hh_243[k];
    }

#pragma omp simd aligned(hh_84, hh_86, hh_87, hh_91, hh_94, hh_100, hh_231, hh_233, hh_234, \
                         hh_238, hh_241, hh_247 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = 19.6875 * hh_86[k]
                  - 118.125 * hh_91[k]
                  + 19.6875 * hh_100[k]
                  - 19.6875 * hh_233[k]
                  + 118.125 * hh_238[k]
                  - 19.6875 * hh_247[k];

        g_21[k] = f_2 * hh_84[k]
                  - f_1 * hh_87[k]
                  + f_0 * hh_94[k]
                  - f_2 * hh_231[k]
                  + f_1 * hh_234[k]
                  - f_0 * hh_241[k];
    }

#pragma omp simd aligned(hh_22, hh_27, hh_36, hh_127, hh_132, hh_141, hh_169, hh_174, hh_183, \
                         hh_316, hh_321, hh_330, hh_358, hh_363, \
                         hh_372 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_3 * hh_22[k]
                  + f_8 * hh_27[k]
                  - f_12 * hh_36[k]
                  - f_4 * hh_127[k]
                  + f_9 * hh_132[k]
                  - f_13 * hh_141[k]
                  + f_5 * hh_169[k]
                  - f_10 * hh_174[k]
                  + f_14 * hh_183[k]
                  + f_6 * hh_316[k]
                  - f_4 * hh_321[k]
                  + f_15 * hh_330[k]
                  - f_7 * hh_358[k]
                  + f_11 * hh_363[k]
                  - f_16 * hh_372[k];
    }

#pragma omp simd aligned(hh_25, hh_32, hh_130, hh_137, hh_172, hh_179, hh_319, hh_326, hh_361, \
                         hh_368 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_51 * hh_25[k]
                  + f_51 * hh_32[k]
                  - f_52 * hh_130[k]
                  + f_52 * hh_137[k]
                  + f_53 * hh_172[k]
                  - f_53 * hh_179[k]
                  + f_54 * hh_319[k]
                  - f_54 * hh_326[k]
                  - f_55 * hh_361[k]
                  + f_55 * hh_368[k];
    }

#pragma omp simd aligned(hh_22, hh_27, hh_29, hh_36, hh_38, hh_127, hh_132, hh_134, hh_141, \
                         hh_143, hh_169, hh_174, hh_176, hh_183, hh_185, hh_316, hh_321, \
                         hh_323, hh_330, hh_332, hh_358, hh_363, hh_365, hh_372, \
                         hh_374 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = 2.4609375 * hh_22[k]
                  + 1.640625 * hh_27[k]
                  - 19.6875 * hh_29[k]
                  - 0.8203125 * hh_36[k]
                  + 6.5625 * hh_38[k]
                  + 1.640625 * hh_127[k]
                  + 1.09375 * hh_132[k]
                  - 13.125 * hh_134[k]
                  - 0.546875 * hh_141[k]
                  + 4.375 * hh_143[k]
                  - 19.6875 * hh_169[k]
                  - 13.125 * hh_174[k]
                  + 157.5 * hh_176[k]
                  + 6.5625 * hh_183[k]
                  - 52.5 * hh_185[k]
                  - 0.8203125 * hh_316[k]
                  - 0.546875 * hh_321[k]
                  + 6.5625 * hh_323[k]
                  + 0.2734375 * hh_330[k]
                  - 2.1875 * hh_332[k]
                  + 6.5625 * hh_358[k]
                  + 4.375 * hh_363[k]
                  - 52.5 * hh_365[k]
                  - 2.1875 * hh_372[k]
                  + 17.5 * hh_374[k];
    }

#pragma omp simd aligned(hh_25, hh_32, hh_34, hh_130, hh_137, hh_139, hh_172, hh_179, hh_181, \
                         hh_319, hh_326, hh_328, hh_361, hh_368, \
                         hh_370 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_67 * hh_25[k]
                  + f_67 * hh_32[k]
                  - f_68 * hh_34[k]
                  + f_69 * hh_130[k]
                  + f_69 * hh_137[k]
                  - f_70 * hh_139[k]
                  - f_71 * hh_172[k]
                  - f_71 * hh_179[k]
                  + f_72 * hh_181[k]
                  - f_73 * hh_319[k]
                  - f_73 * hh_326[k]
                  + f_69 * hh_328[k]
                  + f_74 * hh_361[k]
                  + f_74 * hh_368[k]
                  - f_75 * hh_370[k];
    }

#pragma omp simd aligned(hh_22, hh_27, hh_29, hh_36, hh_38, hh_40, hh_127, hh_132, hh_134, \
                         hh_141, hh_143, hh_145, hh_169, hh_174, hh_176, hh_183, hh_185, \
                         hh_187, hh_316, hh_321, hh_323, hh_330, hh_332, hh_334, hh_358, \
                         hh_363, hh_365, hh_372, hh_374, hh_376 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_76 * hh_22[k]
                  - f_77 * hh_27[k]
                  + f_78 * hh_29[k]
                  - f_76 * hh_36[k]
                  + f_78 * hh_38[k]
                  - f_79 * hh_40[k]
                  - f_80 * hh_127[k]
                  - f_81 * hh_132[k]
                  + f_79 * hh_134[k]
                  - f_80 * hh_141[k]
                  + f_79 * hh_143[k]
                  - f_82 * hh_145[k]
                  + f_79 * hh_169[k]
                  + f_83 * hh_174[k]
                  - f_84 * hh_176[k]
                  + f_79 * hh_183[k]
                  - f_84 * hh_185[k]
                  + f_85 * hh_187[k]
                  + f_86 * hh_316[k]
                  + f_80 * hh_321[k]
                  - f_87 * hh_323[k]
                  + f_86 * hh_330[k]
                  - f_87 * hh_332[k]
                  + f_88 * hh_334[k]
                  - f_88 * hh_358[k]
                  - f_82 * hh_363[k]
                  + f_89 * hh_365[k]
                  - f_88 * hh_372[k]
                  + f_89 * hh_374[k]
                  - f_90 * hh_376[k];
    }

#pragma omp simd aligned(hh_23, hh_28, hh_30, hh_37, hh_39, hh_41, hh_128, hh_133, hh_135, \
                         hh_142, hh_144, hh_146, hh_170, hh_175, hh_177, hh_184, hh_186, \
                         hh_188, hh_317, hh_322, hh_324, hh_331, hh_333, hh_335, hh_359, \
                         hh_364, hh_366, hh_373, hh_375, hh_377 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_91 * hh_23[k]
                  - f_92 * hh_28[k]
                  + f_93 * hh_30[k]
                  - f_91 * hh_37[k]
                  + f_93 * hh_39[k]
                  - f_94 * hh_41[k]
                  - f_95 * hh_128[k]
                  - f_96 * hh_133[k]
                  + f_97 * hh_135[k]
                  - f_95 * hh_142[k]
                  + f_97 * hh_144[k]
                  - f_98 * hh_146[k]
                  + f_99 * hh_170[k]
                  + f_100 * hh_175[k]
                  - f_101 * hh_177[k]
                  + f_99 * hh_184[k]
                  - f_101 * hh_186[k]
                  + f_102 * hh_188[k]
                  + f_103 * hh_317[k]
                  + f_95 * hh_322[k]
                  - f_104 * hh_324[k]
                  + f_103 * hh_331[k]
                  - f_104 * hh_333[k]
                  + f_105 * hh_335[k]
                  - f_93 * hh_359[k]
                  - f_106 * hh_364[k]
                  + f_107 * hh_366[k]
                  - f_93 * hh_373[k]
                  + f_107 * hh_375[k]
                  - f_108 * hh_377[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_26, hh_31, hh_33, hh_35, hh_126, hh_129, hh_131, \
                         hh_136, hh_138, hh_140, hh_168, hh_171, hh_173, hh_178, hh_180, \
                         hh_182, hh_315, hh_318, hh_320, hh_325, hh_327, hh_329, hh_357, \
                         hh_360, hh_362, hh_367, hh_369, hh_371 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_76 * hh_21[k]
                  - f_77 * hh_24[k]
                  + f_78 * hh_26[k]
                  - f_76 * hh_31[k]
                  + f_78 * hh_33[k]
                  - f_79 * hh_35[k]
                  - f_80 * hh_126[k]
                  - f_81 * hh_129[k]
                  + f_79 * hh_131[k]
                  - f_80 * hh_136[k]
                  + f_79 * hh_138[k]
                  - f_82 * hh_140[k]
                  + f_79 * hh_168[k]
                  + f_83 * hh_171[k]
                  - f_84 * hh_173[k]
                  + f_79 * hh_178[k]
                  - f_84 * hh_180[k]
                  + f_85 * hh_182[k]
                  + f_86 * hh_315[k]
                  + f_80 * hh_318[k]
                  - f_87 * hh_320[k]
                  + f_86 * hh_325[k]
                  - f_87 * hh_327[k]
                  + f_88 * hh_329[k]
                  - f_88 * hh_357[k]
                  - f_82 * hh_360[k]
                  + f_89 * hh_362[k]
                  - f_88 * hh_367[k]
                  + f_89 * hh_369[k]
                  - f_90 * hh_371[k];
    }

#pragma omp simd aligned(hh_23, hh_30, hh_37, hh_39, hh_128, hh_135, hh_142, hh_144, hh_170, \
                         hh_177, hh_184, hh_186, hh_317, hh_324, hh_331, hh_333, hh_359, \
                         hh_366, hh_373, hh_375 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_109 * hh_23[k]
                  - f_67 * hh_30[k]
                  - f_109 * hh_37[k]
                  + f_67 * hh_39[k]
                  + f_73 * hh_128[k]
                  - f_69 * hh_135[k]
                  - f_73 * hh_142[k]
                  + f_69 * hh_144[k]
                  - f_110 * hh_170[k]
                  + f_71 * hh_177[k]
                  + f_110 * hh_184[k]
                  - f_71 * hh_186[k]
                  - f_111 * hh_317[k]
                  + f_73 * hh_324[k]
                  + f_111 * hh_331[k]
                  - f_73 * hh_333[k]
                  + f_70 * hh_359[k]
                  - f_74 * hh_366[k]
                  - f_70 * hh_373[k]
                  + f_74 * hh_375[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_26, hh_31, hh_33, hh_126, hh_129, hh_131, hh_136, \
                         hh_138, hh_168, hh_171, hh_173, hh_178, hh_180, hh_315, hh_318, \
                         hh_320, hh_325, hh_327, hh_357, hh_360, hh_362, hh_367, \
                         hh_369 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = 0.8203125 * hh_21[k]
                  - 1.640625 * hh_24[k]
                  - 6.5625 * hh_26[k]
                  - 2.4609375 * hh_31[k]
                  + 19.6875 * hh_33[k]
                  + 0.546875 * hh_126[k]
                  - 1.09375 * hh_129[k]
                  - 4.375 * hh_131[k]
                  - 1.640625 * hh_136[k]
                  + 13.125 * hh_138[k]
                  - 6.5625 * hh_168[k]
                  + 13.125 * hh_171[k]
                  + 52.5 * hh_173[k]
                  + 19.6875 * hh_178[k]
                  - 157.5 * hh_180[k]
                  - 0.2734375 * hh_315[k]
                  + 0.546875 * hh_318[k]
                  + 2.1875 * hh_320[k]
                  + 0.8203125 * hh_325[k]
                  - 6.5625 * hh_327[k]
                  + 2.1875 * hh_357[k]
                  - 4.375 * hh_360[k]
                  - 17.5 * hh_362[k]
                  - 6.5625 * hh_367[k]
                  + 52.5 * hh_369[k];
    }

#pragma omp simd aligned(hh_23, hh_28, hh_37, hh_128, hh_133, hh_142, hh_170, hh_175, hh_184, \
                         hh_317, hh_322, hh_331, hh_359, hh_364, \
                         hh_373 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_112 * hh_23[k]
                  + f_113 * hh_28[k]
                  - f_112 * hh_37[k]
                  - f_114 * hh_128[k]
                  + f_51 * hh_133[k]
                  - f_114 * hh_142[k]
                  + f_115 * hh_170[k]
                  - f_116 * hh_175[k]
                  + f_115 * hh_184[k]
                  + f_117 * hh_317[k]
                  - f_118 * hh_322[k]
                  + f_117 * hh_331[k]
                  - f_52 * hh_359[k]
                  + f_119 * hh_364[k]
                  - f_52 * hh_373[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_31, hh_126, hh_129, hh_136, hh_168, hh_171, hh_178, \
                         hh_315, hh_318, hh_325, hh_357, hh_360, \
                         hh_367 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_12 * hh_21[k]
                  + f_8 * hh_24[k]
                  - f_3 * hh_31[k]
                  - f_13 * hh_126[k]
                  + f_9 * hh_129[k]
                  - f_4 * hh_136[k]
                  + f_14 * hh_168[k]
                  - f_10 * hh_171[k]
                  + f_5 * hh_178[k]
                  + f_15 * hh_315[k]
                  - f_4 * hh_318[k]
                  + f_6 * hh_325[k]
                  - f_16 * hh_357[k]
                  + f_11 * hh_360[k]
                  - f_7 * hh_367[k];
    }

#pragma omp simd aligned(hh_85, hh_90, hh_99, hh_232, hh_237, hh_246, hh_274, hh_279, \
                         hh_288 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_17 * hh_85[k]
                  + f_18 * hh_90[k]
                  - f_20 * hh_99[k]
                  - f_17 * hh_232[k]
                  + f_18 * hh_237[k]
                  - f_20 * hh_246[k]
                  + f_18 * hh_274[k]
                  - f_19 * hh_279[k]
                  + f_21 * hh_288[k];
    }

#pragma omp simd aligned(hh_88, hh_95, hh_235, hh_242, hh_277, hh_284 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_56 * hh_88[k]
                  + f_56 * hh_95[k]
                  - f_56 * hh_235[k]
                  + f_56 * hh_242[k]
                  + f_57 * hh_277[k]
                  - f_57 * hh_284[k];
    }

#pragma omp simd aligned(hh_85, hh_90, hh_92, hh_99, hh_101, hh_232, hh_237, hh_239, hh_246, \
                         hh_248, hh_274, hh_279, hh_281, hh_288, \
                         hh_290 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_67 * hh_85[k]
                  + f_69 * hh_90[k]
                  - f_71 * hh_92[k]
                  - f_73 * hh_99[k]
                  + f_74 * hh_101[k]
                  + f_67 * hh_232[k]
                  + f_69 * hh_237[k]
                  - f_71 * hh_239[k]
                  - f_73 * hh_246[k]
                  + f_74 * hh_248[k]
                  - f_68 * hh_274[k]
                  - f_70 * hh_279[k]
                  + f_72 * hh_281[k]
                  + f_69 * hh_288[k]
                  - f_75 * hh_290[k];
    }

#pragma omp simd aligned(hh_88, hh_95, hh_97, hh_235, hh_242, hh_244, hh_277, hh_284, \
                         hh_286 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = 26.25 * hh_88[k]
                  + 26.25 * hh_95[k]
                  - 52.5 * hh_97[k]
                  + 26.25 * hh_235[k]
                  + 26.25 * hh_242[k]
                  - 52.5 * hh_244[k]
                  - 52.5 * hh_277[k]
                  - 52.5 * hh_284[k]
                  + 105.0 * hh_286[k];
    }

#pragma omp simd aligned(hh_85, hh_90, hh_92, hh_99, hh_101, hh_103, hh_232, hh_237, hh_239, \
                         hh_246, hh_248, hh_250, hh_274, hh_279, hh_281, hh_288, hh_290, \
                         hh_292 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_120 * hh_85[k]
                  - f_121 * hh_90[k]
                  + f_122 * hh_92[k]
                  - f_120 * hh_99[k]
                  + f_122 * hh_101[k]
                  - f_123 * hh_103[k]
                  - f_120 * hh_232[k]
                  - f_121 * hh_237[k]
                  + f_122 * hh_239[k]
                  - f_120 * hh_246[k]
                  + f_122 * hh_248[k]
                  - f_123 * hh_250[k]
                  + f_121 * hh_274[k]
                  + f_124 * hh_279[k]
                  - f_125 * hh_281[k]
                  + f_121 * hh_288[k]
                  - f_125 * hh_290[k]
                  + f_126 * hh_292[k];
    }

#pragma omp simd aligned(hh_86, hh_91, hh_93, hh_100, hh_102, hh_104, hh_233, hh_238, hh_240, \
                         hh_247, hh_249, hh_251, hh_275, hh_280, hh_282, hh_289, hh_291, \
                         hh_293 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_127 * hh_86[k]
                  - f_128 * hh_91[k]
                  + f_129 * hh_93[k]
                  - f_127 * hh_100[k]
                  + f_129 * hh_102[k]
                  - f_130 * hh_104[k]
                  - f_127 * hh_233[k]
                  - f_128 * hh_238[k]
                  + f_129 * hh_240[k]
                  - f_127 * hh_247[k]
                  + f_129 * hh_249[k]
                  - f_130 * hh_251[k]
                  + f_128 * hh_275[k]
                  + f_131 * hh_280[k]
                  - f_132 * hh_282[k]
                  + f_128 * hh_289[k]
                  - f_132 * hh_291[k]
                  + f_133 * hh_293[k];
    }

#pragma omp simd aligned(hh_84, hh_87, hh_89, hh_94, hh_96, hh_98, hh_231, hh_234, hh_236, \
                         hh_241, hh_243, hh_245, hh_273, hh_276, hh_278, hh_283, hh_285, \
                         hh_287 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_120 * hh_84[k]
                  - f_121 * hh_87[k]
                  + f_122 * hh_89[k]
                  - f_120 * hh_94[k]
                  + f_122 * hh_96[k]
                  - f_123 * hh_98[k]
                  - f_120 * hh_231[k]
                  - f_121 * hh_234[k]
                  + f_122 * hh_236[k]
                  - f_120 * hh_241[k]
                  + f_122 * hh_243[k]
                  - f_123 * hh_245[k]
                  + f_121 * hh_273[k]
                  + f_124 * hh_276[k]
                  - f_125 * hh_278[k]
                  + f_121 * hh_283[k]
                  - f_125 * hh_285[k]
                  + f_126 * hh_287[k];
    }

#pragma omp simd aligned(hh_86, hh_93, hh_100, hh_102, hh_233, hh_240, hh_247, hh_249, hh_275, \
                         hh_282, hh_289, hh_291 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = 13.125 * hh_86[k]
                  - 26.25 * hh_93[k]
                  - 13.125 * hh_100[k]
                  + 26.25 * hh_102[k]
                  + 13.125 * hh_233[k]
                  - 26.25 * hh_240[k]
                  - 13.125 * hh_247[k]
                  + 26.25 * hh_249[k]
                  - 26.25 * hh_275[k]
                  + 52.5 * hh_282[k]
                  + 26.25 * hh_289[k]
                  - 52.5 * hh_291[k];
    }

#pragma omp simd aligned(hh_84, hh_87, hh_89, hh_94, hh_96, hh_231, hh_234, hh_236, hh_241, \
                         hh_243, hh_273, hh_276, hh_278, hh_283, \
                         hh_285 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_73 * hh_84[k]
                  - f_69 * hh_87[k]
                  - f_74 * hh_89[k]
                  - f_67 * hh_94[k]
                  + f_71 * hh_96[k]
                  + f_73 * hh_231[k]
                  - f_69 * hh_234[k]
                  - f_74 * hh_236[k]
                  - f_67 * hh_241[k]
                  + f_71 * hh_243[k]
                  - f_69 * hh_273[k]
                  + f_70 * hh_276[k]
                  + f_75 * hh_278[k]
                  + f_68 * hh_283[k]
                  - f_72 * hh_285[k];
    }

#pragma omp simd aligned(hh_86, hh_91, hh_100, hh_233, hh_238, hh_247, hh_275, hh_280, \
                         hh_289 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_134 * hh_86[k]
                  + f_135 * hh_91[k]
                  - f_134 * hh_100[k]
                  - f_134 * hh_233[k]
                  + f_135 * hh_238[k]
                  - f_134 * hh_247[k]
                  + f_66 * hh_275[k]
                  - f_136 * hh_280[k]
                  + f_66 * hh_289[k];
    }

#pragma omp simd aligned(hh_84, hh_87, hh_94, hh_231, hh_234, hh_241, hh_273, hh_276, \
                         hh_283 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_20 * hh_84[k]
                  + f_18 * hh_87[k]
                  - f_17 * hh_94[k]
                  - f_20 * hh_231[k]
                  + f_18 * hh_234[k]
                  - f_17 * hh_241[k]
                  + f_21 * hh_273[k]
                  - f_19 * hh_276[k]
                  + f_18 * hh_283[k];
    }

#pragma omp simd aligned(hh_22, hh_27, hh_36, hh_127, hh_132, hh_141, hh_169, hh_174, hh_183, \
                         hh_316, hh_321, hh_330, hh_358, hh_363, hh_372, hh_400, hh_405, \
                         hh_414 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_22 * hh_22[k]
                  - f_23 * hh_27[k]
                  + f_29 * hh_36[k]
                  + f_23 * hh_127[k]
                  - f_26 * hh_132[k]
                  + f_30 * hh_141[k]
                  - f_24 * hh_169[k]
                  + f_27 * hh_174[k]
                  - f_31 * hh_183[k]
                  + f_22 * hh_316[k]
                  - f_23 * hh_321[k]
                  + f_29 * hh_330[k]
                  - f_24 * hh_358[k]
                  + f_27 * hh_363[k]
                  - f_31 * hh_372[k]
                  + f_25 * hh_400[k]
                  - f_28 * hh_405[k]
                  + f_32 * hh_414[k];
    }

#pragma omp simd aligned(hh_25, hh_32, hh_130, hh_137, hh_172, hh_179, hh_319, hh_326, hh_361, \
                         hh_368, hh_403, hh_410 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_58 * hh_25[k]
                  - f_58 * hh_32[k]
                  + f_59 * hh_130[k]
                  - f_59 * hh_137[k]
                  - f_60 * hh_172[k]
                  + f_60 * hh_179[k]
                  + f_58 * hh_319[k]
                  - f_58 * hh_326[k]
                  - f_60 * hh_361[k]
                  + f_60 * hh_368[k]
                  + f_61 * hh_403[k]
                  - f_61 * hh_410[k];
    }

#pragma omp simd aligned(hh_22, hh_27, hh_29, hh_36, hh_38, hh_127, hh_132, hh_134, hh_141, \
                         hh_143, hh_169, hh_174, hh_176, hh_183, hh_185, hh_316, hh_321, \
                         hh_323, hh_330, hh_332, hh_358, hh_363, hh_365, hh_372, hh_374, \
                         hh_400, hh_405, hh_407, hh_414, hh_416 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_76 * hh_22[k]
                  - f_80 * hh_27[k]
                  + f_79 * hh_29[k]
                  + f_86 * hh_36[k]
                  - f_88 * hh_38[k]
                  - f_77 * hh_127[k]
                  - f_81 * hh_132[k]
                  + f_83 * hh_134[k]
                  + f_80 * hh_141[k]
                  - f_82 * hh_143[k]
                  + f_78 * hh_169[k]
                  + f_79 * hh_174[k]
                  - f_84 * hh_176[k]
                  - f_87 * hh_183[k]
                  + f_89 * hh_185[k]
                  - f_76 * hh_316[k]
                  - f_80 * hh_321[k]
                  + f_79 * hh_323[k]
                  + f_86 * hh_330[k]
                  - f_88 * hh_332[k]
                  + f_78 * hh_358[k]
                  + f_79 * hh_363[k]
                  - f_84 * hh_365[k]
                  - f_87 * hh_372[k]
                  + f_89 * hh_374[k]
                  - f_79 * hh_400[k]
                  - f_82 * hh_405[k]
                  + f_85 * hh_407[k]
                  + f_88 * hh_414[k]
                  - f_90 * hh_416[k];
    }

#pragma omp simd aligned(hh_25, hh_32, hh_34, hh_130, hh_137, hh_139, hh_172, hh_179, hh_181, \
                         hh_319, hh_326, hh_328, hh_361, hh_368, hh_370, hh_403, hh_410, \
                         hh_412 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_120 * hh_25[k]
                  - f_120 * hh_32[k]
                  + f_121 * hh_34[k]
                  - f_121 * hh_130[k]
                  - f_121 * hh_137[k]
                  + f_124 * hh_139[k]
                  + f_122 * hh_172[k]
                  + f_122 * hh_179[k]
                  - f_125 * hh_181[k]
                  - f_120 * hh_319[k]
                  - f_120 * hh_326[k]
                  + f_121 * hh_328[k]
                  + f_122 * hh_361[k]
                  + f_122 * hh_368[k]
                  - f_125 * hh_370[k]
                  - f_123 * hh_403[k]
                  - f_123 * hh_410[k]
                  + f_126 * hh_412[k];
    }

#pragma omp simd aligned(hh_22, hh_27, hh_29, hh_36, hh_38, hh_40, hh_127, hh_132, hh_134, \
                         hh_141, hh_143, hh_145, hh_169, hh_174, hh_176, hh_183, hh_185, \
                         hh_187, hh_316, hh_321, hh_323, hh_330, hh_332, hh_334, hh_358, \
                         hh_363, hh_365, hh_372, hh_374, hh_376, hh_400, hh_405, hh_407, \
                         hh_414, hh_416, hh_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = 0.234375 * hh_22[k]
                  + 0.46875 * hh_27[k]
                  - 2.8125 * hh_29[k]
                  + 0.234375 * hh_36[k]
                  - 2.8125 * hh_38[k]
                  + 1.875 * hh_40[k]
                  + 0.46875 * hh_127[k]
                  + 0.9375 * hh_132[k]
                  - 5.625 * hh_134[k]
                  + 0.46875 * hh_141[k]
                  - 5.625 * hh_143[k]
                  + 3.75 * hh_145[k]
                  - 2.8125 * hh_169[k]
                  - 5.625 * hh_174[k]
                  + 33.75 * hh_176[k]
                  - 2.8125 * hh_183[k]
                  + 33.75 * hh_185[k]
                  - 22.5 * hh_187[k]
                  + 0.234375 * hh_316[k]
                  + 0.46875 * hh_321[k]
                  - 2.8125 * hh_323[k]
                  + 0.234375 * hh_330[k]
                  - 2.8125 * hh_332[k]
                  + 1.875 * hh_334[k]
                  - 2.8125 * hh_358[k]
                  - 5.625 * hh_363[k]
                  + 33.75 * hh_365[k]
                  - 2.8125 * hh_372[k]
                  + 33.75 * hh_374[k]
                  - 22.5 * hh_376[k]
                  + 1.875 * hh_400[k]
                  + 3.75 * hh_405[k]
                  - 22.5 * hh_407[k]
                  + 1.875 * hh_414[k]
                  - 22.5 * hh_416[k]
                  + 15.0 * hh_418[k];
    }

#pragma omp simd aligned(hh_23, hh_28, hh_30, hh_37, hh_39, hh_41, hh_128, hh_133, hh_135, \
                         hh_142, hh_144, hh_146, hh_170, hh_175, hh_177, hh_184, hh_186, \
                         hh_188, hh_317, hh_322, hh_324, hh_331, hh_333, hh_335, hh_359, \
                         hh_364, hh_366, hh_373, hh_375, hh_377, hh_401, hh_406, hh_408, \
                         hh_415, hh_417, hh_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_137 * hh_23[k]
                  + f_138 * hh_28[k]
                  - f_139 * hh_30[k]
                  + f_137 * hh_37[k]
                  - f_139 * hh_39[k]
                  + f_140 * hh_41[k]
                  + f_138 * hh_128[k]
                  + f_141 * hh_133[k]
                  - f_142 * hh_135[k]
                  + f_138 * hh_142[k]
                  - f_142 * hh_144[k]
                  + f_143 * hh_146[k]
                  - f_144 * hh_170[k]
                  - f_145 * hh_175[k]
                  + f_146 * hh_177[k]
                  - f_144 * hh_184[k]
                  + f_146 * hh_186[k]
                  - f_147 * hh_188[k]
                  + f_137 * hh_317[k]
                  + f_138 * hh_322[k]
                  - f_139 * hh_324[k]
                  + f_137 * hh_331[k]
                  - f_139 * hh_333[k]
                  + f_140 * hh_335[k]
                  - f_144 * hh_359[k]
                  - f_145 * hh_364[k]
                  + f_146 * hh_366[k]
                  - f_144 * hh_373[k]
                  + f_146 * hh_375[k]
                  - f_147 * hh_377[k]
                  + f_148 * hh_401[k]
                  + f_149 * hh_406[k]
                  - f_150 * hh_408[k]
                  + f_148 * hh_415[k]
                  - f_150 * hh_417[k]
                  + f_151 * hh_419[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_26, hh_31, hh_33, hh_35, hh_126, hh_129, hh_131, \
                         hh_136, hh_138, hh_140, hh_168, hh_171, hh_173, hh_178, hh_180, \
                         hh_182, hh_315, hh_318, hh_320, hh_325, hh_327, hh_329, hh_357, \
                         hh_360, hh_362, hh_367, hh_369, hh_371, hh_399, hh_402, hh_404, \
                         hh_409, hh_411, hh_413 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = 0.234375 * hh_21[k]
                  + 0.46875 * hh_24[k]
                  - 2.8125 * hh_26[k]
                  + 0.234375 * hh_31[k]
                  - 2.8125 * hh_33[k]
                  + 1.875 * hh_35[k]
                  + 0.46875 * hh_126[k]
                  + 0.9375 * hh_129[k]
                  - 5.625 * hh_131[k]
                  + 0.46875 * hh_136[k]
                  - 5.625 * hh_138[k]
                  + 3.75 * hh_140[k]
                  - 2.8125 * hh_168[k]
                  - 5.625 * hh_171[k]
                  + 33.75 * hh_173[k]
                  - 2.8125 * hh_178[k]
                  + 33.75 * hh_180[k]
                  - 22.5 * hh_182[k]
                  + 0.234375 * hh_315[k]
                  + 0.46875 * hh_318[k]
                  - 2.8125 * hh_320[k]
                  + 0.234375 * hh_325[k]
                  - 2.8125 * hh_327[k]
                  + 1.875 * hh_329[k]
                  - 2.8125 * hh_357[k]
                  - 5.625 * hh_360[k]
                  + 33.75 * hh_362[k]
                  - 2.8125 * hh_367[k]
                  + 33.75 * hh_369[k]
                  - 22.5 * hh_371[k]
                  + 1.875 * hh_399[k]
                  + 3.75 * hh_402[k]
                  - 22.5 * hh_404[k]
                  + 1.875 * hh_409[k]
                  - 22.5 * hh_411[k]
                  + 15.0 * hh_413[k];
    }

#pragma omp simd aligned(hh_23, hh_30, hh_37, hh_39, hh_128, hh_135, hh_142, hh_144, hh_170, \
                         hh_177, hh_184, hh_186, hh_317, hh_324, hh_331, hh_333, hh_359, \
                         hh_366, hh_373, hh_375, hh_401, hh_408, hh_415, \
                         hh_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_152 * hh_23[k]
                  + f_120 * hh_30[k]
                  + f_152 * hh_37[k]
                  - f_120 * hh_39[k]
                  - f_120 * hh_128[k]
                  + f_121 * hh_135[k]
                  + f_120 * hh_142[k]
                  - f_121 * hh_144[k]
                  + f_153 * hh_170[k]
                  - f_122 * hh_177[k]
                  - f_153 * hh_184[k]
                  + f_122 * hh_186[k]
                  - f_152 * hh_317[k]
                  + f_120 * hh_324[k]
                  + f_152 * hh_331[k]
                  - f_120 * hh_333[k]
                  + f_153 * hh_359[k]
                  - f_122 * hh_366[k]
                  - f_153 * hh_373[k]
                  + f_122 * hh_375[k]
                  - f_124 * hh_401[k]
                  + f_123 * hh_408[k]
                  + f_124 * hh_415[k]
                  - f_123 * hh_417[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_26, hh_31, hh_33, hh_126, hh_129, hh_131, hh_136, \
                         hh_138, hh_168, hh_171, hh_173, hh_178, hh_180, hh_315, hh_318, \
                         hh_320, hh_325, hh_327, hh_357, hh_360, hh_362, hh_367, hh_369, \
                         hh_399, hh_402, hh_404, hh_409, hh_411 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_86 * hh_21[k]
                  + f_80 * hh_24[k]
                  + f_88 * hh_26[k]
                  + f_76 * hh_31[k]
                  - f_79 * hh_33[k]
                  - f_80 * hh_126[k]
                  + f_81 * hh_129[k]
                  + f_82 * hh_131[k]
                  + f_77 * hh_136[k]
                  - f_83 * hh_138[k]
                  + f_87 * hh_168[k]
                  - f_79 * hh_171[k]
                  - f_89 * hh_173[k]
                  - f_78 * hh_178[k]
                  + f_84 * hh_180[k]
                  - f_86 * hh_315[k]
                  + f_80 * hh_318[k]
                  + f_88 * hh_320[k]
                  + f_76 * hh_325[k]
                  - f_79 * hh_327[k]
                  + f_87 * hh_357[k]
                  - f_79 * hh_360[k]
                  - f_89 * hh_362[k]
                  - f_78 * hh_367[k]
                  + f_84 * hh_369[k]
                  - f_88 * hh_399[k]
                  + f_82 * hh_402[k]
                  + f_90 * hh_404[k]
                  + f_79 * hh_409[k]
                  - f_85 * hh_411[k];
    }

#pragma omp simd aligned(hh_23, hh_28, hh_37, hh_128, hh_133, hh_142, hh_170, hh_175, hh_184, \
                         hh_317, hh_322, hh_331, hh_359, hh_364, hh_373, hh_401, hh_406, \
                         hh_415 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_154 * hh_23[k]
                  - f_155 * hh_28[k]
                  + f_154 * hh_37[k]
                  + f_156 * hh_128[k]
                  - f_157 * hh_133[k]
                  + f_156 * hh_142[k]
                  - f_157 * hh_170[k]
                  + f_158 * hh_175[k]
                  - f_157 * hh_184[k]
                  + f_154 * hh_317[k]
                  - f_155 * hh_322[k]
                  + f_154 * hh_331[k]
                  - f_157 * hh_359[k]
                  + f_158 * hh_364[k]
                  - f_157 * hh_373[k]
                  + f_59 * hh_401[k]
                  - f_60 * hh_406[k]
                  + f_59 * hh_415[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_31, hh_126, hh_129, hh_136, hh_168, hh_171, hh_178, \
                         hh_315, hh_318, hh_325, hh_357, hh_360, hh_367, hh_399, hh_402, \
                         hh_409 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_29 * hh_21[k]
                  - f_23 * hh_24[k]
                  + f_22 * hh_31[k]
                  + f_30 * hh_126[k]
                  - f_26 * hh_129[k]
                  + f_23 * hh_136[k]
                  - f_31 * hh_168[k]
                  + f_27 * hh_171[k]
                  - f_24 * hh_178[k]
                  + f_29 * hh_315[k]
                  - f_23 * hh_318[k]
                  + f_22 * hh_325[k]
                  - f_31 * hh_357[k]
                  + f_27 * hh_360[k]
                  - f_24 * hh_367[k]
                  + f_32 * hh_399[k]
                  - f_28 * hh_402[k]
                  + f_25 * hh_409[k];
    }

#pragma omp simd aligned(hh_43, hh_48, hh_57, hh_148, hh_153, hh_162, hh_190, hh_195, hh_204, \
                         hh_337, hh_342, hh_351, hh_379, hh_384, hh_393, hh_421, hh_426, \
                         hh_435 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_33 * hh_43[k]
                  - f_34 * hh_48[k]
                  + f_40 * hh_57[k]
                  + f_34 * hh_148[k]
                  - f_37 * hh_153[k]
                  + f_41 * hh_162[k]
                  - f_35 * hh_190[k]
                  + f_38 * hh_195[k]
                  - f_36 * hh_204[k]
                  + f_33 * hh_337[k]
                  - f_34 * hh_342[k]
                  + f_40 * hh_351[k]
                  - f_35 * hh_379[k]
                  + f_38 * hh_384[k]
                  - f_36 * hh_393[k]
                  + f_36 * hh_421[k]
                  - f_39 * hh_426[k]
                  + f_42 * hh_435[k];
    }

#pragma omp simd aligned(hh_46, hh_53, hh_151, hh_158, hh_193, hh_200, hh_340, hh_347, hh_382, \
                         hh_389, hh_424, hh_431 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_62 * hh_46[k]
                  - f_62 * hh_53[k]
                  + f_63 * hh_151[k]
                  - f_63 * hh_158[k]
                  - f_64 * hh_193[k]
                  + f_64 * hh_200[k]
                  + f_62 * hh_340[k]
                  - f_62 * hh_347[k]
                  - f_64 * hh_382[k]
                  + f_64 * hh_389[k]
                  + f_65 * hh_424[k]
                  - f_65 * hh_431[k];
    }

#pragma omp simd aligned(hh_43, hh_48, hh_50, hh_57, hh_59, hh_148, hh_153, hh_155, hh_162, \
                         hh_164, hh_190, hh_195, hh_197, hh_204, hh_206, hh_337, hh_342, \
                         hh_344, hh_351, hh_353, hh_379, hh_384, hh_386, hh_393, hh_395, \
                         hh_421, hh_426, hh_428, hh_435, hh_437 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = -f_91 * hh_43[k]
                  - f_95 * hh_48[k]
                  + f_99 * hh_50[k]
                  + f_103 * hh_57[k]
                  - f_93 * hh_59[k]
                  - f_92 * hh_148[k]
                  - f_96 * hh_153[k]
                  + f_100 * hh_155[k]
                  + f_95 * hh_162[k]
                  - f_106 * hh_164[k]
                  + f_93 * hh_190[k]
                  + f_97 * hh_195[k]
                  - f_101 * hh_197[k]
                  - f_104 * hh_204[k]
                  + f_107 * hh_206[k]
                  - f_91 * hh_337[k]
                  - f_95 * hh_342[k]
                  + f_99 * hh_344[k]
                  + f_103 * hh_351[k]
                  - f_93 * hh_353[k]
                  + f_93 * hh_379[k]
                  + f_97 * hh_384[k]
                  - f_101 * hh_386[k]
                  - f_104 * hh_393[k]
                  + f_107 * hh_395[k]
                  - f_94 * hh_421[k]
                  - f_98 * hh_426[k]
                  + f_102 * hh_428[k]
                  + f_105 * hh_435[k]
                  - f_108 * hh_437[k];
    }

#pragma omp simd aligned(hh_46, hh_53, hh_55, hh_151, hh_158, hh_160, hh_193, hh_200, hh_202, \
                         hh_340, hh_347, hh_349, hh_382, hh_389, hh_391, hh_424, hh_431, \
                         hh_433 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_127 * hh_46[k]
                  - f_127 * hh_53[k]
                  + f_128 * hh_55[k]
                  - f_128 * hh_151[k]
                  - f_128 * hh_158[k]
                  + f_131 * hh_160[k]
                  + f_129 * hh_193[k]
                  + f_129 * hh_200[k]
                  - f_132 * hh_202[k]
                  - f_127 * hh_340[k]
                  - f_127 * hh_347[k]
                  + f_128 * hh_349[k]
                  + f_129 * hh_382[k]
                  + f_129 * hh_389[k]
                  - f_132 * hh_391[k]
                  - f_130 * hh_424[k]
                  - f_130 * hh_431[k]
                  + f_133 * hh_433[k];
    }

#pragma omp simd aligned(hh_43, hh_48, hh_50, hh_57, hh_59, hh_61, hh_148, hh_153, hh_155, \
                         hh_162, hh_164, hh_166, hh_190, hh_195, hh_197, hh_204, hh_206, \
                         hh_208, hh_337, hh_342, hh_344, hh_351, hh_353, hh_355, hh_379, \
                         hh_384, hh_386, hh_393, hh_395, hh_397, hh_421, hh_426, hh_428, \
                         hh_435, hh_437, hh_439 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_137 * hh_43[k]
                  + f_138 * hh_48[k]
                  - f_144 * hh_50[k]
                  + f_137 * hh_57[k]
                  - f_144 * hh_59[k]
                  + f_148 * hh_61[k]
                  + f_138 * hh_148[k]
                  + f_141 * hh_153[k]
                  - f_145 * hh_155[k]
                  + f_138 * hh_162[k]
                  - f_145 * hh_164[k]
                  + f_149 * hh_166[k]
                  - f_139 * hh_190[k]
                  - f_142 * hh_195[k]
                  + f_146 * hh_197[k]
                  - f_139 * hh_204[k]
                  + f_146 * hh_206[k]
                  - f_150 * hh_208[k]
                  + f_137 * hh_337[k]
                  + f_138 * hh_342[k]
                  - f_144 * hh_344[k]
                  + f_137 * hh_351[k]
                  - f_144 * hh_353[k]
                  + f_148 * hh_355[k]
                  - f_139 * hh_379[k]
                  - f_142 * hh_384[k]
                  + f_146 * hh_386[k]
                  - f_139 * hh_393[k]
                  + f_146 * hh_395[k]
                  - f_150 * hh_397[k]
                  + f_140 * hh_421[k]
                  + f_143 * hh_426[k]
                  - f_147 * hh_428[k]
                  + f_140 * hh_435[k]
                  - f_147 * hh_437[k]
                  + f_151 * hh_439[k];
    }

#pragma omp simd aligned(hh_44, hh_49, hh_51, hh_58, hh_60, hh_62, hh_149, hh_154, hh_156, \
                         hh_163, hh_165, hh_167, hh_191, hh_196, hh_198, hh_205, hh_207, \
                         hh_209, hh_338, hh_343, hh_345, hh_352, hh_354, hh_356, hh_380, \
                         hh_385, hh_387, hh_394, hh_396, hh_398, hh_422, hh_427, hh_429, \
                         hh_436, hh_438, hh_440 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = 3.515625 * hh_44[k]
                  + 7.03125 * hh_49[k]
                  - 9.375 * hh_51[k]
                  + 3.515625 * hh_58[k]
                  - 9.375 * hh_60[k]
                  + 1.875 * hh_62[k]
                  + 7.03125 * hh_149[k]
                  + 14.0625 * hh_154[k]
                  - 18.75 * hh_156[k]
                  + 7.03125 * hh_163[k]
                  - 18.75 * hh_165[k]
                  + 3.75 * hh_167[k]
                  - 9.375 * hh_191[k]
                  - 18.75 * hh_196[k]
                  + 25.0 * hh_198[k]
                  - 9.375 * hh_205[k]
                  + 25.0 * hh_207[k]
                  - 5.0 * hh_209[k]
                  + 3.515625 * hh_338[k]
                  + 7.03125 * hh_343[k]
                  - 9.375 * hh_345[k]
                  + 3.515625 * hh_352[k]
                  - 9.375 * hh_354[k]
                  + 1.875 * hh_356[k]
                  - 9.375 * hh_380[k]
                  - 18.75 * hh_385[k]
                  + 25.0 * hh_387[k]
                  - 9.375 * hh_394[k]
                  + 25.0 * hh_396[k]
                  - 5.0 * hh_398[k]
                  + 1.875 * hh_422[k]
                  + 3.75 * hh_427[k]
                  - 5.0 * hh_429[k]
                  + 1.875 * hh_436[k]
                  - 5.0 * hh_438[k]
                  + hh_440[k];
    }

#pragma omp simd aligned(hh_42, hh_45, hh_47, hh_52, hh_54, hh_56, hh_147, hh_150, hh_152, \
                         hh_157, hh_159, hh_161, hh_189, hh_192, hh_194, hh_199, hh_201, \
                         hh_203, hh_336, hh_339, hh_341, hh_346, hh_348, hh_350, hh_378, \
                         hh_381, hh_383, hh_388, hh_390, hh_392, hh_420, hh_423, hh_425, \
                         hh_430, hh_432, hh_434 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_137 * hh_42[k]
                  + f_138 * hh_45[k]
                  - f_144 * hh_47[k]
                  + f_137 * hh_52[k]
                  - f_144 * hh_54[k]
                  + f_148 * hh_56[k]
                  + f_138 * hh_147[k]
                  + f_141 * hh_150[k]
                  - f_145 * hh_152[k]
                  + f_138 * hh_157[k]
                  - f_145 * hh_159[k]
                  + f_149 * hh_161[k]
                  - f_139 * hh_189[k]
                  - f_142 * hh_192[k]
                  + f_146 * hh_194[k]
                  - f_139 * hh_199[k]
                  + f_146 * hh_201[k]
                  - f_150 * hh_203[k]
                  + f_137 * hh_336[k]
                  + f_138 * hh_339[k]
                  - f_144 * hh_341[k]
                  + f_137 * hh_346[k]
                  - f_144 * hh_348[k]
                  + f_148 * hh_350[k]
                  - f_139 * hh_378[k]
                  - f_142 * hh_381[k]
                  + f_146 * hh_383[k]
                  - f_139 * hh_388[k]
                  + f_146 * hh_390[k]
                  - f_150 * hh_392[k]
                  + f_140 * hh_420[k]
                  + f_143 * hh_423[k]
                  - f_147 * hh_425[k]
                  + f_140 * hh_430[k]
                  - f_147 * hh_432[k]
                  + f_151 * hh_434[k];
    }

#pragma omp simd aligned(hh_44, hh_51, hh_58, hh_60, hh_149, hh_156, hh_163, hh_165, hh_191, \
                         hh_198, hh_205, hh_207, hh_338, hh_345, hh_352, hh_354, hh_380, \
                         hh_387, hh_394, hh_396, hh_422, hh_429, hh_436, \
                         hh_438 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_159 * hh_44[k]
                  + f_127 * hh_51[k]
                  + f_159 * hh_58[k]
                  - f_127 * hh_60[k]
                  - f_127 * hh_149[k]
                  + f_128 * hh_156[k]
                  + f_127 * hh_163[k]
                  - f_128 * hh_165[k]
                  + f_160 * hh_191[k]
                  - f_129 * hh_198[k]
                  - f_160 * hh_205[k]
                  + f_129 * hh_207[k]
                  - f_159 * hh_338[k]
                  + f_127 * hh_345[k]
                  + f_159 * hh_352[k]
                  - f_127 * hh_354[k]
                  + f_160 * hh_380[k]
                  - f_129 * hh_387[k]
                  - f_160 * hh_394[k]
                  + f_129 * hh_396[k]
                  - f_161 * hh_422[k]
                  + f_130 * hh_429[k]
                  + f_161 * hh_436[k]
                  - f_130 * hh_438[k];
    }

#pragma omp simd aligned(hh_42, hh_45, hh_47, hh_52, hh_54, hh_147, hh_150, hh_152, hh_157, \
                         hh_159, hh_189, hh_192, hh_194, hh_199, hh_201, hh_336, hh_339, \
                         hh_341, hh_346, hh_348, hh_378, hh_381, hh_383, hh_388, hh_390, \
                         hh_420, hh_423, hh_425, hh_430, hh_432 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_103 * hh_42[k]
                  + f_95 * hh_45[k]
                  + f_93 * hh_47[k]
                  + f_91 * hh_52[k]
                  - f_99 * hh_54[k]
                  - f_95 * hh_147[k]
                  + f_96 * hh_150[k]
                  + f_106 * hh_152[k]
                  + f_92 * hh_157[k]
                  - f_100 * hh_159[k]
                  + f_104 * hh_189[k]
                  - f_97 * hh_192[k]
                  - f_107 * hh_194[k]
                  - f_93 * hh_199[k]
                  + f_101 * hh_201[k]
                  - f_103 * hh_336[k]
                  + f_95 * hh_339[k]
                  + f_93 * hh_341[k]
                  + f_91 * hh_346[k]
                  - f_99 * hh_348[k]
                  + f_104 * hh_378[k]
                  - f_97 * hh_381[k]
                  - f_107 * hh_383[k]
                  - f_93 * hh_388[k]
                  + f_101 * hh_390[k]
                  - f_105 * hh_420[k]
                  + f_98 * hh_423[k]
                  + f_108 * hh_425[k]
                  + f_94 * hh_430[k]
                  - f_102 * hh_432[k];
    }

#pragma omp simd aligned(hh_44, hh_49, hh_58, hh_149, hh_154, hh_163, hh_191, hh_196, hh_205, \
                         hh_338, hh_343, hh_352, hh_380, hh_385, hh_394, hh_422, hh_427, \
                         hh_436 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_162 * hh_44[k]
                  - f_163 * hh_49[k]
                  + f_162 * hh_58[k]
                  + f_164 * hh_149[k]
                  - f_165 * hh_154[k]
                  + f_164 * hh_163[k]
                  - f_166 * hh_191[k]
                  + f_167 * hh_196[k]
                  - f_166 * hh_205[k]
                  + f_162 * hh_338[k]
                  - f_163 * hh_343[k]
                  + f_162 * hh_352[k]
                  - f_166 * hh_380[k]
                  + f_167 * hh_385[k]
                  - f_166 * hh_394[k]
                  + f_168 * hh_422[k]
                  - f_169 * hh_427[k]
                  + f_168 * hh_436[k];
    }

#pragma omp simd aligned(hh_42, hh_45, hh_52, hh_147, hh_150, hh_157, hh_189, hh_192, hh_199, \
                         hh_336, hh_339, hh_346, hh_378, hh_381, hh_388, hh_420, hh_423, \
                         hh_430 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_40 * hh_42[k]
                  - f_34 * hh_45[k]
                  + f_33 * hh_52[k]
                  + f_41 * hh_147[k]
                  - f_37 * hh_150[k]
                  + f_34 * hh_157[k]
                  - f_36 * hh_189[k]
                  + f_38 * hh_192[k]
                  - f_35 * hh_199[k]
                  + f_40 * hh_336[k]
                  - f_34 * hh_339[k]
                  + f_33 * hh_346[k]
                  - f_36 * hh_378[k]
                  + f_38 * hh_381[k]
                  - f_35 * hh_388[k]
                  + f_42 * hh_420[k]
                  - f_39 * hh_423[k]
                  + f_36 * hh_430[k];
    }

#pragma omp simd aligned(hh_1, hh_6, hh_15, hh_64, hh_69, hh_78, hh_106, hh_111, hh_120, \
                         hh_211, hh_216, hh_225, hh_253, hh_258, hh_267, hh_295, hh_300, \
                         hh_309 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = f_22 * hh_1[k]
                  - f_23 * hh_6[k]
                  + f_29 * hh_15[k]
                  + f_23 * hh_64[k]
                  - f_26 * hh_69[k]
                  + f_30 * hh_78[k]
                  - f_24 * hh_106[k]
                  + f_27 * hh_111[k]
                  - f_31 * hh_120[k]
                  + f_22 * hh_211[k]
                  - f_23 * hh_216[k]
                  + f_29 * hh_225[k]
                  - f_24 * hh_253[k]
                  + f_27 * hh_258[k]
                  - f_31 * hh_267[k]
                  + f_25 * hh_295[k]
                  - f_28 * hh_300[k]
                  + f_32 * hh_309[k];
    }

#pragma omp simd aligned(hh_4, hh_11, hh_67, hh_74, hh_109, hh_116, hh_214, hh_221, hh_256, \
                         hh_263, hh_298, hh_305 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = f_58 * hh_4[k]
                  - f_58 * hh_11[k]
                  + f_59 * hh_67[k]
                  - f_59 * hh_74[k]
                  - f_60 * hh_109[k]
                  + f_60 * hh_116[k]
                  + f_58 * hh_214[k]
                  - f_58 * hh_221[k]
                  - f_60 * hh_256[k]
                  + f_60 * hh_263[k]
                  + f_61 * hh_298[k]
                  - f_61 * hh_305[k];
    }

#pragma omp simd aligned(hh_1, hh_6, hh_8, hh_15, hh_17, hh_64, hh_69, hh_71, hh_78, hh_80, \
                         hh_106, hh_111, hh_113, hh_120, hh_122, hh_211, hh_216, hh_218, \
                         hh_225, hh_227, hh_253, hh_258, hh_260, hh_267, hh_269, hh_295, \
                         hh_300, hh_302, hh_309, hh_311 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = -f_76 * hh_1[k]
                  - f_80 * hh_6[k]
                  + f_79 * hh_8[k]
                  + f_86 * hh_15[k]
                  - f_88 * hh_17[k]
                  - f_77 * hh_64[k]
                  - f_81 * hh_69[k]
                  + f_83 * hh_71[k]
                  + f_80 * hh_78[k]
                  - f_82 * hh_80[k]
                  + f_78 * hh_106[k]
                  + f_79 * hh_111[k]
                  - f_84 * hh_113[k]
                  - f_87 * hh_120[k]
                  + f_89 * hh_122[k]
                  - f_76 * hh_211[k]
                  - f_80 * hh_216[k]
                  + f_79 * hh_218[k]
                  + f_86 * hh_225[k]
                  - f_88 * hh_227[k]
                  + f_78 * hh_253[k]
                  + f_79 * hh_258[k]
                  - f_84 * hh_260[k]
                  - f_87 * hh_267[k]
                  + f_89 * hh_269[k]
                  - f_79 * hh_295[k]
                  - f_82 * hh_300[k]
                  + f_85 * hh_302[k]
                  + f_88 * hh_309[k]
                  - f_90 * hh_311[k];
    }

#pragma omp simd aligned(hh_4, hh_11, hh_13, hh_67, hh_74, hh_76, hh_109, hh_116, hh_118, \
                         hh_214, hh_221, hh_223, hh_256, hh_263, hh_265, hh_298, hh_305, \
                         hh_307 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = -f_120 * hh_4[k]
                  - f_120 * hh_11[k]
                  + f_121 * hh_13[k]
                  - f_121 * hh_67[k]
                  - f_121 * hh_74[k]
                  + f_124 * hh_76[k]
                  + f_122 * hh_109[k]
                  + f_122 * hh_116[k]
                  - f_125 * hh_118[k]
                  - f_120 * hh_214[k]
                  - f_120 * hh_221[k]
                  + f_121 * hh_223[k]
                  + f_122 * hh_256[k]
                  + f_122 * hh_263[k]
                  - f_125 * hh_265[k]
                  - f_123 * hh_298[k]
                  - f_123 * hh_305[k]
                  + f_126 * hh_307[k];
    }

#pragma omp simd aligned(hh_1, hh_6, hh_8, hh_15, hh_17, hh_19, hh_64, hh_69, hh_71, hh_78, \
                         hh_80, hh_82, hh_106, hh_111, hh_113, hh_120, hh_122, hh_124, hh_211, \
                         hh_216, hh_218, hh_225, hh_227, hh_229, hh_253, hh_258, hh_260, \
                         hh_267, hh_269, hh_271, hh_295, hh_300, hh_302, hh_309, hh_311, \
                         hh_313 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = 0.234375 * hh_1[k]
                  + 0.46875 * hh_6[k]
                  - 2.8125 * hh_8[k]
                  + 0.234375 * hh_15[k]
                  - 2.8125 * hh_17[k]
                  + 1.875 * hh_19[k]
                  + 0.46875 * hh_64[k]
                  + 0.9375 * hh_69[k]
                  - 5.625 * hh_71[k]
                  + 0.46875 * hh_78[k]
                  - 5.625 * hh_80[k]
                  + 3.75 * hh_82[k]
                  - 2.8125 * hh_106[k]
                  - 5.625 * hh_111[k]
                  + 33.75 * hh_113[k]
                  - 2.8125 * hh_120[k]
                  + 33.75 * hh_122[k]
                  - 22.5 * hh_124[k]
                  + 0.234375 * hh_211[k]
                  + 0.46875 * hh_216[k]
                  - 2.8125 * hh_218[k]
                  + 0.234375 * hh_225[k]
                  - 2.8125 * hh_227[k]
                  + 1.875 * hh_229[k]
                  - 2.8125 * hh_253[k]
                  - 5.625 * hh_258[k]
                  + 33.75 * hh_260[k]
                  - 2.8125 * hh_267[k]
                  + 33.75 * hh_269[k]
                  - 22.5 * hh_271[k]
                  + 1.875 * hh_295[k]
                  + 3.75 * hh_300[k]
                  - 22.5 * hh_302[k]
                  + 1.875 * hh_309[k]
                  - 22.5 * hh_311[k]
                  + 15.0 * hh_313[k];
    }

#pragma omp simd aligned(hh_2, hh_7, hh_9, hh_16, hh_18, hh_20, hh_65, hh_70, hh_72, hh_79, \
                         hh_81, hh_83, hh_107, hh_112, hh_114, hh_121, hh_123, hh_125, hh_212, \
                         hh_217, hh_219, hh_226, hh_228, hh_230, hh_254, hh_259, hh_261, \
                         hh_268, hh_270, hh_272, hh_296, hh_301, hh_303, hh_310, hh_312, \
                         hh_314 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = f_137 * hh_2[k]
                  + f_138 * hh_7[k]
                  - f_139 * hh_9[k]
                  + f_137 * hh_16[k]
                  - f_139 * hh_18[k]
                  + f_140 * hh_20[k]
                  + f_138 * hh_65[k]
                  + f_141 * hh_70[k]
                  - f_142 * hh_72[k]
                  + f_138 * hh_79[k]
                  - f_142 * hh_81[k]
                  + f_143 * hh_83[k]
                  - f_144 * hh_107[k]
                  - f_145 * hh_112[k]
                  + f_146 * hh_114[k]
                  - f_144 * hh_121[k]
                  + f_146 * hh_123[k]
                  - f_147 * hh_125[k]
                  + f_137 * hh_212[k]
                  + f_138 * hh_217[k]
                  - f_139 * hh_219[k]
                  + f_137 * hh_226[k]
                  - f_139 * hh_228[k]
                  + f_140 * hh_230[k]
                  - f_144 * hh_254[k]
                  - f_145 * hh_259[k]
                  + f_146 * hh_261[k]
                  - f_144 * hh_268[k]
                  + f_146 * hh_270[k]
                  - f_147 * hh_272[k]
                  + f_148 * hh_296[k]
                  + f_149 * hh_301[k]
                  - f_150 * hh_303[k]
                  + f_148 * hh_310[k]
                  - f_150 * hh_312[k]
                  + f_151 * hh_314[k];
    }

#pragma omp simd aligned(hh_0, hh_3, hh_5, hh_10, hh_12, hh_14, hh_63, hh_66, hh_68, hh_73, \
                         hh_75, hh_77, hh_105, hh_108, hh_110, hh_115, hh_117, hh_119, hh_210, \
                         hh_213, hh_215, hh_220, hh_222, hh_224, hh_252, hh_255, hh_257, \
                         hh_262, hh_264, hh_266, hh_294, hh_297, hh_299, hh_304, hh_306, \
                         hh_308 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = 0.234375 * hh_0[k]
                  + 0.46875 * hh_3[k]
                  - 2.8125 * hh_5[k]
                  + 0.234375 * hh_10[k]
                  - 2.8125 * hh_12[k]
                  + 1.875 * hh_14[k]
                  + 0.46875 * hh_63[k]
                  + 0.9375 * hh_66[k]
                  - 5.625 * hh_68[k]
                  + 0.46875 * hh_73[k]
                  - 5.625 * hh_75[k]
                  + 3.75 * hh_77[k]
                  - 2.8125 * hh_105[k]
                  - 5.625 * hh_108[k]
                  + 33.75 * hh_110[k]
                  - 2.8125 * hh_115[k]
                  + 33.75 * hh_117[k]
                  - 22.5 * hh_119[k]
                  + 0.234375 * hh_210[k]
                  + 0.46875 * hh_213[k]
                  - 2.8125 * hh_215[k]
                  + 0.234375 * hh_220[k]
                  - 2.8125 * hh_222[k]
                  + 1.875 * hh_224[k]
                  - 2.8125 * hh_252[k]
                  - 5.625 * hh_255[k]
                  + 33.75 * hh_257[k]
                  - 2.8125 * hh_262[k]
                  + 33.75 * hh_264[k]
                  - 22.5 * hh_266[k]
                  + 1.875 * hh_294[k]
                  + 3.75 * hh_297[k]
                  - 22.5 * hh_299[k]
                  + 1.875 * hh_304[k]
                  - 22.5 * hh_306[k]
                  + 15.0 * hh_308[k];
    }

#pragma omp simd aligned(hh_2, hh_9, hh_16, hh_18, hh_65, hh_72, hh_79, hh_81, hh_107, hh_114, \
                         hh_121, hh_123, hh_212, hh_219, hh_226, hh_228, hh_254, hh_261, \
                         hh_268, hh_270, hh_296, hh_303, hh_310, \
                         hh_312 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = -f_152 * hh_2[k]
                  + f_120 * hh_9[k]
                  + f_152 * hh_16[k]
                  - f_120 * hh_18[k]
                  - f_120 * hh_65[k]
                  + f_121 * hh_72[k]
                  + f_120 * hh_79[k]
                  - f_121 * hh_81[k]
                  + f_153 * hh_107[k]
                  - f_122 * hh_114[k]
                  - f_153 * hh_121[k]
                  + f_122 * hh_123[k]
                  - f_152 * hh_212[k]
                  + f_120 * hh_219[k]
                  + f_152 * hh_226[k]
                  - f_120 * hh_228[k]
                  + f_153 * hh_254[k]
                  - f_122 * hh_261[k]
                  - f_153 * hh_268[k]
                  + f_122 * hh_270[k]
                  - f_124 * hh_296[k]
                  + f_123 * hh_303[k]
                  + f_124 * hh_310[k]
                  - f_123 * hh_312[k];
    }

#pragma omp simd aligned(hh_0, hh_3, hh_5, hh_10, hh_12, hh_63, hh_66, hh_68, hh_73, hh_75, \
                         hh_105, hh_108, hh_110, hh_115, hh_117, hh_210, hh_213, hh_215, \
                         hh_220, hh_222, hh_252, hh_255, hh_257, hh_262, hh_264, hh_294, \
                         hh_297, hh_299, hh_304, hh_306 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_86 * hh_0[k]
                  + f_80 * hh_3[k]
                  + f_88 * hh_5[k]
                  + f_76 * hh_10[k]
                  - f_79 * hh_12[k]
                  - f_80 * hh_63[k]
                  + f_81 * hh_66[k]
                  + f_82 * hh_68[k]
                  + f_77 * hh_73[k]
                  - f_83 * hh_75[k]
                  + f_87 * hh_105[k]
                  - f_79 * hh_108[k]
                  - f_89 * hh_110[k]
                  - f_78 * hh_115[k]
                  + f_84 * hh_117[k]
                  - f_86 * hh_210[k]
                  + f_80 * hh_213[k]
                  + f_88 * hh_215[k]
                  + f_76 * hh_220[k]
                  - f_79 * hh_222[k]
                  + f_87 * hh_252[k]
                  - f_79 * hh_255[k]
                  - f_89 * hh_257[k]
                  - f_78 * hh_262[k]
                  + f_84 * hh_264[k]
                  - f_88 * hh_294[k]
                  + f_82 * hh_297[k]
                  + f_90 * hh_299[k]
                  + f_79 * hh_304[k]
                  - f_85 * hh_306[k];
    }

#pragma omp simd aligned(hh_2, hh_7, hh_16, hh_65, hh_70, hh_79, hh_107, hh_112, hh_121, \
                         hh_212, hh_217, hh_226, hh_254, hh_259, hh_268, hh_296, hh_301, \
                         hh_310 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = f_154 * hh_2[k]
                  - f_155 * hh_7[k]
                  + f_154 * hh_16[k]
                  + f_156 * hh_65[k]
                  - f_157 * hh_70[k]
                  + f_156 * hh_79[k]
                  - f_157 * hh_107[k]
                  + f_158 * hh_112[k]
                  - f_157 * hh_121[k]
                  + f_154 * hh_212[k]
                  - f_155 * hh_217[k]
                  + f_154 * hh_226[k]
                  - f_157 * hh_254[k]
                  + f_158 * hh_259[k]
                  - f_157 * hh_268[k]
                  + f_59 * hh_296[k]
                  - f_60 * hh_301[k]
                  + f_59 * hh_310[k];
    }

#pragma omp simd aligned(hh_0, hh_3, hh_10, hh_63, hh_66, hh_73, hh_105, hh_108, hh_115, \
                         hh_210, hh_213, hh_220, hh_252, hh_255, hh_262, hh_294, hh_297, \
                         hh_304 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = f_29 * hh_0[k]
                  - f_23 * hh_3[k]
                  + f_22 * hh_10[k]
                  + f_30 * hh_63[k]
                  - f_26 * hh_66[k]
                  + f_23 * hh_73[k]
                  - f_31 * hh_105[k]
                  + f_27 * hh_108[k]
                  - f_24 * hh_115[k]
                  + f_29 * hh_210[k]
                  - f_23 * hh_213[k]
                  + f_22 * hh_220[k]
                  - f_31 * hh_252[k]
                  + f_27 * hh_255[k]
                  - f_24 * hh_262[k]
                  + f_32 * hh_294[k]
                  - f_28 * hh_297[k]
                  + f_25 * hh_304[k];
    }

#pragma omp simd aligned(hh_43, hh_48, hh_57, hh_190, hh_195, hh_204, hh_337, hh_342, hh_351, \
                         hh_379, hh_384, hh_393 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = -f_43 * hh_43[k]
                  + f_17 * hh_48[k]
                  - f_44 * hh_57[k]
                  + f_17 * hh_190[k]
                  - f_18 * hh_195[k]
                  + f_20 * hh_204[k]
                  + f_43 * hh_337[k]
                  - f_17 * hh_342[k]
                  + f_44 * hh_351[k]
                  - f_17 * hh_379[k]
                  + f_18 * hh_384[k]
                  - f_20 * hh_393[k];
    }

#pragma omp simd aligned(hh_46, hh_53, hh_193, hh_200, hh_340, hh_347, hh_382, \
                         hh_389 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = -f_66 * hh_46[k]
                  + f_66 * hh_53[k]
                  + f_56 * hh_193[k]
                  - f_56 * hh_200[k]
                  + f_66 * hh_340[k]
                  - f_66 * hh_347[k]
                  - f_56 * hh_382[k]
                  + f_56 * hh_389[k];
    }

#pragma omp simd aligned(hh_43, hh_48, hh_50, hh_57, hh_59, hh_190, hh_195, hh_197, hh_204, \
                         hh_206, hh_337, hh_342, hh_344, hh_351, hh_353, hh_379, hh_384, \
                         hh_386, hh_393, hh_395 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = f_109 * hh_43[k]
                  + f_73 * hh_48[k]
                  - f_110 * hh_50[k]
                  - f_111 * hh_57[k]
                  + f_70 * hh_59[k]
                  - f_67 * hh_190[k]
                  - f_69 * hh_195[k]
                  + f_71 * hh_197[k]
                  + f_73 * hh_204[k]
                  - f_74 * hh_206[k]
                  - f_109 * hh_337[k]
                  - f_73 * hh_342[k]
                  + f_110 * hh_344[k]
                  + f_111 * hh_351[k]
                  - f_70 * hh_353[k]
                  + f_67 * hh_379[k]
                  + f_69 * hh_384[k]
                  - f_71 * hh_386[k]
                  - f_73 * hh_393[k]
                  + f_74 * hh_395[k];
    }

#pragma omp simd aligned(hh_46, hh_53, hh_55, hh_193, hh_200, hh_202, hh_340, hh_347, hh_349, \
                         hh_382, hh_389, hh_391 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = 13.125 * hh_46[k]
                  + 13.125 * hh_53[k]
                  - 26.25 * hh_55[k]
                  - 26.25 * hh_193[k]
                  - 26.25 * hh_200[k]
                  + 52.5 * hh_202[k]
                  - 13.125 * hh_340[k]
                  - 13.125 * hh_347[k]
                  + 26.25 * hh_349[k]
                  + 26.25 * hh_382[k]
                  + 26.25 * hh_389[k]
                  - 52.5 * hh_391[k];
    }

#pragma omp simd aligned(hh_43, hh_48, hh_50, hh_57, hh_59, hh_61, hh_190, hh_195, hh_197, \
                         hh_204, hh_206, hh_208, hh_337, hh_342, hh_344, hh_351, hh_353, \
                         hh_355, hh_379, hh_384, hh_386, hh_393, hh_395, \
                         hh_397 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = -f_152 * hh_43[k]
                  - f_120 * hh_48[k]
                  + f_153 * hh_50[k]
                  - f_152 * hh_57[k]
                  + f_153 * hh_59[k]
                  - f_124 * hh_61[k]
                  + f_120 * hh_190[k]
                  + f_121 * hh_195[k]
                  - f_122 * hh_197[k]
                  + f_120 * hh_204[k]
                  - f_122 * hh_206[k]
                  + f_123 * hh_208[k]
                  + f_152 * hh_337[k]
                  + f_120 * hh_342[k]
                  - f_153 * hh_344[k]
                  + f_152 * hh_351[k]
                  - f_153 * hh_353[k]
                  + f_124 * hh_355[k]
                  - f_120 * hh_379[k]
                  - f_121 * hh_384[k]
                  + f_122 * hh_386[k]
                  - f_120 * hh_393[k]
                  + f_122 * hh_395[k]
                  - f_123 * hh_397[k];
    }

#pragma omp simd aligned(hh_44, hh_49, hh_51, hh_58, hh_60, hh_62, hh_191, hh_196, hh_198, \
                         hh_205, hh_207, hh_209, hh_338, hh_343, hh_345, hh_352, hh_354, \
                         hh_356, hh_380, hh_385, hh_387, hh_394, hh_396, \
                         hh_398 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = -f_159 * hh_44[k]
                  - f_127 * hh_49[k]
                  + f_160 * hh_51[k]
                  - f_159 * hh_58[k]
                  + f_160 * hh_60[k]
                  - f_161 * hh_62[k]
                  + f_127 * hh_191[k]
                  + f_128 * hh_196[k]
                  - f_129 * hh_198[k]
                  + f_127 * hh_205[k]
                  - f_129 * hh_207[k]
                  + f_130 * hh_209[k]
                  + f_159 * hh_338[k]
                  + f_127 * hh_343[k]
                  - f_160 * hh_345[k]
                  + f_159 * hh_352[k]
                  - f_160 * hh_354[k]
                  + f_161 * hh_356[k]
                  - f_127 * hh_380[k]
                  - f_128 * hh_385[k]
                  + f_129 * hh_387[k]
                  - f_127 * hh_394[k]
                  + f_129 * hh_396[k]
                  - f_130 * hh_398[k];
    }

#pragma omp simd aligned(hh_42, hh_45, hh_47, hh_52, hh_54, hh_56, hh_189, hh_192, hh_194, \
                         hh_199, hh_201, hh_203, hh_336, hh_339, hh_341, hh_346, hh_348, \
                         hh_350, hh_378, hh_381, hh_383, hh_388, hh_390, \
                         hh_392 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = -f_152 * hh_42[k]
                  - f_120 * hh_45[k]
                  + f_153 * hh_47[k]
                  - f_152 * hh_52[k]
                  + f_153 * hh_54[k]
                  - f_124 * hh_56[k]
                  + f_120 * hh_189[k]
                  + f_121 * hh_192[k]
                  - f_122 * hh_194[k]
                  + f_120 * hh_199[k]
                  - f_122 * hh_201[k]
                  + f_123 * hh_203[k]
                  + f_152 * hh_336[k]
                  + f_120 * hh_339[k]
                  - f_153 * hh_341[k]
                  + f_152 * hh_346[k]
                  - f_153 * hh_348[k]
                  + f_124 * hh_350[k]
                  - f_120 * hh_378[k]
                  - f_121 * hh_381[k]
                  + f_122 * hh_383[k]
                  - f_120 * hh_388[k]
                  + f_122 * hh_390[k]
                  - f_123 * hh_392[k];
    }

#pragma omp simd aligned(hh_44, hh_51, hh_58, hh_60, hh_191, hh_198, hh_205, hh_207, hh_338, \
                         hh_345, hh_352, hh_354, hh_380, hh_387, hh_394, \
                         hh_396 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = 6.5625 * hh_44[k]
                  - 13.125 * hh_51[k]
                  - 6.5625 * hh_58[k]
                  + 13.125 * hh_60[k]
                  - 13.125 * hh_191[k]
                  + 26.25 * hh_198[k]
                  + 13.125 * hh_205[k]
                  - 26.25 * hh_207[k]
                  - 6.5625 * hh_338[k]
                  + 13.125 * hh_345[k]
                  + 6.5625 * hh_352[k]
                  - 13.125 * hh_354[k]
                  + 13.125 * hh_380[k]
                  - 26.25 * hh_387[k]
                  - 13.125 * hh_394[k]
                  + 26.25 * hh_396[k];
    }

#pragma omp simd aligned(hh_42, hh_45, hh_47, hh_52, hh_54, hh_189, hh_192, hh_194, hh_199, \
                         hh_201, hh_336, hh_339, hh_341, hh_346, hh_348, hh_378, hh_381, \
                         hh_383, hh_388, hh_390 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = f_111 * hh_42[k]
                  - f_73 * hh_45[k]
                  - f_70 * hh_47[k]
                  - f_109 * hh_52[k]
                  + f_110 * hh_54[k]
                  - f_73 * hh_189[k]
                  + f_69 * hh_192[k]
                  + f_74 * hh_194[k]
                  + f_67 * hh_199[k]
                  - f_71 * hh_201[k]
                  - f_111 * hh_336[k]
                  + f_73 * hh_339[k]
                  + f_70 * hh_341[k]
                  + f_109 * hh_346[k]
                  - f_110 * hh_348[k]
                  + f_73 * hh_378[k]
                  - f_69 * hh_381[k]
                  - f_74 * hh_383[k]
                  - f_67 * hh_388[k]
                  + f_71 * hh_390[k];
    }

#pragma omp simd aligned(hh_44, hh_49, hh_58, hh_191, hh_196, hh_205, hh_338, hh_343, hh_352, \
                         hh_380, hh_385, hh_394 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_170 * hh_44[k]
                  + f_171 * hh_49[k]
                  - f_170 * hh_58[k]
                  + f_134 * hh_191[k]
                  - f_135 * hh_196[k]
                  + f_134 * hh_205[k]
                  + f_170 * hh_338[k]
                  - f_171 * hh_343[k]
                  + f_170 * hh_352[k]
                  - f_134 * hh_380[k]
                  + f_135 * hh_385[k]
                  - f_134 * hh_394[k];
    }

#pragma omp simd aligned(hh_42, hh_45, hh_52, hh_189, hh_192, hh_199, hh_336, hh_339, hh_346, \
                         hh_378, hh_381, hh_388 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = -f_44 * hh_42[k]
                  + f_17 * hh_45[k]
                  - f_43 * hh_52[k]
                  + f_20 * hh_189[k]
                  - f_18 * hh_192[k]
                  + f_17 * hh_199[k]
                  + f_44 * hh_336[k]
                  - f_17 * hh_339[k]
                  + f_43 * hh_346[k]
                  - f_20 * hh_378[k]
                  + f_18 * hh_381[k]
                  - f_17 * hh_388[k];
    }

#pragma omp simd aligned(hh_1, hh_6, hh_15, hh_64, hh_69, hh_78, hh_106, hh_111, hh_120, \
                         hh_211, hh_216, hh_225, hh_253, hh_258, \
                         hh_267 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = -f_6 * hh_1[k]
                  + f_4 * hh_6[k]
                  - f_15 * hh_15[k]
                  + f_4 * hh_64[k]
                  - f_9 * hh_69[k]
                  + f_13 * hh_78[k]
                  + f_7 * hh_106[k]
                  - f_11 * hh_111[k]
                  + f_16 * hh_120[k]
                  + f_3 * hh_211[k]
                  - f_8 * hh_216[k]
                  + f_12 * hh_225[k]
                  - f_5 * hh_253[k]
                  + f_10 * hh_258[k]
                  - f_14 * hh_267[k];
    }

#pragma omp simd aligned(hh_4, hh_11, hh_67, hh_74, hh_109, hh_116, hh_214, hh_221, hh_256, \
                         hh_263 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = -f_54 * hh_4[k]
                  + f_54 * hh_11[k]
                  + f_52 * hh_67[k]
                  - f_52 * hh_74[k]
                  + f_55 * hh_109[k]
                  - f_55 * hh_116[k]
                  + f_51 * hh_214[k]
                  - f_51 * hh_221[k]
                  - f_53 * hh_256[k]
                  + f_53 * hh_263[k];
    }

#pragma omp simd aligned(hh_1, hh_6, hh_8, hh_15, hh_17, hh_64, hh_69, hh_71, hh_78, hh_80, \
                         hh_106, hh_111, hh_113, hh_120, hh_122, hh_211, hh_216, hh_218, \
                         hh_225, hh_227, hh_253, hh_258, hh_260, hh_267, \
                         hh_269 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = 0.8203125 * hh_1[k]
                  + 0.546875 * hh_6[k]
                  - 6.5625 * hh_8[k]
                  - 0.2734375 * hh_15[k]
                  + 2.1875 * hh_17[k]
                  - 1.640625 * hh_64[k]
                  - 1.09375 * hh_69[k]
                  + 13.125 * hh_71[k]
                  + 0.546875 * hh_78[k]
                  - 4.375 * hh_80[k]
                  - 6.5625 * hh_106[k]
                  - 4.375 * hh_111[k]
                  + 52.5 * hh_113[k]
                  + 2.1875 * hh_120[k]
                  - 17.5 * hh_122[k]
                  - 2.4609375 * hh_211[k]
                  - 1.640625 * hh_216[k]
                  + 19.6875 * hh_218[k]
                  + 0.8203125 * hh_225[k]
                  - 6.5625 * hh_227[k]
                  + 19.6875 * hh_253[k]
                  + 13.125 * hh_258[k]
                  - 157.5 * hh_260[k]
                  - 6.5625 * hh_267[k]
                  + 52.5 * hh_269[k];
    }

#pragma omp simd aligned(hh_4, hh_11, hh_13, hh_67, hh_74, hh_76, hh_109, hh_116, hh_118, \
                         hh_214, hh_221, hh_223, hh_256, hh_263, \
                         hh_265 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = f_73 * hh_4[k]
                  + f_73 * hh_11[k]
                  - f_69 * hh_13[k]
                  - f_69 * hh_67[k]
                  - f_69 * hh_74[k]
                  + f_70 * hh_76[k]
                  - f_74 * hh_109[k]
                  - f_74 * hh_116[k]
                  + f_75 * hh_118[k]
                  - f_67 * hh_214[k]
                  - f_67 * hh_221[k]
                  + f_68 * hh_223[k]
                  + f_71 * hh_256[k]
                  + f_71 * hh_263[k]
                  - f_72 * hh_265[k];
    }

#pragma omp simd aligned(hh_1, hh_6, hh_8, hh_15, hh_17, hh_19, hh_64, hh_69, hh_71, hh_78, \
                         hh_80, hh_82, hh_106, hh_111, hh_113, hh_120, hh_122, hh_124, hh_211, \
                         hh_216, hh_218, hh_225, hh_227, hh_229, hh_253, hh_258, hh_260, \
                         hh_267, hh_269, hh_271 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = -f_86 * hh_1[k]
                  - f_80 * hh_6[k]
                  + f_87 * hh_8[k]
                  - f_86 * hh_15[k]
                  + f_87 * hh_17[k]
                  - f_88 * hh_19[k]
                  + f_80 * hh_64[k]
                  + f_81 * hh_69[k]
                  - f_79 * hh_71[k]
                  + f_80 * hh_78[k]
                  - f_79 * hh_80[k]
                  + f_82 * hh_82[k]
                  + f_88 * hh_106[k]
                  + f_82 * hh_111[k]
                  - f_89 * hh_113[k]
                  + f_88 * hh_120[k]
                  - f_89 * hh_122[k]
                  + f_90 * hh_124[k]
                  + f_76 * hh_211[k]
                  + f_77 * hh_216[k]
                  - f_78 * hh_218[k]
                  + f_76 * hh_225[k]
                  - f_78 * hh_227[k]
                  + f_79 * hh_229[k]
                  - f_79 * hh_253[k]
                  - f_83 * hh_258[k]
                  + f_84 * hh_260[k]
                  - f_79 * hh_267[k]
                  + f_84 * hh_269[k]
                  - f_85 * hh_271[k];
    }

#pragma omp simd aligned(hh_2, hh_7, hh_9, hh_16, hh_18, hh_20, hh_65, hh_70, hh_72, hh_79, \
                         hh_81, hh_83, hh_107, hh_112, hh_114, hh_121, hh_123, hh_125, hh_212, \
                         hh_217, hh_219, hh_226, hh_228, hh_230, hh_254, hh_259, hh_261, \
                         hh_268, hh_270, hh_272 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = -f_103 * hh_2[k]
                  - f_95 * hh_7[k]
                  + f_104 * hh_9[k]
                  - f_103 * hh_16[k]
                  + f_104 * hh_18[k]
                  - f_105 * hh_20[k]
                  + f_95 * hh_65[k]
                  + f_96 * hh_70[k]
                  - f_97 * hh_72[k]
                  + f_95 * hh_79[k]
                  - f_97 * hh_81[k]
                  + f_98 * hh_83[k]
                  + f_93 * hh_107[k]
                  + f_106 * hh_112[k]
                  - f_107 * hh_114[k]
                  + f_93 * hh_121[k]
                  - f_107 * hh_123[k]
                  + f_108 * hh_125[k]
                  + f_91 * hh_212[k]
                  + f_92 * hh_217[k]
                  - f_93 * hh_219[k]
                  + f_91 * hh_226[k]
                  - f_93 * hh_228[k]
                  + f_94 * hh_230[k]
                  - f_99 * hh_254[k]
                  - f_100 * hh_259[k]
                  + f_101 * hh_261[k]
                  - f_99 * hh_268[k]
                  + f_101 * hh_270[k]
                  - f_102 * hh_272[k];
    }

#pragma omp simd aligned(hh_0, hh_3, hh_5, hh_10, hh_12, hh_14, hh_63, hh_66, hh_68, hh_73, \
                         hh_75, hh_77, hh_105, hh_108, hh_110, hh_115, hh_117, hh_119, hh_210, \
                         hh_213, hh_215, hh_220, hh_222, hh_224, hh_252, hh_255, hh_257, \
                         hh_262, hh_264, hh_266 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = -f_86 * hh_0[k]
                  - f_80 * hh_3[k]
                  + f_87 * hh_5[k]
                  - f_86 * hh_10[k]
                  + f_87 * hh_12[k]
                  - f_88 * hh_14[k]
                  + f_80 * hh_63[k]
                  + f_81 * hh_66[k]
                  - f_79 * hh_68[k]
                  + f_80 * hh_73[k]
                  - f_79 * hh_75[k]
                  + f_82 * hh_77[k]
                  + f_88 * hh_105[k]
                  + f_82 * hh_108[k]
                  - f_89 * hh_110[k]
                  + f_88 * hh_115[k]
                  - f_89 * hh_117[k]
                  + f_90 * hh_119[k]
                  + f_76 * hh_210[k]
                  + f_77 * hh_213[k]
                  - f_78 * hh_215[k]
                  + f_76 * hh_220[k]
                  - f_78 * hh_222[k]
                  + f_79 * hh_224[k]
                  - f_79 * hh_252[k]
                  - f_83 * hh_255[k]
                  + f_84 * hh_257[k]
                  - f_79 * hh_262[k]
                  + f_84 * hh_264[k]
                  - f_85 * hh_266[k];
    }

#pragma omp simd aligned(hh_2, hh_9, hh_16, hh_18, hh_65, hh_72, hh_79, hh_81, hh_107, hh_114, \
                         hh_121, hh_123, hh_212, hh_219, hh_226, hh_228, hh_254, hh_261, \
                         hh_268, hh_270 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = f_111 * hh_2[k]
                  - f_73 * hh_9[k]
                  - f_111 * hh_16[k]
                  + f_73 * hh_18[k]
                  - f_73 * hh_65[k]
                  + f_69 * hh_72[k]
                  + f_73 * hh_79[k]
                  - f_69 * hh_81[k]
                  - f_70 * hh_107[k]
                  + f_74 * hh_114[k]
                  + f_70 * hh_121[k]
                  - f_74 * hh_123[k]
                  - f_109 * hh_212[k]
                  + f_67 * hh_219[k]
                  + f_109 * hh_226[k]
                  - f_67 * hh_228[k]
                  + f_110 * hh_254[k]
                  - f_71 * hh_261[k]
                  - f_110 * hh_268[k]
                  + f_71 * hh_270[k];
    }

#pragma omp simd aligned(hh_0, hh_3, hh_5, hh_10, hh_12, hh_63, hh_66, hh_68, hh_73, hh_75, \
                         hh_105, hh_108, hh_110, hh_115, hh_117, hh_210, hh_213, hh_215, \
                         hh_220, hh_222, hh_252, hh_255, hh_257, hh_262, \
                         hh_264 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = 0.2734375 * hh_0[k]
                  - 0.546875 * hh_3[k]
                  - 2.1875 * hh_5[k]
                  - 0.8203125 * hh_10[k]
                  + 6.5625 * hh_12[k]
                  - 0.546875 * hh_63[k]
                  + 1.09375 * hh_66[k]
                  + 4.375 * hh_68[k]
                  + 1.640625 * hh_73[k]
                  - 13.125 * hh_75[k]
                  - 2.1875 * hh_105[k]
                  + 4.375 * hh_108[k]
                  + 17.5 * hh_110[k]
                  + 6.5625 * hh_115[k]
                  - 52.5 * hh_117[k]
                  - 0.8203125 * hh_210[k]
                  + 1.640625 * hh_213[k]
                  + 6.5625 * hh_215[k]
                  + 2.4609375 * hh_220[k]
                  - 19.6875 * hh_222[k]
                  + 6.5625 * hh_252[k]
                  - 13.125 * hh_255[k]
                  - 52.5 * hh_257[k]
                  - 19.6875 * hh_262[k]
                  + 157.5 * hh_264[k];
    }

#pragma omp simd aligned(hh_2, hh_7, hh_16, hh_65, hh_70, hh_79, hh_107, hh_112, hh_121, \
                         hh_212, hh_217, hh_226, hh_254, hh_259, \
                         hh_268 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = -f_117 * hh_2[k]
                  + f_118 * hh_7[k]
                  - f_117 * hh_16[k]
                  + f_114 * hh_65[k]
                  - f_51 * hh_70[k]
                  + f_114 * hh_79[k]
                  + f_52 * hh_107[k]
                  - f_119 * hh_112[k]
                  + f_52 * hh_121[k]
                  + f_112 * hh_212[k]
                  - f_113 * hh_217[k]
                  + f_112 * hh_226[k]
                  - f_115 * hh_254[k]
                  + f_116 * hh_259[k]
                  - f_115 * hh_268[k];
    }

#pragma omp simd aligned(hh_0, hh_3, hh_10, hh_63, hh_66, hh_73, hh_105, hh_108, hh_115, \
                         hh_210, hh_213, hh_220, hh_252, hh_255, \
                         hh_262 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = -f_15 * hh_0[k]
                  + f_4 * hh_3[k]
                  - f_6 * hh_10[k]
                  + f_13 * hh_63[k]
                  - f_9 * hh_66[k]
                  + f_4 * hh_73[k]
                  + f_16 * hh_105[k]
                  - f_11 * hh_108[k]
                  + f_7 * hh_115[k]
                  + f_12 * hh_210[k]
                  - f_8 * hh_213[k]
                  + f_3 * hh_220[k]
                  - f_14 * hh_252[k]
                  + f_10 * hh_255[k]
                  - f_5 * hh_262[k];
    }

#pragma omp simd aligned(hh_43, hh_48, hh_57, hh_148, hh_153, hh_162, hh_337, hh_342, \
                         hh_351 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = f_45 * hh_43[k]
                  - f_47 * hh_48[k]
                  + f_49 * hh_57[k]
                  - f_46 * hh_148[k]
                  + f_48 * hh_153[k]
                  - f_50 * hh_162[k]
                  + f_45 * hh_337[k]
                  - f_47 * hh_342[k]
                  + f_49 * hh_351[k];
    }

#pragma omp simd aligned(hh_46, hh_53, hh_151, hh_158, hh_340, hh_347 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = 19.6875 * hh_46[k]
                   - 19.6875 * hh_53[k]
                   - 118.125 * hh_151[k]
                   + 118.125 * hh_158[k]
                   + 19.6875 * hh_340[k]
                   - 19.6875 * hh_347[k];
    }

#pragma omp simd aligned(hh_43, hh_48, hh_50, hh_57, hh_59, hh_148, hh_153, hh_155, hh_162, \
                         hh_164, hh_337, hh_342, hh_344, hh_351, \
                         hh_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = -f_112 * hh_43[k]
                   - f_114 * hh_48[k]
                   + f_115 * hh_50[k]
                   + f_117 * hh_57[k]
                   - f_52 * hh_59[k]
                   + f_113 * hh_148[k]
                   + f_51 * hh_153[k]
                   - f_116 * hh_155[k]
                   - f_118 * hh_162[k]
                   + f_119 * hh_164[k]
                   - f_112 * hh_337[k]
                   - f_114 * hh_342[k]
                   + f_115 * hh_344[k]
                   + f_117 * hh_351[k]
                   - f_52 * hh_353[k];
    }

#pragma omp simd aligned(hh_46, hh_53, hh_55, hh_151, hh_158, hh_160, hh_340, hh_347, \
                         hh_349 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = -f_134 * hh_46[k]
                   - f_134 * hh_53[k]
                   + f_66 * hh_55[k]
                   + f_135 * hh_151[k]
                   + f_135 * hh_158[k]
                   - f_136 * hh_160[k]
                   - f_134 * hh_340[k]
                   - f_134 * hh_347[k]
                   + f_66 * hh_349[k];
    }

#pragma omp simd aligned(hh_43, hh_48, hh_50, hh_57, hh_59, hh_61, hh_148, hh_153, hh_155, \
                         hh_162, hh_164, hh_166, hh_337, hh_342, hh_344, hh_351, hh_353, \
                         hh_355 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = f_154 * hh_43[k]
                   + f_156 * hh_48[k]
                   - f_157 * hh_50[k]
                   + f_154 * hh_57[k]
                   - f_157 * hh_59[k]
                   + f_59 * hh_61[k]
                   - f_155 * hh_148[k]
                   - f_157 * hh_153[k]
                   + f_158 * hh_155[k]
                   - f_155 * hh_162[k]
                   + f_158 * hh_164[k]
                   - f_60 * hh_166[k]
                   + f_154 * hh_337[k]
                   + f_156 * hh_342[k]
                   - f_157 * hh_344[k]
                   + f_154 * hh_351[k]
                   - f_157 * hh_353[k]
                   + f_59 * hh_355[k];
    }

#pragma omp simd aligned(hh_44, hh_49, hh_51, hh_58, hh_60, hh_62, hh_149, hh_154, hh_156, \
                         hh_163, hh_165, hh_167, hh_338, hh_343, hh_345, hh_352, hh_354, \
                         hh_356 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = f_162 * hh_44[k]
                   + f_164 * hh_49[k]
                   - f_166 * hh_51[k]
                   + f_162 * hh_58[k]
                   - f_166 * hh_60[k]
                   + f_168 * hh_62[k]
                   - f_163 * hh_149[k]
                   - f_165 * hh_154[k]
                   + f_167 * hh_156[k]
                   - f_163 * hh_163[k]
                   + f_167 * hh_165[k]
                   - f_169 * hh_167[k]
                   + f_162 * hh_338[k]
                   + f_164 * hh_343[k]
                   - f_166 * hh_345[k]
                   + f_162 * hh_352[k]
                   - f_166 * hh_354[k]
                   + f_168 * hh_356[k];
    }

#pragma omp simd aligned(hh_42, hh_45, hh_47, hh_52, hh_54, hh_56, hh_147, hh_150, hh_152, \
                         hh_157, hh_159, hh_161, hh_336, hh_339, hh_341, hh_346, hh_348, \
                         hh_350 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = f_154 * hh_42[k]
                   + f_156 * hh_45[k]
                   - f_157 * hh_47[k]
                   + f_154 * hh_52[k]
                   - f_157 * hh_54[k]
                   + f_59 * hh_56[k]
                   - f_155 * hh_147[k]
                   - f_157 * hh_150[k]
                   + f_158 * hh_152[k]
                   - f_155 * hh_157[k]
                   + f_158 * hh_159[k]
                   - f_60 * hh_161[k]
                   + f_154 * hh_336[k]
                   + f_156 * hh_339[k]
                   - f_157 * hh_341[k]
                   + f_154 * hh_346[k]
                   - f_157 * hh_348[k]
                   + f_59 * hh_350[k];
    }

#pragma omp simd aligned(hh_44, hh_51, hh_58, hh_60, hh_149, hh_156, hh_163, hh_165, hh_338, \
                         hh_345, hh_352, hh_354 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_106[k] = -f_170 * hh_44[k]
                   + f_134 * hh_51[k]
                   + f_170 * hh_58[k]
                   - f_134 * hh_60[k]
                   + f_171 * hh_149[k]
                   - f_135 * hh_156[k]
                   - f_171 * hh_163[k]
                   + f_135 * hh_165[k]
                   - f_170 * hh_338[k]
                   + f_134 * hh_345[k]
                   + f_170 * hh_352[k]
                   - f_134 * hh_354[k];
    }

#pragma omp simd aligned(hh_42, hh_45, hh_47, hh_52, hh_54, hh_147, hh_150, hh_152, hh_157, \
                         hh_159, hh_336, hh_339, hh_341, hh_346, \
                         hh_348 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = -f_117 * hh_42[k]
                   + f_114 * hh_45[k]
                   + f_52 * hh_47[k]
                   + f_112 * hh_52[k]
                   - f_115 * hh_54[k]
                   + f_118 * hh_147[k]
                   - f_51 * hh_150[k]
                   - f_119 * hh_152[k]
                   - f_113 * hh_157[k]
                   + f_116 * hh_159[k]
                   - f_117 * hh_336[k]
                   + f_114 * hh_339[k]
                   + f_52 * hh_341[k]
                   + f_112 * hh_346[k]
                   - f_115 * hh_348[k];
    }

#pragma omp simd aligned(hh_44, hh_49, hh_58, hh_149, hh_154, hh_163, hh_338, hh_343, \
                         hh_352 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = 4.921875 * hh_44[k]
                   - 29.53125 * hh_49[k]
                   + 4.921875 * hh_58[k]
                   - 29.53125 * hh_149[k]
                   + 177.1875 * hh_154[k]
                   - 29.53125 * hh_163[k]
                   + 4.921875 * hh_338[k]
                   - 29.53125 * hh_343[k]
                   + 4.921875 * hh_352[k];
    }

#pragma omp simd aligned(hh_42, hh_45, hh_52, hh_147, hh_150, hh_157, hh_336, hh_339, \
                         hh_346 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = f_49 * hh_42[k]
                   - f_47 * hh_45[k]
                   + f_45 * hh_52[k]
                   - f_50 * hh_147[k]
                   + f_48 * hh_150[k]
                   - f_46 * hh_157[k]
                   + f_49 * hh_336[k]
                   - f_47 * hh_339[k]
                   + f_45 * hh_346[k];
    }

#pragma omp simd aligned(hh_1, hh_6, hh_15, hh_64, hh_69, hh_78, hh_211, hh_216, \
                         hh_225 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = 2.4609375 * hh_1[k]
                   - 4.921875 * hh_6[k]
                   + 0.4921875 * hh_15[k]
                   - 24.609375 * hh_64[k]
                   + 49.21875 * hh_69[k]
                   - 4.921875 * hh_78[k]
                   + 12.3046875 * hh_211[k]
                   - 24.609375 * hh_216[k]
                   + 2.4609375 * hh_225[k];
    }

#pragma omp simd aligned(hh_4, hh_11, hh_67, hh_74, hh_214, hh_221 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = f_2 * hh_4[k]
                   - f_2 * hh_11[k]
                   - f_1 * hh_67[k]
                   + f_1 * hh_74[k]
                   + f_0 * hh_214[k]
                   - f_0 * hh_221[k];
    }

#pragma omp simd aligned(hh_1, hh_6, hh_8, hh_15, hh_17, hh_64, hh_69, hh_71, hh_78, hh_80, \
                         hh_211, hh_216, hh_218, hh_225, hh_227 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = -f_12 * hh_1[k]
                   - f_13 * hh_6[k]
                   + f_14 * hh_8[k]
                   + f_15 * hh_15[k]
                   - f_16 * hh_17[k]
                   + f_8 * hh_64[k]
                   + f_9 * hh_69[k]
                   - f_10 * hh_71[k]
                   - f_4 * hh_78[k]
                   + f_11 * hh_80[k]
                   - f_3 * hh_211[k]
                   - f_4 * hh_216[k]
                   + f_5 * hh_218[k]
                   + f_6 * hh_225[k]
                   - f_7 * hh_227[k];
    }

#pragma omp simd aligned(hh_4, hh_11, hh_13, hh_67, hh_74, hh_76, hh_214, hh_221, \
                         hh_223 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = -f_20 * hh_4[k]
                   - f_20 * hh_11[k]
                   + f_21 * hh_13[k]
                   + f_18 * hh_67[k]
                   + f_18 * hh_74[k]
                   - f_19 * hh_76[k]
                   - f_17 * hh_214[k]
                   - f_17 * hh_221[k]
                   + f_18 * hh_223[k];
    }

#pragma omp simd aligned(hh_1, hh_6, hh_8, hh_15, hh_17, hh_19, hh_64, hh_69, hh_71, hh_78, \
                         hh_80, hh_82, hh_211, hh_216, hh_218, hh_225, hh_227, \
                         hh_229 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = f_29 * hh_1[k]
                   + f_30 * hh_6[k]
                   - f_31 * hh_8[k]
                   + f_29 * hh_15[k]
                   - f_31 * hh_17[k]
                   + f_32 * hh_19[k]
                   - f_23 * hh_64[k]
                   - f_26 * hh_69[k]
                   + f_27 * hh_71[k]
                   - f_23 * hh_78[k]
                   + f_27 * hh_80[k]
                   - f_28 * hh_82[k]
                   + f_22 * hh_211[k]
                   + f_23 * hh_216[k]
                   - f_24 * hh_218[k]
                   + f_22 * hh_225[k]
                   - f_24 * hh_227[k]
                   + f_25 * hh_229[k];
    }

#pragma omp simd aligned(hh_2, hh_7, hh_9, hh_16, hh_18, hh_20, hh_65, hh_70, hh_72, hh_79, \
                         hh_81, hh_83, hh_212, hh_217, hh_219, hh_226, hh_228, \
                         hh_230 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = f_40 * hh_2[k]
                   + f_41 * hh_7[k]
                   - f_36 * hh_9[k]
                   + f_40 * hh_16[k]
                   - f_36 * hh_18[k]
                   + f_42 * hh_20[k]
                   - f_34 * hh_65[k]
                   - f_37 * hh_70[k]
                   + f_38 * hh_72[k]
                   - f_34 * hh_79[k]
                   + f_38 * hh_81[k]
                   - f_39 * hh_83[k]
                   + f_33 * hh_212[k]
                   + f_34 * hh_217[k]
                   - f_35 * hh_219[k]
                   + f_33 * hh_226[k]
                   - f_35 * hh_228[k]
                   + f_36 * hh_230[k];
    }

#pragma omp simd aligned(hh_0, hh_3, hh_5, hh_10, hh_12, hh_14, hh_63, hh_66, hh_68, hh_73, \
                         hh_75, hh_77, hh_210, hh_213, hh_215, hh_220, hh_222, \
                         hh_224 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = f_29 * hh_0[k]
                   + f_30 * hh_3[k]
                   - f_31 * hh_5[k]
                   + f_29 * hh_10[k]
                   - f_31 * hh_12[k]
                   + f_32 * hh_14[k]
                   - f_23 * hh_63[k]
                   - f_26 * hh_66[k]
                   + f_27 * hh_68[k]
                   - f_23 * hh_73[k]
                   + f_27 * hh_75[k]
                   - f_28 * hh_77[k]
                   + f_22 * hh_210[k]
                   + f_23 * hh_213[k]
                   - f_24 * hh_215[k]
                   + f_22 * hh_220[k]
                   - f_24 * hh_222[k]
                   + f_25 * hh_224[k];
    }

#pragma omp simd aligned(hh_2, hh_9, hh_16, hh_18, hh_65, hh_72, hh_79, hh_81, hh_212, hh_219, \
                         hh_226, hh_228 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_117[k] = -f_44 * hh_2[k]
                   + f_20 * hh_9[k]
                   + f_44 * hh_16[k]
                   - f_20 * hh_18[k]
                   + f_17 * hh_65[k]
                   - f_18 * hh_72[k]
                   - f_17 * hh_79[k]
                   + f_18 * hh_81[k]
                   - f_43 * hh_212[k]
                   + f_17 * hh_219[k]
                   + f_43 * hh_226[k]
                   - f_17 * hh_228[k];
    }

#pragma omp simd aligned(hh_0, hh_3, hh_5, hh_10, hh_12, hh_63, hh_66, hh_68, hh_73, hh_75, \
                         hh_210, hh_213, hh_215, hh_220, hh_222 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_118[k] = -f_15 * hh_0[k]
                   + f_13 * hh_3[k]
                   + f_16 * hh_5[k]
                   + f_12 * hh_10[k]
                   - f_14 * hh_12[k]
                   + f_4 * hh_63[k]
                   - f_9 * hh_66[k]
                   - f_11 * hh_68[k]
                   - f_8 * hh_73[k]
                   + f_10 * hh_75[k]
                   - f_6 * hh_210[k]
                   + f_4 * hh_213[k]
                   + f_7 * hh_215[k]
                   + f_3 * hh_220[k]
                   - f_5 * hh_222[k];
    }

#pragma omp simd aligned(hh_2, hh_7, hh_16, hh_65, hh_70, hh_79, hh_212, hh_217, \
                         hh_226 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_119[k] = f_49 * hh_2[k]
                   - f_50 * hh_7[k]
                   + f_49 * hh_16[k]
                   - f_47 * hh_65[k]
                   + f_48 * hh_70[k]
                   - f_47 * hh_79[k]
                   + f_45 * hh_212[k]
                   - f_46 * hh_217[k]
                   + f_45 * hh_226[k];
    }

#pragma omp simd aligned(hh_0, hh_3, hh_10, hh_63, hh_66, hh_73, hh_210, hh_213, \
                         hh_220 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_120[k] = 0.4921875 * hh_0[k]
                   - 4.921875 * hh_3[k]
                   + 2.4609375 * hh_10[k]
                   - 4.921875 * hh_63[k]
                   + 49.21875 * hh_66[k]
                   - 24.609375 * hh_73[k]
                   + 2.4609375 * hh_210[k]
                   - 24.609375 * hh_213[k]
                   + 12.3046875 * hh_220[k];
    }
}

auto
transform_hh_tri(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t hh,
                 const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 9.84375 * std::sqrt(10.0);
    const auto f_1 = 19.6875 * std::sqrt(10.0);
    const auto f_2 = 1.96875 * std::sqrt(10.0);
    const auto f_3 = 2.4609375 * std::sqrt(5.0);
    const auto f_4 = 1.640625 * std::sqrt(5.0);
    const auto f_5 = 19.6875 * std::sqrt(5.0);
    const auto f_6 = 0.8203125 * std::sqrt(5.0);
    const auto f_7 = 6.5625 * std::sqrt(5.0);
    const auto f_8 = 4.921875 * std::sqrt(5.0);
    const auto f_9 = 3.28125 * std::sqrt(5.0);
    const auto f_10 = 39.375 * std::sqrt(5.0);
    const auto f_11 = 13.125 * std::sqrt(5.0);
    const auto f_12 = 0.4921875 * std::sqrt(5.0);
    const auto f_13 = 0.328125 * std::sqrt(5.0);
    const auto f_14 = 3.9375 * std::sqrt(5.0);
    const auto f_15 = 0.1640625 * std::sqrt(5.0);
    const auto f_16 = 1.3125 * std::sqrt(5.0);
    const auto f_17 = 3.28125 * std::sqrt(30.0);
    const auto f_18 = 6.5625 * std::sqrt(30.0);
    const auto f_19 = 13.125 * std::sqrt(30.0);
    const auto f_20 = 0.65625 * std::sqrt(30.0);
    const auto f_21 = 1.3125 * std::sqrt(30.0);
    const auto f_22 = 0.1171875 * std::sqrt(210.0);
    const auto f_23 = 0.234375 * std::sqrt(210.0);
    const auto f_24 = 1.40625 * std::sqrt(210.0);
    const auto f_25 = 0.9375 * std::sqrt(210.0);
    const auto f_26 = 0.46875 * std::sqrt(210.0);
    const auto f_27 = 2.8125 * std::sqrt(210.0);
    const auto f_28 = 1.875 * std::sqrt(210.0);
    const auto f_29 = 0.0234375 * std::sqrt(210.0);
    const auto f_30 = 0.046875 * std::sqrt(210.0);
    const auto f_31 = 0.28125 * std::sqrt(210.0);
    const auto f_32 = 0.1875 * std::sqrt(210.0);
    const auto f_33 = 1.7578125 * std::sqrt(14.0);
    const auto f_34 = 3.515625 * std::sqrt(14.0);
    const auto f_35 = 4.6875 * std::sqrt(14.0);
    const auto f_36 = 0.9375 * std::sqrt(14.0);
    const auto f_37 = 7.03125 * std::sqrt(14.0);
    const auto f_38 = 9.375 * std::sqrt(14.0);
    const auto f_39 = 1.875 * std::sqrt(14.0);
    const auto f_40 = 0.3515625 * std::sqrt(14.0);
    const auto f_41 = 0.703125 * std::sqrt(14.0);
    const auto f_42 = 0.1875 * std::sqrt(14.0);
    const auto f_43 = 1.640625 * std::sqrt(30.0);
    const auto f_44 = 0.328125 * std::sqrt(30.0);
    const auto f_45 = 2.4609375 * std::sqrt(10.0);
    const auto f_46 = 14.765625 * std::sqrt(10.0);
    const auto f_47 = 4.921875 * std::sqrt(10.0);
    const auto f_48 = 29.53125 * std::sqrt(10.0);
    const auto f_49 = 0.4921875 * std::sqrt(10.0);
    const auto f_50 = 2.953125 * std::sqrt(10.0);
    const auto f_51 = 9.84375 * std::sqrt(2.0);
    const auto f_52 = 6.5625 * std::sqrt(2.0);
    const auto f_53 = 78.75 * std::sqrt(2.0);
    const auto f_54 = 3.28125 * std::sqrt(2.0);
    const auto f_55 = 26.25 * std::sqrt(2.0);
    const auto f_56 = 26.25 * std::sqrt(3.0);
    const auto f_57 = 52.5 * std::sqrt(3.0);
    const auto f_58 = 0.9375 * std::sqrt(21.0);
    const auto f_59 = 1.875 * std::sqrt(21.0);
    const auto f_60 = 11.25 * std::sqrt(21.0);
    const auto f_61 = 7.5 * std::sqrt(21.0);
    const auto f_62 = 2.8125 * std::sqrt(35.0);
    const auto f_63 = 5.625 * std::sqrt(35.0);
    const auto f_64 = 7.5 * std::sqrt(35.0);
    const auto f_65 = 1.5 * std::sqrt(35.0);
    const auto f_66 = 13.125 * std::sqrt(3.0);
    const auto f_67 = 3.28125 * std::sqrt(6.0);
    const auto f_68 = 6.5625 * std::sqrt(6.0);
    const auto f_69 = 2.1875 * std::sqrt(6.0);
    const auto f_70 = 4.375 * std::sqrt(6.0);
    const auto f_71 = 26.25 * std::sqrt(6.0);
    const auto f_72 = 52.5 * std::sqrt(6.0);
    const auto f_73 = 1.09375 * std::sqrt(6.0);
    const auto f_74 = 8.75 * std::sqrt(6.0);
    const auto f_75 = 17.5 * std::sqrt(6.0);
    const auto f_76 = 0.1171875 * std::sqrt(42.0);
    const auto f_77 = 0.234375 * std::sqrt(42.0);
    const auto f_78 = 1.40625 * std::sqrt(42.0);
    const auto f_79 = 0.9375 * std::sqrt(42.0);
    const auto f_80 = 0.078125 * std::sqrt(42.0);
    const auto f_81 = 0.15625 * std::sqrt(42.0);
    const auto f_82 = 0.625 * std::sqrt(42.0);
    const auto f_83 = 1.875 * std::sqrt(42.0);
    const auto f_84 = 11.25 * std::sqrt(42.0);
    const auto f_85 = 7.5 * std::sqrt(42.0);
    const auto f_86 = 0.0390625 * std::sqrt(42.0);
    const auto f_87 = 0.46875 * std::sqrt(42.0);
    const auto f_88 = 0.3125 * std::sqrt(42.0);
    const auto f_89 = 3.75 * std::sqrt(42.0);
    const auto f_90 = 2.5 * std::sqrt(42.0);
    const auto f_91 = 0.3515625 * std::sqrt(70.0);
    const auto f_92 = 0.703125 * std::sqrt(70.0);
    const auto f_93 = 0.9375 * std::sqrt(70.0);
    const auto f_94 = 0.1875 * std::sqrt(70.0);
    const auto f_95 = 0.234375 * std::sqrt(70.0);
    const auto f_96 = 0.46875 * std::sqrt(70.0);
    const auto f_97 = 0.625 * std::sqrt(70.0);
    const auto f_98 = 0.125 * std::sqrt(70.0);
    const auto f_99 = 2.8125 * std::sqrt(70.0);
    const auto f_100 = 5.625 * std::sqrt(70.0);
    const auto f_101 = 7.5 * std::sqrt(70.0);
    const auto f_102 = 1.5 * std::sqrt(70.0);
    const auto f_103 = 0.1171875 * std::sqrt(70.0);
    const auto f_104 = 0.3125 * std::sqrt(70.0);
    const auto f_105 = 0.0625 * std::sqrt(70.0);
    const auto f_106 = 1.875 * std::sqrt(70.0);
    const auto f_107 = 2.5 * std::sqrt(70.0);
    const auto f_108 = 0.5 * std::sqrt(70.0);
    const auto f_109 = 1.640625 * std::sqrt(6.0);
    const auto f_110 = 13.125 * std::sqrt(6.0);
    const auto f_111 = 0.546875 * std::sqrt(6.0);
    const auto f_112 = 2.4609375 * std::sqrt(2.0);
    const auto f_113 = 14.765625 * std::sqrt(2.0);
    const auto f_114 = 1.640625 * std::sqrt(2.0);
    const auto f_115 = 19.6875 * std::sqrt(2.0);
    const auto f_116 = 118.125 * std::sqrt(2.0);
    const auto f_117 = 0.8203125 * std::sqrt(2.0);
    const auto f_118 = 4.921875 * std::sqrt(2.0);
    const auto f_119 = 39.375 * std::sqrt(2.0);
    const auto f_120 = 0.9375 * std::sqrt(7.0);
    const auto f_121 = 1.875 * std::sqrt(7.0);
    const auto f_122 = 11.25 * std::sqrt(7.0);
    const auto f_123 = 7.5 * std::sqrt(7.0);
    const auto f_124 = 3.75 * std::sqrt(7.0);
    const auto f_125 = 22.5 * std::sqrt(7.0);
    const auto f_126 = 15.0 * std::sqrt(7.0);
    const auto f_127 = 0.9375 * std::sqrt(105.0);
    const auto f_128 = 1.875 * std::sqrt(105.0);
    const auto f_129 = 2.5 * std::sqrt(105.0);
    const auto f_130 = 0.5 * std::sqrt(105.0);
    const auto f_131 = 3.75 * std::sqrt(105.0);
    const auto f_132 = 5.0 * std::sqrt(105.0);
    const auto f_133 = std::sqrt(105.0);
    const auto f_134 = 6.5625 * std::sqrt(3.0);
    const auto f_135 = 39.375 * std::sqrt(3.0);
    const auto f_136 = 78.75 * std::sqrt(3.0);
    const auto f_137 = 0.234375 * std::sqrt(15.0);
    const auto f_138 = 0.46875 * std::sqrt(15.0);
    const auto f_139 = 0.625 * std::sqrt(15.0);
    const auto f_140 = 0.125 * std::sqrt(15.0);
    const auto f_141 = 0.9375 * std::sqrt(15.0);
    const auto f_142 = 1.25 * std::sqrt(15.0);
    const auto f_143 = 0.25 * std::sqrt(15.0);
    const auto f_144 = 2.8125 * std::sqrt(15.0);
    const auto f_145 = 5.625 * std::sqrt(15.0);
    const auto f_146 = 7.5 * std::sqrt(15.0);
    const auto f_147 = 1.5 * std::sqrt(15.0);
    const auto f_148 = 1.875 * std::sqrt(15.0);
    const auto f_149 = 3.75 * std::sqrt(15.0);
    const auto f_150 = 5.0 * std::sqrt(15.0);
    const auto f_151 = std::sqrt(15.0);
    const auto f_152 = 0.46875 * std::sqrt(7.0);
    const auto f_153 = 5.625 * std::sqrt(7.0);
    const auto f_154 = 0.234375 * std::sqrt(21.0);
    const auto f_155 = 1.40625 * std::sqrt(21.0);
    const auto f_156 = 0.46875 * std::sqrt(21.0);
    const auto f_157 = 2.8125 * std::sqrt(21.0);
    const auto f_158 = 16.875 * std::sqrt(21.0);
    const auto f_159 = 0.46875 * std::sqrt(105.0);
    const auto f_160 = 1.25 * std::sqrt(105.0);
    const auto f_161 = 0.25 * std::sqrt(105.0);
    const auto f_162 = 0.703125 * std::sqrt(35.0);
    const auto f_163 = 4.21875 * std::sqrt(35.0);
    const auto f_164 = 1.40625 * std::sqrt(35.0);
    const auto f_165 = 8.4375 * std::sqrt(35.0);
    const auto f_166 = 1.875 * std::sqrt(35.0);
    const auto f_167 = 11.25 * std::sqrt(35.0);
    const auto f_168 = 0.375 * std::sqrt(35.0);
    const auto f_169 = 2.25 * std::sqrt(35.0);
    const auto f_170 = 3.28125 * std::sqrt(3.0);
    const auto f_171 = 19.6875 * std::sqrt(3.0);

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
    auto *g_85 = values + 85 * nvalues;
    auto *g_86 = values + 86 * nvalues;
    auto *g_87 = values + 87 * nvalues;
    auto *g_88 = values + 88 * nvalues;
    auto *g_89 = values + 89 * nvalues;
    auto *g_90 = values + 90 * nvalues;
    auto *g_91 = values + 91 * nvalues;
    auto *g_92 = values + 92 * nvalues;
    auto *g_93 = values + 93 * nvalues;
    auto *g_94 = values + 94 * nvalues;
    auto *g_95 = values + 95 * nvalues;
    auto *g_96 = values + 96 * nvalues;
    auto *g_97 = values + 97 * nvalues;
    auto *g_98 = values + 98 * nvalues;
    auto *g_99 = values + 99 * nvalues;
    auto *g_100 = values + 100 * nvalues;
    auto *g_101 = values + 101 * nvalues;
    auto *g_102 = values + 102 * nvalues;
    auto *g_103 = values + 103 * nvalues;
    auto *g_104 = values + 104 * nvalues;
    auto *g_105 = values + 105 * nvalues;
    auto *g_106 = values + 106 * nvalues;
    auto *g_107 = values + 107 * nvalues;
    auto *g_108 = values + 108 * nvalues;
    auto *g_109 = values + 109 * nvalues;
    auto *g_110 = values + 110 * nvalues;
    auto *g_111 = values + 111 * nvalues;
    auto *g_112 = values + 112 * nvalues;
    auto *g_113 = values + 113 * nvalues;
    auto *g_114 = values + 114 * nvalues;
    auto *g_115 = values + 115 * nvalues;
    auto *g_116 = values + 116 * nvalues;
    auto *g_117 = values + 117 * nvalues;
    auto *g_118 = values + 118 * nvalues;
    auto *g_119 = values + 119 * nvalues;
    auto *g_120 = values + 120 * nvalues;

    const auto *hh_0 = buffer.data(hh + 0);
    const auto *hh_2 = buffer.data(hh + 2);
    const auto *hh_3 = buffer.data(hh + 3);
    const auto *hh_5 = buffer.data(hh + 5);
    const auto *hh_7 = buffer.data(hh + 7);
    const auto *hh_9 = buffer.data(hh + 9);
    const auto *hh_10 = buffer.data(hh + 10);
    const auto *hh_12 = buffer.data(hh + 12);
    const auto *hh_14 = buffer.data(hh + 14);
    const auto *hh_16 = buffer.data(hh + 16);
    const auto *hh_18 = buffer.data(hh + 18);
    const auto *hh_21 = buffer.data(hh + 21);
    const auto *hh_22 = buffer.data(hh + 22);
    const auto *hh_23 = buffer.data(hh + 23);
    const auto *hh_24 = buffer.data(hh + 24);
    const auto *hh_25 = buffer.data(hh + 25);
    const auto *hh_26 = buffer.data(hh + 26);
    const auto *hh_27 = buffer.data(hh + 27);
    const auto *hh_28 = buffer.data(hh + 28);
    const auto *hh_29 = buffer.data(hh + 29);
    const auto *hh_30 = buffer.data(hh + 30);
    const auto *hh_31 = buffer.data(hh + 31);
    const auto *hh_32 = buffer.data(hh + 32);
    const auto *hh_33 = buffer.data(hh + 33);
    const auto *hh_34 = buffer.data(hh + 34);
    const auto *hh_35 = buffer.data(hh + 35);
    const auto *hh_36 = buffer.data(hh + 36);
    const auto *hh_37 = buffer.data(hh + 37);
    const auto *hh_38 = buffer.data(hh + 38);
    const auto *hh_39 = buffer.data(hh + 39);
    const auto *hh_40 = buffer.data(hh + 40);
    const auto *hh_41 = buffer.data(hh + 41);
    const auto *hh_42 = buffer.data(hh + 42);
    const auto *hh_44 = buffer.data(hh + 44);
    const auto *hh_45 = buffer.data(hh + 45);
    const auto *hh_47 = buffer.data(hh + 47);
    const auto *hh_49 = buffer.data(hh + 49);
    const auto *hh_51 = buffer.data(hh + 51);
    const auto *hh_52 = buffer.data(hh + 52);
    const auto *hh_54 = buffer.data(hh + 54);
    const auto *hh_56 = buffer.data(hh + 56);
    const auto *hh_58 = buffer.data(hh + 58);
    const auto *hh_60 = buffer.data(hh + 60);
    const auto *hh_62 = buffer.data(hh + 62);
    const auto *hh_63 = buffer.data(hh + 63);
    const auto *hh_65 = buffer.data(hh + 65);
    const auto *hh_66 = buffer.data(hh + 66);
    const auto *hh_68 = buffer.data(hh + 68);
    const auto *hh_70 = buffer.data(hh + 70);
    const auto *hh_72 = buffer.data(hh + 72);
    const auto *hh_73 = buffer.data(hh + 73);
    const auto *hh_75 = buffer.data(hh + 75);
    const auto *hh_77 = buffer.data(hh + 77);
    const auto *hh_79 = buffer.data(hh + 79);
    const auto *hh_81 = buffer.data(hh + 81);
    const auto *hh_84 = buffer.data(hh + 84);
    const auto *hh_85 = buffer.data(hh + 85);
    const auto *hh_86 = buffer.data(hh + 86);
    const auto *hh_87 = buffer.data(hh + 87);
    const auto *hh_88 = buffer.data(hh + 88);
    const auto *hh_89 = buffer.data(hh + 89);
    const auto *hh_90 = buffer.data(hh + 90);
    const auto *hh_91 = buffer.data(hh + 91);
    const auto *hh_92 = buffer.data(hh + 92);
    const auto *hh_93 = buffer.data(hh + 93);
    const auto *hh_94 = buffer.data(hh + 94);
    const auto *hh_95 = buffer.data(hh + 95);
    const auto *hh_96 = buffer.data(hh + 96);
    const auto *hh_97 = buffer.data(hh + 97);
    const auto *hh_98 = buffer.data(hh + 98);
    const auto *hh_99 = buffer.data(hh + 99);
    const auto *hh_100 = buffer.data(hh + 100);
    const auto *hh_101 = buffer.data(hh + 101);
    const auto *hh_102 = buffer.data(hh + 102);
    const auto *hh_103 = buffer.data(hh + 103);
    const auto *hh_104 = buffer.data(hh + 104);
    const auto *hh_105 = buffer.data(hh + 105);
    const auto *hh_107 = buffer.data(hh + 107);
    const auto *hh_108 = buffer.data(hh + 108);
    const auto *hh_110 = buffer.data(hh + 110);
    const auto *hh_112 = buffer.data(hh + 112);
    const auto *hh_114 = buffer.data(hh + 114);
    const auto *hh_115 = buffer.data(hh + 115);
    const auto *hh_117 = buffer.data(hh + 117);
    const auto *hh_119 = buffer.data(hh + 119);
    const auto *hh_121 = buffer.data(hh + 121);
    const auto *hh_123 = buffer.data(hh + 123);
    const auto *hh_126 = buffer.data(hh + 126);
    const auto *hh_127 = buffer.data(hh + 127);
    const auto *hh_128 = buffer.data(hh + 128);
    const auto *hh_129 = buffer.data(hh + 129);
    const auto *hh_130 = buffer.data(hh + 130);
    const auto *hh_131 = buffer.data(hh + 131);
    const auto *hh_132 = buffer.data(hh + 132);
    const auto *hh_133 = buffer.data(hh + 133);
    const auto *hh_134 = buffer.data(hh + 134);
    const auto *hh_135 = buffer.data(hh + 135);
    const auto *hh_136 = buffer.data(hh + 136);
    const auto *hh_137 = buffer.data(hh + 137);
    const auto *hh_138 = buffer.data(hh + 138);
    const auto *hh_139 = buffer.data(hh + 139);
    const auto *hh_140 = buffer.data(hh + 140);
    const auto *hh_141 = buffer.data(hh + 141);
    const auto *hh_142 = buffer.data(hh + 142);
    const auto *hh_143 = buffer.data(hh + 143);
    const auto *hh_144 = buffer.data(hh + 144);
    const auto *hh_145 = buffer.data(hh + 145);
    const auto *hh_146 = buffer.data(hh + 146);
    const auto *hh_147 = buffer.data(hh + 147);
    const auto *hh_149 = buffer.data(hh + 149);
    const auto *hh_150 = buffer.data(hh + 150);
    const auto *hh_152 = buffer.data(hh + 152);
    const auto *hh_154 = buffer.data(hh + 154);
    const auto *hh_156 = buffer.data(hh + 156);
    const auto *hh_157 = buffer.data(hh + 157);
    const auto *hh_159 = buffer.data(hh + 159);
    const auto *hh_161 = buffer.data(hh + 161);
    const auto *hh_163 = buffer.data(hh + 163);
    const auto *hh_165 = buffer.data(hh + 165);
    const auto *hh_167 = buffer.data(hh + 167);
    const auto *hh_168 = buffer.data(hh + 168);
    const auto *hh_169 = buffer.data(hh + 169);
    const auto *hh_170 = buffer.data(hh + 170);
    const auto *hh_171 = buffer.data(hh + 171);
    const auto *hh_172 = buffer.data(hh + 172);
    const auto *hh_173 = buffer.data(hh + 173);
    const auto *hh_174 = buffer.data(hh + 174);
    const auto *hh_175 = buffer.data(hh + 175);
    const auto *hh_176 = buffer.data(hh + 176);
    const auto *hh_177 = buffer.data(hh + 177);
    const auto *hh_178 = buffer.data(hh + 178);
    const auto *hh_179 = buffer.data(hh + 179);
    const auto *hh_180 = buffer.data(hh + 180);
    const auto *hh_181 = buffer.data(hh + 181);
    const auto *hh_182 = buffer.data(hh + 182);
    const auto *hh_183 = buffer.data(hh + 183);
    const auto *hh_184 = buffer.data(hh + 184);
    const auto *hh_185 = buffer.data(hh + 185);
    const auto *hh_186 = buffer.data(hh + 186);
    const auto *hh_187 = buffer.data(hh + 187);
    const auto *hh_188 = buffer.data(hh + 188);
    const auto *hh_189 = buffer.data(hh + 189);
    const auto *hh_191 = buffer.data(hh + 191);
    const auto *hh_192 = buffer.data(hh + 192);
    const auto *hh_194 = buffer.data(hh + 194);
    const auto *hh_196 = buffer.data(hh + 196);
    const auto *hh_198 = buffer.data(hh + 198);
    const auto *hh_199 = buffer.data(hh + 199);
    const auto *hh_201 = buffer.data(hh + 201);
    const auto *hh_203 = buffer.data(hh + 203);
    const auto *hh_205 = buffer.data(hh + 205);
    const auto *hh_207 = buffer.data(hh + 207);
    const auto *hh_209 = buffer.data(hh + 209);
    const auto *hh_210 = buffer.data(hh + 210);
    const auto *hh_212 = buffer.data(hh + 212);
    const auto *hh_213 = buffer.data(hh + 213);
    const auto *hh_215 = buffer.data(hh + 215);
    const auto *hh_217 = buffer.data(hh + 217);
    const auto *hh_219 = buffer.data(hh + 219);
    const auto *hh_220 = buffer.data(hh + 220);
    const auto *hh_222 = buffer.data(hh + 222);
    const auto *hh_224 = buffer.data(hh + 224);
    const auto *hh_226 = buffer.data(hh + 226);
    const auto *hh_228 = buffer.data(hh + 228);
    const auto *hh_231 = buffer.data(hh + 231);
    const auto *hh_232 = buffer.data(hh + 232);
    const auto *hh_233 = buffer.data(hh + 233);
    const auto *hh_234 = buffer.data(hh + 234);
    const auto *hh_235 = buffer.data(hh + 235);
    const auto *hh_236 = buffer.data(hh + 236);
    const auto *hh_237 = buffer.data(hh + 237);
    const auto *hh_238 = buffer.data(hh + 238);
    const auto *hh_239 = buffer.data(hh + 239);
    const auto *hh_240 = buffer.data(hh + 240);
    const auto *hh_241 = buffer.data(hh + 241);
    const auto *hh_242 = buffer.data(hh + 242);
    const auto *hh_243 = buffer.data(hh + 243);
    const auto *hh_244 = buffer.data(hh + 244);
    const auto *hh_245 = buffer.data(hh + 245);
    const auto *hh_246 = buffer.data(hh + 246);
    const auto *hh_247 = buffer.data(hh + 247);
    const auto *hh_248 = buffer.data(hh + 248);
    const auto *hh_249 = buffer.data(hh + 249);
    const auto *hh_250 = buffer.data(hh + 250);
    const auto *hh_251 = buffer.data(hh + 251);
    const auto *hh_252 = buffer.data(hh + 252);
    const auto *hh_254 = buffer.data(hh + 254);
    const auto *hh_255 = buffer.data(hh + 255);
    const auto *hh_257 = buffer.data(hh + 257);
    const auto *hh_259 = buffer.data(hh + 259);
    const auto *hh_261 = buffer.data(hh + 261);
    const auto *hh_262 = buffer.data(hh + 262);
    const auto *hh_264 = buffer.data(hh + 264);
    const auto *hh_266 = buffer.data(hh + 266);
    const auto *hh_268 = buffer.data(hh + 268);
    const auto *hh_270 = buffer.data(hh + 270);
    const auto *hh_273 = buffer.data(hh + 273);
    const auto *hh_274 = buffer.data(hh + 274);
    const auto *hh_275 = buffer.data(hh + 275);
    const auto *hh_276 = buffer.data(hh + 276);
    const auto *hh_277 = buffer.data(hh + 277);
    const auto *hh_278 = buffer.data(hh + 278);
    const auto *hh_279 = buffer.data(hh + 279);
    const auto *hh_280 = buffer.data(hh + 280);
    const auto *hh_281 = buffer.data(hh + 281);
    const auto *hh_282 = buffer.data(hh + 282);
    const auto *hh_283 = buffer.data(hh + 283);
    const auto *hh_284 = buffer.data(hh + 284);
    const auto *hh_285 = buffer.data(hh + 285);
    const auto *hh_286 = buffer.data(hh + 286);
    const auto *hh_287 = buffer.data(hh + 287);
    const auto *hh_288 = buffer.data(hh + 288);
    const auto *hh_289 = buffer.data(hh + 289);
    const auto *hh_290 = buffer.data(hh + 290);
    const auto *hh_291 = buffer.data(hh + 291);
    const auto *hh_292 = buffer.data(hh + 292);
    const auto *hh_293 = buffer.data(hh + 293);
    const auto *hh_294 = buffer.data(hh + 294);
    const auto *hh_296 = buffer.data(hh + 296);
    const auto *hh_297 = buffer.data(hh + 297);
    const auto *hh_299 = buffer.data(hh + 299);
    const auto *hh_301 = buffer.data(hh + 301);
    const auto *hh_303 = buffer.data(hh + 303);
    const auto *hh_304 = buffer.data(hh + 304);
    const auto *hh_306 = buffer.data(hh + 306);
    const auto *hh_308 = buffer.data(hh + 308);
    const auto *hh_310 = buffer.data(hh + 310);
    const auto *hh_312 = buffer.data(hh + 312);
    const auto *hh_315 = buffer.data(hh + 315);
    const auto *hh_316 = buffer.data(hh + 316);
    const auto *hh_317 = buffer.data(hh + 317);
    const auto *hh_318 = buffer.data(hh + 318);
    const auto *hh_319 = buffer.data(hh + 319);
    const auto *hh_320 = buffer.data(hh + 320);
    const auto *hh_321 = buffer.data(hh + 321);
    const auto *hh_322 = buffer.data(hh + 322);
    const auto *hh_323 = buffer.data(hh + 323);
    const auto *hh_324 = buffer.data(hh + 324);
    const auto *hh_325 = buffer.data(hh + 325);
    const auto *hh_326 = buffer.data(hh + 326);
    const auto *hh_327 = buffer.data(hh + 327);
    const auto *hh_328 = buffer.data(hh + 328);
    const auto *hh_329 = buffer.data(hh + 329);
    const auto *hh_330 = buffer.data(hh + 330);
    const auto *hh_331 = buffer.data(hh + 331);
    const auto *hh_332 = buffer.data(hh + 332);
    const auto *hh_333 = buffer.data(hh + 333);
    const auto *hh_334 = buffer.data(hh + 334);
    const auto *hh_335 = buffer.data(hh + 335);
    const auto *hh_336 = buffer.data(hh + 336);
    const auto *hh_338 = buffer.data(hh + 338);
    const auto *hh_339 = buffer.data(hh + 339);
    const auto *hh_341 = buffer.data(hh + 341);
    const auto *hh_343 = buffer.data(hh + 343);
    const auto *hh_345 = buffer.data(hh + 345);
    const auto *hh_346 = buffer.data(hh + 346);
    const auto *hh_348 = buffer.data(hh + 348);
    const auto *hh_350 = buffer.data(hh + 350);
    const auto *hh_352 = buffer.data(hh + 352);
    const auto *hh_354 = buffer.data(hh + 354);
    const auto *hh_356 = buffer.data(hh + 356);
    const auto *hh_357 = buffer.data(hh + 357);
    const auto *hh_358 = buffer.data(hh + 358);
    const auto *hh_359 = buffer.data(hh + 359);
    const auto *hh_360 = buffer.data(hh + 360);
    const auto *hh_361 = buffer.data(hh + 361);
    const auto *hh_362 = buffer.data(hh + 362);
    const auto *hh_363 = buffer.data(hh + 363);
    const auto *hh_364 = buffer.data(hh + 364);
    const auto *hh_365 = buffer.data(hh + 365);
    const auto *hh_366 = buffer.data(hh + 366);
    const auto *hh_367 = buffer.data(hh + 367);
    const auto *hh_368 = buffer.data(hh + 368);
    const auto *hh_369 = buffer.data(hh + 369);
    const auto *hh_370 = buffer.data(hh + 370);
    const auto *hh_371 = buffer.data(hh + 371);
    const auto *hh_372 = buffer.data(hh + 372);
    const auto *hh_373 = buffer.data(hh + 373);
    const auto *hh_374 = buffer.data(hh + 374);
    const auto *hh_375 = buffer.data(hh + 375);
    const auto *hh_376 = buffer.data(hh + 376);
    const auto *hh_377 = buffer.data(hh + 377);
    const auto *hh_378 = buffer.data(hh + 378);
    const auto *hh_380 = buffer.data(hh + 380);
    const auto *hh_381 = buffer.data(hh + 381);
    const auto *hh_383 = buffer.data(hh + 383);
    const auto *hh_385 = buffer.data(hh + 385);
    const auto *hh_387 = buffer.data(hh + 387);
    const auto *hh_388 = buffer.data(hh + 388);
    const auto *hh_390 = buffer.data(hh + 390);
    const auto *hh_392 = buffer.data(hh + 392);
    const auto *hh_394 = buffer.data(hh + 394);
    const auto *hh_396 = buffer.data(hh + 396);
    const auto *hh_398 = buffer.data(hh + 398);
    const auto *hh_399 = buffer.data(hh + 399);
    const auto *hh_400 = buffer.data(hh + 400);
    const auto *hh_401 = buffer.data(hh + 401);
    const auto *hh_402 = buffer.data(hh + 402);
    const auto *hh_404 = buffer.data(hh + 404);
    const auto *hh_405 = buffer.data(hh + 405);
    const auto *hh_406 = buffer.data(hh + 406);
    const auto *hh_407 = buffer.data(hh + 407);
    const auto *hh_408 = buffer.data(hh + 408);
    const auto *hh_409 = buffer.data(hh + 409);
    const auto *hh_411 = buffer.data(hh + 411);
    const auto *hh_413 = buffer.data(hh + 413);
    const auto *hh_414 = buffer.data(hh + 414);
    const auto *hh_415 = buffer.data(hh + 415);
    const auto *hh_416 = buffer.data(hh + 416);
    const auto *hh_417 = buffer.data(hh + 417);
    const auto *hh_418 = buffer.data(hh + 418);
    const auto *hh_419 = buffer.data(hh + 419);
    const auto *hh_420 = buffer.data(hh + 420);
    const auto *hh_422 = buffer.data(hh + 422);
    const auto *hh_423 = buffer.data(hh + 423);
    const auto *hh_425 = buffer.data(hh + 425);
    const auto *hh_427 = buffer.data(hh + 427);
    const auto *hh_429 = buffer.data(hh + 429);
    const auto *hh_430 = buffer.data(hh + 430);
    const auto *hh_432 = buffer.data(hh + 432);
    const auto *hh_434 = buffer.data(hh + 434);
    const auto *hh_436 = buffer.data(hh + 436);
    const auto *hh_438 = buffer.data(hh + 438);
    const auto *hh_440 = buffer.data(hh + 440);

#pragma omp simd aligned(hh_22, hh_27, hh_36, hh_127, hh_132, hh_141, hh_316, hh_321, \
                         hh_330 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = 12.3046875 * hh_22[k]
                 - 24.609375 * hh_27[k]
                 + 2.4609375 * hh_36[k]
                 - 24.609375 * hh_127[k]
                 + 49.21875 * hh_132[k]
                 - 4.921875 * hh_141[k]
                 + 2.4609375 * hh_316[k]
                 - 4.921875 * hh_321[k]
                 + 0.4921875 * hh_330[k];
    }

#pragma omp simd aligned(hh_25, hh_32, hh_130, hh_137, hh_319, hh_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = f_0 * hh_25[k]
                 - f_0 * hh_32[k]
                 - f_1 * hh_130[k]
                 + f_1 * hh_137[k]
                 + f_2 * hh_319[k]
                 - f_2 * hh_326[k];
        g_11[k] = g_1[k];
    }

#pragma omp simd aligned(hh_22, hh_27, hh_29, hh_36, hh_38, hh_127, hh_132, hh_134, hh_141, \
                         hh_143, hh_316, hh_321, hh_323, hh_330, \
                         hh_332 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_3 * hh_22[k]
                 - f_4 * hh_27[k]
                 + f_5 * hh_29[k]
                 + f_6 * hh_36[k]
                 - f_7 * hh_38[k]
                 + f_8 * hh_127[k]
                 + f_9 * hh_132[k]
                 - f_10 * hh_134[k]
                 - f_4 * hh_141[k]
                 + f_11 * hh_143[k]
                 - f_12 * hh_316[k]
                 - f_13 * hh_321[k]
                 + f_14 * hh_323[k]
                 + f_15 * hh_330[k]
                 - f_16 * hh_332[k];
        g_22[k] = g_2[k];
    }

#pragma omp simd aligned(hh_25, hh_32, hh_34, hh_130, hh_137, hh_139, hh_319, hh_326, \
                         hh_328 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_17 * hh_25[k]
                 - f_17 * hh_32[k]
                 + f_18 * hh_34[k]
                 + f_18 * hh_130[k]
                 + f_18 * hh_137[k]
                 - f_19 * hh_139[k]
                 - f_20 * hh_319[k]
                 - f_20 * hh_326[k]
                 + f_21 * hh_328[k];
        g_33[k] = g_3[k];
    }

#pragma omp simd aligned(hh_22, hh_27, hh_29, hh_36, hh_38, hh_40, hh_127, hh_132, hh_134, \
                         hh_141, hh_143, hh_145, hh_316, hh_321, hh_323, hh_330, hh_332, \
                         hh_334 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_22 * hh_22[k]
                 + f_23 * hh_27[k]
                 - f_24 * hh_29[k]
                 + f_22 * hh_36[k]
                 - f_24 * hh_38[k]
                 + f_25 * hh_40[k]
                 - f_23 * hh_127[k]
                 - f_26 * hh_132[k]
                 + f_27 * hh_134[k]
                 - f_23 * hh_141[k]
                 + f_27 * hh_143[k]
                 - f_28 * hh_145[k]
                 + f_29 * hh_316[k]
                 + f_30 * hh_321[k]
                 - f_31 * hh_323[k]
                 + f_29 * hh_330[k]
                 - f_31 * hh_332[k]
                 + f_32 * hh_334[k];
        g_44[k] = g_4[k];
    }

#pragma omp simd aligned(hh_23, hh_28, hh_30, hh_37, hh_39, hh_41, hh_128, hh_133, hh_135, \
                         hh_142, hh_144, hh_146, hh_317, hh_322, hh_324, hh_331, hh_333, \
                         hh_335 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_33 * hh_23[k]
                 + f_34 * hh_28[k]
                 - f_35 * hh_30[k]
                 + f_33 * hh_37[k]
                 - f_35 * hh_39[k]
                 + f_36 * hh_41[k]
                 - f_34 * hh_128[k]
                 - f_37 * hh_133[k]
                 + f_38 * hh_135[k]
                 - f_34 * hh_142[k]
                 + f_38 * hh_144[k]
                 - f_39 * hh_146[k]
                 + f_40 * hh_317[k]
                 + f_41 * hh_322[k]
                 - f_36 * hh_324[k]
                 + f_40 * hh_331[k]
                 - f_36 * hh_333[k]
                 + f_42 * hh_335[k];
        g_55[k] = g_5[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_26, hh_31, hh_33, hh_35, hh_126, hh_129, hh_131, \
                         hh_136, hh_138, hh_140, hh_315, hh_318, hh_320, hh_325, hh_327, \
                         hh_329 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_22 * hh_21[k]
                 + f_23 * hh_24[k]
                 - f_24 * hh_26[k]
                 + f_22 * hh_31[k]
                 - f_24 * hh_33[k]
                 + f_25 * hh_35[k]
                 - f_23 * hh_126[k]
                 - f_26 * hh_129[k]
                 + f_27 * hh_131[k]
                 - f_23 * hh_136[k]
                 + f_27 * hh_138[k]
                 - f_28 * hh_140[k]
                 + f_29 * hh_315[k]
                 + f_30 * hh_318[k]
                 - f_31 * hh_320[k]
                 + f_29 * hh_325[k]
                 - f_31 * hh_327[k]
                 + f_32 * hh_329[k];
        g_66[k] = g_6[k];
    }

#pragma omp simd aligned(hh_23, hh_30, hh_37, hh_39, hh_128, hh_135, hh_142, hh_144, hh_317, \
                         hh_324, hh_331, hh_333 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_43 * hh_23[k]
                 + f_17 * hh_30[k]
                 + f_43 * hh_37[k]
                 - f_17 * hh_39[k]
                 + f_17 * hh_128[k]
                 - f_18 * hh_135[k]
                 - f_17 * hh_142[k]
                 + f_18 * hh_144[k]
                 - f_44 * hh_317[k]
                 + f_20 * hh_324[k]
                 + f_44 * hh_331[k]
                 - f_20 * hh_333[k];
        g_77[k] = g_7[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_26, hh_31, hh_33, hh_126, hh_129, hh_131, hh_136, \
                         hh_138, hh_315, hh_318, hh_320, hh_325, \
                         hh_327 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_6 * hh_21[k]
                 + f_4 * hh_24[k]
                 + f_7 * hh_26[k]
                 + f_3 * hh_31[k]
                 - f_5 * hh_33[k]
                 + f_4 * hh_126[k]
                 - f_9 * hh_129[k]
                 - f_11 * hh_131[k]
                 - f_8 * hh_136[k]
                 + f_10 * hh_138[k]
                 - f_15 * hh_315[k]
                 + f_13 * hh_318[k]
                 + f_16 * hh_320[k]
                 + f_12 * hh_325[k]
                 - f_14 * hh_327[k];
        g_88[k] = g_8[k];
    }

#pragma omp simd aligned(hh_23, hh_28, hh_37, hh_128, hh_133, hh_142, hh_317, hh_322, \
                         hh_331 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_45 * hh_23[k]
                 - f_46 * hh_28[k]
                 + f_45 * hh_37[k]
                 - f_47 * hh_128[k]
                 + f_48 * hh_133[k]
                 - f_47 * hh_142[k]
                 + f_49 * hh_317[k]
                 - f_50 * hh_322[k]
                 + f_49 * hh_331[k];
        g_99[k] = g_9[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_31, hh_88, hh_95, hh_126, hh_129, hh_136, hh_235, \
                         hh_242, hh_315, hh_318, hh_325 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = 2.4609375 * hh_21[k]
                  - 24.609375 * hh_24[k]
                  + 12.3046875 * hh_31[k]
                  - 4.921875 * hh_126[k]
                  + 49.21875 * hh_129[k]
                  - 24.609375 * hh_136[k]
                  + 0.4921875 * hh_315[k]
                  - 4.921875 * hh_318[k]
                  + 2.4609375 * hh_325[k];
        g_110[k] = g_10[k];

        g_12[k] = 78.75 * hh_88[k]
                  - 78.75 * hh_95[k]
                  - 78.75 * hh_235[k]
                  + 78.75 * hh_242[k];
    }

#pragma omp simd aligned(hh_85, hh_90, hh_92, hh_99, hh_101, hh_232, hh_237, hh_239, hh_246, \
                         hh_248 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_51 * hh_85[k]
                  - f_52 * hh_90[k]
                  + f_53 * hh_92[k]
                  + f_54 * hh_99[k]
                  - f_55 * hh_101[k]
                  + f_51 * hh_232[k]
                  + f_52 * hh_237[k]
                  - f_53 * hh_239[k]
                  - f_54 * hh_246[k]
                  + f_55 * hh_248[k];
        g_23[k] = g_13[k];
    }

#pragma omp simd aligned(hh_88, hh_95, hh_97, hh_235, hh_242, hh_244 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_56 * hh_88[k]
                  - f_56 * hh_95[k]
                  + f_57 * hh_97[k]
                  + f_56 * hh_235[k]
                  + f_56 * hh_242[k]
                  - f_57 * hh_244[k];
        g_34[k] = g_14[k];
    }

#pragma omp simd aligned(hh_85, hh_90, hh_92, hh_99, hh_101, hh_103, hh_232, hh_237, hh_239, \
                         hh_246, hh_248, hh_250 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_58 * hh_85[k]
                  + f_59 * hh_90[k]
                  - f_60 * hh_92[k]
                  + f_58 * hh_99[k]
                  - f_60 * hh_101[k]
                  + f_61 * hh_103[k]
                  - f_58 * hh_232[k]
                  - f_59 * hh_237[k]
                  + f_60 * hh_239[k]
                  - f_58 * hh_246[k]
                  + f_60 * hh_248[k]
                  - f_61 * hh_250[k];
        g_45[k] = g_15[k];
    }

#pragma omp simd aligned(hh_86, hh_91, hh_93, hh_100, hh_102, hh_104, hh_233, hh_238, hh_240, \
                         hh_247, hh_249, hh_251 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_62 * hh_86[k]
                  + f_63 * hh_91[k]
                  - f_64 * hh_93[k]
                  + f_62 * hh_100[k]
                  - f_64 * hh_102[k]
                  + f_65 * hh_104[k]
                  - f_62 * hh_233[k]
                  - f_63 * hh_238[k]
                  + f_64 * hh_240[k]
                  - f_62 * hh_247[k]
                  + f_64 * hh_249[k]
                  - f_65 * hh_251[k];
        g_56[k] = g_16[k];
    }

#pragma omp simd aligned(hh_84, hh_87, hh_89, hh_94, hh_96, hh_98, hh_231, hh_234, hh_236, \
                         hh_241, hh_243, hh_245 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_58 * hh_84[k]
                  + f_59 * hh_87[k]
                  - f_60 * hh_89[k]
                  + f_58 * hh_94[k]
                  - f_60 * hh_96[k]
                  + f_61 * hh_98[k]
                  - f_58 * hh_231[k]
                  - f_59 * hh_234[k]
                  + f_60 * hh_236[k]
                  - f_58 * hh_241[k]
                  + f_60 * hh_243[k]
                  - f_61 * hh_245[k];
        g_67[k] = g_17[k];
    }

#pragma omp simd aligned(hh_86, hh_93, hh_100, hh_102, hh_233, hh_240, hh_247, \
                         hh_249 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_66 * hh_86[k]
                  + f_56 * hh_93[k]
                  + f_66 * hh_100[k]
                  - f_56 * hh_102[k]
                  + f_66 * hh_233[k]
                  - f_56 * hh_240[k]
                  - f_66 * hh_247[k]
                  + f_56 * hh_249[k];
        g_78[k] = g_18[k];
    }

#pragma omp simd aligned(hh_84, hh_87, hh_89, hh_94, hh_96, hh_231, hh_234, hh_236, hh_241, \
                         hh_243 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_54 * hh_84[k]
                  + f_52 * hh_87[k]
                  + f_55 * hh_89[k]
                  + f_51 * hh_94[k]
                  - f_53 * hh_96[k]
                  + f_54 * hh_231[k]
                  - f_52 * hh_234[k]
                  - f_55 * hh_236[k]
                  - f_51 * hh_241[k]
                  + f_53 * hh_243[k];
        g_89[k] = g_19[k];
    }

#pragma omp simd aligned(hh_84, hh_86, hh_87, hh_91, hh_94, hh_100, hh_231, hh_233, hh_234, \
                         hh_238, hh_241, hh_247 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = 19.6875 * hh_86[k]
                  - 118.125 * hh_91[k]
                  + 19.6875 * hh_100[k]
                  - 19.6875 * hh_233[k]
                  + 118.125 * hh_238[k]
                  - 19.6875 * hh_247[k];
        g_100[k] = g_20[k];

        g_21[k] = f_2 * hh_84[k]
                  - f_1 * hh_87[k]
                  + f_0 * hh_94[k]
                  - f_2 * hh_231[k]
                  + f_1 * hh_234[k]
                  - f_0 * hh_241[k];
        g_111[k] = g_21[k];
    }

#pragma omp simd aligned(hh_22, hh_27, hh_29, hh_36, hh_38, hh_127, hh_132, hh_134, hh_141, \
                         hh_143, hh_169, hh_174, hh_176, hh_183, hh_185, hh_316, hh_321, \
                         hh_323, hh_330, hh_332, hh_358, hh_363, hh_365, hh_372, \
                         hh_374 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = 2.4609375 * hh_22[k]
                  + 1.640625 * hh_27[k]
                  - 19.6875 * hh_29[k]
                  - 0.8203125 * hh_36[k]
                  + 6.5625 * hh_38[k]
                  + 1.640625 * hh_127[k]
                  + 1.09375 * hh_132[k]
                  - 13.125 * hh_134[k]
                  - 0.546875 * hh_141[k]
                  + 4.375 * hh_143[k]
                  - 19.6875 * hh_169[k]
                  - 13.125 * hh_174[k]
                  + 157.5 * hh_176[k]
                  + 6.5625 * hh_183[k]
                  - 52.5 * hh_185[k]
                  - 0.8203125 * hh_316[k]
                  - 0.546875 * hh_321[k]
                  + 6.5625 * hh_323[k]
                  + 0.2734375 * hh_330[k]
                  - 2.1875 * hh_332[k]
                  + 6.5625 * hh_358[k]
                  + 4.375 * hh_363[k]
                  - 52.5 * hh_365[k]
                  - 2.1875 * hh_372[k]
                  + 17.5 * hh_374[k];
    }

#pragma omp simd aligned(hh_25, hh_32, hh_34, hh_130, hh_137, hh_139, hh_172, hh_179, hh_181, \
                         hh_319, hh_326, hh_328, hh_361, hh_368, \
                         hh_370 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_67 * hh_25[k]
                  + f_67 * hh_32[k]
                  - f_68 * hh_34[k]
                  + f_69 * hh_130[k]
                  + f_69 * hh_137[k]
                  - f_70 * hh_139[k]
                  - f_71 * hh_172[k]
                  - f_71 * hh_179[k]
                  + f_72 * hh_181[k]
                  - f_73 * hh_319[k]
                  - f_73 * hh_326[k]
                  + f_69 * hh_328[k]
                  + f_74 * hh_361[k]
                  + f_74 * hh_368[k]
                  - f_75 * hh_370[k];
        g_35[k] = g_25[k];
    }

#pragma omp simd aligned(hh_22, hh_27, hh_29, hh_36, hh_38, hh_40, hh_127, hh_132, hh_134, \
                         hh_141, hh_143, hh_145, hh_169, hh_174, hh_176, hh_183, hh_185, \
                         hh_187, hh_316, hh_321, hh_323, hh_330, hh_332, hh_334, hh_358, \
                         hh_363, hh_365, hh_372, hh_374, hh_376 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_76 * hh_22[k]
                  - f_77 * hh_27[k]
                  + f_78 * hh_29[k]
                  - f_76 * hh_36[k]
                  + f_78 * hh_38[k]
                  - f_79 * hh_40[k]
                  - f_80 * hh_127[k]
                  - f_81 * hh_132[k]
                  + f_79 * hh_134[k]
                  - f_80 * hh_141[k]
                  + f_79 * hh_143[k]
                  - f_82 * hh_145[k]
                  + f_79 * hh_169[k]
                  + f_83 * hh_174[k]
                  - f_84 * hh_176[k]
                  + f_79 * hh_183[k]
                  - f_84 * hh_185[k]
                  + f_85 * hh_187[k]
                  + f_86 * hh_316[k]
                  + f_80 * hh_321[k]
                  - f_87 * hh_323[k]
                  + f_86 * hh_330[k]
                  - f_87 * hh_332[k]
                  + f_88 * hh_334[k]
                  - f_88 * hh_358[k]
                  - f_82 * hh_363[k]
                  + f_89 * hh_365[k]
                  - f_88 * hh_372[k]
                  + f_89 * hh_374[k]
                  - f_90 * hh_376[k];
        g_46[k] = g_26[k];
    }

#pragma omp simd aligned(hh_23, hh_28, hh_30, hh_37, hh_39, hh_41, hh_128, hh_133, hh_135, \
                         hh_142, hh_144, hh_146, hh_170, hh_175, hh_177, hh_184, hh_186, \
                         hh_188, hh_317, hh_322, hh_324, hh_331, hh_333, hh_335, hh_359, \
                         hh_364, hh_366, hh_373, hh_375, hh_377 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_91 * hh_23[k]
                  - f_92 * hh_28[k]
                  + f_93 * hh_30[k]
                  - f_91 * hh_37[k]
                  + f_93 * hh_39[k]
                  - f_94 * hh_41[k]
                  - f_95 * hh_128[k]
                  - f_96 * hh_133[k]
                  + f_97 * hh_135[k]
                  - f_95 * hh_142[k]
                  + f_97 * hh_144[k]
                  - f_98 * hh_146[k]
                  + f_99 * hh_170[k]
                  + f_100 * hh_175[k]
                  - f_101 * hh_177[k]
                  + f_99 * hh_184[k]
                  - f_101 * hh_186[k]
                  + f_102 * hh_188[k]
                  + f_103 * hh_317[k]
                  + f_95 * hh_322[k]
                  - f_104 * hh_324[k]
                  + f_103 * hh_331[k]
                  - f_104 * hh_333[k]
                  + f_105 * hh_335[k]
                  - f_93 * hh_359[k]
                  - f_106 * hh_364[k]
                  + f_107 * hh_366[k]
                  - f_93 * hh_373[k]
                  + f_107 * hh_375[k]
                  - f_108 * hh_377[k];
        g_57[k] = g_27[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_26, hh_31, hh_33, hh_35, hh_126, hh_129, hh_131, \
                         hh_136, hh_138, hh_140, hh_168, hh_171, hh_173, hh_178, hh_180, \
                         hh_182, hh_315, hh_318, hh_320, hh_325, hh_327, hh_329, hh_357, \
                         hh_360, hh_362, hh_367, hh_369, hh_371 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_76 * hh_21[k]
                  - f_77 * hh_24[k]
                  + f_78 * hh_26[k]
                  - f_76 * hh_31[k]
                  + f_78 * hh_33[k]
                  - f_79 * hh_35[k]
                  - f_80 * hh_126[k]
                  - f_81 * hh_129[k]
                  + f_79 * hh_131[k]
                  - f_80 * hh_136[k]
                  + f_79 * hh_138[k]
                  - f_82 * hh_140[k]
                  + f_79 * hh_168[k]
                  + f_83 * hh_171[k]
                  - f_84 * hh_173[k]
                  + f_79 * hh_178[k]
                  - f_84 * hh_180[k]
                  + f_85 * hh_182[k]
                  + f_86 * hh_315[k]
                  + f_80 * hh_318[k]
                  - f_87 * hh_320[k]
                  + f_86 * hh_325[k]
                  - f_87 * hh_327[k]
                  + f_88 * hh_329[k]
                  - f_88 * hh_357[k]
                  - f_82 * hh_360[k]
                  + f_89 * hh_362[k]
                  - f_88 * hh_367[k]
                  + f_89 * hh_369[k]
                  - f_90 * hh_371[k];
        g_68[k] = g_28[k];
    }

#pragma omp simd aligned(hh_23, hh_30, hh_37, hh_39, hh_128, hh_135, hh_142, hh_144, hh_170, \
                         hh_177, hh_184, hh_186, hh_317, hh_324, hh_331, hh_333, hh_359, \
                         hh_366, hh_373, hh_375 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_109 * hh_23[k]
                  - f_67 * hh_30[k]
                  - f_109 * hh_37[k]
                  + f_67 * hh_39[k]
                  + f_73 * hh_128[k]
                  - f_69 * hh_135[k]
                  - f_73 * hh_142[k]
                  + f_69 * hh_144[k]
                  - f_110 * hh_170[k]
                  + f_71 * hh_177[k]
                  + f_110 * hh_184[k]
                  - f_71 * hh_186[k]
                  - f_111 * hh_317[k]
                  + f_73 * hh_324[k]
                  + f_111 * hh_331[k]
                  - f_73 * hh_333[k]
                  + f_70 * hh_359[k]
                  - f_74 * hh_366[k]
                  - f_70 * hh_373[k]
                  + f_74 * hh_375[k];
        g_79[k] = g_29[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_26, hh_31, hh_33, hh_126, hh_129, hh_131, hh_136, \
                         hh_138, hh_168, hh_171, hh_173, hh_178, hh_180, hh_315, hh_318, \
                         hh_320, hh_325, hh_327, hh_357, hh_360, hh_362, hh_367, \
                         hh_369 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = 0.8203125 * hh_21[k]
                  - 1.640625 * hh_24[k]
                  - 6.5625 * hh_26[k]
                  - 2.4609375 * hh_31[k]
                  + 19.6875 * hh_33[k]
                  + 0.546875 * hh_126[k]
                  - 1.09375 * hh_129[k]
                  - 4.375 * hh_131[k]
                  - 1.640625 * hh_136[k]
                  + 13.125 * hh_138[k]
                  - 6.5625 * hh_168[k]
                  + 13.125 * hh_171[k]
                  + 52.5 * hh_173[k]
                  + 19.6875 * hh_178[k]
                  - 157.5 * hh_180[k]
                  - 0.2734375 * hh_315[k]
                  + 0.546875 * hh_318[k]
                  + 2.1875 * hh_320[k]
                  + 0.8203125 * hh_325[k]
                  - 6.5625 * hh_327[k]
                  + 2.1875 * hh_357[k]
                  - 4.375 * hh_360[k]
                  - 17.5 * hh_362[k]
                  - 6.5625 * hh_367[k]
                  + 52.5 * hh_369[k];
        g_90[k] = g_30[k];
    }

#pragma omp simd aligned(hh_23, hh_28, hh_37, hh_128, hh_133, hh_142, hh_170, hh_175, hh_184, \
                         hh_317, hh_322, hh_331, hh_359, hh_364, \
                         hh_373 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_112 * hh_23[k]
                  + f_113 * hh_28[k]
                  - f_112 * hh_37[k]
                  - f_114 * hh_128[k]
                  + f_51 * hh_133[k]
                  - f_114 * hh_142[k]
                  + f_115 * hh_170[k]
                  - f_116 * hh_175[k]
                  + f_115 * hh_184[k]
                  + f_117 * hh_317[k]
                  - f_118 * hh_322[k]
                  + f_117 * hh_331[k]
                  - f_52 * hh_359[k]
                  + f_119 * hh_364[k]
                  - f_52 * hh_373[k];
        g_101[k] = g_31[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_31, hh_126, hh_129, hh_136, hh_168, hh_171, hh_178, \
                         hh_315, hh_318, hh_325, hh_357, hh_360, \
                         hh_367 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_12 * hh_21[k]
                  + f_8 * hh_24[k]
                  - f_3 * hh_31[k]
                  - f_13 * hh_126[k]
                  + f_9 * hh_129[k]
                  - f_4 * hh_136[k]
                  + f_14 * hh_168[k]
                  - f_10 * hh_171[k]
                  + f_5 * hh_178[k]
                  + f_15 * hh_315[k]
                  - f_4 * hh_318[k]
                  + f_6 * hh_325[k]
                  - f_16 * hh_357[k]
                  + f_11 * hh_360[k]
                  - f_7 * hh_367[k];
        g_112[k] = g_32[k];
    }

#pragma omp simd aligned(hh_88, hh_95, hh_97, hh_235, hh_242, hh_244, hh_277, hh_284, \
                         hh_286 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = 26.25 * hh_88[k]
                  + 26.25 * hh_95[k]
                  - 52.5 * hh_97[k]
                  + 26.25 * hh_235[k]
                  + 26.25 * hh_242[k]
                  - 52.5 * hh_244[k]
                  - 52.5 * hh_277[k]
                  - 52.5 * hh_284[k]
                  + 105.0 * hh_286[k];
    }

#pragma omp simd aligned(hh_85, hh_90, hh_92, hh_99, hh_101, hh_103, hh_232, hh_237, hh_239, \
                         hh_246, hh_248, hh_250, hh_274, hh_279, hh_281, hh_288, hh_290, \
                         hh_292 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_120 * hh_85[k]
                  - f_121 * hh_90[k]
                  + f_122 * hh_92[k]
                  - f_120 * hh_99[k]
                  + f_122 * hh_101[k]
                  - f_123 * hh_103[k]
                  - f_120 * hh_232[k]
                  - f_121 * hh_237[k]
                  + f_122 * hh_239[k]
                  - f_120 * hh_246[k]
                  + f_122 * hh_248[k]
                  - f_123 * hh_250[k]
                  + f_121 * hh_274[k]
                  + f_124 * hh_279[k]
                  - f_125 * hh_281[k]
                  + f_121 * hh_288[k]
                  - f_125 * hh_290[k]
                  + f_126 * hh_292[k];
        g_47[k] = g_37[k];
    }

#pragma omp simd aligned(hh_86, hh_91, hh_93, hh_100, hh_102, hh_104, hh_233, hh_238, hh_240, \
                         hh_247, hh_249, hh_251, hh_275, hh_280, hh_282, hh_289, hh_291, \
                         hh_293 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_127 * hh_86[k]
                  - f_128 * hh_91[k]
                  + f_129 * hh_93[k]
                  - f_127 * hh_100[k]
                  + f_129 * hh_102[k]
                  - f_130 * hh_104[k]
                  - f_127 * hh_233[k]
                  - f_128 * hh_238[k]
                  + f_129 * hh_240[k]
                  - f_127 * hh_247[k]
                  + f_129 * hh_249[k]
                  - f_130 * hh_251[k]
                  + f_128 * hh_275[k]
                  + f_131 * hh_280[k]
                  - f_132 * hh_282[k]
                  + f_128 * hh_289[k]
                  - f_132 * hh_291[k]
                  + f_133 * hh_293[k];
        g_58[k] = g_38[k];
    }

#pragma omp simd aligned(hh_84, hh_87, hh_89, hh_94, hh_96, hh_98, hh_231, hh_234, hh_236, \
                         hh_241, hh_243, hh_245, hh_273, hh_276, hh_278, hh_283, hh_285, \
                         hh_287 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_120 * hh_84[k]
                  - f_121 * hh_87[k]
                  + f_122 * hh_89[k]
                  - f_120 * hh_94[k]
                  + f_122 * hh_96[k]
                  - f_123 * hh_98[k]
                  - f_120 * hh_231[k]
                  - f_121 * hh_234[k]
                  + f_122 * hh_236[k]
                  - f_120 * hh_241[k]
                  + f_122 * hh_243[k]
                  - f_123 * hh_245[k]
                  + f_121 * hh_273[k]
                  + f_124 * hh_276[k]
                  - f_125 * hh_278[k]
                  + f_121 * hh_283[k]
                  - f_125 * hh_285[k]
                  + f_126 * hh_287[k];
        g_69[k] = g_39[k];
    }

#pragma omp simd aligned(hh_86, hh_93, hh_100, hh_102, hh_233, hh_240, hh_247, hh_249, hh_275, \
                         hh_282, hh_289, hh_291 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = 13.125 * hh_86[k]
                  - 26.25 * hh_93[k]
                  - 13.125 * hh_100[k]
                  + 26.25 * hh_102[k]
                  + 13.125 * hh_233[k]
                  - 26.25 * hh_240[k]
                  - 13.125 * hh_247[k]
                  + 26.25 * hh_249[k]
                  - 26.25 * hh_275[k]
                  + 52.5 * hh_282[k]
                  + 26.25 * hh_289[k]
                  - 52.5 * hh_291[k];
        g_80[k] = g_40[k];
    }

#pragma omp simd aligned(hh_84, hh_87, hh_89, hh_94, hh_96, hh_231, hh_234, hh_236, hh_241, \
                         hh_243, hh_273, hh_276, hh_278, hh_283, \
                         hh_285 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_73 * hh_84[k]
                  - f_69 * hh_87[k]
                  - f_74 * hh_89[k]
                  - f_67 * hh_94[k]
                  + f_71 * hh_96[k]
                  + f_73 * hh_231[k]
                  - f_69 * hh_234[k]
                  - f_74 * hh_236[k]
                  - f_67 * hh_241[k]
                  + f_71 * hh_243[k]
                  - f_69 * hh_273[k]
                  + f_70 * hh_276[k]
                  + f_75 * hh_278[k]
                  + f_68 * hh_283[k]
                  - f_72 * hh_285[k];
        g_91[k] = g_41[k];
    }

#pragma omp simd aligned(hh_86, hh_91, hh_100, hh_233, hh_238, hh_247, hh_275, hh_280, \
                         hh_289 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_134 * hh_86[k]
                  + f_135 * hh_91[k]
                  - f_134 * hh_100[k]
                  - f_134 * hh_233[k]
                  + f_135 * hh_238[k]
                  - f_134 * hh_247[k]
                  + f_66 * hh_275[k]
                  - f_136 * hh_280[k]
                  + f_66 * hh_289[k];
        g_102[k] = g_42[k];
    }

#pragma omp simd aligned(hh_84, hh_87, hh_94, hh_231, hh_234, hh_241, hh_273, hh_276, \
                         hh_283 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_20 * hh_84[k]
                  + f_18 * hh_87[k]
                  - f_17 * hh_94[k]
                  - f_20 * hh_231[k]
                  + f_18 * hh_234[k]
                  - f_17 * hh_241[k]
                  + f_21 * hh_273[k]
                  - f_19 * hh_276[k]
                  + f_18 * hh_283[k];
        g_113[k] = g_43[k];
    }

#pragma omp simd aligned(hh_22, hh_27, hh_29, hh_36, hh_38, hh_40, hh_127, hh_132, hh_134, \
                         hh_141, hh_143, hh_145, hh_169, hh_174, hh_176, hh_183, hh_185, \
                         hh_187, hh_316, hh_321, hh_323, hh_330, hh_332, hh_334, hh_358, \
                         hh_363, hh_365, hh_372, hh_374, hh_376, hh_400, hh_405, hh_407, \
                         hh_414, hh_416, hh_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = 0.234375 * hh_22[k]
                  + 0.46875 * hh_27[k]
                  - 2.8125 * hh_29[k]
                  + 0.234375 * hh_36[k]
                  - 2.8125 * hh_38[k]
                  + 1.875 * hh_40[k]
                  + 0.46875 * hh_127[k]
                  + 0.9375 * hh_132[k]
                  - 5.625 * hh_134[k]
                  + 0.46875 * hh_141[k]
                  - 5.625 * hh_143[k]
                  + 3.75 * hh_145[k]
                  - 2.8125 * hh_169[k]
                  - 5.625 * hh_174[k]
                  + 33.75 * hh_176[k]
                  - 2.8125 * hh_183[k]
                  + 33.75 * hh_185[k]
                  - 22.5 * hh_187[k]
                  + 0.234375 * hh_316[k]
                  + 0.46875 * hh_321[k]
                  - 2.8125 * hh_323[k]
                  + 0.234375 * hh_330[k]
                  - 2.8125 * hh_332[k]
                  + 1.875 * hh_334[k]
                  - 2.8125 * hh_358[k]
                  - 5.625 * hh_363[k]
                  + 33.75 * hh_365[k]
                  - 2.8125 * hh_372[k]
                  + 33.75 * hh_374[k]
                  - 22.5 * hh_376[k]
                  + 1.875 * hh_400[k]
                  + 3.75 * hh_405[k]
                  - 22.5 * hh_407[k]
                  + 1.875 * hh_414[k]
                  - 22.5 * hh_416[k]
                  + 15.0 * hh_418[k];
    }

#pragma omp simd aligned(hh_23, hh_28, hh_30, hh_37, hh_39, hh_41, hh_128, hh_133, hh_135, \
                         hh_142, hh_144, hh_146, hh_170, hh_175, hh_177, hh_184, hh_186, \
                         hh_188, hh_317, hh_322, hh_324, hh_331, hh_333, hh_335, hh_359, \
                         hh_364, hh_366, hh_373, hh_375, hh_377, hh_401, hh_406, hh_408, \
                         hh_415, hh_417, hh_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_137 * hh_23[k]
                  + f_138 * hh_28[k]
                  - f_139 * hh_30[k]
                  + f_137 * hh_37[k]
                  - f_139 * hh_39[k]
                  + f_140 * hh_41[k]
                  + f_138 * hh_128[k]
                  + f_141 * hh_133[k]
                  - f_142 * hh_135[k]
                  + f_138 * hh_142[k]
                  - f_142 * hh_144[k]
                  + f_143 * hh_146[k]
                  - f_144 * hh_170[k]
                  - f_145 * hh_175[k]
                  + f_146 * hh_177[k]
                  - f_144 * hh_184[k]
                  + f_146 * hh_186[k]
                  - f_147 * hh_188[k]
                  + f_137 * hh_317[k]
                  + f_138 * hh_322[k]
                  - f_139 * hh_324[k]
                  + f_137 * hh_331[k]
                  - f_139 * hh_333[k]
                  + f_140 * hh_335[k]
                  - f_144 * hh_359[k]
                  - f_145 * hh_364[k]
                  + f_146 * hh_366[k]
                  - f_144 * hh_373[k]
                  + f_146 * hh_375[k]
                  - f_147 * hh_377[k]
                  + f_148 * hh_401[k]
                  + f_149 * hh_406[k]
                  - f_150 * hh_408[k]
                  + f_148 * hh_415[k]
                  - f_150 * hh_417[k]
                  + f_151 * hh_419[k];
        g_59[k] = g_49[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_26, hh_31, hh_33, hh_35, hh_126, hh_129, hh_131, \
                         hh_136, hh_138, hh_140, hh_168, hh_171, hh_173, hh_178, hh_180, \
                         hh_182, hh_315, hh_318, hh_320, hh_325, hh_327, hh_329, hh_357, \
                         hh_360, hh_362, hh_367, hh_369, hh_371, hh_399, hh_402, hh_404, \
                         hh_409, hh_411, hh_413 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = 0.234375 * hh_21[k]
                  + 0.46875 * hh_24[k]
                  - 2.8125 * hh_26[k]
                  + 0.234375 * hh_31[k]
                  - 2.8125 * hh_33[k]
                  + 1.875 * hh_35[k]
                  + 0.46875 * hh_126[k]
                  + 0.9375 * hh_129[k]
                  - 5.625 * hh_131[k]
                  + 0.46875 * hh_136[k]
                  - 5.625 * hh_138[k]
                  + 3.75 * hh_140[k]
                  - 2.8125 * hh_168[k]
                  - 5.625 * hh_171[k]
                  + 33.75 * hh_173[k]
                  - 2.8125 * hh_178[k]
                  + 33.75 * hh_180[k]
                  - 22.5 * hh_182[k]
                  + 0.234375 * hh_315[k]
                  + 0.46875 * hh_318[k]
                  - 2.8125 * hh_320[k]
                  + 0.234375 * hh_325[k]
                  - 2.8125 * hh_327[k]
                  + 1.875 * hh_329[k]
                  - 2.8125 * hh_357[k]
                  - 5.625 * hh_360[k]
                  + 33.75 * hh_362[k]
                  - 2.8125 * hh_367[k]
                  + 33.75 * hh_369[k]
                  - 22.5 * hh_371[k]
                  + 1.875 * hh_399[k]
                  + 3.75 * hh_402[k]
                  - 22.5 * hh_404[k]
                  + 1.875 * hh_409[k]
                  - 22.5 * hh_411[k]
                  + 15.0 * hh_413[k];
        g_70[k] = g_50[k];
    }

#pragma omp simd aligned(hh_23, hh_30, hh_37, hh_39, hh_128, hh_135, hh_142, hh_144, hh_170, \
                         hh_177, hh_184, hh_186, hh_317, hh_324, hh_331, hh_333, hh_359, \
                         hh_366, hh_373, hh_375, hh_401, hh_408, hh_415, \
                         hh_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_152 * hh_23[k]
                  + f_120 * hh_30[k]
                  + f_152 * hh_37[k]
                  - f_120 * hh_39[k]
                  - f_120 * hh_128[k]
                  + f_121 * hh_135[k]
                  + f_120 * hh_142[k]
                  - f_121 * hh_144[k]
                  + f_153 * hh_170[k]
                  - f_122 * hh_177[k]
                  - f_153 * hh_184[k]
                  + f_122 * hh_186[k]
                  - f_152 * hh_317[k]
                  + f_120 * hh_324[k]
                  + f_152 * hh_331[k]
                  - f_120 * hh_333[k]
                  + f_153 * hh_359[k]
                  - f_122 * hh_366[k]
                  - f_153 * hh_373[k]
                  + f_122 * hh_375[k]
                  - f_124 * hh_401[k]
                  + f_123 * hh_408[k]
                  + f_124 * hh_415[k]
                  - f_123 * hh_417[k];
        g_81[k] = g_51[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_26, hh_31, hh_33, hh_126, hh_129, hh_131, hh_136, \
                         hh_138, hh_168, hh_171, hh_173, hh_178, hh_180, hh_315, hh_318, \
                         hh_320, hh_325, hh_327, hh_357, hh_360, hh_362, hh_367, hh_369, \
                         hh_399, hh_402, hh_404, hh_409, hh_411 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_86 * hh_21[k]
                  + f_80 * hh_24[k]
                  + f_88 * hh_26[k]
                  + f_76 * hh_31[k]
                  - f_79 * hh_33[k]
                  - f_80 * hh_126[k]
                  + f_81 * hh_129[k]
                  + f_82 * hh_131[k]
                  + f_77 * hh_136[k]
                  - f_83 * hh_138[k]
                  + f_87 * hh_168[k]
                  - f_79 * hh_171[k]
                  - f_89 * hh_173[k]
                  - f_78 * hh_178[k]
                  + f_84 * hh_180[k]
                  - f_86 * hh_315[k]
                  + f_80 * hh_318[k]
                  + f_88 * hh_320[k]
                  + f_76 * hh_325[k]
                  - f_79 * hh_327[k]
                  + f_87 * hh_357[k]
                  - f_79 * hh_360[k]
                  - f_89 * hh_362[k]
                  - f_78 * hh_367[k]
                  + f_84 * hh_369[k]
                  - f_88 * hh_399[k]
                  + f_82 * hh_402[k]
                  + f_90 * hh_404[k]
                  + f_79 * hh_409[k]
                  - f_85 * hh_411[k];
        g_92[k] = g_52[k];
    }

#pragma omp simd aligned(hh_23, hh_28, hh_37, hh_128, hh_133, hh_142, hh_170, hh_175, hh_184, \
                         hh_317, hh_322, hh_331, hh_359, hh_364, hh_373, hh_401, hh_406, \
                         hh_415 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_154 * hh_23[k]
                  - f_155 * hh_28[k]
                  + f_154 * hh_37[k]
                  + f_156 * hh_128[k]
                  - f_157 * hh_133[k]
                  + f_156 * hh_142[k]
                  - f_157 * hh_170[k]
                  + f_158 * hh_175[k]
                  - f_157 * hh_184[k]
                  + f_154 * hh_317[k]
                  - f_155 * hh_322[k]
                  + f_154 * hh_331[k]
                  - f_157 * hh_359[k]
                  + f_158 * hh_364[k]
                  - f_157 * hh_373[k]
                  + f_59 * hh_401[k]
                  - f_60 * hh_406[k]
                  + f_59 * hh_415[k];
        g_103[k] = g_53[k];
    }

#pragma omp simd aligned(hh_21, hh_24, hh_31, hh_126, hh_129, hh_136, hh_168, hh_171, hh_178, \
                         hh_315, hh_318, hh_325, hh_357, hh_360, hh_367, hh_399, hh_402, \
                         hh_409 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_29 * hh_21[k]
                  - f_23 * hh_24[k]
                  + f_22 * hh_31[k]
                  + f_30 * hh_126[k]
                  - f_26 * hh_129[k]
                  + f_23 * hh_136[k]
                  - f_31 * hh_168[k]
                  + f_27 * hh_171[k]
                  - f_24 * hh_178[k]
                  + f_29 * hh_315[k]
                  - f_23 * hh_318[k]
                  + f_22 * hh_325[k]
                  - f_31 * hh_357[k]
                  + f_27 * hh_360[k]
                  - f_24 * hh_367[k]
                  + f_32 * hh_399[k]
                  - f_28 * hh_402[k]
                  + f_25 * hh_409[k];
        g_114[k] = g_54[k];
    }

#pragma omp simd aligned(hh_44, hh_49, hh_51, hh_58, hh_60, hh_62, hh_149, hh_154, hh_156, \
                         hh_163, hh_165, hh_167, hh_191, hh_196, hh_198, hh_205, hh_207, \
                         hh_209, hh_338, hh_343, hh_345, hh_352, hh_354, hh_356, hh_380, \
                         hh_385, hh_387, hh_394, hh_396, hh_398, hh_422, hh_427, hh_429, \
                         hh_436, hh_438, hh_440 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = 3.515625 * hh_44[k]
                  + 7.03125 * hh_49[k]
                  - 9.375 * hh_51[k]
                  + 3.515625 * hh_58[k]
                  - 9.375 * hh_60[k]
                  + 1.875 * hh_62[k]
                  + 7.03125 * hh_149[k]
                  + 14.0625 * hh_154[k]
                  - 18.75 * hh_156[k]
                  + 7.03125 * hh_163[k]
                  - 18.75 * hh_165[k]
                  + 3.75 * hh_167[k]
                  - 9.375 * hh_191[k]
                  - 18.75 * hh_196[k]
                  + 25.0 * hh_198[k]
                  - 9.375 * hh_205[k]
                  + 25.0 * hh_207[k]
                  - 5.0 * hh_209[k]
                  + 3.515625 * hh_338[k]
                  + 7.03125 * hh_343[k]
                  - 9.375 * hh_345[k]
                  + 3.515625 * hh_352[k]
                  - 9.375 * hh_354[k]
                  + 1.875 * hh_356[k]
                  - 9.375 * hh_380[k]
                  - 18.75 * hh_385[k]
                  + 25.0 * hh_387[k]
                  - 9.375 * hh_394[k]
                  + 25.0 * hh_396[k]
                  - 5.0 * hh_398[k]
                  + 1.875 * hh_422[k]
                  + 3.75 * hh_427[k]
                  - 5.0 * hh_429[k]
                  + 1.875 * hh_436[k]
                  - 5.0 * hh_438[k]
                  + hh_440[k];
    }

#pragma omp simd aligned(hh_42, hh_45, hh_47, hh_52, hh_54, hh_56, hh_147, hh_150, hh_152, \
                         hh_157, hh_159, hh_161, hh_189, hh_192, hh_194, hh_199, hh_201, \
                         hh_203, hh_336, hh_339, hh_341, hh_346, hh_348, hh_350, hh_378, \
                         hh_381, hh_383, hh_388, hh_390, hh_392, hh_420, hh_423, hh_425, \
                         hh_430, hh_432, hh_434 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_137 * hh_42[k]
                  + f_138 * hh_45[k]
                  - f_144 * hh_47[k]
                  + f_137 * hh_52[k]
                  - f_144 * hh_54[k]
                  + f_148 * hh_56[k]
                  + f_138 * hh_147[k]
                  + f_141 * hh_150[k]
                  - f_145 * hh_152[k]
                  + f_138 * hh_157[k]
                  - f_145 * hh_159[k]
                  + f_149 * hh_161[k]
                  - f_139 * hh_189[k]
                  - f_142 * hh_192[k]
                  + f_146 * hh_194[k]
                  - f_139 * hh_199[k]
                  + f_146 * hh_201[k]
                  - f_150 * hh_203[k]
                  + f_137 * hh_336[k]
                  + f_138 * hh_339[k]
                  - f_144 * hh_341[k]
                  + f_137 * hh_346[k]
                  - f_144 * hh_348[k]
                  + f_148 * hh_350[k]
                  - f_139 * hh_378[k]
                  - f_142 * hh_381[k]
                  + f_146 * hh_383[k]
                  - f_139 * hh_388[k]
                  + f_146 * hh_390[k]
                  - f_150 * hh_392[k]
                  + f_140 * hh_420[k]
                  + f_143 * hh_423[k]
                  - f_147 * hh_425[k]
                  + f_140 * hh_430[k]
                  - f_147 * hh_432[k]
                  + f_151 * hh_434[k];
        g_71[k] = g_61[k];
    }

#pragma omp simd aligned(hh_44, hh_51, hh_58, hh_60, hh_149, hh_156, hh_163, hh_165, hh_191, \
                         hh_198, hh_205, hh_207, hh_338, hh_345, hh_352, hh_354, hh_380, \
                         hh_387, hh_394, hh_396, hh_422, hh_429, hh_436, \
                         hh_438 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_159 * hh_44[k]
                  + f_127 * hh_51[k]
                  + f_159 * hh_58[k]
                  - f_127 * hh_60[k]
                  - f_127 * hh_149[k]
                  + f_128 * hh_156[k]
                  + f_127 * hh_163[k]
                  - f_128 * hh_165[k]
                  + f_160 * hh_191[k]
                  - f_129 * hh_198[k]
                  - f_160 * hh_205[k]
                  + f_129 * hh_207[k]
                  - f_159 * hh_338[k]
                  + f_127 * hh_345[k]
                  + f_159 * hh_352[k]
                  - f_127 * hh_354[k]
                  + f_160 * hh_380[k]
                  - f_129 * hh_387[k]
                  - f_160 * hh_394[k]
                  + f_129 * hh_396[k]
                  - f_161 * hh_422[k]
                  + f_130 * hh_429[k]
                  + f_161 * hh_436[k]
                  - f_130 * hh_438[k];
        g_82[k] = g_62[k];
    }

#pragma omp simd aligned(hh_42, hh_45, hh_47, hh_52, hh_54, hh_147, hh_150, hh_152, hh_157, \
                         hh_159, hh_189, hh_192, hh_194, hh_199, hh_201, hh_336, hh_339, \
                         hh_341, hh_346, hh_348, hh_378, hh_381, hh_383, hh_388, hh_390, \
                         hh_420, hh_423, hh_425, hh_430, hh_432 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_103 * hh_42[k]
                  + f_95 * hh_45[k]
                  + f_93 * hh_47[k]
                  + f_91 * hh_52[k]
                  - f_99 * hh_54[k]
                  - f_95 * hh_147[k]
                  + f_96 * hh_150[k]
                  + f_106 * hh_152[k]
                  + f_92 * hh_157[k]
                  - f_100 * hh_159[k]
                  + f_104 * hh_189[k]
                  - f_97 * hh_192[k]
                  - f_107 * hh_194[k]
                  - f_93 * hh_199[k]
                  + f_101 * hh_201[k]
                  - f_103 * hh_336[k]
                  + f_95 * hh_339[k]
                  + f_93 * hh_341[k]
                  + f_91 * hh_346[k]
                  - f_99 * hh_348[k]
                  + f_104 * hh_378[k]
                  - f_97 * hh_381[k]
                  - f_107 * hh_383[k]
                  - f_93 * hh_388[k]
                  + f_101 * hh_390[k]
                  - f_105 * hh_420[k]
                  + f_98 * hh_423[k]
                  + f_108 * hh_425[k]
                  + f_94 * hh_430[k]
                  - f_102 * hh_432[k];
        g_93[k] = g_63[k];
    }

#pragma omp simd aligned(hh_44, hh_49, hh_58, hh_149, hh_154, hh_163, hh_191, hh_196, hh_205, \
                         hh_338, hh_343, hh_352, hh_380, hh_385, hh_394, hh_422, hh_427, \
                         hh_436 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_162 * hh_44[k]
                  - f_163 * hh_49[k]
                  + f_162 * hh_58[k]
                  + f_164 * hh_149[k]
                  - f_165 * hh_154[k]
                  + f_164 * hh_163[k]
                  - f_166 * hh_191[k]
                  + f_167 * hh_196[k]
                  - f_166 * hh_205[k]
                  + f_162 * hh_338[k]
                  - f_163 * hh_343[k]
                  + f_162 * hh_352[k]
                  - f_166 * hh_380[k]
                  + f_167 * hh_385[k]
                  - f_166 * hh_394[k]
                  + f_168 * hh_422[k]
                  - f_169 * hh_427[k]
                  + f_168 * hh_436[k];
        g_104[k] = g_64[k];
    }

#pragma omp simd aligned(hh_42, hh_45, hh_52, hh_147, hh_150, hh_157, hh_189, hh_192, hh_199, \
                         hh_336, hh_339, hh_346, hh_378, hh_381, hh_388, hh_420, hh_423, \
                         hh_430 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_40 * hh_42[k]
                  - f_34 * hh_45[k]
                  + f_33 * hh_52[k]
                  + f_41 * hh_147[k]
                  - f_37 * hh_150[k]
                  + f_34 * hh_157[k]
                  - f_36 * hh_189[k]
                  + f_38 * hh_192[k]
                  - f_35 * hh_199[k]
                  + f_40 * hh_336[k]
                  - f_34 * hh_339[k]
                  + f_33 * hh_346[k]
                  - f_36 * hh_378[k]
                  + f_38 * hh_381[k]
                  - f_35 * hh_388[k]
                  + f_42 * hh_420[k]
                  - f_39 * hh_423[k]
                  + f_36 * hh_430[k];
        g_115[k] = g_65[k];
    }

#pragma omp simd aligned(hh_0, hh_3, hh_5, hh_10, hh_12, hh_14, hh_63, hh_66, hh_68, hh_73, \
                         hh_75, hh_77, hh_105, hh_108, hh_110, hh_115, hh_117, hh_119, hh_210, \
                         hh_213, hh_215, hh_220, hh_222, hh_224, hh_252, hh_255, hh_257, \
                         hh_262, hh_264, hh_266, hh_294, hh_297, hh_299, hh_304, hh_306, \
                         hh_308 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = 0.234375 * hh_0[k]
                  + 0.46875 * hh_3[k]
                  - 2.8125 * hh_5[k]
                  + 0.234375 * hh_10[k]
                  - 2.8125 * hh_12[k]
                  + 1.875 * hh_14[k]
                  + 0.46875 * hh_63[k]
                  + 0.9375 * hh_66[k]
                  - 5.625 * hh_68[k]
                  + 0.46875 * hh_73[k]
                  - 5.625 * hh_75[k]
                  + 3.75 * hh_77[k]
                  - 2.8125 * hh_105[k]
                  - 5.625 * hh_108[k]
                  + 33.75 * hh_110[k]
                  - 2.8125 * hh_115[k]
                  + 33.75 * hh_117[k]
                  - 22.5 * hh_119[k]
                  + 0.234375 * hh_210[k]
                  + 0.46875 * hh_213[k]
                  - 2.8125 * hh_215[k]
                  + 0.234375 * hh_220[k]
                  - 2.8125 * hh_222[k]
                  + 1.875 * hh_224[k]
                  - 2.8125 * hh_252[k]
                  - 5.625 * hh_255[k]
                  + 33.75 * hh_257[k]
                  - 2.8125 * hh_262[k]
                  + 33.75 * hh_264[k]
                  - 22.5 * hh_266[k]
                  + 1.875 * hh_294[k]
                  + 3.75 * hh_297[k]
                  - 22.5 * hh_299[k]
                  + 1.875 * hh_304[k]
                  - 22.5 * hh_306[k]
                  + 15.0 * hh_308[k];
    }

#pragma omp simd aligned(hh_2, hh_9, hh_16, hh_18, hh_65, hh_72, hh_79, hh_81, hh_107, hh_114, \
                         hh_121, hh_123, hh_212, hh_219, hh_226, hh_228, hh_254, hh_261, \
                         hh_268, hh_270, hh_296, hh_303, hh_310, \
                         hh_312 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = -f_152 * hh_2[k]
                  + f_120 * hh_9[k]
                  + f_152 * hh_16[k]
                  - f_120 * hh_18[k]
                  - f_120 * hh_65[k]
                  + f_121 * hh_72[k]
                  + f_120 * hh_79[k]
                  - f_121 * hh_81[k]
                  + f_153 * hh_107[k]
                  - f_122 * hh_114[k]
                  - f_153 * hh_121[k]
                  + f_122 * hh_123[k]
                  - f_152 * hh_212[k]
                  + f_120 * hh_219[k]
                  + f_152 * hh_226[k]
                  - f_120 * hh_228[k]
                  + f_153 * hh_254[k]
                  - f_122 * hh_261[k]
                  - f_153 * hh_268[k]
                  + f_122 * hh_270[k]
                  - f_124 * hh_296[k]
                  + f_123 * hh_303[k]
                  + f_124 * hh_310[k]
                  - f_123 * hh_312[k];
        g_83[k] = g_73[k];
    }

#pragma omp simd aligned(hh_0, hh_3, hh_5, hh_10, hh_12, hh_63, hh_66, hh_68, hh_73, hh_75, \
                         hh_105, hh_108, hh_110, hh_115, hh_117, hh_210, hh_213, hh_215, \
                         hh_220, hh_222, hh_252, hh_255, hh_257, hh_262, hh_264, hh_294, \
                         hh_297, hh_299, hh_304, hh_306 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_86 * hh_0[k]
                  + f_80 * hh_3[k]
                  + f_88 * hh_5[k]
                  + f_76 * hh_10[k]
                  - f_79 * hh_12[k]
                  - f_80 * hh_63[k]
                  + f_81 * hh_66[k]
                  + f_82 * hh_68[k]
                  + f_77 * hh_73[k]
                  - f_83 * hh_75[k]
                  + f_87 * hh_105[k]
                  - f_79 * hh_108[k]
                  - f_89 * hh_110[k]
                  - f_78 * hh_115[k]
                  + f_84 * hh_117[k]
                  - f_86 * hh_210[k]
                  + f_80 * hh_213[k]
                  + f_88 * hh_215[k]
                  + f_76 * hh_220[k]
                  - f_79 * hh_222[k]
                  + f_87 * hh_252[k]
                  - f_79 * hh_255[k]
                  - f_89 * hh_257[k]
                  - f_78 * hh_262[k]
                  + f_84 * hh_264[k]
                  - f_88 * hh_294[k]
                  + f_82 * hh_297[k]
                  + f_90 * hh_299[k]
                  + f_79 * hh_304[k]
                  - f_85 * hh_306[k];
        g_94[k] = g_74[k];
    }

#pragma omp simd aligned(hh_2, hh_7, hh_16, hh_65, hh_70, hh_79, hh_107, hh_112, hh_121, \
                         hh_212, hh_217, hh_226, hh_254, hh_259, hh_268, hh_296, hh_301, \
                         hh_310 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = f_154 * hh_2[k]
                  - f_155 * hh_7[k]
                  + f_154 * hh_16[k]
                  + f_156 * hh_65[k]
                  - f_157 * hh_70[k]
                  + f_156 * hh_79[k]
                  - f_157 * hh_107[k]
                  + f_158 * hh_112[k]
                  - f_157 * hh_121[k]
                  + f_154 * hh_212[k]
                  - f_155 * hh_217[k]
                  + f_154 * hh_226[k]
                  - f_157 * hh_254[k]
                  + f_158 * hh_259[k]
                  - f_157 * hh_268[k]
                  + f_59 * hh_296[k]
                  - f_60 * hh_301[k]
                  + f_59 * hh_310[k];
        g_105[k] = g_75[k];
    }

#pragma omp simd aligned(hh_0, hh_3, hh_10, hh_63, hh_66, hh_73, hh_105, hh_108, hh_115, \
                         hh_210, hh_213, hh_220, hh_252, hh_255, hh_262, hh_294, hh_297, \
                         hh_304 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = f_29 * hh_0[k]
                  - f_23 * hh_3[k]
                  + f_22 * hh_10[k]
                  + f_30 * hh_63[k]
                  - f_26 * hh_66[k]
                  + f_23 * hh_73[k]
                  - f_31 * hh_105[k]
                  + f_27 * hh_108[k]
                  - f_24 * hh_115[k]
                  + f_29 * hh_210[k]
                  - f_23 * hh_213[k]
                  + f_22 * hh_220[k]
                  - f_31 * hh_252[k]
                  + f_27 * hh_255[k]
                  - f_24 * hh_262[k]
                  + f_32 * hh_294[k]
                  - f_28 * hh_297[k]
                  + f_25 * hh_304[k];
        g_116[k] = g_76[k];
    }

#pragma omp simd aligned(hh_44, hh_51, hh_58, hh_60, hh_191, hh_198, hh_205, hh_207, hh_338, \
                         hh_345, hh_352, hh_354, hh_380, hh_387, hh_394, \
                         hh_396 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = 6.5625 * hh_44[k]
                  - 13.125 * hh_51[k]
                  - 6.5625 * hh_58[k]
                  + 13.125 * hh_60[k]
                  - 13.125 * hh_191[k]
                  + 26.25 * hh_198[k]
                  + 13.125 * hh_205[k]
                  - 26.25 * hh_207[k]
                  - 6.5625 * hh_338[k]
                  + 13.125 * hh_345[k]
                  + 6.5625 * hh_352[k]
                  - 13.125 * hh_354[k]
                  + 13.125 * hh_380[k]
                  - 26.25 * hh_387[k]
                  - 13.125 * hh_394[k]
                  + 26.25 * hh_396[k];
    }

#pragma omp simd aligned(hh_42, hh_45, hh_47, hh_52, hh_54, hh_189, hh_192, hh_194, hh_199, \
                         hh_201, hh_336, hh_339, hh_341, hh_346, hh_348, hh_378, hh_381, \
                         hh_383, hh_388, hh_390 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = f_111 * hh_42[k]
                  - f_73 * hh_45[k]
                  - f_70 * hh_47[k]
                  - f_109 * hh_52[k]
                  + f_110 * hh_54[k]
                  - f_73 * hh_189[k]
                  + f_69 * hh_192[k]
                  + f_74 * hh_194[k]
                  + f_67 * hh_199[k]
                  - f_71 * hh_201[k]
                  - f_111 * hh_336[k]
                  + f_73 * hh_339[k]
                  + f_70 * hh_341[k]
                  + f_109 * hh_346[k]
                  - f_110 * hh_348[k]
                  + f_73 * hh_378[k]
                  - f_69 * hh_381[k]
                  - f_74 * hh_383[k]
                  - f_67 * hh_388[k]
                  + f_71 * hh_390[k];
        g_95[k] = g_85[k];
    }

#pragma omp simd aligned(hh_44, hh_49, hh_58, hh_191, hh_196, hh_205, hh_338, hh_343, hh_352, \
                         hh_380, hh_385, hh_394 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_170 * hh_44[k]
                  + f_171 * hh_49[k]
                  - f_170 * hh_58[k]
                  + f_134 * hh_191[k]
                  - f_135 * hh_196[k]
                  + f_134 * hh_205[k]
                  + f_170 * hh_338[k]
                  - f_171 * hh_343[k]
                  + f_170 * hh_352[k]
                  - f_134 * hh_380[k]
                  + f_135 * hh_385[k]
                  - f_134 * hh_394[k];
        g_106[k] = g_86[k];
    }

#pragma omp simd aligned(hh_42, hh_45, hh_52, hh_189, hh_192, hh_199, hh_336, hh_339, hh_346, \
                         hh_378, hh_381, hh_388 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = -f_44 * hh_42[k]
                  + f_17 * hh_45[k]
                  - f_43 * hh_52[k]
                  + f_20 * hh_189[k]
                  - f_18 * hh_192[k]
                  + f_17 * hh_199[k]
                  + f_44 * hh_336[k]
                  - f_17 * hh_339[k]
                  + f_43 * hh_346[k]
                  - f_20 * hh_378[k]
                  + f_18 * hh_381[k]
                  - f_17 * hh_388[k];
        g_117[k] = g_87[k];
    }

#pragma omp simd aligned(hh_0, hh_3, hh_5, hh_10, hh_12, hh_63, hh_66, hh_68, hh_73, hh_75, \
                         hh_105, hh_108, hh_110, hh_115, hh_117, hh_210, hh_213, hh_215, \
                         hh_220, hh_222, hh_252, hh_255, hh_257, hh_262, \
                         hh_264 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = 0.2734375 * hh_0[k]
                  - 0.546875 * hh_3[k]
                  - 2.1875 * hh_5[k]
                  - 0.8203125 * hh_10[k]
                  + 6.5625 * hh_12[k]
                  - 0.546875 * hh_63[k]
                  + 1.09375 * hh_66[k]
                  + 4.375 * hh_68[k]
                  + 1.640625 * hh_73[k]
                  - 13.125 * hh_75[k]
                  - 2.1875 * hh_105[k]
                  + 4.375 * hh_108[k]
                  + 17.5 * hh_110[k]
                  + 6.5625 * hh_115[k]
                  - 52.5 * hh_117[k]
                  - 0.8203125 * hh_210[k]
                  + 1.640625 * hh_213[k]
                  + 6.5625 * hh_215[k]
                  + 2.4609375 * hh_220[k]
                  - 19.6875 * hh_222[k]
                  + 6.5625 * hh_252[k]
                  - 13.125 * hh_255[k]
                  - 52.5 * hh_257[k]
                  - 19.6875 * hh_262[k]
                  + 157.5 * hh_264[k];
    }

#pragma omp simd aligned(hh_2, hh_7, hh_16, hh_65, hh_70, hh_79, hh_107, hh_112, hh_121, \
                         hh_212, hh_217, hh_226, hh_254, hh_259, \
                         hh_268 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = -f_117 * hh_2[k]
                  + f_118 * hh_7[k]
                  - f_117 * hh_16[k]
                  + f_114 * hh_65[k]
                  - f_51 * hh_70[k]
                  + f_114 * hh_79[k]
                  + f_52 * hh_107[k]
                  - f_119 * hh_112[k]
                  + f_52 * hh_121[k]
                  + f_112 * hh_212[k]
                  - f_113 * hh_217[k]
                  + f_112 * hh_226[k]
                  - f_115 * hh_254[k]
                  + f_116 * hh_259[k]
                  - f_115 * hh_268[k];
        g_107[k] = g_97[k];
    }

#pragma omp simd aligned(hh_0, hh_3, hh_10, hh_63, hh_66, hh_73, hh_105, hh_108, hh_115, \
                         hh_210, hh_213, hh_220, hh_252, hh_255, \
                         hh_262 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = -f_15 * hh_0[k]
                  + f_4 * hh_3[k]
                  - f_6 * hh_10[k]
                  + f_13 * hh_63[k]
                  - f_9 * hh_66[k]
                  + f_4 * hh_73[k]
                  + f_16 * hh_105[k]
                  - f_11 * hh_108[k]
                  + f_7 * hh_115[k]
                  + f_12 * hh_210[k]
                  - f_8 * hh_213[k]
                  + f_3 * hh_220[k]
                  - f_14 * hh_252[k]
                  + f_10 * hh_255[k]
                  - f_5 * hh_262[k];
        g_118[k] = g_98[k];
    }

#pragma omp simd aligned(hh_44, hh_49, hh_58, hh_149, hh_154, hh_163, hh_338, hh_343, \
                         hh_352 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = 4.921875 * hh_44[k]
                   - 29.53125 * hh_49[k]
                   + 4.921875 * hh_58[k]
                   - 29.53125 * hh_149[k]
                   + 177.1875 * hh_154[k]
                   - 29.53125 * hh_163[k]
                   + 4.921875 * hh_338[k]
                   - 29.53125 * hh_343[k]
                   + 4.921875 * hh_352[k];
    }

#pragma omp simd aligned(hh_42, hh_45, hh_52, hh_147, hh_150, hh_157, hh_336, hh_339, \
                         hh_346 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = f_49 * hh_42[k]
                   - f_47 * hh_45[k]
                   + f_45 * hh_52[k]
                   - f_50 * hh_147[k]
                   + f_48 * hh_150[k]
                   - f_46 * hh_157[k]
                   + f_49 * hh_336[k]
                   - f_47 * hh_339[k]
                   + f_45 * hh_346[k];
        g_119[k] = g_109[k];
    }

#pragma omp simd aligned(hh_0, hh_3, hh_10, hh_63, hh_66, hh_73, hh_210, hh_213, \
                         hh_220 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_120[k] = 0.4921875 * hh_0[k]
                   - 4.921875 * hh_3[k]
                   + 2.4609375 * hh_10[k]
                   - 4.921875 * hh_63[k]
                   + 49.21875 * hh_66[k]
                   - 24.609375 * hh_73[k]
                   + 2.4609375 * hh_210[k]
                   - 24.609375 * hh_213[k]
                   + 12.3046875 * hh_220[k];
    }
}

}  // namespace simdtrf
