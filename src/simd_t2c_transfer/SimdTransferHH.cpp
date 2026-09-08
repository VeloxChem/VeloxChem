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


#include "SimdTransferHH.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_hh_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t gh, const size_t gi,
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_1 = buffer.data(gh + 1);
    const auto *gh_2 = buffer.data(gh + 2);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_4 = buffer.data(gh + 4);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_6 = buffer.data(gh + 6);
    const auto *gh_7 = buffer.data(gh + 7);
    const auto *gh_8 = buffer.data(gh + 8);
    const auto *gh_9 = buffer.data(gh + 9);
    const auto *gh_10 = buffer.data(gh + 10);
    const auto *gh_11 = buffer.data(gh + 11);
    const auto *gh_12 = buffer.data(gh + 12);
    const auto *gh_13 = buffer.data(gh + 13);
    const auto *gh_14 = buffer.data(gh + 14);
    const auto *gh_15 = buffer.data(gh + 15);
    const auto *gh_16 = buffer.data(gh + 16);
    const auto *gh_17 = buffer.data(gh + 17);
    const auto *gh_18 = buffer.data(gh + 18);
    const auto *gh_19 = buffer.data(gh + 19);
    const auto *gh_20 = buffer.data(gh + 20);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_23 = buffer.data(gh + 23);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_25 = buffer.data(gh + 25);
    const auto *gh_26 = buffer.data(gh + 26);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_28 = buffer.data(gh + 28);
    const auto *gh_29 = buffer.data(gh + 29);
    const auto *gh_30 = buffer.data(gh + 30);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_32 = buffer.data(gh + 32);
    const auto *gh_33 = buffer.data(gh + 33);
    const auto *gh_34 = buffer.data(gh + 34);
    const auto *gh_35 = buffer.data(gh + 35);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_37 = buffer.data(gh + 37);
    const auto *gh_38 = buffer.data(gh + 38);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_40 = buffer.data(gh + 40);
    const auto *gh_41 = buffer.data(gh + 41);
    const auto *gh_42 = buffer.data(gh + 42);
    const auto *gh_43 = buffer.data(gh + 43);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_45 = buffer.data(gh + 45);
    const auto *gh_46 = buffer.data(gh + 46);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_48 = buffer.data(gh + 48);
    const auto *gh_49 = buffer.data(gh + 49);
    const auto *gh_50 = buffer.data(gh + 50);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_52 = buffer.data(gh + 52);
    const auto *gh_53 = buffer.data(gh + 53);
    const auto *gh_54 = buffer.data(gh + 54);
    const auto *gh_55 = buffer.data(gh + 55);
    const auto *gh_56 = buffer.data(gh + 56);
    const auto *gh_57 = buffer.data(gh + 57);
    const auto *gh_58 = buffer.data(gh + 58);
    const auto *gh_59 = buffer.data(gh + 59);
    const auto *gh_60 = buffer.data(gh + 60);
    const auto *gh_61 = buffer.data(gh + 61);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_64 = buffer.data(gh + 64);
    const auto *gh_65 = buffer.data(gh + 65);
    const auto *gh_66 = buffer.data(gh + 66);
    const auto *gh_67 = buffer.data(gh + 67);
    const auto *gh_68 = buffer.data(gh + 68);
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_70 = buffer.data(gh + 70);
    const auto *gh_71 = buffer.data(gh + 71);
    const auto *gh_72 = buffer.data(gh + 72);
    const auto *gh_73 = buffer.data(gh + 73);
    const auto *gh_74 = buffer.data(gh + 74);
    const auto *gh_75 = buffer.data(gh + 75);
    const auto *gh_76 = buffer.data(gh + 76);
    const auto *gh_77 = buffer.data(gh + 77);
    const auto *gh_78 = buffer.data(gh + 78);
    const auto *gh_79 = buffer.data(gh + 79);
    const auto *gh_80 = buffer.data(gh + 80);
    const auto *gh_81 = buffer.data(gh + 81);
    const auto *gh_82 = buffer.data(gh + 82);
    const auto *gh_83 = buffer.data(gh + 83);
    const auto *gh_84 = buffer.data(gh + 84);
    const auto *gh_85 = buffer.data(gh + 85);
    const auto *gh_86 = buffer.data(gh + 86);
    const auto *gh_87 = buffer.data(gh + 87);
    const auto *gh_88 = buffer.data(gh + 88);
    const auto *gh_89 = buffer.data(gh + 89);
    const auto *gh_90 = buffer.data(gh + 90);
    const auto *gh_91 = buffer.data(gh + 91);
    const auto *gh_92 = buffer.data(gh + 92);
    const auto *gh_93 = buffer.data(gh + 93);
    const auto *gh_94 = buffer.data(gh + 94);
    const auto *gh_95 = buffer.data(gh + 95);
    const auto *gh_96 = buffer.data(gh + 96);
    const auto *gh_97 = buffer.data(gh + 97);
    const auto *gh_98 = buffer.data(gh + 98);
    const auto *gh_99 = buffer.data(gh + 99);
    const auto *gh_100 = buffer.data(gh + 100);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_103 = buffer.data(gh + 103);
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_106 = buffer.data(gh + 106);
    const auto *gh_107 = buffer.data(gh + 107);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_109 = buffer.data(gh + 109);
    const auto *gh_110 = buffer.data(gh + 110);
    const auto *gh_111 = buffer.data(gh + 111);
    const auto *gh_112 = buffer.data(gh + 112);
    const auto *gh_113 = buffer.data(gh + 113);
    const auto *gh_114 = buffer.data(gh + 114);
    const auto *gh_115 = buffer.data(gh + 115);
    const auto *gh_116 = buffer.data(gh + 116);
    const auto *gh_117 = buffer.data(gh + 117);
    const auto *gh_118 = buffer.data(gh + 118);
    const auto *gh_119 = buffer.data(gh + 119);
    const auto *gh_120 = buffer.data(gh + 120);
    const auto *gh_121 = buffer.data(gh + 121);
    const auto *gh_122 = buffer.data(gh + 122);
    const auto *gh_123 = buffer.data(gh + 123);
    const auto *gh_124 = buffer.data(gh + 124);
    const auto *gh_125 = buffer.data(gh + 125);
    const auto *gh_126 = buffer.data(gh + 126);
    const auto *gh_127 = buffer.data(gh + 127);
    const auto *gh_128 = buffer.data(gh + 128);
    const auto *gh_129 = buffer.data(gh + 129);
    const auto *gh_130 = buffer.data(gh + 130);
    const auto *gh_131 = buffer.data(gh + 131);
    const auto *gh_132 = buffer.data(gh + 132);
    const auto *gh_133 = buffer.data(gh + 133);
    const auto *gh_134 = buffer.data(gh + 134);
    const auto *gh_135 = buffer.data(gh + 135);
    const auto *gh_136 = buffer.data(gh + 136);
    const auto *gh_137 = buffer.data(gh + 137);
    const auto *gh_138 = buffer.data(gh + 138);
    const auto *gh_139 = buffer.data(gh + 139);
    const auto *gh_140 = buffer.data(gh + 140);
    const auto *gh_141 = buffer.data(gh + 141);
    const auto *gh_142 = buffer.data(gh + 142);
    const auto *gh_143 = buffer.data(gh + 143);
    const auto *gh_144 = buffer.data(gh + 144);
    const auto *gh_145 = buffer.data(gh + 145);
    const auto *gh_146 = buffer.data(gh + 146);
    const auto *gh_147 = buffer.data(gh + 147);
    const auto *gh_148 = buffer.data(gh + 148);
    const auto *gh_149 = buffer.data(gh + 149);
    const auto *gh_150 = buffer.data(gh + 150);
    const auto *gh_151 = buffer.data(gh + 151);
    const auto *gh_152 = buffer.data(gh + 152);
    const auto *gh_153 = buffer.data(gh + 153);
    const auto *gh_154 = buffer.data(gh + 154);
    const auto *gh_155 = buffer.data(gh + 155);
    const auto *gh_156 = buffer.data(gh + 156);
    const auto *gh_157 = buffer.data(gh + 157);
    const auto *gh_158 = buffer.data(gh + 158);
    const auto *gh_159 = buffer.data(gh + 159);
    const auto *gh_160 = buffer.data(gh + 160);
    const auto *gh_161 = buffer.data(gh + 161);
    const auto *gh_162 = buffer.data(gh + 162);
    const auto *gh_163 = buffer.data(gh + 163);
    const auto *gh_164 = buffer.data(gh + 164);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_166 = buffer.data(gh + 166);
    const auto *gh_167 = buffer.data(gh + 167);
    const auto *gh_168 = buffer.data(gh + 168);
    const auto *gh_169 = buffer.data(gh + 169);
    const auto *gh_170 = buffer.data(gh + 170);
    const auto *gh_171 = buffer.data(gh + 171);
    const auto *gh_172 = buffer.data(gh + 172);
    const auto *gh_173 = buffer.data(gh + 173);
    const auto *gh_174 = buffer.data(gh + 174);
    const auto *gh_175 = buffer.data(gh + 175);
    const auto *gh_176 = buffer.data(gh + 176);
    const auto *gh_177 = buffer.data(gh + 177);
    const auto *gh_178 = buffer.data(gh + 178);
    const auto *gh_179 = buffer.data(gh + 179);
    const auto *gh_180 = buffer.data(gh + 180);
    const auto *gh_181 = buffer.data(gh + 181);
    const auto *gh_182 = buffer.data(gh + 182);
    const auto *gh_183 = buffer.data(gh + 183);
    const auto *gh_184 = buffer.data(gh + 184);
    const auto *gh_185 = buffer.data(gh + 185);
    const auto *gh_186 = buffer.data(gh + 186);
    const auto *gh_187 = buffer.data(gh + 187);
    const auto *gh_188 = buffer.data(gh + 188);
    const auto *gh_189 = buffer.data(gh + 189);
    const auto *gh_190 = buffer.data(gh + 190);
    const auto *gh_191 = buffer.data(gh + 191);
    const auto *gh_192 = buffer.data(gh + 192);
    const auto *gh_193 = buffer.data(gh + 193);
    const auto *gh_194 = buffer.data(gh + 194);
    const auto *gh_195 = buffer.data(gh + 195);
    const auto *gh_196 = buffer.data(gh + 196);
    const auto *gh_197 = buffer.data(gh + 197);
    const auto *gh_198 = buffer.data(gh + 198);
    const auto *gh_199 = buffer.data(gh + 199);
    const auto *gh_200 = buffer.data(gh + 200);
    const auto *gh_201 = buffer.data(gh + 201);
    const auto *gh_202 = buffer.data(gh + 202);
    const auto *gh_203 = buffer.data(gh + 203);
    const auto *gh_204 = buffer.data(gh + 204);
    const auto *gh_205 = buffer.data(gh + 205);
    const auto *gh_206 = buffer.data(gh + 206);
    const auto *gh_207 = buffer.data(gh + 207);
    const auto *gh_208 = buffer.data(gh + 208);
    const auto *gh_209 = buffer.data(gh + 209);
    const auto *gh_210 = buffer.data(gh + 210);
    const auto *gh_211 = buffer.data(gh + 211);
    const auto *gh_212 = buffer.data(gh + 212);
    const auto *gh_213 = buffer.data(gh + 213);
    const auto *gh_214 = buffer.data(gh + 214);
    const auto *gh_215 = buffer.data(gh + 215);
    const auto *gh_216 = buffer.data(gh + 216);
    const auto *gh_217 = buffer.data(gh + 217);
    const auto *gh_218 = buffer.data(gh + 218);
    const auto *gh_219 = buffer.data(gh + 219);
    const auto *gh_220 = buffer.data(gh + 220);
    const auto *gh_221 = buffer.data(gh + 221);
    const auto *gh_222 = buffer.data(gh + 222);
    const auto *gh_223 = buffer.data(gh + 223);
    const auto *gh_224 = buffer.data(gh + 224);
    const auto *gh_225 = buffer.data(gh + 225);
    const auto *gh_226 = buffer.data(gh + 226);
    const auto *gh_227 = buffer.data(gh + 227);
    const auto *gh_228 = buffer.data(gh + 228);
    const auto *gh_229 = buffer.data(gh + 229);
    const auto *gh_230 = buffer.data(gh + 230);
    const auto *gh_231 = buffer.data(gh + 231);
    const auto *gh_232 = buffer.data(gh + 232);
    const auto *gh_233 = buffer.data(gh + 233);
    const auto *gh_234 = buffer.data(gh + 234);
    const auto *gh_235 = buffer.data(gh + 235);
    const auto *gh_236 = buffer.data(gh + 236);
    const auto *gh_237 = buffer.data(gh + 237);
    const auto *gh_238 = buffer.data(gh + 238);
    const auto *gh_239 = buffer.data(gh + 239);
    const auto *gh_240 = buffer.data(gh + 240);
    const auto *gh_241 = buffer.data(gh + 241);
    const auto *gh_242 = buffer.data(gh + 242);
    const auto *gh_243 = buffer.data(gh + 243);
    const auto *gh_244 = buffer.data(gh + 244);
    const auto *gh_245 = buffer.data(gh + 245);
    const auto *gh_246 = buffer.data(gh + 246);
    const auto *gh_247 = buffer.data(gh + 247);
    const auto *gh_248 = buffer.data(gh + 248);
    const auto *gh_249 = buffer.data(gh + 249);
    const auto *gh_250 = buffer.data(gh + 250);
    const auto *gh_251 = buffer.data(gh + 251);
    const auto *gh_252 = buffer.data(gh + 252);
    const auto *gh_253 = buffer.data(gh + 253);
    const auto *gh_254 = buffer.data(gh + 254);
    const auto *gh_255 = buffer.data(gh + 255);
    const auto *gh_256 = buffer.data(gh + 256);
    const auto *gh_257 = buffer.data(gh + 257);
    const auto *gh_258 = buffer.data(gh + 258);
    const auto *gh_259 = buffer.data(gh + 259);
    const auto *gh_260 = buffer.data(gh + 260);
    const auto *gh_261 = buffer.data(gh + 261);
    const auto *gh_262 = buffer.data(gh + 262);
    const auto *gh_263 = buffer.data(gh + 263);
    const auto *gh_264 = buffer.data(gh + 264);
    const auto *gh_265 = buffer.data(gh + 265);
    const auto *gh_266 = buffer.data(gh + 266);
    const auto *gh_267 = buffer.data(gh + 267);
    const auto *gh_268 = buffer.data(gh + 268);
    const auto *gh_269 = buffer.data(gh + 269);
    const auto *gh_270 = buffer.data(gh + 270);
    const auto *gh_271 = buffer.data(gh + 271);
    const auto *gh_272 = buffer.data(gh + 272);
    const auto *gh_273 = buffer.data(gh + 273);
    const auto *gh_274 = buffer.data(gh + 274);
    const auto *gh_275 = buffer.data(gh + 275);
    const auto *gh_276 = buffer.data(gh + 276);
    const auto *gh_277 = buffer.data(gh + 277);
    const auto *gh_278 = buffer.data(gh + 278);
    const auto *gh_279 = buffer.data(gh + 279);
    const auto *gh_280 = buffer.data(gh + 280);
    const auto *gh_281 = buffer.data(gh + 281);
    const auto *gh_282 = buffer.data(gh + 282);
    const auto *gh_283 = buffer.data(gh + 283);
    const auto *gh_284 = buffer.data(gh + 284);
    const auto *gh_285 = buffer.data(gh + 285);
    const auto *gh_286 = buffer.data(gh + 286);
    const auto *gh_287 = buffer.data(gh + 287);
    const auto *gh_288 = buffer.data(gh + 288);
    const auto *gh_289 = buffer.data(gh + 289);
    const auto *gh_290 = buffer.data(gh + 290);
    const auto *gh_291 = buffer.data(gh + 291);
    const auto *gh_292 = buffer.data(gh + 292);
    const auto *gh_293 = buffer.data(gh + 293);
    const auto *gh_294 = buffer.data(gh + 294);
    const auto *gh_295 = buffer.data(gh + 295);
    const auto *gh_296 = buffer.data(gh + 296);
    const auto *gh_297 = buffer.data(gh + 297);
    const auto *gh_298 = buffer.data(gh + 298);
    const auto *gh_299 = buffer.data(gh + 299);
    const auto *gh_300 = buffer.data(gh + 300);
    const auto *gh_301 = buffer.data(gh + 301);
    const auto *gh_302 = buffer.data(gh + 302);
    const auto *gh_303 = buffer.data(gh + 303);
    const auto *gh_304 = buffer.data(gh + 304);
    const auto *gh_305 = buffer.data(gh + 305);
    const auto *gh_306 = buffer.data(gh + 306);
    const auto *gh_307 = buffer.data(gh + 307);
    const auto *gh_308 = buffer.data(gh + 308);
    const auto *gh_309 = buffer.data(gh + 309);
    const auto *gh_310 = buffer.data(gh + 310);
    const auto *gh_311 = buffer.data(gh + 311);
    const auto *gh_312 = buffer.data(gh + 312);
    const auto *gh_313 = buffer.data(gh + 313);
    const auto *gh_314 = buffer.data(gh + 314);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_1 = buffer.data(gi + 1);
    const auto *gi_2 = buffer.data(gi + 2);
    const auto *gi_3 = buffer.data(gi + 3);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_8 = buffer.data(gi + 8);
    const auto *gi_9 = buffer.data(gi + 9);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_15 = buffer.data(gi + 15);
    const auto *gi_16 = buffer.data(gi + 16);
    const auto *gi_17 = buffer.data(gi + 17);
    const auto *gi_18 = buffer.data(gi + 18);
    const auto *gi_19 = buffer.data(gi + 19);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_30 = buffer.data(gi + 30);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_35 = buffer.data(gi + 35);
    const auto *gi_36 = buffer.data(gi + 36);
    const auto *gi_37 = buffer.data(gi + 37);
    const auto *gi_38 = buffer.data(gi + 38);
    const auto *gi_39 = buffer.data(gi + 39);
    const auto *gi_40 = buffer.data(gi + 40);
    const auto *gi_41 = buffer.data(gi + 41);
    const auto *gi_42 = buffer.data(gi + 42);
    const auto *gi_43 = buffer.data(gi + 43);
    const auto *gi_44 = buffer.data(gi + 44);
    const auto *gi_45 = buffer.data(gi + 45);
    const auto *gi_46 = buffer.data(gi + 46);
    const auto *gi_47 = buffer.data(gi + 47);
    const auto *gi_48 = buffer.data(gi + 48);
    const auto *gi_56 = buffer.data(gi + 56);
    const auto *gi_57 = buffer.data(gi + 57);
    const auto *gi_58 = buffer.data(gi + 58);
    const auto *gi_59 = buffer.data(gi + 59);
    const auto *gi_60 = buffer.data(gi + 60);
    const auto *gi_61 = buffer.data(gi + 61);
    const auto *gi_62 = buffer.data(gi + 62);
    const auto *gi_63 = buffer.data(gi + 63);
    const auto *gi_64 = buffer.data(gi + 64);
    const auto *gi_65 = buffer.data(gi + 65);
    const auto *gi_66 = buffer.data(gi + 66);
    const auto *gi_67 = buffer.data(gi + 67);
    const auto *gi_68 = buffer.data(gi + 68);
    const auto *gi_69 = buffer.data(gi + 69);
    const auto *gi_70 = buffer.data(gi + 70);
    const auto *gi_71 = buffer.data(gi + 71);
    const auto *gi_72 = buffer.data(gi + 72);
    const auto *gi_73 = buffer.data(gi + 73);
    const auto *gi_74 = buffer.data(gi + 74);
    const auto *gi_75 = buffer.data(gi + 75);
    const auto *gi_76 = buffer.data(gi + 76);
    const auto *gi_84 = buffer.data(gi + 84);
    const auto *gi_85 = buffer.data(gi + 85);
    const auto *gi_86 = buffer.data(gi + 86);
    const auto *gi_87 = buffer.data(gi + 87);
    const auto *gi_88 = buffer.data(gi + 88);
    const auto *gi_89 = buffer.data(gi + 89);
    const auto *gi_90 = buffer.data(gi + 90);
    const auto *gi_91 = buffer.data(gi + 91);
    const auto *gi_92 = buffer.data(gi + 92);
    const auto *gi_93 = buffer.data(gi + 93);
    const auto *gi_94 = buffer.data(gi + 94);
    const auto *gi_95 = buffer.data(gi + 95);
    const auto *gi_96 = buffer.data(gi + 96);
    const auto *gi_97 = buffer.data(gi + 97);
    const auto *gi_98 = buffer.data(gi + 98);
    const auto *gi_99 = buffer.data(gi + 99);
    const auto *gi_100 = buffer.data(gi + 100);
    const auto *gi_101 = buffer.data(gi + 101);
    const auto *gi_102 = buffer.data(gi + 102);
    const auto *gi_103 = buffer.data(gi + 103);
    const auto *gi_104 = buffer.data(gi + 104);
    const auto *gi_112 = buffer.data(gi + 112);
    const auto *gi_113 = buffer.data(gi + 113);
    const auto *gi_114 = buffer.data(gi + 114);
    const auto *gi_115 = buffer.data(gi + 115);
    const auto *gi_116 = buffer.data(gi + 116);
    const auto *gi_117 = buffer.data(gi + 117);
    const auto *gi_118 = buffer.data(gi + 118);
    const auto *gi_119 = buffer.data(gi + 119);
    const auto *gi_120 = buffer.data(gi + 120);
    const auto *gi_121 = buffer.data(gi + 121);
    const auto *gi_122 = buffer.data(gi + 122);
    const auto *gi_123 = buffer.data(gi + 123);
    const auto *gi_124 = buffer.data(gi + 124);
    const auto *gi_125 = buffer.data(gi + 125);
    const auto *gi_126 = buffer.data(gi + 126);
    const auto *gi_127 = buffer.data(gi + 127);
    const auto *gi_128 = buffer.data(gi + 128);
    const auto *gi_129 = buffer.data(gi + 129);
    const auto *gi_130 = buffer.data(gi + 130);
    const auto *gi_131 = buffer.data(gi + 131);
    const auto *gi_132 = buffer.data(gi + 132);
    const auto *gi_140 = buffer.data(gi + 140);
    const auto *gi_141 = buffer.data(gi + 141);
    const auto *gi_142 = buffer.data(gi + 142);
    const auto *gi_143 = buffer.data(gi + 143);
    const auto *gi_144 = buffer.data(gi + 144);
    const auto *gi_145 = buffer.data(gi + 145);
    const auto *gi_146 = buffer.data(gi + 146);
    const auto *gi_147 = buffer.data(gi + 147);
    const auto *gi_148 = buffer.data(gi + 148);
    const auto *gi_149 = buffer.data(gi + 149);
    const auto *gi_150 = buffer.data(gi + 150);
    const auto *gi_151 = buffer.data(gi + 151);
    const auto *gi_152 = buffer.data(gi + 152);
    const auto *gi_153 = buffer.data(gi + 153);
    const auto *gi_154 = buffer.data(gi + 154);
    const auto *gi_155 = buffer.data(gi + 155);
    const auto *gi_156 = buffer.data(gi + 156);
    const auto *gi_157 = buffer.data(gi + 157);
    const auto *gi_158 = buffer.data(gi + 158);
    const auto *gi_159 = buffer.data(gi + 159);
    const auto *gi_160 = buffer.data(gi + 160);
    const auto *gi_168 = buffer.data(gi + 168);
    const auto *gi_169 = buffer.data(gi + 169);
    const auto *gi_170 = buffer.data(gi + 170);
    const auto *gi_171 = buffer.data(gi + 171);
    const auto *gi_172 = buffer.data(gi + 172);
    const auto *gi_173 = buffer.data(gi + 173);
    const auto *gi_174 = buffer.data(gi + 174);
    const auto *gi_175 = buffer.data(gi + 175);
    const auto *gi_176 = buffer.data(gi + 176);
    const auto *gi_177 = buffer.data(gi + 177);
    const auto *gi_178 = buffer.data(gi + 178);
    const auto *gi_179 = buffer.data(gi + 179);
    const auto *gi_180 = buffer.data(gi + 180);
    const auto *gi_181 = buffer.data(gi + 181);
    const auto *gi_182 = buffer.data(gi + 182);
    const auto *gi_183 = buffer.data(gi + 183);
    const auto *gi_184 = buffer.data(gi + 184);
    const auto *gi_185 = buffer.data(gi + 185);
    const auto *gi_186 = buffer.data(gi + 186);
    const auto *gi_187 = buffer.data(gi + 187);
    const auto *gi_188 = buffer.data(gi + 188);
    const auto *gi_196 = buffer.data(gi + 196);
    const auto *gi_197 = buffer.data(gi + 197);
    const auto *gi_198 = buffer.data(gi + 198);
    const auto *gi_199 = buffer.data(gi + 199);
    const auto *gi_200 = buffer.data(gi + 200);
    const auto *gi_201 = buffer.data(gi + 201);
    const auto *gi_202 = buffer.data(gi + 202);
    const auto *gi_203 = buffer.data(gi + 203);
    const auto *gi_204 = buffer.data(gi + 204);
    const auto *gi_205 = buffer.data(gi + 205);
    const auto *gi_206 = buffer.data(gi + 206);
    const auto *gi_207 = buffer.data(gi + 207);
    const auto *gi_208 = buffer.data(gi + 208);
    const auto *gi_209 = buffer.data(gi + 209);
    const auto *gi_210 = buffer.data(gi + 210);
    const auto *gi_211 = buffer.data(gi + 211);
    const auto *gi_212 = buffer.data(gi + 212);
    const auto *gi_213 = buffer.data(gi + 213);
    const auto *gi_214 = buffer.data(gi + 214);
    const auto *gi_215 = buffer.data(gi + 215);
    const auto *gi_216 = buffer.data(gi + 216);
    const auto *gi_224 = buffer.data(gi + 224);
    const auto *gi_225 = buffer.data(gi + 225);
    const auto *gi_226 = buffer.data(gi + 226);
    const auto *gi_227 = buffer.data(gi + 227);
    const auto *gi_228 = buffer.data(gi + 228);
    const auto *gi_229 = buffer.data(gi + 229);
    const auto *gi_230 = buffer.data(gi + 230);
    const auto *gi_231 = buffer.data(gi + 231);
    const auto *gi_232 = buffer.data(gi + 232);
    const auto *gi_233 = buffer.data(gi + 233);
    const auto *gi_234 = buffer.data(gi + 234);
    const auto *gi_235 = buffer.data(gi + 235);
    const auto *gi_236 = buffer.data(gi + 236);
    const auto *gi_237 = buffer.data(gi + 237);
    const auto *gi_238 = buffer.data(gi + 238);
    const auto *gi_239 = buffer.data(gi + 239);
    const auto *gi_240 = buffer.data(gi + 240);
    const auto *gi_241 = buffer.data(gi + 241);
    const auto *gi_242 = buffer.data(gi + 242);
    const auto *gi_243 = buffer.data(gi + 243);
    const auto *gi_244 = buffer.data(gi + 244);
    const auto *gi_252 = buffer.data(gi + 252);
    const auto *gi_253 = buffer.data(gi + 253);
    const auto *gi_254 = buffer.data(gi + 254);
    const auto *gi_255 = buffer.data(gi + 255);
    const auto *gi_256 = buffer.data(gi + 256);
    const auto *gi_257 = buffer.data(gi + 257);
    const auto *gi_258 = buffer.data(gi + 258);
    const auto *gi_259 = buffer.data(gi + 259);
    const auto *gi_260 = buffer.data(gi + 260);
    const auto *gi_261 = buffer.data(gi + 261);
    const auto *gi_262 = buffer.data(gi + 262);
    const auto *gi_263 = buffer.data(gi + 263);
    const auto *gi_264 = buffer.data(gi + 264);
    const auto *gi_265 = buffer.data(gi + 265);
    const auto *gi_266 = buffer.data(gi + 266);
    const auto *gi_267 = buffer.data(gi + 267);
    const auto *gi_268 = buffer.data(gi + 268);
    const auto *gi_269 = buffer.data(gi + 269);
    const auto *gi_270 = buffer.data(gi + 270);
    const auto *gi_271 = buffer.data(gi + 271);
    const auto *gi_272 = buffer.data(gi + 272);
    const auto *gi_280 = buffer.data(gi + 280);
    const auto *gi_281 = buffer.data(gi + 281);
    const auto *gi_282 = buffer.data(gi + 282);
    const auto *gi_283 = buffer.data(gi + 283);
    const auto *gi_284 = buffer.data(gi + 284);
    const auto *gi_285 = buffer.data(gi + 285);
    const auto *gi_286 = buffer.data(gi + 286);
    const auto *gi_287 = buffer.data(gi + 287);
    const auto *gi_288 = buffer.data(gi + 288);
    const auto *gi_289 = buffer.data(gi + 289);
    const auto *gi_290 = buffer.data(gi + 290);
    const auto *gi_291 = buffer.data(gi + 291);
    const auto *gi_292 = buffer.data(gi + 292);
    const auto *gi_293 = buffer.data(gi + 293);
    const auto *gi_294 = buffer.data(gi + 294);
    const auto *gi_295 = buffer.data(gi + 295);
    const auto *gi_296 = buffer.data(gi + 296);
    const auto *gi_297 = buffer.data(gi + 297);
    const auto *gi_298 = buffer.data(gi + 298);
    const auto *gi_299 = buffer.data(gi + 299);
    const auto *gi_300 = buffer.data(gi + 300);
    const auto *gi_301 = buffer.data(gi + 301);
    const auto *gi_302 = buffer.data(gi + 302);
    const auto *gi_303 = buffer.data(gi + 303);
    const auto *gi_304 = buffer.data(gi + 304);
    const auto *gi_305 = buffer.data(gi + 305);
    const auto *gi_306 = buffer.data(gi + 306);
    const auto *gi_308 = buffer.data(gi + 308);
    const auto *gi_309 = buffer.data(gi + 309);
    const auto *gi_310 = buffer.data(gi + 310);
    const auto *gi_311 = buffer.data(gi + 311);
    const auto *gi_312 = buffer.data(gi + 312);
    const auto *gi_313 = buffer.data(gi + 313);
    const auto *gi_314 = buffer.data(gi + 314);
    const auto *gi_315 = buffer.data(gi + 315);
    const auto *gi_316 = buffer.data(gi + 316);
    const auto *gi_317 = buffer.data(gi + 317);
    const auto *gi_318 = buffer.data(gi + 318);
    const auto *gi_319 = buffer.data(gi + 319);
    const auto *gi_320 = buffer.data(gi + 320);
    const auto *gi_321 = buffer.data(gi + 321);
    const auto *gi_322 = buffer.data(gi + 322);
    const auto *gi_323 = buffer.data(gi + 323);
    const auto *gi_324 = buffer.data(gi + 324);
    const auto *gi_325 = buffer.data(gi + 325);
    const auto *gi_326 = buffer.data(gi + 326);
    const auto *gi_327 = buffer.data(gi + 327);
    const auto *gi_328 = buffer.data(gi + 328);
    const auto *gi_329 = buffer.data(gi + 329);
    const auto *gi_330 = buffer.data(gi + 330);
    const auto *gi_331 = buffer.data(gi + 331);
    const auto *gi_332 = buffer.data(gi + 332);
    const auto *gi_333 = buffer.data(gi + 333);
    const auto *gi_334 = buffer.data(gi + 334);
    const auto *gi_336 = buffer.data(gi + 336);
    const auto *gi_337 = buffer.data(gi + 337);
    const auto *gi_338 = buffer.data(gi + 338);
    const auto *gi_339 = buffer.data(gi + 339);
    const auto *gi_340 = buffer.data(gi + 340);
    const auto *gi_341 = buffer.data(gi + 341);
    const auto *gi_342 = buffer.data(gi + 342);
    const auto *gi_343 = buffer.data(gi + 343);
    const auto *gi_344 = buffer.data(gi + 344);
    const auto *gi_345 = buffer.data(gi + 345);
    const auto *gi_346 = buffer.data(gi + 346);
    const auto *gi_347 = buffer.data(gi + 347);
    const auto *gi_348 = buffer.data(gi + 348);
    const auto *gi_349 = buffer.data(gi + 349);
    const auto *gi_350 = buffer.data(gi + 350);
    const auto *gi_351 = buffer.data(gi + 351);
    const auto *gi_352 = buffer.data(gi + 352);
    const auto *gi_353 = buffer.data(gi + 353);
    const auto *gi_354 = buffer.data(gi + 354);
    const auto *gi_355 = buffer.data(gi + 355);
    const auto *gi_356 = buffer.data(gi + 356);
    const auto *gi_357 = buffer.data(gi + 357);
    const auto *gi_358 = buffer.data(gi + 358);
    const auto *gi_359 = buffer.data(gi + 359);
    const auto *gi_360 = buffer.data(gi + 360);
    const auto *gi_361 = buffer.data(gi + 361);
    const auto *gi_362 = buffer.data(gi + 362);
    const auto *gi_364 = buffer.data(gi + 364);
    const auto *gi_365 = buffer.data(gi + 365);
    const auto *gi_366 = buffer.data(gi + 366);
    const auto *gi_367 = buffer.data(gi + 367);
    const auto *gi_368 = buffer.data(gi + 368);
    const auto *gi_369 = buffer.data(gi + 369);
    const auto *gi_370 = buffer.data(gi + 370);
    const auto *gi_371 = buffer.data(gi + 371);
    const auto *gi_372 = buffer.data(gi + 372);
    const auto *gi_373 = buffer.data(gi + 373);
    const auto *gi_374 = buffer.data(gi + 374);
    const auto *gi_375 = buffer.data(gi + 375);
    const auto *gi_376 = buffer.data(gi + 376);
    const auto *gi_377 = buffer.data(gi + 377);
    const auto *gi_378 = buffer.data(gi + 378);
    const auto *gi_379 = buffer.data(gi + 379);
    const auto *gi_380 = buffer.data(gi + 380);
    const auto *gi_381 = buffer.data(gi + 381);
    const auto *gi_382 = buffer.data(gi + 382);
    const auto *gi_383 = buffer.data(gi + 383);
    const auto *gi_384 = buffer.data(gi + 384);
    const auto *gi_385 = buffer.data(gi + 385);
    const auto *gi_386 = buffer.data(gi + 386);
    const auto *gi_387 = buffer.data(gi + 387);
    const auto *gi_388 = buffer.data(gi + 388);
    const auto *gi_389 = buffer.data(gi + 389);
    const auto *gi_390 = buffer.data(gi + 390);
    const auto *gi_392 = buffer.data(gi + 392);
    const auto *gi_393 = buffer.data(gi + 393);
    const auto *gi_394 = buffer.data(gi + 394);
    const auto *gi_395 = buffer.data(gi + 395);
    const auto *gi_396 = buffer.data(gi + 396);
    const auto *gi_397 = buffer.data(gi + 397);
    const auto *gi_398 = buffer.data(gi + 398);
    const auto *gi_399 = buffer.data(gi + 399);
    const auto *gi_400 = buffer.data(gi + 400);
    const auto *gi_401 = buffer.data(gi + 401);
    const auto *gi_402 = buffer.data(gi + 402);
    const auto *gi_403 = buffer.data(gi + 403);
    const auto *gi_404 = buffer.data(gi + 404);
    const auto *gi_405 = buffer.data(gi + 405);
    const auto *gi_406 = buffer.data(gi + 406);
    const auto *gi_407 = buffer.data(gi + 407);
    const auto *gi_408 = buffer.data(gi + 408);
    const auto *gi_409 = buffer.data(gi + 409);
    const auto *gi_410 = buffer.data(gi + 410);
    const auto *gi_411 = buffer.data(gi + 411);
    const auto *gi_412 = buffer.data(gi + 412);
    const auto *gi_413 = buffer.data(gi + 413);
    const auto *gi_414 = buffer.data(gi + 414);
    const auto *gi_415 = buffer.data(gi + 415);
    const auto *gi_416 = buffer.data(gi + 416);
    const auto *gi_417 = buffer.data(gi + 417);
    const auto *gi_418 = buffer.data(gi + 418);
    const auto *gi_419 = buffer.data(gi + 419);

#pragma omp simd aligned(ab_x, ab_y, gh_22, gh_27, gh_36, gh_127, gh_132, gh_141, gh_211, \
                         gh_216, gh_225, gi_29, gi_34, gi_43, gi_169, gi_174, gi_183, gi_283, \
                         gi_290, gi_301 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -12.3046875 * ab_x[k] * gh_22[k]
                 + 24.609375 * ab_x[k] * gh_27[k]
                 - 2.4609375 * ab_x[k] * gh_36[k]
                 + 24.609375 * ab_x[k] * gh_127[k]
                 - 49.21875 * ab_x[k] * gh_132[k]
                 + 4.921875 * ab_x[k] * gh_141[k]
                 - 2.4609375 * ab_y[k] * gh_211[k]
                 + 4.921875 * ab_y[k] * gh_216[k]
                 - 0.4921875 * ab_y[k] * gh_225[k]
                 + 12.3046875 * gi_29[k]
                 - 24.609375 * gi_34[k]
                 + 2.4609375 * gi_43[k]
                 - 24.609375 * gi_169[k]
                 + 49.21875 * gi_174[k]
                 - 4.921875 * gi_183[k]
                 + 2.4609375 * gi_283[k]
                 - 4.921875 * gi_290[k]
                 + 0.4921875 * gi_301[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_25, gh_32, gh_130, gh_137, gh_214, gh_221, gi_32, \
                         gi_39, gi_172, gi_179, gi_287, gi_296 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = -f_0 * ab_x[k] * gh_25[k]
                 + f_0 * ab_x[k] * gh_32[k]
                 + f_1 * ab_x[k] * gh_130[k]
                 - f_1 * ab_x[k] * gh_137[k]
                 - f_2 * ab_y[k] * gh_214[k]
                 + f_2 * ab_y[k] * gh_221[k]
                 + f_0 * gi_32[k]
                 - f_0 * gi_39[k]
                 - f_1 * gi_172[k]
                 + f_1 * gi_179[k]
                 + f_2 * gi_287[k]
                 - f_2 * gi_296[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_22, gh_27, gh_29, gh_36, gh_38, gh_127, gh_132, \
                         gh_134, gh_141, gh_143, gh_211, gh_216, gh_218, gh_225, gh_227, \
                         gi_29, gi_34, gi_36, gi_43, gi_45, gi_169, gi_174, gi_176, gi_183, \
                         gi_185, gi_283, gi_290, gi_292, gi_301, \
                         gi_303 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = f_3 * ab_x[k] * gh_22[k]
                 + f_4 * ab_x[k] * gh_27[k]
                 - f_5 * ab_x[k] * gh_29[k]
                 - f_6 * ab_x[k] * gh_36[k]
                 + f_7 * ab_x[k] * gh_38[k]
                 - f_8 * ab_x[k] * gh_127[k]
                 - f_9 * ab_x[k] * gh_132[k]
                 + f_10 * ab_x[k] * gh_134[k]
                 + f_4 * ab_x[k] * gh_141[k]
                 - f_11 * ab_x[k] * gh_143[k]
                 + f_12 * ab_y[k] * gh_211[k]
                 + f_13 * ab_y[k] * gh_216[k]
                 - f_14 * ab_y[k] * gh_218[k]
                 - f_15 * ab_y[k] * gh_225[k]
                 + f_16 * ab_y[k] * gh_227[k]
                 - f_3 * gi_29[k]
                 - f_4 * gi_34[k]
                 + f_5 * gi_36[k]
                 + f_6 * gi_43[k]
                 - f_7 * gi_45[k]
                 + f_8 * gi_169[k]
                 + f_9 * gi_174[k]
                 - f_10 * gi_176[k]
                 - f_4 * gi_183[k]
                 + f_11 * gi_185[k]
                 - f_12 * gi_283[k]
                 - f_13 * gi_290[k]
                 + f_14 * gi_292[k]
                 + f_15 * gi_301[k]
                 - f_16 * gi_303[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_25, gh_32, gh_34, gh_130, gh_137, gh_139, gh_214, \
                         gh_221, gh_223, gi_32, gi_39, gi_41, gi_172, gi_179, gi_181, gi_287, \
                         gi_296, gi_298 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_17 * ab_x[k] * gh_25[k]
                 + f_17 * ab_x[k] * gh_32[k]
                 - f_18 * ab_x[k] * gh_34[k]
                 - f_18 * ab_x[k] * gh_130[k]
                 - f_18 * ab_x[k] * gh_137[k]
                 + f_19 * ab_x[k] * gh_139[k]
                 + f_20 * ab_y[k] * gh_214[k]
                 + f_20 * ab_y[k] * gh_221[k]
                 - f_21 * ab_y[k] * gh_223[k]
                 - f_17 * gi_32[k]
                 - f_17 * gi_39[k]
                 + f_18 * gi_41[k]
                 + f_18 * gi_172[k]
                 + f_18 * gi_179[k]
                 - f_19 * gi_181[k]
                 - f_20 * gi_287[k]
                 - f_20 * gi_296[k]
                 + f_21 * gi_298[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_22, gh_27, gh_29, gh_36, gh_38, gh_40, gh_127, gh_132, \
                         gh_134, gh_141, gh_143, gh_145, gh_211, gh_216, gh_218, gh_225, \
                         gh_227, gh_229, gi_29, gi_34, gi_36, gi_43, gi_45, gi_47, gi_169, \
                         gi_174, gi_176, gi_183, gi_185, gi_187, gi_283, gi_290, gi_292, \
                         gi_301, gi_303, gi_305 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_22 * ab_x[k] * gh_22[k]
                 - f_23 * ab_x[k] * gh_27[k]
                 + f_24 * ab_x[k] * gh_29[k]
                 - f_22 * ab_x[k] * gh_36[k]
                 + f_24 * ab_x[k] * gh_38[k]
                 - f_25 * ab_x[k] * gh_40[k]
                 + f_23 * ab_x[k] * gh_127[k]
                 + f_26 * ab_x[k] * gh_132[k]
                 - f_27 * ab_x[k] * gh_134[k]
                 + f_23 * ab_x[k] * gh_141[k]
                 - f_27 * ab_x[k] * gh_143[k]
                 + f_28 * ab_x[k] * gh_145[k]
                 - f_29 * ab_y[k] * gh_211[k]
                 - f_30 * ab_y[k] * gh_216[k]
                 + f_31 * ab_y[k] * gh_218[k]
                 - f_29 * ab_y[k] * gh_225[k]
                 + f_31 * ab_y[k] * gh_227[k]
                 - f_32 * ab_y[k] * gh_229[k]
                 + f_22 * gi_29[k]
                 + f_23 * gi_34[k]
                 - f_24 * gi_36[k]
                 + f_22 * gi_43[k]
                 - f_24 * gi_45[k]
                 + f_25 * gi_47[k]
                 - f_23 * gi_169[k]
                 - f_26 * gi_174[k]
                 + f_27 * gi_176[k]
                 - f_23 * gi_183[k]
                 + f_27 * gi_185[k]
                 - f_28 * gi_187[k]
                 + f_29 * gi_283[k]
                 + f_30 * gi_290[k]
                 - f_31 * gi_292[k]
                 + f_29 * gi_301[k]
                 - f_31 * gi_303[k]
                 + f_32 * gi_305[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_28, gh_30, gh_37, gh_39, gh_41, gh_128, gh_133, \
                         gh_135, gh_142, gh_144, gh_146, gh_212, gh_217, gh_219, gh_226, \
                         gh_228, gh_230, gi_30, gi_35, gi_37, gi_44, gi_46, gi_48, gi_170, \
                         gi_175, gi_177, gi_184, gi_186, gi_188, gi_284, gi_291, gi_293, \
                         gi_302, gi_304, gi_306 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -f_33 * ab_x[k] * gh_23[k]
                 - f_34 * ab_x[k] * gh_28[k]
                 + f_35 * ab_x[k] * gh_30[k]
                 - f_33 * ab_x[k] * gh_37[k]
                 + f_35 * ab_x[k] * gh_39[k]
                 - f_36 * ab_x[k] * gh_41[k]
                 + f_34 * ab_x[k] * gh_128[k]
                 + f_37 * ab_x[k] * gh_133[k]
                 - f_38 * ab_x[k] * gh_135[k]
                 + f_34 * ab_x[k] * gh_142[k]
                 - f_38 * ab_x[k] * gh_144[k]
                 + f_39 * ab_x[k] * gh_146[k]
                 - f_40 * ab_y[k] * gh_212[k]
                 - f_41 * ab_y[k] * gh_217[k]
                 + f_36 * ab_y[k] * gh_219[k]
                 - f_40 * ab_y[k] * gh_226[k]
                 + f_36 * ab_y[k] * gh_228[k]
                 - f_42 * ab_y[k] * gh_230[k]
                 + f_33 * gi_30[k]
                 + f_34 * gi_35[k]
                 - f_35 * gi_37[k]
                 + f_33 * gi_44[k]
                 - f_35 * gi_46[k]
                 + f_36 * gi_48[k]
                 - f_34 * gi_170[k]
                 - f_37 * gi_175[k]
                 + f_38 * gi_177[k]
                 - f_34 * gi_184[k]
                 + f_38 * gi_186[k]
                 - f_39 * gi_188[k]
                 + f_40 * gi_284[k]
                 + f_41 * gi_291[k]
                 - f_36 * gi_293[k]
                 + f_40 * gi_302[k]
                 - f_36 * gi_304[k]
                 + f_42 * gi_306[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_26, gh_31, gh_33, gh_35, gh_126, gh_129, \
                         gh_131, gh_136, gh_138, gh_140, gh_210, gh_213, gh_215, gh_220, \
                         gh_222, gh_224, gi_28, gi_31, gi_33, gi_38, gi_40, gi_42, gi_168, \
                         gi_171, gi_173, gi_178, gi_180, gi_182, gi_281, gi_286, gi_288, \
                         gi_295, gi_297, gi_299 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_22 * ab_x[k] * gh_21[k]
                 - f_23 * ab_x[k] * gh_24[k]
                 + f_24 * ab_x[k] * gh_26[k]
                 - f_22 * ab_x[k] * gh_31[k]
                 + f_24 * ab_x[k] * gh_33[k]
                 - f_25 * ab_x[k] * gh_35[k]
                 + f_23 * ab_x[k] * gh_126[k]
                 + f_26 * ab_x[k] * gh_129[k]
                 - f_27 * ab_x[k] * gh_131[k]
                 + f_23 * ab_x[k] * gh_136[k]
                 - f_27 * ab_x[k] * gh_138[k]
                 + f_28 * ab_x[k] * gh_140[k]
                 - f_29 * ab_y[k] * gh_210[k]
                 - f_30 * ab_y[k] * gh_213[k]
                 + f_31 * ab_y[k] * gh_215[k]
                 - f_29 * ab_y[k] * gh_220[k]
                 + f_31 * ab_y[k] * gh_222[k]
                 - f_32 * ab_y[k] * gh_224[k]
                 + f_22 * gi_28[k]
                 + f_23 * gi_31[k]
                 - f_24 * gi_33[k]
                 + f_22 * gi_38[k]
                 - f_24 * gi_40[k]
                 + f_25 * gi_42[k]
                 - f_23 * gi_168[k]
                 - f_26 * gi_171[k]
                 + f_27 * gi_173[k]
                 - f_23 * gi_178[k]
                 + f_27 * gi_180[k]
                 - f_28 * gi_182[k]
                 + f_29 * gi_281[k]
                 + f_30 * gi_286[k]
                 - f_31 * gi_288[k]
                 + f_29 * gi_295[k]
                 - f_31 * gi_297[k]
                 + f_32 * gi_299[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_30, gh_37, gh_39, gh_128, gh_135, gh_142, \
                         gh_144, gh_212, gh_219, gh_226, gh_228, gi_30, gi_37, gi_44, gi_46, \
                         gi_170, gi_177, gi_184, gi_186, gi_284, gi_293, gi_302, \
                         gi_304 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_43 * ab_x[k] * gh_23[k]
                 - f_17 * ab_x[k] * gh_30[k]
                 - f_43 * ab_x[k] * gh_37[k]
                 + f_17 * ab_x[k] * gh_39[k]
                 - f_17 * ab_x[k] * gh_128[k]
                 + f_18 * ab_x[k] * gh_135[k]
                 + f_17 * ab_x[k] * gh_142[k]
                 - f_18 * ab_x[k] * gh_144[k]
                 + f_44 * ab_y[k] * gh_212[k]
                 - f_20 * ab_y[k] * gh_219[k]
                 - f_44 * ab_y[k] * gh_226[k]
                 + f_20 * ab_y[k] * gh_228[k]
                 - f_43 * gi_30[k]
                 + f_17 * gi_37[k]
                 + f_43 * gi_44[k]
                 - f_17 * gi_46[k]
                 + f_17 * gi_170[k]
                 - f_18 * gi_177[k]
                 - f_17 * gi_184[k]
                 + f_18 * gi_186[k]
                 - f_44 * gi_284[k]
                 + f_20 * gi_293[k]
                 + f_44 * gi_302[k]
                 - f_20 * gi_304[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_26, gh_31, gh_33, gh_126, gh_129, \
                         gh_131, gh_136, gh_138, gh_210, gh_213, gh_215, gh_220, gh_222, \
                         gi_28, gi_31, gi_33, gi_38, gi_40, gi_168, gi_171, gi_173, gi_178, \
                         gi_180, gi_281, gi_286, gi_288, gi_295, \
                         gi_297 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_6 * ab_x[k] * gh_21[k]
                 - f_4 * ab_x[k] * gh_24[k]
                 - f_7 * ab_x[k] * gh_26[k]
                 - f_3 * ab_x[k] * gh_31[k]
                 + f_5 * ab_x[k] * gh_33[k]
                 - f_4 * ab_x[k] * gh_126[k]
                 + f_9 * ab_x[k] * gh_129[k]
                 + f_11 * ab_x[k] * gh_131[k]
                 + f_8 * ab_x[k] * gh_136[k]
                 - f_10 * ab_x[k] * gh_138[k]
                 + f_15 * ab_y[k] * gh_210[k]
                 - f_13 * ab_y[k] * gh_213[k]
                 - f_16 * ab_y[k] * gh_215[k]
                 - f_12 * ab_y[k] * gh_220[k]
                 + f_14 * ab_y[k] * gh_222[k]
                 - f_6 * gi_28[k]
                 + f_4 * gi_31[k]
                 + f_7 * gi_33[k]
                 + f_3 * gi_38[k]
                 - f_5 * gi_40[k]
                 + f_4 * gi_168[k]
                 - f_9 * gi_171[k]
                 - f_11 * gi_173[k]
                 - f_8 * gi_178[k]
                 + f_10 * gi_180[k]
                 - f_15 * gi_281[k]
                 + f_13 * gi_286[k]
                 + f_16 * gi_288[k]
                 + f_12 * gi_295[k]
                 - f_14 * gi_297[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_28, gh_37, gh_128, gh_133, gh_142, gh_212, \
                         gh_217, gh_226, gi_30, gi_35, gi_44, gi_170, gi_175, gi_184, gi_284, \
                         gi_291, gi_302 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_45 * ab_x[k] * gh_23[k]
                 + f_46 * ab_x[k] * gh_28[k]
                 - f_45 * ab_x[k] * gh_37[k]
                 + f_47 * ab_x[k] * gh_128[k]
                 - f_48 * ab_x[k] * gh_133[k]
                 + f_47 * ab_x[k] * gh_142[k]
                 - f_49 * ab_y[k] * gh_212[k]
                 + f_50 * ab_y[k] * gh_217[k]
                 - f_49 * ab_y[k] * gh_226[k]
                 + f_45 * gi_30[k]
                 - f_46 * gi_35[k]
                 + f_45 * gi_44[k]
                 - f_47 * gi_170[k]
                 + f_48 * gi_175[k]
                 - f_47 * gi_184[k]
                 + f_49 * gi_284[k]
                 - f_50 * gi_291[k]
                 + f_49 * gi_302[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_31, gh_126, gh_129, gh_136, gh_210, \
                         gh_213, gh_220, gi_28, gi_31, gi_38, gi_168, gi_171, gi_178, gi_281, \
                         gi_286, gi_295 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -2.4609375 * ab_x[k] * gh_21[k]
                  + 24.609375 * ab_x[k] * gh_24[k]
                  - 12.3046875 * ab_x[k] * gh_31[k]
                  + 4.921875 * ab_x[k] * gh_126[k]
                  - 49.21875 * ab_x[k] * gh_129[k]
                  + 24.609375 * ab_x[k] * gh_136[k]
                  - 0.4921875 * ab_y[k] * gh_210[k]
                  + 4.921875 * ab_y[k] * gh_213[k]
                  - 2.4609375 * ab_y[k] * gh_220[k]
                  + 2.4609375 * gi_28[k]
                  - 24.609375 * gi_31[k]
                  + 12.3046875 * gi_38[k]
                  - 4.921875 * gi_168[k]
                  + 49.21875 * gi_171[k]
                  - 24.609375 * gi_178[k]
                  + 0.4921875 * gi_281[k]
                  - 4.921875 * gi_286[k]
                  + 2.4609375 * gi_295[k];
    }

#pragma omp simd aligned(ab_x, gh_85, gh_90, gh_99, gh_232, gh_237, gh_246, gi_113, gi_118, \
                         gi_127, gi_309, gi_314, gi_323 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_0 * ab_x[k] * gh_85[k]
                  + f_1 * ab_x[k] * gh_90[k]
                  - f_2 * ab_x[k] * gh_99[k]
                  + f_0 * ab_x[k] * gh_232[k]
                  - f_1 * ab_x[k] * gh_237[k]
                  + f_2 * ab_x[k] * gh_246[k]
                  + f_0 * gi_113[k]
                  - f_1 * gi_118[k]
                  + f_2 * gi_127[k]
                  - f_0 * gi_309[k]
                  + f_1 * gi_314[k]
                  - f_2 * gi_323[k];
    }

#pragma omp simd aligned(ab_x, gh_88, gh_95, gh_235, gh_242, gi_116, gi_123, gi_312, \
                         gi_319 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = -78.75 * ab_x[k] * gh_88[k]
                  + 78.75 * ab_x[k] * gh_95[k]
                  + 78.75 * ab_x[k] * gh_235[k]
                  - 78.75 * ab_x[k] * gh_242[k]
                  + 78.75 * gi_116[k]
                  - 78.75 * gi_123[k]
                  - 78.75 * gi_312[k]
                  + 78.75 * gi_319[k];
    }

#pragma omp simd aligned(ab_x, gh_85, gh_90, gh_92, gh_99, gh_101, gh_232, gh_237, gh_239, \
                         gh_246, gh_248, gi_113, gi_118, gi_120, gi_127, gi_129, gi_309, \
                         gi_314, gi_316, gi_323, gi_325 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_51 * ab_x[k] * gh_85[k]
                  + f_52 * ab_x[k] * gh_90[k]
                  - f_53 * ab_x[k] * gh_92[k]
                  - f_54 * ab_x[k] * gh_99[k]
                  + f_55 * ab_x[k] * gh_101[k]
                  - f_51 * ab_x[k] * gh_232[k]
                  - f_52 * ab_x[k] * gh_237[k]
                  + f_53 * ab_x[k] * gh_239[k]
                  + f_54 * ab_x[k] * gh_246[k]
                  - f_55 * ab_x[k] * gh_248[k]
                  - f_51 * gi_113[k]
                  - f_52 * gi_118[k]
                  + f_53 * gi_120[k]
                  + f_54 * gi_127[k]
                  - f_55 * gi_129[k]
                  + f_51 * gi_309[k]
                  + f_52 * gi_314[k]
                  - f_53 * gi_316[k]
                  - f_54 * gi_323[k]
                  + f_55 * gi_325[k];
    }

#pragma omp simd aligned(ab_x, gh_88, gh_95, gh_97, gh_235, gh_242, gh_244, gi_116, gi_123, \
                         gi_125, gi_312, gi_319, gi_321 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_56 * ab_x[k] * gh_88[k]
                  + f_56 * ab_x[k] * gh_95[k]
                  - f_57 * ab_x[k] * gh_97[k]
                  - f_56 * ab_x[k] * gh_235[k]
                  - f_56 * ab_x[k] * gh_242[k]
                  + f_57 * ab_x[k] * gh_244[k]
                  - f_56 * gi_116[k]
                  - f_56 * gi_123[k]
                  + f_57 * gi_125[k]
                  + f_56 * gi_312[k]
                  + f_56 * gi_319[k]
                  - f_57 * gi_321[k];
    }

#pragma omp simd aligned(ab_x, gh_85, gh_90, gh_92, gh_99, gh_101, gh_103, gh_232, gh_237, \
                         gh_239, gh_246, gh_248, gh_250, gi_113, gi_118, gi_120, gi_127, \
                         gi_129, gi_131, gi_309, gi_314, gi_316, gi_323, gi_325, \
                         gi_327 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_58 * ab_x[k] * gh_85[k]
                  - f_59 * ab_x[k] * gh_90[k]
                  + f_60 * ab_x[k] * gh_92[k]
                  - f_58 * ab_x[k] * gh_99[k]
                  + f_60 * ab_x[k] * gh_101[k]
                  - f_61 * ab_x[k] * gh_103[k]
                  + f_58 * ab_x[k] * gh_232[k]
                  + f_59 * ab_x[k] * gh_237[k]
                  - f_60 * ab_x[k] * gh_239[k]
                  + f_58 * ab_x[k] * gh_246[k]
                  - f_60 * ab_x[k] * gh_248[k]
                  + f_61 * ab_x[k] * gh_250[k]
                  + f_58 * gi_113[k]
                  + f_59 * gi_118[k]
                  - f_60 * gi_120[k]
                  + f_58 * gi_127[k]
                  - f_60 * gi_129[k]
                  + f_61 * gi_131[k]
                  - f_58 * gi_309[k]
                  - f_59 * gi_314[k]
                  + f_60 * gi_316[k]
                  - f_58 * gi_323[k]
                  + f_60 * gi_325[k]
                  - f_61 * gi_327[k];
    }

#pragma omp simd aligned(ab_x, gh_86, gh_91, gh_93, gh_100, gh_102, gh_104, gh_233, gh_238, \
                         gh_240, gh_247, gh_249, gh_251, gi_114, gi_119, gi_121, gi_128, \
                         gi_130, gi_132, gi_310, gi_315, gi_317, gi_324, gi_326, \
                         gi_328 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = -f_62 * ab_x[k] * gh_86[k]
                  - f_63 * ab_x[k] * gh_91[k]
                  + f_64 * ab_x[k] * gh_93[k]
                  - f_62 * ab_x[k] * gh_100[k]
                  + f_64 * ab_x[k] * gh_102[k]
                  - f_65 * ab_x[k] * gh_104[k]
                  + f_62 * ab_x[k] * gh_233[k]
                  + f_63 * ab_x[k] * gh_238[k]
                  - f_64 * ab_x[k] * gh_240[k]
                  + f_62 * ab_x[k] * gh_247[k]
                  - f_64 * ab_x[k] * gh_249[k]
                  + f_65 * ab_x[k] * gh_251[k]
                  + f_62 * gi_114[k]
                  + f_63 * gi_119[k]
                  - f_64 * gi_121[k]
                  + f_62 * gi_128[k]
                  - f_64 * gi_130[k]
                  + f_65 * gi_132[k]
                  - f_62 * gi_310[k]
                  - f_63 * gi_315[k]
                  + f_64 * gi_317[k]
                  - f_62 * gi_324[k]
                  + f_64 * gi_326[k]
                  - f_65 * gi_328[k];
    }

#pragma omp simd aligned(ab_x, gh_84, gh_87, gh_89, gh_94, gh_96, gh_98, gh_231, gh_234, \
                         gh_236, gh_241, gh_243, gh_245, gi_112, gi_115, gi_117, gi_122, \
                         gi_124, gi_126, gi_308, gi_311, gi_313, gi_318, gi_320, \
                         gi_322 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_58 * ab_x[k] * gh_84[k]
                  - f_59 * ab_x[k] * gh_87[k]
                  + f_60 * ab_x[k] * gh_89[k]
                  - f_58 * ab_x[k] * gh_94[k]
                  + f_60 * ab_x[k] * gh_96[k]
                  - f_61 * ab_x[k] * gh_98[k]
                  + f_58 * ab_x[k] * gh_231[k]
                  + f_59 * ab_x[k] * gh_234[k]
                  - f_60 * ab_x[k] * gh_236[k]
                  + f_58 * ab_x[k] * gh_241[k]
                  - f_60 * ab_x[k] * gh_243[k]
                  + f_61 * ab_x[k] * gh_245[k]
                  + f_58 * gi_112[k]
                  + f_59 * gi_115[k]
                  - f_60 * gi_117[k]
                  + f_58 * gi_122[k]
                  - f_60 * gi_124[k]
                  + f_61 * gi_126[k]
                  - f_58 * gi_308[k]
                  - f_59 * gi_311[k]
                  + f_60 * gi_313[k]
                  - f_58 * gi_318[k]
                  + f_60 * gi_320[k]
                  - f_61 * gi_322[k];
    }

#pragma omp simd aligned(ab_x, gh_86, gh_93, gh_100, gh_102, gh_233, gh_240, gh_247, gh_249, \
                         gi_114, gi_121, gi_128, gi_130, gi_310, gi_317, gi_324, \
                         gi_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_66 * ab_x[k] * gh_86[k]
                  - f_56 * ab_x[k] * gh_93[k]
                  - f_66 * ab_x[k] * gh_100[k]
                  + f_56 * ab_x[k] * gh_102[k]
                  - f_66 * ab_x[k] * gh_233[k]
                  + f_56 * ab_x[k] * gh_240[k]
                  + f_66 * ab_x[k] * gh_247[k]
                  - f_56 * ab_x[k] * gh_249[k]
                  - f_66 * gi_114[k]
                  + f_56 * gi_121[k]
                  + f_66 * gi_128[k]
                  - f_56 * gi_130[k]
                  + f_66 * gi_310[k]
                  - f_56 * gi_317[k]
                  - f_66 * gi_324[k]
                  + f_56 * gi_326[k];
    }

#pragma omp simd aligned(ab_x, gh_84, gh_87, gh_89, gh_94, gh_96, gh_231, gh_234, gh_236, \
                         gh_241, gh_243, gi_112, gi_115, gi_117, gi_122, gi_124, gi_308, \
                         gi_311, gi_313, gi_318, gi_320 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_54 * ab_x[k] * gh_84[k]
                  - f_52 * ab_x[k] * gh_87[k]
                  - f_55 * ab_x[k] * gh_89[k]
                  - f_51 * ab_x[k] * gh_94[k]
                  + f_53 * ab_x[k] * gh_96[k]
                  - f_54 * ab_x[k] * gh_231[k]
                  + f_52 * ab_x[k] * gh_234[k]
                  + f_55 * ab_x[k] * gh_236[k]
                  + f_51 * ab_x[k] * gh_241[k]
                  - f_53 * ab_x[k] * gh_243[k]
                  - f_54 * gi_112[k]
                  + f_52 * gi_115[k]
                  + f_55 * gi_117[k]
                  + f_51 * gi_122[k]
                  - f_53 * gi_124[k]
                  + f_54 * gi_308[k]
                  - f_52 * gi_311[k]
                  - f_55 * gi_313[k]
                  - f_51 * gi_318[k]
                  + f_53 * gi_320[k];
    }

#pragma omp simd aligned(ab_x, gh_86, gh_91, gh_100, gh_233, gh_238, gh_247, gi_114, gi_119, \
                         gi_128, gi_310, gi_315, gi_324 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -19.6875 * ab_x[k] * gh_86[k]
                  + 118.125 * ab_x[k] * gh_91[k]
                  - 19.6875 * ab_x[k] * gh_100[k]
                  + 19.6875 * ab_x[k] * gh_233[k]
                  - 118.125 * ab_x[k] * gh_238[k]
                  + 19.6875 * ab_x[k] * gh_247[k]
                  + 19.6875 * gi_114[k]
                  - 118.125 * gi_119[k]
                  + 19.6875 * gi_128[k]
                  - 19.6875 * gi_310[k]
                  + 118.125 * gi_315[k]
                  - 19.6875 * gi_324[k];
    }

#pragma omp simd aligned(ab_x, gh_84, gh_87, gh_94, gh_231, gh_234, gh_241, gi_112, gi_115, \
                         gi_122, gi_308, gi_311, gi_318 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_2 * ab_x[k] * gh_84[k]
                  + f_1 * ab_x[k] * gh_87[k]
                  - f_0 * ab_x[k] * gh_94[k]
                  + f_2 * ab_x[k] * gh_231[k]
                  - f_1 * ab_x[k] * gh_234[k]
                  + f_0 * ab_x[k] * gh_241[k]
                  + f_2 * gi_112[k]
                  - f_1 * gi_115[k]
                  + f_0 * gi_122[k]
                  - f_2 * gi_308[k]
                  + f_1 * gi_311[k]
                  - f_0 * gi_318[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_22, gh_27, gh_36, gh_127, gh_132, gh_141, gh_169, \
                         gh_174, gh_183, gh_211, gh_216, gh_225, gh_253, gh_258, gh_267, \
                         gi_29, gi_34, gi_43, gi_169, gi_174, gi_183, gi_225, gi_230, gi_239, \
                         gi_283, gi_290, gi_301, gi_339, gi_346, \
                         gi_357 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_3 * ab_x[k] * gh_22[k]
                  - f_8 * ab_x[k] * gh_27[k]
                  + f_12 * ab_x[k] * gh_36[k]
                  + f_4 * ab_x[k] * gh_127[k]
                  - f_9 * ab_x[k] * gh_132[k]
                  + f_13 * ab_x[k] * gh_141[k]
                  - f_5 * ab_x[k] * gh_169[k]
                  + f_10 * ab_x[k] * gh_174[k]
                  - f_14 * ab_x[k] * gh_183[k]
                  - f_6 * ab_y[k] * gh_211[k]
                  + f_4 * ab_y[k] * gh_216[k]
                  - f_15 * ab_y[k] * gh_225[k]
                  + f_7 * ab_y[k] * gh_253[k]
                  - f_11 * ab_y[k] * gh_258[k]
                  + f_16 * ab_y[k] * gh_267[k]
                  - f_3 * gi_29[k]
                  + f_8 * gi_34[k]
                  - f_12 * gi_43[k]
                  - f_4 * gi_169[k]
                  + f_9 * gi_174[k]
                  - f_13 * gi_183[k]
                  + f_5 * gi_225[k]
                  - f_10 * gi_230[k]
                  + f_14 * gi_239[k]
                  + f_6 * gi_283[k]
                  - f_4 * gi_290[k]
                  + f_15 * gi_301[k]
                  - f_7 * gi_339[k]
                  + f_11 * gi_346[k]
                  - f_16 * gi_357[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_25, gh_32, gh_130, gh_137, gh_172, gh_179, gh_214, \
                         gh_221, gh_256, gh_263, gi_32, gi_39, gi_172, gi_179, gi_228, gi_235, \
                         gi_287, gi_296, gi_343, gi_352 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_51 * ab_x[k] * gh_25[k]
                  - f_51 * ab_x[k] * gh_32[k]
                  + f_52 * ab_x[k] * gh_130[k]
                  - f_52 * ab_x[k] * gh_137[k]
                  - f_53 * ab_x[k] * gh_172[k]
                  + f_53 * ab_x[k] * gh_179[k]
                  - f_54 * ab_y[k] * gh_214[k]
                  + f_54 * ab_y[k] * gh_221[k]
                  + f_55 * ab_y[k] * gh_256[k]
                  - f_55 * ab_y[k] * gh_263[k]
                  - f_51 * gi_32[k]
                  + f_51 * gi_39[k]
                  - f_52 * gi_172[k]
                  + f_52 * gi_179[k]
                  + f_53 * gi_228[k]
                  - f_53 * gi_235[k]
                  + f_54 * gi_287[k]
                  - f_54 * gi_296[k]
                  - f_55 * gi_343[k]
                  + f_55 * gi_352[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_22, gh_27, gh_29, gh_36, gh_38, gh_127, gh_132, \
                         gh_134, gh_141, gh_143, gh_169, gh_174, gh_176, gh_183, gh_185, \
                         gh_211, gh_216, gh_218, gh_225, gh_227, gh_253, gh_258, gh_260, \
                         gh_267, gh_269, gi_29, gi_34, gi_36, gi_43, gi_45, gi_169, gi_174, \
                         gi_176, gi_183, gi_185, gi_225, gi_230, gi_232, gi_239, gi_241, \
                         gi_283, gi_290, gi_292, gi_301, gi_303, gi_339, gi_346, gi_348, \
                         gi_357, gi_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -2.4609375 * ab_x[k] * gh_22[k]
                  - 1.640625 * ab_x[k] * gh_27[k]
                  + 19.6875 * ab_x[k] * gh_29[k]
                  + 0.8203125 * ab_x[k] * gh_36[k]
                  - 6.5625 * ab_x[k] * gh_38[k]
                  - 1.640625 * ab_x[k] * gh_127[k]
                  - 1.09375 * ab_x[k] * gh_132[k]
                  + 13.125 * ab_x[k] * gh_134[k]
                  + 0.546875 * ab_x[k] * gh_141[k]
                  - 4.375 * ab_x[k] * gh_143[k]
                  + 19.6875 * ab_x[k] * gh_169[k]
                  + 13.125 * ab_x[k] * gh_174[k]
                  - 157.5 * ab_x[k] * gh_176[k]
                  - 6.5625 * ab_x[k] * gh_183[k]
                  + 52.5 * ab_x[k] * gh_185[k]
                  + 0.8203125 * ab_y[k] * gh_211[k]
                  + 0.546875 * ab_y[k] * gh_216[k]
                  - 6.5625 * ab_y[k] * gh_218[k]
                  - 0.2734375 * ab_y[k] * gh_225[k]
                  + 2.1875 * ab_y[k] * gh_227[k]
                  - 6.5625 * ab_y[k] * gh_253[k]
                  - 4.375 * ab_y[k] * gh_258[k]
                  + 52.5 * ab_y[k] * gh_260[k]
                  + 2.1875 * ab_y[k] * gh_267[k]
                  - 17.5 * ab_y[k] * gh_269[k]
                  + 2.4609375 * gi_29[k]
                  + 1.640625 * gi_34[k]
                  - 19.6875 * gi_36[k]
                  - 0.8203125 * gi_43[k]
                  + 6.5625 * gi_45[k]
                  + 1.640625 * gi_169[k]
                  + 1.09375 * gi_174[k]
                  - 13.125 * gi_176[k]
                  - 0.546875 * gi_183[k]
                  + 4.375 * gi_185[k]
                  - 19.6875 * gi_225[k]
                  - 13.125 * gi_230[k]
                  + 157.5 * gi_232[k]
                  + 6.5625 * gi_239[k]
                  - 52.5 * gi_241[k]
                  - 0.8203125 * gi_283[k]
                  - 0.546875 * gi_290[k]
                  + 6.5625 * gi_292[k]
                  + 0.2734375 * gi_301[k]
                  - 2.1875 * gi_303[k]
                  + 6.5625 * gi_339[k]
                  + 4.375 * gi_346[k]
                  - 52.5 * gi_348[k]
                  - 2.1875 * gi_357[k]
                  + 17.5 * gi_359[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_25, gh_32, gh_34, gh_130, gh_137, gh_139, gh_172, \
                         gh_179, gh_181, gh_214, gh_221, gh_223, gh_256, gh_263, gh_265, \
                         gi_32, gi_39, gi_41, gi_172, gi_179, gi_181, gi_228, gi_235, gi_237, \
                         gi_287, gi_296, gi_298, gi_343, gi_352, \
                         gi_354 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_67 * ab_x[k] * gh_25[k]
                  - f_67 * ab_x[k] * gh_32[k]
                  + f_68 * ab_x[k] * gh_34[k]
                  - f_69 * ab_x[k] * gh_130[k]
                  - f_69 * ab_x[k] * gh_137[k]
                  + f_70 * ab_x[k] * gh_139[k]
                  + f_71 * ab_x[k] * gh_172[k]
                  + f_71 * ab_x[k] * gh_179[k]
                  - f_72 * ab_x[k] * gh_181[k]
                  + f_73 * ab_y[k] * gh_214[k]
                  + f_73 * ab_y[k] * gh_221[k]
                  - f_69 * ab_y[k] * gh_223[k]
                  - f_74 * ab_y[k] * gh_256[k]
                  - f_74 * ab_y[k] * gh_263[k]
                  + f_75 * ab_y[k] * gh_265[k]
                  + f_67 * gi_32[k]
                  + f_67 * gi_39[k]
                  - f_68 * gi_41[k]
                  + f_69 * gi_172[k]
                  + f_69 * gi_179[k]
                  - f_70 * gi_181[k]
                  - f_71 * gi_228[k]
                  - f_71 * gi_235[k]
                  + f_72 * gi_237[k]
                  - f_73 * gi_287[k]
                  - f_73 * gi_296[k]
                  + f_69 * gi_298[k]
                  + f_74 * gi_343[k]
                  + f_74 * gi_352[k]
                  - f_75 * gi_354[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_22, gh_27, gh_29, gh_36, gh_38, gh_40, gh_127, gh_132, \
                         gh_134, gh_141, gh_143, gh_145, gh_169, gh_174, gh_176, gh_183, \
                         gh_185, gh_187, gh_211, gh_216, gh_218, gh_225, gh_227, gh_229, \
                         gh_253, gh_258, gh_260, gh_267, gh_269, gh_271, gi_29, gi_34, gi_36, \
                         gi_43, gi_45, gi_47, gi_169, gi_174, gi_176, gi_183, gi_185, gi_187, \
                         gi_225, gi_230, gi_232, gi_239, gi_241, gi_243, gi_283, gi_290, \
                         gi_292, gi_301, gi_303, gi_305, gi_339, gi_346, gi_348, gi_357, \
                         gi_359, gi_361 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = f_76 * ab_x[k] * gh_22[k]
                  + f_77 * ab_x[k] * gh_27[k]
                  - f_78 * ab_x[k] * gh_29[k]
                  + f_76 * ab_x[k] * gh_36[k]
                  - f_78 * ab_x[k] * gh_38[k]
                  + f_79 * ab_x[k] * gh_40[k]
                  + f_80 * ab_x[k] * gh_127[k]
                  + f_81 * ab_x[k] * gh_132[k]
                  - f_79 * ab_x[k] * gh_134[k]
                  + f_80 * ab_x[k] * gh_141[k]
                  - f_79 * ab_x[k] * gh_143[k]
                  + f_82 * ab_x[k] * gh_145[k]
                  - f_79 * ab_x[k] * gh_169[k]
                  - f_83 * ab_x[k] * gh_174[k]
                  + f_84 * ab_x[k] * gh_176[k]
                  - f_79 * ab_x[k] * gh_183[k]
                  + f_84 * ab_x[k] * gh_185[k]
                  - f_85 * ab_x[k] * gh_187[k]
                  - f_86 * ab_y[k] * gh_211[k]
                  - f_80 * ab_y[k] * gh_216[k]
                  + f_87 * ab_y[k] * gh_218[k]
                  - f_86 * ab_y[k] * gh_225[k]
                  + f_87 * ab_y[k] * gh_227[k]
                  - f_88 * ab_y[k] * gh_229[k]
                  + f_88 * ab_y[k] * gh_253[k]
                  + f_82 * ab_y[k] * gh_258[k]
                  - f_89 * ab_y[k] * gh_260[k]
                  + f_88 * ab_y[k] * gh_267[k]
                  - f_89 * ab_y[k] * gh_269[k]
                  + f_90 * ab_y[k] * gh_271[k]
                  - f_76 * gi_29[k]
                  - f_77 * gi_34[k]
                  + f_78 * gi_36[k]
                  - f_76 * gi_43[k]
                  + f_78 * gi_45[k]
                  - f_79 * gi_47[k]
                  - f_80 * gi_169[k]
                  - f_81 * gi_174[k]
                  + f_79 * gi_176[k]
                  - f_80 * gi_183[k]
                  + f_79 * gi_185[k]
                  - f_82 * gi_187[k]
                  + f_79 * gi_225[k]
                  + f_83 * gi_230[k]
                  - f_84 * gi_232[k]
                  + f_79 * gi_239[k]
                  - f_84 * gi_241[k]
                  + f_85 * gi_243[k]
                  + f_86 * gi_283[k]
                  + f_80 * gi_290[k]
                  - f_87 * gi_292[k]
                  + f_86 * gi_301[k]
                  - f_87 * gi_303[k]
                  + f_88 * gi_305[k]
                  - f_88 * gi_339[k]
                  - f_82 * gi_346[k]
                  + f_89 * gi_348[k]
                  - f_88 * gi_357[k]
                  + f_89 * gi_359[k]
                  - f_90 * gi_361[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_28, gh_30, gh_37, gh_39, gh_41, gh_128, gh_133, \
                         gh_135, gh_142, gh_144, gh_146, gh_170, gh_175, gh_177, gh_184, \
                         gh_186, gh_188, gh_212, gh_217, gh_219, gh_226, gh_228, gh_230, \
                         gh_254, gh_259, gh_261, gh_268, gh_270, gh_272, gi_30, gi_35, gi_37, \
                         gi_44, gi_46, gi_48, gi_170, gi_175, gi_177, gi_184, gi_186, gi_188, \
                         gi_226, gi_231, gi_233, gi_240, gi_242, gi_244, gi_284, gi_291, \
                         gi_293, gi_302, gi_304, gi_306, gi_340, gi_347, gi_349, gi_358, \
                         gi_360, gi_362 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = f_91 * ab_x[k] * gh_23[k]
                  + f_92 * ab_x[k] * gh_28[k]
                  - f_93 * ab_x[k] * gh_30[k]
                  + f_91 * ab_x[k] * gh_37[k]
                  - f_93 * ab_x[k] * gh_39[k]
                  + f_94 * ab_x[k] * gh_41[k]
                  + f_95 * ab_x[k] * gh_128[k]
                  + f_96 * ab_x[k] * gh_133[k]
                  - f_97 * ab_x[k] * gh_135[k]
                  + f_95 * ab_x[k] * gh_142[k]
                  - f_97 * ab_x[k] * gh_144[k]
                  + f_98 * ab_x[k] * gh_146[k]
                  - f_99 * ab_x[k] * gh_170[k]
                  - f_100 * ab_x[k] * gh_175[k]
                  + f_101 * ab_x[k] * gh_177[k]
                  - f_99 * ab_x[k] * gh_184[k]
                  + f_101 * ab_x[k] * gh_186[k]
                  - f_102 * ab_x[k] * gh_188[k]
                  - f_103 * ab_y[k] * gh_212[k]
                  - f_95 * ab_y[k] * gh_217[k]
                  + f_104 * ab_y[k] * gh_219[k]
                  - f_103 * ab_y[k] * gh_226[k]
                  + f_104 * ab_y[k] * gh_228[k]
                  - f_105 * ab_y[k] * gh_230[k]
                  + f_93 * ab_y[k] * gh_254[k]
                  + f_106 * ab_y[k] * gh_259[k]
                  - f_107 * ab_y[k] * gh_261[k]
                  + f_93 * ab_y[k] * gh_268[k]
                  - f_107 * ab_y[k] * gh_270[k]
                  + f_108 * ab_y[k] * gh_272[k]
                  - f_91 * gi_30[k]
                  - f_92 * gi_35[k]
                  + f_93 * gi_37[k]
                  - f_91 * gi_44[k]
                  + f_93 * gi_46[k]
                  - f_94 * gi_48[k]
                  - f_95 * gi_170[k]
                  - f_96 * gi_175[k]
                  + f_97 * gi_177[k]
                  - f_95 * gi_184[k]
                  + f_97 * gi_186[k]
                  - f_98 * gi_188[k]
                  + f_99 * gi_226[k]
                  + f_100 * gi_231[k]
                  - f_101 * gi_233[k]
                  + f_99 * gi_240[k]
                  - f_101 * gi_242[k]
                  + f_102 * gi_244[k]
                  + f_103 * gi_284[k]
                  + f_95 * gi_291[k]
                  - f_104 * gi_293[k]
                  + f_103 * gi_302[k]
                  - f_104 * gi_304[k]
                  + f_105 * gi_306[k]
                  - f_93 * gi_340[k]
                  - f_106 * gi_347[k]
                  + f_107 * gi_349[k]
                  - f_93 * gi_358[k]
                  + f_107 * gi_360[k]
                  - f_108 * gi_362[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_26, gh_31, gh_33, gh_35, gh_126, gh_129, \
                         gh_131, gh_136, gh_138, gh_140, gh_168, gh_171, gh_173, gh_178, \
                         gh_180, gh_182, gh_210, gh_213, gh_215, gh_220, gh_222, gh_224, \
                         gh_252, gh_255, gh_257, gh_262, gh_264, gh_266, gi_28, gi_31, gi_33, \
                         gi_38, gi_40, gi_42, gi_168, gi_171, gi_173, gi_178, gi_180, gi_182, \
                         gi_224, gi_227, gi_229, gi_234, gi_236, gi_238, gi_281, gi_286, \
                         gi_288, gi_295, gi_297, gi_299, gi_337, gi_342, gi_344, gi_351, \
                         gi_353, gi_355 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_76 * ab_x[k] * gh_21[k]
                  + f_77 * ab_x[k] * gh_24[k]
                  - f_78 * ab_x[k] * gh_26[k]
                  + f_76 * ab_x[k] * gh_31[k]
                  - f_78 * ab_x[k] * gh_33[k]
                  + f_79 * ab_x[k] * gh_35[k]
                  + f_80 * ab_x[k] * gh_126[k]
                  + f_81 * ab_x[k] * gh_129[k]
                  - f_79 * ab_x[k] * gh_131[k]
                  + f_80 * ab_x[k] * gh_136[k]
                  - f_79 * ab_x[k] * gh_138[k]
                  + f_82 * ab_x[k] * gh_140[k]
                  - f_79 * ab_x[k] * gh_168[k]
                  - f_83 * ab_x[k] * gh_171[k]
                  + f_84 * ab_x[k] * gh_173[k]
                  - f_79 * ab_x[k] * gh_178[k]
                  + f_84 * ab_x[k] * gh_180[k]
                  - f_85 * ab_x[k] * gh_182[k]
                  - f_86 * ab_y[k] * gh_210[k]
                  - f_80 * ab_y[k] * gh_213[k]
                  + f_87 * ab_y[k] * gh_215[k]
                  - f_86 * ab_y[k] * gh_220[k]
                  + f_87 * ab_y[k] * gh_222[k]
                  - f_88 * ab_y[k] * gh_224[k]
                  + f_88 * ab_y[k] * gh_252[k]
                  + f_82 * ab_y[k] * gh_255[k]
                  - f_89 * ab_y[k] * gh_257[k]
                  + f_88 * ab_y[k] * gh_262[k]
                  - f_89 * ab_y[k] * gh_264[k]
                  + f_90 * ab_y[k] * gh_266[k]
                  - f_76 * gi_28[k]
                  - f_77 * gi_31[k]
                  + f_78 * gi_33[k]
                  - f_76 * gi_38[k]
                  + f_78 * gi_40[k]
                  - f_79 * gi_42[k]
                  - f_80 * gi_168[k]
                  - f_81 * gi_171[k]
                  + f_79 * gi_173[k]
                  - f_80 * gi_178[k]
                  + f_79 * gi_180[k]
                  - f_82 * gi_182[k]
                  + f_79 * gi_224[k]
                  + f_83 * gi_227[k]
                  - f_84 * gi_229[k]
                  + f_79 * gi_234[k]
                  - f_84 * gi_236[k]
                  + f_85 * gi_238[k]
                  + f_86 * gi_281[k]
                  + f_80 * gi_286[k]
                  - f_87 * gi_288[k]
                  + f_86 * gi_295[k]
                  - f_87 * gi_297[k]
                  + f_88 * gi_299[k]
                  - f_88 * gi_337[k]
                  - f_82 * gi_342[k]
                  + f_89 * gi_344[k]
                  - f_88 * gi_351[k]
                  + f_89 * gi_353[k]
                  - f_90 * gi_355[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_30, gh_37, gh_39, gh_128, gh_135, gh_142, \
                         gh_144, gh_170, gh_177, gh_184, gh_186, gh_212, gh_219, gh_226, \
                         gh_228, gh_254, gh_261, gh_268, gh_270, gi_30, gi_37, gi_44, gi_46, \
                         gi_170, gi_177, gi_184, gi_186, gi_226, gi_233, gi_240, gi_242, \
                         gi_284, gi_293, gi_302, gi_304, gi_340, gi_349, gi_358, \
                         gi_360 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = -f_109 * ab_x[k] * gh_23[k]
                  + f_67 * ab_x[k] * gh_30[k]
                  + f_109 * ab_x[k] * gh_37[k]
                  - f_67 * ab_x[k] * gh_39[k]
                  - f_73 * ab_x[k] * gh_128[k]
                  + f_69 * ab_x[k] * gh_135[k]
                  + f_73 * ab_x[k] * gh_142[k]
                  - f_69 * ab_x[k] * gh_144[k]
                  + f_110 * ab_x[k] * gh_170[k]
                  - f_71 * ab_x[k] * gh_177[k]
                  - f_110 * ab_x[k] * gh_184[k]
                  + f_71 * ab_x[k] * gh_186[k]
                  + f_111 * ab_y[k] * gh_212[k]
                  - f_73 * ab_y[k] * gh_219[k]
                  - f_111 * ab_y[k] * gh_226[k]
                  + f_73 * ab_y[k] * gh_228[k]
                  - f_70 * ab_y[k] * gh_254[k]
                  + f_74 * ab_y[k] * gh_261[k]
                  + f_70 * ab_y[k] * gh_268[k]
                  - f_74 * ab_y[k] * gh_270[k]
                  + f_109 * gi_30[k]
                  - f_67 * gi_37[k]
                  - f_109 * gi_44[k]
                  + f_67 * gi_46[k]
                  + f_73 * gi_170[k]
                  - f_69 * gi_177[k]
                  - f_73 * gi_184[k]
                  + f_69 * gi_186[k]
                  - f_110 * gi_226[k]
                  + f_71 * gi_233[k]
                  + f_110 * gi_240[k]
                  - f_71 * gi_242[k]
                  - f_111 * gi_284[k]
                  + f_73 * gi_293[k]
                  + f_111 * gi_302[k]
                  - f_73 * gi_304[k]
                  + f_70 * gi_340[k]
                  - f_74 * gi_349[k]
                  - f_70 * gi_358[k]
                  + f_74 * gi_360[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_26, gh_31, gh_33, gh_126, gh_129, \
                         gh_131, gh_136, gh_138, gh_168, gh_171, gh_173, gh_178, gh_180, \
                         gh_210, gh_213, gh_215, gh_220, gh_222, gh_252, gh_255, gh_257, \
                         gh_262, gh_264, gi_28, gi_31, gi_33, gi_38, gi_40, gi_168, gi_171, \
                         gi_173, gi_178, gi_180, gi_224, gi_227, gi_229, gi_234, gi_236, \
                         gi_281, gi_286, gi_288, gi_295, gi_297, gi_337, gi_342, gi_344, \
                         gi_351, gi_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -0.8203125 * ab_x[k] * gh_21[k]
                  + 1.640625 * ab_x[k] * gh_24[k]
                  + 6.5625 * ab_x[k] * gh_26[k]
                  + 2.4609375 * ab_x[k] * gh_31[k]
                  - 19.6875 * ab_x[k] * gh_33[k]
                  - 0.546875 * ab_x[k] * gh_126[k]
                  + 1.09375 * ab_x[k] * gh_129[k]
                  + 4.375 * ab_x[k] * gh_131[k]
                  + 1.640625 * ab_x[k] * gh_136[k]
                  - 13.125 * ab_x[k] * gh_138[k]
                  + 6.5625 * ab_x[k] * gh_168[k]
                  - 13.125 * ab_x[k] * gh_171[k]
                  - 52.5 * ab_x[k] * gh_173[k]
                  - 19.6875 * ab_x[k] * gh_178[k]
                  + 157.5 * ab_x[k] * gh_180[k]
                  + 0.2734375 * ab_y[k] * gh_210[k]
                  - 0.546875 * ab_y[k] * gh_213[k]
                  - 2.1875 * ab_y[k] * gh_215[k]
                  - 0.8203125 * ab_y[k] * gh_220[k]
                  + 6.5625 * ab_y[k] * gh_222[k]
                  - 2.1875 * ab_y[k] * gh_252[k]
                  + 4.375 * ab_y[k] * gh_255[k]
                  + 17.5 * ab_y[k] * gh_257[k]
                  + 6.5625 * ab_y[k] * gh_262[k]
                  - 52.5 * ab_y[k] * gh_264[k]
                  + 0.8203125 * gi_28[k]
                  - 1.640625 * gi_31[k]
                  - 6.5625 * gi_33[k]
                  - 2.4609375 * gi_38[k]
                  + 19.6875 * gi_40[k]
                  + 0.546875 * gi_168[k]
                  - 1.09375 * gi_171[k]
                  - 4.375 * gi_173[k]
                  - 1.640625 * gi_178[k]
                  + 13.125 * gi_180[k]
                  - 6.5625 * gi_224[k]
                  + 13.125 * gi_227[k]
                  + 52.5 * gi_229[k]
                  + 19.6875 * gi_234[k]
                  - 157.5 * gi_236[k]
                  - 0.2734375 * gi_281[k]
                  + 0.546875 * gi_286[k]
                  + 2.1875 * gi_288[k]
                  + 0.8203125 * gi_295[k]
                  - 6.5625 * gi_297[k]
                  + 2.1875 * gi_337[k]
                  - 4.375 * gi_342[k]
                  - 17.5 * gi_344[k]
                  - 6.5625 * gi_351[k]
                  + 52.5 * gi_353[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_28, gh_37, gh_128, gh_133, gh_142, gh_170, \
                         gh_175, gh_184, gh_212, gh_217, gh_226, gh_254, gh_259, gh_268, \
                         gi_30, gi_35, gi_44, gi_170, gi_175, gi_184, gi_226, gi_231, gi_240, \
                         gi_284, gi_291, gi_302, gi_340, gi_347, \
                         gi_358 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = f_112 * ab_x[k] * gh_23[k]
                  - f_113 * ab_x[k] * gh_28[k]
                  + f_112 * ab_x[k] * gh_37[k]
                  + f_114 * ab_x[k] * gh_128[k]
                  - f_51 * ab_x[k] * gh_133[k]
                  + f_114 * ab_x[k] * gh_142[k]
                  - f_115 * ab_x[k] * gh_170[k]
                  + f_116 * ab_x[k] * gh_175[k]
                  - f_115 * ab_x[k] * gh_184[k]
                  - f_117 * ab_y[k] * gh_212[k]
                  + f_118 * ab_y[k] * gh_217[k]
                  - f_117 * ab_y[k] * gh_226[k]
                  + f_52 * ab_y[k] * gh_254[k]
                  - f_119 * ab_y[k] * gh_259[k]
                  + f_52 * ab_y[k] * gh_268[k]
                  - f_112 * gi_30[k]
                  + f_113 * gi_35[k]
                  - f_112 * gi_44[k]
                  - f_114 * gi_170[k]
                  + f_51 * gi_175[k]
                  - f_114 * gi_184[k]
                  + f_115 * gi_226[k]
                  - f_116 * gi_231[k]
                  + f_115 * gi_240[k]
                  + f_117 * gi_284[k]
                  - f_118 * gi_291[k]
                  + f_117 * gi_302[k]
                  - f_52 * gi_340[k]
                  + f_119 * gi_347[k]
                  - f_52 * gi_358[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_31, gh_126, gh_129, gh_136, gh_168, \
                         gh_171, gh_178, gh_210, gh_213, gh_220, gh_252, gh_255, gh_262, \
                         gi_28, gi_31, gi_38, gi_168, gi_171, gi_178, gi_224, gi_227, gi_234, \
                         gi_281, gi_286, gi_295, gi_337, gi_342, \
                         gi_351 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_12 * ab_x[k] * gh_21[k]
                  - f_8 * ab_x[k] * gh_24[k]
                  + f_3 * ab_x[k] * gh_31[k]
                  + f_13 * ab_x[k] * gh_126[k]
                  - f_9 * ab_x[k] * gh_129[k]
                  + f_4 * ab_x[k] * gh_136[k]
                  - f_14 * ab_x[k] * gh_168[k]
                  + f_10 * ab_x[k] * gh_171[k]
                  - f_5 * ab_x[k] * gh_178[k]
                  - f_15 * ab_y[k] * gh_210[k]
                  + f_4 * ab_y[k] * gh_213[k]
                  - f_6 * ab_y[k] * gh_220[k]
                  + f_16 * ab_y[k] * gh_252[k]
                  - f_11 * ab_y[k] * gh_255[k]
                  + f_7 * ab_y[k] * gh_262[k]
                  - f_12 * gi_28[k]
                  + f_8 * gi_31[k]
                  - f_3 * gi_38[k]
                  - f_13 * gi_168[k]
                  + f_9 * gi_171[k]
                  - f_4 * gi_178[k]
                  + f_14 * gi_224[k]
                  - f_10 * gi_227[k]
                  + f_5 * gi_234[k]
                  + f_15 * gi_281[k]
                  - f_4 * gi_286[k]
                  + f_6 * gi_295[k]
                  - f_16 * gi_337[k]
                  + f_11 * gi_342[k]
                  - f_7 * gi_351[k];
    }

#pragma omp simd aligned(ab_x, gh_85, gh_90, gh_99, gh_232, gh_237, gh_246, gh_274, gh_279, \
                         gh_288, gi_113, gi_118, gi_127, gi_309, gi_314, gi_323, gi_365, \
                         gi_370, gi_379 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_17 * ab_x[k] * gh_85[k]
                  - f_18 * ab_x[k] * gh_90[k]
                  + f_20 * ab_x[k] * gh_99[k]
                  + f_17 * ab_x[k] * gh_232[k]
                  - f_18 * ab_x[k] * gh_237[k]
                  + f_20 * ab_x[k] * gh_246[k]
                  - f_18 * ab_x[k] * gh_274[k]
                  + f_19 * ab_x[k] * gh_279[k]
                  - f_21 * ab_x[k] * gh_288[k]
                  - f_17 * gi_113[k]
                  + f_18 * gi_118[k]
                  - f_20 * gi_127[k]
                  - f_17 * gi_309[k]
                  + f_18 * gi_314[k]
                  - f_20 * gi_323[k]
                  + f_18 * gi_365[k]
                  - f_19 * gi_370[k]
                  + f_21 * gi_379[k];
    }

#pragma omp simd aligned(ab_x, gh_88, gh_95, gh_235, gh_242, gh_277, gh_284, gi_116, gi_123, \
                         gi_312, gi_319, gi_368, gi_375 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_56 * ab_x[k] * gh_88[k]
                  - f_56 * ab_x[k] * gh_95[k]
                  + f_56 * ab_x[k] * gh_235[k]
                  - f_56 * ab_x[k] * gh_242[k]
                  - f_57 * ab_x[k] * gh_277[k]
                  + f_57 * ab_x[k] * gh_284[k]
                  - f_56 * gi_116[k]
                  + f_56 * gi_123[k]
                  - f_56 * gi_312[k]
                  + f_56 * gi_319[k]
                  + f_57 * gi_368[k]
                  - f_57 * gi_375[k];
    }

#pragma omp simd aligned(ab_x, gh_85, gh_90, gh_92, gh_99, gh_101, gh_232, gh_237, gh_239, \
                         gh_246, gh_248, gh_274, gh_279, gh_281, gh_288, gh_290, gi_113, \
                         gi_118, gi_120, gi_127, gi_129, gi_309, gi_314, gi_316, gi_323, \
                         gi_325, gi_365, gi_370, gi_372, gi_379, \
                         gi_381 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_67 * ab_x[k] * gh_85[k]
                  - f_69 * ab_x[k] * gh_90[k]
                  + f_71 * ab_x[k] * gh_92[k]
                  + f_73 * ab_x[k] * gh_99[k]
                  - f_74 * ab_x[k] * gh_101[k]
                  - f_67 * ab_x[k] * gh_232[k]
                  - f_69 * ab_x[k] * gh_237[k]
                  + f_71 * ab_x[k] * gh_239[k]
                  + f_73 * ab_x[k] * gh_246[k]
                  - f_74 * ab_x[k] * gh_248[k]
                  + f_68 * ab_x[k] * gh_274[k]
                  + f_70 * ab_x[k] * gh_279[k]
                  - f_72 * ab_x[k] * gh_281[k]
                  - f_69 * ab_x[k] * gh_288[k]
                  + f_75 * ab_x[k] * gh_290[k]
                  + f_67 * gi_113[k]
                  + f_69 * gi_118[k]
                  - f_71 * gi_120[k]
                  - f_73 * gi_127[k]
                  + f_74 * gi_129[k]
                  + f_67 * gi_309[k]
                  + f_69 * gi_314[k]
                  - f_71 * gi_316[k]
                  - f_73 * gi_323[k]
                  + f_74 * gi_325[k]
                  - f_68 * gi_365[k]
                  - f_70 * gi_370[k]
                  + f_72 * gi_372[k]
                  + f_69 * gi_379[k]
                  - f_75 * gi_381[k];
    }

#pragma omp simd aligned(ab_x, gh_88, gh_95, gh_97, gh_235, gh_242, gh_244, gh_277, gh_284, \
                         gh_286, gi_116, gi_123, gi_125, gi_312, gi_319, gi_321, gi_368, \
                         gi_375, gi_377 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = -26.25 * ab_x[k] * gh_88[k]
                  - 26.25 * ab_x[k] * gh_95[k]
                  + 52.5 * ab_x[k] * gh_97[k]
                  - 26.25 * ab_x[k] * gh_235[k]
                  - 26.25 * ab_x[k] * gh_242[k]
                  + 52.5 * ab_x[k] * gh_244[k]
                  + 52.5 * ab_x[k] * gh_277[k]
                  + 52.5 * ab_x[k] * gh_284[k]
                  - 105.0 * ab_x[k] * gh_286[k]
                  + 26.25 * gi_116[k]
                  + 26.25 * gi_123[k]
                  - 52.5 * gi_125[k]
                  + 26.25 * gi_312[k]
                  + 26.25 * gi_319[k]
                  - 52.5 * gi_321[k]
                  - 52.5 * gi_368[k]
                  - 52.5 * gi_375[k]
                  + 105.0 * gi_377[k];
    }

#pragma omp simd aligned(ab_x, gh_85, gh_90, gh_92, gh_99, gh_101, gh_103, gh_232, gh_237, \
                         gh_239, gh_246, gh_248, gh_250, gh_274, gh_279, gh_281, gh_288, \
                         gh_290, gh_292, gi_113, gi_118, gi_120, gi_127, gi_129, gi_131, \
                         gi_309, gi_314, gi_316, gi_323, gi_325, gi_327, gi_365, gi_370, \
                         gi_372, gi_379, gi_381, gi_383 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_120 * ab_x[k] * gh_85[k]
                  + f_121 * ab_x[k] * gh_90[k]
                  - f_122 * ab_x[k] * gh_92[k]
                  + f_120 * ab_x[k] * gh_99[k]
                  - f_122 * ab_x[k] * gh_101[k]
                  + f_123 * ab_x[k] * gh_103[k]
                  + f_120 * ab_x[k] * gh_232[k]
                  + f_121 * ab_x[k] * gh_237[k]
                  - f_122 * ab_x[k] * gh_239[k]
                  + f_120 * ab_x[k] * gh_246[k]
                  - f_122 * ab_x[k] * gh_248[k]
                  + f_123 * ab_x[k] * gh_250[k]
                  - f_121 * ab_x[k] * gh_274[k]
                  - f_124 * ab_x[k] * gh_279[k]
                  + f_125 * ab_x[k] * gh_281[k]
                  - f_121 * ab_x[k] * gh_288[k]
                  + f_125 * ab_x[k] * gh_290[k]
                  - f_126 * ab_x[k] * gh_292[k]
                  - f_120 * gi_113[k]
                  - f_121 * gi_118[k]
                  + f_122 * gi_120[k]
                  - f_120 * gi_127[k]
                  + f_122 * gi_129[k]
                  - f_123 * gi_131[k]
                  - f_120 * gi_309[k]
                  - f_121 * gi_314[k]
                  + f_122 * gi_316[k]
                  - f_120 * gi_323[k]
                  + f_122 * gi_325[k]
                  - f_123 * gi_327[k]
                  + f_121 * gi_365[k]
                  + f_124 * gi_370[k]
                  - f_125 * gi_372[k]
                  + f_121 * gi_379[k]
                  - f_125 * gi_381[k]
                  + f_126 * gi_383[k];
    }

#pragma omp simd aligned(ab_x, gh_86, gh_91, gh_93, gh_100, gh_102, gh_104, gh_233, gh_238, \
                         gh_240, gh_247, gh_249, gh_251, gh_275, gh_280, gh_282, gh_289, \
                         gh_291, gh_293, gi_114, gi_119, gi_121, gi_128, gi_130, gi_132, \
                         gi_310, gi_315, gi_317, gi_324, gi_326, gi_328, gi_366, gi_371, \
                         gi_373, gi_380, gi_382, gi_384 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_127 * ab_x[k] * gh_86[k]
                  + f_128 * ab_x[k] * gh_91[k]
                  - f_129 * ab_x[k] * gh_93[k]
                  + f_127 * ab_x[k] * gh_100[k]
                  - f_129 * ab_x[k] * gh_102[k]
                  + f_130 * ab_x[k] * gh_104[k]
                  + f_127 * ab_x[k] * gh_233[k]
                  + f_128 * ab_x[k] * gh_238[k]
                  - f_129 * ab_x[k] * gh_240[k]
                  + f_127 * ab_x[k] * gh_247[k]
                  - f_129 * ab_x[k] * gh_249[k]
                  + f_130 * ab_x[k] * gh_251[k]
                  - f_128 * ab_x[k] * gh_275[k]
                  - f_131 * ab_x[k] * gh_280[k]
                  + f_132 * ab_x[k] * gh_282[k]
                  - f_128 * ab_x[k] * gh_289[k]
                  + f_132 * ab_x[k] * gh_291[k]
                  - f_133 * ab_x[k] * gh_293[k]
                  - f_127 * gi_114[k]
                  - f_128 * gi_119[k]
                  + f_129 * gi_121[k]
                  - f_127 * gi_128[k]
                  + f_129 * gi_130[k]
                  - f_130 * gi_132[k]
                  - f_127 * gi_310[k]
                  - f_128 * gi_315[k]
                  + f_129 * gi_317[k]
                  - f_127 * gi_324[k]
                  + f_129 * gi_326[k]
                  - f_130 * gi_328[k]
                  + f_128 * gi_366[k]
                  + f_131 * gi_371[k]
                  - f_132 * gi_373[k]
                  + f_128 * gi_380[k]
                  - f_132 * gi_382[k]
                  + f_133 * gi_384[k];
    }

#pragma omp simd aligned(ab_x, gh_84, gh_87, gh_89, gh_94, gh_96, gh_98, gh_231, gh_234, \
                         gh_236, gh_241, gh_243, gh_245, gh_273, gh_276, gh_278, gh_283, \
                         gh_285, gh_287, gi_112, gi_115, gi_117, gi_122, gi_124, gi_126, \
                         gi_308, gi_311, gi_313, gi_318, gi_320, gi_322, gi_364, gi_367, \
                         gi_369, gi_374, gi_376, gi_378 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = f_120 * ab_x[k] * gh_84[k]
                  + f_121 * ab_x[k] * gh_87[k]
                  - f_122 * ab_x[k] * gh_89[k]
                  + f_120 * ab_x[k] * gh_94[k]
                  - f_122 * ab_x[k] * gh_96[k]
                  + f_123 * ab_x[k] * gh_98[k]
                  + f_120 * ab_x[k] * gh_231[k]
                  + f_121 * ab_x[k] * gh_234[k]
                  - f_122 * ab_x[k] * gh_236[k]
                  + f_120 * ab_x[k] * gh_241[k]
                  - f_122 * ab_x[k] * gh_243[k]
                  + f_123 * ab_x[k] * gh_245[k]
                  - f_121 * ab_x[k] * gh_273[k]
                  - f_124 * ab_x[k] * gh_276[k]
                  + f_125 * ab_x[k] * gh_278[k]
                  - f_121 * ab_x[k] * gh_283[k]
                  + f_125 * ab_x[k] * gh_285[k]
                  - f_126 * ab_x[k] * gh_287[k]
                  - f_120 * gi_112[k]
                  - f_121 * gi_115[k]
                  + f_122 * gi_117[k]
                  - f_120 * gi_122[k]
                  + f_122 * gi_124[k]
                  - f_123 * gi_126[k]
                  - f_120 * gi_308[k]
                  - f_121 * gi_311[k]
                  + f_122 * gi_313[k]
                  - f_120 * gi_318[k]
                  + f_122 * gi_320[k]
                  - f_123 * gi_322[k]
                  + f_121 * gi_364[k]
                  + f_124 * gi_367[k]
                  - f_125 * gi_369[k]
                  + f_121 * gi_374[k]
                  - f_125 * gi_376[k]
                  + f_126 * gi_378[k];
    }

#pragma omp simd aligned(ab_x, gh_86, gh_93, gh_100, gh_102, gh_233, gh_240, gh_247, gh_249, \
                         gh_275, gh_282, gh_289, gh_291, gi_114, gi_121, gi_128, gi_130, \
                         gi_310, gi_317, gi_324, gi_326, gi_366, gi_373, gi_380, \
                         gi_382 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -13.125 * ab_x[k] * gh_86[k]
                  + 26.25 * ab_x[k] * gh_93[k]
                  + 13.125 * ab_x[k] * gh_100[k]
                  - 26.25 * ab_x[k] * gh_102[k]
                  - 13.125 * ab_x[k] * gh_233[k]
                  + 26.25 * ab_x[k] * gh_240[k]
                  + 13.125 * ab_x[k] * gh_247[k]
                  - 26.25 * ab_x[k] * gh_249[k]
                  + 26.25 * ab_x[k] * gh_275[k]
                  - 52.5 * ab_x[k] * gh_282[k]
                  - 26.25 * ab_x[k] * gh_289[k]
                  + 52.5 * ab_x[k] * gh_291[k]
                  + 13.125 * gi_114[k]
                  - 26.25 * gi_121[k]
                  - 13.125 * gi_128[k]
                  + 26.25 * gi_130[k]
                  + 13.125 * gi_310[k]
                  - 26.25 * gi_317[k]
                  - 13.125 * gi_324[k]
                  + 26.25 * gi_326[k]
                  - 26.25 * gi_366[k]
                  + 52.5 * gi_373[k]
                  + 26.25 * gi_380[k]
                  - 52.5 * gi_382[k];
    }

#pragma omp simd aligned(ab_x, gh_84, gh_87, gh_89, gh_94, gh_96, gh_231, gh_234, gh_236, \
                         gh_241, gh_243, gh_273, gh_276, gh_278, gh_283, gh_285, gi_112, \
                         gi_115, gi_117, gi_122, gi_124, gi_308, gi_311, gi_313, gi_318, \
                         gi_320, gi_364, gi_367, gi_369, gi_374, \
                         gi_376 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_73 * ab_x[k] * gh_84[k]
                  + f_69 * ab_x[k] * gh_87[k]
                  + f_74 * ab_x[k] * gh_89[k]
                  + f_67 * ab_x[k] * gh_94[k]
                  - f_71 * ab_x[k] * gh_96[k]
                  - f_73 * ab_x[k] * gh_231[k]
                  + f_69 * ab_x[k] * gh_234[k]
                  + f_74 * ab_x[k] * gh_236[k]
                  + f_67 * ab_x[k] * gh_241[k]
                  - f_71 * ab_x[k] * gh_243[k]
                  + f_69 * ab_x[k] * gh_273[k]
                  - f_70 * ab_x[k] * gh_276[k]
                  - f_75 * ab_x[k] * gh_278[k]
                  - f_68 * ab_x[k] * gh_283[k]
                  + f_72 * ab_x[k] * gh_285[k]
                  + f_73 * gi_112[k]
                  - f_69 * gi_115[k]
                  - f_74 * gi_117[k]
                  - f_67 * gi_122[k]
                  + f_71 * gi_124[k]
                  + f_73 * gi_308[k]
                  - f_69 * gi_311[k]
                  - f_74 * gi_313[k]
                  - f_67 * gi_318[k]
                  + f_71 * gi_320[k]
                  - f_69 * gi_364[k]
                  + f_70 * gi_367[k]
                  + f_75 * gi_369[k]
                  + f_68 * gi_374[k]
                  - f_72 * gi_376[k];
    }

#pragma omp simd aligned(ab_x, gh_86, gh_91, gh_100, gh_233, gh_238, gh_247, gh_275, gh_280, \
                         gh_289, gi_114, gi_119, gi_128, gi_310, gi_315, gi_324, gi_366, \
                         gi_371, gi_380 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_134 * ab_x[k] * gh_86[k]
                  - f_135 * ab_x[k] * gh_91[k]
                  + f_134 * ab_x[k] * gh_100[k]
                  + f_134 * ab_x[k] * gh_233[k]
                  - f_135 * ab_x[k] * gh_238[k]
                  + f_134 * ab_x[k] * gh_247[k]
                  - f_66 * ab_x[k] * gh_275[k]
                  + f_136 * ab_x[k] * gh_280[k]
                  - f_66 * ab_x[k] * gh_289[k]
                  - f_134 * gi_114[k]
                  + f_135 * gi_119[k]
                  - f_134 * gi_128[k]
                  - f_134 * gi_310[k]
                  + f_135 * gi_315[k]
                  - f_134 * gi_324[k]
                  + f_66 * gi_366[k]
                  - f_136 * gi_371[k]
                  + f_66 * gi_380[k];
    }

#pragma omp simd aligned(ab_x, gh_84, gh_87, gh_94, gh_231, gh_234, gh_241, gh_273, gh_276, \
                         gh_283, gi_112, gi_115, gi_122, gi_308, gi_311, gi_318, gi_364, \
                         gi_367, gi_374 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_20 * ab_x[k] * gh_84[k]
                  - f_18 * ab_x[k] * gh_87[k]
                  + f_17 * ab_x[k] * gh_94[k]
                  + f_20 * ab_x[k] * gh_231[k]
                  - f_18 * ab_x[k] * gh_234[k]
                  + f_17 * ab_x[k] * gh_241[k]
                  - f_21 * ab_x[k] * gh_273[k]
                  + f_19 * ab_x[k] * gh_276[k]
                  - f_18 * ab_x[k] * gh_283[k]
                  - f_20 * gi_112[k]
                  + f_18 * gi_115[k]
                  - f_17 * gi_122[k]
                  - f_20 * gi_308[k]
                  + f_18 * gi_311[k]
                  - f_17 * gi_318[k]
                  + f_21 * gi_364[k]
                  - f_19 * gi_367[k]
                  + f_18 * gi_374[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_22, gh_27, gh_36, gh_127, gh_132, gh_141, gh_169, \
                         gh_174, gh_183, gh_211, gh_216, gh_225, gh_253, gh_258, gh_267, \
                         gh_295, gh_300, gh_309, gi_29, gi_34, gi_43, gi_169, gi_174, gi_183, \
                         gi_225, gi_230, gi_239, gi_283, gi_290, gi_301, gi_339, gi_346, \
                         gi_357, gi_395, gi_402, gi_413 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_22 * ab_x[k] * gh_22[k]
                  + f_23 * ab_x[k] * gh_27[k]
                  - f_29 * ab_x[k] * gh_36[k]
                  - f_23 * ab_x[k] * gh_127[k]
                  + f_26 * ab_x[k] * gh_132[k]
                  - f_30 * ab_x[k] * gh_141[k]
                  + f_24 * ab_x[k] * gh_169[k]
                  - f_27 * ab_x[k] * gh_174[k]
                  + f_31 * ab_x[k] * gh_183[k]
                  - f_22 * ab_y[k] * gh_211[k]
                  + f_23 * ab_y[k] * gh_216[k]
                  - f_29 * ab_y[k] * gh_225[k]
                  + f_24 * ab_y[k] * gh_253[k]
                  - f_27 * ab_y[k] * gh_258[k]
                  + f_31 * ab_y[k] * gh_267[k]
                  - f_25 * ab_y[k] * gh_295[k]
                  + f_28 * ab_y[k] * gh_300[k]
                  - f_32 * ab_y[k] * gh_309[k]
                  + f_22 * gi_29[k]
                  - f_23 * gi_34[k]
                  + f_29 * gi_43[k]
                  + f_23 * gi_169[k]
                  - f_26 * gi_174[k]
                  + f_30 * gi_183[k]
                  - f_24 * gi_225[k]
                  + f_27 * gi_230[k]
                  - f_31 * gi_239[k]
                  + f_22 * gi_283[k]
                  - f_23 * gi_290[k]
                  + f_29 * gi_301[k]
                  - f_24 * gi_339[k]
                  + f_27 * gi_346[k]
                  - f_31 * gi_357[k]
                  + f_25 * gi_395[k]
                  - f_28 * gi_402[k]
                  + f_32 * gi_413[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_25, gh_32, gh_130, gh_137, gh_172, gh_179, gh_214, \
                         gh_221, gh_256, gh_263, gh_298, gh_305, gi_32, gi_39, gi_172, gi_179, \
                         gi_228, gi_235, gi_287, gi_296, gi_343, gi_352, gi_399, \
                         gi_408 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_58 * ab_x[k] * gh_25[k]
                  + f_58 * ab_x[k] * gh_32[k]
                  - f_59 * ab_x[k] * gh_130[k]
                  + f_59 * ab_x[k] * gh_137[k]
                  + f_60 * ab_x[k] * gh_172[k]
                  - f_60 * ab_x[k] * gh_179[k]
                  - f_58 * ab_y[k] * gh_214[k]
                  + f_58 * ab_y[k] * gh_221[k]
                  + f_60 * ab_y[k] * gh_256[k]
                  - f_60 * ab_y[k] * gh_263[k]
                  - f_61 * ab_y[k] * gh_298[k]
                  + f_61 * ab_y[k] * gh_305[k]
                  + f_58 * gi_32[k]
                  - f_58 * gi_39[k]
                  + f_59 * gi_172[k]
                  - f_59 * gi_179[k]
                  - f_60 * gi_228[k]
                  + f_60 * gi_235[k]
                  + f_58 * gi_287[k]
                  - f_58 * gi_296[k]
                  - f_60 * gi_343[k]
                  + f_60 * gi_352[k]
                  + f_61 * gi_399[k]
                  - f_61 * gi_408[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_22, gh_27, gh_29, gh_36, gh_38, gh_127, gh_132, \
                         gh_134, gh_141, gh_143, gh_169, gh_174, gh_176, gh_183, gh_185, \
                         gh_211, gh_216, gh_218, gh_225, gh_227, gh_253, gh_258, gh_260, \
                         gh_267, gh_269, gh_295, gh_300, gh_302, gh_309, gh_311, gi_29, gi_34, \
                         gi_36, gi_43, gi_45, gi_169, gi_174, gi_176, gi_183, gi_185, gi_225, \
                         gi_230, gi_232, gi_239, gi_241, gi_283, gi_290, gi_292, gi_301, \
                         gi_303, gi_339, gi_346, gi_348, gi_357, gi_359, gi_395, gi_402, \
                         gi_404, gi_413, gi_415 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_76 * ab_x[k] * gh_22[k]
                  + f_80 * ab_x[k] * gh_27[k]
                  - f_79 * ab_x[k] * gh_29[k]
                  - f_86 * ab_x[k] * gh_36[k]
                  + f_88 * ab_x[k] * gh_38[k]
                  + f_77 * ab_x[k] * gh_127[k]
                  + f_81 * ab_x[k] * gh_132[k]
                  - f_83 * ab_x[k] * gh_134[k]
                  - f_80 * ab_x[k] * gh_141[k]
                  + f_82 * ab_x[k] * gh_143[k]
                  - f_78 * ab_x[k] * gh_169[k]
                  - f_79 * ab_x[k] * gh_174[k]
                  + f_84 * ab_x[k] * gh_176[k]
                  + f_87 * ab_x[k] * gh_183[k]
                  - f_89 * ab_x[k] * gh_185[k]
                  + f_76 * ab_y[k] * gh_211[k]
                  + f_80 * ab_y[k] * gh_216[k]
                  - f_79 * ab_y[k] * gh_218[k]
                  - f_86 * ab_y[k] * gh_225[k]
                  + f_88 * ab_y[k] * gh_227[k]
                  - f_78 * ab_y[k] * gh_253[k]
                  - f_79 * ab_y[k] * gh_258[k]
                  + f_84 * ab_y[k] * gh_260[k]
                  + f_87 * ab_y[k] * gh_267[k]
                  - f_89 * ab_y[k] * gh_269[k]
                  + f_79 * ab_y[k] * gh_295[k]
                  + f_82 * ab_y[k] * gh_300[k]
                  - f_85 * ab_y[k] * gh_302[k]
                  - f_88 * ab_y[k] * gh_309[k]
                  + f_90 * ab_y[k] * gh_311[k]
                  - f_76 * gi_29[k]
                  - f_80 * gi_34[k]
                  + f_79 * gi_36[k]
                  + f_86 * gi_43[k]
                  - f_88 * gi_45[k]
                  - f_77 * gi_169[k]
                  - f_81 * gi_174[k]
                  + f_83 * gi_176[k]
                  + f_80 * gi_183[k]
                  - f_82 * gi_185[k]
                  + f_78 * gi_225[k]
                  + f_79 * gi_230[k]
                  - f_84 * gi_232[k]
                  - f_87 * gi_239[k]
                  + f_89 * gi_241[k]
                  - f_76 * gi_283[k]
                  - f_80 * gi_290[k]
                  + f_79 * gi_292[k]
                  + f_86 * gi_301[k]
                  - f_88 * gi_303[k]
                  + f_78 * gi_339[k]
                  + f_79 * gi_346[k]
                  - f_84 * gi_348[k]
                  - f_87 * gi_357[k]
                  + f_89 * gi_359[k]
                  - f_79 * gi_395[k]
                  - f_82 * gi_402[k]
                  + f_85 * gi_404[k]
                  + f_88 * gi_413[k]
                  - f_90 * gi_415[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_25, gh_32, gh_34, gh_130, gh_137, gh_139, gh_172, \
                         gh_179, gh_181, gh_214, gh_221, gh_223, gh_256, gh_263, gh_265, \
                         gh_298, gh_305, gh_307, gi_32, gi_39, gi_41, gi_172, gi_179, gi_181, \
                         gi_228, gi_235, gi_237, gi_287, gi_296, gi_298, gi_343, gi_352, \
                         gi_354, gi_399, gi_408, gi_410 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = f_120 * ab_x[k] * gh_25[k]
                  + f_120 * ab_x[k] * gh_32[k]
                  - f_121 * ab_x[k] * gh_34[k]
                  + f_121 * ab_x[k] * gh_130[k]
                  + f_121 * ab_x[k] * gh_137[k]
                  - f_124 * ab_x[k] * gh_139[k]
                  - f_122 * ab_x[k] * gh_172[k]
                  - f_122 * ab_x[k] * gh_179[k]
                  + f_125 * ab_x[k] * gh_181[k]
                  + f_120 * ab_y[k] * gh_214[k]
                  + f_120 * ab_y[k] * gh_221[k]
                  - f_121 * ab_y[k] * gh_223[k]
                  - f_122 * ab_y[k] * gh_256[k]
                  - f_122 * ab_y[k] * gh_263[k]
                  + f_125 * ab_y[k] * gh_265[k]
                  + f_123 * ab_y[k] * gh_298[k]
                  + f_123 * ab_y[k] * gh_305[k]
                  - f_126 * ab_y[k] * gh_307[k]
                  - f_120 * gi_32[k]
                  - f_120 * gi_39[k]
                  + f_121 * gi_41[k]
                  - f_121 * gi_172[k]
                  - f_121 * gi_179[k]
                  + f_124 * gi_181[k]
                  + f_122 * gi_228[k]
                  + f_122 * gi_235[k]
                  - f_125 * gi_237[k]
                  - f_120 * gi_287[k]
                  - f_120 * gi_296[k]
                  + f_121 * gi_298[k]
                  + f_122 * gi_343[k]
                  + f_122 * gi_352[k]
                  - f_125 * gi_354[k]
                  - f_123 * gi_399[k]
                  - f_123 * gi_408[k]
                  + f_126 * gi_410[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_22, gh_27, gh_29, gh_36, gh_38, gh_40, gh_127, gh_132, \
                         gh_134, gh_141, gh_143, gh_145, gh_169, gh_174, gh_176, gh_183, \
                         gh_185, gh_187, gh_211, gh_216, gh_218, gh_225, gh_227, gh_229, \
                         gh_253, gh_258, gh_260, gh_267, gh_269, gh_271, gh_295, gh_300, \
                         gh_302, gh_309, gh_311, gh_313, gi_29, gi_34, gi_36, gi_43, gi_45, \
                         gi_47, gi_169, gi_174, gi_176, gi_183, gi_185, gi_187, gi_225, \
                         gi_230, gi_232, gi_239, gi_241, gi_243, gi_283, gi_290, gi_292, \
                         gi_301, gi_303, gi_305, gi_339, gi_346, gi_348, gi_357, gi_359, \
                         gi_361, gi_395, gi_402, gi_404, gi_413, gi_415, \
                         gi_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -0.234375 * ab_x[k] * gh_22[k]
                  - 0.46875 * ab_x[k] * gh_27[k]
                  + 2.8125 * ab_x[k] * gh_29[k]
                  - 0.234375 * ab_x[k] * gh_36[k]
                  + 2.8125 * ab_x[k] * gh_38[k]
                  - 1.875 * ab_x[k] * gh_40[k]
                  - 0.46875 * ab_x[k] * gh_127[k]
                  - 0.9375 * ab_x[k] * gh_132[k]
                  + 5.625 * ab_x[k] * gh_134[k]
                  - 0.46875 * ab_x[k] * gh_141[k]
                  + 5.625 * ab_x[k] * gh_143[k]
                  - 3.75 * ab_x[k] * gh_145[k]
                  + 2.8125 * ab_x[k] * gh_169[k]
                  + 5.625 * ab_x[k] * gh_174[k]
                  - 33.75 * ab_x[k] * gh_176[k]
                  + 2.8125 * ab_x[k] * gh_183[k]
                  - 33.75 * ab_x[k] * gh_185[k]
                  + 22.5 * ab_x[k] * gh_187[k]
                  - 0.234375 * ab_y[k] * gh_211[k]
                  - 0.46875 * ab_y[k] * gh_216[k]
                  + 2.8125 * ab_y[k] * gh_218[k]
                  - 0.234375 * ab_y[k] * gh_225[k]
                  + 2.8125 * ab_y[k] * gh_227[k]
                  - 1.875 * ab_y[k] * gh_229[k]
                  + 2.8125 * ab_y[k] * gh_253[k]
                  + 5.625 * ab_y[k] * gh_258[k]
                  - 33.75 * ab_y[k] * gh_260[k]
                  + 2.8125 * ab_y[k] * gh_267[k]
                  - 33.75 * ab_y[k] * gh_269[k]
                  + 22.5 * ab_y[k] * gh_271[k]
                  - 1.875 * ab_y[k] * gh_295[k]
                  - 3.75 * ab_y[k] * gh_300[k]
                  + 22.5 * ab_y[k] * gh_302[k]
                  - 1.875 * ab_y[k] * gh_309[k]
                  + 22.5 * ab_y[k] * gh_311[k]
                  - 15.0 * ab_y[k] * gh_313[k]
                  + 0.234375 * gi_29[k]
                  + 0.46875 * gi_34[k]
                  - 2.8125 * gi_36[k]
                  + 0.234375 * gi_43[k]
                  - 2.8125 * gi_45[k]
                  + 1.875 * gi_47[k]
                  + 0.46875 * gi_169[k]
                  + 0.9375 * gi_174[k]
                  - 5.625 * gi_176[k]
                  + 0.46875 * gi_183[k]
                  - 5.625 * gi_185[k]
                  + 3.75 * gi_187[k]
                  - 2.8125 * gi_225[k]
                  - 5.625 * gi_230[k]
                  + 33.75 * gi_232[k]
                  - 2.8125 * gi_239[k]
                  + 33.75 * gi_241[k]
                  - 22.5 * gi_243[k]
                  + 0.234375 * gi_283[k]
                  + 0.46875 * gi_290[k]
                  - 2.8125 * gi_292[k]
                  + 0.234375 * gi_301[k]
                  - 2.8125 * gi_303[k]
                  + 1.875 * gi_305[k]
                  - 2.8125 * gi_339[k]
                  - 5.625 * gi_346[k]
                  + 33.75 * gi_348[k]
                  - 2.8125 * gi_357[k]
                  + 33.75 * gi_359[k]
                  - 22.5 * gi_361[k]
                  + 1.875 * gi_395[k]
                  + 3.75 * gi_402[k]
                  - 22.5 * gi_404[k]
                  + 1.875 * gi_413[k]
                  - 22.5 * gi_415[k]
                  + 15.0 * gi_417[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_28, gh_30, gh_37, gh_39, gh_41, gh_128, gh_133, \
                         gh_135, gh_142, gh_144, gh_146, gh_170, gh_175, gh_177, gh_184, \
                         gh_186, gh_188, gh_212, gh_217, gh_219, gh_226, gh_228, gh_230, \
                         gh_254, gh_259, gh_261, gh_268, gh_270, gh_272, gh_296, gh_301, \
                         gh_303, gh_310, gh_312, gh_314, gi_30, gi_35, gi_37, gi_44, gi_46, \
                         gi_48, gi_170, gi_175, gi_177, gi_184, gi_186, gi_188, gi_226, \
                         gi_231, gi_233, gi_240, gi_242, gi_244, gi_284, gi_291, gi_293, \
                         gi_302, gi_304, gi_306, gi_340, gi_347, gi_349, gi_358, gi_360, \
                         gi_362, gi_396, gi_403, gi_405, gi_414, gi_416, \
                         gi_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_137 * ab_x[k] * gh_23[k]
                  - f_138 * ab_x[k] * gh_28[k]
                  + f_139 * ab_x[k] * gh_30[k]
                  - f_137 * ab_x[k] * gh_37[k]
                  + f_139 * ab_x[k] * gh_39[k]
                  - f_140 * ab_x[k] * gh_41[k]
                  - f_138 * ab_x[k] * gh_128[k]
                  - f_141 * ab_x[k] * gh_133[k]
                  + f_142 * ab_x[k] * gh_135[k]
                  - f_138 * ab_x[k] * gh_142[k]
                  + f_142 * ab_x[k] * gh_144[k]
                  - f_143 * ab_x[k] * gh_146[k]
                  + f_144 * ab_x[k] * gh_170[k]
                  + f_145 * ab_x[k] * gh_175[k]
                  - f_146 * ab_x[k] * gh_177[k]
                  + f_144 * ab_x[k] * gh_184[k]
                  - f_146 * ab_x[k] * gh_186[k]
                  + f_147 * ab_x[k] * gh_188[k]
                  - f_137 * ab_y[k] * gh_212[k]
                  - f_138 * ab_y[k] * gh_217[k]
                  + f_139 * ab_y[k] * gh_219[k]
                  - f_137 * ab_y[k] * gh_226[k]
                  + f_139 * ab_y[k] * gh_228[k]
                  - f_140 * ab_y[k] * gh_230[k]
                  + f_144 * ab_y[k] * gh_254[k]
                  + f_145 * ab_y[k] * gh_259[k]
                  - f_146 * ab_y[k] * gh_261[k]
                  + f_144 * ab_y[k] * gh_268[k]
                  - f_146 * ab_y[k] * gh_270[k]
                  + f_147 * ab_y[k] * gh_272[k]
                  - f_148 * ab_y[k] * gh_296[k]
                  - f_149 * ab_y[k] * gh_301[k]
                  + f_150 * ab_y[k] * gh_303[k]
                  - f_148 * ab_y[k] * gh_310[k]
                  + f_150 * ab_y[k] * gh_312[k]
                  - f_151 * ab_y[k] * gh_314[k]
                  + f_137 * gi_30[k]
                  + f_138 * gi_35[k]
                  - f_139 * gi_37[k]
                  + f_137 * gi_44[k]
                  - f_139 * gi_46[k]
                  + f_140 * gi_48[k]
                  + f_138 * gi_170[k]
                  + f_141 * gi_175[k]
                  - f_142 * gi_177[k]
                  + f_138 * gi_184[k]
                  - f_142 * gi_186[k]
                  + f_143 * gi_188[k]
                  - f_144 * gi_226[k]
                  - f_145 * gi_231[k]
                  + f_146 * gi_233[k]
                  - f_144 * gi_240[k]
                  + f_146 * gi_242[k]
                  - f_147 * gi_244[k]
                  + f_137 * gi_284[k]
                  + f_138 * gi_291[k]
                  - f_139 * gi_293[k]
                  + f_137 * gi_302[k]
                  - f_139 * gi_304[k]
                  + f_140 * gi_306[k]
                  - f_144 * gi_340[k]
                  - f_145 * gi_347[k]
                  + f_146 * gi_349[k]
                  - f_144 * gi_358[k]
                  + f_146 * gi_360[k]
                  - f_147 * gi_362[k]
                  + f_148 * gi_396[k]
                  + f_149 * gi_403[k]
                  - f_150 * gi_405[k]
                  + f_148 * gi_414[k]
                  - f_150 * gi_416[k]
                  + f_151 * gi_418[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_26, gh_31, gh_33, gh_35, gh_126, gh_129, \
                         gh_131, gh_136, gh_138, gh_140, gh_168, gh_171, gh_173, gh_178, \
                         gh_180, gh_182, gh_210, gh_213, gh_215, gh_220, gh_222, gh_224, \
                         gh_252, gh_255, gh_257, gh_262, gh_264, gh_266, gh_294, gh_297, \
                         gh_299, gh_304, gh_306, gh_308, gi_28, gi_31, gi_33, gi_38, gi_40, \
                         gi_42, gi_168, gi_171, gi_173, gi_178, gi_180, gi_182, gi_224, \
                         gi_227, gi_229, gi_234, gi_236, gi_238, gi_281, gi_286, gi_288, \
                         gi_295, gi_297, gi_299, gi_337, gi_342, gi_344, gi_351, gi_353, \
                         gi_355, gi_393, gi_398, gi_400, gi_407, gi_409, \
                         gi_411 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -0.234375 * ab_x[k] * gh_21[k]
                  - 0.46875 * ab_x[k] * gh_24[k]
                  + 2.8125 * ab_x[k] * gh_26[k]
                  - 0.234375 * ab_x[k] * gh_31[k]
                  + 2.8125 * ab_x[k] * gh_33[k]
                  - 1.875 * ab_x[k] * gh_35[k]
                  - 0.46875 * ab_x[k] * gh_126[k]
                  - 0.9375 * ab_x[k] * gh_129[k]
                  + 5.625 * ab_x[k] * gh_131[k]
                  - 0.46875 * ab_x[k] * gh_136[k]
                  + 5.625 * ab_x[k] * gh_138[k]
                  - 3.75 * ab_x[k] * gh_140[k]
                  + 2.8125 * ab_x[k] * gh_168[k]
                  + 5.625 * ab_x[k] * gh_171[k]
                  - 33.75 * ab_x[k] * gh_173[k]
                  + 2.8125 * ab_x[k] * gh_178[k]
                  - 33.75 * ab_x[k] * gh_180[k]
                  + 22.5 * ab_x[k] * gh_182[k]
                  - 0.234375 * ab_y[k] * gh_210[k]
                  - 0.46875 * ab_y[k] * gh_213[k]
                  + 2.8125 * ab_y[k] * gh_215[k]
                  - 0.234375 * ab_y[k] * gh_220[k]
                  + 2.8125 * ab_y[k] * gh_222[k]
                  - 1.875 * ab_y[k] * gh_224[k]
                  + 2.8125 * ab_y[k] * gh_252[k]
                  + 5.625 * ab_y[k] * gh_255[k]
                  - 33.75 * ab_y[k] * gh_257[k]
                  + 2.8125 * ab_y[k] * gh_262[k]
                  - 33.75 * ab_y[k] * gh_264[k]
                  + 22.5 * ab_y[k] * gh_266[k]
                  - 1.875 * ab_y[k] * gh_294[k]
                  - 3.75 * ab_y[k] * gh_297[k]
                  + 22.5 * ab_y[k] * gh_299[k]
                  - 1.875 * ab_y[k] * gh_304[k]
                  + 22.5 * ab_y[k] * gh_306[k]
                  - 15.0 * ab_y[k] * gh_308[k]
                  + 0.234375 * gi_28[k]
                  + 0.46875 * gi_31[k]
                  - 2.8125 * gi_33[k]
                  + 0.234375 * gi_38[k]
                  - 2.8125 * gi_40[k]
                  + 1.875 * gi_42[k]
                  + 0.46875 * gi_168[k]
                  + 0.9375 * gi_171[k]
                  - 5.625 * gi_173[k]
                  + 0.46875 * gi_178[k]
                  - 5.625 * gi_180[k]
                  + 3.75 * gi_182[k]
                  - 2.8125 * gi_224[k]
                  - 5.625 * gi_227[k]
                  + 33.75 * gi_229[k]
                  - 2.8125 * gi_234[k]
                  + 33.75 * gi_236[k]
                  - 22.5 * gi_238[k]
                  + 0.234375 * gi_281[k]
                  + 0.46875 * gi_286[k]
                  - 2.8125 * gi_288[k]
                  + 0.234375 * gi_295[k]
                  - 2.8125 * gi_297[k]
                  + 1.875 * gi_299[k]
                  - 2.8125 * gi_337[k]
                  - 5.625 * gi_342[k]
                  + 33.75 * gi_344[k]
                  - 2.8125 * gi_351[k]
                  + 33.75 * gi_353[k]
                  - 22.5 * gi_355[k]
                  + 1.875 * gi_393[k]
                  + 3.75 * gi_398[k]
                  - 22.5 * gi_400[k]
                  + 1.875 * gi_407[k]
                  - 22.5 * gi_409[k]
                  + 15.0 * gi_411[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_30, gh_37, gh_39, gh_128, gh_135, gh_142, \
                         gh_144, gh_170, gh_177, gh_184, gh_186, gh_212, gh_219, gh_226, \
                         gh_228, gh_254, gh_261, gh_268, gh_270, gh_296, gh_303, gh_310, \
                         gh_312, gi_30, gi_37, gi_44, gi_46, gi_170, gi_177, gi_184, gi_186, \
                         gi_226, gi_233, gi_240, gi_242, gi_284, gi_293, gi_302, gi_304, \
                         gi_340, gi_349, gi_358, gi_360, gi_396, gi_405, gi_414, \
                         gi_416 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_152 * ab_x[k] * gh_23[k]
                  - f_120 * ab_x[k] * gh_30[k]
                  - f_152 * ab_x[k] * gh_37[k]
                  + f_120 * ab_x[k] * gh_39[k]
                  + f_120 * ab_x[k] * gh_128[k]
                  - f_121 * ab_x[k] * gh_135[k]
                  - f_120 * ab_x[k] * gh_142[k]
                  + f_121 * ab_x[k] * gh_144[k]
                  - f_153 * ab_x[k] * gh_170[k]
                  + f_122 * ab_x[k] * gh_177[k]
                  + f_153 * ab_x[k] * gh_184[k]
                  - f_122 * ab_x[k] * gh_186[k]
                  + f_152 * ab_y[k] * gh_212[k]
                  - f_120 * ab_y[k] * gh_219[k]
                  - f_152 * ab_y[k] * gh_226[k]
                  + f_120 * ab_y[k] * gh_228[k]
                  - f_153 * ab_y[k] * gh_254[k]
                  + f_122 * ab_y[k] * gh_261[k]
                  + f_153 * ab_y[k] * gh_268[k]
                  - f_122 * ab_y[k] * gh_270[k]
                  + f_124 * ab_y[k] * gh_296[k]
                  - f_123 * ab_y[k] * gh_303[k]
                  - f_124 * ab_y[k] * gh_310[k]
                  + f_123 * ab_y[k] * gh_312[k]
                  - f_152 * gi_30[k]
                  + f_120 * gi_37[k]
                  + f_152 * gi_44[k]
                  - f_120 * gi_46[k]
                  - f_120 * gi_170[k]
                  + f_121 * gi_177[k]
                  + f_120 * gi_184[k]
                  - f_121 * gi_186[k]
                  + f_153 * gi_226[k]
                  - f_122 * gi_233[k]
                  - f_153 * gi_240[k]
                  + f_122 * gi_242[k]
                  - f_152 * gi_284[k]
                  + f_120 * gi_293[k]
                  + f_152 * gi_302[k]
                  - f_120 * gi_304[k]
                  + f_153 * gi_340[k]
                  - f_122 * gi_349[k]
                  - f_153 * gi_358[k]
                  + f_122 * gi_360[k]
                  - f_124 * gi_396[k]
                  + f_123 * gi_405[k]
                  + f_124 * gi_414[k]
                  - f_123 * gi_416[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_26, gh_31, gh_33, gh_126, gh_129, \
                         gh_131, gh_136, gh_138, gh_168, gh_171, gh_173, gh_178, gh_180, \
                         gh_210, gh_213, gh_215, gh_220, gh_222, gh_252, gh_255, gh_257, \
                         gh_262, gh_264, gh_294, gh_297, gh_299, gh_304, gh_306, gi_28, gi_31, \
                         gi_33, gi_38, gi_40, gi_168, gi_171, gi_173, gi_178, gi_180, gi_224, \
                         gi_227, gi_229, gi_234, gi_236, gi_281, gi_286, gi_288, gi_295, \
                         gi_297, gi_337, gi_342, gi_344, gi_351, gi_353, gi_393, gi_398, \
                         gi_400, gi_407, gi_409 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_86 * ab_x[k] * gh_21[k]
                  - f_80 * ab_x[k] * gh_24[k]
                  - f_88 * ab_x[k] * gh_26[k]
                  - f_76 * ab_x[k] * gh_31[k]
                  + f_79 * ab_x[k] * gh_33[k]
                  + f_80 * ab_x[k] * gh_126[k]
                  - f_81 * ab_x[k] * gh_129[k]
                  - f_82 * ab_x[k] * gh_131[k]
                  - f_77 * ab_x[k] * gh_136[k]
                  + f_83 * ab_x[k] * gh_138[k]
                  - f_87 * ab_x[k] * gh_168[k]
                  + f_79 * ab_x[k] * gh_171[k]
                  + f_89 * ab_x[k] * gh_173[k]
                  + f_78 * ab_x[k] * gh_178[k]
                  - f_84 * ab_x[k] * gh_180[k]
                  + f_86 * ab_y[k] * gh_210[k]
                  - f_80 * ab_y[k] * gh_213[k]
                  - f_88 * ab_y[k] * gh_215[k]
                  - f_76 * ab_y[k] * gh_220[k]
                  + f_79 * ab_y[k] * gh_222[k]
                  - f_87 * ab_y[k] * gh_252[k]
                  + f_79 * ab_y[k] * gh_255[k]
                  + f_89 * ab_y[k] * gh_257[k]
                  + f_78 * ab_y[k] * gh_262[k]
                  - f_84 * ab_y[k] * gh_264[k]
                  + f_88 * ab_y[k] * gh_294[k]
                  - f_82 * ab_y[k] * gh_297[k]
                  - f_90 * ab_y[k] * gh_299[k]
                  - f_79 * ab_y[k] * gh_304[k]
                  + f_85 * ab_y[k] * gh_306[k]
                  - f_86 * gi_28[k]
                  + f_80 * gi_31[k]
                  + f_88 * gi_33[k]
                  + f_76 * gi_38[k]
                  - f_79 * gi_40[k]
                  - f_80 * gi_168[k]
                  + f_81 * gi_171[k]
                  + f_82 * gi_173[k]
                  + f_77 * gi_178[k]
                  - f_83 * gi_180[k]
                  + f_87 * gi_224[k]
                  - f_79 * gi_227[k]
                  - f_89 * gi_229[k]
                  - f_78 * gi_234[k]
                  + f_84 * gi_236[k]
                  - f_86 * gi_281[k]
                  + f_80 * gi_286[k]
                  + f_88 * gi_288[k]
                  + f_76 * gi_295[k]
                  - f_79 * gi_297[k]
                  + f_87 * gi_337[k]
                  - f_79 * gi_342[k]
                  - f_89 * gi_344[k]
                  - f_78 * gi_351[k]
                  + f_84 * gi_353[k]
                  - f_88 * gi_393[k]
                  + f_82 * gi_398[k]
                  + f_90 * gi_400[k]
                  + f_79 * gi_407[k]
                  - f_85 * gi_409[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_28, gh_37, gh_128, gh_133, gh_142, gh_170, \
                         gh_175, gh_184, gh_212, gh_217, gh_226, gh_254, gh_259, gh_268, \
                         gh_296, gh_301, gh_310, gi_30, gi_35, gi_44, gi_170, gi_175, gi_184, \
                         gi_226, gi_231, gi_240, gi_284, gi_291, gi_302, gi_340, gi_347, \
                         gi_358, gi_396, gi_403, gi_414 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = -f_154 * ab_x[k] * gh_23[k]
                  + f_155 * ab_x[k] * gh_28[k]
                  - f_154 * ab_x[k] * gh_37[k]
                  - f_156 * ab_x[k] * gh_128[k]
                  + f_157 * ab_x[k] * gh_133[k]
                  - f_156 * ab_x[k] * gh_142[k]
                  + f_157 * ab_x[k] * gh_170[k]
                  - f_158 * ab_x[k] * gh_175[k]
                  + f_157 * ab_x[k] * gh_184[k]
                  - f_154 * ab_y[k] * gh_212[k]
                  + f_155 * ab_y[k] * gh_217[k]
                  - f_154 * ab_y[k] * gh_226[k]
                  + f_157 * ab_y[k] * gh_254[k]
                  - f_158 * ab_y[k] * gh_259[k]
                  + f_157 * ab_y[k] * gh_268[k]
                  - f_59 * ab_y[k] * gh_296[k]
                  + f_60 * ab_y[k] * gh_301[k]
                  - f_59 * ab_y[k] * gh_310[k]
                  + f_154 * gi_30[k]
                  - f_155 * gi_35[k]
                  + f_154 * gi_44[k]
                  + f_156 * gi_170[k]
                  - f_157 * gi_175[k]
                  + f_156 * gi_184[k]
                  - f_157 * gi_226[k]
                  + f_158 * gi_231[k]
                  - f_157 * gi_240[k]
                  + f_154 * gi_284[k]
                  - f_155 * gi_291[k]
                  + f_154 * gi_302[k]
                  - f_157 * gi_340[k]
                  + f_158 * gi_347[k]
                  - f_157 * gi_358[k]
                  + f_59 * gi_396[k]
                  - f_60 * gi_403[k]
                  + f_59 * gi_414[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_31, gh_126, gh_129, gh_136, gh_168, \
                         gh_171, gh_178, gh_210, gh_213, gh_220, gh_252, gh_255, gh_262, \
                         gh_294, gh_297, gh_304, gi_28, gi_31, gi_38, gi_168, gi_171, gi_178, \
                         gi_224, gi_227, gi_234, gi_281, gi_286, gi_295, gi_337, gi_342, \
                         gi_351, gi_393, gi_398, gi_407 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_29 * ab_x[k] * gh_21[k]
                  + f_23 * ab_x[k] * gh_24[k]
                  - f_22 * ab_x[k] * gh_31[k]
                  - f_30 * ab_x[k] * gh_126[k]
                  + f_26 * ab_x[k] * gh_129[k]
                  - f_23 * ab_x[k] * gh_136[k]
                  + f_31 * ab_x[k] * gh_168[k]
                  - f_27 * ab_x[k] * gh_171[k]
                  + f_24 * ab_x[k] * gh_178[k]
                  - f_29 * ab_y[k] * gh_210[k]
                  + f_23 * ab_y[k] * gh_213[k]
                  - f_22 * ab_y[k] * gh_220[k]
                  + f_31 * ab_y[k] * gh_252[k]
                  - f_27 * ab_y[k] * gh_255[k]
                  + f_24 * ab_y[k] * gh_262[k]
                  - f_32 * ab_y[k] * gh_294[k]
                  + f_28 * ab_y[k] * gh_297[k]
                  - f_25 * ab_y[k] * gh_304[k]
                  + f_29 * gi_28[k]
                  - f_23 * gi_31[k]
                  + f_22 * gi_38[k]
                  + f_30 * gi_168[k]
                  - f_26 * gi_171[k]
                  + f_23 * gi_178[k]
                  - f_31 * gi_224[k]
                  + f_27 * gi_227[k]
                  - f_24 * gi_234[k]
                  + f_29 * gi_281[k]
                  - f_23 * gi_286[k]
                  + f_22 * gi_295[k]
                  - f_31 * gi_337[k]
                  + f_27 * gi_342[k]
                  - f_24 * gi_351[k]
                  + f_32 * gi_393[k]
                  - f_28 * gi_398[k]
                  + f_25 * gi_407[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_43, gh_48, gh_57, gh_148, gh_153, gh_162, \
                         gh_190, gh_195, gh_204, gh_232, gh_237, gh_246, gh_274, gh_279, \
                         gh_288, gh_295, gh_300, gh_309, gi_57, gi_62, gi_71, gi_197, gi_202, \
                         gi_211, gi_253, gi_258, gi_267, gi_311, gi_318, gi_329, gi_367, \
                         gi_374, gi_385, gi_396, gi_403, gi_414 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_33 * ab_x[k] * gh_43[k]
                  + f_34 * ab_x[k] * gh_48[k]
                  - f_40 * ab_x[k] * gh_57[k]
                  - f_34 * ab_x[k] * gh_148[k]
                  + f_37 * ab_x[k] * gh_153[k]
                  - f_41 * ab_x[k] * gh_162[k]
                  + f_35 * ab_x[k] * gh_190[k]
                  - f_38 * ab_x[k] * gh_195[k]
                  + f_36 * ab_x[k] * gh_204[k]
                  - f_33 * ab_y[k] * gh_232[k]
                  + f_34 * ab_y[k] * gh_237[k]
                  - f_40 * ab_y[k] * gh_246[k]
                  + f_35 * ab_y[k] * gh_274[k]
                  - f_38 * ab_y[k] * gh_279[k]
                  + f_36 * ab_y[k] * gh_288[k]
                  - f_36 * ab_z[k] * gh_295[k]
                  + f_39 * ab_z[k] * gh_300[k]
                  - f_42 * ab_z[k] * gh_309[k]
                  + f_33 * gi_57[k]
                  - f_34 * gi_62[k]
                  + f_40 * gi_71[k]
                  + f_34 * gi_197[k]
                  - f_37 * gi_202[k]
                  + f_41 * gi_211[k]
                  - f_35 * gi_253[k]
                  + f_38 * gi_258[k]
                  - f_36 * gi_267[k]
                  + f_33 * gi_311[k]
                  - f_34 * gi_318[k]
                  + f_40 * gi_329[k]
                  - f_35 * gi_367[k]
                  + f_38 * gi_374[k]
                  - f_36 * gi_385[k]
                  + f_36 * gi_396[k]
                  - f_39 * gi_403[k]
                  + f_42 * gi_414[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_46, gh_53, gh_151, gh_158, gh_193, gh_200, \
                         gh_235, gh_242, gh_277, gh_284, gh_298, gh_305, gi_60, gi_67, gi_200, \
                         gi_207, gi_256, gi_263, gi_315, gi_324, gi_371, gi_380, gi_400, \
                         gi_409 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = -f_62 * ab_x[k] * gh_46[k]
                  + f_62 * ab_x[k] * gh_53[k]
                  - f_63 * ab_x[k] * gh_151[k]
                  + f_63 * ab_x[k] * gh_158[k]
                  + f_64 * ab_x[k] * gh_193[k]
                  - f_64 * ab_x[k] * gh_200[k]
                  - f_62 * ab_y[k] * gh_235[k]
                  + f_62 * ab_y[k] * gh_242[k]
                  + f_64 * ab_y[k] * gh_277[k]
                  - f_64 * ab_y[k] * gh_284[k]
                  - f_65 * ab_z[k] * gh_298[k]
                  + f_65 * ab_z[k] * gh_305[k]
                  + f_62 * gi_60[k]
                  - f_62 * gi_67[k]
                  + f_63 * gi_200[k]
                  - f_63 * gi_207[k]
                  - f_64 * gi_256[k]
                  + f_64 * gi_263[k]
                  + f_62 * gi_315[k]
                  - f_62 * gi_324[k]
                  - f_64 * gi_371[k]
                  + f_64 * gi_380[k]
                  + f_65 * gi_400[k]
                  - f_65 * gi_409[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_43, gh_48, gh_50, gh_57, gh_59, gh_148, gh_153, \
                         gh_155, gh_162, gh_164, gh_190, gh_195, gh_197, gh_204, gh_206, \
                         gh_232, gh_237, gh_239, gh_246, gh_248, gh_274, gh_279, gh_281, \
                         gh_288, gh_290, gh_295, gh_300, gh_302, gh_309, gh_311, gi_57, gi_62, \
                         gi_64, gi_71, gi_73, gi_197, gi_202, gi_204, gi_211, gi_213, gi_253, \
                         gi_258, gi_260, gi_267, gi_269, gi_311, gi_318, gi_320, gi_329, \
                         gi_331, gi_367, gi_374, gi_376, gi_385, gi_387, gi_396, gi_403, \
                         gi_405, gi_414, gi_416 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_91 * ab_x[k] * gh_43[k]
                  + f_95 * ab_x[k] * gh_48[k]
                  - f_99 * ab_x[k] * gh_50[k]
                  - f_103 * ab_x[k] * gh_57[k]
                  + f_93 * ab_x[k] * gh_59[k]
                  + f_92 * ab_x[k] * gh_148[k]
                  + f_96 * ab_x[k] * gh_153[k]
                  - f_100 * ab_x[k] * gh_155[k]
                  - f_95 * ab_x[k] * gh_162[k]
                  + f_106 * ab_x[k] * gh_164[k]
                  - f_93 * ab_x[k] * gh_190[k]
                  - f_97 * ab_x[k] * gh_195[k]
                  + f_101 * ab_x[k] * gh_197[k]
                  + f_104 * ab_x[k] * gh_204[k]
                  - f_107 * ab_x[k] * gh_206[k]
                  + f_91 * ab_y[k] * gh_232[k]
                  + f_95 * ab_y[k] * gh_237[k]
                  - f_99 * ab_y[k] * gh_239[k]
                  - f_103 * ab_y[k] * gh_246[k]
                  + f_93 * ab_y[k] * gh_248[k]
                  - f_93 * ab_y[k] * gh_274[k]
                  - f_97 * ab_y[k] * gh_279[k]
                  + f_101 * ab_y[k] * gh_281[k]
                  + f_104 * ab_y[k] * gh_288[k]
                  - f_107 * ab_y[k] * gh_290[k]
                  + f_94 * ab_z[k] * gh_295[k]
                  + f_98 * ab_z[k] * gh_300[k]
                  - f_102 * ab_z[k] * gh_302[k]
                  - f_105 * ab_z[k] * gh_309[k]
                  + f_108 * ab_z[k] * gh_311[k]
                  - f_91 * gi_57[k]
                  - f_95 * gi_62[k]
                  + f_99 * gi_64[k]
                  + f_103 * gi_71[k]
                  - f_93 * gi_73[k]
                  - f_92 * gi_197[k]
                  - f_96 * gi_202[k]
                  + f_100 * gi_204[k]
                  + f_95 * gi_211[k]
                  - f_106 * gi_213[k]
                  + f_93 * gi_253[k]
                  + f_97 * gi_258[k]
                  - f_101 * gi_260[k]
                  - f_104 * gi_267[k]
                  + f_107 * gi_269[k]
                  - f_91 * gi_311[k]
                  - f_95 * gi_318[k]
                  + f_99 * gi_320[k]
                  + f_103 * gi_329[k]
                  - f_93 * gi_331[k]
                  + f_93 * gi_367[k]
                  + f_97 * gi_374[k]
                  - f_101 * gi_376[k]
                  - f_104 * gi_385[k]
                  + f_107 * gi_387[k]
                  - f_94 * gi_396[k]
                  - f_98 * gi_403[k]
                  + f_102 * gi_405[k]
                  + f_105 * gi_414[k]
                  - f_108 * gi_416[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_46, gh_53, gh_55, gh_151, gh_158, gh_160, \
                         gh_193, gh_200, gh_202, gh_235, gh_242, gh_244, gh_277, gh_284, \
                         gh_286, gh_298, gh_305, gh_307, gi_60, gi_67, gi_69, gi_200, gi_207, \
                         gi_209, gi_256, gi_263, gi_265, gi_315, gi_324, gi_326, gi_371, \
                         gi_380, gi_382, gi_400, gi_409, gi_411 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = f_127 * ab_x[k] * gh_46[k]
                  + f_127 * ab_x[k] * gh_53[k]
                  - f_128 * ab_x[k] * gh_55[k]
                  + f_128 * ab_x[k] * gh_151[k]
                  + f_128 * ab_x[k] * gh_158[k]
                  - f_131 * ab_x[k] * gh_160[k]
                  - f_129 * ab_x[k] * gh_193[k]
                  - f_129 * ab_x[k] * gh_200[k]
                  + f_132 * ab_x[k] * gh_202[k]
                  + f_127 * ab_y[k] * gh_235[k]
                  + f_127 * ab_y[k] * gh_242[k]
                  - f_128 * ab_y[k] * gh_244[k]
                  - f_129 * ab_y[k] * gh_277[k]
                  - f_129 * ab_y[k] * gh_284[k]
                  + f_132 * ab_y[k] * gh_286[k]
                  + f_130 * ab_z[k] * gh_298[k]
                  + f_130 * ab_z[k] * gh_305[k]
                  - f_133 * ab_z[k] * gh_307[k]
                  - f_127 * gi_60[k]
                  - f_127 * gi_67[k]
                  + f_128 * gi_69[k]
                  - f_128 * gi_200[k]
                  - f_128 * gi_207[k]
                  + f_131 * gi_209[k]
                  + f_129 * gi_256[k]
                  + f_129 * gi_263[k]
                  - f_132 * gi_265[k]
                  - f_127 * gi_315[k]
                  - f_127 * gi_324[k]
                  + f_128 * gi_326[k]
                  + f_129 * gi_371[k]
                  + f_129 * gi_380[k]
                  - f_132 * gi_382[k]
                  - f_130 * gi_400[k]
                  - f_130 * gi_409[k]
                  + f_133 * gi_411[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_43, gh_48, gh_50, gh_57, gh_59, gh_61, gh_148, \
                         gh_153, gh_155, gh_162, gh_164, gh_166, gh_190, gh_195, gh_197, \
                         gh_204, gh_206, gh_208, gh_232, gh_237, gh_239, gh_246, gh_248, \
                         gh_250, gh_274, gh_279, gh_281, gh_288, gh_290, gh_292, gh_295, \
                         gh_300, gh_302, gh_309, gh_311, gh_313, gi_57, gi_62, gi_64, gi_71, \
                         gi_73, gi_75, gi_197, gi_202, gi_204, gi_211, gi_213, gi_215, gi_253, \
                         gi_258, gi_260, gi_267, gi_269, gi_271, gi_311, gi_318, gi_320, \
                         gi_329, gi_331, gi_333, gi_367, gi_374, gi_376, gi_385, gi_387, \
                         gi_389, gi_396, gi_403, gi_405, gi_414, gi_416, \
                         gi_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_137 * ab_x[k] * gh_43[k]
                  - f_138 * ab_x[k] * gh_48[k]
                  + f_144 * ab_x[k] * gh_50[k]
                  - f_137 * ab_x[k] * gh_57[k]
                  + f_144 * ab_x[k] * gh_59[k]
                  - f_148 * ab_x[k] * gh_61[k]
                  - f_138 * ab_x[k] * gh_148[k]
                  - f_141 * ab_x[k] * gh_153[k]
                  + f_145 * ab_x[k] * gh_155[k]
                  - f_138 * ab_x[k] * gh_162[k]
                  + f_145 * ab_x[k] * gh_164[k]
                  - f_149 * ab_x[k] * gh_166[k]
                  + f_139 * ab_x[k] * gh_190[k]
                  + f_142 * ab_x[k] * gh_195[k]
                  - f_146 * ab_x[k] * gh_197[k]
                  + f_139 * ab_x[k] * gh_204[k]
                  - f_146 * ab_x[k] * gh_206[k]
                  + f_150 * ab_x[k] * gh_208[k]
                  - f_137 * ab_y[k] * gh_232[k]
                  - f_138 * ab_y[k] * gh_237[k]
                  + f_144 * ab_y[k] * gh_239[k]
                  - f_137 * ab_y[k] * gh_246[k]
                  + f_144 * ab_y[k] * gh_248[k]
                  - f_148 * ab_y[k] * gh_250[k]
                  + f_139 * ab_y[k] * gh_274[k]
                  + f_142 * ab_y[k] * gh_279[k]
                  - f_146 * ab_y[k] * gh_281[k]
                  + f_139 * ab_y[k] * gh_288[k]
                  - f_146 * ab_y[k] * gh_290[k]
                  + f_150 * ab_y[k] * gh_292[k]
                  - f_140 * ab_z[k] * gh_295[k]
                  - f_143 * ab_z[k] * gh_300[k]
                  + f_147 * ab_z[k] * gh_302[k]
                  - f_140 * ab_z[k] * gh_309[k]
                  + f_147 * ab_z[k] * gh_311[k]
                  - f_151 * ab_z[k] * gh_313[k]
                  + f_137 * gi_57[k]
                  + f_138 * gi_62[k]
                  - f_144 * gi_64[k]
                  + f_137 * gi_71[k]
                  - f_144 * gi_73[k]
                  + f_148 * gi_75[k]
                  + f_138 * gi_197[k]
                  + f_141 * gi_202[k]
                  - f_145 * gi_204[k]
                  + f_138 * gi_211[k]
                  - f_145 * gi_213[k]
                  + f_149 * gi_215[k]
                  - f_139 * gi_253[k]
                  - f_142 * gi_258[k]
                  + f_146 * gi_260[k]
                  - f_139 * gi_267[k]
                  + f_146 * gi_269[k]
                  - f_150 * gi_271[k]
                  + f_137 * gi_311[k]
                  + f_138 * gi_318[k]
                  - f_144 * gi_320[k]
                  + f_137 * gi_329[k]
                  - f_144 * gi_331[k]
                  + f_148 * gi_333[k]
                  - f_139 * gi_367[k]
                  - f_142 * gi_374[k]
                  + f_146 * gi_376[k]
                  - f_139 * gi_385[k]
                  + f_146 * gi_387[k]
                  - f_150 * gi_389[k]
                  + f_140 * gi_396[k]
                  + f_143 * gi_403[k]
                  - f_147 * gi_405[k]
                  + f_140 * gi_414[k]
                  - f_147 * gi_416[k]
                  + f_151 * gi_418[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_44, gh_49, gh_51, gh_58, gh_60, gh_62, gh_149, \
                         gh_154, gh_156, gh_163, gh_165, gh_167, gh_191, gh_196, gh_198, \
                         gh_205, gh_207, gh_209, gh_233, gh_238, gh_240, gh_247, gh_249, \
                         gh_251, gh_275, gh_280, gh_282, gh_289, gh_291, gh_293, gh_296, \
                         gh_301, gh_303, gh_310, gh_312, gh_314, gi_58, gi_63, gi_65, gi_72, \
                         gi_74, gi_76, gi_198, gi_203, gi_205, gi_212, gi_214, gi_216, gi_254, \
                         gi_259, gi_261, gi_268, gi_270, gi_272, gi_312, gi_319, gi_321, \
                         gi_330, gi_332, gi_334, gi_368, gi_375, gi_377, gi_386, gi_388, \
                         gi_390, gi_397, gi_404, gi_406, gi_415, gi_417, \
                         gi_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -3.515625 * ab_x[k] * gh_44[k]
                  - 7.03125 * ab_x[k] * gh_49[k]
                  + 9.375 * ab_x[k] * gh_51[k]
                  - 3.515625 * ab_x[k] * gh_58[k]
                  + 9.375 * ab_x[k] * gh_60[k]
                  - 1.875 * ab_x[k] * gh_62[k]
                  - 7.03125 * ab_x[k] * gh_149[k]
                  - 14.0625 * ab_x[k] * gh_154[k]
                  + 18.75 * ab_x[k] * gh_156[k]
                  - 7.03125 * ab_x[k] * gh_163[k]
                  + 18.75 * ab_x[k] * gh_165[k]
                  - 3.75 * ab_x[k] * gh_167[k]
                  + 9.375 * ab_x[k] * gh_191[k]
                  + 18.75 * ab_x[k] * gh_196[k]
                  - 25.0 * ab_x[k] * gh_198[k]
                  + 9.375 * ab_x[k] * gh_205[k]
                  - 25.0 * ab_x[k] * gh_207[k]
                  + 5.0 * ab_x[k] * gh_209[k]
                  - 3.515625 * ab_y[k] * gh_233[k]
                  - 7.03125 * ab_y[k] * gh_238[k]
                  + 9.375 * ab_y[k] * gh_240[k]
                  - 3.515625 * ab_y[k] * gh_247[k]
                  + 9.375 * ab_y[k] * gh_249[k]
                  - 1.875 * ab_y[k] * gh_251[k]
                  + 9.375 * ab_y[k] * gh_275[k]
                  + 18.75 * ab_y[k] * gh_280[k]
                  - 25.0 * ab_y[k] * gh_282[k]
                  + 9.375 * ab_y[k] * gh_289[k]
                  - 25.0 * ab_y[k] * gh_291[k]
                  + 5.0 * ab_y[k] * gh_293[k]
                  - 1.875 * ab_z[k] * gh_296[k]
                  - 3.75 * ab_z[k] * gh_301[k]
                  + 5.0 * ab_z[k] * gh_303[k]
                  - 1.875 * ab_z[k] * gh_310[k]
                  + 5.0 * ab_z[k] * gh_312[k]
                  - ab_z[k] * gh_314[k]
                  + 3.515625 * gi_58[k]
                  + 7.03125 * gi_63[k]
                  - 9.375 * gi_65[k]
                  + 3.515625 * gi_72[k]
                  - 9.375 * gi_74[k]
                  + 1.875 * gi_76[k]
                  + 7.03125 * gi_198[k]
                  + 14.0625 * gi_203[k]
                  - 18.75 * gi_205[k]
                  + 7.03125 * gi_212[k]
                  - 18.75 * gi_214[k]
                  + 3.75 * gi_216[k]
                  - 9.375 * gi_254[k]
                  - 18.75 * gi_259[k]
                  + 25.0 * gi_261[k]
                  - 9.375 * gi_268[k]
                  + 25.0 * gi_270[k]
                  - 5.0 * gi_272[k]
                  + 3.515625 * gi_312[k]
                  + 7.03125 * gi_319[k]
                  - 9.375 * gi_321[k]
                  + 3.515625 * gi_330[k]
                  - 9.375 * gi_332[k]
                  + 1.875 * gi_334[k]
                  - 9.375 * gi_368[k]
                  - 18.75 * gi_375[k]
                  + 25.0 * gi_377[k]
                  - 9.375 * gi_386[k]
                  + 25.0 * gi_388[k]
                  - 5.0 * gi_390[k]
                  + 1.875 * gi_397[k]
                  + 3.75 * gi_404[k]
                  - 5.0 * gi_406[k]
                  + 1.875 * gi_415[k]
                  - 5.0 * gi_417[k]
                  + gi_419[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_42, gh_45, gh_47, gh_52, gh_54, gh_56, gh_147, \
                         gh_150, gh_152, gh_157, gh_159, gh_161, gh_189, gh_192, gh_194, \
                         gh_199, gh_201, gh_203, gh_231, gh_234, gh_236, gh_241, gh_243, \
                         gh_245, gh_273, gh_276, gh_278, gh_283, gh_285, gh_287, gh_294, \
                         gh_297, gh_299, gh_304, gh_306, gh_308, gi_56, gi_59, gi_61, gi_66, \
                         gi_68, gi_70, gi_196, gi_199, gi_201, gi_206, gi_208, gi_210, gi_252, \
                         gi_255, gi_257, gi_262, gi_264, gi_266, gi_309, gi_314, gi_316, \
                         gi_323, gi_325, gi_327, gi_365, gi_370, gi_372, gi_379, gi_381, \
                         gi_383, gi_394, gi_399, gi_401, gi_408, gi_410, \
                         gi_412 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_137 * ab_x[k] * gh_42[k]
                  - f_138 * ab_x[k] * gh_45[k]
                  + f_144 * ab_x[k] * gh_47[k]
                  - f_137 * ab_x[k] * gh_52[k]
                  + f_144 * ab_x[k] * gh_54[k]
                  - f_148 * ab_x[k] * gh_56[k]
                  - f_138 * ab_x[k] * gh_147[k]
                  - f_141 * ab_x[k] * gh_150[k]
                  + f_145 * ab_x[k] * gh_152[k]
                  - f_138 * ab_x[k] * gh_157[k]
                  + f_145 * ab_x[k] * gh_159[k]
                  - f_149 * ab_x[k] * gh_161[k]
                  + f_139 * ab_x[k] * gh_189[k]
                  + f_142 * ab_x[k] * gh_192[k]
                  - f_146 * ab_x[k] * gh_194[k]
                  + f_139 * ab_x[k] * gh_199[k]
                  - f_146 * ab_x[k] * gh_201[k]
                  + f_150 * ab_x[k] * gh_203[k]
                  - f_137 * ab_y[k] * gh_231[k]
                  - f_138 * ab_y[k] * gh_234[k]
                  + f_144 * ab_y[k] * gh_236[k]
                  - f_137 * ab_y[k] * gh_241[k]
                  + f_144 * ab_y[k] * gh_243[k]
                  - f_148 * ab_y[k] * gh_245[k]
                  + f_139 * ab_y[k] * gh_273[k]
                  + f_142 * ab_y[k] * gh_276[k]
                  - f_146 * ab_y[k] * gh_278[k]
                  + f_139 * ab_y[k] * gh_283[k]
                  - f_146 * ab_y[k] * gh_285[k]
                  + f_150 * ab_y[k] * gh_287[k]
                  - f_140 * ab_z[k] * gh_294[k]
                  - f_143 * ab_z[k] * gh_297[k]
                  + f_147 * ab_z[k] * gh_299[k]
                  - f_140 * ab_z[k] * gh_304[k]
                  + f_147 * ab_z[k] * gh_306[k]
                  - f_151 * ab_z[k] * gh_308[k]
                  + f_137 * gi_56[k]
                  + f_138 * gi_59[k]
                  - f_144 * gi_61[k]
                  + f_137 * gi_66[k]
                  - f_144 * gi_68[k]
                  + f_148 * gi_70[k]
                  + f_138 * gi_196[k]
                  + f_141 * gi_199[k]
                  - f_145 * gi_201[k]
                  + f_138 * gi_206[k]
                  - f_145 * gi_208[k]
                  + f_149 * gi_210[k]
                  - f_139 * gi_252[k]
                  - f_142 * gi_255[k]
                  + f_146 * gi_257[k]
                  - f_139 * gi_262[k]
                  + f_146 * gi_264[k]
                  - f_150 * gi_266[k]
                  + f_137 * gi_309[k]
                  + f_138 * gi_314[k]
                  - f_144 * gi_316[k]
                  + f_137 * gi_323[k]
                  - f_144 * gi_325[k]
                  + f_148 * gi_327[k]
                  - f_139 * gi_365[k]
                  - f_142 * gi_370[k]
                  + f_146 * gi_372[k]
                  - f_139 * gi_379[k]
                  + f_146 * gi_381[k]
                  - f_150 * gi_383[k]
                  + f_140 * gi_394[k]
                  + f_143 * gi_399[k]
                  - f_147 * gi_401[k]
                  + f_140 * gi_408[k]
                  - f_147 * gi_410[k]
                  + f_151 * gi_412[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_44, gh_51, gh_58, gh_60, gh_149, gh_156, gh_163, \
                         gh_165, gh_191, gh_198, gh_205, gh_207, gh_233, gh_240, gh_247, \
                         gh_249, gh_275, gh_282, gh_289, gh_291, gh_296, gh_303, gh_310, \
                         gh_312, gi_58, gi_65, gi_72, gi_74, gi_198, gi_205, gi_212, gi_214, \
                         gi_254, gi_261, gi_268, gi_270, gi_312, gi_321, gi_330, gi_332, \
                         gi_368, gi_377, gi_386, gi_388, gi_397, gi_406, gi_415, \
                         gi_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_159 * ab_x[k] * gh_44[k]
                  - f_127 * ab_x[k] * gh_51[k]
                  - f_159 * ab_x[k] * gh_58[k]
                  + f_127 * ab_x[k] * gh_60[k]
                  + f_127 * ab_x[k] * gh_149[k]
                  - f_128 * ab_x[k] * gh_156[k]
                  - f_127 * ab_x[k] * gh_163[k]
                  + f_128 * ab_x[k] * gh_165[k]
                  - f_160 * ab_x[k] * gh_191[k]
                  + f_129 * ab_x[k] * gh_198[k]
                  + f_160 * ab_x[k] * gh_205[k]
                  - f_129 * ab_x[k] * gh_207[k]
                  + f_159 * ab_y[k] * gh_233[k]
                  - f_127 * ab_y[k] * gh_240[k]
                  - f_159 * ab_y[k] * gh_247[k]
                  + f_127 * ab_y[k] * gh_249[k]
                  - f_160 * ab_y[k] * gh_275[k]
                  + f_129 * ab_y[k] * gh_282[k]
                  + f_160 * ab_y[k] * gh_289[k]
                  - f_129 * ab_y[k] * gh_291[k]
                  + f_161 * ab_z[k] * gh_296[k]
                  - f_130 * ab_z[k] * gh_303[k]
                  - f_161 * ab_z[k] * gh_310[k]
                  + f_130 * ab_z[k] * gh_312[k]
                  - f_159 * gi_58[k]
                  + f_127 * gi_65[k]
                  + f_159 * gi_72[k]
                  - f_127 * gi_74[k]
                  - f_127 * gi_198[k]
                  + f_128 * gi_205[k]
                  + f_127 * gi_212[k]
                  - f_128 * gi_214[k]
                  + f_160 * gi_254[k]
                  - f_129 * gi_261[k]
                  - f_160 * gi_268[k]
                  + f_129 * gi_270[k]
                  - f_159 * gi_312[k]
                  + f_127 * gi_321[k]
                  + f_159 * gi_330[k]
                  - f_127 * gi_332[k]
                  + f_160 * gi_368[k]
                  - f_129 * gi_377[k]
                  - f_160 * gi_386[k]
                  + f_129 * gi_388[k]
                  - f_161 * gi_397[k]
                  + f_130 * gi_406[k]
                  + f_161 * gi_415[k]
                  - f_130 * gi_417[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_42, gh_45, gh_47, gh_52, gh_54, gh_147, gh_150, \
                         gh_152, gh_157, gh_159, gh_189, gh_192, gh_194, gh_199, gh_201, \
                         gh_231, gh_234, gh_236, gh_241, gh_243, gh_273, gh_276, gh_278, \
                         gh_283, gh_285, gh_294, gh_297, gh_299, gh_304, gh_306, gi_56, gi_59, \
                         gi_61, gi_66, gi_68, gi_196, gi_199, gi_201, gi_206, gi_208, gi_252, \
                         gi_255, gi_257, gi_262, gi_264, gi_309, gi_314, gi_316, gi_323, \
                         gi_325, gi_365, gi_370, gi_372, gi_379, gi_381, gi_394, gi_399, \
                         gi_401, gi_408, gi_410 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_103 * ab_x[k] * gh_42[k]
                  - f_95 * ab_x[k] * gh_45[k]
                  - f_93 * ab_x[k] * gh_47[k]
                  - f_91 * ab_x[k] * gh_52[k]
                  + f_99 * ab_x[k] * gh_54[k]
                  + f_95 * ab_x[k] * gh_147[k]
                  - f_96 * ab_x[k] * gh_150[k]
                  - f_106 * ab_x[k] * gh_152[k]
                  - f_92 * ab_x[k] * gh_157[k]
                  + f_100 * ab_x[k] * gh_159[k]
                  - f_104 * ab_x[k] * gh_189[k]
                  + f_97 * ab_x[k] * gh_192[k]
                  + f_107 * ab_x[k] * gh_194[k]
                  + f_93 * ab_x[k] * gh_199[k]
                  - f_101 * ab_x[k] * gh_201[k]
                  + f_103 * ab_y[k] * gh_231[k]
                  - f_95 * ab_y[k] * gh_234[k]
                  - f_93 * ab_y[k] * gh_236[k]
                  - f_91 * ab_y[k] * gh_241[k]
                  + f_99 * ab_y[k] * gh_243[k]
                  - f_104 * ab_y[k] * gh_273[k]
                  + f_97 * ab_y[k] * gh_276[k]
                  + f_107 * ab_y[k] * gh_278[k]
                  + f_93 * ab_y[k] * gh_283[k]
                  - f_101 * ab_y[k] * gh_285[k]
                  + f_105 * ab_z[k] * gh_294[k]
                  - f_98 * ab_z[k] * gh_297[k]
                  - f_108 * ab_z[k] * gh_299[k]
                  - f_94 * ab_z[k] * gh_304[k]
                  + f_102 * ab_z[k] * gh_306[k]
                  - f_103 * gi_56[k]
                  + f_95 * gi_59[k]
                  + f_93 * gi_61[k]
                  + f_91 * gi_66[k]
                  - f_99 * gi_68[k]
                  - f_95 * gi_196[k]
                  + f_96 * gi_199[k]
                  + f_106 * gi_201[k]
                  + f_92 * gi_206[k]
                  - f_100 * gi_208[k]
                  + f_104 * gi_252[k]
                  - f_97 * gi_255[k]
                  - f_107 * gi_257[k]
                  - f_93 * gi_262[k]
                  + f_101 * gi_264[k]
                  - f_103 * gi_309[k]
                  + f_95 * gi_314[k]
                  + f_93 * gi_316[k]
                  + f_91 * gi_323[k]
                  - f_99 * gi_325[k]
                  + f_104 * gi_365[k]
                  - f_97 * gi_370[k]
                  - f_107 * gi_372[k]
                  - f_93 * gi_379[k]
                  + f_101 * gi_381[k]
                  - f_105 * gi_394[k]
                  + f_98 * gi_399[k]
                  + f_108 * gi_401[k]
                  + f_94 * gi_408[k]
                  - f_102 * gi_410[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_44, gh_49, gh_58, gh_149, gh_154, gh_163, \
                         gh_191, gh_196, gh_205, gh_233, gh_238, gh_247, gh_275, gh_280, \
                         gh_289, gh_296, gh_301, gh_310, gi_58, gi_63, gi_72, gi_198, gi_203, \
                         gi_212, gi_254, gi_259, gi_268, gi_312, gi_319, gi_330, gi_368, \
                         gi_375, gi_386, gi_397, gi_404, gi_415 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_162 * ab_x[k] * gh_44[k]
                  + f_163 * ab_x[k] * gh_49[k]
                  - f_162 * ab_x[k] * gh_58[k]
                  - f_164 * ab_x[k] * gh_149[k]
                  + f_165 * ab_x[k] * gh_154[k]
                  - f_164 * ab_x[k] * gh_163[k]
                  + f_166 * ab_x[k] * gh_191[k]
                  - f_167 * ab_x[k] * gh_196[k]
                  + f_166 * ab_x[k] * gh_205[k]
                  - f_162 * ab_y[k] * gh_233[k]
                  + f_163 * ab_y[k] * gh_238[k]
                  - f_162 * ab_y[k] * gh_247[k]
                  + f_166 * ab_y[k] * gh_275[k]
                  - f_167 * ab_y[k] * gh_280[k]
                  + f_166 * ab_y[k] * gh_289[k]
                  - f_168 * ab_z[k] * gh_296[k]
                  + f_169 * ab_z[k] * gh_301[k]
                  - f_168 * ab_z[k] * gh_310[k]
                  + f_162 * gi_58[k]
                  - f_163 * gi_63[k]
                  + f_162 * gi_72[k]
                  + f_164 * gi_198[k]
                  - f_165 * gi_203[k]
                  + f_164 * gi_212[k]
                  - f_166 * gi_254[k]
                  + f_167 * gi_259[k]
                  - f_166 * gi_268[k]
                  + f_162 * gi_312[k]
                  - f_163 * gi_319[k]
                  + f_162 * gi_330[k]
                  - f_166 * gi_368[k]
                  + f_167 * gi_375[k]
                  - f_166 * gi_386[k]
                  + f_168 * gi_397[k]
                  - f_169 * gi_404[k]
                  + f_168 * gi_415[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_42, gh_45, gh_52, gh_147, gh_150, gh_157, \
                         gh_189, gh_192, gh_199, gh_231, gh_234, gh_241, gh_273, gh_276, \
                         gh_283, gh_294, gh_297, gh_304, gi_56, gi_59, gi_66, gi_196, gi_199, \
                         gi_206, gi_252, gi_255, gi_262, gi_309, gi_314, gi_323, gi_365, \
                         gi_370, gi_379, gi_394, gi_399, gi_408 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = -f_40 * ab_x[k] * gh_42[k]
                  + f_34 * ab_x[k] * gh_45[k]
                  - f_33 * ab_x[k] * gh_52[k]
                  - f_41 * ab_x[k] * gh_147[k]
                  + f_37 * ab_x[k] * gh_150[k]
                  - f_34 * ab_x[k] * gh_157[k]
                  + f_36 * ab_x[k] * gh_189[k]
                  - f_38 * ab_x[k] * gh_192[k]
                  + f_35 * ab_x[k] * gh_199[k]
                  - f_40 * ab_y[k] * gh_231[k]
                  + f_34 * ab_y[k] * gh_234[k]
                  - f_33 * ab_y[k] * gh_241[k]
                  + f_36 * ab_y[k] * gh_273[k]
                  - f_38 * ab_y[k] * gh_276[k]
                  + f_35 * ab_y[k] * gh_283[k]
                  - f_42 * ab_z[k] * gh_294[k]
                  + f_39 * ab_z[k] * gh_297[k]
                  - f_36 * ab_z[k] * gh_304[k]
                  + f_40 * gi_56[k]
                  - f_34 * gi_59[k]
                  + f_33 * gi_66[k]
                  + f_41 * gi_196[k]
                  - f_37 * gi_199[k]
                  + f_34 * gi_206[k]
                  - f_36 * gi_252[k]
                  + f_38 * gi_255[k]
                  - f_35 * gi_262[k]
                  + f_40 * gi_309[k]
                  - f_34 * gi_314[k]
                  + f_33 * gi_323[k]
                  - f_36 * gi_365[k]
                  + f_38 * gi_370[k]
                  - f_35 * gi_379[k]
                  + f_42 * gi_394[k]
                  - f_39 * gi_399[k]
                  + f_36 * gi_408[k];
    }

#pragma omp simd aligned(ab_x, gh_1, gh_6, gh_15, gh_64, gh_69, gh_78, gh_106, gh_111, gh_120, \
                         gh_211, gh_216, gh_225, gh_253, gh_258, gh_267, gh_295, gh_300, \
                         gh_309, gi_1, gi_6, gi_15, gi_85, gi_90, gi_99, gi_141, gi_146, \
                         gi_155, gi_281, gi_286, gi_295, gi_337, gi_342, gi_351, gi_393, \
                         gi_398, gi_407 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = -f_22 * ab_x[k] * gh_1[k]
                  + f_23 * ab_x[k] * gh_6[k]
                  - f_29 * ab_x[k] * gh_15[k]
                  - f_23 * ab_x[k] * gh_64[k]
                  + f_26 * ab_x[k] * gh_69[k]
                  - f_30 * ab_x[k] * gh_78[k]
                  + f_24 * ab_x[k] * gh_106[k]
                  - f_27 * ab_x[k] * gh_111[k]
                  + f_31 * ab_x[k] * gh_120[k]
                  - f_22 * ab_x[k] * gh_211[k]
                  + f_23 * ab_x[k] * gh_216[k]
                  - f_29 * ab_x[k] * gh_225[k]
                  + f_24 * ab_x[k] * gh_253[k]
                  - f_27 * ab_x[k] * gh_258[k]
                  + f_31 * ab_x[k] * gh_267[k]
                  - f_25 * ab_x[k] * gh_295[k]
                  + f_28 * ab_x[k] * gh_300[k]
                  - f_32 * ab_x[k] * gh_309[k]
                  + f_22 * gi_1[k]
                  - f_23 * gi_6[k]
                  + f_29 * gi_15[k]
                  + f_23 * gi_85[k]
                  - f_26 * gi_90[k]
                  + f_30 * gi_99[k]
                  - f_24 * gi_141[k]
                  + f_27 * gi_146[k]
                  - f_31 * gi_155[k]
                  + f_22 * gi_281[k]
                  - f_23 * gi_286[k]
                  + f_29 * gi_295[k]
                  - f_24 * gi_337[k]
                  + f_27 * gi_342[k]
                  - f_31 * gi_351[k]
                  + f_25 * gi_393[k]
                  - f_28 * gi_398[k]
                  + f_32 * gi_407[k];
    }

#pragma omp simd aligned(ab_x, gh_4, gh_11, gh_67, gh_74, gh_109, gh_116, gh_214, gh_221, \
                         gh_256, gh_263, gh_298, gh_305, gi_4, gi_11, gi_88, gi_95, gi_144, \
                         gi_151, gi_284, gi_291, gi_340, gi_347, gi_396, \
                         gi_403 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_58 * ab_x[k] * gh_4[k]
                  + f_58 * ab_x[k] * gh_11[k]
                  - f_59 * ab_x[k] * gh_67[k]
                  + f_59 * ab_x[k] * gh_74[k]
                  + f_60 * ab_x[k] * gh_109[k]
                  - f_60 * ab_x[k] * gh_116[k]
                  - f_58 * ab_x[k] * gh_214[k]
                  + f_58 * ab_x[k] * gh_221[k]
                  + f_60 * ab_x[k] * gh_256[k]
                  - f_60 * ab_x[k] * gh_263[k]
                  - f_61 * ab_x[k] * gh_298[k]
                  + f_61 * ab_x[k] * gh_305[k]
                  + f_58 * gi_4[k]
                  - f_58 * gi_11[k]
                  + f_59 * gi_88[k]
                  - f_59 * gi_95[k]
                  - f_60 * gi_144[k]
                  + f_60 * gi_151[k]
                  + f_58 * gi_284[k]
                  - f_58 * gi_291[k]
                  - f_60 * gi_340[k]
                  + f_60 * gi_347[k]
                  + f_61 * gi_396[k]
                  - f_61 * gi_403[k];
    }

#pragma omp simd aligned(ab_x, gh_1, gh_6, gh_8, gh_15, gh_17, gh_64, gh_69, gh_71, gh_78, \
                         gh_80, gh_106, gh_111, gh_113, gh_120, gh_122, gh_211, gh_216, \
                         gh_218, gh_225, gh_227, gh_253, gh_258, gh_260, gh_267, gh_269, \
                         gh_295, gh_300, gh_302, gh_309, gh_311, gi_1, gi_6, gi_8, gi_15, \
                         gi_17, gi_85, gi_90, gi_92, gi_99, gi_101, gi_141, gi_146, gi_148, \
                         gi_155, gi_157, gi_281, gi_286, gi_288, gi_295, gi_297, gi_337, \
                         gi_342, gi_344, gi_351, gi_353, gi_393, gi_398, gi_400, gi_407, \
                         gi_409 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_76 * ab_x[k] * gh_1[k]
                  + f_80 * ab_x[k] * gh_6[k]
                  - f_79 * ab_x[k] * gh_8[k]
                  - f_86 * ab_x[k] * gh_15[k]
                  + f_88 * ab_x[k] * gh_17[k]
                  + f_77 * ab_x[k] * gh_64[k]
                  + f_81 * ab_x[k] * gh_69[k]
                  - f_83 * ab_x[k] * gh_71[k]
                  - f_80 * ab_x[k] * gh_78[k]
                  + f_82 * ab_x[k] * gh_80[k]
                  - f_78 * ab_x[k] * gh_106[k]
                  - f_79 * ab_x[k] * gh_111[k]
                  + f_84 * ab_x[k] * gh_113[k]
                  + f_87 * ab_x[k] * gh_120[k]
                  - f_89 * ab_x[k] * gh_122[k]
                  + f_76 * ab_x[k] * gh_211[k]
                  + f_80 * ab_x[k] * gh_216[k]
                  - f_79 * ab_x[k] * gh_218[k]
                  - f_86 * ab_x[k] * gh_225[k]
                  + f_88 * ab_x[k] * gh_227[k]
                  - f_78 * ab_x[k] * gh_253[k]
                  - f_79 * ab_x[k] * gh_258[k]
                  + f_84 * ab_x[k] * gh_260[k]
                  + f_87 * ab_x[k] * gh_267[k]
                  - f_89 * ab_x[k] * gh_269[k]
                  + f_79 * ab_x[k] * gh_295[k]
                  + f_82 * ab_x[k] * gh_300[k]
                  - f_85 * ab_x[k] * gh_302[k]
                  - f_88 * ab_x[k] * gh_309[k]
                  + f_90 * ab_x[k] * gh_311[k]
                  - f_76 * gi_1[k]
                  - f_80 * gi_6[k]
                  + f_79 * gi_8[k]
                  + f_86 * gi_15[k]
                  - f_88 * gi_17[k]
                  - f_77 * gi_85[k]
                  - f_81 * gi_90[k]
                  + f_83 * gi_92[k]
                  + f_80 * gi_99[k]
                  - f_82 * gi_101[k]
                  + f_78 * gi_141[k]
                  + f_79 * gi_146[k]
                  - f_84 * gi_148[k]
                  - f_87 * gi_155[k]
                  + f_89 * gi_157[k]
                  - f_76 * gi_281[k]
                  - f_80 * gi_286[k]
                  + f_79 * gi_288[k]
                  + f_86 * gi_295[k]
                  - f_88 * gi_297[k]
                  + f_78 * gi_337[k]
                  + f_79 * gi_342[k]
                  - f_84 * gi_344[k]
                  - f_87 * gi_351[k]
                  + f_89 * gi_353[k]
                  - f_79 * gi_393[k]
                  - f_82 * gi_398[k]
                  + f_85 * gi_400[k]
                  + f_88 * gi_407[k]
                  - f_90 * gi_409[k];
    }

#pragma omp simd aligned(ab_x, gh_4, gh_11, gh_13, gh_67, gh_74, gh_76, gh_109, gh_116, \
                         gh_118, gh_214, gh_221, gh_223, gh_256, gh_263, gh_265, gh_298, \
                         gh_305, gh_307, gi_4, gi_11, gi_13, gi_88, gi_95, gi_97, gi_144, \
                         gi_151, gi_153, gi_284, gi_291, gi_293, gi_340, gi_347, gi_349, \
                         gi_396, gi_403, gi_405 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_120 * ab_x[k] * gh_4[k]
                  + f_120 * ab_x[k] * gh_11[k]
                  - f_121 * ab_x[k] * gh_13[k]
                  + f_121 * ab_x[k] * gh_67[k]
                  + f_121 * ab_x[k] * gh_74[k]
                  - f_124 * ab_x[k] * gh_76[k]
                  - f_122 * ab_x[k] * gh_109[k]
                  - f_122 * ab_x[k] * gh_116[k]
                  + f_125 * ab_x[k] * gh_118[k]
                  + f_120 * ab_x[k] * gh_214[k]
                  + f_120 * ab_x[k] * gh_221[k]
                  - f_121 * ab_x[k] * gh_223[k]
                  - f_122 * ab_x[k] * gh_256[k]
                  - f_122 * ab_x[k] * gh_263[k]
                  + f_125 * ab_x[k] * gh_265[k]
                  + f_123 * ab_x[k] * gh_298[k]
                  + f_123 * ab_x[k] * gh_305[k]
                  - f_126 * ab_x[k] * gh_307[k]
                  - f_120 * gi_4[k]
                  - f_120 * gi_11[k]
                  + f_121 * gi_13[k]
                  - f_121 * gi_88[k]
                  - f_121 * gi_95[k]
                  + f_124 * gi_97[k]
                  + f_122 * gi_144[k]
                  + f_122 * gi_151[k]
                  - f_125 * gi_153[k]
                  - f_120 * gi_284[k]
                  - f_120 * gi_291[k]
                  + f_121 * gi_293[k]
                  + f_122 * gi_340[k]
                  + f_122 * gi_347[k]
                  - f_125 * gi_349[k]
                  - f_123 * gi_396[k]
                  - f_123 * gi_403[k]
                  + f_126 * gi_405[k];
    }

#pragma omp simd aligned(ab_x, gh_1, gh_6, gh_8, gh_15, gh_17, gh_19, gh_64, gh_69, gh_71, \
                         gh_78, gh_80, gh_82, gh_106, gh_111, gh_113, gh_120, gh_122, gh_124, \
                         gh_211, gh_216, gh_218, gh_225, gh_227, gh_229, gh_253, gh_258, \
                         gh_260, gh_267, gh_269, gh_271, gh_295, gh_300, gh_302, gh_309, \
                         gh_311, gh_313, gi_1, gi_6, gi_8, gi_15, gi_17, gi_19, gi_85, gi_90, \
                         gi_92, gi_99, gi_101, gi_103, gi_141, gi_146, gi_148, gi_155, gi_157, \
                         gi_159, gi_281, gi_286, gi_288, gi_295, gi_297, gi_299, gi_337, \
                         gi_342, gi_344, gi_351, gi_353, gi_355, gi_393, gi_398, gi_400, \
                         gi_407, gi_409, gi_411 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -0.234375 * ab_x[k] * gh_1[k]
                  - 0.46875 * ab_x[k] * gh_6[k]
                  + 2.8125 * ab_x[k] * gh_8[k]
                  - 0.234375 * ab_x[k] * gh_15[k]
                  + 2.8125 * ab_x[k] * gh_17[k]
                  - 1.875 * ab_x[k] * gh_19[k]
                  - 0.46875 * ab_x[k] * gh_64[k]
                  - 0.9375 * ab_x[k] * gh_69[k]
                  + 5.625 * ab_x[k] * gh_71[k]
                  - 0.46875 * ab_x[k] * gh_78[k]
                  + 5.625 * ab_x[k] * gh_80[k]
                  - 3.75 * ab_x[k] * gh_82[k]
                  + 2.8125 * ab_x[k] * gh_106[k]
                  + 5.625 * ab_x[k] * gh_111[k]
                  - 33.75 * ab_x[k] * gh_113[k]
                  + 2.8125 * ab_x[k] * gh_120[k]
                  - 33.75 * ab_x[k] * gh_122[k]
                  + 22.5 * ab_x[k] * gh_124[k]
                  - 0.234375 * ab_x[k] * gh_211[k]
                  - 0.46875 * ab_x[k] * gh_216[k]
                  + 2.8125 * ab_x[k] * gh_218[k]
                  - 0.234375 * ab_x[k] * gh_225[k]
                  + 2.8125 * ab_x[k] * gh_227[k]
                  - 1.875 * ab_x[k] * gh_229[k]
                  + 2.8125 * ab_x[k] * gh_253[k]
                  + 5.625 * ab_x[k] * gh_258[k]
                  - 33.75 * ab_x[k] * gh_260[k]
                  + 2.8125 * ab_x[k] * gh_267[k]
                  - 33.75 * ab_x[k] * gh_269[k]
                  + 22.5 * ab_x[k] * gh_271[k]
                  - 1.875 * ab_x[k] * gh_295[k]
                  - 3.75 * ab_x[k] * gh_300[k]
                  + 22.5 * ab_x[k] * gh_302[k]
                  - 1.875 * ab_x[k] * gh_309[k]
                  + 22.5 * ab_x[k] * gh_311[k]
                  - 15.0 * ab_x[k] * gh_313[k]
                  + 0.234375 * gi_1[k]
                  + 0.46875 * gi_6[k]
                  - 2.8125 * gi_8[k]
                  + 0.234375 * gi_15[k]
                  - 2.8125 * gi_17[k]
                  + 1.875 * gi_19[k]
                  + 0.46875 * gi_85[k]
                  + 0.9375 * gi_90[k]
                  - 5.625 * gi_92[k]
                  + 0.46875 * gi_99[k]
                  - 5.625 * gi_101[k]
                  + 3.75 * gi_103[k]
                  - 2.8125 * gi_141[k]
                  - 5.625 * gi_146[k]
                  + 33.75 * gi_148[k]
                  - 2.8125 * gi_155[k]
                  + 33.75 * gi_157[k]
                  - 22.5 * gi_159[k]
                  + 0.234375 * gi_281[k]
                  + 0.46875 * gi_286[k]
                  - 2.8125 * gi_288[k]
                  + 0.234375 * gi_295[k]
                  - 2.8125 * gi_297[k]
                  + 1.875 * gi_299[k]
                  - 2.8125 * gi_337[k]
                  - 5.625 * gi_342[k]
                  + 33.75 * gi_344[k]
                  - 2.8125 * gi_351[k]
                  + 33.75 * gi_353[k]
                  - 22.5 * gi_355[k]
                  + 1.875 * gi_393[k]
                  + 3.75 * gi_398[k]
                  - 22.5 * gi_400[k]
                  + 1.875 * gi_407[k]
                  - 22.5 * gi_409[k]
                  + 15.0 * gi_411[k];
    }

#pragma omp simd aligned(ab_x, gh_2, gh_7, gh_9, gh_16, gh_18, gh_20, gh_65, gh_70, gh_72, \
                         gh_79, gh_81, gh_83, gh_107, gh_112, gh_114, gh_121, gh_123, gh_125, \
                         gh_212, gh_217, gh_219, gh_226, gh_228, gh_230, gh_254, gh_259, \
                         gh_261, gh_268, gh_270, gh_272, gh_296, gh_301, gh_303, gh_310, \
                         gh_312, gh_314, gi_2, gi_7, gi_9, gi_16, gi_18, gi_20, gi_86, gi_91, \
                         gi_93, gi_100, gi_102, gi_104, gi_142, gi_147, gi_149, gi_156, \
                         gi_158, gi_160, gi_282, gi_287, gi_289, gi_296, gi_298, gi_300, \
                         gi_338, gi_343, gi_345, gi_352, gi_354, gi_356, gi_394, gi_399, \
                         gi_401, gi_408, gi_410, gi_412 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_137 * ab_x[k] * gh_2[k]
                  - f_138 * ab_x[k] * gh_7[k]
                  + f_139 * ab_x[k] * gh_9[k]
                  - f_137 * ab_x[k] * gh_16[k]
                  + f_139 * ab_x[k] * gh_18[k]
                  - f_140 * ab_x[k] * gh_20[k]
                  - f_138 * ab_x[k] * gh_65[k]
                  - f_141 * ab_x[k] * gh_70[k]
                  + f_142 * ab_x[k] * gh_72[k]
                  - f_138 * ab_x[k] * gh_79[k]
                  + f_142 * ab_x[k] * gh_81[k]
                  - f_143 * ab_x[k] * gh_83[k]
                  + f_144 * ab_x[k] * gh_107[k]
                  + f_145 * ab_x[k] * gh_112[k]
                  - f_146 * ab_x[k] * gh_114[k]
                  + f_144 * ab_x[k] * gh_121[k]
                  - f_146 * ab_x[k] * gh_123[k]
                  + f_147 * ab_x[k] * gh_125[k]
                  - f_137 * ab_x[k] * gh_212[k]
                  - f_138 * ab_x[k] * gh_217[k]
                  + f_139 * ab_x[k] * gh_219[k]
                  - f_137 * ab_x[k] * gh_226[k]
                  + f_139 * ab_x[k] * gh_228[k]
                  - f_140 * ab_x[k] * gh_230[k]
                  + f_144 * ab_x[k] * gh_254[k]
                  + f_145 * ab_x[k] * gh_259[k]
                  - f_146 * ab_x[k] * gh_261[k]
                  + f_144 * ab_x[k] * gh_268[k]
                  - f_146 * ab_x[k] * gh_270[k]
                  + f_147 * ab_x[k] * gh_272[k]
                  - f_148 * ab_x[k] * gh_296[k]
                  - f_149 * ab_x[k] * gh_301[k]
                  + f_150 * ab_x[k] * gh_303[k]
                  - f_148 * ab_x[k] * gh_310[k]
                  + f_150 * ab_x[k] * gh_312[k]
                  - f_151 * ab_x[k] * gh_314[k]
                  + f_137 * gi_2[k]
                  + f_138 * gi_7[k]
                  - f_139 * gi_9[k]
                  + f_137 * gi_16[k]
                  - f_139 * gi_18[k]
                  + f_140 * gi_20[k]
                  + f_138 * gi_86[k]
                  + f_141 * gi_91[k]
                  - f_142 * gi_93[k]
                  + f_138 * gi_100[k]
                  - f_142 * gi_102[k]
                  + f_143 * gi_104[k]
                  - f_144 * gi_142[k]
                  - f_145 * gi_147[k]
                  + f_146 * gi_149[k]
                  - f_144 * gi_156[k]
                  + f_146 * gi_158[k]
                  - f_147 * gi_160[k]
                  + f_137 * gi_282[k]
                  + f_138 * gi_287[k]
                  - f_139 * gi_289[k]
                  + f_137 * gi_296[k]
                  - f_139 * gi_298[k]
                  + f_140 * gi_300[k]
                  - f_144 * gi_338[k]
                  - f_145 * gi_343[k]
                  + f_146 * gi_345[k]
                  - f_144 * gi_352[k]
                  + f_146 * gi_354[k]
                  - f_147 * gi_356[k]
                  + f_148 * gi_394[k]
                  + f_149 * gi_399[k]
                  - f_150 * gi_401[k]
                  + f_148 * gi_408[k]
                  - f_150 * gi_410[k]
                  + f_151 * gi_412[k];
    }

#pragma omp simd aligned(ab_x, gh_0, gh_3, gh_5, gh_10, gh_12, gh_14, gh_63, gh_66, gh_68, \
                         gh_73, gh_75, gh_77, gh_105, gh_108, gh_110, gh_115, gh_117, gh_119, \
                         gh_210, gh_213, gh_215, gh_220, gh_222, gh_224, gh_252, gh_255, \
                         gh_257, gh_262, gh_264, gh_266, gh_294, gh_297, gh_299, gh_304, \
                         gh_306, gh_308, gi_0, gi_3, gi_5, gi_10, gi_12, gi_14, gi_84, gi_87, \
                         gi_89, gi_94, gi_96, gi_98, gi_140, gi_143, gi_145, gi_150, gi_152, \
                         gi_154, gi_280, gi_283, gi_285, gi_290, gi_292, gi_294, gi_336, \
                         gi_339, gi_341, gi_346, gi_348, gi_350, gi_392, gi_395, gi_397, \
                         gi_402, gi_404, gi_406 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -0.234375 * ab_x[k] * gh_0[k]
                  - 0.46875 * ab_x[k] * gh_3[k]
                  + 2.8125 * ab_x[k] * gh_5[k]
                  - 0.234375 * ab_x[k] * gh_10[k]
                  + 2.8125 * ab_x[k] * gh_12[k]
                  - 1.875 * ab_x[k] * gh_14[k]
                  - 0.46875 * ab_x[k] * gh_63[k]
                  - 0.9375 * ab_x[k] * gh_66[k]
                  + 5.625 * ab_x[k] * gh_68[k]
                  - 0.46875 * ab_x[k] * gh_73[k]
                  + 5.625 * ab_x[k] * gh_75[k]
                  - 3.75 * ab_x[k] * gh_77[k]
                  + 2.8125 * ab_x[k] * gh_105[k]
                  + 5.625 * ab_x[k] * gh_108[k]
                  - 33.75 * ab_x[k] * gh_110[k]
                  + 2.8125 * ab_x[k] * gh_115[k]
                  - 33.75 * ab_x[k] * gh_117[k]
                  + 22.5 * ab_x[k] * gh_119[k]
                  - 0.234375 * ab_x[k] * gh_210[k]
                  - 0.46875 * ab_x[k] * gh_213[k]
                  + 2.8125 * ab_x[k] * gh_215[k]
                  - 0.234375 * ab_x[k] * gh_220[k]
                  + 2.8125 * ab_x[k] * gh_222[k]
                  - 1.875 * ab_x[k] * gh_224[k]
                  + 2.8125 * ab_x[k] * gh_252[k]
                  + 5.625 * ab_x[k] * gh_255[k]
                  - 33.75 * ab_x[k] * gh_257[k]
                  + 2.8125 * ab_x[k] * gh_262[k]
                  - 33.75 * ab_x[k] * gh_264[k]
                  + 22.5 * ab_x[k] * gh_266[k]
                  - 1.875 * ab_x[k] * gh_294[k]
                  - 3.75 * ab_x[k] * gh_297[k]
                  + 22.5 * ab_x[k] * gh_299[k]
                  - 1.875 * ab_x[k] * gh_304[k]
                  + 22.5 * ab_x[k] * gh_306[k]
                  - 15.0 * ab_x[k] * gh_308[k]
                  + 0.234375 * gi_0[k]
                  + 0.46875 * gi_3[k]
                  - 2.8125 * gi_5[k]
                  + 0.234375 * gi_10[k]
                  - 2.8125 * gi_12[k]
                  + 1.875 * gi_14[k]
                  + 0.46875 * gi_84[k]
                  + 0.9375 * gi_87[k]
                  - 5.625 * gi_89[k]
                  + 0.46875 * gi_94[k]
                  - 5.625 * gi_96[k]
                  + 3.75 * gi_98[k]
                  - 2.8125 * gi_140[k]
                  - 5.625 * gi_143[k]
                  + 33.75 * gi_145[k]
                  - 2.8125 * gi_150[k]
                  + 33.75 * gi_152[k]
                  - 22.5 * gi_154[k]
                  + 0.234375 * gi_280[k]
                  + 0.46875 * gi_283[k]
                  - 2.8125 * gi_285[k]
                  + 0.234375 * gi_290[k]
                  - 2.8125 * gi_292[k]
                  + 1.875 * gi_294[k]
                  - 2.8125 * gi_336[k]
                  - 5.625 * gi_339[k]
                  + 33.75 * gi_341[k]
                  - 2.8125 * gi_346[k]
                  + 33.75 * gi_348[k]
                  - 22.5 * gi_350[k]
                  + 1.875 * gi_392[k]
                  + 3.75 * gi_395[k]
                  - 22.5 * gi_397[k]
                  + 1.875 * gi_402[k]
                  - 22.5 * gi_404[k]
                  + 15.0 * gi_406[k];
    }

#pragma omp simd aligned(ab_x, gh_2, gh_9, gh_16, gh_18, gh_65, gh_72, gh_79, gh_81, gh_107, \
                         gh_114, gh_121, gh_123, gh_212, gh_219, gh_226, gh_228, gh_254, \
                         gh_261, gh_268, gh_270, gh_296, gh_303, gh_310, gh_312, gi_2, gi_9, \
                         gi_16, gi_18, gi_86, gi_93, gi_100, gi_102, gi_142, gi_149, gi_156, \
                         gi_158, gi_282, gi_289, gi_296, gi_298, gi_338, gi_345, gi_352, \
                         gi_354, gi_394, gi_401, gi_408, gi_410 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_152 * ab_x[k] * gh_2[k]
                  - f_120 * ab_x[k] * gh_9[k]
                  - f_152 * ab_x[k] * gh_16[k]
                  + f_120 * ab_x[k] * gh_18[k]
                  + f_120 * ab_x[k] * gh_65[k]
                  - f_121 * ab_x[k] * gh_72[k]
                  - f_120 * ab_x[k] * gh_79[k]
                  + f_121 * ab_x[k] * gh_81[k]
                  - f_153 * ab_x[k] * gh_107[k]
                  + f_122 * ab_x[k] * gh_114[k]
                  + f_153 * ab_x[k] * gh_121[k]
                  - f_122 * ab_x[k] * gh_123[k]
                  + f_152 * ab_x[k] * gh_212[k]
                  - f_120 * ab_x[k] * gh_219[k]
                  - f_152 * ab_x[k] * gh_226[k]
                  + f_120 * ab_x[k] * gh_228[k]
                  - f_153 * ab_x[k] * gh_254[k]
                  + f_122 * ab_x[k] * gh_261[k]
                  + f_153 * ab_x[k] * gh_268[k]
                  - f_122 * ab_x[k] * gh_270[k]
                  + f_124 * ab_x[k] * gh_296[k]
                  - f_123 * ab_x[k] * gh_303[k]
                  - f_124 * ab_x[k] * gh_310[k]
                  + f_123 * ab_x[k] * gh_312[k]
                  - f_152 * gi_2[k]
                  + f_120 * gi_9[k]
                  + f_152 * gi_16[k]
                  - f_120 * gi_18[k]
                  - f_120 * gi_86[k]
                  + f_121 * gi_93[k]
                  + f_120 * gi_100[k]
                  - f_121 * gi_102[k]
                  + f_153 * gi_142[k]
                  - f_122 * gi_149[k]
                  - f_153 * gi_156[k]
                  + f_122 * gi_158[k]
                  - f_152 * gi_282[k]
                  + f_120 * gi_289[k]
                  + f_152 * gi_296[k]
                  - f_120 * gi_298[k]
                  + f_153 * gi_338[k]
                  - f_122 * gi_345[k]
                  - f_153 * gi_352[k]
                  + f_122 * gi_354[k]
                  - f_124 * gi_394[k]
                  + f_123 * gi_401[k]
                  + f_124 * gi_408[k]
                  - f_123 * gi_410[k];
    }

#pragma omp simd aligned(ab_x, gh_0, gh_3, gh_5, gh_10, gh_12, gh_63, gh_66, gh_68, gh_73, \
                         gh_75, gh_105, gh_108, gh_110, gh_115, gh_117, gh_210, gh_213, \
                         gh_215, gh_220, gh_222, gh_252, gh_255, gh_257, gh_262, gh_264, \
                         gh_294, gh_297, gh_299, gh_304, gh_306, gi_0, gi_3, gi_5, gi_10, \
                         gi_12, gi_84, gi_87, gi_89, gi_94, gi_96, gi_140, gi_143, gi_145, \
                         gi_150, gi_152, gi_280, gi_283, gi_285, gi_290, gi_292, gi_336, \
                         gi_339, gi_341, gi_346, gi_348, gi_392, gi_395, gi_397, gi_402, \
                         gi_404 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_86 * ab_x[k] * gh_0[k]
                  - f_80 * ab_x[k] * gh_3[k]
                  - f_88 * ab_x[k] * gh_5[k]
                  - f_76 * ab_x[k] * gh_10[k]
                  + f_79 * ab_x[k] * gh_12[k]
                  + f_80 * ab_x[k] * gh_63[k]
                  - f_81 * ab_x[k] * gh_66[k]
                  - f_82 * ab_x[k] * gh_68[k]
                  - f_77 * ab_x[k] * gh_73[k]
                  + f_83 * ab_x[k] * gh_75[k]
                  - f_87 * ab_x[k] * gh_105[k]
                  + f_79 * ab_x[k] * gh_108[k]
                  + f_89 * ab_x[k] * gh_110[k]
                  + f_78 * ab_x[k] * gh_115[k]
                  - f_84 * ab_x[k] * gh_117[k]
                  + f_86 * ab_x[k] * gh_210[k]
                  - f_80 * ab_x[k] * gh_213[k]
                  - f_88 * ab_x[k] * gh_215[k]
                  - f_76 * ab_x[k] * gh_220[k]
                  + f_79 * ab_x[k] * gh_222[k]
                  - f_87 * ab_x[k] * gh_252[k]
                  + f_79 * ab_x[k] * gh_255[k]
                  + f_89 * ab_x[k] * gh_257[k]
                  + f_78 * ab_x[k] * gh_262[k]
                  - f_84 * ab_x[k] * gh_264[k]
                  + f_88 * ab_x[k] * gh_294[k]
                  - f_82 * ab_x[k] * gh_297[k]
                  - f_90 * ab_x[k] * gh_299[k]
                  - f_79 * ab_x[k] * gh_304[k]
                  + f_85 * ab_x[k] * gh_306[k]
                  - f_86 * gi_0[k]
                  + f_80 * gi_3[k]
                  + f_88 * gi_5[k]
                  + f_76 * gi_10[k]
                  - f_79 * gi_12[k]
                  - f_80 * gi_84[k]
                  + f_81 * gi_87[k]
                  + f_82 * gi_89[k]
                  + f_77 * gi_94[k]
                  - f_83 * gi_96[k]
                  + f_87 * gi_140[k]
                  - f_79 * gi_143[k]
                  - f_89 * gi_145[k]
                  - f_78 * gi_150[k]
                  + f_84 * gi_152[k]
                  - f_86 * gi_280[k]
                  + f_80 * gi_283[k]
                  + f_88 * gi_285[k]
                  + f_76 * gi_290[k]
                  - f_79 * gi_292[k]
                  + f_87 * gi_336[k]
                  - f_79 * gi_339[k]
                  - f_89 * gi_341[k]
                  - f_78 * gi_346[k]
                  + f_84 * gi_348[k]
                  - f_88 * gi_392[k]
                  + f_82 * gi_395[k]
                  + f_90 * gi_397[k]
                  + f_79 * gi_402[k]
                  - f_85 * gi_404[k];
    }

#pragma omp simd aligned(ab_x, gh_2, gh_7, gh_16, gh_65, gh_70, gh_79, gh_107, gh_112, gh_121, \
                         gh_212, gh_217, gh_226, gh_254, gh_259, gh_268, gh_296, gh_301, \
                         gh_310, gi_2, gi_7, gi_16, gi_86, gi_91, gi_100, gi_142, gi_147, \
                         gi_156, gi_282, gi_287, gi_296, gi_338, gi_343, gi_352, gi_394, \
                         gi_399, gi_408 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_154 * ab_x[k] * gh_2[k]
                  + f_155 * ab_x[k] * gh_7[k]
                  - f_154 * ab_x[k] * gh_16[k]
                  - f_156 * ab_x[k] * gh_65[k]
                  + f_157 * ab_x[k] * gh_70[k]
                  - f_156 * ab_x[k] * gh_79[k]
                  + f_157 * ab_x[k] * gh_107[k]
                  - f_158 * ab_x[k] * gh_112[k]
                  + f_157 * ab_x[k] * gh_121[k]
                  - f_154 * ab_x[k] * gh_212[k]
                  + f_155 * ab_x[k] * gh_217[k]
                  - f_154 * ab_x[k] * gh_226[k]
                  + f_157 * ab_x[k] * gh_254[k]
                  - f_158 * ab_x[k] * gh_259[k]
                  + f_157 * ab_x[k] * gh_268[k]
                  - f_59 * ab_x[k] * gh_296[k]
                  + f_60 * ab_x[k] * gh_301[k]
                  - f_59 * ab_x[k] * gh_310[k]
                  + f_154 * gi_2[k]
                  - f_155 * gi_7[k]
                  + f_154 * gi_16[k]
                  + f_156 * gi_86[k]
                  - f_157 * gi_91[k]
                  + f_156 * gi_100[k]
                  - f_157 * gi_142[k]
                  + f_158 * gi_147[k]
                  - f_157 * gi_156[k]
                  + f_154 * gi_282[k]
                  - f_155 * gi_287[k]
                  + f_154 * gi_296[k]
                  - f_157 * gi_338[k]
                  + f_158 * gi_343[k]
                  - f_157 * gi_352[k]
                  + f_59 * gi_394[k]
                  - f_60 * gi_399[k]
                  + f_59 * gi_408[k];
    }

#pragma omp simd aligned(ab_x, gh_0, gh_3, gh_10, gh_63, gh_66, gh_73, gh_105, gh_108, gh_115, \
                         gh_210, gh_213, gh_220, gh_252, gh_255, gh_262, gh_294, gh_297, \
                         gh_304, gi_0, gi_3, gi_10, gi_84, gi_87, gi_94, gi_140, gi_143, \
                         gi_150, gi_280, gi_283, gi_290, gi_336, gi_339, gi_346, gi_392, \
                         gi_395, gi_402 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_29 * ab_x[k] * gh_0[k]
                  + f_23 * ab_x[k] * gh_3[k]
                  - f_22 * ab_x[k] * gh_10[k]
                  - f_30 * ab_x[k] * gh_63[k]
                  + f_26 * ab_x[k] * gh_66[k]
                  - f_23 * ab_x[k] * gh_73[k]
                  + f_31 * ab_x[k] * gh_105[k]
                  - f_27 * ab_x[k] * gh_108[k]
                  + f_24 * ab_x[k] * gh_115[k]
                  - f_29 * ab_x[k] * gh_210[k]
                  + f_23 * ab_x[k] * gh_213[k]
                  - f_22 * ab_x[k] * gh_220[k]
                  + f_31 * ab_x[k] * gh_252[k]
                  - f_27 * ab_x[k] * gh_255[k]
                  + f_24 * ab_x[k] * gh_262[k]
                  - f_32 * ab_x[k] * gh_294[k]
                  + f_28 * ab_x[k] * gh_297[k]
                  - f_25 * ab_x[k] * gh_304[k]
                  + f_29 * gi_0[k]
                  - f_23 * gi_3[k]
                  + f_22 * gi_10[k]
                  + f_30 * gi_84[k]
                  - f_26 * gi_87[k]
                  + f_23 * gi_94[k]
                  - f_31 * gi_140[k]
                  + f_27 * gi_143[k]
                  - f_24 * gi_150[k]
                  + f_29 * gi_280[k]
                  - f_23 * gi_283[k]
                  + f_22 * gi_290[k]
                  - f_31 * gi_336[k]
                  + f_27 * gi_339[k]
                  - f_24 * gi_346[k]
                  + f_32 * gi_392[k]
                  - f_28 * gi_395[k]
                  + f_25 * gi_402[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_43, gh_48, gh_57, gh_190, gh_195, gh_204, gh_232, \
                         gh_237, gh_246, gh_274, gh_279, gh_288, gi_57, gi_62, gi_71, gi_253, \
                         gi_258, gi_267, gi_311, gi_318, gi_329, gi_367, gi_374, \
                         gi_385 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = f_43 * ab_x[k] * gh_43[k]
                  - f_17 * ab_x[k] * gh_48[k]
                  + f_44 * ab_x[k] * gh_57[k]
                  - f_17 * ab_x[k] * gh_190[k]
                  + f_18 * ab_x[k] * gh_195[k]
                  - f_20 * ab_x[k] * gh_204[k]
                  - f_43 * ab_y[k] * gh_232[k]
                  + f_17 * ab_y[k] * gh_237[k]
                  - f_44 * ab_y[k] * gh_246[k]
                  + f_17 * ab_y[k] * gh_274[k]
                  - f_18 * ab_y[k] * gh_279[k]
                  + f_20 * ab_y[k] * gh_288[k]
                  - f_43 * gi_57[k]
                  + f_17 * gi_62[k]
                  - f_44 * gi_71[k]
                  + f_17 * gi_253[k]
                  - f_18 * gi_258[k]
                  + f_20 * gi_267[k]
                  + f_43 * gi_311[k]
                  - f_17 * gi_318[k]
                  + f_44 * gi_329[k]
                  - f_17 * gi_367[k]
                  + f_18 * gi_374[k]
                  - f_20 * gi_385[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_46, gh_53, gh_193, gh_200, gh_235, gh_242, gh_277, \
                         gh_284, gi_60, gi_67, gi_256, gi_263, gi_315, gi_324, gi_371, \
                         gi_380 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = f_66 * ab_x[k] * gh_46[k]
                  - f_66 * ab_x[k] * gh_53[k]
                  - f_56 * ab_x[k] * gh_193[k]
                  + f_56 * ab_x[k] * gh_200[k]
                  - f_66 * ab_y[k] * gh_235[k]
                  + f_66 * ab_y[k] * gh_242[k]
                  + f_56 * ab_y[k] * gh_277[k]
                  - f_56 * ab_y[k] * gh_284[k]
                  - f_66 * gi_60[k]
                  + f_66 * gi_67[k]
                  + f_56 * gi_256[k]
                  - f_56 * gi_263[k]
                  + f_66 * gi_315[k]
                  - f_66 * gi_324[k]
                  - f_56 * gi_371[k]
                  + f_56 * gi_380[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_43, gh_48, gh_50, gh_57, gh_59, gh_190, gh_195, \
                         gh_197, gh_204, gh_206, gh_232, gh_237, gh_239, gh_246, gh_248, \
                         gh_274, gh_279, gh_281, gh_288, gh_290, gi_57, gi_62, gi_64, gi_71, \
                         gi_73, gi_253, gi_258, gi_260, gi_267, gi_269, gi_311, gi_318, \
                         gi_320, gi_329, gi_331, gi_367, gi_374, gi_376, gi_385, \
                         gi_387 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = -f_109 * ab_x[k] * gh_43[k]
                  - f_73 * ab_x[k] * gh_48[k]
                  + f_110 * ab_x[k] * gh_50[k]
                  + f_111 * ab_x[k] * gh_57[k]
                  - f_70 * ab_x[k] * gh_59[k]
                  + f_67 * ab_x[k] * gh_190[k]
                  + f_69 * ab_x[k] * gh_195[k]
                  - f_71 * ab_x[k] * gh_197[k]
                  - f_73 * ab_x[k] * gh_204[k]
                  + f_74 * ab_x[k] * gh_206[k]
                  + f_109 * ab_y[k] * gh_232[k]
                  + f_73 * ab_y[k] * gh_237[k]
                  - f_110 * ab_y[k] * gh_239[k]
                  - f_111 * ab_y[k] * gh_246[k]
                  + f_70 * ab_y[k] * gh_248[k]
                  - f_67 * ab_y[k] * gh_274[k]
                  - f_69 * ab_y[k] * gh_279[k]
                  + f_71 * ab_y[k] * gh_281[k]
                  + f_73 * ab_y[k] * gh_288[k]
                  - f_74 * ab_y[k] * gh_290[k]
                  + f_109 * gi_57[k]
                  + f_73 * gi_62[k]
                  - f_110 * gi_64[k]
                  - f_111 * gi_71[k]
                  + f_70 * gi_73[k]
                  - f_67 * gi_253[k]
                  - f_69 * gi_258[k]
                  + f_71 * gi_260[k]
                  + f_73 * gi_267[k]
                  - f_74 * gi_269[k]
                  - f_109 * gi_311[k]
                  - f_73 * gi_318[k]
                  + f_110 * gi_320[k]
                  + f_111 * gi_329[k]
                  - f_70 * gi_331[k]
                  + f_67 * gi_367[k]
                  + f_69 * gi_374[k]
                  - f_71 * gi_376[k]
                  - f_73 * gi_385[k]
                  + f_74 * gi_387[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_46, gh_53, gh_55, gh_193, gh_200, gh_202, gh_235, \
                         gh_242, gh_244, gh_277, gh_284, gh_286, gi_60, gi_67, gi_69, gi_256, \
                         gi_263, gi_265, gi_315, gi_324, gi_326, gi_371, gi_380, \
                         gi_382 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = -13.125 * ab_x[k] * gh_46[k]
                  - 13.125 * ab_x[k] * gh_53[k]
                  + 26.25 * ab_x[k] * gh_55[k]
                  + 26.25 * ab_x[k] * gh_193[k]
                  + 26.25 * ab_x[k] * gh_200[k]
                  - 52.5 * ab_x[k] * gh_202[k]
                  + 13.125 * ab_y[k] * gh_235[k]
                  + 13.125 * ab_y[k] * gh_242[k]
                  - 26.25 * ab_y[k] * gh_244[k]
                  - 26.25 * ab_y[k] * gh_277[k]
                  - 26.25 * ab_y[k] * gh_284[k]
                  + 52.5 * ab_y[k] * gh_286[k]
                  + 13.125 * gi_60[k]
                  + 13.125 * gi_67[k]
                  - 26.25 * gi_69[k]
                  - 26.25 * gi_256[k]
                  - 26.25 * gi_263[k]
                  + 52.5 * gi_265[k]
                  - 13.125 * gi_315[k]
                  - 13.125 * gi_324[k]
                  + 26.25 * gi_326[k]
                  + 26.25 * gi_371[k]
                  + 26.25 * gi_380[k]
                  - 52.5 * gi_382[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_43, gh_48, gh_50, gh_57, gh_59, gh_61, gh_190, gh_195, \
                         gh_197, gh_204, gh_206, gh_208, gh_232, gh_237, gh_239, gh_246, \
                         gh_248, gh_250, gh_274, gh_279, gh_281, gh_288, gh_290, gh_292, \
                         gi_57, gi_62, gi_64, gi_71, gi_73, gi_75, gi_253, gi_258, gi_260, \
                         gi_267, gi_269, gi_271, gi_311, gi_318, gi_320, gi_329, gi_331, \
                         gi_333, gi_367, gi_374, gi_376, gi_385, gi_387, \
                         gi_389 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = f_152 * ab_x[k] * gh_43[k]
                  + f_120 * ab_x[k] * gh_48[k]
                  - f_153 * ab_x[k] * gh_50[k]
                  + f_152 * ab_x[k] * gh_57[k]
                  - f_153 * ab_x[k] * gh_59[k]
                  + f_124 * ab_x[k] * gh_61[k]
                  - f_120 * ab_x[k] * gh_190[k]
                  - f_121 * ab_x[k] * gh_195[k]
                  + f_122 * ab_x[k] * gh_197[k]
                  - f_120 * ab_x[k] * gh_204[k]
                  + f_122 * ab_x[k] * gh_206[k]
                  - f_123 * ab_x[k] * gh_208[k]
                  - f_152 * ab_y[k] * gh_232[k]
                  - f_120 * ab_y[k] * gh_237[k]
                  + f_153 * ab_y[k] * gh_239[k]
                  - f_152 * ab_y[k] * gh_246[k]
                  + f_153 * ab_y[k] * gh_248[k]
                  - f_124 * ab_y[k] * gh_250[k]
                  + f_120 * ab_y[k] * gh_274[k]
                  + f_121 * ab_y[k] * gh_279[k]
                  - f_122 * ab_y[k] * gh_281[k]
                  + f_120 * ab_y[k] * gh_288[k]
                  - f_122 * ab_y[k] * gh_290[k]
                  + f_123 * ab_y[k] * gh_292[k]
                  - f_152 * gi_57[k]
                  - f_120 * gi_62[k]
                  + f_153 * gi_64[k]
                  - f_152 * gi_71[k]
                  + f_153 * gi_73[k]
                  - f_124 * gi_75[k]
                  + f_120 * gi_253[k]
                  + f_121 * gi_258[k]
                  - f_122 * gi_260[k]
                  + f_120 * gi_267[k]
                  - f_122 * gi_269[k]
                  + f_123 * gi_271[k]
                  + f_152 * gi_311[k]
                  + f_120 * gi_318[k]
                  - f_153 * gi_320[k]
                  + f_152 * gi_329[k]
                  - f_153 * gi_331[k]
                  + f_124 * gi_333[k]
                  - f_120 * gi_367[k]
                  - f_121 * gi_374[k]
                  + f_122 * gi_376[k]
                  - f_120 * gi_385[k]
                  + f_122 * gi_387[k]
                  - f_123 * gi_389[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_44, gh_49, gh_51, gh_58, gh_60, gh_62, gh_191, gh_196, \
                         gh_198, gh_205, gh_207, gh_209, gh_233, gh_238, gh_240, gh_247, \
                         gh_249, gh_251, gh_275, gh_280, gh_282, gh_289, gh_291, gh_293, \
                         gi_58, gi_63, gi_65, gi_72, gi_74, gi_76, gi_254, gi_259, gi_261, \
                         gi_268, gi_270, gi_272, gi_312, gi_319, gi_321, gi_330, gi_332, \
                         gi_334, gi_368, gi_375, gi_377, gi_386, gi_388, \
                         gi_390 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = f_159 * ab_x[k] * gh_44[k]
                  + f_127 * ab_x[k] * gh_49[k]
                  - f_160 * ab_x[k] * gh_51[k]
                  + f_159 * ab_x[k] * gh_58[k]
                  - f_160 * ab_x[k] * gh_60[k]
                  + f_161 * ab_x[k] * gh_62[k]
                  - f_127 * ab_x[k] * gh_191[k]
                  - f_128 * ab_x[k] * gh_196[k]
                  + f_129 * ab_x[k] * gh_198[k]
                  - f_127 * ab_x[k] * gh_205[k]
                  + f_129 * ab_x[k] * gh_207[k]
                  - f_130 * ab_x[k] * gh_209[k]
                  - f_159 * ab_y[k] * gh_233[k]
                  - f_127 * ab_y[k] * gh_238[k]
                  + f_160 * ab_y[k] * gh_240[k]
                  - f_159 * ab_y[k] * gh_247[k]
                  + f_160 * ab_y[k] * gh_249[k]
                  - f_161 * ab_y[k] * gh_251[k]
                  + f_127 * ab_y[k] * gh_275[k]
                  + f_128 * ab_y[k] * gh_280[k]
                  - f_129 * ab_y[k] * gh_282[k]
                  + f_127 * ab_y[k] * gh_289[k]
                  - f_129 * ab_y[k] * gh_291[k]
                  + f_130 * ab_y[k] * gh_293[k]
                  - f_159 * gi_58[k]
                  - f_127 * gi_63[k]
                  + f_160 * gi_65[k]
                  - f_159 * gi_72[k]
                  + f_160 * gi_74[k]
                  - f_161 * gi_76[k]
                  + f_127 * gi_254[k]
                  + f_128 * gi_259[k]
                  - f_129 * gi_261[k]
                  + f_127 * gi_268[k]
                  - f_129 * gi_270[k]
                  + f_130 * gi_272[k]
                  + f_159 * gi_312[k]
                  + f_127 * gi_319[k]
                  - f_160 * gi_321[k]
                  + f_159 * gi_330[k]
                  - f_160 * gi_332[k]
                  + f_161 * gi_334[k]
                  - f_127 * gi_368[k]
                  - f_128 * gi_375[k]
                  + f_129 * gi_377[k]
                  - f_127 * gi_386[k]
                  + f_129 * gi_388[k]
                  - f_130 * gi_390[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_42, gh_45, gh_47, gh_52, gh_54, gh_56, gh_189, gh_192, \
                         gh_194, gh_199, gh_201, gh_203, gh_231, gh_234, gh_236, gh_241, \
                         gh_243, gh_245, gh_273, gh_276, gh_278, gh_283, gh_285, gh_287, \
                         gi_56, gi_59, gi_61, gi_66, gi_68, gi_70, gi_252, gi_255, gi_257, \
                         gi_262, gi_264, gi_266, gi_309, gi_314, gi_316, gi_323, gi_325, \
                         gi_327, gi_365, gi_370, gi_372, gi_379, gi_381, \
                         gi_383 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = f_152 * ab_x[k] * gh_42[k]
                  + f_120 * ab_x[k] * gh_45[k]
                  - f_153 * ab_x[k] * gh_47[k]
                  + f_152 * ab_x[k] * gh_52[k]
                  - f_153 * ab_x[k] * gh_54[k]
                  + f_124 * ab_x[k] * gh_56[k]
                  - f_120 * ab_x[k] * gh_189[k]
                  - f_121 * ab_x[k] * gh_192[k]
                  + f_122 * ab_x[k] * gh_194[k]
                  - f_120 * ab_x[k] * gh_199[k]
                  + f_122 * ab_x[k] * gh_201[k]
                  - f_123 * ab_x[k] * gh_203[k]
                  - f_152 * ab_y[k] * gh_231[k]
                  - f_120 * ab_y[k] * gh_234[k]
                  + f_153 * ab_y[k] * gh_236[k]
                  - f_152 * ab_y[k] * gh_241[k]
                  + f_153 * ab_y[k] * gh_243[k]
                  - f_124 * ab_y[k] * gh_245[k]
                  + f_120 * ab_y[k] * gh_273[k]
                  + f_121 * ab_y[k] * gh_276[k]
                  - f_122 * ab_y[k] * gh_278[k]
                  + f_120 * ab_y[k] * gh_283[k]
                  - f_122 * ab_y[k] * gh_285[k]
                  + f_123 * ab_y[k] * gh_287[k]
                  - f_152 * gi_56[k]
                  - f_120 * gi_59[k]
                  + f_153 * gi_61[k]
                  - f_152 * gi_66[k]
                  + f_153 * gi_68[k]
                  - f_124 * gi_70[k]
                  + f_120 * gi_252[k]
                  + f_121 * gi_255[k]
                  - f_122 * gi_257[k]
                  + f_120 * gi_262[k]
                  - f_122 * gi_264[k]
                  + f_123 * gi_266[k]
                  + f_152 * gi_309[k]
                  + f_120 * gi_314[k]
                  - f_153 * gi_316[k]
                  + f_152 * gi_323[k]
                  - f_153 * gi_325[k]
                  + f_124 * gi_327[k]
                  - f_120 * gi_365[k]
                  - f_121 * gi_370[k]
                  + f_122 * gi_372[k]
                  - f_120 * gi_379[k]
                  + f_122 * gi_381[k]
                  - f_123 * gi_383[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_44, gh_51, gh_58, gh_60, gh_191, gh_198, gh_205, \
                         gh_207, gh_233, gh_240, gh_247, gh_249, gh_275, gh_282, gh_289, \
                         gh_291, gi_58, gi_65, gi_72, gi_74, gi_254, gi_261, gi_268, gi_270, \
                         gi_312, gi_321, gi_330, gi_332, gi_368, gi_377, gi_386, \
                         gi_388 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = -6.5625 * ab_x[k] * gh_44[k]
                  + 13.125 * ab_x[k] * gh_51[k]
                  + 6.5625 * ab_x[k] * gh_58[k]
                  - 13.125 * ab_x[k] * gh_60[k]
                  + 13.125 * ab_x[k] * gh_191[k]
                  - 26.25 * ab_x[k] * gh_198[k]
                  - 13.125 * ab_x[k] * gh_205[k]
                  + 26.25 * ab_x[k] * gh_207[k]
                  + 6.5625 * ab_y[k] * gh_233[k]
                  - 13.125 * ab_y[k] * gh_240[k]
                  - 6.5625 * ab_y[k] * gh_247[k]
                  + 13.125 * ab_y[k] * gh_249[k]
                  - 13.125 * ab_y[k] * gh_275[k]
                  + 26.25 * ab_y[k] * gh_282[k]
                  + 13.125 * ab_y[k] * gh_289[k]
                  - 26.25 * ab_y[k] * gh_291[k]
                  + 6.5625 * gi_58[k]
                  - 13.125 * gi_65[k]
                  - 6.5625 * gi_72[k]
                  + 13.125 * gi_74[k]
                  - 13.125 * gi_254[k]
                  + 26.25 * gi_261[k]
                  + 13.125 * gi_268[k]
                  - 26.25 * gi_270[k]
                  - 6.5625 * gi_312[k]
                  + 13.125 * gi_321[k]
                  + 6.5625 * gi_330[k]
                  - 13.125 * gi_332[k]
                  + 13.125 * gi_368[k]
                  - 26.25 * gi_377[k]
                  - 13.125 * gi_386[k]
                  + 26.25 * gi_388[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_42, gh_45, gh_47, gh_52, gh_54, gh_189, gh_192, \
                         gh_194, gh_199, gh_201, gh_231, gh_234, gh_236, gh_241, gh_243, \
                         gh_273, gh_276, gh_278, gh_283, gh_285, gi_56, gi_59, gi_61, gi_66, \
                         gi_68, gi_252, gi_255, gi_257, gi_262, gi_264, gi_309, gi_314, \
                         gi_316, gi_323, gi_325, gi_365, gi_370, gi_372, gi_379, \
                         gi_381 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = -f_111 * ab_x[k] * gh_42[k]
                  + f_73 * ab_x[k] * gh_45[k]
                  + f_70 * ab_x[k] * gh_47[k]
                  + f_109 * ab_x[k] * gh_52[k]
                  - f_110 * ab_x[k] * gh_54[k]
                  + f_73 * ab_x[k] * gh_189[k]
                  - f_69 * ab_x[k] * gh_192[k]
                  - f_74 * ab_x[k] * gh_194[k]
                  - f_67 * ab_x[k] * gh_199[k]
                  + f_71 * ab_x[k] * gh_201[k]
                  + f_111 * ab_y[k] * gh_231[k]
                  - f_73 * ab_y[k] * gh_234[k]
                  - f_70 * ab_y[k] * gh_236[k]
                  - f_109 * ab_y[k] * gh_241[k]
                  + f_110 * ab_y[k] * gh_243[k]
                  - f_73 * ab_y[k] * gh_273[k]
                  + f_69 * ab_y[k] * gh_276[k]
                  + f_74 * ab_y[k] * gh_278[k]
                  + f_67 * ab_y[k] * gh_283[k]
                  - f_71 * ab_y[k] * gh_285[k]
                  + f_111 * gi_56[k]
                  - f_73 * gi_59[k]
                  - f_70 * gi_61[k]
                  - f_109 * gi_66[k]
                  + f_110 * gi_68[k]
                  - f_73 * gi_252[k]
                  + f_69 * gi_255[k]
                  + f_74 * gi_257[k]
                  + f_67 * gi_262[k]
                  - f_71 * gi_264[k]
                  - f_111 * gi_309[k]
                  + f_73 * gi_314[k]
                  + f_70 * gi_316[k]
                  + f_109 * gi_323[k]
                  - f_110 * gi_325[k]
                  + f_73 * gi_365[k]
                  - f_69 * gi_370[k]
                  - f_74 * gi_372[k]
                  - f_67 * gi_379[k]
                  + f_71 * gi_381[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_44, gh_49, gh_58, gh_191, gh_196, gh_205, gh_233, \
                         gh_238, gh_247, gh_275, gh_280, gh_289, gi_58, gi_63, gi_72, gi_254, \
                         gi_259, gi_268, gi_312, gi_319, gi_330, gi_368, gi_375, \
                         gi_386 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = f_170 * ab_x[k] * gh_44[k]
                  - f_171 * ab_x[k] * gh_49[k]
                  + f_170 * ab_x[k] * gh_58[k]
                  - f_134 * ab_x[k] * gh_191[k]
                  + f_135 * ab_x[k] * gh_196[k]
                  - f_134 * ab_x[k] * gh_205[k]
                  - f_170 * ab_y[k] * gh_233[k]
                  + f_171 * ab_y[k] * gh_238[k]
                  - f_170 * ab_y[k] * gh_247[k]
                  + f_134 * ab_y[k] * gh_275[k]
                  - f_135 * ab_y[k] * gh_280[k]
                  + f_134 * ab_y[k] * gh_289[k]
                  - f_170 * gi_58[k]
                  + f_171 * gi_63[k]
                  - f_170 * gi_72[k]
                  + f_134 * gi_254[k]
                  - f_135 * gi_259[k]
                  + f_134 * gi_268[k]
                  + f_170 * gi_312[k]
                  - f_171 * gi_319[k]
                  + f_170 * gi_330[k]
                  - f_134 * gi_368[k]
                  + f_135 * gi_375[k]
                  - f_134 * gi_386[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_42, gh_45, gh_52, gh_189, gh_192, gh_199, gh_231, \
                         gh_234, gh_241, gh_273, gh_276, gh_283, gi_56, gi_59, gi_66, gi_252, \
                         gi_255, gi_262, gi_309, gi_314, gi_323, gi_365, gi_370, \
                         gi_379 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = f_44 * ab_x[k] * gh_42[k]
                  - f_17 * ab_x[k] * gh_45[k]
                  + f_43 * ab_x[k] * gh_52[k]
                  - f_20 * ab_x[k] * gh_189[k]
                  + f_18 * ab_x[k] * gh_192[k]
                  - f_17 * ab_x[k] * gh_199[k]
                  - f_44 * ab_y[k] * gh_231[k]
                  + f_17 * ab_y[k] * gh_234[k]
                  - f_43 * ab_y[k] * gh_241[k]
                  + f_20 * ab_y[k] * gh_273[k]
                  - f_18 * ab_y[k] * gh_276[k]
                  + f_17 * ab_y[k] * gh_283[k]
                  - f_44 * gi_56[k]
                  + f_17 * gi_59[k]
                  - f_43 * gi_66[k]
                  + f_20 * gi_252[k]
                  - f_18 * gi_255[k]
                  + f_17 * gi_262[k]
                  + f_44 * gi_309[k]
                  - f_17 * gi_314[k]
                  + f_43 * gi_323[k]
                  - f_20 * gi_365[k]
                  + f_18 * gi_370[k]
                  - f_17 * gi_379[k];
    }

#pragma omp simd aligned(ab_x, gh_1, gh_6, gh_15, gh_64, gh_69, gh_78, gh_106, gh_111, gh_120, \
                         gh_211, gh_216, gh_225, gh_253, gh_258, gh_267, gi_1, gi_6, gi_15, \
                         gi_85, gi_90, gi_99, gi_141, gi_146, gi_155, gi_281, gi_286, gi_295, \
                         gi_337, gi_342, gi_351 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = f_6 * ab_x[k] * gh_1[k]
                  - f_4 * ab_x[k] * gh_6[k]
                  + f_15 * ab_x[k] * gh_15[k]
                  - f_4 * ab_x[k] * gh_64[k]
                  + f_9 * ab_x[k] * gh_69[k]
                  - f_13 * ab_x[k] * gh_78[k]
                  - f_7 * ab_x[k] * gh_106[k]
                  + f_11 * ab_x[k] * gh_111[k]
                  - f_16 * ab_x[k] * gh_120[k]
                  - f_3 * ab_x[k] * gh_211[k]
                  + f_8 * ab_x[k] * gh_216[k]
                  - f_12 * ab_x[k] * gh_225[k]
                  + f_5 * ab_x[k] * gh_253[k]
                  - f_10 * ab_x[k] * gh_258[k]
                  + f_14 * ab_x[k] * gh_267[k]
                  - f_6 * gi_1[k]
                  + f_4 * gi_6[k]
                  - f_15 * gi_15[k]
                  + f_4 * gi_85[k]
                  - f_9 * gi_90[k]
                  + f_13 * gi_99[k]
                  + f_7 * gi_141[k]
                  - f_11 * gi_146[k]
                  + f_16 * gi_155[k]
                  + f_3 * gi_281[k]
                  - f_8 * gi_286[k]
                  + f_12 * gi_295[k]
                  - f_5 * gi_337[k]
                  + f_10 * gi_342[k]
                  - f_14 * gi_351[k];
    }

#pragma omp simd aligned(ab_x, gh_4, gh_11, gh_67, gh_74, gh_109, gh_116, gh_214, gh_221, \
                         gh_256, gh_263, gi_4, gi_11, gi_88, gi_95, gi_144, gi_151, gi_284, \
                         gi_291, gi_340, gi_347 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = f_54 * ab_x[k] * gh_4[k]
                  - f_54 * ab_x[k] * gh_11[k]
                  - f_52 * ab_x[k] * gh_67[k]
                  + f_52 * ab_x[k] * gh_74[k]
                  - f_55 * ab_x[k] * gh_109[k]
                  + f_55 * ab_x[k] * gh_116[k]
                  - f_51 * ab_x[k] * gh_214[k]
                  + f_51 * ab_x[k] * gh_221[k]
                  + f_53 * ab_x[k] * gh_256[k]
                  - f_53 * ab_x[k] * gh_263[k]
                  - f_54 * gi_4[k]
                  + f_54 * gi_11[k]
                  + f_52 * gi_88[k]
                  - f_52 * gi_95[k]
                  + f_55 * gi_144[k]
                  - f_55 * gi_151[k]
                  + f_51 * gi_284[k]
                  - f_51 * gi_291[k]
                  - f_53 * gi_340[k]
                  + f_53 * gi_347[k];
    }

#pragma omp simd aligned(ab_x, gh_1, gh_6, gh_8, gh_15, gh_17, gh_64, gh_69, gh_71, gh_78, \
                         gh_80, gh_106, gh_111, gh_113, gh_120, gh_122, gh_211, gh_216, \
                         gh_218, gh_225, gh_227, gh_253, gh_258, gh_260, gh_267, gh_269, gi_1, \
                         gi_6, gi_8, gi_15, gi_17, gi_85, gi_90, gi_92, gi_99, gi_101, gi_141, \
                         gi_146, gi_148, gi_155, gi_157, gi_281, gi_286, gi_288, gi_295, \
                         gi_297, gi_337, gi_342, gi_344, gi_351, \
                         gi_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = -0.8203125 * ab_x[k] * gh_1[k]
                  - 0.546875 * ab_x[k] * gh_6[k]
                  + 6.5625 * ab_x[k] * gh_8[k]
                  + 0.2734375 * ab_x[k] * gh_15[k]
                  - 2.1875 * ab_x[k] * gh_17[k]
                  + 1.640625 * ab_x[k] * gh_64[k]
                  + 1.09375 * ab_x[k] * gh_69[k]
                  - 13.125 * ab_x[k] * gh_71[k]
                  - 0.546875 * ab_x[k] * gh_78[k]
                  + 4.375 * ab_x[k] * gh_80[k]
                  + 6.5625 * ab_x[k] * gh_106[k]
                  + 4.375 * ab_x[k] * gh_111[k]
                  - 52.5 * ab_x[k] * gh_113[k]
                  - 2.1875 * ab_x[k] * gh_120[k]
                  + 17.5 * ab_x[k] * gh_122[k]
                  + 2.4609375 * ab_x[k] * gh_211[k]
                  + 1.640625 * ab_x[k] * gh_216[k]
                  - 19.6875 * ab_x[k] * gh_218[k]
                  - 0.8203125 * ab_x[k] * gh_225[k]
                  + 6.5625 * ab_x[k] * gh_227[k]
                  - 19.6875 * ab_x[k] * gh_253[k]
                  - 13.125 * ab_x[k] * gh_258[k]
                  + 157.5 * ab_x[k] * gh_260[k]
                  + 6.5625 * ab_x[k] * gh_267[k]
                  - 52.5 * ab_x[k] * gh_269[k]
                  + 0.8203125 * gi_1[k]
                  + 0.546875 * gi_6[k]
                  - 6.5625 * gi_8[k]
                  - 0.2734375 * gi_15[k]
                  + 2.1875 * gi_17[k]
                  - 1.640625 * gi_85[k]
                  - 1.09375 * gi_90[k]
                  + 13.125 * gi_92[k]
                  + 0.546875 * gi_99[k]
                  - 4.375 * gi_101[k]
                  - 6.5625 * gi_141[k]
                  - 4.375 * gi_146[k]
                  + 52.5 * gi_148[k]
                  + 2.1875 * gi_155[k]
                  - 17.5 * gi_157[k]
                  - 2.4609375 * gi_281[k]
                  - 1.640625 * gi_286[k]
                  + 19.6875 * gi_288[k]
                  + 0.8203125 * gi_295[k]
                  - 6.5625 * gi_297[k]
                  + 19.6875 * gi_337[k]
                  + 13.125 * gi_342[k]
                  - 157.5 * gi_344[k]
                  - 6.5625 * gi_351[k]
                  + 52.5 * gi_353[k];
    }

#pragma omp simd aligned(ab_x, gh_4, gh_11, gh_13, gh_67, gh_74, gh_76, gh_109, gh_116, \
                         gh_118, gh_214, gh_221, gh_223, gh_256, gh_263, gh_265, gi_4, gi_11, \
                         gi_13, gi_88, gi_95, gi_97, gi_144, gi_151, gi_153, gi_284, gi_291, \
                         gi_293, gi_340, gi_347, gi_349 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = -f_73 * ab_x[k] * gh_4[k]
                  - f_73 * ab_x[k] * gh_11[k]
                  + f_69 * ab_x[k] * gh_13[k]
                  + f_69 * ab_x[k] * gh_67[k]
                  + f_69 * ab_x[k] * gh_74[k]
                  - f_70 * ab_x[k] * gh_76[k]
                  + f_74 * ab_x[k] * gh_109[k]
                  + f_74 * ab_x[k] * gh_116[k]
                  - f_75 * ab_x[k] * gh_118[k]
                  + f_67 * ab_x[k] * gh_214[k]
                  + f_67 * ab_x[k] * gh_221[k]
                  - f_68 * ab_x[k] * gh_223[k]
                  - f_71 * ab_x[k] * gh_256[k]
                  - f_71 * ab_x[k] * gh_263[k]
                  + f_72 * ab_x[k] * gh_265[k]
                  + f_73 * gi_4[k]
                  + f_73 * gi_11[k]
                  - f_69 * gi_13[k]
                  - f_69 * gi_88[k]
                  - f_69 * gi_95[k]
                  + f_70 * gi_97[k]
                  - f_74 * gi_144[k]
                  - f_74 * gi_151[k]
                  + f_75 * gi_153[k]
                  - f_67 * gi_284[k]
                  - f_67 * gi_291[k]
                  + f_68 * gi_293[k]
                  + f_71 * gi_340[k]
                  + f_71 * gi_347[k]
                  - f_72 * gi_349[k];
    }

#pragma omp simd aligned(ab_x, gh_1, gh_6, gh_8, gh_15, gh_17, gh_19, gh_64, gh_69, gh_71, \
                         gh_78, gh_80, gh_82, gh_106, gh_111, gh_113, gh_120, gh_122, gh_124, \
                         gh_211, gh_216, gh_218, gh_225, gh_227, gh_229, gh_253, gh_258, \
                         gh_260, gh_267, gh_269, gh_271, gi_1, gi_6, gi_8, gi_15, gi_17, \
                         gi_19, gi_85, gi_90, gi_92, gi_99, gi_101, gi_103, gi_141, gi_146, \
                         gi_148, gi_155, gi_157, gi_159, gi_281, gi_286, gi_288, gi_295, \
                         gi_297, gi_299, gi_337, gi_342, gi_344, gi_351, gi_353, \
                         gi_355 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = f_86 * ab_x[k] * gh_1[k]
                  + f_80 * ab_x[k] * gh_6[k]
                  - f_87 * ab_x[k] * gh_8[k]
                  + f_86 * ab_x[k] * gh_15[k]
                  - f_87 * ab_x[k] * gh_17[k]
                  + f_88 * ab_x[k] * gh_19[k]
                  - f_80 * ab_x[k] * gh_64[k]
                  - f_81 * ab_x[k] * gh_69[k]
                  + f_79 * ab_x[k] * gh_71[k]
                  - f_80 * ab_x[k] * gh_78[k]
                  + f_79 * ab_x[k] * gh_80[k]
                  - f_82 * ab_x[k] * gh_82[k]
                  - f_88 * ab_x[k] * gh_106[k]
                  - f_82 * ab_x[k] * gh_111[k]
                  + f_89 * ab_x[k] * gh_113[k]
                  - f_88 * ab_x[k] * gh_120[k]
                  + f_89 * ab_x[k] * gh_122[k]
                  - f_90 * ab_x[k] * gh_124[k]
                  - f_76 * ab_x[k] * gh_211[k]
                  - f_77 * ab_x[k] * gh_216[k]
                  + f_78 * ab_x[k] * gh_218[k]
                  - f_76 * ab_x[k] * gh_225[k]
                  + f_78 * ab_x[k] * gh_227[k]
                  - f_79 * ab_x[k] * gh_229[k]
                  + f_79 * ab_x[k] * gh_253[k]
                  + f_83 * ab_x[k] * gh_258[k]
                  - f_84 * ab_x[k] * gh_260[k]
                  + f_79 * ab_x[k] * gh_267[k]
                  - f_84 * ab_x[k] * gh_269[k]
                  + f_85 * ab_x[k] * gh_271[k]
                  - f_86 * gi_1[k]
                  - f_80 * gi_6[k]
                  + f_87 * gi_8[k]
                  - f_86 * gi_15[k]
                  + f_87 * gi_17[k]
                  - f_88 * gi_19[k]
                  + f_80 * gi_85[k]
                  + f_81 * gi_90[k]
                  - f_79 * gi_92[k]
                  + f_80 * gi_99[k]
                  - f_79 * gi_101[k]
                  + f_82 * gi_103[k]
                  + f_88 * gi_141[k]
                  + f_82 * gi_146[k]
                  - f_89 * gi_148[k]
                  + f_88 * gi_155[k]
                  - f_89 * gi_157[k]
                  + f_90 * gi_159[k]
                  + f_76 * gi_281[k]
                  + f_77 * gi_286[k]
                  - f_78 * gi_288[k]
                  + f_76 * gi_295[k]
                  - f_78 * gi_297[k]
                  + f_79 * gi_299[k]
                  - f_79 * gi_337[k]
                  - f_83 * gi_342[k]
                  + f_84 * gi_344[k]
                  - f_79 * gi_351[k]
                  + f_84 * gi_353[k]
                  - f_85 * gi_355[k];
    }

#pragma omp simd aligned(ab_x, gh_2, gh_7, gh_9, gh_16, gh_18, gh_20, gh_65, gh_70, gh_72, \
                         gh_79, gh_81, gh_83, gh_107, gh_112, gh_114, gh_121, gh_123, gh_125, \
                         gh_212, gh_217, gh_219, gh_226, gh_228, gh_230, gh_254, gh_259, \
                         gh_261, gh_268, gh_270, gh_272, gi_2, gi_7, gi_9, gi_16, gi_18, \
                         gi_20, gi_86, gi_91, gi_93, gi_100, gi_102, gi_104, gi_142, gi_147, \
                         gi_149, gi_156, gi_158, gi_160, gi_282, gi_287, gi_289, gi_296, \
                         gi_298, gi_300, gi_338, gi_343, gi_345, gi_352, gi_354, \
                         gi_356 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = f_103 * ab_x[k] * gh_2[k]
                  + f_95 * ab_x[k] * gh_7[k]
                  - f_104 * ab_x[k] * gh_9[k]
                  + f_103 * ab_x[k] * gh_16[k]
                  - f_104 * ab_x[k] * gh_18[k]
                  + f_105 * ab_x[k] * gh_20[k]
                  - f_95 * ab_x[k] * gh_65[k]
                  - f_96 * ab_x[k] * gh_70[k]
                  + f_97 * ab_x[k] * gh_72[k]
                  - f_95 * ab_x[k] * gh_79[k]
                  + f_97 * ab_x[k] * gh_81[k]
                  - f_98 * ab_x[k] * gh_83[k]
                  - f_93 * ab_x[k] * gh_107[k]
                  - f_106 * ab_x[k] * gh_112[k]
                  + f_107 * ab_x[k] * gh_114[k]
                  - f_93 * ab_x[k] * gh_121[k]
                  + f_107 * ab_x[k] * gh_123[k]
                  - f_108 * ab_x[k] * gh_125[k]
                  - f_91 * ab_x[k] * gh_212[k]
                  - f_92 * ab_x[k] * gh_217[k]
                  + f_93 * ab_x[k] * gh_219[k]
                  - f_91 * ab_x[k] * gh_226[k]
                  + f_93 * ab_x[k] * gh_228[k]
                  - f_94 * ab_x[k] * gh_230[k]
                  + f_99 * ab_x[k] * gh_254[k]
                  + f_100 * ab_x[k] * gh_259[k]
                  - f_101 * ab_x[k] * gh_261[k]
                  + f_99 * ab_x[k] * gh_268[k]
                  - f_101 * ab_x[k] * gh_270[k]
                  + f_102 * ab_x[k] * gh_272[k]
                  - f_103 * gi_2[k]
                  - f_95 * gi_7[k]
                  + f_104 * gi_9[k]
                  - f_103 * gi_16[k]
                  + f_104 * gi_18[k]
                  - f_105 * gi_20[k]
                  + f_95 * gi_86[k]
                  + f_96 * gi_91[k]
                  - f_97 * gi_93[k]
                  + f_95 * gi_100[k]
                  - f_97 * gi_102[k]
                  + f_98 * gi_104[k]
                  + f_93 * gi_142[k]
                  + f_106 * gi_147[k]
                  - f_107 * gi_149[k]
                  + f_93 * gi_156[k]
                  - f_107 * gi_158[k]
                  + f_108 * gi_160[k]
                  + f_91 * gi_282[k]
                  + f_92 * gi_287[k]
                  - f_93 * gi_289[k]
                  + f_91 * gi_296[k]
                  - f_93 * gi_298[k]
                  + f_94 * gi_300[k]
                  - f_99 * gi_338[k]
                  - f_100 * gi_343[k]
                  + f_101 * gi_345[k]
                  - f_99 * gi_352[k]
                  + f_101 * gi_354[k]
                  - f_102 * gi_356[k];
    }

#pragma omp simd aligned(ab_x, gh_0, gh_3, gh_5, gh_10, gh_12, gh_14, gh_63, gh_66, gh_68, \
                         gh_73, gh_75, gh_77, gh_105, gh_108, gh_110, gh_115, gh_117, gh_119, \
                         gh_210, gh_213, gh_215, gh_220, gh_222, gh_224, gh_252, gh_255, \
                         gh_257, gh_262, gh_264, gh_266, gi_0, gi_3, gi_5, gi_10, gi_12, \
                         gi_14, gi_84, gi_87, gi_89, gi_94, gi_96, gi_98, gi_140, gi_143, \
                         gi_145, gi_150, gi_152, gi_154, gi_280, gi_283, gi_285, gi_290, \
                         gi_292, gi_294, gi_336, gi_339, gi_341, gi_346, gi_348, \
                         gi_350 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = f_86 * ab_x[k] * gh_0[k]
                  + f_80 * ab_x[k] * gh_3[k]
                  - f_87 * ab_x[k] * gh_5[k]
                  + f_86 * ab_x[k] * gh_10[k]
                  - f_87 * ab_x[k] * gh_12[k]
                  + f_88 * ab_x[k] * gh_14[k]
                  - f_80 * ab_x[k] * gh_63[k]
                  - f_81 * ab_x[k] * gh_66[k]
                  + f_79 * ab_x[k] * gh_68[k]
                  - f_80 * ab_x[k] * gh_73[k]
                  + f_79 * ab_x[k] * gh_75[k]
                  - f_82 * ab_x[k] * gh_77[k]
                  - f_88 * ab_x[k] * gh_105[k]
                  - f_82 * ab_x[k] * gh_108[k]
                  + f_89 * ab_x[k] * gh_110[k]
                  - f_88 * ab_x[k] * gh_115[k]
                  + f_89 * ab_x[k] * gh_117[k]
                  - f_90 * ab_x[k] * gh_119[k]
                  - f_76 * ab_x[k] * gh_210[k]
                  - f_77 * ab_x[k] * gh_213[k]
                  + f_78 * ab_x[k] * gh_215[k]
                  - f_76 * ab_x[k] * gh_220[k]
                  + f_78 * ab_x[k] * gh_222[k]
                  - f_79 * ab_x[k] * gh_224[k]
                  + f_79 * ab_x[k] * gh_252[k]
                  + f_83 * ab_x[k] * gh_255[k]
                  - f_84 * ab_x[k] * gh_257[k]
                  + f_79 * ab_x[k] * gh_262[k]
                  - f_84 * ab_x[k] * gh_264[k]
                  + f_85 * ab_x[k] * gh_266[k]
                  - f_86 * gi_0[k]
                  - f_80 * gi_3[k]
                  + f_87 * gi_5[k]
                  - f_86 * gi_10[k]
                  + f_87 * gi_12[k]
                  - f_88 * gi_14[k]
                  + f_80 * gi_84[k]
                  + f_81 * gi_87[k]
                  - f_79 * gi_89[k]
                  + f_80 * gi_94[k]
                  - f_79 * gi_96[k]
                  + f_82 * gi_98[k]
                  + f_88 * gi_140[k]
                  + f_82 * gi_143[k]
                  - f_89 * gi_145[k]
                  + f_88 * gi_150[k]
                  - f_89 * gi_152[k]
                  + f_90 * gi_154[k]
                  + f_76 * gi_280[k]
                  + f_77 * gi_283[k]
                  - f_78 * gi_285[k]
                  + f_76 * gi_290[k]
                  - f_78 * gi_292[k]
                  + f_79 * gi_294[k]
                  - f_79 * gi_336[k]
                  - f_83 * gi_339[k]
                  + f_84 * gi_341[k]
                  - f_79 * gi_346[k]
                  + f_84 * gi_348[k]
                  - f_85 * gi_350[k];
    }

#pragma omp simd aligned(ab_x, gh_2, gh_9, gh_16, gh_18, gh_65, gh_72, gh_79, gh_81, gh_107, \
                         gh_114, gh_121, gh_123, gh_212, gh_219, gh_226, gh_228, gh_254, \
                         gh_261, gh_268, gh_270, gi_2, gi_9, gi_16, gi_18, gi_86, gi_93, \
                         gi_100, gi_102, gi_142, gi_149, gi_156, gi_158, gi_282, gi_289, \
                         gi_296, gi_298, gi_338, gi_345, gi_352, \
                         gi_354 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = -f_111 * ab_x[k] * gh_2[k]
                  + f_73 * ab_x[k] * gh_9[k]
                  + f_111 * ab_x[k] * gh_16[k]
                  - f_73 * ab_x[k] * gh_18[k]
                  + f_73 * ab_x[k] * gh_65[k]
                  - f_69 * ab_x[k] * gh_72[k]
                  - f_73 * ab_x[k] * gh_79[k]
                  + f_69 * ab_x[k] * gh_81[k]
                  + f_70 * ab_x[k] * gh_107[k]
                  - f_74 * ab_x[k] * gh_114[k]
                  - f_70 * ab_x[k] * gh_121[k]
                  + f_74 * ab_x[k] * gh_123[k]
                  + f_109 * ab_x[k] * gh_212[k]
                  - f_67 * ab_x[k] * gh_219[k]
                  - f_109 * ab_x[k] * gh_226[k]
                  + f_67 * ab_x[k] * gh_228[k]
                  - f_110 * ab_x[k] * gh_254[k]
                  + f_71 * ab_x[k] * gh_261[k]
                  + f_110 * ab_x[k] * gh_268[k]
                  - f_71 * ab_x[k] * gh_270[k]
                  + f_111 * gi_2[k]
                  - f_73 * gi_9[k]
                  - f_111 * gi_16[k]
                  + f_73 * gi_18[k]
                  - f_73 * gi_86[k]
                  + f_69 * gi_93[k]
                  + f_73 * gi_100[k]
                  - f_69 * gi_102[k]
                  - f_70 * gi_142[k]
                  + f_74 * gi_149[k]
                  + f_70 * gi_156[k]
                  - f_74 * gi_158[k]
                  - f_109 * gi_282[k]
                  + f_67 * gi_289[k]
                  + f_109 * gi_296[k]
                  - f_67 * gi_298[k]
                  + f_110 * gi_338[k]
                  - f_71 * gi_345[k]
                  - f_110 * gi_352[k]
                  + f_71 * gi_354[k];
    }

#pragma omp simd aligned(ab_x, gh_0, gh_3, gh_5, gh_10, gh_12, gh_63, gh_66, gh_68, gh_73, \
                         gh_75, gh_105, gh_108, gh_110, gh_115, gh_117, gh_210, gh_213, \
                         gh_215, gh_220, gh_222, gh_252, gh_255, gh_257, gh_262, gh_264, gi_0, \
                         gi_3, gi_5, gi_10, gi_12, gi_84, gi_87, gi_89, gi_94, gi_96, gi_140, \
                         gi_143, gi_145, gi_150, gi_152, gi_280, gi_283, gi_285, gi_290, \
                         gi_292, gi_336, gi_339, gi_341, gi_346, \
                         gi_348 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = -0.2734375 * ab_x[k] * gh_0[k]
                  + 0.546875 * ab_x[k] * gh_3[k]
                  + 2.1875 * ab_x[k] * gh_5[k]
                  + 0.8203125 * ab_x[k] * gh_10[k]
                  - 6.5625 * ab_x[k] * gh_12[k]
                  + 0.546875 * ab_x[k] * gh_63[k]
                  - 1.09375 * ab_x[k] * gh_66[k]
                  - 4.375 * ab_x[k] * gh_68[k]
                  - 1.640625 * ab_x[k] * gh_73[k]
                  + 13.125 * ab_x[k] * gh_75[k]
                  + 2.1875 * ab_x[k] * gh_105[k]
                  - 4.375 * ab_x[k] * gh_108[k]
                  - 17.5 * ab_x[k] * gh_110[k]
                  - 6.5625 * ab_x[k] * gh_115[k]
                  + 52.5 * ab_x[k] * gh_117[k]
                  + 0.8203125 * ab_x[k] * gh_210[k]
                  - 1.640625 * ab_x[k] * gh_213[k]
                  - 6.5625 * ab_x[k] * gh_215[k]
                  - 2.4609375 * ab_x[k] * gh_220[k]
                  + 19.6875 * ab_x[k] * gh_222[k]
                  - 6.5625 * ab_x[k] * gh_252[k]
                  + 13.125 * ab_x[k] * gh_255[k]
                  + 52.5 * ab_x[k] * gh_257[k]
                  + 19.6875 * ab_x[k] * gh_262[k]
                  - 157.5 * ab_x[k] * gh_264[k]
                  + 0.2734375 * gi_0[k]
                  - 0.546875 * gi_3[k]
                  - 2.1875 * gi_5[k]
                  - 0.8203125 * gi_10[k]
                  + 6.5625 * gi_12[k]
                  - 0.546875 * gi_84[k]
                  + 1.09375 * gi_87[k]
                  + 4.375 * gi_89[k]
                  + 1.640625 * gi_94[k]
                  - 13.125 * gi_96[k]
                  - 2.1875 * gi_140[k]
                  + 4.375 * gi_143[k]
                  + 17.5 * gi_145[k]
                  + 6.5625 * gi_150[k]
                  - 52.5 * gi_152[k]
                  - 0.8203125 * gi_280[k]
                  + 1.640625 * gi_283[k]
                  + 6.5625 * gi_285[k]
                  + 2.4609375 * gi_290[k]
                  - 19.6875 * gi_292[k]
                  + 6.5625 * gi_336[k]
                  - 13.125 * gi_339[k]
                  - 52.5 * gi_341[k]
                  - 19.6875 * gi_346[k]
                  + 157.5 * gi_348[k];
    }

#pragma omp simd aligned(ab_x, gh_2, gh_7, gh_16, gh_65, gh_70, gh_79, gh_107, gh_112, gh_121, \
                         gh_212, gh_217, gh_226, gh_254, gh_259, gh_268, gi_2, gi_7, gi_16, \
                         gi_86, gi_91, gi_100, gi_142, gi_147, gi_156, gi_282, gi_287, gi_296, \
                         gi_338, gi_343, gi_352 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = f_117 * ab_x[k] * gh_2[k]
                  - f_118 * ab_x[k] * gh_7[k]
                  + f_117 * ab_x[k] * gh_16[k]
                  - f_114 * ab_x[k] * gh_65[k]
                  + f_51 * ab_x[k] * gh_70[k]
                  - f_114 * ab_x[k] * gh_79[k]
                  - f_52 * ab_x[k] * gh_107[k]
                  + f_119 * ab_x[k] * gh_112[k]
                  - f_52 * ab_x[k] * gh_121[k]
                  - f_112 * ab_x[k] * gh_212[k]
                  + f_113 * ab_x[k] * gh_217[k]
                  - f_112 * ab_x[k] * gh_226[k]
                  + f_115 * ab_x[k] * gh_254[k]
                  - f_116 * ab_x[k] * gh_259[k]
                  + f_115 * ab_x[k] * gh_268[k]
                  - f_117 * gi_2[k]
                  + f_118 * gi_7[k]
                  - f_117 * gi_16[k]
                  + f_114 * gi_86[k]
                  - f_51 * gi_91[k]
                  + f_114 * gi_100[k]
                  + f_52 * gi_142[k]
                  - f_119 * gi_147[k]
                  + f_52 * gi_156[k]
                  + f_112 * gi_282[k]
                  - f_113 * gi_287[k]
                  + f_112 * gi_296[k]
                  - f_115 * gi_338[k]
                  + f_116 * gi_343[k]
                  - f_115 * gi_352[k];
    }

#pragma omp simd aligned(ab_x, gh_0, gh_3, gh_10, gh_63, gh_66, gh_73, gh_105, gh_108, gh_115, \
                         gh_210, gh_213, gh_220, gh_252, gh_255, gh_262, gi_0, gi_3, gi_10, \
                         gi_84, gi_87, gi_94, gi_140, gi_143, gi_150, gi_280, gi_283, gi_290, \
                         gi_336, gi_339, gi_346 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = f_15 * ab_x[k] * gh_0[k]
                  - f_4 * ab_x[k] * gh_3[k]
                  + f_6 * ab_x[k] * gh_10[k]
                  - f_13 * ab_x[k] * gh_63[k]
                  + f_9 * ab_x[k] * gh_66[k]
                  - f_4 * ab_x[k] * gh_73[k]
                  - f_16 * ab_x[k] * gh_105[k]
                  + f_11 * ab_x[k] * gh_108[k]
                  - f_7 * ab_x[k] * gh_115[k]
                  - f_12 * ab_x[k] * gh_210[k]
                  + f_8 * ab_x[k] * gh_213[k]
                  - f_3 * ab_x[k] * gh_220[k]
                  + f_14 * ab_x[k] * gh_252[k]
                  - f_10 * ab_x[k] * gh_255[k]
                  + f_5 * ab_x[k] * gh_262[k]
                  - f_15 * gi_0[k]
                  + f_4 * gi_3[k]
                  - f_6 * gi_10[k]
                  + f_13 * gi_84[k]
                  - f_9 * gi_87[k]
                  + f_4 * gi_94[k]
                  + f_16 * gi_140[k]
                  - f_11 * gi_143[k]
                  + f_7 * gi_150[k]
                  + f_12 * gi_280[k]
                  - f_8 * gi_283[k]
                  + f_3 * gi_290[k]
                  - f_14 * gi_336[k]
                  + f_10 * gi_339[k]
                  - f_5 * gi_346[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_43, gh_48, gh_57, gh_148, gh_153, gh_162, gh_232, \
                         gh_237, gh_246, gi_57, gi_62, gi_71, gi_197, gi_202, gi_211, gi_311, \
                         gi_318, gi_329 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = -f_45 * ab_x[k] * gh_43[k]
                  + f_47 * ab_x[k] * gh_48[k]
                  - f_49 * ab_x[k] * gh_57[k]
                  + f_46 * ab_x[k] * gh_148[k]
                  - f_48 * ab_x[k] * gh_153[k]
                  + f_50 * ab_x[k] * gh_162[k]
                  - f_45 * ab_y[k] * gh_232[k]
                  + f_47 * ab_y[k] * gh_237[k]
                  - f_49 * ab_y[k] * gh_246[k]
                  + f_45 * gi_57[k]
                  - f_47 * gi_62[k]
                  + f_49 * gi_71[k]
                  - f_46 * gi_197[k]
                  + f_48 * gi_202[k]
                  - f_50 * gi_211[k]
                  + f_45 * gi_311[k]
                  - f_47 * gi_318[k]
                  + f_49 * gi_329[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_46, gh_53, gh_151, gh_158, gh_235, gh_242, gi_60, \
                         gi_67, gi_200, gi_207, gi_315, gi_324 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = -19.6875 * ab_x[k] * gh_46[k]
                   + 19.6875 * ab_x[k] * gh_53[k]
                   + 118.125 * ab_x[k] * gh_151[k]
                   - 118.125 * ab_x[k] * gh_158[k]
                   - 19.6875 * ab_y[k] * gh_235[k]
                   + 19.6875 * ab_y[k] * gh_242[k]
                   + 19.6875 * gi_60[k]
                   - 19.6875 * gi_67[k]
                   - 118.125 * gi_200[k]
                   + 118.125 * gi_207[k]
                   + 19.6875 * gi_315[k]
                   - 19.6875 * gi_324[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_43, gh_48, gh_50, gh_57, gh_59, gh_148, gh_153, \
                         gh_155, gh_162, gh_164, gh_232, gh_237, gh_239, gh_246, gh_248, \
                         gi_57, gi_62, gi_64, gi_71, gi_73, gi_197, gi_202, gi_204, gi_211, \
                         gi_213, gi_311, gi_318, gi_320, gi_329, \
                         gi_331 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = f_112 * ab_x[k] * gh_43[k]
                   + f_114 * ab_x[k] * gh_48[k]
                   - f_115 * ab_x[k] * gh_50[k]
                   - f_117 * ab_x[k] * gh_57[k]
                   + f_52 * ab_x[k] * gh_59[k]
                   - f_113 * ab_x[k] * gh_148[k]
                   - f_51 * ab_x[k] * gh_153[k]
                   + f_116 * ab_x[k] * gh_155[k]
                   + f_118 * ab_x[k] * gh_162[k]
                   - f_119 * ab_x[k] * gh_164[k]
                   + f_112 * ab_y[k] * gh_232[k]
                   + f_114 * ab_y[k] * gh_237[k]
                   - f_115 * ab_y[k] * gh_239[k]
                   - f_117 * ab_y[k] * gh_246[k]
                   + f_52 * ab_y[k] * gh_248[k]
                   - f_112 * gi_57[k]
                   - f_114 * gi_62[k]
                   + f_115 * gi_64[k]
                   + f_117 * gi_71[k]
                   - f_52 * gi_73[k]
                   + f_113 * gi_197[k]
                   + f_51 * gi_202[k]
                   - f_116 * gi_204[k]
                   - f_118 * gi_211[k]
                   + f_119 * gi_213[k]
                   - f_112 * gi_311[k]
                   - f_114 * gi_318[k]
                   + f_115 * gi_320[k]
                   + f_117 * gi_329[k]
                   - f_52 * gi_331[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_46, gh_53, gh_55, gh_151, gh_158, gh_160, gh_235, \
                         gh_242, gh_244, gi_60, gi_67, gi_69, gi_200, gi_207, gi_209, gi_315, \
                         gi_324, gi_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = f_134 * ab_x[k] * gh_46[k]
                   + f_134 * ab_x[k] * gh_53[k]
                   - f_66 * ab_x[k] * gh_55[k]
                   - f_135 * ab_x[k] * gh_151[k]
                   - f_135 * ab_x[k] * gh_158[k]
                   + f_136 * ab_x[k] * gh_160[k]
                   + f_134 * ab_y[k] * gh_235[k]
                   + f_134 * ab_y[k] * gh_242[k]
                   - f_66 * ab_y[k] * gh_244[k]
                   - f_134 * gi_60[k]
                   - f_134 * gi_67[k]
                   + f_66 * gi_69[k]
                   + f_135 * gi_200[k]
                   + f_135 * gi_207[k]
                   - f_136 * gi_209[k]
                   - f_134 * gi_315[k]
                   - f_134 * gi_324[k]
                   + f_66 * gi_326[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_43, gh_48, gh_50, gh_57, gh_59, gh_61, gh_148, gh_153, \
                         gh_155, gh_162, gh_164, gh_166, gh_232, gh_237, gh_239, gh_246, \
                         gh_248, gh_250, gi_57, gi_62, gi_64, gi_71, gi_73, gi_75, gi_197, \
                         gi_202, gi_204, gi_211, gi_213, gi_215, gi_311, gi_318, gi_320, \
                         gi_329, gi_331, gi_333 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = -f_154 * ab_x[k] * gh_43[k]
                   - f_156 * ab_x[k] * gh_48[k]
                   + f_157 * ab_x[k] * gh_50[k]
                   - f_154 * ab_x[k] * gh_57[k]
                   + f_157 * ab_x[k] * gh_59[k]
                   - f_59 * ab_x[k] * gh_61[k]
                   + f_155 * ab_x[k] * gh_148[k]
                   + f_157 * ab_x[k] * gh_153[k]
                   - f_158 * ab_x[k] * gh_155[k]
                   + f_155 * ab_x[k] * gh_162[k]
                   - f_158 * ab_x[k] * gh_164[k]
                   + f_60 * ab_x[k] * gh_166[k]
                   - f_154 * ab_y[k] * gh_232[k]
                   - f_156 * ab_y[k] * gh_237[k]
                   + f_157 * ab_y[k] * gh_239[k]
                   - f_154 * ab_y[k] * gh_246[k]
                   + f_157 * ab_y[k] * gh_248[k]
                   - f_59 * ab_y[k] * gh_250[k]
                   + f_154 * gi_57[k]
                   + f_156 * gi_62[k]
                   - f_157 * gi_64[k]
                   + f_154 * gi_71[k]
                   - f_157 * gi_73[k]
                   + f_59 * gi_75[k]
                   - f_155 * gi_197[k]
                   - f_157 * gi_202[k]
                   + f_158 * gi_204[k]
                   - f_155 * gi_211[k]
                   + f_158 * gi_213[k]
                   - f_60 * gi_215[k]
                   + f_154 * gi_311[k]
                   + f_156 * gi_318[k]
                   - f_157 * gi_320[k]
                   + f_154 * gi_329[k]
                   - f_157 * gi_331[k]
                   + f_59 * gi_333[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_44, gh_49, gh_51, gh_58, gh_60, gh_62, gh_149, gh_154, \
                         gh_156, gh_163, gh_165, gh_167, gh_233, gh_238, gh_240, gh_247, \
                         gh_249, gh_251, gi_58, gi_63, gi_65, gi_72, gi_74, gi_76, gi_198, \
                         gi_203, gi_205, gi_212, gi_214, gi_216, gi_312, gi_319, gi_321, \
                         gi_330, gi_332, gi_334 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = -f_162 * ab_x[k] * gh_44[k]
                   - f_164 * ab_x[k] * gh_49[k]
                   + f_166 * ab_x[k] * gh_51[k]
                   - f_162 * ab_x[k] * gh_58[k]
                   + f_166 * ab_x[k] * gh_60[k]
                   - f_168 * ab_x[k] * gh_62[k]
                   + f_163 * ab_x[k] * gh_149[k]
                   + f_165 * ab_x[k] * gh_154[k]
                   - f_167 * ab_x[k] * gh_156[k]
                   + f_163 * ab_x[k] * gh_163[k]
                   - f_167 * ab_x[k] * gh_165[k]
                   + f_169 * ab_x[k] * gh_167[k]
                   - f_162 * ab_y[k] * gh_233[k]
                   - f_164 * ab_y[k] * gh_238[k]
                   + f_166 * ab_y[k] * gh_240[k]
                   - f_162 * ab_y[k] * gh_247[k]
                   + f_166 * ab_y[k] * gh_249[k]
                   - f_168 * ab_y[k] * gh_251[k]
                   + f_162 * gi_58[k]
                   + f_164 * gi_63[k]
                   - f_166 * gi_65[k]
                   + f_162 * gi_72[k]
                   - f_166 * gi_74[k]
                   + f_168 * gi_76[k]
                   - f_163 * gi_198[k]
                   - f_165 * gi_203[k]
                   + f_167 * gi_205[k]
                   - f_163 * gi_212[k]
                   + f_167 * gi_214[k]
                   - f_169 * gi_216[k]
                   + f_162 * gi_312[k]
                   + f_164 * gi_319[k]
                   - f_166 * gi_321[k]
                   + f_162 * gi_330[k]
                   - f_166 * gi_332[k]
                   + f_168 * gi_334[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_42, gh_45, gh_47, gh_52, gh_54, gh_56, gh_147, gh_150, \
                         gh_152, gh_157, gh_159, gh_161, gh_231, gh_234, gh_236, gh_241, \
                         gh_243, gh_245, gi_56, gi_59, gi_61, gi_66, gi_68, gi_70, gi_196, \
                         gi_199, gi_201, gi_206, gi_208, gi_210, gi_309, gi_314, gi_316, \
                         gi_323, gi_325, gi_327 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = -f_154 * ab_x[k] * gh_42[k]
                   - f_156 * ab_x[k] * gh_45[k]
                   + f_157 * ab_x[k] * gh_47[k]
                   - f_154 * ab_x[k] * gh_52[k]
                   + f_157 * ab_x[k] * gh_54[k]
                   - f_59 * ab_x[k] * gh_56[k]
                   + f_155 * ab_x[k] * gh_147[k]
                   + f_157 * ab_x[k] * gh_150[k]
                   - f_158 * ab_x[k] * gh_152[k]
                   + f_155 * ab_x[k] * gh_157[k]
                   - f_158 * ab_x[k] * gh_159[k]
                   + f_60 * ab_x[k] * gh_161[k]
                   - f_154 * ab_y[k] * gh_231[k]
                   - f_156 * ab_y[k] * gh_234[k]
                   + f_157 * ab_y[k] * gh_236[k]
                   - f_154 * ab_y[k] * gh_241[k]
                   + f_157 * ab_y[k] * gh_243[k]
                   - f_59 * ab_y[k] * gh_245[k]
                   + f_154 * gi_56[k]
                   + f_156 * gi_59[k]
                   - f_157 * gi_61[k]
                   + f_154 * gi_66[k]
                   - f_157 * gi_68[k]
                   + f_59 * gi_70[k]
                   - f_155 * gi_196[k]
                   - f_157 * gi_199[k]
                   + f_158 * gi_201[k]
                   - f_155 * gi_206[k]
                   + f_158 * gi_208[k]
                   - f_60 * gi_210[k]
                   + f_154 * gi_309[k]
                   + f_156 * gi_314[k]
                   - f_157 * gi_316[k]
                   + f_154 * gi_323[k]
                   - f_157 * gi_325[k]
                   + f_59 * gi_327[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_44, gh_51, gh_58, gh_60, gh_149, gh_156, gh_163, \
                         gh_165, gh_233, gh_240, gh_247, gh_249, gi_58, gi_65, gi_72, gi_74, \
                         gi_198, gi_205, gi_212, gi_214, gi_312, gi_321, gi_330, \
                         gi_332 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_106[k] = f_170 * ab_x[k] * gh_44[k]
                   - f_134 * ab_x[k] * gh_51[k]
                   - f_170 * ab_x[k] * gh_58[k]
                   + f_134 * ab_x[k] * gh_60[k]
                   - f_171 * ab_x[k] * gh_149[k]
                   + f_135 * ab_x[k] * gh_156[k]
                   + f_171 * ab_x[k] * gh_163[k]
                   - f_135 * ab_x[k] * gh_165[k]
                   + f_170 * ab_y[k] * gh_233[k]
                   - f_134 * ab_y[k] * gh_240[k]
                   - f_170 * ab_y[k] * gh_247[k]
                   + f_134 * ab_y[k] * gh_249[k]
                   - f_170 * gi_58[k]
                   + f_134 * gi_65[k]
                   + f_170 * gi_72[k]
                   - f_134 * gi_74[k]
                   + f_171 * gi_198[k]
                   - f_135 * gi_205[k]
                   - f_171 * gi_212[k]
                   + f_135 * gi_214[k]
                   - f_170 * gi_312[k]
                   + f_134 * gi_321[k]
                   + f_170 * gi_330[k]
                   - f_134 * gi_332[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_42, gh_45, gh_47, gh_52, gh_54, gh_147, gh_150, \
                         gh_152, gh_157, gh_159, gh_231, gh_234, gh_236, gh_241, gh_243, \
                         gi_56, gi_59, gi_61, gi_66, gi_68, gi_196, gi_199, gi_201, gi_206, \
                         gi_208, gi_309, gi_314, gi_316, gi_323, \
                         gi_325 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = f_117 * ab_x[k] * gh_42[k]
                   - f_114 * ab_x[k] * gh_45[k]
                   - f_52 * ab_x[k] * gh_47[k]
                   - f_112 * ab_x[k] * gh_52[k]
                   + f_115 * ab_x[k] * gh_54[k]
                   - f_118 * ab_x[k] * gh_147[k]
                   + f_51 * ab_x[k] * gh_150[k]
                   + f_119 * ab_x[k] * gh_152[k]
                   + f_113 * ab_x[k] * gh_157[k]
                   - f_116 * ab_x[k] * gh_159[k]
                   + f_117 * ab_y[k] * gh_231[k]
                   - f_114 * ab_y[k] * gh_234[k]
                   - f_52 * ab_y[k] * gh_236[k]
                   - f_112 * ab_y[k] * gh_241[k]
                   + f_115 * ab_y[k] * gh_243[k]
                   - f_117 * gi_56[k]
                   + f_114 * gi_59[k]
                   + f_52 * gi_61[k]
                   + f_112 * gi_66[k]
                   - f_115 * gi_68[k]
                   + f_118 * gi_196[k]
                   - f_51 * gi_199[k]
                   - f_119 * gi_201[k]
                   - f_113 * gi_206[k]
                   + f_116 * gi_208[k]
                   - f_117 * gi_309[k]
                   + f_114 * gi_314[k]
                   + f_52 * gi_316[k]
                   + f_112 * gi_323[k]
                   - f_115 * gi_325[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_44, gh_49, gh_58, gh_149, gh_154, gh_163, gh_233, \
                         gh_238, gh_247, gi_58, gi_63, gi_72, gi_198, gi_203, gi_212, gi_312, \
                         gi_319, gi_330 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = -4.921875 * ab_x[k] * gh_44[k]
                   + 29.53125 * ab_x[k] * gh_49[k]
                   - 4.921875 * ab_x[k] * gh_58[k]
                   + 29.53125 * ab_x[k] * gh_149[k]
                   - 177.1875 * ab_x[k] * gh_154[k]
                   + 29.53125 * ab_x[k] * gh_163[k]
                   - 4.921875 * ab_y[k] * gh_233[k]
                   + 29.53125 * ab_y[k] * gh_238[k]
                   - 4.921875 * ab_y[k] * gh_247[k]
                   + 4.921875 * gi_58[k]
                   - 29.53125 * gi_63[k]
                   + 4.921875 * gi_72[k]
                   - 29.53125 * gi_198[k]
                   + 177.1875 * gi_203[k]
                   - 29.53125 * gi_212[k]
                   + 4.921875 * gi_312[k]
                   - 29.53125 * gi_319[k]
                   + 4.921875 * gi_330[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_42, gh_45, gh_52, gh_147, gh_150, gh_157, gh_231, \
                         gh_234, gh_241, gi_56, gi_59, gi_66, gi_196, gi_199, gi_206, gi_309, \
                         gi_314, gi_323 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = -f_49 * ab_x[k] * gh_42[k]
                   + f_47 * ab_x[k] * gh_45[k]
                   - f_45 * ab_x[k] * gh_52[k]
                   + f_50 * ab_x[k] * gh_147[k]
                   - f_48 * ab_x[k] * gh_150[k]
                   + f_46 * ab_x[k] * gh_157[k]
                   - f_49 * ab_y[k] * gh_231[k]
                   + f_47 * ab_y[k] * gh_234[k]
                   - f_45 * ab_y[k] * gh_241[k]
                   + f_49 * gi_56[k]
                   - f_47 * gi_59[k]
                   + f_45 * gi_66[k]
                   - f_50 * gi_196[k]
                   + f_48 * gi_199[k]
                   - f_46 * gi_206[k]
                   + f_49 * gi_309[k]
                   - f_47 * gi_314[k]
                   + f_45 * gi_323[k];
    }

#pragma omp simd aligned(ab_x, gh_1, gh_6, gh_15, gh_64, gh_69, gh_78, gh_211, gh_216, gh_225, \
                         gi_1, gi_6, gi_15, gi_85, gi_90, gi_99, gi_281, gi_286, \
                         gi_295 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = -2.4609375 * ab_x[k] * gh_1[k]
                   + 4.921875 * ab_x[k] * gh_6[k]
                   - 0.4921875 * ab_x[k] * gh_15[k]
                   + 24.609375 * ab_x[k] * gh_64[k]
                   - 49.21875 * ab_x[k] * gh_69[k]
                   + 4.921875 * ab_x[k] * gh_78[k]
                   - 12.3046875 * ab_x[k] * gh_211[k]
                   + 24.609375 * ab_x[k] * gh_216[k]
                   - 2.4609375 * ab_x[k] * gh_225[k]
                   + 2.4609375 * gi_1[k]
                   - 4.921875 * gi_6[k]
                   + 0.4921875 * gi_15[k]
                   - 24.609375 * gi_85[k]
                   + 49.21875 * gi_90[k]
                   - 4.921875 * gi_99[k]
                   + 12.3046875 * gi_281[k]
                   - 24.609375 * gi_286[k]
                   + 2.4609375 * gi_295[k];
    }

#pragma omp simd aligned(ab_x, gh_4, gh_11, gh_67, gh_74, gh_214, gh_221, gi_4, gi_11, gi_88, \
                         gi_95, gi_284, gi_291 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = -f_2 * ab_x[k] * gh_4[k]
                   + f_2 * ab_x[k] * gh_11[k]
                   + f_1 * ab_x[k] * gh_67[k]
                   - f_1 * ab_x[k] * gh_74[k]
                   - f_0 * ab_x[k] * gh_214[k]
                   + f_0 * ab_x[k] * gh_221[k]
                   + f_2 * gi_4[k]
                   - f_2 * gi_11[k]
                   - f_1 * gi_88[k]
                   + f_1 * gi_95[k]
                   + f_0 * gi_284[k]
                   - f_0 * gi_291[k];
    }

#pragma omp simd aligned(ab_x, gh_1, gh_6, gh_8, gh_15, gh_17, gh_64, gh_69, gh_71, gh_78, \
                         gh_80, gh_211, gh_216, gh_218, gh_225, gh_227, gi_1, gi_6, gi_8, \
                         gi_15, gi_17, gi_85, gi_90, gi_92, gi_99, gi_101, gi_281, gi_286, \
                         gi_288, gi_295, gi_297 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = f_12 * ab_x[k] * gh_1[k]
                   + f_13 * ab_x[k] * gh_6[k]
                   - f_14 * ab_x[k] * gh_8[k]
                   - f_15 * ab_x[k] * gh_15[k]
                   + f_16 * ab_x[k] * gh_17[k]
                   - f_8 * ab_x[k] * gh_64[k]
                   - f_9 * ab_x[k] * gh_69[k]
                   + f_10 * ab_x[k] * gh_71[k]
                   + f_4 * ab_x[k] * gh_78[k]
                   - f_11 * ab_x[k] * gh_80[k]
                   + f_3 * ab_x[k] * gh_211[k]
                   + f_4 * ab_x[k] * gh_216[k]
                   - f_5 * ab_x[k] * gh_218[k]
                   - f_6 * ab_x[k] * gh_225[k]
                   + f_7 * ab_x[k] * gh_227[k]
                   - f_12 * gi_1[k]
                   - f_13 * gi_6[k]
                   + f_14 * gi_8[k]
                   + f_15 * gi_15[k]
                   - f_16 * gi_17[k]
                   + f_8 * gi_85[k]
                   + f_9 * gi_90[k]
                   - f_10 * gi_92[k]
                   - f_4 * gi_99[k]
                   + f_11 * gi_101[k]
                   - f_3 * gi_281[k]
                   - f_4 * gi_286[k]
                   + f_5 * gi_288[k]
                   + f_6 * gi_295[k]
                   - f_7 * gi_297[k];
    }

#pragma omp simd aligned(ab_x, gh_4, gh_11, gh_13, gh_67, gh_74, gh_76, gh_214, gh_221, \
                         gh_223, gi_4, gi_11, gi_13, gi_88, gi_95, gi_97, gi_284, gi_291, \
                         gi_293 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = f_20 * ab_x[k] * gh_4[k]
                   + f_20 * ab_x[k] * gh_11[k]
                   - f_21 * ab_x[k] * gh_13[k]
                   - f_18 * ab_x[k] * gh_67[k]
                   - f_18 * ab_x[k] * gh_74[k]
                   + f_19 * ab_x[k] * gh_76[k]
                   + f_17 * ab_x[k] * gh_214[k]
                   + f_17 * ab_x[k] * gh_221[k]
                   - f_18 * ab_x[k] * gh_223[k]
                   - f_20 * gi_4[k]
                   - f_20 * gi_11[k]
                   + f_21 * gi_13[k]
                   + f_18 * gi_88[k]
                   + f_18 * gi_95[k]
                   - f_19 * gi_97[k]
                   - f_17 * gi_284[k]
                   - f_17 * gi_291[k]
                   + f_18 * gi_293[k];
    }

#pragma omp simd aligned(ab_x, gh_1, gh_6, gh_8, gh_15, gh_17, gh_19, gh_64, gh_69, gh_71, \
                         gh_78, gh_80, gh_82, gh_211, gh_216, gh_218, gh_225, gh_227, gh_229, \
                         gi_1, gi_6, gi_8, gi_15, gi_17, gi_19, gi_85, gi_90, gi_92, gi_99, \
                         gi_101, gi_103, gi_281, gi_286, gi_288, gi_295, gi_297, \
                         gi_299 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = -f_29 * ab_x[k] * gh_1[k]
                   - f_30 * ab_x[k] * gh_6[k]
                   + f_31 * ab_x[k] * gh_8[k]
                   - f_29 * ab_x[k] * gh_15[k]
                   + f_31 * ab_x[k] * gh_17[k]
                   - f_32 * ab_x[k] * gh_19[k]
                   + f_23 * ab_x[k] * gh_64[k]
                   + f_26 * ab_x[k] * gh_69[k]
                   - f_27 * ab_x[k] * gh_71[k]
                   + f_23 * ab_x[k] * gh_78[k]
                   - f_27 * ab_x[k] * gh_80[k]
                   + f_28 * ab_x[k] * gh_82[k]
                   - f_22 * ab_x[k] * gh_211[k]
                   - f_23 * ab_x[k] * gh_216[k]
                   + f_24 * ab_x[k] * gh_218[k]
                   - f_22 * ab_x[k] * gh_225[k]
                   + f_24 * ab_x[k] * gh_227[k]
                   - f_25 * ab_x[k] * gh_229[k]
                   + f_29 * gi_1[k]
                   + f_30 * gi_6[k]
                   - f_31 * gi_8[k]
                   + f_29 * gi_15[k]
                   - f_31 * gi_17[k]
                   + f_32 * gi_19[k]
                   - f_23 * gi_85[k]
                   - f_26 * gi_90[k]
                   + f_27 * gi_92[k]
                   - f_23 * gi_99[k]
                   + f_27 * gi_101[k]
                   - f_28 * gi_103[k]
                   + f_22 * gi_281[k]
                   + f_23 * gi_286[k]
                   - f_24 * gi_288[k]
                   + f_22 * gi_295[k]
                   - f_24 * gi_297[k]
                   + f_25 * gi_299[k];
    }

#pragma omp simd aligned(ab_x, gh_2, gh_7, gh_9, gh_16, gh_18, gh_20, gh_65, gh_70, gh_72, \
                         gh_79, gh_81, gh_83, gh_212, gh_217, gh_219, gh_226, gh_228, gh_230, \
                         gi_2, gi_7, gi_9, gi_16, gi_18, gi_20, gi_86, gi_91, gi_93, gi_100, \
                         gi_102, gi_104, gi_282, gi_287, gi_289, gi_296, gi_298, \
                         gi_300 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = -f_40 * ab_x[k] * gh_2[k]
                   - f_41 * ab_x[k] * gh_7[k]
                   + f_36 * ab_x[k] * gh_9[k]
                   - f_40 * ab_x[k] * gh_16[k]
                   + f_36 * ab_x[k] * gh_18[k]
                   - f_42 * ab_x[k] * gh_20[k]
                   + f_34 * ab_x[k] * gh_65[k]
                   + f_37 * ab_x[k] * gh_70[k]
                   - f_38 * ab_x[k] * gh_72[k]
                   + f_34 * ab_x[k] * gh_79[k]
                   - f_38 * ab_x[k] * gh_81[k]
                   + f_39 * ab_x[k] * gh_83[k]
                   - f_33 * ab_x[k] * gh_212[k]
                   - f_34 * ab_x[k] * gh_217[k]
                   + f_35 * ab_x[k] * gh_219[k]
                   - f_33 * ab_x[k] * gh_226[k]
                   + f_35 * ab_x[k] * gh_228[k]
                   - f_36 * ab_x[k] * gh_230[k]
                   + f_40 * gi_2[k]
                   + f_41 * gi_7[k]
                   - f_36 * gi_9[k]
                   + f_40 * gi_16[k]
                   - f_36 * gi_18[k]
                   + f_42 * gi_20[k]
                   - f_34 * gi_86[k]
                   - f_37 * gi_91[k]
                   + f_38 * gi_93[k]
                   - f_34 * gi_100[k]
                   + f_38 * gi_102[k]
                   - f_39 * gi_104[k]
                   + f_33 * gi_282[k]
                   + f_34 * gi_287[k]
                   - f_35 * gi_289[k]
                   + f_33 * gi_296[k]
                   - f_35 * gi_298[k]
                   + f_36 * gi_300[k];
    }

#pragma omp simd aligned(ab_x, gh_0, gh_3, gh_5, gh_10, gh_12, gh_14, gh_63, gh_66, gh_68, \
                         gh_73, gh_75, gh_77, gh_210, gh_213, gh_215, gh_220, gh_222, gh_224, \
                         gi_0, gi_3, gi_5, gi_10, gi_12, gi_14, gi_84, gi_87, gi_89, gi_94, \
                         gi_96, gi_98, gi_280, gi_283, gi_285, gi_290, gi_292, \
                         gi_294 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = -f_29 * ab_x[k] * gh_0[k]
                   - f_30 * ab_x[k] * gh_3[k]
                   + f_31 * ab_x[k] * gh_5[k]
                   - f_29 * ab_x[k] * gh_10[k]
                   + f_31 * ab_x[k] * gh_12[k]
                   - f_32 * ab_x[k] * gh_14[k]
                   + f_23 * ab_x[k] * gh_63[k]
                   + f_26 * ab_x[k] * gh_66[k]
                   - f_27 * ab_x[k] * gh_68[k]
                   + f_23 * ab_x[k] * gh_73[k]
                   - f_27 * ab_x[k] * gh_75[k]
                   + f_28 * ab_x[k] * gh_77[k]
                   - f_22 * ab_x[k] * gh_210[k]
                   - f_23 * ab_x[k] * gh_213[k]
                   + f_24 * ab_x[k] * gh_215[k]
                   - f_22 * ab_x[k] * gh_220[k]
                   + f_24 * ab_x[k] * gh_222[k]
                   - f_25 * ab_x[k] * gh_224[k]
                   + f_29 * gi_0[k]
                   + f_30 * gi_3[k]
                   - f_31 * gi_5[k]
                   + f_29 * gi_10[k]
                   - f_31 * gi_12[k]
                   + f_32 * gi_14[k]
                   - f_23 * gi_84[k]
                   - f_26 * gi_87[k]
                   + f_27 * gi_89[k]
                   - f_23 * gi_94[k]
                   + f_27 * gi_96[k]
                   - f_28 * gi_98[k]
                   + f_22 * gi_280[k]
                   + f_23 * gi_283[k]
                   - f_24 * gi_285[k]
                   + f_22 * gi_290[k]
                   - f_24 * gi_292[k]
                   + f_25 * gi_294[k];
    }

#pragma omp simd aligned(ab_x, gh_2, gh_9, gh_16, gh_18, gh_65, gh_72, gh_79, gh_81, gh_212, \
                         gh_219, gh_226, gh_228, gi_2, gi_9, gi_16, gi_18, gi_86, gi_93, \
                         gi_100, gi_102, gi_282, gi_289, gi_296, \
                         gi_298 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_117[k] = f_44 * ab_x[k] * gh_2[k]
                   - f_20 * ab_x[k] * gh_9[k]
                   - f_44 * ab_x[k] * gh_16[k]
                   + f_20 * ab_x[k] * gh_18[k]
                   - f_17 * ab_x[k] * gh_65[k]
                   + f_18 * ab_x[k] * gh_72[k]
                   + f_17 * ab_x[k] * gh_79[k]
                   - f_18 * ab_x[k] * gh_81[k]
                   + f_43 * ab_x[k] * gh_212[k]
                   - f_17 * ab_x[k] * gh_219[k]
                   - f_43 * ab_x[k] * gh_226[k]
                   + f_17 * ab_x[k] * gh_228[k]
                   - f_44 * gi_2[k]
                   + f_20 * gi_9[k]
                   + f_44 * gi_16[k]
                   - f_20 * gi_18[k]
                   + f_17 * gi_86[k]
                   - f_18 * gi_93[k]
                   - f_17 * gi_100[k]
                   + f_18 * gi_102[k]
                   - f_43 * gi_282[k]
                   + f_17 * gi_289[k]
                   + f_43 * gi_296[k]
                   - f_17 * gi_298[k];
    }

#pragma omp simd aligned(ab_x, gh_0, gh_3, gh_5, gh_10, gh_12, gh_63, gh_66, gh_68, gh_73, \
                         gh_75, gh_210, gh_213, gh_215, gh_220, gh_222, gi_0, gi_3, gi_5, \
                         gi_10, gi_12, gi_84, gi_87, gi_89, gi_94, gi_96, gi_280, gi_283, \
                         gi_285, gi_290, gi_292 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_118[k] = f_15 * ab_x[k] * gh_0[k]
                   - f_13 * ab_x[k] * gh_3[k]
                   - f_16 * ab_x[k] * gh_5[k]
                   - f_12 * ab_x[k] * gh_10[k]
                   + f_14 * ab_x[k] * gh_12[k]
                   - f_4 * ab_x[k] * gh_63[k]
                   + f_9 * ab_x[k] * gh_66[k]
                   + f_11 * ab_x[k] * gh_68[k]
                   + f_8 * ab_x[k] * gh_73[k]
                   - f_10 * ab_x[k] * gh_75[k]
                   + f_6 * ab_x[k] * gh_210[k]
                   - f_4 * ab_x[k] * gh_213[k]
                   - f_7 * ab_x[k] * gh_215[k]
                   - f_3 * ab_x[k] * gh_220[k]
                   + f_5 * ab_x[k] * gh_222[k]
                   - f_15 * gi_0[k]
                   + f_13 * gi_3[k]
                   + f_16 * gi_5[k]
                   + f_12 * gi_10[k]
                   - f_14 * gi_12[k]
                   + f_4 * gi_84[k]
                   - f_9 * gi_87[k]
                   - f_11 * gi_89[k]
                   - f_8 * gi_94[k]
                   + f_10 * gi_96[k]
                   - f_6 * gi_280[k]
                   + f_4 * gi_283[k]
                   + f_7 * gi_285[k]
                   + f_3 * gi_290[k]
                   - f_5 * gi_292[k];
    }

#pragma omp simd aligned(ab_x, gh_2, gh_7, gh_16, gh_65, gh_70, gh_79, gh_212, gh_217, gh_226, \
                         gi_2, gi_7, gi_16, gi_86, gi_91, gi_100, gi_282, gi_287, \
                         gi_296 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_119[k] = -f_49 * ab_x[k] * gh_2[k]
                   + f_50 * ab_x[k] * gh_7[k]
                   - f_49 * ab_x[k] * gh_16[k]
                   + f_47 * ab_x[k] * gh_65[k]
                   - f_48 * ab_x[k] * gh_70[k]
                   + f_47 * ab_x[k] * gh_79[k]
                   - f_45 * ab_x[k] * gh_212[k]
                   + f_46 * ab_x[k] * gh_217[k]
                   - f_45 * ab_x[k] * gh_226[k]
                   + f_49 * gi_2[k]
                   - f_50 * gi_7[k]
                   + f_49 * gi_16[k]
                   - f_47 * gi_86[k]
                   + f_48 * gi_91[k]
                   - f_47 * gi_100[k]
                   + f_45 * gi_282[k]
                   - f_46 * gi_287[k]
                   + f_45 * gi_296[k];
    }

#pragma omp simd aligned(ab_x, gh_0, gh_3, gh_10, gh_63, gh_66, gh_73, gh_210, gh_213, gh_220, \
                         gi_0, gi_3, gi_10, gi_84, gi_87, gi_94, gi_280, gi_283, \
                         gi_290 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_120[k] = -0.4921875 * ab_x[k] * gh_0[k]
                   + 4.921875 * ab_x[k] * gh_3[k]
                   - 2.4609375 * ab_x[k] * gh_10[k]
                   + 4.921875 * ab_x[k] * gh_63[k]
                   - 49.21875 * ab_x[k] * gh_66[k]
                   + 24.609375 * ab_x[k] * gh_73[k]
                   - 2.4609375 * ab_x[k] * gh_210[k]
                   + 24.609375 * ab_x[k] * gh_213[k]
                   - 12.3046875 * ab_x[k] * gh_220[k]
                   + 0.4921875 * gi_0[k]
                   - 4.921875 * gi_3[k]
                   + 2.4609375 * gi_10[k]
                   - 4.921875 * gi_84[k]
                   + 49.21875 * gi_87[k]
                   - 24.609375 * gi_94[k]
                   + 2.4609375 * gi_280[k]
                   - 24.609375 * gi_283[k]
                   + 12.3046875 * gi_290[k];
    }
}

auto
compute_hrr_hh_sph_tri(double *values, const size_t nvalues, CSimdMatrix &buffer,
                       const CSimdMatrix &coordinates, const size_t gh, const size_t gi,
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_2 = buffer.data(gh + 2);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_7 = buffer.data(gh + 7);
    const auto *gh_9 = buffer.data(gh + 9);
    const auto *gh_10 = buffer.data(gh + 10);
    const auto *gh_12 = buffer.data(gh + 12);
    const auto *gh_14 = buffer.data(gh + 14);
    const auto *gh_16 = buffer.data(gh + 16);
    const auto *gh_18 = buffer.data(gh + 18);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_23 = buffer.data(gh + 23);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_25 = buffer.data(gh + 25);
    const auto *gh_26 = buffer.data(gh + 26);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_28 = buffer.data(gh + 28);
    const auto *gh_29 = buffer.data(gh + 29);
    const auto *gh_30 = buffer.data(gh + 30);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_32 = buffer.data(gh + 32);
    const auto *gh_33 = buffer.data(gh + 33);
    const auto *gh_34 = buffer.data(gh + 34);
    const auto *gh_35 = buffer.data(gh + 35);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_37 = buffer.data(gh + 37);
    const auto *gh_38 = buffer.data(gh + 38);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_40 = buffer.data(gh + 40);
    const auto *gh_41 = buffer.data(gh + 41);
    const auto *gh_42 = buffer.data(gh + 42);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_45 = buffer.data(gh + 45);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_49 = buffer.data(gh + 49);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_52 = buffer.data(gh + 52);
    const auto *gh_54 = buffer.data(gh + 54);
    const auto *gh_56 = buffer.data(gh + 56);
    const auto *gh_58 = buffer.data(gh + 58);
    const auto *gh_60 = buffer.data(gh + 60);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_65 = buffer.data(gh + 65);
    const auto *gh_66 = buffer.data(gh + 66);
    const auto *gh_68 = buffer.data(gh + 68);
    const auto *gh_70 = buffer.data(gh + 70);
    const auto *gh_72 = buffer.data(gh + 72);
    const auto *gh_73 = buffer.data(gh + 73);
    const auto *gh_75 = buffer.data(gh + 75);
    const auto *gh_77 = buffer.data(gh + 77);
    const auto *gh_79 = buffer.data(gh + 79);
    const auto *gh_81 = buffer.data(gh + 81);
    const auto *gh_84 = buffer.data(gh + 84);
    const auto *gh_85 = buffer.data(gh + 85);
    const auto *gh_86 = buffer.data(gh + 86);
    const auto *gh_87 = buffer.data(gh + 87);
    const auto *gh_88 = buffer.data(gh + 88);
    const auto *gh_89 = buffer.data(gh + 89);
    const auto *gh_90 = buffer.data(gh + 90);
    const auto *gh_91 = buffer.data(gh + 91);
    const auto *gh_92 = buffer.data(gh + 92);
    const auto *gh_93 = buffer.data(gh + 93);
    const auto *gh_94 = buffer.data(gh + 94);
    const auto *gh_95 = buffer.data(gh + 95);
    const auto *gh_96 = buffer.data(gh + 96);
    const auto *gh_97 = buffer.data(gh + 97);
    const auto *gh_98 = buffer.data(gh + 98);
    const auto *gh_99 = buffer.data(gh + 99);
    const auto *gh_100 = buffer.data(gh + 100);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_103 = buffer.data(gh + 103);
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_107 = buffer.data(gh + 107);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_110 = buffer.data(gh + 110);
    const auto *gh_112 = buffer.data(gh + 112);
    const auto *gh_114 = buffer.data(gh + 114);
    const auto *gh_115 = buffer.data(gh + 115);
    const auto *gh_117 = buffer.data(gh + 117);
    const auto *gh_119 = buffer.data(gh + 119);
    const auto *gh_121 = buffer.data(gh + 121);
    const auto *gh_123 = buffer.data(gh + 123);
    const auto *gh_126 = buffer.data(gh + 126);
    const auto *gh_127 = buffer.data(gh + 127);
    const auto *gh_128 = buffer.data(gh + 128);
    const auto *gh_129 = buffer.data(gh + 129);
    const auto *gh_130 = buffer.data(gh + 130);
    const auto *gh_131 = buffer.data(gh + 131);
    const auto *gh_132 = buffer.data(gh + 132);
    const auto *gh_133 = buffer.data(gh + 133);
    const auto *gh_134 = buffer.data(gh + 134);
    const auto *gh_135 = buffer.data(gh + 135);
    const auto *gh_136 = buffer.data(gh + 136);
    const auto *gh_137 = buffer.data(gh + 137);
    const auto *gh_138 = buffer.data(gh + 138);
    const auto *gh_139 = buffer.data(gh + 139);
    const auto *gh_140 = buffer.data(gh + 140);
    const auto *gh_141 = buffer.data(gh + 141);
    const auto *gh_142 = buffer.data(gh + 142);
    const auto *gh_143 = buffer.data(gh + 143);
    const auto *gh_144 = buffer.data(gh + 144);
    const auto *gh_145 = buffer.data(gh + 145);
    const auto *gh_146 = buffer.data(gh + 146);
    const auto *gh_147 = buffer.data(gh + 147);
    const auto *gh_149 = buffer.data(gh + 149);
    const auto *gh_150 = buffer.data(gh + 150);
    const auto *gh_152 = buffer.data(gh + 152);
    const auto *gh_154 = buffer.data(gh + 154);
    const auto *gh_156 = buffer.data(gh + 156);
    const auto *gh_157 = buffer.data(gh + 157);
    const auto *gh_159 = buffer.data(gh + 159);
    const auto *gh_161 = buffer.data(gh + 161);
    const auto *gh_163 = buffer.data(gh + 163);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_167 = buffer.data(gh + 167);
    const auto *gh_168 = buffer.data(gh + 168);
    const auto *gh_169 = buffer.data(gh + 169);
    const auto *gh_170 = buffer.data(gh + 170);
    const auto *gh_171 = buffer.data(gh + 171);
    const auto *gh_172 = buffer.data(gh + 172);
    const auto *gh_173 = buffer.data(gh + 173);
    const auto *gh_174 = buffer.data(gh + 174);
    const auto *gh_175 = buffer.data(gh + 175);
    const auto *gh_176 = buffer.data(gh + 176);
    const auto *gh_177 = buffer.data(gh + 177);
    const auto *gh_178 = buffer.data(gh + 178);
    const auto *gh_179 = buffer.data(gh + 179);
    const auto *gh_180 = buffer.data(gh + 180);
    const auto *gh_181 = buffer.data(gh + 181);
    const auto *gh_182 = buffer.data(gh + 182);
    const auto *gh_183 = buffer.data(gh + 183);
    const auto *gh_184 = buffer.data(gh + 184);
    const auto *gh_185 = buffer.data(gh + 185);
    const auto *gh_186 = buffer.data(gh + 186);
    const auto *gh_187 = buffer.data(gh + 187);
    const auto *gh_188 = buffer.data(gh + 188);
    const auto *gh_189 = buffer.data(gh + 189);
    const auto *gh_191 = buffer.data(gh + 191);
    const auto *gh_192 = buffer.data(gh + 192);
    const auto *gh_194 = buffer.data(gh + 194);
    const auto *gh_196 = buffer.data(gh + 196);
    const auto *gh_198 = buffer.data(gh + 198);
    const auto *gh_199 = buffer.data(gh + 199);
    const auto *gh_201 = buffer.data(gh + 201);
    const auto *gh_203 = buffer.data(gh + 203);
    const auto *gh_205 = buffer.data(gh + 205);
    const auto *gh_207 = buffer.data(gh + 207);
    const auto *gh_209 = buffer.data(gh + 209);
    const auto *gh_210 = buffer.data(gh + 210);
    const auto *gh_211 = buffer.data(gh + 211);
    const auto *gh_212 = buffer.data(gh + 212);
    const auto *gh_213 = buffer.data(gh + 213);
    const auto *gh_214 = buffer.data(gh + 214);
    const auto *gh_215 = buffer.data(gh + 215);
    const auto *gh_216 = buffer.data(gh + 216);
    const auto *gh_217 = buffer.data(gh + 217);
    const auto *gh_218 = buffer.data(gh + 218);
    const auto *gh_219 = buffer.data(gh + 219);
    const auto *gh_220 = buffer.data(gh + 220);
    const auto *gh_221 = buffer.data(gh + 221);
    const auto *gh_222 = buffer.data(gh + 222);
    const auto *gh_223 = buffer.data(gh + 223);
    const auto *gh_224 = buffer.data(gh + 224);
    const auto *gh_225 = buffer.data(gh + 225);
    const auto *gh_226 = buffer.data(gh + 226);
    const auto *gh_227 = buffer.data(gh + 227);
    const auto *gh_228 = buffer.data(gh + 228);
    const auto *gh_229 = buffer.data(gh + 229);
    const auto *gh_230 = buffer.data(gh + 230);
    const auto *gh_231 = buffer.data(gh + 231);
    const auto *gh_232 = buffer.data(gh + 232);
    const auto *gh_233 = buffer.data(gh + 233);
    const auto *gh_234 = buffer.data(gh + 234);
    const auto *gh_235 = buffer.data(gh + 235);
    const auto *gh_236 = buffer.data(gh + 236);
    const auto *gh_237 = buffer.data(gh + 237);
    const auto *gh_238 = buffer.data(gh + 238);
    const auto *gh_239 = buffer.data(gh + 239);
    const auto *gh_240 = buffer.data(gh + 240);
    const auto *gh_241 = buffer.data(gh + 241);
    const auto *gh_242 = buffer.data(gh + 242);
    const auto *gh_243 = buffer.data(gh + 243);
    const auto *gh_244 = buffer.data(gh + 244);
    const auto *gh_245 = buffer.data(gh + 245);
    const auto *gh_246 = buffer.data(gh + 246);
    const auto *gh_247 = buffer.data(gh + 247);
    const auto *gh_248 = buffer.data(gh + 248);
    const auto *gh_249 = buffer.data(gh + 249);
    const auto *gh_250 = buffer.data(gh + 250);
    const auto *gh_251 = buffer.data(gh + 251);
    const auto *gh_252 = buffer.data(gh + 252);
    const auto *gh_253 = buffer.data(gh + 253);
    const auto *gh_254 = buffer.data(gh + 254);
    const auto *gh_255 = buffer.data(gh + 255);
    const auto *gh_256 = buffer.data(gh + 256);
    const auto *gh_257 = buffer.data(gh + 257);
    const auto *gh_258 = buffer.data(gh + 258);
    const auto *gh_259 = buffer.data(gh + 259);
    const auto *gh_260 = buffer.data(gh + 260);
    const auto *gh_261 = buffer.data(gh + 261);
    const auto *gh_262 = buffer.data(gh + 262);
    const auto *gh_263 = buffer.data(gh + 263);
    const auto *gh_264 = buffer.data(gh + 264);
    const auto *gh_265 = buffer.data(gh + 265);
    const auto *gh_266 = buffer.data(gh + 266);
    const auto *gh_267 = buffer.data(gh + 267);
    const auto *gh_268 = buffer.data(gh + 268);
    const auto *gh_269 = buffer.data(gh + 269);
    const auto *gh_270 = buffer.data(gh + 270);
    const auto *gh_271 = buffer.data(gh + 271);
    const auto *gh_272 = buffer.data(gh + 272);
    const auto *gh_273 = buffer.data(gh + 273);
    const auto *gh_274 = buffer.data(gh + 274);
    const auto *gh_275 = buffer.data(gh + 275);
    const auto *gh_276 = buffer.data(gh + 276);
    const auto *gh_277 = buffer.data(gh + 277);
    const auto *gh_278 = buffer.data(gh + 278);
    const auto *gh_279 = buffer.data(gh + 279);
    const auto *gh_280 = buffer.data(gh + 280);
    const auto *gh_281 = buffer.data(gh + 281);
    const auto *gh_282 = buffer.data(gh + 282);
    const auto *gh_283 = buffer.data(gh + 283);
    const auto *gh_284 = buffer.data(gh + 284);
    const auto *gh_285 = buffer.data(gh + 285);
    const auto *gh_286 = buffer.data(gh + 286);
    const auto *gh_287 = buffer.data(gh + 287);
    const auto *gh_288 = buffer.data(gh + 288);
    const auto *gh_289 = buffer.data(gh + 289);
    const auto *gh_290 = buffer.data(gh + 290);
    const auto *gh_291 = buffer.data(gh + 291);
    const auto *gh_292 = buffer.data(gh + 292);
    const auto *gh_293 = buffer.data(gh + 293);
    const auto *gh_294 = buffer.data(gh + 294);
    const auto *gh_295 = buffer.data(gh + 295);
    const auto *gh_296 = buffer.data(gh + 296);
    const auto *gh_297 = buffer.data(gh + 297);
    const auto *gh_299 = buffer.data(gh + 299);
    const auto *gh_300 = buffer.data(gh + 300);
    const auto *gh_301 = buffer.data(gh + 301);
    const auto *gh_302 = buffer.data(gh + 302);
    const auto *gh_303 = buffer.data(gh + 303);
    const auto *gh_304 = buffer.data(gh + 304);
    const auto *gh_306 = buffer.data(gh + 306);
    const auto *gh_308 = buffer.data(gh + 308);
    const auto *gh_309 = buffer.data(gh + 309);
    const auto *gh_310 = buffer.data(gh + 310);
    const auto *gh_311 = buffer.data(gh + 311);
    const auto *gh_312 = buffer.data(gh + 312);
    const auto *gh_313 = buffer.data(gh + 313);
    const auto *gh_314 = buffer.data(gh + 314);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_2 = buffer.data(gi + 2);
    const auto *gi_3 = buffer.data(gi + 3);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_9 = buffer.data(gi + 9);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_16 = buffer.data(gi + 16);
    const auto *gi_18 = buffer.data(gi + 18);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_30 = buffer.data(gi + 30);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_35 = buffer.data(gi + 35);
    const auto *gi_36 = buffer.data(gi + 36);
    const auto *gi_37 = buffer.data(gi + 37);
    const auto *gi_38 = buffer.data(gi + 38);
    const auto *gi_39 = buffer.data(gi + 39);
    const auto *gi_40 = buffer.data(gi + 40);
    const auto *gi_41 = buffer.data(gi + 41);
    const auto *gi_42 = buffer.data(gi + 42);
    const auto *gi_43 = buffer.data(gi + 43);
    const auto *gi_44 = buffer.data(gi + 44);
    const auto *gi_45 = buffer.data(gi + 45);
    const auto *gi_46 = buffer.data(gi + 46);
    const auto *gi_47 = buffer.data(gi + 47);
    const auto *gi_48 = buffer.data(gi + 48);
    const auto *gi_56 = buffer.data(gi + 56);
    const auto *gi_58 = buffer.data(gi + 58);
    const auto *gi_59 = buffer.data(gi + 59);
    const auto *gi_61 = buffer.data(gi + 61);
    const auto *gi_63 = buffer.data(gi + 63);
    const auto *gi_65 = buffer.data(gi + 65);
    const auto *gi_66 = buffer.data(gi + 66);
    const auto *gi_68 = buffer.data(gi + 68);
    const auto *gi_70 = buffer.data(gi + 70);
    const auto *gi_72 = buffer.data(gi + 72);
    const auto *gi_74 = buffer.data(gi + 74);
    const auto *gi_76 = buffer.data(gi + 76);
    const auto *gi_84 = buffer.data(gi + 84);
    const auto *gi_86 = buffer.data(gi + 86);
    const auto *gi_87 = buffer.data(gi + 87);
    const auto *gi_89 = buffer.data(gi + 89);
    const auto *gi_91 = buffer.data(gi + 91);
    const auto *gi_93 = buffer.data(gi + 93);
    const auto *gi_94 = buffer.data(gi + 94);
    const auto *gi_96 = buffer.data(gi + 96);
    const auto *gi_98 = buffer.data(gi + 98);
    const auto *gi_100 = buffer.data(gi + 100);
    const auto *gi_102 = buffer.data(gi + 102);
    const auto *gi_112 = buffer.data(gi + 112);
    const auto *gi_113 = buffer.data(gi + 113);
    const auto *gi_114 = buffer.data(gi + 114);
    const auto *gi_115 = buffer.data(gi + 115);
    const auto *gi_116 = buffer.data(gi + 116);
    const auto *gi_117 = buffer.data(gi + 117);
    const auto *gi_118 = buffer.data(gi + 118);
    const auto *gi_119 = buffer.data(gi + 119);
    const auto *gi_120 = buffer.data(gi + 120);
    const auto *gi_121 = buffer.data(gi + 121);
    const auto *gi_122 = buffer.data(gi + 122);
    const auto *gi_123 = buffer.data(gi + 123);
    const auto *gi_124 = buffer.data(gi + 124);
    const auto *gi_125 = buffer.data(gi + 125);
    const auto *gi_126 = buffer.data(gi + 126);
    const auto *gi_127 = buffer.data(gi + 127);
    const auto *gi_128 = buffer.data(gi + 128);
    const auto *gi_129 = buffer.data(gi + 129);
    const auto *gi_130 = buffer.data(gi + 130);
    const auto *gi_131 = buffer.data(gi + 131);
    const auto *gi_132 = buffer.data(gi + 132);
    const auto *gi_140 = buffer.data(gi + 140);
    const auto *gi_142 = buffer.data(gi + 142);
    const auto *gi_143 = buffer.data(gi + 143);
    const auto *gi_145 = buffer.data(gi + 145);
    const auto *gi_147 = buffer.data(gi + 147);
    const auto *gi_149 = buffer.data(gi + 149);
    const auto *gi_150 = buffer.data(gi + 150);
    const auto *gi_152 = buffer.data(gi + 152);
    const auto *gi_154 = buffer.data(gi + 154);
    const auto *gi_156 = buffer.data(gi + 156);
    const auto *gi_158 = buffer.data(gi + 158);
    const auto *gi_168 = buffer.data(gi + 168);
    const auto *gi_169 = buffer.data(gi + 169);
    const auto *gi_170 = buffer.data(gi + 170);
    const auto *gi_171 = buffer.data(gi + 171);
    const auto *gi_172 = buffer.data(gi + 172);
    const auto *gi_173 = buffer.data(gi + 173);
    const auto *gi_174 = buffer.data(gi + 174);
    const auto *gi_175 = buffer.data(gi + 175);
    const auto *gi_176 = buffer.data(gi + 176);
    const auto *gi_177 = buffer.data(gi + 177);
    const auto *gi_178 = buffer.data(gi + 178);
    const auto *gi_179 = buffer.data(gi + 179);
    const auto *gi_180 = buffer.data(gi + 180);
    const auto *gi_181 = buffer.data(gi + 181);
    const auto *gi_182 = buffer.data(gi + 182);
    const auto *gi_183 = buffer.data(gi + 183);
    const auto *gi_184 = buffer.data(gi + 184);
    const auto *gi_185 = buffer.data(gi + 185);
    const auto *gi_186 = buffer.data(gi + 186);
    const auto *gi_187 = buffer.data(gi + 187);
    const auto *gi_188 = buffer.data(gi + 188);
    const auto *gi_196 = buffer.data(gi + 196);
    const auto *gi_198 = buffer.data(gi + 198);
    const auto *gi_199 = buffer.data(gi + 199);
    const auto *gi_201 = buffer.data(gi + 201);
    const auto *gi_203 = buffer.data(gi + 203);
    const auto *gi_205 = buffer.data(gi + 205);
    const auto *gi_206 = buffer.data(gi + 206);
    const auto *gi_208 = buffer.data(gi + 208);
    const auto *gi_210 = buffer.data(gi + 210);
    const auto *gi_212 = buffer.data(gi + 212);
    const auto *gi_214 = buffer.data(gi + 214);
    const auto *gi_216 = buffer.data(gi + 216);
    const auto *gi_224 = buffer.data(gi + 224);
    const auto *gi_225 = buffer.data(gi + 225);
    const auto *gi_226 = buffer.data(gi + 226);
    const auto *gi_227 = buffer.data(gi + 227);
    const auto *gi_228 = buffer.data(gi + 228);
    const auto *gi_229 = buffer.data(gi + 229);
    const auto *gi_230 = buffer.data(gi + 230);
    const auto *gi_231 = buffer.data(gi + 231);
    const auto *gi_232 = buffer.data(gi + 232);
    const auto *gi_233 = buffer.data(gi + 233);
    const auto *gi_234 = buffer.data(gi + 234);
    const auto *gi_235 = buffer.data(gi + 235);
    const auto *gi_236 = buffer.data(gi + 236);
    const auto *gi_237 = buffer.data(gi + 237);
    const auto *gi_238 = buffer.data(gi + 238);
    const auto *gi_239 = buffer.data(gi + 239);
    const auto *gi_240 = buffer.data(gi + 240);
    const auto *gi_241 = buffer.data(gi + 241);
    const auto *gi_242 = buffer.data(gi + 242);
    const auto *gi_243 = buffer.data(gi + 243);
    const auto *gi_244 = buffer.data(gi + 244);
    const auto *gi_252 = buffer.data(gi + 252);
    const auto *gi_254 = buffer.data(gi + 254);
    const auto *gi_255 = buffer.data(gi + 255);
    const auto *gi_257 = buffer.data(gi + 257);
    const auto *gi_259 = buffer.data(gi + 259);
    const auto *gi_261 = buffer.data(gi + 261);
    const auto *gi_262 = buffer.data(gi + 262);
    const auto *gi_264 = buffer.data(gi + 264);
    const auto *gi_266 = buffer.data(gi + 266);
    const auto *gi_268 = buffer.data(gi + 268);
    const auto *gi_270 = buffer.data(gi + 270);
    const auto *gi_272 = buffer.data(gi + 272);
    const auto *gi_280 = buffer.data(gi + 280);
    const auto *gi_281 = buffer.data(gi + 281);
    const auto *gi_282 = buffer.data(gi + 282);
    const auto *gi_283 = buffer.data(gi + 283);
    const auto *gi_284 = buffer.data(gi + 284);
    const auto *gi_285 = buffer.data(gi + 285);
    const auto *gi_286 = buffer.data(gi + 286);
    const auto *gi_287 = buffer.data(gi + 287);
    const auto *gi_288 = buffer.data(gi + 288);
    const auto *gi_289 = buffer.data(gi + 289);
    const auto *gi_290 = buffer.data(gi + 290);
    const auto *gi_291 = buffer.data(gi + 291);
    const auto *gi_292 = buffer.data(gi + 292);
    const auto *gi_293 = buffer.data(gi + 293);
    const auto *gi_294 = buffer.data(gi + 294);
    const auto *gi_295 = buffer.data(gi + 295);
    const auto *gi_296 = buffer.data(gi + 296);
    const auto *gi_297 = buffer.data(gi + 297);
    const auto *gi_298 = buffer.data(gi + 298);
    const auto *gi_299 = buffer.data(gi + 299);
    const auto *gi_301 = buffer.data(gi + 301);
    const auto *gi_302 = buffer.data(gi + 302);
    const auto *gi_303 = buffer.data(gi + 303);
    const auto *gi_304 = buffer.data(gi + 304);
    const auto *gi_305 = buffer.data(gi + 305);
    const auto *gi_306 = buffer.data(gi + 306);
    const auto *gi_308 = buffer.data(gi + 308);
    const auto *gi_309 = buffer.data(gi + 309);
    const auto *gi_310 = buffer.data(gi + 310);
    const auto *gi_311 = buffer.data(gi + 311);
    const auto *gi_312 = buffer.data(gi + 312);
    const auto *gi_313 = buffer.data(gi + 313);
    const auto *gi_314 = buffer.data(gi + 314);
    const auto *gi_315 = buffer.data(gi + 315);
    const auto *gi_316 = buffer.data(gi + 316);
    const auto *gi_317 = buffer.data(gi + 317);
    const auto *gi_318 = buffer.data(gi + 318);
    const auto *gi_319 = buffer.data(gi + 319);
    const auto *gi_320 = buffer.data(gi + 320);
    const auto *gi_321 = buffer.data(gi + 321);
    const auto *gi_322 = buffer.data(gi + 322);
    const auto *gi_323 = buffer.data(gi + 323);
    const auto *gi_324 = buffer.data(gi + 324);
    const auto *gi_325 = buffer.data(gi + 325);
    const auto *gi_326 = buffer.data(gi + 326);
    const auto *gi_327 = buffer.data(gi + 327);
    const auto *gi_328 = buffer.data(gi + 328);
    const auto *gi_330 = buffer.data(gi + 330);
    const auto *gi_332 = buffer.data(gi + 332);
    const auto *gi_334 = buffer.data(gi + 334);
    const auto *gi_336 = buffer.data(gi + 336);
    const auto *gi_337 = buffer.data(gi + 337);
    const auto *gi_338 = buffer.data(gi + 338);
    const auto *gi_339 = buffer.data(gi + 339);
    const auto *gi_340 = buffer.data(gi + 340);
    const auto *gi_341 = buffer.data(gi + 341);
    const auto *gi_342 = buffer.data(gi + 342);
    const auto *gi_343 = buffer.data(gi + 343);
    const auto *gi_344 = buffer.data(gi + 344);
    const auto *gi_345 = buffer.data(gi + 345);
    const auto *gi_346 = buffer.data(gi + 346);
    const auto *gi_347 = buffer.data(gi + 347);
    const auto *gi_348 = buffer.data(gi + 348);
    const auto *gi_349 = buffer.data(gi + 349);
    const auto *gi_350 = buffer.data(gi + 350);
    const auto *gi_351 = buffer.data(gi + 351);
    const auto *gi_352 = buffer.data(gi + 352);
    const auto *gi_353 = buffer.data(gi + 353);
    const auto *gi_354 = buffer.data(gi + 354);
    const auto *gi_355 = buffer.data(gi + 355);
    const auto *gi_357 = buffer.data(gi + 357);
    const auto *gi_358 = buffer.data(gi + 358);
    const auto *gi_359 = buffer.data(gi + 359);
    const auto *gi_360 = buffer.data(gi + 360);
    const auto *gi_361 = buffer.data(gi + 361);
    const auto *gi_362 = buffer.data(gi + 362);
    const auto *gi_364 = buffer.data(gi + 364);
    const auto *gi_365 = buffer.data(gi + 365);
    const auto *gi_366 = buffer.data(gi + 366);
    const auto *gi_367 = buffer.data(gi + 367);
    const auto *gi_368 = buffer.data(gi + 368);
    const auto *gi_369 = buffer.data(gi + 369);
    const auto *gi_370 = buffer.data(gi + 370);
    const auto *gi_371 = buffer.data(gi + 371);
    const auto *gi_372 = buffer.data(gi + 372);
    const auto *gi_373 = buffer.data(gi + 373);
    const auto *gi_374 = buffer.data(gi + 374);
    const auto *gi_375 = buffer.data(gi + 375);
    const auto *gi_376 = buffer.data(gi + 376);
    const auto *gi_377 = buffer.data(gi + 377);
    const auto *gi_378 = buffer.data(gi + 378);
    const auto *gi_379 = buffer.data(gi + 379);
    const auto *gi_380 = buffer.data(gi + 380);
    const auto *gi_381 = buffer.data(gi + 381);
    const auto *gi_382 = buffer.data(gi + 382);
    const auto *gi_383 = buffer.data(gi + 383);
    const auto *gi_384 = buffer.data(gi + 384);
    const auto *gi_386 = buffer.data(gi + 386);
    const auto *gi_388 = buffer.data(gi + 388);
    const auto *gi_390 = buffer.data(gi + 390);
    const auto *gi_392 = buffer.data(gi + 392);
    const auto *gi_393 = buffer.data(gi + 393);
    const auto *gi_394 = buffer.data(gi + 394);
    const auto *gi_395 = buffer.data(gi + 395);
    const auto *gi_396 = buffer.data(gi + 396);
    const auto *gi_397 = buffer.data(gi + 397);
    const auto *gi_398 = buffer.data(gi + 398);
    const auto *gi_399 = buffer.data(gi + 399);
    const auto *gi_400 = buffer.data(gi + 400);
    const auto *gi_401 = buffer.data(gi + 401);
    const auto *gi_402 = buffer.data(gi + 402);
    const auto *gi_403 = buffer.data(gi + 403);
    const auto *gi_404 = buffer.data(gi + 404);
    const auto *gi_405 = buffer.data(gi + 405);
    const auto *gi_406 = buffer.data(gi + 406);
    const auto *gi_407 = buffer.data(gi + 407);
    const auto *gi_408 = buffer.data(gi + 408);
    const auto *gi_409 = buffer.data(gi + 409);
    const auto *gi_410 = buffer.data(gi + 410);
    const auto *gi_411 = buffer.data(gi + 411);
    const auto *gi_412 = buffer.data(gi + 412);
    const auto *gi_413 = buffer.data(gi + 413);
    const auto *gi_414 = buffer.data(gi + 414);
    const auto *gi_415 = buffer.data(gi + 415);
    const auto *gi_416 = buffer.data(gi + 416);
    const auto *gi_417 = buffer.data(gi + 417);
    const auto *gi_418 = buffer.data(gi + 418);
    const auto *gi_419 = buffer.data(gi + 419);

#pragma omp simd aligned(ab_x, ab_y, gh_22, gh_27, gh_36, gh_127, gh_132, gh_141, gh_211, \
                         gh_216, gh_225, gi_29, gi_34, gi_43, gi_169, gi_174, gi_183, gi_283, \
                         gi_290, gi_301 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -12.3046875 * ab_x[k] * gh_22[k]
                 + 24.609375 * ab_x[k] * gh_27[k]
                 - 2.4609375 * ab_x[k] * gh_36[k]
                 + 24.609375 * ab_x[k] * gh_127[k]
                 - 49.21875 * ab_x[k] * gh_132[k]
                 + 4.921875 * ab_x[k] * gh_141[k]
                 - 2.4609375 * ab_y[k] * gh_211[k]
                 + 4.921875 * ab_y[k] * gh_216[k]
                 - 0.4921875 * ab_y[k] * gh_225[k]
                 + 12.3046875 * gi_29[k]
                 - 24.609375 * gi_34[k]
                 + 2.4609375 * gi_43[k]
                 - 24.609375 * gi_169[k]
                 + 49.21875 * gi_174[k]
                 - 4.921875 * gi_183[k]
                 + 2.4609375 * gi_283[k]
                 - 4.921875 * gi_290[k]
                 + 0.4921875 * gi_301[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_25, gh_32, gh_130, gh_137, gh_214, gh_221, gi_32, \
                         gi_39, gi_172, gi_179, gi_287, gi_296 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = -f_0 * ab_x[k] * gh_25[k]
                 + f_0 * ab_x[k] * gh_32[k]
                 + f_1 * ab_x[k] * gh_130[k]
                 - f_1 * ab_x[k] * gh_137[k]
                 - f_2 * ab_y[k] * gh_214[k]
                 + f_2 * ab_y[k] * gh_221[k]
                 + f_0 * gi_32[k]
                 - f_0 * gi_39[k]
                 - f_1 * gi_172[k]
                 + f_1 * gi_179[k]
                 + f_2 * gi_287[k]
                 - f_2 * gi_296[k];
        g_11[k] = g_1[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_22, gh_27, gh_29, gh_36, gh_38, gh_127, gh_132, \
                         gh_134, gh_141, gh_143, gh_211, gh_216, gh_218, gh_225, gh_227, \
                         gi_29, gi_34, gi_36, gi_43, gi_45, gi_169, gi_174, gi_176, gi_183, \
                         gi_185, gi_283, gi_290, gi_292, gi_301, \
                         gi_303 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = f_3 * ab_x[k] * gh_22[k]
                 + f_4 * ab_x[k] * gh_27[k]
                 - f_5 * ab_x[k] * gh_29[k]
                 - f_6 * ab_x[k] * gh_36[k]
                 + f_7 * ab_x[k] * gh_38[k]
                 - f_8 * ab_x[k] * gh_127[k]
                 - f_9 * ab_x[k] * gh_132[k]
                 + f_10 * ab_x[k] * gh_134[k]
                 + f_4 * ab_x[k] * gh_141[k]
                 - f_11 * ab_x[k] * gh_143[k]
                 + f_12 * ab_y[k] * gh_211[k]
                 + f_13 * ab_y[k] * gh_216[k]
                 - f_14 * ab_y[k] * gh_218[k]
                 - f_15 * ab_y[k] * gh_225[k]
                 + f_16 * ab_y[k] * gh_227[k]
                 - f_3 * gi_29[k]
                 - f_4 * gi_34[k]
                 + f_5 * gi_36[k]
                 + f_6 * gi_43[k]
                 - f_7 * gi_45[k]
                 + f_8 * gi_169[k]
                 + f_9 * gi_174[k]
                 - f_10 * gi_176[k]
                 - f_4 * gi_183[k]
                 + f_11 * gi_185[k]
                 - f_12 * gi_283[k]
                 - f_13 * gi_290[k]
                 + f_14 * gi_292[k]
                 + f_15 * gi_301[k]
                 - f_16 * gi_303[k];
        g_22[k] = g_2[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_25, gh_32, gh_34, gh_130, gh_137, gh_139, gh_214, \
                         gh_221, gh_223, gi_32, gi_39, gi_41, gi_172, gi_179, gi_181, gi_287, \
                         gi_296, gi_298 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_17 * ab_x[k] * gh_25[k]
                 + f_17 * ab_x[k] * gh_32[k]
                 - f_18 * ab_x[k] * gh_34[k]
                 - f_18 * ab_x[k] * gh_130[k]
                 - f_18 * ab_x[k] * gh_137[k]
                 + f_19 * ab_x[k] * gh_139[k]
                 + f_20 * ab_y[k] * gh_214[k]
                 + f_20 * ab_y[k] * gh_221[k]
                 - f_21 * ab_y[k] * gh_223[k]
                 - f_17 * gi_32[k]
                 - f_17 * gi_39[k]
                 + f_18 * gi_41[k]
                 + f_18 * gi_172[k]
                 + f_18 * gi_179[k]
                 - f_19 * gi_181[k]
                 - f_20 * gi_287[k]
                 - f_20 * gi_296[k]
                 + f_21 * gi_298[k];
        g_33[k] = g_3[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_22, gh_27, gh_29, gh_36, gh_38, gh_40, gh_127, gh_132, \
                         gh_134, gh_141, gh_143, gh_145, gh_211, gh_216, gh_218, gh_225, \
                         gh_227, gh_229, gi_29, gi_34, gi_36, gi_43, gi_45, gi_47, gi_169, \
                         gi_174, gi_176, gi_183, gi_185, gi_187, gi_283, gi_290, gi_292, \
                         gi_301, gi_303, gi_305 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_22 * ab_x[k] * gh_22[k]
                 - f_23 * ab_x[k] * gh_27[k]
                 + f_24 * ab_x[k] * gh_29[k]
                 - f_22 * ab_x[k] * gh_36[k]
                 + f_24 * ab_x[k] * gh_38[k]
                 - f_25 * ab_x[k] * gh_40[k]
                 + f_23 * ab_x[k] * gh_127[k]
                 + f_26 * ab_x[k] * gh_132[k]
                 - f_27 * ab_x[k] * gh_134[k]
                 + f_23 * ab_x[k] * gh_141[k]
                 - f_27 * ab_x[k] * gh_143[k]
                 + f_28 * ab_x[k] * gh_145[k]
                 - f_29 * ab_y[k] * gh_211[k]
                 - f_30 * ab_y[k] * gh_216[k]
                 + f_31 * ab_y[k] * gh_218[k]
                 - f_29 * ab_y[k] * gh_225[k]
                 + f_31 * ab_y[k] * gh_227[k]
                 - f_32 * ab_y[k] * gh_229[k]
                 + f_22 * gi_29[k]
                 + f_23 * gi_34[k]
                 - f_24 * gi_36[k]
                 + f_22 * gi_43[k]
                 - f_24 * gi_45[k]
                 + f_25 * gi_47[k]
                 - f_23 * gi_169[k]
                 - f_26 * gi_174[k]
                 + f_27 * gi_176[k]
                 - f_23 * gi_183[k]
                 + f_27 * gi_185[k]
                 - f_28 * gi_187[k]
                 + f_29 * gi_283[k]
                 + f_30 * gi_290[k]
                 - f_31 * gi_292[k]
                 + f_29 * gi_301[k]
                 - f_31 * gi_303[k]
                 + f_32 * gi_305[k];
        g_44[k] = g_4[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_28, gh_30, gh_37, gh_39, gh_41, gh_128, gh_133, \
                         gh_135, gh_142, gh_144, gh_146, gh_212, gh_217, gh_219, gh_226, \
                         gh_228, gh_230, gi_30, gi_35, gi_37, gi_44, gi_46, gi_48, gi_170, \
                         gi_175, gi_177, gi_184, gi_186, gi_188, gi_284, gi_291, gi_293, \
                         gi_302, gi_304, gi_306 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -f_33 * ab_x[k] * gh_23[k]
                 - f_34 * ab_x[k] * gh_28[k]
                 + f_35 * ab_x[k] * gh_30[k]
                 - f_33 * ab_x[k] * gh_37[k]
                 + f_35 * ab_x[k] * gh_39[k]
                 - f_36 * ab_x[k] * gh_41[k]
                 + f_34 * ab_x[k] * gh_128[k]
                 + f_37 * ab_x[k] * gh_133[k]
                 - f_38 * ab_x[k] * gh_135[k]
                 + f_34 * ab_x[k] * gh_142[k]
                 - f_38 * ab_x[k] * gh_144[k]
                 + f_39 * ab_x[k] * gh_146[k]
                 - f_40 * ab_y[k] * gh_212[k]
                 - f_41 * ab_y[k] * gh_217[k]
                 + f_36 * ab_y[k] * gh_219[k]
                 - f_40 * ab_y[k] * gh_226[k]
                 + f_36 * ab_y[k] * gh_228[k]
                 - f_42 * ab_y[k] * gh_230[k]
                 + f_33 * gi_30[k]
                 + f_34 * gi_35[k]
                 - f_35 * gi_37[k]
                 + f_33 * gi_44[k]
                 - f_35 * gi_46[k]
                 + f_36 * gi_48[k]
                 - f_34 * gi_170[k]
                 - f_37 * gi_175[k]
                 + f_38 * gi_177[k]
                 - f_34 * gi_184[k]
                 + f_38 * gi_186[k]
                 - f_39 * gi_188[k]
                 + f_40 * gi_284[k]
                 + f_41 * gi_291[k]
                 - f_36 * gi_293[k]
                 + f_40 * gi_302[k]
                 - f_36 * gi_304[k]
                 + f_42 * gi_306[k];
        g_55[k] = g_5[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_26, gh_31, gh_33, gh_35, gh_126, gh_129, \
                         gh_131, gh_136, gh_138, gh_140, gh_210, gh_213, gh_215, gh_220, \
                         gh_222, gh_224, gi_28, gi_31, gi_33, gi_38, gi_40, gi_42, gi_168, \
                         gi_171, gi_173, gi_178, gi_180, gi_182, gi_281, gi_286, gi_288, \
                         gi_295, gi_297, gi_299 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_22 * ab_x[k] * gh_21[k]
                 - f_23 * ab_x[k] * gh_24[k]
                 + f_24 * ab_x[k] * gh_26[k]
                 - f_22 * ab_x[k] * gh_31[k]
                 + f_24 * ab_x[k] * gh_33[k]
                 - f_25 * ab_x[k] * gh_35[k]
                 + f_23 * ab_x[k] * gh_126[k]
                 + f_26 * ab_x[k] * gh_129[k]
                 - f_27 * ab_x[k] * gh_131[k]
                 + f_23 * ab_x[k] * gh_136[k]
                 - f_27 * ab_x[k] * gh_138[k]
                 + f_28 * ab_x[k] * gh_140[k]
                 - f_29 * ab_y[k] * gh_210[k]
                 - f_30 * ab_y[k] * gh_213[k]
                 + f_31 * ab_y[k] * gh_215[k]
                 - f_29 * ab_y[k] * gh_220[k]
                 + f_31 * ab_y[k] * gh_222[k]
                 - f_32 * ab_y[k] * gh_224[k]
                 + f_22 * gi_28[k]
                 + f_23 * gi_31[k]
                 - f_24 * gi_33[k]
                 + f_22 * gi_38[k]
                 - f_24 * gi_40[k]
                 + f_25 * gi_42[k]
                 - f_23 * gi_168[k]
                 - f_26 * gi_171[k]
                 + f_27 * gi_173[k]
                 - f_23 * gi_178[k]
                 + f_27 * gi_180[k]
                 - f_28 * gi_182[k]
                 + f_29 * gi_281[k]
                 + f_30 * gi_286[k]
                 - f_31 * gi_288[k]
                 + f_29 * gi_295[k]
                 - f_31 * gi_297[k]
                 + f_32 * gi_299[k];
        g_66[k] = g_6[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_30, gh_37, gh_39, gh_128, gh_135, gh_142, \
                         gh_144, gh_212, gh_219, gh_226, gh_228, gi_30, gi_37, gi_44, gi_46, \
                         gi_170, gi_177, gi_184, gi_186, gi_284, gi_293, gi_302, \
                         gi_304 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_43 * ab_x[k] * gh_23[k]
                 - f_17 * ab_x[k] * gh_30[k]
                 - f_43 * ab_x[k] * gh_37[k]
                 + f_17 * ab_x[k] * gh_39[k]
                 - f_17 * ab_x[k] * gh_128[k]
                 + f_18 * ab_x[k] * gh_135[k]
                 + f_17 * ab_x[k] * gh_142[k]
                 - f_18 * ab_x[k] * gh_144[k]
                 + f_44 * ab_y[k] * gh_212[k]
                 - f_20 * ab_y[k] * gh_219[k]
                 - f_44 * ab_y[k] * gh_226[k]
                 + f_20 * ab_y[k] * gh_228[k]
                 - f_43 * gi_30[k]
                 + f_17 * gi_37[k]
                 + f_43 * gi_44[k]
                 - f_17 * gi_46[k]
                 + f_17 * gi_170[k]
                 - f_18 * gi_177[k]
                 - f_17 * gi_184[k]
                 + f_18 * gi_186[k]
                 - f_44 * gi_284[k]
                 + f_20 * gi_293[k]
                 + f_44 * gi_302[k]
                 - f_20 * gi_304[k];
        g_77[k] = g_7[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_26, gh_31, gh_33, gh_126, gh_129, \
                         gh_131, gh_136, gh_138, gh_210, gh_213, gh_215, gh_220, gh_222, \
                         gi_28, gi_31, gi_33, gi_38, gi_40, gi_168, gi_171, gi_173, gi_178, \
                         gi_180, gi_281, gi_286, gi_288, gi_295, \
                         gi_297 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_6 * ab_x[k] * gh_21[k]
                 - f_4 * ab_x[k] * gh_24[k]
                 - f_7 * ab_x[k] * gh_26[k]
                 - f_3 * ab_x[k] * gh_31[k]
                 + f_5 * ab_x[k] * gh_33[k]
                 - f_4 * ab_x[k] * gh_126[k]
                 + f_9 * ab_x[k] * gh_129[k]
                 + f_11 * ab_x[k] * gh_131[k]
                 + f_8 * ab_x[k] * gh_136[k]
                 - f_10 * ab_x[k] * gh_138[k]
                 + f_15 * ab_y[k] * gh_210[k]
                 - f_13 * ab_y[k] * gh_213[k]
                 - f_16 * ab_y[k] * gh_215[k]
                 - f_12 * ab_y[k] * gh_220[k]
                 + f_14 * ab_y[k] * gh_222[k]
                 - f_6 * gi_28[k]
                 + f_4 * gi_31[k]
                 + f_7 * gi_33[k]
                 + f_3 * gi_38[k]
                 - f_5 * gi_40[k]
                 + f_4 * gi_168[k]
                 - f_9 * gi_171[k]
                 - f_11 * gi_173[k]
                 - f_8 * gi_178[k]
                 + f_10 * gi_180[k]
                 - f_15 * gi_281[k]
                 + f_13 * gi_286[k]
                 + f_16 * gi_288[k]
                 + f_12 * gi_295[k]
                 - f_14 * gi_297[k];
        g_88[k] = g_8[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_28, gh_37, gh_128, gh_133, gh_142, gh_212, \
                         gh_217, gh_226, gi_30, gi_35, gi_44, gi_170, gi_175, gi_184, gi_284, \
                         gi_291, gi_302 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_45 * ab_x[k] * gh_23[k]
                 + f_46 * ab_x[k] * gh_28[k]
                 - f_45 * ab_x[k] * gh_37[k]
                 + f_47 * ab_x[k] * gh_128[k]
                 - f_48 * ab_x[k] * gh_133[k]
                 + f_47 * ab_x[k] * gh_142[k]
                 - f_49 * ab_y[k] * gh_212[k]
                 + f_50 * ab_y[k] * gh_217[k]
                 - f_49 * ab_y[k] * gh_226[k]
                 + f_45 * gi_30[k]
                 - f_46 * gi_35[k]
                 + f_45 * gi_44[k]
                 - f_47 * gi_170[k]
                 + f_48 * gi_175[k]
                 - f_47 * gi_184[k]
                 + f_49 * gi_284[k]
                 - f_50 * gi_291[k]
                 + f_49 * gi_302[k];
        g_99[k] = g_9[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_31, gh_126, gh_129, gh_136, gh_210, \
                         gh_213, gh_220, gi_28, gi_31, gi_38, gi_168, gi_171, gi_178, gi_281, \
                         gi_286, gi_295 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -2.4609375 * ab_x[k] * gh_21[k]
                  + 24.609375 * ab_x[k] * gh_24[k]
                  - 12.3046875 * ab_x[k] * gh_31[k]
                  + 4.921875 * ab_x[k] * gh_126[k]
                  - 49.21875 * ab_x[k] * gh_129[k]
                  + 24.609375 * ab_x[k] * gh_136[k]
                  - 0.4921875 * ab_y[k] * gh_210[k]
                  + 4.921875 * ab_y[k] * gh_213[k]
                  - 2.4609375 * ab_y[k] * gh_220[k]
                  + 2.4609375 * gi_28[k]
                  - 24.609375 * gi_31[k]
                  + 12.3046875 * gi_38[k]
                  - 4.921875 * gi_168[k]
                  + 49.21875 * gi_171[k]
                  - 24.609375 * gi_178[k]
                  + 0.4921875 * gi_281[k]
                  - 4.921875 * gi_286[k]
                  + 2.4609375 * gi_295[k];
        g_110[k] = g_10[k];
    }

#pragma omp simd aligned(ab_x, gh_88, gh_95, gh_235, gh_242, gi_116, gi_123, gi_312, \
                         gi_319 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = -78.75 * ab_x[k] * gh_88[k]
                  + 78.75 * ab_x[k] * gh_95[k]
                  + 78.75 * ab_x[k] * gh_235[k]
                  - 78.75 * ab_x[k] * gh_242[k]
                  + 78.75 * gi_116[k]
                  - 78.75 * gi_123[k]
                  - 78.75 * gi_312[k]
                  + 78.75 * gi_319[k];
    }

#pragma omp simd aligned(ab_x, gh_85, gh_90, gh_92, gh_99, gh_101, gh_232, gh_237, gh_239, \
                         gh_246, gh_248, gi_113, gi_118, gi_120, gi_127, gi_129, gi_309, \
                         gi_314, gi_316, gi_323, gi_325 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_51 * ab_x[k] * gh_85[k]
                  + f_52 * ab_x[k] * gh_90[k]
                  - f_53 * ab_x[k] * gh_92[k]
                  - f_54 * ab_x[k] * gh_99[k]
                  + f_55 * ab_x[k] * gh_101[k]
                  - f_51 * ab_x[k] * gh_232[k]
                  - f_52 * ab_x[k] * gh_237[k]
                  + f_53 * ab_x[k] * gh_239[k]
                  + f_54 * ab_x[k] * gh_246[k]
                  - f_55 * ab_x[k] * gh_248[k]
                  - f_51 * gi_113[k]
                  - f_52 * gi_118[k]
                  + f_53 * gi_120[k]
                  + f_54 * gi_127[k]
                  - f_55 * gi_129[k]
                  + f_51 * gi_309[k]
                  + f_52 * gi_314[k]
                  - f_53 * gi_316[k]
                  - f_54 * gi_323[k]
                  + f_55 * gi_325[k];
        g_23[k] = g_13[k];
    }

#pragma omp simd aligned(ab_x, gh_88, gh_95, gh_97, gh_235, gh_242, gh_244, gi_116, gi_123, \
                         gi_125, gi_312, gi_319, gi_321 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_56 * ab_x[k] * gh_88[k]
                  + f_56 * ab_x[k] * gh_95[k]
                  - f_57 * ab_x[k] * gh_97[k]
                  - f_56 * ab_x[k] * gh_235[k]
                  - f_56 * ab_x[k] * gh_242[k]
                  + f_57 * ab_x[k] * gh_244[k]
                  - f_56 * gi_116[k]
                  - f_56 * gi_123[k]
                  + f_57 * gi_125[k]
                  + f_56 * gi_312[k]
                  + f_56 * gi_319[k]
                  - f_57 * gi_321[k];
        g_34[k] = g_14[k];
    }

#pragma omp simd aligned(ab_x, gh_85, gh_90, gh_92, gh_99, gh_101, gh_103, gh_232, gh_237, \
                         gh_239, gh_246, gh_248, gh_250, gi_113, gi_118, gi_120, gi_127, \
                         gi_129, gi_131, gi_309, gi_314, gi_316, gi_323, gi_325, \
                         gi_327 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_58 * ab_x[k] * gh_85[k]
                  - f_59 * ab_x[k] * gh_90[k]
                  + f_60 * ab_x[k] * gh_92[k]
                  - f_58 * ab_x[k] * gh_99[k]
                  + f_60 * ab_x[k] * gh_101[k]
                  - f_61 * ab_x[k] * gh_103[k]
                  + f_58 * ab_x[k] * gh_232[k]
                  + f_59 * ab_x[k] * gh_237[k]
                  - f_60 * ab_x[k] * gh_239[k]
                  + f_58 * ab_x[k] * gh_246[k]
                  - f_60 * ab_x[k] * gh_248[k]
                  + f_61 * ab_x[k] * gh_250[k]
                  + f_58 * gi_113[k]
                  + f_59 * gi_118[k]
                  - f_60 * gi_120[k]
                  + f_58 * gi_127[k]
                  - f_60 * gi_129[k]
                  + f_61 * gi_131[k]
                  - f_58 * gi_309[k]
                  - f_59 * gi_314[k]
                  + f_60 * gi_316[k]
                  - f_58 * gi_323[k]
                  + f_60 * gi_325[k]
                  - f_61 * gi_327[k];
        g_45[k] = g_15[k];
    }

#pragma omp simd aligned(ab_x, gh_86, gh_91, gh_93, gh_100, gh_102, gh_104, gh_233, gh_238, \
                         gh_240, gh_247, gh_249, gh_251, gi_114, gi_119, gi_121, gi_128, \
                         gi_130, gi_132, gi_310, gi_315, gi_317, gi_324, gi_326, \
                         gi_328 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = -f_62 * ab_x[k] * gh_86[k]
                  - f_63 * ab_x[k] * gh_91[k]
                  + f_64 * ab_x[k] * gh_93[k]
                  - f_62 * ab_x[k] * gh_100[k]
                  + f_64 * ab_x[k] * gh_102[k]
                  - f_65 * ab_x[k] * gh_104[k]
                  + f_62 * ab_x[k] * gh_233[k]
                  + f_63 * ab_x[k] * gh_238[k]
                  - f_64 * ab_x[k] * gh_240[k]
                  + f_62 * ab_x[k] * gh_247[k]
                  - f_64 * ab_x[k] * gh_249[k]
                  + f_65 * ab_x[k] * gh_251[k]
                  + f_62 * gi_114[k]
                  + f_63 * gi_119[k]
                  - f_64 * gi_121[k]
                  + f_62 * gi_128[k]
                  - f_64 * gi_130[k]
                  + f_65 * gi_132[k]
                  - f_62 * gi_310[k]
                  - f_63 * gi_315[k]
                  + f_64 * gi_317[k]
                  - f_62 * gi_324[k]
                  + f_64 * gi_326[k]
                  - f_65 * gi_328[k];
        g_56[k] = g_16[k];
    }

#pragma omp simd aligned(ab_x, gh_84, gh_87, gh_89, gh_94, gh_96, gh_98, gh_231, gh_234, \
                         gh_236, gh_241, gh_243, gh_245, gi_112, gi_115, gi_117, gi_122, \
                         gi_124, gi_126, gi_308, gi_311, gi_313, gi_318, gi_320, \
                         gi_322 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_58 * ab_x[k] * gh_84[k]
                  - f_59 * ab_x[k] * gh_87[k]
                  + f_60 * ab_x[k] * gh_89[k]
                  - f_58 * ab_x[k] * gh_94[k]
                  + f_60 * ab_x[k] * gh_96[k]
                  - f_61 * ab_x[k] * gh_98[k]
                  + f_58 * ab_x[k] * gh_231[k]
                  + f_59 * ab_x[k] * gh_234[k]
                  - f_60 * ab_x[k] * gh_236[k]
                  + f_58 * ab_x[k] * gh_241[k]
                  - f_60 * ab_x[k] * gh_243[k]
                  + f_61 * ab_x[k] * gh_245[k]
                  + f_58 * gi_112[k]
                  + f_59 * gi_115[k]
                  - f_60 * gi_117[k]
                  + f_58 * gi_122[k]
                  - f_60 * gi_124[k]
                  + f_61 * gi_126[k]
                  - f_58 * gi_308[k]
                  - f_59 * gi_311[k]
                  + f_60 * gi_313[k]
                  - f_58 * gi_318[k]
                  + f_60 * gi_320[k]
                  - f_61 * gi_322[k];
        g_67[k] = g_17[k];
    }

#pragma omp simd aligned(ab_x, gh_86, gh_93, gh_100, gh_102, gh_233, gh_240, gh_247, gh_249, \
                         gi_114, gi_121, gi_128, gi_130, gi_310, gi_317, gi_324, \
                         gi_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_66 * ab_x[k] * gh_86[k]
                  - f_56 * ab_x[k] * gh_93[k]
                  - f_66 * ab_x[k] * gh_100[k]
                  + f_56 * ab_x[k] * gh_102[k]
                  - f_66 * ab_x[k] * gh_233[k]
                  + f_56 * ab_x[k] * gh_240[k]
                  + f_66 * ab_x[k] * gh_247[k]
                  - f_56 * ab_x[k] * gh_249[k]
                  - f_66 * gi_114[k]
                  + f_56 * gi_121[k]
                  + f_66 * gi_128[k]
                  - f_56 * gi_130[k]
                  + f_66 * gi_310[k]
                  - f_56 * gi_317[k]
                  - f_66 * gi_324[k]
                  + f_56 * gi_326[k];
        g_78[k] = g_18[k];
    }

#pragma omp simd aligned(ab_x, gh_84, gh_87, gh_89, gh_94, gh_96, gh_231, gh_234, gh_236, \
                         gh_241, gh_243, gi_112, gi_115, gi_117, gi_122, gi_124, gi_308, \
                         gi_311, gi_313, gi_318, gi_320 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_54 * ab_x[k] * gh_84[k]
                  - f_52 * ab_x[k] * gh_87[k]
                  - f_55 * ab_x[k] * gh_89[k]
                  - f_51 * ab_x[k] * gh_94[k]
                  + f_53 * ab_x[k] * gh_96[k]
                  - f_54 * ab_x[k] * gh_231[k]
                  + f_52 * ab_x[k] * gh_234[k]
                  + f_55 * ab_x[k] * gh_236[k]
                  + f_51 * ab_x[k] * gh_241[k]
                  - f_53 * ab_x[k] * gh_243[k]
                  - f_54 * gi_112[k]
                  + f_52 * gi_115[k]
                  + f_55 * gi_117[k]
                  + f_51 * gi_122[k]
                  - f_53 * gi_124[k]
                  + f_54 * gi_308[k]
                  - f_52 * gi_311[k]
                  - f_55 * gi_313[k]
                  - f_51 * gi_318[k]
                  + f_53 * gi_320[k];
        g_89[k] = g_19[k];
    }

#pragma omp simd aligned(ab_x, gh_86, gh_91, gh_100, gh_233, gh_238, gh_247, gi_114, gi_119, \
                         gi_128, gi_310, gi_315, gi_324 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -19.6875 * ab_x[k] * gh_86[k]
                  + 118.125 * ab_x[k] * gh_91[k]
                  - 19.6875 * ab_x[k] * gh_100[k]
                  + 19.6875 * ab_x[k] * gh_233[k]
                  - 118.125 * ab_x[k] * gh_238[k]
                  + 19.6875 * ab_x[k] * gh_247[k]
                  + 19.6875 * gi_114[k]
                  - 118.125 * gi_119[k]
                  + 19.6875 * gi_128[k]
                  - 19.6875 * gi_310[k]
                  + 118.125 * gi_315[k]
                  - 19.6875 * gi_324[k];
        g_100[k] = g_20[k];
    }

#pragma omp simd aligned(ab_x, gh_84, gh_87, gh_94, gh_231, gh_234, gh_241, gi_112, gi_115, \
                         gi_122, gi_308, gi_311, gi_318 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_2 * ab_x[k] * gh_84[k]
                  + f_1 * ab_x[k] * gh_87[k]
                  - f_0 * ab_x[k] * gh_94[k]
                  + f_2 * ab_x[k] * gh_231[k]
                  - f_1 * ab_x[k] * gh_234[k]
                  + f_0 * ab_x[k] * gh_241[k]
                  + f_2 * gi_112[k]
                  - f_1 * gi_115[k]
                  + f_0 * gi_122[k]
                  - f_2 * gi_308[k]
                  + f_1 * gi_311[k]
                  - f_0 * gi_318[k];
        g_111[k] = g_21[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_22, gh_27, gh_29, gh_36, gh_38, gh_127, gh_132, \
                         gh_134, gh_141, gh_143, gh_169, gh_174, gh_176, gh_183, gh_185, \
                         gh_211, gh_216, gh_218, gh_225, gh_227, gh_253, gh_258, gh_260, \
                         gh_267, gh_269, gi_29, gi_34, gi_36, gi_43, gi_45, gi_169, gi_174, \
                         gi_176, gi_183, gi_185, gi_225, gi_230, gi_232, gi_239, gi_241, \
                         gi_283, gi_290, gi_292, gi_301, gi_303, gi_339, gi_346, gi_348, \
                         gi_357, gi_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -2.4609375 * ab_x[k] * gh_22[k]
                  - 1.640625 * ab_x[k] * gh_27[k]
                  + 19.6875 * ab_x[k] * gh_29[k]
                  + 0.8203125 * ab_x[k] * gh_36[k]
                  - 6.5625 * ab_x[k] * gh_38[k]
                  - 1.640625 * ab_x[k] * gh_127[k]
                  - 1.09375 * ab_x[k] * gh_132[k]
                  + 13.125 * ab_x[k] * gh_134[k]
                  + 0.546875 * ab_x[k] * gh_141[k]
                  - 4.375 * ab_x[k] * gh_143[k]
                  + 19.6875 * ab_x[k] * gh_169[k]
                  + 13.125 * ab_x[k] * gh_174[k]
                  - 157.5 * ab_x[k] * gh_176[k]
                  - 6.5625 * ab_x[k] * gh_183[k]
                  + 52.5 * ab_x[k] * gh_185[k]
                  + 0.8203125 * ab_y[k] * gh_211[k]
                  + 0.546875 * ab_y[k] * gh_216[k]
                  - 6.5625 * ab_y[k] * gh_218[k]
                  - 0.2734375 * ab_y[k] * gh_225[k]
                  + 2.1875 * ab_y[k] * gh_227[k]
                  - 6.5625 * ab_y[k] * gh_253[k]
                  - 4.375 * ab_y[k] * gh_258[k]
                  + 52.5 * ab_y[k] * gh_260[k]
                  + 2.1875 * ab_y[k] * gh_267[k]
                  - 17.5 * ab_y[k] * gh_269[k]
                  + 2.4609375 * gi_29[k]
                  + 1.640625 * gi_34[k]
                  - 19.6875 * gi_36[k]
                  - 0.8203125 * gi_43[k]
                  + 6.5625 * gi_45[k]
                  + 1.640625 * gi_169[k]
                  + 1.09375 * gi_174[k]
                  - 13.125 * gi_176[k]
                  - 0.546875 * gi_183[k]
                  + 4.375 * gi_185[k]
                  - 19.6875 * gi_225[k]
                  - 13.125 * gi_230[k]
                  + 157.5 * gi_232[k]
                  + 6.5625 * gi_239[k]
                  - 52.5 * gi_241[k]
                  - 0.8203125 * gi_283[k]
                  - 0.546875 * gi_290[k]
                  + 6.5625 * gi_292[k]
                  + 0.2734375 * gi_301[k]
                  - 2.1875 * gi_303[k]
                  + 6.5625 * gi_339[k]
                  + 4.375 * gi_346[k]
                  - 52.5 * gi_348[k]
                  - 2.1875 * gi_357[k]
                  + 17.5 * gi_359[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_25, gh_32, gh_34, gh_130, gh_137, gh_139, gh_172, \
                         gh_179, gh_181, gh_214, gh_221, gh_223, gh_256, gh_263, gh_265, \
                         gi_32, gi_39, gi_41, gi_172, gi_179, gi_181, gi_228, gi_235, gi_237, \
                         gi_287, gi_296, gi_298, gi_343, gi_352, \
                         gi_354 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_67 * ab_x[k] * gh_25[k]
                  - f_67 * ab_x[k] * gh_32[k]
                  + f_68 * ab_x[k] * gh_34[k]
                  - f_69 * ab_x[k] * gh_130[k]
                  - f_69 * ab_x[k] * gh_137[k]
                  + f_70 * ab_x[k] * gh_139[k]
                  + f_71 * ab_x[k] * gh_172[k]
                  + f_71 * ab_x[k] * gh_179[k]
                  - f_72 * ab_x[k] * gh_181[k]
                  + f_73 * ab_y[k] * gh_214[k]
                  + f_73 * ab_y[k] * gh_221[k]
                  - f_69 * ab_y[k] * gh_223[k]
                  - f_74 * ab_y[k] * gh_256[k]
                  - f_74 * ab_y[k] * gh_263[k]
                  + f_75 * ab_y[k] * gh_265[k]
                  + f_67 * gi_32[k]
                  + f_67 * gi_39[k]
                  - f_68 * gi_41[k]
                  + f_69 * gi_172[k]
                  + f_69 * gi_179[k]
                  - f_70 * gi_181[k]
                  - f_71 * gi_228[k]
                  - f_71 * gi_235[k]
                  + f_72 * gi_237[k]
                  - f_73 * gi_287[k]
                  - f_73 * gi_296[k]
                  + f_69 * gi_298[k]
                  + f_74 * gi_343[k]
                  + f_74 * gi_352[k]
                  - f_75 * gi_354[k];
        g_35[k] = g_25[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_22, gh_27, gh_29, gh_36, gh_38, gh_40, gh_127, gh_132, \
                         gh_134, gh_141, gh_143, gh_145, gh_169, gh_174, gh_176, gh_183, \
                         gh_185, gh_187, gh_211, gh_216, gh_218, gh_225, gh_227, gh_229, \
                         gh_253, gh_258, gh_260, gh_267, gh_269, gh_271, gi_29, gi_34, gi_36, \
                         gi_43, gi_45, gi_47, gi_169, gi_174, gi_176, gi_183, gi_185, gi_187, \
                         gi_225, gi_230, gi_232, gi_239, gi_241, gi_243, gi_283, gi_290, \
                         gi_292, gi_301, gi_303, gi_305, gi_339, gi_346, gi_348, gi_357, \
                         gi_359, gi_361 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = f_76 * ab_x[k] * gh_22[k]
                  + f_77 * ab_x[k] * gh_27[k]
                  - f_78 * ab_x[k] * gh_29[k]
                  + f_76 * ab_x[k] * gh_36[k]
                  - f_78 * ab_x[k] * gh_38[k]
                  + f_79 * ab_x[k] * gh_40[k]
                  + f_80 * ab_x[k] * gh_127[k]
                  + f_81 * ab_x[k] * gh_132[k]
                  - f_79 * ab_x[k] * gh_134[k]
                  + f_80 * ab_x[k] * gh_141[k]
                  - f_79 * ab_x[k] * gh_143[k]
                  + f_82 * ab_x[k] * gh_145[k]
                  - f_79 * ab_x[k] * gh_169[k]
                  - f_83 * ab_x[k] * gh_174[k]
                  + f_84 * ab_x[k] * gh_176[k]
                  - f_79 * ab_x[k] * gh_183[k]
                  + f_84 * ab_x[k] * gh_185[k]
                  - f_85 * ab_x[k] * gh_187[k]
                  - f_86 * ab_y[k] * gh_211[k]
                  - f_80 * ab_y[k] * gh_216[k]
                  + f_87 * ab_y[k] * gh_218[k]
                  - f_86 * ab_y[k] * gh_225[k]
                  + f_87 * ab_y[k] * gh_227[k]
                  - f_88 * ab_y[k] * gh_229[k]
                  + f_88 * ab_y[k] * gh_253[k]
                  + f_82 * ab_y[k] * gh_258[k]
                  - f_89 * ab_y[k] * gh_260[k]
                  + f_88 * ab_y[k] * gh_267[k]
                  - f_89 * ab_y[k] * gh_269[k]
                  + f_90 * ab_y[k] * gh_271[k]
                  - f_76 * gi_29[k]
                  - f_77 * gi_34[k]
                  + f_78 * gi_36[k]
                  - f_76 * gi_43[k]
                  + f_78 * gi_45[k]
                  - f_79 * gi_47[k]
                  - f_80 * gi_169[k]
                  - f_81 * gi_174[k]
                  + f_79 * gi_176[k]
                  - f_80 * gi_183[k]
                  + f_79 * gi_185[k]
                  - f_82 * gi_187[k]
                  + f_79 * gi_225[k]
                  + f_83 * gi_230[k]
                  - f_84 * gi_232[k]
                  + f_79 * gi_239[k]
                  - f_84 * gi_241[k]
                  + f_85 * gi_243[k]
                  + f_86 * gi_283[k]
                  + f_80 * gi_290[k]
                  - f_87 * gi_292[k]
                  + f_86 * gi_301[k]
                  - f_87 * gi_303[k]
                  + f_88 * gi_305[k]
                  - f_88 * gi_339[k]
                  - f_82 * gi_346[k]
                  + f_89 * gi_348[k]
                  - f_88 * gi_357[k]
                  + f_89 * gi_359[k]
                  - f_90 * gi_361[k];
        g_46[k] = g_26[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_28, gh_30, gh_37, gh_39, gh_41, gh_128, gh_133, \
                         gh_135, gh_142, gh_144, gh_146, gh_170, gh_175, gh_177, gh_184, \
                         gh_186, gh_188, gh_212, gh_217, gh_219, gh_226, gh_228, gh_230, \
                         gh_254, gh_259, gh_261, gh_268, gh_270, gh_272, gi_30, gi_35, gi_37, \
                         gi_44, gi_46, gi_48, gi_170, gi_175, gi_177, gi_184, gi_186, gi_188, \
                         gi_226, gi_231, gi_233, gi_240, gi_242, gi_244, gi_284, gi_291, \
                         gi_293, gi_302, gi_304, gi_306, gi_340, gi_347, gi_349, gi_358, \
                         gi_360, gi_362 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = f_91 * ab_x[k] * gh_23[k]
                  + f_92 * ab_x[k] * gh_28[k]
                  - f_93 * ab_x[k] * gh_30[k]
                  + f_91 * ab_x[k] * gh_37[k]
                  - f_93 * ab_x[k] * gh_39[k]
                  + f_94 * ab_x[k] * gh_41[k]
                  + f_95 * ab_x[k] * gh_128[k]
                  + f_96 * ab_x[k] * gh_133[k]
                  - f_97 * ab_x[k] * gh_135[k]
                  + f_95 * ab_x[k] * gh_142[k]
                  - f_97 * ab_x[k] * gh_144[k]
                  + f_98 * ab_x[k] * gh_146[k]
                  - f_99 * ab_x[k] * gh_170[k]
                  - f_100 * ab_x[k] * gh_175[k]
                  + f_101 * ab_x[k] * gh_177[k]
                  - f_99 * ab_x[k] * gh_184[k]
                  + f_101 * ab_x[k] * gh_186[k]
                  - f_102 * ab_x[k] * gh_188[k]
                  - f_103 * ab_y[k] * gh_212[k]
                  - f_95 * ab_y[k] * gh_217[k]
                  + f_104 * ab_y[k] * gh_219[k]
                  - f_103 * ab_y[k] * gh_226[k]
                  + f_104 * ab_y[k] * gh_228[k]
                  - f_105 * ab_y[k] * gh_230[k]
                  + f_93 * ab_y[k] * gh_254[k]
                  + f_106 * ab_y[k] * gh_259[k]
                  - f_107 * ab_y[k] * gh_261[k]
                  + f_93 * ab_y[k] * gh_268[k]
                  - f_107 * ab_y[k] * gh_270[k]
                  + f_108 * ab_y[k] * gh_272[k]
                  - f_91 * gi_30[k]
                  - f_92 * gi_35[k]
                  + f_93 * gi_37[k]
                  - f_91 * gi_44[k]
                  + f_93 * gi_46[k]
                  - f_94 * gi_48[k]
                  - f_95 * gi_170[k]
                  - f_96 * gi_175[k]
                  + f_97 * gi_177[k]
                  - f_95 * gi_184[k]
                  + f_97 * gi_186[k]
                  - f_98 * gi_188[k]
                  + f_99 * gi_226[k]
                  + f_100 * gi_231[k]
                  - f_101 * gi_233[k]
                  + f_99 * gi_240[k]
                  - f_101 * gi_242[k]
                  + f_102 * gi_244[k]
                  + f_103 * gi_284[k]
                  + f_95 * gi_291[k]
                  - f_104 * gi_293[k]
                  + f_103 * gi_302[k]
                  - f_104 * gi_304[k]
                  + f_105 * gi_306[k]
                  - f_93 * gi_340[k]
                  - f_106 * gi_347[k]
                  + f_107 * gi_349[k]
                  - f_93 * gi_358[k]
                  + f_107 * gi_360[k]
                  - f_108 * gi_362[k];
        g_57[k] = g_27[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_26, gh_31, gh_33, gh_35, gh_126, gh_129, \
                         gh_131, gh_136, gh_138, gh_140, gh_168, gh_171, gh_173, gh_178, \
                         gh_180, gh_182, gh_210, gh_213, gh_215, gh_220, gh_222, gh_224, \
                         gh_252, gh_255, gh_257, gh_262, gh_264, gh_266, gi_28, gi_31, gi_33, \
                         gi_38, gi_40, gi_42, gi_168, gi_171, gi_173, gi_178, gi_180, gi_182, \
                         gi_224, gi_227, gi_229, gi_234, gi_236, gi_238, gi_281, gi_286, \
                         gi_288, gi_295, gi_297, gi_299, gi_337, gi_342, gi_344, gi_351, \
                         gi_353, gi_355 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_76 * ab_x[k] * gh_21[k]
                  + f_77 * ab_x[k] * gh_24[k]
                  - f_78 * ab_x[k] * gh_26[k]
                  + f_76 * ab_x[k] * gh_31[k]
                  - f_78 * ab_x[k] * gh_33[k]
                  + f_79 * ab_x[k] * gh_35[k]
                  + f_80 * ab_x[k] * gh_126[k]
                  + f_81 * ab_x[k] * gh_129[k]
                  - f_79 * ab_x[k] * gh_131[k]
                  + f_80 * ab_x[k] * gh_136[k]
                  - f_79 * ab_x[k] * gh_138[k]
                  + f_82 * ab_x[k] * gh_140[k]
                  - f_79 * ab_x[k] * gh_168[k]
                  - f_83 * ab_x[k] * gh_171[k]
                  + f_84 * ab_x[k] * gh_173[k]
                  - f_79 * ab_x[k] * gh_178[k]
                  + f_84 * ab_x[k] * gh_180[k]
                  - f_85 * ab_x[k] * gh_182[k]
                  - f_86 * ab_y[k] * gh_210[k]
                  - f_80 * ab_y[k] * gh_213[k]
                  + f_87 * ab_y[k] * gh_215[k]
                  - f_86 * ab_y[k] * gh_220[k]
                  + f_87 * ab_y[k] * gh_222[k]
                  - f_88 * ab_y[k] * gh_224[k]
                  + f_88 * ab_y[k] * gh_252[k]
                  + f_82 * ab_y[k] * gh_255[k]
                  - f_89 * ab_y[k] * gh_257[k]
                  + f_88 * ab_y[k] * gh_262[k]
                  - f_89 * ab_y[k] * gh_264[k]
                  + f_90 * ab_y[k] * gh_266[k]
                  - f_76 * gi_28[k]
                  - f_77 * gi_31[k]
                  + f_78 * gi_33[k]
                  - f_76 * gi_38[k]
                  + f_78 * gi_40[k]
                  - f_79 * gi_42[k]
                  - f_80 * gi_168[k]
                  - f_81 * gi_171[k]
                  + f_79 * gi_173[k]
                  - f_80 * gi_178[k]
                  + f_79 * gi_180[k]
                  - f_82 * gi_182[k]
                  + f_79 * gi_224[k]
                  + f_83 * gi_227[k]
                  - f_84 * gi_229[k]
                  + f_79 * gi_234[k]
                  - f_84 * gi_236[k]
                  + f_85 * gi_238[k]
                  + f_86 * gi_281[k]
                  + f_80 * gi_286[k]
                  - f_87 * gi_288[k]
                  + f_86 * gi_295[k]
                  - f_87 * gi_297[k]
                  + f_88 * gi_299[k]
                  - f_88 * gi_337[k]
                  - f_82 * gi_342[k]
                  + f_89 * gi_344[k]
                  - f_88 * gi_351[k]
                  + f_89 * gi_353[k]
                  - f_90 * gi_355[k];
        g_68[k] = g_28[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_30, gh_37, gh_39, gh_128, gh_135, gh_142, \
                         gh_144, gh_170, gh_177, gh_184, gh_186, gh_212, gh_219, gh_226, \
                         gh_228, gh_254, gh_261, gh_268, gh_270, gi_30, gi_37, gi_44, gi_46, \
                         gi_170, gi_177, gi_184, gi_186, gi_226, gi_233, gi_240, gi_242, \
                         gi_284, gi_293, gi_302, gi_304, gi_340, gi_349, gi_358, \
                         gi_360 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = -f_109 * ab_x[k] * gh_23[k]
                  + f_67 * ab_x[k] * gh_30[k]
                  + f_109 * ab_x[k] * gh_37[k]
                  - f_67 * ab_x[k] * gh_39[k]
                  - f_73 * ab_x[k] * gh_128[k]
                  + f_69 * ab_x[k] * gh_135[k]
                  + f_73 * ab_x[k] * gh_142[k]
                  - f_69 * ab_x[k] * gh_144[k]
                  + f_110 * ab_x[k] * gh_170[k]
                  - f_71 * ab_x[k] * gh_177[k]
                  - f_110 * ab_x[k] * gh_184[k]
                  + f_71 * ab_x[k] * gh_186[k]
                  + f_111 * ab_y[k] * gh_212[k]
                  - f_73 * ab_y[k] * gh_219[k]
                  - f_111 * ab_y[k] * gh_226[k]
                  + f_73 * ab_y[k] * gh_228[k]
                  - f_70 * ab_y[k] * gh_254[k]
                  + f_74 * ab_y[k] * gh_261[k]
                  + f_70 * ab_y[k] * gh_268[k]
                  - f_74 * ab_y[k] * gh_270[k]
                  + f_109 * gi_30[k]
                  - f_67 * gi_37[k]
                  - f_109 * gi_44[k]
                  + f_67 * gi_46[k]
                  + f_73 * gi_170[k]
                  - f_69 * gi_177[k]
                  - f_73 * gi_184[k]
                  + f_69 * gi_186[k]
                  - f_110 * gi_226[k]
                  + f_71 * gi_233[k]
                  + f_110 * gi_240[k]
                  - f_71 * gi_242[k]
                  - f_111 * gi_284[k]
                  + f_73 * gi_293[k]
                  + f_111 * gi_302[k]
                  - f_73 * gi_304[k]
                  + f_70 * gi_340[k]
                  - f_74 * gi_349[k]
                  - f_70 * gi_358[k]
                  + f_74 * gi_360[k];
        g_79[k] = g_29[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_26, gh_31, gh_33, gh_126, gh_129, \
                         gh_131, gh_136, gh_138, gh_168, gh_171, gh_173, gh_178, gh_180, \
                         gh_210, gh_213, gh_215, gh_220, gh_222, gh_252, gh_255, gh_257, \
                         gh_262, gh_264, gi_28, gi_31, gi_33, gi_38, gi_40, gi_168, gi_171, \
                         gi_173, gi_178, gi_180, gi_224, gi_227, gi_229, gi_234, gi_236, \
                         gi_281, gi_286, gi_288, gi_295, gi_297, gi_337, gi_342, gi_344, \
                         gi_351, gi_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -0.8203125 * ab_x[k] * gh_21[k]
                  + 1.640625 * ab_x[k] * gh_24[k]
                  + 6.5625 * ab_x[k] * gh_26[k]
                  + 2.4609375 * ab_x[k] * gh_31[k]
                  - 19.6875 * ab_x[k] * gh_33[k]
                  - 0.546875 * ab_x[k] * gh_126[k]
                  + 1.09375 * ab_x[k] * gh_129[k]
                  + 4.375 * ab_x[k] * gh_131[k]
                  + 1.640625 * ab_x[k] * gh_136[k]
                  - 13.125 * ab_x[k] * gh_138[k]
                  + 6.5625 * ab_x[k] * gh_168[k]
                  - 13.125 * ab_x[k] * gh_171[k]
                  - 52.5 * ab_x[k] * gh_173[k]
                  - 19.6875 * ab_x[k] * gh_178[k]
                  + 157.5 * ab_x[k] * gh_180[k]
                  + 0.2734375 * ab_y[k] * gh_210[k]
                  - 0.546875 * ab_y[k] * gh_213[k]
                  - 2.1875 * ab_y[k] * gh_215[k]
                  - 0.8203125 * ab_y[k] * gh_220[k]
                  + 6.5625 * ab_y[k] * gh_222[k]
                  - 2.1875 * ab_y[k] * gh_252[k]
                  + 4.375 * ab_y[k] * gh_255[k]
                  + 17.5 * ab_y[k] * gh_257[k]
                  + 6.5625 * ab_y[k] * gh_262[k]
                  - 52.5 * ab_y[k] * gh_264[k]
                  + 0.8203125 * gi_28[k]
                  - 1.640625 * gi_31[k]
                  - 6.5625 * gi_33[k]
                  - 2.4609375 * gi_38[k]
                  + 19.6875 * gi_40[k]
                  + 0.546875 * gi_168[k]
                  - 1.09375 * gi_171[k]
                  - 4.375 * gi_173[k]
                  - 1.640625 * gi_178[k]
                  + 13.125 * gi_180[k]
                  - 6.5625 * gi_224[k]
                  + 13.125 * gi_227[k]
                  + 52.5 * gi_229[k]
                  + 19.6875 * gi_234[k]
                  - 157.5 * gi_236[k]
                  - 0.2734375 * gi_281[k]
                  + 0.546875 * gi_286[k]
                  + 2.1875 * gi_288[k]
                  + 0.8203125 * gi_295[k]
                  - 6.5625 * gi_297[k]
                  + 2.1875 * gi_337[k]
                  - 4.375 * gi_342[k]
                  - 17.5 * gi_344[k]
                  - 6.5625 * gi_351[k]
                  + 52.5 * gi_353[k];
        g_90[k] = g_30[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_28, gh_37, gh_128, gh_133, gh_142, gh_170, \
                         gh_175, gh_184, gh_212, gh_217, gh_226, gh_254, gh_259, gh_268, \
                         gi_30, gi_35, gi_44, gi_170, gi_175, gi_184, gi_226, gi_231, gi_240, \
                         gi_284, gi_291, gi_302, gi_340, gi_347, \
                         gi_358 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = f_112 * ab_x[k] * gh_23[k]
                  - f_113 * ab_x[k] * gh_28[k]
                  + f_112 * ab_x[k] * gh_37[k]
                  + f_114 * ab_x[k] * gh_128[k]
                  - f_51 * ab_x[k] * gh_133[k]
                  + f_114 * ab_x[k] * gh_142[k]
                  - f_115 * ab_x[k] * gh_170[k]
                  + f_116 * ab_x[k] * gh_175[k]
                  - f_115 * ab_x[k] * gh_184[k]
                  - f_117 * ab_y[k] * gh_212[k]
                  + f_118 * ab_y[k] * gh_217[k]
                  - f_117 * ab_y[k] * gh_226[k]
                  + f_52 * ab_y[k] * gh_254[k]
                  - f_119 * ab_y[k] * gh_259[k]
                  + f_52 * ab_y[k] * gh_268[k]
                  - f_112 * gi_30[k]
                  + f_113 * gi_35[k]
                  - f_112 * gi_44[k]
                  - f_114 * gi_170[k]
                  + f_51 * gi_175[k]
                  - f_114 * gi_184[k]
                  + f_115 * gi_226[k]
                  - f_116 * gi_231[k]
                  + f_115 * gi_240[k]
                  + f_117 * gi_284[k]
                  - f_118 * gi_291[k]
                  + f_117 * gi_302[k]
                  - f_52 * gi_340[k]
                  + f_119 * gi_347[k]
                  - f_52 * gi_358[k];
        g_101[k] = g_31[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_31, gh_126, gh_129, gh_136, gh_168, \
                         gh_171, gh_178, gh_210, gh_213, gh_220, gh_252, gh_255, gh_262, \
                         gi_28, gi_31, gi_38, gi_168, gi_171, gi_178, gi_224, gi_227, gi_234, \
                         gi_281, gi_286, gi_295, gi_337, gi_342, \
                         gi_351 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_12 * ab_x[k] * gh_21[k]
                  - f_8 * ab_x[k] * gh_24[k]
                  + f_3 * ab_x[k] * gh_31[k]
                  + f_13 * ab_x[k] * gh_126[k]
                  - f_9 * ab_x[k] * gh_129[k]
                  + f_4 * ab_x[k] * gh_136[k]
                  - f_14 * ab_x[k] * gh_168[k]
                  + f_10 * ab_x[k] * gh_171[k]
                  - f_5 * ab_x[k] * gh_178[k]
                  - f_15 * ab_y[k] * gh_210[k]
                  + f_4 * ab_y[k] * gh_213[k]
                  - f_6 * ab_y[k] * gh_220[k]
                  + f_16 * ab_y[k] * gh_252[k]
                  - f_11 * ab_y[k] * gh_255[k]
                  + f_7 * ab_y[k] * gh_262[k]
                  - f_12 * gi_28[k]
                  + f_8 * gi_31[k]
                  - f_3 * gi_38[k]
                  - f_13 * gi_168[k]
                  + f_9 * gi_171[k]
                  - f_4 * gi_178[k]
                  + f_14 * gi_224[k]
                  - f_10 * gi_227[k]
                  + f_5 * gi_234[k]
                  + f_15 * gi_281[k]
                  - f_4 * gi_286[k]
                  + f_6 * gi_295[k]
                  - f_16 * gi_337[k]
                  + f_11 * gi_342[k]
                  - f_7 * gi_351[k];
        g_112[k] = g_32[k];
    }

#pragma omp simd aligned(ab_x, gh_88, gh_95, gh_97, gh_235, gh_242, gh_244, gh_277, gh_284, \
                         gh_286, gi_116, gi_123, gi_125, gi_312, gi_319, gi_321, gi_368, \
                         gi_375, gi_377 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = -26.25 * ab_x[k] * gh_88[k]
                  - 26.25 * ab_x[k] * gh_95[k]
                  + 52.5 * ab_x[k] * gh_97[k]
                  - 26.25 * ab_x[k] * gh_235[k]
                  - 26.25 * ab_x[k] * gh_242[k]
                  + 52.5 * ab_x[k] * gh_244[k]
                  + 52.5 * ab_x[k] * gh_277[k]
                  + 52.5 * ab_x[k] * gh_284[k]
                  - 105.0 * ab_x[k] * gh_286[k]
                  + 26.25 * gi_116[k]
                  + 26.25 * gi_123[k]
                  - 52.5 * gi_125[k]
                  + 26.25 * gi_312[k]
                  + 26.25 * gi_319[k]
                  - 52.5 * gi_321[k]
                  - 52.5 * gi_368[k]
                  - 52.5 * gi_375[k]
                  + 105.0 * gi_377[k];
    }

#pragma omp simd aligned(ab_x, gh_85, gh_90, gh_92, gh_99, gh_101, gh_103, gh_232, gh_237, \
                         gh_239, gh_246, gh_248, gh_250, gh_274, gh_279, gh_281, gh_288, \
                         gh_290, gh_292, gi_113, gi_118, gi_120, gi_127, gi_129, gi_131, \
                         gi_309, gi_314, gi_316, gi_323, gi_325, gi_327, gi_365, gi_370, \
                         gi_372, gi_379, gi_381, gi_383 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_120 * ab_x[k] * gh_85[k]
                  + f_121 * ab_x[k] * gh_90[k]
                  - f_122 * ab_x[k] * gh_92[k]
                  + f_120 * ab_x[k] * gh_99[k]
                  - f_122 * ab_x[k] * gh_101[k]
                  + f_123 * ab_x[k] * gh_103[k]
                  + f_120 * ab_x[k] * gh_232[k]
                  + f_121 * ab_x[k] * gh_237[k]
                  - f_122 * ab_x[k] * gh_239[k]
                  + f_120 * ab_x[k] * gh_246[k]
                  - f_122 * ab_x[k] * gh_248[k]
                  + f_123 * ab_x[k] * gh_250[k]
                  - f_121 * ab_x[k] * gh_274[k]
                  - f_124 * ab_x[k] * gh_279[k]
                  + f_125 * ab_x[k] * gh_281[k]
                  - f_121 * ab_x[k] * gh_288[k]
                  + f_125 * ab_x[k] * gh_290[k]
                  - f_126 * ab_x[k] * gh_292[k]
                  - f_120 * gi_113[k]
                  - f_121 * gi_118[k]
                  + f_122 * gi_120[k]
                  - f_120 * gi_127[k]
                  + f_122 * gi_129[k]
                  - f_123 * gi_131[k]
                  - f_120 * gi_309[k]
                  - f_121 * gi_314[k]
                  + f_122 * gi_316[k]
                  - f_120 * gi_323[k]
                  + f_122 * gi_325[k]
                  - f_123 * gi_327[k]
                  + f_121 * gi_365[k]
                  + f_124 * gi_370[k]
                  - f_125 * gi_372[k]
                  + f_121 * gi_379[k]
                  - f_125 * gi_381[k]
                  + f_126 * gi_383[k];
        g_47[k] = g_37[k];
    }

#pragma omp simd aligned(ab_x, gh_86, gh_91, gh_93, gh_100, gh_102, gh_104, gh_233, gh_238, \
                         gh_240, gh_247, gh_249, gh_251, gh_275, gh_280, gh_282, gh_289, \
                         gh_291, gh_293, gi_114, gi_119, gi_121, gi_128, gi_130, gi_132, \
                         gi_310, gi_315, gi_317, gi_324, gi_326, gi_328, gi_366, gi_371, \
                         gi_373, gi_380, gi_382, gi_384 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_127 * ab_x[k] * gh_86[k]
                  + f_128 * ab_x[k] * gh_91[k]
                  - f_129 * ab_x[k] * gh_93[k]
                  + f_127 * ab_x[k] * gh_100[k]
                  - f_129 * ab_x[k] * gh_102[k]
                  + f_130 * ab_x[k] * gh_104[k]
                  + f_127 * ab_x[k] * gh_233[k]
                  + f_128 * ab_x[k] * gh_238[k]
                  - f_129 * ab_x[k] * gh_240[k]
                  + f_127 * ab_x[k] * gh_247[k]
                  - f_129 * ab_x[k] * gh_249[k]
                  + f_130 * ab_x[k] * gh_251[k]
                  - f_128 * ab_x[k] * gh_275[k]
                  - f_131 * ab_x[k] * gh_280[k]
                  + f_132 * ab_x[k] * gh_282[k]
                  - f_128 * ab_x[k] * gh_289[k]
                  + f_132 * ab_x[k] * gh_291[k]
                  - f_133 * ab_x[k] * gh_293[k]
                  - f_127 * gi_114[k]
                  - f_128 * gi_119[k]
                  + f_129 * gi_121[k]
                  - f_127 * gi_128[k]
                  + f_129 * gi_130[k]
                  - f_130 * gi_132[k]
                  - f_127 * gi_310[k]
                  - f_128 * gi_315[k]
                  + f_129 * gi_317[k]
                  - f_127 * gi_324[k]
                  + f_129 * gi_326[k]
                  - f_130 * gi_328[k]
                  + f_128 * gi_366[k]
                  + f_131 * gi_371[k]
                  - f_132 * gi_373[k]
                  + f_128 * gi_380[k]
                  - f_132 * gi_382[k]
                  + f_133 * gi_384[k];
        g_58[k] = g_38[k];
    }

#pragma omp simd aligned(ab_x, gh_84, gh_87, gh_89, gh_94, gh_96, gh_98, gh_231, gh_234, \
                         gh_236, gh_241, gh_243, gh_245, gh_273, gh_276, gh_278, gh_283, \
                         gh_285, gh_287, gi_112, gi_115, gi_117, gi_122, gi_124, gi_126, \
                         gi_308, gi_311, gi_313, gi_318, gi_320, gi_322, gi_364, gi_367, \
                         gi_369, gi_374, gi_376, gi_378 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = f_120 * ab_x[k] * gh_84[k]
                  + f_121 * ab_x[k] * gh_87[k]
                  - f_122 * ab_x[k] * gh_89[k]
                  + f_120 * ab_x[k] * gh_94[k]
                  - f_122 * ab_x[k] * gh_96[k]
                  + f_123 * ab_x[k] * gh_98[k]
                  + f_120 * ab_x[k] * gh_231[k]
                  + f_121 * ab_x[k] * gh_234[k]
                  - f_122 * ab_x[k] * gh_236[k]
                  + f_120 * ab_x[k] * gh_241[k]
                  - f_122 * ab_x[k] * gh_243[k]
                  + f_123 * ab_x[k] * gh_245[k]
                  - f_121 * ab_x[k] * gh_273[k]
                  - f_124 * ab_x[k] * gh_276[k]
                  + f_125 * ab_x[k] * gh_278[k]
                  - f_121 * ab_x[k] * gh_283[k]
                  + f_125 * ab_x[k] * gh_285[k]
                  - f_126 * ab_x[k] * gh_287[k]
                  - f_120 * gi_112[k]
                  - f_121 * gi_115[k]
                  + f_122 * gi_117[k]
                  - f_120 * gi_122[k]
                  + f_122 * gi_124[k]
                  - f_123 * gi_126[k]
                  - f_120 * gi_308[k]
                  - f_121 * gi_311[k]
                  + f_122 * gi_313[k]
                  - f_120 * gi_318[k]
                  + f_122 * gi_320[k]
                  - f_123 * gi_322[k]
                  + f_121 * gi_364[k]
                  + f_124 * gi_367[k]
                  - f_125 * gi_369[k]
                  + f_121 * gi_374[k]
                  - f_125 * gi_376[k]
                  + f_126 * gi_378[k];
        g_69[k] = g_39[k];
    }

#pragma omp simd aligned(ab_x, gh_86, gh_93, gh_100, gh_102, gh_233, gh_240, gh_247, gh_249, \
                         gh_275, gh_282, gh_289, gh_291, gi_114, gi_121, gi_128, gi_130, \
                         gi_310, gi_317, gi_324, gi_326, gi_366, gi_373, gi_380, \
                         gi_382 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -13.125 * ab_x[k] * gh_86[k]
                  + 26.25 * ab_x[k] * gh_93[k]
                  + 13.125 * ab_x[k] * gh_100[k]
                  - 26.25 * ab_x[k] * gh_102[k]
                  - 13.125 * ab_x[k] * gh_233[k]
                  + 26.25 * ab_x[k] * gh_240[k]
                  + 13.125 * ab_x[k] * gh_247[k]
                  - 26.25 * ab_x[k] * gh_249[k]
                  + 26.25 * ab_x[k] * gh_275[k]
                  - 52.5 * ab_x[k] * gh_282[k]
                  - 26.25 * ab_x[k] * gh_289[k]
                  + 52.5 * ab_x[k] * gh_291[k]
                  + 13.125 * gi_114[k]
                  - 26.25 * gi_121[k]
                  - 13.125 * gi_128[k]
                  + 26.25 * gi_130[k]
                  + 13.125 * gi_310[k]
                  - 26.25 * gi_317[k]
                  - 13.125 * gi_324[k]
                  + 26.25 * gi_326[k]
                  - 26.25 * gi_366[k]
                  + 52.5 * gi_373[k]
                  + 26.25 * gi_380[k]
                  - 52.5 * gi_382[k];
        g_80[k] = g_40[k];
    }

#pragma omp simd aligned(ab_x, gh_84, gh_87, gh_89, gh_94, gh_96, gh_231, gh_234, gh_236, \
                         gh_241, gh_243, gh_273, gh_276, gh_278, gh_283, gh_285, gi_112, \
                         gi_115, gi_117, gi_122, gi_124, gi_308, gi_311, gi_313, gi_318, \
                         gi_320, gi_364, gi_367, gi_369, gi_374, \
                         gi_376 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_73 * ab_x[k] * gh_84[k]
                  + f_69 * ab_x[k] * gh_87[k]
                  + f_74 * ab_x[k] * gh_89[k]
                  + f_67 * ab_x[k] * gh_94[k]
                  - f_71 * ab_x[k] * gh_96[k]
                  - f_73 * ab_x[k] * gh_231[k]
                  + f_69 * ab_x[k] * gh_234[k]
                  + f_74 * ab_x[k] * gh_236[k]
                  + f_67 * ab_x[k] * gh_241[k]
                  - f_71 * ab_x[k] * gh_243[k]
                  + f_69 * ab_x[k] * gh_273[k]
                  - f_70 * ab_x[k] * gh_276[k]
                  - f_75 * ab_x[k] * gh_278[k]
                  - f_68 * ab_x[k] * gh_283[k]
                  + f_72 * ab_x[k] * gh_285[k]
                  + f_73 * gi_112[k]
                  - f_69 * gi_115[k]
                  - f_74 * gi_117[k]
                  - f_67 * gi_122[k]
                  + f_71 * gi_124[k]
                  + f_73 * gi_308[k]
                  - f_69 * gi_311[k]
                  - f_74 * gi_313[k]
                  - f_67 * gi_318[k]
                  + f_71 * gi_320[k]
                  - f_69 * gi_364[k]
                  + f_70 * gi_367[k]
                  + f_75 * gi_369[k]
                  + f_68 * gi_374[k]
                  - f_72 * gi_376[k];
        g_91[k] = g_41[k];
    }

#pragma omp simd aligned(ab_x, gh_86, gh_91, gh_100, gh_233, gh_238, gh_247, gh_275, gh_280, \
                         gh_289, gi_114, gi_119, gi_128, gi_310, gi_315, gi_324, gi_366, \
                         gi_371, gi_380 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_134 * ab_x[k] * gh_86[k]
                  - f_135 * ab_x[k] * gh_91[k]
                  + f_134 * ab_x[k] * gh_100[k]
                  + f_134 * ab_x[k] * gh_233[k]
                  - f_135 * ab_x[k] * gh_238[k]
                  + f_134 * ab_x[k] * gh_247[k]
                  - f_66 * ab_x[k] * gh_275[k]
                  + f_136 * ab_x[k] * gh_280[k]
                  - f_66 * ab_x[k] * gh_289[k]
                  - f_134 * gi_114[k]
                  + f_135 * gi_119[k]
                  - f_134 * gi_128[k]
                  - f_134 * gi_310[k]
                  + f_135 * gi_315[k]
                  - f_134 * gi_324[k]
                  + f_66 * gi_366[k]
                  - f_136 * gi_371[k]
                  + f_66 * gi_380[k];
        g_102[k] = g_42[k];
    }

#pragma omp simd aligned(ab_x, gh_84, gh_87, gh_94, gh_231, gh_234, gh_241, gh_273, gh_276, \
                         gh_283, gi_112, gi_115, gi_122, gi_308, gi_311, gi_318, gi_364, \
                         gi_367, gi_374 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_20 * ab_x[k] * gh_84[k]
                  - f_18 * ab_x[k] * gh_87[k]
                  + f_17 * ab_x[k] * gh_94[k]
                  + f_20 * ab_x[k] * gh_231[k]
                  - f_18 * ab_x[k] * gh_234[k]
                  + f_17 * ab_x[k] * gh_241[k]
                  - f_21 * ab_x[k] * gh_273[k]
                  + f_19 * ab_x[k] * gh_276[k]
                  - f_18 * ab_x[k] * gh_283[k]
                  - f_20 * gi_112[k]
                  + f_18 * gi_115[k]
                  - f_17 * gi_122[k]
                  - f_20 * gi_308[k]
                  + f_18 * gi_311[k]
                  - f_17 * gi_318[k]
                  + f_21 * gi_364[k]
                  - f_19 * gi_367[k]
                  + f_18 * gi_374[k];
        g_113[k] = g_43[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_22, gh_27, gh_29, gh_36, gh_38, gh_40, gh_127, gh_132, \
                         gh_134, gh_141, gh_143, gh_145, gh_169, gh_174, gh_176, gh_183, \
                         gh_185, gh_187, gh_211, gh_216, gh_218, gh_225, gh_227, gh_229, \
                         gh_253, gh_258, gh_260, gh_267, gh_269, gh_271, gh_295, gh_300, \
                         gh_302, gh_309, gh_311, gh_313, gi_29, gi_34, gi_36, gi_43, gi_45, \
                         gi_47, gi_169, gi_174, gi_176, gi_183, gi_185, gi_187, gi_225, \
                         gi_230, gi_232, gi_239, gi_241, gi_243, gi_283, gi_290, gi_292, \
                         gi_301, gi_303, gi_305, gi_339, gi_346, gi_348, gi_357, gi_359, \
                         gi_361, gi_395, gi_402, gi_404, gi_413, gi_415, \
                         gi_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -0.234375 * ab_x[k] * gh_22[k]
                  - 0.46875 * ab_x[k] * gh_27[k]
                  + 2.8125 * ab_x[k] * gh_29[k]
                  - 0.234375 * ab_x[k] * gh_36[k]
                  + 2.8125 * ab_x[k] * gh_38[k]
                  - 1.875 * ab_x[k] * gh_40[k]
                  - 0.46875 * ab_x[k] * gh_127[k]
                  - 0.9375 * ab_x[k] * gh_132[k]
                  + 5.625 * ab_x[k] * gh_134[k]
                  - 0.46875 * ab_x[k] * gh_141[k]
                  + 5.625 * ab_x[k] * gh_143[k]
                  - 3.75 * ab_x[k] * gh_145[k]
                  + 2.8125 * ab_x[k] * gh_169[k]
                  + 5.625 * ab_x[k] * gh_174[k]
                  - 33.75 * ab_x[k] * gh_176[k]
                  + 2.8125 * ab_x[k] * gh_183[k]
                  - 33.75 * ab_x[k] * gh_185[k]
                  + 22.5 * ab_x[k] * gh_187[k]
                  - 0.234375 * ab_y[k] * gh_211[k]
                  - 0.46875 * ab_y[k] * gh_216[k]
                  + 2.8125 * ab_y[k] * gh_218[k]
                  - 0.234375 * ab_y[k] * gh_225[k]
                  + 2.8125 * ab_y[k] * gh_227[k]
                  - 1.875 * ab_y[k] * gh_229[k]
                  + 2.8125 * ab_y[k] * gh_253[k]
                  + 5.625 * ab_y[k] * gh_258[k]
                  - 33.75 * ab_y[k] * gh_260[k]
                  + 2.8125 * ab_y[k] * gh_267[k]
                  - 33.75 * ab_y[k] * gh_269[k]
                  + 22.5 * ab_y[k] * gh_271[k]
                  - 1.875 * ab_y[k] * gh_295[k]
                  - 3.75 * ab_y[k] * gh_300[k]
                  + 22.5 * ab_y[k] * gh_302[k]
                  - 1.875 * ab_y[k] * gh_309[k]
                  + 22.5 * ab_y[k] * gh_311[k]
                  - 15.0 * ab_y[k] * gh_313[k]
                  + 0.234375 * gi_29[k]
                  + 0.46875 * gi_34[k]
                  - 2.8125 * gi_36[k]
                  + 0.234375 * gi_43[k]
                  - 2.8125 * gi_45[k]
                  + 1.875 * gi_47[k]
                  + 0.46875 * gi_169[k]
                  + 0.9375 * gi_174[k]
                  - 5.625 * gi_176[k]
                  + 0.46875 * gi_183[k]
                  - 5.625 * gi_185[k]
                  + 3.75 * gi_187[k]
                  - 2.8125 * gi_225[k]
                  - 5.625 * gi_230[k]
                  + 33.75 * gi_232[k]
                  - 2.8125 * gi_239[k]
                  + 33.75 * gi_241[k]
                  - 22.5 * gi_243[k]
                  + 0.234375 * gi_283[k]
                  + 0.46875 * gi_290[k]
                  - 2.8125 * gi_292[k]
                  + 0.234375 * gi_301[k]
                  - 2.8125 * gi_303[k]
                  + 1.875 * gi_305[k]
                  - 2.8125 * gi_339[k]
                  - 5.625 * gi_346[k]
                  + 33.75 * gi_348[k]
                  - 2.8125 * gi_357[k]
                  + 33.75 * gi_359[k]
                  - 22.5 * gi_361[k]
                  + 1.875 * gi_395[k]
                  + 3.75 * gi_402[k]
                  - 22.5 * gi_404[k]
                  + 1.875 * gi_413[k]
                  - 22.5 * gi_415[k]
                  + 15.0 * gi_417[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_28, gh_30, gh_37, gh_39, gh_41, gh_128, gh_133, \
                         gh_135, gh_142, gh_144, gh_146, gh_170, gh_175, gh_177, gh_184, \
                         gh_186, gh_188, gh_212, gh_217, gh_219, gh_226, gh_228, gh_230, \
                         gh_254, gh_259, gh_261, gh_268, gh_270, gh_272, gh_296, gh_301, \
                         gh_303, gh_310, gh_312, gh_314, gi_30, gi_35, gi_37, gi_44, gi_46, \
                         gi_48, gi_170, gi_175, gi_177, gi_184, gi_186, gi_188, gi_226, \
                         gi_231, gi_233, gi_240, gi_242, gi_244, gi_284, gi_291, gi_293, \
                         gi_302, gi_304, gi_306, gi_340, gi_347, gi_349, gi_358, gi_360, \
                         gi_362, gi_396, gi_403, gi_405, gi_414, gi_416, \
                         gi_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_137 * ab_x[k] * gh_23[k]
                  - f_138 * ab_x[k] * gh_28[k]
                  + f_139 * ab_x[k] * gh_30[k]
                  - f_137 * ab_x[k] * gh_37[k]
                  + f_139 * ab_x[k] * gh_39[k]
                  - f_140 * ab_x[k] * gh_41[k]
                  - f_138 * ab_x[k] * gh_128[k]
                  - f_141 * ab_x[k] * gh_133[k]
                  + f_142 * ab_x[k] * gh_135[k]
                  - f_138 * ab_x[k] * gh_142[k]
                  + f_142 * ab_x[k] * gh_144[k]
                  - f_143 * ab_x[k] * gh_146[k]
                  + f_144 * ab_x[k] * gh_170[k]
                  + f_145 * ab_x[k] * gh_175[k]
                  - f_146 * ab_x[k] * gh_177[k]
                  + f_144 * ab_x[k] * gh_184[k]
                  - f_146 * ab_x[k] * gh_186[k]
                  + f_147 * ab_x[k] * gh_188[k]
                  - f_137 * ab_y[k] * gh_212[k]
                  - f_138 * ab_y[k] * gh_217[k]
                  + f_139 * ab_y[k] * gh_219[k]
                  - f_137 * ab_y[k] * gh_226[k]
                  + f_139 * ab_y[k] * gh_228[k]
                  - f_140 * ab_y[k] * gh_230[k]
                  + f_144 * ab_y[k] * gh_254[k]
                  + f_145 * ab_y[k] * gh_259[k]
                  - f_146 * ab_y[k] * gh_261[k]
                  + f_144 * ab_y[k] * gh_268[k]
                  - f_146 * ab_y[k] * gh_270[k]
                  + f_147 * ab_y[k] * gh_272[k]
                  - f_148 * ab_y[k] * gh_296[k]
                  - f_149 * ab_y[k] * gh_301[k]
                  + f_150 * ab_y[k] * gh_303[k]
                  - f_148 * ab_y[k] * gh_310[k]
                  + f_150 * ab_y[k] * gh_312[k]
                  - f_151 * ab_y[k] * gh_314[k]
                  + f_137 * gi_30[k]
                  + f_138 * gi_35[k]
                  - f_139 * gi_37[k]
                  + f_137 * gi_44[k]
                  - f_139 * gi_46[k]
                  + f_140 * gi_48[k]
                  + f_138 * gi_170[k]
                  + f_141 * gi_175[k]
                  - f_142 * gi_177[k]
                  + f_138 * gi_184[k]
                  - f_142 * gi_186[k]
                  + f_143 * gi_188[k]
                  - f_144 * gi_226[k]
                  - f_145 * gi_231[k]
                  + f_146 * gi_233[k]
                  - f_144 * gi_240[k]
                  + f_146 * gi_242[k]
                  - f_147 * gi_244[k]
                  + f_137 * gi_284[k]
                  + f_138 * gi_291[k]
                  - f_139 * gi_293[k]
                  + f_137 * gi_302[k]
                  - f_139 * gi_304[k]
                  + f_140 * gi_306[k]
                  - f_144 * gi_340[k]
                  - f_145 * gi_347[k]
                  + f_146 * gi_349[k]
                  - f_144 * gi_358[k]
                  + f_146 * gi_360[k]
                  - f_147 * gi_362[k]
                  + f_148 * gi_396[k]
                  + f_149 * gi_403[k]
                  - f_150 * gi_405[k]
                  + f_148 * gi_414[k]
                  - f_150 * gi_416[k]
                  + f_151 * gi_418[k];
        g_59[k] = g_49[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_26, gh_31, gh_33, gh_35, gh_126, gh_129, \
                         gh_131, gh_136, gh_138, gh_140, gh_168, gh_171, gh_173, gh_178, \
                         gh_180, gh_182, gh_210, gh_213, gh_215, gh_220, gh_222, gh_224, \
                         gh_252, gh_255, gh_257, gh_262, gh_264, gh_266, gh_294, gh_297, \
                         gh_299, gh_304, gh_306, gh_308, gi_28, gi_31, gi_33, gi_38, gi_40, \
                         gi_42, gi_168, gi_171, gi_173, gi_178, gi_180, gi_182, gi_224, \
                         gi_227, gi_229, gi_234, gi_236, gi_238, gi_281, gi_286, gi_288, \
                         gi_295, gi_297, gi_299, gi_337, gi_342, gi_344, gi_351, gi_353, \
                         gi_355, gi_393, gi_398, gi_400, gi_407, gi_409, \
                         gi_411 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -0.234375 * ab_x[k] * gh_21[k]
                  - 0.46875 * ab_x[k] * gh_24[k]
                  + 2.8125 * ab_x[k] * gh_26[k]
                  - 0.234375 * ab_x[k] * gh_31[k]
                  + 2.8125 * ab_x[k] * gh_33[k]
                  - 1.875 * ab_x[k] * gh_35[k]
                  - 0.46875 * ab_x[k] * gh_126[k]
                  - 0.9375 * ab_x[k] * gh_129[k]
                  + 5.625 * ab_x[k] * gh_131[k]
                  - 0.46875 * ab_x[k] * gh_136[k]
                  + 5.625 * ab_x[k] * gh_138[k]
                  - 3.75 * ab_x[k] * gh_140[k]
                  + 2.8125 * ab_x[k] * gh_168[k]
                  + 5.625 * ab_x[k] * gh_171[k]
                  - 33.75 * ab_x[k] * gh_173[k]
                  + 2.8125 * ab_x[k] * gh_178[k]
                  - 33.75 * ab_x[k] * gh_180[k]
                  + 22.5 * ab_x[k] * gh_182[k]
                  - 0.234375 * ab_y[k] * gh_210[k]
                  - 0.46875 * ab_y[k] * gh_213[k]
                  + 2.8125 * ab_y[k] * gh_215[k]
                  - 0.234375 * ab_y[k] * gh_220[k]
                  + 2.8125 * ab_y[k] * gh_222[k]
                  - 1.875 * ab_y[k] * gh_224[k]
                  + 2.8125 * ab_y[k] * gh_252[k]
                  + 5.625 * ab_y[k] * gh_255[k]
                  - 33.75 * ab_y[k] * gh_257[k]
                  + 2.8125 * ab_y[k] * gh_262[k]
                  - 33.75 * ab_y[k] * gh_264[k]
                  + 22.5 * ab_y[k] * gh_266[k]
                  - 1.875 * ab_y[k] * gh_294[k]
                  - 3.75 * ab_y[k] * gh_297[k]
                  + 22.5 * ab_y[k] * gh_299[k]
                  - 1.875 * ab_y[k] * gh_304[k]
                  + 22.5 * ab_y[k] * gh_306[k]
                  - 15.0 * ab_y[k] * gh_308[k]
                  + 0.234375 * gi_28[k]
                  + 0.46875 * gi_31[k]
                  - 2.8125 * gi_33[k]
                  + 0.234375 * gi_38[k]
                  - 2.8125 * gi_40[k]
                  + 1.875 * gi_42[k]
                  + 0.46875 * gi_168[k]
                  + 0.9375 * gi_171[k]
                  - 5.625 * gi_173[k]
                  + 0.46875 * gi_178[k]
                  - 5.625 * gi_180[k]
                  + 3.75 * gi_182[k]
                  - 2.8125 * gi_224[k]
                  - 5.625 * gi_227[k]
                  + 33.75 * gi_229[k]
                  - 2.8125 * gi_234[k]
                  + 33.75 * gi_236[k]
                  - 22.5 * gi_238[k]
                  + 0.234375 * gi_281[k]
                  + 0.46875 * gi_286[k]
                  - 2.8125 * gi_288[k]
                  + 0.234375 * gi_295[k]
                  - 2.8125 * gi_297[k]
                  + 1.875 * gi_299[k]
                  - 2.8125 * gi_337[k]
                  - 5.625 * gi_342[k]
                  + 33.75 * gi_344[k]
                  - 2.8125 * gi_351[k]
                  + 33.75 * gi_353[k]
                  - 22.5 * gi_355[k]
                  + 1.875 * gi_393[k]
                  + 3.75 * gi_398[k]
                  - 22.5 * gi_400[k]
                  + 1.875 * gi_407[k]
                  - 22.5 * gi_409[k]
                  + 15.0 * gi_411[k];
        g_70[k] = g_50[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_30, gh_37, gh_39, gh_128, gh_135, gh_142, \
                         gh_144, gh_170, gh_177, gh_184, gh_186, gh_212, gh_219, gh_226, \
                         gh_228, gh_254, gh_261, gh_268, gh_270, gh_296, gh_303, gh_310, \
                         gh_312, gi_30, gi_37, gi_44, gi_46, gi_170, gi_177, gi_184, gi_186, \
                         gi_226, gi_233, gi_240, gi_242, gi_284, gi_293, gi_302, gi_304, \
                         gi_340, gi_349, gi_358, gi_360, gi_396, gi_405, gi_414, \
                         gi_416 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_152 * ab_x[k] * gh_23[k]
                  - f_120 * ab_x[k] * gh_30[k]
                  - f_152 * ab_x[k] * gh_37[k]
                  + f_120 * ab_x[k] * gh_39[k]
                  + f_120 * ab_x[k] * gh_128[k]
                  - f_121 * ab_x[k] * gh_135[k]
                  - f_120 * ab_x[k] * gh_142[k]
                  + f_121 * ab_x[k] * gh_144[k]
                  - f_153 * ab_x[k] * gh_170[k]
                  + f_122 * ab_x[k] * gh_177[k]
                  + f_153 * ab_x[k] * gh_184[k]
                  - f_122 * ab_x[k] * gh_186[k]
                  + f_152 * ab_y[k] * gh_212[k]
                  - f_120 * ab_y[k] * gh_219[k]
                  - f_152 * ab_y[k] * gh_226[k]
                  + f_120 * ab_y[k] * gh_228[k]
                  - f_153 * ab_y[k] * gh_254[k]
                  + f_122 * ab_y[k] * gh_261[k]
                  + f_153 * ab_y[k] * gh_268[k]
                  - f_122 * ab_y[k] * gh_270[k]
                  + f_124 * ab_y[k] * gh_296[k]
                  - f_123 * ab_y[k] * gh_303[k]
                  - f_124 * ab_y[k] * gh_310[k]
                  + f_123 * ab_y[k] * gh_312[k]
                  - f_152 * gi_30[k]
                  + f_120 * gi_37[k]
                  + f_152 * gi_44[k]
                  - f_120 * gi_46[k]
                  - f_120 * gi_170[k]
                  + f_121 * gi_177[k]
                  + f_120 * gi_184[k]
                  - f_121 * gi_186[k]
                  + f_153 * gi_226[k]
                  - f_122 * gi_233[k]
                  - f_153 * gi_240[k]
                  + f_122 * gi_242[k]
                  - f_152 * gi_284[k]
                  + f_120 * gi_293[k]
                  + f_152 * gi_302[k]
                  - f_120 * gi_304[k]
                  + f_153 * gi_340[k]
                  - f_122 * gi_349[k]
                  - f_153 * gi_358[k]
                  + f_122 * gi_360[k]
                  - f_124 * gi_396[k]
                  + f_123 * gi_405[k]
                  + f_124 * gi_414[k]
                  - f_123 * gi_416[k];
        g_81[k] = g_51[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_26, gh_31, gh_33, gh_126, gh_129, \
                         gh_131, gh_136, gh_138, gh_168, gh_171, gh_173, gh_178, gh_180, \
                         gh_210, gh_213, gh_215, gh_220, gh_222, gh_252, gh_255, gh_257, \
                         gh_262, gh_264, gh_294, gh_297, gh_299, gh_304, gh_306, gi_28, gi_31, \
                         gi_33, gi_38, gi_40, gi_168, gi_171, gi_173, gi_178, gi_180, gi_224, \
                         gi_227, gi_229, gi_234, gi_236, gi_281, gi_286, gi_288, gi_295, \
                         gi_297, gi_337, gi_342, gi_344, gi_351, gi_353, gi_393, gi_398, \
                         gi_400, gi_407, gi_409 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_86 * ab_x[k] * gh_21[k]
                  - f_80 * ab_x[k] * gh_24[k]
                  - f_88 * ab_x[k] * gh_26[k]
                  - f_76 * ab_x[k] * gh_31[k]
                  + f_79 * ab_x[k] * gh_33[k]
                  + f_80 * ab_x[k] * gh_126[k]
                  - f_81 * ab_x[k] * gh_129[k]
                  - f_82 * ab_x[k] * gh_131[k]
                  - f_77 * ab_x[k] * gh_136[k]
                  + f_83 * ab_x[k] * gh_138[k]
                  - f_87 * ab_x[k] * gh_168[k]
                  + f_79 * ab_x[k] * gh_171[k]
                  + f_89 * ab_x[k] * gh_173[k]
                  + f_78 * ab_x[k] * gh_178[k]
                  - f_84 * ab_x[k] * gh_180[k]
                  + f_86 * ab_y[k] * gh_210[k]
                  - f_80 * ab_y[k] * gh_213[k]
                  - f_88 * ab_y[k] * gh_215[k]
                  - f_76 * ab_y[k] * gh_220[k]
                  + f_79 * ab_y[k] * gh_222[k]
                  - f_87 * ab_y[k] * gh_252[k]
                  + f_79 * ab_y[k] * gh_255[k]
                  + f_89 * ab_y[k] * gh_257[k]
                  + f_78 * ab_y[k] * gh_262[k]
                  - f_84 * ab_y[k] * gh_264[k]
                  + f_88 * ab_y[k] * gh_294[k]
                  - f_82 * ab_y[k] * gh_297[k]
                  - f_90 * ab_y[k] * gh_299[k]
                  - f_79 * ab_y[k] * gh_304[k]
                  + f_85 * ab_y[k] * gh_306[k]
                  - f_86 * gi_28[k]
                  + f_80 * gi_31[k]
                  + f_88 * gi_33[k]
                  + f_76 * gi_38[k]
                  - f_79 * gi_40[k]
                  - f_80 * gi_168[k]
                  + f_81 * gi_171[k]
                  + f_82 * gi_173[k]
                  + f_77 * gi_178[k]
                  - f_83 * gi_180[k]
                  + f_87 * gi_224[k]
                  - f_79 * gi_227[k]
                  - f_89 * gi_229[k]
                  - f_78 * gi_234[k]
                  + f_84 * gi_236[k]
                  - f_86 * gi_281[k]
                  + f_80 * gi_286[k]
                  + f_88 * gi_288[k]
                  + f_76 * gi_295[k]
                  - f_79 * gi_297[k]
                  + f_87 * gi_337[k]
                  - f_79 * gi_342[k]
                  - f_89 * gi_344[k]
                  - f_78 * gi_351[k]
                  + f_84 * gi_353[k]
                  - f_88 * gi_393[k]
                  + f_82 * gi_398[k]
                  + f_90 * gi_400[k]
                  + f_79 * gi_407[k]
                  - f_85 * gi_409[k];
        g_92[k] = g_52[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_23, gh_28, gh_37, gh_128, gh_133, gh_142, gh_170, \
                         gh_175, gh_184, gh_212, gh_217, gh_226, gh_254, gh_259, gh_268, \
                         gh_296, gh_301, gh_310, gi_30, gi_35, gi_44, gi_170, gi_175, gi_184, \
                         gi_226, gi_231, gi_240, gi_284, gi_291, gi_302, gi_340, gi_347, \
                         gi_358, gi_396, gi_403, gi_414 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = -f_154 * ab_x[k] * gh_23[k]
                  + f_155 * ab_x[k] * gh_28[k]
                  - f_154 * ab_x[k] * gh_37[k]
                  - f_156 * ab_x[k] * gh_128[k]
                  + f_157 * ab_x[k] * gh_133[k]
                  - f_156 * ab_x[k] * gh_142[k]
                  + f_157 * ab_x[k] * gh_170[k]
                  - f_158 * ab_x[k] * gh_175[k]
                  + f_157 * ab_x[k] * gh_184[k]
                  - f_154 * ab_y[k] * gh_212[k]
                  + f_155 * ab_y[k] * gh_217[k]
                  - f_154 * ab_y[k] * gh_226[k]
                  + f_157 * ab_y[k] * gh_254[k]
                  - f_158 * ab_y[k] * gh_259[k]
                  + f_157 * ab_y[k] * gh_268[k]
                  - f_59 * ab_y[k] * gh_296[k]
                  + f_60 * ab_y[k] * gh_301[k]
                  - f_59 * ab_y[k] * gh_310[k]
                  + f_154 * gi_30[k]
                  - f_155 * gi_35[k]
                  + f_154 * gi_44[k]
                  + f_156 * gi_170[k]
                  - f_157 * gi_175[k]
                  + f_156 * gi_184[k]
                  - f_157 * gi_226[k]
                  + f_158 * gi_231[k]
                  - f_157 * gi_240[k]
                  + f_154 * gi_284[k]
                  - f_155 * gi_291[k]
                  + f_154 * gi_302[k]
                  - f_157 * gi_340[k]
                  + f_158 * gi_347[k]
                  - f_157 * gi_358[k]
                  + f_59 * gi_396[k]
                  - f_60 * gi_403[k]
                  + f_59 * gi_414[k];
        g_103[k] = g_53[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_21, gh_24, gh_31, gh_126, gh_129, gh_136, gh_168, \
                         gh_171, gh_178, gh_210, gh_213, gh_220, gh_252, gh_255, gh_262, \
                         gh_294, gh_297, gh_304, gi_28, gi_31, gi_38, gi_168, gi_171, gi_178, \
                         gi_224, gi_227, gi_234, gi_281, gi_286, gi_295, gi_337, gi_342, \
                         gi_351, gi_393, gi_398, gi_407 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_29 * ab_x[k] * gh_21[k]
                  + f_23 * ab_x[k] * gh_24[k]
                  - f_22 * ab_x[k] * gh_31[k]
                  - f_30 * ab_x[k] * gh_126[k]
                  + f_26 * ab_x[k] * gh_129[k]
                  - f_23 * ab_x[k] * gh_136[k]
                  + f_31 * ab_x[k] * gh_168[k]
                  - f_27 * ab_x[k] * gh_171[k]
                  + f_24 * ab_x[k] * gh_178[k]
                  - f_29 * ab_y[k] * gh_210[k]
                  + f_23 * ab_y[k] * gh_213[k]
                  - f_22 * ab_y[k] * gh_220[k]
                  + f_31 * ab_y[k] * gh_252[k]
                  - f_27 * ab_y[k] * gh_255[k]
                  + f_24 * ab_y[k] * gh_262[k]
                  - f_32 * ab_y[k] * gh_294[k]
                  + f_28 * ab_y[k] * gh_297[k]
                  - f_25 * ab_y[k] * gh_304[k]
                  + f_29 * gi_28[k]
                  - f_23 * gi_31[k]
                  + f_22 * gi_38[k]
                  + f_30 * gi_168[k]
                  - f_26 * gi_171[k]
                  + f_23 * gi_178[k]
                  - f_31 * gi_224[k]
                  + f_27 * gi_227[k]
                  - f_24 * gi_234[k]
                  + f_29 * gi_281[k]
                  - f_23 * gi_286[k]
                  + f_22 * gi_295[k]
                  - f_31 * gi_337[k]
                  + f_27 * gi_342[k]
                  - f_24 * gi_351[k]
                  + f_32 * gi_393[k]
                  - f_28 * gi_398[k]
                  + f_25 * gi_407[k];
        g_114[k] = g_54[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_44, gh_49, gh_51, gh_58, gh_60, gh_62, gh_149, \
                         gh_154, gh_156, gh_163, gh_165, gh_167, gh_191, gh_196, gh_198, \
                         gh_205, gh_207, gh_209, gh_233, gh_238, gh_240, gh_247, gh_249, \
                         gh_251, gh_275, gh_280, gh_282, gh_289, gh_291, gh_293, gh_296, \
                         gh_301, gh_303, gh_310, gh_312, gh_314, gi_58, gi_63, gi_65, gi_72, \
                         gi_74, gi_76, gi_198, gi_203, gi_205, gi_212, gi_214, gi_216, gi_254, \
                         gi_259, gi_261, gi_268, gi_270, gi_272, gi_312, gi_319, gi_321, \
                         gi_330, gi_332, gi_334, gi_368, gi_375, gi_377, gi_386, gi_388, \
                         gi_390, gi_397, gi_404, gi_406, gi_415, gi_417, \
                         gi_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -3.515625 * ab_x[k] * gh_44[k]
                  - 7.03125 * ab_x[k] * gh_49[k]
                  + 9.375 * ab_x[k] * gh_51[k]
                  - 3.515625 * ab_x[k] * gh_58[k]
                  + 9.375 * ab_x[k] * gh_60[k]
                  - 1.875 * ab_x[k] * gh_62[k]
                  - 7.03125 * ab_x[k] * gh_149[k]
                  - 14.0625 * ab_x[k] * gh_154[k]
                  + 18.75 * ab_x[k] * gh_156[k]
                  - 7.03125 * ab_x[k] * gh_163[k]
                  + 18.75 * ab_x[k] * gh_165[k]
                  - 3.75 * ab_x[k] * gh_167[k]
                  + 9.375 * ab_x[k] * gh_191[k]
                  + 18.75 * ab_x[k] * gh_196[k]
                  - 25.0 * ab_x[k] * gh_198[k]
                  + 9.375 * ab_x[k] * gh_205[k]
                  - 25.0 * ab_x[k] * gh_207[k]
                  + 5.0 * ab_x[k] * gh_209[k]
                  - 3.515625 * ab_y[k] * gh_233[k]
                  - 7.03125 * ab_y[k] * gh_238[k]
                  + 9.375 * ab_y[k] * gh_240[k]
                  - 3.515625 * ab_y[k] * gh_247[k]
                  + 9.375 * ab_y[k] * gh_249[k]
                  - 1.875 * ab_y[k] * gh_251[k]
                  + 9.375 * ab_y[k] * gh_275[k]
                  + 18.75 * ab_y[k] * gh_280[k]
                  - 25.0 * ab_y[k] * gh_282[k]
                  + 9.375 * ab_y[k] * gh_289[k]
                  - 25.0 * ab_y[k] * gh_291[k]
                  + 5.0 * ab_y[k] * gh_293[k]
                  - 1.875 * ab_z[k] * gh_296[k]
                  - 3.75 * ab_z[k] * gh_301[k]
                  + 5.0 * ab_z[k] * gh_303[k]
                  - 1.875 * ab_z[k] * gh_310[k]
                  + 5.0 * ab_z[k] * gh_312[k]
                  - ab_z[k] * gh_314[k]
                  + 3.515625 * gi_58[k]
                  + 7.03125 * gi_63[k]
                  - 9.375 * gi_65[k]
                  + 3.515625 * gi_72[k]
                  - 9.375 * gi_74[k]
                  + 1.875 * gi_76[k]
                  + 7.03125 * gi_198[k]
                  + 14.0625 * gi_203[k]
                  - 18.75 * gi_205[k]
                  + 7.03125 * gi_212[k]
                  - 18.75 * gi_214[k]
                  + 3.75 * gi_216[k]
                  - 9.375 * gi_254[k]
                  - 18.75 * gi_259[k]
                  + 25.0 * gi_261[k]
                  - 9.375 * gi_268[k]
                  + 25.0 * gi_270[k]
                  - 5.0 * gi_272[k]
                  + 3.515625 * gi_312[k]
                  + 7.03125 * gi_319[k]
                  - 9.375 * gi_321[k]
                  + 3.515625 * gi_330[k]
                  - 9.375 * gi_332[k]
                  + 1.875 * gi_334[k]
                  - 9.375 * gi_368[k]
                  - 18.75 * gi_375[k]
                  + 25.0 * gi_377[k]
                  - 9.375 * gi_386[k]
                  + 25.0 * gi_388[k]
                  - 5.0 * gi_390[k]
                  + 1.875 * gi_397[k]
                  + 3.75 * gi_404[k]
                  - 5.0 * gi_406[k]
                  + 1.875 * gi_415[k]
                  - 5.0 * gi_417[k]
                  + gi_419[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_42, gh_45, gh_47, gh_52, gh_54, gh_56, gh_147, \
                         gh_150, gh_152, gh_157, gh_159, gh_161, gh_189, gh_192, gh_194, \
                         gh_199, gh_201, gh_203, gh_231, gh_234, gh_236, gh_241, gh_243, \
                         gh_245, gh_273, gh_276, gh_278, gh_283, gh_285, gh_287, gh_294, \
                         gh_297, gh_299, gh_304, gh_306, gh_308, gi_56, gi_59, gi_61, gi_66, \
                         gi_68, gi_70, gi_196, gi_199, gi_201, gi_206, gi_208, gi_210, gi_252, \
                         gi_255, gi_257, gi_262, gi_264, gi_266, gi_309, gi_314, gi_316, \
                         gi_323, gi_325, gi_327, gi_365, gi_370, gi_372, gi_379, gi_381, \
                         gi_383, gi_394, gi_399, gi_401, gi_408, gi_410, \
                         gi_412 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_137 * ab_x[k] * gh_42[k]
                  - f_138 * ab_x[k] * gh_45[k]
                  + f_144 * ab_x[k] * gh_47[k]
                  - f_137 * ab_x[k] * gh_52[k]
                  + f_144 * ab_x[k] * gh_54[k]
                  - f_148 * ab_x[k] * gh_56[k]
                  - f_138 * ab_x[k] * gh_147[k]
                  - f_141 * ab_x[k] * gh_150[k]
                  + f_145 * ab_x[k] * gh_152[k]
                  - f_138 * ab_x[k] * gh_157[k]
                  + f_145 * ab_x[k] * gh_159[k]
                  - f_149 * ab_x[k] * gh_161[k]
                  + f_139 * ab_x[k] * gh_189[k]
                  + f_142 * ab_x[k] * gh_192[k]
                  - f_146 * ab_x[k] * gh_194[k]
                  + f_139 * ab_x[k] * gh_199[k]
                  - f_146 * ab_x[k] * gh_201[k]
                  + f_150 * ab_x[k] * gh_203[k]
                  - f_137 * ab_y[k] * gh_231[k]
                  - f_138 * ab_y[k] * gh_234[k]
                  + f_144 * ab_y[k] * gh_236[k]
                  - f_137 * ab_y[k] * gh_241[k]
                  + f_144 * ab_y[k] * gh_243[k]
                  - f_148 * ab_y[k] * gh_245[k]
                  + f_139 * ab_y[k] * gh_273[k]
                  + f_142 * ab_y[k] * gh_276[k]
                  - f_146 * ab_y[k] * gh_278[k]
                  + f_139 * ab_y[k] * gh_283[k]
                  - f_146 * ab_y[k] * gh_285[k]
                  + f_150 * ab_y[k] * gh_287[k]
                  - f_140 * ab_z[k] * gh_294[k]
                  - f_143 * ab_z[k] * gh_297[k]
                  + f_147 * ab_z[k] * gh_299[k]
                  - f_140 * ab_z[k] * gh_304[k]
                  + f_147 * ab_z[k] * gh_306[k]
                  - f_151 * ab_z[k] * gh_308[k]
                  + f_137 * gi_56[k]
                  + f_138 * gi_59[k]
                  - f_144 * gi_61[k]
                  + f_137 * gi_66[k]
                  - f_144 * gi_68[k]
                  + f_148 * gi_70[k]
                  + f_138 * gi_196[k]
                  + f_141 * gi_199[k]
                  - f_145 * gi_201[k]
                  + f_138 * gi_206[k]
                  - f_145 * gi_208[k]
                  + f_149 * gi_210[k]
                  - f_139 * gi_252[k]
                  - f_142 * gi_255[k]
                  + f_146 * gi_257[k]
                  - f_139 * gi_262[k]
                  + f_146 * gi_264[k]
                  - f_150 * gi_266[k]
                  + f_137 * gi_309[k]
                  + f_138 * gi_314[k]
                  - f_144 * gi_316[k]
                  + f_137 * gi_323[k]
                  - f_144 * gi_325[k]
                  + f_148 * gi_327[k]
                  - f_139 * gi_365[k]
                  - f_142 * gi_370[k]
                  + f_146 * gi_372[k]
                  - f_139 * gi_379[k]
                  + f_146 * gi_381[k]
                  - f_150 * gi_383[k]
                  + f_140 * gi_394[k]
                  + f_143 * gi_399[k]
                  - f_147 * gi_401[k]
                  + f_140 * gi_408[k]
                  - f_147 * gi_410[k]
                  + f_151 * gi_412[k];
        g_71[k] = g_61[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_44, gh_51, gh_58, gh_60, gh_149, gh_156, gh_163, \
                         gh_165, gh_191, gh_198, gh_205, gh_207, gh_233, gh_240, gh_247, \
                         gh_249, gh_275, gh_282, gh_289, gh_291, gh_296, gh_303, gh_310, \
                         gh_312, gi_58, gi_65, gi_72, gi_74, gi_198, gi_205, gi_212, gi_214, \
                         gi_254, gi_261, gi_268, gi_270, gi_312, gi_321, gi_330, gi_332, \
                         gi_368, gi_377, gi_386, gi_388, gi_397, gi_406, gi_415, \
                         gi_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_159 * ab_x[k] * gh_44[k]
                  - f_127 * ab_x[k] * gh_51[k]
                  - f_159 * ab_x[k] * gh_58[k]
                  + f_127 * ab_x[k] * gh_60[k]
                  + f_127 * ab_x[k] * gh_149[k]
                  - f_128 * ab_x[k] * gh_156[k]
                  - f_127 * ab_x[k] * gh_163[k]
                  + f_128 * ab_x[k] * gh_165[k]
                  - f_160 * ab_x[k] * gh_191[k]
                  + f_129 * ab_x[k] * gh_198[k]
                  + f_160 * ab_x[k] * gh_205[k]
                  - f_129 * ab_x[k] * gh_207[k]
                  + f_159 * ab_y[k] * gh_233[k]
                  - f_127 * ab_y[k] * gh_240[k]
                  - f_159 * ab_y[k] * gh_247[k]
                  + f_127 * ab_y[k] * gh_249[k]
                  - f_160 * ab_y[k] * gh_275[k]
                  + f_129 * ab_y[k] * gh_282[k]
                  + f_160 * ab_y[k] * gh_289[k]
                  - f_129 * ab_y[k] * gh_291[k]
                  + f_161 * ab_z[k] * gh_296[k]
                  - f_130 * ab_z[k] * gh_303[k]
                  - f_161 * ab_z[k] * gh_310[k]
                  + f_130 * ab_z[k] * gh_312[k]
                  - f_159 * gi_58[k]
                  + f_127 * gi_65[k]
                  + f_159 * gi_72[k]
                  - f_127 * gi_74[k]
                  - f_127 * gi_198[k]
                  + f_128 * gi_205[k]
                  + f_127 * gi_212[k]
                  - f_128 * gi_214[k]
                  + f_160 * gi_254[k]
                  - f_129 * gi_261[k]
                  - f_160 * gi_268[k]
                  + f_129 * gi_270[k]
                  - f_159 * gi_312[k]
                  + f_127 * gi_321[k]
                  + f_159 * gi_330[k]
                  - f_127 * gi_332[k]
                  + f_160 * gi_368[k]
                  - f_129 * gi_377[k]
                  - f_160 * gi_386[k]
                  + f_129 * gi_388[k]
                  - f_161 * gi_397[k]
                  + f_130 * gi_406[k]
                  + f_161 * gi_415[k]
                  - f_130 * gi_417[k];
        g_82[k] = g_62[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_42, gh_45, gh_47, gh_52, gh_54, gh_147, gh_150, \
                         gh_152, gh_157, gh_159, gh_189, gh_192, gh_194, gh_199, gh_201, \
                         gh_231, gh_234, gh_236, gh_241, gh_243, gh_273, gh_276, gh_278, \
                         gh_283, gh_285, gh_294, gh_297, gh_299, gh_304, gh_306, gi_56, gi_59, \
                         gi_61, gi_66, gi_68, gi_196, gi_199, gi_201, gi_206, gi_208, gi_252, \
                         gi_255, gi_257, gi_262, gi_264, gi_309, gi_314, gi_316, gi_323, \
                         gi_325, gi_365, gi_370, gi_372, gi_379, gi_381, gi_394, gi_399, \
                         gi_401, gi_408, gi_410 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_103 * ab_x[k] * gh_42[k]
                  - f_95 * ab_x[k] * gh_45[k]
                  - f_93 * ab_x[k] * gh_47[k]
                  - f_91 * ab_x[k] * gh_52[k]
                  + f_99 * ab_x[k] * gh_54[k]
                  + f_95 * ab_x[k] * gh_147[k]
                  - f_96 * ab_x[k] * gh_150[k]
                  - f_106 * ab_x[k] * gh_152[k]
                  - f_92 * ab_x[k] * gh_157[k]
                  + f_100 * ab_x[k] * gh_159[k]
                  - f_104 * ab_x[k] * gh_189[k]
                  + f_97 * ab_x[k] * gh_192[k]
                  + f_107 * ab_x[k] * gh_194[k]
                  + f_93 * ab_x[k] * gh_199[k]
                  - f_101 * ab_x[k] * gh_201[k]
                  + f_103 * ab_y[k] * gh_231[k]
                  - f_95 * ab_y[k] * gh_234[k]
                  - f_93 * ab_y[k] * gh_236[k]
                  - f_91 * ab_y[k] * gh_241[k]
                  + f_99 * ab_y[k] * gh_243[k]
                  - f_104 * ab_y[k] * gh_273[k]
                  + f_97 * ab_y[k] * gh_276[k]
                  + f_107 * ab_y[k] * gh_278[k]
                  + f_93 * ab_y[k] * gh_283[k]
                  - f_101 * ab_y[k] * gh_285[k]
                  + f_105 * ab_z[k] * gh_294[k]
                  - f_98 * ab_z[k] * gh_297[k]
                  - f_108 * ab_z[k] * gh_299[k]
                  - f_94 * ab_z[k] * gh_304[k]
                  + f_102 * ab_z[k] * gh_306[k]
                  - f_103 * gi_56[k]
                  + f_95 * gi_59[k]
                  + f_93 * gi_61[k]
                  + f_91 * gi_66[k]
                  - f_99 * gi_68[k]
                  - f_95 * gi_196[k]
                  + f_96 * gi_199[k]
                  + f_106 * gi_201[k]
                  + f_92 * gi_206[k]
                  - f_100 * gi_208[k]
                  + f_104 * gi_252[k]
                  - f_97 * gi_255[k]
                  - f_107 * gi_257[k]
                  - f_93 * gi_262[k]
                  + f_101 * gi_264[k]
                  - f_103 * gi_309[k]
                  + f_95 * gi_314[k]
                  + f_93 * gi_316[k]
                  + f_91 * gi_323[k]
                  - f_99 * gi_325[k]
                  + f_104 * gi_365[k]
                  - f_97 * gi_370[k]
                  - f_107 * gi_372[k]
                  - f_93 * gi_379[k]
                  + f_101 * gi_381[k]
                  - f_105 * gi_394[k]
                  + f_98 * gi_399[k]
                  + f_108 * gi_401[k]
                  + f_94 * gi_408[k]
                  - f_102 * gi_410[k];
        g_93[k] = g_63[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_44, gh_49, gh_58, gh_149, gh_154, gh_163, \
                         gh_191, gh_196, gh_205, gh_233, gh_238, gh_247, gh_275, gh_280, \
                         gh_289, gh_296, gh_301, gh_310, gi_58, gi_63, gi_72, gi_198, gi_203, \
                         gi_212, gi_254, gi_259, gi_268, gi_312, gi_319, gi_330, gi_368, \
                         gi_375, gi_386, gi_397, gi_404, gi_415 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_162 * ab_x[k] * gh_44[k]
                  + f_163 * ab_x[k] * gh_49[k]
                  - f_162 * ab_x[k] * gh_58[k]
                  - f_164 * ab_x[k] * gh_149[k]
                  + f_165 * ab_x[k] * gh_154[k]
                  - f_164 * ab_x[k] * gh_163[k]
                  + f_166 * ab_x[k] * gh_191[k]
                  - f_167 * ab_x[k] * gh_196[k]
                  + f_166 * ab_x[k] * gh_205[k]
                  - f_162 * ab_y[k] * gh_233[k]
                  + f_163 * ab_y[k] * gh_238[k]
                  - f_162 * ab_y[k] * gh_247[k]
                  + f_166 * ab_y[k] * gh_275[k]
                  - f_167 * ab_y[k] * gh_280[k]
                  + f_166 * ab_y[k] * gh_289[k]
                  - f_168 * ab_z[k] * gh_296[k]
                  + f_169 * ab_z[k] * gh_301[k]
                  - f_168 * ab_z[k] * gh_310[k]
                  + f_162 * gi_58[k]
                  - f_163 * gi_63[k]
                  + f_162 * gi_72[k]
                  + f_164 * gi_198[k]
                  - f_165 * gi_203[k]
                  + f_164 * gi_212[k]
                  - f_166 * gi_254[k]
                  + f_167 * gi_259[k]
                  - f_166 * gi_268[k]
                  + f_162 * gi_312[k]
                  - f_163 * gi_319[k]
                  + f_162 * gi_330[k]
                  - f_166 * gi_368[k]
                  + f_167 * gi_375[k]
                  - f_166 * gi_386[k]
                  + f_168 * gi_397[k]
                  - f_169 * gi_404[k]
                  + f_168 * gi_415[k];
        g_104[k] = g_64[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gh_42, gh_45, gh_52, gh_147, gh_150, gh_157, \
                         gh_189, gh_192, gh_199, gh_231, gh_234, gh_241, gh_273, gh_276, \
                         gh_283, gh_294, gh_297, gh_304, gi_56, gi_59, gi_66, gi_196, gi_199, \
                         gi_206, gi_252, gi_255, gi_262, gi_309, gi_314, gi_323, gi_365, \
                         gi_370, gi_379, gi_394, gi_399, gi_408 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = -f_40 * ab_x[k] * gh_42[k]
                  + f_34 * ab_x[k] * gh_45[k]
                  - f_33 * ab_x[k] * gh_52[k]
                  - f_41 * ab_x[k] * gh_147[k]
                  + f_37 * ab_x[k] * gh_150[k]
                  - f_34 * ab_x[k] * gh_157[k]
                  + f_36 * ab_x[k] * gh_189[k]
                  - f_38 * ab_x[k] * gh_192[k]
                  + f_35 * ab_x[k] * gh_199[k]
                  - f_40 * ab_y[k] * gh_231[k]
                  + f_34 * ab_y[k] * gh_234[k]
                  - f_33 * ab_y[k] * gh_241[k]
                  + f_36 * ab_y[k] * gh_273[k]
                  - f_38 * ab_y[k] * gh_276[k]
                  + f_35 * ab_y[k] * gh_283[k]
                  - f_42 * ab_z[k] * gh_294[k]
                  + f_39 * ab_z[k] * gh_297[k]
                  - f_36 * ab_z[k] * gh_304[k]
                  + f_40 * gi_56[k]
                  - f_34 * gi_59[k]
                  + f_33 * gi_66[k]
                  + f_41 * gi_196[k]
                  - f_37 * gi_199[k]
                  + f_34 * gi_206[k]
                  - f_36 * gi_252[k]
                  + f_38 * gi_255[k]
                  - f_35 * gi_262[k]
                  + f_40 * gi_309[k]
                  - f_34 * gi_314[k]
                  + f_33 * gi_323[k]
                  - f_36 * gi_365[k]
                  + f_38 * gi_370[k]
                  - f_35 * gi_379[k]
                  + f_42 * gi_394[k]
                  - f_39 * gi_399[k]
                  + f_36 * gi_408[k];
        g_115[k] = g_65[k];
    }

#pragma omp simd aligned(ab_x, gh_0, gh_3, gh_5, gh_10, gh_12, gh_14, gh_63, gh_66, gh_68, \
                         gh_73, gh_75, gh_77, gh_105, gh_108, gh_110, gh_115, gh_117, gh_119, \
                         gh_210, gh_213, gh_215, gh_220, gh_222, gh_224, gh_252, gh_255, \
                         gh_257, gh_262, gh_264, gh_266, gh_294, gh_297, gh_299, gh_304, \
                         gh_306, gh_308, gi_0, gi_3, gi_5, gi_10, gi_12, gi_14, gi_84, gi_87, \
                         gi_89, gi_94, gi_96, gi_98, gi_140, gi_143, gi_145, gi_150, gi_152, \
                         gi_154, gi_280, gi_283, gi_285, gi_290, gi_292, gi_294, gi_336, \
                         gi_339, gi_341, gi_346, gi_348, gi_350, gi_392, gi_395, gi_397, \
                         gi_402, gi_404, gi_406 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -0.234375 * ab_x[k] * gh_0[k]
                  - 0.46875 * ab_x[k] * gh_3[k]
                  + 2.8125 * ab_x[k] * gh_5[k]
                  - 0.234375 * ab_x[k] * gh_10[k]
                  + 2.8125 * ab_x[k] * gh_12[k]
                  - 1.875 * ab_x[k] * gh_14[k]
                  - 0.46875 * ab_x[k] * gh_63[k]
                  - 0.9375 * ab_x[k] * gh_66[k]
                  + 5.625 * ab_x[k] * gh_68[k]
                  - 0.46875 * ab_x[k] * gh_73[k]
                  + 5.625 * ab_x[k] * gh_75[k]
                  - 3.75 * ab_x[k] * gh_77[k]
                  + 2.8125 * ab_x[k] * gh_105[k]
                  + 5.625 * ab_x[k] * gh_108[k]
                  - 33.75 * ab_x[k] * gh_110[k]
                  + 2.8125 * ab_x[k] * gh_115[k]
                  - 33.75 * ab_x[k] * gh_117[k]
                  + 22.5 * ab_x[k] * gh_119[k]
                  - 0.234375 * ab_x[k] * gh_210[k]
                  - 0.46875 * ab_x[k] * gh_213[k]
                  + 2.8125 * ab_x[k] * gh_215[k]
                  - 0.234375 * ab_x[k] * gh_220[k]
                  + 2.8125 * ab_x[k] * gh_222[k]
                  - 1.875 * ab_x[k] * gh_224[k]
                  + 2.8125 * ab_x[k] * gh_252[k]
                  + 5.625 * ab_x[k] * gh_255[k]
                  - 33.75 * ab_x[k] * gh_257[k]
                  + 2.8125 * ab_x[k] * gh_262[k]
                  - 33.75 * ab_x[k] * gh_264[k]
                  + 22.5 * ab_x[k] * gh_266[k]
                  - 1.875 * ab_x[k] * gh_294[k]
                  - 3.75 * ab_x[k] * gh_297[k]
                  + 22.5 * ab_x[k] * gh_299[k]
                  - 1.875 * ab_x[k] * gh_304[k]
                  + 22.5 * ab_x[k] * gh_306[k]
                  - 15.0 * ab_x[k] * gh_308[k]
                  + 0.234375 * gi_0[k]
                  + 0.46875 * gi_3[k]
                  - 2.8125 * gi_5[k]
                  + 0.234375 * gi_10[k]
                  - 2.8125 * gi_12[k]
                  + 1.875 * gi_14[k]
                  + 0.46875 * gi_84[k]
                  + 0.9375 * gi_87[k]
                  - 5.625 * gi_89[k]
                  + 0.46875 * gi_94[k]
                  - 5.625 * gi_96[k]
                  + 3.75 * gi_98[k]
                  - 2.8125 * gi_140[k]
                  - 5.625 * gi_143[k]
                  + 33.75 * gi_145[k]
                  - 2.8125 * gi_150[k]
                  + 33.75 * gi_152[k]
                  - 22.5 * gi_154[k]
                  + 0.234375 * gi_280[k]
                  + 0.46875 * gi_283[k]
                  - 2.8125 * gi_285[k]
                  + 0.234375 * gi_290[k]
                  - 2.8125 * gi_292[k]
                  + 1.875 * gi_294[k]
                  - 2.8125 * gi_336[k]
                  - 5.625 * gi_339[k]
                  + 33.75 * gi_341[k]
                  - 2.8125 * gi_346[k]
                  + 33.75 * gi_348[k]
                  - 22.5 * gi_350[k]
                  + 1.875 * gi_392[k]
                  + 3.75 * gi_395[k]
                  - 22.5 * gi_397[k]
                  + 1.875 * gi_402[k]
                  - 22.5 * gi_404[k]
                  + 15.0 * gi_406[k];
    }

#pragma omp simd aligned(ab_x, gh_2, gh_9, gh_16, gh_18, gh_65, gh_72, gh_79, gh_81, gh_107, \
                         gh_114, gh_121, gh_123, gh_212, gh_219, gh_226, gh_228, gh_254, \
                         gh_261, gh_268, gh_270, gh_296, gh_303, gh_310, gh_312, gi_2, gi_9, \
                         gi_16, gi_18, gi_86, gi_93, gi_100, gi_102, gi_142, gi_149, gi_156, \
                         gi_158, gi_282, gi_289, gi_296, gi_298, gi_338, gi_345, gi_352, \
                         gi_354, gi_394, gi_401, gi_408, gi_410 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_152 * ab_x[k] * gh_2[k]
                  - f_120 * ab_x[k] * gh_9[k]
                  - f_152 * ab_x[k] * gh_16[k]
                  + f_120 * ab_x[k] * gh_18[k]
                  + f_120 * ab_x[k] * gh_65[k]
                  - f_121 * ab_x[k] * gh_72[k]
                  - f_120 * ab_x[k] * gh_79[k]
                  + f_121 * ab_x[k] * gh_81[k]
                  - f_153 * ab_x[k] * gh_107[k]
                  + f_122 * ab_x[k] * gh_114[k]
                  + f_153 * ab_x[k] * gh_121[k]
                  - f_122 * ab_x[k] * gh_123[k]
                  + f_152 * ab_x[k] * gh_212[k]
                  - f_120 * ab_x[k] * gh_219[k]
                  - f_152 * ab_x[k] * gh_226[k]
                  + f_120 * ab_x[k] * gh_228[k]
                  - f_153 * ab_x[k] * gh_254[k]
                  + f_122 * ab_x[k] * gh_261[k]
                  + f_153 * ab_x[k] * gh_268[k]
                  - f_122 * ab_x[k] * gh_270[k]
                  + f_124 * ab_x[k] * gh_296[k]
                  - f_123 * ab_x[k] * gh_303[k]
                  - f_124 * ab_x[k] * gh_310[k]
                  + f_123 * ab_x[k] * gh_312[k]
                  - f_152 * gi_2[k]
                  + f_120 * gi_9[k]
                  + f_152 * gi_16[k]
                  - f_120 * gi_18[k]
                  - f_120 * gi_86[k]
                  + f_121 * gi_93[k]
                  + f_120 * gi_100[k]
                  - f_121 * gi_102[k]
                  + f_153 * gi_142[k]
                  - f_122 * gi_149[k]
                  - f_153 * gi_156[k]
                  + f_122 * gi_158[k]
                  - f_152 * gi_282[k]
                  + f_120 * gi_289[k]
                  + f_152 * gi_296[k]
                  - f_120 * gi_298[k]
                  + f_153 * gi_338[k]
                  - f_122 * gi_345[k]
                  - f_153 * gi_352[k]
                  + f_122 * gi_354[k]
                  - f_124 * gi_394[k]
                  + f_123 * gi_401[k]
                  + f_124 * gi_408[k]
                  - f_123 * gi_410[k];
        g_83[k] = g_73[k];
    }

#pragma omp simd aligned(ab_x, gh_0, gh_3, gh_5, gh_10, gh_12, gh_63, gh_66, gh_68, gh_73, \
                         gh_75, gh_105, gh_108, gh_110, gh_115, gh_117, gh_210, gh_213, \
                         gh_215, gh_220, gh_222, gh_252, gh_255, gh_257, gh_262, gh_264, \
                         gh_294, gh_297, gh_299, gh_304, gh_306, gi_0, gi_3, gi_5, gi_10, \
                         gi_12, gi_84, gi_87, gi_89, gi_94, gi_96, gi_140, gi_143, gi_145, \
                         gi_150, gi_152, gi_280, gi_283, gi_285, gi_290, gi_292, gi_336, \
                         gi_339, gi_341, gi_346, gi_348, gi_392, gi_395, gi_397, gi_402, \
                         gi_404 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_86 * ab_x[k] * gh_0[k]
                  - f_80 * ab_x[k] * gh_3[k]
                  - f_88 * ab_x[k] * gh_5[k]
                  - f_76 * ab_x[k] * gh_10[k]
                  + f_79 * ab_x[k] * gh_12[k]
                  + f_80 * ab_x[k] * gh_63[k]
                  - f_81 * ab_x[k] * gh_66[k]
                  - f_82 * ab_x[k] * gh_68[k]
                  - f_77 * ab_x[k] * gh_73[k]
                  + f_83 * ab_x[k] * gh_75[k]
                  - f_87 * ab_x[k] * gh_105[k]
                  + f_79 * ab_x[k] * gh_108[k]
                  + f_89 * ab_x[k] * gh_110[k]
                  + f_78 * ab_x[k] * gh_115[k]
                  - f_84 * ab_x[k] * gh_117[k]
                  + f_86 * ab_x[k] * gh_210[k]
                  - f_80 * ab_x[k] * gh_213[k]
                  - f_88 * ab_x[k] * gh_215[k]
                  - f_76 * ab_x[k] * gh_220[k]
                  + f_79 * ab_x[k] * gh_222[k]
                  - f_87 * ab_x[k] * gh_252[k]
                  + f_79 * ab_x[k] * gh_255[k]
                  + f_89 * ab_x[k] * gh_257[k]
                  + f_78 * ab_x[k] * gh_262[k]
                  - f_84 * ab_x[k] * gh_264[k]
                  + f_88 * ab_x[k] * gh_294[k]
                  - f_82 * ab_x[k] * gh_297[k]
                  - f_90 * ab_x[k] * gh_299[k]
                  - f_79 * ab_x[k] * gh_304[k]
                  + f_85 * ab_x[k] * gh_306[k]
                  - f_86 * gi_0[k]
                  + f_80 * gi_3[k]
                  + f_88 * gi_5[k]
                  + f_76 * gi_10[k]
                  - f_79 * gi_12[k]
                  - f_80 * gi_84[k]
                  + f_81 * gi_87[k]
                  + f_82 * gi_89[k]
                  + f_77 * gi_94[k]
                  - f_83 * gi_96[k]
                  + f_87 * gi_140[k]
                  - f_79 * gi_143[k]
                  - f_89 * gi_145[k]
                  - f_78 * gi_150[k]
                  + f_84 * gi_152[k]
                  - f_86 * gi_280[k]
                  + f_80 * gi_283[k]
                  + f_88 * gi_285[k]
                  + f_76 * gi_290[k]
                  - f_79 * gi_292[k]
                  + f_87 * gi_336[k]
                  - f_79 * gi_339[k]
                  - f_89 * gi_341[k]
                  - f_78 * gi_346[k]
                  + f_84 * gi_348[k]
                  - f_88 * gi_392[k]
                  + f_82 * gi_395[k]
                  + f_90 * gi_397[k]
                  + f_79 * gi_402[k]
                  - f_85 * gi_404[k];
        g_94[k] = g_74[k];
    }

#pragma omp simd aligned(ab_x, gh_2, gh_7, gh_16, gh_65, gh_70, gh_79, gh_107, gh_112, gh_121, \
                         gh_212, gh_217, gh_226, gh_254, gh_259, gh_268, gh_296, gh_301, \
                         gh_310, gi_2, gi_7, gi_16, gi_86, gi_91, gi_100, gi_142, gi_147, \
                         gi_156, gi_282, gi_287, gi_296, gi_338, gi_343, gi_352, gi_394, \
                         gi_399, gi_408 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_154 * ab_x[k] * gh_2[k]
                  + f_155 * ab_x[k] * gh_7[k]
                  - f_154 * ab_x[k] * gh_16[k]
                  - f_156 * ab_x[k] * gh_65[k]
                  + f_157 * ab_x[k] * gh_70[k]
                  - f_156 * ab_x[k] * gh_79[k]
                  + f_157 * ab_x[k] * gh_107[k]
                  - f_158 * ab_x[k] * gh_112[k]
                  + f_157 * ab_x[k] * gh_121[k]
                  - f_154 * ab_x[k] * gh_212[k]
                  + f_155 * ab_x[k] * gh_217[k]
                  - f_154 * ab_x[k] * gh_226[k]
                  + f_157 * ab_x[k] * gh_254[k]
                  - f_158 * ab_x[k] * gh_259[k]
                  + f_157 * ab_x[k] * gh_268[k]
                  - f_59 * ab_x[k] * gh_296[k]
                  + f_60 * ab_x[k] * gh_301[k]
                  - f_59 * ab_x[k] * gh_310[k]
                  + f_154 * gi_2[k]
                  - f_155 * gi_7[k]
                  + f_154 * gi_16[k]
                  + f_156 * gi_86[k]
                  - f_157 * gi_91[k]
                  + f_156 * gi_100[k]
                  - f_157 * gi_142[k]
                  + f_158 * gi_147[k]
                  - f_157 * gi_156[k]
                  + f_154 * gi_282[k]
                  - f_155 * gi_287[k]
                  + f_154 * gi_296[k]
                  - f_157 * gi_338[k]
                  + f_158 * gi_343[k]
                  - f_157 * gi_352[k]
                  + f_59 * gi_394[k]
                  - f_60 * gi_399[k]
                  + f_59 * gi_408[k];
        g_105[k] = g_75[k];
    }

#pragma omp simd aligned(ab_x, gh_0, gh_3, gh_10, gh_63, gh_66, gh_73, gh_105, gh_108, gh_115, \
                         gh_210, gh_213, gh_220, gh_252, gh_255, gh_262, gh_294, gh_297, \
                         gh_304, gi_0, gi_3, gi_10, gi_84, gi_87, gi_94, gi_140, gi_143, \
                         gi_150, gi_280, gi_283, gi_290, gi_336, gi_339, gi_346, gi_392, \
                         gi_395, gi_402 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_29 * ab_x[k] * gh_0[k]
                  + f_23 * ab_x[k] * gh_3[k]
                  - f_22 * ab_x[k] * gh_10[k]
                  - f_30 * ab_x[k] * gh_63[k]
                  + f_26 * ab_x[k] * gh_66[k]
                  - f_23 * ab_x[k] * gh_73[k]
                  + f_31 * ab_x[k] * gh_105[k]
                  - f_27 * ab_x[k] * gh_108[k]
                  + f_24 * ab_x[k] * gh_115[k]
                  - f_29 * ab_x[k] * gh_210[k]
                  + f_23 * ab_x[k] * gh_213[k]
                  - f_22 * ab_x[k] * gh_220[k]
                  + f_31 * ab_x[k] * gh_252[k]
                  - f_27 * ab_x[k] * gh_255[k]
                  + f_24 * ab_x[k] * gh_262[k]
                  - f_32 * ab_x[k] * gh_294[k]
                  + f_28 * ab_x[k] * gh_297[k]
                  - f_25 * ab_x[k] * gh_304[k]
                  + f_29 * gi_0[k]
                  - f_23 * gi_3[k]
                  + f_22 * gi_10[k]
                  + f_30 * gi_84[k]
                  - f_26 * gi_87[k]
                  + f_23 * gi_94[k]
                  - f_31 * gi_140[k]
                  + f_27 * gi_143[k]
                  - f_24 * gi_150[k]
                  + f_29 * gi_280[k]
                  - f_23 * gi_283[k]
                  + f_22 * gi_290[k]
                  - f_31 * gi_336[k]
                  + f_27 * gi_339[k]
                  - f_24 * gi_346[k]
                  + f_32 * gi_392[k]
                  - f_28 * gi_395[k]
                  + f_25 * gi_402[k];
        g_116[k] = g_76[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_44, gh_51, gh_58, gh_60, gh_191, gh_198, gh_205, \
                         gh_207, gh_233, gh_240, gh_247, gh_249, gh_275, gh_282, gh_289, \
                         gh_291, gi_58, gi_65, gi_72, gi_74, gi_254, gi_261, gi_268, gi_270, \
                         gi_312, gi_321, gi_330, gi_332, gi_368, gi_377, gi_386, \
                         gi_388 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = -6.5625 * ab_x[k] * gh_44[k]
                  + 13.125 * ab_x[k] * gh_51[k]
                  + 6.5625 * ab_x[k] * gh_58[k]
                  - 13.125 * ab_x[k] * gh_60[k]
                  + 13.125 * ab_x[k] * gh_191[k]
                  - 26.25 * ab_x[k] * gh_198[k]
                  - 13.125 * ab_x[k] * gh_205[k]
                  + 26.25 * ab_x[k] * gh_207[k]
                  + 6.5625 * ab_y[k] * gh_233[k]
                  - 13.125 * ab_y[k] * gh_240[k]
                  - 6.5625 * ab_y[k] * gh_247[k]
                  + 13.125 * ab_y[k] * gh_249[k]
                  - 13.125 * ab_y[k] * gh_275[k]
                  + 26.25 * ab_y[k] * gh_282[k]
                  + 13.125 * ab_y[k] * gh_289[k]
                  - 26.25 * ab_y[k] * gh_291[k]
                  + 6.5625 * gi_58[k]
                  - 13.125 * gi_65[k]
                  - 6.5625 * gi_72[k]
                  + 13.125 * gi_74[k]
                  - 13.125 * gi_254[k]
                  + 26.25 * gi_261[k]
                  + 13.125 * gi_268[k]
                  - 26.25 * gi_270[k]
                  - 6.5625 * gi_312[k]
                  + 13.125 * gi_321[k]
                  + 6.5625 * gi_330[k]
                  - 13.125 * gi_332[k]
                  + 13.125 * gi_368[k]
                  - 26.25 * gi_377[k]
                  - 13.125 * gi_386[k]
                  + 26.25 * gi_388[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_42, gh_45, gh_47, gh_52, gh_54, gh_189, gh_192, \
                         gh_194, gh_199, gh_201, gh_231, gh_234, gh_236, gh_241, gh_243, \
                         gh_273, gh_276, gh_278, gh_283, gh_285, gi_56, gi_59, gi_61, gi_66, \
                         gi_68, gi_252, gi_255, gi_257, gi_262, gi_264, gi_309, gi_314, \
                         gi_316, gi_323, gi_325, gi_365, gi_370, gi_372, gi_379, \
                         gi_381 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = -f_111 * ab_x[k] * gh_42[k]
                  + f_73 * ab_x[k] * gh_45[k]
                  + f_70 * ab_x[k] * gh_47[k]
                  + f_109 * ab_x[k] * gh_52[k]
                  - f_110 * ab_x[k] * gh_54[k]
                  + f_73 * ab_x[k] * gh_189[k]
                  - f_69 * ab_x[k] * gh_192[k]
                  - f_74 * ab_x[k] * gh_194[k]
                  - f_67 * ab_x[k] * gh_199[k]
                  + f_71 * ab_x[k] * gh_201[k]
                  + f_111 * ab_y[k] * gh_231[k]
                  - f_73 * ab_y[k] * gh_234[k]
                  - f_70 * ab_y[k] * gh_236[k]
                  - f_109 * ab_y[k] * gh_241[k]
                  + f_110 * ab_y[k] * gh_243[k]
                  - f_73 * ab_y[k] * gh_273[k]
                  + f_69 * ab_y[k] * gh_276[k]
                  + f_74 * ab_y[k] * gh_278[k]
                  + f_67 * ab_y[k] * gh_283[k]
                  - f_71 * ab_y[k] * gh_285[k]
                  + f_111 * gi_56[k]
                  - f_73 * gi_59[k]
                  - f_70 * gi_61[k]
                  - f_109 * gi_66[k]
                  + f_110 * gi_68[k]
                  - f_73 * gi_252[k]
                  + f_69 * gi_255[k]
                  + f_74 * gi_257[k]
                  + f_67 * gi_262[k]
                  - f_71 * gi_264[k]
                  - f_111 * gi_309[k]
                  + f_73 * gi_314[k]
                  + f_70 * gi_316[k]
                  + f_109 * gi_323[k]
                  - f_110 * gi_325[k]
                  + f_73 * gi_365[k]
                  - f_69 * gi_370[k]
                  - f_74 * gi_372[k]
                  - f_67 * gi_379[k]
                  + f_71 * gi_381[k];
        g_95[k] = g_85[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_44, gh_49, gh_58, gh_191, gh_196, gh_205, gh_233, \
                         gh_238, gh_247, gh_275, gh_280, gh_289, gi_58, gi_63, gi_72, gi_254, \
                         gi_259, gi_268, gi_312, gi_319, gi_330, gi_368, gi_375, \
                         gi_386 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = f_170 * ab_x[k] * gh_44[k]
                  - f_171 * ab_x[k] * gh_49[k]
                  + f_170 * ab_x[k] * gh_58[k]
                  - f_134 * ab_x[k] * gh_191[k]
                  + f_135 * ab_x[k] * gh_196[k]
                  - f_134 * ab_x[k] * gh_205[k]
                  - f_170 * ab_y[k] * gh_233[k]
                  + f_171 * ab_y[k] * gh_238[k]
                  - f_170 * ab_y[k] * gh_247[k]
                  + f_134 * ab_y[k] * gh_275[k]
                  - f_135 * ab_y[k] * gh_280[k]
                  + f_134 * ab_y[k] * gh_289[k]
                  - f_170 * gi_58[k]
                  + f_171 * gi_63[k]
                  - f_170 * gi_72[k]
                  + f_134 * gi_254[k]
                  - f_135 * gi_259[k]
                  + f_134 * gi_268[k]
                  + f_170 * gi_312[k]
                  - f_171 * gi_319[k]
                  + f_170 * gi_330[k]
                  - f_134 * gi_368[k]
                  + f_135 * gi_375[k]
                  - f_134 * gi_386[k];
        g_106[k] = g_86[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_42, gh_45, gh_52, gh_189, gh_192, gh_199, gh_231, \
                         gh_234, gh_241, gh_273, gh_276, gh_283, gi_56, gi_59, gi_66, gi_252, \
                         gi_255, gi_262, gi_309, gi_314, gi_323, gi_365, gi_370, \
                         gi_379 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = f_44 * ab_x[k] * gh_42[k]
                  - f_17 * ab_x[k] * gh_45[k]
                  + f_43 * ab_x[k] * gh_52[k]
                  - f_20 * ab_x[k] * gh_189[k]
                  + f_18 * ab_x[k] * gh_192[k]
                  - f_17 * ab_x[k] * gh_199[k]
                  - f_44 * ab_y[k] * gh_231[k]
                  + f_17 * ab_y[k] * gh_234[k]
                  - f_43 * ab_y[k] * gh_241[k]
                  + f_20 * ab_y[k] * gh_273[k]
                  - f_18 * ab_y[k] * gh_276[k]
                  + f_17 * ab_y[k] * gh_283[k]
                  - f_44 * gi_56[k]
                  + f_17 * gi_59[k]
                  - f_43 * gi_66[k]
                  + f_20 * gi_252[k]
                  - f_18 * gi_255[k]
                  + f_17 * gi_262[k]
                  + f_44 * gi_309[k]
                  - f_17 * gi_314[k]
                  + f_43 * gi_323[k]
                  - f_20 * gi_365[k]
                  + f_18 * gi_370[k]
                  - f_17 * gi_379[k];
        g_117[k] = g_87[k];
    }

#pragma omp simd aligned(ab_x, gh_0, gh_3, gh_5, gh_10, gh_12, gh_63, gh_66, gh_68, gh_73, \
                         gh_75, gh_105, gh_108, gh_110, gh_115, gh_117, gh_210, gh_213, \
                         gh_215, gh_220, gh_222, gh_252, gh_255, gh_257, gh_262, gh_264, gi_0, \
                         gi_3, gi_5, gi_10, gi_12, gi_84, gi_87, gi_89, gi_94, gi_96, gi_140, \
                         gi_143, gi_145, gi_150, gi_152, gi_280, gi_283, gi_285, gi_290, \
                         gi_292, gi_336, gi_339, gi_341, gi_346, \
                         gi_348 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = -0.2734375 * ab_x[k] * gh_0[k]
                  + 0.546875 * ab_x[k] * gh_3[k]
                  + 2.1875 * ab_x[k] * gh_5[k]
                  + 0.8203125 * ab_x[k] * gh_10[k]
                  - 6.5625 * ab_x[k] * gh_12[k]
                  + 0.546875 * ab_x[k] * gh_63[k]
                  - 1.09375 * ab_x[k] * gh_66[k]
                  - 4.375 * ab_x[k] * gh_68[k]
                  - 1.640625 * ab_x[k] * gh_73[k]
                  + 13.125 * ab_x[k] * gh_75[k]
                  + 2.1875 * ab_x[k] * gh_105[k]
                  - 4.375 * ab_x[k] * gh_108[k]
                  - 17.5 * ab_x[k] * gh_110[k]
                  - 6.5625 * ab_x[k] * gh_115[k]
                  + 52.5 * ab_x[k] * gh_117[k]
                  + 0.8203125 * ab_x[k] * gh_210[k]
                  - 1.640625 * ab_x[k] * gh_213[k]
                  - 6.5625 * ab_x[k] * gh_215[k]
                  - 2.4609375 * ab_x[k] * gh_220[k]
                  + 19.6875 * ab_x[k] * gh_222[k]
                  - 6.5625 * ab_x[k] * gh_252[k]
                  + 13.125 * ab_x[k] * gh_255[k]
                  + 52.5 * ab_x[k] * gh_257[k]
                  + 19.6875 * ab_x[k] * gh_262[k]
                  - 157.5 * ab_x[k] * gh_264[k]
                  + 0.2734375 * gi_0[k]
                  - 0.546875 * gi_3[k]
                  - 2.1875 * gi_5[k]
                  - 0.8203125 * gi_10[k]
                  + 6.5625 * gi_12[k]
                  - 0.546875 * gi_84[k]
                  + 1.09375 * gi_87[k]
                  + 4.375 * gi_89[k]
                  + 1.640625 * gi_94[k]
                  - 13.125 * gi_96[k]
                  - 2.1875 * gi_140[k]
                  + 4.375 * gi_143[k]
                  + 17.5 * gi_145[k]
                  + 6.5625 * gi_150[k]
                  - 52.5 * gi_152[k]
                  - 0.8203125 * gi_280[k]
                  + 1.640625 * gi_283[k]
                  + 6.5625 * gi_285[k]
                  + 2.4609375 * gi_290[k]
                  - 19.6875 * gi_292[k]
                  + 6.5625 * gi_336[k]
                  - 13.125 * gi_339[k]
                  - 52.5 * gi_341[k]
                  - 19.6875 * gi_346[k]
                  + 157.5 * gi_348[k];
    }

#pragma omp simd aligned(ab_x, gh_2, gh_7, gh_16, gh_65, gh_70, gh_79, gh_107, gh_112, gh_121, \
                         gh_212, gh_217, gh_226, gh_254, gh_259, gh_268, gi_2, gi_7, gi_16, \
                         gi_86, gi_91, gi_100, gi_142, gi_147, gi_156, gi_282, gi_287, gi_296, \
                         gi_338, gi_343, gi_352 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = f_117 * ab_x[k] * gh_2[k]
                  - f_118 * ab_x[k] * gh_7[k]
                  + f_117 * ab_x[k] * gh_16[k]
                  - f_114 * ab_x[k] * gh_65[k]
                  + f_51 * ab_x[k] * gh_70[k]
                  - f_114 * ab_x[k] * gh_79[k]
                  - f_52 * ab_x[k] * gh_107[k]
                  + f_119 * ab_x[k] * gh_112[k]
                  - f_52 * ab_x[k] * gh_121[k]
                  - f_112 * ab_x[k] * gh_212[k]
                  + f_113 * ab_x[k] * gh_217[k]
                  - f_112 * ab_x[k] * gh_226[k]
                  + f_115 * ab_x[k] * gh_254[k]
                  - f_116 * ab_x[k] * gh_259[k]
                  + f_115 * ab_x[k] * gh_268[k]
                  - f_117 * gi_2[k]
                  + f_118 * gi_7[k]
                  - f_117 * gi_16[k]
                  + f_114 * gi_86[k]
                  - f_51 * gi_91[k]
                  + f_114 * gi_100[k]
                  + f_52 * gi_142[k]
                  - f_119 * gi_147[k]
                  + f_52 * gi_156[k]
                  + f_112 * gi_282[k]
                  - f_113 * gi_287[k]
                  + f_112 * gi_296[k]
                  - f_115 * gi_338[k]
                  + f_116 * gi_343[k]
                  - f_115 * gi_352[k];
        g_107[k] = g_97[k];
    }

#pragma omp simd aligned(ab_x, gh_0, gh_3, gh_10, gh_63, gh_66, gh_73, gh_105, gh_108, gh_115, \
                         gh_210, gh_213, gh_220, gh_252, gh_255, gh_262, gi_0, gi_3, gi_10, \
                         gi_84, gi_87, gi_94, gi_140, gi_143, gi_150, gi_280, gi_283, gi_290, \
                         gi_336, gi_339, gi_346 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = f_15 * ab_x[k] * gh_0[k]
                  - f_4 * ab_x[k] * gh_3[k]
                  + f_6 * ab_x[k] * gh_10[k]
                  - f_13 * ab_x[k] * gh_63[k]
                  + f_9 * ab_x[k] * gh_66[k]
                  - f_4 * ab_x[k] * gh_73[k]
                  - f_16 * ab_x[k] * gh_105[k]
                  + f_11 * ab_x[k] * gh_108[k]
                  - f_7 * ab_x[k] * gh_115[k]
                  - f_12 * ab_x[k] * gh_210[k]
                  + f_8 * ab_x[k] * gh_213[k]
                  - f_3 * ab_x[k] * gh_220[k]
                  + f_14 * ab_x[k] * gh_252[k]
                  - f_10 * ab_x[k] * gh_255[k]
                  + f_5 * ab_x[k] * gh_262[k]
                  - f_15 * gi_0[k]
                  + f_4 * gi_3[k]
                  - f_6 * gi_10[k]
                  + f_13 * gi_84[k]
                  - f_9 * gi_87[k]
                  + f_4 * gi_94[k]
                  + f_16 * gi_140[k]
                  - f_11 * gi_143[k]
                  + f_7 * gi_150[k]
                  + f_12 * gi_280[k]
                  - f_8 * gi_283[k]
                  + f_3 * gi_290[k]
                  - f_14 * gi_336[k]
                  + f_10 * gi_339[k]
                  - f_5 * gi_346[k];
        g_118[k] = g_98[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_44, gh_49, gh_58, gh_149, gh_154, gh_163, gh_233, \
                         gh_238, gh_247, gi_58, gi_63, gi_72, gi_198, gi_203, gi_212, gi_312, \
                         gi_319, gi_330 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = -4.921875 * ab_x[k] * gh_44[k]
                   + 29.53125 * ab_x[k] * gh_49[k]
                   - 4.921875 * ab_x[k] * gh_58[k]
                   + 29.53125 * ab_x[k] * gh_149[k]
                   - 177.1875 * ab_x[k] * gh_154[k]
                   + 29.53125 * ab_x[k] * gh_163[k]
                   - 4.921875 * ab_y[k] * gh_233[k]
                   + 29.53125 * ab_y[k] * gh_238[k]
                   - 4.921875 * ab_y[k] * gh_247[k]
                   + 4.921875 * gi_58[k]
                   - 29.53125 * gi_63[k]
                   + 4.921875 * gi_72[k]
                   - 29.53125 * gi_198[k]
                   + 177.1875 * gi_203[k]
                   - 29.53125 * gi_212[k]
                   + 4.921875 * gi_312[k]
                   - 29.53125 * gi_319[k]
                   + 4.921875 * gi_330[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gh_42, gh_45, gh_52, gh_147, gh_150, gh_157, gh_231, \
                         gh_234, gh_241, gi_56, gi_59, gi_66, gi_196, gi_199, gi_206, gi_309, \
                         gi_314, gi_323 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = -f_49 * ab_x[k] * gh_42[k]
                   + f_47 * ab_x[k] * gh_45[k]
                   - f_45 * ab_x[k] * gh_52[k]
                   + f_50 * ab_x[k] * gh_147[k]
                   - f_48 * ab_x[k] * gh_150[k]
                   + f_46 * ab_x[k] * gh_157[k]
                   - f_49 * ab_y[k] * gh_231[k]
                   + f_47 * ab_y[k] * gh_234[k]
                   - f_45 * ab_y[k] * gh_241[k]
                   + f_49 * gi_56[k]
                   - f_47 * gi_59[k]
                   + f_45 * gi_66[k]
                   - f_50 * gi_196[k]
                   + f_48 * gi_199[k]
                   - f_46 * gi_206[k]
                   + f_49 * gi_309[k]
                   - f_47 * gi_314[k]
                   + f_45 * gi_323[k];
        g_119[k] = g_109[k];
    }

#pragma omp simd aligned(ab_x, gh_0, gh_3, gh_10, gh_63, gh_66, gh_73, gh_210, gh_213, gh_220, \
                         gi_0, gi_3, gi_10, gi_84, gi_87, gi_94, gi_280, gi_283, \
                         gi_290 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_120[k] = -0.4921875 * ab_x[k] * gh_0[k]
                   + 4.921875 * ab_x[k] * gh_3[k]
                   - 2.4609375 * ab_x[k] * gh_10[k]
                   + 4.921875 * ab_x[k] * gh_63[k]
                   - 49.21875 * ab_x[k] * gh_66[k]
                   + 24.609375 * ab_x[k] * gh_73[k]
                   - 2.4609375 * ab_x[k] * gh_210[k]
                   + 24.609375 * ab_x[k] * gh_213[k]
                   - 12.3046875 * ab_x[k] * gh_220[k]
                   + 0.4921875 * gi_0[k]
                   - 4.921875 * gi_3[k]
                   + 2.4609375 * gi_10[k]
                   - 4.921875 * gi_84[k]
                   + 49.21875 * gi_87[k]
                   - 24.609375 * gi_94[k]
                   + 2.4609375 * gi_280[k]
                   - 24.609375 * gi_283[k]
                   + 12.3046875 * gi_290[k];
    }
}

}  // namespace simdtrf
