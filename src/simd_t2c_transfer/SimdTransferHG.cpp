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


#include "SimdTransferHG.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_hg_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t hf, const size_t if_,
                   const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 3.28125 * std::sqrt(10.0);
    const auto f_1 = 6.5625 * std::sqrt(10.0);
    const auto f_2 = 0.65625 * std::sqrt(10.0);
    const auto f_3 = 9.84375 * std::sqrt(5.0);
    const auto f_4 = 3.28125 * std::sqrt(5.0);
    const auto f_5 = 19.6875 * std::sqrt(5.0);
    const auto f_6 = 6.5625 * std::sqrt(5.0);
    const auto f_7 = 1.96875 * std::sqrt(5.0);
    const auto f_8 = 0.65625 * std::sqrt(5.0);
    const auto f_9 = 0.46875 * std::sqrt(70.0);
    const auto f_10 = 2.8125 * std::sqrt(70.0);
    const auto f_11 = 0.9375 * std::sqrt(70.0);
    const auto f_12 = 5.625 * std::sqrt(70.0);
    const auto f_13 = 0.09375 * std::sqrt(70.0);
    const auto f_14 = 0.5625 * std::sqrt(70.0);
    const auto f_15 = 1.40625 * std::sqrt(35.0);
    const auto f_16 = 1.875 * std::sqrt(35.0);
    const auto f_17 = 2.8125 * std::sqrt(35.0);
    const auto f_18 = 3.75 * std::sqrt(35.0);
    const auto f_19 = 0.28125 * std::sqrt(35.0);
    const auto f_20 = 0.375 * std::sqrt(35.0);
    const auto f_21 = 0.3515625 * std::sqrt(14.0);
    const auto f_22 = 0.703125 * std::sqrt(14.0);
    const auto f_23 = 2.8125 * std::sqrt(14.0);
    const auto f_24 = 0.9375 * std::sqrt(14.0);
    const auto f_25 = 1.40625 * std::sqrt(14.0);
    const auto f_26 = 5.625 * std::sqrt(14.0);
    const auto f_27 = 1.875 * std::sqrt(14.0);
    const auto f_28 = 0.0703125 * std::sqrt(14.0);
    const auto f_29 = 0.140625 * std::sqrt(14.0);
    const auto f_30 = 0.5625 * std::sqrt(14.0);
    const auto f_31 = 0.1875 * std::sqrt(14.0);
    const auto f_32 = 0.234375 * std::sqrt(70.0);
    const auto f_33 = 1.40625 * std::sqrt(70.0);
    const auto f_34 = 0.046875 * std::sqrt(70.0);
    const auto f_35 = 0.28125 * std::sqrt(70.0);
    const auto f_36 = 0.8203125 * std::sqrt(10.0);
    const auto f_37 = 4.921875 * std::sqrt(10.0);
    const auto f_38 = 1.640625 * std::sqrt(10.0);
    const auto f_39 = 9.84375 * std::sqrt(10.0);
    const auto f_40 = 0.1640625 * std::sqrt(10.0);
    const auto f_41 = 0.984375 * std::sqrt(10.0);
    const auto f_42 = 39.375 * std::sqrt(2.0);
    const auto f_43 = 13.125 * std::sqrt(2.0);
    const auto f_44 = 3.75 * std::sqrt(7.0);
    const auto f_45 = 22.5 * std::sqrt(7.0);
    const auto f_46 = 7.5 * std::sqrt(14.0);
    const auto f_47 = 0.5625 * std::sqrt(35.0);
    const auto f_48 = 1.125 * std::sqrt(35.0);
    const auto f_49 = 4.5 * std::sqrt(35.0);
    const auto f_50 = 1.5 * std::sqrt(35.0);
    const auto f_51 = 1.875 * std::sqrt(7.0);
    const auto f_52 = 11.25 * std::sqrt(7.0);
    const auto f_53 = 3.28125 * std::sqrt(2.0);
    const auto f_54 = 2.1875 * std::sqrt(2.0);
    const auto f_55 = 26.25 * std::sqrt(2.0);
    const auto f_56 = 1.09375 * std::sqrt(2.0);
    const auto f_57 = 8.75 * std::sqrt(2.0);
    const auto f_58 = 0.46875 * std::sqrt(14.0);
    const auto f_59 = 0.3125 * std::sqrt(14.0);
    const auto f_60 = 3.75 * std::sqrt(14.0);
    const auto f_61 = 22.5 * std::sqrt(14.0);
    const auto f_62 = 0.15625 * std::sqrt(14.0);
    const auto f_63 = 1.25 * std::sqrt(14.0);
    const auto f_64 = 1.40625 * std::sqrt(7.0);
    const auto f_65 = 0.9375 * std::sqrt(7.0);
    const auto f_66 = 1.25 * std::sqrt(7.0);
    const auto f_67 = 15.0 * std::sqrt(7.0);
    const auto f_68 = 0.46875 * std::sqrt(7.0);
    const auto f_69 = 0.625 * std::sqrt(7.0);
    const auto f_70 = 5.0 * std::sqrt(7.0);
    const auto f_71 = 0.0703125 * std::sqrt(70.0);
    const auto f_72 = 0.140625 * std::sqrt(70.0);
    const auto f_73 = 0.1875 * std::sqrt(70.0);
    const auto f_74 = 0.375 * std::sqrt(70.0);
    const auto f_75 = 0.125 * std::sqrt(70.0);
    const auto f_76 = 1.125 * std::sqrt(70.0);
    const auto f_77 = 4.5 * std::sqrt(70.0);
    const auto f_78 = 1.5 * std::sqrt(70.0);
    const auto f_79 = 0.0234375 * std::sqrt(70.0);
    const auto f_80 = 0.0625 * std::sqrt(70.0);
    const auto f_81 = 0.5 * std::sqrt(70.0);
    const auto f_82 = 0.234375 * std::sqrt(14.0);
    const auto f_83 = 11.25 * std::sqrt(14.0);
    const auto f_84 = 0.078125 * std::sqrt(14.0);
    const auto f_85 = 0.625 * std::sqrt(14.0);
    const auto f_86 = 0.8203125 * std::sqrt(2.0);
    const auto f_87 = 4.921875 * std::sqrt(2.0);
    const auto f_88 = 0.546875 * std::sqrt(2.0);
    const auto f_89 = 6.5625 * std::sqrt(2.0);
    const auto f_90 = 0.2734375 * std::sqrt(2.0);
    const auto f_91 = 1.640625 * std::sqrt(2.0);
    const auto f_92 = 8.75 * std::sqrt(3.0);
    const auto f_93 = 17.5 * std::sqrt(3.0);
    const auto f_94 = 13.125 * std::sqrt(6.0);
    const auto f_95 = 4.375 * std::sqrt(6.0);
    const auto f_96 = 26.25 * std::sqrt(6.0);
    const auto f_97 = 8.75 * std::sqrt(6.0);
    const auto f_98 = 1.25 * std::sqrt(21.0);
    const auto f_99 = 7.5 * std::sqrt(21.0);
    const auto f_100 = 2.5 * std::sqrt(21.0);
    const auto f_101 = 15.0 * std::sqrt(21.0);
    const auto f_102 = 1.875 * std::sqrt(42.0);
    const auto f_103 = 2.5 * std::sqrt(42.0);
    const auto f_104 = 3.75 * std::sqrt(42.0);
    const auto f_105 = 5.0 * std::sqrt(42.0);
    const auto f_106 = 0.1875 * std::sqrt(105.0);
    const auto f_107 = 0.375 * std::sqrt(105.0);
    const auto f_108 = 1.5 * std::sqrt(105.0);
    const auto f_109 = 0.5 * std::sqrt(105.0);
    const auto f_110 = 0.75 * std::sqrt(105.0);
    const auto f_111 = 3.0 * std::sqrt(105.0);
    const auto f_112 = std::sqrt(105.0);
    const auto f_113 = 0.625 * std::sqrt(21.0);
    const auto f_114 = 3.75 * std::sqrt(21.0);
    const auto f_115 = 2.1875 * std::sqrt(3.0);
    const auto f_116 = 13.125 * std::sqrt(3.0);
    const auto f_117 = 4.375 * std::sqrt(3.0);
    const auto f_118 = 26.25 * std::sqrt(3.0);
    const auto f_119 = 0.3125 * std::sqrt(21.0);
    const auto f_120 = 0.46875 * std::sqrt(42.0);
    const auto f_121 = 0.15625 * std::sqrt(42.0);
    const auto f_122 = 0.9375 * std::sqrt(42.0);
    const auto f_123 = 0.3125 * std::sqrt(42.0);
    const auto f_124 = 5.625 * std::sqrt(42.0);
    const auto f_125 = 1.25 * std::sqrt(42.0);
    const auto f_126 = 0.3125 * std::sqrt(3.0);
    const auto f_127 = 1.875 * std::sqrt(3.0);
    const auto f_128 = 0.625 * std::sqrt(3.0);
    const auto f_129 = 3.75 * std::sqrt(3.0);
    const auto f_130 = 22.5 * std::sqrt(3.0);
    const auto f_131 = 2.5 * std::sqrt(3.0);
    const auto f_132 = 15.0 * std::sqrt(3.0);
    const auto f_133 = 0.46875 * std::sqrt(6.0);
    const auto f_134 = 0.625 * std::sqrt(6.0);
    const auto f_135 = 0.9375 * std::sqrt(6.0);
    const auto f_136 = 1.25 * std::sqrt(6.0);
    const auto f_137 = 5.625 * std::sqrt(6.0);
    const auto f_138 = 7.5 * std::sqrt(6.0);
    const auto f_139 = 3.75 * std::sqrt(6.0);
    const auto f_140 = 5.0 * std::sqrt(6.0);
    const auto f_141 = 0.046875 * std::sqrt(15.0);
    const auto f_142 = 0.09375 * std::sqrt(15.0);
    const auto f_143 = 0.375 * std::sqrt(15.0);
    const auto f_144 = 0.125 * std::sqrt(15.0);
    const auto f_145 = 0.1875 * std::sqrt(15.0);
    const auto f_146 = 0.75 * std::sqrt(15.0);
    const auto f_147 = 0.25 * std::sqrt(15.0);
    const auto f_148 = 0.5625 * std::sqrt(15.0);
    const auto f_149 = 1.125 * std::sqrt(15.0);
    const auto f_150 = 4.5 * std::sqrt(15.0);
    const auto f_151 = 1.5 * std::sqrt(15.0);
    const auto f_152 = 3.0 * std::sqrt(15.0);
    const auto f_153 = std::sqrt(15.0);
    const auto f_154 = 0.15625 * std::sqrt(3.0);
    const auto f_155 = 0.9375 * std::sqrt(3.0);
    const auto f_156 = 11.25 * std::sqrt(3.0);
    const auto f_157 = 1.25 * std::sqrt(3.0);
    const auto f_158 = 7.5 * std::sqrt(3.0);
    const auto f_159 = 0.078125 * std::sqrt(21.0);
    const auto f_160 = 0.46875 * std::sqrt(21.0);
    const auto f_161 = 0.15625 * std::sqrt(21.0);
    const auto f_162 = 0.9375 * std::sqrt(21.0);
    const auto f_163 = 5.625 * std::sqrt(21.0);
    const auto f_164 = 0.9375 * std::sqrt(35.0);
    const auto f_165 = 2.5 * std::sqrt(35.0);
    const auto f_166 = 0.5 * std::sqrt(35.0);
    const auto f_167 = 3.75 * std::sqrt(70.0);
    const auto f_168 = 1.25 * std::sqrt(70.0);
    const auto f_169 = 0.75 * std::sqrt(70.0);
    const auto f_170 = 0.25 * std::sqrt(70.0);
    const auto f_171 = 0.9375 * std::sqrt(5.0);
    const auto f_172 = 5.625 * std::sqrt(5.0);
    const auto f_173 = 1.875 * std::sqrt(5.0);
    const auto f_174 = 11.25 * std::sqrt(5.0);
    const auto f_175 = 2.5 * std::sqrt(5.0);
    const auto f_176 = 15.0 * std::sqrt(5.0);
    const auto f_177 = 0.5 * std::sqrt(5.0);
    const auto f_178 = 3.0 * std::sqrt(5.0);
    const auto f_179 = 1.40625 * std::sqrt(10.0);
    const auto f_180 = 1.875 * std::sqrt(10.0);
    const auto f_181 = 2.8125 * std::sqrt(10.0);
    const auto f_182 = 3.75 * std::sqrt(10.0);
    const auto f_183 = 5.0 * std::sqrt(10.0);
    const auto f_184 = 0.75 * std::sqrt(10.0);
    const auto f_185 = std::sqrt(10.0);
    const auto f_186 = 0.46875 * std::sqrt(5.0);
    const auto f_187 = 2.8125 * std::sqrt(5.0);
    const auto f_188 = 1.25 * std::sqrt(5.0);
    const auto f_189 = 7.5 * std::sqrt(5.0);
    const auto f_190 = 0.25 * std::sqrt(5.0);
    const auto f_191 = 1.5 * std::sqrt(5.0);
    const auto f_192 = 0.234375 * std::sqrt(35.0);
    const auto f_193 = 0.46875 * std::sqrt(35.0);
    const auto f_194 = 0.625 * std::sqrt(35.0);
    const auto f_195 = 0.125 * std::sqrt(35.0);
    const auto f_196 = 0.75 * std::sqrt(35.0);
    const auto f_197 = 6.5625 * std::sqrt(6.0);
    const auto f_198 = 2.1875 * std::sqrt(6.0);
    const auto f_199 = 0.09375 * std::sqrt(105.0);
    const auto f_200 = 0.25 * std::sqrt(105.0);
    const auto f_201 = 1.875 * std::sqrt(21.0);
    const auto f_202 = 1.09375 * std::sqrt(3.0);
    const auto f_203 = 6.5625 * std::sqrt(3.0);
    const auto f_204 = 9.84375 * std::sqrt(2.0);
    const auto f_205 = 59.0625 * std::sqrt(2.0);
    const auto f_206 = 19.6875 * std::sqrt(2.0);
    const auto f_207 = 5.625 * std::sqrt(7.0);
    const auto f_208 = 33.75 * std::sqrt(7.0);
    const auto f_209 = 8.4375 * std::sqrt(14.0);
    const auto f_210 = 0.140625 * std::sqrt(35.0);
    const auto f_211 = 0.84375 * std::sqrt(35.0);
    const auto f_212 = 1.6875 * std::sqrt(35.0);
    const auto f_213 = 6.75 * std::sqrt(35.0);
    const auto f_214 = 2.25 * std::sqrt(35.0);
    const auto f_215 = 2.8125 * std::sqrt(7.0);
    const auto f_216 = 16.875 * std::sqrt(7.0);

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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_70 = buffer.data(hf + 70);
    const auto *hf_71 = buffer.data(hf + 71);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_74 = buffer.data(hf + 74);
    const auto *hf_75 = buffer.data(hf + 75);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_81 = buffer.data(hf + 81);
    const auto *hf_82 = buffer.data(hf + 82);
    const auto *hf_83 = buffer.data(hf + 83);
    const auto *hf_84 = buffer.data(hf + 84);
    const auto *hf_85 = buffer.data(hf + 85);
    const auto *hf_86 = buffer.data(hf + 86);
    const auto *hf_87 = buffer.data(hf + 87);
    const auto *hf_88 = buffer.data(hf + 88);
    const auto *hf_89 = buffer.data(hf + 89);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_93 = buffer.data(hf + 93);
    const auto *hf_94 = buffer.data(hf + 94);
    const auto *hf_95 = buffer.data(hf + 95);
    const auto *hf_96 = buffer.data(hf + 96);
    const auto *hf_97 = buffer.data(hf + 97);
    const auto *hf_98 = buffer.data(hf + 98);
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_100 = buffer.data(hf + 100);
    const auto *hf_101 = buffer.data(hf + 101);
    const auto *hf_102 = buffer.data(hf + 102);
    const auto *hf_103 = buffer.data(hf + 103);
    const auto *hf_104 = buffer.data(hf + 104);
    const auto *hf_105 = buffer.data(hf + 105);
    const auto *hf_106 = buffer.data(hf + 106);
    const auto *hf_107 = buffer.data(hf + 107);
    const auto *hf_108 = buffer.data(hf + 108);
    const auto *hf_109 = buffer.data(hf + 109);
    const auto *hf_110 = buffer.data(hf + 110);
    const auto *hf_111 = buffer.data(hf + 111);
    const auto *hf_112 = buffer.data(hf + 112);
    const auto *hf_113 = buffer.data(hf + 113);
    const auto *hf_114 = buffer.data(hf + 114);
    const auto *hf_115 = buffer.data(hf + 115);
    const auto *hf_116 = buffer.data(hf + 116);
    const auto *hf_117 = buffer.data(hf + 117);
    const auto *hf_118 = buffer.data(hf + 118);
    const auto *hf_119 = buffer.data(hf + 119);
    const auto *hf_120 = buffer.data(hf + 120);
    const auto *hf_121 = buffer.data(hf + 121);
    const auto *hf_122 = buffer.data(hf + 122);
    const auto *hf_123 = buffer.data(hf + 123);
    const auto *hf_124 = buffer.data(hf + 124);
    const auto *hf_125 = buffer.data(hf + 125);
    const auto *hf_126 = buffer.data(hf + 126);
    const auto *hf_127 = buffer.data(hf + 127);
    const auto *hf_128 = buffer.data(hf + 128);
    const auto *hf_129 = buffer.data(hf + 129);
    const auto *hf_130 = buffer.data(hf + 130);
    const auto *hf_131 = buffer.data(hf + 131);
    const auto *hf_132 = buffer.data(hf + 132);
    const auto *hf_133 = buffer.data(hf + 133);
    const auto *hf_134 = buffer.data(hf + 134);
    const auto *hf_135 = buffer.data(hf + 135);
    const auto *hf_136 = buffer.data(hf + 136);
    const auto *hf_137 = buffer.data(hf + 137);
    const auto *hf_138 = buffer.data(hf + 138);
    const auto *hf_139 = buffer.data(hf + 139);
    const auto *hf_140 = buffer.data(hf + 140);
    const auto *hf_141 = buffer.data(hf + 141);
    const auto *hf_142 = buffer.data(hf + 142);
    const auto *hf_143 = buffer.data(hf + 143);
    const auto *hf_144 = buffer.data(hf + 144);
    const auto *hf_145 = buffer.data(hf + 145);
    const auto *hf_146 = buffer.data(hf + 146);
    const auto *hf_147 = buffer.data(hf + 147);
    const auto *hf_148 = buffer.data(hf + 148);
    const auto *hf_149 = buffer.data(hf + 149);
    const auto *hf_150 = buffer.data(hf + 150);
    const auto *hf_151 = buffer.data(hf + 151);
    const auto *hf_152 = buffer.data(hf + 152);
    const auto *hf_153 = buffer.data(hf + 153);
    const auto *hf_154 = buffer.data(hf + 154);
    const auto *hf_155 = buffer.data(hf + 155);
    const auto *hf_156 = buffer.data(hf + 156);
    const auto *hf_157 = buffer.data(hf + 157);
    const auto *hf_158 = buffer.data(hf + 158);
    const auto *hf_159 = buffer.data(hf + 159);
    const auto *hf_160 = buffer.data(hf + 160);
    const auto *hf_161 = buffer.data(hf + 161);
    const auto *hf_162 = buffer.data(hf + 162);
    const auto *hf_163 = buffer.data(hf + 163);
    const auto *hf_164 = buffer.data(hf + 164);
    const auto *hf_165 = buffer.data(hf + 165);
    const auto *hf_166 = buffer.data(hf + 166);
    const auto *hf_167 = buffer.data(hf + 167);
    const auto *hf_168 = buffer.data(hf + 168);
    const auto *hf_169 = buffer.data(hf + 169);
    const auto *hf_170 = buffer.data(hf + 170);
    const auto *hf_171 = buffer.data(hf + 171);
    const auto *hf_172 = buffer.data(hf + 172);
    const auto *hf_173 = buffer.data(hf + 173);
    const auto *hf_174 = buffer.data(hf + 174);
    const auto *hf_175 = buffer.data(hf + 175);
    const auto *hf_176 = buffer.data(hf + 176);
    const auto *hf_177 = buffer.data(hf + 177);
    const auto *hf_178 = buffer.data(hf + 178);
    const auto *hf_179 = buffer.data(hf + 179);
    const auto *hf_180 = buffer.data(hf + 180);
    const auto *hf_181 = buffer.data(hf + 181);
    const auto *hf_182 = buffer.data(hf + 182);
    const auto *hf_183 = buffer.data(hf + 183);
    const auto *hf_184 = buffer.data(hf + 184);
    const auto *hf_185 = buffer.data(hf + 185);
    const auto *hf_186 = buffer.data(hf + 186);
    const auto *hf_187 = buffer.data(hf + 187);
    const auto *hf_188 = buffer.data(hf + 188);
    const auto *hf_189 = buffer.data(hf + 189);
    const auto *hf_190 = buffer.data(hf + 190);
    const auto *hf_191 = buffer.data(hf + 191);
    const auto *hf_192 = buffer.data(hf + 192);
    const auto *hf_193 = buffer.data(hf + 193);
    const auto *hf_194 = buffer.data(hf + 194);
    const auto *hf_195 = buffer.data(hf + 195);
    const auto *hf_196 = buffer.data(hf + 196);
    const auto *hf_197 = buffer.data(hf + 197);
    const auto *hf_198 = buffer.data(hf + 198);
    const auto *hf_199 = buffer.data(hf + 199);
    const auto *hf_200 = buffer.data(hf + 200);
    const auto *hf_201 = buffer.data(hf + 201);
    const auto *hf_202 = buffer.data(hf + 202);
    const auto *hf_203 = buffer.data(hf + 203);
    const auto *hf_204 = buffer.data(hf + 204);
    const auto *hf_205 = buffer.data(hf + 205);
    const auto *hf_206 = buffer.data(hf + 206);
    const auto *hf_207 = buffer.data(hf + 207);
    const auto *hf_208 = buffer.data(hf + 208);
    const auto *hf_209 = buffer.data(hf + 209);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__1 = buffer.data(if_ + 1);
    const auto *if__2 = buffer.data(if_ + 2);
    const auto *if__3 = buffer.data(if_ + 3);
    const auto *if__4 = buffer.data(if_ + 4);
    const auto *if__5 = buffer.data(if_ + 5);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__9 = buffer.data(if_ + 9);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__12 = buffer.data(if_ + 12);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__15 = buffer.data(if_ + 15);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__17 = buffer.data(if_ + 17);
    const auto *if__18 = buffer.data(if_ + 18);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__21 = buffer.data(if_ + 21);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__24 = buffer.data(if_ + 24);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__29 = buffer.data(if_ + 29);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__34 = buffer.data(if_ + 34);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__39 = buffer.data(if_ + 39);
    const auto *if__40 = buffer.data(if_ + 40);
    const auto *if__41 = buffer.data(if_ + 41);
    const auto *if__42 = buffer.data(if_ + 42);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__48 = buffer.data(if_ + 48);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__50 = buffer.data(if_ + 50);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__53 = buffer.data(if_ + 53);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__55 = buffer.data(if_ + 55);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__57 = buffer.data(if_ + 57);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__59 = buffer.data(if_ + 59);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__61 = buffer.data(if_ + 61);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__63 = buffer.data(if_ + 63);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__65 = buffer.data(if_ + 65);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__67 = buffer.data(if_ + 67);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__69 = buffer.data(if_ + 69);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__71 = buffer.data(if_ + 71);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__73 = buffer.data(if_ + 73);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__75 = buffer.data(if_ + 75);
    const auto *if__76 = buffer.data(if_ + 76);
    const auto *if__77 = buffer.data(if_ + 77);
    const auto *if__78 = buffer.data(if_ + 78);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__80 = buffer.data(if_ + 80);
    const auto *if__81 = buffer.data(if_ + 81);
    const auto *if__82 = buffer.data(if_ + 82);
    const auto *if__83 = buffer.data(if_ + 83);
    const auto *if__84 = buffer.data(if_ + 84);
    const auto *if__85 = buffer.data(if_ + 85);
    const auto *if__86 = buffer.data(if_ + 86);
    const auto *if__87 = buffer.data(if_ + 87);
    const auto *if__88 = buffer.data(if_ + 88);
    const auto *if__89 = buffer.data(if_ + 89);
    const auto *if__90 = buffer.data(if_ + 90);
    const auto *if__91 = buffer.data(if_ + 91);
    const auto *if__92 = buffer.data(if_ + 92);
    const auto *if__93 = buffer.data(if_ + 93);
    const auto *if__94 = buffer.data(if_ + 94);
    const auto *if__95 = buffer.data(if_ + 95);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__97 = buffer.data(if_ + 97);
    const auto *if__98 = buffer.data(if_ + 98);
    const auto *if__99 = buffer.data(if_ + 99);
    const auto *if__100 = buffer.data(if_ + 100);
    const auto *if__101 = buffer.data(if_ + 101);
    const auto *if__102 = buffer.data(if_ + 102);
    const auto *if__103 = buffer.data(if_ + 103);
    const auto *if__104 = buffer.data(if_ + 104);
    const auto *if__105 = buffer.data(if_ + 105);
    const auto *if__106 = buffer.data(if_ + 106);
    const auto *if__107 = buffer.data(if_ + 107);
    const auto *if__108 = buffer.data(if_ + 108);
    const auto *if__109 = buffer.data(if_ + 109);
    const auto *if__110 = buffer.data(if_ + 110);
    const auto *if__111 = buffer.data(if_ + 111);
    const auto *if__112 = buffer.data(if_ + 112);
    const auto *if__113 = buffer.data(if_ + 113);
    const auto *if__114 = buffer.data(if_ + 114);
    const auto *if__115 = buffer.data(if_ + 115);
    const auto *if__116 = buffer.data(if_ + 116);
    const auto *if__117 = buffer.data(if_ + 117);
    const auto *if__118 = buffer.data(if_ + 118);
    const auto *if__119 = buffer.data(if_ + 119);
    const auto *if__120 = buffer.data(if_ + 120);
    const auto *if__121 = buffer.data(if_ + 121);
    const auto *if__122 = buffer.data(if_ + 122);
    const auto *if__123 = buffer.data(if_ + 123);
    const auto *if__124 = buffer.data(if_ + 124);
    const auto *if__125 = buffer.data(if_ + 125);
    const auto *if__126 = buffer.data(if_ + 126);
    const auto *if__127 = buffer.data(if_ + 127);
    const auto *if__128 = buffer.data(if_ + 128);
    const auto *if__129 = buffer.data(if_ + 129);
    const auto *if__130 = buffer.data(if_ + 130);
    const auto *if__131 = buffer.data(if_ + 131);
    const auto *if__132 = buffer.data(if_ + 132);
    const auto *if__133 = buffer.data(if_ + 133);
    const auto *if__134 = buffer.data(if_ + 134);
    const auto *if__135 = buffer.data(if_ + 135);
    const auto *if__136 = buffer.data(if_ + 136);
    const auto *if__137 = buffer.data(if_ + 137);
    const auto *if__138 = buffer.data(if_ + 138);
    const auto *if__139 = buffer.data(if_ + 139);
    const auto *if__140 = buffer.data(if_ + 140);
    const auto *if__141 = buffer.data(if_ + 141);
    const auto *if__142 = buffer.data(if_ + 142);
    const auto *if__143 = buffer.data(if_ + 143);
    const auto *if__144 = buffer.data(if_ + 144);
    const auto *if__145 = buffer.data(if_ + 145);
    const auto *if__146 = buffer.data(if_ + 146);
    const auto *if__147 = buffer.data(if_ + 147);
    const auto *if__148 = buffer.data(if_ + 148);
    const auto *if__149 = buffer.data(if_ + 149);
    const auto *if__150 = buffer.data(if_ + 150);
    const auto *if__151 = buffer.data(if_ + 151);
    const auto *if__152 = buffer.data(if_ + 152);
    const auto *if__153 = buffer.data(if_ + 153);
    const auto *if__154 = buffer.data(if_ + 154);
    const auto *if__155 = buffer.data(if_ + 155);
    const auto *if__156 = buffer.data(if_ + 156);
    const auto *if__157 = buffer.data(if_ + 157);
    const auto *if__158 = buffer.data(if_ + 158);
    const auto *if__159 = buffer.data(if_ + 159);
    const auto *if__160 = buffer.data(if_ + 160);
    const auto *if__161 = buffer.data(if_ + 161);
    const auto *if__162 = buffer.data(if_ + 162);
    const auto *if__163 = buffer.data(if_ + 163);
    const auto *if__164 = buffer.data(if_ + 164);
    const auto *if__165 = buffer.data(if_ + 165);
    const auto *if__166 = buffer.data(if_ + 166);
    const auto *if__167 = buffer.data(if_ + 167);
    const auto *if__168 = buffer.data(if_ + 168);
    const auto *if__169 = buffer.data(if_ + 169);
    const auto *if__170 = buffer.data(if_ + 170);
    const auto *if__171 = buffer.data(if_ + 171);
    const auto *if__172 = buffer.data(if_ + 172);
    const auto *if__173 = buffer.data(if_ + 173);
    const auto *if__174 = buffer.data(if_ + 174);
    const auto *if__175 = buffer.data(if_ + 175);
    const auto *if__176 = buffer.data(if_ + 176);
    const auto *if__177 = buffer.data(if_ + 177);
    const auto *if__178 = buffer.data(if_ + 178);
    const auto *if__179 = buffer.data(if_ + 179);
    const auto *if__180 = buffer.data(if_ + 180);
    const auto *if__181 = buffer.data(if_ + 181);
    const auto *if__182 = buffer.data(if_ + 182);
    const auto *if__183 = buffer.data(if_ + 183);
    const auto *if__184 = buffer.data(if_ + 184);
    const auto *if__185 = buffer.data(if_ + 185);
    const auto *if__186 = buffer.data(if_ + 186);
    const auto *if__187 = buffer.data(if_ + 187);
    const auto *if__188 = buffer.data(if_ + 188);
    const auto *if__189 = buffer.data(if_ + 189);
    const auto *if__190 = buffer.data(if_ + 190);
    const auto *if__191 = buffer.data(if_ + 191);
    const auto *if__192 = buffer.data(if_ + 192);
    const auto *if__193 = buffer.data(if_ + 193);
    const auto *if__194 = buffer.data(if_ + 194);
    const auto *if__195 = buffer.data(if_ + 195);
    const auto *if__196 = buffer.data(if_ + 196);
    const auto *if__197 = buffer.data(if_ + 197);
    const auto *if__198 = buffer.data(if_ + 198);
    const auto *if__199 = buffer.data(if_ + 199);
    const auto *if__200 = buffer.data(if_ + 200);
    const auto *if__201 = buffer.data(if_ + 201);
    const auto *if__202 = buffer.data(if_ + 202);
    const auto *if__203 = buffer.data(if_ + 203);
    const auto *if__204 = buffer.data(if_ + 204);
    const auto *if__205 = buffer.data(if_ + 205);
    const auto *if__206 = buffer.data(if_ + 206);
    const auto *if__207 = buffer.data(if_ + 207);
    const auto *if__208 = buffer.data(if_ + 208);
    const auto *if__209 = buffer.data(if_ + 209);
    const auto *if__216 = buffer.data(if_ + 216);
    const auto *if__217 = buffer.data(if_ + 217);
    const auto *if__218 = buffer.data(if_ + 218);
    const auto *if__219 = buffer.data(if_ + 219);
    const auto *if__226 = buffer.data(if_ + 226);
    const auto *if__227 = buffer.data(if_ + 227);
    const auto *if__228 = buffer.data(if_ + 228);
    const auto *if__229 = buffer.data(if_ + 229);
    const auto *if__236 = buffer.data(if_ + 236);
    const auto *if__237 = buffer.data(if_ + 237);
    const auto *if__238 = buffer.data(if_ + 238);
    const auto *if__239 = buffer.data(if_ + 239);
    const auto *if__246 = buffer.data(if_ + 246);
    const auto *if__247 = buffer.data(if_ + 247);
    const auto *if__248 = buffer.data(if_ + 248);
    const auto *if__249 = buffer.data(if_ + 249);
    const auto *if__256 = buffer.data(if_ + 256);
    const auto *if__257 = buffer.data(if_ + 257);
    const auto *if__258 = buffer.data(if_ + 258);
    const auto *if__259 = buffer.data(if_ + 259);
    const auto *if__266 = buffer.data(if_ + 266);
    const auto *if__267 = buffer.data(if_ + 267);
    const auto *if__268 = buffer.data(if_ + 268);
    const auto *if__269 = buffer.data(if_ + 269);
    const auto *if__279 = buffer.data(if_ + 279);

#pragma omp simd aligned(ab_x, hf_11, hf_16, hf_61, hf_66, hf_151, hf_156, if__11, if__16, \
                         if__61, if__66, if__151, if__156 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * ab_x[k] * hf_11[k]
                 - f_0 * ab_x[k] * hf_16[k]
                 - f_1 * ab_x[k] * hf_61[k]
                 + f_1 * ab_x[k] * hf_66[k]
                 + f_2 * ab_x[k] * hf_151[k]
                 - f_2 * ab_x[k] * hf_156[k]
                 + f_0 * if__11[k]
                 - f_0 * if__16[k]
                 - f_1 * if__61[k]
                 + f_1 * if__66[k]
                 + f_2 * if__151[k]
                 - f_2 * if__156[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_14, hf_17, hf_64, hf_67, hf_154, hf_157, if__14, \
                         if__37, if__64, if__107, if__154, if__217 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = f_3 * ab_x[k] * hf_14[k]
                 - f_4 * ab_y[k] * hf_17[k]
                 - f_5 * ab_x[k] * hf_64[k]
                 + f_6 * ab_y[k] * hf_67[k]
                 + f_7 * ab_x[k] * hf_154[k]
                 - f_8 * ab_y[k] * hf_157[k]
                 + f_3 * if__14[k]
                 - f_4 * if__37[k]
                 - f_5 * if__64[k]
                 + f_6 * if__107[k]
                 + f_7 * if__154[k]
                 - f_8 * if__217[k];
    }

#pragma omp simd aligned(ab_x, hf_11, hf_16, hf_18, hf_61, hf_66, hf_68, hf_151, hf_156, \
                         hf_158, if__11, if__16, if__18, if__61, if__66, if__68, if__151, \
                         if__156, if__158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_9 * ab_x[k] * hf_11[k]
                 - f_9 * ab_x[k] * hf_16[k]
                 + f_10 * ab_x[k] * hf_18[k]
                 + f_11 * ab_x[k] * hf_61[k]
                 + f_11 * ab_x[k] * hf_66[k]
                 - f_12 * ab_x[k] * hf_68[k]
                 - f_13 * ab_x[k] * hf_151[k]
                 - f_13 * ab_x[k] * hf_156[k]
                 + f_14 * ab_x[k] * hf_158[k]
                 - f_9 * if__11[k]
                 - f_9 * if__16[k]
                 + f_10 * if__18[k]
                 + f_11 * if__61[k]
                 + f_11 * if__66[k]
                 - f_12 * if__68[k]
                 - f_13 * if__151[k]
                 - f_13 * if__156[k]
                 + f_14 * if__158[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_14, hf_17, hf_19, hf_64, hf_67, hf_69, hf_154, hf_157, \
                         hf_159, if__14, if__37, if__39, if__64, if__107, if__109, if__154, \
                         if__217, if__219 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_15 * ab_x[k] * hf_14[k]
                 - f_15 * ab_y[k] * hf_17[k]
                 + f_16 * ab_y[k] * hf_19[k]
                 + f_17 * ab_x[k] * hf_64[k]
                 + f_17 * ab_y[k] * hf_67[k]
                 - f_18 * ab_y[k] * hf_69[k]
                 - f_19 * ab_x[k] * hf_154[k]
                 - f_19 * ab_y[k] * hf_157[k]
                 + f_20 * ab_y[k] * hf_159[k]
                 - f_15 * if__14[k]
                 - f_15 * if__37[k]
                 + f_16 * if__39[k]
                 + f_17 * if__64[k]
                 + f_17 * if__107[k]
                 - f_18 * if__109[k]
                 - f_19 * if__154[k]
                 - f_19 * if__217[k]
                 + f_20 * if__219[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hf_10, hf_13, hf_15, hf_16, hf_18, hf_19, hf_60, \
                         hf_63, hf_65, hf_66, hf_68, hf_69, hf_150, hf_153, hf_155, hf_156, \
                         hf_158, hf_159, if__10, if__13, if__15, if__36, if__38, if__49, \
                         if__60, if__63, if__65, if__106, if__108, if__119, if__150, if__153, \
                         if__155, if__216, if__218, if__229 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_21 * ab_x[k] * hf_10[k]
                 + f_22 * ab_x[k] * hf_13[k]
                 - f_23 * ab_x[k] * hf_15[k]
                 + f_21 * ab_y[k] * hf_16[k]
                 - f_23 * ab_y[k] * hf_18[k]
                 + f_24 * ab_z[k] * hf_19[k]
                 - f_22 * ab_x[k] * hf_60[k]
                 - f_25 * ab_x[k] * hf_63[k]
                 + f_26 * ab_x[k] * hf_65[k]
                 - f_22 * ab_y[k] * hf_66[k]
                 + f_26 * ab_y[k] * hf_68[k]
                 - f_27 * ab_z[k] * hf_69[k]
                 + f_28 * ab_x[k] * hf_150[k]
                 + f_29 * ab_x[k] * hf_153[k]
                 - f_30 * ab_x[k] * hf_155[k]
                 + f_28 * ab_y[k] * hf_156[k]
                 - f_30 * ab_y[k] * hf_158[k]
                 + f_31 * ab_z[k] * hf_159[k]
                 + f_21 * if__10[k]
                 + f_22 * if__13[k]
                 - f_23 * if__15[k]
                 + f_21 * if__36[k]
                 - f_23 * if__38[k]
                 + f_24 * if__49[k]
                 - f_22 * if__60[k]
                 - f_25 * if__63[k]
                 + f_26 * if__65[k]
                 - f_22 * if__106[k]
                 + f_26 * if__108[k]
                 - f_27 * if__119[k]
                 + f_28 * if__150[k]
                 + f_29 * if__153[k]
                 - f_30 * if__155[k]
                 + f_28 * if__216[k]
                 - f_30 * if__218[k]
                 + f_31 * if__229[k];
    }

#pragma omp simd aligned(ab_x, hf_12, hf_17, hf_19, hf_62, hf_67, hf_69, hf_152, hf_157, \
                         hf_159, if__12, if__17, if__19, if__62, if__67, if__69, if__152, \
                         if__157, if__159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -f_15 * ab_x[k] * hf_12[k]
                 - f_15 * ab_x[k] * hf_17[k]
                 + f_16 * ab_x[k] * hf_19[k]
                 + f_17 * ab_x[k] * hf_62[k]
                 + f_17 * ab_x[k] * hf_67[k]
                 - f_18 * ab_x[k] * hf_69[k]
                 - f_19 * ab_x[k] * hf_152[k]
                 - f_19 * ab_x[k] * hf_157[k]
                 + f_20 * ab_x[k] * hf_159[k]
                 - f_15 * if__12[k]
                 - f_15 * if__17[k]
                 + f_16 * if__19[k]
                 + f_17 * if__62[k]
                 + f_17 * if__67[k]
                 - f_18 * if__69[k]
                 - f_19 * if__152[k]
                 - f_19 * if__157[k]
                 + f_20 * if__159[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_10, hf_15, hf_16, hf_18, hf_60, hf_65, hf_66, hf_68, \
                         hf_150, hf_155, hf_156, hf_158, if__10, if__15, if__36, if__38, \
                         if__60, if__65, if__106, if__108, if__150, if__155, if__216, \
                         if__218 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_32 * ab_x[k] * hf_10[k]
                 + f_33 * ab_x[k] * hf_15[k]
                 + f_32 * ab_y[k] * hf_16[k]
                 - f_33 * ab_y[k] * hf_18[k]
                 + f_9 * ab_x[k] * hf_60[k]
                 - f_10 * ab_x[k] * hf_65[k]
                 - f_9 * ab_y[k] * hf_66[k]
                 + f_10 * ab_y[k] * hf_68[k]
                 - f_34 * ab_x[k] * hf_150[k]
                 + f_35 * ab_x[k] * hf_155[k]
                 + f_34 * ab_y[k] * hf_156[k]
                 - f_35 * ab_y[k] * hf_158[k]
                 - f_32 * if__10[k]
                 + f_33 * if__15[k]
                 + f_32 * if__36[k]
                 - f_33 * if__38[k]
                 + f_9 * if__60[k]
                 - f_10 * if__65[k]
                 - f_9 * if__106[k]
                 + f_10 * if__108[k]
                 - f_34 * if__150[k]
                 + f_35 * if__155[k]
                 + f_34 * if__216[k]
                 - f_35 * if__218[k];
    }

#pragma omp simd aligned(ab_x, hf_12, hf_17, hf_62, hf_67, hf_152, hf_157, if__12, if__17, \
                         if__62, if__67, if__152, if__157 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_4 * ab_x[k] * hf_12[k]
                 - f_3 * ab_x[k] * hf_17[k]
                 - f_6 * ab_x[k] * hf_62[k]
                 + f_5 * ab_x[k] * hf_67[k]
                 + f_8 * ab_x[k] * hf_152[k]
                 - f_7 * ab_x[k] * hf_157[k]
                 + f_4 * if__12[k]
                 - f_3 * if__17[k]
                 - f_6 * if__62[k]
                 + f_5 * if__67[k]
                 + f_8 * if__152[k]
                 - f_7 * if__157[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_10, hf_13, hf_16, hf_60, hf_63, hf_66, hf_150, hf_153, \
                         hf_156, if__10, if__13, if__36, if__60, if__63, if__106, if__150, \
                         if__153, if__216 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_36 * ab_x[k] * hf_10[k]
                 - f_37 * ab_x[k] * hf_13[k]
                 + f_36 * ab_y[k] * hf_16[k]
                 - f_38 * ab_x[k] * hf_60[k]
                 + f_39 * ab_x[k] * hf_63[k]
                 - f_38 * ab_y[k] * hf_66[k]
                 + f_40 * ab_x[k] * hf_150[k]
                 - f_41 * ab_x[k] * hf_153[k]
                 + f_40 * ab_y[k] * hf_156[k]
                 + f_36 * if__10[k]
                 - f_37 * if__13[k]
                 + f_36 * if__36[k]
                 - f_38 * if__60[k]
                 + f_39 * if__63[k]
                 - f_38 * if__106[k]
                 + f_40 * if__150[k]
                 - f_41 * if__153[k]
                 + f_40 * if__216[k];
    }

#pragma omp simd aligned(ab_x, hf_41, hf_46, hf_111, hf_116, if__41, if__46, if__111, \
                         if__116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = 26.25 * ab_x[k] * hf_41[k]
                 - 26.25 * ab_x[k] * hf_46[k]
                 - 26.25 * ab_x[k] * hf_111[k]
                 + 26.25 * ab_x[k] * hf_116[k]
                 + 26.25 * if__41[k]
                 - 26.25 * if__46[k]
                 - 26.25 * if__111[k]
                 + 26.25 * if__116[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_44, hf_47, hf_114, hf_117, if__44, if__77, if__114, \
                         if__167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = f_42 * ab_x[k] * hf_44[k]
                  - f_43 * ab_y[k] * hf_47[k]
                  - f_42 * ab_x[k] * hf_114[k]
                  + f_43 * ab_y[k] * hf_117[k]
                  + f_42 * if__44[k]
                  - f_43 * if__77[k]
                  - f_42 * if__114[k]
                  + f_43 * if__167[k];
    }

#pragma omp simd aligned(ab_x, hf_41, hf_46, hf_48, hf_111, hf_116, hf_118, if__41, if__46, \
                         if__48, if__111, if__116, if__118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_44 * ab_x[k] * hf_41[k]
                  - f_44 * ab_x[k] * hf_46[k]
                  + f_45 * ab_x[k] * hf_48[k]
                  + f_44 * ab_x[k] * hf_111[k]
                  + f_44 * ab_x[k] * hf_116[k]
                  - f_45 * ab_x[k] * hf_118[k]
                  - f_44 * if__41[k]
                  - f_44 * if__46[k]
                  + f_45 * if__48[k]
                  + f_44 * if__111[k]
                  + f_44 * if__116[k]
                  - f_45 * if__118[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_44, hf_47, hf_49, hf_114, hf_117, hf_119, if__44, \
                         if__77, if__79, if__114, if__167, if__169 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = -f_26 * ab_x[k] * hf_44[k]
                  - f_26 * ab_y[k] * hf_47[k]
                  + f_46 * ab_y[k] * hf_49[k]
                  + f_26 * ab_x[k] * hf_114[k]
                  + f_26 * ab_y[k] * hf_117[k]
                  - f_46 * ab_y[k] * hf_119[k]
                  - f_26 * if__44[k]
                  - f_26 * if__77[k]
                  + f_46 * if__79[k]
                  + f_26 * if__114[k]
                  + f_26 * if__167[k]
                  - f_46 * if__169[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hf_40, hf_43, hf_45, hf_46, hf_48, hf_49, hf_110, \
                         hf_113, hf_115, hf_116, hf_118, hf_119, if__40, if__43, if__45, \
                         if__76, if__78, if__89, if__110, if__113, if__115, if__166, if__168, \
                         if__179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_47 * ab_x[k] * hf_40[k]
                  + f_48 * ab_x[k] * hf_43[k]
                  - f_49 * ab_x[k] * hf_45[k]
                  + f_47 * ab_y[k] * hf_46[k]
                  - f_49 * ab_y[k] * hf_48[k]
                  + f_50 * ab_z[k] * hf_49[k]
                  - f_47 * ab_x[k] * hf_110[k]
                  - f_48 * ab_x[k] * hf_113[k]
                  + f_49 * ab_x[k] * hf_115[k]
                  - f_47 * ab_y[k] * hf_116[k]
                  + f_49 * ab_y[k] * hf_118[k]
                  - f_50 * ab_z[k] * hf_119[k]
                  + f_47 * if__40[k]
                  + f_48 * if__43[k]
                  - f_49 * if__45[k]
                  + f_47 * if__76[k]
                  - f_49 * if__78[k]
                  + f_50 * if__89[k]
                  - f_47 * if__110[k]
                  - f_48 * if__113[k]
                  + f_49 * if__115[k]
                  - f_47 * if__166[k]
                  + f_49 * if__168[k]
                  - f_50 * if__179[k];
    }

#pragma omp simd aligned(ab_x, hf_42, hf_47, hf_49, hf_112, hf_117, hf_119, if__42, if__47, \
                         if__49, if__112, if__117, if__119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_26 * ab_x[k] * hf_42[k]
                  - f_26 * ab_x[k] * hf_47[k]
                  + f_46 * ab_x[k] * hf_49[k]
                  + f_26 * ab_x[k] * hf_112[k]
                  + f_26 * ab_x[k] * hf_117[k]
                  - f_46 * ab_x[k] * hf_119[k]
                  - f_26 * if__42[k]
                  - f_26 * if__47[k]
                  + f_46 * if__49[k]
                  + f_26 * if__112[k]
                  + f_26 * if__117[k]
                  - f_46 * if__119[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_40, hf_45, hf_46, hf_48, hf_110, hf_115, hf_116, \
                         hf_118, if__40, if__45, if__76, if__78, if__110, if__115, if__166, \
                         if__168 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_51 * ab_x[k] * hf_40[k]
                  + f_52 * ab_x[k] * hf_45[k]
                  + f_51 * ab_y[k] * hf_46[k]
                  - f_52 * ab_y[k] * hf_48[k]
                  + f_51 * ab_x[k] * hf_110[k]
                  - f_52 * ab_x[k] * hf_115[k]
                  - f_51 * ab_y[k] * hf_116[k]
                  + f_52 * ab_y[k] * hf_118[k]
                  - f_51 * if__40[k]
                  + f_52 * if__45[k]
                  + f_51 * if__76[k]
                  - f_52 * if__78[k]
                  + f_51 * if__110[k]
                  - f_52 * if__115[k]
                  - f_51 * if__166[k]
                  + f_52 * if__168[k];
    }

#pragma omp simd aligned(ab_x, hf_42, hf_47, hf_112, hf_117, if__42, if__47, if__112, \
                         if__117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_43 * ab_x[k] * hf_42[k]
                  - f_42 * ab_x[k] * hf_47[k]
                  - f_43 * ab_x[k] * hf_112[k]
                  + f_42 * ab_x[k] * hf_117[k]
                  + f_43 * if__42[k]
                  - f_42 * if__47[k]
                  - f_43 * if__112[k]
                  + f_42 * if__117[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_40, hf_43, hf_46, hf_110, hf_113, hf_116, if__40, \
                         if__43, if__76, if__110, if__113, if__166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = 6.5625 * ab_x[k] * hf_40[k]
                  - 39.375 * ab_x[k] * hf_43[k]
                  + 6.5625 * ab_y[k] * hf_46[k]
                  - 6.5625 * ab_x[k] * hf_110[k]
                  + 39.375 * ab_x[k] * hf_113[k]
                  - 6.5625 * ab_y[k] * hf_116[k]
                  + 6.5625 * if__40[k]
                  - 39.375 * if__43[k]
                  + 6.5625 * if__76[k]
                  - 6.5625 * if__110[k]
                  + 39.375 * if__113[k]
                  - 6.5625 * if__166[k];
    }

#pragma omp simd aligned(ab_x, hf_11, hf_16, hf_61, hf_66, hf_81, hf_86, hf_151, hf_156, \
                         hf_171, hf_176, if__11, if__16, if__61, if__66, if__81, if__86, \
                         if__151, if__156, if__171, if__176 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_53 * ab_x[k] * hf_11[k]
                  + f_53 * ab_x[k] * hf_16[k]
                  - f_54 * ab_x[k] * hf_61[k]
                  + f_54 * ab_x[k] * hf_66[k]
                  + f_55 * ab_x[k] * hf_81[k]
                  - f_55 * ab_x[k] * hf_86[k]
                  + f_56 * ab_x[k] * hf_151[k]
                  - f_56 * ab_x[k] * hf_156[k]
                  - f_57 * ab_x[k] * hf_171[k]
                  + f_57 * ab_x[k] * hf_176[k]
                  - f_53 * if__11[k]
                  + f_53 * if__16[k]
                  - f_54 * if__61[k]
                  + f_54 * if__66[k]
                  + f_55 * if__81[k]
                  - f_55 * if__86[k]
                  + f_56 * if__151[k]
                  - f_56 * if__156[k]
                  - f_57 * if__171[k]
                  + f_57 * if__176[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_14, hf_17, hf_64, hf_67, hf_84, hf_87, hf_154, hf_157, \
                         hf_174, hf_177, if__14, if__37, if__64, if__84, if__107, if__127, \
                         if__154, if__174, if__217, if__237 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -9.84375 * ab_x[k] * hf_14[k]
                  + 3.28125 * ab_y[k] * hf_17[k]
                  - 6.5625 * ab_x[k] * hf_64[k]
                  + 2.1875 * ab_y[k] * hf_67[k]
                  + 78.75 * ab_x[k] * hf_84[k]
                  - 26.25 * ab_y[k] * hf_87[k]
                  + 3.28125 * ab_x[k] * hf_154[k]
                  - 1.09375 * ab_y[k] * hf_157[k]
                  - 26.25 * ab_x[k] * hf_174[k]
                  + 8.75 * ab_y[k] * hf_177[k]
                  - 9.84375 * if__14[k]
                  + 3.28125 * if__37[k]
                  - 6.5625 * if__64[k]
                  + 78.75 * if__84[k]
                  + 2.1875 * if__107[k]
                  - 26.25 * if__127[k]
                  + 3.28125 * if__154[k]
                  - 26.25 * if__174[k]
                  - 1.09375 * if__217[k]
                  + 8.75 * if__237[k];
    }

#pragma omp simd aligned(ab_x, hf_11, hf_16, hf_18, hf_61, hf_66, hf_68, hf_81, hf_86, hf_88, \
                         hf_151, hf_156, hf_158, hf_171, hf_176, hf_178, if__11, if__16, \
                         if__18, if__61, if__66, if__68, if__81, if__86, if__88, if__151, \
                         if__156, if__158, if__171, if__176, if__178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_58 * ab_x[k] * hf_11[k]
                  + f_58 * ab_x[k] * hf_16[k]
                  - f_23 * ab_x[k] * hf_18[k]
                  + f_59 * ab_x[k] * hf_61[k]
                  + f_59 * ab_x[k] * hf_66[k]
                  - f_27 * ab_x[k] * hf_68[k]
                  - f_60 * ab_x[k] * hf_81[k]
                  - f_60 * ab_x[k] * hf_86[k]
                  + f_61 * ab_x[k] * hf_88[k]
                  - f_62 * ab_x[k] * hf_151[k]
                  - f_62 * ab_x[k] * hf_156[k]
                  + f_24 * ab_x[k] * hf_158[k]
                  + f_63 * ab_x[k] * hf_171[k]
                  + f_63 * ab_x[k] * hf_176[k]
                  - f_46 * ab_x[k] * hf_178[k]
                  + f_58 * if__11[k]
                  + f_58 * if__16[k]
                  - f_23 * if__18[k]
                  + f_59 * if__61[k]
                  + f_59 * if__66[k]
                  - f_27 * if__68[k]
                  - f_60 * if__81[k]
                  - f_60 * if__86[k]
                  + f_61 * if__88[k]
                  - f_62 * if__151[k]
                  - f_62 * if__156[k]
                  + f_24 * if__158[k]
                  + f_63 * if__171[k]
                  + f_63 * if__176[k]
                  - f_46 * if__178[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_14, hf_17, hf_19, hf_64, hf_67, hf_69, hf_84, hf_87, \
                         hf_89, hf_154, hf_157, hf_159, hf_174, hf_177, hf_179, if__14, \
                         if__37, if__39, if__64, if__84, if__107, if__109, if__127, if__129, \
                         if__154, if__174, if__217, if__219, if__237, \
                         if__239 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_64 * ab_x[k] * hf_14[k]
                  + f_64 * ab_y[k] * hf_17[k]
                  - f_51 * ab_y[k] * hf_19[k]
                  + f_65 * ab_x[k] * hf_64[k]
                  + f_65 * ab_y[k] * hf_67[k]
                  - f_66 * ab_y[k] * hf_69[k]
                  - f_52 * ab_x[k] * hf_84[k]
                  - f_52 * ab_y[k] * hf_87[k]
                  + f_67 * ab_y[k] * hf_89[k]
                  - f_68 * ab_x[k] * hf_154[k]
                  - f_68 * ab_y[k] * hf_157[k]
                  + f_69 * ab_y[k] * hf_159[k]
                  + f_44 * ab_x[k] * hf_174[k]
                  + f_44 * ab_y[k] * hf_177[k]
                  - f_70 * ab_y[k] * hf_179[k]
                  + f_64 * if__14[k]
                  + f_64 * if__37[k]
                  - f_51 * if__39[k]
                  + f_65 * if__64[k]
                  - f_52 * if__84[k]
                  + f_65 * if__107[k]
                  - f_66 * if__109[k]
                  - f_52 * if__127[k]
                  + f_67 * if__129[k]
                  - f_68 * if__154[k]
                  + f_44 * if__174[k]
                  - f_68 * if__217[k]
                  + f_69 * if__219[k]
                  + f_44 * if__237[k]
                  - f_70 * if__239[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hf_10, hf_13, hf_15, hf_16, hf_18, hf_19, hf_60, \
                         hf_63, hf_65, hf_66, hf_68, hf_69, hf_80, hf_83, hf_85, hf_86, hf_88, \
                         hf_89, hf_150, hf_153, hf_155, hf_156, hf_158, hf_159, hf_170, \
                         hf_173, hf_175, hf_176, hf_178, hf_179, if__10, if__13, if__15, \
                         if__36, if__38, if__49, if__60, if__63, if__65, if__80, if__83, \
                         if__85, if__106, if__108, if__119, if__126, if__128, if__139, \
                         if__150, if__153, if__155, if__170, if__173, if__175, if__216, \
                         if__218, if__229, if__236, if__238, if__249 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_71 * ab_x[k] * hf_10[k]
                  - f_72 * ab_x[k] * hf_13[k]
                  + f_14 * ab_x[k] * hf_15[k]
                  - f_71 * ab_y[k] * hf_16[k]
                  + f_14 * ab_y[k] * hf_18[k]
                  - f_73 * ab_z[k] * hf_19[k]
                  - f_34 * ab_x[k] * hf_60[k]
                  - f_13 * ab_x[k] * hf_63[k]
                  + f_74 * ab_x[k] * hf_65[k]
                  - f_34 * ab_y[k] * hf_66[k]
                  + f_74 * ab_y[k] * hf_68[k]
                  - f_75 * ab_z[k] * hf_69[k]
                  + f_14 * ab_x[k] * hf_80[k]
                  + f_76 * ab_x[k] * hf_83[k]
                  - f_77 * ab_x[k] * hf_85[k]
                  + f_14 * ab_y[k] * hf_86[k]
                  - f_77 * ab_y[k] * hf_88[k]
                  + f_78 * ab_z[k] * hf_89[k]
                  + f_79 * ab_x[k] * hf_150[k]
                  + f_34 * ab_x[k] * hf_153[k]
                  - f_73 * ab_x[k] * hf_155[k]
                  + f_79 * ab_y[k] * hf_156[k]
                  - f_73 * ab_y[k] * hf_158[k]
                  + f_80 * ab_z[k] * hf_159[k]
                  - f_73 * ab_x[k] * hf_170[k]
                  - f_74 * ab_x[k] * hf_173[k]
                  + f_78 * ab_x[k] * hf_175[k]
                  - f_73 * ab_y[k] * hf_176[k]
                  + f_78 * ab_y[k] * hf_178[k]
                  - f_81 * ab_z[k] * hf_179[k]
                  - f_71 * if__10[k]
                  - f_72 * if__13[k]
                  + f_14 * if__15[k]
                  - f_71 * if__36[k]
                  + f_14 * if__38[k]
                  - f_73 * if__49[k]
                  - f_34 * if__60[k]
                  - f_13 * if__63[k]
                  + f_74 * if__65[k]
                  + f_14 * if__80[k]
                  + f_76 * if__83[k]
                  - f_77 * if__85[k]
                  - f_34 * if__106[k]
                  + f_74 * if__108[k]
                  - f_75 * if__119[k]
                  + f_14 * if__126[k]
                  - f_77 * if__128[k]
                  + f_78 * if__139[k]
                  + f_79 * if__150[k]
                  + f_34 * if__153[k]
                  - f_73 * if__155[k]
                  - f_73 * if__170[k]
                  - f_74 * if__173[k]
                  + f_78 * if__175[k]
                  + f_79 * if__216[k]
                  - f_73 * if__218[k]
                  + f_80 * if__229[k]
                  - f_73 * if__236[k]
                  + f_78 * if__238[k]
                  - f_81 * if__249[k];
    }

#pragma omp simd aligned(ab_x, hf_12, hf_17, hf_19, hf_62, hf_67, hf_69, hf_82, hf_87, hf_89, \
                         hf_152, hf_157, hf_159, hf_172, hf_177, hf_179, if__12, if__17, \
                         if__19, if__62, if__67, if__69, if__82, if__87, if__89, if__152, \
                         if__157, if__159, if__172, if__177, if__179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_64 * ab_x[k] * hf_12[k]
                  + f_64 * ab_x[k] * hf_17[k]
                  - f_51 * ab_x[k] * hf_19[k]
                  + f_65 * ab_x[k] * hf_62[k]
                  + f_65 * ab_x[k] * hf_67[k]
                  - f_66 * ab_x[k] * hf_69[k]
                  - f_52 * ab_x[k] * hf_82[k]
                  - f_52 * ab_x[k] * hf_87[k]
                  + f_67 * ab_x[k] * hf_89[k]
                  - f_68 * ab_x[k] * hf_152[k]
                  - f_68 * ab_x[k] * hf_157[k]
                  + f_69 * ab_x[k] * hf_159[k]
                  + f_44 * ab_x[k] * hf_172[k]
                  + f_44 * ab_x[k] * hf_177[k]
                  - f_70 * ab_x[k] * hf_179[k]
                  + f_64 * if__12[k]
                  + f_64 * if__17[k]
                  - f_51 * if__19[k]
                  + f_65 * if__62[k]
                  + f_65 * if__67[k]
                  - f_66 * if__69[k]
                  - f_52 * if__82[k]
                  - f_52 * if__87[k]
                  + f_67 * if__89[k]
                  - f_68 * if__152[k]
                  - f_68 * if__157[k]
                  + f_69 * if__159[k]
                  + f_44 * if__172[k]
                  + f_44 * if__177[k]
                  - f_70 * if__179[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_10, hf_15, hf_16, hf_18, hf_60, hf_65, hf_66, hf_68, \
                         hf_80, hf_85, hf_86, hf_88, hf_150, hf_155, hf_156, hf_158, hf_170, \
                         hf_175, hf_176, hf_178, if__10, if__15, if__36, if__38, if__60, \
                         if__65, if__80, if__85, if__106, if__108, if__126, if__128, if__150, \
                         if__155, if__170, if__175, if__216, if__218, if__236, \
                         if__238 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_82 * ab_x[k] * hf_10[k]
                  - f_25 * ab_x[k] * hf_15[k]
                  - f_82 * ab_y[k] * hf_16[k]
                  + f_25 * ab_y[k] * hf_18[k]
                  + f_62 * ab_x[k] * hf_60[k]
                  - f_24 * ab_x[k] * hf_65[k]
                  - f_62 * ab_y[k] * hf_66[k]
                  + f_24 * ab_y[k] * hf_68[k]
                  - f_27 * ab_x[k] * hf_80[k]
                  + f_83 * ab_x[k] * hf_85[k]
                  + f_27 * ab_y[k] * hf_86[k]
                  - f_83 * ab_y[k] * hf_88[k]
                  - f_84 * ab_x[k] * hf_150[k]
                  + f_58 * ab_x[k] * hf_155[k]
                  + f_84 * ab_y[k] * hf_156[k]
                  - f_58 * ab_y[k] * hf_158[k]
                  + f_85 * ab_x[k] * hf_170[k]
                  - f_60 * ab_x[k] * hf_175[k]
                  - f_85 * ab_y[k] * hf_176[k]
                  + f_60 * ab_y[k] * hf_178[k]
                  + f_82 * if__10[k]
                  - f_25 * if__15[k]
                  - f_82 * if__36[k]
                  + f_25 * if__38[k]
                  + f_62 * if__60[k]
                  - f_24 * if__65[k]
                  - f_27 * if__80[k]
                  + f_83 * if__85[k]
                  - f_62 * if__106[k]
                  + f_24 * if__108[k]
                  + f_27 * if__126[k]
                  - f_83 * if__128[k]
                  - f_84 * if__150[k]
                  + f_58 * if__155[k]
                  + f_85 * if__170[k]
                  - f_60 * if__175[k]
                  + f_84 * if__216[k]
                  - f_58 * if__218[k]
                  - f_85 * if__236[k]
                  + f_60 * if__238[k];
    }

#pragma omp simd aligned(ab_x, hf_12, hf_17, hf_62, hf_67, hf_82, hf_87, hf_152, hf_157, \
                         hf_172, hf_177, if__12, if__17, if__62, if__67, if__82, if__87, \
                         if__152, if__157, if__172, if__177 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -3.28125 * ab_x[k] * hf_12[k]
                  + 9.84375 * ab_x[k] * hf_17[k]
                  - 2.1875 * ab_x[k] * hf_62[k]
                  + 6.5625 * ab_x[k] * hf_67[k]
                  + 26.25 * ab_x[k] * hf_82[k]
                  - 78.75 * ab_x[k] * hf_87[k]
                  + 1.09375 * ab_x[k] * hf_152[k]
                  - 3.28125 * ab_x[k] * hf_157[k]
                  - 8.75 * ab_x[k] * hf_172[k]
                  + 26.25 * ab_x[k] * hf_177[k]
                  - 3.28125 * if__12[k]
                  + 9.84375 * if__17[k]
                  - 2.1875 * if__62[k]
                  + 6.5625 * if__67[k]
                  + 26.25 * if__82[k]
                  - 78.75 * if__87[k]
                  + 1.09375 * if__152[k]
                  - 3.28125 * if__157[k]
                  - 8.75 * if__172[k]
                  + 26.25 * if__177[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_10, hf_13, hf_16, hf_60, hf_63, hf_66, hf_80, hf_83, \
                         hf_86, hf_150, hf_153, hf_156, hf_170, hf_173, hf_176, if__10, \
                         if__13, if__36, if__60, if__63, if__80, if__83, if__106, if__126, \
                         if__150, if__153, if__170, if__173, if__216, \
                         if__236 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_86 * ab_x[k] * hf_10[k]
                  + f_87 * ab_x[k] * hf_13[k]
                  - f_86 * ab_y[k] * hf_16[k]
                  - f_88 * ab_x[k] * hf_60[k]
                  + f_53 * ab_x[k] * hf_63[k]
                  - f_88 * ab_y[k] * hf_66[k]
                  + f_89 * ab_x[k] * hf_80[k]
                  - f_42 * ab_x[k] * hf_83[k]
                  + f_89 * ab_y[k] * hf_86[k]
                  + f_90 * ab_x[k] * hf_150[k]
                  - f_91 * ab_x[k] * hf_153[k]
                  + f_90 * ab_y[k] * hf_156[k]
                  - f_54 * ab_x[k] * hf_170[k]
                  + f_43 * ab_x[k] * hf_173[k]
                  - f_54 * ab_y[k] * hf_176[k]
                  - f_86 * if__10[k]
                  + f_87 * if__13[k]
                  - f_86 * if__36[k]
                  - f_88 * if__60[k]
                  + f_53 * if__63[k]
                  + f_89 * if__80[k]
                  - f_42 * if__83[k]
                  - f_88 * if__106[k]
                  + f_89 * if__126[k]
                  + f_90 * if__150[k]
                  - f_91 * if__153[k]
                  - f_54 * if__170[k]
                  + f_43 * if__173[k]
                  + f_90 * if__216[k]
                  - f_54 * if__236[k];
    }

#pragma omp simd aligned(ab_x, hf_41, hf_46, hf_111, hf_116, hf_131, hf_136, if__41, if__46, \
                         if__111, if__116, if__131, if__136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_92 * ab_x[k] * hf_41[k]
                  + f_92 * ab_x[k] * hf_46[k]
                  - f_92 * ab_x[k] * hf_111[k]
                  + f_92 * ab_x[k] * hf_116[k]
                  + f_93 * ab_x[k] * hf_131[k]
                  - f_93 * ab_x[k] * hf_136[k]
                  - f_92 * if__41[k]
                  + f_92 * if__46[k]
                  - f_92 * if__111[k]
                  + f_92 * if__116[k]
                  + f_93 * if__131[k]
                  - f_93 * if__136[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_44, hf_47, hf_114, hf_117, hf_134, hf_137, if__44, \
                         if__77, if__114, if__134, if__167, if__187 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_94 * ab_x[k] * hf_44[k]
                  + f_95 * ab_y[k] * hf_47[k]
                  - f_94 * ab_x[k] * hf_114[k]
                  + f_95 * ab_y[k] * hf_117[k]
                  + f_96 * ab_x[k] * hf_134[k]
                  - f_97 * ab_y[k] * hf_137[k]
                  - f_94 * if__44[k]
                  + f_95 * if__77[k]
                  - f_94 * if__114[k]
                  + f_96 * if__134[k]
                  + f_95 * if__167[k]
                  - f_97 * if__187[k];
    }

#pragma omp simd aligned(ab_x, hf_41, hf_46, hf_48, hf_111, hf_116, hf_118, hf_131, hf_136, \
                         hf_138, if__41, if__46, if__48, if__111, if__116, if__118, if__131, \
                         if__136, if__138 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_98 * ab_x[k] * hf_41[k]
                  + f_98 * ab_x[k] * hf_46[k]
                  - f_99 * ab_x[k] * hf_48[k]
                  + f_98 * ab_x[k] * hf_111[k]
                  + f_98 * ab_x[k] * hf_116[k]
                  - f_99 * ab_x[k] * hf_118[k]
                  - f_100 * ab_x[k] * hf_131[k]
                  - f_100 * ab_x[k] * hf_136[k]
                  + f_101 * ab_x[k] * hf_138[k]
                  + f_98 * if__41[k]
                  + f_98 * if__46[k]
                  - f_99 * if__48[k]
                  + f_98 * if__111[k]
                  + f_98 * if__116[k]
                  - f_99 * if__118[k]
                  - f_100 * if__131[k]
                  - f_100 * if__136[k]
                  + f_101 * if__138[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_44, hf_47, hf_49, hf_114, hf_117, hf_119, hf_134, \
                         hf_137, hf_139, if__44, if__77, if__79, if__114, if__134, if__167, \
                         if__169, if__187, if__189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_102 * ab_x[k] * hf_44[k]
                  + f_102 * ab_y[k] * hf_47[k]
                  - f_103 * ab_y[k] * hf_49[k]
                  + f_102 * ab_x[k] * hf_114[k]
                  + f_102 * ab_y[k] * hf_117[k]
                  - f_103 * ab_y[k] * hf_119[k]
                  - f_104 * ab_x[k] * hf_134[k]
                  - f_104 * ab_y[k] * hf_137[k]
                  + f_105 * ab_y[k] * hf_139[k]
                  + f_102 * if__44[k]
                  + f_102 * if__77[k]
                  - f_103 * if__79[k]
                  + f_102 * if__114[k]
                  - f_104 * if__134[k]
                  + f_102 * if__167[k]
                  - f_103 * if__169[k]
                  - f_104 * if__187[k]
                  + f_105 * if__189[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hf_40, hf_43, hf_45, hf_46, hf_48, hf_49, hf_110, \
                         hf_113, hf_115, hf_116, hf_118, hf_119, hf_130, hf_133, hf_135, \
                         hf_136, hf_138, hf_139, if__40, if__43, if__45, if__76, if__78, \
                         if__89, if__110, if__113, if__115, if__130, if__133, if__135, \
                         if__166, if__168, if__179, if__186, if__188, \
                         if__199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_106 * ab_x[k] * hf_40[k]
                  - f_107 * ab_x[k] * hf_43[k]
                  + f_108 * ab_x[k] * hf_45[k]
                  - f_106 * ab_y[k] * hf_46[k]
                  + f_108 * ab_y[k] * hf_48[k]
                  - f_109 * ab_z[k] * hf_49[k]
                  - f_106 * ab_x[k] * hf_110[k]
                  - f_107 * ab_x[k] * hf_113[k]
                  + f_108 * ab_x[k] * hf_115[k]
                  - f_106 * ab_y[k] * hf_116[k]
                  + f_108 * ab_y[k] * hf_118[k]
                  - f_109 * ab_z[k] * hf_119[k]
                  + f_107 * ab_x[k] * hf_130[k]
                  + f_110 * ab_x[k] * hf_133[k]
                  - f_111 * ab_x[k] * hf_135[k]
                  + f_107 * ab_y[k] * hf_136[k]
                  - f_111 * ab_y[k] * hf_138[k]
                  + f_112 * ab_z[k] * hf_139[k]
                  - f_106 * if__40[k]
                  - f_107 * if__43[k]
                  + f_108 * if__45[k]
                  - f_106 * if__76[k]
                  + f_108 * if__78[k]
                  - f_109 * if__89[k]
                  - f_106 * if__110[k]
                  - f_107 * if__113[k]
                  + f_108 * if__115[k]
                  + f_107 * if__130[k]
                  + f_110 * if__133[k]
                  - f_111 * if__135[k]
                  - f_106 * if__166[k]
                  + f_108 * if__168[k]
                  - f_109 * if__179[k]
                  + f_107 * if__186[k]
                  - f_111 * if__188[k]
                  + f_112 * if__199[k];
    }

#pragma omp simd aligned(ab_x, hf_42, hf_47, hf_49, hf_112, hf_117, hf_119, hf_132, hf_137, \
                         hf_139, if__42, if__47, if__49, if__112, if__117, if__119, if__132, \
                         if__137, if__139 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_102 * ab_x[k] * hf_42[k]
                  + f_102 * ab_x[k] * hf_47[k]
                  - f_103 * ab_x[k] * hf_49[k]
                  + f_102 * ab_x[k] * hf_112[k]
                  + f_102 * ab_x[k] * hf_117[k]
                  - f_103 * ab_x[k] * hf_119[k]
                  - f_104 * ab_x[k] * hf_132[k]
                  - f_104 * ab_x[k] * hf_137[k]
                  + f_105 * ab_x[k] * hf_139[k]
                  + f_102 * if__42[k]
                  + f_102 * if__47[k]
                  - f_103 * if__49[k]
                  + f_102 * if__112[k]
                  + f_102 * if__117[k]
                  - f_103 * if__119[k]
                  - f_104 * if__132[k]
                  - f_104 * if__137[k]
                  + f_105 * if__139[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_40, hf_45, hf_46, hf_48, hf_110, hf_115, hf_116, \
                         hf_118, hf_130, hf_135, hf_136, hf_138, if__40, if__45, if__76, \
                         if__78, if__110, if__115, if__130, if__135, if__166, if__168, \
                         if__186, if__188 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_113 * ab_x[k] * hf_40[k]
                  - f_114 * ab_x[k] * hf_45[k]
                  - f_113 * ab_y[k] * hf_46[k]
                  + f_114 * ab_y[k] * hf_48[k]
                  + f_113 * ab_x[k] * hf_110[k]
                  - f_114 * ab_x[k] * hf_115[k]
                  - f_113 * ab_y[k] * hf_116[k]
                  + f_114 * ab_y[k] * hf_118[k]
                  - f_98 * ab_x[k] * hf_130[k]
                  + f_99 * ab_x[k] * hf_135[k]
                  + f_98 * ab_y[k] * hf_136[k]
                  - f_99 * ab_y[k] * hf_138[k]
                  + f_113 * if__40[k]
                  - f_114 * if__45[k]
                  - f_113 * if__76[k]
                  + f_114 * if__78[k]
                  + f_113 * if__110[k]
                  - f_114 * if__115[k]
                  - f_98 * if__130[k]
                  + f_99 * if__135[k]
                  - f_113 * if__166[k]
                  + f_114 * if__168[k]
                  + f_98 * if__186[k]
                  - f_99 * if__188[k];
    }

#pragma omp simd aligned(ab_x, hf_42, hf_47, hf_112, hf_117, hf_132, hf_137, if__42, if__47, \
                         if__112, if__117, if__132, if__137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_95 * ab_x[k] * hf_42[k]
                  + f_94 * ab_x[k] * hf_47[k]
                  - f_95 * ab_x[k] * hf_112[k]
                  + f_94 * ab_x[k] * hf_117[k]
                  + f_97 * ab_x[k] * hf_132[k]
                  - f_96 * ab_x[k] * hf_137[k]
                  - f_95 * if__42[k]
                  + f_94 * if__47[k]
                  - f_95 * if__112[k]
                  + f_94 * if__117[k]
                  + f_97 * if__132[k]
                  - f_96 * if__137[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_40, hf_43, hf_46, hf_110, hf_113, hf_116, hf_130, \
                         hf_133, hf_136, if__40, if__43, if__76, if__110, if__113, if__130, \
                         if__133, if__166, if__186 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_115 * ab_x[k] * hf_40[k]
                  + f_116 * ab_x[k] * hf_43[k]
                  - f_115 * ab_y[k] * hf_46[k]
                  - f_115 * ab_x[k] * hf_110[k]
                  + f_116 * ab_x[k] * hf_113[k]
                  - f_115 * ab_y[k] * hf_116[k]
                  + f_117 * ab_x[k] * hf_130[k]
                  - f_118 * ab_x[k] * hf_133[k]
                  + f_117 * ab_y[k] * hf_136[k]
                  - f_115 * if__40[k]
                  + f_116 * if__43[k]
                  - f_115 * if__76[k]
                  - f_115 * if__110[k]
                  + f_116 * if__113[k]
                  + f_117 * if__130[k]
                  - f_118 * if__133[k]
                  - f_115 * if__166[k]
                  + f_117 * if__186[k];
    }

#pragma omp simd aligned(ab_x, hf_11, hf_16, hf_61, hf_66, hf_81, hf_86, hf_151, hf_156, \
                         hf_171, hf_176, hf_191, hf_196, if__11, if__16, if__61, if__66, \
                         if__81, if__86, if__151, if__156, if__171, if__176, if__191, \
                         if__196 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_119 * ab_x[k] * hf_11[k]
                  - f_119 * ab_x[k] * hf_16[k]
                  + f_113 * ab_x[k] * hf_61[k]
                  - f_113 * ab_x[k] * hf_66[k]
                  - f_114 * ab_x[k] * hf_81[k]
                  + f_114 * ab_x[k] * hf_86[k]
                  + f_119 * ab_x[k] * hf_151[k]
                  - f_119 * ab_x[k] * hf_156[k]
                  - f_114 * ab_x[k] * hf_171[k]
                  + f_114 * ab_x[k] * hf_176[k]
                  + f_100 * ab_x[k] * hf_191[k]
                  - f_100 * ab_x[k] * hf_196[k]
                  + f_119 * if__11[k]
                  - f_119 * if__16[k]
                  + f_113 * if__61[k]
                  - f_113 * if__66[k]
                  - f_114 * if__81[k]
                  + f_114 * if__86[k]
                  + f_119 * if__151[k]
                  - f_119 * if__156[k]
                  - f_114 * if__171[k]
                  + f_114 * if__176[k]
                  + f_100 * if__191[k]
                  - f_100 * if__196[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_14, hf_17, hf_64, hf_67, hf_84, hf_87, hf_154, hf_157, \
                         hf_174, hf_177, hf_194, hf_197, if__14, if__37, if__64, if__84, \
                         if__107, if__127, if__154, if__174, if__194, if__217, if__237, \
                         if__257 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_120 * ab_x[k] * hf_14[k]
                  - f_121 * ab_y[k] * hf_17[k]
                  + f_122 * ab_x[k] * hf_64[k]
                  - f_123 * ab_y[k] * hf_67[k]
                  - f_124 * ab_x[k] * hf_84[k]
                  + f_102 * ab_y[k] * hf_87[k]
                  + f_120 * ab_x[k] * hf_154[k]
                  - f_121 * ab_y[k] * hf_157[k]
                  - f_124 * ab_x[k] * hf_174[k]
                  + f_102 * ab_y[k] * hf_177[k]
                  + f_104 * ab_x[k] * hf_194[k]
                  - f_125 * ab_y[k] * hf_197[k]
                  + f_120 * if__14[k]
                  - f_121 * if__37[k]
                  + f_122 * if__64[k]
                  - f_124 * if__84[k]
                  - f_123 * if__107[k]
                  + f_102 * if__127[k]
                  + f_120 * if__154[k]
                  - f_124 * if__174[k]
                  + f_104 * if__194[k]
                  - f_121 * if__217[k]
                  + f_102 * if__237[k]
                  - f_125 * if__257[k];
    }

#pragma omp simd aligned(ab_x, hf_11, hf_16, hf_18, hf_61, hf_66, hf_68, hf_81, hf_86, hf_88, \
                         hf_151, hf_156, hf_158, hf_171, hf_176, hf_178, hf_191, hf_196, \
                         hf_198, if__11, if__16, if__18, if__61, if__66, if__68, if__81, \
                         if__86, if__88, if__151, if__156, if__158, if__171, if__176, if__178, \
                         if__191, if__196, if__198 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_126 * ab_x[k] * hf_11[k]
                  - f_126 * ab_x[k] * hf_16[k]
                  + f_127 * ab_x[k] * hf_18[k]
                  - f_128 * ab_x[k] * hf_61[k]
                  - f_128 * ab_x[k] * hf_66[k]
                  + f_129 * ab_x[k] * hf_68[k]
                  + f_129 * ab_x[k] * hf_81[k]
                  + f_129 * ab_x[k] * hf_86[k]
                  - f_130 * ab_x[k] * hf_88[k]
                  - f_126 * ab_x[k] * hf_151[k]
                  - f_126 * ab_x[k] * hf_156[k]
                  + f_127 * ab_x[k] * hf_158[k]
                  + f_129 * ab_x[k] * hf_171[k]
                  + f_129 * ab_x[k] * hf_176[k]
                  - f_130 * ab_x[k] * hf_178[k]
                  - f_131 * ab_x[k] * hf_191[k]
                  - f_131 * ab_x[k] * hf_196[k]
                  + f_132 * ab_x[k] * hf_198[k]
                  - f_126 * if__11[k]
                  - f_126 * if__16[k]
                  + f_127 * if__18[k]
                  - f_128 * if__61[k]
                  - f_128 * if__66[k]
                  + f_129 * if__68[k]
                  + f_129 * if__81[k]
                  + f_129 * if__86[k]
                  - f_130 * if__88[k]
                  - f_126 * if__151[k]
                  - f_126 * if__156[k]
                  + f_127 * if__158[k]
                  + f_129 * if__171[k]
                  + f_129 * if__176[k]
                  - f_130 * if__178[k]
                  - f_131 * if__191[k]
                  - f_131 * if__196[k]
                  + f_132 * if__198[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_14, hf_17, hf_19, hf_64, hf_67, hf_69, hf_84, hf_87, \
                         hf_89, hf_154, hf_157, hf_159, hf_174, hf_177, hf_179, hf_194, \
                         hf_197, hf_199, if__14, if__37, if__39, if__64, if__84, if__107, \
                         if__109, if__127, if__129, if__154, if__174, if__194, if__217, \
                         if__219, if__237, if__239, if__257, if__259 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_133 * ab_x[k] * hf_14[k]
                  - f_133 * ab_y[k] * hf_17[k]
                  + f_134 * ab_y[k] * hf_19[k]
                  - f_135 * ab_x[k] * hf_64[k]
                  - f_135 * ab_y[k] * hf_67[k]
                  + f_136 * ab_y[k] * hf_69[k]
                  + f_137 * ab_x[k] * hf_84[k]
                  + f_137 * ab_y[k] * hf_87[k]
                  - f_138 * ab_y[k] * hf_89[k]
                  - f_133 * ab_x[k] * hf_154[k]
                  - f_133 * ab_y[k] * hf_157[k]
                  + f_134 * ab_y[k] * hf_159[k]
                  + f_137 * ab_x[k] * hf_174[k]
                  + f_137 * ab_y[k] * hf_177[k]
                  - f_138 * ab_y[k] * hf_179[k]
                  - f_139 * ab_x[k] * hf_194[k]
                  - f_139 * ab_y[k] * hf_197[k]
                  + f_140 * ab_y[k] * hf_199[k]
                  - f_133 * if__14[k]
                  - f_133 * if__37[k]
                  + f_134 * if__39[k]
                  - f_135 * if__64[k]
                  + f_137 * if__84[k]
                  - f_135 * if__107[k]
                  + f_136 * if__109[k]
                  + f_137 * if__127[k]
                  - f_138 * if__129[k]
                  - f_133 * if__154[k]
                  + f_137 * if__174[k]
                  - f_139 * if__194[k]
                  - f_133 * if__217[k]
                  + f_134 * if__219[k]
                  + f_137 * if__237[k]
                  - f_138 * if__239[k]
                  - f_139 * if__257[k]
                  + f_140 * if__259[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hf_10, hf_13, hf_15, hf_16, hf_18, hf_19, hf_60, \
                         hf_63, hf_65, hf_66, hf_68, hf_69, hf_80, hf_83, hf_85, hf_86, hf_88, \
                         hf_89, hf_150, hf_153, hf_155, hf_156, hf_158, hf_159, hf_170, \
                         hf_173, hf_175, hf_176, hf_178, hf_179, hf_190, hf_193, hf_195, \
                         hf_196, hf_198, hf_199, if__10, if__13, if__15, if__36, if__38, \
                         if__49, if__60, if__63, if__65, if__80, if__83, if__85, if__106, \
                         if__108, if__119, if__126, if__128, if__139, if__150, if__153, \
                         if__155, if__170, if__173, if__175, if__190, if__193, if__195, \
                         if__216, if__218, if__229, if__236, if__238, if__249, if__256, \
                         if__258, if__269 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_141 * ab_x[k] * hf_10[k]
                  + f_142 * ab_x[k] * hf_13[k]
                  - f_143 * ab_x[k] * hf_15[k]
                  + f_141 * ab_y[k] * hf_16[k]
                  - f_143 * ab_y[k] * hf_18[k]
                  + f_144 * ab_z[k] * hf_19[k]
                  + f_142 * ab_x[k] * hf_60[k]
                  + f_145 * ab_x[k] * hf_63[k]
                  - f_146 * ab_x[k] * hf_65[k]
                  + f_142 * ab_y[k] * hf_66[k]
                  - f_146 * ab_y[k] * hf_68[k]
                  + f_147 * ab_z[k] * hf_69[k]
                  - f_148 * ab_x[k] * hf_80[k]
                  - f_149 * ab_x[k] * hf_83[k]
                  + f_150 * ab_x[k] * hf_85[k]
                  - f_148 * ab_y[k] * hf_86[k]
                  + f_150 * ab_y[k] * hf_88[k]
                  - f_151 * ab_z[k] * hf_89[k]
                  + f_141 * ab_x[k] * hf_150[k]
                  + f_142 * ab_x[k] * hf_153[k]
                  - f_143 * ab_x[k] * hf_155[k]
                  + f_141 * ab_y[k] * hf_156[k]
                  - f_143 * ab_y[k] * hf_158[k]
                  + f_144 * ab_z[k] * hf_159[k]
                  - f_148 * ab_x[k] * hf_170[k]
                  - f_149 * ab_x[k] * hf_173[k]
                  + f_150 * ab_x[k] * hf_175[k]
                  - f_148 * ab_y[k] * hf_176[k]
                  + f_150 * ab_y[k] * hf_178[k]
                  - f_151 * ab_z[k] * hf_179[k]
                  + f_143 * ab_x[k] * hf_190[k]
                  + f_146 * ab_x[k] * hf_193[k]
                  - f_152 * ab_x[k] * hf_195[k]
                  + f_143 * ab_y[k] * hf_196[k]
                  - f_152 * ab_y[k] * hf_198[k]
                  + f_153 * ab_z[k] * hf_199[k]
                  + f_141 * if__10[k]
                  + f_142 * if__13[k]
                  - f_143 * if__15[k]
                  + f_141 * if__36[k]
                  - f_143 * if__38[k]
                  + f_144 * if__49[k]
                  + f_142 * if__60[k]
                  + f_145 * if__63[k]
                  - f_146 * if__65[k]
                  - f_148 * if__80[k]
                  - f_149 * if__83[k]
                  + f_150 * if__85[k]
                  + f_142 * if__106[k]
                  - f_146 * if__108[k]
                  + f_147 * if__119[k]
                  - f_148 * if__126[k]
                  + f_150 * if__128[k]
                  - f_151 * if__139[k]
                  + f_141 * if__150[k]
                  + f_142 * if__153[k]
                  - f_143 * if__155[k]
                  - f_148 * if__170[k]
                  - f_149 * if__173[k]
                  + f_150 * if__175[k]
                  + f_143 * if__190[k]
                  + f_146 * if__193[k]
                  - f_152 * if__195[k]
                  + f_141 * if__216[k]
                  - f_143 * if__218[k]
                  + f_144 * if__229[k]
                  - f_148 * if__236[k]
                  + f_150 * if__238[k]
                  - f_151 * if__249[k]
                  + f_143 * if__256[k]
                  - f_152 * if__258[k]
                  + f_153 * if__269[k];
    }

#pragma omp simd aligned(ab_x, hf_12, hf_17, hf_19, hf_62, hf_67, hf_69, hf_82, hf_87, hf_89, \
                         hf_152, hf_157, hf_159, hf_172, hf_177, hf_179, hf_192, hf_197, \
                         hf_199, if__12, if__17, if__19, if__62, if__67, if__69, if__82, \
                         if__87, if__89, if__152, if__157, if__159, if__172, if__177, if__179, \
                         if__192, if__197, if__199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_133 * ab_x[k] * hf_12[k]
                  - f_133 * ab_x[k] * hf_17[k]
                  + f_134 * ab_x[k] * hf_19[k]
                  - f_135 * ab_x[k] * hf_62[k]
                  - f_135 * ab_x[k] * hf_67[k]
                  + f_136 * ab_x[k] * hf_69[k]
                  + f_137 * ab_x[k] * hf_82[k]
                  + f_137 * ab_x[k] * hf_87[k]
                  - f_138 * ab_x[k] * hf_89[k]
                  - f_133 * ab_x[k] * hf_152[k]
                  - f_133 * ab_x[k] * hf_157[k]
                  + f_134 * ab_x[k] * hf_159[k]
                  + f_137 * ab_x[k] * hf_172[k]
                  + f_137 * ab_x[k] * hf_177[k]
                  - f_138 * ab_x[k] * hf_179[k]
                  - f_139 * ab_x[k] * hf_192[k]
                  - f_139 * ab_x[k] * hf_197[k]
                  + f_140 * ab_x[k] * hf_199[k]
                  - f_133 * if__12[k]
                  - f_133 * if__17[k]
                  + f_134 * if__19[k]
                  - f_135 * if__62[k]
                  - f_135 * if__67[k]
                  + f_136 * if__69[k]
                  + f_137 * if__82[k]
                  + f_137 * if__87[k]
                  - f_138 * if__89[k]
                  - f_133 * if__152[k]
                  - f_133 * if__157[k]
                  + f_134 * if__159[k]
                  + f_137 * if__172[k]
                  + f_137 * if__177[k]
                  - f_138 * if__179[k]
                  - f_139 * if__192[k]
                  - f_139 * if__197[k]
                  + f_140 * if__199[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_10, hf_15, hf_16, hf_18, hf_60, hf_65, hf_66, hf_68, \
                         hf_80, hf_85, hf_86, hf_88, hf_150, hf_155, hf_156, hf_158, hf_170, \
                         hf_175, hf_176, hf_178, hf_190, hf_195, hf_196, hf_198, if__10, \
                         if__15, if__36, if__38, if__60, if__65, if__80, if__85, if__106, \
                         if__108, if__126, if__128, if__150, if__155, if__170, if__175, \
                         if__190, if__195, if__216, if__218, if__236, if__238, if__256, \
                         if__258 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_154 * ab_x[k] * hf_10[k]
                  + f_155 * ab_x[k] * hf_15[k]
                  + f_154 * ab_y[k] * hf_16[k]
                  - f_155 * ab_y[k] * hf_18[k]
                  - f_126 * ab_x[k] * hf_60[k]
                  + f_127 * ab_x[k] * hf_65[k]
                  + f_126 * ab_y[k] * hf_66[k]
                  - f_127 * ab_y[k] * hf_68[k]
                  + f_127 * ab_x[k] * hf_80[k]
                  - f_156 * ab_x[k] * hf_85[k]
                  - f_127 * ab_y[k] * hf_86[k]
                  + f_156 * ab_y[k] * hf_88[k]
                  - f_154 * ab_x[k] * hf_150[k]
                  + f_155 * ab_x[k] * hf_155[k]
                  + f_154 * ab_y[k] * hf_156[k]
                  - f_155 * ab_y[k] * hf_158[k]
                  + f_127 * ab_x[k] * hf_170[k]
                  - f_156 * ab_x[k] * hf_175[k]
                  - f_127 * ab_y[k] * hf_176[k]
                  + f_156 * ab_y[k] * hf_178[k]
                  - f_157 * ab_x[k] * hf_190[k]
                  + f_158 * ab_x[k] * hf_195[k]
                  + f_157 * ab_y[k] * hf_196[k]
                  - f_158 * ab_y[k] * hf_198[k]
                  - f_154 * if__10[k]
                  + f_155 * if__15[k]
                  + f_154 * if__36[k]
                  - f_155 * if__38[k]
                  - f_126 * if__60[k]
                  + f_127 * if__65[k]
                  + f_127 * if__80[k]
                  - f_156 * if__85[k]
                  + f_126 * if__106[k]
                  - f_127 * if__108[k]
                  - f_127 * if__126[k]
                  + f_156 * if__128[k]
                  - f_154 * if__150[k]
                  + f_155 * if__155[k]
                  + f_127 * if__170[k]
                  - f_156 * if__175[k]
                  - f_157 * if__190[k]
                  + f_158 * if__195[k]
                  + f_154 * if__216[k]
                  - f_155 * if__218[k]
                  - f_127 * if__236[k]
                  + f_156 * if__238[k]
                  + f_157 * if__256[k]
                  - f_158 * if__258[k];
    }

#pragma omp simd aligned(ab_x, hf_12, hf_17, hf_62, hf_67, hf_82, hf_87, hf_152, hf_157, \
                         hf_172, hf_177, hf_192, hf_197, if__12, if__17, if__62, if__67, \
                         if__82, if__87, if__152, if__157, if__172, if__177, if__192, \
                         if__197 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_121 * ab_x[k] * hf_12[k]
                  - f_120 * ab_x[k] * hf_17[k]
                  + f_123 * ab_x[k] * hf_62[k]
                  - f_122 * ab_x[k] * hf_67[k]
                  - f_102 * ab_x[k] * hf_82[k]
                  + f_124 * ab_x[k] * hf_87[k]
                  + f_121 * ab_x[k] * hf_152[k]
                  - f_120 * ab_x[k] * hf_157[k]
                  - f_102 * ab_x[k] * hf_172[k]
                  + f_124 * ab_x[k] * hf_177[k]
                  + f_125 * ab_x[k] * hf_192[k]
                  - f_104 * ab_x[k] * hf_197[k]
                  + f_121 * if__12[k]
                  - f_120 * if__17[k]
                  + f_123 * if__62[k]
                  - f_122 * if__67[k]
                  - f_102 * if__82[k]
                  + f_124 * if__87[k]
                  + f_121 * if__152[k]
                  - f_120 * if__157[k]
                  - f_102 * if__172[k]
                  + f_124 * if__177[k]
                  + f_125 * if__192[k]
                  - f_104 * if__197[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_10, hf_13, hf_16, hf_60, hf_63, hf_66, hf_80, hf_83, \
                         hf_86, hf_150, hf_153, hf_156, hf_170, hf_173, hf_176, hf_190, \
                         hf_193, hf_196, if__10, if__13, if__36, if__60, if__63, if__80, \
                         if__83, if__106, if__126, if__150, if__153, if__170, if__173, \
                         if__190, if__193, if__216, if__236, if__256 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_159 * ab_x[k] * hf_10[k]
                  - f_160 * ab_x[k] * hf_13[k]
                  + f_159 * ab_y[k] * hf_16[k]
                  + f_161 * ab_x[k] * hf_60[k]
                  - f_162 * ab_x[k] * hf_63[k]
                  + f_161 * ab_y[k] * hf_66[k]
                  - f_162 * ab_x[k] * hf_80[k]
                  + f_163 * ab_x[k] * hf_83[k]
                  - f_162 * ab_y[k] * hf_86[k]
                  + f_159 * ab_x[k] * hf_150[k]
                  - f_160 * ab_x[k] * hf_153[k]
                  + f_159 * ab_y[k] * hf_156[k]
                  - f_162 * ab_x[k] * hf_170[k]
                  + f_163 * ab_x[k] * hf_173[k]
                  - f_162 * ab_y[k] * hf_176[k]
                  + f_113 * ab_x[k] * hf_190[k]
                  - f_114 * ab_x[k] * hf_193[k]
                  + f_113 * ab_y[k] * hf_196[k]
                  + f_159 * if__10[k]
                  - f_160 * if__13[k]
                  + f_159 * if__36[k]
                  + f_161 * if__60[k]
                  - f_162 * if__63[k]
                  - f_162 * if__80[k]
                  + f_163 * if__83[k]
                  + f_161 * if__106[k]
                  - f_162 * if__126[k]
                  + f_159 * if__150[k]
                  - f_160 * if__153[k]
                  - f_162 * if__170[k]
                  + f_163 * if__173[k]
                  + f_113 * if__190[k]
                  - f_114 * if__193[k]
                  + f_159 * if__216[k]
                  - f_162 * if__236[k]
                  + f_113 * if__256[k];
    }

#pragma omp simd aligned(ab_x, hf_21, hf_26, hf_71, hf_76, hf_91, hf_96, hf_161, hf_166, \
                         hf_181, hf_186, hf_201, hf_206, if__21, if__26, if__71, if__76, \
                         if__91, if__96, if__161, if__166, if__181, if__186, if__201, \
                         if__206 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_164 * ab_x[k] * hf_21[k]
                  - f_164 * ab_x[k] * hf_26[k]
                  + f_16 * ab_x[k] * hf_71[k]
                  - f_16 * ab_x[k] * hf_76[k]
                  - f_165 * ab_x[k] * hf_91[k]
                  + f_165 * ab_x[k] * hf_96[k]
                  + f_164 * ab_x[k] * hf_161[k]
                  - f_164 * ab_x[k] * hf_166[k]
                  - f_165 * ab_x[k] * hf_181[k]
                  + f_165 * ab_x[k] * hf_186[k]
                  + f_166 * ab_x[k] * hf_201[k]
                  - f_166 * ab_x[k] * hf_206[k]
                  + f_164 * if__21[k]
                  - f_164 * if__26[k]
                  + f_16 * if__71[k]
                  - f_16 * if__76[k]
                  - f_165 * if__91[k]
                  + f_165 * if__96[k]
                  + f_164 * if__161[k]
                  - f_164 * if__166[k]
                  - f_165 * if__181[k]
                  + f_165 * if__186[k]
                  + f_166 * if__201[k]
                  - f_166 * if__206[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_24, hf_27, hf_74, hf_77, hf_94, hf_97, hf_164, hf_167, \
                         hf_184, hf_187, hf_204, hf_207, if__24, if__47, if__74, if__94, \
                         if__117, if__137, if__164, if__184, if__204, if__227, if__247, \
                         if__267 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_33 * ab_x[k] * hf_24[k]
                  - f_9 * ab_y[k] * hf_27[k]
                  + f_10 * ab_x[k] * hf_74[k]
                  - f_11 * ab_y[k] * hf_77[k]
                  - f_167 * ab_x[k] * hf_94[k]
                  + f_168 * ab_y[k] * hf_97[k]
                  + f_33 * ab_x[k] * hf_164[k]
                  - f_9 * ab_y[k] * hf_167[k]
                  - f_167 * ab_x[k] * hf_184[k]
                  + f_168 * ab_y[k] * hf_187[k]
                  + f_169 * ab_x[k] * hf_204[k]
                  - f_170 * ab_y[k] * hf_207[k]
                  + f_33 * if__24[k]
                  - f_9 * if__47[k]
                  + f_10 * if__74[k]
                  - f_167 * if__94[k]
                  - f_11 * if__117[k]
                  + f_168 * if__137[k]
                  + f_33 * if__164[k]
                  - f_167 * if__184[k]
                  + f_169 * if__204[k]
                  - f_9 * if__227[k]
                  + f_168 * if__247[k]
                  - f_170 * if__267[k];
    }

#pragma omp simd aligned(ab_x, hf_21, hf_26, hf_28, hf_71, hf_76, hf_78, hf_91, hf_96, hf_98, \
                         hf_161, hf_166, hf_168, hf_181, hf_186, hf_188, hf_201, hf_206, \
                         hf_208, if__21, if__26, if__28, if__71, if__76, if__78, if__91, \
                         if__96, if__98, if__161, if__166, if__168, if__181, if__186, if__188, \
                         if__201, if__206, if__208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_171 * ab_x[k] * hf_21[k]
                  - f_171 * ab_x[k] * hf_26[k]
                  + f_172 * ab_x[k] * hf_28[k]
                  - f_173 * ab_x[k] * hf_71[k]
                  - f_173 * ab_x[k] * hf_76[k]
                  + f_174 * ab_x[k] * hf_78[k]
                  + f_175 * ab_x[k] * hf_91[k]
                  + f_175 * ab_x[k] * hf_96[k]
                  - f_176 * ab_x[k] * hf_98[k]
                  - f_171 * ab_x[k] * hf_161[k]
                  - f_171 * ab_x[k] * hf_166[k]
                  + f_172 * ab_x[k] * hf_168[k]
                  + f_175 * ab_x[k] * hf_181[k]
                  + f_175 * ab_x[k] * hf_186[k]
                  - f_176 * ab_x[k] * hf_188[k]
                  - f_177 * ab_x[k] * hf_201[k]
                  - f_177 * ab_x[k] * hf_206[k]
                  + f_178 * ab_x[k] * hf_208[k]
                  - f_171 * if__21[k]
                  - f_171 * if__26[k]
                  + f_172 * if__28[k]
                  - f_173 * if__71[k]
                  - f_173 * if__76[k]
                  + f_174 * if__78[k]
                  + f_175 * if__91[k]
                  + f_175 * if__96[k]
                  - f_176 * if__98[k]
                  - f_171 * if__161[k]
                  - f_171 * if__166[k]
                  + f_172 * if__168[k]
                  + f_175 * if__181[k]
                  + f_175 * if__186[k]
                  - f_176 * if__188[k]
                  - f_177 * if__201[k]
                  - f_177 * if__206[k]
                  + f_178 * if__208[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_24, hf_27, hf_29, hf_74, hf_77, hf_79, hf_94, hf_97, \
                         hf_99, hf_164, hf_167, hf_169, hf_184, hf_187, hf_189, hf_204, \
                         hf_207, hf_209, if__24, if__47, if__49, if__74, if__94, if__117, \
                         if__119, if__137, if__139, if__164, if__184, if__204, if__227, \
                         if__229, if__247, if__249, if__267, if__269 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_179 * ab_x[k] * hf_24[k]
                  - f_179 * ab_y[k] * hf_27[k]
                  + f_180 * ab_y[k] * hf_29[k]
                  - f_181 * ab_x[k] * hf_74[k]
                  - f_181 * ab_y[k] * hf_77[k]
                  + f_182 * ab_y[k] * hf_79[k]
                  + f_182 * ab_x[k] * hf_94[k]
                  + f_182 * ab_y[k] * hf_97[k]
                  - f_183 * ab_y[k] * hf_99[k]
                  - f_179 * ab_x[k] * hf_164[k]
                  - f_179 * ab_y[k] * hf_167[k]
                  + f_180 * ab_y[k] * hf_169[k]
                  + f_182 * ab_x[k] * hf_184[k]
                  + f_182 * ab_y[k] * hf_187[k]
                  - f_183 * ab_y[k] * hf_189[k]
                  - f_184 * ab_x[k] * hf_204[k]
                  - f_184 * ab_y[k] * hf_207[k]
                  + f_185 * ab_y[k] * hf_209[k]
                  - f_179 * if__24[k]
                  - f_179 * if__47[k]
                  + f_180 * if__49[k]
                  - f_181 * if__74[k]
                  + f_182 * if__94[k]
                  - f_181 * if__117[k]
                  + f_182 * if__119[k]
                  + f_182 * if__137[k]
                  - f_183 * if__139[k]
                  - f_179 * if__164[k]
                  + f_182 * if__184[k]
                  - f_184 * if__204[k]
                  - f_179 * if__227[k]
                  + f_180 * if__229[k]
                  + f_182 * if__247[k]
                  - f_183 * if__249[k]
                  - f_184 * if__267[k]
                  + f_185 * if__269[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hf_20, hf_23, hf_25, hf_26, hf_28, hf_29, hf_70, \
                         hf_73, hf_75, hf_76, hf_78, hf_79, hf_90, hf_93, hf_95, hf_96, hf_98, \
                         hf_99, hf_160, hf_163, hf_165, hf_166, hf_168, hf_169, hf_180, \
                         hf_183, hf_185, hf_186, hf_188, hf_189, hf_200, hf_203, hf_205, \
                         hf_206, hf_208, hf_209, if__20, if__23, if__25, if__46, if__48, \
                         if__59, if__70, if__73, if__75, if__90, if__93, if__95, if__116, \
                         if__118, if__129, if__136, if__138, if__149, if__160, if__163, \
                         if__165, if__180, if__183, if__185, if__200, if__203, if__205, \
                         if__226, if__228, if__239, if__246, if__248, if__259, if__266, \
                         if__268, if__279 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = 0.703125 * ab_x[k] * hf_20[k]
                  + 1.40625 * ab_x[k] * hf_23[k]
                  - 5.625 * ab_x[k] * hf_25[k]
                  + 0.703125 * ab_y[k] * hf_26[k]
                  - 5.625 * ab_y[k] * hf_28[k]
                  + 1.875 * ab_z[k] * hf_29[k]
                  + 1.40625 * ab_x[k] * hf_70[k]
                  + 2.8125 * ab_x[k] * hf_73[k]
                  - 11.25 * ab_x[k] * hf_75[k]
                  + 1.40625 * ab_y[k] * hf_76[k]
                  - 11.25 * ab_y[k] * hf_78[k]
                  + 3.75 * ab_z[k] * hf_79[k]
                  - 1.875 * ab_x[k] * hf_90[k]
                  - 3.75 * ab_x[k] * hf_93[k]
                  + 15.0 * ab_x[k] * hf_95[k]
                  - 1.875 * ab_y[k] * hf_96[k]
                  + 15.0 * ab_y[k] * hf_98[k]
                  - 5.0 * ab_z[k] * hf_99[k]
                  + 0.703125 * ab_x[k] * hf_160[k]
                  + 1.40625 * ab_x[k] * hf_163[k]
                  - 5.625 * ab_x[k] * hf_165[k]
                  + 0.703125 * ab_y[k] * hf_166[k]
                  - 5.625 * ab_y[k] * hf_168[k]
                  + 1.875 * ab_z[k] * hf_169[k]
                  - 1.875 * ab_x[k] * hf_180[k]
                  - 3.75 * ab_x[k] * hf_183[k]
                  + 15.0 * ab_x[k] * hf_185[k]
                  - 1.875 * ab_y[k] * hf_186[k]
                  + 15.0 * ab_y[k] * hf_188[k]
                  - 5.0 * ab_z[k] * hf_189[k]
                  + 0.375 * ab_x[k] * hf_200[k]
                  + 0.75 * ab_x[k] * hf_203[k]
                  - 3.0 * ab_x[k] * hf_205[k]
                  + 0.375 * ab_y[k] * hf_206[k]
                  - 3.0 * ab_y[k] * hf_208[k]
                  + ab_z[k] * hf_209[k]
                  + 0.703125 * if__20[k]
                  + 1.40625 * if__23[k]
                  - 5.625 * if__25[k]
                  + 0.703125 * if__46[k]
                  - 5.625 * if__48[k]
                  + 1.875 * if__59[k]
                  + 1.40625 * if__70[k]
                  + 2.8125 * if__73[k]
                  - 11.25 * if__75[k]
                  - 1.875 * if__90[k]
                  - 3.75 * if__93[k]
                  + 15.0 * if__95[k]
                  + 1.40625 * if__116[k]
                  - 11.25 * if__118[k]
                  + 3.75 * if__129[k]
                  - 1.875 * if__136[k]
                  + 15.0 * if__138[k]
                  - 5.0 * if__149[k]
                  + 0.703125 * if__160[k]
                  + 1.40625 * if__163[k]
                  - 5.625 * if__165[k]
                  - 1.875 * if__180[k]
                  - 3.75 * if__183[k]
                  + 15.0 * if__185[k]
                  + 0.375 * if__200[k]
                  + 0.75 * if__203[k]
                  - 3.0 * if__205[k]
                  + 0.703125 * if__226[k]
                  - 5.625 * if__228[k]
                  + 1.875 * if__239[k]
                  - 1.875 * if__246[k]
                  + 15.0 * if__248[k]
                  - 5.0 * if__259[k]
                  + 0.375 * if__266[k]
                  - 3.0 * if__268[k]
                  + if__279[k];
    }

#pragma omp simd aligned(ab_x, hf_22, hf_27, hf_29, hf_72, hf_77, hf_79, hf_92, hf_97, hf_99, \
                         hf_162, hf_167, hf_169, hf_182, hf_187, hf_189, hf_202, hf_207, \
                         hf_209, if__22, if__27, if__29, if__72, if__77, if__79, if__92, \
                         if__97, if__99, if__162, if__167, if__169, if__182, if__187, if__189, \
                         if__202, if__207, if__209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_179 * ab_x[k] * hf_22[k]
                  - f_179 * ab_x[k] * hf_27[k]
                  + f_180 * ab_x[k] * hf_29[k]
                  - f_181 * ab_x[k] * hf_72[k]
                  - f_181 * ab_x[k] * hf_77[k]
                  + f_182 * ab_x[k] * hf_79[k]
                  + f_182 * ab_x[k] * hf_92[k]
                  + f_182 * ab_x[k] * hf_97[k]
                  - f_183 * ab_x[k] * hf_99[k]
                  - f_179 * ab_x[k] * hf_162[k]
                  - f_179 * ab_x[k] * hf_167[k]
                  + f_180 * ab_x[k] * hf_169[k]
                  + f_182 * ab_x[k] * hf_182[k]
                  + f_182 * ab_x[k] * hf_187[k]
                  - f_183 * ab_x[k] * hf_189[k]
                  - f_184 * ab_x[k] * hf_202[k]
                  - f_184 * ab_x[k] * hf_207[k]
                  + f_185 * ab_x[k] * hf_209[k]
                  - f_179 * if__22[k]
                  - f_179 * if__27[k]
                  + f_180 * if__29[k]
                  - f_181 * if__72[k]
                  - f_181 * if__77[k]
                  + f_182 * if__79[k]
                  + f_182 * if__92[k]
                  + f_182 * if__97[k]
                  - f_183 * if__99[k]
                  - f_179 * if__162[k]
                  - f_179 * if__167[k]
                  + f_180 * if__169[k]
                  + f_182 * if__182[k]
                  + f_182 * if__187[k]
                  - f_183 * if__189[k]
                  - f_184 * if__202[k]
                  - f_184 * if__207[k]
                  + f_185 * if__209[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_20, hf_25, hf_26, hf_28, hf_70, hf_75, hf_76, hf_78, \
                         hf_90, hf_95, hf_96, hf_98, hf_160, hf_165, hf_166, hf_168, hf_180, \
                         hf_185, hf_186, hf_188, hf_200, hf_205, hf_206, hf_208, if__20, \
                         if__25, if__46, if__48, if__70, if__75, if__90, if__95, if__116, \
                         if__118, if__136, if__138, if__160, if__165, if__180, if__185, \
                         if__200, if__205, if__226, if__228, if__246, if__248, if__266, \
                         if__268 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_186 * ab_x[k] * hf_20[k]
                  + f_187 * ab_x[k] * hf_25[k]
                  + f_186 * ab_y[k] * hf_26[k]
                  - f_187 * ab_y[k] * hf_28[k]
                  - f_171 * ab_x[k] * hf_70[k]
                  + f_172 * ab_x[k] * hf_75[k]
                  + f_171 * ab_y[k] * hf_76[k]
                  - f_172 * ab_y[k] * hf_78[k]
                  + f_188 * ab_x[k] * hf_90[k]
                  - f_189 * ab_x[k] * hf_95[k]
                  - f_188 * ab_y[k] * hf_96[k]
                  + f_189 * ab_y[k] * hf_98[k]
                  - f_186 * ab_x[k] * hf_160[k]
                  + f_187 * ab_x[k] * hf_165[k]
                  + f_186 * ab_y[k] * hf_166[k]
                  - f_187 * ab_y[k] * hf_168[k]
                  + f_188 * ab_x[k] * hf_180[k]
                  - f_189 * ab_x[k] * hf_185[k]
                  - f_188 * ab_y[k] * hf_186[k]
                  + f_189 * ab_y[k] * hf_188[k]
                  - f_190 * ab_x[k] * hf_200[k]
                  + f_191 * ab_x[k] * hf_205[k]
                  + f_190 * ab_y[k] * hf_206[k]
                  - f_191 * ab_y[k] * hf_208[k]
                  - f_186 * if__20[k]
                  + f_187 * if__25[k]
                  + f_186 * if__46[k]
                  - f_187 * if__48[k]
                  - f_171 * if__70[k]
                  + f_172 * if__75[k]
                  + f_188 * if__90[k]
                  - f_189 * if__95[k]
                  + f_171 * if__116[k]
                  - f_172 * if__118[k]
                  - f_188 * if__136[k]
                  + f_189 * if__138[k]
                  - f_186 * if__160[k]
                  + f_187 * if__165[k]
                  + f_188 * if__180[k]
                  - f_189 * if__185[k]
                  - f_190 * if__200[k]
                  + f_191 * if__205[k]
                  + f_186 * if__226[k]
                  - f_187 * if__228[k]
                  - f_188 * if__246[k]
                  + f_189 * if__248[k]
                  + f_190 * if__266[k]
                  - f_191 * if__268[k];
    }

#pragma omp simd aligned(ab_x, hf_22, hf_27, hf_72, hf_77, hf_92, hf_97, hf_162, hf_167, \
                         hf_182, hf_187, hf_202, hf_207, if__22, if__27, if__72, if__77, \
                         if__92, if__97, if__162, if__167, if__182, if__187, if__202, \
                         if__207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_9 * ab_x[k] * hf_22[k]
                  - f_33 * ab_x[k] * hf_27[k]
                  + f_11 * ab_x[k] * hf_72[k]
                  - f_10 * ab_x[k] * hf_77[k]
                  - f_168 * ab_x[k] * hf_92[k]
                  + f_167 * ab_x[k] * hf_97[k]
                  + f_9 * ab_x[k] * hf_162[k]
                  - f_33 * ab_x[k] * hf_167[k]
                  - f_168 * ab_x[k] * hf_182[k]
                  + f_167 * ab_x[k] * hf_187[k]
                  + f_170 * ab_x[k] * hf_202[k]
                  - f_169 * ab_x[k] * hf_207[k]
                  + f_9 * if__22[k]
                  - f_33 * if__27[k]
                  + f_11 * if__72[k]
                  - f_10 * if__77[k]
                  - f_168 * if__92[k]
                  + f_167 * if__97[k]
                  + f_9 * if__162[k]
                  - f_33 * if__167[k]
                  - f_168 * if__182[k]
                  + f_167 * if__187[k]
                  + f_170 * if__202[k]
                  - f_169 * if__207[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_20, hf_23, hf_26, hf_70, hf_73, hf_76, hf_90, hf_93, \
                         hf_96, hf_160, hf_163, hf_166, hf_180, hf_183, hf_186, hf_200, \
                         hf_203, hf_206, if__20, if__23, if__46, if__70, if__73, if__90, \
                         if__93, if__116, if__136, if__160, if__163, if__180, if__183, \
                         if__200, if__203, if__226, if__246, if__266 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_192 * ab_x[k] * hf_20[k]
                  - f_15 * ab_x[k] * hf_23[k]
                  + f_192 * ab_y[k] * hf_26[k]
                  + f_193 * ab_x[k] * hf_70[k]
                  - f_17 * ab_x[k] * hf_73[k]
                  + f_193 * ab_y[k] * hf_76[k]
                  - f_194 * ab_x[k] * hf_90[k]
                  + f_18 * ab_x[k] * hf_93[k]
                  - f_194 * ab_y[k] * hf_96[k]
                  + f_192 * ab_x[k] * hf_160[k]
                  - f_15 * ab_x[k] * hf_163[k]
                  + f_192 * ab_y[k] * hf_166[k]
                  - f_194 * ab_x[k] * hf_180[k]
                  + f_18 * ab_x[k] * hf_183[k]
                  - f_194 * ab_y[k] * hf_186[k]
                  + f_195 * ab_x[k] * hf_200[k]
                  - f_196 * ab_x[k] * hf_203[k]
                  + f_195 * ab_y[k] * hf_206[k]
                  + f_192 * if__20[k]
                  - f_15 * if__23[k]
                  + f_192 * if__46[k]
                  + f_193 * if__70[k]
                  - f_17 * if__73[k]
                  - f_194 * if__90[k]
                  + f_18 * if__93[k]
                  + f_193 * if__116[k]
                  - f_194 * if__136[k]
                  + f_192 * if__160[k]
                  - f_15 * if__163[k]
                  - f_194 * if__180[k]
                  + f_18 * if__183[k]
                  + f_195 * if__200[k]
                  - f_196 * if__203[k]
                  + f_192 * if__226[k]
                  - f_194 * if__246[k]
                  + f_195 * if__266[k];
    }

#pragma omp simd aligned(ab_x, hf_1, hf_6, hf_31, hf_36, hf_51, hf_56, hf_101, hf_106, hf_121, \
                         hf_126, hf_141, hf_146, if__1, if__6, if__31, if__36, if__51, if__56, \
                         if__101, if__106, if__121, if__126, if__141, \
                         if__146 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_119 * ab_x[k] * hf_1[k]
                  - f_119 * ab_x[k] * hf_6[k]
                  + f_113 * ab_x[k] * hf_31[k]
                  - f_113 * ab_x[k] * hf_36[k]
                  - f_114 * ab_x[k] * hf_51[k]
                  + f_114 * ab_x[k] * hf_56[k]
                  + f_119 * ab_x[k] * hf_101[k]
                  - f_119 * ab_x[k] * hf_106[k]
                  - f_114 * ab_x[k] * hf_121[k]
                  + f_114 * ab_x[k] * hf_126[k]
                  + f_100 * ab_x[k] * hf_141[k]
                  - f_100 * ab_x[k] * hf_146[k]
                  + f_119 * if__1[k]
                  - f_119 * if__6[k]
                  + f_113 * if__31[k]
                  - f_113 * if__36[k]
                  - f_114 * if__51[k]
                  + f_114 * if__56[k]
                  + f_119 * if__101[k]
                  - f_119 * if__106[k]
                  - f_114 * if__121[k]
                  + f_114 * if__126[k]
                  + f_100 * if__141[k]
                  - f_100 * if__146[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_4, hf_7, hf_34, hf_37, hf_54, hf_57, hf_104, hf_107, \
                         hf_124, hf_127, hf_144, hf_147, if__4, if__17, if__34, if__54, \
                         if__67, if__87, if__104, if__124, if__144, if__157, if__177, \
                         if__197 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_120 * ab_x[k] * hf_4[k]
                  - f_121 * ab_y[k] * hf_7[k]
                  + f_122 * ab_x[k] * hf_34[k]
                  - f_123 * ab_y[k] * hf_37[k]
                  - f_124 * ab_x[k] * hf_54[k]
                  + f_102 * ab_y[k] * hf_57[k]
                  + f_120 * ab_x[k] * hf_104[k]
                  - f_121 * ab_y[k] * hf_107[k]
                  - f_124 * ab_x[k] * hf_124[k]
                  + f_102 * ab_y[k] * hf_127[k]
                  + f_104 * ab_x[k] * hf_144[k]
                  - f_125 * ab_y[k] * hf_147[k]
                  + f_120 * if__4[k]
                  - f_121 * if__17[k]
                  + f_122 * if__34[k]
                  - f_124 * if__54[k]
                  - f_123 * if__67[k]
                  + f_102 * if__87[k]
                  + f_120 * if__104[k]
                  - f_124 * if__124[k]
                  + f_104 * if__144[k]
                  - f_121 * if__157[k]
                  + f_102 * if__177[k]
                  - f_125 * if__197[k];
    }

#pragma omp simd aligned(ab_x, hf_1, hf_6, hf_8, hf_31, hf_36, hf_38, hf_51, hf_56, hf_58, \
                         hf_101, hf_106, hf_108, hf_121, hf_126, hf_128, hf_141, hf_146, \
                         hf_148, if__1, if__6, if__8, if__31, if__36, if__38, if__51, if__56, \
                         if__58, if__101, if__106, if__108, if__121, if__126, if__128, \
                         if__141, if__146, if__148 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = -f_126 * ab_x[k] * hf_1[k]
                  - f_126 * ab_x[k] * hf_6[k]
                  + f_127 * ab_x[k] * hf_8[k]
                  - f_128 * ab_x[k] * hf_31[k]
                  - f_128 * ab_x[k] * hf_36[k]
                  + f_129 * ab_x[k] * hf_38[k]
                  + f_129 * ab_x[k] * hf_51[k]
                  + f_129 * ab_x[k] * hf_56[k]
                  - f_130 * ab_x[k] * hf_58[k]
                  - f_126 * ab_x[k] * hf_101[k]
                  - f_126 * ab_x[k] * hf_106[k]
                  + f_127 * ab_x[k] * hf_108[k]
                  + f_129 * ab_x[k] * hf_121[k]
                  + f_129 * ab_x[k] * hf_126[k]
                  - f_130 * ab_x[k] * hf_128[k]
                  - f_131 * ab_x[k] * hf_141[k]
                  - f_131 * ab_x[k] * hf_146[k]
                  + f_132 * ab_x[k] * hf_148[k]
                  - f_126 * if__1[k]
                  - f_126 * if__6[k]
                  + f_127 * if__8[k]
                  - f_128 * if__31[k]
                  - f_128 * if__36[k]
                  + f_129 * if__38[k]
                  + f_129 * if__51[k]
                  + f_129 * if__56[k]
                  - f_130 * if__58[k]
                  - f_126 * if__101[k]
                  - f_126 * if__106[k]
                  + f_127 * if__108[k]
                  + f_129 * if__121[k]
                  + f_129 * if__126[k]
                  - f_130 * if__128[k]
                  - f_131 * if__141[k]
                  - f_131 * if__146[k]
                  + f_132 * if__148[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_4, hf_7, hf_9, hf_34, hf_37, hf_39, hf_54, hf_57, \
                         hf_59, hf_104, hf_107, hf_109, hf_124, hf_127, hf_129, hf_144, \
                         hf_147, hf_149, if__4, if__17, if__19, if__34, if__54, if__67, \
                         if__69, if__87, if__89, if__104, if__124, if__144, if__157, if__159, \
                         if__177, if__179, if__197, if__199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = -f_133 * ab_x[k] * hf_4[k]
                  - f_133 * ab_y[k] * hf_7[k]
                  + f_134 * ab_y[k] * hf_9[k]
                  - f_135 * ab_x[k] * hf_34[k]
                  - f_135 * ab_y[k] * hf_37[k]
                  + f_136 * ab_y[k] * hf_39[k]
                  + f_137 * ab_x[k] * hf_54[k]
                  + f_137 * ab_y[k] * hf_57[k]
                  - f_138 * ab_y[k] * hf_59[k]
                  - f_133 * ab_x[k] * hf_104[k]
                  - f_133 * ab_y[k] * hf_107[k]
                  + f_134 * ab_y[k] * hf_109[k]
                  + f_137 * ab_x[k] * hf_124[k]
                  + f_137 * ab_y[k] * hf_127[k]
                  - f_138 * ab_y[k] * hf_129[k]
                  - f_139 * ab_x[k] * hf_144[k]
                  - f_139 * ab_y[k] * hf_147[k]
                  + f_140 * ab_y[k] * hf_149[k]
                  - f_133 * if__4[k]
                  - f_133 * if__17[k]
                  + f_134 * if__19[k]
                  - f_135 * if__34[k]
                  + f_137 * if__54[k]
                  - f_135 * if__67[k]
                  + f_136 * if__69[k]
                  + f_137 * if__87[k]
                  - f_138 * if__89[k]
                  - f_133 * if__104[k]
                  + f_137 * if__124[k]
                  - f_139 * if__144[k]
                  - f_133 * if__157[k]
                  + f_134 * if__159[k]
                  + f_137 * if__177[k]
                  - f_138 * if__179[k]
                  - f_139 * if__197[k]
                  + f_140 * if__199[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hf_0, hf_3, hf_5, hf_6, hf_8, hf_9, hf_30, hf_33, \
                         hf_35, hf_36, hf_38, hf_39, hf_50, hf_53, hf_55, hf_56, hf_58, hf_59, \
                         hf_100, hf_103, hf_105, hf_106, hf_108, hf_109, hf_120, hf_123, \
                         hf_125, hf_126, hf_128, hf_129, hf_140, hf_143, hf_145, hf_146, \
                         hf_148, hf_149, if__0, if__3, if__5, if__16, if__18, if__29, if__30, \
                         if__33, if__35, if__50, if__53, if__55, if__66, if__68, if__79, \
                         if__86, if__88, if__99, if__100, if__103, if__105, if__120, if__123, \
                         if__125, if__140, if__143, if__145, if__156, if__158, if__169, \
                         if__176, if__178, if__189, if__196, if__198, \
                         if__209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = f_141 * ab_x[k] * hf_0[k]
                  + f_142 * ab_x[k] * hf_3[k]
                  - f_143 * ab_x[k] * hf_5[k]
                  + f_141 * ab_y[k] * hf_6[k]
                  - f_143 * ab_y[k] * hf_8[k]
                  + f_144 * ab_z[k] * hf_9[k]
                  + f_142 * ab_x[k] * hf_30[k]
                  + f_145 * ab_x[k] * hf_33[k]
                  - f_146 * ab_x[k] * hf_35[k]
                  + f_142 * ab_y[k] * hf_36[k]
                  - f_146 * ab_y[k] * hf_38[k]
                  + f_147 * ab_z[k] * hf_39[k]
                  - f_148 * ab_x[k] * hf_50[k]
                  - f_149 * ab_x[k] * hf_53[k]
                  + f_150 * ab_x[k] * hf_55[k]
                  - f_148 * ab_y[k] * hf_56[k]
                  + f_150 * ab_y[k] * hf_58[k]
                  - f_151 * ab_z[k] * hf_59[k]
                  + f_141 * ab_x[k] * hf_100[k]
                  + f_142 * ab_x[k] * hf_103[k]
                  - f_143 * ab_x[k] * hf_105[k]
                  + f_141 * ab_y[k] * hf_106[k]
                  - f_143 * ab_y[k] * hf_108[k]
                  + f_144 * ab_z[k] * hf_109[k]
                  - f_148 * ab_x[k] * hf_120[k]
                  - f_149 * ab_x[k] * hf_123[k]
                  + f_150 * ab_x[k] * hf_125[k]
                  - f_148 * ab_y[k] * hf_126[k]
                  + f_150 * ab_y[k] * hf_128[k]
                  - f_151 * ab_z[k] * hf_129[k]
                  + f_143 * ab_x[k] * hf_140[k]
                  + f_146 * ab_x[k] * hf_143[k]
                  - f_152 * ab_x[k] * hf_145[k]
                  + f_143 * ab_y[k] * hf_146[k]
                  - f_152 * ab_y[k] * hf_148[k]
                  + f_153 * ab_z[k] * hf_149[k]
                  + f_141 * if__0[k]
                  + f_142 * if__3[k]
                  - f_143 * if__5[k]
                  + f_141 * if__16[k]
                  - f_143 * if__18[k]
                  + f_144 * if__29[k]
                  + f_142 * if__30[k]
                  + f_145 * if__33[k]
                  - f_146 * if__35[k]
                  - f_148 * if__50[k]
                  - f_149 * if__53[k]
                  + f_150 * if__55[k]
                  + f_142 * if__66[k]
                  - f_146 * if__68[k]
                  + f_147 * if__79[k]
                  - f_148 * if__86[k]
                  + f_150 * if__88[k]
                  - f_151 * if__99[k]
                  + f_141 * if__100[k]
                  + f_142 * if__103[k]
                  - f_143 * if__105[k]
                  - f_148 * if__120[k]
                  - f_149 * if__123[k]
                  + f_150 * if__125[k]
                  + f_143 * if__140[k]
                  + f_146 * if__143[k]
                  - f_152 * if__145[k]
                  + f_141 * if__156[k]
                  - f_143 * if__158[k]
                  + f_144 * if__169[k]
                  - f_148 * if__176[k]
                  + f_150 * if__178[k]
                  - f_151 * if__189[k]
                  + f_143 * if__196[k]
                  - f_152 * if__198[k]
                  + f_153 * if__209[k];
    }

#pragma omp simd aligned(ab_x, hf_2, hf_7, hf_9, hf_32, hf_37, hf_39, hf_52, hf_57, hf_59, \
                         hf_102, hf_107, hf_109, hf_122, hf_127, hf_129, hf_142, hf_147, \
                         hf_149, if__2, if__7, if__9, if__32, if__37, if__39, if__52, if__57, \
                         if__59, if__102, if__107, if__109, if__122, if__127, if__129, \
                         if__142, if__147, if__149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_133 * ab_x[k] * hf_2[k]
                  - f_133 * ab_x[k] * hf_7[k]
                  + f_134 * ab_x[k] * hf_9[k]
                  - f_135 * ab_x[k] * hf_32[k]
                  - f_135 * ab_x[k] * hf_37[k]
                  + f_136 * ab_x[k] * hf_39[k]
                  + f_137 * ab_x[k] * hf_52[k]
                  + f_137 * ab_x[k] * hf_57[k]
                  - f_138 * ab_x[k] * hf_59[k]
                  - f_133 * ab_x[k] * hf_102[k]
                  - f_133 * ab_x[k] * hf_107[k]
                  + f_134 * ab_x[k] * hf_109[k]
                  + f_137 * ab_x[k] * hf_122[k]
                  + f_137 * ab_x[k] * hf_127[k]
                  - f_138 * ab_x[k] * hf_129[k]
                  - f_139 * ab_x[k] * hf_142[k]
                  - f_139 * ab_x[k] * hf_147[k]
                  + f_140 * ab_x[k] * hf_149[k]
                  - f_133 * if__2[k]
                  - f_133 * if__7[k]
                  + f_134 * if__9[k]
                  - f_135 * if__32[k]
                  - f_135 * if__37[k]
                  + f_136 * if__39[k]
                  + f_137 * if__52[k]
                  + f_137 * if__57[k]
                  - f_138 * if__59[k]
                  - f_133 * if__102[k]
                  - f_133 * if__107[k]
                  + f_134 * if__109[k]
                  + f_137 * if__122[k]
                  + f_137 * if__127[k]
                  - f_138 * if__129[k]
                  - f_139 * if__142[k]
                  - f_139 * if__147[k]
                  + f_140 * if__149[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_0, hf_5, hf_6, hf_8, hf_30, hf_35, hf_36, hf_38, \
                         hf_50, hf_55, hf_56, hf_58, hf_100, hf_105, hf_106, hf_108, hf_120, \
                         hf_125, hf_126, hf_128, hf_140, hf_145, hf_146, hf_148, if__0, if__5, \
                         if__16, if__18, if__30, if__35, if__50, if__55, if__66, if__68, \
                         if__86, if__88, if__100, if__105, if__120, if__125, if__140, if__145, \
                         if__156, if__158, if__176, if__178, if__196, \
                         if__198 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -f_154 * ab_x[k] * hf_0[k]
                  + f_155 * ab_x[k] * hf_5[k]
                  + f_154 * ab_y[k] * hf_6[k]
                  - f_155 * ab_y[k] * hf_8[k]
                  - f_126 * ab_x[k] * hf_30[k]
                  + f_127 * ab_x[k] * hf_35[k]
                  + f_126 * ab_y[k] * hf_36[k]
                  - f_127 * ab_y[k] * hf_38[k]
                  + f_127 * ab_x[k] * hf_50[k]
                  - f_156 * ab_x[k] * hf_55[k]
                  - f_127 * ab_y[k] * hf_56[k]
                  + f_156 * ab_y[k] * hf_58[k]
                  - f_154 * ab_x[k] * hf_100[k]
                  + f_155 * ab_x[k] * hf_105[k]
                  + f_154 * ab_y[k] * hf_106[k]
                  - f_155 * ab_y[k] * hf_108[k]
                  + f_127 * ab_x[k] * hf_120[k]
                  - f_156 * ab_x[k] * hf_125[k]
                  - f_127 * ab_y[k] * hf_126[k]
                  + f_156 * ab_y[k] * hf_128[k]
                  - f_157 * ab_x[k] * hf_140[k]
                  + f_158 * ab_x[k] * hf_145[k]
                  + f_157 * ab_y[k] * hf_146[k]
                  - f_158 * ab_y[k] * hf_148[k]
                  - f_154 * if__0[k]
                  + f_155 * if__5[k]
                  + f_154 * if__16[k]
                  - f_155 * if__18[k]
                  - f_126 * if__30[k]
                  + f_127 * if__35[k]
                  + f_127 * if__50[k]
                  - f_156 * if__55[k]
                  + f_126 * if__66[k]
                  - f_127 * if__68[k]
                  - f_127 * if__86[k]
                  + f_156 * if__88[k]
                  - f_154 * if__100[k]
                  + f_155 * if__105[k]
                  + f_127 * if__120[k]
                  - f_156 * if__125[k]
                  - f_157 * if__140[k]
                  + f_158 * if__145[k]
                  + f_154 * if__156[k]
                  - f_155 * if__158[k]
                  - f_127 * if__176[k]
                  + f_156 * if__178[k]
                  + f_157 * if__196[k]
                  - f_158 * if__198[k];
    }

#pragma omp simd aligned(ab_x, hf_2, hf_7, hf_32, hf_37, hf_52, hf_57, hf_102, hf_107, hf_122, \
                         hf_127, hf_142, hf_147, if__2, if__7, if__32, if__37, if__52, if__57, \
                         if__102, if__107, if__122, if__127, if__142, \
                         if__147 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_121 * ab_x[k] * hf_2[k]
                  - f_120 * ab_x[k] * hf_7[k]
                  + f_123 * ab_x[k] * hf_32[k]
                  - f_122 * ab_x[k] * hf_37[k]
                  - f_102 * ab_x[k] * hf_52[k]
                  + f_124 * ab_x[k] * hf_57[k]
                  + f_121 * ab_x[k] * hf_102[k]
                  - f_120 * ab_x[k] * hf_107[k]
                  - f_102 * ab_x[k] * hf_122[k]
                  + f_124 * ab_x[k] * hf_127[k]
                  + f_125 * ab_x[k] * hf_142[k]
                  - f_104 * ab_x[k] * hf_147[k]
                  + f_121 * if__2[k]
                  - f_120 * if__7[k]
                  + f_123 * if__32[k]
                  - f_122 * if__37[k]
                  - f_102 * if__52[k]
                  + f_124 * if__57[k]
                  + f_121 * if__102[k]
                  - f_120 * if__107[k]
                  - f_102 * if__122[k]
                  + f_124 * if__127[k]
                  + f_125 * if__142[k]
                  - f_104 * if__147[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_0, hf_3, hf_6, hf_30, hf_33, hf_36, hf_50, hf_53, \
                         hf_56, hf_100, hf_103, hf_106, hf_120, hf_123, hf_126, hf_140, \
                         hf_143, hf_146, if__0, if__3, if__16, if__30, if__33, if__50, if__53, \
                         if__66, if__86, if__100, if__103, if__120, if__123, if__140, if__143, \
                         if__156, if__176, if__196 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_159 * ab_x[k] * hf_0[k]
                  - f_160 * ab_x[k] * hf_3[k]
                  + f_159 * ab_y[k] * hf_6[k]
                  + f_161 * ab_x[k] * hf_30[k]
                  - f_162 * ab_x[k] * hf_33[k]
                  + f_161 * ab_y[k] * hf_36[k]
                  - f_162 * ab_x[k] * hf_50[k]
                  + f_163 * ab_x[k] * hf_53[k]
                  - f_162 * ab_y[k] * hf_56[k]
                  + f_159 * ab_x[k] * hf_100[k]
                  - f_160 * ab_x[k] * hf_103[k]
                  + f_159 * ab_y[k] * hf_106[k]
                  - f_162 * ab_x[k] * hf_120[k]
                  + f_163 * ab_x[k] * hf_123[k]
                  - f_162 * ab_y[k] * hf_126[k]
                  + f_113 * ab_x[k] * hf_140[k]
                  - f_114 * ab_x[k] * hf_143[k]
                  + f_113 * ab_y[k] * hf_146[k]
                  + f_159 * if__0[k]
                  - f_160 * if__3[k]
                  + f_159 * if__16[k]
                  + f_161 * if__30[k]
                  - f_162 * if__33[k]
                  - f_162 * if__50[k]
                  + f_163 * if__53[k]
                  + f_161 * if__66[k]
                  - f_162 * if__86[k]
                  + f_159 * if__100[k]
                  - f_160 * if__103[k]
                  - f_162 * if__120[k]
                  + f_163 * if__123[k]
                  + f_113 * if__140[k]
                  - f_114 * if__143[k]
                  + f_159 * if__156[k]
                  - f_162 * if__176[k]
                  + f_113 * if__196[k];
    }

#pragma omp simd aligned(ab_x, hf_21, hf_26, hf_91, hf_96, hf_161, hf_166, hf_181, hf_186, \
                         if__21, if__26, if__91, if__96, if__161, if__166, if__181, \
                         if__186 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_117 * ab_x[k] * hf_21[k]
                  + f_117 * ab_x[k] * hf_26[k]
                  + f_92 * ab_x[k] * hf_91[k]
                  - f_92 * ab_x[k] * hf_96[k]
                  + f_117 * ab_x[k] * hf_161[k]
                  - f_117 * ab_x[k] * hf_166[k]
                  - f_92 * ab_x[k] * hf_181[k]
                  + f_92 * ab_x[k] * hf_186[k]
                  - f_117 * if__21[k]
                  + f_117 * if__26[k]
                  + f_92 * if__91[k]
                  - f_92 * if__96[k]
                  + f_117 * if__161[k]
                  - f_117 * if__166[k]
                  - f_92 * if__181[k]
                  + f_92 * if__186[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_24, hf_27, hf_94, hf_97, hf_164, hf_167, hf_184, \
                         hf_187, if__24, if__47, if__94, if__137, if__164, if__184, if__227, \
                         if__247 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_197 * ab_x[k] * hf_24[k]
                  + f_198 * ab_y[k] * hf_27[k]
                  + f_94 * ab_x[k] * hf_94[k]
                  - f_95 * ab_y[k] * hf_97[k]
                  + f_197 * ab_x[k] * hf_164[k]
                  - f_198 * ab_y[k] * hf_167[k]
                  - f_94 * ab_x[k] * hf_184[k]
                  + f_95 * ab_y[k] * hf_187[k]
                  - f_197 * if__24[k]
                  + f_198 * if__47[k]
                  + f_94 * if__94[k]
                  - f_95 * if__137[k]
                  + f_197 * if__164[k]
                  - f_94 * if__184[k]
                  - f_198 * if__227[k]
                  + f_95 * if__247[k];
    }

#pragma omp simd aligned(ab_x, hf_21, hf_26, hf_28, hf_91, hf_96, hf_98, hf_161, hf_166, \
                         hf_168, hf_181, hf_186, hf_188, if__21, if__26, if__28, if__91, \
                         if__96, if__98, if__161, if__166, if__168, if__181, if__186, \
                         if__188 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_113 * ab_x[k] * hf_21[k]
                  + f_113 * ab_x[k] * hf_26[k]
                  - f_114 * ab_x[k] * hf_28[k]
                  - f_98 * ab_x[k] * hf_91[k]
                  - f_98 * ab_x[k] * hf_96[k]
                  + f_99 * ab_x[k] * hf_98[k]
                  - f_113 * ab_x[k] * hf_161[k]
                  - f_113 * ab_x[k] * hf_166[k]
                  + f_114 * ab_x[k] * hf_168[k]
                  + f_98 * ab_x[k] * hf_181[k]
                  + f_98 * ab_x[k] * hf_186[k]
                  - f_99 * ab_x[k] * hf_188[k]
                  + f_113 * if__21[k]
                  + f_113 * if__26[k]
                  - f_114 * if__28[k]
                  - f_98 * if__91[k]
                  - f_98 * if__96[k]
                  + f_99 * if__98[k]
                  - f_113 * if__161[k]
                  - f_113 * if__166[k]
                  + f_114 * if__168[k]
                  + f_98 * if__181[k]
                  + f_98 * if__186[k]
                  - f_99 * if__188[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_24, hf_27, hf_29, hf_94, hf_97, hf_99, hf_164, hf_167, \
                         hf_169, hf_184, hf_187, hf_189, if__24, if__47, if__49, if__94, \
                         if__137, if__139, if__164, if__184, if__227, if__229, if__247, \
                         if__249 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = f_122 * ab_x[k] * hf_24[k]
                  + f_122 * ab_y[k] * hf_27[k]
                  - f_125 * ab_y[k] * hf_29[k]
                  - f_102 * ab_x[k] * hf_94[k]
                  - f_102 * ab_y[k] * hf_97[k]
                  + f_103 * ab_y[k] * hf_99[k]
                  - f_122 * ab_x[k] * hf_164[k]
                  - f_122 * ab_y[k] * hf_167[k]
                  + f_125 * ab_y[k] * hf_169[k]
                  + f_102 * ab_x[k] * hf_184[k]
                  + f_102 * ab_y[k] * hf_187[k]
                  - f_103 * ab_y[k] * hf_189[k]
                  + f_122 * if__24[k]
                  + f_122 * if__47[k]
                  - f_125 * if__49[k]
                  - f_102 * if__94[k]
                  - f_102 * if__137[k]
                  + f_103 * if__139[k]
                  - f_122 * if__164[k]
                  + f_102 * if__184[k]
                  - f_122 * if__227[k]
                  + f_125 * if__229[k]
                  + f_102 * if__247[k]
                  - f_103 * if__249[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hf_20, hf_23, hf_25, hf_26, hf_28, hf_29, hf_90, \
                         hf_93, hf_95, hf_96, hf_98, hf_99, hf_160, hf_163, hf_165, hf_166, \
                         hf_168, hf_169, hf_180, hf_183, hf_185, hf_186, hf_188, hf_189, \
                         if__20, if__23, if__25, if__46, if__48, if__59, if__90, if__93, \
                         if__95, if__136, if__138, if__149, if__160, if__163, if__165, \
                         if__180, if__183, if__185, if__226, if__228, if__239, if__246, \
                         if__248, if__259 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_199 * ab_x[k] * hf_20[k]
                  - f_106 * ab_x[k] * hf_23[k]
                  + f_110 * ab_x[k] * hf_25[k]
                  - f_199 * ab_y[k] * hf_26[k]
                  + f_110 * ab_y[k] * hf_28[k]
                  - f_200 * ab_z[k] * hf_29[k]
                  + f_106 * ab_x[k] * hf_90[k]
                  + f_107 * ab_x[k] * hf_93[k]
                  - f_108 * ab_x[k] * hf_95[k]
                  + f_106 * ab_y[k] * hf_96[k]
                  - f_108 * ab_y[k] * hf_98[k]
                  + f_109 * ab_z[k] * hf_99[k]
                  + f_199 * ab_x[k] * hf_160[k]
                  + f_106 * ab_x[k] * hf_163[k]
                  - f_110 * ab_x[k] * hf_165[k]
                  + f_199 * ab_y[k] * hf_166[k]
                  - f_110 * ab_y[k] * hf_168[k]
                  + f_200 * ab_z[k] * hf_169[k]
                  - f_106 * ab_x[k] * hf_180[k]
                  - f_107 * ab_x[k] * hf_183[k]
                  + f_108 * ab_x[k] * hf_185[k]
                  - f_106 * ab_y[k] * hf_186[k]
                  + f_108 * ab_y[k] * hf_188[k]
                  - f_109 * ab_z[k] * hf_189[k]
                  - f_199 * if__20[k]
                  - f_106 * if__23[k]
                  + f_110 * if__25[k]
                  - f_199 * if__46[k]
                  + f_110 * if__48[k]
                  - f_200 * if__59[k]
                  + f_106 * if__90[k]
                  + f_107 * if__93[k]
                  - f_108 * if__95[k]
                  + f_106 * if__136[k]
                  - f_108 * if__138[k]
                  + f_109 * if__149[k]
                  + f_199 * if__160[k]
                  + f_106 * if__163[k]
                  - f_110 * if__165[k]
                  - f_106 * if__180[k]
                  - f_107 * if__183[k]
                  + f_108 * if__185[k]
                  + f_199 * if__226[k]
                  - f_110 * if__228[k]
                  + f_200 * if__239[k]
                  - f_106 * if__246[k]
                  + f_108 * if__248[k]
                  - f_109 * if__259[k];
    }

#pragma omp simd aligned(ab_x, hf_22, hf_27, hf_29, hf_92, hf_97, hf_99, hf_162, hf_167, \
                         hf_169, hf_182, hf_187, hf_189, if__22, if__27, if__29, if__92, \
                         if__97, if__99, if__162, if__167, if__169, if__182, if__187, \
                         if__189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_122 * ab_x[k] * hf_22[k]
                  + f_122 * ab_x[k] * hf_27[k]
                  - f_125 * ab_x[k] * hf_29[k]
                  - f_102 * ab_x[k] * hf_92[k]
                  - f_102 * ab_x[k] * hf_97[k]
                  + f_103 * ab_x[k] * hf_99[k]
                  - f_122 * ab_x[k] * hf_162[k]
                  - f_122 * ab_x[k] * hf_167[k]
                  + f_125 * ab_x[k] * hf_169[k]
                  + f_102 * ab_x[k] * hf_182[k]
                  + f_102 * ab_x[k] * hf_187[k]
                  - f_103 * ab_x[k] * hf_189[k]
                  + f_122 * if__22[k]
                  + f_122 * if__27[k]
                  - f_125 * if__29[k]
                  - f_102 * if__92[k]
                  - f_102 * if__97[k]
                  + f_103 * if__99[k]
                  - f_122 * if__162[k]
                  - f_122 * if__167[k]
                  + f_125 * if__169[k]
                  + f_102 * if__182[k]
                  + f_102 * if__187[k]
                  - f_103 * if__189[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_20, hf_25, hf_26, hf_28, hf_90, hf_95, hf_96, hf_98, \
                         hf_160, hf_165, hf_166, hf_168, hf_180, hf_185, hf_186, hf_188, \
                         if__20, if__25, if__46, if__48, if__90, if__95, if__136, if__138, \
                         if__160, if__165, if__180, if__185, if__226, if__228, if__246, \
                         if__248 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_119 * ab_x[k] * hf_20[k]
                  - f_201 * ab_x[k] * hf_25[k]
                  - f_119 * ab_y[k] * hf_26[k]
                  + f_201 * ab_y[k] * hf_28[k]
                  - f_113 * ab_x[k] * hf_90[k]
                  + f_114 * ab_x[k] * hf_95[k]
                  + f_113 * ab_y[k] * hf_96[k]
                  - f_114 * ab_y[k] * hf_98[k]
                  - f_119 * ab_x[k] * hf_160[k]
                  + f_201 * ab_x[k] * hf_165[k]
                  + f_119 * ab_y[k] * hf_166[k]
                  - f_201 * ab_y[k] * hf_168[k]
                  + f_113 * ab_x[k] * hf_180[k]
                  - f_114 * ab_x[k] * hf_185[k]
                  - f_113 * ab_y[k] * hf_186[k]
                  + f_114 * ab_y[k] * hf_188[k]
                  + f_119 * if__20[k]
                  - f_201 * if__25[k]
                  - f_119 * if__46[k]
                  + f_201 * if__48[k]
                  - f_113 * if__90[k]
                  + f_114 * if__95[k]
                  + f_113 * if__136[k]
                  - f_114 * if__138[k]
                  - f_119 * if__160[k]
                  + f_201 * if__165[k]
                  + f_113 * if__180[k]
                  - f_114 * if__185[k]
                  + f_119 * if__226[k]
                  - f_201 * if__228[k]
                  - f_113 * if__246[k]
                  + f_114 * if__248[k];
    }

#pragma omp simd aligned(ab_x, hf_22, hf_27, hf_92, hf_97, hf_162, hf_167, hf_182, hf_187, \
                         if__22, if__27, if__92, if__97, if__162, if__167, if__182, \
                         if__187 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_198 * ab_x[k] * hf_22[k]
                  + f_197 * ab_x[k] * hf_27[k]
                  + f_95 * ab_x[k] * hf_92[k]
                  - f_94 * ab_x[k] * hf_97[k]
                  + f_198 * ab_x[k] * hf_162[k]
                  - f_197 * ab_x[k] * hf_167[k]
                  - f_95 * ab_x[k] * hf_182[k]
                  + f_94 * ab_x[k] * hf_187[k]
                  - f_198 * if__22[k]
                  + f_197 * if__27[k]
                  + f_95 * if__92[k]
                  - f_94 * if__97[k]
                  + f_198 * if__162[k]
                  - f_197 * if__167[k]
                  - f_95 * if__182[k]
                  + f_94 * if__187[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_20, hf_23, hf_26, hf_90, hf_93, hf_96, hf_160, hf_163, \
                         hf_166, hf_180, hf_183, hf_186, if__20, if__23, if__46, if__90, \
                         if__93, if__136, if__160, if__163, if__180, if__183, if__226, \
                         if__246 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_202 * ab_x[k] * hf_20[k]
                  + f_203 * ab_x[k] * hf_23[k]
                  - f_202 * ab_y[k] * hf_26[k]
                  + f_115 * ab_x[k] * hf_90[k]
                  - f_116 * ab_x[k] * hf_93[k]
                  + f_115 * ab_y[k] * hf_96[k]
                  + f_202 * ab_x[k] * hf_160[k]
                  - f_203 * ab_x[k] * hf_163[k]
                  + f_202 * ab_y[k] * hf_166[k]
                  - f_115 * ab_x[k] * hf_180[k]
                  + f_116 * ab_x[k] * hf_183[k]
                  - f_115 * ab_y[k] * hf_186[k]
                  - f_202 * if__20[k]
                  + f_203 * if__23[k]
                  - f_202 * if__46[k]
                  + f_115 * if__90[k]
                  - f_116 * if__93[k]
                  + f_115 * if__136[k]
                  + f_202 * if__160[k]
                  - f_203 * if__163[k]
                  - f_115 * if__180[k]
                  + f_116 * if__183[k]
                  + f_202 * if__226[k]
                  - f_115 * if__246[k];
    }

#pragma omp simd aligned(ab_x, hf_1, hf_6, hf_31, hf_36, hf_51, hf_56, hf_101, hf_106, hf_121, \
                         hf_126, if__1, if__6, if__31, if__36, if__51, if__56, if__101, \
                         if__106, if__121, if__126 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -f_56 * ab_x[k] * hf_1[k]
                  + f_56 * ab_x[k] * hf_6[k]
                  + f_54 * ab_x[k] * hf_31[k]
                  - f_54 * ab_x[k] * hf_36[k]
                  + f_57 * ab_x[k] * hf_51[k]
                  - f_57 * ab_x[k] * hf_56[k]
                  + f_53 * ab_x[k] * hf_101[k]
                  - f_53 * ab_x[k] * hf_106[k]
                  - f_55 * ab_x[k] * hf_121[k]
                  + f_55 * ab_x[k] * hf_126[k]
                  - f_56 * if__1[k]
                  + f_56 * if__6[k]
                  + f_54 * if__31[k]
                  - f_54 * if__36[k]
                  + f_57 * if__51[k]
                  - f_57 * if__56[k]
                  + f_53 * if__101[k]
                  - f_53 * if__106[k]
                  - f_55 * if__121[k]
                  + f_55 * if__126[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_4, hf_7, hf_34, hf_37, hf_54, hf_57, hf_104, hf_107, \
                         hf_124, hf_127, if__4, if__17, if__34, if__54, if__67, if__87, \
                         if__104, if__124, if__157, if__177 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = -3.28125 * ab_x[k] * hf_4[k]
                  + 1.09375 * ab_y[k] * hf_7[k]
                  + 6.5625 * ab_x[k] * hf_34[k]
                  - 2.1875 * ab_y[k] * hf_37[k]
                  + 26.25 * ab_x[k] * hf_54[k]
                  - 8.75 * ab_y[k] * hf_57[k]
                  + 9.84375 * ab_x[k] * hf_104[k]
                  - 3.28125 * ab_y[k] * hf_107[k]
                  - 78.75 * ab_x[k] * hf_124[k]
                  + 26.25 * ab_y[k] * hf_127[k]
                  - 3.28125 * if__4[k]
                  + 1.09375 * if__17[k]
                  + 6.5625 * if__34[k]
                  + 26.25 * if__54[k]
                  - 2.1875 * if__67[k]
                  - 8.75 * if__87[k]
                  + 9.84375 * if__104[k]
                  - 78.75 * if__124[k]
                  - 3.28125 * if__157[k]
                  + 26.25 * if__177[k];
    }

#pragma omp simd aligned(ab_x, hf_1, hf_6, hf_8, hf_31, hf_36, hf_38, hf_51, hf_56, hf_58, \
                         hf_101, hf_106, hf_108, hf_121, hf_126, hf_128, if__1, if__6, if__8, \
                         if__31, if__36, if__38, if__51, if__56, if__58, if__101, if__106, \
                         if__108, if__121, if__126, if__128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_62 * ab_x[k] * hf_1[k]
                  + f_62 * ab_x[k] * hf_6[k]
                  - f_24 * ab_x[k] * hf_8[k]
                  - f_59 * ab_x[k] * hf_31[k]
                  - f_59 * ab_x[k] * hf_36[k]
                  + f_27 * ab_x[k] * hf_38[k]
                  - f_63 * ab_x[k] * hf_51[k]
                  - f_63 * ab_x[k] * hf_56[k]
                  + f_46 * ab_x[k] * hf_58[k]
                  - f_58 * ab_x[k] * hf_101[k]
                  - f_58 * ab_x[k] * hf_106[k]
                  + f_23 * ab_x[k] * hf_108[k]
                  + f_60 * ab_x[k] * hf_121[k]
                  + f_60 * ab_x[k] * hf_126[k]
                  - f_61 * ab_x[k] * hf_128[k]
                  + f_62 * if__1[k]
                  + f_62 * if__6[k]
                  - f_24 * if__8[k]
                  - f_59 * if__31[k]
                  - f_59 * if__36[k]
                  + f_27 * if__38[k]
                  - f_63 * if__51[k]
                  - f_63 * if__56[k]
                  + f_46 * if__58[k]
                  - f_58 * if__101[k]
                  - f_58 * if__106[k]
                  + f_23 * if__108[k]
                  + f_60 * if__121[k]
                  + f_60 * if__126[k]
                  - f_61 * if__128[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_4, hf_7, hf_9, hf_34, hf_37, hf_39, hf_54, hf_57, \
                         hf_59, hf_104, hf_107, hf_109, hf_124, hf_127, hf_129, if__4, if__17, \
                         if__19, if__34, if__54, if__67, if__69, if__87, if__89, if__104, \
                         if__124, if__157, if__159, if__177, if__179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = f_68 * ab_x[k] * hf_4[k]
                  + f_68 * ab_y[k] * hf_7[k]
                  - f_69 * ab_y[k] * hf_9[k]
                  - f_65 * ab_x[k] * hf_34[k]
                  - f_65 * ab_y[k] * hf_37[k]
                  + f_66 * ab_y[k] * hf_39[k]
                  - f_44 * ab_x[k] * hf_54[k]
                  - f_44 * ab_y[k] * hf_57[k]
                  + f_70 * ab_y[k] * hf_59[k]
                  - f_64 * ab_x[k] * hf_104[k]
                  - f_64 * ab_y[k] * hf_107[k]
                  + f_51 * ab_y[k] * hf_109[k]
                  + f_52 * ab_x[k] * hf_124[k]
                  + f_52 * ab_y[k] * hf_127[k]
                  - f_67 * ab_y[k] * hf_129[k]
                  + f_68 * if__4[k]
                  + f_68 * if__17[k]
                  - f_69 * if__19[k]
                  - f_65 * if__34[k]
                  - f_44 * if__54[k]
                  - f_65 * if__67[k]
                  + f_66 * if__69[k]
                  - f_44 * if__87[k]
                  + f_70 * if__89[k]
                  - f_64 * if__104[k]
                  + f_52 * if__124[k]
                  - f_64 * if__157[k]
                  + f_51 * if__159[k]
                  + f_52 * if__177[k]
                  - f_67 * if__179[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hf_0, hf_3, hf_5, hf_6, hf_8, hf_9, hf_30, hf_33, \
                         hf_35, hf_36, hf_38, hf_39, hf_50, hf_53, hf_55, hf_56, hf_58, hf_59, \
                         hf_100, hf_103, hf_105, hf_106, hf_108, hf_109, hf_120, hf_123, \
                         hf_125, hf_126, hf_128, hf_129, if__0, if__3, if__5, if__16, if__18, \
                         if__29, if__30, if__33, if__35, if__50, if__53, if__55, if__66, \
                         if__68, if__79, if__86, if__88, if__99, if__100, if__103, if__105, \
                         if__120, if__123, if__125, if__156, if__158, if__169, if__176, \
                         if__178, if__189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_79 * ab_x[k] * hf_0[k]
                  - f_34 * ab_x[k] * hf_3[k]
                  + f_73 * ab_x[k] * hf_5[k]
                  - f_79 * ab_y[k] * hf_6[k]
                  + f_73 * ab_y[k] * hf_8[k]
                  - f_80 * ab_z[k] * hf_9[k]
                  + f_34 * ab_x[k] * hf_30[k]
                  + f_13 * ab_x[k] * hf_33[k]
                  - f_74 * ab_x[k] * hf_35[k]
                  + f_34 * ab_y[k] * hf_36[k]
                  - f_74 * ab_y[k] * hf_38[k]
                  + f_75 * ab_z[k] * hf_39[k]
                  + f_73 * ab_x[k] * hf_50[k]
                  + f_74 * ab_x[k] * hf_53[k]
                  - f_78 * ab_x[k] * hf_55[k]
                  + f_73 * ab_y[k] * hf_56[k]
                  - f_78 * ab_y[k] * hf_58[k]
                  + f_81 * ab_z[k] * hf_59[k]
                  + f_71 * ab_x[k] * hf_100[k]
                  + f_72 * ab_x[k] * hf_103[k]
                  - f_14 * ab_x[k] * hf_105[k]
                  + f_71 * ab_y[k] * hf_106[k]
                  - f_14 * ab_y[k] * hf_108[k]
                  + f_73 * ab_z[k] * hf_109[k]
                  - f_14 * ab_x[k] * hf_120[k]
                  - f_76 * ab_x[k] * hf_123[k]
                  + f_77 * ab_x[k] * hf_125[k]
                  - f_14 * ab_y[k] * hf_126[k]
                  + f_77 * ab_y[k] * hf_128[k]
                  - f_78 * ab_z[k] * hf_129[k]
                  - f_79 * if__0[k]
                  - f_34 * if__3[k]
                  + f_73 * if__5[k]
                  - f_79 * if__16[k]
                  + f_73 * if__18[k]
                  - f_80 * if__29[k]
                  + f_34 * if__30[k]
                  + f_13 * if__33[k]
                  - f_74 * if__35[k]
                  + f_73 * if__50[k]
                  + f_74 * if__53[k]
                  - f_78 * if__55[k]
                  + f_34 * if__66[k]
                  - f_74 * if__68[k]
                  + f_75 * if__79[k]
                  + f_73 * if__86[k]
                  - f_78 * if__88[k]
                  + f_81 * if__99[k]
                  + f_71 * if__100[k]
                  + f_72 * if__103[k]
                  - f_14 * if__105[k]
                  - f_14 * if__120[k]
                  - f_76 * if__123[k]
                  + f_77 * if__125[k]
                  + f_71 * if__156[k]
                  - f_14 * if__158[k]
                  + f_73 * if__169[k]
                  - f_14 * if__176[k]
                  + f_77 * if__178[k]
                  - f_78 * if__189[k];
    }

#pragma omp simd aligned(ab_x, hf_2, hf_7, hf_9, hf_32, hf_37, hf_39, hf_52, hf_57, hf_59, \
                         hf_102, hf_107, hf_109, hf_122, hf_127, hf_129, if__2, if__7, if__9, \
                         if__32, if__37, if__39, if__52, if__57, if__59, if__102, if__107, \
                         if__109, if__122, if__127, if__129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = f_68 * ab_x[k] * hf_2[k]
                  + f_68 * ab_x[k] * hf_7[k]
                  - f_69 * ab_x[k] * hf_9[k]
                  - f_65 * ab_x[k] * hf_32[k]
                  - f_65 * ab_x[k] * hf_37[k]
                  + f_66 * ab_x[k] * hf_39[k]
                  - f_44 * ab_x[k] * hf_52[k]
                  - f_44 * ab_x[k] * hf_57[k]
                  + f_70 * ab_x[k] * hf_59[k]
                  - f_64 * ab_x[k] * hf_102[k]
                  - f_64 * ab_x[k] * hf_107[k]
                  + f_51 * ab_x[k] * hf_109[k]
                  + f_52 * ab_x[k] * hf_122[k]
                  + f_52 * ab_x[k] * hf_127[k]
                  - f_67 * ab_x[k] * hf_129[k]
                  + f_68 * if__2[k]
                  + f_68 * if__7[k]
                  - f_69 * if__9[k]
                  - f_65 * if__32[k]
                  - f_65 * if__37[k]
                  + f_66 * if__39[k]
                  - f_44 * if__52[k]
                  - f_44 * if__57[k]
                  + f_70 * if__59[k]
                  - f_64 * if__102[k]
                  - f_64 * if__107[k]
                  + f_51 * if__109[k]
                  + f_52 * if__122[k]
                  + f_52 * if__127[k]
                  - f_67 * if__129[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_0, hf_5, hf_6, hf_8, hf_30, hf_35, hf_36, hf_38, \
                         hf_50, hf_55, hf_56, hf_58, hf_100, hf_105, hf_106, hf_108, hf_120, \
                         hf_125, hf_126, hf_128, if__0, if__5, if__16, if__18, if__30, if__35, \
                         if__50, if__55, if__66, if__68, if__86, if__88, if__100, if__105, \
                         if__120, if__125, if__156, if__158, if__176, \
                         if__178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = f_84 * ab_x[k] * hf_0[k]
                  - f_58 * ab_x[k] * hf_5[k]
                  - f_84 * ab_y[k] * hf_6[k]
                  + f_58 * ab_y[k] * hf_8[k]
                  - f_62 * ab_x[k] * hf_30[k]
                  + f_24 * ab_x[k] * hf_35[k]
                  + f_62 * ab_y[k] * hf_36[k]
                  - f_24 * ab_y[k] * hf_38[k]
                  - f_85 * ab_x[k] * hf_50[k]
                  + f_60 * ab_x[k] * hf_55[k]
                  + f_85 * ab_y[k] * hf_56[k]
                  - f_60 * ab_y[k] * hf_58[k]
                  - f_82 * ab_x[k] * hf_100[k]
                  + f_25 * ab_x[k] * hf_105[k]
                  + f_82 * ab_y[k] * hf_106[k]
                  - f_25 * ab_y[k] * hf_108[k]
                  + f_27 * ab_x[k] * hf_120[k]
                  - f_83 * ab_x[k] * hf_125[k]
                  - f_27 * ab_y[k] * hf_126[k]
                  + f_83 * ab_y[k] * hf_128[k]
                  + f_84 * if__0[k]
                  - f_58 * if__5[k]
                  - f_84 * if__16[k]
                  + f_58 * if__18[k]
                  - f_62 * if__30[k]
                  + f_24 * if__35[k]
                  - f_85 * if__50[k]
                  + f_60 * if__55[k]
                  + f_62 * if__66[k]
                  - f_24 * if__68[k]
                  + f_85 * if__86[k]
                  - f_60 * if__88[k]
                  - f_82 * if__100[k]
                  + f_25 * if__105[k]
                  + f_27 * if__120[k]
                  - f_83 * if__125[k]
                  + f_82 * if__156[k]
                  - f_25 * if__158[k]
                  - f_27 * if__176[k]
                  + f_83 * if__178[k];
    }

#pragma omp simd aligned(ab_x, hf_2, hf_7, hf_32, hf_37, hf_52, hf_57, hf_102, hf_107, hf_122, \
                         hf_127, if__2, if__7, if__32, if__37, if__52, if__57, if__102, \
                         if__107, if__122, if__127 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = -1.09375 * ab_x[k] * hf_2[k]
                  + 3.28125 * ab_x[k] * hf_7[k]
                  + 2.1875 * ab_x[k] * hf_32[k]
                  - 6.5625 * ab_x[k] * hf_37[k]
                  + 8.75 * ab_x[k] * hf_52[k]
                  - 26.25 * ab_x[k] * hf_57[k]
                  + 3.28125 * ab_x[k] * hf_102[k]
                  - 9.84375 * ab_x[k] * hf_107[k]
                  - 26.25 * ab_x[k] * hf_122[k]
                  + 78.75 * ab_x[k] * hf_127[k]
                  - 1.09375 * if__2[k]
                  + 3.28125 * if__7[k]
                  + 2.1875 * if__32[k]
                  - 6.5625 * if__37[k]
                  + 8.75 * if__52[k]
                  - 26.25 * if__57[k]
                  + 3.28125 * if__102[k]
                  - 9.84375 * if__107[k]
                  - 26.25 * if__122[k]
                  + 78.75 * if__127[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_0, hf_3, hf_6, hf_30, hf_33, hf_36, hf_50, hf_53, \
                         hf_56, hf_100, hf_103, hf_106, hf_120, hf_123, hf_126, if__0, if__3, \
                         if__16, if__30, if__33, if__50, if__53, if__66, if__86, if__100, \
                         if__103, if__120, if__123, if__156, if__176 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = -f_90 * ab_x[k] * hf_0[k]
                  + f_91 * ab_x[k] * hf_3[k]
                  - f_90 * ab_y[k] * hf_6[k]
                  + f_88 * ab_x[k] * hf_30[k]
                  - f_53 * ab_x[k] * hf_33[k]
                  + f_88 * ab_y[k] * hf_36[k]
                  + f_54 * ab_x[k] * hf_50[k]
                  - f_43 * ab_x[k] * hf_53[k]
                  + f_54 * ab_y[k] * hf_56[k]
                  + f_86 * ab_x[k] * hf_100[k]
                  - f_87 * ab_x[k] * hf_103[k]
                  + f_86 * ab_y[k] * hf_106[k]
                  - f_89 * ab_x[k] * hf_120[k]
                  + f_42 * ab_x[k] * hf_123[k]
                  - f_89 * ab_y[k] * hf_126[k]
                  - f_90 * if__0[k]
                  + f_91 * if__3[k]
                  - f_90 * if__16[k]
                  + f_88 * if__30[k]
                  - f_53 * if__33[k]
                  + f_54 * if__50[k]
                  - f_43 * if__53[k]
                  + f_88 * if__66[k]
                  + f_54 * if__86[k]
                  + f_86 * if__100[k]
                  - f_87 * if__103[k]
                  - f_89 * if__120[k]
                  + f_42 * if__123[k]
                  + f_86 * if__156[k]
                  - f_89 * if__176[k];
    }

#pragma omp simd aligned(ab_x, hf_21, hf_26, hf_71, hf_76, hf_161, hf_166, if__21, if__26, \
                         if__71, if__76, if__161, if__166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = 6.5625 * ab_x[k] * hf_21[k]
                  - 6.5625 * ab_x[k] * hf_26[k]
                  - 39.375 * ab_x[k] * hf_71[k]
                  + 39.375 * ab_x[k] * hf_76[k]
                  + 6.5625 * ab_x[k] * hf_161[k]
                  - 6.5625 * ab_x[k] * hf_166[k]
                  + 6.5625 * if__21[k]
                  - 6.5625 * if__26[k]
                  - 39.375 * if__71[k]
                  + 39.375 * if__76[k]
                  + 6.5625 * if__161[k]
                  - 6.5625 * if__166[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_24, hf_27, hf_74, hf_77, hf_164, hf_167, if__24, \
                         if__47, if__74, if__117, if__164, if__227 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = f_204 * ab_x[k] * hf_24[k]
                  - f_53 * ab_y[k] * hf_27[k]
                  - f_205 * ab_x[k] * hf_74[k]
                  + f_206 * ab_y[k] * hf_77[k]
                  + f_204 * ab_x[k] * hf_164[k]
                  - f_53 * ab_y[k] * hf_167[k]
                  + f_204 * if__24[k]
                  - f_53 * if__47[k]
                  - f_205 * if__74[k]
                  + f_206 * if__117[k]
                  + f_204 * if__164[k]
                  - f_53 * if__227[k];
    }

#pragma omp simd aligned(ab_x, hf_21, hf_26, hf_28, hf_71, hf_76, hf_78, hf_161, hf_166, \
                         hf_168, if__21, if__26, if__28, if__71, if__76, if__78, if__161, \
                         if__166, if__168 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = -f_65 * ab_x[k] * hf_21[k]
                  - f_65 * ab_x[k] * hf_26[k]
                  + f_207 * ab_x[k] * hf_28[k]
                  + f_207 * ab_x[k] * hf_71[k]
                  + f_207 * ab_x[k] * hf_76[k]
                  - f_208 * ab_x[k] * hf_78[k]
                  - f_65 * ab_x[k] * hf_161[k]
                  - f_65 * ab_x[k] * hf_166[k]
                  + f_207 * ab_x[k] * hf_168[k]
                  - f_65 * if__21[k]
                  - f_65 * if__26[k]
                  + f_207 * if__28[k]
                  + f_207 * if__71[k]
                  + f_207 * if__76[k]
                  - f_208 * if__78[k]
                  - f_65 * if__161[k]
                  - f_65 * if__166[k]
                  + f_207 * if__168[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_24, hf_27, hf_29, hf_74, hf_77, hf_79, hf_164, hf_167, \
                         hf_169, if__24, if__47, if__49, if__74, if__117, if__119, if__164, \
                         if__227, if__229 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = -f_25 * ab_x[k] * hf_24[k]
                  - f_25 * ab_y[k] * hf_27[k]
                  + f_27 * ab_y[k] * hf_29[k]
                  + f_209 * ab_x[k] * hf_74[k]
                  + f_209 * ab_y[k] * hf_77[k]
                  - f_83 * ab_y[k] * hf_79[k]
                  - f_25 * ab_x[k] * hf_164[k]
                  - f_25 * ab_y[k] * hf_167[k]
                  + f_27 * ab_y[k] * hf_169[k]
                  - f_25 * if__24[k]
                  - f_25 * if__47[k]
                  + f_27 * if__49[k]
                  + f_209 * if__74[k]
                  + f_209 * if__117[k]
                  - f_83 * if__119[k]
                  - f_25 * if__164[k]
                  - f_25 * if__227[k]
                  + f_27 * if__229[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hf_20, hf_23, hf_25, hf_26, hf_28, hf_29, hf_70, \
                         hf_73, hf_75, hf_76, hf_78, hf_79, hf_160, hf_163, hf_165, hf_166, \
                         hf_168, hf_169, if__20, if__23, if__25, if__46, if__48, if__59, \
                         if__70, if__73, if__75, if__116, if__118, if__129, if__160, if__163, \
                         if__165, if__226, if__228, if__239 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = f_210 * ab_x[k] * hf_20[k]
                  + f_19 * ab_x[k] * hf_23[k]
                  - f_48 * ab_x[k] * hf_25[k]
                  + f_210 * ab_y[k] * hf_26[k]
                  - f_48 * ab_y[k] * hf_28[k]
                  + f_20 * ab_z[k] * hf_29[k]
                  - f_211 * ab_x[k] * hf_70[k]
                  - f_212 * ab_x[k] * hf_73[k]
                  + f_213 * ab_x[k] * hf_75[k]
                  - f_211 * ab_y[k] * hf_76[k]
                  + f_213 * ab_y[k] * hf_78[k]
                  - f_214 * ab_z[k] * hf_79[k]
                  + f_210 * ab_x[k] * hf_160[k]
                  + f_19 * ab_x[k] * hf_163[k]
                  - f_48 * ab_x[k] * hf_165[k]
                  + f_210 * ab_y[k] * hf_166[k]
                  - f_48 * ab_y[k] * hf_168[k]
                  + f_20 * ab_z[k] * hf_169[k]
                  + f_210 * if__20[k]
                  + f_19 * if__23[k]
                  - f_48 * if__25[k]
                  + f_210 * if__46[k]
                  - f_48 * if__48[k]
                  + f_20 * if__59[k]
                  - f_211 * if__70[k]
                  - f_212 * if__73[k]
                  + f_213 * if__75[k]
                  - f_211 * if__116[k]
                  + f_213 * if__118[k]
                  - f_214 * if__129[k]
                  + f_210 * if__160[k]
                  + f_19 * if__163[k]
                  - f_48 * if__165[k]
                  + f_210 * if__226[k]
                  - f_48 * if__228[k]
                  + f_20 * if__239[k];
    }

#pragma omp simd aligned(ab_x, hf_22, hf_27, hf_29, hf_72, hf_77, hf_79, hf_162, hf_167, \
                         hf_169, if__22, if__27, if__29, if__72, if__77, if__79, if__162, \
                         if__167, if__169 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_25 * ab_x[k] * hf_22[k]
                  - f_25 * ab_x[k] * hf_27[k]
                  + f_27 * ab_x[k] * hf_29[k]
                  + f_209 * ab_x[k] * hf_72[k]
                  + f_209 * ab_x[k] * hf_77[k]
                  - f_83 * ab_x[k] * hf_79[k]
                  - f_25 * ab_x[k] * hf_162[k]
                  - f_25 * ab_x[k] * hf_167[k]
                  + f_27 * ab_x[k] * hf_169[k]
                  - f_25 * if__22[k]
                  - f_25 * if__27[k]
                  + f_27 * if__29[k]
                  + f_209 * if__72[k]
                  + f_209 * if__77[k]
                  - f_83 * if__79[k]
                  - f_25 * if__162[k]
                  - f_25 * if__167[k]
                  + f_27 * if__169[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_20, hf_25, hf_26, hf_28, hf_70, hf_75, hf_76, hf_78, \
                         hf_160, hf_165, hf_166, hf_168, if__20, if__25, if__46, if__48, \
                         if__70, if__75, if__116, if__118, if__160, if__165, if__226, \
                         if__228 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = -f_68 * ab_x[k] * hf_20[k]
                  + f_215 * ab_x[k] * hf_25[k]
                  + f_68 * ab_y[k] * hf_26[k]
                  - f_215 * ab_y[k] * hf_28[k]
                  + f_215 * ab_x[k] * hf_70[k]
                  - f_216 * ab_x[k] * hf_75[k]
                  - f_215 * ab_y[k] * hf_76[k]
                  + f_216 * ab_y[k] * hf_78[k]
                  - f_68 * ab_x[k] * hf_160[k]
                  + f_215 * ab_x[k] * hf_165[k]
                  + f_68 * ab_y[k] * hf_166[k]
                  - f_215 * ab_y[k] * hf_168[k]
                  - f_68 * if__20[k]
                  + f_215 * if__25[k]
                  + f_68 * if__46[k]
                  - f_215 * if__48[k]
                  + f_215 * if__70[k]
                  - f_216 * if__75[k]
                  - f_215 * if__116[k]
                  + f_216 * if__118[k]
                  - f_68 * if__160[k]
                  + f_215 * if__165[k]
                  + f_68 * if__226[k]
                  - f_215 * if__228[k];
    }

#pragma omp simd aligned(ab_x, hf_22, hf_27, hf_72, hf_77, hf_162, hf_167, if__22, if__27, \
                         if__72, if__77, if__162, if__167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = f_53 * ab_x[k] * hf_22[k]
                  - f_204 * ab_x[k] * hf_27[k]
                  - f_206 * ab_x[k] * hf_72[k]
                  + f_205 * ab_x[k] * hf_77[k]
                  + f_53 * ab_x[k] * hf_162[k]
                  - f_204 * ab_x[k] * hf_167[k]
                  + f_53 * if__22[k]
                  - f_204 * if__27[k]
                  - f_206 * if__72[k]
                  + f_205 * if__77[k]
                  + f_53 * if__162[k]
                  - f_204 * if__167[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_20, hf_23, hf_26, hf_70, hf_73, hf_76, hf_160, hf_163, \
                         hf_166, if__20, if__23, if__46, if__70, if__73, if__116, if__160, \
                         if__163, if__226 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = 1.640625 * ab_x[k] * hf_20[k]
                  - 9.84375 * ab_x[k] * hf_23[k]
                  + 1.640625 * ab_y[k] * hf_26[k]
                  - 9.84375 * ab_x[k] * hf_70[k]
                  + 59.0625 * ab_x[k] * hf_73[k]
                  - 9.84375 * ab_y[k] * hf_76[k]
                  + 1.640625 * ab_x[k] * hf_160[k]
                  - 9.84375 * ab_x[k] * hf_163[k]
                  + 1.640625 * ab_y[k] * hf_166[k]
                  + 1.640625 * if__20[k]
                  - 9.84375 * if__23[k]
                  + 1.640625 * if__46[k]
                  - 9.84375 * if__70[k]
                  + 59.0625 * if__73[k]
                  - 9.84375 * if__116[k]
                  + 1.640625 * if__160[k]
                  - 9.84375 * if__163[k]
                  + 1.640625 * if__226[k];
    }

#pragma omp simd aligned(ab_x, hf_1, hf_6, hf_31, hf_36, hf_101, hf_106, if__1, if__6, if__31, \
                         if__36, if__101, if__106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = f_2 * ab_x[k] * hf_1[k]
                  - f_2 * ab_x[k] * hf_6[k]
                  - f_1 * ab_x[k] * hf_31[k]
                  + f_1 * ab_x[k] * hf_36[k]
                  + f_0 * ab_x[k] * hf_101[k]
                  - f_0 * ab_x[k] * hf_106[k]
                  + f_2 * if__1[k]
                  - f_2 * if__6[k]
                  - f_1 * if__31[k]
                  + f_1 * if__36[k]
                  + f_0 * if__101[k]
                  - f_0 * if__106[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_4, hf_7, hf_34, hf_37, hf_104, hf_107, if__4, if__17, \
                         if__34, if__67, if__104, if__157 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = f_7 * ab_x[k] * hf_4[k]
                  - f_8 * ab_y[k] * hf_7[k]
                  - f_5 * ab_x[k] * hf_34[k]
                  + f_6 * ab_y[k] * hf_37[k]
                  + f_3 * ab_x[k] * hf_104[k]
                  - f_4 * ab_y[k] * hf_107[k]
                  + f_7 * if__4[k]
                  - f_8 * if__17[k]
                  - f_5 * if__34[k]
                  + f_6 * if__67[k]
                  + f_3 * if__104[k]
                  - f_4 * if__157[k];
    }

#pragma omp simd aligned(ab_x, hf_1, hf_6, hf_8, hf_31, hf_36, hf_38, hf_101, hf_106, hf_108, \
                         if__1, if__6, if__8, if__31, if__36, if__38, if__101, if__106, \
                         if__108 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = -f_13 * ab_x[k] * hf_1[k]
                  - f_13 * ab_x[k] * hf_6[k]
                  + f_14 * ab_x[k] * hf_8[k]
                  + f_11 * ab_x[k] * hf_31[k]
                  + f_11 * ab_x[k] * hf_36[k]
                  - f_12 * ab_x[k] * hf_38[k]
                  - f_9 * ab_x[k] * hf_101[k]
                  - f_9 * ab_x[k] * hf_106[k]
                  + f_10 * ab_x[k] * hf_108[k]
                  - f_13 * if__1[k]
                  - f_13 * if__6[k]
                  + f_14 * if__8[k]
                  + f_11 * if__31[k]
                  + f_11 * if__36[k]
                  - f_12 * if__38[k]
                  - f_9 * if__101[k]
                  - f_9 * if__106[k]
                  + f_10 * if__108[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_4, hf_7, hf_9, hf_34, hf_37, hf_39, hf_104, hf_107, \
                         hf_109, if__4, if__17, if__19, if__34, if__67, if__69, if__104, \
                         if__157, if__159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = -f_19 * ab_x[k] * hf_4[k]
                  - f_19 * ab_y[k] * hf_7[k]
                  + f_20 * ab_y[k] * hf_9[k]
                  + f_17 * ab_x[k] * hf_34[k]
                  + f_17 * ab_y[k] * hf_37[k]
                  - f_18 * ab_y[k] * hf_39[k]
                  - f_15 * ab_x[k] * hf_104[k]
                  - f_15 * ab_y[k] * hf_107[k]
                  + f_16 * ab_y[k] * hf_109[k]
                  - f_19 * if__4[k]
                  - f_19 * if__17[k]
                  + f_20 * if__19[k]
                  + f_17 * if__34[k]
                  + f_17 * if__67[k]
                  - f_18 * if__69[k]
                  - f_15 * if__104[k]
                  - f_15 * if__157[k]
                  + f_16 * if__159[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hf_0, hf_3, hf_5, hf_6, hf_8, hf_9, hf_30, hf_33, \
                         hf_35, hf_36, hf_38, hf_39, hf_100, hf_103, hf_105, hf_106, hf_108, \
                         hf_109, if__0, if__3, if__5, if__16, if__18, if__29, if__30, if__33, \
                         if__35, if__66, if__68, if__79, if__100, if__103, if__105, if__156, \
                         if__158, if__169 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = f_28 * ab_x[k] * hf_0[k]
                  + f_29 * ab_x[k] * hf_3[k]
                  - f_30 * ab_x[k] * hf_5[k]
                  + f_28 * ab_y[k] * hf_6[k]
                  - f_30 * ab_y[k] * hf_8[k]
                  + f_31 * ab_z[k] * hf_9[k]
                  - f_22 * ab_x[k] * hf_30[k]
                  - f_25 * ab_x[k] * hf_33[k]
                  + f_26 * ab_x[k] * hf_35[k]
                  - f_22 * ab_y[k] * hf_36[k]
                  + f_26 * ab_y[k] * hf_38[k]
                  - f_27 * ab_z[k] * hf_39[k]
                  + f_21 * ab_x[k] * hf_100[k]
                  + f_22 * ab_x[k] * hf_103[k]
                  - f_23 * ab_x[k] * hf_105[k]
                  + f_21 * ab_y[k] * hf_106[k]
                  - f_23 * ab_y[k] * hf_108[k]
                  + f_24 * ab_z[k] * hf_109[k]
                  + f_28 * if__0[k]
                  + f_29 * if__3[k]
                  - f_30 * if__5[k]
                  + f_28 * if__16[k]
                  - f_30 * if__18[k]
                  + f_31 * if__29[k]
                  - f_22 * if__30[k]
                  - f_25 * if__33[k]
                  + f_26 * if__35[k]
                  - f_22 * if__66[k]
                  + f_26 * if__68[k]
                  - f_27 * if__79[k]
                  + f_21 * if__100[k]
                  + f_22 * if__103[k]
                  - f_23 * if__105[k]
                  + f_21 * if__156[k]
                  - f_23 * if__158[k]
                  + f_24 * if__169[k];
    }

#pragma omp simd aligned(ab_x, hf_2, hf_7, hf_9, hf_32, hf_37, hf_39, hf_102, hf_107, hf_109, \
                         if__2, if__7, if__9, if__32, if__37, if__39, if__102, if__107, \
                         if__109 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = -f_19 * ab_x[k] * hf_2[k]
                  - f_19 * ab_x[k] * hf_7[k]
                  + f_20 * ab_x[k] * hf_9[k]
                  + f_17 * ab_x[k] * hf_32[k]
                  + f_17 * ab_x[k] * hf_37[k]
                  - f_18 * ab_x[k] * hf_39[k]
                  - f_15 * ab_x[k] * hf_102[k]
                  - f_15 * ab_x[k] * hf_107[k]
                  + f_16 * ab_x[k] * hf_109[k]
                  - f_19 * if__2[k]
                  - f_19 * if__7[k]
                  + f_20 * if__9[k]
                  + f_17 * if__32[k]
                  + f_17 * if__37[k]
                  - f_18 * if__39[k]
                  - f_15 * if__102[k]
                  - f_15 * if__107[k]
                  + f_16 * if__109[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_0, hf_5, hf_6, hf_8, hf_30, hf_35, hf_36, hf_38, \
                         hf_100, hf_105, hf_106, hf_108, if__0, if__5, if__16, if__18, if__30, \
                         if__35, if__66, if__68, if__100, if__105, if__156, \
                         if__158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = -f_34 * ab_x[k] * hf_0[k]
                  + f_35 * ab_x[k] * hf_5[k]
                  + f_34 * ab_y[k] * hf_6[k]
                  - f_35 * ab_y[k] * hf_8[k]
                  + f_9 * ab_x[k] * hf_30[k]
                  - f_10 * ab_x[k] * hf_35[k]
                  - f_9 * ab_y[k] * hf_36[k]
                  + f_10 * ab_y[k] * hf_38[k]
                  - f_32 * ab_x[k] * hf_100[k]
                  + f_33 * ab_x[k] * hf_105[k]
                  + f_32 * ab_y[k] * hf_106[k]
                  - f_33 * ab_y[k] * hf_108[k]
                  - f_34 * if__0[k]
                  + f_35 * if__5[k]
                  + f_34 * if__16[k]
                  - f_35 * if__18[k]
                  + f_9 * if__30[k]
                  - f_10 * if__35[k]
                  - f_9 * if__66[k]
                  + f_10 * if__68[k]
                  - f_32 * if__100[k]
                  + f_33 * if__105[k]
                  + f_32 * if__156[k]
                  - f_33 * if__158[k];
    }

#pragma omp simd aligned(ab_x, hf_2, hf_7, hf_32, hf_37, hf_102, hf_107, if__2, if__7, if__32, \
                         if__37, if__102, if__107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = f_8 * ab_x[k] * hf_2[k]
                  - f_7 * ab_x[k] * hf_7[k]
                  - f_6 * ab_x[k] * hf_32[k]
                  + f_5 * ab_x[k] * hf_37[k]
                  + f_4 * ab_x[k] * hf_102[k]
                  - f_3 * ab_x[k] * hf_107[k]
                  + f_8 * if__2[k]
                  - f_7 * if__7[k]
                  - f_6 * if__32[k]
                  + f_5 * if__37[k]
                  + f_4 * if__102[k]
                  - f_3 * if__107[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hf_0, hf_3, hf_6, hf_30, hf_33, hf_36, hf_100, hf_103, \
                         hf_106, if__0, if__3, if__16, if__30, if__33, if__66, if__100, \
                         if__103, if__156 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = f_40 * ab_x[k] * hf_0[k]
                  - f_41 * ab_x[k] * hf_3[k]
                  + f_40 * ab_y[k] * hf_6[k]
                  - f_38 * ab_x[k] * hf_30[k]
                  + f_39 * ab_x[k] * hf_33[k]
                  - f_38 * ab_y[k] * hf_36[k]
                  + f_36 * ab_x[k] * hf_100[k]
                  - f_37 * ab_x[k] * hf_103[k]
                  + f_36 * ab_y[k] * hf_106[k]
                  + f_40 * if__0[k]
                  - f_41 * if__3[k]
                  + f_40 * if__16[k]
                  - f_38 * if__30[k]
                  + f_39 * if__33[k]
                  - f_38 * if__66[k]
                  + f_36 * if__100[k]
                  - f_37 * if__103[k]
                  + f_36 * if__156[k];
    }
}

}  // namespace simdtrf
