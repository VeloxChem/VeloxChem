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


#include "SimdTransformHF.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_hf(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t hf,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 1.40625 * std::sqrt(35.0);
    const auto f_1 = 0.46875 * std::sqrt(35.0);
    const auto f_2 = 2.8125 * std::sqrt(35.0);
    const auto f_3 = 0.9375 * std::sqrt(35.0);
    const auto f_4 = 0.28125 * std::sqrt(35.0);
    const auto f_5 = 0.09375 * std::sqrt(35.0);
    const auto f_6 = 0.9375 * std::sqrt(210.0);
    const auto f_7 = 1.875 * std::sqrt(210.0);
    const auto f_8 = 0.1875 * std::sqrt(210.0);
    const auto f_9 = 0.46875 * std::sqrt(21.0);
    const auto f_10 = 1.875 * std::sqrt(21.0);
    const auto f_11 = 0.9375 * std::sqrt(21.0);
    const auto f_12 = 3.75 * std::sqrt(21.0);
    const auto f_13 = 0.09375 * std::sqrt(21.0);
    const auto f_14 = 0.375 * std::sqrt(21.0);
    const auto f_15 = 1.40625 * std::sqrt(14.0);
    const auto f_16 = 0.9375 * std::sqrt(14.0);
    const auto f_17 = 2.8125 * std::sqrt(14.0);
    const auto f_18 = 1.875 * std::sqrt(14.0);
    const auto f_19 = 0.28125 * std::sqrt(14.0);
    const auto f_20 = 0.1875 * std::sqrt(14.0);
    const auto f_21 = 0.46875 * std::sqrt(210.0);
    const auto f_22 = 0.09375 * std::sqrt(210.0);
    const auto f_23 = 5.625 * std::sqrt(14.0);
    const auto f_24 = 7.5 * std::sqrt(21.0);
    const auto f_25 = 0.375 * std::sqrt(210.0);
    const auto f_26 = 1.5 * std::sqrt(210.0);
    const auto f_27 = 2.25 * std::sqrt(35.0);
    const auto f_28 = 1.5 * std::sqrt(35.0);
    const auto f_29 = 1.40625 * std::sqrt(7.0);
    const auto f_30 = 0.46875 * std::sqrt(7.0);
    const auto f_31 = 0.9375 * std::sqrt(7.0);
    const auto f_32 = 0.3125 * std::sqrt(7.0);
    const auto f_33 = 11.25 * std::sqrt(7.0);
    const auto f_34 = 3.75 * std::sqrt(7.0);
    const auto f_35 = 0.15625 * std::sqrt(7.0);
    const auto f_36 = 1.25 * std::sqrt(7.0);
    const auto f_37 = 0.9375 * std::sqrt(42.0);
    const auto f_38 = 0.625 * std::sqrt(42.0);
    const auto f_39 = 7.5 * std::sqrt(42.0);
    const auto f_40 = 0.3125 * std::sqrt(42.0);
    const auto f_41 = 2.5 * std::sqrt(42.0);
    const auto f_42 = 0.09375 * std::sqrt(105.0);
    const auto f_43 = 0.375 * std::sqrt(105.0);
    const auto f_44 = 0.0625 * std::sqrt(105.0);
    const auto f_45 = 0.25 * std::sqrt(105.0);
    const auto f_46 = 0.75 * std::sqrt(105.0);
    const auto f_47 = 3.0 * std::sqrt(105.0);
    const auto f_48 = 0.03125 * std::sqrt(105.0);
    const auto f_49 = 0.125 * std::sqrt(105.0);
    const auto f_50 = std::sqrt(105.0);
    const auto f_51 = 0.28125 * std::sqrt(70.0);
    const auto f_52 = 0.1875 * std::sqrt(70.0);
    const auto f_53 = 0.125 * std::sqrt(70.0);
    const auto f_54 = 2.25 * std::sqrt(70.0);
    const auto f_55 = 1.5 * std::sqrt(70.0);
    const auto f_56 = 0.09375 * std::sqrt(70.0);
    const auto f_57 = 0.0625 * std::sqrt(70.0);
    const auto f_58 = 0.75 * std::sqrt(70.0);
    const auto f_59 = 0.5 * std::sqrt(70.0);
    const auto f_60 = 0.46875 * std::sqrt(42.0);
    const auto f_61 = 3.75 * std::sqrt(42.0);
    const auto f_62 = 0.15625 * std::sqrt(42.0);
    const auto f_63 = 1.25 * std::sqrt(42.0);
    const auto f_64 = 1.875 * std::sqrt(42.0);
    const auto f_65 = 7.5 * std::sqrt(7.0);
    const auto f_66 = 15.0 * std::sqrt(7.0);
    const auto f_67 = 0.375 * std::sqrt(70.0);
    const auto f_68 = 3.0 * std::sqrt(70.0);
    const auto f_69 = 0.5 * std::sqrt(105.0);
    const auto f_70 = 1.5 * std::sqrt(105.0);
    const auto f_71 = 0.46875 * std::sqrt(6.0);
    const auto f_72 = 0.15625 * std::sqrt(6.0);
    const auto f_73 = 0.9375 * std::sqrt(6.0);
    const auto f_74 = 0.3125 * std::sqrt(6.0);
    const auto f_75 = 5.625 * std::sqrt(6.0);
    const auto f_76 = 1.875 * std::sqrt(6.0);
    const auto f_77 = 3.75 * std::sqrt(6.0);
    const auto f_78 = 1.25 * std::sqrt(6.0);
    const auto f_79 = 0.09375 * std::sqrt(10.0);
    const auto f_80 = 0.375 * std::sqrt(10.0);
    const auto f_81 = 0.1875 * std::sqrt(10.0);
    const auto f_82 = 0.75 * std::sqrt(10.0);
    const auto f_83 = 1.125 * std::sqrt(10.0);
    const auto f_84 = 4.5 * std::sqrt(10.0);
    const auto f_85 = 3.0 * std::sqrt(10.0);
    const auto f_86 = 0.1875 * std::sqrt(15.0);
    const auto f_87 = 0.125 * std::sqrt(15.0);
    const auto f_88 = 0.375 * std::sqrt(15.0);
    const auto f_89 = 0.25 * std::sqrt(15.0);
    const auto f_90 = 2.25 * std::sqrt(15.0);
    const auto f_91 = 1.5 * std::sqrt(15.0);
    const auto f_92 = std::sqrt(15.0);
    const auto f_93 = 1.40625 * std::sqrt(10.0);
    const auto f_94 = 0.46875 * std::sqrt(10.0);
    const auto f_95 = 2.8125 * std::sqrt(10.0);
    const auto f_96 = 0.9375 * std::sqrt(10.0);
    const auto f_97 = 3.75 * std::sqrt(10.0);
    const auto f_98 = 1.25 * std::sqrt(10.0);
    const auto f_99 = 0.25 * std::sqrt(10.0);
    const auto f_100 = 1.875 * std::sqrt(15.0);
    const auto f_101 = 3.75 * std::sqrt(15.0);
    const auto f_102 = 5.0 * std::sqrt(15.0);
    const auto f_103 = 5.0 * std::sqrt(6.0);
    const auto f_104 = 0.25 * std::sqrt(6.0);
    const auto f_105 = std::sqrt(6.0);
    const auto f_106 = 0.9375 * std::sqrt(15.0);
    const auto f_107 = 2.5 * std::sqrt(15.0);
    const auto f_108 = 0.5 * std::sqrt(15.0);
    const auto f_109 = 1.875 * std::sqrt(7.0);
    const auto f_110 = 0.46875 * std::sqrt(14.0);
    const auto f_111 = 8.4375 * std::sqrt(14.0);
    const auto f_112 = 11.25 * std::sqrt(21.0);
    const auto f_113 = 0.5625 * std::sqrt(210.0);
    const auto f_114 = 2.25 * std::sqrt(210.0);
    const auto f_115 = 0.5625 * std::sqrt(35.0);
    const auto f_116 = 0.375 * std::sqrt(35.0);
    const auto f_117 = 3.375 * std::sqrt(35.0);
    const auto f_118 = 5.625 * std::sqrt(21.0);

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

#pragma omp simd aligned(hf_11, hf_14, hf_16, hf_18, hf_61, hf_64, hf_66, hf_68, hf_151, \
                         hf_154, hf_156, hf_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * hf_11[k]
                 - f_1 * hf_16[k]
                 - f_2 * hf_61[k]
                 + f_3 * hf_66[k]
                 + f_4 * hf_151[k]
                 - f_5 * hf_156[k];

        g_1[k] = f_6 * hf_14[k]
                 - f_7 * hf_64[k]
                 + f_8 * hf_154[k];

        g_2[k] = -f_9 * hf_11[k]
                 - f_9 * hf_16[k]
                 + f_10 * hf_18[k]
                 + f_11 * hf_61[k]
                 + f_11 * hf_66[k]
                 - f_12 * hf_68[k]
                 - f_13 * hf_151[k]
                 - f_13 * hf_156[k]
                 + f_14 * hf_158[k];
    }

#pragma omp simd aligned(hf_12, hf_17, hf_19, hf_62, hf_67, hf_69, hf_152, hf_157, \
                         hf_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_15 * hf_12[k]
                 - f_15 * hf_17[k]
                 + f_16 * hf_19[k]
                 + f_17 * hf_62[k]
                 + f_17 * hf_67[k]
                 - f_18 * hf_69[k]
                 - f_19 * hf_152[k]
                 - f_19 * hf_157[k]
                 + f_20 * hf_159[k];
    }

#pragma omp simd aligned(hf_10, hf_13, hf_15, hf_60, hf_63, hf_65, hf_150, hf_153, \
                         hf_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_9 * hf_10[k]
                 - f_9 * hf_13[k]
                 + f_10 * hf_15[k]
                 + f_11 * hf_60[k]
                 + f_11 * hf_63[k]
                 - f_12 * hf_65[k]
                 - f_13 * hf_150[k]
                 - f_13 * hf_153[k]
                 + f_14 * hf_155[k];
    }

#pragma omp simd aligned(hf_10, hf_12, hf_13, hf_17, hf_60, hf_62, hf_63, hf_67, hf_150, \
                         hf_152, hf_153, hf_157 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_21 * hf_12[k]
                 - f_21 * hf_17[k]
                 - f_6 * hf_62[k]
                 + f_6 * hf_67[k]
                 + f_22 * hf_152[k]
                 - f_22 * hf_157[k];

        g_6[k] = f_1 * hf_10[k]
                 - f_0 * hf_13[k]
                 - f_3 * hf_60[k]
                 + f_2 * hf_63[k]
                 + f_5 * hf_150[k]
                 - f_4 * hf_153[k];
    }

#pragma omp simd aligned(hf_41, hf_44, hf_46, hf_48, hf_111, hf_114, hf_116, \
                         hf_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_23 * hf_41[k]
                 - f_18 * hf_46[k]
                 - f_23 * hf_111[k]
                 + f_18 * hf_116[k];

        g_8[k] = f_24 * hf_44[k]
                 - f_24 * hf_114[k];

        g_9[k] = -f_25 * hf_41[k]
                 - f_25 * hf_46[k]
                 + f_26 * hf_48[k]
                 + f_25 * hf_111[k]
                 + f_25 * hf_116[k]
                 - f_26 * hf_118[k];
    }

#pragma omp simd aligned(hf_40, hf_42, hf_43, hf_45, hf_47, hf_49, hf_110, hf_112, hf_113, \
                         hf_115, hf_117, hf_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_27 * hf_42[k]
                  - f_27 * hf_47[k]
                  + f_28 * hf_49[k]
                  + f_27 * hf_112[k]
                  + f_27 * hf_117[k]
                  - f_28 * hf_119[k];

        g_11[k] = -f_25 * hf_40[k]
                  - f_25 * hf_43[k]
                  + f_26 * hf_45[k]
                  + f_25 * hf_110[k]
                  + f_25 * hf_113[k]
                  - f_26 * hf_115[k];

        g_12[k] = f_12 * hf_42[k]
                  - f_12 * hf_47[k]
                  - f_12 * hf_112[k]
                  + f_12 * hf_117[k];

        g_13[k] = f_18 * hf_40[k]
                  - f_23 * hf_43[k]
                  - f_18 * hf_110[k]
                  + f_23 * hf_113[k];
    }

#pragma omp simd aligned(hf_11, hf_16, hf_61, hf_66, hf_81, hf_86, hf_151, hf_156, hf_171, \
                         hf_176 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_29 * hf_11[k]
                  + f_30 * hf_16[k]
                  - f_31 * hf_61[k]
                  + f_32 * hf_66[k]
                  + f_33 * hf_81[k]
                  - f_34 * hf_86[k]
                  + f_30 * hf_151[k]
                  - f_35 * hf_156[k]
                  - f_34 * hf_171[k]
                  + f_36 * hf_176[k];
    }

#pragma omp simd aligned(hf_14, hf_64, hf_84, hf_154, hf_174 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_37 * hf_14[k]
                  - f_38 * hf_64[k]
                  + f_39 * hf_84[k]
                  + f_40 * hf_154[k]
                  - f_41 * hf_174[k];
    }

#pragma omp simd aligned(hf_11, hf_16, hf_18, hf_61, hf_66, hf_68, hf_81, hf_86, hf_88, \
                         hf_151, hf_156, hf_158, hf_171, hf_176, \
                         hf_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_42 * hf_11[k]
                  + f_42 * hf_16[k]
                  - f_43 * hf_18[k]
                  + f_44 * hf_61[k]
                  + f_44 * hf_66[k]
                  - f_45 * hf_68[k]
                  - f_46 * hf_81[k]
                  - f_46 * hf_86[k]
                  + f_47 * hf_88[k]
                  - f_48 * hf_151[k]
                  - f_48 * hf_156[k]
                  + f_49 * hf_158[k]
                  + f_45 * hf_171[k]
                  + f_45 * hf_176[k]
                  - f_50 * hf_178[k];
    }

#pragma omp simd aligned(hf_12, hf_17, hf_19, hf_62, hf_67, hf_69, hf_82, hf_87, hf_89, \
                         hf_152, hf_157, hf_159, hf_172, hf_177, \
                         hf_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_51 * hf_12[k]
                  + f_51 * hf_17[k]
                  - f_52 * hf_19[k]
                  + f_52 * hf_62[k]
                  + f_52 * hf_67[k]
                  - f_53 * hf_69[k]
                  - f_54 * hf_82[k]
                  - f_54 * hf_87[k]
                  + f_55 * hf_89[k]
                  - f_56 * hf_152[k]
                  - f_56 * hf_157[k]
                  + f_57 * hf_159[k]
                  + f_58 * hf_172[k]
                  + f_58 * hf_177[k]
                  - f_59 * hf_179[k];
    }

#pragma omp simd aligned(hf_10, hf_13, hf_15, hf_60, hf_63, hf_65, hf_80, hf_83, hf_85, \
                         hf_150, hf_153, hf_155, hf_170, hf_173, \
                         hf_175 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_42 * hf_10[k]
                  + f_42 * hf_13[k]
                  - f_43 * hf_15[k]
                  + f_44 * hf_60[k]
                  + f_44 * hf_63[k]
                  - f_45 * hf_65[k]
                  - f_46 * hf_80[k]
                  - f_46 * hf_83[k]
                  + f_47 * hf_85[k]
                  - f_48 * hf_150[k]
                  - f_48 * hf_153[k]
                  + f_49 * hf_155[k]
                  + f_45 * hf_170[k]
                  + f_45 * hf_173[k]
                  - f_50 * hf_175[k];
    }

#pragma omp simd aligned(hf_12, hf_17, hf_62, hf_67, hf_82, hf_87, hf_152, hf_157, hf_172, \
                         hf_177 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_60 * hf_12[k]
                  + f_60 * hf_17[k]
                  - f_40 * hf_62[k]
                  + f_40 * hf_67[k]
                  + f_61 * hf_82[k]
                  - f_61 * hf_87[k]
                  + f_62 * hf_152[k]
                  - f_62 * hf_157[k]
                  - f_63 * hf_172[k]
                  + f_63 * hf_177[k];
    }

#pragma omp simd aligned(hf_10, hf_13, hf_60, hf_63, hf_80, hf_83, hf_150, hf_153, hf_170, \
                         hf_173 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_30 * hf_10[k]
                  + f_29 * hf_13[k]
                  - f_32 * hf_60[k]
                  + f_31 * hf_63[k]
                  + f_34 * hf_80[k]
                  - f_33 * hf_83[k]
                  + f_35 * hf_150[k]
                  - f_30 * hf_153[k]
                  - f_36 * hf_170[k]
                  + f_34 * hf_173[k];
    }

#pragma omp simd aligned(hf_41, hf_44, hf_46, hf_48, hf_111, hf_114, hf_116, hf_118, hf_131, \
                         hf_134, hf_136, hf_138 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_64 * hf_41[k]
                  + f_38 * hf_46[k]
                  - f_64 * hf_111[k]
                  + f_38 * hf_116[k]
                  + f_61 * hf_131[k]
                  - f_63 * hf_136[k];

        g_22[k] = -f_65 * hf_44[k]
                  - f_65 * hf_114[k]
                  + f_66 * hf_134[k];

        g_23[k] = f_67 * hf_41[k]
                  + f_67 * hf_46[k]
                  - f_55 * hf_48[k]
                  + f_67 * hf_111[k]
                  + f_67 * hf_116[k]
                  - f_55 * hf_118[k]
                  - f_58 * hf_131[k]
                  - f_58 * hf_136[k]
                  + f_68 * hf_138[k];
    }

#pragma omp simd aligned(hf_42, hf_47, hf_49, hf_112, hf_117, hf_119, hf_132, hf_137, \
                         hf_139 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_46 * hf_42[k]
                  + f_46 * hf_47[k]
                  - f_69 * hf_49[k]
                  + f_46 * hf_112[k]
                  + f_46 * hf_117[k]
                  - f_69 * hf_119[k]
                  - f_70 * hf_132[k]
                  - f_70 * hf_137[k]
                  + f_50 * hf_139[k];
    }

#pragma omp simd aligned(hf_40, hf_43, hf_45, hf_110, hf_113, hf_115, hf_130, hf_133, \
                         hf_135 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_67 * hf_40[k]
                  + f_67 * hf_43[k]
                  - f_55 * hf_45[k]
                  + f_67 * hf_110[k]
                  + f_67 * hf_113[k]
                  - f_55 * hf_115[k]
                  - f_58 * hf_130[k]
                  - f_58 * hf_133[k]
                  + f_68 * hf_135[k];
    }

#pragma omp simd aligned(hf_40, hf_42, hf_43, hf_47, hf_110, hf_112, hf_113, hf_117, hf_130, \
                         hf_132, hf_133, hf_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_34 * hf_42[k]
                  + f_34 * hf_47[k]
                  - f_34 * hf_112[k]
                  + f_34 * hf_117[k]
                  + f_65 * hf_132[k]
                  - f_65 * hf_137[k];

        g_27[k] = -f_38 * hf_40[k]
                  + f_64 * hf_43[k]
                  - f_38 * hf_110[k]
                  + f_64 * hf_113[k]
                  + f_63 * hf_130[k]
                  - f_61 * hf_133[k];
    }

#pragma omp simd aligned(hf_11, hf_16, hf_61, hf_66, hf_81, hf_86, hf_151, hf_156, hf_171, \
                         hf_176, hf_191, hf_196 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_71 * hf_11[k]
                  - f_72 * hf_16[k]
                  + f_73 * hf_61[k]
                  - f_74 * hf_66[k]
                  - f_75 * hf_81[k]
                  + f_76 * hf_86[k]
                  + f_71 * hf_151[k]
                  - f_72 * hf_156[k]
                  - f_75 * hf_171[k]
                  + f_76 * hf_176[k]
                  + f_77 * hf_191[k]
                  - f_78 * hf_196[k];
    }

#pragma omp simd aligned(hf_14, hf_64, hf_84, hf_154, hf_174, hf_194 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = 1.875 * hf_14[k]
                  + 3.75 * hf_64[k]
                  - 22.5 * hf_84[k]
                  + 1.875 * hf_154[k]
                  - 22.5 * hf_174[k]
                  + 15.0 * hf_194[k];
    }

#pragma omp simd aligned(hf_11, hf_16, hf_18, hf_61, hf_66, hf_68, hf_81, hf_86, hf_88, \
                         hf_151, hf_156, hf_158, hf_171, hf_176, hf_178, hf_191, hf_196, \
                         hf_198 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_79 * hf_11[k]
                  - f_79 * hf_16[k]
                  + f_80 * hf_18[k]
                  - f_81 * hf_61[k]
                  - f_81 * hf_66[k]
                  + f_82 * hf_68[k]
                  + f_83 * hf_81[k]
                  + f_83 * hf_86[k]
                  - f_84 * hf_88[k]
                  - f_79 * hf_151[k]
                  - f_79 * hf_156[k]
                  + f_80 * hf_158[k]
                  + f_83 * hf_171[k]
                  + f_83 * hf_176[k]
                  - f_84 * hf_178[k]
                  - f_82 * hf_191[k]
                  - f_82 * hf_196[k]
                  + f_85 * hf_198[k];
    }

#pragma omp simd aligned(hf_12, hf_17, hf_19, hf_62, hf_67, hf_69, hf_82, hf_87, hf_89, \
                         hf_152, hf_157, hf_159, hf_172, hf_177, hf_179, hf_192, hf_197, \
                         hf_199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_86 * hf_12[k]
                  - f_86 * hf_17[k]
                  + f_87 * hf_19[k]
                  - f_88 * hf_62[k]
                  - f_88 * hf_67[k]
                  + f_89 * hf_69[k]
                  + f_90 * hf_82[k]
                  + f_90 * hf_87[k]
                  - f_91 * hf_89[k]
                  - f_86 * hf_152[k]
                  - f_86 * hf_157[k]
                  + f_87 * hf_159[k]
                  + f_90 * hf_172[k]
                  + f_90 * hf_177[k]
                  - f_91 * hf_179[k]
                  - f_91 * hf_192[k]
                  - f_91 * hf_197[k]
                  + f_92 * hf_199[k];
    }

#pragma omp simd aligned(hf_10, hf_13, hf_15, hf_60, hf_63, hf_65, hf_80, hf_83, hf_85, \
                         hf_150, hf_153, hf_155, hf_170, hf_173, hf_175, hf_190, hf_193, \
                         hf_195 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_79 * hf_10[k]
                  - f_79 * hf_13[k]
                  + f_80 * hf_15[k]
                  - f_81 * hf_60[k]
                  - f_81 * hf_63[k]
                  + f_82 * hf_65[k]
                  + f_83 * hf_80[k]
                  + f_83 * hf_83[k]
                  - f_84 * hf_85[k]
                  - f_79 * hf_150[k]
                  - f_79 * hf_153[k]
                  + f_80 * hf_155[k]
                  + f_83 * hf_170[k]
                  + f_83 * hf_173[k]
                  - f_84 * hf_175[k]
                  - f_82 * hf_190[k]
                  - f_82 * hf_193[k]
                  + f_85 * hf_195[k];
    }

#pragma omp simd aligned(hf_12, hf_17, hf_62, hf_67, hf_82, hf_87, hf_152, hf_157, hf_172, \
                         hf_177, hf_192, hf_197 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = 0.9375 * hf_12[k]
                  - 0.9375 * hf_17[k]
                  + 1.875 * hf_62[k]
                  - 1.875 * hf_67[k]
                  - 11.25 * hf_82[k]
                  + 11.25 * hf_87[k]
                  + 0.9375 * hf_152[k]
                  - 0.9375 * hf_157[k]
                  - 11.25 * hf_172[k]
                  + 11.25 * hf_177[k]
                  + 7.5 * hf_192[k]
                  - 7.5 * hf_197[k];
    }

#pragma omp simd aligned(hf_10, hf_13, hf_60, hf_63, hf_80, hf_83, hf_150, hf_153, hf_170, \
                         hf_173, hf_190, hf_193 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_72 * hf_10[k]
                  - f_71 * hf_13[k]
                  + f_74 * hf_60[k]
                  - f_73 * hf_63[k]
                  - f_76 * hf_80[k]
                  + f_75 * hf_83[k]
                  + f_72 * hf_150[k]
                  - f_71 * hf_153[k]
                  - f_76 * hf_170[k]
                  + f_75 * hf_173[k]
                  + f_78 * hf_190[k]
                  - f_77 * hf_193[k];
    }

#pragma omp simd aligned(hf_21, hf_26, hf_71, hf_76, hf_91, hf_96, hf_161, hf_166, hf_181, \
                         hf_186, hf_201, hf_206 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_93 * hf_21[k]
                  - f_94 * hf_26[k]
                  + f_95 * hf_71[k]
                  - f_96 * hf_76[k]
                  - f_97 * hf_91[k]
                  + f_98 * hf_96[k]
                  + f_93 * hf_161[k]
                  - f_94 * hf_166[k]
                  - f_97 * hf_181[k]
                  + f_98 * hf_186[k]
                  + f_82 * hf_201[k]
                  - f_99 * hf_206[k];
    }

#pragma omp simd aligned(hf_24, hf_74, hf_94, hf_164, hf_184, hf_204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_100 * hf_24[k]
                  + f_101 * hf_74[k]
                  - f_102 * hf_94[k]
                  + f_100 * hf_164[k]
                  - f_102 * hf_184[k]
                  + f_92 * hf_204[k];
    }

#pragma omp simd aligned(hf_21, hf_26, hf_28, hf_71, hf_76, hf_78, hf_91, hf_96, hf_98, \
                         hf_161, hf_166, hf_168, hf_181, hf_186, hf_188, hf_201, hf_206, \
                         hf_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_71 * hf_21[k]
                  - f_71 * hf_26[k]
                  + f_76 * hf_28[k]
                  - f_73 * hf_71[k]
                  - f_73 * hf_76[k]
                  + f_77 * hf_78[k]
                  + f_78 * hf_91[k]
                  + f_78 * hf_96[k]
                  - f_103 * hf_98[k]
                  - f_71 * hf_161[k]
                  - f_71 * hf_166[k]
                  + f_76 * hf_168[k]
                  + f_78 * hf_181[k]
                  + f_78 * hf_186[k]
                  - f_103 * hf_188[k]
                  - f_104 * hf_201[k]
                  - f_104 * hf_206[k]
                  + f_105 * hf_208[k];
    }

#pragma omp simd aligned(hf_22, hf_27, hf_29, hf_72, hf_77, hf_79, hf_92, hf_97, hf_99, \
                         hf_162, hf_167, hf_169, hf_182, hf_187, hf_189, hf_202, hf_207, \
                         hf_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -2.8125 * hf_22[k]
                  - 2.8125 * hf_27[k]
                  + 1.875 * hf_29[k]
                  - 5.625 * hf_72[k]
                  - 5.625 * hf_77[k]
                  + 3.75 * hf_79[k]
                  + 7.5 * hf_92[k]
                  + 7.5 * hf_97[k]
                  - 5.0 * hf_99[k]
                  - 2.8125 * hf_162[k]
                  - 2.8125 * hf_167[k]
                  + 1.875 * hf_169[k]
                  + 7.5 * hf_182[k]
                  + 7.5 * hf_187[k]
                  - 5.0 * hf_189[k]
                  - 1.5 * hf_202[k]
                  - 1.5 * hf_207[k]
                  + hf_209[k];
    }

#pragma omp simd aligned(hf_20, hf_23, hf_25, hf_70, hf_73, hf_75, hf_90, hf_93, hf_95, \
                         hf_160, hf_163, hf_165, hf_180, hf_183, hf_185, hf_200, hf_203, \
                         hf_205 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_71 * hf_20[k]
                  - f_71 * hf_23[k]
                  + f_76 * hf_25[k]
                  - f_73 * hf_70[k]
                  - f_73 * hf_73[k]
                  + f_77 * hf_75[k]
                  + f_78 * hf_90[k]
                  + f_78 * hf_93[k]
                  - f_103 * hf_95[k]
                  - f_71 * hf_160[k]
                  - f_71 * hf_163[k]
                  + f_76 * hf_165[k]
                  + f_78 * hf_180[k]
                  + f_78 * hf_183[k]
                  - f_103 * hf_185[k]
                  - f_104 * hf_200[k]
                  - f_104 * hf_203[k]
                  + f_105 * hf_205[k];
    }

#pragma omp simd aligned(hf_22, hf_27, hf_72, hf_77, hf_92, hf_97, hf_162, hf_167, hf_182, \
                         hf_187, hf_202, hf_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_106 * hf_22[k]
                  - f_106 * hf_27[k]
                  + f_100 * hf_72[k]
                  - f_100 * hf_77[k]
                  - f_107 * hf_92[k]
                  + f_107 * hf_97[k]
                  + f_106 * hf_162[k]
                  - f_106 * hf_167[k]
                  - f_107 * hf_182[k]
                  + f_107 * hf_187[k]
                  + f_108 * hf_202[k]
                  - f_108 * hf_207[k];
    }

#pragma omp simd aligned(hf_20, hf_23, hf_70, hf_73, hf_90, hf_93, hf_160, hf_163, hf_180, \
                         hf_183, hf_200, hf_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_94 * hf_20[k]
                  - f_93 * hf_23[k]
                  + f_96 * hf_70[k]
                  - f_95 * hf_73[k]
                  - f_98 * hf_90[k]
                  + f_97 * hf_93[k]
                  + f_94 * hf_160[k]
                  - f_93 * hf_163[k]
                  - f_98 * hf_180[k]
                  + f_97 * hf_183[k]
                  + f_99 * hf_200[k]
                  - f_82 * hf_203[k];
    }

#pragma omp simd aligned(hf_1, hf_6, hf_31, hf_36, hf_51, hf_56, hf_101, hf_106, hf_121, \
                         hf_126, hf_141, hf_146 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_71 * hf_1[k]
                  - f_72 * hf_6[k]
                  + f_73 * hf_31[k]
                  - f_74 * hf_36[k]
                  - f_75 * hf_51[k]
                  + f_76 * hf_56[k]
                  + f_71 * hf_101[k]
                  - f_72 * hf_106[k]
                  - f_75 * hf_121[k]
                  + f_76 * hf_126[k]
                  + f_77 * hf_141[k]
                  - f_78 * hf_146[k];
    }

#pragma omp simd aligned(hf_4, hf_34, hf_54, hf_104, hf_124, hf_144 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = 1.875 * hf_4[k]
                  + 3.75 * hf_34[k]
                  - 22.5 * hf_54[k]
                  + 1.875 * hf_104[k]
                  - 22.5 * hf_124[k]
                  + 15.0 * hf_144[k];
    }

#pragma omp simd aligned(hf_1, hf_6, hf_8, hf_31, hf_36, hf_38, hf_51, hf_56, hf_58, hf_101, \
                         hf_106, hf_108, hf_121, hf_126, hf_128, hf_141, hf_146, \
                         hf_148 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_79 * hf_1[k]
                  - f_79 * hf_6[k]
                  + f_80 * hf_8[k]
                  - f_81 * hf_31[k]
                  - f_81 * hf_36[k]
                  + f_82 * hf_38[k]
                  + f_83 * hf_51[k]
                  + f_83 * hf_56[k]
                  - f_84 * hf_58[k]
                  - f_79 * hf_101[k]
                  - f_79 * hf_106[k]
                  + f_80 * hf_108[k]
                  + f_83 * hf_121[k]
                  + f_83 * hf_126[k]
                  - f_84 * hf_128[k]
                  - f_82 * hf_141[k]
                  - f_82 * hf_146[k]
                  + f_85 * hf_148[k];
    }

#pragma omp simd aligned(hf_2, hf_7, hf_9, hf_32, hf_37, hf_39, hf_52, hf_57, hf_59, hf_102, \
                         hf_107, hf_109, hf_122, hf_127, hf_129, hf_142, hf_147, \
                         hf_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_86 * hf_2[k]
                  - f_86 * hf_7[k]
                  + f_87 * hf_9[k]
                  - f_88 * hf_32[k]
                  - f_88 * hf_37[k]
                  + f_89 * hf_39[k]
                  + f_90 * hf_52[k]
                  + f_90 * hf_57[k]
                  - f_91 * hf_59[k]
                  - f_86 * hf_102[k]
                  - f_86 * hf_107[k]
                  + f_87 * hf_109[k]
                  + f_90 * hf_122[k]
                  + f_90 * hf_127[k]
                  - f_91 * hf_129[k]
                  - f_91 * hf_142[k]
                  - f_91 * hf_147[k]
                  + f_92 * hf_149[k];
    }

#pragma omp simd aligned(hf_0, hf_3, hf_5, hf_30, hf_33, hf_35, hf_50, hf_53, hf_55, hf_100, \
                         hf_103, hf_105, hf_120, hf_123, hf_125, hf_140, hf_143, \
                         hf_145 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_79 * hf_0[k]
                  - f_79 * hf_3[k]
                  + f_80 * hf_5[k]
                  - f_81 * hf_30[k]
                  - f_81 * hf_33[k]
                  + f_82 * hf_35[k]
                  + f_83 * hf_50[k]
                  + f_83 * hf_53[k]
                  - f_84 * hf_55[k]
                  - f_79 * hf_100[k]
                  - f_79 * hf_103[k]
                  + f_80 * hf_105[k]
                  + f_83 * hf_120[k]
                  + f_83 * hf_123[k]
                  - f_84 * hf_125[k]
                  - f_82 * hf_140[k]
                  - f_82 * hf_143[k]
                  + f_85 * hf_145[k];
    }

#pragma omp simd aligned(hf_2, hf_7, hf_32, hf_37, hf_52, hf_57, hf_102, hf_107, hf_122, \
                         hf_127, hf_142, hf_147 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = 0.9375 * hf_2[k]
                  - 0.9375 * hf_7[k]
                  + 1.875 * hf_32[k]
                  - 1.875 * hf_37[k]
                  - 11.25 * hf_52[k]
                  + 11.25 * hf_57[k]
                  + 0.9375 * hf_102[k]
                  - 0.9375 * hf_107[k]
                  - 11.25 * hf_122[k]
                  + 11.25 * hf_127[k]
                  + 7.5 * hf_142[k]
                  - 7.5 * hf_147[k];
    }

#pragma omp simd aligned(hf_0, hf_3, hf_30, hf_33, hf_50, hf_53, hf_100, hf_103, hf_120, \
                         hf_123, hf_140, hf_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_72 * hf_0[k]
                  - f_71 * hf_3[k]
                  + f_74 * hf_30[k]
                  - f_73 * hf_33[k]
                  - f_76 * hf_50[k]
                  + f_75 * hf_53[k]
                  + f_72 * hf_100[k]
                  - f_71 * hf_103[k]
                  - f_76 * hf_120[k]
                  + f_75 * hf_123[k]
                  + f_78 * hf_140[k]
                  - f_77 * hf_143[k];
    }

#pragma omp simd aligned(hf_21, hf_24, hf_26, hf_91, hf_94, hf_96, hf_161, hf_164, hf_166, \
                         hf_181, hf_184, hf_186 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_37 * hf_21[k]
                  + f_40 * hf_26[k]
                  + f_64 * hf_91[k]
                  - f_38 * hf_96[k]
                  + f_37 * hf_161[k]
                  - f_40 * hf_166[k]
                  - f_64 * hf_181[k]
                  + f_38 * hf_186[k];

        g_50[k] = -f_34 * hf_24[k]
                  + f_65 * hf_94[k]
                  + f_34 * hf_164[k]
                  - f_65 * hf_184[k];
    }

#pragma omp simd aligned(hf_21, hf_26, hf_28, hf_91, hf_96, hf_98, hf_161, hf_166, hf_168, \
                         hf_181, hf_186, hf_188 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_52 * hf_21[k]
                  + f_52 * hf_26[k]
                  - f_58 * hf_28[k]
                  - f_67 * hf_91[k]
                  - f_67 * hf_96[k]
                  + f_55 * hf_98[k]
                  - f_52 * hf_161[k]
                  - f_52 * hf_166[k]
                  + f_58 * hf_168[k]
                  + f_67 * hf_181[k]
                  + f_67 * hf_186[k]
                  - f_55 * hf_188[k];
    }

#pragma omp simd aligned(hf_22, hf_27, hf_29, hf_92, hf_97, hf_99, hf_162, hf_167, hf_169, \
                         hf_182, hf_187, hf_189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_43 * hf_22[k]
                  + f_43 * hf_27[k]
                  - f_45 * hf_29[k]
                  - f_46 * hf_92[k]
                  - f_46 * hf_97[k]
                  + f_69 * hf_99[k]
                  - f_43 * hf_162[k]
                  - f_43 * hf_167[k]
                  + f_45 * hf_169[k]
                  + f_46 * hf_182[k]
                  + f_46 * hf_187[k]
                  - f_69 * hf_189[k];
    }

#pragma omp simd aligned(hf_20, hf_23, hf_25, hf_90, hf_93, hf_95, hf_160, hf_163, hf_165, \
                         hf_180, hf_183, hf_185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_52 * hf_20[k]
                  + f_52 * hf_23[k]
                  - f_58 * hf_25[k]
                  - f_67 * hf_90[k]
                  - f_67 * hf_93[k]
                  + f_55 * hf_95[k]
                  - f_52 * hf_160[k]
                  - f_52 * hf_163[k]
                  + f_58 * hf_165[k]
                  + f_67 * hf_180[k]
                  + f_67 * hf_183[k]
                  - f_55 * hf_185[k];
    }

#pragma omp simd aligned(hf_22, hf_27, hf_92, hf_97, hf_162, hf_167, hf_182, \
                         hf_187 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_109 * hf_22[k]
                  + f_109 * hf_27[k]
                  + f_34 * hf_92[k]
                  - f_34 * hf_97[k]
                  + f_109 * hf_162[k]
                  - f_109 * hf_167[k]
                  - f_34 * hf_182[k]
                  + f_34 * hf_187[k];
    }

#pragma omp simd aligned(hf_20, hf_23, hf_90, hf_93, hf_160, hf_163, hf_180, \
                         hf_183 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_40 * hf_20[k]
                  + f_37 * hf_23[k]
                  + f_38 * hf_90[k]
                  - f_64 * hf_93[k]
                  + f_40 * hf_160[k]
                  - f_37 * hf_163[k]
                  - f_38 * hf_180[k]
                  + f_64 * hf_183[k];
    }

#pragma omp simd aligned(hf_1, hf_6, hf_31, hf_36, hf_51, hf_56, hf_101, hf_106, hf_121, \
                         hf_126 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = -f_30 * hf_1[k]
                  + f_35 * hf_6[k]
                  + f_31 * hf_31[k]
                  - f_32 * hf_36[k]
                  + f_34 * hf_51[k]
                  - f_36 * hf_56[k]
                  + f_29 * hf_101[k]
                  - f_30 * hf_106[k]
                  - f_33 * hf_121[k]
                  + f_34 * hf_126[k];
    }

#pragma omp simd aligned(hf_4, hf_34, hf_54, hf_104, hf_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = -f_40 * hf_4[k]
                  + f_38 * hf_34[k]
                  + f_41 * hf_54[k]
                  + f_37 * hf_104[k]
                  - f_39 * hf_124[k];
    }

#pragma omp simd aligned(hf_1, hf_6, hf_8, hf_31, hf_36, hf_38, hf_51, hf_56, hf_58, hf_101, \
                         hf_106, hf_108, hf_121, hf_126, hf_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = f_48 * hf_1[k]
                  + f_48 * hf_6[k]
                  - f_49 * hf_8[k]
                  - f_44 * hf_31[k]
                  - f_44 * hf_36[k]
                  + f_45 * hf_38[k]
                  - f_45 * hf_51[k]
                  - f_45 * hf_56[k]
                  + f_50 * hf_58[k]
                  - f_42 * hf_101[k]
                  - f_42 * hf_106[k]
                  + f_43 * hf_108[k]
                  + f_46 * hf_121[k]
                  + f_46 * hf_126[k]
                  - f_47 * hf_128[k];
    }

#pragma omp simd aligned(hf_2, hf_7, hf_9, hf_32, hf_37, hf_39, hf_52, hf_57, hf_59, hf_102, \
                         hf_107, hf_109, hf_122, hf_127, hf_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_56 * hf_2[k]
                  + f_56 * hf_7[k]
                  - f_57 * hf_9[k]
                  - f_52 * hf_32[k]
                  - f_52 * hf_37[k]
                  + f_53 * hf_39[k]
                  - f_58 * hf_52[k]
                  - f_58 * hf_57[k]
                  + f_59 * hf_59[k]
                  - f_51 * hf_102[k]
                  - f_51 * hf_107[k]
                  + f_52 * hf_109[k]
                  + f_54 * hf_122[k]
                  + f_54 * hf_127[k]
                  - f_55 * hf_129[k];
    }

#pragma omp simd aligned(hf_0, hf_3, hf_5, hf_30, hf_33, hf_35, hf_50, hf_53, hf_55, hf_100, \
                         hf_103, hf_105, hf_120, hf_123, hf_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_48 * hf_0[k]
                  + f_48 * hf_3[k]
                  - f_49 * hf_5[k]
                  - f_44 * hf_30[k]
                  - f_44 * hf_33[k]
                  + f_45 * hf_35[k]
                  - f_45 * hf_50[k]
                  - f_45 * hf_53[k]
                  + f_50 * hf_55[k]
                  - f_42 * hf_100[k]
                  - f_42 * hf_103[k]
                  + f_43 * hf_105[k]
                  + f_46 * hf_120[k]
                  + f_46 * hf_123[k]
                  - f_47 * hf_125[k];
    }

#pragma omp simd aligned(hf_2, hf_7, hf_32, hf_37, hf_52, hf_57, hf_102, hf_107, hf_122, \
                         hf_127 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_62 * hf_2[k]
                  + f_62 * hf_7[k]
                  + f_40 * hf_32[k]
                  - f_40 * hf_37[k]
                  + f_63 * hf_52[k]
                  - f_63 * hf_57[k]
                  + f_60 * hf_102[k]
                  - f_60 * hf_107[k]
                  - f_61 * hf_122[k]
                  + f_61 * hf_127[k];
    }

#pragma omp simd aligned(hf_0, hf_3, hf_30, hf_33, hf_50, hf_53, hf_100, hf_103, hf_120, \
                         hf_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_35 * hf_0[k]
                  + f_30 * hf_3[k]
                  + f_32 * hf_30[k]
                  - f_31 * hf_33[k]
                  + f_36 * hf_50[k]
                  - f_34 * hf_53[k]
                  + f_30 * hf_100[k]
                  - f_29 * hf_103[k]
                  - f_34 * hf_120[k]
                  + f_33 * hf_123[k];
    }

#pragma omp simd aligned(hf_21, hf_24, hf_26, hf_28, hf_71, hf_74, hf_76, hf_78, hf_161, \
                         hf_164, hf_166, hf_168 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_15 * hf_21[k]
                  - f_110 * hf_26[k]
                  - f_111 * hf_71[k]
                  + f_17 * hf_76[k]
                  + f_15 * hf_161[k]
                  - f_110 * hf_166[k];

        g_64[k] = f_10 * hf_24[k]
                  - f_112 * hf_74[k]
                  + f_10 * hf_164[k];

        g_65[k] = -f_22 * hf_21[k]
                  - f_22 * hf_26[k]
                  + f_25 * hf_28[k]
                  + f_113 * hf_71[k]
                  + f_113 * hf_76[k]
                  - f_114 * hf_78[k]
                  - f_22 * hf_161[k]
                  - f_22 * hf_166[k]
                  + f_25 * hf_168[k];
    }

#pragma omp simd aligned(hf_22, hf_27, hf_29, hf_72, hf_77, hf_79, hf_162, hf_167, \
                         hf_169 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = -f_115 * hf_22[k]
                  - f_115 * hf_27[k]
                  + f_116 * hf_29[k]
                  + f_117 * hf_72[k]
                  + f_117 * hf_77[k]
                  - f_27 * hf_79[k]
                  - f_115 * hf_162[k]
                  - f_115 * hf_167[k]
                  + f_116 * hf_169[k];
    }

#pragma omp simd aligned(hf_20, hf_23, hf_25, hf_70, hf_73, hf_75, hf_160, hf_163, \
                         hf_165 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_22 * hf_20[k]
                  - f_22 * hf_23[k]
                  + f_25 * hf_25[k]
                  + f_113 * hf_70[k]
                  + f_113 * hf_73[k]
                  - f_114 * hf_75[k]
                  - f_22 * hf_160[k]
                  - f_22 * hf_163[k]
                  + f_25 * hf_165[k];
    }

#pragma omp simd aligned(hf_20, hf_22, hf_23, hf_27, hf_70, hf_72, hf_73, hf_77, hf_160, \
                         hf_162, hf_163, hf_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_11 * hf_22[k]
                  - f_11 * hf_27[k]
                  - f_118 * hf_72[k]
                  + f_118 * hf_77[k]
                  + f_11 * hf_162[k]
                  - f_11 * hf_167[k];

        g_69[k] = f_110 * hf_20[k]
                  - f_15 * hf_23[k]
                  - f_17 * hf_70[k]
                  + f_111 * hf_73[k]
                  + f_110 * hf_160[k]
                  - f_15 * hf_163[k];
    }

#pragma omp simd aligned(hf_1, hf_4, hf_6, hf_8, hf_31, hf_34, hf_36, hf_38, hf_101, hf_104, \
                         hf_106, hf_108 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = f_4 * hf_1[k]
                  - f_5 * hf_6[k]
                  - f_2 * hf_31[k]
                  + f_3 * hf_36[k]
                  + f_0 * hf_101[k]
                  - f_1 * hf_106[k];

        g_71[k] = f_8 * hf_4[k]
                  - f_7 * hf_34[k]
                  + f_6 * hf_104[k];

        g_72[k] = -f_13 * hf_1[k]
                  - f_13 * hf_6[k]
                  + f_14 * hf_8[k]
                  + f_11 * hf_31[k]
                  + f_11 * hf_36[k]
                  - f_12 * hf_38[k]
                  - f_9 * hf_101[k]
                  - f_9 * hf_106[k]
                  + f_10 * hf_108[k];
    }

#pragma omp simd aligned(hf_2, hf_7, hf_9, hf_32, hf_37, hf_39, hf_102, hf_107, \
                         hf_109 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = -f_19 * hf_2[k]
                  - f_19 * hf_7[k]
                  + f_20 * hf_9[k]
                  + f_17 * hf_32[k]
                  + f_17 * hf_37[k]
                  - f_18 * hf_39[k]
                  - f_15 * hf_102[k]
                  - f_15 * hf_107[k]
                  + f_16 * hf_109[k];
    }

#pragma omp simd aligned(hf_0, hf_3, hf_5, hf_30, hf_33, hf_35, hf_100, hf_103, \
                         hf_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_13 * hf_0[k]
                  - f_13 * hf_3[k]
                  + f_14 * hf_5[k]
                  + f_11 * hf_30[k]
                  + f_11 * hf_33[k]
                  - f_12 * hf_35[k]
                  - f_9 * hf_100[k]
                  - f_9 * hf_103[k]
                  + f_10 * hf_105[k];
    }

#pragma omp simd aligned(hf_0, hf_2, hf_3, hf_7, hf_30, hf_32, hf_33, hf_37, hf_100, hf_102, \
                         hf_103, hf_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = f_22 * hf_2[k]
                  - f_22 * hf_7[k]
                  - f_6 * hf_32[k]
                  + f_6 * hf_37[k]
                  + f_21 * hf_102[k]
                  - f_21 * hf_107[k];

        g_76[k] = f_5 * hf_0[k]
                  - f_4 * hf_3[k]
                  - f_3 * hf_30[k]
                  + f_2 * hf_33[k]
                  + f_1 * hf_100[k]
                  - f_0 * hf_103[k];
    }
}

}  // namespace simdtrf
