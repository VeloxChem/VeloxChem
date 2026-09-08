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


#include "SimdTransformFH.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_fh(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t fh,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 1.40625 * std::sqrt(35.0);
    const auto f_1 = 2.8125 * std::sqrt(35.0);
    const auto f_2 = 0.28125 * std::sqrt(35.0);
    const auto f_3 = 0.46875 * std::sqrt(35.0);
    const auto f_4 = 0.9375 * std::sqrt(35.0);
    const auto f_5 = 0.09375 * std::sqrt(35.0);
    const auto f_6 = 5.625 * std::sqrt(14.0);
    const auto f_7 = 1.875 * std::sqrt(14.0);
    const auto f_8 = 1.40625 * std::sqrt(7.0);
    const auto f_9 = 0.9375 * std::sqrt(7.0);
    const auto f_10 = 11.25 * std::sqrt(7.0);
    const auto f_11 = 0.46875 * std::sqrt(7.0);
    const auto f_12 = 3.75 * std::sqrt(7.0);
    const auto f_13 = 0.3125 * std::sqrt(7.0);
    const auto f_14 = 0.15625 * std::sqrt(7.0);
    const auto f_15 = 1.25 * std::sqrt(7.0);
    const auto f_16 = 1.875 * std::sqrt(42.0);
    const auto f_17 = 3.75 * std::sqrt(42.0);
    const auto f_18 = 0.625 * std::sqrt(42.0);
    const auto f_19 = 1.25 * std::sqrt(42.0);
    const auto f_20 = 0.46875 * std::sqrt(6.0);
    const auto f_21 = 0.9375 * std::sqrt(6.0);
    const auto f_22 = 5.625 * std::sqrt(6.0);
    const auto f_23 = 3.75 * std::sqrt(6.0);
    const auto f_24 = 0.15625 * std::sqrt(6.0);
    const auto f_25 = 0.3125 * std::sqrt(6.0);
    const auto f_26 = 1.875 * std::sqrt(6.0);
    const auto f_27 = 1.25 * std::sqrt(6.0);
    const auto f_28 = 1.40625 * std::sqrt(10.0);
    const auto f_29 = 2.8125 * std::sqrt(10.0);
    const auto f_30 = 3.75 * std::sqrt(10.0);
    const auto f_31 = 0.75 * std::sqrt(10.0);
    const auto f_32 = 0.46875 * std::sqrt(10.0);
    const auto f_33 = 0.9375 * std::sqrt(10.0);
    const auto f_34 = 1.25 * std::sqrt(10.0);
    const auto f_35 = 0.25 * std::sqrt(10.0);
    const auto f_36 = 0.9375 * std::sqrt(42.0);
    const auto f_37 = 0.3125 * std::sqrt(42.0);
    const auto f_38 = 1.40625 * std::sqrt(14.0);
    const auto f_39 = 8.4375 * std::sqrt(14.0);
    const auto f_40 = 0.46875 * std::sqrt(14.0);
    const auto f_41 = 2.8125 * std::sqrt(14.0);
    const auto f_42 = 0.9375 * std::sqrt(210.0);
    const auto f_43 = 1.875 * std::sqrt(210.0);
    const auto f_44 = 0.1875 * std::sqrt(210.0);
    const auto f_45 = 7.5 * std::sqrt(21.0);
    const auto f_46 = 7.5 * std::sqrt(42.0);
    const auto f_47 = 2.5 * std::sqrt(42.0);
    const auto f_48 = 7.5 * std::sqrt(7.0);
    const auto f_49 = 15.0 * std::sqrt(7.0);
    const auto f_50 = 1.875 * std::sqrt(15.0);
    const auto f_51 = 3.75 * std::sqrt(15.0);
    const auto f_52 = 5.0 * std::sqrt(15.0);
    const auto f_53 = std::sqrt(15.0);
    const auto f_54 = 1.875 * std::sqrt(21.0);
    const auto f_55 = 11.25 * std::sqrt(21.0);
    const auto f_56 = 0.46875 * std::sqrt(21.0);
    const auto f_57 = 0.9375 * std::sqrt(21.0);
    const auto f_58 = 0.09375 * std::sqrt(21.0);
    const auto f_59 = 3.75 * std::sqrt(21.0);
    const auto f_60 = 0.375 * std::sqrt(21.0);
    const auto f_61 = 0.375 * std::sqrt(210.0);
    const auto f_62 = 1.5 * std::sqrt(210.0);
    const auto f_63 = 0.09375 * std::sqrt(105.0);
    const auto f_64 = 0.0625 * std::sqrt(105.0);
    const auto f_65 = 0.75 * std::sqrt(105.0);
    const auto f_66 = 0.03125 * std::sqrt(105.0);
    const auto f_67 = 0.25 * std::sqrt(105.0);
    const auto f_68 = 0.375 * std::sqrt(105.0);
    const auto f_69 = 3.0 * std::sqrt(105.0);
    const auto f_70 = 0.125 * std::sqrt(105.0);
    const auto f_71 = std::sqrt(105.0);
    const auto f_72 = 0.375 * std::sqrt(70.0);
    const auto f_73 = 0.75 * std::sqrt(70.0);
    const auto f_74 = 1.5 * std::sqrt(70.0);
    const auto f_75 = 3.0 * std::sqrt(70.0);
    const auto f_76 = 0.09375 * std::sqrt(10.0);
    const auto f_77 = 0.1875 * std::sqrt(10.0);
    const auto f_78 = 1.125 * std::sqrt(10.0);
    const auto f_79 = 0.375 * std::sqrt(10.0);
    const auto f_80 = 4.5 * std::sqrt(10.0);
    const auto f_81 = 3.0 * std::sqrt(10.0);
    const auto f_82 = 0.25 * std::sqrt(6.0);
    const auto f_83 = 5.0 * std::sqrt(6.0);
    const auto f_84 = std::sqrt(6.0);
    const auto f_85 = 0.1875 * std::sqrt(70.0);
    const auto f_86 = 0.09375 * std::sqrt(210.0);
    const auto f_87 = 0.5625 * std::sqrt(210.0);
    const auto f_88 = 2.25 * std::sqrt(210.0);
    const auto f_89 = 0.28125 * std::sqrt(14.0);
    const auto f_90 = 0.9375 * std::sqrt(14.0);
    const auto f_91 = 0.1875 * std::sqrt(14.0);
    const auto f_92 = 2.25 * std::sqrt(35.0);
    const auto f_93 = 1.5 * std::sqrt(35.0);
    const auto f_94 = 0.28125 * std::sqrt(70.0);
    const auto f_95 = 2.25 * std::sqrt(70.0);
    const auto f_96 = 0.09375 * std::sqrt(70.0);
    const auto f_97 = 0.125 * std::sqrt(70.0);
    const auto f_98 = 0.0625 * std::sqrt(70.0);
    const auto f_99 = 0.5 * std::sqrt(70.0);
    const auto f_100 = 1.5 * std::sqrt(105.0);
    const auto f_101 = 0.5 * std::sqrt(105.0);
    const auto f_102 = 0.1875 * std::sqrt(15.0);
    const auto f_103 = 0.375 * std::sqrt(15.0);
    const auto f_104 = 2.25 * std::sqrt(15.0);
    const auto f_105 = 1.5 * std::sqrt(15.0);
    const auto f_106 = 0.125 * std::sqrt(15.0);
    const auto f_107 = 0.25 * std::sqrt(15.0);
    const auto f_108 = 0.5625 * std::sqrt(35.0);
    const auto f_109 = 3.375 * std::sqrt(35.0);
    const auto f_110 = 0.375 * std::sqrt(35.0);
    const auto f_111 = 0.46875 * std::sqrt(210.0);
    const auto f_112 = 0.46875 * std::sqrt(42.0);
    const auto f_113 = 0.15625 * std::sqrt(42.0);
    const auto f_114 = 0.9375 * std::sqrt(15.0);
    const auto f_115 = 2.5 * std::sqrt(15.0);
    const auto f_116 = 0.5 * std::sqrt(15.0);
    const auto f_117 = 1.875 * std::sqrt(7.0);
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

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_37 = buffer.data(fh + 37);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_40 = buffer.data(fh + 40);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_45 = buffer.data(fh + 45);
    const auto *fh_46 = buffer.data(fh + 46);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_48 = buffer.data(fh + 48);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_50 = buffer.data(fh + 50);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_52 = buffer.data(fh + 52);
    const auto *fh_53 = buffer.data(fh + 53);
    const auto *fh_54 = buffer.data(fh + 54);
    const auto *fh_55 = buffer.data(fh + 55);
    const auto *fh_56 = buffer.data(fh + 56);
    const auto *fh_57 = buffer.data(fh + 57);
    const auto *fh_58 = buffer.data(fh + 58);
    const auto *fh_59 = buffer.data(fh + 59);
    const auto *fh_60 = buffer.data(fh + 60);
    const auto *fh_61 = buffer.data(fh + 61);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_65 = buffer.data(fh + 65);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_67 = buffer.data(fh + 67);
    const auto *fh_68 = buffer.data(fh + 68);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_70 = buffer.data(fh + 70);
    const auto *fh_71 = buffer.data(fh + 71);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_74 = buffer.data(fh + 74);
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_76 = buffer.data(fh + 76);
    const auto *fh_77 = buffer.data(fh + 77);
    const auto *fh_78 = buffer.data(fh + 78);
    const auto *fh_79 = buffer.data(fh + 79);
    const auto *fh_80 = buffer.data(fh + 80);
    const auto *fh_81 = buffer.data(fh + 81);
    const auto *fh_82 = buffer.data(fh + 82);
    const auto *fh_83 = buffer.data(fh + 83);
    const auto *fh_84 = buffer.data(fh + 84);
    const auto *fh_85 = buffer.data(fh + 85);
    const auto *fh_86 = buffer.data(fh + 86);
    const auto *fh_87 = buffer.data(fh + 87);
    const auto *fh_88 = buffer.data(fh + 88);
    const auto *fh_89 = buffer.data(fh + 89);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_91 = buffer.data(fh + 91);
    const auto *fh_92 = buffer.data(fh + 92);
    const auto *fh_93 = buffer.data(fh + 93);
    const auto *fh_94 = buffer.data(fh + 94);
    const auto *fh_95 = buffer.data(fh + 95);
    const auto *fh_96 = buffer.data(fh + 96);
    const auto *fh_97 = buffer.data(fh + 97);
    const auto *fh_98 = buffer.data(fh + 98);
    const auto *fh_99 = buffer.data(fh + 99);
    const auto *fh_100 = buffer.data(fh + 100);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_102 = buffer.data(fh + 102);
    const auto *fh_103 = buffer.data(fh + 103);
    const auto *fh_104 = buffer.data(fh + 104);
    const auto *fh_105 = buffer.data(fh + 105);
    const auto *fh_106 = buffer.data(fh + 106);
    const auto *fh_107 = buffer.data(fh + 107);
    const auto *fh_108 = buffer.data(fh + 108);
    const auto *fh_109 = buffer.data(fh + 109);
    const auto *fh_110 = buffer.data(fh + 110);
    const auto *fh_111 = buffer.data(fh + 111);
    const auto *fh_112 = buffer.data(fh + 112);
    const auto *fh_113 = buffer.data(fh + 113);
    const auto *fh_114 = buffer.data(fh + 114);
    const auto *fh_115 = buffer.data(fh + 115);
    const auto *fh_116 = buffer.data(fh + 116);
    const auto *fh_117 = buffer.data(fh + 117);
    const auto *fh_118 = buffer.data(fh + 118);
    const auto *fh_119 = buffer.data(fh + 119);
    const auto *fh_120 = buffer.data(fh + 120);
    const auto *fh_121 = buffer.data(fh + 121);
    const auto *fh_122 = buffer.data(fh + 122);
    const auto *fh_123 = buffer.data(fh + 123);
    const auto *fh_124 = buffer.data(fh + 124);
    const auto *fh_125 = buffer.data(fh + 125);
    const auto *fh_126 = buffer.data(fh + 126);
    const auto *fh_127 = buffer.data(fh + 127);
    const auto *fh_128 = buffer.data(fh + 128);
    const auto *fh_129 = buffer.data(fh + 129);
    const auto *fh_130 = buffer.data(fh + 130);
    const auto *fh_131 = buffer.data(fh + 131);
    const auto *fh_132 = buffer.data(fh + 132);
    const auto *fh_133 = buffer.data(fh + 133);
    const auto *fh_134 = buffer.data(fh + 134);
    const auto *fh_135 = buffer.data(fh + 135);
    const auto *fh_136 = buffer.data(fh + 136);
    const auto *fh_137 = buffer.data(fh + 137);
    const auto *fh_138 = buffer.data(fh + 138);
    const auto *fh_139 = buffer.data(fh + 139);
    const auto *fh_140 = buffer.data(fh + 140);
    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_142 = buffer.data(fh + 142);
    const auto *fh_143 = buffer.data(fh + 143);
    const auto *fh_144 = buffer.data(fh + 144);
    const auto *fh_145 = buffer.data(fh + 145);
    const auto *fh_146 = buffer.data(fh + 146);
    const auto *fh_147 = buffer.data(fh + 147);
    const auto *fh_148 = buffer.data(fh + 148);
    const auto *fh_149 = buffer.data(fh + 149);
    const auto *fh_150 = buffer.data(fh + 150);
    const auto *fh_151 = buffer.data(fh + 151);
    const auto *fh_152 = buffer.data(fh + 152);
    const auto *fh_153 = buffer.data(fh + 153);
    const auto *fh_154 = buffer.data(fh + 154);
    const auto *fh_155 = buffer.data(fh + 155);
    const auto *fh_156 = buffer.data(fh + 156);
    const auto *fh_157 = buffer.data(fh + 157);
    const auto *fh_158 = buffer.data(fh + 158);
    const auto *fh_159 = buffer.data(fh + 159);
    const auto *fh_160 = buffer.data(fh + 160);
    const auto *fh_161 = buffer.data(fh + 161);
    const auto *fh_162 = buffer.data(fh + 162);
    const auto *fh_163 = buffer.data(fh + 163);
    const auto *fh_164 = buffer.data(fh + 164);
    const auto *fh_165 = buffer.data(fh + 165);
    const auto *fh_166 = buffer.data(fh + 166);
    const auto *fh_167 = buffer.data(fh + 167);
    const auto *fh_168 = buffer.data(fh + 168);
    const auto *fh_169 = buffer.data(fh + 169);
    const auto *fh_170 = buffer.data(fh + 170);
    const auto *fh_171 = buffer.data(fh + 171);
    const auto *fh_172 = buffer.data(fh + 172);
    const auto *fh_173 = buffer.data(fh + 173);
    const auto *fh_174 = buffer.data(fh + 174);
    const auto *fh_175 = buffer.data(fh + 175);
    const auto *fh_176 = buffer.data(fh + 176);
    const auto *fh_177 = buffer.data(fh + 177);
    const auto *fh_178 = buffer.data(fh + 178);
    const auto *fh_179 = buffer.data(fh + 179);
    const auto *fh_180 = buffer.data(fh + 180);
    const auto *fh_181 = buffer.data(fh + 181);
    const auto *fh_182 = buffer.data(fh + 182);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_184 = buffer.data(fh + 184);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_186 = buffer.data(fh + 186);
    const auto *fh_187 = buffer.data(fh + 187);
    const auto *fh_188 = buffer.data(fh + 188);
    const auto *fh_189 = buffer.data(fh + 189);
    const auto *fh_190 = buffer.data(fh + 190);
    const auto *fh_191 = buffer.data(fh + 191);
    const auto *fh_192 = buffer.data(fh + 192);
    const auto *fh_193 = buffer.data(fh + 193);
    const auto *fh_194 = buffer.data(fh + 194);
    const auto *fh_195 = buffer.data(fh + 195);
    const auto *fh_196 = buffer.data(fh + 196);
    const auto *fh_197 = buffer.data(fh + 197);
    const auto *fh_198 = buffer.data(fh + 198);
    const auto *fh_199 = buffer.data(fh + 199);
    const auto *fh_200 = buffer.data(fh + 200);
    const auto *fh_201 = buffer.data(fh + 201);
    const auto *fh_202 = buffer.data(fh + 202);
    const auto *fh_203 = buffer.data(fh + 203);
    const auto *fh_204 = buffer.data(fh + 204);
    const auto *fh_205 = buffer.data(fh + 205);
    const auto *fh_206 = buffer.data(fh + 206);
    const auto *fh_207 = buffer.data(fh + 207);
    const auto *fh_208 = buffer.data(fh + 208);
    const auto *fh_209 = buffer.data(fh + 209);

#pragma omp simd aligned(fh_22, fh_25, fh_27, fh_32, fh_36, fh_127, fh_130, fh_132, fh_137, \
                         fh_141 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * fh_22[k]
                 - f_1 * fh_27[k]
                 + f_2 * fh_36[k]
                 - f_3 * fh_127[k]
                 + f_4 * fh_132[k]
                 - f_5 * fh_141[k];

        g_1[k] = f_6 * fh_25[k]
                 - f_6 * fh_32[k]
                 - f_7 * fh_130[k]
                 + f_7 * fh_137[k];
    }

#pragma omp simd aligned(fh_22, fh_27, fh_29, fh_36, fh_38, fh_127, fh_132, fh_134, fh_141, \
                         fh_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_8 * fh_22[k]
                 - f_9 * fh_27[k]
                 + f_10 * fh_29[k]
                 + f_11 * fh_36[k]
                 - f_12 * fh_38[k]
                 + f_11 * fh_127[k]
                 + f_13 * fh_132[k]
                 - f_12 * fh_134[k]
                 - f_14 * fh_141[k]
                 + f_15 * fh_143[k];
    }

#pragma omp simd aligned(fh_25, fh_32, fh_34, fh_130, fh_137, fh_139 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_16 * fh_25[k]
                 - f_16 * fh_32[k]
                 + f_17 * fh_34[k]
                 + f_18 * fh_130[k]
                 + f_18 * fh_137[k]
                 - f_19 * fh_139[k];
    }

#pragma omp simd aligned(fh_22, fh_27, fh_29, fh_36, fh_38, fh_40, fh_127, fh_132, fh_134, \
                         fh_141, fh_143, fh_145 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_20 * fh_22[k]
                 + f_21 * fh_27[k]
                 - f_22 * fh_29[k]
                 + f_20 * fh_36[k]
                 - f_22 * fh_38[k]
                 + f_23 * fh_40[k]
                 - f_24 * fh_127[k]
                 - f_25 * fh_132[k]
                 + f_26 * fh_134[k]
                 - f_24 * fh_141[k]
                 + f_26 * fh_143[k]
                 - f_27 * fh_145[k];
    }

#pragma omp simd aligned(fh_23, fh_28, fh_30, fh_37, fh_39, fh_41, fh_128, fh_133, fh_135, \
                         fh_142, fh_144, fh_146 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_28 * fh_23[k]
                 + f_29 * fh_28[k]
                 - f_30 * fh_30[k]
                 + f_28 * fh_37[k]
                 - f_30 * fh_39[k]
                 + f_31 * fh_41[k]
                 - f_32 * fh_128[k]
                 - f_33 * fh_133[k]
                 + f_34 * fh_135[k]
                 - f_32 * fh_142[k]
                 + f_34 * fh_144[k]
                 - f_35 * fh_146[k];
    }

#pragma omp simd aligned(fh_21, fh_24, fh_26, fh_31, fh_33, fh_35, fh_126, fh_129, fh_131, \
                         fh_136, fh_138, fh_140 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_20 * fh_21[k]
                 + f_21 * fh_24[k]
                 - f_22 * fh_26[k]
                 + f_20 * fh_31[k]
                 - f_22 * fh_33[k]
                 + f_23 * fh_35[k]
                 - f_24 * fh_126[k]
                 - f_25 * fh_129[k]
                 + f_26 * fh_131[k]
                 - f_24 * fh_136[k]
                 + f_26 * fh_138[k]
                 - f_27 * fh_140[k];
    }

#pragma omp simd aligned(fh_23, fh_30, fh_37, fh_39, fh_128, fh_135, fh_142, \
                         fh_144 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_36 * fh_23[k]
                 + f_16 * fh_30[k]
                 + f_36 * fh_37[k]
                 - f_16 * fh_39[k]
                 + f_37 * fh_128[k]
                 - f_18 * fh_135[k]
                 - f_37 * fh_142[k]
                 + f_18 * fh_144[k];
    }

#pragma omp simd aligned(fh_21, fh_24, fh_26, fh_31, fh_33, fh_126, fh_129, fh_131, fh_136, \
                         fh_138 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_11 * fh_21[k]
                 + f_9 * fh_24[k]
                 + f_12 * fh_26[k]
                 + f_8 * fh_31[k]
                 - f_10 * fh_33[k]
                 + f_14 * fh_126[k]
                 - f_13 * fh_129[k]
                 - f_15 * fh_131[k]
                 - f_11 * fh_136[k]
                 + f_12 * fh_138[k];
    }

#pragma omp simd aligned(fh_21, fh_23, fh_24, fh_28, fh_31, fh_37, fh_126, fh_128, fh_129, \
                         fh_133, fh_136, fh_142 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_38 * fh_23[k]
                 - f_39 * fh_28[k]
                 + f_38 * fh_37[k]
                 - f_40 * fh_128[k]
                 + f_41 * fh_133[k]
                 - f_40 * fh_142[k];

        g_10[k] = f_2 * fh_21[k]
                  - f_1 * fh_24[k]
                  + f_0 * fh_31[k]
                  - f_5 * fh_126[k]
                  + f_4 * fh_129[k]
                  - f_3 * fh_136[k];
    }

#pragma omp simd aligned(fh_85, fh_88, fh_90, fh_92, fh_95, fh_97, fh_99, fh_101, \
                         fh_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_42 * fh_85[k]
                  - f_43 * fh_90[k]
                  + f_44 * fh_99[k];

        g_12[k] = f_45 * fh_88[k]
                  - f_45 * fh_95[k];

        g_13[k] = -f_36 * fh_85[k]
                  - f_18 * fh_90[k]
                  + f_46 * fh_92[k]
                  + f_37 * fh_99[k]
                  - f_47 * fh_101[k];

        g_14[k] = -f_48 * fh_88[k]
                  - f_48 * fh_95[k]
                  + f_49 * fh_97[k];

        g_15[k] = 1.875 * fh_85[k]
                  + 3.75 * fh_90[k]
                  - 22.5 * fh_92[k]
                  + 1.875 * fh_99[k]
                  - 22.5 * fh_101[k]
                  + 15.0 * fh_103[k];
    }

#pragma omp simd aligned(fh_84, fh_86, fh_87, fh_89, fh_91, fh_93, fh_94, fh_96, fh_98, \
                         fh_100, fh_102, fh_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_50 * fh_86[k]
                  + f_51 * fh_91[k]
                  - f_52 * fh_93[k]
                  + f_50 * fh_100[k]
                  - f_52 * fh_102[k]
                  + f_53 * fh_104[k];

        g_17[k] = 1.875 * fh_84[k]
                  + 3.75 * fh_87[k]
                  - 22.5 * fh_89[k]
                  + 1.875 * fh_94[k]
                  - 22.5 * fh_96[k]
                  + 15.0 * fh_98[k];

        g_18[k] = -f_12 * fh_86[k]
                  + f_48 * fh_93[k]
                  + f_12 * fh_100[k]
                  - f_48 * fh_102[k];

        g_19[k] = -f_37 * fh_84[k]
                  + f_18 * fh_87[k]
                  + f_47 * fh_89[k]
                  + f_36 * fh_94[k]
                  - f_46 * fh_96[k];
    }

#pragma omp simd aligned(fh_84, fh_86, fh_87, fh_91, fh_94, fh_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_54 * fh_86[k]
                  - f_55 * fh_91[k]
                  + f_54 * fh_100[k];

        g_21[k] = f_44 * fh_84[k]
                  - f_43 * fh_87[k]
                  + f_42 * fh_94[k];
    }

#pragma omp simd aligned(fh_22, fh_27, fh_36, fh_127, fh_132, fh_141, fh_169, fh_174, \
                         fh_183 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_56 * fh_22[k]
                  + f_57 * fh_27[k]
                  - f_58 * fh_36[k]
                  - f_56 * fh_127[k]
                  + f_57 * fh_132[k]
                  - f_58 * fh_141[k]
                  + f_54 * fh_169[k]
                  - f_59 * fh_174[k]
                  + f_60 * fh_183[k];
    }

#pragma omp simd aligned(fh_25, fh_32, fh_130, fh_137, fh_172, fh_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_61 * fh_25[k]
                  + f_61 * fh_32[k]
                  - f_61 * fh_130[k]
                  + f_61 * fh_137[k]
                  + f_62 * fh_172[k]
                  - f_62 * fh_179[k];
    }

#pragma omp simd aligned(fh_22, fh_27, fh_29, fh_36, fh_38, fh_127, fh_132, fh_134, fh_141, \
                         fh_143, fh_169, fh_174, fh_176, fh_183, \
                         fh_185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_63 * fh_22[k]
                  + f_64 * fh_27[k]
                  - f_65 * fh_29[k]
                  - f_66 * fh_36[k]
                  + f_67 * fh_38[k]
                  + f_63 * fh_127[k]
                  + f_64 * fh_132[k]
                  - f_65 * fh_134[k]
                  - f_66 * fh_141[k]
                  + f_67 * fh_143[k]
                  - f_68 * fh_169[k]
                  - f_67 * fh_174[k]
                  + f_69 * fh_176[k]
                  + f_70 * fh_183[k]
                  - f_71 * fh_185[k];
    }

#pragma omp simd aligned(fh_25, fh_32, fh_34, fh_130, fh_137, fh_139, fh_172, fh_179, \
                         fh_181 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_72 * fh_25[k]
                  + f_72 * fh_32[k]
                  - f_73 * fh_34[k]
                  + f_72 * fh_130[k]
                  + f_72 * fh_137[k]
                  - f_73 * fh_139[k]
                  - f_74 * fh_172[k]
                  - f_74 * fh_179[k]
                  + f_75 * fh_181[k];
    }

#pragma omp simd aligned(fh_22, fh_27, fh_29, fh_36, fh_38, fh_40, fh_127, fh_132, fh_134, \
                         fh_141, fh_143, fh_145, fh_169, fh_174, fh_176, fh_183, fh_185, \
                         fh_187 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_76 * fh_22[k]
                  - f_77 * fh_27[k]
                  + f_78 * fh_29[k]
                  - f_76 * fh_36[k]
                  + f_78 * fh_38[k]
                  - f_31 * fh_40[k]
                  - f_76 * fh_127[k]
                  - f_77 * fh_132[k]
                  + f_78 * fh_134[k]
                  - f_76 * fh_141[k]
                  + f_78 * fh_143[k]
                  - f_31 * fh_145[k]
                  + f_79 * fh_169[k]
                  + f_31 * fh_174[k]
                  - f_80 * fh_176[k]
                  + f_79 * fh_183[k]
                  - f_80 * fh_185[k]
                  + f_81 * fh_187[k];
    }

#pragma omp simd aligned(fh_23, fh_28, fh_30, fh_37, fh_39, fh_41, fh_128, fh_133, fh_135, \
                         fh_142, fh_144, fh_146, fh_170, fh_175, fh_177, fh_184, fh_186, \
                         fh_188 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_20 * fh_23[k]
                  - f_21 * fh_28[k]
                  + f_27 * fh_30[k]
                  - f_20 * fh_37[k]
                  + f_27 * fh_39[k]
                  - f_82 * fh_41[k]
                  - f_20 * fh_128[k]
                  - f_21 * fh_133[k]
                  + f_27 * fh_135[k]
                  - f_20 * fh_142[k]
                  + f_27 * fh_144[k]
                  - f_82 * fh_146[k]
                  + f_26 * fh_170[k]
                  + f_23 * fh_175[k]
                  - f_83 * fh_177[k]
                  + f_26 * fh_184[k]
                  - f_83 * fh_186[k]
                  + f_84 * fh_188[k];
    }

#pragma omp simd aligned(fh_21, fh_24, fh_26, fh_31, fh_33, fh_35, fh_126, fh_129, fh_131, \
                         fh_136, fh_138, fh_140, fh_168, fh_171, fh_173, fh_178, fh_180, \
                         fh_182 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_76 * fh_21[k]
                  - f_77 * fh_24[k]
                  + f_78 * fh_26[k]
                  - f_76 * fh_31[k]
                  + f_78 * fh_33[k]
                  - f_31 * fh_35[k]
                  - f_76 * fh_126[k]
                  - f_77 * fh_129[k]
                  + f_78 * fh_131[k]
                  - f_76 * fh_136[k]
                  + f_78 * fh_138[k]
                  - f_31 * fh_140[k]
                  + f_79 * fh_168[k]
                  + f_31 * fh_171[k]
                  - f_80 * fh_173[k]
                  + f_79 * fh_178[k]
                  - f_80 * fh_180[k]
                  + f_81 * fh_182[k];
    }

#pragma omp simd aligned(fh_23, fh_30, fh_37, fh_39, fh_128, fh_135, fh_142, fh_144, fh_170, \
                         fh_177, fh_184, fh_186 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_85 * fh_23[k]
                  - f_72 * fh_30[k]
                  - f_85 * fh_37[k]
                  + f_72 * fh_39[k]
                  + f_85 * fh_128[k]
                  - f_72 * fh_135[k]
                  - f_85 * fh_142[k]
                  + f_72 * fh_144[k]
                  - f_73 * fh_170[k]
                  + f_74 * fh_177[k]
                  + f_73 * fh_184[k]
                  - f_74 * fh_186[k];
    }

#pragma omp simd aligned(fh_21, fh_24, fh_26, fh_31, fh_33, fh_126, fh_129, fh_131, fh_136, \
                         fh_138, fh_168, fh_171, fh_173, fh_178, \
                         fh_180 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_66 * fh_21[k]
                  - f_64 * fh_24[k]
                  - f_67 * fh_26[k]
                  - f_63 * fh_31[k]
                  + f_65 * fh_33[k]
                  + f_66 * fh_126[k]
                  - f_64 * fh_129[k]
                  - f_67 * fh_131[k]
                  - f_63 * fh_136[k]
                  + f_65 * fh_138[k]
                  - f_70 * fh_168[k]
                  + f_67 * fh_171[k]
                  + f_71 * fh_173[k]
                  + f_68 * fh_178[k]
                  - f_69 * fh_180[k];
    }

#pragma omp simd aligned(fh_23, fh_28, fh_37, fh_128, fh_133, fh_142, fh_170, fh_175, \
                         fh_184 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_86 * fh_23[k]
                  + f_87 * fh_28[k]
                  - f_86 * fh_37[k]
                  - f_86 * fh_128[k]
                  + f_87 * fh_133[k]
                  - f_86 * fh_142[k]
                  + f_61 * fh_170[k]
                  - f_88 * fh_175[k]
                  + f_61 * fh_184[k];
    }

#pragma omp simd aligned(fh_21, fh_24, fh_31, fh_126, fh_129, fh_136, fh_168, fh_171, \
                         fh_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_58 * fh_21[k]
                  + f_57 * fh_24[k]
                  - f_56 * fh_31[k]
                  - f_58 * fh_126[k]
                  + f_57 * fh_129[k]
                  - f_56 * fh_136[k]
                  + f_60 * fh_168[k]
                  - f_59 * fh_171[k]
                  + f_54 * fh_178[k];
    }

#pragma omp simd aligned(fh_43, fh_48, fh_57, fh_148, fh_153, fh_162, fh_190, fh_195, \
                         fh_204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_38 * fh_43[k]
                  + f_41 * fh_48[k]
                  - f_89 * fh_57[k]
                  - f_38 * fh_148[k]
                  + f_41 * fh_153[k]
                  - f_89 * fh_162[k]
                  + f_90 * fh_190[k]
                  - f_7 * fh_195[k]
                  + f_91 * fh_204[k];
    }

#pragma omp simd aligned(fh_46, fh_53, fh_151, fh_158, fh_193, fh_200 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_92 * fh_46[k]
                  + f_92 * fh_53[k]
                  - f_92 * fh_151[k]
                  + f_92 * fh_158[k]
                  + f_93 * fh_193[k]
                  - f_93 * fh_200[k];
    }

#pragma omp simd aligned(fh_43, fh_48, fh_50, fh_57, fh_59, fh_148, fh_153, fh_155, fh_162, \
                         fh_164, fh_190, fh_195, fh_197, fh_204, \
                         fh_206 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_94 * fh_43[k]
                  + f_85 * fh_48[k]
                  - f_95 * fh_50[k]
                  - f_96 * fh_57[k]
                  + f_73 * fh_59[k]
                  + f_94 * fh_148[k]
                  + f_85 * fh_153[k]
                  - f_95 * fh_155[k]
                  - f_96 * fh_162[k]
                  + f_73 * fh_164[k]
                  - f_85 * fh_190[k]
                  - f_97 * fh_195[k]
                  + f_74 * fh_197[k]
                  + f_98 * fh_204[k]
                  - f_99 * fh_206[k];
    }

#pragma omp simd aligned(fh_46, fh_53, fh_55, fh_151, fh_158, fh_160, fh_193, fh_200, \
                         fh_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_65 * fh_46[k]
                  + f_65 * fh_53[k]
                  - f_100 * fh_55[k]
                  + f_65 * fh_151[k]
                  + f_65 * fh_158[k]
                  - f_100 * fh_160[k]
                  - f_101 * fh_193[k]
                  - f_101 * fh_200[k]
                  + f_71 * fh_202[k];
    }

#pragma omp simd aligned(fh_43, fh_48, fh_50, fh_57, fh_59, fh_61, fh_148, fh_153, fh_155, \
                         fh_162, fh_164, fh_166, fh_190, fh_195, fh_197, fh_204, fh_206, \
                         fh_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_102 * fh_43[k]
                  - f_103 * fh_48[k]
                  + f_104 * fh_50[k]
                  - f_102 * fh_57[k]
                  + f_104 * fh_59[k]
                  - f_105 * fh_61[k]
                  - f_102 * fh_148[k]
                  - f_103 * fh_153[k]
                  + f_104 * fh_155[k]
                  - f_102 * fh_162[k]
                  + f_104 * fh_164[k]
                  - f_105 * fh_166[k]
                  + f_106 * fh_190[k]
                  + f_107 * fh_195[k]
                  - f_105 * fh_197[k]
                  + f_106 * fh_204[k]
                  - f_105 * fh_206[k]
                  + f_53 * fh_208[k];
    }

#pragma omp simd aligned(fh_44, fh_49, fh_51, fh_58, fh_60, fh_62, fh_149, fh_154, fh_156, \
                         fh_163, fh_165, fh_167, fh_191, fh_196, fh_198, fh_205, fh_207, \
                         fh_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -2.8125 * fh_44[k]
                  - 5.625 * fh_49[k]
                  + 7.5 * fh_51[k]
                  - 2.8125 * fh_58[k]
                  + 7.5 * fh_60[k]
                  - 1.5 * fh_62[k]
                  - 2.8125 * fh_149[k]
                  - 5.625 * fh_154[k]
                  + 7.5 * fh_156[k]
                  - 2.8125 * fh_163[k]
                  + 7.5 * fh_165[k]
                  - 1.5 * fh_167[k]
                  + 1.875 * fh_191[k]
                  + 3.75 * fh_196[k]
                  - 5.0 * fh_198[k]
                  + 1.875 * fh_205[k]
                  - 5.0 * fh_207[k]
                  + fh_209[k];
    }

#pragma omp simd aligned(fh_42, fh_45, fh_47, fh_52, fh_54, fh_56, fh_147, fh_150, fh_152, \
                         fh_157, fh_159, fh_161, fh_189, fh_192, fh_194, fh_199, fh_201, \
                         fh_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_102 * fh_42[k]
                  - f_103 * fh_45[k]
                  + f_104 * fh_47[k]
                  - f_102 * fh_52[k]
                  + f_104 * fh_54[k]
                  - f_105 * fh_56[k]
                  - f_102 * fh_147[k]
                  - f_103 * fh_150[k]
                  + f_104 * fh_152[k]
                  - f_102 * fh_157[k]
                  + f_104 * fh_159[k]
                  - f_105 * fh_161[k]
                  + f_106 * fh_189[k]
                  + f_107 * fh_192[k]
                  - f_105 * fh_194[k]
                  + f_106 * fh_199[k]
                  - f_105 * fh_201[k]
                  + f_53 * fh_203[k];
    }

#pragma omp simd aligned(fh_44, fh_51, fh_58, fh_60, fh_149, fh_156, fh_163, fh_165, fh_191, \
                         fh_198, fh_205, fh_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_68 * fh_44[k]
                  - f_65 * fh_51[k]
                  - f_68 * fh_58[k]
                  + f_65 * fh_60[k]
                  + f_68 * fh_149[k]
                  - f_65 * fh_156[k]
                  - f_68 * fh_163[k]
                  + f_65 * fh_165[k]
                  - f_67 * fh_191[k]
                  + f_101 * fh_198[k]
                  + f_67 * fh_205[k]
                  - f_101 * fh_207[k];
    }

#pragma omp simd aligned(fh_42, fh_45, fh_47, fh_52, fh_54, fh_147, fh_150, fh_152, fh_157, \
                         fh_159, fh_189, fh_192, fh_194, fh_199, \
                         fh_201 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_96 * fh_42[k]
                  - f_85 * fh_45[k]
                  - f_73 * fh_47[k]
                  - f_94 * fh_52[k]
                  + f_95 * fh_54[k]
                  + f_96 * fh_147[k]
                  - f_85 * fh_150[k]
                  - f_73 * fh_152[k]
                  - f_94 * fh_157[k]
                  + f_95 * fh_159[k]
                  - f_98 * fh_189[k]
                  + f_97 * fh_192[k]
                  + f_99 * fh_194[k]
                  + f_85 * fh_199[k]
                  - f_74 * fh_201[k];
    }

#pragma omp simd aligned(fh_44, fh_49, fh_58, fh_149, fh_154, fh_163, fh_191, fh_196, \
                         fh_205 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_108 * fh_44[k]
                  + f_109 * fh_49[k]
                  - f_108 * fh_58[k]
                  - f_108 * fh_149[k]
                  + f_109 * fh_154[k]
                  - f_108 * fh_163[k]
                  + f_110 * fh_191[k]
                  - f_92 * fh_196[k]
                  + f_110 * fh_205[k];
    }

#pragma omp simd aligned(fh_42, fh_45, fh_52, fh_147, fh_150, fh_157, fh_189, fh_192, \
                         fh_199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_89 * fh_42[k]
                  + f_41 * fh_45[k]
                  - f_38 * fh_52[k]
                  - f_89 * fh_147[k]
                  + f_41 * fh_150[k]
                  - f_38 * fh_157[k]
                  + f_91 * fh_189[k]
                  - f_7 * fh_192[k]
                  + f_90 * fh_199[k];
    }

#pragma omp simd aligned(fh_1, fh_6, fh_15, fh_64, fh_69, fh_78, fh_106, fh_111, \
                         fh_120 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_56 * fh_1[k]
                  + f_57 * fh_6[k]
                  - f_58 * fh_15[k]
                  - f_56 * fh_64[k]
                  + f_57 * fh_69[k]
                  - f_58 * fh_78[k]
                  + f_54 * fh_106[k]
                  - f_59 * fh_111[k]
                  + f_60 * fh_120[k];
    }

#pragma omp simd aligned(fh_4, fh_11, fh_67, fh_74, fh_109, fh_116 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_61 * fh_4[k]
                  + f_61 * fh_11[k]
                  - f_61 * fh_67[k]
                  + f_61 * fh_74[k]
                  + f_62 * fh_109[k]
                  - f_62 * fh_116[k];
    }

#pragma omp simd aligned(fh_1, fh_6, fh_8, fh_15, fh_17, fh_64, fh_69, fh_71, fh_78, fh_80, \
                         fh_106, fh_111, fh_113, fh_120, fh_122 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_63 * fh_1[k]
                  + f_64 * fh_6[k]
                  - f_65 * fh_8[k]
                  - f_66 * fh_15[k]
                  + f_67 * fh_17[k]
                  + f_63 * fh_64[k]
                  + f_64 * fh_69[k]
                  - f_65 * fh_71[k]
                  - f_66 * fh_78[k]
                  + f_67 * fh_80[k]
                  - f_68 * fh_106[k]
                  - f_67 * fh_111[k]
                  + f_69 * fh_113[k]
                  + f_70 * fh_120[k]
                  - f_71 * fh_122[k];
    }

#pragma omp simd aligned(fh_4, fh_11, fh_13, fh_67, fh_74, fh_76, fh_109, fh_116, \
                         fh_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = f_72 * fh_4[k]
                  + f_72 * fh_11[k]
                  - f_73 * fh_13[k]
                  + f_72 * fh_67[k]
                  + f_72 * fh_74[k]
                  - f_73 * fh_76[k]
                  - f_74 * fh_109[k]
                  - f_74 * fh_116[k]
                  + f_75 * fh_118[k];
    }

#pragma omp simd aligned(fh_1, fh_6, fh_8, fh_15, fh_17, fh_19, fh_64, fh_69, fh_71, fh_78, \
                         fh_80, fh_82, fh_106, fh_111, fh_113, fh_120, fh_122, \
                         fh_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_76 * fh_1[k]
                  - f_77 * fh_6[k]
                  + f_78 * fh_8[k]
                  - f_76 * fh_15[k]
                  + f_78 * fh_17[k]
                  - f_31 * fh_19[k]
                  - f_76 * fh_64[k]
                  - f_77 * fh_69[k]
                  + f_78 * fh_71[k]
                  - f_76 * fh_78[k]
                  + f_78 * fh_80[k]
                  - f_31 * fh_82[k]
                  + f_79 * fh_106[k]
                  + f_31 * fh_111[k]
                  - f_80 * fh_113[k]
                  + f_79 * fh_120[k]
                  - f_80 * fh_122[k]
                  + f_81 * fh_124[k];
    }

#pragma omp simd aligned(fh_2, fh_7, fh_9, fh_16, fh_18, fh_20, fh_65, fh_70, fh_72, fh_79, \
                         fh_81, fh_83, fh_107, fh_112, fh_114, fh_121, fh_123, \
                         fh_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_20 * fh_2[k]
                  - f_21 * fh_7[k]
                  + f_27 * fh_9[k]
                  - f_20 * fh_16[k]
                  + f_27 * fh_18[k]
                  - f_82 * fh_20[k]
                  - f_20 * fh_65[k]
                  - f_21 * fh_70[k]
                  + f_27 * fh_72[k]
                  - f_20 * fh_79[k]
                  + f_27 * fh_81[k]
                  - f_82 * fh_83[k]
                  + f_26 * fh_107[k]
                  + f_23 * fh_112[k]
                  - f_83 * fh_114[k]
                  + f_26 * fh_121[k]
                  - f_83 * fh_123[k]
                  + f_84 * fh_125[k];
    }

#pragma omp simd aligned(fh_0, fh_3, fh_5, fh_10, fh_12, fh_14, fh_63, fh_66, fh_68, fh_73, \
                         fh_75, fh_77, fh_105, fh_108, fh_110, fh_115, fh_117, \
                         fh_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_76 * fh_0[k]
                  - f_77 * fh_3[k]
                  + f_78 * fh_5[k]
                  - f_76 * fh_10[k]
                  + f_78 * fh_12[k]
                  - f_31 * fh_14[k]
                  - f_76 * fh_63[k]
                  - f_77 * fh_66[k]
                  + f_78 * fh_68[k]
                  - f_76 * fh_73[k]
                  + f_78 * fh_75[k]
                  - f_31 * fh_77[k]
                  + f_79 * fh_105[k]
                  + f_31 * fh_108[k]
                  - f_80 * fh_110[k]
                  + f_79 * fh_115[k]
                  - f_80 * fh_117[k]
                  + f_81 * fh_119[k];
    }

#pragma omp simd aligned(fh_2, fh_9, fh_16, fh_18, fh_65, fh_72, fh_79, fh_81, fh_107, fh_114, \
                         fh_121, fh_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_85 * fh_2[k]
                  - f_72 * fh_9[k]
                  - f_85 * fh_16[k]
                  + f_72 * fh_18[k]
                  + f_85 * fh_65[k]
                  - f_72 * fh_72[k]
                  - f_85 * fh_79[k]
                  + f_72 * fh_81[k]
                  - f_73 * fh_107[k]
                  + f_74 * fh_114[k]
                  + f_73 * fh_121[k]
                  - f_74 * fh_123[k];
    }

#pragma omp simd aligned(fh_0, fh_3, fh_5, fh_10, fh_12, fh_63, fh_66, fh_68, fh_73, fh_75, \
                         fh_105, fh_108, fh_110, fh_115, fh_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_66 * fh_0[k]
                  - f_64 * fh_3[k]
                  - f_67 * fh_5[k]
                  - f_63 * fh_10[k]
                  + f_65 * fh_12[k]
                  + f_66 * fh_63[k]
                  - f_64 * fh_66[k]
                  - f_67 * fh_68[k]
                  - f_63 * fh_73[k]
                  + f_65 * fh_75[k]
                  - f_70 * fh_105[k]
                  + f_67 * fh_108[k]
                  + f_71 * fh_110[k]
                  + f_68 * fh_115[k]
                  - f_69 * fh_117[k];
    }

#pragma omp simd aligned(fh_2, fh_7, fh_16, fh_65, fh_70, fh_79, fh_107, fh_112, \
                         fh_121 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = -f_86 * fh_2[k]
                  + f_87 * fh_7[k]
                  - f_86 * fh_16[k]
                  - f_86 * fh_65[k]
                  + f_87 * fh_70[k]
                  - f_86 * fh_79[k]
                  + f_61 * fh_107[k]
                  - f_88 * fh_112[k]
                  + f_61 * fh_121[k];
    }

#pragma omp simd aligned(fh_0, fh_3, fh_10, fh_63, fh_66, fh_73, fh_105, fh_108, \
                         fh_115 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_58 * fh_0[k]
                  + f_57 * fh_3[k]
                  - f_56 * fh_10[k]
                  - f_58 * fh_63[k]
                  + f_57 * fh_66[k]
                  - f_56 * fh_73[k]
                  + f_60 * fh_105[k]
                  - f_59 * fh_108[k]
                  + f_54 * fh_115[k];
    }

#pragma omp simd aligned(fh_43, fh_46, fh_48, fh_53, fh_57, fh_148, fh_151, fh_153, fh_158, \
                         fh_162 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_111 * fh_43[k]
                  - f_42 * fh_48[k]
                  + f_86 * fh_57[k]
                  - f_111 * fh_148[k]
                  + f_42 * fh_153[k]
                  - f_86 * fh_162[k];

        g_56[k] = f_59 * fh_46[k]
                  - f_59 * fh_53[k]
                  - f_59 * fh_151[k]
                  + f_59 * fh_158[k];
    }

#pragma omp simd aligned(fh_43, fh_48, fh_50, fh_57, fh_59, fh_148, fh_153, fh_155, fh_162, \
                         fh_164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = -f_112 * fh_43[k]
                  - f_37 * fh_48[k]
                  + f_17 * fh_50[k]
                  + f_113 * fh_57[k]
                  - f_19 * fh_59[k]
                  + f_112 * fh_148[k]
                  + f_37 * fh_153[k]
                  - f_17 * fh_155[k]
                  - f_113 * fh_162[k]
                  + f_19 * fh_164[k];
    }

#pragma omp simd aligned(fh_46, fh_53, fh_55, fh_151, fh_158, fh_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_12 * fh_46[k]
                  - f_12 * fh_53[k]
                  + f_48 * fh_55[k]
                  + f_12 * fh_151[k]
                  + f_12 * fh_158[k]
                  - f_48 * fh_160[k];
    }

#pragma omp simd aligned(fh_43, fh_48, fh_50, fh_57, fh_59, fh_61, fh_148, fh_153, fh_155, \
                         fh_162, fh_164, fh_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = 0.9375 * fh_43[k]
                  + 1.875 * fh_48[k]
                  - 11.25 * fh_50[k]
                  + 0.9375 * fh_57[k]
                  - 11.25 * fh_59[k]
                  + 7.5 * fh_61[k]
                  - 0.9375 * fh_148[k]
                  - 1.875 * fh_153[k]
                  + 11.25 * fh_155[k]
                  - 0.9375 * fh_162[k]
                  + 11.25 * fh_164[k]
                  - 7.5 * fh_166[k];
    }

#pragma omp simd aligned(fh_44, fh_49, fh_51, fh_58, fh_60, fh_62, fh_149, fh_154, fh_156, \
                         fh_163, fh_165, fh_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_114 * fh_44[k]
                  + f_50 * fh_49[k]
                  - f_115 * fh_51[k]
                  + f_114 * fh_58[k]
                  - f_115 * fh_60[k]
                  + f_116 * fh_62[k]
                  - f_114 * fh_149[k]
                  - f_50 * fh_154[k]
                  + f_115 * fh_156[k]
                  - f_114 * fh_163[k]
                  + f_115 * fh_165[k]
                  - f_116 * fh_167[k];
    }

#pragma omp simd aligned(fh_42, fh_45, fh_47, fh_52, fh_54, fh_56, fh_147, fh_150, fh_152, \
                         fh_157, fh_159, fh_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = 0.9375 * fh_42[k]
                  + 1.875 * fh_45[k]
                  - 11.25 * fh_47[k]
                  + 0.9375 * fh_52[k]
                  - 11.25 * fh_54[k]
                  + 7.5 * fh_56[k]
                  - 0.9375 * fh_147[k]
                  - 1.875 * fh_150[k]
                  + 11.25 * fh_152[k]
                  - 0.9375 * fh_157[k]
                  + 11.25 * fh_159[k]
                  - 7.5 * fh_161[k];
    }

#pragma omp simd aligned(fh_44, fh_51, fh_58, fh_60, fh_149, fh_156, fh_163, \
                         fh_165 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_117 * fh_44[k]
                  + f_12 * fh_51[k]
                  + f_117 * fh_58[k]
                  - f_12 * fh_60[k]
                  + f_117 * fh_149[k]
                  - f_12 * fh_156[k]
                  - f_117 * fh_163[k]
                  + f_12 * fh_165[k];
    }

#pragma omp simd aligned(fh_42, fh_45, fh_47, fh_52, fh_54, fh_147, fh_150, fh_152, fh_157, \
                         fh_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_113 * fh_42[k]
                  + f_37 * fh_45[k]
                  + f_19 * fh_47[k]
                  + f_112 * fh_52[k]
                  - f_17 * fh_54[k]
                  + f_113 * fh_147[k]
                  - f_37 * fh_150[k]
                  - f_19 * fh_152[k]
                  - f_112 * fh_157[k]
                  + f_17 * fh_159[k];
    }

#pragma omp simd aligned(fh_42, fh_44, fh_45, fh_49, fh_52, fh_58, fh_147, fh_149, fh_150, \
                         fh_154, fh_157, fh_163 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_57 * fh_44[k]
                  - f_118 * fh_49[k]
                  + f_57 * fh_58[k]
                  - f_57 * fh_149[k]
                  + f_118 * fh_154[k]
                  - f_57 * fh_163[k];

        g_65[k] = f_86 * fh_42[k]
                  - f_42 * fh_45[k]
                  + f_111 * fh_52[k]
                  - f_86 * fh_147[k]
                  + f_42 * fh_150[k]
                  - f_111 * fh_157[k];
    }

#pragma omp simd aligned(fh_1, fh_4, fh_6, fh_11, fh_15, fh_64, fh_67, fh_69, fh_74, \
                         fh_78 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = f_3 * fh_1[k]
                  - f_4 * fh_6[k]
                  + f_5 * fh_15[k]
                  - f_0 * fh_64[k]
                  + f_1 * fh_69[k]
                  - f_2 * fh_78[k];

        g_67[k] = f_7 * fh_4[k]
                  - f_7 * fh_11[k]
                  - f_6 * fh_67[k]
                  + f_6 * fh_74[k];
    }

#pragma omp simd aligned(fh_1, fh_6, fh_8, fh_15, fh_17, fh_64, fh_69, fh_71, fh_78, \
                         fh_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = -f_11 * fh_1[k]
                  - f_13 * fh_6[k]
                  + f_12 * fh_8[k]
                  + f_14 * fh_15[k]
                  - f_15 * fh_17[k]
                  + f_8 * fh_64[k]
                  + f_9 * fh_69[k]
                  - f_10 * fh_71[k]
                  - f_11 * fh_78[k]
                  + f_12 * fh_80[k];
    }

#pragma omp simd aligned(fh_4, fh_11, fh_13, fh_67, fh_74, fh_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = -f_18 * fh_4[k]
                  - f_18 * fh_11[k]
                  + f_19 * fh_13[k]
                  + f_16 * fh_67[k]
                  + f_16 * fh_74[k]
                  - f_17 * fh_76[k];
    }

#pragma omp simd aligned(fh_1, fh_6, fh_8, fh_15, fh_17, fh_19, fh_64, fh_69, fh_71, fh_78, \
                         fh_80, fh_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = f_24 * fh_1[k]
                  + f_25 * fh_6[k]
                  - f_26 * fh_8[k]
                  + f_24 * fh_15[k]
                  - f_26 * fh_17[k]
                  + f_27 * fh_19[k]
                  - f_20 * fh_64[k]
                  - f_21 * fh_69[k]
                  + f_22 * fh_71[k]
                  - f_20 * fh_78[k]
                  + f_22 * fh_80[k]
                  - f_23 * fh_82[k];
    }

#pragma omp simd aligned(fh_2, fh_7, fh_9, fh_16, fh_18, fh_20, fh_65, fh_70, fh_72, fh_79, \
                         fh_81, fh_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = f_32 * fh_2[k]
                  + f_33 * fh_7[k]
                  - f_34 * fh_9[k]
                  + f_32 * fh_16[k]
                  - f_34 * fh_18[k]
                  + f_35 * fh_20[k]
                  - f_28 * fh_65[k]
                  - f_29 * fh_70[k]
                  + f_30 * fh_72[k]
                  - f_28 * fh_79[k]
                  + f_30 * fh_81[k]
                  - f_31 * fh_83[k];
    }

#pragma omp simd aligned(fh_0, fh_3, fh_5, fh_10, fh_12, fh_14, fh_63, fh_66, fh_68, fh_73, \
                         fh_75, fh_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = f_24 * fh_0[k]
                  + f_25 * fh_3[k]
                  - f_26 * fh_5[k]
                  + f_24 * fh_10[k]
                  - f_26 * fh_12[k]
                  + f_27 * fh_14[k]
                  - f_20 * fh_63[k]
                  - f_21 * fh_66[k]
                  + f_22 * fh_68[k]
                  - f_20 * fh_73[k]
                  + f_22 * fh_75[k]
                  - f_23 * fh_77[k];
    }

#pragma omp simd aligned(fh_2, fh_9, fh_16, fh_18, fh_65, fh_72, fh_79, \
                         fh_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = -f_37 * fh_2[k]
                  + f_18 * fh_9[k]
                  + f_37 * fh_16[k]
                  - f_18 * fh_18[k]
                  + f_36 * fh_65[k]
                  - f_16 * fh_72[k]
                  - f_36 * fh_79[k]
                  + f_16 * fh_81[k];
    }

#pragma omp simd aligned(fh_0, fh_3, fh_5, fh_10, fh_12, fh_63, fh_66, fh_68, fh_73, \
                         fh_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_14 * fh_0[k]
                  + f_13 * fh_3[k]
                  + f_15 * fh_5[k]
                  + f_11 * fh_10[k]
                  - f_12 * fh_12[k]
                  + f_11 * fh_63[k]
                  - f_9 * fh_66[k]
                  - f_12 * fh_68[k]
                  - f_8 * fh_73[k]
                  + f_10 * fh_75[k];
    }

#pragma omp simd aligned(fh_0, fh_2, fh_3, fh_7, fh_10, fh_16, fh_63, fh_65, fh_66, fh_70, \
                         fh_73, fh_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = f_40 * fh_2[k]
                  - f_41 * fh_7[k]
                  + f_40 * fh_16[k]
                  - f_38 * fh_65[k]
                  + f_39 * fh_70[k]
                  - f_38 * fh_79[k];

        g_76[k] = f_5 * fh_0[k]
                  - f_4 * fh_3[k]
                  + f_3 * fh_10[k]
                  - f_2 * fh_63[k]
                  + f_1 * fh_66[k]
                  - f_0 * fh_73[k];
    }
}

}  // namespace simdtrf
