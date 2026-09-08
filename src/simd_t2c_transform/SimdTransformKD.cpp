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


#include "SimdTransformKD.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_kd(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t kd,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.65625 * std::sqrt(143.0);
    const auto f_1 = 3.28125 * std::sqrt(143.0);
    const auto f_2 = 1.96875 * std::sqrt(143.0);
    const auto f_3 = 0.09375 * std::sqrt(143.0);
    const auto f_4 = 0.109375 * std::sqrt(429.0);
    const auto f_5 = 0.21875 * std::sqrt(429.0);
    const auto f_6 = 0.546875 * std::sqrt(429.0);
    const auto f_7 = 1.09375 * std::sqrt(429.0);
    const auto f_8 = 0.328125 * std::sqrt(429.0);
    const auto f_9 = 0.65625 * std::sqrt(429.0);
    const auto f_10 = 0.015625 * std::sqrt(429.0);
    const auto f_11 = 0.03125 * std::sqrt(429.0);
    const auto f_12 = 0.328125 * std::sqrt(143.0);
    const auto f_13 = 1.640625 * std::sqrt(143.0);
    const auto f_14 = 0.984375 * std::sqrt(143.0);
    const auto f_15 = 0.046875 * std::sqrt(143.0);
    const auto f_16 = 0.5625 * std::sqrt(2002.0);
    const auto f_17 = 1.875 * std::sqrt(2002.0);
    const auto f_18 = 0.09375 * std::sqrt(6006.0);
    const auto f_19 = 0.1875 * std::sqrt(6006.0);
    const auto f_20 = 0.3125 * std::sqrt(6006.0);
    const auto f_21 = 0.625 * std::sqrt(6006.0);
    const auto f_22 = 0.28125 * std::sqrt(2002.0);
    const auto f_23 = 0.9375 * std::sqrt(2002.0);
    const auto f_24 = 0.46875 * std::sqrt(77.0);
    const auto f_25 = 5.625 * std::sqrt(77.0);
    const auto f_26 = 0.84375 * std::sqrt(77.0);
    const auto f_27 = 11.25 * std::sqrt(77.0);
    const auto f_28 = 0.09375 * std::sqrt(77.0);
    const auto f_29 = 1.125 * std::sqrt(77.0);
    const auto f_30 = 0.078125 * std::sqrt(231.0);
    const auto f_31 = 0.15625 * std::sqrt(231.0);
    const auto f_32 = 0.9375 * std::sqrt(231.0);
    const auto f_33 = 1.875 * std::sqrt(231.0);
    const auto f_34 = 0.140625 * std::sqrt(231.0);
    const auto f_35 = 0.28125 * std::sqrt(231.0);
    const auto f_36 = 3.75 * std::sqrt(231.0);
    const auto f_37 = 0.015625 * std::sqrt(231.0);
    const auto f_38 = 0.03125 * std::sqrt(231.0);
    const auto f_39 = 0.1875 * std::sqrt(231.0);
    const auto f_40 = 0.375 * std::sqrt(231.0);
    const auto f_41 = 0.234375 * std::sqrt(77.0);
    const auto f_42 = 2.8125 * std::sqrt(77.0);
    const auto f_43 = 0.421875 * std::sqrt(77.0);
    const auto f_44 = 0.046875 * std::sqrt(77.0);
    const auto f_45 = 0.5625 * std::sqrt(77.0);
    const auto f_46 = 2.25 * std::sqrt(77.0);
    const auto f_47 = 7.5 * std::sqrt(77.0);
    const auto f_48 = 0.75 * std::sqrt(231.0);
    const auto f_49 = 1.25 * std::sqrt(231.0);
    const auto f_50 = 2.5 * std::sqrt(231.0);
    const auto f_51 = 3.75 * std::sqrt(77.0);
    const auto f_52 = 0.84375 * std::sqrt(7.0);
    const auto f_53 = 1.40625 * std::sqrt(7.0);
    const auto f_54 = 16.875 * std::sqrt(7.0);
    const auto f_55 = 0.28125 * std::sqrt(7.0);
    const auto f_56 = 11.25 * std::sqrt(7.0);
    const auto f_57 = 22.5 * std::sqrt(7.0);
    const auto f_58 = 5.625 * std::sqrt(7.0);
    const auto f_59 = 7.5 * std::sqrt(7.0);
    const auto f_60 = 0.140625 * std::sqrt(21.0);
    const auto f_61 = 0.28125 * std::sqrt(21.0);
    const auto f_62 = 0.234375 * std::sqrt(21.0);
    const auto f_63 = 0.46875 * std::sqrt(21.0);
    const auto f_64 = 2.8125 * std::sqrt(21.0);
    const auto f_65 = 5.625 * std::sqrt(21.0);
    const auto f_66 = 0.046875 * std::sqrt(21.0);
    const auto f_67 = 0.09375 * std::sqrt(21.0);
    const auto f_68 = 1.875 * std::sqrt(21.0);
    const auto f_69 = 3.75 * std::sqrt(21.0);
    const auto f_70 = 7.5 * std::sqrt(21.0);
    const auto f_71 = 0.9375 * std::sqrt(21.0);
    const auto f_72 = 1.25 * std::sqrt(21.0);
    const auto f_73 = 2.5 * std::sqrt(21.0);
    const auto f_74 = 0.421875 * std::sqrt(7.0);
    const auto f_75 = 0.703125 * std::sqrt(7.0);
    const auto f_76 = 8.4375 * std::sqrt(7.0);
    const auto f_77 = 0.140625 * std::sqrt(7.0);
    const auto f_78 = 2.8125 * std::sqrt(7.0);
    const auto f_79 = 3.75 * std::sqrt(7.0);
    const auto f_80 = 2.8125 * std::sqrt(14.0);
    const auto f_81 = 5.625 * std::sqrt(14.0);
    const auto f_82 = 15.0 * std::sqrt(14.0);
    const auto f_83 = 9.0 * std::sqrt(14.0);
    const auto f_84 = 0.46875 * std::sqrt(42.0);
    const auto f_85 = 0.9375 * std::sqrt(42.0);
    const auto f_86 = 1.875 * std::sqrt(42.0);
    const auto f_87 = 2.5 * std::sqrt(42.0);
    const auto f_88 = 5.0 * std::sqrt(42.0);
    const auto f_89 = 1.5 * std::sqrt(42.0);
    const auto f_90 = 3.0 * std::sqrt(42.0);
    const auto f_91 = 1.40625 * std::sqrt(14.0);
    const auto f_92 = 7.5 * std::sqrt(14.0);
    const auto f_93 = 4.5 * std::sqrt(14.0);
    const auto f_94 = 0.15625 * std::sqrt(21.0);
    const auto f_95 = 2.0 * std::sqrt(21.0);
    const auto f_96 = 0.078125 * std::sqrt(7.0);
    const auto f_97 = 0.15625 * std::sqrt(7.0);
    const auto f_98 = 0.234375 * std::sqrt(7.0);
    const auto f_99 = 0.46875 * std::sqrt(7.0);
    const auto f_100 = 1.875 * std::sqrt(7.0);
    const auto f_101 = std::sqrt(7.0);
    const auto f_102 = 2.0 * std::sqrt(7.0);
    const auto f_103 = 0.078125 * std::sqrt(21.0);
    const auto f_104 = std::sqrt(21.0);
    const auto f_105 = 2.1875 * std::sqrt(3.0);
    const auto f_106 = 6.5625 * std::sqrt(3.0);
    const auto f_107 = 13.125 * std::sqrt(3.0);
    const auto f_108 = 26.25 * std::sqrt(3.0);
    const auto f_109 = 10.5 * std::sqrt(3.0);
    const auto f_110 = std::sqrt(3.0);
    const auto f_111 = 1.09375 * std::sqrt(3.0);
    const auto f_112 = 3.28125 * std::sqrt(3.0);
    const auto f_113 = 5.25 * std::sqrt(3.0);
    const auto f_114 = 0.5 * std::sqrt(3.0);
    const auto f_115 = 0.234375 * std::sqrt(42.0);
    const auto f_116 = 1.25 * std::sqrt(42.0);
    const auto f_117 = 0.75 * std::sqrt(42.0);
    const auto f_118 = 0.703125 * std::sqrt(14.0);
    const auto f_119 = 3.75 * std::sqrt(14.0);
    const auto f_120 = 2.25 * std::sqrt(14.0);
    const auto f_121 = 1.875 * std::sqrt(77.0);
    const auto f_122 = 0.09375 * std::sqrt(231.0);
    const auto f_123 = 0.46875 * std::sqrt(231.0);
    const auto f_124 = 0.3125 * std::sqrt(231.0);
    const auto f_125 = 0.625 * std::sqrt(231.0);
    const auto f_126 = 0.28125 * std::sqrt(77.0);
    const auto f_127 = 1.40625 * std::sqrt(77.0);
    const auto f_128 = 0.9375 * std::sqrt(77.0);
    const auto f_129 = 0.09375 * std::sqrt(2002.0);
    const auto f_130 = 1.40625 * std::sqrt(2002.0);
    const auto f_131 = 0.015625 * std::sqrt(6006.0);
    const auto f_132 = 0.03125 * std::sqrt(6006.0);
    const auto f_133 = 0.234375 * std::sqrt(6006.0);
    const auto f_134 = 0.46875 * std::sqrt(6006.0);
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

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_79 = buffer.data(kd + 79);
    const auto *kd_80 = buffer.data(kd + 80);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_82 = buffer.data(kd + 82);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_85 = buffer.data(kd + 85);
    const auto *kd_86 = buffer.data(kd + 86);
    const auto *kd_87 = buffer.data(kd + 87);
    const auto *kd_88 = buffer.data(kd + 88);
    const auto *kd_89 = buffer.data(kd + 89);
    const auto *kd_90 = buffer.data(kd + 90);
    const auto *kd_91 = buffer.data(kd + 91);
    const auto *kd_92 = buffer.data(kd + 92);
    const auto *kd_93 = buffer.data(kd + 93);
    const auto *kd_94 = buffer.data(kd + 94);
    const auto *kd_95 = buffer.data(kd + 95);
    const auto *kd_96 = buffer.data(kd + 96);
    const auto *kd_97 = buffer.data(kd + 97);
    const auto *kd_98 = buffer.data(kd + 98);
    const auto *kd_99 = buffer.data(kd + 99);
    const auto *kd_100 = buffer.data(kd + 100);
    const auto *kd_101 = buffer.data(kd + 101);
    const auto *kd_102 = buffer.data(kd + 102);
    const auto *kd_103 = buffer.data(kd + 103);
    const auto *kd_104 = buffer.data(kd + 104);
    const auto *kd_105 = buffer.data(kd + 105);
    const auto *kd_106 = buffer.data(kd + 106);
    const auto *kd_107 = buffer.data(kd + 107);
    const auto *kd_108 = buffer.data(kd + 108);
    const auto *kd_109 = buffer.data(kd + 109);
    const auto *kd_110 = buffer.data(kd + 110);
    const auto *kd_111 = buffer.data(kd + 111);
    const auto *kd_112 = buffer.data(kd + 112);
    const auto *kd_113 = buffer.data(kd + 113);
    const auto *kd_114 = buffer.data(kd + 114);
    const auto *kd_115 = buffer.data(kd + 115);
    const auto *kd_116 = buffer.data(kd + 116);
    const auto *kd_117 = buffer.data(kd + 117);
    const auto *kd_118 = buffer.data(kd + 118);
    const auto *kd_119 = buffer.data(kd + 119);
    const auto *kd_120 = buffer.data(kd + 120);
    const auto *kd_121 = buffer.data(kd + 121);
    const auto *kd_122 = buffer.data(kd + 122);
    const auto *kd_123 = buffer.data(kd + 123);
    const auto *kd_124 = buffer.data(kd + 124);
    const auto *kd_125 = buffer.data(kd + 125);
    const auto *kd_126 = buffer.data(kd + 126);
    const auto *kd_127 = buffer.data(kd + 127);
    const auto *kd_128 = buffer.data(kd + 128);
    const auto *kd_129 = buffer.data(kd + 129);
    const auto *kd_130 = buffer.data(kd + 130);
    const auto *kd_131 = buffer.data(kd + 131);
    const auto *kd_132 = buffer.data(kd + 132);
    const auto *kd_133 = buffer.data(kd + 133);
    const auto *kd_134 = buffer.data(kd + 134);
    const auto *kd_135 = buffer.data(kd + 135);
    const auto *kd_136 = buffer.data(kd + 136);
    const auto *kd_137 = buffer.data(kd + 137);
    const auto *kd_138 = buffer.data(kd + 138);
    const auto *kd_139 = buffer.data(kd + 139);
    const auto *kd_140 = buffer.data(kd + 140);
    const auto *kd_141 = buffer.data(kd + 141);
    const auto *kd_142 = buffer.data(kd + 142);
    const auto *kd_143 = buffer.data(kd + 143);
    const auto *kd_144 = buffer.data(kd + 144);
    const auto *kd_145 = buffer.data(kd + 145);
    const auto *kd_146 = buffer.data(kd + 146);
    const auto *kd_147 = buffer.data(kd + 147);
    const auto *kd_148 = buffer.data(kd + 148);
    const auto *kd_149 = buffer.data(kd + 149);
    const auto *kd_150 = buffer.data(kd + 150);
    const auto *kd_151 = buffer.data(kd + 151);
    const auto *kd_152 = buffer.data(kd + 152);
    const auto *kd_153 = buffer.data(kd + 153);
    const auto *kd_154 = buffer.data(kd + 154);
    const auto *kd_155 = buffer.data(kd + 155);
    const auto *kd_156 = buffer.data(kd + 156);
    const auto *kd_157 = buffer.data(kd + 157);
    const auto *kd_158 = buffer.data(kd + 158);
    const auto *kd_159 = buffer.data(kd + 159);
    const auto *kd_160 = buffer.data(kd + 160);
    const auto *kd_161 = buffer.data(kd + 161);
    const auto *kd_162 = buffer.data(kd + 162);
    const auto *kd_163 = buffer.data(kd + 163);
    const auto *kd_164 = buffer.data(kd + 164);
    const auto *kd_165 = buffer.data(kd + 165);
    const auto *kd_166 = buffer.data(kd + 166);
    const auto *kd_167 = buffer.data(kd + 167);
    const auto *kd_168 = buffer.data(kd + 168);
    const auto *kd_169 = buffer.data(kd + 169);
    const auto *kd_170 = buffer.data(kd + 170);
    const auto *kd_171 = buffer.data(kd + 171);
    const auto *kd_172 = buffer.data(kd + 172);
    const auto *kd_173 = buffer.data(kd + 173);
    const auto *kd_174 = buffer.data(kd + 174);
    const auto *kd_175 = buffer.data(kd + 175);
    const auto *kd_176 = buffer.data(kd + 176);
    const auto *kd_177 = buffer.data(kd + 177);
    const auto *kd_178 = buffer.data(kd + 178);
    const auto *kd_179 = buffer.data(kd + 179);
    const auto *kd_180 = buffer.data(kd + 180);
    const auto *kd_181 = buffer.data(kd + 181);
    const auto *kd_182 = buffer.data(kd + 182);
    const auto *kd_183 = buffer.data(kd + 183);
    const auto *kd_184 = buffer.data(kd + 184);
    const auto *kd_185 = buffer.data(kd + 185);
    const auto *kd_186 = buffer.data(kd + 186);
    const auto *kd_187 = buffer.data(kd + 187);
    const auto *kd_188 = buffer.data(kd + 188);
    const auto *kd_189 = buffer.data(kd + 189);
    const auto *kd_190 = buffer.data(kd + 190);
    const auto *kd_191 = buffer.data(kd + 191);
    const auto *kd_192 = buffer.data(kd + 192);
    const auto *kd_193 = buffer.data(kd + 193);
    const auto *kd_194 = buffer.data(kd + 194);
    const auto *kd_195 = buffer.data(kd + 195);
    const auto *kd_196 = buffer.data(kd + 196);
    const auto *kd_197 = buffer.data(kd + 197);
    const auto *kd_198 = buffer.data(kd + 198);
    const auto *kd_199 = buffer.data(kd + 199);
    const auto *kd_200 = buffer.data(kd + 200);
    const auto *kd_201 = buffer.data(kd + 201);
    const auto *kd_202 = buffer.data(kd + 202);
    const auto *kd_203 = buffer.data(kd + 203);
    const auto *kd_204 = buffer.data(kd + 204);
    const auto *kd_205 = buffer.data(kd + 205);
    const auto *kd_206 = buffer.data(kd + 206);
    const auto *kd_207 = buffer.data(kd + 207);
    const auto *kd_208 = buffer.data(kd + 208);
    const auto *kd_209 = buffer.data(kd + 209);
    const auto *kd_210 = buffer.data(kd + 210);
    const auto *kd_211 = buffer.data(kd + 211);
    const auto *kd_212 = buffer.data(kd + 212);
    const auto *kd_213 = buffer.data(kd + 213);
    const auto *kd_214 = buffer.data(kd + 214);
    const auto *kd_215 = buffer.data(kd + 215);

#pragma omp simd aligned(kd_7, kd_10, kd_37, kd_40, kd_91, kd_94, kd_169, \
                         kd_172 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * kd_7[k]
                 - f_1 * kd_37[k]
                 + f_2 * kd_91[k]
                 - f_3 * kd_169[k];

        g_1[k] = f_0 * kd_10[k]
                 - f_1 * kd_40[k]
                 + f_2 * kd_94[k]
                 - f_3 * kd_172[k];
    }

#pragma omp simd aligned(kd_6, kd_9, kd_11, kd_36, kd_39, kd_41, kd_90, kd_93, kd_95, kd_168, \
                         kd_171, kd_173 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_4 * kd_6[k]
                 - f_4 * kd_9[k]
                 + f_5 * kd_11[k]
                 + f_6 * kd_36[k]
                 + f_6 * kd_39[k]
                 - f_7 * kd_41[k]
                 - f_8 * kd_90[k]
                 - f_8 * kd_93[k]
                 + f_9 * kd_95[k]
                 + f_10 * kd_168[k]
                 + f_10 * kd_171[k]
                 - f_11 * kd_173[k];
    }

#pragma omp simd aligned(kd_6, kd_8, kd_9, kd_36, kd_38, kd_39, kd_90, kd_92, kd_93, kd_168, \
                         kd_170, kd_171 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_0 * kd_8[k]
                 - f_1 * kd_38[k]
                 + f_2 * kd_92[k]
                 - f_3 * kd_170[k];

        g_4[k] = f_12 * kd_6[k]
                 - f_12 * kd_9[k]
                 - f_13 * kd_36[k]
                 + f_13 * kd_39[k]
                 + f_14 * kd_90[k]
                 - f_14 * kd_93[k]
                 - f_15 * kd_168[k]
                 + f_15 * kd_171[k];
    }

#pragma omp simd aligned(kd_25, kd_28, kd_67, kd_70, kd_133, kd_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_16 * kd_25[k]
                 - f_17 * kd_67[k]
                 + f_16 * kd_133[k];

        g_6[k] = f_16 * kd_28[k]
                 - f_17 * kd_70[k]
                 + f_16 * kd_136[k];
    }

#pragma omp simd aligned(kd_24, kd_26, kd_27, kd_29, kd_66, kd_68, kd_69, kd_71, kd_132, \
                         kd_134, kd_135, kd_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_18 * kd_24[k]
                 - f_18 * kd_27[k]
                 + f_19 * kd_29[k]
                 + f_20 * kd_66[k]
                 + f_20 * kd_69[k]
                 - f_21 * kd_71[k]
                 - f_18 * kd_132[k]
                 - f_18 * kd_135[k]
                 + f_19 * kd_137[k];

        g_8[k] = f_16 * kd_26[k]
                 - f_17 * kd_68[k]
                 + f_16 * kd_134[k];

        g_9[k] = f_22 * kd_24[k]
                 - f_22 * kd_27[k]
                 - f_23 * kd_66[k]
                 + f_23 * kd_69[k]
                 + f_22 * kd_132[k]
                 - f_22 * kd_135[k];
    }

#pragma omp simd aligned(kd_7, kd_10, kd_37, kd_40, kd_49, kd_52, kd_91, kd_94, kd_103, \
                         kd_106, kd_169, kd_172, kd_181, kd_184 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_24 * kd_7[k]
                  + f_24 * kd_37[k]
                  + f_25 * kd_49[k]
                  + f_26 * kd_91[k]
                  - f_27 * kd_103[k]
                  - f_28 * kd_169[k]
                  + f_29 * kd_181[k];

        g_11[k] = -f_24 * kd_10[k]
                  + f_24 * kd_40[k]
                  + f_25 * kd_52[k]
                  + f_26 * kd_94[k]
                  - f_27 * kd_106[k]
                  - f_28 * kd_172[k]
                  + f_29 * kd_184[k];
    }

#pragma omp simd aligned(kd_6, kd_9, kd_11, kd_36, kd_39, kd_41, kd_48, kd_51, kd_53, kd_90, \
                         kd_93, kd_95, kd_102, kd_105, kd_107, kd_168, kd_171, kd_173, kd_180, \
                         kd_183, kd_185 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_30 * kd_6[k]
                  + f_30 * kd_9[k]
                  - f_31 * kd_11[k]
                  - f_30 * kd_36[k]
                  - f_30 * kd_39[k]
                  + f_31 * kd_41[k]
                  - f_32 * kd_48[k]
                  - f_32 * kd_51[k]
                  + f_33 * kd_53[k]
                  - f_34 * kd_90[k]
                  - f_34 * kd_93[k]
                  + f_35 * kd_95[k]
                  + f_33 * kd_102[k]
                  + f_33 * kd_105[k]
                  - f_36 * kd_107[k]
                  + f_37 * kd_168[k]
                  + f_37 * kd_171[k]
                  - f_38 * kd_173[k]
                  - f_39 * kd_180[k]
                  - f_39 * kd_183[k]
                  + f_40 * kd_185[k];
    }

#pragma omp simd aligned(kd_8, kd_38, kd_50, kd_92, kd_104, kd_170, \
                         kd_182 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_24 * kd_8[k]
                  + f_24 * kd_38[k]
                  + f_25 * kd_50[k]
                  + f_26 * kd_92[k]
                  - f_27 * kd_104[k]
                  - f_28 * kd_170[k]
                  + f_29 * kd_182[k];
    }

#pragma omp simd aligned(kd_6, kd_9, kd_36, kd_39, kd_48, kd_51, kd_90, kd_93, kd_102, kd_105, \
                         kd_168, kd_171, kd_180, kd_183 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_41 * kd_6[k]
                  + f_41 * kd_9[k]
                  + f_41 * kd_36[k]
                  - f_41 * kd_39[k]
                  + f_42 * kd_48[k]
                  - f_42 * kd_51[k]
                  + f_43 * kd_90[k]
                  - f_43 * kd_93[k]
                  - f_25 * kd_102[k]
                  + f_25 * kd_105[k]
                  - f_44 * kd_168[k]
                  + f_44 * kd_171[k]
                  + f_45 * kd_180[k]
                  - f_45 * kd_183[k];
    }

#pragma omp simd aligned(kd_25, kd_28, kd_79, kd_82, kd_133, kd_136, kd_145, \
                         kd_148 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_46 * kd_25[k]
                  + f_47 * kd_79[k]
                  + f_46 * kd_133[k]
                  - f_47 * kd_145[k];

        g_16[k] = -f_46 * kd_28[k]
                  + f_47 * kd_82[k]
                  + f_46 * kd_136[k]
                  - f_47 * kd_148[k];
    }

#pragma omp simd aligned(kd_24, kd_27, kd_29, kd_78, kd_81, kd_83, kd_132, kd_135, kd_137, \
                         kd_144, kd_147, kd_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_40 * kd_24[k]
                  + f_40 * kd_27[k]
                  - f_48 * kd_29[k]
                  - f_49 * kd_78[k]
                  - f_49 * kd_81[k]
                  + f_50 * kd_83[k]
                  - f_40 * kd_132[k]
                  - f_40 * kd_135[k]
                  + f_48 * kd_137[k]
                  + f_49 * kd_144[k]
                  + f_49 * kd_147[k]
                  - f_50 * kd_149[k];
    }

#pragma omp simd aligned(kd_24, kd_26, kd_27, kd_78, kd_80, kd_81, kd_132, kd_134, kd_135, \
                         kd_144, kd_146, kd_147 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_46 * kd_26[k]
                  + f_47 * kd_80[k]
                  + f_46 * kd_134[k]
                  - f_47 * kd_146[k];

        g_19[k] = -f_29 * kd_24[k]
                  + f_29 * kd_27[k]
                  + f_51 * kd_78[k]
                  - f_51 * kd_81[k]
                  + f_29 * kd_132[k]
                  - f_29 * kd_135[k]
                  - f_51 * kd_144[k]
                  + f_51 * kd_147[k];
    }

#pragma omp simd aligned(kd_7, kd_37, kd_49, kd_91, kd_103, kd_115, kd_169, kd_181, \
                         kd_193 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_52 * kd_7[k]
                  + f_53 * kd_37[k]
                  - f_54 * kd_49[k]
                  + f_55 * kd_91[k]
                  - f_56 * kd_103[k]
                  + f_57 * kd_115[k]
                  - f_55 * kd_169[k]
                  + f_58 * kd_181[k]
                  - f_59 * kd_193[k];
    }

#pragma omp simd aligned(kd_10, kd_40, kd_52, kd_94, kd_106, kd_118, kd_172, kd_184, \
                         kd_196 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_52 * kd_10[k]
                  + f_53 * kd_40[k]
                  - f_54 * kd_52[k]
                  + f_55 * kd_94[k]
                  - f_56 * kd_106[k]
                  + f_57 * kd_118[k]
                  - f_55 * kd_172[k]
                  + f_58 * kd_184[k]
                  - f_59 * kd_196[k];
    }

#pragma omp simd aligned(kd_6, kd_9, kd_11, kd_36, kd_39, kd_41, kd_48, kd_51, kd_53, kd_90, \
                         kd_93, kd_95, kd_102, kd_105, kd_107, kd_114, kd_117, kd_119, kd_168, \
                         kd_171, kd_173, kd_180, kd_183, kd_185, kd_192, kd_195, \
                         kd_197 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_60 * kd_6[k]
                  - f_60 * kd_9[k]
                  + f_61 * kd_11[k]
                  - f_62 * kd_36[k]
                  - f_62 * kd_39[k]
                  + f_63 * kd_41[k]
                  + f_64 * kd_48[k]
                  + f_64 * kd_51[k]
                  - f_65 * kd_53[k]
                  - f_66 * kd_90[k]
                  - f_66 * kd_93[k]
                  + f_67 * kd_95[k]
                  + f_68 * kd_102[k]
                  + f_68 * kd_105[k]
                  - f_69 * kd_107[k]
                  - f_69 * kd_114[k]
                  - f_69 * kd_117[k]
                  + f_70 * kd_119[k]
                  + f_66 * kd_168[k]
                  + f_66 * kd_171[k]
                  - f_67 * kd_173[k]
                  - f_71 * kd_180[k]
                  - f_71 * kd_183[k]
                  + f_68 * kd_185[k]
                  + f_72 * kd_192[k]
                  + f_72 * kd_195[k]
                  - f_73 * kd_197[k];
    }

#pragma omp simd aligned(kd_8, kd_38, kd_50, kd_92, kd_104, kd_116, kd_170, kd_182, \
                         kd_194 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_52 * kd_8[k]
                  + f_53 * kd_38[k]
                  - f_54 * kd_50[k]
                  + f_55 * kd_92[k]
                  - f_56 * kd_104[k]
                  + f_57 * kd_116[k]
                  - f_55 * kd_170[k]
                  + f_58 * kd_182[k]
                  - f_59 * kd_194[k];
    }

#pragma omp simd aligned(kd_6, kd_9, kd_36, kd_39, kd_48, kd_51, kd_90, kd_93, kd_102, kd_105, \
                         kd_114, kd_117, kd_168, kd_171, kd_180, kd_183, kd_192, \
                         kd_195 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_74 * kd_6[k]
                  - f_74 * kd_9[k]
                  + f_75 * kd_36[k]
                  - f_75 * kd_39[k]
                  - f_76 * kd_48[k]
                  + f_76 * kd_51[k]
                  + f_77 * kd_90[k]
                  - f_77 * kd_93[k]
                  - f_58 * kd_102[k]
                  + f_58 * kd_105[k]
                  + f_56 * kd_114[k]
                  - f_56 * kd_117[k]
                  - f_77 * kd_168[k]
                  + f_77 * kd_171[k]
                  + f_78 * kd_180[k]
                  - f_78 * kd_183[k]
                  - f_79 * kd_192[k]
                  + f_79 * kd_195[k];
    }

#pragma omp simd aligned(kd_25, kd_28, kd_67, kd_70, kd_79, kd_82, kd_133, kd_136, kd_145, \
                         kd_148, kd_157, kd_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_80 * kd_25[k]
                  + f_81 * kd_67[k]
                  - f_82 * kd_79[k]
                  + f_80 * kd_133[k]
                  - f_82 * kd_145[k]
                  + f_83 * kd_157[k];

        g_26[k] = f_80 * kd_28[k]
                  + f_81 * kd_70[k]
                  - f_82 * kd_82[k]
                  + f_80 * kd_136[k]
                  - f_82 * kd_148[k]
                  + f_83 * kd_160[k];
    }

#pragma omp simd aligned(kd_24, kd_27, kd_29, kd_66, kd_69, kd_71, kd_78, kd_81, kd_83, \
                         kd_132, kd_135, kd_137, kd_144, kd_147, kd_149, kd_156, kd_159, \
                         kd_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_84 * kd_24[k]
                  - f_84 * kd_27[k]
                  + f_85 * kd_29[k]
                  - f_85 * kd_66[k]
                  - f_85 * kd_69[k]
                  + f_86 * kd_71[k]
                  + f_87 * kd_78[k]
                  + f_87 * kd_81[k]
                  - f_88 * kd_83[k]
                  - f_84 * kd_132[k]
                  - f_84 * kd_135[k]
                  + f_85 * kd_137[k]
                  + f_87 * kd_144[k]
                  + f_87 * kd_147[k]
                  - f_88 * kd_149[k]
                  - f_89 * kd_156[k]
                  - f_89 * kd_159[k]
                  + f_90 * kd_161[k];
    }

#pragma omp simd aligned(kd_26, kd_68, kd_80, kd_134, kd_146, kd_158 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_80 * kd_26[k]
                  + f_81 * kd_68[k]
                  - f_82 * kd_80[k]
                  + f_80 * kd_134[k]
                  - f_82 * kd_146[k]
                  + f_83 * kd_158[k];
    }

#pragma omp simd aligned(kd_24, kd_27, kd_66, kd_69, kd_78, kd_81, kd_132, kd_135, kd_144, \
                         kd_147, kd_156, kd_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_91 * kd_24[k]
                  - f_91 * kd_27[k]
                  + f_80 * kd_66[k]
                  - f_80 * kd_69[k]
                  - f_92 * kd_78[k]
                  + f_92 * kd_81[k]
                  + f_91 * kd_132[k]
                  - f_91 * kd_135[k]
                  - f_92 * kd_144[k]
                  + f_92 * kd_147[k]
                  + f_93 * kd_156[k]
                  - f_93 * kd_159[k];
    }

#pragma omp simd aligned(kd_7, kd_37, kd_49, kd_91, kd_103, kd_115, kd_169, kd_181, kd_193, \
                         kd_205 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_94 * kd_7[k]
                  - f_63 * kd_37[k]
                  + f_69 * kd_49[k]
                  - f_63 * kd_91[k]
                  + f_70 * kd_103[k]
                  - f_70 * kd_115[k]
                  - f_94 * kd_169[k]
                  + f_69 * kd_181[k]
                  - f_70 * kd_193[k]
                  + f_95 * kd_205[k];
    }

#pragma omp simd aligned(kd_10, kd_40, kd_52, kd_94, kd_106, kd_118, kd_172, kd_184, kd_196, \
                         kd_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_94 * kd_10[k]
                  - f_63 * kd_40[k]
                  + f_69 * kd_52[k]
                  - f_63 * kd_94[k]
                  + f_70 * kd_106[k]
                  - f_70 * kd_118[k]
                  - f_94 * kd_172[k]
                  + f_69 * kd_184[k]
                  - f_70 * kd_196[k]
                  + f_95 * kd_208[k];
    }

#pragma omp simd aligned(kd_6, kd_9, kd_11, kd_36, kd_39, kd_41, kd_48, kd_51, kd_53, kd_90, \
                         kd_93, kd_95, kd_102, kd_105, kd_107, kd_114, kd_117, kd_119, kd_168, \
                         kd_171, kd_173, kd_180, kd_183, kd_185, kd_192, kd_195, kd_197, \
                         kd_204, kd_207, kd_209 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_96 * kd_6[k]
                  + f_96 * kd_9[k]
                  - f_97 * kd_11[k]
                  + f_98 * kd_36[k]
                  + f_98 * kd_39[k]
                  - f_99 * kd_41[k]
                  - f_100 * kd_48[k]
                  - f_100 * kd_51[k]
                  + f_79 * kd_53[k]
                  + f_98 * kd_90[k]
                  + f_98 * kd_93[k]
                  - f_99 * kd_95[k]
                  - f_79 * kd_102[k]
                  - f_79 * kd_105[k]
                  + f_59 * kd_107[k]
                  + f_79 * kd_114[k]
                  + f_79 * kd_117[k]
                  - f_59 * kd_119[k]
                  + f_96 * kd_168[k]
                  + f_96 * kd_171[k]
                  - f_97 * kd_173[k]
                  - f_100 * kd_180[k]
                  - f_100 * kd_183[k]
                  + f_79 * kd_185[k]
                  + f_79 * kd_192[k]
                  + f_79 * kd_195[k]
                  - f_59 * kd_197[k]
                  - f_101 * kd_204[k]
                  - f_101 * kd_207[k]
                  + f_102 * kd_209[k];
    }

#pragma omp simd aligned(kd_8, kd_38, kd_50, kd_92, kd_104, kd_116, kd_170, kd_182, kd_194, \
                         kd_206 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_94 * kd_8[k]
                  - f_63 * kd_38[k]
                  + f_69 * kd_50[k]
                  - f_63 * kd_92[k]
                  + f_70 * kd_104[k]
                  - f_70 * kd_116[k]
                  - f_94 * kd_170[k]
                  + f_69 * kd_182[k]
                  - f_70 * kd_194[k]
                  + f_95 * kd_206[k];
    }

#pragma omp simd aligned(kd_6, kd_9, kd_36, kd_39, kd_48, kd_51, kd_90, kd_93, kd_102, kd_105, \
                         kd_114, kd_117, kd_168, kd_171, kd_180, kd_183, kd_192, kd_195, \
                         kd_204, kd_207 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_103 * kd_6[k]
                  + f_103 * kd_9[k]
                  - f_62 * kd_36[k]
                  + f_62 * kd_39[k]
                  + f_68 * kd_48[k]
                  - f_68 * kd_51[k]
                  - f_62 * kd_90[k]
                  + f_62 * kd_93[k]
                  + f_69 * kd_102[k]
                  - f_69 * kd_105[k]
                  - f_69 * kd_114[k]
                  + f_69 * kd_117[k]
                  - f_103 * kd_168[k]
                  + f_103 * kd_171[k]
                  + f_68 * kd_180[k]
                  - f_68 * kd_183[k]
                  - f_69 * kd_192[k]
                  + f_69 * kd_195[k]
                  + f_104 * kd_204[k]
                  - f_104 * kd_207[k];
    }

#pragma omp simd aligned(kd_13, kd_43, kd_55, kd_97, kd_109, kd_121, kd_175, kd_187, kd_199, \
                         kd_211 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_105 * kd_13[k]
                  - f_106 * kd_43[k]
                  + f_107 * kd_55[k]
                  - f_106 * kd_97[k]
                  + f_108 * kd_109[k]
                  - f_109 * kd_121[k]
                  - f_105 * kd_175[k]
                  + f_107 * kd_187[k]
                  - f_109 * kd_199[k]
                  + f_110 * kd_211[k];
    }

#pragma omp simd aligned(kd_16, kd_46, kd_58, kd_100, kd_112, kd_124, kd_178, kd_190, kd_202, \
                         kd_214 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = -f_105 * kd_16[k]
                  - f_106 * kd_46[k]
                  + f_107 * kd_58[k]
                  - f_106 * kd_100[k]
                  + f_108 * kd_112[k]
                  - f_109 * kd_124[k]
                  - f_105 * kd_178[k]
                  + f_107 * kd_190[k]
                  - f_109 * kd_202[k]
                  + f_110 * kd_214[k];
    }

#pragma omp simd aligned(kd_12, kd_15, kd_17, kd_42, kd_45, kd_47, kd_54, kd_57, kd_59, kd_96, \
                         kd_99, kd_101, kd_108, kd_111, kd_113, kd_120, kd_123, kd_125, \
                         kd_174, kd_177, kd_179, kd_186, kd_189, kd_191, kd_198, kd_201, \
                         kd_203, kd_210, kd_213, kd_215 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = 1.09375 * kd_12[k]
                  + 1.09375 * kd_15[k]
                  - 2.1875 * kd_17[k]
                  + 3.28125 * kd_42[k]
                  + 3.28125 * kd_45[k]
                  - 6.5625 * kd_47[k]
                  - 6.5625 * kd_54[k]
                  - 6.5625 * kd_57[k]
                  + 13.125 * kd_59[k]
                  + 3.28125 * kd_96[k]
                  + 3.28125 * kd_99[k]
                  - 6.5625 * kd_101[k]
                  - 13.125 * kd_108[k]
                  - 13.125 * kd_111[k]
                  + 26.25 * kd_113[k]
                  + 5.25 * kd_120[k]
                  + 5.25 * kd_123[k]
                  - 10.5 * kd_125[k]
                  + 1.09375 * kd_174[k]
                  + 1.09375 * kd_177[k]
                  - 2.1875 * kd_179[k]
                  - 6.5625 * kd_186[k]
                  - 6.5625 * kd_189[k]
                  + 13.125 * kd_191[k]
                  + 5.25 * kd_198[k]
                  + 5.25 * kd_201[k]
                  - 10.5 * kd_203[k]
                  - 0.5 * kd_210[k]
                  - 0.5 * kd_213[k]
                  + kd_215[k];
    }

#pragma omp simd aligned(kd_14, kd_44, kd_56, kd_98, kd_110, kd_122, kd_176, kd_188, kd_200, \
                         kd_212 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_105 * kd_14[k]
                  - f_106 * kd_44[k]
                  + f_107 * kd_56[k]
                  - f_106 * kd_98[k]
                  + f_108 * kd_110[k]
                  - f_109 * kd_122[k]
                  - f_105 * kd_176[k]
                  + f_107 * kd_188[k]
                  - f_109 * kd_200[k]
                  + f_110 * kd_212[k];
    }

#pragma omp simd aligned(kd_12, kd_15, kd_42, kd_45, kd_54, kd_57, kd_96, kd_99, kd_108, \
                         kd_111, kd_120, kd_123, kd_174, kd_177, kd_186, kd_189, kd_198, \
                         kd_201, kd_210, kd_213 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_111 * kd_12[k]
                  + f_111 * kd_15[k]
                  - f_112 * kd_42[k]
                  + f_112 * kd_45[k]
                  + f_106 * kd_54[k]
                  - f_106 * kd_57[k]
                  - f_112 * kd_96[k]
                  + f_112 * kd_99[k]
                  + f_107 * kd_108[k]
                  - f_107 * kd_111[k]
                  - f_113 * kd_120[k]
                  + f_113 * kd_123[k]
                  - f_111 * kd_174[k]
                  + f_111 * kd_177[k]
                  + f_106 * kd_186[k]
                  - f_106 * kd_189[k]
                  - f_113 * kd_198[k]
                  + f_113 * kd_201[k]
                  + f_114 * kd_210[k]
                  - f_114 * kd_213[k];
    }

#pragma omp simd aligned(kd_1, kd_19, kd_31, kd_61, kd_73, kd_85, kd_127, kd_139, kd_151, \
                         kd_163 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_94 * kd_1[k]
                  - f_63 * kd_19[k]
                  + f_69 * kd_31[k]
                  - f_63 * kd_61[k]
                  + f_70 * kd_73[k]
                  - f_70 * kd_85[k]
                  - f_94 * kd_127[k]
                  + f_69 * kd_139[k]
                  - f_70 * kd_151[k]
                  + f_95 * kd_163[k];
    }

#pragma omp simd aligned(kd_4, kd_22, kd_34, kd_64, kd_76, kd_88, kd_130, kd_142, kd_154, \
                         kd_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_94 * kd_4[k]
                  - f_63 * kd_22[k]
                  + f_69 * kd_34[k]
                  - f_63 * kd_64[k]
                  + f_70 * kd_76[k]
                  - f_70 * kd_88[k]
                  - f_94 * kd_130[k]
                  + f_69 * kd_142[k]
                  - f_70 * kd_154[k]
                  + f_95 * kd_166[k];
    }

#pragma omp simd aligned(kd_0, kd_3, kd_5, kd_18, kd_21, kd_23, kd_30, kd_33, kd_35, kd_60, \
                         kd_63, kd_65, kd_72, kd_75, kd_77, kd_84, kd_87, kd_89, kd_126, \
                         kd_129, kd_131, kd_138, kd_141, kd_143, kd_150, kd_153, kd_155, \
                         kd_162, kd_165, kd_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_96 * kd_0[k]
                  + f_96 * kd_3[k]
                  - f_97 * kd_5[k]
                  + f_98 * kd_18[k]
                  + f_98 * kd_21[k]
                  - f_99 * kd_23[k]
                  - f_100 * kd_30[k]
                  - f_100 * kd_33[k]
                  + f_79 * kd_35[k]
                  + f_98 * kd_60[k]
                  + f_98 * kd_63[k]
                  - f_99 * kd_65[k]
                  - f_79 * kd_72[k]
                  - f_79 * kd_75[k]
                  + f_59 * kd_77[k]
                  + f_79 * kd_84[k]
                  + f_79 * kd_87[k]
                  - f_59 * kd_89[k]
                  + f_96 * kd_126[k]
                  + f_96 * kd_129[k]
                  - f_97 * kd_131[k]
                  - f_100 * kd_138[k]
                  - f_100 * kd_141[k]
                  + f_79 * kd_143[k]
                  + f_79 * kd_150[k]
                  + f_79 * kd_153[k]
                  - f_59 * kd_155[k]
                  - f_101 * kd_162[k]
                  - f_101 * kd_165[k]
                  + f_102 * kd_167[k];
    }

#pragma omp simd aligned(kd_2, kd_20, kd_32, kd_62, kd_74, kd_86, kd_128, kd_140, kd_152, \
                         kd_164 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_94 * kd_2[k]
                  - f_63 * kd_20[k]
                  + f_69 * kd_32[k]
                  - f_63 * kd_62[k]
                  + f_70 * kd_74[k]
                  - f_70 * kd_86[k]
                  - f_94 * kd_128[k]
                  + f_69 * kd_140[k]
                  - f_70 * kd_152[k]
                  + f_95 * kd_164[k];
    }

#pragma omp simd aligned(kd_0, kd_3, kd_18, kd_21, kd_30, kd_33, kd_60, kd_63, kd_72, kd_75, \
                         kd_84, kd_87, kd_126, kd_129, kd_138, kd_141, kd_150, kd_153, kd_162, \
                         kd_165 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_103 * kd_0[k]
                  + f_103 * kd_3[k]
                  - f_62 * kd_18[k]
                  + f_62 * kd_21[k]
                  + f_68 * kd_30[k]
                  - f_68 * kd_33[k]
                  - f_62 * kd_60[k]
                  + f_62 * kd_63[k]
                  + f_69 * kd_72[k]
                  - f_69 * kd_75[k]
                  - f_69 * kd_84[k]
                  + f_69 * kd_87[k]
                  - f_103 * kd_126[k]
                  + f_103 * kd_129[k]
                  + f_68 * kd_138[k]
                  - f_68 * kd_141[k]
                  - f_69 * kd_150[k]
                  + f_69 * kd_153[k]
                  + f_104 * kd_162[k]
                  - f_104 * kd_165[k];
    }

#pragma omp simd aligned(kd_13, kd_43, kd_55, kd_97, kd_121, kd_175, kd_187, \
                         kd_199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_91 * kd_13[k]
                  + f_91 * kd_43[k]
                  - f_92 * kd_55[k]
                  - f_91 * kd_97[k]
                  + f_93 * kd_121[k]
                  - f_91 * kd_175[k]
                  + f_92 * kd_187[k]
                  - f_93 * kd_199[k];
    }

#pragma omp simd aligned(kd_16, kd_46, kd_58, kd_100, kd_124, kd_178, kd_190, \
                         kd_202 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_91 * kd_16[k]
                  + f_91 * kd_46[k]
                  - f_92 * kd_58[k]
                  - f_91 * kd_100[k]
                  + f_93 * kd_124[k]
                  - f_91 * kd_178[k]
                  + f_92 * kd_190[k]
                  - f_93 * kd_202[k];
    }

#pragma omp simd aligned(kd_12, kd_15, kd_17, kd_42, kd_45, kd_47, kd_54, kd_57, kd_59, kd_96, \
                         kd_99, kd_101, kd_120, kd_123, kd_125, kd_174, kd_177, kd_179, \
                         kd_186, kd_189, kd_191, kd_198, kd_201, \
                         kd_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_115 * kd_12[k]
                  - f_115 * kd_15[k]
                  + f_84 * kd_17[k]
                  - f_115 * kd_42[k]
                  - f_115 * kd_45[k]
                  + f_84 * kd_47[k]
                  + f_116 * kd_54[k]
                  + f_116 * kd_57[k]
                  - f_87 * kd_59[k]
                  + f_115 * kd_96[k]
                  + f_115 * kd_99[k]
                  - f_84 * kd_101[k]
                  - f_117 * kd_120[k]
                  - f_117 * kd_123[k]
                  + f_89 * kd_125[k]
                  + f_115 * kd_174[k]
                  + f_115 * kd_177[k]
                  - f_84 * kd_179[k]
                  - f_116 * kd_186[k]
                  - f_116 * kd_189[k]
                  + f_87 * kd_191[k]
                  + f_117 * kd_198[k]
                  + f_117 * kd_201[k]
                  - f_89 * kd_203[k];
    }

#pragma omp simd aligned(kd_14, kd_44, kd_56, kd_98, kd_122, kd_176, kd_188, \
                         kd_200 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_91 * kd_14[k]
                  + f_91 * kd_44[k]
                  - f_92 * kd_56[k]
                  - f_91 * kd_98[k]
                  + f_93 * kd_122[k]
                  - f_91 * kd_176[k]
                  + f_92 * kd_188[k]
                  - f_93 * kd_200[k];
    }

#pragma omp simd aligned(kd_12, kd_15, kd_42, kd_45, kd_54, kd_57, kd_96, kd_99, kd_120, \
                         kd_123, kd_174, kd_177, kd_186, kd_189, kd_198, \
                         kd_201 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_118 * kd_12[k]
                  - f_118 * kd_15[k]
                  + f_118 * kd_42[k]
                  - f_118 * kd_45[k]
                  - f_119 * kd_54[k]
                  + f_119 * kd_57[k]
                  - f_118 * kd_96[k]
                  + f_118 * kd_99[k]
                  + f_120 * kd_120[k]
                  - f_120 * kd_123[k]
                  - f_118 * kd_174[k]
                  + f_118 * kd_177[k]
                  + f_119 * kd_186[k]
                  - f_119 * kd_189[k]
                  - f_120 * kd_198[k]
                  + f_120 * kd_201[k];
    }

#pragma omp simd aligned(kd_1, kd_19, kd_31, kd_61, kd_73, kd_85, kd_127, kd_139, \
                         kd_151 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = f_55 * kd_1[k]
                  - f_55 * kd_19[k]
                  - f_58 * kd_31[k]
                  - f_53 * kd_61[k]
                  + f_56 * kd_73[k]
                  + f_59 * kd_85[k]
                  - f_52 * kd_127[k]
                  + f_54 * kd_139[k]
                  - f_57 * kd_151[k];
    }

#pragma omp simd aligned(kd_4, kd_22, kd_34, kd_64, kd_76, kd_88, kd_130, kd_142, \
                         kd_154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_55 * kd_4[k]
                  - f_55 * kd_22[k]
                  - f_58 * kd_34[k]
                  - f_53 * kd_64[k]
                  + f_56 * kd_76[k]
                  + f_59 * kd_88[k]
                  - f_52 * kd_130[k]
                  + f_54 * kd_142[k]
                  - f_57 * kd_154[k];
    }

#pragma omp simd aligned(kd_0, kd_3, kd_5, kd_18, kd_21, kd_23, kd_30, kd_33, kd_35, kd_60, \
                         kd_63, kd_65, kd_72, kd_75, kd_77, kd_84, kd_87, kd_89, kd_126, \
                         kd_129, kd_131, kd_138, kd_141, kd_143, kd_150, kd_153, \
                         kd_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_66 * kd_0[k]
                  - f_66 * kd_3[k]
                  + f_67 * kd_5[k]
                  + f_66 * kd_18[k]
                  + f_66 * kd_21[k]
                  - f_67 * kd_23[k]
                  + f_71 * kd_30[k]
                  + f_71 * kd_33[k]
                  - f_68 * kd_35[k]
                  + f_62 * kd_60[k]
                  + f_62 * kd_63[k]
                  - f_63 * kd_65[k]
                  - f_68 * kd_72[k]
                  - f_68 * kd_75[k]
                  + f_69 * kd_77[k]
                  - f_72 * kd_84[k]
                  - f_72 * kd_87[k]
                  + f_73 * kd_89[k]
                  + f_60 * kd_126[k]
                  + f_60 * kd_129[k]
                  - f_61 * kd_131[k]
                  - f_64 * kd_138[k]
                  - f_64 * kd_141[k]
                  + f_65 * kd_143[k]
                  + f_69 * kd_150[k]
                  + f_69 * kd_153[k]
                  - f_70 * kd_155[k];
    }

#pragma omp simd aligned(kd_2, kd_20, kd_32, kd_62, kd_74, kd_86, kd_128, kd_140, \
                         kd_152 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_55 * kd_2[k]
                  - f_55 * kd_20[k]
                  - f_58 * kd_32[k]
                  - f_53 * kd_62[k]
                  + f_56 * kd_74[k]
                  + f_59 * kd_86[k]
                  - f_52 * kd_128[k]
                  + f_54 * kd_140[k]
                  - f_57 * kd_152[k];
    }

#pragma omp simd aligned(kd_0, kd_3, kd_18, kd_21, kd_30, kd_33, kd_60, kd_63, kd_72, kd_75, \
                         kd_84, kd_87, kd_126, kd_129, kd_138, kd_141, kd_150, \
                         kd_153 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_77 * kd_0[k]
                  - f_77 * kd_3[k]
                  - f_77 * kd_18[k]
                  + f_77 * kd_21[k]
                  - f_78 * kd_30[k]
                  + f_78 * kd_33[k]
                  - f_75 * kd_60[k]
                  + f_75 * kd_63[k]
                  + f_58 * kd_72[k]
                  - f_58 * kd_75[k]
                  + f_79 * kd_84[k]
                  - f_79 * kd_87[k]
                  - f_74 * kd_126[k]
                  + f_74 * kd_129[k]
                  + f_76 * kd_138[k]
                  - f_76 * kd_141[k]
                  - f_56 * kd_150[k]
                  + f_56 * kd_153[k];
    }

#pragma omp simd aligned(kd_13, kd_16, kd_43, kd_46, kd_55, kd_58, kd_97, kd_100, kd_109, \
                         kd_112, kd_175, kd_178, kd_187, kd_190 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_45 * kd_13[k]
                  + f_42 * kd_43[k]
                  + f_121 * kd_55[k]
                  + f_42 * kd_97[k]
                  - f_27 * kd_109[k]
                  - f_45 * kd_175[k]
                  + f_121 * kd_187[k];

        g_56[k] = -f_45 * kd_16[k]
                  + f_42 * kd_46[k]
                  + f_121 * kd_58[k]
                  + f_42 * kd_100[k]
                  - f_27 * kd_112[k]
                  - f_45 * kd_178[k]
                  + f_121 * kd_190[k];
    }

#pragma omp simd aligned(kd_12, kd_15, kd_17, kd_42, kd_45, kd_47, kd_54, kd_57, kd_59, kd_96, \
                         kd_99, kd_101, kd_108, kd_111, kd_113, kd_174, kd_177, kd_179, \
                         kd_186, kd_189, kd_191 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_122 * kd_12[k]
                  + f_122 * kd_15[k]
                  - f_39 * kd_17[k]
                  - f_123 * kd_42[k]
                  - f_123 * kd_45[k]
                  + f_32 * kd_47[k]
                  - f_124 * kd_54[k]
                  - f_124 * kd_57[k]
                  + f_125 * kd_59[k]
                  - f_123 * kd_96[k]
                  - f_123 * kd_99[k]
                  + f_32 * kd_101[k]
                  + f_33 * kd_108[k]
                  + f_33 * kd_111[k]
                  - f_36 * kd_113[k]
                  + f_122 * kd_174[k]
                  + f_122 * kd_177[k]
                  - f_39 * kd_179[k]
                  - f_124 * kd_186[k]
                  - f_124 * kd_189[k]
                  + f_125 * kd_191[k];
    }

#pragma omp simd aligned(kd_14, kd_44, kd_56, kd_98, kd_110, kd_176, \
                         kd_188 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_45 * kd_14[k]
                  + f_42 * kd_44[k]
                  + f_121 * kd_56[k]
                  + f_42 * kd_98[k]
                  - f_27 * kd_110[k]
                  - f_45 * kd_176[k]
                  + f_121 * kd_188[k];
    }

#pragma omp simd aligned(kd_12, kd_15, kd_42, kd_45, kd_54, kd_57, kd_96, kd_99, kd_108, \
                         kd_111, kd_174, kd_177, kd_186, kd_189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_126 * kd_12[k]
                  + f_126 * kd_15[k]
                  + f_127 * kd_42[k]
                  - f_127 * kd_45[k]
                  + f_128 * kd_54[k]
                  - f_128 * kd_57[k]
                  + f_127 * kd_96[k]
                  - f_127 * kd_99[k]
                  - f_25 * kd_108[k]
                  + f_25 * kd_111[k]
                  - f_126 * kd_174[k]
                  + f_126 * kd_177[k]
                  + f_128 * kd_186[k]
                  - f_128 * kd_189[k];
    }

#pragma omp simd aligned(kd_1, kd_4, kd_19, kd_22, kd_31, kd_34, kd_61, kd_64, kd_73, kd_76, \
                         kd_127, kd_130, kd_139, kd_142 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -f_28 * kd_1[k]
                  + f_26 * kd_19[k]
                  + f_29 * kd_31[k]
                  + f_24 * kd_61[k]
                  - f_27 * kd_73[k]
                  - f_24 * kd_127[k]
                  + f_25 * kd_139[k];

        g_61[k] = -f_28 * kd_4[k]
                  + f_26 * kd_22[k]
                  + f_29 * kd_34[k]
                  + f_24 * kd_64[k]
                  - f_27 * kd_76[k]
                  - f_24 * kd_130[k]
                  + f_25 * kd_142[k];
    }

#pragma omp simd aligned(kd_0, kd_3, kd_5, kd_18, kd_21, kd_23, kd_30, kd_33, kd_35, kd_60, \
                         kd_63, kd_65, kd_72, kd_75, kd_77, kd_126, kd_129, kd_131, kd_138, \
                         kd_141, kd_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_37 * kd_0[k]
                  + f_37 * kd_3[k]
                  - f_38 * kd_5[k]
                  - f_34 * kd_18[k]
                  - f_34 * kd_21[k]
                  + f_35 * kd_23[k]
                  - f_39 * kd_30[k]
                  - f_39 * kd_33[k]
                  + f_40 * kd_35[k]
                  - f_30 * kd_60[k]
                  - f_30 * kd_63[k]
                  + f_31 * kd_65[k]
                  + f_33 * kd_72[k]
                  + f_33 * kd_75[k]
                  - f_36 * kd_77[k]
                  + f_30 * kd_126[k]
                  + f_30 * kd_129[k]
                  - f_31 * kd_131[k]
                  - f_32 * kd_138[k]
                  - f_32 * kd_141[k]
                  + f_33 * kd_143[k];
    }

#pragma omp simd aligned(kd_2, kd_20, kd_32, kd_62, kd_74, kd_128, \
                         kd_140 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_28 * kd_2[k]
                  + f_26 * kd_20[k]
                  + f_29 * kd_32[k]
                  + f_24 * kd_62[k]
                  - f_27 * kd_74[k]
                  - f_24 * kd_128[k]
                  + f_25 * kd_140[k];
    }

#pragma omp simd aligned(kd_0, kd_3, kd_18, kd_21, kd_30, kd_33, kd_60, kd_63, kd_72, kd_75, \
                         kd_126, kd_129, kd_138, kd_141 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_44 * kd_0[k]
                  + f_44 * kd_3[k]
                  + f_43 * kd_18[k]
                  - f_43 * kd_21[k]
                  + f_45 * kd_30[k]
                  - f_45 * kd_33[k]
                  + f_41 * kd_60[k]
                  - f_41 * kd_63[k]
                  - f_25 * kd_72[k]
                  + f_25 * kd_75[k]
                  - f_41 * kd_126[k]
                  + f_41 * kd_129[k]
                  + f_42 * kd_138[k]
                  - f_42 * kd_141[k];
    }

#pragma omp simd aligned(kd_13, kd_16, kd_43, kd_46, kd_97, kd_100, kd_175, \
                         kd_178 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_129 * kd_13[k]
                  - f_130 * kd_43[k]
                  + f_130 * kd_97[k]
                  - f_129 * kd_175[k];

        g_66[k] = f_129 * kd_16[k]
                  - f_130 * kd_46[k]
                  + f_130 * kd_100[k]
                  - f_129 * kd_178[k];
    }

#pragma omp simd aligned(kd_12, kd_15, kd_17, kd_42, kd_45, kd_47, kd_96, kd_99, kd_101, \
                         kd_174, kd_177, kd_179 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_131 * kd_12[k]
                  - f_131 * kd_15[k]
                  + f_132 * kd_17[k]
                  + f_133 * kd_42[k]
                  + f_133 * kd_45[k]
                  - f_134 * kd_47[k]
                  - f_133 * kd_96[k]
                  - f_133 * kd_99[k]
                  + f_134 * kd_101[k]
                  + f_131 * kd_174[k]
                  + f_131 * kd_177[k]
                  - f_132 * kd_179[k];
    }

#pragma omp simd aligned(kd_12, kd_14, kd_15, kd_42, kd_44, kd_45, kd_96, kd_98, kd_99, \
                         kd_174, kd_176, kd_177 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_129 * kd_14[k]
                  - f_130 * kd_44[k]
                  + f_130 * kd_98[k]
                  - f_129 * kd_176[k];

        g_69[k] = f_135 * kd_12[k]
                  - f_135 * kd_15[k]
                  - f_136 * kd_42[k]
                  + f_136 * kd_45[k]
                  + f_136 * kd_96[k]
                  - f_136 * kd_99[k]
                  - f_135 * kd_174[k]
                  + f_135 * kd_177[k];
    }

#pragma omp simd aligned(kd_1, kd_4, kd_19, kd_22, kd_61, kd_64, kd_127, \
                         kd_130 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = f_3 * kd_1[k]
                  - f_2 * kd_19[k]
                  + f_1 * kd_61[k]
                  - f_0 * kd_127[k];

        g_71[k] = f_3 * kd_4[k]
                  - f_2 * kd_22[k]
                  + f_1 * kd_64[k]
                  - f_0 * kd_130[k];
    }

#pragma omp simd aligned(kd_0, kd_3, kd_5, kd_18, kd_21, kd_23, kd_60, kd_63, kd_65, kd_126, \
                         kd_129, kd_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -f_10 * kd_0[k]
                  - f_10 * kd_3[k]
                  + f_11 * kd_5[k]
                  + f_8 * kd_18[k]
                  + f_8 * kd_21[k]
                  - f_9 * kd_23[k]
                  - f_6 * kd_60[k]
                  - f_6 * kd_63[k]
                  + f_7 * kd_65[k]
                  + f_4 * kd_126[k]
                  + f_4 * kd_129[k]
                  - f_5 * kd_131[k];
    }

#pragma omp simd aligned(kd_0, kd_2, kd_3, kd_18, kd_20, kd_21, kd_60, kd_62, kd_63, kd_126, \
                         kd_128, kd_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_3 * kd_2[k]
                  - f_2 * kd_20[k]
                  + f_1 * kd_62[k]
                  - f_0 * kd_128[k];

        g_74[k] = f_15 * kd_0[k]
                  - f_15 * kd_3[k]
                  - f_14 * kd_18[k]
                  + f_14 * kd_21[k]
                  + f_13 * kd_60[k]
                  - f_13 * kd_63[k]
                  - f_12 * kd_126[k]
                  + f_12 * kd_129[k];
    }
}

}  // namespace simdtrf
