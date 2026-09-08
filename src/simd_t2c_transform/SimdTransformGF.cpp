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


#include "SimdTransformGF.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_gf(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t gf,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 1.875 * std::sqrt(14.0);
    const auto f_1 = 0.625 * std::sqrt(14.0);
    const auto f_2 = 2.5 * std::sqrt(21.0);
    const auto f_3 = 0.125 * std::sqrt(210.0);
    const auto f_4 = 0.5 * std::sqrt(210.0);
    const auto f_5 = 0.75 * std::sqrt(35.0);
    const auto f_6 = 0.5 * std::sqrt(35.0);
    const auto f_7 = 1.25 * std::sqrt(21.0);
    const auto f_8 = 5.625 * std::sqrt(7.0);
    const auto f_9 = 1.875 * std::sqrt(7.0);
    const auto f_10 = 0.625 * std::sqrt(7.0);
    const auto f_11 = 3.75 * std::sqrt(42.0);
    const auto f_12 = 1.25 * std::sqrt(42.0);
    const auto f_13 = 0.375 * std::sqrt(105.0);
    const auto f_14 = 1.5 * std::sqrt(105.0);
    const auto f_15 = 0.125 * std::sqrt(105.0);
    const auto f_16 = 0.5 * std::sqrt(105.0);
    const auto f_17 = 1.125 * std::sqrt(70.0);
    const auto f_18 = 0.75 * std::sqrt(70.0);
    const auto f_19 = 0.375 * std::sqrt(70.0);
    const auto f_20 = 0.25 * std::sqrt(70.0);
    const auto f_21 = 1.875 * std::sqrt(42.0);
    const auto f_22 = 0.625 * std::sqrt(42.0);
    const auto f_23 = 1.875 * std::sqrt(2.0);
    const auto f_24 = 0.625 * std::sqrt(2.0);
    const auto f_25 = 11.25 * std::sqrt(2.0);
    const auto f_26 = 3.75 * std::sqrt(2.0);
    const auto f_27 = 2.5 * std::sqrt(3.0);
    const auto f_28 = 15.0 * std::sqrt(3.0);
    const auto f_29 = 0.125 * std::sqrt(30.0);
    const auto f_30 = 0.5 * std::sqrt(30.0);
    const auto f_31 = 0.75 * std::sqrt(30.0);
    const auto f_32 = 3.0 * std::sqrt(30.0);
    const auto f_33 = 0.75 * std::sqrt(5.0);
    const auto f_34 = 0.5 * std::sqrt(5.0);
    const auto f_35 = 4.5 * std::sqrt(5.0);
    const auto f_36 = 3.0 * std::sqrt(5.0);
    const auto f_37 = 1.25 * std::sqrt(3.0);
    const auto f_38 = 7.5 * std::sqrt(3.0);
    const auto f_39 = 3.75 * std::sqrt(6.0);
    const auto f_40 = 5.0 * std::sqrt(6.0);
    const auto f_41 = 0.375 * std::sqrt(15.0);
    const auto f_42 = 1.5 * std::sqrt(15.0);
    const auto f_43 = 0.5 * std::sqrt(15.0);
    const auto f_44 = 2.0 * std::sqrt(15.0);
    const auto f_45 = 1.125 * std::sqrt(10.0);
    const auto f_46 = 0.75 * std::sqrt(10.0);
    const auto f_47 = 1.5 * std::sqrt(10.0);
    const auto f_48 = std::sqrt(10.0);
    const auto f_49 = 1.875 * std::sqrt(6.0);
    const auto f_50 = 2.5 * std::sqrt(6.0);
    const auto f_51 = 0.28125 * std::sqrt(10.0);
    const auto f_52 = 0.09375 * std::sqrt(10.0);
    const auto f_53 = 0.5625 * std::sqrt(10.0);
    const auto f_54 = 0.1875 * std::sqrt(10.0);
    const auto f_55 = 2.25 * std::sqrt(10.0);
    const auto f_56 = 0.25 * std::sqrt(10.0);
    const auto f_57 = 0.75 * std::sqrt(15.0);
    const auto f_58 = 3.0 * std::sqrt(15.0);
    const auto f_59 = std::sqrt(15.0);
    const auto f_60 = 0.09375 * std::sqrt(6.0);
    const auto f_61 = 0.375 * std::sqrt(6.0);
    const auto f_62 = 0.1875 * std::sqrt(6.0);
    const auto f_63 = 0.75 * std::sqrt(6.0);
    const auto f_64 = 3.0 * std::sqrt(6.0);
    const auto f_65 = 0.25 * std::sqrt(6.0);
    const auto f_66 = std::sqrt(6.0);
    const auto f_67 = 0.1875 * std::sqrt(15.0);
    const auto f_68 = 0.9375 * std::sqrt(2.0);
    const auto f_69 = 0.3125 * std::sqrt(2.0);
    const auto f_70 = 5.625 * std::sqrt(2.0);
    const auto f_71 = 0.0625 * std::sqrt(30.0);
    const auto f_72 = 0.25 * std::sqrt(30.0);
    const auto f_73 = 0.375 * std::sqrt(30.0);
    const auto f_74 = 1.5 * std::sqrt(30.0);
    const auto f_75 = 0.375 * std::sqrt(5.0);
    const auto f_76 = 0.25 * std::sqrt(5.0);
    const auto f_77 = 2.25 * std::sqrt(5.0);
    const auto f_78 = 1.5 * std::sqrt(5.0);
    const auto f_79 = 0.625 * std::sqrt(3.0);
    const auto f_80 = 3.75 * std::sqrt(3.0);
    const auto f_81 = 0.46875 * std::sqrt(14.0);
    const auto f_82 = 0.15625 * std::sqrt(14.0);
    const auto f_83 = 2.8125 * std::sqrt(14.0);
    const auto f_84 = 0.9375 * std::sqrt(14.0);
    const auto f_85 = 0.625 * std::sqrt(21.0);
    const auto f_86 = 3.75 * std::sqrt(21.0);
    const auto f_87 = 0.03125 * std::sqrt(210.0);
    const auto f_88 = 0.1875 * std::sqrt(210.0);
    const auto f_89 = 0.75 * std::sqrt(210.0);
    const auto f_90 = 0.1875 * std::sqrt(35.0);
    const auto f_91 = 0.125 * std::sqrt(35.0);
    const auto f_92 = 1.125 * std::sqrt(35.0);
    const auto f_93 = 0.3125 * std::sqrt(21.0);
    const auto f_94 = 1.875 * std::sqrt(21.0);

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

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_62 = buffer.data(gf + 62);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_64 = buffer.data(gf + 64);
    const auto *gf_65 = buffer.data(gf + 65);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_67 = buffer.data(gf + 67);
    const auto *gf_68 = buffer.data(gf + 68);
    const auto *gf_69 = buffer.data(gf + 69);
    const auto *gf_70 = buffer.data(gf + 70);
    const auto *gf_71 = buffer.data(gf + 71);
    const auto *gf_72 = buffer.data(gf + 72);
    const auto *gf_73 = buffer.data(gf + 73);
    const auto *gf_74 = buffer.data(gf + 74);
    const auto *gf_75 = buffer.data(gf + 75);
    const auto *gf_76 = buffer.data(gf + 76);
    const auto *gf_77 = buffer.data(gf + 77);
    const auto *gf_78 = buffer.data(gf + 78);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_80 = buffer.data(gf + 80);
    const auto *gf_81 = buffer.data(gf + 81);
    const auto *gf_82 = buffer.data(gf + 82);
    const auto *gf_83 = buffer.data(gf + 83);
    const auto *gf_84 = buffer.data(gf + 84);
    const auto *gf_85 = buffer.data(gf + 85);
    const auto *gf_86 = buffer.data(gf + 86);
    const auto *gf_87 = buffer.data(gf + 87);
    const auto *gf_88 = buffer.data(gf + 88);
    const auto *gf_89 = buffer.data(gf + 89);
    const auto *gf_90 = buffer.data(gf + 90);
    const auto *gf_91 = buffer.data(gf + 91);
    const auto *gf_92 = buffer.data(gf + 92);
    const auto *gf_93 = buffer.data(gf + 93);
    const auto *gf_94 = buffer.data(gf + 94);
    const auto *gf_95 = buffer.data(gf + 95);
    const auto *gf_96 = buffer.data(gf + 96);
    const auto *gf_97 = buffer.data(gf + 97);
    const auto *gf_98 = buffer.data(gf + 98);
    const auto *gf_99 = buffer.data(gf + 99);
    const auto *gf_100 = buffer.data(gf + 100);
    const auto *gf_101 = buffer.data(gf + 101);
    const auto *gf_102 = buffer.data(gf + 102);
    const auto *gf_103 = buffer.data(gf + 103);
    const auto *gf_104 = buffer.data(gf + 104);
    const auto *gf_105 = buffer.data(gf + 105);
    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_107 = buffer.data(gf + 107);
    const auto *gf_108 = buffer.data(gf + 108);
    const auto *gf_109 = buffer.data(gf + 109);
    const auto *gf_110 = buffer.data(gf + 110);
    const auto *gf_111 = buffer.data(gf + 111);
    const auto *gf_112 = buffer.data(gf + 112);
    const auto *gf_113 = buffer.data(gf + 113);
    const auto *gf_114 = buffer.data(gf + 114);
    const auto *gf_115 = buffer.data(gf + 115);
    const auto *gf_116 = buffer.data(gf + 116);
    const auto *gf_117 = buffer.data(gf + 117);
    const auto *gf_118 = buffer.data(gf + 118);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_120 = buffer.data(gf + 120);
    const auto *gf_121 = buffer.data(gf + 121);
    const auto *gf_122 = buffer.data(gf + 122);
    const auto *gf_123 = buffer.data(gf + 123);
    const auto *gf_124 = buffer.data(gf + 124);
    const auto *gf_125 = buffer.data(gf + 125);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_127 = buffer.data(gf + 127);
    const auto *gf_128 = buffer.data(gf + 128);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_130 = buffer.data(gf + 130);
    const auto *gf_131 = buffer.data(gf + 131);
    const auto *gf_132 = buffer.data(gf + 132);
    const auto *gf_133 = buffer.data(gf + 133);
    const auto *gf_134 = buffer.data(gf + 134);
    const auto *gf_135 = buffer.data(gf + 135);
    const auto *gf_136 = buffer.data(gf + 136);
    const auto *gf_137 = buffer.data(gf + 137);
    const auto *gf_138 = buffer.data(gf + 138);
    const auto *gf_139 = buffer.data(gf + 139);
    const auto *gf_140 = buffer.data(gf + 140);
    const auto *gf_141 = buffer.data(gf + 141);
    const auto *gf_142 = buffer.data(gf + 142);
    const auto *gf_143 = buffer.data(gf + 143);
    const auto *gf_144 = buffer.data(gf + 144);
    const auto *gf_145 = buffer.data(gf + 145);
    const auto *gf_146 = buffer.data(gf + 146);
    const auto *gf_147 = buffer.data(gf + 147);
    const auto *gf_148 = buffer.data(gf + 148);
    const auto *gf_149 = buffer.data(gf + 149);

#pragma omp simd aligned(gf_11, gf_14, gf_16, gf_18, gf_61, gf_64, gf_66, \
                         gf_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * gf_11[k]
                 - f_1 * gf_16[k]
                 - f_0 * gf_61[k]
                 + f_1 * gf_66[k];

        g_1[k] = f_2 * gf_14[k]
                 - f_2 * gf_64[k];

        g_2[k] = -f_3 * gf_11[k]
                 - f_3 * gf_16[k]
                 + f_4 * gf_18[k]
                 + f_3 * gf_61[k]
                 + f_3 * gf_66[k]
                 - f_4 * gf_68[k];
    }

#pragma omp simd aligned(gf_10, gf_12, gf_13, gf_15, gf_17, gf_19, gf_60, gf_62, gf_63, gf_65, \
                         gf_67, gf_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_5 * gf_12[k]
                 - f_5 * gf_17[k]
                 + f_6 * gf_19[k]
                 + f_5 * gf_62[k]
                 + f_5 * gf_67[k]
                 - f_6 * gf_69[k];

        g_4[k] = -f_3 * gf_10[k]
                 - f_3 * gf_13[k]
                 + f_4 * gf_15[k]
                 + f_3 * gf_60[k]
                 + f_3 * gf_63[k]
                 - f_4 * gf_65[k];

        g_5[k] = f_7 * gf_12[k]
                 - f_7 * gf_17[k]
                 - f_7 * gf_62[k]
                 + f_7 * gf_67[k];

        g_6[k] = f_1 * gf_10[k]
                 - f_0 * gf_13[k]
                 - f_1 * gf_60[k]
                 + f_0 * gf_63[k];
    }

#pragma omp simd aligned(gf_41, gf_44, gf_46, gf_48, gf_111, gf_114, gf_116, \
                         gf_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_8 * gf_41[k]
                 - f_9 * gf_46[k]
                 - f_9 * gf_111[k]
                 + f_10 * gf_116[k];

        g_8[k] = f_11 * gf_44[k]
                 - f_12 * gf_114[k];

        g_9[k] = -f_13 * gf_41[k]
                 - f_13 * gf_46[k]
                 + f_14 * gf_48[k]
                 + f_15 * gf_111[k]
                 + f_15 * gf_116[k]
                 - f_16 * gf_118[k];
    }

#pragma omp simd aligned(gf_40, gf_42, gf_43, gf_45, gf_47, gf_49, gf_110, gf_112, gf_113, \
                         gf_115, gf_117, gf_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_17 * gf_42[k]
                  - f_17 * gf_47[k]
                  + f_18 * gf_49[k]
                  + f_19 * gf_112[k]
                  + f_19 * gf_117[k]
                  - f_20 * gf_119[k];

        g_11[k] = -f_13 * gf_40[k]
                  - f_13 * gf_43[k]
                  + f_14 * gf_45[k]
                  + f_15 * gf_110[k]
                  + f_15 * gf_113[k]
                  - f_16 * gf_115[k];

        g_12[k] = f_21 * gf_42[k]
                  - f_21 * gf_47[k]
                  - f_22 * gf_112[k]
                  + f_22 * gf_117[k];

        g_13[k] = f_9 * gf_40[k]
                  - f_8 * gf_43[k]
                  - f_10 * gf_110[k]
                  + f_9 * gf_113[k];
    }

#pragma omp simd aligned(gf_11, gf_14, gf_16, gf_18, gf_61, gf_64, gf_66, gf_68, gf_81, gf_84, \
                         gf_86, gf_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_23 * gf_11[k]
                  + f_24 * gf_16[k]
                  - f_23 * gf_61[k]
                  + f_24 * gf_66[k]
                  + f_25 * gf_81[k]
                  - f_26 * gf_86[k];

        g_15[k] = -f_27 * gf_14[k]
                  - f_27 * gf_64[k]
                  + f_28 * gf_84[k];

        g_16[k] = f_29 * gf_11[k]
                  + f_29 * gf_16[k]
                  - f_30 * gf_18[k]
                  + f_29 * gf_61[k]
                  + f_29 * gf_66[k]
                  - f_30 * gf_68[k]
                  - f_31 * gf_81[k]
                  - f_31 * gf_86[k]
                  + f_32 * gf_88[k];
    }

#pragma omp simd aligned(gf_12, gf_17, gf_19, gf_62, gf_67, gf_69, gf_82, gf_87, \
                         gf_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_33 * gf_12[k]
                  + f_33 * gf_17[k]
                  - f_34 * gf_19[k]
                  + f_33 * gf_62[k]
                  + f_33 * gf_67[k]
                  - f_34 * gf_69[k]
                  - f_35 * gf_82[k]
                  - f_35 * gf_87[k]
                  + f_36 * gf_89[k];
    }

#pragma omp simd aligned(gf_10, gf_13, gf_15, gf_60, gf_63, gf_65, gf_80, gf_83, \
                         gf_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_29 * gf_10[k]
                  + f_29 * gf_13[k]
                  - f_30 * gf_15[k]
                  + f_29 * gf_60[k]
                  + f_29 * gf_63[k]
                  - f_30 * gf_65[k]
                  - f_31 * gf_80[k]
                  - f_31 * gf_83[k]
                  + f_32 * gf_85[k];
    }

#pragma omp simd aligned(gf_10, gf_12, gf_13, gf_17, gf_60, gf_62, gf_63, gf_67, gf_80, gf_82, \
                         gf_83, gf_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_37 * gf_12[k]
                  + f_37 * gf_17[k]
                  - f_37 * gf_62[k]
                  + f_37 * gf_67[k]
                  + f_38 * gf_82[k]
                  - f_38 * gf_87[k];

        g_20[k] = -f_24 * gf_10[k]
                  + f_23 * gf_13[k]
                  - f_24 * gf_60[k]
                  + f_23 * gf_63[k]
                  + f_26 * gf_80[k]
                  - f_25 * gf_83[k];
    }

#pragma omp simd aligned(gf_41, gf_44, gf_46, gf_48, gf_111, gf_114, gf_116, gf_118, gf_131, \
                         gf_134, gf_136, gf_138 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -5.625 * gf_41[k]
                  + 1.875 * gf_46[k]
                  - 5.625 * gf_111[k]
                  + 1.875 * gf_116[k]
                  + 7.5 * gf_131[k]
                  - 2.5 * gf_136[k];

        g_22[k] = -f_39 * gf_44[k]
                  - f_39 * gf_114[k]
                  + f_40 * gf_134[k];

        g_23[k] = f_41 * gf_41[k]
                  + f_41 * gf_46[k]
                  - f_42 * gf_48[k]
                  + f_41 * gf_111[k]
                  + f_41 * gf_116[k]
                  - f_42 * gf_118[k]
                  - f_43 * gf_131[k]
                  - f_43 * gf_136[k]
                  + f_44 * gf_138[k];
    }

#pragma omp simd aligned(gf_42, gf_47, gf_49, gf_112, gf_117, gf_119, gf_132, gf_137, \
                         gf_139 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_45 * gf_42[k]
                  + f_45 * gf_47[k]
                  - f_46 * gf_49[k]
                  + f_45 * gf_112[k]
                  + f_45 * gf_117[k]
                  - f_46 * gf_119[k]
                  - f_47 * gf_132[k]
                  - f_47 * gf_137[k]
                  + f_48 * gf_139[k];
    }

#pragma omp simd aligned(gf_40, gf_43, gf_45, gf_110, gf_113, gf_115, gf_130, gf_133, \
                         gf_135 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_41 * gf_40[k]
                  + f_41 * gf_43[k]
                  - f_42 * gf_45[k]
                  + f_41 * gf_110[k]
                  + f_41 * gf_113[k]
                  - f_42 * gf_115[k]
                  - f_43 * gf_130[k]
                  - f_43 * gf_133[k]
                  + f_44 * gf_135[k];
    }

#pragma omp simd aligned(gf_40, gf_42, gf_43, gf_47, gf_110, gf_112, gf_113, gf_117, gf_130, \
                         gf_132, gf_133, gf_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_49 * gf_42[k]
                  + f_49 * gf_47[k]
                  - f_49 * gf_112[k]
                  + f_49 * gf_117[k]
                  + f_50 * gf_132[k]
                  - f_50 * gf_137[k];

        g_27[k] = -1.875 * gf_40[k]
                  + 5.625 * gf_43[k]
                  - 1.875 * gf_110[k]
                  + 5.625 * gf_113[k]
                  + 2.5 * gf_130[k]
                  - 7.5 * gf_133[k];
    }

#pragma omp simd aligned(gf_1, gf_6, gf_31, gf_36, gf_51, gf_56, gf_101, gf_106, gf_121, \
                         gf_126, gf_141, gf_146 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_51 * gf_1[k]
                  - f_52 * gf_6[k]
                  + f_53 * gf_31[k]
                  - f_54 * gf_36[k]
                  - f_55 * gf_51[k]
                  + f_46 * gf_56[k]
                  + f_51 * gf_101[k]
                  - f_52 * gf_106[k]
                  - f_55 * gf_121[k]
                  + f_46 * gf_126[k]
                  + f_46 * gf_141[k]
                  - f_56 * gf_146[k];
    }

#pragma omp simd aligned(gf_4, gf_34, gf_54, gf_104, gf_124, gf_144 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_41 * gf_4[k]
                  + f_57 * gf_34[k]
                  - f_58 * gf_54[k]
                  + f_41 * gf_104[k]
                  - f_58 * gf_124[k]
                  + f_59 * gf_144[k];
    }

#pragma omp simd aligned(gf_1, gf_6, gf_8, gf_31, gf_36, gf_38, gf_51, gf_56, gf_58, gf_101, \
                         gf_106, gf_108, gf_121, gf_126, gf_128, gf_141, gf_146, \
                         gf_148 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_60 * gf_1[k]
                  - f_60 * gf_6[k]
                  + f_61 * gf_8[k]
                  - f_62 * gf_31[k]
                  - f_62 * gf_36[k]
                  + f_63 * gf_38[k]
                  + f_63 * gf_51[k]
                  + f_63 * gf_56[k]
                  - f_64 * gf_58[k]
                  - f_60 * gf_101[k]
                  - f_60 * gf_106[k]
                  + f_61 * gf_108[k]
                  + f_63 * gf_121[k]
                  + f_63 * gf_126[k]
                  - f_64 * gf_128[k]
                  - f_65 * gf_141[k]
                  - f_65 * gf_146[k]
                  + f_66 * gf_148[k];
    }

#pragma omp simd aligned(gf_2, gf_7, gf_9, gf_32, gf_37, gf_39, gf_52, gf_57, gf_59, gf_102, \
                         gf_107, gf_109, gf_122, gf_127, gf_129, gf_142, gf_147, \
                         gf_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -0.5625 * gf_2[k]
                  - 0.5625 * gf_7[k]
                  + 0.375 * gf_9[k]
                  - 1.125 * gf_32[k]
                  - 1.125 * gf_37[k]
                  + 0.75 * gf_39[k]
                  + 4.5 * gf_52[k]
                  + 4.5 * gf_57[k]
                  - 3.0 * gf_59[k]
                  - 0.5625 * gf_102[k]
                  - 0.5625 * gf_107[k]
                  + 0.375 * gf_109[k]
                  + 4.5 * gf_122[k]
                  + 4.5 * gf_127[k]
                  - 3.0 * gf_129[k]
                  - 1.5 * gf_142[k]
                  - 1.5 * gf_147[k]
                  + gf_149[k];
    }

#pragma omp simd aligned(gf_0, gf_3, gf_5, gf_30, gf_33, gf_35, gf_50, gf_53, gf_55, gf_100, \
                         gf_103, gf_105, gf_120, gf_123, gf_125, gf_140, gf_143, \
                         gf_145 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_60 * gf_0[k]
                  - f_60 * gf_3[k]
                  + f_61 * gf_5[k]
                  - f_62 * gf_30[k]
                  - f_62 * gf_33[k]
                  + f_63 * gf_35[k]
                  + f_63 * gf_50[k]
                  + f_63 * gf_53[k]
                  - f_64 * gf_55[k]
                  - f_60 * gf_100[k]
                  - f_60 * gf_103[k]
                  + f_61 * gf_105[k]
                  + f_63 * gf_120[k]
                  + f_63 * gf_123[k]
                  - f_64 * gf_125[k]
                  - f_65 * gf_140[k]
                  - f_65 * gf_143[k]
                  + f_66 * gf_145[k];
    }

#pragma omp simd aligned(gf_2, gf_7, gf_32, gf_37, gf_52, gf_57, gf_102, gf_107, gf_122, \
                         gf_127, gf_142, gf_147 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_67 * gf_2[k]
                  - f_67 * gf_7[k]
                  + f_41 * gf_32[k]
                  - f_41 * gf_37[k]
                  - f_42 * gf_52[k]
                  + f_42 * gf_57[k]
                  + f_67 * gf_102[k]
                  - f_67 * gf_107[k]
                  - f_42 * gf_122[k]
                  + f_42 * gf_127[k]
                  + f_43 * gf_142[k]
                  - f_43 * gf_147[k];
    }

#pragma omp simd aligned(gf_0, gf_3, gf_30, gf_33, gf_50, gf_53, gf_100, gf_103, gf_120, \
                         gf_123, gf_140, gf_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_52 * gf_0[k]
                  - f_51 * gf_3[k]
                  + f_54 * gf_30[k]
                  - f_53 * gf_33[k]
                  - f_46 * gf_50[k]
                  + f_55 * gf_53[k]
                  + f_52 * gf_100[k]
                  - f_51 * gf_103[k]
                  - f_46 * gf_120[k]
                  + f_55 * gf_123[k]
                  + f_56 * gf_140[k]
                  - f_46 * gf_143[k];
    }

#pragma omp simd aligned(gf_21, gf_24, gf_26, gf_28, gf_71, gf_74, gf_76, gf_78, gf_91, gf_94, \
                         gf_96, gf_98 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -5.625 * gf_21[k]
                  + 1.875 * gf_26[k]
                  - 5.625 * gf_71[k]
                  + 1.875 * gf_76[k]
                  + 7.5 * gf_91[k]
                  - 2.5 * gf_96[k];

        g_36[k] = -f_39 * gf_24[k]
                  - f_39 * gf_74[k]
                  + f_40 * gf_94[k];

        g_37[k] = f_41 * gf_21[k]
                  + f_41 * gf_26[k]
                  - f_42 * gf_28[k]
                  + f_41 * gf_71[k]
                  + f_41 * gf_76[k]
                  - f_42 * gf_78[k]
                  - f_43 * gf_91[k]
                  - f_43 * gf_96[k]
                  + f_44 * gf_98[k];
    }

#pragma omp simd aligned(gf_22, gf_27, gf_29, gf_72, gf_77, gf_79, gf_92, gf_97, \
                         gf_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_45 * gf_22[k]
                  + f_45 * gf_27[k]
                  - f_46 * gf_29[k]
                  + f_45 * gf_72[k]
                  + f_45 * gf_77[k]
                  - f_46 * gf_79[k]
                  - f_47 * gf_92[k]
                  - f_47 * gf_97[k]
                  + f_48 * gf_99[k];
    }

#pragma omp simd aligned(gf_20, gf_23, gf_25, gf_70, gf_73, gf_75, gf_90, gf_93, \
                         gf_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = f_41 * gf_20[k]
                  + f_41 * gf_23[k]
                  - f_42 * gf_25[k]
                  + f_41 * gf_70[k]
                  + f_41 * gf_73[k]
                  - f_42 * gf_75[k]
                  - f_43 * gf_90[k]
                  - f_43 * gf_93[k]
                  + f_44 * gf_95[k];
    }

#pragma omp simd aligned(gf_20, gf_22, gf_23, gf_27, gf_70, gf_72, gf_73, gf_77, gf_90, gf_92, \
                         gf_93, gf_97 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_49 * gf_22[k]
                  + f_49 * gf_27[k]
                  - f_49 * gf_72[k]
                  + f_49 * gf_77[k]
                  + f_50 * gf_92[k]
                  - f_50 * gf_97[k];

        g_41[k] = -1.875 * gf_20[k]
                  + 5.625 * gf_23[k]
                  - 1.875 * gf_70[k]
                  + 5.625 * gf_73[k]
                  + 2.5 * gf_90[k]
                  - 7.5 * gf_93[k];
    }

#pragma omp simd aligned(gf_1, gf_4, gf_6, gf_51, gf_54, gf_56, gf_101, gf_104, gf_106, \
                         gf_121, gf_124, gf_126 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_68 * gf_1[k]
                  + f_69 * gf_6[k]
                  + f_70 * gf_51[k]
                  - f_23 * gf_56[k]
                  + f_68 * gf_101[k]
                  - f_69 * gf_106[k]
                  - f_70 * gf_121[k]
                  + f_23 * gf_126[k];

        g_43[k] = -f_37 * gf_4[k]
                  + f_38 * gf_54[k]
                  + f_37 * gf_104[k]
                  - f_38 * gf_124[k];
    }

#pragma omp simd aligned(gf_1, gf_6, gf_8, gf_51, gf_56, gf_58, gf_101, gf_106, gf_108, \
                         gf_121, gf_126, gf_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_71 * gf_1[k]
                  + f_71 * gf_6[k]
                  - f_72 * gf_8[k]
                  - f_73 * gf_51[k]
                  - f_73 * gf_56[k]
                  + f_74 * gf_58[k]
                  - f_71 * gf_101[k]
                  - f_71 * gf_106[k]
                  + f_72 * gf_108[k]
                  + f_73 * gf_121[k]
                  + f_73 * gf_126[k]
                  - f_74 * gf_128[k];
    }

#pragma omp simd aligned(gf_2, gf_7, gf_9, gf_52, gf_57, gf_59, gf_102, gf_107, gf_109, \
                         gf_122, gf_127, gf_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_75 * gf_2[k]
                  + f_75 * gf_7[k]
                  - f_76 * gf_9[k]
                  - f_77 * gf_52[k]
                  - f_77 * gf_57[k]
                  + f_78 * gf_59[k]
                  - f_75 * gf_102[k]
                  - f_75 * gf_107[k]
                  + f_76 * gf_109[k]
                  + f_77 * gf_122[k]
                  + f_77 * gf_127[k]
                  - f_78 * gf_129[k];
    }

#pragma omp simd aligned(gf_0, gf_3, gf_5, gf_50, gf_53, gf_55, gf_100, gf_103, gf_105, \
                         gf_120, gf_123, gf_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_71 * gf_0[k]
                  + f_71 * gf_3[k]
                  - f_72 * gf_5[k]
                  - f_73 * gf_50[k]
                  - f_73 * gf_53[k]
                  + f_74 * gf_55[k]
                  - f_71 * gf_100[k]
                  - f_71 * gf_103[k]
                  + f_72 * gf_105[k]
                  + f_73 * gf_120[k]
                  + f_73 * gf_123[k]
                  - f_74 * gf_125[k];
    }

#pragma omp simd aligned(gf_2, gf_7, gf_52, gf_57, gf_102, gf_107, gf_122, \
                         gf_127 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_79 * gf_2[k]
                  + f_79 * gf_7[k]
                  + f_80 * gf_52[k]
                  - f_80 * gf_57[k]
                  + f_79 * gf_102[k]
                  - f_79 * gf_107[k]
                  - f_80 * gf_122[k]
                  + f_80 * gf_127[k];
    }

#pragma omp simd aligned(gf_0, gf_3, gf_21, gf_26, gf_50, gf_53, gf_71, gf_76, gf_100, gf_103, \
                         gf_120, gf_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_69 * gf_0[k]
                  + f_68 * gf_3[k]
                  + f_23 * gf_50[k]
                  - f_70 * gf_53[k]
                  + f_69 * gf_100[k]
                  - f_68 * gf_103[k]
                  - f_23 * gf_120[k]
                  + f_70 * gf_123[k];

        g_49[k] = f_9 * gf_21[k]
                  - f_10 * gf_26[k]
                  - f_8 * gf_71[k]
                  + f_9 * gf_76[k];
    }

#pragma omp simd aligned(gf_21, gf_24, gf_26, gf_28, gf_71, gf_74, gf_76, \
                         gf_78 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = f_12 * gf_24[k]
                  - f_11 * gf_74[k];

        g_51[k] = -f_15 * gf_21[k]
                  - f_15 * gf_26[k]
                  + f_16 * gf_28[k]
                  + f_13 * gf_71[k]
                  + f_13 * gf_76[k]
                  - f_14 * gf_78[k];
    }

#pragma omp simd aligned(gf_20, gf_22, gf_23, gf_25, gf_27, gf_29, gf_70, gf_72, gf_73, gf_75, \
                         gf_77, gf_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_19 * gf_22[k]
                  - f_19 * gf_27[k]
                  + f_20 * gf_29[k]
                  + f_17 * gf_72[k]
                  + f_17 * gf_77[k]
                  - f_18 * gf_79[k];

        g_53[k] = -f_15 * gf_20[k]
                  - f_15 * gf_23[k]
                  + f_16 * gf_25[k]
                  + f_13 * gf_70[k]
                  + f_13 * gf_73[k]
                  - f_14 * gf_75[k];

        g_54[k] = f_22 * gf_22[k]
                  - f_22 * gf_27[k]
                  - f_21 * gf_72[k]
                  + f_21 * gf_77[k];

        g_55[k] = f_10 * gf_20[k]
                  - f_9 * gf_23[k]
                  - f_9 * gf_70[k]
                  + f_8 * gf_73[k];
    }

#pragma omp simd aligned(gf_1, gf_4, gf_6, gf_8, gf_31, gf_34, gf_36, gf_38, gf_101, gf_104, \
                         gf_106, gf_108 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_81 * gf_1[k]
                  - f_82 * gf_6[k]
                  - f_83 * gf_31[k]
                  + f_84 * gf_36[k]
                  + f_81 * gf_101[k]
                  - f_82 * gf_106[k];

        g_57[k] = f_85 * gf_4[k]
                  - f_86 * gf_34[k]
                  + f_85 * gf_104[k];

        g_58[k] = -f_87 * gf_1[k]
                  - f_87 * gf_6[k]
                  + f_3 * gf_8[k]
                  + f_88 * gf_31[k]
                  + f_88 * gf_36[k]
                  - f_89 * gf_38[k]
                  - f_87 * gf_101[k]
                  - f_87 * gf_106[k]
                  + f_3 * gf_108[k];
    }

#pragma omp simd aligned(gf_2, gf_7, gf_9, gf_32, gf_37, gf_39, gf_102, gf_107, \
                         gf_109 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_90 * gf_2[k]
                  - f_90 * gf_7[k]
                  + f_91 * gf_9[k]
                  + f_92 * gf_32[k]
                  + f_92 * gf_37[k]
                  - f_5 * gf_39[k]
                  - f_90 * gf_102[k]
                  - f_90 * gf_107[k]
                  + f_91 * gf_109[k];
    }

#pragma omp simd aligned(gf_0, gf_3, gf_5, gf_30, gf_33, gf_35, gf_100, gf_103, \
                         gf_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -f_87 * gf_0[k]
                  - f_87 * gf_3[k]
                  + f_3 * gf_5[k]
                  + f_88 * gf_30[k]
                  + f_88 * gf_33[k]
                  - f_89 * gf_35[k]
                  - f_87 * gf_100[k]
                  - f_87 * gf_103[k]
                  + f_3 * gf_105[k];
    }

#pragma omp simd aligned(gf_0, gf_2, gf_3, gf_7, gf_30, gf_32, gf_33, gf_37, gf_100, gf_102, \
                         gf_103, gf_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_93 * gf_2[k]
                  - f_93 * gf_7[k]
                  - f_94 * gf_32[k]
                  + f_94 * gf_37[k]
                  + f_93 * gf_102[k]
                  - f_93 * gf_107[k];

        g_62[k] = f_82 * gf_0[k]
                  - f_81 * gf_3[k]
                  - f_84 * gf_30[k]
                  + f_83 * gf_33[k]
                  + f_82 * gf_100[k]
                  - f_81 * gf_103[k];
    }
}

}  // namespace simdtrf
