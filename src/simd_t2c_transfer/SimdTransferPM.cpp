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


#include "SimdTransferPM.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_hrr_pm(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t sm, const size_t sn, const size_t nmax) -> void
{
    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);
    auto *t_9 = buffer.data(target + 9);
    auto *t_10 = buffer.data(target + 10);
    auto *t_11 = buffer.data(target + 11);
    auto *t_12 = buffer.data(target + 12);
    auto *t_13 = buffer.data(target + 13);
    auto *t_14 = buffer.data(target + 14);
    auto *t_15 = buffer.data(target + 15);
    auto *t_16 = buffer.data(target + 16);
    auto *t_17 = buffer.data(target + 17);
    auto *t_18 = buffer.data(target + 18);
    auto *t_19 = buffer.data(target + 19);
    auto *t_20 = buffer.data(target + 20);
    auto *t_21 = buffer.data(target + 21);
    auto *t_22 = buffer.data(target + 22);
    auto *t_23 = buffer.data(target + 23);
    auto *t_24 = buffer.data(target + 24);
    auto *t_25 = buffer.data(target + 25);
    auto *t_26 = buffer.data(target + 26);
    auto *t_27 = buffer.data(target + 27);
    auto *t_28 = buffer.data(target + 28);
    auto *t_29 = buffer.data(target + 29);
    auto *t_30 = buffer.data(target + 30);
    auto *t_31 = buffer.data(target + 31);
    auto *t_32 = buffer.data(target + 32);
    auto *t_33 = buffer.data(target + 33);
    auto *t_34 = buffer.data(target + 34);
    auto *t_35 = buffer.data(target + 35);
    auto *t_36 = buffer.data(target + 36);
    auto *t_37 = buffer.data(target + 37);
    auto *t_38 = buffer.data(target + 38);
    auto *t_39 = buffer.data(target + 39);
    auto *t_40 = buffer.data(target + 40);
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);
    auto *t_45 = buffer.data(target + 45);
    auto *t_46 = buffer.data(target + 46);
    auto *t_47 = buffer.data(target + 47);
    auto *t_48 = buffer.data(target + 48);
    auto *t_49 = buffer.data(target + 49);
    auto *t_50 = buffer.data(target + 50);
    auto *t_51 = buffer.data(target + 51);
    auto *t_52 = buffer.data(target + 52);
    auto *t_53 = buffer.data(target + 53);
    auto *t_54 = buffer.data(target + 54);
    auto *t_55 = buffer.data(target + 55);
    auto *t_56 = buffer.data(target + 56);
    auto *t_57 = buffer.data(target + 57);
    auto *t_58 = buffer.data(target + 58);
    auto *t_59 = buffer.data(target + 59);
    auto *t_60 = buffer.data(target + 60);
    auto *t_61 = buffer.data(target + 61);
    auto *t_62 = buffer.data(target + 62);
    auto *t_63 = buffer.data(target + 63);
    auto *t_64 = buffer.data(target + 64);
    auto *t_65 = buffer.data(target + 65);
    auto *t_66 = buffer.data(target + 66);
    auto *t_67 = buffer.data(target + 67);
    auto *t_68 = buffer.data(target + 68);
    auto *t_69 = buffer.data(target + 69);
    auto *t_70 = buffer.data(target + 70);
    auto *t_71 = buffer.data(target + 71);
    auto *t_72 = buffer.data(target + 72);
    auto *t_73 = buffer.data(target + 73);
    auto *t_74 = buffer.data(target + 74);
    auto *t_75 = buffer.data(target + 75);
    auto *t_76 = buffer.data(target + 76);
    auto *t_77 = buffer.data(target + 77);
    auto *t_78 = buffer.data(target + 78);
    auto *t_79 = buffer.data(target + 79);
    auto *t_80 = buffer.data(target + 80);
    auto *t_81 = buffer.data(target + 81);
    auto *t_82 = buffer.data(target + 82);
    auto *t_83 = buffer.data(target + 83);
    auto *t_84 = buffer.data(target + 84);
    auto *t_85 = buffer.data(target + 85);
    auto *t_86 = buffer.data(target + 86);
    auto *t_87 = buffer.data(target + 87);
    auto *t_88 = buffer.data(target + 88);
    auto *t_89 = buffer.data(target + 89);
    auto *t_90 = buffer.data(target + 90);
    auto *t_91 = buffer.data(target + 91);
    auto *t_92 = buffer.data(target + 92);
    auto *t_93 = buffer.data(target + 93);
    auto *t_94 = buffer.data(target + 94);
    auto *t_95 = buffer.data(target + 95);
    auto *t_96 = buffer.data(target + 96);
    auto *t_97 = buffer.data(target + 97);
    auto *t_98 = buffer.data(target + 98);
    auto *t_99 = buffer.data(target + 99);
    auto *t_100 = buffer.data(target + 100);
    auto *t_101 = buffer.data(target + 101);
    auto *t_102 = buffer.data(target + 102);
    auto *t_103 = buffer.data(target + 103);
    auto *t_104 = buffer.data(target + 104);
    auto *t_105 = buffer.data(target + 105);
    auto *t_106 = buffer.data(target + 106);
    auto *t_107 = buffer.data(target + 107);
    auto *t_108 = buffer.data(target + 108);
    auto *t_109 = buffer.data(target + 109);
    auto *t_110 = buffer.data(target + 110);
    auto *t_111 = buffer.data(target + 111);
    auto *t_112 = buffer.data(target + 112);
    auto *t_113 = buffer.data(target + 113);
    auto *t_114 = buffer.data(target + 114);
    auto *t_115 = buffer.data(target + 115);
    auto *t_116 = buffer.data(target + 116);
    auto *t_117 = buffer.data(target + 117);
    auto *t_118 = buffer.data(target + 118);
    auto *t_119 = buffer.data(target + 119);
    auto *t_120 = buffer.data(target + 120);
    auto *t_121 = buffer.data(target + 121);
    auto *t_122 = buffer.data(target + 122);
    auto *t_123 = buffer.data(target + 123);
    auto *t_124 = buffer.data(target + 124);
    auto *t_125 = buffer.data(target + 125);
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);
    auto *t_135 = buffer.data(target + 135);
    auto *t_136 = buffer.data(target + 136);
    auto *t_137 = buffer.data(target + 137);
    auto *t_138 = buffer.data(target + 138);
    auto *t_139 = buffer.data(target + 139);
    auto *t_140 = buffer.data(target + 140);
    auto *t_141 = buffer.data(target + 141);
    auto *t_142 = buffer.data(target + 142);
    auto *t_143 = buffer.data(target + 143);
    auto *t_144 = buffer.data(target + 144);
    auto *t_145 = buffer.data(target + 145);
    auto *t_146 = buffer.data(target + 146);
    auto *t_147 = buffer.data(target + 147);
    auto *t_148 = buffer.data(target + 148);
    auto *t_149 = buffer.data(target + 149);
    auto *t_150 = buffer.data(target + 150);
    auto *t_151 = buffer.data(target + 151);
    auto *t_152 = buffer.data(target + 152);
    auto *t_153 = buffer.data(target + 153);
    auto *t_154 = buffer.data(target + 154);
    auto *t_155 = buffer.data(target + 155);
    auto *t_156 = buffer.data(target + 156);
    auto *t_157 = buffer.data(target + 157);
    auto *t_158 = buffer.data(target + 158);
    auto *t_159 = buffer.data(target + 159);
    auto *t_160 = buffer.data(target + 160);
    auto *t_161 = buffer.data(target + 161);
    auto *t_162 = buffer.data(target + 162);
    auto *t_163 = buffer.data(target + 163);
    auto *t_164 = buffer.data(target + 164);

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *sm_0 = buffer.data(sm + 0);
    const auto *sm_1 = buffer.data(sm + 1);
    const auto *sm_2 = buffer.data(sm + 2);
    const auto *sm_3 = buffer.data(sm + 3);
    const auto *sm_4 = buffer.data(sm + 4);
    const auto *sm_5 = buffer.data(sm + 5);
    const auto *sm_6 = buffer.data(sm + 6);
    const auto *sm_7 = buffer.data(sm + 7);
    const auto *sm_8 = buffer.data(sm + 8);
    const auto *sm_9 = buffer.data(sm + 9);
    const auto *sm_10 = buffer.data(sm + 10);
    const auto *sm_11 = buffer.data(sm + 11);
    const auto *sm_12 = buffer.data(sm + 12);
    const auto *sm_13 = buffer.data(sm + 13);
    const auto *sm_14 = buffer.data(sm + 14);
    const auto *sm_15 = buffer.data(sm + 15);
    const auto *sm_16 = buffer.data(sm + 16);
    const auto *sm_17 = buffer.data(sm + 17);
    const auto *sm_18 = buffer.data(sm + 18);
    const auto *sm_19 = buffer.data(sm + 19);
    const auto *sm_20 = buffer.data(sm + 20);
    const auto *sm_21 = buffer.data(sm + 21);
    const auto *sm_22 = buffer.data(sm + 22);
    const auto *sm_23 = buffer.data(sm + 23);
    const auto *sm_24 = buffer.data(sm + 24);
    const auto *sm_25 = buffer.data(sm + 25);
    const auto *sm_26 = buffer.data(sm + 26);
    const auto *sm_27 = buffer.data(sm + 27);
    const auto *sm_28 = buffer.data(sm + 28);
    const auto *sm_29 = buffer.data(sm + 29);
    const auto *sm_30 = buffer.data(sm + 30);
    const auto *sm_31 = buffer.data(sm + 31);
    const auto *sm_32 = buffer.data(sm + 32);
    const auto *sm_33 = buffer.data(sm + 33);
    const auto *sm_34 = buffer.data(sm + 34);
    const auto *sm_35 = buffer.data(sm + 35);
    const auto *sm_36 = buffer.data(sm + 36);
    const auto *sm_37 = buffer.data(sm + 37);
    const auto *sm_38 = buffer.data(sm + 38);
    const auto *sm_39 = buffer.data(sm + 39);
    const auto *sm_40 = buffer.data(sm + 40);
    const auto *sm_41 = buffer.data(sm + 41);
    const auto *sm_42 = buffer.data(sm + 42);
    const auto *sm_43 = buffer.data(sm + 43);
    const auto *sm_44 = buffer.data(sm + 44);
    const auto *sm_45 = buffer.data(sm + 45);
    const auto *sm_46 = buffer.data(sm + 46);
    const auto *sm_47 = buffer.data(sm + 47);
    const auto *sm_48 = buffer.data(sm + 48);
    const auto *sm_49 = buffer.data(sm + 49);
    const auto *sm_50 = buffer.data(sm + 50);
    const auto *sm_51 = buffer.data(sm + 51);
    const auto *sm_52 = buffer.data(sm + 52);
    const auto *sm_53 = buffer.data(sm + 53);
    const auto *sm_54 = buffer.data(sm + 54);

    const auto *sn_0 = buffer.data(sn + 0);
    const auto *sn_1 = buffer.data(sn + 1);
    const auto *sn_2 = buffer.data(sn + 2);
    const auto *sn_3 = buffer.data(sn + 3);
    const auto *sn_4 = buffer.data(sn + 4);
    const auto *sn_5 = buffer.data(sn + 5);
    const auto *sn_6 = buffer.data(sn + 6);
    const auto *sn_7 = buffer.data(sn + 7);
    const auto *sn_8 = buffer.data(sn + 8);
    const auto *sn_9 = buffer.data(sn + 9);
    const auto *sn_10 = buffer.data(sn + 10);
    const auto *sn_11 = buffer.data(sn + 11);
    const auto *sn_12 = buffer.data(sn + 12);
    const auto *sn_13 = buffer.data(sn + 13);
    const auto *sn_14 = buffer.data(sn + 14);
    const auto *sn_15 = buffer.data(sn + 15);
    const auto *sn_16 = buffer.data(sn + 16);
    const auto *sn_17 = buffer.data(sn + 17);
    const auto *sn_18 = buffer.data(sn + 18);
    const auto *sn_19 = buffer.data(sn + 19);
    const auto *sn_20 = buffer.data(sn + 20);
    const auto *sn_21 = buffer.data(sn + 21);
    const auto *sn_22 = buffer.data(sn + 22);
    const auto *sn_23 = buffer.data(sn + 23);
    const auto *sn_24 = buffer.data(sn + 24);
    const auto *sn_25 = buffer.data(sn + 25);
    const auto *sn_26 = buffer.data(sn + 26);
    const auto *sn_27 = buffer.data(sn + 27);
    const auto *sn_28 = buffer.data(sn + 28);
    const auto *sn_29 = buffer.data(sn + 29);
    const auto *sn_30 = buffer.data(sn + 30);
    const auto *sn_31 = buffer.data(sn + 31);
    const auto *sn_32 = buffer.data(sn + 32);
    const auto *sn_33 = buffer.data(sn + 33);
    const auto *sn_34 = buffer.data(sn + 34);
    const auto *sn_35 = buffer.data(sn + 35);
    const auto *sn_36 = buffer.data(sn + 36);
    const auto *sn_37 = buffer.data(sn + 37);
    const auto *sn_38 = buffer.data(sn + 38);
    const auto *sn_39 = buffer.data(sn + 39);
    const auto *sn_40 = buffer.data(sn + 40);
    const auto *sn_41 = buffer.data(sn + 41);
    const auto *sn_42 = buffer.data(sn + 42);
    const auto *sn_43 = buffer.data(sn + 43);
    const auto *sn_44 = buffer.data(sn + 44);
    const auto *sn_45 = buffer.data(sn + 45);
    const auto *sn_46 = buffer.data(sn + 46);
    const auto *sn_47 = buffer.data(sn + 47);
    const auto *sn_48 = buffer.data(sn + 48);
    const auto *sn_49 = buffer.data(sn + 49);
    const auto *sn_50 = buffer.data(sn + 50);
    const auto *sn_51 = buffer.data(sn + 51);
    const auto *sn_52 = buffer.data(sn + 52);
    const auto *sn_53 = buffer.data(sn + 53);
    const auto *sn_54 = buffer.data(sn + 54);
    const auto *sn_55 = buffer.data(sn + 55);
    const auto *sn_56 = buffer.data(sn + 56);
    const auto *sn_57 = buffer.data(sn + 57);
    const auto *sn_58 = buffer.data(sn + 58);
    const auto *sn_59 = buffer.data(sn + 59);
    const auto *sn_60 = buffer.data(sn + 60);
    const auto *sn_61 = buffer.data(sn + 61);
    const auto *sn_62 = buffer.data(sn + 62);
    const auto *sn_63 = buffer.data(sn + 63);
    const auto *sn_64 = buffer.data(sn + 64);
    const auto *sn_65 = buffer.data(sn + 65);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, sm_0, sm_1, sm_2, sm_3, sm_4, sn_0, \
                         sn_1, sn_2, sn_3, sn_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * sm_0[k]
                 + sn_0[k];

        t_1[k] = -ab_x[k] * sm_1[k]
                 + sn_1[k];

        t_2[k] = -ab_x[k] * sm_2[k]
                 + sn_2[k];

        t_3[k] = -ab_x[k] * sm_3[k]
                 + sn_3[k];

        t_4[k] = -ab_x[k] * sm_4[k]
                 + sn_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, sm_5, sm_6, sm_7, sm_8, sm_9, sn_5, \
                         sn_6, sn_7, sn_8, sn_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * sm_5[k]
                 + sn_5[k];

        t_6[k] = -ab_x[k] * sm_6[k]
                 + sn_6[k];

        t_7[k] = -ab_x[k] * sm_7[k]
                 + sn_7[k];

        t_8[k] = -ab_x[k] * sm_8[k]
                 + sn_8[k];

        t_9[k] = -ab_x[k] * sm_9[k]
                 + sn_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, sm_10, sm_11, sm_12, sm_13, \
                         sm_14, sn_10, sn_11, sn_12, sn_13, sn_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * sm_10[k]
                  + sn_10[k];

        t_11[k] = -ab_x[k] * sm_11[k]
                  + sn_11[k];

        t_12[k] = -ab_x[k] * sm_12[k]
                  + sn_12[k];

        t_13[k] = -ab_x[k] * sm_13[k]
                  + sn_13[k];

        t_14[k] = -ab_x[k] * sm_14[k]
                  + sn_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, sm_15, sm_16, sm_17, sm_18, \
                         sm_19, sn_15, sn_16, sn_17, sn_18, sn_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * sm_15[k]
                  + sn_15[k];

        t_16[k] = -ab_x[k] * sm_16[k]
                  + sn_16[k];

        t_17[k] = -ab_x[k] * sm_17[k]
                  + sn_17[k];

        t_18[k] = -ab_x[k] * sm_18[k]
                  + sn_18[k];

        t_19[k] = -ab_x[k] * sm_19[k]
                  + sn_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, sm_20, sm_21, sm_22, sm_23, \
                         sm_24, sn_20, sn_21, sn_22, sn_23, sn_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * sm_20[k]
                  + sn_20[k];

        t_21[k] = -ab_x[k] * sm_21[k]
                  + sn_21[k];

        t_22[k] = -ab_x[k] * sm_22[k]
                  + sn_22[k];

        t_23[k] = -ab_x[k] * sm_23[k]
                  + sn_23[k];

        t_24[k] = -ab_x[k] * sm_24[k]
                  + sn_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, sm_25, sm_26, sm_27, sm_28, \
                         sm_29, sn_25, sn_26, sn_27, sn_28, sn_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * sm_25[k]
                  + sn_25[k];

        t_26[k] = -ab_x[k] * sm_26[k]
                  + sn_26[k];

        t_27[k] = -ab_x[k] * sm_27[k]
                  + sn_27[k];

        t_28[k] = -ab_x[k] * sm_28[k]
                  + sn_28[k];

        t_29[k] = -ab_x[k] * sm_29[k]
                  + sn_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, sm_30, sm_31, sm_32, sm_33, \
                         sm_34, sn_30, sn_31, sn_32, sn_33, sn_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_x[k] * sm_30[k]
                  + sn_30[k];

        t_31[k] = -ab_x[k] * sm_31[k]
                  + sn_31[k];

        t_32[k] = -ab_x[k] * sm_32[k]
                  + sn_32[k];

        t_33[k] = -ab_x[k] * sm_33[k]
                  + sn_33[k];

        t_34[k] = -ab_x[k] * sm_34[k]
                  + sn_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, sm_35, sm_36, sm_37, sm_38, \
                         sm_39, sn_35, sn_36, sn_37, sn_38, sn_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_x[k] * sm_35[k]
                  + sn_35[k];

        t_36[k] = -ab_x[k] * sm_36[k]
                  + sn_36[k];

        t_37[k] = -ab_x[k] * sm_37[k]
                  + sn_37[k];

        t_38[k] = -ab_x[k] * sm_38[k]
                  + sn_38[k];

        t_39[k] = -ab_x[k] * sm_39[k]
                  + sn_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, sm_40, sm_41, sm_42, sm_43, \
                         sm_44, sn_40, sn_41, sn_42, sn_43, sn_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_x[k] * sm_40[k]
                  + sn_40[k];

        t_41[k] = -ab_x[k] * sm_41[k]
                  + sn_41[k];

        t_42[k] = -ab_x[k] * sm_42[k]
                  + sn_42[k];

        t_43[k] = -ab_x[k] * sm_43[k]
                  + sn_43[k];

        t_44[k] = -ab_x[k] * sm_44[k]
                  + sn_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, sm_45, sm_46, sm_47, sm_48, \
                         sm_49, sn_45, sn_46, sn_47, sn_48, sn_49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = -ab_x[k] * sm_45[k]
                  + sn_45[k];

        t_46[k] = -ab_x[k] * sm_46[k]
                  + sn_46[k];

        t_47[k] = -ab_x[k] * sm_47[k]
                  + sn_47[k];

        t_48[k] = -ab_x[k] * sm_48[k]
                  + sn_48[k];

        t_49[k] = -ab_x[k] * sm_49[k]
                  + sn_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, sm_50, sm_51, sm_52, sm_53, \
                         sm_54, sn_50, sn_51, sn_52, sn_53, sn_54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = -ab_x[k] * sm_50[k]
                  + sn_50[k];

        t_51[k] = -ab_x[k] * sm_51[k]
                  + sn_51[k];

        t_52[k] = -ab_x[k] * sm_52[k]
                  + sn_52[k];

        t_53[k] = -ab_x[k] * sm_53[k]
                  + sn_53[k];

        t_54[k] = -ab_x[k] * sm_54[k]
                  + sn_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_y, sm_0, sm_1, sm_2, sm_3, sm_4, \
                         sn_1, sn_3, sn_4, sn_6, sn_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = -ab_y[k] * sm_0[k]
                  + sn_1[k];

        t_56[k] = -ab_y[k] * sm_1[k]
                  + sn_3[k];

        t_57[k] = -ab_y[k] * sm_2[k]
                  + sn_4[k];

        t_58[k] = -ab_y[k] * sm_3[k]
                  + sn_6[k];

        t_59[k] = -ab_y[k] * sm_4[k]
                  + sn_7[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_y, sm_5, sm_6, sm_7, sm_8, sm_9, \
                         sn_8, sn_10, sn_11, sn_12, sn_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = -ab_y[k] * sm_5[k]
                  + sn_8[k];

        t_61[k] = -ab_y[k] * sm_6[k]
                  + sn_10[k];

        t_62[k] = -ab_y[k] * sm_7[k]
                  + sn_11[k];

        t_63[k] = -ab_y[k] * sm_8[k]
                  + sn_12[k];

        t_64[k] = -ab_y[k] * sm_9[k]
                  + sn_13[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_y, sm_10, sm_11, sm_12, sm_13, \
                         sm_14, sn_15, sn_16, sn_17, sn_18, sn_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = -ab_y[k] * sm_10[k]
                  + sn_15[k];

        t_66[k] = -ab_y[k] * sm_11[k]
                  + sn_16[k];

        t_67[k] = -ab_y[k] * sm_12[k]
                  + sn_17[k];

        t_68[k] = -ab_y[k] * sm_13[k]
                  + sn_18[k];

        t_69[k] = -ab_y[k] * sm_14[k]
                  + sn_19[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_y, sm_15, sm_16, sm_17, sm_18, \
                         sm_19, sn_21, sn_22, sn_23, sn_24, sn_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = -ab_y[k] * sm_15[k]
                  + sn_21[k];

        t_71[k] = -ab_y[k] * sm_16[k]
                  + sn_22[k];

        t_72[k] = -ab_y[k] * sm_17[k]
                  + sn_23[k];

        t_73[k] = -ab_y[k] * sm_18[k]
                  + sn_24[k];

        t_74[k] = -ab_y[k] * sm_19[k]
                  + sn_25[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_y, sm_20, sm_21, sm_22, sm_23, \
                         sm_24, sn_26, sn_28, sn_29, sn_30, sn_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = -ab_y[k] * sm_20[k]
                  + sn_26[k];

        t_76[k] = -ab_y[k] * sm_21[k]
                  + sn_28[k];

        t_77[k] = -ab_y[k] * sm_22[k]
                  + sn_29[k];

        t_78[k] = -ab_y[k] * sm_23[k]
                  + sn_30[k];

        t_79[k] = -ab_y[k] * sm_24[k]
                  + sn_31[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_y, sm_25, sm_26, sm_27, sm_28, \
                         sm_29, sn_32, sn_33, sn_34, sn_36, sn_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = -ab_y[k] * sm_25[k]
                  + sn_32[k];

        t_81[k] = -ab_y[k] * sm_26[k]
                  + sn_33[k];

        t_82[k] = -ab_y[k] * sm_27[k]
                  + sn_34[k];

        t_83[k] = -ab_y[k] * sm_28[k]
                  + sn_36[k];

        t_84[k] = -ab_y[k] * sm_29[k]
                  + sn_37[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_y, sm_30, sm_31, sm_32, sm_33, \
                         sm_34, sn_38, sn_39, sn_40, sn_41, sn_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = -ab_y[k] * sm_30[k]
                  + sn_38[k];

        t_86[k] = -ab_y[k] * sm_31[k]
                  + sn_39[k];

        t_87[k] = -ab_y[k] * sm_32[k]
                  + sn_40[k];

        t_88[k] = -ab_y[k] * sm_33[k]
                  + sn_41[k];

        t_89[k] = -ab_y[k] * sm_34[k]
                  + sn_42[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_y, sm_35, sm_36, sm_37, sm_38, \
                         sm_39, sn_43, sn_45, sn_46, sn_47, sn_48 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = -ab_y[k] * sm_35[k]
                  + sn_43[k];

        t_91[k] = -ab_y[k] * sm_36[k]
                  + sn_45[k];

        t_92[k] = -ab_y[k] * sm_37[k]
                  + sn_46[k];

        t_93[k] = -ab_y[k] * sm_38[k]
                  + sn_47[k];

        t_94[k] = -ab_y[k] * sm_39[k]
                  + sn_48[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_y, sm_40, sm_41, sm_42, sm_43, \
                         sm_44, sn_49, sn_50, sn_51, sn_52, sn_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = -ab_y[k] * sm_40[k]
                  + sn_49[k];

        t_96[k] = -ab_y[k] * sm_41[k]
                  + sn_50[k];

        t_97[k] = -ab_y[k] * sm_42[k]
                  + sn_51[k];

        t_98[k] = -ab_y[k] * sm_43[k]
                  + sn_52[k];

        t_99[k] = -ab_y[k] * sm_44[k]
                  + sn_53[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_y, sm_45, sm_46, sm_47, sm_48, \
                         sm_49, sn_55, sn_56, sn_57, sn_58, sn_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_100[k] = -ab_y[k] * sm_45[k]
                   + sn_55[k];

        t_101[k] = -ab_y[k] * sm_46[k]
                   + sn_56[k];

        t_102[k] = -ab_y[k] * sm_47[k]
                   + sn_57[k];

        t_103[k] = -ab_y[k] * sm_48[k]
                   + sn_58[k];

        t_104[k] = -ab_y[k] * sm_49[k]
                   + sn_59[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_y, sm_50, sm_51, sm_52, sm_53, \
                         sm_54, sn_60, sn_61, sn_62, sn_63, sn_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_105[k] = -ab_y[k] * sm_50[k]
                   + sn_60[k];

        t_106[k] = -ab_y[k] * sm_51[k]
                   + sn_61[k];

        t_107[k] = -ab_y[k] * sm_52[k]
                   + sn_62[k];

        t_108[k] = -ab_y[k] * sm_53[k]
                   + sn_63[k];

        t_109[k] = -ab_y[k] * sm_54[k]
                   + sn_64[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_z, sm_0, sm_1, sm_2, sm_3, \
                         sm_4, sn_2, sn_4, sn_5, sn_7, sn_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_110[k] = -ab_z[k] * sm_0[k]
                   + sn_2[k];

        t_111[k] = -ab_z[k] * sm_1[k]
                   + sn_4[k];

        t_112[k] = -ab_z[k] * sm_2[k]
                   + sn_5[k];

        t_113[k] = -ab_z[k] * sm_3[k]
                   + sn_7[k];

        t_114[k] = -ab_z[k] * sm_4[k]
                   + sn_8[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_z, sm_5, sm_6, sm_7, sm_8, \
                         sm_9, sn_9, sn_11, sn_12, sn_13, sn_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_115[k] = -ab_z[k] * sm_5[k]
                   + sn_9[k];

        t_116[k] = -ab_z[k] * sm_6[k]
                   + sn_11[k];

        t_117[k] = -ab_z[k] * sm_7[k]
                   + sn_12[k];

        t_118[k] = -ab_z[k] * sm_8[k]
                   + sn_13[k];

        t_119[k] = -ab_z[k] * sm_9[k]
                   + sn_14[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_z, sm_10, sm_11, sm_12, sm_13, \
                         sm_14, sn_16, sn_17, sn_18, sn_19, sn_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_120[k] = -ab_z[k] * sm_10[k]
                   + sn_16[k];

        t_121[k] = -ab_z[k] * sm_11[k]
                   + sn_17[k];

        t_122[k] = -ab_z[k] * sm_12[k]
                   + sn_18[k];

        t_123[k] = -ab_z[k] * sm_13[k]
                   + sn_19[k];

        t_124[k] = -ab_z[k] * sm_14[k]
                   + sn_20[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_z, sm_15, sm_16, sm_17, sm_18, \
                         sm_19, sn_22, sn_23, sn_24, sn_25, sn_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_125[k] = -ab_z[k] * sm_15[k]
                   + sn_22[k];

        t_126[k] = -ab_z[k] * sm_16[k]
                   + sn_23[k];

        t_127[k] = -ab_z[k] * sm_17[k]
                   + sn_24[k];

        t_128[k] = -ab_z[k] * sm_18[k]
                   + sn_25[k];

        t_129[k] = -ab_z[k] * sm_19[k]
                   + sn_26[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_z, sm_20, sm_21, sm_22, sm_23, \
                         sm_24, sn_27, sn_29, sn_30, sn_31, sn_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_130[k] = -ab_z[k] * sm_20[k]
                   + sn_27[k];

        t_131[k] = -ab_z[k] * sm_21[k]
                   + sn_29[k];

        t_132[k] = -ab_z[k] * sm_22[k]
                   + sn_30[k];

        t_133[k] = -ab_z[k] * sm_23[k]
                   + sn_31[k];

        t_134[k] = -ab_z[k] * sm_24[k]
                   + sn_32[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_z, sm_25, sm_26, sm_27, sm_28, \
                         sm_29, sn_33, sn_34, sn_35, sn_37, sn_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = -ab_z[k] * sm_25[k]
                   + sn_33[k];

        t_136[k] = -ab_z[k] * sm_26[k]
                   + sn_34[k];

        t_137[k] = -ab_z[k] * sm_27[k]
                   + sn_35[k];

        t_138[k] = -ab_z[k] * sm_28[k]
                   + sn_37[k];

        t_139[k] = -ab_z[k] * sm_29[k]
                   + sn_38[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_z, sm_30, sm_31, sm_32, sm_33, \
                         sm_34, sn_39, sn_40, sn_41, sn_42, sn_43 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_140[k] = -ab_z[k] * sm_30[k]
                   + sn_39[k];

        t_141[k] = -ab_z[k] * sm_31[k]
                   + sn_40[k];

        t_142[k] = -ab_z[k] * sm_32[k]
                   + sn_41[k];

        t_143[k] = -ab_z[k] * sm_33[k]
                   + sn_42[k];

        t_144[k] = -ab_z[k] * sm_34[k]
                   + sn_43[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_z, sm_35, sm_36, sm_37, sm_38, \
                         sm_39, sn_44, sn_46, sn_47, sn_48, sn_49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_145[k] = -ab_z[k] * sm_35[k]
                   + sn_44[k];

        t_146[k] = -ab_z[k] * sm_36[k]
                   + sn_46[k];

        t_147[k] = -ab_z[k] * sm_37[k]
                   + sn_47[k];

        t_148[k] = -ab_z[k] * sm_38[k]
                   + sn_48[k];

        t_149[k] = -ab_z[k] * sm_39[k]
                   + sn_49[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_z, sm_40, sm_41, sm_42, sm_43, \
                         sm_44, sn_50, sn_51, sn_52, sn_53, sn_54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_150[k] = -ab_z[k] * sm_40[k]
                   + sn_50[k];

        t_151[k] = -ab_z[k] * sm_41[k]
                   + sn_51[k];

        t_152[k] = -ab_z[k] * sm_42[k]
                   + sn_52[k];

        t_153[k] = -ab_z[k] * sm_43[k]
                   + sn_53[k];

        t_154[k] = -ab_z[k] * sm_44[k]
                   + sn_54[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_z, sm_45, sm_46, sm_47, sm_48, \
                         sm_49, sn_56, sn_57, sn_58, sn_59, sn_60 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_155[k] = -ab_z[k] * sm_45[k]
                   + sn_56[k];

        t_156[k] = -ab_z[k] * sm_46[k]
                   + sn_57[k];

        t_157[k] = -ab_z[k] * sm_47[k]
                   + sn_58[k];

        t_158[k] = -ab_z[k] * sm_48[k]
                   + sn_59[k];

        t_159[k] = -ab_z[k] * sm_49[k]
                   + sn_60[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_z, sm_50, sm_51, sm_52, sm_53, \
                         sm_54, sn_61, sn_62, sn_63, sn_64, sn_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_160[k] = -ab_z[k] * sm_50[k]
                   + sn_61[k];

        t_161[k] = -ab_z[k] * sm_51[k]
                   + sn_62[k];

        t_162[k] = -ab_z[k] * sm_52[k]
                   + sn_63[k];

        t_163[k] = -ab_z[k] * sm_53[k]
                   + sn_64[k];

        t_164[k] = -ab_z[k] * sm_54[k]
                   + sn_65[k];
    }
}

}  // namespace simdovl
