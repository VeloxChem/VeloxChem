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


#include "SimdThreeCenterElectronRepulsionVrrRecSOP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sop_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t sns, const size_t sos,
                                                          const size_t ncols, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = p / q;
    const auto f_2 = 5.0 / q;
    const auto f_3 = 0.5 / q;
    const auto f_4 = 4.5 / q;
    const auto f_5 = 1.0 / q;
    const auto f_6 = 4.0 / q;
    const auto f_7 = 1.5 / q;
    const auto f_8 = 3.5 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 3.0 / q;
    const auto f_11 = 2.5 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sns_0 = buffer.data(sns + 0);
    const auto *sns_1 = buffer.data(sns + 1);
    const auto *sns_2 = buffer.data(sns + 2);
    const auto *sns_3 = buffer.data(sns + 3);
    const auto *sns_4 = buffer.data(sns + 4);
    const auto *sns_5 = buffer.data(sns + 5);
    const auto *sns_6 = buffer.data(sns + 6);
    const auto *sns_7 = buffer.data(sns + 7);
    const auto *sns_8 = buffer.data(sns + 8);
    const auto *sns_9 = buffer.data(sns + 9);
    const auto *sns_10 = buffer.data(sns + 10);
    const auto *sns_11 = buffer.data(sns + 11);
    const auto *sns_12 = buffer.data(sns + 12);
    const auto *sns_13 = buffer.data(sns + 13);
    const auto *sns_14 = buffer.data(sns + 14);
    const auto *sns_15 = buffer.data(sns + 15);
    const auto *sns_16 = buffer.data(sns + 16);
    const auto *sns_17 = buffer.data(sns + 17);
    const auto *sns_18 = buffer.data(sns + 18);
    const auto *sns_19 = buffer.data(sns + 19);
    const auto *sns_20 = buffer.data(sns + 20);
    const auto *sns_21 = buffer.data(sns + 21);
    const auto *sns_22 = buffer.data(sns + 22);
    const auto *sns_23 = buffer.data(sns + 23);
    const auto *sns_24 = buffer.data(sns + 24);
    const auto *sns_25 = buffer.data(sns + 25);
    const auto *sns_26 = buffer.data(sns + 26);
    const auto *sns_27 = buffer.data(sns + 27);
    const auto *sns_28 = buffer.data(sns + 28);
    const auto *sns_29 = buffer.data(sns + 29);
    const auto *sns_30 = buffer.data(sns + 30);
    const auto *sns_31 = buffer.data(sns + 31);
    const auto *sns_32 = buffer.data(sns + 32);
    const auto *sns_33 = buffer.data(sns + 33);
    const auto *sns_34 = buffer.data(sns + 34);
    const auto *sns_35 = buffer.data(sns + 35);
    const auto *sns_36 = buffer.data(sns + 36);
    const auto *sns_37 = buffer.data(sns + 37);
    const auto *sns_38 = buffer.data(sns + 38);
    const auto *sns_39 = buffer.data(sns + 39);
    const auto *sns_40 = buffer.data(sns + 40);
    const auto *sns_41 = buffer.data(sns + 41);
    const auto *sns_42 = buffer.data(sns + 42);
    const auto *sns_43 = buffer.data(sns + 43);
    const auto *sns_44 = buffer.data(sns + 44);
    const auto *sns_45 = buffer.data(sns + 45);
    const auto *sns_46 = buffer.data(sns + 46);
    const auto *sns_47 = buffer.data(sns + 47);
    const auto *sns_48 = buffer.data(sns + 48);
    const auto *sns_49 = buffer.data(sns + 49);
    const auto *sns_50 = buffer.data(sns + 50);
    const auto *sns_51 = buffer.data(sns + 51);
    const auto *sns_52 = buffer.data(sns + 52);
    const auto *sns_53 = buffer.data(sns + 53);

    const auto *sos_0 = buffer.data(sos + 0);
    const auto *sos_1 = buffer.data(sos + 1);
    const auto *sos_2 = buffer.data(sos + 2);
    const auto *sos_3 = buffer.data(sos + 3);
    const auto *sos_4 = buffer.data(sos + 4);
    const auto *sos_5 = buffer.data(sos + 5);
    const auto *sos_6 = buffer.data(sos + 6);
    const auto *sos_7 = buffer.data(sos + 7);
    const auto *sos_8 = buffer.data(sos + 8);
    const auto *sos_9 = buffer.data(sos + 9);
    const auto *sos_10 = buffer.data(sos + 10);
    const auto *sos_11 = buffer.data(sos + 11);
    const auto *sos_12 = buffer.data(sos + 12);
    const auto *sos_13 = buffer.data(sos + 13);
    const auto *sos_14 = buffer.data(sos + 14);
    const auto *sos_15 = buffer.data(sos + 15);
    const auto *sos_16 = buffer.data(sos + 16);
    const auto *sos_17 = buffer.data(sos + 17);
    const auto *sos_18 = buffer.data(sos + 18);
    const auto *sos_19 = buffer.data(sos + 19);
    const auto *sos_20 = buffer.data(sos + 20);
    const auto *sos_21 = buffer.data(sos + 21);
    const auto *sos_22 = buffer.data(sos + 22);
    const auto *sos_23 = buffer.data(sos + 23);
    const auto *sos_24 = buffer.data(sos + 24);
    const auto *sos_25 = buffer.data(sos + 25);
    const auto *sos_26 = buffer.data(sos + 26);
    const auto *sos_27 = buffer.data(sos + 27);
    const auto *sos_28 = buffer.data(sos + 28);
    const auto *sos_29 = buffer.data(sos + 29);
    const auto *sos_30 = buffer.data(sos + 30);
    const auto *sos_31 = buffer.data(sos + 31);
    const auto *sos_32 = buffer.data(sos + 32);
    const auto *sos_33 = buffer.data(sos + 33);
    const auto *sos_34 = buffer.data(sos + 34);
    const auto *sos_35 = buffer.data(sos + 35);
    const auto *sos_36 = buffer.data(sos + 36);
    const auto *sos_37 = buffer.data(sos + 37);
    const auto *sos_38 = buffer.data(sos + 38);
    const auto *sos_39 = buffer.data(sos + 39);
    const auto *sos_40 = buffer.data(sos + 40);
    const auto *sos_41 = buffer.data(sos + 41);
    const auto *sos_42 = buffer.data(sos + 42);
    const auto *sos_43 = buffer.data(sos + 43);
    const auto *sos_44 = buffer.data(sos + 44);
    const auto *sos_45 = buffer.data(sos + 45);
    const auto *sos_46 = buffer.data(sos + 46);
    const auto *sos_47 = buffer.data(sos + 47);
    const auto *sos_48 = buffer.data(sos + 48);
    const auto *sos_49 = buffer.data(sos + 49);
    const auto *sos_50 = buffer.data(sos + 50);
    const auto *sos_51 = buffer.data(sos + 51);
    const auto *sos_52 = buffer.data(sos + 52);
    const auto *sos_53 = buffer.data(sos + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pc_x, pc_y, pc_z, sns_0, sns_1, \
                         sns_2, sos_0, sos_1, sos_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sns_0[k]
                 + f_1 * pc_x[k] * sos_0[k];

        t_1[k] = f_1 * pc_y[k] * sos_0[k];

        t_2[k] = f_1 * pc_z[k] * sos_0[k];

        t_3[k] = f_2 * sns_1[k]
                 + f_1 * pc_x[k] * sos_1[k];

        t_4[k] = f_3 * sns_0[k]
                 + f_1 * pc_y[k] * sos_1[k];

        t_5[k] = f_1 * pc_z[k] * sos_1[k];

        t_6[k] = f_2 * sns_2[k]
                 + f_1 * pc_x[k] * sos_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pc_x, pc_y, pc_z, sns_0, sns_1, \
                         sns_3, sns_4, sos_2, sos_3, sos_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * pc_y[k] * sos_2[k];

        t_8[k] = f_3 * sns_0[k]
                 + f_1 * pc_z[k] * sos_2[k];

        t_9[k] = f_4 * sns_3[k]
                 + f_1 * pc_x[k] * sos_3[k];

        t_10[k] = f_5 * sns_1[k]
                  + f_1 * pc_y[k] * sos_3[k];

        t_11[k] = f_1 * pc_z[k] * sos_3[k];

        t_12[k] = f_4 * sns_4[k]
                  + f_1 * pc_x[k] * sos_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, sns_1, sns_2, \
                         sns_5, sns_6, sos_4, sos_5, sos_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * sns_2[k]
                  + f_1 * pc_y[k] * sos_4[k];

        t_14[k] = f_3 * sns_1[k]
                  + f_1 * pc_z[k] * sos_4[k];

        t_15[k] = f_4 * sns_5[k]
                  + f_1 * pc_x[k] * sos_5[k];

        t_16[k] = f_1 * pc_y[k] * sos_5[k];

        t_17[k] = f_5 * sns_2[k]
                  + f_1 * pc_z[k] * sos_5[k];

        t_18[k] = f_6 * sns_6[k]
                  + f_1 * pc_x[k] * sos_6[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, pc_x, pc_y, pc_z, sns_3, sns_4, \
                         sns_7, sns_8, sos_6, sos_7, sos_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * sns_3[k]
                  + f_1 * pc_y[k] * sos_6[k];

        t_20[k] = f_1 * pc_z[k] * sos_6[k];

        t_21[k] = f_6 * sns_7[k]
                  + f_1 * pc_x[k] * sos_7[k];

        t_22[k] = f_5 * sns_4[k]
                  + f_1 * pc_y[k] * sos_7[k];

        t_23[k] = f_3 * sns_3[k]
                  + f_1 * pc_z[k] * sos_7[k];

        t_24[k] = f_6 * sns_8[k]
                  + f_1 * pc_x[k] * sos_8[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, sns_4, sns_5, \
                         sns_9, sns_10, sos_8, sos_9, sos_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * sns_5[k]
                  + f_1 * pc_y[k] * sos_8[k];

        t_26[k] = f_5 * sns_4[k]
                  + f_1 * pc_z[k] * sos_8[k];

        t_27[k] = f_6 * sns_9[k]
                  + f_1 * pc_x[k] * sos_9[k];

        t_28[k] = f_1 * pc_y[k] * sos_9[k];

        t_29[k] = f_7 * sns_5[k]
                  + f_1 * pc_z[k] * sos_9[k];

        t_30[k] = f_8 * sns_10[k]
                  + f_1 * pc_x[k] * sos_10[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pc_x, pc_y, pc_z, sns_6, sns_7, \
                         sns_11, sns_12, sos_10, sos_11, sos_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_9 * sns_6[k]
                  + f_1 * pc_y[k] * sos_10[k];

        t_32[k] = f_1 * pc_z[k] * sos_10[k];

        t_33[k] = f_8 * sns_11[k]
                  + f_1 * pc_x[k] * sos_11[k];

        t_34[k] = f_7 * sns_7[k]
                  + f_1 * pc_y[k] * sos_11[k];

        t_35[k] = f_3 * sns_6[k]
                  + f_1 * pc_z[k] * sos_11[k];

        t_36[k] = f_8 * sns_12[k]
                  + f_1 * pc_x[k] * sos_12[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, sns_7, sns_8, sns_9, \
                         sns_13, sos_12, sos_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * sns_8[k]
                  + f_1 * pc_y[k] * sos_12[k];

        t_38[k] = f_5 * sns_7[k]
                  + f_1 * pc_z[k] * sos_12[k];

        t_39[k] = f_8 * sns_13[k]
                  + f_1 * pc_x[k] * sos_13[k];

        t_40[k] = f_3 * sns_9[k]
                  + f_1 * pc_y[k] * sos_13[k];

        t_41[k] = f_7 * sns_8[k]
                  + f_1 * pc_z[k] * sos_13[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pc_x, pc_y, pc_z, sns_9, sns_10, \
                         sns_14, sns_15, sos_14, sos_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_8 * sns_14[k]
                  + f_1 * pc_x[k] * sos_14[k];

        t_43[k] = f_1 * pc_y[k] * sos_14[k];

        t_44[k] = f_9 * sns_9[k]
                  + f_1 * pc_z[k] * sos_14[k];

        t_45[k] = f_10 * sns_15[k]
                  + f_1 * pc_x[k] * sos_15[k];

        t_46[k] = f_11 * sns_10[k]
                  + f_1 * pc_y[k] * sos_15[k];

        t_47[k] = f_1 * pc_z[k] * sos_15[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, pc_x, pc_y, pc_z, sns_10, sns_11, \
                         sns_12, sns_16, sns_17, sos_16, sos_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_10 * sns_16[k]
                  + f_1 * pc_x[k] * sos_16[k];

        t_49[k] = f_9 * sns_11[k]
                  + f_1 * pc_y[k] * sos_16[k];

        t_50[k] = f_3 * sns_10[k]
                  + f_1 * pc_z[k] * sos_16[k];

        t_51[k] = f_10 * sns_17[k]
                  + f_1 * pc_x[k] * sos_17[k];

        t_52[k] = f_7 * sns_12[k]
                  + f_1 * pc_y[k] * sos_17[k];

        t_53[k] = f_5 * sns_11[k]
                  + f_1 * pc_z[k] * sos_17[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, sns_12, sns_13, \
                         sns_14, sns_18, sns_19, sos_18, sos_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_10 * sns_18[k]
                  + f_1 * pc_x[k] * sos_18[k];

        t_55[k] = f_5 * sns_13[k]
                  + f_1 * pc_y[k] * sos_18[k];

        t_56[k] = f_7 * sns_12[k]
                  + f_1 * pc_z[k] * sos_18[k];

        t_57[k] = f_10 * sns_19[k]
                  + f_1 * pc_x[k] * sos_19[k];

        t_58[k] = f_3 * sns_14[k]
                  + f_1 * pc_y[k] * sos_19[k];

        t_59[k] = f_9 * sns_13[k]
                  + f_1 * pc_z[k] * sos_19[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, sns_14, sns_15, \
                         sns_20, sns_21, sos_20, sos_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_10 * sns_20[k]
                  + f_1 * pc_x[k] * sos_20[k];

        t_61[k] = f_1 * pc_y[k] * sos_20[k];

        t_62[k] = f_11 * sns_14[k]
                  + f_1 * pc_z[k] * sos_20[k];

        t_63[k] = f_11 * sns_21[k]
                  + f_1 * pc_x[k] * sos_21[k];

        t_64[k] = f_10 * sns_15[k]
                  + f_1 * pc_y[k] * sos_21[k];

        t_65[k] = f_1 * pc_z[k] * sos_21[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, pc_x, pc_y, pc_z, sns_15, sns_16, \
                         sns_17, sns_22, sns_23, sos_22, sos_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_11 * sns_22[k]
                  + f_1 * pc_x[k] * sos_22[k];

        t_67[k] = f_11 * sns_16[k]
                  + f_1 * pc_y[k] * sos_22[k];

        t_68[k] = f_3 * sns_15[k]
                  + f_1 * pc_z[k] * sos_22[k];

        t_69[k] = f_11 * sns_23[k]
                  + f_1 * pc_x[k] * sos_23[k];

        t_70[k] = f_9 * sns_17[k]
                  + f_1 * pc_y[k] * sos_23[k];

        t_71[k] = f_5 * sns_16[k]
                  + f_1 * pc_z[k] * sos_23[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, pc_x, pc_y, pc_z, sns_17, sns_18, \
                         sns_19, sns_24, sns_25, sos_24, sos_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_11 * sns_24[k]
                  + f_1 * pc_x[k] * sos_24[k];

        t_73[k] = f_7 * sns_18[k]
                  + f_1 * pc_y[k] * sos_24[k];

        t_74[k] = f_7 * sns_17[k]
                  + f_1 * pc_z[k] * sos_24[k];

        t_75[k] = f_11 * sns_25[k]
                  + f_1 * pc_x[k] * sos_25[k];

        t_76[k] = f_5 * sns_19[k]
                  + f_1 * pc_y[k] * sos_25[k];

        t_77[k] = f_9 * sns_18[k]
                  + f_1 * pc_z[k] * sos_25[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, pc_x, pc_y, pc_z, sns_19, sns_20, \
                         sns_26, sns_27, sos_26, sos_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_11 * sns_26[k]
                  + f_1 * pc_x[k] * sos_26[k];

        t_79[k] = f_3 * sns_20[k]
                  + f_1 * pc_y[k] * sos_26[k];

        t_80[k] = f_11 * sns_19[k]
                  + f_1 * pc_z[k] * sos_26[k];

        t_81[k] = f_11 * sns_27[k]
                  + f_1 * pc_x[k] * sos_27[k];

        t_82[k] = f_1 * pc_y[k] * sos_27[k];

        t_83[k] = f_10 * sns_20[k]
                  + f_1 * pc_z[k] * sos_27[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, sns_21, sns_22, \
                         sns_28, sns_29, sos_28, sos_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_9 * sns_28[k]
                  + f_1 * pc_x[k] * sos_28[k];

        t_85[k] = f_8 * sns_21[k]
                  + f_1 * pc_y[k] * sos_28[k];

        t_86[k] = f_1 * pc_z[k] * sos_28[k];

        t_87[k] = f_9 * sns_29[k]
                  + f_1 * pc_x[k] * sos_29[k];

        t_88[k] = f_10 * sns_22[k]
                  + f_1 * pc_y[k] * sos_29[k];

        t_89[k] = f_3 * sns_21[k]
                  + f_1 * pc_z[k] * sos_29[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, sns_22, sns_23, \
                         sns_24, sns_30, sns_31, sos_30, sos_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_9 * sns_30[k]
                  + f_1 * pc_x[k] * sos_30[k];

        t_91[k] = f_11 * sns_23[k]
                  + f_1 * pc_y[k] * sos_30[k];

        t_92[k] = f_5 * sns_22[k]
                  + f_1 * pc_z[k] * sos_30[k];

        t_93[k] = f_9 * sns_31[k]
                  + f_1 * pc_x[k] * sos_31[k];

        t_94[k] = f_9 * sns_24[k]
                  + f_1 * pc_y[k] * sos_31[k];

        t_95[k] = f_7 * sns_23[k]
                  + f_1 * pc_z[k] * sos_31[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, t_101, pc_x, pc_y, pc_z, sns_24, \
                         sns_25, sns_26, sns_32, sns_33, sos_32, \
                         sos_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_9 * sns_32[k]
                  + f_1 * pc_x[k] * sos_32[k];

        t_97[k] = f_7 * sns_25[k]
                  + f_1 * pc_y[k] * sos_32[k];

        t_98[k] = f_9 * sns_24[k]
                  + f_1 * pc_z[k] * sos_32[k];

        t_99[k] = f_9 * sns_33[k]
                  + f_1 * pc_x[k] * sos_33[k];

        t_100[k] = f_5 * sns_26[k]
                   + f_1 * pc_y[k] * sos_33[k];

        t_101[k] = f_11 * sns_25[k]
                   + f_1 * pc_z[k] * sos_33[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, t_107, pc_x, pc_y, pc_z, sns_26, \
                         sns_27, sns_34, sns_35, sos_34, sos_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_9 * sns_34[k]
                   + f_1 * pc_x[k] * sos_34[k];

        t_103[k] = f_3 * sns_27[k]
                   + f_1 * pc_y[k] * sos_34[k];

        t_104[k] = f_10 * sns_26[k]
                   + f_1 * pc_z[k] * sos_34[k];

        t_105[k] = f_9 * sns_35[k]
                   + f_1 * pc_x[k] * sos_35[k];

        t_106[k] = f_1 * pc_y[k] * sos_35[k];

        t_107[k] = f_8 * sns_27[k]
                   + f_1 * pc_z[k] * sos_35[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, t_113, pc_x, pc_y, pc_z, sns_28, \
                         sns_29, sns_36, sns_37, sos_36, sos_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_7 * sns_36[k]
                   + f_1 * pc_x[k] * sos_36[k];

        t_109[k] = f_6 * sns_28[k]
                   + f_1 * pc_y[k] * sos_36[k];

        t_110[k] = f_1 * pc_z[k] * sos_36[k];

        t_111[k] = f_7 * sns_37[k]
                   + f_1 * pc_x[k] * sos_37[k];

        t_112[k] = f_8 * sns_29[k]
                   + f_1 * pc_y[k] * sos_37[k];

        t_113[k] = f_3 * sns_28[k]
                   + f_1 * pc_z[k] * sos_37[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, pc_x, pc_y, pc_z, sns_29, \
                         sns_30, sns_31, sns_38, sns_39, sos_38, \
                         sos_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_7 * sns_38[k]
                   + f_1 * pc_x[k] * sos_38[k];

        t_115[k] = f_10 * sns_30[k]
                   + f_1 * pc_y[k] * sos_38[k];

        t_116[k] = f_5 * sns_29[k]
                   + f_1 * pc_z[k] * sos_38[k];

        t_117[k] = f_7 * sns_39[k]
                   + f_1 * pc_x[k] * sos_39[k];

        t_118[k] = f_11 * sns_31[k]
                   + f_1 * pc_y[k] * sos_39[k];

        t_119[k] = f_7 * sns_30[k]
                   + f_1 * pc_z[k] * sos_39[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, t_125, pc_x, pc_y, pc_z, sns_31, \
                         sns_32, sns_33, sns_40, sns_41, sos_40, \
                         sos_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_7 * sns_40[k]
                   + f_1 * pc_x[k] * sos_40[k];

        t_121[k] = f_9 * sns_32[k]
                   + f_1 * pc_y[k] * sos_40[k];

        t_122[k] = f_9 * sns_31[k]
                   + f_1 * pc_z[k] * sos_40[k];

        t_123[k] = f_7 * sns_41[k]
                   + f_1 * pc_x[k] * sos_41[k];

        t_124[k] = f_7 * sns_33[k]
                   + f_1 * pc_y[k] * sos_41[k];

        t_125[k] = f_11 * sns_32[k]
                   + f_1 * pc_z[k] * sos_41[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, t_131, pc_x, pc_y, pc_z, sns_33, \
                         sns_34, sns_35, sns_42, sns_43, sos_42, \
                         sos_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_7 * sns_42[k]
                   + f_1 * pc_x[k] * sos_42[k];

        t_127[k] = f_5 * sns_34[k]
                   + f_1 * pc_y[k] * sos_42[k];

        t_128[k] = f_10 * sns_33[k]
                   + f_1 * pc_z[k] * sos_42[k];

        t_129[k] = f_7 * sns_43[k]
                   + f_1 * pc_x[k] * sos_43[k];

        t_130[k] = f_3 * sns_35[k]
                   + f_1 * pc_y[k] * sos_43[k];

        t_131[k] = f_8 * sns_34[k]
                   + f_1 * pc_z[k] * sos_43[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, sns_35, \
                         sns_36, sns_44, sns_45, sos_44, sos_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_7 * sns_44[k]
                   + f_1 * pc_x[k] * sos_44[k];

        t_133[k] = f_1 * pc_y[k] * sos_44[k];

        t_134[k] = f_6 * sns_35[k]
                   + f_1 * pc_z[k] * sos_44[k];

        t_135[k] = f_5 * sns_45[k]
                   + f_1 * pc_x[k] * sos_45[k];

        t_136[k] = f_4 * sns_36[k]
                   + f_1 * pc_y[k] * sos_45[k];

        t_137[k] = f_1 * pc_z[k] * sos_45[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, t_143, pc_x, pc_y, pc_z, sns_36, \
                         sns_37, sns_38, sns_46, sns_47, sos_46, \
                         sos_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_5 * sns_46[k]
                   + f_1 * pc_x[k] * sos_46[k];

        t_139[k] = f_6 * sns_37[k]
                   + f_1 * pc_y[k] * sos_46[k];

        t_140[k] = f_3 * sns_36[k]
                   + f_1 * pc_z[k] * sos_46[k];

        t_141[k] = f_5 * sns_47[k]
                   + f_1 * pc_x[k] * sos_47[k];

        t_142[k] = f_8 * sns_38[k]
                   + f_1 * pc_y[k] * sos_47[k];

        t_143[k] = f_5 * sns_37[k]
                   + f_1 * pc_z[k] * sos_47[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, t_149, pc_x, pc_y, pc_z, sns_38, \
                         sns_39, sns_40, sns_48, sns_49, sos_48, \
                         sos_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_5 * sns_48[k]
                   + f_1 * pc_x[k] * sos_48[k];

        t_145[k] = f_10 * sns_39[k]
                   + f_1 * pc_y[k] * sos_48[k];

        t_146[k] = f_7 * sns_38[k]
                   + f_1 * pc_z[k] * sos_48[k];

        t_147[k] = f_5 * sns_49[k]
                   + f_1 * pc_x[k] * sos_49[k];

        t_148[k] = f_11 * sns_40[k]
                   + f_1 * pc_y[k] * sos_49[k];

        t_149[k] = f_9 * sns_39[k]
                   + f_1 * pc_z[k] * sos_49[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, t_155, pc_x, pc_y, pc_z, sns_40, \
                         sns_41, sns_42, sns_50, sns_51, sos_50, \
                         sos_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_5 * sns_50[k]
                   + f_1 * pc_x[k] * sos_50[k];

        t_151[k] = f_9 * sns_41[k]
                   + f_1 * pc_y[k] * sos_50[k];

        t_152[k] = f_11 * sns_40[k]
                   + f_1 * pc_z[k] * sos_50[k];

        t_153[k] = f_5 * sns_51[k]
                   + f_1 * pc_x[k] * sos_51[k];

        t_154[k] = f_7 * sns_42[k]
                   + f_1 * pc_y[k] * sos_51[k];

        t_155[k] = f_10 * sns_41[k]
                   + f_1 * pc_z[k] * sos_51[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, pc_x, pc_y, pc_z, sns_42, \
                         sns_43, sns_44, sns_52, sns_53, sos_52, \
                         sos_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_5 * sns_52[k]
                   + f_1 * pc_x[k] * sos_52[k];

        t_157[k] = f_5 * sns_43[k]
                   + f_1 * pc_y[k] * sos_52[k];

        t_158[k] = f_8 * sns_42[k]
                   + f_1 * pc_z[k] * sos_52[k];

        t_159[k] = f_5 * sns_53[k]
                   + f_1 * pc_x[k] * sos_53[k];

        t_160[k] = f_3 * sns_44[k]
                   + f_1 * pc_y[k] * sos_53[k];

        t_161[k] = f_6 * sns_43[k]
                   + f_1 * pc_z[k] * sos_53[k];
    }
}

static auto
compute_prim_sop_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t sns, const size_t sos,
                                                          const size_t ncols, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = p / q;
    const auto f_2 = 5.0 / q;
    const auto f_3 = 0.5 / q;
    const auto f_4 = 4.5 / q;
    const auto f_5 = 1.0 / q;
    const auto f_6 = 4.0 / q;
    const auto f_7 = 1.5 / q;
    const auto f_8 = 3.5 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 3.0 / q;
    const auto f_11 = 2.5 / q;

    auto *t_162 = buffer.data(target + 162);
    auto *t_163 = buffer.data(target + 163);
    auto *t_164 = buffer.data(target + 164);
    auto *t_165 = buffer.data(target + 165);
    auto *t_166 = buffer.data(target + 166);
    auto *t_167 = buffer.data(target + 167);
    auto *t_168 = buffer.data(target + 168);
    auto *t_169 = buffer.data(target + 169);
    auto *t_170 = buffer.data(target + 170);
    auto *t_171 = buffer.data(target + 171);
    auto *t_172 = buffer.data(target + 172);
    auto *t_173 = buffer.data(target + 173);
    auto *t_174 = buffer.data(target + 174);
    auto *t_175 = buffer.data(target + 175);
    auto *t_176 = buffer.data(target + 176);
    auto *t_177 = buffer.data(target + 177);
    auto *t_178 = buffer.data(target + 178);
    auto *t_179 = buffer.data(target + 179);
    auto *t_180 = buffer.data(target + 180);
    auto *t_181 = buffer.data(target + 181);
    auto *t_182 = buffer.data(target + 182);
    auto *t_183 = buffer.data(target + 183);
    auto *t_184 = buffer.data(target + 184);
    auto *t_185 = buffer.data(target + 185);
    auto *t_186 = buffer.data(target + 186);
    auto *t_187 = buffer.data(target + 187);
    auto *t_188 = buffer.data(target + 188);
    auto *t_189 = buffer.data(target + 189);
    auto *t_190 = buffer.data(target + 190);
    auto *t_191 = buffer.data(target + 191);
    auto *t_192 = buffer.data(target + 192);
    auto *t_193 = buffer.data(target + 193);
    auto *t_194 = buffer.data(target + 194);
    auto *t_195 = buffer.data(target + 195);
    auto *t_196 = buffer.data(target + 196);
    auto *t_197 = buffer.data(target + 197);
    auto *t_198 = buffer.data(target + 198);
    auto *t_199 = buffer.data(target + 199);
    auto *t_200 = buffer.data(target + 200);
    auto *t_201 = buffer.data(target + 201);
    auto *t_202 = buffer.data(target + 202);
    auto *t_203 = buffer.data(target + 203);
    auto *t_204 = buffer.data(target + 204);
    auto *t_205 = buffer.data(target + 205);
    auto *t_206 = buffer.data(target + 206);
    auto *t_207 = buffer.data(target + 207);
    auto *t_208 = buffer.data(target + 208);
    auto *t_209 = buffer.data(target + 209);
    auto *t_210 = buffer.data(target + 210);
    auto *t_211 = buffer.data(target + 211);
    auto *t_212 = buffer.data(target + 212);
    auto *t_213 = buffer.data(target + 213);
    auto *t_214 = buffer.data(target + 214);
    auto *t_215 = buffer.data(target + 215);
    auto *t_216 = buffer.data(target + 216);
    auto *t_217 = buffer.data(target + 217);
    auto *t_218 = buffer.data(target + 218);
    auto *t_219 = buffer.data(target + 219);
    auto *t_220 = buffer.data(target + 220);
    auto *t_221 = buffer.data(target + 221);
    auto *t_222 = buffer.data(target + 222);
    auto *t_223 = buffer.data(target + 223);
    auto *t_224 = buffer.data(target + 224);
    auto *t_225 = buffer.data(target + 225);
    auto *t_226 = buffer.data(target + 226);
    auto *t_227 = buffer.data(target + 227);
    auto *t_228 = buffer.data(target + 228);
    auto *t_229 = buffer.data(target + 229);
    auto *t_230 = buffer.data(target + 230);
    auto *t_231 = buffer.data(target + 231);
    auto *t_232 = buffer.data(target + 232);
    auto *t_233 = buffer.data(target + 233);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sns_44 = buffer.data(sns + 44);
    const auto *sns_45 = buffer.data(sns + 45);
    const auto *sns_46 = buffer.data(sns + 46);
    const auto *sns_47 = buffer.data(sns + 47);
    const auto *sns_48 = buffer.data(sns + 48);
    const auto *sns_49 = buffer.data(sns + 49);
    const auto *sns_50 = buffer.data(sns + 50);
    const auto *sns_51 = buffer.data(sns + 51);
    const auto *sns_52 = buffer.data(sns + 52);
    const auto *sns_53 = buffer.data(sns + 53);
    const auto *sns_54 = buffer.data(sns + 54);
    const auto *sns_55 = buffer.data(sns + 55);
    const auto *sns_56 = buffer.data(sns + 56);
    const auto *sns_57 = buffer.data(sns + 57);
    const auto *sns_58 = buffer.data(sns + 58);
    const auto *sns_59 = buffer.data(sns + 59);
    const auto *sns_60 = buffer.data(sns + 60);
    const auto *sns_61 = buffer.data(sns + 61);
    const auto *sns_62 = buffer.data(sns + 62);
    const auto *sns_63 = buffer.data(sns + 63);
    const auto *sns_64 = buffer.data(sns + 64);
    const auto *sns_65 = buffer.data(sns + 65);

    const auto *sos_54 = buffer.data(sos + 54);
    const auto *sos_55 = buffer.data(sos + 55);
    const auto *sos_56 = buffer.data(sos + 56);
    const auto *sos_57 = buffer.data(sos + 57);
    const auto *sos_58 = buffer.data(sos + 58);
    const auto *sos_59 = buffer.data(sos + 59);
    const auto *sos_60 = buffer.data(sos + 60);
    const auto *sos_61 = buffer.data(sos + 61);
    const auto *sos_62 = buffer.data(sos + 62);
    const auto *sos_63 = buffer.data(sos + 63);
    const auto *sos_64 = buffer.data(sos + 64);
    const auto *sos_65 = buffer.data(sos + 65);
    const auto *sos_66 = buffer.data(sos + 66);
    const auto *sos_67 = buffer.data(sos + 67);
    const auto *sos_68 = buffer.data(sos + 68);
    const auto *sos_69 = buffer.data(sos + 69);
    const auto *sos_70 = buffer.data(sos + 70);
    const auto *sos_71 = buffer.data(sos + 71);
    const auto *sos_72 = buffer.data(sos + 72);
    const auto *sos_73 = buffer.data(sos + 73);
    const auto *sos_74 = buffer.data(sos + 74);
    const auto *sos_75 = buffer.data(sos + 75);
    const auto *sos_76 = buffer.data(sos + 76);
    const auto *sos_77 = buffer.data(sos + 77);

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, t_167, pc_x, pc_y, pc_z, sns_44, \
                         sns_45, sns_54, sns_55, sos_54, sos_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_5 * sns_54[k]
                   + f_1 * pc_x[k] * sos_54[k];

        t_163[k] = f_1 * pc_y[k] * sos_54[k];

        t_164[k] = f_4 * sns_44[k]
                   + f_1 * pc_z[k] * sos_54[k];

        t_165[k] = f_3 * sns_55[k]
                   + f_1 * pc_x[k] * sos_55[k];

        t_166[k] = f_2 * sns_45[k]
                   + f_1 * pc_y[k] * sos_55[k];

        t_167[k] = f_1 * pc_z[k] * sos_55[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, t_173, pc_x, pc_y, pc_z, sns_45, \
                         sns_46, sns_47, sns_56, sns_57, sos_56, \
                         sos_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_3 * sns_56[k]
                   + f_1 * pc_x[k] * sos_56[k];

        t_169[k] = f_4 * sns_46[k]
                   + f_1 * pc_y[k] * sos_56[k];

        t_170[k] = f_3 * sns_45[k]
                   + f_1 * pc_z[k] * sos_56[k];

        t_171[k] = f_3 * sns_57[k]
                   + f_1 * pc_x[k] * sos_57[k];

        t_172[k] = f_6 * sns_47[k]
                   + f_1 * pc_y[k] * sos_57[k];

        t_173[k] = f_5 * sns_46[k]
                   + f_1 * pc_z[k] * sos_57[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, t_179, pc_x, pc_y, pc_z, sns_47, \
                         sns_48, sns_49, sns_58, sns_59, sos_58, \
                         sos_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_3 * sns_58[k]
                   + f_1 * pc_x[k] * sos_58[k];

        t_175[k] = f_8 * sns_48[k]
                   + f_1 * pc_y[k] * sos_58[k];

        t_176[k] = f_7 * sns_47[k]
                   + f_1 * pc_z[k] * sos_58[k];

        t_177[k] = f_3 * sns_59[k]
                   + f_1 * pc_x[k] * sos_59[k];

        t_178[k] = f_10 * sns_49[k]
                   + f_1 * pc_y[k] * sos_59[k];

        t_179[k] = f_9 * sns_48[k]
                   + f_1 * pc_z[k] * sos_59[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, t_185, pc_x, pc_y, pc_z, sns_49, \
                         sns_50, sns_51, sns_60, sns_61, sos_60, \
                         sos_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_3 * sns_60[k]
                   + f_1 * pc_x[k] * sos_60[k];

        t_181[k] = f_11 * sns_50[k]
                   + f_1 * pc_y[k] * sos_60[k];

        t_182[k] = f_11 * sns_49[k]
                   + f_1 * pc_z[k] * sos_60[k];

        t_183[k] = f_3 * sns_61[k]
                   + f_1 * pc_x[k] * sos_61[k];

        t_184[k] = f_9 * sns_51[k]
                   + f_1 * pc_y[k] * sos_61[k];

        t_185[k] = f_10 * sns_50[k]
                   + f_1 * pc_z[k] * sos_61[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, t_191, pc_x, pc_y, pc_z, sns_51, \
                         sns_52, sns_53, sns_62, sns_63, sos_62, \
                         sos_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_3 * sns_62[k]
                   + f_1 * pc_x[k] * sos_62[k];

        t_187[k] = f_7 * sns_52[k]
                   + f_1 * pc_y[k] * sos_62[k];

        t_188[k] = f_8 * sns_51[k]
                   + f_1 * pc_z[k] * sos_62[k];

        t_189[k] = f_3 * sns_63[k]
                   + f_1 * pc_x[k] * sos_63[k];

        t_190[k] = f_5 * sns_53[k]
                   + f_1 * pc_y[k] * sos_63[k];

        t_191[k] = f_6 * sns_52[k]
                   + f_1 * pc_z[k] * sos_63[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, t_197, pc_x, pc_y, pc_z, sns_53, \
                         sns_54, sns_64, sns_65, sos_64, sos_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_3 * sns_64[k]
                   + f_1 * pc_x[k] * sos_64[k];

        t_193[k] = f_3 * sns_54[k]
                   + f_1 * pc_y[k] * sos_64[k];

        t_194[k] = f_4 * sns_53[k]
                   + f_1 * pc_z[k] * sos_64[k];

        t_195[k] = f_3 * sns_65[k]
                   + f_1 * pc_x[k] * sos_65[k];

        t_196[k] = f_1 * pc_y[k] * sos_65[k];

        t_197[k] = f_2 * sns_54[k]
                   + f_1 * pc_z[k] * sos_65[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, t_203, t_204, pc_x, pc_y, pc_z, \
                         sns_55, sns_56, sos_66, sos_67, sos_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_1 * pc_x[k] * sos_66[k];

        t_199[k] = f_0 * sns_55[k]
                   + f_1 * pc_y[k] * sos_66[k];

        t_200[k] = f_1 * pc_z[k] * sos_66[k];

        t_201[k] = f_1 * pc_x[k] * sos_67[k];

        t_202[k] = f_2 * sns_56[k]
                   + f_1 * pc_y[k] * sos_67[k];

        t_203[k] = f_3 * sns_55[k]
                   + f_1 * pc_z[k] * sos_67[k];

        t_204[k] = f_1 * pc_x[k] * sos_68[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, t_210, pc_x, pc_y, pc_z, sns_56, \
                         sns_57, sns_58, sos_68, sos_69, sos_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_4 * sns_57[k]
                   + f_1 * pc_y[k] * sos_68[k];

        t_206[k] = f_5 * sns_56[k]
                   + f_1 * pc_z[k] * sos_68[k];

        t_207[k] = f_1 * pc_x[k] * sos_69[k];

        t_208[k] = f_6 * sns_58[k]
                   + f_1 * pc_y[k] * sos_69[k];

        t_209[k] = f_7 * sns_57[k]
                   + f_1 * pc_z[k] * sos_69[k];

        t_210[k] = f_1 * pc_x[k] * sos_70[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, t_216, pc_x, pc_y, pc_z, sns_58, \
                         sns_59, sns_60, sos_70, sos_71, sos_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_8 * sns_59[k]
                   + f_1 * pc_y[k] * sos_70[k];

        t_212[k] = f_9 * sns_58[k]
                   + f_1 * pc_z[k] * sos_70[k];

        t_213[k] = f_1 * pc_x[k] * sos_71[k];

        t_214[k] = f_10 * sns_60[k]
                   + f_1 * pc_y[k] * sos_71[k];

        t_215[k] = f_11 * sns_59[k]
                   + f_1 * pc_z[k] * sos_71[k];

        t_216[k] = f_1 * pc_x[k] * sos_72[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, t_222, pc_x, pc_y, pc_z, sns_60, \
                         sns_61, sns_62, sos_72, sos_73, sos_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_11 * sns_61[k]
                   + f_1 * pc_y[k] * sos_72[k];

        t_218[k] = f_10 * sns_60[k]
                   + f_1 * pc_z[k] * sos_72[k];

        t_219[k] = f_1 * pc_x[k] * sos_73[k];

        t_220[k] = f_9 * sns_62[k]
                   + f_1 * pc_y[k] * sos_73[k];

        t_221[k] = f_8 * sns_61[k]
                   + f_1 * pc_z[k] * sos_73[k];

        t_222[k] = f_1 * pc_x[k] * sos_74[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, t_227, t_228, pc_x, pc_y, pc_z, sns_62, \
                         sns_63, sns_64, sos_74, sos_75, sos_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_7 * sns_63[k]
                   + f_1 * pc_y[k] * sos_74[k];

        t_224[k] = f_6 * sns_62[k]
                   + f_1 * pc_z[k] * sos_74[k];

        t_225[k] = f_1 * pc_x[k] * sos_75[k];

        t_226[k] = f_5 * sns_64[k]
                   + f_1 * pc_y[k] * sos_75[k];

        t_227[k] = f_4 * sns_63[k]
                   + f_1 * pc_z[k] * sos_75[k];

        t_228[k] = f_1 * pc_x[k] * sos_76[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, pc_x, pc_y, pc_z, sns_64, sns_65, \
                         sos_76, sos_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_3 * sns_65[k]
                   + f_1 * pc_y[k] * sos_76[k];

        t_230[k] = f_2 * sns_64[k]
                   + f_1 * pc_z[k] * sos_76[k];

        t_231[k] = f_1 * pc_x[k] * sos_77[k];

        t_232[k] = f_1 * pc_y[k] * sos_77[k];

        t_233[k] = f_0 * sns_65[k]
                   + f_1 * pc_z[k] * sos_77[k];
    }
}

auto
compute_prim_sop_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t sns,
                                                   const size_t sos, const size_t ncols,
                                                   const double p, const double q) -> void
{
    compute_prim_sop_three_center_electron_repulsion_0_piece0(buffer, target, pc, sns, sos,
                                                              ncols, p, q);

    compute_prim_sop_three_center_electron_repulsion_0_piece1(buffer, target, pc, sns, sos,
                                                              ncols, p, q);
}

}  // namespace simdt3ceri
