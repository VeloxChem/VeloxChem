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


#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_qsp_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t oss, const size_t qss,
                                                          const size_t ncols, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
    const auto f_1 = p / q;
    const auto f_2 = 5.5 / q;
    const auto f_3 = 0.5 / q;
    const auto f_4 = 5.0 / q;
    const auto f_5 = 1.0 / q;
    const auto f_6 = 4.5 / q;
    const auto f_7 = 1.5 / q;
    const auto f_8 = 4.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 3.5 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 3.0 / q;

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

    const auto *oss_0 = buffer.data(oss + 0);
    const auto *oss_1 = buffer.data(oss + 1);
    const auto *oss_2 = buffer.data(oss + 2);
    const auto *oss_3 = buffer.data(oss + 3);
    const auto *oss_4 = buffer.data(oss + 4);
    const auto *oss_5 = buffer.data(oss + 5);
    const auto *oss_6 = buffer.data(oss + 6);
    const auto *oss_7 = buffer.data(oss + 7);
    const auto *oss_8 = buffer.data(oss + 8);
    const auto *oss_9 = buffer.data(oss + 9);
    const auto *oss_10 = buffer.data(oss + 10);
    const auto *oss_11 = buffer.data(oss + 11);
    const auto *oss_12 = buffer.data(oss + 12);
    const auto *oss_13 = buffer.data(oss + 13);
    const auto *oss_14 = buffer.data(oss + 14);
    const auto *oss_15 = buffer.data(oss + 15);
    const auto *oss_16 = buffer.data(oss + 16);
    const auto *oss_17 = buffer.data(oss + 17);
    const auto *oss_18 = buffer.data(oss + 18);
    const auto *oss_19 = buffer.data(oss + 19);
    const auto *oss_20 = buffer.data(oss + 20);
    const auto *oss_21 = buffer.data(oss + 21);
    const auto *oss_22 = buffer.data(oss + 22);
    const auto *oss_23 = buffer.data(oss + 23);
    const auto *oss_24 = buffer.data(oss + 24);
    const auto *oss_25 = buffer.data(oss + 25);
    const auto *oss_26 = buffer.data(oss + 26);
    const auto *oss_27 = buffer.data(oss + 27);
    const auto *oss_28 = buffer.data(oss + 28);
    const auto *oss_29 = buffer.data(oss + 29);
    const auto *oss_30 = buffer.data(oss + 30);
    const auto *oss_31 = buffer.data(oss + 31);
    const auto *oss_32 = buffer.data(oss + 32);
    const auto *oss_33 = buffer.data(oss + 33);
    const auto *oss_34 = buffer.data(oss + 34);
    const auto *oss_35 = buffer.data(oss + 35);
    const auto *oss_36 = buffer.data(oss + 36);
    const auto *oss_37 = buffer.data(oss + 37);
    const auto *oss_38 = buffer.data(oss + 38);
    const auto *oss_39 = buffer.data(oss + 39);
    const auto *oss_40 = buffer.data(oss + 40);
    const auto *oss_41 = buffer.data(oss + 41);
    const auto *oss_42 = buffer.data(oss + 42);
    const auto *oss_43 = buffer.data(oss + 43);
    const auto *oss_44 = buffer.data(oss + 44);
    const auto *oss_45 = buffer.data(oss + 45);
    const auto *oss_46 = buffer.data(oss + 46);
    const auto *oss_47 = buffer.data(oss + 47);
    const auto *oss_48 = buffer.data(oss + 48);
    const auto *oss_49 = buffer.data(oss + 49);
    const auto *oss_50 = buffer.data(oss + 50);
    const auto *oss_51 = buffer.data(oss + 51);
    const auto *oss_52 = buffer.data(oss + 52);
    const auto *oss_53 = buffer.data(oss + 53);

    const auto *qss_0 = buffer.data(qss + 0);
    const auto *qss_1 = buffer.data(qss + 1);
    const auto *qss_2 = buffer.data(qss + 2);
    const auto *qss_3 = buffer.data(qss + 3);
    const auto *qss_4 = buffer.data(qss + 4);
    const auto *qss_5 = buffer.data(qss + 5);
    const auto *qss_6 = buffer.data(qss + 6);
    const auto *qss_7 = buffer.data(qss + 7);
    const auto *qss_8 = buffer.data(qss + 8);
    const auto *qss_9 = buffer.data(qss + 9);
    const auto *qss_10 = buffer.data(qss + 10);
    const auto *qss_11 = buffer.data(qss + 11);
    const auto *qss_12 = buffer.data(qss + 12);
    const auto *qss_13 = buffer.data(qss + 13);
    const auto *qss_14 = buffer.data(qss + 14);
    const auto *qss_15 = buffer.data(qss + 15);
    const auto *qss_16 = buffer.data(qss + 16);
    const auto *qss_17 = buffer.data(qss + 17);
    const auto *qss_18 = buffer.data(qss + 18);
    const auto *qss_19 = buffer.data(qss + 19);
    const auto *qss_20 = buffer.data(qss + 20);
    const auto *qss_21 = buffer.data(qss + 21);
    const auto *qss_22 = buffer.data(qss + 22);
    const auto *qss_23 = buffer.data(qss + 23);
    const auto *qss_24 = buffer.data(qss + 24);
    const auto *qss_25 = buffer.data(qss + 25);
    const auto *qss_26 = buffer.data(qss + 26);
    const auto *qss_27 = buffer.data(qss + 27);
    const auto *qss_28 = buffer.data(qss + 28);
    const auto *qss_29 = buffer.data(qss + 29);
    const auto *qss_30 = buffer.data(qss + 30);
    const auto *qss_31 = buffer.data(qss + 31);
    const auto *qss_32 = buffer.data(qss + 32);
    const auto *qss_33 = buffer.data(qss + 33);
    const auto *qss_34 = buffer.data(qss + 34);
    const auto *qss_35 = buffer.data(qss + 35);
    const auto *qss_36 = buffer.data(qss + 36);
    const auto *qss_37 = buffer.data(qss + 37);
    const auto *qss_38 = buffer.data(qss + 38);
    const auto *qss_39 = buffer.data(qss + 39);
    const auto *qss_40 = buffer.data(qss + 40);
    const auto *qss_41 = buffer.data(qss + 41);
    const auto *qss_42 = buffer.data(qss + 42);
    const auto *qss_43 = buffer.data(qss + 43);
    const auto *qss_44 = buffer.data(qss + 44);
    const auto *qss_45 = buffer.data(qss + 45);
    const auto *qss_46 = buffer.data(qss + 46);
    const auto *qss_47 = buffer.data(qss + 47);
    const auto *qss_48 = buffer.data(qss + 48);
    const auto *qss_49 = buffer.data(qss + 49);
    const auto *qss_50 = buffer.data(qss + 50);
    const auto *qss_51 = buffer.data(qss + 51);
    const auto *qss_52 = buffer.data(qss + 52);
    const auto *qss_53 = buffer.data(qss + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pc_x, pc_y, pc_z, oss_0, oss_1, \
                         oss_2, qss_0, qss_1, qss_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * oss_0[k]
                 + f_1 * pc_x[k] * qss_0[k];

        t_1[k] = f_1 * pc_y[k] * qss_0[k];

        t_2[k] = f_1 * pc_z[k] * qss_0[k];

        t_3[k] = f_2 * oss_1[k]
                 + f_1 * pc_x[k] * qss_1[k];

        t_4[k] = f_3 * oss_0[k]
                 + f_1 * pc_y[k] * qss_1[k];

        t_5[k] = f_1 * pc_z[k] * qss_1[k];

        t_6[k] = f_2 * oss_2[k]
                 + f_1 * pc_x[k] * qss_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pc_x, pc_y, pc_z, oss_0, oss_1, \
                         oss_3, oss_4, qss_2, qss_3, qss_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * pc_y[k] * qss_2[k];

        t_8[k] = f_3 * oss_0[k]
                 + f_1 * pc_z[k] * qss_2[k];

        t_9[k] = f_4 * oss_3[k]
                 + f_1 * pc_x[k] * qss_3[k];

        t_10[k] = f_5 * oss_1[k]
                  + f_1 * pc_y[k] * qss_3[k];

        t_11[k] = f_1 * pc_z[k] * qss_3[k];

        t_12[k] = f_4 * oss_4[k]
                  + f_1 * pc_x[k] * qss_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, oss_1, oss_2, \
                         oss_5, oss_6, qss_4, qss_5, qss_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * oss_2[k]
                  + f_1 * pc_y[k] * qss_4[k];

        t_14[k] = f_3 * oss_1[k]
                  + f_1 * pc_z[k] * qss_4[k];

        t_15[k] = f_4 * oss_5[k]
                  + f_1 * pc_x[k] * qss_5[k];

        t_16[k] = f_1 * pc_y[k] * qss_5[k];

        t_17[k] = f_5 * oss_2[k]
                  + f_1 * pc_z[k] * qss_5[k];

        t_18[k] = f_6 * oss_6[k]
                  + f_1 * pc_x[k] * qss_6[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, pc_x, pc_y, pc_z, oss_3, oss_4, \
                         oss_7, oss_8, qss_6, qss_7, qss_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * oss_3[k]
                  + f_1 * pc_y[k] * qss_6[k];

        t_20[k] = f_1 * pc_z[k] * qss_6[k];

        t_21[k] = f_6 * oss_7[k]
                  + f_1 * pc_x[k] * qss_7[k];

        t_22[k] = f_5 * oss_4[k]
                  + f_1 * pc_y[k] * qss_7[k];

        t_23[k] = f_3 * oss_3[k]
                  + f_1 * pc_z[k] * qss_7[k];

        t_24[k] = f_6 * oss_8[k]
                  + f_1 * pc_x[k] * qss_8[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, oss_4, oss_5, \
                         oss_9, oss_10, qss_8, qss_9, qss_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * oss_5[k]
                  + f_1 * pc_y[k] * qss_8[k];

        t_26[k] = f_5 * oss_4[k]
                  + f_1 * pc_z[k] * qss_8[k];

        t_27[k] = f_6 * oss_9[k]
                  + f_1 * pc_x[k] * qss_9[k];

        t_28[k] = f_1 * pc_y[k] * qss_9[k];

        t_29[k] = f_7 * oss_5[k]
                  + f_1 * pc_z[k] * qss_9[k];

        t_30[k] = f_8 * oss_10[k]
                  + f_1 * pc_x[k] * qss_10[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pc_x, pc_y, pc_z, oss_6, oss_7, \
                         oss_11, oss_12, qss_10, qss_11, qss_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_9 * oss_6[k]
                  + f_1 * pc_y[k] * qss_10[k];

        t_32[k] = f_1 * pc_z[k] * qss_10[k];

        t_33[k] = f_8 * oss_11[k]
                  + f_1 * pc_x[k] * qss_11[k];

        t_34[k] = f_7 * oss_7[k]
                  + f_1 * pc_y[k] * qss_11[k];

        t_35[k] = f_3 * oss_6[k]
                  + f_1 * pc_z[k] * qss_11[k];

        t_36[k] = f_8 * oss_12[k]
                  + f_1 * pc_x[k] * qss_12[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, oss_7, oss_8, oss_9, \
                         oss_13, qss_12, qss_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * oss_8[k]
                  + f_1 * pc_y[k] * qss_12[k];

        t_38[k] = f_5 * oss_7[k]
                  + f_1 * pc_z[k] * qss_12[k];

        t_39[k] = f_8 * oss_13[k]
                  + f_1 * pc_x[k] * qss_13[k];

        t_40[k] = f_3 * oss_9[k]
                  + f_1 * pc_y[k] * qss_13[k];

        t_41[k] = f_7 * oss_8[k]
                  + f_1 * pc_z[k] * qss_13[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pc_x, pc_y, pc_z, oss_9, oss_10, \
                         oss_14, oss_15, qss_14, qss_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_8 * oss_14[k]
                  + f_1 * pc_x[k] * qss_14[k];

        t_43[k] = f_1 * pc_y[k] * qss_14[k];

        t_44[k] = f_9 * oss_9[k]
                  + f_1 * pc_z[k] * qss_14[k];

        t_45[k] = f_10 * oss_15[k]
                  + f_1 * pc_x[k] * qss_15[k];

        t_46[k] = f_11 * oss_10[k]
                  + f_1 * pc_y[k] * qss_15[k];

        t_47[k] = f_1 * pc_z[k] * qss_15[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, pc_x, pc_y, pc_z, oss_10, oss_11, \
                         oss_12, oss_16, oss_17, qss_16, qss_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_10 * oss_16[k]
                  + f_1 * pc_x[k] * qss_16[k];

        t_49[k] = f_9 * oss_11[k]
                  + f_1 * pc_y[k] * qss_16[k];

        t_50[k] = f_3 * oss_10[k]
                  + f_1 * pc_z[k] * qss_16[k];

        t_51[k] = f_10 * oss_17[k]
                  + f_1 * pc_x[k] * qss_17[k];

        t_52[k] = f_7 * oss_12[k]
                  + f_1 * pc_y[k] * qss_17[k];

        t_53[k] = f_5 * oss_11[k]
                  + f_1 * pc_z[k] * qss_17[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, oss_12, oss_13, \
                         oss_14, oss_18, oss_19, qss_18, qss_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_10 * oss_18[k]
                  + f_1 * pc_x[k] * qss_18[k];

        t_55[k] = f_5 * oss_13[k]
                  + f_1 * pc_y[k] * qss_18[k];

        t_56[k] = f_7 * oss_12[k]
                  + f_1 * pc_z[k] * qss_18[k];

        t_57[k] = f_10 * oss_19[k]
                  + f_1 * pc_x[k] * qss_19[k];

        t_58[k] = f_3 * oss_14[k]
                  + f_1 * pc_y[k] * qss_19[k];

        t_59[k] = f_9 * oss_13[k]
                  + f_1 * pc_z[k] * qss_19[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, oss_14, oss_15, \
                         oss_20, oss_21, qss_20, qss_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_10 * oss_20[k]
                  + f_1 * pc_x[k] * qss_20[k];

        t_61[k] = f_1 * pc_y[k] * qss_20[k];

        t_62[k] = f_11 * oss_14[k]
                  + f_1 * pc_z[k] * qss_20[k];

        t_63[k] = f_12 * oss_21[k]
                  + f_1 * pc_x[k] * qss_21[k];

        t_64[k] = f_12 * oss_15[k]
                  + f_1 * pc_y[k] * qss_21[k];

        t_65[k] = f_1 * pc_z[k] * qss_21[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, pc_x, pc_y, pc_z, oss_15, oss_16, \
                         oss_17, oss_22, oss_23, qss_22, qss_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_12 * oss_22[k]
                  + f_1 * pc_x[k] * qss_22[k];

        t_67[k] = f_11 * oss_16[k]
                  + f_1 * pc_y[k] * qss_22[k];

        t_68[k] = f_3 * oss_15[k]
                  + f_1 * pc_z[k] * qss_22[k];

        t_69[k] = f_12 * oss_23[k]
                  + f_1 * pc_x[k] * qss_23[k];

        t_70[k] = f_9 * oss_17[k]
                  + f_1 * pc_y[k] * qss_23[k];

        t_71[k] = f_5 * oss_16[k]
                  + f_1 * pc_z[k] * qss_23[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, pc_x, pc_y, pc_z, oss_17, oss_18, \
                         oss_19, oss_24, oss_25, qss_24, qss_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_12 * oss_24[k]
                  + f_1 * pc_x[k] * qss_24[k];

        t_73[k] = f_7 * oss_18[k]
                  + f_1 * pc_y[k] * qss_24[k];

        t_74[k] = f_7 * oss_17[k]
                  + f_1 * pc_z[k] * qss_24[k];

        t_75[k] = f_12 * oss_25[k]
                  + f_1 * pc_x[k] * qss_25[k];

        t_76[k] = f_5 * oss_19[k]
                  + f_1 * pc_y[k] * qss_25[k];

        t_77[k] = f_9 * oss_18[k]
                  + f_1 * pc_z[k] * qss_25[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, pc_x, pc_y, pc_z, oss_19, oss_20, \
                         oss_26, oss_27, qss_26, qss_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_12 * oss_26[k]
                  + f_1 * pc_x[k] * qss_26[k];

        t_79[k] = f_3 * oss_20[k]
                  + f_1 * pc_y[k] * qss_26[k];

        t_80[k] = f_11 * oss_19[k]
                  + f_1 * pc_z[k] * qss_26[k];

        t_81[k] = f_12 * oss_27[k]
                  + f_1 * pc_x[k] * qss_27[k];

        t_82[k] = f_1 * pc_y[k] * qss_27[k];

        t_83[k] = f_12 * oss_20[k]
                  + f_1 * pc_z[k] * qss_27[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, oss_21, oss_22, \
                         oss_28, oss_29, qss_28, qss_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_11 * oss_28[k]
                  + f_1 * pc_x[k] * qss_28[k];

        t_85[k] = f_10 * oss_21[k]
                  + f_1 * pc_y[k] * qss_28[k];

        t_86[k] = f_1 * pc_z[k] * qss_28[k];

        t_87[k] = f_11 * oss_29[k]
                  + f_1 * pc_x[k] * qss_29[k];

        t_88[k] = f_12 * oss_22[k]
                  + f_1 * pc_y[k] * qss_29[k];

        t_89[k] = f_3 * oss_21[k]
                  + f_1 * pc_z[k] * qss_29[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, oss_22, oss_23, \
                         oss_24, oss_30, oss_31, qss_30, qss_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_11 * oss_30[k]
                  + f_1 * pc_x[k] * qss_30[k];

        t_91[k] = f_11 * oss_23[k]
                  + f_1 * pc_y[k] * qss_30[k];

        t_92[k] = f_5 * oss_22[k]
                  + f_1 * pc_z[k] * qss_30[k];

        t_93[k] = f_11 * oss_31[k]
                  + f_1 * pc_x[k] * qss_31[k];

        t_94[k] = f_9 * oss_24[k]
                  + f_1 * pc_y[k] * qss_31[k];

        t_95[k] = f_7 * oss_23[k]
                  + f_1 * pc_z[k] * qss_31[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, t_101, pc_x, pc_y, pc_z, oss_24, \
                         oss_25, oss_26, oss_32, oss_33, qss_32, \
                         qss_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_11 * oss_32[k]
                  + f_1 * pc_x[k] * qss_32[k];

        t_97[k] = f_7 * oss_25[k]
                  + f_1 * pc_y[k] * qss_32[k];

        t_98[k] = f_9 * oss_24[k]
                  + f_1 * pc_z[k] * qss_32[k];

        t_99[k] = f_11 * oss_33[k]
                  + f_1 * pc_x[k] * qss_33[k];

        t_100[k] = f_5 * oss_26[k]
                   + f_1 * pc_y[k] * qss_33[k];

        t_101[k] = f_11 * oss_25[k]
                   + f_1 * pc_z[k] * qss_33[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, t_107, pc_x, pc_y, pc_z, oss_26, \
                         oss_27, oss_34, oss_35, qss_34, qss_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_11 * oss_34[k]
                   + f_1 * pc_x[k] * qss_34[k];

        t_103[k] = f_3 * oss_27[k]
                   + f_1 * pc_y[k] * qss_34[k];

        t_104[k] = f_12 * oss_26[k]
                   + f_1 * pc_z[k] * qss_34[k];

        t_105[k] = f_11 * oss_35[k]
                   + f_1 * pc_x[k] * qss_35[k];

        t_106[k] = f_1 * pc_y[k] * qss_35[k];

        t_107[k] = f_10 * oss_27[k]
                   + f_1 * pc_z[k] * qss_35[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, t_113, pc_x, pc_y, pc_z, oss_28, \
                         oss_29, oss_36, oss_37, qss_36, qss_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_9 * oss_36[k]
                   + f_1 * pc_x[k] * qss_36[k];

        t_109[k] = f_8 * oss_28[k]
                   + f_1 * pc_y[k] * qss_36[k];

        t_110[k] = f_1 * pc_z[k] * qss_36[k];

        t_111[k] = f_9 * oss_37[k]
                   + f_1 * pc_x[k] * qss_37[k];

        t_112[k] = f_10 * oss_29[k]
                   + f_1 * pc_y[k] * qss_37[k];

        t_113[k] = f_3 * oss_28[k]
                   + f_1 * pc_z[k] * qss_37[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, pc_x, pc_y, pc_z, oss_29, \
                         oss_30, oss_31, oss_38, oss_39, qss_38, \
                         qss_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_9 * oss_38[k]
                   + f_1 * pc_x[k] * qss_38[k];

        t_115[k] = f_12 * oss_30[k]
                   + f_1 * pc_y[k] * qss_38[k];

        t_116[k] = f_5 * oss_29[k]
                   + f_1 * pc_z[k] * qss_38[k];

        t_117[k] = f_9 * oss_39[k]
                   + f_1 * pc_x[k] * qss_39[k];

        t_118[k] = f_11 * oss_31[k]
                   + f_1 * pc_y[k] * qss_39[k];

        t_119[k] = f_7 * oss_30[k]
                   + f_1 * pc_z[k] * qss_39[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, t_125, pc_x, pc_y, pc_z, oss_31, \
                         oss_32, oss_33, oss_40, oss_41, qss_40, \
                         qss_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_9 * oss_40[k]
                   + f_1 * pc_x[k] * qss_40[k];

        t_121[k] = f_9 * oss_32[k]
                   + f_1 * pc_y[k] * qss_40[k];

        t_122[k] = f_9 * oss_31[k]
                   + f_1 * pc_z[k] * qss_40[k];

        t_123[k] = f_9 * oss_41[k]
                   + f_1 * pc_x[k] * qss_41[k];

        t_124[k] = f_7 * oss_33[k]
                   + f_1 * pc_y[k] * qss_41[k];

        t_125[k] = f_11 * oss_32[k]
                   + f_1 * pc_z[k] * qss_41[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, t_131, pc_x, pc_y, pc_z, oss_33, \
                         oss_34, oss_35, oss_42, oss_43, qss_42, \
                         qss_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_9 * oss_42[k]
                   + f_1 * pc_x[k] * qss_42[k];

        t_127[k] = f_5 * oss_34[k]
                   + f_1 * pc_y[k] * qss_42[k];

        t_128[k] = f_12 * oss_33[k]
                   + f_1 * pc_z[k] * qss_42[k];

        t_129[k] = f_9 * oss_43[k]
                   + f_1 * pc_x[k] * qss_43[k];

        t_130[k] = f_3 * oss_35[k]
                   + f_1 * pc_y[k] * qss_43[k];

        t_131[k] = f_10 * oss_34[k]
                   + f_1 * pc_z[k] * qss_43[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, oss_35, \
                         oss_36, oss_44, oss_45, qss_44, qss_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_9 * oss_44[k]
                   + f_1 * pc_x[k] * qss_44[k];

        t_133[k] = f_1 * pc_y[k] * qss_44[k];

        t_134[k] = f_8 * oss_35[k]
                   + f_1 * pc_z[k] * qss_44[k];

        t_135[k] = f_7 * oss_45[k]
                   + f_1 * pc_x[k] * qss_45[k];

        t_136[k] = f_6 * oss_36[k]
                   + f_1 * pc_y[k] * qss_45[k];

        t_137[k] = f_1 * pc_z[k] * qss_45[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, t_143, pc_x, pc_y, pc_z, oss_36, \
                         oss_37, oss_38, oss_46, oss_47, qss_46, \
                         qss_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_7 * oss_46[k]
                   + f_1 * pc_x[k] * qss_46[k];

        t_139[k] = f_8 * oss_37[k]
                   + f_1 * pc_y[k] * qss_46[k];

        t_140[k] = f_3 * oss_36[k]
                   + f_1 * pc_z[k] * qss_46[k];

        t_141[k] = f_7 * oss_47[k]
                   + f_1 * pc_x[k] * qss_47[k];

        t_142[k] = f_10 * oss_38[k]
                   + f_1 * pc_y[k] * qss_47[k];

        t_143[k] = f_5 * oss_37[k]
                   + f_1 * pc_z[k] * qss_47[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, t_149, pc_x, pc_y, pc_z, oss_38, \
                         oss_39, oss_40, oss_48, oss_49, qss_48, \
                         qss_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_7 * oss_48[k]
                   + f_1 * pc_x[k] * qss_48[k];

        t_145[k] = f_12 * oss_39[k]
                   + f_1 * pc_y[k] * qss_48[k];

        t_146[k] = f_7 * oss_38[k]
                   + f_1 * pc_z[k] * qss_48[k];

        t_147[k] = f_7 * oss_49[k]
                   + f_1 * pc_x[k] * qss_49[k];

        t_148[k] = f_11 * oss_40[k]
                   + f_1 * pc_y[k] * qss_49[k];

        t_149[k] = f_9 * oss_39[k]
                   + f_1 * pc_z[k] * qss_49[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, t_155, pc_x, pc_y, pc_z, oss_40, \
                         oss_41, oss_42, oss_50, oss_51, qss_50, \
                         qss_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_7 * oss_50[k]
                   + f_1 * pc_x[k] * qss_50[k];

        t_151[k] = f_9 * oss_41[k]
                   + f_1 * pc_y[k] * qss_50[k];

        t_152[k] = f_11 * oss_40[k]
                   + f_1 * pc_z[k] * qss_50[k];

        t_153[k] = f_7 * oss_51[k]
                   + f_1 * pc_x[k] * qss_51[k];

        t_154[k] = f_7 * oss_42[k]
                   + f_1 * pc_y[k] * qss_51[k];

        t_155[k] = f_12 * oss_41[k]
                   + f_1 * pc_z[k] * qss_51[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, pc_x, pc_y, pc_z, oss_42, \
                         oss_43, oss_44, oss_52, oss_53, qss_52, \
                         qss_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_7 * oss_52[k]
                   + f_1 * pc_x[k] * qss_52[k];

        t_157[k] = f_5 * oss_43[k]
                   + f_1 * pc_y[k] * qss_52[k];

        t_158[k] = f_10 * oss_42[k]
                   + f_1 * pc_z[k] * qss_52[k];

        t_159[k] = f_7 * oss_53[k]
                   + f_1 * pc_x[k] * qss_53[k];

        t_160[k] = f_3 * oss_44[k]
                   + f_1 * pc_y[k] * qss_53[k];

        t_161[k] = f_8 * oss_43[k]
                   + f_1 * pc_z[k] * qss_53[k];
    }
}

static auto
compute_prim_qsp_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t oss, const size_t qss,
                                                          const size_t ncols, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 6.0 / q;
    const auto f_1 = p / q;
    const auto f_2 = 5.5 / q;
    const auto f_3 = 0.5 / q;
    const auto f_4 = 5.0 / q;
    const auto f_5 = 1.0 / q;
    const auto f_6 = 4.5 / q;
    const auto f_7 = 1.5 / q;
    const auto f_8 = 4.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 3.5 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 3.0 / q;

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
    auto *t_234 = buffer.data(target + 234);
    auto *t_235 = buffer.data(target + 235);
    auto *t_236 = buffer.data(target + 236);
    auto *t_237 = buffer.data(target + 237);
    auto *t_238 = buffer.data(target + 238);
    auto *t_239 = buffer.data(target + 239);
    auto *t_240 = buffer.data(target + 240);
    auto *t_241 = buffer.data(target + 241);
    auto *t_242 = buffer.data(target + 242);
    auto *t_243 = buffer.data(target + 243);
    auto *t_244 = buffer.data(target + 244);
    auto *t_245 = buffer.data(target + 245);
    auto *t_246 = buffer.data(target + 246);
    auto *t_247 = buffer.data(target + 247);
    auto *t_248 = buffer.data(target + 248);
    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);
    auto *t_255 = buffer.data(target + 255);
    auto *t_256 = buffer.data(target + 256);
    auto *t_257 = buffer.data(target + 257);
    auto *t_258 = buffer.data(target + 258);
    auto *t_259 = buffer.data(target + 259);
    auto *t_260 = buffer.data(target + 260);
    auto *t_261 = buffer.data(target + 261);
    auto *t_262 = buffer.data(target + 262);
    auto *t_263 = buffer.data(target + 263);
    auto *t_264 = buffer.data(target + 264);
    auto *t_265 = buffer.data(target + 265);
    auto *t_266 = buffer.data(target + 266);
    auto *t_267 = buffer.data(target + 267);
    auto *t_268 = buffer.data(target + 268);
    auto *t_269 = buffer.data(target + 269);
    auto *t_270 = buffer.data(target + 270);
    auto *t_271 = buffer.data(target + 271);
    auto *t_272 = buffer.data(target + 272);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *oss_44 = buffer.data(oss + 44);
    const auto *oss_45 = buffer.data(oss + 45);
    const auto *oss_46 = buffer.data(oss + 46);
    const auto *oss_47 = buffer.data(oss + 47);
    const auto *oss_48 = buffer.data(oss + 48);
    const auto *oss_49 = buffer.data(oss + 49);
    const auto *oss_50 = buffer.data(oss + 50);
    const auto *oss_51 = buffer.data(oss + 51);
    const auto *oss_52 = buffer.data(oss + 52);
    const auto *oss_53 = buffer.data(oss + 53);
    const auto *oss_54 = buffer.data(oss + 54);
    const auto *oss_55 = buffer.data(oss + 55);
    const auto *oss_56 = buffer.data(oss + 56);
    const auto *oss_57 = buffer.data(oss + 57);
    const auto *oss_58 = buffer.data(oss + 58);
    const auto *oss_59 = buffer.data(oss + 59);
    const auto *oss_60 = buffer.data(oss + 60);
    const auto *oss_61 = buffer.data(oss + 61);
    const auto *oss_62 = buffer.data(oss + 62);
    const auto *oss_63 = buffer.data(oss + 63);
    const auto *oss_64 = buffer.data(oss + 64);
    const auto *oss_65 = buffer.data(oss + 65);
    const auto *oss_66 = buffer.data(oss + 66);
    const auto *oss_67 = buffer.data(oss + 67);
    const auto *oss_68 = buffer.data(oss + 68);
    const auto *oss_69 = buffer.data(oss + 69);
    const auto *oss_70 = buffer.data(oss + 70);
    const auto *oss_71 = buffer.data(oss + 71);
    const auto *oss_72 = buffer.data(oss + 72);
    const auto *oss_73 = buffer.data(oss + 73);
    const auto *oss_74 = buffer.data(oss + 74);
    const auto *oss_75 = buffer.data(oss + 75);
    const auto *oss_76 = buffer.data(oss + 76);
    const auto *oss_77 = buffer.data(oss + 77);

    const auto *qss_54 = buffer.data(qss + 54);
    const auto *qss_55 = buffer.data(qss + 55);
    const auto *qss_56 = buffer.data(qss + 56);
    const auto *qss_57 = buffer.data(qss + 57);
    const auto *qss_58 = buffer.data(qss + 58);
    const auto *qss_59 = buffer.data(qss + 59);
    const auto *qss_60 = buffer.data(qss + 60);
    const auto *qss_61 = buffer.data(qss + 61);
    const auto *qss_62 = buffer.data(qss + 62);
    const auto *qss_63 = buffer.data(qss + 63);
    const auto *qss_64 = buffer.data(qss + 64);
    const auto *qss_65 = buffer.data(qss + 65);
    const auto *qss_66 = buffer.data(qss + 66);
    const auto *qss_67 = buffer.data(qss + 67);
    const auto *qss_68 = buffer.data(qss + 68);
    const auto *qss_69 = buffer.data(qss + 69);
    const auto *qss_70 = buffer.data(qss + 70);
    const auto *qss_71 = buffer.data(qss + 71);
    const auto *qss_72 = buffer.data(qss + 72);
    const auto *qss_73 = buffer.data(qss + 73);
    const auto *qss_74 = buffer.data(qss + 74);
    const auto *qss_75 = buffer.data(qss + 75);
    const auto *qss_76 = buffer.data(qss + 76);
    const auto *qss_77 = buffer.data(qss + 77);
    const auto *qss_78 = buffer.data(qss + 78);
    const auto *qss_79 = buffer.data(qss + 79);
    const auto *qss_80 = buffer.data(qss + 80);
    const auto *qss_81 = buffer.data(qss + 81);
    const auto *qss_82 = buffer.data(qss + 82);
    const auto *qss_83 = buffer.data(qss + 83);
    const auto *qss_84 = buffer.data(qss + 84);
    const auto *qss_85 = buffer.data(qss + 85);
    const auto *qss_86 = buffer.data(qss + 86);
    const auto *qss_87 = buffer.data(qss + 87);
    const auto *qss_88 = buffer.data(qss + 88);
    const auto *qss_89 = buffer.data(qss + 89);
    const auto *qss_90 = buffer.data(qss + 90);

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, t_167, pc_x, pc_y, pc_z, oss_44, \
                         oss_45, oss_54, oss_55, qss_54, qss_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_7 * oss_54[k]
                   + f_1 * pc_x[k] * qss_54[k];

        t_163[k] = f_1 * pc_y[k] * qss_54[k];

        t_164[k] = f_6 * oss_44[k]
                   + f_1 * pc_z[k] * qss_54[k];

        t_165[k] = f_5 * oss_55[k]
                   + f_1 * pc_x[k] * qss_55[k];

        t_166[k] = f_4 * oss_45[k]
                   + f_1 * pc_y[k] * qss_55[k];

        t_167[k] = f_1 * pc_z[k] * qss_55[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, t_173, pc_x, pc_y, pc_z, oss_45, \
                         oss_46, oss_47, oss_56, oss_57, qss_56, \
                         qss_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_5 * oss_56[k]
                   + f_1 * pc_x[k] * qss_56[k];

        t_169[k] = f_6 * oss_46[k]
                   + f_1 * pc_y[k] * qss_56[k];

        t_170[k] = f_3 * oss_45[k]
                   + f_1 * pc_z[k] * qss_56[k];

        t_171[k] = f_5 * oss_57[k]
                   + f_1 * pc_x[k] * qss_57[k];

        t_172[k] = f_8 * oss_47[k]
                   + f_1 * pc_y[k] * qss_57[k];

        t_173[k] = f_5 * oss_46[k]
                   + f_1 * pc_z[k] * qss_57[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, t_179, pc_x, pc_y, pc_z, oss_47, \
                         oss_48, oss_49, oss_58, oss_59, qss_58, \
                         qss_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_5 * oss_58[k]
                   + f_1 * pc_x[k] * qss_58[k];

        t_175[k] = f_10 * oss_48[k]
                   + f_1 * pc_y[k] * qss_58[k];

        t_176[k] = f_7 * oss_47[k]
                   + f_1 * pc_z[k] * qss_58[k];

        t_177[k] = f_5 * oss_59[k]
                   + f_1 * pc_x[k] * qss_59[k];

        t_178[k] = f_12 * oss_49[k]
                   + f_1 * pc_y[k] * qss_59[k];

        t_179[k] = f_9 * oss_48[k]
                   + f_1 * pc_z[k] * qss_59[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, t_185, pc_x, pc_y, pc_z, oss_49, \
                         oss_50, oss_51, oss_60, oss_61, qss_60, \
                         qss_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_5 * oss_60[k]
                   + f_1 * pc_x[k] * qss_60[k];

        t_181[k] = f_11 * oss_50[k]
                   + f_1 * pc_y[k] * qss_60[k];

        t_182[k] = f_11 * oss_49[k]
                   + f_1 * pc_z[k] * qss_60[k];

        t_183[k] = f_5 * oss_61[k]
                   + f_1 * pc_x[k] * qss_61[k];

        t_184[k] = f_9 * oss_51[k]
                   + f_1 * pc_y[k] * qss_61[k];

        t_185[k] = f_12 * oss_50[k]
                   + f_1 * pc_z[k] * qss_61[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, t_191, pc_x, pc_y, pc_z, oss_51, \
                         oss_52, oss_53, oss_62, oss_63, qss_62, \
                         qss_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_5 * oss_62[k]
                   + f_1 * pc_x[k] * qss_62[k];

        t_187[k] = f_7 * oss_52[k]
                   + f_1 * pc_y[k] * qss_62[k];

        t_188[k] = f_10 * oss_51[k]
                   + f_1 * pc_z[k] * qss_62[k];

        t_189[k] = f_5 * oss_63[k]
                   + f_1 * pc_x[k] * qss_63[k];

        t_190[k] = f_5 * oss_53[k]
                   + f_1 * pc_y[k] * qss_63[k];

        t_191[k] = f_8 * oss_52[k]
                   + f_1 * pc_z[k] * qss_63[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, t_197, pc_x, pc_y, pc_z, oss_53, \
                         oss_54, oss_64, oss_65, qss_64, qss_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_5 * oss_64[k]
                   + f_1 * pc_x[k] * qss_64[k];

        t_193[k] = f_3 * oss_54[k]
                   + f_1 * pc_y[k] * qss_64[k];

        t_194[k] = f_6 * oss_53[k]
                   + f_1 * pc_z[k] * qss_64[k];

        t_195[k] = f_5 * oss_65[k]
                   + f_1 * pc_x[k] * qss_65[k];

        t_196[k] = f_1 * pc_y[k] * qss_65[k];

        t_197[k] = f_4 * oss_54[k]
                   + f_1 * pc_z[k] * qss_65[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, t_203, pc_x, pc_y, pc_z, oss_55, \
                         oss_56, oss_66, oss_67, qss_66, qss_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_3 * oss_66[k]
                   + f_1 * pc_x[k] * qss_66[k];

        t_199[k] = f_2 * oss_55[k]
                   + f_1 * pc_y[k] * qss_66[k];

        t_200[k] = f_1 * pc_z[k] * qss_66[k];

        t_201[k] = f_3 * oss_67[k]
                   + f_1 * pc_x[k] * qss_67[k];

        t_202[k] = f_4 * oss_56[k]
                   + f_1 * pc_y[k] * qss_67[k];

        t_203[k] = f_3 * oss_55[k]
                   + f_1 * pc_z[k] * qss_67[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, t_209, pc_x, pc_y, pc_z, oss_56, \
                         oss_57, oss_58, oss_68, oss_69, qss_68, \
                         qss_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_3 * oss_68[k]
                   + f_1 * pc_x[k] * qss_68[k];

        t_205[k] = f_6 * oss_57[k]
                   + f_1 * pc_y[k] * qss_68[k];

        t_206[k] = f_5 * oss_56[k]
                   + f_1 * pc_z[k] * qss_68[k];

        t_207[k] = f_3 * oss_69[k]
                   + f_1 * pc_x[k] * qss_69[k];

        t_208[k] = f_8 * oss_58[k]
                   + f_1 * pc_y[k] * qss_69[k];

        t_209[k] = f_7 * oss_57[k]
                   + f_1 * pc_z[k] * qss_69[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, pc_x, pc_y, pc_z, oss_58, \
                         oss_59, oss_60, oss_70, oss_71, qss_70, \
                         qss_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_3 * oss_70[k]
                   + f_1 * pc_x[k] * qss_70[k];

        t_211[k] = f_10 * oss_59[k]
                   + f_1 * pc_y[k] * qss_70[k];

        t_212[k] = f_9 * oss_58[k]
                   + f_1 * pc_z[k] * qss_70[k];

        t_213[k] = f_3 * oss_71[k]
                   + f_1 * pc_x[k] * qss_71[k];

        t_214[k] = f_12 * oss_60[k]
                   + f_1 * pc_y[k] * qss_71[k];

        t_215[k] = f_11 * oss_59[k]
                   + f_1 * pc_z[k] * qss_71[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, t_221, pc_x, pc_y, pc_z, oss_60, \
                         oss_61, oss_62, oss_72, oss_73, qss_72, \
                         qss_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_3 * oss_72[k]
                   + f_1 * pc_x[k] * qss_72[k];

        t_217[k] = f_11 * oss_61[k]
                   + f_1 * pc_y[k] * qss_72[k];

        t_218[k] = f_12 * oss_60[k]
                   + f_1 * pc_z[k] * qss_72[k];

        t_219[k] = f_3 * oss_73[k]
                   + f_1 * pc_x[k] * qss_73[k];

        t_220[k] = f_9 * oss_62[k]
                   + f_1 * pc_y[k] * qss_73[k];

        t_221[k] = f_10 * oss_61[k]
                   + f_1 * pc_z[k] * qss_73[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, t_227, pc_x, pc_y, pc_z, oss_62, \
                         oss_63, oss_64, oss_74, oss_75, qss_74, \
                         qss_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_3 * oss_74[k]
                   + f_1 * pc_x[k] * qss_74[k];

        t_223[k] = f_7 * oss_63[k]
                   + f_1 * pc_y[k] * qss_74[k];

        t_224[k] = f_8 * oss_62[k]
                   + f_1 * pc_z[k] * qss_74[k];

        t_225[k] = f_3 * oss_75[k]
                   + f_1 * pc_x[k] * qss_75[k];

        t_226[k] = f_5 * oss_64[k]
                   + f_1 * pc_y[k] * qss_75[k];

        t_227[k] = f_6 * oss_63[k]
                   + f_1 * pc_z[k] * qss_75[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, t_232, t_233, pc_x, pc_y, pc_z, oss_64, \
                         oss_65, oss_76, oss_77, qss_76, qss_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_3 * oss_76[k]
                   + f_1 * pc_x[k] * qss_76[k];

        t_229[k] = f_3 * oss_65[k]
                   + f_1 * pc_y[k] * qss_76[k];

        t_230[k] = f_4 * oss_64[k]
                   + f_1 * pc_z[k] * qss_76[k];

        t_231[k] = f_3 * oss_77[k]
                   + f_1 * pc_x[k] * qss_77[k];

        t_232[k] = f_1 * pc_y[k] * qss_77[k];

        t_233[k] = f_2 * oss_65[k]
                   + f_1 * pc_z[k] * qss_77[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, t_239, t_240, pc_x, pc_y, pc_z, \
                         oss_66, oss_67, qss_78, qss_79, qss_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_1 * pc_x[k] * qss_78[k];

        t_235[k] = f_0 * oss_66[k]
                   + f_1 * pc_y[k] * qss_78[k];

        t_236[k] = f_1 * pc_z[k] * qss_78[k];

        t_237[k] = f_1 * pc_x[k] * qss_79[k];

        t_238[k] = f_2 * oss_67[k]
                   + f_1 * pc_y[k] * qss_79[k];

        t_239[k] = f_3 * oss_66[k]
                   + f_1 * pc_z[k] * qss_79[k];

        t_240[k] = f_1 * pc_x[k] * qss_80[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, t_246, pc_x, pc_y, pc_z, oss_67, \
                         oss_68, oss_69, qss_80, qss_81, qss_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_4 * oss_68[k]
                   + f_1 * pc_y[k] * qss_80[k];

        t_242[k] = f_5 * oss_67[k]
                   + f_1 * pc_z[k] * qss_80[k];

        t_243[k] = f_1 * pc_x[k] * qss_81[k];

        t_244[k] = f_6 * oss_69[k]
                   + f_1 * pc_y[k] * qss_81[k];

        t_245[k] = f_7 * oss_68[k]
                   + f_1 * pc_z[k] * qss_81[k];

        t_246[k] = f_1 * pc_x[k] * qss_82[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, t_252, pc_x, pc_y, pc_z, oss_69, \
                         oss_70, oss_71, qss_82, qss_83, qss_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_8 * oss_70[k]
                   + f_1 * pc_y[k] * qss_82[k];

        t_248[k] = f_9 * oss_69[k]
                   + f_1 * pc_z[k] * qss_82[k];

        t_249[k] = f_1 * pc_x[k] * qss_83[k];

        t_250[k] = f_10 * oss_71[k]
                   + f_1 * pc_y[k] * qss_83[k];

        t_251[k] = f_11 * oss_70[k]
                   + f_1 * pc_z[k] * qss_83[k];

        t_252[k] = f_1 * pc_x[k] * qss_84[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, t_257, t_258, pc_x, pc_y, pc_z, oss_71, \
                         oss_72, oss_73, qss_84, qss_85, qss_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_12 * oss_72[k]
                   + f_1 * pc_y[k] * qss_84[k];

        t_254[k] = f_12 * oss_71[k]
                   + f_1 * pc_z[k] * qss_84[k];

        t_255[k] = f_1 * pc_x[k] * qss_85[k];

        t_256[k] = f_11 * oss_73[k]
                   + f_1 * pc_y[k] * qss_85[k];

        t_257[k] = f_10 * oss_72[k]
                   + f_1 * pc_z[k] * qss_85[k];

        t_258[k] = f_1 * pc_x[k] * qss_86[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, t_263, t_264, pc_x, pc_y, pc_z, oss_73, \
                         oss_74, oss_75, qss_86, qss_87, qss_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_9 * oss_74[k]
                   + f_1 * pc_y[k] * qss_86[k];

        t_260[k] = f_8 * oss_73[k]
                   + f_1 * pc_z[k] * qss_86[k];

        t_261[k] = f_1 * pc_x[k] * qss_87[k];

        t_262[k] = f_7 * oss_75[k]
                   + f_1 * pc_y[k] * qss_87[k];

        t_263[k] = f_6 * oss_74[k]
                   + f_1 * pc_z[k] * qss_87[k];

        t_264[k] = f_1 * pc_x[k] * qss_88[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, t_270, t_271, pc_x, pc_y, pc_z, \
                         oss_75, oss_76, oss_77, qss_88, qss_89, \
                         qss_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_5 * oss_76[k]
                   + f_1 * pc_y[k] * qss_88[k];

        t_266[k] = f_4 * oss_75[k]
                   + f_1 * pc_z[k] * qss_88[k];

        t_267[k] = f_1 * pc_x[k] * qss_89[k];

        t_268[k] = f_3 * oss_77[k]
                   + f_1 * pc_y[k] * qss_89[k];

        t_269[k] = f_2 * oss_76[k]
                   + f_1 * pc_z[k] * qss_89[k];

        t_270[k] = f_1 * pc_x[k] * qss_90[k];

        t_271[k] = f_1 * pc_y[k] * qss_90[k];
    }

#pragma omp simd aligned(t_272, pc_z, oss_77, qss_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_0 * oss_77[k]
                   + f_1 * pc_z[k] * qss_90[k];
    }
}

auto
compute_prim_qsp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t oss,
                                                   const size_t qss, const size_t ncols,
                                                   const double p, const double q) -> void
{
    compute_prim_qsp_three_center_electron_repulsion_0_piece0(buffer, target, pc, oss, qss,
                                                              ncols, p, q);

    compute_prim_qsp_three_center_electron_repulsion_0_piece1(buffer, target, pc, oss, qss,
                                                              ncols, p, q);
}

}  // namespace simdt3ceri
