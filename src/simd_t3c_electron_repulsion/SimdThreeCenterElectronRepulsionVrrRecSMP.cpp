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


#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_smp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t sls,
                                                   const size_t sms, const size_t ncols,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = p / q;
    const auto f_2 = 4.0 / q;
    const auto f_3 = 0.5 / q;
    const auto f_4 = 3.5 / q;
    const auto f_5 = 1.0 / q;
    const auto f_6 = 3.0 / q;
    const auto f_7 = 1.5 / q;
    const auto f_8 = 2.5 / q;
    const auto f_9 = 2.0 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sls_0 = buffer.data(sls + 0);
    const auto *sls_1 = buffer.data(sls + 1);
    const auto *sls_2 = buffer.data(sls + 2);
    const auto *sls_3 = buffer.data(sls + 3);
    const auto *sls_4 = buffer.data(sls + 4);
    const auto *sls_5 = buffer.data(sls + 5);
    const auto *sls_6 = buffer.data(sls + 6);
    const auto *sls_7 = buffer.data(sls + 7);
    const auto *sls_8 = buffer.data(sls + 8);
    const auto *sls_9 = buffer.data(sls + 9);
    const auto *sls_10 = buffer.data(sls + 10);
    const auto *sls_11 = buffer.data(sls + 11);
    const auto *sls_12 = buffer.data(sls + 12);
    const auto *sls_13 = buffer.data(sls + 13);
    const auto *sls_14 = buffer.data(sls + 14);
    const auto *sls_15 = buffer.data(sls + 15);
    const auto *sls_16 = buffer.data(sls + 16);
    const auto *sls_17 = buffer.data(sls + 17);
    const auto *sls_18 = buffer.data(sls + 18);
    const auto *sls_19 = buffer.data(sls + 19);
    const auto *sls_20 = buffer.data(sls + 20);
    const auto *sls_21 = buffer.data(sls + 21);
    const auto *sls_22 = buffer.data(sls + 22);
    const auto *sls_23 = buffer.data(sls + 23);
    const auto *sls_24 = buffer.data(sls + 24);
    const auto *sls_25 = buffer.data(sls + 25);
    const auto *sls_26 = buffer.data(sls + 26);
    const auto *sls_27 = buffer.data(sls + 27);
    const auto *sls_28 = buffer.data(sls + 28);
    const auto *sls_29 = buffer.data(sls + 29);
    const auto *sls_30 = buffer.data(sls + 30);
    const auto *sls_31 = buffer.data(sls + 31);
    const auto *sls_32 = buffer.data(sls + 32);
    const auto *sls_33 = buffer.data(sls + 33);
    const auto *sls_34 = buffer.data(sls + 34);
    const auto *sls_35 = buffer.data(sls + 35);
    const auto *sls_36 = buffer.data(sls + 36);
    const auto *sls_37 = buffer.data(sls + 37);
    const auto *sls_38 = buffer.data(sls + 38);
    const auto *sls_39 = buffer.data(sls + 39);
    const auto *sls_40 = buffer.data(sls + 40);
    const auto *sls_41 = buffer.data(sls + 41);
    const auto *sls_42 = buffer.data(sls + 42);
    const auto *sls_43 = buffer.data(sls + 43);
    const auto *sls_44 = buffer.data(sls + 44);

    const auto *sms_0 = buffer.data(sms + 0);
    const auto *sms_1 = buffer.data(sms + 1);
    const auto *sms_2 = buffer.data(sms + 2);
    const auto *sms_3 = buffer.data(sms + 3);
    const auto *sms_4 = buffer.data(sms + 4);
    const auto *sms_5 = buffer.data(sms + 5);
    const auto *sms_6 = buffer.data(sms + 6);
    const auto *sms_7 = buffer.data(sms + 7);
    const auto *sms_8 = buffer.data(sms + 8);
    const auto *sms_9 = buffer.data(sms + 9);
    const auto *sms_10 = buffer.data(sms + 10);
    const auto *sms_11 = buffer.data(sms + 11);
    const auto *sms_12 = buffer.data(sms + 12);
    const auto *sms_13 = buffer.data(sms + 13);
    const auto *sms_14 = buffer.data(sms + 14);
    const auto *sms_15 = buffer.data(sms + 15);
    const auto *sms_16 = buffer.data(sms + 16);
    const auto *sms_17 = buffer.data(sms + 17);
    const auto *sms_18 = buffer.data(sms + 18);
    const auto *sms_19 = buffer.data(sms + 19);
    const auto *sms_20 = buffer.data(sms + 20);
    const auto *sms_21 = buffer.data(sms + 21);
    const auto *sms_22 = buffer.data(sms + 22);
    const auto *sms_23 = buffer.data(sms + 23);
    const auto *sms_24 = buffer.data(sms + 24);
    const auto *sms_25 = buffer.data(sms + 25);
    const auto *sms_26 = buffer.data(sms + 26);
    const auto *sms_27 = buffer.data(sms + 27);
    const auto *sms_28 = buffer.data(sms + 28);
    const auto *sms_29 = buffer.data(sms + 29);
    const auto *sms_30 = buffer.data(sms + 30);
    const auto *sms_31 = buffer.data(sms + 31);
    const auto *sms_32 = buffer.data(sms + 32);
    const auto *sms_33 = buffer.data(sms + 33);
    const auto *sms_34 = buffer.data(sms + 34);
    const auto *sms_35 = buffer.data(sms + 35);
    const auto *sms_36 = buffer.data(sms + 36);
    const auto *sms_37 = buffer.data(sms + 37);
    const auto *sms_38 = buffer.data(sms + 38);
    const auto *sms_39 = buffer.data(sms + 39);
    const auto *sms_40 = buffer.data(sms + 40);
    const auto *sms_41 = buffer.data(sms + 41);
    const auto *sms_42 = buffer.data(sms + 42);
    const auto *sms_43 = buffer.data(sms + 43);
    const auto *sms_44 = buffer.data(sms + 44);
    const auto *sms_45 = buffer.data(sms + 45);
    const auto *sms_46 = buffer.data(sms + 46);
    const auto *sms_47 = buffer.data(sms + 47);
    const auto *sms_48 = buffer.data(sms + 48);
    const auto *sms_49 = buffer.data(sms + 49);
    const auto *sms_50 = buffer.data(sms + 50);
    const auto *sms_51 = buffer.data(sms + 51);
    const auto *sms_52 = buffer.data(sms + 52);
    const auto *sms_53 = buffer.data(sms + 53);
    const auto *sms_54 = buffer.data(sms + 54);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pc_x, pc_y, pc_z, sls_0, sls_1, \
                         sls_2, sms_0, sms_1, sms_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sls_0[k]
                 + f_1 * pc_x[k] * sms_0[k];

        t_1[k] = f_1 * pc_y[k] * sms_0[k];

        t_2[k] = f_1 * pc_z[k] * sms_0[k];

        t_3[k] = f_2 * sls_1[k]
                 + f_1 * pc_x[k] * sms_1[k];

        t_4[k] = f_3 * sls_0[k]
                 + f_1 * pc_y[k] * sms_1[k];

        t_5[k] = f_1 * pc_z[k] * sms_1[k];

        t_6[k] = f_2 * sls_2[k]
                 + f_1 * pc_x[k] * sms_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pc_x, pc_y, pc_z, sls_0, sls_1, \
                         sls_3, sls_4, sms_2, sms_3, sms_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * pc_y[k] * sms_2[k];

        t_8[k] = f_3 * sls_0[k]
                 + f_1 * pc_z[k] * sms_2[k];

        t_9[k] = f_4 * sls_3[k]
                 + f_1 * pc_x[k] * sms_3[k];

        t_10[k] = f_5 * sls_1[k]
                  + f_1 * pc_y[k] * sms_3[k];

        t_11[k] = f_1 * pc_z[k] * sms_3[k];

        t_12[k] = f_4 * sls_4[k]
                  + f_1 * pc_x[k] * sms_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, sls_1, sls_2, \
                         sls_5, sls_6, sms_4, sms_5, sms_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * sls_2[k]
                  + f_1 * pc_y[k] * sms_4[k];

        t_14[k] = f_3 * sls_1[k]
                  + f_1 * pc_z[k] * sms_4[k];

        t_15[k] = f_4 * sls_5[k]
                  + f_1 * pc_x[k] * sms_5[k];

        t_16[k] = f_1 * pc_y[k] * sms_5[k];

        t_17[k] = f_5 * sls_2[k]
                  + f_1 * pc_z[k] * sms_5[k];

        t_18[k] = f_6 * sls_6[k]
                  + f_1 * pc_x[k] * sms_6[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, pc_x, pc_y, pc_z, sls_3, sls_4, \
                         sls_7, sls_8, sms_6, sms_7, sms_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * sls_3[k]
                  + f_1 * pc_y[k] * sms_6[k];

        t_20[k] = f_1 * pc_z[k] * sms_6[k];

        t_21[k] = f_6 * sls_7[k]
                  + f_1 * pc_x[k] * sms_7[k];

        t_22[k] = f_5 * sls_4[k]
                  + f_1 * pc_y[k] * sms_7[k];

        t_23[k] = f_3 * sls_3[k]
                  + f_1 * pc_z[k] * sms_7[k];

        t_24[k] = f_6 * sls_8[k]
                  + f_1 * pc_x[k] * sms_8[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, sls_4, sls_5, \
                         sls_9, sls_10, sms_8, sms_9, sms_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * sls_5[k]
                  + f_1 * pc_y[k] * sms_8[k];

        t_26[k] = f_5 * sls_4[k]
                  + f_1 * pc_z[k] * sms_8[k];

        t_27[k] = f_6 * sls_9[k]
                  + f_1 * pc_x[k] * sms_9[k];

        t_28[k] = f_1 * pc_y[k] * sms_9[k];

        t_29[k] = f_7 * sls_5[k]
                  + f_1 * pc_z[k] * sms_9[k];

        t_30[k] = f_8 * sls_10[k]
                  + f_1 * pc_x[k] * sms_10[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pc_x, pc_y, pc_z, sls_6, sls_7, \
                         sls_11, sls_12, sms_10, sms_11, sms_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_9 * sls_6[k]
                  + f_1 * pc_y[k] * sms_10[k];

        t_32[k] = f_1 * pc_z[k] * sms_10[k];

        t_33[k] = f_8 * sls_11[k]
                  + f_1 * pc_x[k] * sms_11[k];

        t_34[k] = f_7 * sls_7[k]
                  + f_1 * pc_y[k] * sms_11[k];

        t_35[k] = f_3 * sls_6[k]
                  + f_1 * pc_z[k] * sms_11[k];

        t_36[k] = f_8 * sls_12[k]
                  + f_1 * pc_x[k] * sms_12[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, sls_7, sls_8, sls_9, \
                         sls_13, sms_12, sms_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * sls_8[k]
                  + f_1 * pc_y[k] * sms_12[k];

        t_38[k] = f_5 * sls_7[k]
                  + f_1 * pc_z[k] * sms_12[k];

        t_39[k] = f_8 * sls_13[k]
                  + f_1 * pc_x[k] * sms_13[k];

        t_40[k] = f_3 * sls_9[k]
                  + f_1 * pc_y[k] * sms_13[k];

        t_41[k] = f_7 * sls_8[k]
                  + f_1 * pc_z[k] * sms_13[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pc_x, pc_y, pc_z, sls_9, sls_10, \
                         sls_14, sls_15, sms_14, sms_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_8 * sls_14[k]
                  + f_1 * pc_x[k] * sms_14[k];

        t_43[k] = f_1 * pc_y[k] * sms_14[k];

        t_44[k] = f_9 * sls_9[k]
                  + f_1 * pc_z[k] * sms_14[k];

        t_45[k] = f_9 * sls_15[k]
                  + f_1 * pc_x[k] * sms_15[k];

        t_46[k] = f_8 * sls_10[k]
                  + f_1 * pc_y[k] * sms_15[k];

        t_47[k] = f_1 * pc_z[k] * sms_15[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, pc_x, pc_y, pc_z, sls_10, sls_11, \
                         sls_12, sls_16, sls_17, sms_16, sms_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_9 * sls_16[k]
                  + f_1 * pc_x[k] * sms_16[k];

        t_49[k] = f_9 * sls_11[k]
                  + f_1 * pc_y[k] * sms_16[k];

        t_50[k] = f_3 * sls_10[k]
                  + f_1 * pc_z[k] * sms_16[k];

        t_51[k] = f_9 * sls_17[k]
                  + f_1 * pc_x[k] * sms_17[k];

        t_52[k] = f_7 * sls_12[k]
                  + f_1 * pc_y[k] * sms_17[k];

        t_53[k] = f_5 * sls_11[k]
                  + f_1 * pc_z[k] * sms_17[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, sls_12, sls_13, \
                         sls_14, sls_18, sls_19, sms_18, sms_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_9 * sls_18[k]
                  + f_1 * pc_x[k] * sms_18[k];

        t_55[k] = f_5 * sls_13[k]
                  + f_1 * pc_y[k] * sms_18[k];

        t_56[k] = f_7 * sls_12[k]
                  + f_1 * pc_z[k] * sms_18[k];

        t_57[k] = f_9 * sls_19[k]
                  + f_1 * pc_x[k] * sms_19[k];

        t_58[k] = f_3 * sls_14[k]
                  + f_1 * pc_y[k] * sms_19[k];

        t_59[k] = f_9 * sls_13[k]
                  + f_1 * pc_z[k] * sms_19[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, sls_14, sls_15, \
                         sls_20, sls_21, sms_20, sms_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_9 * sls_20[k]
                  + f_1 * pc_x[k] * sms_20[k];

        t_61[k] = f_1 * pc_y[k] * sms_20[k];

        t_62[k] = f_8 * sls_14[k]
                  + f_1 * pc_z[k] * sms_20[k];

        t_63[k] = f_7 * sls_21[k]
                  + f_1 * pc_x[k] * sms_21[k];

        t_64[k] = f_6 * sls_15[k]
                  + f_1 * pc_y[k] * sms_21[k];

        t_65[k] = f_1 * pc_z[k] * sms_21[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, pc_x, pc_y, pc_z, sls_15, sls_16, \
                         sls_17, sls_22, sls_23, sms_22, sms_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_7 * sls_22[k]
                  + f_1 * pc_x[k] * sms_22[k];

        t_67[k] = f_8 * sls_16[k]
                  + f_1 * pc_y[k] * sms_22[k];

        t_68[k] = f_3 * sls_15[k]
                  + f_1 * pc_z[k] * sms_22[k];

        t_69[k] = f_7 * sls_23[k]
                  + f_1 * pc_x[k] * sms_23[k];

        t_70[k] = f_9 * sls_17[k]
                  + f_1 * pc_y[k] * sms_23[k];

        t_71[k] = f_5 * sls_16[k]
                  + f_1 * pc_z[k] * sms_23[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, pc_x, pc_y, pc_z, sls_17, sls_18, \
                         sls_19, sls_24, sls_25, sms_24, sms_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_7 * sls_24[k]
                  + f_1 * pc_x[k] * sms_24[k];

        t_73[k] = f_7 * sls_18[k]
                  + f_1 * pc_y[k] * sms_24[k];

        t_74[k] = f_7 * sls_17[k]
                  + f_1 * pc_z[k] * sms_24[k];

        t_75[k] = f_7 * sls_25[k]
                  + f_1 * pc_x[k] * sms_25[k];

        t_76[k] = f_5 * sls_19[k]
                  + f_1 * pc_y[k] * sms_25[k];

        t_77[k] = f_9 * sls_18[k]
                  + f_1 * pc_z[k] * sms_25[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, pc_x, pc_y, pc_z, sls_19, sls_20, \
                         sls_26, sls_27, sms_26, sms_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_7 * sls_26[k]
                  + f_1 * pc_x[k] * sms_26[k];

        t_79[k] = f_3 * sls_20[k]
                  + f_1 * pc_y[k] * sms_26[k];

        t_80[k] = f_8 * sls_19[k]
                  + f_1 * pc_z[k] * sms_26[k];

        t_81[k] = f_7 * sls_27[k]
                  + f_1 * pc_x[k] * sms_27[k];

        t_82[k] = f_1 * pc_y[k] * sms_27[k];

        t_83[k] = f_6 * sls_20[k]
                  + f_1 * pc_z[k] * sms_27[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, sls_21, sls_22, \
                         sls_28, sls_29, sms_28, sms_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_5 * sls_28[k]
                  + f_1 * pc_x[k] * sms_28[k];

        t_85[k] = f_4 * sls_21[k]
                  + f_1 * pc_y[k] * sms_28[k];

        t_86[k] = f_1 * pc_z[k] * sms_28[k];

        t_87[k] = f_5 * sls_29[k]
                  + f_1 * pc_x[k] * sms_29[k];

        t_88[k] = f_6 * sls_22[k]
                  + f_1 * pc_y[k] * sms_29[k];

        t_89[k] = f_3 * sls_21[k]
                  + f_1 * pc_z[k] * sms_29[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, sls_22, sls_23, \
                         sls_24, sls_30, sls_31, sms_30, sms_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_5 * sls_30[k]
                  + f_1 * pc_x[k] * sms_30[k];

        t_91[k] = f_8 * sls_23[k]
                  + f_1 * pc_y[k] * sms_30[k];

        t_92[k] = f_5 * sls_22[k]
                  + f_1 * pc_z[k] * sms_30[k];

        t_93[k] = f_5 * sls_31[k]
                  + f_1 * pc_x[k] * sms_31[k];

        t_94[k] = f_9 * sls_24[k]
                  + f_1 * pc_y[k] * sms_31[k];

        t_95[k] = f_7 * sls_23[k]
                  + f_1 * pc_z[k] * sms_31[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, t_101, pc_x, pc_y, pc_z, sls_24, \
                         sls_25, sls_26, sls_32, sls_33, sms_32, \
                         sms_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_5 * sls_32[k]
                  + f_1 * pc_x[k] * sms_32[k];

        t_97[k] = f_7 * sls_25[k]
                  + f_1 * pc_y[k] * sms_32[k];

        t_98[k] = f_9 * sls_24[k]
                  + f_1 * pc_z[k] * sms_32[k];

        t_99[k] = f_5 * sls_33[k]
                  + f_1 * pc_x[k] * sms_33[k];

        t_100[k] = f_5 * sls_26[k]
                   + f_1 * pc_y[k] * sms_33[k];

        t_101[k] = f_8 * sls_25[k]
                   + f_1 * pc_z[k] * sms_33[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, t_107, pc_x, pc_y, pc_z, sls_26, \
                         sls_27, sls_34, sls_35, sms_34, sms_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_5 * sls_34[k]
                   + f_1 * pc_x[k] * sms_34[k];

        t_103[k] = f_3 * sls_27[k]
                   + f_1 * pc_y[k] * sms_34[k];

        t_104[k] = f_6 * sls_26[k]
                   + f_1 * pc_z[k] * sms_34[k];

        t_105[k] = f_5 * sls_35[k]
                   + f_1 * pc_x[k] * sms_35[k];

        t_106[k] = f_1 * pc_y[k] * sms_35[k];

        t_107[k] = f_4 * sls_27[k]
                   + f_1 * pc_z[k] * sms_35[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, t_113, pc_x, pc_y, pc_z, sls_28, \
                         sls_29, sls_36, sls_37, sms_36, sms_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_3 * sls_36[k]
                   + f_1 * pc_x[k] * sms_36[k];

        t_109[k] = f_2 * sls_28[k]
                   + f_1 * pc_y[k] * sms_36[k];

        t_110[k] = f_1 * pc_z[k] * sms_36[k];

        t_111[k] = f_3 * sls_37[k]
                   + f_1 * pc_x[k] * sms_37[k];

        t_112[k] = f_4 * sls_29[k]
                   + f_1 * pc_y[k] * sms_37[k];

        t_113[k] = f_3 * sls_28[k]
                   + f_1 * pc_z[k] * sms_37[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, pc_x, pc_y, pc_z, sls_29, \
                         sls_30, sls_31, sls_38, sls_39, sms_38, \
                         sms_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_3 * sls_38[k]
                   + f_1 * pc_x[k] * sms_38[k];

        t_115[k] = f_6 * sls_30[k]
                   + f_1 * pc_y[k] * sms_38[k];

        t_116[k] = f_5 * sls_29[k]
                   + f_1 * pc_z[k] * sms_38[k];

        t_117[k] = f_3 * sls_39[k]
                   + f_1 * pc_x[k] * sms_39[k];

        t_118[k] = f_8 * sls_31[k]
                   + f_1 * pc_y[k] * sms_39[k];

        t_119[k] = f_7 * sls_30[k]
                   + f_1 * pc_z[k] * sms_39[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, t_125, pc_x, pc_y, pc_z, sls_31, \
                         sls_32, sls_33, sls_40, sls_41, sms_40, \
                         sms_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_3 * sls_40[k]
                   + f_1 * pc_x[k] * sms_40[k];

        t_121[k] = f_9 * sls_32[k]
                   + f_1 * pc_y[k] * sms_40[k];

        t_122[k] = f_9 * sls_31[k]
                   + f_1 * pc_z[k] * sms_40[k];

        t_123[k] = f_3 * sls_41[k]
                   + f_1 * pc_x[k] * sms_41[k];

        t_124[k] = f_7 * sls_33[k]
                   + f_1 * pc_y[k] * sms_41[k];

        t_125[k] = f_8 * sls_32[k]
                   + f_1 * pc_z[k] * sms_41[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, t_131, pc_x, pc_y, pc_z, sls_33, \
                         sls_34, sls_35, sls_42, sls_43, sms_42, \
                         sms_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_3 * sls_42[k]
                   + f_1 * pc_x[k] * sms_42[k];

        t_127[k] = f_5 * sls_34[k]
                   + f_1 * pc_y[k] * sms_42[k];

        t_128[k] = f_6 * sls_33[k]
                   + f_1 * pc_z[k] * sms_42[k];

        t_129[k] = f_3 * sls_43[k]
                   + f_1 * pc_x[k] * sms_43[k];

        t_130[k] = f_3 * sls_35[k]
                   + f_1 * pc_y[k] * sms_43[k];

        t_131[k] = f_4 * sls_34[k]
                   + f_1 * pc_z[k] * sms_43[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, t_137, t_138, pc_x, pc_y, pc_z, \
                         sls_35, sls_36, sls_44, sms_44, sms_45, \
                         sms_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_3 * sls_44[k]
                   + f_1 * pc_x[k] * sms_44[k];

        t_133[k] = f_1 * pc_y[k] * sms_44[k];

        t_134[k] = f_2 * sls_35[k]
                   + f_1 * pc_z[k] * sms_44[k];

        t_135[k] = f_1 * pc_x[k] * sms_45[k];

        t_136[k] = f_0 * sls_36[k]
                   + f_1 * pc_y[k] * sms_45[k];

        t_137[k] = f_1 * pc_z[k] * sms_45[k];

        t_138[k] = f_1 * pc_x[k] * sms_46[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, t_144, pc_x, pc_y, pc_z, sls_36, \
                         sls_37, sls_38, sms_46, sms_47, sms_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_2 * sls_37[k]
                   + f_1 * pc_y[k] * sms_46[k];

        t_140[k] = f_3 * sls_36[k]
                   + f_1 * pc_z[k] * sms_46[k];

        t_141[k] = f_1 * pc_x[k] * sms_47[k];

        t_142[k] = f_4 * sls_38[k]
                   + f_1 * pc_y[k] * sms_47[k];

        t_143[k] = f_5 * sls_37[k]
                   + f_1 * pc_z[k] * sms_47[k];

        t_144[k] = f_1 * pc_x[k] * sms_48[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, t_150, pc_x, pc_y, pc_z, sls_38, \
                         sls_39, sls_40, sms_48, sms_49, sms_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_6 * sls_39[k]
                   + f_1 * pc_y[k] * sms_48[k];

        t_146[k] = f_7 * sls_38[k]
                   + f_1 * pc_z[k] * sms_48[k];

        t_147[k] = f_1 * pc_x[k] * sms_49[k];

        t_148[k] = f_8 * sls_40[k]
                   + f_1 * pc_y[k] * sms_49[k];

        t_149[k] = f_9 * sls_39[k]
                   + f_1 * pc_z[k] * sms_49[k];

        t_150[k] = f_1 * pc_x[k] * sms_50[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, t_156, pc_x, pc_y, pc_z, sls_40, \
                         sls_41, sls_42, sms_50, sms_51, sms_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_9 * sls_41[k]
                   + f_1 * pc_y[k] * sms_50[k];

        t_152[k] = f_8 * sls_40[k]
                   + f_1 * pc_z[k] * sms_50[k];

        t_153[k] = f_1 * pc_x[k] * sms_51[k];

        t_154[k] = f_7 * sls_42[k]
                   + f_1 * pc_y[k] * sms_51[k];

        t_155[k] = f_6 * sls_41[k]
                   + f_1 * pc_z[k] * sms_51[k];

        t_156[k] = f_1 * pc_x[k] * sms_52[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, \
                         sls_42, sls_43, sls_44, sms_52, sms_53, \
                         sms_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_5 * sls_43[k]
                   + f_1 * pc_y[k] * sms_52[k];

        t_158[k] = f_4 * sls_42[k]
                   + f_1 * pc_z[k] * sms_52[k];

        t_159[k] = f_1 * pc_x[k] * sms_53[k];

        t_160[k] = f_3 * sls_44[k]
                   + f_1 * pc_y[k] * sms_53[k];

        t_161[k] = f_2 * sls_43[k]
                   + f_1 * pc_z[k] * sms_53[k];

        t_162[k] = f_1 * pc_x[k] * sms_54[k];

        t_163[k] = f_1 * pc_y[k] * sms_54[k];
    }

#pragma omp simd aligned(t_164, pc_z, sls_44, sms_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_0 * sls_44[k]
                   + f_1 * pc_z[k] * sms_54[k];
    }
}

}  // namespace simdt3ceri
