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


#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_nsp_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t mss, const size_t nss,
                                                          const size_t ncols, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = p / q;
    const auto f_2 = 4.5 / q;
    const auto f_3 = 0.5 / q;
    const auto f_4 = 4.0 / q;
    const auto f_5 = 1.0 / q;
    const auto f_6 = 3.5 / q;
    const auto f_7 = 1.5 / q;
    const auto f_8 = 3.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 2.5 / q;

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

    const auto *mss_0 = buffer.data(mss + 0);
    const auto *mss_1 = buffer.data(mss + 1);
    const auto *mss_2 = buffer.data(mss + 2);
    const auto *mss_3 = buffer.data(mss + 3);
    const auto *mss_4 = buffer.data(mss + 4);
    const auto *mss_5 = buffer.data(mss + 5);
    const auto *mss_6 = buffer.data(mss + 6);
    const auto *mss_7 = buffer.data(mss + 7);
    const auto *mss_8 = buffer.data(mss + 8);
    const auto *mss_9 = buffer.data(mss + 9);
    const auto *mss_10 = buffer.data(mss + 10);
    const auto *mss_11 = buffer.data(mss + 11);
    const auto *mss_12 = buffer.data(mss + 12);
    const auto *mss_13 = buffer.data(mss + 13);
    const auto *mss_14 = buffer.data(mss + 14);
    const auto *mss_15 = buffer.data(mss + 15);
    const auto *mss_16 = buffer.data(mss + 16);
    const auto *mss_17 = buffer.data(mss + 17);
    const auto *mss_18 = buffer.data(mss + 18);
    const auto *mss_19 = buffer.data(mss + 19);
    const auto *mss_20 = buffer.data(mss + 20);
    const auto *mss_21 = buffer.data(mss + 21);
    const auto *mss_22 = buffer.data(mss + 22);
    const auto *mss_23 = buffer.data(mss + 23);
    const auto *mss_24 = buffer.data(mss + 24);
    const auto *mss_25 = buffer.data(mss + 25);
    const auto *mss_26 = buffer.data(mss + 26);
    const auto *mss_27 = buffer.data(mss + 27);
    const auto *mss_28 = buffer.data(mss + 28);
    const auto *mss_29 = buffer.data(mss + 29);
    const auto *mss_30 = buffer.data(mss + 30);
    const auto *mss_31 = buffer.data(mss + 31);
    const auto *mss_32 = buffer.data(mss + 32);
    const auto *mss_33 = buffer.data(mss + 33);
    const auto *mss_34 = buffer.data(mss + 34);
    const auto *mss_35 = buffer.data(mss + 35);
    const auto *mss_36 = buffer.data(mss + 36);
    const auto *mss_37 = buffer.data(mss + 37);
    const auto *mss_38 = buffer.data(mss + 38);
    const auto *mss_39 = buffer.data(mss + 39);
    const auto *mss_40 = buffer.data(mss + 40);
    const auto *mss_41 = buffer.data(mss + 41);
    const auto *mss_42 = buffer.data(mss + 42);
    const auto *mss_43 = buffer.data(mss + 43);
    const auto *mss_44 = buffer.data(mss + 44);
    const auto *mss_45 = buffer.data(mss + 45);
    const auto *mss_46 = buffer.data(mss + 46);
    const auto *mss_47 = buffer.data(mss + 47);
    const auto *mss_48 = buffer.data(mss + 48);
    const auto *mss_49 = buffer.data(mss + 49);
    const auto *mss_50 = buffer.data(mss + 50);
    const auto *mss_51 = buffer.data(mss + 51);
    const auto *mss_52 = buffer.data(mss + 52);
    const auto *mss_53 = buffer.data(mss + 53);

    const auto *nss_0 = buffer.data(nss + 0);
    const auto *nss_1 = buffer.data(nss + 1);
    const auto *nss_2 = buffer.data(nss + 2);
    const auto *nss_3 = buffer.data(nss + 3);
    const auto *nss_4 = buffer.data(nss + 4);
    const auto *nss_5 = buffer.data(nss + 5);
    const auto *nss_6 = buffer.data(nss + 6);
    const auto *nss_7 = buffer.data(nss + 7);
    const auto *nss_8 = buffer.data(nss + 8);
    const auto *nss_9 = buffer.data(nss + 9);
    const auto *nss_10 = buffer.data(nss + 10);
    const auto *nss_11 = buffer.data(nss + 11);
    const auto *nss_12 = buffer.data(nss + 12);
    const auto *nss_13 = buffer.data(nss + 13);
    const auto *nss_14 = buffer.data(nss + 14);
    const auto *nss_15 = buffer.data(nss + 15);
    const auto *nss_16 = buffer.data(nss + 16);
    const auto *nss_17 = buffer.data(nss + 17);
    const auto *nss_18 = buffer.data(nss + 18);
    const auto *nss_19 = buffer.data(nss + 19);
    const auto *nss_20 = buffer.data(nss + 20);
    const auto *nss_21 = buffer.data(nss + 21);
    const auto *nss_22 = buffer.data(nss + 22);
    const auto *nss_23 = buffer.data(nss + 23);
    const auto *nss_24 = buffer.data(nss + 24);
    const auto *nss_25 = buffer.data(nss + 25);
    const auto *nss_26 = buffer.data(nss + 26);
    const auto *nss_27 = buffer.data(nss + 27);
    const auto *nss_28 = buffer.data(nss + 28);
    const auto *nss_29 = buffer.data(nss + 29);
    const auto *nss_30 = buffer.data(nss + 30);
    const auto *nss_31 = buffer.data(nss + 31);
    const auto *nss_32 = buffer.data(nss + 32);
    const auto *nss_33 = buffer.data(nss + 33);
    const auto *nss_34 = buffer.data(nss + 34);
    const auto *nss_35 = buffer.data(nss + 35);
    const auto *nss_36 = buffer.data(nss + 36);
    const auto *nss_37 = buffer.data(nss + 37);
    const auto *nss_38 = buffer.data(nss + 38);
    const auto *nss_39 = buffer.data(nss + 39);
    const auto *nss_40 = buffer.data(nss + 40);
    const auto *nss_41 = buffer.data(nss + 41);
    const auto *nss_42 = buffer.data(nss + 42);
    const auto *nss_43 = buffer.data(nss + 43);
    const auto *nss_44 = buffer.data(nss + 44);
    const auto *nss_45 = buffer.data(nss + 45);
    const auto *nss_46 = buffer.data(nss + 46);
    const auto *nss_47 = buffer.data(nss + 47);
    const auto *nss_48 = buffer.data(nss + 48);
    const auto *nss_49 = buffer.data(nss + 49);
    const auto *nss_50 = buffer.data(nss + 50);
    const auto *nss_51 = buffer.data(nss + 51);
    const auto *nss_52 = buffer.data(nss + 52);
    const auto *nss_53 = buffer.data(nss + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pc_x, pc_y, pc_z, mss_0, mss_1, \
                         mss_2, nss_0, nss_1, nss_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * mss_0[k]
                 + f_1 * pc_x[k] * nss_0[k];

        t_1[k] = f_1 * pc_y[k] * nss_0[k];

        t_2[k] = f_1 * pc_z[k] * nss_0[k];

        t_3[k] = f_2 * mss_1[k]
                 + f_1 * pc_x[k] * nss_1[k];

        t_4[k] = f_3 * mss_0[k]
                 + f_1 * pc_y[k] * nss_1[k];

        t_5[k] = f_1 * pc_z[k] * nss_1[k];

        t_6[k] = f_2 * mss_2[k]
                 + f_1 * pc_x[k] * nss_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pc_x, pc_y, pc_z, mss_0, mss_1, \
                         mss_3, mss_4, nss_2, nss_3, nss_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * pc_y[k] * nss_2[k];

        t_8[k] = f_3 * mss_0[k]
                 + f_1 * pc_z[k] * nss_2[k];

        t_9[k] = f_4 * mss_3[k]
                 + f_1 * pc_x[k] * nss_3[k];

        t_10[k] = f_5 * mss_1[k]
                  + f_1 * pc_y[k] * nss_3[k];

        t_11[k] = f_1 * pc_z[k] * nss_3[k];

        t_12[k] = f_4 * mss_4[k]
                  + f_1 * pc_x[k] * nss_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, mss_1, mss_2, \
                         mss_5, mss_6, nss_4, nss_5, nss_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * mss_2[k]
                  + f_1 * pc_y[k] * nss_4[k];

        t_14[k] = f_3 * mss_1[k]
                  + f_1 * pc_z[k] * nss_4[k];

        t_15[k] = f_4 * mss_5[k]
                  + f_1 * pc_x[k] * nss_5[k];

        t_16[k] = f_1 * pc_y[k] * nss_5[k];

        t_17[k] = f_5 * mss_2[k]
                  + f_1 * pc_z[k] * nss_5[k];

        t_18[k] = f_6 * mss_6[k]
                  + f_1 * pc_x[k] * nss_6[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, pc_x, pc_y, pc_z, mss_3, mss_4, \
                         mss_7, mss_8, nss_6, nss_7, nss_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * mss_3[k]
                  + f_1 * pc_y[k] * nss_6[k];

        t_20[k] = f_1 * pc_z[k] * nss_6[k];

        t_21[k] = f_6 * mss_7[k]
                  + f_1 * pc_x[k] * nss_7[k];

        t_22[k] = f_5 * mss_4[k]
                  + f_1 * pc_y[k] * nss_7[k];

        t_23[k] = f_3 * mss_3[k]
                  + f_1 * pc_z[k] * nss_7[k];

        t_24[k] = f_6 * mss_8[k]
                  + f_1 * pc_x[k] * nss_8[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, mss_4, mss_5, \
                         mss_9, mss_10, nss_8, nss_9, nss_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * mss_5[k]
                  + f_1 * pc_y[k] * nss_8[k];

        t_26[k] = f_5 * mss_4[k]
                  + f_1 * pc_z[k] * nss_8[k];

        t_27[k] = f_6 * mss_9[k]
                  + f_1 * pc_x[k] * nss_9[k];

        t_28[k] = f_1 * pc_y[k] * nss_9[k];

        t_29[k] = f_7 * mss_5[k]
                  + f_1 * pc_z[k] * nss_9[k];

        t_30[k] = f_8 * mss_10[k]
                  + f_1 * pc_x[k] * nss_10[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pc_x, pc_y, pc_z, mss_6, mss_7, \
                         mss_11, mss_12, nss_10, nss_11, nss_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_9 * mss_6[k]
                  + f_1 * pc_y[k] * nss_10[k];

        t_32[k] = f_1 * pc_z[k] * nss_10[k];

        t_33[k] = f_8 * mss_11[k]
                  + f_1 * pc_x[k] * nss_11[k];

        t_34[k] = f_7 * mss_7[k]
                  + f_1 * pc_y[k] * nss_11[k];

        t_35[k] = f_3 * mss_6[k]
                  + f_1 * pc_z[k] * nss_11[k];

        t_36[k] = f_8 * mss_12[k]
                  + f_1 * pc_x[k] * nss_12[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, mss_7, mss_8, mss_9, \
                         mss_13, nss_12, nss_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * mss_8[k]
                  + f_1 * pc_y[k] * nss_12[k];

        t_38[k] = f_5 * mss_7[k]
                  + f_1 * pc_z[k] * nss_12[k];

        t_39[k] = f_8 * mss_13[k]
                  + f_1 * pc_x[k] * nss_13[k];

        t_40[k] = f_3 * mss_9[k]
                  + f_1 * pc_y[k] * nss_13[k];

        t_41[k] = f_7 * mss_8[k]
                  + f_1 * pc_z[k] * nss_13[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pc_x, pc_y, pc_z, mss_9, mss_10, \
                         mss_14, mss_15, nss_14, nss_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_8 * mss_14[k]
                  + f_1 * pc_x[k] * nss_14[k];

        t_43[k] = f_1 * pc_y[k] * nss_14[k];

        t_44[k] = f_9 * mss_9[k]
                  + f_1 * pc_z[k] * nss_14[k];

        t_45[k] = f_10 * mss_15[k]
                  + f_1 * pc_x[k] * nss_15[k];

        t_46[k] = f_10 * mss_10[k]
                  + f_1 * pc_y[k] * nss_15[k];

        t_47[k] = f_1 * pc_z[k] * nss_15[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, pc_x, pc_y, pc_z, mss_10, mss_11, \
                         mss_12, mss_16, mss_17, nss_16, nss_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_10 * mss_16[k]
                  + f_1 * pc_x[k] * nss_16[k];

        t_49[k] = f_9 * mss_11[k]
                  + f_1 * pc_y[k] * nss_16[k];

        t_50[k] = f_3 * mss_10[k]
                  + f_1 * pc_z[k] * nss_16[k];

        t_51[k] = f_10 * mss_17[k]
                  + f_1 * pc_x[k] * nss_17[k];

        t_52[k] = f_7 * mss_12[k]
                  + f_1 * pc_y[k] * nss_17[k];

        t_53[k] = f_5 * mss_11[k]
                  + f_1 * pc_z[k] * nss_17[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, mss_12, mss_13, \
                         mss_14, mss_18, mss_19, nss_18, nss_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_10 * mss_18[k]
                  + f_1 * pc_x[k] * nss_18[k];

        t_55[k] = f_5 * mss_13[k]
                  + f_1 * pc_y[k] * nss_18[k];

        t_56[k] = f_7 * mss_12[k]
                  + f_1 * pc_z[k] * nss_18[k];

        t_57[k] = f_10 * mss_19[k]
                  + f_1 * pc_x[k] * nss_19[k];

        t_58[k] = f_3 * mss_14[k]
                  + f_1 * pc_y[k] * nss_19[k];

        t_59[k] = f_9 * mss_13[k]
                  + f_1 * pc_z[k] * nss_19[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, mss_14, mss_15, \
                         mss_20, mss_21, nss_20, nss_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_10 * mss_20[k]
                  + f_1 * pc_x[k] * nss_20[k];

        t_61[k] = f_1 * pc_y[k] * nss_20[k];

        t_62[k] = f_10 * mss_14[k]
                  + f_1 * pc_z[k] * nss_20[k];

        t_63[k] = f_9 * mss_21[k]
                  + f_1 * pc_x[k] * nss_21[k];

        t_64[k] = f_8 * mss_15[k]
                  + f_1 * pc_y[k] * nss_21[k];

        t_65[k] = f_1 * pc_z[k] * nss_21[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, pc_x, pc_y, pc_z, mss_15, mss_16, \
                         mss_17, mss_22, mss_23, nss_22, nss_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_9 * mss_22[k]
                  + f_1 * pc_x[k] * nss_22[k];

        t_67[k] = f_10 * mss_16[k]
                  + f_1 * pc_y[k] * nss_22[k];

        t_68[k] = f_3 * mss_15[k]
                  + f_1 * pc_z[k] * nss_22[k];

        t_69[k] = f_9 * mss_23[k]
                  + f_1 * pc_x[k] * nss_23[k];

        t_70[k] = f_9 * mss_17[k]
                  + f_1 * pc_y[k] * nss_23[k];

        t_71[k] = f_5 * mss_16[k]
                  + f_1 * pc_z[k] * nss_23[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, pc_x, pc_y, pc_z, mss_17, mss_18, \
                         mss_19, mss_24, mss_25, nss_24, nss_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_9 * mss_24[k]
                  + f_1 * pc_x[k] * nss_24[k];

        t_73[k] = f_7 * mss_18[k]
                  + f_1 * pc_y[k] * nss_24[k];

        t_74[k] = f_7 * mss_17[k]
                  + f_1 * pc_z[k] * nss_24[k];

        t_75[k] = f_9 * mss_25[k]
                  + f_1 * pc_x[k] * nss_25[k];

        t_76[k] = f_5 * mss_19[k]
                  + f_1 * pc_y[k] * nss_25[k];

        t_77[k] = f_9 * mss_18[k]
                  + f_1 * pc_z[k] * nss_25[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, pc_x, pc_y, pc_z, mss_19, mss_20, \
                         mss_26, mss_27, nss_26, nss_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_9 * mss_26[k]
                  + f_1 * pc_x[k] * nss_26[k];

        t_79[k] = f_3 * mss_20[k]
                  + f_1 * pc_y[k] * nss_26[k];

        t_80[k] = f_10 * mss_19[k]
                  + f_1 * pc_z[k] * nss_26[k];

        t_81[k] = f_9 * mss_27[k]
                  + f_1 * pc_x[k] * nss_27[k];

        t_82[k] = f_1 * pc_y[k] * nss_27[k];

        t_83[k] = f_8 * mss_20[k]
                  + f_1 * pc_z[k] * nss_27[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, mss_21, mss_22, \
                         mss_28, mss_29, nss_28, nss_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_7 * mss_28[k]
                  + f_1 * pc_x[k] * nss_28[k];

        t_85[k] = f_6 * mss_21[k]
                  + f_1 * pc_y[k] * nss_28[k];

        t_86[k] = f_1 * pc_z[k] * nss_28[k];

        t_87[k] = f_7 * mss_29[k]
                  + f_1 * pc_x[k] * nss_29[k];

        t_88[k] = f_8 * mss_22[k]
                  + f_1 * pc_y[k] * nss_29[k];

        t_89[k] = f_3 * mss_21[k]
                  + f_1 * pc_z[k] * nss_29[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, mss_22, mss_23, \
                         mss_24, mss_30, mss_31, nss_30, nss_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_7 * mss_30[k]
                  + f_1 * pc_x[k] * nss_30[k];

        t_91[k] = f_10 * mss_23[k]
                  + f_1 * pc_y[k] * nss_30[k];

        t_92[k] = f_5 * mss_22[k]
                  + f_1 * pc_z[k] * nss_30[k];

        t_93[k] = f_7 * mss_31[k]
                  + f_1 * pc_x[k] * nss_31[k];

        t_94[k] = f_9 * mss_24[k]
                  + f_1 * pc_y[k] * nss_31[k];

        t_95[k] = f_7 * mss_23[k]
                  + f_1 * pc_z[k] * nss_31[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, t_101, pc_x, pc_y, pc_z, mss_24, \
                         mss_25, mss_26, mss_32, mss_33, nss_32, \
                         nss_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_7 * mss_32[k]
                  + f_1 * pc_x[k] * nss_32[k];

        t_97[k] = f_7 * mss_25[k]
                  + f_1 * pc_y[k] * nss_32[k];

        t_98[k] = f_9 * mss_24[k]
                  + f_1 * pc_z[k] * nss_32[k];

        t_99[k] = f_7 * mss_33[k]
                  + f_1 * pc_x[k] * nss_33[k];

        t_100[k] = f_5 * mss_26[k]
                   + f_1 * pc_y[k] * nss_33[k];

        t_101[k] = f_10 * mss_25[k]
                   + f_1 * pc_z[k] * nss_33[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, t_107, pc_x, pc_y, pc_z, mss_26, \
                         mss_27, mss_34, mss_35, nss_34, nss_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_7 * mss_34[k]
                   + f_1 * pc_x[k] * nss_34[k];

        t_103[k] = f_3 * mss_27[k]
                   + f_1 * pc_y[k] * nss_34[k];

        t_104[k] = f_8 * mss_26[k]
                   + f_1 * pc_z[k] * nss_34[k];

        t_105[k] = f_7 * mss_35[k]
                   + f_1 * pc_x[k] * nss_35[k];

        t_106[k] = f_1 * pc_y[k] * nss_35[k];

        t_107[k] = f_6 * mss_27[k]
                   + f_1 * pc_z[k] * nss_35[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, t_113, pc_x, pc_y, pc_z, mss_28, \
                         mss_29, mss_36, mss_37, nss_36, nss_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_5 * mss_36[k]
                   + f_1 * pc_x[k] * nss_36[k];

        t_109[k] = f_4 * mss_28[k]
                   + f_1 * pc_y[k] * nss_36[k];

        t_110[k] = f_1 * pc_z[k] * nss_36[k];

        t_111[k] = f_5 * mss_37[k]
                   + f_1 * pc_x[k] * nss_37[k];

        t_112[k] = f_6 * mss_29[k]
                   + f_1 * pc_y[k] * nss_37[k];

        t_113[k] = f_3 * mss_28[k]
                   + f_1 * pc_z[k] * nss_37[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, pc_x, pc_y, pc_z, mss_29, \
                         mss_30, mss_31, mss_38, mss_39, nss_38, \
                         nss_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_5 * mss_38[k]
                   + f_1 * pc_x[k] * nss_38[k];

        t_115[k] = f_8 * mss_30[k]
                   + f_1 * pc_y[k] * nss_38[k];

        t_116[k] = f_5 * mss_29[k]
                   + f_1 * pc_z[k] * nss_38[k];

        t_117[k] = f_5 * mss_39[k]
                   + f_1 * pc_x[k] * nss_39[k];

        t_118[k] = f_10 * mss_31[k]
                   + f_1 * pc_y[k] * nss_39[k];

        t_119[k] = f_7 * mss_30[k]
                   + f_1 * pc_z[k] * nss_39[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, t_125, pc_x, pc_y, pc_z, mss_31, \
                         mss_32, mss_33, mss_40, mss_41, nss_40, \
                         nss_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_5 * mss_40[k]
                   + f_1 * pc_x[k] * nss_40[k];

        t_121[k] = f_9 * mss_32[k]
                   + f_1 * pc_y[k] * nss_40[k];

        t_122[k] = f_9 * mss_31[k]
                   + f_1 * pc_z[k] * nss_40[k];

        t_123[k] = f_5 * mss_41[k]
                   + f_1 * pc_x[k] * nss_41[k];

        t_124[k] = f_7 * mss_33[k]
                   + f_1 * pc_y[k] * nss_41[k];

        t_125[k] = f_10 * mss_32[k]
                   + f_1 * pc_z[k] * nss_41[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, t_131, pc_x, pc_y, pc_z, mss_33, \
                         mss_34, mss_35, mss_42, mss_43, nss_42, \
                         nss_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_5 * mss_42[k]
                   + f_1 * pc_x[k] * nss_42[k];

        t_127[k] = f_5 * mss_34[k]
                   + f_1 * pc_y[k] * nss_42[k];

        t_128[k] = f_8 * mss_33[k]
                   + f_1 * pc_z[k] * nss_42[k];

        t_129[k] = f_5 * mss_43[k]
                   + f_1 * pc_x[k] * nss_43[k];

        t_130[k] = f_3 * mss_35[k]
                   + f_1 * pc_y[k] * nss_43[k];

        t_131[k] = f_6 * mss_34[k]
                   + f_1 * pc_z[k] * nss_43[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, mss_35, \
                         mss_36, mss_44, mss_45, nss_44, nss_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_5 * mss_44[k]
                   + f_1 * pc_x[k] * nss_44[k];

        t_133[k] = f_1 * pc_y[k] * nss_44[k];

        t_134[k] = f_4 * mss_35[k]
                   + f_1 * pc_z[k] * nss_44[k];

        t_135[k] = f_3 * mss_45[k]
                   + f_1 * pc_x[k] * nss_45[k];

        t_136[k] = f_2 * mss_36[k]
                   + f_1 * pc_y[k] * nss_45[k];

        t_137[k] = f_1 * pc_z[k] * nss_45[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, t_143, pc_x, pc_y, pc_z, mss_36, \
                         mss_37, mss_38, mss_46, mss_47, nss_46, \
                         nss_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_3 * mss_46[k]
                   + f_1 * pc_x[k] * nss_46[k];

        t_139[k] = f_4 * mss_37[k]
                   + f_1 * pc_y[k] * nss_46[k];

        t_140[k] = f_3 * mss_36[k]
                   + f_1 * pc_z[k] * nss_46[k];

        t_141[k] = f_3 * mss_47[k]
                   + f_1 * pc_x[k] * nss_47[k];

        t_142[k] = f_6 * mss_38[k]
                   + f_1 * pc_y[k] * nss_47[k];

        t_143[k] = f_5 * mss_37[k]
                   + f_1 * pc_z[k] * nss_47[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, t_149, pc_x, pc_y, pc_z, mss_38, \
                         mss_39, mss_40, mss_48, mss_49, nss_48, \
                         nss_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_3 * mss_48[k]
                   + f_1 * pc_x[k] * nss_48[k];

        t_145[k] = f_8 * mss_39[k]
                   + f_1 * pc_y[k] * nss_48[k];

        t_146[k] = f_7 * mss_38[k]
                   + f_1 * pc_z[k] * nss_48[k];

        t_147[k] = f_3 * mss_49[k]
                   + f_1 * pc_x[k] * nss_49[k];

        t_148[k] = f_10 * mss_40[k]
                   + f_1 * pc_y[k] * nss_49[k];

        t_149[k] = f_9 * mss_39[k]
                   + f_1 * pc_z[k] * nss_49[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, t_155, pc_x, pc_y, pc_z, mss_40, \
                         mss_41, mss_42, mss_50, mss_51, nss_50, \
                         nss_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_3 * mss_50[k]
                   + f_1 * pc_x[k] * nss_50[k];

        t_151[k] = f_9 * mss_41[k]
                   + f_1 * pc_y[k] * nss_50[k];

        t_152[k] = f_10 * mss_40[k]
                   + f_1 * pc_z[k] * nss_50[k];

        t_153[k] = f_3 * mss_51[k]
                   + f_1 * pc_x[k] * nss_51[k];

        t_154[k] = f_7 * mss_42[k]
                   + f_1 * pc_y[k] * nss_51[k];

        t_155[k] = f_8 * mss_41[k]
                   + f_1 * pc_z[k] * nss_51[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, pc_x, pc_y, pc_z, mss_42, \
                         mss_43, mss_44, mss_52, mss_53, nss_52, \
                         nss_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_3 * mss_52[k]
                   + f_1 * pc_x[k] * nss_52[k];

        t_157[k] = f_5 * mss_43[k]
                   + f_1 * pc_y[k] * nss_52[k];

        t_158[k] = f_6 * mss_42[k]
                   + f_1 * pc_z[k] * nss_52[k];

        t_159[k] = f_3 * mss_53[k]
                   + f_1 * pc_x[k] * nss_53[k];

        t_160[k] = f_3 * mss_44[k]
                   + f_1 * pc_y[k] * nss_53[k];

        t_161[k] = f_4 * mss_43[k]
                   + f_1 * pc_z[k] * nss_53[k];
    }
}

static auto
compute_prim_nsp_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t mss, const size_t nss,
                                                          const size_t ncols, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = p / q;
    const auto f_2 = 4.5 / q;
    const auto f_3 = 0.5 / q;
    const auto f_4 = 4.0 / q;
    const auto f_5 = 1.0 / q;
    const auto f_6 = 3.5 / q;
    const auto f_7 = 1.5 / q;
    const auto f_8 = 3.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 2.5 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *mss_44 = buffer.data(mss + 44);
    const auto *mss_45 = buffer.data(mss + 45);
    const auto *mss_46 = buffer.data(mss + 46);
    const auto *mss_47 = buffer.data(mss + 47);
    const auto *mss_48 = buffer.data(mss + 48);
    const auto *mss_49 = buffer.data(mss + 49);
    const auto *mss_50 = buffer.data(mss + 50);
    const auto *mss_51 = buffer.data(mss + 51);
    const auto *mss_52 = buffer.data(mss + 52);
    const auto *mss_53 = buffer.data(mss + 53);
    const auto *mss_54 = buffer.data(mss + 54);

    const auto *nss_54 = buffer.data(nss + 54);
    const auto *nss_55 = buffer.data(nss + 55);
    const auto *nss_56 = buffer.data(nss + 56);
    const auto *nss_57 = buffer.data(nss + 57);
    const auto *nss_58 = buffer.data(nss + 58);
    const auto *nss_59 = buffer.data(nss + 59);
    const auto *nss_60 = buffer.data(nss + 60);
    const auto *nss_61 = buffer.data(nss + 61);
    const auto *nss_62 = buffer.data(nss + 62);
    const auto *nss_63 = buffer.data(nss + 63);
    const auto *nss_64 = buffer.data(nss + 64);
    const auto *nss_65 = buffer.data(nss + 65);

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, t_167, t_168, pc_x, pc_y, pc_z, \
                         mss_44, mss_45, mss_54, nss_54, nss_55, \
                         nss_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_3 * mss_54[k]
                   + f_1 * pc_x[k] * nss_54[k];

        t_163[k] = f_1 * pc_y[k] * nss_54[k];

        t_164[k] = f_2 * mss_44[k]
                   + f_1 * pc_z[k] * nss_54[k];

        t_165[k] = f_1 * pc_x[k] * nss_55[k];

        t_166[k] = f_0 * mss_45[k]
                   + f_1 * pc_y[k] * nss_55[k];

        t_167[k] = f_1 * pc_z[k] * nss_55[k];

        t_168[k] = f_1 * pc_x[k] * nss_56[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, t_174, pc_x, pc_y, pc_z, mss_45, \
                         mss_46, mss_47, nss_56, nss_57, nss_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_2 * mss_46[k]
                   + f_1 * pc_y[k] * nss_56[k];

        t_170[k] = f_3 * mss_45[k]
                   + f_1 * pc_z[k] * nss_56[k];

        t_171[k] = f_1 * pc_x[k] * nss_57[k];

        t_172[k] = f_4 * mss_47[k]
                   + f_1 * pc_y[k] * nss_57[k];

        t_173[k] = f_5 * mss_46[k]
                   + f_1 * pc_z[k] * nss_57[k];

        t_174[k] = f_1 * pc_x[k] * nss_58[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, t_180, pc_x, pc_y, pc_z, mss_47, \
                         mss_48, mss_49, nss_58, nss_59, nss_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_6 * mss_48[k]
                   + f_1 * pc_y[k] * nss_58[k];

        t_176[k] = f_7 * mss_47[k]
                   + f_1 * pc_z[k] * nss_58[k];

        t_177[k] = f_1 * pc_x[k] * nss_59[k];

        t_178[k] = f_8 * mss_49[k]
                   + f_1 * pc_y[k] * nss_59[k];

        t_179[k] = f_9 * mss_48[k]
                   + f_1 * pc_z[k] * nss_59[k];

        t_180[k] = f_1 * pc_x[k] * nss_60[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, t_185, t_186, pc_x, pc_y, pc_z, mss_49, \
                         mss_50, mss_51, nss_60, nss_61, nss_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_10 * mss_50[k]
                   + f_1 * pc_y[k] * nss_60[k];

        t_182[k] = f_10 * mss_49[k]
                   + f_1 * pc_z[k] * nss_60[k];

        t_183[k] = f_1 * pc_x[k] * nss_61[k];

        t_184[k] = f_9 * mss_51[k]
                   + f_1 * pc_y[k] * nss_61[k];

        t_185[k] = f_8 * mss_50[k]
                   + f_1 * pc_z[k] * nss_61[k];

        t_186[k] = f_1 * pc_x[k] * nss_62[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, t_192, pc_x, pc_y, pc_z, mss_51, \
                         mss_52, mss_53, nss_62, nss_63, nss_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_7 * mss_52[k]
                   + f_1 * pc_y[k] * nss_62[k];

        t_188[k] = f_6 * mss_51[k]
                   + f_1 * pc_z[k] * nss_62[k];

        t_189[k] = f_1 * pc_x[k] * nss_63[k];

        t_190[k] = f_5 * mss_53[k]
                   + f_1 * pc_y[k] * nss_63[k];

        t_191[k] = f_4 * mss_52[k]
                   + f_1 * pc_z[k] * nss_63[k];

        t_192[k] = f_1 * pc_x[k] * nss_64[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, pc_x, pc_y, pc_z, mss_53, mss_54, \
                         nss_64, nss_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * mss_54[k]
                   + f_1 * pc_y[k] * nss_64[k];

        t_194[k] = f_2 * mss_53[k]
                   + f_1 * pc_z[k] * nss_64[k];

        t_195[k] = f_1 * pc_x[k] * nss_65[k];

        t_196[k] = f_1 * pc_y[k] * nss_65[k];

        t_197[k] = f_0 * mss_54[k]
                   + f_1 * pc_z[k] * nss_65[k];
    }
}

auto
compute_prim_nsp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t mss,
                                                   const size_t nss, const size_t ncols,
                                                   const double p, const double q) -> void
{
    compute_prim_nsp_three_center_electron_repulsion_0_piece0(buffer, target, pc, mss, nss,
                                                              ncols, p, q);

    compute_prim_nsp_three_center_electron_repulsion_0_piece1(buffer, target, pc, mss, nss,
                                                              ncols, p, q);
}

}  // namespace simdt3ceri
