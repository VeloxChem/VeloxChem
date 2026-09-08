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


#include "SimdElectronRepulsionVrrRecFG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_fg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t df, const size_t dg,
                                     const size_t fd0, const size_t fd1, const size_t ff,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_26 = buffer.data(df + 26);
    const auto *df_27 = buffer.data(df + 27);
    const auto *df_28 = buffer.data(df + 28);
    const auto *df_29 = buffer.data(df + 29);
    const auto *df_30 = buffer.data(df + 30);
    const auto *df_31 = buffer.data(df + 31);
    const auto *df_32 = buffer.data(df + 32);
    const auto *df_33 = buffer.data(df + 33);
    const auto *df_34 = buffer.data(df + 34);
    const auto *df_35 = buffer.data(df + 35);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_3 = buffer.data(fd0 + 3);
    const auto *fd0_4 = buffer.data(fd0 + 4);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_6 = buffer.data(fd0 + 6);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_4 = buffer.data(fd1 + 4);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_8 = buffer.data(fd1 + 8);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_41 = buffer.data(ff + 41);
    const auto *ff_42 = buffer.data(ff + 42);
    const auto *ff_43 = buffer.data(ff + 43);
    const auto *ff_44 = buffer.data(ff + 44);
    const auto *ff_45 = buffer.data(ff + 45);
    const auto *ff_46 = buffer.data(ff + 46);
    const auto *ff_47 = buffer.data(ff + 47);
    const auto *ff_48 = buffer.data(ff + 48);
    const auto *ff_49 = buffer.data(ff + 49);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_51 = buffer.data(ff + 51);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_53 = buffer.data(ff + 53);
    const auto *ff_54 = buffer.data(ff + 54);
    const auto *ff_55 = buffer.data(ff + 55);
    const auto *ff_56 = buffer.data(ff + 56);
    const auto *ff_57 = buffer.data(ff + 57);
    const auto *ff_58 = buffer.data(ff + 58);
    const auto *ff_59 = buffer.data(ff + 59);
    const auto *ff_60 = buffer.data(ff + 60);
    const auto *ff_61 = buffer.data(ff + 61);
    const auto *ff_62 = buffer.data(ff + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, \
                         ff_0, ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = pb_y[k] * ff_2[k];

        t_5[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, df_3, df_6, fd0_1, fd1_1, \
                         ff_3, ff_4, ff_5, ff_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * df_3[k]
                 + pb_x[k] * ff_5[k];

        t_7[k] = pb_z[k] * ff_3[k];

        t_8[k] = pb_y[k] * ff_4[k];

        t_9[k] = f_0 * df_6[k]
                 + pb_x[k] * ff_7[k];

        t_10[k] = f_1 * fd0_1[k]
                  - f_2 * fd1_1[k]
                  + pb_y[k] * ff_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, dg_0, fd0_2, fd1_2, \
                         ff_5, ff_6, ff_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * ff_5[k];

        t_12[k] = f_3 * fd0_2[k]
                  - f_4 * fd1_2[k]
                  + pb_y[k] * ff_6[k];

        t_13[k] = pb_y[k] * ff_7[k];

        t_14[k] = f_1 * fd0_2[k]
                  - f_2 * fd1_2[k]
                  + pb_z[k] * ff_7[k];

        t_15[k] = pa_y[k] * dg_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pb_y, pb_z, df_0, df_1, dg_1, \
                         dg_2, ff_8, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * df_0[k]
                  + pb_y[k] * ff_8[k];

        t_17[k] = pb_z[k] * ff_8[k];

        t_18[k] = f_6 * df_1[k]
                  + pa_y[k] * dg_1[k];

        t_19[k] = pb_z[k] * ff_9[k];

        t_20[k] = pa_y[k] * dg_2[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pb_x, pb_z, df_3, df_8, df_9, \
                         dg_4, dg_5, ff_10, ff_11, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_6 * df_8[k]
                  + pb_x[k] * ff_11[k];

        t_22[k] = pb_z[k] * ff_10[k];

        t_23[k] = f_6 * df_9[k]
                  + pb_x[k] * ff_12[k];

        t_24[k] = pa_y[k] * dg_4[k];

        t_25[k] = f_7 * df_3[k]
                  + pa_y[k] * dg_5[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, pb_z, df_5, df_6, \
                         dg_0, dg_6, dg_7, ff_11, ff_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * ff_11[k];

        t_27[k] = f_6 * df_5[k]
                  + pa_y[k] * dg_6[k];

        t_28[k] = f_5 * df_6[k]
                  + pb_y[k] * ff_13[k];

        t_29[k] = pa_y[k] * dg_7[k];

        t_30[k] = pa_z[k] * dg_0[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, df_0, df_2, \
                         dg_1, dg_2, dg_3, ff_14, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * ff_14[k];

        t_32[k] = f_5 * df_0[k]
                  + pb_z[k] * ff_14[k];

        t_33[k] = pa_z[k] * dg_1[k];

        t_34[k] = pb_y[k] * ff_15[k];

        t_35[k] = f_6 * df_2[k]
                  + pa_z[k] * dg_2[k];

        t_36[k] = pa_z[k] * dg_3[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, pb_y, df_12, df_13, dg_5, ff_16, \
                         ff_18, ff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_6 * df_12[k]
                  + pb_x[k] * ff_18[k];

        t_38[k] = pb_y[k] * ff_16[k];

        t_39[k] = f_6 * df_13[k]
                  + pb_x[k] * ff_19[k];

        t_40[k] = pa_z[k] * dg_5[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_z, pb_y, pb_z, df_3, df_4, df_6, dg_6, \
                         dg_7, ff_17, ff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * df_3[k]
                  + pb_z[k] * ff_17[k];

        t_42[k] = f_6 * df_4[k]
                  + pa_z[k] * dg_6[k];

        t_43[k] = pb_y[k] * ff_19[k];

        t_44[k] = f_7 * df_6[k]
                  + pa_z[k] * dg_7[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_x, pb_y, pb_z, df_7, df_14, df_16, \
                         dg_15, dg_17, ff_20, ff_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_7 * df_14[k]
                  + pa_x[k] * dg_15[k];

        t_46[k] = f_6 * df_7[k]
                  + pb_y[k] * ff_20[k];

        t_47[k] = pb_z[k] * ff_20[k];

        t_48[k] = f_6 * df_16[k]
                  + pa_x[k] * dg_17[k];

        t_49[k] = pb_z[k] * ff_21[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_x, pb_x, pb_z, df_17, df_18, df_20, dg_18, \
                         ff_22, ff_23, ff_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_6 * df_17[k]
                  + pa_x[k] * dg_18[k];

        t_51[k] = f_5 * df_18[k]
                  + pb_x[k] * ff_23[k];

        t_52[k] = pb_z[k] * ff_22[k];

        t_53[k] = f_5 * df_20[k]
                  + pb_x[k] * ff_24[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, pa_x, pb_x, pb_z, df_21, dg_19, \
                         dg_20, dg_21, dg_22, ff_23, ff_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_5 * df_21[k]
                  + pb_x[k] * ff_25[k];

        t_55[k] = pa_x[k] * dg_19[k];

        t_56[k] = pb_z[k] * ff_23[k];

        t_57[k] = pa_x[k] * dg_20[k];

        t_58[k] = pa_x[k] * dg_21[k];

        t_59[k] = pa_x[k] * dg_22[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, df_11, dg_8, \
                         dg_9, dg_11, dg_12, dg_13, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * dg_11[k];

        t_61[k] = pa_z[k] * dg_8[k];

        t_62[k] = pa_y[k] * dg_12[k];

        t_63[k] = pa_z[k] * dg_9[k];

        t_64[k] = f_5 * df_11[k]
                  + pb_y[k] * ff_26[k];

        t_65[k] = pa_y[k] * dg_13[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_x, pa_y, pa_z, pb_x, df_24, df_25, \
                         dg_10, dg_14, dg_23, ff_27, ff_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * dg_10[k];

        t_67[k] = f_5 * df_24[k]
                  + pb_x[k] * ff_27[k];

        t_68[k] = f_5 * df_25[k]
                  + pb_x[k] * ff_28[k];

        t_69[k] = pa_y[k] * dg_14[k];

        t_70[k] = pa_x[k] * dg_23[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, t_76, pa_x, pb_y, df_27, dg_24, dg_25, \
                         dg_26, dg_27, dg_28, ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pa_x[k] * dg_24[k];

        t_72[k] = pa_x[k] * dg_25[k];

        t_73[k] = pa_x[k] * dg_26[k];

        t_74[k] = pa_x[k] * dg_27[k];

        t_75[k] = f_7 * df_27[k]
                  + pa_x[k] * dg_28[k];

        t_76[k] = pb_y[k] * ff_29[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_x, pb_y, pb_z, df_10, df_30, df_31, dg_30, \
                         dg_31, ff_29, ff_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_6 * df_10[k]
                  + pb_z[k] * ff_29[k];

        t_78[k] = f_6 * df_30[k]
                  + pa_x[k] * dg_30[k];

        t_79[k] = pb_y[k] * ff_30[k];

        t_80[k] = f_6 * df_31[k]
                  + pa_x[k] * dg_31[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, pa_x, pb_x, pb_y, df_32, df_33, df_35, \
                         dg_32, ff_31, ff_32, ff_33, ff_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * df_32[k]
                  + pb_x[k] * ff_32[k];

        t_82[k] = f_5 * df_33[k]
                  + pb_x[k] * ff_33[k];

        t_83[k] = pb_y[k] * ff_31[k];

        t_84[k] = f_5 * df_35[k]
                  + pb_x[k] * ff_34[k];

        t_85[k] = pa_x[k] * dg_32[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, pa_x, pb_x, pb_y, dg_33, dg_34, dg_35, \
                         fd0_3, fd1_3, ff_34, ff_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pa_x[k] * dg_33[k];

        t_87[k] = pa_x[k] * dg_34[k];

        t_88[k] = pb_y[k] * ff_34[k];

        t_89[k] = pa_x[k] * dg_35[k];

        t_90[k] = f_1 * fd0_3[k]
                  - f_2 * fd1_3[k]
                  + pb_x[k] * ff_35[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pb_x, pb_y, pb_z, df_14, fd0_4, fd1_4, ff_35, \
                         ff_36, ff_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_0 * df_14[k]
                  + pb_y[k] * ff_35[k];

        t_92[k] = pb_z[k] * ff_35[k];

        t_93[k] = f_3 * fd0_4[k]
                  - f_4 * fd1_4[k]
                  + pb_x[k] * ff_37[k];

        t_94[k] = pb_z[k] * ff_36[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pb_x, fd0_5, fd1_5, ff_38, ff_39, \
                         ff_40, ff_41, ff_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_5[k]
                  + pb_x[k] * ff_38[k];

        t_96[k] = pb_x[k] * ff_39[k];

        t_97[k] = pb_x[k] * ff_40[k];

        t_98[k] = pb_x[k] * ff_41[k];

        t_99[k] = pb_x[k] * ff_42[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pb_y, pb_z, df_18, df_21, fd0_4, \
                         fd0_5, fd1_4, fd1_5, ff_39, ff_40, ff_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_0 * df_18[k]
                   + f_1 * fd0_4[k]
                   - f_2 * fd1_4[k]
                   + pb_y[k] * ff_39[k];

        t_101[k] = pb_z[k] * ff_39[k];

        t_102[k] = f_3 * fd0_4[k]
                   - f_4 * fd1_4[k]
                   + pb_z[k] * ff_40[k];

        t_103[k] = f_0 * df_21[k]
                   + pb_y[k] * ff_42[k];

        t_104[k] = f_1 * fd0_5[k]
                   - f_2 * fd1_5[k]
                   + pb_z[k] * ff_42[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pa_z, pb_y, pb_z, df_14, df_22, \
                         dg_15, dg_16, dg_17, ff_43, ff_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pa_z[k] * dg_15[k];

        t_106[k] = pa_z[k] * dg_16[k];

        t_107[k] = f_5 * df_14[k]
                   + pb_z[k] * ff_43[k];

        t_108[k] = pa_z[k] * dg_17[k];

        t_109[k] = f_6 * df_22[k]
                   + pb_y[k] * ff_44[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, t_115, pa_z, pb_x, df_15, dg_18, \
                         dg_19, ff_45, ff_46, ff_47, ff_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_6 * df_15[k]
                   + pa_z[k] * dg_18[k];

        t_111[k] = pb_x[k] * ff_45[k];

        t_112[k] = pb_x[k] * ff_46[k];

        t_113[k] = pb_x[k] * ff_47[k];

        t_114[k] = pb_x[k] * ff_48[k];

        t_115[k] = pa_z[k] * dg_19[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_z, pb_y, pb_z, df_18, df_19, df_21, \
                         df_26, dg_20, dg_22, ff_45, ff_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_5 * df_18[k]
                   + pb_z[k] * ff_45[k];

        t_117[k] = f_6 * df_19[k]
                   + pa_z[k] * dg_20[k];

        t_118[k] = f_6 * df_26[k]
                   + pb_y[k] * ff_48[k];

        t_119[k] = f_7 * df_21[k]
                   + pa_z[k] * dg_22[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pa_y, pb_y, df_27, df_28, df_29, \
                         dg_28, dg_29, dg_30, ff_49, ff_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_y[k] * dg_28[k];

        t_121[k] = f_5 * df_27[k]
                   + pb_y[k] * ff_49[k];

        t_122[k] = pa_y[k] * dg_29[k];

        t_123[k] = f_6 * df_28[k]
                   + pa_y[k] * dg_30[k];

        t_124[k] = f_5 * df_29[k]
                   + pb_y[k] * ff_50[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, t_130, pa_y, pb_x, df_32, dg_31, \
                         dg_32, ff_51, ff_52, ff_53, ff_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pa_y[k] * dg_31[k];

        t_126[k] = pb_x[k] * ff_51[k];

        t_127[k] = pb_x[k] * ff_52[k];

        t_128[k] = pb_x[k] * ff_53[k];

        t_129[k] = pb_x[k] * ff_54[k];

        t_130[k] = f_7 * df_32[k]
                   + pa_y[k] * dg_32[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_y, pb_y, pb_z, df_23, df_34, df_35, \
                         dg_34, dg_35, ff_51, ff_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_6 * df_23[k]
                   + pb_z[k] * ff_51[k];

        t_132[k] = f_6 * df_34[k]
                   + pa_y[k] * dg_34[k];

        t_133[k] = f_5 * df_35[k]
                   + pb_y[k] * ff_54[k];

        t_134[k] = pa_y[k] * dg_35[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, pb_x, pb_y, pb_z, df_27, fd0_6, \
                         fd0_7, fd1_6, fd1_7, ff_55, ff_56, ff_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_1 * fd0_6[k]
                   - f_2 * fd1_6[k]
                   + pb_x[k] * ff_55[k];

        t_136[k] = pb_y[k] * ff_55[k];

        t_137[k] = f_0 * df_27[k]
                   + pb_z[k] * ff_55[k];

        t_138[k] = f_3 * fd0_7[k]
                   - f_4 * fd1_7[k]
                   + pb_x[k] * ff_57[k];

        t_139[k] = pb_y[k] * ff_56[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pb_x, fd0_8, fd1_8, ff_58, ff_59, \
                         ff_60, ff_61, ff_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_3 * fd0_8[k]
                   - f_4 * fd1_8[k]
                   + pb_x[k] * ff_58[k];

        t_141[k] = pb_x[k] * ff_59[k];

        t_142[k] = pb_x[k] * ff_60[k];

        t_143[k] = pb_x[k] * ff_61[k];

        t_144[k] = pb_x[k] * ff_62[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pb_y, pb_z, df_32, df_35, fd0_7, \
                         fd0_8, fd1_7, fd1_8, ff_59, ff_61, ff_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_1 * fd0_7[k]
                   - f_2 * fd1_7[k]
                   + pb_y[k] * ff_59[k];

        t_146[k] = f_0 * df_32[k]
                   + pb_z[k] * ff_59[k];

        t_147[k] = f_3 * fd0_8[k]
                   - f_4 * fd1_8[k]
                   + pb_y[k] * ff_61[k];

        t_148[k] = pb_y[k] * ff_62[k];

        t_149[k] = f_0 * df_35[k]
                   + f_1 * fd0_8[k]
                   - f_2 * fd1_8[k]
                   + pb_z[k] * ff_62[k];
    }
}

auto
compute_prim_fg_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t df, const size_t dg,
                                     const size_t fd0, const size_t fd1, const size_t ff,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_26 = buffer.data(df + 26);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_41 = buffer.data(dg + 41);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_3 = buffer.data(fd0 + 3);
    const auto *fd0_4 = buffer.data(fd0 + 4);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_6 = buffer.data(fd0 + 6);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_4 = buffer.data(fd1 + 4);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_8 = buffer.data(fd1 + 8);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, df_3, df_6, fd0_1, fd0_2, fd1_1, \
                         fd1_2, ff_3, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * df_3[k]
                 + pb_x[k] * ff_3[k];

        t_6[k] = f_0 * df_6[k]
                 + pb_x[k] * ff_5[k];

        t_7[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_8[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_4[k];

        t_9[k] = pb_y[k] * ff_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, df_0, df_1, dg_0, dg_3, \
                         fd0_2, fd1_2, ff_5, ff_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * fd0_2[k]
                  - f_2 * fd1_2[k]
                  + pb_z[k] * ff_5[k];

        t_11[k] = pa_y[k] * dg_0[k];

        t_12[k] = f_5 * df_0[k]
                  + pb_y[k] * ff_6[k];

        t_13[k] = f_6 * df_1[k]
                  + pa_y[k] * dg_3[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pb_x, df_3, df_5, df_8, dg_4, dg_5, \
                         dg_6, ff_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_y[k] * dg_4[k];

        t_15[k] = f_6 * df_8[k]
                  + pb_x[k] * ff_7[k];

        t_16[k] = f_7 * df_3[k]
                  + pa_y[k] * dg_5[k];

        t_17[k] = f_6 * df_5[k]
                  + pa_y[k] * dg_6[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_y, pa_z, pb_y, pb_z, df_0, df_6, \
                         dg_0, dg_3, dg_8, ff_8, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * df_6[k]
                  + pb_y[k] * ff_8[k];

        t_19[k] = pa_y[k] * dg_8[k];

        t_20[k] = pa_z[k] * dg_0[k];

        t_21[k] = f_5 * df_0[k]
                  + pb_z[k] * ff_9[k];

        t_22[k] = pa_z[k] * dg_3[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pb_x, pb_z, df_2, df_3, df_10, dg_4, \
                         dg_5, ff_10, ff_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * df_2[k]
                  + pa_z[k] * dg_4[k];

        t_24[k] = f_6 * df_10[k]
                  + pb_x[k] * ff_11[k];

        t_25[k] = pa_z[k] * dg_5[k];

        t_26[k] = f_5 * df_3[k]
                  + pb_z[k] * ff_10[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pa_z, pb_y, df_4, df_6, df_7, df_11, \
                         dg_6, dg_8, dg_16, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_6 * df_4[k]
                  + pa_z[k] * dg_6[k];

        t_28[k] = f_7 * df_6[k]
                  + pa_z[k] * dg_8[k];

        t_29[k] = f_7 * df_11[k]
                  + pa_x[k] * dg_16[k];

        t_30[k] = f_6 * df_7[k]
                  + pb_y[k] * ff_12[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_x, pb_x, df_13, df_14, df_15, dg_17, \
                         dg_18, dg_21, dg_23, ff_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_6 * df_13[k]
                  + pa_x[k] * dg_17[k];

        t_32[k] = f_6 * df_14[k]
                  + pa_x[k] * dg_18[k];

        t_33[k] = f_5 * df_15[k]
                  + pb_x[k] * ff_13[k];

        t_34[k] = pa_x[k] * dg_21[k];

        t_35[k] = pa_x[k] * dg_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, pa_x, pa_y, pa_z, dg_10, dg_13, \
                         dg_14, dg_24, dg_25, dg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_x[k] * dg_24[k];

        t_37[k] = pa_x[k] * dg_25[k];

        t_38[k] = pa_y[k] * dg_13[k];

        t_39[k] = pa_z[k] * dg_10[k];

        t_40[k] = pa_y[k] * dg_14[k];

        t_41[k] = pa_x[k] * dg_27[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_x, pb_z, df_9, df_20, df_22, dg_28, \
                         dg_29, dg_31, dg_33, ff_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_x[k] * dg_28[k];

        t_43[k] = pa_x[k] * dg_29[k];

        t_44[k] = f_7 * df_20[k]
                  + pa_x[k] * dg_31[k];

        t_45[k] = f_6 * df_9[k]
                  + pb_z[k] * ff_14[k];

        t_46[k] = f_6 * df_22[k]
                  + pa_x[k] * dg_33[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, pa_x, pb_x, df_23, df_26, dg_34, \
                         dg_37, dg_38, dg_39, dg_41, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * df_23[k]
                  + pa_x[k] * dg_34[k];

        t_48[k] = f_5 * df_26[k]
                  + pb_x[k] * ff_15[k];

        t_49[k] = pa_x[k] * dg_37[k];

        t_50[k] = pa_x[k] * dg_38[k];

        t_51[k] = pa_x[k] * dg_39[k];

        t_52[k] = pa_x[k] * dg_41[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_x, pb_y, df_11, fd0_3, fd0_4, fd0_5, \
                         fd1_3, fd1_4, fd1_5, ff_16, ff_17, ff_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_1 * fd0_3[k]
                  - f_2 * fd1_3[k]
                  + pb_x[k] * ff_16[k];

        t_54[k] = f_0 * df_11[k]
                  + pb_y[k] * ff_16[k];

        t_55[k] = f_3 * fd0_4[k]
                  - f_4 * fd1_4[k]
                  + pb_x[k] * ff_17[k];

        t_56[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_5[k]
                  + pb_x[k] * ff_18[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, t_62, pb_x, pb_y, pb_z, df_15, df_17, \
                         fd0_4, fd1_4, ff_19, ff_20, ff_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_x[k] * ff_19[k];

        t_58[k] = pb_x[k] * ff_21[k];

        t_59[k] = f_0 * df_15[k]
                  + f_1 * fd0_4[k]
                  - f_2 * fd1_4[k]
                  + pb_y[k] * ff_19[k];

        t_60[k] = pb_z[k] * ff_19[k];

        t_61[k] = f_3 * fd0_4[k]
                  - f_4 * fd1_4[k]
                  + pb_z[k] * ff_20[k];

        t_62[k] = f_0 * df_17[k]
                  + pb_y[k] * ff_21[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_z, pb_z, df_11, df_12, dg_16, dg_17, \
                         dg_18, fd0_5, fd1_5, ff_21, ff_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_1 * fd0_5[k]
                  - f_2 * fd1_5[k]
                  + pb_z[k] * ff_21[k];

        t_64[k] = pa_z[k] * dg_16[k];

        t_65[k] = f_5 * df_11[k]
                  + pb_z[k] * ff_22[k];

        t_66[k] = pa_z[k] * dg_17[k];

        t_67[k] = f_6 * df_12[k]
                  + pa_z[k] * dg_18[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_z, pb_y, pb_z, df_15, df_16, df_19, dg_21, \
                         dg_23, ff_23, ff_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pa_z[k] * dg_21[k];

        t_69[k] = f_5 * df_15[k]
                  + pb_z[k] * ff_23[k];

        t_70[k] = f_6 * df_16[k]
                  + pa_z[k] * dg_23[k];

        t_71[k] = f_6 * df_19[k]
                  + pb_y[k] * ff_24[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pa_y, pa_z, df_17, df_21, dg_25, dg_31, \
                         dg_32, dg_33, dg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_7 * df_17[k]
                  + pa_z[k] * dg_25[k];

        t_73[k] = pa_y[k] * dg_31[k];

        t_74[k] = pa_y[k] * dg_32[k];

        t_75[k] = f_6 * df_21[k]
                  + pa_y[k] * dg_33[k];

        t_76[k] = pa_y[k] * dg_34[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pb_y, pb_z, df_18, df_24, df_25, df_26, \
                         dg_37, dg_39, ff_25, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_7 * df_24[k]
                  + pa_y[k] * dg_37[k];

        t_78[k] = f_6 * df_18[k]
                  + pb_z[k] * ff_25[k];

        t_79[k] = f_6 * df_25[k]
                  + pa_y[k] * dg_39[k];

        t_80[k] = f_5 * df_26[k]
                  + pb_y[k] * ff_26[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_y, pb_x, pb_z, df_20, dg_41, fd0_6, fd0_7, \
                         fd1_6, fd1_7, ff_27, ff_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pa_y[k] * dg_41[k];

        t_82[k] = f_1 * fd0_6[k]
                  - f_2 * fd1_6[k]
                  + pb_x[k] * ff_27[k];

        t_83[k] = f_0 * df_20[k]
                  + pb_z[k] * ff_27[k];

        t_84[k] = f_3 * fd0_7[k]
                  - f_4 * fd1_7[k]
                  + pb_x[k] * ff_28[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pb_x, pb_y, pb_z, df_24, fd0_7, fd0_8, \
                         fd1_7, fd1_8, ff_29, ff_30, ff_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_x[k] * ff_29[k];

        t_86[k] = pb_x[k] * ff_30[k];

        t_87[k] = pb_x[k] * ff_32[k];

        t_88[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_7[k]
                  + pb_y[k] * ff_30[k];

        t_89[k] = f_0 * df_24[k]
                  + pb_z[k] * ff_30[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_y, pb_z, df_26, fd0_8, fd1_8, ff_31, \
                         ff_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_y[k] * ff_31[k];

        t_91[k] = pb_y[k] * ff_32[k];

        t_92[k] = f_0 * df_26[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_8[k]
                  + pb_z[k] * ff_32[k];
    }
}

auto
compute_prim_fg_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t df, const size_t fd0, const size_t fd1,
                                     const size_t ff, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_17 = buffer.data(df + 17);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_3 = buffer.data(fd0 + 3);
    const auto *fd0_4 = buffer.data(fd0 + 4);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_6 = buffer.data(fd0 + 6);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_4 = buffer.data(fd1 + 4);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_8 = buffer.data(fd1 + 8);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, fd0_1, fd0_2, fd1_1, fd1_2, ff_3, \
                         ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_6[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_4[k];

        t_7[k] = pb_y[k] * ff_5[k];

        t_8[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_x, fd0_3, fd0_4, fd0_5, fd1_3, fd1_4, \
                         fd1_5, ff_6, ff_7, ff_8, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_1 * fd0_3[k]
                 - f_2 * fd1_3[k]
                 + pb_x[k] * ff_6[k];

        t_10[k] = f_3 * fd0_4[k]
                  - f_4 * fd1_4[k]
                  + pb_x[k] * ff_7[k];

        t_11[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_5[k]
                  + pb_x[k] * ff_8[k];

        t_12[k] = pb_x[k] * ff_9[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pb_x, pb_y, pb_z, df_9, fd0_4, fd0_5, \
                         fd1_4, fd1_5, ff_9, ff_10, ff_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_x[k] * ff_11[k];

        t_14[k] = f_0 * df_9[k]
                  + f_1 * fd0_4[k]
                  - f_2 * fd1_4[k]
                  + pb_y[k] * ff_9[k];

        t_15[k] = pb_z[k] * ff_9[k];

        t_16[k] = f_3 * fd0_4[k]
                  - f_4 * fd1_4[k]
                  + pb_z[k] * ff_10[k];

        t_17[k] = f_1 * fd0_5[k]
                  - f_2 * fd1_5[k]
                  + pb_z[k] * ff_11[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_x, fd0_6, fd0_7, fd0_8, fd1_6, fd1_7, \
                         fd1_8, ff_12, ff_13, ff_14, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_1 * fd0_6[k]
                  - f_2 * fd1_6[k]
                  + pb_x[k] * ff_12[k];

        t_19[k] = f_3 * fd0_7[k]
                  - f_4 * fd1_7[k]
                  + pb_x[k] * ff_13[k];

        t_20[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_x[k] * ff_14[k];

        t_21[k] = pb_x[k] * ff_15[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pb_x, pb_y, pb_z, df_17, fd0_7, fd0_8, \
                         fd1_7, fd1_8, ff_15, ff_16, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_x[k] * ff_17[k];

        t_23[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_7[k]
                  + pb_y[k] * ff_15[k];

        t_24[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_y[k] * ff_16[k];

        t_25[k] = pb_y[k] * ff_17[k];

        t_26[k] = f_0 * df_17[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_8[k]
                  + pb_z[k] * ff_17[k];
    }
}

auto
compute_prim_fg_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t df, const size_t dg,
                                     const size_t fd0, const size_t fd1, const size_t ff,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_17 = buffer.data(df + 17);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_29 = buffer.data(dg + 29);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_3 = buffer.data(fd0 + 3);
    const auto *fd0_4 = buffer.data(fd0 + 4);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_6 = buffer.data(fd0 + 6);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_4 = buffer.data(fd1 + 4);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_8 = buffer.data(fd1 + 8);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, dg_0, fd0_1, fd0_2, fd1_1, \
                         fd1_2, ff_3, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_6[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_4[k];

        t_7[k] = pb_y[k] * ff_5[k];

        t_8[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_5[k];

        t_9[k] = pa_y[k] * dg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pa_x, pa_y, pa_z, df_3, df_5, \
                         df_6, dg_0, dg_5, dg_8, dg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * df_3[k]
                  + pa_y[k] * dg_5[k];

        t_11[k] = pa_y[k] * dg_8[k];

        t_12[k] = pa_z[k] * dg_0[k];

        t_13[k] = pa_z[k] * dg_5[k];

        t_14[k] = f_5 * df_5[k]
                  + pa_z[k] * dg_8[k];

        t_15[k] = f_5 * df_6[k]
                  + pa_x[k] * dg_10[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_x, df_12, dg_15, dg_21, dg_29, \
                         fd0_3, fd1_3, ff_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_x[k] * dg_15[k];

        t_17[k] = f_5 * df_12[k]
                  + pa_x[k] * dg_21[k];

        t_18[k] = pa_x[k] * dg_29[k];

        t_19[k] = f_1 * fd0_3[k]
                  - f_2 * fd1_3[k]
                  + pb_x[k] * ff_6[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pb_x, pb_y, df_9, fd0_4, fd0_5, fd1_4, \
                         fd1_5, ff_7, ff_8, ff_9, ff_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_3 * fd0_4[k]
                  - f_4 * fd1_4[k]
                  + pb_x[k] * ff_7[k];

        t_21[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_5[k]
                  + pb_x[k] * ff_8[k];

        t_22[k] = pb_x[k] * ff_9[k];

        t_23[k] = pb_x[k] * ff_11[k];

        t_24[k] = f_0 * df_9[k]
                  + f_1 * fd0_4[k]
                  - f_2 * fd1_4[k]
                  + pb_y[k] * ff_9[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_z, pb_z, dg_10, dg_15, fd0_4, fd0_5, \
                         fd1_4, fd1_5, ff_9, ff_10, ff_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pb_z[k] * ff_9[k];

        t_26[k] = f_3 * fd0_4[k]
                  - f_4 * fd1_4[k]
                  + pb_z[k] * ff_10[k];

        t_27[k] = f_1 * fd0_5[k]
                  - f_2 * fd1_5[k]
                  + pb_z[k] * ff_11[k];

        t_28[k] = pa_z[k] * dg_10[k];

        t_29[k] = pa_z[k] * dg_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pa_z, pb_x, df_11, df_15, dg_18, dg_26, \
                         dg_29, fd0_6, fd1_6, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * df_11[k]
                  + pa_z[k] * dg_18[k];

        t_31[k] = f_5 * df_15[k]
                  + pa_y[k] * dg_26[k];

        t_32[k] = pa_y[k] * dg_29[k];

        t_33[k] = f_1 * fd0_6[k]
                  - f_2 * fd1_6[k]
                  + pb_x[k] * ff_12[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pb_x, pb_y, fd0_7, fd0_8, fd1_7, fd1_8, \
                         ff_13, ff_14, ff_15, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * fd0_7[k]
                  - f_4 * fd1_7[k]
                  + pb_x[k] * ff_13[k];

        t_35[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_x[k] * ff_14[k];

        t_36[k] = pb_x[k] * ff_15[k];

        t_37[k] = pb_x[k] * ff_17[k];

        t_38[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_7[k]
                  + pb_y[k] * ff_15[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, pb_z, df_17, fd0_8, fd1_8, ff_16, \
                         ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_y[k] * ff_16[k];

        t_40[k] = pb_y[k] * ff_17[k];

        t_41[k] = f_0 * df_17[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_8[k]
                  + pb_z[k] * ff_17[k];
    }
}

auto
compute_prim_fg_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t df, const size_t dg,
                                     const size_t fd0, const size_t fd1, const size_t ff,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_26 = buffer.data(df + 26);
    const auto *df_28 = buffer.data(df + 28);
    const auto *df_29 = buffer.data(df + 29);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_9 = buffer.data(fd0 + 9);
    const auto *fd0_12 = buffer.data(fd0 + 12);
    const auto *fd0_13 = buffer.data(fd0 + 13);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_13 = buffer.data(fd1 + 13);
    const auto *fd1_14 = buffer.data(fd1 + 14);
    const auto *fd1_15 = buffer.data(fd1 + 15);
    const auto *fd1_21 = buffer.data(fd1 + 21);
    const auto *fd1_22 = buffer.data(fd1 + 22);
    const auto *fd1_23 = buffer.data(fd1 + 23);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, df_0, df_3, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2, ff_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_2[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];

        t_3[k] = f_0 * df_3[k]
                 + pb_x[k] * ff_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, df_5, fd0_1, fd0_2, fd1_1, \
                         fd1_2, ff_3, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_0 * df_5[k]
                 + pb_x[k] * ff_5[k];

        t_5[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_6[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_4[k];

        t_7[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_5[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_y, pb_x, pb_y, df_0, df_1, df_7, dg_0, dg_1, \
                         ff_6, ff_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = pa_y[k] * dg_0[k];

        t_9[k] = f_5 * df_0[k]
                 + pb_y[k] * ff_6[k];

        t_10[k] = f_6 * df_1[k]
                  + pa_y[k] * dg_1[k];

        t_11[k] = f_6 * df_7[k]
                  + pb_x[k] * ff_7[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_y, pa_z, pb_z, df_0, df_2, df_3, dg_0, \
                         dg_2, dg_3, ff_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * df_3[k]
                  + pa_y[k] * dg_3[k];

        t_13[k] = pa_z[k] * dg_0[k];

        t_14[k] = f_5 * df_0[k]
                  + pb_z[k] * ff_8[k];

        t_15[k] = f_6 * df_2[k]
                  + pa_z[k] * dg_2[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_z, pb_x, df_4, df_5, df_10, df_11, \
                         dg_4, dg_5, dg_6, ff_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_6 * df_10[k]
                  + pb_x[k] * ff_10[k];

        t_17[k] = f_6 * df_4[k]
                  + pa_z[k] * dg_4[k];

        t_18[k] = f_7 * df_5[k]
                  + pa_z[k] * dg_5[k];

        t_19[k] = f_7 * df_11[k]
                  + pa_x[k] * dg_6[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_y, df_6, df_13, df_14, dg_7, \
                         dg_9, ff_11, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_6 * df_6[k]
                  + pb_y[k] * ff_11[k];

        t_21[k] = f_6 * df_13[k]
                  + pa_x[k] * dg_7[k];

        t_22[k] = f_5 * df_14[k]
                  + pb_x[k] * ff_12[k];

        t_23[k] = pa_x[k] * dg_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_x, pb_x, pb_z, df_8, df_22, df_25, df_29, \
                         dg_12, dg_14, ff_13, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_7 * df_22[k]
                  + pa_x[k] * dg_12[k];

        t_25[k] = f_6 * df_8[k]
                  + pb_z[k] * ff_13[k];

        t_26[k] = f_6 * df_25[k]
                  + pa_x[k] * dg_14[k];

        t_27[k] = f_5 * df_29[k]
                  + pb_x[k] * ff_15[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_x, pb_y, df_11, dg_17, fd0_7, fd0_8, \
                         fd1_13, fd1_14, ff_16, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_x[k] * dg_17[k];

        t_29[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_13[k]
                  + pb_x[k] * ff_16[k];

        t_30[k] = f_0 * df_11[k]
                  + pb_y[k] * ff_16[k];

        t_31[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_14[k]
                  + pb_x[k] * ff_17[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, pb_z, df_14, fd0_8, fd0_9, fd1_14, \
                         fd1_15, ff_18, ff_19, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_15[k]
                  + pb_x[k] * ff_18[k];

        t_33[k] = f_0 * df_14[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_14[k]
                  + pb_y[k] * ff_19[k];

        t_34[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_14[k]
                  + pb_z[k] * ff_20[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_z, pb_y, pb_z, df_12, df_17, dg_8, dg_9, \
                         fd0_9, fd1_15, ff_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * df_17[k]
                  + pb_y[k] * ff_22[k];

        t_36[k] = f_1 * fd0_9[k]
                  - f_2 * fd1_15[k]
                  + pb_z[k] * ff_22[k];

        t_37[k] = f_6 * df_12[k]
                  + pa_z[k] * dg_8[k];

        t_38[k] = pa_z[k] * dg_9[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pb_y, pb_z, df_14, df_15, df_17, df_21, \
                         dg_10, dg_11, ff_23, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_5 * df_14[k]
                  + pb_z[k] * ff_23[k];

        t_40[k] = f_6 * df_15[k]
                  + pa_z[k] * dg_10[k];

        t_41[k] = f_6 * df_21[k]
                  + pb_y[k] * ff_26[k];

        t_42[k] = f_7 * df_17[k]
                  + pa_z[k] * dg_11[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_y, pb_z, df_18, df_23, df_26, df_28, \
                         dg_13, dg_15, dg_16, ff_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_6 * df_23[k]
                  + pa_y[k] * dg_13[k];

        t_44[k] = f_7 * df_26[k]
                  + pa_y[k] * dg_15[k];

        t_45[k] = f_6 * df_18[k]
                  + pb_z[k] * ff_27[k];

        t_46[k] = f_6 * df_28[k]
                  + pa_y[k] * dg_16[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pb_x, pb_y, pb_z, df_22, df_29, dg_17, \
                         fd0_12, fd1_21, ff_30, ff_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_5 * df_29[k]
                  + pb_y[k] * ff_30[k];

        t_48[k] = pa_y[k] * dg_17[k];

        t_49[k] = f_1 * fd0_12[k]
                  - f_2 * fd1_21[k]
                  + pb_x[k] * ff_31[k];

        t_50[k] = f_0 * df_22[k]
                  + pb_z[k] * ff_31[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_x, pb_y, pb_z, df_26, fd0_13, fd0_14, \
                         fd1_22, fd1_23, ff_33, ff_34, ff_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * fd0_13[k]
                  - f_4 * fd1_22[k]
                  + pb_x[k] * ff_33[k];

        t_52[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_23[k]
                  + pb_x[k] * ff_34[k];

        t_53[k] = f_1 * fd0_13[k]
                  - f_2 * fd1_22[k]
                  + pb_y[k] * ff_35[k];

        t_54[k] = f_0 * df_26[k]
                  + pb_z[k] * ff_35[k];
    }

#pragma omp simd aligned(t_55, t_56, pb_y, pb_z, df_29, fd0_14, fd1_23, ff_37, \
                         ff_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_23[k]
                  + pb_y[k] * ff_37[k];

        t_56[k] = f_0 * df_29[k]
                  + f_1 * fd0_14[k]
                  - f_2 * fd1_23[k]
                  + pb_z[k] * ff_38[k];
    }
}

auto
compute_prim_fg_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t df, const size_t dg,
                                     const size_t fd0, const size_t fd1, const size_t ff,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / p;
    const auto f_7 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_41 = buffer.data(dg + 41);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_6 = buffer.data(fd0 + 6);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_9 = buffer.data(fd0 + 9);
    const auto *fd0_10 = buffer.data(fd0 + 10);
    const auto *fd0_11 = buffer.data(fd0 + 11);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_9 = buffer.data(fd1 + 9);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_13 = buffer.data(fd1 + 13);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_41 = buffer.data(ff + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, fd0_1, fd0_2, fd1_1, fd1_2, \
                         ff_3, ff_5, ff_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_6[k] = pb_z[k] * ff_3[k];

        t_7[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_5[k];

        t_8[k] = pb_y[k] * ff_6[k];

        t_9[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, df_1, df_3, df_5, dg_0, dg_3, \
                         dg_4, dg_5, dg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * dg_0[k];

        t_11[k] = f_5 * df_1[k]
                  + pa_y[k] * dg_3[k];

        t_12[k] = pa_y[k] * dg_4[k];

        t_13[k] = f_6 * df_3[k]
                  + pa_y[k] * dg_5[k];

        t_14[k] = f_5 * df_5[k]
                  + pa_y[k] * dg_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pa_z, pb_y, pb_z, df_0, df_6, \
                         dg_0, dg_3, dg_8, ff_9, ff_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_7 * df_6[k]
                  + pb_y[k] * ff_9[k];

        t_16[k] = pa_y[k] * dg_8[k];

        t_17[k] = pa_z[k] * dg_0[k];

        t_18[k] = f_7 * df_0[k]
                  + pb_z[k] * ff_10[k];

        t_19[k] = pa_z[k] * dg_3[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, df_2, df_3, df_4, \
                         dg_4, dg_5, dg_7, ff_11, ff_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * df_2[k]
                  + pa_z[k] * dg_4[k];

        t_21[k] = pa_z[k] * dg_5[k];

        t_22[k] = f_7 * df_3[k]
                  + pb_z[k] * ff_11[k];

        t_23[k] = f_5 * df_4[k]
                  + pa_z[k] * dg_7[k];

        t_24[k] = pb_y[k] * ff_13[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pa_z, df_6, df_8, df_10, df_11, dg_8, \
                         dg_16, dg_18, dg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_6 * df_6[k]
                  + pa_z[k] * dg_8[k];

        t_26[k] = f_6 * df_8[k]
                  + pa_x[k] * dg_16[k];

        t_27[k] = f_5 * df_10[k]
                  + pa_x[k] * dg_18[k];

        t_28[k] = f_5 * df_11[k]
                  + pa_x[k] * dg_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, pa_x, pa_y, pb_x, df_12, dg_13, \
                         dg_21, dg_23, dg_24, dg_25, ff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_7 * df_12[k]
                  + pb_x[k] * ff_16[k];

        t_30[k] = pa_x[k] * dg_21[k];

        t_31[k] = pa_x[k] * dg_23[k];

        t_32[k] = pa_x[k] * dg_24[k];

        t_33[k] = pa_x[k] * dg_25[k];

        t_34[k] = pa_y[k] * dg_13[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, t_40, pa_x, pa_y, pa_z, df_17, dg_10, \
                         dg_14, dg_27, dg_28, dg_29, dg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_z[k] * dg_10[k];

        t_36[k] = pa_y[k] * dg_14[k];

        t_37[k] = pa_x[k] * dg_27[k];

        t_38[k] = pa_x[k] * dg_28[k];

        t_39[k] = pa_x[k] * dg_29[k];

        t_40[k] = f_6 * df_17[k]
                  + pa_x[k] * dg_31[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pb_x, pb_z, df_7, df_19, df_20, df_23, \
                         dg_34, dg_35, ff_17, ff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * df_7[k]
                  + pb_z[k] * ff_17[k];

        t_42[k] = f_5 * df_19[k]
                  + pa_x[k] * dg_34[k];

        t_43[k] = f_5 * df_20[k]
                  + pa_x[k] * dg_35[k];

        t_44[k] = f_7 * df_23[k]
                  + pb_x[k] * ff_19[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, pa_x, pb_x, pb_z, dg_37, dg_38, \
                         dg_39, dg_41, fd0_5, fd1_7, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_x[k] * dg_37[k];

        t_46[k] = pa_x[k] * dg_38[k];

        t_47[k] = pa_x[k] * dg_39[k];

        t_48[k] = pa_x[k] * dg_41[k];

        t_49[k] = f_1 * fd0_5[k]
                  - f_2 * fd1_7[k]
                  + pb_x[k] * ff_20[k];

        t_50[k] = pb_z[k] * ff_20[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pb_x, pb_y, df_12, fd0_6, fd0_7, fd1_8, \
                         fd1_9, ff_22, ff_23, ff_24, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * fd0_6[k]
                  - f_4 * fd1_8[k]
                  + pb_x[k] * ff_22[k];

        t_52[k] = f_3 * fd0_7[k]
                  - f_4 * fd1_9[k]
                  + pb_x[k] * ff_23[k];

        t_53[k] = pb_x[k] * ff_24[k];

        t_54[k] = pb_x[k] * ff_26[k];

        t_55[k] = f_0 * df_12[k]
                  + f_1 * fd0_6[k]
                  - f_2 * fd1_8[k]
                  + pb_y[k] * ff_24[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pb_z, df_14, fd0_6, fd0_7, fd1_8, \
                         fd1_9, ff_24, ff_25, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * ff_24[k];

        t_57[k] = f_3 * fd0_6[k]
                  - f_4 * fd1_8[k]
                  + pb_z[k] * ff_25[k];

        t_58[k] = f_0 * df_14[k]
                  + pb_y[k] * ff_26[k];

        t_59[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_9[k]
                  + pb_z[k] * ff_26[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pa_z, pb_x, pb_z, df_8, df_9, dg_16, \
                         dg_18, dg_19, ff_27, ff_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_z[k] * dg_16[k];

        t_61[k] = f_7 * df_8[k]
                  + pb_z[k] * ff_27[k];

        t_62[k] = pa_z[k] * dg_18[k];

        t_63[k] = f_5 * df_9[k]
                  + pa_z[k] * dg_19[k];

        t_64[k] = pb_x[k] * ff_30[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_z, pb_y, pb_z, df_12, df_13, df_16, dg_21, \
                         dg_23, ff_29, ff_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pa_z[k] * dg_21[k];

        t_66[k] = f_7 * df_12[k]
                  + pb_z[k] * ff_29[k];

        t_67[k] = f_5 * df_13[k]
                  + pa_z[k] * dg_23[k];

        t_68[k] = f_5 * df_16[k]
                  + pb_y[k] * ff_30[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pa_y, pa_z, df_14, df_18, dg_25, dg_31, \
                         dg_33, dg_34, dg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_6 * df_14[k]
                  + pa_z[k] * dg_25[k];

        t_70[k] = pa_y[k] * dg_31[k];

        t_71[k] = pa_y[k] * dg_33[k];

        t_72[k] = f_5 * df_18[k]
                  + pa_y[k] * dg_34[k];

        t_73[k] = pa_y[k] * dg_35[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_y, pb_x, pb_z, df_15, df_21, df_22, dg_37, \
                         dg_39, ff_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pb_x[k] * ff_32[k];

        t_75[k] = f_6 * df_21[k]
                  + pa_y[k] * dg_37[k];

        t_76[k] = f_5 * df_15[k]
                  + pb_z[k] * ff_32[k];

        t_77[k] = f_5 * df_22[k]
                  + pa_y[k] * dg_39[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pb_x, pb_y, pb_z, df_17, df_23, \
                         dg_41, fd0_9, fd1_12, ff_34, ff_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_7 * df_23[k]
                  + pb_y[k] * ff_34[k];

        t_79[k] = pa_y[k] * dg_41[k];

        t_80[k] = f_1 * fd0_9[k]
                  - f_2 * fd1_12[k]
                  + pb_x[k] * ff_35[k];

        t_81[k] = pb_y[k] * ff_35[k];

        t_82[k] = f_0 * df_17[k]
                  + pb_z[k] * ff_35[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pb_x, pb_y, fd0_10, fd0_11, fd1_13, \
                         fd1_14, ff_37, ff_38, ff_39, ff_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_3 * fd0_10[k]
                  - f_4 * fd1_13[k]
                  + pb_x[k] * ff_37[k];

        t_84[k] = f_3 * fd0_11[k]
                  - f_4 * fd1_14[k]
                  + pb_x[k] * ff_38[k];

        t_85[k] = pb_x[k] * ff_39[k];

        t_86[k] = pb_x[k] * ff_41[k];

        t_87[k] = f_1 * fd0_10[k]
                  - f_2 * fd1_13[k]
                  + pb_y[k] * ff_39[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pb_y, pb_z, df_21, df_23, fd0_11, fd1_14, \
                         ff_39, ff_40, ff_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_0 * df_21[k]
                  + pb_z[k] * ff_39[k];

        t_89[k] = f_3 * fd0_11[k]
                  - f_4 * fd1_14[k]
                  + pb_y[k] * ff_40[k];

        t_90[k] = pb_y[k] * ff_41[k];

        t_91[k] = f_0 * df_23[k]
                  + f_1 * fd0_11[k]
                  - f_2 * fd1_14[k]
                  + pb_z[k] * ff_41[k];
    }
}

auto
compute_prim_fg_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t df, const size_t dg,
                                     const size_t fd0, const size_t fd1, const size_t ff,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / p;
    const auto f_7 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_9 = buffer.data(fd0 + 9);
    const auto *fd0_12 = buffer.data(fd0 + 12);
    const auto *fd0_13 = buffer.data(fd0 + 13);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_9 = buffer.data(fd1 + 9);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_15 = buffer.data(fd1 + 15);
    const auto *fd1_16 = buffer.data(fd1 + 16);
    const auto *fd1_17 = buffer.data(fd1 + 17);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, dg_0, fd0_1, fd0_2, fd1_1, \
                         fd1_2, ff_3, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_6[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_4[k];

        t_7[k] = pb_y[k] * ff_5[k];

        t_8[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_5[k];

        t_9[k] = pa_y[k] * dg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_z, df_0, df_1, df_3, dg_0, \
                         dg_3, dg_5, ff_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * df_1[k]
                  + pa_y[k] * dg_3[k];

        t_11[k] = f_6 * df_3[k]
                  + pa_y[k] * dg_5[k];

        t_12[k] = pa_z[k] * dg_0[k];

        t_13[k] = f_7 * df_0[k]
                  + pb_z[k] * ff_8[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_z, df_2, df_4, df_6, df_8, dg_4, \
                         dg_6, dg_8, dg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * df_2[k]
                  + pa_z[k] * dg_4[k];

        t_15[k] = f_5 * df_4[k]
                  + pa_z[k] * dg_6[k];

        t_16[k] = f_6 * df_6[k]
                  + pa_z[k] * dg_8[k];

        t_17[k] = f_6 * df_8[k]
                  + pa_x[k] * dg_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_x, pb_z, df_7, df_10, df_17, df_20, \
                         dg_10, dg_14, dg_18, dg_20, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * df_10[k]
                  + pa_x[k] * dg_10[k];

        t_19[k] = pa_x[k] * dg_14[k];

        t_20[k] = f_6 * df_17[k]
                  + pa_x[k] * dg_18[k];

        t_21[k] = f_5 * df_7[k]
                  + pb_z[k] * ff_12[k];

        t_22[k] = f_5 * df_20[k]
                  + pa_x[k] * dg_20[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_x, dg_26, fd0_7, fd0_8, fd0_9, \
                         fd1_8, fd1_9, fd1_10, ff_14, ff_15, ff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_x[k] * dg_26[k];

        t_24[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_8[k]
                  + pb_x[k] * ff_14[k];

        t_25[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_9[k]
                  + pb_x[k] * ff_15[k];

        t_26[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_10[k]
                  + pb_x[k] * ff_16[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, t_32, pb_x, pb_y, pb_z, df_12, df_14, \
                         fd0_8, fd1_9, ff_17, ff_18, ff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_x[k] * ff_17[k];

        t_28[k] = pb_x[k] * ff_19[k];

        t_29[k] = f_0 * df_12[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_9[k]
                  + pb_y[k] * ff_17[k];

        t_30[k] = pb_z[k] * ff_17[k];

        t_31[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_9[k]
                  + pb_z[k] * ff_18[k];

        t_32[k] = f_0 * df_14[k]
                  + pb_y[k] * ff_19[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pb_z, df_9, df_12, dg_11, dg_14, fd0_9, \
                         fd1_10, ff_19, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_1 * fd0_9[k]
                  - f_2 * fd1_10[k]
                  + pb_z[k] * ff_19[k];

        t_34[k] = f_5 * df_9[k]
                  + pa_z[k] * dg_11[k];

        t_35[k] = pa_z[k] * dg_14[k];

        t_36[k] = f_7 * df_12[k]
                  + pb_z[k] * ff_20[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_y, pa_z, pb_y, df_13, df_14, df_16, df_18, \
                         dg_16, dg_17, dg_19, ff_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * df_13[k]
                  + pa_z[k] * dg_16[k];

        t_38[k] = f_5 * df_16[k]
                  + pb_y[k] * ff_21[k];

        t_39[k] = f_6 * df_14[k]
                  + pa_z[k] * dg_17[k];

        t_40[k] = f_5 * df_18[k]
                  + pa_y[k] * dg_19[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_y, pb_y, pb_z, df_15, df_21, df_22, df_23, \
                         dg_23, dg_24, ff_22, ff_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_6 * df_21[k]
                  + pa_y[k] * dg_23[k];

        t_42[k] = f_5 * df_15[k]
                  + pb_z[k] * ff_22[k];

        t_43[k] = f_5 * df_22[k]
                  + pa_y[k] * dg_24[k];

        t_44[k] = f_7 * df_23[k]
                  + pb_y[k] * ff_23[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_y, pb_x, pb_z, df_17, dg_26, fd0_12, \
                         fd0_13, fd1_15, fd1_16, ff_24, ff_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_y[k] * dg_26[k];

        t_46[k] = f_1 * fd0_12[k]
                  - f_2 * fd1_15[k]
                  + pb_x[k] * ff_24[k];

        t_47[k] = f_0 * df_17[k]
                  + pb_z[k] * ff_24[k];

        t_48[k] = f_3 * fd0_13[k]
                  - f_4 * fd1_16[k]
                  + pb_x[k] * ff_25[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pb_x, pb_y, pb_z, df_21, fd0_13, \
                         fd0_14, fd1_16, fd1_17, ff_26, ff_27, ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_17[k]
                  + pb_x[k] * ff_26[k];

        t_50[k] = pb_x[k] * ff_27[k];

        t_51[k] = pb_x[k] * ff_29[k];

        t_52[k] = f_1 * fd0_13[k]
                  - f_2 * fd1_16[k]
                  + pb_y[k] * ff_27[k];

        t_53[k] = f_0 * df_21[k]
                  + pb_z[k] * ff_27[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_y, pb_z, df_23, fd0_14, fd1_17, ff_28, \
                         ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_17[k]
                  + pb_y[k] * ff_28[k];

        t_55[k] = pb_y[k] * ff_29[k];

        t_56[k] = f_0 * df_23[k]
                  + f_1 * fd0_14[k]
                  - f_2 * fd1_17[k]
                  + pb_z[k] * ff_29[k];
    }
}

auto
compute_prim_fg_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t df, const size_t dg,
                                     const size_t fd0, const size_t fd1, const size_t ff,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_17 = buffer.data(df + 17);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_32 = buffer.data(dg + 32);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_9 = buffer.data(fd0 + 9);
    const auto *fd0_12 = buffer.data(fd0 + 12);
    const auto *fd0_13 = buffer.data(fd0 + 13);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_9 = buffer.data(fd1 + 9);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_13 = buffer.data(fd1 + 13);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, dg_0, fd0_1, fd0_2, fd1_1, \
                         fd1_2, ff_3, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_6[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_4[k];

        t_7[k] = pb_y[k] * ff_5[k];

        t_8[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_5[k];

        t_9[k] = pa_y[k] * dg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_z, pb_x, dg_0, dg_17, dg_32, fd0_7, \
                         fd1_7, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * dg_0[k];

        t_11[k] = pa_x[k] * dg_17[k];

        t_12[k] = pa_x[k] * dg_32[k];

        t_13[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_7[k]
                  + pb_x[k] * ff_12[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pb_x, pb_y, df_9, fd0_8, fd0_9, fd1_8, \
                         fd1_9, ff_13, ff_14, ff_15, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_x[k] * ff_13[k];

        t_15[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_9[k]
                  + pb_x[k] * ff_14[k];

        t_16[k] = pb_x[k] * ff_15[k];

        t_17[k] = pb_x[k] * ff_17[k];

        t_18[k] = f_0 * df_9[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_8[k]
                  + pb_y[k] * ff_15[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_z, pb_z, dg_17, fd0_8, fd0_9, fd1_8, \
                         fd1_9, ff_15, ff_16, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_z[k] * ff_15[k];

        t_20[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_z[k] * ff_16[k];

        t_21[k] = f_1 * fd0_9[k]
                  - f_2 * fd1_9[k]
                  + pb_z[k] * ff_17[k];

        t_22[k] = pa_z[k] * dg_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pb_x, dg_32, fd0_12, fd0_13, fd0_14, \
                         fd1_12, fd1_13, fd1_14, ff_21, ff_22, ff_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_y[k] * dg_32[k];

        t_24[k] = f_1 * fd0_12[k]
                  - f_2 * fd1_12[k]
                  + pb_x[k] * ff_21[k];

        t_25[k] = f_3 * fd0_13[k]
                  - f_4 * fd1_13[k]
                  + pb_x[k] * ff_22[k];

        t_26[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_x[k] * ff_23[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pb_x, pb_y, fd0_13, fd0_14, fd1_13, \
                         fd1_14, ff_24, ff_25, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_x[k] * ff_24[k];

        t_28[k] = pb_x[k] * ff_26[k];

        t_29[k] = f_1 * fd0_13[k]
                  - f_2 * fd1_13[k]
                  + pb_y[k] * ff_24[k];

        t_30[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_y[k] * ff_25[k];

        t_31[k] = pb_y[k] * ff_26[k];
    }

#pragma omp simd aligned(t_32, pb_z, df_17, fd0_14, fd1_14, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * df_17[k]
                  + f_1 * fd0_14[k]
                  - f_2 * fd1_14[k]
                  + pb_z[k] * ff_26[k];
    }
}

auto
compute_prim_fg_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t df, const size_t dg,
                                     const size_t fd0, const size_t fd1, const size_t ff,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_32 = buffer.data(dg + 32);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_9 = buffer.data(fd0 + 9);
    const auto *fd0_12 = buffer.data(fd0 + 12);
    const auto *fd0_13 = buffer.data(fd0 + 13);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_9 = buffer.data(fd1 + 9);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_13 = buffer.data(fd1 + 13);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, fd0_1, fd0_2, fd1_1, fd1_2, \
                         ff_3, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_6[k] = pb_z[k] * ff_3[k];

        t_7[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_4[k];

        t_8[k] = pb_y[k] * ff_5[k];

        t_9[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pa_y, pa_z, df_2, df_3, dg_0, \
                         dg_4, dg_5, dg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * dg_0[k];

        t_11[k] = f_5 * df_3[k]
                  + pa_y[k] * dg_5[k];

        t_12[k] = pa_y[k] * dg_9[k];

        t_13[k] = pa_z[k] * dg_0[k];

        t_14[k] = f_6 * df_2[k]
                  + pa_z[k] * dg_4[k];

        t_15[k] = pa_z[k] * dg_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_z, pb_y, df_4, df_6, df_7, dg_7, \
                         dg_9, dg_11, ff_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_6 * df_4[k]
                  + pa_z[k] * dg_7[k];

        t_17[k] = pb_y[k] * ff_8[k];

        t_18[k] = f_5 * df_6[k]
                  + pa_z[k] * dg_9[k];

        t_19[k] = f_5 * df_7[k]
                  + pa_x[k] * dg_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_x, pb_x, df_9, df_11, df_14, df_17, \
                         dg_13, dg_17, dg_23, dg_26, ff_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_6 * df_9[k]
                  + pa_x[k] * dg_13[k];

        t_21[k] = f_7 * df_11[k]
                  + pb_x[k] * ff_10[k];

        t_22[k] = pa_x[k] * dg_17[k];

        t_23[k] = f_5 * df_14[k]
                  + pa_x[k] * dg_23[k];

        t_24[k] = f_6 * df_17[k]
                  + pa_x[k] * dg_26[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pb_x, pb_z, df_20, dg_32, fd0_7, fd1_7, \
                         ff_12, ff_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * df_20[k]
                  + pb_x[k] * ff_12[k];

        t_26[k] = pa_x[k] * dg_32[k];

        t_27[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_7[k]
                  + pb_x[k] * ff_13[k];

        t_28[k] = pb_z[k] * ff_13[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pb_x, pb_y, df_11, fd0_8, fd0_9, fd1_8, \
                         fd1_9, ff_14, ff_15, ff_16, ff_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_x[k] * ff_14[k];

        t_30[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_9[k]
                  + pb_x[k] * ff_15[k];

        t_31[k] = pb_x[k] * ff_16[k];

        t_32[k] = pb_x[k] * ff_18[k];

        t_33[k] = f_0 * df_11[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_8[k]
                  + pb_y[k] * ff_16[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_z, pb_z, dg_11, fd0_8, fd0_9, fd1_8, \
                         fd1_9, ff_16, ff_17, ff_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_z[k] * ff_16[k];

        t_35[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_z[k] * ff_17[k];

        t_36[k] = f_1 * fd0_9[k]
                  - f_2 * fd1_9[k]
                  + pb_z[k] * ff_18[k];

        t_37[k] = pa_z[k] * dg_11[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_z, pb_x, df_8, df_12, df_13, dg_14, \
                         dg_17, dg_19, dg_20, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * df_8[k]
                  + pa_z[k] * dg_14[k];

        t_39[k] = pb_x[k] * ff_20[k];

        t_40[k] = pa_z[k] * dg_17[k];

        t_41[k] = f_6 * df_12[k]
                  + pa_z[k] * dg_19[k];

        t_42[k] = f_5 * df_13[k]
                  + pa_z[k] * dg_20[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_y, pb_x, df_15, df_18, df_19, dg_25, \
                         dg_29, dg_30, ff_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_6 * df_15[k]
                  + pa_y[k] * dg_25[k];

        t_44[k] = pb_x[k] * ff_21[k];

        t_45[k] = f_5 * df_18[k]
                  + pa_y[k] * dg_29[k];

        t_46[k] = f_6 * df_19[k]
                  + pa_y[k] * dg_30[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pb_x, pb_y, df_20, dg_32, fd0_12, \
                         fd1_12, ff_23, ff_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_7 * df_20[k]
                  + pb_y[k] * ff_23[k];

        t_48[k] = pa_y[k] * dg_32[k];

        t_49[k] = f_1 * fd0_12[k]
                  - f_2 * fd1_12[k]
                  + pb_x[k] * ff_24[k];

        t_50[k] = pb_y[k] * ff_24[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pb_x, pb_y, fd0_13, fd0_14, fd1_13, \
                         fd1_14, ff_25, ff_26, ff_27, ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * fd0_13[k]
                  - f_4 * fd1_13[k]
                  + pb_x[k] * ff_25[k];

        t_52[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_x[k] * ff_26[k];

        t_53[k] = pb_x[k] * ff_27[k];

        t_54[k] = pb_x[k] * ff_29[k];

        t_55[k] = f_1 * fd0_13[k]
                  - f_2 * fd1_13[k]
                  + pb_y[k] * ff_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pb_y, pb_z, df_20, fd0_14, fd1_14, ff_28, \
                         ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_y[k] * ff_28[k];

        t_57[k] = pb_y[k] * ff_29[k];

        t_58[k] = f_0 * df_20[k]
                  + f_1 * fd0_14[k]
                  - f_2 * fd1_14[k]
                  + pb_z[k] * ff_29[k];
    }
}

auto
compute_prim_fg_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t df, const size_t dg,
                                     const size_t fd0, const size_t fd1, const size_t ff,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_17 = buffer.data(df + 17);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_9 = buffer.data(fd0 + 9);
    const auto *fd0_12 = buffer.data(fd0 + 12);
    const auto *fd0_13 = buffer.data(fd0 + 13);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_9 = buffer.data(fd1 + 9);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_13 = buffer.data(fd1 + 13);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, dg_0, fd0_1, fd0_2, fd1_1, \
                         fd1_2, ff_3, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_6[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_4[k];

        t_7[k] = pb_y[k] * ff_5[k];

        t_8[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_5[k];

        t_9[k] = pa_y[k] * dg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pa_y, pa_z, df_3, df_5, df_6, \
                         dg_0, dg_5, dg_8, dg_9, dg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * df_3[k]
                  + pa_y[k] * dg_5[k];

        t_11[k] = pa_z[k] * dg_0[k];

        t_12[k] = f_5 * df_5[k]
                  + pa_z[k] * dg_8[k];

        t_13[k] = f_5 * df_6[k]
                  + pa_x[k] * dg_9[k];

        t_14[k] = pa_x[k] * dg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pb_x, df_12, dg_18, dg_26, fd0_7, \
                         fd0_8, fd1_7, fd1_8, ff_10, ff_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * df_12[k]
                  + pa_x[k] * dg_18[k];

        t_16[k] = pa_x[k] * dg_26[k];

        t_17[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_7[k]
                  + pb_x[k] * ff_10[k];

        t_18[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_x[k] * ff_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pb_x, pb_y, pb_z, df_9, fd0_8, fd0_9, \
                         fd1_8, fd1_9, ff_12, ff_13, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_9[k]
                  + pb_x[k] * ff_12[k];

        t_20[k] = pb_x[k] * ff_13[k];

        t_21[k] = pb_x[k] * ff_15[k];

        t_22[k] = f_0 * df_9[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_8[k]
                  + pb_y[k] * ff_13[k];

        t_23[k] = pb_z[k] * ff_13[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_z, df_11, dg_14, dg_17, fd0_8, \
                         fd0_9, fd1_8, fd1_9, ff_14, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_z[k] * ff_14[k];

        t_25[k] = f_1 * fd0_9[k]
                  - f_2 * fd1_9[k]
                  + pb_z[k] * ff_15[k];

        t_26[k] = pa_z[k] * dg_14[k];

        t_27[k] = f_5 * df_11[k]
                  + pa_z[k] * dg_17[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pb_x, df_15, dg_23, dg_26, fd0_12, \
                         fd0_13, fd1_12, fd1_13, ff_18, ff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_5 * df_15[k]
                  + pa_y[k] * dg_23[k];

        t_29[k] = pa_y[k] * dg_26[k];

        t_30[k] = f_1 * fd0_12[k]
                  - f_2 * fd1_12[k]
                  + pb_x[k] * ff_18[k];

        t_31[k] = f_3 * fd0_13[k]
                  - f_4 * fd1_13[k]
                  + pb_x[k] * ff_19[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, pb_x, pb_y, fd0_13, fd0_14, \
                         fd1_13, fd1_14, ff_20, ff_21, ff_22, ff_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_x[k] * ff_20[k];

        t_33[k] = pb_x[k] * ff_21[k];

        t_34[k] = pb_x[k] * ff_23[k];

        t_35[k] = f_1 * fd0_13[k]
                  - f_2 * fd1_13[k]
                  + pb_y[k] * ff_21[k];

        t_36[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_y[k] * ff_22[k];

        t_37[k] = pb_y[k] * ff_23[k];
    }

#pragma omp simd aligned(t_38, pb_z, df_17, fd0_14, fd1_14, ff_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * df_17[k]
                  + f_1 * fd0_14[k]
                  - f_2 * fd1_14[k]
                  + pb_z[k] * ff_23[k];
    }
}

auto
compute_prim_fg_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_8 = buffer.data(fd0 + 8);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_8 = buffer.data(fd1 + 8);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_8 = buffer.data(ff + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, df_0, dg_0, dg_1, \
                         dg_2, fd0_0, fd1_0, ff_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pa_y[k] * dg_0[k];

        t_2[k] = pa_z[k] * dg_0[k];

        t_3[k] = pa_x[k] * dg_1[k];

        t_4[k] = pa_x[k] * dg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, df_1, dg_1, dg_2, fd0_5, fd1_5, \
                         ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * df_1[k]
                 + f_1 * fd0_5[k]
                 - f_2 * fd1_5[k]
                 + pb_y[k] * ff_5[k];

        t_6[k] = pa_z[k] * dg_1[k];

        t_7[k] = pa_y[k] * dg_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, df_2, fd0_8, fd1_8, ff_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * df_2[k]
                 + f_1 * fd0_8[k]
                 - f_2 * fd1_8[k]
                 + pb_z[k] * ff_8[k];
    }
}

auto
compute_prim_fg_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_8 = buffer.data(fd0 + 8);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_14 = buffer.data(fd1 + 14);
    const auto *fd1_23 = buffer.data(fd1 + 23);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_40 = buffer.data(ff + 40);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, df_0, dg_0, dg_1, \
                         dg_2, fd0_0, fd1_0, ff_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pa_y[k] * dg_0[k];

        t_2[k] = pa_z[k] * dg_0[k];

        t_3[k] = pa_x[k] * dg_1[k];

        t_4[k] = pa_x[k] * dg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, df_1, dg_1, dg_2, fd0_5, fd1_14, \
                         ff_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * df_1[k]
                 + f_1 * fd0_5[k]
                 - f_2 * fd1_14[k]
                 + pb_y[k] * ff_21[k];

        t_6[k] = pa_z[k] * dg_1[k];

        t_7[k] = pa_y[k] * dg_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, df_2, fd0_8, fd1_23, ff_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * df_2[k]
                 + f_1 * fd0_8[k]
                 - f_2 * fd1_23[k]
                 + pb_z[k] * ff_40[k];
    }
}

auto
compute_prim_fg_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_8 = buffer.data(df + 8);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_13 = buffer.data(fd0 + 13);
    const auto *fd0_14 = buffer.data(fd0 + 14);
    const auto *fd0_15 = buffer.data(fd0 + 15);
    const auto *fd0_21 = buffer.data(fd0 + 21);
    const auto *fd0_22 = buffer.data(fd0 + 22);
    const auto *fd0_23 = buffer.data(fd0 + 23);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_11 = buffer.data(fd1 + 11);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_13 = buffer.data(fd1 + 13);
    const auto *fd1_18 = buffer.data(fd1 + 18);
    const auto *fd1_19 = buffer.data(fd1 + 19);
    const auto *fd1_20 = buffer.data(fd1 + 20);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, df_0, fd0_0, fd0_1, fd1_0, \
                         fd1_1, ff_0, ff_1, ff_2, ff_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_2[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];

        t_3[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, pb_y, pb_z, dg_0, dg_1, \
                         fd0_2, fd1_2, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_4[k];

        t_5[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_5[k];

        t_6[k] = pa_y[k] * dg_0[k];

        t_7[k] = pa_z[k] * dg_0[k];

        t_8[k] = pa_x[k] * dg_1[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pb_x, dg_2, fd0_13, fd0_14, fd0_15, \
                         fd1_11, fd1_12, fd1_13, ff_10, ff_11, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_x[k] * dg_2[k];

        t_10[k] = f_1 * fd0_13[k]
                  - f_2 * fd1_11[k]
                  + pb_x[k] * ff_10[k];

        t_11[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_12[k]
                  + pb_x[k] * ff_11[k];

        t_12[k] = f_3 * fd0_15[k]
                  - f_4 * fd1_13[k]
                  + pb_x[k] * ff_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_z, pb_y, pb_z, df_4, dg_1, fd0_14, fd0_15, \
                         fd1_12, fd1_13, ff_13, ff_14, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * df_4[k]
                  + f_1 * fd0_14[k]
                  - f_2 * fd1_12[k]
                  + pb_y[k] * ff_13[k];

        t_14[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_12[k]
                  + pb_z[k] * ff_14[k];

        t_15[k] = f_1 * fd0_15[k]
                  - f_2 * fd1_13[k]
                  + pb_z[k] * ff_15[k];

        t_16[k] = pa_z[k] * dg_1[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pb_x, dg_2, fd0_21, fd0_22, fd0_23, \
                         fd1_18, fd1_19, fd1_20, ff_18, ff_19, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = pa_y[k] * dg_2[k];

        t_18[k] = f_1 * fd0_21[k]
                  - f_2 * fd1_18[k]
                  + pb_x[k] * ff_18[k];

        t_19[k] = f_3 * fd0_22[k]
                  - f_4 * fd1_19[k]
                  + pb_x[k] * ff_19[k];

        t_20[k] = f_3 * fd0_23[k]
                  - f_4 * fd1_20[k]
                  + pb_x[k] * ff_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_y, pb_z, df_8, fd0_22, fd0_23, fd1_19, fd1_20, \
                         ff_21, ff_22, ff_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * fd0_22[k]
                  - f_2 * fd1_19[k]
                  + pb_y[k] * ff_21[k];

        t_22[k] = f_3 * fd0_23[k]
                  - f_4 * fd1_20[k]
                  + pb_y[k] * ff_22[k];

        t_23[k] = f_0 * df_8[k]
                  + f_1 * fd0_23[k]
                  - f_2 * fd1_20[k]
                  + pb_z[k] * ff_23[k];
    }
}

auto
compute_prim_fg_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_8 = buffer.data(fd0 + 8);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_41 = buffer.data(ff + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, df_0, dg_0, dg_1, \
                         dg_2, fd0_0, fd1_0, ff_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pa_y[k] * dg_0[k];

        t_2[k] = pa_z[k] * dg_0[k];

        t_3[k] = pa_x[k] * dg_1[k];

        t_4[k] = pa_x[k] * dg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, df_1, dg_1, dg_2, fd0_5, fd1_8, \
                         ff_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * df_1[k]
                 + f_1 * fd0_5[k]
                 - f_2 * fd1_8[k]
                 + pb_y[k] * ff_24[k];

        t_6[k] = pa_z[k] * dg_1[k];

        t_7[k] = pa_y[k] * dg_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, df_2, fd0_8, fd1_14, ff_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * df_2[k]
                 + f_1 * fd0_8[k]
                 - f_2 * fd1_14[k]
                 + pb_z[k] * ff_41[k];
    }
}

auto
compute_prim_fg_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 2.0 / p;
    const auto f_7 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_9 = buffer.data(fd0 + 9);
    const auto *fd0_12 = buffer.data(fd0 + 12);
    const auto *fd0_13 = buffer.data(fd0 + 13);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_9 = buffer.data(fd1 + 9);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_15 = buffer.data(fd1 + 15);
    const auto *fd1_16 = buffer.data(fd1 + 16);
    const auto *fd1_17 = buffer.data(fd1 + 17);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, dg_0, fd0_1, fd0_2, fd1_1, \
                         fd1_2, ff_3, ff_5, ff_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_6[k] = pb_z[k] * ff_3[k];

        t_7[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_5[k];

        t_8[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_6[k];

        t_9[k] = pa_y[k] * dg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_z, df_0, df_1, df_3, dg_0, \
                         dg_1, dg_3, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * df_1[k]
                  + pa_y[k] * dg_1[k];

        t_11[k] = f_6 * df_3[k]
                  + pa_y[k] * dg_3[k];

        t_12[k] = pa_z[k] * dg_0[k];

        t_13[k] = f_7 * df_0[k]
                  + pb_z[k] * ff_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_z, df_2, df_4, df_5, df_7, dg_2, \
                         dg_4, dg_5, dg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * df_2[k]
                  + pa_z[k] * dg_2[k];

        t_15[k] = f_5 * df_4[k]
                  + pa_z[k] * dg_4[k];

        t_16[k] = f_6 * df_5[k]
                  + pa_z[k] * dg_5[k];

        t_17[k] = f_6 * df_7[k]
                  + pa_x[k] * dg_6[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pb_x, df_9, df_10, df_15, dg_7, dg_9, \
                         dg_12, ff_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * df_9[k]
                  + pa_x[k] * dg_7[k];

        t_19[k] = f_7 * df_10[k]
                  + pb_x[k] * ff_13[k];

        t_20[k] = pa_x[k] * dg_9[k];

        t_21[k] = f_6 * df_15[k]
                  + pa_x[k] * dg_12[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, df_6, df_17, df_20, dg_14, \
                         dg_17, ff_14, ff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * df_6[k]
                  + pb_z[k] * ff_14[k];

        t_23[k] = f_5 * df_17[k]
                  + pa_x[k] * dg_14[k];

        t_24[k] = f_7 * df_20[k]
                  + pb_x[k] * ff_16[k];

        t_25[k] = pa_x[k] * dg_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_x, pb_z, fd0_7, fd0_8, fd0_9, fd1_8, \
                         fd1_9, fd1_10, ff_17, ff_19, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_8[k]
                  + pb_x[k] * ff_17[k];

        t_27[k] = pb_z[k] * ff_17[k];

        t_28[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_9[k]
                  + pb_x[k] * ff_19[k];

        t_29[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_10[k]
                  + pb_x[k] * ff_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_x, pb_y, pb_z, df_10, df_12, fd0_8, \
                         fd1_9, ff_21, ff_22, ff_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_x[k] * ff_21[k];

        t_31[k] = f_0 * df_10[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_9[k]
                  + pb_y[k] * ff_21[k];

        t_32[k] = pb_z[k] * ff_21[k];

        t_33[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_9[k]
                  + pb_z[k] * ff_22[k];

        t_34[k] = f_0 * df_12[k]
                  + pb_y[k] * ff_23[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_z, pb_z, df_8, df_10, dg_8, dg_9, fd0_9, \
                         fd1_10, ff_23, ff_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * fd0_9[k]
                  - f_2 * fd1_10[k]
                  + pb_z[k] * ff_23[k];

        t_36[k] = f_5 * df_8[k]
                  + pa_z[k] * dg_8[k];

        t_37[k] = pa_z[k] * dg_9[k];

        t_38[k] = f_7 * df_10[k]
                  + pb_z[k] * ff_24[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, df_11, df_12, df_14, df_16, \
                         dg_10, dg_11, dg_13, ff_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_5 * df_11[k]
                  + pa_z[k] * dg_10[k];

        t_40[k] = f_5 * df_14[k]
                  + pb_y[k] * ff_25[k];

        t_41[k] = f_6 * df_12[k]
                  + pa_z[k] * dg_11[k];

        t_42[k] = f_5 * df_16[k]
                  + pa_y[k] * dg_13[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_y, pb_y, pb_z, df_13, df_18, df_19, df_20, \
                         dg_15, dg_16, ff_26, ff_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_6 * df_18[k]
                  + pa_y[k] * dg_15[k];

        t_44[k] = f_5 * df_13[k]
                  + pb_z[k] * ff_26[k];

        t_45[k] = f_5 * df_19[k]
                  + pa_y[k] * dg_16[k];

        t_46[k] = f_7 * df_20[k]
                  + pb_y[k] * ff_28[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pb_x, pb_y, pb_z, df_15, dg_17, fd0_12, \
                         fd1_15, ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pa_y[k] * dg_17[k];

        t_48[k] = f_1 * fd0_12[k]
                  - f_2 * fd1_15[k]
                  + pb_x[k] * ff_29[k];

        t_49[k] = pb_y[k] * ff_29[k];

        t_50[k] = f_0 * df_15[k]
                  + pb_z[k] * ff_29[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_x, pb_y, fd0_13, fd0_14, fd1_16, fd1_17, \
                         ff_31, ff_32, ff_33, ff_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * fd0_13[k]
                  - f_4 * fd1_16[k]
                  + pb_x[k] * ff_31[k];

        t_52[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_17[k]
                  + pb_x[k] * ff_32[k];

        t_53[k] = pb_x[k] * ff_35[k];

        t_54[k] = f_1 * fd0_13[k]
                  - f_2 * fd1_16[k]
                  + pb_y[k] * ff_33[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_y, pb_z, df_18, df_20, fd0_14, fd1_17, \
                         ff_33, ff_34, ff_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_0 * df_18[k]
                  + pb_z[k] * ff_33[k];

        t_56[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_17[k]
                  + pb_y[k] * ff_34[k];

        t_57[k] = pb_y[k] * ff_35[k];

        t_58[k] = f_0 * df_20[k]
                  + f_1 * fd0_14[k]
                  - f_2 * fd1_17[k]
                  + pb_z[k] * ff_35[k];
    }
}

auto
compute_prim_fg_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_17 = buffer.data(df + 17);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_9 = buffer.data(fd0 + 9);
    const auto *fd0_10 = buffer.data(fd0 + 10);
    const auto *fd0_15 = buffer.data(fd0 + 15);
    const auto *fd0_16 = buffer.data(fd0 + 16);
    const auto *fd0_17 = buffer.data(fd0 + 17);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_9 = buffer.data(fd1 + 9);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_13 = buffer.data(fd1 + 13);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, dg_0, fd0_1, fd0_2, fd1_1, \
                         fd1_2, ff_3, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_6[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_4[k];

        t_7[k] = pb_y[k] * ff_5[k];

        t_8[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_5[k];

        t_9[k] = pa_y[k] * dg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_z, pb_x, dg_0, dg_1, dg_2, fd0_8, \
                         fd1_7, ff_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * dg_0[k];

        t_11[k] = pa_x[k] * dg_1[k];

        t_12[k] = pa_x[k] * dg_2[k];

        t_13[k] = f_1 * fd0_8[k]
                  - f_2 * fd1_7[k]
                  + pb_x[k] * ff_10[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pb_x, pb_y, df_9, fd0_9, fd0_10, fd1_8, \
                         fd1_9, ff_11, ff_12, ff_13, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_8[k]
                  + pb_x[k] * ff_11[k];

        t_15[k] = f_3 * fd0_10[k]
                  - f_4 * fd1_9[k]
                  + pb_x[k] * ff_12[k];

        t_16[k] = pb_x[k] * ff_13[k];

        t_17[k] = pb_x[k] * ff_15[k];

        t_18[k] = f_0 * df_9[k]
                  + f_1 * fd0_9[k]
                  - f_2 * fd1_8[k]
                  + pb_y[k] * ff_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_z, pb_z, dg_1, fd0_9, fd0_10, fd1_8, \
                         fd1_9, ff_13, ff_14, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_z[k] * ff_13[k];

        t_20[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_8[k]
                  + pb_z[k] * ff_14[k];

        t_21[k] = f_1 * fd0_10[k]
                  - f_2 * fd1_9[k]
                  + pb_z[k] * ff_15[k];

        t_22[k] = pa_z[k] * dg_1[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pb_x, dg_2, fd0_15, fd0_16, fd0_17, \
                         fd1_12, fd1_13, fd1_14, ff_18, ff_19, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_y[k] * dg_2[k];

        t_24[k] = f_1 * fd0_15[k]
                  - f_2 * fd1_12[k]
                  + pb_x[k] * ff_18[k];

        t_25[k] = f_3 * fd0_16[k]
                  - f_4 * fd1_13[k]
                  + pb_x[k] * ff_19[k];

        t_26[k] = f_3 * fd0_17[k]
                  - f_4 * fd1_14[k]
                  + pb_x[k] * ff_20[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pb_x, pb_y, fd0_16, fd0_17, fd1_13, \
                         fd1_14, ff_21, ff_22, ff_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_x[k] * ff_21[k];

        t_28[k] = pb_x[k] * ff_23[k];

        t_29[k] = f_1 * fd0_16[k]
                  - f_2 * fd1_13[k]
                  + pb_y[k] * ff_21[k];

        t_30[k] = f_3 * fd0_17[k]
                  - f_4 * fd1_14[k]
                  + pb_y[k] * ff_22[k];

        t_31[k] = pb_y[k] * ff_23[k];
    }

#pragma omp simd aligned(t_32, pb_z, df_17, fd0_17, fd1_14, ff_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * df_17[k]
                  + f_1 * fd0_17[k]
                  - f_2 * fd1_14[k]
                  + pb_z[k] * ff_23[k];
    }
}

auto
compute_prim_fg_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_8 = buffer.data(fd0 + 8);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_26 = buffer.data(ff + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, df_0, dg_0, dg_1, \
                         dg_2, fd0_0, fd1_0, ff_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pa_y[k] * dg_0[k];

        t_2[k] = pa_z[k] * dg_0[k];

        t_3[k] = pa_x[k] * dg_1[k];

        t_4[k] = pa_x[k] * dg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, df_1, dg_1, dg_2, fd0_5, fd1_8, \
                         ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * df_1[k]
                 + f_1 * fd0_5[k]
                 - f_2 * fd1_8[k]
                 + pb_y[k] * ff_15[k];

        t_6[k] = pa_z[k] * dg_1[k];

        t_7[k] = pa_y[k] * dg_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, df_2, fd0_8, fd1_14, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * df_2[k]
                 + f_1 * fd0_8[k]
                 - f_2 * fd1_14[k]
                 + pb_z[k] * ff_26[k];
    }
}

auto
compute_prim_fg_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_29 = buffer.data(dg + 29);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_9 = buffer.data(fd0 + 9);
    const auto *fd0_12 = buffer.data(fd0 + 12);
    const auto *fd0_13 = buffer.data(fd0 + 13);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_9 = buffer.data(fd1 + 9);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_13 = buffer.data(fd1 + 13);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, fd0_1, fd0_2, fd1_1, fd1_2, \
                         ff_3, ff_5, ff_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_6[k] = pb_z[k] * ff_3[k];

        t_7[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_5[k];

        t_8[k] = pb_y[k] * ff_6[k];

        t_9[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pa_y, pa_z, df_3, df_4, dg_0, \
                         dg_5, dg_7, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * dg_0[k];

        t_11[k] = f_5 * df_3[k]
                  + pa_y[k] * dg_5[k];

        t_12[k] = pa_y[k] * dg_8[k];

        t_13[k] = pa_z[k] * dg_0[k];

        t_14[k] = pa_z[k] * dg_5[k];

        t_15[k] = f_6 * df_4[k]
                  + pa_z[k] * dg_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_z, pb_y, df_6, df_7, df_9, dg_8, \
                         dg_10, dg_12, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * ff_9[k];

        t_17[k] = f_5 * df_6[k]
                  + pa_z[k] * dg_8[k];

        t_18[k] = f_5 * df_7[k]
                  + pa_x[k] * dg_10[k];

        t_19[k] = f_6 * df_9[k]
                  + pa_x[k] * dg_12[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_x, pb_x, df_11, df_14, df_17, df_20, \
                         dg_15, dg_21, dg_24, ff_11, ff_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_7 * df_11[k]
                  + pb_x[k] * ff_11[k];

        t_21[k] = pa_x[k] * dg_15[k];

        t_22[k] = f_5 * df_14[k]
                  + pa_x[k] * dg_21[k];

        t_23[k] = f_6 * df_17[k]
                  + pa_x[k] * dg_24[k];

        t_24[k] = f_7 * df_20[k]
                  + pb_x[k] * ff_13[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pb_x, pb_z, dg_29, fd0_7, fd0_8, fd1_7, \
                         fd1_8, ff_14, ff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_x[k] * dg_29[k];

        t_26[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_7[k]
                  + pb_x[k] * ff_14[k];

        t_27[k] = pb_z[k] * ff_14[k];

        t_28[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_x[k] * ff_16[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pb_x, pb_y, pb_z, df_11, fd0_8, fd0_9, \
                         fd1_8, fd1_9, ff_17, ff_18, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_9[k]
                  + pb_x[k] * ff_17[k];

        t_30[k] = pb_x[k] * ff_18[k];

        t_31[k] = pb_x[k] * ff_20[k];

        t_32[k] = f_0 * df_11[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_8[k]
                  + pb_y[k] * ff_18[k];

        t_33[k] = pb_z[k] * ff_18[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_z, pb_z, df_8, dg_10, dg_13, fd0_8, fd0_9, \
                         fd1_8, fd1_9, ff_19, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_z[k] * ff_19[k];

        t_35[k] = f_1 * fd0_9[k]
                  - f_2 * fd1_9[k]
                  + pb_z[k] * ff_20[k];

        t_36[k] = pa_z[k] * dg_10[k];

        t_37[k] = f_6 * df_8[k]
                  + pa_z[k] * dg_13[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_x, df_13, df_15, dg_15, \
                         dg_18, dg_23, ff_22, ff_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_x[k] * ff_22[k];

        t_39[k] = pa_z[k] * dg_15[k];

        t_40[k] = f_5 * df_13[k]
                  + pa_z[k] * dg_18[k];

        t_41[k] = f_6 * df_15[k]
                  + pa_y[k] * dg_23[k];

        t_42[k] = pb_x[k] * ff_23[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_y, pb_y, df_18, df_19, df_20, dg_26, \
                         dg_27, dg_29, ff_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_5 * df_18[k]
                  + pa_y[k] * dg_26[k];

        t_44[k] = f_6 * df_19[k]
                  + pa_y[k] * dg_27[k];

        t_45[k] = f_7 * df_20[k]
                  + pb_y[k] * ff_25[k];

        t_46[k] = pa_y[k] * dg_29[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_x, pb_y, fd0_12, fd0_13, fd0_14, fd1_12, \
                         fd1_13, fd1_14, ff_26, ff_28, ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_1 * fd0_12[k]
                  - f_2 * fd1_12[k]
                  + pb_x[k] * ff_26[k];

        t_48[k] = pb_y[k] * ff_26[k];

        t_49[k] = f_3 * fd0_13[k]
                  - f_4 * fd1_13[k]
                  + pb_x[k] * ff_28[k];

        t_50[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_x[k] * ff_29[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pb_x, pb_y, fd0_13, fd0_14, fd1_13, \
                         fd1_14, ff_30, ff_31, ff_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_x[k] * ff_30[k];

        t_52[k] = pb_x[k] * ff_32[k];

        t_53[k] = f_1 * fd0_13[k]
                  - f_2 * fd1_13[k]
                  + pb_y[k] * ff_30[k];

        t_54[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_y[k] * ff_31[k];

        t_55[k] = pb_y[k] * ff_32[k];
    }

#pragma omp simd aligned(t_56, pb_z, df_20, fd0_14, fd1_14, ff_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * df_20[k]
                  + f_1 * fd0_14[k]
                  - f_2 * fd1_14[k]
                  + pb_z[k] * ff_32[k];
    }
}

auto
compute_prim_fg_electron_repulsion_18(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_26 = buffer.data(dg + 26);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_9 = buffer.data(fd0 + 9);
    const auto *fd0_12 = buffer.data(fd0 + 12);
    const auto *fd0_13 = buffer.data(fd0 + 13);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_9 = buffer.data(fd1 + 9);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_13 = buffer.data(fd1 + 13);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, fd0_1, fd0_2, fd1_1, fd1_2, \
                         ff_3, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_6[k] = pb_z[k] * ff_3[k];

        t_7[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_4[k];

        t_8[k] = pb_y[k] * ff_5[k];

        t_9[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pa_y, pa_z, df_3, df_5, df_6, \
                         dg_0, dg_5, dg_8, dg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * dg_0[k];

        t_11[k] = f_5 * df_3[k]
                  + pa_y[k] * dg_5[k];

        t_12[k] = pa_z[k] * dg_0[k];

        t_13[k] = f_5 * df_5[k]
                  + pa_z[k] * dg_8[k];

        t_14[k] = f_5 * df_6[k]
                  + pa_x[k] * dg_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_x, pb_x, df_7, df_9, df_12, df_14, \
                         dg_10, dg_14, dg_18, dg_20, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * df_7[k]
                  + pa_x[k] * dg_10[k];

        t_16[k] = f_7 * df_9[k]
                  + pb_x[k] * ff_9[k];

        t_17[k] = pa_x[k] * dg_14[k];

        t_18[k] = f_5 * df_12[k]
                  + pa_x[k] * dg_18[k];

        t_19[k] = f_6 * df_14[k]
                  + pa_x[k] * dg_20[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, df_17, dg_26, fd0_7, fd1_7, \
                         ff_11, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_7 * df_17[k]
                  + pb_x[k] * ff_11[k];

        t_21[k] = pa_x[k] * dg_26[k];

        t_22[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_7[k]
                  + pb_x[k] * ff_12[k];

        t_23[k] = pb_z[k] * ff_12[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pb_x, pb_y, df_9, fd0_8, fd0_9, fd1_8, \
                         fd1_9, ff_13, ff_14, ff_15, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_x[k] * ff_13[k];

        t_25[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_9[k]
                  + pb_x[k] * ff_14[k];

        t_26[k] = pb_x[k] * ff_15[k];

        t_27[k] = pb_x[k] * ff_17[k];

        t_28[k] = f_0 * df_9[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_8[k]
                  + pb_y[k] * ff_15[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_z, pb_z, dg_14, fd0_8, fd0_9, fd1_8, \
                         fd1_9, ff_15, ff_16, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_z[k] * ff_15[k];

        t_30[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_z[k] * ff_16[k];

        t_31[k] = f_1 * fd0_9[k]
                  - f_2 * fd1_9[k]
                  + pb_z[k] * ff_17[k];

        t_32[k] = pa_z[k] * dg_14[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_y, pa_z, pb_y, df_11, df_15, df_16, df_17, \
                         dg_17, dg_23, dg_24, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_5 * df_11[k]
                  + pa_z[k] * dg_17[k];

        t_34[k] = f_5 * df_15[k]
                  + pa_y[k] * dg_23[k];

        t_35[k] = f_6 * df_16[k]
                  + pa_y[k] * dg_24[k];

        t_36[k] = f_7 * df_17[k]
                  + pb_y[k] * ff_20[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_y, pb_x, pb_y, dg_26, fd0_12, fd0_13, \
                         fd1_12, fd1_13, ff_21, ff_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pa_y[k] * dg_26[k];

        t_38[k] = f_1 * fd0_12[k]
                  - f_2 * fd1_12[k]
                  + pb_x[k] * ff_21[k];

        t_39[k] = pb_y[k] * ff_21[k];

        t_40[k] = f_3 * fd0_13[k]
                  - f_4 * fd1_13[k]
                  + pb_x[k] * ff_22[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, t_46, pb_x, pb_y, fd0_13, fd0_14, \
                         fd1_13, fd1_14, ff_23, ff_24, ff_25, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_x[k] * ff_23[k];

        t_42[k] = pb_x[k] * ff_24[k];

        t_43[k] = pb_x[k] * ff_26[k];

        t_44[k] = f_1 * fd0_13[k]
                  - f_2 * fd1_13[k]
                  + pb_y[k] * ff_24[k];

        t_45[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_y[k] * ff_25[k];

        t_46[k] = pb_y[k] * ff_26[k];
    }

#pragma omp simd aligned(t_47, pb_z, df_17, fd0_14, fd1_14, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * df_17[k]
                  + f_1 * fd0_14[k]
                  - f_2 * fd1_14[k]
                  + pb_z[k] * ff_26[k];
    }
}

auto
compute_prim_fg_electron_repulsion_19(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_17 = buffer.data(df + 17);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_9 = buffer.data(fd0 + 9);
    const auto *fd0_12 = buffer.data(fd0 + 12);
    const auto *fd0_13 = buffer.data(fd0 + 13);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_9 = buffer.data(fd1 + 9);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_13 = buffer.data(fd1 + 13);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, dg_0, fd0_1, fd0_2, fd1_1, \
                         fd1_2, ff_3, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_6[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_4[k];

        t_7[k] = pb_y[k] * ff_5[k];

        t_8[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_5[k];

        t_9[k] = pa_y[k] * dg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_z, pb_x, dg_0, dg_1, dg_2, fd0_7, \
                         fd1_7, ff_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * dg_0[k];

        t_11[k] = pa_x[k] * dg_1[k];

        t_12[k] = pa_x[k] * dg_2[k];

        t_13[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_7[k]
                  + pb_x[k] * ff_10[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pb_x, pb_y, df_9, fd0_8, fd0_9, fd1_8, \
                         fd1_9, ff_11, ff_12, ff_13, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_x[k] * ff_11[k];

        t_15[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_9[k]
                  + pb_x[k] * ff_12[k];

        t_16[k] = pb_x[k] * ff_13[k];

        t_17[k] = pb_x[k] * ff_15[k];

        t_18[k] = f_0 * df_9[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_8[k]
                  + pb_y[k] * ff_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_z, pb_z, dg_1, fd0_8, fd0_9, fd1_8, fd1_9, \
                         ff_13, ff_14, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_z[k] * ff_13[k];

        t_20[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_z[k] * ff_14[k];

        t_21[k] = f_1 * fd0_9[k]
                  - f_2 * fd1_9[k]
                  + pb_z[k] * ff_15[k];

        t_22[k] = pa_z[k] * dg_1[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pb_x, dg_2, fd0_12, fd0_13, fd0_14, \
                         fd1_12, fd1_13, fd1_14, ff_18, ff_19, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_y[k] * dg_2[k];

        t_24[k] = f_1 * fd0_12[k]
                  - f_2 * fd1_12[k]
                  + pb_x[k] * ff_18[k];

        t_25[k] = f_3 * fd0_13[k]
                  - f_4 * fd1_13[k]
                  + pb_x[k] * ff_19[k];

        t_26[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_x[k] * ff_20[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pb_x, pb_y, fd0_13, fd0_14, fd1_13, \
                         fd1_14, ff_21, ff_22, ff_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_x[k] * ff_21[k];

        t_28[k] = pb_x[k] * ff_23[k];

        t_29[k] = f_1 * fd0_13[k]
                  - f_2 * fd1_13[k]
                  + pb_y[k] * ff_21[k];

        t_30[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_y[k] * ff_22[k];

        t_31[k] = pb_y[k] * ff_23[k];
    }

#pragma omp simd aligned(t_32, pb_z, df_17, fd0_14, fd1_14, ff_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * df_17[k]
                  + f_1 * fd0_14[k]
                  - f_2 * fd1_14[k]
                  + pb_z[k] * ff_23[k];
    }
}

auto
compute_prim_fg_electron_repulsion_20(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_8 = buffer.data(fd0 + 8);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_8 = buffer.data(fd1 + 8);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_11 = buffer.data(ff + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, df_0, dg_0, dg_1, \
                         dg_2, fd0_0, fd1_0, ff_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pa_y[k] * dg_0[k];

        t_2[k] = pa_z[k] * dg_0[k];

        t_3[k] = pa_x[k] * dg_1[k];

        t_4[k] = pa_x[k] * dg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, df_1, dg_1, dg_2, fd0_5, fd1_5, \
                         ff_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * df_1[k]
                 + f_1 * fd0_5[k]
                 - f_2 * fd1_5[k]
                 + pb_y[k] * ff_7[k];

        t_6[k] = pa_z[k] * dg_1[k];

        t_7[k] = pa_y[k] * dg_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, df_2, fd0_8, fd1_8, ff_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * df_2[k]
                 + f_1 * fd0_8[k]
                 - f_2 * fd1_8[k]
                 + pb_z[k] * ff_11[k];
    }
}

auto
compute_prim_fg_electron_repulsion_21(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_8 = buffer.data(fd0 + 8);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_20 = buffer.data(fd1 + 20);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_17 = buffer.data(ff + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, df_0, dg_0, dg_1, \
                         dg_2, fd0_0, fd1_0, ff_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pa_y[k] * dg_0[k];

        t_2[k] = pa_z[k] * dg_0[k];

        t_3[k] = pa_x[k] * dg_1[k];

        t_4[k] = pa_x[k] * dg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, df_1, dg_1, dg_2, fd0_5, fd1_12, \
                         ff_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * df_1[k]
                 + f_1 * fd0_5[k]
                 - f_2 * fd1_12[k]
                 + pb_y[k] * ff_10[k];

        t_6[k] = pa_z[k] * dg_1[k];

        t_7[k] = pa_y[k] * dg_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, df_2, fd0_8, fd1_20, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * df_2[k]
                 + f_1 * fd0_8[k]
                 - f_2 * fd1_20[k]
                 + pb_z[k] * ff_17[k];
    }
}

auto
compute_prim_fg_electron_repulsion_22(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_12 = buffer.data(fd0 + 12);
    const auto *fd0_20 = buffer.data(fd0 + 20);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_8 = buffer.data(ff + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, df_0, dg_0, dg_1, \
                         dg_2, fd0_0, fd1_0, ff_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pa_y[k] * dg_0[k];

        t_2[k] = pa_z[k] * dg_0[k];

        t_3[k] = pa_x[k] * dg_1[k];

        t_4[k] = pa_x[k] * dg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, df_1, dg_1, dg_2, fd0_12, fd1_8, \
                         ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * df_1[k]
                 + f_1 * fd0_12[k]
                 - f_2 * fd1_8[k]
                 + pb_y[k] * ff_5[k];

        t_6[k] = pa_z[k] * dg_1[k];

        t_7[k] = pa_y[k] * dg_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, df_2, fd0_20, fd1_14, ff_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * df_2[k]
                 + f_1 * fd0_20[k]
                 - f_2 * fd1_14[k]
                 + pb_z[k] * ff_8[k];
    }
}

auto
compute_prim_fg_electron_repulsion_23(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_8 = buffer.data(fd0 + 8);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_8 = buffer.data(fd1 + 8);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, df_0, df_1, dg_0, dg_1, \
                         fd0_0, fd1_0, ff_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pa_y[k] * dg_0[k];

        t_2[k] = pa_z[k] * dg_0[k];

        t_3[k] = f_3 * df_1[k]
                 + pa_x[k] * dg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, df_2, df_3, df_5, dg_2, dg_3, \
                         dg_5, ff_4, ff_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * df_2[k]
                 + pb_x[k] * ff_4[k];

        t_5[k] = pa_x[k] * dg_2[k];

        t_6[k] = f_3 * df_3[k]
                 + pa_x[k] * dg_3[k];

        t_7[k] = f_4 * df_5[k]
                 + pb_x[k] * ff_6[k];

        t_8[k] = pa_x[k] * dg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_y, pa_z, pb_y, df_2, df_4, df_5, dg_2, \
                         dg_4, fd0_5, fd1_5, ff_7, ff_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * df_2[k]
                 + f_1 * fd0_5[k]
                 - f_2 * fd1_5[k]
                 + pb_y[k] * ff_7[k];

        t_10[k] = pa_z[k] * dg_2[k];

        t_11[k] = f_3 * df_4[k]
                  + pa_y[k] * dg_4[k];

        t_12[k] = f_4 * df_5[k]
                  + pb_y[k] * ff_10[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_y, pb_z, df_5, dg_5, fd0_8, fd1_8, \
                         ff_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * dg_5[k];

        t_14[k] = f_0 * df_5[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_8[k]
                  + pb_z[k] * ff_11[k];
    }
}

auto
compute_prim_fg_electron_repulsion_24(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_8 = buffer.data(fd0 + 8);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_32 = buffer.data(ff + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, df_0, df_1, dg_0, dg_1, \
                         fd0_0, fd1_0, ff_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pa_y[k] * dg_0[k];

        t_2[k] = pa_z[k] * dg_0[k];

        t_3[k] = f_3 * df_1[k]
                 + pa_x[k] * dg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, df_2, df_3, df_5, dg_2, dg_3, \
                         dg_5, ff_12, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * df_2[k]
                 + pb_x[k] * ff_12[k];

        t_5[k] = pa_x[k] * dg_2[k];

        t_6[k] = f_3 * df_3[k]
                 + pa_x[k] * dg_3[k];

        t_7[k] = f_4 * df_5[k]
                 + pb_x[k] * ff_15[k];

        t_8[k] = pa_x[k] * dg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_y, pa_z, pb_y, df_2, df_4, df_5, dg_2, \
                         dg_4, fd0_5, fd1_8, ff_19, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * df_2[k]
                 + f_1 * fd0_5[k]
                 - f_2 * fd1_8[k]
                 + pb_y[k] * ff_19[k];

        t_10[k] = pa_z[k] * dg_2[k];

        t_11[k] = f_3 * df_4[k]
                  + pa_y[k] * dg_4[k];

        t_12[k] = f_4 * df_5[k]
                  + pb_y[k] * ff_26[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_y, pb_z, df_5, dg_5, fd0_8, fd1_14, \
                         ff_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * dg_5[k];

        t_14[k] = f_0 * df_5[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_14[k]
                  + pb_z[k] * ff_32[k];
    }
}

auto
compute_prim_fg_electron_repulsion_25(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_9 = buffer.data(fd0 + 9);
    const auto *fd0_12 = buffer.data(fd0 + 12);
    const auto *fd0_13 = buffer.data(fd0 + 13);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_9 = buffer.data(fd1 + 9);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_13 = buffer.data(fd1 + 13);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, df_0, fd0_0, fd0_1, fd1_0, \
                         fd1_1, ff_0, ff_1, ff_2, ff_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_2[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];

        t_3[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, pb_y, pb_z, dg_0, fd0_2, fd1_2, ff_4, \
                         ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_4[k];

        t_5[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_5[k];

        t_6[k] = pa_y[k] * dg_0[k];

        t_7[k] = pa_z[k] * dg_0[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, pa_x, pb_x, df_4, df_5, df_8, df_11, \
                         dg_1, dg_2, dg_3, ff_9, ff_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * df_4[k]
                 + pa_x[k] * dg_1[k];

        t_9[k] = f_6 * df_5[k]
                 + pb_x[k] * ff_9[k];

        t_10[k] = pa_x[k] * dg_2[k];

        t_11[k] = f_5 * df_8[k]
                  + pa_x[k] * dg_3[k];

        t_12[k] = f_6 * df_11[k]
                  + pb_x[k] * ff_11[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pb_x, dg_5, fd0_7, fd0_8, fd0_9, fd1_7, \
                         fd1_8, fd1_9, ff_12, ff_13, ff_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_x[k] * dg_5[k];

        t_14[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_7[k]
                  + pb_x[k] * ff_12[k];

        t_15[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_x[k] * ff_13[k];

        t_16[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_9[k]
                  + pb_x[k] * ff_14[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_z, pb_y, pb_z, df_5, dg_2, fd0_8, fd0_9, \
                         fd1_8, fd1_9, ff_15, ff_16, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * df_5[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_8[k]
                  + pb_y[k] * ff_15[k];

        t_18[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_z[k] * ff_16[k];

        t_19[k] = f_1 * fd0_9[k]
                  - f_2 * fd1_9[k]
                  + pb_z[k] * ff_17[k];

        t_20[k] = pa_z[k] * dg_2[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_y, pb_x, pb_y, df_10, df_11, dg_4, dg_5, \
                         fd0_12, fd1_12, ff_20, ff_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * df_10[k]
                  + pa_y[k] * dg_4[k];

        t_22[k] = f_6 * df_11[k]
                  + pb_y[k] * ff_20[k];

        t_23[k] = pa_y[k] * dg_5[k];

        t_24[k] = f_1 * fd0_12[k]
                  - f_2 * fd1_12[k]
                  + pb_x[k] * ff_21[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pb_x, pb_y, fd0_13, fd0_14, fd1_13, fd1_14, \
                         ff_22, ff_23, ff_24, ff_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * fd0_13[k]
                  - f_4 * fd1_13[k]
                  + pb_x[k] * ff_22[k];

        t_26[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_x[k] * ff_23[k];

        t_27[k] = f_1 * fd0_13[k]
                  - f_2 * fd1_13[k]
                  + pb_y[k] * ff_24[k];

        t_28[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_y[k] * ff_25[k];
    }

#pragma omp simd aligned(t_29, pb_z, df_11, fd0_14, fd1_14, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * df_11[k]
                  + f_1 * fd0_14[k]
                  - f_2 * fd1_14[k]
                  + pb_z[k] * ff_26[k];
    }
}

auto
compute_prim_fg_electron_repulsion_26(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_5 = buffer.data(df + 5);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_8 = buffer.data(ff + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, df_0, dg_0, dg_1, \
                         dg_2, fd0_0, fd1_0, ff_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pa_y[k] * dg_0[k];

        t_2[k] = pa_z[k] * dg_0[k];

        t_3[k] = pa_x[k] * dg_1[k];

        t_4[k] = pa_x[k] * dg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, df_2, dg_1, dg_2, fd0_8, fd1_8, \
                         ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * df_2[k]
                 + f_1 * fd0_8[k]
                 - f_2 * fd1_8[k]
                 + pb_y[k] * ff_5[k];

        t_6[k] = pa_z[k] * dg_1[k];

        t_7[k] = pa_y[k] * dg_2[k];
    }

#pragma omp simd aligned(t_8, pb_z, df_5, fd0_14, fd1_14, ff_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * df_5[k]
                 + f_1 * fd0_14[k]
                 - f_2 * fd1_14[k]
                 + pb_z[k] * ff_8[k];
    }
}

auto
compute_prim_fg_electron_repulsion_27(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_5 = buffer.data(df + 5);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_8 = buffer.data(fd0 + 8);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_8 = buffer.data(fd1 + 8);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_11 = buffer.data(ff + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, df_0, dg_0, dg_3, \
                         dg_8, fd0_0, fd1_0, ff_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pa_y[k] * dg_0[k];

        t_2[k] = pa_z[k] * dg_0[k];

        t_3[k] = pa_x[k] * dg_3[k];

        t_4[k] = pa_x[k] * dg_8[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_y, pa_z, pb_y, df_2, dg_3, dg_8, fd0_5, fd1_5, \
                         ff_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * df_2[k]
                 + f_1 * fd0_5[k]
                 - f_2 * fd1_5[k]
                 + pb_y[k] * ff_7[k];

        t_6[k] = pa_z[k] * dg_3[k];

        t_7[k] = pa_y[k] * dg_8[k];
    }

#pragma omp simd aligned(t_8, pb_z, df_5, fd0_8, fd1_8, ff_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * df_5[k]
                 + f_1 * fd0_8[k]
                 - f_2 * fd1_8[k]
                 + pb_z[k] * ff_11[k];
    }
}

auto
compute_prim_fg_electron_repulsion_28(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_8 = buffer.data(fd0 + 8);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_29 = buffer.data(ff + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, df_0, df_1, dg_0, dg_1, \
                         fd0_0, fd1_0, ff_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pa_y[k] * dg_0[k];

        t_2[k] = pa_z[k] * dg_0[k];

        t_3[k] = f_3 * df_1[k]
                 + pa_x[k] * dg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pb_x, df_2, df_3, df_5, dg_3, dg_4, \
                         dg_8, ff_10, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * df_2[k]
                 + pb_x[k] * ff_10[k];

        t_5[k] = pa_x[k] * dg_3[k];

        t_6[k] = f_3 * df_3[k]
                 + pa_x[k] * dg_4[k];

        t_7[k] = f_4 * df_5[k]
                 + pb_x[k] * ff_12[k];

        t_8[k] = pa_x[k] * dg_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_y, pa_z, pb_y, df_2, df_4, df_5, dg_3, \
                         dg_6, fd0_5, fd1_8, ff_17, ff_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * df_2[k]
                 + f_1 * fd0_5[k]
                 - f_2 * fd1_8[k]
                 + pb_y[k] * ff_17[k];

        t_10[k] = pa_z[k] * dg_3[k];

        t_11[k] = f_3 * df_4[k]
                  + pa_y[k] * dg_6[k];

        t_12[k] = f_4 * df_5[k]
                  + pb_y[k] * ff_22[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_y, pb_z, df_5, dg_8, fd0_8, fd1_14, \
                         ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * dg_8[k];

        t_14[k] = f_0 * df_5[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_14[k]
                  + pb_z[k] * ff_29[k];
    }
}

auto
compute_prim_fg_electron_repulsion_29(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_14 = buffer.data(dg + 14);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_9 = buffer.data(fd0 + 9);
    const auto *fd0_12 = buffer.data(fd0 + 12);
    const auto *fd0_13 = buffer.data(fd0 + 13);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_9 = buffer.data(fd1 + 9);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_13 = buffer.data(fd1 + 13);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, dg_0, fd0_1, fd0_2, fd1_1, \
                         fd1_2, ff_3, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_6[k] = pb_z[k] * ff_3[k];

        t_7[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_4[k];

        t_8[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_5[k];

        t_9[k] = pa_y[k] * dg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, pa_z, df_1, df_2, df_3, dg_0, \
                         dg_1, dg_2, dg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * df_1[k]
                  + pa_y[k] * dg_1[k];

        t_11[k] = pa_z[k] * dg_0[k];

        t_12[k] = f_5 * df_2[k]
                  + pa_z[k] * dg_2[k];

        t_13[k] = f_5 * df_3[k]
                  + pa_x[k] * dg_3[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_x, pb_x, df_4, df_5, df_7, df_8, \
                         dg_4, dg_6, dg_8, dg_9, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_6 * df_4[k]
                  + pa_x[k] * dg_4[k];

        t_15[k] = f_7 * df_5[k]
                  + pb_x[k] * ff_9[k];

        t_16[k] = pa_x[k] * dg_6[k];

        t_17[k] = f_5 * df_7[k]
                  + pa_x[k] * dg_8[k];

        t_18[k] = f_6 * df_8[k]
                  + pa_x[k] * dg_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, pb_z, df_11, dg_14, fd0_7, fd1_7, \
                         ff_11, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * df_11[k]
                  + pb_x[k] * ff_11[k];

        t_20[k] = pa_x[k] * dg_14[k];

        t_21[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_7[k]
                  + pb_x[k] * ff_12[k];

        t_22[k] = pb_z[k] * ff_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pb_x, pb_y, pb_z, df_5, fd0_8, fd0_9, \
                         fd1_8, fd1_9, ff_13, ff_14, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_x[k] * ff_13[k];

        t_24[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_9[k]
                  + pb_x[k] * ff_14[k];

        t_25[k] = pb_x[k] * ff_15[k];

        t_26[k] = f_0 * df_5[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_8[k]
                  + pb_y[k] * ff_15[k];

        t_27[k] = pb_z[k] * ff_15[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_z, pb_z, df_6, dg_6, dg_7, fd0_8, fd0_9, \
                         fd1_8, fd1_9, ff_16, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_z[k] * ff_16[k];

        t_29[k] = f_1 * fd0_9[k]
                  - f_2 * fd1_9[k]
                  + pb_z[k] * ff_17[k];

        t_30[k] = pa_z[k] * dg_6[k];

        t_31[k] = f_5 * df_6[k]
                  + pa_z[k] * dg_7[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_y, pb_y, df_9, df_10, df_11, dg_11, dg_12, \
                         dg_14, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_5 * df_9[k]
                  + pa_y[k] * dg_11[k];

        t_33[k] = f_6 * df_10[k]
                  + pa_y[k] * dg_12[k];

        t_34[k] = f_7 * df_11[k]
                  + pb_y[k] * ff_20[k];

        t_35[k] = pa_y[k] * dg_14[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_x, pb_y, fd0_12, fd0_13, fd0_14, fd1_12, \
                         fd1_13, fd1_14, ff_21, ff_22, ff_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * fd0_12[k]
                  - f_2 * fd1_12[k]
                  + pb_x[k] * ff_21[k];

        t_37[k] = pb_y[k] * ff_21[k];

        t_38[k] = f_3 * fd0_13[k]
                  - f_4 * fd1_13[k]
                  + pb_x[k] * ff_22[k];

        t_39[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_x[k] * ff_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pb_x, pb_y, pb_z, df_11, fd0_13, \
                         fd0_14, fd1_13, fd1_14, ff_24, ff_25, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_x[k] * ff_26[k];

        t_41[k] = f_1 * fd0_13[k]
                  - f_2 * fd1_13[k]
                  + pb_y[k] * ff_24[k];

        t_42[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_y[k] * ff_25[k];

        t_43[k] = pb_y[k] * ff_26[k];

        t_44[k] = f_0 * df_11[k]
                  + f_1 * fd0_14[k]
                  - f_2 * fd1_14[k]
                  + pb_z[k] * ff_26[k];
    }
}

auto
compute_prim_fg_electron_repulsion_30(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t df, const size_t dg,
                                      const size_t fd0, const size_t fd1, const size_t ff,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_8 = buffer.data(dg + 8);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_1 = buffer.data(fd0 + 1);
    const auto *fd0_2 = buffer.data(fd0 + 2);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_9 = buffer.data(fd0 + 9);
    const auto *fd0_12 = buffer.data(fd0 + 12);
    const auto *fd0_13 = buffer.data(fd0 + 13);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_1 = buffer.data(fd1 + 1);
    const auto *fd1_2 = buffer.data(fd1 + 2);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_9 = buffer.data(fd1 + 9);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_13 = buffer.data(fd1 + 13);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, df_0, fd0_0, fd1_0, ff_0, \
                         ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 + f_1 * fd0_0[k]
                 - f_2 * fd1_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = pb_y[k] * ff_0[k];

        t_2[k] = pb_z[k] * ff_0[k];

        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_y[k] * ff_1[k];

        t_4[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, dg_0, fd0_1, fd0_2, fd1_1, \
                         fd1_2, ff_3, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * fd0_1[k]
                 - f_2 * fd1_1[k]
                 + pb_y[k] * ff_3[k];

        t_6[k] = f_3 * fd0_2[k]
                 - f_4 * fd1_2[k]
                 + pb_y[k] * ff_4[k];

        t_7[k] = pb_y[k] * ff_5[k];

        t_8[k] = f_1 * fd0_2[k]
                 - f_2 * fd1_2[k]
                 + pb_z[k] * ff_5[k];

        t_9[k] = pa_y[k] * dg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pa_z, pb_x, df_4, df_5, df_8, \
                         dg_0, dg_1, dg_3, dg_4, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * dg_0[k];

        t_11[k] = f_5 * df_4[k]
                  + pa_x[k] * dg_1[k];

        t_12[k] = f_6 * df_5[k]
                  + pb_x[k] * ff_9[k];

        t_13[k] = pa_x[k] * dg_3[k];

        t_14[k] = f_5 * df_8[k]
                  + pa_x[k] * dg_4[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pb_x, df_11, dg_8, fd0_7, fd0_8, fd1_7, \
                         fd1_8, ff_11, ff_12, ff_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * df_11[k]
                  + pb_x[k] * ff_11[k];

        t_16[k] = pa_x[k] * dg_8[k];

        t_17[k] = f_1 * fd0_7[k]
                  - f_2 * fd1_7[k]
                  + pb_x[k] * ff_12[k];

        t_18[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_x[k] * ff_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pb_x, pb_y, pb_z, df_5, fd0_8, fd0_9, \
                         fd1_8, fd1_9, ff_14, ff_15, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_9[k]
                  + pb_x[k] * ff_14[k];

        t_20[k] = pb_x[k] * ff_15[k];

        t_21[k] = pb_x[k] * ff_17[k];

        t_22[k] = f_0 * df_5[k]
                  + f_1 * fd0_8[k]
                  - f_2 * fd1_8[k]
                  + pb_y[k] * ff_15[k];

        t_23[k] = pb_z[k] * ff_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pa_z, pb_z, df_10, dg_3, dg_6, fd0_8, \
                         fd0_9, fd1_8, fd1_9, ff_16, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pb_z[k] * ff_16[k];

        t_25[k] = f_1 * fd0_9[k]
                  - f_2 * fd1_9[k]
                  + pb_z[k] * ff_17[k];

        t_26[k] = pa_z[k] * dg_3[k];

        t_27[k] = f_5 * df_10[k]
                  + pa_y[k] * dg_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pb_x, pb_y, df_11, dg_8, fd0_12, \
                         fd0_13, fd1_12, fd1_13, ff_20, ff_21, ff_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_6 * df_11[k]
                  + pb_y[k] * ff_20[k];

        t_29[k] = pa_y[k] * dg_8[k];

        t_30[k] = f_1 * fd0_12[k]
                  - f_2 * fd1_12[k]
                  + pb_x[k] * ff_21[k];

        t_31[k] = f_3 * fd0_13[k]
                  - f_4 * fd1_13[k]
                  + pb_x[k] * ff_22[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, pb_x, pb_y, fd0_13, fd0_14, \
                         fd1_13, fd1_14, ff_23, ff_24, ff_25, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_x[k] * ff_23[k];

        t_33[k] = pb_x[k] * ff_24[k];

        t_34[k] = pb_x[k] * ff_26[k];

        t_35[k] = f_1 * fd0_13[k]
                  - f_2 * fd1_13[k]
                  + pb_y[k] * ff_24[k];

        t_36[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pb_y[k] * ff_25[k];

        t_37[k] = pb_y[k] * ff_26[k];
    }

#pragma omp simd aligned(t_38, pb_z, df_11, fd0_14, fd1_14, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * df_11[k]
                  + f_1 * fd0_14[k]
                  - f_2 * fd1_14[k]
                  + pb_z[k] * ff_26[k];
    }
}

}  // namespace simdt2ceri
