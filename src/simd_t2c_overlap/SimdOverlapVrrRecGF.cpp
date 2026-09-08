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


#include "SimdOverlapVrrRecGF.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_gf_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t df, const size_t fd, const size_t ff,
                          const size_t gp, const size_t gd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;

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
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_8 = buffer.data(df + 8);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_22 = buffer.data(fd + 22);
    const auto *fd_23 = buffer.data(fd + 23);
    const auto *fd_24 = buffer.data(fd + 24);
    const auto *fd_25 = buffer.data(fd + 25);
    const auto *fd_27 = buffer.data(fd + 27);
    const auto *fd_29 = buffer.data(fd + 29);

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
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_21 = buffer.data(gp + 21);
    const auto *gp_22 = buffer.data(gp + 22);
    const auto *gp_23 = buffer.data(gp + 23);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_40 = buffer.data(gd + 40);
    const auto *gd_41 = buffer.data(gd + 41);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_43 = buffer.data(gd + 43);
    const auto *gd_44 = buffer.data(gd + 44);
    const auto *gd_45 = buffer.data(gd + 45);
    const auto *gd_46 = buffer.data(gd + 46);
    const auto *gd_47 = buffer.data(gd + 47);
    const auto *gd_48 = buffer.data(gd + 48);
    const auto *gd_49 = buffer.data(gd + 49);
    const auto *gd_50 = buffer.data(gd + 50);
    const auto *gd_51 = buffer.data(gd + 51);
    const auto *gd_52 = buffer.data(gd + 52);
    const auto *gd_53 = buffer.data(gd + 53);
    const auto *gd_54 = buffer.data(gd + 54);
    const auto *gd_55 = buffer.data(gd + 55);
    const auto *gd_56 = buffer.data(gd + 56);
    const auto *gd_57 = buffer.data(gd + 57);
    const auto *gd_58 = buffer.data(gd + 58);
    const auto *gd_59 = buffer.data(gd + 59);
    const auto *gd_60 = buffer.data(gd + 60);
    const auto *gd_61 = buffer.data(gd + 61);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fd_0, fd_1, gp_0, gd_0, \
                         gd_1, gd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_0[k]
                 + f_1 * gp_0[k]
                 + pb_x[k] * gd_0[k];

        t_1[k] = pb_y[k] * gd_0[k];

        t_2[k] = pb_z[k] * gd_0[k];

        t_3[k] = f_0 * fd_1[k]
                 + pb_x[k] * gd_2[k];

        t_4[k] = pb_y[k] * gd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, t_10, pa_y, pb_x, pb_y, pb_z, fd_2, ff_0, \
                         gp_1, gp_2, gd_2, gd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * fd_2[k]
                 + pb_x[k] * gd_3[k];

        t_6[k] = f_1 * gp_1[k]
                 + pb_y[k] * gd_2[k];

        t_7[k] = pb_z[k] * gd_2[k];

        t_8[k] = pb_y[k] * gd_3[k];

        t_9[k] = f_1 * gp_2[k]
                 + pb_z[k] * gd_3[k];

        t_10[k] = pa_y[k] * ff_0[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_x, pb_y, pb_z, fd_0, fd_4, \
                         ff_2, gd_4, gd_5, gd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * fd_0[k]
                  + pb_y[k] * gd_4[k];

        t_12[k] = pb_z[k] * gd_4[k];

        t_13[k] = f_3 * fd_4[k]
                  + pb_x[k] * gd_6[k];

        t_14[k] = pb_z[k] * gd_5[k];

        t_15[k] = pa_y[k] * ff_2[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, df_1, fd_2, ff_4, \
                         ff_8, gd_6, gd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * df_1[k]
                  + pa_x[k] * ff_8[k];

        t_17[k] = pb_z[k] * gd_6[k];

        t_18[k] = f_2 * fd_2[k]
                  + pb_y[k] * gd_7[k];

        t_19[k] = pa_y[k] * ff_4[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, fd_0, ff_0, ff_1, \
                         gd_8, gd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * ff_0[k];

        t_21[k] = pb_y[k] * gd_8[k];

        t_22[k] = f_2 * fd_0[k]
                  + pb_z[k] * gd_8[k];

        t_23[k] = pa_z[k] * ff_1[k];

        t_24[k] = pb_y[k] * gd_9[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_x, pa_z, pb_x, pb_y, df_2, fd_7, \
                         ff_3, ff_12, gp_4, gd_10, gd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * fd_7[k]
                  + pb_x[k] * gd_11[k];

        t_26[k] = pa_z[k] * ff_3[k];

        t_27[k] = f_2 * gp_4[k]
                  + pb_y[k] * gd_10[k];

        t_28[k] = pb_y[k] * gd_11[k];

        t_29[k] = f_1 * df_2[k]
                  + pa_x[k] * ff_12[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pa_y, pb_x, pb_y, pb_z, df_0, fd_3, \
                         fd_9, ff_5, gd_12, gd_13, gd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * df_0[k]
                  + pa_y[k] * ff_5[k];

        t_31[k] = f_1 * fd_3[k]
                  + pb_y[k] * gd_12[k];

        t_32[k] = pb_z[k] * gd_12[k];

        t_33[k] = f_1 * fd_9[k]
                  + pb_x[k] * gd_14[k];

        t_34[k] = pb_z[k] * gd_13[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_x, pb_x, pb_y, pb_z, df_4, fd_5, \
                         fd_10, ff_16, gp_6, gd_14, gd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * fd_10[k]
                  + pb_x[k] * gd_15[k];

        t_36[k] = f_2 * df_4[k]
                  + pa_x[k] * ff_16[k];

        t_37[k] = pb_z[k] * gd_14[k];

        t_38[k] = f_1 * fd_5[k]
                  + pb_y[k] * gd_15[k];

        t_39[k] = f_1 * gp_6[k]
                  + pb_z[k] * gd_15[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, pa_y, pa_z, pb_x, fd_11, ff_6, \
                         ff_7, ff_9, ff_10, ff_11, gd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * ff_9[k];

        t_41[k] = pa_z[k] * ff_6[k];

        t_42[k] = pa_y[k] * ff_10[k];

        t_43[k] = pa_z[k] * ff_7[k];

        t_44[k] = f_1 * fd_11[k]
                  + pb_x[k] * gd_17[k];

        t_45[k] = pa_y[k] * ff_11[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pa_z, pb_y, pb_z, fd_4, fd_7, ff_8, \
                         ff_12, gd_16, gd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_z[k] * ff_8[k];

        t_47[k] = f_2 * fd_4[k]
                  + pb_z[k] * gd_16[k];

        t_48[k] = f_2 * fd_7[k]
                  + pb_y[k] * gd_18[k];

        t_49[k] = pa_y[k] * ff_12[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_z, pb_x, pb_y, pb_z, df_0, fd_6, \
                         fd_13, ff_9, gd_19, gd_20, gd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_2 * df_0[k]
                  + pa_z[k] * ff_9[k];

        t_51[k] = pb_y[k] * gd_19[k];

        t_52[k] = f_1 * fd_6[k]
                  + pb_z[k] * gd_19[k];

        t_53[k] = f_1 * fd_13[k]
                  + pb_x[k] * gd_21[k];

        t_54[k] = pb_y[k] * gd_20[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_x, pb_x, pb_y, df_8, fd_14, ff_20, \
                         gp_8, gp_9, gd_21, gd_22, gd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * fd_14[k]
                  + pb_x[k] * gd_23[k];

        t_56[k] = f_1 * gp_8[k]
                  + pb_y[k] * gd_21[k];

        t_57[k] = f_2 * gp_9[k]
                  + pb_y[k] * gd_22[k];

        t_58[k] = pb_y[k] * gd_23[k];

        t_59[k] = f_2 * df_8[k]
                  + pa_x[k] * ff_20[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pa_x, pb_x, pb_y, pb_z, fd_8, fd_15, \
                         fd_17, ff_21, gd_24, gd_25, gd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_3 * fd_15[k]
                  + pa_x[k] * ff_21[k];

        t_61[k] = f_3 * fd_8[k]
                  + pb_y[k] * gd_24[k];

        t_62[k] = pb_z[k] * gd_24[k];

        t_63[k] = f_2 * fd_17[k]
                  + pb_x[k] * gd_26[k];

        t_64[k] = pb_z[k] * gd_25[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pa_x, pb_x, pb_z, fd_18, ff_23, ff_24, \
                         ff_25, gd_26, gd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_2 * fd_18[k]
                  + pb_x[k] * gd_27[k];

        t_66[k] = pa_x[k] * ff_23[k];

        t_67[k] = pb_z[k] * gd_26[k];

        t_68[k] = pa_x[k] * ff_24[k];

        t_69[k] = pa_x[k] * ff_25[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pa_z, pb_x, pb_z, fd_8, fd_20, ff_13, \
                         ff_14, ff_15, gd_28, gd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * ff_13[k];

        t_71[k] = pa_z[k] * ff_14[k];

        t_72[k] = f_2 * fd_8[k]
                  + pb_z[k] * gd_28[k];

        t_73[k] = pa_z[k] * ff_15[k];

        t_74[k] = f_2 * fd_20[k]
                  + pb_x[k] * gd_29[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, t_80, pa_x, pa_y, pb_x, fd_21, ff_17, \
                         ff_26, ff_27, ff_28, ff_29, gd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_2 * fd_21[k]
                  + pb_x[k] * gd_30[k];

        t_76[k] = pa_x[k] * ff_26[k];

        t_77[k] = pa_x[k] * ff_27[k];

        t_78[k] = pa_x[k] * ff_28[k];

        t_79[k] = pa_x[k] * ff_29[k];

        t_80[k] = pa_y[k] * ff_17[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, pa_y, pb_x, pb_y, fd_12, fd_22, fd_23, \
                         ff_18, ff_19, gd_31, gd_32, gd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_2 * fd_12[k]
                  + pb_y[k] * gd_31[k];

        t_82[k] = pa_y[k] * ff_18[k];

        t_83[k] = f_2 * fd_22[k]
                  + pb_x[k] * gd_32[k];

        t_84[k] = f_2 * fd_23[k]
                  + pb_x[k] * gd_33[k];

        t_85[k] = pa_y[k] * ff_19[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, t_91, pa_x, pb_y, fd_25, ff_30, ff_31, \
                         ff_32, ff_33, ff_34, gd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pa_x[k] * ff_30[k];

        t_87[k] = pa_x[k] * ff_31[k];

        t_88[k] = pa_x[k] * ff_32[k];

        t_89[k] = pa_x[k] * ff_33[k];

        t_90[k] = f_3 * fd_25[k]
                  + pa_x[k] * ff_34[k];

        t_91[k] = pb_y[k] * gd_34[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, fd_12, fd_27, fd_29, gd_34, \
                         gd_35, gd_36, gd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_3 * fd_12[k]
                  + pb_z[k] * gd_34[k];

        t_93[k] = f_2 * fd_27[k]
                  + pb_x[k] * gd_36[k];

        t_94[k] = pb_y[k] * gd_35[k];

        t_95[k] = f_2 * fd_29[k]
                  + pb_x[k] * gd_37[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pa_x, pb_x, pb_y, ff_37, ff_38, ff_39, \
                         gp_12, gd_37, gd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_x[k] * ff_37[k];

        t_97[k] = pa_x[k] * ff_38[k];

        t_98[k] = pb_y[k] * gd_37[k];

        t_99[k] = pa_x[k] * ff_39[k];

        t_100[k] = f_1 * gp_12[k]
                   + pb_x[k] * gd_38[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, pb_x, pb_y, pb_z, fd_17, \
                         gp_13, gd_38, gd_39, gd_40, gd_41, gd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_2 * gp_13[k]
                   + pb_x[k] * gd_39[k];

        t_102[k] = pb_z[k] * gd_38[k];

        t_103[k] = pb_x[k] * gd_40[k];

        t_104[k] = pb_x[k] * gd_41[k];

        t_105[k] = pb_x[k] * gd_42[k];

        t_106[k] = f_0 * fd_17[k]
                   + f_1 * gp_13[k]
                   + pb_y[k] * gd_40[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, pa_z, pb_y, pb_z, fd_18, ff_21, \
                         ff_22, gp_14, gd_40, gd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pb_z[k] * gd_40[k];

        t_108[k] = f_0 * fd_18[k]
                   + pb_y[k] * gd_42[k];

        t_109[k] = f_1 * gp_14[k]
                   + pb_z[k] * gd_42[k];

        t_110[k] = pa_z[k] * ff_21[k];

        t_111[k] = pa_z[k] * ff_22[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, t_117, pa_z, pb_x, pb_z, fd_17, \
                         ff_23, gp_15, gd_43, gd_44, gd_45, gd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_2 * gp_15[k]
                   + pb_x[k] * gd_43[k];

        t_113[k] = pb_x[k] * gd_44[k];

        t_114[k] = pb_x[k] * gd_45[k];

        t_115[k] = pb_x[k] * gd_46[k];

        t_116[k] = pa_z[k] * ff_23[k];

        t_117[k] = f_2 * fd_17[k]
                   + pb_z[k] * gd_44[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pb_x, pb_y, df_5, fd_21, ff_29, \
                         gp_16, gp_17, gd_46, gd_47, gd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_3 * fd_21[k]
                   + pb_y[k] * gd_46[k];

        t_119[k] = f_1 * df_5[k]
                   + pa_y[k] * ff_29[k];

        t_120[k] = f_1 * gp_16[k]
                   + pb_x[k] * gd_47[k];

        t_121[k] = f_2 * gp_17[k]
                   + pb_x[k] * gd_48[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, pa_z, pb_x, df_4, ff_26, gp_18, \
                         gd_49, gd_50, gd_51, gd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_2 * gp_18[k]
                   + pb_x[k] * gd_49[k];

        t_123[k] = pb_x[k] * gd_50[k];

        t_124[k] = pb_x[k] * gd_51[k];

        t_125[k] = pb_x[k] * gd_52[k];

        t_126[k] = f_2 * df_4[k]
                   + pa_z[k] * ff_26[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pa_y, pb_y, pb_z, df_8, fd_19, fd_24, \
                         ff_33, ff_34, gd_50, gd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_1 * fd_19[k]
                   + pb_z[k] * gd_50[k];

        t_128[k] = f_1 * fd_24[k]
                   + pb_y[k] * gd_52[k];

        t_129[k] = f_2 * df_8[k]
                   + pa_y[k] * ff_33[k];

        t_130[k] = pa_y[k] * ff_34[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, t_136, pa_y, pb_x, fd_27, ff_36, \
                         ff_37, gp_19, gd_53, gd_54, gd_55, gd_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_2 * gp_19[k]
                   + pb_x[k] * gd_53[k];

        t_132[k] = pa_y[k] * ff_36[k];

        t_133[k] = pb_x[k] * gd_54[k];

        t_134[k] = pb_x[k] * gd_55[k];

        t_135[k] = pb_x[k] * gd_56[k];

        t_136[k] = f_3 * fd_27[k]
                   + pa_y[k] * ff_37[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, pa_y, pb_x, pb_y, pb_z, fd_22, \
                         fd_29, ff_39, gp_21, gd_54, gd_56, gd_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_3 * fd_22[k]
                   + pb_z[k] * gd_54[k];

        t_138[k] = f_2 * fd_29[k]
                   + pb_y[k] * gd_56[k];

        t_139[k] = pa_y[k] * ff_39[k];

        t_140[k] = f_1 * gp_21[k]
                   + pb_x[k] * gd_57[k];

        t_141[k] = pb_y[k] * gd_57[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, t_147, t_148, pb_x, pb_y, gp_22, \
                         gp_23, gd_58, gd_59, gd_60, gd_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_2 * gp_23[k]
                   + pb_x[k] * gd_58[k];

        t_143[k] = pb_x[k] * gd_59[k];

        t_144[k] = pb_x[k] * gd_60[k];

        t_145[k] = pb_x[k] * gd_61[k];

        t_146[k] = f_1 * gp_22[k]
                   + pb_y[k] * gd_59[k];

        t_147[k] = f_2 * gp_23[k]
                   + pb_y[k] * gd_60[k];

        t_148[k] = pb_y[k] * gd_61[k];
    }

#pragma omp simd aligned(t_149, pb_z, fd_29, gp_23, gd_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_0 * fd_29[k]
                   + f_1 * gp_23[k]
                   + pb_z[k] * gd_61[k];
    }
}

auto
compute_prim_gf_overlap_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t df, const size_t fd, const size_t ff,
                          const size_t gp, const size_t gd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_27 = buffer.data(df + 27);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_23 = buffer.data(fd + 23);
    const auto *fd_25 = buffer.data(fd + 25);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_41 = buffer.data(ff + 41);
    const auto *ff_42 = buffer.data(ff + 42);
    const auto *ff_43 = buffer.data(ff + 43);
    const auto *ff_44 = buffer.data(ff + 44);
    const auto *ff_45 = buffer.data(ff + 45);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_51 = buffer.data(ff + 51);
    const auto *ff_53 = buffer.data(ff + 53);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_40 = buffer.data(gd + 40);
    const auto *gd_41 = buffer.data(gd + 41);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_43 = buffer.data(gd + 43);
    const auto *gd_44 = buffer.data(gd + 44);
    const auto *gd_45 = buffer.data(gd + 45);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fd_0, fd_1, fd_2, gp_0, \
                         gd_0, gd_1, gd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_0[k]
                 + f_1 * gp_0[k]
                 + pb_x[k] * gd_0[k];

        t_1[k] = pb_y[k] * gd_0[k];

        t_2[k] = pb_z[k] * gd_0[k];

        t_3[k] = f_0 * fd_1[k]
                 + pb_x[k] * gd_1[k];

        t_4[k] = f_0 * fd_2[k]
                 + pb_x[k] * gd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, fd_0, ff_0, gp_1, gp_2, \
                         gd_1, gd_2, gd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * gp_1[k]
                 + pb_y[k] * gd_1[k];

        t_6[k] = pb_y[k] * gd_2[k];

        t_7[k] = f_1 * gp_2[k]
                 + pb_z[k] * gd_2[k];

        t_8[k] = pa_y[k] * ff_0[k];

        t_9[k] = f_2 * fd_0[k]
                 + pb_y[k] * gd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_y, pb_z, df_6, fd_2, fd_4, \
                         ff_6, gd_4, gd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * fd_4[k]
                  + pb_x[k] * gd_4[k];

        t_11[k] = f_1 * df_6[k]
                  + pa_x[k] * ff_6[k];

        t_12[k] = pb_z[k] * gd_4[k];

        t_13[k] = f_2 * fd_2[k]
                  + pb_y[k] * gd_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pa_z, pb_x, pb_z, fd_0, fd_8, ff_0, \
                         ff_4, gd_6, gd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_y[k] * ff_4[k];

        t_15[k] = pa_z[k] * ff_0[k];

        t_16[k] = f_2 * fd_0[k]
                  + pb_z[k] * gd_6[k];

        t_17[k] = f_3 * fd_8[k]
                  + pb_x[k] * gd_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, df_0, df_9, ff_5, ff_12, \
                         gp_3, gd_7, gd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * gp_3[k]
                  + pb_y[k] * gd_7[k];

        t_19[k] = pb_y[k] * gd_8[k];

        t_20[k] = f_1 * df_9[k]
                  + pa_x[k] * ff_12[k];

        t_21[k] = f_2 * df_0[k]
                  + pa_y[k] * ff_5[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pa_x, pb_x, pb_y, pb_z, df_12, fd_3, \
                         fd_10, ff_15, gd_9, gd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_1 * fd_3[k]
                  + pb_y[k] * gd_9[k];

        t_23[k] = pb_z[k] * gd_9[k];

        t_24[k] = f_1 * fd_10[k]
                  + pb_x[k] * gd_10[k];

        t_25[k] = f_2 * df_12[k]
                  + pa_x[k] * ff_15[k];

        t_26[k] = pb_z[k] * gd_10[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_y, pa_z, pb_y, pb_z, fd_4, fd_5, \
                         ff_6, ff_10, gp_4, gd_11, gd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * fd_5[k]
                  + pb_y[k] * gd_11[k];

        t_28[k] = f_1 * gp_4[k]
                  + pb_z[k] * gd_11[k];

        t_29[k] = pa_y[k] * ff_10[k];

        t_30[k] = pa_z[k] * ff_6[k];

        t_31[k] = f_2 * fd_4[k]
                  + pb_z[k] * gd_12[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, pa_y, pa_z, pb_y, pb_z, df_0, fd_6, \
                         fd_8, ff_9, ff_12, gd_13, gd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_2 * fd_8[k]
                  + pb_y[k] * gd_13[k];

        t_33[k] = pa_y[k] * ff_12[k];

        t_34[k] = f_2 * df_0[k]
                  + pa_z[k] * ff_9[k];

        t_35[k] = pb_y[k] * gd_14[k];

        t_36[k] = f_1 * fd_6[k]
                  + pb_z[k] * gd_14[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pa_x, pb_x, pb_y, df_27, fd_12, ff_26, \
                         gp_5, gp_6, gd_15, gd_16, gd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_1 * fd_12[k]
                  + pb_x[k] * gd_17[k];

        t_38[k] = f_1 * gp_5[k]
                  + pb_y[k] * gd_15[k];

        t_39[k] = f_2 * gp_6[k]
                  + pb_y[k] * gd_16[k];

        t_40[k] = pb_y[k] * gd_17[k];

        t_41[k] = f_2 * df_27[k]
                  + pa_x[k] * ff_26[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_x, pb_x, pb_y, pb_z, fd_9, fd_13, \
                         fd_15, ff_27, ff_31, gd_18, gd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * fd_13[k]
                  + pa_x[k] * ff_27[k];

        t_43[k] = f_3 * fd_9[k]
                  + pb_y[k] * gd_18[k];

        t_44[k] = pb_z[k] * gd_18[k];

        t_45[k] = f_2 * fd_15[k]
                  + pb_x[k] * gd_19[k];

        t_46[k] = pa_x[k] * ff_31[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, pa_x, pa_z, pb_z, fd_9, ff_13, \
                         ff_33, ff_34, ff_37, ff_38, gd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pa_x[k] * ff_33[k];

        t_48[k] = pa_x[k] * ff_34[k];

        t_49[k] = pa_z[k] * ff_13[k];

        t_50[k] = f_2 * fd_9[k]
                  + pb_z[k] * gd_20[k];

        t_51[k] = pa_x[k] * ff_37[k];

        t_52[k] = pa_x[k] * ff_38[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, t_58, pa_x, pa_y, ff_21, ff_22, ff_39, \
                         ff_41, ff_42, ff_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pa_x[k] * ff_39[k];

        t_54[k] = pa_y[k] * ff_21[k];

        t_55[k] = pa_y[k] * ff_22[k];

        t_56[k] = pa_x[k] * ff_41[k];

        t_57[k] = pa_x[k] * ff_42[k];

        t_58[k] = pa_x[k] * ff_43[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, t_63, pa_x, pb_x, pb_y, pb_z, fd_11, fd_21, \
                         fd_25, ff_45, ff_50, gd_21, gd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_3 * fd_21[k]
                  + pa_x[k] * ff_45[k];

        t_60[k] = pb_y[k] * gd_21[k];

        t_61[k] = f_3 * fd_11[k]
                  + pb_z[k] * gd_21[k];

        t_62[k] = f_2 * fd_25[k]
                  + pb_x[k] * gd_22[k];

        t_63[k] = pa_x[k] * ff_50[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, t_69, pa_x, pb_x, ff_51, ff_53, gp_7, \
                         gp_8, gd_23, gd_24, gd_25, gd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pa_x[k] * ff_51[k];

        t_65[k] = pa_x[k] * ff_53[k];

        t_66[k] = f_1 * gp_7[k]
                  + pb_x[k] * gd_23[k];

        t_67[k] = f_2 * gp_8[k]
                  + pb_x[k] * gd_24[k];

        t_68[k] = pb_x[k] * gd_25[k];

        t_69[k] = pb_x[k] * gd_26[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pb_x, pb_y, pb_z, fd_15, fd_16, gp_8, \
                         gp_9, gp_10, gd_25, gd_26, gd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_0 * fd_15[k]
                  + f_1 * gp_8[k]
                  + pb_y[k] * gd_25[k];

        t_71[k] = pb_z[k] * gd_25[k];

        t_72[k] = f_0 * fd_16[k]
                  + pb_y[k] * gd_26[k];

        t_73[k] = f_1 * gp_9[k]
                  + pb_z[k] * gd_26[k];

        t_74[k] = f_2 * gp_10[k]
                  + pb_x[k] * gd_27[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pa_z, pb_x, pb_y, pb_z, fd_15, fd_18, \
                         ff_31, gd_28, gd_29, gd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = pb_x[k] * gd_29[k];

        t_76[k] = pb_x[k] * gd_30[k];

        t_77[k] = pa_z[k] * ff_31[k];

        t_78[k] = f_2 * fd_15[k]
                  + pb_z[k] * gd_28[k];

        t_79[k] = f_3 * fd_18[k]
                  + pb_y[k] * gd_30[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pa_y, pb_x, df_19, ff_39, gp_11, gp_12, \
                         gp_13, gd_31, gd_32, gd_33, gd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_1 * df_19[k]
                  + pa_y[k] * ff_39[k];

        t_81[k] = f_1 * gp_11[k]
                  + pb_x[k] * gd_31[k];

        t_82[k] = f_2 * gp_12[k]
                  + pb_x[k] * gd_32[k];

        t_83[k] = f_2 * gp_13[k]
                  + pb_x[k] * gd_33[k];

        t_84[k] = pb_x[k] * gd_34[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pa_z, pb_x, pb_y, pb_z, df_12, fd_17, \
                         fd_20, ff_36, gd_34, gd_35, gd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pb_x[k] * gd_35[k];

        t_86[k] = pb_x[k] * gd_36[k];

        t_87[k] = f_2 * df_12[k]
                  + pa_z[k] * ff_36[k];

        t_88[k] = f_1 * fd_17[k]
                  + pb_z[k] * gd_34[k];

        t_89[k] = f_1 * fd_20[k]
                  + pb_y[k] * gd_36[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pa_y, pb_x, df_27, fd_23, ff_44, ff_50, \
                         gp_14, gd_37, gd_38, gd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_2 * df_27[k]
                  + pa_y[k] * ff_44[k];

        t_91[k] = f_2 * gp_14[k]
                  + pb_x[k] * gd_37[k];

        t_92[k] = pb_x[k] * gd_38[k];

        t_93[k] = pb_x[k] * gd_39[k];

        t_94[k] = f_3 * fd_23[k]
                  + pa_y[k] * ff_50[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_y, pb_x, pb_y, pb_z, fd_19, fd_25, ff_53, \
                         gp_16, gd_38, gd_40, gd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_3 * fd_19[k]
                  + pb_z[k] * gd_38[k];

        t_96[k] = f_2 * fd_25[k]
                  + pb_y[k] * gd_40[k];

        t_97[k] = pa_y[k] * ff_53[k];

        t_98[k] = f_1 * gp_16[k]
                  + pb_x[k] * gd_41[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, pb_x, pb_y, gp_17, gp_18, \
                         gd_42, gd_43, gd_44, gd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_2 * gp_18[k]
                  + pb_x[k] * gd_42[k];

        t_100[k] = pb_x[k] * gd_43[k];

        t_101[k] = pb_x[k] * gd_45[k];

        t_102[k] = f_1 * gp_17[k]
                   + pb_y[k] * gd_43[k];

        t_103[k] = f_2 * gp_18[k]
                   + pb_y[k] * gd_44[k];

        t_104[k] = pb_y[k] * gd_45[k];
    }

#pragma omp simd aligned(t_105, pb_z, fd_25, gp_18, gd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_0 * fd_25[k]
                   + f_1 * gp_18[k]
                   + pb_z[k] * gd_45[k];
    }
}

auto
compute_prim_gf_overlap_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t df, const size_t fd, const size_t ff,
                          const size_t gp, const size_t gd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_23 = buffer.data(df + 23);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_22 = buffer.data(fd + 22);
    const auto *fd_24 = buffer.data(fd + 24);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_43 = buffer.data(ff + 43);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, fd_0, gp_0, gp_1, \
                         gp_2, gd_0, gd_1, gd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_0[k]
                 + f_1 * gp_0[k]
                 + pb_x[k] * gd_0[k];

        t_1[k] = pb_y[k] * gd_0[k];

        t_2[k] = pb_z[k] * gd_0[k];

        t_3[k] = f_1 * gp_1[k]
                 + pb_y[k] * gd_1[k];

        t_4[k] = pb_y[k] * gd_2[k];

        t_5[k] = f_1 * gp_2[k]
                 + pb_z[k] * gd_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, pb_z, df_5, ff_0, ff_5, \
                         ff_7, gd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * ff_0[k];

        t_7[k] = f_1 * df_5[k]
                 + pa_x[k] * ff_7[k];

        t_8[k] = pb_z[k] * gd_3[k];

        t_9[k] = pa_y[k] * ff_5[k];

        t_10[k] = pa_z[k] * ff_0[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pa_x, pb_y, pb_z, df_7, fd_0, ff_13, gp_3, \
                         gd_4, gd_5, gd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * fd_0[k]
                  + pb_z[k] * gd_4[k];

        t_12[k] = f_2 * gp_3[k]
                  + pb_y[k] * gd_5[k];

        t_13[k] = pb_y[k] * gd_6[k];

        t_14[k] = f_1 * df_7[k]
                  + pa_x[k] * ff_13[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_x, pa_y, pb_x, pb_z, df_0, df_11, \
                         fd_9, ff_6, ff_17, gd_7, gd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_2 * df_0[k]
                  + pa_y[k] * ff_6[k];

        t_16[k] = pb_z[k] * gd_7[k];

        t_17[k] = f_1 * fd_9[k]
                  + pb_x[k] * gd_8[k];

        t_18[k] = f_2 * df_11[k]
                  + pa_x[k] * ff_17[k];

        t_19[k] = pb_z[k] * gd_8[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pa_z, pb_y, pb_z, df_0, ff_7, \
                         ff_10, ff_13, gp_4, gd_9, gd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * gp_4[k]
                  + pb_z[k] * gd_9[k];

        t_21[k] = pa_z[k] * ff_7[k];

        t_22[k] = pa_y[k] * ff_13[k];

        t_23[k] = f_2 * df_0[k]
                  + pa_z[k] * ff_10[k];

        t_24[k] = pb_y[k] * gd_10[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, fd_5, fd_11, gp_5, \
                         gp_6, gd_10, gd_11, gd_12, gd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * fd_5[k]
                  + pb_z[k] * gd_10[k];

        t_26[k] = f_1 * fd_11[k]
                  + pb_x[k] * gd_13[k];

        t_27[k] = f_1 * gp_5[k]
                  + pb_y[k] * gd_11[k];

        t_28[k] = f_2 * gp_6[k]
                  + pb_y[k] * gd_12[k];

        t_29[k] = pb_y[k] * gd_13[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pa_x, pa_z, pb_z, df_23, fd_12, ff_14, \
                         ff_20, ff_21, ff_25, gd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * df_23[k]
                  + pa_x[k] * ff_20[k];

        t_31[k] = f_3 * fd_12[k]
                  + pa_x[k] * ff_21[k];

        t_32[k] = pb_z[k] * gd_14[k];

        t_33[k] = pa_x[k] * ff_25[k];

        t_34[k] = pa_z[k] * ff_14[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_x, pb_x, pb_y, pb_z, fd_10, fd_20, \
                         ff_36, ff_43, gp_7, gd_15, gd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_3 * fd_20[k]
                  + pa_x[k] * ff_36[k];

        t_36[k] = pb_y[k] * gd_15[k];

        t_37[k] = f_3 * fd_10[k]
                  + pb_z[k] * gd_15[k];

        t_38[k] = pa_x[k] * ff_43[k];

        t_39[k] = f_1 * gp_7[k]
                  + pb_x[k] * gd_16[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, pb_x, pb_y, pb_z, fd_14, fd_15, \
                         gp_8, gd_17, gd_18, gd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_2 * gp_8[k]
                  + pb_x[k] * gd_17[k];

        t_41[k] = pb_x[k] * gd_18[k];

        t_42[k] = pb_x[k] * gd_19[k];

        t_43[k] = f_0 * fd_14[k]
                  + f_1 * gp_8[k]
                  + pb_y[k] * gd_18[k];

        t_44[k] = pb_z[k] * gd_18[k];

        t_45[k] = f_0 * fd_15[k]
                  + pb_y[k] * gd_19[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_z, pb_x, pb_z, ff_25, gp_9, gp_10, \
                         gd_19, gd_20, gd_22, gd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_1 * gp_9[k]
                  + pb_z[k] * gd_19[k];

        t_47[k] = f_2 * gp_10[k]
                  + pb_x[k] * gd_20[k];

        t_48[k] = pb_x[k] * gd_22[k];

        t_49[k] = pb_x[k] * gd_23[k];

        t_50[k] = pa_z[k] * ff_25[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pb_x, pb_y, pb_z, df_15, fd_14, fd_17, \
                         ff_31, gp_11, gd_21, gd_23, gd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_2 * fd_14[k]
                  + pb_z[k] * gd_21[k];

        t_52[k] = f_3 * fd_17[k]
                  + pb_y[k] * gd_23[k];

        t_53[k] = f_1 * df_15[k]
                  + pa_y[k] * ff_31[k];

        t_54[k] = f_1 * gp_11[k]
                  + pb_x[k] * gd_24[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pb_x, gp_12, gp_13, gd_25, gd_26, \
                         gd_27, gd_28, gd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_2 * gp_12[k]
                  + pb_x[k] * gd_25[k];

        t_56[k] = f_2 * gp_13[k]
                  + pb_x[k] * gd_26[k];

        t_57[k] = pb_x[k] * gd_27[k];

        t_58[k] = pb_x[k] * gd_28[k];

        t_59[k] = pb_x[k] * gd_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pa_z, pb_y, pb_z, df_11, df_23, fd_16, \
                         fd_19, ff_30, ff_35, gd_27, gd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_2 * df_11[k]
                  + pa_z[k] * ff_30[k];

        t_61[k] = f_1 * fd_16[k]
                  + pb_z[k] * gd_27[k];

        t_62[k] = f_1 * fd_19[k]
                  + pb_y[k] * gd_29[k];

        t_63[k] = f_2 * df_23[k]
                  + pa_y[k] * ff_35[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, pa_y, pb_x, pb_z, fd_18, fd_22, ff_40, \
                         gp_14, gd_30, gd_31, gd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_2 * gp_14[k]
                  + pb_x[k] * gd_30[k];

        t_65[k] = pb_x[k] * gd_31[k];

        t_66[k] = pb_x[k] * gd_32[k];

        t_67[k] = f_3 * fd_22[k]
                  + pa_y[k] * ff_40[k];

        t_68[k] = f_3 * fd_18[k]
                  + pb_z[k] * gd_31[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pa_y, pb_x, pb_y, fd_24, ff_43, gp_16, \
                         gp_18, gd_33, gd_34, gd_35, gd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_2 * fd_24[k]
                  + pb_y[k] * gd_33[k];

        t_70[k] = pa_y[k] * ff_43[k];

        t_71[k] = f_1 * gp_16[k]
                  + pb_x[k] * gd_34[k];

        t_72[k] = f_2 * gp_18[k]
                  + pb_x[k] * gd_35[k];

        t_73[k] = pb_x[k] * gd_36[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, pb_x, pb_y, pb_z, fd_24, gp_17, gp_18, \
                         gd_36, gd_37, gd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pb_x[k] * gd_38[k];

        t_75[k] = f_1 * gp_17[k]
                  + pb_y[k] * gd_36[k];

        t_76[k] = f_2 * gp_18[k]
                  + pb_y[k] * gd_37[k];

        t_77[k] = pb_y[k] * gd_38[k];

        t_78[k] = f_0 * fd_24[k]
                  + f_1 * gp_18[k]
                  + pb_z[k] * gd_38[k];
    }
}

auto
compute_prim_gf_overlap_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t df, const size_t fd, const size_t ff,
                          const size_t gp, const size_t gd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_21 = buffer.data(df + 21);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_22 = buffer.data(fd + 22);
    const auto *fd_24 = buffer.data(fd + 24);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_42 = buffer.data(ff + 42);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, fd_0, gp_0, gp_1, \
                         gp_2, gd_0, gd_1, gd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_0[k]
                 + f_1 * gp_0[k]
                 + pb_x[k] * gd_0[k];

        t_1[k] = pb_y[k] * gd_0[k];

        t_2[k] = pb_z[k] * gd_0[k];

        t_3[k] = f_1 * gp_1[k]
                 + pb_y[k] * gd_1[k];

        t_4[k] = pb_y[k] * gd_2[k];

        t_5[k] = f_1 * gp_2[k]
                 + pb_z[k] * gd_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, pa_z, pb_z, df_5, fd_0, ff_0, ff_7, gd_3, \
                         gd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * df_5[k]
                 + pa_x[k] * ff_7[k];

        t_7[k] = pb_z[k] * gd_3[k];

        t_8[k] = pa_z[k] * ff_0[k];

        t_9[k] = f_2 * fd_0[k]
                 + pb_z[k] * gd_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, pb_y, df_0, df_6, ff_6, ff_12, \
                         gp_3, gd_5, gd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * gp_3[k]
                  + pb_y[k] * gd_5[k];

        t_11[k] = pb_y[k] * gd_6[k];

        t_12[k] = f_1 * df_6[k]
                  + pa_x[k] * ff_12[k];

        t_13[k] = f_2 * df_0[k]
                  + pa_y[k] * ff_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_x, pb_x, pb_z, df_10, fd_9, ff_16, \
                         gp_4, gd_7, gd_8, gd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pb_z[k] * gd_7[k];

        t_15[k] = f_1 * fd_9[k]
                  + pb_x[k] * gd_8[k];

        t_16[k] = f_2 * df_10[k]
                  + pa_x[k] * ff_16[k];

        t_17[k] = pb_z[k] * gd_8[k];

        t_18[k] = f_1 * gp_4[k]
                  + pb_z[k] * gd_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_z, pb_x, pb_y, pb_z, df_0, fd_5, fd_11, \
                         ff_9, gd_10, gd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * df_0[k]
                  + pa_z[k] * ff_9[k];

        t_20[k] = pb_y[k] * gd_10[k];

        t_21[k] = f_1 * fd_5[k]
                  + pb_z[k] * gd_10[k];

        t_22[k] = f_1 * fd_11[k]
                  + pb_x[k] * gd_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pa_x, pb_y, df_21, fd_12, ff_19, ff_20, \
                         gp_5, gp_6, gd_11, gd_12, gd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * gp_5[k]
                  + pb_y[k] * gd_11[k];

        t_24[k] = f_2 * gp_6[k]
                  + pb_y[k] * gd_12[k];

        t_25[k] = pb_y[k] * gd_13[k];

        t_26[k] = f_2 * df_21[k]
                  + pa_x[k] * ff_19[k];

        t_27[k] = f_3 * fd_12[k]
                  + pa_x[k] * ff_20[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_x, pb_x, pb_y, pb_z, fd_10, fd_20, \
                         ff_35, gp_7, gd_14, gd_15, gd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_z[k] * gd_14[k];

        t_29[k] = f_3 * fd_20[k]
                  + pa_x[k] * ff_35[k];

        t_30[k] = pb_y[k] * gd_15[k];

        t_31[k] = f_3 * fd_10[k]
                  + pb_z[k] * gd_15[k];

        t_32[k] = f_1 * gp_7[k]
                  + pb_x[k] * gd_16[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, fd_14, fd_15, \
                         gp_8, gd_17, gd_18, gd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_2 * gp_8[k]
                  + pb_x[k] * gd_17[k];

        t_34[k] = pb_x[k] * gd_18[k];

        t_35[k] = pb_x[k] * gd_19[k];

        t_36[k] = f_0 * fd_14[k]
                  + f_1 * gp_8[k]
                  + pb_y[k] * gd_18[k];

        t_37[k] = pb_z[k] * gd_18[k];

        t_38[k] = f_0 * fd_15[k]
                  + pb_y[k] * gd_19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_z, pb_x, pb_z, ff_24, gp_9, gp_10, \
                         gd_19, gd_20, gd_22, gd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_1 * gp_9[k]
                  + pb_z[k] * gd_19[k];

        t_40[k] = f_2 * gp_10[k]
                  + pb_x[k] * gd_20[k];

        t_41[k] = pb_x[k] * gd_22[k];

        t_42[k] = pb_x[k] * gd_23[k];

        t_43[k] = pa_z[k] * ff_24[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_y, pb_x, pb_y, pb_z, df_13, fd_14, fd_17, \
                         ff_30, gp_11, gd_21, gd_23, gd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_2 * fd_14[k]
                  + pb_z[k] * gd_21[k];

        t_45[k] = f_3 * fd_17[k]
                  + pb_y[k] * gd_23[k];

        t_46[k] = f_1 * df_13[k]
                  + pa_y[k] * ff_30[k];

        t_47[k] = f_1 * gp_11[k]
                  + pb_x[k] * gd_24[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pb_x, gp_12, gp_13, gd_25, gd_26, \
                         gd_27, gd_28, gd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_2 * gp_12[k]
                  + pb_x[k] * gd_25[k];

        t_49[k] = f_2 * gp_13[k]
                  + pb_x[k] * gd_26[k];

        t_50[k] = pb_x[k] * gd_27[k];

        t_51[k] = pb_x[k] * gd_28[k];

        t_52[k] = pb_x[k] * gd_29[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_y, pa_z, pb_y, pb_z, df_10, df_21, fd_16, \
                         fd_19, ff_29, ff_34, gd_27, gd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_2 * df_10[k]
                  + pa_z[k] * ff_29[k];

        t_54[k] = f_1 * fd_16[k]
                  + pb_z[k] * gd_27[k];

        t_55[k] = f_1 * fd_19[k]
                  + pb_y[k] * gd_29[k];

        t_56[k] = f_2 * df_21[k]
                  + pa_y[k] * ff_34[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pa_y, pb_x, pb_z, fd_18, fd_22, ff_39, \
                         gp_14, gd_30, gd_31, gd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_2 * gp_14[k]
                  + pb_x[k] * gd_30[k];

        t_58[k] = pb_x[k] * gd_31[k];

        t_59[k] = pb_x[k] * gd_32[k];

        t_60[k] = f_3 * fd_22[k]
                  + pa_y[k] * ff_39[k];

        t_61[k] = f_3 * fd_18[k]
                  + pb_z[k] * gd_31[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, pa_y, pb_x, pb_y, fd_24, ff_42, gp_16, \
                         gp_18, gd_33, gd_34, gd_35, gd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_2 * fd_24[k]
                  + pb_y[k] * gd_33[k];

        t_63[k] = pa_y[k] * ff_42[k];

        t_64[k] = f_1 * gp_16[k]
                  + pb_x[k] * gd_34[k];

        t_65[k] = f_2 * gp_18[k]
                  + pb_x[k] * gd_35[k];

        t_66[k] = pb_x[k] * gd_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, fd_24, gp_17, gp_18, \
                         gd_36, gd_37, gd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pb_x[k] * gd_38[k];

        t_68[k] = f_1 * gp_17[k]
                  + pb_y[k] * gd_36[k];

        t_69[k] = f_2 * gp_18[k]
                  + pb_y[k] * gd_37[k];

        t_70[k] = pb_y[k] * gd_38[k];

        t_71[k] = f_0 * fd_24[k]
                  + f_1 * gp_18[k]
                  + pb_z[k] * gd_38[k];
    }
}

auto
compute_prim_gf_overlap_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t df, const size_t fd, const size_t ff,
                          const size_t gp, const size_t gd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;

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

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_22 = buffer.data(fd + 22);

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

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_40 = buffer.data(gd + 40);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_43 = buffer.data(gd + 43);
    const auto *gd_44 = buffer.data(gd + 44);
    const auto *gd_45 = buffer.data(gd + 45);
    const auto *gd_46 = buffer.data(gd + 46);
    const auto *gd_47 = buffer.data(gd + 47);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, fd_0, fd_1, fd_2, gp_0, gp_1, gd_0, \
                         gd_1, gd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_0[k]
                 + f_1 * gp_0[k]
                 + pb_x[k] * gd_0[k];

        t_1[k] = f_0 * fd_1[k]
                 + pb_x[k] * gd_1[k];

        t_2[k] = f_0 * fd_2[k]
                 + pb_x[k] * gd_2[k];

        t_3[k] = f_1 * gp_1[k]
                 + pb_y[k] * gd_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pb_x, pb_y, pb_z, fd_0, fd_4, ff_0, gp_2, \
                         gd_2, gd_3, gd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * gp_2[k]
                 + pb_z[k] * gd_2[k];

        t_5[k] = pa_y[k] * ff_0[k];

        t_6[k] = f_2 * fd_0[k]
                 + pb_y[k] * gd_3[k];

        t_7[k] = f_3 * fd_4[k]
                 + pb_x[k] * gd_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_z, pb_x, pb_z, df_1, fd_0, fd_6, ff_0, \
                         ff_2, gd_6, gd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * df_1[k]
                 + pa_x[k] * ff_2[k];

        t_9[k] = pa_z[k] * ff_0[k];

        t_10[k] = f_2 * fd_0[k]
                  + pb_z[k] * gd_6[k];

        t_11[k] = f_3 * fd_6[k]
                  + pb_x[k] * gd_7[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pa_y, pb_x, pb_y, df_0, df_2, fd_3, \
                         fd_8, ff_1, ff_4, gd_8, gd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_1 * df_2[k]
                  + pa_x[k] * ff_4[k];

        t_13[k] = f_2 * df_0[k]
                  + pa_y[k] * ff_1[k];

        t_14[k] = f_1 * fd_3[k]
                  + pb_y[k] * gd_8[k];

        t_15[k] = f_1 * fd_8[k]
                  + pb_x[k] * gd_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_z, pb_x, pb_z, df_0, df_3, fd_5, \
                         fd_10, ff_3, ff_5, gd_14, gd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * df_3[k]
                  + pa_x[k] * ff_5[k];

        t_17[k] = f_2 * df_0[k]
                  + pa_z[k] * ff_3[k];

        t_18[k] = f_1 * fd_5[k]
                  + pb_z[k] * gd_14[k];

        t_19[k] = f_1 * fd_10[k]
                  + pb_x[k] * gd_16[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_y, df_6, fd_7, fd_11, fd_12, \
                         ff_6, ff_7, gd_17, gd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_2 * df_6[k]
                  + pa_x[k] * ff_6[k];

        t_21[k] = f_3 * fd_11[k]
                  + pa_x[k] * ff_7[k];

        t_22[k] = f_3 * fd_7[k]
                  + pb_y[k] * gd_17[k];

        t_23[k] = f_2 * fd_12[k]
                  + pb_x[k] * gd_18[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_x, pb_z, fd_9, fd_20, ff_8, ff_10, \
                         ff_11, ff_13, gd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_x[k] * ff_8[k];

        t_25[k] = pa_x[k] * ff_10[k];

        t_26[k] = pa_x[k] * ff_11[k];

        t_27[k] = f_3 * fd_20[k]
                  + pa_x[k] * ff_13[k];

        t_28[k] = f_3 * fd_9[k]
                  + pb_z[k] * gd_26[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_x, pb_x, fd_22, ff_15, gp_12, gp_13, \
                         gd_28, gd_29, gd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_2 * fd_22[k]
                  + pb_x[k] * gd_28[k];

        t_30[k] = pa_x[k] * ff_15[k];

        t_31[k] = f_1 * gp_12[k]
                  + pb_x[k] * gd_29[k];

        t_32[k] = f_2 * gp_13[k]
                  + pb_x[k] * gd_30[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_z, pb_y, pb_z, fd_12, fd_13, ff_8, \
                         gp_13, gp_14, gd_31, gd_32, gd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * fd_12[k]
                  + f_1 * gp_13[k]
                  + pb_y[k] * gd_31[k];

        t_34[k] = f_0 * fd_13[k]
                  + pb_y[k] * gd_32[k];

        t_35[k] = f_1 * gp_14[k]
                  + pb_z[k] * gd_32[k];

        t_36[k] = pa_z[k] * ff_8[k];

        t_37[k] = f_2 * fd_12[k]
                  + pb_z[k] * gd_33[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pa_z, pb_x, pb_y, df_3, df_4, fd_16, \
                         ff_9, ff_10, gp_16, gd_35, gd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * fd_16[k]
                  + pb_y[k] * gd_35[k];

        t_39[k] = f_1 * df_4[k]
                  + pa_y[k] * ff_10[k];

        t_40[k] = f_1 * gp_16[k]
                  + pb_x[k] * gd_36[k];

        t_41[k] = f_2 * df_3[k]
                  + pa_z[k] * ff_9[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pb_y, pb_z, df_6, fd_14, fd_19, fd_21, \
                         ff_12, ff_14, gd_37, gd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_1 * fd_14[k]
                  + pb_z[k] * gd_37[k];

        t_43[k] = f_1 * fd_19[k]
                  + pb_y[k] * gd_39[k];

        t_44[k] = f_2 * df_6[k]
                  + pa_y[k] * ff_12[k];

        t_45[k] = f_3 * fd_21[k]
                  + pa_y[k] * ff_14[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_y, pb_x, pb_y, pb_z, fd_17, fd_22, \
                         ff_15, gp_19, gd_40, gd_42, gd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_3 * fd_17[k]
                  + pb_z[k] * gd_40[k];

        t_47[k] = f_2 * fd_22[k]
                  + pb_y[k] * gd_42[k];

        t_48[k] = pa_y[k] * ff_15[k];

        t_49[k] = f_1 * gp_19[k]
                  + pb_x[k] * gd_43[k];

        t_50[k] = pb_y[k] * gd_43[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_x, pb_y, pb_z, fd_22, gp_20, gp_21, gd_44, \
                         gd_45, gd_46, gd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_2 * gp_21[k]
                  + pb_x[k] * gd_44[k];

        t_52[k] = f_1 * gp_20[k]
                  + pb_y[k] * gd_45[k];

        t_53[k] = f_2 * gp_21[k]
                  + pb_y[k] * gd_46[k];

        t_54[k] = f_0 * fd_22[k]
                  + f_1 * gp_21[k]
                  + pb_z[k] * gd_47[k];
    }
}

auto
compute_prim_gf_overlap_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t df, const size_t fd, const size_t ff,
                          const size_t gp, const size_t gd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_11 = buffer.data(df + 11);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_22 = buffer.data(fd + 22);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_19 = buffer.data(ff + 19);
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
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_38 = buffer.data(ff + 38);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_40 = buffer.data(gd + 40);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fd_0, gp_0, gp_1, gp_2, \
                         gd_0, gd_1, gd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_0[k]
                 + f_1 * gp_0[k]
                 + pb_x[k] * gd_0[k];

        t_1[k] = pb_y[k] * gd_0[k];

        t_2[k] = pb_z[k] * gd_0[k];

        t_3[k] = f_1 * gp_1[k]
                 + pb_y[k] * gd_1[k];

        t_4[k] = f_1 * gp_2[k]
                 + pb_z[k] * gd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pa_z, pb_y, df_2, fd_2, ff_0, \
                         ff_4, ff_6, gd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_y[k] * ff_0[k];

        t_6[k] = f_1 * df_2[k]
                 + pa_x[k] * ff_6[k];

        t_7[k] = f_2 * fd_2[k]
                 + pb_y[k] * gd_5[k];

        t_8[k] = pa_y[k] * ff_4[k];

        t_9[k] = pa_z[k] * ff_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, pb_y, pb_z, df_0, df_3, fd_0, \
                         ff_5, ff_9, gp_3, gd_6, gd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * fd_0[k]
                  + pb_z[k] * gd_6[k];

        t_11[k] = f_2 * gp_3[k]
                  + pb_y[k] * gd_7[k];

        t_12[k] = f_1 * df_3[k]
                  + pa_x[k] * ff_9[k];

        t_13[k] = f_2 * df_0[k]
                  + pa_y[k] * ff_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pb_x, pb_y, pb_z, df_5, fd_5, fd_9, \
                         ff_12, gp_4, gd_10, gd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * fd_9[k]
                  + pb_x[k] * gd_10[k];

        t_15[k] = f_2 * df_5[k]
                  + pa_x[k] * ff_12[k];

        t_16[k] = f_1 * fd_5[k]
                  + pb_y[k] * gd_11[k];

        t_17[k] = f_1 * gp_4[k]
                  + pb_z[k] * gd_11[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_y, pa_z, pb_y, pb_z, fd_4, fd_7, \
                         ff_6, ff_8, ff_9, gd_12, gd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pa_y[k] * ff_8[k];

        t_19[k] = pa_z[k] * ff_6[k];

        t_20[k] = f_2 * fd_4[k]
                  + pb_z[k] * gd_12[k];

        t_21[k] = f_2 * fd_7[k]
                  + pb_y[k] * gd_13[k];

        t_22[k] = pa_y[k] * ff_9[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pb_x, pb_y, pb_z, df_0, fd_6, fd_11, \
                         ff_7, gd_14, gd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_2 * df_0[k]
                  + pa_z[k] * ff_7[k];

        t_24[k] = pb_y[k] * gd_14[k];

        t_25[k] = f_1 * fd_6[k]
                  + pb_z[k] * gd_14[k];

        t_26[k] = f_1 * fd_11[k]
                  + pb_x[k] * gd_17[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pb_y, df_11, fd_12, ff_16, ff_17, gp_5, \
                         gp_6, gd_15, gd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * gp_5[k]
                  + pb_y[k] * gd_15[k];

        t_28[k] = f_2 * gp_6[k]
                  + pb_y[k] * gd_16[k];

        t_29[k] = f_2 * df_11[k]
                  + pa_x[k] * ff_16[k];

        t_30[k] = f_3 * fd_12[k]
                  + pa_x[k] * ff_17[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_x, pa_z, pb_x, fd_13, ff_10, ff_19, \
                         ff_21, ff_22, gd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_2 * fd_13[k]
                  + pb_x[k] * gd_19[k];

        t_32[k] = pa_x[k] * ff_19[k];

        t_33[k] = pa_x[k] * ff_21[k];

        t_34[k] = pa_x[k] * ff_22[k];

        t_35[k] = pa_z[k] * ff_10[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, pa_x, pa_y, pb_z, fd_8, ff_13, \
                         ff_14, ff_24, ff_25, ff_26, gd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_2 * fd_8[k]
                  + pb_z[k] * gd_20[k];

        t_37[k] = pa_x[k] * ff_24[k];

        t_38[k] = pa_x[k] * ff_25[k];

        t_39[k] = pa_x[k] * ff_26[k];

        t_40[k] = pa_y[k] * ff_13[k];

        t_41[k] = pa_y[k] * ff_14[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_x, pb_z, fd_10, fd_19, ff_27, ff_28, \
                         ff_29, ff_31, gd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_x[k] * ff_27[k];

        t_43[k] = pa_x[k] * ff_28[k];

        t_44[k] = pa_x[k] * ff_29[k];

        t_45[k] = f_3 * fd_19[k]
                  + pa_x[k] * ff_31[k];

        t_46[k] = f_3 * fd_10[k]
                  + pb_z[k] * gd_23[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pa_x, pb_x, fd_22, ff_35, ff_36, ff_38, \
                         gp_7, gd_24, gd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_2 * fd_22[k]
                  + pb_x[k] * gd_24[k];

        t_48[k] = pa_x[k] * ff_35[k];

        t_49[k] = pa_x[k] * ff_36[k];

        t_50[k] = pa_x[k] * ff_38[k];

        t_51[k] = f_1 * gp_7[k]
                  + pb_x[k] * gd_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, t_57, pb_x, pb_y, pb_z, fd_13, fd_14, \
                         gp_8, gd_26, gd_27, gd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_2 * gp_8[k]
                  + pb_x[k] * gd_26[k];

        t_53[k] = pb_x[k] * gd_27[k];

        t_54[k] = pb_x[k] * gd_28[k];

        t_55[k] = f_0 * fd_13[k]
                  + f_1 * gp_8[k]
                  + pb_y[k] * gd_27[k];

        t_56[k] = pb_z[k] * gd_27[k];

        t_57[k] = f_0 * fd_14[k]
                  + pb_y[k] * gd_28[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_z, pb_x, pb_y, pb_z, fd_13, fd_16, \
                         ff_19, gp_9, gd_28, gd_29, gd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * gp_9[k]
                  + pb_z[k] * gd_28[k];

        t_59[k] = pb_x[k] * gd_30[k];

        t_60[k] = pa_z[k] * ff_19[k];

        t_61[k] = f_2 * fd_13[k]
                  + pb_z[k] * gd_29[k];

        t_62[k] = f_3 * fd_16[k]
                  + pb_y[k] * gd_30[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pa_z, pb_x, df_5, df_7, ff_23, \
                         ff_26, gp_11, gd_31, gd_32, gd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_1 * df_7[k]
                  + pa_y[k] * ff_26[k];

        t_64[k] = f_1 * gp_11[k]
                  + pb_x[k] * gd_31[k];

        t_65[k] = pb_x[k] * gd_32[k];

        t_66[k] = pb_x[k] * gd_33[k];

        t_67[k] = f_2 * df_5[k]
                  + pa_z[k] * ff_23[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_y, pb_x, pb_y, pb_z, df_11, fd_15, fd_18, \
                         ff_30, gd_32, gd_33, gd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_1 * fd_15[k]
                  + pb_z[k] * gd_32[k];

        t_69[k] = f_1 * fd_18[k]
                  + pb_y[k] * gd_33[k];

        t_70[k] = f_2 * df_11[k]
                  + pa_y[k] * ff_30[k];

        t_71[k] = pb_x[k] * gd_34[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_y, pb_y, pb_z, fd_17, fd_20, fd_22, ff_35, \
                         ff_38, gd_34, gd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_3 * fd_20[k]
                  + pa_y[k] * ff_35[k];

        t_73[k] = f_3 * fd_17[k]
                  + pb_z[k] * gd_34[k];

        t_74[k] = f_2 * fd_22[k]
                  + pb_y[k] * gd_35[k];

        t_75[k] = pa_y[k] * ff_38[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, t_81, pb_x, pb_y, gp_14, gp_15, gp_16, \
                         gd_36, gd_37, gd_38, gd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_1 * gp_14[k]
                  + pb_x[k] * gd_36[k];

        t_77[k] = pb_y[k] * gd_36[k];

        t_78[k] = f_2 * gp_16[k]
                  + pb_x[k] * gd_37[k];

        t_79[k] = pb_x[k] * gd_38[k];

        t_80[k] = pb_x[k] * gd_40[k];

        t_81[k] = f_1 * gp_15[k]
                  + pb_y[k] * gd_38[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pb_y, pb_z, fd_22, gp_16, gd_39, \
                         gd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_2 * gp_16[k]
                  + pb_y[k] * gd_39[k];

        t_83[k] = pb_y[k] * gd_40[k];

        t_84[k] = f_0 * fd_22[k]
                  + f_1 * gp_16[k]
                  + pb_z[k] * gd_40[k];
    }
}

auto
compute_prim_gf_overlap_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t df, const size_t fd, const size_t ff,
                          const size_t gp, const size_t gd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_20 = buffer.data(df + 20);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_19 = buffer.data(fd + 19);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_36 = buffer.data(ff + 36);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, fd_0, gp_0, gp_1, \
                         gp_2, gd_0, gd_1, gd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_0[k]
                 + f_1 * gp_0[k]
                 + pb_x[k] * gd_0[k];

        t_1[k] = pb_y[k] * gd_0[k];

        t_2[k] = pb_z[k] * gd_0[k];

        t_3[k] = f_1 * gp_1[k]
                 + pb_y[k] * gd_1[k];

        t_4[k] = pb_y[k] * gd_2[k];

        t_5[k] = f_1 * gp_2[k]
                 + pb_z[k] * gd_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, pb_z, df_4, ff_0, ff_4, \
                         ff_6, gd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * ff_0[k];

        t_7[k] = f_1 * df_4[k]
                 + pa_x[k] * ff_6[k];

        t_8[k] = pb_z[k] * gd_4[k];

        t_9[k] = pa_y[k] * ff_4[k];

        t_10[k] = pa_z[k] * ff_0[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pa_x, pa_y, pb_y, df_0, df_6, ff_5, ff_9, \
                         gp_3, gd_6, gd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * gp_3[k]
                  + pb_y[k] * gd_6[k];

        t_12[k] = pb_y[k] * gd_7[k];

        t_13[k] = f_1 * df_6[k]
                  + pa_x[k] * ff_9[k];

        t_14[k] = f_2 * df_0[k]
                  + pa_y[k] * ff_5[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_x, pb_x, pb_z, df_9, fd_7, ff_12, \
                         gp_4, gd_8, gd_9, gd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_z[k] * gd_8[k];

        t_16[k] = f_1 * fd_7[k]
                  + pb_x[k] * gd_9[k];

        t_17[k] = f_2 * df_9[k]
                  + pa_x[k] * ff_12[k];

        t_18[k] = pb_z[k] * gd_9[k];

        t_19[k] = f_1 * gp_4[k]
                  + pb_z[k] * gd_10[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pa_z, pb_x, pb_y, df_0, fd_8, \
                         ff_6, ff_8, ff_9, gd_11, gd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * ff_6[k];

        t_21[k] = pa_y[k] * ff_9[k];

        t_22[k] = f_2 * df_0[k]
                  + pa_z[k] * ff_8[k];

        t_23[k] = pb_y[k] * gd_11[k];

        t_24[k] = f_1 * fd_8[k]
                  + pb_x[k] * gd_14[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_x, pb_y, df_20, fd_9, ff_15, ff_16, \
                         gp_5, gp_6, gd_12, gd_13, gd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * gp_5[k]
                  + pb_y[k] * gd_12[k];

        t_26[k] = f_2 * gp_6[k]
                  + pb_y[k] * gd_13[k];

        t_27[k] = pb_y[k] * gd_14[k];

        t_28[k] = f_2 * df_20[k]
                  + pa_x[k] * ff_15[k];

        t_29[k] = f_3 * fd_9[k]
                  + pa_x[k] * ff_16[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pa_x, pa_z, pb_x, pb_z, fd_10, ff_10, \
                         ff_19, ff_24, gd_15, gd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_z[k] * gd_15[k];

        t_31[k] = f_2 * fd_10[k]
                  + pb_x[k] * gd_16[k];

        t_32[k] = pa_x[k] * ff_19[k];

        t_33[k] = pa_z[k] * ff_10[k];

        t_34[k] = pa_x[k] * ff_24[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_x, pb_x, fd_16, fd_19, ff_26, ff_29, \
                         ff_36, gp_7, gd_18, gd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_x[k] * ff_26[k];

        t_36[k] = f_3 * fd_16[k]
                  + pa_x[k] * ff_29[k];

        t_37[k] = f_2 * fd_19[k]
                  + pb_x[k] * gd_18[k];

        t_38[k] = pa_x[k] * ff_36[k];

        t_39[k] = f_1 * gp_7[k]
                  + pb_x[k] * gd_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, pb_x, pb_y, pb_z, fd_10, fd_11, \
                         gp_8, gd_20, gd_21, gd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_2 * gp_8[k]
                  + pb_x[k] * gd_20[k];

        t_41[k] = pb_x[k] * gd_21[k];

        t_42[k] = pb_x[k] * gd_22[k];

        t_43[k] = f_0 * fd_10[k]
                  + f_1 * gp_8[k]
                  + pb_y[k] * gd_21[k];

        t_44[k] = pb_z[k] * gd_21[k];

        t_45[k] = f_0 * fd_11[k]
                  + pb_y[k] * gd_22[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pa_z, pb_x, pb_z, df_13, ff_19, ff_24, \
                         gp_9, gd_22, gd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_1 * gp_9[k]
                  + pb_z[k] * gd_22[k];

        t_47[k] = pb_x[k] * gd_24[k];

        t_48[k] = pa_z[k] * ff_19[k];

        t_49[k] = f_1 * df_13[k]
                  + pa_y[k] * ff_24[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_z, pb_x, pb_y, df_9, fd_15, ff_23, \
                         gp_11, gd_25, gd_26, gd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_1 * gp_11[k]
                  + pb_x[k] * gd_25[k];

        t_51[k] = pb_x[k] * gd_26[k];

        t_52[k] = pb_x[k] * gd_27[k];

        t_53[k] = f_2 * df_9[k]
                  + pa_z[k] * ff_23[k];

        t_54[k] = f_1 * fd_15[k]
                  + pb_y[k] * gd_27[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_y, pb_x, pb_y, df_20, fd_17, fd_19, \
                         ff_28, ff_33, ff_36, gd_28, gd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_2 * df_20[k]
                  + pa_y[k] * ff_28[k];

        t_56[k] = pb_x[k] * gd_28[k];

        t_57[k] = f_3 * fd_17[k]
                  + pa_y[k] * ff_33[k];

        t_58[k] = f_2 * fd_19[k]
                  + pb_y[k] * gd_29[k];

        t_59[k] = pa_y[k] * ff_36[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pb_x, pb_y, gp_14, gp_15, gp_16, \
                         gd_30, gd_31, gd_32, gd_33, gd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * gp_14[k]
                  + pb_x[k] * gd_30[k];

        t_61[k] = f_2 * gp_16[k]
                  + pb_x[k] * gd_31[k];

        t_62[k] = pb_x[k] * gd_32[k];

        t_63[k] = pb_x[k] * gd_34[k];

        t_64[k] = f_1 * gp_15[k]
                  + pb_y[k] * gd_32[k];

        t_65[k] = f_2 * gp_16[k]
                  + pb_y[k] * gd_33[k];
    }

#pragma omp simd aligned(t_66, t_67, pb_y, pb_z, fd_19, gp_16, gd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pb_y[k] * gd_34[k];

        t_67[k] = f_0 * fd_19[k]
                  + f_1 * gp_16[k]
                  + pb_z[k] * gd_34[k];
    }
}

auto
compute_prim_gf_overlap_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t df, const size_t fd, const size_t ff,
                          const size_t gp, const size_t gd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;

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
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_18 = buffer.data(df + 18);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_19 = buffer.data(fd + 19);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_33 = buffer.data(ff + 33);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, fd_0, gp_0, gp_1, \
                         gp_2, gd_0, gd_1, gd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_0[k]
                 + f_1 * gp_0[k]
                 + pb_x[k] * gd_0[k];

        t_1[k] = pb_y[k] * gd_0[k];

        t_2[k] = pb_z[k] * gd_0[k];

        t_3[k] = f_1 * gp_1[k]
                 + pb_y[k] * gd_1[k];

        t_4[k] = pb_y[k] * gd_2[k];

        t_5[k] = f_1 * gp_2[k]
                 + pb_z[k] * gd_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, pb_y, pb_z, df_4, ff_0, \
                         ff_6, gp_3, gd_4, gd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * ff_0[k];

        t_7[k] = f_1 * df_4[k]
                 + pa_x[k] * ff_6[k];

        t_8[k] = pb_z[k] * gd_4[k];

        t_9[k] = pa_z[k] * ff_0[k];

        t_10[k] = f_2 * gp_3[k]
                  + pb_y[k] * gd_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pa_x, pa_y, pb_y, pb_z, df_0, df_5, ff_5, \
                         ff_8, gd_7, gd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_y[k] * gd_7[k];

        t_12[k] = f_1 * df_5[k]
                  + pa_x[k] * ff_8[k];

        t_13[k] = f_2 * df_0[k]
                  + pa_y[k] * ff_5[k];

        t_14[k] = pb_z[k] * gd_8[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pb_x, pb_z, df_8, fd_7, ff_10, gp_4, \
                         gd_9, gd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * fd_7[k]
                  + pb_x[k] * gd_9[k];

        t_16[k] = f_2 * df_8[k]
                  + pa_x[k] * ff_10[k];

        t_17[k] = pb_z[k] * gd_9[k];

        t_18[k] = f_1 * gp_4[k]
                  + pb_z[k] * gd_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_z, pb_x, pb_y, df_0, fd_8, ff_7, gp_5, \
                         gd_11, gd_12, gd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * df_0[k]
                  + pa_z[k] * ff_7[k];

        t_20[k] = pb_y[k] * gd_11[k];

        t_21[k] = f_1 * fd_8[k]
                  + pb_x[k] * gd_14[k];

        t_22[k] = f_1 * gp_5[k]
                  + pb_y[k] * gd_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pa_x, pb_y, pb_z, df_18, fd_9, ff_12, \
                         ff_13, gp_6, gd_13, gd_14, gd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_2 * gp_6[k]
                  + pb_y[k] * gd_13[k];

        t_24[k] = pb_y[k] * gd_14[k];

        t_25[k] = f_2 * df_18[k]
                  + pa_x[k] * ff_12[k];

        t_26[k] = f_3 * fd_9[k]
                  + pa_x[k] * ff_13[k];

        t_27[k] = pb_z[k] * gd_15[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_x, pb_x, fd_10, fd_16, fd_19, ff_16, \
                         ff_26, ff_33, gd_16, gd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_2 * fd_10[k]
                  + pb_x[k] * gd_16[k];

        t_29[k] = pa_x[k] * ff_16[k];

        t_30[k] = f_3 * fd_16[k]
                  + pa_x[k] * ff_26[k];

        t_31[k] = f_2 * fd_19[k]
                  + pb_x[k] * gd_18[k];

        t_32[k] = pa_x[k] * ff_33[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, fd_10, gp_7, \
                         gp_8, gd_19, gd_20, gd_21, gd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_1 * gp_7[k]
                  + pb_x[k] * gd_19[k];

        t_34[k] = f_2 * gp_8[k]
                  + pb_x[k] * gd_20[k];

        t_35[k] = pb_x[k] * gd_21[k];

        t_36[k] = pb_x[k] * gd_22[k];

        t_37[k] = f_0 * fd_10[k]
                  + f_1 * gp_8[k]
                  + pb_y[k] * gd_21[k];

        t_38[k] = pb_z[k] * gd_21[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pb_x, pb_y, pb_z, fd_11, ff_16, gp_9, \
                         gd_22, gd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * fd_11[k]
                  + pb_y[k] * gd_22[k];

        t_40[k] = f_1 * gp_9[k]
                  + pb_z[k] * gd_22[k];

        t_41[k] = pb_x[k] * gd_24[k];

        t_42[k] = pa_z[k] * ff_16[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pa_y, pa_z, pb_x, df_8, df_11, ff_20, \
                         ff_21, gp_11, gd_25, gd_26, gd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * df_11[k]
                  + pa_y[k] * ff_21[k];

        t_44[k] = f_1 * gp_11[k]
                  + pb_x[k] * gd_25[k];

        t_45[k] = pb_x[k] * gd_26[k];

        t_46[k] = pb_x[k] * gd_27[k];

        t_47[k] = f_2 * df_8[k]
                  + pa_z[k] * ff_20[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pb_x, pb_y, df_18, fd_15, fd_17, ff_25, \
                         ff_30, gd_27, gd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * fd_15[k]
                  + pb_y[k] * gd_27[k];

        t_49[k] = f_2 * df_18[k]
                  + pa_y[k] * ff_25[k];

        t_50[k] = pb_x[k] * gd_28[k];

        t_51[k] = f_3 * fd_17[k]
                  + pa_y[k] * ff_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pa_y, pb_x, pb_y, fd_19, ff_33, gp_14, \
                         gp_16, gd_29, gd_30, gd_31, gd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_2 * fd_19[k]
                  + pb_y[k] * gd_29[k];

        t_53[k] = pa_y[k] * ff_33[k];

        t_54[k] = f_1 * gp_14[k]
                  + pb_x[k] * gd_30[k];

        t_55[k] = f_2 * gp_16[k]
                  + pb_x[k] * gd_31[k];

        t_56[k] = pb_x[k] * gd_32[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pb_x, pb_y, pb_z, fd_19, gp_15, gp_16, \
                         gd_32, gd_33, gd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_x[k] * gd_34[k];

        t_58[k] = f_1 * gp_15[k]
                  + pb_y[k] * gd_32[k];

        t_59[k] = f_2 * gp_16[k]
                  + pb_y[k] * gd_33[k];

        t_60[k] = pb_y[k] * gd_34[k];

        t_61[k] = f_0 * fd_19[k]
                  + f_1 * gp_16[k]
                  + pb_z[k] * gd_34[k];
    }
}

auto
compute_prim_gf_overlap_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t df, const size_t fd, const size_t ff,
                          const size_t gp, const size_t gd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;

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
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_6 = buffer.data(df + 6);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);

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

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, df_1, fd_0, ff_0, ff_2, \
                         gp_0, gd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_0[k]
                 + f_1 * gp_0[k]
                 + pb_x[k] * gd_0[k];

        t_1[k] = pa_y[k] * ff_0[k];

        t_2[k] = f_1 * df_1[k]
                 + pa_x[k] * ff_2[k];

        t_3[k] = pa_z[k] * ff_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, df_0, df_2, df_3, ff_1, ff_3, \
                         ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * df_2[k]
                 + pa_x[k] * ff_4[k];

        t_5[k] = f_2 * df_0[k]
                 + pa_y[k] * ff_1[k];

        t_6[k] = f_2 * df_3[k]
                 + pa_x[k] * ff_5[k];

        t_7[k] = f_2 * df_0[k]
                 + pa_z[k] * ff_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, pa_x, pb_x, df_6, ff_6, ff_7, ff_9, \
                         ff_10, ff_13, gp_6, gd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * df_6[k]
                 + pa_x[k] * ff_6[k];

        t_9[k] = pa_x[k] * ff_7[k];

        t_10[k] = pa_x[k] * ff_9[k];

        t_11[k] = pa_x[k] * ff_10[k];

        t_12[k] = pa_x[k] * ff_13[k];

        t_13[k] = f_1 * gp_6[k]
                  + pb_x[k] * gd_17[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pa_z, pb_y, df_3, df_4, fd_8, ff_7, \
                         ff_8, ff_9, gp_7, gd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * fd_8[k]
                  + f_1 * gp_7[k]
                  + pb_y[k] * gd_18[k];

        t_15[k] = pa_z[k] * ff_7[k];

        t_16[k] = f_1 * df_4[k]
                  + pa_y[k] * ff_9[k];

        t_17[k] = f_2 * df_3[k]
                  + pa_z[k] * ff_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_y, pb_x, df_6, fd_14, ff_11, ff_12, ff_13, \
                         gp_13, gd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * df_6[k]
                  + pa_y[k] * ff_11[k];

        t_19[k] = f_3 * fd_14[k]
                  + pa_y[k] * ff_12[k];

        t_20[k] = pa_y[k] * ff_13[k];

        t_21[k] = f_1 * gp_13[k]
                  + pb_x[k] * gd_30[k];
    }

#pragma omp simd aligned(t_22, t_23, pb_y, pb_z, fd_15, gp_14, gp_15, gd_31, \
                         gd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_1 * gp_14[k]
                  + pb_y[k] * gd_31[k];

        t_23[k] = f_0 * fd_15[k]
                  + f_1 * gp_15[k]
                  + pb_z[k] * gd_32[k];
    }
}

auto
compute_prim_gf_overlap_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t df, const size_t fd, const size_t ff,
                          const size_t gp, const size_t gd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;

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
    const auto *df_9 = buffer.data(df + 9);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_21 = buffer.data(ff + 21);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fd_0, gp_0, gp_1, gp_2, \
                         gd_0, gd_1, gd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_0[k]
                 + f_1 * gp_0[k]
                 + pb_x[k] * gd_0[k];

        t_1[k] = pb_y[k] * gd_0[k];

        t_2[k] = pb_z[k] * gd_0[k];

        t_3[k] = f_1 * gp_1[k]
                 + pb_y[k] * gd_1[k];

        t_4[k] = f_1 * gp_2[k]
                 + pb_z[k] * gd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pa_z, pb_z, df_2, df_3, fd_0, \
                         ff_0, ff_3, ff_5, gd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_y[k] * ff_0[k];

        t_6[k] = f_1 * df_2[k]
                 + pa_x[k] * ff_3[k];

        t_7[k] = pa_z[k] * ff_0[k];

        t_8[k] = f_2 * fd_0[k]
                 + pb_z[k] * gd_6[k];

        t_9[k] = f_1 * df_3[k]
                 + pa_x[k] * ff_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, pa_z, pb_x, df_0, df_4, fd_7, \
                         ff_2, ff_4, ff_7, gd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * df_0[k]
                  + pa_y[k] * ff_2[k];

        t_11[k] = f_1 * fd_7[k]
                  + pb_x[k] * gd_9[k];

        t_12[k] = f_2 * df_4[k]
                  + pa_x[k] * ff_7[k];

        t_13[k] = f_2 * df_0[k]
                  + pa_z[k] * ff_4[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pb_x, pb_z, df_9, fd_4, fd_9, fd_10, \
                         ff_9, ff_10, gd_13, gd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * fd_4[k]
                  + pb_z[k] * gd_13[k];

        t_15[k] = f_1 * fd_9[k]
                  + pb_x[k] * gd_14[k];

        t_16[k] = f_2 * df_9[k]
                  + pa_x[k] * ff_9[k];

        t_17[k] = f_3 * fd_10[k]
                  + pa_x[k] * ff_10[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_x, pb_x, fd_11, fd_17, ff_11, ff_14, \
                         ff_15, ff_18, gd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * fd_11[k]
                  + pb_x[k] * gd_16[k];

        t_19[k] = pa_x[k] * ff_11[k];

        t_20[k] = pa_x[k] * ff_14[k];

        t_21[k] = pa_x[k] * ff_15[k];

        t_22[k] = f_3 * fd_17[k]
                  + pa_x[k] * ff_18[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pa_x, pb_x, pb_z, fd_8, fd_19, ff_21, \
                         gp_4, gd_20, gd_21, gd_22, gd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * fd_8[k]
                  + pb_z[k] * gd_20[k];

        t_24[k] = f_2 * fd_19[k]
                  + pb_x[k] * gd_21[k];

        t_25[k] = pa_x[k] * ff_21[k];

        t_26[k] = f_1 * gp_4[k]
                  + pb_x[k] * gd_22[k];

        t_27[k] = pb_x[k] * gd_23[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_z, pb_y, pb_z, fd_11, fd_12, ff_11, \
                         gp_5, gp_6, gd_23, gd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * fd_11[k]
                  + f_1 * gp_5[k]
                  + pb_y[k] * gd_23[k];

        t_29[k] = pb_z[k] * gd_23[k];

        t_30[k] = f_0 * fd_12[k]
                  + pb_y[k] * gd_24[k];

        t_31[k] = f_1 * gp_6[k]
                  + pb_z[k] * gd_24[k];

        t_32[k] = pa_z[k] * ff_11[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_y, pb_x, pb_y, pb_z, df_6, fd_11, fd_14, \
                         ff_14, gp_7, gd_25, gd_26, gd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_2 * fd_11[k]
                  + pb_z[k] * gd_25[k];

        t_34[k] = f_3 * fd_14[k]
                  + pb_y[k] * gd_26[k];

        t_35[k] = f_1 * df_6[k]
                  + pa_y[k] * ff_14[k];

        t_36[k] = f_1 * gp_7[k]
                  + pb_x[k] * gd_27[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_y, pa_z, pb_y, pb_z, df_4, df_9, fd_13, \
                         fd_16, ff_13, ff_17, gd_28, gd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_2 * df_4[k]
                  + pa_z[k] * ff_13[k];

        t_38[k] = f_1 * fd_13[k]
                  + pb_z[k] * gd_28[k];

        t_39[k] = f_1 * fd_16[k]
                  + pb_y[k] * gd_29[k];

        t_40[k] = f_2 * df_9[k]
                  + pa_y[k] * ff_17[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_y, pb_y, pb_z, fd_15, fd_18, fd_19, ff_19, \
                         ff_21, gd_30, gd_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_3 * fd_18[k]
                  + pa_y[k] * ff_19[k];

        t_42[k] = f_3 * fd_15[k]
                  + pb_z[k] * gd_30[k];

        t_43[k] = f_2 * fd_19[k]
                  + pb_y[k] * gd_31[k];

        t_44[k] = pa_y[k] * ff_21[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, t_51, pb_x, pb_y, gp_8, gp_9, \
                         gp_10, gd_32, gd_33, gd_34, gd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_1 * gp_8[k]
                  + pb_x[k] * gd_32[k];

        t_46[k] = pb_y[k] * gd_32[k];

        t_47[k] = pb_x[k] * gd_33[k];

        t_48[k] = pb_x[k] * gd_35[k];

        t_49[k] = f_1 * gp_9[k]
                  + pb_y[k] * gd_33[k];

        t_50[k] = f_2 * gp_10[k]
                  + pb_y[k] * gd_34[k];

        t_51[k] = pb_y[k] * gd_35[k];
    }

#pragma omp simd aligned(t_52, pb_z, fd_19, gp_10, gd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_0 * fd_19[k]
                  + f_1 * gp_10[k]
                  + pb_z[k] * gd_35[k];
    }
}

auto
compute_prim_gf_overlap_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                           const size_t pb, const size_t df, const size_t fd, const size_t ff,
                           const size_t gp, const size_t gd, const size_t ncols,
                           const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_13 = buffer.data(df + 13);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_30 = buffer.data(ff + 30);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fd_0, gp_0, gp_1, gp_2, \
                         gd_0, gd_1, gd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_0[k]
                 + f_1 * gp_0[k]
                 + pb_x[k] * gd_0[k];

        t_1[k] = pb_y[k] * gd_0[k];

        t_2[k] = pb_z[k] * gd_0[k];

        t_3[k] = f_1 * gp_1[k]
                 + pb_y[k] * gd_1[k];

        t_4[k] = f_1 * gp_2[k]
                 + pb_z[k] * gd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pa_z, df_2, df_3, ff_0, ff_3, \
                         ff_5, ff_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_y[k] * ff_0[k];

        t_6[k] = f_1 * df_2[k]
                 + pa_x[k] * ff_5[k];

        t_7[k] = pa_y[k] * ff_3[k];

        t_8[k] = pa_z[k] * ff_0[k];

        t_9[k] = f_1 * df_3[k]
                 + pa_x[k] * ff_7[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, pb_x, pb_z, df_0, df_6, fd_6, \
                         ff_4, ff_10, gp_3, gd_8, gd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * df_0[k]
                  + pa_y[k] * ff_4[k];

        t_11[k] = f_1 * fd_6[k]
                  + pb_x[k] * gd_8[k];

        t_12[k] = f_2 * df_6[k]
                  + pa_x[k] * ff_10[k];

        t_13[k] = f_1 * gp_3[k]
                  + pb_z[k] * gd_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_y, pa_z, pb_x, pb_y, df_0, fd_7, \
                         ff_5, ff_6, ff_7, gd_10, gd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_z[k] * ff_5[k];

        t_15[k] = pa_y[k] * ff_7[k];

        t_16[k] = f_2 * df_0[k]
                  + pa_z[k] * ff_6[k];

        t_17[k] = pb_y[k] * gd_10[k];

        t_18[k] = f_1 * fd_7[k]
                  + pb_x[k] * gd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_x, pa_z, pb_x, df_13, fd_8, fd_9, \
                         ff_8, ff_13, ff_14, ff_16, gd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * df_13[k]
                  + pa_x[k] * ff_13[k];

        t_20[k] = f_3 * fd_8[k]
                  + pa_x[k] * ff_14[k];

        t_21[k] = f_2 * fd_9[k]
                  + pb_x[k] * gd_12[k];

        t_22[k] = pa_x[k] * ff_16[k];

        t_23[k] = pa_z[k] * ff_8[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_x, pb_x, fd_14, fd_16, ff_20, ff_21, \
                         ff_24, ff_30, gd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_x[k] * ff_20[k];

        t_25[k] = pa_x[k] * ff_21[k];

        t_26[k] = f_3 * fd_14[k]
                  + pa_x[k] * ff_24[k];

        t_27[k] = f_2 * fd_16[k]
                  + pb_x[k] * gd_13[k];

        t_28[k] = pa_x[k] * ff_30[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, pb_x, pb_y, pb_z, fd_9, gp_4, \
                         gp_5, gp_6, gd_14, gd_15, gd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * gp_4[k]
                  + pb_x[k] * gd_14[k];

        t_30[k] = pb_x[k] * gd_15[k];

        t_31[k] = pb_x[k] * gd_16[k];

        t_32[k] = f_0 * fd_9[k]
                  + f_1 * gp_5[k]
                  + pb_y[k] * gd_15[k];

        t_33[k] = pb_z[k] * gd_15[k];

        t_34[k] = f_1 * gp_6[k]
                  + pb_z[k] * gd_16[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_y, pa_z, pb_x, df_8, ff_16, ff_20, \
                         gp_7, gd_18, gd_19, gd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pb_x[k] * gd_18[k];

        t_36[k] = pa_z[k] * ff_16[k];

        t_37[k] = f_1 * df_8[k]
                  + pa_y[k] * ff_20[k];

        t_38[k] = f_1 * gp_7[k]
                  + pb_x[k] * gd_19[k];

        t_39[k] = pb_x[k] * gd_20[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pa_y, pa_z, pb_x, pb_y, df_6, df_13, \
                         fd_13, ff_19, ff_23, gd_21, gd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_x[k] * gd_21[k];

        t_41[k] = f_2 * df_6[k]
                  + pa_z[k] * ff_19[k];

        t_42[k] = f_1 * fd_13[k]
                  + pb_y[k] * gd_21[k];

        t_43[k] = f_2 * df_13[k]
                  + pa_y[k] * ff_23[k];

        t_44[k] = pb_x[k] * gd_22[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_y, pb_x, pb_y, fd_15, fd_16, ff_28, \
                         ff_30, gp_8, gd_23, gd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_3 * fd_15[k]
                  + pa_y[k] * ff_28[k];

        t_46[k] = f_2 * fd_16[k]
                  + pb_y[k] * gd_23[k];

        t_47[k] = pa_y[k] * ff_30[k];

        t_48[k] = f_1 * gp_8[k]
                  + pb_x[k] * gd_24[k];

        t_49[k] = pb_y[k] * gd_24[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, t_55, pb_x, pb_y, pb_z, fd_16, gp_9, \
                         gp_10, gd_25, gd_26, gd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pb_x[k] * gd_25[k];

        t_51[k] = pb_x[k] * gd_27[k];

        t_52[k] = f_1 * gp_9[k]
                  + pb_y[k] * gd_25[k];

        t_53[k] = f_2 * gp_10[k]
                  + pb_y[k] * gd_26[k];

        t_54[k] = pb_y[k] * gd_27[k];

        t_55[k] = f_0 * fd_16[k]
                  + f_1 * gp_10[k]
                  + pb_z[k] * gd_27[k];
    }
}

auto
compute_prim_gf_overlap_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                           const size_t pb, const size_t df, const size_t fd, const size_t ff,
                           const size_t gp, const size_t gd, const size_t ncols,
                           const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;

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
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_13 = buffer.data(df + 13);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_27 = buffer.data(ff + 27);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fd_0, gp_0, gp_1, gp_2, \
                         gd_0, gd_1, gd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_0[k]
                 + f_1 * gp_0[k]
                 + pb_x[k] * gd_0[k];

        t_1[k] = pb_y[k] * gd_0[k];

        t_2[k] = pb_z[k] * gd_0[k];

        t_3[k] = f_1 * gp_1[k]
                 + pb_y[k] * gd_1[k];

        t_4[k] = f_1 * gp_2[k]
                 + pb_z[k] * gd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pa_z, df_0, df_2, df_3, ff_0, \
                         ff_4, ff_5, ff_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_y[k] * ff_0[k];

        t_6[k] = f_1 * df_2[k]
                 + pa_x[k] * ff_5[k];

        t_7[k] = pa_z[k] * ff_0[k];

        t_8[k] = f_1 * df_3[k]
                 + pa_x[k] * ff_7[k];

        t_9[k] = f_2 * df_0[k]
                 + pa_y[k] * ff_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_z, pb_x, pb_z, df_0, df_6, fd_6, \
                         ff_6, ff_9, gp_3, gd_8, gd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * fd_6[k]
                  + pb_x[k] * gd_8[k];

        t_11[k] = f_2 * df_6[k]
                  + pa_x[k] * ff_9[k];

        t_12[k] = f_1 * gp_3[k]
                  + pb_z[k] * gd_9[k];

        t_13[k] = f_2 * df_0[k]
                  + pa_z[k] * ff_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_x, pb_x, pb_y, df_13, fd_7, fd_9, \
                         ff_11, ff_14, gd_10, gd_11, gd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pb_y[k] * gd_10[k];

        t_15[k] = f_1 * fd_7[k]
                  + pb_x[k] * gd_11[k];

        t_16[k] = f_2 * df_13[k]
                  + pa_x[k] * ff_11[k];

        t_17[k] = f_2 * fd_9[k]
                  + pb_x[k] * gd_12[k];

        t_18[k] = pa_x[k] * ff_14[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_x, pb_x, fd_16, ff_27, gp_4, gd_13, \
                         gd_14, gd_15, gd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * fd_16[k]
                  + pb_x[k] * gd_13[k];

        t_20[k] = pa_x[k] * ff_27[k];

        t_21[k] = f_1 * gp_4[k]
                  + pb_x[k] * gd_14[k];

        t_22[k] = pb_x[k] * gd_15[k];

        t_23[k] = pb_x[k] * gd_16[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_z, pb_x, pb_y, pb_z, fd_9, ff_14, \
                         gp_5, gp_6, gd_15, gd_16, gd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * fd_9[k]
                  + f_1 * gp_5[k]
                  + pb_y[k] * gd_15[k];

        t_25[k] = pb_z[k] * gd_15[k];

        t_26[k] = f_1 * gp_6[k]
                  + pb_z[k] * gd_16[k];

        t_27[k] = pb_x[k] * gd_18[k];

        t_28[k] = pa_z[k] * ff_14[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pa_y, pa_z, pb_x, df_6, df_8, ff_17, \
                         ff_18, gp_7, gd_19, gd_20, gd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * df_8[k]
                  + pa_y[k] * ff_18[k];

        t_30[k] = f_1 * gp_7[k]
                  + pb_x[k] * gd_19[k];

        t_31[k] = pb_x[k] * gd_20[k];

        t_32[k] = pb_x[k] * gd_21[k];

        t_33[k] = f_2 * df_6[k]
                  + pa_z[k] * ff_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pb_x, pb_y, df_13, fd_13, fd_15, ff_20, \
                         ff_25, gd_21, gd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_1 * fd_13[k]
                  + pb_y[k] * gd_21[k];

        t_35[k] = f_2 * df_13[k]
                  + pa_y[k] * ff_20[k];

        t_36[k] = pb_x[k] * gd_22[k];

        t_37[k] = f_3 * fd_15[k]
                  + pa_y[k] * ff_25[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, pa_y, pb_x, pb_y, fd_16, ff_27, \
                         gp_8, gd_23, gd_24, gd_25, gd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_2 * fd_16[k]
                  + pb_y[k] * gd_23[k];

        t_39[k] = pa_y[k] * ff_27[k];

        t_40[k] = f_1 * gp_8[k]
                  + pb_x[k] * gd_24[k];

        t_41[k] = pb_y[k] * gd_24[k];

        t_42[k] = pb_x[k] * gd_25[k];

        t_43[k] = pb_x[k] * gd_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_y, pb_z, fd_16, gp_9, gp_10, gd_25, gd_26, \
                         gd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_1 * gp_9[k]
                  + pb_y[k] * gd_25[k];

        t_45[k] = f_2 * gp_10[k]
                  + pb_y[k] * gd_26[k];

        t_46[k] = pb_y[k] * gd_27[k];

        t_47[k] = f_0 * fd_16[k]
                  + f_1 * gp_10[k]
                  + pb_z[k] * gd_27[k];
    }
}

}  // namespace simdovl
