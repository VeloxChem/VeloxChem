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


#include "SimdElectronRepulsionVrrRecGF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_gf_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t df0, const size_t df1,
                                     const size_t fd, const size_t ff, const size_t gp0,
                                     const size_t gp1, const size_t gd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_4 = 1.5 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / p;

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

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_36 = buffer.data(df0 + 36);
    const auto *df0_59 = buffer.data(df0 + 59);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_36 = buffer.data(df1 + 36);
    const auto *df1_59 = buffer.data(df1 + 59);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_23 = buffer.data(fd + 23);
    const auto *fd_28 = buffer.data(fd + 28);
    const auto *fd_30 = buffer.data(fd + 30);
    const auto *fd_33 = buffer.data(fd + 33);
    const auto *fd_35 = buffer.data(fd + 35);
    const auto *fd_36 = buffer.data(fd + 36);
    const auto *fd_39 = buffer.data(fd + 39);
    const auto *fd_41 = buffer.data(fd + 41);
    const auto *fd_42 = buffer.data(fd + 42);
    const auto *fd_45 = buffer.data(fd + 45);
    const auto *fd_46 = buffer.data(fd + 46);
    const auto *fd_47 = buffer.data(fd + 47);
    const auto *fd_48 = buffer.data(fd + 48);
    const auto *fd_51 = buffer.data(fd + 51);
    const auto *fd_52 = buffer.data(fd + 52);
    const auto *fd_53 = buffer.data(fd + 53);
    const auto *fd_54 = buffer.data(fd + 54);
    const auto *fd_57 = buffer.data(fd + 57);
    const auto *fd_59 = buffer.data(fd + 59);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_55 = buffer.data(ff + 55);
    const auto *ff_59 = buffer.data(ff + 59);
    const auto *ff_60 = buffer.data(ff + 60);
    const auto *ff_61 = buffer.data(ff + 61);
    const auto *ff_66 = buffer.data(ff + 66);
    const auto *ff_68 = buffer.data(ff + 68);
    const auto *ff_69 = buffer.data(ff + 69);
    const auto *ff_76 = buffer.data(ff + 76);
    const auto *ff_77 = buffer.data(ff + 77);
    const auto *ff_78 = buffer.data(ff + 78);
    const auto *ff_79 = buffer.data(ff + 79);
    const auto *ff_86 = buffer.data(ff + 86);
    const auto *ff_87 = buffer.data(ff + 87);
    const auto *ff_88 = buffer.data(ff + 88);
    const auto *ff_89 = buffer.data(ff + 89);
    const auto *ff_90 = buffer.data(ff + 90);
    const auto *ff_92 = buffer.data(ff + 92);
    const auto *ff_96 = buffer.data(ff + 96);
    const auto *ff_97 = buffer.data(ff + 97);
    const auto *ff_99 = buffer.data(ff + 99);

    const auto *gp0_0 = buffer.data(gp0 + 0);
    const auto *gp0_1 = buffer.data(gp0 + 1);
    const auto *gp0_2 = buffer.data(gp0 + 2);
    const auto *gp0_11 = buffer.data(gp0 + 11);
    const auto *gp0_16 = buffer.data(gp0 + 16);
    const auto *gp0_30 = buffer.data(gp0 + 30);
    const auto *gp0_31 = buffer.data(gp0 + 31);
    const auto *gp0_32 = buffer.data(gp0 + 32);
    const auto *gp0_36 = buffer.data(gp0 + 36);
    const auto *gp0_42 = buffer.data(gp0 + 42);
    const auto *gp0_43 = buffer.data(gp0 + 43);
    const auto *gp0_44 = buffer.data(gp0 + 44);

    const auto *gp1_0 = buffer.data(gp1 + 0);
    const auto *gp1_1 = buffer.data(gp1 + 1);
    const auto *gp1_2 = buffer.data(gp1 + 2);
    const auto *gp1_11 = buffer.data(gp1 + 11);
    const auto *gp1_16 = buffer.data(gp1 + 16);
    const auto *gp1_30 = buffer.data(gp1 + 30);
    const auto *gp1_31 = buffer.data(gp1 + 31);
    const auto *gp1_32 = buffer.data(gp1 + 32);
    const auto *gp1_36 = buffer.data(gp1 + 36);
    const auto *gp1_42 = buffer.data(gp1 + 42);
    const auto *gp1_43 = buffer.data(gp1 + 43);
    const auto *gp1_44 = buffer.data(gp1 + 44);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_41 = buffer.data(gd + 41);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_46 = buffer.data(gd + 46);
    const auto *gd_47 = buffer.data(gd + 47);
    const auto *gd_48 = buffer.data(gd + 48);
    const auto *gd_51 = buffer.data(gd + 51);
    const auto *gd_52 = buffer.data(gd + 52);
    const auto *gd_54 = buffer.data(gd + 54);
    const auto *gd_56 = buffer.data(gd + 56);
    const auto *gd_57 = buffer.data(gd + 57);
    const auto *gd_59 = buffer.data(gd + 59);
    const auto *gd_60 = buffer.data(gd + 60);
    const auto *gd_63 = buffer.data(gd + 63);
    const auto *gd_64 = buffer.data(gd + 64);
    const auto *gd_65 = buffer.data(gd + 65);
    const auto *gd_66 = buffer.data(gd + 66);
    const auto *gd_69 = buffer.data(gd + 69);
    const auto *gd_70 = buffer.data(gd + 70);
    const auto *gd_71 = buffer.data(gd + 71);
    const auto *gd_72 = buffer.data(gd + 72);
    const auto *gd_75 = buffer.data(gd + 75);
    const auto *gd_76 = buffer.data(gd + 76);
    const auto *gd_77 = buffer.data(gd + 77);
    const auto *gd_78 = buffer.data(gd + 78);
    const auto *gd_81 = buffer.data(gd + 81);
    const auto *gd_82 = buffer.data(gd + 82);
    const auto *gd_83 = buffer.data(gd + 83);
    const auto *gd_84 = buffer.data(gd + 84);
    const auto *gd_87 = buffer.data(gd + 87);
    const auto *gd_88 = buffer.data(gd + 88);
    const auto *gd_89 = buffer.data(gd + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fd_0, fd_3, gp0_0, gp1_0, \
                         gd_0, gd_2, gd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_0[k]
                 + f_1 * gp0_0[k]
                 - f_2 * gp1_0[k]
                 + pb_x[k] * gd_0[k];

        t_1[k] = pb_y[k] * gd_0[k];

        t_2[k] = pb_z[k] * gd_0[k];

        t_3[k] = f_0 * fd_3[k]
                 + pb_x[k] * gd_3[k];

        t_4[k] = pb_y[k] * gd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, fd_5, gp0_1, gp0_2, gp1_1, \
                         gp1_2, gd_3, gd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * fd_5[k]
                 + pb_x[k] * gd_5[k];

        t_6[k] = f_1 * gp0_1[k]
                 - f_2 * gp1_1[k]
                 + pb_y[k] * gd_3[k];

        t_7[k] = pb_z[k] * gd_3[k];

        t_8[k] = pb_y[k] * gd_5[k];

        t_9[k] = f_1 * gp0_2[k]
                 - f_2 * gp1_2[k]
                 + pb_z[k] * gd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, pb_y, pb_z, fd_0, fd_9, \
                         ff_0, gd_6, gd_7, gd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * ff_0[k];

        t_11[k] = f_3 * fd_0[k]
                  + pb_y[k] * gd_6[k];

        t_12[k] = pb_z[k] * gd_6[k];

        t_13[k] = f_4 * fd_9[k]
                  + pb_x[k] * gd_9[k];

        t_14[k] = pb_z[k] * gd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, fd_3, fd_5, ff_5, \
                         ff_6, ff_9, gd_9, gd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * ff_5[k];

        t_16[k] = f_4 * fd_3[k]
                  + pa_y[k] * ff_6[k];

        t_17[k] = pb_z[k] * gd_9[k];

        t_18[k] = f_3 * fd_5[k]
                  + pb_y[k] * gd_11[k];

        t_19[k] = pa_y[k] * ff_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, fd_0, ff_0, ff_3, \
                         gd_12, gd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * ff_0[k];

        t_21[k] = pb_y[k] * gd_12[k];

        t_22[k] = f_3 * fd_0[k]
                  + pb_z[k] * gd_12[k];

        t_23[k] = pa_z[k] * ff_3[k];

        t_24[k] = pb_y[k] * gd_14[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, pb_z, fd_3, fd_5, \
                         fd_17, ff_6, ff_9, gd_15, gd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * fd_17[k]
                  + pb_x[k] * gd_17[k];

        t_26[k] = pa_z[k] * ff_6[k];

        t_27[k] = f_3 * fd_3[k]
                  + pb_z[k] * gd_15[k];

        t_28[k] = pb_y[k] * gd_17[k];

        t_29[k] = f_4 * fd_5[k]
                  + pa_z[k] * ff_9[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_x, pb_y, pb_z, df0_0, df1_0, fd_6, \
                         fd_21, ff_10, gd_18, gd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * df0_0[k]
                  - f_6 * df1_0[k]
                  + pa_y[k] * ff_10[k];

        t_31[k] = f_7 * fd_6[k]
                  + pb_y[k] * gd_18[k];

        t_32[k] = pb_z[k] * gd_18[k];

        t_33[k] = f_7 * fd_21[k]
                  + pb_x[k] * gd_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, df0_36, df1_36, fd_23, \
                         ff_36, gd_19, gd_21, gd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_z[k] * gd_19[k];

        t_35[k] = f_7 * fd_23[k]
                  + pb_x[k] * gd_23[k];

        t_36[k] = f_5 * df0_36[k]
                  - f_6 * df1_36[k]
                  + pa_x[k] * ff_36[k];

        t_37[k] = pb_z[k] * gd_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, pb_z, fd_11, ff_11, \
                         ff_20, ff_22, gp0_11, gp1_11, gd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_7 * fd_11[k]
                  + pb_y[k] * gd_23[k];

        t_39[k] = f_1 * gp0_11[k]
                  - f_2 * gp1_11[k]
                  + pb_z[k] * gd_23[k];

        t_40[k] = pa_y[k] * ff_20[k];

        t_41[k] = pa_z[k] * ff_11[k];

        t_42[k] = pa_y[k] * ff_22[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pa_y, pa_z, pb_x, pb_z, fd_9, fd_28, \
                         ff_13, ff_16, ff_25, gd_27, gd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pa_z[k] * ff_13[k];

        t_44[k] = f_7 * fd_28[k]
                  + pb_x[k] * gd_28[k];

        t_45[k] = pa_y[k] * ff_25[k];

        t_46[k] = pa_z[k] * ff_16[k];

        t_47[k] = f_3 * fd_9[k]
                  + pb_z[k] * gd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pa_z, pb_y, df0_0, df1_0, fd_17, ff_20, \
                         ff_29, gd_29, gd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * fd_17[k]
                  + pb_y[k] * gd_29[k];

        t_49[k] = pa_y[k] * ff_29[k];

        t_50[k] = f_5 * df0_0[k]
                  - f_6 * df1_0[k]
                  + pa_z[k] * ff_20[k];

        t_51[k] = pb_y[k] * gd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_x, pb_y, pb_z, fd_12, fd_33, fd_35, gd_30, \
                         gd_32, gd_33, gd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_7 * fd_12[k]
                  + pb_z[k] * gd_30[k];

        t_53[k] = f_7 * fd_33[k]
                  + pb_x[k] * gd_33[k];

        t_54[k] = pb_y[k] * gd_32[k];

        t_55[k] = f_7 * fd_35[k]
                  + pb_x[k] * gd_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_y, pb_z, df0_59, df1_59, fd_15, \
                         ff_59, gp0_16, gp1_16, gd_33, gd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * gp0_16[k]
                  - f_2 * gp1_16[k]
                  + pb_y[k] * gd_33[k];

        t_57[k] = f_7 * fd_15[k]
                  + pb_z[k] * gd_33[k];

        t_58[k] = pb_y[k] * gd_35[k];

        t_59[k] = f_5 * df0_59[k]
                  - f_6 * df1_59[k]
                  + pa_x[k] * ff_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pa_x, pb_x, pb_y, pb_z, fd_18, fd_36, \
                         fd_39, ff_60, gd_36, gd_37, gd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_4 * fd_36[k]
                  + pa_x[k] * ff_60[k];

        t_61[k] = f_4 * fd_18[k]
                  + pb_y[k] * gd_36[k];

        t_62[k] = pb_z[k] * gd_36[k];

        t_63[k] = f_3 * fd_39[k]
                  + pb_x[k] * gd_39[k];

        t_64[k] = pb_z[k] * gd_37[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pa_x, pb_x, pb_z, fd_41, ff_66, ff_68, \
                         ff_69, gd_39, gd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_3 * fd_41[k]
                  + pb_x[k] * gd_41[k];

        t_66[k] = pa_x[k] * ff_66[k];

        t_67[k] = pb_z[k] * gd_39[k];

        t_68[k] = pa_x[k] * ff_68[k];

        t_69[k] = pa_x[k] * ff_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pa_z, pb_x, pb_z, fd_18, fd_46, ff_30, \
                         ff_31, ff_33, gd_42, gd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * ff_30[k];

        t_71[k] = pa_z[k] * ff_31[k];

        t_72[k] = f_3 * fd_18[k]
                  + pb_z[k] * gd_42[k];

        t_73[k] = pa_z[k] * ff_33[k];

        t_74[k] = f_3 * fd_46[k]
                  + pb_x[k] * gd_46[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, t_80, pa_x, pa_y, pb_x, fd_47, ff_50, \
                         ff_76, ff_77, ff_78, ff_79, gd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_3 * fd_47[k]
                  + pb_x[k] * gd_47[k];

        t_76[k] = pa_x[k] * ff_76[k];

        t_77[k] = pa_x[k] * ff_77[k];

        t_78[k] = pa_x[k] * ff_78[k];

        t_79[k] = pa_x[k] * ff_79[k];

        t_80[k] = pa_y[k] * ff_50[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, pa_y, pb_x, pb_y, fd_30, fd_51, fd_52, \
                         ff_52, ff_55, gd_48, gd_51, gd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_3 * fd_30[k]
                  + pb_y[k] * gd_48[k];

        t_82[k] = pa_y[k] * ff_52[k];

        t_83[k] = f_3 * fd_51[k]
                  + pb_x[k] * gd_51[k];

        t_84[k] = f_3 * fd_52[k]
                  + pb_x[k] * gd_52[k];

        t_85[k] = pa_y[k] * ff_55[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, t_91, pa_x, pb_y, fd_54, ff_86, ff_87, \
                         ff_88, ff_89, ff_90, gd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pa_x[k] * ff_86[k];

        t_87[k] = pa_x[k] * ff_87[k];

        t_88[k] = pa_x[k] * ff_88[k];

        t_89[k] = pa_x[k] * ff_89[k];

        t_90[k] = f_4 * fd_54[k]
                  + pa_x[k] * ff_90[k];

        t_91[k] = pb_y[k] * gd_54[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, fd_30, fd_57, fd_59, gd_54, \
                         gd_56, gd_57, gd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_4 * fd_30[k]
                  + pb_z[k] * gd_54[k];

        t_93[k] = f_3 * fd_57[k]
                  + pb_x[k] * gd_57[k];

        t_94[k] = pb_y[k] * gd_56[k];

        t_95[k] = f_3 * fd_59[k]
                  + pb_x[k] * gd_59[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pa_x, pb_x, pb_y, ff_96, ff_97, ff_99, \
                         gp0_30, gp1_30, gd_59, gd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_x[k] * ff_96[k];

        t_97[k] = pa_x[k] * ff_97[k];

        t_98[k] = pb_y[k] * gd_59[k];

        t_99[k] = pa_x[k] * ff_99[k];

        t_100[k] = f_1 * gp0_30[k]
                   - f_2 * gp1_30[k]
                   + pb_x[k] * gd_60[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, pb_x, pb_y, pb_z, fd_36, gd_60, \
                         gd_63, gd_64, gd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * fd_36[k]
                   + pb_y[k] * gd_60[k];

        t_102[k] = pb_z[k] * gd_60[k];

        t_103[k] = pb_x[k] * gd_63[k];

        t_104[k] = pb_x[k] * gd_64[k];

        t_105[k] = pb_x[k] * gd_65[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pb_y, pb_z, fd_39, fd_41, gp0_31, gp0_32, \
                         gp1_31, gp1_32, gd_63, gd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * fd_39[k]
                   + f_1 * gp0_31[k]
                   - f_2 * gp1_31[k]
                   + pb_y[k] * gd_63[k];

        t_107[k] = pb_z[k] * gd_63[k];

        t_108[k] = f_0 * fd_41[k]
                   + pb_y[k] * gd_65[k];

        t_109[k] = f_1 * gp0_32[k]
                   - f_2 * gp1_32[k]
                   + pb_z[k] * gd_65[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, t_115, pa_z, pb_x, pb_z, fd_36, \
                         ff_60, ff_61, gd_66, gd_69, gd_70, gd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * ff_60[k];

        t_111[k] = pa_z[k] * ff_61[k];

        t_112[k] = f_3 * fd_36[k]
                   + pb_z[k] * gd_66[k];

        t_113[k] = pb_x[k] * gd_69[k];

        t_114[k] = pb_x[k] * gd_70[k];

        t_115[k] = pb_x[k] * gd_71[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_z, pb_y, pb_z, fd_39, fd_41, fd_47, \
                         ff_66, ff_69, gd_69, gd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = pa_z[k] * ff_66[k];

        t_117[k] = f_3 * fd_39[k]
                   + pb_z[k] * gd_69[k];

        t_118[k] = f_4 * fd_47[k]
                   + pb_y[k] * gd_71[k];

        t_119[k] = f_4 * fd_41[k]
                   + pa_z[k] * ff_69[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pb_x, pb_y, pb_z, fd_42, fd_48, \
                         gp0_36, gp1_36, gd_72, gd_75, gd_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_1 * gp0_36[k]
                   - f_2 * gp1_36[k]
                   + pb_x[k] * gd_72[k];

        t_121[k] = f_7 * fd_48[k]
                   + pb_y[k] * gd_72[k];

        t_122[k] = f_7 * fd_42[k]
                   + pb_z[k] * gd_72[k];

        t_123[k] = pb_x[k] * gd_75[k];

        t_124[k] = pb_x[k] * gd_76[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_z, pb_x, pb_y, pb_z, df0_36, df1_36, \
                         fd_45, fd_53, ff_76, gd_75, gd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pb_x[k] * gd_77[k];

        t_126[k] = f_5 * df0_36[k]
                   - f_6 * df1_36[k]
                   + pa_z[k] * ff_76[k];

        t_127[k] = f_7 * fd_45[k]
                   + pb_z[k] * gd_75[k];

        t_128[k] = f_7 * fd_53[k]
                   + pb_y[k] * gd_77[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pa_y, pb_x, pb_y, df0_59, df1_59, \
                         fd_54, ff_89, ff_90, ff_92, gd_78, gd_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_5 * df0_59[k]
                   - f_6 * df1_59[k]
                   + pa_y[k] * ff_89[k];

        t_130[k] = pa_y[k] * ff_90[k];

        t_131[k] = f_3 * fd_54[k]
                   + pb_y[k] * gd_78[k];

        t_132[k] = pa_y[k] * ff_92[k];

        t_133[k] = pb_x[k] * gd_81[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, pa_y, pb_x, pb_y, pb_z, fd_51, \
                         fd_57, fd_59, ff_96, gd_81, gd_82, gd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = pb_x[k] * gd_82[k];

        t_135[k] = pb_x[k] * gd_83[k];

        t_136[k] = f_4 * fd_57[k]
                   + pa_y[k] * ff_96[k];

        t_137[k] = f_4 * fd_51[k]
                   + pb_z[k] * gd_81[k];

        t_138[k] = f_3 * fd_59[k]
                   + pb_y[k] * gd_83[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pa_y, pb_x, pb_y, pb_z, fd_54, \
                         ff_99, gp0_42, gp1_42, gd_84, gd_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_y[k] * ff_99[k];

        t_140[k] = f_1 * gp0_42[k]
                   - f_2 * gp1_42[k]
                   + pb_x[k] * gd_84[k];

        t_141[k] = pb_y[k] * gd_84[k];

        t_142[k] = f_0 * fd_54[k]
                   + pb_z[k] * gd_84[k];

        t_143[k] = pb_x[k] * gd_87[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, pb_x, pb_y, pb_z, fd_57, gp0_43, \
                         gp1_43, gd_87, gd_88, gd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = pb_x[k] * gd_88[k];

        t_145[k] = pb_x[k] * gd_89[k];

        t_146[k] = f_1 * gp0_43[k]
                   - f_2 * gp1_43[k]
                   + pb_y[k] * gd_87[k];

        t_147[k] = f_0 * fd_57[k]
                   + pb_z[k] * gd_87[k];

        t_148[k] = pb_y[k] * gd_89[k];
    }

#pragma omp simd aligned(t_149, pb_z, fd_59, gp0_44, gp1_44, gd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_0 * fd_59[k]
                   + f_1 * gp0_44[k]
                   - f_2 * gp1_44[k]
                   + pb_z[k] * gd_89[k];
    }
}

}  // namespace simdt2ceri
