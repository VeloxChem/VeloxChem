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


#include "SimdOverlapVrrRecHD.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_hd_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fd, const size_t gp, const size_t gd,
                          const size_t hs, const size_t hp, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_35 = buffer.data(fd + 35);
    const auto *fd_39 = buffer.data(fd + 39);
    const auto *fd_45 = buffer.data(fd + 45);
    const auto *fd_47 = buffer.data(fd + 47);
    const auto *fd_51 = buffer.data(fd + 51);
    const auto *fd_53 = buffer.data(fd + 53);
    const auto *fd_59 = buffer.data(fd + 59);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_23 = buffer.data(gp + 23);
    const auto *gp_25 = buffer.data(gp + 25);
    const auto *gp_29 = buffer.data(gp + 29);
    const auto *gp_30 = buffer.data(gp + 30);
    const auto *gp_31 = buffer.data(gp + 31);
    const auto *gp_35 = buffer.data(gp + 35);
    const auto *gp_36 = buffer.data(gp + 36);
    const auto *gp_37 = buffer.data(gp + 37);
    const auto *gp_38 = buffer.data(gp + 38);
    const auto *gp_40 = buffer.data(gp + 40);
    const auto *gp_41 = buffer.data(gp + 41);
    const auto *gp_42 = buffer.data(gp + 42);
    const auto *gp_43 = buffer.data(gp + 43);
    const auto *gp_44 = buffer.data(gp + 44);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_47 = buffer.data(gd + 47);
    const auto *gd_51 = buffer.data(gd + 51);
    const auto *gd_54 = buffer.data(gd + 54);
    const auto *gd_56 = buffer.data(gd + 56);
    const auto *gd_59 = buffer.data(gd + 59);
    const auto *gd_60 = buffer.data(gd + 60);
    const auto *gd_63 = buffer.data(gd + 63);
    const auto *gd_65 = buffer.data(gd + 65);
    const auto *gd_69 = buffer.data(gd + 69);
    const auto *gd_70 = buffer.data(gd + 70);
    const auto *gd_71 = buffer.data(gd + 71);
    const auto *gd_72 = buffer.data(gd + 72);
    const auto *gd_75 = buffer.data(gd + 75);
    const auto *gd_76 = buffer.data(gd + 76);
    const auto *gd_77 = buffer.data(gd + 77);
    const auto *gd_81 = buffer.data(gd + 81);
    const auto *gd_82 = buffer.data(gd + 82);
    const auto *gd_83 = buffer.data(gd + 83);
    const auto *gd_84 = buffer.data(gd + 84);
    const auto *gd_87 = buffer.data(gd + 87);
    const auto *gd_89 = buffer.data(gd + 89);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_15 = buffer.data(hs + 15);
    const auto *hs_17 = buffer.data(hs + 17);
    const auto *hs_18 = buffer.data(hs + 18);
    const auto *hs_20 = buffer.data(hs + 20);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);
    const auto *hp_30 = buffer.data(hp + 30);
    const auto *hp_31 = buffer.data(hp + 31);
    const auto *hp_35 = buffer.data(hp + 35);
    const auto *hp_37 = buffer.data(hp + 37);
    const auto *hp_38 = buffer.data(hp + 38);
    const auto *hp_40 = buffer.data(hp + 40);
    const auto *hp_42 = buffer.data(hp + 42);
    const auto *hp_44 = buffer.data(hp + 44);
    const auto *hp_45 = buffer.data(hp + 45);
    const auto *hp_46 = buffer.data(hp + 46);
    const auto *hp_47 = buffer.data(hp + 47);
    const auto *hp_49 = buffer.data(hp + 49);
    const auto *hp_50 = buffer.data(hp + 50);
    const auto *hp_51 = buffer.data(hp + 51);
    const auto *hp_52 = buffer.data(hp + 52);
    const auto *hp_53 = buffer.data(hp + 53);
    const auto *hp_54 = buffer.data(hp + 54);
    const auto *hp_55 = buffer.data(hp + 55);
    const auto *hp_56 = buffer.data(hp + 56);
    const auto *hp_58 = buffer.data(hp + 58);
    const auto *hp_59 = buffer.data(hp + 59);
    const auto *hp_60 = buffer.data(hp + 60);
    const auto *hp_61 = buffer.data(hp + 61);
    const auto *hp_62 = buffer.data(hp + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gp_0, hs_0, hp_0, \
                         hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = pb_y[k] * hp_0[k];

        t_2[k] = pb_z[k] * hp_0[k];

        t_3[k] = f_1 * hs_0[k]
                 + pb_y[k] * hp_1[k];

        t_4[k] = pb_y[k] * hp_2[k];

        t_5[k] = f_1 * hs_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pb_x, pb_z, fd_9, gp_4, gd_0, \
                         gd_9, hp_3, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * gd_0[k];

        t_7[k] = f_2 * gp_4[k]
                 + pb_x[k] * hp_4[k];

        t_8[k] = pb_z[k] * hp_3[k];

        t_9[k] = f_3 * fd_9[k]
                 + pa_x[k] * gd_9[k];

        t_10[k] = pb_z[k] * hp_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pa_y, pa_z, pb_x, pb_y, gp_8, \
                         gd_0, gd_3, gd_5, hp_6, hp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_y[k] * gd_5[k];

        t_12[k] = pa_z[k] * gd_0[k];

        t_13[k] = pb_y[k] * hp_6[k];

        t_14[k] = f_2 * gp_8[k]
                  + pb_x[k] * hp_8[k];

        t_15[k] = pa_z[k] * gd_3[k];

        t_16[k] = pb_y[k] * hp_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_x, pa_y, pb_x, pb_z, fd_0, fd_17, gp_10, \
                         gd_6, gd_17, hp_9, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * fd_17[k]
                  + pa_x[k] * gd_17[k];

        t_18[k] = f_1 * fd_0[k]
                  + pa_y[k] * gd_6[k];

        t_19[k] = f_3 * gp_10[k]
                  + pb_x[k] * hp_10[k];

        t_20[k] = pb_z[k] * hp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_x, pa_y, pa_z, pb_z, fd_21, gd_7, \
                         gd_12, gd_21, hs_3, hp_10, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_4 * fd_21[k]
                  + pa_x[k] * gd_21[k];

        t_22[k] = pb_z[k] * hp_10[k];

        t_23[k] = f_1 * hs_3[k]
                  + pb_z[k] * hp_11[k];

        t_24[k] = pa_y[k] * gd_12[k];

        t_25[k] = pa_z[k] * gd_7[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, fd_0, gp_8, gd_9, \
                         gd_12, gd_14, gd_17, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_y[k] * gd_14[k];

        t_27[k] = pa_z[k] * gd_9[k];

        t_28[k] = f_1 * gp_8[k]
                  + pb_y[k] * hp_14[k];

        t_29[k] = pa_y[k] * gd_17[k];

        t_30[k] = f_1 * fd_0[k]
                  + pa_z[k] * gd_12[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_x, pb_x, pb_y, fd_35, gp_17, gd_35, \
                         hs_5, hp_15, hp_16, hp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * hp_15[k];

        t_32[k] = f_3 * gp_17[k]
                  + pb_x[k] * hp_17[k];

        t_33[k] = f_1 * hs_5[k]
                  + pb_y[k] * hp_16[k];

        t_34[k] = pb_y[k] * hp_17[k];

        t_35[k] = f_4 * fd_35[k]
                  + pa_x[k] * gd_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pa_x, pa_y, pb_x, pb_z, fd_6, fd_39, \
                         gp_19, gd_18, gd_39, hp_18, hp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_4 * fd_6[k]
                  + pa_y[k] * gd_18[k];

        t_37[k] = f_4 * gp_19[k]
                  + pb_x[k] * hp_19[k];

        t_38[k] = pb_z[k] * hp_18[k];

        t_39[k] = f_1 * fd_39[k]
                  + pa_x[k] * gd_39[k];

        t_40[k] = pb_z[k] * hp_19[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_z, pb_x, pb_z, gp_23, gd_18, gd_19, \
                         gd_21, hs_6, hp_20, hp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * hs_6[k]
                  + pb_z[k] * hp_20[k];

        t_42[k] = pa_z[k] * gd_18[k];

        t_43[k] = pa_z[k] * gd_19[k];

        t_44[k] = f_4 * gp_23[k]
                  + pb_x[k] * hp_23[k];

        t_45[k] = pa_z[k] * gd_21[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pa_y, pb_x, pb_y, fd_47, gp_14, gp_25, \
                         gd_30, gd_47, hp_23, hp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_4 * gp_14[k]
                  + pb_y[k] * hp_23[k];

        t_47[k] = f_1 * fd_47[k]
                  + pa_x[k] * gd_47[k];

        t_48[k] = pa_y[k] * gd_30[k];

        t_49[k] = f_4 * gp_25[k]
                  + pb_x[k] * hp_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_x, pa_y, pb_y, fd_51, gp_17, gd_32, gd_35, \
                         gd_51, hp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * gd_32[k];

        t_51[k] = f_1 * fd_51[k]
                  + pa_x[k] * gd_51[k];

        t_52[k] = f_1 * gp_17[k]
                  + pb_y[k] * hp_26[k];

        t_53[k] = pa_y[k] * gd_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pa_z, pb_x, pb_y, fd_12, gp_29, gd_30, \
                         hs_9, hp_27, hp_28, hp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_4 * fd_12[k]
                  + pa_z[k] * gd_30[k];

        t_55[k] = pb_y[k] * hp_27[k];

        t_56[k] = f_4 * gp_29[k]
                  + pb_x[k] * hp_29[k];

        t_57[k] = f_1 * hs_9[k]
                  + pb_y[k] * hp_28[k];

        t_58[k] = pb_y[k] * hp_29[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, t_63, pa_x, pb_x, pb_z, fd_59, gp_30, gp_31, \
                         gd_59, gd_60, gd_63, hp_30, hp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * fd_59[k]
                  + pa_x[k] * gd_59[k];

        t_60[k] = f_4 * gp_30[k]
                  + pa_x[k] * gd_60[k];

        t_61[k] = f_1 * gp_31[k]
                  + pb_x[k] * hp_31[k];

        t_62[k] = pb_z[k] * hp_30[k];

        t_63[k] = pa_x[k] * gd_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, pa_x, pa_z, pb_x, pb_z, gp_35, gd_36, \
                         gd_37, gd_65, hp_31, hp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_z[k] * hp_31[k];

        t_65[k] = pa_x[k] * gd_65[k];

        t_66[k] = pa_z[k] * gd_36[k];

        t_67[k] = pa_z[k] * gd_37[k];

        t_68[k] = f_1 * gp_35[k]
                  + pb_x[k] * hp_35[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pa_x, pb_x, gp_36, gp_37, gd_69, gd_70, \
                         gd_71, gd_72, hp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pa_x[k] * gd_69[k];

        t_70[k] = pa_x[k] * gd_70[k];

        t_71[k] = pa_x[k] * gd_71[k];

        t_72[k] = f_4 * gp_36[k]
                  + pa_x[k] * gd_72[k];

        t_73[k] = f_1 * gp_37[k]
                  + pb_x[k] * hp_37[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, pa_x, pa_y, pb_x, gp_38, gd_54, gd_75, \
                         gd_76, gd_77, hp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_1 * gp_38[k]
                  + pb_x[k] * hp_38[k];

        t_75[k] = pa_x[k] * gd_75[k];

        t_76[k] = pa_x[k] * gd_76[k];

        t_77[k] = pa_x[k] * gd_77[k];

        t_78[k] = pa_y[k] * gd_54[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pa_x, pa_y, pb_x, gp_40, gd_56, gd_81, \
                         gd_82, gd_83, hp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_1 * gp_40[k]
                  + pb_x[k] * hp_40[k];

        t_80[k] = pa_y[k] * gd_56[k];

        t_81[k] = pa_x[k] * gd_81[k];

        t_82[k] = pa_x[k] * gd_82[k];

        t_83[k] = pa_x[k] * gd_83[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, pa_x, pb_x, pb_y, gp_42, gp_44, \
                         gd_84, gd_87, gd_89, hp_42, hp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * gp_42[k]
                  + pa_x[k] * gd_84[k];

        t_85[k] = pb_y[k] * hp_42[k];

        t_86[k] = f_1 * gp_44[k]
                  + pb_x[k] * hp_44[k];

        t_87[k] = pa_x[k] * gd_87[k];

        t_88[k] = pb_y[k] * hp_44[k];

        t_89[k] = pa_x[k] * gd_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, gp_31, hs_15, \
                         hp_45, hp_46, hp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_1 * hs_15[k]
                  + pb_x[k] * hp_45[k];

        t_91[k] = pb_x[k] * hp_46[k];

        t_92[k] = pb_x[k] * hp_47[k];

        t_93[k] = f_0 * gp_31[k]
                  + f_1 * hs_15[k]
                  + pb_y[k] * hp_46[k];

        t_94[k] = pb_z[k] * hp_46[k];

        t_95[k] = f_1 * hs_15[k]
                  + pb_z[k] * hp_47[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pa_z, pb_x, pb_y, gp_35, gd_60, gd_63, \
                         hp_49, hp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * gd_60[k];

        t_97[k] = pb_x[k] * hp_49[k];

        t_98[k] = pb_x[k] * hp_50[k];

        t_99[k] = pa_z[k] * gd_63[k];

        t_100[k] = f_2 * gp_35[k]
                   + pb_y[k] * hp_50[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, pa_y, pa_z, pb_x, fd_39, fd_47, \
                         gd_69, gd_71, hs_17, hp_51, hp_52, hp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_3 * fd_47[k]
                   + pa_y[k] * gd_71[k];

        t_102[k] = f_1 * hs_17[k]
                   + pb_x[k] * hp_51[k];

        t_103[k] = pb_x[k] * hp_52[k];

        t_104[k] = pb_x[k] * hp_53[k];

        t_105[k] = f_1 * fd_39[k]
                   + pa_z[k] * gd_69[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, pa_y, pb_x, pb_y, fd_53, gp_38, \
                         gd_77, hs_18, hp_53, hp_54, hp_55, hp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_3 * gp_38[k]
                   + pb_y[k] * hp_53[k];

        t_107[k] = f_4 * fd_53[k]
                   + pa_y[k] * gd_77[k];

        t_108[k] = f_1 * hs_18[k]
                   + pb_x[k] * hp_54[k];

        t_109[k] = pb_x[k] * hp_55[k];

        t_110[k] = pb_x[k] * hp_56[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_y, pa_z, pb_y, fd_45, fd_59, gp_41, \
                         gd_75, gd_83, gd_84, hp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_4 * fd_45[k]
                   + pa_z[k] * gd_75[k];

        t_112[k] = f_4 * gp_41[k]
                   + pb_y[k] * hp_56[k];

        t_113[k] = f_1 * fd_59[k]
                   + pa_y[k] * gd_83[k];

        t_114[k] = pa_y[k] * gd_84[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pa_y, pb_x, pb_y, gp_43, gp_44, \
                         gd_87, gd_89, hp_58, hp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pb_x[k] * hp_58[k];

        t_116[k] = pb_x[k] * hp_59[k];

        t_117[k] = f_4 * gp_43[k]
                   + pa_y[k] * gd_87[k];

        t_118[k] = f_1 * gp_44[k]
                   + pb_y[k] * hp_59[k];

        t_119[k] = pa_y[k] * gd_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, t_125, pb_x, pb_y, pb_z, gp_44, \
                         hs_20, hp_60, hp_61, hp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_1 * hs_20[k]
                   + pb_x[k] * hp_60[k];

        t_121[k] = pb_x[k] * hp_61[k];

        t_122[k] = pb_x[k] * hp_62[k];

        t_123[k] = f_1 * hs_20[k]
                   + pb_y[k] * hp_61[k];

        t_124[k] = pb_y[k] * hp_62[k];

        t_125[k] = f_0 * gp_44[k]
                   + f_1 * hs_20[k]
                   + pb_z[k] * hp_62[k];
    }
}

}  // namespace simdovl
