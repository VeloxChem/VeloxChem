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


#include "SimdElectronRepulsionVrrRecHD.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_hd_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fd0, const size_t fd1,
                                     const size_t gp, const size_t gd, const size_t hs0,
                                     const size_t hs1, const size_t hp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 2.0 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.5 / p;
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 0.5 / p;

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
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);
    const auto *gp_22 = buffer.data(gp + 22);
    const auto *gp_23 = buffer.data(gp + 23);
    const auto *gp_24 = buffer.data(gp + 24);
    const auto *gp_25 = buffer.data(gp + 25);

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

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);
    const auto *hp_30 = buffer.data(hp + 30);
    const auto *hp_31 = buffer.data(hp + 31);
    const auto *hp_32 = buffer.data(hp + 32);
    const auto *hp_33 = buffer.data(hp + 33);
    const auto *hp_34 = buffer.data(hp + 34);
    const auto *hp_35 = buffer.data(hp + 35);
    const auto *hp_36 = buffer.data(hp + 36);
    const auto *hp_37 = buffer.data(hp + 37);
    const auto *hp_38 = buffer.data(hp + 38);
    const auto *hp_39 = buffer.data(hp + 39);
    const auto *hp_40 = buffer.data(hp + 40);
    const auto *hp_41 = buffer.data(hp + 41);
    const auto *hp_42 = buffer.data(hp + 42);
    const auto *hp_43 = buffer.data(hp + 43);
    const auto *hp_44 = buffer.data(hp + 44);
    const auto *hp_45 = buffer.data(hp + 45);
    const auto *hp_46 = buffer.data(hp + 46);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gp_0, hs0_0, hs1_0, \
                         hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = pb_y[k] * hp_0[k];

        t_2[k] = pb_z[k] * hp_0[k];

        t_3[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_4[k] = pb_y[k] * hp_2[k];

        t_5[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pa_y, pb_x, pb_z, gp_1, gp_3, gd_0, \
                         gd_1, gd_2, hp_3, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * gd_0[k];

        t_7[k] = f_3 * gp_3[k]
                 + pb_x[k] * hp_4[k];

        t_8[k] = pb_z[k] * hp_3[k];

        t_9[k] = f_4 * gp_1[k]
                 + pa_y[k] * gd_1[k];

        t_10[k] = pb_z[k] * hp_4[k];

        t_11[k] = pa_y[k] * gd_2[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, gp_2, gp_4, \
                         gd_0, gd_1, gd_2, hp_5, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_z[k] * gd_0[k];

        t_13[k] = pb_y[k] * hp_5[k];

        t_14[k] = f_3 * gp_4[k]
                  + pb_x[k] * hp_6[k];

        t_15[k] = pa_z[k] * gd_1[k];

        t_16[k] = pb_y[k] * hp_6[k];

        t_17[k] = f_4 * gp_2[k]
                  + pa_z[k] * gd_2[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pb_x, pb_z, fd0_0, fd1_0, gp_5, gd_3, hp_7, \
                         hp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * fd0_0[k]
                  - f_6 * fd1_0[k]
                  + pa_y[k] * gd_3[k];

        t_19[k] = f_7 * gp_5[k]
                  + pb_x[k] * hp_8[k];

        t_20[k] = pb_z[k] * hp_7[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_x, pa_y, pb_z, fd0_3, fd1_3, gd_6, gd_11, \
                         hs0_1, hs1_1, hp_8, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * fd0_3[k]
                  - f_9 * fd1_3[k]
                  + pa_x[k] * gd_11[k];

        t_22[k] = pb_z[k] * hp_8[k];

        t_23[k] = f_1 * hs0_1[k]
                  - f_2 * hs1_1[k]
                  + pb_z[k] * hp_9[k];

        t_24[k] = pa_y[k] * gd_6[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pa_z, pb_y, gp_4, gd_4, gd_5, \
                         gd_7, gd_8, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_z[k] * gd_4[k];

        t_26[k] = pa_y[k] * gd_7[k];

        t_27[k] = pa_z[k] * gd_5[k];

        t_28[k] = f_10 * gp_4[k]
                  + pb_y[k] * hp_10[k];

        t_29[k] = pa_y[k] * gd_8[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_z, pb_x, pb_y, fd0_0, fd1_0, gp_9, gd_6, \
                         hs0_2, hs1_2, hp_11, hp_12, hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * fd0_0[k]
                  - f_6 * fd1_0[k]
                  + pa_z[k] * gd_6[k];

        t_31[k] = pb_y[k] * hp_11[k];

        t_32[k] = f_7 * gp_9[k]
                  + pb_x[k] * hp_13[k];

        t_33[k] = f_1 * hs0_2[k]
                  - f_2 * hs1_2[k]
                  + pb_y[k] * hp_12[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_x, pa_y, pb_y, fd0_1, fd0_4, fd1_1, fd1_4, gd_9, \
                         gd_16, hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * hp_13[k];

        t_35[k] = f_8 * fd0_4[k]
                  - f_9 * fd1_4[k]
                  + pa_x[k] * gd_16[k];

        t_36[k] = f_8 * fd0_1[k]
                  - f_9 * fd1_1[k]
                  + pa_y[k] * gd_9[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_x, pb_x, pb_z, fd0_5, fd1_5, gp_10, gd_19, \
                         hp_14, hp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_4 * gp_10[k]
                  + pb_x[k] * hp_15[k];

        t_38[k] = pb_z[k] * hp_14[k];

        t_39[k] = f_5 * fd0_5[k]
                  - f_6 * fd1_5[k]
                  + pa_x[k] * gd_19[k];

        t_40[k] = pb_z[k] * hp_15[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_z, pb_x, pb_z, gp_11, gd_9, gd_10, \
                         gd_11, hs0_3, hs1_3, hp_16, hp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_16[k];

        t_42[k] = pa_z[k] * gd_9[k];

        t_43[k] = pa_z[k] * gd_10[k];

        t_44[k] = f_4 * gp_11[k]
                  + pb_x[k] * hp_17[k];

        t_45[k] = pa_z[k] * gd_11[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pa_z, pb_x, pb_y, gp_6, gp_7, gp_12, \
                         gd_12, gd_13, hp_17, hp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_4 * gp_7[k]
                  + pb_y[k] * hp_17[k];

        t_47[k] = f_4 * gp_6[k]
                  + pa_z[k] * gd_12[k];

        t_48[k] = pa_y[k] * gd_13[k];

        t_49[k] = f_4 * gp_12[k]
                  + pb_x[k] * hp_18[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pb_y, gp_8, gp_9, gd_14, gd_15, gd_16, \
                         hp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * gd_14[k];

        t_51[k] = f_4 * gp_8[k]
                  + pa_y[k] * gd_15[k];

        t_52[k] = f_10 * gp_9[k]
                  + pb_y[k] * hp_19[k];

        t_53[k] = pa_y[k] * gd_16[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_z, pb_x, pb_y, fd0_2, fd1_2, gp_13, gd_13, \
                         hs0_4, hs1_4, hp_20, hp_21, hp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_8 * fd0_2[k]
                  - f_9 * fd1_2[k]
                  + pa_z[k] * gd_13[k];

        t_55[k] = pb_y[k] * hp_20[k];

        t_56[k] = f_4 * gp_13[k]
                  + pb_x[k] * hp_22[k];

        t_57[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_21[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_x, pb_x, pb_y, fd0_8, fd1_8, gp_14, gp_15, \
                         gd_22, gd_23, hp_22, hp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_y[k] * hp_22[k];

        t_59[k] = f_5 * fd0_8[k]
                  - f_6 * fd1_8[k]
                  + pa_x[k] * gd_22[k];

        t_60[k] = f_4 * gp_14[k]
                  + pa_x[k] * gd_23[k];

        t_61[k] = f_10 * gp_15[k]
                  + pb_x[k] * hp_24[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, t_67, pa_x, pa_z, pb_z, gd_17, gd_18, \
                         gd_24, gd_25, hp_23, hp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_z[k] * hp_23[k];

        t_63[k] = pa_x[k] * gd_24[k];

        t_64[k] = pb_z[k] * hp_24[k];

        t_65[k] = pa_x[k] * gd_25[k];

        t_66[k] = pa_z[k] * gd_17[k];

        t_67[k] = pa_z[k] * gd_18[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_x, pb_x, gp_17, gp_18, gd_26, gd_27, \
                         gd_28, gd_29, hp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_10 * gp_17[k]
                  + pb_x[k] * hp_25[k];

        t_69[k] = pa_x[k] * gd_26[k];

        t_70[k] = pa_x[k] * gd_27[k];

        t_71[k] = pa_x[k] * gd_28[k];

        t_72[k] = f_4 * gp_18[k]
                  + pa_x[k] * gd_29[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_x, pb_x, gp_19, gp_20, gd_30, gd_31, \
                         gd_32, hp_26, hp_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_10 * gp_19[k]
                  + pb_x[k] * hp_26[k];

        t_74[k] = f_10 * gp_20[k]
                  + pb_x[k] * hp_27[k];

        t_75[k] = pa_x[k] * gd_30[k];

        t_76[k] = pa_x[k] * gd_31[k];

        t_77[k] = pa_x[k] * gd_32[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, pa_x, pa_y, pb_x, gp_21, gd_20, \
                         gd_21, gd_33, gd_34, gd_35, hp_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_y[k] * gd_20[k];

        t_79[k] = f_10 * gp_21[k]
                  + pb_x[k] * hp_28[k];

        t_80[k] = pa_y[k] * gd_21[k];

        t_81[k] = pa_x[k] * gd_33[k];

        t_82[k] = pa_x[k] * gd_34[k];

        t_83[k] = pa_x[k] * gd_35[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, pa_x, pb_x, pb_y, gp_23, gp_25, \
                         gd_36, gd_37, gd_38, hp_29, hp_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * gp_23[k]
                  + pa_x[k] * gd_36[k];

        t_85[k] = pb_y[k] * hp_29[k];

        t_86[k] = f_10 * gp_25[k]
                  + pb_x[k] * hp_30[k];

        t_87[k] = pa_x[k] * gd_37[k];

        t_88[k] = pb_y[k] * hp_30[k];

        t_89[k] = pa_x[k] * gd_38[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, gp_15, hs0_5, \
                         hs1_5, hp_31, hp_32, hp_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_x[k] * hp_31[k];

        t_91[k] = pb_x[k] * hp_32[k];

        t_92[k] = pb_x[k] * hp_33[k];

        t_93[k] = f_0 * gp_15[k]
                  + f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_y[k] * hp_32[k];

        t_94[k] = pb_z[k] * hp_32[k];

        t_95[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_z[k] * hp_33[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, t_101, pa_z, pb_x, pb_y, gp_16, gp_17, \
                         gd_23, gd_24, gd_25, hp_34, hp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * gd_23[k];

        t_97[k] = pb_x[k] * hp_34[k];

        t_98[k] = pb_x[k] * hp_35[k];

        t_99[k] = pa_z[k] * gd_24[k];

        t_100[k] = f_3 * gp_17[k]
                   + pb_y[k] * hp_35[k];

        t_101[k] = f_4 * gp_16[k]
                   + pa_z[k] * gd_25[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pa_z, pb_x, fd0_5, fd1_5, gd_26, hs0_6, \
                         hs1_6, hp_36, hp_37, hp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_1 * hs0_6[k]
                   - f_2 * hs1_6[k]
                   + pb_x[k] * hp_36[k];

        t_103[k] = pb_x[k] * hp_37[k];

        t_104[k] = pb_x[k] * hp_38[k];

        t_105[k] = f_5 * fd0_5[k]
                   - f_6 * fd1_5[k]
                   + pa_z[k] * gd_26[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pa_y, pb_x, pb_y, fd0_7, fd1_7, gp_20, \
                         gd_32, hs0_7, hs1_7, hp_38, hp_39, hp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_7 * gp_20[k]
                   + pb_y[k] * hp_38[k];

        t_107[k] = f_8 * fd0_7[k]
                   - f_9 * fd1_7[k]
                   + pa_y[k] * gd_32[k];

        t_108[k] = f_1 * hs0_7[k]
                   - f_2 * hs1_7[k]
                   + pb_x[k] * hp_39[k];

        t_109[k] = pb_x[k] * hp_40[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_y, pa_z, pb_x, pb_y, fd0_6, fd0_8, \
                         fd1_6, fd1_8, gp_22, gd_30, gd_35, hp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pb_x[k] * hp_41[k];

        t_111[k] = f_8 * fd0_6[k]
                   - f_9 * fd1_6[k]
                   + pa_z[k] * gd_30[k];

        t_112[k] = f_4 * gp_22[k]
                   + pb_y[k] * hp_41[k];

        t_113[k] = f_5 * fd0_8[k]
                   - f_6 * fd1_8[k]
                   + pa_y[k] * gd_35[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, pa_y, pb_x, pb_y, gp_24, \
                         gp_25, gd_36, gd_37, gd_38, hp_42, hp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = pa_y[k] * gd_36[k];

        t_115[k] = pb_x[k] * hp_42[k];

        t_116[k] = pb_x[k] * hp_43[k];

        t_117[k] = f_4 * gp_24[k]
                   + pa_y[k] * gd_37[k];

        t_118[k] = f_10 * gp_25[k]
                   + pb_y[k] * hp_43[k];

        t_119[k] = pa_y[k] * gd_38[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, t_125, pb_x, pb_y, pb_z, gp_25, \
                         hs0_8, hs1_8, hp_44, hp_45, hp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_1 * hs0_8[k]
                   - f_2 * hs1_8[k]
                   + pb_x[k] * hp_44[k];

        t_121[k] = pb_x[k] * hp_45[k];

        t_122[k] = pb_x[k] * hp_46[k];

        t_123[k] = f_1 * hs0_8[k]
                   - f_2 * hs1_8[k]
                   + pb_y[k] * hp_45[k];

        t_124[k] = pb_y[k] * hp_46[k];

        t_125[k] = f_0 * gp_25[k]
                   + f_1 * hs0_8[k]
                   - f_2 * hs1_8[k]
                   + pb_z[k] * hp_46[k];
    }
}

auto
compute_prim_hd_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fd0, const size_t fd1,
                                     const size_t gp, const size_t gd, const size_t hs0,
                                     const size_t hs1, const size_t hp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 2.0 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.5 / p;
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_14 = buffer.data(fd1 + 14);
    const auto *fd1_19 = buffer.data(fd1 + 19);
    const auto *fd1_22 = buffer.data(fd1 + 22);

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
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);

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

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_y, pb_x, pb_y, pb_z, gp_0, gd_0, hs0_0, \
                         hs1_0, hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = pb_z[k] * hp_0[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_3[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];

        t_4[k] = pa_y[k] * gd_0[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_x, gp_1, gp_3, gp_4, gd_0, \
                         gd_1, gd_2, hp_3, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * gp_3[k]
                 + pb_x[k] * hp_3[k];

        t_6[k] = f_4 * gp_1[k]
                 + pa_y[k] * gd_1[k];

        t_7[k] = pa_y[k] * gd_2[k];

        t_8[k] = pa_z[k] * gd_0[k];

        t_9[k] = f_3 * gp_4[k]
                 + pb_x[k] * hp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_x, fd0_0, fd1_0, gp_2, gp_5, \
                         gd_1, gd_2, gd_3, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * gd_1[k];

        t_11[k] = f_4 * gp_2[k]
                  + pa_z[k] * gd_2[k];

        t_12[k] = f_5 * fd0_0[k]
                  - f_6 * fd1_0[k]
                  + pa_y[k] * gd_3[k];

        t_13[k] = f_7 * gp_5[k]
                  + pb_x[k] * hp_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_z, fd0_3, fd1_8, gd_4, gd_8, hs0_1, \
                         hs1_1, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_8 * fd0_3[k]
                  - f_9 * fd1_8[k]
                  + pa_x[k] * gd_8[k];

        t_15[k] = f_1 * hs0_1[k]
                  - f_2 * hs1_1[k]
                  + pb_z[k] * hp_6[k];

        t_16[k] = pa_z[k] * gd_4[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pa_z, pb_x, pb_y, fd0_0, fd1_0, gp_4, \
                         gp_9, gd_5, gd_6, hp_7, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_10 * gp_4[k]
                  + pb_y[k] * hp_7[k];

        t_18[k] = pa_y[k] * gd_6[k];

        t_19[k] = f_5 * fd0_0[k]
                  - f_6 * fd1_0[k]
                  + pa_z[k] * gd_5[k];

        t_20[k] = f_7 * gp_9[k]
                  + pb_x[k] * hp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_y, fd0_1, fd0_4, fd1_3, fd1_10, \
                         gd_7, gd_12, hs0_2, hs1_2, hp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * hs0_2[k]
                  - f_2 * hs1_2[k]
                  + pb_y[k] * hp_8[k];

        t_22[k] = f_8 * fd0_4[k]
                  - f_9 * fd1_10[k]
                  + pa_x[k] * gd_12[k];

        t_23[k] = f_8 * fd0_1[k]
                  - f_9 * fd1_3[k]
                  + pa_y[k] * gd_7[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pb_x, pb_z, fd0_5, fd1_12, gp_10, gd_14, \
                         hs0_3, hs1_3, hp_10, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_4 * gp_10[k]
                  + pb_x[k] * hp_10[k];

        t_25[k] = f_5 * fd0_5[k]
                  - f_6 * fd1_12[k]
                  + pa_x[k] * gd_14[k];

        t_26[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_11[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_y, pa_z, pb_y, gp_6, gp_7, gd_7, \
                         gd_8, gd_9, gd_10, hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_z[k] * gd_7[k];

        t_28[k] = pa_z[k] * gd_8[k];

        t_29[k] = f_4 * gp_7[k]
                  + pb_y[k] * hp_12[k];

        t_30[k] = f_4 * gp_6[k]
                  + pa_z[k] * gd_9[k];

        t_31[k] = pa_y[k] * gd_10[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_y, pa_z, pb_y, fd0_2, fd1_5, gp_8, gp_9, \
                         gd_10, gd_11, gd_12, hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_4 * gp_8[k]
                  + pa_y[k] * gd_11[k];

        t_33[k] = f_10 * gp_9[k]
                  + pb_y[k] * hp_13[k];

        t_34[k] = pa_y[k] * gd_12[k];

        t_35[k] = f_8 * fd0_2[k]
                  - f_9 * fd1_5[k]
                  + pa_z[k] * gd_10[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_x, pb_x, pb_y, fd0_8, fd1_22, gp_11, gd_16, \
                         hs0_4, hs1_4, hp_14, hp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_4 * gp_11[k]
                  + pb_x[k] * hp_15[k];

        t_37[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_14[k];

        t_38[k] = f_5 * fd0_8[k]
                  - f_6 * fd1_22[k]
                  + pa_x[k] * gd_16[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_x, pa_z, pb_x, gp_12, gp_13, gd_13, \
                         gd_17, gd_18, gd_19, hp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_4 * gp_12[k]
                  + pa_x[k] * gd_17[k];

        t_40[k] = f_10 * gp_13[k]
                  + pb_x[k] * hp_16[k];

        t_41[k] = pa_x[k] * gd_18[k];

        t_42[k] = pa_x[k] * gd_19[k];

        t_43[k] = pa_z[k] * gd_13[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, pa_x, gp_16, gd_21, gd_22, gd_23, \
                         gd_24, gd_25, gd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_x[k] * gd_21[k];

        t_45[k] = pa_x[k] * gd_22[k];

        t_46[k] = f_4 * gp_16[k]
                  + pa_x[k] * gd_23[k];

        t_47[k] = pa_x[k] * gd_24[k];

        t_48[k] = pa_x[k] * gd_25[k];

        t_49[k] = pa_x[k] * gd_26[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_x, pa_y, pb_x, gp_19, gp_21, gd_15, \
                         gd_27, gd_28, gd_30, hp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * gd_15[k];

        t_51[k] = pa_x[k] * gd_27[k];

        t_52[k] = pa_x[k] * gd_28[k];

        t_53[k] = f_4 * gp_19[k]
                  + pa_x[k] * gd_30[k];

        t_54[k] = f_10 * gp_21[k]
                  + pb_x[k] * hp_17[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_x, pb_x, pb_y, pb_z, gp_13, gd_31, \
                         gd_32, hs0_5, hs1_5, hp_18, hp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_x[k] * gd_31[k];

        t_56[k] = pa_x[k] * gd_32[k];

        t_57[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_x[k] * hp_18[k];

        t_58[k] = f_0 * gp_13[k]
                  + f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_y[k] * hp_19[k];

        t_59[k] = pb_z[k] * hp_19[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_y, pb_z, gp_15, gd_17, gd_18, hs0_5, \
                         hs1_5, hp_20, hp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_z[k] * hp_20[k];

        t_61[k] = pa_z[k] * gd_17[k];

        t_62[k] = pa_z[k] * gd_18[k];

        t_63[k] = f_3 * gp_15[k]
                  + pb_y[k] * hp_21[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_z, pb_x, fd0_5, fd1_12, gp_14, gd_19, gd_20, \
                         hs0_6, hs1_6, hp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_4 * gp_14[k]
                  + pa_z[k] * gd_19[k];

        t_65[k] = f_1 * hs0_6[k]
                  - f_2 * hs1_6[k]
                  + pb_x[k] * hp_22[k];

        t_66[k] = f_5 * fd0_5[k]
                  - f_6 * fd1_12[k]
                  + pa_z[k] * gd_20[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pa_y, pb_x, pb_y, fd0_7, fd1_19, gp_17, gd_26, \
                         hs0_7, hs1_7, hp_23, hp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_7 * gp_17[k]
                  + pb_y[k] * hp_23[k];

        t_68[k] = f_8 * fd0_7[k]
                  - f_9 * fd1_19[k]
                  + pa_y[k] * gd_26[k];

        t_69[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_24[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_y, pa_z, pb_y, fd0_6, fd0_8, fd1_14, \
                         fd1_22, gp_18, gd_24, gd_29, gd_30, hp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_8 * fd0_6[k]
                  - f_9 * fd1_14[k]
                  + pa_z[k] * gd_24[k];

        t_71[k] = f_4 * gp_18[k]
                  + pb_y[k] * hp_25[k];

        t_72[k] = f_5 * fd0_8[k]
                  - f_6 * fd1_22[k]
                  + pa_y[k] * gd_29[k];

        t_73[k] = pa_y[k] * gd_30[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_y, pb_x, pb_y, gp_20, gp_21, gd_31, gd_32, \
                         hs0_8, hs1_8, hp_26, hp_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_4 * gp_20[k]
                  + pa_y[k] * gd_31[k];

        t_75[k] = f_10 * gp_21[k]
                  + pb_y[k] * hp_26[k];

        t_76[k] = pa_y[k] * gd_32[k];

        t_77[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_27[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pb_y, pb_z, gp_21, hs0_8, hs1_8, hp_28, \
                         hp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_y[k] * hp_28[k];

        t_79[k] = pb_y[k] * hp_29[k];

        t_80[k] = f_0 * gp_21[k]
                  + f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_z[k] * hp_29[k];
    }
}

auto
compute_prim_hd_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fd0, const size_t fd1,
                                     const size_t gp, const size_t gd, const size_t hs0,
                                     const size_t hs1, const size_t hp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / alpha;
    const auto f_6 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_11 = buffer.data(gp + 11);

    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs0_0, hs1_0, hp_0, hp_1, \
                         hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, fd0_0, fd0_3, fd1_0, fd1_3, gd_3, \
                         gd_6, hs0_1, hs1_1, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_y[k] * gd_3[k];

        t_4[k] = f_5 * fd0_3[k]
                 - f_6 * fd1_3[k]
                 + pa_x[k] * gd_6[k];

        t_5[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_1[k]
                 + pb_z[k] * hp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, fd0_0, fd0_4, fd1_0, fd1_4, gd_4, \
                         gd_10, hs0_2, hs1_2, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_z[k] * gd_4[k];

        t_7[k] = f_1 * hs0_2[k]
                 - f_2 * hs1_2[k]
                 + pb_y[k] * hp_4[k];

        t_8[k] = f_5 * fd0_4[k]
                 - f_6 * fd1_4[k]
                 + pa_x[k] * gd_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, fd0_1, fd0_5, fd1_1, fd1_5, gd_5, \
                         gd_11, hs0_3, hs1_3, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * fd0_1[k]
                 - f_6 * fd1_1[k]
                 + pa_y[k] * gd_5[k];

        t_10[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_5[k]
                  + pa_x[k] * gd_11[k];

        t_11[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, fd0_2, fd0_8, fd1_2, fd1_8, gd_8, \
                         gd_12, hs0_4, hs1_4, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * fd0_2[k]
                  - f_6 * fd1_2[k]
                  + pa_z[k] * gd_8[k];

        t_13[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];

        t_14[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pa_x[k] * gd_12[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, gp_6, hs0_5, hs0_6, hs1_5, \
                         hs1_6, hp_7, hp_8, hp_9, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_x[k] * hp_7[k];

        t_16[k] = f_0 * gp_6[k]
                  + f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_y[k] * hp_8[k];

        t_17[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_z[k] * hp_9[k];

        t_18[k] = f_1 * hs0_6[k]
                  - f_2 * hs1_6[k]
                  + pb_x[k] * hp_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pa_z, pb_x, fd0_5, fd0_7, fd1_5, fd1_7, \
                         gd_16, gd_19, hs0_7, hs1_7, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_5[k]
                  + pa_z[k] * gd_16[k];

        t_20[k] = f_5 * fd0_7[k]
                  - f_6 * fd1_7[k]
                  + pa_y[k] * gd_19[k];

        t_21[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pa_z, pb_x, fd0_6, fd0_8, fd1_6, fd1_8, \
                         gd_18, gd_20, hs0_8, hs1_8, hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * fd0_6[k]
                  - f_6 * fd1_6[k]
                  + pa_z[k] * gd_18[k];

        t_23[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pa_y[k] * gd_20[k];

        t_24[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_12[k];
    }

#pragma omp simd aligned(t_25, t_26, pb_y, pb_z, gp_11, hs0_8, hs1_8, hp_13, \
                         hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_y[k] * hp_13[k];

        t_26[k] = f_0 * gp_11[k]
                  + f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fd0, const size_t fd1,
                                     const size_t gp, const size_t gd, const size_t hs0,
                                     const size_t hs1, const size_t hp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / alpha;
    const auto f_6 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_14 = buffer.data(fd1 + 14);
    const auto *fd1_17 = buffer.data(fd1 + 17);
    const auto *fd1_20 = buffer.data(fd1 + 20);
    const auto *fd1_23 = buffer.data(fd1 + 23);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_11 = buffer.data(gp + 11);

    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_32 = buffer.data(gd + 32);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs0_0, hs1_0, hp_0, hp_1, \
                         hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, fd0_0, fd0_3, fd1_0, fd1_10, gd_3, \
                         gd_10, hs0_1, hs1_1, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_y[k] * gd_3[k];

        t_4[k] = f_5 * fd0_3[k]
                 - f_6 * fd1_10[k]
                 + pa_x[k] * gd_10[k];

        t_5[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_1[k]
                 + pb_z[k] * hp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, fd0_0, fd0_4, fd1_0, fd1_12, gd_6, \
                         gd_16, hs0_2, hs1_2, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_z[k] * gd_6[k];

        t_7[k] = f_1 * hs0_2[k]
                 - f_2 * hs1_2[k]
                 + pb_y[k] * hp_4[k];

        t_8[k] = f_5 * fd0_4[k]
                 - f_6 * fd1_12[k]
                 + pa_x[k] * gd_16[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, fd0_1, fd0_5, fd1_3, fd1_14, gd_9, \
                         gd_18, hs0_3, hs1_3, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * fd0_1[k]
                 - f_6 * fd1_3[k]
                 + pa_y[k] * gd_9[k];

        t_10[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_14[k]
                  + pa_x[k] * gd_18[k];

        t_11[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, fd0_2, fd0_8, fd1_6, fd1_23, \
                         gd_14, gd_21, hs0_4, hs1_4, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * fd0_2[k]
                  - f_6 * fd1_6[k]
                  + pa_z[k] * gd_14[k];

        t_13[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];

        t_14[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_23[k]
                  + pa_x[k] * gd_21[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, gp_6, hs0_5, hs0_6, hs1_5, \
                         hs1_6, hp_7, hp_8, hp_9, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_x[k] * hp_7[k];

        t_16[k] = f_0 * gp_6[k]
                  + f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_y[k] * hp_8[k];

        t_17[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_z[k] * hp_9[k];

        t_18[k] = f_1 * hs0_6[k]
                  - f_2 * hs1_6[k]
                  + pb_x[k] * hp_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pa_z, pb_x, fd0_5, fd0_7, fd1_14, fd1_20, \
                         gd_26, gd_30, hs0_7, hs1_7, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_14[k]
                  + pa_z[k] * gd_26[k];

        t_20[k] = f_5 * fd0_7[k]
                  - f_6 * fd1_20[k]
                  + pa_y[k] * gd_30[k];

        t_21[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pa_z, pb_x, fd0_6, fd0_8, fd1_17, fd1_23, \
                         gd_29, gd_32, hs0_8, hs1_8, hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * fd0_6[k]
                  - f_6 * fd1_17[k]
                  + pa_z[k] * gd_29[k];

        t_23[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_23[k]
                  + pa_y[k] * gd_32[k];

        t_24[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_12[k];
    }

#pragma omp simd aligned(t_25, t_26, pb_y, pb_z, gp_11, hs0_8, hs1_8, hp_13, \
                         hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_y[k] * hp_13[k];

        t_26[k] = f_0 * gp_11[k]
                  + f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fd0, const size_t fd1,
                                     const size_t gp, const size_t gd, const size_t hs0,
                                     const size_t hs1, const size_t hp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_3 = buffer.data(fd0 + 3);
    const auto *fd0_6 = buffer.data(fd0 + 6);
    const auto *fd0_10 = buffer.data(fd0 + 10);
    const auto *fd0_12 = buffer.data(fd0 + 12);
    const auto *fd0_14 = buffer.data(fd0 + 14);
    const auto *fd0_17 = buffer.data(fd0 + 17);
    const auto *fd0_20 = buffer.data(fd0 + 20);
    const auto *fd0_23 = buffer.data(fd0 + 23);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_14 = buffer.data(fd1 + 14);
    const auto *fd1_17 = buffer.data(fd1 + 17);
    const auto *fd1_20 = buffer.data(fd1 + 20);

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
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, gp_0, gd_0, hs0_0, hs1_0, \
                         hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];

        t_3[k] = pa_y[k] * gd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, fd0_0, fd1_0, gp_1, gp_2, \
                         gd_0, gd_1, gd_2, gd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * gp_1[k]
                 + pa_y[k] * gd_1[k];

        t_5[k] = pa_y[k] * gd_2[k];

        t_6[k] = pa_z[k] * gd_0[k];

        t_7[k] = pa_z[k] * gd_1[k];

        t_8[k] = f_3 * gp_2[k]
                 + pa_z[k] * gd_2[k];

        t_9[k] = f_4 * fd0_0[k]
                 - f_5 * fd1_0[k]
                 + pa_y[k] * gd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, pa_z, pb_z, fd0_10, fd1_8, gd_4, \
                         gd_6, gd_8, hs0_1, hs1_1, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_6 * fd0_10[k]
                  - f_7 * fd1_8[k]
                  + pa_x[k] * gd_8[k];

        t_11[k] = f_1 * hs0_1[k]
                  - f_2 * hs1_1[k]
                  + pb_z[k] * hp_3[k];

        t_12[k] = pa_z[k] * gd_4[k];

        t_13[k] = pa_y[k] * gd_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_y, fd0_0, fd0_12, fd1_0, fd1_10, \
                         gd_5, gd_12, hs0_2, hs1_2, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_4 * fd0_0[k]
                  - f_5 * fd1_0[k]
                  + pa_z[k] * gd_5[k];

        t_15[k] = f_1 * hs0_2[k]
                  - f_2 * hs1_2[k]
                  + pb_y[k] * hp_4[k];

        t_16[k] = f_6 * fd0_12[k]
                  - f_7 * fd1_10[k]
                  + pa_x[k] * gd_12[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, pb_z, fd0_3, fd0_14, fd1_3, fd1_12, \
                         gd_7, gd_14, hs0_3, hs1_3, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * fd0_3[k]
                  - f_7 * fd1_3[k]
                  + pa_y[k] * gd_7[k];

        t_18[k] = f_4 * fd0_14[k]
                  - f_5 * fd1_12[k]
                  + pa_x[k] * gd_14[k];

        t_19[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pa_z, gp_3, gp_4, gd_7, gd_8, \
                         gd_9, gd_11, gd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * gd_7[k];

        t_21[k] = pa_z[k] * gd_8[k];

        t_22[k] = f_3 * gp_3[k]
                  + pa_z[k] * gd_9[k];

        t_23[k] = f_3 * gp_4[k]
                  + pa_y[k] * gd_11[k];

        t_24[k] = pa_y[k] * gd_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_x, pa_z, pb_y, fd0_6, fd0_23, fd1_5, fd1_20, \
                         gd_10, gd_16, hs0_4, hs1_4, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_6 * fd0_6[k]
                  - f_7 * fd1_5[k]
                  + pa_z[k] * gd_10[k];

        t_26[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];

        t_27[k] = f_4 * fd0_23[k]
                  - f_5 * fd1_20[k]
                  + pa_x[k] * gd_16[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_x, pa_z, gp_5, gp_8, gp_9, gd_13, \
                         gd_17, gd_18, gd_22, gd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * gp_5[k]
                  + pa_x[k] * gd_17[k];

        t_29[k] = pa_x[k] * gd_18[k];

        t_30[k] = pa_z[k] * gd_13[k];

        t_31[k] = f_3 * gp_8[k]
                  + pa_x[k] * gd_22[k];

        t_32[k] = f_3 * gp_9[k]
                  + pa_x[k] * gd_27[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_x, pb_x, pb_y, pb_z, gp_6, gd_29, hs0_5, \
                         hs1_5, hp_7, hp_8, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_x[k] * gd_29[k];

        t_34[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_x[k] * hp_7[k];

        t_35[k] = f_0 * gp_6[k]
                  + f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_y[k] * hp_8[k];

        t_36[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_z[k] * hp_9[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, gp_7, gd_17, gd_18, gd_19, hs0_6, \
                         hs1_6, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pa_z[k] * gd_17[k];

        t_38[k] = pa_z[k] * gd_18[k];

        t_39[k] = f_3 * gp_7[k]
                  + pa_z[k] * gd_19[k];

        t_40[k] = f_1 * hs0_6[k]
                  - f_2 * hs1_6[k]
                  + pb_x[k] * hp_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_y, pa_z, pb_x, fd0_14, fd0_20, fd1_12, fd1_17, \
                         gd_20, gd_24, hs0_7, hs1_7, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_4 * fd0_14[k]
                  - f_5 * fd1_12[k]
                  + pa_z[k] * gd_20[k];

        t_42[k] = f_6 * fd0_20[k]
                  - f_7 * fd1_17[k]
                  + pa_y[k] * gd_24[k];

        t_43[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_11[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_y, pa_z, fd0_17, fd0_23, fd1_14, fd1_20, \
                         gp_10, gd_23, gd_26, gd_28, gd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_6 * fd0_17[k]
                  - f_7 * fd1_14[k]
                  + pa_z[k] * gd_23[k];

        t_45[k] = f_4 * fd0_23[k]
                  - f_5 * fd1_20[k]
                  + pa_y[k] * gd_26[k];

        t_46[k] = f_3 * gp_10[k]
                  + pa_y[k] * gd_28[k];

        t_47[k] = pa_y[k] * gd_29[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_x, pb_y, pb_z, gp_11, hs0_8, hs1_8, hp_12, \
                         hp_13, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_12[k];

        t_49[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_y[k] * hp_13[k];

        t_50[k] = f_0 * gp_11[k]
                  + f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fd0, const size_t fd1,
                                     const size_t gp, const size_t gd, const size_t hs0,
                                     const size_t hs1, const size_t hp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / alpha;
    const auto f_6 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_4 = buffer.data(fd1 + 4);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_11 = buffer.data(fd1 + 11);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_11 = buffer.data(gp + 11);

    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs0_0, hs1_0, hp_0, hp_1, \
                         hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, fd0_0, fd0_3, fd1_0, fd1_5, gd_3, \
                         gd_6, hs0_1, hs1_1, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_y[k] * gd_3[k];

        t_4[k] = f_5 * fd0_3[k]
                 - f_6 * fd1_5[k]
                 + pa_x[k] * gd_6[k];

        t_5[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_1[k]
                 + pb_z[k] * hp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, fd0_0, fd0_4, fd1_0, fd1_6, gd_4, \
                         gd_10, hs0_2, hs1_2, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_z[k] * gd_4[k];

        t_7[k] = f_1 * hs0_2[k]
                 - f_2 * hs1_2[k]
                 + pb_y[k] * hp_4[k];

        t_8[k] = f_5 * fd0_4[k]
                 - f_6 * fd1_6[k]
                 + pa_x[k] * gd_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, fd0_1, fd0_5, fd1_3, fd1_8, gd_5, \
                         gd_11, hs0_3, hs1_3, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * fd0_1[k]
                 - f_6 * fd1_3[k]
                 + pa_y[k] * gd_5[k];

        t_10[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_8[k]
                  + pa_x[k] * gd_11[k];

        t_11[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, fd0_2, fd0_8, fd1_4, fd1_14, \
                         gd_8, gd_12, hs0_4, hs1_4, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * fd0_2[k]
                  - f_6 * fd1_4[k]
                  + pa_z[k] * gd_8[k];

        t_13[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];

        t_14[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_14[k]
                  + pa_x[k] * gd_12[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, gp_6, hs0_5, hs0_6, hs1_5, \
                         hs1_6, hp_7, hp_8, hp_9, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_x[k] * hp_7[k];

        t_16[k] = f_0 * gp_6[k]
                  + f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_y[k] * hp_8[k];

        t_17[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_z[k] * hp_9[k];

        t_18[k] = f_1 * hs0_6[k]
                  - f_2 * hs1_6[k]
                  + pb_x[k] * hp_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pa_z, pb_x, fd0_5, fd0_7, fd1_8, fd1_11, \
                         gd_16, gd_19, hs0_7, hs1_7, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_8[k]
                  + pa_z[k] * gd_16[k];

        t_20[k] = f_5 * fd0_7[k]
                  - f_6 * fd1_11[k]
                  + pa_y[k] * gd_19[k];

        t_21[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pa_z, pb_x, fd0_6, fd0_8, fd1_10, fd1_14, \
                         gd_18, gd_20, hs0_8, hs1_8, hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * fd0_6[k]
                  - f_6 * fd1_10[k]
                  + pa_z[k] * gd_18[k];

        t_23[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_14[k]
                  + pa_y[k] * gd_20[k];

        t_24[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_12[k];
    }

#pragma omp simd aligned(t_25, t_26, pb_y, pb_z, gp_11, hs0_8, hs1_8, hp_13, \
                         hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_y[k] * hp_13[k];

        t_26[k] = f_0 * gp_11[k]
                  + f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fd0, const size_t fd1,
                                     const size_t gp, const size_t gd, const size_t hs0,
                                     const size_t hs1, const size_t hp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / alpha;
    const auto f_6 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_3 = buffer.data(fd0 + 3);
    const auto *fd0_4 = buffer.data(fd0 + 4);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_6 = buffer.data(fd0 + 6);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_10 = buffer.data(fd0 + 10);
    const auto *fd0_11 = buffer.data(fd0 + 11);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_4 = buffer.data(fd1 + 4);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_7 = buffer.data(fd1 + 7);
    const auto *fd1_9 = buffer.data(fd1 + 9);
    const auto *fd1_11 = buffer.data(fd1 + 11);
    const auto *fd1_14 = buffer.data(fd1 + 14);
    const auto *fd1_17 = buffer.data(fd1 + 17);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_11 = buffer.data(gp + 11);

    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs0_0, hs1_0, hp_0, hp_1, \
                         hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, fd0_0, fd0_5, fd1_0, fd1_6, gd_3, \
                         gd_7, hs0_1, hs1_1, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_y[k] * gd_3[k];

        t_4[k] = f_5 * fd0_5[k]
                 - f_6 * fd1_6[k]
                 + pa_x[k] * gd_7[k];

        t_5[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_1[k]
                 + pb_z[k] * hp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, fd0_0, fd0_6, fd1_0, fd1_7, gd_4, \
                         gd_11, hs0_2, hs1_2, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_z[k] * gd_4[k];

        t_7[k] = f_1 * hs0_2[k]
                 - f_2 * hs1_2[k]
                 + pb_y[k] * hp_4[k];

        t_8[k] = f_5 * fd0_6[k]
                 - f_6 * fd1_7[k]
                 + pa_x[k] * gd_11[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, fd0_3, fd0_8, fd1_3, fd1_9, gd_6, \
                         gd_12, hs0_3, hs1_3, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * fd0_3[k]
                 - f_6 * fd1_3[k]
                 + pa_y[k] * gd_6[k];

        t_10[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_9[k]
                  + pa_x[k] * gd_12[k];

        t_11[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, fd0_4, fd0_14, fd1_4, fd1_17, \
                         gd_9, gd_13, hs0_4, hs1_4, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * fd0_4[k]
                  - f_6 * fd1_4[k]
                  + pa_z[k] * gd_9[k];

        t_13[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];

        t_14[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_17[k]
                  + pa_x[k] * gd_13[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, gp_6, hs0_5, hs0_6, hs1_5, \
                         hs1_6, hp_7, hp_8, hp_9, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_x[k] * hp_7[k];

        t_16[k] = f_0 * gp_6[k]
                  + f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_y[k] * hp_8[k];

        t_17[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_z[k] * hp_9[k];

        t_18[k] = f_1 * hs0_6[k]
                  - f_2 * hs1_6[k]
                  + pb_x[k] * hp_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pa_z, pb_x, fd0_8, fd0_11, fd1_9, fd1_14, \
                         gd_17, gd_21, hs0_7, hs1_7, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_9[k]
                  + pa_z[k] * gd_17[k];

        t_20[k] = f_5 * fd0_11[k]
                  - f_6 * fd1_14[k]
                  + pa_y[k] * gd_21[k];

        t_21[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pa_z, pb_x, fd0_10, fd0_14, fd1_11, fd1_17, \
                         gd_20, gd_23, hs0_8, hs1_8, hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * fd0_10[k]
                  - f_6 * fd1_11[k]
                  + pa_z[k] * gd_20[k];

        t_23[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_17[k]
                  + pa_y[k] * gd_23[k];

        t_24[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_12[k];
    }

#pragma omp simd aligned(t_25, t_26, pb_y, pb_z, gp_11, hs0_8, hs1_8, hp_13, \
                         hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_y[k] * hp_13[k];

        t_26[k] = f_0 * gp_11[k]
                  + f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fd0, const size_t fd1,
                                     const size_t gp, const size_t gd, const size_t hs0,
                                     const size_t hs1, const size_t hp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / alpha;
    const auto f_6 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_3 = buffer.data(fd0 + 3);
    const auto *fd0_4 = buffer.data(fd0 + 4);
    const auto *fd0_6 = buffer.data(fd0 + 6);
    const auto *fd0_7 = buffer.data(fd0 + 7);
    const auto *fd0_9 = buffer.data(fd0 + 9);
    const auto *fd0_11 = buffer.data(fd0 + 11);
    const auto *fd0_14 = buffer.data(fd0 + 14);
    const auto *fd0_17 = buffer.data(fd0 + 17);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_4 = buffer.data(fd1 + 4);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_11 = buffer.data(fd1 + 11);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_11 = buffer.data(gp + 11);

    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs0_0, hs1_0, hp_0, hp_1, \
                         hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, fd0_0, fd0_6, fd1_0, fd1_5, gd_3, \
                         gd_6, hs0_1, hs1_1, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_y[k] * gd_3[k];

        t_4[k] = f_5 * fd0_6[k]
                 - f_6 * fd1_5[k]
                 + pa_x[k] * gd_6[k];

        t_5[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_1[k]
                 + pb_z[k] * hp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, fd0_0, fd0_7, fd1_0, fd1_6, gd_4, \
                         gd_10, hs0_2, hs1_2, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_z[k] * gd_4[k];

        t_7[k] = f_1 * hs0_2[k]
                 - f_2 * hs1_2[k]
                 + pb_y[k] * hp_4[k];

        t_8[k] = f_5 * fd0_7[k]
                 - f_6 * fd1_6[k]
                 + pa_x[k] * gd_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, fd0_3, fd0_9, fd1_3, fd1_8, gd_5, \
                         gd_11, hs0_3, hs1_3, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * fd0_3[k]
                 - f_6 * fd1_3[k]
                 + pa_y[k] * gd_5[k];

        t_10[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_8[k]
                  + pa_x[k] * gd_11[k];

        t_11[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, fd0_4, fd0_17, fd1_4, fd1_14, \
                         gd_8, gd_12, hs0_4, hs1_4, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * fd0_4[k]
                  - f_6 * fd1_4[k]
                  + pa_z[k] * gd_8[k];

        t_13[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];

        t_14[k] = f_3 * fd0_17[k]
                  - f_4 * fd1_14[k]
                  + pa_x[k] * gd_12[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, gp_6, hs0_5, hs0_6, hs1_5, \
                         hs1_6, hp_7, hp_8, hp_9, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_x[k] * hp_7[k];

        t_16[k] = f_0 * gp_6[k]
                  + f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_y[k] * hp_8[k];

        t_17[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_z[k] * hp_9[k];

        t_18[k] = f_1 * hs0_6[k]
                  - f_2 * hs1_6[k]
                  + pb_x[k] * hp_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pa_z, pb_x, fd0_9, fd0_14, fd1_8, fd1_11, \
                         gd_16, gd_19, hs0_7, hs1_7, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * fd0_9[k]
                  - f_4 * fd1_8[k]
                  + pa_z[k] * gd_16[k];

        t_20[k] = f_5 * fd0_14[k]
                  - f_6 * fd1_11[k]
                  + pa_y[k] * gd_19[k];

        t_21[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pa_z, pb_x, fd0_11, fd0_17, fd1_10, fd1_14, \
                         gd_18, gd_20, hs0_8, hs1_8, hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * fd0_11[k]
                  - f_6 * fd1_10[k]
                  + pa_z[k] * gd_18[k];

        t_23[k] = f_3 * fd0_17[k]
                  - f_4 * fd1_14[k]
                  + pa_y[k] * gd_20[k];

        t_24[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_12[k];
    }

#pragma omp simd aligned(t_25, t_26, pb_y, pb_z, gp_11, hs0_8, hs1_8, hp_13, \
                         hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_y[k] * hp_13[k];

        t_26[k] = f_0 * gp_11[k]
                  + f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fd0, const size_t fd1,
                                     const size_t gp, const size_t gd, const size_t hs0,
                                     const size_t hs1, const size_t hp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / alpha;
    const auto f_6 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_3 = buffer.data(fd0 + 3);
    const auto *fd0_4 = buffer.data(fd0 + 4);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_6 = buffer.data(fd0 + 6);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_10 = buffer.data(fd0 + 10);
    const auto *fd0_11 = buffer.data(fd0 + 11);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_4 = buffer.data(fd1 + 4);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_11 = buffer.data(fd1 + 11);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_11 = buffer.data(gp + 11);

    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs0_0, hs1_0, hp_0, hp_1, \
                         hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, fd0_0, fd0_5, fd1_0, fd1_5, gd_3, \
                         gd_6, hs0_1, hs1_1, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_y[k] * gd_3[k];

        t_4[k] = f_5 * fd0_5[k]
                 - f_6 * fd1_5[k]
                 + pa_x[k] * gd_6[k];

        t_5[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_1[k]
                 + pb_z[k] * hp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, fd0_0, fd0_6, fd1_0, fd1_6, gd_4, \
                         gd_10, hs0_2, hs1_2, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_z[k] * gd_4[k];

        t_7[k] = f_1 * hs0_2[k]
                 - f_2 * hs1_2[k]
                 + pb_y[k] * hp_4[k];

        t_8[k] = f_5 * fd0_6[k]
                 - f_6 * fd1_6[k]
                 + pa_x[k] * gd_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, fd0_3, fd0_8, fd1_3, fd1_8, gd_5, \
                         gd_11, hs0_3, hs1_3, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * fd0_3[k]
                 - f_6 * fd1_3[k]
                 + pa_y[k] * gd_5[k];

        t_10[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pa_x[k] * gd_11[k];

        t_11[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, fd0_4, fd0_14, fd1_4, fd1_14, \
                         gd_8, gd_12, hs0_4, hs1_4, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * fd0_4[k]
                  - f_6 * fd1_4[k]
                  + pa_z[k] * gd_8[k];

        t_13[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];

        t_14[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pa_x[k] * gd_12[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, gp_6, hs0_5, hs0_6, hs1_5, \
                         hs1_6, hp_7, hp_8, hp_9, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_x[k] * hp_7[k];

        t_16[k] = f_0 * gp_6[k]
                  + f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_y[k] * hp_8[k];

        t_17[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_5[k]
                  + pb_z[k] * hp_9[k];

        t_18[k] = f_1 * hs0_6[k]
                  - f_2 * hs1_6[k]
                  + pb_x[k] * hp_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pa_z, pb_x, fd0_8, fd0_11, fd1_8, fd1_11, \
                         gd_16, gd_19, hs0_7, hs1_7, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pa_z[k] * gd_16[k];

        t_20[k] = f_5 * fd0_11[k]
                  - f_6 * fd1_11[k]
                  + pa_y[k] * gd_19[k];

        t_21[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pa_z, pb_x, fd0_10, fd0_14, fd1_10, fd1_14, \
                         gd_18, gd_20, hs0_8, hs1_8, hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * fd0_10[k]
                  - f_6 * fd1_10[k]
                  + pa_z[k] * gd_18[k];

        t_23[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pa_y[k] * gd_20[k];

        t_24[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_12[k];
    }

#pragma omp simd aligned(t_25, t_26, pb_y, pb_z, gp_11, hs0_8, hs1_8, hp_13, \
                         hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_y[k] * hp_13[k];

        t_26[k] = f_0 * gp_11[k]
                  + f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fd0, const size_t fd1,
                                     const size_t gp, const size_t gd, const size_t hs0,
                                     const size_t hs1, const size_t hp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 2.0 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 1.5 / p;
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);

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

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_11 = buffer.data(hs1 + 11);
    const auto *hs1_13 = buffer.data(hs1 + 13);
    const auto *hs1_14 = buffer.data(hs1 + 14);
    const auto *hs1_16 = buffer.data(hs1 + 16);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, gp_0, gd_0, hs0_0, hs1_0, \
                         hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];

        t_3[k] = pa_y[k] * gd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, pb_x, gp_1, gp_3, gp_4, gd_0, gd_1, \
                         hp_3, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * gp_3[k]
                 + pb_x[k] * hp_3[k];

        t_5[k] = f_4 * gp_1[k]
                 + pa_y[k] * gd_1[k];

        t_6[k] = pa_z[k] * gd_0[k];

        t_7[k] = f_3 * gp_4[k]
                 + pb_x[k] * hp_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_y, pa_z, pb_x, fd0_0, fd1_0, gp_2, gp_5, gd_2, \
                         gd_3, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * gp_2[k]
                 + pa_z[k] * gd_2[k];

        t_9[k] = f_5 * fd0_0[k]
                 - f_6 * fd1_0[k]
                 + pa_y[k] * gd_3[k];

        t_10[k] = f_7 * gp_5[k]
                  + pb_x[k] * hp_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_z, pb_z, fd0_0, fd0_3, fd1_0, fd1_3, gd_4, \
                         gd_6, hs0_1, hs1_3, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_8 * fd0_3[k]
                  - f_9 * fd1_3[k]
                  + pa_x[k] * gd_6[k];

        t_12[k] = f_1 * hs0_1[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_6[k];

        t_13[k] = f_5 * fd0_0[k]
                  - f_6 * fd1_0[k]
                  + pa_z[k] * gd_4[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pb_x, pb_y, fd0_4, fd1_4, gp_6, gd_8, hs0_2, \
                         hs1_4, hp_7, hp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_7 * gp_6[k]
                  + pb_x[k] * hp_8[k];

        t_15[k] = f_1 * hs0_2[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_7[k];

        t_16[k] = f_8 * fd0_4[k]
                  - f_9 * fd1_4[k]
                  + pa_x[k] * gd_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, pb_x, fd0_1, fd0_5, fd1_1, fd1_5, gp_7, \
                         gd_5, gd_9, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * fd0_1[k]
                  - f_9 * fd1_1[k]
                  + pa_y[k] * gd_5[k];

        t_18[k] = f_4 * gp_7[k]
                  + pb_x[k] * hp_9[k];

        t_19[k] = f_5 * fd0_5[k]
                  - f_6 * fd1_5[k]
                  + pa_x[k] * gd_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_y, pa_z, pb_x, pb_z, fd0_2, fd1_2, gp_8, \
                         gd_7, hs0_3, hs1_5, hp_10, hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_5[k]
                  + pb_z[k] * hp_10[k];

        t_21[k] = pa_y[k] * gd_7[k];

        t_22[k] = f_8 * fd0_2[k]
                  - f_9 * fd1_2[k]
                  + pa_z[k] * gd_7[k];

        t_23[k] = f_4 * gp_8[k]
                  + pb_x[k] * hp_12[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pb_y, fd0_8, fd1_8, gp_9, gd_10, gd_11, \
                         hs0_4, hs1_8, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_8[k]
                  + pb_y[k] * hp_11[k];

        t_25[k] = f_5 * fd0_8[k]
                  - f_6 * fd1_8[k]
                  + pa_x[k] * gd_10[k];

        t_26[k] = f_4 * gp_9[k]
                  + pa_x[k] * gd_11[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_x, pb_x, gp_10, gp_17, gd_12, gd_15, \
                         gd_16, gd_18, hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_10 * gp_10[k]
                  + pb_x[k] * hp_13[k];

        t_28[k] = pa_x[k] * gd_12[k];

        t_29[k] = pa_x[k] * gd_15[k];

        t_30[k] = pa_x[k] * gd_16[k];

        t_31[k] = f_4 * gp_17[k]
                  + pa_x[k] * gd_18[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_y, gp_10, gp_19, gd_20, hs0_7, \
                         hs1_11, hp_14, hp_15, hp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_10 * gp_19[k]
                  + pb_x[k] * hp_14[k];

        t_33[k] = pa_x[k] * gd_20[k];

        t_34[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_11[k]
                  + pb_x[k] * hp_15[k];

        t_35[k] = f_0 * gp_10[k]
                  + f_1 * hs0_7[k]
                  - f_2 * hs1_11[k]
                  + pb_y[k] * hp_16[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_y, pb_z, gp_11, gp_12, gd_12, gd_13, \
                         hs0_7, hs1_11, hp_17, hp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_11[k]
                  + pb_z[k] * hp_17[k];

        t_37[k] = pa_z[k] * gd_12[k];

        t_38[k] = f_3 * gp_12[k]
                  + pb_y[k] * hp_18[k];

        t_39[k] = f_4 * gp_11[k]
                  + pa_z[k] * gd_13[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_z, pb_x, pb_y, fd0_5, fd1_5, gp_14, gd_14, \
                         hs0_8, hs1_13, hp_19, hp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_13[k]
                  + pb_x[k] * hp_19[k];

        t_41[k] = f_5 * fd0_5[k]
                  - f_6 * fd1_5[k]
                  + pa_z[k] * gd_14[k];

        t_42[k] = f_7 * gp_14[k]
                  + pb_y[k] * hp_21[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_y, pa_z, pb_x, fd0_6, fd0_7, fd1_6, fd1_7, \
                         gd_15, gd_16, hs0_9, hs1_14, hp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_8 * fd0_7[k]
                  - f_9 * fd1_7[k]
                  + pa_y[k] * gd_16[k];

        t_44[k] = f_1 * hs0_9[k]
                  - f_2 * hs1_14[k]
                  + pb_x[k] * hp_22[k];

        t_45[k] = f_8 * fd0_6[k]
                  - f_9 * fd1_6[k]
                  + pa_z[k] * gd_15[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pb_y, fd0_8, fd1_8, gp_16, gp_18, \
                         gp_19, gd_17, gd_19, hp_24, hp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_4 * gp_16[k]
                  + pb_y[k] * hp_24[k];

        t_47[k] = f_5 * fd0_8[k]
                  - f_6 * fd1_8[k]
                  + pa_y[k] * gd_17[k];

        t_48[k] = f_4 * gp_18[k]
                  + pa_y[k] * gd_19[k];

        t_49[k] = f_10 * gp_19[k]
                  + pb_y[k] * hp_26[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pb_x, pb_y, pb_z, gp_19, gd_20, hs0_11, \
                         hs1_16, hp_27, hp_28, hp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * gd_20[k];

        t_51[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_16[k]
                  + pb_x[k] * hp_27[k];

        t_52[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_16[k]
                  + pb_y[k] * hp_28[k];

        t_53[k] = f_0 * gp_19[k]
                  + f_1 * hs0_11[k]
                  - f_2 * hs1_16[k]
                  + pb_z[k] * hp_29[k];
    }
}

auto
compute_prim_hd_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t fd0, const size_t fd1,
                                      const size_t gp, const size_t gd, const size_t hs0,
                                      const size_t hs1, const size_t hp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / p;
    const auto f_9 = 2.0 / p;
    const auto f_10 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);

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

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_11 = buffer.data(hs1 + 11);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, gp_0, gd_0, hs0_0, hs1_0, \
                         hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];

        t_3[k] = pa_y[k] * gd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, fd0_0, fd1_0, gp_1, gp_2, \
                         gd_0, gd_1, gd_2, gd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * gp_1[k]
                 + pa_y[k] * gd_1[k];

        t_5[k] = pa_y[k] * gd_2[k];

        t_6[k] = pa_z[k] * gd_0[k];

        t_7[k] = pa_z[k] * gd_1[k];

        t_8[k] = f_3 * gp_2[k]
                 + pa_z[k] * gd_2[k];

        t_9[k] = f_4 * fd0_0[k]
                 - f_5 * fd1_0[k]
                 + pa_y[k] * gd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_z, pb_z, fd0_3, fd1_3, gd_4, gd_8, hs0_1, \
                         hs1_1, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_6 * fd0_3[k]
                  - f_7 * fd1_3[k]
                  + pa_x[k] * gd_8[k];

        t_11[k] = f_1 * hs0_1[k]
                  - f_2 * hs1_1[k]
                  + pb_z[k] * hp_6[k];

        t_12[k] = pa_z[k] * gd_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pa_z, pb_y, fd0_0, fd1_0, gp_3, gd_5, \
                         gd_6, hs0_2, hs1_2, hp_7, hp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_8 * gp_3[k]
                  + pb_y[k] * hp_7[k];

        t_14[k] = pa_y[k] * gd_6[k];

        t_15[k] = f_4 * fd0_0[k]
                  - f_5 * fd1_0[k]
                  + pa_z[k] * gd_5[k];

        t_16[k] = f_1 * hs0_2[k]
                  - f_2 * hs1_2[k]
                  + pb_y[k] * hp_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, fd0_1, fd0_4, fd0_5, fd1_1, fd1_4, \
                         fd1_5, gd_7, gd_12, gd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * fd0_4[k]
                  - f_7 * fd1_4[k]
                  + pa_x[k] * gd_12[k];

        t_18[k] = f_6 * fd0_1[k]
                  - f_7 * fd1_1[k]
                  + pa_y[k] * gd_7[k];

        t_19[k] = f_4 * fd0_5[k]
                  - f_5 * fd1_5[k]
                  + pa_x[k] * gd_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_z, pb_y, pb_z, gp_5, gd_7, gd_8, hs0_3, \
                         hs1_3, hp_11, hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_11[k];

        t_21[k] = pa_z[k] * gd_7[k];

        t_22[k] = pa_z[k] * gd_8[k];

        t_23[k] = f_3 * gp_5[k]
                  + pb_y[k] * hp_12[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_y, pa_z, pb_y, gp_4, gp_6, gp_7, \
                         gd_9, gd_10, gd_11, gd_12, hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_3 * gp_4[k]
                  + pa_z[k] * gd_9[k];

        t_25[k] = pa_y[k] * gd_10[k];

        t_26[k] = f_3 * gp_6[k]
                  + pa_y[k] * gd_11[k];

        t_27[k] = f_8 * gp_7[k]
                  + pb_y[k] * hp_13[k];

        t_28[k] = pa_y[k] * gd_12[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_x, pa_z, pb_y, fd0_2, fd0_8, fd1_2, fd1_8, \
                         gd_10, gd_16, hs0_4, hs1_4, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_6 * fd0_2[k]
                  - f_7 * fd1_2[k]
                  + pa_z[k] * gd_10[k];

        t_30[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_14[k];

        t_31[k] = f_4 * fd0_8[k]
                  - f_5 * fd1_8[k]
                  + pa_x[k] * gd_16[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, pa_x, pa_z, gp_8, gd_13, gd_17, \
                         gd_18, gd_19, gd_21, gd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * gp_8[k]
                  + pa_x[k] * gd_17[k];

        t_33[k] = pa_x[k] * gd_18[k];

        t_34[k] = pa_x[k] * gd_19[k];

        t_35[k] = pa_z[k] * gd_13[k];

        t_36[k] = pa_x[k] * gd_21[k];

        t_37[k] = pa_x[k] * gd_22[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, pa_x, pa_y, gp_12, gd_15, gd_23, \
                         gd_24, gd_25, gd_26, gd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_3 * gp_12[k]
                  + pa_x[k] * gd_23[k];

        t_39[k] = pa_x[k] * gd_24[k];

        t_40[k] = pa_x[k] * gd_25[k];

        t_41[k] = pa_x[k] * gd_26[k];

        t_42[k] = pa_y[k] * gd_15[k];

        t_43[k] = pa_x[k] * gd_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pa_x, pb_x, gp_15, gd_28, gd_30, gd_31, \
                         gd_32, hs0_5, hs1_7, hp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_x[k] * gd_28[k];

        t_45[k] = f_3 * gp_15[k]
                  + pa_x[k] * gd_30[k];

        t_46[k] = pa_x[k] * gd_31[k];

        t_47[k] = pa_x[k] * gd_32[k];

        t_48[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_18[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_z, pb_y, pb_z, gp_9, gd_17, gd_18, hs0_5, \
                         hs1_7, hp_19, hp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_0 * gp_9[k]
                  + f_1 * hs0_5[k]
                  - f_2 * hs1_7[k]
                  + pb_y[k] * hp_19[k];

        t_50[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_7[k]
                  + pb_z[k] * hp_20[k];

        t_51[k] = pa_z[k] * gd_17[k];

        t_52[k] = pa_z[k] * gd_18[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_z, pb_x, pb_y, gp_10, gp_11, gd_19, hs0_6, \
                         hs1_8, hp_21, hp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_9 * gp_11[k]
                  + pb_y[k] * hp_21[k];

        t_54[k] = f_3 * gp_10[k]
                  + pa_z[k] * gd_19[k];

        t_55[k] = f_1 * hs0_6[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_22[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pa_y, pa_z, pb_y, fd0_5, fd0_7, fd1_5, fd1_7, \
                         gp_13, gd_20, gd_26, hp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_4 * fd0_5[k]
                  - f_5 * fd1_5[k]
                  + pa_z[k] * gd_20[k];

        t_57[k] = f_10 * gp_13[k]
                  + pb_y[k] * hp_23[k];

        t_58[k] = f_6 * fd0_7[k]
                  - f_7 * fd1_7[k]
                  + pa_y[k] * gd_26[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pa_z, pb_x, pb_y, fd0_6, fd1_6, gp_14, gd_24, \
                         hs0_7, hs1_9, hp_24, hp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_9[k]
                  + pb_x[k] * hp_24[k];

        t_60[k] = f_6 * fd0_6[k]
                  - f_7 * fd1_6[k]
                  + pa_z[k] * gd_24[k];

        t_61[k] = f_3 * gp_14[k]
                  + pb_y[k] * hp_25[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, pa_y, pb_y, fd0_8, fd1_8, gp_16, gp_17, \
                         gd_29, gd_30, gd_31, gd_32, hp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_4 * fd0_8[k]
                  - f_5 * fd1_8[k]
                  + pa_y[k] * gd_29[k];

        t_63[k] = pa_y[k] * gd_30[k];

        t_64[k] = f_3 * gp_16[k]
                  + pa_y[k] * gd_31[k];

        t_65[k] = f_8 * gp_17[k]
                  + pb_y[k] * hp_26[k];

        t_66[k] = pa_y[k] * gd_32[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pb_x, pb_y, pb_z, gp_17, hs0_8, hs1_11, hp_27, \
                         hp_28, hp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_11[k]
                  + pb_x[k] * hp_27[k];

        t_68[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_11[k]
                  + pb_y[k] * hp_28[k];

        t_69[k] = f_0 * gp_17[k]
                  + f_1 * hs0_8[k]
                  - f_2 * hs1_11[k]
                  + pb_z[k] * hp_29[k];
    }
}

auto
compute_prim_hd_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t fd0, const size_t fd1,
                                      const size_t gp, const size_t gd, const size_t hs0,
                                      const size_t hs1, const size_t hp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 2.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_10 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_4 = buffer.data(fd1 + 4);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_11 = buffer.data(fd1 + 11);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);

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

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_12 = buffer.data(hs1 + 12);
    const auto *hs1_13 = buffer.data(hs1 + 13);
    const auto *hs1_15 = buffer.data(hs1 + 15);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_21 = buffer.data(hp + 21);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, gp_0, gd_0, hs0_0, hs1_0, \
                         hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];

        t_3[k] = pa_y[k] * gd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, fd0_0, fd1_0, gp_1, gp_2, gd_0, gd_1, \
                         gd_2, gd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * gp_1[k]
                 + pa_y[k] * gd_1[k];

        t_5[k] = pa_z[k] * gd_0[k];

        t_6[k] = f_3 * gp_2[k]
                 + pa_z[k] * gd_2[k];

        t_7[k] = f_4 * fd0_0[k]
                 - f_5 * fd1_0[k]
                 + pa_y[k] * gd_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_x, pa_z, pb_z, fd0_0, fd0_3, fd1_0, fd1_5, gd_4, \
                         gd_6, hs0_1, hs1_3, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_6 * fd0_3[k]
                 - f_7 * fd1_5[k]
                 + pa_x[k] * gd_6[k];

        t_9[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_3[k]
                 + pb_z[k] * hp_4[k];

        t_10[k] = f_4 * fd0_0[k]
                  - f_5 * fd1_0[k]
                  + pa_z[k] * gd_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, pb_y, fd0_1, fd0_4, fd1_3, fd1_6, gd_5, \
                         gd_8, hs0_2, hs1_4, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_1 * hs0_2[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_5[k];

        t_12[k] = f_6 * fd0_4[k]
                  - f_7 * fd1_6[k]
                  + pa_x[k] * gd_8[k];

        t_13[k] = f_6 * fd0_1[k]
                  - f_7 * fd1_3[k]
                  + pa_y[k] * gd_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_y, pb_z, fd0_5, fd1_8, gd_7, gd_9, hs0_3, \
                         hs1_5, hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_4 * fd0_5[k]
                  - f_5 * fd1_8[k]
                  + pa_x[k] * gd_9[k];

        t_15[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_5[k]
                  + pb_z[k] * hp_7[k];

        t_16[k] = pa_y[k] * gd_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_z, pb_y, fd0_2, fd0_8, fd1_4, fd1_14, \
                         gd_7, gd_10, hs0_4, hs1_7, hp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * fd0_2[k]
                  - f_7 * fd1_4[k]
                  + pa_z[k] * gd_7[k];

        t_18[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_7[k]
                  + pb_y[k] * hp_8[k];

        t_19[k] = f_4 * fd0_8[k]
                  - f_5 * fd1_14[k]
                  + pa_x[k] * gd_10[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, pa_x, gp_7, gp_14, gd_11, gd_12, \
                         gd_15, gd_16, gd_18, gd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_3 * gp_7[k]
                  + pa_x[k] * gd_11[k];

        t_21[k] = pa_x[k] * gd_12[k];

        t_22[k] = pa_x[k] * gd_15[k];

        t_23[k] = pa_x[k] * gd_16[k];

        t_24[k] = f_3 * gp_14[k]
                  + pa_x[k] * gd_18[k];

        t_25[k] = pa_x[k] * gd_20[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, pb_z, gp_8, gd_12, hs0_7, \
                         hs1_10, hp_10, hp_11, hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_10[k]
                  + pb_x[k] * hp_10[k];

        t_27[k] = f_0 * gp_8[k]
                  + f_1 * hs0_7[k]
                  - f_2 * hs1_10[k]
                  + pb_y[k] * hp_11[k];

        t_28[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_10[k]
                  + pb_z[k] * hp_12[k];

        t_29[k] = pa_z[k] * gd_12[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_z, pb_x, pb_y, gp_9, gp_10, gd_13, hs0_8, \
                         hs1_12, hp_13, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_8 * gp_10[k]
                  + pb_y[k] * hp_13[k];

        t_31[k] = f_3 * gp_9[k]
                  + pa_z[k] * gd_13[k];

        t_32[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_12[k]
                  + pb_x[k] * hp_14[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pa_z, pb_y, fd0_5, fd0_7, fd1_8, fd1_11, \
                         gp_12, gd_14, gd_16, hp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_4 * fd0_5[k]
                  - f_5 * fd1_8[k]
                  + pa_z[k] * gd_14[k];

        t_34[k] = f_9 * gp_12[k]
                  + pb_y[k] * hp_15[k];

        t_35[k] = f_6 * fd0_7[k]
                  - f_7 * fd1_11[k]
                  + pa_y[k] * gd_16[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_z, pb_x, pb_y, fd0_6, fd1_10, gp_13, gd_15, \
                         hs0_9, hs1_13, hp_16, hp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * hs0_9[k]
                  - f_2 * hs1_13[k]
                  + pb_x[k] * hp_16[k];

        t_37[k] = f_6 * fd0_6[k]
                  - f_7 * fd1_10[k]
                  + pa_z[k] * gd_15[k];

        t_38[k] = f_3 * gp_13[k]
                  + pb_y[k] * hp_17[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_y, pb_y, fd0_8, fd1_14, gp_15, gp_16, \
                         gd_17, gd_19, gd_20, hp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_4 * fd0_8[k]
                  - f_5 * fd1_14[k]
                  + pa_y[k] * gd_17[k];

        t_40[k] = f_3 * gp_15[k]
                  + pa_y[k] * gd_19[k];

        t_41[k] = f_10 * gp_16[k]
                  + pb_y[k] * hp_18[k];

        t_42[k] = pa_y[k] * gd_20[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pb_x, pb_y, pb_z, gp_16, hs0_11, hs1_15, hp_19, \
                         hp_20, hp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_15[k]
                  + pb_x[k] * hp_19[k];

        t_44[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_15[k]
                  + pb_y[k] * hp_20[k];

        t_45[k] = f_0 * gp_16[k]
                  + f_1 * hs0_11[k]
                  - f_2 * hs1_15[k]
                  + pb_z[k] * hp_21[k];
    }
}

auto
compute_prim_hd_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t fd0, const size_t fd1,
                                      const size_t gp, const size_t gd, const size_t hs0,
                                      const size_t hs1, const size_t hp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / alpha;
    const auto f_6 = beta / (alpha * p);

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

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_11 = buffer.data(gp + 11);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_35 = buffer.data(gd + 35);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_11 = buffer.data(hs1 + 11);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, gp_0, gd_0, hs0_0, hs1_0, \
                         hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];

        t_3[k] = pa_y[k] * gd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, fd0_0, fd0_3, fd1_0, fd1_3, gd_0, \
                         gd_3, gd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * gd_0[k];

        t_5[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_y[k] * gd_3[k];

        t_6[k] = f_5 * fd0_3[k]
                 - f_6 * fd1_3[k]
                 + pa_x[k] * gd_10[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, fd0_0, fd1_0, gd_6, hs0_1, hs0_2, \
                         hs1_1, hs1_2, hp_3, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_1[k]
                 + pb_z[k] * hp_3[k];

        t_8[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_z[k] * gd_6[k];

        t_9[k] = f_1 * hs0_2[k]
                 - f_2 * hs1_2[k]
                 + pb_y[k] * hp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, fd0_1, fd0_4, fd0_5, fd1_1, fd1_4, \
                         fd1_5, gd_9, gd_16, gd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fd0_4[k]
                  - f_6 * fd1_4[k]
                  + pa_x[k] * gd_16[k];

        t_11[k] = f_5 * fd0_1[k]
                  - f_6 * fd1_1[k]
                  + pa_y[k] * gd_9[k];

        t_12[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_5[k]
                  + pa_x[k] * gd_18[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, fd0_2, fd1_2, gd_14, hs0_3, \
                         hs0_4, hs1_3, hs1_4, hp_5, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];

        t_14[k] = f_5 * fd0_2[k]
                  - f_6 * fd1_2[k]
                  + pa_z[k] * gd_14[k];

        t_15[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_x, fd0_8, fd1_8, gd_21, gd_23, \
                         gd_35, hs0_7, hs1_7, hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pa_x[k] * gd_21[k];

        t_17[k] = pa_x[k] * gd_23[k];

        t_18[k] = pa_x[k] * gd_35[k];

        t_19[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_7[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, pb_z, gp_6, gd_23, hs0_7, hs1_7, hp_8, \
                         hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * gp_6[k]
                  + f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_y[k] * hp_8[k];

        t_21[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_z[k] * hp_9[k];

        t_22[k] = pa_z[k] * gd_23[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pa_z, pb_x, fd0_5, fd0_7, fd1_5, fd1_7, \
                         gd_26, gd_30, hs0_8, hs1_8, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_10[k];

        t_24[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_5[k]
                  + pa_z[k] * gd_26[k];

        t_25[k] = f_5 * fd0_7[k]
                  - f_6 * fd1_7[k]
                  + pa_y[k] * gd_30[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pa_z, pb_x, fd0_6, fd0_8, fd1_6, fd1_8, \
                         gd_29, gd_32, hs0_9, hs1_9, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * hs0_9[k]
                  - f_2 * hs1_9[k]
                  + pb_x[k] * hp_11[k];

        t_27[k] = f_5 * fd0_6[k]
                  - f_6 * fd1_6[k]
                  + pa_z[k] * gd_29[k];

        t_28[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pa_y[k] * gd_32[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, gp_11, gd_35, hs0_11, \
                         hs1_11, hp_12, hp_13, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * gd_35[k];

        t_30[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_x[k] * hp_12[k];

        t_31[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_y[k] * hp_13[k];

        t_32[k] = f_0 * gp_11[k]
                  + f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t fd0, const size_t fd1,
                                      const size_t gp, const size_t gd, const size_t hs0,
                                      const size_t hs1, const size_t hp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_12 = buffer.data(fd1 + 12);
    const auto *fd1_14 = buffer.data(fd1 + 14);
    const auto *fd1_17 = buffer.data(fd1 + 17);
    const auto *fd1_20 = buffer.data(fd1 + 20);

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
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_11 = buffer.data(hs1 + 11);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, gp_0, gd_0, hs0_0, hs1_0, \
                         hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];

        t_3[k] = pa_y[k] * gd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, fd0_0, fd1_0, gp_1, gp_2, \
                         gd_0, gd_1, gd_2, gd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * gp_1[k]
                 + pa_y[k] * gd_1[k];

        t_5[k] = pa_y[k] * gd_2[k];

        t_6[k] = pa_z[k] * gd_0[k];

        t_7[k] = pa_z[k] * gd_1[k];

        t_8[k] = f_3 * gp_2[k]
                 + pa_z[k] * gd_2[k];

        t_9[k] = f_4 * fd0_0[k]
                 - f_5 * fd1_0[k]
                 + pa_y[k] * gd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, pa_z, pb_z, fd0_3, fd1_8, gd_4, \
                         gd_6, gd_8, hs0_1, hs1_1, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_6 * fd0_3[k]
                  - f_7 * fd1_8[k]
                  + pa_x[k] * gd_8[k];

        t_11[k] = f_1 * hs0_1[k]
                  - f_2 * hs1_1[k]
                  + pb_z[k] * hp_3[k];

        t_12[k] = pa_z[k] * gd_4[k];

        t_13[k] = pa_y[k] * gd_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_y, fd0_0, fd0_4, fd1_0, fd1_10, \
                         gd_5, gd_12, hs0_2, hs1_2, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_4 * fd0_0[k]
                  - f_5 * fd1_0[k]
                  + pa_z[k] * gd_5[k];

        t_15[k] = f_1 * hs0_2[k]
                  - f_2 * hs1_2[k]
                  + pb_y[k] * hp_4[k];

        t_16[k] = f_6 * fd0_4[k]
                  - f_7 * fd1_10[k]
                  + pa_x[k] * gd_12[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, pb_z, fd0_1, fd0_5, fd1_3, fd1_12, \
                         gd_7, gd_14, hs0_3, hs1_3, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * fd0_1[k]
                  - f_7 * fd1_3[k]
                  + pa_y[k] * gd_7[k];

        t_18[k] = f_4 * fd0_5[k]
                  - f_5 * fd1_12[k]
                  + pa_x[k] * gd_14[k];

        t_19[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pa_z, gp_3, gp_4, gd_7, gd_8, \
                         gd_9, gd_11, gd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * gd_7[k];

        t_21[k] = pa_z[k] * gd_8[k];

        t_22[k] = f_3 * gp_3[k]
                  + pa_z[k] * gd_9[k];

        t_23[k] = f_3 * gp_4[k]
                  + pa_y[k] * gd_11[k];

        t_24[k] = pa_y[k] * gd_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_x, pa_z, pb_y, fd0_2, fd0_8, fd1_5, fd1_20, \
                         gd_10, gd_16, hs0_4, hs1_4, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_6 * fd0_2[k]
                  - f_7 * fd1_5[k]
                  + pa_z[k] * gd_10[k];

        t_26[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];

        t_27[k] = f_4 * fd0_8[k]
                  - f_5 * fd1_20[k]
                  + pa_x[k] * gd_16[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_x, pa_z, gp_5, gp_8, gp_9, gd_13, \
                         gd_17, gd_18, gd_22, gd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * gp_5[k]
                  + pa_x[k] * gd_17[k];

        t_29[k] = pa_x[k] * gd_18[k];

        t_30[k] = pa_z[k] * gd_13[k];

        t_31[k] = f_3 * gp_8[k]
                  + pa_x[k] * gd_22[k];

        t_32[k] = f_3 * gp_9[k]
                  + pa_x[k] * gd_27[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_x, pb_x, pb_y, pb_z, gp_6, gd_29, hs0_7, \
                         hs1_7, hp_7, hp_8, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_x[k] * gd_29[k];

        t_34[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_7[k];

        t_35[k] = f_0 * gp_6[k]
                  + f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_y[k] * hp_8[k];

        t_36[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_z[k] * hp_9[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, gp_7, gd_17, gd_18, gd_19, hs0_8, \
                         hs1_8, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pa_z[k] * gd_17[k];

        t_38[k] = pa_z[k] * gd_18[k];

        t_39[k] = f_3 * gp_7[k]
                  + pa_z[k] * gd_19[k];

        t_40[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_y, pa_z, pb_x, fd0_5, fd0_7, fd1_12, fd1_17, \
                         gd_20, gd_24, hs0_9, hs1_9, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_4 * fd0_5[k]
                  - f_5 * fd1_12[k]
                  + pa_z[k] * gd_20[k];

        t_42[k] = f_6 * fd0_7[k]
                  - f_7 * fd1_17[k]
                  + pa_y[k] * gd_24[k];

        t_43[k] = f_1 * hs0_9[k]
                  - f_2 * hs1_9[k]
                  + pb_x[k] * hp_11[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_y, pa_z, fd0_6, fd0_8, fd1_14, fd1_20, \
                         gp_10, gd_23, gd_26, gd_28, gd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_6 * fd0_6[k]
                  - f_7 * fd1_14[k]
                  + pa_z[k] * gd_23[k];

        t_45[k] = f_4 * fd0_8[k]
                  - f_5 * fd1_20[k]
                  + pa_y[k] * gd_26[k];

        t_46[k] = f_3 * gp_10[k]
                  + pa_y[k] * gd_28[k];

        t_47[k] = pa_y[k] * gd_29[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_x, pb_y, pb_z, gp_11, hs0_11, hs1_11, hp_12, \
                         hp_13, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_x[k] * hp_12[k];

        t_49[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_y[k] * hp_13[k];

        t_50[k] = f_0 * gp_11[k]
                  + f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t fd0, const size_t fd1,
                                      const size_t gp, const size_t gd, const size_t hs0,
                                      const size_t hs1, const size_t hp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);

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

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_3 = buffer.data(fd0 + 3);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_10 = buffer.data(fd0 + 10);
    const auto *fd0_12 = buffer.data(fd0 + 12);
    const auto *fd0_14 = buffer.data(fd0 + 14);
    const auto *fd0_17 = buffer.data(fd0 + 17);
    const auto *fd0_20 = buffer.data(fd0 + 20);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_4 = buffer.data(fd1 + 4);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_11 = buffer.data(fd1 + 11);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);

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

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_11 = buffer.data(hs1 + 11);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, gp_0, gd_0, hs0_0, hs1_0, \
                         hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];

        t_3[k] = pa_y[k] * gd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, fd0_0, fd1_0, gp_1, gp_2, gd_0, gd_1, \
                         gd_2, gd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * gp_1[k]
                 + pa_y[k] * gd_1[k];

        t_5[k] = pa_z[k] * gd_0[k];

        t_6[k] = f_3 * gp_2[k]
                 + pa_z[k] * gd_2[k];

        t_7[k] = f_4 * fd0_0[k]
                 - f_5 * fd1_0[k]
                 + pa_y[k] * gd_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_x, pa_z, pb_z, fd0_0, fd0_8, fd1_0, fd1_5, gd_4, \
                         gd_6, hs0_1, hs1_1, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_6 * fd0_8[k]
                 - f_7 * fd1_5[k]
                 + pa_x[k] * gd_6[k];

        t_9[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_1[k]
                 + pb_z[k] * hp_3[k];

        t_10[k] = f_4 * fd0_0[k]
                  - f_5 * fd1_0[k]
                  + pa_z[k] * gd_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, pb_y, fd0_3, fd0_10, fd1_3, fd1_6, \
                         gd_5, gd_8, hs0_2, hs1_2, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_1 * hs0_2[k]
                  - f_2 * hs1_2[k]
                  + pb_y[k] * hp_4[k];

        t_12[k] = f_6 * fd0_10[k]
                  - f_7 * fd1_6[k]
                  + pa_x[k] * gd_8[k];

        t_13[k] = f_6 * fd0_3[k]
                  - f_7 * fd1_3[k]
                  + pa_y[k] * gd_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_y, pb_z, fd0_12, fd1_8, gd_7, gd_9, hs0_3, \
                         hs1_3, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_4 * fd0_12[k]
                  - f_5 * fd1_8[k]
                  + pa_x[k] * gd_9[k];

        t_15[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];

        t_16[k] = pa_y[k] * gd_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_z, pb_y, fd0_5, fd0_20, fd1_4, fd1_14, \
                         gd_7, gd_10, hs0_4, hs1_4, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * fd0_5[k]
                  - f_7 * fd1_4[k]
                  + pa_z[k] * gd_7[k];

        t_18[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];

        t_19[k] = f_4 * fd0_20[k]
                  - f_5 * fd1_14[k]
                  + pa_x[k] * gd_10[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, pa_x, gp_5, gp_9, gd_11, gd_12, \
                         gd_15, gd_16, gd_18, gd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_3 * gp_5[k]
                  + pa_x[k] * gd_11[k];

        t_21[k] = pa_x[k] * gd_12[k];

        t_22[k] = pa_x[k] * gd_15[k];

        t_23[k] = pa_x[k] * gd_16[k];

        t_24[k] = f_3 * gp_9[k]
                  + pa_x[k] * gd_18[k];

        t_25[k] = pa_x[k] * gd_20[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, pb_z, gp_6, gd_12, hs0_7, \
                         hs1_7, hp_7, hp_8, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_7[k];

        t_27[k] = f_0 * gp_6[k]
                  + f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_y[k] * hp_8[k];

        t_28[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_z[k] * hp_9[k];

        t_29[k] = pa_z[k] * gd_12[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_z, pb_x, fd0_12, fd1_8, gp_7, gd_13, gd_14, \
                         hs0_8, hs1_8, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * gp_7[k]
                  + pa_z[k] * gd_13[k];

        t_31[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_10[k];

        t_32[k] = f_4 * fd0_12[k]
                  - f_5 * fd1_8[k]
                  + pa_z[k] * gd_14[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pa_z, pb_x, fd0_14, fd0_17, fd1_10, fd1_11, \
                         gd_15, gd_16, hs0_9, hs1_9, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_6 * fd0_17[k]
                  - f_7 * fd1_11[k]
                  + pa_y[k] * gd_16[k];

        t_34[k] = f_1 * hs0_9[k]
                  - f_2 * hs1_9[k]
                  + pb_x[k] * hp_11[k];

        t_35[k] = f_6 * fd0_14[k]
                  - f_7 * fd1_10[k]
                  + pa_z[k] * gd_15[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_y, pb_x, fd0_20, fd1_14, gp_10, gd_17, \
                         gd_19, gd_20, hs0_11, hs1_11, hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_4 * fd0_20[k]
                  - f_5 * fd1_14[k]
                  + pa_y[k] * gd_17[k];

        t_37[k] = f_3 * gp_10[k]
                  + pa_y[k] * gd_19[k];

        t_38[k] = pa_y[k] * gd_20[k];

        t_39[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_x[k] * hp_12[k];
    }

#pragma omp simd aligned(t_40, t_41, pb_y, pb_z, gp_11, hs0_11, hs1_11, hp_13, \
                         hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_y[k] * hp_13[k];

        t_41[k] = f_0 * gp_11[k]
                  + f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t fd0, const size_t fd1,
                                      const size_t gp, const size_t gd, const size_t hs0,
                                      const size_t hs1, const size_t hp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / alpha;
    const auto f_6 = beta / (alpha * p);

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

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_11 = buffer.data(gp + 11);

    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_11 = buffer.data(hs1 + 11);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs0_0, hs1_0, hp_0, hp_1, \
                         hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, fd0_0, fd0_3, fd1_0, fd1_3, gd_3, \
                         gd_6, hs0_1, hs1_1, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_y[k] * gd_3[k];

        t_4[k] = f_5 * fd0_3[k]
                 - f_6 * fd1_3[k]
                 + pa_x[k] * gd_6[k];

        t_5[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_1[k]
                 + pb_z[k] * hp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, fd0_0, fd0_4, fd1_0, fd1_4, gd_4, \
                         gd_10, hs0_2, hs1_2, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_z[k] * gd_4[k];

        t_7[k] = f_1 * hs0_2[k]
                 - f_2 * hs1_2[k]
                 + pb_y[k] * hp_4[k];

        t_8[k] = f_5 * fd0_4[k]
                 - f_6 * fd1_4[k]
                 + pa_x[k] * gd_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, fd0_1, fd0_5, fd1_1, fd1_5, gd_5, \
                         gd_11, hs0_3, hs1_3, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * fd0_1[k]
                 - f_6 * fd1_1[k]
                 + pa_y[k] * gd_5[k];

        t_10[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_5[k]
                  + pa_x[k] * gd_11[k];

        t_11[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, fd0_2, fd0_8, fd1_2, fd1_8, gd_8, \
                         gd_12, hs0_4, hs1_4, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * fd0_2[k]
                  - f_6 * fd1_2[k]
                  + pa_z[k] * gd_8[k];

        t_13[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];

        t_14[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pa_x[k] * gd_12[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pb_x, pb_y, gp_6, gd_14, gd_23, hs0_5, \
                         hs1_7, hp_7, hp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_x[k] * gd_14[k];

        t_16[k] = pa_x[k] * gd_23[k];

        t_17[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_7[k];

        t_18[k] = f_0 * gp_6[k]
                  + f_1 * hs0_5[k]
                  - f_2 * hs1_7[k]
                  + pb_y[k] * hp_8[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_x, pb_z, fd0_5, fd1_5, gd_16, hs0_5, \
                         hs0_6, hs1_7, hs1_8, hp_9, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * hs0_5[k]
                  - f_2 * hs1_7[k]
                  + pb_z[k] * hp_9[k];

        t_20[k] = f_1 * hs0_6[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_10[k];

        t_21[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_5[k]
                  + pa_z[k] * gd_16[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pa_z, pb_x, fd0_6, fd0_7, fd1_6, fd1_7, \
                         gd_18, gd_19, hs0_7, hs1_9, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * fd0_7[k]
                  - f_6 * fd1_7[k]
                  + pa_y[k] * gd_19[k];

        t_23[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_9[k]
                  + pb_x[k] * hp_11[k];

        t_24[k] = f_5 * fd0_6[k]
                  - f_6 * fd1_6[k]
                  + pa_z[k] * gd_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, pb_x, pb_y, fd0_8, fd1_8, gd_20, gd_23, \
                         hs0_8, hs1_11, hp_12, hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pa_y[k] * gd_20[k];

        t_26[k] = pa_y[k] * gd_23[k];

        t_27[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_11[k]
                  + pb_x[k] * hp_12[k];

        t_28[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_11[k]
                  + pb_y[k] * hp_13[k];
    }

#pragma omp simd aligned(t_29, pb_z, gp_11, hs0_8, hs1_11, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * gp_11[k]
                  + f_1 * hs0_8[k]
                  - f_2 * hs1_11[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t fd0, const size_t fd1,
                                      const size_t gp, const size_t gd, const size_t hs0,
                                      const size_t hs1, const size_t hp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / alpha;
    const auto f_6 = beta / (alpha * p);

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
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_4 = buffer.data(fd1 + 4);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_11 = buffer.data(fd1 + 11);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_11 = buffer.data(gp + 11);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_26 = buffer.data(gd + 26);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_11 = buffer.data(hs1 + 11);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, gp_0, gd_0, hs0_0, hs1_0, \
                         hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];

        t_3[k] = pa_y[k] * gd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, fd0_0, fd0_3, fd1_0, fd1_5, gd_0, \
                         gd_3, gd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * gd_0[k];

        t_5[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_y[k] * gd_3[k];

        t_6[k] = f_5 * fd0_3[k]
                 - f_6 * fd1_5[k]
                 + pa_x[k] * gd_7[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, fd0_0, fd1_0, gd_4, hs0_1, hs0_2, \
                         hs1_1, hs1_2, hp_3, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_1[k]
                 + pb_z[k] * hp_3[k];

        t_8[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_z[k] * gd_4[k];

        t_9[k] = f_1 * hs0_2[k]
                 - f_2 * hs1_2[k]
                 + pb_y[k] * hp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, fd0_1, fd0_4, fd0_5, fd1_3, fd1_6, \
                         fd1_8, gd_6, gd_11, gd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fd0_4[k]
                  - f_6 * fd1_6[k]
                  + pa_x[k] * gd_11[k];

        t_11[k] = f_5 * fd0_1[k]
                  - f_6 * fd1_3[k]
                  + pa_y[k] * gd_6[k];

        t_12[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_8[k]
                  + pa_x[k] * gd_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, fd0_2, fd1_4, gd_9, hs0_3, hs0_4, \
                         hs1_3, hs1_4, hp_5, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];

        t_14[k] = f_5 * fd0_2[k]
                  - f_6 * fd1_4[k]
                  + pa_z[k] * gd_9[k];

        t_15[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_x, fd0_8, fd1_14, gd_13, gd_15, \
                         gd_26, hs0_7, hs1_7, hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_14[k]
                  + pa_x[k] * gd_13[k];

        t_17[k] = pa_x[k] * gd_15[k];

        t_18[k] = pa_x[k] * gd_26[k];

        t_19[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_7[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, pb_z, gp_6, gd_15, hs0_7, hs1_7, hp_8, \
                         hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * gp_6[k]
                  + f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_y[k] * hp_8[k];

        t_21[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_z[k] * hp_9[k];

        t_22[k] = pa_z[k] * gd_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pa_z, pb_x, fd0_5, fd0_7, fd1_8, fd1_11, \
                         gd_17, gd_21, hs0_8, hs1_8, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_10[k];

        t_24[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_8[k]
                  + pa_z[k] * gd_17[k];

        t_25[k] = f_5 * fd0_7[k]
                  - f_6 * fd1_11[k]
                  + pa_y[k] * gd_21[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pa_z, pb_x, fd0_6, fd0_8, fd1_10, fd1_14, \
                         gd_20, gd_23, hs0_9, hs1_9, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * hs0_9[k]
                  - f_2 * hs1_9[k]
                  + pb_x[k] * hp_11[k];

        t_27[k] = f_5 * fd0_6[k]
                  - f_6 * fd1_10[k]
                  + pa_z[k] * gd_20[k];

        t_28[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_14[k]
                  + pa_y[k] * gd_23[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, gp_11, gd_26, hs0_11, \
                         hs1_11, hp_12, hp_13, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * gd_26[k];

        t_30[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_x[k] * hp_12[k];

        t_31[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_y[k] * hp_13[k];

        t_32[k] = f_0 * gp_11[k]
                  + f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t fd0, const size_t fd1,
                                      const size_t gp, const size_t gd, const size_t hs0,
                                      const size_t hs1, const size_t hp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_3 = buffer.data(fd0 + 3);
    const auto *fd0_4 = buffer.data(fd0 + 4);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_6 = buffer.data(fd0 + 6);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_10 = buffer.data(fd0 + 10);
    const auto *fd0_11 = buffer.data(fd0 + 11);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_4 = buffer.data(fd1 + 4);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_11 = buffer.data(fd1 + 11);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_11 = buffer.data(hs1 + 11);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, gp_0, gd_0, hs0_0, hs1_0, \
                         hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];

        t_3[k] = pa_y[k] * gd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, fd0_0, fd0_5, fd1_0, fd1_5, \
                         gp_2, gd_0, gd_2, gd_3, gd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * gd_0[k];

        t_5[k] = f_3 * gp_2[k]
                 + pa_z[k] * gd_2[k];

        t_6[k] = f_4 * fd0_0[k]
                 - f_5 * fd1_0[k]
                 + pa_y[k] * gd_3[k];

        t_7[k] = f_6 * fd0_5[k]
                 - f_7 * fd1_5[k]
                 + pa_x[k] * gd_6[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_z, pb_y, pb_z, fd0_0, fd1_0, gd_4, hs0_1, hs0_2, \
                         hs1_1, hs1_2, hp_3, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_1[k]
                 + pb_z[k] * hp_3[k];

        t_9[k] = f_4 * fd0_0[k]
                 - f_5 * fd1_0[k]
                 + pa_z[k] * gd_4[k];

        t_10[k] = f_1 * hs0_2[k]
                  - f_2 * hs1_2[k]
                  + pb_y[k] * hp_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, fd0_3, fd0_6, fd0_8, fd1_3, fd1_6, \
                         fd1_8, gd_5, gd_10, gd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_6 * fd0_6[k]
                  - f_7 * fd1_6[k]
                  + pa_x[k] * gd_10[k];

        t_12[k] = f_6 * fd0_3[k]
                  - f_7 * fd1_3[k]
                  + pa_y[k] * gd_5[k];

        t_13[k] = f_4 * fd0_8[k]
                  - f_5 * fd1_8[k]
                  + pa_x[k] * gd_11[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_z, pb_y, pb_z, fd0_4, fd1_4, gd_8, hs0_3, hs0_4, \
                         hs1_3, hs1_4, hp_5, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];

        t_15[k] = f_6 * fd0_4[k]
                  - f_7 * fd1_4[k]
                  + pa_z[k] * gd_8[k];

        t_16[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_x, pb_x, fd0_14, fd1_14, gd_12, gd_14, \
                         gd_23, hs0_7, hs1_7, hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_4 * fd0_14[k]
                  - f_5 * fd1_14[k]
                  + pa_x[k] * gd_12[k];

        t_18[k] = pa_x[k] * gd_14[k];

        t_19[k] = pa_x[k] * gd_23[k];

        t_20[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_7[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, gp_6, gp_7, gd_14, gd_15, \
                         hs0_7, hs1_7, hp_8, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * gp_6[k]
                  + f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_y[k] * hp_8[k];

        t_22[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_z[k] * hp_9[k];

        t_23[k] = pa_z[k] * gd_14[k];

        t_24[k] = f_3 * gp_7[k]
                  + pa_z[k] * gd_15[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pb_x, fd0_8, fd0_11, fd1_8, fd1_11, \
                         gd_16, gd_19, hs0_8, hs1_8, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_10[k];

        t_26[k] = f_4 * fd0_8[k]
                  - f_5 * fd1_8[k]
                  + pa_z[k] * gd_16[k];

        t_27[k] = f_6 * fd0_11[k]
                  - f_7 * fd1_11[k]
                  + pa_y[k] * gd_19[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pa_z, pb_x, fd0_10, fd0_14, fd1_10, fd1_14, \
                         gd_18, gd_20, hs0_9, hs1_9, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * hs0_9[k]
                  - f_2 * hs1_9[k]
                  + pb_x[k] * hp_11[k];

        t_29[k] = f_6 * fd0_10[k]
                  - f_7 * fd1_10[k]
                  + pa_z[k] * gd_18[k];

        t_30[k] = f_4 * fd0_14[k]
                  - f_5 * fd1_14[k]
                  + pa_y[k] * gd_20[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pb_x, pb_y, gp_10, gd_22, gd_23, \
                         hs0_11, hs1_11, hp_12, hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_3 * gp_10[k]
                  + pa_y[k] * gd_22[k];

        t_32[k] = pa_y[k] * gd_23[k];

        t_33[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_x[k] * hp_12[k];

        t_34[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_y[k] * hp_13[k];
    }

#pragma omp simd aligned(t_35, pb_z, gp_11, hs0_11, hs1_11, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * gp_11[k]
                  + f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_18(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t fd0, const size_t fd1,
                                      const size_t gp, const size_t gd, const size_t hs0,
                                      const size_t hs1, const size_t hp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / alpha;
    const auto f_6 = beta / (alpha * p);

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

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_3 = buffer.data(fd0 + 3);
    const auto *fd0_4 = buffer.data(fd0 + 4);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_6 = buffer.data(fd0 + 6);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_10 = buffer.data(fd0 + 10);
    const auto *fd0_11 = buffer.data(fd0 + 11);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_4 = buffer.data(fd1 + 4);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_11 = buffer.data(fd1 + 11);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_11 = buffer.data(gp + 11);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_20 = buffer.data(gd + 20);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_11 = buffer.data(hs1 + 11);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, gp_0, gd_0, hs0_0, hs1_0, \
                         hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];

        t_3[k] = pa_y[k] * gd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, fd0_0, fd0_5, fd1_0, fd1_5, gd_0, \
                         gd_3, gd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * gd_0[k];

        t_5[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_y[k] * gd_3[k];

        t_6[k] = f_5 * fd0_5[k]
                 - f_6 * fd1_5[k]
                 + pa_x[k] * gd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, fd0_0, fd1_0, gd_4, hs0_1, hs0_2, \
                         hs1_1, hs1_2, hp_3, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_1[k]
                 + pb_z[k] * hp_3[k];

        t_8[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_z[k] * gd_4[k];

        t_9[k] = f_1 * hs0_2[k]
                 - f_2 * hs1_2[k]
                 + pb_y[k] * hp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, fd0_3, fd0_6, fd0_8, fd1_3, fd1_6, \
                         fd1_8, gd_5, gd_8, gd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fd0_6[k]
                  - f_6 * fd1_6[k]
                  + pa_x[k] * gd_8[k];

        t_11[k] = f_5 * fd0_3[k]
                  - f_6 * fd1_3[k]
                  + pa_y[k] * gd_5[k];

        t_12[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pa_x[k] * gd_9[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, fd0_4, fd1_4, gd_7, hs0_3, hs0_4, \
                         hs1_3, hs1_4, hp_5, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];

        t_14[k] = f_5 * fd0_4[k]
                  - f_6 * fd1_4[k]
                  + pa_z[k] * gd_7[k];

        t_15[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_x, fd0_14, fd1_14, gd_10, gd_12, \
                         gd_20, hs0_7, hs1_7, hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pa_x[k] * gd_10[k];

        t_17[k] = pa_x[k] * gd_12[k];

        t_18[k] = pa_x[k] * gd_20[k];

        t_19[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_7[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, pb_z, gp_6, gd_12, hs0_7, hs1_7, hp_8, \
                         hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * gp_6[k]
                  + f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_y[k] * hp_8[k];

        t_21[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_z[k] * hp_9[k];

        t_22[k] = pa_z[k] * gd_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pa_z, pb_x, fd0_8, fd0_11, fd1_8, fd1_11, \
                         gd_14, gd_16, hs0_8, hs1_8, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_10[k];

        t_24[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pa_z[k] * gd_14[k];

        t_25[k] = f_5 * fd0_11[k]
                  - f_6 * fd1_11[k]
                  + pa_y[k] * gd_16[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pa_z, pb_x, fd0_10, fd0_14, fd1_10, fd1_14, \
                         gd_15, gd_17, hs0_9, hs1_9, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * hs0_9[k]
                  - f_2 * hs1_9[k]
                  + pb_x[k] * hp_11[k];

        t_27[k] = f_5 * fd0_10[k]
                  - f_6 * fd1_10[k]
                  + pa_z[k] * gd_15[k];

        t_28[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pa_y[k] * gd_17[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, gp_11, gd_20, hs0_11, \
                         hs1_11, hp_12, hp_13, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * gd_20[k];

        t_30[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_x[k] * hp_12[k];

        t_31[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_y[k] * hp_13[k];

        t_32[k] = f_0 * gp_11[k]
                  + f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_19(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t fd0, const size_t fd1,
                                      const size_t gp, const size_t gd, const size_t hs0,
                                      const size_t hs1, const size_t hp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / alpha;
    const auto f_6 = beta / (alpha * p);

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

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_11 = buffer.data(gp + 11);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_11 = buffer.data(hs1 + 11);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, gp_0, gd_0, hs0_0, hs1_0, \
                         hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];

        t_3[k] = pa_y[k] * gd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, fd0_0, fd0_3, fd1_0, fd1_3, gd_0, \
                         gd_3, gd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * gd_0[k];

        t_5[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_y[k] * gd_3[k];

        t_6[k] = f_5 * fd0_3[k]
                 - f_6 * fd1_3[k]
                 + pa_x[k] * gd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, fd0_0, fd1_0, gd_4, hs0_1, hs0_2, \
                         hs1_1, hs1_2, hp_3, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_1[k]
                 + pb_z[k] * hp_3[k];

        t_8[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_z[k] * gd_4[k];

        t_9[k] = f_1 * hs0_2[k]
                 - f_2 * hs1_2[k]
                 + pb_y[k] * hp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, fd0_1, fd0_4, fd0_5, fd1_1, fd1_4, \
                         fd1_5, gd_5, gd_10, gd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fd0_4[k]
                  - f_6 * fd1_4[k]
                  + pa_x[k] * gd_10[k];

        t_11[k] = f_5 * fd0_1[k]
                  - f_6 * fd1_1[k]
                  + pa_y[k] * gd_5[k];

        t_12[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_5[k]
                  + pa_x[k] * gd_11[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, fd0_2, fd1_2, gd_8, hs0_3, hs0_4, \
                         hs1_3, hs1_4, hp_5, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];

        t_14[k] = f_5 * fd0_2[k]
                  - f_6 * fd1_2[k]
                  + pa_z[k] * gd_8[k];

        t_15[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_x, fd0_8, fd1_8, gd_12, gd_14, \
                         gd_23, hs0_7, hs1_7, hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pa_x[k] * gd_12[k];

        t_17[k] = pa_x[k] * gd_14[k];

        t_18[k] = pa_x[k] * gd_23[k];

        t_19[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_7[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, pb_z, gp_6, gd_14, hs0_7, hs1_7, hp_8, \
                         hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * gp_6[k]
                  + f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_y[k] * hp_8[k];

        t_21[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_z[k] * hp_9[k];

        t_22[k] = pa_z[k] * gd_14[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pa_z, pb_x, fd0_5, fd0_7, fd1_5, fd1_7, \
                         gd_16, gd_19, hs0_8, hs1_8, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_10[k];

        t_24[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_5[k]
                  + pa_z[k] * gd_16[k];

        t_25[k] = f_5 * fd0_7[k]
                  - f_6 * fd1_7[k]
                  + pa_y[k] * gd_19[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pa_z, pb_x, fd0_6, fd0_8, fd1_6, fd1_8, \
                         gd_18, gd_20, hs0_9, hs1_9, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * hs0_9[k]
                  - f_2 * hs1_9[k]
                  + pb_x[k] * hp_11[k];

        t_27[k] = f_5 * fd0_6[k]
                  - f_6 * fd1_6[k]
                  + pa_z[k] * gd_18[k];

        t_28[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pa_y[k] * gd_20[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, gp_11, gd_23, hs0_11, \
                         hs1_11, hp_12, hp_13, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * gd_23[k];

        t_30[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_x[k] * hp_12[k];

        t_31[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_y[k] * hp_13[k];

        t_32[k] = f_0 * gp_11[k]
                  + f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_20(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t fd0, const size_t fd1,
                                      const size_t gp, const size_t gd, const size_t hs0,
                                      const size_t hs1, const size_t hp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / alpha;
    const auto f_6 = beta / (alpha * p);

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
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_4 = buffer.data(fd1 + 4);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_11 = buffer.data(fd1 + 11);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_11 = buffer.data(gp + 11);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_11 = buffer.data(hs1 + 11);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, gp_0, gd_0, hs0_0, hs1_0, \
                         hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];

        t_3[k] = pa_y[k] * gd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, fd0_0, fd0_3, fd1_0, fd1_5, gd_0, \
                         gd_3, gd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * gd_0[k];

        t_5[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_y[k] * gd_3[k];

        t_6[k] = f_5 * fd0_3[k]
                 - f_6 * fd1_5[k]
                 + pa_x[k] * gd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, fd0_0, fd1_0, gd_4, hs0_1, hs0_2, \
                         hs1_1, hs1_2, hp_3, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_1[k]
                 + pb_z[k] * hp_3[k];

        t_8[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_z[k] * gd_4[k];

        t_9[k] = f_1 * hs0_2[k]
                 - f_2 * hs1_2[k]
                 + pb_y[k] * hp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, fd0_1, fd0_4, fd0_5, fd1_3, fd1_6, \
                         fd1_8, gd_5, gd_10, gd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fd0_4[k]
                  - f_6 * fd1_6[k]
                  + pa_x[k] * gd_10[k];

        t_11[k] = f_5 * fd0_1[k]
                  - f_6 * fd1_3[k]
                  + pa_y[k] * gd_5[k];

        t_12[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_8[k]
                  + pa_x[k] * gd_11[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, fd0_2, fd1_4, gd_8, hs0_3, hs0_4, \
                         hs1_3, hs1_4, hp_5, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];

        t_14[k] = f_5 * fd0_2[k]
                  - f_6 * fd1_4[k]
                  + pa_z[k] * gd_8[k];

        t_15[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_x, fd0_8, fd1_14, gd_12, gd_14, \
                         gd_23, hs0_7, hs1_7, hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_14[k]
                  + pa_x[k] * gd_12[k];

        t_17[k] = pa_x[k] * gd_14[k];

        t_18[k] = pa_x[k] * gd_23[k];

        t_19[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_7[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, pb_z, gp_6, gd_14, hs0_7, hs1_7, hp_8, \
                         hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * gp_6[k]
                  + f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_y[k] * hp_8[k];

        t_21[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_z[k] * hp_9[k];

        t_22[k] = pa_z[k] * gd_14[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pa_z, pb_x, fd0_5, fd0_7, fd1_8, fd1_11, \
                         gd_16, gd_19, hs0_8, hs1_8, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_10[k];

        t_24[k] = f_3 * fd0_5[k]
                  - f_4 * fd1_8[k]
                  + pa_z[k] * gd_16[k];

        t_25[k] = f_5 * fd0_7[k]
                  - f_6 * fd1_11[k]
                  + pa_y[k] * gd_19[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pa_z, pb_x, fd0_6, fd0_8, fd1_10, fd1_14, \
                         gd_18, gd_20, hs0_9, hs1_9, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * hs0_9[k]
                  - f_2 * hs1_9[k]
                  + pb_x[k] * hp_11[k];

        t_27[k] = f_5 * fd0_6[k]
                  - f_6 * fd1_10[k]
                  + pa_z[k] * gd_18[k];

        t_28[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_14[k]
                  + pa_y[k] * gd_20[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, gp_11, gd_23, hs0_11, \
                         hs1_11, hp_12, hp_13, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * gd_23[k];

        t_30[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_x[k] * hp_12[k];

        t_31[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_y[k] * hp_13[k];

        t_32[k] = f_0 * gp_11[k]
                  + f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_z[k] * hp_14[k];
    }
}

auto
compute_prim_hd_electron_repulsion_21(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t fd0, const size_t fd1,
                                      const size_t gp, const size_t gd, const size_t hs0,
                                      const size_t hs1, const size_t hp, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.0 / alpha;
    const auto f_6 = beta / (alpha * p);

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

    const auto *fd0_0 = buffer.data(fd0 + 0);
    const auto *fd0_3 = buffer.data(fd0 + 3);
    const auto *fd0_4 = buffer.data(fd0 + 4);
    const auto *fd0_5 = buffer.data(fd0 + 5);
    const auto *fd0_6 = buffer.data(fd0 + 6);
    const auto *fd0_8 = buffer.data(fd0 + 8);
    const auto *fd0_10 = buffer.data(fd0 + 10);
    const auto *fd0_11 = buffer.data(fd0 + 11);
    const auto *fd0_14 = buffer.data(fd0 + 14);

    const auto *fd1_0 = buffer.data(fd1 + 0);
    const auto *fd1_3 = buffer.data(fd1 + 3);
    const auto *fd1_4 = buffer.data(fd1 + 4);
    const auto *fd1_5 = buffer.data(fd1 + 5);
    const auto *fd1_6 = buffer.data(fd1 + 6);
    const auto *fd1_8 = buffer.data(fd1 + 8);
    const auto *fd1_10 = buffer.data(fd1 + 10);
    const auto *fd1_11 = buffer.data(fd1 + 11);
    const auto *fd1_14 = buffer.data(fd1 + 14);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_11 = buffer.data(gp + 11);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_11 = buffer.data(hs1 + 11);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, gp_0, gd_0, hs0_0, hs1_0, \
                         hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 + f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = f_1 * hs0_0[k]
                 - f_2 * hs1_0[k]
                 + pb_z[k] * hp_2[k];

        t_3[k] = pa_y[k] * gd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, fd0_0, fd0_5, fd1_0, fd1_5, gd_0, \
                         gd_3, gd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * gd_0[k];

        t_5[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_y[k] * gd_3[k];

        t_6[k] = f_5 * fd0_5[k]
                 - f_6 * fd1_5[k]
                 + pa_x[k] * gd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, fd0_0, fd1_0, gd_4, hs0_1, hs0_2, \
                         hs1_1, hs1_2, hp_3, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * hs0_1[k]
                 - f_2 * hs1_1[k]
                 + pb_z[k] * hp_3[k];

        t_8[k] = f_3 * fd0_0[k]
                 - f_4 * fd1_0[k]
                 + pa_z[k] * gd_4[k];

        t_9[k] = f_1 * hs0_2[k]
                 - f_2 * hs1_2[k]
                 + pb_y[k] * hp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, fd0_3, fd0_6, fd0_8, fd1_3, fd1_6, \
                         fd1_8, gd_5, gd_10, gd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * fd0_6[k]
                  - f_6 * fd1_6[k]
                  + pa_x[k] * gd_10[k];

        t_11[k] = f_5 * fd0_3[k]
                  - f_6 * fd1_3[k]
                  + pa_y[k] * gd_5[k];

        t_12[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pa_x[k] * gd_11[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, fd0_4, fd1_4, gd_8, hs0_3, hs0_4, \
                         hs1_3, hs1_4, hp_5, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * hs0_3[k]
                  - f_2 * hs1_3[k]
                  + pb_z[k] * hp_5[k];

        t_14[k] = f_5 * fd0_4[k]
                  - f_6 * fd1_4[k]
                  + pa_z[k] * gd_8[k];

        t_15[k] = f_1 * hs0_4[k]
                  - f_2 * hs1_4[k]
                  + pb_y[k] * hp_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_x, fd0_14, fd1_14, gd_12, gd_14, \
                         gd_23, hs0_7, hs1_7, hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pa_x[k] * gd_12[k];

        t_17[k] = pa_x[k] * gd_14[k];

        t_18[k] = pa_x[k] * gd_23[k];

        t_19[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_x[k] * hp_7[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, pb_z, gp_6, gd_14, hs0_7, hs1_7, hp_8, \
                         hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * gp_6[k]
                  + f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_y[k] * hp_8[k];

        t_21[k] = f_1 * hs0_7[k]
                  - f_2 * hs1_7[k]
                  + pb_z[k] * hp_9[k];

        t_22[k] = pa_z[k] * gd_14[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pa_z, pb_x, fd0_8, fd0_11, fd1_8, fd1_11, \
                         gd_16, gd_19, hs0_8, hs1_8, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * hs0_8[k]
                  - f_2 * hs1_8[k]
                  + pb_x[k] * hp_10[k];

        t_24[k] = f_3 * fd0_8[k]
                  - f_4 * fd1_8[k]
                  + pa_z[k] * gd_16[k];

        t_25[k] = f_5 * fd0_11[k]
                  - f_6 * fd1_11[k]
                  + pa_y[k] * gd_19[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pa_z, pb_x, fd0_10, fd0_14, fd1_10, fd1_14, \
                         gd_18, gd_20, hs0_9, hs1_9, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * hs0_9[k]
                  - f_2 * hs1_9[k]
                  + pb_x[k] * hp_11[k];

        t_27[k] = f_5 * fd0_10[k]
                  - f_6 * fd1_10[k]
                  + pa_z[k] * gd_18[k];

        t_28[k] = f_3 * fd0_14[k]
                  - f_4 * fd1_14[k]
                  + pa_y[k] * gd_20[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, gp_11, gd_23, hs0_11, \
                         hs1_11, hp_12, hp_13, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * gd_23[k];

        t_30[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_x[k] * hp_12[k];

        t_31[k] = f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_y[k] * hp_13[k];

        t_32[k] = f_0 * gp_11[k]
                  + f_1 * hs0_11[k]
                  - f_2 * hs1_11[k]
                  + pb_z[k] * hp_14[k];
    }
}

}  // namespace simdt2ceri
