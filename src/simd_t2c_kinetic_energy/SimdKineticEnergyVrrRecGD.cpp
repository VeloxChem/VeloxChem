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


#include "SimdKineticEnergyVrrRecGD.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_gd_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dd_s, const size_t dd,
                                 const size_t fp, const size_t fd, const size_t gs_s,
                                 const size_t gd_s, const size_t gs, const size_t gp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 1.5 / p;
    const auto f_5 = 2.0 * beta / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_5 = buffer.data(dd_s + 5);
    const auto *dd_s_8 = buffer.data(dd_s + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);
    const auto *fp_15 = buffer.data(fp + 15);

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
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_22 = buffer.data(fd + 22);
    const auto *fd_23 = buffer.data(fd + 23);
    const auto *fd_24 = buffer.data(fd + 24);
    const auto *fd_25 = buffer.data(fd + 25);
    const auto *fd_26 = buffer.data(fd + 26);

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_3 = buffer.data(gs_s + 3);
    const auto *gs_s_4 = buffer.data(gs_s + 4);
    const auto *gs_s_7 = buffer.data(gs_s + 7);
    const auto *gs_s_9 = buffer.data(gs_s + 9);
    const auto *gs_s_11 = buffer.data(gs_s + 11);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_15 = buffer.data(gd_s + 15);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_22 = buffer.data(gd_s + 22);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_24 = buffer.data(gd_s + 24);
    const auto *gd_s_25 = buffer.data(gd_s + 25);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_27 = buffer.data(gd_s + 27);
    const auto *gd_s_28 = buffer.data(gd_s + 28);
    const auto *gd_s_29 = buffer.data(gd_s + 29);
    const auto *gd_s_30 = buffer.data(gd_s + 30);
    const auto *gd_s_31 = buffer.data(gd_s + 31);
    const auto *gd_s_32 = buffer.data(gd_s + 32);
    const auto *gd_s_33 = buffer.data(gd_s + 33);
    const auto *gd_s_34 = buffer.data(gd_s + 34);
    const auto *gd_s_35 = buffer.data(gd_s + 35);
    const auto *gd_s_36 = buffer.data(gd_s + 36);
    const auto *gd_s_37 = buffer.data(gd_s + 37);
    const auto *gd_s_38 = buffer.data(gd_s + 38);
    const auto *gd_s_39 = buffer.data(gd_s + 39);
    const auto *gd_s_40 = buffer.data(gd_s + 40);
    const auto *gd_s_41 = buffer.data(gd_s + 41);
    const auto *gd_s_42 = buffer.data(gd_s + 42);
    const auto *gd_s_43 = buffer.data(gd_s + 43);
    const auto *gd_s_44 = buffer.data(gd_s + 44);
    const auto *gd_s_45 = buffer.data(gd_s + 45);
    const auto *gd_s_46 = buffer.data(gd_s + 46);
    const auto *gd_s_47 = buffer.data(gd_s + 47);
    const auto *gd_s_48 = buffer.data(gd_s + 48);
    const auto *gd_s_49 = buffer.data(gd_s + 49);
    const auto *gd_s_50 = buffer.data(gd_s + 50);
    const auto *gd_s_51 = buffer.data(gd_s + 51);
    const auto *gd_s_52 = buffer.data(gd_s + 52);
    const auto *gd_s_53 = buffer.data(gd_s + 53);
    const auto *gd_s_54 = buffer.data(gd_s + 54);
    const auto *gd_s_55 = buffer.data(gd_s + 55);
    const auto *gd_s_56 = buffer.data(gd_s + 56);
    const auto *gd_s_57 = buffer.data(gd_s + 57);
    const auto *gd_s_58 = buffer.data(gd_s + 58);
    const auto *gd_s_59 = buffer.data(gd_s + 59);
    const auto *gd_s_60 = buffer.data(gd_s + 60);
    const auto *gd_s_61 = buffer.data(gd_s + 61);
    const auto *gd_s_62 = buffer.data(gd_s + 62);
    const auto *gd_s_63 = buffer.data(gd_s + 63);
    const auto *gd_s_64 = buffer.data(gd_s + 64);
    const auto *gd_s_65 = buffer.data(gd_s + 65);
    const auto *gd_s_66 = buffer.data(gd_s + 66);
    const auto *gd_s_67 = buffer.data(gd_s + 67);
    const auto *gd_s_68 = buffer.data(gd_s + 68);
    const auto *gd_s_69 = buffer.data(gd_s + 69);
    const auto *gd_s_70 = buffer.data(gd_s + 70);
    const auto *gd_s_71 = buffer.data(gd_s + 71);
    const auto *gd_s_72 = buffer.data(gd_s + 72);
    const auto *gd_s_73 = buffer.data(gd_s + 73);
    const auto *gd_s_74 = buffer.data(gd_s + 74);
    const auto *gd_s_75 = buffer.data(gd_s + 75);
    const auto *gd_s_76 = buffer.data(gd_s + 76);
    const auto *gd_s_77 = buffer.data(gd_s + 77);
    const auto *gd_s_78 = buffer.data(gd_s + 78);
    const auto *gd_s_79 = buffer.data(gd_s + 79);
    const auto *gd_s_80 = buffer.data(gd_s + 80);
    const auto *gd_s_81 = buffer.data(gd_s + 81);
    const auto *gd_s_82 = buffer.data(gd_s + 82);
    const auto *gd_s_83 = buffer.data(gd_s + 83);
    const auto *gd_s_84 = buffer.data(gd_s + 84);
    const auto *gd_s_85 = buffer.data(gd_s + 85);
    const auto *gd_s_86 = buffer.data(gd_s + 86);
    const auto *gd_s_87 = buffer.data(gd_s + 87);
    const auto *gd_s_88 = buffer.data(gd_s + 88);
    const auto *gd_s_89 = buffer.data(gd_s + 89);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_11 = buffer.data(gs + 11);

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
    const auto *gp_26 = buffer.data(gp + 26);
    const auto *gp_27 = buffer.data(gp + 27);
    const auto *gp_28 = buffer.data(gp + 28);
    const auto *gp_29 = buffer.data(gp + 29);
    const auto *gp_30 = buffer.data(gp + 30);
    const auto *gp_31 = buffer.data(gp + 31);
    const auto *gp_32 = buffer.data(gp + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, fp_0, gs_s_0, gd_s_0, gd_s_1, \
                         gd_s_2, gd_s_3, gs_0, gp_0, gp_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 - f_1 * gs_s_0[k]
                 + f_2 * gd_s_0[k]
                 + f_3 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = f_2 * gd_s_1[k]
                 + pb_y[k] * gp_0[k];

        t_2[k] = f_2 * gd_s_2[k]
                 + pb_z[k] * gp_0[k];

        t_3[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_3[k]
                 + f_3 * gs_0[k]
                 + pb_y[k] * gp_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_y, pb_y, pb_z, fd_0, gs_s_0, gd_s_4, gd_s_5, \
                         gd_s_6, gs_0, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * gd_s_4[k]
                 + pb_y[k] * gp_2[k];

        t_5[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_5[k]
                 + f_3 * gs_0[k]
                 + pb_z[k] * gp_2[k];

        t_6[k] = pa_y[k] * fd_0[k]
                 + f_2 * gd_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, pb_z, dd_s_1, dd_1, fp_3, fd_5, gd_s_7, \
                         gd_s_8, gd_s_9, gp_3, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_4 * fp_3[k]
                 + f_2 * gd_s_7[k]
                 + pb_x[k] * gp_4[k];

        t_8[k] = f_2 * gd_s_8[k]
                 + pb_z[k] * gp_3[k];

        t_9[k] = -f_5 * dd_s_1[k]
                 + f_6 * dd_1[k]
                 + pa_x[k] * fd_5[k]
                 + f_2 * gd_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_y, pb_z, fd_0, fd_2, gd_s_10, \
                         gd_s_11, gd_s_12, gd_s_13, gp_4, gp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * gd_s_10[k]
                  + pb_z[k] * gp_4[k];

        t_11[k] = pa_y[k] * fd_2[k]
                  + f_2 * gd_s_11[k];

        t_12[k] = pa_z[k] * fd_0[k]
                  + f_2 * gd_s_12[k];

        t_13[k] = f_2 * gd_s_13[k]
                  + pb_y[k] * gp_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_z, pb_x, pb_y, fp_4, fd_1, gd_s_14, gd_s_15, \
                         gd_s_16, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_4 * fp_4[k]
                  + f_2 * gd_s_14[k]
                  + pb_x[k] * gp_6[k];

        t_15[k] = pa_z[k] * fd_1[k]
                  + f_2 * gd_s_15[k];

        t_16[k] = f_2 * gd_s_16[k]
                  + pb_y[k] * gp_6[k];
    }

#pragma omp simd aligned(t_17, t_18, pa_x, pa_y, dd_s_0, dd_s_2, dd_0, dd_2, fd_3, fd_8, \
                         gd_s_17, gd_s_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_5 * dd_s_2[k]
                  + f_6 * dd_2[k]
                  + pa_x[k] * fd_8[k]
                  + f_2 * gd_s_17[k];

        t_18[k] = -f_7 * dd_s_0[k]
                  + f_3 * dd_0[k]
                  + pa_y[k] * fd_3[k]
                  + f_2 * gd_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pb_x, pb_z, dd_s_4, dd_4, fp_5, fd_11, \
                         gd_s_19, gd_s_20, gd_s_21, gp_7, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_6 * fp_5[k]
                  + f_2 * gd_s_19[k]
                  + pb_x[k] * gp_8[k];

        t_20[k] = f_2 * gd_s_20[k]
                  + pb_z[k] * gp_7[k];

        t_21[k] = -f_7 * dd_s_4[k]
                  + f_3 * dd_4[k]
                  + pa_x[k] * fd_11[k]
                  + f_2 * gd_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pb_z, fd_6, gs_s_3, gd_s_22, gd_s_23, \
                         gd_s_24, gs_3, gp_8, gp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_2 * gd_s_22[k]
                  + pb_z[k] * gp_8[k];

        t_23[k] = -f_1 * gs_s_3[k]
                  + f_2 * gd_s_23[k]
                  + f_3 * gs_3[k]
                  + pb_z[k] * gp_9[k];

        t_24[k] = pa_y[k] * fd_6[k]
                  + f_2 * gd_s_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, pa_z, pb_y, fp_4, fd_4, fd_5, fd_7, \
                         gd_s_25, gd_s_26, gd_s_27, gd_s_28, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_z[k] * fd_4[k]
                  + f_2 * gd_s_25[k];

        t_26[k] = pa_y[k] * fd_7[k]
                  + f_2 * gd_s_26[k];

        t_27[k] = pa_z[k] * fd_5[k]
                  + f_2 * gd_s_27[k];

        t_28[k] = f_3 * fp_4[k]
                  + f_2 * gd_s_28[k]
                  + pb_y[k] * gp_10[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pa_z, pb_y, dd_s_0, dd_0, fd_6, fd_8, \
                         gd_s_29, gd_s_30, gd_s_31, gp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * fd_8[k]
                  + f_2 * gd_s_29[k];

        t_30[k] = -f_7 * dd_s_0[k]
                  + f_3 * dd_0[k]
                  + pa_z[k] * fd_6[k]
                  + f_2 * gd_s_30[k];

        t_31[k] = f_2 * gd_s_31[k]
                  + pb_y[k] * gp_11[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, fp_6, gs_s_4, gd_s_32, gd_s_33, \
                         gd_s_34, gs_4, gp_12, gp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_6 * fp_6[k]
                  + f_2 * gd_s_32[k]
                  + pb_x[k] * gp_13[k];

        t_33[k] = -f_1 * gs_s_4[k]
                  + f_2 * gd_s_33[k]
                  + f_3 * gs_4[k]
                  + pb_y[k] * gp_12[k];

        t_34[k] = f_2 * gd_s_34[k]
                  + pb_y[k] * gp_13[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_x, pb_x, dd_s_8, dd_8, fp_7, fp_8, fd_14, fd_15, \
                         gd_s_35, gd_s_36, gd_s_37, gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -f_7 * dd_s_8[k]
                  + f_3 * dd_8[k]
                  + pa_x[k] * fd_14[k]
                  + f_2 * gd_s_35[k];

        t_36[k] = f_6 * fp_7[k]
                  + pa_x[k] * fd_15[k]
                  + f_2 * gd_s_36[k];

        t_37[k] = f_3 * fp_8[k]
                  + f_2 * gd_s_37[k]
                  + pb_x[k] * gp_15[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_x, pb_z, fd_16, fd_17, gd_s_38, gd_s_39, \
                         gd_s_40, gd_s_41, gp_14, gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_2 * gd_s_38[k]
                  + pb_z[k] * gp_14[k];

        t_39[k] = pa_x[k] * fd_16[k]
                  + f_2 * gd_s_39[k];

        t_40[k] = f_2 * gd_s_40[k]
                  + pb_z[k] * gp_15[k];

        t_41[k] = pa_x[k] * fd_17[k]
                  + f_2 * gd_s_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pa_z, pb_x, fp_10, fd_9, fd_10, fd_18, \
                         gd_s_42, gd_s_43, gd_s_44, gd_s_45, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * fd_9[k]
                  + f_2 * gd_s_42[k];

        t_43[k] = pa_z[k] * fd_10[k]
                  + f_2 * gd_s_43[k];

        t_44[k] = f_3 * fp_10[k]
                  + f_2 * gd_s_44[k]
                  + pb_x[k] * gp_16[k];

        t_45[k] = pa_x[k] * fd_18[k]
                  + f_2 * gd_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pa_y, pb_x, fp_11, fd_12, fd_19, fd_20, \
                         gd_s_46, gd_s_47, gd_s_48, gd_s_49, gp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_x[k] * fd_19[k]
                  + f_2 * gd_s_46[k];

        t_47[k] = pa_x[k] * fd_20[k]
                  + f_2 * gd_s_47[k];

        t_48[k] = pa_y[k] * fd_12[k]
                  + f_2 * gd_s_48[k];

        t_49[k] = f_3 * fp_11[k]
                  + f_2 * gd_s_49[k]
                  + pb_x[k] * gp_17[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_x, pa_y, fd_13, fd_21, fd_22, fd_23, \
                         gd_s_50, gd_s_51, gd_s_52, gd_s_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * fd_13[k]
                  + f_2 * gd_s_50[k];

        t_51[k] = pa_x[k] * fd_21[k]
                  + f_2 * gd_s_51[k];

        t_52[k] = pa_x[k] * fd_22[k]
                  + f_2 * gd_s_52[k];

        t_53[k] = pa_x[k] * fd_23[k]
                  + f_2 * gd_s_53[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_x, pb_x, pb_y, fp_13, fp_15, fd_24, gd_s_54, \
                         gd_s_55, gd_s_56, gp_18, gp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_6 * fp_13[k]
                  + pa_x[k] * fd_24[k]
                  + f_2 * gd_s_54[k];

        t_55[k] = f_2 * gd_s_55[k]
                  + pb_y[k] * gp_18[k];

        t_56[k] = f_3 * fp_15[k]
                  + f_2 * gd_s_56[k]
                  + pb_x[k] * gp_19[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_x, pb_y, fd_25, fd_26, gd_s_57, gd_s_58, \
                         gd_s_59, gp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_x[k] * fd_25[k]
                  + f_2 * gd_s_57[k];

        t_58[k] = f_2 * gd_s_58[k]
                  + pb_y[k] * gp_19[k];

        t_59[k] = pa_x[k] * fd_26[k]
                  + f_2 * gd_s_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pb_x, pb_y, fp_8, gs_s_7, gd_s_60, gd_s_61, \
                         gd_s_62, gd_s_63, gs_7, gp_20, gp_21, gp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -f_1 * gs_s_7[k]
                  + f_2 * gd_s_60[k]
                  + f_3 * gs_7[k]
                  + pb_x[k] * gp_20[k];

        t_61[k] = f_2 * gd_s_61[k]
                  + pb_x[k] * gp_21[k];

        t_62[k] = f_2 * gd_s_62[k]
                  + pb_x[k] * gp_22[k];

        t_63[k] = f_0 * fp_8[k]
                  - f_1 * gs_s_7[k]
                  + f_2 * gd_s_63[k]
                  + f_3 * gs_7[k]
                  + pb_y[k] * gp_21[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_z, pb_z, fd_15, gs_s_7, gd_s_64, gd_s_65, \
                         gd_s_66, gs_7, gp_21, gp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_2 * gd_s_64[k]
                  + pb_z[k] * gp_21[k];

        t_65[k] = -f_1 * gs_s_7[k]
                  + f_2 * gd_s_65[k]
                  + f_3 * gs_7[k]
                  + pb_z[k] * gp_22[k];

        t_66[k] = pa_z[k] * fd_15[k]
                  + f_2 * gd_s_66[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pb_x, pb_y, fp_10, fd_16, gd_s_67, \
                         gd_s_68, gd_s_69, gd_s_70, gp_23, gp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_2 * gd_s_67[k]
                  + pb_x[k] * gp_23[k];

        t_68[k] = f_2 * gd_s_68[k]
                  + pb_x[k] * gp_24[k];

        t_69[k] = pa_z[k] * fd_16[k]
                  + f_2 * gd_s_69[k];

        t_70[k] = f_4 * fp_10[k]
                  + f_2 * gd_s_70[k]
                  + pb_y[k] * gp_24[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pa_y, pb_x, dd_s_5, dd_5, fd_20, gs_s_9, gd_s_71, \
                         gd_s_72, gd_s_73, gs_9, gp_25, gp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -f_5 * dd_s_5[k]
                  + f_6 * dd_5[k]
                  + pa_y[k] * fd_20[k]
                  + f_2 * gd_s_71[k];

        t_72[k] = -f_1 * gs_s_9[k]
                  + f_2 * gd_s_72[k]
                  + f_3 * gs_9[k]
                  + pb_x[k] * gp_25[k];

        t_73[k] = f_2 * gd_s_73[k]
                  + pb_x[k] * gp_26[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pa_z, pb_x, pb_y, dd_s_4, dd_4, fp_12, fd_18, \
                         gd_s_74, gd_s_75, gd_s_76, gp_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_2 * gd_s_74[k]
                  + pb_x[k] * gp_27[k];

        t_75[k] = -f_7 * dd_s_4[k]
                  + f_3 * dd_4[k]
                  + pa_z[k] * fd_18[k]
                  + f_2 * gd_s_75[k];

        t_76[k] = f_6 * fp_12[k]
                  + f_2 * gd_s_76[k]
                  + pb_y[k] * gp_27[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pb_x, dd_s_8, dd_8, fd_23, fd_24, \
                         gd_s_77, gd_s_78, gd_s_79, gd_s_80, gp_28, \
                         gp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -f_7 * dd_s_8[k]
                  + f_3 * dd_8[k]
                  + pa_y[k] * fd_23[k]
                  + f_2 * gd_s_77[k];

        t_78[k] = pa_y[k] * fd_24[k]
                  + f_2 * gd_s_78[k];

        t_79[k] = f_2 * gd_s_79[k]
                  + pb_x[k] * gp_28[k];

        t_80[k] = f_2 * gd_s_80[k]
                  + pb_x[k] * gp_29[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_y, pb_y, fp_14, fp_15, fd_25, fd_26, gd_s_81, \
                         gd_s_82, gd_s_83, gp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_6 * fp_14[k]
                  + pa_y[k] * fd_25[k]
                  + f_2 * gd_s_81[k];

        t_82[k] = f_3 * fp_15[k]
                  + f_2 * gd_s_82[k]
                  + pb_y[k] * gp_29[k];

        t_83[k] = pa_y[k] * fd_26[k]
                  + f_2 * gd_s_83[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_x, pb_y, gs_s_11, gd_s_84, gd_s_85, \
                         gd_s_86, gd_s_87, gs_11, gp_30, gp_31, gp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -f_1 * gs_s_11[k]
                  + f_2 * gd_s_84[k]
                  + f_3 * gs_11[k]
                  + pb_x[k] * gp_30[k];

        t_85[k] = f_2 * gd_s_85[k]
                  + pb_x[k] * gp_31[k];

        t_86[k] = f_2 * gd_s_86[k]
                  + pb_x[k] * gp_32[k];

        t_87[k] = -f_1 * gs_s_11[k]
                  + f_2 * gd_s_87[k]
                  + f_3 * gs_11[k]
                  + pb_y[k] * gp_31[k];
    }

#pragma omp simd aligned(t_88, t_89, pb_y, pb_z, fp_15, gs_s_11, gd_s_88, gd_s_89, gs_11, \
                         gp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_2 * gd_s_88[k]
                  + pb_y[k] * gp_32[k];

        t_89[k] = f_0 * fp_15[k]
                  - f_1 * gs_s_11[k]
                  + f_2 * gd_s_89[k]
                  + f_3 * gs_11[k]
                  + pb_z[k] * gp_32[k];
    }
}

auto
compute_prim_gd_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dd_s, const size_t dd,
                                 const size_t fp, const size_t fd, const size_t gs_s,
                                 const size_t gd_s, const size_t gs, const size_t gp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 1.5 / p;
    const auto f_5 = 2.0 * beta / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_6 = buffer.data(dd_s + 6);
    const auto *dd_s_8 = buffer.data(dd_s + 8);
    const auto *dd_s_12 = buffer.data(dd_s + 12);
    const auto *dd_s_15 = buffer.data(dd_s + 15);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_15 = buffer.data(dd + 15);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_12 = buffer.data(fd + 12);
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

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_3 = buffer.data(gs_s + 3);
    const auto *gs_s_4 = buffer.data(gs_s + 4);
    const auto *gs_s_5 = buffer.data(gs_s + 5);
    const auto *gs_s_7 = buffer.data(gs_s + 7);
    const auto *gs_s_9 = buffer.data(gs_s + 9);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_15 = buffer.data(gd_s + 15);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_22 = buffer.data(gd_s + 22);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_24 = buffer.data(gd_s + 24);
    const auto *gd_s_25 = buffer.data(gd_s + 25);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_27 = buffer.data(gd_s + 27);
    const auto *gd_s_28 = buffer.data(gd_s + 28);
    const auto *gd_s_29 = buffer.data(gd_s + 29);
    const auto *gd_s_30 = buffer.data(gd_s + 30);
    const auto *gd_s_31 = buffer.data(gd_s + 31);
    const auto *gd_s_32 = buffer.data(gd_s + 32);
    const auto *gd_s_33 = buffer.data(gd_s + 33);
    const auto *gd_s_34 = buffer.data(gd_s + 34);
    const auto *gd_s_35 = buffer.data(gd_s + 35);
    const auto *gd_s_36 = buffer.data(gd_s + 36);
    const auto *gd_s_37 = buffer.data(gd_s + 37);
    const auto *gd_s_38 = buffer.data(gd_s + 38);
    const auto *gd_s_39 = buffer.data(gd_s + 39);
    const auto *gd_s_40 = buffer.data(gd_s + 40);
    const auto *gd_s_41 = buffer.data(gd_s + 41);
    const auto *gd_s_42 = buffer.data(gd_s + 42);
    const auto *gd_s_43 = buffer.data(gd_s + 43);
    const auto *gd_s_44 = buffer.data(gd_s + 44);
    const auto *gd_s_45 = buffer.data(gd_s + 45);
    const auto *gd_s_46 = buffer.data(gd_s + 46);
    const auto *gd_s_47 = buffer.data(gd_s + 47);
    const auto *gd_s_48 = buffer.data(gd_s + 48);
    const auto *gd_s_49 = buffer.data(gd_s + 49);
    const auto *gd_s_50 = buffer.data(gd_s + 50);
    const auto *gd_s_51 = buffer.data(gd_s + 51);
    const auto *gd_s_52 = buffer.data(gd_s + 52);
    const auto *gd_s_54 = buffer.data(gd_s + 54);
    const auto *gd_s_55 = buffer.data(gd_s + 55);
    const auto *gd_s_56 = buffer.data(gd_s + 56);
    const auto *gd_s_57 = buffer.data(gd_s + 57);
    const auto *gd_s_58 = buffer.data(gd_s + 58);
    const auto *gd_s_59 = buffer.data(gd_s + 59);
    const auto *gd_s_60 = buffer.data(gd_s + 60);
    const auto *gd_s_61 = buffer.data(gd_s + 61);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fp_0, gs_s_0, gd_s_0, gd_s_1, \
                         gd_s_2, gs_0, gp_0, gp_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 - f_1 * gs_s_0[k]
                 + f_2 * gd_s_0[k]
                 + f_3 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = f_2 * gd_s_1[k]
                 + pb_z[k] * gp_0[k];

        t_2[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_2[k]
                 + f_3 * gs_0[k]
                 + pb_y[k] * gp_1[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_y, pb_x, pb_z, fp_3, fd_0, gs_s_0, gd_s_3, gd_s_4, \
                         gd_s_5, gs_0, gp_2, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_3[k]
                 + f_3 * gs_0[k]
                 + pb_z[k] * gp_2[k];

        t_4[k] = pa_y[k] * fd_0[k]
                 + f_2 * gd_s_4[k];

        t_5[k] = f_4 * fp_3[k]
                 + f_2 * gd_s_5[k]
                 + pb_x[k] * gp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, pa_z, dd_s_4, dd_4, fd_0, fd_2, fd_4, \
                         gd_s_6, gd_s_7, gd_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_5 * dd_s_4[k]
                 + f_6 * dd_4[k]
                 + pa_x[k] * fd_4[k]
                 + f_2 * gd_s_6[k];

        t_7[k] = pa_y[k] * fd_2[k]
                 + f_2 * gd_s_7[k];

        t_8[k] = pa_z[k] * fd_0[k]
                 + f_2 * gd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_x, pb_y, dd_s_6, dd_6, fp_4, fd_7, gd_s_9, \
                         gd_s_10, gd_s_11, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * fp_4[k]
                 + f_2 * gd_s_9[k]
                 + pb_x[k] * gp_4[k];

        t_10[k] = f_2 * gd_s_10[k]
                  + pb_y[k] * gp_4[k];

        t_11[k] = -f_5 * dd_s_6[k]
                  + f_6 * dd_6[k]
                  + pa_x[k] * fd_7[k]
                  + f_2 * gd_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, pa_y, pb_x, dd_s_0, dd_0, fp_5, fd_3, gd_s_12, gd_s_13, \
                         gp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_7 * dd_s_0[k]
                  + f_3 * dd_0[k]
                  + pa_y[k] * fd_3[k]
                  + f_2 * gd_s_12[k];

        t_13[k] = f_6 * fp_5[k]
                  + f_2 * gd_s_13[k]
                  + pb_x[k] * gp_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_z, dd_s_8, dd_8, fd_4, fd_9, gs_s_3, \
                         gd_s_14, gd_s_15, gd_s_16, gs_3, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -f_7 * dd_s_8[k]
                  + f_3 * dd_8[k]
                  + pa_x[k] * fd_9[k]
                  + f_2 * gd_s_14[k];

        t_15[k] = -f_1 * gs_s_3[k]
                  + f_2 * gd_s_15[k]
                  + f_3 * gs_3[k]
                  + pb_z[k] * gp_6[k];

        t_16[k] = pa_z[k] * fd_4[k]
                  + f_2 * gd_s_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_y, pa_z, pb_y, dd_s_0, dd_0, fp_4, fd_6, fd_7, \
                         gd_s_17, gd_s_18, gd_s_19, gp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * fp_4[k]
                  + f_2 * gd_s_17[k]
                  + pb_y[k] * gp_7[k];

        t_18[k] = pa_y[k] * fd_7[k]
                  + f_2 * gd_s_18[k];

        t_19[k] = -f_7 * dd_s_0[k]
                  + f_3 * dd_0[k]
                  + pa_z[k] * fd_6[k]
                  + f_2 * gd_s_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_x, pb_y, fp_6, gs_s_4, gd_s_20, gd_s_21, \
                         gd_s_22, gs_4, gp_8, gp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_6 * fp_6[k]
                  + f_2 * gd_s_20[k]
                  + pb_x[k] * gp_9[k];

        t_21[k] = -f_1 * gs_s_4[k]
                  + f_2 * gd_s_21[k]
                  + f_3 * gs_4[k]
                  + pb_y[k] * gp_8[k];

        t_22[k] = f_2 * gd_s_22[k]
                  + pb_y[k] * gp_9[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pb_x, dd_s_15, dd_15, fp_7, fp_8, fd_14, \
                         fd_15, gd_s_23, gd_s_24, gd_s_25, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -f_7 * dd_s_15[k]
                  + f_3 * dd_15[k]
                  + pa_x[k] * fd_14[k]
                  + f_2 * gd_s_23[k];

        t_24[k] = f_6 * fp_7[k]
                  + pa_x[k] * fd_15[k]
                  + f_2 * gd_s_24[k];

        t_25[k] = f_3 * fp_8[k]
                  + f_2 * gd_s_25[k]
                  + pb_x[k] * gp_10[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_x, pa_z, fd_8, fd_17, fd_18, fd_20, \
                         gd_s_26, gd_s_27, gd_s_28, gd_s_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_x[k] * fd_17[k]
                  + f_2 * gd_s_26[k];

        t_27[k] = pa_x[k] * fd_18[k]
                  + f_2 * gd_s_27[k];

        t_28[k] = pa_z[k] * fd_8[k]
                  + f_2 * gd_s_28[k];

        t_29[k] = pa_x[k] * fd_20[k]
                  + f_2 * gd_s_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, fd_12, fd_21, fd_22, fd_23, \
                         gd_s_30, gd_s_31, gd_s_32, gd_s_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_x[k] * fd_21[k]
                  + f_2 * gd_s_30[k];

        t_31[k] = pa_y[k] * fd_12[k]
                  + f_2 * gd_s_31[k];

        t_32[k] = pa_x[k] * fd_22[k]
                  + f_2 * gd_s_32[k];

        t_33[k] = pa_x[k] * fd_23[k]
                  + f_2 * gd_s_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, fp_12, fp_14, fd_25, fd_27, \
                         fd_29, gd_s_34, gd_s_35, gd_s_36, gd_s_37, \
                         gp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_6 * fp_12[k]
                  + pa_x[k] * fd_25[k]
                  + f_2 * gd_s_34[k];

        t_35[k] = f_3 * fp_14[k]
                  + f_2 * gd_s_35[k]
                  + pb_x[k] * gp_11[k];

        t_36[k] = pa_x[k] * fd_27[k]
                  + f_2 * gd_s_36[k];

        t_37[k] = pa_x[k] * fd_29[k]
                  + f_2 * gd_s_37[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_y, pb_z, fp_8, gs_s_5, gd_s_38, \
                         gd_s_39, gd_s_40, gd_s_41, gs_5, gp_12, \
                         gp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -f_1 * gs_s_5[k]
                  + f_2 * gd_s_38[k]
                  + f_3 * gs_5[k]
                  + pb_x[k] * gp_12[k];

        t_39[k] = f_2 * gd_s_39[k]
                  + pb_x[k] * gp_13[k];

        t_40[k] = f_0 * fp_8[k]
                  - f_1 * gs_s_5[k]
                  + f_2 * gd_s_40[k]
                  + f_3 * gs_5[k]
                  + pb_y[k] * gp_13[k];

        t_41[k] = f_2 * gd_s_41[k]
                  + pb_z[k] * gp_13[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_z, pb_x, pb_z, fd_17, gs_s_5, gd_s_42, gd_s_43, \
                         gd_s_44, gs_5, gp_14, gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = -f_1 * gs_s_5[k]
                  + f_2 * gd_s_42[k]
                  + f_3 * gs_5[k]
                  + pb_z[k] * gp_14[k];

        t_43[k] = f_2 * gd_s_43[k]
                  + pb_x[k] * gp_15[k];

        t_44[k] = pa_z[k] * fd_17[k]
                  + f_2 * gd_s_44[k];
    }

#pragma omp simd aligned(t_45, t_46, pa_y, pb_y, dd_s_12, dd_12, fp_10, fd_21, gd_s_45, \
                         gd_s_46, gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_4 * fp_10[k]
                  + f_2 * gd_s_45[k]
                  + pb_y[k] * gp_15[k];

        t_46[k] = -f_5 * dd_s_12[k]
                  + f_6 * dd_12[k]
                  + pa_y[k] * fd_21[k]
                  + f_2 * gd_s_46[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pb_x, gs_s_7, gd_s_47, gd_s_48, gd_s_49, gs_7, \
                         gp_16, gp_17, gp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -f_1 * gs_s_7[k]
                  + f_2 * gd_s_47[k]
                  + f_3 * gs_7[k]
                  + pb_x[k] * gp_16[k];

        t_48[k] = f_2 * gd_s_48[k]
                  + pb_x[k] * gp_17[k];

        t_49[k] = f_2 * gd_s_49[k]
                  + pb_x[k] * gp_18[k];
    }

#pragma omp simd aligned(t_50, t_51, pa_z, pb_y, dd_s_8, dd_8, fp_11, fd_19, gd_s_50, gd_s_51, \
                         gp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -f_7 * dd_s_8[k]
                  + f_3 * dd_8[k]
                  + pa_z[k] * fd_19[k]
                  + f_2 * gd_s_50[k];

        t_51[k] = f_6 * fp_11[k]
                  + f_2 * gd_s_51[k]
                  + pb_y[k] * gp_18[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_y, pb_y, dd_s_15, dd_15, fp_13, fp_14, fd_24, \
                         fd_27, gd_s_52, gd_s_54, gd_s_55, gp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = -f_7 * dd_s_15[k]
                  + f_3 * dd_15[k]
                  + pa_y[k] * fd_24[k]
                  + f_2 * gd_s_52[k];

        t_53[k] = f_6 * fp_13[k]
                  + pa_y[k] * fd_27[k]
                  + f_2 * gd_s_54[k];

        t_54[k] = f_3 * fp_14[k]
                  + f_2 * gd_s_55[k]
                  + pb_y[k] * gp_19[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_y, pb_x, fd_29, gs_s_9, gd_s_56, gd_s_57, \
                         gd_s_58, gs_8, gp_20, gp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_y[k] * fd_29[k]
                  + f_2 * gd_s_56[k];

        t_56[k] = -f_1 * gs_s_9[k]
                  + f_2 * gd_s_57[k]
                  + f_3 * gs_8[k]
                  + pb_x[k] * gp_20[k];

        t_57[k] = f_2 * gd_s_58[k]
                  + pb_x[k] * gp_22[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pb_y, pb_z, fp_14, gs_s_9, gd_s_59, gd_s_60, \
                         gd_s_61, gs_8, gp_21, gp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -f_1 * gs_s_9[k]
                  + f_2 * gd_s_59[k]
                  + f_3 * gs_8[k]
                  + pb_y[k] * gp_21[k];

        t_59[k] = f_2 * gd_s_60[k]
                  + pb_y[k] * gp_22[k];

        t_60[k] = f_0 * fp_14[k]
                  - f_1 * gs_s_9[k]
                  + f_2 * gd_s_61[k]
                  + f_3 * gs_8[k]
                  + pb_z[k] * gp_22[k];
    }
}

auto
compute_prim_gd_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dd_s, const size_t dd,
                                 const size_t fp, const size_t fd, const size_t gs_s,
                                 const size_t gd_s, const size_t gs, const size_t gp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 * beta / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_3 = buffer.data(dd_s + 3);
    const auto *dd_s_5 = buffer.data(dd_s + 5);
    const auto *dd_s_7 = buffer.data(dd_s + 7);
    const auto *dd_s_10 = buffer.data(dd_s + 10);
    const auto *dd_s_14 = buffer.data(dd_s + 14);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_14 = buffer.data(dd + 14);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_23 = buffer.data(fd + 23);
    const auto *fd_25 = buffer.data(fd + 25);

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_2 = buffer.data(gs_s + 2);
    const auto *gs_s_3 = buffer.data(gs_s + 3);
    const auto *gs_s_4 = buffer.data(gs_s + 4);
    const auto *gs_s_6 = buffer.data(gs_s + 6);
    const auto *gs_s_8 = buffer.data(gs_s + 8);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_15 = buffer.data(gd_s + 15);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_22 = buffer.data(gd_s + 22);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_24 = buffer.data(gd_s + 24);
    const auto *gd_s_25 = buffer.data(gd_s + 25);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_27 = buffer.data(gd_s + 27);
    const auto *gd_s_28 = buffer.data(gd_s + 28);
    const auto *gd_s_30 = buffer.data(gd_s + 30);
    const auto *gd_s_31 = buffer.data(gd_s + 31);
    const auto *gd_s_32 = buffer.data(gd_s + 32);
    const auto *gd_s_33 = buffer.data(gd_s + 33);
    const auto *gd_s_34 = buffer.data(gd_s + 34);
    const auto *gd_s_35 = buffer.data(gd_s + 35);
    const auto *gd_s_36 = buffer.data(gd_s + 36);
    const auto *gd_s_38 = buffer.data(gd_s + 38);
    const auto *gd_s_40 = buffer.data(gd_s + 40);
    const auto *gd_s_41 = buffer.data(gd_s + 41);
    const auto *gd_s_42 = buffer.data(gd_s + 42);
    const auto *gd_s_43 = buffer.data(gd_s + 43);
    const auto *gd_s_44 = buffer.data(gd_s + 44);
    const auto *gd_s_45 = buffer.data(gd_s + 45);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fp_0, gs_s_0, gd_s_0, gd_s_1, \
                         gd_s_2, gs_0, gp_0, gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 - f_1 * gs_s_0[k]
                 + f_2 * gd_s_0[k]
                 + f_3 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_1[k]
                 + f_3 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_2[k]
                 + f_3 * gs_0[k]
                 + pb_z[k] * gp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pa_x, pa_y, pa_z, dd_s_3, dd_3, fd_0, fd_2, fd_4, \
                         gd_s_3, gd_s_4, gd_s_5, gd_s_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * fd_0[k]
                 + f_2 * gd_s_3[k];

        t_4[k] = -f_4 * dd_s_3[k]
                 + f_5 * dd_3[k]
                 + pa_x[k] * fd_4[k]
                 + f_2 * gd_s_4[k];

        t_5[k] = pa_y[k] * fd_2[k]
                 + f_2 * gd_s_5[k];

        t_6[k] = pa_z[k] * fd_0[k]
                 + f_2 * gd_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_y, pb_y, dd_s_0, dd_s_5, dd_0, dd_5, fd_3, \
                         fd_8, gd_s_7, gd_s_8, gd_s_9, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * gd_s_7[k]
                 + pb_y[k] * gp_3[k];

        t_8[k] = -f_4 * dd_s_5[k]
                 + f_5 * dd_5[k]
                 + pa_x[k] * fd_8[k]
                 + f_2 * gd_s_8[k];

        t_9[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_y[k] * fd_3[k]
                 + f_2 * gd_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_z, pb_z, dd_s_7, dd_7, fd_4, fd_10, \
                         gs_s_2, gd_s_10, gd_s_11, gd_s_12, gs_2, \
                         gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_6 * dd_s_7[k]
                  + f_3 * dd_7[k]
                  + pa_x[k] * fd_10[k]
                  + f_2 * gd_s_10[k];

        t_11[k] = -f_1 * gs_s_2[k]
                  + f_2 * gd_s_11[k]
                  + f_3 * gs_2[k]
                  + pb_z[k] * gp_4[k];

        t_12[k] = pa_z[k] * fd_4[k]
                  + f_2 * gd_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pa_z, pb_y, dd_s_0, dd_0, fd_6, fd_8, gs_s_3, \
                         gd_s_13, gd_s_14, gd_s_15, gs_3, gp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * fd_8[k]
                  + f_2 * gd_s_13[k];

        t_14[k] = -f_6 * dd_s_0[k]
                  + f_3 * dd_0[k]
                  + pa_z[k] * fd_6[k]
                  + f_2 * gd_s_14[k];

        t_15[k] = -f_1 * gs_s_3[k]
                  + f_2 * gd_s_15[k]
                  + f_3 * gs_3[k]
                  + pb_y[k] * gp_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pb_y, dd_s_14, dd_14, fp_4, fd_12, fd_13, \
                         gd_s_16, gd_s_17, gd_s_18, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * gd_s_16[k]
                  + pb_y[k] * gp_6[k];

        t_17[k] = -f_6 * dd_s_14[k]
                  + f_3 * dd_14[k]
                  + pa_x[k] * fd_12[k]
                  + f_2 * gd_s_17[k];

        t_18[k] = f_5 * fp_4[k]
                  + pa_x[k] * fd_13[k]
                  + f_2 * gd_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pa_z, fp_8, fd_9, fd_15, fd_21, fd_25, \
                         gd_s_19, gd_s_20, gd_s_21, gd_s_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_x[k] * fd_15[k]
                  + f_2 * gd_s_19[k];

        t_20[k] = pa_z[k] * fd_9[k]
                  + f_2 * gd_s_20[k];

        t_21[k] = f_5 * fp_8[k]
                  + pa_x[k] * fd_21[k]
                  + f_2 * gd_s_21[k];

        t_22[k] = pa_x[k] * fd_25[k]
                  + f_2 * gd_s_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_x, pb_y, fp_5, gs_s_4, gd_s_23, gd_s_24, \
                         gd_s_25, gs_4, gp_7, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -f_1 * gs_s_4[k]
                  + f_2 * gd_s_23[k]
                  + f_3 * gs_4[k]
                  + pb_x[k] * gp_7[k];

        t_24[k] = f_2 * gd_s_24[k]
                  + pb_x[k] * gp_8[k];

        t_25[k] = f_0 * fp_5[k]
                  - f_1 * gs_s_4[k]
                  + f_2 * gd_s_25[k]
                  + f_3 * gs_4[k]
                  + pb_y[k] * gp_8[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_x, pb_z, fd_15, gs_s_4, gd_s_26, gd_s_27, \
                         gd_s_28, gs_4, gp_9, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -f_1 * gs_s_4[k]
                  + f_2 * gd_s_26[k]
                  + f_3 * gs_4[k]
                  + pb_z[k] * gp_9[k];

        t_27[k] = f_2 * gd_s_27[k]
                  + pb_x[k] * gp_10[k];

        t_28[k] = pa_z[k] * fd_15[k]
                  + f_2 * gd_s_28[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pb_x, dd_s_10, dd_10, fd_18, gs_s_6, gd_s_30, \
                         gd_s_31, gd_s_32, gs_6, gp_11, gp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -f_4 * dd_s_10[k]
                  + f_5 * dd_10[k]
                  + pa_y[k] * fd_18[k]
                  + f_2 * gd_s_30[k];

        t_30[k] = -f_1 * gs_s_6[k]
                  + f_2 * gd_s_31[k]
                  + f_3 * gs_6[k]
                  + pb_x[k] * gp_11[k];

        t_31[k] = f_2 * gd_s_32[k]
                  + pb_x[k] * gp_12[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_z, pb_x, pb_y, dd_s_7, dd_7, fp_7, fd_17, \
                         gd_s_33, gd_s_34, gd_s_35, gp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_2 * gd_s_33[k]
                  + pb_x[k] * gp_13[k];

        t_33[k] = -f_6 * dd_s_7[k]
                  + f_3 * dd_7[k]
                  + pa_z[k] * fd_17[k]
                  + f_2 * gd_s_34[k];

        t_34[k] = f_5 * fp_7[k]
                  + f_2 * gd_s_35[k]
                  + pb_y[k] * gp_13[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, dd_s_14, dd_14, fp_9, fd_20, fd_23, fd_25, \
                         gd_s_36, gd_s_38, gd_s_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -f_6 * dd_s_14[k]
                  + f_3 * dd_14[k]
                  + pa_y[k] * fd_20[k]
                  + f_2 * gd_s_36[k];

        t_36[k] = f_5 * fp_9[k]
                  + pa_y[k] * fd_23[k]
                  + f_2 * gd_s_38[k];

        t_37[k] = pa_y[k] * fd_25[k]
                  + f_2 * gd_s_40[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_y, gs_s_8, gd_s_41, gd_s_42, \
                         gd_s_43, gd_s_44, gs_7, gp_14, gp_15, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -f_1 * gs_s_8[k]
                  + f_2 * gd_s_41[k]
                  + f_3 * gs_7[k]
                  + pb_x[k] * gp_14[k];

        t_39[k] = f_2 * gd_s_42[k]
                  + pb_x[k] * gp_16[k];

        t_40[k] = -f_1 * gs_s_8[k]
                  + f_2 * gd_s_43[k]
                  + f_3 * gs_7[k]
                  + pb_y[k] * gp_15[k];

        t_41[k] = f_2 * gd_s_44[k]
                  + pb_y[k] * gp_16[k];
    }

#pragma omp simd aligned(t_42, pb_z, fp_10, gs_s_8, gd_s_45, gs_7, \
                         gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_0 * fp_10[k]
                  - f_1 * gs_s_8[k]
                  + f_2 * gd_s_45[k]
                  + f_3 * gs_7[k]
                  + pb_z[k] * gp_16[k];
    }
}

auto
compute_prim_gd_kinetic_energy_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dd_s, const size_t dd,
                                 const size_t fp, const size_t fd, const size_t gs_s,
                                 const size_t gd_s, const size_t gs, const size_t gp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 * beta / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;

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

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_3 = buffer.data(dd_s + 3);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_6 = buffer.data(dd_s + 6);
    const auto *dd_s_8 = buffer.data(dd_s + 8);
    const auto *dd_s_12 = buffer.data(dd_s + 12);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_12 = buffer.data(dd + 12);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_10 = buffer.data(fp + 10);

    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_17 = buffer.data(fd + 17);

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_2 = buffer.data(gs_s + 2);
    const auto *gs_s_3 = buffer.data(gs_s + 3);
    const auto *gs_s_4 = buffer.data(gs_s + 4);
    const auto *gs_s_6 = buffer.data(gs_s + 6);
    const auto *gs_s_8 = buffer.data(gs_s + 8);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_24 = buffer.data(gd_s + 24);
    const auto *gd_s_25 = buffer.data(gd_s + 25);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_27 = buffer.data(gd_s + 27);
    const auto *gd_s_28 = buffer.data(gd_s + 28);
    const auto *gd_s_29 = buffer.data(gd_s + 29);
    const auto *gd_s_34 = buffer.data(gd_s + 34);
    const auto *gd_s_35 = buffer.data(gd_s + 35);
    const auto *gd_s_36 = buffer.data(gd_s + 36);
    const auto *gd_s_37 = buffer.data(gd_s + 37);
    const auto *gd_s_38 = buffer.data(gd_s + 38);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fp_0, gs_s_0, gd_s_0, gd_s_1, \
                         gd_s_2, gs_0, gp_0, gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 - f_1 * gs_s_0[k]
                 + f_2 * gd_s_0[k]
                 + f_3 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_1[k]
                 + f_3 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_2[k]
                 + f_3 * gs_0[k]
                 + pb_z[k] * gp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_y, dd_s_3, dd_s_4, dd_3, dd_4, fd_4, fd_7, \
                         gd_s_3, gd_s_5, gd_s_6, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_4 * dd_s_3[k]
                 + f_5 * dd_3[k]
                 + pa_x[k] * fd_4[k]
                 + f_2 * gd_s_3[k];

        t_4[k] = f_2 * gd_s_5[k]
                 + pb_y[k] * gp_3[k];

        t_5[k] = -f_4 * dd_s_4[k]
                 + f_5 * dd_4[k]
                 + pa_x[k] * fd_7[k]
                 + f_2 * gd_s_6[k];
    }

#pragma omp simd aligned(t_6, t_7, pa_x, pa_y, dd_s_0, dd_s_6, dd_0, dd_6, fd_3, fd_8, gd_s_7, \
                         gd_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_y[k] * fd_3[k]
                 + f_2 * gd_s_7[k];

        t_7[k] = -f_6 * dd_s_6[k]
                 + f_3 * dd_6[k]
                 + pa_x[k] * fd_8[k]
                 + f_2 * gd_s_8[k];
    }

#pragma omp simd aligned(t_8, t_9, pa_z, pb_z, dd_s_0, dd_0, fd_5, gs_s_2, gd_s_9, gd_s_10, \
                         gs_2, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = -f_1 * gs_s_2[k]
                 + f_2 * gd_s_9[k]
                 + f_3 * gs_2[k]
                 + pb_z[k] * gp_4[k];

        t_9[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_z[k] * fd_5[k]
                 + f_2 * gd_s_10[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pb_y, dd_s_12, dd_12, fd_9, gs_s_3, gd_s_11, \
                         gd_s_12, gd_s_13, gs_3, gp_5, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * gs_s_3[k]
                  + f_2 * gd_s_11[k]
                  + f_3 * gs_3[k]
                  + pb_y[k] * gp_5[k];

        t_11[k] = f_2 * gd_s_12[k]
                  + pb_y[k] * gp_6[k];

        t_12[k] = -f_6 * dd_s_12[k]
                  + f_3 * dd_12[k]
                  + pa_x[k] * fd_9[k]
                  + f_2 * gd_s_13[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, pb_y, fp_5, gs_s_4, gd_s_16, gd_s_17, \
                         gd_s_18, gs_4, gp_7, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_1 * gs_s_4[k]
                  + f_2 * gd_s_16[k]
                  + f_3 * gs_4[k]
                  + pb_x[k] * gp_7[k];

        t_14[k] = f_2 * gd_s_17[k]
                  + pb_x[k] * gp_8[k];

        t_15[k] = f_0 * fp_5[k]
                  - f_1 * gs_s_4[k]
                  + f_2 * gd_s_18[k]
                  + f_3 * gs_4[k]
                  + pb_y[k] * gp_8[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_y, pb_x, pb_z, dd_s_8, dd_8, fd_15, gs_s_4, \
                         gd_s_19, gd_s_20, gd_s_23, gs_4, gp_9, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_1 * gs_s_4[k]
                  + f_2 * gd_s_19[k]
                  + f_3 * gs_4[k]
                  + pb_z[k] * gp_9[k];

        t_17[k] = f_2 * gd_s_20[k]
                  + pb_x[k] * gp_10[k];

        t_18[k] = -f_4 * dd_s_8[k]
                  + f_5 * dd_8[k]
                  + pa_y[k] * fd_15[k]
                  + f_2 * gd_s_23[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_x, gs_s_6, gd_s_24, gd_s_25, gd_s_26, gs_6, \
                         gp_11, gp_12, gp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_1 * gs_s_6[k]
                  + f_2 * gd_s_24[k]
                  + f_3 * gs_6[k]
                  + pb_x[k] * gp_11[k];

        t_20[k] = f_2 * gd_s_25[k]
                  + pb_x[k] * gp_12[k];

        t_21[k] = f_2 * gd_s_26[k]
                  + pb_x[k] * gp_13[k];
    }

#pragma omp simd aligned(t_22, t_23, pa_z, pb_y, dd_s_6, dd_6, fp_7, fd_14, gd_s_27, gd_s_28, \
                         gp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -f_6 * dd_s_6[k]
                  + f_3 * dd_6[k]
                  + pa_z[k] * fd_14[k]
                  + f_2 * gd_s_27[k];

        t_23[k] = f_5 * fp_7[k]
                  + f_2 * gd_s_28[k]
                  + pb_y[k] * gp_13[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_y, pb_x, dd_s_12, dd_12, fd_17, gs_s_8, gd_s_29, \
                         gd_s_34, gd_s_35, gs_7, gp_14, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_6 * dd_s_12[k]
                  + f_3 * dd_12[k]
                  + pa_y[k] * fd_17[k]
                  + f_2 * gd_s_29[k];

        t_25[k] = -f_1 * gs_s_8[k]
                  + f_2 * gd_s_34[k]
                  + f_3 * gs_7[k]
                  + pb_x[k] * gp_14[k];

        t_26[k] = f_2 * gd_s_35[k]
                  + pb_x[k] * gp_16[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pb_z, fp_10, gs_s_8, gd_s_36, gd_s_37, \
                         gd_s_38, gs_7, gp_15, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_1 * gs_s_8[k]
                  + f_2 * gd_s_36[k]
                  + f_3 * gs_7[k]
                  + pb_y[k] * gp_15[k];

        t_28[k] = f_2 * gd_s_37[k]
                  + pb_y[k] * gp_16[k];

        t_29[k] = f_0 * fp_10[k]
                  - f_1 * gs_s_8[k]
                  + f_2 * gd_s_38[k]
                  + f_3 * gs_7[k]
                  + pb_z[k] * gp_16[k];
    }
}

auto
compute_prim_gd_kinetic_energy_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dd_s, const size_t dd,
                                 const size_t fp, const size_t fd, const size_t gs_s,
                                 const size_t gd_s, const size_t gs, const size_t gp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 * beta / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;

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

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_3 = buffer.data(dd_s + 3);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_6 = buffer.data(dd_s + 6);
    const auto *dd_s_8 = buffer.data(dd_s + 8);
    const auto *dd_s_12 = buffer.data(dd_s + 12);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_12 = buffer.data(dd + 12);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_10 = buffer.data(fp + 10);

    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_2 = buffer.data(gs_s + 2);
    const auto *gs_s_3 = buffer.data(gs_s + 3);
    const auto *gs_s_4 = buffer.data(gs_s + 4);
    const auto *gs_s_6 = buffer.data(gs_s + 6);
    const auto *gs_s_8 = buffer.data(gs_s + 8);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_24 = buffer.data(gd_s + 24);
    const auto *gd_s_25 = buffer.data(gd_s + 25);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_27 = buffer.data(gd_s + 27);
    const auto *gd_s_28 = buffer.data(gd_s + 28);
    const auto *gd_s_29 = buffer.data(gd_s + 29);
    const auto *gd_s_34 = buffer.data(gd_s + 34);
    const auto *gd_s_35 = buffer.data(gd_s + 35);
    const auto *gd_s_36 = buffer.data(gd_s + 36);
    const auto *gd_s_37 = buffer.data(gd_s + 37);
    const auto *gd_s_38 = buffer.data(gd_s + 38);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fp_0, gs_s_0, gd_s_0, gd_s_1, \
                         gd_s_2, gs_0, gp_0, gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 - f_1 * gs_s_0[k]
                 + f_2 * gd_s_0[k]
                 + f_3 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_1[k]
                 + f_3 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_2[k]
                 + f_3 * gs_0[k]
                 + pb_z[k] * gp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_y, dd_s_3, dd_s_4, dd_3, dd_4, fd_4, fd_7, \
                         gd_s_3, gd_s_5, gd_s_6, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_4 * dd_s_3[k]
                 + f_5 * dd_3[k]
                 + pa_x[k] * fd_4[k]
                 + f_2 * gd_s_3[k];

        t_4[k] = f_2 * gd_s_5[k]
                 + pb_y[k] * gp_3[k];

        t_5[k] = -f_4 * dd_s_4[k]
                 + f_5 * dd_4[k]
                 + pa_x[k] * fd_7[k]
                 + f_2 * gd_s_6[k];
    }

#pragma omp simd aligned(t_6, t_7, pa_x, pa_y, dd_s_0, dd_s_6, dd_0, dd_6, fd_3, fd_8, gd_s_7, \
                         gd_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_y[k] * fd_3[k]
                 + f_2 * gd_s_7[k];

        t_7[k] = -f_6 * dd_s_6[k]
                 + f_3 * dd_6[k]
                 + pa_x[k] * fd_8[k]
                 + f_2 * gd_s_8[k];
    }

#pragma omp simd aligned(t_8, t_9, pa_z, pb_z, dd_s_0, dd_0, fd_5, gs_s_2, gd_s_9, gd_s_10, \
                         gs_2, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = -f_1 * gs_s_2[k]
                 + f_2 * gd_s_9[k]
                 + f_3 * gs_2[k]
                 + pb_z[k] * gp_4[k];

        t_9[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_z[k] * fd_5[k]
                 + f_2 * gd_s_10[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pb_y, dd_s_12, dd_12, fd_9, gs_s_3, gd_s_11, \
                         gd_s_12, gd_s_13, gs_3, gp_5, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * gs_s_3[k]
                  + f_2 * gd_s_11[k]
                  + f_3 * gs_3[k]
                  + pb_y[k] * gp_5[k];

        t_11[k] = f_2 * gd_s_12[k]
                  + pb_y[k] * gp_6[k];

        t_12[k] = -f_6 * dd_s_12[k]
                  + f_3 * dd_12[k]
                  + pa_x[k] * fd_9[k]
                  + f_2 * gd_s_13[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, pb_y, fp_5, gs_s_4, gd_s_16, gd_s_17, \
                         gd_s_18, gs_4, gp_7, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_1 * gs_s_4[k]
                  + f_2 * gd_s_16[k]
                  + f_3 * gs_4[k]
                  + pb_x[k] * gp_7[k];

        t_14[k] = f_2 * gd_s_17[k]
                  + pb_x[k] * gp_8[k];

        t_15[k] = f_0 * fp_5[k]
                  - f_1 * gs_s_4[k]
                  + f_2 * gd_s_18[k]
                  + f_3 * gs_4[k]
                  + pb_y[k] * gp_8[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_y, pb_x, pb_z, dd_s_8, dd_8, fd_15, gs_s_4, \
                         gd_s_19, gd_s_20, gd_s_23, gs_4, gp_9, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_1 * gs_s_4[k]
                  + f_2 * gd_s_19[k]
                  + f_3 * gs_4[k]
                  + pb_z[k] * gp_9[k];

        t_17[k] = f_2 * gd_s_20[k]
                  + pb_x[k] * gp_10[k];

        t_18[k] = -f_4 * dd_s_8[k]
                  + f_5 * dd_8[k]
                  + pa_y[k] * fd_15[k]
                  + f_2 * gd_s_23[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_x, gs_s_6, gd_s_24, gd_s_25, gd_s_26, gs_6, \
                         gp_11, gp_12, gp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_1 * gs_s_6[k]
                  + f_2 * gd_s_24[k]
                  + f_3 * gs_6[k]
                  + pb_x[k] * gp_11[k];

        t_20[k] = f_2 * gd_s_25[k]
                  + pb_x[k] * gp_12[k];

        t_21[k] = f_2 * gd_s_26[k]
                  + pb_x[k] * gp_13[k];
    }

#pragma omp simd aligned(t_22, t_23, pa_z, pb_y, dd_s_6, dd_6, fp_7, fd_14, gd_s_27, gd_s_28, \
                         gp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -f_6 * dd_s_6[k]
                  + f_3 * dd_6[k]
                  + pa_z[k] * fd_14[k]
                  + f_2 * gd_s_27[k];

        t_23[k] = f_5 * fp_7[k]
                  + f_2 * gd_s_28[k]
                  + pb_y[k] * gp_13[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_y, pb_x, dd_s_12, dd_12, fd_16, gs_s_8, gd_s_29, \
                         gd_s_34, gd_s_35, gs_7, gp_14, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_6 * dd_s_12[k]
                  + f_3 * dd_12[k]
                  + pa_y[k] * fd_16[k]
                  + f_2 * gd_s_29[k];

        t_25[k] = -f_1 * gs_s_8[k]
                  + f_2 * gd_s_34[k]
                  + f_3 * gs_7[k]
                  + pb_x[k] * gp_14[k];

        t_26[k] = f_2 * gd_s_35[k]
                  + pb_x[k] * gp_16[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pb_z, fp_10, gs_s_8, gd_s_36, gd_s_37, \
                         gd_s_38, gs_7, gp_15, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_1 * gs_s_8[k]
                  + f_2 * gd_s_36[k]
                  + f_3 * gs_7[k]
                  + pb_y[k] * gp_15[k];

        t_28[k] = f_2 * gd_s_37[k]
                  + pb_y[k] * gp_16[k];

        t_29[k] = f_0 * fp_10[k]
                  - f_1 * gs_s_8[k]
                  + f_2 * gd_s_38[k]
                  + f_3 * gs_7[k]
                  + pb_z[k] * gp_16[k];
    }
}

auto
compute_prim_gd_kinetic_energy_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dd_s, const size_t dd,
                                 const size_t fp, const size_t fd, const size_t gs_s,
                                 const size_t gd_s, const size_t gs, const size_t gp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 1.5 / p;
    const auto f_5 = 2.0 * beta / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = beta / p;

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

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_3 = buffer.data(dd_s + 3);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_6 = buffer.data(dd_s + 6);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_12 = buffer.data(fp + 12);

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

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_7 = buffer.data(gs_s + 7);
    const auto *gs_s_9 = buffer.data(gs_s + 9);
    const auto *gs_s_11 = buffer.data(gs_s + 11);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_15 = buffer.data(gd_s + 15);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_22 = buffer.data(gd_s + 22);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_24 = buffer.data(gd_s + 24);
    const auto *gd_s_25 = buffer.data(gd_s + 25);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_27 = buffer.data(gd_s + 27);
    const auto *gd_s_28 = buffer.data(gd_s + 28);
    const auto *gd_s_29 = buffer.data(gd_s + 29);
    const auto *gd_s_30 = buffer.data(gd_s + 30);
    const auto *gd_s_31 = buffer.data(gd_s + 31);
    const auto *gd_s_32 = buffer.data(gd_s + 32);
    const auto *gd_s_33 = buffer.data(gd_s + 33);
    const auto *gd_s_34 = buffer.data(gd_s + 34);
    const auto *gd_s_35 = buffer.data(gd_s + 35);
    const auto *gd_s_36 = buffer.data(gd_s + 36);
    const auto *gd_s_37 = buffer.data(gd_s + 37);
    const auto *gd_s_38 = buffer.data(gd_s + 38);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_11 = buffer.data(gs + 11);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);
    const auto *gp_22 = buffer.data(gp + 22);
    const auto *gp_23 = buffer.data(gp + 23);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fp_0, gs_s_0, gd_s_0, gd_s_1, \
                         gd_s_2, gs_0, gp_0, gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 - f_1 * gs_s_0[k]
                 + f_2 * gd_s_0[k]
                 + f_3 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_1[k]
                 + f_3 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_2[k]
                 + f_3 * gs_0[k]
                 + pb_z[k] * gp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_x, dd_s_1, dd_1, fp_1, fd_0, fd_2, \
                         gd_s_3, gd_s_4, gd_s_5, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * fd_0[k]
                 + f_2 * gd_s_3[k];

        t_4[k] = f_4 * fp_1[k]
                 + f_2 * gd_s_4[k]
                 + pb_x[k] * gp_3[k];

        t_5[k] = -f_5 * dd_s_1[k]
                 + f_6 * dd_1[k]
                 + pa_x[k] * fd_2[k]
                 + f_2 * gd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_x, dd_s_2, dd_2, fp_2, fd_0, fd_4, \
                         gd_s_6, gd_s_7, gd_s_8, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_z[k] * fd_0[k]
                 + f_2 * gd_s_6[k];

        t_7[k] = f_4 * fp_2[k]
                 + f_2 * gd_s_7[k]
                 + pb_x[k] * gp_4[k];

        t_8[k] = -f_5 * dd_s_2[k]
                 + f_6 * dd_2[k]
                 + pa_x[k] * fd_4[k]
                 + f_2 * gd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, pa_y, pb_x, dd_s_0, dd_0, fp_3, fd_1, gd_s_9, gd_s_10, \
                         gp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_7 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_y[k] * fd_1[k]
                 + f_2 * gd_s_9[k];

        t_10[k] = f_6 * fp_3[k]
                  + f_2 * gd_s_10[k]
                  + pb_x[k] * gp_5[k];
    }

#pragma omp simd aligned(t_11, t_12, pa_x, pa_z, dd_s_0, dd_s_3, dd_0, dd_3, fd_3, fd_5, \
                         gd_s_11, gd_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = -f_7 * dd_s_3[k]
                  + f_3 * dd_3[k]
                  + pa_x[k] * fd_5[k]
                  + f_2 * gd_s_11[k];

        t_12[k] = -f_7 * dd_s_0[k]
                  + f_3 * dd_0[k]
                  + pa_z[k] * fd_3[k]
                  + f_2 * gd_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pb_x, dd_s_6, dd_6, fp_4, fp_5, fd_6, fd_7, \
                         gd_s_13, gd_s_14, gd_s_15, gp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_6 * fp_4[k]
                  + f_2 * gd_s_13[k]
                  + pb_x[k] * gp_7[k];

        t_14[k] = -f_7 * dd_s_6[k]
                  + f_3 * dd_6[k]
                  + pa_x[k] * fd_6[k]
                  + f_2 * gd_s_14[k];

        t_15[k] = f_6 * fp_5[k]
                  + pa_x[k] * fd_7[k]
                  + f_2 * gd_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_x, fp_6, fd_8, fd_10, fd_11, \
                         gd_s_16, gd_s_17, gd_s_18, gd_s_19, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * fp_6[k]
                  + f_2 * gd_s_16[k]
                  + pb_x[k] * gp_8[k];

        t_17[k] = pa_x[k] * fd_8[k]
                  + f_2 * gd_s_17[k];

        t_18[k] = pa_x[k] * fd_10[k]
                  + f_2 * gd_s_18[k];

        t_19[k] = pa_x[k] * fd_11[k]
                  + f_2 * gd_s_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_x, pb_x, fp_10, fp_12, fd_13, fd_15, gd_s_20, \
                         gd_s_21, gd_s_22, gp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_6 * fp_10[k]
                  + pa_x[k] * fd_13[k]
                  + f_2 * gd_s_20[k];

        t_21[k] = f_3 * fp_12[k]
                  + f_2 * gd_s_21[k]
                  + pb_x[k] * gp_11[k];

        t_22[k] = pa_x[k] * fd_15[k]
                  + f_2 * gd_s_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_x, pb_y, pb_z, fp_6, gs_s_7, gd_s_23, gd_s_24, \
                         gd_s_25, gs_7, gp_12, gp_13, gp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -f_1 * gs_s_7[k]
                  + f_2 * gd_s_23[k]
                  + f_3 * gs_7[k]
                  + pb_x[k] * gp_12[k];

        t_24[k] = f_0 * fp_6[k]
                  - f_1 * gs_s_7[k]
                  + f_2 * gd_s_24[k]
                  + f_3 * gs_7[k]
                  + pb_y[k] * gp_13[k];

        t_25[k] = -f_1 * gs_s_7[k]
                  + f_2 * gd_s_25[k]
                  + f_3 * gs_7[k]
                  + pb_z[k] * gp_14[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pa_z, pb_y, dd_s_4, dd_4, fp_7, fd_8, fd_10, \
                         gd_s_26, gd_s_27, gd_s_28, gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_z[k] * fd_8[k]
                  + f_2 * gd_s_26[k];

        t_27[k] = f_4 * fp_7[k]
                  + f_2 * gd_s_27[k]
                  + pb_y[k] * gp_15[k];

        t_28[k] = -f_5 * dd_s_4[k]
                  + f_6 * dd_4[k]
                  + pa_y[k] * fd_10[k]
                  + f_2 * gd_s_28[k];
    }

#pragma omp simd aligned(t_29, t_30, pa_z, pb_x, dd_s_3, dd_3, fd_9, gs_s_9, gd_s_29, gd_s_30, \
                         gs_9, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -f_1 * gs_s_9[k]
                  + f_2 * gd_s_29[k]
                  + f_3 * gs_9[k]
                  + pb_x[k] * gp_16[k];

        t_30[k] = -f_7 * dd_s_3[k]
                  + f_3 * dd_3[k]
                  + pa_z[k] * fd_9[k]
                  + f_2 * gd_s_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pb_y, dd_s_6, dd_6, fp_9, fp_11, fd_12, \
                         fd_14, gd_s_31, gd_s_32, gd_s_33, gp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_6 * fp_9[k]
                  + f_2 * gd_s_31[k]
                  + pb_y[k] * gp_18[k];

        t_32[k] = -f_7 * dd_s_6[k]
                  + f_3 * dd_6[k]
                  + pa_y[k] * fd_12[k]
                  + f_2 * gd_s_32[k];

        t_33[k] = f_6 * fp_11[k]
                  + pa_y[k] * fd_14[k]
                  + f_2 * gd_s_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pb_x, pb_y, fp_12, fd_15, gs_s_11, gd_s_34, \
                         gd_s_35, gd_s_36, gs_11, gp_20, gp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * fp_12[k]
                  + f_2 * gd_s_34[k]
                  + pb_y[k] * gp_20[k];

        t_35[k] = pa_y[k] * fd_15[k]
                  + f_2 * gd_s_35[k];

        t_36[k] = -f_1 * gs_s_11[k]
                  + f_2 * gd_s_36[k]
                  + f_3 * gs_11[k]
                  + pb_x[k] * gp_21[k];
    }

#pragma omp simd aligned(t_37, t_38, pb_y, pb_z, fp_12, gs_s_11, gd_s_37, gd_s_38, gs_11, \
                         gp_22, gp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -f_1 * gs_s_11[k]
                  + f_2 * gd_s_37[k]
                  + f_3 * gs_11[k]
                  + pb_y[k] * gp_22[k];

        t_38[k] = f_0 * fp_12[k]
                  - f_1 * gs_s_11[k]
                  + f_2 * gd_s_38[k]
                  + f_3 * gs_11[k]
                  + pb_z[k] * gp_23[k];
    }
}

auto
compute_prim_gd_kinetic_energy_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dd_s, const size_t dd,
                                 const size_t fp, const size_t fd, const size_t gs_s,
                                 const size_t gd_s, const size_t gs, const size_t gp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 * beta / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 1.5 / p;

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

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_5 = buffer.data(dd_s + 5);
    const auto *dd_s_8 = buffer.data(dd_s + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);

    const auto *fd_0 = buffer.data(fd + 0);
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
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_22 = buffer.data(fd + 22);

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_3 = buffer.data(gs_s + 3);
    const auto *gs_s_4 = buffer.data(gs_s + 4);
    const auto *gs_s_7 = buffer.data(gs_s + 7);
    const auto *gs_s_9 = buffer.data(gs_s + 9);
    const auto *gs_s_11 = buffer.data(gs_s + 11);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_15 = buffer.data(gd_s + 15);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_22 = buffer.data(gd_s + 22);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_24 = buffer.data(gd_s + 24);
    const auto *gd_s_25 = buffer.data(gd_s + 25);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_27 = buffer.data(gd_s + 27);
    const auto *gd_s_28 = buffer.data(gd_s + 28);
    const auto *gd_s_29 = buffer.data(gd_s + 29);
    const auto *gd_s_30 = buffer.data(gd_s + 30);
    const auto *gd_s_31 = buffer.data(gd_s + 31);
    const auto *gd_s_32 = buffer.data(gd_s + 32);
    const auto *gd_s_33 = buffer.data(gd_s + 33);
    const auto *gd_s_34 = buffer.data(gd_s + 34);
    const auto *gd_s_35 = buffer.data(gd_s + 35);
    const auto *gd_s_36 = buffer.data(gd_s + 36);
    const auto *gd_s_37 = buffer.data(gd_s + 37);
    const auto *gd_s_38 = buffer.data(gd_s + 38);
    const auto *gd_s_39 = buffer.data(gd_s + 39);
    const auto *gd_s_40 = buffer.data(gd_s + 40);
    const auto *gd_s_41 = buffer.data(gd_s + 41);
    const auto *gd_s_42 = buffer.data(gd_s + 42);
    const auto *gd_s_43 = buffer.data(gd_s + 43);
    const auto *gd_s_44 = buffer.data(gd_s + 44);
    const auto *gd_s_45 = buffer.data(gd_s + 45);
    const auto *gd_s_46 = buffer.data(gd_s + 46);
    const auto *gd_s_47 = buffer.data(gd_s + 47);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_11 = buffer.data(gs + 11);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fp_0, gs_s_0, gd_s_0, gd_s_1, \
                         gd_s_2, gs_0, gp_0, gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 - f_1 * gs_s_0[k]
                 + f_2 * gd_s_0[k]
                 + f_3 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_1[k]
                 + f_3 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_2[k]
                 + f_3 * gs_0[k]
                 + pb_z[k] * gp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pa_x, pa_y, pa_z, dd_s_1, dd_1, fd_0, fd_2, fd_4, \
                         gd_s_3, gd_s_4, gd_s_5, gd_s_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * fd_0[k]
                 + f_2 * gd_s_3[k];

        t_4[k] = -f_4 * dd_s_1[k]
                 + f_5 * dd_1[k]
                 + pa_x[k] * fd_4[k]
                 + f_2 * gd_s_4[k];

        t_5[k] = pa_y[k] * fd_2[k]
                 + f_2 * gd_s_5[k];

        t_6[k] = pa_z[k] * fd_0[k]
                 + f_2 * gd_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pa_x, pa_y, dd_s_0, dd_s_2, dd_0, dd_2, fd_3, fd_6, gd_s_7, \
                         gd_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -f_4 * dd_s_2[k]
                 + f_5 * dd_2[k]
                 + pa_x[k] * fd_6[k]
                 + f_2 * gd_s_7[k];

        t_8[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_y[k] * fd_3[k]
                 + f_2 * gd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_z, pb_z, dd_s_4, dd_4, fd_4, fd_8, gs_s_3, \
                         gd_s_9, gd_s_10, gd_s_11, gs_3, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_6 * dd_s_4[k]
                 + f_3 * dd_4[k]
                 + pa_x[k] * fd_8[k]
                 + f_2 * gd_s_9[k];

        t_10[k] = -f_1 * gs_s_3[k]
                  + f_2 * gd_s_10[k]
                  + f_3 * gs_3[k]
                  + pb_z[k] * gp_6[k];

        t_11[k] = pa_z[k] * fd_4[k]
                  + f_2 * gd_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pa_z, pb_y, dd_s_0, dd_0, fp_3, fd_5, fd_6, \
                         gd_s_12, gd_s_13, gd_s_14, gp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * fp_3[k]
                  + f_2 * gd_s_12[k]
                  + pb_y[k] * gp_7[k];

        t_13[k] = pa_y[k] * fd_6[k]
                  + f_2 * gd_s_13[k];

        t_14[k] = -f_6 * dd_s_0[k]
                  + f_3 * dd_0[k]
                  + pa_z[k] * fd_5[k]
                  + f_2 * gd_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pb_y, dd_s_8, dd_8, fp_4, fd_10, fd_11, \
                         gs_s_4, gd_s_15, gd_s_16, gd_s_17, gs_4, \
                         gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_1 * gs_s_4[k]
                  + f_2 * gd_s_15[k]
                  + f_3 * gs_4[k]
                  + pb_y[k] * gp_8[k];

        t_16[k] = -f_6 * dd_s_8[k]
                  + f_3 * dd_8[k]
                  + pa_x[k] * fd_10[k]
                  + f_2 * gd_s_16[k];

        t_17[k] = f_5 * fp_4[k]
                  + pa_x[k] * fd_11[k]
                  + f_2 * gd_s_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_z, fd_7, fd_12, fd_13, fd_15, \
                         gd_s_18, gd_s_19, gd_s_20, gd_s_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pa_x[k] * fd_12[k]
                  + f_2 * gd_s_18[k];

        t_19[k] = pa_x[k] * fd_13[k]
                  + f_2 * gd_s_19[k];

        t_20[k] = pa_z[k] * fd_7[k]
                  + f_2 * gd_s_20[k];

        t_21[k] = pa_x[k] * fd_15[k]
                  + f_2 * gd_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pa_y, fd_9, fd_16, fd_17, fd_18, \
                         gd_s_22, gd_s_23, gd_s_24, gd_s_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pa_x[k] * fd_16[k]
                  + f_2 * gd_s_22[k];

        t_23[k] = pa_y[k] * fd_9[k]
                  + f_2 * gd_s_23[k];

        t_24[k] = pa_x[k] * fd_17[k]
                  + f_2 * gd_s_24[k];

        t_25[k] = pa_x[k] * fd_18[k]
                  + f_2 * gd_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_x, fp_9, fd_20, fd_21, fd_22, gd_s_26, gd_s_27, \
                         gd_s_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_5 * fp_9[k]
                  + pa_x[k] * fd_20[k]
                  + f_2 * gd_s_26[k];

        t_27[k] = pa_x[k] * fd_21[k]
                  + f_2 * gd_s_27[k];

        t_28[k] = pa_x[k] * fd_22[k]
                  + f_2 * gd_s_28[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, pb_y, fp_5, gs_s_7, gd_s_29, gd_s_30, \
                         gd_s_31, gs_7, gp_12, gp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -f_1 * gs_s_7[k]
                  + f_2 * gd_s_29[k]
                  + f_3 * gs_7[k]
                  + pb_x[k] * gp_12[k];

        t_30[k] = f_2 * gd_s_30[k]
                  + pb_x[k] * gp_13[k];

        t_31[k] = f_0 * fp_5[k]
                  - f_1 * gs_s_7[k]
                  + f_2 * gd_s_31[k]
                  + f_3 * gs_7[k]
                  + pb_y[k] * gp_13[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_z, pb_y, pb_z, fp_7, fd_12, gs_s_7, gd_s_32, \
                         gd_s_33, gd_s_34, gs_7, gp_14, gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -f_1 * gs_s_7[k]
                  + f_2 * gd_s_32[k]
                  + f_3 * gs_7[k]
                  + pb_z[k] * gp_14[k];

        t_33[k] = pa_z[k] * fd_12[k]
                  + f_2 * gd_s_33[k];

        t_34[k] = f_7 * fp_7[k]
                  + f_2 * gd_s_34[k]
                  + pb_y[k] * gp_15[k];
    }

#pragma omp simd aligned(t_35, t_36, pa_y, pb_x, dd_s_5, dd_5, fd_16, gs_s_9, gd_s_35, \
                         gd_s_36, gs_9, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -f_4 * dd_s_5[k]
                  + f_5 * dd_5[k]
                  + pa_y[k] * fd_16[k]
                  + f_2 * gd_s_35[k];

        t_36[k] = -f_1 * gs_s_9[k]
                  + f_2 * gd_s_36[k]
                  + f_3 * gs_9[k]
                  + pb_x[k] * gp_16[k];
    }

#pragma omp simd aligned(t_37, t_38, pa_z, pb_y, dd_s_4, dd_4, fp_8, fd_14, gd_s_37, gd_s_38, \
                         gp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -f_6 * dd_s_4[k]
                  + f_3 * dd_4[k]
                  + pa_z[k] * fd_14[k]
                  + f_2 * gd_s_37[k];

        t_38[k] = f_5 * fp_8[k]
                  + f_2 * gd_s_38[k]
                  + pb_y[k] * gp_17[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_y, dd_s_8, dd_8, fp_10, fp_11, fd_19, \
                         fd_21, gd_s_39, gd_s_40, gd_s_41, gp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -f_6 * dd_s_8[k]
                  + f_3 * dd_8[k]
                  + pa_y[k] * fd_19[k]
                  + f_2 * gd_s_39[k];

        t_40[k] = f_5 * fp_10[k]
                  + pa_y[k] * fd_21[k]
                  + f_2 * gd_s_40[k];

        t_41[k] = f_3 * fp_11[k]
                  + f_2 * gd_s_41[k]
                  + pb_y[k] * gp_18[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_y, pb_x, fd_22, gs_s_11, gd_s_42, gd_s_43, \
                         gd_s_44, gs_11, gp_19, gp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * fd_22[k]
                  + f_2 * gd_s_42[k];

        t_43[k] = -f_1 * gs_s_11[k]
                  + f_2 * gd_s_43[k]
                  + f_3 * gs_11[k]
                  + pb_x[k] * gp_19[k];

        t_44[k] = f_2 * gd_s_44[k]
                  + pb_x[k] * gp_21[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pb_y, pb_z, fp_11, gs_s_11, gd_s_45, gd_s_46, \
                         gd_s_47, gs_11, gp_20, gp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -f_1 * gs_s_11[k]
                  + f_2 * gd_s_45[k]
                  + f_3 * gs_11[k]
                  + pb_y[k] * gp_20[k];

        t_46[k] = f_2 * gd_s_46[k]
                  + pb_y[k] * gp_21[k];

        t_47[k] = f_0 * fp_11[k]
                  - f_1 * gs_s_11[k]
                  + f_2 * gd_s_47[k]
                  + f_3 * gs_11[k]
                  + pb_z[k] * gp_21[k];
    }
}

auto
compute_prim_gd_kinetic_energy_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dd_s, const size_t dd,
                                 const size_t fp, const size_t fd, const size_t gs_s,
                                 const size_t gd_s, const size_t gs, const size_t gp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 * beta / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_6 = buffer.data(dd_s + 6);
    const auto *dd_s_9 = buffer.data(dd_s + 9);
    const auto *dd_s_12 = buffer.data(dd_s + 12);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_12 = buffer.data(dd + 12);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_22 = buffer.data(fd + 22);

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_2 = buffer.data(gs_s + 2);
    const auto *gs_s_3 = buffer.data(gs_s + 3);
    const auto *gs_s_6 = buffer.data(gs_s + 6);
    const auto *gs_s_8 = buffer.data(gs_s + 8);
    const auto *gs_s_10 = buffer.data(gs_s + 10);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_15 = buffer.data(gd_s + 15);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_22 = buffer.data(gd_s + 22);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_24 = buffer.data(gd_s + 24);
    const auto *gd_s_25 = buffer.data(gd_s + 25);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_27 = buffer.data(gd_s + 27);
    const auto *gd_s_28 = buffer.data(gd_s + 28);
    const auto *gd_s_29 = buffer.data(gd_s + 29);
    const auto *gd_s_30 = buffer.data(gd_s + 30);
    const auto *gd_s_31 = buffer.data(gd_s + 31);
    const auto *gd_s_32 = buffer.data(gd_s + 32);
    const auto *gd_s_33 = buffer.data(gd_s + 33);
    const auto *gd_s_34 = buffer.data(gd_s + 34);
    const auto *gd_s_35 = buffer.data(gd_s + 35);
    const auto *gd_s_36 = buffer.data(gd_s + 36);
    const auto *gd_s_37 = buffer.data(gd_s + 37);
    const auto *gd_s_38 = buffer.data(gd_s + 38);
    const auto *gd_s_39 = buffer.data(gd_s + 39);
    const auto *gd_s_40 = buffer.data(gd_s + 40);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_10 = buffer.data(gs + 10);

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
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fp_0, gs_s_0, gd_s_0, gd_s_1, \
                         gd_s_2, gs_0, gp_0, gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 - f_1 * gs_s_0[k]
                 + f_2 * gd_s_0[k]
                 + f_3 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_1[k]
                 + f_3 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_2[k]
                 + f_3 * gs_0[k]
                 + pb_z[k] * gp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pa_x, pa_y, pa_z, dd_s_2, dd_2, fd_0, fd_2, fd_4, \
                         gd_s_3, gd_s_4, gd_s_5, gd_s_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * fd_0[k]
                 + f_2 * gd_s_3[k];

        t_4[k] = -f_4 * dd_s_2[k]
                 + f_5 * dd_2[k]
                 + pa_x[k] * fd_4[k]
                 + f_2 * gd_s_4[k];

        t_5[k] = pa_y[k] * fd_2[k]
                 + f_2 * gd_s_5[k];

        t_6[k] = pa_z[k] * fd_0[k]
                 + f_2 * gd_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_y, pb_y, dd_s_0, dd_s_4, dd_0, dd_4, fd_3, \
                         fd_7, gd_s_7, gd_s_8, gd_s_9, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * gd_s_7[k]
                 + pb_y[k] * gp_3[k];

        t_8[k] = -f_4 * dd_s_4[k]
                 + f_5 * dd_4[k]
                 + pa_x[k] * fd_7[k]
                 + f_2 * gd_s_8[k];

        t_9[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_y[k] * fd_3[k]
                 + f_2 * gd_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_z, pb_z, dd_s_6, dd_6, fd_4, fd_9, gs_s_2, \
                         gd_s_10, gd_s_11, gd_s_12, gs_2, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_6 * dd_s_6[k]
                  + f_3 * dd_6[k]
                  + pa_x[k] * fd_9[k]
                  + f_2 * gd_s_10[k];

        t_11[k] = -f_1 * gs_s_2[k]
                  + f_2 * gd_s_11[k]
                  + f_3 * gs_2[k]
                  + pb_z[k] * gp_4[k];

        t_12[k] = pa_z[k] * fd_4[k]
                  + f_2 * gd_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pa_z, pb_y, dd_s_0, dd_0, fd_6, fd_7, gs_s_3, \
                         gd_s_13, gd_s_14, gd_s_15, gs_3, gp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * fd_7[k]
                  + f_2 * gd_s_13[k];

        t_14[k] = -f_6 * dd_s_0[k]
                  + f_3 * dd_0[k]
                  + pa_z[k] * fd_6[k]
                  + f_2 * gd_s_14[k];

        t_15[k] = -f_1 * gs_s_3[k]
                  + f_2 * gd_s_15[k]
                  + f_3 * gs_3[k]
                  + pb_y[k] * gp_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pb_y, dd_s_12, dd_12, fp_3, fd_11, fd_12, \
                         gd_s_16, gd_s_17, gd_s_18, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * gd_s_16[k]
                  + pb_y[k] * gp_6[k];

        t_17[k] = -f_6 * dd_s_12[k]
                  + f_3 * dd_12[k]
                  + pa_x[k] * fd_11[k]
                  + f_2 * gd_s_17[k];

        t_18[k] = f_5 * fp_3[k]
                  + pa_x[k] * fd_12[k]
                  + f_2 * gd_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pa_z, fd_8, fd_13, fd_16, fd_17, \
                         gd_s_19, gd_s_20, gd_s_21, gd_s_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_x[k] * fd_13[k]
                  + f_2 * gd_s_19[k];

        t_20[k] = pa_z[k] * fd_8[k]
                  + f_2 * gd_s_20[k];

        t_21[k] = pa_x[k] * fd_16[k]
                  + f_2 * gd_s_21[k];

        t_22[k] = pa_x[k] * fd_17[k]
                  + f_2 * gd_s_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pb_x, fp_6, fd_19, fd_22, gs_s_6, gd_s_23, \
                         gd_s_24, gd_s_25, gs_6, gp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * fp_6[k]
                  + pa_x[k] * fd_19[k]
                  + f_2 * gd_s_23[k];

        t_24[k] = pa_x[k] * fd_22[k]
                  + f_2 * gd_s_24[k];

        t_25[k] = -f_1 * gs_s_6[k]
                  + f_2 * gd_s_25[k]
                  + f_3 * gs_6[k]
                  + pb_x[k] * gp_7[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pb_x, pb_y, pb_z, fp_4, gs_s_6, gd_s_26, gd_s_27, \
                         gd_s_28, gs_6, gp_8, gp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * gd_s_26[k]
                  + pb_x[k] * gp_8[k];

        t_27[k] = f_0 * fp_4[k]
                  - f_1 * gs_s_6[k]
                  + f_2 * gd_s_27[k]
                  + f_3 * gs_6[k]
                  + pb_y[k] * gp_8[k];

        t_28[k] = -f_1 * gs_s_6[k]
                  + f_2 * gd_s_28[k]
                  + f_3 * gs_6[k]
                  + pb_z[k] * gp_9[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pa_z, pb_x, dd_s_9, dd_9, fd_13, fd_16, \
                         gs_s_8, gd_s_29, gd_s_30, gd_s_31, gs_8, \
                         gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_z[k] * fd_13[k]
                  + f_2 * gd_s_29[k];

        t_30[k] = -f_4 * dd_s_9[k]
                  + f_5 * dd_9[k]
                  + pa_y[k] * fd_16[k]
                  + f_2 * gd_s_30[k];

        t_31[k] = -f_1 * gs_s_8[k]
                  + f_2 * gd_s_31[k]
                  + f_3 * gs_8[k]
                  + pb_x[k] * gp_10[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pa_z, dd_s_6, dd_s_12, dd_6, dd_12, fp_7, \
                         fd_15, fd_18, fd_20, gd_s_32, gd_s_33, \
                         gd_s_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -f_6 * dd_s_6[k]
                  + f_3 * dd_6[k]
                  + pa_z[k] * fd_15[k]
                  + f_2 * gd_s_32[k];

        t_33[k] = -f_6 * dd_s_12[k]
                  + f_3 * dd_12[k]
                  + pa_y[k] * fd_18[k]
                  + f_2 * gd_s_33[k];

        t_34[k] = f_5 * fp_7[k]
                  + pa_y[k] * fd_20[k]
                  + f_2 * gd_s_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pb_x, fd_22, gs_s_10, gd_s_35, gd_s_36, \
                         gd_s_37, gs_10, gp_13, gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_y[k] * fd_22[k]
                  + f_2 * gd_s_35[k];

        t_36[k] = -f_1 * gs_s_10[k]
                  + f_2 * gd_s_36[k]
                  + f_3 * gs_10[k]
                  + pb_x[k] * gp_13[k];

        t_37[k] = f_2 * gd_s_37[k]
                  + pb_x[k] * gp_15[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pb_y, pb_z, fp_8, gs_s_10, gd_s_38, gd_s_39, \
                         gd_s_40, gs_10, gp_14, gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -f_1 * gs_s_10[k]
                  + f_2 * gd_s_38[k]
                  + f_3 * gs_10[k]
                  + pb_y[k] * gp_14[k];

        t_39[k] = f_2 * gd_s_39[k]
                  + pb_y[k] * gp_15[k];

        t_40[k] = f_0 * fp_8[k]
                  - f_1 * gs_s_10[k]
                  + f_2 * gd_s_40[k]
                  + f_3 * gs_10[k]
                  + pb_z[k] * gp_15[k];
    }
}

auto
compute_prim_gd_kinetic_energy_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dd_s, const size_t dd,
                                 const size_t fp, const size_t fd, const size_t gs_s,
                                 const size_t gd_s, const size_t gs, const size_t gp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 * beta / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;

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

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_3 = buffer.data(dd_s + 3);
    const auto *dd_s_5 = buffer.data(dd_s + 5);
    const auto *dd_s_7 = buffer.data(dd_s + 7);
    const auto *dd_s_10 = buffer.data(dd_s + 10);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_10 = buffer.data(dd + 10);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_19 = buffer.data(fd + 19);

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_2 = buffer.data(gs_s + 2);
    const auto *gs_s_3 = buffer.data(gs_s + 3);
    const auto *gs_s_6 = buffer.data(gs_s + 6);
    const auto *gs_s_8 = buffer.data(gs_s + 8);
    const auto *gs_s_10 = buffer.data(gs_s + 10);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_22 = buffer.data(gd_s + 22);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_24 = buffer.data(gd_s + 24);
    const auto *gd_s_25 = buffer.data(gd_s + 25);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_27 = buffer.data(gd_s + 27);
    const auto *gd_s_28 = buffer.data(gd_s + 28);
    const auto *gd_s_29 = buffer.data(gd_s + 29);
    const auto *gd_s_30 = buffer.data(gd_s + 30);
    const auto *gd_s_31 = buffer.data(gd_s + 31);
    const auto *gd_s_32 = buffer.data(gd_s + 32);
    const auto *gd_s_33 = buffer.data(gd_s + 33);
    const auto *gd_s_34 = buffer.data(gd_s + 34);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_10 = buffer.data(gs + 10);

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
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fp_0, gs_s_0, gd_s_0, gd_s_1, \
                         gd_s_2, gs_0, gp_0, gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 - f_1 * gs_s_0[k]
                 + f_2 * gd_s_0[k]
                 + f_3 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_1[k]
                 + f_3 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_2[k]
                 + f_3 * gs_0[k]
                 + pb_z[k] * gp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pa_z, dd_s_2, dd_2, fd_0, fd_4, gd_s_3, \
                         gd_s_4, gd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * fd_0[k]
                 + f_2 * gd_s_3[k];

        t_4[k] = -f_4 * dd_s_2[k]
                 + f_5 * dd_2[k]
                 + pa_x[k] * fd_4[k]
                 + f_2 * gd_s_4[k];

        t_5[k] = pa_z[k] * fd_0[k]
                 + f_2 * gd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, pb_y, dd_s_0, dd_s_3, dd_0, dd_3, fd_3, \
                         fd_6, gd_s_6, gd_s_7, gd_s_8, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * gd_s_6[k]
                 + pb_y[k] * gp_3[k];

        t_7[k] = -f_4 * dd_s_3[k]
                 + f_5 * dd_3[k]
                 + pa_x[k] * fd_6[k]
                 + f_2 * gd_s_7[k];

        t_8[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_y[k] * fd_3[k]
                 + f_2 * gd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, pa_x, pb_z, dd_s_5, dd_5, fd_7, gs_s_2, gd_s_9, gd_s_10, \
                         gs_2, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_6 * dd_s_5[k]
                 + f_3 * dd_5[k]
                 + pa_x[k] * fd_7[k]
                 + f_2 * gd_s_9[k];

        t_10[k] = -f_1 * gs_s_2[k]
                  + f_2 * gd_s_10[k]
                  + f_3 * gs_2[k]
                  + pb_z[k] * gp_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_z, pb_y, dd_s_0, dd_0, fd_5, gs_s_3, gd_s_11, \
                         gd_s_12, gd_s_13, gs_3, gp_5, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = -f_6 * dd_s_0[k]
                  + f_3 * dd_0[k]
                  + pa_z[k] * fd_5[k]
                  + f_2 * gd_s_11[k];

        t_12[k] = -f_1 * gs_s_3[k]
                  + f_2 * gd_s_12[k]
                  + f_3 * gs_3[k]
                  + pb_y[k] * gp_5[k];

        t_13[k] = f_2 * gd_s_13[k]
                  + pb_y[k] * gp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, dd_s_10, dd_10, fd_8, fd_10, fd_19, gd_s_14, \
                         gd_s_16, gd_s_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -f_6 * dd_s_10[k]
                  + f_3 * dd_10[k]
                  + pa_x[k] * fd_8[k]
                  + f_2 * gd_s_14[k];

        t_15[k] = pa_x[k] * fd_10[k]
                  + f_2 * gd_s_16[k];

        t_16[k] = pa_x[k] * fd_19[k]
                  + f_2 * gd_s_18[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_x, pb_y, fp_4, gs_s_6, gd_s_19, gd_s_20, \
                         gd_s_21, gs_6, gp_7, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_1 * gs_s_6[k]
                  + f_2 * gd_s_19[k]
                  + f_3 * gs_6[k]
                  + pb_x[k] * gp_7[k];

        t_18[k] = f_2 * gd_s_20[k]
                  + pb_x[k] * gp_8[k];

        t_19[k] = f_0 * fp_4[k]
                  - f_1 * gs_s_6[k]
                  + f_2 * gd_s_21[k]
                  + f_3 * gs_6[k]
                  + pb_y[k] * gp_8[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pa_z, pb_z, dd_s_7, dd_7, fd_10, fd_13, \
                         gs_s_6, gd_s_22, gd_s_23, gd_s_24, gs_6, \
                         gp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_1 * gs_s_6[k]
                  + f_2 * gd_s_22[k]
                  + f_3 * gs_6[k]
                  + pb_z[k] * gp_9[k];

        t_21[k] = pa_z[k] * fd_10[k]
                  + f_2 * gd_s_23[k];

        t_22[k] = -f_4 * dd_s_7[k]
                  + f_5 * dd_7[k]
                  + pa_y[k] * fd_13[k]
                  + f_2 * gd_s_24[k];
    }

#pragma omp simd aligned(t_23, t_24, pa_z, pb_x, dd_s_5, dd_5, fd_12, gs_s_8, gd_s_25, \
                         gd_s_26, gs_8, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -f_1 * gs_s_8[k]
                  + f_2 * gd_s_25[k]
                  + f_3 * gs_8[k]
                  + pb_x[k] * gp_10[k];

        t_24[k] = -f_6 * dd_s_5[k]
                  + f_3 * dd_5[k]
                  + pa_z[k] * fd_12[k]
                  + f_2 * gd_s_26[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, dd_s_10, dd_10, fp_7, fd_15, fd_17, fd_19, \
                         gd_s_27, gd_s_28, gd_s_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -f_6 * dd_s_10[k]
                  + f_3 * dd_10[k]
                  + pa_y[k] * fd_15[k]
                  + f_2 * gd_s_27[k];

        t_26[k] = f_5 * fp_7[k]
                  + pa_y[k] * fd_17[k]
                  + f_2 * gd_s_28[k];

        t_27[k] = pa_y[k] * fd_19[k]
                  + f_2 * gd_s_29[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_x, pb_y, gs_s_10, gd_s_30, gd_s_31, \
                         gd_s_32, gd_s_33, gs_10, gp_13, gp_14, gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -f_1 * gs_s_10[k]
                  + f_2 * gd_s_30[k]
                  + f_3 * gs_10[k]
                  + pb_x[k] * gp_13[k];

        t_29[k] = f_2 * gd_s_31[k]
                  + pb_x[k] * gp_15[k];

        t_30[k] = -f_1 * gs_s_10[k]
                  + f_2 * gd_s_32[k]
                  + f_3 * gs_10[k]
                  + pb_y[k] * gp_14[k];

        t_31[k] = f_2 * gd_s_33[k]
                  + pb_y[k] * gp_15[k];
    }

#pragma omp simd aligned(t_32, pb_z, fp_8, gs_s_10, gd_s_34, gs_10, \
                         gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * fp_8[k]
                  - f_1 * gs_s_10[k]
                  + f_2 * gd_s_34[k]
                  + f_3 * gs_10[k]
                  + pb_z[k] * gp_15[k];
    }
}

auto
compute_prim_gd_kinetic_energy_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dd_s, const size_t dd,
                                 const size_t fp, const size_t fd, const size_t gs_s,
                                 const size_t gd_s, const size_t gs, const size_t gp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 * beta / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_3 = buffer.data(dd_s + 3);
    const auto *dd_s_5 = buffer.data(dd_s + 5);
    const auto *dd_s_7 = buffer.data(dd_s + 7);
    const auto *dd_s_10 = buffer.data(dd_s + 10);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_10 = buffer.data(dd + 10);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_8 = buffer.data(fp + 8);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_18 = buffer.data(fd + 18);

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_2 = buffer.data(gs_s + 2);
    const auto *gs_s_3 = buffer.data(gs_s + 3);
    const auto *gs_s_6 = buffer.data(gs_s + 6);
    const auto *gs_s_8 = buffer.data(gs_s + 8);
    const auto *gs_s_10 = buffer.data(gs_s + 10);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_22 = buffer.data(gd_s + 22);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_24 = buffer.data(gd_s + 24);
    const auto *gd_s_25 = buffer.data(gd_s + 25);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_27 = buffer.data(gd_s + 27);
    const auto *gd_s_29 = buffer.data(gd_s + 29);
    const auto *gd_s_30 = buffer.data(gd_s + 30);
    const auto *gd_s_31 = buffer.data(gd_s + 31);
    const auto *gd_s_32 = buffer.data(gd_s + 32);
    const auto *gd_s_33 = buffer.data(gd_s + 33);
    const auto *gd_s_34 = buffer.data(gd_s + 34);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_10 = buffer.data(gs + 10);

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
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fp_0, gs_s_0, gd_s_0, gd_s_1, \
                         gd_s_2, gs_0, gp_0, gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 - f_1 * gs_s_0[k]
                 + f_2 * gd_s_0[k]
                 + f_3 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_1[k]
                 + f_3 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_2[k]
                 + f_3 * gs_0[k]
                 + pb_z[k] * gp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pa_z, dd_s_2, dd_2, fd_0, fd_4, gd_s_3, \
                         gd_s_4, gd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * fd_0[k]
                 + f_2 * gd_s_3[k];

        t_4[k] = -f_4 * dd_s_2[k]
                 + f_5 * dd_2[k]
                 + pa_x[k] * fd_4[k]
                 + f_2 * gd_s_4[k];

        t_5[k] = pa_z[k] * fd_0[k]
                 + f_2 * gd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, pb_y, dd_s_0, dd_s_3, dd_0, dd_3, fd_3, \
                         fd_6, gd_s_6, gd_s_7, gd_s_8, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * gd_s_6[k]
                 + pb_y[k] * gp_3[k];

        t_7[k] = -f_4 * dd_s_3[k]
                 + f_5 * dd_3[k]
                 + pa_x[k] * fd_6[k]
                 + f_2 * gd_s_7[k];

        t_8[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_y[k] * fd_3[k]
                 + f_2 * gd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, pa_x, pb_z, dd_s_5, dd_5, fd_7, gs_s_2, gd_s_9, gd_s_10, \
                         gs_2, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_6 * dd_s_5[k]
                 + f_3 * dd_5[k]
                 + pa_x[k] * fd_7[k]
                 + f_2 * gd_s_9[k];

        t_10[k] = -f_1 * gs_s_2[k]
                  + f_2 * gd_s_10[k]
                  + f_3 * gs_2[k]
                  + pb_z[k] * gp_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_z, pb_y, dd_s_0, dd_0, fd_5, gs_s_3, gd_s_11, \
                         gd_s_12, gd_s_13, gs_3, gp_5, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = -f_6 * dd_s_0[k]
                  + f_3 * dd_0[k]
                  + pa_z[k] * fd_5[k]
                  + f_2 * gd_s_11[k];

        t_12[k] = -f_1 * gs_s_3[k]
                  + f_2 * gd_s_12[k]
                  + f_3 * gs_3[k]
                  + pb_y[k] * gp_5[k];

        t_13[k] = f_2 * gd_s_13[k]
                  + pb_y[k] * gp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, dd_s_10, dd_10, fd_8, fd_10, fd_18, gd_s_14, \
                         gd_s_16, gd_s_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -f_6 * dd_s_10[k]
                  + f_3 * dd_10[k]
                  + pa_x[k] * fd_8[k]
                  + f_2 * gd_s_14[k];

        t_15[k] = pa_x[k] * fd_10[k]
                  + f_2 * gd_s_16[k];

        t_16[k] = pa_x[k] * fd_18[k]
                  + f_2 * gd_s_18[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_x, pb_y, fp_4, gs_s_6, gd_s_19, gd_s_20, \
                         gd_s_21, gs_6, gp_7, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_1 * gs_s_6[k]
                  + f_2 * gd_s_19[k]
                  + f_3 * gs_6[k]
                  + pb_x[k] * gp_7[k];

        t_18[k] = f_2 * gd_s_20[k]
                  + pb_x[k] * gp_8[k];

        t_19[k] = f_0 * fp_4[k]
                  - f_1 * gs_s_6[k]
                  + f_2 * gd_s_21[k]
                  + f_3 * gs_6[k]
                  + pb_y[k] * gp_8[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pa_z, pb_z, dd_s_7, dd_7, fd_10, fd_13, \
                         gs_s_6, gd_s_22, gd_s_23, gd_s_24, gs_6, \
                         gp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_1 * gs_s_6[k]
                  + f_2 * gd_s_22[k]
                  + f_3 * gs_6[k]
                  + pb_z[k] * gp_9[k];

        t_21[k] = pa_z[k] * fd_10[k]
                  + f_2 * gd_s_23[k];

        t_22[k] = -f_4 * dd_s_7[k]
                  + f_5 * dd_7[k]
                  + pa_y[k] * fd_13[k]
                  + f_2 * gd_s_24[k];
    }

#pragma omp simd aligned(t_23, t_24, pa_z, pb_x, dd_s_5, dd_5, fd_12, gs_s_8, gd_s_25, \
                         gd_s_26, gs_8, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -f_1 * gs_s_8[k]
                  + f_2 * gd_s_25[k]
                  + f_3 * gs_8[k]
                  + pb_x[k] * gp_10[k];

        t_24[k] = -f_6 * dd_s_5[k]
                  + f_3 * dd_5[k]
                  + pa_z[k] * fd_12[k]
                  + f_2 * gd_s_26[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pb_x, dd_s_10, dd_10, fd_14, fd_18, gs_s_10, \
                         gd_s_27, gd_s_29, gd_s_30, gs_10, gp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -f_6 * dd_s_10[k]
                  + f_3 * dd_10[k]
                  + pa_y[k] * fd_14[k]
                  + f_2 * gd_s_27[k];

        t_26[k] = pa_y[k] * fd_18[k]
                  + f_2 * gd_s_29[k];

        t_27[k] = -f_1 * gs_s_10[k]
                  + f_2 * gd_s_30[k]
                  + f_3 * gs_10[k]
                  + pb_x[k] * gp_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_x, pb_y, pb_z, fp_8, gs_s_10, gd_s_31, \
                         gd_s_32, gd_s_33, gd_s_34, gs_10, gp_14, \
                         gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_2 * gd_s_31[k]
                  + pb_x[k] * gp_15[k];

        t_29[k] = -f_1 * gs_s_10[k]
                  + f_2 * gd_s_32[k]
                  + f_3 * gs_10[k]
                  + pb_y[k] * gp_14[k];

        t_30[k] = f_2 * gd_s_33[k]
                  + pb_y[k] * gp_15[k];

        t_31[k] = f_0 * fp_8[k]
                  - f_1 * gs_s_10[k]
                  + f_2 * gd_s_34[k]
                  + f_3 * gs_10[k]
                  + pb_z[k] * gp_15[k];
    }
}

auto
compute_prim_gd_kinetic_energy_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t dd_s, const size_t dd,
                                  const size_t fp, const size_t fd, const size_t gs_s,
                                  const size_t gd_s, const size_t gs, const size_t gp,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 * beta / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;

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

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_3 = buffer.data(dd_s + 3);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_6 = buffer.data(dd_s + 6);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);

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

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_7 = buffer.data(gs_s + 7);
    const auto *gs_s_11 = buffer.data(gs_s + 11);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_15 = buffer.data(gd_s + 15);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_22 = buffer.data(gd_s + 22);
    const auto *gd_s_23 = buffer.data(gd_s + 23);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_11 = buffer.data(gs + 11);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);

#pragma omp simd aligned(t_0, t_1, pa_y, pb_x, fp_0, fd_0, gs_s_0, gd_s_0, gd_s_1, gs_0, \
                         gp_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 - f_1 * gs_s_0[k]
                 + f_2 * gd_s_0[k]
                 + f_3 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = pa_y[k] * fd_0[k]
                 + f_2 * gd_s_1[k];
    }

#pragma omp simd aligned(t_2, t_3, t_4, pa_x, pa_z, dd_s_1, dd_s_2, dd_1, dd_2, fd_0, fd_2, \
                         fd_4, gd_s_2, gd_s_3, gd_s_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = -f_4 * dd_s_1[k]
                 + f_5 * dd_1[k]
                 + pa_x[k] * fd_2[k]
                 + f_2 * gd_s_2[k];

        t_3[k] = pa_z[k] * fd_0[k]
                 + f_2 * gd_s_3[k];

        t_4[k] = -f_4 * dd_s_2[k]
                 + f_5 * dd_2[k]
                 + pa_x[k] * fd_4[k]
                 + f_2 * gd_s_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_x, pa_y, pa_z, dd_s_0, dd_s_3, dd_0, dd_3, fd_1, \
                         fd_3, fd_5, gd_s_5, gd_s_6, gd_s_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_y[k] * fd_1[k]
                 + f_2 * gd_s_5[k];

        t_6[k] = -f_6 * dd_s_3[k]
                 + f_3 * dd_3[k]
                 + pa_x[k] * fd_5[k]
                 + f_2 * gd_s_6[k];

        t_7[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_z[k] * fd_3[k]
                 + f_2 * gd_s_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, dd_s_6, dd_6, fd_6, fd_7, fd_9, fd_10, \
                         gd_s_8, gd_s_9, gd_s_10, gd_s_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = -f_6 * dd_s_6[k]
                 + f_3 * dd_6[k]
                 + pa_x[k] * fd_6[k]
                 + f_2 * gd_s_8[k];

        t_9[k] = pa_x[k] * fd_7[k]
                 + f_2 * gd_s_9[k];

        t_10[k] = pa_x[k] * fd_9[k]
                  + f_2 * gd_s_10[k];

        t_11[k] = pa_x[k] * fd_10[k]
                  + f_2 * gd_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_x, pb_y, fp_1, fd_13, gs_s_7, gd_s_12, \
                         gd_s_13, gd_s_14, gs_7, gp_7, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_x[k] * fd_13[k]
                  + f_2 * gd_s_12[k];

        t_13[k] = -f_1 * gs_s_7[k]
                  + f_2 * gd_s_13[k]
                  + f_3 * gs_7[k]
                  + pb_x[k] * gp_7[k];

        t_14[k] = f_0 * fp_1[k]
                  - f_1 * gs_s_7[k]
                  + f_2 * gd_s_14[k]
                  + f_3 * gs_7[k]
                  + pb_y[k] * gp_8[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_y, pa_z, dd_s_3, dd_s_4, dd_3, dd_4, fd_7, fd_8, \
                         fd_9, gd_s_15, gd_s_16, gd_s_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_z[k] * fd_7[k]
                  + f_2 * gd_s_15[k];

        t_16[k] = -f_4 * dd_s_4[k]
                  + f_5 * dd_4[k]
                  + pa_y[k] * fd_9[k]
                  + f_2 * gd_s_16[k];

        t_17[k] = -f_6 * dd_s_3[k]
                  + f_3 * dd_3[k]
                  + pa_z[k] * fd_8[k]
                  + f_2 * gd_s_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, dd_s_6, dd_6, fp_2, fd_11, fd_12, fd_13, \
                         gd_s_18, gd_s_19, gd_s_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -f_6 * dd_s_6[k]
                  + f_3 * dd_6[k]
                  + pa_y[k] * fd_11[k]
                  + f_2 * gd_s_18[k];

        t_19[k] = f_5 * fp_2[k]
                  + pa_y[k] * fd_12[k]
                  + f_2 * gd_s_19[k];

        t_20[k] = pa_y[k] * fd_13[k]
                  + f_2 * gd_s_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_x, pb_y, pb_z, fp_3, gs_s_11, gd_s_21, gd_s_22, \
                         gd_s_23, gs_11, gp_14, gp_15, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_1 * gs_s_11[k]
                  + f_2 * gd_s_21[k]
                  + f_3 * gs_11[k]
                  + pb_x[k] * gp_14[k];

        t_22[k] = -f_1 * gs_s_11[k]
                  + f_2 * gd_s_22[k]
                  + f_3 * gs_11[k]
                  + pb_y[k] * gp_15[k];

        t_23[k] = f_0 * fp_3[k]
                  - f_1 * gs_s_11[k]
                  + f_2 * gd_s_23[k]
                  + f_3 * gs_11[k]
                  + pb_z[k] * gp_16[k];
    }
}

auto
compute_prim_gd_kinetic_energy_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t dd_s, const size_t dd,
                                  const size_t fp, const size_t fd, const size_t gs_s,
                                  const size_t gd_s, const size_t gs, const size_t gp,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 * beta / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 1.5 / p;

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

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_3 = buffer.data(dd_s + 3);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_6 = buffer.data(dd_s + 6);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);

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

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_7 = buffer.data(gs_s + 7);
    const auto *gs_s_9 = buffer.data(gs_s + 9);
    const auto *gs_s_11 = buffer.data(gs_s + 11);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_15 = buffer.data(gd_s + 15);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_22 = buffer.data(gd_s + 22);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_24 = buffer.data(gd_s + 24);
    const auto *gd_s_25 = buffer.data(gd_s + 25);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_27 = buffer.data(gd_s + 27);
    const auto *gd_s_28 = buffer.data(gd_s + 28);
    const auto *gd_s_29 = buffer.data(gd_s + 29);
    const auto *gd_s_30 = buffer.data(gd_s + 30);
    const auto *gd_s_31 = buffer.data(gd_s + 31);
    const auto *gd_s_32 = buffer.data(gd_s + 32);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_11 = buffer.data(gs + 11);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fp_0, gs_s_0, gd_s_0, gd_s_1, \
                         gd_s_2, gs_0, gp_0, gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 - f_1 * gs_s_0[k]
                 + f_2 * gd_s_0[k]
                 + f_3 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_1[k]
                 + f_3 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_2[k]
                 + f_3 * gs_0[k]
                 + pb_z[k] * gp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pa_z, dd_s_1, dd_1, fd_0, fd_2, gd_s_3, \
                         gd_s_4, gd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * fd_0[k]
                 + f_2 * gd_s_3[k];

        t_4[k] = -f_4 * dd_s_1[k]
                 + f_5 * dd_1[k]
                 + pa_x[k] * fd_2[k]
                 + f_2 * gd_s_4[k];

        t_5[k] = pa_z[k] * fd_0[k]
                 + f_2 * gd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, pa_x, pa_y, dd_s_0, dd_s_2, dd_0, dd_2, fd_1, fd_4, gd_s_6, \
                         gd_s_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * dd_s_2[k]
                 + f_5 * dd_2[k]
                 + pa_x[k] * fd_4[k]
                 + f_2 * gd_s_6[k];

        t_7[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_y[k] * fd_1[k]
                 + f_2 * gd_s_7[k];
    }

#pragma omp simd aligned(t_8, t_9, pa_x, pa_z, dd_s_0, dd_s_3, dd_0, dd_3, fd_3, fd_5, gd_s_8, \
                         gd_s_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = -f_6 * dd_s_3[k]
                 + f_3 * dd_3[k]
                 + pa_x[k] * fd_5[k]
                 + f_2 * gd_s_8[k];

        t_9[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_z[k] * fd_3[k]
                 + f_2 * gd_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, dd_s_6, dd_6, fp_2, fd_6, fd_7, fd_8, \
                         fd_10, gd_s_10, gd_s_11, gd_s_12, gd_s_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_6 * dd_s_6[k]
                  + f_3 * dd_6[k]
                  + pa_x[k] * fd_6[k]
                  + f_2 * gd_s_10[k];

        t_11[k] = f_5 * fp_2[k]
                  + pa_x[k] * fd_7[k]
                  + f_2 * gd_s_11[k];

        t_12[k] = pa_x[k] * fd_8[k]
                  + f_2 * gd_s_12[k];

        t_13[k] = pa_x[k] * fd_10[k]
                  + f_2 * gd_s_13[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, fp_6, fd_11, fd_13, fd_15, gd_s_14, gd_s_15, \
                         gd_s_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_x[k] * fd_11[k]
                  + f_2 * gd_s_14[k];

        t_15[k] = f_5 * fp_6[k]
                  + pa_x[k] * fd_13[k]
                  + f_2 * gd_s_15[k];

        t_16[k] = pa_x[k] * fd_15[k]
                  + f_2 * gd_s_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_x, pb_y, pb_z, fp_3, gs_s_7, gd_s_17, gd_s_18, \
                         gd_s_19, gs_7, gp_6, gp_7, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_1 * gs_s_7[k]
                  + f_2 * gd_s_17[k]
                  + f_3 * gs_7[k]
                  + pb_x[k] * gp_6[k];

        t_18[k] = f_0 * fp_3[k]
                  - f_1 * gs_s_7[k]
                  + f_2 * gd_s_18[k]
                  + f_3 * gs_7[k]
                  + pb_y[k] * gp_7[k];

        t_19[k] = -f_1 * gs_s_7[k]
                  + f_2 * gd_s_19[k]
                  + f_3 * gs_7[k]
                  + pb_z[k] * gp_8[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pa_z, pb_y, dd_s_4, dd_4, fp_4, fd_8, fd_10, \
                         gd_s_20, gd_s_21, gd_s_22, gp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * fd_8[k]
                  + f_2 * gd_s_20[k];

        t_21[k] = f_7 * fp_4[k]
                  + f_2 * gd_s_21[k]
                  + pb_y[k] * gp_9[k];

        t_22[k] = -f_4 * dd_s_4[k]
                  + f_5 * dd_4[k]
                  + pa_y[k] * fd_10[k]
                  + f_2 * gd_s_22[k];
    }

#pragma omp simd aligned(t_23, t_24, pa_z, pb_x, dd_s_3, dd_3, fd_9, gs_s_9, gd_s_23, gd_s_24, \
                         gs_9, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -f_1 * gs_s_9[k]
                  + f_2 * gd_s_23[k]
                  + f_3 * gs_9[k]
                  + pb_x[k] * gp_10[k];

        t_24[k] = -f_6 * dd_s_3[k]
                  + f_3 * dd_3[k]
                  + pa_z[k] * fd_9[k]
                  + f_2 * gd_s_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pb_y, dd_s_6, dd_6, fp_5, fp_7, fd_12, fd_14, \
                         gd_s_25, gd_s_26, gd_s_27, gp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_5 * fp_5[k]
                  + f_2 * gd_s_25[k]
                  + pb_y[k] * gp_11[k];

        t_26[k] = -f_6 * dd_s_6[k]
                  + f_3 * dd_6[k]
                  + pa_y[k] * fd_12[k]
                  + f_2 * gd_s_26[k];

        t_27[k] = f_5 * fp_7[k]
                  + pa_y[k] * fd_14[k]
                  + f_2 * gd_s_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pb_x, pb_y, fp_8, fd_15, gs_s_11, gd_s_28, \
                         gd_s_29, gd_s_30, gs_11, gp_12, gp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * fp_8[k]
                  + f_2 * gd_s_28[k]
                  + pb_y[k] * gp_12[k];

        t_29[k] = pa_y[k] * fd_15[k]
                  + f_2 * gd_s_29[k];

        t_30[k] = -f_1 * gs_s_11[k]
                  + f_2 * gd_s_30[k]
                  + f_3 * gs_11[k]
                  + pb_x[k] * gp_13[k];
    }

#pragma omp simd aligned(t_31, t_32, pb_y, pb_z, fp_8, gs_s_11, gd_s_31, gd_s_32, gs_11, \
                         gp_14, gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -f_1 * gs_s_11[k]
                  + f_2 * gd_s_31[k]
                  + f_3 * gs_11[k]
                  + pb_y[k] * gp_14[k];

        t_32[k] = f_0 * fp_8[k]
                  - f_1 * gs_s_11[k]
                  + f_2 * gd_s_32[k]
                  + f_3 * gs_11[k]
                  + pb_z[k] * gp_15[k];
    }
}

auto
compute_prim_gd_kinetic_energy_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t dd_s, const size_t dd,
                                  const size_t fp, const size_t fd, const size_t gs_s,
                                  const size_t gd_s, const size_t gs, const size_t gp,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 * beta / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;

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

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_5 = buffer.data(dd_s + 5);
    const auto *dd_s_8 = buffer.data(dd_s + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_2 = buffer.data(gs_s + 2);
    const auto *gs_s_6 = buffer.data(gs_s + 6);
    const auto *gs_s_8 = buffer.data(gs_s + 8);
    const auto *gs_s_10 = buffer.data(gs_s + 10);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_15 = buffer.data(gd_s + 15);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_22 = buffer.data(gd_s + 22);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_24 = buffer.data(gd_s + 24);
    const auto *gd_s_25 = buffer.data(gd_s + 25);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_27 = buffer.data(gd_s + 27);
    const auto *gd_s_28 = buffer.data(gd_s + 28);
    const auto *gd_s_29 = buffer.data(gd_s + 29);
    const auto *gd_s_30 = buffer.data(gd_s + 30);
    const auto *gd_s_31 = buffer.data(gd_s + 31);
    const auto *gd_s_32 = buffer.data(gd_s + 32);
    const auto *gd_s_33 = buffer.data(gd_s + 33);
    const auto *gd_s_34 = buffer.data(gd_s + 34);
    const auto *gd_s_35 = buffer.data(gd_s + 35);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_10 = buffer.data(gs + 10);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fp_0, gs_s_0, gd_s_0, gd_s_1, \
                         gd_s_2, gs_0, gp_0, gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 - f_1 * gs_s_0[k]
                 + f_2 * gd_s_0[k]
                 + f_3 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_1[k]
                 + f_3 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_2[k]
                 + f_3 * gs_0[k]
                 + pb_z[k] * gp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pa_x, pa_y, pa_z, dd_s_1, dd_1, fd_0, fd_1, fd_3, \
                         gd_s_3, gd_s_4, gd_s_5, gd_s_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * fd_0[k]
                 + f_2 * gd_s_3[k];

        t_4[k] = -f_4 * dd_s_1[k]
                 + f_5 * dd_1[k]
                 + pa_x[k] * fd_3[k]
                 + f_2 * gd_s_4[k];

        t_5[k] = pa_y[k] * fd_1[k]
                 + f_2 * gd_s_5[k];

        t_6[k] = pa_z[k] * fd_0[k]
                 + f_2 * gd_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pa_x, pa_y, dd_s_0, dd_s_2, dd_0, dd_2, fd_2, fd_5, gd_s_7, \
                         gd_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -f_4 * dd_s_2[k]
                 + f_5 * dd_2[k]
                 + pa_x[k] * fd_5[k]
                 + f_2 * gd_s_7[k];

        t_8[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_y[k] * fd_2[k]
                 + f_2 * gd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_z, pb_z, dd_s_4, dd_4, fd_3, fd_7, gs_s_2, \
                         gd_s_9, gd_s_10, gd_s_11, gs_2, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_6 * dd_s_4[k]
                 + f_3 * dd_4[k]
                 + pa_x[k] * fd_7[k]
                 + f_2 * gd_s_9[k];

        t_10[k] = -f_1 * gs_s_2[k]
                  + f_2 * gd_s_10[k]
                  + f_3 * gs_2[k]
                  + pb_z[k] * gp_3[k];

        t_11[k] = pa_z[k] * fd_3[k]
                  + f_2 * gd_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_y, pa_z, dd_s_0, dd_s_8, dd_0, dd_8, fd_4, \
                         fd_5, fd_9, gd_s_12, gd_s_13, gd_s_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_y[k] * fd_5[k]
                  + f_2 * gd_s_12[k];

        t_13[k] = -f_6 * dd_s_0[k]
                  + f_3 * dd_0[k]
                  + pa_z[k] * fd_4[k]
                  + f_2 * gd_s_13[k];

        t_14[k] = -f_6 * dd_s_8[k]
                  + f_3 * dd_8[k]
                  + pa_x[k] * fd_9[k]
                  + f_2 * gd_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pa_z, fp_2, fd_6, fd_10, fd_11, fd_14, \
                         gd_s_15, gd_s_16, gd_s_17, gd_s_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * fp_2[k]
                  + pa_x[k] * fd_10[k]
                  + f_2 * gd_s_15[k];

        t_16[k] = pa_x[k] * fd_11[k]
                  + f_2 * gd_s_16[k];

        t_17[k] = pa_z[k] * fd_6[k]
                  + f_2 * gd_s_17[k];

        t_18[k] = pa_x[k] * fd_14[k]
                  + f_2 * gd_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, fp_5, fd_15, fd_17, fd_19, gd_s_19, gd_s_20, \
                         gd_s_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_x[k] * fd_15[k]
                  + f_2 * gd_s_19[k];

        t_20[k] = f_5 * fp_5[k]
                  + pa_x[k] * fd_17[k]
                  + f_2 * gd_s_20[k];

        t_21[k] = pa_x[k] * fd_19[k]
                  + f_2 * gd_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pb_x, pb_y, pb_z, fp_3, gs_s_6, gd_s_22, gd_s_23, \
                         gd_s_24, gs_6, gp_4, gp_5, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -f_1 * gs_s_6[k]
                  + f_2 * gd_s_22[k]
                  + f_3 * gs_6[k]
                  + pb_x[k] * gp_4[k];

        t_23[k] = f_0 * fp_3[k]
                  - f_1 * gs_s_6[k]
                  + f_2 * gd_s_23[k]
                  + f_3 * gs_6[k]
                  + pb_y[k] * gp_5[k];

        t_24[k] = -f_1 * gs_s_6[k]
                  + f_2 * gd_s_24[k]
                  + f_3 * gs_6[k]
                  + pb_z[k] * gp_6[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pb_x, dd_s_5, dd_5, fd_11, fd_14, \
                         gs_s_8, gd_s_25, gd_s_26, gd_s_27, gs_8, \
                         gp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_z[k] * fd_11[k]
                  + f_2 * gd_s_25[k];

        t_26[k] = -f_4 * dd_s_5[k]
                  + f_5 * dd_5[k]
                  + pa_y[k] * fd_14[k]
                  + f_2 * gd_s_26[k];

        t_27[k] = -f_1 * gs_s_8[k]
                  + f_2 * gd_s_27[k]
                  + f_3 * gs_8[k]
                  + pb_x[k] * gp_7[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pa_z, dd_s_4, dd_s_8, dd_4, dd_8, fp_6, \
                         fd_13, fd_16, fd_18, gd_s_28, gd_s_29, \
                         gd_s_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -f_6 * dd_s_4[k]
                  + f_3 * dd_4[k]
                  + pa_z[k] * fd_13[k]
                  + f_2 * gd_s_28[k];

        t_29[k] = -f_6 * dd_s_8[k]
                  + f_3 * dd_8[k]
                  + pa_y[k] * fd_16[k]
                  + f_2 * gd_s_29[k];

        t_30[k] = f_5 * fp_6[k]
                  + pa_y[k] * fd_18[k]
                  + f_2 * gd_s_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pb_x, pb_y, fd_19, gs_s_10, gd_s_31, gd_s_32, \
                         gd_s_33, gs_10, gp_8, gp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pa_y[k] * fd_19[k]
                  + f_2 * gd_s_31[k];

        t_32[k] = -f_1 * gs_s_10[k]
                  + f_2 * gd_s_32[k]
                  + f_3 * gs_10[k]
                  + pb_x[k] * gp_8[k];

        t_33[k] = -f_1 * gs_s_10[k]
                  + f_2 * gd_s_33[k]
                  + f_3 * gs_10[k]
                  + pb_y[k] * gp_9[k];
    }

#pragma omp simd aligned(t_34, t_35, pb_y, pb_z, fp_7, gs_s_10, gd_s_34, gd_s_35, gs_10, \
                         gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_2 * gd_s_34[k]
                  + pb_y[k] * gp_10[k];

        t_35[k] = f_0 * fp_7[k]
                  - f_1 * gs_s_10[k]
                  + f_2 * gd_s_35[k]
                  + f_3 * gs_10[k]
                  + pb_z[k] * gp_10[k];
    }
}

auto
compute_prim_gd_kinetic_energy_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t dd_s, const size_t dd,
                                  const size_t fp, const size_t fd, const size_t gs_s,
                                  const size_t gd_s, const size_t gs, const size_t gp,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 * beta / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_5 = buffer.data(dd_s + 5);
    const auto *dd_s_8 = buffer.data(dd_s + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_2 = buffer.data(gs_s + 2);
    const auto *gs_s_6 = buffer.data(gs_s + 6);
    const auto *gs_s_8 = buffer.data(gs_s + 8);
    const auto *gs_s_10 = buffer.data(gs_s + 10);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_15 = buffer.data(gd_s + 15);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_22 = buffer.data(gd_s + 22);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_24 = buffer.data(gd_s + 24);
    const auto *gd_s_25 = buffer.data(gd_s + 25);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_27 = buffer.data(gd_s + 27);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_10 = buffer.data(gs + 10);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fp_0, gs_s_0, gd_s_0, gd_s_1, \
                         gd_s_2, gs_0, gp_0, gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 - f_1 * gs_s_0[k]
                 + f_2 * gd_s_0[k]
                 + f_3 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_1[k]
                 + f_3 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_2[k]
                 + f_3 * gs_0[k]
                 + pb_z[k] * gp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pa_z, dd_s_1, dd_1, fd_0, fd_3, gd_s_3, \
                         gd_s_4, gd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * fd_0[k]
                 + f_2 * gd_s_3[k];

        t_4[k] = -f_4 * dd_s_1[k]
                 + f_5 * dd_1[k]
                 + pa_x[k] * fd_3[k]
                 + f_2 * gd_s_4[k];

        t_5[k] = pa_z[k] * fd_0[k]
                 + f_2 * gd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, pa_x, pa_y, dd_s_0, dd_s_2, dd_0, dd_2, fd_2, fd_5, gd_s_6, \
                         gd_s_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * dd_s_2[k]
                 + f_5 * dd_2[k]
                 + pa_x[k] * fd_5[k]
                 + f_2 * gd_s_6[k];

        t_7[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_y[k] * fd_2[k]
                 + f_2 * gd_s_7[k];
    }

#pragma omp simd aligned(t_8, t_9, pa_x, pb_z, dd_s_4, dd_4, fd_6, gs_s_2, gd_s_8, gd_s_9, \
                         gs_2, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = -f_6 * dd_s_4[k]
                 + f_3 * dd_4[k]
                 + pa_x[k] * fd_6[k]
                 + f_2 * gd_s_8[k];

        t_9[k] = -f_1 * gs_s_2[k]
                 + f_2 * gd_s_9[k]
                 + f_3 * gs_2[k]
                 + pb_z[k] * gp_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_z, dd_s_0, dd_s_8, dd_0, dd_8, fd_4, fd_7, \
                         fd_9, gd_s_10, gd_s_11, gd_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_6 * dd_s_0[k]
                  + f_3 * dd_0[k]
                  + pa_z[k] * fd_4[k]
                  + f_2 * gd_s_10[k];

        t_11[k] = -f_6 * dd_s_8[k]
                  + f_3 * dd_8[k]
                  + pa_x[k] * fd_7[k]
                  + f_2 * gd_s_11[k];

        t_12[k] = pa_x[k] * fd_9[k]
                  + f_2 * gd_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pb_x, pb_y, fp_3, fd_16, gs_s_6, gd_s_13, \
                         gd_s_14, gd_s_15, gs_6, gp_4, gp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_x[k] * fd_16[k]
                  + f_2 * gd_s_13[k];

        t_14[k] = -f_1 * gs_s_6[k]
                  + f_2 * gd_s_14[k]
                  + f_3 * gs_6[k]
                  + pb_x[k] * gp_4[k];

        t_15[k] = f_0 * fp_3[k]
                  - f_1 * gs_s_6[k]
                  + f_2 * gd_s_15[k]
                  + f_3 * gs_6[k]
                  + pb_y[k] * gp_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_y, pa_z, pb_z, dd_s_5, dd_5, fd_9, fd_12, \
                         gs_s_6, gd_s_16, gd_s_17, gd_s_18, gs_6, \
                         gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_1 * gs_s_6[k]
                  + f_2 * gd_s_16[k]
                  + f_3 * gs_6[k]
                  + pb_z[k] * gp_6[k];

        t_17[k] = pa_z[k] * fd_9[k]
                  + f_2 * gd_s_17[k];

        t_18[k] = -f_4 * dd_s_5[k]
                  + f_5 * dd_5[k]
                  + pa_y[k] * fd_12[k]
                  + f_2 * gd_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, pa_z, pb_x, dd_s_4, dd_4, fd_11, gs_s_8, gd_s_19, \
                         gd_s_20, gs_8, gp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_1 * gs_s_8[k]
                  + f_2 * gd_s_19[k]
                  + f_3 * gs_8[k]
                  + pb_x[k] * gp_7[k];

        t_20[k] = -f_6 * dd_s_4[k]
                  + f_3 * dd_4[k]
                  + pa_z[k] * fd_11[k]
                  + f_2 * gd_s_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_y, dd_s_8, dd_8, fp_6, fd_13, fd_15, fd_16, \
                         gd_s_21, gd_s_22, gd_s_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_6 * dd_s_8[k]
                  + f_3 * dd_8[k]
                  + pa_y[k] * fd_13[k]
                  + f_2 * gd_s_21[k];

        t_22[k] = f_5 * fp_6[k]
                  + pa_y[k] * fd_15[k]
                  + f_2 * gd_s_22[k];

        t_23[k] = pa_y[k] * fd_16[k]
                  + f_2 * gd_s_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_x, pb_y, gs_s_10, gd_s_24, gd_s_25, gd_s_26, \
                         gs_10, gp_8, gp_9, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_1 * gs_s_10[k]
                  + f_2 * gd_s_24[k]
                  + f_3 * gs_10[k]
                  + pb_x[k] * gp_8[k];

        t_25[k] = -f_1 * gs_s_10[k]
                  + f_2 * gd_s_25[k]
                  + f_3 * gs_10[k]
                  + pb_y[k] * gp_9[k];

        t_26[k] = f_2 * gd_s_26[k]
                  + pb_y[k] * gp_10[k];
    }

#pragma omp simd aligned(t_27, pb_z, fp_7, gs_s_10, gd_s_27, gs_10, \
                         gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * fp_7[k]
                  - f_1 * gs_s_10[k]
                  + f_2 * gd_s_27[k]
                  + f_3 * gs_10[k]
                  + pb_z[k] * gp_10[k];
    }
}

auto
compute_prim_gd_kinetic_energy_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t dd_s, const size_t dd,
                                  const size_t fp, const size_t fd, const size_t gs_s,
                                  const size_t gd_s, const size_t gs, const size_t gp,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 * beta / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;

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

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_5 = buffer.data(dd_s + 5);
    const auto *dd_s_8 = buffer.data(dd_s + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_7 = buffer.data(fp + 7);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_16 = buffer.data(fd + 16);

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_2 = buffer.data(gs_s + 2);
    const auto *gs_s_6 = buffer.data(gs_s + 6);
    const auto *gs_s_8 = buffer.data(gs_s + 8);
    const auto *gs_s_10 = buffer.data(gs_s + 10);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_15 = buffer.data(gd_s + 15);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_24 = buffer.data(gd_s + 24);
    const auto *gd_s_25 = buffer.data(gd_s + 25);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_27 = buffer.data(gd_s + 27);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_10 = buffer.data(gs + 10);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fp_0, gs_s_0, gd_s_0, gd_s_1, \
                         gd_s_2, gs_0, gp_0, gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 - f_1 * gs_s_0[k]
                 + f_2 * gd_s_0[k]
                 + f_3 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_1[k]
                 + f_3 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = -f_1 * gs_s_0[k]
                 + f_2 * gd_s_2[k]
                 + f_3 * gs_0[k]
                 + pb_z[k] * gp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pa_z, dd_s_1, dd_1, fd_0, fd_3, gd_s_3, \
                         gd_s_4, gd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * fd_0[k]
                 + f_2 * gd_s_3[k];

        t_4[k] = -f_4 * dd_s_1[k]
                 + f_5 * dd_1[k]
                 + pa_x[k] * fd_3[k]
                 + f_2 * gd_s_4[k];

        t_5[k] = pa_z[k] * fd_0[k]
                 + f_2 * gd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, pa_x, pa_y, dd_s_0, dd_s_2, dd_0, dd_2, fd_2, fd_5, gd_s_6, \
                         gd_s_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * dd_s_2[k]
                 + f_5 * dd_2[k]
                 + pa_x[k] * fd_5[k]
                 + f_2 * gd_s_6[k];

        t_7[k] = -f_6 * dd_s_0[k]
                 + f_3 * dd_0[k]
                 + pa_y[k] * fd_2[k]
                 + f_2 * gd_s_7[k];
    }

#pragma omp simd aligned(t_8, t_9, pa_x, pb_z, dd_s_4, dd_4, fd_6, gs_s_2, gd_s_8, gd_s_9, \
                         gs_2, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = -f_6 * dd_s_4[k]
                 + f_3 * dd_4[k]
                 + pa_x[k] * fd_6[k]
                 + f_2 * gd_s_8[k];

        t_9[k] = -f_1 * gs_s_2[k]
                 + f_2 * gd_s_9[k]
                 + f_3 * gs_2[k]
                 + pb_z[k] * gp_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_z, dd_s_0, dd_s_8, dd_0, dd_8, fd_4, fd_7, \
                         fd_9, gd_s_10, gd_s_11, gd_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_6 * dd_s_0[k]
                  + f_3 * dd_0[k]
                  + pa_z[k] * fd_4[k]
                  + f_2 * gd_s_10[k];

        t_11[k] = -f_6 * dd_s_8[k]
                  + f_3 * dd_8[k]
                  + pa_x[k] * fd_7[k]
                  + f_2 * gd_s_11[k];

        t_12[k] = pa_x[k] * fd_9[k]
                  + f_2 * gd_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pb_x, pb_y, fp_3, fd_16, gs_s_6, gd_s_13, \
                         gd_s_14, gd_s_15, gs_6, gp_4, gp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_x[k] * fd_16[k]
                  + f_2 * gd_s_13[k];

        t_14[k] = -f_1 * gs_s_6[k]
                  + f_2 * gd_s_14[k]
                  + f_3 * gs_6[k]
                  + pb_x[k] * gp_4[k];

        t_15[k] = f_0 * fp_3[k]
                  - f_1 * gs_s_6[k]
                  + f_2 * gd_s_15[k]
                  + f_3 * gs_6[k]
                  + pb_y[k] * gp_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_y, pa_z, pb_z, dd_s_5, dd_5, fd_9, fd_12, \
                         gs_s_6, gd_s_16, gd_s_17, gd_s_18, gs_6, \
                         gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_1 * gs_s_6[k]
                  + f_2 * gd_s_16[k]
                  + f_3 * gs_6[k]
                  + pb_z[k] * gp_6[k];

        t_17[k] = pa_z[k] * fd_9[k]
                  + f_2 * gd_s_17[k];

        t_18[k] = -f_4 * dd_s_5[k]
                  + f_5 * dd_5[k]
                  + pa_y[k] * fd_12[k]
                  + f_2 * gd_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, pa_z, pb_x, dd_s_4, dd_4, fd_11, gs_s_8, gd_s_19, \
                         gd_s_20, gs_8, gp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_1 * gs_s_8[k]
                  + f_2 * gd_s_19[k]
                  + f_3 * gs_8[k]
                  + pb_x[k] * gp_7[k];

        t_20[k] = -f_6 * dd_s_4[k]
                  + f_3 * dd_4[k]
                  + pa_z[k] * fd_11[k]
                  + f_2 * gd_s_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_y, pb_x, dd_s_8, dd_8, fd_13, fd_16, gs_s_10, \
                         gd_s_21, gd_s_23, gd_s_24, gs_10, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_6 * dd_s_8[k]
                  + f_3 * dd_8[k]
                  + pa_y[k] * fd_13[k]
                  + f_2 * gd_s_21[k];

        t_22[k] = pa_y[k] * fd_16[k]
                  + f_2 * gd_s_23[k];

        t_23[k] = -f_1 * gs_s_10[k]
                  + f_2 * gd_s_24[k]
                  + f_3 * gs_10[k]
                  + pb_x[k] * gp_8[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_y, pb_z, fp_7, gs_s_10, gd_s_25, gd_s_26, \
                         gd_s_27, gs_10, gp_9, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_1 * gs_s_10[k]
                  + f_2 * gd_s_25[k]
                  + f_3 * gs_10[k]
                  + pb_y[k] * gp_9[k];

        t_25[k] = f_2 * gd_s_26[k]
                  + pb_y[k] * gp_10[k];

        t_26[k] = f_0 * fp_7[k]
                  - f_1 * gs_s_10[k]
                  + f_2 * gd_s_27[k]
                  + f_3 * gs_10[k]
                  + pb_z[k] * gp_10[k];
    }
}

}  // namespace simdkin
