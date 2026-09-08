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


#include "SimdOverlapVrrRecGD.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_gd_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dd, const size_t fp, const size_t fd,
                          const size_t gs, const size_t gp, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = 1.0 / p;

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, fp_0, gs_0, gp_0, \
                         gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 + f_1 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = pb_y[k] * gp_0[k];

        t_2[k] = pb_z[k] * gp_0[k];

        t_3[k] = f_1 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_4[k] = pb_y[k] * gp_2[k];

        t_5[k] = f_1 * gs_0[k]
                 + pb_z[k] * gp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pb_x, pb_z, dd_1, fp_3, fd_0, \
                         fd_5, gp_3, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * fd_0[k];

        t_7[k] = f_2 * fp_3[k]
                 + pb_x[k] * gp_4[k];

        t_8[k] = pb_z[k] * gp_3[k];

        t_9[k] = f_3 * dd_1[k]
                 + pa_x[k] * fd_5[k];

        t_10[k] = pb_z[k] * gp_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pa_y, pa_z, pb_x, pb_y, fp_4, \
                         fd_0, fd_1, fd_2, gp_5, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_y[k] * fd_2[k];

        t_12[k] = pa_z[k] * fd_0[k];

        t_13[k] = pb_y[k] * gp_5[k];

        t_14[k] = f_2 * fp_4[k]
                  + pb_x[k] * gp_6[k];

        t_15[k] = pa_z[k] * fd_1[k];

        t_16[k] = pb_y[k] * gp_6[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_x, pa_y, pb_x, pb_z, dd_0, dd_2, fp_5, \
                         fd_3, fd_8, gp_7, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * dd_2[k]
                  + pa_x[k] * fd_8[k];

        t_18[k] = f_1 * dd_0[k]
                  + pa_y[k] * fd_3[k];

        t_19[k] = f_3 * fp_5[k]
                  + pb_x[k] * gp_8[k];

        t_20[k] = pb_z[k] * gp_7[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_x, pa_y, pa_z, pb_z, dd_4, fd_4, \
                         fd_6, fd_11, gs_3, gp_8, gp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * dd_4[k]
                  + pa_x[k] * fd_11[k];

        t_22[k] = pb_z[k] * gp_8[k];

        t_23[k] = f_1 * gs_3[k]
                  + pb_z[k] * gp_9[k];

        t_24[k] = pa_y[k] * fd_6[k];

        t_25[k] = pa_z[k] * fd_4[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, dd_0, fp_4, fd_5, \
                         fd_6, fd_7, fd_8, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_y[k] * fd_7[k];

        t_27[k] = pa_z[k] * fd_5[k];

        t_28[k] = f_1 * fp_4[k]
                  + pb_y[k] * gp_10[k];

        t_29[k] = pa_y[k] * fd_8[k];

        t_30[k] = f_1 * dd_0[k]
                  + pa_z[k] * fd_6[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_x, pb_x, pb_y, dd_8, fp_6, fd_14, \
                         gs_4, gp_11, gp_12, gp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * gp_11[k];

        t_32[k] = f_3 * fp_6[k]
                  + pb_x[k] * gp_13[k];

        t_33[k] = f_1 * gs_4[k]
                  + pb_y[k] * gp_12[k];

        t_34[k] = pb_y[k] * gp_13[k];

        t_35[k] = f_1 * dd_8[k]
                  + pa_x[k] * fd_14[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, pa_x, pb_x, pb_z, fp_7, fp_8, \
                         fd_15, fd_16, fd_17, gp_14, gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_3 * fp_7[k]
                  + pa_x[k] * fd_15[k];

        t_37[k] = f_1 * fp_8[k]
                  + pb_x[k] * gp_15[k];

        t_38[k] = pb_z[k] * gp_14[k];

        t_39[k] = pa_x[k] * fd_16[k];

        t_40[k] = pb_z[k] * gp_15[k];

        t_41[k] = pa_x[k] * fd_17[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pa_x, pa_z, pb_x, fp_10, fd_9, \
                         fd_10, fd_18, fd_19, fd_20, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * fd_9[k];

        t_43[k] = pa_z[k] * fd_10[k];

        t_44[k] = f_1 * fp_10[k]
                  + pb_x[k] * gp_16[k];

        t_45[k] = pa_x[k] * fd_18[k];

        t_46[k] = pa_x[k] * fd_19[k];

        t_47[k] = pa_x[k] * fd_20[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, pa_x, pa_y, pb_x, fp_11, fd_12, \
                         fd_13, fd_21, fd_22, fd_23, gp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_y[k] * fd_12[k];

        t_49[k] = f_1 * fp_11[k]
                  + pb_x[k] * gp_17[k];

        t_50[k] = pa_y[k] * fd_13[k];

        t_51[k] = pa_x[k] * fd_21[k];

        t_52[k] = pa_x[k] * fd_22[k];

        t_53[k] = pa_x[k] * fd_23[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, pa_x, pb_x, pb_y, fp_13, fp_15, \
                         fd_24, fd_25, fd_26, gp_18, gp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_3 * fp_13[k]
                  + pa_x[k] * fd_24[k];

        t_55[k] = pb_y[k] * gp_18[k];

        t_56[k] = f_1 * fp_15[k]
                  + pb_x[k] * gp_19[k];

        t_57[k] = pa_x[k] * fd_25[k];

        t_58[k] = pb_y[k] * gp_19[k];

        t_59[k] = pa_x[k] * fd_26[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pb_x, pb_y, pb_z, fp_8, gs_7, \
                         gp_20, gp_21, gp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * gs_7[k]
                  + pb_x[k] * gp_20[k];

        t_61[k] = pb_x[k] * gp_21[k];

        t_62[k] = pb_x[k] * gp_22[k];

        t_63[k] = f_0 * fp_8[k]
                  + f_1 * gs_7[k]
                  + pb_y[k] * gp_21[k];

        t_64[k] = pb_z[k] * gp_21[k];

        t_65[k] = f_1 * gs_7[k]
                  + pb_z[k] * gp_22[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_z, pb_x, pb_y, fp_10, fd_15, fd_16, \
                         gp_23, gp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * fd_15[k];

        t_67[k] = pb_x[k] * gp_23[k];

        t_68[k] = pb_x[k] * gp_24[k];

        t_69[k] = pa_z[k] * fd_16[k];

        t_70[k] = f_2 * fp_10[k]
                  + pb_y[k] * gp_24[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, pa_y, pa_z, pb_x, dd_4, dd_5, fd_18, \
                         fd_20, gs_9, gp_25, gp_26, gp_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_3 * dd_5[k]
                  + pa_y[k] * fd_20[k];

        t_72[k] = f_1 * gs_9[k]
                  + pb_x[k] * gp_25[k];

        t_73[k] = pb_x[k] * gp_26[k];

        t_74[k] = pb_x[k] * gp_27[k];

        t_75[k] = f_1 * dd_4[k]
                  + pa_z[k] * fd_18[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pa_y, pb_x, pb_y, dd_8, fp_12, fd_23, \
                         fd_24, gp_27, gp_28, gp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_3 * fp_12[k]
                  + pb_y[k] * gp_27[k];

        t_77[k] = f_1 * dd_8[k]
                  + pa_y[k] * fd_23[k];

        t_78[k] = pa_y[k] * fd_24[k];

        t_79[k] = pb_x[k] * gp_28[k];

        t_80[k] = pb_x[k] * gp_29[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, pa_y, pb_x, pb_y, fp_14, fp_15, fd_25, \
                         fd_26, gs_11, gp_29, gp_30, gp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_3 * fp_14[k]
                  + pa_y[k] * fd_25[k];

        t_82[k] = f_1 * fp_15[k]
                  + pb_y[k] * gp_29[k];

        t_83[k] = pa_y[k] * fd_26[k];

        t_84[k] = f_1 * gs_11[k]
                  + pb_x[k] * gp_30[k];

        t_85[k] = pb_x[k] * gp_31[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_x, pb_y, pb_z, fp_15, gs_11, gp_31, \
                         gp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_x[k] * gp_32[k];

        t_87[k] = f_1 * gs_11[k]
                  + pb_y[k] * gp_31[k];

        t_88[k] = pb_y[k] * gp_32[k];

        t_89[k] = f_0 * fp_15[k]
                  + f_1 * gs_11[k]
                  + pb_z[k] * gp_32[k];
    }
}

auto
compute_prim_gd_overlap_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dd, const size_t fp, const size_t fd,
                          const size_t gs, const size_t gp, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = 1.0 / p;

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

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_9 = buffer.data(gs + 9);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_y, pb_x, pb_y, pb_z, fp_0, fd_0, gs_0, \
                         gp_0, gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 + f_1 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = pb_z[k] * gp_0[k];

        t_2[k] = f_1 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_3[k] = f_1 * gs_0[k]
                 + pb_z[k] * gp_2[k];

        t_4[k] = pa_y[k] * fd_0[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, pb_x, dd_4, fp_3, fd_0, fd_2, \
                         fd_4, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_2 * fp_3[k]
                 + pb_x[k] * gp_3[k];

        t_6[k] = f_3 * dd_4[k]
                 + pa_x[k] * fd_4[k];

        t_7[k] = pa_y[k] * fd_2[k];

        t_8[k] = pa_z[k] * fd_0[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pa_y, pb_x, pb_y, dd_0, dd_6, fp_4, \
                         fd_3, fd_7, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * fp_4[k]
                 + pb_x[k] * gp_4[k];

        t_10[k] = pb_y[k] * gp_4[k];

        t_11[k] = f_3 * dd_6[k]
                  + pa_x[k] * fd_7[k];

        t_12[k] = f_1 * dd_0[k]
                  + pa_y[k] * fd_3[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pa_z, pb_x, pb_z, dd_8, fp_5, fd_4, \
                         fd_9, gs_3, gp_5, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * fp_5[k]
                  + pb_x[k] * gp_5[k];

        t_14[k] = f_1 * dd_8[k]
                  + pa_x[k] * fd_9[k];

        t_15[k] = f_1 * gs_3[k]
                  + pb_z[k] * gp_6[k];

        t_16[k] = pa_z[k] * fd_4[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pa_z, pb_x, pb_y, dd_0, fp_4, fp_6, \
                         fd_6, fd_7, gp_7, gp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * fp_4[k]
                  + pb_y[k] * gp_7[k];

        t_18[k] = pa_y[k] * fd_7[k];

        t_19[k] = f_1 * dd_0[k]
                  + pa_z[k] * fd_6[k];

        t_20[k] = f_3 * fp_6[k]
                  + pb_x[k] * gp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_x, pb_y, dd_15, fp_7, fd_14, fd_15, gs_4, \
                         gp_8, gp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * gs_4[k]
                  + pb_y[k] * gp_8[k];

        t_22[k] = pb_y[k] * gp_9[k];

        t_23[k] = f_1 * dd_15[k]
                  + pa_x[k] * fd_14[k];

        t_24[k] = f_3 * fp_7[k]
                  + pa_x[k] * fd_15[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, t_30, pa_x, pa_z, pb_x, fp_8, fd_8, \
                         fd_17, fd_18, fd_20, fd_21, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * fp_8[k]
                  + pb_x[k] * gp_10[k];

        t_26[k] = pa_x[k] * fd_17[k];

        t_27[k] = pa_x[k] * fd_18[k];

        t_28[k] = pa_z[k] * fd_8[k];

        t_29[k] = pa_x[k] * fd_20[k];

        t_30[k] = pa_x[k] * fd_21[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_x, pa_y, pb_x, fp_12, fp_14, fd_12, \
                         fd_22, fd_23, fd_25, gp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pa_y[k] * fd_12[k];

        t_32[k] = pa_x[k] * fd_22[k];

        t_33[k] = pa_x[k] * fd_23[k];

        t_34[k] = f_3 * fp_12[k]
                  + pa_x[k] * fd_25[k];

        t_35[k] = f_1 * fp_14[k]
                  + pb_x[k] * gp_11[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, pa_x, pb_x, pb_y, pb_z, fp_8, \
                         fd_27, fd_29, gs_5, gp_12, gp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_x[k] * fd_27[k];

        t_37[k] = pa_x[k] * fd_29[k];

        t_38[k] = f_1 * gs_5[k]
                  + pb_x[k] * gp_12[k];

        t_39[k] = pb_x[k] * gp_13[k];

        t_40[k] = f_0 * fp_8[k]
                  + f_1 * gs_5[k]
                  + pb_y[k] * gp_13[k];

        t_41[k] = pb_z[k] * gp_13[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_z, pb_x, pb_y, pb_z, fp_10, fd_17, gs_5, \
                         gp_14, gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_1 * gs_5[k]
                  + pb_z[k] * gp_14[k];

        t_43[k] = pb_x[k] * gp_15[k];

        t_44[k] = pa_z[k] * fd_17[k];

        t_45[k] = f_2 * fp_10[k]
                  + pb_y[k] * gp_15[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_y, pa_z, pb_x, dd_8, dd_12, fd_19, \
                         fd_21, gs_7, gp_16, gp_17, gp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_3 * dd_12[k]
                  + pa_y[k] * fd_21[k];

        t_47[k] = f_1 * gs_7[k]
                  + pb_x[k] * gp_16[k];

        t_48[k] = pb_x[k] * gp_17[k];

        t_49[k] = pb_x[k] * gp_18[k];

        t_50[k] = f_1 * dd_8[k]
                  + pa_z[k] * fd_19[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pb_x, pb_y, dd_15, fp_11, fp_13, fd_24, \
                         fd_27, gp_18, gp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * fp_11[k]
                  + pb_y[k] * gp_18[k];

        t_52[k] = f_1 * dd_15[k]
                  + pa_y[k] * fd_24[k];

        t_53[k] = pb_x[k] * gp_19[k];

        t_54[k] = f_3 * fp_13[k]
                  + pa_y[k] * fd_27[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, pa_y, pb_x, pb_y, fp_14, fd_29, \
                         gs_9, gp_20, gp_21, gp_22, gp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * fp_14[k]
                  + pb_y[k] * gp_20[k];

        t_56[k] = pa_y[k] * fd_29[k];

        t_57[k] = f_1 * gs_9[k]
                  + pb_x[k] * gp_21[k];

        t_58[k] = pb_x[k] * gp_23[k];

        t_59[k] = f_1 * gs_9[k]
                  + pb_y[k] * gp_22[k];

        t_60[k] = pb_y[k] * gp_23[k];
    }

#pragma omp simd aligned(t_61, pb_z, fp_14, gs_9, gp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_0 * fp_14[k]
                  + f_1 * gs_9[k]
                  + pb_z[k] * gp_23[k];
    }
}

auto
compute_prim_gd_overlap_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dd, const size_t fp, const size_t fd,
                          const size_t gs, const size_t gp, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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
    const auto *fp_11 = buffer.data(fp + 11);

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

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_6 = buffer.data(gs + 6);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, fp_0, fd_0, gs_0, gp_0, \
                         gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 + f_1 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = f_1 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = f_1 * gs_0[k]
                 + pb_z[k] * gp_2[k];

        t_3[k] = pa_y[k] * fd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, pb_y, dd_3, dd_5, fd_0, \
                         fd_2, fd_4, fd_8, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * dd_3[k]
                 + pa_x[k] * fd_4[k];

        t_5[k] = pa_y[k] * fd_2[k];

        t_6[k] = pa_z[k] * fd_0[k];

        t_7[k] = pb_y[k] * gp_3[k];

        t_8[k] = f_2 * dd_5[k]
                 + pa_x[k] * fd_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pa_y, pa_z, pb_z, dd_0, dd_7, fd_3, \
                         fd_4, fd_10, gs_2, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_1 * dd_0[k]
                 + pa_y[k] * fd_3[k];

        t_10[k] = f_1 * dd_7[k]
                  + pa_x[k] * fd_10[k];

        t_11[k] = f_1 * gs_2[k]
                  + pb_z[k] * gp_4[k];

        t_12[k] = pa_z[k] * fd_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pa_z, pb_y, dd_0, fd_6, fd_8, gs_3, \
                         gp_5, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * fd_8[k];

        t_14[k] = f_1 * dd_0[k]
                  + pa_z[k] * fd_6[k];

        t_15[k] = f_1 * gs_3[k]
                  + pb_y[k] * gp_5[k];

        t_16[k] = pb_y[k] * gp_6[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_x, pa_z, dd_14, fp_4, fp_9, fd_9, \
                         fd_12, fd_13, fd_15, fd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * dd_14[k]
                  + pa_x[k] * fd_12[k];

        t_18[k] = f_2 * fp_4[k]
                  + pa_x[k] * fd_13[k];

        t_19[k] = pa_x[k] * fd_15[k];

        t_20[k] = pa_z[k] * fd_9[k];

        t_21[k] = f_2 * fp_9[k]
                  + pa_x[k] * fd_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pa_x, pb_x, pb_y, pb_z, fp_5, fd_25, \
                         gs_4, gp_7, gp_8, gp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pa_x[k] * fd_25[k];

        t_23[k] = f_1 * gs_4[k]
                  + pb_x[k] * gp_7[k];

        t_24[k] = pb_x[k] * gp_8[k];

        t_25[k] = f_0 * fp_5[k]
                  + f_1 * gs_4[k]
                  + pb_y[k] * gp_8[k];

        t_26[k] = f_1 * gs_4[k]
                  + pb_z[k] * gp_9[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_y, pa_z, pb_x, pb_y, dd_10, fp_7, \
                         fd_15, fd_18, gs_6, gp_10, gp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_x[k] * gp_10[k];

        t_28[k] = pa_z[k] * fd_15[k];

        t_29[k] = f_3 * fp_7[k]
                  + pb_y[k] * gp_10[k];

        t_30[k] = f_2 * dd_10[k]
                  + pa_y[k] * fd_18[k];

        t_31[k] = f_1 * gs_6[k]
                  + pb_x[k] * gp_11[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, pa_y, pa_z, pb_x, pb_y, dd_7, dd_14, \
                         fp_8, fd_17, fd_20, gp_12, gp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_x[k] * gp_12[k];

        t_33[k] = pb_x[k] * gp_13[k];

        t_34[k] = f_1 * dd_7[k]
                  + pa_z[k] * fd_17[k];

        t_35[k] = f_2 * fp_8[k]
                  + pb_y[k] * gp_13[k];

        t_36[k] = f_1 * dd_14[k]
                  + pa_y[k] * fd_20[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pa_y, pb_x, pb_y, fp_10, fp_11, fd_23, \
                         fd_25, gs_8, gp_14, gp_15, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_x[k] * gp_14[k];

        t_38[k] = f_2 * fp_10[k]
                  + pa_y[k] * fd_23[k];

        t_39[k] = f_1 * fp_11[k]
                  + pb_y[k] * gp_15[k];

        t_40[k] = pa_y[k] * fd_25[k];

        t_41[k] = f_1 * gs_8[k]
                  + pb_x[k] * gp_16[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_x, pb_y, pb_z, fp_11, gs_8, gp_17, \
                         gp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_x[k] * gp_18[k];

        t_43[k] = f_1 * gs_8[k]
                  + pb_y[k] * gp_17[k];

        t_44[k] = pb_y[k] * gp_18[k];

        t_45[k] = f_0 * fp_11[k]
                  + f_1 * gs_8[k]
                  + pb_z[k] * gp_18[k];
    }
}

auto
compute_prim_gd_overlap_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dd, const size_t fp, const size_t fd,
                          const size_t gs, const size_t gp, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_12 = buffer.data(dd + 12);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_22 = buffer.data(fd + 22);
    const auto *fd_24 = buffer.data(fd + 24);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_6 = buffer.data(gs + 6);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pb_x, pb_y, pb_z, dd_3, fp_0, fd_4, gs_0, \
                         gp_0, gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 + f_1 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = f_1 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = f_1 * gs_0[k]
                 + pb_z[k] * gp_2[k];

        t_3[k] = f_2 * dd_3[k]
                 + pa_x[k] * fd_4[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, pb_y, dd_0, dd_4, fd_0, fd_3, \
                         fd_7, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * fd_0[k];

        t_5[k] = pb_y[k] * gp_3[k];

        t_6[k] = f_2 * dd_4[k]
                 + pa_x[k] * fd_7[k];

        t_7[k] = f_1 * dd_0[k]
                 + pa_y[k] * fd_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_z, pb_y, pb_z, dd_0, dd_6, fd_5, fd_9, \
                         gs_2, gs_3, gp_4, gp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * dd_6[k]
                 + pa_x[k] * fd_9[k];

        t_9[k] = f_1 * gs_2[k]
                 + pb_z[k] * gp_4[k];

        t_10[k] = f_1 * dd_0[k]
                  + pa_z[k] * fd_5[k];

        t_11[k] = f_1 * gs_3[k]
                  + pb_y[k] * gp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_y, dd_12, fp_4, fp_9, fd_11, fd_12, \
                         fd_20, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_y[k] * gp_6[k];

        t_13[k] = f_1 * dd_12[k]
                  + pa_x[k] * fd_11[k];

        t_14[k] = f_2 * fp_4[k]
                  + pa_x[k] * fd_12[k];

        t_15[k] = f_2 * fp_9[k]
                  + pa_x[k] * fd_20[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pb_x, pb_y, pb_z, fp_5, gs_4, gp_7, \
                         gp_8, gp_9, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * gs_4[k]
                  + pb_x[k] * gp_7[k];

        t_17[k] = pb_x[k] * gp_8[k];

        t_18[k] = f_0 * fp_5[k]
                  + f_1 * gs_4[k]
                  + pb_y[k] * gp_8[k];

        t_19[k] = f_1 * gs_4[k]
                  + pb_z[k] * gp_9[k];

        t_20[k] = pb_x[k] * gp_10[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_y, pa_z, pb_x, pb_y, dd_8, fp_7, fd_14, \
                         fd_17, gs_6, gp_10, gp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_z[k] * fd_14[k];

        t_22[k] = f_3 * fp_7[k]
                  + pb_y[k] * gp_10[k];

        t_23[k] = f_2 * dd_8[k]
                  + pa_y[k] * fd_17[k];

        t_24[k] = f_1 * gs_6[k]
                  + pb_x[k] * gp_11[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pa_z, pb_x, pb_y, dd_6, dd_12, \
                         fp_8, fd_16, fd_19, gp_12, gp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pb_x[k] * gp_12[k];

        t_26[k] = pb_x[k] * gp_13[k];

        t_27[k] = f_1 * dd_6[k]
                  + pa_z[k] * fd_16[k];

        t_28[k] = f_2 * fp_8[k]
                  + pb_y[k] * gp_13[k];

        t_29[k] = f_1 * dd_12[k]
                  + pa_y[k] * fd_19[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pa_y, pb_x, pb_y, fp_10, fp_11, fd_22, \
                         fd_24, gs_8, gp_14, gp_15, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_x[k] * gp_14[k];

        t_31[k] = f_2 * fp_10[k]
                  + pa_y[k] * fd_22[k];

        t_32[k] = f_1 * fp_11[k]
                  + pb_y[k] * gp_15[k];

        t_33[k] = pa_y[k] * fd_24[k];

        t_34[k] = f_1 * gs_8[k]
                  + pb_x[k] * gp_16[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, fp_11, gs_8, gp_17, \
                         gp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pb_x[k] * gp_18[k];

        t_36[k] = f_1 * gs_8[k]
                  + pb_y[k] * gp_17[k];

        t_37[k] = pb_y[k] * gp_18[k];

        t_38[k] = f_0 * fp_11[k]
                  + f_1 * gs_8[k]
                  + pb_z[k] * gp_18[k];
    }
}

auto
compute_prim_gd_overlap_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dd, const size_t fp, const size_t fd,
                          const size_t gs, const size_t gp, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = 1.0 / p;

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, fp_0, fd_0, gs_0, gp_0, \
                         gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 + f_1 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = f_1 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = f_1 * gs_0[k]
                 + pb_z[k] * gp_2[k];

        t_3[k] = pa_y[k] * fd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_z, pb_x, dd_1, fp_1, fp_2, fd_0, fd_2, \
                         gp_3, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * fp_1[k]
                 + pb_x[k] * gp_3[k];

        t_5[k] = f_3 * dd_1[k]
                 + pa_x[k] * fd_2[k];

        t_6[k] = pa_z[k] * fd_0[k];

        t_7[k] = f_2 * fp_2[k]
                 + pb_x[k] * gp_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_y, pb_x, dd_0, dd_2, dd_3, fp_3, fd_1, \
                         fd_4, fd_5, gp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_3 * dd_2[k]
                 + pa_x[k] * fd_4[k];

        t_9[k] = f_1 * dd_0[k]
                 + pa_y[k] * fd_1[k];

        t_10[k] = f_3 * fp_3[k]
                  + pb_x[k] * gp_5[k];

        t_11[k] = f_1 * dd_3[k]
                  + pa_x[k] * fd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pa_z, pb_x, dd_0, dd_6, fp_4, fp_5, \
                         fd_3, fd_6, fd_7, gp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_1 * dd_0[k]
                  + pa_z[k] * fd_3[k];

        t_13[k] = f_3 * fp_4[k]
                  + pb_x[k] * gp_7[k];

        t_14[k] = f_1 * dd_6[k]
                  + pa_x[k] * fd_6[k];

        t_15[k] = f_3 * fp_5[k]
                  + pa_x[k] * fd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_x, pb_x, fp_6, fp_10, fd_8, fd_10, \
                         fd_11, fd_13, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * fp_6[k]
                  + pb_x[k] * gp_8[k];

        t_17[k] = pa_x[k] * fd_8[k];

        t_18[k] = pa_x[k] * fd_10[k];

        t_19[k] = pa_x[k] * fd_11[k];

        t_20[k] = f_3 * fp_10[k]
                  + pa_x[k] * fd_13[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_x, pb_x, pb_y, fp_6, fp_12, fd_15, gs_7, \
                         gp_11, gp_12, gp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * fp_12[k]
                  + pb_x[k] * gp_11[k];

        t_22[k] = pa_x[k] * fd_15[k];

        t_23[k] = f_1 * gs_7[k]
                  + pb_x[k] * gp_12[k];

        t_24[k] = f_0 * fp_6[k]
                  + f_1 * gs_7[k]
                  + pb_y[k] * gp_13[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, pa_z, pb_y, pb_z, dd_4, fp_7, fd_8, \
                         fd_10, gs_7, gp_14, gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * gs_7[k]
                  + pb_z[k] * gp_14[k];

        t_26[k] = pa_z[k] * fd_8[k];

        t_27[k] = f_2 * fp_7[k]
                  + pb_y[k] * gp_15[k];

        t_28[k] = f_3 * dd_4[k]
                  + pa_y[k] * fd_10[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pb_x, pb_y, dd_3, dd_6, fp_9, \
                         fd_9, fd_12, gs_9, gp_16, gp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * gs_9[k]
                  + pb_x[k] * gp_16[k];

        t_30[k] = f_1 * dd_3[k]
                  + pa_z[k] * fd_9[k];

        t_31[k] = f_3 * fp_9[k]
                  + pb_y[k] * gp_18[k];

        t_32[k] = f_1 * dd_6[k]
                  + pa_y[k] * fd_12[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pb_x, pb_y, fp_11, fp_12, fd_14, \
                         fd_15, gs_11, gp_20, gp_21, gp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * fp_11[k]
                  + pa_y[k] * fd_14[k];

        t_34[k] = f_1 * fp_12[k]
                  + pb_y[k] * gp_20[k];

        t_35[k] = pa_y[k] * fd_15[k];

        t_36[k] = f_1 * gs_11[k]
                  + pb_x[k] * gp_21[k];

        t_37[k] = f_1 * gs_11[k]
                  + pb_y[k] * gp_22[k];
    }

#pragma omp simd aligned(t_38, pb_z, fp_12, gs_11, gp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * fp_12[k]
                  + f_1 * gs_11[k]
                  + pb_z[k] * gp_23[k];
    }
}

auto
compute_prim_gd_overlap_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dd, const size_t fp, const size_t fd,
                          const size_t gs, const size_t gp, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, fp_0, fd_0, gs_0, gp_0, \
                         gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 + f_1 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = f_1 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = f_1 * gs_0[k]
                 + pb_z[k] * gp_2[k];

        t_3[k] = pa_y[k] * fd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, dd_0, dd_1, dd_2, fd_0, \
                         fd_2, fd_3, fd_4, fd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * dd_1[k]
                 + pa_x[k] * fd_4[k];

        t_5[k] = pa_y[k] * fd_2[k];

        t_6[k] = pa_z[k] * fd_0[k];

        t_7[k] = f_2 * dd_2[k]
                 + pa_x[k] * fd_6[k];

        t_8[k] = f_1 * dd_0[k]
                 + pa_y[k] * fd_3[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pa_z, pb_y, pb_z, dd_4, fp_3, fd_4, \
                         fd_8, gs_3, gp_6, gp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_1 * dd_4[k]
                 + pa_x[k] * fd_8[k];

        t_10[k] = f_1 * gs_3[k]
                  + pb_z[k] * gp_6[k];

        t_11[k] = pa_z[k] * fd_4[k];

        t_12[k] = f_1 * fp_3[k]
                  + pb_y[k] * gp_7[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pa_y, pa_z, pb_y, dd_0, dd_8, fd_5, \
                         fd_6, fd_10, gs_4, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * fd_6[k];

        t_14[k] = f_1 * dd_0[k]
                  + pa_z[k] * fd_5[k];

        t_15[k] = f_1 * gs_4[k]
                  + pb_y[k] * gp_8[k];

        t_16[k] = f_1 * dd_8[k]
                  + pa_x[k] * fd_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, t_22, pa_x, pa_z, fp_4, fd_7, fd_11, \
                         fd_12, fd_13, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_2 * fp_4[k]
                  + pa_x[k] * fd_11[k];

        t_18[k] = pa_x[k] * fd_12[k];

        t_19[k] = pa_x[k] * fd_13[k];

        t_20[k] = pa_z[k] * fd_7[k];

        t_21[k] = pa_x[k] * fd_15[k];

        t_22[k] = pa_x[k] * fd_16[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, t_28, pa_x, pa_y, fp_9, fd_9, fd_17, \
                         fd_18, fd_20, fd_21, fd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_y[k] * fd_9[k];

        t_24[k] = pa_x[k] * fd_17[k];

        t_25[k] = pa_x[k] * fd_18[k];

        t_26[k] = f_2 * fp_9[k]
                  + pa_x[k] * fd_20[k];

        t_27[k] = pa_x[k] * fd_21[k];

        t_28[k] = pa_x[k] * fd_22[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pa_z, pb_x, pb_y, pb_z, fp_5, fd_12, \
                         gs_7, gp_12, gp_13, gp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * gs_7[k]
                  + pb_x[k] * gp_12[k];

        t_30[k] = pb_x[k] * gp_13[k];

        t_31[k] = f_0 * fp_5[k]
                  + f_1 * gs_7[k]
                  + pb_y[k] * gp_13[k];

        t_32[k] = f_1 * gs_7[k]
                  + pb_z[k] * gp_14[k];

        t_33[k] = pa_z[k] * fd_12[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pa_z, pb_x, pb_y, dd_4, dd_5, fp_7, \
                         fd_14, fd_16, gs_9, gp_15, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * fp_7[k]
                  + pb_y[k] * gp_15[k];

        t_35[k] = f_2 * dd_5[k]
                  + pa_y[k] * fd_16[k];

        t_36[k] = f_1 * gs_9[k]
                  + pb_x[k] * gp_16[k];

        t_37[k] = f_1 * dd_4[k]
                  + pa_z[k] * fd_14[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pb_y, dd_8, fp_8, fp_10, fp_11, \
                         fd_19, fd_21, fd_22, gp_17, gp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_2 * fp_8[k]
                  + pb_y[k] * gp_17[k];

        t_39[k] = f_1 * dd_8[k]
                  + pa_y[k] * fd_19[k];

        t_40[k] = f_2 * fp_10[k]
                  + pa_y[k] * fd_21[k];

        t_41[k] = f_1 * fp_11[k]
                  + pb_y[k] * gp_18[k];

        t_42[k] = pa_y[k] * fd_22[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pb_x, pb_y, pb_z, fp_11, gs_11, gp_19, \
                         gp_20, gp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * gs_11[k]
                  + pb_x[k] * gp_19[k];

        t_44[k] = pb_x[k] * gp_21[k];

        t_45[k] = f_1 * gs_11[k]
                  + pb_y[k] * gp_20[k];

        t_46[k] = pb_y[k] * gp_21[k];

        t_47[k] = f_0 * fp_11[k]
                  + f_1 * gs_11[k]
                  + pb_z[k] * gp_21[k];
    }
}

auto
compute_prim_gd_overlap_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dd, const size_t fp, const size_t fd,
                          const size_t gs, const size_t gp, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;

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
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, fp_0, fd_0, gs_0, gp_0, \
                         gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 + f_1 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = f_1 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = f_1 * gs_0[k]
                 + pb_z[k] * gp_2[k];

        t_3[k] = pa_y[k] * fd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, pb_y, dd_2, dd_4, fd_0, \
                         fd_2, fd_4, fd_7, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * dd_2[k]
                 + pa_x[k] * fd_4[k];

        t_5[k] = pa_y[k] * fd_2[k];

        t_6[k] = pa_z[k] * fd_0[k];

        t_7[k] = pb_y[k] * gp_3[k];

        t_8[k] = f_2 * dd_4[k]
                 + pa_x[k] * fd_7[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pa_y, pa_z, pb_z, dd_0, dd_6, fd_3, \
                         fd_4, fd_9, gs_2, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_1 * dd_0[k]
                 + pa_y[k] * fd_3[k];

        t_10[k] = f_1 * dd_6[k]
                  + pa_x[k] * fd_9[k];

        t_11[k] = f_1 * gs_2[k]
                  + pb_z[k] * gp_4[k];

        t_12[k] = pa_z[k] * fd_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pa_z, pb_y, dd_0, fd_6, fd_7, gs_3, \
                         gp_5, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * fd_7[k];

        t_14[k] = f_1 * dd_0[k]
                  + pa_z[k] * fd_6[k];

        t_15[k] = f_1 * gs_3[k]
                  + pb_y[k] * gp_5[k];

        t_16[k] = pb_y[k] * gp_6[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, t_22, pa_x, pa_z, dd_12, fp_3, fd_8, \
                         fd_11, fd_12, fd_13, fd_16, fd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * dd_12[k]
                  + pa_x[k] * fd_11[k];

        t_18[k] = f_2 * fp_3[k]
                  + pa_x[k] * fd_12[k];

        t_19[k] = pa_x[k] * fd_13[k];

        t_20[k] = pa_z[k] * fd_8[k];

        t_21[k] = pa_x[k] * fd_16[k];

        t_22[k] = pa_x[k] * fd_17[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pa_x, pb_x, pb_y, fp_4, fp_6, fd_19, \
                         fd_22, gs_6, gp_7, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_2 * fp_6[k]
                  + pa_x[k] * fd_19[k];

        t_24[k] = pa_x[k] * fd_22[k];

        t_25[k] = f_1 * gs_6[k]
                  + pb_x[k] * gp_7[k];

        t_26[k] = pb_x[k] * gp_8[k];

        t_27[k] = f_0 * fp_4[k]
                  + f_1 * gs_6[k]
                  + pb_y[k] * gp_8[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pa_z, pb_x, pb_z, dd_9, fd_13, fd_16, \
                         gs_6, gs_8, gp_9, gp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * gs_6[k]
                  + pb_z[k] * gp_9[k];

        t_29[k] = pa_z[k] * fd_13[k];

        t_30[k] = f_2 * dd_9[k]
                  + pa_y[k] * fd_16[k];

        t_31[k] = f_1 * gs_8[k]
                  + pb_x[k] * gp_11[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_y, pa_z, dd_6, dd_12, fp_7, fd_15, fd_18, \
                         fd_20, fd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * dd_6[k]
                  + pa_z[k] * fd_15[k];

        t_33[k] = f_1 * dd_12[k]
                  + pa_y[k] * fd_18[k];

        t_34[k] = f_2 * fp_7[k]
                  + pa_y[k] * fd_20[k];

        t_35[k] = pa_y[k] * fd_22[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pb_x, pb_y, pb_z, fp_8, gs_10, gp_14, \
                         gp_15, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * gs_10[k]
                  + pb_x[k] * gp_14[k];

        t_37[k] = pb_x[k] * gp_16[k];

        t_38[k] = f_1 * gs_10[k]
                  + pb_y[k] * gp_15[k];

        t_39[k] = pb_y[k] * gp_16[k];

        t_40[k] = f_0 * fp_8[k]
                  + f_1 * gs_10[k]
                  + pb_z[k] * gp_16[k];
    }
}

auto
compute_prim_gd_overlap_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dd, const size_t fp, const size_t fd,
                          const size_t gs, const size_t gp, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_10 = buffer.data(dd + 10);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_19 = buffer.data(fd + 19);

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
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, fp_0, fd_0, gs_0, gp_0, \
                         gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 + f_1 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = f_1 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = f_1 * gs_0[k]
                 + pb_z[k] * gp_2[k];

        t_3[k] = pa_y[k] * fd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_z, pb_y, dd_2, dd_3, fd_0, fd_4, fd_6, \
                         gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * dd_2[k]
                 + pa_x[k] * fd_4[k];

        t_5[k] = pa_z[k] * fd_0[k];

        t_6[k] = pb_y[k] * gp_3[k];

        t_7[k] = f_2 * dd_3[k]
                 + pa_x[k] * fd_6[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_y, pa_z, pb_z, dd_0, dd_5, fd_3, fd_5, \
                         fd_7, gs_2, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * dd_0[k]
                 + pa_y[k] * fd_3[k];

        t_9[k] = f_1 * dd_5[k]
                 + pa_x[k] * fd_7[k];

        t_10[k] = f_1 * gs_2[k]
                  + pb_z[k] * gp_4[k];

        t_11[k] = f_1 * dd_0[k]
                  + pa_z[k] * fd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pa_x, pb_y, dd_10, fp_3, fd_8, fd_9, \
                         fd_10, gs_3, gp_5, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_1 * gs_3[k]
                  + pb_y[k] * gp_5[k];

        t_13[k] = pb_y[k] * gp_6[k];

        t_14[k] = f_1 * dd_10[k]
                  + pa_x[k] * fd_8[k];

        t_15[k] = f_2 * fp_3[k]
                  + pa_x[k] * fd_9[k];

        t_16[k] = pa_x[k] * fd_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_x, pb_x, pb_y, fp_4, fp_6, fd_16, \
                         fd_19, gs_6, gp_7, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_2 * fp_6[k]
                  + pa_x[k] * fd_16[k];

        t_18[k] = pa_x[k] * fd_19[k];

        t_19[k] = f_1 * gs_6[k]
                  + pb_x[k] * gp_7[k];

        t_20[k] = pb_x[k] * gp_8[k];

        t_21[k] = f_0 * fp_4[k]
                  + f_1 * gs_6[k]
                  + pb_y[k] * gp_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_y, pa_z, pb_x, pb_z, dd_7, fd_10, fd_13, \
                         gs_6, gs_8, gp_9, gp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_1 * gs_6[k]
                  + pb_z[k] * gp_9[k];

        t_23[k] = pa_z[k] * fd_10[k];

        t_24[k] = f_2 * dd_7[k]
                  + pa_y[k] * fd_13[k];

        t_25[k] = f_1 * gs_8[k]
                  + pb_x[k] * gp_11[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, dd_5, dd_10, fp_7, fd_12, fd_15, \
                         fd_17, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * dd_5[k]
                  + pa_z[k] * fd_12[k];

        t_27[k] = f_1 * dd_10[k]
                  + pa_y[k] * fd_15[k];

        t_28[k] = f_2 * fp_7[k]
                  + pa_y[k] * fd_17[k];

        t_29[k] = pa_y[k] * fd_19[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_x, pb_y, pb_z, fp_8, gs_10, gp_14, \
                         gp_15, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_1 * gs_10[k]
                  + pb_x[k] * gp_14[k];

        t_31[k] = pb_x[k] * gp_16[k];

        t_32[k] = f_1 * gs_10[k]
                  + pb_y[k] * gp_15[k];

        t_33[k] = pb_y[k] * gp_16[k];

        t_34[k] = f_0 * fp_8[k]
                  + f_1 * gs_10[k]
                  + pb_z[k] * gp_16[k];
    }
}

auto
compute_prim_gd_overlap_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dd, const size_t fp, const size_t fd,
                          const size_t gs, const size_t gp, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;

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

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_11 = buffer.data(gs + 11);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, dd_1, fp_0, fd_0, fd_2, \
                         gs_0, gp_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 + f_1 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = pa_y[k] * fd_0[k];

        t_2[k] = f_2 * dd_1[k]
                 + pa_x[k] * fd_2[k];

        t_3[k] = pa_z[k] * fd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, dd_0, dd_2, dd_3, fd_1, fd_3, \
                         fd_4, fd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * dd_2[k]
                 + pa_x[k] * fd_4[k];

        t_5[k] = f_1 * dd_0[k]
                 + pa_y[k] * fd_1[k];

        t_6[k] = f_1 * dd_3[k]
                 + pa_x[k] * fd_5[k];

        t_7[k] = f_1 * dd_0[k]
                 + pa_z[k] * fd_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, pa_x, pb_x, dd_6, fd_6, fd_7, fd_9, \
                         fd_10, fd_13, gs_7, gp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * dd_6[k]
                 + pa_x[k] * fd_6[k];

        t_9[k] = pa_x[k] * fd_7[k];

        t_10[k] = pa_x[k] * fd_9[k];

        t_11[k] = pa_x[k] * fd_10[k];

        t_12[k] = pa_x[k] * fd_13[k];

        t_13[k] = f_1 * gs_7[k]
                  + pb_x[k] * gp_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pa_z, pb_y, dd_3, dd_4, fp_1, fd_7, \
                         fd_8, fd_9, gs_7, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * fp_1[k]
                  + f_1 * gs_7[k]
                  + pb_y[k] * gp_8[k];

        t_15[k] = pa_z[k] * fd_7[k];

        t_16[k] = f_2 * dd_4[k]
                  + pa_y[k] * fd_9[k];

        t_17[k] = f_1 * dd_3[k]
                  + pa_z[k] * fd_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_y, pb_x, pb_y, dd_6, fp_2, fd_11, \
                         fd_12, fd_13, gs_11, gp_14, gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_1 * dd_6[k]
                  + pa_y[k] * fd_11[k];

        t_19[k] = f_2 * fp_2[k]
                  + pa_y[k] * fd_12[k];

        t_20[k] = pa_y[k] * fd_13[k];

        t_21[k] = f_1 * gs_11[k]
                  + pb_x[k] * gp_14[k];

        t_22[k] = f_1 * gs_11[k]
                  + pb_y[k] * gp_15[k];
    }

#pragma omp simd aligned(t_23, pb_z, fp_3, gs_11, gp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * fp_3[k]
                  + f_1 * gs_11[k]
                  + pb_z[k] * gp_16[k];
    }
}

auto
compute_prim_gd_overlap_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t dd, const size_t fp, const size_t fd,
                          const size_t gs, const size_t gp, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, fp_0, fd_0, gs_0, gp_0, \
                         gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 + f_1 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = f_1 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = f_1 * gs_0[k]
                 + pb_z[k] * gp_2[k];

        t_3[k] = pa_y[k] * fd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, dd_0, dd_1, dd_2, fd_0, fd_1, \
                         fd_2, fd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * dd_1[k]
                 + pa_x[k] * fd_2[k];

        t_5[k] = pa_z[k] * fd_0[k];

        t_6[k] = f_2 * dd_2[k]
                 + pa_x[k] * fd_4[k];

        t_7[k] = f_1 * dd_0[k]
                 + pa_y[k] * fd_1[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, pa_x, pa_z, dd_0, dd_3, dd_6, fp_2, fd_3, \
                         fd_5, fd_6, fd_7, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * dd_3[k]
                 + pa_x[k] * fd_5[k];

        t_9[k] = f_1 * dd_0[k]
                 + pa_z[k] * fd_3[k];

        t_10[k] = f_1 * dd_6[k]
                  + pa_x[k] * fd_6[k];

        t_11[k] = f_2 * fp_2[k]
                  + pa_x[k] * fd_7[k];

        t_12[k] = pa_x[k] * fd_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_x, pb_x, fp_6, fd_10, fd_11, fd_13, \
                         fd_15, gs_7, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_x[k] * fd_10[k];

        t_14[k] = pa_x[k] * fd_11[k];

        t_15[k] = f_2 * fp_6[k]
                  + pa_x[k] * fd_13[k];

        t_16[k] = pa_x[k] * fd_15[k];

        t_17[k] = f_1 * gs_7[k]
                  + pb_x[k] * gp_6[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_z, pb_y, pb_z, fp_3, fp_4, fd_8, gs_7, \
                         gp_7, gp_8, gp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * fp_3[k]
                  + f_1 * gs_7[k]
                  + pb_y[k] * gp_7[k];

        t_19[k] = f_1 * gs_7[k]
                  + pb_z[k] * gp_8[k];

        t_20[k] = pa_z[k] * fd_8[k];

        t_21[k] = f_3 * fp_4[k]
                  + pb_y[k] * gp_9[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_y, pa_z, pb_x, pb_y, dd_3, dd_4, fp_5, \
                         fd_9, fd_10, gs_9, gp_10, gp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_2 * dd_4[k]
                  + pa_y[k] * fd_10[k];

        t_23[k] = f_1 * gs_9[k]
                  + pb_x[k] * gp_10[k];

        t_24[k] = f_1 * dd_3[k]
                  + pa_z[k] * fd_9[k];

        t_25[k] = f_2 * fp_5[k]
                  + pb_y[k] * gp_11[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pb_y, dd_6, fp_7, fp_8, fd_12, fd_14, \
                         fd_15, gp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * dd_6[k]
                  + pa_y[k] * fd_12[k];

        t_27[k] = f_2 * fp_7[k]
                  + pa_y[k] * fd_14[k];

        t_28[k] = f_1 * fp_8[k]
                  + pb_y[k] * gp_12[k];

        t_29[k] = pa_y[k] * fd_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_x, pb_y, pb_z, fp_8, gs_11, gp_13, gp_14, \
                         gp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_1 * gs_11[k]
                  + pb_x[k] * gp_13[k];

        t_31[k] = f_1 * gs_11[k]
                  + pb_y[k] * gp_14[k];

        t_32[k] = f_0 * fp_8[k]
                  + f_1 * gs_11[k]
                  + pb_z[k] * gp_15[k];
    }
}

auto
compute_prim_gd_overlap_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                           const size_t pb, const size_t dd, const size_t fp, const size_t fd,
                           const size_t gs, const size_t gp, const size_t ncols,
                           const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, fp_0, fd_0, gs_0, gp_0, \
                         gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 + f_1 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = f_1 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = f_1 * gs_0[k]
                 + pb_z[k] * gp_2[k];

        t_3[k] = pa_y[k] * fd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, dd_0, dd_1, dd_2, fd_0, \
                         fd_1, fd_2, fd_3, fd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * dd_1[k]
                 + pa_x[k] * fd_3[k];

        t_5[k] = pa_y[k] * fd_1[k];

        t_6[k] = pa_z[k] * fd_0[k];

        t_7[k] = f_2 * dd_2[k]
                 + pa_x[k] * fd_5[k];

        t_8[k] = f_1 * dd_0[k]
                 + pa_y[k] * fd_2[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pa_y, pa_z, pb_z, dd_4, fd_3, fd_5, \
                         fd_7, gs_2, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_1 * dd_4[k]
                 + pa_x[k] * fd_7[k];

        t_10[k] = f_1 * gs_2[k]
                  + pb_z[k] * gp_3[k];

        t_11[k] = pa_z[k] * fd_3[k];

        t_12[k] = pa_y[k] * fd_5[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_x, pa_z, dd_0, dd_8, fp_2, fd_4, \
                         fd_6, fd_9, fd_10, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * dd_0[k]
                  + pa_z[k] * fd_4[k];

        t_14[k] = f_1 * dd_8[k]
                  + pa_x[k] * fd_9[k];

        t_15[k] = f_2 * fp_2[k]
                  + pa_x[k] * fd_10[k];

        t_16[k] = pa_x[k] * fd_11[k];

        t_17[k] = pa_z[k] * fd_6[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_x, pb_x, fp_5, fd_14, fd_15, fd_17, \
                         fd_19, gs_6, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pa_x[k] * fd_14[k];

        t_19[k] = pa_x[k] * fd_15[k];

        t_20[k] = f_2 * fp_5[k]
                  + pa_x[k] * fd_17[k];

        t_21[k] = pa_x[k] * fd_19[k];

        t_22[k] = f_1 * gs_6[k]
                  + pb_x[k] * gp_4[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pa_z, pb_y, pb_z, dd_5, fp_3, fd_11, \
                         fd_14, gs_6, gp_5, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * fp_3[k]
                  + f_1 * gs_6[k]
                  + pb_y[k] * gp_5[k];

        t_24[k] = f_1 * gs_6[k]
                  + pb_z[k] * gp_6[k];

        t_25[k] = pa_z[k] * fd_11[k];

        t_26[k] = f_2 * dd_5[k]
                  + pa_y[k] * fd_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_y, pa_z, pb_x, dd_4, dd_8, fp_6, fd_13, \
                         fd_16, fd_18, gs_8, gp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * gs_8[k]
                  + pb_x[k] * gp_7[k];

        t_28[k] = f_1 * dd_4[k]
                  + pa_z[k] * fd_13[k];

        t_29[k] = f_1 * dd_8[k]
                  + pa_y[k] * fd_16[k];

        t_30[k] = f_2 * fp_6[k]
                  + pa_y[k] * fd_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_y, pb_x, pb_y, pb_z, fp_7, fd_19, \
                         gs_10, gp_8, gp_9, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pa_y[k] * fd_19[k];

        t_32[k] = f_1 * gs_10[k]
                  + pb_x[k] * gp_8[k];

        t_33[k] = f_1 * gs_10[k]
                  + pb_y[k] * gp_9[k];

        t_34[k] = pb_y[k] * gp_10[k];

        t_35[k] = f_0 * fp_7[k]
                  + f_1 * gs_10[k]
                  + pb_z[k] * gp_10[k];
    }
}

auto
compute_prim_gd_overlap_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                           const size_t pb, const size_t dd, const size_t fp, const size_t fd,
                           const size_t gs, const size_t gp, const size_t ncols,
                           const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, fp_0, fd_0, gs_0, gp_0, \
                         gp_1, gp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fp_0[k]
                 + f_1 * gs_0[k]
                 + pb_x[k] * gp_0[k];

        t_1[k] = f_1 * gs_0[k]
                 + pb_y[k] * gp_1[k];

        t_2[k] = f_1 * gs_0[k]
                 + pb_z[k] * gp_2[k];

        t_3[k] = pa_y[k] * fd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, dd_0, dd_1, dd_2, fd_0, fd_2, \
                         fd_3, fd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * dd_1[k]
                 + pa_x[k] * fd_3[k];

        t_5[k] = pa_z[k] * fd_0[k];

        t_6[k] = f_2 * dd_2[k]
                 + pa_x[k] * fd_5[k];

        t_7[k] = f_1 * dd_0[k]
                 + pa_y[k] * fd_2[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_z, pb_z, dd_0, dd_4, dd_8, fd_4, fd_6, \
                         fd_7, gs_2, gp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * dd_4[k]
                 + pa_x[k] * fd_6[k];

        t_9[k] = f_1 * gs_2[k]
                 + pb_z[k] * gp_3[k];

        t_10[k] = f_1 * dd_0[k]
                  + pa_z[k] * fd_4[k];

        t_11[k] = f_1 * dd_8[k]
                  + pa_x[k] * fd_7[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pa_x, pb_x, pb_y, pb_z, fp_3, fd_9, \
                         fd_16, gs_6, gp_4, gp_5, gp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_x[k] * fd_9[k];

        t_13[k] = pa_x[k] * fd_16[k];

        t_14[k] = f_1 * gs_6[k]
                  + pb_x[k] * gp_4[k];

        t_15[k] = f_0 * fp_3[k]
                  + f_1 * gs_6[k]
                  + pb_y[k] * gp_5[k];

        t_16[k] = f_1 * gs_6[k]
                  + pb_z[k] * gp_6[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pa_z, pb_x, dd_4, dd_5, fd_9, fd_11, \
                         fd_12, gs_8, gp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = pa_z[k] * fd_9[k];

        t_18[k] = f_2 * dd_5[k]
                  + pa_y[k] * fd_12[k];

        t_19[k] = f_1 * gs_8[k]
                  + pb_x[k] * gp_7[k];

        t_20[k] = f_1 * dd_4[k]
                  + pa_z[k] * fd_11[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pb_x, pb_y, dd_8, fp_6, fd_13, \
                         fd_15, fd_16, gs_10, gp_8, gp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * dd_8[k]
                  + pa_y[k] * fd_13[k];

        t_22[k] = f_2 * fp_6[k]
                  + pa_y[k] * fd_15[k];

        t_23[k] = pa_y[k] * fd_16[k];

        t_24[k] = f_1 * gs_10[k]
                  + pb_x[k] * gp_8[k];

        t_25[k] = f_1 * gs_10[k]
                  + pb_y[k] * gp_9[k];
    }

#pragma omp simd aligned(t_26, t_27, pb_y, pb_z, fp_7, gs_10, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_y[k] * gp_10[k];

        t_27[k] = f_0 * fp_7[k]
                  + f_1 * gs_10[k]
                  + pb_z[k] * gp_10[k];
    }
}

}  // namespace simdovl
