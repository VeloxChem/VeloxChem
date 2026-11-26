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
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_17 = buffer.data(dd + 17);
    const auto *dd_21 = buffer.data(dd + 21);
    const auto *dd_29 = buffer.data(dd + 29);
    const auto *dd_35 = buffer.data(dd + 35);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_17 = buffer.data(fp + 17);
    const auto *fp_18 = buffer.data(fp + 18);
    const auto *fp_19 = buffer.data(fp + 19);
    const auto *fp_23 = buffer.data(fp + 23);
    const auto *fp_25 = buffer.data(fp + 25);
    const auto *fp_26 = buffer.data(fp + 26);
    const auto *fp_27 = buffer.data(fp + 27);
    const auto *fp_28 = buffer.data(fp + 28);
    const auto *fp_29 = buffer.data(fp + 29);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_30 = buffer.data(fd + 30);
    const auto *fd_32 = buffer.data(fd + 32);
    const auto *fd_35 = buffer.data(fd + 35);
    const auto *fd_36 = buffer.data(fd + 36);
    const auto *fd_39 = buffer.data(fd + 39);
    const auto *fd_41 = buffer.data(fd + 41);
    const auto *fd_45 = buffer.data(fd + 45);
    const auto *fd_46 = buffer.data(fd + 46);
    const auto *fd_47 = buffer.data(fd + 47);
    const auto *fd_51 = buffer.data(fd + 51);
    const auto *fd_52 = buffer.data(fd + 52);
    const auto *fd_53 = buffer.data(fd + 53);
    const auto *fd_54 = buffer.data(fd + 54);
    const auto *fd_57 = buffer.data(fd + 57);
    const auto *fd_59 = buffer.data(fd + 59);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_10 = buffer.data(gs + 10);
    const auto *gs_12 = buffer.data(gs + 12);
    const auto *gs_14 = buffer.data(gs + 14);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_23 = buffer.data(gp + 23);
    const auto *gp_25 = buffer.data(gp + 25);
    const auto *gp_27 = buffer.data(gp + 27);
    const auto *gp_29 = buffer.data(gp + 29);
    const auto *gp_30 = buffer.data(gp + 30);
    const auto *gp_31 = buffer.data(gp + 31);
    const auto *gp_32 = buffer.data(gp + 32);
    const auto *gp_34 = buffer.data(gp + 34);
    const auto *gp_35 = buffer.data(gp + 35);
    const auto *gp_36 = buffer.data(gp + 36);
    const auto *gp_37 = buffer.data(gp + 37);
    const auto *gp_38 = buffer.data(gp + 38);
    const auto *gp_40 = buffer.data(gp + 40);
    const auto *gp_41 = buffer.data(gp + 41);
    const auto *gp_42 = buffer.data(gp + 42);
    const auto *gp_43 = buffer.data(gp + 43);
    const auto *gp_44 = buffer.data(gp + 44);

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

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pb_x, pb_z, dd_9, fp_4, fd_0, \
                         fd_9, gp_3, gp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * fd_0[k];

        t_7[k] = f_2 * fp_4[k]
                 + pb_x[k] * gp_4[k];

        t_8[k] = pb_z[k] * gp_3[k];

        t_9[k] = f_3 * dd_9[k]
                 + pa_x[k] * fd_9[k];

        t_10[k] = pb_z[k] * gp_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pa_y, pa_z, pb_x, pb_y, fp_8, \
                         fd_0, fd_3, fd_5, gp_6, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_y[k] * fd_5[k];

        t_12[k] = pa_z[k] * fd_0[k];

        t_13[k] = pb_y[k] * gp_6[k];

        t_14[k] = f_2 * fp_8[k]
                  + pb_x[k] * gp_8[k];

        t_15[k] = pa_z[k] * fd_3[k];

        t_16[k] = pb_y[k] * gp_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_x, pa_y, pb_x, pb_z, dd_0, dd_17, fp_10, \
                         fd_6, fd_17, gp_9, gp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * dd_17[k]
                  + pa_x[k] * fd_17[k];

        t_18[k] = f_1 * dd_0[k]
                  + pa_y[k] * fd_6[k];

        t_19[k] = f_3 * fp_10[k]
                  + pb_x[k] * gp_10[k];

        t_20[k] = pb_z[k] * gp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_x, pa_y, pa_z, pb_z, dd_21, fd_7, \
                         fd_12, fd_21, gs_3, gp_10, gp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * dd_21[k]
                  + pa_x[k] * fd_21[k];

        t_22[k] = pb_z[k] * gp_10[k];

        t_23[k] = f_1 * gs_3[k]
                  + pb_z[k] * gp_11[k];

        t_24[k] = pa_y[k] * fd_12[k];

        t_25[k] = pa_z[k] * fd_7[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, dd_0, fp_8, fd_9, \
                         fd_12, fd_14, fd_17, gp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_y[k] * fd_14[k];

        t_27[k] = pa_z[k] * fd_9[k];

        t_28[k] = f_1 * fp_8[k]
                  + pb_y[k] * gp_14[k];

        t_29[k] = pa_y[k] * fd_17[k];

        t_30[k] = f_1 * dd_0[k]
                  + pa_z[k] * fd_12[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_x, pb_x, pb_y, dd_35, fp_17, fd_35, \
                         gs_5, gp_15, gp_16, gp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * gp_15[k];

        t_32[k] = f_3 * fp_17[k]
                  + pb_x[k] * gp_17[k];

        t_33[k] = f_1 * gs_5[k]
                  + pb_y[k] * gp_16[k];

        t_34[k] = pb_y[k] * gp_17[k];

        t_35[k] = f_1 * dd_35[k]
                  + pa_x[k] * fd_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, pa_x, pb_x, pb_z, fp_18, fp_19, \
                         fd_36, fd_39, fd_41, gp_18, gp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_3 * fp_18[k]
                  + pa_x[k] * fd_36[k];

        t_37[k] = f_1 * fp_19[k]
                  + pb_x[k] * gp_19[k];

        t_38[k] = pb_z[k] * gp_18[k];

        t_39[k] = pa_x[k] * fd_39[k];

        t_40[k] = pb_z[k] * gp_19[k];

        t_41[k] = pa_x[k] * fd_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pa_x, pa_z, pb_x, fp_23, fd_18, \
                         fd_19, fd_45, fd_46, fd_47, gp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * fd_18[k];

        t_43[k] = pa_z[k] * fd_19[k];

        t_44[k] = f_1 * fp_23[k]
                  + pb_x[k] * gp_23[k];

        t_45[k] = pa_x[k] * fd_45[k];

        t_46[k] = pa_x[k] * fd_46[k];

        t_47[k] = pa_x[k] * fd_47[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, pa_x, pa_y, pb_x, fp_25, fd_30, \
                         fd_32, fd_51, fd_52, fd_53, gp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_y[k] * fd_30[k];

        t_49[k] = f_1 * fp_25[k]
                  + pb_x[k] * gp_25[k];

        t_50[k] = pa_y[k] * fd_32[k];

        t_51[k] = pa_x[k] * fd_51[k];

        t_52[k] = pa_x[k] * fd_52[k];

        t_53[k] = pa_x[k] * fd_53[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, pa_x, pb_x, pb_y, fp_27, fp_29, \
                         fd_54, fd_57, fd_59, gp_27, gp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_3 * fp_27[k]
                  + pa_x[k] * fd_54[k];

        t_55[k] = pb_y[k] * gp_27[k];

        t_56[k] = f_1 * fp_29[k]
                  + pb_x[k] * gp_29[k];

        t_57[k] = pa_x[k] * fd_57[k];

        t_58[k] = pb_y[k] * gp_29[k];

        t_59[k] = pa_x[k] * fd_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pb_x, pb_y, pb_z, fp_19, gs_10, \
                         gp_30, gp_31, gp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * gs_10[k]
                  + pb_x[k] * gp_30[k];

        t_61[k] = pb_x[k] * gp_31[k];

        t_62[k] = pb_x[k] * gp_32[k];

        t_63[k] = f_0 * fp_19[k]
                  + f_1 * gs_10[k]
                  + pb_y[k] * gp_31[k];

        t_64[k] = pb_z[k] * gp_31[k];

        t_65[k] = f_1 * gs_10[k]
                  + pb_z[k] * gp_32[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_z, pb_x, pb_y, fp_23, fd_36, fd_39, \
                         gp_34, gp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * fd_36[k];

        t_67[k] = pb_x[k] * gp_34[k];

        t_68[k] = pb_x[k] * gp_35[k];

        t_69[k] = pa_z[k] * fd_39[k];

        t_70[k] = f_2 * fp_23[k]
                  + pb_y[k] * gp_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, pa_y, pa_z, pb_x, dd_21, dd_29, fd_45, \
                         fd_47, gs_12, gp_36, gp_37, gp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_3 * dd_29[k]
                  + pa_y[k] * fd_47[k];

        t_72[k] = f_1 * gs_12[k]
                  + pb_x[k] * gp_36[k];

        t_73[k] = pb_x[k] * gp_37[k];

        t_74[k] = pb_x[k] * gp_38[k];

        t_75[k] = f_1 * dd_21[k]
                  + pa_z[k] * fd_45[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pa_y, pb_x, pb_y, dd_35, fp_26, fd_53, \
                         fd_54, gp_38, gp_40, gp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_3 * fp_26[k]
                  + pb_y[k] * gp_38[k];

        t_77[k] = f_1 * dd_35[k]
                  + pa_y[k] * fd_53[k];

        t_78[k] = pa_y[k] * fd_54[k];

        t_79[k] = pb_x[k] * gp_40[k];

        t_80[k] = pb_x[k] * gp_41[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, pa_y, pb_x, pb_y, fp_28, fp_29, fd_57, \
                         fd_59, gs_14, gp_41, gp_42, gp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_3 * fp_28[k]
                  + pa_y[k] * fd_57[k];

        t_82[k] = f_1 * fp_29[k]
                  + pb_y[k] * gp_41[k];

        t_83[k] = pa_y[k] * fd_59[k];

        t_84[k] = f_1 * gs_14[k]
                  + pb_x[k] * gp_42[k];

        t_85[k] = pb_x[k] * gp_43[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_x, pb_y, pb_z, fp_29, gs_14, gp_43, \
                         gp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_x[k] * gp_44[k];

        t_87[k] = f_1 * gs_14[k]
                  + pb_y[k] * gp_43[k];

        t_88[k] = pb_y[k] * gp_44[k];

        t_89[k] = f_0 * fp_29[k]
                  + f_1 * gs_14[k]
                  + pb_z[k] * gp_44[k];
    }
}

}  // namespace simdovl
