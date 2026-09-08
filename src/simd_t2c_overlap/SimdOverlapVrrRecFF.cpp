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


#include "SimdOverlapVrrRecFF.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_ff_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pf, const size_t dd, const size_t df,
                          const size_t fp, const size_t fd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_16 = buffer.data(dd + 16);
    const auto *dd_18 = buffer.data(dd + 18);

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
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_26 = buffer.data(df + 26);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);
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
    const auto *fd_27 = buffer.data(fd + 27);
    const auto *fd_28 = buffer.data(fd + 28);
    const auto *fd_29 = buffer.data(fd + 29);
    const auto *fd_30 = buffer.data(fd + 30);
    const auto *fd_31 = buffer.data(fd + 31);
    const auto *fd_32 = buffer.data(fd + 32);
    const auto *fd_33 = buffer.data(fd + 33);
    const auto *fd_34 = buffer.data(fd + 34);
    const auto *fd_35 = buffer.data(fd + 35);
    const auto *fd_36 = buffer.data(fd + 36);
    const auto *fd_37 = buffer.data(fd + 37);
    const auto *fd_38 = buffer.data(fd + 38);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dd_0, dd_1, fp_0, fd_0, \
                         fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pb_y[k] * fd_0[k];

        t_2[k] = pb_z[k] * fd_0[k];

        t_3[k] = f_0 * dd_1[k]
                 + pb_x[k] * fd_2[k];

        t_4[k] = pb_y[k] * fd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, t_10, pa_y, pb_x, pb_y, pb_z, dd_2, df_0, \
                         fp_1, fp_2, fd_2, fd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dd_2[k]
                 + pb_x[k] * fd_3[k];

        t_6[k] = f_1 * fp_1[k]
                 + pb_y[k] * fd_2[k];

        t_7[k] = pb_z[k] * fd_2[k];

        t_8[k] = pb_y[k] * fd_3[k];

        t_9[k] = f_1 * fp_2[k]
                 + pb_z[k] * fd_3[k];

        t_10[k] = pa_y[k] * df_0[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_x, pb_y, pb_z, dd_0, dd_4, \
                         df_2, fd_4, fd_5, fd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * dd_0[k]
                  + pb_y[k] * fd_4[k];

        t_12[k] = pb_z[k] * fd_4[k];

        t_13[k] = f_1 * dd_4[k]
                  + pb_x[k] * fd_6[k];

        t_14[k] = pb_z[k] * fd_5[k];

        t_15[k] = pa_y[k] * df_2[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, pf_1, dd_2, df_4, \
                         df_7, fd_6, fd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * pf_1[k]
                  + pa_x[k] * df_7[k];

        t_17[k] = pb_z[k] * fd_6[k];

        t_18[k] = f_2 * dd_2[k]
                  + pb_y[k] * fd_7[k];

        t_19[k] = pa_y[k] * df_4[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, dd_0, df_0, df_1, \
                         fd_8, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * df_0[k];

        t_21[k] = pb_y[k] * fd_8[k];

        t_22[k] = f_2 * dd_0[k]
                  + pb_z[k] * fd_8[k];

        t_23[k] = pa_z[k] * df_1[k];

        t_24[k] = pb_y[k] * fd_9[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_x, pa_z, pb_x, pb_y, pf_2, dd_6, \
                         df_3, df_11, fp_4, fd_10, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * dd_6[k]
                  + pb_x[k] * fd_11[k];

        t_26[k] = pa_z[k] * df_3[k];

        t_27[k] = f_2 * fp_4[k]
                  + pb_y[k] * fd_10[k];

        t_28[k] = pb_y[k] * fd_11[k];

        t_29[k] = f_2 * pf_2[k]
                  + pa_x[k] * df_11[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pa_x, pb_x, pb_y, pb_z, dd_3, dd_7, \
                         dd_9, df_12, fd_12, fd_13, fd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * dd_7[k]
                  + pa_x[k] * df_12[k];

        t_31[k] = f_1 * dd_3[k]
                  + pb_y[k] * fd_12[k];

        t_32[k] = pb_z[k] * fd_12[k];

        t_33[k] = f_2 * dd_9[k]
                  + pb_x[k] * fd_14[k];

        t_34[k] = pb_z[k] * fd_13[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_x, pb_x, pb_z, dd_10, df_14, df_15, \
                         df_16, fd_14, fd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_2 * dd_10[k]
                  + pb_x[k] * fd_15[k];

        t_36[k] = pa_x[k] * df_14[k];

        t_37[k] = pb_z[k] * fd_14[k];

        t_38[k] = pa_x[k] * df_15[k];

        t_39[k] = pa_x[k] * df_16[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, pa_y, pa_z, pb_x, dd_12, df_5, \
                         df_6, df_8, df_9, df_10, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * df_8[k];

        t_41[k] = pa_z[k] * df_5[k];

        t_42[k] = pa_y[k] * df_9[k];

        t_43[k] = pa_z[k] * df_6[k];

        t_44[k] = f_2 * dd_12[k]
                  + pb_x[k] * fd_16[k];

        t_45[k] = pa_y[k] * df_10[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, t_51, pa_x, pb_y, dd_14, df_17, df_18, \
                         df_19, df_20, df_21, fd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_x[k] * df_17[k];

        t_47[k] = pa_x[k] * df_18[k];

        t_48[k] = pa_x[k] * df_19[k];

        t_49[k] = pa_x[k] * df_20[k];

        t_50[k] = f_0 * dd_14[k]
                  + pa_x[k] * df_21[k];

        t_51[k] = pb_y[k] * fd_17[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_x, pb_y, pb_z, dd_5, dd_16, dd_18, fd_17, \
                         fd_18, fd_19, fd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * dd_5[k]
                  + pb_z[k] * fd_17[k];

        t_53[k] = f_2 * dd_16[k]
                  + pb_x[k] * fd_19[k];

        t_54[k] = pb_y[k] * fd_18[k];

        t_55[k] = f_2 * dd_18[k]
                  + pb_x[k] * fd_20[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_x, pb_x, pb_y, df_24, df_25, df_26, \
                         fp_7, fd_20, fd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_x[k] * df_24[k];

        t_57[k] = pa_x[k] * df_25[k];

        t_58[k] = pb_y[k] * fd_20[k];

        t_59[k] = pa_x[k] * df_26[k];

        t_60[k] = f_1 * fp_7[k]
                  + pb_x[k] * fd_21[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, t_66, pb_x, pb_y, pb_z, dd_9, fp_8, \
                         fd_21, fd_22, fd_23, fd_24, fd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_2 * fp_8[k]
                  + pb_x[k] * fd_22[k];

        t_62[k] = pb_z[k] * fd_21[k];

        t_63[k] = pb_x[k] * fd_23[k];

        t_64[k] = pb_x[k] * fd_24[k];

        t_65[k] = pb_x[k] * fd_25[k];

        t_66[k] = f_0 * dd_9[k]
                  + f_1 * fp_8[k]
                  + pb_y[k] * fd_23[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, pa_z, pb_y, pb_z, dd_10, df_12, df_13, \
                         fp_9, fd_23, fd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pb_z[k] * fd_23[k];

        t_68[k] = f_0 * dd_10[k]
                  + pb_y[k] * fd_25[k];

        t_69[k] = f_1 * fp_9[k]
                  + pb_z[k] * fd_25[k];

        t_70[k] = pa_z[k] * df_12[k];

        t_71[k] = pa_z[k] * df_13[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, pa_z, pb_x, pb_z, dd_9, df_14, \
                         fp_10, fd_26, fd_27, fd_28, fd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_2 * fp_10[k]
                  + pb_x[k] * fd_26[k];

        t_73[k] = pb_x[k] * fd_27[k];

        t_74[k] = pb_x[k] * fd_28[k];

        t_75[k] = pb_x[k] * fd_29[k];

        t_76[k] = pa_z[k] * df_14[k];

        t_77[k] = f_2 * dd_9[k]
                  + pb_z[k] * fd_27[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pb_x, pb_y, pf_2, dd_13, df_20, \
                         df_21, df_23, fp_11, fd_29, fd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_1 * dd_13[k]
                  + pb_y[k] * fd_29[k];

        t_79[k] = f_2 * pf_2[k]
                  + pa_y[k] * df_20[k];

        t_80[k] = pa_y[k] * df_21[k];

        t_81[k] = f_2 * fp_11[k]
                  + pb_x[k] * fd_30[k];

        t_82[k] = pa_y[k] * df_23[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, pb_x, pb_z, dd_11, dd_16, df_24, \
                         fd_31, fd_32, fd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pb_x[k] * fd_31[k];

        t_84[k] = pb_x[k] * fd_32[k];

        t_85[k] = pb_x[k] * fd_33[k];

        t_86[k] = f_0 * dd_16[k]
                  + pa_y[k] * df_24[k];

        t_87[k] = f_1 * dd_11[k]
                  + pb_z[k] * fd_31[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_y, pb_x, pb_y, dd_18, df_26, fp_13, \
                         fp_15, fd_33, fd_34, fd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_2 * dd_18[k]
                  + pb_y[k] * fd_33[k];

        t_89[k] = pa_y[k] * df_26[k];

        t_90[k] = f_1 * fp_13[k]
                  + pb_x[k] * fd_34[k];

        t_91[k] = pb_y[k] * fd_34[k];

        t_92[k] = f_2 * fp_15[k]
                  + pb_x[k] * fd_35[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, pb_x, pb_y, pb_z, dd_18, \
                         fp_14, fp_15, fd_36, fd_37, fd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pb_x[k] * fd_36[k];

        t_94[k] = pb_x[k] * fd_37[k];

        t_95[k] = pb_x[k] * fd_38[k];

        t_96[k] = f_1 * fp_14[k]
                  + pb_y[k] * fd_36[k];

        t_97[k] = f_2 * fp_15[k]
                  + pb_y[k] * fd_37[k];

        t_98[k] = pb_y[k] * fd_38[k];

        t_99[k] = f_0 * dd_18[k]
                  + f_1 * fp_15[k]
                  + pb_z[k] * fd_38[k];
    }
}

auto
compute_prim_ff_overlap_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pf, const size_t dd, const size_t df,
                          const size_t fp, const size_t fd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_10 = buffer.data(pf + 10);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_15 = buffer.data(dd + 15);
    const auto *dd_17 = buffer.data(dd + 17);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_29 = buffer.data(df + 29);
    const auto *df_30 = buffer.data(df + 30);
    const auto *df_32 = buffer.data(df + 32);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
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
    const auto *fd_27 = buffer.data(fd + 27);
    const auto *fd_28 = buffer.data(fd + 28);
    const auto *fd_29 = buffer.data(fd + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dd_0, dd_1, dd_2, fp_0, \
                         fd_0, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pb_y[k] * fd_0[k];

        t_2[k] = pb_z[k] * fd_0[k];

        t_3[k] = f_0 * dd_1[k]
                 + pb_x[k] * fd_1[k];

        t_4[k] = f_0 * dd_2[k]
                 + pb_x[k] * fd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, dd_0, df_0, fp_1, fp_2, \
                         fd_1, fd_2, fd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_6[k] = pb_y[k] * fd_2[k];

        t_7[k] = f_1 * fp_2[k]
                 + pb_z[k] * fd_2[k];

        t_8[k] = pa_y[k] * df_0[k];

        t_9[k] = f_2 * dd_0[k]
                 + pb_y[k] * fd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_y, pb_z, pf_4, dd_2, dd_4, \
                         df_6, fd_4, fd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * dd_4[k]
                  + pb_x[k] * fd_4[k];

        t_11[k] = f_2 * pf_4[k]
                  + pa_x[k] * df_6[k];

        t_12[k] = pb_z[k] * fd_4[k];

        t_13[k] = f_2 * dd_2[k]
                  + pb_y[k] * fd_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pa_z, pb_x, pb_z, dd_0, dd_6, df_0, \
                         df_4, fd_6, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_y[k] * df_4[k];

        t_15[k] = pa_z[k] * df_0[k];

        t_16[k] = f_2 * dd_0[k]
                  + pb_z[k] * fd_6[k];

        t_17[k] = f_1 * dd_6[k]
                  + pb_x[k] * fd_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_x, pb_y, pf_10, dd_3, dd_7, df_11, \
                         df_12, fp_3, fd_7, fd_8, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * fp_3[k]
                  + pb_y[k] * fd_7[k];

        t_19[k] = pb_y[k] * fd_8[k];

        t_20[k] = f_2 * pf_10[k]
                  + pa_x[k] * df_11[k];

        t_21[k] = f_0 * dd_7[k]
                  + pa_x[k] * df_12[k];

        t_22[k] = f_1 * dd_3[k]
                  + pb_y[k] * fd_9[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pa_x, pb_x, pb_z, dd_9, df_16, df_18, \
                         df_19, fd_9, fd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pb_z[k] * fd_9[k];

        t_24[k] = f_2 * dd_9[k]
                  + pb_x[k] * fd_10[k];

        t_25[k] = pa_x[k] * df_16[k];

        t_26[k] = pa_x[k] * df_18[k];

        t_27[k] = pa_x[k] * df_19[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_x, pa_y, pb_y, dd_13, df_9, df_21, \
                         df_22, df_24, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_y[k] * df_9[k];

        t_29[k] = pa_x[k] * df_21[k];

        t_30[k] = pa_x[k] * df_22[k];

        t_31[k] = f_0 * dd_13[k]
                  + pa_x[k] * df_24[k];

        t_32[k] = pb_y[k] * fd_11[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, dd_5, dd_17, df_29, \
                         df_30, df_32, fd_11, fd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_1 * dd_5[k]
                  + pb_z[k] * fd_11[k];

        t_34[k] = f_2 * dd_17[k]
                  + pb_x[k] * fd_12[k];

        t_35[k] = pa_x[k] * df_29[k];

        t_36[k] = pa_x[k] * df_30[k];

        t_37[k] = pa_x[k] * df_32[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, pb_x, pb_y, pb_z, dd_9, fp_4, \
                         fp_5, fd_13, fd_14, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_1 * fp_4[k]
                  + pb_x[k] * fd_13[k];

        t_39[k] = f_2 * fp_5[k]
                  + pb_x[k] * fd_14[k];

        t_40[k] = pb_x[k] * fd_15[k];

        t_41[k] = pb_x[k] * fd_16[k];

        t_42[k] = f_0 * dd_9[k]
                  + f_1 * fp_5[k]
                  + pb_y[k] * fd_15[k];

        t_43[k] = pb_z[k] * fd_15[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pb_x, pb_y, pb_z, dd_10, fp_6, fp_7, \
                         fd_16, fd_17, fd_19, fd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * dd_10[k]
                  + pb_y[k] * fd_16[k];

        t_45[k] = f_1 * fp_6[k]
                  + pb_z[k] * fd_16[k];

        t_46[k] = f_2 * fp_7[k]
                  + pb_x[k] * fd_17[k];

        t_47[k] = pb_x[k] * fd_19[k];

        t_48[k] = pb_x[k] * fd_20[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_y, pa_z, pb_y, pb_z, pf_10, dd_9, dd_12, \
                         df_16, df_23, fd_18, fd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pa_z[k] * df_16[k];

        t_50[k] = f_2 * dd_9[k]
                  + pb_z[k] * fd_18[k];

        t_51[k] = f_1 * dd_12[k]
                  + pb_y[k] * fd_20[k];

        t_52[k] = f_2 * pf_10[k]
                  + pa_y[k] * df_23[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_y, pb_x, pb_z, dd_11, dd_15, df_29, \
                         fp_8, fd_21, fd_22, fd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_2 * fp_8[k]
                  + pb_x[k] * fd_21[k];

        t_54[k] = pb_x[k] * fd_22[k];

        t_55[k] = pb_x[k] * fd_23[k];

        t_56[k] = f_0 * dd_15[k]
                  + pa_y[k] * df_29[k];

        t_57[k] = f_1 * dd_11[k]
                  + pb_z[k] * fd_22[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_y, pb_x, pb_y, dd_17, df_32, fp_10, \
                         fp_12, fd_24, fd_25, fd_26, fd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_2 * dd_17[k]
                  + pb_y[k] * fd_24[k];

        t_59[k] = pa_y[k] * df_32[k];

        t_60[k] = f_1 * fp_10[k]
                  + pb_x[k] * fd_25[k];

        t_61[k] = f_2 * fp_12[k]
                  + pb_x[k] * fd_26[k];

        t_62[k] = pb_x[k] * fd_27[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pb_x, pb_y, pb_z, dd_17, fp_11, fp_12, \
                         fd_27, fd_28, fd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_x[k] * fd_29[k];

        t_64[k] = f_1 * fp_11[k]
                  + pb_y[k] * fd_27[k];

        t_65[k] = f_2 * fp_12[k]
                  + pb_y[k] * fd_28[k];

        t_66[k] = pb_y[k] * fd_29[k];

        t_67[k] = f_0 * dd_17[k]
                  + f_1 * fp_12[k]
                  + pb_z[k] * fd_29[k];
    }
}

auto
compute_prim_ff_overlap_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pf, const size_t dd, const size_t df,
                          const size_t fp, const size_t fd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_9 = buffer.data(pf + 9);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_16 = buffer.data(dd + 16);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_27 = buffer.data(df + 27);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, dd_0, fp_0, fp_1, \
                         fp_2, fd_0, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pb_y[k] * fd_0[k];

        t_2[k] = pb_z[k] * fd_0[k];

        t_3[k] = f_1 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = pb_y[k] * fd_2[k];

        t_5[k] = f_1 * fp_2[k]
                 + pb_z[k] * fd_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, pb_z, pf_3, df_0, df_5, \
                         df_7, fd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * df_0[k];

        t_7[k] = f_2 * pf_3[k]
                 + pa_x[k] * df_7[k];

        t_8[k] = pb_z[k] * fd_3[k];

        t_9[k] = pa_y[k] * df_5[k];

        t_10[k] = pa_z[k] * df_0[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pa_x, pb_y, pb_z, pf_9, dd_0, df_9, fp_3, \
                         fd_4, fd_5, fd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * dd_0[k]
                  + pb_z[k] * fd_4[k];

        t_12[k] = f_2 * fp_3[k]
                  + pb_y[k] * fd_5[k];

        t_13[k] = pb_y[k] * fd_6[k];

        t_14[k] = f_2 * pf_9[k]
                  + pa_x[k] * df_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_x, pb_y, pb_z, dd_6, dd_12, df_10, \
                         df_14, df_20, fd_7, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * dd_6[k]
                  + pa_x[k] * df_10[k];

        t_16[k] = pb_z[k] * fd_7[k];

        t_17[k] = pa_x[k] * df_14[k];

        t_18[k] = f_0 * dd_12[k]
                  + pa_x[k] * df_20[k];

        t_19[k] = pb_y[k] * fd_8[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_x, pb_x, pb_z, dd_4, df_27, fp_4, \
                         fp_5, fd_8, fd_9, fd_10, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * dd_4[k]
                  + pb_z[k] * fd_8[k];

        t_21[k] = pa_x[k] * df_27[k];

        t_22[k] = f_1 * fp_4[k]
                  + pb_x[k] * fd_9[k];

        t_23[k] = f_2 * fp_5[k]
                  + pb_x[k] * fd_10[k];

        t_24[k] = pb_x[k] * fd_11[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, dd_8, dd_9, fp_5, \
                         fp_6, fd_11, fd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pb_x[k] * fd_12[k];

        t_26[k] = f_0 * dd_8[k]
                  + f_1 * fp_5[k]
                  + pb_y[k] * fd_11[k];

        t_27[k] = pb_z[k] * fd_11[k];

        t_28[k] = f_0 * dd_9[k]
                  + pb_y[k] * fd_12[k];

        t_29[k] = f_1 * fp_6[k]
                  + pb_z[k] * fd_12[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pa_z, pb_x, pb_z, dd_8, df_14, fp_7, \
                         fd_13, fd_14, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * fp_7[k]
                  + pb_x[k] * fd_13[k];

        t_31[k] = pb_x[k] * fd_15[k];

        t_32[k] = pb_x[k] * fd_16[k];

        t_33[k] = pa_z[k] * df_14[k];

        t_34[k] = f_2 * dd_8[k]
                  + pb_z[k] * fd_14[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_y, pb_x, pb_y, pf_9, dd_11, df_19, \
                         fp_8, fd_16, fd_17, fd_18, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * dd_11[k]
                  + pb_y[k] * fd_16[k];

        t_36[k] = f_2 * pf_9[k]
                  + pa_y[k] * df_19[k];

        t_37[k] = f_2 * fp_8[k]
                  + pb_x[k] * fd_17[k];

        t_38[k] = pb_x[k] * fd_18[k];

        t_39[k] = pb_x[k] * fd_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pb_y, pb_z, dd_10, dd_14, dd_16, df_24, \
                         df_27, fd_18, fd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * dd_14[k]
                  + pa_y[k] * df_24[k];

        t_41[k] = f_1 * dd_10[k]
                  + pb_z[k] * fd_18[k];

        t_42[k] = f_2 * dd_16[k]
                  + pb_y[k] * fd_20[k];

        t_43[k] = pa_y[k] * df_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, pb_x, pb_y, fp_10, fp_11, fp_12, \
                         fd_21, fd_22, fd_23, fd_24, fd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_1 * fp_10[k]
                  + pb_x[k] * fd_21[k];

        t_45[k] = f_2 * fp_12[k]
                  + pb_x[k] * fd_22[k];

        t_46[k] = pb_x[k] * fd_23[k];

        t_47[k] = pb_x[k] * fd_25[k];

        t_48[k] = f_1 * fp_11[k]
                  + pb_y[k] * fd_23[k];

        t_49[k] = f_2 * fp_12[k]
                  + pb_y[k] * fd_24[k];
    }

#pragma omp simd aligned(t_50, t_51, pb_y, pb_z, dd_16, fp_12, fd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pb_y[k] * fd_25[k];

        t_51[k] = f_0 * dd_16[k]
                  + f_1 * fp_12[k]
                  + pb_z[k] * fd_25[k];
    }
}

auto
compute_prim_ff_overlap_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pf, const size_t dd, const size_t df,
                          const size_t fp, const size_t fd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_9 = buffer.data(pf + 9);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_16 = buffer.data(dd + 16);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_26 = buffer.data(df + 26);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, dd_0, fp_0, fp_1, \
                         fp_2, fd_0, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pb_y[k] * fd_0[k];

        t_2[k] = pb_z[k] * fd_0[k];

        t_3[k] = f_1 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = pb_y[k] * fd_2[k];

        t_5[k] = f_1 * fp_2[k]
                 + pb_z[k] * fd_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, pa_z, pb_z, pf_3, dd_0, df_0, df_6, fd_3, \
                         fd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * pf_3[k]
                 + pa_x[k] * df_6[k];

        t_7[k] = pb_z[k] * fd_3[k];

        t_8[k] = pa_z[k] * df_0[k];

        t_9[k] = f_2 * dd_0[k]
                 + pb_z[k] * fd_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pb_y, pb_z, pf_9, dd_6, df_8, \
                         df_9, fp_3, fd_5, fd_6, fd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * fp_3[k]
                  + pb_y[k] * fd_5[k];

        t_11[k] = pb_y[k] * fd_6[k];

        t_12[k] = f_2 * pf_9[k]
                  + pa_x[k] * df_8[k];

        t_13[k] = f_0 * dd_6[k]
                  + pa_x[k] * df_9[k];

        t_14[k] = pb_z[k] * fd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pb_x, pb_y, pb_z, dd_4, dd_12, df_19, \
                         fp_4, fd_8, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * dd_12[k]
                  + pa_x[k] * df_19[k];

        t_16[k] = pb_y[k] * fd_8[k];

        t_17[k] = f_1 * dd_4[k]
                  + pb_z[k] * fd_8[k];

        t_18[k] = f_1 * fp_4[k]
                  + pb_x[k] * fd_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, pb_x, pb_y, pb_z, dd_8, dd_9, \
                         fp_5, fd_10, fd_11, fd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * fp_5[k]
                  + pb_x[k] * fd_10[k];

        t_20[k] = pb_x[k] * fd_11[k];

        t_21[k] = pb_x[k] * fd_12[k];

        t_22[k] = f_0 * dd_8[k]
                  + f_1 * fp_5[k]
                  + pb_y[k] * fd_11[k];

        t_23[k] = pb_z[k] * fd_11[k];

        t_24[k] = f_0 * dd_9[k]
                  + pb_y[k] * fd_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_z, pb_x, pb_z, df_13, fp_6, fp_7, \
                         fd_12, fd_13, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * fp_6[k]
                  + pb_z[k] * fd_12[k];

        t_26[k] = f_2 * fp_7[k]
                  + pb_x[k] * fd_13[k];

        t_27[k] = pb_x[k] * fd_15[k];

        t_28[k] = pb_x[k] * fd_16[k];

        t_29[k] = pa_z[k] * df_13[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_x, pb_y, pb_z, pf_9, dd_8, dd_11, \
                         df_18, fp_8, fd_14, fd_16, fd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * dd_8[k]
                  + pb_z[k] * fd_14[k];

        t_31[k] = f_1 * dd_11[k]
                  + pb_y[k] * fd_16[k];

        t_32[k] = f_2 * pf_9[k]
                  + pa_y[k] * df_18[k];

        t_33[k] = f_2 * fp_8[k]
                  + pb_x[k] * fd_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_y, pb_x, pb_y, pb_z, dd_10, dd_14, \
                         dd_16, df_23, fd_18, fd_19, fd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_x[k] * fd_18[k];

        t_35[k] = pb_x[k] * fd_19[k];

        t_36[k] = f_0 * dd_14[k]
                  + pa_y[k] * df_23[k];

        t_37[k] = f_1 * dd_10[k]
                  + pb_z[k] * fd_18[k];

        t_38[k] = f_2 * dd_16[k]
                  + pb_y[k] * fd_20[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_y, pb_x, df_26, fp_10, fp_12, fd_21, \
                         fd_22, fd_23, fd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pa_y[k] * df_26[k];

        t_40[k] = f_1 * fp_10[k]
                  + pb_x[k] * fd_21[k];

        t_41[k] = f_2 * fp_12[k]
                  + pb_x[k] * fd_22[k];

        t_42[k] = pb_x[k] * fd_23[k];

        t_43[k] = pb_x[k] * fd_25[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_y, pb_z, dd_16, fp_11, fp_12, fd_23, \
                         fd_24, fd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_1 * fp_11[k]
                  + pb_y[k] * fd_23[k];

        t_45[k] = f_2 * fp_12[k]
                  + pb_y[k] * fd_24[k];

        t_46[k] = pb_y[k] * fd_25[k];

        t_47[k] = f_0 * dd_16[k]
                  + f_1 * fp_12[k]
                  + pb_z[k] * fd_25[k];
    }
}

auto
compute_prim_ff_overlap_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pf, const size_t dd, const size_t df,
                          const size_t fp, const size_t fd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_15 = buffer.data(dd + 15);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
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
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_22 = buffer.data(fd + 22);
    const auto *fd_24 = buffer.data(fd + 24);
    const auto *fd_25 = buffer.data(fd + 25);
    const auto *fd_26 = buffer.data(fd + 26);
    const auto *fd_27 = buffer.data(fd + 27);
    const auto *fd_28 = buffer.data(fd + 28);
    const auto *fd_29 = buffer.data(fd + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, dd_0, dd_1, dd_2, fp_0, fp_1, fd_0, \
                         fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = f_0 * dd_1[k]
                 + pb_x[k] * fd_1[k];

        t_2[k] = f_0 * dd_2[k]
                 + pb_x[k] * fd_2[k];

        t_3[k] = f_1 * fp_1[k]
                 + pb_y[k] * fd_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pb_x, pb_y, pb_z, dd_0, dd_4, df_0, fp_2, \
                         fd_2, fd_3, fd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * fp_2[k]
                 + pb_z[k] * fd_2[k];

        t_5[k] = pa_y[k] * df_0[k];

        t_6[k] = f_2 * dd_0[k]
                 + pb_y[k] * fd_3[k];

        t_7[k] = f_1 * dd_4[k]
                 + pb_x[k] * fd_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_z, pb_x, pb_z, pf_0, dd_0, dd_6, df_0, \
                         df_1, fd_6, fd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * pf_0[k]
                 + pa_x[k] * df_1[k];

        t_9[k] = pa_z[k] * df_0[k];

        t_10[k] = f_2 * dd_0[k]
                  + pb_z[k] * fd_6[k];

        t_11[k] = f_1 * dd_6[k]
                  + pb_x[k] * fd_7[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, pb_y, pf_1, dd_3, dd_7, dd_8, \
                         df_2, df_3, fd_8, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_2 * pf_1[k]
                  + pa_x[k] * df_2[k];

        t_13[k] = f_0 * dd_7[k]
                  + pa_x[k] * df_3[k];

        t_14[k] = f_1 * dd_3[k]
                  + pb_y[k] * fd_8[k];

        t_15[k] = f_2 * dd_8[k]
                  + pb_x[k] * fd_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_x, pb_x, pb_z, dd_5, dd_13, dd_15, \
                         df_4, df_6, df_8, fd_12, fd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_x[k] * df_4[k];

        t_17[k] = f_0 * dd_13[k]
                  + pa_x[k] * df_6[k];

        t_18[k] = f_1 * dd_5[k]
                  + pb_z[k] * fd_12[k];

        t_19[k] = f_2 * dd_15[k]
                  + pb_x[k] * fd_14[k];

        t_20[k] = pa_x[k] * df_8[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_x, pb_y, dd_8, dd_9, fp_7, fp_8, fd_15, \
                         fd_16, fd_17, fd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * fp_7[k]
                  + pb_x[k] * fd_15[k];

        t_22[k] = f_2 * fp_8[k]
                  + pb_x[k] * fd_16[k];

        t_23[k] = f_0 * dd_8[k]
                  + f_1 * fp_8[k]
                  + pb_y[k] * fd_17[k];

        t_24[k] = f_0 * dd_9[k]
                  + pb_y[k] * fd_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_z, pb_y, pb_z, dd_8, dd_12, df_4, fp_9, \
                         fd_18, fd_19, fd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * fp_9[k]
                  + pb_z[k] * fd_18[k];

        t_26[k] = pa_z[k] * df_4[k];

        t_27[k] = f_2 * dd_8[k]
                  + pb_z[k] * fd_19[k];

        t_28[k] = f_1 * dd_12[k]
                  + pb_y[k] * fd_21[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pb_y, pb_z, pf_1, dd_10, dd_14, dd_15, \
                         df_5, df_7, fd_22, fd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_2 * pf_1[k]
                  + pa_y[k] * df_5[k];

        t_30[k] = f_0 * dd_14[k]
                  + pa_y[k] * df_7[k];

        t_31[k] = f_1 * dd_10[k]
                  + pb_z[k] * fd_22[k];

        t_32[k] = f_2 * dd_15[k]
                  + pb_y[k] * fd_24[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pb_x, pb_y, df_8, fp_12, fp_13, \
                         fp_14, fd_25, fd_26, fd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_y[k] * df_8[k];

        t_34[k] = f_1 * fp_12[k]
                  + pb_x[k] * fd_25[k];

        t_35[k] = pb_y[k] * fd_25[k];

        t_36[k] = f_2 * fp_14[k]
                  + pb_x[k] * fd_26[k];

        t_37[k] = f_1 * fp_13[k]
                  + pb_y[k] * fd_27[k];
    }

#pragma omp simd aligned(t_38, t_39, pb_y, pb_z, dd_15, fp_14, fd_28, \
                         fd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_2 * fp_14[k]
                  + pb_y[k] * fd_28[k];

        t_39[k] = f_0 * dd_15[k]
                  + f_1 * fp_14[k]
                  + pb_z[k] * fd_29[k];
    }
}

auto
compute_prim_ff_overlap_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pf, const size_t dd, const size_t df,
                          const size_t fp, const size_t fd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_3 = buffer.data(pf + 3);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_14 = buffer.data(dd + 14);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_27 = buffer.data(df + 27);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dd_0, fp_0, fp_1, fp_2, \
                         fd_0, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pb_y[k] * fd_0[k];

        t_2[k] = pb_z[k] * fd_0[k];

        t_3[k] = f_1 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = f_1 * fp_2[k]
                 + pb_z[k] * fd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pa_z, pb_y, pf_1, dd_2, df_0, \
                         df_4, df_6, fd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_y[k] * df_0[k];

        t_6[k] = f_2 * pf_1[k]
                 + pa_x[k] * df_6[k];

        t_7[k] = f_2 * dd_2[k]
                 + pb_y[k] * fd_5[k];

        t_8[k] = pa_y[k] * df_4[k];

        t_9[k] = pa_z[k] * df_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_y, pb_z, pf_3, dd_0, dd_6, df_9, \
                         df_10, fp_3, fd_6, fd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * dd_0[k]
                  + pb_z[k] * fd_6[k];

        t_11[k] = f_2 * fp_3[k]
                  + pb_y[k] * fd_7[k];

        t_12[k] = f_2 * pf_3[k]
                  + pa_x[k] * df_9[k];

        t_13[k] = f_0 * dd_6[k]
                  + pa_x[k] * df_10[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, t_19, pa_x, pa_y, pb_x, dd_7, df_8, \
                         df_12, df_14, df_15, df_17, fd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_2 * dd_7[k]
                  + pb_x[k] * fd_10[k];

        t_15[k] = pa_x[k] * df_12[k];

        t_16[k] = pa_x[k] * df_14[k];

        t_17[k] = pa_x[k] * df_15[k];

        t_18[k] = pa_y[k] * df_8[k];

        t_19[k] = pa_x[k] * df_17[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_x, pb_x, pb_z, dd_4, dd_11, dd_14, \
                         df_18, df_20, df_24, fd_11, fd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_x[k] * df_18[k];

        t_21[k] = f_0 * dd_11[k]
                  + pa_x[k] * df_20[k];

        t_22[k] = f_1 * dd_4[k]
                  + pb_z[k] * fd_11[k];

        t_23[k] = f_2 * dd_14[k]
                  + pb_x[k] * fd_12[k];

        t_24[k] = pa_x[k] * df_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, t_30, pa_x, pb_x, df_25, df_27, fp_4, \
                         fp_5, fd_13, fd_14, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_x[k] * df_25[k];

        t_26[k] = pa_x[k] * df_27[k];

        t_27[k] = f_1 * fp_4[k]
                  + pb_x[k] * fd_13[k];

        t_28[k] = f_2 * fp_5[k]
                  + pb_x[k] * fd_14[k];

        t_29[k] = pb_x[k] * fd_15[k];

        t_30[k] = pb_x[k] * fd_16[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pb_x, pb_y, pb_z, dd_7, dd_8, fp_5, \
                         fp_6, fd_15, fd_16, fd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * dd_7[k]
                  + f_1 * fp_5[k]
                  + pb_y[k] * fd_15[k];

        t_32[k] = pb_z[k] * fd_15[k];

        t_33[k] = f_0 * dd_8[k]
                  + pb_y[k] * fd_16[k];

        t_34[k] = f_1 * fp_6[k]
                  + pb_z[k] * fd_16[k];

        t_35[k] = pb_x[k] * fd_18[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_y, pa_z, pb_y, pb_z, pf_3, dd_7, dd_10, \
                         df_12, df_19, fd_17, fd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_z[k] * df_12[k];

        t_37[k] = f_2 * dd_7[k]
                  + pb_z[k] * fd_17[k];

        t_38[k] = f_1 * dd_10[k]
                  + pb_y[k] * fd_18[k];

        t_39[k] = f_2 * pf_3[k]
                  + pa_y[k] * df_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pa_y, pb_x, pb_y, pb_z, dd_9, dd_12, \
                         dd_14, df_24, df_27, fd_19, fd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_x[k] * fd_19[k];

        t_41[k] = f_0 * dd_12[k]
                  + pa_y[k] * df_24[k];

        t_42[k] = f_1 * dd_9[k]
                  + pb_z[k] * fd_19[k];

        t_43[k] = f_2 * dd_14[k]
                  + pb_y[k] * fd_20[k];

        t_44[k] = pa_y[k] * df_27[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, pb_x, pb_y, fp_9, fp_10, fp_11, \
                         fd_21, fd_22, fd_23, fd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_1 * fp_9[k]
                  + pb_x[k] * fd_21[k];

        t_46[k] = pb_y[k] * fd_21[k];

        t_47[k] = f_2 * fp_11[k]
                  + pb_x[k] * fd_22[k];

        t_48[k] = pb_x[k] * fd_23[k];

        t_49[k] = pb_x[k] * fd_25[k];

        t_50[k] = f_1 * fp_10[k]
                  + pb_y[k] * fd_23[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_y, pb_z, dd_14, fp_11, fd_24, \
                         fd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_2 * fp_11[k]
                  + pb_y[k] * fd_24[k];

        t_52[k] = pb_y[k] * fd_25[k];

        t_53[k] = f_0 * dd_14[k]
                  + f_1 * fp_11[k]
                  + pb_z[k] * fd_25[k];
    }
}

auto
compute_prim_ff_overlap_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pf, const size_t dd, const size_t df,
                          const size_t fp, const size_t fd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_6 = buffer.data(pf + 6);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_12 = buffer.data(dd + 12);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_23 = buffer.data(df + 23);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, dd_0, fp_0, fp_1, \
                         fp_2, fd_0, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pb_y[k] * fd_0[k];

        t_2[k] = pb_z[k] * fd_0[k];

        t_3[k] = f_1 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = pb_y[k] * fd_2[k];

        t_5[k] = f_1 * fp_2[k]
                 + pb_z[k] * fd_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, pb_z, pf_2, df_0, df_4, \
                         df_5, fd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * df_0[k];

        t_7[k] = f_2 * pf_2[k]
                 + pa_x[k] * df_5[k];

        t_8[k] = pb_z[k] * fd_4[k];

        t_9[k] = pa_y[k] * df_4[k];

        t_10[k] = pa_z[k] * df_0[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_x, pb_y, pb_z, pf_6, dd_5, df_7, \
                         df_8, fp_3, fd_6, fd_7, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * fp_3[k]
                  + pb_y[k] * fd_6[k];

        t_12[k] = pb_y[k] * fd_7[k];

        t_13[k] = f_2 * pf_6[k]
                  + pa_x[k] * df_7[k];

        t_14[k] = f_0 * dd_5[k]
                  + pa_x[k] * df_8[k];

        t_15[k] = pb_z[k] * fd_8[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_x, pb_x, dd_6, dd_9, dd_12, df_11, \
                         df_16, df_23, fd_9, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * dd_6[k]
                  + pb_x[k] * fd_9[k];

        t_17[k] = pa_x[k] * df_11[k];

        t_18[k] = f_0 * dd_9[k]
                  + pa_x[k] * df_16[k];

        t_19[k] = f_2 * dd_12[k]
                  + pb_x[k] * fd_11[k];

        t_20[k] = pa_x[k] * df_23[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, t_26, pb_x, pb_y, pb_z, dd_6, fp_4, \
                         fp_5, fd_12, fd_13, fd_14, fd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * fp_4[k]
                  + pb_x[k] * fd_12[k];

        t_22[k] = f_2 * fp_5[k]
                  + pb_x[k] * fd_13[k];

        t_23[k] = pb_x[k] * fd_14[k];

        t_24[k] = pb_x[k] * fd_15[k];

        t_25[k] = f_0 * dd_6[k]
                  + f_1 * fp_5[k]
                  + pb_y[k] * fd_14[k];

        t_26[k] = pb_z[k] * fd_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_z, pb_x, pb_y, pb_z, dd_7, df_11, fp_6, \
                         fd_15, fd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * dd_7[k]
                  + pb_y[k] * fd_15[k];

        t_28[k] = f_1 * fp_6[k]
                  + pb_z[k] * fd_15[k];

        t_29[k] = pb_x[k] * fd_17[k];

        t_30[k] = pa_z[k] * df_11[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_y, pb_x, pb_y, pf_6, dd_10, dd_12, \
                         df_15, df_20, df_23, fd_18, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_2 * pf_6[k]
                  + pa_y[k] * df_15[k];

        t_32[k] = pb_x[k] * fd_18[k];

        t_33[k] = f_0 * dd_10[k]
                  + pa_y[k] * df_20[k];

        t_34[k] = f_2 * dd_12[k]
                  + pb_y[k] * fd_19[k];

        t_35[k] = pa_y[k] * df_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, pb_x, pb_y, fp_9, fp_10, fp_11, \
                         fd_20, fd_21, fd_22, fd_23, fd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * fp_9[k]
                  + pb_x[k] * fd_20[k];

        t_37[k] = f_2 * fp_11[k]
                  + pb_x[k] * fd_21[k];

        t_38[k] = pb_x[k] * fd_22[k];

        t_39[k] = pb_x[k] * fd_24[k];

        t_40[k] = f_1 * fp_10[k]
                  + pb_y[k] * fd_22[k];

        t_41[k] = f_2 * fp_11[k]
                  + pb_y[k] * fd_23[k];
    }

#pragma omp simd aligned(t_42, t_43, pb_y, pb_z, dd_12, fp_11, fd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * fd_24[k];

        t_43[k] = f_0 * dd_12[k]
                  + f_1 * fp_11[k]
                  + pb_z[k] * fd_24[k];
    }
}

auto
compute_prim_ff_overlap_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pf, const size_t dd, const size_t df,
                          const size_t fp, const size_t fd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_5 = buffer.data(pf + 5);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_12 = buffer.data(dd + 12);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_21 = buffer.data(df + 21);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, dd_0, fp_0, fp_1, \
                         fp_2, fd_0, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pb_y[k] * fd_0[k];

        t_2[k] = pb_z[k] * fd_0[k];

        t_3[k] = f_1 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = pb_y[k] * fd_2[k];

        t_5[k] = f_1 * fp_2[k]
                 + pb_z[k] * fd_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, pb_y, pb_z, pf_2, df_0, \
                         df_5, fp_3, fd_4, fd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * df_0[k];

        t_7[k] = f_2 * pf_2[k]
                 + pa_x[k] * df_5[k];

        t_8[k] = pb_z[k] * fd_4[k];

        t_9[k] = pa_z[k] * df_0[k];

        t_10[k] = f_2 * fp_3[k]
                  + pb_y[k] * fd_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pa_x, pb_y, pb_z, pf_5, dd_5, df_6, df_7, \
                         fd_7, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_y[k] * fd_7[k];

        t_12[k] = f_2 * pf_5[k]
                  + pa_x[k] * df_6[k];

        t_13[k] = f_0 * dd_5[k]
                  + pa_x[k] * df_7[k];

        t_14[k] = pb_z[k] * fd_8[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_x, pb_x, dd_6, dd_9, dd_12, df_10, \
                         df_14, df_21, fd_9, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_2 * dd_6[k]
                  + pb_x[k] * fd_9[k];

        t_16[k] = pa_x[k] * df_10[k];

        t_17[k] = f_0 * dd_9[k]
                  + pa_x[k] * df_14[k];

        t_18[k] = f_2 * dd_12[k]
                  + pb_x[k] * fd_11[k];

        t_19[k] = pa_x[k] * df_21[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, pb_x, pb_y, pb_z, dd_6, fp_4, \
                         fp_5, fd_12, fd_13, fd_14, fd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * fp_4[k]
                  + pb_x[k] * fd_12[k];

        t_21[k] = f_2 * fp_5[k]
                  + pb_x[k] * fd_13[k];

        t_22[k] = pb_x[k] * fd_14[k];

        t_23[k] = pb_x[k] * fd_15[k];

        t_24[k] = f_0 * dd_6[k]
                  + f_1 * fp_5[k]
                  + pb_y[k] * fd_14[k];

        t_25[k] = pb_z[k] * fd_14[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, pb_z, dd_7, df_10, fp_6, \
                         fd_15, fd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * dd_7[k]
                  + pb_y[k] * fd_15[k];

        t_27[k] = f_1 * fp_6[k]
                  + pb_z[k] * fd_15[k];

        t_28[k] = pb_x[k] * fd_17[k];

        t_29[k] = pa_z[k] * df_10[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pa_y, pb_x, pb_y, pf_5, dd_10, dd_12, \
                         df_13, df_18, df_21, fd_18, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * pf_5[k]
                  + pa_y[k] * df_13[k];

        t_31[k] = pb_x[k] * fd_18[k];

        t_32[k] = f_0 * dd_10[k]
                  + pa_y[k] * df_18[k];

        t_33[k] = f_2 * dd_12[k]
                  + pb_y[k] * fd_19[k];

        t_34[k] = pa_y[k] * df_21[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, t_40, pb_x, pb_y, fp_9, fp_10, fp_11, \
                         fd_20, fd_21, fd_22, fd_23, fd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * fp_9[k]
                  + pb_x[k] * fd_20[k];

        t_36[k] = f_2 * fp_11[k]
                  + pb_x[k] * fd_21[k];

        t_37[k] = pb_x[k] * fd_22[k];

        t_38[k] = pb_x[k] * fd_24[k];

        t_39[k] = f_1 * fp_10[k]
                  + pb_y[k] * fd_22[k];

        t_40[k] = f_2 * fp_11[k]
                  + pb_y[k] * fd_23[k];
    }

#pragma omp simd aligned(t_41, t_42, pb_y, pb_z, dd_12, fp_11, fd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pb_y[k] * fd_24[k];

        t_42[k] = f_0 * dd_12[k]
                  + f_1 * fp_11[k]
                  + pb_z[k] * fd_24[k];
    }
}

auto
compute_prim_ff_overlap_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pf, const size_t dd, const size_t df,
                          const size_t fp, const size_t fd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_22 = buffer.data(fd + 22);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, pf_0, dd_0, df_0, df_1, \
                         fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pa_y[k] * df_0[k];

        t_2[k] = f_2 * pf_0[k]
                 + pa_x[k] * df_1[k];

        t_3[k] = pa_z[k] * df_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pb_x, pf_1, df_2, df_3, df_6, fp_4, \
                         fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * pf_1[k]
                 + pa_x[k] * df_2[k];

        t_5[k] = pa_x[k] * df_3[k];

        t_6[k] = pa_x[k] * df_6[k];

        t_7[k] = f_1 * fp_4[k]
                 + pb_x[k] * fd_11[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_y, pa_z, pb_y, pf_1, dd_4, dd_7, df_3, df_4, \
                         df_5, fp_5, fd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dd_4[k]
                 + f_1 * fp_5[k]
                 + pb_y[k] * fd_12[k];

        t_9[k] = pa_z[k] * df_3[k];

        t_10[k] = f_2 * pf_1[k]
                  + pa_y[k] * df_4[k];

        t_11[k] = f_0 * dd_7[k]
                  + pa_y[k] * df_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_y, pb_x, pb_y, pb_z, dd_8, df_6, fp_9, \
                         fp_10, fp_11, fd_20, fd_21, fd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_y[k] * df_6[k];

        t_13[k] = f_1 * fp_9[k]
                  + pb_x[k] * fd_20[k];

        t_14[k] = f_1 * fp_10[k]
                  + pb_y[k] * fd_21[k];

        t_15[k] = f_0 * dd_8[k]
                  + f_1 * fp_11[k]
                  + pb_z[k] * fd_22[k];
    }
}

auto
compute_prim_ff_overlap_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pf, const size_t dd, const size_t df,
                          const size_t fp, const size_t fd, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_3 = buffer.data(pf + 3);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_11 = buffer.data(df + 11);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
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
    const auto *fd_6 = buffer.data(fd + 6);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dd_0, fp_0, fp_1, fp_2, \
                         fd_0, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pb_y[k] * fd_0[k];

        t_2[k] = pb_z[k] * fd_0[k];

        t_3[k] = f_1 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = f_1 * fp_2[k]
                 + pb_z[k] * fd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pa_z, pb_z, pf_1, pf_3, dd_0, \
                         df_0, df_2, df_3, fd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_y[k] * df_0[k];

        t_6[k] = f_2 * pf_1[k]
                 + pa_x[k] * df_2[k];

        t_7[k] = pa_z[k] * df_0[k];

        t_8[k] = f_2 * dd_0[k]
                 + pb_z[k] * fd_6[k];

        t_9[k] = f_2 * pf_3[k]
                 + pa_x[k] * df_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, dd_5, dd_6, dd_10, df_4, df_5, \
                         df_8, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * dd_5[k]
                  + pa_x[k] * df_4[k];

        t_11[k] = f_2 * dd_6[k]
                  + pb_x[k] * fd_9[k];

        t_12[k] = pa_x[k] * df_5[k];

        t_13[k] = f_0 * dd_10[k]
                  + pa_x[k] * df_8[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_x, pb_x, pb_z, dd_3, dd_12, df_11, \
                         fp_3, fd_10, fd_11, fd_12, fd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * dd_3[k]
                  + pb_z[k] * fd_10[k];

        t_15[k] = f_2 * dd_12[k]
                  + pb_x[k] * fd_11[k];

        t_16[k] = pa_x[k] * df_11[k];

        t_17[k] = f_1 * fp_3[k]
                  + pb_x[k] * fd_12[k];

        t_18[k] = pb_x[k] * fd_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_z, pb_y, pb_z, dd_6, dd_7, df_5, \
                         fp_4, fp_5, fd_13, fd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * dd_6[k]
                  + f_1 * fp_4[k]
                  + pb_y[k] * fd_13[k];

        t_20[k] = pb_z[k] * fd_13[k];

        t_21[k] = f_0 * dd_7[k]
                  + pb_y[k] * fd_14[k];

        t_22[k] = f_1 * fp_5[k]
                  + pb_z[k] * fd_14[k];

        t_23[k] = pa_z[k] * df_5[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pb_y, pb_z, pf_3, dd_6, dd_9, dd_11, \
                         df_7, df_9, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_2 * dd_6[k]
                  + pb_z[k] * fd_15[k];

        t_25[k] = f_1 * dd_9[k]
                  + pb_y[k] * fd_16[k];

        t_26[k] = f_2 * pf_3[k]
                  + pa_y[k] * df_7[k];

        t_27[k] = f_0 * dd_11[k]
                  + pa_y[k] * df_9[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, dd_8, dd_12, \
                         df_11, fp_6, fd_17, fd_18, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * dd_8[k]
                  + pb_z[k] * fd_17[k];

        t_29[k] = f_2 * dd_12[k]
                  + pb_y[k] * fd_18[k];

        t_30[k] = pa_y[k] * df_11[k];

        t_31[k] = f_1 * fp_6[k]
                  + pb_x[k] * fd_19[k];

        t_32[k] = pb_y[k] * fd_19[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, dd_12, fp_7, \
                         fp_8, fd_20, fd_21, fd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_x[k] * fd_20[k];

        t_34[k] = pb_x[k] * fd_22[k];

        t_35[k] = f_1 * fp_7[k]
                  + pb_y[k] * fd_20[k];

        t_36[k] = f_2 * fp_8[k]
                  + pb_y[k] * fd_21[k];

        t_37[k] = pb_y[k] * fd_22[k];

        t_38[k] = f_0 * dd_12[k]
                  + f_1 * fp_8[k]
                  + pb_z[k] * fd_22[k];
    }
}

auto
compute_prim_ff_overlap_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                           const size_t pb, const size_t pf, const size_t dd, const size_t df,
                           const size_t fp, const size_t fd, const size_t ncols,
                           const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_3 = buffer.data(pf + 3);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_20 = buffer.data(df + 20);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
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
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dd_0, fp_0, fp_1, fp_2, \
                         fd_0, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pb_y[k] * fd_0[k];

        t_2[k] = pb_z[k] * fd_0[k];

        t_3[k] = f_1 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = f_1 * fp_2[k]
                 + pb_z[k] * fd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pa_z, pf_1, pf_3, df_0, df_3, \
                         df_4, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_y[k] * df_0[k];

        t_6[k] = f_2 * pf_1[k]
                 + pa_x[k] * df_4[k];

        t_7[k] = pa_y[k] * df_3[k];

        t_8[k] = pa_z[k] * df_0[k];

        t_9[k] = f_2 * pf_3[k]
                 + pa_x[k] * df_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pb_x, dd_4, dd_5, dd_8, dd_10, \
                         df_7, df_9, df_14, fd_7, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * dd_4[k]
                  + pa_x[k] * df_7[k];

        t_11[k] = f_2 * dd_5[k]
                  + pb_x[k] * fd_7[k];

        t_12[k] = pa_x[k] * df_9[k];

        t_13[k] = f_0 * dd_8[k]
                  + pa_x[k] * df_14[k];

        t_14[k] = f_2 * dd_10[k]
                  + pb_x[k] * fd_8[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_x, pb_x, pb_y, dd_5, df_20, fp_3, \
                         fp_4, fd_9, fd_10, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_x[k] * df_20[k];

        t_16[k] = f_1 * fp_3[k]
                  + pb_x[k] * fd_9[k];

        t_17[k] = pb_x[k] * fd_10[k];

        t_18[k] = pb_x[k] * fd_11[k];

        t_19[k] = f_0 * dd_5[k]
                  + f_1 * fp_4[k]
                  + pb_y[k] * fd_10[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pa_z, pb_x, pb_z, pf_3, df_9, \
                         df_13, fp_5, fd_10, fd_11, fd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_z[k] * fd_10[k];

        t_21[k] = f_1 * fp_5[k]
                  + pb_z[k] * fd_11[k];

        t_22[k] = pb_x[k] * fd_13[k];

        t_23[k] = pa_z[k] * df_9[k];

        t_24[k] = f_2 * pf_3[k]
                  + pa_y[k] * df_13[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pb_x, pb_y, dd_9, dd_10, df_18, \
                         df_20, fp_6, fd_14, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pb_x[k] * fd_14[k];

        t_26[k] = f_0 * dd_9[k]
                  + pa_y[k] * df_18[k];

        t_27[k] = f_2 * dd_10[k]
                  + pb_y[k] * fd_15[k];

        t_28[k] = pa_y[k] * df_20[k];

        t_29[k] = f_1 * fp_6[k]
                  + pb_x[k] * fd_16[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pb_x, pb_y, fp_7, fp_8, fd_16, \
                         fd_17, fd_18, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * fd_16[k];

        t_31[k] = pb_x[k] * fd_17[k];

        t_32[k] = pb_x[k] * fd_19[k];

        t_33[k] = f_1 * fp_7[k]
                  + pb_y[k] * fd_17[k];

        t_34[k] = f_2 * fp_8[k]
                  + pb_y[k] * fd_18[k];

        t_35[k] = pb_y[k] * fd_19[k];
    }

#pragma omp simd aligned(t_36, pb_z, dd_10, fp_8, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * dd_10[k]
                  + f_1 * fp_8[k]
                  + pb_z[k] * fd_19[k];
    }
}

auto
compute_prim_ff_overlap_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                           const size_t pb, const size_t pf, const size_t dd, const size_t df,
                           const size_t fp, const size_t fd, const size_t ncols,
                           const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_3 = buffer.data(pf + 3);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_18 = buffer.data(df + 18);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
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
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dd_0, fp_0, fp_1, fp_2, \
                         fd_0, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pb_y[k] * fd_0[k];

        t_2[k] = pb_z[k] * fd_0[k];

        t_3[k] = f_1 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = f_1 * fp_2[k]
                 + pb_z[k] * fd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pa_z, pb_x, pf_1, pf_3, dd_5, \
                         df_0, df_4, df_5, fd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_y[k] * df_0[k];

        t_6[k] = f_2 * pf_1[k]
                 + pa_x[k] * df_4[k];

        t_7[k] = pa_z[k] * df_0[k];

        t_8[k] = f_2 * pf_3[k]
                 + pa_x[k] * df_5[k];

        t_9[k] = f_2 * dd_5[k]
                 + pb_x[k] * fd_7[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pa_x, pb_x, dd_10, df_8, df_18, \
                         fp_3, fd_8, fd_9, fd_10, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * df_8[k];

        t_11[k] = f_2 * dd_10[k]
                  + pb_x[k] * fd_8[k];

        t_12[k] = pa_x[k] * df_18[k];

        t_13[k] = f_1 * fp_3[k]
                  + pb_x[k] * fd_9[k];

        t_14[k] = pb_x[k] * fd_10[k];

        t_15[k] = pb_x[k] * fd_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_z, pb_x, pb_y, pb_z, dd_5, df_8, \
                         fp_4, fp_5, fd_10, fd_11, fd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * dd_5[k]
                  + f_1 * fp_4[k]
                  + pb_y[k] * fd_10[k];

        t_17[k] = pb_z[k] * fd_10[k];

        t_18[k] = f_1 * fp_5[k]
                  + pb_z[k] * fd_11[k];

        t_19[k] = pb_x[k] * fd_13[k];

        t_20[k] = pa_z[k] * df_8[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pb_x, pb_y, pf_3, dd_9, dd_10, \
                         df_11, df_16, df_18, fd_14, fd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_2 * pf_3[k]
                  + pa_y[k] * df_11[k];

        t_22[k] = pb_x[k] * fd_14[k];

        t_23[k] = f_0 * dd_9[k]
                  + pa_y[k] * df_16[k];

        t_24[k] = f_2 * dd_10[k]
                  + pb_y[k] * fd_15[k];

        t_25[k] = pa_y[k] * df_18[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, t_32, pb_x, pb_y, fp_6, fp_7, \
                         fp_8, fd_16, fd_17, fd_18, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * fp_6[k]
                  + pb_x[k] * fd_16[k];

        t_27[k] = pb_y[k] * fd_16[k];

        t_28[k] = pb_x[k] * fd_17[k];

        t_29[k] = pb_x[k] * fd_19[k];

        t_30[k] = f_1 * fp_7[k]
                  + pb_y[k] * fd_17[k];

        t_31[k] = f_2 * fp_8[k]
                  + pb_y[k] * fd_18[k];

        t_32[k] = pb_y[k] * fd_19[k];
    }

#pragma omp simd aligned(t_33, pb_z, dd_10, fp_8, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * dd_10[k]
                  + f_1 * fp_8[k]
                  + pb_z[k] * fd_19[k];
    }
}

auto
compute_prim_ff_overlap_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                           const size_t pb, const size_t pf, const size_t dd, const size_t df,
                           const size_t fp, const size_t fd, const size_t ncols,
                           const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, pb_x, pf_0, dd_0, df_0, df_1, \
                         fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pa_y[k] * df_0[k];

        t_2[k] = f_2 * pf_0[k]
                 + pa_x[k] * df_1[k];

        t_3[k] = pa_z[k] * df_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_x, pa_z, pb_y, pf_1, dd_3, df_2, df_3, \
                         df_6, fp_3, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * pf_1[k]
                 + pa_x[k] * df_2[k];

        t_5[k] = pa_x[k] * df_3[k];

        t_6[k] = pa_x[k] * df_6[k];

        t_7[k] = f_0 * dd_3[k]
                 + f_1 * fp_3[k]
                 + pb_y[k] * fd_8[k];

        t_8[k] = pa_z[k] * df_3[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_y, pb_y, pf_1, dd_5, df_4, df_5, df_6, \
                         fp_7, fd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * pf_1[k]
                 + pa_y[k] * df_4[k];

        t_10[k] = f_0 * dd_5[k]
                  + pa_y[k] * df_5[k];

        t_11[k] = pa_y[k] * df_6[k];

        t_12[k] = f_1 * fp_7[k]
                  + pb_y[k] * fd_14[k];
    }

#pragma omp simd aligned(t_13, pb_z, dd_6, fp_8, fd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dd_6[k]
                  + f_1 * fp_8[k]
                  + pb_z[k] * fd_15[k];
    }
}

auto
compute_prim_ff_overlap_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                           const size_t pb, const size_t pf, const size_t dd, const size_t df,
                           const size_t fp, const size_t fd, const size_t ncols,
                           const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_3 = buffer.data(pf + 3);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_9 = buffer.data(df + 9);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, pb_x, pb_z, pf_1, dd_0, \
                         df_0, df_2, fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pb_z[k] * fd_0[k];

        t_2[k] = pa_y[k] * df_0[k];

        t_3[k] = f_2 * pf_1[k]
                 + pa_x[k] * df_2[k];

        t_4[k] = pa_z[k] * df_0[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pb_x, pf_3, dd_4, dd_8, df_3, df_4, \
                         df_9, fd_7, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_2 * pf_3[k]
                 + pa_x[k] * df_3[k];

        t_6[k] = f_2 * dd_4[k]
                 + pb_x[k] * fd_7[k];

        t_7[k] = pa_x[k] * df_4[k];

        t_8[k] = f_2 * dd_8[k]
                 + pb_x[k] * fd_9[k];

        t_9[k] = pa_x[k] * df_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_z, pb_x, pb_y, pb_z, dd_4, df_4, fp_2, \
                         fp_3, fd_10, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * fp_2[k]
                  + pb_x[k] * fd_10[k];

        t_11[k] = f_0 * dd_4[k]
                  + f_1 * fp_3[k]
                  + pb_y[k] * fd_11[k];

        t_12[k] = pb_z[k] * fd_11[k];

        t_13[k] = pa_z[k] * df_4[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pb_y, pf_3, dd_7, dd_8, df_6, df_7, \
                         df_9, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_2 * pf_3[k]
                  + pa_y[k] * df_6[k];

        t_15[k] = f_0 * dd_7[k]
                  + pa_y[k] * df_7[k];

        t_16[k] = f_2 * dd_8[k]
                  + pb_y[k] * fd_16[k];

        t_17[k] = pa_y[k] * df_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_x, pb_y, pb_z, dd_8, fp_5, fp_6, fp_7, \
                         fd_17, fd_18, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_1 * fp_5[k]
                  + pb_x[k] * fd_17[k];

        t_19[k] = f_1 * fp_6[k]
                  + pb_y[k] * fd_18[k];

        t_20[k] = pb_y[k] * fd_19[k];

        t_21[k] = f_0 * dd_8[k]
                  + f_1 * fp_7[k]
                  + pb_z[k] * fd_19[k];
    }
}

auto
compute_prim_ff_overlap_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                           const size_t pb, const size_t pf, const size_t dd, const size_t df,
                           const size_t fp, const size_t fd, const size_t ncols,
                           const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_3 = buffer.data(pf + 3);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_13 = buffer.data(df + 13);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_y, pb_x, pb_y, pb_z, dd_0, df_0, fp_0, \
                         fp_1, fd_0, fd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pb_y[k] * fd_0[k];

        t_2[k] = pb_z[k] * fd_0[k];

        t_3[k] = f_1 * fp_1[k]
                 + pb_z[k] * fd_1[k];

        t_4[k] = pa_y[k] * df_0[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_x, pa_z, pf_1, pf_3, dd_3, df_0, df_2, df_3, \
                         df_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_2 * pf_1[k]
                 + pa_x[k] * df_2[k];

        t_6[k] = pa_z[k] * df_0[k];

        t_7[k] = f_2 * pf_3[k]
                 + pa_x[k] * df_3[k];

        t_8[k] = f_0 * dd_3[k]
                 + pa_x[k] * df_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_x, pb_x, dd_4, dd_6, dd_8, df_6, \
                         df_9, df_13, fd_6, fd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * dd_4[k]
                 + pb_x[k] * fd_6[k];

        t_10[k] = pa_x[k] * df_6[k];

        t_11[k] = f_0 * dd_6[k]
                  + pa_x[k] * df_9[k];

        t_12[k] = f_2 * dd_8[k]
                  + pb_x[k] * fd_7[k];

        t_13[k] = pa_x[k] * df_13[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, dd_4, fp_2, fp_3, \
                         fp_4, fd_8, fd_9, fd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * fp_2[k]
                  + pb_x[k] * fd_8[k];

        t_15[k] = pb_x[k] * fd_9[k];

        t_16[k] = f_0 * dd_4[k]
                  + f_1 * fp_3[k]
                  + pb_y[k] * fd_9[k];

        t_17[k] = pb_z[k] * fd_9[k];

        t_18[k] = f_1 * fp_4[k]
                  + pb_z[k] * fd_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pa_z, pb_y, pf_3, dd_7, dd_8, \
                         df_6, df_8, df_11, df_13, fd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_z[k] * df_6[k];

        t_20[k] = f_2 * pf_3[k]
                  + pa_y[k] * df_8[k];

        t_21[k] = f_0 * dd_7[k]
                  + pa_y[k] * df_11[k];

        t_22[k] = f_2 * dd_8[k]
                  + pb_y[k] * fd_13[k];

        t_23[k] = pa_y[k] * df_13[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, pb_x, pb_y, fp_5, fp_6, fd_14, \
                         fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_1 * fp_5[k]
                  + pb_x[k] * fd_14[k];

        t_25[k] = pb_y[k] * fd_14[k];

        t_26[k] = pb_x[k] * fd_15[k];

        t_27[k] = pb_x[k] * fd_16[k];

        t_28[k] = f_1 * fp_6[k]
                  + pb_y[k] * fd_15[k];

        t_29[k] = pb_y[k] * fd_16[k];
    }

#pragma omp simd aligned(t_30, pb_z, dd_8, fp_7, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * dd_8[k]
                  + f_1 * fp_7[k]
                  + pb_z[k] * fd_16[k];
    }
}

auto
compute_prim_ff_overlap_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                           const size_t pb, const size_t pf, const size_t dd, const size_t df,
                           const size_t fp, const size_t fd, const size_t ncols,
                           const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;

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

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_3 = buffer.data(pf + 3);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_13 = buffer.data(df + 13);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_y, pb_x, pb_y, pb_z, dd_0, df_0, fp_0, \
                         fp_1, fd_0, fd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pb_y[k] * fd_0[k];

        t_2[k] = pb_z[k] * fd_0[k];

        t_3[k] = f_1 * fp_1[k]
                 + pb_z[k] * fd_1[k];

        t_4[k] = pa_y[k] * df_0[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_z, pb_x, pf_1, pf_3, dd_4, df_0, \
                         df_2, df_3, df_6, fd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_2 * pf_1[k]
                 + pa_x[k] * df_2[k];

        t_6[k] = pa_z[k] * df_0[k];

        t_7[k] = f_2 * pf_3[k]
                 + pa_x[k] * df_3[k];

        t_8[k] = f_2 * dd_4[k]
                 + pb_x[k] * fd_6[k];

        t_9[k] = pa_x[k] * df_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pb_x, pb_y, dd_4, dd_8, df_13, \
                         fp_2, fp_3, fd_7, fd_8, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * dd_8[k]
                  + pb_x[k] * fd_7[k];

        t_11[k] = pa_x[k] * df_13[k];

        t_12[k] = f_1 * fp_2[k]
                  + pb_x[k] * fd_8[k];

        t_13[k] = pb_x[k] * fd_9[k];

        t_14[k] = f_0 * dd_4[k]
                  + f_1 * fp_3[k]
                  + pb_y[k] * fd_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_y, pa_z, pb_z, pf_3, df_6, df_8, fp_4, \
                         fd_9, fd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_z[k] * fd_9[k];

        t_16[k] = f_1 * fp_4[k]
                  + pb_z[k] * fd_10[k];

        t_17[k] = pa_z[k] * df_6[k];

        t_18[k] = f_2 * pf_3[k]
                  + pa_y[k] * df_8[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, pa_y, pb_x, pb_y, dd_8, df_13, \
                         fp_5, fd_13, fd_14, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * dd_8[k]
                  + pb_y[k] * fd_13[k];

        t_20[k] = pa_y[k] * df_13[k];

        t_21[k] = f_1 * fp_5[k]
                  + pb_x[k] * fd_14[k];

        t_22[k] = pb_y[k] * fd_14[k];

        t_23[k] = pb_x[k] * fd_15[k];

        t_24[k] = pb_x[k] * fd_16[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pb_y, pb_z, dd_8, fp_6, fp_7, fd_15, \
                         fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * fp_6[k]
                  + pb_y[k] * fd_15[k];

        t_26[k] = pb_y[k] * fd_16[k];

        t_27[k] = f_0 * dd_8[k]
                  + f_1 * fp_7[k]
                  + pb_z[k] * fd_16[k];
    }
}

}  // namespace simdovl
