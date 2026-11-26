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


#include "SimdElectronRepulsionVrrRecFF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_ff_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dd, const size_t df,
                                     const size_t fp0, const size_t fp1, const size_t fd,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_17 = buffer.data(dd + 17);
    const auto *dd_18 = buffer.data(dd + 18);
    const auto *dd_21 = buffer.data(dd + 21);
    const auto *dd_23 = buffer.data(dd + 23);
    const auto *dd_27 = buffer.data(dd + 27);
    const auto *dd_28 = buffer.data(dd + 28);
    const auto *dd_29 = buffer.data(dd + 29);
    const auto *dd_30 = buffer.data(dd + 30);
    const auto *dd_33 = buffer.data(dd + 33);
    const auto *dd_35 = buffer.data(dd + 35);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_30 = buffer.data(df + 30);
    const auto *df_31 = buffer.data(df + 31);
    const auto *df_36 = buffer.data(df + 36);
    const auto *df_38 = buffer.data(df + 38);
    const auto *df_39 = buffer.data(df + 39);
    const auto *df_46 = buffer.data(df + 46);
    const auto *df_47 = buffer.data(df + 47);
    const auto *df_48 = buffer.data(df + 48);
    const auto *df_49 = buffer.data(df + 49);
    const auto *df_50 = buffer.data(df + 50);
    const auto *df_52 = buffer.data(df + 52);
    const auto *df_56 = buffer.data(df + 56);
    const auto *df_57 = buffer.data(df + 57);
    const auto *df_59 = buffer.data(df + 59);

    const auto *fp0_0 = buffer.data(fp0 + 0);
    const auto *fp0_1 = buffer.data(fp0 + 1);
    const auto *fp0_2 = buffer.data(fp0 + 2);
    const auto *fp0_18 = buffer.data(fp0 + 18);
    const auto *fp0_19 = buffer.data(fp0 + 19);
    const auto *fp0_20 = buffer.data(fp0 + 20);
    const auto *fp0_27 = buffer.data(fp0 + 27);
    const auto *fp0_28 = buffer.data(fp0 + 28);
    const auto *fp0_29 = buffer.data(fp0 + 29);

    const auto *fp1_0 = buffer.data(fp1 + 0);
    const auto *fp1_1 = buffer.data(fp1 + 1);
    const auto *fp1_2 = buffer.data(fp1 + 2);
    const auto *fp1_18 = buffer.data(fp1 + 18);
    const auto *fp1_19 = buffer.data(fp1 + 19);
    const auto *fp1_20 = buffer.data(fp1 + 20);
    const auto *fp1_27 = buffer.data(fp1 + 27);
    const auto *fp1_28 = buffer.data(fp1 + 28);
    const auto *fp1_29 = buffer.data(fp1 + 29);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_23 = buffer.data(fd + 23);
    const auto *fd_28 = buffer.data(fd + 28);
    const auto *fd_30 = buffer.data(fd + 30);
    const auto *fd_32 = buffer.data(fd + 32);
    const auto *fd_33 = buffer.data(fd + 33);
    const auto *fd_35 = buffer.data(fd + 35);
    const auto *fd_36 = buffer.data(fd + 36);
    const auto *fd_39 = buffer.data(fd + 39);
    const auto *fd_40 = buffer.data(fd + 40);
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
    const auto *fd_58 = buffer.data(fd + 58);
    const auto *fd_59 = buffer.data(fd + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, dd_0, dd_3, fp0_0, fp1_0, \
                         fd_0, fd_2, fd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 + f_1 * fp0_0[k]
                 - f_2 * fp1_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pb_y[k] * fd_0[k];

        t_2[k] = pb_z[k] * fd_0[k];

        t_3[k] = f_0 * dd_3[k]
                 + pb_x[k] * fd_3[k];

        t_4[k] = pb_y[k] * fd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, dd_5, fp0_1, fp0_2, fp1_1, \
                         fp1_2, fd_3, fd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * dd_5[k]
                 + pb_x[k] * fd_5[k];

        t_6[k] = f_1 * fp0_1[k]
                 - f_2 * fp1_1[k]
                 + pb_y[k] * fd_3[k];

        t_7[k] = pb_z[k] * fd_3[k];

        t_8[k] = pb_y[k] * fd_5[k];

        t_9[k] = f_1 * fp0_2[k]
                 - f_2 * fp1_2[k]
                 + pb_z[k] * fd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, pb_y, pb_z, dd_0, dd_9, \
                         df_0, fd_6, fd_7, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * df_0[k];

        t_11[k] = f_3 * dd_0[k]
                  + pb_y[k] * fd_6[k];

        t_12[k] = pb_z[k] * fd_6[k];

        t_13[k] = f_4 * dd_9[k]
                  + pb_x[k] * fd_9[k];

        t_14[k] = pb_z[k] * fd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, dd_3, dd_5, df_5, \
                         df_6, df_9, fd_9, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * df_5[k];

        t_16[k] = f_0 * dd_3[k]
                  + pa_y[k] * df_6[k];

        t_17[k] = pb_z[k] * fd_9[k];

        t_18[k] = f_3 * dd_5[k]
                  + pb_y[k] * fd_11[k];

        t_19[k] = pa_y[k] * df_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, dd_0, df_0, df_3, \
                         fd_12, fd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * df_0[k];

        t_21[k] = pb_y[k] * fd_12[k];

        t_22[k] = f_3 * dd_0[k]
                  + pb_z[k] * fd_12[k];

        t_23[k] = pa_z[k] * df_3[k];

        t_24[k] = pb_y[k] * fd_14[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, pb_z, dd_3, dd_5, \
                         dd_17, df_6, df_9, fd_15, fd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * dd_17[k]
                  + pb_x[k] * fd_17[k];

        t_26[k] = pa_z[k] * df_6[k];

        t_27[k] = f_3 * dd_3[k]
                  + pb_z[k] * fd_15[k];

        t_28[k] = pb_y[k] * fd_17[k];

        t_29[k] = f_0 * dd_5[k]
                  + pa_z[k] * df_9[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pa_x, pb_x, pb_y, pb_z, dd_6, dd_18, \
                         dd_21, df_30, fd_18, fd_19, fd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * dd_18[k]
                  + pa_x[k] * df_30[k];

        t_31[k] = f_4 * dd_6[k]
                  + pb_y[k] * fd_18[k];

        t_32[k] = pb_z[k] * fd_18[k];

        t_33[k] = f_3 * dd_21[k]
                  + pb_x[k] * fd_21[k];

        t_34[k] = pb_z[k] * fd_19[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_x, pb_x, pb_z, dd_23, df_36, df_38, \
                         df_39, fd_21, fd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_3 * dd_23[k]
                  + pb_x[k] * fd_23[k];

        t_36[k] = pa_x[k] * df_36[k];

        t_37[k] = pb_z[k] * fd_21[k];

        t_38[k] = pa_x[k] * df_38[k];

        t_39[k] = pa_x[k] * df_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, pa_y, pa_z, pb_x, dd_28, df_11, \
                         df_13, df_20, df_22, df_25, fd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * df_20[k];

        t_41[k] = pa_z[k] * df_11[k];

        t_42[k] = pa_y[k] * df_22[k];

        t_43[k] = pa_z[k] * df_13[k];

        t_44[k] = f_3 * dd_28[k]
                  + pb_x[k] * fd_28[k];

        t_45[k] = pa_y[k] * df_25[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, t_51, pa_x, pb_y, dd_30, df_46, df_47, \
                         df_48, df_49, df_50, fd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_x[k] * df_46[k];

        t_47[k] = pa_x[k] * df_47[k];

        t_48[k] = pa_x[k] * df_48[k];

        t_49[k] = pa_x[k] * df_49[k];

        t_50[k] = f_0 * dd_30[k]
                  + pa_x[k] * df_50[k];

        t_51[k] = pb_y[k] * fd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_x, pb_y, pb_z, dd_12, dd_33, dd_35, fd_30, \
                         fd_32, fd_33, fd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_4 * dd_12[k]
                  + pb_z[k] * fd_30[k];

        t_53[k] = f_3 * dd_33[k]
                  + pb_x[k] * fd_33[k];

        t_54[k] = pb_y[k] * fd_32[k];

        t_55[k] = f_3 * dd_35[k]
                  + pb_x[k] * fd_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_x, pb_x, pb_y, df_56, df_57, df_59, \
                         fp0_18, fp1_18, fd_35, fd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_x[k] * df_56[k];

        t_57[k] = pa_x[k] * df_57[k];

        t_58[k] = pb_y[k] * fd_35[k];

        t_59[k] = pa_x[k] * df_59[k];

        t_60[k] = f_1 * fp0_18[k]
                  - f_2 * fp1_18[k]
                  + pb_x[k] * fd_36[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pb_x, pb_y, pb_z, dd_18, fd_36, fd_39, \
                         fd_40, fd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_0 * dd_18[k]
                  + pb_y[k] * fd_36[k];

        t_62[k] = pb_z[k] * fd_36[k];

        t_63[k] = pb_x[k] * fd_39[k];

        t_64[k] = pb_x[k] * fd_40[k];

        t_65[k] = pb_x[k] * fd_41[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pb_y, pb_z, dd_21, dd_23, fp0_19, fp0_20, \
                         fp1_19, fp1_20, fd_39, fd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_0 * dd_21[k]
                  + f_1 * fp0_19[k]
                  - f_2 * fp1_19[k]
                  + pb_y[k] * fd_39[k];

        t_67[k] = pb_z[k] * fd_39[k];

        t_68[k] = f_0 * dd_23[k]
                  + pb_y[k] * fd_41[k];

        t_69[k] = f_1 * fp0_20[k]
                  - f_2 * fp1_20[k]
                  + pb_z[k] * fd_41[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, t_75, pa_z, pb_x, pb_z, dd_18, df_30, \
                         df_31, fd_42, fd_45, fd_46, fd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * df_30[k];

        t_71[k] = pa_z[k] * df_31[k];

        t_72[k] = f_3 * dd_18[k]
                  + pb_z[k] * fd_42[k];

        t_73[k] = pb_x[k] * fd_45[k];

        t_74[k] = pb_x[k] * fd_46[k];

        t_75[k] = pb_x[k] * fd_47[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_z, pb_y, pb_z, dd_21, dd_23, dd_29, df_36, \
                         df_39, fd_45, fd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pa_z[k] * df_36[k];

        t_77[k] = f_3 * dd_21[k]
                  + pb_z[k] * fd_45[k];

        t_78[k] = f_4 * dd_29[k]
                  + pb_y[k] * fd_47[k];

        t_79[k] = f_0 * dd_23[k]
                  + pa_z[k] * df_39[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, t_85, pa_y, pb_x, pb_y, dd_30, df_50, \
                         df_52, fd_48, fd_51, fd_52, fd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pa_y[k] * df_50[k];

        t_81[k] = f_3 * dd_30[k]
                  + pb_y[k] * fd_48[k];

        t_82[k] = pa_y[k] * df_52[k];

        t_83[k] = pb_x[k] * fd_51[k];

        t_84[k] = pb_x[k] * fd_52[k];

        t_85[k] = pb_x[k] * fd_53[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_y, pb_y, pb_z, dd_27, dd_33, dd_35, df_56, \
                         df_59, fd_51, fd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_0 * dd_33[k]
                  + pa_y[k] * df_56[k];

        t_87[k] = f_4 * dd_27[k]
                  + pb_z[k] * fd_51[k];

        t_88[k] = f_3 * dd_35[k]
                  + pb_y[k] * fd_53[k];

        t_89[k] = pa_y[k] * df_59[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, dd_30, fp0_27, \
                         fp1_27, fd_54, fd_57, fd_58, fd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_1 * fp0_27[k]
                  - f_2 * fp1_27[k]
                  + pb_x[k] * fd_54[k];

        t_91[k] = pb_y[k] * fd_54[k];

        t_92[k] = f_0 * dd_30[k]
                  + pb_z[k] * fd_54[k];

        t_93[k] = pb_x[k] * fd_57[k];

        t_94[k] = pb_x[k] * fd_58[k];

        t_95[k] = pb_x[k] * fd_59[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pb_y, pb_z, dd_33, dd_35, fp0_28, fp0_29, \
                         fp1_28, fp1_29, fd_57, fd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * fp0_28[k]
                  - f_2 * fp1_28[k]
                  + pb_y[k] * fd_57[k];

        t_97[k] = f_0 * dd_33[k]
                  + pb_z[k] * fd_57[k];

        t_98[k] = pb_y[k] * fd_59[k];

        t_99[k] = f_0 * dd_35[k]
                  + f_1 * fp0_29[k]
                  - f_2 * fp1_29[k]
                  + pb_z[k] * fd_59[k];
    }
}

}  // namespace simdt2ceri
