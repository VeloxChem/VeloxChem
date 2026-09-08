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


#include "SimdElectronRepulsionVrrRecDG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_dg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pf, const size_t pg,
                                     const size_t dd0, const size_t dd1, const size_t df,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;

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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_12 = buffer.data(pf + 12);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_14 = buffer.data(pf + 14);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_11 = buffer.data(pg + 11);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_14 = buffer.data(pg + 14);
    const auto *pg_15 = buffer.data(pg + 15);
    const auto *pg_16 = buffer.data(pg + 16);
    const auto *pg_17 = buffer.data(pg + 17);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_3 = buffer.data(dd1 + 3);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);

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
    const auto *df_36 = buffer.data(df + 36);
    const auto *df_37 = buffer.data(df + 37);
    const auto *df_38 = buffer.data(df + 38);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, pf_0, dd0_0, dd1_0, \
                         df_0, df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = pb_y[k] * df_0[k];

        t_2[k] = pb_z[k] * df_0[k];

        t_3[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_4[k] = pb_y[k] * df_2[k];

        t_5[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, pf_1, pf_2, dd0_1, dd1_1, \
                         df_3, df_4, df_5, df_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * pf_1[k]
                 + pb_x[k] * df_5[k];

        t_7[k] = pb_z[k] * df_3[k];

        t_8[k] = pb_y[k] * df_4[k];

        t_9[k] = f_0 * pf_2[k]
                 + pb_x[k] * df_7[k];

        t_10[k] = f_1 * dd0_1[k]
                  - f_2 * dd1_1[k]
                  + pb_y[k] * df_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, pg_0, dd0_2, dd1_2, \
                         df_5, df_6, df_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * df_5[k];

        t_12[k] = f_3 * dd0_2[k]
                  - f_4 * dd1_2[k]
                  + pb_y[k] * df_6[k];

        t_13[k] = pb_y[k] * df_7[k];

        t_14[k] = f_1 * dd0_2[k]
                  - f_2 * dd1_2[k]
                  + pb_z[k] * df_7[k];

        t_15[k] = pa_y[k] * pg_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_x, pa_y, pb_y, pb_z, pf_0, pf_4, \
                         pg_2, pg_6, df_8, df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * pf_0[k]
                  + pb_y[k] * df_8[k];

        t_17[k] = pb_z[k] * df_8[k];

        t_18[k] = f_0 * pf_4[k]
                  + pa_x[k] * pg_6[k];

        t_19[k] = pb_z[k] * df_9[k];

        t_20[k] = pa_y[k] * pg_2[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_x, pa_y, pb_x, pb_z, pf_5, pf_6, \
                         pg_4, pg_7, df_10, df_11, df_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * pf_5[k]
                  + pb_x[k] * df_11[k];

        t_22[k] = pb_z[k] * df_10[k];

        t_23[k] = f_5 * pf_6[k]
                  + pb_x[k] * df_12[k];

        t_24[k] = pa_y[k] * pg_4[k];

        t_25[k] = pa_x[k] * pg_7[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, pa_x, pa_z, pb_y, pb_z, pg_0, \
                         pg_8, pg_9, pg_10, df_11, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * df_11[k];

        t_27[k] = pa_x[k] * pg_8[k];

        t_28[k] = pa_x[k] * pg_9[k];

        t_29[k] = pa_x[k] * pg_10[k];

        t_30[k] = pa_z[k] * pg_0[k];

        t_31[k] = pb_y[k] * df_13[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, pa_x, pa_z, pb_y, pb_z, pf_0, pf_10, \
                         pg_1, pg_3, pg_13, df_13, df_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_5 * pf_0[k]
                  + pb_z[k] * df_13[k];

        t_33[k] = pa_z[k] * pg_1[k];

        t_34[k] = pb_y[k] * df_14[k];

        t_35[k] = f_0 * pf_10[k]
                  + pa_x[k] * pg_13[k];

        t_36[k] = pa_z[k] * pg_3[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pa_x, pb_x, pb_y, pf_12, pf_14, pg_14, \
                         pg_15, df_15, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * pf_12[k]
                  + pb_x[k] * df_16[k];

        t_38[k] = pb_y[k] * df_15[k];

        t_39[k] = f_5 * pf_14[k]
                  + pb_x[k] * df_17[k];

        t_40[k] = pa_x[k] * pg_14[k];

        t_41[k] = pa_x[k] * pg_15[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_x, pb_x, pb_y, pf_3, pg_16, pg_17, \
                         dd0_3, dd1_3, df_17, df_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_x[k] * pg_16[k];

        t_43[k] = pb_y[k] * df_17[k];

        t_44[k] = pa_x[k] * pg_17[k];

        t_45[k] = f_1 * dd0_3[k]
                  - f_2 * dd1_3[k]
                  + pb_x[k] * df_18[k];

        t_46[k] = f_0 * pf_3[k]
                  + pb_y[k] * df_18[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pb_x, pb_z, dd0_4, dd0_5, dd1_4, dd1_5, \
                         df_18, df_19, df_20, df_21, df_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_z[k] * df_18[k];

        t_48[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_x[k] * df_20[k];

        t_49[k] = pb_z[k] * df_19[k];

        t_50[k] = f_3 * dd0_5[k]
                  - f_4 * dd1_5[k]
                  + pb_x[k] * df_21[k];

        t_51[k] = pb_x[k] * df_22[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, t_57, pb_x, pb_y, pb_z, pf_5, dd0_4, \
                         dd1_4, df_22, df_23, df_24, df_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_x[k] * df_23[k];

        t_53[k] = pb_x[k] * df_24[k];

        t_54[k] = pb_x[k] * df_25[k];

        t_55[k] = f_0 * pf_5[k]
                  + f_1 * dd0_4[k]
                  - f_2 * dd1_4[k]
                  + pb_y[k] * df_22[k];

        t_56[k] = pb_z[k] * df_22[k];

        t_57[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_z[k] * df_23[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_y, pa_z, pb_y, pb_z, pf_7, pg_5, \
                         pg_11, pg_12, dd0_5, dd1_5, df_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_0 * pf_7[k]
                  + pb_y[k] * df_25[k];

        t_59[k] = f_1 * dd0_5[k]
                  - f_2 * dd1_5[k]
                  + pb_z[k] * df_25[k];

        t_60[k] = pa_y[k] * pg_11[k];

        t_61[k] = pa_z[k] * pg_5[k];

        t_62[k] = pa_y[k] * pg_12[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pa_z, pb_x, pb_y, pf_9, pg_6, \
                         pg_13, df_26, df_27, df_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_z[k] * pg_6[k];

        t_64[k] = f_5 * pf_9[k]
                  + pb_y[k] * df_26[k];

        t_65[k] = pa_y[k] * pg_13[k];

        t_66[k] = pb_x[k] * df_27[k];

        t_67[k] = pb_x[k] * df_28[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_y, pa_z, pb_x, pb_z, pf_5, pf_13, \
                         pg_7, pg_16, df_27, df_29, df_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * df_29[k];

        t_69[k] = pb_x[k] * df_30[k];

        t_70[k] = pa_z[k] * pg_7[k];

        t_71[k] = f_5 * pf_5[k]
                  + pb_z[k] * df_27[k];

        t_72[k] = f_0 * pf_13[k]
                  + pa_y[k] * pg_16[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_y, pb_x, pb_y, pb_z, pf_8, pf_14, \
                         pg_17, dd0_6, dd1_6, df_30, df_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_5 * pf_14[k]
                  + pb_y[k] * df_30[k];

        t_74[k] = pa_y[k] * pg_17[k];

        t_75[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_6[k]
                  + pb_x[k] * df_31[k];

        t_76[k] = pb_y[k] * df_31[k];

        t_77[k] = f_0 * pf_8[k]
                  + pb_z[k] * df_31[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pb_x, pb_y, dd0_7, dd0_8, dd1_7, dd1_8, \
                         df_32, df_33, df_34, df_35, df_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_3 * dd0_7[k]
                  - f_4 * dd1_7[k]
                  + pb_x[k] * df_33[k];

        t_79[k] = pb_y[k] * df_32[k];

        t_80[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_x[k] * df_34[k];

        t_81[k] = pb_x[k] * df_35[k];

        t_82[k] = pb_x[k] * df_36[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pb_x, pb_y, pb_z, pf_11, dd0_7, dd0_8, \
                         dd1_7, dd1_8, df_35, df_37, df_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pb_x[k] * df_37[k];

        t_84[k] = pb_x[k] * df_38[k];

        t_85[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_7[k]
                  + pb_y[k] * df_35[k];

        t_86[k] = f_0 * pf_11[k]
                  + pb_z[k] * df_35[k];

        t_87[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_y[k] * df_37[k];
    }

#pragma omp simd aligned(t_88, t_89, pb_y, pb_z, pf_14, dd0_8, dd1_8, \
                         df_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pb_y[k] * df_38[k];

        t_89[k] = f_0 * pf_14[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_8[k]
                  + pb_z[k] * df_38[k];
    }
}

auto
compute_prim_dg_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pf, const size_t pg,
                                     const size_t dd0, const size_t dd1, const size_t df,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;

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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_15 = buffer.data(pg + 15);
    const auto *pg_16 = buffer.data(pg + 16);
    const auto *pg_17 = buffer.data(pg + 17);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_3 = buffer.data(dd1 + 3);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pf_0, dd0_0, dd1_0, df_0, \
                         df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = pb_y[k] * df_0[k];

        t_2[k] = pb_z[k] * df_0[k];

        t_3[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_4[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pf_1, pf_2, dd0_1, dd0_2, dd1_1, \
                         dd1_2, df_3, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * pf_1[k]
                 + pb_x[k] * df_3[k];

        t_6[k] = f_0 * pf_2[k]
                 + pb_x[k] * df_5[k];

        t_7[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_3[k];

        t_8[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_4[k];

        t_9[k] = pb_y[k] * df_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, pb_y, pb_z, pf_0, pf_4, pg_0, \
                         pg_6, dd0_2, dd1_2, df_5, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * dd0_2[k]
                  - f_2 * dd1_2[k]
                  + pb_z[k] * df_5[k];

        t_11[k] = pa_y[k] * pg_0[k];

        t_12[k] = f_5 * pf_0[k]
                  + pb_y[k] * df_6[k];

        t_13[k] = f_0 * pf_4[k]
                  + pa_x[k] * pg_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_x, pa_y, pb_x, pf_5, pg_2, pg_7, \
                         pg_8, pg_9, df_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_y[k] * pg_2[k];

        t_15[k] = f_5 * pf_5[k]
                  + pb_x[k] * df_7[k];

        t_16[k] = pa_x[k] * pg_7[k];

        t_17[k] = pa_x[k] * pg_8[k];

        t_18[k] = pa_x[k] * pg_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pa_z, pb_z, pf_0, pf_8, pg_0, pg_1, \
                         pg_13, df_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_z[k] * pg_0[k];

        t_20[k] = f_5 * pf_0[k]
                  + pb_z[k] * df_8[k];

        t_21[k] = pa_z[k] * pg_1[k];

        t_22[k] = f_0 * pf_8[k]
                  + pa_x[k] * pg_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pa_x, pb_x, pf_11, pg_15, pg_16, pg_17, \
                         dd0_3, dd1_3, df_9, df_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * pf_11[k]
                  + pb_x[k] * df_9[k];

        t_24[k] = pa_x[k] * pg_15[k];

        t_25[k] = pa_x[k] * pg_16[k];

        t_26[k] = pa_x[k] * pg_17[k];

        t_27[k] = f_1 * dd0_3[k]
                  - f_2 * dd1_3[k]
                  + pb_x[k] * df_10[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_x, pb_y, pf_3, dd0_4, dd0_5, dd1_4, dd1_5, \
                         df_10, df_11, df_12, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * pf_3[k]
                  + pb_y[k] * df_10[k];

        t_29[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_x[k] * df_11[k];

        t_30[k] = f_3 * dd0_5[k]
                  - f_4 * dd1_5[k]
                  + pb_x[k] * df_12[k];

        t_31[k] = pb_x[k] * df_13[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, pb_x, pb_y, pb_z, pf_5, pf_6, dd0_4, \
                         dd1_4, df_13, df_14, df_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_x[k] * df_15[k];

        t_33[k] = f_0 * pf_5[k]
                  + f_1 * dd0_4[k]
                  - f_2 * dd1_4[k]
                  + pb_y[k] * df_13[k];

        t_34[k] = pb_z[k] * df_13[k];

        t_35[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_z[k] * df_14[k];

        t_36[k] = f_0 * pf_6[k]
                  + pb_y[k] * df_15[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pa_y, pa_z, pb_z, pg_6, pg_7, pg_12, \
                         pg_13, dd0_5, dd1_5, df_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_1 * dd0_5[k]
                  - f_2 * dd1_5[k]
                  + pb_z[k] * df_15[k];

        t_38[k] = pa_y[k] * pg_12[k];

        t_39[k] = pa_z[k] * pg_6[k];

        t_40[k] = pa_y[k] * pg_13[k];

        t_41[k] = pa_z[k] * pg_7[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pb_y, pb_z, pf_5, pf_10, pf_11, pg_16, \
                         pg_17, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * pf_5[k]
                  + pb_z[k] * df_16[k];

        t_43[k] = f_0 * pf_10[k]
                  + pa_y[k] * pg_16[k];

        t_44[k] = f_5 * pf_11[k]
                  + pb_y[k] * df_17[k];

        t_45[k] = pa_y[k] * pg_17[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_x, pb_z, pf_7, dd0_6, dd0_7, dd0_8, dd1_6, \
                         dd1_7, dd1_8, df_18, df_19, df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_6[k]
                  + pb_x[k] * df_18[k];

        t_47[k] = f_0 * pf_7[k]
                  + pb_z[k] * df_18[k];

        t_48[k] = f_3 * dd0_7[k]
                  - f_4 * dd1_7[k]
                  + pb_x[k] * df_19[k];

        t_49[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_x[k] * df_20[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pb_x, pb_y, pb_z, pf_9, dd0_7, dd0_8, \
                         dd1_7, dd1_8, df_21, df_22, df_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pb_x[k] * df_21[k];

        t_51[k] = pb_x[k] * df_23[k];

        t_52[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_7[k]
                  + pb_y[k] * df_21[k];

        t_53[k] = f_0 * pf_9[k]
                  + pb_z[k] * df_21[k];

        t_54[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_y[k] * df_22[k];
    }

#pragma omp simd aligned(t_55, t_56, pb_y, pb_z, pf_11, dd0_8, dd1_8, \
                         df_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pb_y[k] * df_23[k];

        t_56[k] = f_0 * pf_11[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_8[k]
                  + pb_z[k] * df_23[k];
    }
}

auto
compute_prim_dg_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pf, const size_t dd0, const size_t dd1,
                                     const size_t df, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_3 = buffer.data(dd1 + 3);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pf_0, dd0_0, dd1_0, df_0, \
                         df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = pb_y[k] * df_0[k];

        t_2[k] = pb_z[k] * df_0[k];

        t_3[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_4[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, dd0_1, dd0_2, dd1_1, dd1_2, df_3, \
                         df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_3[k];

        t_6[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_4[k];

        t_7[k] = pb_y[k] * df_5[k];

        t_8[k] = f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_x, dd0_3, dd0_4, dd0_5, dd1_3, dd1_4, \
                         dd1_5, df_6, df_7, df_8, df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_1 * dd0_3[k]
                 - f_2 * dd1_3[k]
                 + pb_x[k] * df_6[k];

        t_10[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_x[k] * df_7[k];

        t_11[k] = f_3 * dd0_5[k]
                  - f_4 * dd1_5[k]
                  + pb_x[k] * df_8[k];

        t_12[k] = pb_x[k] * df_9[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pb_x, pb_y, pb_z, pf_1, dd0_4, dd0_5, \
                         dd1_4, dd1_5, df_9, df_10, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_x[k] * df_11[k];

        t_14[k] = f_0 * pf_1[k]
                  + f_1 * dd0_4[k]
                  - f_2 * dd1_4[k]
                  + pb_y[k] * df_9[k];

        t_15[k] = pb_z[k] * df_9[k];

        t_16[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_z[k] * df_10[k];

        t_17[k] = f_1 * dd0_5[k]
                  - f_2 * dd1_5[k]
                  + pb_z[k] * df_11[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_x, dd0_6, dd0_7, dd0_8, dd1_6, dd1_7, \
                         dd1_8, df_12, df_13, df_14, df_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_6[k]
                  + pb_x[k] * df_12[k];

        t_19[k] = f_3 * dd0_7[k]
                  - f_4 * dd1_7[k]
                  + pb_x[k] * df_13[k];

        t_20[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_x[k] * df_14[k];

        t_21[k] = pb_x[k] * df_15[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pb_x, pb_y, pb_z, pf_2, dd0_7, dd0_8, \
                         dd1_7, dd1_8, df_15, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_x[k] * df_17[k];

        t_23[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_7[k]
                  + pb_y[k] * df_15[k];

        t_24[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_y[k] * df_16[k];

        t_25[k] = pb_y[k] * df_17[k];

        t_26[k] = f_0 * pf_2[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_8[k]
                  + pb_z[k] * df_17[k];
    }
}

auto
compute_prim_dg_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pf, const size_t pg,
                                     const size_t dd0, const size_t dd1, const size_t df,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_5 = buffer.data(pg + 5);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_3 = buffer.data(dd1 + 3);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pf_0, dd0_0, dd1_0, df_0, \
                         df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = pb_y[k] * df_0[k];

        t_2[k] = pb_z[k] * df_0[k];

        t_3[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_4[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, pg_0, dd0_1, dd0_2, dd1_1, \
                         dd1_2, df_3, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_3[k];

        t_6[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_4[k];

        t_7[k] = pb_y[k] * df_5[k];

        t_8[k] = f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_5[k];

        t_9[k] = pa_y[k] * pg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_z, pb_x, pg_0, pg_1, pg_5, dd0_3, \
                         dd1_3, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * pg_1[k];

        t_11[k] = pa_z[k] * pg_0[k];

        t_12[k] = pa_x[k] * pg_5[k];

        t_13[k] = f_1 * dd0_3[k]
                  - f_2 * dd1_3[k]
                  + pb_x[k] * df_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pb_x, pb_y, pf_1, dd0_4, dd0_5, dd1_4, \
                         dd1_5, df_7, df_8, df_9, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_x[k] * df_7[k];

        t_15[k] = f_3 * dd0_5[k]
                  - f_4 * dd1_5[k]
                  + pb_x[k] * df_8[k];

        t_16[k] = pb_x[k] * df_9[k];

        t_17[k] = pb_x[k] * df_11[k];

        t_18[k] = f_0 * pf_1[k]
                  + f_1 * dd0_4[k]
                  - f_2 * dd1_4[k]
                  + pb_y[k] * df_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_z, pb_z, pg_1, dd0_4, dd0_5, dd1_4, dd1_5, \
                         df_9, df_10, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_z[k] * df_9[k];

        t_20[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_z[k] * df_10[k];

        t_21[k] = f_1 * dd0_5[k]
                  - f_2 * dd1_5[k]
                  + pb_z[k] * df_11[k];

        t_22[k] = pa_z[k] * pg_1[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pb_x, pg_5, dd0_6, dd0_7, dd0_8, dd1_6, \
                         dd1_7, dd1_8, df_12, df_13, df_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_y[k] * pg_5[k];

        t_24[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_6[k]
                  + pb_x[k] * df_12[k];

        t_25[k] = f_3 * dd0_7[k]
                  - f_4 * dd1_7[k]
                  + pb_x[k] * df_13[k];

        t_26[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_x[k] * df_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pb_x, pb_y, dd0_7, dd0_8, dd1_7, dd1_8, \
                         df_15, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_x[k] * df_15[k];

        t_28[k] = pb_x[k] * df_17[k];

        t_29[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_7[k]
                  + pb_y[k] * df_15[k];

        t_30[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_y[k] * df_16[k];

        t_31[k] = pb_y[k] * df_17[k];
    }

#pragma omp simd aligned(t_32, pb_z, pf_2, dd0_8, dd1_8, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * pf_2[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_8[k]
                  + pb_z[k] * df_17[k];
    }
}

auto
compute_prim_dg_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pf, const size_t pg,
                                     const size_t dd0, const size_t dd1, const size_t df,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;

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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_14 = buffer.data(pf + 14);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);
    const auto *dd1_9 = buffer.data(dd1 + 9);
    const auto *dd1_12 = buffer.data(dd1 + 12);
    const auto *dd1_13 = buffer.data(dd1 + 13);
    const auto *dd1_14 = buffer.data(dd1 + 14);

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
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_26 = buffer.data(df + 26);
    const auto *df_28 = buffer.data(df + 28);
    const auto *df_29 = buffer.data(df + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pf_0, pf_1, dd0_0, dd1_0, df_0, \
                         df_1, df_2, df_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_2[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];

        t_3[k] = f_0 * pf_1[k]
                 + pb_x[k] * df_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, pf_2, dd0_1, dd0_2, dd1_1, \
                         dd1_2, df_3, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_0 * pf_2[k]
                 + pb_x[k] * df_5[k];

        t_5[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_3[k];

        t_6[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_4[k];

        t_7[k] = f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_5[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_z, pb_x, pb_y, pf_0, pf_4, pf_5, pg_0, \
                         pg_1, df_6, df_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * pf_0[k]
                 + pb_y[k] * df_6[k];

        t_9[k] = f_0 * pf_4[k]
                 + pa_x[k] * pg_1[k];

        t_10[k] = f_5 * pf_5[k]
                  + pb_x[k] * df_7[k];

        t_11[k] = pa_z[k] * pg_0[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, pb_z, pf_0, pf_10, pf_14, pg_3, \
                         dd0_3, dd1_7, df_8, df_10, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * pf_0[k]
                  + pb_z[k] * df_8[k];

        t_13[k] = f_0 * pf_10[k]
                  + pa_x[k] * pg_3[k];

        t_14[k] = f_5 * pf_14[k]
                  + pb_x[k] * df_10[k];

        t_15[k] = f_1 * dd0_3[k]
                  - f_2 * dd1_7[k]
                  + pb_x[k] * df_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_x, pb_y, pf_3, pf_5, dd0_4, dd0_5, dd1_8, \
                         dd1_9, df_11, df_12, df_13, df_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * pf_3[k]
                  + pb_y[k] * df_11[k];

        t_17[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_8[k]
                  + pb_x[k] * df_12[k];

        t_18[k] = f_3 * dd0_5[k]
                  - f_4 * dd1_9[k]
                  + pb_x[k] * df_13[k];

        t_19[k] = f_0 * pf_5[k]
                  + f_1 * dd0_4[k]
                  - f_2 * dd1_8[k]
                  + pb_y[k] * df_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_z, pb_y, pb_z, pf_7, pg_2, dd0_4, dd0_5, \
                         dd1_8, dd1_9, df_15, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_8[k]
                  + pb_z[k] * df_15[k];

        t_21[k] = f_0 * pf_7[k]
                  + pb_y[k] * df_17[k];

        t_22[k] = f_1 * dd0_5[k]
                  - f_2 * dd1_9[k]
                  + pb_z[k] * df_17[k];

        t_23[k] = pa_z[k] * pg_2[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pb_y, pb_z, pf_5, pf_13, pf_14, pg_4, \
                         pg_5, df_18, df_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * pf_5[k]
                  + pb_z[k] * df_18[k];

        t_25[k] = f_0 * pf_13[k]
                  + pa_y[k] * pg_4[k];

        t_26[k] = f_5 * pf_14[k]
                  + pb_y[k] * df_21[k];

        t_27[k] = pa_y[k] * pg_5[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_x, pb_z, pf_8, dd0_6, dd0_7, dd0_8, \
                         dd1_12, dd1_13, dd1_14, df_22, df_24, df_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_12[k]
                  + pb_x[k] * df_22[k];

        t_29[k] = f_0 * pf_8[k]
                  + pb_z[k] * df_22[k];

        t_30[k] = f_3 * dd0_7[k]
                  - f_4 * dd1_13[k]
                  + pb_x[k] * df_24[k];

        t_31[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_14[k]
                  + pb_x[k] * df_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_y, pb_z, pf_11, pf_14, dd0_7, dd0_8, \
                         dd1_13, dd1_14, df_26, df_28, df_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_13[k]
                  + pb_y[k] * df_26[k];

        t_33[k] = f_0 * pf_11[k]
                  + pb_z[k] * df_26[k];

        t_34[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_14[k]
                  + pb_y[k] * df_28[k];

        t_35[k] = f_0 * pf_14[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_14[k]
                  + pb_z[k] * df_29[k];
    }
}

auto
compute_prim_dg_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pf, const size_t pg,
                                     const size_t dd0, const size_t dd1, const size_t df,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;

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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_14 = buffer.data(pg + 14);
    const auto *pg_17 = buffer.data(pg + 17);
    const auto *pg_18 = buffer.data(pg + 18);
    const auto *pg_20 = buffer.data(pg + 20);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_3 = buffer.data(dd1 + 3);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_26 = buffer.data(df + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pf_0, dd0_0, dd1_0, df_0, \
                         df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = pb_y[k] * df_0[k];

        t_2[k] = pb_z[k] * df_0[k];

        t_3[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_4[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, dd0_1, dd0_2, dd1_1, dd1_2, \
                         df_3, df_5, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_3[k];

        t_6[k] = pb_z[k] * df_3[k];

        t_7[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_5[k];

        t_8[k] = pb_y[k] * df_6[k];

        t_9[k] = f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, pa_x, pa_y, pf_1, pg_0, pg_2, \
                         pg_6, pg_8, pg_9, pg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * pg_0[k];

        t_11[k] = f_0 * pf_1[k]
                  + pa_x[k] * pg_6[k];

        t_12[k] = pa_y[k] * pg_2[k];

        t_13[k] = pa_x[k] * pg_8[k];

        t_14[k] = pa_x[k] * pg_9[k];

        t_15[k] = pa_x[k] * pg_10[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_x, pa_z, pb_z, pf_0, pf_5, pg_0, \
                         pg_1, pg_14, pg_17, df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_z[k] * pg_0[k];

        t_17[k] = f_5 * pf_0[k]
                  + pb_z[k] * df_9[k];

        t_18[k] = pa_z[k] * pg_1[k];

        t_19[k] = f_0 * pf_5[k]
                  + pa_x[k] * pg_14[k];

        t_20[k] = pa_x[k] * pg_17[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, pg_18, pg_20, dd0_3, \
                         dd0_4, dd1_3, dd1_4, df_11, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_x[k] * pg_18[k];

        t_22[k] = pa_x[k] * pg_20[k];

        t_23[k] = f_1 * dd0_3[k]
                  - f_2 * dd1_3[k]
                  + pb_x[k] * df_11[k];

        t_24[k] = pb_z[k] * df_11[k];

        t_25[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_x[k] * df_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pb_x, pb_y, pb_z, pf_2, dd0_4, dd0_5, \
                         dd1_4, dd1_5, df_14, df_15, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * dd0_5[k]
                  - f_4 * dd1_5[k]
                  + pb_x[k] * df_14[k];

        t_27[k] = pb_x[k] * df_15[k];

        t_28[k] = pb_x[k] * df_17[k];

        t_29[k] = f_0 * pf_2[k]
                  + f_1 * dd0_4[k]
                  - f_2 * dd1_4[k]
                  + pb_y[k] * df_15[k];

        t_30[k] = pb_z[k] * df_15[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pb_y, pb_z, pf_3, pg_13, dd0_4, dd0_5, \
                         dd1_4, dd1_5, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_z[k] * df_16[k];

        t_32[k] = f_0 * pf_3[k]
                  + pb_y[k] * df_17[k];

        t_33[k] = f_1 * dd0_5[k]
                  - f_2 * dd1_5[k]
                  + pb_z[k] * df_17[k];

        t_34[k] = pa_y[k] * pg_13[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_y, pa_z, pb_z, pf_2, pf_7, pg_6, \
                         pg_8, pg_14, pg_18, df_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_z[k] * pg_6[k];

        t_36[k] = pa_y[k] * pg_14[k];

        t_37[k] = pa_z[k] * pg_8[k];

        t_38[k] = f_5 * pf_2[k]
                  + pb_z[k] * df_18[k];

        t_39[k] = f_0 * pf_7[k]
                  + pa_y[k] * pg_18[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pa_y, pb_x, pb_y, pb_z, pf_4, pf_8, \
                         pg_20, dd0_6, dd1_6, df_19, df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * pf_8[k]
                  + pb_y[k] * df_19[k];

        t_41[k] = pa_y[k] * pg_20[k];

        t_42[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_6[k]
                  + pb_x[k] * df_20[k];

        t_43[k] = pb_y[k] * df_20[k];

        t_44[k] = f_0 * pf_4[k]
                  + pb_z[k] * df_20[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pb_x, pb_y, dd0_7, dd0_8, dd1_7, dd1_8, \
                         df_22, df_23, df_24, df_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_3 * dd0_7[k]
                  - f_4 * dd1_7[k]
                  + pb_x[k] * df_22[k];

        t_46[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_x[k] * df_23[k];

        t_47[k] = pb_x[k] * df_24[k];

        t_48[k] = pb_x[k] * df_26[k];

        t_49[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_7[k]
                  + pb_y[k] * df_24[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pb_y, pb_z, pf_6, pf_8, dd0_8, dd1_8, df_24, \
                         df_25, df_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_0 * pf_6[k]
                  + pb_z[k] * df_24[k];

        t_51[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_y[k] * df_25[k];

        t_52[k] = pb_y[k] * df_26[k];

        t_53[k] = f_0 * pf_8[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_8[k]
                  + pb_z[k] * df_26[k];
    }
}

auto
compute_prim_dg_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pf, const size_t pg,
                                     const size_t dd0, const size_t dd1, const size_t df,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;

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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_3 = buffer.data(dd1 + 3);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pf_0, dd0_0, dd1_0, df_0, \
                         df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = pb_y[k] * df_0[k];

        t_2[k] = pb_z[k] * df_0[k];

        t_3[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_4[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, pg_0, dd0_1, dd0_2, dd1_1, \
                         dd1_2, df_3, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_3[k];

        t_6[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_4[k];

        t_7[k] = pb_y[k] * df_5[k];

        t_8[k] = f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_5[k];

        t_9[k] = pa_y[k] * pg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pa_z, pb_z, pf_0, pf_1, pf_5, \
                         pg_0, pg_1, pg_2, pg_3, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * pf_1[k]
                  + pa_x[k] * pg_1[k];

        t_11[k] = pa_x[k] * pg_2[k];

        t_12[k] = pa_z[k] * pg_0[k];

        t_13[k] = f_5 * pf_0[k]
                  + pb_z[k] * df_6[k];

        t_14[k] = f_0 * pf_5[k]
                  + pa_x[k] * pg_3[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pb_x, pg_5, dd0_3, dd0_4, dd0_5, dd1_3, \
                         dd1_4, dd1_5, df_7, df_8, df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_x[k] * pg_5[k];

        t_16[k] = f_1 * dd0_3[k]
                  - f_2 * dd1_3[k]
                  + pb_x[k] * df_7[k];

        t_17[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_x[k] * df_8[k];

        t_18[k] = f_3 * dd0_5[k]
                  - f_4 * dd1_5[k]
                  + pb_x[k] * df_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, pb_x, pb_y, pb_z, pf_2, pf_3, \
                         dd0_4, dd1_4, df_10, df_11, df_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_x[k] * df_10[k];

        t_20[k] = pb_x[k] * df_12[k];

        t_21[k] = f_0 * pf_2[k]
                  + f_1 * dd0_4[k]
                  - f_2 * dd1_4[k]
                  + pb_y[k] * df_10[k];

        t_22[k] = pb_z[k] * df_10[k];

        t_23[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_z[k] * df_11[k];

        t_24[k] = f_0 * pf_3[k]
                  + pb_y[k] * df_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, pa_z, pb_z, pf_2, pf_7, pg_2, pg_4, \
                         dd0_5, dd1_5, df_12, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * dd0_5[k]
                  - f_2 * dd1_5[k]
                  + pb_z[k] * df_12[k];

        t_26[k] = pa_z[k] * pg_2[k];

        t_27[k] = f_5 * pf_2[k]
                  + pb_z[k] * df_13[k];

        t_28[k] = f_0 * pf_7[k]
                  + pa_y[k] * pg_4[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, pf_4, pf_8, pg_5, \
                         dd0_6, dd1_6, df_14, df_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_5 * pf_8[k]
                  + pb_y[k] * df_14[k];

        t_30[k] = pa_y[k] * pg_5[k];

        t_31[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_6[k]
                  + pb_x[k] * df_15[k];

        t_32[k] = f_0 * pf_4[k]
                  + pb_z[k] * df_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pb_x, pb_y, dd0_7, dd0_8, dd1_7, dd1_8, \
                         df_16, df_17, df_18, df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * dd0_7[k]
                  - f_4 * dd1_7[k]
                  + pb_x[k] * df_16[k];

        t_34[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_x[k] * df_17[k];

        t_35[k] = pb_x[k] * df_18[k];

        t_36[k] = pb_x[k] * df_20[k];

        t_37[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_7[k]
                  + pb_y[k] * df_18[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pb_z, pf_6, pf_8, dd0_8, dd1_8, df_18, \
                         df_19, df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * pf_6[k]
                  + pb_z[k] * df_18[k];

        t_39[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_y[k] * df_19[k];

        t_40[k] = pb_y[k] * df_20[k];

        t_41[k] = f_0 * pf_8[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_8[k]
                  + pb_z[k] * df_20[k];
    }
}

auto
compute_prim_dg_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pf, const size_t pg,
                                     const size_t dd0, const size_t dd1, const size_t df,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_11 = buffer.data(pg + 11);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_3 = buffer.data(dd1 + 3);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pf_0, dd0_0, dd1_0, df_0, \
                         df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = pb_y[k] * df_0[k];

        t_2[k] = pb_z[k] * df_0[k];

        t_3[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_4[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, dd0_1, dd0_2, dd1_1, dd1_2, \
                         df_3, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_3[k];

        t_6[k] = pb_z[k] * df_3[k];

        t_7[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_4[k];

        t_8[k] = pb_y[k] * df_5[k];

        t_9[k] = f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pa_y, pa_z, pb_x, pg_0, pg_3, \
                         pg_11, dd0_3, dd1_3, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * pg_0[k];

        t_11[k] = pa_x[k] * pg_3[k];

        t_12[k] = pa_z[k] * pg_0[k];

        t_13[k] = pa_x[k] * pg_11[k];

        t_14[k] = f_1 * dd0_3[k]
                  - f_2 * dd1_3[k]
                  + pb_x[k] * df_6[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pb_x, pb_z, dd0_4, dd0_5, dd1_4, dd1_5, \
                         df_6, df_7, df_8, df_9, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_z[k] * df_6[k];

        t_16[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_x[k] * df_7[k];

        t_17[k] = f_3 * dd0_5[k]
                  - f_4 * dd1_5[k]
                  + pb_x[k] * df_8[k];

        t_18[k] = pb_x[k] * df_9[k];

        t_19[k] = pb_x[k] * df_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pb_y, pb_z, pf_1, dd0_4, dd0_5, dd1_4, dd1_5, \
                         df_9, df_10, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * pf_1[k]
                  + f_1 * dd0_4[k]
                  - f_2 * dd1_4[k]
                  + pb_y[k] * df_9[k];

        t_21[k] = pb_z[k] * df_9[k];

        t_22[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_z[k] * df_10[k];

        t_23[k] = f_1 * dd0_5[k]
                  - f_2 * dd1_5[k]
                  + pb_z[k] * df_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_y, pa_z, pb_x, pb_y, pf_2, pg_3, \
                         pg_9, pg_11, dd0_6, dd1_6, df_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_z[k] * pg_3[k];

        t_25[k] = f_0 * pf_2[k]
                  + pa_y[k] * pg_9[k];

        t_26[k] = pa_y[k] * pg_11[k];

        t_27[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_6[k]
                  + pb_x[k] * df_12[k];

        t_28[k] = pb_y[k] * df_12[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pb_x, pb_y, dd0_7, dd0_8, dd1_7, dd1_8, \
                         df_13, df_14, df_15, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_3 * dd0_7[k]
                  - f_4 * dd1_7[k]
                  + pb_x[k] * df_13[k];

        t_30[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_x[k] * df_14[k];

        t_31[k] = pb_x[k] * df_15[k];

        t_32[k] = pb_x[k] * df_17[k];

        t_33[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_7[k]
                  + pb_y[k] * df_15[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pb_y, pb_z, pf_3, dd0_8, dd1_8, df_16, \
                         df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_y[k] * df_16[k];

        t_35[k] = pb_y[k] * df_17[k];

        t_36[k] = f_0 * pf_3[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_8[k]
                  + pb_z[k] * df_17[k];
    }
}

auto
compute_prim_dg_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pf, const size_t pg,
                                     const size_t dd0, const size_t dd1, const size_t df,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_3 = buffer.data(dd1 + 3);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pf_0, dd0_0, dd1_0, df_0, \
                         df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = pb_y[k] * df_0[k];

        t_2[k] = pb_z[k] * df_0[k];

        t_3[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_4[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_z, pb_y, pb_z, pg_0, dd0_1, dd0_2, dd1_1, \
                         dd1_2, df_3, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_3[k];

        t_6[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_4[k];

        t_7[k] = pb_y[k] * df_5[k];

        t_8[k] = f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_5[k];

        t_9[k] = pa_z[k] * pg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_x, dd0_3, dd0_4, dd0_5, dd1_3, dd1_4, \
                         dd1_5, df_6, df_7, df_8, df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * dd0_3[k]
                  - f_2 * dd1_3[k]
                  + pb_x[k] * df_6[k];

        t_11[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_x[k] * df_7[k];

        t_12[k] = f_3 * dd0_5[k]
                  - f_4 * dd1_5[k]
                  + pb_x[k] * df_8[k];

        t_13[k] = pb_x[k] * df_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, pf_1, dd0_4, dd0_5, \
                         dd1_4, dd1_5, df_9, df_10, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pb_x[k] * df_11[k];

        t_15[k] = f_0 * pf_1[k]
                  + f_1 * dd0_4[k]
                  - f_2 * dd1_4[k]
                  + pb_y[k] * df_9[k];

        t_16[k] = pb_z[k] * df_9[k];

        t_17[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_z[k] * df_10[k];

        t_18[k] = f_1 * dd0_5[k]
                  - f_2 * dd1_5[k]
                  + pb_z[k] * df_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_x, pg_1, pg_2, dd0_6, dd0_7, \
                         dd1_6, dd1_7, df_12, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_z[k] * pg_1[k];

        t_20[k] = pa_y[k] * pg_2[k];

        t_21[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_6[k]
                  + pb_x[k] * df_12[k];

        t_22[k] = f_3 * dd0_7[k]
                  - f_4 * dd1_7[k]
                  + pb_x[k] * df_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, t_28, pb_x, pb_y, dd0_7, dd0_8, dd1_7, \
                         dd1_8, df_14, df_15, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_x[k] * df_14[k];

        t_24[k] = pb_x[k] * df_15[k];

        t_25[k] = pb_x[k] * df_17[k];

        t_26[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_7[k]
                  + pb_y[k] * df_15[k];

        t_27[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_y[k] * df_16[k];

        t_28[k] = pb_y[k] * df_17[k];
    }

#pragma omp simd aligned(t_29, pb_z, pf_2, dd0_8, dd1_8, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * pf_2[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_8[k]
                  + pb_z[k] * df_17[k];
    }
}

auto
compute_prim_dg_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pf, const size_t dd0, const size_t dd1,
                                     const size_t df, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_1, dd1_0, dd1_1, df_0, \
                         df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_0 * pf_1[k]
                 + f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pf_2, dd0_2, dd1_2, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pf_2[k]
                 + f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_2[k];
    }
}

auto
compute_prim_dg_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_11 = buffer.data(dd1 + 11);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_29 = buffer.data(df + 29);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_1, dd1_0, dd1_5, df_0, \
                         df_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_0 * pf_1[k]
                 + f_1 * dd0_1[k]
                 - f_2 * dd1_5[k]
                 + pb_y[k] * df_14[k];
    }

#pragma omp simd aligned(t_2, pb_z, pf_2, dd0_2, dd1_11, df_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pf_2[k]
                 + f_1 * dd0_2[k]
                 - f_2 * dd1_11[k]
                 + pb_z[k] * df_29[k];
    }
}

auto
compute_prim_dg_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_9 = buffer.data(dd0 + 9);
    const auto *dd0_10 = buffer.data(dd0 + 10);
    const auto *dd0_11 = buffer.data(dd0 + 11);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_9 = buffer.data(dd1 + 9);
    const auto *dd1_10 = buffer.data(dd1 + 10);
    const auto *dd1_11 = buffer.data(dd1 + 11);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pf_0, dd0_0, dd0_1, dd1_0, \
                         dd1_1, df_0, df_1, df_2, df_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_2[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];

        t_3[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, pb_y, pb_z, dd0_2, dd0_4, dd1_2, dd1_4, df_4, \
                         df_5, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_4[k];

        t_5[k] = f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_5[k];

        t_6[k] = f_1 * dd0_4[k]
                 - f_2 * dd1_4[k]
                 + pb_x[k] * df_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, pf_1, dd0_5, dd0_6, dd1_5, \
                         dd1_6, df_7, df_8, df_9, df_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * dd0_5[k]
                 - f_4 * dd1_5[k]
                 + pb_x[k] * df_7[k];

        t_8[k] = f_3 * dd0_6[k]
                 - f_4 * dd1_6[k]
                 + pb_x[k] * df_8[k];

        t_9[k] = f_0 * pf_1[k]
                 + f_1 * dd0_5[k]
                 - f_2 * dd1_5[k]
                 + pb_y[k] * df_9[k];

        t_10[k] = f_3 * dd0_5[k]
                  - f_4 * dd1_5[k]
                  + pb_z[k] * df_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pb_x, pb_z, dd0_6, dd0_9, dd0_10, dd1_6, dd1_9, \
                         dd1_10, df_11, df_12, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_6[k]
                  + pb_z[k] * df_11[k];

        t_12[k] = f_1 * dd0_9[k]
                  - f_2 * dd1_9[k]
                  + pb_x[k] * df_12[k];

        t_13[k] = f_3 * dd0_10[k]
                  - f_4 * dd1_10[k]
                  + pb_x[k] * df_13[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_x, pb_y, pb_z, pf_2, dd0_10, dd0_11, \
                         dd1_10, dd1_11, df_14, df_15, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * dd0_11[k]
                  - f_4 * dd1_11[k]
                  + pb_x[k] * df_14[k];

        t_15[k] = f_1 * dd0_10[k]
                  - f_2 * dd1_10[k]
                  + pb_y[k] * df_15[k];

        t_16[k] = f_3 * dd0_11[k]
                  - f_4 * dd1_11[k]
                  + pb_y[k] * df_16[k];

        t_17[k] = f_0 * pf_2[k]
                  + f_1 * dd0_11[k]
                  - f_2 * dd1_11[k]
                  + pb_z[k] * df_17[k];
    }
}

auto
compute_prim_dg_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_8 = buffer.data(dd1 + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_23 = buffer.data(df + 23);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_1, dd1_0, dd1_4, df_0, \
                         df_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_0 * pf_1[k]
                 + f_1 * dd0_1[k]
                 - f_2 * dd1_4[k]
                 + pb_y[k] * df_12[k];
    }

#pragma omp simd aligned(t_2, pb_z, pf_2, dd0_2, dd1_8, df_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pf_2[k]
                 + f_1 * dd0_2[k]
                 - f_2 * dd1_8[k]
                 + pb_z[k] * df_23[k];
    }
}

auto
compute_prim_dg_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t pf, const size_t pg,
                                      const size_t dd0, const size_t dd1, const size_t df,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;

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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_3 = buffer.data(dd1 + 3);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
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
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pf_0, dd0_0, dd1_0, df_0, \
                         df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = pb_y[k] * df_0[k];

        t_2[k] = pb_z[k] * df_0[k];

        t_3[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_4[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, pg_0, dd0_1, dd0_2, dd1_1, \
                         dd1_2, df_3, df_5, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_3[k];

        t_6[k] = pb_z[k] * df_3[k];

        t_7[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_5[k];

        t_8[k] = f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_6[k];

        t_9[k] = pa_y[k] * pg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pa_z, pb_z, pf_0, pf_1, pf_5, \
                         pg_0, pg_1, pg_2, pg_3, df_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * pf_1[k]
                  + pa_x[k] * pg_1[k];

        t_11[k] = pa_x[k] * pg_2[k];

        t_12[k] = pa_z[k] * pg_0[k];

        t_13[k] = f_5 * pf_0[k]
                  + pb_z[k] * df_7[k];

        t_14[k] = f_0 * pf_5[k]
                  + pa_x[k] * pg_3[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pb_x, pb_z, pg_5, dd0_3, dd0_4, dd1_3, \
                         dd1_4, df_8, df_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_x[k] * pg_5[k];

        t_16[k] = f_1 * dd0_3[k]
                  - f_2 * dd1_3[k]
                  + pb_x[k] * df_8[k];

        t_17[k] = pb_z[k] * df_8[k];

        t_18[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_x[k] * df_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pb_x, pb_y, pb_z, pf_2, dd0_4, dd0_5, \
                         dd1_4, dd1_5, df_11, df_12, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * dd0_5[k]
                  - f_4 * dd1_5[k]
                  + pb_x[k] * df_11[k];

        t_20[k] = pb_x[k] * df_12[k];

        t_21[k] = f_0 * pf_2[k]
                  + f_1 * dd0_4[k]
                  - f_2 * dd1_4[k]
                  + pb_y[k] * df_12[k];

        t_22[k] = pb_z[k] * df_12[k];

        t_23[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_z[k] * df_13[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_y, pb_z, pf_2, pf_3, pg_2, dd0_5, \
                         dd1_5, df_14, df_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * pf_3[k]
                  + pb_y[k] * df_14[k];

        t_25[k] = f_1 * dd0_5[k]
                  - f_2 * dd1_5[k]
                  + pb_z[k] * df_14[k];

        t_26[k] = pa_z[k] * pg_2[k];

        t_27[k] = f_5 * pf_2[k]
                  + pb_z[k] * df_15[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pf_7, pf_8, pg_4, \
                         pg_5, dd0_6, dd1_6, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * pf_7[k]
                  + pa_y[k] * pg_4[k];

        t_29[k] = f_5 * pf_8[k]
                  + pb_y[k] * df_16[k];

        t_30[k] = pa_y[k] * pg_5[k];

        t_31[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_6[k]
                  + pb_x[k] * df_17[k];

        t_32[k] = pb_y[k] * df_17[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_x, pb_z, pf_4, dd0_7, dd0_8, dd1_7, dd1_8, \
                         df_17, df_19, df_20, df_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * pf_4[k]
                  + pb_z[k] * df_17[k];

        t_34[k] = f_3 * dd0_7[k]
                  - f_4 * dd1_7[k]
                  + pb_x[k] * df_19[k];

        t_35[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_x[k] * df_20[k];

        t_36[k] = pb_x[k] * df_23[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pb_y, pb_z, pf_6, pf_8, dd0_7, dd0_8, \
                         dd1_7, dd1_8, df_21, df_22, df_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_7[k]
                  + pb_y[k] * df_21[k];

        t_38[k] = f_0 * pf_6[k]
                  + pb_z[k] * df_21[k];

        t_39[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_y[k] * df_22[k];

        t_40[k] = pb_y[k] * df_23[k];

        t_41[k] = f_0 * pf_8[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_8[k]
                  + pb_z[k] * df_23[k];
    }
}

auto
compute_prim_dg_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_8 = buffer.data(dd1 + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_17 = buffer.data(df + 17);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_1, dd1_0, dd1_4, df_0, \
                         df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_0 * pf_1[k]
                 + f_1 * dd0_1[k]
                 - f_2 * dd1_4[k]
                 + pb_y[k] * df_9[k];
    }

#pragma omp simd aligned(t_2, pb_z, pf_2, dd0_2, dd1_8, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pf_2[k]
                 + f_1 * dd0_2[k]
                 - f_2 * dd1_8[k]
                 + pb_z[k] * df_17[k];
    }
}

auto
compute_prim_dg_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t pf, const size_t pg,
                                      const size_t dd0, const size_t dd1, const size_t df,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_11 = buffer.data(pg + 11);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_3 = buffer.data(dd1 + 3);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pf_0, dd0_0, dd1_0, df_0, \
                         df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = pb_y[k] * df_0[k];

        t_2[k] = pb_z[k] * df_0[k];

        t_3[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_4[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, dd0_1, dd0_2, dd1_1, dd1_2, \
                         df_3, df_5, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_3[k];

        t_6[k] = pb_z[k] * df_3[k];

        t_7[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_5[k];

        t_8[k] = pb_y[k] * df_6[k];

        t_9[k] = f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_z, pb_x, pb_z, pg_0, dd0_3, dd0_4, dd1_3, \
                         dd1_4, df_7, df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * pg_0[k];

        t_11[k] = f_1 * dd0_3[k]
                  - f_2 * dd1_3[k]
                  + pb_x[k] * df_7[k];

        t_12[k] = pb_z[k] * df_7[k];

        t_13[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_x[k] * df_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, pf_1, dd0_4, dd0_5, \
                         dd1_4, dd1_5, df_10, df_11, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * dd0_5[k]
                  - f_4 * dd1_5[k]
                  + pb_x[k] * df_10[k];

        t_15[k] = pb_x[k] * df_11[k];

        t_16[k] = pb_x[k] * df_13[k];

        t_17[k] = f_0 * pf_1[k]
                  + f_1 * dd0_4[k]
                  - f_2 * dd1_4[k]
                  + pb_y[k] * df_11[k];

        t_18[k] = pb_z[k] * df_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_z, pg_3, pg_11, dd0_4, dd0_5, \
                         dd1_4, dd1_5, df_12, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_z[k] * df_12[k];

        t_20[k] = f_1 * dd0_5[k]
                  - f_2 * dd1_5[k]
                  + pb_z[k] * df_13[k];

        t_21[k] = pa_z[k] * pg_3[k];

        t_22[k] = pa_y[k] * pg_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_x, pb_y, dd0_6, dd0_7, dd0_8, dd1_6, \
                         dd1_7, dd1_8, df_14, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_6[k]
                  + pb_x[k] * df_14[k];

        t_24[k] = pb_y[k] * df_14[k];

        t_25[k] = f_3 * dd0_7[k]
                  - f_4 * dd1_7[k]
                  + pb_x[k] * df_16[k];

        t_26[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_x[k] * df_17[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pb_x, pb_y, dd0_7, dd0_8, dd1_7, dd1_8, \
                         df_18, df_19, df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_x[k] * df_18[k];

        t_28[k] = pb_x[k] * df_20[k];

        t_29[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_7[k]
                  + pb_y[k] * df_18[k];

        t_30[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_y[k] * df_19[k];

        t_31[k] = pb_y[k] * df_20[k];
    }

#pragma omp simd aligned(t_32, pb_z, pf_2, dd0_8, dd1_8, df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * pf_2[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_8[k]
                  + pb_z[k] * df_20[k];
    }
}

auto
compute_prim_dg_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t pf, const size_t pg,
                                      const size_t dd0, const size_t dd1, const size_t df,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_3 = buffer.data(dd1 + 3);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pf_0, dd0_0, dd1_0, df_0, \
                         df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = pb_y[k] * df_0[k];

        t_2[k] = pb_z[k] * df_0[k];

        t_3[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_4[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, dd0_1, dd0_2, dd1_1, dd1_2, \
                         df_3, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_3[k];

        t_6[k] = pb_z[k] * df_3[k];

        t_7[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_4[k];

        t_8[k] = pb_y[k] * df_5[k];

        t_9[k] = f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_z, pb_x, pb_z, pg_0, dd0_3, dd0_4, dd1_3, \
                         dd1_4, df_6, df_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * pg_0[k];

        t_11[k] = f_1 * dd0_3[k]
                  - f_2 * dd1_3[k]
                  + pb_x[k] * df_6[k];

        t_12[k] = pb_z[k] * df_6[k];

        t_13[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_x[k] * df_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, pf_1, dd0_4, dd0_5, \
                         dd1_4, dd1_5, df_8, df_9, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * dd0_5[k]
                  - f_4 * dd1_5[k]
                  + pb_x[k] * df_8[k];

        t_15[k] = pb_x[k] * df_9[k];

        t_16[k] = pb_x[k] * df_11[k];

        t_17[k] = f_0 * pf_1[k]
                  + f_1 * dd0_4[k]
                  - f_2 * dd1_4[k]
                  + pb_y[k] * df_9[k];

        t_18[k] = pb_z[k] * df_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_z, pg_1, pg_2, dd0_4, dd0_5, \
                         dd1_4, dd1_5, df_10, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_z[k] * df_10[k];

        t_20[k] = f_1 * dd0_5[k]
                  - f_2 * dd1_5[k]
                  + pb_z[k] * df_11[k];

        t_21[k] = pa_z[k] * pg_1[k];

        t_22[k] = pa_y[k] * pg_2[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_x, pb_y, dd0_6, dd0_7, dd0_8, dd1_6, \
                         dd1_7, dd1_8, df_12, df_13, df_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_6[k]
                  + pb_x[k] * df_12[k];

        t_24[k] = pb_y[k] * df_12[k];

        t_25[k] = f_3 * dd0_7[k]
                  - f_4 * dd1_7[k]
                  + pb_x[k] * df_13[k];

        t_26[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_x[k] * df_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pb_x, pb_y, dd0_7, dd0_8, dd1_7, dd1_8, \
                         df_15, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_x[k] * df_15[k];

        t_28[k] = pb_x[k] * df_17[k];

        t_29[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_7[k]
                  + pb_y[k] * df_15[k];

        t_30[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_y[k] * df_16[k];

        t_31[k] = pb_y[k] * df_17[k];
    }

#pragma omp simd aligned(t_32, pb_z, pf_2, dd0_8, dd1_8, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * pf_2[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_8[k]
                  + pb_z[k] * df_17[k];
    }
}

auto
compute_prim_dg_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_3 = buffer.data(dd1 + 3);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pf_0, dd0_0, dd1_0, df_0, \
                         df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = pb_y[k] * df_0[k];

        t_2[k] = pb_z[k] * df_0[k];

        t_3[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_4[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, dd0_1, dd0_2, dd1_1, dd1_2, \
                         df_3, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_3[k];

        t_6[k] = pb_z[k] * df_3[k];

        t_7[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_4[k];

        t_8[k] = pb_y[k] * df_5[k];

        t_9[k] = f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_x, pb_z, dd0_3, dd0_4, dd0_5, dd1_3, \
                         dd1_4, dd1_5, df_6, df_7, df_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * dd0_3[k]
                  - f_2 * dd1_3[k]
                  + pb_x[k] * df_6[k];

        t_11[k] = pb_z[k] * df_6[k];

        t_12[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_x[k] * df_7[k];

        t_13[k] = f_3 * dd0_5[k]
                  - f_4 * dd1_5[k]
                  + pb_x[k] * df_8[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, pf_1, dd0_4, dd1_4, \
                         df_9, df_10, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pb_x[k] * df_9[k];

        t_15[k] = pb_x[k] * df_11[k];

        t_16[k] = f_0 * pf_1[k]
                  + f_1 * dd0_4[k]
                  - f_2 * dd1_4[k]
                  + pb_y[k] * df_9[k];

        t_17[k] = pb_z[k] * df_9[k];

        t_18[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_z[k] * df_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_x, pb_y, pb_z, dd0_5, dd0_6, dd0_7, dd1_5, \
                         dd1_6, dd1_7, df_11, df_12, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * dd0_5[k]
                  - f_2 * dd1_5[k]
                  + pb_z[k] * df_11[k];

        t_20[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_6[k]
                  + pb_x[k] * df_12[k];

        t_21[k] = pb_y[k] * df_12[k];

        t_22[k] = f_3 * dd0_7[k]
                  - f_4 * dd1_7[k]
                  + pb_x[k] * df_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, t_28, pb_x, pb_y, dd0_7, dd0_8, dd1_7, \
                         dd1_8, df_14, df_15, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_x[k] * df_14[k];

        t_24[k] = pb_x[k] * df_15[k];

        t_25[k] = pb_x[k] * df_17[k];

        t_26[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_7[k]
                  + pb_y[k] * df_15[k];

        t_27[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_y[k] * df_16[k];

        t_28[k] = pb_y[k] * df_17[k];
    }

#pragma omp simd aligned(t_29, pb_z, pf_2, dd0_8, dd1_8, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * pf_2[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_8[k]
                  + pb_z[k] * df_17[k];
    }
}

auto
compute_prim_dg_electron_repulsion_18(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_11 = buffer.data(dd1 + 11);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_8 = buffer.data(df + 8);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_1, dd1_0, dd1_5, df_0, \
                         df_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_0 * pf_1[k]
                 + f_1 * dd0_1[k]
                 - f_2 * dd1_5[k]
                 + pb_y[k] * df_4[k];
    }

#pragma omp simd aligned(t_2, pb_z, pf_2, dd0_2, dd1_11, df_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pf_2[k]
                 + f_1 * dd0_2[k]
                 - f_2 * dd1_11[k]
                 + pb_z[k] * df_8[k];
    }
}

auto
compute_prim_dg_electron_repulsion_19(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_11 = buffer.data(dd0 + 11);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_8 = buffer.data(dd1 + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_5, dd1_0, dd1_4, df_0, \
                         df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_0 * pf_1[k]
                 + f_1 * dd0_5[k]
                 - f_2 * dd1_4[k]
                 + pb_y[k] * df_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pf_2, dd0_11, dd1_8, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pf_2[k]
                 + f_1 * dd0_11[k]
                 - f_2 * dd1_8[k]
                 + pb_z[k] * df_2[k];
    }
}

auto
compute_prim_dg_electron_repulsion_20(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_8 = buffer.data(dd1 + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_20 = buffer.data(df + 20);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_1, dd1_0, dd1_4, df_0, \
                         df_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_0 * pf_1[k]
                 + f_1 * dd0_1[k]
                 - f_2 * dd1_4[k]
                 + pb_y[k] * df_10[k];
    }

#pragma omp simd aligned(t_2, pb_z, pf_2, dd0_2, dd1_8, df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pf_2[k]
                 + f_1 * dd0_2[k]
                 - f_2 * dd1_8[k]
                 + pb_z[k] * df_20[k];
    }
}

auto
compute_prim_dg_electron_repulsion_21(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_3 = buffer.data(dd1 + 3);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pf_0, dd0_0, dd0_1, dd1_0, \
                         dd1_1, df_0, df_1, df_2, df_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_2[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];

        t_3[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, pb_y, pb_z, dd0_2, dd0_3, dd1_2, dd1_3, df_4, \
                         df_5, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_4[k];

        t_5[k] = f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_5[k];

        t_6[k] = f_1 * dd0_3[k]
                 - f_2 * dd1_3[k]
                 + pb_x[k] * df_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, pf_1, dd0_4, dd0_5, dd1_4, \
                         dd1_5, df_7, df_8, df_9, df_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * dd0_4[k]
                 - f_4 * dd1_4[k]
                 + pb_x[k] * df_7[k];

        t_8[k] = f_3 * dd0_5[k]
                 - f_4 * dd1_5[k]
                 + pb_x[k] * df_8[k];

        t_9[k] = f_0 * pf_1[k]
                 + f_1 * dd0_4[k]
                 - f_2 * dd1_4[k]
                 + pb_y[k] * df_9[k];

        t_10[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_z[k] * df_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pb_x, pb_z, dd0_5, dd0_6, dd0_7, dd1_5, dd1_6, \
                         dd1_7, df_11, df_12, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_1 * dd0_5[k]
                  - f_2 * dd1_5[k]
                  + pb_z[k] * df_11[k];

        t_12[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_6[k]
                  + pb_x[k] * df_12[k];

        t_13[k] = f_3 * dd0_7[k]
                  - f_4 * dd1_7[k]
                  + pb_x[k] * df_13[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_x, pb_y, pb_z, pf_2, dd0_7, dd0_8, dd1_7, \
                         dd1_8, df_14, df_15, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_x[k] * df_14[k];

        t_15[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_7[k]
                  + pb_y[k] * df_15[k];

        t_16[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_y[k] * df_16[k];

        t_17[k] = f_0 * pf_2[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_8[k]
                  + pb_z[k] * df_17[k];
    }
}

auto
compute_prim_dg_electron_repulsion_22(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_8 = buffer.data(dd1 + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_4, dd1_0, dd1_4, df_0, \
                         df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_0 * pf_1[k]
                 + f_1 * dd0_4[k]
                 - f_2 * dd1_4[k]
                 + pb_y[k] * df_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pf_2, dd0_8, dd1_8, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pf_2[k]
                 + f_1 * dd0_8[k]
                 - f_2 * dd1_8[k]
                 + pb_z[k] * df_2[k];
    }
}

auto
compute_prim_dg_electron_repulsion_23(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_8 = buffer.data(dd1 + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_20 = buffer.data(df + 20);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_1, dd1_0, dd1_4, df_0, \
                         df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_0 * pf_1[k]
                 + f_1 * dd0_1[k]
                 - f_2 * dd1_4[k]
                 + pb_y[k] * df_11[k];
    }

#pragma omp simd aligned(t_2, pb_z, pf_2, dd0_2, dd1_8, df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pf_2[k]
                 + f_1 * dd0_2[k]
                 - f_2 * dd1_8[k]
                 + pb_z[k] * df_20[k];
    }
}

auto
compute_prim_dg_electron_repulsion_24(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t pf, const size_t pg,
                                      const size_t dd0, const size_t dd1, const size_t df,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_3 = buffer.data(dd1 + 3);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pf_0, dd0_0, dd1_0, df_0, \
                         df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = pb_y[k] * df_0[k];

        t_2[k] = pb_z[k] * df_0[k];

        t_3[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_y[k] * df_1[k];

        t_4[k] = f_3 * dd0_0[k]
                 - f_4 * dd1_0[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_z, pb_y, pb_z, pg_0, dd0_1, dd0_2, dd1_1, \
                         dd1_2, df_3, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_3[k];

        t_6[k] = pb_z[k] * df_3[k];

        t_7[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_4[k];

        t_8[k] = f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_5[k];

        t_9[k] = pa_z[k] * pg_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_x, pb_z, dd0_3, dd0_4, dd0_5, dd1_3, \
                         dd1_4, dd1_5, df_6, df_7, df_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * dd0_3[k]
                  - f_2 * dd1_3[k]
                  + pb_x[k] * df_6[k];

        t_11[k] = pb_z[k] * df_6[k];

        t_12[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_x[k] * df_7[k];

        t_13[k] = f_3 * dd0_5[k]
                  - f_4 * dd1_5[k]
                  + pb_x[k] * df_8[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, pf_1, dd0_4, dd0_5, \
                         dd1_4, dd1_5, df_9, df_10, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pb_x[k] * df_9[k];

        t_15[k] = f_0 * pf_1[k]
                  + f_1 * dd0_4[k]
                  - f_2 * dd1_4[k]
                  + pb_y[k] * df_9[k];

        t_16[k] = pb_z[k] * df_9[k];

        t_17[k] = f_3 * dd0_4[k]
                  - f_4 * dd1_4[k]
                  + pb_z[k] * df_10[k];

        t_18[k] = f_1 * dd0_5[k]
                  - f_2 * dd1_5[k]
                  + pb_z[k] * df_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_x, pb_y, pg_1, pg_2, dd0_6, \
                         dd1_6, df_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_z[k] * pg_1[k];

        t_20[k] = pa_y[k] * pg_2[k];

        t_21[k] = f_1 * dd0_6[k]
                  - f_2 * dd1_6[k]
                  + pb_x[k] * df_12[k];

        t_22[k] = pb_y[k] * df_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pb_x, pb_y, dd0_7, dd0_8, dd1_7, dd1_8, \
                         df_13, df_14, df_15, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * dd0_7[k]
                  - f_4 * dd1_7[k]
                  + pb_x[k] * df_13[k];

        t_24[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_x[k] * df_14[k];

        t_25[k] = pb_x[k] * df_17[k];

        t_26[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_7[k]
                  + pb_y[k] * df_15[k];

        t_27[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_y[k] * df_16[k];
    }

#pragma omp simd aligned(t_28, t_29, pb_y, pb_z, pf_2, dd0_8, dd1_8, \
                         df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * df_17[k];

        t_29[k] = f_0 * pf_2[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_8[k]
                  + pb_z[k] * df_17[k];
    }
}

auto
compute_prim_dg_electron_repulsion_25(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_8 = buffer.data(dd1 + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_1, dd1_0, dd1_4, df_0, \
                         df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_0 * pf_1[k]
                 + f_1 * dd0_1[k]
                 - f_2 * dd1_4[k]
                 + pb_y[k] * df_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pf_2, dd0_2, dd1_8, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pf_2[k]
                 + f_1 * dd0_2[k]
                 - f_2 * dd1_8[k]
                 + pb_z[k] * df_2[k];
    }
}

auto
compute_prim_dg_electron_repulsion_26(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_4, dd1_0, dd1_1, df_0, \
                         df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_0 * pf_1[k]
                 + f_1 * dd0_4[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pf_2, dd0_8, dd1_2, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pf_2[k]
                 + f_1 * dd0_8[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_2[k];
    }
}

auto
compute_prim_dg_electron_repulsion_27(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_5 = buffer.data(df + 5);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_1, dd1_0, dd1_1, df_0, \
                         df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_0 * pf_1[k]
                 + f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_2[k];
    }

#pragma omp simd aligned(t_2, pb_z, pf_2, dd0_2, dd1_2, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pf_2[k]
                 + f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_5[k];
    }
}

auto
compute_prim_dg_electron_repulsion_28(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_1, dd1_0, dd1_1, \
                         df_0, df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_3 * dd0_1[k]
                 - f_4 * dd1_1[k]
                 + pb_x[k] * df_1[k];

        t_2[k] = f_0 * pf_1[k]
                 + f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, pb_z, pf_2, dd0_2, dd1_2, df_3, df_4, \
                         df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_x[k] * df_3[k];

        t_4[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_4[k];

        t_5[k] = f_0 * pf_2[k]
                 + f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_5[k];
    }
}

auto
compute_prim_dg_electron_repulsion_29(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_8 = buffer.data(dd1 + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_1, dd1_0, dd1_4, \
                         df_0, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_3 * dd0_1[k]
                 - f_4 * dd1_4[k]
                 + pb_x[k] * df_4[k];

        t_2[k] = f_0 * pf_1[k]
                 + f_1 * dd0_1[k]
                 - f_2 * dd1_4[k]
                 + pb_y[k] * df_5[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, pb_z, pf_2, dd0_2, dd1_8, df_8, df_10, \
                         df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_8[k]
                 + pb_x[k] * df_8[k];

        t_4[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_8[k]
                 + pb_y[k] * df_10[k];

        t_5[k] = f_0 * pf_2[k]
                 + f_1 * dd0_2[k]
                 - f_2 * dd1_8[k]
                 + pb_z[k] * df_11[k];
    }
}

auto
compute_prim_dg_electron_repulsion_30(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_8 = buffer.data(dd1 + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_4, dd1_0, dd1_4, \
                         df_0, df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_3 * dd0_4[k]
                 - f_4 * dd1_4[k]
                 + pb_x[k] * df_1[k];

        t_2[k] = f_0 * pf_1[k]
                 + f_1 * dd0_4[k]
                 - f_2 * dd1_4[k]
                 + pb_y[k] * df_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, pb_z, pf_2, dd0_8, dd1_8, df_3, df_4, \
                         df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * dd0_8[k]
                 - f_4 * dd1_8[k]
                 + pb_x[k] * df_3[k];

        t_4[k] = f_3 * dd0_8[k]
                 - f_4 * dd1_8[k]
                 + pb_y[k] * df_4[k];

        t_5[k] = f_0 * pf_2[k]
                 + f_1 * dd0_8[k]
                 - f_2 * dd1_8[k]
                 + pb_z[k] * df_5[k];
    }
}

auto
compute_prim_dg_electron_repulsion_31(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_1, dd1_0, \
                         dd1_1, df_0, df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_3 * dd0_1[k]
                 - f_4 * dd1_1[k]
                 + pb_x[k] * df_1[k];

        t_2[k] = pb_x[k] * df_2[k];

        t_3[k] = f_0 * pf_1[k]
                 + f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_2[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, pf_2, dd0_2, dd1_2, df_3, \
                         df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_x[k] * df_3[k];

        t_5[k] = pb_x[k] * df_5[k];

        t_6[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_2[k]
                 + pb_y[k] * df_4[k];

        t_7[k] = pb_y[k] * df_5[k];

        t_8[k] = f_0 * pf_2[k]
                 + f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_5[k];
    }
}

auto
compute_prim_dg_electron_repulsion_32(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_8 = buffer.data(dd1 + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_1, dd1_0, \
                         dd1_4, df_0, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_3 * dd0_1[k]
                 - f_4 * dd1_4[k]
                 + pb_x[k] * df_4[k];

        t_2[k] = pb_x[k] * df_5[k];

        t_3[k] = f_0 * pf_1[k]
                 + f_1 * dd0_1[k]
                 - f_2 * dd1_4[k]
                 + pb_y[k] * df_5[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, pf_2, dd0_2, dd1_8, df_8, \
                         df_10, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_8[k]
                 + pb_x[k] * df_8[k];

        t_5[k] = pb_x[k] * df_11[k];

        t_6[k] = f_3 * dd0_2[k]
                 - f_4 * dd1_8[k]
                 + pb_y[k] * df_10[k];

        t_7[k] = pb_y[k] * df_11[k];

        t_8[k] = f_0 * pf_2[k]
                 + f_1 * dd0_2[k]
                 - f_2 * dd1_8[k]
                 + pb_z[k] * df_11[k];
    }
}

auto
compute_prim_dg_electron_repulsion_33(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_1 = buffer.data(dd0 + 1);
    const auto *dd0_2 = buffer.data(dd0 + 2);
    const auto *dd0_3 = buffer.data(dd0 + 3);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_5 = buffer.data(dd0 + 5);
    const auto *dd0_6 = buffer.data(dd0 + 6);
    const auto *dd0_7 = buffer.data(dd0 + 7);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_1 = buffer.data(dd1 + 1);
    const auto *dd1_2 = buffer.data(dd1 + 2);
    const auto *dd1_3 = buffer.data(dd1 + 3);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_5 = buffer.data(dd1 + 5);
    const auto *dd1_6 = buffer.data(dd1 + 6);
    const auto *dd1_7 = buffer.data(dd1 + 7);
    const auto *dd1_8 = buffer.data(dd1 + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pf_0, dd0_0, dd0_1, dd0_2, dd1_0, \
                         dd1_1, dd1_2, df_0, df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_1 * dd0_1[k]
                 - f_2 * dd1_1[k]
                 + pb_y[k] * df_1[k];

        t_2[k] = f_1 * dd0_2[k]
                 - f_2 * dd1_2[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pb_x, pb_y, pf_1, dd0_3, dd0_4, dd1_3, dd1_4, \
                         df_3, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_1 * dd0_3[k]
                 - f_2 * dd1_3[k]
                 + pb_x[k] * df_3[k];

        t_4[k] = f_3 * dd0_4[k]
                 - f_4 * dd1_4[k]
                 + pb_x[k] * df_4[k];

        t_5[k] = pb_x[k] * df_5[k];

        t_6[k] = f_0 * pf_1[k]
                 + f_1 * dd0_4[k]
                 - f_2 * dd1_4[k]
                 + pb_y[k] * df_5[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pb_x, pb_z, dd0_5, dd0_6, dd0_8, dd1_5, dd1_6, \
                         dd1_8, df_6, df_7, df_8, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * dd0_5[k]
                 - f_2 * dd1_5[k]
                 + pb_z[k] * df_6[k];

        t_8[k] = f_1 * dd0_6[k]
                 - f_2 * dd1_6[k]
                 + pb_x[k] * df_7[k];

        t_9[k] = f_3 * dd0_8[k]
                 - f_4 * dd1_8[k]
                 + pb_x[k] * df_8[k];

        t_10[k] = pb_x[k] * df_11[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pb_y, pb_z, pf_2, dd0_7, dd0_8, dd1_7, dd1_8, \
                         df_9, df_10, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_1 * dd0_7[k]
                  - f_2 * dd1_7[k]
                  + pb_y[k] * df_9[k];

        t_12[k] = f_3 * dd0_8[k]
                  - f_4 * dd1_8[k]
                  + pb_y[k] * df_10[k];

        t_13[k] = pb_y[k] * df_11[k];

        t_14[k] = f_0 * pf_2[k]
                  + f_1 * dd0_8[k]
                  - f_2 * dd1_8[k]
                  + pb_z[k] * df_11[k];
    }
}

auto
compute_prim_dg_electron_repulsion_34(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pf, const size_t dd0, const size_t dd1,
                                      const size_t df, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd0_0 = buffer.data(dd0 + 0);
    const auto *dd0_4 = buffer.data(dd0 + 4);
    const auto *dd0_8 = buffer.data(dd0 + 8);

    const auto *dd1_0 = buffer.data(dd1 + 0);
    const auto *dd1_4 = buffer.data(dd1 + 4);
    const auto *dd1_8 = buffer.data(dd1 + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pf_0, pf_1, dd0_0, dd0_4, dd1_0, \
                         dd1_4, df_0, df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 + f_1 * dd0_0[k]
                 - f_2 * dd1_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_3 * dd0_4[k]
                 - f_4 * dd1_4[k]
                 + pb_x[k] * df_1[k];

        t_2[k] = pb_x[k] * df_2[k];

        t_3[k] = f_0 * pf_1[k]
                 + f_1 * dd0_4[k]
                 - f_2 * dd1_4[k]
                 + pb_y[k] * df_2[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, pf_2, dd0_8, dd1_8, df_3, \
                         df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * dd0_8[k]
                 - f_4 * dd1_8[k]
                 + pb_x[k] * df_3[k];

        t_5[k] = pb_x[k] * df_5[k];

        t_6[k] = f_3 * dd0_8[k]
                 - f_4 * dd1_8[k]
                 + pb_y[k] * df_4[k];

        t_7[k] = pb_y[k] * df_5[k];

        t_8[k] = f_0 * pf_2[k]
                 + f_1 * dd0_8[k]
                 - f_2 * dd1_8[k]
                 + pb_z[k] * df_5[k];
    }
}

}  // namespace simdt2ceri
