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


#include "SimdElectronRepulsionVrrRecDH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_dh_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pg, const size_t ph,
                                     const size_t df0, const size_t df1, const size_t dg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.5 / p;

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
    const auto *pg_18 = buffer.data(pg + 18);
    const auto *pg_19 = buffer.data(pg + 19);
    const auto *pg_20 = buffer.data(pg + 20);
    const auto *pg_21 = buffer.data(pg + 21);
    const auto *pg_22 = buffer.data(pg + 22);
    const auto *pg_23 = buffer.data(pg + 23);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_8 = buffer.data(ph + 8);
    const auto *ph_9 = buffer.data(ph + 9);
    const auto *ph_10 = buffer.data(ph + 10);
    const auto *ph_11 = buffer.data(ph + 11);
    const auto *ph_12 = buffer.data(ph + 12);
    const auto *ph_13 = buffer.data(ph + 13);
    const auto *ph_14 = buffer.data(ph + 14);
    const auto *ph_15 = buffer.data(ph + 15);
    const auto *ph_16 = buffer.data(ph + 16);
    const auto *ph_17 = buffer.data(ph + 17);
    const auto *ph_18 = buffer.data(ph + 18);
    const auto *ph_19 = buffer.data(ph + 19);
    const auto *ph_20 = buffer.data(ph + 20);
    const auto *ph_21 = buffer.data(ph + 21);
    const auto *ph_22 = buffer.data(ph + 22);
    const auto *ph_23 = buffer.data(ph + 23);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);
    const auto *df0_3 = buffer.data(df0 + 3);
    const auto *df0_4 = buffer.data(df0 + 4);
    const auto *df0_5 = buffer.data(df0 + 5);
    const auto *df0_6 = buffer.data(df0 + 6);
    const auto *df0_7 = buffer.data(df0 + 7);
    const auto *df0_8 = buffer.data(df0 + 8);
    const auto *df0_9 = buffer.data(df0 + 9);
    const auto *df0_10 = buffer.data(df0 + 10);
    const auto *df0_11 = buffer.data(df0 + 11);
    const auto *df0_12 = buffer.data(df0 + 12);
    const auto *df0_13 = buffer.data(df0 + 13);
    const auto *df0_14 = buffer.data(df0 + 14);
    const auto *df0_15 = buffer.data(df0 + 15);
    const auto *df0_16 = buffer.data(df0 + 16);
    const auto *df0_17 = buffer.data(df0 + 17);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_1 = buffer.data(df1 + 1);
    const auto *df1_2 = buffer.data(df1 + 2);
    const auto *df1_3 = buffer.data(df1 + 3);
    const auto *df1_4 = buffer.data(df1 + 4);
    const auto *df1_5 = buffer.data(df1 + 5);
    const auto *df1_6 = buffer.data(df1 + 6);
    const auto *df1_7 = buffer.data(df1 + 7);
    const auto *df1_8 = buffer.data(df1 + 8);
    const auto *df1_9 = buffer.data(df1 + 9);
    const auto *df1_10 = buffer.data(df1 + 10);
    const auto *df1_11 = buffer.data(df1 + 11);
    const auto *df1_12 = buffer.data(df1 + 12);
    const auto *df1_13 = buffer.data(df1 + 13);
    const auto *df1_14 = buffer.data(df1 + 14);
    const auto *df1_15 = buffer.data(df1 + 15);
    const auto *df1_16 = buffer.data(df1 + 16);
    const auto *df1_17 = buffer.data(df1 + 17);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_40 = buffer.data(dg + 40);
    const auto *dg_41 = buffer.data(dg + 41);
    const auto *dg_42 = buffer.data(dg + 42);
    const auto *dg_43 = buffer.data(dg + 43);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_46 = buffer.data(dg + 46);
    const auto *dg_47 = buffer.data(dg + 47);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_49 = buffer.data(dg + 49);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_51 = buffer.data(dg + 51);
    const auto *dg_52 = buffer.data(dg + 52);
    const auto *dg_53 = buffer.data(dg + 53);
    const auto *dg_54 = buffer.data(dg + 54);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_56 = buffer.data(dg + 56);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, pg_0, df0_0, df1_0, \
                         dg_0, dg_1, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = pb_y[k] * dg_0[k];

        t_2[k] = pb_z[k] * dg_0[k];

        t_3[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_y[k] * dg_1[k];

        t_4[k] = pb_y[k] * dg_2[k];

        t_5[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_z[k] * dg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, pg_3, df0_1, df0_2, \
                         df1_1, df1_2, dg_3, dg_4, dg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * df0_1[k]
                 - f_6 * df1_1[k]
                 + pb_y[k] * dg_3[k];

        t_7[k] = pb_z[k] * dg_3[k];

        t_8[k] = pb_y[k] * dg_4[k];

        t_9[k] = f_5 * df0_2[k]
                 - f_6 * df1_2[k]
                 + pb_z[k] * dg_4[k];

        t_10[k] = f_0 * pg_3[k]
                  + pb_x[k] * dg_7[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pb_x, pb_y, pb_z, pg_4, pg_5, dg_5, dg_6, \
                         dg_8, dg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * dg_5[k];

        t_12[k] = f_0 * pg_4[k]
                  + pb_x[k] * dg_8[k];

        t_13[k] = pb_y[k] * dg_6[k];

        t_14[k] = f_0 * pg_5[k]
                  + pb_x[k] * dg_10[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pb_z, df0_3, df0_4, df0_5, df1_3, \
                         df1_4, df1_5, dg_7, dg_8, dg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * df0_3[k]
                  - f_2 * df1_3[k]
                  + pb_y[k] * dg_7[k];

        t_16[k] = pb_z[k] * dg_7[k];

        t_17[k] = f_5 * df0_4[k]
                  - f_6 * df1_4[k]
                  + pb_y[k] * dg_8[k];

        t_18[k] = f_3 * df0_5[k]
                  - f_4 * df1_5[k]
                  + pb_y[k] * dg_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pb_y, pb_z, pg_0, ph_0, df0_5, \
                         df1_5, dg_10, dg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * dg_10[k];

        t_20[k] = f_1 * df0_5[k]
                  - f_2 * df1_5[k]
                  + pb_z[k] * dg_10[k];

        t_21[k] = pa_y[k] * ph_0[k];

        t_22[k] = f_7 * pg_0[k]
                  + pb_y[k] * dg_11[k];

        t_23[k] = pb_z[k] * dg_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_x, pa_y, pb_z, pg_7, pg_9, ph_2, \
                         ph_8, ph_9, dg_12, dg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * pg_7[k]
                  + pa_x[k] * ph_8[k];

        t_25[k] = pb_z[k] * dg_12[k];

        t_26[k] = pa_y[k] * ph_2[k];

        t_27[k] = f_0 * pg_9[k]
                  + pa_x[k] * ph_9[k];

        t_28[k] = pb_z[k] * dg_13[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, pg_2, pg_10, ph_4, \
                         dg_14, dg_15, dg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_7 * pg_2[k]
                  + pb_y[k] * dg_14[k];

        t_30[k] = pa_y[k] * ph_4[k];

        t_31[k] = f_7 * pg_10[k]
                  + pb_x[k] * dg_16[k];

        t_32[k] = pb_z[k] * dg_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_x, pa_y, pb_x, pb_z, pg_11, pg_12, \
                         ph_6, ph_10, dg_16, dg_17, dg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_7 * pg_11[k]
                  + pb_x[k] * dg_17[k];

        t_34[k] = f_7 * pg_12[k]
                  + pb_x[k] * dg_18[k];

        t_35[k] = pa_y[k] * ph_6[k];

        t_36[k] = pa_x[k] * ph_10[k];

        t_37[k] = pb_z[k] * dg_16[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, pa_x, pa_z, pb_y, ph_0, ph_11, \
                         ph_12, ph_13, ph_14, dg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_x[k] * ph_11[k];

        t_39[k] = pa_x[k] * ph_12[k];

        t_40[k] = pa_x[k] * ph_13[k];

        t_41[k] = pa_x[k] * ph_14[k];

        t_42[k] = pa_z[k] * ph_0[k];

        t_43[k] = pb_y[k] * dg_19[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pa_x, pa_z, pb_y, pb_z, pg_0, pg_17, \
                         ph_1, ph_3, ph_17, dg_19, dg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_7 * pg_0[k]
                  + pb_z[k] * dg_19[k];

        t_45[k] = pa_z[k] * ph_1[k];

        t_46[k] = pb_y[k] * dg_20[k];

        t_47[k] = f_8 * pg_17[k]
                  + pa_x[k] * ph_17[k];

        t_48[k] = pa_z[k] * ph_3[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_x, pa_z, pb_y, pb_z, pg_1, pg_18, ph_5, \
                         ph_18, dg_21, dg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_7 * pg_1[k]
                  + pb_z[k] * dg_21[k];

        t_50[k] = pb_y[k] * dg_22[k];

        t_51[k] = f_0 * pg_18[k]
                  + pa_x[k] * ph_18[k];

        t_52[k] = pa_z[k] * ph_5[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_x, pb_x, pb_y, pg_20, pg_21, pg_23, \
                         ph_19, dg_23, dg_24, dg_25, dg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_7 * pg_20[k]
                  + pb_x[k] * dg_24[k];

        t_54[k] = f_7 * pg_21[k]
                  + pb_x[k] * dg_25[k];

        t_55[k] = pb_y[k] * dg_23[k];

        t_56[k] = f_7 * pg_23[k]
                  + pb_x[k] * dg_26[k];

        t_57[k] = pa_x[k] * ph_19[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_x, pb_y, ph_20, ph_21, ph_22, ph_23, \
                         dg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pa_x[k] * ph_20[k];

        t_59[k] = pa_x[k] * ph_21[k];

        t_60[k] = pa_x[k] * ph_22[k];

        t_61[k] = pb_y[k] * dg_26[k];

        t_62[k] = pa_x[k] * ph_23[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pb_x, pb_y, pb_z, pg_6, df0_6, df0_7, \
                         df1_6, df1_7, dg_27, dg_28, dg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_1 * df0_6[k]
                  - f_2 * df1_6[k]
                  + pb_x[k] * dg_27[k];

        t_64[k] = f_0 * pg_6[k]
                  + pb_y[k] * dg_27[k];

        t_65[k] = pb_z[k] * dg_27[k];

        t_66[k] = f_5 * df0_7[k]
                  - f_6 * df1_7[k]
                  + pb_x[k] * dg_29[k];

        t_67[k] = pb_z[k] * dg_28[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pb_y, pb_z, pg_8, df0_8, df0_9, df1_8, \
                         df1_9, dg_29, dg_30, dg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_5 * df0_8[k]
                  - f_6 * df1_8[k]
                  + pb_x[k] * dg_30[k];

        t_69[k] = f_3 * df0_9[k]
                  - f_4 * df1_9[k]
                  + pb_x[k] * dg_31[k];

        t_70[k] = pb_z[k] * dg_29[k];

        t_71[k] = f_0 * pg_8[k]
                  + pb_y[k] * dg_30[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, pb_x, df0_11, df1_11, dg_32, \
                         dg_33, dg_34, dg_35, dg_36, dg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_3 * df0_11[k]
                  - f_4 * df1_11[k]
                  + pb_x[k] * dg_32[k];

        t_73[k] = pb_x[k] * dg_33[k];

        t_74[k] = pb_x[k] * dg_34[k];

        t_75[k] = pb_x[k] * dg_35[k];

        t_76[k] = pb_x[k] * dg_36[k];

        t_77[k] = pb_x[k] * dg_37[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_y, pb_z, pg_10, df0_9, df0_10, df1_9, \
                         df1_10, dg_33, dg_34, dg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * pg_10[k]
                  + f_1 * df0_9[k]
                  - f_2 * df1_9[k]
                  + pb_y[k] * dg_33[k];

        t_79[k] = pb_z[k] * dg_33[k];

        t_80[k] = f_3 * df0_9[k]
                  - f_4 * df1_9[k]
                  + pb_z[k] * dg_34[k];

        t_81[k] = f_5 * df0_10[k]
                  - f_6 * df1_10[k]
                  + pb_z[k] * dg_35[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pa_y, pa_z, pb_y, pb_z, pg_13, ph_7, \
                         ph_15, ph_16, df0_11, df1_11, dg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * pg_13[k]
                  + pb_y[k] * dg_37[k];

        t_83[k] = f_1 * df0_11[k]
                  - f_2 * df1_11[k]
                  + pb_z[k] * dg_37[k];

        t_84[k] = pa_y[k] * ph_15[k];

        t_85[k] = pa_z[k] * ph_7[k];

        t_86[k] = pa_y[k] * ph_16[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pa_y, pa_z, pb_y, pb_z, pg_7, pg_15, \
                         ph_8, ph_9, ph_17, dg_38, dg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * ph_8[k];

        t_88[k] = f_7 * pg_15[k]
                  + pb_y[k] * dg_38[k];

        t_89[k] = pa_y[k] * ph_17[k];

        t_90[k] = pa_z[k] * ph_9[k];

        t_91[k] = f_7 * pg_7[k]
                  + pb_z[k] * dg_39[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, t_97, pa_y, pb_x, pb_y, pg_17, ph_18, \
                         dg_40, dg_41, dg_42, dg_43, dg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_7 * pg_17[k]
                  + pb_y[k] * dg_40[k];

        t_93[k] = pa_y[k] * ph_18[k];

        t_94[k] = pb_x[k] * dg_41[k];

        t_95[k] = pb_x[k] * dg_42[k];

        t_96[k] = pb_x[k] * dg_43[k];

        t_97[k] = pb_x[k] * dg_44[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pa_y, pa_z, pb_x, pb_z, pg_10, pg_21, \
                         ph_10, ph_21, dg_41, dg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = pb_x[k] * dg_45[k];

        t_99[k] = pa_z[k] * ph_10[k];

        t_100[k] = f_7 * pg_10[k]
                   + pb_z[k] * dg_41[k];

        t_101[k] = f_8 * pg_21[k]
                   + pa_y[k] * ph_21[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, pa_y, pb_x, pb_y, pg_22, pg_23, \
                         ph_22, ph_23, df0_12, df1_12, dg_45, dg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_0 * pg_22[k]
                   + pa_y[k] * ph_22[k];

        t_103[k] = f_7 * pg_23[k]
                   + pb_y[k] * dg_45[k];

        t_104[k] = pa_y[k] * ph_23[k];

        t_105[k] = f_1 * df0_12[k]
                   - f_2 * df1_12[k]
                   + pb_x[k] * dg_46[k];

        t_106[k] = pb_y[k] * dg_46[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pb_x, pb_y, pb_z, pg_14, df0_13, df0_14, \
                         df1_13, df1_14, dg_46, dg_47, dg_48, dg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_0 * pg_14[k]
                   + pb_z[k] * dg_46[k];

        t_108[k] = f_5 * df0_13[k]
                   - f_6 * df1_13[k]
                   + pb_x[k] * dg_48[k];

        t_109[k] = pb_y[k] * dg_47[k];

        t_110[k] = f_5 * df0_14[k]
                   - f_6 * df1_14[k]
                   + pb_x[k] * dg_49[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pb_x, pb_y, pb_z, pg_16, df0_15, df0_17, \
                         df1_15, df1_17, dg_48, dg_49, dg_50, dg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_3 * df0_15[k]
                   - f_4 * df1_15[k]
                   + pb_x[k] * dg_50[k];

        t_112[k] = f_0 * pg_16[k]
                   + pb_z[k] * dg_48[k];

        t_113[k] = pb_y[k] * dg_49[k];

        t_114[k] = f_3 * df0_17[k]
                   - f_4 * df1_17[k]
                   + pb_x[k] * dg_51[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, t_120, pb_x, pb_y, df0_15, df1_15, \
                         dg_52, dg_53, dg_54, dg_55, dg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pb_x[k] * dg_52[k];

        t_116[k] = pb_x[k] * dg_53[k];

        t_117[k] = pb_x[k] * dg_54[k];

        t_118[k] = pb_x[k] * dg_55[k];

        t_119[k] = pb_x[k] * dg_56[k];

        t_120[k] = f_1 * df0_15[k]
                   - f_2 * df1_15[k]
                   + pb_y[k] * dg_52[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_y, pb_z, pg_19, df0_16, df0_17, \
                         df1_16, df1_17, dg_52, dg_54, dg_55, dg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_0 * pg_19[k]
                   + pb_z[k] * dg_52[k];

        t_122[k] = f_5 * df0_16[k]
                   - f_6 * df1_16[k]
                   + pb_y[k] * dg_54[k];

        t_123[k] = f_3 * df0_17[k]
                   - f_4 * df1_17[k]
                   + pb_y[k] * dg_55[k];

        t_124[k] = pb_y[k] * dg_56[k];
    }

#pragma omp simd aligned(t_125, pb_z, pg_23, df0_17, df1_17, dg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_0 * pg_23[k]
                   + f_1 * df0_17[k]
                   - f_2 * df1_17[k]
                   + pb_z[k] * dg_56[k];
    }
}

auto
compute_prim_dh_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pg, const size_t ph,
                                     const size_t df0, const size_t df1, const size_t dg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_11 = buffer.data(pg + 11);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_14 = buffer.data(pg + 14);
    const auto *pg_15 = buffer.data(pg + 15);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_8 = buffer.data(ph + 8);
    const auto *ph_9 = buffer.data(ph + 9);
    const auto *ph_10 = buffer.data(ph + 10);
    const auto *ph_11 = buffer.data(ph + 11);
    const auto *ph_12 = buffer.data(ph + 12);
    const auto *ph_13 = buffer.data(ph + 13);
    const auto *ph_16 = buffer.data(ph + 16);
    const auto *ph_17 = buffer.data(ph + 17);
    const auto *ph_18 = buffer.data(ph + 18);
    const auto *ph_20 = buffer.data(ph + 20);
    const auto *ph_21 = buffer.data(ph + 21);
    const auto *ph_22 = buffer.data(ph + 22);
    const auto *ph_23 = buffer.data(ph + 23);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);
    const auto *df0_3 = buffer.data(df0 + 3);
    const auto *df0_4 = buffer.data(df0 + 4);
    const auto *df0_5 = buffer.data(df0 + 5);
    const auto *df0_6 = buffer.data(df0 + 6);
    const auto *df0_7 = buffer.data(df0 + 7);
    const auto *df0_8 = buffer.data(df0 + 8);
    const auto *df0_9 = buffer.data(df0 + 9);
    const auto *df0_10 = buffer.data(df0 + 10);
    const auto *df0_11 = buffer.data(df0 + 11);
    const auto *df0_12 = buffer.data(df0 + 12);
    const auto *df0_13 = buffer.data(df0 + 13);
    const auto *df0_14 = buffer.data(df0 + 14);
    const auto *df0_15 = buffer.data(df0 + 15);
    const auto *df0_16 = buffer.data(df0 + 16);
    const auto *df0_17 = buffer.data(df0 + 17);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_1 = buffer.data(df1 + 1);
    const auto *df1_2 = buffer.data(df1 + 2);
    const auto *df1_3 = buffer.data(df1 + 3);
    const auto *df1_4 = buffer.data(df1 + 4);
    const auto *df1_5 = buffer.data(df1 + 5);
    const auto *df1_6 = buffer.data(df1 + 6);
    const auto *df1_7 = buffer.data(df1 + 7);
    const auto *df1_8 = buffer.data(df1 + 8);
    const auto *df1_9 = buffer.data(df1 + 9);
    const auto *df1_10 = buffer.data(df1 + 10);
    const auto *df1_11 = buffer.data(df1 + 11);
    const auto *df1_12 = buffer.data(df1 + 12);
    const auto *df1_13 = buffer.data(df1 + 13);
    const auto *df1_14 = buffer.data(df1 + 14);
    const auto *df1_15 = buffer.data(df1 + 15);
    const auto *df1_16 = buffer.data(df1 + 16);
    const auto *df1_17 = buffer.data(df1 + 17);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pg_0, df0_0, df1_0, dg_0, \
                         dg_1, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = pb_y[k] * dg_0[k];

        t_2[k] = pb_z[k] * dg_0[k];

        t_3[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_y[k] * dg_1[k];

        t_4[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_z[k] * dg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, pg_1, df0_1, df0_2, df1_1, \
                         df1_2, dg_3, dg_4, dg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * df0_1[k]
                 - f_6 * df1_1[k]
                 + pb_y[k] * dg_3[k];

        t_6[k] = pb_y[k] * dg_4[k];

        t_7[k] = f_5 * df0_2[k]
                 - f_6 * df1_2[k]
                 + pb_z[k] * dg_4[k];

        t_8[k] = f_0 * pg_1[k]
                 + pb_x[k] * dg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, pb_y, pg_2, df0_3, df0_4, df1_3, df1_4, dg_5, \
                         dg_6, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * pg_2[k]
                 + pb_x[k] * dg_8[k];

        t_10[k] = f_1 * df0_3[k]
                  - f_2 * df1_3[k]
                  + pb_y[k] * dg_5[k];

        t_11[k] = f_5 * df0_4[k]
                  - f_6 * df1_4[k]
                  + pb_y[k] * dg_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pa_y, pb_y, pb_z, pg_0, ph_0, df0_5, \
                         df1_5, dg_7, dg_8, dg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * df0_5[k]
                  - f_4 * df1_5[k]
                  + pb_y[k] * dg_7[k];

        t_13[k] = pb_y[k] * dg_8[k];

        t_14[k] = f_1 * df0_5[k]
                  - f_2 * df1_5[k]
                  + pb_z[k] * dg_8[k];

        t_15[k] = pa_y[k] * ph_0[k];

        t_16[k] = f_7 * pg_0[k]
                  + pb_y[k] * dg_9[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_x, pa_y, pb_x, pg_4, pg_5, pg_6, \
                         ph_2, ph_4, ph_8, ph_9, dg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * pg_4[k]
                  + pa_x[k] * ph_8[k];

        t_18[k] = pa_y[k] * ph_2[k];

        t_19[k] = f_0 * pg_5[k]
                  + pa_x[k] * ph_9[k];

        t_20[k] = pa_y[k] * ph_4[k];

        t_21[k] = f_7 * pg_6[k]
                  + pb_x[k] * dg_10[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, t_27, pa_x, pa_z, pb_z, pg_0, ph_0, \
                         ph_10, ph_11, ph_12, ph_13, dg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pa_x[k] * ph_10[k];

        t_23[k] = pa_x[k] * ph_11[k];

        t_24[k] = pa_x[k] * ph_12[k];

        t_25[k] = pa_x[k] * ph_13[k];

        t_26[k] = pa_z[k] * ph_0[k];

        t_27[k] = f_7 * pg_0[k]
                  + pb_z[k] * dg_11[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_x, pa_z, pb_x, pg_10, pg_11, pg_15, \
                         ph_1, ph_3, ph_17, ph_18, dg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_z[k] * ph_1[k];

        t_29[k] = f_8 * pg_10[k]
                  + pa_x[k] * ph_17[k];

        t_30[k] = pa_z[k] * ph_3[k];

        t_31[k] = f_0 * pg_11[k]
                  + pa_x[k] * ph_18[k];

        t_32[k] = f_7 * pg_15[k]
                  + pb_x[k] * dg_12[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_x, pb_x, ph_20, ph_21, ph_22, ph_23, \
                         df0_6, df1_6, dg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_x[k] * ph_20[k];

        t_34[k] = pa_x[k] * ph_21[k];

        t_35[k] = pa_x[k] * ph_22[k];

        t_36[k] = pa_x[k] * ph_23[k];

        t_37[k] = f_1 * df0_6[k]
                  - f_2 * df1_6[k]
                  + pb_x[k] * dg_13[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pb_x, pb_y, pg_3, df0_7, df0_8, df1_7, df1_8, \
                         dg_13, dg_14, dg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * pg_3[k]
                  + pb_y[k] * dg_13[k];

        t_39[k] = f_5 * df0_7[k]
                  - f_6 * df1_7[k]
                  + pb_x[k] * dg_14[k];

        t_40[k] = f_5 * df0_8[k]
                  - f_6 * df1_8[k]
                  + pb_x[k] * dg_15[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pb_x, df0_9, df0_11, df1_9, df1_11, \
                         dg_16, dg_17, dg_18, dg_20, dg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_3 * df0_9[k]
                  - f_4 * df1_9[k]
                  + pb_x[k] * dg_16[k];

        t_42[k] = f_3 * df0_11[k]
                  - f_4 * df1_11[k]
                  + pb_x[k] * dg_17[k];

        t_43[k] = pb_x[k] * dg_18[k];

        t_44[k] = pb_x[k] * dg_20[k];

        t_45[k] = pb_x[k] * dg_21[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_y, pb_z, pg_6, df0_9, df0_10, df1_9, \
                         df1_10, dg_18, dg_19, dg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_0 * pg_6[k]
                  + f_1 * df0_9[k]
                  - f_2 * df1_9[k]
                  + pb_y[k] * dg_18[k];

        t_47[k] = pb_z[k] * dg_18[k];

        t_48[k] = f_3 * df0_9[k]
                  - f_4 * df1_9[k]
                  + pb_z[k] * dg_19[k];

        t_49[k] = f_5 * df0_10[k]
                  - f_6 * df1_10[k]
                  + pb_z[k] * dg_20[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_y, pa_z, pb_y, pb_z, pg_8, ph_8, \
                         ph_16, ph_17, df0_11, df1_11, dg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_0 * pg_8[k]
                  + pb_y[k] * dg_21[k];

        t_51[k] = f_1 * df0_11[k]
                  - f_2 * df1_11[k]
                  + pb_z[k] * dg_21[k];

        t_52[k] = pa_y[k] * ph_16[k];

        t_53[k] = pa_z[k] * ph_8[k];

        t_54[k] = pa_y[k] * ph_17[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_y, pa_z, pb_z, pg_6, pg_13, ph_9, \
                         ph_10, ph_18, ph_21, dg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_z[k] * ph_9[k];

        t_56[k] = pa_y[k] * ph_18[k];

        t_57[k] = pa_z[k] * ph_10[k];

        t_58[k] = f_7 * pg_6[k]
                  + pb_z[k] * dg_22[k];

        t_59[k] = f_8 * pg_13[k]
                  + pa_y[k] * ph_21[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pb_x, pb_y, pg_14, pg_15, ph_22, ph_23, \
                         df0_12, df1_12, dg_23, dg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * pg_14[k]
                  + pa_y[k] * ph_22[k];

        t_61[k] = f_7 * pg_15[k]
                  + pb_y[k] * dg_23[k];

        t_62[k] = pa_y[k] * ph_23[k];

        t_63[k] = f_1 * df0_12[k]
                  - f_2 * df1_12[k]
                  + pb_x[k] * dg_24[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_x, pb_z, pg_9, df0_13, df0_14, df1_13, df1_14, \
                         dg_24, dg_25, dg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * pg_9[k]
                  + pb_z[k] * dg_24[k];

        t_65[k] = f_5 * df0_13[k]
                  - f_6 * df1_13[k]
                  + pb_x[k] * dg_25[k];

        t_66[k] = f_5 * df0_14[k]
                  - f_6 * df1_14[k]
                  + pb_x[k] * dg_26[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, pb_x, df0_15, df0_17, df1_15, df1_17, \
                         dg_27, dg_28, dg_29, dg_30, dg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_3 * df0_15[k]
                  - f_4 * df1_15[k]
                  + pb_x[k] * dg_27[k];

        t_68[k] = f_3 * df0_17[k]
                  - f_4 * df1_17[k]
                  + pb_x[k] * dg_28[k];

        t_69[k] = pb_x[k] * dg_29[k];

        t_70[k] = pb_x[k] * dg_30[k];

        t_71[k] = pb_x[k] * dg_32[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pb_y, pb_z, pg_12, df0_15, df0_16, df0_17, \
                         df1_15, df1_16, df1_17, dg_29, dg_30, dg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * df0_15[k]
                  - f_2 * df1_15[k]
                  + pb_y[k] * dg_29[k];

        t_73[k] = f_0 * pg_12[k]
                  + pb_z[k] * dg_29[k];

        t_74[k] = f_5 * df0_16[k]
                  - f_6 * df1_16[k]
                  + pb_y[k] * dg_30[k];

        t_75[k] = f_3 * df0_17[k]
                  - f_4 * df1_17[k]
                  + pb_y[k] * dg_31[k];
    }

#pragma omp simd aligned(t_76, t_77, pb_y, pb_z, pg_15, df0_17, df1_17, \
                         dg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pb_y[k] * dg_32[k];

        t_77[k] = f_0 * pg_15[k]
                  + f_1 * df0_17[k]
                  - f_2 * df1_17[k]
                  + pb_z[k] * dg_32[k];
    }
}

auto
compute_prim_dh_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pg, const size_t ph,
                                     const size_t df0, const size_t df1, const size_t dg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.5 / p;

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

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_11 = buffer.data(pg + 11);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_14 = buffer.data(pg + 14);
    const auto *pg_15 = buffer.data(pg + 15);
    const auto *pg_16 = buffer.data(pg + 16);
    const auto *pg_18 = buffer.data(pg + 18);
    const auto *pg_19 = buffer.data(pg + 19);
    const auto *pg_20 = buffer.data(pg + 20);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_8 = buffer.data(ph + 8);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);
    const auto *df0_3 = buffer.data(df0 + 3);
    const auto *df0_4 = buffer.data(df0 + 4);
    const auto *df0_5 = buffer.data(df0 + 5);
    const auto *df0_6 = buffer.data(df0 + 6);
    const auto *df0_7 = buffer.data(df0 + 7);
    const auto *df0_8 = buffer.data(df0 + 8);
    const auto *df0_9 = buffer.data(df0 + 9);
    const auto *df0_10 = buffer.data(df0 + 10);
    const auto *df0_11 = buffer.data(df0 + 11);
    const auto *df0_12 = buffer.data(df0 + 12);
    const auto *df0_13 = buffer.data(df0 + 13);
    const auto *df0_14 = buffer.data(df0 + 14);
    const auto *df0_15 = buffer.data(df0 + 15);
    const auto *df0_16 = buffer.data(df0 + 16);
    const auto *df0_17 = buffer.data(df0 + 17);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_1 = buffer.data(df1 + 1);
    const auto *df1_2 = buffer.data(df1 + 2);
    const auto *df1_3 = buffer.data(df1 + 3);
    const auto *df1_5 = buffer.data(df1 + 5);
    const auto *df1_6 = buffer.data(df1 + 6);
    const auto *df1_11 = buffer.data(df1 + 11);
    const auto *df1_13 = buffer.data(df1 + 13);
    const auto *df1_14 = buffer.data(df1 + 14);
    const auto *df1_15 = buffer.data(df1 + 15);
    const auto *df1_16 = buffer.data(df1 + 16);
    const auto *df1_17 = buffer.data(df1 + 17);
    const auto *df1_20 = buffer.data(df1 + 20);
    const auto *df1_22 = buffer.data(df1 + 22);
    const auto *df1_23 = buffer.data(df1 + 23);
    const auto *df1_24 = buffer.data(df1 + 24);
    const auto *df1_25 = buffer.data(df1 + 25);
    const auto *df1_26 = buffer.data(df1 + 26);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_40 = buffer.data(dg + 40);
    const auto *dg_41 = buffer.data(dg + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pg_0, df0_0, df0_1, df1_0, \
                         df1_1, dg_0, dg_1, dg_2, dg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_y[k] * dg_1[k];

        t_2[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_z[k] * dg_2[k];

        t_3[k] = f_5 * df0_1[k]
                 - f_6 * df1_1[k]
                 + pb_y[k] * dg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, pg_3, pg_4, df0_2, df0_3, \
                         df1_2, df1_3, dg_4, dg_5, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * df0_2[k]
                 - f_6 * df1_2[k]
                 + pb_z[k] * dg_4[k];

        t_5[k] = f_0 * pg_3[k]
                 + pb_x[k] * dg_5[k];

        t_6[k] = f_0 * pg_4[k]
                 + pb_x[k] * dg_8[k];

        t_7[k] = f_1 * df0_3[k]
                 - f_2 * df1_3[k]
                 + pb_y[k] * dg_5[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, pg_0, df0_4, df0_5, df1_5, df1_6, \
                         dg_6, dg_7, dg_8, dg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * df0_4[k]
                 - f_6 * df1_5[k]
                 + pb_y[k] * dg_6[k];

        t_9[k] = f_3 * df0_5[k]
                 - f_4 * df1_6[k]
                 + pb_y[k] * dg_7[k];

        t_10[k] = f_1 * df0_5[k]
                  - f_2 * df1_6[k]
                  + pb_z[k] * dg_8[k];

        t_11[k] = f_7 * pg_0[k]
                  + pb_y[k] * dg_9[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pa_z, pb_x, pg_6, pg_7, pg_8, ph_0, \
                         ph_1, ph_2, dg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_8 * pg_6[k]
                  + pa_x[k] * ph_1[k];

        t_13[k] = f_0 * pg_7[k]
                  + pa_x[k] * ph_2[k];

        t_14[k] = f_7 * pg_8[k]
                  + pb_x[k] * dg_11[k];

        t_15[k] = pa_z[k] * ph_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_x, pb_z, pg_0, pg_14, pg_15, pg_20, \
                         ph_4, ph_5, dg_12, dg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * pg_0[k]
                  + pb_z[k] * dg_12[k];

        t_17[k] = f_8 * pg_14[k]
                  + pa_x[k] * ph_4[k];

        t_18[k] = f_0 * pg_15[k]
                  + pa_x[k] * ph_5[k];

        t_19[k] = f_7 * pg_20[k]
                  + pb_x[k] * dg_15[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pb_x, pb_y, pg_5, df0_6, df0_7, df0_8, \
                         df1_11, df1_13, df1_14, dg_16, dg_17, dg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * df0_6[k]
                  - f_2 * df1_11[k]
                  + pb_x[k] * dg_16[k];

        t_21[k] = f_0 * pg_5[k]
                  + pb_y[k] * dg_16[k];

        t_22[k] = f_5 * df0_7[k]
                  - f_6 * df1_13[k]
                  + pb_x[k] * dg_17[k];

        t_23[k] = f_5 * df0_8[k]
                  - f_6 * df1_14[k]
                  + pb_x[k] * dg_18[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_x, pb_y, pb_z, pg_8, df0_9, df0_11, \
                         df1_15, df1_17, dg_19, dg_20, dg_21, dg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_3 * df0_9[k]
                  - f_4 * df1_15[k]
                  + pb_x[k] * dg_19[k];

        t_25[k] = f_3 * df0_11[k]
                  - f_4 * df1_17[k]
                  + pb_x[k] * dg_20[k];

        t_26[k] = f_0 * pg_8[k]
                  + f_1 * df0_9[k]
                  - f_2 * df1_15[k]
                  + pb_y[k] * dg_21[k];

        t_27[k] = f_3 * df0_9[k]
                  - f_4 * df1_15[k]
                  + pb_z[k] * dg_22[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_z, pb_y, pb_z, pg_11, ph_3, df0_10, \
                         df0_11, df1_16, df1_17, dg_23, dg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_5 * df0_10[k]
                  - f_6 * df1_16[k]
                  + pb_z[k] * dg_23[k];

        t_29[k] = f_0 * pg_11[k]
                  + pb_y[k] * dg_25[k];

        t_30[k] = f_1 * df0_11[k]
                  - f_2 * df1_17[k]
                  + pb_z[k] * dg_25[k];

        t_31[k] = pa_z[k] * ph_3[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_y, pb_y, pb_z, pg_8, pg_18, pg_19, pg_20, \
                         ph_6, ph_7, dg_26, dg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_7 * pg_8[k]
                  + pb_z[k] * dg_26[k];

        t_33[k] = f_8 * pg_18[k]
                  + pa_y[k] * ph_6[k];

        t_34[k] = f_0 * pg_19[k]
                  + pa_y[k] * ph_7[k];

        t_35[k] = f_7 * pg_20[k]
                  + pb_y[k] * dg_30[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_y, pb_x, pb_z, pg_12, ph_8, df0_12, \
                         df0_13, df1_20, df1_22, dg_31, dg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_y[k] * ph_8[k];

        t_37[k] = f_1 * df0_12[k]
                  - f_2 * df1_20[k]
                  + pb_x[k] * dg_31[k];

        t_38[k] = f_0 * pg_12[k]
                  + pb_z[k] * dg_31[k];

        t_39[k] = f_5 * df0_13[k]
                  - f_6 * df1_22[k]
                  + pb_x[k] * dg_33[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_x, pb_y, df0_14, df0_15, df0_17, df1_23, \
                         df1_24, df1_26, dg_34, dg_35, dg_36, dg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_5 * df0_14[k]
                  - f_6 * df1_23[k]
                  + pb_x[k] * dg_34[k];

        t_41[k] = f_3 * df0_15[k]
                  - f_4 * df1_24[k]
                  + pb_x[k] * dg_35[k];

        t_42[k] = f_3 * df0_17[k]
                  - f_4 * df1_26[k]
                  + pb_x[k] * dg_36[k];

        t_43[k] = f_1 * df0_15[k]
                  - f_2 * df1_24[k]
                  + pb_y[k] * dg_37[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_y, pb_z, pg_16, pg_20, df0_16, df0_17, \
                         df1_25, df1_26, dg_37, dg_39, dg_40, dg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * pg_16[k]
                  + pb_z[k] * dg_37[k];

        t_45[k] = f_5 * df0_16[k]
                  - f_6 * df1_25[k]
                  + pb_y[k] * dg_39[k];

        t_46[k] = f_3 * df0_17[k]
                  - f_4 * df1_26[k]
                  + pb_y[k] * dg_40[k];

        t_47[k] = f_0 * pg_20[k]
                  + f_1 * df0_17[k]
                  - f_2 * df1_26[k]
                  + pb_z[k] * dg_41[k];
    }
}

auto
compute_prim_dh_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pg, const size_t ph,
                                     const size_t df0, const size_t df1, const size_t dg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / p;
    const auto f_8 = 0.5 / p;

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

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_8 = buffer.data(ph + 8);
    const auto *ph_9 = buffer.data(ph + 9);
    const auto *ph_11 = buffer.data(ph + 11);
    const auto *ph_12 = buffer.data(ph + 12);
    const auto *ph_13 = buffer.data(ph + 13);
    const auto *ph_14 = buffer.data(ph + 14);
    const auto *ph_17 = buffer.data(ph + 17);
    const auto *ph_18 = buffer.data(ph + 18);
    const auto *ph_20 = buffer.data(ph + 20);
    const auto *ph_23 = buffer.data(ph + 23);
    const auto *ph_24 = buffer.data(ph + 24);
    const auto *ph_25 = buffer.data(ph + 25);
    const auto *ph_27 = buffer.data(ph + 27);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);
    const auto *df0_3 = buffer.data(df0 + 3);
    const auto *df0_4 = buffer.data(df0 + 4);
    const auto *df0_5 = buffer.data(df0 + 5);
    const auto *df0_6 = buffer.data(df0 + 6);
    const auto *df0_7 = buffer.data(df0 + 7);
    const auto *df0_8 = buffer.data(df0 + 8);
    const auto *df0_9 = buffer.data(df0 + 9);
    const auto *df0_10 = buffer.data(df0 + 10);
    const auto *df0_11 = buffer.data(df0 + 11);
    const auto *df0_12 = buffer.data(df0 + 12);
    const auto *df0_13 = buffer.data(df0 + 13);
    const auto *df0_14 = buffer.data(df0 + 14);
    const auto *df0_15 = buffer.data(df0 + 15);
    const auto *df0_16 = buffer.data(df0 + 16);
    const auto *df0_17 = buffer.data(df0 + 17);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_1 = buffer.data(df1 + 1);
    const auto *df1_2 = buffer.data(df1 + 2);
    const auto *df1_3 = buffer.data(df1 + 3);
    const auto *df1_4 = buffer.data(df1 + 4);
    const auto *df1_5 = buffer.data(df1 + 5);
    const auto *df1_6 = buffer.data(df1 + 6);
    const auto *df1_7 = buffer.data(df1 + 7);
    const auto *df1_8 = buffer.data(df1 + 8);
    const auto *df1_9 = buffer.data(df1 + 9);
    const auto *df1_10 = buffer.data(df1 + 10);
    const auto *df1_11 = buffer.data(df1 + 11);
    const auto *df1_12 = buffer.data(df1 + 12);
    const auto *df1_13 = buffer.data(df1 + 13);
    const auto *df1_14 = buffer.data(df1 + 14);
    const auto *df1_15 = buffer.data(df1 + 15);
    const auto *df1_16 = buffer.data(df1 + 16);
    const auto *df1_17 = buffer.data(df1 + 17);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pg_0, df0_0, df1_0, dg_0, \
                         dg_1, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = pb_y[k] * dg_0[k];

        t_2[k] = pb_z[k] * dg_0[k];

        t_3[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_y[k] * dg_1[k];

        t_4[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_z[k] * dg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, df0_1, df0_2, df0_3, df1_1, \
                         df1_2, df1_3, dg_3, dg_4, dg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * df0_1[k]
                 - f_6 * df1_1[k]
                 + pb_y[k] * dg_3[k];

        t_6[k] = pb_z[k] * dg_3[k];

        t_7[k] = pb_y[k] * dg_4[k];

        t_8[k] = f_5 * df0_2[k]
                 - f_6 * df1_2[k]
                 + pb_z[k] * dg_4[k];

        t_9[k] = f_1 * df0_3[k]
                 - f_2 * df1_3[k]
                 + pb_y[k] * dg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pb_y, pb_z, df0_4, df0_5, df1_4, df1_5, \
                         dg_5, dg_7, dg_8, dg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_z[k] * dg_5[k];

        t_11[k] = f_5 * df0_4[k]
                  - f_6 * df1_4[k]
                  + pb_y[k] * dg_7[k];

        t_12[k] = f_3 * df0_5[k]
                  - f_4 * df1_5[k]
                  + pb_y[k] * dg_8[k];

        t_13[k] = pb_y[k] * dg_9[k];

        t_14[k] = f_1 * df0_5[k]
                  - f_2 * df1_5[k]
                  + pb_z[k] * dg_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, t_20, pa_x, pa_y, pg_1, pg_2, ph_0, \
                         ph_2, ph_4, ph_8, ph_9, ph_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * ph_0[k];

        t_16[k] = f_7 * pg_1[k]
                  + pa_x[k] * ph_8[k];

        t_17[k] = pa_y[k] * ph_2[k];

        t_18[k] = f_0 * pg_2[k]
                  + pa_x[k] * ph_9[k];

        t_19[k] = pa_y[k] * ph_4[k];

        t_20[k] = pa_x[k] * ph_11[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, t_26, pa_x, pa_z, pb_z, pg_0, ph_0, \
                         ph_1, ph_12, ph_13, ph_14, dg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_x[k] * ph_12[k];

        t_22[k] = pa_x[k] * ph_13[k];

        t_23[k] = pa_x[k] * ph_14[k];

        t_24[k] = pa_z[k] * ph_0[k];

        t_25[k] = f_8 * pg_0[k]
                  + pb_z[k] * dg_12[k];

        t_26[k] = pa_z[k] * ph_1[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, t_32, pa_x, pa_z, pg_6, pg_7, ph_3, \
                         ph_18, ph_20, ph_23, ph_24, ph_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * pg_6[k]
                  + pa_x[k] * ph_18[k];

        t_28[k] = pa_z[k] * ph_3[k];

        t_29[k] = f_0 * pg_7[k]
                  + pa_x[k] * ph_20[k];

        t_30[k] = pa_x[k] * ph_23[k];

        t_31[k] = pa_x[k] * ph_24[k];

        t_32[k] = pa_x[k] * ph_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_x, pb_x, pb_z, ph_27, df0_6, df0_7, df1_6, \
                         df1_7, dg_14, dg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_x[k] * ph_27[k];

        t_34[k] = f_1 * df0_6[k]
                  - f_2 * df1_6[k]
                  + pb_x[k] * dg_14[k];

        t_35[k] = pb_z[k] * dg_14[k];

        t_36[k] = f_5 * df0_7[k]
                  - f_6 * df1_7[k]
                  + pb_x[k] * dg_16[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_x, pb_z, df0_8, df0_9, df0_11, df1_8, \
                         df1_9, df1_11, dg_16, dg_17, dg_18, dg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * df0_8[k]
                  - f_6 * df1_8[k]
                  + pb_x[k] * dg_17[k];

        t_38[k] = f_3 * df0_9[k]
                  - f_4 * df1_9[k]
                  + pb_x[k] * dg_18[k];

        t_39[k] = pb_z[k] * dg_16[k];

        t_40[k] = f_3 * df0_11[k]
                  - f_4 * df1_11[k]
                  + pb_x[k] * dg_19[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, t_46, pb_x, pb_y, pb_z, pg_3, df0_9, \
                         df1_9, dg_20, dg_21, dg_22, dg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pb_x[k] * dg_20[k];

        t_42[k] = pb_x[k] * dg_22[k];

        t_43[k] = pb_x[k] * dg_23[k];

        t_44[k] = f_0 * pg_3[k]
                  + f_1 * df0_9[k]
                  - f_2 * df1_9[k]
                  + pb_y[k] * dg_20[k];

        t_45[k] = pb_z[k] * dg_20[k];

        t_46[k] = f_3 * df0_9[k]
                  - f_4 * df1_9[k]
                  + pb_z[k] * dg_21[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pb_y, pb_z, pg_4, ph_17, df0_10, \
                         df0_11, df1_10, df1_11, dg_22, dg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_5 * df0_10[k]
                  - f_6 * df1_10[k]
                  + pb_z[k] * dg_22[k];

        t_48[k] = f_0 * pg_4[k]
                  + pb_y[k] * dg_23[k];

        t_49[k] = f_1 * df0_11[k]
                  - f_2 * df1_11[k]
                  + pb_z[k] * dg_23[k];

        t_50[k] = pa_y[k] * ph_17[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, pa_y, pa_z, pb_x, ph_8, ph_9, \
                         ph_11, ph_18, ph_20, dg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_z[k] * ph_8[k];

        t_52[k] = pa_y[k] * ph_18[k];

        t_53[k] = pa_z[k] * ph_9[k];

        t_54[k] = pa_y[k] * ph_20[k];

        t_55[k] = pb_x[k] * dg_25[k];

        t_56[k] = pa_z[k] * ph_11[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_y, pb_y, pb_z, pg_3, pg_9, pg_10, pg_11, \
                         ph_24, ph_25, dg_24, dg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_8 * pg_3[k]
                  + pb_z[k] * dg_24[k];

        t_58[k] = f_7 * pg_9[k]
                  + pa_y[k] * ph_24[k];

        t_59[k] = f_0 * pg_10[k]
                  + pa_y[k] * ph_25[k];

        t_60[k] = f_8 * pg_11[k]
                  + pb_y[k] * dg_26[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_y, pb_x, pb_y, pb_z, pg_5, ph_27, df0_12, \
                         df1_12, dg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pa_y[k] * ph_27[k];

        t_62[k] = f_1 * df0_12[k]
                  - f_2 * df1_12[k]
                  + pb_x[k] * dg_27[k];

        t_63[k] = pb_y[k] * dg_27[k];

        t_64[k] = f_0 * pg_5[k]
                  + pb_z[k] * dg_27[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_x, pb_y, df0_13, df0_14, df0_15, df1_13, \
                         df1_14, df1_15, dg_29, dg_30, dg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_5 * df0_13[k]
                  - f_6 * df1_13[k]
                  + pb_x[k] * dg_29[k];

        t_66[k] = f_5 * df0_14[k]
                  - f_6 * df1_14[k]
                  + pb_x[k] * dg_30[k];

        t_67[k] = f_3 * df0_15[k]
                  - f_4 * df1_15[k]
                  + pb_x[k] * dg_31[k];

        t_68[k] = pb_y[k] * dg_30[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pb_x, pb_y, df0_15, df0_17, df1_15, \
                         df1_17, dg_32, dg_33, dg_34, dg_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * df0_17[k]
                  - f_4 * df1_17[k]
                  + pb_x[k] * dg_32[k];

        t_70[k] = pb_x[k] * dg_33[k];

        t_71[k] = pb_x[k] * dg_34[k];

        t_72[k] = pb_x[k] * dg_36[k];

        t_73[k] = f_1 * df0_15[k]
                  - f_2 * df1_15[k]
                  + pb_y[k] * dg_33[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_y, pb_z, pg_8, df0_16, df0_17, df1_16, \
                         df1_17, dg_33, dg_34, dg_35, dg_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * pg_8[k]
                  + pb_z[k] * dg_33[k];

        t_75[k] = f_5 * df0_16[k]
                  - f_6 * df1_16[k]
                  + pb_y[k] * dg_34[k];

        t_76[k] = f_3 * df0_17[k]
                  - f_4 * df1_17[k]
                  + pb_y[k] * dg_35[k];

        t_77[k] = pb_y[k] * dg_36[k];
    }

#pragma omp simd aligned(t_78, pb_z, pg_11, df0_17, df1_17, dg_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * pg_11[k]
                  + f_1 * df0_17[k]
                  - f_2 * df1_17[k]
                  + pb_z[k] * dg_36[k];
    }
}

auto
compute_prim_dh_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pg, const size_t ph,
                                     const size_t df0, const size_t df1, const size_t dg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / p;
    const auto f_8 = 0.5 / p;

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

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_8 = buffer.data(ph + 8);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);
    const auto *df0_3 = buffer.data(df0 + 3);
    const auto *df0_4 = buffer.data(df0 + 4);
    const auto *df0_5 = buffer.data(df0 + 5);
    const auto *df0_6 = buffer.data(df0 + 6);
    const auto *df0_7 = buffer.data(df0 + 7);
    const auto *df0_8 = buffer.data(df0 + 8);
    const auto *df0_9 = buffer.data(df0 + 9);
    const auto *df0_10 = buffer.data(df0 + 10);
    const auto *df0_11 = buffer.data(df0 + 11);
    const auto *df0_12 = buffer.data(df0 + 12);
    const auto *df0_13 = buffer.data(df0 + 13);
    const auto *df0_14 = buffer.data(df0 + 14);
    const auto *df0_15 = buffer.data(df0 + 15);
    const auto *df0_16 = buffer.data(df0 + 16);
    const auto *df0_17 = buffer.data(df0 + 17);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_1 = buffer.data(df1 + 1);
    const auto *df1_2 = buffer.data(df1 + 2);
    const auto *df1_3 = buffer.data(df1 + 3);
    const auto *df1_4 = buffer.data(df1 + 4);
    const auto *df1_5 = buffer.data(df1 + 5);
    const auto *df1_6 = buffer.data(df1 + 6);
    const auto *df1_7 = buffer.data(df1 + 7);
    const auto *df1_8 = buffer.data(df1 + 8);
    const auto *df1_9 = buffer.data(df1 + 9);
    const auto *df1_10 = buffer.data(df1 + 10);
    const auto *df1_11 = buffer.data(df1 + 11);
    const auto *df1_12 = buffer.data(df1 + 12);
    const auto *df1_13 = buffer.data(df1 + 13);
    const auto *df1_14 = buffer.data(df1 + 14);
    const auto *df1_15 = buffer.data(df1 + 15);
    const auto *df1_16 = buffer.data(df1 + 16);
    const auto *df1_17 = buffer.data(df1 + 17);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pg_0, df0_0, df1_0, dg_0, \
                         dg_1, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = pb_y[k] * dg_0[k];

        t_2[k] = pb_z[k] * dg_0[k];

        t_3[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_y[k] * dg_1[k];

        t_4[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_z[k] * dg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, df0_1, df0_2, df0_3, df1_1, df1_2, \
                         df1_3, dg_3, dg_4, dg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * df0_1[k]
                 - f_6 * df1_1[k]
                 + pb_y[k] * dg_3[k];

        t_6[k] = pb_y[k] * dg_4[k];

        t_7[k] = f_5 * df0_2[k]
                 - f_6 * df1_2[k]
                 + pb_z[k] * dg_4[k];

        t_8[k] = f_1 * df0_3[k]
                 - f_2 * df1_3[k]
                 + pb_y[k] * dg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, ph_0, df0_4, df0_5, \
                         df1_4, df1_5, dg_6, dg_7, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * df0_4[k]
                 - f_6 * df1_4[k]
                 + pb_y[k] * dg_6[k];

        t_10[k] = f_3 * df0_5[k]
                  - f_4 * df1_5[k]
                  + pb_y[k] * dg_7[k];

        t_11[k] = pb_y[k] * dg_8[k];

        t_12[k] = f_1 * df0_5[k]
                  - f_2 * df1_5[k]
                  + pb_z[k] * dg_8[k];

        t_13[k] = pa_y[k] * ph_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_x, pa_z, pb_z, pg_0, pg_1, pg_2, \
                         ph_0, ph_1, ph_2, ph_3, dg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_7 * pg_1[k]
                  + pa_x[k] * ph_1[k];

        t_15[k] = f_0 * pg_2[k]
                  + pa_x[k] * ph_2[k];

        t_16[k] = pa_x[k] * ph_3[k];

        t_17[k] = pa_z[k] * ph_0[k];

        t_18[k] = f_8 * pg_0[k]
                  + pb_z[k] * dg_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, pg_6, pg_7, ph_4, ph_5, ph_8, \
                         df0_6, df1_6, dg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * pg_6[k]
                  + pa_x[k] * ph_4[k];

        t_20[k] = f_0 * pg_7[k]
                  + pa_x[k] * ph_5[k];

        t_21[k] = pa_x[k] * ph_8[k];

        t_22[k] = f_1 * df0_6[k]
                  - f_2 * df1_6[k]
                  + pb_x[k] * dg_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_x, df0_7, df0_8, df0_9, df1_7, df1_8, df1_9, \
                         dg_11, dg_12, dg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * df0_7[k]
                  - f_6 * df1_7[k]
                  + pb_x[k] * dg_11[k];

        t_24[k] = f_5 * df0_8[k]
                  - f_6 * df1_8[k]
                  + pb_x[k] * dg_12[k];

        t_25[k] = f_3 * df0_9[k]
                  - f_4 * df1_9[k]
                  + pb_x[k] * dg_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pb_x, pb_y, pg_3, df0_9, df0_11, df1_9, \
                         df1_11, dg_14, dg_15, dg_17, dg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * df0_11[k]
                  - f_4 * df1_11[k]
                  + pb_x[k] * dg_14[k];

        t_27[k] = pb_x[k] * dg_15[k];

        t_28[k] = pb_x[k] * dg_17[k];

        t_29[k] = pb_x[k] * dg_18[k];

        t_30[k] = f_0 * pg_3[k]
                  + f_1 * df0_9[k]
                  - f_2 * df1_9[k]
                  + pb_y[k] * dg_15[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pb_y, pb_z, pg_4, df0_9, df0_10, df1_9, \
                         df1_10, dg_15, dg_16, dg_17, dg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_z[k] * dg_15[k];

        t_32[k] = f_3 * df0_9[k]
                  - f_4 * df1_9[k]
                  + pb_z[k] * dg_16[k];

        t_33[k] = f_5 * df0_10[k]
                  - f_6 * df1_10[k]
                  + pb_z[k] * dg_17[k];

        t_34[k] = f_0 * pg_4[k]
                  + pb_y[k] * dg_18[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_y, pa_z, pb_z, pg_3, pg_9, ph_3, ph_6, \
                         df0_11, df1_11, dg_18, dg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * df0_11[k]
                  - f_2 * df1_11[k]
                  + pb_z[k] * dg_18[k];

        t_36[k] = pa_z[k] * ph_3[k];

        t_37[k] = f_8 * pg_3[k]
                  + pb_z[k] * dg_19[k];

        t_38[k] = f_7 * pg_9[k]
                  + pa_y[k] * ph_6[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_y, pb_x, pb_y, pg_10, pg_11, ph_7, ph_8, \
                         df0_12, df1_12, dg_20, dg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * pg_10[k]
                  + pa_y[k] * ph_7[k];

        t_40[k] = f_8 * pg_11[k]
                  + pb_y[k] * dg_20[k];

        t_41[k] = pa_y[k] * ph_8[k];

        t_42[k] = f_1 * df0_12[k]
                  - f_2 * df1_12[k]
                  + pb_x[k] * dg_21[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pb_x, pb_z, pg_5, df0_13, df0_14, df1_13, df1_14, \
                         dg_21, dg_22, dg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_0 * pg_5[k]
                  + pb_z[k] * dg_21[k];

        t_44[k] = f_5 * df0_13[k]
                  - f_6 * df1_13[k]
                  + pb_x[k] * dg_22[k];

        t_45[k] = f_5 * df0_14[k]
                  - f_6 * df1_14[k]
                  + pb_x[k] * dg_23[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pb_x, df0_15, df0_17, df1_15, df1_17, \
                         dg_24, dg_25, dg_26, dg_27, dg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_3 * df0_15[k]
                  - f_4 * df1_15[k]
                  + pb_x[k] * dg_24[k];

        t_47[k] = f_3 * df0_17[k]
                  - f_4 * df1_17[k]
                  + pb_x[k] * dg_25[k];

        t_48[k] = pb_x[k] * dg_26[k];

        t_49[k] = pb_x[k] * dg_27[k];

        t_50[k] = pb_x[k] * dg_29[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_y, pb_z, pg_8, df0_15, df0_16, df0_17, \
                         df1_15, df1_16, df1_17, dg_26, dg_27, dg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * df0_15[k]
                  - f_2 * df1_15[k]
                  + pb_y[k] * dg_26[k];

        t_52[k] = f_0 * pg_8[k]
                  + pb_z[k] * dg_26[k];

        t_53[k] = f_5 * df0_16[k]
                  - f_6 * df1_16[k]
                  + pb_y[k] * dg_27[k];

        t_54[k] = f_3 * df0_17[k]
                  - f_4 * df1_17[k]
                  + pb_y[k] * dg_28[k];
    }

#pragma omp simd aligned(t_55, t_56, pb_y, pb_z, pg_11, df0_17, df1_17, \
                         dg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pb_y[k] * dg_29[k];

        t_56[k] = f_0 * pg_11[k]
                  + f_1 * df0_17[k]
                  - f_2 * df1_17[k]
                  + pb_z[k] * dg_29[k];
    }
}

auto
compute_prim_dh_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pg, const size_t df0, const size_t df1,
                                     const size_t dg, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_1 = buffer.data(df1 + 1);
    const auto *df1_2 = buffer.data(df1 + 2);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pg_0, pg_1, df0_0, df0_1, df1_0, df1_1, dg_0, \
                         dg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_0 * pg_1[k]
                 + f_1 * df0_1[k]
                 - f_2 * df1_1[k]
                 + pb_y[k] * dg_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pg_2, df0_2, df1_2, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pg_2[k]
                 + f_1 * df0_2[k]
                 - f_2 * df1_2[k]
                 + pb_z[k] * dg_2[k];
    }
}

auto
compute_prim_dh_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pg, const size_t df0, const size_t df1,
                                     const size_t dg, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_12 = buffer.data(df1 + 12);
    const auto *df1_23 = buffer.data(df1 + 23);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_41 = buffer.data(dg + 41);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pg_0, pg_1, df0_0, df0_1, df1_0, df1_12, dg_0, \
                         dg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_0 * pg_1[k]
                 + f_1 * df0_1[k]
                 - f_2 * df1_12[k]
                 + pb_y[k] * dg_21[k];
    }

#pragma omp simd aligned(t_2, pb_z, pg_2, df0_2, df1_23, dg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pg_2[k]
                 + f_1 * df0_2[k]
                 - f_2 * df1_23[k]
                 + pb_z[k] * dg_41[k];
    }
}

auto
compute_prim_dh_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pg, const size_t df0, const size_t df1,
                                     const size_t dg, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);

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

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_5 = buffer.data(pg + 5);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);
    const auto *df0_3 = buffer.data(df0 + 3);
    const auto *df0_5 = buffer.data(df0 + 5);
    const auto *df0_6 = buffer.data(df0 + 6);
    const auto *df0_8 = buffer.data(df0 + 8);
    const auto *df0_10 = buffer.data(df0 + 10);
    const auto *df0_11 = buffer.data(df0 + 11);
    const auto *df0_12 = buffer.data(df0 + 12);
    const auto *df0_13 = buffer.data(df0 + 13);
    const auto *df0_14 = buffer.data(df0 + 14);
    const auto *df0_17 = buffer.data(df0 + 17);
    const auto *df0_19 = buffer.data(df0 + 19);
    const auto *df0_20 = buffer.data(df0 + 20);
    const auto *df0_21 = buffer.data(df0 + 21);
    const auto *df0_22 = buffer.data(df0 + 22);
    const auto *df0_23 = buffer.data(df0 + 23);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_1 = buffer.data(df1 + 1);
    const auto *df1_2 = buffer.data(df1 + 2);
    const auto *df1_3 = buffer.data(df1 + 3);
    const auto *df1_5 = buffer.data(df1 + 5);
    const auto *df1_6 = buffer.data(df1 + 6);
    const auto *df1_8 = buffer.data(df1 + 8);
    const auto *df1_10 = buffer.data(df1 + 10);
    const auto *df1_11 = buffer.data(df1 + 11);
    const auto *df1_12 = buffer.data(df1 + 12);
    const auto *df1_13 = buffer.data(df1 + 13);
    const auto *df1_14 = buffer.data(df1 + 14);
    const auto *df1_17 = buffer.data(df1 + 17);
    const auto *df1_19 = buffer.data(df1 + 19);
    const auto *df1_20 = buffer.data(df1 + 20);
    const auto *df1_21 = buffer.data(df1 + 21);
    const auto *df1_22 = buffer.data(df1 + 22);
    const auto *df1_23 = buffer.data(df1 + 23);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pg_0, df0_0, df0_1, df1_0, \
                         df1_1, dg_0, dg_1, dg_2, dg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_y[k] * dg_1[k];

        t_2[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_z[k] * dg_2[k];

        t_3[k] = f_5 * df0_1[k]
                 - f_6 * df1_1[k]
                 + pb_y[k] * dg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, df0_2, df0_3, df0_5, df1_2, df1_3, df1_5, \
                         dg_4, dg_5, dg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * df0_2[k]
                 - f_6 * df1_2[k]
                 + pb_z[k] * dg_4[k];

        t_5[k] = f_1 * df0_3[k]
                 - f_2 * df1_3[k]
                 + pb_y[k] * dg_5[k];

        t_6[k] = f_5 * df0_5[k]
                 - f_6 * df1_5[k]
                 + pb_y[k] * dg_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_x, pb_y, pb_z, df0_6, df0_8, df1_6, df1_8, dg_7, \
                         dg_8, dg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * df0_6[k]
                 - f_4 * df1_6[k]
                 + pb_y[k] * dg_7[k];

        t_8[k] = f_1 * df0_6[k]
                 - f_2 * df1_6[k]
                 + pb_z[k] * dg_8[k];

        t_9[k] = f_1 * df0_8[k]
                 - f_2 * df1_8[k]
                 + pb_x[k] * dg_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, df0_10, df0_11, df0_12, df1_10, df1_11, \
                         df1_12, dg_10, dg_11, dg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * df0_10[k]
                  - f_6 * df1_10[k]
                  + pb_x[k] * dg_10[k];

        t_11[k] = f_5 * df0_11[k]
                  - f_6 * df1_11[k]
                  + pb_x[k] * dg_11[k];

        t_12[k] = f_3 * df0_12[k]
                  - f_4 * df1_12[k]
                  + pb_x[k] * dg_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, pb_y, pb_z, pg_2, df0_12, df0_14, df1_12, \
                         df1_14, dg_13, dg_14, dg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * df0_14[k]
                  - f_4 * df1_14[k]
                  + pb_x[k] * dg_13[k];

        t_14[k] = f_0 * pg_2[k]
                  + f_1 * df0_12[k]
                  - f_2 * df1_12[k]
                  + pb_y[k] * dg_14[k];

        t_15[k] = f_3 * df0_12[k]
                  - f_4 * df1_12[k]
                  + pb_z[k] * dg_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, pb_z, df0_13, df0_14, df0_17, df1_13, df1_14, \
                         df1_17, dg_16, dg_17, dg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * df0_13[k]
                  - f_6 * df1_13[k]
                  + pb_z[k] * dg_16[k];

        t_17[k] = f_1 * df0_14[k]
                  - f_2 * df1_14[k]
                  + pb_z[k] * dg_17[k];

        t_18[k] = f_1 * df0_17[k]
                  - f_2 * df1_17[k]
                  + pb_x[k] * dg_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_x, df0_19, df0_20, df0_21, df1_19, df1_20, \
                         df1_21, dg_19, dg_20, dg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_5 * df0_19[k]
                  - f_6 * df1_19[k]
                  + pb_x[k] * dg_19[k];

        t_20[k] = f_5 * df0_20[k]
                  - f_6 * df1_20[k]
                  + pb_x[k] * dg_20[k];

        t_21[k] = f_3 * df0_21[k]
                  - f_4 * df1_21[k]
                  + pb_x[k] * dg_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pb_x, pb_y, df0_21, df0_22, df0_23, df1_21, \
                         df1_22, df1_23, dg_22, dg_23, dg_24, dg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * df0_23[k]
                  - f_4 * df1_23[k]
                  + pb_x[k] * dg_22[k];

        t_23[k] = f_1 * df0_21[k]
                  - f_2 * df1_21[k]
                  + pb_y[k] * dg_23[k];

        t_24[k] = f_5 * df0_22[k]
                  - f_6 * df1_22[k]
                  + pb_y[k] * dg_24[k];

        t_25[k] = f_3 * df0_23[k]
                  - f_4 * df1_23[k]
                  + pb_y[k] * dg_25[k];
    }

#pragma omp simd aligned(t_26, pb_z, pg_5, df0_23, df1_23, dg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * pg_5[k]
                  + f_1 * df0_23[k]
                  - f_2 * df1_23[k]
                  + pb_z[k] * dg_26[k];
    }
}

auto
compute_prim_dh_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pg, const size_t df0, const size_t df1,
                                     const size_t dg, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_9 = buffer.data(df1 + 9);
    const auto *df1_17 = buffer.data(df1 + 17);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_32 = buffer.data(dg + 32);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pg_0, pg_1, df0_0, df0_1, df1_0, df1_9, dg_0, \
                         dg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_0 * pg_1[k]
                 + f_1 * df0_1[k]
                 - f_2 * df1_9[k]
                 + pb_y[k] * dg_17[k];
    }

#pragma omp simd aligned(t_2, pb_z, pg_2, df0_2, df1_17, dg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pg_2[k]
                 + f_1 * df0_2[k]
                 - f_2 * df1_17[k]
                 + pb_z[k] * dg_32[k];
    }
}

auto
compute_prim_dh_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pg, const size_t ph,
                                     const size_t df0, const size_t df1, const size_t dg,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / p;
    const auto f_8 = 0.5 / p;

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

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_8 = buffer.data(ph + 8);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);
    const auto *df0_3 = buffer.data(df0 + 3);
    const auto *df0_4 = buffer.data(df0 + 4);
    const auto *df0_5 = buffer.data(df0 + 5);
    const auto *df0_6 = buffer.data(df0 + 6);
    const auto *df0_7 = buffer.data(df0 + 7);
    const auto *df0_8 = buffer.data(df0 + 8);
    const auto *df0_9 = buffer.data(df0 + 9);
    const auto *df0_10 = buffer.data(df0 + 10);
    const auto *df0_11 = buffer.data(df0 + 11);
    const auto *df0_12 = buffer.data(df0 + 12);
    const auto *df0_13 = buffer.data(df0 + 13);
    const auto *df0_14 = buffer.data(df0 + 14);
    const auto *df0_15 = buffer.data(df0 + 15);
    const auto *df0_16 = buffer.data(df0 + 16);
    const auto *df0_17 = buffer.data(df0 + 17);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_1 = buffer.data(df1 + 1);
    const auto *df1_2 = buffer.data(df1 + 2);
    const auto *df1_3 = buffer.data(df1 + 3);
    const auto *df1_5 = buffer.data(df1 + 5);
    const auto *df1_6 = buffer.data(df1 + 6);
    const auto *df1_7 = buffer.data(df1 + 7);
    const auto *df1_9 = buffer.data(df1 + 9);
    const auto *df1_10 = buffer.data(df1 + 10);
    const auto *df1_11 = buffer.data(df1 + 11);
    const auto *df1_12 = buffer.data(df1 + 12);
    const auto *df1_13 = buffer.data(df1 + 13);
    const auto *df1_14 = buffer.data(df1 + 14);
    const auto *df1_16 = buffer.data(df1 + 16);
    const auto *df1_17 = buffer.data(df1 + 17);
    const auto *df1_18 = buffer.data(df1 + 18);
    const auto *df1_19 = buffer.data(df1 + 19);
    const auto *df1_20 = buffer.data(df1 + 20);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pg_0, df0_0, df1_0, dg_0, \
                         dg_1, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = pb_y[k] * dg_0[k];

        t_2[k] = pb_z[k] * dg_0[k];

        t_3[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_y[k] * dg_1[k];

        t_4[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_z[k] * dg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, df0_1, df0_2, df0_3, df1_1, \
                         df1_2, df1_3, dg_3, dg_4, dg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * df0_1[k]
                 - f_6 * df1_1[k]
                 + pb_y[k] * dg_3[k];

        t_6[k] = pb_z[k] * dg_3[k];

        t_7[k] = f_5 * df0_2[k]
                 - f_6 * df1_2[k]
                 + pb_z[k] * dg_4[k];

        t_8[k] = f_1 * df0_3[k]
                 - f_2 * df1_3[k]
                 + pb_y[k] * dg_5[k];

        t_9[k] = pb_z[k] * dg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, ph_0, df0_4, df0_5, df1_5, \
                         df1_6, dg_7, dg_8, dg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * df0_4[k]
                  - f_6 * df1_5[k]
                  + pb_y[k] * dg_7[k];

        t_11[k] = f_3 * df0_5[k]
                  - f_4 * df1_6[k]
                  + pb_y[k] * dg_8[k];

        t_12[k] = f_1 * df0_5[k]
                  - f_2 * df1_6[k]
                  + pb_z[k] * dg_9[k];

        t_13[k] = pa_y[k] * ph_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_x, pa_z, pb_z, pg_0, pg_1, pg_2, \
                         ph_0, ph_1, ph_2, ph_3, dg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_7 * pg_1[k]
                  + pa_x[k] * ph_1[k];

        t_15[k] = f_0 * pg_2[k]
                  + pa_x[k] * ph_2[k];

        t_16[k] = pa_x[k] * ph_3[k];

        t_17[k] = pa_z[k] * ph_0[k];

        t_18[k] = f_8 * pg_0[k]
                  + pb_z[k] * dg_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, pg_6, pg_7, ph_4, \
                         ph_5, ph_8, df0_6, df1_7, dg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * pg_6[k]
                  + pa_x[k] * ph_4[k];

        t_20[k] = f_0 * pg_7[k]
                  + pa_x[k] * ph_5[k];

        t_21[k] = pa_x[k] * ph_8[k];

        t_22[k] = f_1 * df0_6[k]
                  - f_2 * df1_7[k]
                  + pb_x[k] * dg_11[k];

        t_23[k] = pb_z[k] * dg_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_x, pb_z, df0_7, df0_8, df0_9, df1_9, \
                         df1_10, df1_11, dg_13, dg_14, dg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * df0_7[k]
                  - f_6 * df1_9[k]
                  + pb_x[k] * dg_13[k];

        t_25[k] = f_5 * df0_8[k]
                  - f_6 * df1_10[k]
                  + pb_x[k] * dg_14[k];

        t_26[k] = f_3 * df0_9[k]
                  - f_4 * df1_11[k]
                  + pb_x[k] * dg_15[k];

        t_27[k] = pb_z[k] * dg_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pb_x, pb_y, pb_z, pg_3, df0_9, df0_11, \
                         df1_11, df1_13, dg_16, dg_17, dg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * df0_11[k]
                  - f_4 * df1_13[k]
                  + pb_x[k] * dg_16[k];

        t_29[k] = pb_x[k] * dg_17[k];

        t_30[k] = f_0 * pg_3[k]
                  + f_1 * df0_9[k]
                  - f_2 * df1_11[k]
                  + pb_y[k] * dg_17[k];

        t_31[k] = pb_z[k] * dg_17[k];

        t_32[k] = f_3 * df0_9[k]
                  - f_4 * df1_11[k]
                  + pb_z[k] * dg_18[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, pg_4, ph_3, df0_10, df0_11, \
                         df1_12, df1_13, dg_19, dg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_5 * df0_10[k]
                  - f_6 * df1_12[k]
                  + pb_z[k] * dg_19[k];

        t_34[k] = f_0 * pg_4[k]
                  + pb_y[k] * dg_20[k];

        t_35[k] = f_1 * df0_11[k]
                  - f_2 * df1_13[k]
                  + pb_z[k] * dg_20[k];

        t_36[k] = pa_z[k] * ph_3[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_y, pb_y, pb_z, pg_3, pg_9, pg_10, pg_11, \
                         ph_6, ph_7, dg_21, dg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_8 * pg_3[k]
                  + pb_z[k] * dg_21[k];

        t_38[k] = f_7 * pg_9[k]
                  + pa_y[k] * ph_6[k];

        t_39[k] = f_0 * pg_10[k]
                  + pa_y[k] * ph_7[k];

        t_40[k] = f_8 * pg_11[k]
                  + pb_y[k] * dg_22[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_y, pb_x, pb_y, pb_z, pg_5, ph_8, df0_12, \
                         df1_14, dg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pa_y[k] * ph_8[k];

        t_42[k] = f_1 * df0_12[k]
                  - f_2 * df1_14[k]
                  + pb_x[k] * dg_23[k];

        t_43[k] = pb_y[k] * dg_23[k];

        t_44[k] = f_0 * pg_5[k]
                  + pb_z[k] * dg_23[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_x, pb_y, df0_13, df0_14, df0_15, df1_16, \
                         df1_17, df1_18, dg_25, dg_26, dg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_5 * df0_13[k]
                  - f_6 * df1_16[k]
                  + pb_x[k] * dg_25[k];

        t_46[k] = f_5 * df0_14[k]
                  - f_6 * df1_17[k]
                  + pb_x[k] * dg_26[k];

        t_47[k] = f_3 * df0_15[k]
                  - f_4 * df1_18[k]
                  + pb_x[k] * dg_27[k];

        t_48[k] = pb_y[k] * dg_26[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pb_x, pb_y, pb_z, pg_8, df0_15, df0_17, \
                         df1_18, df1_20, dg_28, dg_29, dg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * df0_17[k]
                  - f_4 * df1_20[k]
                  + pb_x[k] * dg_28[k];

        t_50[k] = pb_x[k] * dg_32[k];

        t_51[k] = f_1 * df0_15[k]
                  - f_2 * df1_18[k]
                  + pb_y[k] * dg_29[k];

        t_52[k] = f_0 * pg_8[k]
                  + pb_z[k] * dg_29[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_y, pb_z, pg_11, df0_16, df0_17, df1_19, \
                         df1_20, dg_30, dg_31, dg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_5 * df0_16[k]
                  - f_6 * df1_19[k]
                  + pb_y[k] * dg_30[k];

        t_54[k] = f_3 * df0_17[k]
                  - f_4 * df1_20[k]
                  + pb_y[k] * dg_31[k];

        t_55[k] = pb_y[k] * dg_32[k];

        t_56[k] = f_0 * pg_11[k]
                  + f_1 * df0_17[k]
                  - f_2 * df1_20[k]
                  + pb_z[k] * dg_32[k];
    }
}

auto
compute_prim_dh_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pg, const size_t df0, const size_t df1,
                                      const size_t dg, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);
    const auto *df0_3 = buffer.data(df0 + 3);
    const auto *df0_5 = buffer.data(df0 + 5);
    const auto *df0_6 = buffer.data(df0 + 6);
    const auto *df0_7 = buffer.data(df0 + 7);
    const auto *df0_9 = buffer.data(df0 + 9);
    const auto *df0_10 = buffer.data(df0 + 10);
    const auto *df0_11 = buffer.data(df0 + 11);
    const auto *df0_12 = buffer.data(df0 + 12);
    const auto *df0_13 = buffer.data(df0 + 13);
    const auto *df0_14 = buffer.data(df0 + 14);
    const auto *df0_16 = buffer.data(df0 + 16);
    const auto *df0_17 = buffer.data(df0 + 17);
    const auto *df0_18 = buffer.data(df0 + 18);
    const auto *df0_19 = buffer.data(df0 + 19);
    const auto *df0_20 = buffer.data(df0 + 20);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_1 = buffer.data(df1 + 1);
    const auto *df1_2 = buffer.data(df1 + 2);
    const auto *df1_3 = buffer.data(df1 + 3);
    const auto *df1_4 = buffer.data(df1 + 4);
    const auto *df1_5 = buffer.data(df1 + 5);
    const auto *df1_6 = buffer.data(df1 + 6);
    const auto *df1_7 = buffer.data(df1 + 7);
    const auto *df1_8 = buffer.data(df1 + 8);
    const auto *df1_9 = buffer.data(df1 + 9);
    const auto *df1_10 = buffer.data(df1 + 10);
    const auto *df1_11 = buffer.data(df1 + 11);
    const auto *df1_12 = buffer.data(df1 + 12);
    const auto *df1_13 = buffer.data(df1 + 13);
    const auto *df1_14 = buffer.data(df1 + 14);
    const auto *df1_15 = buffer.data(df1 + 15);
    const auto *df1_16 = buffer.data(df1 + 16);
    const auto *df1_17 = buffer.data(df1 + 17);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pg_0, df0_0, df1_0, dg_0, \
                         dg_1, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = pb_y[k] * dg_0[k];

        t_2[k] = pb_z[k] * dg_0[k];

        t_3[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_y[k] * dg_1[k];

        t_4[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_z[k] * dg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, df0_1, df0_2, df0_3, df1_1, df1_2, \
                         df1_3, dg_3, dg_4, dg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * df0_1[k]
                 - f_6 * df1_1[k]
                 + pb_y[k] * dg_3[k];

        t_6[k] = pb_y[k] * dg_4[k];

        t_7[k] = f_5 * df0_2[k]
                 - f_6 * df1_2[k]
                 + pb_z[k] * dg_4[k];

        t_8[k] = f_1 * df0_3[k]
                 - f_2 * df1_3[k]
                 + pb_y[k] * dg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, df0_5, df0_6, df1_4, df1_5, dg_6, \
                         dg_7, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * df0_5[k]
                 - f_6 * df1_4[k]
                 + pb_y[k] * dg_6[k];

        t_10[k] = f_3 * df0_6[k]
                  - f_4 * df1_5[k]
                  + pb_y[k] * dg_7[k];

        t_11[k] = pb_y[k] * dg_8[k];

        t_12[k] = f_1 * df0_6[k]
                  - f_2 * df1_5[k]
                  + pb_z[k] * dg_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, df0_7, df0_9, df0_10, df1_6, df1_7, df1_8, \
                         dg_9, dg_10, dg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * df0_7[k]
                  - f_2 * df1_6[k]
                  + pb_x[k] * dg_9[k];

        t_14[k] = f_5 * df0_9[k]
                  - f_6 * df1_7[k]
                  + pb_x[k] * dg_10[k];

        t_15[k] = f_5 * df0_10[k]
                  - f_6 * df1_8[k]
                  + pb_x[k] * dg_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pb_x, df0_11, df0_13, df1_9, df1_11, \
                         dg_12, dg_13, dg_14, dg_16, dg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * df0_11[k]
                  - f_4 * df1_9[k]
                  + pb_x[k] * dg_12[k];

        t_17[k] = f_3 * df0_13[k]
                  - f_4 * df1_11[k]
                  + pb_x[k] * dg_13[k];

        t_18[k] = pb_x[k] * dg_14[k];

        t_19[k] = pb_x[k] * dg_16[k];

        t_20[k] = pb_x[k] * dg_17[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, pb_z, pg_1, df0_11, df0_12, df1_9, \
                         df1_10, dg_14, dg_15, dg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * pg_1[k]
                  + f_1 * df0_11[k]
                  - f_2 * df1_9[k]
                  + pb_y[k] * dg_14[k];

        t_22[k] = pb_z[k] * dg_14[k];

        t_23[k] = f_3 * df0_11[k]
                  - f_4 * df1_9[k]
                  + pb_z[k] * dg_15[k];

        t_24[k] = f_5 * df0_12[k]
                  - f_6 * df1_10[k]
                  + pb_z[k] * dg_16[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pb_x, pb_z, df0_13, df0_14, df0_16, df1_11, df1_12, \
                         df1_13, dg_17, dg_18, dg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * df0_13[k]
                  - f_2 * df1_11[k]
                  + pb_z[k] * dg_17[k];

        t_26[k] = f_1 * df0_14[k]
                  - f_2 * df1_12[k]
                  + pb_x[k] * dg_18[k];

        t_27[k] = f_5 * df0_16[k]
                  - f_6 * df1_13[k]
                  + pb_x[k] * dg_19[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_x, df0_17, df0_18, df0_20, df1_14, df1_15, \
                         df1_17, dg_20, dg_21, dg_22, dg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_5 * df0_17[k]
                  - f_6 * df1_14[k]
                  + pb_x[k] * dg_20[k];

        t_29[k] = f_3 * df0_18[k]
                  - f_4 * df1_15[k]
                  + pb_x[k] * dg_21[k];

        t_30[k] = f_3 * df0_20[k]
                  - f_4 * df1_17[k]
                  + pb_x[k] * dg_22[k];

        t_31[k] = pb_x[k] * dg_23[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_x, pb_y, df0_18, df0_19, df1_15, df1_16, \
                         dg_23, dg_24, dg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_x[k] * dg_24[k];

        t_33[k] = pb_x[k] * dg_26[k];

        t_34[k] = f_1 * df0_18[k]
                  - f_2 * df1_15[k]
                  + pb_y[k] * dg_23[k];

        t_35[k] = f_5 * df0_19[k]
                  - f_6 * df1_16[k]
                  + pb_y[k] * dg_24[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_y, pb_z, pg_2, df0_20, df1_17, dg_25, \
                         dg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_3 * df0_20[k]
                  - f_4 * df1_17[k]
                  + pb_y[k] * dg_25[k];

        t_37[k] = pb_y[k] * dg_26[k];

        t_38[k] = f_0 * pg_2[k]
                  + f_1 * df0_20[k]
                  - f_2 * df1_17[k]
                  + pb_z[k] * dg_26[k];
    }
}

auto
compute_prim_dh_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pg, const size_t df0, const size_t df1,
                                      const size_t dg, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_10 = buffer.data(df1 + 10);
    const auto *df1_20 = buffer.data(df1 + 20);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_17 = buffer.data(dg + 17);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pg_0, pg_1, df0_0, df0_1, df1_0, df1_10, dg_0, \
                         dg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_0 * pg_1[k]
                 + f_1 * df0_1[k]
                 - f_2 * df1_10[k]
                 + pb_y[k] * dg_9[k];
    }

#pragma omp simd aligned(t_2, pb_z, pg_2, df0_2, df1_20, dg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pg_2[k]
                 + f_1 * df0_2[k]
                 - f_2 * df1_20[k]
                 + pb_z[k] * dg_17[k];
    }
}

auto
compute_prim_dh_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pg, const size_t df0, const size_t df1,
                                      const size_t dg, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_10 = buffer.data(df0 + 10);
    const auto *df0_20 = buffer.data(df0 + 20);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_9 = buffer.data(df1 + 9);
    const auto *df1_17 = buffer.data(df1 + 17);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pg_0, pg_1, df0_0, df0_10, df1_0, df1_9, dg_0, \
                         dg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_0 * pg_1[k]
                 + f_1 * df0_10[k]
                 - f_2 * df1_9[k]
                 + pb_y[k] * dg_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pg_2, df0_20, df1_17, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pg_2[k]
                 + f_1 * df0_20[k]
                 - f_2 * df1_17[k]
                 + pb_z[k] * dg_2[k];
    }
}

auto
compute_prim_dh_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pg, const size_t df0, const size_t df1,
                                      const size_t dg, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_11 = buffer.data(df1 + 11);
    const auto *df1_20 = buffer.data(df1 + 20);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_29 = buffer.data(dg + 29);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pg_0, pg_1, df0_0, df0_1, df1_0, df1_11, dg_0, \
                         dg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_0 * pg_1[k]
                 + f_1 * df0_1[k]
                 - f_2 * df1_11[k]
                 + pb_y[k] * dg_15[k];
    }

#pragma omp simd aligned(t_2, pb_z, pg_2, df0_2, df1_20, dg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pg_2[k]
                 + f_1 * df0_2[k]
                 - f_2 * df1_20[k]
                 + pb_z[k] * dg_29[k];
    }
}

auto
compute_prim_dh_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pg, const size_t df0, const size_t df1,
                                      const size_t dg, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);

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

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);
    const auto *df0_3 = buffer.data(df0 + 3);
    const auto *df0_5 = buffer.data(df0 + 5);
    const auto *df0_6 = buffer.data(df0 + 6);
    const auto *df0_7 = buffer.data(df0 + 7);
    const auto *df0_9 = buffer.data(df0 + 9);
    const auto *df0_10 = buffer.data(df0 + 10);
    const auto *df0_11 = buffer.data(df0 + 11);
    const auto *df0_12 = buffer.data(df0 + 12);
    const auto *df0_13 = buffer.data(df0 + 13);
    const auto *df0_14 = buffer.data(df0 + 14);
    const auto *df0_16 = buffer.data(df0 + 16);
    const auto *df0_17 = buffer.data(df0 + 17);
    const auto *df0_18 = buffer.data(df0 + 18);
    const auto *df0_19 = buffer.data(df0 + 19);
    const auto *df0_20 = buffer.data(df0 + 20);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_1 = buffer.data(df1 + 1);
    const auto *df1_2 = buffer.data(df1 + 2);
    const auto *df1_3 = buffer.data(df1 + 3);
    const auto *df1_4 = buffer.data(df1 + 4);
    const auto *df1_5 = buffer.data(df1 + 5);
    const auto *df1_6 = buffer.data(df1 + 6);
    const auto *df1_7 = buffer.data(df1 + 7);
    const auto *df1_8 = buffer.data(df1 + 8);
    const auto *df1_9 = buffer.data(df1 + 9);
    const auto *df1_10 = buffer.data(df1 + 10);
    const auto *df1_11 = buffer.data(df1 + 11);
    const auto *df1_12 = buffer.data(df1 + 12);
    const auto *df1_13 = buffer.data(df1 + 13);
    const auto *df1_14 = buffer.data(df1 + 14);
    const auto *df1_15 = buffer.data(df1 + 15);
    const auto *df1_16 = buffer.data(df1 + 16);
    const auto *df1_17 = buffer.data(df1 + 17);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pg_0, df0_0, df0_1, df1_0, \
                         df1_1, dg_0, dg_1, dg_2, dg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_y[k] * dg_1[k];

        t_2[k] = f_3 * df0_0[k]
                 - f_4 * df1_0[k]
                 + pb_z[k] * dg_2[k];

        t_3[k] = f_5 * df0_1[k]
                 - f_6 * df1_1[k]
                 + pb_y[k] * dg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, df0_2, df0_3, df0_5, df1_2, df1_3, df1_4, \
                         dg_4, dg_5, dg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * df0_2[k]
                 - f_6 * df1_2[k]
                 + pb_z[k] * dg_4[k];

        t_5[k] = f_1 * df0_3[k]
                 - f_2 * df1_3[k]
                 + pb_y[k] * dg_5[k];

        t_6[k] = f_5 * df0_5[k]
                 - f_6 * df1_4[k]
                 + pb_y[k] * dg_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_x, pb_y, pb_z, df0_6, df0_7, df1_5, df1_6, dg_7, \
                         dg_8, dg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * df0_6[k]
                 - f_4 * df1_5[k]
                 + pb_y[k] * dg_7[k];

        t_8[k] = f_1 * df0_6[k]
                 - f_2 * df1_5[k]
                 + pb_z[k] * dg_8[k];

        t_9[k] = f_1 * df0_7[k]
                 - f_2 * df1_6[k]
                 + pb_x[k] * dg_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, df0_9, df0_10, df0_11, df1_7, df1_8, df1_9, \
                         dg_10, dg_11, dg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * df0_9[k]
                  - f_6 * df1_7[k]
                  + pb_x[k] * dg_10[k];

        t_11[k] = f_5 * df0_10[k]
                  - f_6 * df1_8[k]
                  + pb_x[k] * dg_11[k];

        t_12[k] = f_3 * df0_11[k]
                  - f_4 * df1_9[k]
                  + pb_x[k] * dg_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, pb_y, pb_z, pg_1, df0_11, df0_13, df1_9, \
                         df1_11, dg_13, dg_14, dg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * df0_13[k]
                  - f_4 * df1_11[k]
                  + pb_x[k] * dg_13[k];

        t_14[k] = f_0 * pg_1[k]
                  + f_1 * df0_11[k]
                  - f_2 * df1_9[k]
                  + pb_y[k] * dg_14[k];

        t_15[k] = f_3 * df0_11[k]
                  - f_4 * df1_9[k]
                  + pb_z[k] * dg_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, pb_z, df0_12, df0_13, df0_14, df1_10, df1_11, \
                         df1_12, dg_16, dg_17, dg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * df0_12[k]
                  - f_6 * df1_10[k]
                  + pb_z[k] * dg_16[k];

        t_17[k] = f_1 * df0_13[k]
                  - f_2 * df1_11[k]
                  + pb_z[k] * dg_17[k];

        t_18[k] = f_1 * df0_14[k]
                  - f_2 * df1_12[k]
                  + pb_x[k] * dg_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_x, df0_16, df0_17, df0_18, df1_13, df1_14, \
                         df1_15, dg_19, dg_20, dg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_5 * df0_16[k]
                  - f_6 * df1_13[k]
                  + pb_x[k] * dg_19[k];

        t_20[k] = f_5 * df0_17[k]
                  - f_6 * df1_14[k]
                  + pb_x[k] * dg_20[k];

        t_21[k] = f_3 * df0_18[k]
                  - f_4 * df1_15[k]
                  + pb_x[k] * dg_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pb_x, pb_y, df0_18, df0_19, df0_20, df1_15, \
                         df1_16, df1_17, dg_22, dg_23, dg_24, dg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * df0_20[k]
                  - f_4 * df1_17[k]
                  + pb_x[k] * dg_22[k];

        t_23[k] = f_1 * df0_18[k]
                  - f_2 * df1_15[k]
                  + pb_y[k] * dg_23[k];

        t_24[k] = f_5 * df0_19[k]
                  - f_6 * df1_16[k]
                  + pb_y[k] * dg_24[k];

        t_25[k] = f_3 * df0_20[k]
                  - f_4 * df1_17[k]
                  + pb_y[k] * dg_25[k];
    }

#pragma omp simd aligned(t_26, pb_z, pg_2, df0_20, df1_17, dg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * pg_2[k]
                  + f_1 * df0_20[k]
                  - f_2 * df1_17[k]
                  + pb_z[k] * dg_26[k];
    }
}

auto
compute_prim_dh_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pg, const size_t df0, const size_t df1,
                                      const size_t dg, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_9 = buffer.data(df0 + 9);
    const auto *df0_17 = buffer.data(df0 + 17);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_9 = buffer.data(df1 + 9);
    const auto *df1_17 = buffer.data(df1 + 17);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pg_0, pg_1, df0_0, df0_9, df1_0, df1_9, dg_0, \
                         dg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_0 * pg_1[k]
                 + f_1 * df0_9[k]
                 - f_2 * df1_9[k]
                 + pb_y[k] * dg_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pg_2, df0_17, df1_17, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pg_2[k]
                 + f_1 * df0_17[k]
                 - f_2 * df1_17[k]
                 + pb_z[k] * dg_2[k];
    }
}

auto
compute_prim_dh_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pg, const size_t df0, const size_t df1,
                                      const size_t dg, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_2 = buffer.data(df1 + 2);
    const auto *df1_5 = buffer.data(df1 + 5);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_5 = buffer.data(dg + 5);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pg_0, pg_1, df0_0, df0_1, df1_0, df1_2, dg_0, \
                         dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_0 * pg_1[k]
                 + f_1 * df0_1[k]
                 - f_2 * df1_2[k]
                 + pb_y[k] * dg_2[k];
    }

#pragma omp simd aligned(t_2, pb_z, pg_2, df0_2, df1_5, dg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pg_2[k]
                 + f_1 * df0_2[k]
                 - f_2 * df1_5[k]
                 + pb_z[k] * dg_5[k];
    }
}

auto
compute_prim_dh_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pg, const size_t df0, const size_t df1,
                                      const size_t dg, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_2 = buffer.data(df0 + 2);
    const auto *df0_5 = buffer.data(df0 + 5);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_2 = buffer.data(df1 + 2);
    const auto *df1_5 = buffer.data(df1 + 5);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_5 = buffer.data(dg + 5);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pg_0, pg_1, df0_0, df0_2, df1_0, df1_2, dg_0, \
                         dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_0 * pg_1[k]
                 + f_1 * df0_2[k]
                 - f_2 * df1_2[k]
                 + pb_y[k] * dg_2[k];
    }

#pragma omp simd aligned(t_2, pb_z, pg_2, df0_5, df1_5, dg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pg_2[k]
                 + f_1 * df0_5[k]
                 - f_2 * df1_5[k]
                 + pb_z[k] * dg_5[k];
    }
}

auto
compute_prim_dh_electron_repulsion_18(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pg, const size_t df0, const size_t df1,
                                      const size_t dg, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_2 = buffer.data(df0 + 2);
    const auto *df0_5 = buffer.data(df0 + 5);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_5 = buffer.data(df1 + 5);
    const auto *df1_11 = buffer.data(df1 + 11);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_5 = buffer.data(dg + 5);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pg_0, pg_1, df0_0, df0_2, df1_0, df1_5, dg_0, \
                         dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_0 * pg_1[k]
                 + f_1 * df0_2[k]
                 - f_2 * df1_5[k]
                 + pb_y[k] * dg_2[k];
    }

#pragma omp simd aligned(t_2, pb_z, pg_2, df0_5, df1_11, dg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pg_2[k]
                 + f_1 * df0_5[k]
                 - f_2 * df1_11[k]
                 + pb_z[k] * dg_5[k];
    }
}

auto
compute_prim_dh_electron_repulsion_19(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pg, const size_t df0, const size_t df1,
                                      const size_t dg, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_5 = buffer.data(df0 + 5);
    const auto *df0_11 = buffer.data(df0 + 11);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_2 = buffer.data(df1 + 2);
    const auto *df1_5 = buffer.data(df1 + 5);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pg_0, pg_1, df0_0, df0_5, df1_0, df1_2, dg_0, \
                         dg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_0 * pg_1[k]
                 + f_1 * df0_5[k]
                 - f_2 * df1_2[k]
                 + pb_y[k] * dg_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pg_2, df0_11, df1_5, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pg_2[k]
                 + f_1 * df0_11[k]
                 - f_2 * df1_5[k]
                 + pb_z[k] * dg_2[k];
    }
}

auto
compute_prim_dh_electron_repulsion_20(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pg, const size_t df0, const size_t df1,
                                      const size_t dg, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_2 = buffer.data(df1 + 2);
    const auto *df1_5 = buffer.data(df1 + 5);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_8 = buffer.data(dg + 8);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pg_0, pg_1, df0_0, df0_1, df1_0, df1_2, dg_0, \
                         dg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_0 * pg_1[k]
                 + f_1 * df0_1[k]
                 - f_2 * df1_2[k]
                 + pb_y[k] * dg_3[k];
    }

#pragma omp simd aligned(t_2, pb_z, pg_2, df0_2, df1_5, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pg_2[k]
                 + f_1 * df0_2[k]
                 - f_2 * df1_5[k]
                 + pb_z[k] * dg_8[k];
    }
}

auto
compute_prim_dh_electron_repulsion_21(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pg, const size_t df0, const size_t df1,
                                      const size_t dg, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / beta;
    const auto f_4 = alpha / (beta * p);
    const auto f_5 = 0.5 / beta;
    const auto f_6 = 0.5 * alpha / (beta * p);

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

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);
    const auto *df0_3 = buffer.data(df0 + 3);
    const auto *df0_4 = buffer.data(df0 + 4);
    const auto *df0_5 = buffer.data(df0 + 5);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_1 = buffer.data(df1 + 1);
    const auto *df1_2 = buffer.data(df1 + 2);
    const auto *df1_3 = buffer.data(df1 + 3);
    const auto *df1_4 = buffer.data(df1 + 4);
    const auto *df1_5 = buffer.data(df1 + 5);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pg_0, df0_0, df0_1, df0_2, df1_0, df1_1, df1_2, \
                         dg_0, dg_1, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_3 * df0_1[k]
                 - f_4 * df1_1[k]
                 + pb_x[k] * dg_1[k];

        t_2[k] = f_5 * df0_2[k]
                 - f_6 * df1_2[k]
                 + pb_x[k] * dg_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, pg_1, df0_2, df0_3, df0_5, df1_2, df1_3, \
                         df1_5, dg_3, dg_4, dg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_0 * pg_1[k]
                 + f_1 * df0_2[k]
                 - f_2 * df1_2[k]
                 + pb_y[k] * dg_3[k];

        t_4[k] = f_3 * df0_3[k]
                 - f_4 * df1_3[k]
                 + pb_x[k] * dg_4[k];

        t_5[k] = f_5 * df0_5[k]
                 - f_6 * df1_5[k]
                 + pb_x[k] * dg_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_y, pb_z, pg_2, df0_4, df0_5, df1_4, df1_5, dg_6, \
                         dg_7, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * df0_4[k]
                 - f_4 * df1_4[k]
                 + pb_y[k] * dg_6[k];

        t_7[k] = f_5 * df0_5[k]
                 - f_6 * df1_5[k]
                 + pb_y[k] * dg_7[k];

        t_8[k] = f_0 * pg_2[k]
                 + f_1 * df0_5[k]
                 - f_2 * df1_5[k]
                 + pb_z[k] * dg_8[k];
    }
}

auto
compute_prim_dh_electron_repulsion_22(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pg, const size_t df0, const size_t df1,
                                      const size_t dg, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / beta;
    const auto f_4 = alpha / (beta * p);
    const auto f_5 = 0.5 / beta;
    const auto f_6 = 0.5 * alpha / (beta * p);

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

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_1 = buffer.data(df0 + 1);
    const auto *df0_2 = buffer.data(df0 + 2);
    const auto *df0_3 = buffer.data(df0 + 3);
    const auto *df0_4 = buffer.data(df0 + 4);
    const auto *df0_5 = buffer.data(df0 + 5);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_4 = buffer.data(df1 + 4);
    const auto *df1_5 = buffer.data(df1 + 5);
    const auto *df1_8 = buffer.data(df1 + 8);
    const auto *df1_10 = buffer.data(df1 + 10);
    const auto *df1_11 = buffer.data(df1 + 11);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pg_0, df0_0, df0_1, df0_2, df1_0, df1_4, df1_5, \
                         dg_0, dg_4, dg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_3 * df0_1[k]
                 - f_4 * df1_4[k]
                 + pb_x[k] * dg_4[k];

        t_2[k] = f_5 * df0_2[k]
                 - f_6 * df1_5[k]
                 + pb_x[k] * dg_5[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, pg_1, df0_2, df0_3, df0_5, df1_5, df1_8, \
                         df1_11, dg_6, dg_9, dg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_0 * pg_1[k]
                 + f_1 * df0_2[k]
                 - f_2 * df1_5[k]
                 + pb_y[k] * dg_6[k];

        t_4[k] = f_3 * df0_3[k]
                 - f_4 * df1_8[k]
                 + pb_x[k] * dg_9[k];

        t_5[k] = f_5 * df0_5[k]
                 - f_6 * df1_11[k]
                 + pb_x[k] * dg_10[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_y, pb_z, pg_2, df0_4, df0_5, df1_10, df1_11, dg_12, \
                         dg_13, dg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * df0_4[k]
                 - f_4 * df1_10[k]
                 + pb_y[k] * dg_12[k];

        t_7[k] = f_5 * df0_5[k]
                 - f_6 * df1_11[k]
                 + pb_y[k] * dg_13[k];

        t_8[k] = f_0 * pg_2[k]
                 + f_1 * df0_5[k]
                 - f_2 * df1_11[k]
                 + pb_z[k] * dg_14[k];
    }
}

auto
compute_prim_dh_electron_repulsion_23(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pg, const size_t df0, const size_t df1,
                                      const size_t dg, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 1.0 / beta;
    const auto f_4 = alpha / (beta * p);
    const auto f_5 = 0.5 / beta;
    const auto f_6 = 0.5 * alpha / (beta * p);

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

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df0_0 = buffer.data(df0 + 0);
    const auto *df0_4 = buffer.data(df0 + 4);
    const auto *df0_5 = buffer.data(df0 + 5);
    const auto *df0_8 = buffer.data(df0 + 8);
    const auto *df0_10 = buffer.data(df0 + 10);
    const auto *df0_11 = buffer.data(df0 + 11);

    const auto *df1_0 = buffer.data(df1 + 0);
    const auto *df1_4 = buffer.data(df1 + 4);
    const auto *df1_5 = buffer.data(df1 + 5);
    const auto *df1_8 = buffer.data(df1 + 8);
    const auto *df1_10 = buffer.data(df1 + 10);
    const auto *df1_11 = buffer.data(df1 + 11);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pg_0, df0_0, df0_4, df0_5, df1_0, df1_4, df1_5, \
                         dg_0, dg_1, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 + f_1 * df0_0[k]
                 - f_2 * df1_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_3 * df0_4[k]
                 - f_4 * df1_4[k]
                 + pb_x[k] * dg_1[k];

        t_2[k] = f_5 * df0_5[k]
                 - f_6 * df1_5[k]
                 + pb_x[k] * dg_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, pg_1, df0_5, df0_8, df0_11, df1_5, df1_8, \
                         df1_11, dg_3, dg_4, dg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_0 * pg_1[k]
                 + f_1 * df0_5[k]
                 - f_2 * df1_5[k]
                 + pb_y[k] * dg_3[k];

        t_4[k] = f_3 * df0_8[k]
                 - f_4 * df1_8[k]
                 + pb_x[k] * dg_4[k];

        t_5[k] = f_5 * df0_11[k]
                 - f_6 * df1_11[k]
                 + pb_x[k] * dg_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_y, pb_z, pg_2, df0_10, df0_11, df1_10, df1_11, \
                         dg_6, dg_7, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * df0_10[k]
                 - f_4 * df1_10[k]
                 + pb_y[k] * dg_6[k];

        t_7[k] = f_5 * df0_11[k]
                 - f_6 * df1_11[k]
                 + pb_y[k] * dg_7[k];

        t_8[k] = f_0 * pg_2[k]
                 + f_1 * df0_11[k]
                 - f_2 * df1_11[k]
                 + pb_z[k] * dg_8[k];
    }
}

}  // namespace simdt2ceri
