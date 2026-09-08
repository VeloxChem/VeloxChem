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


#include "SimdKineticEnergyVrrRecDG.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_dg_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf, const size_t pg,
                                 const size_t dd_s, const size_t dg_s, const size_t dd,
                                 const size_t df, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;

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

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_7 = buffer.data(dd_s + 7);
    const auto *dd_s_8 = buffer.data(dd_s + 8);
    const auto *dd_s_9 = buffer.data(dd_s + 9);
    const auto *dd_s_10 = buffer.data(dd_s + 10);
    const auto *dd_s_14 = buffer.data(dd_s + 14);
    const auto *dd_s_15 = buffer.data(dd_s + 15);
    const auto *dd_s_16 = buffer.data(dd_s + 16);
    const auto *dd_s_17 = buffer.data(dd_s + 17);
    const auto *dd_s_18 = buffer.data(dd_s + 18);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_1 = buffer.data(dg_s + 1);
    const auto *dg_s_2 = buffer.data(dg_s + 2);
    const auto *dg_s_3 = buffer.data(dg_s + 3);
    const auto *dg_s_4 = buffer.data(dg_s + 4);
    const auto *dg_s_5 = buffer.data(dg_s + 5);
    const auto *dg_s_6 = buffer.data(dg_s + 6);
    const auto *dg_s_7 = buffer.data(dg_s + 7);
    const auto *dg_s_8 = buffer.data(dg_s + 8);
    const auto *dg_s_9 = buffer.data(dg_s + 9);
    const auto *dg_s_10 = buffer.data(dg_s + 10);
    const auto *dg_s_11 = buffer.data(dg_s + 11);
    const auto *dg_s_12 = buffer.data(dg_s + 12);
    const auto *dg_s_13 = buffer.data(dg_s + 13);
    const auto *dg_s_14 = buffer.data(dg_s + 14);
    const auto *dg_s_15 = buffer.data(dg_s + 15);
    const auto *dg_s_16 = buffer.data(dg_s + 16);
    const auto *dg_s_17 = buffer.data(dg_s + 17);
    const auto *dg_s_18 = buffer.data(dg_s + 18);
    const auto *dg_s_19 = buffer.data(dg_s + 19);
    const auto *dg_s_20 = buffer.data(dg_s + 20);
    const auto *dg_s_21 = buffer.data(dg_s + 21);
    const auto *dg_s_22 = buffer.data(dg_s + 22);
    const auto *dg_s_23 = buffer.data(dg_s + 23);
    const auto *dg_s_24 = buffer.data(dg_s + 24);
    const auto *dg_s_25 = buffer.data(dg_s + 25);
    const auto *dg_s_26 = buffer.data(dg_s + 26);
    const auto *dg_s_27 = buffer.data(dg_s + 27);
    const auto *dg_s_28 = buffer.data(dg_s + 28);
    const auto *dg_s_29 = buffer.data(dg_s + 29);
    const auto *dg_s_30 = buffer.data(dg_s + 30);
    const auto *dg_s_31 = buffer.data(dg_s + 31);
    const auto *dg_s_32 = buffer.data(dg_s + 32);
    const auto *dg_s_33 = buffer.data(dg_s + 33);
    const auto *dg_s_34 = buffer.data(dg_s + 34);
    const auto *dg_s_35 = buffer.data(dg_s + 35);
    const auto *dg_s_36 = buffer.data(dg_s + 36);
    const auto *dg_s_37 = buffer.data(dg_s + 37);
    const auto *dg_s_38 = buffer.data(dg_s + 38);
    const auto *dg_s_39 = buffer.data(dg_s + 39);
    const auto *dg_s_40 = buffer.data(dg_s + 40);
    const auto *dg_s_41 = buffer.data(dg_s + 41);
    const auto *dg_s_42 = buffer.data(dg_s + 42);
    const auto *dg_s_43 = buffer.data(dg_s + 43);
    const auto *dg_s_44 = buffer.data(dg_s + 44);
    const auto *dg_s_45 = buffer.data(dg_s + 45);
    const auto *dg_s_46 = buffer.data(dg_s + 46);
    const auto *dg_s_47 = buffer.data(dg_s + 47);
    const auto *dg_s_48 = buffer.data(dg_s + 48);
    const auto *dg_s_49 = buffer.data(dg_s + 49);
    const auto *dg_s_50 = buffer.data(dg_s + 50);
    const auto *dg_s_51 = buffer.data(dg_s + 51);
    const auto *dg_s_52 = buffer.data(dg_s + 52);
    const auto *dg_s_53 = buffer.data(dg_s + 53);
    const auto *dg_s_54 = buffer.data(dg_s + 54);
    const auto *dg_s_55 = buffer.data(dg_s + 55);
    const auto *dg_s_56 = buffer.data(dg_s + 56);
    const auto *dg_s_57 = buffer.data(dg_s + 57);
    const auto *dg_s_58 = buffer.data(dg_s + 58);
    const auto *dg_s_59 = buffer.data(dg_s + 59);
    const auto *dg_s_60 = buffer.data(dg_s + 60);
    const auto *dg_s_61 = buffer.data(dg_s + 61);
    const auto *dg_s_62 = buffer.data(dg_s + 62);
    const auto *dg_s_63 = buffer.data(dg_s + 63);
    const auto *dg_s_64 = buffer.data(dg_s + 64);
    const auto *dg_s_65 = buffer.data(dg_s + 65);
    const auto *dg_s_66 = buffer.data(dg_s + 66);
    const auto *dg_s_67 = buffer.data(dg_s + 67);
    const auto *dg_s_68 = buffer.data(dg_s + 68);
    const auto *dg_s_69 = buffer.data(dg_s + 69);
    const auto *dg_s_70 = buffer.data(dg_s + 70);
    const auto *dg_s_71 = buffer.data(dg_s + 71);
    const auto *dg_s_72 = buffer.data(dg_s + 72);
    const auto *dg_s_73 = buffer.data(dg_s + 73);
    const auto *dg_s_74 = buffer.data(dg_s + 74);
    const auto *dg_s_75 = buffer.data(dg_s + 75);
    const auto *dg_s_76 = buffer.data(dg_s + 76);
    const auto *dg_s_77 = buffer.data(dg_s + 77);
    const auto *dg_s_78 = buffer.data(dg_s + 78);
    const auto *dg_s_79 = buffer.data(dg_s + 79);
    const auto *dg_s_80 = buffer.data(dg_s + 80);
    const auto *dg_s_81 = buffer.data(dg_s + 81);
    const auto *dg_s_82 = buffer.data(dg_s + 82);
    const auto *dg_s_83 = buffer.data(dg_s + 83);
    const auto *dg_s_84 = buffer.data(dg_s + 84);
    const auto *dg_s_85 = buffer.data(dg_s + 85);
    const auto *dg_s_86 = buffer.data(dg_s + 86);
    const auto *dg_s_87 = buffer.data(dg_s + 87);
    const auto *dg_s_88 = buffer.data(dg_s + 88);
    const auto *dg_s_89 = buffer.data(dg_s + 89);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_15 = buffer.data(dd + 15);
    const auto *dd_16 = buffer.data(dd + 16);
    const auto *dd_17 = buffer.data(dd + 17);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pf_0, dd_s_0, dg_s_0, dg_s_1, \
                         dg_s_2, dg_s_3, dd_0, df_0, df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 - f_1 * dd_s_0[k]
                 + f_2 * dg_s_0[k]
                 + f_3 * dd_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_2 * dg_s_1[k]
                 + pb_y[k] * df_0[k];

        t_2[k] = f_2 * dg_s_2[k]
                 + pb_z[k] * df_0[k];

        t_3[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_3[k]
                 + f_5 * dd_0[k]
                 + pb_y[k] * df_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, pb_y, pb_z, pf_1, dd_s_0, dg_s_4, dg_s_5, \
                         dg_s_6, dd_0, df_2, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * dg_s_4[k]
                 + pb_y[k] * df_2[k];

        t_5[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_5[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_2[k];

        t_6[k] = f_0 * pf_1[k]
                 + f_2 * dg_s_6[k]
                 + pb_x[k] * df_5[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_x, pb_y, pb_z, pf_2, dg_s_7, dg_s_8, dg_s_9, df_3, \
                         df_4, df_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * dg_s_7[k]
                 + pb_z[k] * df_3[k];

        t_8[k] = f_2 * dg_s_8[k]
                 + pb_y[k] * df_4[k];

        t_9[k] = f_0 * pf_2[k]
                 + f_2 * dg_s_9[k]
                 + pb_x[k] * df_7[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, pb_z, dd_s_1, dd_s_2, dg_s_10, dg_s_11, \
                         dg_s_12, dd_1, dd_2, df_5, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * dd_s_1[k]
                  + f_2 * dg_s_10[k]
                  + f_3 * dd_1[k]
                  + pb_y[k] * df_5[k];

        t_11[k] = f_2 * dg_s_11[k]
                  + pb_z[k] * df_5[k];

        t_12[k] = -f_4 * dd_s_2[k]
                  + f_2 * dg_s_12[k]
                  + f_5 * dd_2[k]
                  + pb_y[k] * df_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_y, pb_z, pg_0, dd_s_2, dg_s_13, dg_s_14, \
                         dg_s_15, dd_2, df_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * dg_s_13[k]
                  + pb_y[k] * df_7[k];

        t_14[k] = -f_1 * dd_s_2[k]
                  + f_2 * dg_s_14[k]
                  + f_3 * dd_2[k]
                  + pb_z[k] * df_7[k];

        t_15[k] = pa_y[k] * pg_0[k]
                  + f_2 * dg_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, pb_z, pf_0, pf_3, pg_6, dg_s_16, \
                         dg_s_17, dg_s_18, dg_s_19, df_8, df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * pf_0[k]
                  + f_2 * dg_s_16[k]
                  + pb_y[k] * df_8[k];

        t_17[k] = f_2 * dg_s_17[k]
                  + pb_z[k] * df_8[k];

        t_18[k] = f_0 * pf_3[k]
                  + pa_x[k] * pg_6[k]
                  + f_2 * dg_s_18[k];

        t_19[k] = f_2 * dg_s_19[k]
                  + pb_z[k] * df_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_x, pb_z, pf_4, pg_2, dg_s_20, dg_s_21, \
                         dg_s_22, df_10, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_y[k] * pg_2[k]
                  + f_2 * dg_s_20[k];

        t_21[k] = f_5 * pf_4[k]
                  + f_2 * dg_s_21[k]
                  + pb_x[k] * df_11[k];

        t_22[k] = f_2 * dg_s_22[k]
                  + pb_z[k] * df_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pa_y, pb_x, pf_5, pg_4, pg_7, dg_s_23, \
                         dg_s_24, dg_s_25, df_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * pf_5[k]
                  + f_2 * dg_s_23[k]
                  + pb_x[k] * df_12[k];

        t_24[k] = pa_y[k] * pg_4[k]
                  + f_2 * dg_s_24[k];

        t_25[k] = pa_x[k] * pg_7[k]
                  + f_2 * dg_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_x, pb_z, pg_8, pg_9, pg_10, dg_s_26, \
                         dg_s_27, dg_s_28, dg_s_29, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * dg_s_26[k]
                  + pb_z[k] * df_11[k];

        t_27[k] = pa_x[k] * pg_8[k]
                  + f_2 * dg_s_27[k];

        t_28[k] = pa_x[k] * pg_9[k]
                  + f_2 * dg_s_28[k];

        t_29[k] = pa_x[k] * pg_10[k]
                  + f_2 * dg_s_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_z, pb_y, pb_z, pf_0, pg_0, pg_1, dg_s_30, \
                         dg_s_31, dg_s_32, dg_s_33, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_z[k] * pg_0[k]
                  + f_2 * dg_s_30[k];

        t_31[k] = f_2 * dg_s_31[k]
                  + pb_y[k] * df_13[k];

        t_32[k] = f_5 * pf_0[k]
                  + f_2 * dg_s_32[k]
                  + pb_z[k] * df_13[k];

        t_33[k] = pa_z[k] * pg_1[k]
                  + f_2 * dg_s_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_x, pa_z, pb_y, pf_8, pg_3, pg_14, dg_s_34, \
                         dg_s_35, dg_s_36, df_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_2 * dg_s_34[k]
                  + pb_y[k] * df_14[k];

        t_35[k] = f_0 * pf_8[k]
                  + pa_x[k] * pg_14[k]
                  + f_2 * dg_s_35[k];

        t_36[k] = pa_z[k] * pg_3[k]
                  + f_2 * dg_s_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, pb_y, pf_9, pf_11, dg_s_37, dg_s_38, dg_s_39, \
                         df_15, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * pf_9[k]
                  + f_2 * dg_s_37[k]
                  + pb_x[k] * df_16[k];

        t_38[k] = f_2 * dg_s_38[k]
                  + pb_y[k] * df_15[k];

        t_39[k] = f_5 * pf_11[k]
                  + f_2 * dg_s_39[k]
                  + pb_x[k] * df_17[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_y, pg_15, pg_16, pg_17, dg_s_40, \
                         dg_s_41, dg_s_42, dg_s_43, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_x[k] * pg_15[k]
                  + f_2 * dg_s_40[k];

        t_41[k] = pa_x[k] * pg_16[k]
                  + f_2 * dg_s_41[k];

        t_42[k] = pa_x[k] * pg_17[k]
                  + f_2 * dg_s_42[k];

        t_43[k] = f_2 * dg_s_43[k]
                  + pb_y[k] * df_17[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_x, pb_x, pg_18, dd_s_7, dd_s_8, dg_s_44, \
                         dg_s_45, dg_s_46, dd_7, dd_8, df_18, df_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_x[k] * pg_18[k]
                  + f_2 * dg_s_44[k];

        t_45[k] = -f_1 * dd_s_7[k]
                  + f_2 * dg_s_45[k]
                  + f_3 * dd_7[k]
                  + pb_x[k] * df_18[k];

        t_46[k] = -f_6 * dd_s_8[k]
                  + f_2 * dg_s_46[k]
                  + f_0 * dd_8[k]
                  + pb_x[k] * df_19[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pb_x, pb_z, dd_s_9, dg_s_47, dg_s_48, dg_s_49, \
                         dd_9, df_18, df_19, df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_2 * dg_s_47[k]
                  + pb_z[k] * df_18[k];

        t_48[k] = -f_4 * dd_s_9[k]
                  + f_2 * dg_s_48[k]
                  + f_5 * dd_9[k]
                  + pb_x[k] * df_20[k];

        t_49[k] = f_2 * dg_s_49[k]
                  + pb_z[k] * df_19[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pb_x, dd_s_10, dg_s_50, dg_s_51, dg_s_52, \
                         dg_s_53, dd_10, df_21, df_22, df_23, df_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -f_4 * dd_s_10[k]
                  + f_2 * dg_s_50[k]
                  + f_5 * dd_10[k]
                  + pb_x[k] * df_21[k];

        t_51[k] = f_2 * dg_s_51[k]
                  + pb_x[k] * df_22[k];

        t_52[k] = f_2 * dg_s_52[k]
                  + pb_x[k] * df_23[k];

        t_53[k] = f_2 * dg_s_53[k]
                  + pb_x[k] * df_24[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_x, pb_y, pb_z, pf_4, dd_s_9, dg_s_54, dg_s_55, \
                         dg_s_56, dd_9, df_22, df_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_2 * dg_s_54[k]
                  + pb_x[k] * df_25[k];

        t_55[k] = f_0 * pf_4[k]
                  - f_1 * dd_s_9[k]
                  + f_2 * dg_s_55[k]
                  + f_3 * dd_9[k]
                  + pb_y[k] * df_22[k];

        t_56[k] = f_2 * dg_s_56[k]
                  + pb_z[k] * df_22[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_y, pb_z, pf_6, dd_s_9, dd_s_10, dg_s_57, \
                         dg_s_58, dg_s_59, dd_9, dd_10, df_23, df_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = -f_4 * dd_s_9[k]
                  + f_2 * dg_s_57[k]
                  + f_5 * dd_9[k]
                  + pb_z[k] * df_23[k];

        t_58[k] = f_0 * pf_6[k]
                  + f_2 * dg_s_58[k]
                  + pb_y[k] * df_25[k];

        t_59[k] = -f_1 * dd_s_10[k]
                  + f_2 * dg_s_59[k]
                  + f_3 * dd_10[k]
                  + pb_z[k] * df_25[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pa_z, pg_5, pg_6, pg_11, pg_12, \
                         dg_s_60, dg_s_61, dg_s_62, dg_s_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * pg_11[k]
                  + f_2 * dg_s_60[k];

        t_61[k] = pa_z[k] * pg_5[k]
                  + f_2 * dg_s_61[k];

        t_62[k] = pa_y[k] * pg_12[k]
                  + f_2 * dg_s_62[k];

        t_63[k] = pa_z[k] * pg_6[k]
                  + f_2 * dg_s_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pf_7, pg_13, pg_14, dg_s_64, \
                         dg_s_65, dg_s_66, dg_s_67, df_26, df_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_5 * pf_7[k]
                  + pa_y[k] * pg_13[k]
                  + f_2 * dg_s_64[k];

        t_65[k] = pa_y[k] * pg_14[k]
                  + f_2 * dg_s_65[k];

        t_66[k] = f_2 * dg_s_66[k]
                  + pb_x[k] * df_26[k];

        t_67[k] = f_2 * dg_s_67[k]
                  + pb_x[k] * df_27[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_z, pb_x, pb_z, pf_4, pg_7, dg_s_68, \
                         dg_s_69, dg_s_70, dg_s_71, df_26, df_28, \
                         df_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_2 * dg_s_68[k]
                  + pb_x[k] * df_28[k];

        t_69[k] = f_2 * dg_s_69[k]
                  + pb_x[k] * df_29[k];

        t_70[k] = pa_z[k] * pg_7[k]
                  + f_2 * dg_s_70[k];

        t_71[k] = f_5 * pf_4[k]
                  + f_2 * dg_s_71[k]
                  + pb_z[k] * df_26[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_y, pb_y, pf_10, pf_11, pg_17, pg_18, dg_s_72, \
                         dg_s_73, dg_s_74, df_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_0 * pf_10[k]
                  + pa_y[k] * pg_17[k]
                  + f_2 * dg_s_72[k];

        t_73[k] = f_5 * pf_11[k]
                  + f_2 * dg_s_73[k]
                  + pb_y[k] * df_29[k];

        t_74[k] = pa_y[k] * pg_18[k]
                  + f_2 * dg_s_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_x, pb_y, dd_s_14, dd_s_15, dg_s_75, dg_s_76, \
                         dg_s_77, dd_13, dd_14, df_30, df_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -f_1 * dd_s_14[k]
                  + f_2 * dg_s_75[k]
                  + f_3 * dd_13[k]
                  + pb_x[k] * df_30[k];

        t_76[k] = f_2 * dg_s_76[k]
                  + pb_y[k] * df_30[k];

        t_77[k] = -f_6 * dd_s_15[k]
                  + f_2 * dg_s_77[k]
                  + f_0 * dd_14[k]
                  + pb_x[k] * df_31[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pb_x, pb_y, dd_s_16, dd_s_18, dg_s_78, dg_s_79, \
                         dg_s_80, dd_15, dd_17, df_31, df_32, df_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = -f_4 * dd_s_16[k]
                  + f_2 * dg_s_78[k]
                  + f_5 * dd_15[k]
                  + pb_x[k] * df_32[k];

        t_79[k] = f_2 * dg_s_79[k]
                  + pb_y[k] * df_31[k];

        t_80[k] = -f_4 * dd_s_18[k]
                  + f_2 * dg_s_80[k]
                  + f_5 * dd_17[k]
                  + pb_x[k] * df_33[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pb_x, dg_s_81, dg_s_82, dg_s_83, dg_s_84, \
                         df_34, df_35, df_36, df_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_2 * dg_s_81[k]
                  + pb_x[k] * df_34[k];

        t_82[k] = f_2 * dg_s_82[k]
                  + pb_x[k] * df_35[k];

        t_83[k] = f_2 * dg_s_83[k]
                  + pb_x[k] * df_36[k];

        t_84[k] = f_2 * dg_s_84[k]
                  + pb_x[k] * df_37[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pb_y, dd_s_16, dd_s_17, dd_s_18, dg_s_85, dg_s_86, \
                         dg_s_87, dd_15, dd_16, dd_17, df_34, df_35, \
                         df_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -f_1 * dd_s_16[k]
                  + f_2 * dg_s_85[k]
                  + f_3 * dd_15[k]
                  + pb_y[k] * df_34[k];

        t_86[k] = -f_6 * dd_s_17[k]
                  + f_2 * dg_s_86[k]
                  + f_0 * dd_16[k]
                  + pb_y[k] * df_35[k];

        t_87[k] = -f_4 * dd_s_18[k]
                  + f_2 * dg_s_87[k]
                  + f_5 * dd_17[k]
                  + pb_y[k] * df_36[k];
    }

#pragma omp simd aligned(t_88, t_89, pb_y, pb_z, pf_11, dd_s_18, dg_s_88, dg_s_89, dd_17, \
                         df_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_2 * dg_s_88[k]
                  + pb_y[k] * df_37[k];

        t_89[k] = f_0 * pf_11[k]
                  - f_1 * dd_s_18[k]
                  + f_2 * dg_s_89[k]
                  + f_3 * dd_17[k]
                  + pb_z[k] * df_37[k];
    }
}

auto
compute_prim_dg_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf, const size_t pg,
                                 const size_t dd_s, const size_t dg_s, const size_t dd,
                                 const size_t df, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;

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
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_10 = buffer.data(pf + 10);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_16 = buffer.data(pg + 16);
    const auto *pg_17 = buffer.data(pg + 17);
    const auto *pg_19 = buffer.data(pg + 19);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_5 = buffer.data(dd_s + 5);
    const auto *dd_s_6 = buffer.data(dd_s + 6);
    const auto *dd_s_7 = buffer.data(dd_s + 7);
    const auto *dd_s_11 = buffer.data(dd_s + 11);
    const auto *dd_s_12 = buffer.data(dd_s + 12);
    const auto *dd_s_13 = buffer.data(dd_s + 13);
    const auto *dd_s_14 = buffer.data(dd_s + 14);
    const auto *dd_s_15 = buffer.data(dd_s + 15);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_1 = buffer.data(dg_s + 1);
    const auto *dg_s_2 = buffer.data(dg_s + 2);
    const auto *dg_s_3 = buffer.data(dg_s + 3);
    const auto *dg_s_4 = buffer.data(dg_s + 4);
    const auto *dg_s_5 = buffer.data(dg_s + 5);
    const auto *dg_s_6 = buffer.data(dg_s + 6);
    const auto *dg_s_7 = buffer.data(dg_s + 7);
    const auto *dg_s_8 = buffer.data(dg_s + 8);
    const auto *dg_s_9 = buffer.data(dg_s + 9);
    const auto *dg_s_10 = buffer.data(dg_s + 10);
    const auto *dg_s_11 = buffer.data(dg_s + 11);
    const auto *dg_s_12 = buffer.data(dg_s + 12);
    const auto *dg_s_13 = buffer.data(dg_s + 13);
    const auto *dg_s_14 = buffer.data(dg_s + 14);
    const auto *dg_s_15 = buffer.data(dg_s + 15);
    const auto *dg_s_16 = buffer.data(dg_s + 16);
    const auto *dg_s_17 = buffer.data(dg_s + 17);
    const auto *dg_s_18 = buffer.data(dg_s + 18);
    const auto *dg_s_19 = buffer.data(dg_s + 19);
    const auto *dg_s_20 = buffer.data(dg_s + 20);
    const auto *dg_s_22 = buffer.data(dg_s + 22);
    const auto *dg_s_23 = buffer.data(dg_s + 23);
    const auto *dg_s_24 = buffer.data(dg_s + 24);
    const auto *dg_s_25 = buffer.data(dg_s + 25);
    const auto *dg_s_26 = buffer.data(dg_s + 26);
    const auto *dg_s_27 = buffer.data(dg_s + 27);
    const auto *dg_s_28 = buffer.data(dg_s + 28);
    const auto *dg_s_29 = buffer.data(dg_s + 29);
    const auto *dg_s_30 = buffer.data(dg_s + 30);
    const auto *dg_s_31 = buffer.data(dg_s + 31);
    const auto *dg_s_32 = buffer.data(dg_s + 32);
    const auto *dg_s_33 = buffer.data(dg_s + 33);
    const auto *dg_s_34 = buffer.data(dg_s + 34);
    const auto *dg_s_35 = buffer.data(dg_s + 35);
    const auto *dg_s_36 = buffer.data(dg_s + 36);
    const auto *dg_s_37 = buffer.data(dg_s + 37);
    const auto *dg_s_38 = buffer.data(dg_s + 38);
    const auto *dg_s_42 = buffer.data(dg_s + 42);
    const auto *dg_s_43 = buffer.data(dg_s + 43);
    const auto *dg_s_44 = buffer.data(dg_s + 44);
    const auto *dg_s_45 = buffer.data(dg_s + 45);
    const auto *dg_s_46 = buffer.data(dg_s + 46);
    const auto *dg_s_47 = buffer.data(dg_s + 47);
    const auto *dg_s_48 = buffer.data(dg_s + 48);
    const auto *dg_s_49 = buffer.data(dg_s + 49);
    const auto *dg_s_50 = buffer.data(dg_s + 50);
    const auto *dg_s_51 = buffer.data(dg_s + 51);
    const auto *dg_s_52 = buffer.data(dg_s + 52);
    const auto *dg_s_53 = buffer.data(dg_s + 53);
    const auto *dg_s_54 = buffer.data(dg_s + 54);
    const auto *dg_s_55 = buffer.data(dg_s + 55);
    const auto *dg_s_56 = buffer.data(dg_s + 56);
    const auto *dg_s_57 = buffer.data(dg_s + 57);
    const auto *dg_s_58 = buffer.data(dg_s + 58);

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
    const auto *dd_11 = buffer.data(dd + 11);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pf_0, dd_s_0, dg_s_0, dg_s_1, \
                         dg_s_2, dg_s_3, dd_0, df_0, df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 - f_1 * dd_s_0[k]
                 + f_2 * dg_s_0[k]
                 + f_3 * dd_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_2 * dg_s_1[k]
                 + pb_y[k] * df_0[k];

        t_2[k] = f_2 * dg_s_2[k]
                 + pb_z[k] * df_0[k];

        t_3[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_3[k]
                 + f_5 * dd_0[k]
                 + pb_y[k] * df_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, pb_z, pf_1, pf_2, dd_s_0, dg_s_4, dg_s_5, \
                         dg_s_6, dd_0, df_2, df_3, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_4[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_2[k];

        t_5[k] = f_0 * pf_1[k]
                 + f_2 * dg_s_5[k]
                 + pb_x[k] * df_3[k];

        t_6[k] = f_0 * pf_2[k]
                 + f_2 * dg_s_6[k]
                 + pb_x[k] * df_5[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_y, dd_s_1, dd_s_2, dg_s_7, dg_s_8, dg_s_9, dd_1, \
                         dd_2, df_3, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -f_1 * dd_s_1[k]
                 + f_2 * dg_s_7[k]
                 + f_3 * dd_1[k]
                 + pb_y[k] * df_3[k];

        t_8[k] = -f_4 * dd_s_2[k]
                 + f_2 * dg_s_8[k]
                 + f_5 * dd_2[k]
                 + pb_y[k] * df_4[k];

        t_9[k] = f_2 * dg_s_9[k]
                 + pb_y[k] * df_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pb_y, pb_z, pf_0, pg_0, dd_s_2, dg_s_10, \
                         dg_s_11, dg_s_12, dd_2, df_5, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * dd_s_2[k]
                  + f_2 * dg_s_10[k]
                  + f_3 * dd_2[k]
                  + pb_z[k] * df_5[k];

        t_11[k] = pa_y[k] * pg_0[k]
                  + f_2 * dg_s_11[k];

        t_12[k] = f_5 * pf_0[k]
                  + f_2 * dg_s_12[k]
                  + pb_y[k] * df_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pa_y, pb_x, pf_3, pf_4, pg_1, pg_5, dg_s_13, \
                         dg_s_14, dg_s_15, df_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * pf_3[k]
                  + pa_x[k] * pg_5[k]
                  + f_2 * dg_s_13[k];

        t_14[k] = pa_y[k] * pg_1[k]
                  + f_2 * dg_s_14[k];

        t_15[k] = f_5 * pf_4[k]
                  + f_2 * dg_s_15[k]
                  + pb_x[k] * df_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_z, pg_0, pg_7, pg_8, pg_9, dg_s_16, \
                         dg_s_17, dg_s_18, dg_s_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_x[k] * pg_7[k]
                  + f_2 * dg_s_16[k];

        t_17[k] = pa_x[k] * pg_8[k]
                  + f_2 * dg_s_17[k];

        t_18[k] = pa_x[k] * pg_9[k]
                  + f_2 * dg_s_18[k];

        t_19[k] = pa_z[k] * pg_0[k]
                  + f_2 * dg_s_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_x, pb_x, pb_z, pf_0, pf_7, pf_10, pg_13, \
                         dg_s_20, dg_s_22, dg_s_23, df_8, df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * pf_0[k]
                  + f_2 * dg_s_20[k]
                  + pb_z[k] * df_8[k];

        t_21[k] = f_0 * pf_7[k]
                  + pa_x[k] * pg_13[k]
                  + f_2 * dg_s_22[k];

        t_22[k] = f_5 * pf_10[k]
                  + f_2 * dg_s_23[k]
                  + pb_x[k] * df_9[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_x, pg_16, pg_17, pg_19, dd_s_4, \
                         dg_s_24, dg_s_25, dg_s_26, dg_s_27, dd_3, \
                         df_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_x[k] * pg_16[k]
                  + f_2 * dg_s_24[k];

        t_24[k] = pa_x[k] * pg_17[k]
                  + f_2 * dg_s_25[k];

        t_25[k] = pa_x[k] * pg_19[k]
                  + f_2 * dg_s_26[k];

        t_26[k] = -f_1 * dd_s_4[k]
                  + f_2 * dg_s_27[k]
                  + f_3 * dd_3[k]
                  + pb_x[k] * df_10[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_x, dd_s_5, dd_s_6, dd_s_7, dg_s_28, dg_s_29, \
                         dg_s_30, dd_4, dd_5, dd_6, df_11, df_12, \
                         df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_6 * dd_s_5[k]
                  + f_2 * dg_s_28[k]
                  + f_0 * dd_4[k]
                  + pb_x[k] * df_11[k];

        t_28[k] = -f_4 * dd_s_6[k]
                  + f_2 * dg_s_29[k]
                  + f_5 * dd_5[k]
                  + pb_x[k] * df_12[k];

        t_29[k] = -f_4 * dd_s_7[k]
                  + f_2 * dg_s_30[k]
                  + f_5 * dd_6[k]
                  + pb_x[k] * df_13[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_x, pb_y, pf_4, dd_s_6, dg_s_31, dg_s_32, \
                         dg_s_33, dg_s_34, dd_5, df_14, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * dg_s_31[k]
                  + pb_x[k] * df_14[k];

        t_31[k] = f_2 * dg_s_32[k]
                  + pb_x[k] * df_16[k];

        t_32[k] = f_2 * dg_s_33[k]
                  + pb_x[k] * df_17[k];

        t_33[k] = f_0 * pf_4[k]
                  - f_1 * dd_s_6[k]
                  + f_2 * dg_s_34[k]
                  + f_3 * dd_5[k]
                  + pb_y[k] * df_14[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pb_y, pb_z, pf_5, dd_s_6, dg_s_35, dg_s_36, \
                         dg_s_37, dd_5, df_14, df_15, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_2 * dg_s_35[k]
                  + pb_z[k] * df_14[k];

        t_35[k] = -f_4 * dd_s_6[k]
                  + f_2 * dg_s_36[k]
                  + f_5 * dd_5[k]
                  + pb_z[k] * df_15[k];

        t_36[k] = f_0 * pf_5[k]
                  + f_2 * dg_s_37[k]
                  + pb_y[k] * df_17[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_z, pb_z, pf_4, pg_7, dd_s_7, dg_s_38, dg_s_42, \
                         dg_s_43, dd_6, df_17, df_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -f_1 * dd_s_7[k]
                  + f_2 * dg_s_38[k]
                  + f_3 * dd_6[k]
                  + pb_z[k] * df_17[k];

        t_38[k] = pa_z[k] * pg_7[k]
                  + f_2 * dg_s_42[k];

        t_39[k] = f_5 * pf_4[k]
                  + f_2 * dg_s_43[k]
                  + pb_z[k] * df_18[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pb_y, pf_9, pf_10, pg_17, pg_19, dg_s_44, \
                         dg_s_45, dg_s_46, df_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * pf_9[k]
                  + pa_y[k] * pg_17[k]
                  + f_2 * dg_s_44[k];

        t_41[k] = f_5 * pf_10[k]
                  + f_2 * dg_s_45[k]
                  + pb_y[k] * df_19[k];

        t_42[k] = pa_y[k] * pg_19[k]
                  + f_2 * dg_s_46[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pb_x, dd_s_11, dd_s_12, dd_s_13, dg_s_47, dg_s_48, \
                         dg_s_49, dd_7, dd_8, dd_9, df_20, df_21, \
                         df_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -f_1 * dd_s_11[k]
                  + f_2 * dg_s_47[k]
                  + f_3 * dd_7[k]
                  + pb_x[k] * df_20[k];

        t_44[k] = -f_6 * dd_s_12[k]
                  + f_2 * dg_s_48[k]
                  + f_0 * dd_8[k]
                  + pb_x[k] * df_21[k];

        t_45[k] = -f_4 * dd_s_13[k]
                  + f_2 * dg_s_49[k]
                  + f_5 * dd_9[k]
                  + pb_x[k] * df_22[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_x, dd_s_15, dg_s_50, dg_s_51, dg_s_52, \
                         dg_s_53, dd_11, df_23, df_24, df_25, df_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -f_4 * dd_s_15[k]
                  + f_2 * dg_s_50[k]
                  + f_5 * dd_11[k]
                  + pb_x[k] * df_23[k];

        t_47[k] = f_2 * dg_s_51[k]
                  + pb_x[k] * df_24[k];

        t_48[k] = f_2 * dg_s_52[k]
                  + pb_x[k] * df_25[k];

        t_49[k] = f_2 * dg_s_53[k]
                  + pb_x[k] * df_27[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pb_y, dd_s_13, dd_s_14, dd_s_15, dg_s_54, dg_s_55, \
                         dg_s_56, dd_9, dd_10, dd_11, df_24, df_25, \
                         df_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -f_1 * dd_s_13[k]
                  + f_2 * dg_s_54[k]
                  + f_3 * dd_9[k]
                  + pb_y[k] * df_24[k];

        t_51[k] = -f_6 * dd_s_14[k]
                  + f_2 * dg_s_55[k]
                  + f_0 * dd_10[k]
                  + pb_y[k] * df_25[k];

        t_52[k] = -f_4 * dd_s_15[k]
                  + f_2 * dg_s_56[k]
                  + f_5 * dd_11[k]
                  + pb_y[k] * df_26[k];
    }

#pragma omp simd aligned(t_53, t_54, pb_y, pb_z, pf_10, dd_s_15, dg_s_57, dg_s_58, dd_11, \
                         df_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_2 * dg_s_57[k]
                  + pb_y[k] * df_27[k];

        t_54[k] = f_0 * pf_10[k]
                  - f_1 * dd_s_15[k]
                  + f_2 * dg_s_58[k]
                  + f_3 * dd_11[k]
                  + pb_z[k] * df_27[k];
    }
}

auto
compute_prim_dg_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf, const size_t pg,
                                 const size_t dd_s, const size_t dg_s, const size_t dd,
                                 const size_t df, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;

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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_14 = buffer.data(pg + 14);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_5 = buffer.data(dd_s + 5);
    const auto *dd_s_6 = buffer.data(dd_s + 6);
    const auto *dd_s_7 = buffer.data(dd_s + 7);
    const auto *dd_s_11 = buffer.data(dd_s + 11);
    const auto *dd_s_12 = buffer.data(dd_s + 12);
    const auto *dd_s_13 = buffer.data(dd_s + 13);
    const auto *dd_s_14 = buffer.data(dd_s + 14);
    const auto *dd_s_15 = buffer.data(dd_s + 15);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_1 = buffer.data(dg_s + 1);
    const auto *dg_s_2 = buffer.data(dg_s + 2);
    const auto *dg_s_3 = buffer.data(dg_s + 3);
    const auto *dg_s_4 = buffer.data(dg_s + 4);
    const auto *dg_s_5 = buffer.data(dg_s + 5);
    const auto *dg_s_6 = buffer.data(dg_s + 6);
    const auto *dg_s_7 = buffer.data(dg_s + 7);
    const auto *dg_s_8 = buffer.data(dg_s + 8);
    const auto *dg_s_9 = buffer.data(dg_s + 9);
    const auto *dg_s_11 = buffer.data(dg_s + 11);
    const auto *dg_s_12 = buffer.data(dg_s + 12);
    const auto *dg_s_16 = buffer.data(dg_s + 16);
    const auto *dg_s_17 = buffer.data(dg_s + 17);
    const auto *dg_s_18 = buffer.data(dg_s + 18);
    const auto *dg_s_19 = buffer.data(dg_s + 19);
    const auto *dg_s_20 = buffer.data(dg_s + 20);
    const auto *dg_s_21 = buffer.data(dg_s + 21);
    const auto *dg_s_22 = buffer.data(dg_s + 22);
    const auto *dg_s_23 = buffer.data(dg_s + 23);
    const auto *dg_s_24 = buffer.data(dg_s + 24);
    const auto *dg_s_25 = buffer.data(dg_s + 25);
    const auto *dg_s_26 = buffer.data(dg_s + 26);
    const auto *dg_s_27 = buffer.data(dg_s + 27);
    const auto *dg_s_28 = buffer.data(dg_s + 28);
    const auto *dg_s_32 = buffer.data(dg_s + 32);
    const auto *dg_s_36 = buffer.data(dg_s + 36);
    const auto *dg_s_37 = buffer.data(dg_s + 37);
    const auto *dg_s_38 = buffer.data(dg_s + 38);
    const auto *dg_s_39 = buffer.data(dg_s + 39);
    const auto *dg_s_40 = buffer.data(dg_s + 40);
    const auto *dg_s_41 = buffer.data(dg_s + 41);
    const auto *dg_s_42 = buffer.data(dg_s + 42);
    const auto *dg_s_43 = buffer.data(dg_s + 43);
    const auto *dg_s_44 = buffer.data(dg_s + 44);
    const auto *dg_s_45 = buffer.data(dg_s + 45);
    const auto *dg_s_46 = buffer.data(dg_s + 46);
    const auto *dg_s_47 = buffer.data(dg_s + 47);
    const auto *dg_s_48 = buffer.data(dg_s + 48);

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
    const auto *dd_11 = buffer.data(dd + 11);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pf_0, dd_s_0, dg_s_0, dg_s_1, \
                         dg_s_2, dg_s_3, dd_0, df_0, df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 - f_1 * dd_s_0[k]
                 + f_2 * dg_s_0[k]
                 + f_3 * dd_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_2 * dg_s_1[k]
                 + pb_y[k] * df_0[k];

        t_2[k] = f_2 * dg_s_2[k]
                 + pb_z[k] * df_0[k];

        t_3[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_3[k]
                 + f_5 * dd_0[k]
                 + pb_y[k] * df_1[k];
    }

#pragma omp simd aligned(t_4, t_5, pb_y, pb_z, dd_s_0, dd_s_1, dg_s_4, dg_s_5, dd_0, dd_1, \
                         df_2, df_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_4[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_2[k];

        t_5[k] = -f_1 * dd_s_1[k]
                 + f_2 * dg_s_5[k]
                 + f_3 * dd_1[k]
                 + pb_y[k] * df_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, pg_0, dd_s_2, dg_s_6, dg_s_7, \
                         dg_s_8, dg_s_9, dd_2, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * dd_s_2[k]
                 + f_2 * dg_s_6[k]
                 + f_5 * dd_2[k]
                 + pb_y[k] * df_4[k];

        t_7[k] = f_2 * dg_s_7[k]
                 + pb_y[k] * df_5[k];

        t_8[k] = -f_1 * dd_s_2[k]
                 + f_2 * dg_s_8[k]
                 + f_3 * dd_2[k]
                 + pb_z[k] * df_5[k];

        t_9[k] = pa_y[k] * pg_0[k]
                 + f_2 * dg_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_z, pg_0, pg_5, pg_14, dg_s_11, dg_s_12, \
                         dg_s_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * pg_5[k]
                  + f_2 * dg_s_11[k];

        t_11[k] = pa_z[k] * pg_0[k]
                  + f_2 * dg_s_12[k];

        t_12[k] = pa_x[k] * pg_14[k]
                  + f_2 * dg_s_16[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, dd_s_4, dd_s_5, dd_s_6, dg_s_17, dg_s_18, \
                         dg_s_19, dd_3, dd_4, dd_5, df_6, df_7, df_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_1 * dd_s_4[k]
                  + f_2 * dg_s_17[k]
                  + f_3 * dd_3[k]
                  + pb_x[k] * df_6[k];

        t_14[k] = -f_6 * dd_s_5[k]
                  + f_2 * dg_s_18[k]
                  + f_0 * dd_4[k]
                  + pb_x[k] * df_7[k];

        t_15[k] = -f_4 * dd_s_6[k]
                  + f_2 * dg_s_19[k]
                  + f_5 * dd_5[k]
                  + pb_x[k] * df_8[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_x, dd_s_7, dg_s_20, dg_s_21, dg_s_22, \
                         dg_s_23, dd_6, df_9, df_10, df_12, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_4 * dd_s_7[k]
                  + f_2 * dg_s_20[k]
                  + f_5 * dd_6[k]
                  + pb_x[k] * df_9[k];

        t_17[k] = f_2 * dg_s_21[k]
                  + pb_x[k] * df_10[k];

        t_18[k] = f_2 * dg_s_22[k]
                  + pb_x[k] * df_12[k];

        t_19[k] = f_2 * dg_s_23[k]
                  + pb_x[k] * df_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_y, pb_z, pf_1, dd_s_6, dg_s_24, dg_s_25, \
                         dg_s_26, dd_5, df_10, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * pf_1[k]
                  - f_1 * dd_s_6[k]
                  + f_2 * dg_s_24[k]
                  + f_3 * dd_5[k]
                  + pb_y[k] * df_10[k];

        t_21[k] = f_2 * dg_s_25[k]
                  + pb_z[k] * df_10[k];

        t_22[k] = -f_4 * dd_s_6[k]
                  + f_2 * dg_s_26[k]
                  + f_5 * dd_5[k]
                  + pb_z[k] * df_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_z, pb_y, pb_z, pf_2, pg_5, dd_s_7, dg_s_27, \
                         dg_s_28, dg_s_32, dd_6, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * pf_2[k]
                  + f_2 * dg_s_27[k]
                  + pb_y[k] * df_13[k];

        t_24[k] = -f_1 * dd_s_7[k]
                  + f_2 * dg_s_28[k]
                  + f_3 * dd_6[k]
                  + pb_z[k] * df_13[k];

        t_25[k] = pa_z[k] * pg_5[k]
                  + f_2 * dg_s_32[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_x, pg_14, dd_s_11, dd_s_12, dg_s_36, \
                         dg_s_37, dg_s_38, dd_7, dd_8, df_14, df_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_y[k] * pg_14[k]
                  + f_2 * dg_s_36[k];

        t_27[k] = -f_1 * dd_s_11[k]
                  + f_2 * dg_s_37[k]
                  + f_3 * dd_7[k]
                  + pb_x[k] * df_14[k];

        t_28[k] = -f_6 * dd_s_12[k]
                  + f_2 * dg_s_38[k]
                  + f_0 * dd_8[k]
                  + pb_x[k] * df_15[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, dd_s_13, dd_s_15, dg_s_39, dg_s_40, dg_s_41, \
                         dd_9, dd_11, df_16, df_17, df_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -f_4 * dd_s_13[k]
                  + f_2 * dg_s_39[k]
                  + f_5 * dd_9[k]
                  + pb_x[k] * df_16[k];

        t_30[k] = -f_4 * dd_s_15[k]
                  + f_2 * dg_s_40[k]
                  + f_5 * dd_11[k]
                  + pb_x[k] * df_17[k];

        t_31[k] = f_2 * dg_s_41[k]
                  + pb_x[k] * df_18[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, dd_s_13, dg_s_42, dg_s_43, dg_s_44, \
                         dd_9, df_18, df_19, df_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_2 * dg_s_42[k]
                  + pb_x[k] * df_19[k];

        t_33[k] = f_2 * dg_s_43[k]
                  + pb_x[k] * df_21[k];

        t_34[k] = -f_1 * dd_s_13[k]
                  + f_2 * dg_s_44[k]
                  + f_3 * dd_9[k]
                  + pb_y[k] * df_18[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_y, dd_s_14, dd_s_15, dg_s_45, dg_s_46, dg_s_47, \
                         dd_10, dd_11, df_19, df_20, df_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -f_6 * dd_s_14[k]
                  + f_2 * dg_s_45[k]
                  + f_0 * dd_10[k]
                  + pb_y[k] * df_19[k];

        t_36[k] = -f_4 * dd_s_15[k]
                  + f_2 * dg_s_46[k]
                  + f_5 * dd_11[k]
                  + pb_y[k] * df_20[k];

        t_37[k] = f_2 * dg_s_47[k]
                  + pb_y[k] * df_21[k];
    }

#pragma omp simd aligned(t_38, pb_z, pf_3, dd_s_15, dg_s_48, dd_11, \
                         df_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * pf_3[k]
                  - f_1 * dd_s_15[k]
                  + f_2 * dg_s_48[k]
                  + f_3 * dd_11[k]
                  + pb_z[k] * df_21[k];
    }
}

auto
compute_prim_dg_kinetic_energy_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf, const size_t pg,
                                 const size_t dd_s, const size_t dg_s, const size_t dd,
                                 const size_t df, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;

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
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_10 = buffer.data(pf + 10);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_7 = buffer.data(dd_s + 7);
    const auto *dd_s_8 = buffer.data(dd_s + 8);
    const auto *dd_s_9 = buffer.data(dd_s + 9);
    const auto *dd_s_10 = buffer.data(dd_s + 10);
    const auto *dd_s_13 = buffer.data(dd_s + 13);
    const auto *dd_s_14 = buffer.data(dd_s + 14);
    const auto *dd_s_15 = buffer.data(dd_s + 15);
    const auto *dd_s_16 = buffer.data(dd_s + 16);
    const auto *dd_s_17 = buffer.data(dd_s + 17);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_1 = buffer.data(dg_s + 1);
    const auto *dg_s_2 = buffer.data(dg_s + 2);
    const auto *dg_s_3 = buffer.data(dg_s + 3);
    const auto *dg_s_4 = buffer.data(dg_s + 4);
    const auto *dg_s_5 = buffer.data(dg_s + 5);
    const auto *dg_s_6 = buffer.data(dg_s + 6);
    const auto *dg_s_7 = buffer.data(dg_s + 7);
    const auto *dg_s_8 = buffer.data(dg_s + 8);
    const auto *dg_s_9 = buffer.data(dg_s + 9);
    const auto *dg_s_10 = buffer.data(dg_s + 10);
    const auto *dg_s_11 = buffer.data(dg_s + 11);
    const auto *dg_s_12 = buffer.data(dg_s + 12);
    const auto *dg_s_13 = buffer.data(dg_s + 13);
    const auto *dg_s_14 = buffer.data(dg_s + 14);
    const auto *dg_s_15 = buffer.data(dg_s + 15);
    const auto *dg_s_16 = buffer.data(dg_s + 16);
    const auto *dg_s_17 = buffer.data(dg_s + 17);
    const auto *dg_s_18 = buffer.data(dg_s + 18);
    const auto *dg_s_19 = buffer.data(dg_s + 19);
    const auto *dg_s_20 = buffer.data(dg_s + 20);
    const auto *dg_s_21 = buffer.data(dg_s + 21);
    const auto *dg_s_22 = buffer.data(dg_s + 22);
    const auto *dg_s_23 = buffer.data(dg_s + 23);
    const auto *dg_s_24 = buffer.data(dg_s + 24);
    const auto *dg_s_25 = buffer.data(dg_s + 25);
    const auto *dg_s_26 = buffer.data(dg_s + 26);
    const auto *dg_s_27 = buffer.data(dg_s + 27);
    const auto *dg_s_28 = buffer.data(dg_s + 28);
    const auto *dg_s_29 = buffer.data(dg_s + 29);
    const auto *dg_s_30 = buffer.data(dg_s + 30);
    const auto *dg_s_31 = buffer.data(dg_s + 31);
    const auto *dg_s_32 = buffer.data(dg_s + 32);
    const auto *dg_s_33 = buffer.data(dg_s + 33);
    const auto *dg_s_34 = buffer.data(dg_s + 34);
    const auto *dg_s_35 = buffer.data(dg_s + 35);
    const auto *dg_s_36 = buffer.data(dg_s + 36);
    const auto *dg_s_37 = buffer.data(dg_s + 37);
    const auto *dg_s_38 = buffer.data(dg_s + 38);
    const auto *dg_s_39 = buffer.data(dg_s + 39);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_15 = buffer.data(dd + 15);
    const auto *dd_16 = buffer.data(dd + 16);
    const auto *dd_17 = buffer.data(dd + 17);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_26 = buffer.data(df + 26);
    const auto *df_27 = buffer.data(df + 27);
    const auto *df_28 = buffer.data(df + 28);
    const auto *df_29 = buffer.data(df + 29);
    const auto *df_30 = buffer.data(df + 30);
    const auto *df_31 = buffer.data(df + 31);
    const auto *df_32 = buffer.data(df + 32);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pf_0, dd_s_0, dg_s_0, dg_s_1, \
                         dg_s_2, dd_0, df_0, df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 - f_1 * dd_s_0[k]
                 + f_2 * dg_s_0[k]
                 + f_3 * dd_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_1[k]
                 + f_5 * dd_0[k]
                 + pb_y[k] * df_1[k];

        t_2[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_2[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, pf_1, pf_2, dd_s_1, dg_s_3, dg_s_4, \
                         dg_s_5, dd_1, df_3, df_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_0 * pf_1[k]
                 + f_2 * dg_s_3[k]
                 + pb_x[k] * df_3[k];

        t_4[k] = f_0 * pf_2[k]
                 + f_2 * dg_s_4[k]
                 + pb_x[k] * df_4[k];

        t_5[k] = -f_1 * dd_s_1[k]
                 + f_2 * dg_s_5[k]
                 + f_3 * dd_1[k]
                 + pb_y[k] * df_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pb_y, pb_z, pf_0, pf_3, pg_1, dd_s_2, dg_s_6, \
                         dg_s_7, dg_s_8, dd_2, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_1 * dd_s_2[k]
                 + f_2 * dg_s_6[k]
                 + f_3 * dd_2[k]
                 + pb_z[k] * df_4[k];

        t_7[k] = f_5 * pf_0[k]
                 + f_2 * dg_s_7[k]
                 + pb_y[k] * df_5[k];

        t_8[k] = f_0 * pf_3[k]
                 + pa_x[k] * pg_1[k]
                 + f_2 * dg_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_z, pb_x, pf_4, pg_0, pg_2, dg_s_9, dg_s_10, \
                         dg_s_11, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * pf_4[k]
                 + f_2 * dg_s_9[k]
                 + pb_x[k] * df_6[k];

        t_10[k] = pa_x[k] * pg_2[k]
                  + f_2 * dg_s_10[k];

        t_11[k] = pa_z[k] * pg_0[k]
                  + f_2 * dg_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_x, pb_z, pf_0, pf_7, pf_10, pg_3, dg_s_12, \
                         dg_s_13, dg_s_14, df_8, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * pf_0[k]
                  + f_2 * dg_s_12[k]
                  + pb_z[k] * df_8[k];

        t_13[k] = f_0 * pf_7[k]
                  + pa_x[k] * pg_3[k]
                  + f_2 * dg_s_13[k];

        t_14[k] = f_5 * pf_10[k]
                  + f_2 * dg_s_14[k]
                  + pb_x[k] * df_11[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pb_x, pg_5, dd_s_7, dd_s_8, dg_s_15, dg_s_16, \
                         dg_s_17, dd_7, dd_8, df_12, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_x[k] * pg_5[k]
                  + f_2 * dg_s_15[k];

        t_16[k] = -f_1 * dd_s_7[k]
                  + f_2 * dg_s_16[k]
                  + f_3 * dd_7[k]
                  + pb_x[k] * df_12[k];

        t_17[k] = -f_6 * dd_s_8[k]
                  + f_2 * dg_s_17[k]
                  + f_0 * dd_8[k]
                  + pb_x[k] * df_13[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_x, pb_z, dd_s_9, dd_s_10, dg_s_18, dg_s_19, \
                         dg_s_20, dd_9, dd_10, df_13, df_14, df_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -f_4 * dd_s_9[k]
                  + f_2 * dg_s_18[k]
                  + f_5 * dd_9[k]
                  + pb_x[k] * df_14[k];

        t_19[k] = f_2 * dg_s_19[k]
                  + pb_z[k] * df_13[k];

        t_20[k] = -f_4 * dd_s_10[k]
                  + f_2 * dg_s_20[k]
                  + f_5 * dd_10[k]
                  + pb_x[k] * df_15[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_y, pb_z, pf_4, pf_6, dd_s_9, dg_s_21, dg_s_22, \
                         dg_s_23, dd_9, df_16, df_17, df_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * pf_4[k]
                  - f_1 * dd_s_9[k]
                  + f_2 * dg_s_21[k]
                  + f_3 * dd_9[k]
                  + pb_y[k] * df_16[k];

        t_22[k] = -f_4 * dd_s_9[k]
                  + f_2 * dg_s_22[k]
                  + f_5 * dd_9[k]
                  + pb_z[k] * df_17[k];

        t_23[k] = f_0 * pf_6[k]
                  + f_2 * dg_s_23[k]
                  + pb_y[k] * df_19[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_z, pb_z, pf_4, pg_2, dd_s_10, dg_s_24, dg_s_25, \
                         dg_s_26, dd_10, df_19, df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_1 * dd_s_10[k]
                  + f_2 * dg_s_24[k]
                  + f_3 * dd_10[k]
                  + pb_z[k] * df_19[k];

        t_25[k] = pa_z[k] * pg_2[k]
                  + f_2 * dg_s_25[k];

        t_26[k] = f_5 * pf_4[k]
                  + f_2 * dg_s_26[k]
                  + pb_z[k] * df_20[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_y, pf_9, pf_10, pg_4, pg_5, dg_s_27, \
                         dg_s_28, dg_s_29, df_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * pf_9[k]
                  + pa_y[k] * pg_4[k]
                  + f_2 * dg_s_27[k];

        t_28[k] = f_5 * pf_10[k]
                  + f_2 * dg_s_28[k]
                  + pb_y[k] * df_23[k];

        t_29[k] = pa_y[k] * pg_5[k]
                  + f_2 * dg_s_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_x, pb_y, dd_s_13, dd_s_14, dg_s_30, dg_s_31, \
                         dg_s_32, dd_13, dd_14, df_24, df_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -f_1 * dd_s_13[k]
                  + f_2 * dg_s_30[k]
                  + f_3 * dd_13[k]
                  + pb_x[k] * df_24[k];

        t_31[k] = f_2 * dg_s_31[k]
                  + pb_y[k] * df_24[k];

        t_32[k] = -f_6 * dd_s_14[k]
                  + f_2 * dg_s_32[k]
                  + f_0 * dd_14[k]
                  + pb_x[k] * df_26[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_x, pb_y, dd_s_15, dd_s_17, dg_s_33, dg_s_34, \
                         dg_s_35, dd_15, dd_17, df_26, df_27, df_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -f_4 * dd_s_15[k]
                  + f_2 * dg_s_33[k]
                  + f_5 * dd_15[k]
                  + pb_x[k] * df_27[k];

        t_34[k] = f_2 * dg_s_34[k]
                  + pb_y[k] * df_26[k];

        t_35[k] = -f_4 * dd_s_17[k]
                  + f_2 * dg_s_35[k]
                  + f_5 * dd_17[k]
                  + pb_x[k] * df_28[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_y, dd_s_15, dd_s_16, dd_s_17, dg_s_36, dg_s_37, \
                         dg_s_38, dd_15, dd_16, dd_17, df_29, df_30, \
                         df_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_1 * dd_s_15[k]
                  + f_2 * dg_s_36[k]
                  + f_3 * dd_15[k]
                  + pb_y[k] * df_29[k];

        t_37[k] = -f_6 * dd_s_16[k]
                  + f_2 * dg_s_37[k]
                  + f_0 * dd_16[k]
                  + pb_y[k] * df_30[k];

        t_38[k] = -f_4 * dd_s_17[k]
                  + f_2 * dg_s_38[k]
                  + f_5 * dd_17[k]
                  + pb_y[k] * df_31[k];
    }

#pragma omp simd aligned(t_39, pb_z, pf_10, dd_s_17, dg_s_39, dd_17, \
                         df_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * pf_10[k]
                  - f_1 * dd_s_17[k]
                  + f_2 * dg_s_39[k]
                  + f_3 * dd_17[k]
                  + pb_z[k] * df_32[k];
    }
}

auto
compute_prim_dg_kinetic_energy_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf, const size_t pg,
                                 const size_t dd_s, const size_t dg_s, const size_t dd,
                                 const size_t df, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;

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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_15 = buffer.data(pg + 15);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_6 = buffer.data(dd_s + 6);
    const auto *dd_s_7 = buffer.data(dd_s + 7);
    const auto *dd_s_8 = buffer.data(dd_s + 8);
    const auto *dd_s_9 = buffer.data(dd_s + 9);
    const auto *dd_s_12 = buffer.data(dd_s + 12);
    const auto *dd_s_13 = buffer.data(dd_s + 13);
    const auto *dd_s_14 = buffer.data(dd_s + 14);
    const auto *dd_s_15 = buffer.data(dd_s + 15);
    const auto *dd_s_16 = buffer.data(dd_s + 16);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_1 = buffer.data(dg_s + 1);
    const auto *dg_s_2 = buffer.data(dg_s + 2);
    const auto *dg_s_3 = buffer.data(dg_s + 3);
    const auto *dg_s_4 = buffer.data(dg_s + 4);
    const auto *dg_s_5 = buffer.data(dg_s + 5);
    const auto *dg_s_6 = buffer.data(dg_s + 6);
    const auto *dg_s_7 = buffer.data(dg_s + 7);
    const auto *dg_s_8 = buffer.data(dg_s + 8);
    const auto *dg_s_9 = buffer.data(dg_s + 9);
    const auto *dg_s_10 = buffer.data(dg_s + 10);
    const auto *dg_s_11 = buffer.data(dg_s + 11);
    const auto *dg_s_12 = buffer.data(dg_s + 12);
    const auto *dg_s_13 = buffer.data(dg_s + 13);
    const auto *dg_s_14 = buffer.data(dg_s + 14);
    const auto *dg_s_15 = buffer.data(dg_s + 15);
    const auto *dg_s_16 = buffer.data(dg_s + 16);
    const auto *dg_s_17 = buffer.data(dg_s + 17);
    const auto *dg_s_18 = buffer.data(dg_s + 18);
    const auto *dg_s_19 = buffer.data(dg_s + 19);
    const auto *dg_s_20 = buffer.data(dg_s + 20);
    const auto *dg_s_21 = buffer.data(dg_s + 21);
    const auto *dg_s_22 = buffer.data(dg_s + 22);
    const auto *dg_s_23 = buffer.data(dg_s + 23);
    const auto *dg_s_24 = buffer.data(dg_s + 24);
    const auto *dg_s_25 = buffer.data(dg_s + 25);
    const auto *dg_s_26 = buffer.data(dg_s + 26);
    const auto *dg_s_27 = buffer.data(dg_s + 27);
    const auto *dg_s_28 = buffer.data(dg_s + 28);
    const auto *dg_s_29 = buffer.data(dg_s + 29);
    const auto *dg_s_30 = buffer.data(dg_s + 30);
    const auto *dg_s_31 = buffer.data(dg_s + 31);
    const auto *dg_s_32 = buffer.data(dg_s + 32);
    const auto *dg_s_33 = buffer.data(dg_s + 33);
    const auto *dg_s_34 = buffer.data(dg_s + 34);
    const auto *dg_s_35 = buffer.data(dg_s + 35);
    const auto *dg_s_36 = buffer.data(dg_s + 36);
    const auto *dg_s_37 = buffer.data(dg_s + 37);
    const auto *dg_s_38 = buffer.data(dg_s + 38);
    const auto *dg_s_39 = buffer.data(dg_s + 39);
    const auto *dg_s_40 = buffer.data(dg_s + 40);
    const auto *dg_s_41 = buffer.data(dg_s + 41);
    const auto *dg_s_42 = buffer.data(dg_s + 42);
    const auto *dg_s_43 = buffer.data(dg_s + 43);
    const auto *dg_s_44 = buffer.data(dg_s + 44);
    const auto *dg_s_45 = buffer.data(dg_s + 45);
    const auto *dg_s_46 = buffer.data(dg_s + 46);
    const auto *dg_s_47 = buffer.data(dg_s + 47);
    const auto *dg_s_48 = buffer.data(dg_s + 48);
    const auto *dg_s_49 = buffer.data(dg_s + 49);
    const auto *dg_s_50 = buffer.data(dg_s + 50);
    const auto *dg_s_51 = buffer.data(dg_s + 51);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_8 = buffer.data(df + 8);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pf_0, dd_s_0, dg_s_0, dg_s_1, \
                         dg_s_2, dg_s_3, dd_0, df_0, df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 - f_1 * dd_s_0[k]
                 + f_2 * dg_s_0[k]
                 + f_3 * dd_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_2 * dg_s_1[k]
                 + pb_y[k] * df_0[k];

        t_2[k] = f_2 * dg_s_2[k]
                 + pb_z[k] * df_0[k];

        t_3[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_3[k]
                 + f_5 * dd_0[k]
                 + pb_y[k] * df_1[k];
    }

#pragma omp simd aligned(t_4, t_5, pb_y, pb_z, dd_s_0, dd_s_1, dg_s_4, dg_s_5, dd_0, dd_1, \
                         df_2, df_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_4[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_2[k];

        t_5[k] = -f_1 * dd_s_1[k]
                 + f_2 * dg_s_5[k]
                 + f_3 * dd_1[k]
                 + pb_y[k] * df_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_y, pb_y, pb_z, pg_0, dd_s_2, dg_s_6, dg_s_7, \
                         dg_s_8, dd_2, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * dd_s_2[k]
                 + f_2 * dg_s_6[k]
                 + f_5 * dd_2[k]
                 + pb_y[k] * df_4[k];

        t_7[k] = -f_1 * dd_s_2[k]
                 + f_2 * dg_s_7[k]
                 + f_3 * dd_2[k]
                 + pb_z[k] * df_5[k];

        t_8[k] = pa_y[k] * pg_0[k]
                 + f_2 * dg_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pa_y, pf_1, pg_1, pg_4, pg_6, pg_7, \
                         dg_s_9, dg_s_10, dg_s_11, dg_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * pf_1[k]
                 + pa_x[k] * pg_4[k]
                 + f_2 * dg_s_9[k];

        t_10[k] = pa_y[k] * pg_1[k]
                  + f_2 * dg_s_10[k];

        t_11[k] = pa_x[k] * pg_6[k]
                  + f_2 * dg_s_11[k];

        t_12[k] = pa_x[k] * pg_7[k]
                  + f_2 * dg_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pa_z, pb_z, pf_0, pg_0, pg_8, dg_s_13, \
                         dg_s_14, dg_s_15, df_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_x[k] * pg_8[k]
                  + f_2 * dg_s_13[k];

        t_14[k] = pa_z[k] * pg_0[k]
                  + f_2 * dg_s_14[k];

        t_15[k] = f_5 * pf_0[k]
                  + f_2 * dg_s_15[k]
                  + pb_z[k] * df_8[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pf_4, pg_10, pg_12, pg_13, pg_15, \
                         dg_s_16, dg_s_17, dg_s_18, dg_s_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * pf_4[k]
                  + pa_x[k] * pg_10[k]
                  + f_2 * dg_s_16[k];

        t_17[k] = pa_x[k] * pg_12[k]
                  + f_2 * dg_s_17[k];

        t_18[k] = pa_x[k] * pg_13[k]
                  + f_2 * dg_s_18[k];

        t_19[k] = pa_x[k] * pg_15[k]
                  + f_2 * dg_s_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_x, dd_s_6, dd_s_7, dd_s_8, dg_s_20, dg_s_21, \
                         dg_s_22, dd_5, dd_6, dd_7, df_10, df_11, \
                         df_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_1 * dd_s_6[k]
                  + f_2 * dg_s_20[k]
                  + f_3 * dd_5[k]
                  + pb_x[k] * df_10[k];

        t_21[k] = -f_6 * dd_s_7[k]
                  + f_2 * dg_s_21[k]
                  + f_0 * dd_6[k]
                  + pb_x[k] * df_11[k];

        t_22[k] = -f_4 * dd_s_8[k]
                  + f_2 * dg_s_22[k]
                  + f_5 * dd_7[k]
                  + pb_x[k] * df_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_x, pb_z, dd_s_9, dg_s_23, dg_s_24, \
                         dg_s_25, dg_s_26, dd_8, df_11, df_13, df_14, \
                         df_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_2 * dg_s_23[k]
                  + pb_z[k] * df_11[k];

        t_24[k] = -f_4 * dd_s_9[k]
                  + f_2 * dg_s_24[k]
                  + f_5 * dd_8[k]
                  + pb_x[k] * df_13[k];

        t_25[k] = f_2 * dg_s_25[k]
                  + pb_x[k] * df_14[k];

        t_26[k] = f_2 * dg_s_26[k]
                  + pb_x[k] * df_16[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_x, pb_y, pb_z, pf_2, dd_s_8, dg_s_27, dg_s_28, \
                         dg_s_29, dd_7, df_14, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_2 * dg_s_27[k]
                  + pb_x[k] * df_17[k];

        t_28[k] = f_0 * pf_2[k]
                  - f_1 * dd_s_8[k]
                  + f_2 * dg_s_28[k]
                  + f_3 * dd_7[k]
                  + pb_y[k] * df_14[k];

        t_29[k] = f_2 * dg_s_29[k]
                  + pb_z[k] * df_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_y, pb_z, pf_3, dd_s_8, dd_s_9, dg_s_30, dg_s_31, \
                         dg_s_32, dd_7, dd_8, df_15, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -f_4 * dd_s_8[k]
                  + f_2 * dg_s_30[k]
                  + f_5 * dd_7[k]
                  + pb_z[k] * df_15[k];

        t_31[k] = f_0 * pf_3[k]
                  + f_2 * dg_s_31[k]
                  + pb_y[k] * df_17[k];

        t_32[k] = -f_1 * dd_s_9[k]
                  + f_2 * dg_s_32[k]
                  + f_3 * dd_8[k]
                  + pb_z[k] * df_17[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pa_z, pb_z, pf_2, pf_5, pg_6, pg_13, dg_s_33, \
                         dg_s_34, dg_s_35, df_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * pg_6[k]
                  + f_2 * dg_s_33[k];

        t_34[k] = f_5 * pf_2[k]
                  + f_2 * dg_s_34[k]
                  + pb_z[k] * df_18[k];

        t_35[k] = f_0 * pf_5[k]
                  + pa_y[k] * pg_13[k]
                  + f_2 * dg_s_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_y, pb_x, pb_y, pf_6, pg_15, dd_s_12, dg_s_36, \
                         dg_s_37, dg_s_38, dd_10, df_19, df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_5 * pf_6[k]
                  + f_2 * dg_s_36[k]
                  + pb_y[k] * df_19[k];

        t_37[k] = pa_y[k] * pg_15[k]
                  + f_2 * dg_s_37[k];

        t_38[k] = -f_1 * dd_s_12[k]
                  + f_2 * dg_s_38[k]
                  + f_3 * dd_10[k]
                  + pb_x[k] * df_20[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_x, pb_y, dd_s_13, dd_s_14, dg_s_39, dg_s_40, \
                         dg_s_41, dd_11, dd_12, df_20, df_21, df_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_2 * dg_s_39[k]
                  + pb_y[k] * df_20[k];

        t_40[k] = -f_6 * dd_s_13[k]
                  + f_2 * dg_s_40[k]
                  + f_0 * dd_11[k]
                  + pb_x[k] * df_21[k];

        t_41[k] = -f_4 * dd_s_14[k]
                  + f_2 * dg_s_41[k]
                  + f_5 * dd_12[k]
                  + pb_x[k] * df_22[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_x, pb_y, dd_s_16, dg_s_42, dg_s_43, \
                         dg_s_44, dg_s_45, dd_14, df_21, df_23, df_24, \
                         df_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_2 * dg_s_42[k]
                  + pb_y[k] * df_21[k];

        t_43[k] = -f_4 * dd_s_16[k]
                  + f_2 * dg_s_43[k]
                  + f_5 * dd_14[k]
                  + pb_x[k] * df_23[k];

        t_44[k] = f_2 * dg_s_44[k]
                  + pb_x[k] * df_24[k];

        t_45[k] = f_2 * dg_s_45[k]
                  + pb_x[k] * df_25[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_x, pb_y, dd_s_14, dd_s_15, dg_s_46, dg_s_47, \
                         dg_s_48, dd_12, dd_13, df_24, df_25, df_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_2 * dg_s_46[k]
                  + pb_x[k] * df_27[k];

        t_47[k] = -f_1 * dd_s_14[k]
                  + f_2 * dg_s_47[k]
                  + f_3 * dd_12[k]
                  + pb_y[k] * df_24[k];

        t_48[k] = -f_6 * dd_s_15[k]
                  + f_2 * dg_s_48[k]
                  + f_0 * dd_13[k]
                  + pb_y[k] * df_25[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_y, pb_z, pf_6, dd_s_16, dg_s_49, dg_s_50, \
                         dg_s_51, dd_14, df_26, df_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = -f_4 * dd_s_16[k]
                  + f_2 * dg_s_49[k]
                  + f_5 * dd_14[k]
                  + pb_y[k] * df_26[k];

        t_50[k] = f_2 * dg_s_50[k]
                  + pb_y[k] * df_27[k];

        t_51[k] = f_0 * pf_6[k]
                  - f_1 * dd_s_16[k]
                  + f_2 * dg_s_51[k]
                  + f_3 * dd_14[k]
                  + pb_z[k] * df_27[k];
    }
}

auto
compute_prim_dg_kinetic_energy_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf, const size_t pg,
                                 const size_t dd_s, const size_t dg_s, const size_t dd,
                                 const size_t df, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;

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
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_9 = buffer.data(pg + 9);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_6 = buffer.data(dd_s + 6);
    const auto *dd_s_7 = buffer.data(dd_s + 7);
    const auto *dd_s_8 = buffer.data(dd_s + 8);
    const auto *dd_s_9 = buffer.data(dd_s + 9);
    const auto *dd_s_12 = buffer.data(dd_s + 12);
    const auto *dd_s_13 = buffer.data(dd_s + 13);
    const auto *dd_s_14 = buffer.data(dd_s + 14);
    const auto *dd_s_15 = buffer.data(dd_s + 15);
    const auto *dd_s_16 = buffer.data(dd_s + 16);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_1 = buffer.data(dg_s + 1);
    const auto *dg_s_2 = buffer.data(dg_s + 2);
    const auto *dg_s_3 = buffer.data(dg_s + 3);
    const auto *dg_s_4 = buffer.data(dg_s + 4);
    const auto *dg_s_5 = buffer.data(dg_s + 5);
    const auto *dg_s_6 = buffer.data(dg_s + 6);
    const auto *dg_s_7 = buffer.data(dg_s + 7);
    const auto *dg_s_8 = buffer.data(dg_s + 8);
    const auto *dg_s_9 = buffer.data(dg_s + 9);
    const auto *dg_s_10 = buffer.data(dg_s + 10);
    const auto *dg_s_11 = buffer.data(dg_s + 11);
    const auto *dg_s_12 = buffer.data(dg_s + 12);
    const auto *dg_s_13 = buffer.data(dg_s + 13);
    const auto *dg_s_14 = buffer.data(dg_s + 14);
    const auto *dg_s_15 = buffer.data(dg_s + 15);
    const auto *dg_s_16 = buffer.data(dg_s + 16);
    const auto *dg_s_17 = buffer.data(dg_s + 17);
    const auto *dg_s_18 = buffer.data(dg_s + 18);
    const auto *dg_s_19 = buffer.data(dg_s + 19);
    const auto *dg_s_20 = buffer.data(dg_s + 20);
    const auto *dg_s_21 = buffer.data(dg_s + 21);
    const auto *dg_s_22 = buffer.data(dg_s + 22);
    const auto *dg_s_23 = buffer.data(dg_s + 23);
    const auto *dg_s_24 = buffer.data(dg_s + 24);
    const auto *dg_s_25 = buffer.data(dg_s + 25);
    const auto *dg_s_26 = buffer.data(dg_s + 26);
    const auto *dg_s_27 = buffer.data(dg_s + 27);
    const auto *dg_s_28 = buffer.data(dg_s + 28);
    const auto *dg_s_29 = buffer.data(dg_s + 29);
    const auto *dg_s_30 = buffer.data(dg_s + 30);
    const auto *dg_s_31 = buffer.data(dg_s + 31);
    const auto *dg_s_32 = buffer.data(dg_s + 32);
    const auto *dg_s_33 = buffer.data(dg_s + 33);
    const auto *dg_s_34 = buffer.data(dg_s + 34);
    const auto *dg_s_35 = buffer.data(dg_s + 35);
    const auto *dg_s_36 = buffer.data(dg_s + 36);
    const auto *dg_s_37 = buffer.data(dg_s + 37);
    const auto *dg_s_38 = buffer.data(dg_s + 38);
    const auto *dg_s_39 = buffer.data(dg_s + 39);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pf_0, dd_s_0, dg_s_0, dg_s_1, \
                         dg_s_2, dg_s_3, dd_0, df_0, df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 - f_1 * dd_s_0[k]
                 + f_2 * dg_s_0[k]
                 + f_3 * dd_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_2 * dg_s_1[k]
                 + pb_y[k] * df_0[k];

        t_2[k] = f_2 * dg_s_2[k]
                 + pb_z[k] * df_0[k];

        t_3[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_3[k]
                 + f_5 * dd_0[k]
                 + pb_y[k] * df_1[k];
    }

#pragma omp simd aligned(t_4, t_5, pb_y, pb_z, dd_s_0, dd_s_1, dg_s_4, dg_s_5, dd_0, dd_1, \
                         df_2, df_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_4[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_2[k];

        t_5[k] = -f_1 * dd_s_1[k]
                 + f_2 * dg_s_5[k]
                 + f_3 * dd_1[k]
                 + pb_y[k] * df_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, pg_0, dd_s_2, dg_s_6, dg_s_7, \
                         dg_s_8, dg_s_9, dd_2, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * dd_s_2[k]
                 + f_2 * dg_s_6[k]
                 + f_5 * dd_2[k]
                 + pb_y[k] * df_4[k];

        t_7[k] = f_2 * dg_s_7[k]
                 + pb_y[k] * df_5[k];

        t_8[k] = -f_1 * dd_s_2[k]
                 + f_2 * dg_s_8[k]
                 + f_3 * dd_2[k]
                 + pb_z[k] * df_5[k];

        t_9[k] = pa_y[k] * pg_0[k]
                 + f_2 * dg_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_z, pg_0, pg_3, pg_9, dg_s_10, dg_s_11, \
                         dg_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * pg_3[k]
                  + f_2 * dg_s_10[k];

        t_11[k] = pa_z[k] * pg_0[k]
                  + f_2 * dg_s_11[k];

        t_12[k] = pa_x[k] * pg_9[k]
                  + f_2 * dg_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, dd_s_6, dd_s_7, dd_s_8, dg_s_13, dg_s_14, \
                         dg_s_15, dd_5, dd_6, dd_7, df_8, df_9, df_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_1 * dd_s_6[k]
                  + f_2 * dg_s_13[k]
                  + f_3 * dd_5[k]
                  + pb_x[k] * df_8[k];

        t_14[k] = -f_6 * dd_s_7[k]
                  + f_2 * dg_s_14[k]
                  + f_0 * dd_6[k]
                  + pb_x[k] * df_9[k];

        t_15[k] = -f_4 * dd_s_8[k]
                  + f_2 * dg_s_15[k]
                  + f_5 * dd_7[k]
                  + pb_x[k] * df_10[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_x, dd_s_9, dg_s_16, dg_s_17, dg_s_18, \
                         dg_s_19, dd_8, df_11, df_12, df_14, df_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_4 * dd_s_9[k]
                  + f_2 * dg_s_16[k]
                  + f_5 * dd_8[k]
                  + pb_x[k] * df_11[k];

        t_17[k] = f_2 * dg_s_17[k]
                  + pb_x[k] * df_12[k];

        t_18[k] = f_2 * dg_s_18[k]
                  + pb_x[k] * df_14[k];

        t_19[k] = f_2 * dg_s_19[k]
                  + pb_x[k] * df_15[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_y, pb_z, pf_1, dd_s_8, dg_s_20, dg_s_21, \
                         dg_s_22, dd_7, df_12, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * pf_1[k]
                  - f_1 * dd_s_8[k]
                  + f_2 * dg_s_20[k]
                  + f_3 * dd_7[k]
                  + pb_y[k] * df_12[k];

        t_21[k] = f_2 * dg_s_21[k]
                  + pb_z[k] * df_12[k];

        t_22[k] = -f_4 * dd_s_8[k]
                  + f_2 * dg_s_22[k]
                  + f_5 * dd_7[k]
                  + pb_z[k] * df_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_z, pb_y, pb_z, pf_2, pg_3, dd_s_9, dg_s_23, \
                         dg_s_24, dg_s_25, dd_8, df_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * pf_2[k]
                  + f_2 * dg_s_23[k]
                  + pb_y[k] * df_15[k];

        t_24[k] = -f_1 * dd_s_9[k]
                  + f_2 * dg_s_24[k]
                  + f_3 * dd_8[k]
                  + pb_z[k] * df_15[k];

        t_25[k] = pa_z[k] * pg_3[k]
                  + f_2 * dg_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_x, pf_3, pg_7, pg_9, dd_s_12, dg_s_26, \
                         dg_s_27, dg_s_28, dd_10, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * pf_3[k]
                  + pa_y[k] * pg_7[k]
                  + f_2 * dg_s_26[k];

        t_27[k] = pa_y[k] * pg_9[k]
                  + f_2 * dg_s_27[k];

        t_28[k] = -f_1 * dd_s_12[k]
                  + f_2 * dg_s_28[k]
                  + f_3 * dd_10[k]
                  + pb_x[k] * df_17[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, dd_s_13, dd_s_14, dd_s_16, dg_s_29, dg_s_30, \
                         dg_s_31, dd_11, dd_12, dd_14, df_18, df_19, \
                         df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -f_6 * dd_s_13[k]
                  + f_2 * dg_s_29[k]
                  + f_0 * dd_11[k]
                  + pb_x[k] * df_18[k];

        t_30[k] = -f_4 * dd_s_14[k]
                  + f_2 * dg_s_30[k]
                  + f_5 * dd_12[k]
                  + pb_x[k] * df_19[k];

        t_31[k] = -f_4 * dd_s_16[k]
                  + f_2 * dg_s_31[k]
                  + f_5 * dd_14[k]
                  + pb_x[k] * df_20[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_x, pb_y, dd_s_14, dg_s_32, dg_s_33, \
                         dg_s_34, dg_s_35, dd_12, df_21, df_22, df_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_2 * dg_s_32[k]
                  + pb_x[k] * df_21[k];

        t_33[k] = f_2 * dg_s_33[k]
                  + pb_x[k] * df_22[k];

        t_34[k] = f_2 * dg_s_34[k]
                  + pb_x[k] * df_24[k];

        t_35[k] = -f_1 * dd_s_14[k]
                  + f_2 * dg_s_35[k]
                  + f_3 * dd_12[k]
                  + pb_y[k] * df_21[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_y, dd_s_15, dd_s_16, dg_s_36, dg_s_37, dg_s_38, \
                         dd_13, dd_14, df_22, df_23, df_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_6 * dd_s_15[k]
                  + f_2 * dg_s_36[k]
                  + f_0 * dd_13[k]
                  + pb_y[k] * df_22[k];

        t_37[k] = -f_4 * dd_s_16[k]
                  + f_2 * dg_s_37[k]
                  + f_5 * dd_14[k]
                  + pb_y[k] * df_23[k];

        t_38[k] = f_2 * dg_s_38[k]
                  + pb_y[k] * df_24[k];
    }

#pragma omp simd aligned(t_39, pb_z, pf_4, dd_s_16, dg_s_39, dd_14, \
                         df_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * pf_4[k]
                  - f_1 * dd_s_16[k]
                  + f_2 * dg_s_39[k]
                  + f_3 * dd_14[k]
                  + pb_z[k] * df_24[k];
    }
}

auto
compute_prim_dg_kinetic_energy_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf, const size_t pg,
                                 const size_t dd_s, const size_t dg_s, const size_t dd,
                                 const size_t df, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_6 = buffer.data(dd_s + 6);
    const auto *dd_s_7 = buffer.data(dd_s + 7);
    const auto *dd_s_11 = buffer.data(dd_s + 11);
    const auto *dd_s_12 = buffer.data(dd_s + 12);
    const auto *dd_s_14 = buffer.data(dd_s + 14);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_1 = buffer.data(dg_s + 1);
    const auto *dg_s_2 = buffer.data(dg_s + 2);
    const auto *dg_s_3 = buffer.data(dg_s + 3);
    const auto *dg_s_4 = buffer.data(dg_s + 4);
    const auto *dg_s_5 = buffer.data(dg_s + 5);
    const auto *dg_s_6 = buffer.data(dg_s + 6);
    const auto *dg_s_7 = buffer.data(dg_s + 7);
    const auto *dg_s_8 = buffer.data(dg_s + 8);
    const auto *dg_s_9 = buffer.data(dg_s + 9);
    const auto *dg_s_10 = buffer.data(dg_s + 10);
    const auto *dg_s_11 = buffer.data(dg_s + 11);
    const auto *dg_s_12 = buffer.data(dg_s + 12);
    const auto *dg_s_13 = buffer.data(dg_s + 13);
    const auto *dg_s_14 = buffer.data(dg_s + 14);
    const auto *dg_s_15 = buffer.data(dg_s + 15);
    const auto *dg_s_16 = buffer.data(dg_s + 16);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_14 = buffer.data(dd + 14);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_26 = buffer.data(df + 26);
    const auto *df_27 = buffer.data(df + 27);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pf_0, dd_s_0, dg_s_0, dg_s_1, \
                         dg_s_2, dd_0, df_0, df_1, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 - f_1 * dd_s_0[k]
                 + f_2 * dg_s_0[k]
                 + f_3 * dd_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_1[k]
                 + f_5 * dd_0[k]
                 + pb_y[k] * df_1[k];

        t_2[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_2[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_x, pg_0, pg_2, dd_s_6, dg_s_3, dg_s_4, \
                         dg_s_5, dd_6, df_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_x[k] * pg_0[k]
                 + f_2 * dg_s_3[k];

        t_4[k] = pa_x[k] * pg_2[k]
                 + f_2 * dg_s_4[k];

        t_5[k] = -f_1 * dd_s_6[k]
                 + f_2 * dg_s_5[k]
                 + f_3 * dd_6[k]
                 + pb_x[k] * df_10[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, pb_z, pf_1, dd_s_7, dg_s_6, dg_s_7, \
                         dg_s_8, dd_7, df_11, df_12, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * dd_s_7[k]
                 + f_2 * dg_s_6[k]
                 + f_5 * dd_7[k]
                 + pb_x[k] * df_11[k];

        t_7[k] = f_0 * pf_1[k]
                 - f_1 * dd_s_7[k]
                 + f_2 * dg_s_7[k]
                 + f_3 * dd_7[k]
                 + pb_y[k] * df_12[k];

        t_8[k] = -f_4 * dd_s_7[k]
                 + f_2 * dg_s_8[k]
                 + f_5 * dd_7[k]
                 + pb_z[k] * df_13[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pf_2, pg_1, pg_2, dd_s_11, dg_s_9, \
                         dg_s_10, dg_s_11, dd_11, df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * pf_2[k]
                 + pa_y[k] * pg_1[k]
                 + f_2 * dg_s_9[k];

        t_10[k] = pa_y[k] * pg_2[k]
                  + f_2 * dg_s_10[k];

        t_11[k] = -f_1 * dd_s_11[k]
                  + f_2 * dg_s_11[k]
                  + f_3 * dd_11[k]
                  + pb_x[k] * df_20[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_x, pb_y, dd_s_12, dd_s_14, dg_s_12, dg_s_13, \
                         dg_s_14, dd_12, dd_14, df_22, df_23, df_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_4 * dd_s_12[k]
                  + f_2 * dg_s_12[k]
                  + f_5 * dd_12[k]
                  + pb_x[k] * df_22[k];

        t_13[k] = -f_4 * dd_s_14[k]
                  + f_2 * dg_s_13[k]
                  + f_5 * dd_14[k]
                  + pb_x[k] * df_23[k];

        t_14[k] = -f_1 * dd_s_12[k]
                  + f_2 * dg_s_14[k]
                  + f_3 * dd_12[k]
                  + pb_y[k] * df_24[k];
    }

#pragma omp simd aligned(t_15, t_16, pb_y, pb_z, pf_3, dd_s_14, dg_s_15, dg_s_16, dd_14, \
                         df_26, df_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_4 * dd_s_14[k]
                  + f_2 * dg_s_15[k]
                  + f_5 * dd_14[k]
                  + pb_y[k] * df_26[k];

        t_16[k] = f_0 * pf_3[k]
                  - f_1 * dd_s_14[k]
                  + f_2 * dg_s_16[k]
                  + f_3 * dd_14[k]
                  + pb_z[k] * df_27[k];
    }
}

auto
compute_prim_dg_kinetic_energy_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf, const size_t pg,
                                 const size_t dd_s, const size_t dg_s, const size_t dd,
                                 const size_t df, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;

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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_6 = buffer.data(pg + 6);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_5 = buffer.data(dd_s + 5);
    const auto *dd_s_6 = buffer.data(dd_s + 6);
    const auto *dd_s_7 = buffer.data(dd_s + 7);
    const auto *dd_s_9 = buffer.data(dd_s + 9);
    const auto *dd_s_10 = buffer.data(dd_s + 10);
    const auto *dd_s_11 = buffer.data(dd_s + 11);
    const auto *dd_s_12 = buffer.data(dd_s + 12);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_1 = buffer.data(dg_s + 1);
    const auto *dg_s_2 = buffer.data(dg_s + 2);
    const auto *dg_s_3 = buffer.data(dg_s + 3);
    const auto *dg_s_4 = buffer.data(dg_s + 4);
    const auto *dg_s_5 = buffer.data(dg_s + 5);
    const auto *dg_s_6 = buffer.data(dg_s + 6);
    const auto *dg_s_7 = buffer.data(dg_s + 7);
    const auto *dg_s_8 = buffer.data(dg_s + 8);
    const auto *dg_s_9 = buffer.data(dg_s + 9);
    const auto *dg_s_10 = buffer.data(dg_s + 10);
    const auto *dg_s_11 = buffer.data(dg_s + 11);
    const auto *dg_s_12 = buffer.data(dg_s + 12);
    const auto *dg_s_13 = buffer.data(dg_s + 13);
    const auto *dg_s_14 = buffer.data(dg_s + 14);
    const auto *dg_s_15 = buffer.data(dg_s + 15);
    const auto *dg_s_16 = buffer.data(dg_s + 16);
    const auto *dg_s_17 = buffer.data(dg_s + 17);
    const auto *dg_s_18 = buffer.data(dg_s + 18);
    const auto *dg_s_19 = buffer.data(dg_s + 19);
    const auto *dg_s_20 = buffer.data(dg_s + 20);
    const auto *dg_s_21 = buffer.data(dg_s + 21);
    const auto *dg_s_22 = buffer.data(dg_s + 22);
    const auto *dg_s_23 = buffer.data(dg_s + 23);
    const auto *dg_s_24 = buffer.data(dg_s + 24);
    const auto *dg_s_25 = buffer.data(dg_s + 25);
    const auto *dg_s_26 = buffer.data(dg_s + 26);
    const auto *dg_s_27 = buffer.data(dg_s + 27);
    const auto *dg_s_28 = buffer.data(dg_s + 28);
    const auto *dg_s_29 = buffer.data(dg_s + 29);
    const auto *dg_s_30 = buffer.data(dg_s + 30);
    const auto *dg_s_31 = buffer.data(dg_s + 31);
    const auto *dg_s_32 = buffer.data(dg_s + 32);
    const auto *dg_s_33 = buffer.data(dg_s + 33);
    const auto *dg_s_34 = buffer.data(dg_s + 34);
    const auto *dg_s_35 = buffer.data(dg_s + 35);
    const auto *dg_s_36 = buffer.data(dg_s + 36);
    const auto *dg_s_37 = buffer.data(dg_s + 37);
    const auto *dg_s_38 = buffer.data(dg_s + 38);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pf_0, dd_s_0, dg_s_0, dg_s_1, \
                         dg_s_2, dg_s_3, dd_0, df_0, df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 - f_1 * dd_s_0[k]
                 + f_2 * dg_s_0[k]
                 + f_3 * dd_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_2 * dg_s_1[k]
                 + pb_y[k] * df_0[k];

        t_2[k] = f_2 * dg_s_2[k]
                 + pb_z[k] * df_0[k];

        t_3[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_3[k]
                 + f_5 * dd_0[k]
                 + pb_y[k] * df_1[k];
    }

#pragma omp simd aligned(t_4, t_5, pb_y, pb_z, dd_s_0, dd_s_1, dg_s_4, dg_s_5, dd_0, dd_1, \
                         df_2, df_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_4[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_2[k];

        t_5[k] = -f_1 * dd_s_1[k]
                 + f_2 * dg_s_5[k]
                 + f_3 * dd_1[k]
                 + pb_y[k] * df_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, pb_z, pf_1, pg_0, pg_1, dd_s_2, dg_s_6, \
                         dg_s_7, dg_s_8, dd_2, df_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_1 * dd_s_2[k]
                 + f_2 * dg_s_6[k]
                 + f_3 * dd_2[k]
                 + pb_z[k] * df_4[k];

        t_7[k] = pa_y[k] * pg_0[k]
                 + f_2 * dg_s_7[k];

        t_8[k] = f_0 * pf_1[k]
                 + pa_x[k] * pg_1[k]
                 + f_2 * dg_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_z, pb_z, pf_0, pg_0, pg_2, dg_s_9, dg_s_10, \
                         dg_s_11, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_x[k] * pg_2[k]
                 + f_2 * dg_s_9[k];

        t_10[k] = pa_z[k] * pg_0[k]
                  + f_2 * dg_s_10[k];

        t_11[k] = f_5 * pf_0[k]
                  + f_2 * dg_s_11[k]
                  + pb_z[k] * df_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_x, pf_4, pg_3, pg_6, dd_s_5, dg_s_12, \
                         dg_s_13, dg_s_14, dd_5, df_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * pf_4[k]
                  + pa_x[k] * pg_3[k]
                  + f_2 * dg_s_12[k];

        t_13[k] = pa_x[k] * pg_6[k]
                  + f_2 * dg_s_13[k];

        t_14[k] = -f_1 * dd_s_5[k]
                  + f_2 * dg_s_14[k]
                  + f_3 * dd_5[k]
                  + pb_x[k] * df_8[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_x, dd_s_6, dd_s_7, dg_s_15, dg_s_16, dg_s_17, \
                         dd_6, dd_7, df_9, df_10, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_4 * dd_s_6[k]
                  + f_2 * dg_s_15[k]
                  + f_5 * dd_6[k]
                  + pb_x[k] * df_9[k];

        t_16[k] = -f_4 * dd_s_7[k]
                  + f_2 * dg_s_16[k]
                  + f_5 * dd_7[k]
                  + pb_x[k] * df_10[k];

        t_17[k] = f_2 * dg_s_17[k]
                  + pb_x[k] * df_11[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_y, pb_z, pf_2, dd_s_6, dg_s_18, dg_s_19, \
                         dg_s_20, dd_6, df_11, df_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * pf_2[k]
                  - f_1 * dd_s_6[k]
                  + f_2 * dg_s_18[k]
                  + f_3 * dd_6[k]
                  + pb_y[k] * df_11[k];

        t_19[k] = f_2 * dg_s_19[k]
                  + pb_z[k] * df_11[k];

        t_20[k] = -f_4 * dd_s_6[k]
                  + f_2 * dg_s_20[k]
                  + f_5 * dd_6[k]
                  + pb_z[k] * df_12[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_z, pb_y, pb_z, pf_3, pg_2, dd_s_7, dg_s_21, \
                         dg_s_22, dg_s_23, dd_7, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * pf_3[k]
                  + f_2 * dg_s_21[k]
                  + pb_y[k] * df_13[k];

        t_22[k] = -f_1 * dd_s_7[k]
                  + f_2 * dg_s_22[k]
                  + f_3 * dd_7[k]
                  + pb_z[k] * df_13[k];

        t_23[k] = pa_z[k] * pg_2[k]
                  + f_2 * dg_s_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_y, pb_y, pb_z, pf_2, pf_5, pf_6, pg_4, dg_s_24, \
                         dg_s_25, dg_s_26, df_14, df_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * pf_2[k]
                  + f_2 * dg_s_24[k]
                  + pb_z[k] * df_14[k];

        t_25[k] = f_0 * pf_5[k]
                  + pa_y[k] * pg_4[k]
                  + f_2 * dg_s_25[k];

        t_26[k] = f_5 * pf_6[k]
                  + f_2 * dg_s_26[k]
                  + pb_y[k] * df_15[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_x, pb_y, pg_6, dd_s_9, dg_s_27, dg_s_28, \
                         dg_s_29, dd_9, df_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_y[k] * pg_6[k]
                  + f_2 * dg_s_27[k];

        t_28[k] = -f_1 * dd_s_9[k]
                  + f_2 * dg_s_28[k]
                  + f_3 * dd_9[k]
                  + pb_x[k] * df_16[k];

        t_29[k] = f_2 * dg_s_29[k]
                  + pb_y[k] * df_16[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_x, dd_s_10, dd_s_12, dg_s_30, dg_s_31, dg_s_32, \
                         dd_10, dd_12, df_18, df_19, df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -f_4 * dd_s_10[k]
                  + f_2 * dg_s_30[k]
                  + f_5 * dd_10[k]
                  + pb_x[k] * df_18[k];

        t_31[k] = -f_4 * dd_s_12[k]
                  + f_2 * dg_s_31[k]
                  + f_5 * dd_12[k]
                  + pb_x[k] * df_19[k];

        t_32[k] = f_2 * dg_s_32[k]
                  + pb_x[k] * df_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_x, pb_y, dd_s_10, dd_s_11, dg_s_33, dg_s_34, \
                         dg_s_35, dd_10, dd_11, df_20, df_21, df_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_2 * dg_s_33[k]
                  + pb_x[k] * df_23[k];

        t_34[k] = -f_1 * dd_s_10[k]
                  + f_2 * dg_s_34[k]
                  + f_3 * dd_10[k]
                  + pb_y[k] * df_20[k];

        t_35[k] = -f_6 * dd_s_11[k]
                  + f_2 * dg_s_35[k]
                  + f_0 * dd_11[k]
                  + pb_y[k] * df_21[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_y, pb_z, pf_6, dd_s_12, dg_s_36, dg_s_37, \
                         dg_s_38, dd_12, df_22, df_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_4 * dd_s_12[k]
                  + f_2 * dg_s_36[k]
                  + f_5 * dd_12[k]
                  + pb_y[k] * df_22[k];

        t_37[k] = f_2 * dg_s_37[k]
                  + pb_y[k] * df_23[k];

        t_38[k] = f_0 * pf_6[k]
                  - f_1 * dd_s_12[k]
                  + f_2 * dg_s_38[k]
                  + f_3 * dd_12[k]
                  + pb_z[k] * df_23[k];
    }
}

auto
compute_prim_dg_kinetic_energy_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf, const size_t pg,
                                 const size_t dd_s, const size_t dg_s, const size_t dd,
                                 const size_t df, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;

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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_9 = buffer.data(pg + 9);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_2 = buffer.data(dd_s + 2);
    const auto *dd_s_5 = buffer.data(dd_s + 5);
    const auto *dd_s_6 = buffer.data(dd_s + 6);
    const auto *dd_s_7 = buffer.data(dd_s + 7);
    const auto *dd_s_9 = buffer.data(dd_s + 9);
    const auto *dd_s_10 = buffer.data(dd_s + 10);
    const auto *dd_s_11 = buffer.data(dd_s + 11);
    const auto *dd_s_12 = buffer.data(dd_s + 12);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_1 = buffer.data(dg_s + 1);
    const auto *dg_s_2 = buffer.data(dg_s + 2);
    const auto *dg_s_3 = buffer.data(dg_s + 3);
    const auto *dg_s_4 = buffer.data(dg_s + 4);
    const auto *dg_s_5 = buffer.data(dg_s + 5);
    const auto *dg_s_6 = buffer.data(dg_s + 6);
    const auto *dg_s_7 = buffer.data(dg_s + 7);
    const auto *dg_s_8 = buffer.data(dg_s + 8);
    const auto *dg_s_9 = buffer.data(dg_s + 9);
    const auto *dg_s_10 = buffer.data(dg_s + 10);
    const auto *dg_s_11 = buffer.data(dg_s + 11);
    const auto *dg_s_12 = buffer.data(dg_s + 12);
    const auto *dg_s_13 = buffer.data(dg_s + 13);
    const auto *dg_s_14 = buffer.data(dg_s + 14);
    const auto *dg_s_15 = buffer.data(dg_s + 15);
    const auto *dg_s_16 = buffer.data(dg_s + 16);
    const auto *dg_s_17 = buffer.data(dg_s + 17);
    const auto *dg_s_18 = buffer.data(dg_s + 18);
    const auto *dg_s_19 = buffer.data(dg_s + 19);
    const auto *dg_s_20 = buffer.data(dg_s + 20);
    const auto *dg_s_21 = buffer.data(dg_s + 21);
    const auto *dg_s_22 = buffer.data(dg_s + 22);
    const auto *dg_s_23 = buffer.data(dg_s + 23);
    const auto *dg_s_24 = buffer.data(dg_s + 24);
    const auto *dg_s_25 = buffer.data(dg_s + 25);
    const auto *dg_s_26 = buffer.data(dg_s + 26);
    const auto *dg_s_27 = buffer.data(dg_s + 27);
    const auto *dg_s_28 = buffer.data(dg_s + 28);
    const auto *dg_s_29 = buffer.data(dg_s + 29);
    const auto *dg_s_30 = buffer.data(dg_s + 30);
    const auto *dg_s_31 = buffer.data(dg_s + 31);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pf_0, dd_s_0, dg_s_0, dg_s_1, \
                         dg_s_2, dg_s_3, dd_0, df_0, df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 - f_1 * dd_s_0[k]
                 + f_2 * dg_s_0[k]
                 + f_3 * dd_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_2 * dg_s_1[k]
                 + pb_y[k] * df_0[k];

        t_2[k] = f_2 * dg_s_2[k]
                 + pb_z[k] * df_0[k];

        t_3[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_3[k]
                 + f_5 * dd_0[k]
                 + pb_y[k] * df_1[k];
    }

#pragma omp simd aligned(t_4, t_5, pb_y, pb_z, dd_s_0, dd_s_1, dg_s_4, dg_s_5, dd_0, dd_1, \
                         df_2, df_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_4[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_2[k];

        t_5[k] = -f_1 * dd_s_1[k]
                 + f_2 * dg_s_5[k]
                 + f_3 * dd_1[k]
                 + pb_y[k] * df_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_z, pg_0, pg_3, dd_s_2, dg_s_6, dg_s_7, \
                         dg_s_8, dd_2, df_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_1 * dd_s_2[k]
                 + f_2 * dg_s_6[k]
                 + f_3 * dd_2[k]
                 + pb_z[k] * df_4[k];

        t_7[k] = pa_x[k] * pg_3[k]
                 + f_2 * dg_s_7[k];

        t_8[k] = pa_z[k] * pg_0[k]
                 + f_2 * dg_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_x, pg_9, dd_s_5, dd_s_6, dg_s_9, dg_s_10, \
                         dg_s_11, dd_5, dd_6, df_7, df_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_x[k] * pg_9[k]
                 + f_2 * dg_s_9[k];

        t_10[k] = -f_1 * dd_s_5[k]
                  + f_2 * dg_s_10[k]
                  + f_3 * dd_5[k]
                  + pb_x[k] * df_7[k];

        t_11[k] = -f_4 * dd_s_6[k]
                  + f_2 * dg_s_11[k]
                  + f_5 * dd_6[k]
                  + pb_x[k] * df_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_x, dd_s_7, dg_s_12, dg_s_13, dg_s_14, dd_7, \
                         df_9, df_10, df_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_4 * dd_s_7[k]
                  + f_2 * dg_s_12[k]
                  + f_5 * dd_7[k]
                  + pb_x[k] * df_9[k];

        t_13[k] = f_2 * dg_s_13[k]
                  + pb_x[k] * df_10[k];

        t_14[k] = f_2 * dg_s_14[k]
                  + pb_x[k] * df_12[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_y, pb_z, pf_1, dd_s_6, dg_s_15, dg_s_16, \
                         dg_s_17, dd_6, df_10, df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * pf_1[k]
                  - f_1 * dd_s_6[k]
                  + f_2 * dg_s_15[k]
                  + f_3 * dd_6[k]
                  + pb_y[k] * df_10[k];

        t_16[k] = f_2 * dg_s_16[k]
                  + pb_z[k] * df_10[k];

        t_17[k] = -f_4 * dd_s_6[k]
                  + f_2 * dg_s_17[k]
                  + f_5 * dd_6[k]
                  + pb_z[k] * df_11[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pa_z, pb_z, pg_3, pg_9, dd_s_7, dg_s_18, \
                         dg_s_19, dg_s_20, dd_7, df_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -f_1 * dd_s_7[k]
                  + f_2 * dg_s_18[k]
                  + f_3 * dd_7[k]
                  + pb_z[k] * df_12[k];

        t_19[k] = pa_z[k] * pg_3[k]
                  + f_2 * dg_s_19[k];

        t_20[k] = pa_y[k] * pg_9[k]
                  + f_2 * dg_s_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_x, pb_y, dd_s_9, dd_s_10, dg_s_21, dg_s_22, \
                         dg_s_23, dd_9, dd_10, df_14, df_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_1 * dd_s_9[k]
                  + f_2 * dg_s_21[k]
                  + f_3 * dd_9[k]
                  + pb_x[k] * df_14[k];

        t_22[k] = f_2 * dg_s_22[k]
                  + pb_y[k] * df_14[k];

        t_23[k] = -f_4 * dd_s_10[k]
                  + f_2 * dg_s_23[k]
                  + f_5 * dd_10[k]
                  + pb_x[k] * df_16[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_x, dd_s_12, dg_s_24, dg_s_25, dg_s_26, dd_12, \
                         df_17, df_18, df_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_4 * dd_s_12[k]
                  + f_2 * dg_s_24[k]
                  + f_5 * dd_12[k]
                  + pb_x[k] * df_17[k];

        t_25[k] = f_2 * dg_s_25[k]
                  + pb_x[k] * df_18[k];

        t_26[k] = f_2 * dg_s_26[k]
                  + pb_x[k] * df_21[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, dd_s_10, dd_s_11, dd_s_12, dg_s_27, dg_s_28, \
                         dg_s_29, dd_10, dd_11, dd_12, df_18, df_19, \
                         df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_1 * dd_s_10[k]
                  + f_2 * dg_s_27[k]
                  + f_3 * dd_10[k]
                  + pb_y[k] * df_18[k];

        t_28[k] = -f_6 * dd_s_11[k]
                  + f_2 * dg_s_28[k]
                  + f_0 * dd_11[k]
                  + pb_y[k] * df_19[k];

        t_29[k] = -f_4 * dd_s_12[k]
                  + f_2 * dg_s_29[k]
                  + f_5 * dd_12[k]
                  + pb_y[k] * df_20[k];
    }

#pragma omp simd aligned(t_30, t_31, pb_y, pb_z, pf_2, dd_s_12, dg_s_30, dg_s_31, dd_12, \
                         df_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * dg_s_30[k]
                  + pb_y[k] * df_21[k];

        t_31[k] = f_0 * pf_2[k]
                  - f_1 * dd_s_12[k]
                  + f_2 * dg_s_31[k]
                  + f_3 * dd_12[k]
                  + pb_z[k] * df_21[k];
    }
}

auto
compute_prim_dg_kinetic_energy_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf, const size_t pg,
                                 const size_t dd_s, const size_t dg_s, const size_t dd,
                                 const size_t df, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_6 = buffer.data(dd_s + 6);
    const auto *dd_s_11 = buffer.data(dd_s + 11);
    const auto *dd_s_12 = buffer.data(dd_s + 12);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_1 = buffer.data(dg_s + 1);
    const auto *dg_s_2 = buffer.data(dg_s + 2);
    const auto *dg_s_3 = buffer.data(dg_s + 3);
    const auto *dg_s_4 = buffer.data(dg_s + 4);
    const auto *dg_s_5 = buffer.data(dg_s + 5);
    const auto *dg_s_6 = buffer.data(dg_s + 6);
    const auto *dg_s_7 = buffer.data(dg_s + 7);
    const auto *dg_s_8 = buffer.data(dg_s + 8);
    const auto *dg_s_9 = buffer.data(dg_s + 9);
    const auto *dg_s_10 = buffer.data(dg_s + 10);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pb_x, pb_z, pf_0, pg_0, dd_s_0, dg_s_0, dg_s_1, \
                         dg_s_2, dd_0, df_0, df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 - f_1 * dd_s_0[k]
                 + f_2 * dg_s_0[k]
                 + f_3 * dd_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_1[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_1[k];

        t_2[k] = pa_x[k] * pg_0[k]
                 + f_2 * dg_s_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_y, pb_z, pf_1, pg_2, dd_s_6, dg_s_3, dg_s_4, \
                         dg_s_5, dd_6, df_5, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_x[k] * pg_2[k]
                 + f_2 * dg_s_3[k];

        t_4[k] = f_0 * pf_1[k]
                 - f_1 * dd_s_6[k]
                 + f_2 * dg_s_4[k]
                 + f_3 * dd_6[k]
                 + pb_y[k] * df_5[k];

        t_5[k] = -f_4 * dd_s_6[k]
                 + f_2 * dg_s_5[k]
                 + f_5 * dd_6[k]
                 + pb_z[k] * df_6[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_y, pb_y, pf_2, pg_1, pg_2, dd_s_11, dg_s_6, dg_s_7, \
                         dg_s_8, dd_11, df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * pf_2[k]
                 + pa_y[k] * pg_1[k]
                 + f_2 * dg_s_6[k];

        t_7[k] = pa_y[k] * pg_2[k]
                 + f_2 * dg_s_7[k];

        t_8[k] = -f_1 * dd_s_11[k]
                 + f_2 * dg_s_8[k]
                 + f_3 * dd_11[k]
                 + pb_y[k] * df_9[k];
    }

#pragma omp simd aligned(t_9, t_10, pb_y, pb_z, pf_3, dd_s_12, dg_s_9, dg_s_10, dd_12, df_10, \
                         df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_4 * dd_s_12[k]
                 + f_2 * dg_s_9[k]
                 + f_5 * dd_12[k]
                 + pb_y[k] * df_10[k];

        t_10[k] = f_0 * pf_3[k]
                  - f_1 * dd_s_12[k]
                  + f_2 * dg_s_10[k]
                  + f_3 * dd_12[k]
                  + pb_z[k] * df_11[k];
    }
}

auto
compute_prim_dg_kinetic_energy_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pf, const size_t pg,
                                  const size_t dd_s, const size_t dg_s, const size_t dd,
                                  const size_t df, const size_t ncols, const double alpha,
                                  const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);

    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_4 = buffer.data(pg + 4);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_5 = buffer.data(dd_s + 5);
    const auto *dd_s_8 = buffer.data(dd_s + 8);
    const auto *dd_s_9 = buffer.data(dd_s + 9);
    const auto *dd_s_10 = buffer.data(dd_s + 10);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_1 = buffer.data(dg_s + 1);
    const auto *dg_s_2 = buffer.data(dg_s + 2);
    const auto *dg_s_3 = buffer.data(dg_s + 3);
    const auto *dg_s_4 = buffer.data(dg_s + 4);
    const auto *dg_s_5 = buffer.data(dg_s + 5);
    const auto *dg_s_6 = buffer.data(dg_s + 6);
    const auto *dg_s_7 = buffer.data(dg_s + 7);
    const auto *dg_s_8 = buffer.data(dg_s + 8);
    const auto *dg_s_9 = buffer.data(dg_s + 9);
    const auto *dg_s_10 = buffer.data(dg_s + 10);
    const auto *dg_s_11 = buffer.data(dg_s + 11);
    const auto *dg_s_12 = buffer.data(dg_s + 12);
    const auto *dg_s_13 = buffer.data(dg_s + 13);
    const auto *dg_s_14 = buffer.data(dg_s + 14);
    const auto *dg_s_15 = buffer.data(dg_s + 15);
    const auto *dg_s_16 = buffer.data(dg_s + 16);
    const auto *dg_s_17 = buffer.data(dg_s + 17);
    const auto *dg_s_18 = buffer.data(dg_s + 18);
    const auto *dg_s_19 = buffer.data(dg_s + 19);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pf_0, dd_s_0, dg_s_0, dg_s_1, \
                         dg_s_2, dd_0, df_0, df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 - f_1 * dd_s_0[k]
                 + f_2 * dg_s_0[k]
                 + f_3 * dd_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_2 * dg_s_1[k]
                 + pb_z[k] * df_0[k];

        t_2[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_2[k]
                 + f_5 * dd_0[k]
                 + pb_y[k] * df_1[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_z, pg_1, pg_4, dd_s_0, dg_s_3, dg_s_4, \
                         dg_s_5, dd_0, df_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_3[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_2[k];

        t_4[k] = pa_x[k] * pg_1[k]
                 + f_2 * dg_s_4[k];

        t_5[k] = pa_x[k] * pg_4[k]
                 + f_2 * dg_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, pf_1, dd_s_4, dd_s_5, dg_s_6, dg_s_7, \
                         dg_s_8, dd_4, dd_5, df_7, df_8, df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_1 * dd_s_4[k]
                 + f_2 * dg_s_6[k]
                 + f_3 * dd_4[k]
                 + pb_x[k] * df_7[k];

        t_7[k] = -f_4 * dd_s_5[k]
                 + f_2 * dg_s_7[k]
                 + f_5 * dd_5[k]
                 + pb_x[k] * df_8[k];

        t_8[k] = f_0 * pf_1[k]
                 - f_1 * dd_s_5[k]
                 + f_2 * dg_s_8[k]
                 + f_3 * dd_5[k]
                 + pb_y[k] * df_9[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_z, pf_2, pg_2, dd_s_5, dg_s_9, dg_s_10, \
                         dg_s_11, dd_5, df_9, df_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * dg_s_9[k]
                 + pb_z[k] * df_9[k];

        t_10[k] = -f_4 * dd_s_5[k]
                  + f_2 * dg_s_10[k]
                  + f_5 * dd_5[k]
                  + pb_z[k] * df_10[k];

        t_11[k] = f_0 * pf_2[k]
                  + pa_y[k] * pg_2[k]
                  + f_2 * dg_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pb_x, pg_4, dd_s_8, dd_s_9, dg_s_12, dg_s_13, \
                         dg_s_14, dd_8, dd_9, df_14, df_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_y[k] * pg_4[k]
                  + f_2 * dg_s_12[k];

        t_13[k] = -f_1 * dd_s_8[k]
                  + f_2 * dg_s_13[k]
                  + f_3 * dd_8[k]
                  + pb_x[k] * df_14[k];

        t_14[k] = -f_4 * dd_s_9[k]
                  + f_2 * dg_s_14[k]
                  + f_5 * dd_9[k]
                  + pb_x[k] * df_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_x, pb_y, dd_s_9, dd_s_10, dg_s_15, dg_s_16, \
                         dg_s_17, dd_9, dd_10, df_17, df_18, df_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_4 * dd_s_10[k]
                  + f_2 * dg_s_15[k]
                  + f_5 * dd_10[k]
                  + pb_x[k] * df_17[k];

        t_16[k] = -f_1 * dd_s_9[k]
                  + f_2 * dg_s_16[k]
                  + f_3 * dd_9[k]
                  + pb_y[k] * df_18[k];

        t_17[k] = -f_4 * dd_s_10[k]
                  + f_2 * dg_s_17[k]
                  + f_5 * dd_10[k]
                  + pb_y[k] * df_19[k];
    }

#pragma omp simd aligned(t_18, t_19, pb_y, pb_z, pf_3, dd_s_10, dg_s_18, dg_s_19, dd_10, \
                         df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * dg_s_18[k]
                  + pb_y[k] * df_20[k];

        t_19[k] = f_0 * pf_3[k]
                  - f_1 * dd_s_10[k]
                  + f_2 * dg_s_19[k]
                  + f_3 * dd_10[k]
                  + pb_z[k] * df_20[k];
    }
}

auto
compute_prim_dg_kinetic_energy_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pf, const size_t pg,
                                  const size_t dd_s, const size_t dg_s, const size_t dd,
                                  const size_t df, const size_t ncols, const double alpha,
                                  const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
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
    const auto *pg_4 = buffer.data(pg + 4);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_1 = buffer.data(dd_s + 1);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_5 = buffer.data(dd_s + 5);
    const auto *dd_s_6 = buffer.data(dd_s + 6);
    const auto *dd_s_8 = buffer.data(dd_s + 8);
    const auto *dd_s_9 = buffer.data(dd_s + 9);
    const auto *dd_s_10 = buffer.data(dd_s + 10);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_1 = buffer.data(dg_s + 1);
    const auto *dg_s_2 = buffer.data(dg_s + 2);
    const auto *dg_s_3 = buffer.data(dg_s + 3);
    const auto *dg_s_4 = buffer.data(dg_s + 4);
    const auto *dg_s_5 = buffer.data(dg_s + 5);
    const auto *dg_s_6 = buffer.data(dg_s + 6);
    const auto *dg_s_7 = buffer.data(dg_s + 7);
    const auto *dg_s_8 = buffer.data(dg_s + 8);
    const auto *dg_s_9 = buffer.data(dg_s + 9);
    const auto *dg_s_10 = buffer.data(dg_s + 10);
    const auto *dg_s_11 = buffer.data(dg_s + 11);
    const auto *dg_s_12 = buffer.data(dg_s + 12);
    const auto *dg_s_13 = buffer.data(dg_s + 13);
    const auto *dg_s_14 = buffer.data(dg_s + 14);
    const auto *dg_s_15 = buffer.data(dg_s + 15);
    const auto *dg_s_16 = buffer.data(dg_s + 16);
    const auto *dg_s_17 = buffer.data(dg_s + 17);
    const auto *dg_s_18 = buffer.data(dg_s + 18);
    const auto *dg_s_19 = buffer.data(dg_s + 19);
    const auto *dg_s_20 = buffer.data(dg_s + 20);
    const auto *dg_s_21 = buffer.data(dg_s + 21);
    const auto *dg_s_22 = buffer.data(dg_s + 22);
    const auto *dg_s_23 = buffer.data(dg_s + 23);
    const auto *dg_s_24 = buffer.data(dg_s + 24);
    const auto *dg_s_25 = buffer.data(dg_s + 25);
    const auto *dg_s_26 = buffer.data(dg_s + 26);
    const auto *dg_s_27 = buffer.data(dg_s + 27);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pf_0, dd_s_0, dg_s_0, dg_s_1, \
                         dg_s_2, dg_s_3, dd_0, df_0, df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 - f_1 * dd_s_0[k]
                 + f_2 * dg_s_0[k]
                 + f_3 * dd_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_2 * dg_s_1[k]
                 + pb_y[k] * df_0[k];

        t_2[k] = f_2 * dg_s_2[k]
                 + pb_z[k] * df_0[k];

        t_3[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_3[k]
                 + f_5 * dd_0[k]
                 + pb_y[k] * df_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pb_z, pg_1, dd_s_0, dd_s_1, dg_s_4, dg_s_5, \
                         dg_s_6, dd_0, dd_1, df_2, df_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_4[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_2[k];

        t_5[k] = -f_1 * dd_s_1[k]
                 + f_2 * dg_s_5[k]
                 + f_3 * dd_1[k]
                 + pb_z[k] * df_3[k];

        t_6[k] = pa_x[k] * pg_1[k]
                 + f_2 * dg_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_z, pb_x, pg_0, pg_4, dd_s_4, dg_s_7, dg_s_8, \
                         dg_s_9, dd_4, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_z[k] * pg_0[k]
                 + f_2 * dg_s_7[k];

        t_8[k] = pa_x[k] * pg_4[k]
                 + f_2 * dg_s_8[k];

        t_9[k] = -f_1 * dd_s_4[k]
                 + f_2 * dg_s_9[k]
                 + f_3 * dd_4[k]
                 + pb_x[k] * df_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_x, pb_y, pb_z, pf_1, dd_s_5, dg_s_10, \
                         dg_s_11, dg_s_12, dg_s_13, dd_5, df_7, df_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_4 * dd_s_5[k]
                  + f_2 * dg_s_10[k]
                  + f_5 * dd_5[k]
                  + pb_x[k] * df_7[k];

        t_11[k] = f_2 * dg_s_11[k]
                  + pb_x[k] * df_8[k];

        t_12[k] = f_0 * pf_1[k]
                  - f_1 * dd_s_5[k]
                  + f_2 * dg_s_12[k]
                  + f_3 * dd_5[k]
                  + pb_y[k] * df_8[k];

        t_13[k] = f_2 * dg_s_13[k]
                  + pb_z[k] * df_8[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_z, pb_z, pg_1, dd_s_5, dd_s_6, dg_s_14, dg_s_15, \
                         dg_s_16, dd_5, dd_6, df_9, df_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -f_4 * dd_s_5[k]
                  + f_2 * dg_s_14[k]
                  + f_5 * dd_5[k]
                  + pb_z[k] * df_9[k];

        t_15[k] = -f_1 * dd_s_6[k]
                  + f_2 * dg_s_15[k]
                  + f_3 * dd_6[k]
                  + pb_z[k] * df_10[k];

        t_16[k] = pa_z[k] * pg_1[k]
                  + f_2 * dg_s_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_y, pb_x, pb_y, pg_4, dd_s_8, dg_s_17, dg_s_18, \
                         dg_s_19, dd_8, df_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = pa_y[k] * pg_4[k]
                  + f_2 * dg_s_17[k];

        t_18[k] = -f_1 * dd_s_8[k]
                  + f_2 * dg_s_18[k]
                  + f_3 * dd_8[k]
                  + pb_x[k] * df_12[k];

        t_19[k] = f_2 * dg_s_19[k]
                  + pb_y[k] * df_12[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_x, dd_s_9, dd_s_10, dg_s_20, dg_s_21, dg_s_22, \
                         dd_9, dd_10, df_14, df_15, df_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_4 * dd_s_9[k]
                  + f_2 * dg_s_20[k]
                  + f_5 * dd_9[k]
                  + pb_x[k] * df_14[k];

        t_21[k] = -f_4 * dd_s_10[k]
                  + f_2 * dg_s_21[k]
                  + f_5 * dd_10[k]
                  + pb_x[k] * df_15[k];

        t_22[k] = f_2 * dg_s_22[k]
                  + pb_x[k] * df_16[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_x, pb_y, dd_s_9, dd_s_10, dg_s_23, dg_s_24, \
                         dg_s_25, dd_9, dd_10, df_16, df_17, df_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_2 * dg_s_23[k]
                  + pb_x[k] * df_18[k];

        t_24[k] = -f_1 * dd_s_9[k]
                  + f_2 * dg_s_24[k]
                  + f_3 * dd_9[k]
                  + pb_y[k] * df_16[k];

        t_25[k] = -f_4 * dd_s_10[k]
                  + f_2 * dg_s_25[k]
                  + f_5 * dd_10[k]
                  + pb_y[k] * df_17[k];
    }

#pragma omp simd aligned(t_26, t_27, pb_y, pb_z, pf_2, dd_s_10, dg_s_26, dg_s_27, dd_10, \
                         df_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * dg_s_26[k]
                  + pb_y[k] * df_18[k];

        t_27[k] = f_0 * pf_2[k]
                  - f_1 * dd_s_10[k]
                  + f_2 * dg_s_27[k]
                  + f_3 * dd_10[k]
                  + pb_z[k] * df_18[k];
    }
}

auto
compute_prim_dg_kinetic_energy_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pf, const size_t pg,
                                  const size_t dd_s, const size_t dg_s, const size_t dd,
                                  const size_t df, const size_t ncols, const double alpha,
                                  const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_7 = buffer.data(dd_s + 7);
    const auto *dd_s_8 = buffer.data(dd_s + 8);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_1 = buffer.data(dg_s + 1);
    const auto *dg_s_2 = buffer.data(dg_s + 2);
    const auto *dg_s_3 = buffer.data(dg_s + 3);
    const auto *dg_s_4 = buffer.data(dg_s + 4);
    const auto *dg_s_5 = buffer.data(dg_s + 5);
    const auto *dg_s_6 = buffer.data(dg_s + 6);
    const auto *dg_s_7 = buffer.data(dg_s + 7);
    const auto *dg_s_8 = buffer.data(dg_s + 8);
    const auto *dg_s_9 = buffer.data(dg_s + 9);
    const auto *dg_s_10 = buffer.data(dg_s + 10);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pb_x, pb_z, pf_0, pg_0, dd_s_0, dg_s_0, dg_s_1, \
                         dg_s_2, dd_0, df_0, df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 - f_1 * dd_s_0[k]
                 + f_2 * dg_s_0[k]
                 + f_3 * dd_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_1[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_1[k];

        t_2[k] = pa_x[k] * pg_0[k]
                 + f_2 * dg_s_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_y, pb_z, pf_1, pg_2, dd_s_4, dg_s_3, dg_s_4, \
                         dg_s_5, dd_4, df_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_x[k] * pg_2[k]
                 + f_2 * dg_s_3[k];

        t_4[k] = f_0 * pf_1[k]
                 - f_1 * dd_s_4[k]
                 + f_2 * dg_s_4[k]
                 + f_3 * dd_4[k]
                 + pb_y[k] * df_4[k];

        t_5[k] = -f_4 * dd_s_4[k]
                 + f_2 * dg_s_5[k]
                 + f_5 * dd_4[k]
                 + pb_z[k] * df_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_y, pb_y, pf_2, pg_1, pg_2, dd_s_7, dg_s_6, dg_s_7, \
                         dg_s_8, dd_7, df_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * pf_2[k]
                 + pa_y[k] * pg_1[k]
                 + f_2 * dg_s_6[k];

        t_7[k] = pa_y[k] * pg_2[k]
                 + f_2 * dg_s_7[k];

        t_8[k] = -f_1 * dd_s_7[k]
                 + f_2 * dg_s_8[k]
                 + f_3 * dd_7[k]
                 + pb_y[k] * df_7[k];
    }

#pragma omp simd aligned(t_9, t_10, pb_y, pb_z, pf_3, dd_s_8, dg_s_9, dg_s_10, dd_8, df_8, \
                         df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_4 * dd_s_8[k]
                 + f_2 * dg_s_9[k]
                 + f_5 * dd_8[k]
                 + pb_y[k] * df_8[k];

        t_10[k] = f_0 * pf_3[k]
                  - f_1 * dd_s_8[k]
                  + f_2 * dg_s_10[k]
                  + f_3 * dd_8[k]
                  + pb_z[k] * df_9[k];
    }
}

auto
compute_prim_dg_kinetic_energy_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pf, const size_t pg,
                                  const size_t dd_s, const size_t dg_s, const size_t dd,
                                  const size_t df, const size_t ncols, const double alpha,
                                  const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);

    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_4 = buffer.data(pg + 4);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_7 = buffer.data(dd_s + 7);
    const auto *dd_s_8 = buffer.data(dd_s + 8);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_1 = buffer.data(dg_s + 1);
    const auto *dg_s_2 = buffer.data(dg_s + 2);
    const auto *dg_s_3 = buffer.data(dg_s + 3);
    const auto *dg_s_4 = buffer.data(dg_s + 4);
    const auto *dg_s_5 = buffer.data(dg_s + 5);
    const auto *dg_s_6 = buffer.data(dg_s + 6);
    const auto *dg_s_7 = buffer.data(dg_s + 7);
    const auto *dg_s_8 = buffer.data(dg_s + 8);
    const auto *dg_s_9 = buffer.data(dg_s + 9);
    const auto *dg_s_10 = buffer.data(dg_s + 10);
    const auto *dg_s_11 = buffer.data(dg_s + 11);
    const auto *dg_s_12 = buffer.data(dg_s + 12);
    const auto *dg_s_13 = buffer.data(dg_s + 13);
    const auto *dg_s_14 = buffer.data(dg_s + 14);
    const auto *dg_s_15 = buffer.data(dg_s + 15);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_z, pf_0, dd_s_0, dg_s_0, dg_s_1, dg_s_2, \
                         dd_0, df_0, df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 - f_1 * dd_s_0[k]
                 + f_2 * dg_s_0[k]
                 + f_3 * dd_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_2 * dg_s_1[k]
                 + pb_z[k] * df_0[k];

        t_2[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_2[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_1[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_x, pg_1, pg_4, dd_s_4, dg_s_3, dg_s_4, \
                         dg_s_5, dd_4, df_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_x[k] * pg_1[k]
                 + f_2 * dg_s_3[k];

        t_4[k] = pa_x[k] * pg_4[k]
                 + f_2 * dg_s_4[k];

        t_5[k] = -f_4 * dd_s_4[k]
                 + f_2 * dg_s_5[k]
                 + f_5 * dd_4[k]
                 + pb_x[k] * df_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_y, pb_z, pf_1, dd_s_4, dg_s_6, dg_s_7, dg_s_8, \
                         dd_4, df_6, df_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * pf_1[k]
                 - f_1 * dd_s_4[k]
                 + f_2 * dg_s_6[k]
                 + f_3 * dd_4[k]
                 + pb_y[k] * df_6[k];

        t_7[k] = f_2 * dg_s_7[k]
                 + pb_z[k] * df_6[k];

        t_8[k] = -f_4 * dd_s_4[k]
                 + f_2 * dg_s_8[k]
                 + f_5 * dd_4[k]
                 + pb_z[k] * df_7[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pf_2, pg_2, pg_4, dd_s_8, dg_s_9, \
                         dg_s_10, dg_s_11, dd_8, df_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * pf_2[k]
                 + pa_y[k] * pg_2[k]
                 + f_2 * dg_s_9[k];

        t_10[k] = pa_y[k] * pg_4[k]
                  + f_2 * dg_s_10[k];

        t_11[k] = -f_4 * dd_s_8[k]
                  + f_2 * dg_s_11[k]
                  + f_5 * dd_8[k]
                  + pb_x[k] * df_10[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, dd_s_7, dd_s_8, dg_s_12, dg_s_13, dg_s_14, \
                         dd_7, dd_8, df_11, df_12, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_1 * dd_s_7[k]
                  + f_2 * dg_s_12[k]
                  + f_3 * dd_7[k]
                  + pb_y[k] * df_11[k];

        t_13[k] = -f_4 * dd_s_8[k]
                  + f_2 * dg_s_13[k]
                  + f_5 * dd_8[k]
                  + pb_y[k] * df_12[k];

        t_14[k] = f_2 * dg_s_14[k]
                  + pb_y[k] * df_13[k];
    }

#pragma omp simd aligned(t_15, pb_z, pf_3, dd_s_8, dg_s_15, dd_8, \
                         df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * pf_3[k]
                  - f_1 * dd_s_8[k]
                  + f_2 * dg_s_15[k]
                  + f_3 * dd_8[k]
                  + pb_z[k] * df_13[k];
    }
}

auto
compute_prim_dg_kinetic_energy_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pf, const size_t pg,
                                  const size_t dd_s, const size_t dg_s, const size_t dd,
                                  const size_t df, const size_t ncols, const double alpha,
                                  const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_4 = buffer.data(pg + 4);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_3 = buffer.data(dd_s + 3);
    const auto *dd_s_4 = buffer.data(dd_s + 4);
    const auto *dd_s_6 = buffer.data(dd_s + 6);
    const auto *dd_s_7 = buffer.data(dd_s + 7);
    const auto *dd_s_8 = buffer.data(dd_s + 8);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_1 = buffer.data(dg_s + 1);
    const auto *dg_s_2 = buffer.data(dg_s + 2);
    const auto *dg_s_3 = buffer.data(dg_s + 3);
    const auto *dg_s_4 = buffer.data(dg_s + 4);
    const auto *dg_s_5 = buffer.data(dg_s + 5);
    const auto *dg_s_6 = buffer.data(dg_s + 6);
    const auto *dg_s_7 = buffer.data(dg_s + 7);
    const auto *dg_s_8 = buffer.data(dg_s + 8);
    const auto *dg_s_9 = buffer.data(dg_s + 9);
    const auto *dg_s_10 = buffer.data(dg_s + 10);
    const auto *dg_s_11 = buffer.data(dg_s + 11);
    const auto *dg_s_12 = buffer.data(dg_s + 12);
    const auto *dg_s_13 = buffer.data(dg_s + 13);
    const auto *dg_s_14 = buffer.data(dg_s + 14);
    const auto *dg_s_15 = buffer.data(dg_s + 15);
    const auto *dg_s_16 = buffer.data(dg_s + 16);
    const auto *dg_s_17 = buffer.data(dg_s + 17);
    const auto *dg_s_18 = buffer.data(dg_s + 18);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_z, pf_0, dd_s_0, dg_s_0, dg_s_1, dg_s_2, \
                         dd_0, df_0, df_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pf_0[k]
                 - f_1 * dd_s_0[k]
                 + f_2 * dg_s_0[k]
                 + f_3 * dd_0[k]
                 + pb_x[k] * df_0[k];

        t_1[k] = f_2 * dg_s_1[k]
                 + pb_z[k] * df_0[k];

        t_2[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_2[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_1[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_x, pg_1, pg_4, dd_s_3, dg_s_3, dg_s_4, \
                         dg_s_5, dd_3, df_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_x[k] * pg_1[k]
                 + f_2 * dg_s_3[k];

        t_4[k] = pa_x[k] * pg_4[k]
                 + f_2 * dg_s_4[k];

        t_5[k] = -f_1 * dd_s_3[k]
                 + f_2 * dg_s_5[k]
                 + f_3 * dd_3[k]
                 + pb_x[k] * df_4[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, pf_1, dd_s_4, dg_s_6, dg_s_7, \
                         dg_s_8, dg_s_9, dd_4, df_5, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * dd_s_4[k]
                 + f_2 * dg_s_6[k]
                 + f_5 * dd_4[k]
                 + pb_x[k] * df_5[k];

        t_7[k] = f_2 * dg_s_7[k]
                 + pb_x[k] * df_6[k];

        t_8[k] = f_0 * pf_1[k]
                 - f_1 * dd_s_4[k]
                 + f_2 * dg_s_8[k]
                 + f_3 * dd_4[k]
                 + pb_y[k] * df_6[k];

        t_9[k] = f_2 * dg_s_9[k]
                 + pb_z[k] * df_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pb_x, pb_z, pg_4, dd_s_4, dd_s_6, dg_s_10, \
                         dg_s_11, dg_s_12, dd_4, dd_6, df_7, df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_4 * dd_s_4[k]
                  + f_2 * dg_s_10[k]
                  + f_5 * dd_4[k]
                  + pb_z[k] * df_7[k];

        t_11[k] = pa_y[k] * pg_4[k]
                  + f_2 * dg_s_11[k];

        t_12[k] = -f_1 * dd_s_6[k]
                  + f_2 * dg_s_12[k]
                  + f_3 * dd_6[k]
                  + pb_x[k] * df_9[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, pb_y, dd_s_7, dd_s_8, dg_s_13, dg_s_14, \
                         dg_s_15, dd_7, dd_8, df_10, df_11, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_4 * dd_s_8[k]
                  + f_2 * dg_s_13[k]
                  + f_5 * dd_8[k]
                  + pb_x[k] * df_10[k];

        t_14[k] = f_2 * dg_s_14[k]
                  + pb_x[k] * df_13[k];

        t_15[k] = -f_1 * dd_s_7[k]
                  + f_2 * dg_s_15[k]
                  + f_3 * dd_7[k]
                  + pb_y[k] * df_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_y, pb_z, pf_2, dd_s_8, dg_s_16, dg_s_17, \
                         dg_s_18, dd_8, df_12, df_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_4 * dd_s_8[k]
                  + f_2 * dg_s_16[k]
                  + f_5 * dd_8[k]
                  + pb_y[k] * df_12[k];

        t_17[k] = f_2 * dg_s_17[k]
                  + pb_y[k] * df_13[k];

        t_18[k] = f_0 * pf_2[k]
                  - f_1 * dd_s_8[k]
                  + f_2 * dg_s_18[k]
                  + f_3 * dd_8[k]
                  + pb_z[k] * df_13[k];
    }
}

}  // namespace simdkin
