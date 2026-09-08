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
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_16 = buffer.data(pf + 16);
    const auto *pf_18 = buffer.data(pf + 18);
    const auto *pf_19 = buffer.data(pf + 19);
    const auto *pf_22 = buffer.data(pf + 22);
    const auto *pf_25 = buffer.data(pf + 25);
    const auto *pf_27 = buffer.data(pf + 27);
    const auto *pf_28 = buffer.data(pf + 28);
    const auto *pf_29 = buffer.data(pf + 29);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_16 = buffer.data(pg + 16);
    const auto *pg_18 = buffer.data(pg + 18);
    const auto *pg_25 = buffer.data(pg + 25);
    const auto *pg_27 = buffer.data(pg + 27);
    const auto *pg_28 = buffer.data(pg + 28);
    const auto *pg_29 = buffer.data(pg + 29);
    const auto *pg_30 = buffer.data(pg + 30);
    const auto *pg_32 = buffer.data(pg + 32);
    const auto *pg_34 = buffer.data(pg + 34);
    const auto *pg_35 = buffer.data(pg + 35);
    const auto *pg_40 = buffer.data(pg + 40);
    const auto *pg_41 = buffer.data(pg + 41);
    const auto *pg_42 = buffer.data(pg + 42);
    const auto *pg_44 = buffer.data(pg + 44);

    const auto *dd_s_0 = buffer.data(dd_s + 0);
    const auto *dd_s_3 = buffer.data(dd_s + 3);
    const auto *dd_s_5 = buffer.data(dd_s + 5);
    const auto *dd_s_18 = buffer.data(dd_s + 18);
    const auto *dd_s_19 = buffer.data(dd_s + 19);
    const auto *dd_s_21 = buffer.data(dd_s + 21);
    const auto *dd_s_23 = buffer.data(dd_s + 23);
    const auto *dd_s_30 = buffer.data(dd_s + 30);
    const auto *dd_s_32 = buffer.data(dd_s + 32);
    const auto *dd_s_33 = buffer.data(dd_s + 33);
    const auto *dd_s_34 = buffer.data(dd_s + 34);
    const auto *dd_s_35 = buffer.data(dd_s + 35);

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
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_18 = buffer.data(dd + 18);
    const auto *dd_19 = buffer.data(dd + 19);
    const auto *dd_21 = buffer.data(dd + 21);
    const auto *dd_23 = buffer.data(dd + 23);
    const auto *dd_30 = buffer.data(dd + 30);
    const auto *dd_32 = buffer.data(dd + 32);
    const auto *dd_33 = buffer.data(dd + 33);
    const auto *dd_34 = buffer.data(dd + 34);
    const auto *dd_35 = buffer.data(dd + 35);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_27 = buffer.data(df + 27);
    const auto *df_29 = buffer.data(df + 29);
    const auto *df_30 = buffer.data(df + 30);
    const auto *df_31 = buffer.data(df + 31);
    const auto *df_33 = buffer.data(df + 33);
    const auto *df_35 = buffer.data(df + 35);
    const auto *df_36 = buffer.data(df + 36);
    const auto *df_37 = buffer.data(df + 37);
    const auto *df_38 = buffer.data(df + 38);
    const auto *df_39 = buffer.data(df + 39);
    const auto *df_46 = buffer.data(df + 46);
    const auto *df_47 = buffer.data(df + 47);
    const auto *df_48 = buffer.data(df + 48);
    const auto *df_49 = buffer.data(df + 49);
    const auto *df_50 = buffer.data(df + 50);
    const auto *df_52 = buffer.data(df + 52);
    const auto *df_53 = buffer.data(df + 53);
    const auto *df_55 = buffer.data(df + 55);
    const auto *df_56 = buffer.data(df + 56);
    const auto *df_57 = buffer.data(df + 57);
    const auto *df_58 = buffer.data(df + 58);
    const auto *df_59 = buffer.data(df + 59);

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

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, pb_y, pb_z, pf_6, dd_s_0, dg_s_4, dg_s_5, \
                         dg_s_6, dd_0, df_2, df_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * dg_s_4[k]
                 + pb_y[k] * df_2[k];

        t_5[k] = -f_4 * dd_s_0[k]
                 + f_2 * dg_s_5[k]
                 + f_5 * dd_0[k]
                 + pb_z[k] * df_2[k];

        t_6[k] = f_0 * pf_6[k]
                 + f_2 * dg_s_6[k]
                 + pb_x[k] * df_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_x, pb_y, pb_z, pf_9, dg_s_7, dg_s_8, dg_s_9, df_3, \
                         df_5, df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * dg_s_7[k]
                 + pb_z[k] * df_3[k];

        t_8[k] = f_2 * dg_s_8[k]
                 + pb_y[k] * df_5[k];

        t_9[k] = f_0 * pf_9[k]
                 + f_2 * dg_s_9[k]
                 + pb_x[k] * df_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, pb_z, dd_s_3, dd_s_5, dg_s_10, dg_s_11, \
                         dg_s_12, dd_3, dd_5, df_6, df_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * dd_s_3[k]
                  + f_2 * dg_s_10[k]
                  + f_3 * dd_3[k]
                  + pb_y[k] * df_6[k];

        t_11[k] = f_2 * dg_s_11[k]
                  + pb_z[k] * df_6[k];

        t_12[k] = -f_4 * dd_s_5[k]
                  + f_2 * dg_s_12[k]
                  + f_5 * dd_5[k]
                  + pb_y[k] * df_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_y, pb_z, pg_0, dd_s_5, dg_s_13, dg_s_14, \
                         dg_s_15, dd_5, df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * dg_s_13[k]
                  + pb_y[k] * df_9[k];

        t_14[k] = -f_1 * dd_s_5[k]
                  + f_2 * dg_s_14[k]
                  + f_3 * dd_5[k]
                  + pb_z[k] * df_9[k];

        t_15[k] = pa_y[k] * pg_0[k]
                  + f_2 * dg_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, pb_z, pf_0, pf_13, pg_18, \
                         dg_s_16, dg_s_17, dg_s_18, dg_s_19, df_10, \
                         df_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * pf_0[k]
                  + f_2 * dg_s_16[k]
                  + pb_y[k] * df_10[k];

        t_17[k] = f_2 * dg_s_17[k]
                  + pb_z[k] * df_10[k];

        t_18[k] = f_0 * pf_13[k]
                  + pa_x[k] * pg_18[k]
                  + f_2 * dg_s_18[k];

        t_19[k] = f_2 * dg_s_19[k]
                  + pb_z[k] * df_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_x, pb_z, pf_16, pg_5, dg_s_20, dg_s_21, \
                         dg_s_22, df_13, df_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_y[k] * pg_5[k]
                  + f_2 * dg_s_20[k];

        t_21[k] = f_5 * pf_16[k]
                  + f_2 * dg_s_21[k]
                  + pb_x[k] * df_16[k];

        t_22[k] = f_2 * dg_s_22[k]
                  + pb_z[k] * df_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pa_y, pb_x, pf_18, pg_9, pg_25, dg_s_23, \
                         dg_s_24, dg_s_25, df_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * pf_18[k]
                  + f_2 * dg_s_23[k]
                  + pb_x[k] * df_18[k];

        t_24[k] = pa_y[k] * pg_9[k]
                  + f_2 * dg_s_24[k];

        t_25[k] = pa_x[k] * pg_25[k]
                  + f_2 * dg_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_x, pb_z, pg_27, pg_28, pg_29, dg_s_26, \
                         dg_s_27, dg_s_28, dg_s_29, df_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * dg_s_26[k]
                  + pb_z[k] * df_16[k];

        t_27[k] = pa_x[k] * pg_27[k]
                  + f_2 * dg_s_27[k];

        t_28[k] = pa_x[k] * pg_28[k]
                  + f_2 * dg_s_28[k];

        t_29[k] = pa_x[k] * pg_29[k]
                  + f_2 * dg_s_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_z, pb_y, pb_z, pf_0, pg_0, pg_3, dg_s_30, \
                         dg_s_31, dg_s_32, dg_s_33, df_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_z[k] * pg_0[k]
                  + f_2 * dg_s_30[k];

        t_31[k] = f_2 * dg_s_31[k]
                  + pb_y[k] * df_20[k];

        t_32[k] = f_5 * pf_0[k]
                  + f_2 * dg_s_32[k]
                  + pb_z[k] * df_20[k];

        t_33[k] = pa_z[k] * pg_3[k]
                  + f_2 * dg_s_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_x, pa_z, pb_y, pf_25, pg_6, pg_35, dg_s_34, \
                         dg_s_35, dg_s_36, df_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_2 * dg_s_34[k]
                  + pb_y[k] * df_22[k];

        t_35[k] = f_0 * pf_25[k]
                  + pa_x[k] * pg_35[k]
                  + f_2 * dg_s_35[k];

        t_36[k] = pa_z[k] * pg_6[k]
                  + f_2 * dg_s_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, pb_y, pf_27, pf_29, dg_s_37, dg_s_38, \
                         dg_s_39, df_25, df_27, df_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * pf_27[k]
                  + f_2 * dg_s_37[k]
                  + pb_x[k] * df_27[k];

        t_38[k] = f_2 * dg_s_38[k]
                  + pb_y[k] * df_25[k];

        t_39[k] = f_5 * pf_29[k]
                  + f_2 * dg_s_39[k]
                  + pb_x[k] * df_29[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pb_y, pg_40, pg_41, pg_42, dg_s_40, \
                         dg_s_41, dg_s_42, dg_s_43, df_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_x[k] * pg_40[k]
                  + f_2 * dg_s_40[k];

        t_41[k] = pa_x[k] * pg_41[k]
                  + f_2 * dg_s_41[k];

        t_42[k] = pa_x[k] * pg_42[k]
                  + f_2 * dg_s_42[k];

        t_43[k] = f_2 * dg_s_43[k]
                  + pb_y[k] * df_29[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_x, pb_x, pg_44, dd_s_18, dd_s_19, dg_s_44, \
                         dg_s_45, dg_s_46, dd_18, dd_19, df_30, df_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_x[k] * pg_44[k]
                  + f_2 * dg_s_44[k];

        t_45[k] = -f_1 * dd_s_18[k]
                  + f_2 * dg_s_45[k]
                  + f_3 * dd_18[k]
                  + pb_x[k] * df_30[k];

        t_46[k] = -f_6 * dd_s_19[k]
                  + f_2 * dg_s_46[k]
                  + f_0 * dd_19[k]
                  + pb_x[k] * df_31[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pb_x, pb_z, dd_s_21, dg_s_47, dg_s_48, dg_s_49, \
                         dd_21, df_30, df_31, df_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_2 * dg_s_47[k]
                  + pb_z[k] * df_30[k];

        t_48[k] = -f_4 * dd_s_21[k]
                  + f_2 * dg_s_48[k]
                  + f_5 * dd_21[k]
                  + pb_x[k] * df_33[k];

        t_49[k] = f_2 * dg_s_49[k]
                  + pb_z[k] * df_31[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pb_x, dd_s_23, dg_s_50, dg_s_51, dg_s_52, \
                         dg_s_53, dd_23, df_35, df_36, df_37, df_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -f_4 * dd_s_23[k]
                  + f_2 * dg_s_50[k]
                  + f_5 * dd_23[k]
                  + pb_x[k] * df_35[k];

        t_51[k] = f_2 * dg_s_51[k]
                  + pb_x[k] * df_36[k];

        t_52[k] = f_2 * dg_s_52[k]
                  + pb_x[k] * df_37[k];

        t_53[k] = f_2 * dg_s_53[k]
                  + pb_x[k] * df_38[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_x, pb_y, pb_z, pf_16, dd_s_21, dg_s_54, dg_s_55, \
                         dg_s_56, dd_21, df_36, df_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_2 * dg_s_54[k]
                  + pb_x[k] * df_39[k];

        t_55[k] = f_0 * pf_16[k]
                  - f_1 * dd_s_21[k]
                  + f_2 * dg_s_55[k]
                  + f_3 * dd_21[k]
                  + pb_y[k] * df_36[k];

        t_56[k] = f_2 * dg_s_56[k]
                  + pb_z[k] * df_36[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_y, pb_z, pf_19, dd_s_21, dd_s_23, dg_s_57, \
                         dg_s_58, dg_s_59, dd_21, dd_23, df_37, df_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = -f_4 * dd_s_21[k]
                  + f_2 * dg_s_57[k]
                  + f_5 * dd_21[k]
                  + pb_z[k] * df_37[k];

        t_58[k] = f_0 * pf_19[k]
                  + f_2 * dg_s_58[k]
                  + pb_y[k] * df_39[k];

        t_59[k] = -f_1 * dd_s_23[k]
                  + f_2 * dg_s_59[k]
                  + f_3 * dd_23[k]
                  + pb_z[k] * df_39[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pa_z, pg_16, pg_18, pg_30, pg_32, \
                         dg_s_60, dg_s_61, dg_s_62, dg_s_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * pg_30[k]
                  + f_2 * dg_s_60[k];

        t_61[k] = pa_z[k] * pg_16[k]
                  + f_2 * dg_s_61[k];

        t_62[k] = pa_y[k] * pg_32[k]
                  + f_2 * dg_s_62[k];

        t_63[k] = pa_z[k] * pg_18[k]
                  + f_2 * dg_s_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pf_22, pg_34, pg_35, dg_s_64, \
                         dg_s_65, dg_s_66, dg_s_67, df_46, df_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_5 * pf_22[k]
                  + pa_y[k] * pg_34[k]
                  + f_2 * dg_s_64[k];

        t_65[k] = pa_y[k] * pg_35[k]
                  + f_2 * dg_s_65[k];

        t_66[k] = f_2 * dg_s_66[k]
                  + pb_x[k] * df_46[k];

        t_67[k] = f_2 * dg_s_67[k]
                  + pb_x[k] * df_47[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_z, pb_x, pb_z, pf_16, pg_25, dg_s_68, \
                         dg_s_69, dg_s_70, dg_s_71, df_46, df_48, \
                         df_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_2 * dg_s_68[k]
                  + pb_x[k] * df_48[k];

        t_69[k] = f_2 * dg_s_69[k]
                  + pb_x[k] * df_49[k];

        t_70[k] = pa_z[k] * pg_25[k]
                  + f_2 * dg_s_70[k];

        t_71[k] = f_5 * pf_16[k]
                  + f_2 * dg_s_71[k]
                  + pb_z[k] * df_46[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_y, pb_y, pf_28, pf_29, pg_42, pg_44, dg_s_72, \
                         dg_s_73, dg_s_74, df_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_0 * pf_28[k]
                  + pa_y[k] * pg_42[k]
                  + f_2 * dg_s_72[k];

        t_73[k] = f_5 * pf_29[k]
                  + f_2 * dg_s_73[k]
                  + pb_y[k] * df_49[k];

        t_74[k] = pa_y[k] * pg_44[k]
                  + f_2 * dg_s_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_x, pb_y, dd_s_30, dd_s_32, dg_s_75, dg_s_76, \
                         dg_s_77, dd_30, dd_32, df_50, df_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -f_1 * dd_s_30[k]
                  + f_2 * dg_s_75[k]
                  + f_3 * dd_30[k]
                  + pb_x[k] * df_50[k];

        t_76[k] = f_2 * dg_s_76[k]
                  + pb_y[k] * df_50[k];

        t_77[k] = -f_6 * dd_s_32[k]
                  + f_2 * dg_s_77[k]
                  + f_0 * dd_32[k]
                  + pb_x[k] * df_52[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pb_x, pb_y, dd_s_33, dd_s_35, dg_s_78, dg_s_79, \
                         dg_s_80, dd_33, dd_35, df_52, df_53, df_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = -f_4 * dd_s_33[k]
                  + f_2 * dg_s_78[k]
                  + f_5 * dd_33[k]
                  + pb_x[k] * df_53[k];

        t_79[k] = f_2 * dg_s_79[k]
                  + pb_y[k] * df_52[k];

        t_80[k] = -f_4 * dd_s_35[k]
                  + f_2 * dg_s_80[k]
                  + f_5 * dd_35[k]
                  + pb_x[k] * df_55[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pb_x, dg_s_81, dg_s_82, dg_s_83, dg_s_84, \
                         df_56, df_57, df_58, df_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_2 * dg_s_81[k]
                  + pb_x[k] * df_56[k];

        t_82[k] = f_2 * dg_s_82[k]
                  + pb_x[k] * df_57[k];

        t_83[k] = f_2 * dg_s_83[k]
                  + pb_x[k] * df_58[k];

        t_84[k] = f_2 * dg_s_84[k]
                  + pb_x[k] * df_59[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pb_y, dd_s_33, dd_s_34, dd_s_35, dg_s_85, dg_s_86, \
                         dg_s_87, dd_33, dd_34, dd_35, df_56, df_57, \
                         df_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -f_1 * dd_s_33[k]
                  + f_2 * dg_s_85[k]
                  + f_3 * dd_33[k]
                  + pb_y[k] * df_56[k];

        t_86[k] = -f_6 * dd_s_34[k]
                  + f_2 * dg_s_86[k]
                  + f_0 * dd_34[k]
                  + pb_y[k] * df_57[k];

        t_87[k] = -f_4 * dd_s_35[k]
                  + f_2 * dg_s_87[k]
                  + f_5 * dd_35[k]
                  + pb_y[k] * df_58[k];
    }

#pragma omp simd aligned(t_88, t_89, pb_y, pb_z, pf_29, dd_s_35, dg_s_88, dg_s_89, dd_35, \
                         df_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_2 * dg_s_88[k]
                  + pb_y[k] * df_59[k];

        t_89[k] = f_0 * pf_29[k]
                  - f_1 * dd_s_35[k]
                  + f_2 * dg_s_89[k]
                  + f_3 * dd_35[k]
                  + pb_z[k] * df_59[k];
    }
}

}  // namespace simdkin
