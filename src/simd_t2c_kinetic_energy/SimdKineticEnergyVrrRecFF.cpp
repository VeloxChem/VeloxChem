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


#include "SimdKineticEnergyVrrRecFF.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_ff_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf_s, const size_t pf,
                                 const size_t dd, const size_t df, const size_t fp_s,
                                 const size_t ff_s, const size_t fp, const size_t fd,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = beta / p;
    const auto f_6 = alpha / p;

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

    const auto *pf_s_1 = buffer.data(pf_s + 1);
    const auto *pf_s_2 = buffer.data(pf_s + 2);

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
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_26 = buffer.data(df + 26);

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_1 = buffer.data(fp_s + 1);
    const auto *fp_s_2 = buffer.data(fp_s + 2);
    const auto *fp_s_4 = buffer.data(fp_s + 4);
    const auto *fp_s_7 = buffer.data(fp_s + 7);
    const auto *fp_s_8 = buffer.data(fp_s + 8);
    const auto *fp_s_9 = buffer.data(fp_s + 9);
    const auto *fp_s_10 = buffer.data(fp_s + 10);
    const auto *fp_s_13 = buffer.data(fp_s + 13);
    const auto *fp_s_14 = buffer.data(fp_s + 14);
    const auto *fp_s_15 = buffer.data(fp_s + 15);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_15 = buffer.data(ff_s + 15);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_30 = buffer.data(ff_s + 30);
    const auto *ff_s_31 = buffer.data(ff_s + 31);
    const auto *ff_s_32 = buffer.data(ff_s + 32);
    const auto *ff_s_33 = buffer.data(ff_s + 33);
    const auto *ff_s_34 = buffer.data(ff_s + 34);
    const auto *ff_s_35 = buffer.data(ff_s + 35);
    const auto *ff_s_36 = buffer.data(ff_s + 36);
    const auto *ff_s_37 = buffer.data(ff_s + 37);
    const auto *ff_s_38 = buffer.data(ff_s + 38);
    const auto *ff_s_39 = buffer.data(ff_s + 39);
    const auto *ff_s_40 = buffer.data(ff_s + 40);
    const auto *ff_s_41 = buffer.data(ff_s + 41);
    const auto *ff_s_42 = buffer.data(ff_s + 42);
    const auto *ff_s_43 = buffer.data(ff_s + 43);
    const auto *ff_s_44 = buffer.data(ff_s + 44);
    const auto *ff_s_45 = buffer.data(ff_s + 45);
    const auto *ff_s_46 = buffer.data(ff_s + 46);
    const auto *ff_s_47 = buffer.data(ff_s + 47);
    const auto *ff_s_48 = buffer.data(ff_s + 48);
    const auto *ff_s_49 = buffer.data(ff_s + 49);
    const auto *ff_s_50 = buffer.data(ff_s + 50);
    const auto *ff_s_51 = buffer.data(ff_s + 51);
    const auto *ff_s_52 = buffer.data(ff_s + 52);
    const auto *ff_s_53 = buffer.data(ff_s + 53);
    const auto *ff_s_54 = buffer.data(ff_s + 54);
    const auto *ff_s_55 = buffer.data(ff_s + 55);
    const auto *ff_s_56 = buffer.data(ff_s + 56);
    const auto *ff_s_57 = buffer.data(ff_s + 57);
    const auto *ff_s_58 = buffer.data(ff_s + 58);
    const auto *ff_s_59 = buffer.data(ff_s + 59);
    const auto *ff_s_60 = buffer.data(ff_s + 60);
    const auto *ff_s_61 = buffer.data(ff_s + 61);
    const auto *ff_s_62 = buffer.data(ff_s + 62);
    const auto *ff_s_63 = buffer.data(ff_s + 63);
    const auto *ff_s_64 = buffer.data(ff_s + 64);
    const auto *ff_s_65 = buffer.data(ff_s + 65);
    const auto *ff_s_66 = buffer.data(ff_s + 66);
    const auto *ff_s_67 = buffer.data(ff_s + 67);
    const auto *ff_s_68 = buffer.data(ff_s + 68);
    const auto *ff_s_69 = buffer.data(ff_s + 69);
    const auto *ff_s_70 = buffer.data(ff_s + 70);
    const auto *ff_s_71 = buffer.data(ff_s + 71);
    const auto *ff_s_72 = buffer.data(ff_s + 72);
    const auto *ff_s_73 = buffer.data(ff_s + 73);
    const auto *ff_s_74 = buffer.data(ff_s + 74);
    const auto *ff_s_75 = buffer.data(ff_s + 75);
    const auto *ff_s_76 = buffer.data(ff_s + 76);
    const auto *ff_s_77 = buffer.data(ff_s + 77);
    const auto *ff_s_78 = buffer.data(ff_s + 78);
    const auto *ff_s_79 = buffer.data(ff_s + 79);
    const auto *ff_s_80 = buffer.data(ff_s + 80);
    const auto *ff_s_81 = buffer.data(ff_s + 81);
    const auto *ff_s_82 = buffer.data(ff_s + 82);
    const auto *ff_s_83 = buffer.data(ff_s + 83);
    const auto *ff_s_84 = buffer.data(ff_s + 84);
    const auto *ff_s_85 = buffer.data(ff_s + 85);
    const auto *ff_s_86 = buffer.data(ff_s + 86);
    const auto *ff_s_87 = buffer.data(ff_s + 87);
    const auto *ff_s_88 = buffer.data(ff_s + 88);
    const auto *ff_s_89 = buffer.data(ff_s + 89);
    const auto *ff_s_90 = buffer.data(ff_s + 90);
    const auto *ff_s_91 = buffer.data(ff_s + 91);
    const auto *ff_s_92 = buffer.data(ff_s + 92);
    const auto *ff_s_93 = buffer.data(ff_s + 93);
    const auto *ff_s_94 = buffer.data(ff_s + 94);
    const auto *ff_s_95 = buffer.data(ff_s + 95);
    const auto *ff_s_96 = buffer.data(ff_s + 96);
    const auto *ff_s_97 = buffer.data(ff_s + 97);
    const auto *ff_s_98 = buffer.data(ff_s + 98);
    const auto *ff_s_99 = buffer.data(ff_s + 99);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dd_0, fp_s_0, ff_s_0, ff_s_1, \
                         ff_s_2, fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = f_2 * ff_s_1[k]
                 + pb_y[k] * fd_0[k];

        t_2[k] = f_2 * ff_s_2[k]
                 + pb_z[k] * fd_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, dd_1, dd_2, ff_s_3, ff_s_4, ff_s_5, fd_1, \
                         fd_2, fd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_0 * dd_1[k]
                 + f_2 * ff_s_3[k]
                 + pb_x[k] * fd_2[k];

        t_4[k] = f_2 * ff_s_4[k]
                 + pb_y[k] * fd_1[k];

        t_5[k] = f_0 * dd_2[k]
                 + f_2 * ff_s_5[k]
                 + pb_x[k] * fd_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_y, pb_z, fp_s_1, fp_s_2, ff_s_6, ff_s_7, \
                         ff_s_8, ff_s_9, fp_1, fp_2, fd_2, fd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_1 * fp_s_1[k]
                 + f_2 * ff_s_6[k]
                 + f_3 * fp_1[k]
                 + pb_y[k] * fd_2[k];

        t_7[k] = f_2 * ff_s_7[k]
                 + pb_z[k] * fd_2[k];

        t_8[k] = f_2 * ff_s_8[k]
                 + pb_y[k] * fd_3[k];

        t_9[k] = -f_1 * fp_s_2[k]
                 + f_2 * ff_s_9[k]
                 + f_3 * fp_2[k]
                 + pb_z[k] * fd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pb_y, pb_z, dd_0, df_0, ff_s_10, ff_s_11, \
                         ff_s_12, fd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * df_0[k]
                  + f_2 * ff_s_10[k];

        t_11[k] = f_4 * dd_0[k]
                  + f_2 * ff_s_11[k]
                  + pb_y[k] * fd_4[k];

        t_12[k] = f_2 * ff_s_12[k]
                  + pb_z[k] * fd_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_x, pb_z, dd_4, df_2, ff_s_13, ff_s_14, \
                         ff_s_15, fd_5, fd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * dd_4[k]
                  + f_2 * ff_s_13[k]
                  + pb_x[k] * fd_6[k];

        t_14[k] = f_2 * ff_s_14[k]
                  + pb_z[k] * fd_5[k];

        t_15[k] = pa_y[k] * df_2[k]
                  + f_2 * ff_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pb_y, pb_z, pf_s_1, pf_1, dd_2, df_7, \
                         ff_s_16, ff_s_17, ff_s_18, fd_6, fd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_5 * pf_s_1[k]
                  + f_4 * pf_1[k]
                  + pa_x[k] * df_7[k]
                  + f_2 * ff_s_16[k];

        t_17[k] = f_2 * ff_s_17[k]
                  + pb_z[k] * fd_6[k];

        t_18[k] = f_4 * dd_2[k]
                  + f_2 * ff_s_18[k]
                  + pb_y[k] * fd_7[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_y, pb_z, dd_0, df_0, df_4, \
                         ff_s_19, ff_s_20, ff_s_21, ff_s_22, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * df_4[k]
                  + f_2 * ff_s_19[k];

        t_20[k] = pa_z[k] * df_0[k]
                  + f_2 * ff_s_20[k];

        t_21[k] = f_2 * ff_s_21[k]
                  + pb_y[k] * fd_8[k];

        t_22[k] = f_4 * dd_0[k]
                  + f_2 * ff_s_22[k]
                  + pb_z[k] * fd_8[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pb_x, pb_y, dd_6, df_1, df_3, ff_s_23, \
                         ff_s_24, ff_s_25, ff_s_26, fd_9, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_z[k] * df_1[k]
                  + f_2 * ff_s_23[k];

        t_24[k] = f_2 * ff_s_24[k]
                  + pb_y[k] * fd_9[k];

        t_25[k] = f_3 * dd_6[k]
                  + f_2 * ff_s_25[k]
                  + pb_x[k] * fd_11[k];

        t_26[k] = pa_z[k] * df_3[k]
                  + f_2 * ff_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pb_y, pf_s_2, pf_2, df_11, fp_s_4, ff_s_27, \
                         ff_s_28, ff_s_29, fp_4, fd_10, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_6 * fp_s_4[k]
                  + f_2 * ff_s_27[k]
                  + f_4 * fp_4[k]
                  + pb_y[k] * fd_10[k];

        t_28[k] = f_2 * ff_s_28[k]
                  + pb_y[k] * fd_11[k];

        t_29[k] = -f_5 * pf_s_2[k]
                  + f_4 * pf_2[k]
                  + pa_x[k] * df_11[k]
                  + f_2 * ff_s_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_x, pb_y, pb_z, dd_3, dd_7, df_12, ff_s_30, \
                         ff_s_31, ff_s_32, fd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * dd_7[k]
                  + pa_x[k] * df_12[k]
                  + f_2 * ff_s_30[k];

        t_31[k] = f_3 * dd_3[k]
                  + f_2 * ff_s_31[k]
                  + pb_y[k] * fd_12[k];

        t_32[k] = f_2 * ff_s_32[k]
                  + pb_z[k] * fd_12[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_x, pb_z, dd_9, dd_10, ff_s_33, ff_s_34, ff_s_35, \
                         fd_13, fd_14, fd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_4 * dd_9[k]
                  + f_2 * ff_s_33[k]
                  + pb_x[k] * fd_14[k];

        t_34[k] = f_2 * ff_s_34[k]
                  + pb_z[k] * fd_13[k];

        t_35[k] = f_4 * dd_10[k]
                  + f_2 * ff_s_35[k]
                  + pb_x[k] * fd_15[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pb_z, df_14, df_15, df_16, ff_s_36, \
                         ff_s_37, ff_s_38, ff_s_39, fd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_x[k] * df_14[k]
                  + f_2 * ff_s_36[k];

        t_37[k] = f_2 * ff_s_37[k]
                  + pb_z[k] * fd_14[k];

        t_38[k] = pa_x[k] * df_15[k]
                  + f_2 * ff_s_38[k];

        t_39[k] = pa_x[k] * df_16[k]
                  + f_2 * ff_s_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pa_z, df_5, df_6, df_8, df_9, ff_s_40, \
                         ff_s_41, ff_s_42, ff_s_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * df_8[k]
                  + f_2 * ff_s_40[k];

        t_41[k] = pa_z[k] * df_5[k]
                  + f_2 * ff_s_41[k];

        t_42[k] = pa_y[k] * df_9[k]
                  + f_2 * ff_s_42[k];

        t_43[k] = pa_z[k] * df_6[k]
                  + f_2 * ff_s_43[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pa_y, pb_x, dd_12, df_10, df_17, df_18, \
                         ff_s_44, ff_s_45, ff_s_46, ff_s_47, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_4 * dd_12[k]
                  + f_2 * ff_s_44[k]
                  + pb_x[k] * fd_16[k];

        t_45[k] = pa_y[k] * df_10[k]
                  + f_2 * ff_s_45[k];

        t_46[k] = pa_x[k] * df_17[k]
                  + f_2 * ff_s_46[k];

        t_47[k] = pa_x[k] * df_18[k]
                  + f_2 * ff_s_47[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pb_y, dd_14, df_19, df_20, df_21, \
                         ff_s_48, ff_s_49, ff_s_50, ff_s_51, fd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_x[k] * df_19[k]
                  + f_2 * ff_s_48[k];

        t_49[k] = pa_x[k] * df_20[k]
                  + f_2 * ff_s_49[k];

        t_50[k] = f_0 * dd_14[k]
                  + pa_x[k] * df_21[k]
                  + f_2 * ff_s_50[k];

        t_51[k] = f_2 * ff_s_51[k]
                  + pb_y[k] * fd_17[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_x, pb_y, pb_z, dd_5, dd_16, ff_s_52, ff_s_53, \
                         ff_s_54, fd_17, fd_18, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * dd_5[k]
                  + f_2 * ff_s_52[k]
                  + pb_z[k] * fd_17[k];

        t_53[k] = f_4 * dd_16[k]
                  + f_2 * ff_s_53[k]
                  + pb_x[k] * fd_19[k];

        t_54[k] = f_2 * ff_s_54[k]
                  + pb_y[k] * fd_18[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_x, pb_x, pb_y, dd_18, df_24, df_25, \
                         ff_s_55, ff_s_56, ff_s_57, ff_s_58, fd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_4 * dd_18[k]
                  + f_2 * ff_s_55[k]
                  + pb_x[k] * fd_20[k];

        t_56[k] = pa_x[k] * df_24[k]
                  + f_2 * ff_s_56[k];

        t_57[k] = pa_x[k] * df_25[k]
                  + f_2 * ff_s_57[k];

        t_58[k] = f_2 * ff_s_58[k]
                  + pb_y[k] * fd_20[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pa_x, pb_x, df_26, fp_s_7, fp_s_8, ff_s_59, \
                         ff_s_60, ff_s_61, fp_7, fp_8, fd_21, fd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_x[k] * df_26[k]
                  + f_2 * ff_s_59[k];

        t_60[k] = -f_1 * fp_s_7[k]
                  + f_2 * ff_s_60[k]
                  + f_3 * fp_7[k]
                  + pb_x[k] * fd_21[k];

        t_61[k] = -f_6 * fp_s_8[k]
                  + f_2 * ff_s_61[k]
                  + f_4 * fp_8[k]
                  + pb_x[k] * fd_22[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pb_x, pb_z, ff_s_62, ff_s_63, ff_s_64, \
                         ff_s_65, fd_21, fd_23, fd_24, fd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_2 * ff_s_62[k]
                  + pb_z[k] * fd_21[k];

        t_63[k] = f_2 * ff_s_63[k]
                  + pb_x[k] * fd_23[k];

        t_64[k] = f_2 * ff_s_64[k]
                  + pb_x[k] * fd_24[k];

        t_65[k] = f_2 * ff_s_65[k]
                  + pb_x[k] * fd_25[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_y, pb_z, dd_9, dd_10, fp_s_8, ff_s_66, ff_s_67, \
                         ff_s_68, fp_8, fd_23, fd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_0 * dd_9[k]
                  - f_1 * fp_s_8[k]
                  + f_2 * ff_s_66[k]
                  + f_3 * fp_8[k]
                  + pb_y[k] * fd_23[k];

        t_67[k] = f_2 * ff_s_67[k]
                  + pb_z[k] * fd_23[k];

        t_68[k] = f_0 * dd_10[k]
                  + f_2 * ff_s_68[k]
                  + pb_y[k] * fd_25[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_z, pb_z, df_12, df_13, fp_s_9, ff_s_69, ff_s_70, \
                         ff_s_71, fp_9, fd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = -f_1 * fp_s_9[k]
                  + f_2 * ff_s_69[k]
                  + f_3 * fp_9[k]
                  + pb_z[k] * fd_25[k];

        t_70[k] = pa_z[k] * df_12[k]
                  + f_2 * ff_s_70[k];

        t_71[k] = pa_z[k] * df_13[k]
                  + f_2 * ff_s_71[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pb_x, fp_s_10, ff_s_72, ff_s_73, ff_s_74, \
                         ff_s_75, fp_10, fd_26, fd_27, fd_28, fd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -f_6 * fp_s_10[k]
                  + f_2 * ff_s_72[k]
                  + f_4 * fp_10[k]
                  + pb_x[k] * fd_26[k];

        t_73[k] = f_2 * ff_s_73[k]
                  + pb_x[k] * fd_27[k];

        t_74[k] = f_2 * ff_s_74[k]
                  + pb_x[k] * fd_28[k];

        t_75[k] = f_2 * ff_s_75[k]
                  + pb_x[k] * fd_29[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pa_z, pb_y, pb_z, dd_9, dd_13, df_14, ff_s_76, \
                         ff_s_77, ff_s_78, fd_27, fd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pa_z[k] * df_14[k]
                  + f_2 * ff_s_76[k];

        t_77[k] = f_4 * dd_9[k]
                  + f_2 * ff_s_77[k]
                  + pb_z[k] * fd_27[k];

        t_78[k] = f_3 * dd_13[k]
                  + f_2 * ff_s_78[k]
                  + pb_y[k] * fd_29[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_y, pf_s_2, pf_2, dd_14, df_20, df_21, \
                         df_22, df_23, ff_s_79, ff_s_80, ff_s_81, \
                         ff_s_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = -f_5 * pf_s_2[k]
                  + f_4 * pf_2[k]
                  + pa_y[k] * df_20[k]
                  + f_2 * ff_s_79[k];

        t_80[k] = pa_y[k] * df_21[k]
                  + f_2 * ff_s_80[k];

        t_81[k] = f_4 * dd_14[k]
                  + pa_y[k] * df_22[k]
                  + f_2 * ff_s_81[k];

        t_82[k] = pa_y[k] * df_23[k]
                  + f_2 * ff_s_82[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pb_x, dd_16, df_24, ff_s_83, ff_s_84, \
                         ff_s_85, ff_s_86, fd_30, fd_31, fd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_2 * ff_s_83[k]
                  + pb_x[k] * fd_30[k];

        t_84[k] = f_2 * ff_s_84[k]
                  + pb_x[k] * fd_31[k];

        t_85[k] = f_2 * ff_s_85[k]
                  + pb_x[k] * fd_32[k];

        t_86[k] = f_0 * dd_16[k]
                  + pa_y[k] * df_24[k]
                  + f_2 * ff_s_86[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pa_y, pb_y, pb_z, dd_11, dd_18, df_26, ff_s_87, \
                         ff_s_88, ff_s_89, fd_30, fd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * dd_11[k]
                  + f_2 * ff_s_87[k]
                  + pb_z[k] * fd_30[k];

        t_88[k] = f_4 * dd_18[k]
                  + f_2 * ff_s_88[k]
                  + pb_y[k] * fd_32[k];

        t_89[k] = pa_y[k] * df_26[k]
                  + f_2 * ff_s_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_x, pb_y, fp_s_13, fp_s_15, ff_s_90, ff_s_91, \
                         ff_s_92, fp_12, fp_14, fd_33, fd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -f_1 * fp_s_13[k]
                  + f_2 * ff_s_90[k]
                  + f_3 * fp_12[k]
                  + pb_x[k] * fd_33[k];

        t_91[k] = f_2 * ff_s_91[k]
                  + pb_y[k] * fd_33[k];

        t_92[k] = -f_6 * fp_s_15[k]
                  + f_2 * ff_s_92[k]
                  + f_4 * fp_14[k]
                  + pb_x[k] * fd_34[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pb_x, pb_y, fp_s_14, ff_s_93, ff_s_94, \
                         ff_s_95, ff_s_96, fp_13, fd_35, fd_36, fd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_2 * ff_s_93[k]
                  + pb_x[k] * fd_35[k];

        t_94[k] = f_2 * ff_s_94[k]
                  + pb_x[k] * fd_36[k];

        t_95[k] = f_2 * ff_s_95[k]
                  + pb_x[k] * fd_37[k];

        t_96[k] = -f_1 * fp_s_14[k]
                  + f_2 * ff_s_96[k]
                  + f_3 * fp_13[k]
                  + pb_y[k] * fd_35[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pb_y, pb_z, dd_18, fp_s_15, ff_s_97, ff_s_98, \
                         ff_s_99, fp_14, fd_36, fd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = -f_6 * fp_s_15[k]
                  + f_2 * ff_s_97[k]
                  + f_4 * fp_14[k]
                  + pb_y[k] * fd_36[k];

        t_98[k] = f_2 * ff_s_98[k]
                  + pb_y[k] * fd_37[k];

        t_99[k] = f_0 * dd_18[k]
                  - f_1 * fp_s_15[k]
                  + f_2 * ff_s_99[k]
                  + f_3 * fp_14[k]
                  + pb_z[k] * fd_37[k];
    }
}

auto
compute_prim_ff_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf_s, const size_t pf,
                                 const size_t dd, const size_t df, const size_t fp_s,
                                 const size_t ff_s, const size_t fp, const size_t fd,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = beta / p;
    const auto f_6 = alpha / p;

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

    const auto *pf_s_4 = buffer.data(pf_s + 4);
    const auto *pf_s_10 = buffer.data(pf_s + 10);

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

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_1 = buffer.data(fp_s + 1);
    const auto *fp_s_2 = buffer.data(fp_s + 2);
    const auto *fp_s_3 = buffer.data(fp_s + 3);
    const auto *fp_s_4 = buffer.data(fp_s + 4);
    const auto *fp_s_5 = buffer.data(fp_s + 5);
    const auto *fp_s_6 = buffer.data(fp_s + 6);
    const auto *fp_s_7 = buffer.data(fp_s + 7);
    const auto *fp_s_10 = buffer.data(fp_s + 10);
    const auto *fp_s_11 = buffer.data(fp_s + 11);
    const auto *fp_s_12 = buffer.data(fp_s + 12);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_15 = buffer.data(ff_s + 15);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_30 = buffer.data(ff_s + 30);
    const auto *ff_s_31 = buffer.data(ff_s + 31);
    const auto *ff_s_33 = buffer.data(ff_s + 33);
    const auto *ff_s_34 = buffer.data(ff_s + 34);
    const auto *ff_s_35 = buffer.data(ff_s + 35);
    const auto *ff_s_36 = buffer.data(ff_s + 36);
    const auto *ff_s_37 = buffer.data(ff_s + 37);
    const auto *ff_s_38 = buffer.data(ff_s + 38);
    const auto *ff_s_39 = buffer.data(ff_s + 39);
    const auto *ff_s_40 = buffer.data(ff_s + 40);
    const auto *ff_s_41 = buffer.data(ff_s + 41);
    const auto *ff_s_42 = buffer.data(ff_s + 42);
    const auto *ff_s_43 = buffer.data(ff_s + 43);
    const auto *ff_s_44 = buffer.data(ff_s + 44);
    const auto *ff_s_45 = buffer.data(ff_s + 45);
    const auto *ff_s_46 = buffer.data(ff_s + 46);
    const auto *ff_s_48 = buffer.data(ff_s + 48);
    const auto *ff_s_49 = buffer.data(ff_s + 49);
    const auto *ff_s_50 = buffer.data(ff_s + 50);
    const auto *ff_s_51 = buffer.data(ff_s + 51);
    const auto *ff_s_52 = buffer.data(ff_s + 52);
    const auto *ff_s_56 = buffer.data(ff_s + 56);
    const auto *ff_s_57 = buffer.data(ff_s + 57);
    const auto *ff_s_58 = buffer.data(ff_s + 58);
    const auto *ff_s_59 = buffer.data(ff_s + 59);
    const auto *ff_s_60 = buffer.data(ff_s + 60);
    const auto *ff_s_61 = buffer.data(ff_s + 61);
    const auto *ff_s_62 = buffer.data(ff_s + 62);
    const auto *ff_s_63 = buffer.data(ff_s + 63);
    const auto *ff_s_64 = buffer.data(ff_s + 64);
    const auto *ff_s_65 = buffer.data(ff_s + 65);
    const auto *ff_s_66 = buffer.data(ff_s + 66);
    const auto *ff_s_67 = buffer.data(ff_s + 67);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dd_0, fp_s_0, ff_s_0, ff_s_1, \
                         ff_s_2, fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = f_2 * ff_s_1[k]
                 + pb_y[k] * fd_0[k];

        t_2[k] = f_2 * ff_s_2[k]
                 + pb_z[k] * fd_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pb_x, pb_y, dd_1, dd_2, fp_s_1, ff_s_3, ff_s_4, \
                         ff_s_5, ff_s_6, fp_1, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_0 * dd_1[k]
                 + f_2 * ff_s_3[k]
                 + pb_x[k] * fd_1[k];

        t_4[k] = f_0 * dd_2[k]
                 + f_2 * ff_s_4[k]
                 + pb_x[k] * fd_2[k];

        t_5[k] = -f_1 * fp_s_1[k]
                 + f_2 * ff_s_5[k]
                 + f_3 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_6[k] = f_2 * ff_s_6[k]
                 + pb_y[k] * fd_2[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_y, pb_y, pb_z, dd_0, df_0, fp_s_2, ff_s_7, ff_s_8, \
                         ff_s_9, fp_2, fd_2, fd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -f_1 * fp_s_2[k]
                 + f_2 * ff_s_7[k]
                 + f_3 * fp_2[k]
                 + pb_z[k] * fd_2[k];

        t_8[k] = pa_y[k] * df_0[k]
                 + f_2 * ff_s_8[k];

        t_9[k] = f_4 * dd_0[k]
                 + f_2 * ff_s_9[k]
                 + pb_y[k] * fd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pb_x, pb_z, pf_s_4, pf_4, dd_4, df_6, \
                         ff_s_10, ff_s_11, ff_s_12, fd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * dd_4[k]
                  + f_2 * ff_s_10[k]
                  + pb_x[k] * fd_4[k];

        t_11[k] = -f_5 * pf_s_4[k]
                  + f_4 * pf_4[k]
                  + pa_x[k] * df_6[k]
                  + f_2 * ff_s_11[k];

        t_12[k] = f_2 * ff_s_12[k]
                  + pb_z[k] * fd_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pa_z, pb_y, dd_2, df_0, df_4, ff_s_13, \
                         ff_s_14, ff_s_15, fd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_4 * dd_2[k]
                  + f_2 * ff_s_13[k]
                  + pb_y[k] * fd_5[k];

        t_14[k] = pa_y[k] * df_4[k]
                  + f_2 * ff_s_14[k];

        t_15[k] = pa_z[k] * df_0[k]
                  + f_2 * ff_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, pb_y, pb_z, dd_0, dd_6, fp_s_3, ff_s_16, \
                         ff_s_17, ff_s_18, fp_3, fd_6, fd_7, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * dd_0[k]
                  + f_2 * ff_s_16[k]
                  + pb_z[k] * fd_6[k];

        t_17[k] = f_3 * dd_6[k]
                  + f_2 * ff_s_17[k]
                  + pb_x[k] * fd_8[k];

        t_18[k] = -f_6 * fp_s_3[k]
                  + f_2 * ff_s_18[k]
                  + f_4 * fp_3[k]
                  + pb_y[k] * fd_7[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pb_y, pf_s_10, pf_10, dd_7, df_11, df_12, \
                         ff_s_19, ff_s_20, ff_s_21, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * ff_s_19[k]
                  + pb_y[k] * fd_8[k];

        t_20[k] = -f_5 * pf_s_10[k]
                  + f_4 * pf_10[k]
                  + pa_x[k] * df_11[k]
                  + f_2 * ff_s_20[k];

        t_21[k] = f_0 * dd_7[k]
                  + pa_x[k] * df_12[k]
                  + f_2 * ff_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pb_x, pb_y, dd_3, dd_9, df_16, ff_s_22, \
                         ff_s_24, ff_s_25, fd_9, fd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * dd_3[k]
                  + f_2 * ff_s_22[k]
                  + pb_y[k] * fd_9[k];

        t_23[k] = f_4 * dd_9[k]
                  + f_2 * ff_s_24[k]
                  + pb_x[k] * fd_10[k];

        t_24[k] = pa_x[k] * df_16[k]
                  + f_2 * ff_s_25[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pa_y, df_9, df_18, df_19, df_21, \
                         ff_s_26, ff_s_27, ff_s_28, ff_s_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_x[k] * df_18[k]
                  + f_2 * ff_s_26[k];

        t_26[k] = pa_x[k] * df_19[k]
                  + f_2 * ff_s_27[k];

        t_27[k] = pa_y[k] * df_9[k]
                  + f_2 * ff_s_28[k];

        t_28[k] = pa_x[k] * df_21[k]
                  + f_2 * ff_s_29[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_x, pb_z, dd_5, dd_13, df_22, df_24, ff_s_30, \
                         ff_s_31, ff_s_33, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_x[k] * df_22[k]
                  + f_2 * ff_s_30[k];

        t_30[k] = f_0 * dd_13[k]
                  + pa_x[k] * df_24[k]
                  + f_2 * ff_s_31[k];

        t_31[k] = f_3 * dd_5[k]
                  + f_2 * ff_s_33[k]
                  + pb_z[k] * fd_11[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, dd_17, df_29, df_30, df_32, \
                         ff_s_34, ff_s_35, ff_s_36, ff_s_37, fd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_4 * dd_17[k]
                  + f_2 * ff_s_34[k]
                  + pb_x[k] * fd_12[k];

        t_33[k] = pa_x[k] * df_29[k]
                  + f_2 * ff_s_35[k];

        t_34[k] = pa_x[k] * df_30[k]
                  + f_2 * ff_s_36[k];

        t_35[k] = pa_x[k] * df_32[k]
                  + f_2 * ff_s_37[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_x, fp_s_4, fp_s_5, ff_s_38, ff_s_39, ff_s_40, \
                         fp_4, fp_5, fd_13, fd_14, fd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_1 * fp_s_4[k]
                  + f_2 * ff_s_38[k]
                  + f_3 * fp_4[k]
                  + pb_x[k] * fd_13[k];

        t_37[k] = -f_6 * fp_s_5[k]
                  + f_2 * ff_s_39[k]
                  + f_4 * fp_5[k]
                  + pb_x[k] * fd_14[k];

        t_38[k] = f_2 * ff_s_40[k]
                  + pb_x[k] * fd_15[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_x, pb_y, pb_z, dd_9, fp_s_5, ff_s_41, ff_s_42, \
                         ff_s_43, fp_5, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_2 * ff_s_41[k]
                  + pb_x[k] * fd_16[k];

        t_40[k] = f_0 * dd_9[k]
                  - f_1 * fp_s_5[k]
                  + f_2 * ff_s_42[k]
                  + f_3 * fp_5[k]
                  + pb_y[k] * fd_15[k];

        t_41[k] = f_2 * ff_s_43[k]
                  + pb_z[k] * fd_15[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, pb_y, pb_z, dd_10, fp_s_6, fp_s_7, ff_s_44, \
                         ff_s_45, ff_s_46, fp_6, fp_7, fd_16, fd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_0 * dd_10[k]
                  + f_2 * ff_s_44[k]
                  + pb_y[k] * fd_16[k];

        t_43[k] = -f_1 * fp_s_6[k]
                  + f_2 * ff_s_45[k]
                  + f_3 * fp_6[k]
                  + pb_z[k] * fd_16[k];

        t_44[k] = -f_6 * fp_s_7[k]
                  + f_2 * ff_s_46[k]
                  + f_4 * fp_7[k]
                  + pb_x[k] * fd_17[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_z, pb_x, pb_z, dd_9, df_16, ff_s_48, ff_s_49, \
                         ff_s_50, fd_18, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_2 * ff_s_48[k]
                  + pb_x[k] * fd_19[k];

        t_46[k] = pa_z[k] * df_16[k]
                  + f_2 * ff_s_49[k];

        t_47[k] = f_4 * dd_9[k]
                  + f_2 * ff_s_50[k]
                  + pb_z[k] * fd_18[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pa_y, pb_y, pf_s_10, pf_10, dd_12, dd_15, df_23, \
                         df_29, ff_s_51, ff_s_52, ff_s_56, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * dd_12[k]
                  + f_2 * ff_s_51[k]
                  + pb_y[k] * fd_19[k];

        t_49[k] = -f_5 * pf_s_10[k]
                  + f_4 * pf_10[k]
                  + pa_y[k] * df_23[k]
                  + f_2 * ff_s_52[k];

        t_50[k] = f_0 * dd_15[k]
                  + pa_y[k] * df_29[k]
                  + f_2 * ff_s_56[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pa_y, pb_y, pb_z, dd_11, dd_17, df_32, ff_s_57, \
                         ff_s_58, ff_s_59, fd_20, fd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * dd_11[k]
                  + f_2 * ff_s_57[k]
                  + pb_z[k] * fd_20[k];

        t_52[k] = f_4 * dd_17[k]
                  + f_2 * ff_s_58[k]
                  + pb_y[k] * fd_21[k];

        t_53[k] = pa_y[k] * df_32[k]
                  + f_2 * ff_s_59[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_x, fp_s_10, fp_s_12, ff_s_60, ff_s_61, ff_s_62, \
                         fp_8, fp_10, fd_22, fd_23, fd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = -f_1 * fp_s_10[k]
                  + f_2 * ff_s_60[k]
                  + f_3 * fp_8[k]
                  + pb_x[k] * fd_22[k];

        t_55[k] = -f_6 * fp_s_12[k]
                  + f_2 * ff_s_61[k]
                  + f_4 * fp_10[k]
                  + pb_x[k] * fd_23[k];

        t_56[k] = f_2 * ff_s_62[k]
                  + pb_x[k] * fd_24[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_x, pb_y, fp_s_11, fp_s_12, ff_s_63, ff_s_64, \
                         ff_s_65, fp_9, fp_10, fd_24, fd_25, fd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_2 * ff_s_63[k]
                  + pb_x[k] * fd_26[k];

        t_58[k] = -f_1 * fp_s_11[k]
                  + f_2 * ff_s_64[k]
                  + f_3 * fp_9[k]
                  + pb_y[k] * fd_24[k];

        t_59[k] = -f_6 * fp_s_12[k]
                  + f_2 * ff_s_65[k]
                  + f_4 * fp_10[k]
                  + pb_y[k] * fd_25[k];
    }

#pragma omp simd aligned(t_60, t_61, pb_y, pb_z, dd_17, fp_s_12, ff_s_66, ff_s_67, fp_10, \
                         fd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_2 * ff_s_66[k]
                  + pb_y[k] * fd_26[k];

        t_61[k] = f_0 * dd_17[k]
                  - f_1 * fp_s_12[k]
                  + f_2 * ff_s_67[k]
                  + f_3 * fp_10[k]
                  + pb_z[k] * fd_26[k];
    }
}

auto
compute_prim_ff_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf_s, const size_t pf,
                                 const size_t dd, const size_t df, const size_t fp_s,
                                 const size_t ff_s, const size_t fp, const size_t fd,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = alpha / p;

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

    const auto *pf_s_3 = buffer.data(pf_s + 3);
    const auto *pf_s_9 = buffer.data(pf_s + 9);

    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_6 = buffer.data(pf + 6);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_14 = buffer.data(dd + 14);

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

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_1 = buffer.data(fp_s + 1);
    const auto *fp_s_2 = buffer.data(fp_s + 2);
    const auto *fp_s_3 = buffer.data(fp_s + 3);
    const auto *fp_s_4 = buffer.data(fp_s + 4);
    const auto *fp_s_5 = buffer.data(fp_s + 5);
    const auto *fp_s_6 = buffer.data(fp_s + 6);
    const auto *fp_s_7 = buffer.data(fp_s + 7);
    const auto *fp_s_10 = buffer.data(fp_s + 10);
    const auto *fp_s_11 = buffer.data(fp_s + 11);
    const auto *fp_s_12 = buffer.data(fp_s + 12);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_15 = buffer.data(ff_s + 15);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_30 = buffer.data(ff_s + 30);
    const auto *ff_s_32 = buffer.data(ff_s + 32);
    const auto *ff_s_33 = buffer.data(ff_s + 33);
    const auto *ff_s_36 = buffer.data(ff_s + 36);
    const auto *ff_s_40 = buffer.data(ff_s + 40);
    const auto *ff_s_43 = buffer.data(ff_s + 43);
    const auto *ff_s_44 = buffer.data(ff_s + 44);
    const auto *ff_s_45 = buffer.data(ff_s + 45);
    const auto *ff_s_46 = buffer.data(ff_s + 46);
    const auto *ff_s_47 = buffer.data(ff_s + 47);
    const auto *ff_s_48 = buffer.data(ff_s + 48);
    const auto *ff_s_49 = buffer.data(ff_s + 49);
    const auto *ff_s_50 = buffer.data(ff_s + 50);
    const auto *ff_s_51 = buffer.data(ff_s + 51);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dd_0, fp_s_0, ff_s_0, ff_s_1, \
                         ff_s_2, fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = f_2 * ff_s_1[k]
                 + pb_y[k] * fd_0[k];

        t_2[k] = f_2 * ff_s_2[k]
                 + pb_z[k] * fd_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_y, pb_z, fp_s_1, fp_s_2, ff_s_3, ff_s_4, ff_s_5, \
                         fp_1, fp_2, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_1 * fp_s_1[k]
                 + f_2 * ff_s_3[k]
                 + f_3 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = f_2 * ff_s_4[k]
                 + pb_y[k] * fd_2[k];

        t_5[k] = -f_1 * fp_s_2[k]
                 + f_2 * ff_s_5[k]
                 + f_3 * fp_2[k]
                 + pb_z[k] * fd_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, pb_z, pf_s_3, pf_2, df_0, df_7, ff_s_6, \
                         ff_s_7, ff_s_8, fd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * df_0[k]
                 + f_2 * ff_s_6[k];

        t_7[k] = -f_4 * pf_s_3[k]
                 + f_5 * pf_2[k]
                 + pa_x[k] * df_7[k]
                 + f_2 * ff_s_7[k];

        t_8[k] = f_2 * ff_s_8[k]
                 + pb_z[k] * fd_3[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pa_z, pb_y, df_0, df_5, fp_s_3, ff_s_9, \
                         ff_s_10, ff_s_12, fp_3, fd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_y[k] * df_5[k]
                 + f_2 * ff_s_9[k];

        t_10[k] = pa_z[k] * df_0[k]
                  + f_2 * ff_s_10[k];

        t_11[k] = -f_6 * fp_s_3[k]
                  + f_2 * ff_s_12[k]
                  + f_5 * fp_3[k]
                  + pb_y[k] * fd_4[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_y, pf_s_9, pf_6, dd_5, df_9, df_10, \
                         ff_s_13, ff_s_14, ff_s_15, fd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_2 * ff_s_13[k]
                  + pb_y[k] * fd_5[k];

        t_13[k] = -f_4 * pf_s_9[k]
                  + f_5 * pf_6[k]
                  + pa_x[k] * df_9[k]
                  + f_2 * ff_s_14[k];

        t_14[k] = f_0 * dd_5[k]
                  + pa_x[k] * df_10[k]
                  + f_2 * ff_s_15[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, dd_10, df_14, df_20, df_27, ff_s_17, ff_s_18, \
                         ff_s_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_x[k] * df_14[k]
                  + f_2 * ff_s_17[k];

        t_16[k] = f_0 * dd_10[k]
                  + pa_x[k] * df_20[k]
                  + f_2 * ff_s_18[k];

        t_17[k] = pa_x[k] * df_27[k]
                  + f_2 * ff_s_21[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_x, fp_s_4, fp_s_5, ff_s_22, ff_s_23, ff_s_24, \
                         fp_4, fp_5, fd_6, fd_7, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -f_1 * fp_s_4[k]
                  + f_2 * ff_s_22[k]
                  + f_3 * fp_4[k]
                  + pb_x[k] * fd_6[k];

        t_19[k] = -f_6 * fp_s_5[k]
                  + f_2 * ff_s_23[k]
                  + f_5 * fp_5[k]
                  + pb_x[k] * fd_7[k];

        t_20[k] = f_2 * ff_s_24[k]
                  + pb_x[k] * fd_8[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_x, pb_y, pb_z, dd_7, fp_s_5, ff_s_25, ff_s_26, \
                         ff_s_27, fp_5, fd_8, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_2 * ff_s_25[k]
                  + pb_x[k] * fd_9[k];

        t_22[k] = f_0 * dd_7[k]
                  - f_1 * fp_s_5[k]
                  + f_2 * ff_s_26[k]
                  + f_3 * fp_5[k]
                  + pb_y[k] * fd_8[k];

        t_23[k] = f_2 * ff_s_27[k]
                  + pb_z[k] * fd_8[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_x, pb_y, pb_z, dd_8, fp_s_6, fp_s_7, ff_s_28, \
                         ff_s_29, ff_s_30, fp_6, fp_7, fd_9, fd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * dd_8[k]
                  + f_2 * ff_s_28[k]
                  + pb_y[k] * fd_9[k];

        t_25[k] = -f_1 * fp_s_6[k]
                  + f_2 * ff_s_29[k]
                  + f_3 * fp_6[k]
                  + pb_z[k] * fd_9[k];

        t_26[k] = -f_6 * fp_s_7[k]
                  + f_2 * ff_s_30[k]
                  + f_5 * fp_7[k]
                  + pb_x[k] * fd_10[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pa_z, pb_x, pf_s_9, pf_6, df_14, df_19, \
                         ff_s_32, ff_s_33, ff_s_36, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_2 * ff_s_32[k]
                  + pb_x[k] * fd_11[k];

        t_28[k] = pa_z[k] * df_14[k]
                  + f_2 * ff_s_33[k];

        t_29[k] = -f_4 * pf_s_9[k]
                  + f_5 * pf_6[k]
                  + pa_y[k] * df_19[k]
                  + f_2 * ff_s_36[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pb_x, dd_12, df_24, df_27, fp_s_10, ff_s_40, \
                         ff_s_43, ff_s_44, fp_8, fd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * dd_12[k]
                  + pa_y[k] * df_24[k]
                  + f_2 * ff_s_40[k];

        t_31[k] = pa_y[k] * df_27[k]
                  + f_2 * ff_s_43[k];

        t_32[k] = -f_1 * fp_s_10[k]
                  + f_2 * ff_s_44[k]
                  + f_3 * fp_8[k]
                  + pb_x[k] * fd_12[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_x, fp_s_12, ff_s_45, ff_s_46, ff_s_47, fp_10, \
                         fd_13, fd_14, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -f_6 * fp_s_12[k]
                  + f_2 * ff_s_45[k]
                  + f_5 * fp_10[k]
                  + pb_x[k] * fd_13[k];

        t_34[k] = f_2 * ff_s_46[k]
                  + pb_x[k] * fd_14[k];

        t_35[k] = f_2 * ff_s_47[k]
                  + pb_x[k] * fd_16[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_y, fp_s_11, fp_s_12, ff_s_48, ff_s_49, ff_s_50, \
                         fp_9, fp_10, fd_14, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_1 * fp_s_11[k]
                  + f_2 * ff_s_48[k]
                  + f_3 * fp_9[k]
                  + pb_y[k] * fd_14[k];

        t_37[k] = -f_6 * fp_s_12[k]
                  + f_2 * ff_s_49[k]
                  + f_5 * fp_10[k]
                  + pb_y[k] * fd_15[k];

        t_38[k] = f_2 * ff_s_50[k]
                  + pb_y[k] * fd_16[k];
    }

#pragma omp simd aligned(t_39, pb_z, dd_14, fp_s_12, ff_s_51, fp_10, \
                         fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * dd_14[k]
                  - f_1 * fp_s_12[k]
                  + f_2 * ff_s_51[k]
                  + f_3 * fp_10[k]
                  + pb_z[k] * fd_16[k];
    }
}

auto
compute_prim_ff_kinetic_energy_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf_s, const size_t pf,
                                 const size_t dd, const size_t df, const size_t fp_s,
                                 const size_t ff_s, const size_t fp, const size_t fd,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = alpha / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_s_3 = buffer.data(pf_s + 3);
    const auto *pf_s_9 = buffer.data(pf_s + 9);

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_4 = buffer.data(pf + 4);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_14 = buffer.data(dd + 14);

    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_16 = buffer.data(df + 16);

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_1 = buffer.data(fp_s + 1);
    const auto *fp_s_2 = buffer.data(fp_s + 2);
    const auto *fp_s_3 = buffer.data(fp_s + 3);
    const auto *fp_s_4 = buffer.data(fp_s + 4);
    const auto *fp_s_5 = buffer.data(fp_s + 5);
    const auto *fp_s_6 = buffer.data(fp_s + 6);
    const auto *fp_s_7 = buffer.data(fp_s + 7);
    const auto *fp_s_10 = buffer.data(fp_s + 10);
    const auto *fp_s_11 = buffer.data(fp_s + 11);
    const auto *fp_s_12 = buffer.data(fp_s + 12);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_32 = buffer.data(ff_s + 32);
    const auto *ff_s_40 = buffer.data(ff_s + 40);
    const auto *ff_s_41 = buffer.data(ff_s + 41);
    const auto *ff_s_42 = buffer.data(ff_s + 42);
    const auto *ff_s_43 = buffer.data(ff_s + 43);
    const auto *ff_s_44 = buffer.data(ff_s + 44);
    const auto *ff_s_45 = buffer.data(ff_s + 45);
    const auto *ff_s_46 = buffer.data(ff_s + 46);
    const auto *ff_s_47 = buffer.data(ff_s + 47);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dd_0, fp_s_0, ff_s_0, ff_s_1, \
                         ff_s_2, fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = f_2 * ff_s_1[k]
                 + pb_y[k] * fd_0[k];

        t_2[k] = f_2 * ff_s_2[k]
                 + pb_z[k] * fd_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_y, pb_z, fp_s_1, fp_s_2, ff_s_3, ff_s_4, ff_s_5, \
                         fp_1, fp_2, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_1 * fp_s_1[k]
                 + f_2 * ff_s_3[k]
                 + f_3 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = f_2 * ff_s_4[k]
                 + pb_y[k] * fd_2[k];

        t_5[k] = -f_1 * fp_s_2[k]
                 + f_2 * ff_s_5[k]
                 + f_3 * fp_2[k]
                 + pb_z[k] * fd_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pb_y, pb_z, pf_s_3, pf_1, df_6, fp_s_3, ff_s_6, \
                         ff_s_7, ff_s_10, fp_3, fd_3, fd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * pf_s_3[k]
                 + f_5 * pf_1[k]
                 + pa_x[k] * df_6[k]
                 + f_2 * ff_s_6[k];

        t_7[k] = f_2 * ff_s_7[k]
                 + pb_z[k] * fd_3[k];

        t_8[k] = -f_6 * fp_s_3[k]
                 + f_2 * ff_s_10[k]
                 + f_5 * fp_3[k]
                 + pb_y[k] * fd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_x, pb_y, pf_s_9, pf_4, df_7, fp_s_4, \
                         ff_s_11, ff_s_12, ff_s_18, fp_4, fd_5, fd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * ff_s_11[k]
                 + pb_y[k] * fd_5[k];

        t_10[k] = -f_4 * pf_s_9[k]
                  + f_5 * pf_4[k]
                  + pa_x[k] * df_7[k]
                  + f_2 * ff_s_12[k];

        t_11[k] = -f_1 * fp_s_4[k]
                  + f_2 * ff_s_18[k]
                  + f_3 * fp_4[k]
                  + pb_x[k] * fd_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_x, pb_y, dd_7, fp_s_5, ff_s_19, ff_s_20, \
                         ff_s_21, ff_s_22, fp_5, fd_7, fd_8, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_6 * fp_s_5[k]
                  + f_2 * ff_s_19[k]
                  + f_5 * fp_5[k]
                  + pb_x[k] * fd_7[k];

        t_13[k] = f_2 * ff_s_20[k]
                  + pb_x[k] * fd_8[k];

        t_14[k] = f_2 * ff_s_21[k]
                  + pb_x[k] * fd_9[k];

        t_15[k] = f_0 * dd_7[k]
                  - f_1 * fp_s_5[k]
                  + f_2 * ff_s_22[k]
                  + f_3 * fp_5[k]
                  + pb_y[k] * fd_8[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_y, pb_z, dd_8, fp_s_6, ff_s_23, ff_s_24, \
                         ff_s_25, fp_6, fd_8, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * ff_s_23[k]
                  + pb_z[k] * fd_8[k];

        t_17[k] = f_0 * dd_8[k]
                  + f_2 * ff_s_24[k]
                  + pb_y[k] * fd_9[k];

        t_18[k] = -f_1 * fp_s_6[k]
                  + f_2 * ff_s_25[k]
                  + f_3 * fp_6[k]
                  + pb_z[k] * fd_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pb_x, pf_s_9, pf_4, df_16, fp_s_7, ff_s_26, \
                         ff_s_28, ff_s_32, fp_7, fd_10, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_6 * fp_s_7[k]
                  + f_2 * ff_s_26[k]
                  + f_5 * fp_7[k]
                  + pb_x[k] * fd_10[k];

        t_20[k] = f_2 * ff_s_28[k]
                  + pb_x[k] * fd_11[k];

        t_21[k] = -f_4 * pf_s_9[k]
                  + f_5 * pf_4[k]
                  + pa_y[k] * df_16[k]
                  + f_2 * ff_s_32[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pb_x, fp_s_10, fp_s_12, ff_s_40, ff_s_41, ff_s_42, \
                         fp_8, fp_10, fd_12, fd_13, fd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -f_1 * fp_s_10[k]
                  + f_2 * ff_s_40[k]
                  + f_3 * fp_8[k]
                  + pb_x[k] * fd_12[k];

        t_23[k] = -f_6 * fp_s_12[k]
                  + f_2 * ff_s_41[k]
                  + f_5 * fp_10[k]
                  + pb_x[k] * fd_13[k];

        t_24[k] = f_2 * ff_s_42[k]
                  + pb_x[k] * fd_14[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pb_x, pb_y, fp_s_11, fp_s_12, ff_s_43, ff_s_44, \
                         ff_s_45, fp_9, fp_10, fd_14, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_2 * ff_s_43[k]
                  + pb_x[k] * fd_16[k];

        t_26[k] = -f_1 * fp_s_11[k]
                  + f_2 * ff_s_44[k]
                  + f_3 * fp_9[k]
                  + pb_y[k] * fd_14[k];

        t_27[k] = -f_6 * fp_s_12[k]
                  + f_2 * ff_s_45[k]
                  + f_5 * fp_10[k]
                  + pb_y[k] * fd_15[k];
    }

#pragma omp simd aligned(t_28, t_29, pb_y, pb_z, dd_14, fp_s_12, ff_s_46, ff_s_47, fp_10, \
                         fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_2 * ff_s_46[k]
                  + pb_y[k] * fd_16[k];

        t_29[k] = f_0 * dd_14[k]
                  - f_1 * fp_s_12[k]
                  + f_2 * ff_s_47[k]
                  + f_3 * fp_10[k]
                  + pb_z[k] * fd_16[k];
    }
}

auto
compute_prim_ff_kinetic_energy_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf_s, const size_t pf,
                                 const size_t dd, const size_t df, const size_t fp_s,
                                 const size_t ff_s, const size_t fp, const size_t fd,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = beta / p;
    const auto f_6 = alpha / p;

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

    const auto *pf_s_0 = buffer.data(pf_s + 0);
    const auto *pf_s_1 = buffer.data(pf_s + 1);

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

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_1 = buffer.data(fp_s + 1);
    const auto *fp_s_2 = buffer.data(fp_s + 2);
    const auto *fp_s_7 = buffer.data(fp_s + 7);
    const auto *fp_s_8 = buffer.data(fp_s + 8);
    const auto *fp_s_9 = buffer.data(fp_s + 9);
    const auto *fp_s_12 = buffer.data(fp_s + 12);
    const auto *fp_s_13 = buffer.data(fp_s + 13);
    const auto *fp_s_14 = buffer.data(fp_s + 14);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_15 = buffer.data(ff_s + 15);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_30 = buffer.data(ff_s + 30);
    const auto *ff_s_31 = buffer.data(ff_s + 31);
    const auto *ff_s_32 = buffer.data(ff_s + 32);
    const auto *ff_s_33 = buffer.data(ff_s + 33);
    const auto *ff_s_34 = buffer.data(ff_s + 34);
    const auto *ff_s_35 = buffer.data(ff_s + 35);
    const auto *ff_s_36 = buffer.data(ff_s + 36);
    const auto *ff_s_37 = buffer.data(ff_s + 37);
    const auto *ff_s_38 = buffer.data(ff_s + 38);
    const auto *ff_s_39 = buffer.data(ff_s + 39);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, dd_0, dd_1, dd_2, fp_s_0, ff_s_0, ff_s_1, \
                         ff_s_2, fp_0, fd_0, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = f_0 * dd_1[k]
                 + f_2 * ff_s_1[k]
                 + pb_x[k] * fd_1[k];

        t_2[k] = f_0 * dd_2[k]
                 + f_2 * ff_s_2[k]
                 + pb_x[k] * fd_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_y, pb_y, pb_z, df_0, fp_s_1, fp_s_2, ff_s_3, \
                         ff_s_4, ff_s_5, fp_1, fp_2, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_1 * fp_s_1[k]
                 + f_2 * ff_s_3[k]
                 + f_3 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = -f_1 * fp_s_2[k]
                 + f_2 * ff_s_4[k]
                 + f_3 * fp_2[k]
                 + pb_z[k] * fd_2[k];

        t_5[k] = pa_y[k] * df_0[k]
                 + f_2 * ff_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pb_x, pb_y, pf_s_0, pf_0, dd_0, dd_4, df_1, \
                         ff_s_6, ff_s_7, ff_s_8, fd_3, fd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_4 * dd_0[k]
                 + f_2 * ff_s_6[k]
                 + pb_y[k] * fd_3[k];

        t_7[k] = f_3 * dd_4[k]
                 + f_2 * ff_s_7[k]
                 + pb_x[k] * fd_4[k];

        t_8[k] = -f_5 * pf_s_0[k]
                 + f_4 * pf_0[k]
                 + pa_x[k] * df_1[k]
                 + f_2 * ff_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_z, pb_x, pb_z, dd_0, dd_6, df_0, ff_s_9, ff_s_10, \
                         ff_s_11, fd_6, fd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_z[k] * df_0[k]
                 + f_2 * ff_s_9[k];

        t_10[k] = f_4 * dd_0[k]
                  + f_2 * ff_s_10[k]
                  + pb_z[k] * fd_6[k];

        t_11[k] = f_3 * dd_6[k]
                  + f_2 * ff_s_11[k]
                  + pb_x[k] * fd_7[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_y, pf_s_1, pf_1, dd_3, dd_7, df_2, df_3, \
                         ff_s_12, ff_s_13, ff_s_14, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_5 * pf_s_1[k]
                  + f_4 * pf_1[k]
                  + pa_x[k] * df_2[k]
                  + f_2 * ff_s_12[k];

        t_13[k] = f_0 * dd_7[k]
                  + pa_x[k] * df_3[k]
                  + f_2 * ff_s_13[k];

        t_14[k] = f_3 * dd_3[k]
                  + f_2 * ff_s_14[k]
                  + pb_y[k] * fd_8[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pb_x, dd_8, dd_13, df_4, df_6, ff_s_15, \
                         ff_s_16, ff_s_17, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * dd_8[k]
                  + f_2 * ff_s_15[k]
                  + pb_x[k] * fd_9[k];

        t_16[k] = pa_x[k] * df_4[k]
                  + f_2 * ff_s_16[k];

        t_17[k] = f_0 * dd_13[k]
                  + pa_x[k] * df_6[k]
                  + f_2 * ff_s_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pb_x, pb_z, dd_5, dd_15, df_8, ff_s_18, \
                         ff_s_19, ff_s_20, fd_12, fd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * dd_5[k]
                  + f_2 * ff_s_18[k]
                  + pb_z[k] * fd_12[k];

        t_19[k] = f_4 * dd_15[k]
                  + f_2 * ff_s_19[k]
                  + pb_x[k] * fd_14[k];

        t_20[k] = pa_x[k] * df_8[k]
                  + f_2 * ff_s_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_x, pb_y, dd_8, fp_s_7, fp_s_8, ff_s_21, ff_s_22, \
                         ff_s_23, fp_7, fp_8, fd_15, fd_16, fd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_1 * fp_s_7[k]
                  + f_2 * ff_s_21[k]
                  + f_3 * fp_7[k]
                  + pb_x[k] * fd_15[k];

        t_22[k] = -f_6 * fp_s_8[k]
                  + f_2 * ff_s_22[k]
                  + f_4 * fp_8[k]
                  + pb_x[k] * fd_16[k];

        t_23[k] = f_0 * dd_8[k]
                  - f_1 * fp_s_8[k]
                  + f_2 * ff_s_23[k]
                  + f_3 * fp_8[k]
                  + pb_y[k] * fd_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_z, pb_y, pb_z, dd_9, df_4, fp_s_9, ff_s_24, \
                         ff_s_25, ff_s_26, fp_9, fd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * dd_9[k]
                  + f_2 * ff_s_24[k]
                  + pb_y[k] * fd_18[k];

        t_25[k] = -f_1 * fp_s_9[k]
                  + f_2 * ff_s_25[k]
                  + f_3 * fp_9[k]
                  + pb_z[k] * fd_18[k];

        t_26[k] = pa_z[k] * df_4[k]
                  + f_2 * ff_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_y, pb_z, pf_s_1, pf_1, dd_8, dd_12, df_5, \
                         ff_s_27, ff_s_28, ff_s_29, fd_19, fd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_4 * dd_8[k]
                  + f_2 * ff_s_27[k]
                  + pb_z[k] * fd_19[k];

        t_28[k] = f_3 * dd_12[k]
                  + f_2 * ff_s_28[k]
                  + pb_y[k] * fd_21[k];

        t_29[k] = -f_5 * pf_s_1[k]
                  + f_4 * pf_1[k]
                  + pa_y[k] * df_5[k]
                  + f_2 * ff_s_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pb_y, pb_z, dd_10, dd_14, dd_15, df_7, \
                         ff_s_30, ff_s_31, ff_s_32, fd_22, fd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * dd_14[k]
                  + pa_y[k] * df_7[k]
                  + f_2 * ff_s_30[k];

        t_31[k] = f_3 * dd_10[k]
                  + f_2 * ff_s_31[k]
                  + pb_z[k] * fd_22[k];

        t_32[k] = f_4 * dd_15[k]
                  + f_2 * ff_s_32[k]
                  + pb_y[k] * fd_24[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pb_x, pb_y, df_8, fp_s_12, ff_s_33, ff_s_34, \
                         ff_s_35, fp_12, fd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_y[k] * df_8[k]
                  + f_2 * ff_s_33[k];

        t_34[k] = -f_1 * fp_s_12[k]
                  + f_2 * ff_s_34[k]
                  + f_3 * fp_12[k]
                  + pb_x[k] * fd_25[k];

        t_35[k] = f_2 * ff_s_35[k]
                  + pb_y[k] * fd_25[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_x, pb_y, fp_s_13, fp_s_14, ff_s_36, ff_s_37, \
                         ff_s_38, fp_13, fp_14, fd_26, fd_27, fd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_6 * fp_s_14[k]
                  + f_2 * ff_s_36[k]
                  + f_4 * fp_14[k]
                  + pb_x[k] * fd_26[k];

        t_37[k] = -f_1 * fp_s_13[k]
                  + f_2 * ff_s_37[k]
                  + f_3 * fp_13[k]
                  + pb_y[k] * fd_27[k];

        t_38[k] = -f_6 * fp_s_14[k]
                  + f_2 * ff_s_38[k]
                  + f_4 * fp_14[k]
                  + pb_y[k] * fd_28[k];
    }

#pragma omp simd aligned(t_39, pb_z, dd_15, fp_s_14, ff_s_39, fp_14, \
                         fd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * dd_15[k]
                  - f_1 * fp_s_14[k]
                  + f_2 * ff_s_39[k]
                  + f_3 * fp_14[k]
                  + pb_z[k] * fd_29[k];
    }
}

auto
compute_prim_ff_kinetic_energy_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf_s, const size_t pf,
                                 const size_t dd, const size_t df, const size_t fp_s,
                                 const size_t ff_s, const size_t fp, const size_t fd,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = alpha / p;

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

    const auto *pf_s_1 = buffer.data(pf_s + 1);
    const auto *pf_s_3 = buffer.data(pf_s + 3);

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

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_1 = buffer.data(fp_s + 1);
    const auto *fp_s_2 = buffer.data(fp_s + 2);
    const auto *fp_s_3 = buffer.data(fp_s + 3);
    const auto *fp_s_4 = buffer.data(fp_s + 4);
    const auto *fp_s_5 = buffer.data(fp_s + 5);
    const auto *fp_s_6 = buffer.data(fp_s + 6);
    const auto *fp_s_9 = buffer.data(fp_s + 9);
    const auto *fp_s_10 = buffer.data(fp_s + 10);
    const auto *fp_s_11 = buffer.data(fp_s + 11);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_15 = buffer.data(ff_s + 15);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_30 = buffer.data(ff_s + 30);
    const auto *ff_s_31 = buffer.data(ff_s + 31);
    const auto *ff_s_32 = buffer.data(ff_s + 32);
    const auto *ff_s_33 = buffer.data(ff_s + 33);
    const auto *ff_s_34 = buffer.data(ff_s + 34);
    const auto *ff_s_35 = buffer.data(ff_s + 35);
    const auto *ff_s_36 = buffer.data(ff_s + 36);
    const auto *ff_s_37 = buffer.data(ff_s + 37);
    const auto *ff_s_38 = buffer.data(ff_s + 38);
    const auto *ff_s_39 = buffer.data(ff_s + 39);
    const auto *ff_s_40 = buffer.data(ff_s + 40);
    const auto *ff_s_41 = buffer.data(ff_s + 41);
    const auto *ff_s_42 = buffer.data(ff_s + 42);
    const auto *ff_s_43 = buffer.data(ff_s + 43);
    const auto *ff_s_44 = buffer.data(ff_s + 44);
    const auto *ff_s_45 = buffer.data(ff_s + 45);
    const auto *ff_s_46 = buffer.data(ff_s + 46);
    const auto *ff_s_47 = buffer.data(ff_s + 47);
    const auto *ff_s_48 = buffer.data(ff_s + 48);
    const auto *ff_s_49 = buffer.data(ff_s + 49);
    const auto *ff_s_50 = buffer.data(ff_s + 50);
    const auto *ff_s_51 = buffer.data(ff_s + 51);
    const auto *ff_s_52 = buffer.data(ff_s + 52);
    const auto *ff_s_53 = buffer.data(ff_s + 53);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dd_0, fp_s_0, ff_s_0, ff_s_1, \
                         ff_s_2, fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = f_2 * ff_s_1[k]
                 + pb_y[k] * fd_0[k];

        t_2[k] = f_2 * ff_s_2[k]
                 + pb_z[k] * fd_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_y, pb_y, pb_z, df_0, fp_s_1, fp_s_2, ff_s_3, \
                         ff_s_4, ff_s_5, fp_1, fp_2, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_1 * fp_s_1[k]
                 + f_2 * ff_s_3[k]
                 + f_3 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = -f_1 * fp_s_2[k]
                 + f_2 * ff_s_4[k]
                 + f_3 * fp_2[k]
                 + pb_z[k] * fd_2[k];

        t_5[k] = pa_y[k] * df_0[k]
                 + f_2 * ff_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, pb_y, pf_s_1, pf_1, dd_2, df_4, df_6, \
                         ff_s_6, ff_s_7, ff_s_8, fd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * pf_s_1[k]
                 + f_5 * pf_1[k]
                 + pa_x[k] * df_6[k]
                 + f_2 * ff_s_6[k];

        t_7[k] = f_5 * dd_2[k]
                 + f_2 * ff_s_7[k]
                 + pb_y[k] * fd_5[k];

        t_8[k] = pa_y[k] * df_4[k]
                 + f_2 * ff_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_z, pb_y, pb_z, dd_0, df_0, fp_s_3, ff_s_9, \
                         ff_s_10, ff_s_11, fp_3, fd_6, fd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_z[k] * df_0[k]
                 + f_2 * ff_s_9[k];

        t_10[k] = f_5 * dd_0[k]
                  + f_2 * ff_s_10[k]
                  + pb_z[k] * fd_6[k];

        t_11[k] = -f_6 * fp_s_3[k]
                  + f_2 * ff_s_11[k]
                  + f_5 * fp_3[k]
                  + pb_y[k] * fd_7[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_x, pf_s_3, pf_3, dd_6, dd_7, df_9, df_10, \
                         ff_s_12, ff_s_13, ff_s_14, fd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_4 * pf_s_3[k]
                  + f_5 * pf_3[k]
                  + pa_x[k] * df_9[k]
                  + f_2 * ff_s_12[k];

        t_13[k] = f_0 * dd_6[k]
                  + pa_x[k] * df_10[k]
                  + f_2 * ff_s_13[k];

        t_14[k] = f_5 * dd_7[k]
                  + f_2 * ff_s_14[k]
                  + pb_x[k] * fd_10[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pa_y, df_8, df_12, df_14, df_15, \
                         ff_s_15, ff_s_16, ff_s_17, ff_s_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_x[k] * df_12[k]
                  + f_2 * ff_s_15[k];

        t_16[k] = pa_x[k] * df_14[k]
                  + f_2 * ff_s_16[k];

        t_17[k] = pa_x[k] * df_15[k]
                  + f_2 * ff_s_17[k];

        t_18[k] = pa_y[k] * df_8[k]
                  + f_2 * ff_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_z, dd_4, dd_11, df_17, df_18, df_20, \
                         ff_s_19, ff_s_20, ff_s_21, ff_s_22, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_x[k] * df_17[k]
                  + f_2 * ff_s_19[k];

        t_20[k] = pa_x[k] * df_18[k]
                  + f_2 * ff_s_20[k];

        t_21[k] = f_0 * dd_11[k]
                  + pa_x[k] * df_20[k]
                  + f_2 * ff_s_21[k];

        t_22[k] = f_3 * dd_4[k]
                  + f_2 * ff_s_22[k]
                  + pb_z[k] * fd_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_x, dd_14, df_24, df_25, df_27, \
                         ff_s_23, ff_s_24, ff_s_25, ff_s_26, fd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * dd_14[k]
                  + f_2 * ff_s_23[k]
                  + pb_x[k] * fd_12[k];

        t_24[k] = pa_x[k] * df_24[k]
                  + f_2 * ff_s_24[k];

        t_25[k] = pa_x[k] * df_25[k]
                  + f_2 * ff_s_25[k];

        t_26[k] = pa_x[k] * df_27[k]
                  + f_2 * ff_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_x, fp_s_4, fp_s_5, ff_s_27, ff_s_28, ff_s_29, \
                         fp_4, fp_5, fd_13, fd_14, fd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_1 * fp_s_4[k]
                  + f_2 * ff_s_27[k]
                  + f_3 * fp_4[k]
                  + pb_x[k] * fd_13[k];

        t_28[k] = -f_6 * fp_s_5[k]
                  + f_2 * ff_s_28[k]
                  + f_5 * fp_5[k]
                  + pb_x[k] * fd_14[k];

        t_29[k] = f_2 * ff_s_29[k]
                  + pb_x[k] * fd_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_x, pb_y, pb_z, dd_7, fp_s_5, ff_s_30, ff_s_31, \
                         ff_s_32, fp_5, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * ff_s_30[k]
                  + pb_x[k] * fd_16[k];

        t_31[k] = f_0 * dd_7[k]
                  - f_1 * fp_s_5[k]
                  + f_2 * ff_s_31[k]
                  + f_3 * fp_5[k]
                  + pb_y[k] * fd_15[k];

        t_32[k] = f_2 * ff_s_32[k]
                  + pb_z[k] * fd_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_x, pb_y, pb_z, dd_8, fp_s_6, ff_s_33, ff_s_34, \
                         ff_s_35, fp_6, fd_16, fd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * dd_8[k]
                  + f_2 * ff_s_33[k]
                  + pb_y[k] * fd_16[k];

        t_34[k] = -f_1 * fp_s_6[k]
                  + f_2 * ff_s_34[k]
                  + f_3 * fp_6[k]
                  + pb_z[k] * fd_16[k];

        t_35[k] = f_2 * ff_s_35[k]
                  + pb_x[k] * fd_18[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_z, pb_y, pb_z, dd_7, dd_10, df_12, ff_s_36, \
                         ff_s_37, ff_s_38, fd_17, fd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_z[k] * df_12[k]
                  + f_2 * ff_s_36[k];

        t_37[k] = f_5 * dd_7[k]
                  + f_2 * ff_s_37[k]
                  + pb_z[k] * fd_17[k];

        t_38[k] = f_3 * dd_10[k]
                  + f_2 * ff_s_38[k]
                  + pb_y[k] * fd_18[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_x, pf_s_3, pf_3, dd_12, df_19, df_24, \
                         ff_s_39, ff_s_40, ff_s_41, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -f_4 * pf_s_3[k]
                  + f_5 * pf_3[k]
                  + pa_y[k] * df_19[k]
                  + f_2 * ff_s_39[k];

        t_40[k] = f_2 * ff_s_40[k]
                  + pb_x[k] * fd_19[k];

        t_41[k] = f_0 * dd_12[k]
                  + pa_y[k] * df_24[k]
                  + f_2 * ff_s_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_y, pb_y, pb_z, dd_9, dd_14, df_27, ff_s_42, \
                         ff_s_43, ff_s_44, fd_19, fd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * dd_9[k]
                  + f_2 * ff_s_42[k]
                  + pb_z[k] * fd_19[k];

        t_43[k] = f_5 * dd_14[k]
                  + f_2 * ff_s_43[k]
                  + pb_y[k] * fd_20[k];

        t_44[k] = pa_y[k] * df_27[k]
                  + f_2 * ff_s_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pb_x, pb_y, fp_s_9, fp_s_11, ff_s_45, ff_s_46, \
                         ff_s_47, fp_8, fp_10, fd_21, fd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -f_1 * fp_s_9[k]
                  + f_2 * ff_s_45[k]
                  + f_3 * fp_8[k]
                  + pb_x[k] * fd_21[k];

        t_46[k] = f_2 * ff_s_46[k]
                  + pb_y[k] * fd_21[k];

        t_47[k] = -f_6 * fp_s_11[k]
                  + f_2 * ff_s_47[k]
                  + f_5 * fp_10[k]
                  + pb_x[k] * fd_22[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_x, pb_y, fp_s_10, ff_s_48, ff_s_49, ff_s_50, \
                         fp_9, fd_23, fd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_2 * ff_s_48[k]
                  + pb_x[k] * fd_23[k];

        t_49[k] = f_2 * ff_s_49[k]
                  + pb_x[k] * fd_25[k];

        t_50[k] = -f_1 * fp_s_10[k]
                  + f_2 * ff_s_50[k]
                  + f_3 * fp_9[k]
                  + pb_y[k] * fd_23[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_y, pb_z, dd_14, fp_s_11, ff_s_51, ff_s_52, \
                         ff_s_53, fp_10, fd_24, fd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -f_6 * fp_s_11[k]
                  + f_2 * ff_s_51[k]
                  + f_5 * fp_10[k]
                  + pb_y[k] * fd_24[k];

        t_52[k] = f_2 * ff_s_52[k]
                  + pb_y[k] * fd_25[k];

        t_53[k] = f_0 * dd_14[k]
                  - f_1 * fp_s_11[k]
                  + f_2 * ff_s_53[k]
                  + f_3 * fp_10[k]
                  + pb_z[k] * fd_25[k];
    }
}

auto
compute_prim_ff_kinetic_energy_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf_s, const size_t pf,
                                 const size_t dd, const size_t df, const size_t fp_s,
                                 const size_t ff_s, const size_t fp, const size_t fd,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = alpha / p;

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

    const auto *pf_s_2 = buffer.data(pf_s + 2);
    const auto *pf_s_6 = buffer.data(pf_s + 6);

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

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_1 = buffer.data(fp_s + 1);
    const auto *fp_s_2 = buffer.data(fp_s + 2);
    const auto *fp_s_3 = buffer.data(fp_s + 3);
    const auto *fp_s_4 = buffer.data(fp_s + 4);
    const auto *fp_s_5 = buffer.data(fp_s + 5);
    const auto *fp_s_6 = buffer.data(fp_s + 6);
    const auto *fp_s_9 = buffer.data(fp_s + 9);
    const auto *fp_s_10 = buffer.data(fp_s + 10);
    const auto *fp_s_11 = buffer.data(fp_s + 11);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_30 = buffer.data(ff_s + 30);
    const auto *ff_s_31 = buffer.data(ff_s + 31);
    const auto *ff_s_32 = buffer.data(ff_s + 32);
    const auto *ff_s_33 = buffer.data(ff_s + 33);
    const auto *ff_s_34 = buffer.data(ff_s + 34);
    const auto *ff_s_35 = buffer.data(ff_s + 35);
    const auto *ff_s_36 = buffer.data(ff_s + 36);
    const auto *ff_s_37 = buffer.data(ff_s + 37);
    const auto *ff_s_38 = buffer.data(ff_s + 38);
    const auto *ff_s_39 = buffer.data(ff_s + 39);
    const auto *ff_s_40 = buffer.data(ff_s + 40);
    const auto *ff_s_41 = buffer.data(ff_s + 41);
    const auto *ff_s_42 = buffer.data(ff_s + 42);
    const auto *ff_s_43 = buffer.data(ff_s + 43);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_22 = buffer.data(fd + 22);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dd_0, fp_s_0, ff_s_0, ff_s_1, \
                         ff_s_2, fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = f_2 * ff_s_1[k]
                 + pb_y[k] * fd_0[k];

        t_2[k] = f_2 * ff_s_2[k]
                 + pb_z[k] * fd_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_y, pb_z, fp_s_1, fp_s_2, ff_s_3, ff_s_4, ff_s_5, \
                         fp_1, fp_2, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_1 * fp_s_1[k]
                 + f_2 * ff_s_3[k]
                 + f_3 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = f_2 * ff_s_4[k]
                 + pb_y[k] * fd_2[k];

        t_5[k] = -f_1 * fp_s_2[k]
                 + f_2 * ff_s_5[k]
                 + f_3 * fp_2[k]
                 + pb_z[k] * fd_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, pb_z, pf_s_2, pf_2, df_0, df_5, ff_s_6, \
                         ff_s_7, ff_s_8, fd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * df_0[k]
                 + f_2 * ff_s_6[k];

        t_7[k] = -f_4 * pf_s_2[k]
                 + f_5 * pf_2[k]
                 + pa_x[k] * df_5[k]
                 + f_2 * ff_s_7[k];

        t_8[k] = f_2 * ff_s_8[k]
                 + pb_z[k] * fd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pa_z, pb_y, df_0, df_4, fp_s_3, ff_s_9, \
                         ff_s_10, ff_s_11, fp_3, fd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_y[k] * df_4[k]
                 + f_2 * ff_s_9[k];

        t_10[k] = pa_z[k] * df_0[k]
                  + f_2 * ff_s_10[k];

        t_11[k] = -f_6 * fp_s_3[k]
                  + f_2 * ff_s_11[k]
                  + f_5 * fp_3[k]
                  + pb_y[k] * fd_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_y, pf_s_6, pf_6, dd_5, df_7, df_8, \
                         ff_s_12, ff_s_13, ff_s_14, fd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_2 * ff_s_12[k]
                  + pb_y[k] * fd_7[k];

        t_13[k] = -f_4 * pf_s_6[k]
                  + f_5 * pf_6[k]
                  + pa_x[k] * df_7[k]
                  + f_2 * ff_s_13[k];

        t_14[k] = f_0 * dd_5[k]
                  + pa_x[k] * df_8[k]
                  + f_2 * ff_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pb_x, dd_6, dd_9, df_11, df_16, ff_s_16, \
                         ff_s_17, ff_s_18, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * dd_6[k]
                  + f_2 * ff_s_16[k]
                  + pb_x[k] * fd_8[k];

        t_16[k] = pa_x[k] * df_11[k]
                  + f_2 * ff_s_17[k];

        t_17[k] = f_0 * dd_9[k]
                  + pa_x[k] * df_16[k]
                  + f_2 * ff_s_18[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pb_x, dd_12, df_23, fp_s_4, ff_s_19, ff_s_20, \
                         ff_s_21, fp_4, fd_9, fd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * dd_12[k]
                  + f_2 * ff_s_19[k]
                  + pb_x[k] * fd_9[k];

        t_19[k] = pa_x[k] * df_23[k]
                  + f_2 * ff_s_20[k];

        t_20[k] = -f_1 * fp_s_4[k]
                  + f_2 * ff_s_21[k]
                  + f_3 * fp_4[k]
                  + pb_x[k] * fd_10[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_x, pb_y, dd_6, fp_s_5, ff_s_22, ff_s_23, \
                         ff_s_24, ff_s_25, fp_5, fd_11, fd_12, fd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_6 * fp_s_5[k]
                  + f_2 * ff_s_22[k]
                  + f_5 * fp_5[k]
                  + pb_x[k] * fd_11[k];

        t_22[k] = f_2 * ff_s_23[k]
                  + pb_x[k] * fd_12[k];

        t_23[k] = f_2 * ff_s_24[k]
                  + pb_x[k] * fd_13[k];

        t_24[k] = f_0 * dd_6[k]
                  - f_1 * fp_s_5[k]
                  + f_2 * ff_s_25[k]
                  + f_3 * fp_5[k]
                  + pb_y[k] * fd_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pb_y, pb_z, dd_7, fp_s_6, ff_s_26, ff_s_27, \
                         ff_s_28, fp_6, fd_12, fd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_2 * ff_s_26[k]
                  + pb_z[k] * fd_12[k];

        t_26[k] = f_0 * dd_7[k]
                  + f_2 * ff_s_27[k]
                  + pb_y[k] * fd_13[k];

        t_27[k] = -f_1 * fp_s_6[k]
                  + f_2 * ff_s_28[k]
                  + f_3 * fp_6[k]
                  + pb_z[k] * fd_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pa_z, pb_x, pf_s_6, pf_6, df_11, df_15, \
                         ff_s_29, ff_s_30, ff_s_31, fd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_2 * ff_s_29[k]
                  + pb_x[k] * fd_15[k];

        t_29[k] = pa_z[k] * df_11[k]
                  + f_2 * ff_s_30[k];

        t_30[k] = -f_4 * pf_s_6[k]
                  + f_5 * pf_6[k]
                  + pa_y[k] * df_15[k]
                  + f_2 * ff_s_31[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pb_x, pb_y, dd_10, dd_12, df_20, ff_s_32, \
                         ff_s_33, ff_s_34, fd_16, fd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_2 * ff_s_32[k]
                  + pb_x[k] * fd_16[k];

        t_32[k] = f_0 * dd_10[k]
                  + pa_y[k] * df_20[k]
                  + f_2 * ff_s_33[k];

        t_33[k] = f_5 * dd_12[k]
                  + f_2 * ff_s_34[k]
                  + pb_y[k] * fd_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pb_x, df_23, fp_s_9, fp_s_11, ff_s_35, \
                         ff_s_36, ff_s_37, fp_8, fp_10, fd_18, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_y[k] * df_23[k]
                  + f_2 * ff_s_35[k];

        t_35[k] = -f_1 * fp_s_9[k]
                  + f_2 * ff_s_36[k]
                  + f_3 * fp_8[k]
                  + pb_x[k] * fd_18[k];

        t_36[k] = -f_6 * fp_s_11[k]
                  + f_2 * ff_s_37[k]
                  + f_5 * fp_10[k]
                  + pb_x[k] * fd_19[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, pb_y, fp_s_10, ff_s_38, ff_s_39, ff_s_40, \
                         fp_9, fd_20, fd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_2 * ff_s_38[k]
                  + pb_x[k] * fd_20[k];

        t_38[k] = f_2 * ff_s_39[k]
                  + pb_x[k] * fd_22[k];

        t_39[k] = -f_1 * fp_s_10[k]
                  + f_2 * ff_s_40[k]
                  + f_3 * fp_9[k]
                  + pb_y[k] * fd_20[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_y, pb_z, dd_12, fp_s_11, ff_s_41, ff_s_42, \
                         ff_s_43, fp_10, fd_21, fd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -f_6 * fp_s_11[k]
                  + f_2 * ff_s_41[k]
                  + f_5 * fp_10[k]
                  + pb_y[k] * fd_21[k];

        t_41[k] = f_2 * ff_s_42[k]
                  + pb_y[k] * fd_22[k];

        t_42[k] = f_0 * dd_12[k]
                  - f_1 * fp_s_11[k]
                  + f_2 * ff_s_43[k]
                  + f_3 * fp_10[k]
                  + pb_z[k] * fd_22[k];
    }
}

auto
compute_prim_ff_kinetic_energy_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf_s, const size_t pf,
                                 const size_t dd, const size_t df, const size_t fp_s,
                                 const size_t ff_s, const size_t fp, const size_t fd,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = alpha / p;

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

    const auto *pf_s_2 = buffer.data(pf_s + 2);
    const auto *pf_s_5 = buffer.data(pf_s + 5);

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_12 = buffer.data(dd + 12);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_21 = buffer.data(df + 21);

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_1 = buffer.data(fp_s + 1);
    const auto *fp_s_2 = buffer.data(fp_s + 2);
    const auto *fp_s_3 = buffer.data(fp_s + 3);
    const auto *fp_s_4 = buffer.data(fp_s + 4);
    const auto *fp_s_5 = buffer.data(fp_s + 5);
    const auto *fp_s_6 = buffer.data(fp_s + 6);
    const auto *fp_s_9 = buffer.data(fp_s + 9);
    const auto *fp_s_10 = buffer.data(fp_s + 10);
    const auto *fp_s_11 = buffer.data(fp_s + 11);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_15 = buffer.data(ff_s + 15);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_30 = buffer.data(ff_s + 30);
    const auto *ff_s_32 = buffer.data(ff_s + 32);
    const auto *ff_s_33 = buffer.data(ff_s + 33);
    const auto *ff_s_34 = buffer.data(ff_s + 34);
    const auto *ff_s_35 = buffer.data(ff_s + 35);
    const auto *ff_s_36 = buffer.data(ff_s + 36);
    const auto *ff_s_37 = buffer.data(ff_s + 37);
    const auto *ff_s_38 = buffer.data(ff_s + 38);
    const auto *ff_s_39 = buffer.data(ff_s + 39);
    const auto *ff_s_40 = buffer.data(ff_s + 40);
    const auto *ff_s_41 = buffer.data(ff_s + 41);
    const auto *ff_s_42 = buffer.data(ff_s + 42);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dd_0, fp_s_0, ff_s_0, ff_s_1, \
                         ff_s_2, fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = f_2 * ff_s_1[k]
                 + pb_y[k] * fd_0[k];

        t_2[k] = f_2 * ff_s_2[k]
                 + pb_z[k] * fd_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_y, pb_z, fp_s_1, fp_s_2, ff_s_3, ff_s_4, ff_s_5, \
                         fp_1, fp_2, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_1 * fp_s_1[k]
                 + f_2 * ff_s_3[k]
                 + f_3 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = f_2 * ff_s_4[k]
                 + pb_y[k] * fd_2[k];

        t_5[k] = -f_1 * fp_s_2[k]
                 + f_2 * ff_s_5[k]
                 + f_3 * fp_2[k]
                 + pb_z[k] * fd_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, pb_z, pf_s_2, pf_1, df_0, df_5, ff_s_6, \
                         ff_s_7, ff_s_8, fd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * df_0[k]
                 + f_2 * ff_s_6[k];

        t_7[k] = -f_4 * pf_s_2[k]
                 + f_5 * pf_1[k]
                 + pa_x[k] * df_5[k]
                 + f_2 * ff_s_7[k];

        t_8[k] = f_2 * ff_s_8[k]
                 + pb_z[k] * fd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_z, pb_y, df_0, fp_s_3, ff_s_9, ff_s_10, ff_s_11, \
                         fp_3, fd_6, fd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_z[k] * df_0[k]
                 + f_2 * ff_s_9[k];

        t_10[k] = -f_6 * fp_s_3[k]
                  + f_2 * ff_s_10[k]
                  + f_5 * fp_3[k]
                  + pb_y[k] * fd_6[k];

        t_11[k] = f_2 * ff_s_11[k]
                  + pb_y[k] * fd_7[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_x, pf_s_5, pf_2, dd_6, df_6, df_10, \
                         ff_s_12, ff_s_15, ff_s_16, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_4 * pf_s_5[k]
                  + f_5 * pf_2[k]
                  + pa_x[k] * df_6[k]
                  + f_2 * ff_s_12[k];

        t_13[k] = f_5 * dd_6[k]
                  + f_2 * ff_s_15[k]
                  + pb_x[k] * fd_8[k];

        t_14[k] = pa_x[k] * df_10[k]
                  + f_2 * ff_s_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pb_x, dd_12, df_21, fp_s_4, ff_s_18, ff_s_19, \
                         ff_s_20, fp_4, fd_9, fd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * dd_12[k]
                  + f_2 * ff_s_18[k]
                  + pb_x[k] * fd_9[k];

        t_16[k] = pa_x[k] * df_21[k]
                  + f_2 * ff_s_19[k];

        t_17[k] = -f_1 * fp_s_4[k]
                  + f_2 * ff_s_20[k]
                  + f_3 * fp_4[k]
                  + pb_x[k] * fd_10[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_x, pb_y, dd_6, fp_s_5, ff_s_21, ff_s_22, \
                         ff_s_23, ff_s_24, fp_5, fd_11, fd_12, fd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -f_6 * fp_s_5[k]
                  + f_2 * ff_s_21[k]
                  + f_5 * fp_5[k]
                  + pb_x[k] * fd_11[k];

        t_19[k] = f_2 * ff_s_22[k]
                  + pb_x[k] * fd_12[k];

        t_20[k] = f_2 * ff_s_23[k]
                  + pb_x[k] * fd_13[k];

        t_21[k] = f_0 * dd_6[k]
                  - f_1 * fp_s_5[k]
                  + f_2 * ff_s_24[k]
                  + f_3 * fp_5[k]
                  + pb_y[k] * fd_12[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pb_y, pb_z, dd_7, fp_s_6, ff_s_25, ff_s_26, \
                         ff_s_27, fp_6, fd_12, fd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_2 * ff_s_25[k]
                  + pb_z[k] * fd_12[k];

        t_23[k] = f_0 * dd_7[k]
                  + f_2 * ff_s_26[k]
                  + pb_y[k] * fd_13[k];

        t_24[k] = -f_1 * fp_s_6[k]
                  + f_2 * ff_s_27[k]
                  + f_3 * fp_6[k]
                  + pb_z[k] * fd_13[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pb_x, pf_s_5, pf_2, df_10, df_13, \
                         ff_s_28, ff_s_29, ff_s_30, fd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_2 * ff_s_28[k]
                  + pb_x[k] * fd_15[k];

        t_26[k] = pa_z[k] * df_10[k]
                  + f_2 * ff_s_29[k];

        t_27[k] = -f_4 * pf_s_5[k]
                  + f_5 * pf_2[k]
                  + pa_y[k] * df_13[k]
                  + f_2 * ff_s_30[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pb_y, dd_10, dd_12, df_18, df_21, ff_s_32, \
                         ff_s_33, ff_s_34, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * dd_10[k]
                  + pa_y[k] * df_18[k]
                  + f_2 * ff_s_32[k];

        t_29[k] = f_5 * dd_12[k]
                  + f_2 * ff_s_33[k]
                  + pb_y[k] * fd_16[k];

        t_30[k] = pa_y[k] * df_21[k]
                  + f_2 * ff_s_34[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pb_x, fp_s_9, fp_s_11, ff_s_35, ff_s_36, ff_s_37, \
                         fp_8, fp_10, fd_17, fd_18, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -f_1 * fp_s_9[k]
                  + f_2 * ff_s_35[k]
                  + f_3 * fp_8[k]
                  + pb_x[k] * fd_17[k];

        t_32[k] = -f_6 * fp_s_11[k]
                  + f_2 * ff_s_36[k]
                  + f_5 * fp_10[k]
                  + pb_x[k] * fd_18[k];

        t_33[k] = f_2 * ff_s_37[k]
                  + pb_x[k] * fd_19[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pb_x, pb_y, fp_s_10, fp_s_11, ff_s_38, ff_s_39, \
                         ff_s_40, fp_9, fp_10, fd_19, fd_20, fd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_2 * ff_s_38[k]
                  + pb_x[k] * fd_21[k];

        t_35[k] = -f_1 * fp_s_10[k]
                  + f_2 * ff_s_39[k]
                  + f_3 * fp_9[k]
                  + pb_y[k] * fd_19[k];

        t_36[k] = -f_6 * fp_s_11[k]
                  + f_2 * ff_s_40[k]
                  + f_5 * fp_10[k]
                  + pb_y[k] * fd_20[k];
    }

#pragma omp simd aligned(t_37, t_38, pb_y, pb_z, dd_12, fp_s_11, ff_s_41, ff_s_42, fp_10, \
                         fd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_2 * ff_s_41[k]
                  + pb_y[k] * fd_21[k];

        t_38[k] = f_0 * dd_12[k]
                  - f_1 * fp_s_11[k]
                  + f_2 * ff_s_42[k]
                  + f_3 * fp_10[k]
                  + pb_z[k] * fd_21[k];
    }
}

auto
compute_prim_ff_kinetic_energy_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf_s, const size_t pf,
                                 const size_t dd, const size_t df, const size_t fp_s,
                                 const size_t ff_s, const size_t fp, const size_t fd,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = beta / p;
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
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_s_0 = buffer.data(pf_s + 0);
    const auto *pf_s_1 = buffer.data(pf_s + 1);

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

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_4 = buffer.data(fp_s + 4);
    const auto *fp_s_5 = buffer.data(fp_s + 5);
    const auto *fp_s_9 = buffer.data(fp_s + 9);
    const auto *fp_s_10 = buffer.data(fp_s + 10);
    const auto *fp_s_11 = buffer.data(fp_s + 11);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_15 = buffer.data(ff_s + 15);

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

#pragma omp simd aligned(t_0, t_1, pa_y, pb_x, dd_0, df_0, fp_s_0, ff_s_0, ff_s_1, fp_0, \
                         fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pa_y[k] * df_0[k]
                 + f_2 * ff_s_1[k];
    }

#pragma omp simd aligned(t_2, t_3, t_4, pa_x, pa_z, pf_s_0, pf_s_1, pf_0, pf_1, df_0, df_1, \
                         df_2, ff_s_2, ff_s_3, ff_s_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = -f_4 * pf_s_0[k]
                 + f_5 * pf_0[k]
                 + pa_x[k] * df_1[k]
                 + f_2 * ff_s_2[k];

        t_3[k] = pa_z[k] * df_0[k]
                 + f_2 * ff_s_3[k];

        t_4[k] = -f_4 * pf_s_1[k]
                 + f_5 * pf_1[k]
                 + pa_x[k] * df_2[k]
                 + f_2 * ff_s_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_x, pb_x, df_3, df_6, fp_s_4, ff_s_5, ff_s_6, \
                         ff_s_7, fp_4, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_x[k] * df_3[k]
                 + f_2 * ff_s_5[k];

        t_6[k] = pa_x[k] * df_6[k]
                 + f_2 * ff_s_6[k];

        t_7[k] = -f_1 * fp_s_4[k]
                 + f_2 * ff_s_7[k]
                 + f_3 * fp_4[k]
                 + pb_x[k] * fd_11[k];
    }

#pragma omp simd aligned(t_8, t_9, pa_z, pb_y, dd_4, df_3, fp_s_5, ff_s_8, ff_s_9, fp_5, \
                         fd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dd_4[k]
                 - f_1 * fp_s_5[k]
                 + f_2 * ff_s_8[k]
                 + f_3 * fp_5[k]
                 + pb_y[k] * fd_12[k];

        t_9[k] = pa_z[k] * df_3[k]
                 + f_2 * ff_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pf_s_1, pf_1, dd_7, df_4, df_5, df_6, \
                         ff_s_10, ff_s_11, ff_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_4 * pf_s_1[k]
                  + f_5 * pf_1[k]
                  + pa_y[k] * df_4[k]
                  + f_2 * ff_s_10[k];

        t_11[k] = f_0 * dd_7[k]
                  + pa_y[k] * df_5[k]
                  + f_2 * ff_s_11[k];

        t_12[k] = pa_y[k] * df_6[k]
                  + f_2 * ff_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pb_x, pb_y, fp_s_9, fp_s_10, ff_s_13, ff_s_14, fp_9, \
                         fp_10, fd_20, fd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_1 * fp_s_9[k]
                  + f_2 * ff_s_13[k]
                  + f_3 * fp_9[k]
                  + pb_x[k] * fd_20[k];

        t_14[k] = -f_1 * fp_s_10[k]
                  + f_2 * ff_s_14[k]
                  + f_3 * fp_10[k]
                  + pb_y[k] * fd_21[k];
    }

#pragma omp simd aligned(t_15, pb_z, dd_8, fp_s_11, ff_s_15, fp_11, \
                         fd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * dd_8[k]
                  - f_1 * fp_s_11[k]
                  + f_2 * ff_s_15[k]
                  + f_3 * fp_11[k]
                  + pb_z[k] * fd_22[k];
    }
}

auto
compute_prim_ff_kinetic_energy_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pf_s, const size_t pf,
                                 const size_t dd, const size_t df, const size_t fp_s,
                                 const size_t ff_s, const size_t fp, const size_t fd,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = alpha / p;

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

    const auto *pf_s_1 = buffer.data(pf_s + 1);
    const auto *pf_s_3 = buffer.data(pf_s + 3);

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

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_1 = buffer.data(fp_s + 1);
    const auto *fp_s_2 = buffer.data(fp_s + 2);
    const auto *fp_s_3 = buffer.data(fp_s + 3);
    const auto *fp_s_4 = buffer.data(fp_s + 4);
    const auto *fp_s_5 = buffer.data(fp_s + 5);
    const auto *fp_s_6 = buffer.data(fp_s + 6);
    const auto *fp_s_7 = buffer.data(fp_s + 7);
    const auto *fp_s_8 = buffer.data(fp_s + 8);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_15 = buffer.data(ff_s + 15);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_30 = buffer.data(ff_s + 30);
    const auto *ff_s_31 = buffer.data(ff_s + 31);
    const auto *ff_s_32 = buffer.data(ff_s + 32);
    const auto *ff_s_33 = buffer.data(ff_s + 33);
    const auto *ff_s_34 = buffer.data(ff_s + 34);
    const auto *ff_s_35 = buffer.data(ff_s + 35);
    const auto *ff_s_36 = buffer.data(ff_s + 36);
    const auto *ff_s_37 = buffer.data(ff_s + 37);
    const auto *ff_s_38 = buffer.data(ff_s + 38);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dd_0, fp_s_0, ff_s_0, ff_s_1, \
                         ff_s_2, fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = f_2 * ff_s_1[k]
                 + pb_y[k] * fd_0[k];

        t_2[k] = f_2 * ff_s_2[k]
                 + pb_z[k] * fd_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_y, pb_y, pb_z, df_0, fp_s_1, fp_s_2, ff_s_3, \
                         ff_s_4, ff_s_5, fp_1, fp_2, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_1 * fp_s_1[k]
                 + f_2 * ff_s_3[k]
                 + f_3 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = -f_1 * fp_s_2[k]
                 + f_2 * ff_s_4[k]
                 + f_3 * fp_2[k]
                 + pb_z[k] * fd_2[k];

        t_5[k] = pa_y[k] * df_0[k]
                 + f_2 * ff_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_z, pf_s_1, pf_1, dd_0, df_0, df_2, \
                         ff_s_6, ff_s_7, ff_s_8, fd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * pf_s_1[k]
                 + f_5 * pf_1[k]
                 + pa_x[k] * df_2[k]
                 + f_2 * ff_s_6[k];

        t_7[k] = pa_z[k] * df_0[k]
                 + f_2 * ff_s_7[k];

        t_8[k] = f_5 * dd_0[k]
                 + f_2 * ff_s_8[k]
                 + pb_z[k] * fd_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_x, pf_s_3, pf_3, dd_5, dd_6, df_3, df_4, \
                         ff_s_9, ff_s_10, ff_s_11, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_4 * pf_s_3[k]
                 + f_5 * pf_3[k]
                 + pa_x[k] * df_3[k]
                 + f_2 * ff_s_9[k];

        t_10[k] = f_0 * dd_5[k]
                  + pa_x[k] * df_4[k]
                  + f_2 * ff_s_10[k];

        t_11[k] = f_5 * dd_6[k]
                  + f_2 * ff_s_11[k]
                  + pb_x[k] * fd_9[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_z, dd_3, dd_10, df_5, df_8, ff_s_12, \
                         ff_s_13, ff_s_14, fd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_x[k] * df_5[k]
                  + f_2 * ff_s_12[k];

        t_13[k] = f_0 * dd_10[k]
                  + pa_x[k] * df_8[k]
                  + f_2 * ff_s_13[k];

        t_14[k] = f_3 * dd_3[k]
                  + f_2 * ff_s_14[k]
                  + pb_z[k] * fd_10[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pb_x, dd_12, df_11, fp_s_3, ff_s_15, ff_s_16, \
                         ff_s_17, fp_3, fd_11, fd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * dd_12[k]
                  + f_2 * ff_s_15[k]
                  + pb_x[k] * fd_11[k];

        t_16[k] = pa_x[k] * df_11[k]
                  + f_2 * ff_s_16[k];

        t_17[k] = -f_1 * fp_s_3[k]
                  + f_2 * ff_s_17[k]
                  + f_3 * fp_3[k]
                  + pb_x[k] * fd_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_x, pb_y, pb_z, dd_6, fp_s_4, ff_s_18, ff_s_19, \
                         ff_s_20, fp_4, fd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * ff_s_18[k]
                  + pb_x[k] * fd_13[k];

        t_19[k] = f_0 * dd_6[k]
                  - f_1 * fp_s_4[k]
                  + f_2 * ff_s_19[k]
                  + f_3 * fp_4[k]
                  + pb_y[k] * fd_13[k];

        t_20[k] = f_2 * ff_s_20[k]
                  + pb_z[k] * fd_13[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_z, pb_y, pb_z, dd_7, df_5, fp_s_5, ff_s_21, \
                         ff_s_22, ff_s_23, fp_5, fd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * dd_7[k]
                  + f_2 * ff_s_21[k]
                  + pb_y[k] * fd_14[k];

        t_22[k] = -f_1 * fp_s_5[k]
                  + f_2 * ff_s_22[k]
                  + f_3 * fp_5[k]
                  + pb_z[k] * fd_14[k];

        t_23[k] = pa_z[k] * df_5[k]
                  + f_2 * ff_s_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_y, pb_y, pb_z, pf_s_3, pf_3, dd_6, dd_9, df_7, \
                         ff_s_24, ff_s_25, ff_s_26, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * dd_6[k]
                  + f_2 * ff_s_24[k]
                  + pb_z[k] * fd_15[k];

        t_25[k] = f_3 * dd_9[k]
                  + f_2 * ff_s_25[k]
                  + pb_y[k] * fd_16[k];

        t_26[k] = -f_4 * pf_s_3[k]
                  + f_5 * pf_3[k]
                  + pa_y[k] * df_7[k]
                  + f_2 * ff_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_y, pb_z, dd_8, dd_11, dd_12, df_9, \
                         ff_s_27, ff_s_28, ff_s_29, fd_17, fd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * dd_11[k]
                  + pa_y[k] * df_9[k]
                  + f_2 * ff_s_27[k];

        t_28[k] = f_3 * dd_8[k]
                  + f_2 * ff_s_28[k]
                  + pb_z[k] * fd_17[k];

        t_29[k] = f_5 * dd_12[k]
                  + f_2 * ff_s_29[k]
                  + pb_y[k] * fd_18[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_x, pb_y, df_11, fp_s_6, ff_s_30, \
                         ff_s_31, ff_s_32, ff_s_33, fp_6, fd_19, \
                         fd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * df_11[k]
                  + f_2 * ff_s_30[k];

        t_31[k] = -f_1 * fp_s_6[k]
                  + f_2 * ff_s_31[k]
                  + f_3 * fp_6[k]
                  + pb_x[k] * fd_19[k];

        t_32[k] = f_2 * ff_s_32[k]
                  + pb_y[k] * fd_19[k];

        t_33[k] = f_2 * ff_s_33[k]
                  + pb_x[k] * fd_20[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pb_x, pb_y, fp_s_7, fp_s_8, ff_s_34, ff_s_35, \
                         ff_s_36, fp_7, fp_8, fd_20, fd_21, fd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_2 * ff_s_34[k]
                  + pb_x[k] * fd_22[k];

        t_35[k] = -f_1 * fp_s_7[k]
                  + f_2 * ff_s_35[k]
                  + f_3 * fp_7[k]
                  + pb_y[k] * fd_20[k];

        t_36[k] = -f_6 * fp_s_8[k]
                  + f_2 * ff_s_36[k]
                  + f_5 * fp_8[k]
                  + pb_y[k] * fd_21[k];
    }

#pragma omp simd aligned(t_37, t_38, pb_y, pb_z, dd_12, fp_s_8, ff_s_37, ff_s_38, fp_8, \
                         fd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_2 * ff_s_37[k]
                  + pb_y[k] * fd_22[k];

        t_38[k] = f_0 * dd_12[k]
                  - f_1 * fp_s_8[k]
                  + f_2 * ff_s_38[k]
                  + f_3 * fp_8[k]
                  + pb_z[k] * fd_22[k];
    }
}

auto
compute_prim_ff_kinetic_energy_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pf_s, const size_t pf,
                                  const size_t dd, const size_t df, const size_t fp_s,
                                  const size_t ff_s, const size_t fp, const size_t fd,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = alpha / p;

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

    const auto *pf_s_1 = buffer.data(pf_s + 1);
    const auto *pf_s_3 = buffer.data(pf_s + 3);

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

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_1 = buffer.data(fp_s + 1);
    const auto *fp_s_2 = buffer.data(fp_s + 2);
    const auto *fp_s_3 = buffer.data(fp_s + 3);
    const auto *fp_s_4 = buffer.data(fp_s + 4);
    const auto *fp_s_5 = buffer.data(fp_s + 5);
    const auto *fp_s_6 = buffer.data(fp_s + 6);
    const auto *fp_s_7 = buffer.data(fp_s + 7);
    const auto *fp_s_8 = buffer.data(fp_s + 8);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_15 = buffer.data(ff_s + 15);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_30 = buffer.data(ff_s + 30);
    const auto *ff_s_31 = buffer.data(ff_s + 31);
    const auto *ff_s_32 = buffer.data(ff_s + 32);
    const auto *ff_s_33 = buffer.data(ff_s + 33);
    const auto *ff_s_34 = buffer.data(ff_s + 34);
    const auto *ff_s_35 = buffer.data(ff_s + 35);
    const auto *ff_s_36 = buffer.data(ff_s + 36);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dd_0, fp_s_0, ff_s_0, ff_s_1, \
                         ff_s_2, fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = f_2 * ff_s_1[k]
                 + pb_y[k] * fd_0[k];

        t_2[k] = f_2 * ff_s_2[k]
                 + pb_z[k] * fd_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_y, pb_y, pb_z, df_0, fp_s_1, fp_s_2, ff_s_3, \
                         ff_s_4, ff_s_5, fp_1, fp_2, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_1 * fp_s_1[k]
                 + f_2 * ff_s_3[k]
                 + f_3 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = -f_1 * fp_s_2[k]
                 + f_2 * ff_s_4[k]
                 + f_3 * fp_2[k]
                 + pb_z[k] * fd_2[k];

        t_5[k] = pa_y[k] * df_0[k]
                 + f_2 * ff_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, pa_z, pf_s_1, pf_1, df_0, df_3, df_4, \
                         ff_s_6, ff_s_7, ff_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * pf_s_1[k]
                 + f_5 * pf_1[k]
                 + pa_x[k] * df_4[k]
                 + f_2 * ff_s_6[k];

        t_7[k] = pa_y[k] * df_3[k]
                 + f_2 * ff_s_7[k];

        t_8[k] = pa_z[k] * df_0[k]
                 + f_2 * ff_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_x, pf_s_3, pf_3, dd_4, dd_5, df_6, df_7, \
                         ff_s_9, ff_s_10, ff_s_11, fd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_4 * pf_s_3[k]
                 + f_5 * pf_3[k]
                 + pa_x[k] * df_6[k]
                 + f_2 * ff_s_9[k];

        t_10[k] = f_0 * dd_4[k]
                  + pa_x[k] * df_7[k]
                  + f_2 * ff_s_10[k];

        t_11[k] = f_5 * dd_5[k]
                  + f_2 * ff_s_11[k]
                  + pb_x[k] * fd_7[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, dd_8, dd_10, df_9, df_14, df_20, \
                         ff_s_12, ff_s_13, ff_s_14, ff_s_15, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_x[k] * df_9[k]
                  + f_2 * ff_s_12[k];

        t_13[k] = f_0 * dd_8[k]
                  + pa_x[k] * df_14[k]
                  + f_2 * ff_s_13[k];

        t_14[k] = f_5 * dd_10[k]
                  + f_2 * ff_s_14[k]
                  + pb_x[k] * fd_8[k];

        t_15[k] = pa_x[k] * df_20[k]
                  + f_2 * ff_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, fp_s_3, ff_s_16, ff_s_17, ff_s_18, fp_3, \
                         fd_9, fd_10, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_1 * fp_s_3[k]
                  + f_2 * ff_s_16[k]
                  + f_3 * fp_3[k]
                  + pb_x[k] * fd_9[k];

        t_17[k] = f_2 * ff_s_17[k]
                  + pb_x[k] * fd_10[k];

        t_18[k] = f_2 * ff_s_18[k]
                  + pb_x[k] * fd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_y, pb_z, dd_5, fp_s_4, fp_s_5, ff_s_19, ff_s_20, \
                         ff_s_21, fp_4, fp_5, fd_10, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * dd_5[k]
                  - f_1 * fp_s_4[k]
                  + f_2 * ff_s_19[k]
                  + f_3 * fp_4[k]
                  + pb_y[k] * fd_10[k];

        t_20[k] = f_2 * ff_s_20[k]
                  + pb_z[k] * fd_10[k];

        t_21[k] = -f_1 * fp_s_5[k]
                  + f_2 * ff_s_21[k]
                  + f_3 * fp_5[k]
                  + pb_z[k] * fd_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pa_z, pb_x, pf_s_3, pf_3, df_9, df_13, \
                         ff_s_22, ff_s_23, ff_s_24, fd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_2 * ff_s_22[k]
                  + pb_x[k] * fd_13[k];

        t_23[k] = pa_z[k] * df_9[k]
                  + f_2 * ff_s_23[k];

        t_24[k] = -f_4 * pf_s_3[k]
                  + f_5 * pf_3[k]
                  + pa_y[k] * df_13[k]
                  + f_2 * ff_s_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pb_x, pb_y, dd_9, dd_10, df_18, ff_s_25, \
                         ff_s_26, ff_s_27, fd_14, fd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_2 * ff_s_25[k]
                  + pb_x[k] * fd_14[k];

        t_26[k] = f_0 * dd_9[k]
                  + pa_y[k] * df_18[k]
                  + f_2 * ff_s_26[k];

        t_27[k] = f_5 * dd_10[k]
                  + f_2 * ff_s_27[k]
                  + pb_y[k] * fd_15[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pb_x, pb_y, df_20, fp_s_6, ff_s_28, \
                         ff_s_29, ff_s_30, ff_s_31, fp_6, fd_16, \
                         fd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_y[k] * df_20[k]
                  + f_2 * ff_s_28[k];

        t_29[k] = -f_1 * fp_s_6[k]
                  + f_2 * ff_s_29[k]
                  + f_3 * fp_6[k]
                  + pb_x[k] * fd_16[k];

        t_30[k] = f_2 * ff_s_30[k]
                  + pb_y[k] * fd_16[k];

        t_31[k] = f_2 * ff_s_31[k]
                  + pb_x[k] * fd_17[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, fp_s_7, fp_s_8, ff_s_32, ff_s_33, \
                         ff_s_34, fp_7, fp_8, fd_17, fd_18, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_2 * ff_s_32[k]
                  + pb_x[k] * fd_19[k];

        t_33[k] = -f_1 * fp_s_7[k]
                  + f_2 * ff_s_33[k]
                  + f_3 * fp_7[k]
                  + pb_y[k] * fd_17[k];

        t_34[k] = -f_6 * fp_s_8[k]
                  + f_2 * ff_s_34[k]
                  + f_5 * fp_8[k]
                  + pb_y[k] * fd_18[k];
    }

#pragma omp simd aligned(t_35, t_36, pb_y, pb_z, dd_10, fp_s_8, ff_s_35, ff_s_36, fp_8, \
                         fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_2 * ff_s_35[k]
                  + pb_y[k] * fd_19[k];

        t_36[k] = f_0 * dd_10[k]
                  - f_1 * fp_s_8[k]
                  + f_2 * ff_s_36[k]
                  + f_3 * fp_8[k]
                  + pb_z[k] * fd_19[k];
    }
}

auto
compute_prim_ff_kinetic_energy_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pf_s, const size_t pf,
                                  const size_t dd, const size_t df, const size_t fp_s,
                                  const size_t ff_s, const size_t fp, const size_t fd,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = alpha / p;

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

    const auto *pf_s_1 = buffer.data(pf_s + 1);
    const auto *pf_s_3 = buffer.data(pf_s + 3);

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

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

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_1 = buffer.data(fp_s + 1);
    const auto *fp_s_2 = buffer.data(fp_s + 2);
    const auto *fp_s_3 = buffer.data(fp_s + 3);
    const auto *fp_s_4 = buffer.data(fp_s + 4);
    const auto *fp_s_5 = buffer.data(fp_s + 5);
    const auto *fp_s_6 = buffer.data(fp_s + 6);
    const auto *fp_s_7 = buffer.data(fp_s + 7);
    const auto *fp_s_8 = buffer.data(fp_s + 8);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_15 = buffer.data(ff_s + 15);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_30 = buffer.data(ff_s + 30);
    const auto *ff_s_31 = buffer.data(ff_s + 31);
    const auto *ff_s_32 = buffer.data(ff_s + 32);
    const auto *ff_s_33 = buffer.data(ff_s + 33);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dd_0, fp_s_0, ff_s_0, ff_s_1, \
                         ff_s_2, fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = f_2 * ff_s_1[k]
                 + pb_y[k] * fd_0[k];

        t_2[k] = f_2 * ff_s_2[k]
                 + pb_z[k] * fd_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_y, pb_y, pb_z, df_0, fp_s_1, fp_s_2, ff_s_3, \
                         ff_s_4, ff_s_5, fp_1, fp_2, fd_1, fd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_1 * fp_s_1[k]
                 + f_2 * ff_s_3[k]
                 + f_3 * fp_1[k]
                 + pb_y[k] * fd_1[k];

        t_4[k] = -f_1 * fp_s_2[k]
                 + f_2 * ff_s_4[k]
                 + f_3 * fp_2[k]
                 + pb_z[k] * fd_2[k];

        t_5[k] = pa_y[k] * df_0[k]
                 + f_2 * ff_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pf_s_1, pf_s_3, pf_1, pf_2, df_0, df_4, \
                         df_5, ff_s_6, ff_s_7, ff_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * pf_s_1[k]
                 + f_5 * pf_1[k]
                 + pa_x[k] * df_4[k]
                 + f_2 * ff_s_6[k];

        t_7[k] = pa_z[k] * df_0[k]
                 + f_2 * ff_s_7[k];

        t_8[k] = -f_4 * pf_s_3[k]
                 + f_5 * pf_2[k]
                 + pa_x[k] * df_5[k]
                 + f_2 * ff_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pb_x, dd_5, dd_10, df_8, df_18, ff_s_9, \
                         ff_s_10, ff_s_11, ff_s_12, fd_7, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * dd_5[k]
                 + f_2 * ff_s_9[k]
                 + pb_x[k] * fd_7[k];

        t_10[k] = pa_x[k] * df_8[k]
                  + f_2 * ff_s_10[k];

        t_11[k] = f_5 * dd_10[k]
                  + f_2 * ff_s_11[k]
                  + pb_x[k] * fd_8[k];

        t_12[k] = pa_x[k] * df_18[k]
                  + f_2 * ff_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, fp_s_3, ff_s_13, ff_s_14, ff_s_15, fp_3, \
                         fd_9, fd_10, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_1 * fp_s_3[k]
                  + f_2 * ff_s_13[k]
                  + f_3 * fp_3[k]
                  + pb_x[k] * fd_9[k];

        t_14[k] = f_2 * ff_s_14[k]
                  + pb_x[k] * fd_10[k];

        t_15[k] = f_2 * ff_s_15[k]
                  + pb_x[k] * fd_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_y, pb_z, dd_5, fp_s_4, fp_s_5, ff_s_16, ff_s_17, \
                         ff_s_18, fp_4, fp_5, fd_10, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * dd_5[k]
                  - f_1 * fp_s_4[k]
                  + f_2 * ff_s_16[k]
                  + f_3 * fp_4[k]
                  + pb_y[k] * fd_10[k];

        t_17[k] = f_2 * ff_s_17[k]
                  + pb_z[k] * fd_10[k];

        t_18[k] = -f_1 * fp_s_5[k]
                  + f_2 * ff_s_18[k]
                  + f_3 * fp_5[k]
                  + pb_z[k] * fd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pa_z, pb_x, pf_s_3, pf_2, df_8, df_11, \
                         ff_s_19, ff_s_20, ff_s_21, fd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * ff_s_19[k]
                  + pb_x[k] * fd_13[k];

        t_20[k] = pa_z[k] * df_8[k]
                  + f_2 * ff_s_20[k];

        t_21[k] = -f_4 * pf_s_3[k]
                  + f_5 * pf_2[k]
                  + pa_y[k] * df_11[k]
                  + f_2 * ff_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pb_y, dd_9, dd_10, df_16, df_18, ff_s_23, \
                         ff_s_24, ff_s_25, fd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_0 * dd_9[k]
                  + pa_y[k] * df_16[k]
                  + f_2 * ff_s_23[k];

        t_23[k] = f_5 * dd_10[k]
                  + f_2 * ff_s_24[k]
                  + pb_y[k] * fd_14[k];

        t_24[k] = pa_y[k] * df_18[k]
                  + f_2 * ff_s_25[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pb_x, pb_y, fp_s_6, ff_s_26, ff_s_27, \
                         ff_s_28, ff_s_29, fp_6, fd_15, fd_16, fd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -f_1 * fp_s_6[k]
                  + f_2 * ff_s_26[k]
                  + f_3 * fp_6[k]
                  + pb_x[k] * fd_15[k];

        t_26[k] = f_2 * ff_s_27[k]
                  + pb_y[k] * fd_15[k];

        t_27[k] = f_2 * ff_s_28[k]
                  + pb_x[k] * fd_16[k];

        t_28[k] = f_2 * ff_s_29[k]
                  + pb_x[k] * fd_18[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_y, fp_s_7, fp_s_8, ff_s_30, ff_s_31, ff_s_32, \
                         fp_7, fp_8, fd_16, fd_17, fd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -f_1 * fp_s_7[k]
                  + f_2 * ff_s_30[k]
                  + f_3 * fp_7[k]
                  + pb_y[k] * fd_16[k];

        t_30[k] = -f_6 * fp_s_8[k]
                  + f_2 * ff_s_31[k]
                  + f_5 * fp_8[k]
                  + pb_y[k] * fd_17[k];

        t_31[k] = f_2 * ff_s_32[k]
                  + pb_y[k] * fd_18[k];
    }

#pragma omp simd aligned(t_32, pb_z, dd_10, fp_s_8, ff_s_33, fp_8, \
                         fd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * dd_10[k]
                  - f_1 * fp_s_8[k]
                  + f_2 * ff_s_33[k]
                  + f_3 * fp_8[k]
                  + pb_z[k] * fd_18[k];
    }
}

auto
compute_prim_ff_kinetic_energy_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pf_s, const size_t pf,
                                  const size_t dd, const size_t df, const size_t fp_s,
                                  const size_t ff_s, const size_t fp, const size_t fd,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = beta / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_s_0 = buffer.data(pf_s + 0);
    const auto *pf_s_1 = buffer.data(pf_s + 1);

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

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_3 = buffer.data(fp_s + 3);
    const auto *fp_s_7 = buffer.data(fp_s + 7);
    const auto *fp_s_8 = buffer.data(fp_s + 8);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);

#pragma omp simd aligned(t_0, t_1, pa_y, pb_x, dd_0, df_0, fp_s_0, ff_s_0, ff_s_1, fp_0, \
                         fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = pa_y[k] * df_0[k]
                 + f_2 * ff_s_1[k];
    }

#pragma omp simd aligned(t_2, t_3, t_4, pa_x, pa_z, pf_s_0, pf_s_1, pf_0, pf_1, df_0, df_1, \
                         df_2, ff_s_2, ff_s_3, ff_s_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = -f_4 * pf_s_0[k]
                 + f_5 * pf_0[k]
                 + pa_x[k] * df_1[k]
                 + f_2 * ff_s_2[k];

        t_3[k] = pa_z[k] * df_0[k]
                 + f_2 * ff_s_3[k];

        t_4[k] = -f_4 * pf_s_1[k]
                 + f_5 * pf_1[k]
                 + pa_x[k] * df_2[k]
                 + f_2 * ff_s_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_x, pb_y, dd_3, df_3, df_6, fp_s_3, ff_s_5, ff_s_6, \
                         ff_s_7, fp_3, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_x[k] * df_3[k]
                 + f_2 * ff_s_5[k];

        t_6[k] = pa_x[k] * df_6[k]
                 + f_2 * ff_s_6[k];

        t_7[k] = f_0 * dd_3[k]
                 - f_1 * fp_s_3[k]
                 + f_2 * ff_s_7[k]
                 + f_3 * fp_3[k]
                 + pb_y[k] * fd_8[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_y, pa_z, pf_s_1, pf_1, dd_5, df_3, df_4, df_5, \
                         ff_s_8, ff_s_9, ff_s_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = pa_z[k] * df_3[k]
                 + f_2 * ff_s_8[k];

        t_9[k] = -f_4 * pf_s_1[k]
                 + f_5 * pf_1[k]
                 + pa_y[k] * df_4[k]
                 + f_2 * ff_s_9[k];

        t_10[k] = f_0 * dd_5[k]
                  + pa_y[k] * df_5[k]
                  + f_2 * ff_s_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pa_y, pb_y, df_6, fp_s_7, ff_s_11, ff_s_12, fp_7, \
                         fd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_y[k] * df_6[k]
                  + f_2 * ff_s_11[k];

        t_12[k] = -f_1 * fp_s_7[k]
                  + f_2 * ff_s_12[k]
                  + f_3 * fp_7[k]
                  + pb_y[k] * fd_14[k];
    }

#pragma omp simd aligned(t_13, pb_z, dd_6, fp_s_8, ff_s_13, fp_8, \
                         fd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dd_6[k]
                  - f_1 * fp_s_8[k]
                  + f_2 * ff_s_13[k]
                  + f_3 * fp_8[k]
                  + pb_z[k] * fd_15[k];
    }
}

auto
compute_prim_ff_kinetic_energy_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pf_s, const size_t pf,
                                  const size_t dd, const size_t df, const size_t fp_s,
                                  const size_t ff_s, const size_t fp, const size_t fd,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = beta / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_s_1 = buffer.data(pf_s + 1);
    const auto *pf_s_3 = buffer.data(pf_s + 3);

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

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_2 = buffer.data(fp_s + 2);
    const auto *fp_s_3 = buffer.data(fp_s + 3);
    const auto *fp_s_5 = buffer.data(fp_s + 5);
    const auto *fp_s_6 = buffer.data(fp_s + 6);
    const auto *fp_s_7 = buffer.data(fp_s + 7);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_15 = buffer.data(ff_s + 15);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_21 = buffer.data(ff_s + 21);

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

#pragma omp simd aligned(t_0, t_1, t_2, pa_y, pb_x, pb_z, dd_0, df_0, fp_s_0, ff_s_0, ff_s_1, \
                         ff_s_2, fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = f_2 * ff_s_1[k]
                 + pb_z[k] * fd_0[k];

        t_2[k] = pa_y[k] * df_0[k]
                 + f_2 * ff_s_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_z, pf_s_1, pf_s_3, pf_1, pf_3, df_0, df_2, \
                         df_3, ff_s_3, ff_s_4, ff_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_4 * pf_s_1[k]
                 + f_5 * pf_1[k]
                 + pa_x[k] * df_2[k]
                 + f_2 * ff_s_3[k];

        t_4[k] = pa_z[k] * df_0[k]
                 + f_2 * ff_s_4[k];

        t_5[k] = -f_4 * pf_s_3[k]
                 + f_5 * pf_3[k]
                 + pa_x[k] * df_3[k]
                 + f_2 * ff_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, pb_x, dd_4, dd_8, df_4, df_9, ff_s_6, \
                         ff_s_7, ff_s_8, ff_s_9, fd_7, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * dd_4[k]
                 + f_2 * ff_s_6[k]
                 + pb_x[k] * fd_7[k];

        t_7[k] = pa_x[k] * df_4[k]
                 + f_2 * ff_s_7[k];

        t_8[k] = f_5 * dd_8[k]
                 + f_2 * ff_s_8[k]
                 + pb_x[k] * fd_9[k];

        t_9[k] = pa_x[k] * df_9[k]
                 + f_2 * ff_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, pb_y, pb_z, dd_4, fp_s_2, fp_s_3, ff_s_10, \
                         ff_s_11, ff_s_12, fp_2, fp_3, fd_10, fd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * fp_s_2[k]
                  + f_2 * ff_s_10[k]
                  + f_3 * fp_2[k]
                  + pb_x[k] * fd_10[k];

        t_11[k] = f_0 * dd_4[k]
                  - f_1 * fp_s_3[k]
                  + f_2 * ff_s_11[k]
                  + f_3 * fp_3[k]
                  + pb_y[k] * fd_11[k];

        t_12[k] = f_2 * ff_s_12[k]
                  + pb_z[k] * fd_11[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pa_z, pf_s_3, pf_3, dd_7, df_4, df_6, df_7, \
                         ff_s_13, ff_s_14, ff_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_z[k] * df_4[k]
                  + f_2 * ff_s_13[k];

        t_14[k] = -f_4 * pf_s_3[k]
                  + f_5 * pf_3[k]
                  + pa_y[k] * df_6[k]
                  + f_2 * ff_s_14[k];

        t_15[k] = f_0 * dd_7[k]
                  + pa_y[k] * df_7[k]
                  + f_2 * ff_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_y, pb_x, pb_y, dd_8, df_9, fp_s_5, ff_s_16, \
                         ff_s_17, ff_s_18, fp_5, fd_16, fd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * dd_8[k]
                  + f_2 * ff_s_16[k]
                  + pb_y[k] * fd_16[k];

        t_17[k] = pa_y[k] * df_9[k]
                  + f_2 * ff_s_17[k];

        t_18[k] = -f_1 * fp_s_5[k]
                  + f_2 * ff_s_18[k]
                  + f_3 * fp_5[k]
                  + pb_x[k] * fd_17[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_y, pb_z, dd_8, fp_s_6, fp_s_7, ff_s_19, ff_s_20, \
                         ff_s_21, fp_6, fp_7, fd_18, fd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_1 * fp_s_6[k]
                  + f_2 * ff_s_19[k]
                  + f_3 * fp_6[k]
                  + pb_y[k] * fd_18[k];

        t_20[k] = f_2 * ff_s_20[k]
                  + pb_y[k] * fd_19[k];

        t_21[k] = f_0 * dd_8[k]
                  - f_1 * fp_s_7[k]
                  + f_2 * ff_s_21[k]
                  + f_3 * fp_7[k]
                  + pb_z[k] * fd_19[k];
    }
}

auto
compute_prim_ff_kinetic_energy_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pf_s, const size_t pf,
                                  const size_t dd, const size_t df, const size_t fp_s,
                                  const size_t ff_s, const size_t fp, const size_t fd,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = beta / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pf_s_1 = buffer.data(pf_s + 1);
    const auto *pf_s_3 = buffer.data(pf_s + 3);

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

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_1 = buffer.data(fp_s + 1);
    const auto *fp_s_2 = buffer.data(fp_s + 2);
    const auto *fp_s_3 = buffer.data(fp_s + 3);
    const auto *fp_s_4 = buffer.data(fp_s + 4);
    const auto *fp_s_5 = buffer.data(fp_s + 5);
    const auto *fp_s_6 = buffer.data(fp_s + 6);
    const auto *fp_s_7 = buffer.data(fp_s + 7);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_15 = buffer.data(ff_s + 15);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_30 = buffer.data(ff_s + 30);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dd_0, fp_s_0, ff_s_0, ff_s_1, \
                         ff_s_2, fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = f_2 * ff_s_1[k]
                 + pb_y[k] * fd_0[k];

        t_2[k] = f_2 * ff_s_2[k]
                 + pb_z[k] * fd_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, pf_s_1, pf_1, df_0, df_2, fp_s_1, \
                         ff_s_3, ff_s_4, ff_s_5, fp_1, fd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_1 * fp_s_1[k]
                 + f_2 * ff_s_3[k]
                 + f_3 * fp_1[k]
                 + pb_z[k] * fd_1[k];

        t_4[k] = pa_y[k] * df_0[k]
                 + f_2 * ff_s_4[k];

        t_5[k] = -f_4 * pf_s_1[k]
                 + f_5 * pf_1[k]
                 + pa_x[k] * df_2[k]
                 + f_2 * ff_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pf_s_3, pf_3, dd_3, df_0, df_3, df_4, \
                         ff_s_6, ff_s_7, ff_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_z[k] * df_0[k]
                 + f_2 * ff_s_6[k];

        t_7[k] = -f_4 * pf_s_3[k]
                 + f_5 * pf_3[k]
                 + pa_x[k] * df_3[k]
                 + f_2 * ff_s_7[k];

        t_8[k] = f_0 * dd_3[k]
                 + pa_x[k] * df_4[k]
                 + f_2 * ff_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_x, dd_4, dd_6, df_6, df_9, ff_s_9, ff_s_10, \
                         ff_s_11, fd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * dd_4[k]
                 + f_2 * ff_s_9[k]
                 + pb_x[k] * fd_6[k];

        t_10[k] = pa_x[k] * df_6[k]
                  + f_2 * ff_s_10[k];

        t_11[k] = f_0 * dd_6[k]
                  + pa_x[k] * df_9[k]
                  + f_2 * ff_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_x, dd_8, df_13, fp_s_2, ff_s_12, ff_s_13, \
                         ff_s_14, fp_2, fd_7, fd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * dd_8[k]
                  + f_2 * ff_s_12[k]
                  + pb_x[k] * fd_7[k];

        t_13[k] = pa_x[k] * df_13[k]
                  + f_2 * ff_s_13[k];

        t_14[k] = -f_1 * fp_s_2[k]
                  + f_2 * ff_s_14[k]
                  + f_3 * fp_2[k]
                  + pb_x[k] * fd_8[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_x, pb_y, pb_z, dd_4, fp_s_3, ff_s_15, ff_s_16, \
                         ff_s_17, fp_3, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_2 * ff_s_15[k]
                  + pb_x[k] * fd_9[k];

        t_16[k] = f_0 * dd_4[k]
                  - f_1 * fp_s_3[k]
                  + f_2 * ff_s_16[k]
                  + f_3 * fp_3[k]
                  + pb_y[k] * fd_9[k];

        t_17[k] = f_2 * ff_s_17[k]
                  + pb_z[k] * fd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pa_z, pb_z, pf_s_3, pf_3, df_6, df_8, fp_s_4, \
                         ff_s_18, ff_s_19, ff_s_20, fp_4, fd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -f_1 * fp_s_4[k]
                  + f_2 * ff_s_18[k]
                  + f_3 * fp_4[k]
                  + pb_z[k] * fd_10[k];

        t_19[k] = pa_z[k] * df_6[k]
                  + f_2 * ff_s_19[k];

        t_20[k] = -f_4 * pf_s_3[k]
                  + f_5 * pf_3[k]
                  + pa_y[k] * df_8[k]
                  + f_2 * ff_s_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_y, pb_y, dd_7, dd_8, df_11, df_13, ff_s_21, \
                         ff_s_22, ff_s_23, fd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * dd_7[k]
                  + pa_y[k] * df_11[k]
                  + f_2 * ff_s_21[k];

        t_22[k] = f_5 * dd_8[k]
                  + f_2 * ff_s_22[k]
                  + pb_y[k] * fd_13[k];

        t_23[k] = pa_y[k] * df_13[k]
                  + f_2 * ff_s_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_x, pb_y, fp_s_5, ff_s_24, ff_s_25, \
                         ff_s_26, ff_s_27, fp_5, fd_14, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_1 * fp_s_5[k]
                  + f_2 * ff_s_24[k]
                  + f_3 * fp_5[k]
                  + pb_x[k] * fd_14[k];

        t_25[k] = f_2 * ff_s_25[k]
                  + pb_y[k] * fd_14[k];

        t_26[k] = f_2 * ff_s_26[k]
                  + pb_x[k] * fd_15[k];

        t_27[k] = f_2 * ff_s_27[k]
                  + pb_x[k] * fd_16[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pb_y, pb_z, dd_8, fp_s_6, fp_s_7, ff_s_28, ff_s_29, \
                         ff_s_30, fp_6, fp_7, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -f_1 * fp_s_6[k]
                  + f_2 * ff_s_28[k]
                  + f_3 * fp_6[k]
                  + pb_y[k] * fd_15[k];

        t_29[k] = f_2 * ff_s_29[k]
                  + pb_y[k] * fd_16[k];

        t_30[k] = f_0 * dd_8[k]
                  - f_1 * fp_s_7[k]
                  + f_2 * ff_s_30[k]
                  + f_3 * fp_7[k]
                  + pb_z[k] * fd_16[k];
    }
}

auto
compute_prim_ff_kinetic_energy_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pf_s, const size_t pf,
                                  const size_t dd, const size_t df, const size_t fp_s,
                                  const size_t ff_s, const size_t fp, const size_t fd,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = beta / p;
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

    const auto *pf_s_1 = buffer.data(pf_s + 1);
    const auto *pf_s_3 = buffer.data(pf_s + 3);

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_13 = buffer.data(df + 13);

    const auto *fp_s_0 = buffer.data(fp_s + 0);
    const auto *fp_s_1 = buffer.data(fp_s + 1);
    const auto *fp_s_2 = buffer.data(fp_s + 2);
    const auto *fp_s_3 = buffer.data(fp_s + 3);
    const auto *fp_s_4 = buffer.data(fp_s + 4);
    const auto *fp_s_5 = buffer.data(fp_s + 5);
    const auto *fp_s_6 = buffer.data(fp_s + 6);
    const auto *fp_s_7 = buffer.data(fp_s + 7);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_3 = buffer.data(ff_s + 3);
    const auto *ff_s_4 = buffer.data(ff_s + 4);
    const auto *ff_s_5 = buffer.data(ff_s + 5);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_7 = buffer.data(ff_s + 7);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_11 = buffer.data(ff_s + 11);
    const auto *ff_s_12 = buffer.data(ff_s + 12);
    const auto *ff_s_13 = buffer.data(ff_s + 13);
    const auto *ff_s_14 = buffer.data(ff_s + 14);
    const auto *ff_s_15 = buffer.data(ff_s + 15);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_18 = buffer.data(ff_s + 18);
    const auto *ff_s_19 = buffer.data(ff_s + 19);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_21 = buffer.data(ff_s + 21);
    const auto *ff_s_22 = buffer.data(ff_s + 22);
    const auto *ff_s_23 = buffer.data(ff_s + 23);
    const auto *ff_s_24 = buffer.data(ff_s + 24);
    const auto *ff_s_25 = buffer.data(ff_s + 25);
    const auto *ff_s_26 = buffer.data(ff_s + 26);
    const auto *ff_s_27 = buffer.data(ff_s + 27);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dd_0, fp_s_0, ff_s_0, ff_s_1, \
                         ff_s_2, fp_0, fd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dd_0[k]
                 - f_1 * fp_s_0[k]
                 + f_2 * ff_s_0[k]
                 + f_3 * fp_0[k]
                 + pb_x[k] * fd_0[k];

        t_1[k] = f_2 * ff_s_1[k]
                 + pb_y[k] * fd_0[k];

        t_2[k] = f_2 * ff_s_2[k]
                 + pb_z[k] * fd_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, pf_s_1, pf_1, df_0, df_2, fp_s_1, \
                         ff_s_3, ff_s_4, ff_s_5, fp_1, fd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_1 * fp_s_1[k]
                 + f_2 * ff_s_3[k]
                 + f_3 * fp_1[k]
                 + pb_z[k] * fd_1[k];

        t_4[k] = pa_y[k] * df_0[k]
                 + f_2 * ff_s_4[k];

        t_5[k] = -f_4 * pf_s_1[k]
                 + f_5 * pf_1[k]
                 + pa_x[k] * df_2[k]
                 + f_2 * ff_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_x, pf_s_3, pf_2, dd_4, df_0, df_3, \
                         ff_s_6, ff_s_7, ff_s_8, fd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_z[k] * df_0[k]
                 + f_2 * ff_s_6[k];

        t_7[k] = -f_4 * pf_s_3[k]
                 + f_5 * pf_2[k]
                 + pa_x[k] * df_3[k]
                 + f_2 * ff_s_7[k];

        t_8[k] = f_5 * dd_4[k]
                 + f_2 * ff_s_8[k]
                 + pb_x[k] * fd_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_x, dd_8, df_6, df_13, ff_s_9, ff_s_10, \
                         ff_s_11, fd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_x[k] * df_6[k]
                 + f_2 * ff_s_9[k];

        t_10[k] = f_5 * dd_8[k]
                  + f_2 * ff_s_10[k]
                  + pb_x[k] * fd_7[k];

        t_11[k] = pa_x[k] * df_13[k]
                  + f_2 * ff_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_x, pb_y, dd_4, fp_s_2, fp_s_3, ff_s_12, ff_s_13, \
                         ff_s_14, fp_2, fp_3, fd_8, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_1 * fp_s_2[k]
                  + f_2 * ff_s_12[k]
                  + f_3 * fp_2[k]
                  + pb_x[k] * fd_8[k];

        t_13[k] = f_2 * ff_s_13[k]
                  + pb_x[k] * fd_9[k];

        t_14[k] = f_0 * dd_4[k]
                  - f_1 * fp_s_3[k]
                  + f_2 * ff_s_14[k]
                  + f_3 * fp_3[k]
                  + pb_y[k] * fd_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_z, pb_z, df_6, fp_s_4, ff_s_15, ff_s_16, \
                         ff_s_17, fp_4, fd_9, fd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_2 * ff_s_15[k]
                  + pb_z[k] * fd_9[k];

        t_16[k] = -f_1 * fp_s_4[k]
                  + f_2 * ff_s_16[k]
                  + f_3 * fp_4[k]
                  + pb_z[k] * fd_10[k];

        t_17[k] = pa_z[k] * df_6[k]
                  + f_2 * ff_s_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pb_y, pf_s_3, pf_2, dd_8, df_8, df_13, \
                         ff_s_18, ff_s_19, ff_s_20, fd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -f_4 * pf_s_3[k]
                  + f_5 * pf_2[k]
                  + pa_y[k] * df_8[k]
                  + f_2 * ff_s_18[k];

        t_19[k] = f_5 * dd_8[k]
                  + f_2 * ff_s_19[k]
                  + pb_y[k] * fd_13[k];

        t_20[k] = pa_y[k] * df_13[k]
                  + f_2 * ff_s_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_x, pb_y, fp_s_5, ff_s_21, ff_s_22, \
                         ff_s_23, ff_s_24, fp_5, fd_14, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_1 * fp_s_5[k]
                  + f_2 * ff_s_21[k]
                  + f_3 * fp_5[k]
                  + pb_x[k] * fd_14[k];

        t_22[k] = f_2 * ff_s_22[k]
                  + pb_y[k] * fd_14[k];

        t_23[k] = f_2 * ff_s_23[k]
                  + pb_x[k] * fd_15[k];

        t_24[k] = f_2 * ff_s_24[k]
                  + pb_x[k] * fd_16[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pb_y, pb_z, dd_8, fp_s_6, fp_s_7, ff_s_25, ff_s_26, \
                         ff_s_27, fp_6, fp_7, fd_15, fd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -f_1 * fp_s_6[k]
                  + f_2 * ff_s_25[k]
                  + f_3 * fp_6[k]
                  + pb_y[k] * fd_15[k];

        t_26[k] = f_2 * ff_s_26[k]
                  + pb_y[k] * fd_16[k];

        t_27[k] = f_0 * dd_8[k]
                  - f_1 * fp_s_7[k]
                  + f_2 * ff_s_27[k]
                  + f_3 * fp_7[k]
                  + pb_z[k] * fd_16[k];
    }
}

}  // namespace simdkin
