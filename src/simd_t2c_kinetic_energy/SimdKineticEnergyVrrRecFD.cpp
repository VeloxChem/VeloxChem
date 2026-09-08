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


#include "SimdKineticEnergyVrrRecFD.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_fd_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pd_s, const size_t pd,
                                 const size_t dp, const size_t dd, const size_t fs_s,
                                 const size_t fd_s, const size_t fs, const size_t fp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_5 = buffer.data(dp + 5);
    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_8 = buffer.data(dp + 8);
    const auto *dp_9 = buffer.data(dp + 9);
    const auto *dp_10 = buffer.data(dp + 10);
    const auto *dp_11 = buffer.data(dp + 11);

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
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_15 = buffer.data(dd + 15);
    const auto *dd_16 = buffer.data(dd + 16);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_5 = buffer.data(fs_s + 5);
    const auto *fs_s_8 = buffer.data(fs_s + 8);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_19 = buffer.data(fd_s + 19);
    const auto *fd_s_20 = buffer.data(fd_s + 20);
    const auto *fd_s_21 = buffer.data(fd_s + 21);
    const auto *fd_s_22 = buffer.data(fd_s + 22);
    const auto *fd_s_23 = buffer.data(fd_s + 23);
    const auto *fd_s_24 = buffer.data(fd_s + 24);
    const auto *fd_s_25 = buffer.data(fd_s + 25);
    const auto *fd_s_26 = buffer.data(fd_s + 26);
    const auto *fd_s_27 = buffer.data(fd_s + 27);
    const auto *fd_s_28 = buffer.data(fd_s + 28);
    const auto *fd_s_29 = buffer.data(fd_s + 29);
    const auto *fd_s_30 = buffer.data(fd_s + 30);
    const auto *fd_s_31 = buffer.data(fd_s + 31);
    const auto *fd_s_32 = buffer.data(fd_s + 32);
    const auto *fd_s_33 = buffer.data(fd_s + 33);
    const auto *fd_s_34 = buffer.data(fd_s + 34);
    const auto *fd_s_35 = buffer.data(fd_s + 35);
    const auto *fd_s_36 = buffer.data(fd_s + 36);
    const auto *fd_s_37 = buffer.data(fd_s + 37);
    const auto *fd_s_38 = buffer.data(fd_s + 38);
    const auto *fd_s_39 = buffer.data(fd_s + 39);
    const auto *fd_s_40 = buffer.data(fd_s + 40);
    const auto *fd_s_41 = buffer.data(fd_s + 41);
    const auto *fd_s_42 = buffer.data(fd_s + 42);
    const auto *fd_s_43 = buffer.data(fd_s + 43);
    const auto *fd_s_44 = buffer.data(fd_s + 44);
    const auto *fd_s_45 = buffer.data(fd_s + 45);
    const auto *fd_s_46 = buffer.data(fd_s + 46);
    const auto *fd_s_47 = buffer.data(fd_s + 47);
    const auto *fd_s_48 = buffer.data(fd_s + 48);
    const auto *fd_s_49 = buffer.data(fd_s + 49);
    const auto *fd_s_50 = buffer.data(fd_s + 50);
    const auto *fd_s_51 = buffer.data(fd_s + 51);
    const auto *fd_s_52 = buffer.data(fd_s + 52);
    const auto *fd_s_53 = buffer.data(fd_s + 53);
    const auto *fd_s_54 = buffer.data(fd_s + 54);
    const auto *fd_s_55 = buffer.data(fd_s + 55);
    const auto *fd_s_56 = buffer.data(fd_s + 56);
    const auto *fd_s_57 = buffer.data(fd_s + 57);
    const auto *fd_s_58 = buffer.data(fd_s + 58);
    const auto *fd_s_59 = buffer.data(fd_s + 59);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_8 = buffer.data(fs + 8);

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
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);
    const auto *fp_15 = buffer.data(fp + 15);
    const auto *fp_16 = buffer.data(fp + 16);
    const auto *fp_17 = buffer.data(fp + 17);
    const auto *fp_18 = buffer.data(fp + 18);
    const auto *fp_19 = buffer.data(fp + 19);
    const auto *fp_20 = buffer.data(fp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, dp_0, fs_s_0, fd_s_0, fd_s_1, \
                         fd_s_2, fd_s_3, fs_0, fp_0, fp_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = f_2 * fd_s_1[k]
                 + pb_y[k] * fp_0[k];

        t_2[k] = f_2 * fd_s_2[k]
                 + pb_z[k] * fp_0[k];

        t_3[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_3[k]
                 + f_3 * fs_0[k]
                 + pb_y[k] * fp_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_y, pb_y, pb_z, dd_0, fs_s_0, fd_s_4, fd_s_5, \
                         fd_s_6, fs_0, fp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * fd_s_4[k]
                 + pb_y[k] * fp_2[k];

        t_5[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_5[k]
                 + f_3 * fs_0[k]
                 + pb_z[k] * fp_2[k];

        t_6[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, pb_z, pd_s_1, pd_1, dp_3, dd_4, fd_s_7, \
                         fd_s_8, fd_s_9, fp_3, fp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_4 * dp_3[k]
                 + f_2 * fd_s_7[k]
                 + pb_x[k] * fp_4[k];

        t_8[k] = f_2 * fd_s_8[k]
                 + pb_z[k] * fp_3[k];

        t_9[k] = -f_5 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_4[k]
                 + f_2 * fd_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_y, pb_z, dd_0, dd_2, fd_s_10, \
                         fd_s_11, fd_s_12, fd_s_13, fp_4, fp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * fd_s_10[k]
                  + pb_z[k] * fp_4[k];

        t_11[k] = pa_y[k] * dd_2[k]
                  + f_2 * fd_s_11[k];

        t_12[k] = pa_z[k] * dd_0[k]
                  + f_2 * fd_s_12[k];

        t_13[k] = f_2 * fd_s_13[k]
                  + pb_y[k] * fp_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_z, pb_x, pb_y, dp_4, dd_1, fd_s_14, fd_s_15, \
                         fd_s_16, fp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_4 * dp_4[k]
                  + f_2 * fd_s_14[k]
                  + pb_x[k] * fp_6[k];

        t_15[k] = pa_z[k] * dd_1[k]
                  + f_2 * fd_s_15[k];

        t_16[k] = f_2 * fd_s_16[k]
                  + pb_y[k] * fp_6[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pb_x, pd_s_2, pd_2, dp_5, dp_6, dd_7, dd_8, \
                         fd_s_17, fd_s_18, fd_s_19, fp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_5 * pd_s_2[k]
                  + f_3 * pd_2[k]
                  + pa_x[k] * dd_7[k]
                  + f_2 * fd_s_17[k];

        t_18[k] = f_4 * dp_5[k]
                  + pa_x[k] * dd_8[k]
                  + f_2 * fd_s_18[k];

        t_19[k] = f_3 * dp_6[k]
                  + f_2 * fd_s_19[k]
                  + pb_x[k] * fp_8[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_z, dd_9, dd_10, fd_s_20, fd_s_21, \
                         fd_s_22, fd_s_23, fp_7, fp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_2 * fd_s_20[k]
                  + pb_z[k] * fp_7[k];

        t_21[k] = pa_x[k] * dd_9[k]
                  + f_2 * fd_s_21[k];

        t_22[k] = f_2 * fd_s_22[k]
                  + pb_z[k] * fp_8[k];

        t_23[k] = pa_x[k] * dd_10[k]
                  + f_2 * fd_s_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_x, pa_y, pa_z, dd_3, dd_5, dd_6, dd_11, \
                         fd_s_24, fd_s_25, fd_s_26, fd_s_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * dd_5[k]
                  + f_2 * fd_s_24[k];

        t_25[k] = pa_z[k] * dd_3[k]
                  + f_2 * fd_s_25[k];

        t_26[k] = pa_y[k] * dd_6[k]
                  + f_2 * fd_s_26[k];

        t_27[k] = pa_x[k] * dd_11[k]
                  + f_2 * fd_s_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_y, dp_9, dd_12, dd_13, dd_14, \
                         fd_s_28, fd_s_29, fd_s_30, fd_s_31, fp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_x[k] * dd_12[k]
                  + f_2 * fd_s_28[k];

        t_29[k] = pa_x[k] * dd_13[k]
                  + f_2 * fd_s_29[k];

        t_30[k] = f_4 * dp_9[k]
                  + pa_x[k] * dd_14[k]
                  + f_2 * fd_s_30[k];

        t_31[k] = f_2 * fd_s_31[k]
                  + pb_y[k] * fp_9[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_y, dp_11, dd_15, dd_16, \
                         fd_s_32, fd_s_33, fd_s_34, fd_s_35, fp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * dp_11[k]
                  + f_2 * fd_s_32[k]
                  + pb_x[k] * fp_10[k];

        t_33[k] = pa_x[k] * dd_15[k]
                  + f_2 * fd_s_33[k];

        t_34[k] = f_2 * fd_s_34[k]
                  + pb_y[k] * fp_10[k];

        t_35[k] = pa_x[k] * dd_16[k]
                  + f_2 * fd_s_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_x, pb_y, dp_6, fs_s_5, fd_s_36, fd_s_37, \
                         fd_s_38, fd_s_39, fs_5, fp_11, fp_12, fp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_36[k]
                  + f_3 * fs_5[k]
                  + pb_x[k] * fp_11[k];

        t_37[k] = f_2 * fd_s_37[k]
                  + pb_x[k] * fp_12[k];

        t_38[k] = f_2 * fd_s_38[k]
                  + pb_x[k] * fp_13[k];

        t_39[k] = f_0 * dp_6[k]
                  - f_1 * fs_s_5[k]
                  + f_2 * fd_s_39[k]
                  + f_3 * fs_5[k]
                  + pb_y[k] * fp_12[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_z, pb_z, dd_8, fs_s_5, fd_s_40, fd_s_41, \
                         fd_s_42, fs_5, fp_12, fp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_2 * fd_s_40[k]
                  + pb_z[k] * fp_12[k];

        t_41[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_41[k]
                  + f_3 * fs_5[k]
                  + pb_z[k] * fp_13[k];

        t_42[k] = pa_z[k] * dd_8[k]
                  + f_2 * fd_s_42[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_z, pb_x, pb_y, dp_8, dd_9, fd_s_43, \
                         fd_s_44, fd_s_45, fd_s_46, fp_14, fp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_2 * fd_s_43[k]
                  + pb_x[k] * fp_14[k];

        t_44[k] = f_2 * fd_s_44[k]
                  + pb_x[k] * fp_15[k];

        t_45[k] = pa_z[k] * dd_9[k]
                  + f_2 * fd_s_45[k];

        t_46[k] = f_4 * dp_8[k]
                  + f_2 * fd_s_46[k]
                  + pb_y[k] * fp_15[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pb_x, pd_s_2, pd_2, dd_13, dd_14, \
                         fd_s_47, fd_s_48, fd_s_49, fd_s_50, fp_16, \
                         fp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -f_5 * pd_s_2[k]
                  + f_3 * pd_2[k]
                  + pa_y[k] * dd_13[k]
                  + f_2 * fd_s_47[k];

        t_48[k] = pa_y[k] * dd_14[k]
                  + f_2 * fd_s_48[k];

        t_49[k] = f_2 * fd_s_49[k]
                  + pb_x[k] * fp_16[k];

        t_50[k] = f_2 * fd_s_50[k]
                  + pb_x[k] * fp_17[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pa_y, pb_y, dp_10, dp_11, dd_15, dd_16, fd_s_51, \
                         fd_s_52, fd_s_53, fp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_4 * dp_10[k]
                  + pa_y[k] * dd_15[k]
                  + f_2 * fd_s_51[k];

        t_52[k] = f_3 * dp_11[k]
                  + f_2 * fd_s_52[k]
                  + pb_y[k] * fp_17[k];

        t_53[k] = pa_y[k] * dd_16[k]
                  + f_2 * fd_s_53[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pb_x, pb_y, fs_s_8, fd_s_54, fd_s_55, \
                         fd_s_56, fd_s_57, fs_8, fp_18, fp_19, fp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = -f_1 * fs_s_8[k]
                  + f_2 * fd_s_54[k]
                  + f_3 * fs_8[k]
                  + pb_x[k] * fp_18[k];

        t_55[k] = f_2 * fd_s_55[k]
                  + pb_x[k] * fp_19[k];

        t_56[k] = f_2 * fd_s_56[k]
                  + pb_x[k] * fp_20[k];

        t_57[k] = -f_1 * fs_s_8[k]
                  + f_2 * fd_s_57[k]
                  + f_3 * fs_8[k]
                  + pb_y[k] * fp_19[k];
    }

#pragma omp simd aligned(t_58, t_59, pb_y, pb_z, dp_11, fs_s_8, fd_s_58, fd_s_59, fs_8, \
                         fp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_2 * fd_s_58[k]
                  + pb_y[k] * fp_20[k];

        t_59[k] = f_0 * dp_11[k]
                  - f_1 * fs_s_8[k]
                  + f_2 * fd_s_59[k]
                  + f_3 * fs_8[k]
                  + pb_z[k] * fp_20[k];
    }
}

auto
compute_prim_fd_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pd_s, const size_t pd,
                                 const size_t dp, const size_t dd, const size_t fs_s,
                                 const size_t fd_s, const size_t fs, const size_t fp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_s_3 = buffer.data(pd_s + 3);
    const auto *pd_s_5 = buffer.data(pd_s + 5);

    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_5 = buffer.data(pd + 5);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_5 = buffer.data(dp + 5);
    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_8 = buffer.data(dp + 8);
    const auto *dp_9 = buffer.data(dp + 9);
    const auto *dp_10 = buffer.data(dp + 10);
    const auto *dp_11 = buffer.data(dp + 11);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_16 = buffer.data(dd + 16);
    const auto *dd_18 = buffer.data(dd + 18);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_3 = buffer.data(fs_s + 3);
    const auto *fs_s_6 = buffer.data(fs_s + 6);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_19 = buffer.data(fd_s + 19);
    const auto *fd_s_20 = buffer.data(fd_s + 20);
    const auto *fd_s_21 = buffer.data(fd_s + 21);
    const auto *fd_s_22 = buffer.data(fd_s + 22);
    const auto *fd_s_23 = buffer.data(fd_s + 23);
    const auto *fd_s_24 = buffer.data(fd_s + 24);
    const auto *fd_s_25 = buffer.data(fd_s + 25);
    const auto *fd_s_26 = buffer.data(fd_s + 26);
    const auto *fd_s_27 = buffer.data(fd_s + 27);
    const auto *fd_s_28 = buffer.data(fd_s + 28);
    const auto *fd_s_29 = buffer.data(fd_s + 29);
    const auto *fd_s_31 = buffer.data(fd_s + 31);
    const auto *fd_s_32 = buffer.data(fd_s + 32);
    const auto *fd_s_33 = buffer.data(fd_s + 33);
    const auto *fd_s_34 = buffer.data(fd_s + 34);
    const auto *fd_s_35 = buffer.data(fd_s + 35);
    const auto *fd_s_36 = buffer.data(fd_s + 36);
    const auto *fd_s_37 = buffer.data(fd_s + 37);
    const auto *fd_s_38 = buffer.data(fd_s + 38);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_5 = buffer.data(fs + 5);

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
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dp_0, fs_s_0, fd_s_0, fd_s_1, \
                         fd_s_2, fs_0, fp_0, fp_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = f_2 * fd_s_1[k]
                 + pb_z[k] * fp_0[k];

        t_2[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_2[k]
                 + f_3 * fs_0[k]
                 + pb_y[k] * fp_1[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_y, pb_x, pb_z, dp_3, dd_0, fs_s_0, fd_s_3, fd_s_4, \
                         fd_s_5, fs_0, fp_2, fp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_3[k]
                 + f_3 * fs_0[k]
                 + pb_z[k] * fp_2[k];

        t_4[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_4[k];

        t_5[k] = f_4 * dp_3[k]
                 + f_2 * fd_s_5[k]
                 + pb_x[k] * fp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, pa_z, pd_s_3, pd_3, dd_0, dd_2, dd_4, \
                         fd_s_6, fd_s_7, fd_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_5 * pd_s_3[k]
                 + f_3 * pd_3[k]
                 + pa_x[k] * dd_4[k]
                 + f_2 * fd_s_6[k];

        t_7[k] = pa_y[k] * dd_2[k]
                 + f_2 * fd_s_7[k];

        t_8[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_x, pb_y, pd_s_5, pd_5, dp_4, dd_6, fd_s_9, \
                         fd_s_10, fd_s_11, fp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dp_4[k]
                 + f_2 * fd_s_9[k]
                 + pb_x[k] * fp_4[k];

        t_10[k] = f_2 * fd_s_10[k]
                  + pb_y[k] * fp_4[k];

        t_11[k] = -f_5 * pd_s_5[k]
                  + f_3 * pd_5[k]
                  + pa_x[k] * dd_6[k]
                  + f_2 * fd_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pb_x, dp_5, dp_6, dd_7, dd_9, dd_10, \
                         fd_s_12, fd_s_13, fd_s_14, fd_s_15, fp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_4 * dp_5[k]
                  + pa_x[k] * dd_7[k]
                  + f_2 * fd_s_12[k];

        t_13[k] = f_3 * dp_6[k]
                  + f_2 * fd_s_13[k]
                  + pb_x[k] * fp_5[k];

        t_14[k] = pa_x[k] * dd_9[k]
                  + f_2 * fd_s_14[k];

        t_15[k] = pa_x[k] * dd_10[k]
                  + f_2 * fd_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_x, dp_9, dp_11, dd_12, dd_14, dd_16, \
                         fd_s_16, fd_s_17, fd_s_18, fd_s_19, fp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_x[k] * dd_12[k]
                  + f_2 * fd_s_16[k];

        t_17[k] = f_4 * dp_9[k]
                  + pa_x[k] * dd_14[k]
                  + f_2 * fd_s_17[k];

        t_18[k] = f_3 * dp_11[k]
                  + f_2 * fd_s_18[k]
                  + pb_x[k] * fp_6[k];

        t_19[k] = pa_x[k] * dd_16[k]
                  + f_2 * fd_s_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_x, pb_x, dd_18, fs_s_3, fd_s_20, fd_s_21, \
                         fd_s_22, fs_3, fp_7, fp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_x[k] * dd_18[k]
                  + f_2 * fd_s_20[k];

        t_21[k] = -f_1 * fs_s_3[k]
                  + f_2 * fd_s_21[k]
                  + f_3 * fs_3[k]
                  + pb_x[k] * fp_7[k];

        t_22[k] = f_2 * fd_s_22[k]
                  + pb_x[k] * fp_8[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, pb_z, dp_6, fs_s_3, fd_s_23, fd_s_24, \
                         fd_s_25, fs_3, fp_8, fp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * dp_6[k]
                  - f_1 * fs_s_3[k]
                  + f_2 * fd_s_23[k]
                  + f_3 * fs_3[k]
                  + pb_y[k] * fp_8[k];

        t_24[k] = f_2 * fd_s_24[k]
                  + pb_z[k] * fp_8[k];

        t_25[k] = -f_1 * fs_s_3[k]
                  + f_2 * fd_s_25[k]
                  + f_3 * fs_3[k]
                  + pb_z[k] * fp_9[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_x, pb_y, dp_8, dd_9, fd_s_26, fd_s_27, \
                         fd_s_28, fp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * fd_s_26[k]
                  + pb_x[k] * fp_10[k];

        t_27[k] = pa_z[k] * dd_9[k]
                  + f_2 * fd_s_27[k];

        t_28[k] = f_4 * dp_8[k]
                  + f_2 * fd_s_28[k]
                  + pb_y[k] * fp_10[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pb_y, pd_s_5, pd_5, dp_10, dp_11, dd_13, \
                         dd_16, fd_s_29, fd_s_31, fd_s_32, fp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -f_5 * pd_s_5[k]
                  + f_3 * pd_5[k]
                  + pa_y[k] * dd_13[k]
                  + f_2 * fd_s_29[k];

        t_30[k] = f_4 * dp_10[k]
                  + pa_y[k] * dd_16[k]
                  + f_2 * fd_s_31[k];

        t_31[k] = f_3 * dp_11[k]
                  + f_2 * fd_s_32[k]
                  + pb_y[k] * fp_11[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pb_x, dd_18, fs_s_6, fd_s_33, fd_s_34, \
                         fd_s_35, fs_5, fp_12, fp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pa_y[k] * dd_18[k]
                  + f_2 * fd_s_33[k];

        t_33[k] = -f_1 * fs_s_6[k]
                  + f_2 * fd_s_34[k]
                  + f_3 * fs_5[k]
                  + pb_x[k] * fp_12[k];

        t_34[k] = f_2 * fd_s_35[k]
                  + pb_x[k] * fp_14[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_y, pb_z, dp_11, fs_s_6, fd_s_36, fd_s_37, \
                         fd_s_38, fs_5, fp_13, fp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -f_1 * fs_s_6[k]
                  + f_2 * fd_s_36[k]
                  + f_3 * fs_5[k]
                  + pb_y[k] * fp_13[k];

        t_36[k] = f_2 * fd_s_37[k]
                  + pb_y[k] * fp_14[k];

        t_37[k] = f_0 * dp_11[k]
                  - f_1 * fs_s_6[k]
                  + f_2 * fd_s_38[k]
                  + f_3 * fs_5[k]
                  + pb_z[k] * fp_14[k];
    }
}

auto
compute_prim_fd_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pd_s, const size_t pd,
                                 const size_t dp, const size_t dd, const size_t fs_s,
                                 const size_t fd_s, const size_t fs, const size_t fp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 1.0 / p;

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

    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_4 = buffer.data(pd_s + 4);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_3 = buffer.data(pd + 3);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_7 = buffer.data(dp + 7);
    const auto *dp_8 = buffer.data(dp + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_15 = buffer.data(dd + 15);
    const auto *dd_17 = buffer.data(dd + 17);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_2 = buffer.data(fs_s + 2);
    const auto *fs_s_5 = buffer.data(fs_s + 5);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_20 = buffer.data(fd_s + 20);
    const auto *fd_s_22 = buffer.data(fd_s + 22);
    const auto *fd_s_24 = buffer.data(fd_s + 24);
    const auto *fd_s_25 = buffer.data(fd_s + 25);
    const auto *fd_s_26 = buffer.data(fd_s + 26);
    const auto *fd_s_27 = buffer.data(fd_s + 27);
    const auto *fd_s_28 = buffer.data(fd_s + 28);
    const auto *fd_s_29 = buffer.data(fd_s + 29);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_4 = buffer.data(fs + 4);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dp_0, fs_s_0, fd_s_0, fd_s_1, \
                         fd_s_2, fs_0, fp_0, fp_1, fp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_1[k]
                 + f_3 * fs_0[k]
                 + pb_y[k] * fp_1[k];

        t_2[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_2[k]
                 + f_3 * fs_0[k]
                 + pb_z[k] * fp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pa_x, pa_y, pa_z, pd_s_1, pd_1, dd_0, dd_2, dd_4, \
                         fd_s_3, fd_s_4, fd_s_5, fd_s_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_4[k]
                 + f_2 * fd_s_4[k];

        t_5[k] = pa_y[k] * dd_2[k]
                 + f_2 * fd_s_5[k];

        t_6[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_y, pd_s_4, pd_3, dp_3, dd_6, dd_7, fd_s_7, \
                         fd_s_8, fd_s_9, fp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * fd_s_7[k]
                 + pb_y[k] * fp_3[k];

        t_8[k] = -f_4 * pd_s_4[k]
                 + f_3 * pd_3[k]
                 + pa_x[k] * dd_6[k]
                 + f_2 * fd_s_8[k];

        t_9[k] = f_5 * dp_3[k]
                 + pa_x[k] * dd_7[k]
                 + f_2 * fd_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, dp_6, dd_9, dd_13, dd_17, fd_s_10, fd_s_11, \
                         fd_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * dd_9[k]
                  + f_2 * fd_s_10[k];

        t_11[k] = f_5 * dp_6[k]
                  + pa_x[k] * dd_13[k]
                  + f_2 * fd_s_11[k];

        t_12[k] = pa_x[k] * dd_17[k]
                  + f_2 * fd_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, pb_y, dp_4, fs_s_2, fd_s_13, fd_s_14, \
                         fd_s_15, fs_2, fp_4, fp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_1 * fs_s_2[k]
                  + f_2 * fd_s_13[k]
                  + f_3 * fs_2[k]
                  + pb_x[k] * fp_4[k];

        t_14[k] = f_2 * fd_s_14[k]
                  + pb_x[k] * fp_5[k];

        t_15[k] = f_0 * dp_4[k]
                  - f_1 * fs_s_2[k]
                  + f_2 * fd_s_15[k]
                  + f_3 * fs_2[k]
                  + pb_y[k] * fp_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_z, pb_x, pb_z, dd_9, fs_s_2, fd_s_16, fd_s_17, \
                         fd_s_18, fs_2, fp_6, fp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_1 * fs_s_2[k]
                  + f_2 * fd_s_16[k]
                  + f_3 * fs_2[k]
                  + pb_z[k] * fp_6[k];

        t_17[k] = f_2 * fd_s_17[k]
                  + pb_x[k] * fp_7[k];

        t_18[k] = pa_z[k] * dd_9[k]
                  + f_2 * fd_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pd_s_4, pd_3, dp_7, dd_12, dd_15, dd_17, \
                         fd_s_20, fd_s_22, fd_s_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_4 * pd_s_4[k]
                  + f_3 * pd_3[k]
                  + pa_y[k] * dd_12[k]
                  + f_2 * fd_s_20[k];

        t_20[k] = f_5 * dp_7[k]
                  + pa_y[k] * dd_15[k]
                  + f_2 * fd_s_22[k];

        t_21[k] = pa_y[k] * dd_17[k]
                  + f_2 * fd_s_24[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pb_x, pb_y, fs_s_5, fd_s_25, fd_s_26, \
                         fd_s_27, fd_s_28, fs_4, fp_8, fp_9, fp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_25[k]
                  + f_3 * fs_4[k]
                  + pb_x[k] * fp_8[k];

        t_23[k] = f_2 * fd_s_26[k]
                  + pb_x[k] * fp_10[k];

        t_24[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_27[k]
                  + f_3 * fs_4[k]
                  + pb_y[k] * fp_9[k];

        t_25[k] = f_2 * fd_s_28[k]
                  + pb_y[k] * fp_10[k];
    }

#pragma omp simd aligned(t_26, pb_z, dp_8, fs_s_5, fd_s_29, fs_4, \
                         fp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * dp_8[k]
                  - f_1 * fs_s_5[k]
                  + f_2 * fd_s_29[k]
                  + f_3 * fs_4[k]
                  + pb_z[k] * fp_10[k];
    }
}

auto
compute_prim_fd_kinetic_energy_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pd_s, const size_t pd,
                                 const size_t dp, const size_t dd, const size_t fs_s,
                                 const size_t fd_s, const size_t fs, const size_t fp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = beta / p;

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

    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_4 = buffer.data(pd_s + 4);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_3 = buffer.data(pd + 3);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_8 = buffer.data(dp + 8);

    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_9 = buffer.data(dd + 9);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_2 = buffer.data(fs_s + 2);
    const auto *fs_s_5 = buffer.data(fs_s + 5);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_21 = buffer.data(fd_s + 21);
    const auto *fd_s_22 = buffer.data(fd_s + 22);
    const auto *fd_s_23 = buffer.data(fd_s + 23);
    const auto *fd_s_24 = buffer.data(fd_s + 24);
    const auto *fd_s_25 = buffer.data(fd_s + 25);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_4 = buffer.data(fs + 4);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dp_0, fs_s_0, fd_s_0, fd_s_1, \
                         fd_s_2, fs_0, fp_0, fp_1, fp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_1[k]
                 + f_3 * fs_0[k]
                 + pb_y[k] * fp_1[k];

        t_2[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_2[k]
                 + f_3 * fs_0[k]
                 + pb_z[k] * fp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_y, pd_s_1, pd_s_4, pd_1, pd_3, dd_3, dd_4, \
                         fd_s_3, fd_s_5, fd_s_6, fp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_3[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = f_2 * fd_s_5[k]
                 + pb_y[k] * fp_3[k];

        t_5[k] = -f_4 * pd_s_4[k]
                 + f_3 * pd_3[k]
                 + pa_x[k] * dd_4[k]
                 + f_2 * fd_s_6[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, dp_4, fs_s_2, fd_s_9, fd_s_10, fd_s_11, \
                         fs_2, fp_4, fp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_1 * fs_s_2[k]
                 + f_2 * fd_s_9[k]
                 + f_3 * fs_2[k]
                 + pb_x[k] * fp_4[k];

        t_7[k] = f_2 * fd_s_10[k]
                 + pb_x[k] * fp_5[k];

        t_8[k] = f_0 * dp_4[k]
                 - f_1 * fs_s_2[k]
                 + f_2 * fd_s_11[k]
                 + f_3 * fs_2[k]
                 + pb_y[k] * fp_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pb_z, pd_s_4, pd_3, dd_9, fs_s_2, \
                         fd_s_12, fd_s_13, fd_s_16, fs_2, fp_6, fp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_1 * fs_s_2[k]
                 + f_2 * fd_s_12[k]
                 + f_3 * fs_2[k]
                 + pb_z[k] * fp_6[k];

        t_10[k] = f_2 * fd_s_13[k]
                  + pb_x[k] * fp_7[k];

        t_11[k] = -f_4 * pd_s_4[k]
                  + f_3 * pd_3[k]
                  + pa_y[k] * dd_9[k]
                  + f_2 * fd_s_16[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_x, pb_y, fs_s_5, fd_s_21, fd_s_22, \
                         fd_s_23, fd_s_24, fs_4, fp_8, fp_9, fp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_21[k]
                  + f_3 * fs_4[k]
                  + pb_x[k] * fp_8[k];

        t_13[k] = f_2 * fd_s_22[k]
                  + pb_x[k] * fp_10[k];

        t_14[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_23[k]
                  + f_3 * fs_4[k]
                  + pb_y[k] * fp_9[k];

        t_15[k] = f_2 * fd_s_24[k]
                  + pb_y[k] * fp_10[k];
    }

#pragma omp simd aligned(t_16, pb_z, dp_8, fs_s_5, fd_s_25, fs_4, \
                         fp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * dp_8[k]
                  - f_1 * fs_s_5[k]
                  + f_2 * fd_s_25[k]
                  + f_3 * fs_4[k]
                  + pb_z[k] * fp_10[k];
    }
}

auto
compute_prim_fd_kinetic_energy_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pd_s, const size_t pd,
                                 const size_t dp, const size_t dd, const size_t fs_s,
                                 const size_t fd_s, const size_t fs, const size_t fp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = beta / p;

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

    const auto *pd_s_0 = buffer.data(pd_s + 0);
    const auto *pd_s_1 = buffer.data(pd_s + 1);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_1 = buffer.data(dp + 1);
    const auto *dp_2 = buffer.data(dp + 2);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_5 = buffer.data(dp + 5);
    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_7 = buffer.data(dp + 7);
    const auto *dp_8 = buffer.data(dp + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_5 = buffer.data(fs_s + 5);
    const auto *fs_s_8 = buffer.data(fs_s + 8);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_19 = buffer.data(fd_s + 19);
    const auto *fd_s_20 = buffer.data(fd_s + 20);
    const auto *fd_s_21 = buffer.data(fd_s + 21);
    const auto *fd_s_22 = buffer.data(fd_s + 22);
    const auto *fd_s_23 = buffer.data(fd_s + 23);
    const auto *fd_s_24 = buffer.data(fd_s + 24);
    const auto *fd_s_25 = buffer.data(fd_s + 25);
    const auto *fd_s_26 = buffer.data(fd_s + 26);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_8 = buffer.data(fs + 8);

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
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);
    const auto *fp_15 = buffer.data(fp + 15);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dp_0, fs_s_0, fd_s_0, fd_s_1, \
                         fd_s_2, fs_0, fp_0, fp_1, fp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_1[k]
                 + f_3 * fs_0[k]
                 + pb_y[k] * fp_1[k];

        t_2[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_2[k]
                 + f_3 * fs_0[k]
                 + pb_z[k] * fp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_x, pd_s_0, pd_0, dp_1, dd_0, dd_1, \
                         fd_s_3, fd_s_4, fd_s_5, fp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = f_4 * dp_1[k]
                 + f_2 * fd_s_4[k]
                 + pb_x[k] * fp_3[k];

        t_5[k] = -f_5 * pd_s_0[k]
                 + f_3 * pd_0[k]
                 + pa_x[k] * dd_1[k]
                 + f_2 * fd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_x, pd_s_1, pd_1, dp_2, dd_0, dd_2, \
                         fd_s_6, fd_s_7, fd_s_8, fp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_6[k];

        t_7[k] = f_4 * dp_2[k]
                 + f_2 * fd_s_7[k]
                 + pb_x[k] * fp_4[k];

        t_8[k] = -f_5 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_2[k]
                 + f_2 * fd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_x, dp_3, dp_4, dd_3, dd_4, fd_s_9, fd_s_10, \
                         fd_s_11, fp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * dp_3[k]
                 + pa_x[k] * dd_3[k]
                 + f_2 * fd_s_9[k];

        t_10[k] = f_3 * dp_4[k]
                  + f_2 * fd_s_10[k]
                  + pb_x[k] * fp_5[k];

        t_11[k] = pa_x[k] * dd_4[k]
                  + f_2 * fd_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_x, dp_6, dp_8, dd_6, dd_8, fd_s_12, \
                         fd_s_13, fd_s_14, fp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_4 * dp_6[k]
                  + pa_x[k] * dd_6[k]
                  + f_2 * fd_s_12[k];

        t_13[k] = f_3 * dp_8[k]
                  + f_2 * fd_s_13[k]
                  + pb_x[k] * fp_6[k];

        t_14[k] = pa_x[k] * dd_8[k]
                  + f_2 * fd_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_x, pb_y, pb_z, dp_4, fs_s_5, fd_s_15, fd_s_16, \
                         fd_s_17, fs_5, fp_7, fp_8, fp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_15[k]
                  + f_3 * fs_5[k]
                  + pb_x[k] * fp_7[k];

        t_16[k] = f_0 * dp_4[k]
                  - f_1 * fs_s_5[k]
                  + f_2 * fd_s_16[k]
                  + f_3 * fs_5[k]
                  + pb_y[k] * fp_8[k];

        t_17[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_17[k]
                  + f_3 * fs_5[k]
                  + pb_z[k] * fp_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pa_z, pb_y, pd_s_1, pd_1, dp_5, dd_4, dd_5, \
                         fd_s_18, fd_s_19, fd_s_20, fp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pa_z[k] * dd_4[k]
                  + f_2 * fd_s_18[k];

        t_19[k] = f_4 * dp_5[k]
                  + f_2 * fd_s_19[k]
                  + pb_y[k] * fp_10[k];

        t_20[k] = -f_5 * pd_s_1[k]
                  + f_3 * pd_1[k]
                  + pa_y[k] * dd_5[k]
                  + f_2 * fd_s_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_y, pb_y, dp_7, dp_8, dd_7, dd_8, fd_s_21, \
                         fd_s_22, fd_s_23, fp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_4 * dp_7[k]
                  + pa_y[k] * dd_7[k]
                  + f_2 * fd_s_21[k];

        t_22[k] = f_3 * dp_8[k]
                  + f_2 * fd_s_22[k]
                  + pb_y[k] * fp_12[k];

        t_23[k] = pa_y[k] * dd_8[k]
                  + f_2 * fd_s_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_x, pb_y, pb_z, dp_8, fs_s_8, fd_s_24, fd_s_25, \
                         fd_s_26, fs_8, fp_13, fp_14, fp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_1 * fs_s_8[k]
                  + f_2 * fd_s_24[k]
                  + f_3 * fs_8[k]
                  + pb_x[k] * fp_13[k];

        t_25[k] = -f_1 * fs_s_8[k]
                  + f_2 * fd_s_25[k]
                  + f_3 * fs_8[k]
                  + pb_y[k] * fp_14[k];

        t_26[k] = f_0 * dp_8[k]
                  - f_1 * fs_s_8[k]
                  + f_2 * fd_s_26[k]
                  + f_3 * fs_8[k]
                  + pb_z[k] * fp_15[k];
    }
}

auto
compute_prim_fd_kinetic_energy_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pd_s, const size_t pd,
                                 const size_t dp, const size_t dd, const size_t fs_s,
                                 const size_t fd_s, const size_t fs, const size_t fp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 1.0 / p;

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

    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_7 = buffer.data(dp + 7);
    const auto *dp_8 = buffer.data(dp + 8);
    const auto *dp_9 = buffer.data(dp + 9);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_15 = buffer.data(dd + 15);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_5 = buffer.data(fs_s + 5);
    const auto *fs_s_8 = buffer.data(fs_s + 8);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_19 = buffer.data(fd_s + 19);
    const auto *fd_s_20 = buffer.data(fd_s + 20);
    const auto *fd_s_21 = buffer.data(fd_s + 21);
    const auto *fd_s_22 = buffer.data(fd_s + 22);
    const auto *fd_s_23 = buffer.data(fd_s + 23);
    const auto *fd_s_24 = buffer.data(fd_s + 24);
    const auto *fd_s_25 = buffer.data(fd_s + 25);
    const auto *fd_s_26 = buffer.data(fd_s + 26);
    const auto *fd_s_27 = buffer.data(fd_s + 27);
    const auto *fd_s_28 = buffer.data(fd_s + 28);
    const auto *fd_s_29 = buffer.data(fd_s + 29);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_8 = buffer.data(fs + 8);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dp_0, fs_s_0, fd_s_0, fd_s_1, \
                         fd_s_2, fs_0, fp_0, fp_1, fp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_1[k]
                 + f_3 * fs_0[k]
                 + pb_y[k] * fp_1[k];

        t_2[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_2[k]
                 + f_3 * fs_0[k]
                 + pb_z[k] * fp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pa_x, pa_y, pa_z, pd_s_1, pd_1, dd_0, dd_2, dd_4, \
                         fd_s_3, fd_s_4, fd_s_5, fd_s_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_4[k]
                 + f_2 * fd_s_4[k];

        t_5[k] = pa_y[k] * dd_2[k]
                 + f_2 * fd_s_5[k];

        t_6[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pa_x, pd_s_2, pd_2, dp_3, dd_6, dd_7, dd_8, \
                         dd_9, fd_s_7, fd_s_8, fd_s_9, fd_s_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -f_4 * pd_s_2[k]
                 + f_3 * pd_2[k]
                 + pa_x[k] * dd_6[k]
                 + f_2 * fd_s_7[k];

        t_8[k] = f_5 * dp_3[k]
                 + pa_x[k] * dd_7[k]
                 + f_2 * fd_s_8[k];

        t_9[k] = pa_x[k] * dd_8[k]
                 + f_2 * fd_s_9[k];

        t_10[k] = pa_x[k] * dd_9[k]
                  + f_2 * fd_s_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pa_x, dp_7, dd_11, dd_13, dd_14, dd_15, \
                         fd_s_11, fd_s_12, fd_s_13, fd_s_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_x[k] * dd_11[k]
                  + f_2 * fd_s_11[k];

        t_12[k] = f_5 * dp_7[k]
                  + pa_x[k] * dd_13[k]
                  + f_2 * fd_s_12[k];

        t_13[k] = pa_x[k] * dd_14[k]
                  + f_2 * fd_s_13[k];

        t_14[k] = pa_x[k] * dd_15[k]
                  + f_2 * fd_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_x, pb_y, dp_4, fs_s_5, fd_s_15, fd_s_16, \
                         fd_s_17, fs_5, fp_7, fp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_15[k]
                  + f_3 * fs_5[k]
                  + pb_x[k] * fp_7[k];

        t_16[k] = f_2 * fd_s_16[k]
                  + pb_x[k] * fp_8[k];

        t_17[k] = f_0 * dp_4[k]
                  - f_1 * fs_s_5[k]
                  + f_2 * fd_s_17[k]
                  + f_3 * fs_5[k]
                  + pb_y[k] * fp_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_z, pb_y, pb_z, dp_6, dd_8, fs_s_5, fd_s_18, \
                         fd_s_19, fd_s_20, fs_5, fp_9, fp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_18[k]
                  + f_3 * fs_5[k]
                  + pb_z[k] * fp_9[k];

        t_19[k] = pa_z[k] * dd_8[k]
                  + f_2 * fd_s_19[k];

        t_20[k] = f_5 * dp_6[k]
                  + f_2 * fd_s_20[k]
                  + pb_y[k] * fp_10[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_y, pb_y, pd_s_2, pd_2, dp_8, dp_9, dd_12, dd_14, \
                         fd_s_21, fd_s_22, fd_s_23, fp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_4 * pd_s_2[k]
                  + f_3 * pd_2[k]
                  + pa_y[k] * dd_12[k]
                  + f_2 * fd_s_21[k];

        t_22[k] = f_5 * dp_8[k]
                  + pa_y[k] * dd_14[k]
                  + f_2 * fd_s_22[k];

        t_23[k] = f_3 * dp_9[k]
                  + f_2 * fd_s_23[k]
                  + pb_y[k] * fp_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_y, pb_x, dd_15, fs_s_8, fd_s_24, fd_s_25, \
                         fd_s_26, fs_8, fp_12, fp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * dd_15[k]
                  + f_2 * fd_s_24[k];

        t_25[k] = -f_1 * fs_s_8[k]
                  + f_2 * fd_s_25[k]
                  + f_3 * fs_8[k]
                  + pb_x[k] * fp_12[k];

        t_26[k] = f_2 * fd_s_26[k]
                  + pb_x[k] * fp_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pb_z, dp_9, fs_s_8, fd_s_27, fd_s_28, \
                         fd_s_29, fs_8, fp_13, fp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_1 * fs_s_8[k]
                  + f_2 * fd_s_27[k]
                  + f_3 * fs_8[k]
                  + pb_y[k] * fp_13[k];

        t_28[k] = f_2 * fd_s_28[k]
                  + pb_y[k] * fp_14[k];

        t_29[k] = f_0 * dp_9[k]
                  - f_1 * fs_s_8[k]
                  + f_2 * fd_s_29[k]
                  + f_3 * fs_8[k]
                  + pb_z[k] * fp_14[k];
    }
}

auto
compute_prim_fd_kinetic_energy_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pd_s, const size_t pd,
                                 const size_t dp, const size_t dd, const size_t fs_s,
                                 const size_t fd_s, const size_t fs, const size_t fp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_3 = buffer.data(pd_s + 3);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_3 = buffer.data(pd + 3);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_7 = buffer.data(dp + 7);
    const auto *dp_8 = buffer.data(dp + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_14 = buffer.data(dd + 14);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_4 = buffer.data(fs_s + 4);
    const auto *fs_s_7 = buffer.data(fs_s + 7);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_19 = buffer.data(fd_s + 19);
    const auto *fd_s_20 = buffer.data(fd_s + 20);
    const auto *fd_s_21 = buffer.data(fd_s + 21);
    const auto *fd_s_22 = buffer.data(fd_s + 22);
    const auto *fd_s_23 = buffer.data(fd_s + 23);
    const auto *fd_s_24 = buffer.data(fd_s + 24);
    const auto *fd_s_25 = buffer.data(fd_s + 25);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_7 = buffer.data(fs + 7);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dp_0, fs_s_0, fd_s_0, fd_s_1, \
                         fd_s_2, fs_0, fp_0, fp_1, fp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_1[k]
                 + f_3 * fs_0[k]
                 + pb_y[k] * fp_1[k];

        t_2[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_2[k]
                 + f_3 * fs_0[k]
                 + pb_z[k] * fp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pa_x, pa_y, pa_z, pd_s_1, pd_1, dd_0, dd_2, dd_3, \
                         fd_s_3, fd_s_4, fd_s_5, fd_s_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_3[k]
                 + f_2 * fd_s_4[k];

        t_5[k] = pa_y[k] * dd_2[k]
                 + f_2 * fd_s_5[k];

        t_6[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_y, pd_s_3, pd_3, dp_3, dd_5, dd_6, fd_s_7, \
                         fd_s_8, fd_s_9, fp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * fd_s_7[k]
                 + pb_y[k] * fp_3[k];

        t_8[k] = -f_4 * pd_s_3[k]
                 + f_3 * pd_3[k]
                 + pa_x[k] * dd_5[k]
                 + f_2 * fd_s_8[k];

        t_9[k] = f_5 * dp_3[k]
                 + pa_x[k] * dd_6[k]
                 + f_2 * fd_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, dp_6, dd_7, dd_11, dd_14, fd_s_10, fd_s_11, \
                         fd_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * dd_7[k]
                  + f_2 * fd_s_10[k];

        t_11[k] = f_5 * dp_6[k]
                  + pa_x[k] * dd_11[k]
                  + f_2 * fd_s_11[k];

        t_12[k] = pa_x[k] * dd_14[k]
                  + f_2 * fd_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, pb_y, dp_4, fs_s_4, fd_s_13, fd_s_14, \
                         fd_s_15, fs_4, fp_4, fp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_1 * fs_s_4[k]
                  + f_2 * fd_s_13[k]
                  + f_3 * fs_4[k]
                  + pb_x[k] * fp_4[k];

        t_14[k] = f_2 * fd_s_14[k]
                  + pb_x[k] * fp_5[k];

        t_15[k] = f_0 * dp_4[k]
                  - f_1 * fs_s_4[k]
                  + f_2 * fd_s_15[k]
                  + f_3 * fs_4[k]
                  + pb_y[k] * fp_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_y, pa_z, pb_z, pd_s_3, pd_3, dd_7, dd_10, \
                         fs_s_4, fd_s_16, fd_s_17, fd_s_18, fs_4, \
                         fp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_1 * fs_s_4[k]
                  + f_2 * fd_s_16[k]
                  + f_3 * fs_4[k]
                  + pb_z[k] * fp_6[k];

        t_17[k] = pa_z[k] * dd_7[k]
                  + f_2 * fd_s_17[k];

        t_18[k] = -f_4 * pd_s_3[k]
                  + f_3 * pd_3[k]
                  + pa_y[k] * dd_10[k]
                  + f_2 * fd_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pb_x, dp_7, dd_12, dd_14, fs_s_7, fd_s_19, \
                         fd_s_20, fd_s_21, fs_7, fp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_5 * dp_7[k]
                  + pa_y[k] * dd_12[k]
                  + f_2 * fd_s_19[k];

        t_20[k] = pa_y[k] * dd_14[k]
                  + f_2 * fd_s_20[k];

        t_21[k] = -f_1 * fs_s_7[k]
                  + f_2 * fd_s_21[k]
                  + f_3 * fs_7[k]
                  + pb_x[k] * fp_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pb_x, pb_y, pb_z, dp_8, fs_s_7, fd_s_22, \
                         fd_s_23, fd_s_24, fd_s_25, fs_7, fp_9, fp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_2 * fd_s_22[k]
                  + pb_x[k] * fp_10[k];

        t_23[k] = -f_1 * fs_s_7[k]
                  + f_2 * fd_s_23[k]
                  + f_3 * fs_7[k]
                  + pb_y[k] * fp_9[k];

        t_24[k] = f_2 * fd_s_24[k]
                  + pb_y[k] * fp_10[k];

        t_25[k] = f_0 * dp_8[k]
                  - f_1 * fs_s_7[k]
                  + f_2 * fd_s_25[k]
                  + f_3 * fs_7[k]
                  + pb_z[k] * fp_10[k];
    }
}

auto
compute_prim_fd_kinetic_energy_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pd_s, const size_t pd,
                                 const size_t dp, const size_t dd, const size_t fs_s,
                                 const size_t fd_s, const size_t fs, const size_t fp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_7 = buffer.data(dp + 7);
    const auto *dp_8 = buffer.data(dp + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_12 = buffer.data(dd + 12);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_4 = buffer.data(fs_s + 4);
    const auto *fs_s_7 = buffer.data(fs_s + 7);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_19 = buffer.data(fd_s + 19);
    const auto *fd_s_20 = buffer.data(fd_s + 20);
    const auto *fd_s_21 = buffer.data(fd_s + 21);
    const auto *fd_s_22 = buffer.data(fd_s + 22);
    const auto *fd_s_23 = buffer.data(fd_s + 23);
    const auto *fd_s_24 = buffer.data(fd_s + 24);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_7 = buffer.data(fs + 7);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dp_0, fs_s_0, fd_s_0, fd_s_1, \
                         fd_s_2, fs_0, fp_0, fp_1, fp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_1[k]
                 + f_3 * fs_0[k]
                 + pb_y[k] * fp_1[k];

        t_2[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_2[k]
                 + f_3 * fs_0[k]
                 + pb_z[k] * fp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pa_z, pd_s_1, pd_1, dd_0, dd_3, fd_s_3, \
                         fd_s_4, fd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_3[k]
                 + f_2 * fd_s_4[k];

        t_5[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, pb_y, pd_s_2, pd_2, dd_4, dd_6, dd_12, \
                         fd_s_6, fd_s_7, fd_s_9, fd_s_11, fp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * fd_s_6[k]
                 + pb_y[k] * fp_3[k];

        t_7[k] = -f_4 * pd_s_2[k]
                 + f_3 * pd_2[k]
                 + pa_x[k] * dd_4[k]
                 + f_2 * fd_s_7[k];

        t_8[k] = pa_x[k] * dd_6[k]
                 + f_2 * fd_s_9[k];

        t_9[k] = pa_x[k] * dd_12[k]
                 + f_2 * fd_s_11[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, pb_y, dp_4, fs_s_4, fd_s_12, fd_s_13, \
                         fd_s_14, fs_4, fp_4, fp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * fs_s_4[k]
                  + f_2 * fd_s_12[k]
                  + f_3 * fs_4[k]
                  + pb_x[k] * fp_4[k];

        t_11[k] = f_2 * fd_s_13[k]
                  + pb_x[k] * fp_5[k];

        t_12[k] = f_0 * dp_4[k]
                  - f_1 * fs_s_4[k]
                  + f_2 * fd_s_14[k]
                  + f_3 * fs_4[k]
                  + pb_y[k] * fp_5[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pa_z, pb_z, pd_s_2, pd_2, dd_6, dd_8, fs_s_4, \
                         fd_s_15, fd_s_16, fd_s_17, fs_4, fp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_1 * fs_s_4[k]
                  + f_2 * fd_s_15[k]
                  + f_3 * fs_4[k]
                  + pb_z[k] * fp_6[k];

        t_14[k] = pa_z[k] * dd_6[k]
                  + f_2 * fd_s_16[k];

        t_15[k] = -f_4 * pd_s_2[k]
                  + f_3 * pd_2[k]
                  + pa_y[k] * dd_8[k]
                  + f_2 * fd_s_17[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_y, pb_x, dp_7, dd_10, dd_12, fs_s_7, fd_s_18, \
                         fd_s_19, fd_s_20, fs_7, fp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * dp_7[k]
                  + pa_y[k] * dd_10[k]
                  + f_2 * fd_s_18[k];

        t_17[k] = pa_y[k] * dd_12[k]
                  + f_2 * fd_s_19[k];

        t_18[k] = -f_1 * fs_s_7[k]
                  + f_2 * fd_s_20[k]
                  + f_3 * fs_7[k]
                  + pb_x[k] * fp_8[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_x, pb_y, pb_z, dp_8, fs_s_7, fd_s_21, \
                         fd_s_22, fd_s_23, fd_s_24, fs_7, fp_9, fp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * fd_s_21[k]
                  + pb_x[k] * fp_10[k];

        t_20[k] = -f_1 * fs_s_7[k]
                  + f_2 * fd_s_22[k]
                  + f_3 * fs_7[k]
                  + pb_y[k] * fp_9[k];

        t_21[k] = f_2 * fd_s_23[k]
                  + pb_y[k] * fp_10[k];

        t_22[k] = f_0 * dp_8[k]
                  - f_1 * fs_s_7[k]
                  + f_2 * fd_s_24[k]
                  + f_3 * fs_7[k]
                  + pb_z[k] * fp_10[k];
    }
}

auto
compute_prim_fd_kinetic_energy_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pd_s, const size_t pd,
                                 const size_t dp, const size_t dd, const size_t fs_s,
                                 const size_t fd_s, const size_t fs, const size_t fp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = beta / p;

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

    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_8 = buffer.data(dp + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_12 = buffer.data(dd + 12);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_4 = buffer.data(fs_s + 4);
    const auto *fs_s_7 = buffer.data(fs_s + 7);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_19 = buffer.data(fd_s + 19);
    const auto *fd_s_20 = buffer.data(fd_s + 20);
    const auto *fd_s_21 = buffer.data(fd_s + 21);
    const auto *fd_s_22 = buffer.data(fd_s + 22);
    const auto *fd_s_23 = buffer.data(fd_s + 23);
    const auto *fd_s_24 = buffer.data(fd_s + 24);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_7 = buffer.data(fs + 7);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dp_0, fs_s_0, fd_s_0, fd_s_1, \
                         fd_s_2, fs_0, fp_0, fp_1, fp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_1[k]
                 + f_3 * fs_0[k]
                 + pb_y[k] * fp_1[k];

        t_2[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_2[k]
                 + f_3 * fs_0[k]
                 + pb_z[k] * fp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pa_z, pd_s_1, pd_1, dd_0, dd_3, fd_s_3, \
                         fd_s_4, fd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_3[k]
                 + f_2 * fd_s_4[k];

        t_5[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, pb_y, pd_s_2, pd_2, dd_4, dd_6, dd_12, \
                         fd_s_6, fd_s_7, fd_s_9, fd_s_11, fp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * fd_s_6[k]
                 + pb_y[k] * fp_3[k];

        t_7[k] = -f_4 * pd_s_2[k]
                 + f_3 * pd_2[k]
                 + pa_x[k] * dd_4[k]
                 + f_2 * fd_s_7[k];

        t_8[k] = pa_x[k] * dd_6[k]
                 + f_2 * fd_s_9[k];

        t_9[k] = pa_x[k] * dd_12[k]
                 + f_2 * fd_s_11[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, pb_y, dp_4, fs_s_4, fd_s_12, fd_s_13, \
                         fd_s_14, fs_4, fp_4, fp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * fs_s_4[k]
                  + f_2 * fd_s_12[k]
                  + f_3 * fs_4[k]
                  + pb_x[k] * fp_4[k];

        t_11[k] = f_2 * fd_s_13[k]
                  + pb_x[k] * fp_5[k];

        t_12[k] = f_0 * dp_4[k]
                  - f_1 * fs_s_4[k]
                  + f_2 * fd_s_14[k]
                  + f_3 * fs_4[k]
                  + pb_y[k] * fp_5[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pa_z, pb_z, pd_s_2, pd_2, dd_6, dd_8, fs_s_4, \
                         fd_s_15, fd_s_16, fd_s_17, fs_4, fp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_1 * fs_s_4[k]
                  + f_2 * fd_s_15[k]
                  + f_3 * fs_4[k]
                  + pb_z[k] * fp_6[k];

        t_14[k] = pa_z[k] * dd_6[k]
                  + f_2 * fd_s_16[k];

        t_15[k] = -f_4 * pd_s_2[k]
                  + f_3 * pd_2[k]
                  + pa_y[k] * dd_8[k]
                  + f_2 * fd_s_17[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_y, pb_x, dd_12, fs_s_7, fd_s_19, fd_s_20, \
                         fd_s_21, fs_7, fp_8, fp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_y[k] * dd_12[k]
                  + f_2 * fd_s_19[k];

        t_17[k] = -f_1 * fs_s_7[k]
                  + f_2 * fd_s_20[k]
                  + f_3 * fs_7[k]
                  + pb_x[k] * fp_8[k];

        t_18[k] = f_2 * fd_s_21[k]
                  + pb_x[k] * fp_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_y, pb_z, dp_8, fs_s_7, fd_s_22, fd_s_23, \
                         fd_s_24, fs_7, fp_9, fp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_1 * fs_s_7[k]
                  + f_2 * fd_s_22[k]
                  + f_3 * fs_7[k]
                  + pb_y[k] * fp_9[k];

        t_20[k] = f_2 * fd_s_23[k]
                  + pb_y[k] * fp_10[k];

        t_21[k] = f_0 * dp_8[k]
                  - f_1 * fs_s_7[k]
                  + f_2 * fd_s_24[k]
                  + f_3 * fs_7[k]
                  + pb_z[k] * fp_10[k];
    }
}

auto
compute_prim_fd_kinetic_energy_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pd_s, const size_t pd,
                                 const size_t dp, const size_t dd, const size_t fs_s,
                                 const size_t fd_s, const size_t fs, const size_t fp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 1.0 / p;

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

    const auto *pd_s_0 = buffer.data(pd_s + 0);
    const auto *pd_s_1 = buffer.data(pd_s + 1);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_1 = buffer.data(dp + 1);
    const auto *dp_2 = buffer.data(dp + 2);
    const auto *dp_3 = buffer.data(dp + 3);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_5 = buffer.data(fs_s + 5);
    const auto *fs_s_8 = buffer.data(fs_s + 8);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_8 = buffer.data(fs + 8);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_12 = buffer.data(fp + 12);

#pragma omp simd aligned(t_0, t_1, pa_y, pb_x, dp_0, dd_0, fs_s_0, fd_s_0, fd_s_1, fs_0, \
                         fp_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_1[k];
    }

#pragma omp simd aligned(t_2, t_3, t_4, pa_x, pa_z, pd_s_0, pd_s_1, pd_0, pd_1, dd_0, dd_1, \
                         dd_2, fd_s_2, fd_s_3, fd_s_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = -f_4 * pd_s_0[k]
                 + f_3 * pd_0[k]
                 + pa_x[k] * dd_1[k]
                 + f_2 * fd_s_2[k];

        t_3[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_2[k]
                 + f_2 * fd_s_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_x, pb_x, dd_3, dd_6, fs_s_5, fd_s_5, fd_s_6, \
                         fd_s_7, fs_5, fp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_x[k] * dd_3[k]
                 + f_2 * fd_s_5[k];

        t_6[k] = pa_x[k] * dd_6[k]
                 + f_2 * fd_s_6[k];

        t_7[k] = -f_1 * fs_s_5[k]
                 + f_2 * fd_s_7[k]
                 + f_3 * fs_5[k]
                 + pb_x[k] * fp_5[k];
    }

#pragma omp simd aligned(t_8, t_9, pa_z, pb_y, dp_1, dd_3, fs_s_5, fd_s_8, fd_s_9, fs_5, \
                         fp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dp_1[k]
                 - f_1 * fs_s_5[k]
                 + f_2 * fd_s_8[k]
                 + f_3 * fs_5[k]
                 + pb_y[k] * fp_6[k];

        t_9[k] = pa_z[k] * dd_3[k]
                 + f_2 * fd_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pd_s_1, pd_1, dp_2, dd_4, dd_5, dd_6, \
                         fd_s_10, fd_s_11, fd_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_4 * pd_s_1[k]
                  + f_3 * pd_1[k]
                  + pa_y[k] * dd_4[k]
                  + f_2 * fd_s_10[k];

        t_11[k] = f_5 * dp_2[k]
                  + pa_y[k] * dd_5[k]
                  + f_2 * fd_s_11[k];

        t_12[k] = pa_y[k] * dd_6[k]
                  + f_2 * fd_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, pb_y, pb_z, dp_3, fs_s_8, fd_s_13, fd_s_14, \
                         fd_s_15, fs_8, fp_10, fp_11, fp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_1 * fs_s_8[k]
                  + f_2 * fd_s_13[k]
                  + f_3 * fs_8[k]
                  + pb_x[k] * fp_10[k];

        t_14[k] = -f_1 * fs_s_8[k]
                  + f_2 * fd_s_14[k]
                  + f_3 * fs_8[k]
                  + pb_y[k] * fp_11[k];

        t_15[k] = f_0 * dp_3[k]
                  - f_1 * fs_s_8[k]
                  + f_2 * fd_s_15[k]
                  + f_3 * fs_8[k]
                  + pb_z[k] * fp_12[k];
    }
}

auto
compute_prim_fd_kinetic_energy_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pd_s, const size_t pd,
                                  const size_t dp, const size_t dd, const size_t fs_s,
                                  const size_t fd_s, const size_t fs, const size_t fp,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_1 = buffer.data(dp + 1);
    const auto *dp_2 = buffer.data(dp + 2);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_5 = buffer.data(dp + 5);
    const auto *dp_6 = buffer.data(dp + 6);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_5 = buffer.data(fs_s + 5);
    const auto *fs_s_8 = buffer.data(fs_s + 8);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_19 = buffer.data(fd_s + 19);
    const auto *fd_s_20 = buffer.data(fd_s + 20);
    const auto *fd_s_21 = buffer.data(fd_s + 21);
    const auto *fd_s_22 = buffer.data(fd_s + 22);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_8 = buffer.data(fs + 8);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dp_0, fs_s_0, fd_s_0, fd_s_1, \
                         fd_s_2, fs_0, fp_0, fp_1, fp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_1[k]
                 + f_3 * fs_0[k]
                 + pb_y[k] * fp_1[k];

        t_2[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_2[k]
                 + f_3 * fs_0[k]
                 + pb_z[k] * fp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pa_z, pd_s_1, pd_1, dd_0, dd_1, fd_s_3, \
                         fd_s_4, fd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_1[k]
                 + f_2 * fd_s_4[k];

        t_5[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pd_s_2, pd_2, dp_1, dd_2, dd_3, dd_4, fd_s_6, \
                         fd_s_7, fd_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * pd_s_2[k]
                 + f_3 * pd_2[k]
                 + pa_x[k] * dd_2[k]
                 + f_2 * fd_s_6[k];

        t_7[k] = f_5 * dp_1[k]
                 + pa_x[k] * dd_3[k]
                 + f_2 * fd_s_7[k];

        t_8[k] = pa_x[k] * dd_4[k]
                 + f_2 * fd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_x, dp_4, dd_6, dd_8, fs_s_5, fd_s_9, \
                         fd_s_10, fd_s_11, fs_5, fp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * dp_4[k]
                 + pa_x[k] * dd_6[k]
                 + f_2 * fd_s_9[k];

        t_10[k] = pa_x[k] * dd_8[k]
                  + f_2 * fd_s_10[k];

        t_11[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_11[k]
                  + f_3 * fs_5[k]
                  + pb_x[k] * fp_4[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_z, pb_y, pb_z, dp_2, dd_4, fs_s_5, fd_s_12, \
                         fd_s_13, fd_s_14, fs_5, fp_5, fp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * dp_2[k]
                  - f_1 * fs_s_5[k]
                  + f_2 * fd_s_12[k]
                  + f_3 * fs_5[k]
                  + pb_y[k] * fp_5[k];

        t_13[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_13[k]
                  + f_3 * fs_5[k]
                  + pb_z[k] * fp_6[k];

        t_14[k] = pa_z[k] * dd_4[k]
                  + f_2 * fd_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_y, pb_y, pd_s_2, pd_2, dp_3, dp_5, dd_5, dd_7, \
                         fd_s_15, fd_s_16, fd_s_17, fp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * dp_3[k]
                  + f_2 * fd_s_15[k]
                  + pb_y[k] * fp_7[k];

        t_16[k] = -f_4 * pd_s_2[k]
                  + f_3 * pd_2[k]
                  + pa_y[k] * dd_5[k]
                  + f_2 * fd_s_16[k];

        t_17[k] = f_5 * dp_5[k]
                  + pa_y[k] * dd_7[k]
                  + f_2 * fd_s_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pb_x, pb_y, dp_6, dd_8, fs_s_8, fd_s_18, \
                         fd_s_19, fd_s_20, fs_8, fp_8, fp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * dp_6[k]
                  + f_2 * fd_s_18[k]
                  + pb_y[k] * fp_8[k];

        t_19[k] = pa_y[k] * dd_8[k]
                  + f_2 * fd_s_19[k];

        t_20[k] = -f_1 * fs_s_8[k]
                  + f_2 * fd_s_20[k]
                  + f_3 * fs_8[k]
                  + pb_x[k] * fp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, pb_y, pb_z, dp_6, fs_s_8, fd_s_21, fd_s_22, fs_8, fp_10, \
                         fp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_1 * fs_s_8[k]
                  + f_2 * fd_s_21[k]
                  + f_3 * fs_8[k]
                  + pb_y[k] * fp_10[k];

        t_22[k] = f_0 * dp_6[k]
                  - f_1 * fs_s_8[k]
                  + f_2 * fd_s_22[k]
                  + f_3 * fs_8[k]
                  + pb_z[k] * fp_11[k];
    }
}

auto
compute_prim_fd_kinetic_energy_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pd_s, const size_t pd,
                                  const size_t dp, const size_t dd, const size_t fs_s,
                                  const size_t fd_s, const size_t fs, const size_t fp,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_2 = buffer.data(dp + 2);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_5 = buffer.data(dp + 5);
    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_7 = buffer.data(dp + 7);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_4 = buffer.data(fs_s + 4);
    const auto *fs_s_7 = buffer.data(fs_s + 7);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_19 = buffer.data(fd_s + 19);
    const auto *fd_s_20 = buffer.data(fd_s + 20);
    const auto *fd_s_21 = buffer.data(fd_s + 21);
    const auto *fd_s_22 = buffer.data(fd_s + 22);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_7 = buffer.data(fs + 7);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dp_0, fs_s_0, fd_s_0, fd_s_1, \
                         fd_s_2, fs_0, fp_0, fp_1, fp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_1[k]
                 + f_3 * fs_0[k]
                 + pb_y[k] * fp_1[k];

        t_2[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_2[k]
                 + f_3 * fs_0[k]
                 + pb_z[k] * fp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pa_x, pa_y, pa_z, pd_s_1, pd_1, dd_0, dd_1, dd_2, \
                         fd_s_3, fd_s_4, fd_s_5, fd_s_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_2[k]
                 + f_2 * fd_s_4[k];

        t_5[k] = pa_y[k] * dd_1[k]
                 + f_2 * fd_s_5[k];

        t_6[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pd_s_2, pd_2, dp_2, dd_4, dd_5, dd_6, fd_s_7, \
                         fd_s_8, fd_s_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -f_4 * pd_s_2[k]
                 + f_3 * pd_2[k]
                 + pa_x[k] * dd_4[k]
                 + f_2 * fd_s_7[k];

        t_8[k] = f_5 * dp_2[k]
                 + pa_x[k] * dd_5[k]
                 + f_2 * fd_s_8[k];

        t_9[k] = pa_x[k] * dd_6[k]
                 + f_2 * fd_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pb_x, dp_5, dd_10, dd_12, fs_s_4, fd_s_10, \
                         fd_s_11, fd_s_12, fs_4, fp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * dp_5[k]
                  + pa_x[k] * dd_10[k]
                  + f_2 * fd_s_10[k];

        t_11[k] = pa_x[k] * dd_12[k]
                  + f_2 * fd_s_11[k];

        t_12[k] = -f_1 * fs_s_4[k]
                  + f_2 * fd_s_12[k]
                  + f_3 * fs_4[k]
                  + pb_x[k] * fp_3[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, dp_3, dd_6, fs_s_4, fd_s_13, \
                         fd_s_14, fd_s_15, fs_4, fp_4, fp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * dp_3[k]
                  - f_1 * fs_s_4[k]
                  + f_2 * fd_s_13[k]
                  + f_3 * fs_4[k]
                  + pb_y[k] * fp_4[k];

        t_14[k] = -f_1 * fs_s_4[k]
                  + f_2 * fd_s_14[k]
                  + f_3 * fs_4[k]
                  + pb_z[k] * fp_5[k];

        t_15[k] = pa_z[k] * dd_6[k]
                  + f_2 * fd_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_y, pd_s_2, pd_2, dp_6, dd_9, dd_11, dd_12, \
                         fd_s_16, fd_s_17, fd_s_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_4 * pd_s_2[k]
                  + f_3 * pd_2[k]
                  + pa_y[k] * dd_9[k]
                  + f_2 * fd_s_16[k];

        t_17[k] = f_5 * dp_6[k]
                  + pa_y[k] * dd_11[k]
                  + f_2 * fd_s_17[k];

        t_18[k] = pa_y[k] * dd_12[k]
                  + f_2 * fd_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pb_x, pb_y, fs_s_7, fd_s_19, fd_s_20, fd_s_21, \
                         fs_7, fp_6, fp_7, fp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_1 * fs_s_7[k]
                  + f_2 * fd_s_19[k]
                  + f_3 * fs_7[k]
                  + pb_x[k] * fp_6[k];

        t_20[k] = -f_1 * fs_s_7[k]
                  + f_2 * fd_s_20[k]
                  + f_3 * fs_7[k]
                  + pb_y[k] * fp_7[k];

        t_21[k] = f_2 * fd_s_21[k]
                  + pb_y[k] * fp_8[k];
    }

#pragma omp simd aligned(t_22, pb_z, dp_7, fs_s_7, fd_s_22, fs_7, \
                         fp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_0 * dp_7[k]
                  - f_1 * fs_s_7[k]
                  + f_2 * fd_s_22[k]
                  + f_3 * fs_7[k]
                  + pb_z[k] * fp_8[k];
    }
}

auto
compute_prim_fd_kinetic_energy_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pd_s, const size_t pd,
                                  const size_t dp, const size_t dd, const size_t fs_s,
                                  const size_t fd_s, const size_t fs, const size_t fp,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 1.0 / p;

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
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_7 = buffer.data(dp + 7);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_4 = buffer.data(fs_s + 4);
    const auto *fs_s_7 = buffer.data(fs_s + 7);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_19 = buffer.data(fd_s + 19);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_7 = buffer.data(fs + 7);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dp_0, fs_s_0, fd_s_0, fd_s_1, \
                         fd_s_2, fs_0, fp_0, fp_1, fp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_1[k]
                 + f_3 * fs_0[k]
                 + pb_y[k] * fp_1[k];

        t_2[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_2[k]
                 + f_3 * fs_0[k]
                 + pb_z[k] * fp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pa_z, pd_s_1, pd_1, dd_0, dd_2, fd_s_3, \
                         fd_s_4, fd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_2[k]
                 + f_2 * fd_s_4[k];

        t_5[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pd_s_2, pd_2, dd_3, dd_5, dd_10, fd_s_6, fd_s_7, \
                         fd_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * pd_s_2[k]
                 + f_3 * pd_2[k]
                 + pa_x[k] * dd_3[k]
                 + f_2 * fd_s_6[k];

        t_7[k] = pa_x[k] * dd_5[k]
                 + f_2 * fd_s_7[k];

        t_8[k] = pa_x[k] * dd_10[k]
                 + f_2 * fd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, pb_y, pb_z, dp_3, fs_s_4, fd_s_9, fd_s_10, \
                         fd_s_11, fs_4, fp_3, fp_4, fp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_1 * fs_s_4[k]
                 + f_2 * fd_s_9[k]
                 + f_3 * fs_4[k]
                 + pb_x[k] * fp_3[k];

        t_10[k] = f_0 * dp_3[k]
                  - f_1 * fs_s_4[k]
                  + f_2 * fd_s_10[k]
                  + f_3 * fs_4[k]
                  + pb_y[k] * fp_4[k];

        t_11[k] = -f_1 * fs_s_4[k]
                  + f_2 * fd_s_11[k]
                  + f_3 * fs_4[k]
                  + pb_z[k] * fp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pa_z, pd_s_2, pd_2, dp_6, dd_5, dd_7, dd_9, \
                         fd_s_12, fd_s_13, fd_s_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_z[k] * dd_5[k]
                  + f_2 * fd_s_12[k];

        t_13[k] = -f_4 * pd_s_2[k]
                  + f_3 * pd_2[k]
                  + pa_y[k] * dd_7[k]
                  + f_2 * fd_s_13[k];

        t_14[k] = f_5 * dp_6[k]
                  + pa_y[k] * dd_9[k]
                  + f_2 * fd_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_y, pb_x, pb_y, dd_10, fs_s_7, fd_s_15, fd_s_16, \
                         fd_s_17, fs_7, fp_6, fp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * dd_10[k]
                  + f_2 * fd_s_15[k];

        t_16[k] = -f_1 * fs_s_7[k]
                  + f_2 * fd_s_16[k]
                  + f_3 * fs_7[k]
                  + pb_x[k] * fp_6[k];

        t_17[k] = -f_1 * fs_s_7[k]
                  + f_2 * fd_s_17[k]
                  + f_3 * fs_7[k]
                  + pb_y[k] * fp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, pb_y, pb_z, dp_7, fs_s_7, fd_s_18, fd_s_19, fs_7, \
                         fp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * fd_s_18[k]
                  + pb_y[k] * fp_8[k];

        t_19[k] = f_0 * dp_7[k]
                  - f_1 * fs_s_7[k]
                  + f_2 * fd_s_19[k]
                  + f_3 * fs_7[k]
                  + pb_z[k] * fp_8[k];
    }
}

auto
compute_prim_fd_kinetic_energy_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pd_s, const size_t pd,
                                  const size_t dp, const size_t dd, const size_t fs_s,
                                  const size_t fd_s, const size_t fs, const size_t fp,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = beta / p;

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
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_7 = buffer.data(dp + 7);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_10 = buffer.data(dd + 10);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_4 = buffer.data(fs_s + 4);
    const auto *fs_s_7 = buffer.data(fs_s + 7);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_19 = buffer.data(fd_s + 19);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_7 = buffer.data(fs + 7);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, dp_0, fs_s_0, fd_s_0, fd_s_1, \
                         fd_s_2, fs_0, fp_0, fp_1, fp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_1[k]
                 + f_3 * fs_0[k]
                 + pb_y[k] * fp_1[k];

        t_2[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_2[k]
                 + f_3 * fs_0[k]
                 + pb_z[k] * fp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pa_z, pd_s_1, pd_1, dd_0, dd_2, fd_s_3, \
                         fd_s_4, fd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_2[k]
                 + f_2 * fd_s_4[k];

        t_5[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pd_s_2, pd_2, dd_3, dd_5, dd_10, fd_s_6, fd_s_7, \
                         fd_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * pd_s_2[k]
                 + f_3 * pd_2[k]
                 + pa_x[k] * dd_3[k]
                 + f_2 * fd_s_6[k];

        t_7[k] = pa_x[k] * dd_5[k]
                 + f_2 * fd_s_7[k];

        t_8[k] = pa_x[k] * dd_10[k]
                 + f_2 * fd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, pb_y, pb_z, dp_3, fs_s_4, fd_s_9, fd_s_10, \
                         fd_s_11, fs_4, fp_3, fp_4, fp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_1 * fs_s_4[k]
                 + f_2 * fd_s_9[k]
                 + f_3 * fs_4[k]
                 + pb_x[k] * fp_3[k];

        t_10[k] = f_0 * dp_3[k]
                  - f_1 * fs_s_4[k]
                  + f_2 * fd_s_10[k]
                  + f_3 * fs_4[k]
                  + pb_y[k] * fp_4[k];

        t_11[k] = -f_1 * fs_s_4[k]
                  + f_2 * fd_s_11[k]
                  + f_3 * fs_4[k]
                  + pb_z[k] * fp_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pa_z, pd_s_2, pd_2, dd_5, dd_7, dd_10, \
                         fd_s_12, fd_s_13, fd_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_z[k] * dd_5[k]
                  + f_2 * fd_s_12[k];

        t_13[k] = -f_4 * pd_s_2[k]
                  + f_3 * pd_2[k]
                  + pa_y[k] * dd_7[k]
                  + f_2 * fd_s_13[k];

        t_14[k] = pa_y[k] * dd_10[k]
                  + f_2 * fd_s_15[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_x, pb_y, fs_s_7, fd_s_16, fd_s_17, fd_s_18, \
                         fs_7, fp_6, fp_7, fp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_1 * fs_s_7[k]
                  + f_2 * fd_s_16[k]
                  + f_3 * fs_7[k]
                  + pb_x[k] * fp_6[k];

        t_16[k] = -f_1 * fs_s_7[k]
                  + f_2 * fd_s_17[k]
                  + f_3 * fs_7[k]
                  + pb_y[k] * fp_7[k];

        t_17[k] = f_2 * fd_s_18[k]
                  + pb_y[k] * fp_8[k];
    }

#pragma omp simd aligned(t_18, pb_z, dp_7, fs_s_7, fd_s_19, fs_7, \
                         fp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * dp_7[k]
                  - f_1 * fs_s_7[k]
                  + f_2 * fd_s_19[k]
                  + f_3 * fs_7[k]
                  + pb_z[k] * fp_8[k];
    }
}

auto
compute_prim_fd_kinetic_energy_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pd_s, const size_t pd,
                                  const size_t dp, const size_t dd, const size_t fs_s,
                                  const size_t fd_s, const size_t fs, const size_t fp,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 1.0 / p;

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

    const auto *pd_s_0 = buffer.data(pd_s + 0);
    const auto *pd_s_1 = buffer.data(pd_s + 1);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_1 = buffer.data(dp + 1);
    const auto *dp_2 = buffer.data(dp + 2);
    const auto *dp_3 = buffer.data(dp + 3);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_5 = buffer.data(fs_s + 5);
    const auto *fs_s_8 = buffer.data(fs_s + 8);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_8 = buffer.data(fs + 8);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);

#pragma omp simd aligned(t_0, t_1, pa_y, pb_x, dp_0, dd_0, fs_s_0, fd_s_0, fd_s_1, fs_0, \
                         fp_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_1[k];
    }

#pragma omp simd aligned(t_2, t_3, t_4, pa_x, pa_z, pd_s_0, pd_s_1, pd_0, pd_1, dd_0, dd_1, \
                         dd_2, fd_s_2, fd_s_3, fd_s_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = -f_4 * pd_s_0[k]
                 + f_3 * pd_0[k]
                 + pa_x[k] * dd_1[k]
                 + f_2 * fd_s_2[k];

        t_3[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_2[k]
                 + f_2 * fd_s_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_x, pb_y, dp_1, dd_3, dd_6, fs_s_5, fd_s_5, fd_s_6, \
                         fd_s_7, fs_5, fp_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_x[k] * dd_3[k]
                 + f_2 * fd_s_5[k];

        t_6[k] = pa_x[k] * dd_6[k]
                 + f_2 * fd_s_6[k];

        t_7[k] = f_0 * dp_1[k]
                 - f_1 * fs_s_5[k]
                 + f_2 * fd_s_7[k]
                 + f_3 * fs_5[k]
                 + pb_y[k] * fp_1[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_y, pa_z, pd_s_1, pd_1, dp_2, dd_3, dd_4, dd_5, \
                         fd_s_8, fd_s_9, fd_s_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = pa_z[k] * dd_3[k]
                 + f_2 * fd_s_8[k];

        t_9[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_y[k] * dd_4[k]
                 + f_2 * fd_s_9[k];

        t_10[k] = f_5 * dp_2[k]
                  + pa_y[k] * dd_5[k]
                  + f_2 * fd_s_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_y, pb_y, pb_z, dp_3, dd_6, fs_s_8, fd_s_11, \
                         fd_s_12, fd_s_13, fs_8, fp_2, fp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_y[k] * dd_6[k]
                  + f_2 * fd_s_11[k];

        t_12[k] = -f_1 * fs_s_8[k]
                  + f_2 * fd_s_12[k]
                  + f_3 * fs_8[k]
                  + pb_y[k] * fp_2[k];

        t_13[k] = f_0 * dp_3[k]
                  - f_1 * fs_s_8[k]
                  + f_2 * fd_s_13[k]
                  + f_3 * fs_8[k]
                  + pb_z[k] * fp_3[k];
    }
}

auto
compute_prim_fd_kinetic_energy_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pd_s, const size_t pd,
                                  const size_t dp, const size_t dd, const size_t fs_s,
                                  const size_t fd_s, const size_t fs, const size_t fp,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 1.0 / p;

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

    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_1 = buffer.data(dp + 1);
    const auto *dp_2 = buffer.data(dp + 2);
    const auto *dp_3 = buffer.data(dp + 3);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_4 = buffer.data(fs_s + 4);
    const auto *fs_s_7 = buffer.data(fs_s + 7);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_7 = buffer.data(fs + 7);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);

#pragma omp simd aligned(t_0, t_1, pa_y, pb_x, dp_0, dd_0, fs_s_0, fd_s_0, fd_s_1, fs_0, \
                         fp_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_1[k];
    }

#pragma omp simd aligned(t_2, t_3, t_4, pa_x, pa_z, pd_s_1, pd_s_2, pd_1, pd_2, dd_0, dd_1, \
                         dd_2, fd_s_2, fd_s_3, fd_s_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_1[k]
                 + f_2 * fd_s_2[k];

        t_3[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = -f_4 * pd_s_2[k]
                 + f_3 * pd_2[k]
                 + pa_x[k] * dd_2[k]
                 + f_2 * fd_s_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_x, pb_x, dd_3, dd_6, fs_s_4, fd_s_5, fd_s_6, \
                         fd_s_7, fs_4, fp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_x[k] * dd_3[k]
                 + f_2 * fd_s_5[k];

        t_6[k] = pa_x[k] * dd_6[k]
                 + f_2 * fd_s_6[k];

        t_7[k] = -f_1 * fs_s_4[k]
                 + f_2 * fd_s_7[k]
                 + f_3 * fs_4[k]
                 + pb_x[k] * fp_2[k];
    }

#pragma omp simd aligned(t_8, t_9, pa_z, pb_y, dp_1, dd_3, fs_s_4, fd_s_8, fd_s_9, fs_4, \
                         fp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dp_1[k]
                 - f_1 * fs_s_4[k]
                 + f_2 * fd_s_8[k]
                 + f_3 * fs_4[k]
                 + pb_y[k] * fp_3[k];

        t_9[k] = pa_z[k] * dd_3[k]
                 + f_2 * fd_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pd_s_2, pd_2, dp_2, dd_4, dd_5, dd_6, \
                         fd_s_10, fd_s_11, fd_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_4 * pd_s_2[k]
                  + f_3 * pd_2[k]
                  + pa_y[k] * dd_4[k]
                  + f_2 * fd_s_10[k];

        t_11[k] = f_5 * dp_2[k]
                  + pa_y[k] * dd_5[k]
                  + f_2 * fd_s_11[k];

        t_12[k] = pa_y[k] * dd_6[k]
                  + f_2 * fd_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_x, pb_y, pb_z, dp_3, fs_s_7, fd_s_13, fd_s_14, \
                         fd_s_15, fs_7, fp_6, fp_7, fp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_1 * fs_s_7[k]
                  + f_2 * fd_s_13[k]
                  + f_3 * fs_7[k]
                  + pb_x[k] * fp_6[k];

        t_14[k] = -f_1 * fs_s_7[k]
                  + f_2 * fd_s_14[k]
                  + f_3 * fs_7[k]
                  + pb_y[k] * fp_7[k];

        t_15[k] = f_0 * dp_3[k]
                  - f_1 * fs_s_7[k]
                  + f_2 * fd_s_15[k]
                  + f_3 * fs_7[k]
                  + pb_z[k] * fp_8[k];
    }
}

auto
compute_prim_fd_kinetic_energy_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pd_s, const size_t pd,
                                  const size_t dp, const size_t dd, const size_t fs_s,
                                  const size_t fd_s, const size_t fs, const size_t fp,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = beta / p;
    const auto f_5 = 1.0 / p;

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
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_1 = buffer.data(dp + 1);
    const auto *dp_2 = buffer.data(dp + 2);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_5 = buffer.data(dp + 5);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_3 = buffer.data(fs_s + 3);
    const auto *fs_s_5 = buffer.data(fs_s + 5);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_19 = buffer.data(fd_s + 19);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_5 = buffer.data(fs + 5);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);

#pragma omp simd aligned(t_0, t_1, t_2, pa_y, pb_x, pb_z, dp_0, dd_0, fs_s_0, fd_s_0, fd_s_1, \
                         fd_s_2, fs_0, fp_0, fp_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_1[k]
                 + f_3 * fs_0[k]
                 + pb_z[k] * fp_1[k];

        t_2[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_z, pd_s_1, pd_s_2, pd_1, pd_2, dd_0, dd_1, \
                         dd_2, fd_s_3, fd_s_4, fd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_1[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_4[k];

        t_5[k] = -f_4 * pd_s_2[k]
                 + f_3 * pd_2[k]
                 + pa_x[k] * dd_2[k]
                 + f_2 * fd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, dp_1, dp_3, dd_3, dd_4, dd_6, dd_8, fd_s_6, \
                         fd_s_7, fd_s_8, fd_s_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * dp_1[k]
                 + pa_x[k] * dd_3[k]
                 + f_2 * fd_s_6[k];

        t_7[k] = pa_x[k] * dd_4[k]
                 + f_2 * fd_s_7[k];

        t_8[k] = f_5 * dp_3[k]
                 + pa_x[k] * dd_6[k]
                 + f_2 * fd_s_8[k];

        t_9[k] = pa_x[k] * dd_8[k]
                 + f_2 * fd_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, pb_y, pb_z, dp_2, fs_s_3, fd_s_10, fd_s_11, \
                         fd_s_12, fs_3, fp_2, fp_3, fp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * fs_s_3[k]
                  + f_2 * fd_s_10[k]
                  + f_3 * fs_3[k]
                  + pb_x[k] * fp_2[k];

        t_11[k] = f_0 * dp_2[k]
                  - f_1 * fs_s_3[k]
                  + f_2 * fd_s_11[k]
                  + f_3 * fs_3[k]
                  + pb_y[k] * fp_3[k];

        t_12[k] = -f_1 * fs_s_3[k]
                  + f_2 * fd_s_12[k]
                  + f_3 * fs_3[k]
                  + pb_z[k] * fp_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pa_z, pd_s_2, pd_2, dp_4, dd_4, dd_5, dd_7, \
                         fd_s_13, fd_s_14, fd_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_z[k] * dd_4[k]
                  + f_2 * fd_s_13[k];

        t_14[k] = -f_4 * pd_s_2[k]
                  + f_3 * pd_2[k]
                  + pa_y[k] * dd_5[k]
                  + f_2 * fd_s_14[k];

        t_15[k] = f_5 * dp_4[k]
                  + pa_y[k] * dd_7[k]
                  + f_2 * fd_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_y, pb_x, pb_y, dd_8, fs_s_5, fd_s_16, fd_s_17, \
                         fd_s_18, fs_5, fp_5, fp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_y[k] * dd_8[k]
                  + f_2 * fd_s_16[k];

        t_17[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_17[k]
                  + f_3 * fs_5[k]
                  + pb_x[k] * fp_5[k];

        t_18[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_18[k]
                  + f_3 * fs_5[k]
                  + pb_y[k] * fp_6[k];
    }

#pragma omp simd aligned(t_19, pb_z, dp_5, fs_s_5, fd_s_19, fs_5, \
                         fp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * dp_5[k]
                  - f_1 * fs_s_5[k]
                  + f_2 * fd_s_19[k]
                  + f_3 * fs_5[k]
                  + pb_z[k] * fp_7[k];
    }
}

auto
compute_prim_fd_kinetic_energy_17(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pd_s, const size_t pd,
                                  const size_t dp, const size_t dd, const size_t fs_s,
                                  const size_t fd_s, const size_t fs, const size_t fp,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = beta / p;

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
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_s_1 = buffer.data(pd_s + 1);
    const auto *pd_s_2 = buffer.data(pd_s + 2);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_2 = buffer.data(dp + 2);
    const auto *dp_5 = buffer.data(dp + 5);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_8 = buffer.data(dd + 8);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_3 = buffer.data(fs_s + 3);
    const auto *fs_s_5 = buffer.data(fs_s + 5);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_5 = buffer.data(fs + 5);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);

#pragma omp simd aligned(t_0, t_1, t_2, pa_y, pb_x, pb_z, dp_0, dd_0, fs_s_0, fd_s_0, fd_s_1, \
                         fd_s_2, fs_0, fp_0, fp_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 - f_1 * fs_s_0[k]
                 + f_2 * fd_s_0[k]
                 + f_3 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = -f_1 * fs_s_0[k]
                 + f_2 * fd_s_1[k]
                 + f_3 * fs_0[k]
                 + pb_z[k] * fp_1[k];

        t_2[k] = pa_y[k] * dd_0[k]
                 + f_2 * fd_s_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_z, pd_s_1, pd_s_2, pd_1, pd_2, dd_0, dd_1, \
                         dd_2, fd_s_3, fd_s_4, fd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_4 * pd_s_1[k]
                 + f_3 * pd_1[k]
                 + pa_x[k] * dd_1[k]
                 + f_2 * fd_s_3[k];

        t_4[k] = pa_z[k] * dd_0[k]
                 + f_2 * fd_s_4[k];

        t_5[k] = -f_4 * pd_s_2[k]
                 + f_3 * pd_2[k]
                 + pa_x[k] * dd_2[k]
                 + f_2 * fd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pb_x, dd_4, dd_8, fs_s_3, fd_s_6, fd_s_7, \
                         fd_s_8, fs_3, fp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_x[k] * dd_4[k]
                 + f_2 * fd_s_6[k];

        t_7[k] = pa_x[k] * dd_8[k]
                 + f_2 * fd_s_7[k];

        t_8[k] = -f_1 * fs_s_3[k]
                 + f_2 * fd_s_8[k]
                 + f_3 * fs_3[k]
                 + pb_x[k] * fp_2[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_z, pb_y, pb_z, dp_2, dd_4, fs_s_3, fd_s_9, \
                         fd_s_10, fd_s_11, fs_3, fp_3, fp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * dp_2[k]
                 - f_1 * fs_s_3[k]
                 + f_2 * fd_s_9[k]
                 + f_3 * fs_3[k]
                 + pb_y[k] * fp_3[k];

        t_10[k] = -f_1 * fs_s_3[k]
                  + f_2 * fd_s_10[k]
                  + f_3 * fs_3[k]
                  + pb_z[k] * fp_4[k];

        t_11[k] = pa_z[k] * dd_4[k]
                  + f_2 * fd_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pb_x, pd_s_2, pd_2, dd_5, dd_8, fs_s_5, \
                         fd_s_12, fd_s_13, fd_s_14, fs_5, fp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_4 * pd_s_2[k]
                  + f_3 * pd_2[k]
                  + pa_y[k] * dd_5[k]
                  + f_2 * fd_s_12[k];

        t_13[k] = pa_y[k] * dd_8[k]
                  + f_2 * fd_s_13[k];

        t_14[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_14[k]
                  + f_3 * fs_5[k]
                  + pb_x[k] * fp_5[k];
    }

#pragma omp simd aligned(t_15, t_16, pb_y, pb_z, dp_5, fs_s_5, fd_s_15, fd_s_16, fs_5, fp_6, \
                         fp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_1 * fs_s_5[k]
                  + f_2 * fd_s_15[k]
                  + f_3 * fs_5[k]
                  + pb_y[k] * fp_6[k];

        t_16[k] = f_0 * dp_5[k]
                  - f_1 * fs_s_5[k]
                  + f_2 * fd_s_16[k]
                  + f_3 * fs_5[k]
                  + pb_z[k] * fp_7[k];
    }
}

}  // namespace simdkin
