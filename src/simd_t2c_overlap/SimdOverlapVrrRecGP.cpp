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


#include "SimdOverlapVrrRecGP.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_gp_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fs, const size_t fp, const size_t gs,
                          const size_t ncols, const double p) -> void
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
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);
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

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_10 = buffer.data(gs + 10);
    const auto *gs_11 = buffer.data(gs + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pa_y, pa_z, pb_x, pb_y, pb_z, \
                         fs_0, fp_0, gs_0, gs_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pb_x[k] * gs_0[k];

        t_1[k] = pb_y[k] * gs_0[k];

        t_2[k] = pb_z[k] * gs_0[k];

        t_3[k] = pa_y[k] * fp_0[k];

        t_4[k] = f_1 * fs_0[k]
                 + pb_y[k] * gs_1[k];

        t_5[k] = pb_z[k] * gs_1[k];

        t_6[k] = pa_z[k] * fp_0[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pa_y, pb_x, pb_y, pb_z, fs_0, fs_1, \
                         fs_3, fp_2, gs_2, gs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_y[k] * gs_2[k];

        t_8[k] = f_1 * fs_0[k]
                 + pb_z[k] * gs_2[k];

        t_9[k] = f_2 * fs_3[k]
                 + pb_x[k] * gs_3[k];

        t_10[k] = f_2 * fs_1[k]
                  + pb_y[k] * gs_3[k];

        t_11[k] = pb_z[k] * gs_3[k];

        t_12[k] = pa_y[k] * fp_2[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_y, pa_z, pb_x, pb_y, pb_z, fs_2, \
                         fs_4, fp_1, fp_3, gs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_z[k] * fp_1[k];

        t_14[k] = pa_y[k] * fp_3[k];

        t_15[k] = f_2 * fs_4[k]
                  + pb_x[k] * gs_4[k];

        t_16[k] = pb_y[k] * gs_4[k];

        t_17[k] = f_2 * fs_2[k]
                  + pb_z[k] * gs_4[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, pa_x, pa_z, pb_x, pb_z, fs_5, \
                         fp_4, fp_6, fp_7, fp_8, gs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_1 * fs_5[k]
                  + pb_x[k] * gs_5[k];

        t_19[k] = pa_x[k] * fp_6[k];

        t_20[k] = pb_z[k] * gs_5[k];

        t_21[k] = pa_z[k] * fp_4[k];

        t_22[k] = pa_x[k] * fp_7[k];

        t_23[k] = pa_x[k] * fp_8[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, pa_x, pa_y, pb_x, pb_y, fs_8, \
                         fp_5, fp_9, fp_10, fp_11, gs_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * fp_5[k];

        t_25[k] = pa_x[k] * fp_9[k];

        t_26[k] = pa_x[k] * fp_10[k];

        t_27[k] = f_1 * fs_8[k]
                  + pb_x[k] * gs_6[k];

        t_28[k] = pb_y[k] * gs_6[k];

        t_29[k] = pa_x[k] * fp_11[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, t_36, pa_z, pb_x, pb_y, pb_z, \
                         fs_5, fp_6, gs_7, gs_8, gs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_x[k] * gs_7[k];

        t_31[k] = f_0 * fs_5[k]
                  + pb_y[k] * gs_7[k];

        t_32[k] = pb_z[k] * gs_7[k];

        t_33[k] = pb_x[k] * gs_8[k];

        t_34[k] = pa_z[k] * fp_6[k];

        t_35[k] = f_1 * fs_5[k]
                  + pb_z[k] * gs_8[k];

        t_36[k] = pb_x[k] * gs_9[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pa_y, pb_x, pb_y, pb_z, fs_6, fs_7, \
                         fs_8, fp_11, gs_9, gs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_2 * fs_7[k]
                  + pb_y[k] * gs_9[k];

        t_38[k] = f_2 * fs_6[k]
                  + pb_z[k] * gs_9[k];

        t_39[k] = pb_x[k] * gs_10[k];

        t_40[k] = f_1 * fs_8[k]
                  + pb_y[k] * gs_10[k];

        t_41[k] = pa_y[k] * fp_11[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, pb_y, pb_z, fs_8, \
                         gs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_x[k] * gs_11[k];

        t_43[k] = pb_y[k] * gs_11[k];

        t_44[k] = f_0 * fs_8[k]
                  + pb_z[k] * gs_11[k];
    }
}

auto
compute_prim_gp_overlap_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fs, const size_t fp, const size_t gs,
                          const size_t ncols, const double p) -> void
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);
    const auto *fs_8 = buffer.data(fs + 8);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_15 = buffer.data(fp + 15);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_10 = buffer.data(gs + 10);
    const auto *gs_11 = buffer.data(gs + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_y, pa_z, pb_x, pb_y, pb_z, fs_0, \
                         fp_0, gs_0, gs_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pb_x[k] * gs_0[k];

        t_1[k] = pb_y[k] * gs_0[k];

        t_2[k] = pb_z[k] * gs_0[k];

        t_3[k] = pa_y[k] * fp_0[k];

        t_4[k] = f_1 * fs_0[k]
                 + pb_y[k] * gs_1[k];

        t_5[k] = pa_z[k] * fp_0[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_y, pb_x, pb_y, pb_z, fs_0, fs_1, fs_3, \
                         fp_4, gs_2, gs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * fs_0[k]
                 + pb_z[k] * gs_2[k];

        t_7[k] = f_2 * fs_3[k]
                 + pb_x[k] * gs_3[k];

        t_8[k] = f_2 * fs_1[k]
                 + pb_y[k] * gs_3[k];

        t_9[k] = pb_z[k] * gs_3[k];

        t_10[k] = pa_y[k] * fp_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_x, pb_x, pb_y, pb_z, fs_2, fs_4, \
                         fs_5, fp_8, gs_4, gs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * fs_4[k]
                  + pb_x[k] * gs_4[k];

        t_12[k] = pb_y[k] * gs_4[k];

        t_13[k] = f_2 * fs_2[k]
                  + pb_z[k] * gs_4[k];

        t_14[k] = f_1 * fs_5[k]
                  + pb_x[k] * gs_5[k];

        t_15[k] = pa_x[k] * fp_8[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, pa_x, pb_x, pb_y, fs_5, fs_8, \
                         fp_10, fp_11, fp_15, gs_6, gs_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_x[k] * fp_10[k];

        t_17[k] = pa_x[k] * fp_11[k];

        t_18[k] = f_1 * fs_8[k]
                  + pb_x[k] * gs_6[k];

        t_19[k] = pa_x[k] * fp_15[k];

        t_20[k] = pb_x[k] * gs_7[k];

        t_21[k] = f_0 * fs_5[k]
                  + pb_y[k] * gs_7[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pa_z, pb_x, pb_y, pb_z, fs_5, fs_7, \
                         fp_8, gs_7, gs_8, gs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_z[k] * gs_7[k];

        t_23[k] = pa_z[k] * fp_8[k];

        t_24[k] = f_1 * fs_5[k]
                  + pb_z[k] * gs_8[k];

        t_25[k] = pb_x[k] * gs_9[k];

        t_26[k] = f_2 * fs_7[k]
                  + pb_y[k] * gs_9[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, fs_6, \
                         fs_8, fp_15, gs_9, gs_10, gs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_2 * fs_6[k]
                  + pb_z[k] * gs_9[k];

        t_28[k] = f_1 * fs_8[k]
                  + pb_y[k] * gs_10[k];

        t_29[k] = pa_y[k] * fp_15[k];

        t_30[k] = pb_x[k] * gs_11[k];

        t_31[k] = pb_y[k] * gs_11[k];

        t_32[k] = f_0 * fs_8[k]
                  + pb_z[k] * gs_11[k];
    }
}

auto
compute_prim_gp_overlap_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fs, const size_t fp, const size_t gs,
                          const size_t ncols, const double p) -> void
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);
    const auto *fs_8 = buffer.data(fs + 8);

    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_14 = buffer.data(fp + 14);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_9 = buffer.data(gs + 9);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pb_x, pb_y, pb_z, fs_0, fs_1, \
                         gs_0, gs_1, gs_2, gs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pb_x[k] * gs_0[k];

        t_1[k] = pb_y[k] * gs_0[k];

        t_2[k] = pb_z[k] * gs_0[k];

        t_3[k] = f_1 * fs_0[k]
                 + pb_y[k] * gs_1[k];

        t_4[k] = f_1 * fs_0[k]
                 + pb_z[k] * gs_2[k];

        t_5[k] = f_2 * fs_1[k]
                 + pb_y[k] * gs_3[k];

        t_6[k] = pb_z[k] * gs_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, pa_x, pa_y, pb_y, pb_z, fs_2, fp_4, fp_8, \
                         fp_14, gs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_y[k] * fp_4[k];

        t_8[k] = pb_y[k] * gs_4[k];

        t_9[k] = f_2 * fs_2[k]
                 + pb_z[k] * gs_4[k];

        t_10[k] = pa_x[k] * fp_8[k];

        t_11[k] = pa_x[k] * fp_14[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, fs_5, \
                         fs_6, fs_7, gs_5, gs_6, gs_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_x[k] * gs_5[k];

        t_13[k] = f_0 * fs_5[k]
                  + pb_y[k] * gs_5[k];

        t_14[k] = pb_z[k] * gs_5[k];

        t_15[k] = f_1 * fs_5[k]
                  + pb_z[k] * gs_6[k];

        t_16[k] = pb_x[k] * gs_7[k];

        t_17[k] = f_2 * fs_7[k]
                  + pb_y[k] * gs_7[k];

        t_18[k] = f_2 * fs_6[k]
                  + pb_z[k] * gs_7[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pb_x, pb_y, pb_z, fs_8, fp_14, \
                         gs_8, gs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * fs_8[k]
                  + pb_y[k] * gs_8[k];

        t_20[k] = pa_y[k] * fp_14[k];

        t_21[k] = pb_x[k] * gs_9[k];

        t_22[k] = pb_y[k] * gs_9[k];

        t_23[k] = f_0 * fs_8[k]
                  + pb_z[k] * gs_9[k];
    }
}

auto
compute_prim_gp_overlap_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fs, const size_t fp, const size_t gs,
                          const size_t ncols, const double p) -> void
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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);

    const auto *fp_11 = buffer.data(fp + 11);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pb_x, pb_y, pb_z, fs_0, fs_1, \
                         gs_0, gs_1, gs_2, gs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pb_x[k] * gs_0[k];

        t_1[k] = pb_y[k] * gs_0[k];

        t_2[k] = pb_z[k] * gs_0[k];

        t_3[k] = f_1 * fs_0[k]
                 + pb_z[k] * gs_1[k];

        t_4[k] = pb_z[k] * gs_2[k];

        t_5[k] = pb_y[k] * gs_3[k];

        t_6[k] = f_2 * fs_1[k]
                 + pb_z[k] * gs_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, t_13, pb_x, pb_y, pb_z, fs_4, fs_5, \
                         fs_6, gs_4, gs_5, gs_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_x[k] * gs_4[k];

        t_8[k] = f_0 * fs_4[k]
                 + pb_y[k] * gs_4[k];

        t_9[k] = pb_z[k] * gs_4[k];

        t_10[k] = f_1 * fs_4[k]
                  + pb_z[k] * gs_5[k];

        t_11[k] = pb_x[k] * gs_6[k];

        t_12[k] = f_2 * fs_6[k]
                  + pb_y[k] * gs_6[k];

        t_13[k] = f_2 * fs_5[k]
                  + pb_z[k] * gs_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_y, pb_x, pb_y, pb_z, fs_7, fp_11, \
                         gs_7, gs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * fs_7[k]
                  + pb_y[k] * gs_7[k];

        t_15[k] = pa_y[k] * fp_11[k];

        t_16[k] = pb_x[k] * gs_8[k];

        t_17[k] = pb_y[k] * gs_8[k];

        t_18[k] = f_0 * fs_7[k]
                  + pb_z[k] * gs_8[k];
    }
}

auto
compute_prim_gp_overlap_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fs, const size_t fp, const size_t gs,
                          const size_t ncols, const double p) -> void
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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);
    const auto *fs_8 = buffer.data(fs + 8);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_10 = buffer.data(gs + 10);
    const auto *gs_11 = buffer.data(gs + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_z, pb_x, pb_y, pb_z, fs_0, fs_3, fp_0, \
                         gs_0, gs_1, gs_2, gs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pb_x[k] * gs_0[k];

        t_1[k] = f_1 * fs_0[k]
                 + pb_y[k] * gs_1[k];

        t_2[k] = pa_z[k] * fp_0[k];

        t_3[k] = f_1 * fs_0[k]
                 + pb_z[k] * gs_2[k];

        t_4[k] = f_2 * fs_3[k]
                 + pb_x[k] * gs_3[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, fs_1, fs_2, fs_4, fs_5, gs_3, \
                         gs_4, gs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_2 * fs_1[k]
                 + pb_y[k] * gs_3[k];

        t_6[k] = f_2 * fs_4[k]
                 + pb_x[k] * gs_4[k];

        t_7[k] = f_2 * fs_2[k]
                 + pb_z[k] * gs_4[k];

        t_8[k] = f_1 * fs_5[k]
                 + pb_x[k] * gs_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_z, pb_x, pb_y, pb_z, fs_5, fs_8, fp_1, \
                         gs_6, gs_7, gs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_1 * fs_8[k]
                 + pb_x[k] * gs_6[k];

        t_10[k] = f_0 * fs_5[k]
                  + pb_y[k] * gs_7[k];

        t_11[k] = pa_z[k] * fp_1[k];

        t_12[k] = f_1 * fs_5[k]
                  + pb_z[k] * gs_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_y, pb_y, pb_z, fs_6, fs_7, fs_8, \
                         fp_2, gs_9, gs_10, gs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * fs_7[k]
                  + pb_y[k] * gs_9[k];

        t_14[k] = f_2 * fs_6[k]
                  + pb_z[k] * gs_9[k];

        t_15[k] = f_1 * fs_8[k]
                  + pb_y[k] * gs_10[k];

        t_16[k] = pa_y[k] * fp_2[k];

        t_17[k] = f_0 * fs_8[k]
                  + pb_z[k] * gs_11[k];
    }
}

auto
compute_prim_gp_overlap_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fs, const size_t fp, const size_t gs,
                          const size_t ncols, const double p) -> void
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);
    const auto *fs_8 = buffer.data(fs + 8);

    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_12 = buffer.data(fp + 12);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_10 = buffer.data(gs + 10);
    const auto *gs_11 = buffer.data(gs + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, fs_0, fs_1, gs_0, \
                         gs_1, gs_2, gs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pb_x[k] * gs_0[k];

        t_1[k] = pb_y[k] * gs_0[k];

        t_2[k] = pb_z[k] * gs_0[k];

        t_3[k] = f_1 * fs_0[k]
                 + pb_y[k] * gs_1[k];

        t_4[k] = f_1 * fs_0[k]
                 + pb_z[k] * gs_2[k];

        t_5[k] = f_2 * fs_1[k]
                 + pb_y[k] * gs_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pa_x, pa_y, pb_z, fs_2, fp_2, fp_6, \
                         fp_7, fp_8, fp_12, gs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * fp_2[k];

        t_7[k] = f_2 * fs_2[k]
                 + pb_z[k] * gs_4[k];

        t_8[k] = pa_x[k] * fp_6[k];

        t_9[k] = pa_x[k] * fp_7[k];

        t_10[k] = pa_x[k] * fp_8[k];

        t_11[k] = pa_x[k] * fp_12[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, fs_5, \
                         fs_6, fs_7, gs_7, gs_8, gs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_x[k] * gs_7[k];

        t_13[k] = f_0 * fs_5[k]
                  + pb_y[k] * gs_7[k];

        t_14[k] = pb_z[k] * gs_7[k];

        t_15[k] = f_1 * fs_5[k]
                  + pb_z[k] * gs_8[k];

        t_16[k] = pb_x[k] * gs_9[k];

        t_17[k] = f_2 * fs_7[k]
                  + pb_y[k] * gs_9[k];

        t_18[k] = f_2 * fs_6[k]
                  + pb_z[k] * gs_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pb_x, pb_y, pb_z, fs_8, fp_12, \
                         gs_10, gs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * fs_8[k]
                  + pb_y[k] * gs_10[k];

        t_20[k] = pa_y[k] * fp_12[k];

        t_21[k] = pb_x[k] * gs_11[k];

        t_22[k] = pb_y[k] * gs_11[k];

        t_23[k] = f_0 * fs_8[k]
                  + pb_z[k] * gs_11[k];
    }
}

auto
compute_prim_gp_overlap_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fs, const size_t fp, const size_t gs,
                          const size_t ncols, const double p) -> void
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_8 = buffer.data(fs + 8);

    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_11 = buffer.data(fp + 11);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_11 = buffer.data(gs + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pb_x, pb_y, pb_z, fs_0, fs_1, \
                         gs_0, gs_1, gs_2, gs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pb_x[k] * gs_0[k];

        t_1[k] = pb_y[k] * gs_0[k];

        t_2[k] = pb_z[k] * gs_0[k];

        t_3[k] = f_1 * fs_0[k]
                 + pb_y[k] * gs_1[k];

        t_4[k] = f_1 * fs_0[k]
                 + pb_z[k] * gs_2[k];

        t_5[k] = f_2 * fs_1[k]
                 + pb_y[k] * gs_3[k];

        t_6[k] = pb_z[k] * gs_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, pa_x, pa_y, pb_y, pb_z, fs_2, fp_3, fp_5, \
                         fp_11, gs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_y[k] * fp_3[k];

        t_8[k] = pb_y[k] * gs_4[k];

        t_9[k] = f_2 * fs_2[k]
                 + pb_z[k] * gs_4[k];

        t_10[k] = pa_x[k] * fp_5[k];

        t_11[k] = pa_x[k] * fp_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pb_x, pb_y, pb_z, fs_5, fs_6, \
                         gs_7, gs_8, gs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_x[k] * gs_7[k];

        t_13[k] = f_0 * fs_5[k]
                  + pb_y[k] * gs_7[k];

        t_14[k] = pb_z[k] * gs_7[k];

        t_15[k] = f_1 * fs_5[k]
                  + pb_z[k] * gs_8[k];

        t_16[k] = pb_x[k] * gs_9[k];

        t_17[k] = f_2 * fs_6[k]
                  + pb_z[k] * gs_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_y, pb_x, pb_y, pb_z, fs_8, fp_11, \
                         gs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pa_y[k] * fp_11[k];

        t_19[k] = pb_x[k] * gs_11[k];

        t_20[k] = pb_y[k] * gs_11[k];

        t_21[k] = f_0 * fs_8[k]
                  + pb_z[k] * gs_11[k];
    }
}

auto
compute_prim_gp_overlap_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fs, const size_t fp, const size_t gs,
                          const size_t ncols, const double p) -> void
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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_7 = buffer.data(fs + 7);

    const auto *fp_8 = buffer.data(fp + 8);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_10 = buffer.data(gs + 10);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pb_x, pb_y, pb_z, fs_0, fs_1, \
                         gs_0, gs_1, gs_2, gs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pb_x[k] * gs_0[k];

        t_1[k] = pb_y[k] * gs_0[k];

        t_2[k] = pb_z[k] * gs_0[k];

        t_3[k] = f_1 * fs_0[k]
                 + pb_z[k] * gs_1[k];

        t_4[k] = pb_z[k] * gs_2[k];

        t_5[k] = pb_y[k] * gs_3[k];

        t_6[k] = f_2 * fs_1[k]
                 + pb_z[k] * gs_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, pb_x, pb_y, pb_z, fs_4, fs_5, gs_6, \
                         gs_7, gs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_x[k] * gs_6[k];

        t_8[k] = f_0 * fs_4[k]
                 + pb_y[k] * gs_6[k];

        t_9[k] = pb_z[k] * gs_6[k];

        t_10[k] = f_1 * fs_4[k]
                  + pb_z[k] * gs_7[k];

        t_11[k] = pb_x[k] * gs_8[k];

        t_12[k] = f_2 * fs_5[k]
                  + pb_z[k] * gs_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pb_x, pb_y, pb_z, fs_7, fp_8, \
                         gs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * fp_8[k];

        t_14[k] = pb_x[k] * gs_10[k];

        t_15[k] = pb_y[k] * gs_10[k];

        t_16[k] = f_0 * fs_7[k]
                  + pb_z[k] * gs_10[k];
    }
}

auto
compute_prim_gp_overlap_8(CSimdMatrix &buffer, const size_t target, const size_t pb,
                          const size_t fs, const size_t gs, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_5 = buffer.data(fs + 5);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_11 = buffer.data(gs + 11);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fs_0, fs_3, fs_5, gs_0, gs_7, \
                         gs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pb_x[k] * gs_0[k];

        t_1[k] = f_0 * fs_3[k]
                 + pb_y[k] * gs_7[k];

        t_2[k] = f_0 * fs_5[k]
                 + pb_z[k] * gs_11[k];
    }
}

auto
compute_prim_gp_overlap_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fs, const size_t fp, const size_t gs,
                          const size_t ncols, const double p) -> void
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);
    const auto *fs_8 = buffer.data(fs + 8);

    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_3 = buffer.data(fp + 3);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_10 = buffer.data(gs + 10);
    const auto *gs_11 = buffer.data(gs + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fs_0, fs_1, fs_2, gs_0, \
                         gs_1, gs_2, gs_3, gs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pb_x[k] * gs_0[k];

        t_1[k] = f_1 * fs_0[k]
                 + pb_y[k] * gs_1[k];

        t_2[k] = f_1 * fs_0[k]
                 + pb_z[k] * gs_2[k];

        t_3[k] = f_2 * fs_1[k]
                 + pb_y[k] * gs_3[k];

        t_4[k] = f_2 * fs_2[k]
                 + pb_z[k] * gs_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pb_x, pb_y, pb_z, fs_5, fp_1, fp_3, \
                         gs_7, gs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_x[k] * fp_1[k];

        t_6[k] = pa_x[k] * fp_3[k];

        t_7[k] = pb_x[k] * gs_7[k];

        t_8[k] = f_0 * fs_5[k]
                 + pb_y[k] * gs_7[k];

        t_9[k] = f_1 * fs_5[k]
                 + pb_z[k] * gs_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, pb_y, pb_z, fs_6, fs_7, \
                         fs_8, fp_3, gs_9, gs_10, gs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * fs_7[k]
                  + pb_y[k] * gs_9[k];

        t_11[k] = f_2 * fs_6[k]
                  + pb_z[k] * gs_9[k];

        t_12[k] = f_1 * fs_8[k]
                  + pb_y[k] * gs_10[k];

        t_13[k] = pa_y[k] * fp_3[k];

        t_14[k] = pb_x[k] * gs_11[k];
    }

#pragma omp simd aligned(t_15, t_16, pb_y, pb_z, fs_8, gs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_y[k] * gs_11[k];

        t_16[k] = f_0 * fs_8[k]
                  + pb_z[k] * gs_11[k];
    }
}

auto
compute_prim_gp_overlap_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                           const size_t pb, const size_t fs, const size_t fp, const size_t gs,
                           const size_t ncols, const double p) -> void
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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_7 = buffer.data(fs + 7);

    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_8 = buffer.data(fp + 8);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_11 = buffer.data(gs + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_y, pb_x, pb_y, pb_z, fs_0, fs_1, \
                         fp_1, gs_0, gs_2, gs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pb_x[k] * gs_0[k];

        t_1[k] = pb_y[k] * gs_0[k];

        t_2[k] = pb_z[k] * gs_0[k];

        t_3[k] = f_1 * fs_0[k]
                 + pb_z[k] * gs_2[k];

        t_4[k] = pa_y[k] * fp_1[k];

        t_5[k] = f_2 * fs_1[k]
                 + pb_z[k] * gs_4[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pb_x, pb_y, pb_z, fs_4, fs_5, gs_7, \
                         gs_8, gs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_x[k] * gs_7[k];

        t_7[k] = f_0 * fs_4[k]
                 + pb_y[k] * gs_7[k];

        t_8[k] = pb_z[k] * gs_7[k];

        t_9[k] = f_1 * fs_4[k]
                 + pb_z[k] * gs_8[k];

        t_10[k] = pb_x[k] * gs_9[k];

        t_11[k] = f_2 * fs_5[k]
                  + pb_z[k] * gs_9[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_y, pb_x, pb_y, pb_z, fs_7, fp_8, \
                         gs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_y[k] * fp_8[k];

        t_13[k] = pb_x[k] * gs_11[k];

        t_14[k] = pb_y[k] * gs_11[k];

        t_15[k] = f_0 * fs_7[k]
                  + pb_z[k] * gs_11[k];
    }
}

auto
compute_prim_gp_overlap_11(CSimdMatrix &buffer, const size_t target, const size_t pb,
                           const size_t fs, const size_t gs, const size_t ncols,
                           const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_5 = buffer.data(fs + 5);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_10 = buffer.data(gs + 10);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pb_x, pb_y, pb_z, fs_0, fs_3, \
                         gs_0, gs_2, gs_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + pb_x[k] * gs_0[k];

        t_1[k] = pb_y[k] * gs_0[k];

        t_2[k] = pb_z[k] * gs_0[k];

        t_3[k] = pb_z[k] * gs_2[k];

        t_4[k] = pb_x[k] * gs_6[k];

        t_5[k] = f_0 * fs_3[k]
                 + pb_y[k] * gs_6[k];

        t_6[k] = pb_z[k] * gs_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, fs_5, gs_8, \
                         gs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_x[k] * gs_8[k];

        t_8[k] = pb_x[k] * gs_10[k];

        t_9[k] = pb_y[k] * gs_10[k];

        t_10[k] = f_0 * fs_5[k]
                  + pb_z[k] * gs_10[k];
    }
}

}  // namespace simdovl
