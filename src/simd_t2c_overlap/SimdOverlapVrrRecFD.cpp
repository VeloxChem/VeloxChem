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


#include "SimdOverlapVrrRecFD.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_fd_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pd, const size_t dp, const size_t dd,
                          const size_t fs, const size_t fp, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
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

    const auto *pd_9 = buffer.data(pd + 9);
    const auto *pd_17 = buffer.data(pd + 17);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_8 = buffer.data(dp + 8);
    const auto *dp_9 = buffer.data(dp + 9);
    const auto *dp_10 = buffer.data(dp + 10);
    const auto *dp_14 = buffer.data(dp + 14);
    const auto *dp_15 = buffer.data(dp + 15);
    const auto *dp_16 = buffer.data(dp + 16);
    const auto *dp_17 = buffer.data(dp + 17);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_14 = buffer.data(dd + 14);
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

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_9 = buffer.data(fs + 9);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_15 = buffer.data(fp + 15);
    const auto *fp_17 = buffer.data(fp + 17);
    const auto *fp_18 = buffer.data(fp + 18);
    const auto *fp_19 = buffer.data(fp + 19);
    const auto *fp_20 = buffer.data(fp + 20);
    const auto *fp_22 = buffer.data(fp + 22);
    const auto *fp_23 = buffer.data(fp + 23);
    const auto *fp_25 = buffer.data(fp + 25);
    const auto *fp_26 = buffer.data(fp + 26);
    const auto *fp_27 = buffer.data(fp + 27);
    const auto *fp_28 = buffer.data(fp + 28);
    const auto *fp_29 = buffer.data(fp + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, dp_0, fs_0, fp_0, \
                         fp_1, fp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dp_0[k]
                 + f_1 * fs_0[k]
                 + pb_x[k] * fp_0[k];

        t_1[k] = pb_y[k] * fp_0[k];

        t_2[k] = pb_z[k] * fp_0[k];

        t_3[k] = f_1 * fs_0[k]
                 + pb_y[k] * fp_1[k];

        t_4[k] = pb_y[k] * fp_2[k];

        t_5[k] = f_1 * fs_0[k]
                 + pb_z[k] * fp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pb_x, pb_z, pd_9, dp_4, dd_0, \
                         dd_9, fp_3, fp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * dd_0[k];

        t_7[k] = f_2 * dp_4[k]
                 + pb_x[k] * fp_4[k];

        t_8[k] = pb_z[k] * fp_3[k];

        t_9[k] = f_1 * pd_9[k]
                 + pa_x[k] * dd_9[k];

        t_10[k] = pb_z[k] * fp_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pa_y, pa_z, pb_x, pb_y, dp_8, \
                         dd_0, dd_3, dd_5, fp_6, fp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_y[k] * dd_5[k];

        t_12[k] = pa_z[k] * dd_0[k];

        t_13[k] = pb_y[k] * fp_6[k];

        t_14[k] = f_2 * dp_8[k]
                  + pb_x[k] * fp_8[k];

        t_15[k] = pa_z[k] * dd_3[k];

        t_16[k] = pb_y[k] * fp_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_x, pb_x, pb_z, pd_17, dp_9, dp_10, \
                         dd_17, dd_18, dd_21, fp_9, fp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * pd_17[k]
                  + pa_x[k] * dd_17[k];

        t_18[k] = f_2 * dp_9[k]
                  + pa_x[k] * dd_18[k];

        t_19[k] = f_1 * dp_10[k]
                  + pb_x[k] * fp_10[k];

        t_20[k] = pb_z[k] * fp_9[k];

        t_21[k] = pa_x[k] * dd_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, t_27, pa_x, pa_y, pa_z, pb_z, dd_7, \
                         dd_12, dd_14, dd_23, dd_27, fp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_z[k] * fp_10[k];

        t_23[k] = pa_x[k] * dd_23[k];

        t_24[k] = pa_y[k] * dd_12[k];

        t_25[k] = pa_z[k] * dd_7[k];

        t_26[k] = pa_y[k] * dd_14[k];

        t_27[k] = pa_x[k] * dd_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_x, pb_x, pb_y, dp_15, dp_17, dd_28, \
                         dd_29, dd_30, fp_15, fp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_x[k] * dd_28[k];

        t_29[k] = pa_x[k] * dd_29[k];

        t_30[k] = f_2 * dp_15[k]
                  + pa_x[k] * dd_30[k];

        t_31[k] = pb_y[k] * fp_15[k];

        t_32[k] = f_1 * dp_17[k]
                  + pb_x[k] * fp_17[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pa_x, pb_x, pb_y, dd_33, dd_35, \
                         fs_6, fp_17, fp_18, fp_19, fp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_x[k] * dd_33[k];

        t_34[k] = pb_y[k] * fp_17[k];

        t_35[k] = pa_x[k] * dd_35[k];

        t_36[k] = f_1 * fs_6[k]
                  + pb_x[k] * fp_18[k];

        t_37[k] = pb_x[k] * fp_19[k];

        t_38[k] = pb_x[k] * fp_20[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_z, pb_x, pb_y, pb_z, dp_10, dd_18, \
                         fs_6, fp_19, fp_20, fp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * dp_10[k]
                  + f_1 * fs_6[k]
                  + pb_y[k] * fp_19[k];

        t_40[k] = pb_z[k] * fp_19[k];

        t_41[k] = f_1 * fs_6[k]
                  + pb_z[k] * fp_20[k];

        t_42[k] = pa_z[k] * dd_18[k];

        t_43[k] = pb_x[k] * fp_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pa_y, pa_z, pb_x, pb_y, pd_17, dp_14, \
                         dd_21, dd_29, dd_30, fp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * fp_23[k];

        t_45[k] = pa_z[k] * dd_21[k];

        t_46[k] = f_2 * dp_14[k]
                  + pb_y[k] * fp_23[k];

        t_47[k] = f_1 * pd_17[k]
                  + pa_y[k] * dd_29[k];

        t_48[k] = pa_y[k] * dd_30[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_y, pb_x, pb_y, dp_16, dp_17, dd_33, \
                         dd_35, fp_25, fp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pb_x[k] * fp_25[k];

        t_50[k] = pb_x[k] * fp_26[k];

        t_51[k] = f_2 * dp_16[k]
                  + pa_y[k] * dd_33[k];

        t_52[k] = f_1 * dp_17[k]
                  + pb_y[k] * fp_26[k];

        t_53[k] = pa_y[k] * dd_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, pb_x, pb_y, pb_z, dp_17, fs_9, \
                         fp_27, fp_28, fp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * fs_9[k]
                  + pb_x[k] * fp_27[k];

        t_55[k] = pb_x[k] * fp_28[k];

        t_56[k] = pb_x[k] * fp_29[k];

        t_57[k] = f_1 * fs_9[k]
                  + pb_y[k] * fp_28[k];

        t_58[k] = pb_y[k] * fp_29[k];

        t_59[k] = f_0 * dp_17[k]
                  + f_1 * fs_9[k]
                  + pb_z[k] * fp_29[k];
    }
}

}  // namespace simdovl
