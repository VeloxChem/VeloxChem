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


#include "SimdKineticEnergyVrrRecGP.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_gp_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fs, const size_t fp,
                                 const size_t gp_s, const size_t gs, const size_t ncols,
                                 const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;

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
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);
    const auto *fs_8 = buffer.data(fs + 8);
    const auto *fs_9 = buffer.data(fs + 9);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_15 = buffer.data(fp + 15);
    const auto *fp_19 = buffer.data(fp + 19);
    const auto *fp_22 = buffer.data(fp + 22);
    const auto *fp_23 = buffer.data(fp + 23);
    const auto *fp_25 = buffer.data(fp + 25);
    const auto *fp_26 = buffer.data(fp + 26);
    const auto *fp_29 = buffer.data(fp + 29);

    const auto *gp_s_0 = buffer.data(gp_s + 0);
    const auto *gp_s_1 = buffer.data(gp_s + 1);
    const auto *gp_s_2 = buffer.data(gp_s + 2);
    const auto *gp_s_3 = buffer.data(gp_s + 3);
    const auto *gp_s_4 = buffer.data(gp_s + 4);
    const auto *gp_s_5 = buffer.data(gp_s + 5);
    const auto *gp_s_6 = buffer.data(gp_s + 6);
    const auto *gp_s_7 = buffer.data(gp_s + 7);
    const auto *gp_s_8 = buffer.data(gp_s + 8);
    const auto *gp_s_9 = buffer.data(gp_s + 9);
    const auto *gp_s_10 = buffer.data(gp_s + 10);
    const auto *gp_s_11 = buffer.data(gp_s + 11);
    const auto *gp_s_12 = buffer.data(gp_s + 12);
    const auto *gp_s_13 = buffer.data(gp_s + 13);
    const auto *gp_s_14 = buffer.data(gp_s + 14);
    const auto *gp_s_15 = buffer.data(gp_s + 15);
    const auto *gp_s_16 = buffer.data(gp_s + 16);
    const auto *gp_s_17 = buffer.data(gp_s + 17);
    const auto *gp_s_18 = buffer.data(gp_s + 18);
    const auto *gp_s_19 = buffer.data(gp_s + 19);
    const auto *gp_s_20 = buffer.data(gp_s + 20);
    const auto *gp_s_21 = buffer.data(gp_s + 21);
    const auto *gp_s_22 = buffer.data(gp_s + 22);
    const auto *gp_s_23 = buffer.data(gp_s + 23);
    const auto *gp_s_24 = buffer.data(gp_s + 24);
    const auto *gp_s_25 = buffer.data(gp_s + 25);
    const auto *gp_s_26 = buffer.data(gp_s + 26);
    const auto *gp_s_27 = buffer.data(gp_s + 27);
    const auto *gp_s_28 = buffer.data(gp_s + 28);
    const auto *gp_s_29 = buffer.data(gp_s + 29);
    const auto *gp_s_30 = buffer.data(gp_s + 30);
    const auto *gp_s_31 = buffer.data(gp_s + 31);
    const auto *gp_s_32 = buffer.data(gp_s + 32);
    const auto *gp_s_33 = buffer.data(gp_s + 33);
    const auto *gp_s_34 = buffer.data(gp_s + 34);
    const auto *gp_s_35 = buffer.data(gp_s + 35);
    const auto *gp_s_36 = buffer.data(gp_s + 36);
    const auto *gp_s_37 = buffer.data(gp_s + 37);
    const auto *gp_s_38 = buffer.data(gp_s + 38);
    const auto *gp_s_39 = buffer.data(gp_s + 39);
    const auto *gp_s_40 = buffer.data(gp_s + 40);
    const auto *gp_s_41 = buffer.data(gp_s + 41);
    const auto *gp_s_42 = buffer.data(gp_s + 42);
    const auto *gp_s_43 = buffer.data(gp_s + 43);
    const auto *gp_s_44 = buffer.data(gp_s + 44);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_10 = buffer.data(gs + 10);
    const auto *gs_11 = buffer.data(gs + 11);
    const auto *gs_12 = buffer.data(gs + 12);
    const auto *gs_13 = buffer.data(gs + 13);
    const auto *gs_14 = buffer.data(gs + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, fs_0, fp_0, gp_s_0, \
                         gp_s_1, gp_s_2, gp_s_3, gs_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fs_0[k]
                 + f_1 * gp_s_0[k]
                 + pb_x[k] * gs_0[k];

        t_1[k] = f_1 * gp_s_1[k]
                 + pb_y[k] * gs_0[k];

        t_2[k] = f_1 * gp_s_2[k]
                 + pb_z[k] * gs_0[k];

        t_3[k] = pa_y[k] * fp_0[k]
                 + f_1 * gp_s_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_z, pb_y, pb_z, fs_0, fp_0, gp_s_4, gp_s_5, \
                         gp_s_6, gp_s_7, gs_1, gs_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * fs_0[k]
                 + f_1 * gp_s_4[k]
                 + pb_y[k] * gs_1[k];

        t_5[k] = f_1 * gp_s_5[k]
                 + pb_z[k] * gs_1[k];

        t_6[k] = pa_z[k] * fp_0[k]
                 + f_1 * gp_s_6[k];

        t_7[k] = f_1 * gp_s_7[k]
                 + pb_y[k] * gs_2[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pb_y, pb_z, fs_0, fs_1, fs_3, gp_s_8, \
                         gp_s_9, gp_s_10, gp_s_11, gs_2, gs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * fs_0[k]
                 + f_1 * gp_s_8[k]
                 + pb_z[k] * gs_2[k];

        t_9[k] = f_3 * fs_3[k]
                 + f_1 * gp_s_9[k]
                 + pb_x[k] * gs_3[k];

        t_10[k] = f_3 * fs_1[k]
                  + f_1 * gp_s_10[k]
                  + pb_y[k] * gs_3[k];

        t_11[k] = f_1 * gp_s_11[k]
                  + pb_z[k] * gs_3[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_y, pa_z, pb_x, fs_5, fp_4, fp_6, fp_8, \
                         gp_s_12, gp_s_13, gp_s_14, gp_s_15, gs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_y[k] * fp_6[k]
                  + f_1 * gp_s_12[k];

        t_13[k] = pa_z[k] * fp_4[k]
                  + f_1 * gp_s_13[k];

        t_14[k] = pa_y[k] * fp_8[k]
                  + f_1 * gp_s_14[k];

        t_15[k] = f_3 * fs_5[k]
                  + f_1 * gp_s_15[k]
                  + pb_x[k] * gs_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, pb_y, pb_z, fs_2, fs_6, gp_s_16, gp_s_17, \
                         gp_s_18, gs_5, gs_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * gp_s_16[k]
                  + pb_y[k] * gs_5[k];

        t_17[k] = f_3 * fs_2[k]
                  + f_1 * gp_s_17[k]
                  + pb_z[k] * gs_5[k];

        t_18[k] = f_2 * fs_6[k]
                  + f_1 * gp_s_18[k]
                  + pb_x[k] * gs_6[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pa_z, pb_z, fp_9, fp_19, fp_22, \
                         gp_s_19, gp_s_20, gp_s_21, gp_s_22, gs_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_x[k] * fp_19[k]
                  + f_1 * gp_s_19[k];

        t_20[k] = f_1 * gp_s_20[k]
                  + pb_z[k] * gs_6[k];

        t_21[k] = pa_z[k] * fp_9[k]
                  + f_1 * gp_s_21[k];

        t_22[k] = pa_x[k] * fp_22[k]
                  + f_1 * gp_s_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pa_y, fp_15, fp_23, fp_25, fp_26, \
                         gp_s_23, gp_s_24, gp_s_25, gp_s_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_x[k] * fp_23[k]
                  + f_1 * gp_s_23[k];

        t_24[k] = pa_y[k] * fp_15[k]
                  + f_1 * gp_s_24[k];

        t_25[k] = pa_x[k] * fp_25[k]
                  + f_1 * gp_s_25[k];

        t_26[k] = pa_x[k] * fp_26[k]
                  + f_1 * gp_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pb_x, pb_y, fs_9, fp_29, gp_s_27, \
                         gp_s_28, gp_s_29, gp_s_30, gs_9, gs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_2 * fs_9[k]
                  + f_1 * gp_s_27[k]
                  + pb_x[k] * gs_9[k];

        t_28[k] = f_1 * gp_s_28[k]
                  + pb_y[k] * gs_9[k];

        t_29[k] = pa_x[k] * fp_29[k]
                  + f_1 * gp_s_29[k];

        t_30[k] = f_1 * gp_s_30[k]
                  + pb_x[k] * gs_10[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_z, pb_x, pb_y, pb_z, fs_6, fp_19, gp_s_31, \
                         gp_s_32, gp_s_33, gp_s_34, gs_10, gs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * fs_6[k]
                  + f_1 * gp_s_31[k]
                  + pb_y[k] * gs_10[k];

        t_32[k] = f_1 * gp_s_32[k]
                  + pb_z[k] * gs_10[k];

        t_33[k] = f_1 * gp_s_33[k]
                  + pb_x[k] * gs_11[k];

        t_34[k] = pa_z[k] * fp_19[k]
                  + f_1 * gp_s_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, fs_6, fs_7, fs_8, gp_s_35, \
                         gp_s_36, gp_s_37, gp_s_38, gs_11, gs_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_2 * fs_6[k]
                  + f_1 * gp_s_35[k]
                  + pb_z[k] * gs_11[k];

        t_36[k] = f_1 * gp_s_36[k]
                  + pb_x[k] * gs_12[k];

        t_37[k] = f_3 * fs_8[k]
                  + f_1 * gp_s_37[k]
                  + pb_y[k] * gs_12[k];

        t_38[k] = f_3 * fs_7[k]
                  + f_1 * gp_s_38[k]
                  + pb_z[k] * gs_12[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_y, pb_x, pb_y, fs_9, fp_29, gp_s_39, \
                         gp_s_40, gp_s_41, gp_s_42, gs_13, gs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_1 * gp_s_39[k]
                  + pb_x[k] * gs_13[k];

        t_40[k] = f_2 * fs_9[k]
                  + f_1 * gp_s_40[k]
                  + pb_y[k] * gs_13[k];

        t_41[k] = pa_y[k] * fp_29[k]
                  + f_1 * gp_s_41[k];

        t_42[k] = f_1 * gp_s_42[k]
                  + pb_x[k] * gs_14[k];
    }

#pragma omp simd aligned(t_43, t_44, pb_y, pb_z, fs_9, gp_s_43, gp_s_44, \
                         gs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * gp_s_43[k]
                  + pb_y[k] * gs_14[k];

        t_44[k] = f_0 * fs_9[k]
                  + f_1 * gp_s_44[k]
                  + pb_z[k] * gs_14[k];
    }
}

}  // namespace simdkin
