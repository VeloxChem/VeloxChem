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


#include "SimdKineticEnergyVrrRecHP.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_hp_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t gs, const size_t gp,
                                 const size_t hp_s, const size_t hs, const size_t ncols,
                                 const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

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

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_27 = buffer.data(gp + 27);
    const auto *gp_31 = buffer.data(gp + 31);
    const auto *gp_34 = buffer.data(gp + 34);
    const auto *gp_35 = buffer.data(gp + 35);
    const auto *gp_37 = buffer.data(gp + 37);
    const auto *gp_38 = buffer.data(gp + 38);
    const auto *gp_40 = buffer.data(gp + 40);
    const auto *gp_41 = buffer.data(gp + 41);
    const auto *gp_44 = buffer.data(gp + 44);

    const auto *hp_s_0 = buffer.data(hp_s + 0);
    const auto *hp_s_1 = buffer.data(hp_s + 1);
    const auto *hp_s_2 = buffer.data(hp_s + 2);
    const auto *hp_s_3 = buffer.data(hp_s + 3);
    const auto *hp_s_4 = buffer.data(hp_s + 4);
    const auto *hp_s_5 = buffer.data(hp_s + 5);
    const auto *hp_s_6 = buffer.data(hp_s + 6);
    const auto *hp_s_7 = buffer.data(hp_s + 7);
    const auto *hp_s_8 = buffer.data(hp_s + 8);
    const auto *hp_s_9 = buffer.data(hp_s + 9);
    const auto *hp_s_10 = buffer.data(hp_s + 10);
    const auto *hp_s_11 = buffer.data(hp_s + 11);
    const auto *hp_s_12 = buffer.data(hp_s + 12);
    const auto *hp_s_13 = buffer.data(hp_s + 13);
    const auto *hp_s_14 = buffer.data(hp_s + 14);
    const auto *hp_s_15 = buffer.data(hp_s + 15);
    const auto *hp_s_16 = buffer.data(hp_s + 16);
    const auto *hp_s_17 = buffer.data(hp_s + 17);
    const auto *hp_s_18 = buffer.data(hp_s + 18);
    const auto *hp_s_19 = buffer.data(hp_s + 19);
    const auto *hp_s_20 = buffer.data(hp_s + 20);
    const auto *hp_s_21 = buffer.data(hp_s + 21);
    const auto *hp_s_22 = buffer.data(hp_s + 22);
    const auto *hp_s_23 = buffer.data(hp_s + 23);
    const auto *hp_s_24 = buffer.data(hp_s + 24);
    const auto *hp_s_25 = buffer.data(hp_s + 25);
    const auto *hp_s_26 = buffer.data(hp_s + 26);
    const auto *hp_s_27 = buffer.data(hp_s + 27);
    const auto *hp_s_28 = buffer.data(hp_s + 28);
    const auto *hp_s_29 = buffer.data(hp_s + 29);
    const auto *hp_s_30 = buffer.data(hp_s + 30);
    const auto *hp_s_31 = buffer.data(hp_s + 31);
    const auto *hp_s_32 = buffer.data(hp_s + 32);
    const auto *hp_s_33 = buffer.data(hp_s + 33);
    const auto *hp_s_34 = buffer.data(hp_s + 34);
    const auto *hp_s_35 = buffer.data(hp_s + 35);
    const auto *hp_s_36 = buffer.data(hp_s + 36);
    const auto *hp_s_37 = buffer.data(hp_s + 37);
    const auto *hp_s_38 = buffer.data(hp_s + 38);
    const auto *hp_s_39 = buffer.data(hp_s + 39);
    const auto *hp_s_40 = buffer.data(hp_s + 40);
    const auto *hp_s_41 = buffer.data(hp_s + 41);
    const auto *hp_s_42 = buffer.data(hp_s + 42);
    const auto *hp_s_43 = buffer.data(hp_s + 43);
    const auto *hp_s_44 = buffer.data(hp_s + 44);
    const auto *hp_s_45 = buffer.data(hp_s + 45);
    const auto *hp_s_46 = buffer.data(hp_s + 46);
    const auto *hp_s_47 = buffer.data(hp_s + 47);
    const auto *hp_s_48 = buffer.data(hp_s + 48);
    const auto *hp_s_49 = buffer.data(hp_s + 49);
    const auto *hp_s_50 = buffer.data(hp_s + 50);
    const auto *hp_s_51 = buffer.data(hp_s + 51);
    const auto *hp_s_52 = buffer.data(hp_s + 52);
    const auto *hp_s_53 = buffer.data(hp_s + 53);
    const auto *hp_s_54 = buffer.data(hp_s + 54);
    const auto *hp_s_55 = buffer.data(hp_s + 55);
    const auto *hp_s_56 = buffer.data(hp_s + 56);
    const auto *hp_s_57 = buffer.data(hp_s + 57);
    const auto *hp_s_58 = buffer.data(hp_s + 58);
    const auto *hp_s_59 = buffer.data(hp_s + 59);
    const auto *hp_s_60 = buffer.data(hp_s + 60);
    const auto *hp_s_61 = buffer.data(hp_s + 61);
    const auto *hp_s_62 = buffer.data(hp_s + 62);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);
    const auto *hs_16 = buffer.data(hs + 16);
    const auto *hs_17 = buffer.data(hs + 17);
    const auto *hs_18 = buffer.data(hs + 18);
    const auto *hs_19 = buffer.data(hs + 19);
    const auto *hs_20 = buffer.data(hs + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, gs_0, gp_0, hp_s_0, \
                         hp_s_1, hp_s_2, hp_s_3, hs_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gs_0[k]
                 + f_1 * hp_s_0[k]
                 + pb_x[k] * hs_0[k];

        t_1[k] = f_1 * hp_s_1[k]
                 + pb_y[k] * hs_0[k];

        t_2[k] = f_1 * hp_s_2[k]
                 + pb_z[k] * hs_0[k];

        t_3[k] = pa_y[k] * gp_0[k]
                 + f_1 * hp_s_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_z, pb_y, pb_z, gs_0, gp_0, hp_s_4, hp_s_5, \
                         hp_s_6, hp_s_7, hs_1, hs_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * gs_0[k]
                 + f_1 * hp_s_4[k]
                 + pb_y[k] * hs_1[k];

        t_5[k] = f_1 * hp_s_5[k]
                 + pb_z[k] * hs_1[k];

        t_6[k] = pa_z[k] * gp_0[k]
                 + f_1 * hp_s_6[k];

        t_7[k] = f_1 * hp_s_7[k]
                 + pb_y[k] * hs_2[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pb_y, pb_z, gs_0, gs_1, gs_3, hp_s_8, \
                         hp_s_9, hp_s_10, hp_s_11, hs_2, hs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * gs_0[k]
                 + f_1 * hp_s_8[k]
                 + pb_z[k] * hs_2[k];

        t_9[k] = f_3 * gs_3[k]
                 + f_1 * hp_s_9[k]
                 + pb_x[k] * hs_3[k];

        t_10[k] = f_4 * gs_1[k]
                  + f_1 * hp_s_10[k]
                  + pb_y[k] * hs_3[k];

        t_11[k] = f_1 * hp_s_11[k]
                  + pb_z[k] * hs_3[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_y, pa_z, pb_x, gs_5, gp_4, gp_6, gp_8, \
                         hp_s_12, hp_s_13, hp_s_14, hp_s_15, hs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_y[k] * gp_6[k]
                  + f_1 * hp_s_12[k];

        t_13[k] = pa_z[k] * gp_4[k]
                  + f_1 * hp_s_13[k];

        t_14[k] = pa_y[k] * gp_8[k]
                  + f_1 * hp_s_14[k];

        t_15[k] = f_3 * gs_5[k]
                  + f_1 * hp_s_15[k]
                  + pb_x[k] * hs_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_x, pb_y, pb_z, gs_2, gs_3, gs_6, hp_s_16, \
                         hp_s_17, hp_s_18, hp_s_19, hs_5, hs_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * hp_s_16[k]
                  + pb_y[k] * hs_5[k];

        t_17[k] = f_4 * gs_2[k]
                  + f_1 * hp_s_17[k]
                  + pb_z[k] * hs_5[k];

        t_18[k] = f_4 * gs_6[k]
                  + f_1 * hp_s_18[k]
                  + pb_x[k] * hs_6[k];

        t_19[k] = f_3 * gs_3[k]
                  + f_1 * hp_s_19[k]
                  + pb_y[k] * hs_6[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_z, pb_z, gs_3, gp_9, gp_10, hp_s_20, \
                         hp_s_21, hp_s_22, hp_s_23, hs_6, hs_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * hp_s_20[k]
                  + pb_z[k] * hs_6[k];

        t_21[k] = pa_z[k] * gp_9[k]
                  + f_1 * hp_s_21[k];

        t_22[k] = pa_z[k] * gp_10[k]
                  + f_1 * hp_s_22[k];

        t_23[k] = f_2 * gs_3[k]
                  + f_1 * hp_s_23[k]
                  + pb_z[k] * hs_7[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_y, pb_y, gs_5, gp_15, gp_17, hp_s_24, hp_s_25, \
                         hp_s_26, hs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * gp_15[k]
                  + f_1 * hp_s_24[k];

        t_25[k] = f_2 * gs_5[k]
                  + f_1 * hp_s_25[k]
                  + pb_y[k] * hs_8[k];

        t_26[k] = pa_y[k] * gp_17[k]
                  + f_1 * hp_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_x, pb_y, pb_z, gs_5, gs_9, gs_10, hp_s_27, \
                         hp_s_28, hp_s_29, hp_s_30, hs_9, hs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_4 * gs_9[k]
                  + f_1 * hp_s_27[k]
                  + pb_x[k] * hs_9[k];

        t_28[k] = f_1 * hp_s_28[k]
                  + pb_y[k] * hs_9[k];

        t_29[k] = f_3 * gs_5[k]
                  + f_1 * hp_s_29[k]
                  + pb_z[k] * hs_9[k];

        t_30[k] = f_2 * gs_10[k]
                  + f_1 * hp_s_30[k]
                  + pb_x[k] * hs_10[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pa_z, pb_z, gp_18, gp_31, gp_34, \
                         hp_s_31, hp_s_32, hp_s_33, hp_s_34, hs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pa_x[k] * gp_31[k]
                  + f_1 * hp_s_31[k];

        t_32[k] = f_1 * hp_s_32[k]
                  + pb_z[k] * hs_10[k];

        t_33[k] = pa_z[k] * gp_18[k]
                  + f_1 * hp_s_33[k];

        t_34[k] = pa_x[k] * gp_34[k]
                  + f_1 * hp_s_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_x, gs_12, gp_35, gp_37, gp_38, \
                         hp_s_35, hp_s_36, hp_s_37, hp_s_38, hs_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_x[k] * gp_35[k]
                  + f_1 * hp_s_35[k];

        t_36[k] = f_2 * gs_12[k]
                  + f_1 * hp_s_36[k]
                  + pb_x[k] * hs_12[k];

        t_37[k] = pa_x[k] * gp_37[k]
                  + f_1 * hp_s_37[k];

        t_38[k] = pa_x[k] * gp_38[k]
                  + f_1 * hp_s_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pa_y, pb_x, gs_14, gp_27, gp_40, gp_41, \
                         hp_s_39, hp_s_40, hp_s_41, hp_s_42, hs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pa_y[k] * gp_27[k]
                  + f_1 * hp_s_39[k];

        t_40[k] = pa_x[k] * gp_40[k]
                  + f_1 * hp_s_40[k];

        t_41[k] = pa_x[k] * gp_41[k]
                  + f_1 * hp_s_41[k];

        t_42[k] = f_2 * gs_14[k]
                  + f_1 * hp_s_42[k]
                  + pb_x[k] * hs_14[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pb_x, pb_y, gs_10, gp_44, hp_s_43, \
                         hp_s_44, hp_s_45, hp_s_46, hs_14, hs_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * hp_s_43[k]
                  + pb_y[k] * hs_14[k];

        t_44[k] = pa_x[k] * gp_44[k]
                  + f_1 * hp_s_44[k];

        t_45[k] = f_1 * hp_s_45[k]
                  + pb_x[k] * hs_15[k];

        t_46[k] = f_0 * gs_10[k]
                  + f_1 * hp_s_46[k]
                  + pb_y[k] * hs_15[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_z, pb_x, pb_z, gs_10, gp_31, hp_s_47, \
                         hp_s_48, hp_s_49, hp_s_50, hs_15, hs_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_1 * hp_s_47[k]
                  + pb_z[k] * hs_15[k];

        t_48[k] = f_1 * hp_s_48[k]
                  + pb_x[k] * hs_16[k];

        t_49[k] = pa_z[k] * gp_31[k]
                  + f_1 * hp_s_49[k];

        t_50[k] = f_2 * gs_10[k]
                  + f_1 * hp_s_50[k]
                  + pb_z[k] * hs_16[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_x, pb_y, pb_z, gs_11, gs_12, hp_s_51, \
                         hp_s_52, hp_s_53, hp_s_54, hs_17, hs_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * hp_s_51[k]
                  + pb_x[k] * hs_17[k];

        t_52[k] = f_3 * gs_12[k]
                  + f_1 * hp_s_52[k]
                  + pb_y[k] * hs_17[k];

        t_53[k] = f_4 * gs_11[k]
                  + f_1 * hp_s_53[k]
                  + pb_z[k] * hs_17[k];

        t_54[k] = f_1 * hp_s_54[k]
                  + pb_x[k] * hs_18[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_x, pb_y, pb_z, gs_12, gs_13, gs_14, \
                         hp_s_55, hp_s_56, hp_s_57, hp_s_58, hs_18, \
                         hs_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_4 * gs_13[k]
                  + f_1 * hp_s_55[k]
                  + pb_y[k] * hs_18[k];

        t_56[k] = f_3 * gs_12[k]
                  + f_1 * hp_s_56[k]
                  + pb_z[k] * hs_18[k];

        t_57[k] = f_1 * hp_s_57[k]
                  + pb_x[k] * hs_19[k];

        t_58[k] = f_2 * gs_14[k]
                  + f_1 * hp_s_58[k]
                  + pb_y[k] * hs_19[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pb_x, pb_y, pb_z, gs_14, gp_44, \
                         hp_s_59, hp_s_60, hp_s_61, hp_s_62, hs_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_y[k] * gp_44[k]
                  + f_1 * hp_s_59[k];

        t_60[k] = f_1 * hp_s_60[k]
                  + pb_x[k] * hs_20[k];

        t_61[k] = f_1 * hp_s_61[k]
                  + pb_y[k] * hs_20[k];

        t_62[k] = f_0 * gs_14[k]
                  + f_1 * hp_s_62[k]
                  + pb_z[k] * hs_20[k];
    }
}

}  // namespace simdkin
