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


#include "SimdKineticEnergyVrrRecPH.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_ph_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t sg, const size_t sh,
                                 const size_t ph_s, const size_t pg, const size_t ncols,
                                 const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 1.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;

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

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);
    const auto *sg_12 = buffer.data(sg + 12);
    const auto *sg_13 = buffer.data(sg + 13);
    const auto *sg_14 = buffer.data(sg + 14);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_15 = buffer.data(sh + 15);
    const auto *sh_16 = buffer.data(sh + 16);
    const auto *sh_17 = buffer.data(sh + 17);
    const auto *sh_18 = buffer.data(sh + 18);
    const auto *sh_20 = buffer.data(sh + 20);

    const auto *ph_s_0 = buffer.data(ph_s + 0);
    const auto *ph_s_1 = buffer.data(ph_s + 1);
    const auto *ph_s_2 = buffer.data(ph_s + 2);
    const auto *ph_s_3 = buffer.data(ph_s + 3);
    const auto *ph_s_4 = buffer.data(ph_s + 4);
    const auto *ph_s_5 = buffer.data(ph_s + 5);
    const auto *ph_s_6 = buffer.data(ph_s + 6);
    const auto *ph_s_7 = buffer.data(ph_s + 7);
    const auto *ph_s_8 = buffer.data(ph_s + 8);
    const auto *ph_s_9 = buffer.data(ph_s + 9);
    const auto *ph_s_10 = buffer.data(ph_s + 10);
    const auto *ph_s_11 = buffer.data(ph_s + 11);
    const auto *ph_s_12 = buffer.data(ph_s + 12);
    const auto *ph_s_13 = buffer.data(ph_s + 13);
    const auto *ph_s_14 = buffer.data(ph_s + 14);
    const auto *ph_s_15 = buffer.data(ph_s + 15);
    const auto *ph_s_16 = buffer.data(ph_s + 16);
    const auto *ph_s_17 = buffer.data(ph_s + 17);
    const auto *ph_s_18 = buffer.data(ph_s + 18);
    const auto *ph_s_19 = buffer.data(ph_s + 19);
    const auto *ph_s_20 = buffer.data(ph_s + 20);
    const auto *ph_s_21 = buffer.data(ph_s + 21);
    const auto *ph_s_22 = buffer.data(ph_s + 22);
    const auto *ph_s_23 = buffer.data(ph_s + 23);
    const auto *ph_s_24 = buffer.data(ph_s + 24);
    const auto *ph_s_25 = buffer.data(ph_s + 25);
    const auto *ph_s_26 = buffer.data(ph_s + 26);
    const auto *ph_s_27 = buffer.data(ph_s + 27);
    const auto *ph_s_28 = buffer.data(ph_s + 28);
    const auto *ph_s_29 = buffer.data(ph_s + 29);
    const auto *ph_s_30 = buffer.data(ph_s + 30);
    const auto *ph_s_31 = buffer.data(ph_s + 31);
    const auto *ph_s_32 = buffer.data(ph_s + 32);
    const auto *ph_s_33 = buffer.data(ph_s + 33);
    const auto *ph_s_34 = buffer.data(ph_s + 34);
    const auto *ph_s_35 = buffer.data(ph_s + 35);
    const auto *ph_s_36 = buffer.data(ph_s + 36);
    const auto *ph_s_37 = buffer.data(ph_s + 37);
    const auto *ph_s_38 = buffer.data(ph_s + 38);
    const auto *ph_s_39 = buffer.data(ph_s + 39);
    const auto *ph_s_40 = buffer.data(ph_s + 40);
    const auto *ph_s_41 = buffer.data(ph_s + 41);
    const auto *ph_s_42 = buffer.data(ph_s + 42);
    const auto *ph_s_43 = buffer.data(ph_s + 43);
    const auto *ph_s_44 = buffer.data(ph_s + 44);
    const auto *ph_s_45 = buffer.data(ph_s + 45);
    const auto *ph_s_46 = buffer.data(ph_s + 46);
    const auto *ph_s_47 = buffer.data(ph_s + 47);
    const auto *ph_s_48 = buffer.data(ph_s + 48);
    const auto *ph_s_49 = buffer.data(ph_s + 49);
    const auto *ph_s_50 = buffer.data(ph_s + 50);
    const auto *ph_s_51 = buffer.data(ph_s + 51);
    const auto *ph_s_52 = buffer.data(ph_s + 52);
    const auto *ph_s_53 = buffer.data(ph_s + 53);
    const auto *ph_s_54 = buffer.data(ph_s + 54);
    const auto *ph_s_55 = buffer.data(ph_s + 55);
    const auto *ph_s_56 = buffer.data(ph_s + 56);
    const auto *ph_s_57 = buffer.data(ph_s + 57);
    const auto *ph_s_58 = buffer.data(ph_s + 58);
    const auto *ph_s_59 = buffer.data(ph_s + 59);
    const auto *ph_s_60 = buffer.data(ph_s + 60);
    const auto *ph_s_61 = buffer.data(ph_s + 61);
    const auto *ph_s_62 = buffer.data(ph_s + 62);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_14 = buffer.data(pg + 14);
    const auto *pg_15 = buffer.data(pg + 15);
    const auto *pg_16 = buffer.data(pg + 16);
    const auto *pg_18 = buffer.data(pg + 18);
    const auto *pg_25 = buffer.data(pg + 25);
    const auto *pg_26 = buffer.data(pg + 26);
    const auto *pg_27 = buffer.data(pg + 27);
    const auto *pg_28 = buffer.data(pg + 28);
    const auto *pg_29 = buffer.data(pg + 29);
    const auto *pg_30 = buffer.data(pg + 30);
    const auto *pg_32 = buffer.data(pg + 32);
    const auto *pg_35 = buffer.data(pg + 35);
    const auto *pg_40 = buffer.data(pg + 40);
    const auto *pg_41 = buffer.data(pg + 41);
    const auto *pg_42 = buffer.data(pg + 42);
    const auto *pg_43 = buffer.data(pg + 43);
    const auto *pg_44 = buffer.data(pg + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pb_y, pb_z, sg_0, sg_3, sh_0, sh_3, ph_s_0, \
                         ph_s_1, ph_s_2, ph_s_3, pg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg_0[k]
                 + pa_x[k] * sh_0[k]
                 + f_1 * ph_s_0[k];

        t_1[k] = f_1 * ph_s_1[k]
                 + pb_y[k] * pg_0[k];

        t_2[k] = f_1 * ph_s_2[k]
                 + pb_z[k] * pg_0[k];

        t_3[k] = f_2 * sg_3[k]
                 + pa_x[k] * sh_3[k]
                 + f_1 * ph_s_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pb_y, sg_5, sg_6, sh_5, sh_6, ph_s_4, ph_s_5, \
                         ph_s_6, pg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * ph_s_4[k]
                 + pb_y[k] * pg_2[k];

        t_5[k] = f_2 * sg_5[k]
                 + pa_x[k] * sh_5[k]
                 + f_1 * ph_s_5[k];

        t_6[k] = f_3 * sg_6[k]
                 + pa_x[k] * sh_6[k]
                 + f_1 * ph_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_y, pb_z, sg_9, sh_9, ph_s_7, ph_s_8, ph_s_9, \
                         pg_3, pg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * ph_s_7[k]
                 + pb_z[k] * pg_3[k];

        t_8[k] = f_1 * ph_s_8[k]
                 + pb_y[k] * pg_5[k];

        t_9[k] = f_3 * sg_9[k]
                 + pa_x[k] * sh_9[k]
                 + f_1 * ph_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_x, pb_z, sg_10, sg_12, ph_s_10, ph_s_11, \
                         ph_s_12, pg_6, pg_10, pg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_4 * sg_10[k]
                  + f_1 * ph_s_10[k]
                  + pb_x[k] * pg_10[k];

        t_11[k] = f_1 * ph_s_11[k]
                  + pb_z[k] * pg_6[k];

        t_12[k] = f_4 * sg_12[k]
                  + f_1 * ph_s_12[k]
                  + pb_x[k] * pg_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pb_x, pb_y, sg_14, sh_15, ph_s_13, ph_s_14, \
                         ph_s_15, pg_9, pg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * ph_s_13[k]
                  + pb_y[k] * pg_9[k];

        t_14[k] = f_4 * sg_14[k]
                  + f_1 * ph_s_14[k]
                  + pb_x[k] * pg_14[k];

        t_15[k] = pa_x[k] * sh_15[k]
                  + f_1 * ph_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, pb_z, sh_17, sh_18, ph_s_16, \
                         ph_s_17, ph_s_18, ph_s_19, pg_10, pg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * ph_s_16[k]
                  + pb_z[k] * pg_10[k];

        t_17[k] = pa_x[k] * sh_17[k]
                  + f_1 * ph_s_17[k];

        t_18[k] = pa_x[k] * sh_18[k]
                  + f_1 * ph_s_18[k];

        t_19[k] = f_1 * ph_s_19[k]
                  + pb_y[k] * pg_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pa_y, pb_z, sg_0, sh_0, sh_1, sh_20, \
                         ph_s_20, ph_s_21, ph_s_22, ph_s_23, pg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_x[k] * sh_20[k]
                  + f_1 * ph_s_20[k];

        t_21[k] = pa_y[k] * sh_0[k]
                  + f_1 * ph_s_21[k];

        t_22[k] = f_4 * sg_0[k]
                  + pa_y[k] * sh_1[k]
                  + f_1 * ph_s_22[k];

        t_23[k] = f_1 * ph_s_23[k]
                  + pb_z[k] * pg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pb_z, sg_1, sg_3, sh_3, sh_5, sh_6, \
                         ph_s_24, ph_s_25, ph_s_26, ph_s_27, pg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_3 * sg_1[k]
                  + pa_y[k] * sh_3[k]
                  + f_1 * ph_s_24[k];

        t_25[k] = f_1 * ph_s_25[k]
                  + pb_z[k] * pg_16[k];

        t_26[k] = pa_y[k] * sh_5[k]
                  + f_1 * ph_s_26[k];

        t_27[k] = f_2 * sg_3[k]
                  + pa_y[k] * sh_6[k]
                  + f_1 * ph_s_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pb_x, pb_z, sg_5, sh_8, sh_9, ph_s_28, \
                         ph_s_29, ph_s_30, ph_s_31, pg_18, pg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * ph_s_28[k]
                  + pb_z[k] * pg_18[k];

        t_29[k] = f_4 * sg_5[k]
                  + pa_y[k] * sh_8[k]
                  + f_1 * ph_s_29[k];

        t_30[k] = pa_y[k] * sh_9[k]
                  + f_1 * ph_s_30[k];

        t_31[k] = f_1 * ph_s_31[k]
                  + pb_x[k] * pg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_x, ph_s_32, ph_s_33, ph_s_34, ph_s_35, \
                         pg_26, pg_27, pg_28, pg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * ph_s_32[k]
                  + pb_x[k] * pg_26[k];

        t_33[k] = f_1 * ph_s_33[k]
                  + pb_x[k] * pg_27[k];

        t_34[k] = f_1 * ph_s_34[k]
                  + pb_x[k] * pg_28[k];

        t_35[k] = f_1 * ph_s_35[k]
                  + pb_x[k] * pg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_y, pb_z, sg_10, sg_12, sh_15, sh_17, ph_s_36, \
                         ph_s_37, ph_s_38, pg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * sg_10[k]
                  + pa_y[k] * sh_15[k]
                  + f_1 * ph_s_36[k];

        t_37[k] = f_1 * ph_s_37[k]
                  + pb_z[k] * pg_25[k];

        t_38[k] = f_2 * sg_12[k]
                  + pa_y[k] * sh_17[k]
                  + f_1 * ph_s_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_y, sg_13, sg_14, sh_18, sh_20, ph_s_39, \
                         ph_s_40, ph_s_41, pg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_3 * sg_13[k]
                  + pa_y[k] * sh_18[k]
                  + f_1 * ph_s_39[k];

        t_40[k] = f_4 * sg_14[k]
                  + f_1 * ph_s_40[k]
                  + pb_y[k] * pg_29[k];

        t_41[k] = pa_y[k] * sh_20[k]
                  + f_1 * ph_s_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_z, pb_y, sg_0, sh_0, sh_2, sh_3, ph_s_42, \
                         ph_s_43, ph_s_44, ph_s_45, pg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * sh_0[k]
                  + f_1 * ph_s_42[k];

        t_43[k] = f_1 * ph_s_43[k]
                  + pb_y[k] * pg_30[k];

        t_44[k] = f_4 * sg_0[k]
                  + pa_z[k] * sh_2[k]
                  + f_1 * ph_s_44[k];

        t_45[k] = pa_z[k] * sh_3[k]
                  + f_1 * ph_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_z, pb_y, sg_2, sg_3, sh_5, sh_6, sh_7, \
                         ph_s_46, ph_s_47, ph_s_48, ph_s_49, pg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_1 * ph_s_46[k]
                  + pb_y[k] * pg_32[k];

        t_47[k] = f_3 * sg_2[k]
                  + pa_z[k] * sh_5[k]
                  + f_1 * ph_s_47[k];

        t_48[k] = pa_z[k] * sh_6[k]
                  + f_1 * ph_s_48[k];

        t_49[k] = f_4 * sg_3[k]
                  + pa_z[k] * sh_7[k]
                  + f_1 * ph_s_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_z, pb_x, pb_y, sg_5, sh_9, ph_s_50, \
                         ph_s_51, ph_s_52, ph_s_53, pg_35, pg_40, \
                         pg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_1 * ph_s_50[k]
                  + pb_y[k] * pg_35[k];

        t_51[k] = f_2 * sg_5[k]
                  + pa_z[k] * sh_9[k]
                  + f_1 * ph_s_51[k];

        t_52[k] = f_1 * ph_s_52[k]
                  + pb_x[k] * pg_40[k];

        t_53[k] = f_1 * ph_s_53[k]
                  + pb_x[k] * pg_41[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_z, pb_x, sh_15, ph_s_54, ph_s_55, ph_s_56, \
                         ph_s_57, pg_42, pg_43, pg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * ph_s_54[k]
                  + pb_x[k] * pg_42[k];

        t_55[k] = f_1 * ph_s_55[k]
                  + pb_x[k] * pg_43[k];

        t_56[k] = f_1 * ph_s_56[k]
                  + pb_x[k] * pg_44[k];

        t_57[k] = pa_z[k] * sh_15[k]
                  + f_1 * ph_s_57[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_z, sg_10, sg_11, sg_12, sh_16, sh_17, sh_18, \
                         ph_s_58, ph_s_59, ph_s_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_4 * sg_10[k]
                  + pa_z[k] * sh_16[k]
                  + f_1 * ph_s_58[k];

        t_59[k] = f_3 * sg_11[k]
                  + pa_z[k] * sh_17[k]
                  + f_1 * ph_s_59[k];

        t_60[k] = f_2 * sg_12[k]
                  + pa_z[k] * sh_18[k]
                  + f_1 * ph_s_60[k];
    }

#pragma omp simd aligned(t_61, t_62, pa_z, pb_y, sg_14, sh_20, ph_s_61, ph_s_62, \
                         pg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_1 * ph_s_61[k]
                  + pb_y[k] * pg_44[k];

        t_62[k] = f_0 * sg_14[k]
                  + pa_z[k] * sh_20[k]
                  + f_1 * ph_s_62[k];
    }
}

}  // namespace simdkin
