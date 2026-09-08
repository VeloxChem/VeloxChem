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


#include "SimdElectronRepulsionVrrRecPH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_ph_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t sg, const size_t sh,
                                     const size_t pg, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 0.5 / p;

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
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);

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
    const auto *pg_24 = buffer.data(pg + 24);
    const auto *pg_25 = buffer.data(pg + 25);
    const auto *pg_26 = buffer.data(pg + 26);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pb_y, pb_z, sg_0, sg_3, sh_0, sh_1, \
                         pg_0, pg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg_0[k]
                 + pa_x[k] * sh_0[k];

        t_1[k] = pb_y[k] * pg_0[k];

        t_2[k] = pb_z[k] * pg_0[k];

        t_3[k] = f_1 * sg_3[k]
                 + pa_x[k] * sh_1[k];

        t_4[k] = pb_y[k] * pg_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pb_y, pb_z, sg_4, sg_5, sg_6, sh_2, \
                         sh_3, sh_4, pg_2, pg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sg_4[k]
                 + pa_x[k] * sh_2[k];

        t_6[k] = f_2 * sg_5[k]
                 + pa_x[k] * sh_3[k];

        t_7[k] = pb_z[k] * pg_2[k];

        t_8[k] = pb_y[k] * pg_3[k];

        t_9[k] = f_2 * sg_6[k]
                 + pa_x[k] * sh_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pb_x, pb_y, pb_z, sg_7, sg_9, sg_11, \
                         pg_4, pg_5, pg_6, pg_7, pg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * sg_7[k]
                  + pb_x[k] * pg_6[k];

        t_11[k] = pb_z[k] * pg_4[k];

        t_12[k] = f_3 * sg_9[k]
                  + pb_x[k] * pg_7[k];

        t_13[k] = pb_y[k] * pg_5[k];

        t_14[k] = f_3 * sg_11[k]
                  + pb_x[k] * pg_8[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, t_20, pa_x, pb_y, pb_z, sh_5, sh_6, \
                         sh_7, sh_8, pg_6, pg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_x[k] * sh_5[k];

        t_16[k] = pb_z[k] * pg_6[k];

        t_17[k] = pa_x[k] * sh_6[k];

        t_18[k] = pa_x[k] * sh_7[k];

        t_19[k] = pb_y[k] * pg_8[k];

        t_20[k] = pa_x[k] * sh_8[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, t_26, pa_y, pb_y, pb_z, sg_0, sg_1, \
                         sh_0, sh_1, sh_2, pg_9, pg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_y[k] * sh_0[k];

        t_22[k] = f_3 * sg_0[k]
                  + pb_y[k] * pg_9[k];

        t_23[k] = pb_z[k] * pg_9[k];

        t_24[k] = f_2 * sg_1[k]
                  + pa_y[k] * sh_1[k];

        t_25[k] = pb_z[k] * pg_10[k];

        t_26[k] = pa_y[k] * sh_2[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_y, pb_x, pb_y, pb_z, sg_3, sg_4, \
                         sh_3, sh_4, pg_11, pg_12, pg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * sg_3[k]
                  + pa_y[k] * sh_3[k];

        t_28[k] = pb_z[k] * pg_11[k];

        t_29[k] = f_3 * sg_4[k]
                  + pb_y[k] * pg_12[k];

        t_30[k] = pa_y[k] * sh_4[k];

        t_31[k] = pb_x[k] * pg_13[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, pa_y, pb_x, pb_z, sg_7, sh_5, \
                         pg_13, pg_14, pg_15, pg_16, pg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_x[k] * pg_14[k];

        t_33[k] = pb_x[k] * pg_15[k];

        t_34[k] = pb_x[k] * pg_16[k];

        t_35[k] = pb_x[k] * pg_17[k];

        t_36[k] = f_0 * sg_7[k]
                  + pa_y[k] * sh_5[k];

        t_37[k] = pb_z[k] * pg_13[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, sg_9, sg_10, sg_11, \
                         sh_0, sh_6, sh_7, sh_8, pg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_1 * sg_9[k]
                  + pa_y[k] * sh_6[k];

        t_39[k] = f_2 * sg_10[k]
                  + pa_y[k] * sh_7[k];

        t_40[k] = f_3 * sg_11[k]
                  + pb_y[k] * pg_17[k];

        t_41[k] = pa_y[k] * sh_8[k];

        t_42[k] = pa_z[k] * sh_0[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, t_48, pa_z, pb_y, pb_z, sg_0, sg_2, \
                         sh_1, sh_2, sh_3, pg_18, pg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * pg_18[k];

        t_44[k] = f_3 * sg_0[k]
                  + pb_z[k] * pg_18[k];

        t_45[k] = pa_z[k] * sh_1[k];

        t_46[k] = pb_y[k] * pg_19[k];

        t_47[k] = f_2 * sg_2[k]
                  + pa_z[k] * sh_2[k];

        t_48[k] = pa_z[k] * sh_3[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_z, pb_x, pb_y, pb_z, sg_3, sg_4, \
                         sh_4, pg_20, pg_21, pg_22, pg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * sg_3[k]
                  + pb_z[k] * pg_20[k];

        t_50[k] = pb_y[k] * pg_21[k];

        t_51[k] = f_1 * sg_4[k]
                  + pa_z[k] * sh_4[k];

        t_52[k] = pb_x[k] * pg_22[k];

        t_53[k] = pb_x[k] * pg_23[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pa_z, pb_x, pb_z, sg_7, sh_5, pg_22, \
                         pg_24, pg_25, pg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pb_x[k] * pg_24[k];

        t_55[k] = pb_x[k] * pg_25[k];

        t_56[k] = pb_x[k] * pg_26[k];

        t_57[k] = pa_z[k] * sh_5[k];

        t_58[k] = f_3 * sg_7[k]
                  + pb_z[k] * pg_22[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_z, pb_y, sg_8, sg_9, sg_11, sh_6, sh_7, \
                         sh_8, pg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_2 * sg_8[k]
                  + pa_z[k] * sh_6[k];

        t_60[k] = f_1 * sg_9[k]
                  + pa_z[k] * sh_7[k];

        t_61[k] = pb_y[k] * pg_26[k];

        t_62[k] = f_0 * sg_11[k]
                  + pa_z[k] * sh_8[k];
    }
}

auto
compute_prim_ph_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t sg, const size_t sh,
                                     const size_t pg, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 0.5 / p;

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
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_12 = buffer.data(sh + 12);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pb_z, sg_0, sg_3, sg_4, sg_5, sh_0, \
                         sh_1, sh_2, sh_3, pg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg_0[k]
                 + pa_x[k] * sh_0[k];

        t_1[k] = pb_z[k] * pg_0[k];

        t_2[k] = f_1 * sg_3[k]
                 + pa_x[k] * sh_1[k];

        t_3[k] = f_1 * sg_4[k]
                 + pa_x[k] * sh_2[k];

        t_4[k] = f_2 * sg_5[k]
                 + pa_x[k] * sh_3[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pb_x, sg_6, sg_7, sg_11, sh_4, sh_8, \
                         sh_9, pg_1, pg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_2 * sg_6[k]
                 + pa_x[k] * sh_4[k];

        t_6[k] = f_3 * sg_7[k]
                 + pb_x[k] * pg_1[k];

        t_7[k] = f_3 * sg_11[k]
                 + pb_x[k] * pg_2[k];

        t_8[k] = pa_x[k] * sh_8[k];

        t_9[k] = pa_x[k] * sh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pa_y, pb_y, sg_0, sg_1, sh_0, \
                         sh_1, sh_10, sh_12, pg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * sh_10[k];

        t_11[k] = pa_x[k] * sh_12[k];

        t_12[k] = pa_y[k] * sh_0[k];

        t_13[k] = f_3 * sg_0[k]
                  + pb_y[k] * pg_3[k];

        t_14[k] = f_2 * sg_1[k]
                  + pa_y[k] * sh_1[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pb_z, sg_3, sg_7, sh_2, sh_3, \
                         sh_4, sh_8, pg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * sh_2[k];

        t_16[k] = f_1 * sg_3[k]
                  + pa_y[k] * sh_3[k];

        t_17[k] = pa_y[k] * sh_4[k];

        t_18[k] = f_0 * sg_7[k]
                  + pa_y[k] * sh_8[k];

        t_19[k] = pb_z[k] * pg_4[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pa_z, pb_y, sg_9, sg_10, sg_11, \
                         sh_0, sh_9, sh_10, sh_12, pg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * sg_9[k]
                  + pa_y[k] * sh_9[k];

        t_21[k] = f_2 * sg_10[k]
                  + pa_y[k] * sh_10[k];

        t_22[k] = f_3 * sg_11[k]
                  + pb_y[k] * pg_5[k];

        t_23[k] = pa_y[k] * sh_12[k];

        t_24[k] = pa_z[k] * sh_0[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_z, pb_z, sg_0, sg_2, sg_4, sh_1, \
                         sh_2, sh_3, sh_4, pg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * sg_0[k]
                  + pb_z[k] * pg_6[k];

        t_26[k] = pa_z[k] * sh_1[k];

        t_27[k] = f_2 * sg_2[k]
                  + pa_z[k] * sh_2[k];

        t_28[k] = pa_z[k] * sh_3[k];

        t_29[k] = f_1 * sg_4[k]
                  + pa_z[k] * sh_4[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pa_z, pb_y, pb_z, sg_7, sg_8, sg_9, \
                         sh_8, sh_9, sh_10, pg_7, pg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_z[k] * sh_8[k];

        t_31[k] = f_3 * sg_7[k]
                  + pb_z[k] * pg_7[k];

        t_32[k] = f_2 * sg_8[k]
                  + pa_z[k] * sh_9[k];

        t_33[k] = f_1 * sg_9[k]
                  + pa_z[k] * sh_10[k];

        t_34[k] = pb_y[k] * pg_8[k];
    }

#pragma omp simd aligned(t_35, pa_z, sg_11, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * sg_11[k]
                  + pa_z[k] * sh_12[k];
    }
}

auto
compute_prim_ph_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t sg, const size_t sh,
                                     const size_t pg, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 0.5 / p;

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
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);

    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_11 = buffer.data(pg + 11);
    const auto *pg_14 = buffer.data(pg + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, sg_0, sg_3, sg_4, sg_5, sg_6, sh_0, \
                         sh_1, sh_2, sh_3, sh_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg_0[k]
                 + pa_x[k] * sh_0[k];

        t_1[k] = f_1 * sg_3[k]
                 + pa_x[k] * sh_1[k];

        t_2[k] = f_1 * sg_4[k]
                 + pa_x[k] * sh_2[k];

        t_3[k] = f_2 * sg_5[k]
                 + pa_x[k] * sh_3[k];

        t_4[k] = f_2 * sg_6[k]
                 + pa_x[k] * sh_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_y, sg_0, sg_1, sg_7, sg_11, sh_1, \
                         pg_3, pg_4, pg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * sg_7[k]
                 + pb_x[k] * pg_3[k];

        t_6[k] = f_3 * sg_11[k]
                 + pb_x[k] * pg_4[k];

        t_7[k] = f_3 * sg_0[k]
                 + pb_y[k] * pg_5[k];

        t_8[k] = f_2 * sg_1[k]
                 + pa_y[k] * sh_1[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_y, sg_3, sg_7, sg_9, sg_10, sh_3, sh_5, \
                         sh_6, sh_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_1 * sg_3[k]
                 + pa_y[k] * sh_3[k];

        t_10[k] = f_0 * sg_7[k]
                  + pa_y[k] * sh_5[k];

        t_11[k] = f_1 * sg_9[k]
                  + pa_y[k] * sh_6[k];

        t_12[k] = f_2 * sg_10[k]
                  + pa_y[k] * sh_7[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pa_z, pb_y, pb_z, sg_0, sg_11, sh_0, \
                         sh_8, pg_10, pg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * sg_11[k]
                  + pb_y[k] * pg_10[k];

        t_14[k] = pa_y[k] * sh_8[k];

        t_15[k] = pa_z[k] * sh_0[k];

        t_16[k] = f_3 * sg_0[k]
                  + pb_z[k] * pg_11[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_z, pb_z, sg_2, sg_4, sg_7, sg_8, \
                         sh_2, sh_4, sh_5, sh_6, pg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_2 * sg_2[k]
                  + pa_z[k] * sh_2[k];

        t_18[k] = f_1 * sg_4[k]
                  + pa_z[k] * sh_4[k];

        t_19[k] = pa_z[k] * sh_5[k];

        t_20[k] = f_3 * sg_7[k]
                  + pb_z[k] * pg_14[k];

        t_21[k] = f_2 * sg_8[k]
                  + pa_z[k] * sh_6[k];
    }

#pragma omp simd aligned(t_22, t_23, pa_z, sg_9, sg_11, sh_7, sh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_1 * sg_9[k]
                  + pa_z[k] * sh_7[k];

        t_23[k] = f_0 * sg_11[k]
                  + pa_z[k] * sh_8[k];
    }
}

auto
compute_prim_ph_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t sg, const size_t sh,
                                     const size_t pg, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 0.5 / p;

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
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_13 = buffer.data(sh + 13);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_16 = buffer.data(sh + 16);

    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_15 = buffer.data(pg + 15);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, sg_0, sg_3, sg_4, sg_5, sg_6, sh_0, \
                         sh_3, sh_4, sh_5, sh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg_0[k]
                 + pa_x[k] * sh_0[k];

        t_1[k] = f_1 * sg_3[k]
                 + pa_x[k] * sh_3[k];

        t_2[k] = f_1 * sg_4[k]
                 + pa_x[k] * sh_4[k];

        t_3[k] = f_2 * sg_5[k]
                 + pa_x[k] * sh_5[k];

        t_4[k] = f_2 * sg_6[k]
                 + pa_x[k] * sh_8[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, sg_1, sh_0, sh_3, sh_11, \
                         sh_13, sh_14, sh_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_x[k] * sh_11[k];

        t_6[k] = pa_x[k] * sh_13[k];

        t_7[k] = pa_x[k] * sh_14[k];

        t_8[k] = pa_x[k] * sh_16[k];

        t_9[k] = pa_y[k] * sh_0[k];

        t_10[k] = f_2 * sg_1[k]
                  + pa_y[k] * sh_3[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pa_y, pb_x, sg_3, sg_7, sh_4, \
                         sh_5, sh_8, sh_11, pg_6, pg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_y[k] * sh_4[k];

        t_12[k] = f_1 * sg_3[k]
                  + pa_y[k] * sh_5[k];

        t_13[k] = pa_y[k] * sh_8[k];

        t_14[k] = pb_x[k] * pg_6[k];

        t_15[k] = pb_x[k] * pg_7[k];

        t_16[k] = f_0 * sg_7[k]
                  + pa_y[k] * sh_11[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_y, pa_z, pb_y, sg_9, sg_10, sg_11, \
                         sh_0, sh_13, sh_14, sh_16, pg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * sg_9[k]
                  + pa_y[k] * sh_13[k];

        t_18[k] = f_2 * sg_10[k]
                  + pa_y[k] * sh_14[k];

        t_19[k] = f_3 * sg_11[k]
                  + pb_y[k] * pg_8[k];

        t_20[k] = pa_y[k] * sh_16[k];

        t_21[k] = pa_z[k] * sh_0[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pa_z, pb_y, pb_z, sg_0, sg_2, sh_3, \
                         sh_4, sh_5, pg_9, pg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * sg_0[k]
                  + pb_z[k] * pg_9[k];

        t_23[k] = pa_z[k] * sh_3[k];

        t_24[k] = f_2 * sg_2[k]
                  + pa_z[k] * sh_4[k];

        t_25[k] = pa_z[k] * sh_5[k];

        t_26[k] = pb_y[k] * pg_10[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_z, pb_x, pb_z, sg_4, sg_7, sh_8, \
                         sh_11, pg_12, pg_13, pg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * sg_4[k]
                  + pa_z[k] * sh_8[k];

        t_28[k] = pb_x[k] * pg_13[k];

        t_29[k] = pb_x[k] * pg_15[k];

        t_30[k] = pa_z[k] * sh_11[k];

        t_31[k] = f_3 * sg_7[k]
                  + pb_z[k] * pg_12[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_z, pb_y, sg_8, sg_9, sg_11, sh_13, sh_14, \
                         sh_16, pg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_2 * sg_8[k]
                  + pa_z[k] * sh_13[k];

        t_33[k] = f_1 * sg_9[k]
                  + pa_z[k] * sh_14[k];

        t_34[k] = pb_y[k] * pg_15[k];

        t_35[k] = f_0 * sg_11[k]
                  + pa_z[k] * sh_16[k];
    }
}

auto
compute_prim_ph_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t sg, const size_t sh,
                                     const size_t pg, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 0.5 / p;

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
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_12 = buffer.data(sh + 12);

    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, sg_0, sg_3, sg_4, sg_5, sg_6, sh_0, \
                         sh_1, sh_2, sh_3, sh_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg_0[k]
                 + pa_x[k] * sh_0[k];

        t_1[k] = f_1 * sg_3[k]
                 + pa_x[k] * sh_1[k];

        t_2[k] = f_1 * sg_4[k]
                 + pa_x[k] * sh_2[k];

        t_3[k] = f_2 * sg_5[k]
                 + pa_x[k] * sh_3[k];

        t_4[k] = f_2 * sg_6[k]
                 + pa_x[k] * sh_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, sg_1, sg_3, sg_7, sh_0, \
                         sh_1, sh_3, sh_8, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_x[k] * sh_8[k];

        t_6[k] = pa_x[k] * sh_12[k];

        t_7[k] = pa_y[k] * sh_0[k];

        t_8[k] = f_2 * sg_1[k]
                 + pa_y[k] * sh_1[k];

        t_9[k] = f_1 * sg_3[k]
                 + pa_y[k] * sh_3[k];

        t_10[k] = f_0 * sg_7[k]
                  + pa_y[k] * sh_8[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pa_z, pb_y, sg_9, sg_10, sg_11, \
                         sh_0, sh_9, sh_10, sh_12, pg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_1 * sg_9[k]
                  + pa_y[k] * sh_9[k];

        t_12[k] = f_2 * sg_10[k]
                  + pa_y[k] * sh_10[k];

        t_13[k] = f_3 * sg_11[k]
                  + pb_y[k] * pg_2[k];

        t_14[k] = pa_y[k] * sh_12[k];

        t_15[k] = pa_z[k] * sh_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_z, pb_z, sg_0, sg_2, sg_4, sg_7, \
                         sh_2, sh_4, sh_8, pg_3, pg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * sg_0[k]
                  + pb_z[k] * pg_3[k];

        t_17[k] = f_2 * sg_2[k]
                  + pa_z[k] * sh_2[k];

        t_18[k] = f_1 * sg_4[k]
                  + pa_z[k] * sh_4[k];

        t_19[k] = pa_z[k] * sh_8[k];

        t_20[k] = f_3 * sg_7[k]
                  + pb_z[k] * pg_4[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_z, sg_8, sg_9, sg_11, sh_9, sh_10, \
                         sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_2 * sg_8[k]
                  + pa_z[k] * sh_9[k];

        t_22[k] = f_1 * sg_9[k]
                  + pa_z[k] * sh_10[k];

        t_23[k] = f_0 * sg_11[k]
                  + pa_z[k] * sh_12[k];
    }
}

auto
compute_prim_ph_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t sg, const size_t sh, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 1.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, sg_0, sg_1, sg_3, sg_5, sh_0, sh_1, \
                         sh_3, sh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg_0[k]
                 + pa_x[k] * sh_0[k];

        t_1[k] = f_1 * sg_1[k]
                 + pa_y[k] * sh_1[k];

        t_2[k] = f_2 * sg_3[k]
                 + pa_y[k] * sh_3[k];

        t_3[k] = f_0 * sg_5[k]
                 + pa_y[k] * sh_5[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_z, sg_2, sg_4, sg_6, sg_7, sg_8, sh_2, \
                         sh_4, sh_6, sh_7, sh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * sg_2[k]
                 + pa_z[k] * sh_2[k];

        t_5[k] = f_2 * sg_4[k]
                 + pa_z[k] * sh_4[k];

        t_6[k] = f_1 * sg_6[k]
                 + pa_z[k] * sh_6[k];

        t_7[k] = f_2 * sg_7[k]
                 + pa_z[k] * sh_7[k];

        t_8[k] = f_0 * sg_8[k]
                 + pa_z[k] * sh_8[k];
    }
}

auto
compute_prim_ph_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t sg, const size_t sh,
                                     const size_t pg, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 0.5 / p;

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

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);

    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_11 = buffer.data(pg + 11);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, sg_0, sg_3, sg_4, sg_5, sg_6, sh_0, \
                         sh_3, sh_4, sh_5, sh_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg_0[k]
                 + pa_x[k] * sh_0[k];

        t_1[k] = f_1 * sg_3[k]
                 + pa_x[k] * sh_3[k];

        t_2[k] = f_1 * sg_4[k]
                 + pa_x[k] * sh_4[k];

        t_3[k] = f_2 * sg_5[k]
                 + pa_x[k] * sh_5[k];

        t_4[k] = f_2 * sg_6[k]
                 + pa_x[k] * sh_7[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, sg_1, sg_3, sh_0, sh_3, sh_5, \
                         sh_8, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pa_x[k] * sh_8[k];

        t_6[k] = pa_x[k] * sh_12[k];

        t_7[k] = pa_y[k] * sh_0[k];

        t_8[k] = f_2 * sg_1[k]
                 + pa_y[k] * sh_3[k];

        t_9[k] = f_1 * sg_3[k]
                 + pa_y[k] * sh_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pb_x, sg_7, sg_9, sg_10, sh_8, sh_10, \
                         sh_11, pg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_x[k] * pg_3[k];

        t_11[k] = f_0 * sg_7[k]
                  + pa_y[k] * sh_8[k];

        t_12[k] = f_1 * sg_9[k]
                  + pa_y[k] * sh_10[k];

        t_13[k] = f_2 * sg_10[k]
                  + pa_y[k] * sh_11[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pa_z, pb_y, pb_z, sg_0, sg_11, sh_0, \
                         sh_12, pg_4, pg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * sg_11[k]
                  + pb_y[k] * pg_4[k];

        t_15[k] = pa_y[k] * sh_12[k];

        t_16[k] = pa_z[k] * sh_0[k];

        t_17[k] = f_3 * sg_0[k]
                  + pb_z[k] * pg_5[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_z, pb_x, pb_y, sg_2, sg_4, sh_4, \
                         sh_7, sh_8, pg_6, pg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * sg_2[k]
                  + pa_z[k] * sh_4[k];

        t_19[k] = pb_y[k] * pg_6[k];

        t_20[k] = f_1 * sg_4[k]
                  + pa_z[k] * sh_7[k];

        t_21[k] = pb_x[k] * pg_11[k];

        t_22[k] = pa_z[k] * sh_8[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pb_y, pb_z, sg_7, sg_8, sg_9, sh_10, \
                         sh_11, pg_8, pg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * sg_7[k]
                  + pb_z[k] * pg_8[k];

        t_24[k] = f_2 * sg_8[k]
                  + pa_z[k] * sh_10[k];

        t_25[k] = f_1 * sg_9[k]
                  + pa_z[k] * sh_11[k];

        t_26[k] = pb_y[k] * pg_11[k];
    }

#pragma omp simd aligned(t_27, pa_z, sg_11, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * sg_11[k]
                  + pa_z[k] * sh_12[k];
    }
}

auto
compute_prim_ph_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t sg, const size_t sh, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 1.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_11 = buffer.data(sg + 11);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_12 = buffer.data(sh + 12);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, sg_0, sg_1, sg_3, sg_7, sh_0, sh_1, \
                         sh_3, sh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg_0[k]
                 + pa_x[k] * sh_0[k];

        t_1[k] = f_1 * sg_1[k]
                 + pa_y[k] * sh_1[k];

        t_2[k] = f_2 * sg_3[k]
                 + pa_y[k] * sh_3[k];

        t_3[k] = f_0 * sg_7[k]
                 + pa_y[k] * sh_8[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_z, sg_2, sg_4, sg_8, sg_9, sg_11, sh_2, \
                         sh_4, sh_9, sh_10, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * sg_2[k]
                 + pa_z[k] * sh_2[k];

        t_5[k] = f_2 * sg_4[k]
                 + pa_z[k] * sh_4[k];

        t_6[k] = f_1 * sg_8[k]
                 + pa_z[k] * sh_9[k];

        t_7[k] = f_2 * sg_9[k]
                 + pa_z[k] * sh_10[k];

        t_8[k] = f_0 * sg_11[k]
                 + pa_z[k] * sh_12[k];
    }
}

auto
compute_prim_ph_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t sg, const size_t sh, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, sg_0, sg_1, sg_2, sh_0, sh_1, \
                         sh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg_0[k]
                 + pa_x[k] * sh_0[k];

        t_1[k] = f_0 * sg_1[k]
                 + pa_y[k] * sh_1[k];

        t_2[k] = f_0 * sg_2[k]
                 + pa_z[k] * sh_2[k];
    }
}

}  // namespace simdt2ceri
