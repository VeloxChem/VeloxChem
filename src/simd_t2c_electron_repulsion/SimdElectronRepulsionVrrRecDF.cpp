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


#include "SimdElectronRepulsionVrrRecDF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_df_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pd, const size_t pf,
                                     const size_t dp0, const size_t dp1, const size_t dd,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_7 = buffer.data(pd + 7);
    const auto *pd_8 = buffer.data(pd + 8);

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

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);

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
    const auto *dd_17 = buffer.data(dd + 17);
    const auto *dd_18 = buffer.data(dd + 18);
    const auto *dd_19 = buffer.data(dd + 19);
    const auto *dd_20 = buffer.data(dd + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pd_0, pd_1, dp0_0, dp1_0, \
                         dd_0, dd_1, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_y[k] * dd_0[k];

        t_2[k] = pb_z[k] * dd_0[k];

        t_3[k] = f_0 * pd_1[k]
                 + pb_x[k] * dd_2[k];

        t_4[k] = pb_y[k] * dd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, pd_2, dp0_1, dp0_2, dp1_1, \
                         dp1_2, dd_2, dd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * pd_2[k]
                 + pb_x[k] * dd_3[k];

        t_6[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_2[k];

        t_7[k] = pb_z[k] * dd_2[k];

        t_8[k] = pb_y[k] * dd_3[k];

        t_9[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, pb_y, pb_z, pd_0, pd_4, \
                         pf_0, dd_4, dd_5, dd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * pf_0[k];

        t_11[k] = f_3 * pd_0[k]
                  + pb_y[k] * dd_4[k];

        t_12[k] = pb_z[k] * dd_4[k];

        t_13[k] = f_3 * pd_4[k]
                  + pb_x[k] * dd_6[k];

        t_14[k] = pb_z[k] * dd_5[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, t_20, pa_x, pa_y, pa_z, pb_z, pf_0, \
                         pf_2, pf_4, pf_5, pf_6, dd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * pf_2[k];

        t_16[k] = pa_x[k] * pf_4[k];

        t_17[k] = pb_z[k] * dd_6[k];

        t_18[k] = pa_x[k] * pf_5[k];

        t_19[k] = pa_x[k] * pf_6[k];

        t_20[k] = pa_z[k] * pf_0[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_z, pb_x, pb_y, pb_z, pd_0, pd_8, \
                         pf_1, dd_7, dd_8, dd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pb_y[k] * dd_7[k];

        t_22[k] = f_3 * pd_0[k]
                  + pb_z[k] * dd_7[k];

        t_23[k] = pa_z[k] * pf_1[k];

        t_24[k] = pb_y[k] * dd_8[k];

        t_25[k] = f_3 * pd_8[k]
                  + pb_x[k] * dd_9[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_x, pb_x, pb_y, pf_9, pf_10, pf_11, \
                         dp0_3, dp1_3, dd_9, dd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_x[k] * pf_9[k];

        t_27[k] = pa_x[k] * pf_10[k];

        t_28[k] = pb_y[k] * dd_9[k];

        t_29[k] = pa_x[k] * pf_11[k];

        t_30[k] = f_1 * dp0_3[k]
                  - f_2 * dp1_3[k]
                  + pb_x[k] * dd_10[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pb_x, pb_y, pb_z, pd_3, dd_10, dd_11, \
                         dd_12, dd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * pd_3[k]
                  + pb_y[k] * dd_10[k];

        t_32[k] = pb_z[k] * dd_10[k];

        t_33[k] = pb_x[k] * dd_11[k];

        t_34[k] = pb_x[k] * dd_12[k];

        t_35[k] = pb_x[k] * dd_13[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pb_z, pd_4, pd_5, dp0_4, dp0_5, dp1_4, \
                         dp1_5, dd_11, dd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * pd_4[k]
                  + f_1 * dp0_4[k]
                  - f_2 * dp1_4[k]
                  + pb_y[k] * dd_11[k];

        t_37[k] = pb_z[k] * dd_11[k];

        t_38[k] = f_0 * pd_5[k]
                  + pb_y[k] * dd_13[k];

        t_39[k] = f_1 * dp0_5[k]
                  - f_2 * dp1_5[k]
                  + pb_z[k] * dd_13[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, pa_y, pa_z, pb_x, pf_3, pf_7, \
                         pf_8, dd_14, dd_15, dd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * pf_7[k];

        t_41[k] = pa_z[k] * pf_3[k];

        t_42[k] = pa_y[k] * pf_8[k];

        t_43[k] = pb_x[k] * dd_14[k];

        t_44[k] = pb_x[k] * dd_15[k];

        t_45[k] = pb_x[k] * dd_16[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pa_z, pb_y, pb_z, pd_4, pd_8, pf_4, \
                         pf_11, dd_14, dd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_z[k] * pf_4[k];

        t_47[k] = f_3 * pd_4[k]
                  + pb_z[k] * dd_14[k];

        t_48[k] = f_3 * pd_8[k]
                  + pb_y[k] * dd_16[k];

        t_49[k] = pa_y[k] * pf_11[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, t_55, pb_x, pb_y, pb_z, pd_6, dp0_6, \
                         dp1_6, dd_17, dd_18, dd_19, dd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_1 * dp0_6[k]
                  - f_2 * dp1_6[k]
                  + pb_x[k] * dd_17[k];

        t_51[k] = pb_y[k] * dd_17[k];

        t_52[k] = f_0 * pd_6[k]
                  + pb_z[k] * dd_17[k];

        t_53[k] = pb_x[k] * dd_18[k];

        t_54[k] = pb_x[k] * dd_19[k];

        t_55[k] = pb_x[k] * dd_20[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pb_z, pd_7, pd_8, dp0_7, dp0_8, dp1_7, \
                         dp1_8, dd_18, dd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * dp0_7[k]
                  - f_2 * dp1_7[k]
                  + pb_y[k] * dd_18[k];

        t_57[k] = f_0 * pd_7[k]
                  + pb_z[k] * dd_18[k];

        t_58[k] = pb_y[k] * dd_20[k];

        t_59[k] = f_0 * pd_8[k]
                  + f_1 * dp0_8[k]
                  - f_2 * dp1_8[k]
                  + pb_z[k] * dd_20[k];
    }
}

auto
compute_prim_df_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pd, const size_t pf,
                                     const size_t dp0, const size_t dp1, const size_t dd,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_7 = buffer.data(pd + 7);
    const auto *pd_8 = buffer.data(pd + 8);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pd_0, pd_1, pd_2, dp0_0, \
                         dp1_0, dd_0, dd_1, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_y[k] * dd_0[k];

        t_2[k] = pb_z[k] * dd_0[k];

        t_3[k] = f_0 * pd_1[k]
                 + pb_x[k] * dd_1[k];

        t_4[k] = f_0 * pd_2[k]
                 + pb_x[k] * dd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_y, pb_z, pf_0, dp0_1, dp0_2, dp1_1, \
                         dp1_2, dd_1, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_6[k] = pb_y[k] * dd_2[k];

        t_7[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];

        t_8[k] = pa_y[k] * pf_0[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_x, pa_z, pb_x, pb_y, pd_0, pd_4, \
                         pf_0, pf_4, pf_5, dd_3, dd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * pd_0[k]
                 + pb_y[k] * dd_3[k];

        t_10[k] = f_3 * pd_4[k]
                  + pb_x[k] * dd_4[k];

        t_11[k] = pa_x[k] * pf_4[k];

        t_12[k] = pa_x[k] * pf_5[k];

        t_13[k] = pa_z[k] * pf_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pb_x, pb_z, pd_0, pd_8, pf_10, pf_11, \
                         dd_5, dd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * pd_0[k]
                  + pb_z[k] * dd_5[k];

        t_15[k] = f_3 * pd_8[k]
                  + pb_x[k] * dd_6[k];

        t_16[k] = pa_x[k] * pf_10[k];

        t_17[k] = pa_x[k] * pf_11[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pb_x, pb_y, pd_3, pd_4, dp0_3, dp0_4, \
                         dp1_3, dp1_4, dd_7, dd_8, dd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_1 * dp0_3[k]
                  - f_2 * dp1_3[k]
                  + pb_x[k] * dd_7[k];

        t_19[k] = f_0 * pd_3[k]
                  + pb_y[k] * dd_7[k];

        t_20[k] = pb_x[k] * dd_8[k];

        t_21[k] = pb_x[k] * dd_9[k];

        t_22[k] = f_0 * pd_4[k]
                  + f_1 * dp0_4[k]
                  - f_2 * dp1_4[k]
                  + pb_y[k] * dd_8[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pa_y, pa_z, pb_y, pb_z, pd_5, pf_4, \
                         pf_8, dp0_5, dp1_5, dd_8, dd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pb_z[k] * dd_8[k];

        t_24[k] = f_0 * pd_5[k]
                  + pb_y[k] * dd_9[k];

        t_25[k] = f_1 * dp0_5[k]
                  - f_2 * dp1_5[k]
                  + pb_z[k] * dd_9[k];

        t_26[k] = pa_y[k] * pf_8[k];

        t_27[k] = pa_z[k] * pf_4[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pb_x, pb_y, pb_z, pd_4, pd_8, pf_11, \
                         dp0_6, dp1_6, dd_10, dd_11, dd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pd_4[k]
                  + pb_z[k] * dd_10[k];

        t_29[k] = f_3 * pd_8[k]
                  + pb_y[k] * dd_11[k];

        t_30[k] = pa_y[k] * pf_11[k];

        t_31[k] = f_1 * dp0_6[k]
                  - f_2 * dp1_6[k]
                  + pb_x[k] * dd_12[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, pb_x, pb_y, pb_z, pd_6, pd_7, \
                         dp0_7, dp1_7, dd_12, dd_13, dd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * pd_6[k]
                  + pb_z[k] * dd_12[k];

        t_33[k] = pb_x[k] * dd_13[k];

        t_34[k] = pb_x[k] * dd_14[k];

        t_35[k] = f_1 * dp0_7[k]
                  - f_2 * dp1_7[k]
                  + pb_y[k] * dd_13[k];

        t_36[k] = f_0 * pd_7[k]
                  + pb_z[k] * dd_13[k];

        t_37[k] = pb_y[k] * dd_14[k];
    }

#pragma omp simd aligned(t_38, pb_z, pd_8, dp0_8, dp1_8, dd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * pd_8[k]
                  + f_1 * dp0_8[k]
                  - f_2 * dp1_8[k]
                  + pb_z[k] * dd_14[k];
    }
}

auto
compute_prim_df_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pd, const size_t dp0, const size_t dp1,
                                     const size_t dd, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pd_0, dp0_0, dp0_1, dp1_0, \
                         dp1_1, dd_0, dd_1, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_y[k] * dd_0[k];

        t_2[k] = pb_z[k] * dd_0[k];

        t_3[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_4[k] = pb_y[k] * dd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_x, pb_z, dp0_2, dp0_3, dp1_2, dp1_3, dd_2, \
                         dd_3, dd_4, dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];

        t_6[k] = f_1 * dp0_3[k]
                 - f_2 * dp1_3[k]
                 + pb_x[k] * dd_3[k];

        t_7[k] = pb_x[k] * dd_4[k];

        t_8[k] = pb_x[k] * dd_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_y, pb_z, pd_1, dp0_4, dp0_5, dp1_4, dp1_5, dd_4, \
                         dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * pd_1[k]
                 + f_1 * dp0_4[k]
                 - f_2 * dp1_4[k]
                 + pb_y[k] * dd_4[k];

        t_10[k] = pb_z[k] * dd_4[k];

        t_11[k] = f_1 * dp0_5[k]
                  - f_2 * dp1_5[k]
                  + pb_z[k] * dd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_x, pb_y, dp0_6, dp0_7, dp1_6, dp1_7, \
                         dd_6, dd_7, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_1 * dp0_6[k]
                  - f_2 * dp1_6[k]
                  + pb_x[k] * dd_6[k];

        t_13[k] = pb_x[k] * dd_7[k];

        t_14[k] = pb_x[k] * dd_8[k];

        t_15[k] = f_1 * dp0_7[k]
                  - f_2 * dp1_7[k]
                  + pb_y[k] * dd_7[k];

        t_16[k] = pb_y[k] * dd_8[k];
    }

#pragma omp simd aligned(t_17, pb_z, pd_2, dp0_8, dp1_8, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * pd_2[k]
                  + f_1 * dp0_8[k]
                  - f_2 * dp1_8[k]
                  + pb_z[k] * dd_8[k];
    }
}

auto
compute_prim_df_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pd, const size_t pf,
                                     const size_t dp0, const size_t dp1, const size_t dd,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

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

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_5 = buffer.data(pf + 5);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pd_0, dp0_0, dp0_1, dp1_0, \
                         dp1_1, dd_0, dd_1, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_y[k] * dd_0[k];

        t_2[k] = pb_z[k] * dd_0[k];

        t_3[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_4[k] = pb_y[k] * dd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pa_z, pb_z, pf_0, pf_1, pf_5, \
                         dp0_2, dp1_2, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];

        t_6[k] = pa_y[k] * pf_0[k];

        t_7[k] = pa_x[k] * pf_1[k];

        t_8[k] = pa_z[k] * pf_0[k];

        t_9[k] = pa_x[k] * pf_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pb_x, pb_y, pb_z, pd_1, dp0_3, dp0_4, \
                         dp1_3, dp1_4, dd_3, dd_4, dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * dp0_3[k]
                  - f_2 * dp1_3[k]
                  + pb_x[k] * dd_3[k];

        t_11[k] = pb_x[k] * dd_4[k];

        t_12[k] = pb_x[k] * dd_5[k];

        t_13[k] = f_0 * pd_1[k]
                  + f_1 * dp0_4[k]
                  - f_2 * dp1_4[k]
                  + pb_y[k] * dd_4[k];

        t_14[k] = pb_z[k] * dd_4[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_y, pa_z, pb_x, pb_z, pf_1, pf_5, dp0_5, \
                         dp0_6, dp1_5, dp1_6, dd_5, dd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * dp0_5[k]
                  - f_2 * dp1_5[k]
                  + pb_z[k] * dd_5[k];

        t_16[k] = pa_z[k] * pf_1[k];

        t_17[k] = pa_y[k] * pf_5[k];

        t_18[k] = f_1 * dp0_6[k]
                  - f_2 * dp1_6[k]
                  + pb_x[k] * dd_6[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pb_x, pb_y, pb_z, pd_2, dp0_7, dp0_8, \
                         dp1_7, dp1_8, dd_7, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_x[k] * dd_7[k];

        t_20[k] = pb_x[k] * dd_8[k];

        t_21[k] = f_1 * dp0_7[k]
                  - f_2 * dp1_7[k]
                  + pb_y[k] * dd_7[k];

        t_22[k] = pb_y[k] * dd_8[k];

        t_23[k] = f_0 * pd_2[k]
                  + f_1 * dp0_8[k]
                  - f_2 * dp1_8[k]
                  + pb_z[k] * dd_8[k];
    }
}

auto
compute_prim_df_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pd, const size_t pf,
                                     const size_t dp0, const size_t dp1, const size_t dd,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_7 = buffer.data(pd + 7);
    const auto *pd_8 = buffer.data(pd + 8);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_9 = buffer.data(dp1 + 9);
    const auto *dp1_10 = buffer.data(dp1 + 10);
    const auto *dp1_11 = buffer.data(dp1 + 11);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pd_0, pd_1, pd_2, dp0_0, dp0_1, \
                         dp1_0, dp1_1, dd_0, dd_1, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = f_0 * pd_1[k]
                 + pb_x[k] * dd_1[k];

        t_2[k] = f_0 * pd_2[k]
                 + pb_x[k] * dd_2[k];

        t_3[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_z, pb_x, pb_y, pb_z, pd_0, pd_4, pf_0, dp0_2, \
                         dp1_2, dd_2, dd_3, dd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];

        t_5[k] = f_3 * pd_0[k]
                 + pb_y[k] * dd_3[k];

        t_6[k] = f_3 * pd_4[k]
                 + pb_x[k] * dd_4[k];

        t_7[k] = pa_z[k] * pf_0[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pb_y, pb_z, pd_0, pd_3, pd_8, dp0_3, \
                         dp1_5, dd_5, dd_6, dd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_3 * pd_0[k]
                 + pb_z[k] * dd_5[k];

        t_9[k] = f_3 * pd_8[k]
                 + pb_x[k] * dd_6[k];

        t_10[k] = f_1 * dp0_3[k]
                  - f_2 * dp1_5[k]
                  + pb_x[k] * dd_7[k];

        t_11[k] = f_0 * pd_3[k]
                  + pb_y[k] * dd_7[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_y, pb_z, pd_4, pd_5, pf_1, dp0_4, \
                         dp0_5, dp1_6, dp1_7, dd_8, dd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * pd_4[k]
                  + f_1 * dp0_4[k]
                  - f_2 * dp1_6[k]
                  + pb_y[k] * dd_8[k];

        t_13[k] = f_0 * pd_5[k]
                  + pb_y[k] * dd_9[k];

        t_14[k] = f_1 * dp0_5[k]
                  - f_2 * dp1_7[k]
                  + pb_z[k] * dd_9[k];

        t_15[k] = pa_z[k] * pf_1[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_y, pb_x, pb_y, pb_z, pd_4, pd_8, pf_2, \
                         dp0_6, dp1_9, dd_10, dd_12, dd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pd_4[k]
                  + pb_z[k] * dd_10[k];

        t_17[k] = f_3 * pd_8[k]
                  + pb_y[k] * dd_12[k];

        t_18[k] = pa_y[k] * pf_2[k];

        t_19[k] = f_1 * dp0_6[k]
                  - f_2 * dp1_9[k]
                  + pb_x[k] * dd_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pb_y, pb_z, pd_6, pd_7, pd_8, dp0_7, dp0_8, \
                         dp1_10, dp1_11, dd_13, dd_14, dd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * pd_6[k]
                  + pb_z[k] * dd_13[k];

        t_21[k] = f_1 * dp0_7[k]
                  - f_2 * dp1_10[k]
                  + pb_y[k] * dd_14[k];

        t_22[k] = f_0 * pd_7[k]
                  + pb_z[k] * dd_14[k];

        t_23[k] = f_0 * pd_8[k]
                  + f_1 * dp0_8[k]
                  - f_2 * dp1_11[k]
                  + pb_z[k] * dd_15[k];
    }
}

auto
compute_prim_df_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pd, const size_t pf,
                                     const size_t dp0, const size_t dp1, const size_t dd,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
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

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_12 = buffer.data(pf + 12);
    const auto *pf_14 = buffer.data(pf + 14);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pd_0, dp0_0, dp0_1, dp1_0, \
                         dp1_1, dd_0, dd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_y[k] * dd_0[k];

        t_2[k] = pb_z[k] * dd_0[k];

        t_3[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_4[k] = pb_z[k] * dd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pb_y, pb_z, pf_0, pf_5, pf_6, \
                         dp0_2, dp1_2, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_y[k] * dd_2[k];

        t_6[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];

        t_7[k] = pa_y[k] * pf_0[k];

        t_8[k] = pa_x[k] * pf_5[k];

        t_9[k] = pa_x[k] * pf_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_z, pb_z, pd_0, pf_0, pf_12, pf_14, \
                         dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * pf_0[k];

        t_11[k] = f_3 * pd_0[k]
                  + pb_z[k] * dd_5[k];

        t_12[k] = pa_x[k] * pf_12[k];

        t_13[k] = pa_x[k] * pf_14[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pb_x, pb_y, pb_z, pd_1, dp0_3, dp0_4, \
                         dp1_3, dp1_4, dd_7, dd_8, dd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * dp0_3[k]
                  - f_2 * dp1_3[k]
                  + pb_x[k] * dd_7[k];

        t_15[k] = pb_z[k] * dd_7[k];

        t_16[k] = pb_x[k] * dd_8[k];

        t_17[k] = pb_x[k] * dd_9[k];

        t_18[k] = f_0 * pd_1[k]
                  + f_1 * dp0_4[k]
                  - f_2 * dp1_4[k]
                  + pb_y[k] * dd_8[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pa_z, pb_y, pb_z, pd_2, pf_5, \
                         pf_9, dp0_5, dp1_5, dd_8, dd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_z[k] * dd_8[k];

        t_20[k] = f_0 * pd_2[k]
                  + pb_y[k] * dd_9[k];

        t_21[k] = f_1 * dp0_5[k]
                  - f_2 * dp1_5[k]
                  + pb_z[k] * dd_9[k];

        t_22[k] = pa_y[k] * pf_9[k];

        t_23[k] = pa_z[k] * pf_5[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pb_x, pb_y, pb_z, pd_1, pd_5, pf_14, \
                         dp0_6, dp1_6, dd_10, dd_11, dd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_3 * pd_1[k]
                  + pb_z[k] * dd_10[k];

        t_25[k] = f_3 * pd_5[k]
                  + pb_y[k] * dd_11[k];

        t_26[k] = pa_y[k] * pf_14[k];

        t_27[k] = f_1 * dp0_6[k]
                  - f_2 * dp1_6[k]
                  + pb_x[k] * dd_12[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, t_33, pb_x, pb_y, pb_z, pd_3, pd_4, \
                         dp0_7, dp1_7, dd_12, dd_13, dd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * dd_12[k];

        t_29[k] = f_0 * pd_3[k]
                  + pb_z[k] * dd_12[k];

        t_30[k] = pb_x[k] * dd_13[k];

        t_31[k] = pb_x[k] * dd_14[k];

        t_32[k] = f_1 * dp0_7[k]
                  - f_2 * dp1_7[k]
                  + pb_y[k] * dd_13[k];

        t_33[k] = f_0 * pd_4[k]
                  + pb_z[k] * dd_13[k];
    }

#pragma omp simd aligned(t_34, t_35, pb_y, pb_z, pd_5, dp0_8, dp1_8, \
                         dd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * dd_14[k];

        t_35[k] = f_0 * pd_5[k]
                  + f_1 * dp0_8[k]
                  - f_2 * dp1_8[k]
                  + pb_z[k] * dd_14[k];
    }
}

auto
compute_prim_df_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pd, const size_t pf,
                                     const size_t dp0, const size_t dp1, const size_t dd,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pd_0, dp0_0, dp0_1, dp1_0, \
                         dp1_1, dd_0, dd_1, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_y[k] * dd_0[k];

        t_2[k] = pb_z[k] * dd_0[k];

        t_3[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_4[k] = pb_y[k] * dd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pa_z, pb_z, pd_0, pf_0, pf_1, \
                         dp0_2, dp1_2, dd_2, dd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];

        t_6[k] = pa_y[k] * pf_0[k];

        t_7[k] = pa_x[k] * pf_1[k];

        t_8[k] = pa_z[k] * pf_0[k];

        t_9[k] = f_3 * pd_0[k]
                 + pb_z[k] * dd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pf_2, dp0_3, dp1_3, dd_4, dd_5, \
                         dd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * pf_2[k];

        t_11[k] = f_1 * dp0_3[k]
                  - f_2 * dp1_3[k]
                  + pb_x[k] * dd_4[k];

        t_12[k] = pb_x[k] * dd_5[k];

        t_13[k] = pb_x[k] * dd_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_y, pb_z, pd_1, pd_2, dp0_4, dp0_5, dp1_4, \
                         dp1_5, dd_5, dd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * pd_1[k]
                  + f_1 * dp0_4[k]
                  - f_2 * dp1_4[k]
                  + pb_y[k] * dd_5[k];

        t_15[k] = pb_z[k] * dd_5[k];

        t_16[k] = f_0 * pd_2[k]
                  + pb_y[k] * dd_6[k];

        t_17[k] = f_1 * dp0_5[k]
                  - f_2 * dp1_5[k]
                  + pb_z[k] * dd_6[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_y, pa_z, pb_y, pb_z, pd_1, pd_5, pf_1, \
                         pf_2, dd_7, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pa_z[k] * pf_1[k];

        t_19[k] = f_3 * pd_1[k]
                  + pb_z[k] * dd_7[k];

        t_20[k] = f_3 * pd_5[k]
                  + pb_y[k] * dd_8[k];

        t_21[k] = pa_y[k] * pf_2[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pb_x, pb_y, pb_z, pd_3, dp0_6, dp0_7, \
                         dp1_6, dp1_7, dd_9, dd_10, dd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_1 * dp0_6[k]
                  - f_2 * dp1_6[k]
                  + pb_x[k] * dd_9[k];

        t_23[k] = f_0 * pd_3[k]
                  + pb_z[k] * dd_9[k];

        t_24[k] = pb_x[k] * dd_10[k];

        t_25[k] = pb_x[k] * dd_11[k];

        t_26[k] = f_1 * dp0_7[k]
                  - f_2 * dp1_7[k]
                  + pb_y[k] * dd_10[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pb_z, pd_4, pd_5, dp0_8, dp1_8, dd_10, \
                         dd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * pd_4[k]
                  + pb_z[k] * dd_10[k];

        t_28[k] = pb_y[k] * dd_11[k];

        t_29[k] = f_0 * pd_5[k]
                  + f_1 * dp0_8[k]
                  - f_2 * dp1_8[k]
                  + pb_z[k] * dd_11[k];
    }
}

auto
compute_prim_df_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pd, const size_t pf,
                                     const size_t dp0, const size_t dp1, const size_t dd,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

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

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_8 = buffer.data(pf + 8);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pd_0, dp0_0, dp0_1, dp1_0, \
                         dp1_1, dd_0, dd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_y[k] * dd_0[k];

        t_2[k] = pb_z[k] * dd_0[k];

        t_3[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_4[k] = pb_z[k] * dd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pa_z, pb_y, pb_z, pf_0, pf_2, \
                         dp0_2, dp1_2, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_y[k] * dd_2[k];

        t_6[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];

        t_7[k] = pa_y[k] * pf_0[k];

        t_8[k] = pa_x[k] * pf_2[k];

        t_9[k] = pa_z[k] * pf_0[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_x, pb_x, pb_z, pf_8, dp0_3, dp1_3, \
                         dd_3, dd_4, dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * pf_8[k];

        t_11[k] = f_1 * dp0_3[k]
                  - f_2 * dp1_3[k]
                  + pb_x[k] * dd_3[k];

        t_12[k] = pb_z[k] * dd_3[k];

        t_13[k] = pb_x[k] * dd_4[k];

        t_14[k] = pb_x[k] * dd_5[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_z, pb_y, pb_z, pd_1, pf_2, dp0_4, dp0_5, \
                         dp1_4, dp1_5, dd_4, dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * pd_1[k]
                  + f_1 * dp0_4[k]
                  - f_2 * dp1_4[k]
                  + pb_y[k] * dd_4[k];

        t_16[k] = pb_z[k] * dd_4[k];

        t_17[k] = f_1 * dp0_5[k]
                  - f_2 * dp1_5[k]
                  + pb_z[k] * dd_5[k];

        t_18[k] = pa_z[k] * pf_2[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pb_x, pb_y, pf_8, dp0_6, dp1_6, \
                         dd_6, dd_7, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * pf_8[k];

        t_20[k] = f_1 * dp0_6[k]
                  - f_2 * dp1_6[k]
                  + pb_x[k] * dd_6[k];

        t_21[k] = pb_y[k] * dd_6[k];

        t_22[k] = pb_x[k] * dd_7[k];

        t_23[k] = pb_x[k] * dd_8[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_y, pb_z, pd_2, dp0_7, dp0_8, dp1_7, dp1_8, dd_7, \
                         dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_1 * dp0_7[k]
                  - f_2 * dp1_7[k]
                  + pb_y[k] * dd_7[k];

        t_25[k] = pb_y[k] * dd_8[k];

        t_26[k] = f_0 * pd_2[k]
                  + f_1 * dp0_8[k]
                  - f_2 * dp1_8[k]
                  + pb_z[k] * dd_8[k];
    }
}

auto
compute_prim_df_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t pd, const size_t pf,
                                     const size_t dp0, const size_t dp1, const size_t dd,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pd_0, dp0_0, dp0_1, dp1_0, \
                         dp1_1, dd_0, dd_1, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_y[k] * dd_0[k];

        t_2[k] = pb_z[k] * dd_0[k];

        t_3[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_4[k] = pb_y[k] * dd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_z, pb_x, pb_z, pf_0, dp0_2, dp0_3, dp1_2, \
                         dp1_3, dd_2, dd_3, dd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];

        t_6[k] = pa_z[k] * pf_0[k];

        t_7[k] = f_1 * dp0_3[k]
                 - f_2 * dp1_3[k]
                 + pb_x[k] * dd_3[k];

        t_8[k] = pb_x[k] * dd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_x, pb_y, pb_z, pd_1, dp0_4, dp0_5, dp1_4, \
                         dp1_5, dd_4, dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pb_x[k] * dd_5[k];

        t_10[k] = f_0 * pd_1[k]
                  + f_1 * dp0_4[k]
                  - f_2 * dp1_4[k]
                  + pb_y[k] * dd_4[k];

        t_11[k] = pb_z[k] * dd_4[k];

        t_12[k] = f_1 * dp0_5[k]
                  - f_2 * dp1_5[k]
                  + pb_z[k] * dd_5[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_y, pa_z, pb_x, pf_1, pf_2, dp0_6, \
                         dp1_6, dd_6, dd_7, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_z[k] * pf_1[k];

        t_14[k] = pa_y[k] * pf_2[k];

        t_15[k] = f_1 * dp0_6[k]
                  - f_2 * dp1_6[k]
                  + pb_x[k] * dd_6[k];

        t_16[k] = pb_x[k] * dd_7[k];

        t_17[k] = pb_x[k] * dd_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_y, pb_z, pd_2, dp0_7, dp0_8, dp1_7, dp1_8, dd_7, \
                         dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_1 * dp0_7[k]
                  - f_2 * dp1_7[k]
                  + pb_y[k] * dd_7[k];

        t_19[k] = pb_y[k] * dd_8[k];

        t_20[k] = f_0 * pd_2[k]
                  + f_1 * dp0_8[k]
                  - f_2 * dp1_8[k]
                  + pb_z[k] * dd_8[k];
    }
}

auto
compute_prim_df_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t pd, const size_t dp0, const size_t dp1,
                                     const size_t dd, const size_t ncols, const double alpha,
                                     const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pd_0, pd_1, dp0_0, dp0_1, dp1_0, dp1_1, dd_0, \
                         dd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = f_0 * pd_1[k]
                 + f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pd_2, dp0_2, dp1_2, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pd_2[k]
                 + f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];
    }
}

auto
compute_prim_df_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pd, const size_t dp0, const size_t dp1,
                                      const size_t dd, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_9 = buffer.data(dp1 + 9);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_15 = buffer.data(dd + 15);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pd_0, pd_1, dp0_0, dp0_1, dp1_0, dp1_4, dd_0, \
                         dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = f_0 * pd_1[k]
                 + f_1 * dp0_1[k]
                 - f_2 * dp1_4[k]
                 + pb_y[k] * dd_8[k];
    }

#pragma omp simd aligned(t_2, pb_z, pd_2, dp0_2, dp1_9, dd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pd_2[k]
                 + f_1 * dp0_2[k]
                 - f_2 * dp1_9[k]
                 + pb_z[k] * dd_15[k];
    }
}

auto
compute_prim_df_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pd, const size_t dp0, const size_t dp1,
                                      const size_t dd, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);
    const auto *dp0_9 = buffer.data(dp0 + 9);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);
    const auto *dp1_9 = buffer.data(dp1 + 9);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pd_0, dp0_0, dp0_1, dp0_2, dp1_0, \
                         dp1_1, dp1_2, dd_0, dd_1, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_2[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, pb_z, pd_1, dp0_3, dp0_4, dp0_5, dp1_3, \
                         dp1_4, dp1_5, dd_3, dd_4, dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_1 * dp0_3[k]
                 - f_2 * dp1_3[k]
                 + pb_x[k] * dd_3[k];

        t_4[k] = f_0 * pd_1[k]
                 + f_1 * dp0_4[k]
                 - f_2 * dp1_4[k]
                 + pb_y[k] * dd_4[k];

        t_5[k] = f_1 * dp0_5[k]
                 - f_2 * dp1_5[k]
                 + pb_z[k] * dd_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, pb_z, pd_2, dp0_7, dp0_8, dp0_9, dp1_7, \
                         dp1_8, dp1_9, dd_6, dd_7, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * dp0_7[k]
                 - f_2 * dp1_7[k]
                 + pb_x[k] * dd_6[k];

        t_7[k] = f_1 * dp0_8[k]
                 - f_2 * dp1_8[k]
                 + pb_y[k] * dd_7[k];

        t_8[k] = f_0 * pd_2[k]
                 + f_1 * dp0_9[k]
                 - f_2 * dp1_9[k]
                 + pb_z[k] * dd_8[k];
    }
}

auto
compute_prim_df_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pd, const size_t dp0, const size_t dp1,
                                      const size_t dd, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_11 = buffer.data(dd + 11);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pd_0, pd_1, dp0_0, dp0_1, dp1_0, dp1_4, dd_0, \
                         dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = f_0 * pd_1[k]
                 + f_1 * dp0_1[k]
                 - f_2 * dp1_4[k]
                 + pb_y[k] * dd_5[k];
    }

#pragma omp simd aligned(t_2, pb_z, pd_2, dp0_2, dp1_8, dd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pd_2[k]
                 + f_1 * dp0_2[k]
                 - f_2 * dp1_8[k]
                 + pb_z[k] * dd_11[k];
    }
}

auto
compute_prim_df_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t pd, const size_t pf,
                                      const size_t dp0, const size_t dp1, const size_t dd,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pd_0, dp0_0, dp0_1, dp1_0, \
                         dp1_1, dd_0, dd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_y[k] * dd_0[k];

        t_2[k] = pb_z[k] * dd_0[k];

        t_3[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_4[k] = pb_z[k] * dd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pa_y, pa_z, pb_z, pd_0, pf_0, pf_1, \
                         dp0_2, dp1_2, dd_2, dd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];

        t_6[k] = pa_y[k] * pf_0[k];

        t_7[k] = pa_x[k] * pf_1[k];

        t_8[k] = pa_z[k] * pf_0[k];

        t_9[k] = f_3 * pd_0[k]
                 + pb_z[k] * dd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, pf_2, dp0_3, dp1_3, dd_4, \
                         dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * pf_2[k];

        t_11[k] = f_1 * dp0_3[k]
                  - f_2 * dp1_3[k]
                  + pb_x[k] * dd_4[k];

        t_12[k] = pb_z[k] * dd_4[k];

        t_13[k] = pb_x[k] * dd_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_y, pb_z, pd_1, pd_2, dp0_4, dp0_5, dp1_4, \
                         dp1_5, dd_5, dd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * pd_1[k]
                  + f_1 * dp0_4[k]
                  - f_2 * dp1_4[k]
                  + pb_y[k] * dd_5[k];

        t_15[k] = pb_z[k] * dd_5[k];

        t_16[k] = f_0 * pd_2[k]
                  + pb_y[k] * dd_6[k];

        t_17[k] = f_1 * dp0_5[k]
                  - f_2 * dp1_5[k]
                  + pb_z[k] * dd_6[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_y, pa_z, pb_y, pb_z, pd_1, pd_5, pf_1, \
                         pf_2, dd_7, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pa_z[k] * pf_1[k];

        t_19[k] = f_3 * pd_1[k]
                  + pb_z[k] * dd_7[k];

        t_20[k] = f_3 * pd_5[k]
                  + pb_y[k] * dd_8[k];

        t_21[k] = pa_y[k] * pf_2[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pb_x, pb_y, pb_z, pd_3, dp0_6, dp0_7, \
                         dp1_6, dp1_7, dd_9, dd_10, dd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_1 * dp0_6[k]
                  - f_2 * dp1_6[k]
                  + pb_x[k] * dd_9[k];

        t_23[k] = pb_y[k] * dd_9[k];

        t_24[k] = f_0 * pd_3[k]
                  + pb_z[k] * dd_9[k];

        t_25[k] = pb_x[k] * dd_11[k];

        t_26[k] = f_1 * dp0_7[k]
                  - f_2 * dp1_7[k]
                  + pb_y[k] * dd_10[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pb_z, pd_4, pd_5, dp0_8, dp1_8, dd_10, \
                         dd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * pd_4[k]
                  + pb_z[k] * dd_10[k];

        t_28[k] = pb_y[k] * dd_11[k];

        t_29[k] = f_0 * pd_5[k]
                  + f_1 * dp0_8[k]
                  - f_2 * dp1_8[k]
                  + pb_z[k] * dd_11[k];
    }
}

auto
compute_prim_df_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pd, const size_t dp0, const size_t dp1,
                                      const size_t dd, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_8 = buffer.data(dd + 8);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pd_0, pd_1, dp0_0, dp0_1, dp1_0, dp1_4, dd_0, \
                         dd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = f_0 * pd_1[k]
                 + f_1 * dp0_1[k]
                 - f_2 * dp1_4[k]
                 + pb_y[k] * dd_4[k];
    }

#pragma omp simd aligned(t_2, pb_z, pd_2, dp0_2, dp1_8, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pd_2[k]
                 + f_1 * dp0_2[k]
                 - f_2 * dp1_8[k]
                 + pb_z[k] * dd_8[k];
    }
}

auto
compute_prim_df_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t pd, const size_t pf,
                                      const size_t dp0, const size_t dp1, const size_t dd,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_8 = buffer.data(pf + 8);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pd_0, dp0_0, dp0_1, dp1_0, \
                         dp1_1, dd_0, dd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_y[k] * dd_0[k];

        t_2[k] = pb_z[k] * dd_0[k];

        t_3[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_4[k] = pb_z[k] * dd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_z, pb_x, pb_y, pb_z, pf_0, dp0_2, dp0_3, \
                         dp1_2, dp1_3, dd_2, dd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_y[k] * dd_2[k];

        t_6[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];

        t_7[k] = pa_z[k] * pf_0[k];

        t_8[k] = f_1 * dp0_3[k]
                 - f_2 * dp1_3[k]
                 + pb_x[k] * dd_3[k];

        t_9[k] = pb_z[k] * dd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pb_x, pb_y, pb_z, pd_1, dp0_4, dp0_5, \
                         dp1_4, dp1_5, dd_4, dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_x[k] * dd_4[k];

        t_11[k] = pb_x[k] * dd_5[k];

        t_12[k] = f_0 * pd_1[k]
                  + f_1 * dp0_4[k]
                  - f_2 * dp1_4[k]
                  + pb_y[k] * dd_4[k];

        t_13[k] = pb_z[k] * dd_4[k];

        t_14[k] = f_1 * dp0_5[k]
                  - f_2 * dp1_5[k]
                  + pb_z[k] * dd_5[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pa_z, pb_x, pb_y, pf_2, pf_8, \
                         dp0_6, dp1_6, dd_6, dd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_z[k] * pf_2[k];

        t_16[k] = pa_y[k] * pf_8[k];

        t_17[k] = f_1 * dp0_6[k]
                  - f_2 * dp1_6[k]
                  + pb_x[k] * dd_6[k];

        t_18[k] = pb_y[k] * dd_6[k];

        t_19[k] = pb_x[k] * dd_7[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pb_x, pb_y, pb_z, pd_2, dp0_7, dp0_8, dp1_7, \
                         dp1_8, dd_7, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_x[k] * dd_8[k];

        t_21[k] = f_1 * dp0_7[k]
                  - f_2 * dp1_7[k]
                  + pb_y[k] * dd_7[k];

        t_22[k] = pb_y[k] * dd_8[k];

        t_23[k] = f_0 * pd_2[k]
                  + f_1 * dp0_8[k]
                  - f_2 * dp1_8[k]
                  + pb_z[k] * dd_8[k];
    }
}

auto
compute_prim_df_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t pd, const size_t pf,
                                      const size_t dp0, const size_t dp1, const size_t dd,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pd_0, dp0_0, dp0_1, dp1_0, \
                         dp1_1, dd_0, dd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_y[k] * dd_0[k];

        t_2[k] = pb_z[k] * dd_0[k];

        t_3[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_4[k] = pb_z[k] * dd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_z, pb_x, pb_y, pb_z, pf_0, dp0_2, dp0_3, \
                         dp1_2, dp1_3, dd_2, dd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_y[k] * dd_2[k];

        t_6[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];

        t_7[k] = pa_z[k] * pf_0[k];

        t_8[k] = f_1 * dp0_3[k]
                 - f_2 * dp1_3[k]
                 + pb_x[k] * dd_3[k];

        t_9[k] = pb_z[k] * dd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pb_x, pb_y, pb_z, pd_1, dp0_4, dp0_5, \
                         dp1_4, dp1_5, dd_4, dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_x[k] * dd_4[k];

        t_11[k] = pb_x[k] * dd_5[k];

        t_12[k] = f_0 * pd_1[k]
                  + f_1 * dp0_4[k]
                  - f_2 * dp1_4[k]
                  + pb_y[k] * dd_4[k];

        t_13[k] = pb_z[k] * dd_4[k];

        t_14[k] = f_1 * dp0_5[k]
                  - f_2 * dp1_5[k]
                  + pb_z[k] * dd_5[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pa_z, pb_x, pb_y, pf_1, pf_2, \
                         dp0_6, dp1_6, dd_6, dd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_z[k] * pf_1[k];

        t_16[k] = pa_y[k] * pf_2[k];

        t_17[k] = f_1 * dp0_6[k]
                  - f_2 * dp1_6[k]
                  + pb_x[k] * dd_6[k];

        t_18[k] = pb_y[k] * dd_6[k];

        t_19[k] = pb_x[k] * dd_7[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pb_x, pb_y, pb_z, pd_2, dp0_7, dp0_8, dp1_7, \
                         dp1_8, dd_7, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_x[k] * dd_8[k];

        t_21[k] = f_1 * dp0_7[k]
                  - f_2 * dp1_7[k]
                  + pb_y[k] * dd_7[k];

        t_22[k] = pb_y[k] * dd_8[k];

        t_23[k] = f_0 * pd_2[k]
                  + f_1 * dp0_8[k]
                  - f_2 * dp1_8[k]
                  + pb_z[k] * dd_8[k];
    }
}

auto
compute_prim_df_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pd, const size_t dp0, const size_t dp1,
                                      const size_t dd, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pd_0, dp0_0, dp0_1, dp1_0, \
                         dp1_1, dd_0, dd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_y[k] * dd_0[k];

        t_2[k] = pb_z[k] * dd_0[k];

        t_3[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_4[k] = pb_z[k] * dd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, dp0_2, dp0_3, dp1_2, \
                         dp1_3, dd_2, dd_3, dd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_y[k] * dd_2[k];

        t_6[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];

        t_7[k] = f_1 * dp0_3[k]
                 - f_2 * dp1_3[k]
                 + pb_x[k] * dd_3[k];

        t_8[k] = pb_z[k] * dd_3[k];

        t_9[k] = pb_x[k] * dd_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_x, pb_y, pb_z, pd_1, dp0_4, dp0_5, dp1_4, \
                         dp1_5, dd_4, dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_x[k] * dd_5[k];

        t_11[k] = f_0 * pd_1[k]
                  + f_1 * dp0_4[k]
                  - f_2 * dp1_4[k]
                  + pb_y[k] * dd_4[k];

        t_12[k] = pb_z[k] * dd_4[k];

        t_13[k] = f_1 * dp0_5[k]
                  - f_2 * dp1_5[k]
                  + pb_z[k] * dd_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, t_19, pb_x, pb_y, dp0_6, dp0_7, dp1_6, \
                         dp1_7, dd_6, dd_7, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * dp0_6[k]
                  - f_2 * dp1_6[k]
                  + pb_x[k] * dd_6[k];

        t_15[k] = pb_y[k] * dd_6[k];

        t_16[k] = pb_x[k] * dd_7[k];

        t_17[k] = pb_x[k] * dd_8[k];

        t_18[k] = f_1 * dp0_7[k]
                  - f_2 * dp1_7[k]
                  + pb_y[k] * dd_7[k];

        t_19[k] = pb_y[k] * dd_8[k];
    }

#pragma omp simd aligned(t_20, pb_z, pd_2, dp0_8, dp1_8, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * pd_2[k]
                  + f_1 * dp0_8[k]
                  - f_2 * dp1_8[k]
                  + pb_z[k] * dd_8[k];
    }
}

auto
compute_prim_df_electron_repulsion_18(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pd, const size_t dp0, const size_t dp1,
                                      const size_t dd, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_9 = buffer.data(dp1 + 9);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_8 = buffer.data(dd + 8);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pd_0, pd_1, dp0_0, dp0_1, dp1_0, dp1_4, dd_0, \
                         dd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = f_0 * pd_1[k]
                 + f_1 * dp0_1[k]
                 - f_2 * dp1_4[k]
                 + pb_y[k] * dd_4[k];
    }

#pragma omp simd aligned(t_2, pb_z, pd_2, dp0_2, dp1_9, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pd_2[k]
                 + f_1 * dp0_2[k]
                 - f_2 * dp1_9[k]
                 + pb_z[k] * dd_8[k];
    }
}

auto
compute_prim_df_electron_repulsion_19(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pd, const size_t dp0, const size_t dp1,
                                      const size_t dd, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_9 = buffer.data(dp0 + 9);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pd_0, pd_1, dp0_0, dp0_4, dp1_0, dp1_4, dd_0, \
                         dd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = f_0 * pd_1[k]
                 + f_1 * dp0_4[k]
                 - f_2 * dp1_4[k]
                 + pb_y[k] * dd_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pd_2, dp0_9, dp1_8, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pd_2[k]
                 + f_1 * dp0_9[k]
                 - f_2 * dp1_8[k]
                 + pb_z[k] * dd_2[k];
    }
}

auto
compute_prim_df_electron_repulsion_20(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pd, const size_t dp0, const size_t dp1,
                                      const size_t dd, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pd_0, dp0_0, dp0_1, dp0_2, dp1_0, \
                         dp1_1, dp1_2, dd_0, dd_1, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_2[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, pb_z, pd_1, dp0_3, dp0_4, dp0_5, dp1_3, \
                         dp1_4, dp1_5, dd_3, dd_4, dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_1 * dp0_3[k]
                 - f_2 * dp1_3[k]
                 + pb_x[k] * dd_3[k];

        t_4[k] = f_0 * pd_1[k]
                 + f_1 * dp0_4[k]
                 - f_2 * dp1_4[k]
                 + pb_y[k] * dd_4[k];

        t_5[k] = f_1 * dp0_5[k]
                 - f_2 * dp1_5[k]
                 + pb_z[k] * dd_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, pb_z, pd_2, dp0_6, dp0_7, dp0_8, dp1_6, \
                         dp1_7, dp1_8, dd_6, dd_7, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * dp0_6[k]
                 - f_2 * dp1_6[k]
                 + pb_x[k] * dd_6[k];

        t_7[k] = f_1 * dp0_7[k]
                 - f_2 * dp1_7[k]
                 + pb_y[k] * dd_7[k];

        t_8[k] = f_0 * pd_2[k]
                 + f_1 * dp0_8[k]
                 - f_2 * dp1_8[k]
                 + pb_z[k] * dd_8[k];
    }
}

auto
compute_prim_df_electron_repulsion_21(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pd, const size_t dp0, const size_t dp1,
                                      const size_t dd, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pd_0, pd_1, dp0_0, dp0_4, dp1_0, dp1_4, dd_0, \
                         dd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = f_0 * pd_1[k]
                 + f_1 * dp0_4[k]
                 - f_2 * dp1_4[k]
                 + pb_y[k] * dd_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pd_2, dp0_8, dp1_8, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pd_2[k]
                 + f_1 * dp0_8[k]
                 - f_2 * dp1_8[k]
                 + pb_z[k] * dd_2[k];
    }
}

auto
compute_prim_df_electron_repulsion_22(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t pd, const size_t pf,
                                      const size_t dp0, const size_t dp1, const size_t dd,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pd_0, dp0_0, dp0_1, dp1_0, \
                         dp1_1, dd_0, dd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_y[k] * dd_0[k];

        t_2[k] = pb_z[k] * dd_0[k];

        t_3[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_4[k] = pb_z[k] * dd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_z, pb_x, pb_z, pf_0, dp0_2, dp0_3, dp1_2, \
                         dp1_3, dd_2, dd_3, dd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];

        t_6[k] = pa_z[k] * pf_0[k];

        t_7[k] = f_1 * dp0_3[k]
                 - f_2 * dp1_3[k]
                 + pb_x[k] * dd_3[k];

        t_8[k] = pb_z[k] * dd_3[k];

        t_9[k] = pb_x[k] * dd_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_z, pb_y, pb_z, pd_1, pf_1, dp0_4, dp0_5, \
                         dp1_4, dp1_5, dd_4, dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * pd_1[k]
                  + f_1 * dp0_4[k]
                  - f_2 * dp1_4[k]
                  + pb_y[k] * dd_4[k];

        t_11[k] = pb_z[k] * dd_4[k];

        t_12[k] = f_1 * dp0_5[k]
                  - f_2 * dp1_5[k]
                  + pb_z[k] * dd_5[k];

        t_13[k] = pa_z[k] * pf_1[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pa_y, pb_x, pb_y, pf_2, dp0_6, dp0_7, \
                         dp1_6, dp1_7, dd_6, dd_7, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_y[k] * pf_2[k];

        t_15[k] = f_1 * dp0_6[k]
                  - f_2 * dp1_6[k]
                  + pb_x[k] * dd_6[k];

        t_16[k] = pb_y[k] * dd_6[k];

        t_17[k] = pb_x[k] * dd_8[k];

        t_18[k] = f_1 * dp0_7[k]
                  - f_2 * dp1_7[k]
                  + pb_y[k] * dd_7[k];
    }

#pragma omp simd aligned(t_19, t_20, pb_y, pb_z, pd_2, dp0_8, dp1_8, \
                         dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * dd_8[k];

        t_20[k] = f_0 * pd_2[k]
                  + f_1 * dp0_8[k]
                  - f_2 * dp1_8[k]
                  + pb_z[k] * dd_8[k];
    }
}

auto
compute_prim_df_electron_repulsion_23(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pd, const size_t dp0, const size_t dp1,
                                      const size_t dd, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);

#pragma omp simd aligned(t_0, t_1, pb_x, pb_y, pd_0, pd_1, dp0_0, dp0_4, dp1_0, dp1_1, dd_0, \
                         dd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = f_0 * pd_1[k]
                 + f_1 * dp0_4[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];
    }

#pragma omp simd aligned(t_2, pb_z, pd_2, dp0_8, dp1_2, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2[k] = f_0 * pd_2[k]
                 + f_1 * dp0_8[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];
    }
}

auto
compute_prim_df_electron_repulsion_24(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pd, const size_t dp0, const size_t dp1,
                                      const size_t dd, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pd_0, pd_1, dp0_0, dp0_1, dp1_0, \
                         dp1_1, dd_0, dd_1, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_x[k] * dd_1[k];

        t_2[k] = f_0 * pd_1[k]
                 + f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_3[k] = pb_x[k] * dd_2[k];

        t_4[k] = pb_y[k] * dd_2[k];
    }

#pragma omp simd aligned(t_5, pb_z, pd_2, dp0_2, dp1_2, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * pd_2[k]
                 + f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];
    }
}

auto
compute_prim_df_electron_repulsion_25(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pd, const size_t dp0, const size_t dp1,
                                      const size_t dd, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_8 = buffer.data(dd + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pd_0, pd_1, dp0_0, dp0_1, dp1_0, \
                         dp1_4, dd_0, dd_4, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_x[k] * dd_4[k];

        t_2[k] = f_0 * pd_1[k]
                 + f_1 * dp0_1[k]
                 - f_2 * dp1_4[k]
                 + pb_y[k] * dd_4[k];

        t_3[k] = pb_x[k] * dd_8[k];

        t_4[k] = pb_y[k] * dd_8[k];
    }

#pragma omp simd aligned(t_5, pb_z, pd_2, dp0_2, dp1_8, dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * pd_2[k]
                 + f_1 * dp0_2[k]
                 - f_2 * dp1_8[k]
                 + pb_z[k] * dd_8[k];
    }
}

auto
compute_prim_df_electron_repulsion_26(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pd, const size_t dp0, const size_t dp1,
                                      const size_t dd, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_3 = buffer.data(dp0 + 3);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_5 = buffer.data(dp0 + 5);
    const auto *dp0_6 = buffer.data(dp0 + 6);
    const auto *dp0_7 = buffer.data(dp0 + 7);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_3 = buffer.data(dp1 + 3);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_5 = buffer.data(dp1 + 5);
    const auto *dp1_6 = buffer.data(dp1 + 6);
    const auto *dp1_7 = buffer.data(dp1 + 7);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pd_0, dp0_0, dp0_1, dp0_2, dp1_0, \
                         dp1_1, dp1_2, dd_0, dd_1, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_1[k];

        t_2[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, pd_1, dp0_3, dp0_4, dp1_3, dp1_4, dd_3, \
                         dd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_1 * dp0_3[k]
                 - f_2 * dp1_3[k]
                 + pb_x[k] * dd_3[k];

        t_4[k] = pb_x[k] * dd_4[k];

        t_5[k] = f_0 * pd_1[k]
                 + f_1 * dp0_4[k]
                 - f_2 * dp1_4[k]
                 + pb_y[k] * dd_4[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_z, dp0_5, dp0_6, dp1_5, dp1_6, dd_5, dd_6, \
                         dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * dp0_5[k]
                 - f_2 * dp1_5[k]
                 + pb_z[k] * dd_5[k];

        t_7[k] = f_1 * dp0_6[k]
                 - f_2 * dp1_6[k]
                 + pb_x[k] * dd_6[k];

        t_8[k] = pb_x[k] * dd_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_y, pb_z, pd_2, dp0_7, dp0_8, dp1_7, dp1_8, dd_7, \
                         dd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_1 * dp0_7[k]
                 - f_2 * dp1_7[k]
                 + pb_y[k] * dd_7[k];

        t_10[k] = pb_y[k] * dd_8[k];

        t_11[k] = f_0 * pd_2[k]
                  + f_1 * dp0_8[k]
                  - f_2 * dp1_8[k]
                  + pb_z[k] * dd_8[k];
    }
}

auto
compute_prim_df_electron_repulsion_27(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t pd, const size_t dp0, const size_t dp1,
                                      const size_t dd, const size_t ncols, const double alpha,
                                      const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_4 = buffer.data(dp0 + 4);
    const auto *dp0_8 = buffer.data(dp0 + 8);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_4 = buffer.data(dp1 + 4);
    const auto *dp1_8 = buffer.data(dp1 + 8);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pd_0, pd_1, dp0_0, dp0_4, dp1_0, \
                         dp1_4, dd_0, dd_1, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_x[k] * dd_1[k];

        t_2[k] = f_0 * pd_1[k]
                 + f_1 * dp0_4[k]
                 - f_2 * dp1_4[k]
                 + pb_y[k] * dd_1[k];

        t_3[k] = pb_x[k] * dd_2[k];

        t_4[k] = pb_y[k] * dd_2[k];
    }

#pragma omp simd aligned(t_5, pb_z, pd_2, dp0_8, dp1_8, dd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * pd_2[k]
                 + f_1 * dp0_8[k]
                 - f_2 * dp1_8[k]
                 + pb_z[k] * dd_2[k];
    }
}

}  // namespace simdt2ceri
