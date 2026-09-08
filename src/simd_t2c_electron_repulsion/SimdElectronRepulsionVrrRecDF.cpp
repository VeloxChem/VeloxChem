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
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_9 = buffer.data(pd + 9);
    const auto *pd_11 = buffer.data(pd + 11);
    const auto *pd_12 = buffer.data(pd + 12);
    const auto *pd_15 = buffer.data(pd + 15);
    const auto *pd_17 = buffer.data(pd + 17);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_16 = buffer.data(pf + 16);
    const auto *pf_18 = buffer.data(pf + 18);
    const auto *pf_19 = buffer.data(pf + 19);
    const auto *pf_20 = buffer.data(pf + 20);
    const auto *pf_22 = buffer.data(pf + 22);
    const auto *pf_26 = buffer.data(pf + 26);
    const auto *pf_27 = buffer.data(pf + 27);
    const auto *pf_29 = buffer.data(pf + 29);

    const auto *dp0_0 = buffer.data(dp0 + 0);
    const auto *dp0_1 = buffer.data(dp0 + 1);
    const auto *dp0_2 = buffer.data(dp0 + 2);
    const auto *dp0_9 = buffer.data(dp0 + 9);
    const auto *dp0_10 = buffer.data(dp0 + 10);
    const auto *dp0_11 = buffer.data(dp0 + 11);
    const auto *dp0_15 = buffer.data(dp0 + 15);
    const auto *dp0_16 = buffer.data(dp0 + 16);
    const auto *dp0_17 = buffer.data(dp0 + 17);

    const auto *dp1_0 = buffer.data(dp1 + 0);
    const auto *dp1_1 = buffer.data(dp1 + 1);
    const auto *dp1_2 = buffer.data(dp1 + 2);
    const auto *dp1_9 = buffer.data(dp1 + 9);
    const auto *dp1_10 = buffer.data(dp1 + 10);
    const auto *dp1_11 = buffer.data(dp1 + 11);
    const auto *dp1_15 = buffer.data(dp1 + 15);
    const auto *dp1_16 = buffer.data(dp1 + 16);
    const auto *dp1_17 = buffer.data(dp1 + 17);

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_17 = buffer.data(dd + 17);
    const auto *dd_18 = buffer.data(dd + 18);
    const auto *dd_21 = buffer.data(dd + 21);
    const auto *dd_22 = buffer.data(dd + 22);
    const auto *dd_23 = buffer.data(dd + 23);
    const auto *dd_27 = buffer.data(dd + 27);
    const auto *dd_28 = buffer.data(dd + 28);
    const auto *dd_29 = buffer.data(dd + 29);
    const auto *dd_30 = buffer.data(dd + 30);
    const auto *dd_33 = buffer.data(dd + 33);
    const auto *dd_34 = buffer.data(dd + 34);
    const auto *dd_35 = buffer.data(dd + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, pd_0, pd_3, dp0_0, dp1_0, \
                         dd_0, dd_2, dd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pd_0[k]
                 + f_1 * dp0_0[k]
                 - f_2 * dp1_0[k]
                 + pb_x[k] * dd_0[k];

        t_1[k] = pb_y[k] * dd_0[k];

        t_2[k] = pb_z[k] * dd_0[k];

        t_3[k] = f_0 * pd_3[k]
                 + pb_x[k] * dd_3[k];

        t_4[k] = pb_y[k] * dd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, pd_5, dp0_1, dp0_2, dp1_1, \
                         dp1_2, dd_3, dd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * pd_5[k]
                 + pb_x[k] * dd_5[k];

        t_6[k] = f_1 * dp0_1[k]
                 - f_2 * dp1_1[k]
                 + pb_y[k] * dd_3[k];

        t_7[k] = pb_z[k] * dd_3[k];

        t_8[k] = pb_y[k] * dd_5[k];

        t_9[k] = f_1 * dp0_2[k]
                 - f_2 * dp1_2[k]
                 + pb_z[k] * dd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, pb_y, pb_z, pd_0, pd_9, \
                         pf_0, dd_6, dd_7, dd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * pf_0[k];

        t_11[k] = f_3 * pd_0[k]
                  + pb_y[k] * dd_6[k];

        t_12[k] = pb_z[k] * dd_6[k];

        t_13[k] = f_3 * pd_9[k]
                  + pb_x[k] * dd_9[k];

        t_14[k] = pb_z[k] * dd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, t_20, pa_x, pa_y, pa_z, pb_z, pf_0, \
                         pf_5, pf_16, pf_18, pf_19, dd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * pf_5[k];

        t_16[k] = pa_x[k] * pf_16[k];

        t_17[k] = pb_z[k] * dd_9[k];

        t_18[k] = pa_x[k] * pf_18[k];

        t_19[k] = pa_x[k] * pf_19[k];

        t_20[k] = pa_z[k] * pf_0[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_z, pb_x, pb_y, pb_z, pd_0, pd_17, \
                         pf_3, dd_12, dd_14, dd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pb_y[k] * dd_12[k];

        t_22[k] = f_3 * pd_0[k]
                  + pb_z[k] * dd_12[k];

        t_23[k] = pa_z[k] * pf_3[k];

        t_24[k] = pb_y[k] * dd_14[k];

        t_25[k] = f_3 * pd_17[k]
                  + pb_x[k] * dd_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_x, pb_x, pb_y, pf_26, pf_27, pf_29, \
                         dp0_9, dp1_9, dd_17, dd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_x[k] * pf_26[k];

        t_27[k] = pa_x[k] * pf_27[k];

        t_28[k] = pb_y[k] * dd_17[k];

        t_29[k] = pa_x[k] * pf_29[k];

        t_30[k] = f_1 * dp0_9[k]
                  - f_2 * dp1_9[k]
                  + pb_x[k] * dd_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pb_x, pb_y, pb_z, pd_6, dd_18, dd_21, \
                         dd_22, dd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * pd_6[k]
                  + pb_y[k] * dd_18[k];

        t_32[k] = pb_z[k] * dd_18[k];

        t_33[k] = pb_x[k] * dd_21[k];

        t_34[k] = pb_x[k] * dd_22[k];

        t_35[k] = pb_x[k] * dd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pb_z, pd_9, pd_11, dp0_10, dp0_11, \
                         dp1_10, dp1_11, dd_21, dd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * pd_9[k]
                  + f_1 * dp0_10[k]
                  - f_2 * dp1_10[k]
                  + pb_y[k] * dd_21[k];

        t_37[k] = pb_z[k] * dd_21[k];

        t_38[k] = f_0 * pd_11[k]
                  + pb_y[k] * dd_23[k];

        t_39[k] = f_1 * dp0_11[k]
                  - f_2 * dp1_11[k]
                  + pb_z[k] * dd_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, pa_y, pa_z, pb_x, pf_11, pf_20, \
                         pf_22, dd_27, dd_28, dd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * pf_20[k];

        t_41[k] = pa_z[k] * pf_11[k];

        t_42[k] = pa_y[k] * pf_22[k];

        t_43[k] = pb_x[k] * dd_27[k];

        t_44[k] = pb_x[k] * dd_28[k];

        t_45[k] = pb_x[k] * dd_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pa_z, pb_y, pb_z, pd_9, pd_17, pf_16, \
                         pf_29, dd_27, dd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_z[k] * pf_16[k];

        t_47[k] = f_3 * pd_9[k]
                  + pb_z[k] * dd_27[k];

        t_48[k] = f_3 * pd_17[k]
                  + pb_y[k] * dd_29[k];

        t_49[k] = pa_y[k] * pf_29[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, t_55, pb_x, pb_y, pb_z, pd_12, dp0_15, \
                         dp1_15, dd_30, dd_33, dd_34, dd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_1 * dp0_15[k]
                  - f_2 * dp1_15[k]
                  + pb_x[k] * dd_30[k];

        t_51[k] = pb_y[k] * dd_30[k];

        t_52[k] = f_0 * pd_12[k]
                  + pb_z[k] * dd_30[k];

        t_53[k] = pb_x[k] * dd_33[k];

        t_54[k] = pb_x[k] * dd_34[k];

        t_55[k] = pb_x[k] * dd_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pb_z, pd_15, pd_17, dp0_16, dp0_17, \
                         dp1_16, dp1_17, dd_33, dd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * dp0_16[k]
                  - f_2 * dp1_16[k]
                  + pb_y[k] * dd_33[k];

        t_57[k] = f_0 * pd_15[k]
                  + pb_z[k] * dd_33[k];

        t_58[k] = pb_y[k] * dd_35[k];

        t_59[k] = f_0 * pd_17[k]
                  + f_1 * dp0_17[k]
                  - f_2 * dp1_17[k]
                  + pb_z[k] * dd_35[k];
    }
}

}  // namespace simdt2ceri
