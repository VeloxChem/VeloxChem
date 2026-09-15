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


#include "SimdElectronRepulsionVrrRecMS.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_ms_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t ks0, const size_t ks1, const size_t ls,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / alpha;
    const auto f_1 = 4.0 * beta / (alpha * p);
    const auto f_2 = 3.0 / alpha;
    const auto f_3 = 3.0 * beta / (alpha * p);
    const auto f_4 = 2.5 / alpha;
    const auto f_5 = 2.5 * beta / (alpha * p);
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.5 / alpha;
    const auto f_9 = 1.5 * beta / (alpha * p);
    const auto f_10 = 1.0 / alpha;
    const auto f_11 = beta / (alpha * p);
    const auto f_12 = 0.5 / alpha;
    const auto f_13 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ks0_0 = buffer.data(ks0 + 0);
    const auto *ks0_3 = buffer.data(ks0 + 3);
    const auto *ks0_5 = buffer.data(ks0 + 5);
    const auto *ks0_6 = buffer.data(ks0 + 6);
    const auto *ks0_9 = buffer.data(ks0 + 9);
    const auto *ks0_10 = buffer.data(ks0 + 10);
    const auto *ks0_12 = buffer.data(ks0 + 12);
    const auto *ks0_14 = buffer.data(ks0 + 14);
    const auto *ks0_15 = buffer.data(ks0 + 15);
    const auto *ks0_17 = buffer.data(ks0 + 17);
    const auto *ks0_18 = buffer.data(ks0 + 18);
    const auto *ks0_20 = buffer.data(ks0 + 20);
    const auto *ks0_21 = buffer.data(ks0 + 21);
    const auto *ks0_23 = buffer.data(ks0 + 23);
    const auto *ks0_24 = buffer.data(ks0 + 24);
    const auto *ks0_25 = buffer.data(ks0 + 25);
    const auto *ks0_27 = buffer.data(ks0 + 27);
    const auto *ks0_28 = buffer.data(ks0 + 28);
    const auto *ks0_30 = buffer.data(ks0 + 30);
    const auto *ks0_31 = buffer.data(ks0 + 31);
    const auto *ks0_32 = buffer.data(ks0 + 32);
    const auto *ks0_33 = buffer.data(ks0 + 33);
    const auto *ks0_34 = buffer.data(ks0 + 34);
    const auto *ks0_35 = buffer.data(ks0 + 35);

    const auto *ks1_0 = buffer.data(ks1 + 0);
    const auto *ks1_3 = buffer.data(ks1 + 3);
    const auto *ks1_5 = buffer.data(ks1 + 5);
    const auto *ks1_6 = buffer.data(ks1 + 6);
    const auto *ks1_9 = buffer.data(ks1 + 9);
    const auto *ks1_10 = buffer.data(ks1 + 10);
    const auto *ks1_12 = buffer.data(ks1 + 12);
    const auto *ks1_14 = buffer.data(ks1 + 14);
    const auto *ks1_15 = buffer.data(ks1 + 15);
    const auto *ks1_17 = buffer.data(ks1 + 17);
    const auto *ks1_18 = buffer.data(ks1 + 18);
    const auto *ks1_20 = buffer.data(ks1 + 20);
    const auto *ks1_21 = buffer.data(ks1 + 21);
    const auto *ks1_23 = buffer.data(ks1 + 23);
    const auto *ks1_24 = buffer.data(ks1 + 24);
    const auto *ks1_25 = buffer.data(ks1 + 25);
    const auto *ks1_27 = buffer.data(ks1 + 27);
    const auto *ks1_28 = buffer.data(ks1 + 28);
    const auto *ks1_30 = buffer.data(ks1 + 30);
    const auto *ks1_31 = buffer.data(ks1 + 31);
    const auto *ks1_32 = buffer.data(ks1 + 32);
    const auto *ks1_33 = buffer.data(ks1 + 33);
    const auto *ks1_34 = buffer.data(ks1 + 34);
    const auto *ks1_35 = buffer.data(ks1 + 35);

    const auto *ls_0 = buffer.data(ls + 0);
    const auto *ls_2 = buffer.data(ls + 2);
    const auto *ls_3 = buffer.data(ls + 3);
    const auto *ls_5 = buffer.data(ls + 5);
    const auto *ls_6 = buffer.data(ls + 6);
    const auto *ls_9 = buffer.data(ls + 9);
    const auto *ls_10 = buffer.data(ls + 10);
    const auto *ls_12 = buffer.data(ls + 12);
    const auto *ls_14 = buffer.data(ls + 14);
    const auto *ls_15 = buffer.data(ls + 15);
    const auto *ls_17 = buffer.data(ls + 17);
    const auto *ls_18 = buffer.data(ls + 18);
    const auto *ls_20 = buffer.data(ls + 20);
    const auto *ls_21 = buffer.data(ls + 21);
    const auto *ls_23 = buffer.data(ls + 23);
    const auto *ls_24 = buffer.data(ls + 24);
    const auto *ls_25 = buffer.data(ls + 25);
    const auto *ls_27 = buffer.data(ls + 27);
    const auto *ls_28 = buffer.data(ls + 28);
    const auto *ls_30 = buffer.data(ls + 30);
    const auto *ls_31 = buffer.data(ls + 31);
    const auto *ls_32 = buffer.data(ls + 32);
    const auto *ls_33 = buffer.data(ls + 33);
    const auto *ls_35 = buffer.data(ls + 35);
    const auto *ls_36 = buffer.data(ls + 36);
    const auto *ls_37 = buffer.data(ls + 37);
    const auto *ls_38 = buffer.data(ls + 38);
    const auto *ls_39 = buffer.data(ls + 39);
    const auto *ls_40 = buffer.data(ls + 40);
    const auto *ls_41 = buffer.data(ls + 41);
    const auto *ls_42 = buffer.data(ls + 42);
    const auto *ls_43 = buffer.data(ls + 43);
    const auto *ls_44 = buffer.data(ls + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, ks0_0, ks0_3, ks1_0, \
                         ks1_3, ls_0, ls_2, ls_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ks0_0[k]
                 - f_1 * ks1_0[k]
                 + pa_x[k] * ls_0[k];

        t_1[k] = pa_y[k] * ls_0[k];

        t_2[k] = pa_z[k] * ls_0[k];

        t_3[k] = f_2 * ks0_3[k]
                 - f_3 * ks1_3[k]
                 + pa_x[k] * ls_3[k];

        t_4[k] = pa_y[k] * ls_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, ks0_5, ks0_6, ks1_5, ks1_6, \
                         ls_3, ls_5, ls_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_2 * ks0_5[k]
                 - f_3 * ks1_5[k]
                 + pa_x[k] * ls_5[k];

        t_6[k] = f_4 * ks0_6[k]
                 - f_5 * ks1_6[k]
                 + pa_x[k] * ls_6[k];

        t_7[k] = pa_z[k] * ls_3[k];

        t_8[k] = pa_y[k] * ls_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pa_z, ks0_9, ks0_10, ks0_12, ks1_9, \
                         ks1_10, ks1_12, ls_6, ls_9, ls_10, ls_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * ks0_9[k]
                 - f_5 * ks1_9[k]
                 + pa_x[k] * ls_9[k];

        t_10[k] = f_6 * ks0_10[k]
                  - f_7 * ks1_10[k]
                  + pa_x[k] * ls_10[k];

        t_11[k] = pa_z[k] * ls_6[k];

        t_12[k] = f_6 * ks0_12[k]
                  - f_7 * ks1_12[k]
                  + pa_x[k] * ls_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pa_y, pa_z, ks0_14, ks0_15, ks1_14, \
                         ks1_15, ls_9, ls_10, ls_14, ls_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * ls_9[k];

        t_14[k] = f_6 * ks0_14[k]
                  - f_7 * ks1_14[k]
                  + pa_x[k] * ls_14[k];

        t_15[k] = f_8 * ks0_15[k]
                  - f_9 * ks1_15[k]
                  + pa_x[k] * ls_15[k];

        t_16[k] = pa_z[k] * ls_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_x, pa_y, ks0_17, ks0_18, ks0_20, ks1_17, \
                         ks1_18, ks1_20, ls_14, ls_17, ls_18, ls_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * ks0_17[k]
                  - f_9 * ks1_17[k]
                  + pa_x[k] * ls_17[k];

        t_18[k] = f_8 * ks0_18[k]
                  - f_9 * ks1_18[k]
                  + pa_x[k] * ls_18[k];

        t_19[k] = pa_y[k] * ls_14[k];

        t_20[k] = f_8 * ks0_20[k]
                  - f_9 * ks1_20[k]
                  + pa_x[k] * ls_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_x, pa_z, ks0_21, ks0_23, ks0_24, ks1_21, \
                         ks1_23, ks1_24, ls_15, ls_21, ls_23, ls_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_10 * ks0_21[k]
                  - f_11 * ks1_21[k]
                  + pa_x[k] * ls_21[k];

        t_22[k] = pa_z[k] * ls_15[k];

        t_23[k] = f_10 * ks0_23[k]
                  - f_11 * ks1_23[k]
                  + pa_x[k] * ls_23[k];

        t_24[k] = f_10 * ks0_24[k]
                  - f_11 * ks1_24[k]
                  + pa_x[k] * ls_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pa_y, ks0_25, ks0_27, ks0_28, ks1_25, \
                         ks1_27, ks1_28, ls_20, ls_25, ls_27, ls_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_10 * ks0_25[k]
                  - f_11 * ks1_25[k]
                  + pa_x[k] * ls_25[k];

        t_26[k] = pa_y[k] * ls_20[k];

        t_27[k] = f_10 * ks0_27[k]
                  - f_11 * ks1_27[k]
                  + pa_x[k] * ls_27[k];

        t_28[k] = f_12 * ks0_28[k]
                  - f_13 * ks1_28[k]
                  + pa_x[k] * ls_28[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_x, pa_z, ks0_30, ks0_31, ks0_32, ks1_30, \
                         ks1_31, ks1_32, ls_21, ls_30, ls_31, ls_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_z[k] * ls_21[k];

        t_30[k] = f_12 * ks0_30[k]
                  - f_13 * ks1_30[k]
                  + pa_x[k] * ls_30[k];

        t_31[k] = f_12 * ks0_31[k]
                  - f_13 * ks1_31[k]
                  + pa_x[k] * ls_31[k];

        t_32[k] = f_12 * ks0_32[k]
                  - f_13 * ks1_32[k]
                  + pa_x[k] * ls_32[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_x, pa_y, ks0_33, ks0_35, ks1_33, \
                         ks1_35, ls_27, ls_33, ls_35, ls_36, ls_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_12 * ks0_33[k]
                  - f_13 * ks1_33[k]
                  + pa_x[k] * ls_33[k];

        t_34[k] = pa_y[k] * ls_27[k];

        t_35[k] = f_12 * ks0_35[k]
                  - f_13 * ks1_35[k]
                  + pa_x[k] * ls_35[k];

        t_36[k] = pa_x[k] * ls_36[k];

        t_37[k] = pa_x[k] * ls_37[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, t_44, pa_x, ls_38, ls_39, ls_40, \
                         ls_41, ls_42, ls_43, ls_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_x[k] * ls_38[k];

        t_39[k] = pa_x[k] * ls_39[k];

        t_40[k] = pa_x[k] * ls_40[k];

        t_41[k] = pa_x[k] * ls_41[k];

        t_42[k] = pa_x[k] * ls_42[k];

        t_43[k] = pa_x[k] * ls_43[k];

        t_44[k] = pa_x[k] * ls_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_y, pa_z, ks0_28, ks0_30, ks0_31, ks1_28, \
                         ks1_30, ks1_31, ls_36, ls_38, ls_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_0 * ks0_28[k]
                  - f_1 * ks1_28[k]
                  + pa_y[k] * ls_36[k];

        t_46[k] = pa_z[k] * ls_36[k];

        t_47[k] = f_2 * ks0_30[k]
                  - f_3 * ks1_30[k]
                  + pa_y[k] * ls_38[k];

        t_48[k] = f_4 * ks0_31[k]
                  - f_5 * ks1_31[k]
                  + pa_y[k] * ls_39[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_y, ks0_32, ks0_33, ks0_34, ks1_32, ks1_33, \
                         ks1_34, ls_40, ls_41, ls_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_6 * ks0_32[k]
                  - f_7 * ks1_32[k]
                  + pa_y[k] * ls_40[k];

        t_50[k] = f_8 * ks0_33[k]
                  - f_9 * ks1_33[k]
                  + pa_y[k] * ls_41[k];

        t_51[k] = f_10 * ks0_34[k]
                  - f_11 * ks1_34[k]
                  + pa_y[k] * ls_42[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_y, pa_z, ks0_35, ks1_35, ls_43, \
                         ls_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_12 * ks0_35[k]
                  - f_13 * ks1_35[k]
                  + pa_y[k] * ls_43[k];

        t_53[k] = pa_y[k] * ls_44[k];

        t_54[k] = f_0 * ks0_35[k]
                  - f_1 * ks1_35[k]
                  + pa_z[k] * ls_44[k];
    }
}

}  // namespace simdt2ceri
