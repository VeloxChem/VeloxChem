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


#include "SimdElectronRepulsionVrrRecLS.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_ls_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t is0, const size_t is1, const size_t ks,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / alpha;
    const auto f_1 = 3.5 * beta / (alpha * p);
    const auto f_2 = 2.5 / alpha;
    const auto f_3 = 2.5 * beta / (alpha * p);
    const auto f_4 = 2.0 / alpha;
    const auto f_5 = 2.0 * beta / (alpha * p);
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 0.5 / alpha;
    const auto f_11 = 0.5 * beta / (alpha * p);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);
    const auto *is0_15 = buffer.data(is0 + 15);
    const auto *is0_17 = buffer.data(is0 + 17);
    const auto *is0_18 = buffer.data(is0 + 18);
    const auto *is0_20 = buffer.data(is0 + 20);
    const auto *is0_21 = buffer.data(is0 + 21);
    const auto *is0_23 = buffer.data(is0 + 23);
    const auto *is0_24 = buffer.data(is0 + 24);
    const auto *is0_25 = buffer.data(is0 + 25);
    const auto *is0_26 = buffer.data(is0 + 26);
    const auto *is0_27 = buffer.data(is0 + 27);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);
    const auto *is1_15 = buffer.data(is1 + 15);
    const auto *is1_17 = buffer.data(is1 + 17);
    const auto *is1_18 = buffer.data(is1 + 18);
    const auto *is1_20 = buffer.data(is1 + 20);
    const auto *is1_21 = buffer.data(is1 + 21);
    const auto *is1_23 = buffer.data(is1 + 23);
    const auto *is1_24 = buffer.data(is1 + 24);
    const auto *is1_25 = buffer.data(is1 + 25);
    const auto *is1_26 = buffer.data(is1 + 26);
    const auto *is1_27 = buffer.data(is1 + 27);

    const auto *ks_0 = buffer.data(ks + 0);
    const auto *ks_2 = buffer.data(ks + 2);
    const auto *ks_3 = buffer.data(ks + 3);
    const auto *ks_5 = buffer.data(ks + 5);
    const auto *ks_6 = buffer.data(ks + 6);
    const auto *ks_9 = buffer.data(ks + 9);
    const auto *ks_10 = buffer.data(ks + 10);
    const auto *ks_12 = buffer.data(ks + 12);
    const auto *ks_14 = buffer.data(ks + 14);
    const auto *ks_15 = buffer.data(ks + 15);
    const auto *ks_17 = buffer.data(ks + 17);
    const auto *ks_18 = buffer.data(ks + 18);
    const auto *ks_20 = buffer.data(ks + 20);
    const auto *ks_21 = buffer.data(ks + 21);
    const auto *ks_23 = buffer.data(ks + 23);
    const auto *ks_24 = buffer.data(ks + 24);
    const auto *ks_25 = buffer.data(ks + 25);
    const auto *ks_27 = buffer.data(ks + 27);
    const auto *ks_28 = buffer.data(ks + 28);
    const auto *ks_29 = buffer.data(ks + 29);
    const auto *ks_30 = buffer.data(ks + 30);
    const auto *ks_31 = buffer.data(ks + 31);
    const auto *ks_32 = buffer.data(ks + 32);
    const auto *ks_33 = buffer.data(ks + 33);
    const auto *ks_34 = buffer.data(ks + 34);
    const auto *ks_35 = buffer.data(ks + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pa_z, is0_0, is0_3, is1_0, \
                         is1_3, ks_0, ks_2, ks_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * is0_0[k]
                 - f_1 * is1_0[k]
                 + pa_x[k] * ks_0[k];

        t_1[k] = pa_y[k] * ks_0[k];

        t_2[k] = pa_z[k] * ks_0[k];

        t_3[k] = f_2 * is0_3[k]
                 - f_3 * is1_3[k]
                 + pa_x[k] * ks_3[k];

        t_4[k] = pa_y[k] * ks_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_x, pa_y, pa_z, is0_5, is0_6, is1_5, is1_6, \
                         ks_3, ks_5, ks_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_2 * is0_5[k]
                 - f_3 * is1_5[k]
                 + pa_x[k] * ks_5[k];

        t_6[k] = f_4 * is0_6[k]
                 - f_5 * is1_6[k]
                 + pa_x[k] * ks_6[k];

        t_7[k] = pa_z[k] * ks_3[k];

        t_8[k] = pa_y[k] * ks_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pa_z, is0_9, is0_10, is0_12, is1_9, \
                         is1_10, is1_12, ks_6, ks_9, ks_10, ks_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * is0_9[k]
                 - f_5 * is1_9[k]
                 + pa_x[k] * ks_9[k];

        t_10[k] = f_6 * is0_10[k]
                  - f_7 * is1_10[k]
                  + pa_x[k] * ks_10[k];

        t_11[k] = pa_z[k] * ks_6[k];

        t_12[k] = f_6 * is0_12[k]
                  - f_7 * is1_12[k]
                  + pa_x[k] * ks_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pa_y, pa_z, is0_14, is0_15, is1_14, \
                         is1_15, ks_9, ks_10, ks_14, ks_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * ks_9[k];

        t_14[k] = f_6 * is0_14[k]
                  - f_7 * is1_14[k]
                  + pa_x[k] * ks_14[k];

        t_15[k] = f_8 * is0_15[k]
                  - f_9 * is1_15[k]
                  + pa_x[k] * ks_15[k];

        t_16[k] = pa_z[k] * ks_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_x, pa_y, is0_17, is0_18, is0_20, is1_17, \
                         is1_18, is1_20, ks_14, ks_17, ks_18, ks_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * is0_17[k]
                  - f_9 * is1_17[k]
                  + pa_x[k] * ks_17[k];

        t_18[k] = f_8 * is0_18[k]
                  - f_9 * is1_18[k]
                  + pa_x[k] * ks_18[k];

        t_19[k] = pa_y[k] * ks_14[k];

        t_20[k] = f_8 * is0_20[k]
                  - f_9 * is1_20[k]
                  + pa_x[k] * ks_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_x, pa_z, is0_21, is0_23, is0_24, is1_21, \
                         is1_23, is1_24, ks_15, ks_21, ks_23, ks_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_10 * is0_21[k]
                  - f_11 * is1_21[k]
                  + pa_x[k] * ks_21[k];

        t_22[k] = pa_z[k] * ks_15[k];

        t_23[k] = f_10 * is0_23[k]
                  - f_11 * is1_23[k]
                  + pa_x[k] * ks_23[k];

        t_24[k] = f_10 * is0_24[k]
                  - f_11 * is1_24[k]
                  + pa_x[k] * ks_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_x, pa_y, is0_25, is0_27, is1_25, \
                         is1_27, ks_20, ks_25, ks_27, ks_28, ks_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_10 * is0_25[k]
                  - f_11 * is1_25[k]
                  + pa_x[k] * ks_25[k];

        t_26[k] = pa_y[k] * ks_20[k];

        t_27[k] = f_10 * is0_27[k]
                  - f_11 * is1_27[k]
                  + pa_x[k] * ks_27[k];

        t_28[k] = pa_x[k] * ks_28[k];

        t_29[k] = pa_x[k] * ks_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pa_x, ks_30, ks_31, ks_32, ks_33, \
                         ks_34, ks_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_x[k] * ks_30[k];

        t_31[k] = pa_x[k] * ks_31[k];

        t_32[k] = pa_x[k] * ks_32[k];

        t_33[k] = pa_x[k] * ks_33[k];

        t_34[k] = pa_x[k] * ks_34[k];

        t_35[k] = pa_x[k] * ks_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_y, pa_z, is0_21, is0_23, is0_24, is1_21, \
                         is1_23, is1_24, ks_28, ks_30, ks_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * is0_21[k]
                  - f_1 * is1_21[k]
                  + pa_y[k] * ks_28[k];

        t_37[k] = pa_z[k] * ks_28[k];

        t_38[k] = f_2 * is0_23[k]
                  - f_3 * is1_23[k]
                  + pa_y[k] * ks_30[k];

        t_39[k] = f_4 * is0_24[k]
                  - f_5 * is1_24[k]
                  + pa_y[k] * ks_31[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, is0_25, is0_26, is0_27, is1_25, is1_26, \
                         is1_27, ks_32, ks_33, ks_34, ks_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_6 * is0_25[k]
                  - f_7 * is1_25[k]
                  + pa_y[k] * ks_32[k];

        t_41[k] = f_8 * is0_26[k]
                  - f_9 * is1_26[k]
                  + pa_y[k] * ks_33[k];

        t_42[k] = f_10 * is0_27[k]
                  - f_11 * is1_27[k]
                  + pa_y[k] * ks_34[k];

        t_43[k] = pa_y[k] * ks_35[k];
    }

#pragma omp simd aligned(t_44, pa_z, is0_27, is1_27, ks_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * is0_27[k]
                  - f_1 * is1_27[k]
                  + pa_z[k] * ks_35[k];
    }
}

}  // namespace simdt2ceri
