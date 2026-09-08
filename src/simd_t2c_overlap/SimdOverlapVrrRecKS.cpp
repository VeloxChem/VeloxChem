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


#include "SimdOverlapVrrRecKS.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_ks_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t hs, const size_t is, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 2.0 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);
    const auto *hs_17 = buffer.data(hs + 17);
    const auto *hs_18 = buffer.data(hs + 18);
    const auto *hs_19 = buffer.data(hs + 19);
    const auto *hs_20 = buffer.data(hs + 20);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_21 = buffer.data(is + 21);
    const auto *is_22 = buffer.data(is + 22);
    const auto *is_23 = buffer.data(is + 23);
    const auto *is_24 = buffer.data(is + 24);
    const auto *is_25 = buffer.data(is + 25);
    const auto *is_26 = buffer.data(is + 26);
    const auto *is_27 = buffer.data(is + 27);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_x, pa_y, pa_z, hs_0, hs_3, hs_5, \
                         is_0, is_2, is_3, is_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = pa_y[k] * is_0[k];

        t_2[k] = pa_z[k] * is_0[k];

        t_3[k] = f_1 * hs_3[k]
                 + pa_x[k] * is_3[k];

        t_4[k] = pa_y[k] * is_2[k];

        t_5[k] = f_1 * hs_5[k]
                 + pa_x[k] * is_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, hs_6, hs_9, hs_10, is_3, \
                         is_5, is_6, is_9, is_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * hs_6[k]
                 + pa_x[k] * is_6[k];

        t_7[k] = pa_z[k] * is_3[k];

        t_8[k] = pa_y[k] * is_5[k];

        t_9[k] = f_2 * hs_9[k]
                 + pa_x[k] * is_9[k];

        t_10[k] = f_3 * hs_10[k]
                  + pa_x[k] * is_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_x, pa_y, pa_z, hs_12, hs_14, hs_15, \
                         is_6, is_9, is_12, is_14, is_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_z[k] * is_6[k];

        t_12[k] = f_3 * hs_12[k]
                  + pa_x[k] * is_12[k];

        t_13[k] = pa_y[k] * is_9[k];

        t_14[k] = f_3 * hs_14[k]
                  + pa_x[k] * is_14[k];

        t_15[k] = f_4 * hs_15[k]
                  + pa_x[k] * is_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_x, pa_y, pa_z, hs_17, hs_18, hs_20, \
                         is_10, is_14, is_17, is_18, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_z[k] * is_10[k];

        t_17[k] = f_4 * hs_17[k]
                  + pa_x[k] * is_17[k];

        t_18[k] = f_4 * hs_18[k]
                  + pa_x[k] * is_18[k];

        t_19[k] = pa_y[k] * is_14[k];

        t_20[k] = f_4 * hs_20[k]
                  + pa_x[k] * is_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, t_26, t_27, pa_x, is_21, is_22, is_23, \
                         is_24, is_25, is_26, is_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_x[k] * is_21[k];

        t_22[k] = pa_x[k] * is_22[k];

        t_23[k] = pa_x[k] * is_23[k];

        t_24[k] = pa_x[k] * is_24[k];

        t_25[k] = pa_x[k] * is_25[k];

        t_26[k] = pa_x[k] * is_26[k];

        t_27[k] = pa_x[k] * is_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_y, pa_z, hs_15, hs_17, hs_18, hs_19, \
                         is_21, is_23, is_24, is_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * hs_15[k]
                  + pa_y[k] * is_21[k];

        t_29[k] = pa_z[k] * is_21[k];

        t_30[k] = f_1 * hs_17[k]
                  + pa_y[k] * is_23[k];

        t_31[k] = f_2 * hs_18[k]
                  + pa_y[k] * is_24[k];

        t_32[k] = f_3 * hs_19[k]
                  + pa_y[k] * is_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pa_z, hs_20, is_26, \
                         is_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_4 * hs_20[k]
                  + pa_y[k] * is_26[k];

        t_34[k] = pa_y[k] * is_27[k];

        t_35[k] = f_0 * hs_20[k]
                  + pa_z[k] * is_27[k];
    }
}

}  // namespace simdovl
