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


#include "SimdOverlapVrrRecDD.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_dd_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pp, const size_t pd, const size_t ds,
                          const size_t dp, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 0.5 / p;

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

    const auto *pp_0 = buffer.data(pp + 0);
    const auto *pp_4 = buffer.data(pp + 4);
    const auto *pp_8 = buffer.data(pp + 8);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_9 = buffer.data(pd + 9);
    const auto *pd_11 = buffer.data(pd + 11);
    const auto *pd_12 = buffer.data(pd + 12);
    const auto *pd_15 = buffer.data(pd + 15);
    const auto *pd_17 = buffer.data(pd + 17);

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_3 = buffer.data(ds + 3);
    const auto *ds_5 = buffer.data(ds + 5);

    const auto *dp_0 = buffer.data(dp + 0);
    const auto *dp_1 = buffer.data(dp + 1);
    const auto *dp_2 = buffer.data(dp + 2);
    const auto *dp_3 = buffer.data(dp + 3);
    const auto *dp_4 = buffer.data(dp + 4);
    const auto *dp_6 = buffer.data(dp + 6);
    const auto *dp_8 = buffer.data(dp + 8);
    const auto *dp_9 = buffer.data(dp + 9);
    const auto *dp_10 = buffer.data(dp + 10);
    const auto *dp_11 = buffer.data(dp + 11);
    const auto *dp_13 = buffer.data(dp + 13);
    const auto *dp_14 = buffer.data(dp + 14);
    const auto *dp_15 = buffer.data(dp + 15);
    const auto *dp_16 = buffer.data(dp + 16);
    const auto *dp_17 = buffer.data(dp + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, pp_0, ds_0, dp_0, \
                         dp_1, dp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pp_0[k]
                 + f_1 * ds_0[k]
                 + pb_x[k] * dp_0[k];

        t_1[k] = pb_y[k] * dp_0[k];

        t_2[k] = pb_z[k] * dp_0[k];

        t_3[k] = f_1 * ds_0[k]
                 + pb_y[k] * dp_1[k];

        t_4[k] = pb_y[k] * dp_2[k];

        t_5[k] = f_1 * ds_0[k]
                 + pb_z[k] * dp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pa_x, pa_y, pb_x, pb_z, pp_4, pd_0, \
                         pd_9, pd_11, dp_3, dp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * pd_0[k];

        t_7[k] = f_1 * pp_4[k]
                 + pb_x[k] * dp_4[k];

        t_8[k] = pb_z[k] * dp_3[k];

        t_9[k] = pa_x[k] * pd_9[k];

        t_10[k] = pb_z[k] * dp_4[k];

        t_11[k] = pa_x[k] * pd_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pa_x, pa_z, pb_x, pb_y, pp_8, \
                         pd_0, pd_15, pd_17, dp_6, dp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_z[k] * pd_0[k];

        t_13[k] = pb_y[k] * dp_6[k];

        t_14[k] = f_1 * pp_8[k]
                  + pb_x[k] * dp_8[k];

        t_15[k] = pa_x[k] * pd_15[k];

        t_16[k] = pb_y[k] * dp_8[k];

        t_17[k] = pa_x[k] * pd_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, pb_x, pb_y, pb_z, pp_4, ds_3, \
                         dp_9, dp_10, dp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_1 * ds_3[k]
                  + pb_x[k] * dp_9[k];

        t_19[k] = pb_x[k] * dp_10[k];

        t_20[k] = pb_x[k] * dp_11[k];

        t_21[k] = f_0 * pp_4[k]
                  + f_1 * ds_3[k]
                  + pb_y[k] * dp_10[k];

        t_22[k] = pb_z[k] * dp_10[k];

        t_23[k] = f_1 * ds_3[k]
                  + pb_z[k] * dp_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, pa_y, pa_z, pb_x, pb_y, pp_8, \
                         pd_9, pd_12, pd_17, dp_13, dp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * pd_12[k];

        t_25[k] = pb_x[k] * dp_13[k];

        t_26[k] = pb_x[k] * dp_14[k];

        t_27[k] = pa_z[k] * pd_9[k];

        t_28[k] = f_1 * pp_8[k]
                  + pb_y[k] * dp_14[k];

        t_29[k] = pa_y[k] * pd_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pb_x, pb_y, pb_z, pp_8, ds_5, \
                         dp_15, dp_16, dp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_1 * ds_5[k]
                  + pb_x[k] * dp_15[k];

        t_31[k] = pb_x[k] * dp_16[k];

        t_32[k] = pb_x[k] * dp_17[k];

        t_33[k] = f_1 * ds_5[k]
                  + pb_y[k] * dp_16[k];

        t_34[k] = pb_y[k] * dp_17[k];

        t_35[k] = f_0 * pp_8[k]
                  + f_1 * ds_5[k]
                  + pb_z[k] * dp_17[k];
    }
}

}  // namespace simdovl
