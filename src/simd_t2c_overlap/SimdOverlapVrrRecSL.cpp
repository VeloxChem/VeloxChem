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


#include "SimdOverlapVrrRecSL.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_sl_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                          const size_t si, const size_t sk, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_10 = buffer.data(si + 10);
    const auto *si_12 = buffer.data(si + 12);
    const auto *si_14 = buffer.data(si + 14);
    const auto *si_15 = buffer.data(si + 15);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_20 = buffer.data(si + 20);
    const auto *si_21 = buffer.data(si + 21);
    const auto *si_23 = buffer.data(si + 23);
    const auto *si_24 = buffer.data(si + 24);
    const auto *si_25 = buffer.data(si + 25);
    const auto *si_26 = buffer.data(si + 26);
    const auto *si_27 = buffer.data(si + 27);

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_2 = buffer.data(sk + 2);
    const auto *sk_3 = buffer.data(sk + 3);
    const auto *sk_5 = buffer.data(sk + 5);
    const auto *sk_6 = buffer.data(sk + 6);
    const auto *sk_9 = buffer.data(sk + 9);
    const auto *sk_10 = buffer.data(sk + 10);
    const auto *sk_12 = buffer.data(sk + 12);
    const auto *sk_14 = buffer.data(sk + 14);
    const auto *sk_15 = buffer.data(sk + 15);
    const auto *sk_17 = buffer.data(sk + 17);
    const auto *sk_18 = buffer.data(sk + 18);
    const auto *sk_20 = buffer.data(sk + 20);
    const auto *sk_21 = buffer.data(sk + 21);
    const auto *sk_23 = buffer.data(sk + 23);
    const auto *sk_24 = buffer.data(sk + 24);
    const auto *sk_25 = buffer.data(sk + 25);
    const auto *sk_27 = buffer.data(sk + 27);
    const auto *sk_28 = buffer.data(sk + 28);
    const auto *sk_29 = buffer.data(sk + 29);
    const auto *sk_30 = buffer.data(sk + 30);
    const auto *sk_31 = buffer.data(sk + 31);
    const auto *sk_32 = buffer.data(sk + 32);
    const auto *sk_33 = buffer.data(sk + 33);
    const auto *sk_34 = buffer.data(sk + 34);
    const auto *sk_35 = buffer.data(sk + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, si_0, si_3, si_5, \
                         sk_0, sk_2, sk_3, sk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * si_0[k]
                 + pb_x[k] * sk_0[k];

        t_1[k] = pb_y[k] * sk_0[k];

        t_2[k] = pb_z[k] * sk_0[k];

        t_3[k] = f_1 * si_3[k]
                 + pb_x[k] * sk_3[k];

        t_4[k] = pb_y[k] * sk_2[k];

        t_5[k] = f_1 * si_5[k]
                 + pb_x[k] * sk_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, si_6, si_9, si_10, sk_3, \
                         sk_5, sk_6, sk_9, sk_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * si_6[k]
                 + pb_x[k] * sk_6[k];

        t_7[k] = pb_z[k] * sk_3[k];

        t_8[k] = pb_y[k] * sk_5[k];

        t_9[k] = f_2 * si_9[k]
                 + pb_x[k] * sk_9[k];

        t_10[k] = f_3 * si_10[k]
                  + pb_x[k] * sk_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pb_x, pb_y, pb_z, si_12, si_14, si_15, \
                         sk_6, sk_9, sk_12, sk_14, sk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * sk_6[k];

        t_12[k] = f_3 * si_12[k]
                  + pb_x[k] * sk_12[k];

        t_13[k] = pb_y[k] * sk_9[k];

        t_14[k] = f_3 * si_14[k]
                  + pb_x[k] * sk_14[k];

        t_15[k] = f_4 * si_15[k]
                  + pb_x[k] * sk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pb_x, pb_y, pb_z, si_17, si_18, si_20, \
                         sk_10, sk_14, sk_17, sk_18, sk_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_z[k] * sk_10[k];

        t_17[k] = f_4 * si_17[k]
                  + pb_x[k] * sk_17[k];

        t_18[k] = f_4 * si_18[k]
                  + pb_x[k] * sk_18[k];

        t_19[k] = pb_y[k] * sk_14[k];

        t_20[k] = f_4 * si_20[k]
                  + pb_x[k] * sk_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_z, si_21, si_23, si_24, si_25, \
                         sk_15, sk_21, sk_23, sk_24, sk_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * si_21[k]
                  + pb_x[k] * sk_21[k];

        t_22[k] = pb_z[k] * sk_15[k];

        t_23[k] = f_5 * si_23[k]
                  + pb_x[k] * sk_23[k];

        t_24[k] = f_5 * si_24[k]
                  + pb_x[k] * sk_24[k];

        t_25[k] = f_5 * si_25[k]
                  + pb_x[k] * sk_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, pb_x, pb_y, si_27, sk_20, sk_27, \
                         sk_28, sk_29, sk_30, sk_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_y[k] * sk_20[k];

        t_27[k] = f_5 * si_27[k]
                  + pb_x[k] * sk_27[k];

        t_28[k] = pb_x[k] * sk_28[k];

        t_29[k] = pb_x[k] * sk_29[k];

        t_30[k] = pb_x[k] * sk_30[k];

        t_31[k] = pb_x[k] * sk_31[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, pb_x, pb_y, pb_z, si_21, sk_28, \
                         sk_32, sk_33, sk_34, sk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_x[k] * sk_32[k];

        t_33[k] = pb_x[k] * sk_33[k];

        t_34[k] = pb_x[k] * sk_34[k];

        t_35[k] = pb_x[k] * sk_35[k];

        t_36[k] = f_0 * si_21[k]
                  + pb_y[k] * sk_28[k];

        t_37[k] = pb_z[k] * sk_28[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pb_y, si_23, si_24, si_25, si_26, \
                         si_27, sk_30, sk_31, sk_32, sk_33, sk_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_1 * si_23[k]
                  + pb_y[k] * sk_30[k];

        t_39[k] = f_2 * si_24[k]
                  + pb_y[k] * sk_31[k];

        t_40[k] = f_3 * si_25[k]
                  + pb_y[k] * sk_32[k];

        t_41[k] = f_4 * si_26[k]
                  + pb_y[k] * sk_33[k];

        t_42[k] = f_5 * si_27[k]
                  + pb_y[k] * sk_34[k];
    }

#pragma omp simd aligned(t_43, t_44, pb_y, pb_z, si_27, sk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * sk_35[k];

        t_44[k] = f_0 * si_27[k]
                  + pb_z[k] * sk_35[k];
    }
}

}  // namespace simdovl
