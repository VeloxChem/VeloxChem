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


#include "SimdOverlapVrrRecSK.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_sk_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                          const size_t sh, const size_t si, const size_t ncols,
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_15 = buffer.data(sh + 15);
    const auto *sh_17 = buffer.data(sh + 17);
    const auto *sh_18 = buffer.data(sh + 18);
    const auto *sh_19 = buffer.data(sh + 19);
    const auto *sh_20 = buffer.data(sh + 20);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_2 = buffer.data(si + 2);
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
    const auto *si_22 = buffer.data(si + 22);
    const auto *si_23 = buffer.data(si + 23);
    const auto *si_24 = buffer.data(si + 24);
    const auto *si_25 = buffer.data(si + 25);
    const auto *si_26 = buffer.data(si + 26);
    const auto *si_27 = buffer.data(si + 27);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, sh_0, sh_3, sh_5, \
                         si_0, si_2, si_3, si_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh_0[k]
                 + pb_x[k] * si_0[k];

        t_1[k] = pb_y[k] * si_0[k];

        t_2[k] = pb_z[k] * si_0[k];

        t_3[k] = f_1 * sh_3[k]
                 + pb_x[k] * si_3[k];

        t_4[k] = pb_y[k] * si_2[k];

        t_5[k] = f_1 * sh_5[k]
                 + pb_x[k] * si_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, sh_6, sh_9, sh_10, si_3, \
                         si_5, si_6, si_9, si_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * sh_6[k]
                 + pb_x[k] * si_6[k];

        t_7[k] = pb_z[k] * si_3[k];

        t_8[k] = pb_y[k] * si_5[k];

        t_9[k] = f_2 * sh_9[k]
                 + pb_x[k] * si_9[k];

        t_10[k] = f_3 * sh_10[k]
                  + pb_x[k] * si_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pb_x, pb_y, pb_z, sh_12, sh_14, sh_15, \
                         si_6, si_9, si_12, si_14, si_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * si_6[k];

        t_12[k] = f_3 * sh_12[k]
                  + pb_x[k] * si_12[k];

        t_13[k] = pb_y[k] * si_9[k];

        t_14[k] = f_3 * sh_14[k]
                  + pb_x[k] * si_14[k];

        t_15[k] = f_4 * sh_15[k]
                  + pb_x[k] * si_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pb_x, pb_y, pb_z, sh_17, sh_18, sh_20, \
                         si_10, si_14, si_17, si_18, si_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_z[k] * si_10[k];

        t_17[k] = f_4 * sh_17[k]
                  + pb_x[k] * si_17[k];

        t_18[k] = f_4 * sh_18[k]
                  + pb_x[k] * si_18[k];

        t_19[k] = pb_y[k] * si_14[k];

        t_20[k] = f_4 * sh_20[k]
                  + pb_x[k] * si_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, t_26, t_27, pb_x, si_21, si_22, si_23, \
                         si_24, si_25, si_26, si_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pb_x[k] * si_21[k];

        t_22[k] = pb_x[k] * si_22[k];

        t_23[k] = pb_x[k] * si_23[k];

        t_24[k] = pb_x[k] * si_24[k];

        t_25[k] = pb_x[k] * si_25[k];

        t_26[k] = pb_x[k] * si_26[k];

        t_27[k] = pb_x[k] * si_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pb_y, pb_z, sh_15, sh_17, sh_18, sh_19, \
                         si_21, si_23, si_24, si_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * sh_15[k]
                  + pb_y[k] * si_21[k];

        t_29[k] = pb_z[k] * si_21[k];

        t_30[k] = f_1 * sh_17[k]
                  + pb_y[k] * si_23[k];

        t_31[k] = f_2 * sh_18[k]
                  + pb_y[k] * si_24[k];

        t_32[k] = f_3 * sh_19[k]
                  + pb_y[k] * si_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_y, pb_z, sh_20, si_26, \
                         si_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_4 * sh_20[k]
                  + pb_y[k] * si_26[k];

        t_34[k] = pb_y[k] * si_27[k];

        t_35[k] = f_0 * sh_20[k]
                  + pb_z[k] * si_27[k];
    }
}

}  // namespace simdovl
