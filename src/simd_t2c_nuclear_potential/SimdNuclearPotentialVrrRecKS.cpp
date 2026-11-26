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


#include "SimdNuclearPotentialVrrRecKS.hpp"

#include "SimdAlign.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_prim_ks_nuclear_potential_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                    const size_t pc, const size_t hs0, const size_t hs1,
                                    const size_t is0, const size_t is1, const size_t ncols,
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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_12 = buffer.data(hs0 + 12);
    const auto *hs0_14 = buffer.data(hs0 + 14);
    const auto *hs0_15 = buffer.data(hs0 + 15);
    const auto *hs0_17 = buffer.data(hs0 + 17);
    const auto *hs0_18 = buffer.data(hs0 + 18);
    const auto *hs0_19 = buffer.data(hs0 + 19);
    const auto *hs0_20 = buffer.data(hs0 + 20);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_12 = buffer.data(hs1 + 12);
    const auto *hs1_14 = buffer.data(hs1 + 14);
    const auto *hs1_15 = buffer.data(hs1 + 15);
    const auto *hs1_17 = buffer.data(hs1 + 17);
    const auto *hs1_18 = buffer.data(hs1 + 18);
    const auto *hs1_19 = buffer.data(hs1 + 19);
    const auto *hs1_20 = buffer.data(hs1 + 20);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_2 = buffer.data(is0 + 2);
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
    const auto *is0_22 = buffer.data(is0 + 22);
    const auto *is0_23 = buffer.data(is0 + 23);
    const auto *is0_24 = buffer.data(is0 + 24);
    const auto *is0_25 = buffer.data(is0 + 25);
    const auto *is0_26 = buffer.data(is0 + 26);
    const auto *is0_27 = buffer.data(is0 + 27);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_2 = buffer.data(is1 + 2);
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
    const auto *is1_22 = buffer.data(is1 + 22);
    const auto *is1_23 = buffer.data(is1 + 23);
    const auto *is1_24 = buffer.data(is1 + 24);
    const auto *is1_25 = buffer.data(is1 + 25);
    const auto *is1_26 = buffer.data(is1 + 26);
    const auto *is1_27 = buffer.data(is1 + 27);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, hs0_0, hs1_0, \
                         is0_0, is1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_0 * hs1_0[k]
                 + pa_x[k] * is0_0[k]
                 - pc_x[k] * is1_0[k];

        t_1[k] = pa_y[k] * is0_0[k]
                 - pc_y[k] * is1_0[k];

        t_2[k] = pa_z[k] * is0_0[k]
                 - pc_z[k] * is1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pa_x, pa_y, pc_x, pc_y, hs0_3, hs1_3, is0_2, is0_3, is1_2, \
                         is1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_1 * hs0_3[k]
                 - f_1 * hs1_3[k]
                 + pa_x[k] * is0_3[k]
                 - pc_x[k] * is1_3[k];

        t_4[k] = pa_y[k] * is0_2[k]
                 - pc_y[k] * is1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pa_x, pc_x, hs0_5, hs0_6, hs1_5, hs1_6, is0_5, is0_6, \
                         is1_5, is1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * hs0_5[k]
                 - f_1 * hs1_5[k]
                 + pa_x[k] * is0_5[k]
                 - pc_x[k] * is1_5[k];

        t_6[k] = f_2 * hs0_6[k]
                 - f_2 * hs1_6[k]
                 + pa_x[k] * is0_6[k]
                 - pc_x[k] * is1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pa_y, pa_z, pc_y, pc_z, is0_3, is0_5, is1_3, \
                         is1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_z[k] * is0_3[k]
                 - pc_z[k] * is1_3[k];

        t_8[k] = pa_y[k] * is0_5[k]
                 - pc_y[k] * is1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, pa_x, pc_x, hs0_9, hs0_10, hs1_9, hs1_10, is0_9, is0_10, \
                         is1_9, is1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * hs0_9[k]
                 - f_2 * hs1_9[k]
                 + pa_x[k] * is0_9[k]
                 - pc_x[k] * is1_9[k];

        t_10[k] = f_3 * hs0_10[k]
                  - f_3 * hs1_10[k]
                  + pa_x[k] * is0_10[k]
                  - pc_x[k] * is1_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pa_x, pa_z, pc_x, pc_z, hs0_12, hs1_12, is0_6, is0_12, \
                         is1_6, is1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_z[k] * is0_6[k]
                  - pc_z[k] * is1_6[k];

        t_12[k] = f_3 * hs0_12[k]
                  - f_3 * hs1_12[k]
                  + pa_x[k] * is0_12[k]
                  - pc_x[k] * is1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_x, pa_y, pc_x, pc_y, hs0_14, hs1_14, is0_9, is0_14, \
                         is1_9, is1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * is0_9[k]
                  - pc_y[k] * is1_9[k];

        t_14[k] = f_3 * hs0_14[k]
                  - f_3 * hs1_14[k]
                  + pa_x[k] * is0_14[k]
                  - pc_x[k] * is1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pa_x, pa_z, pc_x, pc_z, hs0_15, hs1_15, is0_10, is0_15, \
                         is1_10, is1_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * hs0_15[k]
                  - f_4 * hs1_15[k]
                  + pa_x[k] * is0_15[k]
                  - pc_x[k] * is1_15[k];

        t_16[k] = pa_z[k] * is0_10[k]
                  - pc_z[k] * is1_10[k];
    }

#pragma omp simd aligned(t_17, t_18, pa_x, pc_x, hs0_17, hs0_18, hs1_17, hs1_18, is0_17, \
                         is0_18, is1_17, is1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_4 * hs0_17[k]
                  - f_4 * hs1_17[k]
                  + pa_x[k] * is0_17[k]
                  - pc_x[k] * is1_17[k];

        t_18[k] = f_4 * hs0_18[k]
                  - f_4 * hs1_18[k]
                  + pa_x[k] * is0_18[k]
                  - pc_x[k] * is1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pa_y, pc_x, pc_y, hs0_20, hs1_20, is0_14, \
                         is0_20, is0_21, is1_14, is1_20, is1_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * is0_14[k]
                  - pc_y[k] * is1_14[k];

        t_20[k] = f_4 * hs0_20[k]
                  - f_4 * hs1_20[k]
                  + pa_x[k] * is0_20[k]
                  - pc_x[k] * is1_20[k];

        t_21[k] = pa_x[k] * is0_21[k]
                  - pc_x[k] * is1_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pc_x, is0_22, is0_23, is0_24, is0_25, \
                         is1_22, is1_23, is1_24, is1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pa_x[k] * is0_22[k]
                  - pc_x[k] * is1_22[k];

        t_23[k] = pa_x[k] * is0_23[k]
                  - pc_x[k] * is1_23[k];

        t_24[k] = pa_x[k] * is0_24[k]
                  - pc_x[k] * is1_24[k];

        t_25[k] = pa_x[k] * is0_25[k]
                  - pc_x[k] * is1_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_x, pa_y, pc_x, pc_y, hs0_15, hs1_15, is0_21, \
                         is0_26, is0_27, is1_21, is1_26, is1_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_x[k] * is0_26[k]
                  - pc_x[k] * is1_26[k];

        t_27[k] = pa_x[k] * is0_27[k]
                  - pc_x[k] * is1_27[k];

        t_28[k] = f_0 * hs0_15[k]
                  - f_0 * hs1_15[k]
                  + pa_y[k] * is0_21[k]
                  - pc_y[k] * is1_21[k];
    }

#pragma omp simd aligned(t_29, t_30, pa_y, pa_z, pc_y, pc_z, hs0_17, hs1_17, is0_21, is0_23, \
                         is1_21, is1_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_z[k] * is0_21[k]
                  - pc_z[k] * is1_21[k];

        t_30[k] = f_1 * hs0_17[k]
                  - f_1 * hs1_17[k]
                  + pa_y[k] * is0_23[k]
                  - pc_y[k] * is1_23[k];
    }

#pragma omp simd aligned(t_31, t_32, pa_y, pc_y, hs0_18, hs0_19, hs1_18, hs1_19, is0_24, \
                         is0_25, is1_24, is1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_2 * hs0_18[k]
                  - f_2 * hs1_18[k]
                  + pa_y[k] * is0_24[k]
                  - pc_y[k] * is1_24[k];

        t_32[k] = f_3 * hs0_19[k]
                  - f_3 * hs1_19[k]
                  + pa_y[k] * is0_25[k]
                  - pc_y[k] * is1_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pa_z, pc_y, pc_z, hs0_20, hs1_20, is0_26, \
                         is0_27, is1_26, is1_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_4 * hs0_20[k]
                  - f_4 * hs1_20[k]
                  + pa_y[k] * is0_26[k]
                  - pc_y[k] * is1_26[k];

        t_34[k] = pa_y[k] * is0_27[k]
                  - pc_y[k] * is1_27[k];

        t_35[k] = f_0 * hs0_20[k]
                  - f_0 * hs1_20[k]
                  + pa_z[k] * is0_27[k]
                  - pc_z[k] * is1_27[k];
    }
}

}  // namespace simdnpot
