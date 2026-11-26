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


#include "SimdNuclearPotentialVrrRecSM.hpp"

#include "SimdAlign.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_prim_sm_nuclear_potential_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                    const size_t pc, const size_t sk0, const size_t sl0,
                                    const size_t sk1, const size_t sl1, const size_t ncols,
                                    const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 3.0 / p;
    const auto f_2 = 2.5 / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = 1.5 / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = 0.5 / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sk0_0 = buffer.data(sk0 + 0);
    const auto *sk0_3 = buffer.data(sk0 + 3);
    const auto *sk0_5 = buffer.data(sk0 + 5);
    const auto *sk0_6 = buffer.data(sk0 + 6);
    const auto *sk0_9 = buffer.data(sk0 + 9);
    const auto *sk0_10 = buffer.data(sk0 + 10);
    const auto *sk0_12 = buffer.data(sk0 + 12);
    const auto *sk0_14 = buffer.data(sk0 + 14);
    const auto *sk0_15 = buffer.data(sk0 + 15);
    const auto *sk0_17 = buffer.data(sk0 + 17);
    const auto *sk0_18 = buffer.data(sk0 + 18);
    const auto *sk0_20 = buffer.data(sk0 + 20);
    const auto *sk0_21 = buffer.data(sk0 + 21);
    const auto *sk0_23 = buffer.data(sk0 + 23);
    const auto *sk0_24 = buffer.data(sk0 + 24);
    const auto *sk0_25 = buffer.data(sk0 + 25);
    const auto *sk0_27 = buffer.data(sk0 + 27);
    const auto *sk0_28 = buffer.data(sk0 + 28);
    const auto *sk0_30 = buffer.data(sk0 + 30);
    const auto *sk0_31 = buffer.data(sk0 + 31);
    const auto *sk0_32 = buffer.data(sk0 + 32);
    const auto *sk0_33 = buffer.data(sk0 + 33);
    const auto *sk0_34 = buffer.data(sk0 + 34);
    const auto *sk0_35 = buffer.data(sk0 + 35);

    const auto *sl0_0 = buffer.data(sl0 + 0);
    const auto *sl0_2 = buffer.data(sl0 + 2);
    const auto *sl0_3 = buffer.data(sl0 + 3);
    const auto *sl0_5 = buffer.data(sl0 + 5);
    const auto *sl0_6 = buffer.data(sl0 + 6);
    const auto *sl0_9 = buffer.data(sl0 + 9);
    const auto *sl0_10 = buffer.data(sl0 + 10);
    const auto *sl0_12 = buffer.data(sl0 + 12);
    const auto *sl0_14 = buffer.data(sl0 + 14);
    const auto *sl0_15 = buffer.data(sl0 + 15);
    const auto *sl0_17 = buffer.data(sl0 + 17);
    const auto *sl0_18 = buffer.data(sl0 + 18);
    const auto *sl0_20 = buffer.data(sl0 + 20);
    const auto *sl0_21 = buffer.data(sl0 + 21);
    const auto *sl0_23 = buffer.data(sl0 + 23);
    const auto *sl0_24 = buffer.data(sl0 + 24);
    const auto *sl0_25 = buffer.data(sl0 + 25);
    const auto *sl0_27 = buffer.data(sl0 + 27);
    const auto *sl0_28 = buffer.data(sl0 + 28);
    const auto *sl0_30 = buffer.data(sl0 + 30);
    const auto *sl0_31 = buffer.data(sl0 + 31);
    const auto *sl0_32 = buffer.data(sl0 + 32);
    const auto *sl0_33 = buffer.data(sl0 + 33);
    const auto *sl0_35 = buffer.data(sl0 + 35);
    const auto *sl0_36 = buffer.data(sl0 + 36);
    const auto *sl0_37 = buffer.data(sl0 + 37);
    const auto *sl0_38 = buffer.data(sl0 + 38);
    const auto *sl0_39 = buffer.data(sl0 + 39);
    const auto *sl0_40 = buffer.data(sl0 + 40);
    const auto *sl0_41 = buffer.data(sl0 + 41);
    const auto *sl0_42 = buffer.data(sl0 + 42);
    const auto *sl0_43 = buffer.data(sl0 + 43);
    const auto *sl0_44 = buffer.data(sl0 + 44);

    const auto *sk1_0 = buffer.data(sk1 + 0);
    const auto *sk1_3 = buffer.data(sk1 + 3);
    const auto *sk1_5 = buffer.data(sk1 + 5);
    const auto *sk1_6 = buffer.data(sk1 + 6);
    const auto *sk1_9 = buffer.data(sk1 + 9);
    const auto *sk1_10 = buffer.data(sk1 + 10);
    const auto *sk1_12 = buffer.data(sk1 + 12);
    const auto *sk1_14 = buffer.data(sk1 + 14);
    const auto *sk1_15 = buffer.data(sk1 + 15);
    const auto *sk1_17 = buffer.data(sk1 + 17);
    const auto *sk1_18 = buffer.data(sk1 + 18);
    const auto *sk1_20 = buffer.data(sk1 + 20);
    const auto *sk1_21 = buffer.data(sk1 + 21);
    const auto *sk1_23 = buffer.data(sk1 + 23);
    const auto *sk1_24 = buffer.data(sk1 + 24);
    const auto *sk1_25 = buffer.data(sk1 + 25);
    const auto *sk1_27 = buffer.data(sk1 + 27);
    const auto *sk1_28 = buffer.data(sk1 + 28);
    const auto *sk1_30 = buffer.data(sk1 + 30);
    const auto *sk1_31 = buffer.data(sk1 + 31);
    const auto *sk1_32 = buffer.data(sk1 + 32);
    const auto *sk1_33 = buffer.data(sk1 + 33);
    const auto *sk1_34 = buffer.data(sk1 + 34);
    const auto *sk1_35 = buffer.data(sk1 + 35);

    const auto *sl1_0 = buffer.data(sl1 + 0);
    const auto *sl1_2 = buffer.data(sl1 + 2);
    const auto *sl1_3 = buffer.data(sl1 + 3);
    const auto *sl1_5 = buffer.data(sl1 + 5);
    const auto *sl1_6 = buffer.data(sl1 + 6);
    const auto *sl1_9 = buffer.data(sl1 + 9);
    const auto *sl1_10 = buffer.data(sl1 + 10);
    const auto *sl1_12 = buffer.data(sl1 + 12);
    const auto *sl1_14 = buffer.data(sl1 + 14);
    const auto *sl1_15 = buffer.data(sl1 + 15);
    const auto *sl1_17 = buffer.data(sl1 + 17);
    const auto *sl1_18 = buffer.data(sl1 + 18);
    const auto *sl1_20 = buffer.data(sl1 + 20);
    const auto *sl1_21 = buffer.data(sl1 + 21);
    const auto *sl1_23 = buffer.data(sl1 + 23);
    const auto *sl1_24 = buffer.data(sl1 + 24);
    const auto *sl1_25 = buffer.data(sl1 + 25);
    const auto *sl1_27 = buffer.data(sl1 + 27);
    const auto *sl1_28 = buffer.data(sl1 + 28);
    const auto *sl1_30 = buffer.data(sl1 + 30);
    const auto *sl1_31 = buffer.data(sl1 + 31);
    const auto *sl1_32 = buffer.data(sl1 + 32);
    const auto *sl1_33 = buffer.data(sl1 + 33);
    const auto *sl1_35 = buffer.data(sl1 + 35);
    const auto *sl1_36 = buffer.data(sl1 + 36);
    const auto *sl1_37 = buffer.data(sl1 + 37);
    const auto *sl1_38 = buffer.data(sl1 + 38);
    const auto *sl1_39 = buffer.data(sl1 + 39);
    const auto *sl1_40 = buffer.data(sl1 + 40);
    const auto *sl1_41 = buffer.data(sl1 + 41);
    const auto *sl1_42 = buffer.data(sl1 + 42);
    const auto *sl1_43 = buffer.data(sl1 + 43);
    const auto *sl1_44 = buffer.data(sl1 + 44);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, sk0_0, sl0_0, \
                         sk1_0, sl1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sk0_0[k]
                 + pb_x[k] * sl0_0[k]
                 - f_0 * sk1_0[k]
                 - pc_x[k] * sl1_0[k];

        t_1[k] = pb_y[k] * sl0_0[k]
                 - pc_y[k] * sl1_0[k];

        t_2[k] = pb_z[k] * sl0_0[k]
                 - pc_z[k] * sl1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pb_x, pb_y, pc_x, pc_y, sk0_3, sl0_2, sl0_3, sk1_3, sl1_2, \
                         sl1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_1 * sk0_3[k]
                 + pb_x[k] * sl0_3[k]
                 - f_1 * sk1_3[k]
                 - pc_x[k] * sl1_3[k];

        t_4[k] = pb_y[k] * sl0_2[k]
                 - pc_y[k] * sl1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pb_x, pc_x, sk0_5, sk0_6, sl0_5, sl0_6, sk1_5, sk1_6, \
                         sl1_5, sl1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sk0_5[k]
                 + pb_x[k] * sl0_5[k]
                 - f_1 * sk1_5[k]
                 - pc_x[k] * sl1_5[k];

        t_6[k] = f_2 * sk0_6[k]
                 + pb_x[k] * sl0_6[k]
                 - f_2 * sk1_6[k]
                 - pc_x[k] * sl1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pb_y, pb_z, pc_y, pc_z, sl0_3, sl0_5, sl1_3, \
                         sl1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_z[k] * sl0_3[k]
                 - pc_z[k] * sl1_3[k];

        t_8[k] = pb_y[k] * sl0_5[k]
                 - pc_y[k] * sl1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, pb_x, pc_x, sk0_9, sk0_10, sl0_9, sl0_10, sk1_9, sk1_10, \
                         sl1_9, sl1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * sk0_9[k]
                 + pb_x[k] * sl0_9[k]
                 - f_2 * sk1_9[k]
                 - pc_x[k] * sl1_9[k];

        t_10[k] = f_3 * sk0_10[k]
                  + pb_x[k] * sl0_10[k]
                  - f_3 * sk1_10[k]
                  - pc_x[k] * sl1_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pb_x, pb_z, pc_x, pc_z, sk0_12, sl0_6, sl0_12, sk1_12, \
                         sl1_6, sl1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * sl0_6[k]
                  - pc_z[k] * sl1_6[k];

        t_12[k] = f_3 * sk0_12[k]
                  + pb_x[k] * sl0_12[k]
                  - f_3 * sk1_12[k]
                  - pc_x[k] * sl1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pb_x, pb_y, pc_x, pc_y, sk0_14, sl0_9, sl0_14, sk1_14, \
                         sl1_9, sl1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_y[k] * sl0_9[k]
                  - pc_y[k] * sl1_9[k];

        t_14[k] = f_3 * sk0_14[k]
                  + pb_x[k] * sl0_14[k]
                  - f_3 * sk1_14[k]
                  - pc_x[k] * sl1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pb_x, pb_z, pc_x, pc_z, sk0_15, sl0_10, sl0_15, sk1_15, \
                         sl1_10, sl1_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * sk0_15[k]
                  + pb_x[k] * sl0_15[k]
                  - f_4 * sk1_15[k]
                  - pc_x[k] * sl1_15[k];

        t_16[k] = pb_z[k] * sl0_10[k]
                  - pc_z[k] * sl1_10[k];
    }

#pragma omp simd aligned(t_17, t_18, pb_x, pc_x, sk0_17, sk0_18, sl0_17, sl0_18, sk1_17, \
                         sk1_18, sl1_17, sl1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_4 * sk0_17[k]
                  + pb_x[k] * sl0_17[k]
                  - f_4 * sk1_17[k]
                  - pc_x[k] * sl1_17[k];

        t_18[k] = f_4 * sk0_18[k]
                  + pb_x[k] * sl0_18[k]
                  - f_4 * sk1_18[k]
                  - pc_x[k] * sl1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, pb_x, pb_y, pc_x, pc_y, sk0_20, sl0_14, sl0_20, sk1_20, \
                         sl1_14, sl1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * sl0_14[k]
                  - pc_y[k] * sl1_14[k];

        t_20[k] = f_4 * sk0_20[k]
                  + pb_x[k] * sl0_20[k]
                  - f_4 * sk1_20[k]
                  - pc_x[k] * sl1_20[k];
    }

#pragma omp simd aligned(t_21, t_22, pb_x, pb_z, pc_x, pc_z, sk0_21, sl0_15, sl0_21, sk1_21, \
                         sl1_15, sl1_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * sk0_21[k]
                  + pb_x[k] * sl0_21[k]
                  - f_5 * sk1_21[k]
                  - pc_x[k] * sl1_21[k];

        t_22[k] = pb_z[k] * sl0_15[k]
                  - pc_z[k] * sl1_15[k];
    }

#pragma omp simd aligned(t_23, t_24, pb_x, pc_x, sk0_23, sk0_24, sl0_23, sl0_24, sk1_23, \
                         sk1_24, sl1_23, sl1_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * sk0_23[k]
                  + pb_x[k] * sl0_23[k]
                  - f_5 * sk1_23[k]
                  - pc_x[k] * sl1_23[k];

        t_24[k] = f_5 * sk0_24[k]
                  + pb_x[k] * sl0_24[k]
                  - f_5 * sk1_24[k]
                  - pc_x[k] * sl1_24[k];
    }

#pragma omp simd aligned(t_25, t_26, pb_x, pb_y, pc_x, pc_y, sk0_25, sl0_20, sl0_25, sk1_25, \
                         sl1_20, sl1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_5 * sk0_25[k]
                  + pb_x[k] * sl0_25[k]
                  - f_5 * sk1_25[k]
                  - pc_x[k] * sl1_25[k];

        t_26[k] = pb_y[k] * sl0_20[k]
                  - pc_y[k] * sl1_20[k];
    }

#pragma omp simd aligned(t_27, t_28, pb_x, pc_x, sk0_27, sk0_28, sl0_27, sl0_28, sk1_27, \
                         sk1_28, sl1_27, sl1_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * sk0_27[k]
                  + pb_x[k] * sl0_27[k]
                  - f_5 * sk1_27[k]
                  - pc_x[k] * sl1_27[k];

        t_28[k] = f_6 * sk0_28[k]
                  + pb_x[k] * sl0_28[k]
                  - f_6 * sk1_28[k]
                  - pc_x[k] * sl1_28[k];
    }

#pragma omp simd aligned(t_29, t_30, pb_x, pb_z, pc_x, pc_z, sk0_30, sl0_21, sl0_30, sk1_30, \
                         sl1_21, sl1_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_z[k] * sl0_21[k]
                  - pc_z[k] * sl1_21[k];

        t_30[k] = f_6 * sk0_30[k]
                  + pb_x[k] * sl0_30[k]
                  - f_6 * sk1_30[k]
                  - pc_x[k] * sl1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, pb_x, pc_x, sk0_31, sk0_32, sl0_31, sl0_32, sk1_31, \
                         sk1_32, sl1_31, sl1_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_6 * sk0_31[k]
                  + pb_x[k] * sl0_31[k]
                  - f_6 * sk1_31[k]
                  - pc_x[k] * sl1_31[k];

        t_32[k] = f_6 * sk0_32[k]
                  + pb_x[k] * sl0_32[k]
                  - f_6 * sk1_32[k]
                  - pc_x[k] * sl1_32[k];
    }

#pragma omp simd aligned(t_33, t_34, pb_x, pb_y, pc_x, pc_y, sk0_33, sl0_27, sl0_33, sk1_33, \
                         sl1_27, sl1_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_6 * sk0_33[k]
                  + pb_x[k] * sl0_33[k]
                  - f_6 * sk1_33[k]
                  - pc_x[k] * sl1_33[k];

        t_34[k] = pb_y[k] * sl0_27[k]
                  - pc_y[k] * sl1_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, pc_x, sk0_35, sl0_35, sl0_36, sl0_37, \
                         sl0_38, sk1_35, sl1_35, sl1_36, sl1_37, \
                         sl1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_6 * sk0_35[k]
                  + pb_x[k] * sl0_35[k]
                  - f_6 * sk1_35[k]
                  - pc_x[k] * sl1_35[k];

        t_36[k] = pb_x[k] * sl0_36[k]
                  - pc_x[k] * sl1_36[k];

        t_37[k] = pb_x[k] * sl0_37[k]
                  - pc_x[k] * sl1_37[k];

        t_38[k] = pb_x[k] * sl0_38[k]
                  - pc_x[k] * sl1_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pb_x, pc_x, sl0_39, sl0_40, sl0_41, sl0_42, \
                         sl1_39, sl1_40, sl1_41, sl1_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pb_x[k] * sl0_39[k]
                  - pc_x[k] * sl1_39[k];

        t_40[k] = pb_x[k] * sl0_40[k]
                  - pc_x[k] * sl1_40[k];

        t_41[k] = pb_x[k] * sl0_41[k]
                  - pc_x[k] * sl1_41[k];

        t_42[k] = pb_x[k] * sl0_42[k]
                  - pc_x[k] * sl1_42[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pb_x, pb_y, pc_x, pc_y, sk0_28, sl0_36, sl0_43, \
                         sl0_44, sk1_28, sl1_36, sl1_43, sl1_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_x[k] * sl0_43[k]
                  - pc_x[k] * sl1_43[k];

        t_44[k] = pb_x[k] * sl0_44[k]
                  - pc_x[k] * sl1_44[k];

        t_45[k] = f_0 * sk0_28[k]
                  + pb_y[k] * sl0_36[k]
                  - f_0 * sk1_28[k]
                  - pc_y[k] * sl1_36[k];
    }

#pragma omp simd aligned(t_46, t_47, pb_y, pb_z, pc_y, pc_z, sk0_30, sl0_36, sl0_38, sk1_30, \
                         sl1_36, sl1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_z[k] * sl0_36[k]
                  - pc_z[k] * sl1_36[k];

        t_47[k] = f_1 * sk0_30[k]
                  + pb_y[k] * sl0_38[k]
                  - f_1 * sk1_30[k]
                  - pc_y[k] * sl1_38[k];
    }

#pragma omp simd aligned(t_48, t_49, pb_y, pc_y, sk0_31, sk0_32, sl0_39, sl0_40, sk1_31, \
                         sk1_32, sl1_39, sl1_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_2 * sk0_31[k]
                  + pb_y[k] * sl0_39[k]
                  - f_2 * sk1_31[k]
                  - pc_y[k] * sl1_39[k];

        t_49[k] = f_3 * sk0_32[k]
                  + pb_y[k] * sl0_40[k]
                  - f_3 * sk1_32[k]
                  - pc_y[k] * sl1_40[k];
    }

#pragma omp simd aligned(t_50, t_51, pb_y, pc_y, sk0_33, sk0_34, sl0_41, sl0_42, sk1_33, \
                         sk1_34, sl1_41, sl1_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_4 * sk0_33[k]
                  + pb_y[k] * sl0_41[k]
                  - f_4 * sk1_33[k]
                  - pc_y[k] * sl1_41[k];

        t_51[k] = f_5 * sk0_34[k]
                  + pb_y[k] * sl0_42[k]
                  - f_5 * sk1_34[k]
                  - pc_y[k] * sl1_42[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_y, pb_z, pc_y, pc_z, sk0_35, sl0_43, sl0_44, \
                         sk1_35, sl1_43, sl1_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_6 * sk0_35[k]
                  + pb_y[k] * sl0_43[k]
                  - f_6 * sk1_35[k]
                  - pc_y[k] * sl1_43[k];

        t_53[k] = pb_y[k] * sl0_44[k]
                  - pc_y[k] * sl1_44[k];

        t_54[k] = f_0 * sk0_35[k]
                  + pb_z[k] * sl0_44[k]
                  - f_0 * sk1_35[k]
                  - pc_z[k] * sl1_44[k];
    }
}

}  // namespace simdnpot
