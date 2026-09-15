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


#include "SimdElectronRepulsionGeom10VrrRecKS.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_ks_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t is, const size_t ls,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_21 = buffer.data(is + 21);
    const auto *is_22 = buffer.data(is + 22);
    const auto *is_23 = buffer.data(is + 23);
    const auto *is_24 = buffer.data(is + 24);
    const auto *is_25 = buffer.data(is + 25);
    const auto *is_26 = buffer.data(is + 26);
    const auto *is_27 = buffer.data(is + 27);

    const auto *ls_0 = buffer.data(ls + 0);
    const auto *ls_1 = buffer.data(ls + 1);
    const auto *ls_2 = buffer.data(ls + 2);
    const auto *ls_3 = buffer.data(ls + 3);
    const auto *ls_4 = buffer.data(ls + 4);
    const auto *ls_5 = buffer.data(ls + 5);
    const auto *ls_6 = buffer.data(ls + 6);
    const auto *ls_7 = buffer.data(ls + 7);
    const auto *ls_8 = buffer.data(ls + 8);
    const auto *ls_9 = buffer.data(ls + 9);
    const auto *ls_10 = buffer.data(ls + 10);
    const auto *ls_11 = buffer.data(ls + 11);
    const auto *ls_12 = buffer.data(ls + 12);
    const auto *ls_13 = buffer.data(ls + 13);
    const auto *ls_14 = buffer.data(ls + 14);
    const auto *ls_15 = buffer.data(ls + 15);
    const auto *ls_16 = buffer.data(ls + 16);
    const auto *ls_17 = buffer.data(ls + 17);
    const auto *ls_18 = buffer.data(ls + 18);
    const auto *ls_19 = buffer.data(ls + 19);
    const auto *ls_20 = buffer.data(ls + 20);
    const auto *ls_21 = buffer.data(ls + 21);
    const auto *ls_22 = buffer.data(ls + 22);
    const auto *ls_23 = buffer.data(ls + 23);
    const auto *ls_24 = buffer.data(ls + 24);
    const auto *ls_25 = buffer.data(ls + 25);
    const auto *ls_26 = buffer.data(ls + 26);
    const auto *ls_27 = buffer.data(ls + 27);
    const auto *ls_28 = buffer.data(ls + 28);
    const auto *ls_29 = buffer.data(ls + 29);
    const auto *ls_30 = buffer.data(ls + 30);
    const auto *ls_31 = buffer.data(ls + 31);
    const auto *ls_32 = buffer.data(ls + 32);
    const auto *ls_33 = buffer.data(ls + 33);
    const auto *ls_34 = buffer.data(ls + 34);
    const auto *ls_35 = buffer.data(ls + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, is_0, is_1, is_2, is_3, is_4, ls_0, ls_1, \
                         ls_2, ls_3, ls_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -7.0 * is_0[k]
                 + f_0 * ls_0[k];

        t_1[k] = -6.0 * is_1[k]
                 + f_0 * ls_1[k];

        t_2[k] = -6.0 * is_2[k]
                 + f_0 * ls_2[k];

        t_3[k] = -5.0 * is_3[k]
                 + f_0 * ls_3[k];

        t_4[k] = -5.0 * is_4[k]
                 + f_0 * ls_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, is_5, is_6, is_7, is_8, is_9, ls_5, ls_6, \
                         ls_7, ls_8, ls_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -5.0 * is_5[k]
                 + f_0 * ls_5[k];

        t_6[k] = -4.0 * is_6[k]
                 + f_0 * ls_6[k];

        t_7[k] = -4.0 * is_7[k]
                 + f_0 * ls_7[k];

        t_8[k] = -4.0 * is_8[k]
                 + f_0 * ls_8[k];

        t_9[k] = -4.0 * is_9[k]
                 + f_0 * ls_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, is_10, is_11, is_12, is_13, is_14, \
                         ls_10, ls_11, ls_12, ls_13, ls_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -3.0 * is_10[k]
                  + f_0 * ls_10[k];

        t_11[k] = -3.0 * is_11[k]
                  + f_0 * ls_11[k];

        t_12[k] = -3.0 * is_12[k]
                  + f_0 * ls_12[k];

        t_13[k] = -3.0 * is_13[k]
                  + f_0 * ls_13[k];

        t_14[k] = -3.0 * is_14[k]
                  + f_0 * ls_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, is_15, is_16, is_17, is_18, is_19, \
                         ls_15, ls_16, ls_17, ls_18, ls_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -2.0 * is_15[k]
                  + f_0 * ls_15[k];

        t_16[k] = -2.0 * is_16[k]
                  + f_0 * ls_16[k];

        t_17[k] = -2.0 * is_17[k]
                  + f_0 * ls_17[k];

        t_18[k] = -2.0 * is_18[k]
                  + f_0 * ls_18[k];

        t_19[k] = -2.0 * is_19[k]
                  + f_0 * ls_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, is_20, is_21, is_22, is_23, is_24, \
                         ls_20, ls_21, ls_22, ls_23, ls_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -2.0 * is_20[k]
                  + f_0 * ls_20[k];

        t_21[k] = -is_21[k]
                  + f_0 * ls_21[k];

        t_22[k] = -is_22[k]
                  + f_0 * ls_22[k];

        t_23[k] = -is_23[k]
                  + f_0 * ls_23[k];

        t_24[k] = -is_24[k]
                  + f_0 * ls_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, t_30, is_25, is_26, is_27, ls_25, \
                         ls_26, ls_27, ls_28, ls_29, ls_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -is_25[k]
                  + f_0 * ls_25[k];

        t_26[k] = -is_26[k]
                  + f_0 * ls_26[k];

        t_27[k] = -is_27[k]
                  + f_0 * ls_27[k];

        t_28[k] = f_0 * ls_28[k];

        t_29[k] = f_0 * ls_29[k];

        t_30[k] = f_0 * ls_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ls_31, ls_32, ls_33, ls_34, \
                         ls_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * ls_31[k];

        t_32[k] = f_0 * ls_32[k];

        t_33[k] = f_0 * ls_33[k];

        t_34[k] = f_0 * ls_34[k];

        t_35[k] = f_0 * ls_35[k];
    }
}

auto
compute_prim_geom_10_ks_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t is, const size_t ls,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_21 = buffer.data(is + 21);
    const auto *is_22 = buffer.data(is + 22);
    const auto *is_23 = buffer.data(is + 23);
    const auto *is_24 = buffer.data(is + 24);
    const auto *is_25 = buffer.data(is + 25);
    const auto *is_26 = buffer.data(is + 26);
    const auto *is_27 = buffer.data(is + 27);

    const auto *ls_1 = buffer.data(ls + 1);
    const auto *ls_3 = buffer.data(ls + 3);
    const auto *ls_4 = buffer.data(ls + 4);
    const auto *ls_6 = buffer.data(ls + 6);
    const auto *ls_7 = buffer.data(ls + 7);
    const auto *ls_8 = buffer.data(ls + 8);
    const auto *ls_10 = buffer.data(ls + 10);
    const auto *ls_11 = buffer.data(ls + 11);
    const auto *ls_12 = buffer.data(ls + 12);
    const auto *ls_13 = buffer.data(ls + 13);
    const auto *ls_15 = buffer.data(ls + 15);
    const auto *ls_16 = buffer.data(ls + 16);
    const auto *ls_17 = buffer.data(ls + 17);
    const auto *ls_18 = buffer.data(ls + 18);
    const auto *ls_19 = buffer.data(ls + 19);
    const auto *ls_21 = buffer.data(ls + 21);
    const auto *ls_22 = buffer.data(ls + 22);
    const auto *ls_23 = buffer.data(ls + 23);
    const auto *ls_24 = buffer.data(ls + 24);
    const auto *ls_25 = buffer.data(ls + 25);
    const auto *ls_26 = buffer.data(ls + 26);
    const auto *ls_28 = buffer.data(ls + 28);
    const auto *ls_29 = buffer.data(ls + 29);
    const auto *ls_30 = buffer.data(ls + 30);
    const auto *ls_31 = buffer.data(ls + 31);
    const auto *ls_32 = buffer.data(ls + 32);
    const auto *ls_33 = buffer.data(ls + 33);
    const auto *ls_34 = buffer.data(ls + 34);
    const auto *ls_36 = buffer.data(ls + 36);
    const auto *ls_37 = buffer.data(ls + 37);
    const auto *ls_38 = buffer.data(ls + 38);
    const auto *ls_39 = buffer.data(ls + 39);
    const auto *ls_40 = buffer.data(ls + 40);
    const auto *ls_41 = buffer.data(ls + 41);
    const auto *ls_42 = buffer.data(ls + 42);
    const auto *ls_43 = buffer.data(ls + 43);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, is_0, is_1, is_2, ls_1, ls_3, ls_4, \
                         ls_6, ls_7, ls_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ls_1[k];

        t_1[k] = -is_0[k]
                 + f_0 * ls_3[k];

        t_2[k] = f_0 * ls_4[k];

        t_3[k] = -2.0 * is_1[k]
                 + f_0 * ls_6[k];

        t_4[k] = -is_2[k]
                 + f_0 * ls_7[k];

        t_5[k] = f_0 * ls_8[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, is_3, is_4, is_5, is_6, ls_10, ls_11, \
                         ls_12, ls_13, ls_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -3.0 * is_3[k]
                 + f_0 * ls_10[k];

        t_7[k] = -2.0 * is_4[k]
                 + f_0 * ls_11[k];

        t_8[k] = -is_5[k]
                 + f_0 * ls_12[k];

        t_9[k] = f_0 * ls_13[k];

        t_10[k] = -4.0 * is_6[k]
                  + f_0 * ls_15[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, is_7, is_8, is_9, is_10, ls_16, ls_17, \
                         ls_18, ls_19, ls_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = -3.0 * is_7[k]
                  + f_0 * ls_16[k];

        t_12[k] = -2.0 * is_8[k]
                  + f_0 * ls_17[k];

        t_13[k] = -is_9[k]
                  + f_0 * ls_18[k];

        t_14[k] = f_0 * ls_19[k];

        t_15[k] = -5.0 * is_10[k]
                  + f_0 * ls_21[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, is_11, is_12, is_13, is_14, ls_22, \
                         ls_23, ls_24, ls_25, ls_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -4.0 * is_11[k]
                  + f_0 * ls_22[k];

        t_17[k] = -3.0 * is_12[k]
                  + f_0 * ls_23[k];

        t_18[k] = -2.0 * is_13[k]
                  + f_0 * ls_24[k];

        t_19[k] = -is_14[k]
                  + f_0 * ls_25[k];

        t_20[k] = f_0 * ls_26[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, is_15, is_16, is_17, is_18, is_19, \
                         ls_28, ls_29, ls_30, ls_31, ls_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -6.0 * is_15[k]
                  + f_0 * ls_28[k];

        t_22[k] = -5.0 * is_16[k]
                  + f_0 * ls_29[k];

        t_23[k] = -4.0 * is_17[k]
                  + f_0 * ls_30[k];

        t_24[k] = -3.0 * is_18[k]
                  + f_0 * ls_31[k];

        t_25[k] = -2.0 * is_19[k]
                  + f_0 * ls_32[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, is_20, is_21, is_22, is_23, ls_33, \
                         ls_34, ls_36, ls_37, ls_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -is_20[k]
                  + f_0 * ls_33[k];

        t_27[k] = f_0 * ls_34[k];

        t_28[k] = -7.0 * is_21[k]
                  + f_0 * ls_36[k];

        t_29[k] = -6.0 * is_22[k]
                  + f_0 * ls_37[k];

        t_30[k] = -5.0 * is_23[k]
                  + f_0 * ls_38[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, is_24, is_25, is_26, is_27, ls_39, \
                         ls_40, ls_41, ls_42, ls_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -4.0 * is_24[k]
                  + f_0 * ls_39[k];

        t_32[k] = -3.0 * is_25[k]
                  + f_0 * ls_40[k];

        t_33[k] = -2.0 * is_26[k]
                  + f_0 * ls_41[k];

        t_34[k] = -is_27[k]
                  + f_0 * ls_42[k];

        t_35[k] = f_0 * ls_43[k];
    }
}

auto
compute_prim_geom_10_ks_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t is, const size_t ls,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_21 = buffer.data(is + 21);
    const auto *is_22 = buffer.data(is + 22);
    const auto *is_23 = buffer.data(is + 23);
    const auto *is_24 = buffer.data(is + 24);
    const auto *is_25 = buffer.data(is + 25);
    const auto *is_26 = buffer.data(is + 26);
    const auto *is_27 = buffer.data(is + 27);

    const auto *ls_2 = buffer.data(ls + 2);
    const auto *ls_4 = buffer.data(ls + 4);
    const auto *ls_5 = buffer.data(ls + 5);
    const auto *ls_7 = buffer.data(ls + 7);
    const auto *ls_8 = buffer.data(ls + 8);
    const auto *ls_9 = buffer.data(ls + 9);
    const auto *ls_11 = buffer.data(ls + 11);
    const auto *ls_12 = buffer.data(ls + 12);
    const auto *ls_13 = buffer.data(ls + 13);
    const auto *ls_14 = buffer.data(ls + 14);
    const auto *ls_16 = buffer.data(ls + 16);
    const auto *ls_17 = buffer.data(ls + 17);
    const auto *ls_18 = buffer.data(ls + 18);
    const auto *ls_19 = buffer.data(ls + 19);
    const auto *ls_20 = buffer.data(ls + 20);
    const auto *ls_22 = buffer.data(ls + 22);
    const auto *ls_23 = buffer.data(ls + 23);
    const auto *ls_24 = buffer.data(ls + 24);
    const auto *ls_25 = buffer.data(ls + 25);
    const auto *ls_26 = buffer.data(ls + 26);
    const auto *ls_27 = buffer.data(ls + 27);
    const auto *ls_29 = buffer.data(ls + 29);
    const auto *ls_30 = buffer.data(ls + 30);
    const auto *ls_31 = buffer.data(ls + 31);
    const auto *ls_32 = buffer.data(ls + 32);
    const auto *ls_33 = buffer.data(ls + 33);
    const auto *ls_34 = buffer.data(ls + 34);
    const auto *ls_35 = buffer.data(ls + 35);
    const auto *ls_37 = buffer.data(ls + 37);
    const auto *ls_38 = buffer.data(ls + 38);
    const auto *ls_39 = buffer.data(ls + 39);
    const auto *ls_40 = buffer.data(ls + 40);
    const auto *ls_41 = buffer.data(ls + 41);
    const auto *ls_42 = buffer.data(ls + 42);
    const auto *ls_43 = buffer.data(ls + 43);
    const auto *ls_44 = buffer.data(ls + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, is_0, is_1, is_2, ls_2, ls_4, ls_5, \
                         ls_7, ls_8, ls_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ls_2[k];

        t_1[k] = f_0 * ls_4[k];

        t_2[k] = -is_0[k]
                 + f_0 * ls_5[k];

        t_3[k] = f_0 * ls_7[k];

        t_4[k] = -is_1[k]
                 + f_0 * ls_8[k];

        t_5[k] = -2.0 * is_2[k]
                 + f_0 * ls_9[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, is_3, is_4, is_5, is_6, ls_11, ls_12, \
                         ls_13, ls_14, ls_16, ls_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * ls_11[k];

        t_7[k] = -is_3[k]
                 + f_0 * ls_12[k];

        t_8[k] = -2.0 * is_4[k]
                 + f_0 * ls_13[k];

        t_9[k] = -3.0 * is_5[k]
                 + f_0 * ls_14[k];

        t_10[k] = f_0 * ls_16[k];

        t_11[k] = -is_6[k]
                  + f_0 * ls_17[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, is_7, is_8, is_9, is_10, ls_18, ls_19, \
                         ls_20, ls_22, ls_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -2.0 * is_7[k]
                  + f_0 * ls_18[k];

        t_13[k] = -3.0 * is_8[k]
                  + f_0 * ls_19[k];

        t_14[k] = -4.0 * is_9[k]
                  + f_0 * ls_20[k];

        t_15[k] = f_0 * ls_22[k];

        t_16[k] = -is_10[k]
                  + f_0 * ls_23[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, is_11, is_12, is_13, is_14, ls_24, \
                         ls_25, ls_26, ls_27, ls_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -2.0 * is_11[k]
                  + f_0 * ls_24[k];

        t_18[k] = -3.0 * is_12[k]
                  + f_0 * ls_25[k];

        t_19[k] = -4.0 * is_13[k]
                  + f_0 * ls_26[k];

        t_20[k] = -5.0 * is_14[k]
                  + f_0 * ls_27[k];

        t_21[k] = f_0 * ls_29[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, is_15, is_16, is_17, is_18, is_19, \
                         ls_30, ls_31, ls_32, ls_33, ls_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -is_15[k]
                  + f_0 * ls_30[k];

        t_23[k] = -2.0 * is_16[k]
                  + f_0 * ls_31[k];

        t_24[k] = -3.0 * is_17[k]
                  + f_0 * ls_32[k];

        t_25[k] = -4.0 * is_18[k]
                  + f_0 * ls_33[k];

        t_26[k] = -5.0 * is_19[k]
                  + f_0 * ls_34[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, is_20, is_21, is_22, is_23, ls_35, \
                         ls_37, ls_38, ls_39, ls_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -6.0 * is_20[k]
                  + f_0 * ls_35[k];

        t_28[k] = f_0 * ls_37[k];

        t_29[k] = -is_21[k]
                  + f_0 * ls_38[k];

        t_30[k] = -2.0 * is_22[k]
                  + f_0 * ls_39[k];

        t_31[k] = -3.0 * is_23[k]
                  + f_0 * ls_40[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, is_24, is_25, is_26, is_27, ls_41, ls_42, \
                         ls_43, ls_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -4.0 * is_24[k]
                  + f_0 * ls_41[k];

        t_33[k] = -5.0 * is_25[k]
                  + f_0 * ls_42[k];

        t_34[k] = -6.0 * is_26[k]
                  + f_0 * ls_43[k];

        t_35[k] = -7.0 * is_27[k]
                  + f_0 * ls_44[k];
    }
}

}  // namespace simdt2ceri
