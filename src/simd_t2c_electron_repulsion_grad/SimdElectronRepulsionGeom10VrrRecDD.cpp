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


#include "SimdElectronRepulsionGeom10VrrRecDD.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_dd_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t pd, const size_t fd,
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

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_7 = buffer.data(pd + 7);
    const auto *pd_8 = buffer.data(pd + 8);
    const auto *pd_9 = buffer.data(pd + 9);
    const auto *pd_10 = buffer.data(pd + 10);
    const auto *pd_11 = buffer.data(pd + 11);
    const auto *pd_12 = buffer.data(pd + 12);
    const auto *pd_13 = buffer.data(pd + 13);
    const auto *pd_14 = buffer.data(pd + 14);
    const auto *pd_15 = buffer.data(pd + 15);
    const auto *pd_16 = buffer.data(pd + 16);
    const auto *pd_17 = buffer.data(pd + 17);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_22 = buffer.data(fd + 22);
    const auto *fd_23 = buffer.data(fd + 23);
    const auto *fd_24 = buffer.data(fd + 24);
    const auto *fd_25 = buffer.data(fd + 25);
    const auto *fd_26 = buffer.data(fd + 26);
    const auto *fd_27 = buffer.data(fd + 27);
    const auto *fd_28 = buffer.data(fd + 28);
    const auto *fd_29 = buffer.data(fd + 29);
    const auto *fd_30 = buffer.data(fd + 30);
    const auto *fd_31 = buffer.data(fd + 31);
    const auto *fd_32 = buffer.data(fd + 32);
    const auto *fd_33 = buffer.data(fd + 33);
    const auto *fd_34 = buffer.data(fd + 34);
    const auto *fd_35 = buffer.data(fd + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pd_0, pd_1, pd_2, pd_3, pd_4, fd_0, fd_1, \
                         fd_2, fd_3, fd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -2.0 * pd_0[k]
                 + f_0 * fd_0[k];

        t_1[k] = -2.0 * pd_1[k]
                 + f_0 * fd_1[k];

        t_2[k] = -2.0 * pd_2[k]
                 + f_0 * fd_2[k];

        t_3[k] = -2.0 * pd_3[k]
                 + f_0 * fd_3[k];

        t_4[k] = -2.0 * pd_4[k]
                 + f_0 * fd_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pd_5, pd_6, pd_7, pd_8, pd_9, fd_5, fd_6, \
                         fd_7, fd_8, fd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -2.0 * pd_5[k]
                 + f_0 * fd_5[k];

        t_6[k] = -pd_6[k]
                 + f_0 * fd_6[k];

        t_7[k] = -pd_7[k]
                 + f_0 * fd_7[k];

        t_8[k] = -pd_8[k]
                 + f_0 * fd_8[k];

        t_9[k] = -pd_9[k]
                 + f_0 * fd_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pd_10, pd_11, pd_12, pd_13, pd_14, \
                         fd_10, fd_11, fd_12, fd_13, fd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -pd_10[k]
                  + f_0 * fd_10[k];

        t_11[k] = -pd_11[k]
                  + f_0 * fd_11[k];

        t_12[k] = -pd_12[k]
                  + f_0 * fd_12[k];

        t_13[k] = -pd_13[k]
                  + f_0 * fd_13[k];

        t_14[k] = -pd_14[k]
                  + f_0 * fd_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, t_20, pd_15, pd_16, pd_17, fd_15, \
                         fd_16, fd_17, fd_18, fd_19, fd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -pd_15[k]
                  + f_0 * fd_15[k];

        t_16[k] = -pd_16[k]
                  + f_0 * fd_16[k];

        t_17[k] = -pd_17[k]
                  + f_0 * fd_17[k];

        t_18[k] = f_0 * fd_18[k];

        t_19[k] = f_0 * fd_19[k];

        t_20[k] = f_0 * fd_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, t_26, t_27, t_28, fd_21, fd_22, fd_23, \
                         fd_24, fd_25, fd_26, fd_27, fd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * fd_21[k];

        t_22[k] = f_0 * fd_22[k];

        t_23[k] = f_0 * fd_23[k];

        t_24[k] = f_0 * fd_24[k];

        t_25[k] = f_0 * fd_25[k];

        t_26[k] = f_0 * fd_26[k];

        t_27[k] = f_0 * fd_27[k];

        t_28[k] = f_0 * fd_28[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, t_35, fd_29, fd_30, fd_31, fd_32, \
                         fd_33, fd_34, fd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * fd_29[k];

        t_30[k] = f_0 * fd_30[k];

        t_31[k] = f_0 * fd_31[k];

        t_32[k] = f_0 * fd_32[k];

        t_33[k] = f_0 * fd_33[k];

        t_34[k] = f_0 * fd_34[k];

        t_35[k] = f_0 * fd_35[k];
    }
}

auto
compute_prim_geom_10_dd_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t pd, const size_t fd,
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

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_7 = buffer.data(pd + 7);
    const auto *pd_8 = buffer.data(pd + 8);
    const auto *pd_9 = buffer.data(pd + 9);
    const auto *pd_10 = buffer.data(pd + 10);
    const auto *pd_11 = buffer.data(pd + 11);
    const auto *pd_12 = buffer.data(pd + 12);
    const auto *pd_13 = buffer.data(pd + 13);
    const auto *pd_14 = buffer.data(pd + 14);
    const auto *pd_15 = buffer.data(pd + 15);
    const auto *pd_16 = buffer.data(pd + 16);
    const auto *pd_17 = buffer.data(pd + 17);

    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_22 = buffer.data(fd + 22);
    const auto *fd_23 = buffer.data(fd + 23);
    const auto *fd_24 = buffer.data(fd + 24);
    const auto *fd_25 = buffer.data(fd + 25);
    const auto *fd_26 = buffer.data(fd + 26);
    const auto *fd_27 = buffer.data(fd + 27);
    const auto *fd_28 = buffer.data(fd + 28);
    const auto *fd_29 = buffer.data(fd + 29);
    const auto *fd_36 = buffer.data(fd + 36);
    const auto *fd_37 = buffer.data(fd + 37);
    const auto *fd_38 = buffer.data(fd + 38);
    const auto *fd_39 = buffer.data(fd + 39);
    const auto *fd_40 = buffer.data(fd + 40);
    const auto *fd_41 = buffer.data(fd + 41);
    const auto *fd_42 = buffer.data(fd + 42);
    const auto *fd_43 = buffer.data(fd + 43);
    const auto *fd_44 = buffer.data(fd + 44);
    const auto *fd_45 = buffer.data(fd + 45);
    const auto *fd_46 = buffer.data(fd + 46);
    const auto *fd_47 = buffer.data(fd + 47);
    const auto *fd_48 = buffer.data(fd + 48);
    const auto *fd_49 = buffer.data(fd + 49);
    const auto *fd_50 = buffer.data(fd + 50);
    const auto *fd_51 = buffer.data(fd + 51);
    const auto *fd_52 = buffer.data(fd + 52);
    const auto *fd_53 = buffer.data(fd + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, pd_0, fd_6, fd_7, fd_8, fd_9, \
                         fd_10, fd_11, fd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_6[k];

        t_1[k] = f_0 * fd_7[k];

        t_2[k] = f_0 * fd_8[k];

        t_3[k] = f_0 * fd_9[k];

        t_4[k] = f_0 * fd_10[k];

        t_5[k] = f_0 * fd_11[k];

        t_6[k] = -pd_0[k]
                 + f_0 * fd_18[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, pd_1, pd_2, pd_3, pd_4, pd_5, fd_19, \
                         fd_20, fd_21, fd_22, fd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -pd_1[k]
                 + f_0 * fd_19[k];

        t_8[k] = -pd_2[k]
                 + f_0 * fd_20[k];

        t_9[k] = -pd_3[k]
                 + f_0 * fd_21[k];

        t_10[k] = -pd_4[k]
                  + f_0 * fd_22[k];

        t_11[k] = -pd_5[k]
                  + f_0 * fd_23[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, t_18, pd_6, fd_24, fd_25, fd_26, \
                         fd_27, fd_28, fd_29, fd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * fd_24[k];

        t_13[k] = f_0 * fd_25[k];

        t_14[k] = f_0 * fd_26[k];

        t_15[k] = f_0 * fd_27[k];

        t_16[k] = f_0 * fd_28[k];

        t_17[k] = f_0 * fd_29[k];

        t_18[k] = -2.0 * pd_6[k]
                  + f_0 * fd_36[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pd_7, pd_8, pd_9, pd_10, pd_11, fd_37, \
                         fd_38, fd_39, fd_40, fd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -2.0 * pd_7[k]
                  + f_0 * fd_37[k];

        t_20[k] = -2.0 * pd_8[k]
                  + f_0 * fd_38[k];

        t_21[k] = -2.0 * pd_9[k]
                  + f_0 * fd_39[k];

        t_22[k] = -2.0 * pd_10[k]
                  + f_0 * fd_40[k];

        t_23[k] = -2.0 * pd_11[k]
                  + f_0 * fd_41[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pd_12, pd_13, pd_14, pd_15, pd_16, \
                         fd_42, fd_43, fd_44, fd_45, fd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -pd_12[k]
                  + f_0 * fd_42[k];

        t_25[k] = -pd_13[k]
                  + f_0 * fd_43[k];

        t_26[k] = -pd_14[k]
                  + f_0 * fd_44[k];

        t_27[k] = -pd_15[k]
                  + f_0 * fd_45[k];

        t_28[k] = -pd_16[k]
                  + f_0 * fd_46[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, t_35, pd_17, fd_47, fd_48, fd_49, \
                         fd_50, fd_51, fd_52, fd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -pd_17[k]
                  + f_0 * fd_47[k];

        t_30[k] = f_0 * fd_48[k];

        t_31[k] = f_0 * fd_49[k];

        t_32[k] = f_0 * fd_50[k];

        t_33[k] = f_0 * fd_51[k];

        t_34[k] = f_0 * fd_52[k];

        t_35[k] = f_0 * fd_53[k];
    }
}

auto
compute_prim_geom_10_dd_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t pd, const size_t fd,
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

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_7 = buffer.data(pd + 7);
    const auto *pd_8 = buffer.data(pd + 8);
    const auto *pd_9 = buffer.data(pd + 9);
    const auto *pd_10 = buffer.data(pd + 10);
    const auto *pd_11 = buffer.data(pd + 11);
    const auto *pd_12 = buffer.data(pd + 12);
    const auto *pd_13 = buffer.data(pd + 13);
    const auto *pd_14 = buffer.data(pd + 14);
    const auto *pd_15 = buffer.data(pd + 15);
    const auto *pd_16 = buffer.data(pd + 16);
    const auto *pd_17 = buffer.data(pd + 17);

    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_24 = buffer.data(fd + 24);
    const auto *fd_25 = buffer.data(fd + 25);
    const auto *fd_26 = buffer.data(fd + 26);
    const auto *fd_27 = buffer.data(fd + 27);
    const auto *fd_28 = buffer.data(fd + 28);
    const auto *fd_29 = buffer.data(fd + 29);
    const auto *fd_30 = buffer.data(fd + 30);
    const auto *fd_31 = buffer.data(fd + 31);
    const auto *fd_32 = buffer.data(fd + 32);
    const auto *fd_33 = buffer.data(fd + 33);
    const auto *fd_34 = buffer.data(fd + 34);
    const auto *fd_35 = buffer.data(fd + 35);
    const auto *fd_42 = buffer.data(fd + 42);
    const auto *fd_43 = buffer.data(fd + 43);
    const auto *fd_44 = buffer.data(fd + 44);
    const auto *fd_45 = buffer.data(fd + 45);
    const auto *fd_46 = buffer.data(fd + 46);
    const auto *fd_47 = buffer.data(fd + 47);
    const auto *fd_48 = buffer.data(fd + 48);
    const auto *fd_49 = buffer.data(fd + 49);
    const auto *fd_50 = buffer.data(fd + 50);
    const auto *fd_51 = buffer.data(fd + 51);
    const auto *fd_52 = buffer.data(fd + 52);
    const auto *fd_53 = buffer.data(fd + 53);
    const auto *fd_54 = buffer.data(fd + 54);
    const auto *fd_55 = buffer.data(fd + 55);
    const auto *fd_56 = buffer.data(fd + 56);
    const auto *fd_57 = buffer.data(fd + 57);
    const auto *fd_58 = buffer.data(fd + 58);
    const auto *fd_59 = buffer.data(fd + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, fd_12, fd_13, fd_14, fd_15, \
                         fd_16, fd_17, fd_24, fd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_12[k];

        t_1[k] = f_0 * fd_13[k];

        t_2[k] = f_0 * fd_14[k];

        t_3[k] = f_0 * fd_15[k];

        t_4[k] = f_0 * fd_16[k];

        t_5[k] = f_0 * fd_17[k];

        t_6[k] = f_0 * fd_24[k];

        t_7[k] = f_0 * fd_25[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, pd_0, pd_1, fd_26, fd_27, fd_28, \
                         fd_29, fd_30, fd_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * fd_26[k];

        t_9[k] = f_0 * fd_27[k];

        t_10[k] = f_0 * fd_28[k];

        t_11[k] = f_0 * fd_29[k];

        t_12[k] = -pd_0[k]
                  + f_0 * fd_30[k];

        t_13[k] = -pd_1[k]
                  + f_0 * fd_31[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, t_19, pd_2, pd_3, pd_4, pd_5, fd_32, \
                         fd_33, fd_34, fd_35, fd_42, fd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -pd_2[k]
                  + f_0 * fd_32[k];

        t_15[k] = -pd_3[k]
                  + f_0 * fd_33[k];

        t_16[k] = -pd_4[k]
                  + f_0 * fd_34[k];

        t_17[k] = -pd_5[k]
                  + f_0 * fd_35[k];

        t_18[k] = f_0 * fd_42[k];

        t_19[k] = f_0 * fd_43[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, pd_6, pd_7, fd_44, fd_45, fd_46, \
                         fd_47, fd_48, fd_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * fd_44[k];

        t_21[k] = f_0 * fd_45[k];

        t_22[k] = f_0 * fd_46[k];

        t_23[k] = f_0 * fd_47[k];

        t_24[k] = -pd_6[k]
                  + f_0 * fd_48[k];

        t_25[k] = -pd_7[k]
                  + f_0 * fd_49[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pd_8, pd_9, pd_10, pd_11, pd_12, fd_50, \
                         fd_51, fd_52, fd_53, fd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -pd_8[k]
                  + f_0 * fd_50[k];

        t_27[k] = -pd_9[k]
                  + f_0 * fd_51[k];

        t_28[k] = -pd_10[k]
                  + f_0 * fd_52[k];

        t_29[k] = -pd_11[k]
                  + f_0 * fd_53[k];

        t_30[k] = -2.0 * pd_12[k]
                  + f_0 * fd_54[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pd_13, pd_14, pd_15, pd_16, pd_17, \
                         fd_55, fd_56, fd_57, fd_58, fd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -2.0 * pd_13[k]
                  + f_0 * fd_55[k];

        t_32[k] = -2.0 * pd_14[k]
                  + f_0 * fd_56[k];

        t_33[k] = -2.0 * pd_15[k]
                  + f_0 * fd_57[k];

        t_34[k] = -2.0 * pd_16[k]
                  + f_0 * fd_58[k];

        t_35[k] = -2.0 * pd_17[k]
                  + f_0 * fd_59[k];
    }
}

}  // namespace simdt2ceri
