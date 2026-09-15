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


#include "SimdElectronRepulsionGeom10VrrRecPF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_pf_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t sf, const size_t df,
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

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_6 = buffer.data(sf + 6);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);
    const auto *sf_9 = buffer.data(sf + 9);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_26 = buffer.data(df + 26);
    const auto *df_27 = buffer.data(df + 27);
    const auto *df_28 = buffer.data(df + 28);
    const auto *df_29 = buffer.data(df + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, sf_0, sf_1, sf_2, sf_3, sf_4, df_0, df_1, \
                         df_2, df_3, df_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -sf_0[k]
                 + f_0 * df_0[k];

        t_1[k] = -sf_1[k]
                 + f_0 * df_1[k];

        t_2[k] = -sf_2[k]
                 + f_0 * df_2[k];

        t_3[k] = -sf_3[k]
                 + f_0 * df_3[k];

        t_4[k] = -sf_4[k]
                 + f_0 * df_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, sf_5, sf_6, sf_7, sf_8, sf_9, df_5, df_6, \
                         df_7, df_8, df_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -sf_5[k]
                 + f_0 * df_5[k];

        t_6[k] = -sf_6[k]
                 + f_0 * df_6[k];

        t_7[k] = -sf_7[k]
                 + f_0 * df_7[k];

        t_8[k] = -sf_8[k]
                 + f_0 * df_8[k];

        t_9[k] = -sf_9[k]
                 + f_0 * df_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, t_15, t_16, t_17, df_10, df_11, df_12, \
                         df_13, df_14, df_15, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * df_10[k];

        t_11[k] = f_0 * df_11[k];

        t_12[k] = f_0 * df_12[k];

        t_13[k] = f_0 * df_13[k];

        t_14[k] = f_0 * df_14[k];

        t_15[k] = f_0 * df_15[k];

        t_16[k] = f_0 * df_16[k];

        t_17[k] = f_0 * df_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, t_24, t_25, df_18, df_19, df_20, \
                         df_21, df_22, df_23, df_24, df_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * df_18[k];

        t_19[k] = f_0 * df_19[k];

        t_20[k] = f_0 * df_20[k];

        t_21[k] = f_0 * df_21[k];

        t_22[k] = f_0 * df_22[k];

        t_23[k] = f_0 * df_23[k];

        t_24[k] = f_0 * df_24[k];

        t_25[k] = f_0 * df_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, df_26, df_27, df_28, \
                         df_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * df_26[k];

        t_27[k] = f_0 * df_27[k];

        t_28[k] = f_0 * df_28[k];

        t_29[k] = f_0 * df_29[k];
    }
}

auto
compute_prim_geom_10_pf_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t sf, const size_t df,
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

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_6 = buffer.data(sf + 6);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);
    const auto *sf_9 = buffer.data(sf + 9);

    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_30 = buffer.data(df + 30);
    const auto *df_31 = buffer.data(df + 31);
    const auto *df_32 = buffer.data(df + 32);
    const auto *df_33 = buffer.data(df + 33);
    const auto *df_34 = buffer.data(df + 34);
    const auto *df_35 = buffer.data(df + 35);
    const auto *df_36 = buffer.data(df + 36);
    const auto *df_37 = buffer.data(df + 37);
    const auto *df_38 = buffer.data(df + 38);
    const auto *df_39 = buffer.data(df + 39);
    const auto *df_40 = buffer.data(df + 40);
    const auto *df_41 = buffer.data(df + 41);
    const auto *df_42 = buffer.data(df + 42);
    const auto *df_43 = buffer.data(df + 43);
    const auto *df_44 = buffer.data(df + 44);
    const auto *df_45 = buffer.data(df + 45);
    const auto *df_46 = buffer.data(df + 46);
    const auto *df_47 = buffer.data(df + 47);
    const auto *df_48 = buffer.data(df + 48);
    const auto *df_49 = buffer.data(df + 49);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, df_10, df_11, df_12, df_13, \
                         df_14, df_15, df_16, df_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_10[k];

        t_1[k] = f_0 * df_11[k];

        t_2[k] = f_0 * df_12[k];

        t_3[k] = f_0 * df_13[k];

        t_4[k] = f_0 * df_14[k];

        t_5[k] = f_0 * df_15[k];

        t_6[k] = f_0 * df_16[k];

        t_7[k] = f_0 * df_17[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, sf_0, sf_1, sf_2, sf_3, df_18, \
                         df_19, df_30, df_31, df_32, df_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * df_18[k];

        t_9[k] = f_0 * df_19[k];

        t_10[k] = -sf_0[k]
                  + f_0 * df_30[k];

        t_11[k] = -sf_1[k]
                  + f_0 * df_31[k];

        t_12[k] = -sf_2[k]
                  + f_0 * df_32[k];

        t_13[k] = -sf_3[k]
                  + f_0 * df_33[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, sf_4, sf_5, sf_6, sf_7, sf_8, df_34, \
                         df_35, df_36, df_37, df_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -sf_4[k]
                  + f_0 * df_34[k];

        t_15[k] = -sf_5[k]
                  + f_0 * df_35[k];

        t_16[k] = -sf_6[k]
                  + f_0 * df_36[k];

        t_17[k] = -sf_7[k]
                  + f_0 * df_37[k];

        t_18[k] = -sf_8[k]
                  + f_0 * df_38[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, t_25, sf_9, df_39, df_40, df_41, \
                         df_42, df_43, df_44, df_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -sf_9[k]
                  + f_0 * df_39[k];

        t_20[k] = f_0 * df_40[k];

        t_21[k] = f_0 * df_41[k];

        t_22[k] = f_0 * df_42[k];

        t_23[k] = f_0 * df_43[k];

        t_24[k] = f_0 * df_44[k];

        t_25[k] = f_0 * df_45[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, df_46, df_47, df_48, \
                         df_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * df_46[k];

        t_27[k] = f_0 * df_47[k];

        t_28[k] = f_0 * df_48[k];

        t_29[k] = f_0 * df_49[k];
    }
}

auto
compute_prim_geom_10_pf_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t sf, const size_t df,
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

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_6 = buffer.data(sf + 6);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);
    const auto *sf_9 = buffer.data(sf + 9);

    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_26 = buffer.data(df + 26);
    const auto *df_27 = buffer.data(df + 27);
    const auto *df_28 = buffer.data(df + 28);
    const auto *df_29 = buffer.data(df + 29);
    const auto *df_40 = buffer.data(df + 40);
    const auto *df_41 = buffer.data(df + 41);
    const auto *df_42 = buffer.data(df + 42);
    const auto *df_43 = buffer.data(df + 43);
    const auto *df_44 = buffer.data(df + 44);
    const auto *df_45 = buffer.data(df + 45);
    const auto *df_46 = buffer.data(df + 46);
    const auto *df_47 = buffer.data(df + 47);
    const auto *df_48 = buffer.data(df + 48);
    const auto *df_49 = buffer.data(df + 49);
    const auto *df_50 = buffer.data(df + 50);
    const auto *df_51 = buffer.data(df + 51);
    const auto *df_52 = buffer.data(df + 52);
    const auto *df_53 = buffer.data(df + 53);
    const auto *df_54 = buffer.data(df + 54);
    const auto *df_55 = buffer.data(df + 55);
    const auto *df_56 = buffer.data(df + 56);
    const auto *df_57 = buffer.data(df + 57);
    const auto *df_58 = buffer.data(df + 58);
    const auto *df_59 = buffer.data(df + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, df_20, df_21, df_22, df_23, \
                         df_24, df_25, df_26, df_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_20[k];

        t_1[k] = f_0 * df_21[k];

        t_2[k] = f_0 * df_22[k];

        t_3[k] = f_0 * df_23[k];

        t_4[k] = f_0 * df_24[k];

        t_5[k] = f_0 * df_25[k];

        t_6[k] = f_0 * df_26[k];

        t_7[k] = f_0 * df_27[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, df_28, df_29, df_40, \
                         df_41, df_42, df_43, df_44, df_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * df_28[k];

        t_9[k] = f_0 * df_29[k];

        t_10[k] = f_0 * df_40[k];

        t_11[k] = f_0 * df_41[k];

        t_12[k] = f_0 * df_42[k];

        t_13[k] = f_0 * df_43[k];

        t_14[k] = f_0 * df_44[k];

        t_15[k] = f_0 * df_45[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, sf_0, sf_1, df_46, df_47, df_48, \
                         df_49, df_50, df_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * df_46[k];

        t_17[k] = f_0 * df_47[k];

        t_18[k] = f_0 * df_48[k];

        t_19[k] = f_0 * df_49[k];

        t_20[k] = -sf_0[k]
                  + f_0 * df_50[k];

        t_21[k] = -sf_1[k]
                  + f_0 * df_51[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, sf_2, sf_3, sf_4, sf_5, sf_6, df_52, \
                         df_53, df_54, df_55, df_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -sf_2[k]
                  + f_0 * df_52[k];

        t_23[k] = -sf_3[k]
                  + f_0 * df_53[k];

        t_24[k] = -sf_4[k]
                  + f_0 * df_54[k];

        t_25[k] = -sf_5[k]
                  + f_0 * df_55[k];

        t_26[k] = -sf_6[k]
                  + f_0 * df_56[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, sf_7, sf_8, sf_9, df_57, df_58, \
                         df_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -sf_7[k]
                  + f_0 * df_57[k];

        t_28[k] = -sf_8[k]
                  + f_0 * df_58[k];

        t_29[k] = -sf_9[k]
                  + f_0 * df_59[k];
    }
}

}  // namespace simdt2ceri
