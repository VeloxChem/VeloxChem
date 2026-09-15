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


#include "SimdElectronRepulsionGeom10VrrRecPG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_pg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t sg, const size_t dg,
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
    auto *t_36 = buffer.data(target + 36);
    auto *t_37 = buffer.data(target + 37);
    auto *t_38 = buffer.data(target + 38);
    auto *t_39 = buffer.data(target + 39);
    auto *t_40 = buffer.data(target + 40);
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);
    const auto *sg_12 = buffer.data(sg + 12);
    const auto *sg_13 = buffer.data(sg + 13);
    const auto *sg_14 = buffer.data(sg + 14);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_40 = buffer.data(dg + 40);
    const auto *dg_41 = buffer.data(dg + 41);
    const auto *dg_42 = buffer.data(dg + 42);
    const auto *dg_43 = buffer.data(dg + 43);
    const auto *dg_44 = buffer.data(dg + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, sg_0, sg_1, sg_2, sg_3, sg_4, dg_0, dg_1, \
                         dg_2, dg_3, dg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -sg_0[k]
                 + f_0 * dg_0[k];

        t_1[k] = -sg_1[k]
                 + f_0 * dg_1[k];

        t_2[k] = -sg_2[k]
                 + f_0 * dg_2[k];

        t_3[k] = -sg_3[k]
                 + f_0 * dg_3[k];

        t_4[k] = -sg_4[k]
                 + f_0 * dg_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, sg_5, sg_6, sg_7, sg_8, sg_9, dg_5, dg_6, \
                         dg_7, dg_8, dg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -sg_5[k]
                 + f_0 * dg_5[k];

        t_6[k] = -sg_6[k]
                 + f_0 * dg_6[k];

        t_7[k] = -sg_7[k]
                 + f_0 * dg_7[k];

        t_8[k] = -sg_8[k]
                 + f_0 * dg_8[k];

        t_9[k] = -sg_9[k]
                 + f_0 * dg_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, sg_10, sg_11, sg_12, sg_13, sg_14, \
                         dg_10, dg_11, dg_12, dg_13, dg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -sg_10[k]
                  + f_0 * dg_10[k];

        t_11[k] = -sg_11[k]
                  + f_0 * dg_11[k];

        t_12[k] = -sg_12[k]
                  + f_0 * dg_12[k];

        t_13[k] = -sg_13[k]
                  + f_0 * dg_13[k];

        t_14[k] = -sg_14[k]
                  + f_0 * dg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, t_20, t_21, t_22, dg_15, dg_16, dg_17, \
                         dg_18, dg_19, dg_20, dg_21, dg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * dg_15[k];

        t_16[k] = f_0 * dg_16[k];

        t_17[k] = f_0 * dg_17[k];

        t_18[k] = f_0 * dg_18[k];

        t_19[k] = f_0 * dg_19[k];

        t_20[k] = f_0 * dg_20[k];

        t_21[k] = f_0 * dg_21[k];

        t_22[k] = f_0 * dg_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, t_28, t_29, t_30, dg_23, dg_24, dg_25, \
                         dg_26, dg_27, dg_28, dg_29, dg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * dg_23[k];

        t_24[k] = f_0 * dg_24[k];

        t_25[k] = f_0 * dg_25[k];

        t_26[k] = f_0 * dg_26[k];

        t_27[k] = f_0 * dg_27[k];

        t_28[k] = f_0 * dg_28[k];

        t_29[k] = f_0 * dg_29[k];

        t_30[k] = f_0 * dg_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, t_37, t_38, dg_31, dg_32, dg_33, \
                         dg_34, dg_35, dg_36, dg_37, dg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * dg_31[k];

        t_32[k] = f_0 * dg_32[k];

        t_33[k] = f_0 * dg_33[k];

        t_34[k] = f_0 * dg_34[k];

        t_35[k] = f_0 * dg_35[k];

        t_36[k] = f_0 * dg_36[k];

        t_37[k] = f_0 * dg_37[k];

        t_38[k] = f_0 * dg_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, t_44, dg_39, dg_40, dg_41, dg_42, \
                         dg_43, dg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * dg_39[k];

        t_40[k] = f_0 * dg_40[k];

        t_41[k] = f_0 * dg_41[k];

        t_42[k] = f_0 * dg_42[k];

        t_43[k] = f_0 * dg_43[k];

        t_44[k] = f_0 * dg_44[k];
    }
}

auto
compute_prim_geom_10_pg_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t sg, const size_t dg,
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
    auto *t_36 = buffer.data(target + 36);
    auto *t_37 = buffer.data(target + 37);
    auto *t_38 = buffer.data(target + 38);
    auto *t_39 = buffer.data(target + 39);
    auto *t_40 = buffer.data(target + 40);
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);
    const auto *sg_12 = buffer.data(sg + 12);
    const auto *sg_13 = buffer.data(sg + 13);
    const auto *sg_14 = buffer.data(sg + 14);

    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_46 = buffer.data(dg + 46);
    const auto *dg_47 = buffer.data(dg + 47);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_49 = buffer.data(dg + 49);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_51 = buffer.data(dg + 51);
    const auto *dg_52 = buffer.data(dg + 52);
    const auto *dg_53 = buffer.data(dg + 53);
    const auto *dg_54 = buffer.data(dg + 54);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_56 = buffer.data(dg + 56);
    const auto *dg_57 = buffer.data(dg + 57);
    const auto *dg_58 = buffer.data(dg + 58);
    const auto *dg_59 = buffer.data(dg + 59);
    const auto *dg_60 = buffer.data(dg + 60);
    const auto *dg_61 = buffer.data(dg + 61);
    const auto *dg_62 = buffer.data(dg + 62);
    const auto *dg_63 = buffer.data(dg + 63);
    const auto *dg_64 = buffer.data(dg + 64);
    const auto *dg_65 = buffer.data(dg + 65);
    const auto *dg_66 = buffer.data(dg + 66);
    const auto *dg_67 = buffer.data(dg + 67);
    const auto *dg_68 = buffer.data(dg + 68);
    const auto *dg_69 = buffer.data(dg + 69);
    const auto *dg_70 = buffer.data(dg + 70);
    const auto *dg_71 = buffer.data(dg + 71);
    const auto *dg_72 = buffer.data(dg + 72);
    const auto *dg_73 = buffer.data(dg + 73);
    const auto *dg_74 = buffer.data(dg + 74);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, dg_15, dg_16, dg_17, dg_18, \
                         dg_19, dg_20, dg_21, dg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_15[k];

        t_1[k] = f_0 * dg_16[k];

        t_2[k] = f_0 * dg_17[k];

        t_3[k] = f_0 * dg_18[k];

        t_4[k] = f_0 * dg_19[k];

        t_5[k] = f_0 * dg_20[k];

        t_6[k] = f_0 * dg_21[k];

        t_7[k] = f_0 * dg_22[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, dg_23, dg_24, dg_25, dg_26, \
                         dg_27, dg_28, dg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_23[k];

        t_9[k] = f_0 * dg_24[k];

        t_10[k] = f_0 * dg_25[k];

        t_11[k] = f_0 * dg_26[k];

        t_12[k] = f_0 * dg_27[k];

        t_13[k] = f_0 * dg_28[k];

        t_14[k] = f_0 * dg_29[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, sg_0, sg_1, sg_2, sg_3, sg_4, dg_45, \
                         dg_46, dg_47, dg_48, dg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -sg_0[k]
                  + f_0 * dg_45[k];

        t_16[k] = -sg_1[k]
                  + f_0 * dg_46[k];

        t_17[k] = -sg_2[k]
                  + f_0 * dg_47[k];

        t_18[k] = -sg_3[k]
                  + f_0 * dg_48[k];

        t_19[k] = -sg_4[k]
                  + f_0 * dg_49[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, sg_5, sg_6, sg_7, sg_8, sg_9, dg_50, \
                         dg_51, dg_52, dg_53, dg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -sg_5[k]
                  + f_0 * dg_50[k];

        t_21[k] = -sg_6[k]
                  + f_0 * dg_51[k];

        t_22[k] = -sg_7[k]
                  + f_0 * dg_52[k];

        t_23[k] = -sg_8[k]
                  + f_0 * dg_53[k];

        t_24[k] = -sg_9[k]
                  + f_0 * dg_54[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, sg_10, sg_11, sg_12, sg_13, sg_14, \
                         dg_55, dg_56, dg_57, dg_58, dg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -sg_10[k]
                  + f_0 * dg_55[k];

        t_26[k] = -sg_11[k]
                  + f_0 * dg_56[k];

        t_27[k] = -sg_12[k]
                  + f_0 * dg_57[k];

        t_28[k] = -sg_13[k]
                  + f_0 * dg_58[k];

        t_29[k] = -sg_14[k]
                  + f_0 * dg_59[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, t_36, t_37, dg_60, dg_61, dg_62, \
                         dg_63, dg_64, dg_65, dg_66, dg_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * dg_60[k];

        t_31[k] = f_0 * dg_61[k];

        t_32[k] = f_0 * dg_62[k];

        t_33[k] = f_0 * dg_63[k];

        t_34[k] = f_0 * dg_64[k];

        t_35[k] = f_0 * dg_65[k];

        t_36[k] = f_0 * dg_66[k];

        t_37[k] = f_0 * dg_67[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, t_44, dg_68, dg_69, dg_70, dg_71, \
                         dg_72, dg_73, dg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * dg_68[k];

        t_39[k] = f_0 * dg_69[k];

        t_40[k] = f_0 * dg_70[k];

        t_41[k] = f_0 * dg_71[k];

        t_42[k] = f_0 * dg_72[k];

        t_43[k] = f_0 * dg_73[k];

        t_44[k] = f_0 * dg_74[k];
    }
}

auto
compute_prim_geom_10_pg_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t sg, const size_t dg,
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
    auto *t_36 = buffer.data(target + 36);
    auto *t_37 = buffer.data(target + 37);
    auto *t_38 = buffer.data(target + 38);
    auto *t_39 = buffer.data(target + 39);
    auto *t_40 = buffer.data(target + 40);
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);
    const auto *sg_12 = buffer.data(sg + 12);
    const auto *sg_13 = buffer.data(sg + 13);
    const auto *sg_14 = buffer.data(sg + 14);

    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_40 = buffer.data(dg + 40);
    const auto *dg_41 = buffer.data(dg + 41);
    const auto *dg_42 = buffer.data(dg + 42);
    const auto *dg_43 = buffer.data(dg + 43);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_60 = buffer.data(dg + 60);
    const auto *dg_61 = buffer.data(dg + 61);
    const auto *dg_62 = buffer.data(dg + 62);
    const auto *dg_63 = buffer.data(dg + 63);
    const auto *dg_64 = buffer.data(dg + 64);
    const auto *dg_65 = buffer.data(dg + 65);
    const auto *dg_66 = buffer.data(dg + 66);
    const auto *dg_67 = buffer.data(dg + 67);
    const auto *dg_68 = buffer.data(dg + 68);
    const auto *dg_69 = buffer.data(dg + 69);
    const auto *dg_70 = buffer.data(dg + 70);
    const auto *dg_71 = buffer.data(dg + 71);
    const auto *dg_72 = buffer.data(dg + 72);
    const auto *dg_73 = buffer.data(dg + 73);
    const auto *dg_74 = buffer.data(dg + 74);
    const auto *dg_75 = buffer.data(dg + 75);
    const auto *dg_76 = buffer.data(dg + 76);
    const auto *dg_77 = buffer.data(dg + 77);
    const auto *dg_78 = buffer.data(dg + 78);
    const auto *dg_79 = buffer.data(dg + 79);
    const auto *dg_80 = buffer.data(dg + 80);
    const auto *dg_81 = buffer.data(dg + 81);
    const auto *dg_82 = buffer.data(dg + 82);
    const auto *dg_83 = buffer.data(dg + 83);
    const auto *dg_84 = buffer.data(dg + 84);
    const auto *dg_85 = buffer.data(dg + 85);
    const auto *dg_86 = buffer.data(dg + 86);
    const auto *dg_87 = buffer.data(dg + 87);
    const auto *dg_88 = buffer.data(dg + 88);
    const auto *dg_89 = buffer.data(dg + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, dg_30, dg_31, dg_32, dg_33, \
                         dg_34, dg_35, dg_36, dg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_30[k];

        t_1[k] = f_0 * dg_31[k];

        t_2[k] = f_0 * dg_32[k];

        t_3[k] = f_0 * dg_33[k];

        t_4[k] = f_0 * dg_34[k];

        t_5[k] = f_0 * dg_35[k];

        t_6[k] = f_0 * dg_36[k];

        t_7[k] = f_0 * dg_37[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, dg_38, dg_39, dg_40, \
                         dg_41, dg_42, dg_43, dg_44, dg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dg_38[k];

        t_9[k] = f_0 * dg_39[k];

        t_10[k] = f_0 * dg_40[k];

        t_11[k] = f_0 * dg_41[k];

        t_12[k] = f_0 * dg_42[k];

        t_13[k] = f_0 * dg_43[k];

        t_14[k] = f_0 * dg_44[k];

        t_15[k] = f_0 * dg_60[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, dg_61, dg_62, dg_63, \
                         dg_64, dg_65, dg_66, dg_67, dg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * dg_61[k];

        t_17[k] = f_0 * dg_62[k];

        t_18[k] = f_0 * dg_63[k];

        t_19[k] = f_0 * dg_64[k];

        t_20[k] = f_0 * dg_65[k];

        t_21[k] = f_0 * dg_66[k];

        t_22[k] = f_0 * dg_67[k];

        t_23[k] = f_0 * dg_68[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, sg_0, dg_69, dg_70, dg_71, \
                         dg_72, dg_73, dg_74, dg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * dg_69[k];

        t_25[k] = f_0 * dg_70[k];

        t_26[k] = f_0 * dg_71[k];

        t_27[k] = f_0 * dg_72[k];

        t_28[k] = f_0 * dg_73[k];

        t_29[k] = f_0 * dg_74[k];

        t_30[k] = -sg_0[k]
                  + f_0 * dg_75[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, sg_1, sg_2, sg_3, sg_4, sg_5, dg_76, \
                         dg_77, dg_78, dg_79, dg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -sg_1[k]
                  + f_0 * dg_76[k];

        t_32[k] = -sg_2[k]
                  + f_0 * dg_77[k];

        t_33[k] = -sg_3[k]
                  + f_0 * dg_78[k];

        t_34[k] = -sg_4[k]
                  + f_0 * dg_79[k];

        t_35[k] = -sg_5[k]
                  + f_0 * dg_80[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, sg_6, sg_7, sg_8, sg_9, sg_10, dg_81, \
                         dg_82, dg_83, dg_84, dg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -sg_6[k]
                  + f_0 * dg_81[k];

        t_37[k] = -sg_7[k]
                  + f_0 * dg_82[k];

        t_38[k] = -sg_8[k]
                  + f_0 * dg_83[k];

        t_39[k] = -sg_9[k]
                  + f_0 * dg_84[k];

        t_40[k] = -sg_10[k]
                  + f_0 * dg_85[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, sg_11, sg_12, sg_13, sg_14, dg_86, dg_87, \
                         dg_88, dg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -sg_11[k]
                  + f_0 * dg_86[k];

        t_42[k] = -sg_12[k]
                  + f_0 * dg_87[k];

        t_43[k] = -sg_13[k]
                  + f_0 * dg_88[k];

        t_44[k] = -sg_14[k]
                  + f_0 * dg_89[k];
    }
}

}  // namespace simdt2ceri
