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


#include "SimdElectronRepulsionGeom10VrrRecPH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_ph_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t sh, const size_t dh,
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
    auto *t_55 = buffer.data(target + 55);
    auto *t_56 = buffer.data(target + 56);
    auto *t_57 = buffer.data(target + 57);
    auto *t_58 = buffer.data(target + 58);
    auto *t_59 = buffer.data(target + 59);
    auto *t_60 = buffer.data(target + 60);
    auto *t_61 = buffer.data(target + 61);
    auto *t_62 = buffer.data(target + 62);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_13 = buffer.data(sh + 13);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_15 = buffer.data(sh + 15);
    const auto *sh_16 = buffer.data(sh + 16);
    const auto *sh_17 = buffer.data(sh + 17);
    const auto *sh_18 = buffer.data(sh + 18);
    const auto *sh_19 = buffer.data(sh + 19);
    const auto *sh_20 = buffer.data(sh + 20);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_4 = buffer.data(dh + 4);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_11 = buffer.data(dh + 11);
    const auto *dh_12 = buffer.data(dh + 12);
    const auto *dh_13 = buffer.data(dh + 13);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_16 = buffer.data(dh + 16);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_37 = buffer.data(dh + 37);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_41 = buffer.data(dh + 41);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_43 = buffer.data(dh + 43);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_45 = buffer.data(dh + 45);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_48 = buffer.data(dh + 48);
    const auto *dh_49 = buffer.data(dh + 49);
    const auto *dh_50 = buffer.data(dh + 50);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_54 = buffer.data(dh + 54);
    const auto *dh_55 = buffer.data(dh + 55);
    const auto *dh_56 = buffer.data(dh + 56);
    const auto *dh_57 = buffer.data(dh + 57);
    const auto *dh_58 = buffer.data(dh + 58);
    const auto *dh_59 = buffer.data(dh + 59);
    const auto *dh_60 = buffer.data(dh + 60);
    const auto *dh_61 = buffer.data(dh + 61);
    const auto *dh_62 = buffer.data(dh + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, sh_0, sh_1, sh_2, sh_3, sh_4, dh_0, dh_1, \
                         dh_2, dh_3, dh_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -sh_0[k]
                 + f_0 * dh_0[k];

        t_1[k] = -sh_1[k]
                 + f_0 * dh_1[k];

        t_2[k] = -sh_2[k]
                 + f_0 * dh_2[k];

        t_3[k] = -sh_3[k]
                 + f_0 * dh_3[k];

        t_4[k] = -sh_4[k]
                 + f_0 * dh_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, sh_5, sh_6, sh_7, sh_8, sh_9, dh_5, dh_6, \
                         dh_7, dh_8, dh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -sh_5[k]
                 + f_0 * dh_5[k];

        t_6[k] = -sh_6[k]
                 + f_0 * dh_6[k];

        t_7[k] = -sh_7[k]
                 + f_0 * dh_7[k];

        t_8[k] = -sh_8[k]
                 + f_0 * dh_8[k];

        t_9[k] = -sh_9[k]
                 + f_0 * dh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, sh_10, sh_11, sh_12, sh_13, sh_14, \
                         dh_10, dh_11, dh_12, dh_13, dh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -sh_10[k]
                  + f_0 * dh_10[k];

        t_11[k] = -sh_11[k]
                  + f_0 * dh_11[k];

        t_12[k] = -sh_12[k]
                  + f_0 * dh_12[k];

        t_13[k] = -sh_13[k]
                  + f_0 * dh_13[k];

        t_14[k] = -sh_14[k]
                  + f_0 * dh_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, sh_15, sh_16, sh_17, sh_18, sh_19, \
                         dh_15, dh_16, dh_17, dh_18, dh_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -sh_15[k]
                  + f_0 * dh_15[k];

        t_16[k] = -sh_16[k]
                  + f_0 * dh_16[k];

        t_17[k] = -sh_17[k]
                  + f_0 * dh_17[k];

        t_18[k] = -sh_18[k]
                  + f_0 * dh_18[k];

        t_19[k] = -sh_19[k]
                  + f_0 * dh_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, t_26, sh_20, dh_20, dh_21, dh_22, \
                         dh_23, dh_24, dh_25, dh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -sh_20[k]
                  + f_0 * dh_20[k];

        t_21[k] = f_0 * dh_21[k];

        t_22[k] = f_0 * dh_22[k];

        t_23[k] = f_0 * dh_23[k];

        t_24[k] = f_0 * dh_24[k];

        t_25[k] = f_0 * dh_25[k];

        t_26[k] = f_0 * dh_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, t_32, t_33, t_34, dh_27, dh_28, dh_29, \
                         dh_30, dh_31, dh_32, dh_33, dh_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * dh_27[k];

        t_28[k] = f_0 * dh_28[k];

        t_29[k] = f_0 * dh_29[k];

        t_30[k] = f_0 * dh_30[k];

        t_31[k] = f_0 * dh_31[k];

        t_32[k] = f_0 * dh_32[k];

        t_33[k] = f_0 * dh_33[k];

        t_34[k] = f_0 * dh_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, t_40, t_41, t_42, dh_35, dh_36, dh_37, \
                         dh_38, dh_39, dh_40, dh_41, dh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * dh_35[k];

        t_36[k] = f_0 * dh_36[k];

        t_37[k] = f_0 * dh_37[k];

        t_38[k] = f_0 * dh_38[k];

        t_39[k] = f_0 * dh_39[k];

        t_40[k] = f_0 * dh_40[k];

        t_41[k] = f_0 * dh_41[k];

        t_42[k] = f_0 * dh_42[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, t_48, t_49, t_50, dh_43, dh_44, dh_45, \
                         dh_46, dh_47, dh_48, dh_49, dh_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_0 * dh_43[k];

        t_44[k] = f_0 * dh_44[k];

        t_45[k] = f_0 * dh_45[k];

        t_46[k] = f_0 * dh_46[k];

        t_47[k] = f_0 * dh_47[k];

        t_48[k] = f_0 * dh_48[k];

        t_49[k] = f_0 * dh_49[k];

        t_50[k] = f_0 * dh_50[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, t_57, t_58, dh_51, dh_52, dh_53, \
                         dh_54, dh_55, dh_56, dh_57, dh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * dh_51[k];

        t_52[k] = f_0 * dh_52[k];

        t_53[k] = f_0 * dh_53[k];

        t_54[k] = f_0 * dh_54[k];

        t_55[k] = f_0 * dh_55[k];

        t_56[k] = f_0 * dh_56[k];

        t_57[k] = f_0 * dh_57[k];

        t_58[k] = f_0 * dh_58[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, dh_59, dh_60, dh_61, \
                         dh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_0 * dh_59[k];

        t_60[k] = f_0 * dh_60[k];

        t_61[k] = f_0 * dh_61[k];

        t_62[k] = f_0 * dh_62[k];
    }
}

auto
compute_prim_geom_10_ph_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t sh, const size_t dh,
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
    auto *t_55 = buffer.data(target + 55);
    auto *t_56 = buffer.data(target + 56);
    auto *t_57 = buffer.data(target + 57);
    auto *t_58 = buffer.data(target + 58);
    auto *t_59 = buffer.data(target + 59);
    auto *t_60 = buffer.data(target + 60);
    auto *t_61 = buffer.data(target + 61);
    auto *t_62 = buffer.data(target + 62);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_13 = buffer.data(sh + 13);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_15 = buffer.data(sh + 15);
    const auto *sh_16 = buffer.data(sh + 16);
    const auto *sh_17 = buffer.data(sh + 17);
    const auto *sh_18 = buffer.data(sh + 18);
    const auto *sh_19 = buffer.data(sh + 19);
    const auto *sh_20 = buffer.data(sh + 20);

    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_23 = buffer.data(dh + 23);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_25 = buffer.data(dh + 25);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_28 = buffer.data(dh + 28);
    const auto *dh_29 = buffer.data(dh + 29);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_32 = buffer.data(dh + 32);
    const auto *dh_33 = buffer.data(dh + 33);
    const auto *dh_34 = buffer.data(dh + 34);
    const auto *dh_35 = buffer.data(dh + 35);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_37 = buffer.data(dh + 37);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_41 = buffer.data(dh + 41);
    const auto *dh_63 = buffer.data(dh + 63);
    const auto *dh_64 = buffer.data(dh + 64);
    const auto *dh_65 = buffer.data(dh + 65);
    const auto *dh_66 = buffer.data(dh + 66);
    const auto *dh_67 = buffer.data(dh + 67);
    const auto *dh_68 = buffer.data(dh + 68);
    const auto *dh_69 = buffer.data(dh + 69);
    const auto *dh_70 = buffer.data(dh + 70);
    const auto *dh_71 = buffer.data(dh + 71);
    const auto *dh_72 = buffer.data(dh + 72);
    const auto *dh_73 = buffer.data(dh + 73);
    const auto *dh_74 = buffer.data(dh + 74);
    const auto *dh_75 = buffer.data(dh + 75);
    const auto *dh_76 = buffer.data(dh + 76);
    const auto *dh_77 = buffer.data(dh + 77);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_79 = buffer.data(dh + 79);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_82 = buffer.data(dh + 82);
    const auto *dh_83 = buffer.data(dh + 83);
    const auto *dh_84 = buffer.data(dh + 84);
    const auto *dh_85 = buffer.data(dh + 85);
    const auto *dh_86 = buffer.data(dh + 86);
    const auto *dh_87 = buffer.data(dh + 87);
    const auto *dh_88 = buffer.data(dh + 88);
    const auto *dh_89 = buffer.data(dh + 89);
    const auto *dh_90 = buffer.data(dh + 90);
    const auto *dh_91 = buffer.data(dh + 91);
    const auto *dh_92 = buffer.data(dh + 92);
    const auto *dh_93 = buffer.data(dh + 93);
    const auto *dh_94 = buffer.data(dh + 94);
    const auto *dh_95 = buffer.data(dh + 95);
    const auto *dh_96 = buffer.data(dh + 96);
    const auto *dh_97 = buffer.data(dh + 97);
    const auto *dh_98 = buffer.data(dh + 98);
    const auto *dh_99 = buffer.data(dh + 99);
    const auto *dh_100 = buffer.data(dh + 100);
    const auto *dh_101 = buffer.data(dh + 101);
    const auto *dh_102 = buffer.data(dh + 102);
    const auto *dh_103 = buffer.data(dh + 103);
    const auto *dh_104 = buffer.data(dh + 104);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, dh_21, dh_22, dh_23, dh_24, \
                         dh_25, dh_26, dh_27, dh_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dh_21[k];

        t_1[k] = f_0 * dh_22[k];

        t_2[k] = f_0 * dh_23[k];

        t_3[k] = f_0 * dh_24[k];

        t_4[k] = f_0 * dh_25[k];

        t_5[k] = f_0 * dh_26[k];

        t_6[k] = f_0 * dh_27[k];

        t_7[k] = f_0 * dh_28[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, dh_29, dh_30, dh_31, \
                         dh_32, dh_33, dh_34, dh_35, dh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dh_29[k];

        t_9[k] = f_0 * dh_30[k];

        t_10[k] = f_0 * dh_31[k];

        t_11[k] = f_0 * dh_32[k];

        t_12[k] = f_0 * dh_33[k];

        t_13[k] = f_0 * dh_34[k];

        t_14[k] = f_0 * dh_35[k];

        t_15[k] = f_0 * dh_36[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, sh_0, sh_1, dh_37, dh_38, \
                         dh_39, dh_40, dh_41, dh_63, dh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * dh_37[k];

        t_17[k] = f_0 * dh_38[k];

        t_18[k] = f_0 * dh_39[k];

        t_19[k] = f_0 * dh_40[k];

        t_20[k] = f_0 * dh_41[k];

        t_21[k] = -sh_0[k]
                  + f_0 * dh_63[k];

        t_22[k] = -sh_1[k]
                  + f_0 * dh_64[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, sh_2, sh_3, sh_4, sh_5, sh_6, dh_65, \
                         dh_66, dh_67, dh_68, dh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -sh_2[k]
                  + f_0 * dh_65[k];

        t_24[k] = -sh_3[k]
                  + f_0 * dh_66[k];

        t_25[k] = -sh_4[k]
                  + f_0 * dh_67[k];

        t_26[k] = -sh_5[k]
                  + f_0 * dh_68[k];

        t_27[k] = -sh_6[k]
                  + f_0 * dh_69[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, sh_7, sh_8, sh_9, sh_10, sh_11, dh_70, \
                         dh_71, dh_72, dh_73, dh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -sh_7[k]
                  + f_0 * dh_70[k];

        t_29[k] = -sh_8[k]
                  + f_0 * dh_71[k];

        t_30[k] = -sh_9[k]
                  + f_0 * dh_72[k];

        t_31[k] = -sh_10[k]
                  + f_0 * dh_73[k];

        t_32[k] = -sh_11[k]
                  + f_0 * dh_74[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, sh_12, sh_13, sh_14, sh_15, sh_16, \
                         dh_75, dh_76, dh_77, dh_78, dh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -sh_12[k]
                  + f_0 * dh_75[k];

        t_34[k] = -sh_13[k]
                  + f_0 * dh_76[k];

        t_35[k] = -sh_14[k]
                  + f_0 * dh_77[k];

        t_36[k] = -sh_15[k]
                  + f_0 * dh_78[k];

        t_37[k] = -sh_16[k]
                  + f_0 * dh_79[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, sh_17, sh_18, sh_19, sh_20, \
                         dh_80, dh_81, dh_82, dh_83, dh_84, dh_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -sh_17[k]
                  + f_0 * dh_80[k];

        t_39[k] = -sh_18[k]
                  + f_0 * dh_81[k];

        t_40[k] = -sh_19[k]
                  + f_0 * dh_82[k];

        t_41[k] = -sh_20[k]
                  + f_0 * dh_83[k];

        t_42[k] = f_0 * dh_84[k];

        t_43[k] = f_0 * dh_85[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, t_50, t_51, dh_86, dh_87, dh_88, \
                         dh_89, dh_90, dh_91, dh_92, dh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * dh_86[k];

        t_45[k] = f_0 * dh_87[k];

        t_46[k] = f_0 * dh_88[k];

        t_47[k] = f_0 * dh_89[k];

        t_48[k] = f_0 * dh_90[k];

        t_49[k] = f_0 * dh_91[k];

        t_50[k] = f_0 * dh_92[k];

        t_51[k] = f_0 * dh_93[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, t_57, t_58, t_59, dh_94, dh_95, dh_96, \
                         dh_97, dh_98, dh_99, dh_100, dh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_0 * dh_94[k];

        t_53[k] = f_0 * dh_95[k];

        t_54[k] = f_0 * dh_96[k];

        t_55[k] = f_0 * dh_97[k];

        t_56[k] = f_0 * dh_98[k];

        t_57[k] = f_0 * dh_99[k];

        t_58[k] = f_0 * dh_100[k];

        t_59[k] = f_0 * dh_101[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, dh_102, dh_103, dh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * dh_102[k];

        t_61[k] = f_0 * dh_103[k];

        t_62[k] = f_0 * dh_104[k];
    }
}

auto
compute_prim_geom_10_ph_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t sh, const size_t dh,
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
    auto *t_55 = buffer.data(target + 55);
    auto *t_56 = buffer.data(target + 56);
    auto *t_57 = buffer.data(target + 57);
    auto *t_58 = buffer.data(target + 58);
    auto *t_59 = buffer.data(target + 59);
    auto *t_60 = buffer.data(target + 60);
    auto *t_61 = buffer.data(target + 61);
    auto *t_62 = buffer.data(target + 62);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_13 = buffer.data(sh + 13);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_15 = buffer.data(sh + 15);
    const auto *sh_16 = buffer.data(sh + 16);
    const auto *sh_17 = buffer.data(sh + 17);
    const auto *sh_18 = buffer.data(sh + 18);
    const auto *sh_19 = buffer.data(sh + 19);
    const auto *sh_20 = buffer.data(sh + 20);

    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_43 = buffer.data(dh + 43);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_45 = buffer.data(dh + 45);
    const auto *dh_46 = buffer.data(dh + 46);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_48 = buffer.data(dh + 48);
    const auto *dh_49 = buffer.data(dh + 49);
    const auto *dh_50 = buffer.data(dh + 50);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_52 = buffer.data(dh + 52);
    const auto *dh_53 = buffer.data(dh + 53);
    const auto *dh_54 = buffer.data(dh + 54);
    const auto *dh_55 = buffer.data(dh + 55);
    const auto *dh_56 = buffer.data(dh + 56);
    const auto *dh_57 = buffer.data(dh + 57);
    const auto *dh_58 = buffer.data(dh + 58);
    const auto *dh_59 = buffer.data(dh + 59);
    const auto *dh_60 = buffer.data(dh + 60);
    const auto *dh_61 = buffer.data(dh + 61);
    const auto *dh_62 = buffer.data(dh + 62);
    const auto *dh_84 = buffer.data(dh + 84);
    const auto *dh_85 = buffer.data(dh + 85);
    const auto *dh_86 = buffer.data(dh + 86);
    const auto *dh_87 = buffer.data(dh + 87);
    const auto *dh_88 = buffer.data(dh + 88);
    const auto *dh_89 = buffer.data(dh + 89);
    const auto *dh_90 = buffer.data(dh + 90);
    const auto *dh_91 = buffer.data(dh + 91);
    const auto *dh_92 = buffer.data(dh + 92);
    const auto *dh_93 = buffer.data(dh + 93);
    const auto *dh_94 = buffer.data(dh + 94);
    const auto *dh_95 = buffer.data(dh + 95);
    const auto *dh_96 = buffer.data(dh + 96);
    const auto *dh_97 = buffer.data(dh + 97);
    const auto *dh_98 = buffer.data(dh + 98);
    const auto *dh_99 = buffer.data(dh + 99);
    const auto *dh_100 = buffer.data(dh + 100);
    const auto *dh_101 = buffer.data(dh + 101);
    const auto *dh_102 = buffer.data(dh + 102);
    const auto *dh_103 = buffer.data(dh + 103);
    const auto *dh_104 = buffer.data(dh + 104);
    const auto *dh_105 = buffer.data(dh + 105);
    const auto *dh_106 = buffer.data(dh + 106);
    const auto *dh_107 = buffer.data(dh + 107);
    const auto *dh_108 = buffer.data(dh + 108);
    const auto *dh_109 = buffer.data(dh + 109);
    const auto *dh_110 = buffer.data(dh + 110);
    const auto *dh_111 = buffer.data(dh + 111);
    const auto *dh_112 = buffer.data(dh + 112);
    const auto *dh_113 = buffer.data(dh + 113);
    const auto *dh_114 = buffer.data(dh + 114);
    const auto *dh_115 = buffer.data(dh + 115);
    const auto *dh_116 = buffer.data(dh + 116);
    const auto *dh_117 = buffer.data(dh + 117);
    const auto *dh_118 = buffer.data(dh + 118);
    const auto *dh_119 = buffer.data(dh + 119);
    const auto *dh_120 = buffer.data(dh + 120);
    const auto *dh_121 = buffer.data(dh + 121);
    const auto *dh_122 = buffer.data(dh + 122);
    const auto *dh_123 = buffer.data(dh + 123);
    const auto *dh_124 = buffer.data(dh + 124);
    const auto *dh_125 = buffer.data(dh + 125);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, dh_42, dh_43, dh_44, dh_45, \
                         dh_46, dh_47, dh_48, dh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dh_42[k];

        t_1[k] = f_0 * dh_43[k];

        t_2[k] = f_0 * dh_44[k];

        t_3[k] = f_0 * dh_45[k];

        t_4[k] = f_0 * dh_46[k];

        t_5[k] = f_0 * dh_47[k];

        t_6[k] = f_0 * dh_48[k];

        t_7[k] = f_0 * dh_49[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, dh_50, dh_51, dh_52, \
                         dh_53, dh_54, dh_55, dh_56, dh_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * dh_50[k];

        t_9[k] = f_0 * dh_51[k];

        t_10[k] = f_0 * dh_52[k];

        t_11[k] = f_0 * dh_53[k];

        t_12[k] = f_0 * dh_54[k];

        t_13[k] = f_0 * dh_55[k];

        t_14[k] = f_0 * dh_56[k];

        t_15[k] = f_0 * dh_57[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, dh_58, dh_59, dh_60, \
                         dh_61, dh_62, dh_84, dh_85, dh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * dh_58[k];

        t_17[k] = f_0 * dh_59[k];

        t_18[k] = f_0 * dh_60[k];

        t_19[k] = f_0 * dh_61[k];

        t_20[k] = f_0 * dh_62[k];

        t_21[k] = f_0 * dh_84[k];

        t_22[k] = f_0 * dh_85[k];

        t_23[k] = f_0 * dh_86[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, dh_87, dh_88, dh_89, \
                         dh_90, dh_91, dh_92, dh_93, dh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * dh_87[k];

        t_25[k] = f_0 * dh_88[k];

        t_26[k] = f_0 * dh_89[k];

        t_27[k] = f_0 * dh_90[k];

        t_28[k] = f_0 * dh_91[k];

        t_29[k] = f_0 * dh_92[k];

        t_30[k] = f_0 * dh_93[k];

        t_31[k] = f_0 * dh_94[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, dh_95, dh_96, dh_97, \
                         dh_98, dh_99, dh_100, dh_101, dh_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * dh_95[k];

        t_33[k] = f_0 * dh_96[k];

        t_34[k] = f_0 * dh_97[k];

        t_35[k] = f_0 * dh_98[k];

        t_36[k] = f_0 * dh_99[k];

        t_37[k] = f_0 * dh_100[k];

        t_38[k] = f_0 * dh_101[k];

        t_39[k] = f_0 * dh_102[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, sh_0, sh_1, sh_2, sh_3, dh_103, \
                         dh_104, dh_105, dh_106, dh_107, dh_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * dh_103[k];

        t_41[k] = f_0 * dh_104[k];

        t_42[k] = -sh_0[k]
                  + f_0 * dh_105[k];

        t_43[k] = -sh_1[k]
                  + f_0 * dh_106[k];

        t_44[k] = -sh_2[k]
                  + f_0 * dh_107[k];

        t_45[k] = -sh_3[k]
                  + f_0 * dh_108[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, sh_4, sh_5, sh_6, sh_7, sh_8, dh_109, \
                         dh_110, dh_111, dh_112, dh_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -sh_4[k]
                  + f_0 * dh_109[k];

        t_47[k] = -sh_5[k]
                  + f_0 * dh_110[k];

        t_48[k] = -sh_6[k]
                  + f_0 * dh_111[k];

        t_49[k] = -sh_7[k]
                  + f_0 * dh_112[k];

        t_50[k] = -sh_8[k]
                  + f_0 * dh_113[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, sh_9, sh_10, sh_11, sh_12, sh_13, \
                         dh_114, dh_115, dh_116, dh_117, dh_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -sh_9[k]
                  + f_0 * dh_114[k];

        t_52[k] = -sh_10[k]
                  + f_0 * dh_115[k];

        t_53[k] = -sh_11[k]
                  + f_0 * dh_116[k];

        t_54[k] = -sh_12[k]
                  + f_0 * dh_117[k];

        t_55[k] = -sh_13[k]
                  + f_0 * dh_118[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, sh_14, sh_15, sh_16, sh_17, sh_18, \
                         dh_119, dh_120, dh_121, dh_122, dh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -sh_14[k]
                  + f_0 * dh_119[k];

        t_57[k] = -sh_15[k]
                  + f_0 * dh_120[k];

        t_58[k] = -sh_16[k]
                  + f_0 * dh_121[k];

        t_59[k] = -sh_17[k]
                  + f_0 * dh_122[k];

        t_60[k] = -sh_18[k]
                  + f_0 * dh_123[k];
    }

#pragma omp simd aligned(t_61, t_62, sh_19, sh_20, dh_124, dh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -sh_19[k]
                  + f_0 * dh_124[k];

        t_62[k] = -sh_20[k]
                  + f_0 * dh_125[k];
    }
}

}  // namespace simdt2ceri
