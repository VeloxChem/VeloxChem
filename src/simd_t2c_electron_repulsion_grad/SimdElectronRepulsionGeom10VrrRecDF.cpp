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


#include "SimdElectronRepulsionGeom10VrrRecDF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_df_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t pf, const size_t ff,
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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_12 = buffer.data(pf + 12);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_14 = buffer.data(pf + 14);
    const auto *pf_15 = buffer.data(pf + 15);
    const auto *pf_16 = buffer.data(pf + 16);
    const auto *pf_17 = buffer.data(pf + 17);
    const auto *pf_18 = buffer.data(pf + 18);
    const auto *pf_19 = buffer.data(pf + 19);
    const auto *pf_20 = buffer.data(pf + 20);
    const auto *pf_21 = buffer.data(pf + 21);
    const auto *pf_22 = buffer.data(pf + 22);
    const auto *pf_23 = buffer.data(pf + 23);
    const auto *pf_24 = buffer.data(pf + 24);
    const auto *pf_25 = buffer.data(pf + 25);
    const auto *pf_26 = buffer.data(pf + 26);
    const auto *pf_27 = buffer.data(pf + 27);
    const auto *pf_28 = buffer.data(pf + 28);
    const auto *pf_29 = buffer.data(pf + 29);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_41 = buffer.data(ff + 41);
    const auto *ff_42 = buffer.data(ff + 42);
    const auto *ff_43 = buffer.data(ff + 43);
    const auto *ff_44 = buffer.data(ff + 44);
    const auto *ff_45 = buffer.data(ff + 45);
    const auto *ff_46 = buffer.data(ff + 46);
    const auto *ff_47 = buffer.data(ff + 47);
    const auto *ff_48 = buffer.data(ff + 48);
    const auto *ff_49 = buffer.data(ff + 49);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_51 = buffer.data(ff + 51);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_53 = buffer.data(ff + 53);
    const auto *ff_54 = buffer.data(ff + 54);
    const auto *ff_55 = buffer.data(ff + 55);
    const auto *ff_56 = buffer.data(ff + 56);
    const auto *ff_57 = buffer.data(ff + 57);
    const auto *ff_58 = buffer.data(ff + 58);
    const auto *ff_59 = buffer.data(ff + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pf_0, pf_1, pf_2, pf_3, pf_4, ff_0, ff_1, \
                         ff_2, ff_3, ff_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -2.0 * pf_0[k]
                 + f_0 * ff_0[k];

        t_1[k] = -2.0 * pf_1[k]
                 + f_0 * ff_1[k];

        t_2[k] = -2.0 * pf_2[k]
                 + f_0 * ff_2[k];

        t_3[k] = -2.0 * pf_3[k]
                 + f_0 * ff_3[k];

        t_4[k] = -2.0 * pf_4[k]
                 + f_0 * ff_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pf_5, pf_6, pf_7, pf_8, pf_9, ff_5, ff_6, \
                         ff_7, ff_8, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -2.0 * pf_5[k]
                 + f_0 * ff_5[k];

        t_6[k] = -2.0 * pf_6[k]
                 + f_0 * ff_6[k];

        t_7[k] = -2.0 * pf_7[k]
                 + f_0 * ff_7[k];

        t_8[k] = -2.0 * pf_8[k]
                 + f_0 * ff_8[k];

        t_9[k] = -2.0 * pf_9[k]
                 + f_0 * ff_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pf_10, pf_11, pf_12, pf_13, pf_14, \
                         ff_10, ff_11, ff_12, ff_13, ff_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -pf_10[k]
                  + f_0 * ff_10[k];

        t_11[k] = -pf_11[k]
                  + f_0 * ff_11[k];

        t_12[k] = -pf_12[k]
                  + f_0 * ff_12[k];

        t_13[k] = -pf_13[k]
                  + f_0 * ff_13[k];

        t_14[k] = -pf_14[k]
                  + f_0 * ff_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pf_15, pf_16, pf_17, pf_18, pf_19, \
                         ff_15, ff_16, ff_17, ff_18, ff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -pf_15[k]
                  + f_0 * ff_15[k];

        t_16[k] = -pf_16[k]
                  + f_0 * ff_16[k];

        t_17[k] = -pf_17[k]
                  + f_0 * ff_17[k];

        t_18[k] = -pf_18[k]
                  + f_0 * ff_18[k];

        t_19[k] = -pf_19[k]
                  + f_0 * ff_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pf_20, pf_21, pf_22, pf_23, pf_24, \
                         ff_20, ff_21, ff_22, ff_23, ff_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -pf_20[k]
                  + f_0 * ff_20[k];

        t_21[k] = -pf_21[k]
                  + f_0 * ff_21[k];

        t_22[k] = -pf_22[k]
                  + f_0 * ff_22[k];

        t_23[k] = -pf_23[k]
                  + f_0 * ff_23[k];

        t_24[k] = -pf_24[k]
                  + f_0 * ff_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pf_25, pf_26, pf_27, pf_28, pf_29, \
                         ff_25, ff_26, ff_27, ff_28, ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -pf_25[k]
                  + f_0 * ff_25[k];

        t_26[k] = -pf_26[k]
                  + f_0 * ff_26[k];

        t_27[k] = -pf_27[k]
                  + f_0 * ff_27[k];

        t_28[k] = -pf_28[k]
                  + f_0 * ff_28[k];

        t_29[k] = -pf_29[k]
                  + f_0 * ff_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, t_36, t_37, ff_30, ff_31, ff_32, \
                         ff_33, ff_34, ff_35, ff_36, ff_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * ff_30[k];

        t_31[k] = f_0 * ff_31[k];

        t_32[k] = f_0 * ff_32[k];

        t_33[k] = f_0 * ff_33[k];

        t_34[k] = f_0 * ff_34[k];

        t_35[k] = f_0 * ff_35[k];

        t_36[k] = f_0 * ff_36[k];

        t_37[k] = f_0 * ff_37[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, t_44, t_45, ff_38, ff_39, ff_40, \
                         ff_41, ff_42, ff_43, ff_44, ff_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * ff_38[k];

        t_39[k] = f_0 * ff_39[k];

        t_40[k] = f_0 * ff_40[k];

        t_41[k] = f_0 * ff_41[k];

        t_42[k] = f_0 * ff_42[k];

        t_43[k] = f_0 * ff_43[k];

        t_44[k] = f_0 * ff_44[k];

        t_45[k] = f_0 * ff_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, t_51, t_52, t_53, ff_46, ff_47, ff_48, \
                         ff_49, ff_50, ff_51, ff_52, ff_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_0 * ff_46[k];

        t_47[k] = f_0 * ff_47[k];

        t_48[k] = f_0 * ff_48[k];

        t_49[k] = f_0 * ff_49[k];

        t_50[k] = f_0 * ff_50[k];

        t_51[k] = f_0 * ff_51[k];

        t_52[k] = f_0 * ff_52[k];

        t_53[k] = f_0 * ff_53[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, ff_54, ff_55, ff_56, ff_57, \
                         ff_58, ff_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_0 * ff_54[k];

        t_55[k] = f_0 * ff_55[k];

        t_56[k] = f_0 * ff_56[k];

        t_57[k] = f_0 * ff_57[k];

        t_58[k] = f_0 * ff_58[k];

        t_59[k] = f_0 * ff_59[k];
    }
}

auto
compute_prim_geom_10_df_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t pf, const size_t ff,
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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_12 = buffer.data(pf + 12);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_14 = buffer.data(pf + 14);
    const auto *pf_15 = buffer.data(pf + 15);
    const auto *pf_16 = buffer.data(pf + 16);
    const auto *pf_17 = buffer.data(pf + 17);
    const auto *pf_18 = buffer.data(pf + 18);
    const auto *pf_19 = buffer.data(pf + 19);
    const auto *pf_20 = buffer.data(pf + 20);
    const auto *pf_21 = buffer.data(pf + 21);
    const auto *pf_22 = buffer.data(pf + 22);
    const auto *pf_23 = buffer.data(pf + 23);
    const auto *pf_24 = buffer.data(pf + 24);
    const auto *pf_25 = buffer.data(pf + 25);
    const auto *pf_26 = buffer.data(pf + 26);
    const auto *pf_27 = buffer.data(pf + 27);
    const auto *pf_28 = buffer.data(pf + 28);
    const auto *pf_29 = buffer.data(pf + 29);

    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_41 = buffer.data(ff + 41);
    const auto *ff_42 = buffer.data(ff + 42);
    const auto *ff_43 = buffer.data(ff + 43);
    const auto *ff_44 = buffer.data(ff + 44);
    const auto *ff_45 = buffer.data(ff + 45);
    const auto *ff_46 = buffer.data(ff + 46);
    const auto *ff_47 = buffer.data(ff + 47);
    const auto *ff_48 = buffer.data(ff + 48);
    const auto *ff_49 = buffer.data(ff + 49);
    const auto *ff_60 = buffer.data(ff + 60);
    const auto *ff_61 = buffer.data(ff + 61);
    const auto *ff_62 = buffer.data(ff + 62);
    const auto *ff_63 = buffer.data(ff + 63);
    const auto *ff_64 = buffer.data(ff + 64);
    const auto *ff_65 = buffer.data(ff + 65);
    const auto *ff_66 = buffer.data(ff + 66);
    const auto *ff_67 = buffer.data(ff + 67);
    const auto *ff_68 = buffer.data(ff + 68);
    const auto *ff_69 = buffer.data(ff + 69);
    const auto *ff_70 = buffer.data(ff + 70);
    const auto *ff_71 = buffer.data(ff + 71);
    const auto *ff_72 = buffer.data(ff + 72);
    const auto *ff_73 = buffer.data(ff + 73);
    const auto *ff_74 = buffer.data(ff + 74);
    const auto *ff_75 = buffer.data(ff + 75);
    const auto *ff_76 = buffer.data(ff + 76);
    const auto *ff_77 = buffer.data(ff + 77);
    const auto *ff_78 = buffer.data(ff + 78);
    const auto *ff_79 = buffer.data(ff + 79);
    const auto *ff_80 = buffer.data(ff + 80);
    const auto *ff_81 = buffer.data(ff + 81);
    const auto *ff_82 = buffer.data(ff + 82);
    const auto *ff_83 = buffer.data(ff + 83);
    const auto *ff_84 = buffer.data(ff + 84);
    const auto *ff_85 = buffer.data(ff + 85);
    const auto *ff_86 = buffer.data(ff + 86);
    const auto *ff_87 = buffer.data(ff + 87);
    const auto *ff_88 = buffer.data(ff + 88);
    const auto *ff_89 = buffer.data(ff + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, ff_10, ff_11, ff_12, ff_13, \
                         ff_14, ff_15, ff_16, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_10[k];

        t_1[k] = f_0 * ff_11[k];

        t_2[k] = f_0 * ff_12[k];

        t_3[k] = f_0 * ff_13[k];

        t_4[k] = f_0 * ff_14[k];

        t_5[k] = f_0 * ff_15[k];

        t_6[k] = f_0 * ff_16[k];

        t_7[k] = f_0 * ff_17[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, pf_0, pf_1, pf_2, pf_3, ff_18, \
                         ff_19, ff_30, ff_31, ff_32, ff_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * ff_18[k];

        t_9[k] = f_0 * ff_19[k];

        t_10[k] = -pf_0[k]
                  + f_0 * ff_30[k];

        t_11[k] = -pf_1[k]
                  + f_0 * ff_31[k];

        t_12[k] = -pf_2[k]
                  + f_0 * ff_32[k];

        t_13[k] = -pf_3[k]
                  + f_0 * ff_33[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, pf_4, pf_5, pf_6, pf_7, pf_8, ff_34, \
                         ff_35, ff_36, ff_37, ff_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -pf_4[k]
                  + f_0 * ff_34[k];

        t_15[k] = -pf_5[k]
                  + f_0 * ff_35[k];

        t_16[k] = -pf_6[k]
                  + f_0 * ff_36[k];

        t_17[k] = -pf_7[k]
                  + f_0 * ff_37[k];

        t_18[k] = -pf_8[k]
                  + f_0 * ff_38[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, t_25, pf_9, ff_39, ff_40, ff_41, \
                         ff_42, ff_43, ff_44, ff_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -pf_9[k]
                  + f_0 * ff_39[k];

        t_20[k] = f_0 * ff_40[k];

        t_21[k] = f_0 * ff_41[k];

        t_22[k] = f_0 * ff_42[k];

        t_23[k] = f_0 * ff_43[k];

        t_24[k] = f_0 * ff_44[k];

        t_25[k] = f_0 * ff_45[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, pf_10, pf_11, ff_46, ff_47, \
                         ff_48, ff_49, ff_60, ff_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * ff_46[k];

        t_27[k] = f_0 * ff_47[k];

        t_28[k] = f_0 * ff_48[k];

        t_29[k] = f_0 * ff_49[k];

        t_30[k] = -2.0 * pf_10[k]
                  + f_0 * ff_60[k];

        t_31[k] = -2.0 * pf_11[k]
                  + f_0 * ff_61[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, pf_12, pf_13, pf_14, pf_15, pf_16, \
                         ff_62, ff_63, ff_64, ff_65, ff_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -2.0 * pf_12[k]
                  + f_0 * ff_62[k];

        t_33[k] = -2.0 * pf_13[k]
                  + f_0 * ff_63[k];

        t_34[k] = -2.0 * pf_14[k]
                  + f_0 * ff_64[k];

        t_35[k] = -2.0 * pf_15[k]
                  + f_0 * ff_65[k];

        t_36[k] = -2.0 * pf_16[k]
                  + f_0 * ff_66[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pf_17, pf_18, pf_19, pf_20, pf_21, \
                         ff_67, ff_68, ff_69, ff_70, ff_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -2.0 * pf_17[k]
                  + f_0 * ff_67[k];

        t_38[k] = -2.0 * pf_18[k]
                  + f_0 * ff_68[k];

        t_39[k] = -2.0 * pf_19[k]
                  + f_0 * ff_69[k];

        t_40[k] = -pf_20[k]
                  + f_0 * ff_70[k];

        t_41[k] = -pf_21[k]
                  + f_0 * ff_71[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pf_22, pf_23, pf_24, pf_25, pf_26, \
                         ff_72, ff_73, ff_74, ff_75, ff_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = -pf_22[k]
                  + f_0 * ff_72[k];

        t_43[k] = -pf_23[k]
                  + f_0 * ff_73[k];

        t_44[k] = -pf_24[k]
                  + f_0 * ff_74[k];

        t_45[k] = -pf_25[k]
                  + f_0 * ff_75[k];

        t_46[k] = -pf_26[k]
                  + f_0 * ff_76[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, pf_27, pf_28, pf_29, ff_77, \
                         ff_78, ff_79, ff_80, ff_81, ff_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -pf_27[k]
                  + f_0 * ff_77[k];

        t_48[k] = -pf_28[k]
                  + f_0 * ff_78[k];

        t_49[k] = -pf_29[k]
                  + f_0 * ff_79[k];

        t_50[k] = f_0 * ff_80[k];

        t_51[k] = f_0 * ff_81[k];

        t_52[k] = f_0 * ff_82[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, t_58, t_59, ff_83, ff_84, ff_85, ff_86, \
                         ff_87, ff_88, ff_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * ff_83[k];

        t_54[k] = f_0 * ff_84[k];

        t_55[k] = f_0 * ff_85[k];

        t_56[k] = f_0 * ff_86[k];

        t_57[k] = f_0 * ff_87[k];

        t_58[k] = f_0 * ff_88[k];

        t_59[k] = f_0 * ff_89[k];
    }
}

auto
compute_prim_geom_10_df_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t pf, const size_t ff,
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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_12 = buffer.data(pf + 12);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_14 = buffer.data(pf + 14);
    const auto *pf_15 = buffer.data(pf + 15);
    const auto *pf_16 = buffer.data(pf + 16);
    const auto *pf_17 = buffer.data(pf + 17);
    const auto *pf_18 = buffer.data(pf + 18);
    const auto *pf_19 = buffer.data(pf + 19);
    const auto *pf_20 = buffer.data(pf + 20);
    const auto *pf_21 = buffer.data(pf + 21);
    const auto *pf_22 = buffer.data(pf + 22);
    const auto *pf_23 = buffer.data(pf + 23);
    const auto *pf_24 = buffer.data(pf + 24);
    const auto *pf_25 = buffer.data(pf + 25);
    const auto *pf_26 = buffer.data(pf + 26);
    const auto *pf_27 = buffer.data(pf + 27);
    const auto *pf_28 = buffer.data(pf + 28);
    const auto *pf_29 = buffer.data(pf + 29);

    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_41 = buffer.data(ff + 41);
    const auto *ff_42 = buffer.data(ff + 42);
    const auto *ff_43 = buffer.data(ff + 43);
    const auto *ff_44 = buffer.data(ff + 44);
    const auto *ff_45 = buffer.data(ff + 45);
    const auto *ff_46 = buffer.data(ff + 46);
    const auto *ff_47 = buffer.data(ff + 47);
    const auto *ff_48 = buffer.data(ff + 48);
    const auto *ff_49 = buffer.data(ff + 49);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_51 = buffer.data(ff + 51);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_53 = buffer.data(ff + 53);
    const auto *ff_54 = buffer.data(ff + 54);
    const auto *ff_55 = buffer.data(ff + 55);
    const auto *ff_56 = buffer.data(ff + 56);
    const auto *ff_57 = buffer.data(ff + 57);
    const auto *ff_58 = buffer.data(ff + 58);
    const auto *ff_59 = buffer.data(ff + 59);
    const auto *ff_70 = buffer.data(ff + 70);
    const auto *ff_71 = buffer.data(ff + 71);
    const auto *ff_72 = buffer.data(ff + 72);
    const auto *ff_73 = buffer.data(ff + 73);
    const auto *ff_74 = buffer.data(ff + 74);
    const auto *ff_75 = buffer.data(ff + 75);
    const auto *ff_76 = buffer.data(ff + 76);
    const auto *ff_77 = buffer.data(ff + 77);
    const auto *ff_78 = buffer.data(ff + 78);
    const auto *ff_79 = buffer.data(ff + 79);
    const auto *ff_80 = buffer.data(ff + 80);
    const auto *ff_81 = buffer.data(ff + 81);
    const auto *ff_82 = buffer.data(ff + 82);
    const auto *ff_83 = buffer.data(ff + 83);
    const auto *ff_84 = buffer.data(ff + 84);
    const auto *ff_85 = buffer.data(ff + 85);
    const auto *ff_86 = buffer.data(ff + 86);
    const auto *ff_87 = buffer.data(ff + 87);
    const auto *ff_88 = buffer.data(ff + 88);
    const auto *ff_89 = buffer.data(ff + 89);
    const auto *ff_90 = buffer.data(ff + 90);
    const auto *ff_91 = buffer.data(ff + 91);
    const auto *ff_92 = buffer.data(ff + 92);
    const auto *ff_93 = buffer.data(ff + 93);
    const auto *ff_94 = buffer.data(ff + 94);
    const auto *ff_95 = buffer.data(ff + 95);
    const auto *ff_96 = buffer.data(ff + 96);
    const auto *ff_97 = buffer.data(ff + 97);
    const auto *ff_98 = buffer.data(ff + 98);
    const auto *ff_99 = buffer.data(ff + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, ff_20, ff_21, ff_22, ff_23, \
                         ff_24, ff_25, ff_26, ff_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_20[k];

        t_1[k] = f_0 * ff_21[k];

        t_2[k] = f_0 * ff_22[k];

        t_3[k] = f_0 * ff_23[k];

        t_4[k] = f_0 * ff_24[k];

        t_5[k] = f_0 * ff_25[k];

        t_6[k] = f_0 * ff_26[k];

        t_7[k] = f_0 * ff_27[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, ff_28, ff_29, ff_40, \
                         ff_41, ff_42, ff_43, ff_44, ff_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * ff_28[k];

        t_9[k] = f_0 * ff_29[k];

        t_10[k] = f_0 * ff_40[k];

        t_11[k] = f_0 * ff_41[k];

        t_12[k] = f_0 * ff_42[k];

        t_13[k] = f_0 * ff_43[k];

        t_14[k] = f_0 * ff_44[k];

        t_15[k] = f_0 * ff_45[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, pf_0, pf_1, ff_46, ff_47, ff_48, \
                         ff_49, ff_50, ff_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * ff_46[k];

        t_17[k] = f_0 * ff_47[k];

        t_18[k] = f_0 * ff_48[k];

        t_19[k] = f_0 * ff_49[k];

        t_20[k] = -pf_0[k]
                  + f_0 * ff_50[k];

        t_21[k] = -pf_1[k]
                  + f_0 * ff_51[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pf_2, pf_3, pf_4, pf_5, pf_6, ff_52, \
                         ff_53, ff_54, ff_55, ff_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -pf_2[k]
                  + f_0 * ff_52[k];

        t_23[k] = -pf_3[k]
                  + f_0 * ff_53[k];

        t_24[k] = -pf_4[k]
                  + f_0 * ff_54[k];

        t_25[k] = -pf_5[k]
                  + f_0 * ff_55[k];

        t_26[k] = -pf_6[k]
                  + f_0 * ff_56[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, t_32, pf_7, pf_8, pf_9, ff_57, ff_58, \
                         ff_59, ff_70, ff_71, ff_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -pf_7[k]
                  + f_0 * ff_57[k];

        t_28[k] = -pf_8[k]
                  + f_0 * ff_58[k];

        t_29[k] = -pf_9[k]
                  + f_0 * ff_59[k];

        t_30[k] = f_0 * ff_70[k];

        t_31[k] = f_0 * ff_71[k];

        t_32[k] = f_0 * ff_72[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, t_39, ff_73, ff_74, ff_75, ff_76, \
                         ff_77, ff_78, ff_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * ff_73[k];

        t_34[k] = f_0 * ff_74[k];

        t_35[k] = f_0 * ff_75[k];

        t_36[k] = f_0 * ff_76[k];

        t_37[k] = f_0 * ff_77[k];

        t_38[k] = f_0 * ff_78[k];

        t_39[k] = f_0 * ff_79[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pf_10, pf_11, pf_12, pf_13, pf_14, \
                         ff_80, ff_81, ff_82, ff_83, ff_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -pf_10[k]
                  + f_0 * ff_80[k];

        t_41[k] = -pf_11[k]
                  + f_0 * ff_81[k];

        t_42[k] = -pf_12[k]
                  + f_0 * ff_82[k];

        t_43[k] = -pf_13[k]
                  + f_0 * ff_83[k];

        t_44[k] = -pf_14[k]
                  + f_0 * ff_84[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pf_15, pf_16, pf_17, pf_18, pf_19, \
                         ff_85, ff_86, ff_87, ff_88, ff_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -pf_15[k]
                  + f_0 * ff_85[k];

        t_46[k] = -pf_16[k]
                  + f_0 * ff_86[k];

        t_47[k] = -pf_17[k]
                  + f_0 * ff_87[k];

        t_48[k] = -pf_18[k]
                  + f_0 * ff_88[k];

        t_49[k] = -pf_19[k]
                  + f_0 * ff_89[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pf_20, pf_21, pf_22, pf_23, pf_24, \
                         ff_90, ff_91, ff_92, ff_93, ff_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * pf_20[k]
                  + f_0 * ff_90[k];

        t_51[k] = -2.0 * pf_21[k]
                  + f_0 * ff_91[k];

        t_52[k] = -2.0 * pf_22[k]
                  + f_0 * ff_92[k];

        t_53[k] = -2.0 * pf_23[k]
                  + f_0 * ff_93[k];

        t_54[k] = -2.0 * pf_24[k]
                  + f_0 * ff_94[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pf_25, pf_26, pf_27, pf_28, pf_29, \
                         ff_95, ff_96, ff_97, ff_98, ff_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * pf_25[k]
                  + f_0 * ff_95[k];

        t_56[k] = -2.0 * pf_26[k]
                  + f_0 * ff_96[k];

        t_57[k] = -2.0 * pf_27[k]
                  + f_0 * ff_97[k];

        t_58[k] = -2.0 * pf_28[k]
                  + f_0 * ff_98[k];

        t_59[k] = -2.0 * pf_29[k]
                  + f_0 * ff_99[k];
    }
}

}  // namespace simdt2ceri
