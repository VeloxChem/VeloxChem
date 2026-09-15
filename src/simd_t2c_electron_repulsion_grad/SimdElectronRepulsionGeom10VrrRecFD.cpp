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


#include "SimdElectronRepulsionGeom10VrrRecFD.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_fd_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t dd, const size_t gd,
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

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_15 = buffer.data(dd + 15);
    const auto *dd_16 = buffer.data(dd + 16);
    const auto *dd_17 = buffer.data(dd + 17);
    const auto *dd_18 = buffer.data(dd + 18);
    const auto *dd_19 = buffer.data(dd + 19);
    const auto *dd_20 = buffer.data(dd + 20);
    const auto *dd_21 = buffer.data(dd + 21);
    const auto *dd_22 = buffer.data(dd + 22);
    const auto *dd_23 = buffer.data(dd + 23);
    const auto *dd_24 = buffer.data(dd + 24);
    const auto *dd_25 = buffer.data(dd + 25);
    const auto *dd_26 = buffer.data(dd + 26);
    const auto *dd_27 = buffer.data(dd + 27);
    const auto *dd_28 = buffer.data(dd + 28);
    const auto *dd_29 = buffer.data(dd + 29);
    const auto *dd_30 = buffer.data(dd + 30);
    const auto *dd_31 = buffer.data(dd + 31);
    const auto *dd_32 = buffer.data(dd + 32);
    const auto *dd_33 = buffer.data(dd + 33);
    const auto *dd_34 = buffer.data(dd + 34);
    const auto *dd_35 = buffer.data(dd + 35);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_40 = buffer.data(gd + 40);
    const auto *gd_41 = buffer.data(gd + 41);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_43 = buffer.data(gd + 43);
    const auto *gd_44 = buffer.data(gd + 44);
    const auto *gd_45 = buffer.data(gd + 45);
    const auto *gd_46 = buffer.data(gd + 46);
    const auto *gd_47 = buffer.data(gd + 47);
    const auto *gd_48 = buffer.data(gd + 48);
    const auto *gd_49 = buffer.data(gd + 49);
    const auto *gd_50 = buffer.data(gd + 50);
    const auto *gd_51 = buffer.data(gd + 51);
    const auto *gd_52 = buffer.data(gd + 52);
    const auto *gd_53 = buffer.data(gd + 53);
    const auto *gd_54 = buffer.data(gd + 54);
    const auto *gd_55 = buffer.data(gd + 55);
    const auto *gd_56 = buffer.data(gd + 56);
    const auto *gd_57 = buffer.data(gd + 57);
    const auto *gd_58 = buffer.data(gd + 58);
    const auto *gd_59 = buffer.data(gd + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, dd_0, dd_1, dd_2, dd_3, dd_4, gd_0, gd_1, \
                         gd_2, gd_3, gd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -3.0 * dd_0[k]
                 + f_0 * gd_0[k];

        t_1[k] = -3.0 * dd_1[k]
                 + f_0 * gd_1[k];

        t_2[k] = -3.0 * dd_2[k]
                 + f_0 * gd_2[k];

        t_3[k] = -3.0 * dd_3[k]
                 + f_0 * gd_3[k];

        t_4[k] = -3.0 * dd_4[k]
                 + f_0 * gd_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, dd_5, dd_6, dd_7, dd_8, dd_9, gd_5, gd_6, \
                         gd_7, gd_8, gd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -3.0 * dd_5[k]
                 + f_0 * gd_5[k];

        t_6[k] = -2.0 * dd_6[k]
                 + f_0 * gd_6[k];

        t_7[k] = -2.0 * dd_7[k]
                 + f_0 * gd_7[k];

        t_8[k] = -2.0 * dd_8[k]
                 + f_0 * gd_8[k];

        t_9[k] = -2.0 * dd_9[k]
                 + f_0 * gd_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, dd_10, dd_11, dd_12, dd_13, dd_14, \
                         gd_10, gd_11, gd_12, gd_13, gd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -2.0 * dd_10[k]
                  + f_0 * gd_10[k];

        t_11[k] = -2.0 * dd_11[k]
                  + f_0 * gd_11[k];

        t_12[k] = -2.0 * dd_12[k]
                  + f_0 * gd_12[k];

        t_13[k] = -2.0 * dd_13[k]
                  + f_0 * gd_13[k];

        t_14[k] = -2.0 * dd_14[k]
                  + f_0 * gd_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, dd_15, dd_16, dd_17, dd_18, dd_19, \
                         gd_15, gd_16, gd_17, gd_18, gd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -2.0 * dd_15[k]
                  + f_0 * gd_15[k];

        t_16[k] = -2.0 * dd_16[k]
                  + f_0 * gd_16[k];

        t_17[k] = -2.0 * dd_17[k]
                  + f_0 * gd_17[k];

        t_18[k] = -dd_18[k]
                  + f_0 * gd_18[k];

        t_19[k] = -dd_19[k]
                  + f_0 * gd_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, dd_20, dd_21, dd_22, dd_23, dd_24, \
                         gd_20, gd_21, gd_22, gd_23, gd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -dd_20[k]
                  + f_0 * gd_20[k];

        t_21[k] = -dd_21[k]
                  + f_0 * gd_21[k];

        t_22[k] = -dd_22[k]
                  + f_0 * gd_22[k];

        t_23[k] = -dd_23[k]
                  + f_0 * gd_23[k];

        t_24[k] = -dd_24[k]
                  + f_0 * gd_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, dd_25, dd_26, dd_27, dd_28, dd_29, \
                         gd_25, gd_26, gd_27, gd_28, gd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -dd_25[k]
                  + f_0 * gd_25[k];

        t_26[k] = -dd_26[k]
                  + f_0 * gd_26[k];

        t_27[k] = -dd_27[k]
                  + f_0 * gd_27[k];

        t_28[k] = -dd_28[k]
                  + f_0 * gd_28[k];

        t_29[k] = -dd_29[k]
                  + f_0 * gd_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, dd_30, dd_31, dd_32, dd_33, dd_34, \
                         gd_30, gd_31, gd_32, gd_33, gd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -dd_30[k]
                  + f_0 * gd_30[k];

        t_31[k] = -dd_31[k]
                  + f_0 * gd_31[k];

        t_32[k] = -dd_32[k]
                  + f_0 * gd_32[k];

        t_33[k] = -dd_33[k]
                  + f_0 * gd_33[k];

        t_34[k] = -dd_34[k]
                  + f_0 * gd_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, t_40, t_41, dd_35, gd_35, gd_36, gd_37, \
                         gd_38, gd_39, gd_40, gd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -dd_35[k]
                  + f_0 * gd_35[k];

        t_36[k] = f_0 * gd_36[k];

        t_37[k] = f_0 * gd_37[k];

        t_38[k] = f_0 * gd_38[k];

        t_39[k] = f_0 * gd_39[k];

        t_40[k] = f_0 * gd_40[k];

        t_41[k] = f_0 * gd_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, t_48, t_49, gd_42, gd_43, gd_44, \
                         gd_45, gd_46, gd_47, gd_48, gd_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_0 * gd_42[k];

        t_43[k] = f_0 * gd_43[k];

        t_44[k] = f_0 * gd_44[k];

        t_45[k] = f_0 * gd_45[k];

        t_46[k] = f_0 * gd_46[k];

        t_47[k] = f_0 * gd_47[k];

        t_48[k] = f_0 * gd_48[k];

        t_49[k] = f_0 * gd_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, t_55, t_56, t_57, gd_50, gd_51, gd_52, \
                         gd_53, gd_54, gd_55, gd_56, gd_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_0 * gd_50[k];

        t_51[k] = f_0 * gd_51[k];

        t_52[k] = f_0 * gd_52[k];

        t_53[k] = f_0 * gd_53[k];

        t_54[k] = f_0 * gd_54[k];

        t_55[k] = f_0 * gd_55[k];

        t_56[k] = f_0 * gd_56[k];

        t_57[k] = f_0 * gd_57[k];
    }

#pragma omp simd aligned(t_58, t_59, gd_58, gd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_0 * gd_58[k];

        t_59[k] = f_0 * gd_59[k];
    }
}

auto
compute_prim_geom_10_fd_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t dd, const size_t gd,
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

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_15 = buffer.data(dd + 15);
    const auto *dd_16 = buffer.data(dd + 16);
    const auto *dd_17 = buffer.data(dd + 17);
    const auto *dd_18 = buffer.data(dd + 18);
    const auto *dd_19 = buffer.data(dd + 19);
    const auto *dd_20 = buffer.data(dd + 20);
    const auto *dd_21 = buffer.data(dd + 21);
    const auto *dd_22 = buffer.data(dd + 22);
    const auto *dd_23 = buffer.data(dd + 23);
    const auto *dd_24 = buffer.data(dd + 24);
    const auto *dd_25 = buffer.data(dd + 25);
    const auto *dd_26 = buffer.data(dd + 26);
    const auto *dd_27 = buffer.data(dd + 27);
    const auto *dd_28 = buffer.data(dd + 28);
    const auto *dd_29 = buffer.data(dd + 29);
    const auto *dd_30 = buffer.data(dd + 30);
    const auto *dd_31 = buffer.data(dd + 31);
    const auto *dd_32 = buffer.data(dd + 32);
    const auto *dd_33 = buffer.data(dd + 33);
    const auto *dd_34 = buffer.data(dd + 34);
    const auto *dd_35 = buffer.data(dd + 35);

    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_40 = buffer.data(gd + 40);
    const auto *gd_41 = buffer.data(gd + 41);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_43 = buffer.data(gd + 43);
    const auto *gd_44 = buffer.data(gd + 44);
    const auto *gd_45 = buffer.data(gd + 45);
    const auto *gd_46 = buffer.data(gd + 46);
    const auto *gd_47 = buffer.data(gd + 47);
    const auto *gd_48 = buffer.data(gd + 48);
    const auto *gd_49 = buffer.data(gd + 49);
    const auto *gd_50 = buffer.data(gd + 50);
    const auto *gd_51 = buffer.data(gd + 51);
    const auto *gd_52 = buffer.data(gd + 52);
    const auto *gd_53 = buffer.data(gd + 53);
    const auto *gd_60 = buffer.data(gd + 60);
    const auto *gd_61 = buffer.data(gd + 61);
    const auto *gd_62 = buffer.data(gd + 62);
    const auto *gd_63 = buffer.data(gd + 63);
    const auto *gd_64 = buffer.data(gd + 64);
    const auto *gd_65 = buffer.data(gd + 65);
    const auto *gd_66 = buffer.data(gd + 66);
    const auto *gd_67 = buffer.data(gd + 67);
    const auto *gd_68 = buffer.data(gd + 68);
    const auto *gd_69 = buffer.data(gd + 69);
    const auto *gd_70 = buffer.data(gd + 70);
    const auto *gd_71 = buffer.data(gd + 71);
    const auto *gd_72 = buffer.data(gd + 72);
    const auto *gd_73 = buffer.data(gd + 73);
    const auto *gd_74 = buffer.data(gd + 74);
    const auto *gd_75 = buffer.data(gd + 75);
    const auto *gd_76 = buffer.data(gd + 76);
    const auto *gd_77 = buffer.data(gd + 77);
    const auto *gd_78 = buffer.data(gd + 78);
    const auto *gd_79 = buffer.data(gd + 79);
    const auto *gd_80 = buffer.data(gd + 80);
    const auto *gd_81 = buffer.data(gd + 81);
    const auto *gd_82 = buffer.data(gd + 82);
    const auto *gd_83 = buffer.data(gd + 83);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, dd_0, gd_6, gd_7, gd_8, gd_9, \
                         gd_10, gd_11, gd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_6[k];

        t_1[k] = f_0 * gd_7[k];

        t_2[k] = f_0 * gd_8[k];

        t_3[k] = f_0 * gd_9[k];

        t_4[k] = f_0 * gd_10[k];

        t_5[k] = f_0 * gd_11[k];

        t_6[k] = -dd_0[k]
                 + f_0 * gd_18[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, dd_1, dd_2, dd_3, dd_4, dd_5, gd_19, \
                         gd_20, gd_21, gd_22, gd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -dd_1[k]
                 + f_0 * gd_19[k];

        t_8[k] = -dd_2[k]
                 + f_0 * gd_20[k];

        t_9[k] = -dd_3[k]
                 + f_0 * gd_21[k];

        t_10[k] = -dd_4[k]
                  + f_0 * gd_22[k];

        t_11[k] = -dd_5[k]
                  + f_0 * gd_23[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, t_18, dd_6, gd_24, gd_25, gd_26, \
                         gd_27, gd_28, gd_29, gd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * gd_24[k];

        t_13[k] = f_0 * gd_25[k];

        t_14[k] = f_0 * gd_26[k];

        t_15[k] = f_0 * gd_27[k];

        t_16[k] = f_0 * gd_28[k];

        t_17[k] = f_0 * gd_29[k];

        t_18[k] = -2.0 * dd_6[k]
                  + f_0 * gd_36[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, dd_7, dd_8, dd_9, dd_10, dd_11, gd_37, \
                         gd_38, gd_39, gd_40, gd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -2.0 * dd_7[k]
                  + f_0 * gd_37[k];

        t_20[k] = -2.0 * dd_8[k]
                  + f_0 * gd_38[k];

        t_21[k] = -2.0 * dd_9[k]
                  + f_0 * gd_39[k];

        t_22[k] = -2.0 * dd_10[k]
                  + f_0 * gd_40[k];

        t_23[k] = -2.0 * dd_11[k]
                  + f_0 * gd_41[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, dd_12, dd_13, dd_14, dd_15, dd_16, \
                         gd_42, gd_43, gd_44, gd_45, gd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -dd_12[k]
                  + f_0 * gd_42[k];

        t_25[k] = -dd_13[k]
                  + f_0 * gd_43[k];

        t_26[k] = -dd_14[k]
                  + f_0 * gd_44[k];

        t_27[k] = -dd_15[k]
                  + f_0 * gd_45[k];

        t_28[k] = -dd_16[k]
                  + f_0 * gd_46[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, t_35, dd_17, gd_47, gd_48, gd_49, \
                         gd_50, gd_51, gd_52, gd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -dd_17[k]
                  + f_0 * gd_47[k];

        t_30[k] = f_0 * gd_48[k];

        t_31[k] = f_0 * gd_49[k];

        t_32[k] = f_0 * gd_50[k];

        t_33[k] = f_0 * gd_51[k];

        t_34[k] = f_0 * gd_52[k];

        t_35[k] = f_0 * gd_53[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, dd_18, dd_19, dd_20, dd_21, dd_22, \
                         gd_60, gd_61, gd_62, gd_63, gd_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -3.0 * dd_18[k]
                  + f_0 * gd_60[k];

        t_37[k] = -3.0 * dd_19[k]
                  + f_0 * gd_61[k];

        t_38[k] = -3.0 * dd_20[k]
                  + f_0 * gd_62[k];

        t_39[k] = -3.0 * dd_21[k]
                  + f_0 * gd_63[k];

        t_40[k] = -3.0 * dd_22[k]
                  + f_0 * gd_64[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, dd_23, dd_24, dd_25, dd_26, dd_27, \
                         gd_65, gd_66, gd_67, gd_68, gd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -3.0 * dd_23[k]
                  + f_0 * gd_65[k];

        t_42[k] = -2.0 * dd_24[k]
                  + f_0 * gd_66[k];

        t_43[k] = -2.0 * dd_25[k]
                  + f_0 * gd_67[k];

        t_44[k] = -2.0 * dd_26[k]
                  + f_0 * gd_68[k];

        t_45[k] = -2.0 * dd_27[k]
                  + f_0 * gd_69[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, dd_28, dd_29, dd_30, dd_31, dd_32, \
                         gd_70, gd_71, gd_72, gd_73, gd_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -2.0 * dd_28[k]
                  + f_0 * gd_70[k];

        t_47[k] = -2.0 * dd_29[k]
                  + f_0 * gd_71[k];

        t_48[k] = -dd_30[k]
                  + f_0 * gd_72[k];

        t_49[k] = -dd_31[k]
                  + f_0 * gd_73[k];

        t_50[k] = -dd_32[k]
                  + f_0 * gd_74[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, dd_33, dd_34, dd_35, gd_75, \
                         gd_76, gd_77, gd_78, gd_79, gd_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -dd_33[k]
                  + f_0 * gd_75[k];

        t_52[k] = -dd_34[k]
                  + f_0 * gd_76[k];

        t_53[k] = -dd_35[k]
                  + f_0 * gd_77[k];

        t_54[k] = f_0 * gd_78[k];

        t_55[k] = f_0 * gd_79[k];

        t_56[k] = f_0 * gd_80[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, gd_81, gd_82, gd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_0 * gd_81[k];

        t_58[k] = f_0 * gd_82[k];

        t_59[k] = f_0 * gd_83[k];
    }
}

auto
compute_prim_geom_10_fd_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t dd, const size_t gd,
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

    const auto *dd_0 = buffer.data(dd + 0);
    const auto *dd_1 = buffer.data(dd + 1);
    const auto *dd_2 = buffer.data(dd + 2);
    const auto *dd_3 = buffer.data(dd + 3);
    const auto *dd_4 = buffer.data(dd + 4);
    const auto *dd_5 = buffer.data(dd + 5);
    const auto *dd_6 = buffer.data(dd + 6);
    const auto *dd_7 = buffer.data(dd + 7);
    const auto *dd_8 = buffer.data(dd + 8);
    const auto *dd_9 = buffer.data(dd + 9);
    const auto *dd_10 = buffer.data(dd + 10);
    const auto *dd_11 = buffer.data(dd + 11);
    const auto *dd_12 = buffer.data(dd + 12);
    const auto *dd_13 = buffer.data(dd + 13);
    const auto *dd_14 = buffer.data(dd + 14);
    const auto *dd_15 = buffer.data(dd + 15);
    const auto *dd_16 = buffer.data(dd + 16);
    const auto *dd_17 = buffer.data(dd + 17);
    const auto *dd_18 = buffer.data(dd + 18);
    const auto *dd_19 = buffer.data(dd + 19);
    const auto *dd_20 = buffer.data(dd + 20);
    const auto *dd_21 = buffer.data(dd + 21);
    const auto *dd_22 = buffer.data(dd + 22);
    const auto *dd_23 = buffer.data(dd + 23);
    const auto *dd_24 = buffer.data(dd + 24);
    const auto *dd_25 = buffer.data(dd + 25);
    const auto *dd_26 = buffer.data(dd + 26);
    const auto *dd_27 = buffer.data(dd + 27);
    const auto *dd_28 = buffer.data(dd + 28);
    const auto *dd_29 = buffer.data(dd + 29);
    const auto *dd_30 = buffer.data(dd + 30);
    const auto *dd_31 = buffer.data(dd + 31);
    const auto *dd_32 = buffer.data(dd + 32);
    const auto *dd_33 = buffer.data(dd + 33);
    const auto *dd_34 = buffer.data(dd + 34);
    const auto *dd_35 = buffer.data(dd + 35);

    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_43 = buffer.data(gd + 43);
    const auto *gd_44 = buffer.data(gd + 44);
    const auto *gd_45 = buffer.data(gd + 45);
    const auto *gd_46 = buffer.data(gd + 46);
    const auto *gd_47 = buffer.data(gd + 47);
    const auto *gd_48 = buffer.data(gd + 48);
    const auto *gd_49 = buffer.data(gd + 49);
    const auto *gd_50 = buffer.data(gd + 50);
    const auto *gd_51 = buffer.data(gd + 51);
    const auto *gd_52 = buffer.data(gd + 52);
    const auto *gd_53 = buffer.data(gd + 53);
    const auto *gd_54 = buffer.data(gd + 54);
    const auto *gd_55 = buffer.data(gd + 55);
    const auto *gd_56 = buffer.data(gd + 56);
    const auto *gd_57 = buffer.data(gd + 57);
    const auto *gd_58 = buffer.data(gd + 58);
    const auto *gd_59 = buffer.data(gd + 59);
    const auto *gd_66 = buffer.data(gd + 66);
    const auto *gd_67 = buffer.data(gd + 67);
    const auto *gd_68 = buffer.data(gd + 68);
    const auto *gd_69 = buffer.data(gd + 69);
    const auto *gd_70 = buffer.data(gd + 70);
    const auto *gd_71 = buffer.data(gd + 71);
    const auto *gd_72 = buffer.data(gd + 72);
    const auto *gd_73 = buffer.data(gd + 73);
    const auto *gd_74 = buffer.data(gd + 74);
    const auto *gd_75 = buffer.data(gd + 75);
    const auto *gd_76 = buffer.data(gd + 76);
    const auto *gd_77 = buffer.data(gd + 77);
    const auto *gd_78 = buffer.data(gd + 78);
    const auto *gd_79 = buffer.data(gd + 79);
    const auto *gd_80 = buffer.data(gd + 80);
    const auto *gd_81 = buffer.data(gd + 81);
    const auto *gd_82 = buffer.data(gd + 82);
    const auto *gd_83 = buffer.data(gd + 83);
    const auto *gd_84 = buffer.data(gd + 84);
    const auto *gd_85 = buffer.data(gd + 85);
    const auto *gd_86 = buffer.data(gd + 86);
    const auto *gd_87 = buffer.data(gd + 87);
    const auto *gd_88 = buffer.data(gd + 88);
    const auto *gd_89 = buffer.data(gd + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, gd_12, gd_13, gd_14, gd_15, \
                         gd_16, gd_17, gd_24, gd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_12[k];

        t_1[k] = f_0 * gd_13[k];

        t_2[k] = f_0 * gd_14[k];

        t_3[k] = f_0 * gd_15[k];

        t_4[k] = f_0 * gd_16[k];

        t_5[k] = f_0 * gd_17[k];

        t_6[k] = f_0 * gd_24[k];

        t_7[k] = f_0 * gd_25[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, dd_0, dd_1, gd_26, gd_27, gd_28, \
                         gd_29, gd_30, gd_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * gd_26[k];

        t_9[k] = f_0 * gd_27[k];

        t_10[k] = f_0 * gd_28[k];

        t_11[k] = f_0 * gd_29[k];

        t_12[k] = -dd_0[k]
                  + f_0 * gd_30[k];

        t_13[k] = -dd_1[k]
                  + f_0 * gd_31[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, t_19, dd_2, dd_3, dd_4, dd_5, gd_32, \
                         gd_33, gd_34, gd_35, gd_42, gd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -dd_2[k]
                  + f_0 * gd_32[k];

        t_15[k] = -dd_3[k]
                  + f_0 * gd_33[k];

        t_16[k] = -dd_4[k]
                  + f_0 * gd_34[k];

        t_17[k] = -dd_5[k]
                  + f_0 * gd_35[k];

        t_18[k] = f_0 * gd_42[k];

        t_19[k] = f_0 * gd_43[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, dd_6, dd_7, gd_44, gd_45, gd_46, \
                         gd_47, gd_48, gd_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * gd_44[k];

        t_21[k] = f_0 * gd_45[k];

        t_22[k] = f_0 * gd_46[k];

        t_23[k] = f_0 * gd_47[k];

        t_24[k] = -dd_6[k]
                  + f_0 * gd_48[k];

        t_25[k] = -dd_7[k]
                  + f_0 * gd_49[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, dd_8, dd_9, dd_10, dd_11, dd_12, gd_50, \
                         gd_51, gd_52, gd_53, gd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -dd_8[k]
                  + f_0 * gd_50[k];

        t_27[k] = -dd_9[k]
                  + f_0 * gd_51[k];

        t_28[k] = -dd_10[k]
                  + f_0 * gd_52[k];

        t_29[k] = -dd_11[k]
                  + f_0 * gd_53[k];

        t_30[k] = -2.0 * dd_12[k]
                  + f_0 * gd_54[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, dd_13, dd_14, dd_15, dd_16, dd_17, \
                         gd_55, gd_56, gd_57, gd_58, gd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -2.0 * dd_13[k]
                  + f_0 * gd_55[k];

        t_32[k] = -2.0 * dd_14[k]
                  + f_0 * gd_56[k];

        t_33[k] = -2.0 * dd_15[k]
                  + f_0 * gd_57[k];

        t_34[k] = -2.0 * dd_16[k]
                  + f_0 * gd_58[k];

        t_35[k] = -2.0 * dd_17[k]
                  + f_0 * gd_59[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, t_42, dd_18, gd_66, gd_67, gd_68, \
                         gd_69, gd_70, gd_71, gd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * gd_66[k];

        t_37[k] = f_0 * gd_67[k];

        t_38[k] = f_0 * gd_68[k];

        t_39[k] = f_0 * gd_69[k];

        t_40[k] = f_0 * gd_70[k];

        t_41[k] = f_0 * gd_71[k];

        t_42[k] = -dd_18[k]
                  + f_0 * gd_72[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, dd_19, dd_20, dd_21, dd_22, dd_23, \
                         gd_73, gd_74, gd_75, gd_76, gd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -dd_19[k]
                  + f_0 * gd_73[k];

        t_44[k] = -dd_20[k]
                  + f_0 * gd_74[k];

        t_45[k] = -dd_21[k]
                  + f_0 * gd_75[k];

        t_46[k] = -dd_22[k]
                  + f_0 * gd_76[k];

        t_47[k] = -dd_23[k]
                  + f_0 * gd_77[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, dd_24, dd_25, dd_26, dd_27, dd_28, \
                         gd_78, gd_79, gd_80, gd_81, gd_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -2.0 * dd_24[k]
                  + f_0 * gd_78[k];

        t_49[k] = -2.0 * dd_25[k]
                  + f_0 * gd_79[k];

        t_50[k] = -2.0 * dd_26[k]
                  + f_0 * gd_80[k];

        t_51[k] = -2.0 * dd_27[k]
                  + f_0 * gd_81[k];

        t_52[k] = -2.0 * dd_28[k]
                  + f_0 * gd_82[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, dd_29, dd_30, dd_31, dd_32, dd_33, \
                         gd_83, gd_84, gd_85, gd_86, gd_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -2.0 * dd_29[k]
                  + f_0 * gd_83[k];

        t_54[k] = -3.0 * dd_30[k]
                  + f_0 * gd_84[k];

        t_55[k] = -3.0 * dd_31[k]
                  + f_0 * gd_85[k];

        t_56[k] = -3.0 * dd_32[k]
                  + f_0 * gd_86[k];

        t_57[k] = -3.0 * dd_33[k]
                  + f_0 * gd_87[k];
    }

#pragma omp simd aligned(t_58, t_59, dd_34, dd_35, gd_88, gd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -3.0 * dd_34[k]
                  + f_0 * gd_88[k];

        t_59[k] = -3.0 * dd_35[k]
                  + f_0 * gd_89[k];
    }
}

}  // namespace simdt2ceri
