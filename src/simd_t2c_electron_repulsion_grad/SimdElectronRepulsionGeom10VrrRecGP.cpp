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


#include "SimdElectronRepulsionGeom10VrrRecGP.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_gp_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t fp, const size_t hp,
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

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);
    const auto *fp_15 = buffer.data(fp + 15);
    const auto *fp_16 = buffer.data(fp + 16);
    const auto *fp_17 = buffer.data(fp + 17);
    const auto *fp_18 = buffer.data(fp + 18);
    const auto *fp_19 = buffer.data(fp + 19);
    const auto *fp_20 = buffer.data(fp + 20);
    const auto *fp_21 = buffer.data(fp + 21);
    const auto *fp_22 = buffer.data(fp + 22);
    const auto *fp_23 = buffer.data(fp + 23);
    const auto *fp_24 = buffer.data(fp + 24);
    const auto *fp_25 = buffer.data(fp + 25);
    const auto *fp_26 = buffer.data(fp + 26);
    const auto *fp_27 = buffer.data(fp + 27);
    const auto *fp_28 = buffer.data(fp + 28);
    const auto *fp_29 = buffer.data(fp + 29);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);
    const auto *hp_30 = buffer.data(hp + 30);
    const auto *hp_31 = buffer.data(hp + 31);
    const auto *hp_32 = buffer.data(hp + 32);
    const auto *hp_33 = buffer.data(hp + 33);
    const auto *hp_34 = buffer.data(hp + 34);
    const auto *hp_35 = buffer.data(hp + 35);
    const auto *hp_36 = buffer.data(hp + 36);
    const auto *hp_37 = buffer.data(hp + 37);
    const auto *hp_38 = buffer.data(hp + 38);
    const auto *hp_39 = buffer.data(hp + 39);
    const auto *hp_40 = buffer.data(hp + 40);
    const auto *hp_41 = buffer.data(hp + 41);
    const auto *hp_42 = buffer.data(hp + 42);
    const auto *hp_43 = buffer.data(hp + 43);
    const auto *hp_44 = buffer.data(hp + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, fp_0, fp_1, fp_2, fp_3, fp_4, hp_0, hp_1, \
                         hp_2, hp_3, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -4.0 * fp_0[k]
                 + f_0 * hp_0[k];

        t_1[k] = -4.0 * fp_1[k]
                 + f_0 * hp_1[k];

        t_2[k] = -4.0 * fp_2[k]
                 + f_0 * hp_2[k];

        t_3[k] = -3.0 * fp_3[k]
                 + f_0 * hp_3[k];

        t_4[k] = -3.0 * fp_4[k]
                 + f_0 * hp_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, fp_5, fp_6, fp_7, fp_8, fp_9, hp_5, hp_6, \
                         hp_7, hp_8, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -3.0 * fp_5[k]
                 + f_0 * hp_5[k];

        t_6[k] = -3.0 * fp_6[k]
                 + f_0 * hp_6[k];

        t_7[k] = -3.0 * fp_7[k]
                 + f_0 * hp_7[k];

        t_8[k] = -3.0 * fp_8[k]
                 + f_0 * hp_8[k];

        t_9[k] = -2.0 * fp_9[k]
                 + f_0 * hp_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, fp_10, fp_11, fp_12, fp_13, fp_14, \
                         hp_10, hp_11, hp_12, hp_13, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -2.0 * fp_10[k]
                  + f_0 * hp_10[k];

        t_11[k] = -2.0 * fp_11[k]
                  + f_0 * hp_11[k];

        t_12[k] = -2.0 * fp_12[k]
                  + f_0 * hp_12[k];

        t_13[k] = -2.0 * fp_13[k]
                  + f_0 * hp_13[k];

        t_14[k] = -2.0 * fp_14[k]
                  + f_0 * hp_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, fp_15, fp_16, fp_17, fp_18, fp_19, \
                         hp_15, hp_16, hp_17, hp_18, hp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -2.0 * fp_15[k]
                  + f_0 * hp_15[k];

        t_16[k] = -2.0 * fp_16[k]
                  + f_0 * hp_16[k];

        t_17[k] = -2.0 * fp_17[k]
                  + f_0 * hp_17[k];

        t_18[k] = -fp_18[k]
                  + f_0 * hp_18[k];

        t_19[k] = -fp_19[k]
                  + f_0 * hp_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, fp_20, fp_21, fp_22, fp_23, fp_24, \
                         hp_20, hp_21, hp_22, hp_23, hp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -fp_20[k]
                  + f_0 * hp_20[k];

        t_21[k] = -fp_21[k]
                  + f_0 * hp_21[k];

        t_22[k] = -fp_22[k]
                  + f_0 * hp_22[k];

        t_23[k] = -fp_23[k]
                  + f_0 * hp_23[k];

        t_24[k] = -fp_24[k]
                  + f_0 * hp_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, fp_25, fp_26, fp_27, fp_28, fp_29, \
                         hp_25, hp_26, hp_27, hp_28, hp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -fp_25[k]
                  + f_0 * hp_25[k];

        t_26[k] = -fp_26[k]
                  + f_0 * hp_26[k];

        t_27[k] = -fp_27[k]
                  + f_0 * hp_27[k];

        t_28[k] = -fp_28[k]
                  + f_0 * hp_28[k];

        t_29[k] = -fp_29[k]
                  + f_0 * hp_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, t_36, t_37, hp_30, hp_31, hp_32, \
                         hp_33, hp_34, hp_35, hp_36, hp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * hp_30[k];

        t_31[k] = f_0 * hp_31[k];

        t_32[k] = f_0 * hp_32[k];

        t_33[k] = f_0 * hp_33[k];

        t_34[k] = f_0 * hp_34[k];

        t_35[k] = f_0 * hp_35[k];

        t_36[k] = f_0 * hp_36[k];

        t_37[k] = f_0 * hp_37[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, t_44, hp_38, hp_39, hp_40, hp_41, \
                         hp_42, hp_43, hp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * hp_38[k];

        t_39[k] = f_0 * hp_39[k];

        t_40[k] = f_0 * hp_40[k];

        t_41[k] = f_0 * hp_41[k];

        t_42[k] = f_0 * hp_42[k];

        t_43[k] = f_0 * hp_43[k];

        t_44[k] = f_0 * hp_44[k];
    }
}

auto
compute_prim_geom_10_gp_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t fp, const size_t hp,
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

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);
    const auto *fp_15 = buffer.data(fp + 15);
    const auto *fp_16 = buffer.data(fp + 16);
    const auto *fp_17 = buffer.data(fp + 17);
    const auto *fp_18 = buffer.data(fp + 18);
    const auto *fp_19 = buffer.data(fp + 19);
    const auto *fp_20 = buffer.data(fp + 20);
    const auto *fp_21 = buffer.data(fp + 21);
    const auto *fp_22 = buffer.data(fp + 22);
    const auto *fp_23 = buffer.data(fp + 23);
    const auto *fp_24 = buffer.data(fp + 24);
    const auto *fp_25 = buffer.data(fp + 25);
    const auto *fp_26 = buffer.data(fp + 26);
    const auto *fp_27 = buffer.data(fp + 27);
    const auto *fp_28 = buffer.data(fp + 28);
    const auto *fp_29 = buffer.data(fp + 29);

    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_30 = buffer.data(hp + 30);
    const auto *hp_31 = buffer.data(hp + 31);
    const auto *hp_32 = buffer.data(hp + 32);
    const auto *hp_33 = buffer.data(hp + 33);
    const auto *hp_34 = buffer.data(hp + 34);
    const auto *hp_35 = buffer.data(hp + 35);
    const auto *hp_36 = buffer.data(hp + 36);
    const auto *hp_37 = buffer.data(hp + 37);
    const auto *hp_38 = buffer.data(hp + 38);
    const auto *hp_39 = buffer.data(hp + 39);
    const auto *hp_40 = buffer.data(hp + 40);
    const auto *hp_41 = buffer.data(hp + 41);
    const auto *hp_45 = buffer.data(hp + 45);
    const auto *hp_46 = buffer.data(hp + 46);
    const auto *hp_47 = buffer.data(hp + 47);
    const auto *hp_48 = buffer.data(hp + 48);
    const auto *hp_49 = buffer.data(hp + 49);
    const auto *hp_50 = buffer.data(hp + 50);
    const auto *hp_51 = buffer.data(hp + 51);
    const auto *hp_52 = buffer.data(hp + 52);
    const auto *hp_53 = buffer.data(hp + 53);
    const auto *hp_54 = buffer.data(hp + 54);
    const auto *hp_55 = buffer.data(hp + 55);
    const auto *hp_56 = buffer.data(hp + 56);
    const auto *hp_57 = buffer.data(hp + 57);
    const auto *hp_58 = buffer.data(hp + 58);
    const auto *hp_59 = buffer.data(hp + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, fp_0, fp_1, fp_2, hp_3, hp_4, hp_5, \
                         hp_9, hp_10, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_3[k];

        t_1[k] = f_0 * hp_4[k];

        t_2[k] = f_0 * hp_5[k];

        t_3[k] = -fp_0[k]
                 + f_0 * hp_9[k];

        t_4[k] = -fp_1[k]
                 + f_0 * hp_10[k];

        t_5[k] = -fp_2[k]
                 + f_0 * hp_11[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, fp_3, fp_4, fp_5, hp_12, hp_13, \
                         hp_14, hp_18, hp_19, hp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * hp_12[k];

        t_7[k] = f_0 * hp_13[k];

        t_8[k] = f_0 * hp_14[k];

        t_9[k] = -2.0 * fp_3[k]
                 + f_0 * hp_18[k];

        t_10[k] = -2.0 * fp_4[k]
                  + f_0 * hp_19[k];

        t_11[k] = -2.0 * fp_5[k]
                  + f_0 * hp_20[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, fp_6, fp_7, fp_8, hp_21, hp_22, \
                         hp_23, hp_24, hp_25, hp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -fp_6[k]
                  + f_0 * hp_21[k];

        t_13[k] = -fp_7[k]
                  + f_0 * hp_22[k];

        t_14[k] = -fp_8[k]
                  + f_0 * hp_23[k];

        t_15[k] = f_0 * hp_24[k];

        t_16[k] = f_0 * hp_25[k];

        t_17[k] = f_0 * hp_26[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, fp_9, fp_10, fp_11, fp_12, fp_13, \
                         hp_30, hp_31, hp_32, hp_33, hp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -3.0 * fp_9[k]
                  + f_0 * hp_30[k];

        t_19[k] = -3.0 * fp_10[k]
                  + f_0 * hp_31[k];

        t_20[k] = -3.0 * fp_11[k]
                  + f_0 * hp_32[k];

        t_21[k] = -2.0 * fp_12[k]
                  + f_0 * hp_33[k];

        t_22[k] = -2.0 * fp_13[k]
                  + f_0 * hp_34[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, t_28, fp_14, fp_15, fp_16, fp_17, \
                         hp_35, hp_36, hp_37, hp_38, hp_39, hp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -2.0 * fp_14[k]
                  + f_0 * hp_35[k];

        t_24[k] = -fp_15[k]
                  + f_0 * hp_36[k];

        t_25[k] = -fp_16[k]
                  + f_0 * hp_37[k];

        t_26[k] = -fp_17[k]
                  + f_0 * hp_38[k];

        t_27[k] = f_0 * hp_39[k];

        t_28[k] = f_0 * hp_40[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, fp_18, fp_19, fp_20, fp_21, hp_41, \
                         hp_45, hp_46, hp_47, hp_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * hp_41[k];

        t_30[k] = -4.0 * fp_18[k]
                  + f_0 * hp_45[k];

        t_31[k] = -4.0 * fp_19[k]
                  + f_0 * hp_46[k];

        t_32[k] = -4.0 * fp_20[k]
                  + f_0 * hp_47[k];

        t_33[k] = -3.0 * fp_21[k]
                  + f_0 * hp_48[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, fp_22, fp_23, fp_24, fp_25, fp_26, \
                         hp_49, hp_50, hp_51, hp_52, hp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = -3.0 * fp_22[k]
                  + f_0 * hp_49[k];

        t_35[k] = -3.0 * fp_23[k]
                  + f_0 * hp_50[k];

        t_36[k] = -2.0 * fp_24[k]
                  + f_0 * hp_51[k];

        t_37[k] = -2.0 * fp_25[k]
                  + f_0 * hp_52[k];

        t_38[k] = -2.0 * fp_26[k]
                  + f_0 * hp_53[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, t_44, fp_27, fp_28, fp_29, hp_54, \
                         hp_55, hp_56, hp_57, hp_58, hp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -fp_27[k]
                  + f_0 * hp_54[k];

        t_40[k] = -fp_28[k]
                  + f_0 * hp_55[k];

        t_41[k] = -fp_29[k]
                  + f_0 * hp_56[k];

        t_42[k] = f_0 * hp_57[k];

        t_43[k] = f_0 * hp_58[k];

        t_44[k] = f_0 * hp_59[k];
    }
}

auto
compute_prim_geom_10_gp_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t fp, const size_t hp,
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

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);
    const auto *fp_15 = buffer.data(fp + 15);
    const auto *fp_16 = buffer.data(fp + 16);
    const auto *fp_17 = buffer.data(fp + 17);
    const auto *fp_18 = buffer.data(fp + 18);
    const auto *fp_19 = buffer.data(fp + 19);
    const auto *fp_20 = buffer.data(fp + 20);
    const auto *fp_21 = buffer.data(fp + 21);
    const auto *fp_22 = buffer.data(fp + 22);
    const auto *fp_23 = buffer.data(fp + 23);
    const auto *fp_24 = buffer.data(fp + 24);
    const auto *fp_25 = buffer.data(fp + 25);
    const auto *fp_26 = buffer.data(fp + 26);
    const auto *fp_27 = buffer.data(fp + 27);
    const auto *fp_28 = buffer.data(fp + 28);
    const auto *fp_29 = buffer.data(fp + 29);

    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);
    const auto *hp_33 = buffer.data(hp + 33);
    const auto *hp_34 = buffer.data(hp + 34);
    const auto *hp_35 = buffer.data(hp + 35);
    const auto *hp_36 = buffer.data(hp + 36);
    const auto *hp_37 = buffer.data(hp + 37);
    const auto *hp_38 = buffer.data(hp + 38);
    const auto *hp_39 = buffer.data(hp + 39);
    const auto *hp_40 = buffer.data(hp + 40);
    const auto *hp_41 = buffer.data(hp + 41);
    const auto *hp_42 = buffer.data(hp + 42);
    const auto *hp_43 = buffer.data(hp + 43);
    const auto *hp_44 = buffer.data(hp + 44);
    const auto *hp_48 = buffer.data(hp + 48);
    const auto *hp_49 = buffer.data(hp + 49);
    const auto *hp_50 = buffer.data(hp + 50);
    const auto *hp_51 = buffer.data(hp + 51);
    const auto *hp_52 = buffer.data(hp + 52);
    const auto *hp_53 = buffer.data(hp + 53);
    const auto *hp_54 = buffer.data(hp + 54);
    const auto *hp_55 = buffer.data(hp + 55);
    const auto *hp_56 = buffer.data(hp + 56);
    const auto *hp_57 = buffer.data(hp + 57);
    const auto *hp_58 = buffer.data(hp + 58);
    const auto *hp_59 = buffer.data(hp + 59);
    const auto *hp_60 = buffer.data(hp + 60);
    const auto *hp_61 = buffer.data(hp + 61);
    const auto *hp_62 = buffer.data(hp + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, fp_0, hp_6, hp_7, hp_8, hp_12, \
                         hp_13, hp_14, hp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_6[k];

        t_1[k] = f_0 * hp_7[k];

        t_2[k] = f_0 * hp_8[k];

        t_3[k] = f_0 * hp_12[k];

        t_4[k] = f_0 * hp_13[k];

        t_5[k] = f_0 * hp_14[k];

        t_6[k] = -fp_0[k]
                 + f_0 * hp_15[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, fp_1, fp_2, fp_3, hp_16, hp_17, \
                         hp_21, hp_22, hp_23, hp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -fp_1[k]
                 + f_0 * hp_16[k];

        t_8[k] = -fp_2[k]
                 + f_0 * hp_17[k];

        t_9[k] = f_0 * hp_21[k];

        t_10[k] = f_0 * hp_22[k];

        t_11[k] = f_0 * hp_23[k];

        t_12[k] = -fp_3[k]
                  + f_0 * hp_24[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, fp_4, fp_5, fp_6, fp_7, fp_8, hp_25, \
                         hp_26, hp_27, hp_28, hp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -fp_4[k]
                  + f_0 * hp_25[k];

        t_14[k] = -fp_5[k]
                  + f_0 * hp_26[k];

        t_15[k] = -2.0 * fp_6[k]
                  + f_0 * hp_27[k];

        t_16[k] = -2.0 * fp_7[k]
                  + f_0 * hp_28[k];

        t_17[k] = -2.0 * fp_8[k]
                  + f_0 * hp_29[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, fp_9, fp_10, fp_11, hp_33, hp_34, \
                         hp_35, hp_36, hp_37, hp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * hp_33[k];

        t_19[k] = f_0 * hp_34[k];

        t_20[k] = f_0 * hp_35[k];

        t_21[k] = -fp_9[k]
                  + f_0 * hp_36[k];

        t_22[k] = -fp_10[k]
                  + f_0 * hp_37[k];

        t_23[k] = -fp_11[k]
                  + f_0 * hp_38[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, fp_12, fp_13, fp_14, fp_15, fp_16, \
                         hp_39, hp_40, hp_41, hp_42, hp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -2.0 * fp_12[k]
                  + f_0 * hp_39[k];

        t_25[k] = -2.0 * fp_13[k]
                  + f_0 * hp_40[k];

        t_26[k] = -2.0 * fp_14[k]
                  + f_0 * hp_41[k];

        t_27[k] = -3.0 * fp_15[k]
                  + f_0 * hp_42[k];

        t_28[k] = -3.0 * fp_16[k]
                  + f_0 * hp_43[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, fp_17, fp_18, fp_19, hp_44, \
                         hp_48, hp_49, hp_50, hp_51, hp_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -3.0 * fp_17[k]
                  + f_0 * hp_44[k];

        t_30[k] = f_0 * hp_48[k];

        t_31[k] = f_0 * hp_49[k];

        t_32[k] = f_0 * hp_50[k];

        t_33[k] = -fp_18[k]
                  + f_0 * hp_51[k];

        t_34[k] = -fp_19[k]
                  + f_0 * hp_52[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, fp_20, fp_21, fp_22, fp_23, fp_24, \
                         hp_53, hp_54, hp_55, hp_56, hp_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -fp_20[k]
                  + f_0 * hp_53[k];

        t_36[k] = -2.0 * fp_21[k]
                  + f_0 * hp_54[k];

        t_37[k] = -2.0 * fp_22[k]
                  + f_0 * hp_55[k];

        t_38[k] = -2.0 * fp_23[k]
                  + f_0 * hp_56[k];

        t_39[k] = -3.0 * fp_24[k]
                  + f_0 * hp_57[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, fp_25, fp_26, fp_27, fp_28, fp_29, \
                         hp_58, hp_59, hp_60, hp_61, hp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -3.0 * fp_25[k]
                  + f_0 * hp_58[k];

        t_41[k] = -3.0 * fp_26[k]
                  + f_0 * hp_59[k];

        t_42[k] = -4.0 * fp_27[k]
                  + f_0 * hp_60[k];

        t_43[k] = -4.0 * fp_28[k]
                  + f_0 * hp_61[k];

        t_44[k] = -4.0 * fp_29[k]
                  + f_0 * hp_62[k];
    }
}

}  // namespace simdt2ceri
