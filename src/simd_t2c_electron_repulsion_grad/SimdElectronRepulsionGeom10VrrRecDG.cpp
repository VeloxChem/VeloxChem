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


#include "SimdElectronRepulsionGeom10VrrRecDG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_dg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t pg, const size_t fg,
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
    auto *t_63 = buffer.data(target + 63);
    auto *t_64 = buffer.data(target + 64);
    auto *t_65 = buffer.data(target + 65);
    auto *t_66 = buffer.data(target + 66);
    auto *t_67 = buffer.data(target + 67);
    auto *t_68 = buffer.data(target + 68);
    auto *t_69 = buffer.data(target + 69);
    auto *t_70 = buffer.data(target + 70);
    auto *t_71 = buffer.data(target + 71);
    auto *t_72 = buffer.data(target + 72);
    auto *t_73 = buffer.data(target + 73);
    auto *t_74 = buffer.data(target + 74);
    auto *t_75 = buffer.data(target + 75);
    auto *t_76 = buffer.data(target + 76);
    auto *t_77 = buffer.data(target + 77);
    auto *t_78 = buffer.data(target + 78);
    auto *t_79 = buffer.data(target + 79);
    auto *t_80 = buffer.data(target + 80);
    auto *t_81 = buffer.data(target + 81);
    auto *t_82 = buffer.data(target + 82);
    auto *t_83 = buffer.data(target + 83);
    auto *t_84 = buffer.data(target + 84);
    auto *t_85 = buffer.data(target + 85);
    auto *t_86 = buffer.data(target + 86);
    auto *t_87 = buffer.data(target + 87);
    auto *t_88 = buffer.data(target + 88);
    auto *t_89 = buffer.data(target + 89);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_11 = buffer.data(pg + 11);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_14 = buffer.data(pg + 14);
    const auto *pg_15 = buffer.data(pg + 15);
    const auto *pg_16 = buffer.data(pg + 16);
    const auto *pg_17 = buffer.data(pg + 17);
    const auto *pg_18 = buffer.data(pg + 18);
    const auto *pg_19 = buffer.data(pg + 19);
    const auto *pg_20 = buffer.data(pg + 20);
    const auto *pg_21 = buffer.data(pg + 21);
    const auto *pg_22 = buffer.data(pg + 22);
    const auto *pg_23 = buffer.data(pg + 23);
    const auto *pg_24 = buffer.data(pg + 24);
    const auto *pg_25 = buffer.data(pg + 25);
    const auto *pg_26 = buffer.data(pg + 26);
    const auto *pg_27 = buffer.data(pg + 27);
    const auto *pg_28 = buffer.data(pg + 28);
    const auto *pg_29 = buffer.data(pg + 29);
    const auto *pg_30 = buffer.data(pg + 30);
    const auto *pg_31 = buffer.data(pg + 31);
    const auto *pg_32 = buffer.data(pg + 32);
    const auto *pg_33 = buffer.data(pg + 33);
    const auto *pg_34 = buffer.data(pg + 34);
    const auto *pg_35 = buffer.data(pg + 35);
    const auto *pg_36 = buffer.data(pg + 36);
    const auto *pg_37 = buffer.data(pg + 37);
    const auto *pg_38 = buffer.data(pg + 38);
    const auto *pg_39 = buffer.data(pg + 39);
    const auto *pg_40 = buffer.data(pg + 40);
    const auto *pg_41 = buffer.data(pg + 41);
    const auto *pg_42 = buffer.data(pg + 42);
    const auto *pg_43 = buffer.data(pg + 43);
    const auto *pg_44 = buffer.data(pg + 44);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_4 = buffer.data(fg + 4);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_7 = buffer.data(fg + 7);
    const auto *fg_8 = buffer.data(fg + 8);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_11 = buffer.data(fg + 11);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);
    const auto *fg_59 = buffer.data(fg + 59);
    const auto *fg_60 = buffer.data(fg + 60);
    const auto *fg_61 = buffer.data(fg + 61);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_63 = buffer.data(fg + 63);
    const auto *fg_64 = buffer.data(fg + 64);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_66 = buffer.data(fg + 66);
    const auto *fg_67 = buffer.data(fg + 67);
    const auto *fg_68 = buffer.data(fg + 68);
    const auto *fg_69 = buffer.data(fg + 69);
    const auto *fg_70 = buffer.data(fg + 70);
    const auto *fg_71 = buffer.data(fg + 71);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_73 = buffer.data(fg + 73);
    const auto *fg_74 = buffer.data(fg + 74);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_76 = buffer.data(fg + 76);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_78 = buffer.data(fg + 78);
    const auto *fg_79 = buffer.data(fg + 79);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_81 = buffer.data(fg + 81);
    const auto *fg_82 = buffer.data(fg + 82);
    const auto *fg_83 = buffer.data(fg + 83);
    const auto *fg_84 = buffer.data(fg + 84);
    const auto *fg_85 = buffer.data(fg + 85);
    const auto *fg_86 = buffer.data(fg + 86);
    const auto *fg_87 = buffer.data(fg + 87);
    const auto *fg_88 = buffer.data(fg + 88);
    const auto *fg_89 = buffer.data(fg + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pg_0, pg_1, pg_2, pg_3, pg_4, fg_0, fg_1, \
                         fg_2, fg_3, fg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -2.0 * pg_0[k]
                 + f_0 * fg_0[k];

        t_1[k] = -2.0 * pg_1[k]
                 + f_0 * fg_1[k];

        t_2[k] = -2.0 * pg_2[k]
                 + f_0 * fg_2[k];

        t_3[k] = -2.0 * pg_3[k]
                 + f_0 * fg_3[k];

        t_4[k] = -2.0 * pg_4[k]
                 + f_0 * fg_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pg_5, pg_6, pg_7, pg_8, pg_9, fg_5, fg_6, \
                         fg_7, fg_8, fg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -2.0 * pg_5[k]
                 + f_0 * fg_5[k];

        t_6[k] = -2.0 * pg_6[k]
                 + f_0 * fg_6[k];

        t_7[k] = -2.0 * pg_7[k]
                 + f_0 * fg_7[k];

        t_8[k] = -2.0 * pg_8[k]
                 + f_0 * fg_8[k];

        t_9[k] = -2.0 * pg_9[k]
                 + f_0 * fg_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pg_10, pg_11, pg_12, pg_13, pg_14, \
                         fg_10, fg_11, fg_12, fg_13, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -2.0 * pg_10[k]
                  + f_0 * fg_10[k];

        t_11[k] = -2.0 * pg_11[k]
                  + f_0 * fg_11[k];

        t_12[k] = -2.0 * pg_12[k]
                  + f_0 * fg_12[k];

        t_13[k] = -2.0 * pg_13[k]
                  + f_0 * fg_13[k];

        t_14[k] = -2.0 * pg_14[k]
                  + f_0 * fg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pg_15, pg_16, pg_17, pg_18, pg_19, \
                         fg_15, fg_16, fg_17, fg_18, fg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -pg_15[k]
                  + f_0 * fg_15[k];

        t_16[k] = -pg_16[k]
                  + f_0 * fg_16[k];

        t_17[k] = -pg_17[k]
                  + f_0 * fg_17[k];

        t_18[k] = -pg_18[k]
                  + f_0 * fg_18[k];

        t_19[k] = -pg_19[k]
                  + f_0 * fg_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pg_20, pg_21, pg_22, pg_23, pg_24, \
                         fg_20, fg_21, fg_22, fg_23, fg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -pg_20[k]
                  + f_0 * fg_20[k];

        t_21[k] = -pg_21[k]
                  + f_0 * fg_21[k];

        t_22[k] = -pg_22[k]
                  + f_0 * fg_22[k];

        t_23[k] = -pg_23[k]
                  + f_0 * fg_23[k];

        t_24[k] = -pg_24[k]
                  + f_0 * fg_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pg_25, pg_26, pg_27, pg_28, pg_29, \
                         fg_25, fg_26, fg_27, fg_28, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -pg_25[k]
                  + f_0 * fg_25[k];

        t_26[k] = -pg_26[k]
                  + f_0 * fg_26[k];

        t_27[k] = -pg_27[k]
                  + f_0 * fg_27[k];

        t_28[k] = -pg_28[k]
                  + f_0 * fg_28[k];

        t_29[k] = -pg_29[k]
                  + f_0 * fg_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pg_30, pg_31, pg_32, pg_33, pg_34, \
                         fg_30, fg_31, fg_32, fg_33, fg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -pg_30[k]
                  + f_0 * fg_30[k];

        t_31[k] = -pg_31[k]
                  + f_0 * fg_31[k];

        t_32[k] = -pg_32[k]
                  + f_0 * fg_32[k];

        t_33[k] = -pg_33[k]
                  + f_0 * fg_33[k];

        t_34[k] = -pg_34[k]
                  + f_0 * fg_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pg_35, pg_36, pg_37, pg_38, pg_39, \
                         fg_35, fg_36, fg_37, fg_38, fg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -pg_35[k]
                  + f_0 * fg_35[k];

        t_36[k] = -pg_36[k]
                  + f_0 * fg_36[k];

        t_37[k] = -pg_37[k]
                  + f_0 * fg_37[k];

        t_38[k] = -pg_38[k]
                  + f_0 * fg_38[k];

        t_39[k] = -pg_39[k]
                  + f_0 * fg_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pg_40, pg_41, pg_42, pg_43, pg_44, \
                         fg_40, fg_41, fg_42, fg_43, fg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -pg_40[k]
                  + f_0 * fg_40[k];

        t_41[k] = -pg_41[k]
                  + f_0 * fg_41[k];

        t_42[k] = -pg_42[k]
                  + f_0 * fg_42[k];

        t_43[k] = -pg_43[k]
                  + f_0 * fg_43[k];

        t_44[k] = -pg_44[k]
                  + f_0 * fg_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, t_51, t_52, fg_45, fg_46, fg_47, \
                         fg_48, fg_49, fg_50, fg_51, fg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_0 * fg_45[k];

        t_46[k] = f_0 * fg_46[k];

        t_47[k] = f_0 * fg_47[k];

        t_48[k] = f_0 * fg_48[k];

        t_49[k] = f_0 * fg_49[k];

        t_50[k] = f_0 * fg_50[k];

        t_51[k] = f_0 * fg_51[k];

        t_52[k] = f_0 * fg_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, t_58, t_59, t_60, fg_53, fg_54, fg_55, \
                         fg_56, fg_57, fg_58, fg_59, fg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * fg_53[k];

        t_54[k] = f_0 * fg_54[k];

        t_55[k] = f_0 * fg_55[k];

        t_56[k] = f_0 * fg_56[k];

        t_57[k] = f_0 * fg_57[k];

        t_58[k] = f_0 * fg_58[k];

        t_59[k] = f_0 * fg_59[k];

        t_60[k] = f_0 * fg_60[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, t_66, t_67, t_68, fg_61, fg_62, fg_63, \
                         fg_64, fg_65, fg_66, fg_67, fg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_0 * fg_61[k];

        t_62[k] = f_0 * fg_62[k];

        t_63[k] = f_0 * fg_63[k];

        t_64[k] = f_0 * fg_64[k];

        t_65[k] = f_0 * fg_65[k];

        t_66[k] = f_0 * fg_66[k];

        t_67[k] = f_0 * fg_67[k];

        t_68[k] = f_0 * fg_68[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, t_74, t_75, t_76, fg_69, fg_70, fg_71, \
                         fg_72, fg_73, fg_74, fg_75, fg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_0 * fg_69[k];

        t_70[k] = f_0 * fg_70[k];

        t_71[k] = f_0 * fg_71[k];

        t_72[k] = f_0 * fg_72[k];

        t_73[k] = f_0 * fg_73[k];

        t_74[k] = f_0 * fg_74[k];

        t_75[k] = f_0 * fg_75[k];

        t_76[k] = f_0 * fg_76[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, t_82, t_83, t_84, fg_77, fg_78, fg_79, \
                         fg_80, fg_81, fg_82, fg_83, fg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_0 * fg_77[k];

        t_78[k] = f_0 * fg_78[k];

        t_79[k] = f_0 * fg_79[k];

        t_80[k] = f_0 * fg_80[k];

        t_81[k] = f_0 * fg_81[k];

        t_82[k] = f_0 * fg_82[k];

        t_83[k] = f_0 * fg_83[k];

        t_84[k] = f_0 * fg_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, fg_85, fg_86, fg_87, fg_88, \
                         fg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_0 * fg_85[k];

        t_86[k] = f_0 * fg_86[k];

        t_87[k] = f_0 * fg_87[k];

        t_88[k] = f_0 * fg_88[k];

        t_89[k] = f_0 * fg_89[k];
    }
}

auto
compute_prim_geom_10_dg_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t pg, const size_t fg,
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
    auto *t_63 = buffer.data(target + 63);
    auto *t_64 = buffer.data(target + 64);
    auto *t_65 = buffer.data(target + 65);
    auto *t_66 = buffer.data(target + 66);
    auto *t_67 = buffer.data(target + 67);
    auto *t_68 = buffer.data(target + 68);
    auto *t_69 = buffer.data(target + 69);
    auto *t_70 = buffer.data(target + 70);
    auto *t_71 = buffer.data(target + 71);
    auto *t_72 = buffer.data(target + 72);
    auto *t_73 = buffer.data(target + 73);
    auto *t_74 = buffer.data(target + 74);
    auto *t_75 = buffer.data(target + 75);
    auto *t_76 = buffer.data(target + 76);
    auto *t_77 = buffer.data(target + 77);
    auto *t_78 = buffer.data(target + 78);
    auto *t_79 = buffer.data(target + 79);
    auto *t_80 = buffer.data(target + 80);
    auto *t_81 = buffer.data(target + 81);
    auto *t_82 = buffer.data(target + 82);
    auto *t_83 = buffer.data(target + 83);
    auto *t_84 = buffer.data(target + 84);
    auto *t_85 = buffer.data(target + 85);
    auto *t_86 = buffer.data(target + 86);
    auto *t_87 = buffer.data(target + 87);
    auto *t_88 = buffer.data(target + 88);
    auto *t_89 = buffer.data(target + 89);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_11 = buffer.data(pg + 11);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_14 = buffer.data(pg + 14);
    const auto *pg_15 = buffer.data(pg + 15);
    const auto *pg_16 = buffer.data(pg + 16);
    const auto *pg_17 = buffer.data(pg + 17);
    const auto *pg_18 = buffer.data(pg + 18);
    const auto *pg_19 = buffer.data(pg + 19);
    const auto *pg_20 = buffer.data(pg + 20);
    const auto *pg_21 = buffer.data(pg + 21);
    const auto *pg_22 = buffer.data(pg + 22);
    const auto *pg_23 = buffer.data(pg + 23);
    const auto *pg_24 = buffer.data(pg + 24);
    const auto *pg_25 = buffer.data(pg + 25);
    const auto *pg_26 = buffer.data(pg + 26);
    const auto *pg_27 = buffer.data(pg + 27);
    const auto *pg_28 = buffer.data(pg + 28);
    const auto *pg_29 = buffer.data(pg + 29);
    const auto *pg_30 = buffer.data(pg + 30);
    const auto *pg_31 = buffer.data(pg + 31);
    const auto *pg_32 = buffer.data(pg + 32);
    const auto *pg_33 = buffer.data(pg + 33);
    const auto *pg_34 = buffer.data(pg + 34);
    const auto *pg_35 = buffer.data(pg + 35);
    const auto *pg_36 = buffer.data(pg + 36);
    const auto *pg_37 = buffer.data(pg + 37);
    const auto *pg_38 = buffer.data(pg + 38);
    const auto *pg_39 = buffer.data(pg + 39);
    const auto *pg_40 = buffer.data(pg + 40);
    const auto *pg_41 = buffer.data(pg + 41);
    const auto *pg_42 = buffer.data(pg + 42);
    const auto *pg_43 = buffer.data(pg + 43);
    const auto *pg_44 = buffer.data(pg + 44);

    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_17 = buffer.data(fg + 17);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_19 = buffer.data(fg + 19);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_22 = buffer.data(fg + 22);
    const auto *fg_23 = buffer.data(fg + 23);
    const auto *fg_24 = buffer.data(fg + 24);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_49 = buffer.data(fg + 49);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_52 = buffer.data(fg + 52);
    const auto *fg_53 = buffer.data(fg + 53);
    const auto *fg_54 = buffer.data(fg + 54);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_56 = buffer.data(fg + 56);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);
    const auto *fg_59 = buffer.data(fg + 59);
    const auto *fg_60 = buffer.data(fg + 60);
    const auto *fg_61 = buffer.data(fg + 61);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_63 = buffer.data(fg + 63);
    const auto *fg_64 = buffer.data(fg + 64);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_66 = buffer.data(fg + 66);
    const auto *fg_67 = buffer.data(fg + 67);
    const auto *fg_68 = buffer.data(fg + 68);
    const auto *fg_69 = buffer.data(fg + 69);
    const auto *fg_70 = buffer.data(fg + 70);
    const auto *fg_71 = buffer.data(fg + 71);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_73 = buffer.data(fg + 73);
    const auto *fg_74 = buffer.data(fg + 74);
    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_91 = buffer.data(fg + 91);
    const auto *fg_92 = buffer.data(fg + 92);
    const auto *fg_93 = buffer.data(fg + 93);
    const auto *fg_94 = buffer.data(fg + 94);
    const auto *fg_95 = buffer.data(fg + 95);
    const auto *fg_96 = buffer.data(fg + 96);
    const auto *fg_97 = buffer.data(fg + 97);
    const auto *fg_98 = buffer.data(fg + 98);
    const auto *fg_99 = buffer.data(fg + 99);
    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_101 = buffer.data(fg + 101);
    const auto *fg_102 = buffer.data(fg + 102);
    const auto *fg_103 = buffer.data(fg + 103);
    const auto *fg_104 = buffer.data(fg + 104);
    const auto *fg_105 = buffer.data(fg + 105);
    const auto *fg_106 = buffer.data(fg + 106);
    const auto *fg_107 = buffer.data(fg + 107);
    const auto *fg_108 = buffer.data(fg + 108);
    const auto *fg_109 = buffer.data(fg + 109);
    const auto *fg_110 = buffer.data(fg + 110);
    const auto *fg_111 = buffer.data(fg + 111);
    const auto *fg_112 = buffer.data(fg + 112);
    const auto *fg_113 = buffer.data(fg + 113);
    const auto *fg_114 = buffer.data(fg + 114);
    const auto *fg_115 = buffer.data(fg + 115);
    const auto *fg_116 = buffer.data(fg + 116);
    const auto *fg_117 = buffer.data(fg + 117);
    const auto *fg_118 = buffer.data(fg + 118);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_120 = buffer.data(fg + 120);
    const auto *fg_121 = buffer.data(fg + 121);
    const auto *fg_122 = buffer.data(fg + 122);
    const auto *fg_123 = buffer.data(fg + 123);
    const auto *fg_124 = buffer.data(fg + 124);
    const auto *fg_125 = buffer.data(fg + 125);
    const auto *fg_126 = buffer.data(fg + 126);
    const auto *fg_127 = buffer.data(fg + 127);
    const auto *fg_128 = buffer.data(fg + 128);
    const auto *fg_129 = buffer.data(fg + 129);
    const auto *fg_130 = buffer.data(fg + 130);
    const auto *fg_131 = buffer.data(fg + 131);
    const auto *fg_132 = buffer.data(fg + 132);
    const auto *fg_133 = buffer.data(fg + 133);
    const auto *fg_134 = buffer.data(fg + 134);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, fg_15, fg_16, fg_17, fg_18, \
                         fg_19, fg_20, fg_21, fg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_15[k];

        t_1[k] = f_0 * fg_16[k];

        t_2[k] = f_0 * fg_17[k];

        t_3[k] = f_0 * fg_18[k];

        t_4[k] = f_0 * fg_19[k];

        t_5[k] = f_0 * fg_20[k];

        t_6[k] = f_0 * fg_21[k];

        t_7[k] = f_0 * fg_22[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, fg_23, fg_24, fg_25, fg_26, \
                         fg_27, fg_28, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * fg_23[k];

        t_9[k] = f_0 * fg_24[k];

        t_10[k] = f_0 * fg_25[k];

        t_11[k] = f_0 * fg_26[k];

        t_12[k] = f_0 * fg_27[k];

        t_13[k] = f_0 * fg_28[k];

        t_14[k] = f_0 * fg_29[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pg_0, pg_1, pg_2, pg_3, pg_4, fg_45, \
                         fg_46, fg_47, fg_48, fg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -pg_0[k]
                  + f_0 * fg_45[k];

        t_16[k] = -pg_1[k]
                  + f_0 * fg_46[k];

        t_17[k] = -pg_2[k]
                  + f_0 * fg_47[k];

        t_18[k] = -pg_3[k]
                  + f_0 * fg_48[k];

        t_19[k] = -pg_4[k]
                  + f_0 * fg_49[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pg_5, pg_6, pg_7, pg_8, pg_9, fg_50, \
                         fg_51, fg_52, fg_53, fg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -pg_5[k]
                  + f_0 * fg_50[k];

        t_21[k] = -pg_6[k]
                  + f_0 * fg_51[k];

        t_22[k] = -pg_7[k]
                  + f_0 * fg_52[k];

        t_23[k] = -pg_8[k]
                  + f_0 * fg_53[k];

        t_24[k] = -pg_9[k]
                  + f_0 * fg_54[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pg_10, pg_11, pg_12, pg_13, pg_14, \
                         fg_55, fg_56, fg_57, fg_58, fg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -pg_10[k]
                  + f_0 * fg_55[k];

        t_26[k] = -pg_11[k]
                  + f_0 * fg_56[k];

        t_27[k] = -pg_12[k]
                  + f_0 * fg_57[k];

        t_28[k] = -pg_13[k]
                  + f_0 * fg_58[k];

        t_29[k] = -pg_14[k]
                  + f_0 * fg_59[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, t_36, t_37, fg_60, fg_61, fg_62, \
                         fg_63, fg_64, fg_65, fg_66, fg_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * fg_60[k];

        t_31[k] = f_0 * fg_61[k];

        t_32[k] = f_0 * fg_62[k];

        t_33[k] = f_0 * fg_63[k];

        t_34[k] = f_0 * fg_64[k];

        t_35[k] = f_0 * fg_65[k];

        t_36[k] = f_0 * fg_66[k];

        t_37[k] = f_0 * fg_67[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, t_44, fg_68, fg_69, fg_70, fg_71, \
                         fg_72, fg_73, fg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * fg_68[k];

        t_39[k] = f_0 * fg_69[k];

        t_40[k] = f_0 * fg_70[k];

        t_41[k] = f_0 * fg_71[k];

        t_42[k] = f_0 * fg_72[k];

        t_43[k] = f_0 * fg_73[k];

        t_44[k] = f_0 * fg_74[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pg_15, pg_16, pg_17, pg_18, pg_19, \
                         fg_90, fg_91, fg_92, fg_93, fg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -2.0 * pg_15[k]
                  + f_0 * fg_90[k];

        t_46[k] = -2.0 * pg_16[k]
                  + f_0 * fg_91[k];

        t_47[k] = -2.0 * pg_17[k]
                  + f_0 * fg_92[k];

        t_48[k] = -2.0 * pg_18[k]
                  + f_0 * fg_93[k];

        t_49[k] = -2.0 * pg_19[k]
                  + f_0 * fg_94[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pg_20, pg_21, pg_22, pg_23, pg_24, \
                         fg_95, fg_96, fg_97, fg_98, fg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * pg_20[k]
                  + f_0 * fg_95[k];

        t_51[k] = -2.0 * pg_21[k]
                  + f_0 * fg_96[k];

        t_52[k] = -2.0 * pg_22[k]
                  + f_0 * fg_97[k];

        t_53[k] = -2.0 * pg_23[k]
                  + f_0 * fg_98[k];

        t_54[k] = -2.0 * pg_24[k]
                  + f_0 * fg_99[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pg_25, pg_26, pg_27, pg_28, pg_29, \
                         fg_100, fg_101, fg_102, fg_103, fg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * pg_25[k]
                  + f_0 * fg_100[k];

        t_56[k] = -2.0 * pg_26[k]
                  + f_0 * fg_101[k];

        t_57[k] = -2.0 * pg_27[k]
                  + f_0 * fg_102[k];

        t_58[k] = -2.0 * pg_28[k]
                  + f_0 * fg_103[k];

        t_59[k] = -2.0 * pg_29[k]
                  + f_0 * fg_104[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pg_30, pg_31, pg_32, pg_33, pg_34, \
                         fg_105, fg_106, fg_107, fg_108, fg_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -pg_30[k]
                  + f_0 * fg_105[k];

        t_61[k] = -pg_31[k]
                  + f_0 * fg_106[k];

        t_62[k] = -pg_32[k]
                  + f_0 * fg_107[k];

        t_63[k] = -pg_33[k]
                  + f_0 * fg_108[k];

        t_64[k] = -pg_34[k]
                  + f_0 * fg_109[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pg_35, pg_36, pg_37, pg_38, pg_39, \
                         fg_110, fg_111, fg_112, fg_113, fg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -pg_35[k]
                  + f_0 * fg_110[k];

        t_66[k] = -pg_36[k]
                  + f_0 * fg_111[k];

        t_67[k] = -pg_37[k]
                  + f_0 * fg_112[k];

        t_68[k] = -pg_38[k]
                  + f_0 * fg_113[k];

        t_69[k] = -pg_39[k]
                  + f_0 * fg_114[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pg_40, pg_41, pg_42, pg_43, pg_44, \
                         fg_115, fg_116, fg_117, fg_118, fg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -pg_40[k]
                  + f_0 * fg_115[k];

        t_71[k] = -pg_41[k]
                  + f_0 * fg_116[k];

        t_72[k] = -pg_42[k]
                  + f_0 * fg_117[k];

        t_73[k] = -pg_43[k]
                  + f_0 * fg_118[k];

        t_74[k] = -pg_44[k]
                  + f_0 * fg_119[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, t_80, t_81, t_82, fg_120, fg_121, \
                         fg_122, fg_123, fg_124, fg_125, fg_126, \
                         fg_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_0 * fg_120[k];

        t_76[k] = f_0 * fg_121[k];

        t_77[k] = f_0 * fg_122[k];

        t_78[k] = f_0 * fg_123[k];

        t_79[k] = f_0 * fg_124[k];

        t_80[k] = f_0 * fg_125[k];

        t_81[k] = f_0 * fg_126[k];

        t_82[k] = f_0 * fg_127[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, t_88, t_89, fg_128, fg_129, fg_130, \
                         fg_131, fg_132, fg_133, fg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_0 * fg_128[k];

        t_84[k] = f_0 * fg_129[k];

        t_85[k] = f_0 * fg_130[k];

        t_86[k] = f_0 * fg_131[k];

        t_87[k] = f_0 * fg_132[k];

        t_88[k] = f_0 * fg_133[k];

        t_89[k] = f_0 * fg_134[k];
    }
}

auto
compute_prim_geom_10_dg_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t pg, const size_t fg,
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
    auto *t_63 = buffer.data(target + 63);
    auto *t_64 = buffer.data(target + 64);
    auto *t_65 = buffer.data(target + 65);
    auto *t_66 = buffer.data(target + 66);
    auto *t_67 = buffer.data(target + 67);
    auto *t_68 = buffer.data(target + 68);
    auto *t_69 = buffer.data(target + 69);
    auto *t_70 = buffer.data(target + 70);
    auto *t_71 = buffer.data(target + 71);
    auto *t_72 = buffer.data(target + 72);
    auto *t_73 = buffer.data(target + 73);
    auto *t_74 = buffer.data(target + 74);
    auto *t_75 = buffer.data(target + 75);
    auto *t_76 = buffer.data(target + 76);
    auto *t_77 = buffer.data(target + 77);
    auto *t_78 = buffer.data(target + 78);
    auto *t_79 = buffer.data(target + 79);
    auto *t_80 = buffer.data(target + 80);
    auto *t_81 = buffer.data(target + 81);
    auto *t_82 = buffer.data(target + 82);
    auto *t_83 = buffer.data(target + 83);
    auto *t_84 = buffer.data(target + 84);
    auto *t_85 = buffer.data(target + 85);
    auto *t_86 = buffer.data(target + 86);
    auto *t_87 = buffer.data(target + 87);
    auto *t_88 = buffer.data(target + 88);
    auto *t_89 = buffer.data(target + 89);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_11 = buffer.data(pg + 11);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_14 = buffer.data(pg + 14);
    const auto *pg_15 = buffer.data(pg + 15);
    const auto *pg_16 = buffer.data(pg + 16);
    const auto *pg_17 = buffer.data(pg + 17);
    const auto *pg_18 = buffer.data(pg + 18);
    const auto *pg_19 = buffer.data(pg + 19);
    const auto *pg_20 = buffer.data(pg + 20);
    const auto *pg_21 = buffer.data(pg + 21);
    const auto *pg_22 = buffer.data(pg + 22);
    const auto *pg_23 = buffer.data(pg + 23);
    const auto *pg_24 = buffer.data(pg + 24);
    const auto *pg_25 = buffer.data(pg + 25);
    const auto *pg_26 = buffer.data(pg + 26);
    const auto *pg_27 = buffer.data(pg + 27);
    const auto *pg_28 = buffer.data(pg + 28);
    const auto *pg_29 = buffer.data(pg + 29);
    const auto *pg_30 = buffer.data(pg + 30);
    const auto *pg_31 = buffer.data(pg + 31);
    const auto *pg_32 = buffer.data(pg + 32);
    const auto *pg_33 = buffer.data(pg + 33);
    const auto *pg_34 = buffer.data(pg + 34);
    const auto *pg_35 = buffer.data(pg + 35);
    const auto *pg_36 = buffer.data(pg + 36);
    const auto *pg_37 = buffer.data(pg + 37);
    const auto *pg_38 = buffer.data(pg + 38);
    const auto *pg_39 = buffer.data(pg + 39);
    const auto *pg_40 = buffer.data(pg + 40);
    const auto *pg_41 = buffer.data(pg + 41);
    const auto *pg_42 = buffer.data(pg + 42);
    const auto *pg_43 = buffer.data(pg + 43);
    const auto *pg_44 = buffer.data(pg + 44);

    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_31 = buffer.data(fg + 31);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_33 = buffer.data(fg + 33);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_36 = buffer.data(fg + 36);
    const auto *fg_37 = buffer.data(fg + 37);
    const auto *fg_38 = buffer.data(fg + 38);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_40 = buffer.data(fg + 40);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_60 = buffer.data(fg + 60);
    const auto *fg_61 = buffer.data(fg + 61);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_63 = buffer.data(fg + 63);
    const auto *fg_64 = buffer.data(fg + 64);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_66 = buffer.data(fg + 66);
    const auto *fg_67 = buffer.data(fg + 67);
    const auto *fg_68 = buffer.data(fg + 68);
    const auto *fg_69 = buffer.data(fg + 69);
    const auto *fg_70 = buffer.data(fg + 70);
    const auto *fg_71 = buffer.data(fg + 71);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_73 = buffer.data(fg + 73);
    const auto *fg_74 = buffer.data(fg + 74);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_76 = buffer.data(fg + 76);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_78 = buffer.data(fg + 78);
    const auto *fg_79 = buffer.data(fg + 79);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_81 = buffer.data(fg + 81);
    const auto *fg_82 = buffer.data(fg + 82);
    const auto *fg_83 = buffer.data(fg + 83);
    const auto *fg_84 = buffer.data(fg + 84);
    const auto *fg_85 = buffer.data(fg + 85);
    const auto *fg_86 = buffer.data(fg + 86);
    const auto *fg_87 = buffer.data(fg + 87);
    const auto *fg_88 = buffer.data(fg + 88);
    const auto *fg_89 = buffer.data(fg + 89);
    const auto *fg_105 = buffer.data(fg + 105);
    const auto *fg_106 = buffer.data(fg + 106);
    const auto *fg_107 = buffer.data(fg + 107);
    const auto *fg_108 = buffer.data(fg + 108);
    const auto *fg_109 = buffer.data(fg + 109);
    const auto *fg_110 = buffer.data(fg + 110);
    const auto *fg_111 = buffer.data(fg + 111);
    const auto *fg_112 = buffer.data(fg + 112);
    const auto *fg_113 = buffer.data(fg + 113);
    const auto *fg_114 = buffer.data(fg + 114);
    const auto *fg_115 = buffer.data(fg + 115);
    const auto *fg_116 = buffer.data(fg + 116);
    const auto *fg_117 = buffer.data(fg + 117);
    const auto *fg_118 = buffer.data(fg + 118);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_120 = buffer.data(fg + 120);
    const auto *fg_121 = buffer.data(fg + 121);
    const auto *fg_122 = buffer.data(fg + 122);
    const auto *fg_123 = buffer.data(fg + 123);
    const auto *fg_124 = buffer.data(fg + 124);
    const auto *fg_125 = buffer.data(fg + 125);
    const auto *fg_126 = buffer.data(fg + 126);
    const auto *fg_127 = buffer.data(fg + 127);
    const auto *fg_128 = buffer.data(fg + 128);
    const auto *fg_129 = buffer.data(fg + 129);
    const auto *fg_130 = buffer.data(fg + 130);
    const auto *fg_131 = buffer.data(fg + 131);
    const auto *fg_132 = buffer.data(fg + 132);
    const auto *fg_133 = buffer.data(fg + 133);
    const auto *fg_134 = buffer.data(fg + 134);
    const auto *fg_135 = buffer.data(fg + 135);
    const auto *fg_136 = buffer.data(fg + 136);
    const auto *fg_137 = buffer.data(fg + 137);
    const auto *fg_138 = buffer.data(fg + 138);
    const auto *fg_139 = buffer.data(fg + 139);
    const auto *fg_140 = buffer.data(fg + 140);
    const auto *fg_141 = buffer.data(fg + 141);
    const auto *fg_142 = buffer.data(fg + 142);
    const auto *fg_143 = buffer.data(fg + 143);
    const auto *fg_144 = buffer.data(fg + 144);
    const auto *fg_145 = buffer.data(fg + 145);
    const auto *fg_146 = buffer.data(fg + 146);
    const auto *fg_147 = buffer.data(fg + 147);
    const auto *fg_148 = buffer.data(fg + 148);
    const auto *fg_149 = buffer.data(fg + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, fg_30, fg_31, fg_32, fg_33, \
                         fg_34, fg_35, fg_36, fg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_30[k];

        t_1[k] = f_0 * fg_31[k];

        t_2[k] = f_0 * fg_32[k];

        t_3[k] = f_0 * fg_33[k];

        t_4[k] = f_0 * fg_34[k];

        t_5[k] = f_0 * fg_35[k];

        t_6[k] = f_0 * fg_36[k];

        t_7[k] = f_0 * fg_37[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, fg_38, fg_39, fg_40, \
                         fg_41, fg_42, fg_43, fg_44, fg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * fg_38[k];

        t_9[k] = f_0 * fg_39[k];

        t_10[k] = f_0 * fg_40[k];

        t_11[k] = f_0 * fg_41[k];

        t_12[k] = f_0 * fg_42[k];

        t_13[k] = f_0 * fg_43[k];

        t_14[k] = f_0 * fg_44[k];

        t_15[k] = f_0 * fg_60[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, fg_61, fg_62, fg_63, \
                         fg_64, fg_65, fg_66, fg_67, fg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * fg_61[k];

        t_17[k] = f_0 * fg_62[k];

        t_18[k] = f_0 * fg_63[k];

        t_19[k] = f_0 * fg_64[k];

        t_20[k] = f_0 * fg_65[k];

        t_21[k] = f_0 * fg_66[k];

        t_22[k] = f_0 * fg_67[k];

        t_23[k] = f_0 * fg_68[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, pg_0, fg_69, fg_70, fg_71, \
                         fg_72, fg_73, fg_74, fg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * fg_69[k];

        t_25[k] = f_0 * fg_70[k];

        t_26[k] = f_0 * fg_71[k];

        t_27[k] = f_0 * fg_72[k];

        t_28[k] = f_0 * fg_73[k];

        t_29[k] = f_0 * fg_74[k];

        t_30[k] = -pg_0[k]
                  + f_0 * fg_75[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pg_1, pg_2, pg_3, pg_4, pg_5, fg_76, \
                         fg_77, fg_78, fg_79, fg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -pg_1[k]
                  + f_0 * fg_76[k];

        t_32[k] = -pg_2[k]
                  + f_0 * fg_77[k];

        t_33[k] = -pg_3[k]
                  + f_0 * fg_78[k];

        t_34[k] = -pg_4[k]
                  + f_0 * fg_79[k];

        t_35[k] = -pg_5[k]
                  + f_0 * fg_80[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pg_6, pg_7, pg_8, pg_9, pg_10, fg_81, \
                         fg_82, fg_83, fg_84, fg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -pg_6[k]
                  + f_0 * fg_81[k];

        t_37[k] = -pg_7[k]
                  + f_0 * fg_82[k];

        t_38[k] = -pg_8[k]
                  + f_0 * fg_83[k];

        t_39[k] = -pg_9[k]
                  + f_0 * fg_84[k];

        t_40[k] = -pg_10[k]
                  + f_0 * fg_85[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, t_46, pg_11, pg_12, pg_13, pg_14, \
                         fg_86, fg_87, fg_88, fg_89, fg_105, fg_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -pg_11[k]
                  + f_0 * fg_86[k];

        t_42[k] = -pg_12[k]
                  + f_0 * fg_87[k];

        t_43[k] = -pg_13[k]
                  + f_0 * fg_88[k];

        t_44[k] = -pg_14[k]
                  + f_0 * fg_89[k];

        t_45[k] = f_0 * fg_105[k];

        t_46[k] = f_0 * fg_106[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, t_53, t_54, fg_107, fg_108, \
                         fg_109, fg_110, fg_111, fg_112, fg_113, \
                         fg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * fg_107[k];

        t_48[k] = f_0 * fg_108[k];

        t_49[k] = f_0 * fg_109[k];

        t_50[k] = f_0 * fg_110[k];

        t_51[k] = f_0 * fg_111[k];

        t_52[k] = f_0 * fg_112[k];

        t_53[k] = f_0 * fg_113[k];

        t_54[k] = f_0 * fg_114[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, pg_15, pg_16, fg_115, \
                         fg_116, fg_117, fg_118, fg_119, fg_120, \
                         fg_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_0 * fg_115[k];

        t_56[k] = f_0 * fg_116[k];

        t_57[k] = f_0 * fg_117[k];

        t_58[k] = f_0 * fg_118[k];

        t_59[k] = f_0 * fg_119[k];

        t_60[k] = -pg_15[k]
                  + f_0 * fg_120[k];

        t_61[k] = -pg_16[k]
                  + f_0 * fg_121[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, pg_17, pg_18, pg_19, pg_20, pg_21, \
                         fg_122, fg_123, fg_124, fg_125, fg_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -pg_17[k]
                  + f_0 * fg_122[k];

        t_63[k] = -pg_18[k]
                  + f_0 * fg_123[k];

        t_64[k] = -pg_19[k]
                  + f_0 * fg_124[k];

        t_65[k] = -pg_20[k]
                  + f_0 * fg_125[k];

        t_66[k] = -pg_21[k]
                  + f_0 * fg_126[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, pg_22, pg_23, pg_24, pg_25, pg_26, \
                         fg_127, fg_128, fg_129, fg_130, fg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -pg_22[k]
                  + f_0 * fg_127[k];

        t_68[k] = -pg_23[k]
                  + f_0 * fg_128[k];

        t_69[k] = -pg_24[k]
                  + f_0 * fg_129[k];

        t_70[k] = -pg_25[k]
                  + f_0 * fg_130[k];

        t_71[k] = -pg_26[k]
                  + f_0 * fg_131[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pg_27, pg_28, pg_29, pg_30, pg_31, \
                         fg_132, fg_133, fg_134, fg_135, fg_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -pg_27[k]
                  + f_0 * fg_132[k];

        t_73[k] = -pg_28[k]
                  + f_0 * fg_133[k];

        t_74[k] = -pg_29[k]
                  + f_0 * fg_134[k];

        t_75[k] = -2.0 * pg_30[k]
                  + f_0 * fg_135[k];

        t_76[k] = -2.0 * pg_31[k]
                  + f_0 * fg_136[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, pg_32, pg_33, pg_34, pg_35, pg_36, \
                         fg_137, fg_138, fg_139, fg_140, fg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -2.0 * pg_32[k]
                  + f_0 * fg_137[k];

        t_78[k] = -2.0 * pg_33[k]
                  + f_0 * fg_138[k];

        t_79[k] = -2.0 * pg_34[k]
                  + f_0 * fg_139[k];

        t_80[k] = -2.0 * pg_35[k]
                  + f_0 * fg_140[k];

        t_81[k] = -2.0 * pg_36[k]
                  + f_0 * fg_141[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pg_37, pg_38, pg_39, pg_40, pg_41, \
                         fg_142, fg_143, fg_144, fg_145, fg_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -2.0 * pg_37[k]
                  + f_0 * fg_142[k];

        t_83[k] = -2.0 * pg_38[k]
                  + f_0 * fg_143[k];

        t_84[k] = -2.0 * pg_39[k]
                  + f_0 * fg_144[k];

        t_85[k] = -2.0 * pg_40[k]
                  + f_0 * fg_145[k];

        t_86[k] = -2.0 * pg_41[k]
                  + f_0 * fg_146[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pg_42, pg_43, pg_44, fg_147, fg_148, \
                         fg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -2.0 * pg_42[k]
                  + f_0 * fg_147[k];

        t_88[k] = -2.0 * pg_43[k]
                  + f_0 * fg_148[k];

        t_89[k] = -2.0 * pg_44[k]
                  + f_0 * fg_149[k];
    }
}

}  // namespace simdt2ceri
