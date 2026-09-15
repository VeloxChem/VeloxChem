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


#include "SimdElectronRepulsionGeom10VrrRecSL.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_sl_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t pl, const size_t ncols,
                                             const double alpha) -> void
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

    const auto *pl_0 = buffer.data(pl + 0);
    const auto *pl_1 = buffer.data(pl + 1);
    const auto *pl_2 = buffer.data(pl + 2);
    const auto *pl_3 = buffer.data(pl + 3);
    const auto *pl_4 = buffer.data(pl + 4);
    const auto *pl_5 = buffer.data(pl + 5);
    const auto *pl_6 = buffer.data(pl + 6);
    const auto *pl_7 = buffer.data(pl + 7);
    const auto *pl_8 = buffer.data(pl + 8);
    const auto *pl_9 = buffer.data(pl + 9);
    const auto *pl_10 = buffer.data(pl + 10);
    const auto *pl_11 = buffer.data(pl + 11);
    const auto *pl_12 = buffer.data(pl + 12);
    const auto *pl_13 = buffer.data(pl + 13);
    const auto *pl_14 = buffer.data(pl + 14);
    const auto *pl_15 = buffer.data(pl + 15);
    const auto *pl_16 = buffer.data(pl + 16);
    const auto *pl_17 = buffer.data(pl + 17);
    const auto *pl_18 = buffer.data(pl + 18);
    const auto *pl_19 = buffer.data(pl + 19);
    const auto *pl_20 = buffer.data(pl + 20);
    const auto *pl_21 = buffer.data(pl + 21);
    const auto *pl_22 = buffer.data(pl + 22);
    const auto *pl_23 = buffer.data(pl + 23);
    const auto *pl_24 = buffer.data(pl + 24);
    const auto *pl_25 = buffer.data(pl + 25);
    const auto *pl_26 = buffer.data(pl + 26);
    const auto *pl_27 = buffer.data(pl + 27);
    const auto *pl_28 = buffer.data(pl + 28);
    const auto *pl_29 = buffer.data(pl + 29);
    const auto *pl_30 = buffer.data(pl + 30);
    const auto *pl_31 = buffer.data(pl + 31);
    const auto *pl_32 = buffer.data(pl + 32);
    const auto *pl_33 = buffer.data(pl + 33);
    const auto *pl_34 = buffer.data(pl + 34);
    const auto *pl_35 = buffer.data(pl + 35);
    const auto *pl_36 = buffer.data(pl + 36);
    const auto *pl_37 = buffer.data(pl + 37);
    const auto *pl_38 = buffer.data(pl + 38);
    const auto *pl_39 = buffer.data(pl + 39);
    const auto *pl_40 = buffer.data(pl + 40);
    const auto *pl_41 = buffer.data(pl + 41);
    const auto *pl_42 = buffer.data(pl + 42);
    const auto *pl_43 = buffer.data(pl + 43);
    const auto *pl_44 = buffer.data(pl + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, pl_0, pl_1, pl_2, pl_3, pl_4, \
                         pl_5, pl_6, pl_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pl_0[k];

        t_1[k] = f_0 * pl_1[k];

        t_2[k] = f_0 * pl_2[k];

        t_3[k] = f_0 * pl_3[k];

        t_4[k] = f_0 * pl_4[k];

        t_5[k] = f_0 * pl_5[k];

        t_6[k] = f_0 * pl_6[k];

        t_7[k] = f_0 * pl_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, pl_8, pl_9, pl_10, \
                         pl_11, pl_12, pl_13, pl_14, pl_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * pl_8[k];

        t_9[k] = f_0 * pl_9[k];

        t_10[k] = f_0 * pl_10[k];

        t_11[k] = f_0 * pl_11[k];

        t_12[k] = f_0 * pl_12[k];

        t_13[k] = f_0 * pl_13[k];

        t_14[k] = f_0 * pl_14[k];

        t_15[k] = f_0 * pl_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, pl_16, pl_17, pl_18, \
                         pl_19, pl_20, pl_21, pl_22, pl_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * pl_16[k];

        t_17[k] = f_0 * pl_17[k];

        t_18[k] = f_0 * pl_18[k];

        t_19[k] = f_0 * pl_19[k];

        t_20[k] = f_0 * pl_20[k];

        t_21[k] = f_0 * pl_21[k];

        t_22[k] = f_0 * pl_22[k];

        t_23[k] = f_0 * pl_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, pl_24, pl_25, pl_26, \
                         pl_27, pl_28, pl_29, pl_30, pl_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * pl_24[k];

        t_25[k] = f_0 * pl_25[k];

        t_26[k] = f_0 * pl_26[k];

        t_27[k] = f_0 * pl_27[k];

        t_28[k] = f_0 * pl_28[k];

        t_29[k] = f_0 * pl_29[k];

        t_30[k] = f_0 * pl_30[k];

        t_31[k] = f_0 * pl_31[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, pl_32, pl_33, pl_34, \
                         pl_35, pl_36, pl_37, pl_38, pl_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * pl_32[k];

        t_33[k] = f_0 * pl_33[k];

        t_34[k] = f_0 * pl_34[k];

        t_35[k] = f_0 * pl_35[k];

        t_36[k] = f_0 * pl_36[k];

        t_37[k] = f_0 * pl_37[k];

        t_38[k] = f_0 * pl_38[k];

        t_39[k] = f_0 * pl_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pl_40, pl_41, pl_42, pl_43, \
                         pl_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * pl_40[k];

        t_41[k] = f_0 * pl_41[k];

        t_42[k] = f_0 * pl_42[k];

        t_43[k] = f_0 * pl_43[k];

        t_44[k] = f_0 * pl_44[k];
    }
}

auto
compute_prim_geom_10_sl_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t pl, const size_t ncols,
                                             const double alpha) -> void
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

    const auto *pl_45 = buffer.data(pl + 45);
    const auto *pl_46 = buffer.data(pl + 46);
    const auto *pl_47 = buffer.data(pl + 47);
    const auto *pl_48 = buffer.data(pl + 48);
    const auto *pl_49 = buffer.data(pl + 49);
    const auto *pl_50 = buffer.data(pl + 50);
    const auto *pl_51 = buffer.data(pl + 51);
    const auto *pl_52 = buffer.data(pl + 52);
    const auto *pl_53 = buffer.data(pl + 53);
    const auto *pl_54 = buffer.data(pl + 54);
    const auto *pl_55 = buffer.data(pl + 55);
    const auto *pl_56 = buffer.data(pl + 56);
    const auto *pl_57 = buffer.data(pl + 57);
    const auto *pl_58 = buffer.data(pl + 58);
    const auto *pl_59 = buffer.data(pl + 59);
    const auto *pl_60 = buffer.data(pl + 60);
    const auto *pl_61 = buffer.data(pl + 61);
    const auto *pl_62 = buffer.data(pl + 62);
    const auto *pl_63 = buffer.data(pl + 63);
    const auto *pl_64 = buffer.data(pl + 64);
    const auto *pl_65 = buffer.data(pl + 65);
    const auto *pl_66 = buffer.data(pl + 66);
    const auto *pl_67 = buffer.data(pl + 67);
    const auto *pl_68 = buffer.data(pl + 68);
    const auto *pl_69 = buffer.data(pl + 69);
    const auto *pl_70 = buffer.data(pl + 70);
    const auto *pl_71 = buffer.data(pl + 71);
    const auto *pl_72 = buffer.data(pl + 72);
    const auto *pl_73 = buffer.data(pl + 73);
    const auto *pl_74 = buffer.data(pl + 74);
    const auto *pl_75 = buffer.data(pl + 75);
    const auto *pl_76 = buffer.data(pl + 76);
    const auto *pl_77 = buffer.data(pl + 77);
    const auto *pl_78 = buffer.data(pl + 78);
    const auto *pl_79 = buffer.data(pl + 79);
    const auto *pl_80 = buffer.data(pl + 80);
    const auto *pl_81 = buffer.data(pl + 81);
    const auto *pl_82 = buffer.data(pl + 82);
    const auto *pl_83 = buffer.data(pl + 83);
    const auto *pl_84 = buffer.data(pl + 84);
    const auto *pl_85 = buffer.data(pl + 85);
    const auto *pl_86 = buffer.data(pl + 86);
    const auto *pl_87 = buffer.data(pl + 87);
    const auto *pl_88 = buffer.data(pl + 88);
    const auto *pl_89 = buffer.data(pl + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, pl_45, pl_46, pl_47, pl_48, \
                         pl_49, pl_50, pl_51, pl_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pl_45[k];

        t_1[k] = f_0 * pl_46[k];

        t_2[k] = f_0 * pl_47[k];

        t_3[k] = f_0 * pl_48[k];

        t_4[k] = f_0 * pl_49[k];

        t_5[k] = f_0 * pl_50[k];

        t_6[k] = f_0 * pl_51[k];

        t_7[k] = f_0 * pl_52[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, pl_53, pl_54, pl_55, \
                         pl_56, pl_57, pl_58, pl_59, pl_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * pl_53[k];

        t_9[k] = f_0 * pl_54[k];

        t_10[k] = f_0 * pl_55[k];

        t_11[k] = f_0 * pl_56[k];

        t_12[k] = f_0 * pl_57[k];

        t_13[k] = f_0 * pl_58[k];

        t_14[k] = f_0 * pl_59[k];

        t_15[k] = f_0 * pl_60[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, pl_61, pl_62, pl_63, \
                         pl_64, pl_65, pl_66, pl_67, pl_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * pl_61[k];

        t_17[k] = f_0 * pl_62[k];

        t_18[k] = f_0 * pl_63[k];

        t_19[k] = f_0 * pl_64[k];

        t_20[k] = f_0 * pl_65[k];

        t_21[k] = f_0 * pl_66[k];

        t_22[k] = f_0 * pl_67[k];

        t_23[k] = f_0 * pl_68[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, pl_69, pl_70, pl_71, \
                         pl_72, pl_73, pl_74, pl_75, pl_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * pl_69[k];

        t_25[k] = f_0 * pl_70[k];

        t_26[k] = f_0 * pl_71[k];

        t_27[k] = f_0 * pl_72[k];

        t_28[k] = f_0 * pl_73[k];

        t_29[k] = f_0 * pl_74[k];

        t_30[k] = f_0 * pl_75[k];

        t_31[k] = f_0 * pl_76[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, pl_77, pl_78, pl_79, \
                         pl_80, pl_81, pl_82, pl_83, pl_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * pl_77[k];

        t_33[k] = f_0 * pl_78[k];

        t_34[k] = f_0 * pl_79[k];

        t_35[k] = f_0 * pl_80[k];

        t_36[k] = f_0 * pl_81[k];

        t_37[k] = f_0 * pl_82[k];

        t_38[k] = f_0 * pl_83[k];

        t_39[k] = f_0 * pl_84[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pl_85, pl_86, pl_87, pl_88, \
                         pl_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * pl_85[k];

        t_41[k] = f_0 * pl_86[k];

        t_42[k] = f_0 * pl_87[k];

        t_43[k] = f_0 * pl_88[k];

        t_44[k] = f_0 * pl_89[k];
    }
}

auto
compute_prim_geom_10_sl_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t pl, const size_t ncols,
                                             const double alpha) -> void
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

    const auto *pl_90 = buffer.data(pl + 90);
    const auto *pl_91 = buffer.data(pl + 91);
    const auto *pl_92 = buffer.data(pl + 92);
    const auto *pl_93 = buffer.data(pl + 93);
    const auto *pl_94 = buffer.data(pl + 94);
    const auto *pl_95 = buffer.data(pl + 95);
    const auto *pl_96 = buffer.data(pl + 96);
    const auto *pl_97 = buffer.data(pl + 97);
    const auto *pl_98 = buffer.data(pl + 98);
    const auto *pl_99 = buffer.data(pl + 99);
    const auto *pl_100 = buffer.data(pl + 100);
    const auto *pl_101 = buffer.data(pl + 101);
    const auto *pl_102 = buffer.data(pl + 102);
    const auto *pl_103 = buffer.data(pl + 103);
    const auto *pl_104 = buffer.data(pl + 104);
    const auto *pl_105 = buffer.data(pl + 105);
    const auto *pl_106 = buffer.data(pl + 106);
    const auto *pl_107 = buffer.data(pl + 107);
    const auto *pl_108 = buffer.data(pl + 108);
    const auto *pl_109 = buffer.data(pl + 109);
    const auto *pl_110 = buffer.data(pl + 110);
    const auto *pl_111 = buffer.data(pl + 111);
    const auto *pl_112 = buffer.data(pl + 112);
    const auto *pl_113 = buffer.data(pl + 113);
    const auto *pl_114 = buffer.data(pl + 114);
    const auto *pl_115 = buffer.data(pl + 115);
    const auto *pl_116 = buffer.data(pl + 116);
    const auto *pl_117 = buffer.data(pl + 117);
    const auto *pl_118 = buffer.data(pl + 118);
    const auto *pl_119 = buffer.data(pl + 119);
    const auto *pl_120 = buffer.data(pl + 120);
    const auto *pl_121 = buffer.data(pl + 121);
    const auto *pl_122 = buffer.data(pl + 122);
    const auto *pl_123 = buffer.data(pl + 123);
    const auto *pl_124 = buffer.data(pl + 124);
    const auto *pl_125 = buffer.data(pl + 125);
    const auto *pl_126 = buffer.data(pl + 126);
    const auto *pl_127 = buffer.data(pl + 127);
    const auto *pl_128 = buffer.data(pl + 128);
    const auto *pl_129 = buffer.data(pl + 129);
    const auto *pl_130 = buffer.data(pl + 130);
    const auto *pl_131 = buffer.data(pl + 131);
    const auto *pl_132 = buffer.data(pl + 132);
    const auto *pl_133 = buffer.data(pl + 133);
    const auto *pl_134 = buffer.data(pl + 134);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, pl_90, pl_91, pl_92, pl_93, \
                         pl_94, pl_95, pl_96, pl_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pl_90[k];

        t_1[k] = f_0 * pl_91[k];

        t_2[k] = f_0 * pl_92[k];

        t_3[k] = f_0 * pl_93[k];

        t_4[k] = f_0 * pl_94[k];

        t_5[k] = f_0 * pl_95[k];

        t_6[k] = f_0 * pl_96[k];

        t_7[k] = f_0 * pl_97[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, pl_98, pl_99, pl_100, \
                         pl_101, pl_102, pl_103, pl_104, pl_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * pl_98[k];

        t_9[k] = f_0 * pl_99[k];

        t_10[k] = f_0 * pl_100[k];

        t_11[k] = f_0 * pl_101[k];

        t_12[k] = f_0 * pl_102[k];

        t_13[k] = f_0 * pl_103[k];

        t_14[k] = f_0 * pl_104[k];

        t_15[k] = f_0 * pl_105[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, pl_106, pl_107, \
                         pl_108, pl_109, pl_110, pl_111, pl_112, \
                         pl_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * pl_106[k];

        t_17[k] = f_0 * pl_107[k];

        t_18[k] = f_0 * pl_108[k];

        t_19[k] = f_0 * pl_109[k];

        t_20[k] = f_0 * pl_110[k];

        t_21[k] = f_0 * pl_111[k];

        t_22[k] = f_0 * pl_112[k];

        t_23[k] = f_0 * pl_113[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, pl_114, pl_115, \
                         pl_116, pl_117, pl_118, pl_119, pl_120, \
                         pl_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * pl_114[k];

        t_25[k] = f_0 * pl_115[k];

        t_26[k] = f_0 * pl_116[k];

        t_27[k] = f_0 * pl_117[k];

        t_28[k] = f_0 * pl_118[k];

        t_29[k] = f_0 * pl_119[k];

        t_30[k] = f_0 * pl_120[k];

        t_31[k] = f_0 * pl_121[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, pl_122, pl_123, \
                         pl_124, pl_125, pl_126, pl_127, pl_128, \
                         pl_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * pl_122[k];

        t_33[k] = f_0 * pl_123[k];

        t_34[k] = f_0 * pl_124[k];

        t_35[k] = f_0 * pl_125[k];

        t_36[k] = f_0 * pl_126[k];

        t_37[k] = f_0 * pl_127[k];

        t_38[k] = f_0 * pl_128[k];

        t_39[k] = f_0 * pl_129[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pl_130, pl_131, pl_132, pl_133, \
                         pl_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * pl_130[k];

        t_41[k] = f_0 * pl_131[k];

        t_42[k] = f_0 * pl_132[k];

        t_43[k] = f_0 * pl_133[k];

        t_44[k] = f_0 * pl_134[k];
    }
}

}  // namespace simdt2ceri
