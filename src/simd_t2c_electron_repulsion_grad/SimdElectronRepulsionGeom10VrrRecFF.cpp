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


#include "SimdElectronRepulsionGeom10VrrRecFF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_ff_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t df, const size_t gf,
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
    auto *t_90 = buffer.data(target + 90);
    auto *t_91 = buffer.data(target + 91);
    auto *t_92 = buffer.data(target + 92);
    auto *t_93 = buffer.data(target + 93);
    auto *t_94 = buffer.data(target + 94);
    auto *t_95 = buffer.data(target + 95);
    auto *t_96 = buffer.data(target + 96);
    auto *t_97 = buffer.data(target + 97);
    auto *t_98 = buffer.data(target + 98);
    auto *t_99 = buffer.data(target + 99);

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

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_4 = buffer.data(gf + 4);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_7 = buffer.data(gf + 7);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_62 = buffer.data(gf + 62);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_64 = buffer.data(gf + 64);
    const auto *gf_65 = buffer.data(gf + 65);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_67 = buffer.data(gf + 67);
    const auto *gf_68 = buffer.data(gf + 68);
    const auto *gf_69 = buffer.data(gf + 69);
    const auto *gf_70 = buffer.data(gf + 70);
    const auto *gf_71 = buffer.data(gf + 71);
    const auto *gf_72 = buffer.data(gf + 72);
    const auto *gf_73 = buffer.data(gf + 73);
    const auto *gf_74 = buffer.data(gf + 74);
    const auto *gf_75 = buffer.data(gf + 75);
    const auto *gf_76 = buffer.data(gf + 76);
    const auto *gf_77 = buffer.data(gf + 77);
    const auto *gf_78 = buffer.data(gf + 78);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_80 = buffer.data(gf + 80);
    const auto *gf_81 = buffer.data(gf + 81);
    const auto *gf_82 = buffer.data(gf + 82);
    const auto *gf_83 = buffer.data(gf + 83);
    const auto *gf_84 = buffer.data(gf + 84);
    const auto *gf_85 = buffer.data(gf + 85);
    const auto *gf_86 = buffer.data(gf + 86);
    const auto *gf_87 = buffer.data(gf + 87);
    const auto *gf_88 = buffer.data(gf + 88);
    const auto *gf_89 = buffer.data(gf + 89);
    const auto *gf_90 = buffer.data(gf + 90);
    const auto *gf_91 = buffer.data(gf + 91);
    const auto *gf_92 = buffer.data(gf + 92);
    const auto *gf_93 = buffer.data(gf + 93);
    const auto *gf_94 = buffer.data(gf + 94);
    const auto *gf_95 = buffer.data(gf + 95);
    const auto *gf_96 = buffer.data(gf + 96);
    const auto *gf_97 = buffer.data(gf + 97);
    const auto *gf_98 = buffer.data(gf + 98);
    const auto *gf_99 = buffer.data(gf + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, df_0, df_1, df_2, df_3, df_4, gf_0, gf_1, \
                         gf_2, gf_3, gf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -3.0 * df_0[k]
                 + f_0 * gf_0[k];

        t_1[k] = -3.0 * df_1[k]
                 + f_0 * gf_1[k];

        t_2[k] = -3.0 * df_2[k]
                 + f_0 * gf_2[k];

        t_3[k] = -3.0 * df_3[k]
                 + f_0 * gf_3[k];

        t_4[k] = -3.0 * df_4[k]
                 + f_0 * gf_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, df_5, df_6, df_7, df_8, df_9, gf_5, gf_6, \
                         gf_7, gf_8, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -3.0 * df_5[k]
                 + f_0 * gf_5[k];

        t_6[k] = -3.0 * df_6[k]
                 + f_0 * gf_6[k];

        t_7[k] = -3.0 * df_7[k]
                 + f_0 * gf_7[k];

        t_8[k] = -3.0 * df_8[k]
                 + f_0 * gf_8[k];

        t_9[k] = -3.0 * df_9[k]
                 + f_0 * gf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, df_10, df_11, df_12, df_13, df_14, \
                         gf_10, gf_11, gf_12, gf_13, gf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -2.0 * df_10[k]
                  + f_0 * gf_10[k];

        t_11[k] = -2.0 * df_11[k]
                  + f_0 * gf_11[k];

        t_12[k] = -2.0 * df_12[k]
                  + f_0 * gf_12[k];

        t_13[k] = -2.0 * df_13[k]
                  + f_0 * gf_13[k];

        t_14[k] = -2.0 * df_14[k]
                  + f_0 * gf_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, df_15, df_16, df_17, df_18, df_19, \
                         gf_15, gf_16, gf_17, gf_18, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -2.0 * df_15[k]
                  + f_0 * gf_15[k];

        t_16[k] = -2.0 * df_16[k]
                  + f_0 * gf_16[k];

        t_17[k] = -2.0 * df_17[k]
                  + f_0 * gf_17[k];

        t_18[k] = -2.0 * df_18[k]
                  + f_0 * gf_18[k];

        t_19[k] = -2.0 * df_19[k]
                  + f_0 * gf_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, df_20, df_21, df_22, df_23, df_24, \
                         gf_20, gf_21, gf_22, gf_23, gf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -2.0 * df_20[k]
                  + f_0 * gf_20[k];

        t_21[k] = -2.0 * df_21[k]
                  + f_0 * gf_21[k];

        t_22[k] = -2.0 * df_22[k]
                  + f_0 * gf_22[k];

        t_23[k] = -2.0 * df_23[k]
                  + f_0 * gf_23[k];

        t_24[k] = -2.0 * df_24[k]
                  + f_0 * gf_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, df_25, df_26, df_27, df_28, df_29, \
                         gf_25, gf_26, gf_27, gf_28, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -2.0 * df_25[k]
                  + f_0 * gf_25[k];

        t_26[k] = -2.0 * df_26[k]
                  + f_0 * gf_26[k];

        t_27[k] = -2.0 * df_27[k]
                  + f_0 * gf_27[k];

        t_28[k] = -2.0 * df_28[k]
                  + f_0 * gf_28[k];

        t_29[k] = -2.0 * df_29[k]
                  + f_0 * gf_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, df_30, df_31, df_32, df_33, df_34, \
                         gf_30, gf_31, gf_32, gf_33, gf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -df_30[k]
                  + f_0 * gf_30[k];

        t_31[k] = -df_31[k]
                  + f_0 * gf_31[k];

        t_32[k] = -df_32[k]
                  + f_0 * gf_32[k];

        t_33[k] = -df_33[k]
                  + f_0 * gf_33[k];

        t_34[k] = -df_34[k]
                  + f_0 * gf_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, df_35, df_36, df_37, df_38, df_39, \
                         gf_35, gf_36, gf_37, gf_38, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -df_35[k]
                  + f_0 * gf_35[k];

        t_36[k] = -df_36[k]
                  + f_0 * gf_36[k];

        t_37[k] = -df_37[k]
                  + f_0 * gf_37[k];

        t_38[k] = -df_38[k]
                  + f_0 * gf_38[k];

        t_39[k] = -df_39[k]
                  + f_0 * gf_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, df_40, df_41, df_42, df_43, df_44, \
                         gf_40, gf_41, gf_42, gf_43, gf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -df_40[k]
                  + f_0 * gf_40[k];

        t_41[k] = -df_41[k]
                  + f_0 * gf_41[k];

        t_42[k] = -df_42[k]
                  + f_0 * gf_42[k];

        t_43[k] = -df_43[k]
                  + f_0 * gf_43[k];

        t_44[k] = -df_44[k]
                  + f_0 * gf_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, df_45, df_46, df_47, df_48, df_49, \
                         gf_45, gf_46, gf_47, gf_48, gf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -df_45[k]
                  + f_0 * gf_45[k];

        t_46[k] = -df_46[k]
                  + f_0 * gf_46[k];

        t_47[k] = -df_47[k]
                  + f_0 * gf_47[k];

        t_48[k] = -df_48[k]
                  + f_0 * gf_48[k];

        t_49[k] = -df_49[k]
                  + f_0 * gf_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, df_50, df_51, df_52, df_53, df_54, \
                         gf_50, gf_51, gf_52, gf_53, gf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -df_50[k]
                  + f_0 * gf_50[k];

        t_51[k] = -df_51[k]
                  + f_0 * gf_51[k];

        t_52[k] = -df_52[k]
                  + f_0 * gf_52[k];

        t_53[k] = -df_53[k]
                  + f_0 * gf_53[k];

        t_54[k] = -df_54[k]
                  + f_0 * gf_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, df_55, df_56, df_57, df_58, df_59, \
                         gf_55, gf_56, gf_57, gf_58, gf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -df_55[k]
                  + f_0 * gf_55[k];

        t_56[k] = -df_56[k]
                  + f_0 * gf_56[k];

        t_57[k] = -df_57[k]
                  + f_0 * gf_57[k];

        t_58[k] = -df_58[k]
                  + f_0 * gf_58[k];

        t_59[k] = -df_59[k]
                  + f_0 * gf_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, t_66, t_67, gf_60, gf_61, gf_62, \
                         gf_63, gf_64, gf_65, gf_66, gf_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * gf_60[k];

        t_61[k] = f_0 * gf_61[k];

        t_62[k] = f_0 * gf_62[k];

        t_63[k] = f_0 * gf_63[k];

        t_64[k] = f_0 * gf_64[k];

        t_65[k] = f_0 * gf_65[k];

        t_66[k] = f_0 * gf_66[k];

        t_67[k] = f_0 * gf_67[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, t_74, t_75, gf_68, gf_69, gf_70, \
                         gf_71, gf_72, gf_73, gf_74, gf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * gf_68[k];

        t_69[k] = f_0 * gf_69[k];

        t_70[k] = f_0 * gf_70[k];

        t_71[k] = f_0 * gf_71[k];

        t_72[k] = f_0 * gf_72[k];

        t_73[k] = f_0 * gf_73[k];

        t_74[k] = f_0 * gf_74[k];

        t_75[k] = f_0 * gf_75[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, t_81, t_82, t_83, gf_76, gf_77, gf_78, \
                         gf_79, gf_80, gf_81, gf_82, gf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_0 * gf_76[k];

        t_77[k] = f_0 * gf_77[k];

        t_78[k] = f_0 * gf_78[k];

        t_79[k] = f_0 * gf_79[k];

        t_80[k] = f_0 * gf_80[k];

        t_81[k] = f_0 * gf_81[k];

        t_82[k] = f_0 * gf_82[k];

        t_83[k] = f_0 * gf_83[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, t_90, t_91, gf_84, gf_85, gf_86, \
                         gf_87, gf_88, gf_89, gf_90, gf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_0 * gf_84[k];

        t_85[k] = f_0 * gf_85[k];

        t_86[k] = f_0 * gf_86[k];

        t_87[k] = f_0 * gf_87[k];

        t_88[k] = f_0 * gf_88[k];

        t_89[k] = f_0 * gf_89[k];

        t_90[k] = f_0 * gf_90[k];

        t_91[k] = f_0 * gf_91[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, t_97, t_98, t_99, gf_92, gf_93, gf_94, \
                         gf_95, gf_96, gf_97, gf_98, gf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_0 * gf_92[k];

        t_93[k] = f_0 * gf_93[k];

        t_94[k] = f_0 * gf_94[k];

        t_95[k] = f_0 * gf_95[k];

        t_96[k] = f_0 * gf_96[k];

        t_97[k] = f_0 * gf_97[k];

        t_98[k] = f_0 * gf_98[k];

        t_99[k] = f_0 * gf_99[k];
    }
}

auto
compute_prim_geom_10_ff_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t df, const size_t gf,
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
    auto *t_90 = buffer.data(target + 90);
    auto *t_91 = buffer.data(target + 91);
    auto *t_92 = buffer.data(target + 92);
    auto *t_93 = buffer.data(target + 93);
    auto *t_94 = buffer.data(target + 94);
    auto *t_95 = buffer.data(target + 95);
    auto *t_96 = buffer.data(target + 96);
    auto *t_97 = buffer.data(target + 97);
    auto *t_98 = buffer.data(target + 98);
    auto *t_99 = buffer.data(target + 99);

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

    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_12 = buffer.data(gf + 12);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_14 = buffer.data(gf + 14);
    const auto *gf_15 = buffer.data(gf + 15);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_34 = buffer.data(gf + 34);
    const auto *gf_35 = buffer.data(gf + 35);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_62 = buffer.data(gf + 62);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_64 = buffer.data(gf + 64);
    const auto *gf_65 = buffer.data(gf + 65);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_67 = buffer.data(gf + 67);
    const auto *gf_68 = buffer.data(gf + 68);
    const auto *gf_69 = buffer.data(gf + 69);
    const auto *gf_70 = buffer.data(gf + 70);
    const auto *gf_71 = buffer.data(gf + 71);
    const auto *gf_72 = buffer.data(gf + 72);
    const auto *gf_73 = buffer.data(gf + 73);
    const auto *gf_74 = buffer.data(gf + 74);
    const auto *gf_75 = buffer.data(gf + 75);
    const auto *gf_76 = buffer.data(gf + 76);
    const auto *gf_77 = buffer.data(gf + 77);
    const auto *gf_78 = buffer.data(gf + 78);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_80 = buffer.data(gf + 80);
    const auto *gf_81 = buffer.data(gf + 81);
    const auto *gf_82 = buffer.data(gf + 82);
    const auto *gf_83 = buffer.data(gf + 83);
    const auto *gf_84 = buffer.data(gf + 84);
    const auto *gf_85 = buffer.data(gf + 85);
    const auto *gf_86 = buffer.data(gf + 86);
    const auto *gf_87 = buffer.data(gf + 87);
    const auto *gf_88 = buffer.data(gf + 88);
    const auto *gf_89 = buffer.data(gf + 89);
    const auto *gf_100 = buffer.data(gf + 100);
    const auto *gf_101 = buffer.data(gf + 101);
    const auto *gf_102 = buffer.data(gf + 102);
    const auto *gf_103 = buffer.data(gf + 103);
    const auto *gf_104 = buffer.data(gf + 104);
    const auto *gf_105 = buffer.data(gf + 105);
    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_107 = buffer.data(gf + 107);
    const auto *gf_108 = buffer.data(gf + 108);
    const auto *gf_109 = buffer.data(gf + 109);
    const auto *gf_110 = buffer.data(gf + 110);
    const auto *gf_111 = buffer.data(gf + 111);
    const auto *gf_112 = buffer.data(gf + 112);
    const auto *gf_113 = buffer.data(gf + 113);
    const auto *gf_114 = buffer.data(gf + 114);
    const auto *gf_115 = buffer.data(gf + 115);
    const auto *gf_116 = buffer.data(gf + 116);
    const auto *gf_117 = buffer.data(gf + 117);
    const auto *gf_118 = buffer.data(gf + 118);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_120 = buffer.data(gf + 120);
    const auto *gf_121 = buffer.data(gf + 121);
    const auto *gf_122 = buffer.data(gf + 122);
    const auto *gf_123 = buffer.data(gf + 123);
    const auto *gf_124 = buffer.data(gf + 124);
    const auto *gf_125 = buffer.data(gf + 125);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_127 = buffer.data(gf + 127);
    const auto *gf_128 = buffer.data(gf + 128);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_130 = buffer.data(gf + 130);
    const auto *gf_131 = buffer.data(gf + 131);
    const auto *gf_132 = buffer.data(gf + 132);
    const auto *gf_133 = buffer.data(gf + 133);
    const auto *gf_134 = buffer.data(gf + 134);
    const auto *gf_135 = buffer.data(gf + 135);
    const auto *gf_136 = buffer.data(gf + 136);
    const auto *gf_137 = buffer.data(gf + 137);
    const auto *gf_138 = buffer.data(gf + 138);
    const auto *gf_139 = buffer.data(gf + 139);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, gf_10, gf_11, gf_12, gf_13, \
                         gf_14, gf_15, gf_16, gf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_10[k];

        t_1[k] = f_0 * gf_11[k];

        t_2[k] = f_0 * gf_12[k];

        t_3[k] = f_0 * gf_13[k];

        t_4[k] = f_0 * gf_14[k];

        t_5[k] = f_0 * gf_15[k];

        t_6[k] = f_0 * gf_16[k];

        t_7[k] = f_0 * gf_17[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, df_0, df_1, df_2, df_3, gf_18, \
                         gf_19, gf_30, gf_31, gf_32, gf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * gf_18[k];

        t_9[k] = f_0 * gf_19[k];

        t_10[k] = -df_0[k]
                  + f_0 * gf_30[k];

        t_11[k] = -df_1[k]
                  + f_0 * gf_31[k];

        t_12[k] = -df_2[k]
                  + f_0 * gf_32[k];

        t_13[k] = -df_3[k]
                  + f_0 * gf_33[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, df_4, df_5, df_6, df_7, df_8, gf_34, \
                         gf_35, gf_36, gf_37, gf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -df_4[k]
                  + f_0 * gf_34[k];

        t_15[k] = -df_5[k]
                  + f_0 * gf_35[k];

        t_16[k] = -df_6[k]
                  + f_0 * gf_36[k];

        t_17[k] = -df_7[k]
                  + f_0 * gf_37[k];

        t_18[k] = -df_8[k]
                  + f_0 * gf_38[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, t_25, df_9, gf_39, gf_40, gf_41, \
                         gf_42, gf_43, gf_44, gf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -df_9[k]
                  + f_0 * gf_39[k];

        t_20[k] = f_0 * gf_40[k];

        t_21[k] = f_0 * gf_41[k];

        t_22[k] = f_0 * gf_42[k];

        t_23[k] = f_0 * gf_43[k];

        t_24[k] = f_0 * gf_44[k];

        t_25[k] = f_0 * gf_45[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, df_10, df_11, gf_46, gf_47, \
                         gf_48, gf_49, gf_60, gf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * gf_46[k];

        t_27[k] = f_0 * gf_47[k];

        t_28[k] = f_0 * gf_48[k];

        t_29[k] = f_0 * gf_49[k];

        t_30[k] = -2.0 * df_10[k]
                  + f_0 * gf_60[k];

        t_31[k] = -2.0 * df_11[k]
                  + f_0 * gf_61[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, df_12, df_13, df_14, df_15, df_16, \
                         gf_62, gf_63, gf_64, gf_65, gf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -2.0 * df_12[k]
                  + f_0 * gf_62[k];

        t_33[k] = -2.0 * df_13[k]
                  + f_0 * gf_63[k];

        t_34[k] = -2.0 * df_14[k]
                  + f_0 * gf_64[k];

        t_35[k] = -2.0 * df_15[k]
                  + f_0 * gf_65[k];

        t_36[k] = -2.0 * df_16[k]
                  + f_0 * gf_66[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, df_17, df_18, df_19, df_20, df_21, \
                         gf_67, gf_68, gf_69, gf_70, gf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -2.0 * df_17[k]
                  + f_0 * gf_67[k];

        t_38[k] = -2.0 * df_18[k]
                  + f_0 * gf_68[k];

        t_39[k] = -2.0 * df_19[k]
                  + f_0 * gf_69[k];

        t_40[k] = -df_20[k]
                  + f_0 * gf_70[k];

        t_41[k] = -df_21[k]
                  + f_0 * gf_71[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, df_22, df_23, df_24, df_25, df_26, \
                         gf_72, gf_73, gf_74, gf_75, gf_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = -df_22[k]
                  + f_0 * gf_72[k];

        t_43[k] = -df_23[k]
                  + f_0 * gf_73[k];

        t_44[k] = -df_24[k]
                  + f_0 * gf_74[k];

        t_45[k] = -df_25[k]
                  + f_0 * gf_75[k];

        t_46[k] = -df_26[k]
                  + f_0 * gf_76[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, df_27, df_28, df_29, gf_77, \
                         gf_78, gf_79, gf_80, gf_81, gf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -df_27[k]
                  + f_0 * gf_77[k];

        t_48[k] = -df_28[k]
                  + f_0 * gf_78[k];

        t_49[k] = -df_29[k]
                  + f_0 * gf_79[k];

        t_50[k] = f_0 * gf_80[k];

        t_51[k] = f_0 * gf_81[k];

        t_52[k] = f_0 * gf_82[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, t_58, t_59, gf_83, gf_84, gf_85, gf_86, \
                         gf_87, gf_88, gf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * gf_83[k];

        t_54[k] = f_0 * gf_84[k];

        t_55[k] = f_0 * gf_85[k];

        t_56[k] = f_0 * gf_86[k];

        t_57[k] = f_0 * gf_87[k];

        t_58[k] = f_0 * gf_88[k];

        t_59[k] = f_0 * gf_89[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, df_30, df_31, df_32, df_33, df_34, \
                         gf_100, gf_101, gf_102, gf_103, gf_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -3.0 * df_30[k]
                  + f_0 * gf_100[k];

        t_61[k] = -3.0 * df_31[k]
                  + f_0 * gf_101[k];

        t_62[k] = -3.0 * df_32[k]
                  + f_0 * gf_102[k];

        t_63[k] = -3.0 * df_33[k]
                  + f_0 * gf_103[k];

        t_64[k] = -3.0 * df_34[k]
                  + f_0 * gf_104[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, df_35, df_36, df_37, df_38, df_39, \
                         gf_105, gf_106, gf_107, gf_108, gf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -3.0 * df_35[k]
                  + f_0 * gf_105[k];

        t_66[k] = -3.0 * df_36[k]
                  + f_0 * gf_106[k];

        t_67[k] = -3.0 * df_37[k]
                  + f_0 * gf_107[k];

        t_68[k] = -3.0 * df_38[k]
                  + f_0 * gf_108[k];

        t_69[k] = -3.0 * df_39[k]
                  + f_0 * gf_109[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, df_40, df_41, df_42, df_43, df_44, \
                         gf_110, gf_111, gf_112, gf_113, gf_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -2.0 * df_40[k]
                  + f_0 * gf_110[k];

        t_71[k] = -2.0 * df_41[k]
                  + f_0 * gf_111[k];

        t_72[k] = -2.0 * df_42[k]
                  + f_0 * gf_112[k];

        t_73[k] = -2.0 * df_43[k]
                  + f_0 * gf_113[k];

        t_74[k] = -2.0 * df_44[k]
                  + f_0 * gf_114[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, df_45, df_46, df_47, df_48, df_49, \
                         gf_115, gf_116, gf_117, gf_118, gf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * df_45[k]
                  + f_0 * gf_115[k];

        t_76[k] = -2.0 * df_46[k]
                  + f_0 * gf_116[k];

        t_77[k] = -2.0 * df_47[k]
                  + f_0 * gf_117[k];

        t_78[k] = -2.0 * df_48[k]
                  + f_0 * gf_118[k];

        t_79[k] = -2.0 * df_49[k]
                  + f_0 * gf_119[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, df_50, df_51, df_52, df_53, df_54, \
                         gf_120, gf_121, gf_122, gf_123, gf_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -df_50[k]
                  + f_0 * gf_120[k];

        t_81[k] = -df_51[k]
                  + f_0 * gf_121[k];

        t_82[k] = -df_52[k]
                  + f_0 * gf_122[k];

        t_83[k] = -df_53[k]
                  + f_0 * gf_123[k];

        t_84[k] = -df_54[k]
                  + f_0 * gf_124[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, df_55, df_56, df_57, df_58, df_59, \
                         gf_125, gf_126, gf_127, gf_128, gf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -df_55[k]
                  + f_0 * gf_125[k];

        t_86[k] = -df_56[k]
                  + f_0 * gf_126[k];

        t_87[k] = -df_57[k]
                  + f_0 * gf_127[k];

        t_88[k] = -df_58[k]
                  + f_0 * gf_128[k];

        t_89[k] = -df_59[k]
                  + f_0 * gf_129[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, gf_130, gf_131, \
                         gf_132, gf_133, gf_134, gf_135, gf_136, \
                         gf_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * gf_130[k];

        t_91[k] = f_0 * gf_131[k];

        t_92[k] = f_0 * gf_132[k];

        t_93[k] = f_0 * gf_133[k];

        t_94[k] = f_0 * gf_134[k];

        t_95[k] = f_0 * gf_135[k];

        t_96[k] = f_0 * gf_136[k];

        t_97[k] = f_0 * gf_137[k];
    }

#pragma omp simd aligned(t_98, t_99, gf_138, gf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * gf_138[k];

        t_99[k] = f_0 * gf_139[k];
    }
}

auto
compute_prim_geom_10_ff_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t df, const size_t gf,
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
    auto *t_90 = buffer.data(target + 90);
    auto *t_91 = buffer.data(target + 91);
    auto *t_92 = buffer.data(target + 92);
    auto *t_93 = buffer.data(target + 93);
    auto *t_94 = buffer.data(target + 94);
    auto *t_95 = buffer.data(target + 95);
    auto *t_96 = buffer.data(target + 96);
    auto *t_97 = buffer.data(target + 97);
    auto *t_98 = buffer.data(target + 98);
    auto *t_99 = buffer.data(target + 99);

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

    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_21 = buffer.data(gf + 21);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_23 = buffer.data(gf + 23);
    const auto *gf_24 = buffer.data(gf + 24);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_26 = buffer.data(gf + 26);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_40 = buffer.data(gf + 40);
    const auto *gf_41 = buffer.data(gf + 41);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_43 = buffer.data(gf + 43);
    const auto *gf_44 = buffer.data(gf + 44);
    const auto *gf_45 = buffer.data(gf + 45);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_53 = buffer.data(gf + 53);
    const auto *gf_54 = buffer.data(gf + 54);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_70 = buffer.data(gf + 70);
    const auto *gf_71 = buffer.data(gf + 71);
    const auto *gf_72 = buffer.data(gf + 72);
    const auto *gf_73 = buffer.data(gf + 73);
    const auto *gf_74 = buffer.data(gf + 74);
    const auto *gf_75 = buffer.data(gf + 75);
    const auto *gf_76 = buffer.data(gf + 76);
    const auto *gf_77 = buffer.data(gf + 77);
    const auto *gf_78 = buffer.data(gf + 78);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_80 = buffer.data(gf + 80);
    const auto *gf_81 = buffer.data(gf + 81);
    const auto *gf_82 = buffer.data(gf + 82);
    const auto *gf_83 = buffer.data(gf + 83);
    const auto *gf_84 = buffer.data(gf + 84);
    const auto *gf_85 = buffer.data(gf + 85);
    const auto *gf_86 = buffer.data(gf + 86);
    const auto *gf_87 = buffer.data(gf + 87);
    const auto *gf_88 = buffer.data(gf + 88);
    const auto *gf_89 = buffer.data(gf + 89);
    const auto *gf_90 = buffer.data(gf + 90);
    const auto *gf_91 = buffer.data(gf + 91);
    const auto *gf_92 = buffer.data(gf + 92);
    const auto *gf_93 = buffer.data(gf + 93);
    const auto *gf_94 = buffer.data(gf + 94);
    const auto *gf_95 = buffer.data(gf + 95);
    const auto *gf_96 = buffer.data(gf + 96);
    const auto *gf_97 = buffer.data(gf + 97);
    const auto *gf_98 = buffer.data(gf + 98);
    const auto *gf_99 = buffer.data(gf + 99);
    const auto *gf_110 = buffer.data(gf + 110);
    const auto *gf_111 = buffer.data(gf + 111);
    const auto *gf_112 = buffer.data(gf + 112);
    const auto *gf_113 = buffer.data(gf + 113);
    const auto *gf_114 = buffer.data(gf + 114);
    const auto *gf_115 = buffer.data(gf + 115);
    const auto *gf_116 = buffer.data(gf + 116);
    const auto *gf_117 = buffer.data(gf + 117);
    const auto *gf_118 = buffer.data(gf + 118);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_120 = buffer.data(gf + 120);
    const auto *gf_121 = buffer.data(gf + 121);
    const auto *gf_122 = buffer.data(gf + 122);
    const auto *gf_123 = buffer.data(gf + 123);
    const auto *gf_124 = buffer.data(gf + 124);
    const auto *gf_125 = buffer.data(gf + 125);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_127 = buffer.data(gf + 127);
    const auto *gf_128 = buffer.data(gf + 128);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_130 = buffer.data(gf + 130);
    const auto *gf_131 = buffer.data(gf + 131);
    const auto *gf_132 = buffer.data(gf + 132);
    const auto *gf_133 = buffer.data(gf + 133);
    const auto *gf_134 = buffer.data(gf + 134);
    const auto *gf_135 = buffer.data(gf + 135);
    const auto *gf_136 = buffer.data(gf + 136);
    const auto *gf_137 = buffer.data(gf + 137);
    const auto *gf_138 = buffer.data(gf + 138);
    const auto *gf_139 = buffer.data(gf + 139);
    const auto *gf_140 = buffer.data(gf + 140);
    const auto *gf_141 = buffer.data(gf + 141);
    const auto *gf_142 = buffer.data(gf + 142);
    const auto *gf_143 = buffer.data(gf + 143);
    const auto *gf_144 = buffer.data(gf + 144);
    const auto *gf_145 = buffer.data(gf + 145);
    const auto *gf_146 = buffer.data(gf + 146);
    const auto *gf_147 = buffer.data(gf + 147);
    const auto *gf_148 = buffer.data(gf + 148);
    const auto *gf_149 = buffer.data(gf + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, gf_20, gf_21, gf_22, gf_23, \
                         gf_24, gf_25, gf_26, gf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_20[k];

        t_1[k] = f_0 * gf_21[k];

        t_2[k] = f_0 * gf_22[k];

        t_3[k] = f_0 * gf_23[k];

        t_4[k] = f_0 * gf_24[k];

        t_5[k] = f_0 * gf_25[k];

        t_6[k] = f_0 * gf_26[k];

        t_7[k] = f_0 * gf_27[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, gf_28, gf_29, gf_40, \
                         gf_41, gf_42, gf_43, gf_44, gf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * gf_28[k];

        t_9[k] = f_0 * gf_29[k];

        t_10[k] = f_0 * gf_40[k];

        t_11[k] = f_0 * gf_41[k];

        t_12[k] = f_0 * gf_42[k];

        t_13[k] = f_0 * gf_43[k];

        t_14[k] = f_0 * gf_44[k];

        t_15[k] = f_0 * gf_45[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, df_0, df_1, gf_46, gf_47, gf_48, \
                         gf_49, gf_50, gf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * gf_46[k];

        t_17[k] = f_0 * gf_47[k];

        t_18[k] = f_0 * gf_48[k];

        t_19[k] = f_0 * gf_49[k];

        t_20[k] = -df_0[k]
                  + f_0 * gf_50[k];

        t_21[k] = -df_1[k]
                  + f_0 * gf_51[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, df_2, df_3, df_4, df_5, df_6, gf_52, \
                         gf_53, gf_54, gf_55, gf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -df_2[k]
                  + f_0 * gf_52[k];

        t_23[k] = -df_3[k]
                  + f_0 * gf_53[k];

        t_24[k] = -df_4[k]
                  + f_0 * gf_54[k];

        t_25[k] = -df_5[k]
                  + f_0 * gf_55[k];

        t_26[k] = -df_6[k]
                  + f_0 * gf_56[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, t_32, df_7, df_8, df_9, gf_57, gf_58, \
                         gf_59, gf_70, gf_71, gf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -df_7[k]
                  + f_0 * gf_57[k];

        t_28[k] = -df_8[k]
                  + f_0 * gf_58[k];

        t_29[k] = -df_9[k]
                  + f_0 * gf_59[k];

        t_30[k] = f_0 * gf_70[k];

        t_31[k] = f_0 * gf_71[k];

        t_32[k] = f_0 * gf_72[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, t_39, gf_73, gf_74, gf_75, gf_76, \
                         gf_77, gf_78, gf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * gf_73[k];

        t_34[k] = f_0 * gf_74[k];

        t_35[k] = f_0 * gf_75[k];

        t_36[k] = f_0 * gf_76[k];

        t_37[k] = f_0 * gf_77[k];

        t_38[k] = f_0 * gf_78[k];

        t_39[k] = f_0 * gf_79[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, df_10, df_11, df_12, df_13, df_14, \
                         gf_80, gf_81, gf_82, gf_83, gf_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -df_10[k]
                  + f_0 * gf_80[k];

        t_41[k] = -df_11[k]
                  + f_0 * gf_81[k];

        t_42[k] = -df_12[k]
                  + f_0 * gf_82[k];

        t_43[k] = -df_13[k]
                  + f_0 * gf_83[k];

        t_44[k] = -df_14[k]
                  + f_0 * gf_84[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, df_15, df_16, df_17, df_18, df_19, \
                         gf_85, gf_86, gf_87, gf_88, gf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -df_15[k]
                  + f_0 * gf_85[k];

        t_46[k] = -df_16[k]
                  + f_0 * gf_86[k];

        t_47[k] = -df_17[k]
                  + f_0 * gf_87[k];

        t_48[k] = -df_18[k]
                  + f_0 * gf_88[k];

        t_49[k] = -df_19[k]
                  + f_0 * gf_89[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, df_20, df_21, df_22, df_23, df_24, \
                         gf_90, gf_91, gf_92, gf_93, gf_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * df_20[k]
                  + f_0 * gf_90[k];

        t_51[k] = -2.0 * df_21[k]
                  + f_0 * gf_91[k];

        t_52[k] = -2.0 * df_22[k]
                  + f_0 * gf_92[k];

        t_53[k] = -2.0 * df_23[k]
                  + f_0 * gf_93[k];

        t_54[k] = -2.0 * df_24[k]
                  + f_0 * gf_94[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, df_25, df_26, df_27, df_28, df_29, \
                         gf_95, gf_96, gf_97, gf_98, gf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * df_25[k]
                  + f_0 * gf_95[k];

        t_56[k] = -2.0 * df_26[k]
                  + f_0 * gf_96[k];

        t_57[k] = -2.0 * df_27[k]
                  + f_0 * gf_97[k];

        t_58[k] = -2.0 * df_28[k]
                  + f_0 * gf_98[k];

        t_59[k] = -2.0 * df_29[k]
                  + f_0 * gf_99[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, t_66, t_67, gf_110, gf_111, \
                         gf_112, gf_113, gf_114, gf_115, gf_116, \
                         gf_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * gf_110[k];

        t_61[k] = f_0 * gf_111[k];

        t_62[k] = f_0 * gf_112[k];

        t_63[k] = f_0 * gf_113[k];

        t_64[k] = f_0 * gf_114[k];

        t_65[k] = f_0 * gf_115[k];

        t_66[k] = f_0 * gf_116[k];

        t_67[k] = f_0 * gf_117[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, df_30, df_31, df_32, df_33, \
                         gf_118, gf_119, gf_120, gf_121, gf_122, \
                         gf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * gf_118[k];

        t_69[k] = f_0 * gf_119[k];

        t_70[k] = -df_30[k]
                  + f_0 * gf_120[k];

        t_71[k] = -df_31[k]
                  + f_0 * gf_121[k];

        t_72[k] = -df_32[k]
                  + f_0 * gf_122[k];

        t_73[k] = -df_33[k]
                  + f_0 * gf_123[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, df_34, df_35, df_36, df_37, df_38, \
                         gf_124, gf_125, gf_126, gf_127, gf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = -df_34[k]
                  + f_0 * gf_124[k];

        t_75[k] = -df_35[k]
                  + f_0 * gf_125[k];

        t_76[k] = -df_36[k]
                  + f_0 * gf_126[k];

        t_77[k] = -df_37[k]
                  + f_0 * gf_127[k];

        t_78[k] = -df_38[k]
                  + f_0 * gf_128[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, df_39, df_40, df_41, df_42, df_43, \
                         gf_129, gf_130, gf_131, gf_132, gf_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = -df_39[k]
                  + f_0 * gf_129[k];

        t_80[k] = -2.0 * df_40[k]
                  + f_0 * gf_130[k];

        t_81[k] = -2.0 * df_41[k]
                  + f_0 * gf_131[k];

        t_82[k] = -2.0 * df_42[k]
                  + f_0 * gf_132[k];

        t_83[k] = -2.0 * df_43[k]
                  + f_0 * gf_133[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, df_44, df_45, df_46, df_47, df_48, \
                         gf_134, gf_135, gf_136, gf_137, gf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -2.0 * df_44[k]
                  + f_0 * gf_134[k];

        t_85[k] = -2.0 * df_45[k]
                  + f_0 * gf_135[k];

        t_86[k] = -2.0 * df_46[k]
                  + f_0 * gf_136[k];

        t_87[k] = -2.0 * df_47[k]
                  + f_0 * gf_137[k];

        t_88[k] = -2.0 * df_48[k]
                  + f_0 * gf_138[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, df_49, df_50, df_51, df_52, df_53, \
                         gf_139, gf_140, gf_141, gf_142, gf_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -2.0 * df_49[k]
                  + f_0 * gf_139[k];

        t_90[k] = -3.0 * df_50[k]
                  + f_0 * gf_140[k];

        t_91[k] = -3.0 * df_51[k]
                  + f_0 * gf_141[k];

        t_92[k] = -3.0 * df_52[k]
                  + f_0 * gf_142[k];

        t_93[k] = -3.0 * df_53[k]
                  + f_0 * gf_143[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, df_54, df_55, df_56, df_57, df_58, \
                         gf_144, gf_145, gf_146, gf_147, gf_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -3.0 * df_54[k]
                  + f_0 * gf_144[k];

        t_95[k] = -3.0 * df_55[k]
                  + f_0 * gf_145[k];

        t_96[k] = -3.0 * df_56[k]
                  + f_0 * gf_146[k];

        t_97[k] = -3.0 * df_57[k]
                  + f_0 * gf_147[k];

        t_98[k] = -3.0 * df_58[k]
                  + f_0 * gf_148[k];
    }

#pragma omp simd aligned(t_99, df_59, gf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -3.0 * df_59[k]
                  + f_0 * gf_149[k];
    }
}

}  // namespace simdt2ceri
