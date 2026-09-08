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


#include "SimdTransferFF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_ff(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t df, const size_t dg, const size_t nmax) -> void
{
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, df_0, df_1, df_2, df_3, df_4, dg_0, \
                         dg_1, dg_2, dg_3, dg_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * df_0[k]
                 + dg_0[k];

        t_1[k] = -ab_x[k] * df_1[k]
                 + dg_1[k];

        t_2[k] = -ab_x[k] * df_2[k]
                 + dg_2[k];

        t_3[k] = -ab_x[k] * df_3[k]
                 + dg_3[k];

        t_4[k] = -ab_x[k] * df_4[k]
                 + dg_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, df_5, df_6, df_7, df_8, df_9, dg_5, \
                         dg_6, dg_7, dg_8, dg_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * df_5[k]
                 + dg_5[k];

        t_6[k] = -ab_x[k] * df_6[k]
                 + dg_6[k];

        t_7[k] = -ab_x[k] * df_7[k]
                 + dg_7[k];

        t_8[k] = -ab_x[k] * df_8[k]
                 + dg_8[k];

        t_9[k] = -ab_x[k] * df_9[k]
                 + dg_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, df_10, df_11, df_12, df_13, \
                         df_14, dg_15, dg_16, dg_17, dg_18, dg_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * df_10[k]
                  + dg_15[k];

        t_11[k] = -ab_x[k] * df_11[k]
                  + dg_16[k];

        t_12[k] = -ab_x[k] * df_12[k]
                  + dg_17[k];

        t_13[k] = -ab_x[k] * df_13[k]
                  + dg_18[k];

        t_14[k] = -ab_x[k] * df_14[k]
                  + dg_19[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, df_15, df_16, df_17, df_18, \
                         df_19, dg_20, dg_21, dg_22, dg_23, dg_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * df_15[k]
                  + dg_20[k];

        t_16[k] = -ab_x[k] * df_16[k]
                  + dg_21[k];

        t_17[k] = -ab_x[k] * df_17[k]
                  + dg_22[k];

        t_18[k] = -ab_x[k] * df_18[k]
                  + dg_23[k];

        t_19[k] = -ab_x[k] * df_19[k]
                  + dg_24[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, df_20, df_21, df_22, df_23, \
                         df_24, dg_30, dg_31, dg_32, dg_33, dg_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * df_20[k]
                  + dg_30[k];

        t_21[k] = -ab_x[k] * df_21[k]
                  + dg_31[k];

        t_22[k] = -ab_x[k] * df_22[k]
                  + dg_32[k];

        t_23[k] = -ab_x[k] * df_23[k]
                  + dg_33[k];

        t_24[k] = -ab_x[k] * df_24[k]
                  + dg_34[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, df_25, df_26, df_27, df_28, \
                         df_29, dg_35, dg_36, dg_37, dg_38, dg_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * df_25[k]
                  + dg_35[k];

        t_26[k] = -ab_x[k] * df_26[k]
                  + dg_36[k];

        t_27[k] = -ab_x[k] * df_27[k]
                  + dg_37[k];

        t_28[k] = -ab_x[k] * df_28[k]
                  + dg_38[k];

        t_29[k] = -ab_x[k] * df_29[k]
                  + dg_39[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, df_30, df_31, df_32, df_33, \
                         df_34, dg_45, dg_46, dg_47, dg_48, dg_49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_x[k] * df_30[k]
                  + dg_45[k];

        t_31[k] = -ab_x[k] * df_31[k]
                  + dg_46[k];

        t_32[k] = -ab_x[k] * df_32[k]
                  + dg_47[k];

        t_33[k] = -ab_x[k] * df_33[k]
                  + dg_48[k];

        t_34[k] = -ab_x[k] * df_34[k]
                  + dg_49[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, df_35, df_36, df_37, df_38, \
                         df_39, dg_50, dg_51, dg_52, dg_53, dg_54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_x[k] * df_35[k]
                  + dg_50[k];

        t_36[k] = -ab_x[k] * df_36[k]
                  + dg_51[k];

        t_37[k] = -ab_x[k] * df_37[k]
                  + dg_52[k];

        t_38[k] = -ab_x[k] * df_38[k]
                  + dg_53[k];

        t_39[k] = -ab_x[k] * df_39[k]
                  + dg_54[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, df_40, df_41, df_42, df_43, \
                         df_44, dg_60, dg_61, dg_62, dg_63, dg_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_x[k] * df_40[k]
                  + dg_60[k];

        t_41[k] = -ab_x[k] * df_41[k]
                  + dg_61[k];

        t_42[k] = -ab_x[k] * df_42[k]
                  + dg_62[k];

        t_43[k] = -ab_x[k] * df_43[k]
                  + dg_63[k];

        t_44[k] = -ab_x[k] * df_44[k]
                  + dg_64[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, df_45, df_46, df_47, df_48, \
                         df_49, dg_65, dg_66, dg_67, dg_68, dg_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = -ab_x[k] * df_45[k]
                  + dg_65[k];

        t_46[k] = -ab_x[k] * df_46[k]
                  + dg_66[k];

        t_47[k] = -ab_x[k] * df_47[k]
                  + dg_67[k];

        t_48[k] = -ab_x[k] * df_48[k]
                  + dg_68[k];

        t_49[k] = -ab_x[k] * df_49[k]
                  + dg_69[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, df_50, df_51, df_52, df_53, \
                         df_54, dg_75, dg_76, dg_77, dg_78, dg_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = -ab_x[k] * df_50[k]
                  + dg_75[k];

        t_51[k] = -ab_x[k] * df_51[k]
                  + dg_76[k];

        t_52[k] = -ab_x[k] * df_52[k]
                  + dg_77[k];

        t_53[k] = -ab_x[k] * df_53[k]
                  + dg_78[k];

        t_54[k] = -ab_x[k] * df_54[k]
                  + dg_79[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, df_55, df_56, df_57, df_58, \
                         df_59, dg_80, dg_81, dg_82, dg_83, dg_84 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = -ab_x[k] * df_55[k]
                  + dg_80[k];

        t_56[k] = -ab_x[k] * df_56[k]
                  + dg_81[k];

        t_57[k] = -ab_x[k] * df_57[k]
                  + dg_82[k];

        t_58[k] = -ab_x[k] * df_58[k]
                  + dg_83[k];

        t_59[k] = -ab_x[k] * df_59[k]
                  + dg_84[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_y, df_30, df_31, df_32, df_33, \
                         df_34, dg_46, dg_48, dg_49, dg_51, dg_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = -ab_y[k] * df_30[k]
                  + dg_46[k];

        t_61[k] = -ab_y[k] * df_31[k]
                  + dg_48[k];

        t_62[k] = -ab_y[k] * df_32[k]
                  + dg_49[k];

        t_63[k] = -ab_y[k] * df_33[k]
                  + dg_51[k];

        t_64[k] = -ab_y[k] * df_34[k]
                  + dg_52[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_y, df_35, df_36, df_37, df_38, \
                         df_39, dg_53, dg_55, dg_56, dg_57, dg_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = -ab_y[k] * df_35[k]
                  + dg_53[k];

        t_66[k] = -ab_y[k] * df_36[k]
                  + dg_55[k];

        t_67[k] = -ab_y[k] * df_37[k]
                  + dg_56[k];

        t_68[k] = -ab_y[k] * df_38[k]
                  + dg_57[k];

        t_69[k] = -ab_y[k] * df_39[k]
                  + dg_58[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_y, df_40, df_41, df_42, df_43, \
                         df_44, dg_61, dg_63, dg_64, dg_66, dg_67 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = -ab_y[k] * df_40[k]
                  + dg_61[k];

        t_71[k] = -ab_y[k] * df_41[k]
                  + dg_63[k];

        t_72[k] = -ab_y[k] * df_42[k]
                  + dg_64[k];

        t_73[k] = -ab_y[k] * df_43[k]
                  + dg_66[k];

        t_74[k] = -ab_y[k] * df_44[k]
                  + dg_67[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_y, df_45, df_46, df_47, df_48, \
                         df_49, dg_68, dg_70, dg_71, dg_72, dg_73 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = -ab_y[k] * df_45[k]
                  + dg_68[k];

        t_76[k] = -ab_y[k] * df_46[k]
                  + dg_70[k];

        t_77[k] = -ab_y[k] * df_47[k]
                  + dg_71[k];

        t_78[k] = -ab_y[k] * df_48[k]
                  + dg_72[k];

        t_79[k] = -ab_y[k] * df_49[k]
                  + dg_73[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_y, df_50, df_51, df_52, df_53, \
                         df_54, dg_76, dg_78, dg_79, dg_81, dg_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = -ab_y[k] * df_50[k]
                  + dg_76[k];

        t_81[k] = -ab_y[k] * df_51[k]
                  + dg_78[k];

        t_82[k] = -ab_y[k] * df_52[k]
                  + dg_79[k];

        t_83[k] = -ab_y[k] * df_53[k]
                  + dg_81[k];

        t_84[k] = -ab_y[k] * df_54[k]
                  + dg_82[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_y, df_55, df_56, df_57, df_58, \
                         df_59, dg_83, dg_85, dg_86, dg_87, dg_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = -ab_y[k] * df_55[k]
                  + dg_83[k];

        t_86[k] = -ab_y[k] * df_56[k]
                  + dg_85[k];

        t_87[k] = -ab_y[k] * df_57[k]
                  + dg_86[k];

        t_88[k] = -ab_y[k] * df_58[k]
                  + dg_87[k];

        t_89[k] = -ab_y[k] * df_59[k]
                  + dg_88[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_z, df_50, df_51, df_52, df_53, \
                         df_54, dg_77, dg_79, dg_80, dg_82, dg_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = -ab_z[k] * df_50[k]
                  + dg_77[k];

        t_91[k] = -ab_z[k] * df_51[k]
                  + dg_79[k];

        t_92[k] = -ab_z[k] * df_52[k]
                  + dg_80[k];

        t_93[k] = -ab_z[k] * df_53[k]
                  + dg_82[k];

        t_94[k] = -ab_z[k] * df_54[k]
                  + dg_83[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_z, df_55, df_56, df_57, df_58, \
                         df_59, dg_84, dg_86, dg_87, dg_88, dg_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = -ab_z[k] * df_55[k]
                  + dg_84[k];

        t_96[k] = -ab_z[k] * df_56[k]
                  + dg_86[k];

        t_97[k] = -ab_z[k] * df_57[k]
                  + dg_87[k];

        t_98[k] = -ab_z[k] * df_58[k]
                  + dg_88[k];

        t_99[k] = -ab_z[k] * df_59[k]
                  + dg_89[k];
    }
}

}  // namespace simdtrf
