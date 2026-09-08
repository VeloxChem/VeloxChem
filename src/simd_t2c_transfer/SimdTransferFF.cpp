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

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_hrr_ff_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t df, const size_t dg,
                   const size_t nmax) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.75 * std::sqrt(6.0);
    const auto f_1 = 1.25 * std::sqrt(6.0);
    const auto f_2 = 0.375 * std::sqrt(15.0);
    const auto f_3 = 1.5 * std::sqrt(15.0);
    const auto f_4 = 0.125 * std::sqrt(15.0);
    const auto f_5 = 0.5 * std::sqrt(15.0);
    const auto f_6 = 1.125 * std::sqrt(10.0);
    const auto f_7 = 0.75 * std::sqrt(10.0);
    const auto f_8 = 0.375 * std::sqrt(10.0);
    const auto f_9 = 0.25 * std::sqrt(10.0);
    const auto f_10 = 1.875 * std::sqrt(6.0);
    const auto f_11 = 0.625 * std::sqrt(6.0);
    const auto f_12 = 3.0 * std::sqrt(10.0);
    const auto f_13 = std::sqrt(15.0);
    const auto f_14 = 0.375 * std::sqrt(6.0);
    const auto f_15 = 0.25 * std::sqrt(6.0);
    const auto f_16 = 1.5 * std::sqrt(6.0);
    const auto f_17 = std::sqrt(6.0);
    const auto f_18 = 1.5 * std::sqrt(10.0);
    const auto f_19 = 0.75 * std::sqrt(15.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;
    auto *g_7 = values + 7 * nvalues;
    auto *g_8 = values + 8 * nvalues;
    auto *g_9 = values + 9 * nvalues;
    auto *g_10 = values + 10 * nvalues;
    auto *g_11 = values + 11 * nvalues;
    auto *g_12 = values + 12 * nvalues;
    auto *g_13 = values + 13 * nvalues;
    auto *g_14 = values + 14 * nvalues;
    auto *g_15 = values + 15 * nvalues;
    auto *g_16 = values + 16 * nvalues;
    auto *g_17 = values + 17 * nvalues;
    auto *g_18 = values + 18 * nvalues;
    auto *g_19 = values + 19 * nvalues;
    auto *g_20 = values + 20 * nvalues;
    auto *g_21 = values + 21 * nvalues;
    auto *g_22 = values + 22 * nvalues;
    auto *g_23 = values + 23 * nvalues;
    auto *g_24 = values + 24 * nvalues;
    auto *g_25 = values + 25 * nvalues;
    auto *g_26 = values + 26 * nvalues;
    auto *g_27 = values + 27 * nvalues;
    auto *g_28 = values + 28 * nvalues;
    auto *g_29 = values + 29 * nvalues;
    auto *g_30 = values + 30 * nvalues;
    auto *g_31 = values + 31 * nvalues;
    auto *g_32 = values + 32 * nvalues;
    auto *g_33 = values + 33 * nvalues;
    auto *g_34 = values + 34 * nvalues;
    auto *g_35 = values + 35 * nvalues;
    auto *g_36 = values + 36 * nvalues;
    auto *g_37 = values + 37 * nvalues;
    auto *g_38 = values + 38 * nvalues;
    auto *g_39 = values + 39 * nvalues;
    auto *g_40 = values + 40 * nvalues;
    auto *g_41 = values + 41 * nvalues;
    auto *g_42 = values + 42 * nvalues;
    auto *g_43 = values + 43 * nvalues;
    auto *g_44 = values + 44 * nvalues;
    auto *g_45 = values + 45 * nvalues;
    auto *g_46 = values + 46 * nvalues;
    auto *g_47 = values + 47 * nvalues;
    auto *g_48 = values + 48 * nvalues;

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

#pragma omp simd aligned(ab_x, ab_y, df_11, df_14, df_16, df_31, df_34, df_36, dg_16, dg_19, \
                         dg_21, dg_48, dg_52, dg_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -5.625 * ab_x[k] * df_11[k]
                 + 1.875 * ab_x[k] * df_16[k]
                 + 1.875 * ab_y[k] * df_31[k]
                 - 0.625 * ab_y[k] * df_36[k]
                 + 5.625 * dg_16[k]
                 - 1.875 * dg_21[k]
                 - 1.875 * dg_48[k]
                 + 0.625 * dg_55[k];

        g_1[k] = -f_0 * ab_x[k] * df_14[k]
                 + f_1 * ab_y[k] * df_34[k]
                 + f_0 * dg_19[k]
                 - f_1 * dg_52[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_11, df_16, df_18, df_31, df_36, df_38, dg_16, dg_21, \
                         dg_23, dg_48, dg_55, dg_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = f_2 * ab_x[k] * df_11[k]
                 + f_2 * ab_x[k] * df_16[k]
                 - f_3 * ab_x[k] * df_18[k]
                 - f_4 * ab_y[k] * df_31[k]
                 - f_4 * ab_y[k] * df_36[k]
                 + f_5 * ab_y[k] * df_38[k]
                 - f_2 * dg_16[k]
                 - f_2 * dg_21[k]
                 + f_3 * dg_23[k]
                 + f_4 * dg_48[k]
                 + f_4 * dg_55[k]
                 - f_5 * dg_57[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_12, df_17, df_19, df_32, df_37, df_39, dg_17, dg_22, \
                         dg_24, dg_49, dg_56, dg_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_6 * ab_x[k] * df_12[k]
                 + f_6 * ab_x[k] * df_17[k]
                 - f_7 * ab_x[k] * df_19[k]
                 - f_8 * ab_y[k] * df_32[k]
                 - f_8 * ab_y[k] * df_37[k]
                 + f_9 * ab_y[k] * df_39[k]
                 - f_6 * dg_17[k]
                 - f_6 * dg_22[k]
                 + f_7 * dg_24[k]
                 + f_8 * dg_49[k]
                 + f_8 * dg_56[k]
                 - f_9 * dg_58[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_10, df_13, df_15, df_30, df_33, df_35, dg_15, dg_18, \
                         dg_20, dg_46, dg_51, dg_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_2 * ab_x[k] * df_10[k]
                 + f_2 * ab_x[k] * df_13[k]
                 - f_3 * ab_x[k] * df_15[k]
                 - f_4 * ab_y[k] * df_30[k]
                 - f_4 * ab_y[k] * df_33[k]
                 + f_5 * ab_y[k] * df_35[k]
                 - f_2 * dg_15[k]
                 - f_2 * dg_18[k]
                 + f_3 * dg_20[k]
                 + f_4 * dg_46[k]
                 + f_4 * dg_51[k]
                 - f_5 * dg_53[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_12, df_17, df_32, df_37, dg_17, dg_22, dg_49, \
                         dg_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -f_10 * ab_x[k] * df_12[k]
                 + f_10 * ab_x[k] * df_17[k]
                 + f_11 * ab_y[k] * df_32[k]
                 - f_11 * ab_y[k] * df_37[k]
                 + f_10 * dg_17[k]
                 - f_10 * dg_22[k]
                 - f_11 * dg_49[k]
                 + f_11 * dg_56[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_10, df_13, df_30, df_33, df_41, df_46, dg_15, dg_18, \
                         dg_46, dg_51, dg_61, dg_66 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -1.875 * ab_x[k] * df_10[k]
                 + 5.625 * ab_x[k] * df_13[k]
                 + 0.625 * ab_y[k] * df_30[k]
                 - 1.875 * ab_y[k] * df_33[k]
                 + 1.875 * dg_15[k]
                 - 5.625 * dg_18[k]
                 - 0.625 * dg_46[k]
                 + 1.875 * dg_51[k];

        g_7[k] = -f_0 * ab_x[k] * df_41[k]
                 + f_1 * ab_x[k] * df_46[k]
                 + f_0 * dg_61[k]
                 - f_1 * dg_66[k];
    }

#pragma omp simd aligned(ab_x, df_41, df_44, df_46, df_48, dg_61, dg_64, dg_66, \
                         dg_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -15.0 * ab_x[k] * df_44[k]
                 + 15.0 * dg_64[k];

        g_9[k] = f_7 * ab_x[k] * df_41[k]
                 + f_7 * ab_x[k] * df_46[k]
                 - f_12 * ab_x[k] * df_48[k]
                 - f_7 * dg_61[k]
                 - f_7 * dg_66[k]
                 + f_12 * dg_68[k];
    }

#pragma omp simd aligned(ab_x, df_40, df_42, df_43, df_45, df_47, df_49, dg_60, dg_62, dg_63, \
                         dg_65, dg_67, dg_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = f_3 * ab_x[k] * df_42[k]
                  + f_3 * ab_x[k] * df_47[k]
                  - f_13 * ab_x[k] * df_49[k]
                  - f_3 * dg_62[k]
                  - f_3 * dg_67[k]
                  + f_13 * dg_69[k];

        g_11[k] = f_7 * ab_x[k] * df_40[k]
                  + f_7 * ab_x[k] * df_43[k]
                  - f_12 * ab_x[k] * df_45[k]
                  - f_7 * dg_60[k]
                  - f_7 * dg_63[k]
                  + f_12 * dg_65[k];

        g_12[k] = -7.5 * ab_x[k] * df_42[k]
                  + 7.5 * ab_x[k] * df_47[k]
                  + 7.5 * dg_62[k]
                  - 7.5 * dg_67[k];
    }

#pragma omp simd aligned(ab_x, df_40, df_43, dg_60, dg_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_1 * ab_x[k] * df_40[k]
                  + f_0 * ab_x[k] * df_43[k]
                  + f_1 * dg_60[k]
                  - f_0 * dg_63[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_11, df_16, df_31, df_36, df_51, df_56, dg_16, dg_21, \
                         dg_48, dg_55, dg_78, dg_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_2 * ab_x[k] * df_11[k]
                  - f_4 * ab_x[k] * df_16[k]
                  + f_2 * ab_y[k] * df_31[k]
                  - f_4 * ab_y[k] * df_36[k]
                  - f_3 * ab_y[k] * df_51[k]
                  + f_5 * ab_y[k] * df_56[k]
                  - f_2 * dg_16[k]
                  + f_4 * dg_21[k]
                  - f_2 * dg_48[k]
                  + f_4 * dg_55[k]
                  + f_3 * dg_78[k]
                  - f_5 * dg_85[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_14, df_34, df_54, dg_19, dg_52, \
                         dg_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_7 * ab_x[k] * df_14[k]
                  + f_7 * ab_y[k] * df_34[k]
                  - f_12 * ab_y[k] * df_54[k]
                  - f_7 * dg_19[k]
                  - f_7 * dg_52[k]
                  + f_12 * dg_82[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_11, df_16, df_18, df_31, df_36, df_38, df_51, df_56, \
                         df_58, dg_16, dg_21, dg_23, dg_48, dg_55, dg_57, dg_78, dg_85, \
                         dg_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = -0.375 * ab_x[k] * df_11[k]
                  - 0.375 * ab_x[k] * df_16[k]
                  + 1.5 * ab_x[k] * df_18[k]
                  - 0.375 * ab_y[k] * df_31[k]
                  - 0.375 * ab_y[k] * df_36[k]
                  + 1.5 * ab_y[k] * df_38[k]
                  + 1.5 * ab_y[k] * df_51[k]
                  + 1.5 * ab_y[k] * df_56[k]
                  - 6.0 * ab_y[k] * df_58[k]
                  + 0.375 * dg_16[k]
                  + 0.375 * dg_21[k]
                  - 1.5 * dg_23[k]
                  + 0.375 * dg_48[k]
                  + 0.375 * dg_55[k]
                  - 1.5 * dg_57[k]
                  - 1.5 * dg_78[k]
                  - 1.5 * dg_85[k]
                  + 6.0 * dg_87[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_12, df_17, df_19, df_32, df_37, df_39, df_52, df_57, \
                         df_59, dg_17, dg_22, dg_24, dg_49, dg_56, dg_58, dg_79, dg_86, \
                         dg_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_14 * ab_x[k] * df_12[k]
                  - f_14 * ab_x[k] * df_17[k]
                  + f_15 * ab_x[k] * df_19[k]
                  - f_14 * ab_y[k] * df_32[k]
                  - f_14 * ab_y[k] * df_37[k]
                  + f_15 * ab_y[k] * df_39[k]
                  + f_16 * ab_y[k] * df_52[k]
                  + f_16 * ab_y[k] * df_57[k]
                  - f_17 * ab_y[k] * df_59[k]
                  + f_14 * dg_17[k]
                  + f_14 * dg_22[k]
                  - f_15 * dg_24[k]
                  + f_14 * dg_49[k]
                  + f_14 * dg_56[k]
                  - f_15 * dg_58[k]
                  - f_16 * dg_79[k]
                  - f_16 * dg_86[k]
                  + f_17 * dg_88[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_10, df_13, df_15, df_30, df_33, df_35, df_50, df_53, \
                         df_55, dg_15, dg_18, dg_20, dg_46, dg_51, dg_53, dg_76, dg_81, \
                         dg_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -0.375 * ab_x[k] * df_10[k]
                  - 0.375 * ab_x[k] * df_13[k]
                  + 1.5 * ab_x[k] * df_15[k]
                  - 0.375 * ab_y[k] * df_30[k]
                  - 0.375 * ab_y[k] * df_33[k]
                  + 1.5 * ab_y[k] * df_35[k]
                  + 1.5 * ab_y[k] * df_50[k]
                  + 1.5 * ab_y[k] * df_53[k]
                  - 6.0 * ab_y[k] * df_55[k]
                  + 0.375 * dg_15[k]
                  + 0.375 * dg_18[k]
                  - 1.5 * dg_20[k]
                  + 0.375 * dg_46[k]
                  + 0.375 * dg_51[k]
                  - 1.5 * dg_53[k]
                  - 1.5 * dg_76[k]
                  - 1.5 * dg_81[k]
                  + 6.0 * dg_83[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_12, df_17, df_32, df_37, df_52, df_57, dg_17, dg_22, \
                         dg_49, dg_56, dg_79, dg_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_8 * ab_x[k] * df_12[k]
                  - f_8 * ab_x[k] * df_17[k]
                  + f_8 * ab_y[k] * df_32[k]
                  - f_8 * ab_y[k] * df_37[k]
                  - f_18 * ab_y[k] * df_52[k]
                  + f_18 * ab_y[k] * df_57[k]
                  - f_8 * dg_17[k]
                  + f_8 * dg_22[k]
                  - f_8 * dg_49[k]
                  + f_8 * dg_56[k]
                  + f_18 * dg_79[k]
                  - f_18 * dg_86[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_10, df_13, df_30, df_33, df_50, df_53, dg_15, dg_18, \
                         dg_46, dg_51, dg_76, dg_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_4 * ab_x[k] * df_10[k]
                  - f_2 * ab_x[k] * df_13[k]
                  + f_4 * ab_y[k] * df_30[k]
                  - f_2 * ab_y[k] * df_33[k]
                  - f_5 * ab_y[k] * df_50[k]
                  + f_3 * ab_y[k] * df_53[k]
                  - f_4 * dg_15[k]
                  + f_2 * dg_18[k]
                  - f_4 * dg_46[k]
                  + f_2 * dg_51[k]
                  + f_5 * dg_76[k]
                  - f_3 * dg_81[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, df_21, df_26, df_41, df_46, df_51, df_56, dg_31, \
                         dg_36, dg_63, dg_70, dg_79, dg_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_6 * ab_x[k] * df_21[k]
                  - f_8 * ab_x[k] * df_26[k]
                  + f_6 * ab_y[k] * df_41[k]
                  - f_8 * ab_y[k] * df_46[k]
                  - f_7 * ab_z[k] * df_51[k]
                  + f_9 * ab_z[k] * df_56[k]
                  - f_6 * dg_31[k]
                  + f_8 * dg_36[k]
                  - f_6 * dg_63[k]
                  + f_8 * dg_70[k]
                  + f_7 * dg_79[k]
                  - f_9 * dg_86[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, df_24, df_44, df_54, dg_34, dg_67, \
                         dg_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_3 * ab_x[k] * df_24[k]
                  + f_3 * ab_y[k] * df_44[k]
                  - f_13 * ab_z[k] * df_54[k]
                  - f_3 * dg_34[k]
                  - f_3 * dg_67[k]
                  + f_13 * dg_83[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, df_21, df_26, df_28, df_41, df_46, df_48, df_51, \
                         df_56, df_58, dg_31, dg_36, dg_38, dg_63, dg_70, dg_72, dg_79, dg_86, \
                         dg_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_14 * ab_x[k] * df_21[k]
                  - f_14 * ab_x[k] * df_26[k]
                  + f_16 * ab_x[k] * df_28[k]
                  - f_14 * ab_y[k] * df_41[k]
                  - f_14 * ab_y[k] * df_46[k]
                  + f_16 * ab_y[k] * df_48[k]
                  + f_15 * ab_z[k] * df_51[k]
                  + f_15 * ab_z[k] * df_56[k]
                  - f_17 * ab_z[k] * df_58[k]
                  + f_14 * dg_31[k]
                  + f_14 * dg_36[k]
                  - f_16 * dg_38[k]
                  + f_14 * dg_63[k]
                  + f_14 * dg_70[k]
                  - f_16 * dg_72[k]
                  - f_15 * dg_79[k]
                  - f_15 * dg_86[k]
                  + f_17 * dg_88[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, df_22, df_27, df_29, df_42, df_47, df_49, df_52, \
                         df_57, df_59, dg_32, dg_37, dg_39, dg_64, dg_71, dg_73, dg_80, dg_87, \
                         dg_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -2.25 * ab_x[k] * df_22[k]
                  - 2.25 * ab_x[k] * df_27[k]
                  + 1.5 * ab_x[k] * df_29[k]
                  - 2.25 * ab_y[k] * df_42[k]
                  - 2.25 * ab_y[k] * df_47[k]
                  + 1.5 * ab_y[k] * df_49[k]
                  + 1.5 * ab_z[k] * df_52[k]
                  + 1.5 * ab_z[k] * df_57[k]
                  - ab_z[k] * df_59[k]
                  + 2.25 * dg_32[k]
                  + 2.25 * dg_37[k]
                  - 1.5 * dg_39[k]
                  + 2.25 * dg_64[k]
                  + 2.25 * dg_71[k]
                  - 1.5 * dg_73[k]
                  - 1.5 * dg_80[k]
                  - 1.5 * dg_87[k]
                  + dg_89[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, df_20, df_23, df_25, df_40, df_43, df_45, df_50, \
                         df_53, df_55, dg_30, dg_33, dg_35, dg_61, dg_66, dg_68, dg_77, dg_82, \
                         dg_84 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_14 * ab_x[k] * df_20[k]
                  - f_14 * ab_x[k] * df_23[k]
                  + f_16 * ab_x[k] * df_25[k]
                  - f_14 * ab_y[k] * df_40[k]
                  - f_14 * ab_y[k] * df_43[k]
                  + f_16 * ab_y[k] * df_45[k]
                  + f_15 * ab_z[k] * df_50[k]
                  + f_15 * ab_z[k] * df_53[k]
                  - f_17 * ab_z[k] * df_55[k]
                  + f_14 * dg_30[k]
                  + f_14 * dg_33[k]
                  - f_16 * dg_35[k]
                  + f_14 * dg_61[k]
                  + f_14 * dg_66[k]
                  - f_16 * dg_68[k]
                  - f_15 * dg_77[k]
                  - f_15 * dg_82[k]
                  + f_17 * dg_84[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, df_22, df_27, df_42, df_47, df_52, df_57, dg_32, \
                         dg_37, dg_64, dg_71, dg_80, dg_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = f_19 * ab_x[k] * df_22[k]
                  - f_19 * ab_x[k] * df_27[k]
                  + f_19 * ab_y[k] * df_42[k]
                  - f_19 * ab_y[k] * df_47[k]
                  - f_5 * ab_z[k] * df_52[k]
                  + f_5 * ab_z[k] * df_57[k]
                  - f_19 * dg_32[k]
                  + f_19 * dg_37[k]
                  - f_19 * dg_64[k]
                  + f_19 * dg_71[k]
                  + f_5 * dg_80[k]
                  - f_5 * dg_87[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, df_20, df_23, df_40, df_43, df_50, df_53, dg_30, \
                         dg_33, dg_61, dg_66, dg_77, dg_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = f_8 * ab_x[k] * df_20[k]
                  - f_6 * ab_x[k] * df_23[k]
                  + f_8 * ab_y[k] * df_40[k]
                  - f_6 * ab_y[k] * df_43[k]
                  - f_9 * ab_z[k] * df_50[k]
                  + f_7 * ab_z[k] * df_53[k]
                  - f_8 * dg_30[k]
                  + f_6 * dg_33[k]
                  - f_8 * dg_61[k]
                  + f_6 * dg_66[k]
                  + f_9 * dg_77[k]
                  - f_7 * dg_82[k];
    }

#pragma omp simd aligned(ab_x, df_1, df_6, df_31, df_36, df_51, df_56, dg_1, dg_6, dg_46, \
                         dg_51, dg_76, dg_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_2 * ab_x[k] * df_1[k]
                  - f_4 * ab_x[k] * df_6[k]
                  + f_2 * ab_x[k] * df_31[k]
                  - f_4 * ab_x[k] * df_36[k]
                  - f_3 * ab_x[k] * df_51[k]
                  + f_5 * ab_x[k] * df_56[k]
                  - f_2 * dg_1[k]
                  + f_4 * dg_6[k]
                  - f_2 * dg_46[k]
                  + f_4 * dg_51[k]
                  + f_3 * dg_76[k]
                  - f_5 * dg_81[k];
    }

#pragma omp simd aligned(ab_x, df_4, df_34, df_54, dg_4, dg_49, dg_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_7 * ab_x[k] * df_4[k]
                  + f_7 * ab_x[k] * df_34[k]
                  - f_12 * ab_x[k] * df_54[k]
                  - f_7 * dg_4[k]
                  - f_7 * dg_49[k]
                  + f_12 * dg_79[k];
    }

#pragma omp simd aligned(ab_x, df_1, df_6, df_8, df_31, df_36, df_38, df_51, df_56, df_58, \
                         dg_1, dg_6, dg_8, dg_46, dg_51, dg_53, dg_76, dg_81, \
                         dg_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -0.375 * ab_x[k] * df_1[k]
                  - 0.375 * ab_x[k] * df_6[k]
                  + 1.5 * ab_x[k] * df_8[k]
                  - 0.375 * ab_x[k] * df_31[k]
                  - 0.375 * ab_x[k] * df_36[k]
                  + 1.5 * ab_x[k] * df_38[k]
                  + 1.5 * ab_x[k] * df_51[k]
                  + 1.5 * ab_x[k] * df_56[k]
                  - 6.0 * ab_x[k] * df_58[k]
                  + 0.375 * dg_1[k]
                  + 0.375 * dg_6[k]
                  - 1.5 * dg_8[k]
                  + 0.375 * dg_46[k]
                  + 0.375 * dg_51[k]
                  - 1.5 * dg_53[k]
                  - 1.5 * dg_76[k]
                  - 1.5 * dg_81[k]
                  + 6.0 * dg_83[k];
    }

#pragma omp simd aligned(ab_x, df_2, df_7, df_9, df_32, df_37, df_39, df_52, df_57, df_59, \
                         dg_2, dg_7, dg_9, dg_47, dg_52, dg_54, dg_77, dg_82, \
                         dg_84 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_14 * ab_x[k] * df_2[k]
                  - f_14 * ab_x[k] * df_7[k]
                  + f_15 * ab_x[k] * df_9[k]
                  - f_14 * ab_x[k] * df_32[k]
                  - f_14 * ab_x[k] * df_37[k]
                  + f_15 * ab_x[k] * df_39[k]
                  + f_16 * ab_x[k] * df_52[k]
                  + f_16 * ab_x[k] * df_57[k]
                  - f_17 * ab_x[k] * df_59[k]
                  + f_14 * dg_2[k]
                  + f_14 * dg_7[k]
                  - f_15 * dg_9[k]
                  + f_14 * dg_47[k]
                  + f_14 * dg_52[k]
                  - f_15 * dg_54[k]
                  - f_16 * dg_77[k]
                  - f_16 * dg_82[k]
                  + f_17 * dg_84[k];
    }

#pragma omp simd aligned(ab_x, df_0, df_3, df_5, df_30, df_33, df_35, df_50, df_53, df_55, \
                         dg_0, dg_3, dg_5, dg_45, dg_48, dg_50, dg_75, dg_78, \
                         dg_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -0.375 * ab_x[k] * df_0[k]
                  - 0.375 * ab_x[k] * df_3[k]
                  + 1.5 * ab_x[k] * df_5[k]
                  - 0.375 * ab_x[k] * df_30[k]
                  - 0.375 * ab_x[k] * df_33[k]
                  + 1.5 * ab_x[k] * df_35[k]
                  + 1.5 * ab_x[k] * df_50[k]
                  + 1.5 * ab_x[k] * df_53[k]
                  - 6.0 * ab_x[k] * df_55[k]
                  + 0.375 * dg_0[k]
                  + 0.375 * dg_3[k]
                  - 1.5 * dg_5[k]
                  + 0.375 * dg_45[k]
                  + 0.375 * dg_48[k]
                  - 1.5 * dg_50[k]
                  - 1.5 * dg_75[k]
                  - 1.5 * dg_78[k]
                  + 6.0 * dg_80[k];
    }

#pragma omp simd aligned(ab_x, df_2, df_7, df_32, df_37, df_52, df_57, dg_2, dg_7, dg_47, \
                         dg_52, dg_77, dg_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_8 * ab_x[k] * df_2[k]
                  - f_8 * ab_x[k] * df_7[k]
                  + f_8 * ab_x[k] * df_32[k]
                  - f_8 * ab_x[k] * df_37[k]
                  - f_18 * ab_x[k] * df_52[k]
                  + f_18 * ab_x[k] * df_57[k]
                  - f_8 * dg_2[k]
                  + f_8 * dg_7[k]
                  - f_8 * dg_47[k]
                  + f_8 * dg_52[k]
                  + f_18 * dg_77[k]
                  - f_18 * dg_82[k];
    }

#pragma omp simd aligned(ab_x, df_0, df_3, df_30, df_33, df_50, df_53, dg_0, dg_3, dg_45, \
                         dg_48, dg_75, dg_78 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_4 * ab_x[k] * df_0[k]
                  - f_2 * ab_x[k] * df_3[k]
                  + f_4 * ab_x[k] * df_30[k]
                  - f_2 * ab_x[k] * df_33[k]
                  - f_5 * ab_x[k] * df_50[k]
                  + f_3 * ab_x[k] * df_53[k]
                  - f_4 * dg_0[k]
                  + f_2 * dg_3[k]
                  - f_4 * dg_45[k]
                  + f_2 * dg_48[k]
                  + f_5 * dg_75[k]
                  - f_3 * dg_78[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_21, df_24, df_26, df_41, df_44, df_46, dg_31, dg_34, \
                         dg_36, dg_63, dg_67, dg_70 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_10 * ab_x[k] * df_21[k]
                  + f_11 * ab_x[k] * df_26[k]
                  + f_10 * ab_y[k] * df_41[k]
                  - f_11 * ab_y[k] * df_46[k]
                  + f_10 * dg_31[k]
                  - f_11 * dg_36[k]
                  - f_10 * dg_63[k]
                  + f_11 * dg_70[k];

        g_36[k] = -7.5 * ab_x[k] * df_24[k]
                  + 7.5 * ab_y[k] * df_44[k]
                  + 7.5 * dg_34[k]
                  - 7.5 * dg_67[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_21, df_26, df_28, df_41, df_46, df_48, dg_31, dg_36, \
                         dg_38, dg_63, dg_70, dg_72 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_8 * ab_x[k] * df_21[k]
                  + f_8 * ab_x[k] * df_26[k]
                  - f_18 * ab_x[k] * df_28[k]
                  - f_8 * ab_y[k] * df_41[k]
                  - f_8 * ab_y[k] * df_46[k]
                  + f_18 * ab_y[k] * df_48[k]
                  - f_8 * dg_31[k]
                  - f_8 * dg_36[k]
                  + f_18 * dg_38[k]
                  + f_8 * dg_63[k]
                  + f_8 * dg_70[k]
                  - f_18 * dg_72[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_22, df_27, df_29, df_42, df_47, df_49, dg_32, dg_37, \
                         dg_39, dg_64, dg_71, dg_73 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_19 * ab_x[k] * df_22[k]
                  + f_19 * ab_x[k] * df_27[k]
                  - f_5 * ab_x[k] * df_29[k]
                  - f_19 * ab_y[k] * df_42[k]
                  - f_19 * ab_y[k] * df_47[k]
                  + f_5 * ab_y[k] * df_49[k]
                  - f_19 * dg_32[k]
                  - f_19 * dg_37[k]
                  + f_5 * dg_39[k]
                  + f_19 * dg_64[k]
                  + f_19 * dg_71[k]
                  - f_5 * dg_73[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_20, df_23, df_25, df_40, df_43, df_45, dg_30, dg_33, \
                         dg_35, dg_61, dg_66, dg_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = f_8 * ab_x[k] * df_20[k]
                  + f_8 * ab_x[k] * df_23[k]
                  - f_18 * ab_x[k] * df_25[k]
                  - f_8 * ab_y[k] * df_40[k]
                  - f_8 * ab_y[k] * df_43[k]
                  + f_18 * ab_y[k] * df_45[k]
                  - f_8 * dg_30[k]
                  - f_8 * dg_33[k]
                  + f_18 * dg_35[k]
                  + f_8 * dg_61[k]
                  + f_8 * dg_66[k]
                  - f_18 * dg_68[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_22, df_27, df_42, df_47, dg_32, dg_37, dg_64, \
                         dg_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -3.75 * ab_x[k] * df_22[k]
                  + 3.75 * ab_x[k] * df_27[k]
                  + 3.75 * ab_y[k] * df_42[k]
                  - 3.75 * ab_y[k] * df_47[k]
                  + 3.75 * dg_32[k]
                  - 3.75 * dg_37[k]
                  - 3.75 * dg_64[k]
                  + 3.75 * dg_71[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_20, df_23, df_40, df_43, dg_30, dg_33, dg_61, \
                         dg_66 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_11 * ab_x[k] * df_20[k]
                  + f_10 * ab_x[k] * df_23[k]
                  + f_11 * ab_y[k] * df_40[k]
                  - f_10 * ab_y[k] * df_43[k]
                  + f_11 * dg_30[k]
                  - f_10 * dg_33[k]
                  - f_11 * dg_61[k]
                  + f_10 * dg_66[k];
    }

#pragma omp simd aligned(ab_x, df_1, df_4, df_6, df_31, df_34, df_36, dg_1, dg_4, dg_6, dg_46, \
                         dg_49, dg_51 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -1.875 * ab_x[k] * df_1[k]
                  + 0.625 * ab_x[k] * df_6[k]
                  + 5.625 * ab_x[k] * df_31[k]
                  - 1.875 * ab_x[k] * df_36[k]
                  + 1.875 * dg_1[k]
                  - 0.625 * dg_6[k]
                  - 5.625 * dg_46[k]
                  + 1.875 * dg_51[k];

        g_43[k] = -f_1 * ab_x[k] * df_4[k]
                  + f_0 * ab_x[k] * df_34[k]
                  + f_1 * dg_4[k]
                  - f_0 * dg_49[k];
    }

#pragma omp simd aligned(ab_x, df_1, df_6, df_8, df_31, df_36, df_38, dg_1, dg_6, dg_8, dg_46, \
                         dg_51, dg_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_4 * ab_x[k] * df_1[k]
                  + f_4 * ab_x[k] * df_6[k]
                  - f_5 * ab_x[k] * df_8[k]
                  - f_2 * ab_x[k] * df_31[k]
                  - f_2 * ab_x[k] * df_36[k]
                  + f_3 * ab_x[k] * df_38[k]
                  - f_4 * dg_1[k]
                  - f_4 * dg_6[k]
                  + f_5 * dg_8[k]
                  + f_2 * dg_46[k]
                  + f_2 * dg_51[k]
                  - f_3 * dg_53[k];
    }

#pragma omp simd aligned(ab_x, df_2, df_7, df_9, df_32, df_37, df_39, dg_2, dg_7, dg_9, dg_47, \
                         dg_52, dg_54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_8 * ab_x[k] * df_2[k]
                  + f_8 * ab_x[k] * df_7[k]
                  - f_9 * ab_x[k] * df_9[k]
                  - f_6 * ab_x[k] * df_32[k]
                  - f_6 * ab_x[k] * df_37[k]
                  + f_7 * ab_x[k] * df_39[k]
                  - f_8 * dg_2[k]
                  - f_8 * dg_7[k]
                  + f_9 * dg_9[k]
                  + f_6 * dg_47[k]
                  + f_6 * dg_52[k]
                  - f_7 * dg_54[k];
    }

#pragma omp simd aligned(ab_x, df_0, df_3, df_5, df_30, df_33, df_35, dg_0, dg_3, dg_5, dg_45, \
                         dg_48, dg_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_4 * ab_x[k] * df_0[k]
                  + f_4 * ab_x[k] * df_3[k]
                  - f_5 * ab_x[k] * df_5[k]
                  - f_2 * ab_x[k] * df_30[k]
                  - f_2 * ab_x[k] * df_33[k]
                  + f_3 * ab_x[k] * df_35[k]
                  - f_4 * dg_0[k]
                  - f_4 * dg_3[k]
                  + f_5 * dg_5[k]
                  + f_2 * dg_45[k]
                  + f_2 * dg_48[k]
                  - f_3 * dg_50[k];
    }

#pragma omp simd aligned(ab_x, df_2, df_7, df_32, df_37, dg_2, dg_7, dg_47, \
                         dg_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_11 * ab_x[k] * df_2[k]
                  + f_11 * ab_x[k] * df_7[k]
                  + f_10 * ab_x[k] * df_32[k]
                  - f_10 * ab_x[k] * df_37[k]
                  + f_11 * dg_2[k]
                  - f_11 * dg_7[k]
                  - f_10 * dg_47[k]
                  + f_10 * dg_52[k];
    }

#pragma omp simd aligned(ab_x, df_0, df_3, df_30, df_33, dg_0, dg_3, dg_45, \
                         dg_48 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -0.625 * ab_x[k] * df_0[k]
                  + 1.875 * ab_x[k] * df_3[k]
                  + 1.875 * ab_x[k] * df_30[k]
                  - 5.625 * ab_x[k] * df_33[k]
                  + 0.625 * dg_0[k]
                  - 1.875 * dg_3[k]
                  - 1.875 * dg_45[k]
                  + 5.625 * dg_48[k];
    }
}

auto
compute_hrr_ff_sph_tri(double *values, const size_t nvalues, CSimdMatrix &buffer,
                       const CSimdMatrix &coordinates, const size_t df, const size_t dg,
                       const size_t nmax) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.75 * std::sqrt(6.0);
    const auto f_1 = 1.25 * std::sqrt(6.0);
    const auto f_2 = 0.375 * std::sqrt(15.0);
    const auto f_3 = 1.5 * std::sqrt(15.0);
    const auto f_4 = 0.125 * std::sqrt(15.0);
    const auto f_5 = 0.5 * std::sqrt(15.0);
    const auto f_6 = 1.125 * std::sqrt(10.0);
    const auto f_7 = 0.75 * std::sqrt(10.0);
    const auto f_8 = 0.375 * std::sqrt(10.0);
    const auto f_9 = 0.25 * std::sqrt(10.0);
    const auto f_10 = 1.875 * std::sqrt(6.0);
    const auto f_11 = 0.625 * std::sqrt(6.0);
    const auto f_12 = 3.0 * std::sqrt(10.0);
    const auto f_13 = std::sqrt(15.0);
    const auto f_14 = 0.375 * std::sqrt(6.0);
    const auto f_15 = 0.25 * std::sqrt(6.0);
    const auto f_16 = 1.5 * std::sqrt(6.0);
    const auto f_17 = std::sqrt(6.0);
    const auto f_18 = 1.5 * std::sqrt(10.0);
    const auto f_19 = 0.75 * std::sqrt(15.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;
    auto *g_7 = values + 7 * nvalues;
    auto *g_8 = values + 8 * nvalues;
    auto *g_9 = values + 9 * nvalues;
    auto *g_10 = values + 10 * nvalues;
    auto *g_11 = values + 11 * nvalues;
    auto *g_12 = values + 12 * nvalues;
    auto *g_13 = values + 13 * nvalues;
    auto *g_14 = values + 14 * nvalues;
    auto *g_15 = values + 15 * nvalues;
    auto *g_16 = values + 16 * nvalues;
    auto *g_17 = values + 17 * nvalues;
    auto *g_18 = values + 18 * nvalues;
    auto *g_19 = values + 19 * nvalues;
    auto *g_20 = values + 20 * nvalues;
    auto *g_21 = values + 21 * nvalues;
    auto *g_22 = values + 22 * nvalues;
    auto *g_23 = values + 23 * nvalues;
    auto *g_24 = values + 24 * nvalues;
    auto *g_25 = values + 25 * nvalues;
    auto *g_26 = values + 26 * nvalues;
    auto *g_27 = values + 27 * nvalues;
    auto *g_28 = values + 28 * nvalues;
    auto *g_29 = values + 29 * nvalues;
    auto *g_30 = values + 30 * nvalues;
    auto *g_31 = values + 31 * nvalues;
    auto *g_32 = values + 32 * nvalues;
    auto *g_33 = values + 33 * nvalues;
    auto *g_34 = values + 34 * nvalues;
    auto *g_35 = values + 35 * nvalues;
    auto *g_36 = values + 36 * nvalues;
    auto *g_37 = values + 37 * nvalues;
    auto *g_38 = values + 38 * nvalues;
    auto *g_39 = values + 39 * nvalues;
    auto *g_40 = values + 40 * nvalues;
    auto *g_41 = values + 41 * nvalues;
    auto *g_42 = values + 42 * nvalues;
    auto *g_43 = values + 43 * nvalues;
    auto *g_44 = values + 44 * nvalues;
    auto *g_45 = values + 45 * nvalues;
    auto *g_46 = values + 46 * nvalues;
    auto *g_47 = values + 47 * nvalues;
    auto *g_48 = values + 48 * nvalues;

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_7 = buffer.data(df + 7);
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
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_27 = buffer.data(df + 27);
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
    const auto *df_55 = buffer.data(df + 55);
    const auto *df_56 = buffer.data(df + 56);
    const auto *df_57 = buffer.data(df + 57);
    const auto *df_58 = buffer.data(df + 58);
    const auto *df_59 = buffer.data(df + 59);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_7 = buffer.data(dg + 7);
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
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_37 = buffer.data(dg + 37);
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
    const auto *dg_71 = buffer.data(dg + 71);
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

#pragma omp simd aligned(ab_x, ab_y, df_11, df_14, df_16, df_31, df_34, df_36, dg_16, dg_19, \
                         dg_21, dg_48, dg_52, dg_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -5.625 * ab_x[k] * df_11[k]
                 + 1.875 * ab_x[k] * df_16[k]
                 + 1.875 * ab_y[k] * df_31[k]
                 - 0.625 * ab_y[k] * df_36[k]
                 + 5.625 * dg_16[k]
                 - 1.875 * dg_21[k]
                 - 1.875 * dg_48[k]
                 + 0.625 * dg_55[k];

        g_1[k] = -f_0 * ab_x[k] * df_14[k]
                 + f_1 * ab_y[k] * df_34[k]
                 + f_0 * dg_19[k]
                 - f_1 * dg_52[k];
        g_7[k] = g_1[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_11, df_16, df_18, df_31, df_36, df_38, dg_16, dg_21, \
                         dg_23, dg_48, dg_55, dg_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = f_2 * ab_x[k] * df_11[k]
                 + f_2 * ab_x[k] * df_16[k]
                 - f_3 * ab_x[k] * df_18[k]
                 - f_4 * ab_y[k] * df_31[k]
                 - f_4 * ab_y[k] * df_36[k]
                 + f_5 * ab_y[k] * df_38[k]
                 - f_2 * dg_16[k]
                 - f_2 * dg_21[k]
                 + f_3 * dg_23[k]
                 + f_4 * dg_48[k]
                 + f_4 * dg_55[k]
                 - f_5 * dg_57[k];
        g_14[k] = g_2[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_12, df_17, df_19, df_32, df_37, df_39, dg_17, dg_22, \
                         dg_24, dg_49, dg_56, dg_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_6 * ab_x[k] * df_12[k]
                 + f_6 * ab_x[k] * df_17[k]
                 - f_7 * ab_x[k] * df_19[k]
                 - f_8 * ab_y[k] * df_32[k]
                 - f_8 * ab_y[k] * df_37[k]
                 + f_9 * ab_y[k] * df_39[k]
                 - f_6 * dg_17[k]
                 - f_6 * dg_22[k]
                 + f_7 * dg_24[k]
                 + f_8 * dg_49[k]
                 + f_8 * dg_56[k]
                 - f_9 * dg_58[k];
        g_21[k] = g_3[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_10, df_13, df_15, df_30, df_33, df_35, dg_15, dg_18, \
                         dg_20, dg_46, dg_51, dg_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_2 * ab_x[k] * df_10[k]
                 + f_2 * ab_x[k] * df_13[k]
                 - f_3 * ab_x[k] * df_15[k]
                 - f_4 * ab_y[k] * df_30[k]
                 - f_4 * ab_y[k] * df_33[k]
                 + f_5 * ab_y[k] * df_35[k]
                 - f_2 * dg_15[k]
                 - f_2 * dg_18[k]
                 + f_3 * dg_20[k]
                 + f_4 * dg_46[k]
                 + f_4 * dg_51[k]
                 - f_5 * dg_53[k];
        g_28[k] = g_4[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_12, df_17, df_32, df_37, dg_17, dg_22, dg_49, \
                         dg_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -f_10 * ab_x[k] * df_12[k]
                 + f_10 * ab_x[k] * df_17[k]
                 + f_11 * ab_y[k] * df_32[k]
                 - f_11 * ab_y[k] * df_37[k]
                 + f_10 * dg_17[k]
                 - f_10 * dg_22[k]
                 - f_11 * dg_49[k]
                 + f_11 * dg_56[k];
        g_35[k] = g_5[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_10, df_13, df_30, df_33, df_44, dg_15, dg_18, dg_46, \
                         dg_51, dg_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -1.875 * ab_x[k] * df_10[k]
                 + 5.625 * ab_x[k] * df_13[k]
                 + 0.625 * ab_y[k] * df_30[k]
                 - 1.875 * ab_y[k] * df_33[k]
                 + 1.875 * dg_15[k]
                 - 5.625 * dg_18[k]
                 - 0.625 * dg_46[k]
                 + 1.875 * dg_51[k];
        g_42[k] = g_6[k];

        g_8[k] = -15.0 * ab_x[k] * df_44[k]
                 + 15.0 * dg_64[k];
    }

#pragma omp simd aligned(ab_x, df_41, df_42, df_46, df_47, df_48, df_49, dg_61, dg_62, dg_66, \
                         dg_67, dg_68, dg_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_7 * ab_x[k] * df_41[k]
                 + f_7 * ab_x[k] * df_46[k]
                 - f_12 * ab_x[k] * df_48[k]
                 - f_7 * dg_61[k]
                 - f_7 * dg_66[k]
                 + f_12 * dg_68[k];
        g_15[k] = g_9[k];

        g_10[k] = f_3 * ab_x[k] * df_42[k]
                  + f_3 * ab_x[k] * df_47[k]
                  - f_13 * ab_x[k] * df_49[k]
                  - f_3 * dg_62[k]
                  - f_3 * dg_67[k]
                  + f_13 * dg_69[k];
        g_22[k] = g_10[k];
    }

#pragma omp simd aligned(ab_x, df_40, df_42, df_43, df_45, df_47, dg_60, dg_62, dg_63, dg_65, \
                         dg_67 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_7 * ab_x[k] * df_40[k]
                  + f_7 * ab_x[k] * df_43[k]
                  - f_12 * ab_x[k] * df_45[k]
                  - f_7 * dg_60[k]
                  - f_7 * dg_63[k]
                  + f_12 * dg_65[k];
        g_29[k] = g_11[k];

        g_12[k] = -7.5 * ab_x[k] * df_42[k]
                  + 7.5 * ab_x[k] * df_47[k]
                  + 7.5 * dg_62[k]
                  - 7.5 * dg_67[k];
        g_36[k] = g_12[k];

        g_13[k] = -f_1 * ab_x[k] * df_40[k]
                  + f_0 * ab_x[k] * df_43[k]
                  + f_1 * dg_60[k]
                  - f_0 * dg_63[k];
        g_43[k] = g_13[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_11, df_16, df_18, df_31, df_36, df_38, df_51, df_56, \
                         df_58, dg_16, dg_21, dg_23, dg_48, dg_55, dg_57, dg_78, dg_85, \
                         dg_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = -0.375 * ab_x[k] * df_11[k]
                  - 0.375 * ab_x[k] * df_16[k]
                  + 1.5 * ab_x[k] * df_18[k]
                  - 0.375 * ab_y[k] * df_31[k]
                  - 0.375 * ab_y[k] * df_36[k]
                  + 1.5 * ab_y[k] * df_38[k]
                  + 1.5 * ab_y[k] * df_51[k]
                  + 1.5 * ab_y[k] * df_56[k]
                  - 6.0 * ab_y[k] * df_58[k]
                  + 0.375 * dg_16[k]
                  + 0.375 * dg_21[k]
                  - 1.5 * dg_23[k]
                  + 0.375 * dg_48[k]
                  + 0.375 * dg_55[k]
                  - 1.5 * dg_57[k]
                  - 1.5 * dg_78[k]
                  - 1.5 * dg_85[k]
                  + 6.0 * dg_87[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_12, df_17, df_19, df_32, df_37, df_39, df_52, df_57, \
                         df_59, dg_17, dg_22, dg_24, dg_49, dg_56, dg_58, dg_79, dg_86, \
                         dg_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_14 * ab_x[k] * df_12[k]
                  - f_14 * ab_x[k] * df_17[k]
                  + f_15 * ab_x[k] * df_19[k]
                  - f_14 * ab_y[k] * df_32[k]
                  - f_14 * ab_y[k] * df_37[k]
                  + f_15 * ab_y[k] * df_39[k]
                  + f_16 * ab_y[k] * df_52[k]
                  + f_16 * ab_y[k] * df_57[k]
                  - f_17 * ab_y[k] * df_59[k]
                  + f_14 * dg_17[k]
                  + f_14 * dg_22[k]
                  - f_15 * dg_24[k]
                  + f_14 * dg_49[k]
                  + f_14 * dg_56[k]
                  - f_15 * dg_58[k]
                  - f_16 * dg_79[k]
                  - f_16 * dg_86[k]
                  + f_17 * dg_88[k];
        g_23[k] = g_17[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_10, df_13, df_15, df_30, df_33, df_35, df_50, df_53, \
                         df_55, dg_15, dg_18, dg_20, dg_46, dg_51, dg_53, dg_76, dg_81, \
                         dg_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -0.375 * ab_x[k] * df_10[k]
                  - 0.375 * ab_x[k] * df_13[k]
                  + 1.5 * ab_x[k] * df_15[k]
                  - 0.375 * ab_y[k] * df_30[k]
                  - 0.375 * ab_y[k] * df_33[k]
                  + 1.5 * ab_y[k] * df_35[k]
                  + 1.5 * ab_y[k] * df_50[k]
                  + 1.5 * ab_y[k] * df_53[k]
                  - 6.0 * ab_y[k] * df_55[k]
                  + 0.375 * dg_15[k]
                  + 0.375 * dg_18[k]
                  - 1.5 * dg_20[k]
                  + 0.375 * dg_46[k]
                  + 0.375 * dg_51[k]
                  - 1.5 * dg_53[k]
                  - 1.5 * dg_76[k]
                  - 1.5 * dg_81[k]
                  + 6.0 * dg_83[k];
        g_30[k] = g_18[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_12, df_17, df_32, df_37, df_52, df_57, dg_17, dg_22, \
                         dg_49, dg_56, dg_79, dg_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_8 * ab_x[k] * df_12[k]
                  - f_8 * ab_x[k] * df_17[k]
                  + f_8 * ab_y[k] * df_32[k]
                  - f_8 * ab_y[k] * df_37[k]
                  - f_18 * ab_y[k] * df_52[k]
                  + f_18 * ab_y[k] * df_57[k]
                  - f_8 * dg_17[k]
                  + f_8 * dg_22[k]
                  - f_8 * dg_49[k]
                  + f_8 * dg_56[k]
                  + f_18 * dg_79[k]
                  - f_18 * dg_86[k];
        g_37[k] = g_19[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_10, df_13, df_30, df_33, df_50, df_53, dg_15, dg_18, \
                         dg_46, dg_51, dg_76, dg_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_4 * ab_x[k] * df_10[k]
                  - f_2 * ab_x[k] * df_13[k]
                  + f_4 * ab_y[k] * df_30[k]
                  - f_2 * ab_y[k] * df_33[k]
                  - f_5 * ab_y[k] * df_50[k]
                  + f_3 * ab_y[k] * df_53[k]
                  - f_4 * dg_15[k]
                  + f_2 * dg_18[k]
                  - f_4 * dg_46[k]
                  + f_2 * dg_51[k]
                  + f_5 * dg_76[k]
                  - f_3 * dg_81[k];
        g_44[k] = g_20[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, df_22, df_27, df_29, df_42, df_47, df_49, df_52, \
                         df_57, df_59, dg_32, dg_37, dg_39, dg_64, dg_71, dg_73, dg_80, dg_87, \
                         dg_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -2.25 * ab_x[k] * df_22[k]
                  - 2.25 * ab_x[k] * df_27[k]
                  + 1.5 * ab_x[k] * df_29[k]
                  - 2.25 * ab_y[k] * df_42[k]
                  - 2.25 * ab_y[k] * df_47[k]
                  + 1.5 * ab_y[k] * df_49[k]
                  + 1.5 * ab_z[k] * df_52[k]
                  + 1.5 * ab_z[k] * df_57[k]
                  - ab_z[k] * df_59[k]
                  + 2.25 * dg_32[k]
                  + 2.25 * dg_37[k]
                  - 1.5 * dg_39[k]
                  + 2.25 * dg_64[k]
                  + 2.25 * dg_71[k]
                  - 1.5 * dg_73[k]
                  - 1.5 * dg_80[k]
                  - 1.5 * dg_87[k]
                  + dg_89[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, df_20, df_23, df_25, df_40, df_43, df_45, df_50, \
                         df_53, df_55, dg_30, dg_33, dg_35, dg_61, dg_66, dg_68, dg_77, dg_82, \
                         dg_84 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_14 * ab_x[k] * df_20[k]
                  - f_14 * ab_x[k] * df_23[k]
                  + f_16 * ab_x[k] * df_25[k]
                  - f_14 * ab_y[k] * df_40[k]
                  - f_14 * ab_y[k] * df_43[k]
                  + f_16 * ab_y[k] * df_45[k]
                  + f_15 * ab_z[k] * df_50[k]
                  + f_15 * ab_z[k] * df_53[k]
                  - f_17 * ab_z[k] * df_55[k]
                  + f_14 * dg_30[k]
                  + f_14 * dg_33[k]
                  - f_16 * dg_35[k]
                  + f_14 * dg_61[k]
                  + f_14 * dg_66[k]
                  - f_16 * dg_68[k]
                  - f_15 * dg_77[k]
                  - f_15 * dg_82[k]
                  + f_17 * dg_84[k];
        g_31[k] = g_25[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, df_22, df_27, df_42, df_47, df_52, df_57, dg_32, \
                         dg_37, dg_64, dg_71, dg_80, dg_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = f_19 * ab_x[k] * df_22[k]
                  - f_19 * ab_x[k] * df_27[k]
                  + f_19 * ab_y[k] * df_42[k]
                  - f_19 * ab_y[k] * df_47[k]
                  - f_5 * ab_z[k] * df_52[k]
                  + f_5 * ab_z[k] * df_57[k]
                  - f_19 * dg_32[k]
                  + f_19 * dg_37[k]
                  - f_19 * dg_64[k]
                  + f_19 * dg_71[k]
                  + f_5 * dg_80[k]
                  - f_5 * dg_87[k];
        g_38[k] = g_26[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, df_20, df_23, df_40, df_43, df_50, df_53, dg_30, \
                         dg_33, dg_61, dg_66, dg_77, dg_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = f_8 * ab_x[k] * df_20[k]
                  - f_6 * ab_x[k] * df_23[k]
                  + f_8 * ab_y[k] * df_40[k]
                  - f_6 * ab_y[k] * df_43[k]
                  - f_9 * ab_z[k] * df_50[k]
                  + f_7 * ab_z[k] * df_53[k]
                  - f_8 * dg_30[k]
                  + f_6 * dg_33[k]
                  - f_8 * dg_61[k]
                  + f_6 * dg_66[k]
                  + f_9 * dg_77[k]
                  - f_7 * dg_82[k];
        g_45[k] = g_27[k];
    }

#pragma omp simd aligned(ab_x, df_0, df_3, df_5, df_30, df_33, df_35, df_50, df_53, df_55, \
                         dg_0, dg_3, dg_5, dg_45, dg_48, dg_50, dg_75, dg_78, \
                         dg_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -0.375 * ab_x[k] * df_0[k]
                  - 0.375 * ab_x[k] * df_3[k]
                  + 1.5 * ab_x[k] * df_5[k]
                  - 0.375 * ab_x[k] * df_30[k]
                  - 0.375 * ab_x[k] * df_33[k]
                  + 1.5 * ab_x[k] * df_35[k]
                  + 1.5 * ab_x[k] * df_50[k]
                  + 1.5 * ab_x[k] * df_53[k]
                  - 6.0 * ab_x[k] * df_55[k]
                  + 0.375 * dg_0[k]
                  + 0.375 * dg_3[k]
                  - 1.5 * dg_5[k]
                  + 0.375 * dg_45[k]
                  + 0.375 * dg_48[k]
                  - 1.5 * dg_50[k]
                  - 1.5 * dg_75[k]
                  - 1.5 * dg_78[k]
                  + 6.0 * dg_80[k];
    }

#pragma omp simd aligned(ab_x, df_2, df_7, df_32, df_37, df_52, df_57, dg_2, dg_7, dg_47, \
                         dg_52, dg_77, dg_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_8 * ab_x[k] * df_2[k]
                  - f_8 * ab_x[k] * df_7[k]
                  + f_8 * ab_x[k] * df_32[k]
                  - f_8 * ab_x[k] * df_37[k]
                  - f_18 * ab_x[k] * df_52[k]
                  + f_18 * ab_x[k] * df_57[k]
                  - f_8 * dg_2[k]
                  + f_8 * dg_7[k]
                  - f_8 * dg_47[k]
                  + f_8 * dg_52[k]
                  + f_18 * dg_77[k]
                  - f_18 * dg_82[k];
        g_39[k] = g_33[k];
    }

#pragma omp simd aligned(ab_x, df_0, df_3, df_30, df_33, df_50, df_53, dg_0, dg_3, dg_45, \
                         dg_48, dg_75, dg_78 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_4 * ab_x[k] * df_0[k]
                  - f_2 * ab_x[k] * df_3[k]
                  + f_4 * ab_x[k] * df_30[k]
                  - f_2 * ab_x[k] * df_33[k]
                  - f_5 * ab_x[k] * df_50[k]
                  + f_3 * ab_x[k] * df_53[k]
                  - f_4 * dg_0[k]
                  + f_2 * dg_3[k]
                  - f_4 * dg_45[k]
                  + f_2 * dg_48[k]
                  + f_5 * dg_75[k]
                  - f_3 * dg_78[k];
        g_46[k] = g_34[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_22, df_27, df_42, df_47, dg_32, dg_37, dg_64, \
                         dg_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -3.75 * ab_x[k] * df_22[k]
                  + 3.75 * ab_x[k] * df_27[k]
                  + 3.75 * ab_y[k] * df_42[k]
                  - 3.75 * ab_y[k] * df_47[k]
                  + 3.75 * dg_32[k]
                  - 3.75 * dg_37[k]
                  - 3.75 * dg_64[k]
                  + 3.75 * dg_71[k];
    }

#pragma omp simd aligned(ab_x, ab_y, df_20, df_23, df_40, df_43, dg_30, dg_33, dg_61, \
                         dg_66 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_11 * ab_x[k] * df_20[k]
                  + f_10 * ab_x[k] * df_23[k]
                  + f_11 * ab_y[k] * df_40[k]
                  - f_10 * ab_y[k] * df_43[k]
                  + f_11 * dg_30[k]
                  - f_10 * dg_33[k]
                  - f_11 * dg_61[k]
                  + f_10 * dg_66[k];
        g_47[k] = g_41[k];
    }

#pragma omp simd aligned(ab_x, df_0, df_3, df_30, df_33, dg_0, dg_3, dg_45, \
                         dg_48 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -0.625 * ab_x[k] * df_0[k]
                  + 1.875 * ab_x[k] * df_3[k]
                  + 1.875 * ab_x[k] * df_30[k]
                  - 5.625 * ab_x[k] * df_33[k]
                  + 0.625 * dg_0[k]
                  - 1.875 * dg_3[k]
                  - 1.875 * dg_45[k]
                  + 5.625 * dg_48[k];
    }
}

}  // namespace simdovl
