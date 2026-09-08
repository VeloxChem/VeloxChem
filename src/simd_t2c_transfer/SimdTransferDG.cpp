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


#include "SimdTransferDG.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_dg(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t pg, const size_t ph, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_8 = buffer.data(ph + 8);
    const auto *ph_9 = buffer.data(ph + 9);
    const auto *ph_10 = buffer.data(ph + 10);
    const auto *ph_11 = buffer.data(ph + 11);
    const auto *ph_12 = buffer.data(ph + 12);
    const auto *ph_13 = buffer.data(ph + 13);
    const auto *ph_14 = buffer.data(ph + 14);
    const auto *ph_21 = buffer.data(ph + 21);
    const auto *ph_22 = buffer.data(ph + 22);
    const auto *ph_23 = buffer.data(ph + 23);
    const auto *ph_24 = buffer.data(ph + 24);
    const auto *ph_25 = buffer.data(ph + 25);
    const auto *ph_26 = buffer.data(ph + 26);
    const auto *ph_27 = buffer.data(ph + 27);
    const auto *ph_28 = buffer.data(ph + 28);
    const auto *ph_29 = buffer.data(ph + 29);
    const auto *ph_30 = buffer.data(ph + 30);
    const auto *ph_31 = buffer.data(ph + 31);
    const auto *ph_32 = buffer.data(ph + 32);
    const auto *ph_33 = buffer.data(ph + 33);
    const auto *ph_34 = buffer.data(ph + 34);
    const auto *ph_35 = buffer.data(ph + 35);
    const auto *ph_36 = buffer.data(ph + 36);
    const auto *ph_37 = buffer.data(ph + 37);
    const auto *ph_38 = buffer.data(ph + 38);
    const auto *ph_39 = buffer.data(ph + 39);
    const auto *ph_40 = buffer.data(ph + 40);
    const auto *ph_42 = buffer.data(ph + 42);
    const auto *ph_43 = buffer.data(ph + 43);
    const auto *ph_44 = buffer.data(ph + 44);
    const auto *ph_45 = buffer.data(ph + 45);
    const auto *ph_46 = buffer.data(ph + 46);
    const auto *ph_47 = buffer.data(ph + 47);
    const auto *ph_48 = buffer.data(ph + 48);
    const auto *ph_49 = buffer.data(ph + 49);
    const auto *ph_50 = buffer.data(ph + 50);
    const auto *ph_51 = buffer.data(ph + 51);
    const auto *ph_52 = buffer.data(ph + 52);
    const auto *ph_53 = buffer.data(ph + 53);
    const auto *ph_54 = buffer.data(ph + 54);
    const auto *ph_55 = buffer.data(ph + 55);
    const auto *ph_56 = buffer.data(ph + 56);
    const auto *ph_57 = buffer.data(ph + 57);
    const auto *ph_58 = buffer.data(ph + 58);
    const auto *ph_59 = buffer.data(ph + 59);
    const auto *ph_60 = buffer.data(ph + 60);
    const auto *ph_61 = buffer.data(ph + 61);
    const auto *ph_62 = buffer.data(ph + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, pg_0, pg_1, pg_2, pg_3, pg_4, ph_0, \
                         ph_1, ph_2, ph_3, ph_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * pg_0[k]
                 + ph_0[k];

        t_1[k] = -ab_x[k] * pg_1[k]
                 + ph_1[k];

        t_2[k] = -ab_x[k] * pg_2[k]
                 + ph_2[k];

        t_3[k] = -ab_x[k] * pg_3[k]
                 + ph_3[k];

        t_4[k] = -ab_x[k] * pg_4[k]
                 + ph_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, pg_5, pg_6, pg_7, pg_8, pg_9, ph_5, \
                         ph_6, ph_7, ph_8, ph_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * pg_5[k]
                 + ph_5[k];

        t_6[k] = -ab_x[k] * pg_6[k]
                 + ph_6[k];

        t_7[k] = -ab_x[k] * pg_7[k]
                 + ph_7[k];

        t_8[k] = -ab_x[k] * pg_8[k]
                 + ph_8[k];

        t_9[k] = -ab_x[k] * pg_9[k]
                 + ph_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, pg_10, pg_11, pg_12, pg_13, \
                         pg_14, ph_10, ph_11, ph_12, ph_13, ph_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * pg_10[k]
                  + ph_10[k];

        t_11[k] = -ab_x[k] * pg_11[k]
                  + ph_11[k];

        t_12[k] = -ab_x[k] * pg_12[k]
                  + ph_12[k];

        t_13[k] = -ab_x[k] * pg_13[k]
                  + ph_13[k];

        t_14[k] = -ab_x[k] * pg_14[k]
                  + ph_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, pg_15, pg_16, pg_17, pg_18, \
                         pg_19, ph_21, ph_22, ph_23, ph_24, ph_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * pg_15[k]
                  + ph_21[k];

        t_16[k] = -ab_x[k] * pg_16[k]
                  + ph_22[k];

        t_17[k] = -ab_x[k] * pg_17[k]
                  + ph_23[k];

        t_18[k] = -ab_x[k] * pg_18[k]
                  + ph_24[k];

        t_19[k] = -ab_x[k] * pg_19[k]
                  + ph_25[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, pg_20, pg_21, pg_22, pg_23, \
                         pg_24, ph_26, ph_27, ph_28, ph_29, ph_30 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * pg_20[k]
                  + ph_26[k];

        t_21[k] = -ab_x[k] * pg_21[k]
                  + ph_27[k];

        t_22[k] = -ab_x[k] * pg_22[k]
                  + ph_28[k];

        t_23[k] = -ab_x[k] * pg_23[k]
                  + ph_29[k];

        t_24[k] = -ab_x[k] * pg_24[k]
                  + ph_30[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, pg_25, pg_26, pg_27, pg_28, \
                         pg_29, ph_31, ph_32, ph_33, ph_34, ph_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * pg_25[k]
                  + ph_31[k];

        t_26[k] = -ab_x[k] * pg_26[k]
                  + ph_32[k];

        t_27[k] = -ab_x[k] * pg_27[k]
                  + ph_33[k];

        t_28[k] = -ab_x[k] * pg_28[k]
                  + ph_34[k];

        t_29[k] = -ab_x[k] * pg_29[k]
                  + ph_35[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, pg_30, pg_31, pg_32, pg_33, \
                         pg_34, ph_42, ph_43, ph_44, ph_45, ph_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_x[k] * pg_30[k]
                  + ph_42[k];

        t_31[k] = -ab_x[k] * pg_31[k]
                  + ph_43[k];

        t_32[k] = -ab_x[k] * pg_32[k]
                  + ph_44[k];

        t_33[k] = -ab_x[k] * pg_33[k]
                  + ph_45[k];

        t_34[k] = -ab_x[k] * pg_34[k]
                  + ph_46[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, pg_35, pg_36, pg_37, pg_38, \
                         pg_39, ph_47, ph_48, ph_49, ph_50, ph_51 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_x[k] * pg_35[k]
                  + ph_47[k];

        t_36[k] = -ab_x[k] * pg_36[k]
                  + ph_48[k];

        t_37[k] = -ab_x[k] * pg_37[k]
                  + ph_49[k];

        t_38[k] = -ab_x[k] * pg_38[k]
                  + ph_50[k];

        t_39[k] = -ab_x[k] * pg_39[k]
                  + ph_51[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, pg_40, pg_41, pg_42, pg_43, \
                         pg_44, ph_52, ph_53, ph_54, ph_55, ph_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_x[k] * pg_40[k]
                  + ph_52[k];

        t_41[k] = -ab_x[k] * pg_41[k]
                  + ph_53[k];

        t_42[k] = -ab_x[k] * pg_42[k]
                  + ph_54[k];

        t_43[k] = -ab_x[k] * pg_43[k]
                  + ph_55[k];

        t_44[k] = -ab_x[k] * pg_44[k]
                  + ph_56[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_y, pg_15, pg_16, pg_17, pg_18, \
                         pg_19, ph_22, ph_24, ph_25, ph_27, ph_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = -ab_y[k] * pg_15[k]
                  + ph_22[k];

        t_46[k] = -ab_y[k] * pg_16[k]
                  + ph_24[k];

        t_47[k] = -ab_y[k] * pg_17[k]
                  + ph_25[k];

        t_48[k] = -ab_y[k] * pg_18[k]
                  + ph_27[k];

        t_49[k] = -ab_y[k] * pg_19[k]
                  + ph_28[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_y, pg_20, pg_21, pg_22, pg_23, \
                         pg_24, ph_29, ph_31, ph_32, ph_33, ph_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = -ab_y[k] * pg_20[k]
                  + ph_29[k];

        t_51[k] = -ab_y[k] * pg_21[k]
                  + ph_31[k];

        t_52[k] = -ab_y[k] * pg_22[k]
                  + ph_32[k];

        t_53[k] = -ab_y[k] * pg_23[k]
                  + ph_33[k];

        t_54[k] = -ab_y[k] * pg_24[k]
                  + ph_34[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_y, pg_25, pg_26, pg_27, pg_28, \
                         pg_29, ph_36, ph_37, ph_38, ph_39, ph_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = -ab_y[k] * pg_25[k]
                  + ph_36[k];

        t_56[k] = -ab_y[k] * pg_26[k]
                  + ph_37[k];

        t_57[k] = -ab_y[k] * pg_27[k]
                  + ph_38[k];

        t_58[k] = -ab_y[k] * pg_28[k]
                  + ph_39[k];

        t_59[k] = -ab_y[k] * pg_29[k]
                  + ph_40[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_y, pg_30, pg_31, pg_32, pg_33, \
                         pg_34, ph_43, ph_45, ph_46, ph_48, ph_49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = -ab_y[k] * pg_30[k]
                  + ph_43[k];

        t_61[k] = -ab_y[k] * pg_31[k]
                  + ph_45[k];

        t_62[k] = -ab_y[k] * pg_32[k]
                  + ph_46[k];

        t_63[k] = -ab_y[k] * pg_33[k]
                  + ph_48[k];

        t_64[k] = -ab_y[k] * pg_34[k]
                  + ph_49[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_y, pg_35, pg_36, pg_37, pg_38, \
                         pg_39, ph_50, ph_52, ph_53, ph_54, ph_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = -ab_y[k] * pg_35[k]
                  + ph_50[k];

        t_66[k] = -ab_y[k] * pg_36[k]
                  + ph_52[k];

        t_67[k] = -ab_y[k] * pg_37[k]
                  + ph_53[k];

        t_68[k] = -ab_y[k] * pg_38[k]
                  + ph_54[k];

        t_69[k] = -ab_y[k] * pg_39[k]
                  + ph_55[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_y, pg_40, pg_41, pg_42, pg_43, \
                         pg_44, ph_57, ph_58, ph_59, ph_60, ph_61 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = -ab_y[k] * pg_40[k]
                  + ph_57[k];

        t_71[k] = -ab_y[k] * pg_41[k]
                  + ph_58[k];

        t_72[k] = -ab_y[k] * pg_42[k]
                  + ph_59[k];

        t_73[k] = -ab_y[k] * pg_43[k]
                  + ph_60[k];

        t_74[k] = -ab_y[k] * pg_44[k]
                  + ph_61[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_z, pg_30, pg_31, pg_32, pg_33, \
                         pg_34, ph_44, ph_46, ph_47, ph_49, ph_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = -ab_z[k] * pg_30[k]
                  + ph_44[k];

        t_76[k] = -ab_z[k] * pg_31[k]
                  + ph_46[k];

        t_77[k] = -ab_z[k] * pg_32[k]
                  + ph_47[k];

        t_78[k] = -ab_z[k] * pg_33[k]
                  + ph_49[k];

        t_79[k] = -ab_z[k] * pg_34[k]
                  + ph_50[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_z, pg_35, pg_36, pg_37, pg_38, \
                         pg_39, ph_51, ph_53, ph_54, ph_55, ph_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = -ab_z[k] * pg_35[k]
                  + ph_51[k];

        t_81[k] = -ab_z[k] * pg_36[k]
                  + ph_53[k];

        t_82[k] = -ab_z[k] * pg_37[k]
                  + ph_54[k];

        t_83[k] = -ab_z[k] * pg_38[k]
                  + ph_55[k];

        t_84[k] = -ab_z[k] * pg_39[k]
                  + ph_56[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_z, pg_40, pg_41, pg_42, pg_43, \
                         pg_44, ph_58, ph_59, ph_60, ph_61, ph_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = -ab_z[k] * pg_40[k]
                  + ph_58[k];

        t_86[k] = -ab_z[k] * pg_41[k]
                  + ph_59[k];

        t_87[k] = -ab_z[k] * pg_42[k]
                  + ph_60[k];

        t_88[k] = -ab_z[k] * pg_43[k]
                  + ph_61[k];

        t_89[k] = -ab_z[k] * pg_44[k]
                  + ph_62[k];
    }
}

}  // namespace simdtrf
