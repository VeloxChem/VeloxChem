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


#include "SimdTransformPK.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_pk(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t pk,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.21875 * std::sqrt(429.0);
    const auto f_1 = 1.09375 * std::sqrt(429.0);
    const auto f_2 = 0.65625 * std::sqrt(429.0);
    const auto f_3 = 0.03125 * std::sqrt(429.0);
    const auto f_4 = 0.1875 * std::sqrt(6006.0);
    const auto f_5 = 0.625 * std::sqrt(6006.0);
    const auto f_6 = 0.15625 * std::sqrt(231.0);
    const auto f_7 = 1.875 * std::sqrt(231.0);
    const auto f_8 = 0.28125 * std::sqrt(231.0);
    const auto f_9 = 3.75 * std::sqrt(231.0);
    const auto f_10 = 0.03125 * std::sqrt(231.0);
    const auto f_11 = 0.375 * std::sqrt(231.0);
    const auto f_12 = 0.75 * std::sqrt(231.0);
    const auto f_13 = 2.5 * std::sqrt(231.0);
    const auto f_14 = 0.28125 * std::sqrt(21.0);
    const auto f_15 = 0.46875 * std::sqrt(21.0);
    const auto f_16 = 5.625 * std::sqrt(21.0);
    const auto f_17 = 0.09375 * std::sqrt(21.0);
    const auto f_18 = 3.75 * std::sqrt(21.0);
    const auto f_19 = 7.5 * std::sqrt(21.0);
    const auto f_20 = 1.875 * std::sqrt(21.0);
    const auto f_21 = 2.5 * std::sqrt(21.0);
    const auto f_22 = 0.9375 * std::sqrt(42.0);
    const auto f_23 = 1.875 * std::sqrt(42.0);
    const auto f_24 = 5.0 * std::sqrt(42.0);
    const auto f_25 = 3.0 * std::sqrt(42.0);
    const auto f_26 = 0.15625 * std::sqrt(7.0);
    const auto f_27 = 0.46875 * std::sqrt(7.0);
    const auto f_28 = 3.75 * std::sqrt(7.0);
    const auto f_29 = 7.5 * std::sqrt(7.0);
    const auto f_30 = 2.0 * std::sqrt(7.0);
    const auto f_31 = 0.46875 * std::sqrt(42.0);
    const auto f_32 = 2.5 * std::sqrt(42.0);
    const auto f_33 = 1.5 * std::sqrt(42.0);
    const auto f_34 = 0.1875 * std::sqrt(231.0);
    const auto f_35 = 0.9375 * std::sqrt(231.0);
    const auto f_36 = 0.625 * std::sqrt(231.0);
    const auto f_37 = 0.03125 * std::sqrt(6006.0);
    const auto f_38 = 0.46875 * std::sqrt(6006.0);

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

    const auto *pk_0 = buffer.data(pk + 0);
    const auto *pk_1 = buffer.data(pk + 1);
    const auto *pk_2 = buffer.data(pk + 2);
    const auto *pk_3 = buffer.data(pk + 3);
    const auto *pk_4 = buffer.data(pk + 4);
    const auto *pk_5 = buffer.data(pk + 5);
    const auto *pk_6 = buffer.data(pk + 6);
    const auto *pk_7 = buffer.data(pk + 7);
    const auto *pk_8 = buffer.data(pk + 8);
    const auto *pk_9 = buffer.data(pk + 9);
    const auto *pk_10 = buffer.data(pk + 10);
    const auto *pk_11 = buffer.data(pk + 11);
    const auto *pk_12 = buffer.data(pk + 12);
    const auto *pk_13 = buffer.data(pk + 13);
    const auto *pk_14 = buffer.data(pk + 14);
    const auto *pk_15 = buffer.data(pk + 15);
    const auto *pk_16 = buffer.data(pk + 16);
    const auto *pk_17 = buffer.data(pk + 17);
    const auto *pk_18 = buffer.data(pk + 18);
    const auto *pk_19 = buffer.data(pk + 19);
    const auto *pk_20 = buffer.data(pk + 20);
    const auto *pk_21 = buffer.data(pk + 21);
    const auto *pk_22 = buffer.data(pk + 22);
    const auto *pk_23 = buffer.data(pk + 23);
    const auto *pk_24 = buffer.data(pk + 24);
    const auto *pk_25 = buffer.data(pk + 25);
    const auto *pk_26 = buffer.data(pk + 26);
    const auto *pk_27 = buffer.data(pk + 27);
    const auto *pk_28 = buffer.data(pk + 28);
    const auto *pk_29 = buffer.data(pk + 29);
    const auto *pk_30 = buffer.data(pk + 30);
    const auto *pk_31 = buffer.data(pk + 31);
    const auto *pk_32 = buffer.data(pk + 32);
    const auto *pk_33 = buffer.data(pk + 33);
    const auto *pk_34 = buffer.data(pk + 34);
    const auto *pk_35 = buffer.data(pk + 35);
    const auto *pk_36 = buffer.data(pk + 36);
    const auto *pk_37 = buffer.data(pk + 37);
    const auto *pk_38 = buffer.data(pk + 38);
    const auto *pk_39 = buffer.data(pk + 39);
    const auto *pk_40 = buffer.data(pk + 40);
    const auto *pk_41 = buffer.data(pk + 41);
    const auto *pk_42 = buffer.data(pk + 42);
    const auto *pk_43 = buffer.data(pk + 43);
    const auto *pk_44 = buffer.data(pk + 44);
    const auto *pk_45 = buffer.data(pk + 45);
    const auto *pk_46 = buffer.data(pk + 46);
    const auto *pk_47 = buffer.data(pk + 47);
    const auto *pk_48 = buffer.data(pk + 48);
    const auto *pk_49 = buffer.data(pk + 49);
    const auto *pk_50 = buffer.data(pk + 50);
    const auto *pk_51 = buffer.data(pk + 51);
    const auto *pk_52 = buffer.data(pk + 52);
    const auto *pk_53 = buffer.data(pk + 53);
    const auto *pk_54 = buffer.data(pk + 54);
    const auto *pk_55 = buffer.data(pk + 55);
    const auto *pk_56 = buffer.data(pk + 56);
    const auto *pk_57 = buffer.data(pk + 57);
    const auto *pk_58 = buffer.data(pk + 58);
    const auto *pk_59 = buffer.data(pk + 59);
    const auto *pk_60 = buffer.data(pk + 60);
    const auto *pk_61 = buffer.data(pk + 61);
    const auto *pk_62 = buffer.data(pk + 62);
    const auto *pk_63 = buffer.data(pk + 63);
    const auto *pk_64 = buffer.data(pk + 64);
    const auto *pk_65 = buffer.data(pk + 65);
    const auto *pk_66 = buffer.data(pk + 66);
    const auto *pk_67 = buffer.data(pk + 67);
    const auto *pk_68 = buffer.data(pk + 68);
    const auto *pk_69 = buffer.data(pk + 69);
    const auto *pk_70 = buffer.data(pk + 70);
    const auto *pk_71 = buffer.data(pk + 71);
    const auto *pk_72 = buffer.data(pk + 72);
    const auto *pk_73 = buffer.data(pk + 73);
    const auto *pk_74 = buffer.data(pk + 74);
    const auto *pk_75 = buffer.data(pk + 75);
    const auto *pk_76 = buffer.data(pk + 76);
    const auto *pk_77 = buffer.data(pk + 77);
    const auto *pk_78 = buffer.data(pk + 78);
    const auto *pk_79 = buffer.data(pk + 79);
    const auto *pk_80 = buffer.data(pk + 80);
    const auto *pk_81 = buffer.data(pk + 81);
    const auto *pk_82 = buffer.data(pk + 82);
    const auto *pk_83 = buffer.data(pk + 83);
    const auto *pk_84 = buffer.data(pk + 84);
    const auto *pk_85 = buffer.data(pk + 85);
    const auto *pk_86 = buffer.data(pk + 86);
    const auto *pk_87 = buffer.data(pk + 87);
    const auto *pk_88 = buffer.data(pk + 88);
    const auto *pk_89 = buffer.data(pk + 89);
    const auto *pk_90 = buffer.data(pk + 90);
    const auto *pk_91 = buffer.data(pk + 91);
    const auto *pk_92 = buffer.data(pk + 92);
    const auto *pk_93 = buffer.data(pk + 93);
    const auto *pk_94 = buffer.data(pk + 94);
    const auto *pk_95 = buffer.data(pk + 95);
    const auto *pk_96 = buffer.data(pk + 96);
    const auto *pk_97 = buffer.data(pk + 97);
    const auto *pk_98 = buffer.data(pk + 98);
    const auto *pk_99 = buffer.data(pk + 99);
    const auto *pk_100 = buffer.data(pk + 100);
    const auto *pk_101 = buffer.data(pk + 101);
    const auto *pk_102 = buffer.data(pk + 102);
    const auto *pk_103 = buffer.data(pk + 103);
    const auto *pk_104 = buffer.data(pk + 104);
    const auto *pk_105 = buffer.data(pk + 105);
    const auto *pk_106 = buffer.data(pk + 106);
    const auto *pk_107 = buffer.data(pk + 107);

#pragma omp simd aligned(pk_37, pk_40, pk_42, pk_44, pk_47, pk_49, pk_51, pk_53, pk_58, pk_60, \
                         pk_64, pk_66 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * pk_37[k]
                 - f_1 * pk_42[k]
                 + f_2 * pk_51[k]
                 - f_3 * pk_64[k];

        g_1[k] = f_4 * pk_40[k]
                 - f_5 * pk_47[k]
                 + f_4 * pk_58[k];

        g_2[k] = -f_6 * pk_37[k]
                 + f_6 * pk_42[k]
                 + f_7 * pk_44[k]
                 + f_8 * pk_51[k]
                 - f_9 * pk_53[k]
                 - f_10 * pk_64[k]
                 + f_11 * pk_66[k];

        g_3[k] = -f_12 * pk_40[k]
                 + f_13 * pk_49[k]
                 + f_12 * pk_58[k]
                 - f_13 * pk_60[k];
    }

#pragma omp simd aligned(pk_37, pk_42, pk_44, pk_51, pk_53, pk_55, pk_64, pk_66, \
                         pk_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_14 * pk_37[k]
                 + f_15 * pk_42[k]
                 - f_16 * pk_44[k]
                 + f_17 * pk_51[k]
                 - f_18 * pk_53[k]
                 + f_19 * pk_55[k]
                 - f_17 * pk_64[k]
                 + f_20 * pk_66[k]
                 - f_21 * pk_68[k];
    }

#pragma omp simd aligned(pk_40, pk_47, pk_49, pk_58, pk_60, pk_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_22 * pk_40[k]
                 + f_23 * pk_47[k]
                 - f_24 * pk_49[k]
                 + f_22 * pk_58[k]
                 - f_24 * pk_60[k]
                 + f_25 * pk_62[k];
    }

#pragma omp simd aligned(pk_37, pk_42, pk_44, pk_51, pk_53, pk_55, pk_64, pk_66, pk_68, \
                         pk_70 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_26 * pk_37[k]
                 - f_27 * pk_42[k]
                 + f_28 * pk_44[k]
                 - f_27 * pk_51[k]
                 + f_29 * pk_53[k]
                 - f_29 * pk_55[k]
                 - f_26 * pk_64[k]
                 + f_28 * pk_66[k]
                 - f_29 * pk_68[k]
                 + f_30 * pk_70[k];
    }

#pragma omp simd aligned(pk_38, pk_43, pk_45, pk_52, pk_54, pk_56, pk_65, pk_67, pk_69, \
                         pk_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -2.1875 * pk_38[k]
                 - 6.5625 * pk_43[k]
                 + 13.125 * pk_45[k]
                 - 6.5625 * pk_52[k]
                 + 26.25 * pk_54[k]
                 - 10.5 * pk_56[k]
                 - 2.1875 * pk_65[k]
                 + 13.125 * pk_67[k]
                 - 10.5 * pk_69[k]
                 + pk_71[k];
    }

#pragma omp simd aligned(pk_36, pk_39, pk_41, pk_46, pk_48, pk_50, pk_57, pk_59, pk_61, \
                         pk_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_26 * pk_36[k]
                 - f_27 * pk_39[k]
                 + f_28 * pk_41[k]
                 - f_27 * pk_46[k]
                 + f_29 * pk_48[k]
                 - f_29 * pk_50[k]
                 - f_26 * pk_57[k]
                 + f_28 * pk_59[k]
                 - f_29 * pk_61[k]
                 + f_30 * pk_63[k];
    }

#pragma omp simd aligned(pk_38, pk_43, pk_45, pk_52, pk_56, pk_65, pk_67, \
                         pk_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_31 * pk_38[k]
                 + f_31 * pk_43[k]
                 - f_32 * pk_45[k]
                 - f_31 * pk_52[k]
                 + f_33 * pk_56[k]
                 - f_31 * pk_65[k]
                 + f_32 * pk_67[k]
                 - f_33 * pk_69[k];
    }

#pragma omp simd aligned(pk_36, pk_39, pk_41, pk_46, pk_48, pk_50, pk_57, pk_59, \
                         pk_61 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = f_17 * pk_36[k]
                  - f_17 * pk_39[k]
                  - f_20 * pk_41[k]
                  - f_15 * pk_46[k]
                  + f_18 * pk_48[k]
                  + f_21 * pk_50[k]
                  - f_14 * pk_57[k]
                  + f_16 * pk_59[k]
                  - f_19 * pk_61[k];
    }

#pragma omp simd aligned(pk_36, pk_38, pk_39, pk_41, pk_43, pk_45, pk_46, pk_48, pk_52, pk_54, \
                         pk_57, pk_59, pk_65, pk_67 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_34 * pk_38[k]
                  + f_35 * pk_43[k]
                  + f_36 * pk_45[k]
                  + f_35 * pk_52[k]
                  - f_9 * pk_54[k]
                  - f_34 * pk_65[k]
                  + f_36 * pk_67[k];

        g_12[k] = -f_10 * pk_36[k]
                  + f_8 * pk_39[k]
                  + f_11 * pk_41[k]
                  + f_6 * pk_46[k]
                  - f_9 * pk_48[k]
                  - f_6 * pk_57[k]
                  + f_7 * pk_59[k];
    }

#pragma omp simd aligned(pk_36, pk_38, pk_39, pk_43, pk_46, pk_52, pk_57, pk_65, pk_73, pk_78, \
                         pk_87, pk_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_37 * pk_38[k]
                  - f_38 * pk_43[k]
                  + f_38 * pk_52[k]
                  - f_37 * pk_65[k];

        g_14[k] = f_3 * pk_36[k]
                  - f_2 * pk_39[k]
                  + f_1 * pk_46[k]
                  - f_0 * pk_57[k];

        g_15[k] = f_0 * pk_73[k]
                  - f_1 * pk_78[k]
                  + f_2 * pk_87[k]
                  - f_3 * pk_100[k];
    }

#pragma omp simd aligned(pk_73, pk_76, pk_78, pk_80, pk_83, pk_85, pk_87, pk_89, pk_94, pk_96, \
                         pk_100, pk_102 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_4 * pk_76[k]
                  - f_5 * pk_83[k]
                  + f_4 * pk_94[k];

        g_17[k] = -f_6 * pk_73[k]
                  + f_6 * pk_78[k]
                  + f_7 * pk_80[k]
                  + f_8 * pk_87[k]
                  - f_9 * pk_89[k]
                  - f_10 * pk_100[k]
                  + f_11 * pk_102[k];

        g_18[k] = -f_12 * pk_76[k]
                  + f_13 * pk_85[k]
                  + f_12 * pk_94[k]
                  - f_13 * pk_96[k];
    }

#pragma omp simd aligned(pk_73, pk_78, pk_80, pk_87, pk_89, pk_91, pk_100, pk_102, \
                         pk_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_14 * pk_73[k]
                  + f_15 * pk_78[k]
                  - f_16 * pk_80[k]
                  + f_17 * pk_87[k]
                  - f_18 * pk_89[k]
                  + f_19 * pk_91[k]
                  - f_17 * pk_100[k]
                  + f_20 * pk_102[k]
                  - f_21 * pk_104[k];
    }

#pragma omp simd aligned(pk_76, pk_83, pk_85, pk_94, pk_96, pk_98 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_22 * pk_76[k]
                  + f_23 * pk_83[k]
                  - f_24 * pk_85[k]
                  + f_22 * pk_94[k]
                  - f_24 * pk_96[k]
                  + f_25 * pk_98[k];
    }

#pragma omp simd aligned(pk_73, pk_78, pk_80, pk_87, pk_89, pk_91, pk_100, pk_102, pk_104, \
                         pk_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_26 * pk_73[k]
                  - f_27 * pk_78[k]
                  + f_28 * pk_80[k]
                  - f_27 * pk_87[k]
                  + f_29 * pk_89[k]
                  - f_29 * pk_91[k]
                  - f_26 * pk_100[k]
                  + f_28 * pk_102[k]
                  - f_29 * pk_104[k]
                  + f_30 * pk_106[k];
    }

#pragma omp simd aligned(pk_74, pk_79, pk_81, pk_88, pk_90, pk_92, pk_101, pk_103, pk_105, \
                         pk_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -2.1875 * pk_74[k]
                  - 6.5625 * pk_79[k]
                  + 13.125 * pk_81[k]
                  - 6.5625 * pk_88[k]
                  + 26.25 * pk_90[k]
                  - 10.5 * pk_92[k]
                  - 2.1875 * pk_101[k]
                  + 13.125 * pk_103[k]
                  - 10.5 * pk_105[k]
                  + pk_107[k];
    }

#pragma omp simd aligned(pk_72, pk_75, pk_77, pk_82, pk_84, pk_86, pk_93, pk_95, pk_97, \
                         pk_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_26 * pk_72[k]
                  - f_27 * pk_75[k]
                  + f_28 * pk_77[k]
                  - f_27 * pk_82[k]
                  + f_29 * pk_84[k]
                  - f_29 * pk_86[k]
                  - f_26 * pk_93[k]
                  + f_28 * pk_95[k]
                  - f_29 * pk_97[k]
                  + f_30 * pk_99[k];
    }

#pragma omp simd aligned(pk_74, pk_79, pk_81, pk_88, pk_92, pk_101, pk_103, \
                         pk_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_31 * pk_74[k]
                  + f_31 * pk_79[k]
                  - f_32 * pk_81[k]
                  - f_31 * pk_88[k]
                  + f_33 * pk_92[k]
                  - f_31 * pk_101[k]
                  + f_32 * pk_103[k]
                  - f_33 * pk_105[k];
    }

#pragma omp simd aligned(pk_72, pk_75, pk_77, pk_82, pk_84, pk_86, pk_93, pk_95, \
                         pk_97 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_17 * pk_72[k]
                  - f_17 * pk_75[k]
                  - f_20 * pk_77[k]
                  - f_15 * pk_82[k]
                  + f_18 * pk_84[k]
                  + f_21 * pk_86[k]
                  - f_14 * pk_93[k]
                  + f_16 * pk_95[k]
                  - f_19 * pk_97[k];
    }

#pragma omp simd aligned(pk_72, pk_74, pk_75, pk_77, pk_79, pk_81, pk_82, pk_84, pk_88, pk_90, \
                         pk_93, pk_95, pk_101, pk_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_34 * pk_74[k]
                  + f_35 * pk_79[k]
                  + f_36 * pk_81[k]
                  + f_35 * pk_88[k]
                  - f_9 * pk_90[k]
                  - f_34 * pk_101[k]
                  + f_36 * pk_103[k];

        g_27[k] = -f_10 * pk_72[k]
                  + f_8 * pk_75[k]
                  + f_11 * pk_77[k]
                  + f_6 * pk_82[k]
                  - f_9 * pk_84[k]
                  - f_6 * pk_93[k]
                  + f_7 * pk_95[k];
    }

#pragma omp simd aligned(pk_1, pk_6, pk_15, pk_28, pk_72, pk_74, pk_75, pk_79, pk_82, pk_88, \
                         pk_93, pk_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_37 * pk_74[k]
                  - f_38 * pk_79[k]
                  + f_38 * pk_88[k]
                  - f_37 * pk_101[k];

        g_29[k] = f_3 * pk_72[k]
                  - f_2 * pk_75[k]
                  + f_1 * pk_82[k]
                  - f_0 * pk_93[k];

        g_30[k] = f_0 * pk_1[k]
                  - f_1 * pk_6[k]
                  + f_2 * pk_15[k]
                  - f_3 * pk_28[k];
    }

#pragma omp simd aligned(pk_1, pk_4, pk_6, pk_8, pk_11, pk_13, pk_15, pk_17, pk_22, pk_24, \
                         pk_28, pk_30 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = f_4 * pk_4[k]
                  - f_5 * pk_11[k]
                  + f_4 * pk_22[k];

        g_32[k] = -f_6 * pk_1[k]
                  + f_6 * pk_6[k]
                  + f_7 * pk_8[k]
                  + f_8 * pk_15[k]
                  - f_9 * pk_17[k]
                  - f_10 * pk_28[k]
                  + f_11 * pk_30[k];

        g_33[k] = -f_12 * pk_4[k]
                  + f_13 * pk_13[k]
                  + f_12 * pk_22[k]
                  - f_13 * pk_24[k];
    }

#pragma omp simd aligned(pk_1, pk_6, pk_8, pk_15, pk_17, pk_19, pk_28, pk_30, \
                         pk_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_14 * pk_1[k]
                  + f_15 * pk_6[k]
                  - f_16 * pk_8[k]
                  + f_17 * pk_15[k]
                  - f_18 * pk_17[k]
                  + f_19 * pk_19[k]
                  - f_17 * pk_28[k]
                  + f_20 * pk_30[k]
                  - f_21 * pk_32[k];
    }

#pragma omp simd aligned(pk_4, pk_11, pk_13, pk_22, pk_24, pk_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_22 * pk_4[k]
                  + f_23 * pk_11[k]
                  - f_24 * pk_13[k]
                  + f_22 * pk_22[k]
                  - f_24 * pk_24[k]
                  + f_25 * pk_26[k];
    }

#pragma omp simd aligned(pk_1, pk_6, pk_8, pk_15, pk_17, pk_19, pk_28, pk_30, pk_32, \
                         pk_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = -f_26 * pk_1[k]
                  - f_27 * pk_6[k]
                  + f_28 * pk_8[k]
                  - f_27 * pk_15[k]
                  + f_29 * pk_17[k]
                  - f_29 * pk_19[k]
                  - f_26 * pk_28[k]
                  + f_28 * pk_30[k]
                  - f_29 * pk_32[k]
                  + f_30 * pk_34[k];
    }

#pragma omp simd aligned(pk_2, pk_7, pk_9, pk_16, pk_18, pk_20, pk_29, pk_31, pk_33, \
                         pk_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -2.1875 * pk_2[k]
                  - 6.5625 * pk_7[k]
                  + 13.125 * pk_9[k]
                  - 6.5625 * pk_16[k]
                  + 26.25 * pk_18[k]
                  - 10.5 * pk_20[k]
                  - 2.1875 * pk_29[k]
                  + 13.125 * pk_31[k]
                  - 10.5 * pk_33[k]
                  + pk_35[k];
    }

#pragma omp simd aligned(pk_0, pk_3, pk_5, pk_10, pk_12, pk_14, pk_21, pk_23, pk_25, \
                         pk_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_26 * pk_0[k]
                  - f_27 * pk_3[k]
                  + f_28 * pk_5[k]
                  - f_27 * pk_10[k]
                  + f_29 * pk_12[k]
                  - f_29 * pk_14[k]
                  - f_26 * pk_21[k]
                  + f_28 * pk_23[k]
                  - f_29 * pk_25[k]
                  + f_30 * pk_27[k];
    }

#pragma omp simd aligned(pk_2, pk_7, pk_9, pk_16, pk_20, pk_29, pk_31, \
                         pk_33 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = f_31 * pk_2[k]
                  + f_31 * pk_7[k]
                  - f_32 * pk_9[k]
                  - f_31 * pk_16[k]
                  + f_33 * pk_20[k]
                  - f_31 * pk_29[k]
                  + f_32 * pk_31[k]
                  - f_33 * pk_33[k];
    }

#pragma omp simd aligned(pk_0, pk_3, pk_5, pk_10, pk_12, pk_14, pk_21, pk_23, \
                         pk_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_17 * pk_0[k]
                  - f_17 * pk_3[k]
                  - f_20 * pk_5[k]
                  - f_15 * pk_10[k]
                  + f_18 * pk_12[k]
                  + f_21 * pk_14[k]
                  - f_14 * pk_21[k]
                  + f_16 * pk_23[k]
                  - f_19 * pk_25[k];
    }

#pragma omp simd aligned(pk_0, pk_2, pk_3, pk_5, pk_7, pk_9, pk_10, pk_12, pk_16, pk_18, \
                         pk_21, pk_23, pk_29, pk_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_34 * pk_2[k]
                  + f_35 * pk_7[k]
                  + f_36 * pk_9[k]
                  + f_35 * pk_16[k]
                  - f_9 * pk_18[k]
                  - f_34 * pk_29[k]
                  + f_36 * pk_31[k];

        g_42[k] = -f_10 * pk_0[k]
                  + f_8 * pk_3[k]
                  + f_11 * pk_5[k]
                  + f_6 * pk_10[k]
                  - f_9 * pk_12[k]
                  - f_6 * pk_21[k]
                  + f_7 * pk_23[k];
    }

#pragma omp simd aligned(pk_0, pk_2, pk_3, pk_7, pk_10, pk_16, pk_21, \
                         pk_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_37 * pk_2[k]
                  - f_38 * pk_7[k]
                  + f_38 * pk_16[k]
                  - f_37 * pk_29[k];

        g_44[k] = f_3 * pk_0[k]
                  - f_2 * pk_3[k]
                  + f_1 * pk_10[k]
                  - f_0 * pk_21[k];
    }
}

}  // namespace simdtrf
