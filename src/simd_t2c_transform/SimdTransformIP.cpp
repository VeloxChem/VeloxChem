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


#include "SimdTransformIP.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_ip(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t ip,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.1875 * std::sqrt(462.0);
    const auto f_1 = 0.625 * std::sqrt(462.0);
    const auto f_2 = 0.9375 * std::sqrt(154.0);
    const auto f_3 = 1.875 * std::sqrt(154.0);
    const auto f_4 = 0.1875 * std::sqrt(154.0);
    const auto f_5 = 0.75 * std::sqrt(7.0);
    const auto f_6 = 7.5 * std::sqrt(7.0);
    const auto f_7 = 0.5625 * std::sqrt(210.0);
    const auto f_8 = 0.375 * std::sqrt(210.0);
    const auto f_9 = 1.5 * std::sqrt(210.0);
    const auto f_10 = 0.1875 * std::sqrt(210.0);
    const auto f_11 = 0.5 * std::sqrt(210.0);
    const auto f_12 = 0.0625 * std::sqrt(210.0);
    const auto f_13 = 0.125 * std::sqrt(210.0);
    const auto f_14 = std::sqrt(210.0);
    const auto f_15 = 0.625 * std::sqrt(21.0);
    const auto f_16 = 1.25 * std::sqrt(21.0);
    const auto f_17 = 2.5 * std::sqrt(21.0);
    const auto f_18 = std::sqrt(21.0);
    const auto f_19 = 0.03125 * std::sqrt(210.0);
    const auto f_20 = 0.1875 * std::sqrt(7.0);
    const auto f_21 = 0.9375 * std::sqrt(7.0);
    const auto f_22 = 1.875 * std::sqrt(7.0);
    const auto f_23 = 11.25 * std::sqrt(7.0);
    const auto f_24 = 0.03125 * std::sqrt(462.0);
    const auto f_25 = 0.46875 * std::sqrt(462.0);

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

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_21 = buffer.data(ip + 21);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_27 = buffer.data(ip + 27);
    const auto *ip_28 = buffer.data(ip + 28);
    const auto *ip_29 = buffer.data(ip + 29);
    const auto *ip_30 = buffer.data(ip + 30);
    const auto *ip_31 = buffer.data(ip + 31);
    const auto *ip_32 = buffer.data(ip + 32);
    const auto *ip_33 = buffer.data(ip + 33);
    const auto *ip_34 = buffer.data(ip + 34);
    const auto *ip_35 = buffer.data(ip + 35);
    const auto *ip_36 = buffer.data(ip + 36);
    const auto *ip_37 = buffer.data(ip + 37);
    const auto *ip_38 = buffer.data(ip + 38);
    const auto *ip_39 = buffer.data(ip + 39);
    const auto *ip_40 = buffer.data(ip + 40);
    const auto *ip_41 = buffer.data(ip + 41);
    const auto *ip_42 = buffer.data(ip + 42);
    const auto *ip_43 = buffer.data(ip + 43);
    const auto *ip_44 = buffer.data(ip + 44);
    const auto *ip_45 = buffer.data(ip + 45);
    const auto *ip_46 = buffer.data(ip + 46);
    const auto *ip_47 = buffer.data(ip + 47);
    const auto *ip_48 = buffer.data(ip + 48);
    const auto *ip_49 = buffer.data(ip + 49);
    const auto *ip_50 = buffer.data(ip + 50);
    const auto *ip_51 = buffer.data(ip + 51);
    const auto *ip_52 = buffer.data(ip + 52);
    const auto *ip_53 = buffer.data(ip + 53);
    const auto *ip_54 = buffer.data(ip + 54);
    const auto *ip_55 = buffer.data(ip + 55);
    const auto *ip_56 = buffer.data(ip + 56);
    const auto *ip_57 = buffer.data(ip + 57);
    const auto *ip_58 = buffer.data(ip + 58);
    const auto *ip_59 = buffer.data(ip + 59);
    const auto *ip_60 = buffer.data(ip + 60);
    const auto *ip_61 = buffer.data(ip + 61);
    const auto *ip_62 = buffer.data(ip + 62);
    const auto *ip_63 = buffer.data(ip + 63);
    const auto *ip_64 = buffer.data(ip + 64);
    const auto *ip_65 = buffer.data(ip + 65);
    const auto *ip_66 = buffer.data(ip + 66);
    const auto *ip_67 = buffer.data(ip + 67);
    const auto *ip_68 = buffer.data(ip + 68);
    const auto *ip_69 = buffer.data(ip + 69);
    const auto *ip_70 = buffer.data(ip + 70);
    const auto *ip_71 = buffer.data(ip + 71);
    const auto *ip_72 = buffer.data(ip + 72);
    const auto *ip_73 = buffer.data(ip + 73);
    const auto *ip_74 = buffer.data(ip + 74);
    const auto *ip_75 = buffer.data(ip + 75);
    const auto *ip_76 = buffer.data(ip + 76);
    const auto *ip_77 = buffer.data(ip + 77);
    const auto *ip_78 = buffer.data(ip + 78);
    const auto *ip_79 = buffer.data(ip + 79);
    const auto *ip_80 = buffer.data(ip + 80);
    const auto *ip_81 = buffer.data(ip + 81);
    const auto *ip_82 = buffer.data(ip + 82);
    const auto *ip_83 = buffer.data(ip + 83);

#pragma omp simd aligned(ip_3, ip_4, ip_5, ip_13, ip_18, ip_19, ip_20, ip_34, ip_45, ip_46, \
                         ip_47, ip_67 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * ip_4[k]
                 - f_1 * ip_19[k]
                 + f_0 * ip_46[k];

        g_1[k] = f_0 * ip_5[k]
                 - f_1 * ip_20[k]
                 + f_0 * ip_47[k];

        g_2[k] = f_0 * ip_3[k]
                 - f_1 * ip_18[k]
                 + f_0 * ip_45[k];

        g_3[k] = f_2 * ip_13[k]
                 - f_3 * ip_34[k]
                 + f_4 * ip_67[k];
    }

#pragma omp simd aligned(ip_4, ip_12, ip_14, ip_25, ip_33, ip_35, ip_46, ip_52, ip_66, \
                         ip_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_2 * ip_14[k]
                 - f_3 * ip_35[k]
                 + f_4 * ip_68[k];

        g_5[k] = f_2 * ip_12[k]
                 - f_3 * ip_33[k]
                 + f_4 * ip_66[k];

        g_6[k] = -f_5 * ip_4[k]
                 + f_6 * ip_25[k]
                 + f_5 * ip_46[k]
                 - f_6 * ip_52[k];
    }

#pragma omp simd aligned(ip_3, ip_5, ip_13, ip_24, ip_26, ip_34, ip_40, ip_45, ip_47, ip_51, \
                         ip_53, ip_67, ip_73 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_5 * ip_5[k]
                 + f_6 * ip_26[k]
                 + f_5 * ip_47[k]
                 - f_6 * ip_53[k];

        g_8[k] = -f_5 * ip_3[k]
                 + f_6 * ip_24[k]
                 + f_5 * ip_45[k]
                 - f_6 * ip_51[k];

        g_9[k] = -f_7 * ip_13[k]
                 - f_8 * ip_34[k]
                 + f_9 * ip_40[k]
                 + f_10 * ip_67[k]
                 - f_11 * ip_73[k];
    }

#pragma omp simd aligned(ip_12, ip_14, ip_33, ip_35, ip_39, ip_41, ip_66, ip_68, ip_72, \
                         ip_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_7 * ip_14[k]
                  - f_8 * ip_35[k]
                  + f_9 * ip_41[k]
                  + f_10 * ip_68[k]
                  - f_11 * ip_74[k];

        g_11[k] = -f_7 * ip_12[k]
                  - f_8 * ip_33[k]
                  + f_9 * ip_39[k]
                  + f_10 * ip_66[k]
                  - f_11 * ip_72[k];
    }

#pragma omp simd aligned(ip_4, ip_5, ip_19, ip_20, ip_25, ip_26, ip_46, ip_47, ip_52, ip_53, \
                         ip_58, ip_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_12 * ip_4[k]
                  + f_13 * ip_19[k]
                  - f_14 * ip_25[k]
                  + f_12 * ip_46[k]
                  - f_14 * ip_52[k]
                  + f_14 * ip_58[k];

        g_13[k] = f_12 * ip_5[k]
                  + f_13 * ip_20[k]
                  - f_14 * ip_26[k]
                  + f_12 * ip_47[k]
                  - f_14 * ip_53[k]
                  + f_14 * ip_59[k];
    }

#pragma omp simd aligned(ip_3, ip_13, ip_18, ip_24, ip_34, ip_40, ip_45, ip_51, ip_57, ip_67, \
                         ip_73, ip_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_12 * ip_3[k]
                  + f_13 * ip_18[k]
                  - f_14 * ip_24[k]
                  + f_12 * ip_45[k]
                  - f_14 * ip_51[k]
                  + f_14 * ip_57[k];

        g_15[k] = f_15 * ip_13[k]
                  + f_16 * ip_34[k]
                  - f_17 * ip_40[k]
                  + f_15 * ip_67[k]
                  - f_17 * ip_73[k]
                  + f_18 * ip_79[k];
    }

#pragma omp simd aligned(ip_12, ip_14, ip_33, ip_35, ip_39, ip_41, ip_66, ip_68, ip_72, ip_74, \
                         ip_78, ip_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_15 * ip_14[k]
                  + f_16 * ip_35[k]
                  - f_17 * ip_41[k]
                  + f_15 * ip_68[k]
                  - f_17 * ip_74[k]
                  + f_18 * ip_80[k];

        g_17[k] = f_15 * ip_12[k]
                  + f_16 * ip_33[k]
                  - f_17 * ip_39[k]
                  + f_15 * ip_66[k]
                  - f_17 * ip_72[k]
                  + f_18 * ip_78[k];
    }

#pragma omp simd aligned(ip_1, ip_10, ip_16, ip_31, ip_37, ip_43, ip_64, ip_70, ip_76, \
                         ip_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -0.3125 * ip_1[k]
                  - 0.9375 * ip_10[k]
                  + 5.625 * ip_16[k]
                  - 0.9375 * ip_31[k]
                  + 11.25 * ip_37[k]
                  - 7.5 * ip_43[k]
                  - 0.3125 * ip_64[k]
                  + 5.625 * ip_70[k]
                  - 7.5 * ip_76[k]
                  + ip_82[k];
    }

#pragma omp simd aligned(ip_2, ip_11, ip_17, ip_32, ip_38, ip_44, ip_65, ip_71, ip_77, \
                         ip_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -0.3125 * ip_2[k]
                  - 0.9375 * ip_11[k]
                  + 5.625 * ip_17[k]
                  - 0.9375 * ip_32[k]
                  + 11.25 * ip_38[k]
                  - 7.5 * ip_44[k]
                  - 0.3125 * ip_65[k]
                  + 5.625 * ip_71[k]
                  - 7.5 * ip_77[k]
                  + ip_83[k];
    }

#pragma omp simd aligned(ip_0, ip_9, ip_15, ip_30, ip_36, ip_42, ip_63, ip_69, ip_75, \
                         ip_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -0.3125 * ip_0[k]
                  - 0.9375 * ip_9[k]
                  + 5.625 * ip_15[k]
                  - 0.9375 * ip_30[k]
                  + 11.25 * ip_36[k]
                  - 7.5 * ip_42[k]
                  - 0.3125 * ip_63[k]
                  + 5.625 * ip_69[k]
                  - 7.5 * ip_75[k]
                  + ip_81[k];
    }

#pragma omp simd aligned(ip_7, ip_8, ip_22, ip_23, ip_28, ip_29, ip_49, ip_50, ip_55, ip_56, \
                         ip_61, ip_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_15 * ip_7[k]
                  + f_16 * ip_22[k]
                  - f_17 * ip_28[k]
                  + f_15 * ip_49[k]
                  - f_17 * ip_55[k]
                  + f_18 * ip_61[k];

        g_22[k] = f_15 * ip_8[k]
                  + f_16 * ip_23[k]
                  - f_17 * ip_29[k]
                  + f_15 * ip_50[k]
                  - f_17 * ip_56[k]
                  + f_18 * ip_62[k];
    }

#pragma omp simd aligned(ip_1, ip_6, ip_10, ip_16, ip_21, ip_27, ip_31, ip_43, ip_48, ip_54, \
                         ip_60, ip_64, ip_70, ip_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_15 * ip_6[k]
                  + f_16 * ip_21[k]
                  - f_17 * ip_27[k]
                  + f_15 * ip_48[k]
                  - f_17 * ip_54[k]
                  + f_18 * ip_60[k];

        g_24[k] = f_19 * ip_1[k]
                  + f_19 * ip_10[k]
                  - f_11 * ip_16[k]
                  - f_19 * ip_31[k]
                  + f_11 * ip_43[k]
                  - f_19 * ip_64[k]
                  + f_11 * ip_70[k]
                  - f_11 * ip_76[k];
    }

#pragma omp simd aligned(ip_2, ip_11, ip_17, ip_32, ip_44, ip_65, ip_71, \
                         ip_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_19 * ip_2[k]
                  + f_19 * ip_11[k]
                  - f_11 * ip_17[k]
                  - f_19 * ip_32[k]
                  + f_11 * ip_44[k]
                  - f_19 * ip_65[k]
                  + f_11 * ip_71[k]
                  - f_11 * ip_77[k];
    }

#pragma omp simd aligned(ip_0, ip_7, ip_9, ip_15, ip_22, ip_28, ip_30, ip_42, ip_49, ip_55, \
                         ip_63, ip_69, ip_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = f_19 * ip_0[k]
                  + f_19 * ip_9[k]
                  - f_11 * ip_15[k]
                  - f_19 * ip_30[k]
                  + f_11 * ip_42[k]
                  - f_19 * ip_63[k]
                  + f_11 * ip_69[k]
                  - f_11 * ip_75[k];

        g_27[k] = -f_10 * ip_7[k]
                  + f_8 * ip_22[k]
                  + f_11 * ip_28[k]
                  + f_7 * ip_49[k]
                  - f_9 * ip_55[k];
    }

#pragma omp simd aligned(ip_6, ip_8, ip_21, ip_23, ip_27, ip_29, ip_48, ip_50, ip_54, \
                         ip_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_10 * ip_8[k]
                  + f_8 * ip_23[k]
                  + f_11 * ip_29[k]
                  + f_7 * ip_50[k]
                  - f_9 * ip_56[k];

        g_29[k] = -f_10 * ip_6[k]
                  + f_8 * ip_21[k]
                  + f_11 * ip_27[k]
                  + f_7 * ip_48[k]
                  - f_9 * ip_54[k];
    }

#pragma omp simd aligned(ip_1, ip_2, ip_10, ip_11, ip_16, ip_17, ip_31, ip_32, ip_37, ip_38, \
                         ip_64, ip_65, ip_70, ip_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_20 * ip_1[k]
                  + f_21 * ip_10[k]
                  + f_22 * ip_16[k]
                  + f_21 * ip_31[k]
                  - f_23 * ip_37[k]
                  - f_20 * ip_64[k]
                  + f_22 * ip_70[k];

        g_31[k] = -f_20 * ip_2[k]
                  + f_21 * ip_11[k]
                  + f_22 * ip_17[k]
                  + f_21 * ip_32[k]
                  - f_23 * ip_38[k]
                  - f_20 * ip_65[k]
                  + f_22 * ip_71[k];
    }

#pragma omp simd aligned(ip_0, ip_7, ip_8, ip_9, ip_15, ip_22, ip_23, ip_30, ip_36, ip_49, \
                         ip_50, ip_63, ip_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_20 * ip_0[k]
                  + f_21 * ip_9[k]
                  + f_22 * ip_15[k]
                  + f_21 * ip_30[k]
                  - f_23 * ip_36[k]
                  - f_20 * ip_63[k]
                  + f_22 * ip_69[k];

        g_33[k] = f_4 * ip_7[k]
                  - f_3 * ip_22[k]
                  + f_2 * ip_49[k];

        g_34[k] = f_4 * ip_8[k]
                  - f_3 * ip_23[k]
                  + f_2 * ip_50[k];
    }

#pragma omp simd aligned(ip_1, ip_2, ip_6, ip_10, ip_11, ip_21, ip_31, ip_32, ip_48, ip_64, \
                         ip_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_4 * ip_6[k]
                  - f_3 * ip_21[k]
                  + f_2 * ip_48[k];

        g_36[k] = f_24 * ip_1[k]
                  - f_25 * ip_10[k]
                  + f_25 * ip_31[k]
                  - f_24 * ip_64[k];

        g_37[k] = f_24 * ip_2[k]
                  - f_25 * ip_11[k]
                  + f_25 * ip_32[k]
                  - f_24 * ip_65[k];
    }

#pragma omp simd aligned(ip_0, ip_9, ip_30, ip_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_24 * ip_0[k]
                  - f_25 * ip_9[k]
                  + f_25 * ip_30[k]
                  - f_24 * ip_63[k];
    }
}

}  // namespace simdtrf
