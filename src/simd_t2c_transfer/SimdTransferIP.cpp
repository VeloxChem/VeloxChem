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


#include "SimdTransferIP.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_ip(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t is, const size_t ks, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_21 = buffer.data(is + 21);
    const auto *is_22 = buffer.data(is + 22);
    const auto *is_23 = buffer.data(is + 23);
    const auto *is_24 = buffer.data(is + 24);
    const auto *is_25 = buffer.data(is + 25);
    const auto *is_26 = buffer.data(is + 26);
    const auto *is_27 = buffer.data(is + 27);

    const auto *ks_0 = buffer.data(ks + 0);
    const auto *ks_1 = buffer.data(ks + 1);
    const auto *ks_2 = buffer.data(ks + 2);
    const auto *ks_3 = buffer.data(ks + 3);
    const auto *ks_4 = buffer.data(ks + 4);
    const auto *ks_5 = buffer.data(ks + 5);
    const auto *ks_6 = buffer.data(ks + 6);
    const auto *ks_7 = buffer.data(ks + 7);
    const auto *ks_8 = buffer.data(ks + 8);
    const auto *ks_9 = buffer.data(ks + 9);
    const auto *ks_10 = buffer.data(ks + 10);
    const auto *ks_11 = buffer.data(ks + 11);
    const auto *ks_12 = buffer.data(ks + 12);
    const auto *ks_13 = buffer.data(ks + 13);
    const auto *ks_14 = buffer.data(ks + 14);
    const auto *ks_15 = buffer.data(ks + 15);
    const auto *ks_16 = buffer.data(ks + 16);
    const auto *ks_17 = buffer.data(ks + 17);
    const auto *ks_18 = buffer.data(ks + 18);
    const auto *ks_19 = buffer.data(ks + 19);
    const auto *ks_20 = buffer.data(ks + 20);
    const auto *ks_21 = buffer.data(ks + 21);
    const auto *ks_22 = buffer.data(ks + 22);
    const auto *ks_23 = buffer.data(ks + 23);
    const auto *ks_24 = buffer.data(ks + 24);
    const auto *ks_25 = buffer.data(ks + 25);
    const auto *ks_26 = buffer.data(ks + 26);
    const auto *ks_27 = buffer.data(ks + 27);
    const auto *ks_28 = buffer.data(ks + 28);
    const auto *ks_29 = buffer.data(ks + 29);
    const auto *ks_30 = buffer.data(ks + 30);
    const auto *ks_31 = buffer.data(ks + 31);
    const auto *ks_32 = buffer.data(ks + 32);
    const auto *ks_33 = buffer.data(ks + 33);
    const auto *ks_34 = buffer.data(ks + 34);
    const auto *ks_35 = buffer.data(ks + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, is_0, is_1, ks_0, \
                         ks_1, ks_2, ks_3, ks_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * is_0[k]
                 + ks_0[k];

        t_1[k] = ab_y[k] * is_0[k]
                 + ks_1[k];

        t_2[k] = ab_z[k] * is_0[k]
                 + ks_2[k];

        t_3[k] = ab_x[k] * is_1[k]
                 + ks_1[k];

        t_4[k] = ab_y[k] * is_1[k]
                 + ks_3[k];

        t_5[k] = ab_z[k] * is_1[k]
                 + ks_4[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, is_2, is_3, ks_2, ks_3, \
                         ks_4, ks_5, ks_6 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_6[k] = ab_x[k] * is_2[k]
                 + ks_2[k];

        t_7[k] = ab_y[k] * is_2[k]
                 + ks_4[k];

        t_8[k] = ab_z[k] * is_2[k]
                 + ks_5[k];

        t_9[k] = ab_x[k] * is_3[k]
                 + ks_3[k];

        t_10[k] = ab_y[k] * is_3[k]
                  + ks_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, is_3, is_4, \
                         is_5, ks_4, ks_5, ks_7, ks_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_11[k] = ab_z[k] * is_3[k]
                  + ks_7[k];

        t_12[k] = ab_x[k] * is_4[k]
                  + ks_4[k];

        t_13[k] = ab_y[k] * is_4[k]
                  + ks_7[k];

        t_14[k] = ab_z[k] * is_4[k]
                  + ks_8[k];

        t_15[k] = ab_x[k] * is_5[k]
                  + ks_5[k];

        t_16[k] = ab_y[k] * is_5[k]
                  + ks_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, is_5, is_6, is_7, \
                         ks_6, ks_7, ks_9, ks_10, ks_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_17[k] = ab_z[k] * is_5[k]
                  + ks_9[k];

        t_18[k] = ab_x[k] * is_6[k]
                  + ks_6[k];

        t_19[k] = ab_y[k] * is_6[k]
                  + ks_10[k];

        t_20[k] = ab_z[k] * is_6[k]
                  + ks_11[k];

        t_21[k] = ab_x[k] * is_7[k]
                  + ks_7[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, is_7, is_8, ks_8, \
                         ks_11, ks_12, ks_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_22[k] = ab_y[k] * is_7[k]
                  + ks_11[k];

        t_23[k] = ab_z[k] * is_7[k]
                  + ks_12[k];

        t_24[k] = ab_x[k] * is_8[k]
                  + ks_8[k];

        t_25[k] = ab_y[k] * is_8[k]
                  + ks_12[k];

        t_26[k] = ab_z[k] * is_8[k]
                  + ks_13[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, is_9, is_10, ks_9, \
                         ks_10, ks_13, ks_14, ks_15 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_27[k] = ab_x[k] * is_9[k]
                  + ks_9[k];

        t_28[k] = ab_y[k] * is_9[k]
                  + ks_13[k];

        t_29[k] = ab_z[k] * is_9[k]
                  + ks_14[k];

        t_30[k] = ab_x[k] * is_10[k]
                  + ks_10[k];

        t_31[k] = ab_y[k] * is_10[k]
                  + ks_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, ab_x, ab_y, ab_z, is_10, is_11, \
                         is_12, ks_11, ks_12, ks_16, ks_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_32[k] = ab_z[k] * is_10[k]
                  + ks_16[k];

        t_33[k] = ab_x[k] * is_11[k]
                  + ks_11[k];

        t_34[k] = ab_y[k] * is_11[k]
                  + ks_16[k];

        t_35[k] = ab_z[k] * is_11[k]
                  + ks_17[k];

        t_36[k] = ab_x[k] * is_12[k]
                  + ks_12[k];

        t_37[k] = ab_y[k] * is_12[k]
                  + ks_17[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, is_12, is_13, \
                         is_14, ks_13, ks_14, ks_18, ks_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_38[k] = ab_z[k] * is_12[k]
                  + ks_18[k];

        t_39[k] = ab_x[k] * is_13[k]
                  + ks_13[k];

        t_40[k] = ab_y[k] * is_13[k]
                  + ks_18[k];

        t_41[k] = ab_z[k] * is_13[k]
                  + ks_19[k];

        t_42[k] = ab_x[k] * is_14[k]
                  + ks_14[k];

        t_43[k] = ab_y[k] * is_14[k]
                  + ks_19[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, is_14, is_15, is_16, \
                         ks_15, ks_16, ks_20, ks_21, ks_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_44[k] = ab_z[k] * is_14[k]
                  + ks_20[k];

        t_45[k] = ab_x[k] * is_15[k]
                  + ks_15[k];

        t_46[k] = ab_y[k] * is_15[k]
                  + ks_21[k];

        t_47[k] = ab_z[k] * is_15[k]
                  + ks_22[k];

        t_48[k] = ab_x[k] * is_16[k]
                  + ks_16[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, is_16, is_17, ks_17, \
                         ks_22, ks_23, ks_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_49[k] = ab_y[k] * is_16[k]
                  + ks_22[k];

        t_50[k] = ab_z[k] * is_16[k]
                  + ks_23[k];

        t_51[k] = ab_x[k] * is_17[k]
                  + ks_17[k];

        t_52[k] = ab_y[k] * is_17[k]
                  + ks_23[k];

        t_53[k] = ab_z[k] * is_17[k]
                  + ks_24[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, is_18, is_19, \
                         ks_18, ks_19, ks_24, ks_25, ks_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_54[k] = ab_x[k] * is_18[k]
                  + ks_18[k];

        t_55[k] = ab_y[k] * is_18[k]
                  + ks_24[k];

        t_56[k] = ab_z[k] * is_18[k]
                  + ks_25[k];

        t_57[k] = ab_x[k] * is_19[k]
                  + ks_19[k];

        t_58[k] = ab_y[k] * is_19[k]
                  + ks_25[k];

        t_59[k] = ab_z[k] * is_19[k]
                  + ks_26[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, ab_y, ab_z, is_20, is_21, ks_20, \
                         ks_21, ks_26, ks_27, ks_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = ab_x[k] * is_20[k]
                  + ks_20[k];

        t_61[k] = ab_y[k] * is_20[k]
                  + ks_26[k];

        t_62[k] = ab_z[k] * is_20[k]
                  + ks_27[k];

        t_63[k] = ab_x[k] * is_21[k]
                  + ks_21[k];

        t_64[k] = ab_y[k] * is_21[k]
                  + ks_28[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, t_70, ab_x, ab_y, ab_z, is_21, is_22, \
                         is_23, ks_22, ks_23, ks_29, ks_30 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = ab_z[k] * is_21[k]
                  + ks_29[k];

        t_66[k] = ab_x[k] * is_22[k]
                  + ks_22[k];

        t_67[k] = ab_y[k] * is_22[k]
                  + ks_29[k];

        t_68[k] = ab_z[k] * is_22[k]
                  + ks_30[k];

        t_69[k] = ab_x[k] * is_23[k]
                  + ks_23[k];

        t_70[k] = ab_y[k] * is_23[k]
                  + ks_30[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, ab_z, is_23, is_24, \
                         is_25, ks_24, ks_25, ks_31, ks_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_71[k] = ab_z[k] * is_23[k]
                  + ks_31[k];

        t_72[k] = ab_x[k] * is_24[k]
                  + ks_24[k];

        t_73[k] = ab_y[k] * is_24[k]
                  + ks_31[k];

        t_74[k] = ab_z[k] * is_24[k]
                  + ks_32[k];

        t_75[k] = ab_x[k] * is_25[k]
                  + ks_25[k];

        t_76[k] = ab_y[k] * is_25[k]
                  + ks_32[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, t_82, ab_x, ab_y, ab_z, is_25, is_26, \
                         is_27, ks_26, ks_27, ks_33, ks_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_77[k] = ab_z[k] * is_25[k]
                  + ks_33[k];

        t_78[k] = ab_x[k] * is_26[k]
                  + ks_26[k];

        t_79[k] = ab_y[k] * is_26[k]
                  + ks_33[k];

        t_80[k] = ab_z[k] * is_26[k]
                  + ks_34[k];

        t_81[k] = ab_x[k] * is_27[k]
                  + ks_27[k];

        t_82[k] = ab_y[k] * is_27[k]
                  + ks_34[k];
    }

#pragma omp simd aligned(t_83, ab_z, is_27, ks_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_83[k] = ab_z[k] * is_27[k]
                  + ks_35[k];
    }
}

auto
compute_hrr_ip_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t is, const size_t ks,
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_21 = buffer.data(is + 21);
    const auto *is_22 = buffer.data(is + 22);
    const auto *is_23 = buffer.data(is + 23);
    const auto *is_24 = buffer.data(is + 24);
    const auto *is_25 = buffer.data(is + 25);
    const auto *is_26 = buffer.data(is + 26);
    const auto *is_27 = buffer.data(is + 27);

    const auto *ks_0 = buffer.data(ks + 0);
    const auto *ks_1 = buffer.data(ks + 1);
    const auto *ks_2 = buffer.data(ks + 2);
    const auto *ks_3 = buffer.data(ks + 3);
    const auto *ks_4 = buffer.data(ks + 4);
    const auto *ks_5 = buffer.data(ks + 5);
    const auto *ks_6 = buffer.data(ks + 6);
    const auto *ks_7 = buffer.data(ks + 7);
    const auto *ks_8 = buffer.data(ks + 8);
    const auto *ks_9 = buffer.data(ks + 9);
    const auto *ks_10 = buffer.data(ks + 10);
    const auto *ks_11 = buffer.data(ks + 11);
    const auto *ks_12 = buffer.data(ks + 12);
    const auto *ks_13 = buffer.data(ks + 13);
    const auto *ks_14 = buffer.data(ks + 14);
    const auto *ks_15 = buffer.data(ks + 15);
    const auto *ks_16 = buffer.data(ks + 16);
    const auto *ks_17 = buffer.data(ks + 17);
    const auto *ks_18 = buffer.data(ks + 18);
    const auto *ks_19 = buffer.data(ks + 19);
    const auto *ks_20 = buffer.data(ks + 20);
    const auto *ks_21 = buffer.data(ks + 21);
    const auto *ks_22 = buffer.data(ks + 22);
    const auto *ks_23 = buffer.data(ks + 23);
    const auto *ks_24 = buffer.data(ks + 24);
    const auto *ks_25 = buffer.data(ks + 25);
    const auto *ks_26 = buffer.data(ks + 26);
    const auto *ks_27 = buffer.data(ks + 27);
    const auto *ks_28 = buffer.data(ks + 28);
    const auto *ks_29 = buffer.data(ks + 29);
    const auto *ks_30 = buffer.data(ks + 30);
    const auto *ks_31 = buffer.data(ks + 31);
    const auto *ks_32 = buffer.data(ks + 32);
    const auto *ks_33 = buffer.data(ks + 33);
    const auto *ks_34 = buffer.data(ks + 34);
    const auto *ks_35 = buffer.data(ks + 35);

#pragma omp simd aligned(ab_y, ab_z, is_1, is_6, is_15, ks_3, ks_4, ks_10, ks_11, ks_21, \
                         ks_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * ab_y[k] * is_1[k]
                 - f_1 * ab_y[k] * is_6[k]
                 + f_0 * ab_y[k] * is_15[k]
                 + f_0 * ks_3[k]
                 - f_1 * ks_10[k]
                 + f_0 * ks_21[k];

        g_1[k] = f_0 * ab_z[k] * is_1[k]
                 - f_1 * ab_z[k] * is_6[k]
                 + f_0 * ab_z[k] * is_15[k]
                 + f_0 * ks_4[k]
                 - f_1 * ks_11[k]
                 + f_0 * ks_22[k];
    }

#pragma omp simd aligned(ab_x, ab_y, is_1, is_4, is_6, is_11, is_15, is_22, ks_1, ks_6, ks_7, \
                         ks_15, ks_16, ks_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = f_0 * ab_x[k] * is_1[k]
                 - f_1 * ab_x[k] * is_6[k]
                 + f_0 * ab_x[k] * is_15[k]
                 + f_0 * ks_1[k]
                 - f_1 * ks_6[k]
                 + f_0 * ks_15[k];

        g_3[k] = f_2 * ab_y[k] * is_4[k]
                 - f_3 * ab_y[k] * is_11[k]
                 + f_4 * ab_y[k] * is_22[k]
                 + f_2 * ks_7[k]
                 - f_3 * ks_16[k]
                 + f_4 * ks_29[k];
    }

#pragma omp simd aligned(ab_x, ab_z, is_4, is_11, is_22, ks_4, ks_8, ks_11, ks_17, ks_22, \
                         ks_30 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_2 * ab_z[k] * is_4[k]
                 - f_3 * ab_z[k] * is_11[k]
                 + f_4 * ab_z[k] * is_22[k]
                 + f_2 * ks_8[k]
                 - f_3 * ks_17[k]
                 + f_4 * ks_30[k];

        g_5[k] = f_2 * ab_x[k] * is_4[k]
                 - f_3 * ab_x[k] * is_11[k]
                 + f_4 * ab_x[k] * is_22[k]
                 + f_2 * ks_4[k]
                 - f_3 * ks_11[k]
                 + f_4 * ks_22[k];
    }

#pragma omp simd aligned(ab_y, ab_z, is_1, is_8, is_15, is_17, ks_3, ks_4, ks_12, ks_13, \
                         ks_21, ks_22, ks_23, ks_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_5 * ab_y[k] * is_1[k]
                 + f_6 * ab_y[k] * is_8[k]
                 + f_5 * ab_y[k] * is_15[k]
                 - f_6 * ab_y[k] * is_17[k]
                 - f_5 * ks_3[k]
                 + f_6 * ks_12[k]
                 + f_5 * ks_21[k]
                 - f_6 * ks_23[k];

        g_7[k] = -f_5 * ab_z[k] * is_1[k]
                 + f_6 * ab_z[k] * is_8[k]
                 + f_5 * ab_z[k] * is_15[k]
                 - f_6 * ab_z[k] * is_17[k]
                 - f_5 * ks_4[k]
                 + f_6 * ks_13[k]
                 + f_5 * ks_22[k]
                 - f_6 * ks_24[k];
    }

#pragma omp simd aligned(ab_x, is_1, is_8, is_15, is_17, ks_1, ks_8, ks_15, \
                         ks_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_5 * ab_x[k] * is_1[k]
                 + f_6 * ab_x[k] * is_8[k]
                 + f_5 * ab_x[k] * is_15[k]
                 - f_6 * ab_x[k] * is_17[k]
                 - f_5 * ks_1[k]
                 + f_6 * ks_8[k]
                 + f_5 * ks_15[k]
                 - f_6 * ks_17[k];
    }

#pragma omp simd aligned(ab_y, is_4, is_11, is_13, is_22, is_24, ks_7, ks_16, ks_18, ks_29, \
                         ks_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_7 * ab_y[k] * is_4[k]
                 - f_8 * ab_y[k] * is_11[k]
                 + f_9 * ab_y[k] * is_13[k]
                 + f_10 * ab_y[k] * is_22[k]
                 - f_11 * ab_y[k] * is_24[k]
                 - f_7 * ks_7[k]
                 - f_8 * ks_16[k]
                 + f_9 * ks_18[k]
                 + f_10 * ks_29[k]
                 - f_11 * ks_31[k];
    }

#pragma omp simd aligned(ab_z, is_4, is_11, is_13, is_22, is_24, ks_8, ks_17, ks_19, ks_30, \
                         ks_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_7 * ab_z[k] * is_4[k]
                  - f_8 * ab_z[k] * is_11[k]
                  + f_9 * ab_z[k] * is_13[k]
                  + f_10 * ab_z[k] * is_22[k]
                  - f_11 * ab_z[k] * is_24[k]
                  - f_7 * ks_8[k]
                  - f_8 * ks_17[k]
                  + f_9 * ks_19[k]
                  + f_10 * ks_30[k]
                  - f_11 * ks_32[k];
    }

#pragma omp simd aligned(ab_x, is_4, is_11, is_13, is_22, is_24, ks_4, ks_11, ks_13, ks_22, \
                         ks_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_7 * ab_x[k] * is_4[k]
                  - f_8 * ab_x[k] * is_11[k]
                  + f_9 * ab_x[k] * is_13[k]
                  + f_10 * ab_x[k] * is_22[k]
                  - f_11 * ab_x[k] * is_24[k]
                  - f_7 * ks_4[k]
                  - f_8 * ks_11[k]
                  + f_9 * ks_13[k]
                  + f_10 * ks_22[k]
                  - f_11 * ks_24[k];
    }

#pragma omp simd aligned(ab_y, is_1, is_6, is_8, is_15, is_17, is_19, ks_3, ks_10, ks_12, \
                         ks_21, ks_23, ks_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_12 * ab_y[k] * is_1[k]
                  + f_13 * ab_y[k] * is_6[k]
                  - f_14 * ab_y[k] * is_8[k]
                  + f_12 * ab_y[k] * is_15[k]
                  - f_14 * ab_y[k] * is_17[k]
                  + f_14 * ab_y[k] * is_19[k]
                  + f_12 * ks_3[k]
                  + f_13 * ks_10[k]
                  - f_14 * ks_12[k]
                  + f_12 * ks_21[k]
                  - f_14 * ks_23[k]
                  + f_14 * ks_25[k];
    }

#pragma omp simd aligned(ab_z, is_1, is_6, is_8, is_15, is_17, is_19, ks_4, ks_11, ks_13, \
                         ks_22, ks_24, ks_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_12 * ab_z[k] * is_1[k]
                  + f_13 * ab_z[k] * is_6[k]
                  - f_14 * ab_z[k] * is_8[k]
                  + f_12 * ab_z[k] * is_15[k]
                  - f_14 * ab_z[k] * is_17[k]
                  + f_14 * ab_z[k] * is_19[k]
                  + f_12 * ks_4[k]
                  + f_13 * ks_11[k]
                  - f_14 * ks_13[k]
                  + f_12 * ks_22[k]
                  - f_14 * ks_24[k]
                  + f_14 * ks_26[k];
    }

#pragma omp simd aligned(ab_x, is_1, is_6, is_8, is_15, is_17, is_19, ks_1, ks_6, ks_8, ks_15, \
                         ks_17, ks_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_12 * ab_x[k] * is_1[k]
                  + f_13 * ab_x[k] * is_6[k]
                  - f_14 * ab_x[k] * is_8[k]
                  + f_12 * ab_x[k] * is_15[k]
                  - f_14 * ab_x[k] * is_17[k]
                  + f_14 * ab_x[k] * is_19[k]
                  + f_12 * ks_1[k]
                  + f_13 * ks_6[k]
                  - f_14 * ks_8[k]
                  + f_12 * ks_15[k]
                  - f_14 * ks_17[k]
                  + f_14 * ks_19[k];
    }

#pragma omp simd aligned(ab_y, is_4, is_11, is_13, is_22, is_24, is_26, ks_7, ks_16, ks_18, \
                         ks_29, ks_31, ks_33 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_15 * ab_y[k] * is_4[k]
                  + f_16 * ab_y[k] * is_11[k]
                  - f_17 * ab_y[k] * is_13[k]
                  + f_15 * ab_y[k] * is_22[k]
                  - f_17 * ab_y[k] * is_24[k]
                  + f_18 * ab_y[k] * is_26[k]
                  + f_15 * ks_7[k]
                  + f_16 * ks_16[k]
                  - f_17 * ks_18[k]
                  + f_15 * ks_29[k]
                  - f_17 * ks_31[k]
                  + f_18 * ks_33[k];
    }

#pragma omp simd aligned(ab_z, is_4, is_11, is_13, is_22, is_24, is_26, ks_8, ks_17, ks_19, \
                         ks_30, ks_32, ks_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_15 * ab_z[k] * is_4[k]
                  + f_16 * ab_z[k] * is_11[k]
                  - f_17 * ab_z[k] * is_13[k]
                  + f_15 * ab_z[k] * is_22[k]
                  - f_17 * ab_z[k] * is_24[k]
                  + f_18 * ab_z[k] * is_26[k]
                  + f_15 * ks_8[k]
                  + f_16 * ks_17[k]
                  - f_17 * ks_19[k]
                  + f_15 * ks_30[k]
                  - f_17 * ks_32[k]
                  + f_18 * ks_34[k];
    }

#pragma omp simd aligned(ab_x, is_4, is_11, is_13, is_22, is_24, is_26, ks_4, ks_11, ks_13, \
                         ks_22, ks_24, ks_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_15 * ab_x[k] * is_4[k]
                  + f_16 * ab_x[k] * is_11[k]
                  - f_17 * ab_x[k] * is_13[k]
                  + f_15 * ab_x[k] * is_22[k]
                  - f_17 * ab_x[k] * is_24[k]
                  + f_18 * ab_x[k] * is_26[k]
                  + f_15 * ks_4[k]
                  + f_16 * ks_11[k]
                  - f_17 * ks_13[k]
                  + f_15 * ks_22[k]
                  - f_17 * ks_24[k]
                  + f_18 * ks_26[k];
    }

#pragma omp simd aligned(ab_y, is_0, is_3, is_5, is_10, is_12, is_14, is_21, is_23, is_25, \
                         is_27, ks_1, ks_6, ks_8, ks_15, ks_17, ks_19, ks_28, ks_30, ks_32, \
                         ks_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -0.3125 * ab_y[k] * is_0[k]
                  - 0.9375 * ab_y[k] * is_3[k]
                  + 5.625 * ab_y[k] * is_5[k]
                  - 0.9375 * ab_y[k] * is_10[k]
                  + 11.25 * ab_y[k] * is_12[k]
                  - 7.5 * ab_y[k] * is_14[k]
                  - 0.3125 * ab_y[k] * is_21[k]
                  + 5.625 * ab_y[k] * is_23[k]
                  - 7.5 * ab_y[k] * is_25[k]
                  + ab_y[k] * is_27[k]
                  - 0.3125 * ks_1[k]
                  - 0.9375 * ks_6[k]
                  + 5.625 * ks_8[k]
                  - 0.9375 * ks_15[k]
                  + 11.25 * ks_17[k]
                  - 7.5 * ks_19[k]
                  - 0.3125 * ks_28[k]
                  + 5.625 * ks_30[k]
                  - 7.5 * ks_32[k]
                  + ks_34[k];
    }

#pragma omp simd aligned(ab_z, is_0, is_3, is_5, is_10, is_12, is_14, is_21, is_23, is_25, \
                         is_27, ks_2, ks_7, ks_9, ks_16, ks_18, ks_20, ks_29, ks_31, ks_33, \
                         ks_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -0.3125 * ab_z[k] * is_0[k]
                  - 0.9375 * ab_z[k] * is_3[k]
                  + 5.625 * ab_z[k] * is_5[k]
                  - 0.9375 * ab_z[k] * is_10[k]
                  + 11.25 * ab_z[k] * is_12[k]
                  - 7.5 * ab_z[k] * is_14[k]
                  - 0.3125 * ab_z[k] * is_21[k]
                  + 5.625 * ab_z[k] * is_23[k]
                  - 7.5 * ab_z[k] * is_25[k]
                  + ab_z[k] * is_27[k]
                  - 0.3125 * ks_2[k]
                  - 0.9375 * ks_7[k]
                  + 5.625 * ks_9[k]
                  - 0.9375 * ks_16[k]
                  + 11.25 * ks_18[k]
                  - 7.5 * ks_20[k]
                  - 0.3125 * ks_29[k]
                  + 5.625 * ks_31[k]
                  - 7.5 * ks_33[k]
                  + ks_35[k];
    }

#pragma omp simd aligned(ab_x, is_0, is_3, is_5, is_10, is_12, is_14, is_21, is_23, is_25, \
                         is_27, ks_0, ks_3, ks_5, ks_10, ks_12, ks_14, ks_21, ks_23, ks_25, \
                         ks_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -0.3125 * ab_x[k] * is_0[k]
                  - 0.9375 * ab_x[k] * is_3[k]
                  + 5.625 * ab_x[k] * is_5[k]
                  - 0.9375 * ab_x[k] * is_10[k]
                  + 11.25 * ab_x[k] * is_12[k]
                  - 7.5 * ab_x[k] * is_14[k]
                  - 0.3125 * ab_x[k] * is_21[k]
                  + 5.625 * ab_x[k] * is_23[k]
                  - 7.5 * ab_x[k] * is_25[k]
                  + ab_x[k] * is_27[k]
                  - 0.3125 * ks_0[k]
                  - 0.9375 * ks_3[k]
                  + 5.625 * ks_5[k]
                  - 0.9375 * ks_10[k]
                  + 11.25 * ks_12[k]
                  - 7.5 * ks_14[k]
                  - 0.3125 * ks_21[k]
                  + 5.625 * ks_23[k]
                  - 7.5 * ks_25[k]
                  + ks_27[k];
    }

#pragma omp simd aligned(ab_y, is_2, is_7, is_9, is_16, is_18, is_20, ks_4, ks_11, ks_13, \
                         ks_22, ks_24, ks_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_15 * ab_y[k] * is_2[k]
                  + f_16 * ab_y[k] * is_7[k]
                  - f_17 * ab_y[k] * is_9[k]
                  + f_15 * ab_y[k] * is_16[k]
                  - f_17 * ab_y[k] * is_18[k]
                  + f_18 * ab_y[k] * is_20[k]
                  + f_15 * ks_4[k]
                  + f_16 * ks_11[k]
                  - f_17 * ks_13[k]
                  + f_15 * ks_22[k]
                  - f_17 * ks_24[k]
                  + f_18 * ks_26[k];
    }

#pragma omp simd aligned(ab_z, is_2, is_7, is_9, is_16, is_18, is_20, ks_5, ks_12, ks_14, \
                         ks_23, ks_25, ks_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_15 * ab_z[k] * is_2[k]
                  + f_16 * ab_z[k] * is_7[k]
                  - f_17 * ab_z[k] * is_9[k]
                  + f_15 * ab_z[k] * is_16[k]
                  - f_17 * ab_z[k] * is_18[k]
                  + f_18 * ab_z[k] * is_20[k]
                  + f_15 * ks_5[k]
                  + f_16 * ks_12[k]
                  - f_17 * ks_14[k]
                  + f_15 * ks_23[k]
                  - f_17 * ks_25[k]
                  + f_18 * ks_27[k];
    }

#pragma omp simd aligned(ab_x, is_2, is_7, is_9, is_16, is_18, is_20, ks_2, ks_7, ks_9, ks_16, \
                         ks_18, ks_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_15 * ab_x[k] * is_2[k]
                  + f_16 * ab_x[k] * is_7[k]
                  - f_17 * ab_x[k] * is_9[k]
                  + f_15 * ab_x[k] * is_16[k]
                  - f_17 * ab_x[k] * is_18[k]
                  + f_18 * ab_x[k] * is_20[k]
                  + f_15 * ks_2[k]
                  + f_16 * ks_7[k]
                  - f_17 * ks_9[k]
                  + f_15 * ks_16[k]
                  - f_17 * ks_18[k]
                  + f_18 * ks_20[k];
    }

#pragma omp simd aligned(ab_y, is_0, is_3, is_5, is_10, is_14, is_21, is_23, is_25, ks_1, \
                         ks_6, ks_8, ks_15, ks_19, ks_28, ks_30, \
                         ks_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_19 * ab_y[k] * is_0[k]
                  + f_19 * ab_y[k] * is_3[k]
                  - f_11 * ab_y[k] * is_5[k]
                  - f_19 * ab_y[k] * is_10[k]
                  + f_11 * ab_y[k] * is_14[k]
                  - f_19 * ab_y[k] * is_21[k]
                  + f_11 * ab_y[k] * is_23[k]
                  - f_11 * ab_y[k] * is_25[k]
                  + f_19 * ks_1[k]
                  + f_19 * ks_6[k]
                  - f_11 * ks_8[k]
                  - f_19 * ks_15[k]
                  + f_11 * ks_19[k]
                  - f_19 * ks_28[k]
                  + f_11 * ks_30[k]
                  - f_11 * ks_32[k];
    }

#pragma omp simd aligned(ab_z, is_0, is_3, is_5, is_10, is_14, is_21, is_23, is_25, ks_2, \
                         ks_7, ks_9, ks_16, ks_20, ks_29, ks_31, \
                         ks_33 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_19 * ab_z[k] * is_0[k]
                  + f_19 * ab_z[k] * is_3[k]
                  - f_11 * ab_z[k] * is_5[k]
                  - f_19 * ab_z[k] * is_10[k]
                  + f_11 * ab_z[k] * is_14[k]
                  - f_19 * ab_z[k] * is_21[k]
                  + f_11 * ab_z[k] * is_23[k]
                  - f_11 * ab_z[k] * is_25[k]
                  + f_19 * ks_2[k]
                  + f_19 * ks_7[k]
                  - f_11 * ks_9[k]
                  - f_19 * ks_16[k]
                  + f_11 * ks_20[k]
                  - f_19 * ks_29[k]
                  + f_11 * ks_31[k]
                  - f_11 * ks_33[k];
    }

#pragma omp simd aligned(ab_x, is_0, is_3, is_5, is_10, is_14, is_21, is_23, is_25, ks_0, \
                         ks_3, ks_5, ks_10, ks_14, ks_21, ks_23, \
                         ks_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = f_19 * ab_x[k] * is_0[k]
                  + f_19 * ab_x[k] * is_3[k]
                  - f_11 * ab_x[k] * is_5[k]
                  - f_19 * ab_x[k] * is_10[k]
                  + f_11 * ab_x[k] * is_14[k]
                  - f_19 * ab_x[k] * is_21[k]
                  + f_11 * ab_x[k] * is_23[k]
                  - f_11 * ab_x[k] * is_25[k]
                  + f_19 * ks_0[k]
                  + f_19 * ks_3[k]
                  - f_11 * ks_5[k]
                  - f_19 * ks_10[k]
                  + f_11 * ks_14[k]
                  - f_19 * ks_21[k]
                  + f_11 * ks_23[k]
                  - f_11 * ks_25[k];
    }

#pragma omp simd aligned(ab_y, is_2, is_7, is_9, is_16, is_18, ks_4, ks_11, ks_13, ks_22, \
                         ks_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_10 * ab_y[k] * is_2[k]
                  + f_8 * ab_y[k] * is_7[k]
                  + f_11 * ab_y[k] * is_9[k]
                  + f_7 * ab_y[k] * is_16[k]
                  - f_9 * ab_y[k] * is_18[k]
                  - f_10 * ks_4[k]
                  + f_8 * ks_11[k]
                  + f_11 * ks_13[k]
                  + f_7 * ks_22[k]
                  - f_9 * ks_24[k];
    }

#pragma omp simd aligned(ab_z, is_2, is_7, is_9, is_16, is_18, ks_5, ks_12, ks_14, ks_23, \
                         ks_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_10 * ab_z[k] * is_2[k]
                  + f_8 * ab_z[k] * is_7[k]
                  + f_11 * ab_z[k] * is_9[k]
                  + f_7 * ab_z[k] * is_16[k]
                  - f_9 * ab_z[k] * is_18[k]
                  - f_10 * ks_5[k]
                  + f_8 * ks_12[k]
                  + f_11 * ks_14[k]
                  + f_7 * ks_23[k]
                  - f_9 * ks_25[k];
    }

#pragma omp simd aligned(ab_x, is_2, is_7, is_9, is_16, is_18, ks_2, ks_7, ks_9, ks_16, \
                         ks_18 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = -f_10 * ab_x[k] * is_2[k]
                  + f_8 * ab_x[k] * is_7[k]
                  + f_11 * ab_x[k] * is_9[k]
                  + f_7 * ab_x[k] * is_16[k]
                  - f_9 * ab_x[k] * is_18[k]
                  - f_10 * ks_2[k]
                  + f_8 * ks_7[k]
                  + f_11 * ks_9[k]
                  + f_7 * ks_16[k]
                  - f_9 * ks_18[k];
    }

#pragma omp simd aligned(ab_y, is_0, is_3, is_5, is_10, is_12, is_21, is_23, ks_1, ks_6, ks_8, \
                         ks_15, ks_17, ks_28, ks_30 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_20 * ab_y[k] * is_0[k]
                  + f_21 * ab_y[k] * is_3[k]
                  + f_22 * ab_y[k] * is_5[k]
                  + f_21 * ab_y[k] * is_10[k]
                  - f_23 * ab_y[k] * is_12[k]
                  - f_20 * ab_y[k] * is_21[k]
                  + f_22 * ab_y[k] * is_23[k]
                  - f_20 * ks_1[k]
                  + f_21 * ks_6[k]
                  + f_22 * ks_8[k]
                  + f_21 * ks_15[k]
                  - f_23 * ks_17[k]
                  - f_20 * ks_28[k]
                  + f_22 * ks_30[k];
    }

#pragma omp simd aligned(ab_z, is_0, is_3, is_5, is_10, is_12, is_21, is_23, ks_2, ks_7, ks_9, \
                         ks_16, ks_18, ks_29, ks_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_20 * ab_z[k] * is_0[k]
                  + f_21 * ab_z[k] * is_3[k]
                  + f_22 * ab_z[k] * is_5[k]
                  + f_21 * ab_z[k] * is_10[k]
                  - f_23 * ab_z[k] * is_12[k]
                  - f_20 * ab_z[k] * is_21[k]
                  + f_22 * ab_z[k] * is_23[k]
                  - f_20 * ks_2[k]
                  + f_21 * ks_7[k]
                  + f_22 * ks_9[k]
                  + f_21 * ks_16[k]
                  - f_23 * ks_18[k]
                  - f_20 * ks_29[k]
                  + f_22 * ks_31[k];
    }

#pragma omp simd aligned(ab_x, is_0, is_3, is_5, is_10, is_12, is_21, is_23, ks_0, ks_3, ks_5, \
                         ks_10, ks_12, ks_21, ks_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_20 * ab_x[k] * is_0[k]
                  + f_21 * ab_x[k] * is_3[k]
                  + f_22 * ab_x[k] * is_5[k]
                  + f_21 * ab_x[k] * is_10[k]
                  - f_23 * ab_x[k] * is_12[k]
                  - f_20 * ab_x[k] * is_21[k]
                  + f_22 * ab_x[k] * is_23[k]
                  - f_20 * ks_0[k]
                  + f_21 * ks_3[k]
                  + f_22 * ks_5[k]
                  + f_21 * ks_10[k]
                  - f_23 * ks_12[k]
                  - f_20 * ks_21[k]
                  + f_22 * ks_23[k];
    }

#pragma omp simd aligned(ab_y, ab_z, is_2, is_7, is_16, ks_4, ks_5, ks_11, ks_12, ks_22, \
                         ks_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_4 * ab_y[k] * is_2[k]
                  - f_3 * ab_y[k] * is_7[k]
                  + f_2 * ab_y[k] * is_16[k]
                  + f_4 * ks_4[k]
                  - f_3 * ks_11[k]
                  + f_2 * ks_22[k];

        g_34[k] = f_4 * ab_z[k] * is_2[k]
                  - f_3 * ab_z[k] * is_7[k]
                  + f_2 * ab_z[k] * is_16[k]
                  + f_4 * ks_5[k]
                  - f_3 * ks_12[k]
                  + f_2 * ks_23[k];
    }

#pragma omp simd aligned(ab_x, is_2, is_7, is_16, ks_2, ks_7, ks_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_4 * ab_x[k] * is_2[k]
                  - f_3 * ab_x[k] * is_7[k]
                  + f_2 * ab_x[k] * is_16[k]
                  + f_4 * ks_2[k]
                  - f_3 * ks_7[k]
                  + f_2 * ks_16[k];
    }

#pragma omp simd aligned(ab_y, ab_z, is_0, is_3, is_10, is_21, ks_1, ks_2, ks_6, ks_7, ks_15, \
                         ks_16, ks_28, ks_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_24 * ab_y[k] * is_0[k]
                  - f_25 * ab_y[k] * is_3[k]
                  + f_25 * ab_y[k] * is_10[k]
                  - f_24 * ab_y[k] * is_21[k]
                  + f_24 * ks_1[k]
                  - f_25 * ks_6[k]
                  + f_25 * ks_15[k]
                  - f_24 * ks_28[k];

        g_37[k] = f_24 * ab_z[k] * is_0[k]
                  - f_25 * ab_z[k] * is_3[k]
                  + f_25 * ab_z[k] * is_10[k]
                  - f_24 * ab_z[k] * is_21[k]
                  + f_24 * ks_2[k]
                  - f_25 * ks_7[k]
                  + f_25 * ks_16[k]
                  - f_24 * ks_29[k];
    }

#pragma omp simd aligned(ab_x, is_0, is_3, is_10, is_21, ks_0, ks_3, ks_10, \
                         ks_21 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_24 * ab_x[k] * is_0[k]
                  - f_25 * ab_x[k] * is_3[k]
                  + f_25 * ab_x[k] * is_10[k]
                  - f_24 * ab_x[k] * is_21[k]
                  + f_24 * ks_0[k]
                  - f_25 * ks_3[k]
                  + f_25 * ks_10[k]
                  - f_24 * ks_21[k];
    }
}

}  // namespace simdtrf
