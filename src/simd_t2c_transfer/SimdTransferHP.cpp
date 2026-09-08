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


#include "SimdTransferHP.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_hp(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t hs, const size_t is, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);
    const auto *hs_16 = buffer.data(hs + 16);
    const auto *hs_17 = buffer.data(hs + 17);
    const auto *hs_18 = buffer.data(hs + 18);
    const auto *hs_19 = buffer.data(hs + 19);
    const auto *hs_20 = buffer.data(hs + 20);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, hs_0, hs_1, is_0, \
                         is_1, is_2, is_3, is_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * hs_0[k]
                 + is_0[k];

        t_1[k] = ab_y[k] * hs_0[k]
                 + is_1[k];

        t_2[k] = ab_z[k] * hs_0[k]
                 + is_2[k];

        t_3[k] = ab_x[k] * hs_1[k]
                 + is_1[k];

        t_4[k] = ab_y[k] * hs_1[k]
                 + is_3[k];

        t_5[k] = ab_z[k] * hs_1[k]
                 + is_4[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, hs_2, hs_3, is_2, is_3, \
                         is_4, is_5, is_6 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_6[k] = ab_x[k] * hs_2[k]
                 + is_2[k];

        t_7[k] = ab_y[k] * hs_2[k]
                 + is_4[k];

        t_8[k] = ab_z[k] * hs_2[k]
                 + is_5[k];

        t_9[k] = ab_x[k] * hs_3[k]
                 + is_3[k];

        t_10[k] = ab_y[k] * hs_3[k]
                  + is_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, hs_3, hs_4, \
                         hs_5, is_4, is_5, is_7, is_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_11[k] = ab_z[k] * hs_3[k]
                  + is_7[k];

        t_12[k] = ab_x[k] * hs_4[k]
                  + is_4[k];

        t_13[k] = ab_y[k] * hs_4[k]
                  + is_7[k];

        t_14[k] = ab_z[k] * hs_4[k]
                  + is_8[k];

        t_15[k] = ab_x[k] * hs_5[k]
                  + is_5[k];

        t_16[k] = ab_y[k] * hs_5[k]
                  + is_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, hs_5, hs_6, hs_7, \
                         is_6, is_7, is_9, is_10, is_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_17[k] = ab_z[k] * hs_5[k]
                  + is_9[k];

        t_18[k] = ab_x[k] * hs_6[k]
                  + is_6[k];

        t_19[k] = ab_y[k] * hs_6[k]
                  + is_10[k];

        t_20[k] = ab_z[k] * hs_6[k]
                  + is_11[k];

        t_21[k] = ab_x[k] * hs_7[k]
                  + is_7[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, hs_7, hs_8, is_8, \
                         is_11, is_12, is_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_22[k] = ab_y[k] * hs_7[k]
                  + is_11[k];

        t_23[k] = ab_z[k] * hs_7[k]
                  + is_12[k];

        t_24[k] = ab_x[k] * hs_8[k]
                  + is_8[k];

        t_25[k] = ab_y[k] * hs_8[k]
                  + is_12[k];

        t_26[k] = ab_z[k] * hs_8[k]
                  + is_13[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, hs_9, hs_10, is_9, \
                         is_10, is_13, is_14, is_15 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_27[k] = ab_x[k] * hs_9[k]
                  + is_9[k];

        t_28[k] = ab_y[k] * hs_9[k]
                  + is_13[k];

        t_29[k] = ab_z[k] * hs_9[k]
                  + is_14[k];

        t_30[k] = ab_x[k] * hs_10[k]
                  + is_10[k];

        t_31[k] = ab_y[k] * hs_10[k]
                  + is_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, ab_x, ab_y, ab_z, hs_10, hs_11, \
                         hs_12, is_11, is_12, is_16, is_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_32[k] = ab_z[k] * hs_10[k]
                  + is_16[k];

        t_33[k] = ab_x[k] * hs_11[k]
                  + is_11[k];

        t_34[k] = ab_y[k] * hs_11[k]
                  + is_16[k];

        t_35[k] = ab_z[k] * hs_11[k]
                  + is_17[k];

        t_36[k] = ab_x[k] * hs_12[k]
                  + is_12[k];

        t_37[k] = ab_y[k] * hs_12[k]
                  + is_17[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, hs_12, hs_13, \
                         hs_14, is_13, is_14, is_18, is_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_38[k] = ab_z[k] * hs_12[k]
                  + is_18[k];

        t_39[k] = ab_x[k] * hs_13[k]
                  + is_13[k];

        t_40[k] = ab_y[k] * hs_13[k]
                  + is_18[k];

        t_41[k] = ab_z[k] * hs_13[k]
                  + is_19[k];

        t_42[k] = ab_x[k] * hs_14[k]
                  + is_14[k];

        t_43[k] = ab_y[k] * hs_14[k]
                  + is_19[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, hs_14, hs_15, hs_16, \
                         is_15, is_16, is_20, is_21, is_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_44[k] = ab_z[k] * hs_14[k]
                  + is_20[k];

        t_45[k] = ab_x[k] * hs_15[k]
                  + is_15[k];

        t_46[k] = ab_y[k] * hs_15[k]
                  + is_21[k];

        t_47[k] = ab_z[k] * hs_15[k]
                  + is_22[k];

        t_48[k] = ab_x[k] * hs_16[k]
                  + is_16[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, hs_16, hs_17, is_17, \
                         is_22, is_23, is_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_49[k] = ab_y[k] * hs_16[k]
                  + is_22[k];

        t_50[k] = ab_z[k] * hs_16[k]
                  + is_23[k];

        t_51[k] = ab_x[k] * hs_17[k]
                  + is_17[k];

        t_52[k] = ab_y[k] * hs_17[k]
                  + is_23[k];

        t_53[k] = ab_z[k] * hs_17[k]
                  + is_24[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, hs_18, hs_19, \
                         is_18, is_19, is_24, is_25, is_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_54[k] = ab_x[k] * hs_18[k]
                  + is_18[k];

        t_55[k] = ab_y[k] * hs_18[k]
                  + is_24[k];

        t_56[k] = ab_z[k] * hs_18[k]
                  + is_25[k];

        t_57[k] = ab_x[k] * hs_19[k]
                  + is_19[k];

        t_58[k] = ab_y[k] * hs_19[k]
                  + is_25[k];

        t_59[k] = ab_z[k] * hs_19[k]
                  + is_26[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, ab_x, ab_y, ab_z, hs_20, is_20, is_26, \
                         is_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = ab_x[k] * hs_20[k]
                  + is_20[k];

        t_61[k] = ab_y[k] * hs_20[k]
                  + is_26[k];

        t_62[k] = ab_z[k] * hs_20[k]
                  + is_27[k];
    }
}

}  // namespace simdtrf
