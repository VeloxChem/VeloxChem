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


#include "SimdTransferFD.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_fd(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t fp, const size_t gp, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *fp_0 = buffer.data(fp + 0);
    const auto *fp_1 = buffer.data(fp + 1);
    const auto *fp_2 = buffer.data(fp + 2);
    const auto *fp_3 = buffer.data(fp + 3);
    const auto *fp_4 = buffer.data(fp + 4);
    const auto *fp_5 = buffer.data(fp + 5);
    const auto *fp_6 = buffer.data(fp + 6);
    const auto *fp_7 = buffer.data(fp + 7);
    const auto *fp_8 = buffer.data(fp + 8);
    const auto *fp_9 = buffer.data(fp + 9);
    const auto *fp_10 = buffer.data(fp + 10);
    const auto *fp_11 = buffer.data(fp + 11);
    const auto *fp_12 = buffer.data(fp + 12);
    const auto *fp_13 = buffer.data(fp + 13);
    const auto *fp_14 = buffer.data(fp + 14);
    const auto *fp_15 = buffer.data(fp + 15);
    const auto *fp_16 = buffer.data(fp + 16);
    const auto *fp_17 = buffer.data(fp + 17);
    const auto *fp_18 = buffer.data(fp + 18);
    const auto *fp_19 = buffer.data(fp + 19);
    const auto *fp_20 = buffer.data(fp + 20);
    const auto *fp_21 = buffer.data(fp + 21);
    const auto *fp_22 = buffer.data(fp + 22);
    const auto *fp_23 = buffer.data(fp + 23);
    const auto *fp_24 = buffer.data(fp + 24);
    const auto *fp_25 = buffer.data(fp + 25);
    const auto *fp_26 = buffer.data(fp + 26);
    const auto *fp_27 = buffer.data(fp + 27);
    const auto *fp_28 = buffer.data(fp + 28);
    const auto *fp_29 = buffer.data(fp + 29);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);
    const auto *gp_22 = buffer.data(gp + 22);
    const auto *gp_23 = buffer.data(gp + 23);
    const auto *gp_24 = buffer.data(gp + 24);
    const auto *gp_25 = buffer.data(gp + 25);
    const auto *gp_26 = buffer.data(gp + 26);
    const auto *gp_27 = buffer.data(gp + 27);
    const auto *gp_28 = buffer.data(gp + 28);
    const auto *gp_29 = buffer.data(gp + 29);
    const auto *gp_31 = buffer.data(gp + 31);
    const auto *gp_32 = buffer.data(gp + 32);
    const auto *gp_34 = buffer.data(gp + 34);
    const auto *gp_35 = buffer.data(gp + 35);
    const auto *gp_37 = buffer.data(gp + 37);
    const auto *gp_38 = buffer.data(gp + 38);
    const auto *gp_40 = buffer.data(gp + 40);
    const auto *gp_41 = buffer.data(gp + 41);
    const auto *gp_44 = buffer.data(gp + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, fp_0, fp_1, fp_2, gp_0, gp_1, \
                         gp_2, gp_4, gp_5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * fp_0[k]
                 + gp_0[k];

        t_1[k] = ab_x[k] * fp_1[k]
                 + gp_1[k];

        t_2[k] = ab_x[k] * fp_2[k]
                 + gp_2[k];

        t_3[k] = ab_y[k] * fp_1[k]
                 + gp_4[k];

        t_4[k] = ab_y[k] * fp_2[k]
                 + gp_5[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_z, fp_2, fp_3, fp_4, fp_5, gp_3, gp_4, \
                         gp_5, gp_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = ab_z[k] * fp_2[k]
                 + gp_8[k];

        t_6[k] = ab_x[k] * fp_3[k]
                 + gp_3[k];

        t_7[k] = ab_x[k] * fp_4[k]
                 + gp_4[k];

        t_8[k] = ab_x[k] * fp_5[k]
                 + gp_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_y, ab_z, fp_4, fp_5, fp_6, gp_6, \
                         gp_10, gp_11, gp_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_9[k] = ab_y[k] * fp_4[k]
                 + gp_10[k];

        t_10[k] = ab_y[k] * fp_5[k]
                  + gp_11[k];

        t_11[k] = ab_z[k] * fp_5[k]
                  + gp_14[k];

        t_12[k] = ab_x[k] * fp_6[k]
                  + gp_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, fp_7, fp_8, gp_7, \
                         gp_8, gp_13, gp_14, gp_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_13[k] = ab_x[k] * fp_7[k]
                  + gp_7[k];

        t_14[k] = ab_x[k] * fp_8[k]
                  + gp_8[k];

        t_15[k] = ab_y[k] * fp_7[k]
                  + gp_13[k];

        t_16[k] = ab_y[k] * fp_8[k]
                  + gp_14[k];

        t_17[k] = ab_z[k] * fp_8[k]
                  + gp_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, fp_9, fp_10, fp_11, gp_9, \
                         gp_10, gp_11, gp_19, gp_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_18[k] = ab_x[k] * fp_9[k]
                  + gp_9[k];

        t_19[k] = ab_x[k] * fp_10[k]
                  + gp_10[k];

        t_20[k] = ab_x[k] * fp_11[k]
                  + gp_11[k];

        t_21[k] = ab_y[k] * fp_10[k]
                  + gp_19[k];

        t_22[k] = ab_y[k] * fp_11[k]
                  + gp_20[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_z, fp_11, fp_12, fp_13, fp_14, \
                         gp_12, gp_13, gp_14, gp_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_23[k] = ab_z[k] * fp_11[k]
                  + gp_23[k];

        t_24[k] = ab_x[k] * fp_12[k]
                  + gp_12[k];

        t_25[k] = ab_x[k] * fp_13[k]
                  + gp_13[k];

        t_26[k] = ab_x[k] * fp_14[k]
                  + gp_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, fp_13, fp_14, fp_15, gp_15, \
                         gp_22, gp_23, gp_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_27[k] = ab_y[k] * fp_13[k]
                  + gp_22[k];

        t_28[k] = ab_y[k] * fp_14[k]
                  + gp_23[k];

        t_29[k] = ab_z[k] * fp_14[k]
                  + gp_26[k];

        t_30[k] = ab_x[k] * fp_15[k]
                  + gp_15[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, fp_16, fp_17, gp_16, \
                         gp_17, gp_25, gp_26, gp_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_31[k] = ab_x[k] * fp_16[k]
                  + gp_16[k];

        t_32[k] = ab_x[k] * fp_17[k]
                  + gp_17[k];

        t_33[k] = ab_y[k] * fp_16[k]
                  + gp_25[k];

        t_34[k] = ab_y[k] * fp_17[k]
                  + gp_26[k];

        t_35[k] = ab_z[k] * fp_17[k]
                  + gp_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_x, ab_y, fp_18, fp_19, fp_20, gp_18, \
                         gp_19, gp_20, gp_31, gp_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_36[k] = ab_x[k] * fp_18[k]
                  + gp_18[k];

        t_37[k] = ab_x[k] * fp_19[k]
                  + gp_19[k];

        t_38[k] = ab_x[k] * fp_20[k]
                  + gp_20[k];

        t_39[k] = ab_y[k] * fp_19[k]
                  + gp_31[k];

        t_40[k] = ab_y[k] * fp_20[k]
                  + gp_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_x, ab_z, fp_20, fp_21, fp_22, fp_23, \
                         gp_21, gp_22, gp_23, gp_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_41[k] = ab_z[k] * fp_20[k]
                  + gp_35[k];

        t_42[k] = ab_x[k] * fp_21[k]
                  + gp_21[k];

        t_43[k] = ab_x[k] * fp_22[k]
                  + gp_22[k];

        t_44[k] = ab_x[k] * fp_23[k]
                  + gp_23[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, fp_22, fp_23, fp_24, gp_24, \
                         gp_34, gp_35, gp_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = ab_y[k] * fp_22[k]
                  + gp_34[k];

        t_46[k] = ab_y[k] * fp_23[k]
                  + gp_35[k];

        t_47[k] = ab_z[k] * fp_23[k]
                  + gp_38[k];

        t_48[k] = ab_x[k] * fp_24[k]
                  + gp_24[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, fp_25, fp_26, gp_25, \
                         gp_26, gp_37, gp_38, gp_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_49[k] = ab_x[k] * fp_25[k]
                  + gp_25[k];

        t_50[k] = ab_x[k] * fp_26[k]
                  + gp_26[k];

        t_51[k] = ab_y[k] * fp_25[k]
                  + gp_37[k];

        t_52[k] = ab_y[k] * fp_26[k]
                  + gp_38[k];

        t_53[k] = ab_z[k] * fp_26[k]
                  + gp_41[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_x, ab_y, fp_27, fp_28, fp_29, gp_27, \
                         gp_28, gp_29, gp_40, gp_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_54[k] = ab_x[k] * fp_27[k]
                  + gp_27[k];

        t_55[k] = ab_x[k] * fp_28[k]
                  + gp_28[k];

        t_56[k] = ab_x[k] * fp_29[k]
                  + gp_29[k];

        t_57[k] = ab_y[k] * fp_28[k]
                  + gp_40[k];

        t_58[k] = ab_y[k] * fp_29[k]
                  + gp_41[k];
    }

#pragma omp simd aligned(t_59, ab_z, fp_29, gp_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_59[k] = ab_z[k] * fp_29[k]
                  + gp_44[k];
    }
}

}  // namespace simdtrf
