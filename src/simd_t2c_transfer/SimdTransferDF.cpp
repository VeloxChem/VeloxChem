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


#include "SimdTransferDF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_df(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t pf, const size_t pg, const size_t nmax) -> void
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

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_12 = buffer.data(pf + 12);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_14 = buffer.data(pf + 14);
    const auto *pf_15 = buffer.data(pf + 15);
    const auto *pf_16 = buffer.data(pf + 16);
    const auto *pf_17 = buffer.data(pf + 17);
    const auto *pf_18 = buffer.data(pf + 18);
    const auto *pf_19 = buffer.data(pf + 19);
    const auto *pf_20 = buffer.data(pf + 20);
    const auto *pf_21 = buffer.data(pf + 21);
    const auto *pf_22 = buffer.data(pf + 22);
    const auto *pf_23 = buffer.data(pf + 23);
    const auto *pf_24 = buffer.data(pf + 24);
    const auto *pf_25 = buffer.data(pf + 25);
    const auto *pf_26 = buffer.data(pf + 26);
    const auto *pf_27 = buffer.data(pf + 27);
    const auto *pf_28 = buffer.data(pf + 28);
    const auto *pf_29 = buffer.data(pf + 29);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, pf_0, pf_1, pf_2, pf_3, pf_4, pg_0, \
                         pg_1, pg_2, pg_3, pg_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * pf_0[k]
                 + pg_0[k];

        t_1[k] = -ab_x[k] * pf_1[k]
                 + pg_1[k];

        t_2[k] = -ab_x[k] * pf_2[k]
                 + pg_2[k];

        t_3[k] = -ab_x[k] * pf_3[k]
                 + pg_3[k];

        t_4[k] = -ab_x[k] * pf_4[k]
                 + pg_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, pf_5, pf_6, pf_7, pf_8, pf_9, pg_5, \
                         pg_6, pg_7, pg_8, pg_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * pf_5[k]
                 + pg_5[k];

        t_6[k] = -ab_x[k] * pf_6[k]
                 + pg_6[k];

        t_7[k] = -ab_x[k] * pf_7[k]
                 + pg_7[k];

        t_8[k] = -ab_x[k] * pf_8[k]
                 + pg_8[k];

        t_9[k] = -ab_x[k] * pf_9[k]
                 + pg_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, pf_10, pf_11, pf_12, pf_13, \
                         pf_14, pg_15, pg_16, pg_17, pg_18, pg_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * pf_10[k]
                  + pg_15[k];

        t_11[k] = -ab_x[k] * pf_11[k]
                  + pg_16[k];

        t_12[k] = -ab_x[k] * pf_12[k]
                  + pg_17[k];

        t_13[k] = -ab_x[k] * pf_13[k]
                  + pg_18[k];

        t_14[k] = -ab_x[k] * pf_14[k]
                  + pg_19[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, pf_15, pf_16, pf_17, pf_18, \
                         pf_19, pg_20, pg_21, pg_22, pg_23, pg_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * pf_15[k]
                  + pg_20[k];

        t_16[k] = -ab_x[k] * pf_16[k]
                  + pg_21[k];

        t_17[k] = -ab_x[k] * pf_17[k]
                  + pg_22[k];

        t_18[k] = -ab_x[k] * pf_18[k]
                  + pg_23[k];

        t_19[k] = -ab_x[k] * pf_19[k]
                  + pg_24[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, pf_20, pf_21, pf_22, pf_23, \
                         pf_24, pg_30, pg_31, pg_32, pg_33, pg_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * pf_20[k]
                  + pg_30[k];

        t_21[k] = -ab_x[k] * pf_21[k]
                  + pg_31[k];

        t_22[k] = -ab_x[k] * pf_22[k]
                  + pg_32[k];

        t_23[k] = -ab_x[k] * pf_23[k]
                  + pg_33[k];

        t_24[k] = -ab_x[k] * pf_24[k]
                  + pg_34[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, pf_25, pf_26, pf_27, pf_28, \
                         pf_29, pg_35, pg_36, pg_37, pg_38, pg_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * pf_25[k]
                  + pg_35[k];

        t_26[k] = -ab_x[k] * pf_26[k]
                  + pg_36[k];

        t_27[k] = -ab_x[k] * pf_27[k]
                  + pg_37[k];

        t_28[k] = -ab_x[k] * pf_28[k]
                  + pg_38[k];

        t_29[k] = -ab_x[k] * pf_29[k]
                  + pg_39[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_y, pf_10, pf_11, pf_12, pf_13, \
                         pf_14, pg_16, pg_18, pg_19, pg_21, pg_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_y[k] * pf_10[k]
                  + pg_16[k];

        t_31[k] = -ab_y[k] * pf_11[k]
                  + pg_18[k];

        t_32[k] = -ab_y[k] * pf_12[k]
                  + pg_19[k];

        t_33[k] = -ab_y[k] * pf_13[k]
                  + pg_21[k];

        t_34[k] = -ab_y[k] * pf_14[k]
                  + pg_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_y, pf_15, pf_16, pf_17, pf_18, \
                         pf_19, pg_23, pg_25, pg_26, pg_27, pg_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_y[k] * pf_15[k]
                  + pg_23[k];

        t_36[k] = -ab_y[k] * pf_16[k]
                  + pg_25[k];

        t_37[k] = -ab_y[k] * pf_17[k]
                  + pg_26[k];

        t_38[k] = -ab_y[k] * pf_18[k]
                  + pg_27[k];

        t_39[k] = -ab_y[k] * pf_19[k]
                  + pg_28[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_y, pf_20, pf_21, pf_22, pf_23, \
                         pf_24, pg_31, pg_33, pg_34, pg_36, pg_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_y[k] * pf_20[k]
                  + pg_31[k];

        t_41[k] = -ab_y[k] * pf_21[k]
                  + pg_33[k];

        t_42[k] = -ab_y[k] * pf_22[k]
                  + pg_34[k];

        t_43[k] = -ab_y[k] * pf_23[k]
                  + pg_36[k];

        t_44[k] = -ab_y[k] * pf_24[k]
                  + pg_37[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_y, pf_25, pf_26, pf_27, pf_28, \
                         pf_29, pg_38, pg_40, pg_41, pg_42, pg_43 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = -ab_y[k] * pf_25[k]
                  + pg_38[k];

        t_46[k] = -ab_y[k] * pf_26[k]
                  + pg_40[k];

        t_47[k] = -ab_y[k] * pf_27[k]
                  + pg_41[k];

        t_48[k] = -ab_y[k] * pf_28[k]
                  + pg_42[k];

        t_49[k] = -ab_y[k] * pf_29[k]
                  + pg_43[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_z, pf_20, pf_21, pf_22, pf_23, \
                         pf_24, pg_32, pg_34, pg_35, pg_37, pg_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = -ab_z[k] * pf_20[k]
                  + pg_32[k];

        t_51[k] = -ab_z[k] * pf_21[k]
                  + pg_34[k];

        t_52[k] = -ab_z[k] * pf_22[k]
                  + pg_35[k];

        t_53[k] = -ab_z[k] * pf_23[k]
                  + pg_37[k];

        t_54[k] = -ab_z[k] * pf_24[k]
                  + pg_38[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_z, pf_25, pf_26, pf_27, pf_28, \
                         pf_29, pg_39, pg_41, pg_42, pg_43, pg_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = -ab_z[k] * pf_25[k]
                  + pg_39[k];

        t_56[k] = -ab_z[k] * pf_26[k]
                  + pg_41[k];

        t_57[k] = -ab_z[k] * pf_27[k]
                  + pg_42[k];

        t_58[k] = -ab_z[k] * pf_28[k]
                  + pg_43[k];

        t_59[k] = -ab_z[k] * pf_29[k]
                  + pg_44[k];
    }
}

}  // namespace simdtrf
