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


#include "SimdTransferDD.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_dd(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t pd, const size_t pf, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_7 = buffer.data(pd + 7);
    const auto *pd_8 = buffer.data(pd + 8);
    const auto *pd_9 = buffer.data(pd + 9);
    const auto *pd_10 = buffer.data(pd + 10);
    const auto *pd_11 = buffer.data(pd + 11);
    const auto *pd_12 = buffer.data(pd + 12);
    const auto *pd_13 = buffer.data(pd + 13);
    const auto *pd_14 = buffer.data(pd + 14);
    const auto *pd_15 = buffer.data(pd + 15);
    const auto *pd_16 = buffer.data(pd + 16);
    const auto *pd_17 = buffer.data(pd + 17);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_12 = buffer.data(pf + 12);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_14 = buffer.data(pf + 14);
    const auto *pf_15 = buffer.data(pf + 15);
    const auto *pf_16 = buffer.data(pf + 16);
    const auto *pf_17 = buffer.data(pf + 17);
    const auto *pf_18 = buffer.data(pf + 18);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, pd_0, pd_1, pd_2, pd_3, pd_4, pf_0, \
                         pf_1, pf_2, pf_3, pf_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * pd_0[k]
                 + pf_0[k];

        t_1[k] = -ab_x[k] * pd_1[k]
                 + pf_1[k];

        t_2[k] = -ab_x[k] * pd_2[k]
                 + pf_2[k];

        t_3[k] = -ab_x[k] * pd_3[k]
                 + pf_3[k];

        t_4[k] = -ab_x[k] * pd_4[k]
                 + pf_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, pd_5, pd_6, pd_7, pd_8, pd_9, pf_5, \
                         pf_10, pf_11, pf_12, pf_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * pd_5[k]
                 + pf_5[k];

        t_6[k] = -ab_x[k] * pd_6[k]
                 + pf_10[k];

        t_7[k] = -ab_x[k] * pd_7[k]
                 + pf_11[k];

        t_8[k] = -ab_x[k] * pd_8[k]
                 + pf_12[k];

        t_9[k] = -ab_x[k] * pd_9[k]
                 + pf_13[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, pd_10, pd_11, pd_12, pd_13, \
                         pd_14, pf_14, pf_15, pf_20, pf_21, pf_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * pd_10[k]
                  + pf_14[k];

        t_11[k] = -ab_x[k] * pd_11[k]
                  + pf_15[k];

        t_12[k] = -ab_x[k] * pd_12[k]
                  + pf_20[k];

        t_13[k] = -ab_x[k] * pd_13[k]
                  + pf_21[k];

        t_14[k] = -ab_x[k] * pd_14[k]
                  + pf_22[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, ab_x, ab_y, pd_6, pd_15, pd_16, pd_17, pf_11, \
                         pf_23, pf_24, pf_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * pd_15[k]
                  + pf_23[k];

        t_16[k] = -ab_x[k] * pd_16[k]
                  + pf_24[k];

        t_17[k] = -ab_x[k] * pd_17[k]
                  + pf_25[k];

        t_18[k] = -ab_y[k] * pd_6[k]
                  + pf_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, ab_y, pd_7, pd_8, pd_9, pd_10, pd_11, \
                         pf_13, pf_14, pf_16, pf_17, pf_18 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_19[k] = -ab_y[k] * pd_7[k]
                  + pf_13[k];

        t_20[k] = -ab_y[k] * pd_8[k]
                  + pf_14[k];

        t_21[k] = -ab_y[k] * pd_9[k]
                  + pf_16[k];

        t_22[k] = -ab_y[k] * pd_10[k]
                  + pf_17[k];

        t_23[k] = -ab_y[k] * pd_11[k]
                  + pf_18[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, ab_y, pd_12, pd_13, pd_14, pd_15, \
                         pd_16, pf_21, pf_23, pf_24, pf_26, pf_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_24[k] = -ab_y[k] * pd_12[k]
                  + pf_21[k];

        t_25[k] = -ab_y[k] * pd_13[k]
                  + pf_23[k];

        t_26[k] = -ab_y[k] * pd_14[k]
                  + pf_24[k];

        t_27[k] = -ab_y[k] * pd_15[k]
                  + pf_26[k];

        t_28[k] = -ab_y[k] * pd_16[k]
                  + pf_27[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, ab_y, ab_z, pd_12, pd_13, pd_14, pd_17, \
                         pf_22, pf_24, pf_25, pf_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_29[k] = -ab_y[k] * pd_17[k]
                  + pf_28[k];

        t_30[k] = -ab_z[k] * pd_12[k]
                  + pf_22[k];

        t_31[k] = -ab_z[k] * pd_13[k]
                  + pf_24[k];

        t_32[k] = -ab_z[k] * pd_14[k]
                  + pf_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, ab_z, pd_15, pd_16, pd_17, pf_27, pf_28, \
                         pf_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_33[k] = -ab_z[k] * pd_15[k]
                  + pf_27[k];

        t_34[k] = -ab_z[k] * pd_16[k]
                  + pf_28[k];

        t_35[k] = -ab_z[k] * pd_17[k]
                  + pf_29[k];
    }
}

}  // namespace simdtrf
