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


#include "SimdOverlapCtrVrrIS.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_ctr_is_overlap_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                         const size_t pa, const size_t gs, const size_t hs, const size_t ncols,
                         const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.625 * std::sqrt(462.0) / p;
    const auto f_1 = 0.1875 * std::sqrt(462.0);
    const auto f_2 = 0.625 * std::sqrt(462.0);
    const auto f_3 = 0.9375 * std::sqrt(154.0);
    const auto f_4 = 1.875 * std::sqrt(154.0);
    const auto f_5 = 0.1875 * std::sqrt(154.0);
    const auto f_6 = 0.75 * std::sqrt(7.0);
    const auto f_7 = 7.5 * std::sqrt(7.0);
    const auto f_8 = 0.5 * std::sqrt(210.0) / p;
    const auto f_9 = 0.5625 * std::sqrt(210.0);
    const auto f_10 = 0.375 * std::sqrt(210.0);
    const auto f_11 = 1.5 * std::sqrt(210.0);
    const auto f_12 = 0.1875 * std::sqrt(210.0);
    const auto f_13 = 0.5 * std::sqrt(210.0);
    const auto f_14 = 0.125 * std::sqrt(210.0) / p;
    const auto f_15 = 0.0625 * std::sqrt(210.0);
    const auto f_16 = std::sqrt(210.0);
    const auto f_17 = 0.125 * std::sqrt(210.0);
    const auto f_18 = 2.5 * std::sqrt(21.0) / p;
    const auto f_19 = 0.625 * std::sqrt(21.0);
    const auto f_20 = 1.25 * std::sqrt(21.0);
    const auto f_21 = 2.5 * std::sqrt(21.0);
    const auto f_22 = std::sqrt(21.0);
    const auto f_23 = 0.78125 / p;
    const auto f_24 = 1.40625 / p;
    const auto f_25 = 8.4375 / p;
    const auto f_26 = 1.25 / p;
    const auto f_27 = 14.0625 / p;
    const auto f_28 = 5.0 / p;
    const auto f_29 = 0.078125 * std::sqrt(210.0) / p;
    const auto f_30 = 0.046875 * std::sqrt(210.0) / p;
    const auto f_31 = 0.75 * std::sqrt(210.0) / p;
    const auto f_32 = 0.09375 * std::sqrt(210.0) / p;
    const auto f_33 = 0.03125 * std::sqrt(210.0);
    const auto f_34 = 0.46875 * std::sqrt(7.0) / p;
    const auto f_35 = 1.40625 * std::sqrt(7.0) / p;
    const auto f_36 = 2.8125 * std::sqrt(7.0) / p;
    const auto f_37 = 0.1875 * std::sqrt(7.0);
    const auto f_38 = 0.9375 * std::sqrt(7.0);
    const auto f_39 = 1.875 * std::sqrt(7.0);
    const auto f_40 = 11.25 * std::sqrt(7.0);
    const auto f_41 = 0.078125 * std::sqrt(462.0) / p;
    const auto f_42 = 0.703125 * std::sqrt(462.0) / p;
    const auto f_43 = 0.15625 * std::sqrt(462.0) / p;
    const auto f_44 = 0.03125 * std::sqrt(462.0);
    const auto f_45 = 0.46875 * std::sqrt(462.0);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);

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

#pragma omp simd aligned(pa_x, pa_y, pa_z, gs_3, hs_0, hs_1, hs_3, hs_4, hs_9, \
                         hs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += -f_0 * gs_3[k]
                  + f_1 * pa_y[k] * hs_0[k]
                  - f_2 * pa_x[k] * hs_4[k]
                  + f_1 * pa_x[k] * hs_9[k];

        g_1[k] += f_3 * pa_y[k] * hs_1[k]
                  - f_4 * pa_z[k] * hs_4[k]
                  + f_5 * pa_z[k] * hs_9[k];

        g_2[k] += -f_6 * pa_y[k] * hs_0[k]
                  + f_7 * pa_y[k] * hs_3[k]
                  + f_6 * pa_x[k] * hs_9[k]
                  - f_7 * pa_x[k] * hs_11[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, gs_3, gs_7, hs_0, hs_1, hs_3, hs_4, hs_5, hs_9, \
                         hs_11, hs_12, hs_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += -f_8 * gs_7[k]
                  - f_9 * pa_y[k] * hs_1[k]
                  - f_10 * pa_z[k] * hs_4[k]
                  + f_11 * pa_y[k] * hs_5[k]
                  + f_12 * pa_z[k] * hs_9[k]
                  - f_13 * pa_y[k] * hs_12[k];

        g_4[k] += f_14 * gs_3[k]
                  + f_15 * pa_y[k] * hs_0[k]
                  - f_16 * pa_y[k] * hs_3[k]
                  + f_17 * pa_x[k] * hs_4[k]
                  + f_15 * pa_x[k] * hs_9[k]
                  - f_16 * pa_x[k] * hs_11[k]
                  + f_16 * pa_x[k] * hs_13[k];
    }

#pragma omp simd aligned(pa_y, pa_z, gs_7, hs_1, hs_4, hs_5, hs_9, hs_12, \
                         hs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += -f_18 * gs_7[k]
                  + f_19 * pa_y[k] * hs_1[k]
                  + f_20 * pa_z[k] * hs_4[k]
                  - f_21 * pa_y[k] * hs_5[k]
                  + f_19 * pa_z[k] * hs_9[k]
                  - f_21 * pa_y[k] * hs_12[k]
                  + f_22 * pa_y[k] * hs_14[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, gs_0, gs_1, gs_2, gs_5, gs_6, gs_8, hs_0, hs_2, \
                         hs_3, hs_6, hs_7, hs_8, hs_9, hs_11, hs_13, \
                         hs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += -f_23 * gs_0[k]
                  - f_24 * gs_1[k]
                  + f_25 * gs_2[k]
                  - f_26 * gs_5[k]
                  + f_27 * gs_6[k]
                  - f_28 * gs_8[k]
                  - 0.3125 * pa_x[k] * hs_0[k]
                  - 0.9375 * pa_x[k] * hs_2[k]
                  + 5.625 * pa_x[k] * hs_3[k]
                  - 0.9375 * pa_x[k] * hs_6[k]
                  + 11.25 * pa_x[k] * hs_7[k]
                  - 7.5 * pa_x[k] * hs_8[k]
                  - 0.3125 * pa_y[k] * hs_9[k]
                  + 5.625 * pa_y[k] * hs_11[k]
                  - 7.5 * pa_y[k] * hs_13[k]
                  + pa_z[k] * hs_14[k];
    }

#pragma omp simd aligned(pa_x, pa_z, gs_4, hs_0, hs_2, hs_5, hs_10, hs_12, \
                         hs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_7[k] += -f_18 * gs_4[k]
                  + f_19 * pa_z[k] * hs_0[k]
                  + f_20 * pa_z[k] * hs_2[k]
                  - f_21 * pa_x[k] * hs_5[k]
                  + f_19 * pa_x[k] * hs_10[k]
                  - f_21 * pa_x[k] * hs_12[k]
                  + f_22 * pa_x[k] * hs_14[k];
    }

#pragma omp simd aligned(pa_x, pa_y, gs_0, gs_1, gs_2, gs_5, gs_6, hs_0, hs_2, hs_3, hs_6, \
                         hs_8, hs_9, hs_11, hs_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += f_29 * gs_0[k]
                  + f_30 * gs_1[k]
                  - f_31 * gs_2[k]
                  - f_32 * gs_5[k]
                  + f_31 * gs_6[k]
                  + f_33 * pa_x[k] * hs_0[k]
                  + f_33 * pa_x[k] * hs_2[k]
                  - f_13 * pa_x[k] * hs_3[k]
                  - f_33 * pa_x[k] * hs_6[k]
                  + f_13 * pa_x[k] * hs_8[k]
                  - f_33 * pa_y[k] * hs_9[k]
                  + f_13 * pa_y[k] * hs_11[k]
                  - f_13 * pa_y[k] * hs_13[k];
    }

#pragma omp simd aligned(pa_x, pa_z, gs_4, hs_0, hs_2, hs_5, hs_10, \
                         hs_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_9[k] += f_8 * gs_4[k]
                  - f_12 * pa_z[k] * hs_0[k]
                  + f_10 * pa_z[k] * hs_2[k]
                  + f_13 * pa_x[k] * hs_5[k]
                  + f_9 * pa_x[k] * hs_10[k]
                  - f_11 * pa_x[k] * hs_12[k];
    }

#pragma omp simd aligned(pa_x, pa_y, gs_0, gs_1, gs_2, gs_6, hs_0, hs_2, hs_3, hs_6, hs_7, \
                         hs_9, hs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_10[k] += -f_34 * gs_0[k]
                   + f_35 * gs_1[k]
                   + f_36 * gs_2[k]
                   - f_36 * gs_6[k]
                   - f_37 * pa_x[k] * hs_0[k]
                   + f_38 * pa_x[k] * hs_2[k]
                   + f_39 * pa_x[k] * hs_3[k]
                   + f_38 * pa_x[k] * hs_6[k]
                   - f_40 * pa_x[k] * hs_7[k]
                   - f_37 * pa_y[k] * hs_9[k]
                   + f_39 * pa_y[k] * hs_11[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, gs_0, gs_1, gs_5, hs_0, hs_2, hs_6, hs_9, \
                         hs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_11[k] += f_5 * pa_z[k] * hs_0[k]
                   - f_4 * pa_z[k] * hs_2[k]
                   + f_3 * pa_x[k] * hs_10[k];

        g_12[k] += f_41 * gs_0[k]
                   - f_42 * gs_1[k]
                   + f_43 * gs_5[k]
                   + f_44 * pa_x[k] * hs_0[k]
                   - f_45 * pa_x[k] * hs_2[k]
                   + f_45 * pa_x[k] * hs_6[k]
                   - f_44 * pa_y[k] * hs_9[k];
    }
}

}  // namespace simdovl
