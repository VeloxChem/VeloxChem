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


#include "SimdElectronRepulsionCtrVrrIS.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ctr_is_electron_repulsion_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                                    const size_t pa, const size_t gs0, const size_t gs1,
                                    const size_t hs, const size_t ncols, const double alpha,
                                    const double beta, const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.625 * std::sqrt(462.0) / alpha;
    const auto f_1 = 0.625 * std::sqrt(462.0) * beta / (alpha * p);
    const auto f_2 = 0.1875 * std::sqrt(462.0);
    const auto f_3 = 0.625 * std::sqrt(462.0);
    const auto f_4 = 0.9375 * std::sqrt(154.0);
    const auto f_5 = 1.875 * std::sqrt(154.0);
    const auto f_6 = 0.1875 * std::sqrt(154.0);
    const auto f_7 = 0.75 * std::sqrt(7.0);
    const auto f_8 = 7.5 * std::sqrt(7.0);
    const auto f_9 = 0.5 * std::sqrt(210.0) / alpha;
    const auto f_10 = 0.5 * std::sqrt(210.0) * beta / (alpha * p);
    const auto f_11 = 0.5625 * std::sqrt(210.0);
    const auto f_12 = 0.375 * std::sqrt(210.0);
    const auto f_13 = 1.5 * std::sqrt(210.0);
    const auto f_14 = 0.1875 * std::sqrt(210.0);
    const auto f_15 = 0.5 * std::sqrt(210.0);
    const auto f_16 = 0.125 * std::sqrt(210.0) / alpha;
    const auto f_17 = 0.125 * std::sqrt(210.0) * beta / (alpha * p);
    const auto f_18 = 0.0625 * std::sqrt(210.0);
    const auto f_19 = std::sqrt(210.0);
    const auto f_20 = 0.125 * std::sqrt(210.0);
    const auto f_21 = 2.5 * std::sqrt(21.0) / alpha;
    const auto f_22 = 2.5 * std::sqrt(21.0) * beta / (alpha * p);
    const auto f_23 = 0.625 * std::sqrt(21.0);
    const auto f_24 = 1.25 * std::sqrt(21.0);
    const auto f_25 = 2.5 * std::sqrt(21.0);
    const auto f_26 = std::sqrt(21.0);
    const auto f_27 = 0.78125 / alpha;
    const auto f_28 = 1.40625 / alpha;
    const auto f_29 = 8.4375 / alpha;
    const auto f_30 = 1.25 / alpha;
    const auto f_31 = 14.0625 / alpha;
    const auto f_32 = 5.0 / alpha;
    const auto f_33 = 0.78125 * beta / (alpha * p);
    const auto f_34 = 1.40625 * beta / (alpha * p);
    const auto f_35 = 8.4375 * beta / (alpha * p);
    const auto f_36 = 1.25 * beta / (alpha * p);
    const auto f_37 = 14.0625 * beta / (alpha * p);
    const auto f_38 = 5.0 * beta / (alpha * p);
    const auto f_39 = 0.078125 * std::sqrt(210.0) / alpha;
    const auto f_40 = 0.046875 * std::sqrt(210.0) / alpha;
    const auto f_41 = 0.75 * std::sqrt(210.0) / alpha;
    const auto f_42 = 0.09375 * std::sqrt(210.0) / alpha;
    const auto f_43 = 0.078125 * std::sqrt(210.0) * beta / (alpha * p);
    const auto f_44 = 0.046875 * std::sqrt(210.0) * beta / (alpha * p);
    const auto f_45 = 0.75 * std::sqrt(210.0) * beta / (alpha * p);
    const auto f_46 = 0.09375 * std::sqrt(210.0) * beta / (alpha * p);
    const auto f_47 = 0.03125 * std::sqrt(210.0);
    const auto f_48 = 0.46875 * std::sqrt(7.0) / alpha;
    const auto f_49 = 1.40625 * std::sqrt(7.0) / alpha;
    const auto f_50 = 2.8125 * std::sqrt(7.0) / alpha;
    const auto f_51 = 0.46875 * std::sqrt(7.0) * beta / (alpha * p);
    const auto f_52 = 1.40625 * std::sqrt(7.0) * beta / (alpha * p);
    const auto f_53 = 2.8125 * std::sqrt(7.0) * beta / (alpha * p);
    const auto f_54 = 0.1875 * std::sqrt(7.0);
    const auto f_55 = 0.9375 * std::sqrt(7.0);
    const auto f_56 = 1.875 * std::sqrt(7.0);
    const auto f_57 = 11.25 * std::sqrt(7.0);
    const auto f_58 = 0.078125 * std::sqrt(462.0) / alpha;
    const auto f_59 = 0.703125 * std::sqrt(462.0) / alpha;
    const auto f_60 = 0.15625 * std::sqrt(462.0) / alpha;
    const auto f_61 = 0.078125 * std::sqrt(462.0) * beta / (alpha * p);
    const auto f_62 = 0.703125 * std::sqrt(462.0) * beta / (alpha * p);
    const auto f_63 = 0.15625 * std::sqrt(462.0) * beta / (alpha * p);
    const auto f_64 = 0.03125 * std::sqrt(462.0);
    const auto f_65 = 0.46875 * std::sqrt(462.0);

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

    const auto *gs0_0 = buffer.data(gs0 + 0);
    const auto *gs0_1 = buffer.data(gs0 + 1);
    const auto *gs0_2 = buffer.data(gs0 + 2);
    const auto *gs0_3 = buffer.data(gs0 + 3);
    const auto *gs0_4 = buffer.data(gs0 + 4);
    const auto *gs0_5 = buffer.data(gs0 + 5);
    const auto *gs0_6 = buffer.data(gs0 + 6);
    const auto *gs0_7 = buffer.data(gs0 + 7);
    const auto *gs0_8 = buffer.data(gs0 + 8);

    const auto *gs1_0 = buffer.data(gs1 + 0);
    const auto *gs1_1 = buffer.data(gs1 + 1);
    const auto *gs1_2 = buffer.data(gs1 + 2);
    const auto *gs1_3 = buffer.data(gs1 + 3);
    const auto *gs1_4 = buffer.data(gs1 + 4);
    const auto *gs1_5 = buffer.data(gs1 + 5);
    const auto *gs1_6 = buffer.data(gs1 + 6);
    const auto *gs1_7 = buffer.data(gs1 + 7);
    const auto *gs1_8 = buffer.data(gs1 + 8);

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

#pragma omp simd aligned(pa_x, pa_y, pa_z, gs0_3, gs1_3, hs_0, hs_1, hs_3, hs_4, hs_9, \
                         hs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += -f_0 * gs0_3[k]
                  + f_1 * gs1_3[k]
                  + f_2 * pa_y[k] * hs_0[k]
                  - f_3 * pa_x[k] * hs_4[k]
                  + f_2 * pa_x[k] * hs_9[k];

        g_1[k] += f_4 * pa_y[k] * hs_1[k]
                  - f_5 * pa_z[k] * hs_4[k]
                  + f_6 * pa_z[k] * hs_9[k];

        g_2[k] += -f_7 * pa_y[k] * hs_0[k]
                  + f_8 * pa_y[k] * hs_3[k]
                  + f_7 * pa_x[k] * hs_9[k]
                  - f_8 * pa_x[k] * hs_11[k];
    }

#pragma omp simd aligned(pa_y, pa_z, gs0_7, gs1_7, hs_1, hs_4, hs_5, hs_9, \
                         hs_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += -f_9 * gs0_7[k]
                  + f_10 * gs1_7[k]
                  - f_11 * pa_y[k] * hs_1[k]
                  - f_12 * pa_z[k] * hs_4[k]
                  + f_13 * pa_y[k] * hs_5[k]
                  + f_14 * pa_z[k] * hs_9[k]
                  - f_15 * pa_y[k] * hs_12[k];
    }

#pragma omp simd aligned(pa_x, pa_y, gs0_3, gs1_3, hs_0, hs_3, hs_4, hs_9, hs_11, \
                         hs_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += f_16 * gs0_3[k]
                  - f_17 * gs1_3[k]
                  + f_18 * pa_y[k] * hs_0[k]
                  - f_19 * pa_y[k] * hs_3[k]
                  + f_20 * pa_x[k] * hs_4[k]
                  + f_18 * pa_x[k] * hs_9[k]
                  - f_19 * pa_x[k] * hs_11[k]
                  + f_19 * pa_x[k] * hs_13[k];
    }

#pragma omp simd aligned(pa_y, pa_z, gs0_7, gs1_7, hs_1, hs_4, hs_5, hs_9, hs_12, \
                         hs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += -f_21 * gs0_7[k]
                  + f_22 * gs1_7[k]
                  + f_23 * pa_y[k] * hs_1[k]
                  + f_24 * pa_z[k] * hs_4[k]
                  - f_25 * pa_y[k] * hs_5[k]
                  + f_23 * pa_z[k] * hs_9[k]
                  - f_25 * pa_y[k] * hs_12[k]
                  + f_26 * pa_y[k] * hs_14[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, gs0_0, gs0_1, gs0_2, gs0_5, gs0_6, gs0_8, gs1_0, \
                         gs1_1, gs1_2, gs1_5, gs1_6, gs1_8, hs_0, hs_2, hs_3, hs_6, hs_7, \
                         hs_8, hs_9, hs_11, hs_13, hs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += -f_27 * gs0_0[k]
                  - f_28 * gs0_1[k]
                  + f_29 * gs0_2[k]
                  - f_30 * gs0_5[k]
                  + f_31 * gs0_6[k]
                  - f_32 * gs0_8[k]
                  + f_33 * gs1_0[k]
                  + f_34 * gs1_1[k]
                  - f_35 * gs1_2[k]
                  + f_36 * gs1_5[k]
                  - f_37 * gs1_6[k]
                  + f_38 * gs1_8[k]
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

#pragma omp simd aligned(pa_x, pa_z, gs0_4, gs1_4, hs_0, hs_2, hs_5, hs_10, hs_12, \
                         hs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_7[k] += -f_21 * gs0_4[k]
                  + f_22 * gs1_4[k]
                  + f_23 * pa_z[k] * hs_0[k]
                  + f_24 * pa_z[k] * hs_2[k]
                  - f_25 * pa_x[k] * hs_5[k]
                  + f_23 * pa_x[k] * hs_10[k]
                  - f_25 * pa_x[k] * hs_12[k]
                  + f_26 * pa_x[k] * hs_14[k];
    }

#pragma omp simd aligned(pa_x, pa_y, gs0_0, gs0_1, gs0_2, gs0_5, gs0_6, gs1_0, gs1_1, gs1_2, \
                         gs1_5, gs1_6, hs_0, hs_2, hs_3, hs_6, hs_8, hs_9, hs_11, \
                         hs_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += f_39 * gs0_0[k]
                  + f_40 * gs0_1[k]
                  - f_41 * gs0_2[k]
                  - f_42 * gs0_5[k]
                  + f_41 * gs0_6[k]
                  - f_43 * gs1_0[k]
                  - f_44 * gs1_1[k]
                  + f_45 * gs1_2[k]
                  + f_46 * gs1_5[k]
                  - f_45 * gs1_6[k]
                  + f_47 * pa_x[k] * hs_0[k]
                  + f_47 * pa_x[k] * hs_2[k]
                  - f_15 * pa_x[k] * hs_3[k]
                  - f_47 * pa_x[k] * hs_6[k]
                  + f_15 * pa_x[k] * hs_8[k]
                  - f_47 * pa_y[k] * hs_9[k]
                  + f_15 * pa_y[k] * hs_11[k]
                  - f_15 * pa_y[k] * hs_13[k];
    }

#pragma omp simd aligned(pa_x, pa_z, gs0_4, gs1_4, hs_0, hs_2, hs_5, hs_10, \
                         hs_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_9[k] += f_9 * gs0_4[k]
                  - f_10 * gs1_4[k]
                  - f_14 * pa_z[k] * hs_0[k]
                  + f_12 * pa_z[k] * hs_2[k]
                  + f_15 * pa_x[k] * hs_5[k]
                  + f_11 * pa_x[k] * hs_10[k]
                  - f_13 * pa_x[k] * hs_12[k];
    }

#pragma omp simd aligned(pa_x, pa_y, gs0_0, gs0_1, gs0_2, gs0_6, gs1_0, gs1_1, gs1_2, gs1_6, \
                         hs_0, hs_2, hs_3, hs_6, hs_7, hs_9, hs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_10[k] += -f_48 * gs0_0[k]
                   + f_49 * gs0_1[k]
                   + f_50 * gs0_2[k]
                   - f_50 * gs0_6[k]
                   + f_51 * gs1_0[k]
                   - f_52 * gs1_1[k]
                   - f_53 * gs1_2[k]
                   + f_53 * gs1_6[k]
                   - f_54 * pa_x[k] * hs_0[k]
                   + f_55 * pa_x[k] * hs_2[k]
                   + f_56 * pa_x[k] * hs_3[k]
                   + f_55 * pa_x[k] * hs_6[k]
                   - f_57 * pa_x[k] * hs_7[k]
                   - f_54 * pa_y[k] * hs_9[k]
                   + f_56 * pa_y[k] * hs_11[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, gs0_0, gs0_1, gs0_5, gs1_0, gs1_1, gs1_5, hs_0, \
                         hs_2, hs_6, hs_9, hs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_11[k] += f_6 * pa_z[k] * hs_0[k]
                   - f_5 * pa_z[k] * hs_2[k]
                   + f_4 * pa_x[k] * hs_10[k];

        g_12[k] += f_58 * gs0_0[k]
                   - f_59 * gs0_1[k]
                   + f_60 * gs0_5[k]
                   - f_61 * gs1_0[k]
                   + f_62 * gs1_1[k]
                   - f_63 * gs1_5[k]
                   + f_64 * pa_x[k] * hs_0[k]
                   - f_65 * pa_x[k] * hs_2[k]
                   + f_65 * pa_x[k] * hs_6[k]
                   - f_64 * pa_y[k] * hs_9[k];
    }
}

}  // namespace simdt2ceri
