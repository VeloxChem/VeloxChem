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


#include "SimdElectronRepulsionCtrVrrHS.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ctr_hs_electron_repulsion_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                                    const size_t pa, const size_t fs0, const size_t fs1,
                                    const size_t gs, const size_t ncols, const double alpha,
                                    const double beta, const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5625 * std::sqrt(14.0) / alpha;
    const auto f_1 = 0.5625 * std::sqrt(14.0) * beta / (alpha * p);
    const auto f_2 = 0.9375 * std::sqrt(14.0);
    const auto f_3 = 1.875 * std::sqrt(14.0);
    const auto f_4 = 0.1875 * std::sqrt(14.0);
    const auto f_5 = 1.5 * std::sqrt(35.0);
    const auto f_6 = 0.0625 * std::sqrt(70.0) / alpha;
    const auto f_7 = 0.5 * std::sqrt(70.0) / alpha;
    const auto f_8 = 0.0625 * std::sqrt(70.0) * beta / (alpha * p);
    const auto f_9 = 0.5 * std::sqrt(70.0) * beta / (alpha * p);
    const auto f_10 = 0.1875 * std::sqrt(70.0);
    const auto f_11 = 1.5 * std::sqrt(70.0);
    const auto f_12 = 0.125 * std::sqrt(70.0);
    const auto f_13 = 0.0625 * std::sqrt(70.0);
    const auto f_14 = 0.5 * std::sqrt(70.0);
    const auto f_15 = 0.5 * std::sqrt(105.0);
    const auto f_16 = std::sqrt(105.0);
    const auto f_17 = 0.375 * std::sqrt(15.0) / alpha;
    const auto f_18 = 1.5 * std::sqrt(15.0) / alpha;
    const auto f_19 = 0.375 * std::sqrt(15.0) * beta / (alpha * p);
    const auto f_20 = 1.5 * std::sqrt(15.0) * beta / (alpha * p);
    const auto f_21 = 0.125 * std::sqrt(15.0);
    const auto f_22 = 1.5 * std::sqrt(15.0);
    const auto f_23 = 0.25 * std::sqrt(15.0);
    const auto f_24 = std::sqrt(15.0);
    const auto f_25 = 3.0 / alpha;
    const auto f_26 = 3.0 * beta / (alpha * p);
    const auto f_27 = 0.25 * std::sqrt(15.0) / alpha;
    const auto f_28 = 0.25 * std::sqrt(15.0) * beta / (alpha * p);
    const auto f_29 = 0.25 * std::sqrt(105.0);
    const auto f_30 = 0.125 * std::sqrt(70.0) / alpha;
    const auto f_31 = 0.125 * std::sqrt(70.0) * beta / (alpha * p);
    const auto f_32 = 0.375 * std::sqrt(35.0);
    const auto f_33 = 2.25 * std::sqrt(35.0);
    const auto f_34 = 0.375 * std::sqrt(14.0) / alpha;
    const auto f_35 = 1.875 * std::sqrt(14.0) / alpha;
    const auto f_36 = 0.375 * std::sqrt(14.0) * beta / (alpha * p);
    const auto f_37 = 1.875 * std::sqrt(14.0) * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *fs0_0 = buffer.data(fs0 + 0);
    const auto *fs0_1 = buffer.data(fs0 + 1);
    const auto *fs0_2 = buffer.data(fs0 + 2);
    const auto *fs0_3 = buffer.data(fs0 + 3);
    const auto *fs0_4 = buffer.data(fs0 + 4);
    const auto *fs0_5 = buffer.data(fs0 + 5);

    const auto *fs1_0 = buffer.data(fs1 + 0);
    const auto *fs1_1 = buffer.data(fs1 + 1);
    const auto *fs1_2 = buffer.data(fs1 + 2);
    const auto *fs1_3 = buffer.data(fs1 + 3);
    const auto *fs1_4 = buffer.data(fs1 + 4);
    const auto *fs1_5 = buffer.data(fs1 + 5);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_10 = buffer.data(gs + 10);

#pragma omp simd aligned(pa_x, pa_y, fs0_3, fs0_4, fs1_3, fs1_4, gs_0, gs_1, gs_3, gs_4, gs_6, \
                         gs_7, gs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += -f_0 * fs0_3[k]
                  + f_1 * fs1_3[k]
                  + f_2 * pa_y[k] * gs_0[k]
                  - f_3 * pa_x[k] * gs_4[k]
                  + f_4 * pa_y[k] * gs_6[k];

        g_1[k] += f_5 * pa_y[k] * gs_1[k]
                  - f_5 * pa_x[k] * gs_7[k];

        g_2[k] += f_6 * fs0_3[k]
                  - f_7 * fs0_4[k]
                  - f_8 * fs1_3[k]
                  + f_9 * fs1_4[k]
                  - f_10 * pa_y[k] * gs_0[k]
                  + f_11 * pa_y[k] * gs_3[k]
                  - f_12 * pa_x[k] * gs_4[k]
                  + f_13 * pa_y[k] * gs_6[k]
                  - f_14 * pa_y[k] * gs_8[k];
    }

#pragma omp simd aligned(pa_x, pa_y, gs_1, gs_7, gs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += -f_15 * pa_y[k] * gs_1[k]
                  - f_15 * pa_x[k] * gs_7[k]
                  + f_16 * pa_x[k] * gs_9[k];
    }

#pragma omp simd aligned(pa_x, pa_y, fs0_3, fs0_4, fs1_3, fs1_4, gs_0, gs_3, gs_4, gs_6, gs_8, \
                         gs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += f_17 * fs0_3[k]
                  - f_18 * fs0_4[k]
                  - f_19 * fs1_3[k]
                  + f_20 * fs1_4[k]
                  + f_21 * pa_y[k] * gs_0[k]
                  - f_22 * pa_y[k] * gs_3[k]
                  + f_23 * pa_x[k] * gs_4[k]
                  + f_21 * pa_y[k] * gs_6[k]
                  - f_22 * pa_y[k] * gs_8[k]
                  + f_24 * pa_y[k] * gs_10[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, fs0_5, fs1_5, gs_0, gs_2, gs_5, gs_6, gs_9, \
                         gs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += -f_25 * fs0_5[k]
                  + f_26 * fs1_5[k]
                  + 1.875 * pa_z[k] * gs_0[k]
                  + 3.75 * pa_z[k] * gs_2[k]
                  - 5.0 * pa_x[k] * gs_5[k]
                  + 1.875 * pa_z[k] * gs_6[k]
                  - 5.0 * pa_y[k] * gs_9[k]
                  + pa_z[k] * gs_10[k];
    }

#pragma omp simd aligned(pa_x, fs0_0, fs0_1, fs0_2, fs1_0, fs1_1, fs1_2, gs_0, gs_2, gs_3, \
                         gs_6, gs_8, gs_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += f_27 * fs0_0[k]
                  + f_27 * fs0_1[k]
                  - f_18 * fs0_2[k]
                  - f_28 * fs1_0[k]
                  - f_28 * fs1_1[k]
                  + f_20 * fs1_2[k]
                  + f_21 * pa_x[k] * gs_0[k]
                  + f_23 * pa_x[k] * gs_2[k]
                  - f_22 * pa_x[k] * gs_3[k]
                  + f_21 * pa_x[k] * gs_6[k]
                  - f_22 * pa_x[k] * gs_8[k]
                  + f_24 * pa_x[k] * gs_10[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, gs_0, gs_5, gs_6, gs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_7[k] += -f_29 * pa_z[k] * gs_0[k]
                  + f_15 * pa_x[k] * gs_5[k]
                  + f_29 * pa_z[k] * gs_6[k]
                  - f_15 * pa_y[k] * gs_9[k];
    }

#pragma omp simd aligned(pa_x, pa_z, fs0_0, fs0_1, fs0_2, fs1_0, fs1_1, fs1_2, gs_0, gs_2, \
                         gs_3, gs_6, gs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += -f_30 * fs0_0[k]
                  + f_30 * fs0_1[k]
                  + f_7 * fs0_2[k]
                  + f_31 * fs1_0[k]
                  - f_31 * fs1_1[k]
                  - f_9 * fs1_2[k]
                  - f_13 * pa_x[k] * gs_0[k]
                  + f_12 * pa_x[k] * gs_2[k]
                  + f_14 * pa_x[k] * gs_3[k]
                  + f_10 * pa_x[k] * gs_6[k]
                  - f_11 * pa_x[k] * gs_8[k];

        g_9[k] += f_32 * pa_z[k] * gs_0[k]
                  - f_33 * pa_z[k] * gs_2[k]
                  + f_32 * pa_z[k] * gs_6[k];

        g_10[k] += f_34 * fs0_0[k]
                   - f_35 * fs0_1[k]
                   - f_36 * fs1_0[k]
                   + f_37 * fs1_1[k]
                   + f_4 * pa_x[k] * gs_0[k]
                   - f_3 * pa_x[k] * gs_2[k]
                   + f_2 * pa_x[k] * gs_6[k];
    }
}

}  // namespace simdt2ceri
