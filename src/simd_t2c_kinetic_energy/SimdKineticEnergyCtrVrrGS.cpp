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


#include "SimdKineticEnergyCtrVrrGS.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_ctr_gs_kinetic_energy_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                                const size_t pa, const size_t ds_s, const size_t ds,
                                const size_t fs, const size_t gs_s, const size_t ncols,
                                const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(35.0);
    const auto f_1 = std::sqrt(35.0) * alpha * beta / p;
    const auto f_2 = 0.75 * std::sqrt(70.0);
    const auto f_3 = 0.25 * std::sqrt(70.0);
    const auto f_4 = 1.5 * std::sqrt(70.0) * alpha * beta / p;
    const auto f_5 = 0.5 * std::sqrt(70.0) * alpha * beta / p;
    const auto f_6 = 0.5 * std::sqrt(5.0);
    const auto f_7 = 3.0 * std::sqrt(5.0);
    const auto f_8 = std::sqrt(5.0) * alpha * beta / p;
    const auto f_9 = 6.0 * std::sqrt(5.0) * alpha * beta / p;
    const auto f_10 = 0.75 * std::sqrt(10.0);
    const auto f_11 = std::sqrt(10.0);
    const auto f_12 = 1.5 * std::sqrt(10.0) * alpha * beta / p;
    const auto f_13 = 2.0 * std::sqrt(10.0) * alpha * beta / p;
    const auto f_14 = 1.125 * beta / p;
    const auto f_15 = 1.875 * beta / p;
    const auto f_16 = 3.0 * beta / p;
    const auto f_17 = 0.5625 / p;
    const auto f_18 = 0.9375 / p;
    const auto f_19 = 1.5 / p;
    const auto f_20 = 0.75 * alpha * beta / p;
    const auto f_21 = 1.5 * alpha * beta / p;
    const auto f_22 = 6.0 * alpha * beta / p;
    const auto f_23 = 2.0 * alpha * beta / p;
    const auto f_24 = 0.75 * std::sqrt(5.0) * beta / p;
    const auto f_25 = 0.375 * std::sqrt(5.0) / p;
    const auto f_26 = 0.25 * std::sqrt(5.0);
    const auto f_27 = 1.5 * std::sqrt(5.0);
    const auto f_28 = 0.5 * std::sqrt(5.0) * alpha * beta / p;
    const auto f_29 = 3.0 * std::sqrt(5.0) * alpha * beta / p;
    const auto f_30 = 0.375 * std::sqrt(35.0) * beta / p;
    const auto f_31 = 0.1875 * std::sqrt(35.0) / p;
    const auto f_32 = 0.125 * std::sqrt(35.0);
    const auto f_33 = 0.75 * std::sqrt(35.0);
    const auto f_34 = 0.25 * std::sqrt(35.0) * alpha * beta / p;
    const auto f_35 = 1.5 * std::sqrt(35.0) * alpha * beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ds_s_0 = buffer.data(ds_s + 0);
    const auto *ds_s_1 = buffer.data(ds_s + 1);
    const auto *ds_s_2 = buffer.data(ds_s + 2);

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_1 = buffer.data(ds + 1);
    const auto *ds_2 = buffer.data(ds + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_1 = buffer.data(gs_s + 1);
    const auto *gs_s_2 = buffer.data(gs_s + 2);
    const auto *gs_s_3 = buffer.data(gs_s + 3);
    const auto *gs_s_4 = buffer.data(gs_s + 4);
    const auto *gs_s_5 = buffer.data(gs_s + 5);
    const auto *gs_s_6 = buffer.data(gs_s + 6);
    const auto *gs_s_7 = buffer.data(gs_s + 7);
    const auto *gs_s_8 = buffer.data(gs_s + 8);
    const auto *gs_s_9 = buffer.data(gs_s + 9);
    const auto *gs_s_10 = buffer.data(gs_s + 10);
    const auto *gs_s_11 = buffer.data(gs_s + 11);
    const auto *gs_s_12 = buffer.data(gs_s + 12);
    const auto *gs_s_13 = buffer.data(gs_s + 13);
    const auto *gs_s_14 = buffer.data(gs_s + 14);

#pragma omp simd aligned(pa_x, pa_y, pa_z, fs_0, fs_1, fs_4, fs_6, gs_s_1, gs_s_4, gs_s_6, \
                         gs_s_8, gs_s_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += f_0 * pa_y[k] * fs_0[k]
                  - f_0 * pa_x[k] * fs_4[k]
                  + f_1 * gs_s_1[k]
                  - f_1 * gs_s_6[k];

        g_1[k] += f_2 * pa_y[k] * fs_1[k]
                  - f_3 * pa_z[k] * fs_4[k]
                  + f_4 * gs_s_4[k]
                  - f_5 * gs_s_11[k];

        g_2[k] += -f_6 * pa_y[k] * fs_0[k]
                  - f_6 * pa_x[k] * fs_4[k]
                  + f_7 * pa_x[k] * fs_6[k]
                  - f_8 * gs_s_1[k]
                  - f_8 * gs_s_6[k]
                  + f_9 * gs_s_8[k];
    }

#pragma omp simd aligned(pa_y, pa_z, fs_1, fs_4, fs_7, gs_s_4, gs_s_11, \
                         gs_s_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += -f_10 * pa_y[k] * fs_1[k]
                  - f_10 * pa_z[k] * fs_4[k]
                  + f_11 * pa_y[k] * fs_7[k]
                  - f_12 * gs_s_4[k]
                  - f_12 * gs_s_11[k]
                  + f_13 * gs_s_13[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, ds_s_0, ds_s_1, ds_s_2, ds_0, ds_1, ds_2, fs_0, \
                         fs_2, fs_3, fs_4, fs_6, fs_7, gs_s_0, gs_s_3, gs_s_5, gs_s_10, \
                         gs_s_12, gs_s_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += -f_14 * ds_s_0[k]
                  - f_15 * ds_s_1[k]
                  + f_16 * ds_s_2[k]
                  + f_17 * ds_0[k]
                  + f_18 * ds_1[k]
                  - f_19 * ds_2[k]
                  + 0.375 * pa_x[k] * fs_0[k]
                  + 0.75 * pa_x[k] * fs_2[k]
                  - 3.0 * pa_x[k] * fs_3[k]
                  + 0.375 * pa_y[k] * fs_4[k]
                  - 3.0 * pa_y[k] * fs_6[k]
                  + pa_z[k] * fs_7[k]
                  + f_20 * gs_s_0[k]
                  + f_21 * gs_s_3[k]
                  - f_22 * gs_s_5[k]
                  + f_20 * gs_s_10[k]
                  - f_22 * gs_s_12[k]
                  + f_23 * gs_s_14[k];
    }

#pragma omp simd aligned(pa_x, pa_z, fs_0, fs_5, fs_7, gs_s_2, gs_s_7, \
                         gs_s_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += -f_10 * pa_z[k] * fs_0[k]
                  - f_10 * pa_x[k] * fs_5[k]
                  + f_11 * pa_x[k] * fs_7[k]
                  - f_12 * gs_s_2[k]
                  - f_12 * gs_s_7[k]
                  + f_13 * gs_s_9[k];
    }

#pragma omp simd aligned(pa_x, pa_y, ds_s_0, ds_s_1, ds_0, ds_1, fs_0, fs_3, fs_4, fs_6, \
                         gs_s_0, gs_s_5, gs_s_10, gs_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += f_24 * ds_s_0[k]
                  - f_24 * ds_s_1[k]
                  - f_25 * ds_0[k]
                  + f_25 * ds_1[k]
                  - f_26 * pa_x[k] * fs_0[k]
                  + f_27 * pa_x[k] * fs_3[k]
                  + f_26 * pa_y[k] * fs_4[k]
                  - f_27 * pa_y[k] * fs_6[k]
                  - f_28 * gs_s_0[k]
                  + f_29 * gs_s_5[k]
                  + f_28 * gs_s_10[k]
                  - f_29 * gs_s_12[k];
    }

#pragma omp simd aligned(pa_x, pa_z, fs_0, fs_5, gs_s_2, gs_s_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_7[k] += f_3 * pa_z[k] * fs_0[k]
                  - f_2 * pa_x[k] * fs_5[k]
                  + f_5 * gs_s_2[k]
                  - f_4 * gs_s_7[k];
    }

#pragma omp simd aligned(pa_x, pa_y, ds_s_0, ds_s_1, ds_0, ds_1, fs_0, fs_2, fs_4, gs_s_0, \
                         gs_s_3, gs_s_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += -f_30 * ds_s_0[k]
                  + f_30 * ds_s_1[k]
                  + f_31 * ds_0[k]
                  - f_31 * ds_1[k]
                  + f_32 * pa_x[k] * fs_0[k]
                  - f_33 * pa_x[k] * fs_2[k]
                  + f_32 * pa_y[k] * fs_4[k]
                  + f_34 * gs_s_0[k]
                  - f_35 * gs_s_3[k]
                  + f_34 * gs_s_10[k];
    }
}

}  // namespace simdkin
