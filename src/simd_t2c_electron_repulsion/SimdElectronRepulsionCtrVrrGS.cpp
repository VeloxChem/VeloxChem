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


#include "SimdElectronRepulsionCtrVrrGS.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ctr_gs_electron_repulsion_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                                    const size_t pa, const size_t ds0, const size_t ds1,
                                    const size_t fs, const size_t ncols, const double alpha,
                                    const double beta, const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(35.0);
    const auto f_1 = 0.75 * std::sqrt(70.0);
    const auto f_2 = 0.25 * std::sqrt(70.0);
    const auto f_3 = 0.5 * std::sqrt(5.0);
    const auto f_4 = 3.0 * std::sqrt(5.0);
    const auto f_5 = 0.75 * std::sqrt(10.0);
    const auto f_6 = std::sqrt(10.0);
    const auto f_7 = 0.5625 / alpha;
    const auto f_8 = 0.9375 / alpha;
    const auto f_9 = 1.5 / alpha;
    const auto f_10 = 0.5625 * beta / (alpha * p);
    const auto f_11 = 0.9375 * beta / (alpha * p);
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 0.375 * std::sqrt(5.0) / alpha;
    const auto f_14 = 0.375 * std::sqrt(5.0) * beta / (alpha * p);
    const auto f_15 = 0.25 * std::sqrt(5.0);
    const auto f_16 = 1.5 * std::sqrt(5.0);
    const auto f_17 = 0.1875 * std::sqrt(35.0) / alpha;
    const auto f_18 = 0.1875 * std::sqrt(35.0) * beta / (alpha * p);
    const auto f_19 = 0.125 * std::sqrt(35.0);
    const auto f_20 = 0.75 * std::sqrt(35.0);

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

    const auto *ds0_0 = buffer.data(ds0 + 0);
    const auto *ds0_1 = buffer.data(ds0 + 1);
    const auto *ds0_2 = buffer.data(ds0 + 2);

    const auto *ds1_0 = buffer.data(ds1 + 0);
    const auto *ds1_1 = buffer.data(ds1 + 1);
    const auto *ds1_2 = buffer.data(ds1 + 2);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);
    const auto *fs_6 = buffer.data(fs + 6);
    const auto *fs_7 = buffer.data(fs + 7);

#pragma omp simd aligned(pa_x, pa_y, pa_z, fs_0, fs_1, fs_4, fs_6, \
                         fs_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += f_0 * pa_y[k] * fs_0[k]
                  - f_0 * pa_x[k] * fs_4[k];

        g_1[k] += f_1 * pa_y[k] * fs_1[k]
                  - f_2 * pa_z[k] * fs_4[k];

        g_2[k] += -f_3 * pa_y[k] * fs_0[k]
                  - f_3 * pa_x[k] * fs_4[k]
                  + f_4 * pa_x[k] * fs_6[k];

        g_3[k] += -f_5 * pa_y[k] * fs_1[k]
                  - f_5 * pa_z[k] * fs_4[k]
                  + f_6 * pa_y[k] * fs_7[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, ds0_0, ds0_1, ds0_2, ds1_0, ds1_1, ds1_2, fs_0, \
                         fs_2, fs_3, fs_4, fs_6, fs_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += f_7 * ds0_0[k]
                  + f_8 * ds0_1[k]
                  - f_9 * ds0_2[k]
                  - f_10 * ds1_0[k]
                  - f_11 * ds1_1[k]
                  + f_12 * ds1_2[k]
                  + 0.375 * pa_x[k] * fs_0[k]
                  + 0.75 * pa_x[k] * fs_2[k]
                  - 3.0 * pa_x[k] * fs_3[k]
                  + 0.375 * pa_y[k] * fs_4[k]
                  - 3.0 * pa_y[k] * fs_6[k]
                  + pa_z[k] * fs_7[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, ds0_0, ds0_1, ds1_0, ds1_1, fs_0, fs_3, fs_4, fs_5, \
                         fs_6, fs_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += -f_5 * pa_z[k] * fs_0[k]
                  - f_5 * pa_x[k] * fs_5[k]
                  + f_6 * pa_x[k] * fs_7[k];

        g_6[k] += -f_13 * ds0_0[k]
                  + f_13 * ds0_1[k]
                  + f_14 * ds1_0[k]
                  - f_14 * ds1_1[k]
                  - f_15 * pa_x[k] * fs_0[k]
                  + f_16 * pa_x[k] * fs_3[k]
                  + f_15 * pa_y[k] * fs_4[k]
                  - f_16 * pa_y[k] * fs_6[k];

        g_7[k] += f_2 * pa_z[k] * fs_0[k]
                  - f_1 * pa_x[k] * fs_5[k];
    }

#pragma omp simd aligned(pa_x, pa_y, ds0_0, ds0_1, ds1_0, ds1_1, fs_0, fs_2, \
                         fs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += f_17 * ds0_0[k]
                  - f_17 * ds0_1[k]
                  - f_18 * ds1_0[k]
                  + f_18 * ds1_1[k]
                  + f_19 * pa_x[k] * fs_0[k]
                  - f_20 * pa_x[k] * fs_2[k]
                  + f_19 * pa_y[k] * fs_4[k];
    }
}

}  // namespace simdt2ceri
