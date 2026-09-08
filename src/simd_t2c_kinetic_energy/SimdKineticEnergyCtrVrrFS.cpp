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


#include "SimdKineticEnergyCtrVrrFS.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_ctr_fs_kinetic_energy_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                                const size_t pa, const size_t ps_s, const size_t ps,
                                const size_t ds, const size_t fs_s, const size_t ncols,
                                const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(10.0) * beta / p;
    const auto f_1 = 0.25 * std::sqrt(10.0) / p;
    const auto f_2 = 0.75 * std::sqrt(10.0);
    const auto f_3 = 0.25 * std::sqrt(10.0);
    const auto f_4 = 1.5 * std::sqrt(10.0) * alpha * beta / p;
    const auto f_5 = 0.5 * std::sqrt(10.0) * alpha * beta / p;
    const auto f_6 = std::sqrt(15.0);
    const auto f_7 = 2.0 * std::sqrt(15.0) * alpha * beta / p;
    const auto f_8 = 0.5 * std::sqrt(6.0) * beta / p;
    const auto f_9 = 0.25 * std::sqrt(6.0) / p;
    const auto f_10 = 0.25 * std::sqrt(6.0);
    const auto f_11 = std::sqrt(6.0);
    const auto f_12 = 0.5 * std::sqrt(6.0) * alpha * beta / p;
    const auto f_13 = 2.0 * std::sqrt(6.0) * alpha * beta / p;
    const auto f_14 = 2.0 * beta / p;
    const auto f_15 = 1.0 / p;
    const auto f_16 = 3.0 * alpha * beta / p;
    const auto f_17 = 2.0 * alpha * beta / p;
    const auto f_18 = 0.5 * std::sqrt(15.0);
    const auto f_19 = std::sqrt(15.0) * alpha * beta / p;

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ps_s_0 = buffer.data(ps_s + 0);
    const auto *ps_s_1 = buffer.data(ps_s + 1);
    const auto *ps_s_2 = buffer.data(ps_s + 2);

    const auto *ps_0 = buffer.data(ps + 0);
    const auto *ps_1 = buffer.data(ps + 1);
    const auto *ps_2 = buffer.data(ps + 2);

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_1 = buffer.data(ds + 1);
    const auto *ds_2 = buffer.data(ds + 2);
    const auto *ds_3 = buffer.data(ds + 3);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_1 = buffer.data(fs_s + 1);
    const auto *fs_s_2 = buffer.data(fs_s + 2);
    const auto *fs_s_3 = buffer.data(fs_s + 3);
    const auto *fs_s_4 = buffer.data(fs_s + 4);
    const auto *fs_s_5 = buffer.data(fs_s + 5);
    const auto *fs_s_6 = buffer.data(fs_s + 6);
    const auto *fs_s_7 = buffer.data(fs_s + 7);
    const auto *fs_s_8 = buffer.data(fs_s + 8);
    const auto *fs_s_9 = buffer.data(fs_s + 9);

#pragma omp simd aligned(pa_x, pa_y, ps_s_1, ps_1, ds_0, ds_1, ds_2, ds_3, fs_s_1, fs_s_4, \
                         fs_s_6, fs_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += f_0 * ps_s_1[k]
                  - f_1 * ps_1[k]
                  + f_2 * pa_y[k] * ds_0[k]
                  - f_3 * pa_y[k] * ds_1[k]
                  + f_4 * fs_s_1[k]
                  - f_5 * fs_s_6[k];

        g_1[k] += f_6 * pa_x[k] * ds_2[k]
                  + f_7 * fs_s_4[k];

        g_2[k] += f_8 * ps_s_1[k]
                  - f_9 * ps_1[k]
                  - f_10 * pa_y[k] * ds_0[k]
                  - f_10 * pa_y[k] * ds_1[k]
                  + f_11 * pa_y[k] * ds_3[k]
                  - f_12 * fs_s_1[k]
                  - f_12 * fs_s_6[k]
                  + f_13 * fs_s_8[k];
    }

#pragma omp simd aligned(pa_z, ps_s_2, ps_2, ds_0, ds_1, ds_3, fs_s_2, fs_s_7, \
                         fs_s_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += -f_14 * ps_s_2[k]
                  + f_15 * ps_2[k]
                  - 1.5 * pa_z[k] * ds_0[k]
                  - 1.5 * pa_z[k] * ds_1[k]
                  + pa_z[k] * ds_3[k]
                  - f_16 * fs_s_2[k]
                  - f_16 * fs_s_7[k]
                  + f_17 * fs_s_9[k];
    }

#pragma omp simd aligned(pa_x, pa_z, ps_s_0, ps_0, ds_0, ds_1, ds_3, fs_s_0, fs_s_2, fs_s_3, \
                         fs_s_5, fs_s_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += f_8 * ps_s_0[k]
                  - f_9 * ps_0[k]
                  - f_10 * pa_x[k] * ds_0[k]
                  - f_10 * pa_x[k] * ds_1[k]
                  + f_11 * pa_x[k] * ds_3[k]
                  - f_12 * fs_s_0[k]
                  - f_12 * fs_s_3[k]
                  + f_13 * fs_s_5[k];

        g_5[k] += f_18 * pa_z[k] * ds_0[k]
                  - f_18 * pa_z[k] * ds_1[k]
                  + f_19 * fs_s_2[k]
                  - f_19 * fs_s_7[k];

        g_6[k] += -f_0 * ps_s_0[k]
                  + f_1 * ps_0[k]
                  + f_3 * pa_x[k] * ds_0[k]
                  - f_2 * pa_x[k] * ds_1[k]
                  + f_5 * fs_s_0[k]
                  - f_4 * fs_s_3[k];
    }
}

}  // namespace simdkin
