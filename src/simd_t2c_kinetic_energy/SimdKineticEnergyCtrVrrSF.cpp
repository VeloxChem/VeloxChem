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


#include "SimdKineticEnergyCtrVrrSF.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_ctr_sf_kinetic_energy_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                                const size_t pb, const size_t sp_s, const size_t sf_s,
                                const size_t sp, const size_t sd, const size_t ncols,
                                const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(10.0) * alpha / p;
    const auto f_1 = 1.5 * std::sqrt(10.0) * alpha * beta / p;
    const auto f_2 = 0.5 * std::sqrt(10.0) * alpha * beta / p;
    const auto f_3 = 0.25 * std::sqrt(10.0) / p;
    const auto f_4 = 0.75 * std::sqrt(10.0);
    const auto f_5 = 0.25 * std::sqrt(10.0);
    const auto f_6 = 2.0 * std::sqrt(15.0) * alpha * beta / p;
    const auto f_7 = std::sqrt(15.0);
    const auto f_8 = 0.5 * std::sqrt(6.0) * alpha / p;
    const auto f_9 = 0.5 * std::sqrt(6.0) * alpha * beta / p;
    const auto f_10 = 2.0 * std::sqrt(6.0) * alpha * beta / p;
    const auto f_11 = 0.25 * std::sqrt(6.0) / p;
    const auto f_12 = 0.25 * std::sqrt(6.0);
    const auto f_13 = std::sqrt(6.0);
    const auto f_14 = 2.0 * alpha / p;
    const auto f_15 = 3.0 * alpha * beta / p;
    const auto f_16 = 2.0 * alpha * beta / p;
    const auto f_17 = 1.0 / p;
    const auto f_18 = std::sqrt(15.0) * alpha * beta / p;
    const auto f_19 = 0.5 * std::sqrt(15.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sp_s_0 = buffer.data(sp_s + 0);
    const auto *sp_s_1 = buffer.data(sp_s + 1);
    const auto *sp_s_2 = buffer.data(sp_s + 2);

    const auto *sf_s_0 = buffer.data(sf_s + 0);
    const auto *sf_s_1 = buffer.data(sf_s + 1);
    const auto *sf_s_2 = buffer.data(sf_s + 2);
    const auto *sf_s_3 = buffer.data(sf_s + 3);
    const auto *sf_s_4 = buffer.data(sf_s + 4);
    const auto *sf_s_5 = buffer.data(sf_s + 5);
    const auto *sf_s_6 = buffer.data(sf_s + 6);
    const auto *sf_s_7 = buffer.data(sf_s + 7);
    const auto *sf_s_8 = buffer.data(sf_s + 8);
    const auto *sf_s_9 = buffer.data(sf_s + 9);

    const auto *sp_0 = buffer.data(sp + 0);
    const auto *sp_1 = buffer.data(sp + 1);
    const auto *sp_2 = buffer.data(sp + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);
    const auto *sd_3 = buffer.data(sd + 3);

#pragma omp simd aligned(pb_x, pb_y, sp_s_1, sf_s_1, sf_s_4, sf_s_6, sf_s_8, sp_1, sd_0, sd_1, \
                         sd_2, sd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += f_0 * sp_s_1[k]
                  + f_1 * sf_s_1[k]
                  - f_2 * sf_s_6[k]
                  - f_3 * sp_1[k]
                  + f_4 * pb_y[k] * sd_0[k]
                  - f_5 * pb_y[k] * sd_1[k];

        g_1[k] += f_6 * sf_s_4[k]
                  + f_7 * pb_x[k] * sd_2[k];

        g_2[k] += f_8 * sp_s_1[k]
                  - f_9 * sf_s_1[k]
                  - f_9 * sf_s_6[k]
                  + f_10 * sf_s_8[k]
                  - f_11 * sp_1[k]
                  - f_12 * pb_y[k] * sd_0[k]
                  - f_12 * pb_y[k] * sd_1[k]
                  + f_13 * pb_y[k] * sd_3[k];
    }

#pragma omp simd aligned(pb_z, sp_s_2, sf_s_2, sf_s_7, sf_s_9, sp_2, sd_0, sd_1, \
                         sd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += -f_14 * sp_s_2[k]
                  - f_15 * sf_s_2[k]
                  - f_15 * sf_s_7[k]
                  + f_16 * sf_s_9[k]
                  + f_17 * sp_2[k]
                  - 1.5 * pb_z[k] * sd_0[k]
                  - 1.5 * pb_z[k] * sd_1[k]
                  + pb_z[k] * sd_3[k];
    }

#pragma omp simd aligned(pb_x, pb_z, sp_s_0, sf_s_0, sf_s_2, sf_s_3, sf_s_5, sf_s_7, sp_0, \
                         sd_0, sd_1, sd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += f_8 * sp_s_0[k]
                  - f_9 * sf_s_0[k]
                  - f_9 * sf_s_3[k]
                  + f_10 * sf_s_5[k]
                  - f_11 * sp_0[k]
                  - f_12 * pb_x[k] * sd_0[k]
                  - f_12 * pb_x[k] * sd_1[k]
                  + f_13 * pb_x[k] * sd_3[k];

        g_5[k] += f_18 * sf_s_2[k]
                  - f_18 * sf_s_7[k]
                  + f_19 * pb_z[k] * sd_0[k]
                  - f_19 * pb_z[k] * sd_1[k];

        g_6[k] += -f_0 * sp_s_0[k]
                  + f_2 * sf_s_0[k]
                  - f_1 * sf_s_3[k]
                  + f_3 * sp_0[k]
                  + f_5 * pb_x[k] * sd_0[k]
                  - f_4 * pb_x[k] * sd_1[k];
    }
}

}  // namespace simdkin
