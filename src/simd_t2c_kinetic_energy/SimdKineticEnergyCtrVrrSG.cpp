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


#include "SimdKineticEnergyCtrVrrSG.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_ctr_sg_kinetic_energy_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                                const size_t pb, const size_t sd_s, const size_t sg_s,
                                const size_t sd, const size_t sf, const size_t ncols,
                                const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = std::sqrt(35.0) * alpha * beta / p;
    const auto f_1 = 0.5 * std::sqrt(35.0);
    const auto f_2 = 1.5 * std::sqrt(70.0) * alpha * beta / p;
    const auto f_3 = 0.5 * std::sqrt(70.0) * alpha * beta / p;
    const auto f_4 = 0.75 * std::sqrt(70.0);
    const auto f_5 = 0.25 * std::sqrt(70.0);
    const auto f_6 = std::sqrt(5.0) * alpha * beta / p;
    const auto f_7 = 6.0 * std::sqrt(5.0) * alpha * beta / p;
    const auto f_8 = 0.5 * std::sqrt(5.0);
    const auto f_9 = 3.0 * std::sqrt(5.0);
    const auto f_10 = 1.5 * std::sqrt(10.0) * alpha * beta / p;
    const auto f_11 = 2.0 * std::sqrt(10.0) * alpha * beta / p;
    const auto f_12 = 0.75 * std::sqrt(10.0);
    const auto f_13 = std::sqrt(10.0);
    const auto f_14 = 1.125 * alpha / p;
    const auto f_15 = 1.875 * alpha / p;
    const auto f_16 = 3.0 * alpha / p;
    const auto f_17 = 0.75 * alpha * beta / p;
    const auto f_18 = 1.5 * alpha * beta / p;
    const auto f_19 = 6.0 * alpha * beta / p;
    const auto f_20 = 2.0 * alpha * beta / p;
    const auto f_21 = 0.5625 / p;
    const auto f_22 = 0.9375 / p;
    const auto f_23 = 1.5 / p;
    const auto f_24 = 0.75 * std::sqrt(5.0) * alpha / p;
    const auto f_25 = 0.5 * std::sqrt(5.0) * alpha * beta / p;
    const auto f_26 = 3.0 * std::sqrt(5.0) * alpha * beta / p;
    const auto f_27 = 0.375 * std::sqrt(5.0) / p;
    const auto f_28 = 0.25 * std::sqrt(5.0);
    const auto f_29 = 1.5 * std::sqrt(5.0);
    const auto f_30 = 0.375 * std::sqrt(35.0) * alpha / p;
    const auto f_31 = 0.25 * std::sqrt(35.0) * alpha * beta / p;
    const auto f_32 = 1.5 * std::sqrt(35.0) * alpha * beta / p;
    const auto f_33 = 0.1875 * std::sqrt(35.0) / p;
    const auto f_34 = 0.125 * std::sqrt(35.0);
    const auto f_35 = 0.75 * std::sqrt(35.0);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd_s_0 = buffer.data(sd_s + 0);
    const auto *sd_s_1 = buffer.data(sd_s + 1);
    const auto *sd_s_2 = buffer.data(sd_s + 2);

    const auto *sg_s_0 = buffer.data(sg_s + 0);
    const auto *sg_s_1 = buffer.data(sg_s + 1);
    const auto *sg_s_2 = buffer.data(sg_s + 2);
    const auto *sg_s_3 = buffer.data(sg_s + 3);
    const auto *sg_s_4 = buffer.data(sg_s + 4);
    const auto *sg_s_5 = buffer.data(sg_s + 5);
    const auto *sg_s_6 = buffer.data(sg_s + 6);
    const auto *sg_s_7 = buffer.data(sg_s + 7);
    const auto *sg_s_8 = buffer.data(sg_s + 8);
    const auto *sg_s_9 = buffer.data(sg_s + 9);
    const auto *sg_s_10 = buffer.data(sg_s + 10);
    const auto *sg_s_11 = buffer.data(sg_s + 11);
    const auto *sg_s_12 = buffer.data(sg_s + 12);
    const auto *sg_s_13 = buffer.data(sg_s + 13);
    const auto *sg_s_14 = buffer.data(sg_s + 14);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_1 = buffer.data(sd + 1);
    const auto *sd_2 = buffer.data(sd + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_6 = buffer.data(sf + 6);
    const auto *sf_7 = buffer.data(sf + 7);

#pragma omp simd aligned(pb_x, pb_y, pb_z, sg_s_1, sg_s_4, sg_s_6, sg_s_8, sg_s_11, sf_0, \
                         sf_1, sf_4, sf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += f_0 * sg_s_1[k]
                  - f_0 * sg_s_6[k]
                  + f_1 * pb_y[k] * sf_0[k]
                  - f_1 * pb_x[k] * sf_4[k];

        g_1[k] += f_2 * sg_s_4[k]
                  - f_3 * sg_s_11[k]
                  + f_4 * pb_y[k] * sf_1[k]
                  - f_5 * pb_z[k] * sf_4[k];

        g_2[k] += -f_6 * sg_s_1[k]
                  - f_6 * sg_s_6[k]
                  + f_7 * sg_s_8[k]
                  - f_8 * pb_y[k] * sf_0[k]
                  - f_8 * pb_x[k] * sf_4[k]
                  + f_9 * pb_x[k] * sf_6[k];
    }

#pragma omp simd aligned(pb_y, pb_z, sg_s_4, sg_s_11, sg_s_13, sf_1, sf_4, \
                         sf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += -f_10 * sg_s_4[k]
                  - f_10 * sg_s_11[k]
                  + f_11 * sg_s_13[k]
                  - f_12 * pb_y[k] * sf_1[k]
                  - f_12 * pb_z[k] * sf_4[k]
                  + f_13 * pb_y[k] * sf_7[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, sd_s_0, sd_s_1, sd_s_2, sg_s_0, sg_s_3, sg_s_5, \
                         sg_s_10, sg_s_12, sg_s_14, sd_0, sd_1, sd_2, sf_0, sf_2, sf_3, sf_4, \
                         sf_6, sf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += -f_14 * sd_s_0[k]
                  - f_15 * sd_s_1[k]
                  + f_16 * sd_s_2[k]
                  + f_17 * sg_s_0[k]
                  + f_18 * sg_s_3[k]
                  - f_19 * sg_s_5[k]
                  + f_17 * sg_s_10[k]
                  - f_19 * sg_s_12[k]
                  + f_20 * sg_s_14[k]
                  + f_21 * sd_0[k]
                  + f_22 * sd_1[k]
                  - f_23 * sd_2[k]
                  + 0.375 * pb_x[k] * sf_0[k]
                  + 0.75 * pb_x[k] * sf_2[k]
                  - 3.0 * pb_x[k] * sf_3[k]
                  + 0.375 * pb_y[k] * sf_4[k]
                  - 3.0 * pb_y[k] * sf_6[k]
                  + pb_z[k] * sf_7[k];
    }

#pragma omp simd aligned(pb_x, pb_z, sg_s_2, sg_s_7, sg_s_9, sf_0, sf_5, \
                         sf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += -f_10 * sg_s_2[k]
                  - f_10 * sg_s_7[k]
                  + f_11 * sg_s_9[k]
                  - f_12 * pb_z[k] * sf_0[k]
                  - f_12 * pb_x[k] * sf_5[k]
                  + f_13 * pb_x[k] * sf_7[k];
    }

#pragma omp simd aligned(pb_x, pb_y, sd_s_0, sd_s_1, sg_s_0, sg_s_5, sg_s_10, sg_s_12, sd_0, \
                         sd_1, sf_0, sf_3, sf_4, sf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += f_24 * sd_s_0[k]
                  - f_24 * sd_s_1[k]
                  - f_25 * sg_s_0[k]
                  + f_26 * sg_s_5[k]
                  + f_25 * sg_s_10[k]
                  - f_26 * sg_s_12[k]
                  - f_27 * sd_0[k]
                  + f_27 * sd_1[k]
                  - f_28 * pb_x[k] * sf_0[k]
                  + f_29 * pb_x[k] * sf_3[k]
                  + f_28 * pb_y[k] * sf_4[k]
                  - f_29 * pb_y[k] * sf_6[k];
    }

#pragma omp simd aligned(pb_x, pb_z, sg_s_2, sg_s_7, sf_0, sf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_7[k] += f_3 * sg_s_2[k]
                  - f_2 * sg_s_7[k]
                  + f_5 * pb_z[k] * sf_0[k]
                  - f_4 * pb_x[k] * sf_5[k];
    }

#pragma omp simd aligned(pb_x, pb_y, sd_s_0, sd_s_1, sg_s_0, sg_s_3, sg_s_10, sd_0, sd_1, \
                         sf_0, sf_2, sf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += -f_30 * sd_s_0[k]
                  + f_30 * sd_s_1[k]
                  + f_31 * sg_s_0[k]
                  - f_32 * sg_s_3[k]
                  + f_31 * sg_s_10[k]
                  + f_33 * sd_0[k]
                  - f_33 * sd_1[k]
                  + f_34 * pb_x[k] * sf_0[k]
                  - f_35 * pb_x[k] * sf_2[k]
                  + f_34 * pb_y[k] * sf_4[k];
    }
}

}  // namespace simdkin
