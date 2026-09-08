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


#include "SimdOverlapCtrVrrSH.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_ctr_sh_overlap_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                         const size_t pb, const size_t sf, const size_t sg, const size_t ncols,
                         const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5625 * std::sqrt(14.0) / p;
    const auto f_1 = 0.9375 * std::sqrt(14.0);
    const auto f_2 = 1.875 * std::sqrt(14.0);
    const auto f_3 = 0.1875 * std::sqrt(14.0);
    const auto f_4 = 1.5 * std::sqrt(35.0);
    const auto f_5 = 0.0625 * std::sqrt(70.0) / p;
    const auto f_6 = 0.5 * std::sqrt(70.0) / p;
    const auto f_7 = 0.1875 * std::sqrt(70.0);
    const auto f_8 = 1.5 * std::sqrt(70.0);
    const auto f_9 = 0.125 * std::sqrt(70.0);
    const auto f_10 = 0.0625 * std::sqrt(70.0);
    const auto f_11 = 0.5 * std::sqrt(70.0);
    const auto f_12 = 0.5 * std::sqrt(105.0);
    const auto f_13 = std::sqrt(105.0);
    const auto f_14 = 0.375 * std::sqrt(15.0) / p;
    const auto f_15 = 1.5 * std::sqrt(15.0) / p;
    const auto f_16 = 0.125 * std::sqrt(15.0);
    const auto f_17 = 1.5 * std::sqrt(15.0);
    const auto f_18 = 0.25 * std::sqrt(15.0);
    const auto f_19 = std::sqrt(15.0);
    const auto f_20 = 3.0 / p;
    const auto f_21 = 0.25 * std::sqrt(15.0) / p;
    const auto f_22 = 0.25 * std::sqrt(105.0);
    const auto f_23 = 0.125 * std::sqrt(70.0) / p;
    const auto f_24 = 0.375 * std::sqrt(35.0);
    const auto f_25 = 2.25 * std::sqrt(35.0);
    const auto f_26 = 0.375 * std::sqrt(14.0) / p;
    const auto f_27 = 1.875 * std::sqrt(14.0) / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_4 = buffer.data(sf + 4);
    const auto *sf_5 = buffer.data(sf + 5);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);

#pragma omp simd aligned(pb_x, pb_y, sf_3, sf_4, sg_0, sg_1, sg_3, sg_4, sg_6, sg_7, sg_8, \
                         sg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += -f_0 * sf_3[k]
                  + f_1 * pb_y[k] * sg_0[k]
                  - f_2 * pb_x[k] * sg_4[k]
                  + f_3 * pb_y[k] * sg_6[k];

        g_1[k] += f_4 * pb_y[k] * sg_1[k]
                  - f_4 * pb_x[k] * sg_7[k];

        g_2[k] += f_5 * sf_3[k]
                  - f_6 * sf_4[k]
                  - f_7 * pb_y[k] * sg_0[k]
                  + f_8 * pb_y[k] * sg_3[k]
                  - f_9 * pb_x[k] * sg_4[k]
                  + f_10 * pb_y[k] * sg_6[k]
                  - f_11 * pb_y[k] * sg_8[k];

        g_3[k] += -f_12 * pb_y[k] * sg_1[k]
                  - f_12 * pb_x[k] * sg_7[k]
                  + f_13 * pb_x[k] * sg_9[k];
    }

#pragma omp simd aligned(pb_x, pb_y, sf_3, sf_4, sg_0, sg_3, sg_4, sg_6, sg_8, \
                         sg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += f_14 * sf_3[k]
                  - f_15 * sf_4[k]
                  + f_16 * pb_y[k] * sg_0[k]
                  - f_17 * pb_y[k] * sg_3[k]
                  + f_18 * pb_x[k] * sg_4[k]
                  + f_16 * pb_y[k] * sg_6[k]
                  - f_17 * pb_y[k] * sg_8[k]
                  + f_19 * pb_y[k] * sg_10[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, sf_5, sg_0, sg_2, sg_5, sg_6, sg_9, \
                         sg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += -f_20 * sf_5[k]
                  + 1.875 * pb_z[k] * sg_0[k]
                  + 3.75 * pb_z[k] * sg_2[k]
                  - 5.0 * pb_x[k] * sg_5[k]
                  + 1.875 * pb_z[k] * sg_6[k]
                  - 5.0 * pb_y[k] * sg_9[k]
                  + pb_z[k] * sg_10[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, sf_0, sf_1, sf_2, sg_0, sg_2, sg_3, sg_5, sg_6, \
                         sg_8, sg_9, sg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += f_21 * sf_0[k]
                  + f_21 * sf_1[k]
                  - f_15 * sf_2[k]
                  + f_16 * pb_x[k] * sg_0[k]
                  + f_18 * pb_x[k] * sg_2[k]
                  - f_17 * pb_x[k] * sg_3[k]
                  + f_16 * pb_x[k] * sg_6[k]
                  - f_17 * pb_x[k] * sg_8[k]
                  + f_19 * pb_x[k] * sg_10[k];

        g_7[k] += -f_22 * pb_z[k] * sg_0[k]
                  + f_12 * pb_x[k] * sg_5[k]
                  + f_22 * pb_z[k] * sg_6[k]
                  - f_12 * pb_y[k] * sg_9[k];
    }

#pragma omp simd aligned(pb_x, pb_z, sf_0, sf_1, sf_2, sg_0, sg_2, sg_3, sg_6, \
                         sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += -f_23 * sf_0[k]
                  + f_23 * sf_1[k]
                  + f_6 * sf_2[k]
                  - f_10 * pb_x[k] * sg_0[k]
                  + f_9 * pb_x[k] * sg_2[k]
                  + f_11 * pb_x[k] * sg_3[k]
                  + f_7 * pb_x[k] * sg_6[k]
                  - f_8 * pb_x[k] * sg_8[k];

        g_9[k] += f_24 * pb_z[k] * sg_0[k]
                  - f_25 * pb_z[k] * sg_2[k]
                  + f_24 * pb_z[k] * sg_6[k];

        g_10[k] += f_26 * sf_0[k]
                   - f_27 * sf_1[k]
                   + f_3 * pb_x[k] * sg_0[k]
                   - f_2 * pb_x[k] * sg_2[k]
                   + f_1 * pb_x[k] * sg_6[k];
    }
}

}  // namespace simdovl
