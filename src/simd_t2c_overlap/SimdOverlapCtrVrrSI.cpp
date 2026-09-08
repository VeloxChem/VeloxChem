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


#include "SimdOverlapCtrVrrSI.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_ctr_si_overlap_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                         const size_t pb, const size_t sg, const size_t sh, const size_t ncols,
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_13 = buffer.data(sh + 13);
    const auto *sh_14 = buffer.data(sh + 14);

#pragma omp simd aligned(pb_x, pb_y, pb_z, sg_3, sh_0, sh_1, sh_3, sh_4, sh_9, \
                         sh_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += -f_0 * sg_3[k]
                  + f_1 * pb_y[k] * sh_0[k]
                  - f_2 * pb_x[k] * sh_4[k]
                  + f_1 * pb_x[k] * sh_9[k];

        g_1[k] += f_3 * pb_y[k] * sh_1[k]
                  - f_4 * pb_z[k] * sh_4[k]
                  + f_5 * pb_z[k] * sh_9[k];

        g_2[k] += -f_6 * pb_y[k] * sh_0[k]
                  + f_7 * pb_y[k] * sh_3[k]
                  + f_6 * pb_x[k] * sh_9[k]
                  - f_7 * pb_x[k] * sh_11[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, sg_3, sg_7, sh_0, sh_1, sh_3, sh_4, sh_5, sh_9, \
                         sh_11, sh_12, sh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += -f_8 * sg_7[k]
                  - f_9 * pb_y[k] * sh_1[k]
                  - f_10 * pb_z[k] * sh_4[k]
                  + f_11 * pb_y[k] * sh_5[k]
                  + f_12 * pb_z[k] * sh_9[k]
                  - f_13 * pb_y[k] * sh_12[k];

        g_4[k] += f_14 * sg_3[k]
                  + f_15 * pb_y[k] * sh_0[k]
                  - f_16 * pb_y[k] * sh_3[k]
                  + f_17 * pb_x[k] * sh_4[k]
                  + f_15 * pb_x[k] * sh_9[k]
                  - f_16 * pb_x[k] * sh_11[k]
                  + f_16 * pb_x[k] * sh_13[k];
    }

#pragma omp simd aligned(pb_y, pb_z, sg_7, sh_1, sh_4, sh_5, sh_9, sh_12, \
                         sh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += -f_18 * sg_7[k]
                  + f_19 * pb_y[k] * sh_1[k]
                  + f_20 * pb_z[k] * sh_4[k]
                  - f_21 * pb_y[k] * sh_5[k]
                  + f_19 * pb_z[k] * sh_9[k]
                  - f_21 * pb_y[k] * sh_12[k]
                  + f_22 * pb_y[k] * sh_14[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, sg_0, sg_1, sg_2, sg_5, sg_6, sg_8, sh_0, sh_2, \
                         sh_3, sh_6, sh_7, sh_8, sh_9, sh_11, sh_13, \
                         sh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += -f_23 * sg_0[k]
                  - f_24 * sg_1[k]
                  + f_25 * sg_2[k]
                  - f_26 * sg_5[k]
                  + f_27 * sg_6[k]
                  - f_28 * sg_8[k]
                  - 0.3125 * pb_x[k] * sh_0[k]
                  - 0.9375 * pb_x[k] * sh_2[k]
                  + 5.625 * pb_x[k] * sh_3[k]
                  - 0.9375 * pb_x[k] * sh_6[k]
                  + 11.25 * pb_x[k] * sh_7[k]
                  - 7.5 * pb_x[k] * sh_8[k]
                  - 0.3125 * pb_y[k] * sh_9[k]
                  + 5.625 * pb_y[k] * sh_11[k]
                  - 7.5 * pb_y[k] * sh_13[k]
                  + pb_z[k] * sh_14[k];
    }

#pragma omp simd aligned(pb_x, pb_z, sg_4, sh_0, sh_2, sh_5, sh_10, sh_12, \
                         sh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_7[k] += -f_18 * sg_4[k]
                  + f_19 * pb_z[k] * sh_0[k]
                  + f_20 * pb_z[k] * sh_2[k]
                  - f_21 * pb_x[k] * sh_5[k]
                  + f_19 * pb_x[k] * sh_10[k]
                  - f_21 * pb_x[k] * sh_12[k]
                  + f_22 * pb_x[k] * sh_14[k];
    }

#pragma omp simd aligned(pb_x, pb_y, sg_0, sg_1, sg_2, sg_5, sg_6, sh_0, sh_2, sh_3, sh_6, \
                         sh_8, sh_9, sh_11, sh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += f_29 * sg_0[k]
                  + f_30 * sg_1[k]
                  - f_31 * sg_2[k]
                  - f_32 * sg_5[k]
                  + f_31 * sg_6[k]
                  + f_33 * pb_x[k] * sh_0[k]
                  + f_33 * pb_x[k] * sh_2[k]
                  - f_13 * pb_x[k] * sh_3[k]
                  - f_33 * pb_x[k] * sh_6[k]
                  + f_13 * pb_x[k] * sh_8[k]
                  - f_33 * pb_y[k] * sh_9[k]
                  + f_13 * pb_y[k] * sh_11[k]
                  - f_13 * pb_y[k] * sh_13[k];
    }

#pragma omp simd aligned(pb_x, pb_z, sg_4, sh_0, sh_2, sh_5, sh_10, \
                         sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_9[k] += f_8 * sg_4[k]
                  - f_12 * pb_z[k] * sh_0[k]
                  + f_10 * pb_z[k] * sh_2[k]
                  + f_13 * pb_x[k] * sh_5[k]
                  + f_9 * pb_x[k] * sh_10[k]
                  - f_11 * pb_x[k] * sh_12[k];
    }

#pragma omp simd aligned(pb_x, pb_y, sg_0, sg_1, sg_2, sg_6, sh_0, sh_2, sh_3, sh_6, sh_7, \
                         sh_9, sh_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_10[k] += -f_34 * sg_0[k]
                   + f_35 * sg_1[k]
                   + f_36 * sg_2[k]
                   - f_36 * sg_6[k]
                   - f_37 * pb_x[k] * sh_0[k]
                   + f_38 * pb_x[k] * sh_2[k]
                   + f_39 * pb_x[k] * sh_3[k]
                   + f_38 * pb_x[k] * sh_6[k]
                   - f_40 * pb_x[k] * sh_7[k]
                   - f_37 * pb_y[k] * sh_9[k]
                   + f_39 * pb_y[k] * sh_11[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, sg_0, sg_1, sg_5, sh_0, sh_2, sh_6, sh_9, \
                         sh_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_11[k] += f_5 * pb_z[k] * sh_0[k]
                   - f_4 * pb_z[k] * sh_2[k]
                   + f_3 * pb_x[k] * sh_10[k];

        g_12[k] += f_41 * sg_0[k]
                   - f_42 * sg_1[k]
                   + f_43 * sg_5[k]
                   + f_44 * pb_x[k] * sh_0[k]
                   - f_45 * pb_x[k] * sh_2[k]
                   + f_45 * pb_x[k] * sh_6[k]
                   - f_44 * pb_y[k] * sh_9[k];
    }
}

}  // namespace simdovl
