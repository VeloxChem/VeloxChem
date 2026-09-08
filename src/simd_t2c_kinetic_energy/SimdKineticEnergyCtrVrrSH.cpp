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


#include "SimdKineticEnergyCtrVrrSH.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_ctr_sh_kinetic_energy_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                                const size_t pb, const size_t sf_s, const size_t sh_s,
                                const size_t sf, const size_t sg, const size_t ncols,
                                const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 1.125 * std::sqrt(14.0) * alpha / p;
    const auto f_1 = 1.875 * std::sqrt(14.0) * alpha * beta / p;
    const auto f_2 = 3.75 * std::sqrt(14.0) * alpha * beta / p;
    const auto f_3 = 0.375 * std::sqrt(14.0) * alpha * beta / p;
    const auto f_4 = 0.5625 * std::sqrt(14.0) / p;
    const auto f_5 = 0.9375 * std::sqrt(14.0);
    const auto f_6 = 1.875 * std::sqrt(14.0);
    const auto f_7 = 0.1875 * std::sqrt(14.0);
    const auto f_8 = 3.0 * std::sqrt(35.0) * alpha * beta / p;
    const auto f_9 = 1.5 * std::sqrt(35.0);
    const auto f_10 = 0.125 * std::sqrt(70.0) * alpha / p;
    const auto f_11 = std::sqrt(70.0) * alpha / p;
    const auto f_12 = 0.375 * std::sqrt(70.0) * alpha * beta / p;
    const auto f_13 = 0.25 * std::sqrt(70.0) * alpha * beta / p;
    const auto f_14 = 3.0 * std::sqrt(70.0) * alpha * beta / p;
    const auto f_15 = 0.125 * std::sqrt(70.0) * alpha * beta / p;
    const auto f_16 = std::sqrt(70.0) * alpha * beta / p;
    const auto f_17 = 0.0625 * std::sqrt(70.0) / p;
    const auto f_18 = 0.5 * std::sqrt(70.0) / p;
    const auto f_19 = 0.1875 * std::sqrt(70.0);
    const auto f_20 = 1.5 * std::sqrt(70.0);
    const auto f_21 = 0.125 * std::sqrt(70.0);
    const auto f_22 = 0.0625 * std::sqrt(70.0);
    const auto f_23 = 0.5 * std::sqrt(70.0);
    const auto f_24 = std::sqrt(105.0) * alpha * beta / p;
    const auto f_25 = 2.0 * std::sqrt(105.0) * alpha * beta / p;
    const auto f_26 = 0.5 * std::sqrt(105.0);
    const auto f_27 = std::sqrt(105.0);
    const auto f_28 = 0.75 * std::sqrt(15.0) * alpha / p;
    const auto f_29 = 3.0 * std::sqrt(15.0) * alpha / p;
    const auto f_30 = 0.25 * std::sqrt(15.0) * alpha * beta / p;
    const auto f_31 = 0.5 * std::sqrt(15.0) * alpha * beta / p;
    const auto f_32 = 3.0 * std::sqrt(15.0) * alpha * beta / p;
    const auto f_33 = 2.0 * std::sqrt(15.0) * alpha * beta / p;
    const auto f_34 = 0.375 * std::sqrt(15.0) / p;
    const auto f_35 = 1.5 * std::sqrt(15.0) / p;
    const auto f_36 = 0.125 * std::sqrt(15.0);
    const auto f_37 = 1.5 * std::sqrt(15.0);
    const auto f_38 = 0.25 * std::sqrt(15.0);
    const auto f_39 = std::sqrt(15.0);
    const auto f_40 = 6.0 * alpha / p;
    const auto f_41 = 3.75 * alpha * beta / p;
    const auto f_42 = 7.5 * alpha * beta / p;
    const auto f_43 = 10.0 * alpha * beta / p;
    const auto f_44 = 2.0 * alpha * beta / p;
    const auto f_45 = 3.0 / p;
    const auto f_46 = 0.5 * std::sqrt(15.0) * alpha / p;
    const auto f_47 = 0.25 * std::sqrt(15.0) / p;
    const auto f_48 = 0.5 * std::sqrt(105.0) * alpha * beta / p;
    const auto f_49 = 0.25 * std::sqrt(105.0);
    const auto f_50 = 0.25 * std::sqrt(70.0) * alpha / p;
    const auto f_51 = 0.125 * std::sqrt(70.0) / p;
    const auto f_52 = 0.75 * std::sqrt(35.0) * alpha * beta / p;
    const auto f_53 = 4.5 * std::sqrt(35.0) * alpha * beta / p;
    const auto f_54 = 0.375 * std::sqrt(35.0);
    const auto f_55 = 2.25 * std::sqrt(35.0);
    const auto f_56 = 0.75 * std::sqrt(14.0) * alpha / p;
    const auto f_57 = 3.75 * std::sqrt(14.0) * alpha / p;
    const auto f_58 = 0.375 * std::sqrt(14.0) / p;
    const auto f_59 = 1.875 * std::sqrt(14.0) / p;

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

    const auto *sf_s_0 = buffer.data(sf_s + 0);
    const auto *sf_s_1 = buffer.data(sf_s + 1);
    const auto *sf_s_2 = buffer.data(sf_s + 2);
    const auto *sf_s_3 = buffer.data(sf_s + 3);
    const auto *sf_s_4 = buffer.data(sf_s + 4);
    const auto *sf_s_5 = buffer.data(sf_s + 5);

    const auto *sh_s_0 = buffer.data(sh_s + 0);
    const auto *sh_s_1 = buffer.data(sh_s + 1);
    const auto *sh_s_2 = buffer.data(sh_s + 2);
    const auto *sh_s_3 = buffer.data(sh_s + 3);
    const auto *sh_s_4 = buffer.data(sh_s + 4);
    const auto *sh_s_5 = buffer.data(sh_s + 5);
    const auto *sh_s_6 = buffer.data(sh_s + 6);
    const auto *sh_s_7 = buffer.data(sh_s + 7);
    const auto *sh_s_8 = buffer.data(sh_s + 8);
    const auto *sh_s_9 = buffer.data(sh_s + 9);
    const auto *sh_s_10 = buffer.data(sh_s + 10);
    const auto *sh_s_11 = buffer.data(sh_s + 11);
    const auto *sh_s_12 = buffer.data(sh_s + 12);
    const auto *sh_s_13 = buffer.data(sh_s + 13);
    const auto *sh_s_14 = buffer.data(sh_s + 14);
    const auto *sh_s_15 = buffer.data(sh_s + 15);
    const auto *sh_s_16 = buffer.data(sh_s + 16);
    const auto *sh_s_17 = buffer.data(sh_s + 17);
    const auto *sh_s_18 = buffer.data(sh_s + 18);
    const auto *sh_s_19 = buffer.data(sh_s + 19);
    const auto *sh_s_20 = buffer.data(sh_s + 20);

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

#pragma omp simd aligned(pb_x, pb_y, sf_s_3, sh_s_1, sh_s_4, sh_s_6, sh_s_11, sh_s_15, sf_3, \
                         sg_0, sg_1, sg_4, sg_6, sg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += f_0 * sf_s_3[k]
                  + f_1 * sh_s_1[k]
                  - f_2 * sh_s_6[k]
                  + f_3 * sh_s_15[k]
                  - f_4 * sf_3[k]
                  + f_5 * pb_y[k] * sg_0[k]
                  - f_6 * pb_x[k] * sg_4[k]
                  + f_7 * pb_y[k] * sg_6[k];

        g_1[k] += f_8 * sh_s_4[k]
                  - f_8 * sh_s_11[k]
                  + f_9 * pb_y[k] * sg_1[k]
                  - f_9 * pb_x[k] * sg_7[k];
    }

#pragma omp simd aligned(pb_x, pb_y, sf_s_3, sf_s_4, sh_s_1, sh_s_6, sh_s_8, sh_s_15, sh_s_17, \
                         sf_3, sf_4, sg_0, sg_3, sg_4, sg_6, sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_2[k] += -f_10 * sf_s_3[k]
                  + f_11 * sf_s_4[k]
                  - f_12 * sh_s_1[k]
                  - f_13 * sh_s_6[k]
                  + f_14 * sh_s_8[k]
                  + f_15 * sh_s_15[k]
                  - f_16 * sh_s_17[k]
                  + f_17 * sf_3[k]
                  - f_18 * sf_4[k]
                  - f_19 * pb_y[k] * sg_0[k]
                  + f_20 * pb_y[k] * sg_3[k]
                  - f_21 * pb_x[k] * sg_4[k]
                  + f_22 * pb_y[k] * sg_6[k]
                  - f_23 * pb_y[k] * sg_8[k];
    }

#pragma omp simd aligned(pb_x, pb_y, sh_s_4, sh_s_11, sh_s_13, sg_1, sg_7, \
                         sg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += -f_24 * sh_s_4[k]
                  - f_24 * sh_s_11[k]
                  + f_25 * sh_s_13[k]
                  - f_26 * pb_y[k] * sg_1[k]
                  - f_26 * pb_x[k] * sg_7[k]
                  + f_27 * pb_x[k] * sg_9[k];
    }

#pragma omp simd aligned(pb_x, pb_y, sf_s_3, sf_s_4, sh_s_1, sh_s_6, sh_s_8, sh_s_15, sh_s_17, \
                         sh_s_19, sf_3, sf_4, sg_0, sg_3, sg_4, sg_6, sg_8, \
                         sg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += -f_28 * sf_s_3[k]
                  + f_29 * sf_s_4[k]
                  + f_30 * sh_s_1[k]
                  + f_31 * sh_s_6[k]
                  - f_32 * sh_s_8[k]
                  + f_30 * sh_s_15[k]
                  - f_32 * sh_s_17[k]
                  + f_33 * sh_s_19[k]
                  + f_34 * sf_3[k]
                  - f_35 * sf_4[k]
                  + f_36 * pb_y[k] * sg_0[k]
                  - f_37 * pb_y[k] * sg_3[k]
                  + f_38 * pb_x[k] * sg_4[k]
                  + f_36 * pb_y[k] * sg_6[k]
                  - f_37 * pb_y[k] * sg_8[k]
                  + f_39 * pb_y[k] * sg_10[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, sf_s_5, sh_s_2, sh_s_7, sh_s_9, sh_s_16, sh_s_18, \
                         sh_s_20, sf_5, sg_0, sg_2, sg_5, sg_6, sg_9, \
                         sg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += f_40 * sf_s_5[k]
                  + f_41 * sh_s_2[k]
                  + f_42 * sh_s_7[k]
                  - f_43 * sh_s_9[k]
                  + f_41 * sh_s_16[k]
                  - f_43 * sh_s_18[k]
                  + f_44 * sh_s_20[k]
                  - f_45 * sf_5[k]
                  + 1.875 * pb_z[k] * sg_0[k]
                  + 3.75 * pb_z[k] * sg_2[k]
                  - 5.0 * pb_x[k] * sg_5[k]
                  + 1.875 * pb_z[k] * sg_6[k]
                  - 5.0 * pb_y[k] * sg_9[k]
                  + pb_z[k] * sg_10[k];
    }

#pragma omp simd aligned(pb_x, sf_s_0, sf_s_1, sf_s_2, sh_s_0, sh_s_3, sh_s_5, sh_s_10, \
                         sh_s_12, sh_s_14, sf_0, sf_1, sf_2, sg_0, sg_2, sg_3, sg_6, sg_8, \
                         sg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += -f_46 * sf_s_0[k]
                  - f_46 * sf_s_1[k]
                  + f_29 * sf_s_2[k]
                  + f_30 * sh_s_0[k]
                  + f_31 * sh_s_3[k]
                  - f_32 * sh_s_5[k]
                  + f_30 * sh_s_10[k]
                  - f_32 * sh_s_12[k]
                  + f_33 * sh_s_14[k]
                  + f_47 * sf_0[k]
                  + f_47 * sf_1[k]
                  - f_35 * sf_2[k]
                  + f_36 * pb_x[k] * sg_0[k]
                  + f_38 * pb_x[k] * sg_2[k]
                  - f_37 * pb_x[k] * sg_3[k]
                  + f_36 * pb_x[k] * sg_6[k]
                  - f_37 * pb_x[k] * sg_8[k]
                  + f_39 * pb_x[k] * sg_10[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, sh_s_2, sh_s_9, sh_s_16, sh_s_18, sg_0, sg_5, sg_6, \
                         sg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_7[k] += -f_48 * sh_s_2[k]
                  + f_24 * sh_s_9[k]
                  + f_48 * sh_s_16[k]
                  - f_24 * sh_s_18[k]
                  - f_49 * pb_z[k] * sg_0[k]
                  + f_26 * pb_x[k] * sg_5[k]
                  + f_49 * pb_z[k] * sg_6[k]
                  - f_26 * pb_y[k] * sg_9[k];
    }

#pragma omp simd aligned(pb_x, sf_s_0, sf_s_1, sf_s_2, sh_s_0, sh_s_3, sh_s_5, sh_s_10, \
                         sh_s_12, sf_0, sf_1, sf_2, sg_0, sg_2, sg_3, sg_6, \
                         sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += f_50 * sf_s_0[k]
                  - f_50 * sf_s_1[k]
                  - f_11 * sf_s_2[k]
                  - f_15 * sh_s_0[k]
                  + f_13 * sh_s_3[k]
                  + f_16 * sh_s_5[k]
                  + f_12 * sh_s_10[k]
                  - f_14 * sh_s_12[k]
                  - f_51 * sf_0[k]
                  + f_51 * sf_1[k]
                  + f_18 * sf_2[k]
                  - f_22 * pb_x[k] * sg_0[k]
                  + f_21 * pb_x[k] * sg_2[k]
                  + f_23 * pb_x[k] * sg_3[k]
                  + f_19 * pb_x[k] * sg_6[k]
                  - f_20 * pb_x[k] * sg_8[k];
    }

#pragma omp simd aligned(pb_z, sh_s_2, sh_s_7, sh_s_16, sg_0, sg_2, \
                         sg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_9[k] += f_52 * sh_s_2[k]
                  - f_53 * sh_s_7[k]
                  + f_52 * sh_s_16[k]
                  + f_54 * pb_z[k] * sg_0[k]
                  - f_55 * pb_z[k] * sg_2[k]
                  + f_54 * pb_z[k] * sg_6[k];
    }

#pragma omp simd aligned(pb_x, sf_s_0, sf_s_1, sh_s_0, sh_s_3, sh_s_10, sf_0, sf_1, sg_0, \
                         sg_2, sg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_10[k] += -f_56 * sf_s_0[k]
                   + f_57 * sf_s_1[k]
                   + f_3 * sh_s_0[k]
                   - f_2 * sh_s_3[k]
                   + f_1 * sh_s_10[k]
                   + f_58 * sf_0[k]
                   - f_59 * sf_1[k]
                   + f_7 * pb_x[k] * sg_0[k]
                   - f_6 * pb_x[k] * sg_2[k]
                   + f_5 * pb_x[k] * sg_6[k];
    }
}

}  // namespace simdkin
