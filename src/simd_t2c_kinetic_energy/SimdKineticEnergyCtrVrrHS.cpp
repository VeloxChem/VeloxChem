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


#include "SimdKineticEnergyCtrVrrHS.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_ctr_hs_kinetic_energy_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                                const size_t pa, const size_t fs_s, const size_t fs,
                                const size_t gs, const size_t hs_s, const size_t ncols,
                                const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 1.125 * std::sqrt(14.0) * beta / p;
    const auto f_1 = 0.5625 * std::sqrt(14.0) / p;
    const auto f_2 = 0.9375 * std::sqrt(14.0);
    const auto f_3 = 1.875 * std::sqrt(14.0);
    const auto f_4 = 0.1875 * std::sqrt(14.0);
    const auto f_5 = 1.875 * std::sqrt(14.0) * alpha * beta / p;
    const auto f_6 = 3.75 * std::sqrt(14.0) * alpha * beta / p;
    const auto f_7 = 0.375 * std::sqrt(14.0) * alpha * beta / p;
    const auto f_8 = 1.5 * std::sqrt(35.0);
    const auto f_9 = 3.0 * std::sqrt(35.0) * alpha * beta / p;
    const auto f_10 = 0.125 * std::sqrt(70.0) * beta / p;
    const auto f_11 = std::sqrt(70.0) * beta / p;
    const auto f_12 = 0.0625 * std::sqrt(70.0) / p;
    const auto f_13 = 0.5 * std::sqrt(70.0) / p;
    const auto f_14 = 0.1875 * std::sqrt(70.0);
    const auto f_15 = 1.5 * std::sqrt(70.0);
    const auto f_16 = 0.125 * std::sqrt(70.0);
    const auto f_17 = 0.0625 * std::sqrt(70.0);
    const auto f_18 = 0.5 * std::sqrt(70.0);
    const auto f_19 = 0.375 * std::sqrt(70.0) * alpha * beta / p;
    const auto f_20 = 0.25 * std::sqrt(70.0) * alpha * beta / p;
    const auto f_21 = 3.0 * std::sqrt(70.0) * alpha * beta / p;
    const auto f_22 = 0.125 * std::sqrt(70.0) * alpha * beta / p;
    const auto f_23 = std::sqrt(70.0) * alpha * beta / p;
    const auto f_24 = 0.5 * std::sqrt(105.0);
    const auto f_25 = std::sqrt(105.0);
    const auto f_26 = std::sqrt(105.0) * alpha * beta / p;
    const auto f_27 = 2.0 * std::sqrt(105.0) * alpha * beta / p;
    const auto f_28 = 0.75 * std::sqrt(15.0) * beta / p;
    const auto f_29 = 3.0 * std::sqrt(15.0) * beta / p;
    const auto f_30 = 0.375 * std::sqrt(15.0) / p;
    const auto f_31 = 1.5 * std::sqrt(15.0) / p;
    const auto f_32 = 0.125 * std::sqrt(15.0);
    const auto f_33 = 1.5 * std::sqrt(15.0);
    const auto f_34 = 0.25 * std::sqrt(15.0);
    const auto f_35 = std::sqrt(15.0);
    const auto f_36 = 0.25 * std::sqrt(15.0) * alpha * beta / p;
    const auto f_37 = 0.5 * std::sqrt(15.0) * alpha * beta / p;
    const auto f_38 = 3.0 * std::sqrt(15.0) * alpha * beta / p;
    const auto f_39 = 2.0 * std::sqrt(15.0) * alpha * beta / p;
    const auto f_40 = 6.0 * beta / p;
    const auto f_41 = 3.0 / p;
    const auto f_42 = 3.75 * alpha * beta / p;
    const auto f_43 = 7.5 * alpha * beta / p;
    const auto f_44 = 10.0 * alpha * beta / p;
    const auto f_45 = 2.0 * alpha * beta / p;
    const auto f_46 = 0.5 * std::sqrt(15.0) * beta / p;
    const auto f_47 = 0.25 * std::sqrt(15.0) / p;
    const auto f_48 = 0.25 * std::sqrt(105.0);
    const auto f_49 = 0.5 * std::sqrt(105.0) * alpha * beta / p;
    const auto f_50 = 0.25 * std::sqrt(70.0) * beta / p;
    const auto f_51 = 0.125 * std::sqrt(70.0) / p;
    const auto f_52 = 0.375 * std::sqrt(35.0);
    const auto f_53 = 2.25 * std::sqrt(35.0);
    const auto f_54 = 0.75 * std::sqrt(35.0) * alpha * beta / p;
    const auto f_55 = 4.5 * std::sqrt(35.0) * alpha * beta / p;
    const auto f_56 = 0.75 * std::sqrt(14.0) * beta / p;
    const auto f_57 = 3.75 * std::sqrt(14.0) * beta / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *fs_s_0 = buffer.data(fs_s + 0);
    const auto *fs_s_1 = buffer.data(fs_s + 1);
    const auto *fs_s_2 = buffer.data(fs_s + 2);
    const auto *fs_s_3 = buffer.data(fs_s + 3);
    const auto *fs_s_4 = buffer.data(fs_s + 4);
    const auto *fs_s_5 = buffer.data(fs_s + 5);

    const auto *fs_0 = buffer.data(fs + 0);
    const auto *fs_1 = buffer.data(fs + 1);
    const auto *fs_2 = buffer.data(fs + 2);
    const auto *fs_3 = buffer.data(fs + 3);
    const auto *fs_4 = buffer.data(fs + 4);
    const auto *fs_5 = buffer.data(fs + 5);

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

    const auto *hs_s_0 = buffer.data(hs_s + 0);
    const auto *hs_s_1 = buffer.data(hs_s + 1);
    const auto *hs_s_2 = buffer.data(hs_s + 2);
    const auto *hs_s_3 = buffer.data(hs_s + 3);
    const auto *hs_s_4 = buffer.data(hs_s + 4);
    const auto *hs_s_5 = buffer.data(hs_s + 5);
    const auto *hs_s_6 = buffer.data(hs_s + 6);
    const auto *hs_s_7 = buffer.data(hs_s + 7);
    const auto *hs_s_8 = buffer.data(hs_s + 8);
    const auto *hs_s_9 = buffer.data(hs_s + 9);
    const auto *hs_s_10 = buffer.data(hs_s + 10);
    const auto *hs_s_11 = buffer.data(hs_s + 11);
    const auto *hs_s_12 = buffer.data(hs_s + 12);
    const auto *hs_s_13 = buffer.data(hs_s + 13);
    const auto *hs_s_14 = buffer.data(hs_s + 14);
    const auto *hs_s_15 = buffer.data(hs_s + 15);
    const auto *hs_s_16 = buffer.data(hs_s + 16);
    const auto *hs_s_17 = buffer.data(hs_s + 17);
    const auto *hs_s_18 = buffer.data(hs_s + 18);
    const auto *hs_s_19 = buffer.data(hs_s + 19);
    const auto *hs_s_20 = buffer.data(hs_s + 20);

#pragma omp simd aligned(pa_x, pa_y, fs_s_3, fs_3, gs_0, gs_1, gs_4, gs_6, gs_7, hs_s_1, \
                         hs_s_4, hs_s_6, hs_s_11, hs_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += f_0 * fs_s_3[k]
                  - f_1 * fs_3[k]
                  + f_2 * pa_y[k] * gs_0[k]
                  - f_3 * pa_x[k] * gs_4[k]
                  + f_4 * pa_y[k] * gs_6[k]
                  + f_5 * hs_s_1[k]
                  - f_6 * hs_s_6[k]
                  + f_7 * hs_s_15[k];

        g_1[k] += f_8 * pa_y[k] * gs_1[k]
                  - f_8 * pa_x[k] * gs_7[k]
                  + f_9 * hs_s_4[k]
                  - f_9 * hs_s_11[k];
    }

#pragma omp simd aligned(pa_x, pa_y, fs_s_3, fs_s_4, fs_3, fs_4, gs_0, gs_3, gs_4, gs_6, gs_8, \
                         hs_s_1, hs_s_6, hs_s_8, hs_s_15, hs_s_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_2[k] += -f_10 * fs_s_3[k]
                  + f_11 * fs_s_4[k]
                  + f_12 * fs_3[k]
                  - f_13 * fs_4[k]
                  - f_14 * pa_y[k] * gs_0[k]
                  + f_15 * pa_y[k] * gs_3[k]
                  - f_16 * pa_x[k] * gs_4[k]
                  + f_17 * pa_y[k] * gs_6[k]
                  - f_18 * pa_y[k] * gs_8[k]
                  - f_19 * hs_s_1[k]
                  - f_20 * hs_s_6[k]
                  + f_21 * hs_s_8[k]
                  + f_22 * hs_s_15[k]
                  - f_23 * hs_s_17[k];
    }

#pragma omp simd aligned(pa_x, pa_y, gs_1, gs_7, gs_9, hs_s_4, hs_s_11, \
                         hs_s_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += -f_24 * pa_y[k] * gs_1[k]
                  - f_24 * pa_x[k] * gs_7[k]
                  + f_25 * pa_x[k] * gs_9[k]
                  - f_26 * hs_s_4[k]
                  - f_26 * hs_s_11[k]
                  + f_27 * hs_s_13[k];
    }

#pragma omp simd aligned(pa_x, pa_y, fs_s_3, fs_s_4, fs_3, fs_4, gs_0, gs_3, gs_4, gs_6, gs_8, \
                         gs_10, hs_s_1, hs_s_6, hs_s_8, hs_s_15, hs_s_17, \
                         hs_s_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += -f_28 * fs_s_3[k]
                  + f_29 * fs_s_4[k]
                  + f_30 * fs_3[k]
                  - f_31 * fs_4[k]
                  + f_32 * pa_y[k] * gs_0[k]
                  - f_33 * pa_y[k] * gs_3[k]
                  + f_34 * pa_x[k] * gs_4[k]
                  + f_32 * pa_y[k] * gs_6[k]
                  - f_33 * pa_y[k] * gs_8[k]
                  + f_35 * pa_y[k] * gs_10[k]
                  + f_36 * hs_s_1[k]
                  + f_37 * hs_s_6[k]
                  - f_38 * hs_s_8[k]
                  + f_36 * hs_s_15[k]
                  - f_38 * hs_s_17[k]
                  + f_39 * hs_s_19[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, fs_s_5, fs_5, gs_0, gs_2, gs_5, gs_6, gs_9, gs_10, \
                         hs_s_2, hs_s_7, hs_s_9, hs_s_16, hs_s_18, \
                         hs_s_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += f_40 * fs_s_5[k]
                  - f_41 * fs_5[k]
                  + 1.875 * pa_z[k] * gs_0[k]
                  + 3.75 * pa_z[k] * gs_2[k]
                  - 5.0 * pa_x[k] * gs_5[k]
                  + 1.875 * pa_z[k] * gs_6[k]
                  - 5.0 * pa_y[k] * gs_9[k]
                  + pa_z[k] * gs_10[k]
                  + f_42 * hs_s_2[k]
                  + f_43 * hs_s_7[k]
                  - f_44 * hs_s_9[k]
                  + f_42 * hs_s_16[k]
                  - f_44 * hs_s_18[k]
                  + f_45 * hs_s_20[k];
    }

#pragma omp simd aligned(pa_x, fs_s_0, fs_s_1, fs_s_2, fs_0, fs_1, fs_2, gs_0, gs_2, gs_3, \
                         gs_6, gs_8, gs_10, hs_s_0, hs_s_3, hs_s_5, hs_s_10, hs_s_12, \
                         hs_s_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += -f_46 * fs_s_0[k]
                  - f_46 * fs_s_1[k]
                  + f_29 * fs_s_2[k]
                  + f_47 * fs_0[k]
                  + f_47 * fs_1[k]
                  - f_31 * fs_2[k]
                  + f_32 * pa_x[k] * gs_0[k]
                  + f_34 * pa_x[k] * gs_2[k]
                  - f_33 * pa_x[k] * gs_3[k]
                  + f_32 * pa_x[k] * gs_6[k]
                  - f_33 * pa_x[k] * gs_8[k]
                  + f_35 * pa_x[k] * gs_10[k]
                  + f_36 * hs_s_0[k]
                  + f_37 * hs_s_3[k]
                  - f_38 * hs_s_5[k]
                  + f_36 * hs_s_10[k]
                  - f_38 * hs_s_12[k]
                  + f_39 * hs_s_14[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, gs_0, gs_5, gs_6, gs_9, hs_s_2, hs_s_9, hs_s_16, \
                         hs_s_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_7[k] += -f_48 * pa_z[k] * gs_0[k]
                  + f_24 * pa_x[k] * gs_5[k]
                  + f_48 * pa_z[k] * gs_6[k]
                  - f_24 * pa_y[k] * gs_9[k]
                  - f_49 * hs_s_2[k]
                  + f_26 * hs_s_9[k]
                  + f_49 * hs_s_16[k]
                  - f_26 * hs_s_18[k];
    }

#pragma omp simd aligned(pa_x, fs_s_0, fs_s_1, fs_s_2, fs_0, fs_1, fs_2, gs_0, gs_2, gs_3, \
                         gs_6, gs_8, hs_s_0, hs_s_3, hs_s_5, hs_s_10, \
                         hs_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += f_50 * fs_s_0[k]
                  - f_50 * fs_s_1[k]
                  - f_11 * fs_s_2[k]
                  - f_51 * fs_0[k]
                  + f_51 * fs_1[k]
                  + f_13 * fs_2[k]
                  - f_17 * pa_x[k] * gs_0[k]
                  + f_16 * pa_x[k] * gs_2[k]
                  + f_18 * pa_x[k] * gs_3[k]
                  + f_14 * pa_x[k] * gs_6[k]
                  - f_15 * pa_x[k] * gs_8[k]
                  - f_22 * hs_s_0[k]
                  + f_20 * hs_s_3[k]
                  + f_23 * hs_s_5[k]
                  + f_19 * hs_s_10[k]
                  - f_21 * hs_s_12[k];
    }

#pragma omp simd aligned(pa_z, gs_0, gs_2, gs_6, hs_s_2, hs_s_7, \
                         hs_s_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_9[k] += f_52 * pa_z[k] * gs_0[k]
                  - f_53 * pa_z[k] * gs_2[k]
                  + f_52 * pa_z[k] * gs_6[k]
                  + f_54 * hs_s_2[k]
                  - f_55 * hs_s_7[k]
                  + f_54 * hs_s_16[k];
    }

#pragma omp simd aligned(pa_x, fs_s_0, fs_s_1, fs_0, fs_1, gs_0, gs_2, gs_6, hs_s_0, hs_s_3, \
                         hs_s_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_10[k] += -f_56 * fs_s_0[k]
                   + f_57 * fs_s_1[k]
                   + f_58 * fs_0[k]
                   - f_59 * fs_1[k]
                   + f_4 * pa_x[k] * gs_0[k]
                   - f_3 * pa_x[k] * gs_2[k]
                   + f_2 * pa_x[k] * gs_6[k]
                   + f_7 * hs_s_0[k]
                   - f_6 * hs_s_3[k]
                   + f_5 * hs_s_10[k];
    }
}

}  // namespace simdkin
