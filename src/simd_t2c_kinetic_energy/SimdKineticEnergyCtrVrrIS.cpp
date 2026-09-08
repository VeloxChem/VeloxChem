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


#include "SimdKineticEnergyCtrVrrIS.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_ctr_is_kinetic_energy_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                                const size_t pa, const size_t gs_s, const size_t gs,
                                const size_t hs, const size_t is_s, const size_t ncols,
                                const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 1.25 * std::sqrt(462.0) * beta / p;
    const auto f_1 = 0.625 * std::sqrt(462.0) / p;
    const auto f_2 = 0.1875 * std::sqrt(462.0);
    const auto f_3 = 0.625 * std::sqrt(462.0);
    const auto f_4 = 0.375 * std::sqrt(462.0) * alpha * beta / p;
    const auto f_5 = 1.25 * std::sqrt(462.0) * alpha * beta / p;
    const auto f_6 = 0.9375 * std::sqrt(154.0);
    const auto f_7 = 1.875 * std::sqrt(154.0);
    const auto f_8 = 0.1875 * std::sqrt(154.0);
    const auto f_9 = 1.875 * std::sqrt(154.0) * alpha * beta / p;
    const auto f_10 = 3.75 * std::sqrt(154.0) * alpha * beta / p;
    const auto f_11 = 0.375 * std::sqrt(154.0) * alpha * beta / p;
    const auto f_12 = 0.75 * std::sqrt(7.0);
    const auto f_13 = 7.5 * std::sqrt(7.0);
    const auto f_14 = 1.5 * std::sqrt(7.0) * alpha * beta / p;
    const auto f_15 = 15.0 * std::sqrt(7.0) * alpha * beta / p;
    const auto f_16 = std::sqrt(210.0) * beta / p;
    const auto f_17 = 0.5 * std::sqrt(210.0) / p;
    const auto f_18 = 0.5625 * std::sqrt(210.0);
    const auto f_19 = 0.375 * std::sqrt(210.0);
    const auto f_20 = 1.5 * std::sqrt(210.0);
    const auto f_21 = 0.1875 * std::sqrt(210.0);
    const auto f_22 = 0.5 * std::sqrt(210.0);
    const auto f_23 = 1.125 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_24 = 0.75 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_25 = 3.0 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_26 = 0.375 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_27 = std::sqrt(210.0) * alpha * beta / p;
    const auto f_28 = 0.25 * std::sqrt(210.0) * beta / p;
    const auto f_29 = 0.125 * std::sqrt(210.0) / p;
    const auto f_30 = 0.0625 * std::sqrt(210.0);
    const auto f_31 = std::sqrt(210.0);
    const auto f_32 = 0.125 * std::sqrt(210.0);
    const auto f_33 = 0.125 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_34 = 0.25 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_35 = 2.0 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_36 = 5.0 * std::sqrt(21.0) * beta / p;
    const auto f_37 = 2.5 * std::sqrt(21.0) / p;
    const auto f_38 = 0.625 * std::sqrt(21.0);
    const auto f_39 = 1.25 * std::sqrt(21.0);
    const auto f_40 = 2.5 * std::sqrt(21.0);
    const auto f_41 = std::sqrt(21.0);
    const auto f_42 = 1.25 * std::sqrt(21.0) * alpha * beta / p;
    const auto f_43 = 2.5 * std::sqrt(21.0) * alpha * beta / p;
    const auto f_44 = 5.0 * std::sqrt(21.0) * alpha * beta / p;
    const auto f_45 = 2.0 * std::sqrt(21.0) * alpha * beta / p;
    const auto f_46 = 1.5625 * beta / p;
    const auto f_47 = 2.8125 * beta / p;
    const auto f_48 = 16.875 * beta / p;
    const auto f_49 = 2.5 * beta / p;
    const auto f_50 = 28.125 * beta / p;
    const auto f_51 = 10.0 * beta / p;
    const auto f_52 = 0.78125 / p;
    const auto f_53 = 1.40625 / p;
    const auto f_54 = 8.4375 / p;
    const auto f_55 = 1.25 / p;
    const auto f_56 = 14.0625 / p;
    const auto f_57 = 5.0 / p;
    const auto f_58 = 0.625 * alpha * beta / p;
    const auto f_59 = 1.875 * alpha * beta / p;
    const auto f_60 = 11.25 * alpha * beta / p;
    const auto f_61 = 22.5 * alpha * beta / p;
    const auto f_62 = 15.0 * alpha * beta / p;
    const auto f_63 = 2.0 * alpha * beta / p;
    const auto f_64 = 0.15625 * std::sqrt(210.0) * beta / p;
    const auto f_65 = 0.09375 * std::sqrt(210.0) * beta / p;
    const auto f_66 = 1.5 * std::sqrt(210.0) * beta / p;
    const auto f_67 = 0.1875 * std::sqrt(210.0) * beta / p;
    const auto f_68 = 0.078125 * std::sqrt(210.0) / p;
    const auto f_69 = 0.046875 * std::sqrt(210.0) / p;
    const auto f_70 = 0.75 * std::sqrt(210.0) / p;
    const auto f_71 = 0.09375 * std::sqrt(210.0) / p;
    const auto f_72 = 0.03125 * std::sqrt(210.0);
    const auto f_73 = 0.0625 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_74 = 0.9375 * std::sqrt(7.0) * beta / p;
    const auto f_75 = 2.8125 * std::sqrt(7.0) * beta / p;
    const auto f_76 = 5.625 * std::sqrt(7.0) * beta / p;
    const auto f_77 = 0.46875 * std::sqrt(7.0) / p;
    const auto f_78 = 1.40625 * std::sqrt(7.0) / p;
    const auto f_79 = 2.8125 * std::sqrt(7.0) / p;
    const auto f_80 = 0.1875 * std::sqrt(7.0);
    const auto f_81 = 0.9375 * std::sqrt(7.0);
    const auto f_82 = 1.875 * std::sqrt(7.0);
    const auto f_83 = 11.25 * std::sqrt(7.0);
    const auto f_84 = 0.375 * std::sqrt(7.0) * alpha * beta / p;
    const auto f_85 = 1.875 * std::sqrt(7.0) * alpha * beta / p;
    const auto f_86 = 3.75 * std::sqrt(7.0) * alpha * beta / p;
    const auto f_87 = 22.5 * std::sqrt(7.0) * alpha * beta / p;
    const auto f_88 = 0.15625 * std::sqrt(462.0) * beta / p;
    const auto f_89 = 1.40625 * std::sqrt(462.0) * beta / p;
    const auto f_90 = 0.3125 * std::sqrt(462.0) * beta / p;
    const auto f_91 = 0.078125 * std::sqrt(462.0) / p;
    const auto f_92 = 0.703125 * std::sqrt(462.0) / p;
    const auto f_93 = 0.15625 * std::sqrt(462.0) / p;
    const auto f_94 = 0.03125 * std::sqrt(462.0);
    const auto f_95 = 0.46875 * std::sqrt(462.0);
    const auto f_96 = 0.0625 * std::sqrt(462.0) * alpha * beta / p;
    const auto f_97 = 0.9375 * std::sqrt(462.0) * alpha * beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *gs_s_0 = buffer.data(gs_s + 0);
    const auto *gs_s_1 = buffer.data(gs_s + 1);
    const auto *gs_s_2 = buffer.data(gs_s + 2);
    const auto *gs_s_3 = buffer.data(gs_s + 3);
    const auto *gs_s_4 = buffer.data(gs_s + 4);
    const auto *gs_s_5 = buffer.data(gs_s + 5);
    const auto *gs_s_6 = buffer.data(gs_s + 6);
    const auto *gs_s_7 = buffer.data(gs_s + 7);
    const auto *gs_s_8 = buffer.data(gs_s + 8);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);

    const auto *is_s_0 = buffer.data(is_s + 0);
    const auto *is_s_1 = buffer.data(is_s + 1);
    const auto *is_s_2 = buffer.data(is_s + 2);
    const auto *is_s_3 = buffer.data(is_s + 3);
    const auto *is_s_4 = buffer.data(is_s + 4);
    const auto *is_s_5 = buffer.data(is_s + 5);
    const auto *is_s_6 = buffer.data(is_s + 6);
    const auto *is_s_7 = buffer.data(is_s + 7);
    const auto *is_s_8 = buffer.data(is_s + 8);
    const auto *is_s_9 = buffer.data(is_s + 9);
    const auto *is_s_10 = buffer.data(is_s + 10);
    const auto *is_s_11 = buffer.data(is_s + 11);
    const auto *is_s_12 = buffer.data(is_s + 12);
    const auto *is_s_13 = buffer.data(is_s + 13);
    const auto *is_s_14 = buffer.data(is_s + 14);
    const auto *is_s_15 = buffer.data(is_s + 15);
    const auto *is_s_16 = buffer.data(is_s + 16);
    const auto *is_s_17 = buffer.data(is_s + 17);
    const auto *is_s_18 = buffer.data(is_s + 18);
    const auto *is_s_19 = buffer.data(is_s + 19);
    const auto *is_s_20 = buffer.data(is_s + 20);
    const auto *is_s_21 = buffer.data(is_s + 21);
    const auto *is_s_22 = buffer.data(is_s + 22);
    const auto *is_s_23 = buffer.data(is_s + 23);
    const auto *is_s_24 = buffer.data(is_s + 24);
    const auto *is_s_25 = buffer.data(is_s + 25);
    const auto *is_s_26 = buffer.data(is_s + 26);
    const auto *is_s_27 = buffer.data(is_s + 27);

#pragma omp simd aligned(pa_x, pa_y, gs_s_3, gs_3, hs_0, hs_4, hs_9, is_s_1, is_s_6, \
                         is_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += f_0 * gs_s_3[k]
                  - f_1 * gs_3[k]
                  + f_2 * pa_y[k] * hs_0[k]
                  - f_3 * pa_x[k] * hs_4[k]
                  + f_2 * pa_x[k] * hs_9[k]
                  + f_4 * is_s_1[k]
                  - f_5 * is_s_6[k]
                  + f_4 * is_s_15[k];
    }

#pragma omp simd aligned(pa_y, pa_z, hs_1, hs_4, hs_9, is_s_4, is_s_11, \
                         is_s_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_1[k] += f_6 * pa_y[k] * hs_1[k]
                  - f_7 * pa_z[k] * hs_4[k]
                  + f_8 * pa_z[k] * hs_9[k]
                  + f_9 * is_s_4[k]
                  - f_10 * is_s_11[k]
                  + f_11 * is_s_22[k];
    }

#pragma omp simd aligned(pa_x, pa_y, hs_0, hs_3, hs_9, hs_11, is_s_1, is_s_8, is_s_15, \
                         is_s_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_2[k] += -f_12 * pa_y[k] * hs_0[k]
                  + f_13 * pa_y[k] * hs_3[k]
                  + f_12 * pa_x[k] * hs_9[k]
                  - f_13 * pa_x[k] * hs_11[k]
                  - f_14 * is_s_1[k]
                  + f_15 * is_s_8[k]
                  + f_14 * is_s_15[k]
                  - f_15 * is_s_17[k];
    }

#pragma omp simd aligned(pa_y, pa_z, gs_s_7, gs_7, hs_1, hs_4, hs_5, hs_9, hs_12, is_s_4, \
                         is_s_11, is_s_13, is_s_22, is_s_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += f_16 * gs_s_7[k]
                  - f_17 * gs_7[k]
                  - f_18 * pa_y[k] * hs_1[k]
                  - f_19 * pa_z[k] * hs_4[k]
                  + f_20 * pa_y[k] * hs_5[k]
                  + f_21 * pa_z[k] * hs_9[k]
                  - f_22 * pa_y[k] * hs_12[k]
                  - f_23 * is_s_4[k]
                  - f_24 * is_s_11[k]
                  + f_25 * is_s_13[k]
                  + f_26 * is_s_22[k]
                  - f_27 * is_s_24[k];
    }

#pragma omp simd aligned(pa_x, pa_y, gs_s_3, gs_3, hs_0, hs_3, hs_4, hs_9, hs_11, hs_13, \
                         is_s_1, is_s_6, is_s_8, is_s_15, is_s_17, \
                         is_s_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += -f_28 * gs_s_3[k]
                  + f_29 * gs_3[k]
                  + f_30 * pa_y[k] * hs_0[k]
                  - f_31 * pa_y[k] * hs_3[k]
                  + f_32 * pa_x[k] * hs_4[k]
                  + f_30 * pa_x[k] * hs_9[k]
                  - f_31 * pa_x[k] * hs_11[k]
                  + f_31 * pa_x[k] * hs_13[k]
                  + f_33 * is_s_1[k]
                  + f_34 * is_s_6[k]
                  - f_35 * is_s_8[k]
                  + f_33 * is_s_15[k]
                  - f_35 * is_s_17[k]
                  + f_35 * is_s_19[k];
    }

#pragma omp simd aligned(pa_y, pa_z, gs_s_7, gs_7, hs_1, hs_4, hs_5, hs_9, hs_12, hs_14, \
                         is_s_4, is_s_11, is_s_13, is_s_22, is_s_24, \
                         is_s_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += f_36 * gs_s_7[k]
                  - f_37 * gs_7[k]
                  + f_38 * pa_y[k] * hs_1[k]
                  + f_39 * pa_z[k] * hs_4[k]
                  - f_40 * pa_y[k] * hs_5[k]
                  + f_38 * pa_z[k] * hs_9[k]
                  - f_40 * pa_y[k] * hs_12[k]
                  + f_41 * pa_y[k] * hs_14[k]
                  + f_42 * is_s_4[k]
                  + f_43 * is_s_11[k]
                  - f_44 * is_s_13[k]
                  + f_42 * is_s_22[k]
                  - f_44 * is_s_24[k]
                  + f_45 * is_s_26[k];
    }

#pragma omp simd aligned(pa_x, pa_y, pa_z, gs_s_0, gs_s_1, gs_s_2, gs_s_5, gs_s_6, gs_s_8, \
                         gs_0, gs_1, gs_2, gs_5, gs_6, gs_8, hs_0, hs_2, hs_3, hs_6, hs_7, \
                         hs_8, hs_9, hs_11, hs_13, hs_14, is_s_0, is_s_3, is_s_5, is_s_10, \
                         is_s_12, is_s_14, is_s_21, is_s_23, is_s_25, \
                         is_s_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += f_46 * gs_s_0[k]
                  + f_47 * gs_s_1[k]
                  - f_48 * gs_s_2[k]
                  + f_49 * gs_s_5[k]
                  - f_50 * gs_s_6[k]
                  + f_51 * gs_s_8[k]
                  - f_52 * gs_0[k]
                  - f_53 * gs_1[k]
                  + f_54 * gs_2[k]
                  - f_55 * gs_5[k]
                  + f_56 * gs_6[k]
                  - f_57 * gs_8[k]
                  - 0.3125 * pa_x[k] * hs_0[k]
                  - 0.9375 * pa_x[k] * hs_2[k]
                  + 5.625 * pa_x[k] * hs_3[k]
                  - 0.9375 * pa_x[k] * hs_6[k]
                  + 11.25 * pa_x[k] * hs_7[k]
                  - 7.5 * pa_x[k] * hs_8[k]
                  - 0.3125 * pa_y[k] * hs_9[k]
                  + 5.625 * pa_y[k] * hs_11[k]
                  - 7.5 * pa_y[k] * hs_13[k]
                  + pa_z[k] * hs_14[k]
                  - f_58 * is_s_0[k]
                  - f_59 * is_s_3[k]
                  + f_60 * is_s_5[k]
                  - f_59 * is_s_10[k]
                  + f_61 * is_s_12[k]
                  - f_62 * is_s_14[k]
                  - f_58 * is_s_21[k]
                  + f_60 * is_s_23[k]
                  - f_62 * is_s_25[k]
                  + f_63 * is_s_27[k];
    }

#pragma omp simd aligned(pa_x, pa_z, gs_s_4, gs_4, hs_0, hs_2, hs_5, hs_10, hs_12, hs_14, \
                         is_s_2, is_s_7, is_s_9, is_s_16, is_s_18, \
                         is_s_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_7[k] += f_36 * gs_s_4[k]
                  - f_37 * gs_4[k]
                  + f_38 * pa_z[k] * hs_0[k]
                  + f_39 * pa_z[k] * hs_2[k]
                  - f_40 * pa_x[k] * hs_5[k]
                  + f_38 * pa_x[k] * hs_10[k]
                  - f_40 * pa_x[k] * hs_12[k]
                  + f_41 * pa_x[k] * hs_14[k]
                  + f_42 * is_s_2[k]
                  + f_43 * is_s_7[k]
                  - f_44 * is_s_9[k]
                  + f_42 * is_s_16[k]
                  - f_44 * is_s_18[k]
                  + f_45 * is_s_20[k];
    }

#pragma omp simd aligned(pa_x, pa_y, gs_s_0, gs_s_1, gs_s_2, gs_s_5, gs_s_6, gs_0, gs_1, gs_2, \
                         gs_5, gs_6, hs_0, hs_2, hs_3, hs_6, hs_8, hs_9, hs_11, hs_13, is_s_0, \
                         is_s_3, is_s_5, is_s_10, is_s_14, is_s_21, is_s_23, \
                         is_s_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += -f_64 * gs_s_0[k]
                  - f_65 * gs_s_1[k]
                  + f_66 * gs_s_2[k]
                  + f_67 * gs_s_5[k]
                  - f_66 * gs_s_6[k]
                  + f_68 * gs_0[k]
                  + f_69 * gs_1[k]
                  - f_70 * gs_2[k]
                  - f_71 * gs_5[k]
                  + f_70 * gs_6[k]
                  + f_72 * pa_x[k] * hs_0[k]
                  + f_72 * pa_x[k] * hs_2[k]
                  - f_22 * pa_x[k] * hs_3[k]
                  - f_72 * pa_x[k] * hs_6[k]
                  + f_22 * pa_x[k] * hs_8[k]
                  - f_72 * pa_y[k] * hs_9[k]
                  + f_22 * pa_y[k] * hs_11[k]
                  - f_22 * pa_y[k] * hs_13[k]
                  + f_73 * is_s_0[k]
                  + f_73 * is_s_3[k]
                  - f_27 * is_s_5[k]
                  - f_73 * is_s_10[k]
                  + f_27 * is_s_14[k]
                  - f_73 * is_s_21[k]
                  + f_27 * is_s_23[k]
                  - f_27 * is_s_25[k];
    }

#pragma omp simd aligned(pa_x, pa_z, gs_s_4, gs_4, hs_0, hs_2, hs_5, hs_10, hs_12, is_s_2, \
                         is_s_7, is_s_9, is_s_16, is_s_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_9[k] += -f_16 * gs_s_4[k]
                  + f_17 * gs_4[k]
                  - f_21 * pa_z[k] * hs_0[k]
                  + f_19 * pa_z[k] * hs_2[k]
                  + f_22 * pa_x[k] * hs_5[k]
                  + f_18 * pa_x[k] * hs_10[k]
                  - f_20 * pa_x[k] * hs_12[k]
                  - f_26 * is_s_2[k]
                  + f_24 * is_s_7[k]
                  + f_27 * is_s_9[k]
                  + f_23 * is_s_16[k]
                  - f_25 * is_s_18[k];
    }

#pragma omp simd aligned(pa_x, pa_y, gs_s_0, gs_s_1, gs_s_2, gs_s_6, gs_0, gs_1, gs_2, gs_6, \
                         hs_0, hs_2, hs_3, hs_6, hs_7, hs_9, hs_11, is_s_0, is_s_3, is_s_5, \
                         is_s_10, is_s_12, is_s_21, is_s_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_10[k] += f_74 * gs_s_0[k]
                   - f_75 * gs_s_1[k]
                   - f_76 * gs_s_2[k]
                   + f_76 * gs_s_6[k]
                   - f_77 * gs_0[k]
                   + f_78 * gs_1[k]
                   + f_79 * gs_2[k]
                   - f_79 * gs_6[k]
                   - f_80 * pa_x[k] * hs_0[k]
                   + f_81 * pa_x[k] * hs_2[k]
                   + f_82 * pa_x[k] * hs_3[k]
                   + f_81 * pa_x[k] * hs_6[k]
                   - f_83 * pa_x[k] * hs_7[k]
                   - f_80 * pa_y[k] * hs_9[k]
                   + f_82 * pa_y[k] * hs_11[k]
                   - f_84 * is_s_0[k]
                   + f_85 * is_s_3[k]
                   + f_86 * is_s_5[k]
                   + f_85 * is_s_10[k]
                   - f_87 * is_s_12[k]
                   - f_84 * is_s_21[k]
                   + f_86 * is_s_23[k];
    }

#pragma omp simd aligned(pa_x, pa_z, hs_0, hs_2, hs_10, is_s_2, is_s_7, \
                         is_s_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_11[k] += f_8 * pa_z[k] * hs_0[k]
                   - f_7 * pa_z[k] * hs_2[k]
                   + f_6 * pa_x[k] * hs_10[k]
                   + f_11 * is_s_2[k]
                   - f_10 * is_s_7[k]
                   + f_9 * is_s_16[k];
    }

#pragma omp simd aligned(pa_x, pa_y, gs_s_0, gs_s_1, gs_s_5, gs_0, gs_1, gs_5, hs_0, hs_2, \
                         hs_6, hs_9, is_s_0, is_s_3, is_s_10, is_s_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_12[k] += -f_88 * gs_s_0[k]
                   + f_89 * gs_s_1[k]
                   - f_90 * gs_s_5[k]
                   + f_91 * gs_0[k]
                   - f_92 * gs_1[k]
                   + f_93 * gs_5[k]
                   + f_94 * pa_x[k] * hs_0[k]
                   - f_95 * pa_x[k] * hs_2[k]
                   + f_95 * pa_x[k] * hs_6[k]
                   - f_94 * pa_y[k] * hs_9[k]
                   + f_96 * is_s_0[k]
                   - f_97 * is_s_3[k]
                   + f_97 * is_s_10[k]
                   - f_96 * is_s_21[k];
    }
}

}  // namespace simdkin
