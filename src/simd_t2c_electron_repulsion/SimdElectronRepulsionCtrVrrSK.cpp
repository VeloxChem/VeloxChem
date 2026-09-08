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


#include "SimdElectronRepulsionCtrVrrSK.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ctr_sk_electron_repulsion_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                                    const size_t pb, const size_t sh0, const size_t sh1,
                                    const size_t si, const size_t ncols, const double alpha,
                                    const double beta, const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 1.640625 * std::sqrt(429.0) / beta;
    const auto f_1 = 0.234375 * std::sqrt(429.0) / beta;
    const auto f_2 = 1.640625 * std::sqrt(429.0) * alpha / (beta * p);
    const auto f_3 = 0.234375 * std::sqrt(429.0) * alpha / (beta * p);
    const auto f_4 = 0.21875 * std::sqrt(429.0);
    const auto f_5 = 1.09375 * std::sqrt(429.0);
    const auto f_6 = 0.65625 * std::sqrt(429.0);
    const auto f_7 = 0.03125 * std::sqrt(429.0);
    const auto f_8 = 0.1875 * std::sqrt(6006.0);
    const auto f_9 = 0.625 * std::sqrt(6006.0);
    const auto f_10 = 0.234375 * std::sqrt(231.0) / beta;
    const auto f_11 = 0.046875 * std::sqrt(231.0) / beta;
    const auto f_12 = 1.125 * std::sqrt(231.0) / beta;
    const auto f_13 = 0.234375 * std::sqrt(231.0) * alpha / (beta * p);
    const auto f_14 = 0.046875 * std::sqrt(231.0) * alpha / (beta * p);
    const auto f_15 = 1.125 * std::sqrt(231.0) * alpha / (beta * p);
    const auto f_16 = 0.15625 * std::sqrt(231.0);
    const auto f_17 = 1.875 * std::sqrt(231.0);
    const auto f_18 = 0.28125 * std::sqrt(231.0);
    const auto f_19 = 3.75 * std::sqrt(231.0);
    const auto f_20 = 0.03125 * std::sqrt(231.0);
    const auto f_21 = 0.375 * std::sqrt(231.0);
    const auto f_22 = 0.75 * std::sqrt(231.0);
    const auto f_23 = 2.5 * std::sqrt(231.0);
    const auto f_24 = 0.703125 * std::sqrt(21.0) / beta;
    const auto f_25 = 0.234375 * std::sqrt(21.0) / beta;
    const auto f_26 = 1.875 * std::sqrt(21.0) / beta;
    const auto f_27 = 2.5 * std::sqrt(21.0) / beta;
    const auto f_28 = 0.703125 * std::sqrt(21.0) * alpha / (beta * p);
    const auto f_29 = 0.234375 * std::sqrt(21.0) * alpha / (beta * p);
    const auto f_30 = 1.875 * std::sqrt(21.0) * alpha / (beta * p);
    const auto f_31 = 2.5 * std::sqrt(21.0) * alpha / (beta * p);
    const auto f_32 = 0.28125 * std::sqrt(21.0);
    const auto f_33 = 5.625 * std::sqrt(21.0);
    const auto f_34 = 0.46875 * std::sqrt(21.0);
    const auto f_35 = 7.5 * std::sqrt(21.0);
    const auto f_36 = 0.09375 * std::sqrt(21.0);
    const auto f_37 = 3.75 * std::sqrt(21.0);
    const auto f_38 = 1.875 * std::sqrt(21.0);
    const auto f_39 = 2.5 * std::sqrt(21.0);
    const auto f_40 = 0.9375 * std::sqrt(42.0);
    const auto f_41 = 1.875 * std::sqrt(42.0);
    const auto f_42 = 5.0 * std::sqrt(42.0);
    const auto f_43 = 3.0 * std::sqrt(42.0);
    const auto f_44 = 0.703125 * std::sqrt(7.0) / beta;
    const auto f_45 = 11.25 * std::sqrt(7.0) / beta;
    const auto f_46 = 7.5 * std::sqrt(7.0) / beta;
    const auto f_47 = 0.703125 * std::sqrt(7.0) * alpha / (beta * p);
    const auto f_48 = 11.25 * std::sqrt(7.0) * alpha / (beta * p);
    const auto f_49 = 7.5 * std::sqrt(7.0) * alpha / (beta * p);
    const auto f_50 = 0.15625 * std::sqrt(7.0);
    const auto f_51 = 3.75 * std::sqrt(7.0);
    const auto f_52 = 0.46875 * std::sqrt(7.0);
    const auto f_53 = 7.5 * std::sqrt(7.0);
    const auto f_54 = 2.0 * std::sqrt(7.0);
    const auto f_55 = 19.6875 / beta;
    const auto f_56 = 32.8125 / beta;
    const auto f_57 = 7.5 / beta;
    const auto f_58 = 19.6875 * alpha / (beta * p);
    const auto f_59 = 32.8125 * alpha / (beta * p);
    const auto f_60 = 7.5 * alpha / (beta * p);
    const auto f_61 = 0.46875 * std::sqrt(7.0) / beta;
    const auto f_62 = 0.9375 * std::sqrt(7.0) / beta;
    const auto f_63 = 0.46875 * std::sqrt(7.0) * alpha / (beta * p);
    const auto f_64 = 0.9375 * std::sqrt(7.0) * alpha / (beta * p);
    const auto f_65 = 3.75 * std::sqrt(42.0) / beta;
    const auto f_66 = 3.75 * std::sqrt(42.0) * alpha / (beta * p);
    const auto f_67 = 0.46875 * std::sqrt(42.0);
    const auto f_68 = 2.5 * std::sqrt(42.0);
    const auto f_69 = 1.5 * std::sqrt(42.0);
    const auto f_70 = 0.28125 * std::sqrt(21.0) / beta;
    const auto f_71 = 0.1875 * std::sqrt(21.0) / beta;
    const auto f_72 = 3.75 * std::sqrt(21.0) / beta;
    const auto f_73 = 0.46875 * std::sqrt(21.0) / beta;
    const auto f_74 = 0.28125 * std::sqrt(21.0) * alpha / (beta * p);
    const auto f_75 = 0.1875 * std::sqrt(21.0) * alpha / (beta * p);
    const auto f_76 = 3.75 * std::sqrt(21.0) * alpha / (beta * p);
    const auto f_77 = 0.46875 * std::sqrt(21.0) * alpha / (beta * p);
    const auto f_78 = 0.9375 * std::sqrt(231.0) / beta;
    const auto f_79 = 0.9375 * std::sqrt(231.0) * alpha / (beta * p);
    const auto f_80 = 0.1875 * std::sqrt(231.0);
    const auto f_81 = 0.9375 * std::sqrt(231.0);
    const auto f_82 = 0.625 * std::sqrt(231.0);
    const auto f_83 = 0.09375 * std::sqrt(231.0) / beta;
    const auto f_84 = 0.5625 * std::sqrt(231.0) / beta;
    const auto f_85 = 0.75 * std::sqrt(231.0) / beta;
    const auto f_86 = 0.15625 * std::sqrt(231.0) / beta;
    const auto f_87 = 3.75 * std::sqrt(231.0) / beta;
    const auto f_88 = 0.09375 * std::sqrt(231.0) * alpha / (beta * p);
    const auto f_89 = 0.5625 * std::sqrt(231.0) * alpha / (beta * p);
    const auto f_90 = 0.75 * std::sqrt(231.0) * alpha / (beta * p);
    const auto f_91 = 0.15625 * std::sqrt(231.0) * alpha / (beta * p);
    const auto f_92 = 3.75 * std::sqrt(231.0) * alpha / (beta * p);
    const auto f_93 = 0.03125 * std::sqrt(6006.0);
    const auto f_94 = 0.46875 * std::sqrt(6006.0);
    const auto f_95 = 0.09375 * std::sqrt(429.0) / beta;
    const auto f_96 = 1.3125 * std::sqrt(429.0) / beta;
    const auto f_97 = 1.09375 * std::sqrt(429.0) / beta;
    const auto f_98 = 0.09375 * std::sqrt(429.0) * alpha / (beta * p);
    const auto f_99 = 1.3125 * std::sqrt(429.0) * alpha / (beta * p);
    const auto f_100 = 1.09375 * std::sqrt(429.0) * alpha / (beta * p);

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
    auto *g_13 = values + 13 * nvalues;
    auto *g_14 = values + 14 * nvalues;

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh0_0 = buffer.data(sh0 + 0);
    const auto *sh0_1 = buffer.data(sh0 + 1);
    const auto *sh0_2 = buffer.data(sh0 + 2);
    const auto *sh0_3 = buffer.data(sh0 + 3);
    const auto *sh0_4 = buffer.data(sh0 + 4);
    const auto *sh0_5 = buffer.data(sh0 + 5);
    const auto *sh0_6 = buffer.data(sh0 + 6);
    const auto *sh0_7 = buffer.data(sh0 + 7);
    const auto *sh0_8 = buffer.data(sh0 + 8);
    const auto *sh0_9 = buffer.data(sh0 + 9);
    const auto *sh0_10 = buffer.data(sh0 + 10);
    const auto *sh0_11 = buffer.data(sh0 + 11);
    const auto *sh0_12 = buffer.data(sh0 + 12);

    const auto *sh1_0 = buffer.data(sh1 + 0);
    const auto *sh1_1 = buffer.data(sh1 + 1);
    const auto *sh1_2 = buffer.data(sh1 + 2);
    const auto *sh1_3 = buffer.data(sh1 + 3);
    const auto *sh1_4 = buffer.data(sh1 + 4);
    const auto *sh1_5 = buffer.data(sh1 + 5);
    const auto *sh1_6 = buffer.data(sh1 + 6);
    const auto *sh1_7 = buffer.data(sh1 + 7);
    const auto *sh1_8 = buffer.data(sh1 + 8);
    const auto *sh1_9 = buffer.data(sh1 + 9);
    const auto *sh1_10 = buffer.data(sh1 + 10);
    const auto *sh1_11 = buffer.data(sh1 + 11);
    const auto *sh1_12 = buffer.data(sh1 + 12);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_7 = buffer.data(si + 7);
    const auto *si_8 = buffer.data(si + 8);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_10 = buffer.data(si + 10);
    const auto *si_11 = buffer.data(si + 11);
    const auto *si_12 = buffer.data(si + 12);
    const auto *si_13 = buffer.data(si + 13);
    const auto *si_14 = buffer.data(si + 14);
    const auto *si_15 = buffer.data(si + 15);
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_19 = buffer.data(si + 19);

#pragma omp simd aligned(pb_x, pb_y, pb_z, sh0_3, sh0_8, sh1_3, sh1_8, si_0, si_1, si_4, si_9, \
                         si_13, si_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += -f_0 * sh0_3[k]
                  + f_1 * sh0_8[k]
                  + f_2 * sh1_3[k]
                  - f_3 * sh1_8[k]
                  + f_4 * pb_y[k] * si_0[k]
                  - f_5 * pb_x[k] * si_4[k]
                  + f_6 * pb_x[k] * si_9[k]
                  - f_7 * pb_y[k] * si_13[k];

        g_1[k] += f_8 * pb_y[k] * si_1[k]
                  - f_9 * pb_z[k] * si_4[k]
                  + f_8 * pb_x[k] * si_14[k];
    }

#pragma omp simd aligned(pb_x, pb_y, sh0_3, sh0_8, sh0_9, sh1_3, sh1_8, sh1_9, si_0, si_3, \
                         si_4, si_9, si_10, si_13, si_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_2[k] += f_10 * sh0_3[k]
                  + f_11 * sh0_8[k]
                  - f_12 * sh0_9[k]
                  - f_13 * sh1_3[k]
                  - f_14 * sh1_8[k]
                  + f_15 * sh1_9[k]
                  - f_16 * pb_y[k] * si_0[k]
                  + f_17 * pb_y[k] * si_3[k]
                  + f_16 * pb_x[k] * si_4[k]
                  + f_18 * pb_x[k] * si_9[k]
                  - f_19 * pb_x[k] * si_10[k]
                  - f_20 * pb_y[k] * si_13[k]
                  + f_21 * pb_y[k] * si_15[k];
    }

#pragma omp simd aligned(pb_x, pb_y, si_1, si_5, si_14, si_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += -f_22 * pb_y[k] * si_1[k]
                  + f_23 * pb_y[k] * si_5[k]
                  + f_22 * pb_x[k] * si_14[k]
                  - f_23 * pb_x[k] * si_16[k];
    }

#pragma omp simd aligned(pb_x, pb_y, sh0_3, sh0_8, sh0_9, sh0_11, sh1_3, sh1_8, sh1_9, sh1_11, \
                         si_0, si_3, si_4, si_8, si_9, si_10, si_13, si_15, \
                         si_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += f_24 * sh0_3[k]
                  - f_25 * sh0_8[k]
                  + f_26 * sh0_9[k]
                  - f_27 * sh0_11[k]
                  - f_28 * sh1_3[k]
                  + f_29 * sh1_8[k]
                  - f_30 * sh1_9[k]
                  + f_31 * sh1_11[k]
                  + f_32 * pb_y[k] * si_0[k]
                  - f_33 * pb_y[k] * si_3[k]
                  + f_34 * pb_x[k] * si_4[k]
                  + f_35 * pb_y[k] * si_8[k]
                  + f_36 * pb_x[k] * si_9[k]
                  - f_37 * pb_x[k] * si_10[k]
                  - f_36 * pb_y[k] * si_13[k]
                  + f_38 * pb_y[k] * si_15[k]
                  - f_39 * pb_y[k] * si_17[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, si_1, si_4, si_5, si_14, si_16, \
                         si_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += f_40 * pb_y[k] * si_1[k]
                  + f_41 * pb_z[k] * si_4[k]
                  - f_42 * pb_y[k] * si_5[k]
                  + f_40 * pb_x[k] * si_14[k]
                  - f_42 * pb_x[k] * si_16[k]
                  + f_43 * pb_x[k] * si_18[k];
    }

#pragma omp simd aligned(pb_x, pb_y, sh0_3, sh0_8, sh0_9, sh0_11, sh1_3, sh1_8, sh1_9, sh1_11, \
                         si_0, si_3, si_4, si_8, si_9, si_10, si_13, si_15, si_17, \
                         si_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += -f_44 * sh0_3[k]
                  - f_44 * sh0_8[k]
                  + f_45 * sh0_9[k]
                  - f_46 * sh0_11[k]
                  + f_47 * sh1_3[k]
                  + f_47 * sh1_8[k]
                  - f_48 * sh1_9[k]
                  + f_49 * sh1_11[k]
                  - f_50 * pb_y[k] * si_0[k]
                  + f_51 * pb_y[k] * si_3[k]
                  - f_52 * pb_x[k] * si_4[k]
                  - f_53 * pb_y[k] * si_8[k]
                  - f_52 * pb_x[k] * si_9[k]
                  + f_53 * pb_x[k] * si_10[k]
                  - f_50 * pb_y[k] * si_13[k]
                  + f_51 * pb_y[k] * si_15[k]
                  - f_53 * pb_y[k] * si_17[k]
                  + f_54 * pb_y[k] * si_19[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, sh0_4, sh0_10, sh0_12, sh1_4, sh1_10, sh1_12, si_0, \
                         si_2, si_5, si_6, si_11, si_12, si_13, si_16, si_18, \
                         si_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_7[k] += f_55 * sh0_4[k]
                  + f_56 * sh0_10[k]
                  - f_57 * sh0_12[k]
                  - f_58 * sh1_4[k]
                  - f_59 * sh1_10[k]
                  + f_60 * sh1_12[k]
                  - 2.1875 * pb_z[k] * si_0[k]
                  - 6.5625 * pb_z[k] * si_2[k]
                  + 13.125 * pb_x[k] * si_5[k]
                  - 6.5625 * pb_z[k] * si_6[k]
                  + 26.25 * pb_x[k] * si_11[k]
                  - 10.5 * pb_x[k] * si_12[k]
                  - 2.1875 * pb_z[k] * si_13[k]
                  + 13.125 * pb_y[k] * si_16[k]
                  - 10.5 * pb_y[k] * si_18[k]
                  + pb_z[k] * si_19[k];
    }

#pragma omp simd aligned(pb_x, sh0_0, sh0_1, sh0_2, sh0_5, sh0_6, sh0_7, sh1_0, sh1_1, sh1_2, \
                         sh1_5, sh1_6, sh1_7, si_0, si_2, si_3, si_6, si_7, si_8, si_13, \
                         si_15, si_17, si_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += -f_61 * sh0_0[k]
                  - f_62 * sh0_1[k]
                  + f_46 * sh0_2[k]
                  - f_61 * sh0_5[k]
                  + f_46 * sh0_6[k]
                  - f_46 * sh0_7[k]
                  + f_63 * sh1_0[k]
                  + f_64 * sh1_1[k]
                  - f_49 * sh1_2[k]
                  + f_63 * sh1_5[k]
                  - f_49 * sh1_6[k]
                  + f_49 * sh1_7[k]
                  - f_50 * pb_x[k] * si_0[k]
                  - f_52 * pb_x[k] * si_2[k]
                  + f_51 * pb_x[k] * si_3[k]
                  - f_52 * pb_x[k] * si_6[k]
                  + f_53 * pb_x[k] * si_7[k]
                  - f_53 * pb_x[k] * si_8[k]
                  - f_50 * pb_x[k] * si_13[k]
                  + f_51 * pb_x[k] * si_15[k]
                  - f_53 * pb_x[k] * si_17[k]
                  + f_54 * pb_x[k] * si_19[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, sh0_4, sh0_10, sh1_4, sh1_10, si_0, si_2, si_5, \
                         si_6, si_12, si_13, si_16, si_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_9[k] += -f_65 * sh0_4[k]
                  + f_65 * sh0_10[k]
                  + f_66 * sh1_4[k]
                  - f_66 * sh1_10[k]
                  + f_67 * pb_z[k] * si_0[k]
                  + f_67 * pb_z[k] * si_2[k]
                  - f_68 * pb_x[k] * si_5[k]
                  - f_67 * pb_z[k] * si_6[k]
                  + f_69 * pb_x[k] * si_12[k]
                  - f_67 * pb_z[k] * si_13[k]
                  + f_68 * pb_y[k] * si_16[k]
                  - f_69 * pb_y[k] * si_18[k];
    }

#pragma omp simd aligned(pb_x, sh0_0, sh0_1, sh0_2, sh0_5, sh0_6, sh0_7, sh1_0, sh1_1, sh1_2, \
                         sh1_5, sh1_6, sh1_7, si_0, si_2, si_3, si_6, si_7, si_8, si_13, \
                         si_15, si_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_10[k] += f_70 * sh0_0[k]
                   - f_71 * sh0_1[k]
                   - f_72 * sh0_2[k]
                   - f_73 * sh0_5[k]
                   + f_72 * sh0_6[k]
                   + f_27 * sh0_7[k]
                   - f_74 * sh1_0[k]
                   + f_75 * sh1_1[k]
                   + f_76 * sh1_2[k]
                   + f_77 * sh1_5[k]
                   - f_76 * sh1_6[k]
                   - f_31 * sh1_7[k]
                   + f_36 * pb_x[k] * si_0[k]
                   - f_36 * pb_x[k] * si_2[k]
                   - f_38 * pb_x[k] * si_3[k]
                   - f_34 * pb_x[k] * si_6[k]
                   + f_37 * pb_x[k] * si_7[k]
                   + f_39 * pb_x[k] * si_8[k]
                   - f_32 * pb_x[k] * si_13[k]
                   + f_33 * pb_x[k] * si_15[k]
                   - f_35 * pb_x[k] * si_17[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, sh0_4, sh0_10, sh1_4, sh1_10, si_0, si_2, si_5, \
                         si_6, si_11, si_13, si_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_11[k] += f_78 * sh0_4[k]
                   - f_78 * sh0_10[k]
                   - f_79 * sh1_4[k]
                   + f_79 * sh1_10[k]
                   - f_80 * pb_z[k] * si_0[k]
                   + f_81 * pb_z[k] * si_2[k]
                   + f_82 * pb_x[k] * si_5[k]
                   + f_81 * pb_z[k] * si_6[k]
                   - f_19 * pb_x[k] * si_11[k]
                   - f_80 * pb_z[k] * si_13[k]
                   + f_82 * pb_y[k] * si_16[k];
    }

#pragma omp simd aligned(pb_x, sh0_0, sh0_1, sh0_2, sh0_5, sh0_6, sh1_0, sh1_1, sh1_2, sh1_5, \
                         sh1_6, si_0, si_2, si_3, si_6, si_7, si_13, \
                         si_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_12[k] += -f_83 * sh0_0[k]
                   + f_84 * sh0_1[k]
                   + f_85 * sh0_2[k]
                   + f_86 * sh0_5[k]
                   - f_87 * sh0_6[k]
                   + f_88 * sh1_0[k]
                   - f_89 * sh1_1[k]
                   - f_90 * sh1_2[k]
                   - f_91 * sh1_5[k]
                   + f_92 * sh1_6[k]
                   - f_20 * pb_x[k] * si_0[k]
                   + f_18 * pb_x[k] * si_2[k]
                   + f_21 * pb_x[k] * si_3[k]
                   + f_16 * pb_x[k] * si_6[k]
                   - f_19 * pb_x[k] * si_7[k]
                   - f_16 * pb_x[k] * si_13[k]
                   + f_17 * pb_x[k] * si_15[k];
    }

#pragma omp simd aligned(pb_x, pb_z, sh0_0, sh0_1, sh0_5, sh1_0, sh1_1, sh1_5, si_0, si_2, \
                         si_6, si_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_13[k] += f_93 * pb_z[k] * si_0[k]
                   - f_94 * pb_z[k] * si_2[k]
                   + f_94 * pb_z[k] * si_6[k]
                   - f_93 * pb_z[k] * si_13[k];

        g_14[k] += f_95 * sh0_0[k]
                   - f_96 * sh0_1[k]
                   + f_97 * sh0_5[k]
                   - f_98 * sh1_0[k]
                   + f_99 * sh1_1[k]
                   - f_100 * sh1_5[k]
                   + f_7 * pb_x[k] * si_0[k]
                   - f_6 * pb_x[k] * si_2[k]
                   + f_5 * pb_x[k] * si_6[k]
                   - f_4 * pb_x[k] * si_13[k];
    }
}

}  // namespace simdt2ceri
