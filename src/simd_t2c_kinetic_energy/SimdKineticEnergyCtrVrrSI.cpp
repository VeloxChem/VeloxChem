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


#include "SimdKineticEnergyCtrVrrSI.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_ctr_si_kinetic_energy_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                                const size_t pb, const size_t sg_s, const size_t si_s,
                                const size_t sg, const size_t sh, const size_t ncols,
                                const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 1.25 * std::sqrt(462.0) * alpha / p;
    const auto f_1 = 0.375 * std::sqrt(462.0) * alpha * beta / p;
    const auto f_2 = 1.25 * std::sqrt(462.0) * alpha * beta / p;
    const auto f_3 = 0.625 * std::sqrt(462.0) / p;
    const auto f_4 = 0.1875 * std::sqrt(462.0);
    const auto f_5 = 0.625 * std::sqrt(462.0);
    const auto f_6 = 1.875 * std::sqrt(154.0) * alpha * beta / p;
    const auto f_7 = 3.75 * std::sqrt(154.0) * alpha * beta / p;
    const auto f_8 = 0.375 * std::sqrt(154.0) * alpha * beta / p;
    const auto f_9 = 0.9375 * std::sqrt(154.0);
    const auto f_10 = 1.875 * std::sqrt(154.0);
    const auto f_11 = 0.1875 * std::sqrt(154.0);
    const auto f_12 = 1.5 * std::sqrt(7.0) * alpha * beta / p;
    const auto f_13 = 15.0 * std::sqrt(7.0) * alpha * beta / p;
    const auto f_14 = 0.75 * std::sqrt(7.0);
    const auto f_15 = 7.5 * std::sqrt(7.0);
    const auto f_16 = std::sqrt(210.0) * alpha / p;
    const auto f_17 = 1.125 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_18 = 0.75 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_19 = 3.0 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_20 = 0.375 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_21 = std::sqrt(210.0) * alpha * beta / p;
    const auto f_22 = 0.5 * std::sqrt(210.0) / p;
    const auto f_23 = 0.5625 * std::sqrt(210.0);
    const auto f_24 = 0.375 * std::sqrt(210.0);
    const auto f_25 = 1.5 * std::sqrt(210.0);
    const auto f_26 = 0.1875 * std::sqrt(210.0);
    const auto f_27 = 0.5 * std::sqrt(210.0);
    const auto f_28 = 0.25 * std::sqrt(210.0) * alpha / p;
    const auto f_29 = 0.125 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_30 = 0.25 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_31 = 2.0 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_32 = 0.125 * std::sqrt(210.0) / p;
    const auto f_33 = 0.0625 * std::sqrt(210.0);
    const auto f_34 = std::sqrt(210.0);
    const auto f_35 = 0.125 * std::sqrt(210.0);
    const auto f_36 = 5.0 * std::sqrt(21.0) * alpha / p;
    const auto f_37 = 1.25 * std::sqrt(21.0) * alpha * beta / p;
    const auto f_38 = 2.5 * std::sqrt(21.0) * alpha * beta / p;
    const auto f_39 = 5.0 * std::sqrt(21.0) * alpha * beta / p;
    const auto f_40 = 2.0 * std::sqrt(21.0) * alpha * beta / p;
    const auto f_41 = 2.5 * std::sqrt(21.0) / p;
    const auto f_42 = 0.625 * std::sqrt(21.0);
    const auto f_43 = 1.25 * std::sqrt(21.0);
    const auto f_44 = 2.5 * std::sqrt(21.0);
    const auto f_45 = std::sqrt(21.0);
    const auto f_46 = 1.5625 * alpha / p;
    const auto f_47 = 2.8125 * alpha / p;
    const auto f_48 = 16.875 * alpha / p;
    const auto f_49 = 2.5 * alpha / p;
    const auto f_50 = 28.125 * alpha / p;
    const auto f_51 = 10.0 * alpha / p;
    const auto f_52 = 0.625 * alpha * beta / p;
    const auto f_53 = 1.875 * alpha * beta / p;
    const auto f_54 = 11.25 * alpha * beta / p;
    const auto f_55 = 22.5 * alpha * beta / p;
    const auto f_56 = 15.0 * alpha * beta / p;
    const auto f_57 = 2.0 * alpha * beta / p;
    const auto f_58 = 0.78125 / p;
    const auto f_59 = 1.40625 / p;
    const auto f_60 = 8.4375 / p;
    const auto f_61 = 1.25 / p;
    const auto f_62 = 14.0625 / p;
    const auto f_63 = 5.0 / p;
    const auto f_64 = 0.15625 * std::sqrt(210.0) * alpha / p;
    const auto f_65 = 0.09375 * std::sqrt(210.0) * alpha / p;
    const auto f_66 = 1.5 * std::sqrt(210.0) * alpha / p;
    const auto f_67 = 0.1875 * std::sqrt(210.0) * alpha / p;
    const auto f_68 = 0.0625 * std::sqrt(210.0) * alpha * beta / p;
    const auto f_69 = 0.078125 * std::sqrt(210.0) / p;
    const auto f_70 = 0.046875 * std::sqrt(210.0) / p;
    const auto f_71 = 0.75 * std::sqrt(210.0) / p;
    const auto f_72 = 0.09375 * std::sqrt(210.0) / p;
    const auto f_73 = 0.03125 * std::sqrt(210.0);
    const auto f_74 = 0.9375 * std::sqrt(7.0) * alpha / p;
    const auto f_75 = 2.8125 * std::sqrt(7.0) * alpha / p;
    const auto f_76 = 5.625 * std::sqrt(7.0) * alpha / p;
    const auto f_77 = 0.375 * std::sqrt(7.0) * alpha * beta / p;
    const auto f_78 = 1.875 * std::sqrt(7.0) * alpha * beta / p;
    const auto f_79 = 3.75 * std::sqrt(7.0) * alpha * beta / p;
    const auto f_80 = 22.5 * std::sqrt(7.0) * alpha * beta / p;
    const auto f_81 = 0.46875 * std::sqrt(7.0) / p;
    const auto f_82 = 1.40625 * std::sqrt(7.0) / p;
    const auto f_83 = 2.8125 * std::sqrt(7.0) / p;
    const auto f_84 = 0.1875 * std::sqrt(7.0);
    const auto f_85 = 0.9375 * std::sqrt(7.0);
    const auto f_86 = 1.875 * std::sqrt(7.0);
    const auto f_87 = 11.25 * std::sqrt(7.0);
    const auto f_88 = 0.15625 * std::sqrt(462.0) * alpha / p;
    const auto f_89 = 1.40625 * std::sqrt(462.0) * alpha / p;
    const auto f_90 = 0.3125 * std::sqrt(462.0) * alpha / p;
    const auto f_91 = 0.0625 * std::sqrt(462.0) * alpha * beta / p;
    const auto f_92 = 0.9375 * std::sqrt(462.0) * alpha * beta / p;
    const auto f_93 = 0.078125 * std::sqrt(462.0) / p;
    const auto f_94 = 0.703125 * std::sqrt(462.0) / p;
    const auto f_95 = 0.15625 * std::sqrt(462.0) / p;
    const auto f_96 = 0.03125 * std::sqrt(462.0);
    const auto f_97 = 0.46875 * std::sqrt(462.0);

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

    const auto *sg_s_0 = buffer.data(sg_s + 0);
    const auto *sg_s_1 = buffer.data(sg_s + 1);
    const auto *sg_s_2 = buffer.data(sg_s + 2);
    const auto *sg_s_3 = buffer.data(sg_s + 3);
    const auto *sg_s_4 = buffer.data(sg_s + 4);
    const auto *sg_s_5 = buffer.data(sg_s + 5);
    const auto *sg_s_6 = buffer.data(sg_s + 6);
    const auto *sg_s_7 = buffer.data(sg_s + 7);
    const auto *sg_s_8 = buffer.data(sg_s + 8);

    const auto *si_s_0 = buffer.data(si_s + 0);
    const auto *si_s_1 = buffer.data(si_s + 1);
    const auto *si_s_2 = buffer.data(si_s + 2);
    const auto *si_s_3 = buffer.data(si_s + 3);
    const auto *si_s_4 = buffer.data(si_s + 4);
    const auto *si_s_5 = buffer.data(si_s + 5);
    const auto *si_s_6 = buffer.data(si_s + 6);
    const auto *si_s_7 = buffer.data(si_s + 7);
    const auto *si_s_8 = buffer.data(si_s + 8);
    const auto *si_s_9 = buffer.data(si_s + 9);
    const auto *si_s_10 = buffer.data(si_s + 10);
    const auto *si_s_11 = buffer.data(si_s + 11);
    const auto *si_s_12 = buffer.data(si_s + 12);
    const auto *si_s_13 = buffer.data(si_s + 13);
    const auto *si_s_14 = buffer.data(si_s + 14);
    const auto *si_s_15 = buffer.data(si_s + 15);
    const auto *si_s_16 = buffer.data(si_s + 16);
    const auto *si_s_17 = buffer.data(si_s + 17);
    const auto *si_s_18 = buffer.data(si_s + 18);
    const auto *si_s_19 = buffer.data(si_s + 19);
    const auto *si_s_20 = buffer.data(si_s + 20);
    const auto *si_s_21 = buffer.data(si_s + 21);
    const auto *si_s_22 = buffer.data(si_s + 22);
    const auto *si_s_23 = buffer.data(si_s + 23);
    const auto *si_s_24 = buffer.data(si_s + 24);
    const auto *si_s_25 = buffer.data(si_s + 25);
    const auto *si_s_26 = buffer.data(si_s + 26);
    const auto *si_s_27 = buffer.data(si_s + 27);

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

#pragma omp simd aligned(pb_x, pb_y, sg_s_3, si_s_1, si_s_6, si_s_15, sg_3, sh_0, sh_4, \
                         sh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += f_0 * sg_s_3[k]
                  + f_1 * si_s_1[k]
                  - f_2 * si_s_6[k]
                  + f_1 * si_s_15[k]
                  - f_3 * sg_3[k]
                  + f_4 * pb_y[k] * sh_0[k]
                  - f_5 * pb_x[k] * sh_4[k]
                  + f_4 * pb_x[k] * sh_9[k];
    }

#pragma omp simd aligned(pb_y, pb_z, si_s_4, si_s_11, si_s_22, sh_1, sh_4, \
                         sh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_1[k] += f_6 * si_s_4[k]
                  - f_7 * si_s_11[k]
                  + f_8 * si_s_22[k]
                  + f_9 * pb_y[k] * sh_1[k]
                  - f_10 * pb_z[k] * sh_4[k]
                  + f_11 * pb_z[k] * sh_9[k];
    }

#pragma omp simd aligned(pb_x, pb_y, si_s_1, si_s_8, si_s_15, si_s_17, sh_0, sh_3, sh_9, \
                         sh_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_2[k] += -f_12 * si_s_1[k]
                  + f_13 * si_s_8[k]
                  + f_12 * si_s_15[k]
                  - f_13 * si_s_17[k]
                  - f_14 * pb_y[k] * sh_0[k]
                  + f_15 * pb_y[k] * sh_3[k]
                  + f_14 * pb_x[k] * sh_9[k]
                  - f_15 * pb_x[k] * sh_11[k];
    }

#pragma omp simd aligned(pb_y, pb_z, sg_s_7, si_s_4, si_s_11, si_s_13, si_s_22, si_s_24, sg_7, \
                         sh_1, sh_4, sh_5, sh_9, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_3[k] += f_16 * sg_s_7[k]
                  - f_17 * si_s_4[k]
                  - f_18 * si_s_11[k]
                  + f_19 * si_s_13[k]
                  + f_20 * si_s_22[k]
                  - f_21 * si_s_24[k]
                  - f_22 * sg_7[k]
                  - f_23 * pb_y[k] * sh_1[k]
                  - f_24 * pb_z[k] * sh_4[k]
                  + f_25 * pb_y[k] * sh_5[k]
                  + f_26 * pb_z[k] * sh_9[k]
                  - f_27 * pb_y[k] * sh_12[k];
    }

#pragma omp simd aligned(pb_x, pb_y, sg_s_3, si_s_1, si_s_6, si_s_8, si_s_15, si_s_17, \
                         si_s_19, sg_3, sh_0, sh_3, sh_4, sh_9, sh_11, \
                         sh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_4[k] += -f_28 * sg_s_3[k]
                  + f_29 * si_s_1[k]
                  + f_30 * si_s_6[k]
                  - f_31 * si_s_8[k]
                  + f_29 * si_s_15[k]
                  - f_31 * si_s_17[k]
                  + f_31 * si_s_19[k]
                  + f_32 * sg_3[k]
                  + f_33 * pb_y[k] * sh_0[k]
                  - f_34 * pb_y[k] * sh_3[k]
                  + f_35 * pb_x[k] * sh_4[k]
                  + f_33 * pb_x[k] * sh_9[k]
                  - f_34 * pb_x[k] * sh_11[k]
                  + f_34 * pb_x[k] * sh_13[k];
    }

#pragma omp simd aligned(pb_y, pb_z, sg_s_7, si_s_4, si_s_11, si_s_13, si_s_22, si_s_24, \
                         si_s_26, sg_7, sh_1, sh_4, sh_5, sh_9, sh_12, \
                         sh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_5[k] += f_36 * sg_s_7[k]
                  + f_37 * si_s_4[k]
                  + f_38 * si_s_11[k]
                  - f_39 * si_s_13[k]
                  + f_37 * si_s_22[k]
                  - f_39 * si_s_24[k]
                  + f_40 * si_s_26[k]
                  - f_41 * sg_7[k]
                  + f_42 * pb_y[k] * sh_1[k]
                  + f_43 * pb_z[k] * sh_4[k]
                  - f_44 * pb_y[k] * sh_5[k]
                  + f_42 * pb_z[k] * sh_9[k]
                  - f_44 * pb_y[k] * sh_12[k]
                  + f_45 * pb_y[k] * sh_14[k];
    }

#pragma omp simd aligned(pb_x, pb_y, pb_z, sg_s_0, sg_s_1, sg_s_2, sg_s_5, sg_s_6, sg_s_8, \
                         si_s_0, si_s_3, si_s_5, si_s_10, si_s_12, si_s_14, si_s_21, si_s_23, \
                         si_s_25, si_s_27, sg_0, sg_1, sg_2, sg_5, sg_6, sg_8, sh_0, sh_2, \
                         sh_3, sh_6, sh_7, sh_8, sh_9, sh_11, sh_13, \
                         sh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += f_46 * sg_s_0[k]
                  + f_47 * sg_s_1[k]
                  - f_48 * sg_s_2[k]
                  + f_49 * sg_s_5[k]
                  - f_50 * sg_s_6[k]
                  + f_51 * sg_s_8[k]
                  - f_52 * si_s_0[k]
                  - f_53 * si_s_3[k]
                  + f_54 * si_s_5[k]
                  - f_53 * si_s_10[k]
                  + f_55 * si_s_12[k]
                  - f_56 * si_s_14[k]
                  - f_52 * si_s_21[k]
                  + f_54 * si_s_23[k]
                  - f_56 * si_s_25[k]
                  + f_57 * si_s_27[k]
                  - f_58 * sg_0[k]
                  - f_59 * sg_1[k]
                  + f_60 * sg_2[k]
                  - f_61 * sg_5[k]
                  + f_62 * sg_6[k]
                  - f_63 * sg_8[k]
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

#pragma omp simd aligned(pb_x, pb_z, sg_s_4, si_s_2, si_s_7, si_s_9, si_s_16, si_s_18, \
                         si_s_20, sg_4, sh_0, sh_2, sh_5, sh_10, sh_12, \
                         sh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_7[k] += f_36 * sg_s_4[k]
                  + f_37 * si_s_2[k]
                  + f_38 * si_s_7[k]
                  - f_39 * si_s_9[k]
                  + f_37 * si_s_16[k]
                  - f_39 * si_s_18[k]
                  + f_40 * si_s_20[k]
                  - f_41 * sg_4[k]
                  + f_42 * pb_z[k] * sh_0[k]
                  + f_43 * pb_z[k] * sh_2[k]
                  - f_44 * pb_x[k] * sh_5[k]
                  + f_42 * pb_x[k] * sh_10[k]
                  - f_44 * pb_x[k] * sh_12[k]
                  + f_45 * pb_x[k] * sh_14[k];
    }

#pragma omp simd aligned(pb_x, pb_y, sg_s_0, sg_s_1, sg_s_2, sg_s_5, sg_s_6, si_s_0, si_s_3, \
                         si_s_5, si_s_10, si_s_14, si_s_21, si_s_23, si_s_25, sg_0, sg_1, \
                         sg_2, sg_5, sg_6, sh_0, sh_2, sh_3, sh_6, sh_8, sh_9, sh_11, \
                         sh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_8[k] += -f_64 * sg_s_0[k]
                  - f_65 * sg_s_1[k]
                  + f_66 * sg_s_2[k]
                  + f_67 * sg_s_5[k]
                  - f_66 * sg_s_6[k]
                  + f_68 * si_s_0[k]
                  + f_68 * si_s_3[k]
                  - f_21 * si_s_5[k]
                  - f_68 * si_s_10[k]
                  + f_21 * si_s_14[k]
                  - f_68 * si_s_21[k]
                  + f_21 * si_s_23[k]
                  - f_21 * si_s_25[k]
                  + f_69 * sg_0[k]
                  + f_70 * sg_1[k]
                  - f_71 * sg_2[k]
                  - f_72 * sg_5[k]
                  + f_71 * sg_6[k]
                  + f_73 * pb_x[k] * sh_0[k]
                  + f_73 * pb_x[k] * sh_2[k]
                  - f_27 * pb_x[k] * sh_3[k]
                  - f_73 * pb_x[k] * sh_6[k]
                  + f_27 * pb_x[k] * sh_8[k]
                  - f_73 * pb_y[k] * sh_9[k]
                  + f_27 * pb_y[k] * sh_11[k]
                  - f_27 * pb_y[k] * sh_13[k];
    }

#pragma omp simd aligned(pb_x, pb_z, sg_s_4, si_s_2, si_s_7, si_s_9, si_s_16, si_s_18, sg_4, \
                         sh_0, sh_2, sh_5, sh_10, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_9[k] += -f_16 * sg_s_4[k]
                  - f_20 * si_s_2[k]
                  + f_18 * si_s_7[k]
                  + f_21 * si_s_9[k]
                  + f_17 * si_s_16[k]
                  - f_19 * si_s_18[k]
                  + f_22 * sg_4[k]
                  - f_26 * pb_z[k] * sh_0[k]
                  + f_24 * pb_z[k] * sh_2[k]
                  + f_27 * pb_x[k] * sh_5[k]
                  + f_23 * pb_x[k] * sh_10[k]
                  - f_25 * pb_x[k] * sh_12[k];
    }

#pragma omp simd aligned(pb_x, pb_y, sg_s_0, sg_s_1, sg_s_2, sg_s_6, si_s_0, si_s_3, si_s_5, \
                         si_s_10, si_s_12, si_s_21, si_s_23, sg_0, sg_1, sg_2, sg_6, sh_0, \
                         sh_2, sh_3, sh_6, sh_7, sh_9, sh_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_10[k] += f_74 * sg_s_0[k]
                   - f_75 * sg_s_1[k]
                   - f_76 * sg_s_2[k]
                   + f_76 * sg_s_6[k]
                   - f_77 * si_s_0[k]
                   + f_78 * si_s_3[k]
                   + f_79 * si_s_5[k]
                   + f_78 * si_s_10[k]
                   - f_80 * si_s_12[k]
                   - f_77 * si_s_21[k]
                   + f_79 * si_s_23[k]
                   - f_81 * sg_0[k]
                   + f_82 * sg_1[k]
                   + f_83 * sg_2[k]
                   - f_83 * sg_6[k]
                   - f_84 * pb_x[k] * sh_0[k]
                   + f_85 * pb_x[k] * sh_2[k]
                   + f_86 * pb_x[k] * sh_3[k]
                   + f_85 * pb_x[k] * sh_6[k]
                   - f_87 * pb_x[k] * sh_7[k]
                   - f_84 * pb_y[k] * sh_9[k]
                   + f_86 * pb_y[k] * sh_11[k];
    }

#pragma omp simd aligned(pb_x, pb_z, si_s_2, si_s_7, si_s_16, sh_0, sh_2, \
                         sh_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_11[k] += f_8 * si_s_2[k]
                   - f_7 * si_s_7[k]
                   + f_6 * si_s_16[k]
                   + f_11 * pb_z[k] * sh_0[k]
                   - f_10 * pb_z[k] * sh_2[k]
                   + f_9 * pb_x[k] * sh_10[k];
    }

#pragma omp simd aligned(pb_x, pb_y, sg_s_0, sg_s_1, sg_s_5, si_s_0, si_s_3, si_s_10, si_s_21, \
                         sg_0, sg_1, sg_5, sh_0, sh_2, sh_6, sh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_12[k] += -f_88 * sg_s_0[k]
                   + f_89 * sg_s_1[k]
                   - f_90 * sg_s_5[k]
                   + f_91 * si_s_0[k]
                   - f_92 * si_s_3[k]
                   + f_92 * si_s_10[k]
                   - f_91 * si_s_21[k]
                   + f_93 * sg_0[k]
                   - f_94 * sg_1[k]
                   + f_95 * sg_5[k]
                   + f_96 * pb_x[k] * sh_0[k]
                   - f_97 * pb_x[k] * sh_2[k]
                   + f_97 * pb_x[k] * sh_6[k]
                   - f_96 * pb_y[k] * sh_9[k];
    }
}

}  // namespace simdkin
