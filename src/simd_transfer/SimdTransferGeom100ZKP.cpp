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


#include "SimdTransferGeom100ZKP.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_100z_kp_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                      const size_t target, const size_t ks_1, const size_t ks_0,
                                      const size_t ls_1, const size_t ncomps,
                                      const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_0 = buffer.data(target + 0 * ncomps + c);
        auto *t_1 = buffer.data(target + 1 * ncomps + c);
        auto *t_2 = buffer.data(target + 2 * ncomps + c);
        auto *t_3 = buffer.data(target + 3 * ncomps + c);
        auto *t_4 = buffer.data(target + 4 * ncomps + c);
        auto *t_5 = buffer.data(target + 5 * ncomps + c);
        auto *t_6 = buffer.data(target + 6 * ncomps + c);
        auto *t_7 = buffer.data(target + 7 * ncomps + c);
        auto *t_8 = buffer.data(target + 8 * ncomps + c);
        auto *t_9 = buffer.data(target + 9 * ncomps + c);
        auto *t_10 = buffer.data(target + 10 * ncomps + c);
        auto *t_11 = buffer.data(target + 11 * ncomps + c);
        auto *t_12 = buffer.data(target + 12 * ncomps + c);
        auto *t_13 = buffer.data(target + 13 * ncomps + c);
        auto *t_14 = buffer.data(target + 14 * ncomps + c);
        auto *t_15 = buffer.data(target + 15 * ncomps + c);
        auto *t_16 = buffer.data(target + 16 * ncomps + c);
        auto *t_17 = buffer.data(target + 17 * ncomps + c);
        auto *t_18 = buffer.data(target + 18 * ncomps + c);
        auto *t_19 = buffer.data(target + 19 * ncomps + c);
        auto *t_20 = buffer.data(target + 20 * ncomps + c);
        auto *t_21 = buffer.data(target + 21 * ncomps + c);
        auto *t_22 = buffer.data(target + 22 * ncomps + c);
        auto *t_23 = buffer.data(target + 23 * ncomps + c);
        auto *t_24 = buffer.data(target + 24 * ncomps + c);
        auto *t_25 = buffer.data(target + 25 * ncomps + c);
        auto *t_26 = buffer.data(target + 26 * ncomps + c);
        auto *t_27 = buffer.data(target + 27 * ncomps + c);
        auto *t_28 = buffer.data(target + 28 * ncomps + c);
        auto *t_29 = buffer.data(target + 29 * ncomps + c);
        auto *t_30 = buffer.data(target + 30 * ncomps + c);
        auto *t_31 = buffer.data(target + 31 * ncomps + c);
        auto *t_32 = buffer.data(target + 32 * ncomps + c);
        auto *t_33 = buffer.data(target + 33 * ncomps + c);
        auto *t_34 = buffer.data(target + 34 * ncomps + c);
        auto *t_35 = buffer.data(target + 35 * ncomps + c);
        auto *t_36 = buffer.data(target + 36 * ncomps + c);
        auto *t_37 = buffer.data(target + 37 * ncomps + c);
        auto *t_38 = buffer.data(target + 38 * ncomps + c);
        auto *t_39 = buffer.data(target + 39 * ncomps + c);
        auto *t_40 = buffer.data(target + 40 * ncomps + c);
        auto *t_41 = buffer.data(target + 41 * ncomps + c);
        auto *t_42 = buffer.data(target + 42 * ncomps + c);
        auto *t_43 = buffer.data(target + 43 * ncomps + c);
        auto *t_44 = buffer.data(target + 44 * ncomps + c);
        auto *t_45 = buffer.data(target + 45 * ncomps + c);
        auto *t_46 = buffer.data(target + 46 * ncomps + c);
        auto *t_47 = buffer.data(target + 47 * ncomps + c);
        auto *t_48 = buffer.data(target + 48 * ncomps + c);
        auto *t_49 = buffer.data(target + 49 * ncomps + c);
        auto *t_50 = buffer.data(target + 50 * ncomps + c);
        auto *t_51 = buffer.data(target + 51 * ncomps + c);
        auto *t_52 = buffer.data(target + 52 * ncomps + c);
        auto *t_53 = buffer.data(target + 53 * ncomps + c);
        auto *t_54 = buffer.data(target + 54 * ncomps + c);
        auto *t_55 = buffer.data(target + 55 * ncomps + c);
        auto *t_56 = buffer.data(target + 56 * ncomps + c);
        auto *t_57 = buffer.data(target + 57 * ncomps + c);
        auto *t_58 = buffer.data(target + 58 * ncomps + c);
        auto *t_59 = buffer.data(target + 59 * ncomps + c);
        auto *t_60 = buffer.data(target + 60 * ncomps + c);
        auto *t_61 = buffer.data(target + 61 * ncomps + c);
        auto *t_62 = buffer.data(target + 62 * ncomps + c);
        auto *t_63 = buffer.data(target + 63 * ncomps + c);
        auto *t_64 = buffer.data(target + 64 * ncomps + c);
        auto *t_65 = buffer.data(target + 65 * ncomps + c);
        auto *t_66 = buffer.data(target + 66 * ncomps + c);
        auto *t_67 = buffer.data(target + 67 * ncomps + c);
        auto *t_68 = buffer.data(target + 68 * ncomps + c);
        auto *t_69 = buffer.data(target + 69 * ncomps + c);
        auto *t_70 = buffer.data(target + 70 * ncomps + c);
        auto *t_71 = buffer.data(target + 71 * ncomps + c);
        auto *t_72 = buffer.data(target + 72 * ncomps + c);
        auto *t_73 = buffer.data(target + 73 * ncomps + c);
        auto *t_74 = buffer.data(target + 74 * ncomps + c);
        auto *t_75 = buffer.data(target + 75 * ncomps + c);
        auto *t_76 = buffer.data(target + 76 * ncomps + c);
        auto *t_77 = buffer.data(target + 77 * ncomps + c);
        auto *t_78 = buffer.data(target + 78 * ncomps + c);
        auto *t_79 = buffer.data(target + 79 * ncomps + c);
        auto *t_80 = buffer.data(target + 80 * ncomps + c);
        auto *t_81 = buffer.data(target + 81 * ncomps + c);
        auto *t_82 = buffer.data(target + 82 * ncomps + c);
        auto *t_83 = buffer.data(target + 83 * ncomps + c);
        auto *t_84 = buffer.data(target + 84 * ncomps + c);
        auto *t_85 = buffer.data(target + 85 * ncomps + c);
        auto *t_86 = buffer.data(target + 86 * ncomps + c);
        auto *t_87 = buffer.data(target + 87 * ncomps + c);
        auto *t_88 = buffer.data(target + 88 * ncomps + c);
        auto *t_89 = buffer.data(target + 89 * ncomps + c);
        auto *t_90 = buffer.data(target + 90 * ncomps + c);
        auto *t_91 = buffer.data(target + 91 * ncomps + c);
        auto *t_92 = buffer.data(target + 92 * ncomps + c);
        auto *t_93 = buffer.data(target + 93 * ncomps + c);
        auto *t_94 = buffer.data(target + 94 * ncomps + c);
        auto *t_95 = buffer.data(target + 95 * ncomps + c);
        auto *t_96 = buffer.data(target + 96 * ncomps + c);
        auto *t_97 = buffer.data(target + 97 * ncomps + c);
        auto *t_98 = buffer.data(target + 98 * ncomps + c);
        auto *t_99 = buffer.data(target + 99 * ncomps + c);
        auto *t_100 = buffer.data(target + 100 * ncomps + c);
        auto *t_101 = buffer.data(target + 101 * ncomps + c);
        auto *t_102 = buffer.data(target + 102 * ncomps + c);
        auto *t_103 = buffer.data(target + 103 * ncomps + c);
        auto *t_104 = buffer.data(target + 104 * ncomps + c);
        auto *t_105 = buffer.data(target + 105 * ncomps + c);
        auto *t_106 = buffer.data(target + 106 * ncomps + c);
        auto *t_107 = buffer.data(target + 107 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ks_1_0 = buffer.data(ks_1 + 0 * ncomps + c);
        const auto *ks_1_1 = buffer.data(ks_1 + 1 * ncomps + c);
        const auto *ks_1_2 = buffer.data(ks_1 + 2 * ncomps + c);
        const auto *ks_1_3 = buffer.data(ks_1 + 3 * ncomps + c);
        const auto *ks_1_4 = buffer.data(ks_1 + 4 * ncomps + c);
        const auto *ks_1_5 = buffer.data(ks_1 + 5 * ncomps + c);
        const auto *ks_1_6 = buffer.data(ks_1 + 6 * ncomps + c);
        const auto *ks_1_7 = buffer.data(ks_1 + 7 * ncomps + c);
        const auto *ks_1_8 = buffer.data(ks_1 + 8 * ncomps + c);
        const auto *ks_1_9 = buffer.data(ks_1 + 9 * ncomps + c);
        const auto *ks_1_10 = buffer.data(ks_1 + 10 * ncomps + c);
        const auto *ks_1_11 = buffer.data(ks_1 + 11 * ncomps + c);
        const auto *ks_1_12 = buffer.data(ks_1 + 12 * ncomps + c);
        const auto *ks_1_13 = buffer.data(ks_1 + 13 * ncomps + c);
        const auto *ks_1_14 = buffer.data(ks_1 + 14 * ncomps + c);
        const auto *ks_1_15 = buffer.data(ks_1 + 15 * ncomps + c);
        const auto *ks_1_16 = buffer.data(ks_1 + 16 * ncomps + c);
        const auto *ks_1_17 = buffer.data(ks_1 + 17 * ncomps + c);
        const auto *ks_1_18 = buffer.data(ks_1 + 18 * ncomps + c);
        const auto *ks_1_19 = buffer.data(ks_1 + 19 * ncomps + c);
        const auto *ks_1_20 = buffer.data(ks_1 + 20 * ncomps + c);
        const auto *ks_1_21 = buffer.data(ks_1 + 21 * ncomps + c);
        const auto *ks_1_22 = buffer.data(ks_1 + 22 * ncomps + c);
        const auto *ks_1_23 = buffer.data(ks_1 + 23 * ncomps + c);
        const auto *ks_1_24 = buffer.data(ks_1 + 24 * ncomps + c);
        const auto *ks_1_25 = buffer.data(ks_1 + 25 * ncomps + c);
        const auto *ks_1_26 = buffer.data(ks_1 + 26 * ncomps + c);
        const auto *ks_1_27 = buffer.data(ks_1 + 27 * ncomps + c);
        const auto *ks_1_28 = buffer.data(ks_1 + 28 * ncomps + c);
        const auto *ks_1_29 = buffer.data(ks_1 + 29 * ncomps + c);
        const auto *ks_1_30 = buffer.data(ks_1 + 30 * ncomps + c);
        const auto *ks_1_31 = buffer.data(ks_1 + 31 * ncomps + c);
        const auto *ks_1_32 = buffer.data(ks_1 + 32 * ncomps + c);
        const auto *ks_1_33 = buffer.data(ks_1 + 33 * ncomps + c);
        const auto *ks_1_34 = buffer.data(ks_1 + 34 * ncomps + c);
        const auto *ks_1_35 = buffer.data(ks_1 + 35 * ncomps + c);

        const auto *ks_0_0 = buffer.data(ks_0 + 0 * ncomps + c);
        const auto *ks_0_1 = buffer.data(ks_0 + 1 * ncomps + c);
        const auto *ks_0_2 = buffer.data(ks_0 + 2 * ncomps + c);
        const auto *ks_0_3 = buffer.data(ks_0 + 3 * ncomps + c);
        const auto *ks_0_4 = buffer.data(ks_0 + 4 * ncomps + c);
        const auto *ks_0_5 = buffer.data(ks_0 + 5 * ncomps + c);
        const auto *ks_0_6 = buffer.data(ks_0 + 6 * ncomps + c);
        const auto *ks_0_7 = buffer.data(ks_0 + 7 * ncomps + c);
        const auto *ks_0_8 = buffer.data(ks_0 + 8 * ncomps + c);
        const auto *ks_0_9 = buffer.data(ks_0 + 9 * ncomps + c);
        const auto *ks_0_10 = buffer.data(ks_0 + 10 * ncomps + c);
        const auto *ks_0_11 = buffer.data(ks_0 + 11 * ncomps + c);
        const auto *ks_0_12 = buffer.data(ks_0 + 12 * ncomps + c);
        const auto *ks_0_13 = buffer.data(ks_0 + 13 * ncomps + c);
        const auto *ks_0_14 = buffer.data(ks_0 + 14 * ncomps + c);
        const auto *ks_0_15 = buffer.data(ks_0 + 15 * ncomps + c);
        const auto *ks_0_16 = buffer.data(ks_0 + 16 * ncomps + c);
        const auto *ks_0_17 = buffer.data(ks_0 + 17 * ncomps + c);
        const auto *ks_0_18 = buffer.data(ks_0 + 18 * ncomps + c);
        const auto *ks_0_19 = buffer.data(ks_0 + 19 * ncomps + c);
        const auto *ks_0_20 = buffer.data(ks_0 + 20 * ncomps + c);
        const auto *ks_0_21 = buffer.data(ks_0 + 21 * ncomps + c);
        const auto *ks_0_22 = buffer.data(ks_0 + 22 * ncomps + c);
        const auto *ks_0_23 = buffer.data(ks_0 + 23 * ncomps + c);
        const auto *ks_0_24 = buffer.data(ks_0 + 24 * ncomps + c);
        const auto *ks_0_25 = buffer.data(ks_0 + 25 * ncomps + c);
        const auto *ks_0_26 = buffer.data(ks_0 + 26 * ncomps + c);
        const auto *ks_0_27 = buffer.data(ks_0 + 27 * ncomps + c);
        const auto *ks_0_28 = buffer.data(ks_0 + 28 * ncomps + c);
        const auto *ks_0_29 = buffer.data(ks_0 + 29 * ncomps + c);
        const auto *ks_0_30 = buffer.data(ks_0 + 30 * ncomps + c);
        const auto *ks_0_31 = buffer.data(ks_0 + 31 * ncomps + c);
        const auto *ks_0_32 = buffer.data(ks_0 + 32 * ncomps + c);
        const auto *ks_0_33 = buffer.data(ks_0 + 33 * ncomps + c);
        const auto *ks_0_34 = buffer.data(ks_0 + 34 * ncomps + c);
        const auto *ks_0_35 = buffer.data(ks_0 + 35 * ncomps + c);

        const auto *ls_1_0 = buffer.data(ls_1 + 0 * ncomps + c);
        const auto *ls_1_1 = buffer.data(ls_1 + 1 * ncomps + c);
        const auto *ls_1_2 = buffer.data(ls_1 + 2 * ncomps + c);
        const auto *ls_1_3 = buffer.data(ls_1 + 3 * ncomps + c);
        const auto *ls_1_4 = buffer.data(ls_1 + 4 * ncomps + c);
        const auto *ls_1_5 = buffer.data(ls_1 + 5 * ncomps + c);
        const auto *ls_1_6 = buffer.data(ls_1 + 6 * ncomps + c);
        const auto *ls_1_7 = buffer.data(ls_1 + 7 * ncomps + c);
        const auto *ls_1_8 = buffer.data(ls_1 + 8 * ncomps + c);
        const auto *ls_1_9 = buffer.data(ls_1 + 9 * ncomps + c);
        const auto *ls_1_10 = buffer.data(ls_1 + 10 * ncomps + c);
        const auto *ls_1_11 = buffer.data(ls_1 + 11 * ncomps + c);
        const auto *ls_1_12 = buffer.data(ls_1 + 12 * ncomps + c);
        const auto *ls_1_13 = buffer.data(ls_1 + 13 * ncomps + c);
        const auto *ls_1_14 = buffer.data(ls_1 + 14 * ncomps + c);
        const auto *ls_1_15 = buffer.data(ls_1 + 15 * ncomps + c);
        const auto *ls_1_16 = buffer.data(ls_1 + 16 * ncomps + c);
        const auto *ls_1_17 = buffer.data(ls_1 + 17 * ncomps + c);
        const auto *ls_1_18 = buffer.data(ls_1 + 18 * ncomps + c);
        const auto *ls_1_19 = buffer.data(ls_1 + 19 * ncomps + c);
        const auto *ls_1_20 = buffer.data(ls_1 + 20 * ncomps + c);
        const auto *ls_1_21 = buffer.data(ls_1 + 21 * ncomps + c);
        const auto *ls_1_22 = buffer.data(ls_1 + 22 * ncomps + c);
        const auto *ls_1_23 = buffer.data(ls_1 + 23 * ncomps + c);
        const auto *ls_1_24 = buffer.data(ls_1 + 24 * ncomps + c);
        const auto *ls_1_25 = buffer.data(ls_1 + 25 * ncomps + c);
        const auto *ls_1_26 = buffer.data(ls_1 + 26 * ncomps + c);
        const auto *ls_1_27 = buffer.data(ls_1 + 27 * ncomps + c);
        const auto *ls_1_28 = buffer.data(ls_1 + 28 * ncomps + c);
        const auto *ls_1_29 = buffer.data(ls_1 + 29 * ncomps + c);
        const auto *ls_1_30 = buffer.data(ls_1 + 30 * ncomps + c);
        const auto *ls_1_31 = buffer.data(ls_1 + 31 * ncomps + c);
        const auto *ls_1_32 = buffer.data(ls_1 + 32 * ncomps + c);
        const auto *ls_1_33 = buffer.data(ls_1 + 33 * ncomps + c);
        const auto *ls_1_34 = buffer.data(ls_1 + 34 * ncomps + c);
        const auto *ls_1_35 = buffer.data(ls_1 + 35 * ncomps + c);
        const auto *ls_1_36 = buffer.data(ls_1 + 36 * ncomps + c);
        const auto *ls_1_37 = buffer.data(ls_1 + 37 * ncomps + c);
        const auto *ls_1_38 = buffer.data(ls_1 + 38 * ncomps + c);
        const auto *ls_1_39 = buffer.data(ls_1 + 39 * ncomps + c);
        const auto *ls_1_40 = buffer.data(ls_1 + 40 * ncomps + c);
        const auto *ls_1_41 = buffer.data(ls_1 + 41 * ncomps + c);
        const auto *ls_1_42 = buffer.data(ls_1 + 42 * ncomps + c);
        const auto *ls_1_43 = buffer.data(ls_1 + 43 * ncomps + c);
        const auto *ls_1_44 = buffer.data(ls_1 + 44 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, ab_z, ks_1_0, ks_1_1, ks_0_0, \
                         ls_1_0, ls_1_1, ls_1_2, ls_1_3 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * ks_1_0[k]
                     + ls_1_0[k];

            t_1[k] = ab_y[k] * ks_1_0[k]
                     + ls_1_1[k];

            t_2[k] = ab_z[k] * ks_1_0[k]
                     + ks_0_0[k]
                     + ls_1_2[k];

            t_3[k] = ab_x[k] * ks_1_1[k]
                     + ls_1_1[k];

            t_4[k] = ab_y[k] * ks_1_1[k]
                     + ls_1_3[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_y, ab_z, ks_1_1, ks_1_2, ks_0_1, ks_0_2, \
                         ls_1_2, ls_1_4, ls_1_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_z[k] * ks_1_1[k]
                     + ks_0_1[k]
                     + ls_1_4[k];

            t_6[k] = ab_x[k] * ks_1_2[k]
                     + ls_1_2[k];

            t_7[k] = ab_y[k] * ks_1_2[k]
                     + ls_1_4[k];

            t_8[k] = ab_z[k] * ks_1_2[k]
                     + ks_0_2[k]
                     + ls_1_5[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, ab_x, ab_y, ab_z, ks_1_3, ks_1_4, \
                         ks_0_3, ls_1_3, ls_1_4, ls_1_6, ls_1_7 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = ab_x[k] * ks_1_3[k]
                     + ls_1_3[k];

            t_10[k] = ab_y[k] * ks_1_3[k]
                      + ls_1_6[k];

            t_11[k] = ab_z[k] * ks_1_3[k]
                      + ks_0_3[k]
                      + ls_1_7[k];

            t_12[k] = ab_x[k] * ks_1_4[k]
                      + ls_1_4[k];

            t_13[k] = ab_y[k] * ks_1_4[k]
                      + ls_1_7[k];
        }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, ks_1_4, ks_1_5, ks_0_4, \
                         ks_0_5, ls_1_5, ls_1_8, ls_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_14[k] = ab_z[k] * ks_1_4[k]
                      + ks_0_4[k]
                      + ls_1_8[k];

            t_15[k] = ab_x[k] * ks_1_5[k]
                      + ls_1_5[k];

            t_16[k] = ab_y[k] * ks_1_5[k]
                      + ls_1_8[k];

            t_17[k] = ab_z[k] * ks_1_5[k]
                      + ks_0_5[k]
                      + ls_1_9[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, ab_z, ks_1_6, ks_1_7, \
                         ks_0_6, ls_1_6, ls_1_7, ls_1_10, ls_1_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = ab_x[k] * ks_1_6[k]
                      + ls_1_6[k];

            t_19[k] = ab_y[k] * ks_1_6[k]
                      + ls_1_10[k];

            t_20[k] = ab_z[k] * ks_1_6[k]
                      + ks_0_6[k]
                      + ls_1_11[k];

            t_21[k] = ab_x[k] * ks_1_7[k]
                      + ls_1_7[k];

            t_22[k] = ab_y[k] * ks_1_7[k]
                      + ls_1_11[k];
        }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, ks_1_7, ks_1_8, ks_0_7, \
                         ks_0_8, ls_1_8, ls_1_12, ls_1_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_23[k] = ab_z[k] * ks_1_7[k]
                      + ks_0_7[k]
                      + ls_1_12[k];

            t_24[k] = ab_x[k] * ks_1_8[k]
                      + ls_1_8[k];

            t_25[k] = ab_y[k] * ks_1_8[k]
                      + ls_1_12[k];

            t_26[k] = ab_z[k] * ks_1_8[k]
                      + ks_0_8[k]
                      + ls_1_13[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, ks_1_9, ks_1_10, \
                         ks_0_9, ls_1_9, ls_1_10, ls_1_13, ls_1_14, \
                         ls_1_15 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_x[k] * ks_1_9[k]
                      + ls_1_9[k];

            t_28[k] = ab_y[k] * ks_1_9[k]
                      + ls_1_13[k];

            t_29[k] = ab_z[k] * ks_1_9[k]
                      + ks_0_9[k]
                      + ls_1_14[k];

            t_30[k] = ab_x[k] * ks_1_10[k]
                      + ls_1_10[k];

            t_31[k] = ab_y[k] * ks_1_10[k]
                      + ls_1_15[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, ks_1_10, ks_1_11, ks_0_10, \
                         ks_0_11, ls_1_11, ls_1_16, ls_1_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_z[k] * ks_1_10[k]
                      + ks_0_10[k]
                      + ls_1_16[k];

            t_33[k] = ab_x[k] * ks_1_11[k]
                      + ls_1_11[k];

            t_34[k] = ab_y[k] * ks_1_11[k]
                      + ls_1_16[k];

            t_35[k] = ab_z[k] * ks_1_11[k]
                      + ks_0_11[k]
                      + ls_1_17[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_x, ab_y, ab_z, ks_1_12, ks_1_13, \
                         ks_0_12, ls_1_12, ls_1_13, ls_1_17, ls_1_18 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = ab_x[k] * ks_1_12[k]
                      + ls_1_12[k];

            t_37[k] = ab_y[k] * ks_1_12[k]
                      + ls_1_17[k];

            t_38[k] = ab_z[k] * ks_1_12[k]
                      + ks_0_12[k]
                      + ls_1_18[k];

            t_39[k] = ab_x[k] * ks_1_13[k]
                      + ls_1_13[k];

            t_40[k] = ab_y[k] * ks_1_13[k]
                      + ls_1_18[k];
        }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_x, ab_y, ab_z, ks_1_13, ks_1_14, ks_0_13, \
                         ks_0_14, ls_1_14, ls_1_19, ls_1_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_41[k] = ab_z[k] * ks_1_13[k]
                      + ks_0_13[k]
                      + ls_1_19[k];

            t_42[k] = ab_x[k] * ks_1_14[k]
                      + ls_1_14[k];

            t_43[k] = ab_y[k] * ks_1_14[k]
                      + ls_1_19[k];

            t_44[k] = ab_z[k] * ks_1_14[k]
                      + ks_0_14[k]
                      + ls_1_20[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, ab_y, ab_z, ks_1_15, ks_1_16, \
                         ks_0_15, ls_1_15, ls_1_16, ls_1_21, ls_1_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * ks_1_15[k]
                      + ls_1_15[k];

            t_46[k] = ab_y[k] * ks_1_15[k]
                      + ls_1_21[k];

            t_47[k] = ab_z[k] * ks_1_15[k]
                      + ks_0_15[k]
                      + ls_1_22[k];

            t_48[k] = ab_x[k] * ks_1_16[k]
                      + ls_1_16[k];

            t_49[k] = ab_y[k] * ks_1_16[k]
                      + ls_1_22[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, ks_1_16, ks_1_17, ks_0_16, \
                         ks_0_17, ls_1_17, ls_1_23, ls_1_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = ab_z[k] * ks_1_16[k]
                      + ks_0_16[k]
                      + ls_1_23[k];

            t_51[k] = ab_x[k] * ks_1_17[k]
                      + ls_1_17[k];

            t_52[k] = ab_y[k] * ks_1_17[k]
                      + ls_1_23[k];

            t_53[k] = ab_z[k] * ks_1_17[k]
                      + ks_0_17[k]
                      + ls_1_24[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_x, ab_y, ab_z, ks_1_18, ks_1_19, \
                         ks_0_18, ls_1_18, ls_1_19, ls_1_24, ls_1_25 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * ks_1_18[k]
                      + ls_1_18[k];

            t_55[k] = ab_y[k] * ks_1_18[k]
                      + ls_1_24[k];

            t_56[k] = ab_z[k] * ks_1_18[k]
                      + ks_0_18[k]
                      + ls_1_25[k];

            t_57[k] = ab_x[k] * ks_1_19[k]
                      + ls_1_19[k];

            t_58[k] = ab_y[k] * ks_1_19[k]
                      + ls_1_25[k];
        }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, ab_x, ab_y, ab_z, ks_1_19, ks_1_20, ks_0_19, \
                         ks_0_20, ls_1_20, ls_1_26, ls_1_27 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_59[k] = ab_z[k] * ks_1_19[k]
                      + ks_0_19[k]
                      + ls_1_26[k];

            t_60[k] = ab_x[k] * ks_1_20[k]
                      + ls_1_20[k];

            t_61[k] = ab_y[k] * ks_1_20[k]
                      + ls_1_26[k];

            t_62[k] = ab_z[k] * ks_1_20[k]
                      + ks_0_20[k]
                      + ls_1_27[k];
        }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, ab_x, ab_y, ab_z, ks_1_21, ks_1_22, \
                         ks_0_21, ls_1_21, ls_1_22, ls_1_28, ls_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_63[k] = ab_x[k] * ks_1_21[k]
                      + ls_1_21[k];

            t_64[k] = ab_y[k] * ks_1_21[k]
                      + ls_1_28[k];

            t_65[k] = ab_z[k] * ks_1_21[k]
                      + ks_0_21[k]
                      + ls_1_29[k];

            t_66[k] = ab_x[k] * ks_1_22[k]
                      + ls_1_22[k];

            t_67[k] = ab_y[k] * ks_1_22[k]
                      + ls_1_29[k];
        }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, ab_x, ab_y, ab_z, ks_1_22, ks_1_23, ks_0_22, \
                         ks_0_23, ls_1_23, ls_1_30, ls_1_31 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_68[k] = ab_z[k] * ks_1_22[k]
                      + ks_0_22[k]
                      + ls_1_30[k];

            t_69[k] = ab_x[k] * ks_1_23[k]
                      + ls_1_23[k];

            t_70[k] = ab_y[k] * ks_1_23[k]
                      + ls_1_30[k];

            t_71[k] = ab_z[k] * ks_1_23[k]
                      + ks_0_23[k]
                      + ls_1_31[k];
        }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, ab_z, ks_1_24, ks_1_25, \
                         ks_0_24, ls_1_24, ls_1_25, ls_1_31, ls_1_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_72[k] = ab_x[k] * ks_1_24[k]
                      + ls_1_24[k];

            t_73[k] = ab_y[k] * ks_1_24[k]
                      + ls_1_31[k];

            t_74[k] = ab_z[k] * ks_1_24[k]
                      + ks_0_24[k]
                      + ls_1_32[k];

            t_75[k] = ab_x[k] * ks_1_25[k]
                      + ls_1_25[k];

            t_76[k] = ab_y[k] * ks_1_25[k]
                      + ls_1_32[k];
        }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, ab_x, ab_y, ab_z, ks_1_25, ks_1_26, ks_0_25, \
                         ks_0_26, ls_1_26, ls_1_33, ls_1_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_77[k] = ab_z[k] * ks_1_25[k]
                      + ks_0_25[k]
                      + ls_1_33[k];

            t_78[k] = ab_x[k] * ks_1_26[k]
                      + ls_1_26[k];

            t_79[k] = ab_y[k] * ks_1_26[k]
                      + ls_1_33[k];

            t_80[k] = ab_z[k] * ks_1_26[k]
                      + ks_0_26[k]
                      + ls_1_34[k];
        }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, ab_x, ab_y, ab_z, ks_1_27, ks_1_28, \
                         ks_0_27, ls_1_27, ls_1_28, ls_1_34, ls_1_35, \
                         ls_1_36 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_81[k] = ab_x[k] * ks_1_27[k]
                      + ls_1_27[k];

            t_82[k] = ab_y[k] * ks_1_27[k]
                      + ls_1_34[k];

            t_83[k] = ab_z[k] * ks_1_27[k]
                      + ks_0_27[k]
                      + ls_1_35[k];

            t_84[k] = ab_x[k] * ks_1_28[k]
                      + ls_1_28[k];

            t_85[k] = ab_y[k] * ks_1_28[k]
                      + ls_1_36[k];
        }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, ks_1_28, ks_1_29, ks_0_28, \
                         ks_0_29, ls_1_29, ls_1_37, ls_1_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_86[k] = ab_z[k] * ks_1_28[k]
                      + ks_0_28[k]
                      + ls_1_37[k];

            t_87[k] = ab_x[k] * ks_1_29[k]
                      + ls_1_29[k];

            t_88[k] = ab_y[k] * ks_1_29[k]
                      + ls_1_37[k];

            t_89[k] = ab_z[k] * ks_1_29[k]
                      + ks_0_29[k]
                      + ls_1_38[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, ab_y, ab_z, ks_1_30, ks_1_31, \
                         ks_0_30, ls_1_30, ls_1_31, ls_1_38, ls_1_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * ks_1_30[k]
                      + ls_1_30[k];

            t_91[k] = ab_y[k] * ks_1_30[k]
                      + ls_1_38[k];

            t_92[k] = ab_z[k] * ks_1_30[k]
                      + ks_0_30[k]
                      + ls_1_39[k];

            t_93[k] = ab_x[k] * ks_1_31[k]
                      + ls_1_31[k];

            t_94[k] = ab_y[k] * ks_1_31[k]
                      + ls_1_39[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, ab_x, ab_y, ab_z, ks_1_31, ks_1_32, ks_0_31, \
                         ks_0_32, ls_1_32, ls_1_40, ls_1_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_z[k] * ks_1_31[k]
                      + ks_0_31[k]
                      + ls_1_40[k];

            t_96[k] = ab_x[k] * ks_1_32[k]
                      + ls_1_32[k];

            t_97[k] = ab_y[k] * ks_1_32[k]
                      + ls_1_40[k];

            t_98[k] = ab_z[k] * ks_1_32[k]
                      + ks_0_32[k]
                      + ls_1_41[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, ab_x, ab_y, ab_z, ks_1_33, ks_1_34, \
                         ks_0_33, ls_1_33, ls_1_34, ls_1_41, ls_1_42 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_x[k] * ks_1_33[k]
                      + ls_1_33[k];

            t_100[k] = ab_y[k] * ks_1_33[k]
                       + ls_1_41[k];

            t_101[k] = ab_z[k] * ks_1_33[k]
                       + ks_0_33[k]
                       + ls_1_42[k];

            t_102[k] = ab_x[k] * ks_1_34[k]
                       + ls_1_34[k];

            t_103[k] = ab_y[k] * ks_1_34[k]
                       + ls_1_42[k];
        }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, ab_x, ab_y, ab_z, ks_1_34, ks_1_35, \
                         ks_0_34, ks_0_35, ls_1_35, ls_1_43, ls_1_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_104[k] = ab_z[k] * ks_1_34[k]
                       + ks_0_34[k]
                       + ls_1_43[k];

            t_105[k] = ab_x[k] * ks_1_35[k]
                       + ls_1_35[k];

            t_106[k] = ab_y[k] * ks_1_35[k]
                       + ls_1_43[k];

            t_107[k] = ab_z[k] * ks_1_35[k]
                       + ks_0_35[k]
                       + ls_1_44[k];
        }
    }
}

}  // namespace simdtrf
