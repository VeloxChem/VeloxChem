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


#include "SimdTransferGeom010ZFF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_010z_ff(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                         const size_t target, const size_t df_1, const size_t df_0,
                         const size_t dg_1, const size_t ncomps, const size_t nmax) -> void
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *df_1_0 = buffer.data(df_1 + 0 * ncomps + c);
        const auto *df_1_1 = buffer.data(df_1 + 1 * ncomps + c);
        const auto *df_1_2 = buffer.data(df_1 + 2 * ncomps + c);
        const auto *df_1_3 = buffer.data(df_1 + 3 * ncomps + c);
        const auto *df_1_4 = buffer.data(df_1 + 4 * ncomps + c);
        const auto *df_1_5 = buffer.data(df_1 + 5 * ncomps + c);
        const auto *df_1_6 = buffer.data(df_1 + 6 * ncomps + c);
        const auto *df_1_7 = buffer.data(df_1 + 7 * ncomps + c);
        const auto *df_1_8 = buffer.data(df_1 + 8 * ncomps + c);
        const auto *df_1_9 = buffer.data(df_1 + 9 * ncomps + c);
        const auto *df_1_10 = buffer.data(df_1 + 10 * ncomps + c);
        const auto *df_1_11 = buffer.data(df_1 + 11 * ncomps + c);
        const auto *df_1_12 = buffer.data(df_1 + 12 * ncomps + c);
        const auto *df_1_13 = buffer.data(df_1 + 13 * ncomps + c);
        const auto *df_1_14 = buffer.data(df_1 + 14 * ncomps + c);
        const auto *df_1_15 = buffer.data(df_1 + 15 * ncomps + c);
        const auto *df_1_16 = buffer.data(df_1 + 16 * ncomps + c);
        const auto *df_1_17 = buffer.data(df_1 + 17 * ncomps + c);
        const auto *df_1_18 = buffer.data(df_1 + 18 * ncomps + c);
        const auto *df_1_19 = buffer.data(df_1 + 19 * ncomps + c);
        const auto *df_1_20 = buffer.data(df_1 + 20 * ncomps + c);
        const auto *df_1_21 = buffer.data(df_1 + 21 * ncomps + c);
        const auto *df_1_22 = buffer.data(df_1 + 22 * ncomps + c);
        const auto *df_1_23 = buffer.data(df_1 + 23 * ncomps + c);
        const auto *df_1_24 = buffer.data(df_1 + 24 * ncomps + c);
        const auto *df_1_25 = buffer.data(df_1 + 25 * ncomps + c);
        const auto *df_1_26 = buffer.data(df_1 + 26 * ncomps + c);
        const auto *df_1_27 = buffer.data(df_1 + 27 * ncomps + c);
        const auto *df_1_28 = buffer.data(df_1 + 28 * ncomps + c);
        const auto *df_1_29 = buffer.data(df_1 + 29 * ncomps + c);
        const auto *df_1_30 = buffer.data(df_1 + 30 * ncomps + c);
        const auto *df_1_31 = buffer.data(df_1 + 31 * ncomps + c);
        const auto *df_1_32 = buffer.data(df_1 + 32 * ncomps + c);
        const auto *df_1_33 = buffer.data(df_1 + 33 * ncomps + c);
        const auto *df_1_34 = buffer.data(df_1 + 34 * ncomps + c);
        const auto *df_1_35 = buffer.data(df_1 + 35 * ncomps + c);
        const auto *df_1_36 = buffer.data(df_1 + 36 * ncomps + c);
        const auto *df_1_37 = buffer.data(df_1 + 37 * ncomps + c);
        const auto *df_1_38 = buffer.data(df_1 + 38 * ncomps + c);
        const auto *df_1_39 = buffer.data(df_1 + 39 * ncomps + c);
        const auto *df_1_40 = buffer.data(df_1 + 40 * ncomps + c);
        const auto *df_1_41 = buffer.data(df_1 + 41 * ncomps + c);
        const auto *df_1_42 = buffer.data(df_1 + 42 * ncomps + c);
        const auto *df_1_43 = buffer.data(df_1 + 43 * ncomps + c);
        const auto *df_1_44 = buffer.data(df_1 + 44 * ncomps + c);
        const auto *df_1_45 = buffer.data(df_1 + 45 * ncomps + c);
        const auto *df_1_46 = buffer.data(df_1 + 46 * ncomps + c);
        const auto *df_1_47 = buffer.data(df_1 + 47 * ncomps + c);
        const auto *df_1_48 = buffer.data(df_1 + 48 * ncomps + c);
        const auto *df_1_49 = buffer.data(df_1 + 49 * ncomps + c);
        const auto *df_1_50 = buffer.data(df_1 + 50 * ncomps + c);
        const auto *df_1_51 = buffer.data(df_1 + 51 * ncomps + c);
        const auto *df_1_52 = buffer.data(df_1 + 52 * ncomps + c);
        const auto *df_1_53 = buffer.data(df_1 + 53 * ncomps + c);
        const auto *df_1_54 = buffer.data(df_1 + 54 * ncomps + c);
        const auto *df_1_55 = buffer.data(df_1 + 55 * ncomps + c);
        const auto *df_1_56 = buffer.data(df_1 + 56 * ncomps + c);
        const auto *df_1_57 = buffer.data(df_1 + 57 * ncomps + c);
        const auto *df_1_58 = buffer.data(df_1 + 58 * ncomps + c);
        const auto *df_1_59 = buffer.data(df_1 + 59 * ncomps + c);

        const auto *df_0_50 = buffer.data(df_0 + 50 * ncomps + c);
        const auto *df_0_51 = buffer.data(df_0 + 51 * ncomps + c);
        const auto *df_0_52 = buffer.data(df_0 + 52 * ncomps + c);
        const auto *df_0_53 = buffer.data(df_0 + 53 * ncomps + c);
        const auto *df_0_54 = buffer.data(df_0 + 54 * ncomps + c);
        const auto *df_0_55 = buffer.data(df_0 + 55 * ncomps + c);
        const auto *df_0_56 = buffer.data(df_0 + 56 * ncomps + c);
        const auto *df_0_57 = buffer.data(df_0 + 57 * ncomps + c);
        const auto *df_0_58 = buffer.data(df_0 + 58 * ncomps + c);
        const auto *df_0_59 = buffer.data(df_0 + 59 * ncomps + c);

        const auto *dg_1_0 = buffer.data(dg_1 + 0 * ncomps + c);
        const auto *dg_1_1 = buffer.data(dg_1 + 1 * ncomps + c);
        const auto *dg_1_2 = buffer.data(dg_1 + 2 * ncomps + c);
        const auto *dg_1_3 = buffer.data(dg_1 + 3 * ncomps + c);
        const auto *dg_1_4 = buffer.data(dg_1 + 4 * ncomps + c);
        const auto *dg_1_5 = buffer.data(dg_1 + 5 * ncomps + c);
        const auto *dg_1_6 = buffer.data(dg_1 + 6 * ncomps + c);
        const auto *dg_1_7 = buffer.data(dg_1 + 7 * ncomps + c);
        const auto *dg_1_8 = buffer.data(dg_1 + 8 * ncomps + c);
        const auto *dg_1_9 = buffer.data(dg_1 + 9 * ncomps + c);
        const auto *dg_1_15 = buffer.data(dg_1 + 15 * ncomps + c);
        const auto *dg_1_16 = buffer.data(dg_1 + 16 * ncomps + c);
        const auto *dg_1_17 = buffer.data(dg_1 + 17 * ncomps + c);
        const auto *dg_1_18 = buffer.data(dg_1 + 18 * ncomps + c);
        const auto *dg_1_19 = buffer.data(dg_1 + 19 * ncomps + c);
        const auto *dg_1_20 = buffer.data(dg_1 + 20 * ncomps + c);
        const auto *dg_1_21 = buffer.data(dg_1 + 21 * ncomps + c);
        const auto *dg_1_22 = buffer.data(dg_1 + 22 * ncomps + c);
        const auto *dg_1_23 = buffer.data(dg_1 + 23 * ncomps + c);
        const auto *dg_1_24 = buffer.data(dg_1 + 24 * ncomps + c);
        const auto *dg_1_30 = buffer.data(dg_1 + 30 * ncomps + c);
        const auto *dg_1_31 = buffer.data(dg_1 + 31 * ncomps + c);
        const auto *dg_1_32 = buffer.data(dg_1 + 32 * ncomps + c);
        const auto *dg_1_33 = buffer.data(dg_1 + 33 * ncomps + c);
        const auto *dg_1_34 = buffer.data(dg_1 + 34 * ncomps + c);
        const auto *dg_1_35 = buffer.data(dg_1 + 35 * ncomps + c);
        const auto *dg_1_36 = buffer.data(dg_1 + 36 * ncomps + c);
        const auto *dg_1_37 = buffer.data(dg_1 + 37 * ncomps + c);
        const auto *dg_1_38 = buffer.data(dg_1 + 38 * ncomps + c);
        const auto *dg_1_39 = buffer.data(dg_1 + 39 * ncomps + c);
        const auto *dg_1_45 = buffer.data(dg_1 + 45 * ncomps + c);
        const auto *dg_1_46 = buffer.data(dg_1 + 46 * ncomps + c);
        const auto *dg_1_47 = buffer.data(dg_1 + 47 * ncomps + c);
        const auto *dg_1_48 = buffer.data(dg_1 + 48 * ncomps + c);
        const auto *dg_1_49 = buffer.data(dg_1 + 49 * ncomps + c);
        const auto *dg_1_50 = buffer.data(dg_1 + 50 * ncomps + c);
        const auto *dg_1_51 = buffer.data(dg_1 + 51 * ncomps + c);
        const auto *dg_1_52 = buffer.data(dg_1 + 52 * ncomps + c);
        const auto *dg_1_53 = buffer.data(dg_1 + 53 * ncomps + c);
        const auto *dg_1_54 = buffer.data(dg_1 + 54 * ncomps + c);
        const auto *dg_1_55 = buffer.data(dg_1 + 55 * ncomps + c);
        const auto *dg_1_56 = buffer.data(dg_1 + 56 * ncomps + c);
        const auto *dg_1_57 = buffer.data(dg_1 + 57 * ncomps + c);
        const auto *dg_1_58 = buffer.data(dg_1 + 58 * ncomps + c);
        const auto *dg_1_60 = buffer.data(dg_1 + 60 * ncomps + c);
        const auto *dg_1_61 = buffer.data(dg_1 + 61 * ncomps + c);
        const auto *dg_1_62 = buffer.data(dg_1 + 62 * ncomps + c);
        const auto *dg_1_63 = buffer.data(dg_1 + 63 * ncomps + c);
        const auto *dg_1_64 = buffer.data(dg_1 + 64 * ncomps + c);
        const auto *dg_1_65 = buffer.data(dg_1 + 65 * ncomps + c);
        const auto *dg_1_66 = buffer.data(dg_1 + 66 * ncomps + c);
        const auto *dg_1_67 = buffer.data(dg_1 + 67 * ncomps + c);
        const auto *dg_1_68 = buffer.data(dg_1 + 68 * ncomps + c);
        const auto *dg_1_69 = buffer.data(dg_1 + 69 * ncomps + c);
        const auto *dg_1_70 = buffer.data(dg_1 + 70 * ncomps + c);
        const auto *dg_1_71 = buffer.data(dg_1 + 71 * ncomps + c);
        const auto *dg_1_72 = buffer.data(dg_1 + 72 * ncomps + c);
        const auto *dg_1_73 = buffer.data(dg_1 + 73 * ncomps + c);
        const auto *dg_1_75 = buffer.data(dg_1 + 75 * ncomps + c);
        const auto *dg_1_76 = buffer.data(dg_1 + 76 * ncomps + c);
        const auto *dg_1_77 = buffer.data(dg_1 + 77 * ncomps + c);
        const auto *dg_1_78 = buffer.data(dg_1 + 78 * ncomps + c);
        const auto *dg_1_79 = buffer.data(dg_1 + 79 * ncomps + c);
        const auto *dg_1_80 = buffer.data(dg_1 + 80 * ncomps + c);
        const auto *dg_1_81 = buffer.data(dg_1 + 81 * ncomps + c);
        const auto *dg_1_82 = buffer.data(dg_1 + 82 * ncomps + c);
        const auto *dg_1_83 = buffer.data(dg_1 + 83 * ncomps + c);
        const auto *dg_1_84 = buffer.data(dg_1 + 84 * ncomps + c);
        const auto *dg_1_85 = buffer.data(dg_1 + 85 * ncomps + c);
        const auto *dg_1_86 = buffer.data(dg_1 + 86 * ncomps + c);
        const auto *dg_1_87 = buffer.data(dg_1 + 87 * ncomps + c);
        const auto *dg_1_88 = buffer.data(dg_1 + 88 * ncomps + c);
        const auto *dg_1_89 = buffer.data(dg_1 + 89 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, df_1_0, df_1_1, df_1_2, df_1_3, \
                         df_1_4, dg_1_0, dg_1_1, dg_1_2, dg_1_3, \
                         dg_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * df_1_0[k]
                     + dg_1_0[k];

            t_1[k] = -ab_x[k] * df_1_1[k]
                     + dg_1_1[k];

            t_2[k] = -ab_x[k] * df_1_2[k]
                     + dg_1_2[k];

            t_3[k] = -ab_x[k] * df_1_3[k]
                     + dg_1_3[k];

            t_4[k] = -ab_x[k] * df_1_4[k]
                     + dg_1_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, df_1_5, df_1_6, df_1_7, df_1_8, \
                         df_1_9, dg_1_5, dg_1_6, dg_1_7, dg_1_8, \
                         dg_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * df_1_5[k]
                     + dg_1_5[k];

            t_6[k] = -ab_x[k] * df_1_6[k]
                     + dg_1_6[k];

            t_7[k] = -ab_x[k] * df_1_7[k]
                     + dg_1_7[k];

            t_8[k] = -ab_x[k] * df_1_8[k]
                     + dg_1_8[k];

            t_9[k] = -ab_x[k] * df_1_9[k]
                     + dg_1_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, df_1_10, df_1_11, df_1_12, \
                         df_1_13, df_1_14, dg_1_15, dg_1_16, dg_1_17, dg_1_18, \
                         dg_1_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * df_1_10[k]
                      + dg_1_15[k];

            t_11[k] = -ab_x[k] * df_1_11[k]
                      + dg_1_16[k];

            t_12[k] = -ab_x[k] * df_1_12[k]
                      + dg_1_17[k];

            t_13[k] = -ab_x[k] * df_1_13[k]
                      + dg_1_18[k];

            t_14[k] = -ab_x[k] * df_1_14[k]
                      + dg_1_19[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, df_1_15, df_1_16, df_1_17, \
                         df_1_18, df_1_19, dg_1_20, dg_1_21, dg_1_22, dg_1_23, \
                         dg_1_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * df_1_15[k]
                      + dg_1_20[k];

            t_16[k] = -ab_x[k] * df_1_16[k]
                      + dg_1_21[k];

            t_17[k] = -ab_x[k] * df_1_17[k]
                      + dg_1_22[k];

            t_18[k] = -ab_x[k] * df_1_18[k]
                      + dg_1_23[k];

            t_19[k] = -ab_x[k] * df_1_19[k]
                      + dg_1_24[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, df_1_20, df_1_21, df_1_22, \
                         df_1_23, df_1_24, dg_1_30, dg_1_31, dg_1_32, dg_1_33, \
                         dg_1_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * df_1_20[k]
                      + dg_1_30[k];

            t_21[k] = -ab_x[k] * df_1_21[k]
                      + dg_1_31[k];

            t_22[k] = -ab_x[k] * df_1_22[k]
                      + dg_1_32[k];

            t_23[k] = -ab_x[k] * df_1_23[k]
                      + dg_1_33[k];

            t_24[k] = -ab_x[k] * df_1_24[k]
                      + dg_1_34[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, df_1_25, df_1_26, df_1_27, \
                         df_1_28, df_1_29, dg_1_35, dg_1_36, dg_1_37, dg_1_38, \
                         dg_1_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * df_1_25[k]
                      + dg_1_35[k];

            t_26[k] = -ab_x[k] * df_1_26[k]
                      + dg_1_36[k];

            t_27[k] = -ab_x[k] * df_1_27[k]
                      + dg_1_37[k];

            t_28[k] = -ab_x[k] * df_1_28[k]
                      + dg_1_38[k];

            t_29[k] = -ab_x[k] * df_1_29[k]
                      + dg_1_39[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, df_1_30, df_1_31, df_1_32, \
                         df_1_33, df_1_34, dg_1_45, dg_1_46, dg_1_47, dg_1_48, \
                         dg_1_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * df_1_30[k]
                      + dg_1_45[k];

            t_31[k] = -ab_x[k] * df_1_31[k]
                      + dg_1_46[k];

            t_32[k] = -ab_x[k] * df_1_32[k]
                      + dg_1_47[k];

            t_33[k] = -ab_x[k] * df_1_33[k]
                      + dg_1_48[k];

            t_34[k] = -ab_x[k] * df_1_34[k]
                      + dg_1_49[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, df_1_35, df_1_36, df_1_37, \
                         df_1_38, df_1_39, dg_1_50, dg_1_51, dg_1_52, dg_1_53, \
                         dg_1_54 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * df_1_35[k]
                      + dg_1_50[k];

            t_36[k] = -ab_x[k] * df_1_36[k]
                      + dg_1_51[k];

            t_37[k] = -ab_x[k] * df_1_37[k]
                      + dg_1_52[k];

            t_38[k] = -ab_x[k] * df_1_38[k]
                      + dg_1_53[k];

            t_39[k] = -ab_x[k] * df_1_39[k]
                      + dg_1_54[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, df_1_40, df_1_41, df_1_42, \
                         df_1_43, df_1_44, dg_1_60, dg_1_61, dg_1_62, dg_1_63, \
                         dg_1_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * df_1_40[k]
                      + dg_1_60[k];

            t_41[k] = -ab_x[k] * df_1_41[k]
                      + dg_1_61[k];

            t_42[k] = -ab_x[k] * df_1_42[k]
                      + dg_1_62[k];

            t_43[k] = -ab_x[k] * df_1_43[k]
                      + dg_1_63[k];

            t_44[k] = -ab_x[k] * df_1_44[k]
                      + dg_1_64[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, df_1_45, df_1_46, df_1_47, \
                         df_1_48, df_1_49, dg_1_65, dg_1_66, dg_1_67, dg_1_68, \
                         dg_1_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * df_1_45[k]
                      + dg_1_65[k];

            t_46[k] = -ab_x[k] * df_1_46[k]
                      + dg_1_66[k];

            t_47[k] = -ab_x[k] * df_1_47[k]
                      + dg_1_67[k];

            t_48[k] = -ab_x[k] * df_1_48[k]
                      + dg_1_68[k];

            t_49[k] = -ab_x[k] * df_1_49[k]
                      + dg_1_69[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, df_1_50, df_1_51, df_1_52, \
                         df_1_53, df_1_54, dg_1_75, dg_1_76, dg_1_77, dg_1_78, \
                         dg_1_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * df_1_50[k]
                      + dg_1_75[k];

            t_51[k] = -ab_x[k] * df_1_51[k]
                      + dg_1_76[k];

            t_52[k] = -ab_x[k] * df_1_52[k]
                      + dg_1_77[k];

            t_53[k] = -ab_x[k] * df_1_53[k]
                      + dg_1_78[k];

            t_54[k] = -ab_x[k] * df_1_54[k]
                      + dg_1_79[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, df_1_55, df_1_56, df_1_57, \
                         df_1_58, df_1_59, dg_1_80, dg_1_81, dg_1_82, dg_1_83, \
                         dg_1_84 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * df_1_55[k]
                      + dg_1_80[k];

            t_56[k] = -ab_x[k] * df_1_56[k]
                      + dg_1_81[k];

            t_57[k] = -ab_x[k] * df_1_57[k]
                      + dg_1_82[k];

            t_58[k] = -ab_x[k] * df_1_58[k]
                      + dg_1_83[k];

            t_59[k] = -ab_x[k] * df_1_59[k]
                      + dg_1_84[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_y, df_1_30, df_1_31, df_1_32, \
                         df_1_33, df_1_34, dg_1_46, dg_1_48, dg_1_49, dg_1_51, \
                         dg_1_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_y[k] * df_1_30[k]
                      + dg_1_46[k];

            t_61[k] = -ab_y[k] * df_1_31[k]
                      + dg_1_48[k];

            t_62[k] = -ab_y[k] * df_1_32[k]
                      + dg_1_49[k];

            t_63[k] = -ab_y[k] * df_1_33[k]
                      + dg_1_51[k];

            t_64[k] = -ab_y[k] * df_1_34[k]
                      + dg_1_52[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_y, df_1_35, df_1_36, df_1_37, \
                         df_1_38, df_1_39, dg_1_53, dg_1_55, dg_1_56, dg_1_57, \
                         dg_1_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_y[k] * df_1_35[k]
                      + dg_1_53[k];

            t_66[k] = -ab_y[k] * df_1_36[k]
                      + dg_1_55[k];

            t_67[k] = -ab_y[k] * df_1_37[k]
                      + dg_1_56[k];

            t_68[k] = -ab_y[k] * df_1_38[k]
                      + dg_1_57[k];

            t_69[k] = -ab_y[k] * df_1_39[k]
                      + dg_1_58[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_y, df_1_40, df_1_41, df_1_42, \
                         df_1_43, df_1_44, dg_1_61, dg_1_63, dg_1_64, dg_1_66, \
                         dg_1_67 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_y[k] * df_1_40[k]
                      + dg_1_61[k];

            t_71[k] = -ab_y[k] * df_1_41[k]
                      + dg_1_63[k];

            t_72[k] = -ab_y[k] * df_1_42[k]
                      + dg_1_64[k];

            t_73[k] = -ab_y[k] * df_1_43[k]
                      + dg_1_66[k];

            t_74[k] = -ab_y[k] * df_1_44[k]
                      + dg_1_67[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_y, df_1_45, df_1_46, df_1_47, \
                         df_1_48, df_1_49, dg_1_68, dg_1_70, dg_1_71, dg_1_72, \
                         dg_1_73 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_y[k] * df_1_45[k]
                      + dg_1_68[k];

            t_76[k] = -ab_y[k] * df_1_46[k]
                      + dg_1_70[k];

            t_77[k] = -ab_y[k] * df_1_47[k]
                      + dg_1_71[k];

            t_78[k] = -ab_y[k] * df_1_48[k]
                      + dg_1_72[k];

            t_79[k] = -ab_y[k] * df_1_49[k]
                      + dg_1_73[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_y, df_1_50, df_1_51, df_1_52, \
                         df_1_53, df_1_54, dg_1_76, dg_1_78, dg_1_79, dg_1_81, \
                         dg_1_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_y[k] * df_1_50[k]
                      + dg_1_76[k];

            t_81[k] = -ab_y[k] * df_1_51[k]
                      + dg_1_78[k];

            t_82[k] = -ab_y[k] * df_1_52[k]
                      + dg_1_79[k];

            t_83[k] = -ab_y[k] * df_1_53[k]
                      + dg_1_81[k];

            t_84[k] = -ab_y[k] * df_1_54[k]
                      + dg_1_82[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_y, df_1_55, df_1_56, df_1_57, \
                         df_1_58, df_1_59, dg_1_83, dg_1_85, dg_1_86, dg_1_87, \
                         dg_1_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_y[k] * df_1_55[k]
                      + dg_1_83[k];

            t_86[k] = -ab_y[k] * df_1_56[k]
                      + dg_1_85[k];

            t_87[k] = -ab_y[k] * df_1_57[k]
                      + dg_1_86[k];

            t_88[k] = -ab_y[k] * df_1_58[k]
                      + dg_1_87[k];

            t_89[k] = -ab_y[k] * df_1_59[k]
                      + dg_1_88[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, ab_z, df_1_50, df_1_51, df_1_52, df_0_50, df_0_51, \
                         df_0_52, dg_1_77, dg_1_79, dg_1_80 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_z[k] * df_1_50[k]
                      + df_0_50[k]
                      + dg_1_77[k];

            t_91[k] = -ab_z[k] * df_1_51[k]
                      + df_0_51[k]
                      + dg_1_79[k];

            t_92[k] = -ab_z[k] * df_1_52[k]
                      + df_0_52[k]
                      + dg_1_80[k];
        }

#pragma omp simd aligned(t_93, t_94, t_95, ab_z, df_1_53, df_1_54, df_1_55, df_0_53, df_0_54, \
                         df_0_55, dg_1_82, dg_1_83, dg_1_84 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_93[k] = -ab_z[k] * df_1_53[k]
                      + df_0_53[k]
                      + dg_1_82[k];

            t_94[k] = -ab_z[k] * df_1_54[k]
                      + df_0_54[k]
                      + dg_1_83[k];

            t_95[k] = -ab_z[k] * df_1_55[k]
                      + df_0_55[k]
                      + dg_1_84[k];
        }

#pragma omp simd aligned(t_96, t_97, t_98, ab_z, df_1_56, df_1_57, df_1_58, df_0_56, df_0_57, \
                         df_0_58, dg_1_86, dg_1_87, dg_1_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_96[k] = -ab_z[k] * df_1_56[k]
                      + df_0_56[k]
                      + dg_1_86[k];

            t_97[k] = -ab_z[k] * df_1_57[k]
                      + df_0_57[k]
                      + dg_1_87[k];

            t_98[k] = -ab_z[k] * df_1_58[k]
                      + df_0_58[k]
                      + dg_1_88[k];
        }

#pragma omp simd aligned(t_99, ab_z, df_1_59, df_0_59, dg_1_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = -ab_z[k] * df_1_59[k]
                      + df_0_59[k]
                      + dg_1_89[k];
        }
    }
}

}  // namespace simdtrf
