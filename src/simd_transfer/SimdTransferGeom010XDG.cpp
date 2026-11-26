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


#include "SimdTransferGeom010XDG.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_010x_dg(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                         const size_t target, const size_t pg_1, const size_t pg_0,
                         const size_t ph_1, const size_t ncomps, const size_t nmax) -> void
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *pg_1_0 = buffer.data(pg_1 + 0 * ncomps + c);
        const auto *pg_1_1 = buffer.data(pg_1 + 1 * ncomps + c);
        const auto *pg_1_2 = buffer.data(pg_1 + 2 * ncomps + c);
        const auto *pg_1_3 = buffer.data(pg_1 + 3 * ncomps + c);
        const auto *pg_1_4 = buffer.data(pg_1 + 4 * ncomps + c);
        const auto *pg_1_5 = buffer.data(pg_1 + 5 * ncomps + c);
        const auto *pg_1_6 = buffer.data(pg_1 + 6 * ncomps + c);
        const auto *pg_1_7 = buffer.data(pg_1 + 7 * ncomps + c);
        const auto *pg_1_8 = buffer.data(pg_1 + 8 * ncomps + c);
        const auto *pg_1_9 = buffer.data(pg_1 + 9 * ncomps + c);
        const auto *pg_1_10 = buffer.data(pg_1 + 10 * ncomps + c);
        const auto *pg_1_11 = buffer.data(pg_1 + 11 * ncomps + c);
        const auto *pg_1_12 = buffer.data(pg_1 + 12 * ncomps + c);
        const auto *pg_1_13 = buffer.data(pg_1 + 13 * ncomps + c);
        const auto *pg_1_14 = buffer.data(pg_1 + 14 * ncomps + c);
        const auto *pg_1_15 = buffer.data(pg_1 + 15 * ncomps + c);
        const auto *pg_1_16 = buffer.data(pg_1 + 16 * ncomps + c);
        const auto *pg_1_17 = buffer.data(pg_1 + 17 * ncomps + c);
        const auto *pg_1_18 = buffer.data(pg_1 + 18 * ncomps + c);
        const auto *pg_1_19 = buffer.data(pg_1 + 19 * ncomps + c);
        const auto *pg_1_20 = buffer.data(pg_1 + 20 * ncomps + c);
        const auto *pg_1_21 = buffer.data(pg_1 + 21 * ncomps + c);
        const auto *pg_1_22 = buffer.data(pg_1 + 22 * ncomps + c);
        const auto *pg_1_23 = buffer.data(pg_1 + 23 * ncomps + c);
        const auto *pg_1_24 = buffer.data(pg_1 + 24 * ncomps + c);
        const auto *pg_1_25 = buffer.data(pg_1 + 25 * ncomps + c);
        const auto *pg_1_26 = buffer.data(pg_1 + 26 * ncomps + c);
        const auto *pg_1_27 = buffer.data(pg_1 + 27 * ncomps + c);
        const auto *pg_1_28 = buffer.data(pg_1 + 28 * ncomps + c);
        const auto *pg_1_29 = buffer.data(pg_1 + 29 * ncomps + c);
        const auto *pg_1_30 = buffer.data(pg_1 + 30 * ncomps + c);
        const auto *pg_1_31 = buffer.data(pg_1 + 31 * ncomps + c);
        const auto *pg_1_32 = buffer.data(pg_1 + 32 * ncomps + c);
        const auto *pg_1_33 = buffer.data(pg_1 + 33 * ncomps + c);
        const auto *pg_1_34 = buffer.data(pg_1 + 34 * ncomps + c);
        const auto *pg_1_35 = buffer.data(pg_1 + 35 * ncomps + c);
        const auto *pg_1_36 = buffer.data(pg_1 + 36 * ncomps + c);
        const auto *pg_1_37 = buffer.data(pg_1 + 37 * ncomps + c);
        const auto *pg_1_38 = buffer.data(pg_1 + 38 * ncomps + c);
        const auto *pg_1_39 = buffer.data(pg_1 + 39 * ncomps + c);
        const auto *pg_1_40 = buffer.data(pg_1 + 40 * ncomps + c);
        const auto *pg_1_41 = buffer.data(pg_1 + 41 * ncomps + c);
        const auto *pg_1_42 = buffer.data(pg_1 + 42 * ncomps + c);
        const auto *pg_1_43 = buffer.data(pg_1 + 43 * ncomps + c);
        const auto *pg_1_44 = buffer.data(pg_1 + 44 * ncomps + c);

        const auto *pg_0_0 = buffer.data(pg_0 + 0 * ncomps + c);
        const auto *pg_0_1 = buffer.data(pg_0 + 1 * ncomps + c);
        const auto *pg_0_2 = buffer.data(pg_0 + 2 * ncomps + c);
        const auto *pg_0_3 = buffer.data(pg_0 + 3 * ncomps + c);
        const auto *pg_0_4 = buffer.data(pg_0 + 4 * ncomps + c);
        const auto *pg_0_5 = buffer.data(pg_0 + 5 * ncomps + c);
        const auto *pg_0_6 = buffer.data(pg_0 + 6 * ncomps + c);
        const auto *pg_0_7 = buffer.data(pg_0 + 7 * ncomps + c);
        const auto *pg_0_8 = buffer.data(pg_0 + 8 * ncomps + c);
        const auto *pg_0_9 = buffer.data(pg_0 + 9 * ncomps + c);
        const auto *pg_0_10 = buffer.data(pg_0 + 10 * ncomps + c);
        const auto *pg_0_11 = buffer.data(pg_0 + 11 * ncomps + c);
        const auto *pg_0_12 = buffer.data(pg_0 + 12 * ncomps + c);
        const auto *pg_0_13 = buffer.data(pg_0 + 13 * ncomps + c);
        const auto *pg_0_14 = buffer.data(pg_0 + 14 * ncomps + c);
        const auto *pg_0_15 = buffer.data(pg_0 + 15 * ncomps + c);
        const auto *pg_0_16 = buffer.data(pg_0 + 16 * ncomps + c);
        const auto *pg_0_17 = buffer.data(pg_0 + 17 * ncomps + c);
        const auto *pg_0_18 = buffer.data(pg_0 + 18 * ncomps + c);
        const auto *pg_0_19 = buffer.data(pg_0 + 19 * ncomps + c);
        const auto *pg_0_20 = buffer.data(pg_0 + 20 * ncomps + c);
        const auto *pg_0_21 = buffer.data(pg_0 + 21 * ncomps + c);
        const auto *pg_0_22 = buffer.data(pg_0 + 22 * ncomps + c);
        const auto *pg_0_23 = buffer.data(pg_0 + 23 * ncomps + c);
        const auto *pg_0_24 = buffer.data(pg_0 + 24 * ncomps + c);
        const auto *pg_0_25 = buffer.data(pg_0 + 25 * ncomps + c);
        const auto *pg_0_26 = buffer.data(pg_0 + 26 * ncomps + c);
        const auto *pg_0_27 = buffer.data(pg_0 + 27 * ncomps + c);
        const auto *pg_0_28 = buffer.data(pg_0 + 28 * ncomps + c);
        const auto *pg_0_29 = buffer.data(pg_0 + 29 * ncomps + c);
        const auto *pg_0_30 = buffer.data(pg_0 + 30 * ncomps + c);
        const auto *pg_0_31 = buffer.data(pg_0 + 31 * ncomps + c);
        const auto *pg_0_32 = buffer.data(pg_0 + 32 * ncomps + c);
        const auto *pg_0_33 = buffer.data(pg_0 + 33 * ncomps + c);
        const auto *pg_0_34 = buffer.data(pg_0 + 34 * ncomps + c);
        const auto *pg_0_35 = buffer.data(pg_0 + 35 * ncomps + c);
        const auto *pg_0_36 = buffer.data(pg_0 + 36 * ncomps + c);
        const auto *pg_0_37 = buffer.data(pg_0 + 37 * ncomps + c);
        const auto *pg_0_38 = buffer.data(pg_0 + 38 * ncomps + c);
        const auto *pg_0_39 = buffer.data(pg_0 + 39 * ncomps + c);
        const auto *pg_0_40 = buffer.data(pg_0 + 40 * ncomps + c);
        const auto *pg_0_41 = buffer.data(pg_0 + 41 * ncomps + c);
        const auto *pg_0_42 = buffer.data(pg_0 + 42 * ncomps + c);
        const auto *pg_0_43 = buffer.data(pg_0 + 43 * ncomps + c);
        const auto *pg_0_44 = buffer.data(pg_0 + 44 * ncomps + c);

        const auto *ph_1_0 = buffer.data(ph_1 + 0 * ncomps + c);
        const auto *ph_1_1 = buffer.data(ph_1 + 1 * ncomps + c);
        const auto *ph_1_2 = buffer.data(ph_1 + 2 * ncomps + c);
        const auto *ph_1_3 = buffer.data(ph_1 + 3 * ncomps + c);
        const auto *ph_1_4 = buffer.data(ph_1 + 4 * ncomps + c);
        const auto *ph_1_5 = buffer.data(ph_1 + 5 * ncomps + c);
        const auto *ph_1_6 = buffer.data(ph_1 + 6 * ncomps + c);
        const auto *ph_1_7 = buffer.data(ph_1 + 7 * ncomps + c);
        const auto *ph_1_8 = buffer.data(ph_1 + 8 * ncomps + c);
        const auto *ph_1_9 = buffer.data(ph_1 + 9 * ncomps + c);
        const auto *ph_1_10 = buffer.data(ph_1 + 10 * ncomps + c);
        const auto *ph_1_11 = buffer.data(ph_1 + 11 * ncomps + c);
        const auto *ph_1_12 = buffer.data(ph_1 + 12 * ncomps + c);
        const auto *ph_1_13 = buffer.data(ph_1 + 13 * ncomps + c);
        const auto *ph_1_14 = buffer.data(ph_1 + 14 * ncomps + c);
        const auto *ph_1_21 = buffer.data(ph_1 + 21 * ncomps + c);
        const auto *ph_1_22 = buffer.data(ph_1 + 22 * ncomps + c);
        const auto *ph_1_23 = buffer.data(ph_1 + 23 * ncomps + c);
        const auto *ph_1_24 = buffer.data(ph_1 + 24 * ncomps + c);
        const auto *ph_1_25 = buffer.data(ph_1 + 25 * ncomps + c);
        const auto *ph_1_26 = buffer.data(ph_1 + 26 * ncomps + c);
        const auto *ph_1_27 = buffer.data(ph_1 + 27 * ncomps + c);
        const auto *ph_1_28 = buffer.data(ph_1 + 28 * ncomps + c);
        const auto *ph_1_29 = buffer.data(ph_1 + 29 * ncomps + c);
        const auto *ph_1_30 = buffer.data(ph_1 + 30 * ncomps + c);
        const auto *ph_1_31 = buffer.data(ph_1 + 31 * ncomps + c);
        const auto *ph_1_32 = buffer.data(ph_1 + 32 * ncomps + c);
        const auto *ph_1_33 = buffer.data(ph_1 + 33 * ncomps + c);
        const auto *ph_1_34 = buffer.data(ph_1 + 34 * ncomps + c);
        const auto *ph_1_35 = buffer.data(ph_1 + 35 * ncomps + c);
        const auto *ph_1_36 = buffer.data(ph_1 + 36 * ncomps + c);
        const auto *ph_1_37 = buffer.data(ph_1 + 37 * ncomps + c);
        const auto *ph_1_38 = buffer.data(ph_1 + 38 * ncomps + c);
        const auto *ph_1_39 = buffer.data(ph_1 + 39 * ncomps + c);
        const auto *ph_1_40 = buffer.data(ph_1 + 40 * ncomps + c);
        const auto *ph_1_42 = buffer.data(ph_1 + 42 * ncomps + c);
        const auto *ph_1_43 = buffer.data(ph_1 + 43 * ncomps + c);
        const auto *ph_1_44 = buffer.data(ph_1 + 44 * ncomps + c);
        const auto *ph_1_45 = buffer.data(ph_1 + 45 * ncomps + c);
        const auto *ph_1_46 = buffer.data(ph_1 + 46 * ncomps + c);
        const auto *ph_1_47 = buffer.data(ph_1 + 47 * ncomps + c);
        const auto *ph_1_48 = buffer.data(ph_1 + 48 * ncomps + c);
        const auto *ph_1_49 = buffer.data(ph_1 + 49 * ncomps + c);
        const auto *ph_1_50 = buffer.data(ph_1 + 50 * ncomps + c);
        const auto *ph_1_51 = buffer.data(ph_1 + 51 * ncomps + c);
        const auto *ph_1_52 = buffer.data(ph_1 + 52 * ncomps + c);
        const auto *ph_1_53 = buffer.data(ph_1 + 53 * ncomps + c);
        const auto *ph_1_54 = buffer.data(ph_1 + 54 * ncomps + c);
        const auto *ph_1_55 = buffer.data(ph_1 + 55 * ncomps + c);
        const auto *ph_1_56 = buffer.data(ph_1 + 56 * ncomps + c);
        const auto *ph_1_57 = buffer.data(ph_1 + 57 * ncomps + c);
        const auto *ph_1_58 = buffer.data(ph_1 + 58 * ncomps + c);
        const auto *ph_1_59 = buffer.data(ph_1 + 59 * ncomps + c);
        const auto *ph_1_60 = buffer.data(ph_1 + 60 * ncomps + c);
        const auto *ph_1_61 = buffer.data(ph_1 + 61 * ncomps + c);
        const auto *ph_1_62 = buffer.data(ph_1 + 62 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, ab_x, pg_1_0, pg_1_1, pg_1_2, pg_0_0, pg_0_1, pg_0_2, \
                         ph_1_0, ph_1_1, ph_1_2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * pg_1_0[k]
                     + pg_0_0[k]
                     + ph_1_0[k];

            t_1[k] = -ab_x[k] * pg_1_1[k]
                     + pg_0_1[k]
                     + ph_1_1[k];

            t_2[k] = -ab_x[k] * pg_1_2[k]
                     + pg_0_2[k]
                     + ph_1_2[k];
        }

#pragma omp simd aligned(t_3, t_4, t_5, ab_x, pg_1_3, pg_1_4, pg_1_5, pg_0_3, pg_0_4, pg_0_5, \
                         ph_1_3, ph_1_4, ph_1_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_3[k] = -ab_x[k] * pg_1_3[k]
                     + pg_0_3[k]
                     + ph_1_3[k];

            t_4[k] = -ab_x[k] * pg_1_4[k]
                     + pg_0_4[k]
                     + ph_1_4[k];

            t_5[k] = -ab_x[k] * pg_1_5[k]
                     + pg_0_5[k]
                     + ph_1_5[k];
        }

#pragma omp simd aligned(t_6, t_7, t_8, ab_x, pg_1_6, pg_1_7, pg_1_8, pg_0_6, pg_0_7, pg_0_8, \
                         ph_1_6, ph_1_7, ph_1_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_6[k] = -ab_x[k] * pg_1_6[k]
                     + pg_0_6[k]
                     + ph_1_6[k];

            t_7[k] = -ab_x[k] * pg_1_7[k]
                     + pg_0_7[k]
                     + ph_1_7[k];

            t_8[k] = -ab_x[k] * pg_1_8[k]
                     + pg_0_8[k]
                     + ph_1_8[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, ab_x, pg_1_9, pg_1_10, pg_1_11, pg_0_9, pg_0_10, \
                         pg_0_11, ph_1_9, ph_1_10, ph_1_11 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = -ab_x[k] * pg_1_9[k]
                     + pg_0_9[k]
                     + ph_1_9[k];

            t_10[k] = -ab_x[k] * pg_1_10[k]
                      + pg_0_10[k]
                      + ph_1_10[k];

            t_11[k] = -ab_x[k] * pg_1_11[k]
                      + pg_0_11[k]
                      + ph_1_11[k];
        }

#pragma omp simd aligned(t_12, t_13, t_14, ab_x, pg_1_12, pg_1_13, pg_1_14, pg_0_12, pg_0_13, \
                         pg_0_14, ph_1_12, ph_1_13, ph_1_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_12[k] = -ab_x[k] * pg_1_12[k]
                      + pg_0_12[k]
                      + ph_1_12[k];

            t_13[k] = -ab_x[k] * pg_1_13[k]
                      + pg_0_13[k]
                      + ph_1_13[k];

            t_14[k] = -ab_x[k] * pg_1_14[k]
                      + pg_0_14[k]
                      + ph_1_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, ab_x, pg_1_15, pg_1_16, pg_1_17, pg_0_15, pg_0_16, \
                         pg_0_17, ph_1_21, ph_1_22, ph_1_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * pg_1_15[k]
                      + pg_0_15[k]
                      + ph_1_21[k];

            t_16[k] = -ab_x[k] * pg_1_16[k]
                      + pg_0_16[k]
                      + ph_1_22[k];

            t_17[k] = -ab_x[k] * pg_1_17[k]
                      + pg_0_17[k]
                      + ph_1_23[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, ab_x, pg_1_18, pg_1_19, pg_1_20, pg_0_18, pg_0_19, \
                         pg_0_20, ph_1_24, ph_1_25, ph_1_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = -ab_x[k] * pg_1_18[k]
                      + pg_0_18[k]
                      + ph_1_24[k];

            t_19[k] = -ab_x[k] * pg_1_19[k]
                      + pg_0_19[k]
                      + ph_1_25[k];

            t_20[k] = -ab_x[k] * pg_1_20[k]
                      + pg_0_20[k]
                      + ph_1_26[k];
        }

#pragma omp simd aligned(t_21, t_22, t_23, ab_x, pg_1_21, pg_1_22, pg_1_23, pg_0_21, pg_0_22, \
                         pg_0_23, ph_1_27, ph_1_28, ph_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_21[k] = -ab_x[k] * pg_1_21[k]
                      + pg_0_21[k]
                      + ph_1_27[k];

            t_22[k] = -ab_x[k] * pg_1_22[k]
                      + pg_0_22[k]
                      + ph_1_28[k];

            t_23[k] = -ab_x[k] * pg_1_23[k]
                      + pg_0_23[k]
                      + ph_1_29[k];
        }

#pragma omp simd aligned(t_24, t_25, t_26, ab_x, pg_1_24, pg_1_25, pg_1_26, pg_0_24, pg_0_25, \
                         pg_0_26, ph_1_30, ph_1_31, ph_1_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_24[k] = -ab_x[k] * pg_1_24[k]
                      + pg_0_24[k]
                      + ph_1_30[k];

            t_25[k] = -ab_x[k] * pg_1_25[k]
                      + pg_0_25[k]
                      + ph_1_31[k];

            t_26[k] = -ab_x[k] * pg_1_26[k]
                      + pg_0_26[k]
                      + ph_1_32[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, ab_x, pg_1_27, pg_1_28, pg_1_29, pg_0_27, pg_0_28, \
                         pg_0_29, ph_1_33, ph_1_34, ph_1_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = -ab_x[k] * pg_1_27[k]
                      + pg_0_27[k]
                      + ph_1_33[k];

            t_28[k] = -ab_x[k] * pg_1_28[k]
                      + pg_0_28[k]
                      + ph_1_34[k];

            t_29[k] = -ab_x[k] * pg_1_29[k]
                      + pg_0_29[k]
                      + ph_1_35[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, ab_x, pg_1_30, pg_1_31, pg_1_32, pg_0_30, pg_0_31, \
                         pg_0_32, ph_1_42, ph_1_43, ph_1_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * pg_1_30[k]
                      + pg_0_30[k]
                      + ph_1_42[k];

            t_31[k] = -ab_x[k] * pg_1_31[k]
                      + pg_0_31[k]
                      + ph_1_43[k];

            t_32[k] = -ab_x[k] * pg_1_32[k]
                      + pg_0_32[k]
                      + ph_1_44[k];
        }

#pragma omp simd aligned(t_33, t_34, t_35, ab_x, pg_1_33, pg_1_34, pg_1_35, pg_0_33, pg_0_34, \
                         pg_0_35, ph_1_45, ph_1_46, ph_1_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_33[k] = -ab_x[k] * pg_1_33[k]
                      + pg_0_33[k]
                      + ph_1_45[k];

            t_34[k] = -ab_x[k] * pg_1_34[k]
                      + pg_0_34[k]
                      + ph_1_46[k];

            t_35[k] = -ab_x[k] * pg_1_35[k]
                      + pg_0_35[k]
                      + ph_1_47[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, ab_x, pg_1_36, pg_1_37, pg_1_38, pg_0_36, pg_0_37, \
                         pg_0_38, ph_1_48, ph_1_49, ph_1_50 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = -ab_x[k] * pg_1_36[k]
                      + pg_0_36[k]
                      + ph_1_48[k];

            t_37[k] = -ab_x[k] * pg_1_37[k]
                      + pg_0_37[k]
                      + ph_1_49[k];

            t_38[k] = -ab_x[k] * pg_1_38[k]
                      + pg_0_38[k]
                      + ph_1_50[k];
        }

#pragma omp simd aligned(t_39, t_40, t_41, ab_x, pg_1_39, pg_1_40, pg_1_41, pg_0_39, pg_0_40, \
                         pg_0_41, ph_1_51, ph_1_52, ph_1_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_39[k] = -ab_x[k] * pg_1_39[k]
                      + pg_0_39[k]
                      + ph_1_51[k];

            t_40[k] = -ab_x[k] * pg_1_40[k]
                      + pg_0_40[k]
                      + ph_1_52[k];

            t_41[k] = -ab_x[k] * pg_1_41[k]
                      + pg_0_41[k]
                      + ph_1_53[k];
        }

#pragma omp simd aligned(t_42, t_43, t_44, ab_x, pg_1_42, pg_1_43, pg_1_44, pg_0_42, pg_0_43, \
                         pg_0_44, ph_1_54, ph_1_55, ph_1_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_42[k] = -ab_x[k] * pg_1_42[k]
                      + pg_0_42[k]
                      + ph_1_54[k];

            t_43[k] = -ab_x[k] * pg_1_43[k]
                      + pg_0_43[k]
                      + ph_1_55[k];

            t_44[k] = -ab_x[k] * pg_1_44[k]
                      + pg_0_44[k]
                      + ph_1_56[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_y, pg_1_15, pg_1_16, pg_1_17, \
                         pg_1_18, pg_1_19, ph_1_22, ph_1_24, ph_1_25, ph_1_27, \
                         ph_1_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_y[k] * pg_1_15[k]
                      + ph_1_22[k];

            t_46[k] = -ab_y[k] * pg_1_16[k]
                      + ph_1_24[k];

            t_47[k] = -ab_y[k] * pg_1_17[k]
                      + ph_1_25[k];

            t_48[k] = -ab_y[k] * pg_1_18[k]
                      + ph_1_27[k];

            t_49[k] = -ab_y[k] * pg_1_19[k]
                      + ph_1_28[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_y, pg_1_20, pg_1_21, pg_1_22, \
                         pg_1_23, pg_1_24, ph_1_29, ph_1_31, ph_1_32, ph_1_33, \
                         ph_1_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_y[k] * pg_1_20[k]
                      + ph_1_29[k];

            t_51[k] = -ab_y[k] * pg_1_21[k]
                      + ph_1_31[k];

            t_52[k] = -ab_y[k] * pg_1_22[k]
                      + ph_1_32[k];

            t_53[k] = -ab_y[k] * pg_1_23[k]
                      + ph_1_33[k];

            t_54[k] = -ab_y[k] * pg_1_24[k]
                      + ph_1_34[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_y, pg_1_25, pg_1_26, pg_1_27, \
                         pg_1_28, pg_1_29, ph_1_36, ph_1_37, ph_1_38, ph_1_39, \
                         ph_1_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_y[k] * pg_1_25[k]
                      + ph_1_36[k];

            t_56[k] = -ab_y[k] * pg_1_26[k]
                      + ph_1_37[k];

            t_57[k] = -ab_y[k] * pg_1_27[k]
                      + ph_1_38[k];

            t_58[k] = -ab_y[k] * pg_1_28[k]
                      + ph_1_39[k];

            t_59[k] = -ab_y[k] * pg_1_29[k]
                      + ph_1_40[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_y, pg_1_30, pg_1_31, pg_1_32, \
                         pg_1_33, pg_1_34, ph_1_43, ph_1_45, ph_1_46, ph_1_48, \
                         ph_1_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_y[k] * pg_1_30[k]
                      + ph_1_43[k];

            t_61[k] = -ab_y[k] * pg_1_31[k]
                      + ph_1_45[k];

            t_62[k] = -ab_y[k] * pg_1_32[k]
                      + ph_1_46[k];

            t_63[k] = -ab_y[k] * pg_1_33[k]
                      + ph_1_48[k];

            t_64[k] = -ab_y[k] * pg_1_34[k]
                      + ph_1_49[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_y, pg_1_35, pg_1_36, pg_1_37, \
                         pg_1_38, pg_1_39, ph_1_50, ph_1_52, ph_1_53, ph_1_54, \
                         ph_1_55 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_y[k] * pg_1_35[k]
                      + ph_1_50[k];

            t_66[k] = -ab_y[k] * pg_1_36[k]
                      + ph_1_52[k];

            t_67[k] = -ab_y[k] * pg_1_37[k]
                      + ph_1_53[k];

            t_68[k] = -ab_y[k] * pg_1_38[k]
                      + ph_1_54[k];

            t_69[k] = -ab_y[k] * pg_1_39[k]
                      + ph_1_55[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_y, pg_1_40, pg_1_41, pg_1_42, \
                         pg_1_43, pg_1_44, ph_1_57, ph_1_58, ph_1_59, ph_1_60, \
                         ph_1_61 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_y[k] * pg_1_40[k]
                      + ph_1_57[k];

            t_71[k] = -ab_y[k] * pg_1_41[k]
                      + ph_1_58[k];

            t_72[k] = -ab_y[k] * pg_1_42[k]
                      + ph_1_59[k];

            t_73[k] = -ab_y[k] * pg_1_43[k]
                      + ph_1_60[k];

            t_74[k] = -ab_y[k] * pg_1_44[k]
                      + ph_1_61[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_z, pg_1_30, pg_1_31, pg_1_32, \
                         pg_1_33, pg_1_34, ph_1_44, ph_1_46, ph_1_47, ph_1_49, \
                         ph_1_50 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_z[k] * pg_1_30[k]
                      + ph_1_44[k];

            t_76[k] = -ab_z[k] * pg_1_31[k]
                      + ph_1_46[k];

            t_77[k] = -ab_z[k] * pg_1_32[k]
                      + ph_1_47[k];

            t_78[k] = -ab_z[k] * pg_1_33[k]
                      + ph_1_49[k];

            t_79[k] = -ab_z[k] * pg_1_34[k]
                      + ph_1_50[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_z, pg_1_35, pg_1_36, pg_1_37, \
                         pg_1_38, pg_1_39, ph_1_51, ph_1_53, ph_1_54, ph_1_55, \
                         ph_1_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_z[k] * pg_1_35[k]
                      + ph_1_51[k];

            t_81[k] = -ab_z[k] * pg_1_36[k]
                      + ph_1_53[k];

            t_82[k] = -ab_z[k] * pg_1_37[k]
                      + ph_1_54[k];

            t_83[k] = -ab_z[k] * pg_1_38[k]
                      + ph_1_55[k];

            t_84[k] = -ab_z[k] * pg_1_39[k]
                      + ph_1_56[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_z, pg_1_40, pg_1_41, pg_1_42, \
                         pg_1_43, pg_1_44, ph_1_58, ph_1_59, ph_1_60, ph_1_61, \
                         ph_1_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_z[k] * pg_1_40[k]
                      + ph_1_58[k];

            t_86[k] = -ab_z[k] * pg_1_41[k]
                      + ph_1_59[k];

            t_87[k] = -ab_z[k] * pg_1_42[k]
                      + ph_1_60[k];

            t_88[k] = -ab_z[k] * pg_1_43[k]
                      + ph_1_61[k];

            t_89[k] = -ab_z[k] * pg_1_44[k]
                      + ph_1_62[k];
        }
    }
}

}  // namespace simdtrf
