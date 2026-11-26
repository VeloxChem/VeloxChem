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


#include "SimdTransferGeom100XGD.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_100x_gd_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                      const size_t target, const size_t gp_1, const size_t gp_0,
                                      const size_t hp_1, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *gp_1_0 = buffer.data(gp_1 + 0 * ncomps + c);
        const auto *gp_1_1 = buffer.data(gp_1 + 1 * ncomps + c);
        const auto *gp_1_2 = buffer.data(gp_1 + 2 * ncomps + c);
        const auto *gp_1_3 = buffer.data(gp_1 + 3 * ncomps + c);
        const auto *gp_1_4 = buffer.data(gp_1 + 4 * ncomps + c);
        const auto *gp_1_5 = buffer.data(gp_1 + 5 * ncomps + c);
        const auto *gp_1_6 = buffer.data(gp_1 + 6 * ncomps + c);
        const auto *gp_1_7 = buffer.data(gp_1 + 7 * ncomps + c);
        const auto *gp_1_8 = buffer.data(gp_1 + 8 * ncomps + c);
        const auto *gp_1_9 = buffer.data(gp_1 + 9 * ncomps + c);
        const auto *gp_1_10 = buffer.data(gp_1 + 10 * ncomps + c);
        const auto *gp_1_11 = buffer.data(gp_1 + 11 * ncomps + c);
        const auto *gp_1_12 = buffer.data(gp_1 + 12 * ncomps + c);
        const auto *gp_1_13 = buffer.data(gp_1 + 13 * ncomps + c);
        const auto *gp_1_14 = buffer.data(gp_1 + 14 * ncomps + c);
        const auto *gp_1_15 = buffer.data(gp_1 + 15 * ncomps + c);
        const auto *gp_1_16 = buffer.data(gp_1 + 16 * ncomps + c);
        const auto *gp_1_17 = buffer.data(gp_1 + 17 * ncomps + c);
        const auto *gp_1_18 = buffer.data(gp_1 + 18 * ncomps + c);
        const auto *gp_1_19 = buffer.data(gp_1 + 19 * ncomps + c);
        const auto *gp_1_20 = buffer.data(gp_1 + 20 * ncomps + c);
        const auto *gp_1_21 = buffer.data(gp_1 + 21 * ncomps + c);
        const auto *gp_1_22 = buffer.data(gp_1 + 22 * ncomps + c);
        const auto *gp_1_23 = buffer.data(gp_1 + 23 * ncomps + c);
        const auto *gp_1_24 = buffer.data(gp_1 + 24 * ncomps + c);
        const auto *gp_1_25 = buffer.data(gp_1 + 25 * ncomps + c);
        const auto *gp_1_26 = buffer.data(gp_1 + 26 * ncomps + c);
        const auto *gp_1_27 = buffer.data(gp_1 + 27 * ncomps + c);
        const auto *gp_1_28 = buffer.data(gp_1 + 28 * ncomps + c);
        const auto *gp_1_29 = buffer.data(gp_1 + 29 * ncomps + c);
        const auto *gp_1_30 = buffer.data(gp_1 + 30 * ncomps + c);
        const auto *gp_1_31 = buffer.data(gp_1 + 31 * ncomps + c);
        const auto *gp_1_32 = buffer.data(gp_1 + 32 * ncomps + c);
        const auto *gp_1_33 = buffer.data(gp_1 + 33 * ncomps + c);
        const auto *gp_1_34 = buffer.data(gp_1 + 34 * ncomps + c);
        const auto *gp_1_35 = buffer.data(gp_1 + 35 * ncomps + c);
        const auto *gp_1_36 = buffer.data(gp_1 + 36 * ncomps + c);
        const auto *gp_1_37 = buffer.data(gp_1 + 37 * ncomps + c);
        const auto *gp_1_38 = buffer.data(gp_1 + 38 * ncomps + c);
        const auto *gp_1_39 = buffer.data(gp_1 + 39 * ncomps + c);
        const auto *gp_1_40 = buffer.data(gp_1 + 40 * ncomps + c);
        const auto *gp_1_41 = buffer.data(gp_1 + 41 * ncomps + c);
        const auto *gp_1_42 = buffer.data(gp_1 + 42 * ncomps + c);
        const auto *gp_1_43 = buffer.data(gp_1 + 43 * ncomps + c);
        const auto *gp_1_44 = buffer.data(gp_1 + 44 * ncomps + c);

        const auto *gp_0_0 = buffer.data(gp_0 + 0 * ncomps + c);
        const auto *gp_0_1 = buffer.data(gp_0 + 1 * ncomps + c);
        const auto *gp_0_2 = buffer.data(gp_0 + 2 * ncomps + c);
        const auto *gp_0_3 = buffer.data(gp_0 + 3 * ncomps + c);
        const auto *gp_0_4 = buffer.data(gp_0 + 4 * ncomps + c);
        const auto *gp_0_5 = buffer.data(gp_0 + 5 * ncomps + c);
        const auto *gp_0_6 = buffer.data(gp_0 + 6 * ncomps + c);
        const auto *gp_0_7 = buffer.data(gp_0 + 7 * ncomps + c);
        const auto *gp_0_8 = buffer.data(gp_0 + 8 * ncomps + c);
        const auto *gp_0_9 = buffer.data(gp_0 + 9 * ncomps + c);
        const auto *gp_0_10 = buffer.data(gp_0 + 10 * ncomps + c);
        const auto *gp_0_11 = buffer.data(gp_0 + 11 * ncomps + c);
        const auto *gp_0_12 = buffer.data(gp_0 + 12 * ncomps + c);
        const auto *gp_0_13 = buffer.data(gp_0 + 13 * ncomps + c);
        const auto *gp_0_14 = buffer.data(gp_0 + 14 * ncomps + c);
        const auto *gp_0_15 = buffer.data(gp_0 + 15 * ncomps + c);
        const auto *gp_0_16 = buffer.data(gp_0 + 16 * ncomps + c);
        const auto *gp_0_17 = buffer.data(gp_0 + 17 * ncomps + c);
        const auto *gp_0_18 = buffer.data(gp_0 + 18 * ncomps + c);
        const auto *gp_0_19 = buffer.data(gp_0 + 19 * ncomps + c);
        const auto *gp_0_20 = buffer.data(gp_0 + 20 * ncomps + c);
        const auto *gp_0_21 = buffer.data(gp_0 + 21 * ncomps + c);
        const auto *gp_0_22 = buffer.data(gp_0 + 22 * ncomps + c);
        const auto *gp_0_23 = buffer.data(gp_0 + 23 * ncomps + c);
        const auto *gp_0_24 = buffer.data(gp_0 + 24 * ncomps + c);
        const auto *gp_0_25 = buffer.data(gp_0 + 25 * ncomps + c);
        const auto *gp_0_26 = buffer.data(gp_0 + 26 * ncomps + c);
        const auto *gp_0_27 = buffer.data(gp_0 + 27 * ncomps + c);
        const auto *gp_0_28 = buffer.data(gp_0 + 28 * ncomps + c);
        const auto *gp_0_29 = buffer.data(gp_0 + 29 * ncomps + c);
        const auto *gp_0_30 = buffer.data(gp_0 + 30 * ncomps + c);
        const auto *gp_0_31 = buffer.data(gp_0 + 31 * ncomps + c);
        const auto *gp_0_32 = buffer.data(gp_0 + 32 * ncomps + c);
        const auto *gp_0_33 = buffer.data(gp_0 + 33 * ncomps + c);
        const auto *gp_0_34 = buffer.data(gp_0 + 34 * ncomps + c);
        const auto *gp_0_35 = buffer.data(gp_0 + 35 * ncomps + c);
        const auto *gp_0_36 = buffer.data(gp_0 + 36 * ncomps + c);
        const auto *gp_0_37 = buffer.data(gp_0 + 37 * ncomps + c);
        const auto *gp_0_38 = buffer.data(gp_0 + 38 * ncomps + c);
        const auto *gp_0_39 = buffer.data(gp_0 + 39 * ncomps + c);
        const auto *gp_0_40 = buffer.data(gp_0 + 40 * ncomps + c);
        const auto *gp_0_41 = buffer.data(gp_0 + 41 * ncomps + c);
        const auto *gp_0_42 = buffer.data(gp_0 + 42 * ncomps + c);
        const auto *gp_0_43 = buffer.data(gp_0 + 43 * ncomps + c);
        const auto *gp_0_44 = buffer.data(gp_0 + 44 * ncomps + c);

        const auto *hp_1_0 = buffer.data(hp_1 + 0 * ncomps + c);
        const auto *hp_1_1 = buffer.data(hp_1 + 1 * ncomps + c);
        const auto *hp_1_2 = buffer.data(hp_1 + 2 * ncomps + c);
        const auto *hp_1_3 = buffer.data(hp_1 + 3 * ncomps + c);
        const auto *hp_1_4 = buffer.data(hp_1 + 4 * ncomps + c);
        const auto *hp_1_5 = buffer.data(hp_1 + 5 * ncomps + c);
        const auto *hp_1_6 = buffer.data(hp_1 + 6 * ncomps + c);
        const auto *hp_1_7 = buffer.data(hp_1 + 7 * ncomps + c);
        const auto *hp_1_8 = buffer.data(hp_1 + 8 * ncomps + c);
        const auto *hp_1_9 = buffer.data(hp_1 + 9 * ncomps + c);
        const auto *hp_1_10 = buffer.data(hp_1 + 10 * ncomps + c);
        const auto *hp_1_11 = buffer.data(hp_1 + 11 * ncomps + c);
        const auto *hp_1_12 = buffer.data(hp_1 + 12 * ncomps + c);
        const auto *hp_1_13 = buffer.data(hp_1 + 13 * ncomps + c);
        const auto *hp_1_14 = buffer.data(hp_1 + 14 * ncomps + c);
        const auto *hp_1_15 = buffer.data(hp_1 + 15 * ncomps + c);
        const auto *hp_1_16 = buffer.data(hp_1 + 16 * ncomps + c);
        const auto *hp_1_17 = buffer.data(hp_1 + 17 * ncomps + c);
        const auto *hp_1_18 = buffer.data(hp_1 + 18 * ncomps + c);
        const auto *hp_1_19 = buffer.data(hp_1 + 19 * ncomps + c);
        const auto *hp_1_20 = buffer.data(hp_1 + 20 * ncomps + c);
        const auto *hp_1_21 = buffer.data(hp_1 + 21 * ncomps + c);
        const auto *hp_1_22 = buffer.data(hp_1 + 22 * ncomps + c);
        const auto *hp_1_23 = buffer.data(hp_1 + 23 * ncomps + c);
        const auto *hp_1_24 = buffer.data(hp_1 + 24 * ncomps + c);
        const auto *hp_1_25 = buffer.data(hp_1 + 25 * ncomps + c);
        const auto *hp_1_26 = buffer.data(hp_1 + 26 * ncomps + c);
        const auto *hp_1_27 = buffer.data(hp_1 + 27 * ncomps + c);
        const auto *hp_1_28 = buffer.data(hp_1 + 28 * ncomps + c);
        const auto *hp_1_29 = buffer.data(hp_1 + 29 * ncomps + c);
        const auto *hp_1_30 = buffer.data(hp_1 + 30 * ncomps + c);
        const auto *hp_1_31 = buffer.data(hp_1 + 31 * ncomps + c);
        const auto *hp_1_32 = buffer.data(hp_1 + 32 * ncomps + c);
        const auto *hp_1_33 = buffer.data(hp_1 + 33 * ncomps + c);
        const auto *hp_1_34 = buffer.data(hp_1 + 34 * ncomps + c);
        const auto *hp_1_35 = buffer.data(hp_1 + 35 * ncomps + c);
        const auto *hp_1_36 = buffer.data(hp_1 + 36 * ncomps + c);
        const auto *hp_1_37 = buffer.data(hp_1 + 37 * ncomps + c);
        const auto *hp_1_38 = buffer.data(hp_1 + 38 * ncomps + c);
        const auto *hp_1_39 = buffer.data(hp_1 + 39 * ncomps + c);
        const auto *hp_1_40 = buffer.data(hp_1 + 40 * ncomps + c);
        const auto *hp_1_41 = buffer.data(hp_1 + 41 * ncomps + c);
        const auto *hp_1_42 = buffer.data(hp_1 + 42 * ncomps + c);
        const auto *hp_1_43 = buffer.data(hp_1 + 43 * ncomps + c);
        const auto *hp_1_44 = buffer.data(hp_1 + 44 * ncomps + c);
        const auto *hp_1_46 = buffer.data(hp_1 + 46 * ncomps + c);
        const auto *hp_1_47 = buffer.data(hp_1 + 47 * ncomps + c);
        const auto *hp_1_49 = buffer.data(hp_1 + 49 * ncomps + c);
        const auto *hp_1_50 = buffer.data(hp_1 + 50 * ncomps + c);
        const auto *hp_1_52 = buffer.data(hp_1 + 52 * ncomps + c);
        const auto *hp_1_53 = buffer.data(hp_1 + 53 * ncomps + c);
        const auto *hp_1_55 = buffer.data(hp_1 + 55 * ncomps + c);
        const auto *hp_1_56 = buffer.data(hp_1 + 56 * ncomps + c);
        const auto *hp_1_58 = buffer.data(hp_1 + 58 * ncomps + c);
        const auto *hp_1_59 = buffer.data(hp_1 + 59 * ncomps + c);
        const auto *hp_1_62 = buffer.data(hp_1 + 62 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, ab_x, ab_y, gp_1_0, gp_1_1, gp_1_2, gp_0_0, \
                         gp_0_1, gp_0_2, hp_1_0, hp_1_1, hp_1_2, \
                         hp_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * gp_1_0[k]
                     + gp_0_0[k]
                     + hp_1_0[k];

            t_1[k] = ab_x[k] * gp_1_1[k]
                     + gp_0_1[k]
                     + hp_1_1[k];

            t_2[k] = ab_x[k] * gp_1_2[k]
                     + gp_0_2[k]
                     + hp_1_2[k];

            t_3[k] = ab_y[k] * gp_1_1[k]
                     + hp_1_4[k];
        }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, ab_x, ab_y, ab_z, gp_1_2, gp_1_3, gp_1_4, gp_0_3, \
                         gp_0_4, hp_1_3, hp_1_4, hp_1_5, hp_1_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_4[k] = ab_y[k] * gp_1_2[k]
                     + hp_1_5[k];

            t_5[k] = ab_z[k] * gp_1_2[k]
                     + hp_1_8[k];

            t_6[k] = ab_x[k] * gp_1_3[k]
                     + gp_0_3[k]
                     + hp_1_3[k];

            t_7[k] = ab_x[k] * gp_1_4[k]
                     + gp_0_4[k]
                     + hp_1_4[k];
        }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, ab_x, ab_y, ab_z, gp_1_4, gp_1_5, gp_0_5, \
                         hp_1_5, hp_1_10, hp_1_11, hp_1_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_8[k] = ab_x[k] * gp_1_5[k]
                     + gp_0_5[k]
                     + hp_1_5[k];

            t_9[k] = ab_y[k] * gp_1_4[k]
                     + hp_1_10[k];

            t_10[k] = ab_y[k] * gp_1_5[k]
                      + hp_1_11[k];

            t_11[k] = ab_z[k] * gp_1_5[k]
                      + hp_1_14[k];
        }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, ab_x, ab_y, gp_1_6, gp_1_7, gp_1_8, gp_0_6, \
                         gp_0_7, gp_0_8, hp_1_6, hp_1_7, hp_1_8, \
                         hp_1_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_12[k] = ab_x[k] * gp_1_6[k]
                      + gp_0_6[k]
                      + hp_1_6[k];

            t_13[k] = ab_x[k] * gp_1_7[k]
                      + gp_0_7[k]
                      + hp_1_7[k];

            t_14[k] = ab_x[k] * gp_1_8[k]
                      + gp_0_8[k]
                      + hp_1_8[k];

            t_15[k] = ab_y[k] * gp_1_7[k]
                      + hp_1_13[k];
        }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, ab_x, ab_y, ab_z, gp_1_8, gp_1_9, gp_1_10, \
                         gp_0_9, gp_0_10, hp_1_9, hp_1_10, hp_1_14, \
                         hp_1_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_16[k] = ab_y[k] * gp_1_8[k]
                      + hp_1_14[k];

            t_17[k] = ab_z[k] * gp_1_8[k]
                      + hp_1_17[k];

            t_18[k] = ab_x[k] * gp_1_9[k]
                      + gp_0_9[k]
                      + hp_1_9[k];

            t_19[k] = ab_x[k] * gp_1_10[k]
                      + gp_0_10[k]
                      + hp_1_10[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, ab_x, ab_y, ab_z, gp_1_10, gp_1_11, gp_0_11, \
                         hp_1_11, hp_1_19, hp_1_20, hp_1_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * gp_1_11[k]
                      + gp_0_11[k]
                      + hp_1_11[k];

            t_21[k] = ab_y[k] * gp_1_10[k]
                      + hp_1_19[k];

            t_22[k] = ab_y[k] * gp_1_11[k]
                      + hp_1_20[k];

            t_23[k] = ab_z[k] * gp_1_11[k]
                      + hp_1_23[k];
        }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, ab_x, ab_y, gp_1_12, gp_1_13, gp_1_14, \
                         gp_0_12, gp_0_13, gp_0_14, hp_1_12, hp_1_13, hp_1_14, \
                         hp_1_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_24[k] = ab_x[k] * gp_1_12[k]
                      + gp_0_12[k]
                      + hp_1_12[k];

            t_25[k] = ab_x[k] * gp_1_13[k]
                      + gp_0_13[k]
                      + hp_1_13[k];

            t_26[k] = ab_x[k] * gp_1_14[k]
                      + gp_0_14[k]
                      + hp_1_14[k];

            t_27[k] = ab_y[k] * gp_1_13[k]
                      + hp_1_22[k];
        }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, gp_1_14, gp_1_15, gp_1_16, \
                         gp_0_15, gp_0_16, hp_1_15, hp_1_16, hp_1_23, \
                         hp_1_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_28[k] = ab_y[k] * gp_1_14[k]
                      + hp_1_23[k];

            t_29[k] = ab_z[k] * gp_1_14[k]
                      + hp_1_26[k];

            t_30[k] = ab_x[k] * gp_1_15[k]
                      + gp_0_15[k]
                      + hp_1_15[k];

            t_31[k] = ab_x[k] * gp_1_16[k]
                      + gp_0_16[k]
                      + hp_1_16[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, gp_1_16, gp_1_17, gp_0_17, \
                         hp_1_17, hp_1_25, hp_1_26, hp_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_x[k] * gp_1_17[k]
                      + gp_0_17[k]
                      + hp_1_17[k];

            t_33[k] = ab_y[k] * gp_1_16[k]
                      + hp_1_25[k];

            t_34[k] = ab_y[k] * gp_1_17[k]
                      + hp_1_26[k];

            t_35[k] = ab_z[k] * gp_1_17[k]
                      + hp_1_29[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, ab_x, ab_y, gp_1_18, gp_1_19, gp_1_20, \
                         gp_0_18, gp_0_19, gp_0_20, hp_1_18, hp_1_19, hp_1_20, \
                         hp_1_31 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = ab_x[k] * gp_1_18[k]
                      + gp_0_18[k]
                      + hp_1_18[k];

            t_37[k] = ab_x[k] * gp_1_19[k]
                      + gp_0_19[k]
                      + hp_1_19[k];

            t_38[k] = ab_x[k] * gp_1_20[k]
                      + gp_0_20[k]
                      + hp_1_20[k];

            t_39[k] = ab_y[k] * gp_1_19[k]
                      + hp_1_31[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, gp_1_20, gp_1_21, gp_1_22, \
                         gp_0_21, gp_0_22, hp_1_21, hp_1_22, hp_1_32, \
                         hp_1_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_y[k] * gp_1_20[k]
                      + hp_1_32[k];

            t_41[k] = ab_z[k] * gp_1_20[k]
                      + hp_1_35[k];

            t_42[k] = ab_x[k] * gp_1_21[k]
                      + gp_0_21[k]
                      + hp_1_21[k];

            t_43[k] = ab_x[k] * gp_1_22[k]
                      + gp_0_22[k]
                      + hp_1_22[k];
        }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, ab_x, ab_y, ab_z, gp_1_22, gp_1_23, gp_0_23, \
                         hp_1_23, hp_1_34, hp_1_35, hp_1_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_44[k] = ab_x[k] * gp_1_23[k]
                      + gp_0_23[k]
                      + hp_1_23[k];

            t_45[k] = ab_y[k] * gp_1_22[k]
                      + hp_1_34[k];

            t_46[k] = ab_y[k] * gp_1_23[k]
                      + hp_1_35[k];

            t_47[k] = ab_z[k] * gp_1_23[k]
                      + hp_1_38[k];
        }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, ab_x, ab_y, gp_1_24, gp_1_25, gp_1_26, \
                         gp_0_24, gp_0_25, gp_0_26, hp_1_24, hp_1_25, hp_1_26, \
                         hp_1_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_48[k] = ab_x[k] * gp_1_24[k]
                      + gp_0_24[k]
                      + hp_1_24[k];

            t_49[k] = ab_x[k] * gp_1_25[k]
                      + gp_0_25[k]
                      + hp_1_25[k];

            t_50[k] = ab_x[k] * gp_1_26[k]
                      + gp_0_26[k]
                      + hp_1_26[k];

            t_51[k] = ab_y[k] * gp_1_25[k]
                      + hp_1_37[k];
        }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, ab_x, ab_y, ab_z, gp_1_26, gp_1_27, gp_1_28, \
                         gp_0_27, gp_0_28, hp_1_27, hp_1_28, hp_1_38, \
                         hp_1_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_52[k] = ab_y[k] * gp_1_26[k]
                      + hp_1_38[k];

            t_53[k] = ab_z[k] * gp_1_26[k]
                      + hp_1_41[k];

            t_54[k] = ab_x[k] * gp_1_27[k]
                      + gp_0_27[k]
                      + hp_1_27[k];

            t_55[k] = ab_x[k] * gp_1_28[k]
                      + gp_0_28[k]
                      + hp_1_28[k];
        }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, gp_1_28, gp_1_29, gp_0_29, \
                         hp_1_29, hp_1_40, hp_1_41, hp_1_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_56[k] = ab_x[k] * gp_1_29[k]
                      + gp_0_29[k]
                      + hp_1_29[k];

            t_57[k] = ab_y[k] * gp_1_28[k]
                      + hp_1_40[k];

            t_58[k] = ab_y[k] * gp_1_29[k]
                      + hp_1_41[k];

            t_59[k] = ab_z[k] * gp_1_29[k]
                      + hp_1_44[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, ab_x, ab_y, gp_1_30, gp_1_31, gp_1_32, \
                         gp_0_30, gp_0_31, gp_0_32, hp_1_30, hp_1_31, hp_1_32, \
                         hp_1_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * gp_1_30[k]
                      + gp_0_30[k]
                      + hp_1_30[k];

            t_61[k] = ab_x[k] * gp_1_31[k]
                      + gp_0_31[k]
                      + hp_1_31[k];

            t_62[k] = ab_x[k] * gp_1_32[k]
                      + gp_0_32[k]
                      + hp_1_32[k];

            t_63[k] = ab_y[k] * gp_1_31[k]
                      + hp_1_46[k];
        }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, ab_x, ab_y, ab_z, gp_1_32, gp_1_33, gp_1_34, \
                         gp_0_33, gp_0_34, hp_1_33, hp_1_34, hp_1_47, \
                         hp_1_50 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_64[k] = ab_y[k] * gp_1_32[k]
                      + hp_1_47[k];

            t_65[k] = ab_z[k] * gp_1_32[k]
                      + hp_1_50[k];

            t_66[k] = ab_x[k] * gp_1_33[k]
                      + gp_0_33[k]
                      + hp_1_33[k];

            t_67[k] = ab_x[k] * gp_1_34[k]
                      + gp_0_34[k]
                      + hp_1_34[k];
        }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, ab_x, ab_y, ab_z, gp_1_34, gp_1_35, gp_0_35, \
                         hp_1_35, hp_1_49, hp_1_50, hp_1_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_68[k] = ab_x[k] * gp_1_35[k]
                      + gp_0_35[k]
                      + hp_1_35[k];

            t_69[k] = ab_y[k] * gp_1_34[k]
                      + hp_1_49[k];

            t_70[k] = ab_y[k] * gp_1_35[k]
                      + hp_1_50[k];

            t_71[k] = ab_z[k] * gp_1_35[k]
                      + hp_1_53[k];
        }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, ab_x, ab_y, gp_1_36, gp_1_37, gp_1_38, \
                         gp_0_36, gp_0_37, gp_0_38, hp_1_36, hp_1_37, hp_1_38, \
                         hp_1_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_72[k] = ab_x[k] * gp_1_36[k]
                      + gp_0_36[k]
                      + hp_1_36[k];

            t_73[k] = ab_x[k] * gp_1_37[k]
                      + gp_0_37[k]
                      + hp_1_37[k];

            t_74[k] = ab_x[k] * gp_1_38[k]
                      + gp_0_38[k]
                      + hp_1_38[k];

            t_75[k] = ab_y[k] * gp_1_37[k]
                      + hp_1_52[k];
        }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, ab_x, ab_y, ab_z, gp_1_38, gp_1_39, gp_1_40, \
                         gp_0_39, gp_0_40, hp_1_39, hp_1_40, hp_1_53, \
                         hp_1_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_76[k] = ab_y[k] * gp_1_38[k]
                      + hp_1_53[k];

            t_77[k] = ab_z[k] * gp_1_38[k]
                      + hp_1_56[k];

            t_78[k] = ab_x[k] * gp_1_39[k]
                      + gp_0_39[k]
                      + hp_1_39[k];

            t_79[k] = ab_x[k] * gp_1_40[k]
                      + gp_0_40[k]
                      + hp_1_40[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, ab_x, ab_y, ab_z, gp_1_40, gp_1_41, gp_0_41, \
                         hp_1_41, hp_1_55, hp_1_56, hp_1_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_x[k] * gp_1_41[k]
                      + gp_0_41[k]
                      + hp_1_41[k];

            t_81[k] = ab_y[k] * gp_1_40[k]
                      + hp_1_55[k];

            t_82[k] = ab_y[k] * gp_1_41[k]
                      + hp_1_56[k];

            t_83[k] = ab_z[k] * gp_1_41[k]
                      + hp_1_59[k];
        }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, ab_x, ab_y, gp_1_42, gp_1_43, gp_1_44, \
                         gp_0_42, gp_0_43, gp_0_44, hp_1_42, hp_1_43, hp_1_44, \
                         hp_1_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_84[k] = ab_x[k] * gp_1_42[k]
                      + gp_0_42[k]
                      + hp_1_42[k];

            t_85[k] = ab_x[k] * gp_1_43[k]
                      + gp_0_43[k]
                      + hp_1_43[k];

            t_86[k] = ab_x[k] * gp_1_44[k]
                      + gp_0_44[k]
                      + hp_1_44[k];

            t_87[k] = ab_y[k] * gp_1_43[k]
                      + hp_1_58[k];
        }

#pragma omp simd aligned(t_88, t_89, ab_y, ab_z, gp_1_44, hp_1_59, \
                         hp_1_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_88[k] = ab_y[k] * gp_1_44[k]
                      + hp_1_59[k];

            t_89[k] = ab_z[k] * gp_1_44[k]
                      + hp_1_62[k];
        }
    }
}

}  // namespace simdtrf
