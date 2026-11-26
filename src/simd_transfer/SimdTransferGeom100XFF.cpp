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


#include "SimdTransferGeom100XFF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_100x_ff_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                      const size_t target, const size_t fd_1, const size_t fd_0,
                                      const size_t gd_1, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *fd_1_0 = buffer.data(fd_1 + 0 * ncomps + c);
        const auto *fd_1_1 = buffer.data(fd_1 + 1 * ncomps + c);
        const auto *fd_1_2 = buffer.data(fd_1 + 2 * ncomps + c);
        const auto *fd_1_3 = buffer.data(fd_1 + 3 * ncomps + c);
        const auto *fd_1_4 = buffer.data(fd_1 + 4 * ncomps + c);
        const auto *fd_1_5 = buffer.data(fd_1 + 5 * ncomps + c);
        const auto *fd_1_6 = buffer.data(fd_1 + 6 * ncomps + c);
        const auto *fd_1_7 = buffer.data(fd_1 + 7 * ncomps + c);
        const auto *fd_1_8 = buffer.data(fd_1 + 8 * ncomps + c);
        const auto *fd_1_9 = buffer.data(fd_1 + 9 * ncomps + c);
        const auto *fd_1_10 = buffer.data(fd_1 + 10 * ncomps + c);
        const auto *fd_1_11 = buffer.data(fd_1 + 11 * ncomps + c);
        const auto *fd_1_12 = buffer.data(fd_1 + 12 * ncomps + c);
        const auto *fd_1_13 = buffer.data(fd_1 + 13 * ncomps + c);
        const auto *fd_1_14 = buffer.data(fd_1 + 14 * ncomps + c);
        const auto *fd_1_15 = buffer.data(fd_1 + 15 * ncomps + c);
        const auto *fd_1_16 = buffer.data(fd_1 + 16 * ncomps + c);
        const auto *fd_1_17 = buffer.data(fd_1 + 17 * ncomps + c);
        const auto *fd_1_18 = buffer.data(fd_1 + 18 * ncomps + c);
        const auto *fd_1_19 = buffer.data(fd_1 + 19 * ncomps + c);
        const auto *fd_1_20 = buffer.data(fd_1 + 20 * ncomps + c);
        const auto *fd_1_21 = buffer.data(fd_1 + 21 * ncomps + c);
        const auto *fd_1_22 = buffer.data(fd_1 + 22 * ncomps + c);
        const auto *fd_1_23 = buffer.data(fd_1 + 23 * ncomps + c);
        const auto *fd_1_24 = buffer.data(fd_1 + 24 * ncomps + c);
        const auto *fd_1_25 = buffer.data(fd_1 + 25 * ncomps + c);
        const auto *fd_1_26 = buffer.data(fd_1 + 26 * ncomps + c);
        const auto *fd_1_27 = buffer.data(fd_1 + 27 * ncomps + c);
        const auto *fd_1_28 = buffer.data(fd_1 + 28 * ncomps + c);
        const auto *fd_1_29 = buffer.data(fd_1 + 29 * ncomps + c);
        const auto *fd_1_30 = buffer.data(fd_1 + 30 * ncomps + c);
        const auto *fd_1_31 = buffer.data(fd_1 + 31 * ncomps + c);
        const auto *fd_1_32 = buffer.data(fd_1 + 32 * ncomps + c);
        const auto *fd_1_33 = buffer.data(fd_1 + 33 * ncomps + c);
        const auto *fd_1_34 = buffer.data(fd_1 + 34 * ncomps + c);
        const auto *fd_1_35 = buffer.data(fd_1 + 35 * ncomps + c);
        const auto *fd_1_36 = buffer.data(fd_1 + 36 * ncomps + c);
        const auto *fd_1_37 = buffer.data(fd_1 + 37 * ncomps + c);
        const auto *fd_1_38 = buffer.data(fd_1 + 38 * ncomps + c);
        const auto *fd_1_39 = buffer.data(fd_1 + 39 * ncomps + c);
        const auto *fd_1_40 = buffer.data(fd_1 + 40 * ncomps + c);
        const auto *fd_1_41 = buffer.data(fd_1 + 41 * ncomps + c);
        const auto *fd_1_42 = buffer.data(fd_1 + 42 * ncomps + c);
        const auto *fd_1_43 = buffer.data(fd_1 + 43 * ncomps + c);
        const auto *fd_1_44 = buffer.data(fd_1 + 44 * ncomps + c);
        const auto *fd_1_45 = buffer.data(fd_1 + 45 * ncomps + c);
        const auto *fd_1_46 = buffer.data(fd_1 + 46 * ncomps + c);
        const auto *fd_1_47 = buffer.data(fd_1 + 47 * ncomps + c);
        const auto *fd_1_48 = buffer.data(fd_1 + 48 * ncomps + c);
        const auto *fd_1_49 = buffer.data(fd_1 + 49 * ncomps + c);
        const auto *fd_1_50 = buffer.data(fd_1 + 50 * ncomps + c);
        const auto *fd_1_51 = buffer.data(fd_1 + 51 * ncomps + c);
        const auto *fd_1_52 = buffer.data(fd_1 + 52 * ncomps + c);
        const auto *fd_1_53 = buffer.data(fd_1 + 53 * ncomps + c);
        const auto *fd_1_54 = buffer.data(fd_1 + 54 * ncomps + c);
        const auto *fd_1_55 = buffer.data(fd_1 + 55 * ncomps + c);
        const auto *fd_1_56 = buffer.data(fd_1 + 56 * ncomps + c);
        const auto *fd_1_57 = buffer.data(fd_1 + 57 * ncomps + c);
        const auto *fd_1_58 = buffer.data(fd_1 + 58 * ncomps + c);
        const auto *fd_1_59 = buffer.data(fd_1 + 59 * ncomps + c);

        const auto *fd_0_0 = buffer.data(fd_0 + 0 * ncomps + c);
        const auto *fd_0_1 = buffer.data(fd_0 + 1 * ncomps + c);
        const auto *fd_0_2 = buffer.data(fd_0 + 2 * ncomps + c);
        const auto *fd_0_3 = buffer.data(fd_0 + 3 * ncomps + c);
        const auto *fd_0_4 = buffer.data(fd_0 + 4 * ncomps + c);
        const auto *fd_0_5 = buffer.data(fd_0 + 5 * ncomps + c);
        const auto *fd_0_6 = buffer.data(fd_0 + 6 * ncomps + c);
        const auto *fd_0_7 = buffer.data(fd_0 + 7 * ncomps + c);
        const auto *fd_0_8 = buffer.data(fd_0 + 8 * ncomps + c);
        const auto *fd_0_9 = buffer.data(fd_0 + 9 * ncomps + c);
        const auto *fd_0_10 = buffer.data(fd_0 + 10 * ncomps + c);
        const auto *fd_0_11 = buffer.data(fd_0 + 11 * ncomps + c);
        const auto *fd_0_12 = buffer.data(fd_0 + 12 * ncomps + c);
        const auto *fd_0_13 = buffer.data(fd_0 + 13 * ncomps + c);
        const auto *fd_0_14 = buffer.data(fd_0 + 14 * ncomps + c);
        const auto *fd_0_15 = buffer.data(fd_0 + 15 * ncomps + c);
        const auto *fd_0_16 = buffer.data(fd_0 + 16 * ncomps + c);
        const auto *fd_0_17 = buffer.data(fd_0 + 17 * ncomps + c);
        const auto *fd_0_18 = buffer.data(fd_0 + 18 * ncomps + c);
        const auto *fd_0_19 = buffer.data(fd_0 + 19 * ncomps + c);
        const auto *fd_0_20 = buffer.data(fd_0 + 20 * ncomps + c);
        const auto *fd_0_21 = buffer.data(fd_0 + 21 * ncomps + c);
        const auto *fd_0_22 = buffer.data(fd_0 + 22 * ncomps + c);
        const auto *fd_0_23 = buffer.data(fd_0 + 23 * ncomps + c);
        const auto *fd_0_24 = buffer.data(fd_0 + 24 * ncomps + c);
        const auto *fd_0_25 = buffer.data(fd_0 + 25 * ncomps + c);
        const auto *fd_0_26 = buffer.data(fd_0 + 26 * ncomps + c);
        const auto *fd_0_27 = buffer.data(fd_0 + 27 * ncomps + c);
        const auto *fd_0_28 = buffer.data(fd_0 + 28 * ncomps + c);
        const auto *fd_0_29 = buffer.data(fd_0 + 29 * ncomps + c);
        const auto *fd_0_30 = buffer.data(fd_0 + 30 * ncomps + c);
        const auto *fd_0_31 = buffer.data(fd_0 + 31 * ncomps + c);
        const auto *fd_0_32 = buffer.data(fd_0 + 32 * ncomps + c);
        const auto *fd_0_33 = buffer.data(fd_0 + 33 * ncomps + c);
        const auto *fd_0_34 = buffer.data(fd_0 + 34 * ncomps + c);
        const auto *fd_0_35 = buffer.data(fd_0 + 35 * ncomps + c);
        const auto *fd_0_36 = buffer.data(fd_0 + 36 * ncomps + c);
        const auto *fd_0_37 = buffer.data(fd_0 + 37 * ncomps + c);
        const auto *fd_0_38 = buffer.data(fd_0 + 38 * ncomps + c);
        const auto *fd_0_39 = buffer.data(fd_0 + 39 * ncomps + c);
        const auto *fd_0_40 = buffer.data(fd_0 + 40 * ncomps + c);
        const auto *fd_0_41 = buffer.data(fd_0 + 41 * ncomps + c);
        const auto *fd_0_42 = buffer.data(fd_0 + 42 * ncomps + c);
        const auto *fd_0_43 = buffer.data(fd_0 + 43 * ncomps + c);
        const auto *fd_0_44 = buffer.data(fd_0 + 44 * ncomps + c);
        const auto *fd_0_45 = buffer.data(fd_0 + 45 * ncomps + c);
        const auto *fd_0_46 = buffer.data(fd_0 + 46 * ncomps + c);
        const auto *fd_0_47 = buffer.data(fd_0 + 47 * ncomps + c);
        const auto *fd_0_48 = buffer.data(fd_0 + 48 * ncomps + c);
        const auto *fd_0_49 = buffer.data(fd_0 + 49 * ncomps + c);
        const auto *fd_0_50 = buffer.data(fd_0 + 50 * ncomps + c);
        const auto *fd_0_51 = buffer.data(fd_0 + 51 * ncomps + c);
        const auto *fd_0_52 = buffer.data(fd_0 + 52 * ncomps + c);
        const auto *fd_0_53 = buffer.data(fd_0 + 53 * ncomps + c);
        const auto *fd_0_54 = buffer.data(fd_0 + 54 * ncomps + c);
        const auto *fd_0_55 = buffer.data(fd_0 + 55 * ncomps + c);
        const auto *fd_0_56 = buffer.data(fd_0 + 56 * ncomps + c);
        const auto *fd_0_57 = buffer.data(fd_0 + 57 * ncomps + c);
        const auto *fd_0_58 = buffer.data(fd_0 + 58 * ncomps + c);
        const auto *fd_0_59 = buffer.data(fd_0 + 59 * ncomps + c);

        const auto *gd_1_0 = buffer.data(gd_1 + 0 * ncomps + c);
        const auto *gd_1_1 = buffer.data(gd_1 + 1 * ncomps + c);
        const auto *gd_1_2 = buffer.data(gd_1 + 2 * ncomps + c);
        const auto *gd_1_3 = buffer.data(gd_1 + 3 * ncomps + c);
        const auto *gd_1_4 = buffer.data(gd_1 + 4 * ncomps + c);
        const auto *gd_1_5 = buffer.data(gd_1 + 5 * ncomps + c);
        const auto *gd_1_6 = buffer.data(gd_1 + 6 * ncomps + c);
        const auto *gd_1_7 = buffer.data(gd_1 + 7 * ncomps + c);
        const auto *gd_1_8 = buffer.data(gd_1 + 8 * ncomps + c);
        const auto *gd_1_9 = buffer.data(gd_1 + 9 * ncomps + c);
        const auto *gd_1_10 = buffer.data(gd_1 + 10 * ncomps + c);
        const auto *gd_1_11 = buffer.data(gd_1 + 11 * ncomps + c);
        const auto *gd_1_12 = buffer.data(gd_1 + 12 * ncomps + c);
        const auto *gd_1_13 = buffer.data(gd_1 + 13 * ncomps + c);
        const auto *gd_1_14 = buffer.data(gd_1 + 14 * ncomps + c);
        const auto *gd_1_15 = buffer.data(gd_1 + 15 * ncomps + c);
        const auto *gd_1_16 = buffer.data(gd_1 + 16 * ncomps + c);
        const auto *gd_1_17 = buffer.data(gd_1 + 17 * ncomps + c);
        const auto *gd_1_18 = buffer.data(gd_1 + 18 * ncomps + c);
        const auto *gd_1_19 = buffer.data(gd_1 + 19 * ncomps + c);
        const auto *gd_1_20 = buffer.data(gd_1 + 20 * ncomps + c);
        const auto *gd_1_21 = buffer.data(gd_1 + 21 * ncomps + c);
        const auto *gd_1_22 = buffer.data(gd_1 + 22 * ncomps + c);
        const auto *gd_1_23 = buffer.data(gd_1 + 23 * ncomps + c);
        const auto *gd_1_24 = buffer.data(gd_1 + 24 * ncomps + c);
        const auto *gd_1_25 = buffer.data(gd_1 + 25 * ncomps + c);
        const auto *gd_1_26 = buffer.data(gd_1 + 26 * ncomps + c);
        const auto *gd_1_27 = buffer.data(gd_1 + 27 * ncomps + c);
        const auto *gd_1_28 = buffer.data(gd_1 + 28 * ncomps + c);
        const auto *gd_1_29 = buffer.data(gd_1 + 29 * ncomps + c);
        const auto *gd_1_30 = buffer.data(gd_1 + 30 * ncomps + c);
        const auto *gd_1_31 = buffer.data(gd_1 + 31 * ncomps + c);
        const auto *gd_1_32 = buffer.data(gd_1 + 32 * ncomps + c);
        const auto *gd_1_33 = buffer.data(gd_1 + 33 * ncomps + c);
        const auto *gd_1_34 = buffer.data(gd_1 + 34 * ncomps + c);
        const auto *gd_1_35 = buffer.data(gd_1 + 35 * ncomps + c);
        const auto *gd_1_36 = buffer.data(gd_1 + 36 * ncomps + c);
        const auto *gd_1_37 = buffer.data(gd_1 + 37 * ncomps + c);
        const auto *gd_1_38 = buffer.data(gd_1 + 38 * ncomps + c);
        const auto *gd_1_39 = buffer.data(gd_1 + 39 * ncomps + c);
        const auto *gd_1_40 = buffer.data(gd_1 + 40 * ncomps + c);
        const auto *gd_1_41 = buffer.data(gd_1 + 41 * ncomps + c);
        const auto *gd_1_42 = buffer.data(gd_1 + 42 * ncomps + c);
        const auto *gd_1_43 = buffer.data(gd_1 + 43 * ncomps + c);
        const auto *gd_1_44 = buffer.data(gd_1 + 44 * ncomps + c);
        const auto *gd_1_45 = buffer.data(gd_1 + 45 * ncomps + c);
        const auto *gd_1_46 = buffer.data(gd_1 + 46 * ncomps + c);
        const auto *gd_1_47 = buffer.data(gd_1 + 47 * ncomps + c);
        const auto *gd_1_48 = buffer.data(gd_1 + 48 * ncomps + c);
        const auto *gd_1_49 = buffer.data(gd_1 + 49 * ncomps + c);
        const auto *gd_1_50 = buffer.data(gd_1 + 50 * ncomps + c);
        const auto *gd_1_51 = buffer.data(gd_1 + 51 * ncomps + c);
        const auto *gd_1_52 = buffer.data(gd_1 + 52 * ncomps + c);
        const auto *gd_1_53 = buffer.data(gd_1 + 53 * ncomps + c);
        const auto *gd_1_54 = buffer.data(gd_1 + 54 * ncomps + c);
        const auto *gd_1_55 = buffer.data(gd_1 + 55 * ncomps + c);
        const auto *gd_1_56 = buffer.data(gd_1 + 56 * ncomps + c);
        const auto *gd_1_57 = buffer.data(gd_1 + 57 * ncomps + c);
        const auto *gd_1_58 = buffer.data(gd_1 + 58 * ncomps + c);
        const auto *gd_1_59 = buffer.data(gd_1 + 59 * ncomps + c);
        const auto *gd_1_63 = buffer.data(gd_1 + 63 * ncomps + c);
        const auto *gd_1_64 = buffer.data(gd_1 + 64 * ncomps + c);
        const auto *gd_1_65 = buffer.data(gd_1 + 65 * ncomps + c);
        const auto *gd_1_69 = buffer.data(gd_1 + 69 * ncomps + c);
        const auto *gd_1_70 = buffer.data(gd_1 + 70 * ncomps + c);
        const auto *gd_1_71 = buffer.data(gd_1 + 71 * ncomps + c);
        const auto *gd_1_75 = buffer.data(gd_1 + 75 * ncomps + c);
        const auto *gd_1_76 = buffer.data(gd_1 + 76 * ncomps + c);
        const auto *gd_1_77 = buffer.data(gd_1 + 77 * ncomps + c);
        const auto *gd_1_81 = buffer.data(gd_1 + 81 * ncomps + c);
        const auto *gd_1_82 = buffer.data(gd_1 + 82 * ncomps + c);
        const auto *gd_1_83 = buffer.data(gd_1 + 83 * ncomps + c);
        const auto *gd_1_89 = buffer.data(gd_1 + 89 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, ab_x, fd_1_0, fd_1_1, fd_1_2, fd_0_0, fd_0_1, fd_0_2, \
                         gd_1_0, gd_1_1, gd_1_2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * fd_1_0[k]
                     + fd_0_0[k]
                     + gd_1_0[k];

            t_1[k] = ab_x[k] * fd_1_1[k]
                     + fd_0_1[k]
                     + gd_1_1[k];

            t_2[k] = ab_x[k] * fd_1_2[k]
                     + fd_0_2[k]
                     + gd_1_2[k];
        }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, ab_x, ab_y, fd_1_3, fd_1_4, fd_1_5, fd_0_3, \
                         fd_0_4, fd_0_5, gd_1_3, gd_1_4, gd_1_5, \
                         gd_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_3[k] = ab_x[k] * fd_1_3[k]
                     + fd_0_3[k]
                     + gd_1_3[k];

            t_4[k] = ab_x[k] * fd_1_4[k]
                     + fd_0_4[k]
                     + gd_1_4[k];

            t_5[k] = ab_x[k] * fd_1_5[k]
                     + fd_0_5[k]
                     + gd_1_5[k];

            t_6[k] = ab_y[k] * fd_1_3[k]
                     + gd_1_9[k];
        }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, fd_1_4, fd_1_5, fd_1_6, \
                         fd_0_6, gd_1_6, gd_1_10, gd_1_11, gd_1_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_7[k] = ab_y[k] * fd_1_4[k]
                     + gd_1_10[k];

            t_8[k] = ab_y[k] * fd_1_5[k]
                     + gd_1_11[k];

            t_9[k] = ab_z[k] * fd_1_5[k]
                     + gd_1_17[k];

            t_10[k] = ab_x[k] * fd_1_6[k]
                      + fd_0_6[k]
                      + gd_1_6[k];
        }

#pragma omp simd aligned(t_11, t_12, t_13, ab_x, fd_1_7, fd_1_8, fd_1_9, fd_0_7, fd_0_8, \
                         fd_0_9, gd_1_7, gd_1_8, gd_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_11[k] = ab_x[k] * fd_1_7[k]
                      + fd_0_7[k]
                      + gd_1_7[k];

            t_12[k] = ab_x[k] * fd_1_8[k]
                      + fd_0_8[k]
                      + gd_1_8[k];

            t_13[k] = ab_x[k] * fd_1_9[k]
                      + fd_0_9[k]
                      + gd_1_9[k];
        }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, ab_x, ab_y, fd_1_9, fd_1_10, fd_1_11, \
                         fd_0_10, fd_0_11, gd_1_10, gd_1_11, gd_1_21, \
                         gd_1_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_14[k] = ab_x[k] * fd_1_10[k]
                      + fd_0_10[k]
                      + gd_1_10[k];

            t_15[k] = ab_x[k] * fd_1_11[k]
                      + fd_0_11[k]
                      + gd_1_11[k];

            t_16[k] = ab_y[k] * fd_1_9[k]
                      + gd_1_21[k];

            t_17[k] = ab_y[k] * fd_1_10[k]
                      + gd_1_22[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, fd_1_11, fd_1_12, fd_1_13, \
                         fd_0_12, fd_0_13, gd_1_12, gd_1_13, gd_1_23, \
                         gd_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = ab_y[k] * fd_1_11[k]
                      + gd_1_23[k];

            t_19[k] = ab_z[k] * fd_1_11[k]
                      + gd_1_29[k];

            t_20[k] = ab_x[k] * fd_1_12[k]
                      + fd_0_12[k]
                      + gd_1_12[k];

            t_21[k] = ab_x[k] * fd_1_13[k]
                      + fd_0_13[k]
                      + gd_1_13[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, ab_x, fd_1_14, fd_1_15, fd_1_16, fd_0_14, fd_0_15, \
                         fd_0_16, gd_1_14, gd_1_15, gd_1_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_x[k] * fd_1_14[k]
                      + fd_0_14[k]
                      + gd_1_14[k];

            t_23[k] = ab_x[k] * fd_1_15[k]
                      + fd_0_15[k]
                      + gd_1_15[k];

            t_24[k] = ab_x[k] * fd_1_16[k]
                      + fd_0_16[k]
                      + gd_1_16[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, ab_x, ab_y, fd_1_15, fd_1_16, fd_1_17, \
                         fd_0_17, gd_1_17, gd_1_27, gd_1_28, gd_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_x[k] * fd_1_17[k]
                      + fd_0_17[k]
                      + gd_1_17[k];

            t_26[k] = ab_y[k] * fd_1_15[k]
                      + gd_1_27[k];

            t_27[k] = ab_y[k] * fd_1_16[k]
                      + gd_1_28[k];

            t_28[k] = ab_y[k] * fd_1_17[k]
                      + gd_1_29[k];
        }

#pragma omp simd aligned(t_29, t_30, t_31, ab_x, ab_z, fd_1_17, fd_1_18, fd_1_19, fd_0_18, \
                         fd_0_19, gd_1_18, gd_1_19, gd_1_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_29[k] = ab_z[k] * fd_1_17[k]
                      + gd_1_35[k];

            t_30[k] = ab_x[k] * fd_1_18[k]
                      + fd_0_18[k]
                      + gd_1_18[k];

            t_31[k] = ab_x[k] * fd_1_19[k]
                      + fd_0_19[k]
                      + gd_1_19[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, ab_x, fd_1_20, fd_1_21, fd_1_22, fd_0_20, fd_0_21, \
                         fd_0_22, gd_1_20, gd_1_21, gd_1_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_x[k] * fd_1_20[k]
                      + fd_0_20[k]
                      + gd_1_20[k];

            t_33[k] = ab_x[k] * fd_1_21[k]
                      + fd_0_21[k]
                      + gd_1_21[k];

            t_34[k] = ab_x[k] * fd_1_22[k]
                      + fd_0_22[k]
                      + gd_1_22[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, ab_x, ab_y, fd_1_21, fd_1_22, fd_1_23, \
                         fd_0_23, gd_1_23, gd_1_39, gd_1_40, gd_1_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * fd_1_23[k]
                      + fd_0_23[k]
                      + gd_1_23[k];

            t_36[k] = ab_y[k] * fd_1_21[k]
                      + gd_1_39[k];

            t_37[k] = ab_y[k] * fd_1_22[k]
                      + gd_1_40[k];

            t_38[k] = ab_y[k] * fd_1_23[k]
                      + gd_1_41[k];
        }

#pragma omp simd aligned(t_39, t_40, t_41, ab_x, ab_z, fd_1_23, fd_1_24, fd_1_25, fd_0_24, \
                         fd_0_25, gd_1_24, gd_1_25, gd_1_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_39[k] = ab_z[k] * fd_1_23[k]
                      + gd_1_47[k];

            t_40[k] = ab_x[k] * fd_1_24[k]
                      + fd_0_24[k]
                      + gd_1_24[k];

            t_41[k] = ab_x[k] * fd_1_25[k]
                      + fd_0_25[k]
                      + gd_1_25[k];
        }

#pragma omp simd aligned(t_42, t_43, t_44, ab_x, fd_1_26, fd_1_27, fd_1_28, fd_0_26, fd_0_27, \
                         fd_0_28, gd_1_26, gd_1_27, gd_1_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_42[k] = ab_x[k] * fd_1_26[k]
                      + fd_0_26[k]
                      + gd_1_26[k];

            t_43[k] = ab_x[k] * fd_1_27[k]
                      + fd_0_27[k]
                      + gd_1_27[k];

            t_44[k] = ab_x[k] * fd_1_28[k]
                      + fd_0_28[k]
                      + gd_1_28[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, fd_1_27, fd_1_28, fd_1_29, \
                         fd_0_29, gd_1_29, gd_1_45, gd_1_46, gd_1_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * fd_1_29[k]
                      + fd_0_29[k]
                      + gd_1_29[k];

            t_46[k] = ab_y[k] * fd_1_27[k]
                      + gd_1_45[k];

            t_47[k] = ab_y[k] * fd_1_28[k]
                      + gd_1_46[k];

            t_48[k] = ab_y[k] * fd_1_29[k]
                      + gd_1_47[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, ab_x, ab_z, fd_1_29, fd_1_30, fd_1_31, fd_0_30, \
                         fd_0_31, gd_1_30, gd_1_31, gd_1_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_z[k] * fd_1_29[k]
                      + gd_1_53[k];

            t_50[k] = ab_x[k] * fd_1_30[k]
                      + fd_0_30[k]
                      + gd_1_30[k];

            t_51[k] = ab_x[k] * fd_1_31[k]
                      + fd_0_31[k]
                      + gd_1_31[k];
        }

#pragma omp simd aligned(t_52, t_53, t_54, ab_x, fd_1_32, fd_1_33, fd_1_34, fd_0_32, fd_0_33, \
                         fd_0_34, gd_1_32, gd_1_33, gd_1_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_52[k] = ab_x[k] * fd_1_32[k]
                      + fd_0_32[k]
                      + gd_1_32[k];

            t_53[k] = ab_x[k] * fd_1_33[k]
                      + fd_0_33[k]
                      + gd_1_33[k];

            t_54[k] = ab_x[k] * fd_1_34[k]
                      + fd_0_34[k]
                      + gd_1_34[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, ab_x, ab_y, fd_1_33, fd_1_34, fd_1_35, \
                         fd_0_35, gd_1_35, gd_1_51, gd_1_52, gd_1_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = ab_x[k] * fd_1_35[k]
                      + fd_0_35[k]
                      + gd_1_35[k];

            t_56[k] = ab_y[k] * fd_1_33[k]
                      + gd_1_51[k];

            t_57[k] = ab_y[k] * fd_1_34[k]
                      + gd_1_52[k];

            t_58[k] = ab_y[k] * fd_1_35[k]
                      + gd_1_53[k];
        }

#pragma omp simd aligned(t_59, t_60, t_61, ab_x, ab_z, fd_1_35, fd_1_36, fd_1_37, fd_0_36, \
                         fd_0_37, gd_1_36, gd_1_37, gd_1_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_59[k] = ab_z[k] * fd_1_35[k]
                      + gd_1_59[k];

            t_60[k] = ab_x[k] * fd_1_36[k]
                      + fd_0_36[k]
                      + gd_1_36[k];

            t_61[k] = ab_x[k] * fd_1_37[k]
                      + fd_0_37[k]
                      + gd_1_37[k];
        }

#pragma omp simd aligned(t_62, t_63, t_64, ab_x, fd_1_38, fd_1_39, fd_1_40, fd_0_38, fd_0_39, \
                         fd_0_40, gd_1_38, gd_1_39, gd_1_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_62[k] = ab_x[k] * fd_1_38[k]
                      + fd_0_38[k]
                      + gd_1_38[k];

            t_63[k] = ab_x[k] * fd_1_39[k]
                      + fd_0_39[k]
                      + gd_1_39[k];

            t_64[k] = ab_x[k] * fd_1_40[k]
                      + fd_0_40[k]
                      + gd_1_40[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, ab_x, ab_y, fd_1_39, fd_1_40, fd_1_41, \
                         fd_0_41, gd_1_41, gd_1_63, gd_1_64, gd_1_65 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * fd_1_41[k]
                      + fd_0_41[k]
                      + gd_1_41[k];

            t_66[k] = ab_y[k] * fd_1_39[k]
                      + gd_1_63[k];

            t_67[k] = ab_y[k] * fd_1_40[k]
                      + gd_1_64[k];

            t_68[k] = ab_y[k] * fd_1_41[k]
                      + gd_1_65[k];
        }

#pragma omp simd aligned(t_69, t_70, t_71, ab_x, ab_z, fd_1_41, fd_1_42, fd_1_43, fd_0_42, \
                         fd_0_43, gd_1_42, gd_1_43, gd_1_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_69[k] = ab_z[k] * fd_1_41[k]
                      + gd_1_71[k];

            t_70[k] = ab_x[k] * fd_1_42[k]
                      + fd_0_42[k]
                      + gd_1_42[k];

            t_71[k] = ab_x[k] * fd_1_43[k]
                      + fd_0_43[k]
                      + gd_1_43[k];
        }

#pragma omp simd aligned(t_72, t_73, t_74, ab_x, fd_1_44, fd_1_45, fd_1_46, fd_0_44, fd_0_45, \
                         fd_0_46, gd_1_44, gd_1_45, gd_1_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_72[k] = ab_x[k] * fd_1_44[k]
                      + fd_0_44[k]
                      + gd_1_44[k];

            t_73[k] = ab_x[k] * fd_1_45[k]
                      + fd_0_45[k]
                      + gd_1_45[k];

            t_74[k] = ab_x[k] * fd_1_46[k]
                      + fd_0_46[k]
                      + gd_1_46[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, ab_x, ab_y, fd_1_45, fd_1_46, fd_1_47, \
                         fd_0_47, gd_1_47, gd_1_69, gd_1_70, gd_1_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = ab_x[k] * fd_1_47[k]
                      + fd_0_47[k]
                      + gd_1_47[k];

            t_76[k] = ab_y[k] * fd_1_45[k]
                      + gd_1_69[k];

            t_77[k] = ab_y[k] * fd_1_46[k]
                      + gd_1_70[k];

            t_78[k] = ab_y[k] * fd_1_47[k]
                      + gd_1_71[k];
        }

#pragma omp simd aligned(t_79, t_80, t_81, ab_x, ab_z, fd_1_47, fd_1_48, fd_1_49, fd_0_48, \
                         fd_0_49, gd_1_48, gd_1_49, gd_1_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_79[k] = ab_z[k] * fd_1_47[k]
                      + gd_1_77[k];

            t_80[k] = ab_x[k] * fd_1_48[k]
                      + fd_0_48[k]
                      + gd_1_48[k];

            t_81[k] = ab_x[k] * fd_1_49[k]
                      + fd_0_49[k]
                      + gd_1_49[k];
        }

#pragma omp simd aligned(t_82, t_83, t_84, ab_x, fd_1_50, fd_1_51, fd_1_52, fd_0_50, fd_0_51, \
                         fd_0_52, gd_1_50, gd_1_51, gd_1_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_82[k] = ab_x[k] * fd_1_50[k]
                      + fd_0_50[k]
                      + gd_1_50[k];

            t_83[k] = ab_x[k] * fd_1_51[k]
                      + fd_0_51[k]
                      + gd_1_51[k];

            t_84[k] = ab_x[k] * fd_1_52[k]
                      + fd_0_52[k]
                      + gd_1_52[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, ab_x, ab_y, fd_1_51, fd_1_52, fd_1_53, \
                         fd_0_53, gd_1_53, gd_1_75, gd_1_76, gd_1_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_x[k] * fd_1_53[k]
                      + fd_0_53[k]
                      + gd_1_53[k];

            t_86[k] = ab_y[k] * fd_1_51[k]
                      + gd_1_75[k];

            t_87[k] = ab_y[k] * fd_1_52[k]
                      + gd_1_76[k];

            t_88[k] = ab_y[k] * fd_1_53[k]
                      + gd_1_77[k];
        }

#pragma omp simd aligned(t_89, t_90, t_91, ab_x, ab_z, fd_1_53, fd_1_54, fd_1_55, fd_0_54, \
                         fd_0_55, gd_1_54, gd_1_55, gd_1_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_89[k] = ab_z[k] * fd_1_53[k]
                      + gd_1_83[k];

            t_90[k] = ab_x[k] * fd_1_54[k]
                      + fd_0_54[k]
                      + gd_1_54[k];

            t_91[k] = ab_x[k] * fd_1_55[k]
                      + fd_0_55[k]
                      + gd_1_55[k];
        }

#pragma omp simd aligned(t_92, t_93, t_94, ab_x, fd_1_56, fd_1_57, fd_1_58, fd_0_56, fd_0_57, \
                         fd_0_58, gd_1_56, gd_1_57, gd_1_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_92[k] = ab_x[k] * fd_1_56[k]
                      + fd_0_56[k]
                      + gd_1_56[k];

            t_93[k] = ab_x[k] * fd_1_57[k]
                      + fd_0_57[k]
                      + gd_1_57[k];

            t_94[k] = ab_x[k] * fd_1_58[k]
                      + fd_0_58[k]
                      + gd_1_58[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, ab_x, ab_y, fd_1_57, fd_1_58, fd_1_59, \
                         fd_0_59, gd_1_59, gd_1_81, gd_1_82, gd_1_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_x[k] * fd_1_59[k]
                      + fd_0_59[k]
                      + gd_1_59[k];

            t_96[k] = ab_y[k] * fd_1_57[k]
                      + gd_1_81[k];

            t_97[k] = ab_y[k] * fd_1_58[k]
                      + gd_1_82[k];

            t_98[k] = ab_y[k] * fd_1_59[k]
                      + gd_1_83[k];
        }

#pragma omp simd aligned(t_99, ab_z, fd_1_59, gd_1_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_z[k] * fd_1_59[k]
                      + gd_1_89[k];
        }
    }
}

}  // namespace simdtrf
