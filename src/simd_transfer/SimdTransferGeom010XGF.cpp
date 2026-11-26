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


#include "SimdTransferGeom010XGF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_geom_010x_gf_out_of_second_piece0(CSimdMatrix &buffer,
                                              const CSimdMatrix &coordinates,
                                              const size_t target, const size_t ff_1,
                                              const size_t ff_0, const size_t fg_1,
                                              const size_t ncomps, const size_t nmax) -> void
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

        const auto *ff_1_0 = buffer.data(ff_1 + 0 * ncomps + c);
        const auto *ff_1_1 = buffer.data(ff_1 + 1 * ncomps + c);
        const auto *ff_1_2 = buffer.data(ff_1 + 2 * ncomps + c);
        const auto *ff_1_3 = buffer.data(ff_1 + 3 * ncomps + c);
        const auto *ff_1_4 = buffer.data(ff_1 + 4 * ncomps + c);
        const auto *ff_1_5 = buffer.data(ff_1 + 5 * ncomps + c);
        const auto *ff_1_6 = buffer.data(ff_1 + 6 * ncomps + c);
        const auto *ff_1_7 = buffer.data(ff_1 + 7 * ncomps + c);
        const auto *ff_1_8 = buffer.data(ff_1 + 8 * ncomps + c);
        const auto *ff_1_9 = buffer.data(ff_1 + 9 * ncomps + c);
        const auto *ff_1_10 = buffer.data(ff_1 + 10 * ncomps + c);
        const auto *ff_1_11 = buffer.data(ff_1 + 11 * ncomps + c);
        const auto *ff_1_12 = buffer.data(ff_1 + 12 * ncomps + c);
        const auto *ff_1_13 = buffer.data(ff_1 + 13 * ncomps + c);
        const auto *ff_1_14 = buffer.data(ff_1 + 14 * ncomps + c);
        const auto *ff_1_15 = buffer.data(ff_1 + 15 * ncomps + c);
        const auto *ff_1_16 = buffer.data(ff_1 + 16 * ncomps + c);
        const auto *ff_1_17 = buffer.data(ff_1 + 17 * ncomps + c);
        const auto *ff_1_18 = buffer.data(ff_1 + 18 * ncomps + c);
        const auto *ff_1_19 = buffer.data(ff_1 + 19 * ncomps + c);
        const auto *ff_1_20 = buffer.data(ff_1 + 20 * ncomps + c);
        const auto *ff_1_21 = buffer.data(ff_1 + 21 * ncomps + c);
        const auto *ff_1_22 = buffer.data(ff_1 + 22 * ncomps + c);
        const auto *ff_1_23 = buffer.data(ff_1 + 23 * ncomps + c);
        const auto *ff_1_24 = buffer.data(ff_1 + 24 * ncomps + c);
        const auto *ff_1_25 = buffer.data(ff_1 + 25 * ncomps + c);
        const auto *ff_1_26 = buffer.data(ff_1 + 26 * ncomps + c);
        const auto *ff_1_27 = buffer.data(ff_1 + 27 * ncomps + c);
        const auto *ff_1_28 = buffer.data(ff_1 + 28 * ncomps + c);
        const auto *ff_1_29 = buffer.data(ff_1 + 29 * ncomps + c);
        const auto *ff_1_30 = buffer.data(ff_1 + 30 * ncomps + c);
        const auto *ff_1_31 = buffer.data(ff_1 + 31 * ncomps + c);
        const auto *ff_1_32 = buffer.data(ff_1 + 32 * ncomps + c);
        const auto *ff_1_33 = buffer.data(ff_1 + 33 * ncomps + c);
        const auto *ff_1_34 = buffer.data(ff_1 + 34 * ncomps + c);
        const auto *ff_1_35 = buffer.data(ff_1 + 35 * ncomps + c);
        const auto *ff_1_36 = buffer.data(ff_1 + 36 * ncomps + c);
        const auto *ff_1_37 = buffer.data(ff_1 + 37 * ncomps + c);
        const auto *ff_1_38 = buffer.data(ff_1 + 38 * ncomps + c);
        const auto *ff_1_39 = buffer.data(ff_1 + 39 * ncomps + c);
        const auto *ff_1_40 = buffer.data(ff_1 + 40 * ncomps + c);
        const auto *ff_1_41 = buffer.data(ff_1 + 41 * ncomps + c);
        const auto *ff_1_42 = buffer.data(ff_1 + 42 * ncomps + c);
        const auto *ff_1_43 = buffer.data(ff_1 + 43 * ncomps + c);
        const auto *ff_1_44 = buffer.data(ff_1 + 44 * ncomps + c);
        const auto *ff_1_45 = buffer.data(ff_1 + 45 * ncomps + c);
        const auto *ff_1_46 = buffer.data(ff_1 + 46 * ncomps + c);
        const auto *ff_1_47 = buffer.data(ff_1 + 47 * ncomps + c);
        const auto *ff_1_48 = buffer.data(ff_1 + 48 * ncomps + c);
        const auto *ff_1_49 = buffer.data(ff_1 + 49 * ncomps + c);
        const auto *ff_1_50 = buffer.data(ff_1 + 50 * ncomps + c);
        const auto *ff_1_51 = buffer.data(ff_1 + 51 * ncomps + c);
        const auto *ff_1_52 = buffer.data(ff_1 + 52 * ncomps + c);
        const auto *ff_1_53 = buffer.data(ff_1 + 53 * ncomps + c);
        const auto *ff_1_54 = buffer.data(ff_1 + 54 * ncomps + c);
        const auto *ff_1_55 = buffer.data(ff_1 + 55 * ncomps + c);
        const auto *ff_1_56 = buffer.data(ff_1 + 56 * ncomps + c);
        const auto *ff_1_57 = buffer.data(ff_1 + 57 * ncomps + c);
        const auto *ff_1_58 = buffer.data(ff_1 + 58 * ncomps + c);
        const auto *ff_1_59 = buffer.data(ff_1 + 59 * ncomps + c);
        const auto *ff_1_60 = buffer.data(ff_1 + 60 * ncomps + c);
        const auto *ff_1_61 = buffer.data(ff_1 + 61 * ncomps + c);
        const auto *ff_1_62 = buffer.data(ff_1 + 62 * ncomps + c);
        const auto *ff_1_63 = buffer.data(ff_1 + 63 * ncomps + c);
        const auto *ff_1_64 = buffer.data(ff_1 + 64 * ncomps + c);
        const auto *ff_1_65 = buffer.data(ff_1 + 65 * ncomps + c);
        const auto *ff_1_66 = buffer.data(ff_1 + 66 * ncomps + c);
        const auto *ff_1_67 = buffer.data(ff_1 + 67 * ncomps + c);
        const auto *ff_1_68 = buffer.data(ff_1 + 68 * ncomps + c);
        const auto *ff_1_69 = buffer.data(ff_1 + 69 * ncomps + c);
        const auto *ff_1_70 = buffer.data(ff_1 + 70 * ncomps + c);
        const auto *ff_1_71 = buffer.data(ff_1 + 71 * ncomps + c);
        const auto *ff_1_72 = buffer.data(ff_1 + 72 * ncomps + c);
        const auto *ff_1_73 = buffer.data(ff_1 + 73 * ncomps + c);
        const auto *ff_1_74 = buffer.data(ff_1 + 74 * ncomps + c);
        const auto *ff_1_75 = buffer.data(ff_1 + 75 * ncomps + c);
        const auto *ff_1_76 = buffer.data(ff_1 + 76 * ncomps + c);
        const auto *ff_1_77 = buffer.data(ff_1 + 77 * ncomps + c);
        const auto *ff_1_78 = buffer.data(ff_1 + 78 * ncomps + c);
        const auto *ff_1_79 = buffer.data(ff_1 + 79 * ncomps + c);
        const auto *ff_1_80 = buffer.data(ff_1 + 80 * ncomps + c);
        const auto *ff_1_81 = buffer.data(ff_1 + 81 * ncomps + c);
        const auto *ff_1_82 = buffer.data(ff_1 + 82 * ncomps + c);
        const auto *ff_1_83 = buffer.data(ff_1 + 83 * ncomps + c);
        const auto *ff_1_84 = buffer.data(ff_1 + 84 * ncomps + c);
        const auto *ff_1_85 = buffer.data(ff_1 + 85 * ncomps + c);
        const auto *ff_1_86 = buffer.data(ff_1 + 86 * ncomps + c);
        const auto *ff_1_87 = buffer.data(ff_1 + 87 * ncomps + c);
        const auto *ff_1_88 = buffer.data(ff_1 + 88 * ncomps + c);
        const auto *ff_1_89 = buffer.data(ff_1 + 89 * ncomps + c);
        const auto *ff_1_90 = buffer.data(ff_1 + 90 * ncomps + c);
        const auto *ff_1_91 = buffer.data(ff_1 + 91 * ncomps + c);
        const auto *ff_1_92 = buffer.data(ff_1 + 92 * ncomps + c);
        const auto *ff_1_93 = buffer.data(ff_1 + 93 * ncomps + c);
        const auto *ff_1_94 = buffer.data(ff_1 + 94 * ncomps + c);
        const auto *ff_1_95 = buffer.data(ff_1 + 95 * ncomps + c);
        const auto *ff_1_96 = buffer.data(ff_1 + 96 * ncomps + c);
        const auto *ff_1_97 = buffer.data(ff_1 + 97 * ncomps + c);
        const auto *ff_1_98 = buffer.data(ff_1 + 98 * ncomps + c);
        const auto *ff_1_99 = buffer.data(ff_1 + 99 * ncomps + c);

        const auto *ff_0_0 = buffer.data(ff_0 + 0 * ncomps + c);
        const auto *ff_0_1 = buffer.data(ff_0 + 1 * ncomps + c);
        const auto *ff_0_2 = buffer.data(ff_0 + 2 * ncomps + c);
        const auto *ff_0_3 = buffer.data(ff_0 + 3 * ncomps + c);
        const auto *ff_0_4 = buffer.data(ff_0 + 4 * ncomps + c);
        const auto *ff_0_5 = buffer.data(ff_0 + 5 * ncomps + c);
        const auto *ff_0_6 = buffer.data(ff_0 + 6 * ncomps + c);
        const auto *ff_0_7 = buffer.data(ff_0 + 7 * ncomps + c);
        const auto *ff_0_8 = buffer.data(ff_0 + 8 * ncomps + c);
        const auto *ff_0_9 = buffer.data(ff_0 + 9 * ncomps + c);
        const auto *ff_0_10 = buffer.data(ff_0 + 10 * ncomps + c);
        const auto *ff_0_11 = buffer.data(ff_0 + 11 * ncomps + c);
        const auto *ff_0_12 = buffer.data(ff_0 + 12 * ncomps + c);
        const auto *ff_0_13 = buffer.data(ff_0 + 13 * ncomps + c);
        const auto *ff_0_14 = buffer.data(ff_0 + 14 * ncomps + c);
        const auto *ff_0_15 = buffer.data(ff_0 + 15 * ncomps + c);
        const auto *ff_0_16 = buffer.data(ff_0 + 16 * ncomps + c);
        const auto *ff_0_17 = buffer.data(ff_0 + 17 * ncomps + c);
        const auto *ff_0_18 = buffer.data(ff_0 + 18 * ncomps + c);
        const auto *ff_0_19 = buffer.data(ff_0 + 19 * ncomps + c);
        const auto *ff_0_20 = buffer.data(ff_0 + 20 * ncomps + c);
        const auto *ff_0_21 = buffer.data(ff_0 + 21 * ncomps + c);
        const auto *ff_0_22 = buffer.data(ff_0 + 22 * ncomps + c);
        const auto *ff_0_23 = buffer.data(ff_0 + 23 * ncomps + c);
        const auto *ff_0_24 = buffer.data(ff_0 + 24 * ncomps + c);
        const auto *ff_0_25 = buffer.data(ff_0 + 25 * ncomps + c);
        const auto *ff_0_26 = buffer.data(ff_0 + 26 * ncomps + c);
        const auto *ff_0_27 = buffer.data(ff_0 + 27 * ncomps + c);
        const auto *ff_0_28 = buffer.data(ff_0 + 28 * ncomps + c);
        const auto *ff_0_29 = buffer.data(ff_0 + 29 * ncomps + c);
        const auto *ff_0_30 = buffer.data(ff_0 + 30 * ncomps + c);
        const auto *ff_0_31 = buffer.data(ff_0 + 31 * ncomps + c);
        const auto *ff_0_32 = buffer.data(ff_0 + 32 * ncomps + c);
        const auto *ff_0_33 = buffer.data(ff_0 + 33 * ncomps + c);
        const auto *ff_0_34 = buffer.data(ff_0 + 34 * ncomps + c);
        const auto *ff_0_35 = buffer.data(ff_0 + 35 * ncomps + c);
        const auto *ff_0_36 = buffer.data(ff_0 + 36 * ncomps + c);
        const auto *ff_0_37 = buffer.data(ff_0 + 37 * ncomps + c);
        const auto *ff_0_38 = buffer.data(ff_0 + 38 * ncomps + c);
        const auto *ff_0_39 = buffer.data(ff_0 + 39 * ncomps + c);
        const auto *ff_0_40 = buffer.data(ff_0 + 40 * ncomps + c);
        const auto *ff_0_41 = buffer.data(ff_0 + 41 * ncomps + c);
        const auto *ff_0_42 = buffer.data(ff_0 + 42 * ncomps + c);
        const auto *ff_0_43 = buffer.data(ff_0 + 43 * ncomps + c);
        const auto *ff_0_44 = buffer.data(ff_0 + 44 * ncomps + c);
        const auto *ff_0_45 = buffer.data(ff_0 + 45 * ncomps + c);
        const auto *ff_0_46 = buffer.data(ff_0 + 46 * ncomps + c);
        const auto *ff_0_47 = buffer.data(ff_0 + 47 * ncomps + c);
        const auto *ff_0_48 = buffer.data(ff_0 + 48 * ncomps + c);
        const auto *ff_0_49 = buffer.data(ff_0 + 49 * ncomps + c);
        const auto *ff_0_50 = buffer.data(ff_0 + 50 * ncomps + c);
        const auto *ff_0_51 = buffer.data(ff_0 + 51 * ncomps + c);
        const auto *ff_0_52 = buffer.data(ff_0 + 52 * ncomps + c);
        const auto *ff_0_53 = buffer.data(ff_0 + 53 * ncomps + c);
        const auto *ff_0_54 = buffer.data(ff_0 + 54 * ncomps + c);
        const auto *ff_0_55 = buffer.data(ff_0 + 55 * ncomps + c);
        const auto *ff_0_56 = buffer.data(ff_0 + 56 * ncomps + c);
        const auto *ff_0_57 = buffer.data(ff_0 + 57 * ncomps + c);
        const auto *ff_0_58 = buffer.data(ff_0 + 58 * ncomps + c);
        const auto *ff_0_59 = buffer.data(ff_0 + 59 * ncomps + c);
        const auto *ff_0_60 = buffer.data(ff_0 + 60 * ncomps + c);
        const auto *ff_0_61 = buffer.data(ff_0 + 61 * ncomps + c);
        const auto *ff_0_62 = buffer.data(ff_0 + 62 * ncomps + c);
        const auto *ff_0_63 = buffer.data(ff_0 + 63 * ncomps + c);
        const auto *ff_0_64 = buffer.data(ff_0 + 64 * ncomps + c);
        const auto *ff_0_65 = buffer.data(ff_0 + 65 * ncomps + c);
        const auto *ff_0_66 = buffer.data(ff_0 + 66 * ncomps + c);
        const auto *ff_0_67 = buffer.data(ff_0 + 67 * ncomps + c);
        const auto *ff_0_68 = buffer.data(ff_0 + 68 * ncomps + c);
        const auto *ff_0_69 = buffer.data(ff_0 + 69 * ncomps + c);
        const auto *ff_0_70 = buffer.data(ff_0 + 70 * ncomps + c);
        const auto *ff_0_71 = buffer.data(ff_0 + 71 * ncomps + c);
        const auto *ff_0_72 = buffer.data(ff_0 + 72 * ncomps + c);
        const auto *ff_0_73 = buffer.data(ff_0 + 73 * ncomps + c);
        const auto *ff_0_74 = buffer.data(ff_0 + 74 * ncomps + c);
        const auto *ff_0_75 = buffer.data(ff_0 + 75 * ncomps + c);
        const auto *ff_0_76 = buffer.data(ff_0 + 76 * ncomps + c);
        const auto *ff_0_77 = buffer.data(ff_0 + 77 * ncomps + c);
        const auto *ff_0_78 = buffer.data(ff_0 + 78 * ncomps + c);
        const auto *ff_0_79 = buffer.data(ff_0 + 79 * ncomps + c);
        const auto *ff_0_80 = buffer.data(ff_0 + 80 * ncomps + c);
        const auto *ff_0_81 = buffer.data(ff_0 + 81 * ncomps + c);
        const auto *ff_0_82 = buffer.data(ff_0 + 82 * ncomps + c);
        const auto *ff_0_83 = buffer.data(ff_0 + 83 * ncomps + c);
        const auto *ff_0_84 = buffer.data(ff_0 + 84 * ncomps + c);
        const auto *ff_0_85 = buffer.data(ff_0 + 85 * ncomps + c);
        const auto *ff_0_86 = buffer.data(ff_0 + 86 * ncomps + c);
        const auto *ff_0_87 = buffer.data(ff_0 + 87 * ncomps + c);
        const auto *ff_0_88 = buffer.data(ff_0 + 88 * ncomps + c);
        const auto *ff_0_89 = buffer.data(ff_0 + 89 * ncomps + c);
        const auto *ff_0_90 = buffer.data(ff_0 + 90 * ncomps + c);
        const auto *ff_0_91 = buffer.data(ff_0 + 91 * ncomps + c);
        const auto *ff_0_92 = buffer.data(ff_0 + 92 * ncomps + c);
        const auto *ff_0_93 = buffer.data(ff_0 + 93 * ncomps + c);
        const auto *ff_0_94 = buffer.data(ff_0 + 94 * ncomps + c);
        const auto *ff_0_95 = buffer.data(ff_0 + 95 * ncomps + c);
        const auto *ff_0_96 = buffer.data(ff_0 + 96 * ncomps + c);
        const auto *ff_0_97 = buffer.data(ff_0 + 97 * ncomps + c);
        const auto *ff_0_98 = buffer.data(ff_0 + 98 * ncomps + c);
        const auto *ff_0_99 = buffer.data(ff_0 + 99 * ncomps + c);

        const auto *fg_1_0 = buffer.data(fg_1 + 0 * ncomps + c);
        const auto *fg_1_1 = buffer.data(fg_1 + 1 * ncomps + c);
        const auto *fg_1_2 = buffer.data(fg_1 + 2 * ncomps + c);
        const auto *fg_1_3 = buffer.data(fg_1 + 3 * ncomps + c);
        const auto *fg_1_4 = buffer.data(fg_1 + 4 * ncomps + c);
        const auto *fg_1_5 = buffer.data(fg_1 + 5 * ncomps + c);
        const auto *fg_1_6 = buffer.data(fg_1 + 6 * ncomps + c);
        const auto *fg_1_7 = buffer.data(fg_1 + 7 * ncomps + c);
        const auto *fg_1_8 = buffer.data(fg_1 + 8 * ncomps + c);
        const auto *fg_1_9 = buffer.data(fg_1 + 9 * ncomps + c);
        const auto *fg_1_15 = buffer.data(fg_1 + 15 * ncomps + c);
        const auto *fg_1_16 = buffer.data(fg_1 + 16 * ncomps + c);
        const auto *fg_1_17 = buffer.data(fg_1 + 17 * ncomps + c);
        const auto *fg_1_18 = buffer.data(fg_1 + 18 * ncomps + c);
        const auto *fg_1_19 = buffer.data(fg_1 + 19 * ncomps + c);
        const auto *fg_1_20 = buffer.data(fg_1 + 20 * ncomps + c);
        const auto *fg_1_21 = buffer.data(fg_1 + 21 * ncomps + c);
        const auto *fg_1_22 = buffer.data(fg_1 + 22 * ncomps + c);
        const auto *fg_1_23 = buffer.data(fg_1 + 23 * ncomps + c);
        const auto *fg_1_24 = buffer.data(fg_1 + 24 * ncomps + c);
        const auto *fg_1_30 = buffer.data(fg_1 + 30 * ncomps + c);
        const auto *fg_1_31 = buffer.data(fg_1 + 31 * ncomps + c);
        const auto *fg_1_32 = buffer.data(fg_1 + 32 * ncomps + c);
        const auto *fg_1_33 = buffer.data(fg_1 + 33 * ncomps + c);
        const auto *fg_1_34 = buffer.data(fg_1 + 34 * ncomps + c);
        const auto *fg_1_35 = buffer.data(fg_1 + 35 * ncomps + c);
        const auto *fg_1_36 = buffer.data(fg_1 + 36 * ncomps + c);
        const auto *fg_1_37 = buffer.data(fg_1 + 37 * ncomps + c);
        const auto *fg_1_38 = buffer.data(fg_1 + 38 * ncomps + c);
        const auto *fg_1_39 = buffer.data(fg_1 + 39 * ncomps + c);
        const auto *fg_1_45 = buffer.data(fg_1 + 45 * ncomps + c);
        const auto *fg_1_46 = buffer.data(fg_1 + 46 * ncomps + c);
        const auto *fg_1_47 = buffer.data(fg_1 + 47 * ncomps + c);
        const auto *fg_1_48 = buffer.data(fg_1 + 48 * ncomps + c);
        const auto *fg_1_49 = buffer.data(fg_1 + 49 * ncomps + c);
        const auto *fg_1_50 = buffer.data(fg_1 + 50 * ncomps + c);
        const auto *fg_1_51 = buffer.data(fg_1 + 51 * ncomps + c);
        const auto *fg_1_52 = buffer.data(fg_1 + 52 * ncomps + c);
        const auto *fg_1_53 = buffer.data(fg_1 + 53 * ncomps + c);
        const auto *fg_1_54 = buffer.data(fg_1 + 54 * ncomps + c);
        const auto *fg_1_60 = buffer.data(fg_1 + 60 * ncomps + c);
        const auto *fg_1_61 = buffer.data(fg_1 + 61 * ncomps + c);
        const auto *fg_1_62 = buffer.data(fg_1 + 62 * ncomps + c);
        const auto *fg_1_63 = buffer.data(fg_1 + 63 * ncomps + c);
        const auto *fg_1_64 = buffer.data(fg_1 + 64 * ncomps + c);
        const auto *fg_1_65 = buffer.data(fg_1 + 65 * ncomps + c);
        const auto *fg_1_66 = buffer.data(fg_1 + 66 * ncomps + c);
        const auto *fg_1_67 = buffer.data(fg_1 + 67 * ncomps + c);
        const auto *fg_1_68 = buffer.data(fg_1 + 68 * ncomps + c);
        const auto *fg_1_69 = buffer.data(fg_1 + 69 * ncomps + c);
        const auto *fg_1_75 = buffer.data(fg_1 + 75 * ncomps + c);
        const auto *fg_1_76 = buffer.data(fg_1 + 76 * ncomps + c);
        const auto *fg_1_77 = buffer.data(fg_1 + 77 * ncomps + c);
        const auto *fg_1_78 = buffer.data(fg_1 + 78 * ncomps + c);
        const auto *fg_1_79 = buffer.data(fg_1 + 79 * ncomps + c);
        const auto *fg_1_80 = buffer.data(fg_1 + 80 * ncomps + c);
        const auto *fg_1_81 = buffer.data(fg_1 + 81 * ncomps + c);
        const auto *fg_1_82 = buffer.data(fg_1 + 82 * ncomps + c);
        const auto *fg_1_83 = buffer.data(fg_1 + 83 * ncomps + c);
        const auto *fg_1_84 = buffer.data(fg_1 + 84 * ncomps + c);
        const auto *fg_1_90 = buffer.data(fg_1 + 90 * ncomps + c);
        const auto *fg_1_91 = buffer.data(fg_1 + 91 * ncomps + c);
        const auto *fg_1_92 = buffer.data(fg_1 + 92 * ncomps + c);
        const auto *fg_1_93 = buffer.data(fg_1 + 93 * ncomps + c);
        const auto *fg_1_94 = buffer.data(fg_1 + 94 * ncomps + c);
        const auto *fg_1_95 = buffer.data(fg_1 + 95 * ncomps + c);
        const auto *fg_1_96 = buffer.data(fg_1 + 96 * ncomps + c);
        const auto *fg_1_97 = buffer.data(fg_1 + 97 * ncomps + c);
        const auto *fg_1_98 = buffer.data(fg_1 + 98 * ncomps + c);
        const auto *fg_1_99 = buffer.data(fg_1 + 99 * ncomps + c);
        const auto *fg_1_100 = buffer.data(fg_1 + 100 * ncomps + c);
        const auto *fg_1_101 = buffer.data(fg_1 + 101 * ncomps + c);
        const auto *fg_1_105 = buffer.data(fg_1 + 105 * ncomps + c);
        const auto *fg_1_106 = buffer.data(fg_1 + 106 * ncomps + c);
        const auto *fg_1_107 = buffer.data(fg_1 + 107 * ncomps + c);
        const auto *fg_1_108 = buffer.data(fg_1 + 108 * ncomps + c);
        const auto *fg_1_109 = buffer.data(fg_1 + 109 * ncomps + c);
        const auto *fg_1_110 = buffer.data(fg_1 + 110 * ncomps + c);
        const auto *fg_1_111 = buffer.data(fg_1 + 111 * ncomps + c);
        const auto *fg_1_112 = buffer.data(fg_1 + 112 * ncomps + c);
        const auto *fg_1_113 = buffer.data(fg_1 + 113 * ncomps + c);
        const auto *fg_1_114 = buffer.data(fg_1 + 114 * ncomps + c);
        const auto *fg_1_120 = buffer.data(fg_1 + 120 * ncomps + c);
        const auto *fg_1_121 = buffer.data(fg_1 + 121 * ncomps + c);
        const auto *fg_1_122 = buffer.data(fg_1 + 122 * ncomps + c);
        const auto *fg_1_123 = buffer.data(fg_1 + 123 * ncomps + c);
        const auto *fg_1_124 = buffer.data(fg_1 + 124 * ncomps + c);
        const auto *fg_1_125 = buffer.data(fg_1 + 125 * ncomps + c);
        const auto *fg_1_126 = buffer.data(fg_1 + 126 * ncomps + c);
        const auto *fg_1_127 = buffer.data(fg_1 + 127 * ncomps + c);
        const auto *fg_1_128 = buffer.data(fg_1 + 128 * ncomps + c);
        const auto *fg_1_129 = buffer.data(fg_1 + 129 * ncomps + c);
        const auto *fg_1_135 = buffer.data(fg_1 + 135 * ncomps + c);
        const auto *fg_1_136 = buffer.data(fg_1 + 136 * ncomps + c);
        const auto *fg_1_137 = buffer.data(fg_1 + 137 * ncomps + c);
        const auto *fg_1_138 = buffer.data(fg_1 + 138 * ncomps + c);
        const auto *fg_1_139 = buffer.data(fg_1 + 139 * ncomps + c);
        const auto *fg_1_140 = buffer.data(fg_1 + 140 * ncomps + c);
        const auto *fg_1_141 = buffer.data(fg_1 + 141 * ncomps + c);
        const auto *fg_1_142 = buffer.data(fg_1 + 142 * ncomps + c);
        const auto *fg_1_143 = buffer.data(fg_1 + 143 * ncomps + c);
        const auto *fg_1_144 = buffer.data(fg_1 + 144 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, ab_x, ff_1_0, ff_1_1, ff_1_2, ff_0_0, ff_0_1, ff_0_2, \
                         fg_1_0, fg_1_1, fg_1_2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * ff_1_0[k]
                     + ff_0_0[k]
                     + fg_1_0[k];

            t_1[k] = -ab_x[k] * ff_1_1[k]
                     + ff_0_1[k]
                     + fg_1_1[k];

            t_2[k] = -ab_x[k] * ff_1_2[k]
                     + ff_0_2[k]
                     + fg_1_2[k];
        }

#pragma omp simd aligned(t_3, t_4, t_5, ab_x, ff_1_3, ff_1_4, ff_1_5, ff_0_3, ff_0_4, ff_0_5, \
                         fg_1_3, fg_1_4, fg_1_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_3[k] = -ab_x[k] * ff_1_3[k]
                     + ff_0_3[k]
                     + fg_1_3[k];

            t_4[k] = -ab_x[k] * ff_1_4[k]
                     + ff_0_4[k]
                     + fg_1_4[k];

            t_5[k] = -ab_x[k] * ff_1_5[k]
                     + ff_0_5[k]
                     + fg_1_5[k];
        }

#pragma omp simd aligned(t_6, t_7, t_8, ab_x, ff_1_6, ff_1_7, ff_1_8, ff_0_6, ff_0_7, ff_0_8, \
                         fg_1_6, fg_1_7, fg_1_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_6[k] = -ab_x[k] * ff_1_6[k]
                     + ff_0_6[k]
                     + fg_1_6[k];

            t_7[k] = -ab_x[k] * ff_1_7[k]
                     + ff_0_7[k]
                     + fg_1_7[k];

            t_8[k] = -ab_x[k] * ff_1_8[k]
                     + ff_0_8[k]
                     + fg_1_8[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, ab_x, ff_1_9, ff_1_10, ff_1_11, ff_0_9, ff_0_10, \
                         ff_0_11, fg_1_9, fg_1_15, fg_1_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = -ab_x[k] * ff_1_9[k]
                     + ff_0_9[k]
                     + fg_1_9[k];

            t_10[k] = -ab_x[k] * ff_1_10[k]
                      + ff_0_10[k]
                      + fg_1_15[k];

            t_11[k] = -ab_x[k] * ff_1_11[k]
                      + ff_0_11[k]
                      + fg_1_16[k];
        }

#pragma omp simd aligned(t_12, t_13, t_14, ab_x, ff_1_12, ff_1_13, ff_1_14, ff_0_12, ff_0_13, \
                         ff_0_14, fg_1_17, fg_1_18, fg_1_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_12[k] = -ab_x[k] * ff_1_12[k]
                      + ff_0_12[k]
                      + fg_1_17[k];

            t_13[k] = -ab_x[k] * ff_1_13[k]
                      + ff_0_13[k]
                      + fg_1_18[k];

            t_14[k] = -ab_x[k] * ff_1_14[k]
                      + ff_0_14[k]
                      + fg_1_19[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, ab_x, ff_1_15, ff_1_16, ff_1_17, ff_0_15, ff_0_16, \
                         ff_0_17, fg_1_20, fg_1_21, fg_1_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * ff_1_15[k]
                      + ff_0_15[k]
                      + fg_1_20[k];

            t_16[k] = -ab_x[k] * ff_1_16[k]
                      + ff_0_16[k]
                      + fg_1_21[k];

            t_17[k] = -ab_x[k] * ff_1_17[k]
                      + ff_0_17[k]
                      + fg_1_22[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, ab_x, ff_1_18, ff_1_19, ff_1_20, ff_0_18, ff_0_19, \
                         ff_0_20, fg_1_23, fg_1_24, fg_1_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = -ab_x[k] * ff_1_18[k]
                      + ff_0_18[k]
                      + fg_1_23[k];

            t_19[k] = -ab_x[k] * ff_1_19[k]
                      + ff_0_19[k]
                      + fg_1_24[k];

            t_20[k] = -ab_x[k] * ff_1_20[k]
                      + ff_0_20[k]
                      + fg_1_30[k];
        }

#pragma omp simd aligned(t_21, t_22, t_23, ab_x, ff_1_21, ff_1_22, ff_1_23, ff_0_21, ff_0_22, \
                         ff_0_23, fg_1_31, fg_1_32, fg_1_33 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_21[k] = -ab_x[k] * ff_1_21[k]
                      + ff_0_21[k]
                      + fg_1_31[k];

            t_22[k] = -ab_x[k] * ff_1_22[k]
                      + ff_0_22[k]
                      + fg_1_32[k];

            t_23[k] = -ab_x[k] * ff_1_23[k]
                      + ff_0_23[k]
                      + fg_1_33[k];
        }

#pragma omp simd aligned(t_24, t_25, t_26, ab_x, ff_1_24, ff_1_25, ff_1_26, ff_0_24, ff_0_25, \
                         ff_0_26, fg_1_34, fg_1_35, fg_1_36 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_24[k] = -ab_x[k] * ff_1_24[k]
                      + ff_0_24[k]
                      + fg_1_34[k];

            t_25[k] = -ab_x[k] * ff_1_25[k]
                      + ff_0_25[k]
                      + fg_1_35[k];

            t_26[k] = -ab_x[k] * ff_1_26[k]
                      + ff_0_26[k]
                      + fg_1_36[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, ab_x, ff_1_27, ff_1_28, ff_1_29, ff_0_27, ff_0_28, \
                         ff_0_29, fg_1_37, fg_1_38, fg_1_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = -ab_x[k] * ff_1_27[k]
                      + ff_0_27[k]
                      + fg_1_37[k];

            t_28[k] = -ab_x[k] * ff_1_28[k]
                      + ff_0_28[k]
                      + fg_1_38[k];

            t_29[k] = -ab_x[k] * ff_1_29[k]
                      + ff_0_29[k]
                      + fg_1_39[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, ab_x, ff_1_30, ff_1_31, ff_1_32, ff_0_30, ff_0_31, \
                         ff_0_32, fg_1_45, fg_1_46, fg_1_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * ff_1_30[k]
                      + ff_0_30[k]
                      + fg_1_45[k];

            t_31[k] = -ab_x[k] * ff_1_31[k]
                      + ff_0_31[k]
                      + fg_1_46[k];

            t_32[k] = -ab_x[k] * ff_1_32[k]
                      + ff_0_32[k]
                      + fg_1_47[k];
        }

#pragma omp simd aligned(t_33, t_34, t_35, ab_x, ff_1_33, ff_1_34, ff_1_35, ff_0_33, ff_0_34, \
                         ff_0_35, fg_1_48, fg_1_49, fg_1_50 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_33[k] = -ab_x[k] * ff_1_33[k]
                      + ff_0_33[k]
                      + fg_1_48[k];

            t_34[k] = -ab_x[k] * ff_1_34[k]
                      + ff_0_34[k]
                      + fg_1_49[k];

            t_35[k] = -ab_x[k] * ff_1_35[k]
                      + ff_0_35[k]
                      + fg_1_50[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, ab_x, ff_1_36, ff_1_37, ff_1_38, ff_0_36, ff_0_37, \
                         ff_0_38, fg_1_51, fg_1_52, fg_1_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = -ab_x[k] * ff_1_36[k]
                      + ff_0_36[k]
                      + fg_1_51[k];

            t_37[k] = -ab_x[k] * ff_1_37[k]
                      + ff_0_37[k]
                      + fg_1_52[k];

            t_38[k] = -ab_x[k] * ff_1_38[k]
                      + ff_0_38[k]
                      + fg_1_53[k];
        }

#pragma omp simd aligned(t_39, t_40, t_41, ab_x, ff_1_39, ff_1_40, ff_1_41, ff_0_39, ff_0_40, \
                         ff_0_41, fg_1_54, fg_1_60, fg_1_61 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_39[k] = -ab_x[k] * ff_1_39[k]
                      + ff_0_39[k]
                      + fg_1_54[k];

            t_40[k] = -ab_x[k] * ff_1_40[k]
                      + ff_0_40[k]
                      + fg_1_60[k];

            t_41[k] = -ab_x[k] * ff_1_41[k]
                      + ff_0_41[k]
                      + fg_1_61[k];
        }

#pragma omp simd aligned(t_42, t_43, t_44, ab_x, ff_1_42, ff_1_43, ff_1_44, ff_0_42, ff_0_43, \
                         ff_0_44, fg_1_62, fg_1_63, fg_1_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_42[k] = -ab_x[k] * ff_1_42[k]
                      + ff_0_42[k]
                      + fg_1_62[k];

            t_43[k] = -ab_x[k] * ff_1_43[k]
                      + ff_0_43[k]
                      + fg_1_63[k];

            t_44[k] = -ab_x[k] * ff_1_44[k]
                      + ff_0_44[k]
                      + fg_1_64[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, ab_x, ff_1_45, ff_1_46, ff_1_47, ff_0_45, ff_0_46, \
                         ff_0_47, fg_1_65, fg_1_66, fg_1_67 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * ff_1_45[k]
                      + ff_0_45[k]
                      + fg_1_65[k];

            t_46[k] = -ab_x[k] * ff_1_46[k]
                      + ff_0_46[k]
                      + fg_1_66[k];

            t_47[k] = -ab_x[k] * ff_1_47[k]
                      + ff_0_47[k]
                      + fg_1_67[k];
        }

#pragma omp simd aligned(t_48, t_49, t_50, ab_x, ff_1_48, ff_1_49, ff_1_50, ff_0_48, ff_0_49, \
                         ff_0_50, fg_1_68, fg_1_69, fg_1_75 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_48[k] = -ab_x[k] * ff_1_48[k]
                      + ff_0_48[k]
                      + fg_1_68[k];

            t_49[k] = -ab_x[k] * ff_1_49[k]
                      + ff_0_49[k]
                      + fg_1_69[k];

            t_50[k] = -ab_x[k] * ff_1_50[k]
                      + ff_0_50[k]
                      + fg_1_75[k];
        }

#pragma omp simd aligned(t_51, t_52, t_53, ab_x, ff_1_51, ff_1_52, ff_1_53, ff_0_51, ff_0_52, \
                         ff_0_53, fg_1_76, fg_1_77, fg_1_78 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_51[k] = -ab_x[k] * ff_1_51[k]
                      + ff_0_51[k]
                      + fg_1_76[k];

            t_52[k] = -ab_x[k] * ff_1_52[k]
                      + ff_0_52[k]
                      + fg_1_77[k];

            t_53[k] = -ab_x[k] * ff_1_53[k]
                      + ff_0_53[k]
                      + fg_1_78[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, ab_x, ff_1_54, ff_1_55, ff_1_56, ff_0_54, ff_0_55, \
                         ff_0_56, fg_1_79, fg_1_80, fg_1_81 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = -ab_x[k] * ff_1_54[k]
                      + ff_0_54[k]
                      + fg_1_79[k];

            t_55[k] = -ab_x[k] * ff_1_55[k]
                      + ff_0_55[k]
                      + fg_1_80[k];

            t_56[k] = -ab_x[k] * ff_1_56[k]
                      + ff_0_56[k]
                      + fg_1_81[k];
        }

#pragma omp simd aligned(t_57, t_58, t_59, ab_x, ff_1_57, ff_1_58, ff_1_59, ff_0_57, ff_0_58, \
                         ff_0_59, fg_1_82, fg_1_83, fg_1_84 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_57[k] = -ab_x[k] * ff_1_57[k]
                      + ff_0_57[k]
                      + fg_1_82[k];

            t_58[k] = -ab_x[k] * ff_1_58[k]
                      + ff_0_58[k]
                      + fg_1_83[k];

            t_59[k] = -ab_x[k] * ff_1_59[k]
                      + ff_0_59[k]
                      + fg_1_84[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, ab_x, ff_1_60, ff_1_61, ff_1_62, ff_0_60, ff_0_61, \
                         ff_0_62, fg_1_90, fg_1_91, fg_1_92 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * ff_1_60[k]
                      + ff_0_60[k]
                      + fg_1_90[k];

            t_61[k] = -ab_x[k] * ff_1_61[k]
                      + ff_0_61[k]
                      + fg_1_91[k];

            t_62[k] = -ab_x[k] * ff_1_62[k]
                      + ff_0_62[k]
                      + fg_1_92[k];
        }

#pragma omp simd aligned(t_63, t_64, t_65, ab_x, ff_1_63, ff_1_64, ff_1_65, ff_0_63, ff_0_64, \
                         ff_0_65, fg_1_93, fg_1_94, fg_1_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_63[k] = -ab_x[k] * ff_1_63[k]
                      + ff_0_63[k]
                      + fg_1_93[k];

            t_64[k] = -ab_x[k] * ff_1_64[k]
                      + ff_0_64[k]
                      + fg_1_94[k];

            t_65[k] = -ab_x[k] * ff_1_65[k]
                      + ff_0_65[k]
                      + fg_1_95[k];
        }

#pragma omp simd aligned(t_66, t_67, t_68, ab_x, ff_1_66, ff_1_67, ff_1_68, ff_0_66, ff_0_67, \
                         ff_0_68, fg_1_96, fg_1_97, fg_1_98 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_66[k] = -ab_x[k] * ff_1_66[k]
                      + ff_0_66[k]
                      + fg_1_96[k];

            t_67[k] = -ab_x[k] * ff_1_67[k]
                      + ff_0_67[k]
                      + fg_1_97[k];

            t_68[k] = -ab_x[k] * ff_1_68[k]
                      + ff_0_68[k]
                      + fg_1_98[k];
        }

#pragma omp simd aligned(t_69, t_70, t_71, ab_x, ff_1_69, ff_1_70, ff_1_71, ff_0_69, ff_0_70, \
                         ff_0_71, fg_1_99, fg_1_105, fg_1_106 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_69[k] = -ab_x[k] * ff_1_69[k]
                      + ff_0_69[k]
                      + fg_1_99[k];

            t_70[k] = -ab_x[k] * ff_1_70[k]
                      + ff_0_70[k]
                      + fg_1_105[k];

            t_71[k] = -ab_x[k] * ff_1_71[k]
                      + ff_0_71[k]
                      + fg_1_106[k];
        }

#pragma omp simd aligned(t_72, t_73, t_74, ab_x, ff_1_72, ff_1_73, ff_1_74, ff_0_72, ff_0_73, \
                         ff_0_74, fg_1_107, fg_1_108, fg_1_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_72[k] = -ab_x[k] * ff_1_72[k]
                      + ff_0_72[k]
                      + fg_1_107[k];

            t_73[k] = -ab_x[k] * ff_1_73[k]
                      + ff_0_73[k]
                      + fg_1_108[k];

            t_74[k] = -ab_x[k] * ff_1_74[k]
                      + ff_0_74[k]
                      + fg_1_109[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, ab_x, ff_1_75, ff_1_76, ff_1_77, ff_0_75, ff_0_76, \
                         ff_0_77, fg_1_110, fg_1_111, fg_1_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * ff_1_75[k]
                      + ff_0_75[k]
                      + fg_1_110[k];

            t_76[k] = -ab_x[k] * ff_1_76[k]
                      + ff_0_76[k]
                      + fg_1_111[k];

            t_77[k] = -ab_x[k] * ff_1_77[k]
                      + ff_0_77[k]
                      + fg_1_112[k];
        }

#pragma omp simd aligned(t_78, t_79, t_80, ab_x, ff_1_78, ff_1_79, ff_1_80, ff_0_78, ff_0_79, \
                         ff_0_80, fg_1_113, fg_1_114, fg_1_120 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_78[k] = -ab_x[k] * ff_1_78[k]
                      + ff_0_78[k]
                      + fg_1_113[k];

            t_79[k] = -ab_x[k] * ff_1_79[k]
                      + ff_0_79[k]
                      + fg_1_114[k];

            t_80[k] = -ab_x[k] * ff_1_80[k]
                      + ff_0_80[k]
                      + fg_1_120[k];
        }

#pragma omp simd aligned(t_81, t_82, t_83, ab_x, ff_1_81, ff_1_82, ff_1_83, ff_0_81, ff_0_82, \
                         ff_0_83, fg_1_121, fg_1_122, fg_1_123 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_81[k] = -ab_x[k] * ff_1_81[k]
                      + ff_0_81[k]
                      + fg_1_121[k];

            t_82[k] = -ab_x[k] * ff_1_82[k]
                      + ff_0_82[k]
                      + fg_1_122[k];

            t_83[k] = -ab_x[k] * ff_1_83[k]
                      + ff_0_83[k]
                      + fg_1_123[k];
        }

#pragma omp simd aligned(t_84, t_85, t_86, ab_x, ff_1_84, ff_1_85, ff_1_86, ff_0_84, ff_0_85, \
                         ff_0_86, fg_1_124, fg_1_125, fg_1_126 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_84[k] = -ab_x[k] * ff_1_84[k]
                      + ff_0_84[k]
                      + fg_1_124[k];

            t_85[k] = -ab_x[k] * ff_1_85[k]
                      + ff_0_85[k]
                      + fg_1_125[k];

            t_86[k] = -ab_x[k] * ff_1_86[k]
                      + ff_0_86[k]
                      + fg_1_126[k];
        }

#pragma omp simd aligned(t_87, t_88, t_89, ab_x, ff_1_87, ff_1_88, ff_1_89, ff_0_87, ff_0_88, \
                         ff_0_89, fg_1_127, fg_1_128, fg_1_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_87[k] = -ab_x[k] * ff_1_87[k]
                      + ff_0_87[k]
                      + fg_1_127[k];

            t_88[k] = -ab_x[k] * ff_1_88[k]
                      + ff_0_88[k]
                      + fg_1_128[k];

            t_89[k] = -ab_x[k] * ff_1_89[k]
                      + ff_0_89[k]
                      + fg_1_129[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, ab_x, ff_1_90, ff_1_91, ff_1_92, ff_0_90, ff_0_91, \
                         ff_0_92, fg_1_135, fg_1_136, fg_1_137 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_x[k] * ff_1_90[k]
                      + ff_0_90[k]
                      + fg_1_135[k];

            t_91[k] = -ab_x[k] * ff_1_91[k]
                      + ff_0_91[k]
                      + fg_1_136[k];

            t_92[k] = -ab_x[k] * ff_1_92[k]
                      + ff_0_92[k]
                      + fg_1_137[k];
        }

#pragma omp simd aligned(t_93, t_94, t_95, ab_x, ff_1_93, ff_1_94, ff_1_95, ff_0_93, ff_0_94, \
                         ff_0_95, fg_1_138, fg_1_139, fg_1_140 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_93[k] = -ab_x[k] * ff_1_93[k]
                      + ff_0_93[k]
                      + fg_1_138[k];

            t_94[k] = -ab_x[k] * ff_1_94[k]
                      + ff_0_94[k]
                      + fg_1_139[k];

            t_95[k] = -ab_x[k] * ff_1_95[k]
                      + ff_0_95[k]
                      + fg_1_140[k];
        }

#pragma omp simd aligned(t_96, t_97, t_98, ab_x, ff_1_96, ff_1_97, ff_1_98, ff_0_96, ff_0_97, \
                         ff_0_98, fg_1_141, fg_1_142, fg_1_143 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_96[k] = -ab_x[k] * ff_1_96[k]
                      + ff_0_96[k]
                      + fg_1_141[k];

            t_97[k] = -ab_x[k] * ff_1_97[k]
                      + ff_0_97[k]
                      + fg_1_142[k];

            t_98[k] = -ab_x[k] * ff_1_98[k]
                      + ff_0_98[k]
                      + fg_1_143[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, ab_x, ab_y, ff_1_60, ff_1_61, ff_1_62, \
                         ff_1_99, ff_0_99, fg_1_91, fg_1_93, fg_1_94, \
                         fg_1_144 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = -ab_x[k] * ff_1_99[k]
                      + ff_0_99[k]
                      + fg_1_144[k];

            t_100[k] = -ab_y[k] * ff_1_60[k]
                       + fg_1_91[k];

            t_101[k] = -ab_y[k] * ff_1_61[k]
                       + fg_1_93[k];

            t_102[k] = -ab_y[k] * ff_1_62[k]
                       + fg_1_94[k];
        }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, ab_y, ff_1_63, ff_1_64, ff_1_65, \
                         ff_1_66, ff_1_67, fg_1_96, fg_1_97, fg_1_98, fg_1_100, \
                         fg_1_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_103[k] = -ab_y[k] * ff_1_63[k]
                       + fg_1_96[k];

            t_104[k] = -ab_y[k] * ff_1_64[k]
                       + fg_1_97[k];

            t_105[k] = -ab_y[k] * ff_1_65[k]
                       + fg_1_98[k];

            t_106[k] = -ab_y[k] * ff_1_66[k]
                       + fg_1_100[k];

            t_107[k] = -ab_y[k] * ff_1_67[k]
                       + fg_1_101[k];
        }
    }
}

static auto
compute_hrr_geom_010x_gf_out_of_second_piece1(CSimdMatrix &buffer,
                                              const CSimdMatrix &coordinates,
                                              const size_t target, const size_t ff_1,
                                              const size_t fg_1, const size_t ncomps,
                                              const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_108 = buffer.data(target + 108 * ncomps + c);
        auto *t_109 = buffer.data(target + 109 * ncomps + c);
        auto *t_110 = buffer.data(target + 110 * ncomps + c);
        auto *t_111 = buffer.data(target + 111 * ncomps + c);
        auto *t_112 = buffer.data(target + 112 * ncomps + c);
        auto *t_113 = buffer.data(target + 113 * ncomps + c);
        auto *t_114 = buffer.data(target + 114 * ncomps + c);
        auto *t_115 = buffer.data(target + 115 * ncomps + c);
        auto *t_116 = buffer.data(target + 116 * ncomps + c);
        auto *t_117 = buffer.data(target + 117 * ncomps + c);
        auto *t_118 = buffer.data(target + 118 * ncomps + c);
        auto *t_119 = buffer.data(target + 119 * ncomps + c);
        auto *t_120 = buffer.data(target + 120 * ncomps + c);
        auto *t_121 = buffer.data(target + 121 * ncomps + c);
        auto *t_122 = buffer.data(target + 122 * ncomps + c);
        auto *t_123 = buffer.data(target + 123 * ncomps + c);
        auto *t_124 = buffer.data(target + 124 * ncomps + c);
        auto *t_125 = buffer.data(target + 125 * ncomps + c);
        auto *t_126 = buffer.data(target + 126 * ncomps + c);
        auto *t_127 = buffer.data(target + 127 * ncomps + c);
        auto *t_128 = buffer.data(target + 128 * ncomps + c);
        auto *t_129 = buffer.data(target + 129 * ncomps + c);
        auto *t_130 = buffer.data(target + 130 * ncomps + c);
        auto *t_131 = buffer.data(target + 131 * ncomps + c);
        auto *t_132 = buffer.data(target + 132 * ncomps + c);
        auto *t_133 = buffer.data(target + 133 * ncomps + c);
        auto *t_134 = buffer.data(target + 134 * ncomps + c);
        auto *t_135 = buffer.data(target + 135 * ncomps + c);
        auto *t_136 = buffer.data(target + 136 * ncomps + c);
        auto *t_137 = buffer.data(target + 137 * ncomps + c);
        auto *t_138 = buffer.data(target + 138 * ncomps + c);
        auto *t_139 = buffer.data(target + 139 * ncomps + c);
        auto *t_140 = buffer.data(target + 140 * ncomps + c);
        auto *t_141 = buffer.data(target + 141 * ncomps + c);
        auto *t_142 = buffer.data(target + 142 * ncomps + c);
        auto *t_143 = buffer.data(target + 143 * ncomps + c);
        auto *t_144 = buffer.data(target + 144 * ncomps + c);
        auto *t_145 = buffer.data(target + 145 * ncomps + c);
        auto *t_146 = buffer.data(target + 146 * ncomps + c);
        auto *t_147 = buffer.data(target + 147 * ncomps + c);
        auto *t_148 = buffer.data(target + 148 * ncomps + c);
        auto *t_149 = buffer.data(target + 149 * ncomps + c);

        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ff_1_68 = buffer.data(ff_1 + 68 * ncomps + c);
        const auto *ff_1_69 = buffer.data(ff_1 + 69 * ncomps + c);
        const auto *ff_1_70 = buffer.data(ff_1 + 70 * ncomps + c);
        const auto *ff_1_71 = buffer.data(ff_1 + 71 * ncomps + c);
        const auto *ff_1_72 = buffer.data(ff_1 + 72 * ncomps + c);
        const auto *ff_1_73 = buffer.data(ff_1 + 73 * ncomps + c);
        const auto *ff_1_74 = buffer.data(ff_1 + 74 * ncomps + c);
        const auto *ff_1_75 = buffer.data(ff_1 + 75 * ncomps + c);
        const auto *ff_1_76 = buffer.data(ff_1 + 76 * ncomps + c);
        const auto *ff_1_77 = buffer.data(ff_1 + 77 * ncomps + c);
        const auto *ff_1_78 = buffer.data(ff_1 + 78 * ncomps + c);
        const auto *ff_1_79 = buffer.data(ff_1 + 79 * ncomps + c);
        const auto *ff_1_80 = buffer.data(ff_1 + 80 * ncomps + c);
        const auto *ff_1_81 = buffer.data(ff_1 + 81 * ncomps + c);
        const auto *ff_1_82 = buffer.data(ff_1 + 82 * ncomps + c);
        const auto *ff_1_83 = buffer.data(ff_1 + 83 * ncomps + c);
        const auto *ff_1_84 = buffer.data(ff_1 + 84 * ncomps + c);
        const auto *ff_1_85 = buffer.data(ff_1 + 85 * ncomps + c);
        const auto *ff_1_86 = buffer.data(ff_1 + 86 * ncomps + c);
        const auto *ff_1_87 = buffer.data(ff_1 + 87 * ncomps + c);
        const auto *ff_1_88 = buffer.data(ff_1 + 88 * ncomps + c);
        const auto *ff_1_89 = buffer.data(ff_1 + 89 * ncomps + c);
        const auto *ff_1_90 = buffer.data(ff_1 + 90 * ncomps + c);
        const auto *ff_1_91 = buffer.data(ff_1 + 91 * ncomps + c);
        const auto *ff_1_92 = buffer.data(ff_1 + 92 * ncomps + c);
        const auto *ff_1_93 = buffer.data(ff_1 + 93 * ncomps + c);
        const auto *ff_1_94 = buffer.data(ff_1 + 94 * ncomps + c);
        const auto *ff_1_95 = buffer.data(ff_1 + 95 * ncomps + c);
        const auto *ff_1_96 = buffer.data(ff_1 + 96 * ncomps + c);
        const auto *ff_1_97 = buffer.data(ff_1 + 97 * ncomps + c);
        const auto *ff_1_98 = buffer.data(ff_1 + 98 * ncomps + c);
        const auto *ff_1_99 = buffer.data(ff_1 + 99 * ncomps + c);

        const auto *fg_1_102 = buffer.data(fg_1 + 102 * ncomps + c);
        const auto *fg_1_103 = buffer.data(fg_1 + 103 * ncomps + c);
        const auto *fg_1_106 = buffer.data(fg_1 + 106 * ncomps + c);
        const auto *fg_1_108 = buffer.data(fg_1 + 108 * ncomps + c);
        const auto *fg_1_109 = buffer.data(fg_1 + 109 * ncomps + c);
        const auto *fg_1_111 = buffer.data(fg_1 + 111 * ncomps + c);
        const auto *fg_1_112 = buffer.data(fg_1 + 112 * ncomps + c);
        const auto *fg_1_113 = buffer.data(fg_1 + 113 * ncomps + c);
        const auto *fg_1_115 = buffer.data(fg_1 + 115 * ncomps + c);
        const auto *fg_1_116 = buffer.data(fg_1 + 116 * ncomps + c);
        const auto *fg_1_117 = buffer.data(fg_1 + 117 * ncomps + c);
        const auto *fg_1_118 = buffer.data(fg_1 + 118 * ncomps + c);
        const auto *fg_1_121 = buffer.data(fg_1 + 121 * ncomps + c);
        const auto *fg_1_123 = buffer.data(fg_1 + 123 * ncomps + c);
        const auto *fg_1_124 = buffer.data(fg_1 + 124 * ncomps + c);
        const auto *fg_1_126 = buffer.data(fg_1 + 126 * ncomps + c);
        const auto *fg_1_127 = buffer.data(fg_1 + 127 * ncomps + c);
        const auto *fg_1_128 = buffer.data(fg_1 + 128 * ncomps + c);
        const auto *fg_1_130 = buffer.data(fg_1 + 130 * ncomps + c);
        const auto *fg_1_131 = buffer.data(fg_1 + 131 * ncomps + c);
        const auto *fg_1_132 = buffer.data(fg_1 + 132 * ncomps + c);
        const auto *fg_1_133 = buffer.data(fg_1 + 133 * ncomps + c);
        const auto *fg_1_136 = buffer.data(fg_1 + 136 * ncomps + c);
        const auto *fg_1_137 = buffer.data(fg_1 + 137 * ncomps + c);
        const auto *fg_1_138 = buffer.data(fg_1 + 138 * ncomps + c);
        const auto *fg_1_139 = buffer.data(fg_1 + 139 * ncomps + c);
        const auto *fg_1_140 = buffer.data(fg_1 + 140 * ncomps + c);
        const auto *fg_1_141 = buffer.data(fg_1 + 141 * ncomps + c);
        const auto *fg_1_142 = buffer.data(fg_1 + 142 * ncomps + c);
        const auto *fg_1_143 = buffer.data(fg_1 + 143 * ncomps + c);
        const auto *fg_1_144 = buffer.data(fg_1 + 144 * ncomps + c);
        const auto *fg_1_145 = buffer.data(fg_1 + 145 * ncomps + c);
        const auto *fg_1_146 = buffer.data(fg_1 + 146 * ncomps + c);
        const auto *fg_1_147 = buffer.data(fg_1 + 147 * ncomps + c);
        const auto *fg_1_148 = buffer.data(fg_1 + 148 * ncomps + c);
        const auto *fg_1_149 = buffer.data(fg_1 + 149 * ncomps + c);

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, ab_y, ff_1_68, ff_1_69, ff_1_70, \
                         ff_1_71, ff_1_72, fg_1_102, fg_1_103, fg_1_106, fg_1_108, \
                         fg_1_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_108[k] = -ab_y[k] * ff_1_68[k]
                       + fg_1_102[k];

            t_109[k] = -ab_y[k] * ff_1_69[k]
                       + fg_1_103[k];

            t_110[k] = -ab_y[k] * ff_1_70[k]
                       + fg_1_106[k];

            t_111[k] = -ab_y[k] * ff_1_71[k]
                       + fg_1_108[k];

            t_112[k] = -ab_y[k] * ff_1_72[k]
                       + fg_1_109[k];
        }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, ab_y, ff_1_73, ff_1_74, ff_1_75, \
                         ff_1_76, ff_1_77, fg_1_111, fg_1_112, fg_1_113, fg_1_115, \
                         fg_1_116 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_113[k] = -ab_y[k] * ff_1_73[k]
                       + fg_1_111[k];

            t_114[k] = -ab_y[k] * ff_1_74[k]
                       + fg_1_112[k];

            t_115[k] = -ab_y[k] * ff_1_75[k]
                       + fg_1_113[k];

            t_116[k] = -ab_y[k] * ff_1_76[k]
                       + fg_1_115[k];

            t_117[k] = -ab_y[k] * ff_1_77[k]
                       + fg_1_116[k];
        }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, ab_y, ff_1_78, ff_1_79, ff_1_80, \
                         ff_1_81, ff_1_82, fg_1_117, fg_1_118, fg_1_121, fg_1_123, \
                         fg_1_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_118[k] = -ab_y[k] * ff_1_78[k]
                       + fg_1_117[k];

            t_119[k] = -ab_y[k] * ff_1_79[k]
                       + fg_1_118[k];

            t_120[k] = -ab_y[k] * ff_1_80[k]
                       + fg_1_121[k];

            t_121[k] = -ab_y[k] * ff_1_81[k]
                       + fg_1_123[k];

            t_122[k] = -ab_y[k] * ff_1_82[k]
                       + fg_1_124[k];
        }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, t_127, ab_y, ff_1_83, ff_1_84, ff_1_85, \
                         ff_1_86, ff_1_87, fg_1_126, fg_1_127, fg_1_128, fg_1_130, \
                         fg_1_131 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_123[k] = -ab_y[k] * ff_1_83[k]
                       + fg_1_126[k];

            t_124[k] = -ab_y[k] * ff_1_84[k]
                       + fg_1_127[k];

            t_125[k] = -ab_y[k] * ff_1_85[k]
                       + fg_1_128[k];

            t_126[k] = -ab_y[k] * ff_1_86[k]
                       + fg_1_130[k];

            t_127[k] = -ab_y[k] * ff_1_87[k]
                       + fg_1_131[k];
        }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, ab_y, ff_1_88, ff_1_89, ff_1_90, \
                         ff_1_91, ff_1_92, fg_1_132, fg_1_133, fg_1_136, fg_1_138, \
                         fg_1_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_128[k] = -ab_y[k] * ff_1_88[k]
                       + fg_1_132[k];

            t_129[k] = -ab_y[k] * ff_1_89[k]
                       + fg_1_133[k];

            t_130[k] = -ab_y[k] * ff_1_90[k]
                       + fg_1_136[k];

            t_131[k] = -ab_y[k] * ff_1_91[k]
                       + fg_1_138[k];

            t_132[k] = -ab_y[k] * ff_1_92[k]
                       + fg_1_139[k];
        }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, ab_y, ff_1_93, ff_1_94, ff_1_95, \
                         ff_1_96, ff_1_97, fg_1_141, fg_1_142, fg_1_143, fg_1_145, \
                         fg_1_146 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_133[k] = -ab_y[k] * ff_1_93[k]
                       + fg_1_141[k];

            t_134[k] = -ab_y[k] * ff_1_94[k]
                       + fg_1_142[k];

            t_135[k] = -ab_y[k] * ff_1_95[k]
                       + fg_1_143[k];

            t_136[k] = -ab_y[k] * ff_1_96[k]
                       + fg_1_145[k];

            t_137[k] = -ab_y[k] * ff_1_97[k]
                       + fg_1_146[k];
        }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, ab_y, ab_z, ff_1_90, ff_1_91, ff_1_98, \
                         ff_1_99, fg_1_137, fg_1_139, fg_1_147, \
                         fg_1_148 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_138[k] = -ab_y[k] * ff_1_98[k]
                       + fg_1_147[k];

            t_139[k] = -ab_y[k] * ff_1_99[k]
                       + fg_1_148[k];

            t_140[k] = -ab_z[k] * ff_1_90[k]
                       + fg_1_137[k];

            t_141[k] = -ab_z[k] * ff_1_91[k]
                       + fg_1_139[k];
        }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, ab_z, ff_1_92, ff_1_93, ff_1_94, \
                         ff_1_95, ff_1_96, fg_1_140, fg_1_142, fg_1_143, fg_1_144, \
                         fg_1_146 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_142[k] = -ab_z[k] * ff_1_92[k]
                       + fg_1_140[k];

            t_143[k] = -ab_z[k] * ff_1_93[k]
                       + fg_1_142[k];

            t_144[k] = -ab_z[k] * ff_1_94[k]
                       + fg_1_143[k];

            t_145[k] = -ab_z[k] * ff_1_95[k]
                       + fg_1_144[k];

            t_146[k] = -ab_z[k] * ff_1_96[k]
                       + fg_1_146[k];
        }

#pragma omp simd aligned(t_147, t_148, t_149, ab_z, ff_1_97, ff_1_98, ff_1_99, fg_1_147, \
                         fg_1_148, fg_1_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_147[k] = -ab_z[k] * ff_1_97[k]
                       + fg_1_147[k];

            t_148[k] = -ab_z[k] * ff_1_98[k]
                       + fg_1_148[k];

            t_149[k] = -ab_z[k] * ff_1_99[k]
                       + fg_1_149[k];
        }
    }
}

auto
compute_hrr_geom_010x_gf_out_of_second(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                       const size_t target, const size_t ff_1, const size_t ff_0,
                                       const size_t fg_1, const size_t ncomps,
                                       const size_t nmax) -> void
{
    compute_hrr_geom_010x_gf_out_of_second_piece0(buffer, coordinates, target, ff_1, ff_0, fg_1,
                                                  ncomps, nmax);

    compute_hrr_geom_010x_gf_out_of_second_piece1(buffer, coordinates, target, ff_1, fg_1,
                                                  ncomps, nmax);
}

}  // namespace simdtrf
