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


#include "SimdTransferGeom100XHF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_geom_100x_hf_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                             const size_t target, const size_t hd_1,
                                             const size_t hd_0, const size_t id_1,
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
        auto *t_108 = buffer.data(target + 108 * ncomps + c);
        auto *t_109 = buffer.data(target + 109 * ncomps + c);
        auto *t_110 = buffer.data(target + 110 * ncomps + c);
        auto *t_111 = buffer.data(target + 111 * ncomps + c);
        auto *t_112 = buffer.data(target + 112 * ncomps + c);
        auto *t_113 = buffer.data(target + 113 * ncomps + c);
        auto *t_114 = buffer.data(target + 114 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *hd_1_0 = buffer.data(hd_1 + 0 * ncomps + c);
        const auto *hd_1_1 = buffer.data(hd_1 + 1 * ncomps + c);
        const auto *hd_1_2 = buffer.data(hd_1 + 2 * ncomps + c);
        const auto *hd_1_3 = buffer.data(hd_1 + 3 * ncomps + c);
        const auto *hd_1_4 = buffer.data(hd_1 + 4 * ncomps + c);
        const auto *hd_1_5 = buffer.data(hd_1 + 5 * ncomps + c);
        const auto *hd_1_6 = buffer.data(hd_1 + 6 * ncomps + c);
        const auto *hd_1_7 = buffer.data(hd_1 + 7 * ncomps + c);
        const auto *hd_1_8 = buffer.data(hd_1 + 8 * ncomps + c);
        const auto *hd_1_9 = buffer.data(hd_1 + 9 * ncomps + c);
        const auto *hd_1_10 = buffer.data(hd_1 + 10 * ncomps + c);
        const auto *hd_1_11 = buffer.data(hd_1 + 11 * ncomps + c);
        const auto *hd_1_12 = buffer.data(hd_1 + 12 * ncomps + c);
        const auto *hd_1_13 = buffer.data(hd_1 + 13 * ncomps + c);
        const auto *hd_1_14 = buffer.data(hd_1 + 14 * ncomps + c);
        const auto *hd_1_15 = buffer.data(hd_1 + 15 * ncomps + c);
        const auto *hd_1_16 = buffer.data(hd_1 + 16 * ncomps + c);
        const auto *hd_1_17 = buffer.data(hd_1 + 17 * ncomps + c);
        const auto *hd_1_18 = buffer.data(hd_1 + 18 * ncomps + c);
        const auto *hd_1_19 = buffer.data(hd_1 + 19 * ncomps + c);
        const auto *hd_1_20 = buffer.data(hd_1 + 20 * ncomps + c);
        const auto *hd_1_21 = buffer.data(hd_1 + 21 * ncomps + c);
        const auto *hd_1_22 = buffer.data(hd_1 + 22 * ncomps + c);
        const auto *hd_1_23 = buffer.data(hd_1 + 23 * ncomps + c);
        const auto *hd_1_24 = buffer.data(hd_1 + 24 * ncomps + c);
        const auto *hd_1_25 = buffer.data(hd_1 + 25 * ncomps + c);
        const auto *hd_1_26 = buffer.data(hd_1 + 26 * ncomps + c);
        const auto *hd_1_27 = buffer.data(hd_1 + 27 * ncomps + c);
        const auto *hd_1_28 = buffer.data(hd_1 + 28 * ncomps + c);
        const auto *hd_1_29 = buffer.data(hd_1 + 29 * ncomps + c);
        const auto *hd_1_30 = buffer.data(hd_1 + 30 * ncomps + c);
        const auto *hd_1_31 = buffer.data(hd_1 + 31 * ncomps + c);
        const auto *hd_1_32 = buffer.data(hd_1 + 32 * ncomps + c);
        const auto *hd_1_33 = buffer.data(hd_1 + 33 * ncomps + c);
        const auto *hd_1_34 = buffer.data(hd_1 + 34 * ncomps + c);
        const auto *hd_1_35 = buffer.data(hd_1 + 35 * ncomps + c);
        const auto *hd_1_36 = buffer.data(hd_1 + 36 * ncomps + c);
        const auto *hd_1_37 = buffer.data(hd_1 + 37 * ncomps + c);
        const auto *hd_1_38 = buffer.data(hd_1 + 38 * ncomps + c);
        const auto *hd_1_39 = buffer.data(hd_1 + 39 * ncomps + c);
        const auto *hd_1_40 = buffer.data(hd_1 + 40 * ncomps + c);
        const auto *hd_1_41 = buffer.data(hd_1 + 41 * ncomps + c);
        const auto *hd_1_42 = buffer.data(hd_1 + 42 * ncomps + c);
        const auto *hd_1_43 = buffer.data(hd_1 + 43 * ncomps + c);
        const auto *hd_1_44 = buffer.data(hd_1 + 44 * ncomps + c);
        const auto *hd_1_45 = buffer.data(hd_1 + 45 * ncomps + c);
        const auto *hd_1_46 = buffer.data(hd_1 + 46 * ncomps + c);
        const auto *hd_1_47 = buffer.data(hd_1 + 47 * ncomps + c);
        const auto *hd_1_48 = buffer.data(hd_1 + 48 * ncomps + c);
        const auto *hd_1_49 = buffer.data(hd_1 + 49 * ncomps + c);
        const auto *hd_1_50 = buffer.data(hd_1 + 50 * ncomps + c);
        const auto *hd_1_51 = buffer.data(hd_1 + 51 * ncomps + c);
        const auto *hd_1_52 = buffer.data(hd_1 + 52 * ncomps + c);
        const auto *hd_1_53 = buffer.data(hd_1 + 53 * ncomps + c);
        const auto *hd_1_54 = buffer.data(hd_1 + 54 * ncomps + c);
        const auto *hd_1_55 = buffer.data(hd_1 + 55 * ncomps + c);
        const auto *hd_1_56 = buffer.data(hd_1 + 56 * ncomps + c);
        const auto *hd_1_57 = buffer.data(hd_1 + 57 * ncomps + c);
        const auto *hd_1_58 = buffer.data(hd_1 + 58 * ncomps + c);
        const auto *hd_1_59 = buffer.data(hd_1 + 59 * ncomps + c);
        const auto *hd_1_60 = buffer.data(hd_1 + 60 * ncomps + c);
        const auto *hd_1_61 = buffer.data(hd_1 + 61 * ncomps + c);
        const auto *hd_1_62 = buffer.data(hd_1 + 62 * ncomps + c);
        const auto *hd_1_63 = buffer.data(hd_1 + 63 * ncomps + c);
        const auto *hd_1_64 = buffer.data(hd_1 + 64 * ncomps + c);
        const auto *hd_1_65 = buffer.data(hd_1 + 65 * ncomps + c);
        const auto *hd_1_66 = buffer.data(hd_1 + 66 * ncomps + c);
        const auto *hd_1_67 = buffer.data(hd_1 + 67 * ncomps + c);
        const auto *hd_1_68 = buffer.data(hd_1 + 68 * ncomps + c);
        const auto *hd_1_69 = buffer.data(hd_1 + 69 * ncomps + c);
        const auto *hd_1_70 = buffer.data(hd_1 + 70 * ncomps + c);

        const auto *hd_0_0 = buffer.data(hd_0 + 0 * ncomps + c);
        const auto *hd_0_1 = buffer.data(hd_0 + 1 * ncomps + c);
        const auto *hd_0_2 = buffer.data(hd_0 + 2 * ncomps + c);
        const auto *hd_0_3 = buffer.data(hd_0 + 3 * ncomps + c);
        const auto *hd_0_4 = buffer.data(hd_0 + 4 * ncomps + c);
        const auto *hd_0_5 = buffer.data(hd_0 + 5 * ncomps + c);
        const auto *hd_0_6 = buffer.data(hd_0 + 6 * ncomps + c);
        const auto *hd_0_7 = buffer.data(hd_0 + 7 * ncomps + c);
        const auto *hd_0_8 = buffer.data(hd_0 + 8 * ncomps + c);
        const auto *hd_0_9 = buffer.data(hd_0 + 9 * ncomps + c);
        const auto *hd_0_10 = buffer.data(hd_0 + 10 * ncomps + c);
        const auto *hd_0_11 = buffer.data(hd_0 + 11 * ncomps + c);
        const auto *hd_0_12 = buffer.data(hd_0 + 12 * ncomps + c);
        const auto *hd_0_13 = buffer.data(hd_0 + 13 * ncomps + c);
        const auto *hd_0_14 = buffer.data(hd_0 + 14 * ncomps + c);
        const auto *hd_0_15 = buffer.data(hd_0 + 15 * ncomps + c);
        const auto *hd_0_16 = buffer.data(hd_0 + 16 * ncomps + c);
        const auto *hd_0_17 = buffer.data(hd_0 + 17 * ncomps + c);
        const auto *hd_0_18 = buffer.data(hd_0 + 18 * ncomps + c);
        const auto *hd_0_19 = buffer.data(hd_0 + 19 * ncomps + c);
        const auto *hd_0_20 = buffer.data(hd_0 + 20 * ncomps + c);
        const auto *hd_0_21 = buffer.data(hd_0 + 21 * ncomps + c);
        const auto *hd_0_22 = buffer.data(hd_0 + 22 * ncomps + c);
        const auto *hd_0_23 = buffer.data(hd_0 + 23 * ncomps + c);
        const auto *hd_0_24 = buffer.data(hd_0 + 24 * ncomps + c);
        const auto *hd_0_25 = buffer.data(hd_0 + 25 * ncomps + c);
        const auto *hd_0_26 = buffer.data(hd_0 + 26 * ncomps + c);
        const auto *hd_0_27 = buffer.data(hd_0 + 27 * ncomps + c);
        const auto *hd_0_28 = buffer.data(hd_0 + 28 * ncomps + c);
        const auto *hd_0_29 = buffer.data(hd_0 + 29 * ncomps + c);
        const auto *hd_0_30 = buffer.data(hd_0 + 30 * ncomps + c);
        const auto *hd_0_31 = buffer.data(hd_0 + 31 * ncomps + c);
        const auto *hd_0_32 = buffer.data(hd_0 + 32 * ncomps + c);
        const auto *hd_0_33 = buffer.data(hd_0 + 33 * ncomps + c);
        const auto *hd_0_34 = buffer.data(hd_0 + 34 * ncomps + c);
        const auto *hd_0_35 = buffer.data(hd_0 + 35 * ncomps + c);
        const auto *hd_0_36 = buffer.data(hd_0 + 36 * ncomps + c);
        const auto *hd_0_37 = buffer.data(hd_0 + 37 * ncomps + c);
        const auto *hd_0_38 = buffer.data(hd_0 + 38 * ncomps + c);
        const auto *hd_0_39 = buffer.data(hd_0 + 39 * ncomps + c);
        const auto *hd_0_40 = buffer.data(hd_0 + 40 * ncomps + c);
        const auto *hd_0_41 = buffer.data(hd_0 + 41 * ncomps + c);
        const auto *hd_0_42 = buffer.data(hd_0 + 42 * ncomps + c);
        const auto *hd_0_43 = buffer.data(hd_0 + 43 * ncomps + c);
        const auto *hd_0_44 = buffer.data(hd_0 + 44 * ncomps + c);
        const auto *hd_0_45 = buffer.data(hd_0 + 45 * ncomps + c);
        const auto *hd_0_46 = buffer.data(hd_0 + 46 * ncomps + c);
        const auto *hd_0_47 = buffer.data(hd_0 + 47 * ncomps + c);
        const auto *hd_0_48 = buffer.data(hd_0 + 48 * ncomps + c);
        const auto *hd_0_49 = buffer.data(hd_0 + 49 * ncomps + c);
        const auto *hd_0_50 = buffer.data(hd_0 + 50 * ncomps + c);
        const auto *hd_0_51 = buffer.data(hd_0 + 51 * ncomps + c);
        const auto *hd_0_52 = buffer.data(hd_0 + 52 * ncomps + c);
        const auto *hd_0_53 = buffer.data(hd_0 + 53 * ncomps + c);
        const auto *hd_0_54 = buffer.data(hd_0 + 54 * ncomps + c);
        const auto *hd_0_55 = buffer.data(hd_0 + 55 * ncomps + c);
        const auto *hd_0_56 = buffer.data(hd_0 + 56 * ncomps + c);
        const auto *hd_0_57 = buffer.data(hd_0 + 57 * ncomps + c);
        const auto *hd_0_58 = buffer.data(hd_0 + 58 * ncomps + c);
        const auto *hd_0_59 = buffer.data(hd_0 + 59 * ncomps + c);
        const auto *hd_0_60 = buffer.data(hd_0 + 60 * ncomps + c);
        const auto *hd_0_61 = buffer.data(hd_0 + 61 * ncomps + c);
        const auto *hd_0_62 = buffer.data(hd_0 + 62 * ncomps + c);
        const auto *hd_0_63 = buffer.data(hd_0 + 63 * ncomps + c);
        const auto *hd_0_64 = buffer.data(hd_0 + 64 * ncomps + c);
        const auto *hd_0_65 = buffer.data(hd_0 + 65 * ncomps + c);
        const auto *hd_0_66 = buffer.data(hd_0 + 66 * ncomps + c);
        const auto *hd_0_67 = buffer.data(hd_0 + 67 * ncomps + c);
        const auto *hd_0_68 = buffer.data(hd_0 + 68 * ncomps + c);
        const auto *hd_0_69 = buffer.data(hd_0 + 69 * ncomps + c);
        const auto *hd_0_70 = buffer.data(hd_0 + 70 * ncomps + c);

        const auto *id_1_0 = buffer.data(id_1 + 0 * ncomps + c);
        const auto *id_1_1 = buffer.data(id_1 + 1 * ncomps + c);
        const auto *id_1_2 = buffer.data(id_1 + 2 * ncomps + c);
        const auto *id_1_3 = buffer.data(id_1 + 3 * ncomps + c);
        const auto *id_1_4 = buffer.data(id_1 + 4 * ncomps + c);
        const auto *id_1_5 = buffer.data(id_1 + 5 * ncomps + c);
        const auto *id_1_6 = buffer.data(id_1 + 6 * ncomps + c);
        const auto *id_1_7 = buffer.data(id_1 + 7 * ncomps + c);
        const auto *id_1_8 = buffer.data(id_1 + 8 * ncomps + c);
        const auto *id_1_9 = buffer.data(id_1 + 9 * ncomps + c);
        const auto *id_1_10 = buffer.data(id_1 + 10 * ncomps + c);
        const auto *id_1_11 = buffer.data(id_1 + 11 * ncomps + c);
        const auto *id_1_12 = buffer.data(id_1 + 12 * ncomps + c);
        const auto *id_1_13 = buffer.data(id_1 + 13 * ncomps + c);
        const auto *id_1_14 = buffer.data(id_1 + 14 * ncomps + c);
        const auto *id_1_15 = buffer.data(id_1 + 15 * ncomps + c);
        const auto *id_1_16 = buffer.data(id_1 + 16 * ncomps + c);
        const auto *id_1_17 = buffer.data(id_1 + 17 * ncomps + c);
        const auto *id_1_18 = buffer.data(id_1 + 18 * ncomps + c);
        const auto *id_1_19 = buffer.data(id_1 + 19 * ncomps + c);
        const auto *id_1_20 = buffer.data(id_1 + 20 * ncomps + c);
        const auto *id_1_21 = buffer.data(id_1 + 21 * ncomps + c);
        const auto *id_1_22 = buffer.data(id_1 + 22 * ncomps + c);
        const auto *id_1_23 = buffer.data(id_1 + 23 * ncomps + c);
        const auto *id_1_24 = buffer.data(id_1 + 24 * ncomps + c);
        const auto *id_1_25 = buffer.data(id_1 + 25 * ncomps + c);
        const auto *id_1_26 = buffer.data(id_1 + 26 * ncomps + c);
        const auto *id_1_27 = buffer.data(id_1 + 27 * ncomps + c);
        const auto *id_1_28 = buffer.data(id_1 + 28 * ncomps + c);
        const auto *id_1_29 = buffer.data(id_1 + 29 * ncomps + c);
        const auto *id_1_30 = buffer.data(id_1 + 30 * ncomps + c);
        const auto *id_1_31 = buffer.data(id_1 + 31 * ncomps + c);
        const auto *id_1_32 = buffer.data(id_1 + 32 * ncomps + c);
        const auto *id_1_33 = buffer.data(id_1 + 33 * ncomps + c);
        const auto *id_1_34 = buffer.data(id_1 + 34 * ncomps + c);
        const auto *id_1_35 = buffer.data(id_1 + 35 * ncomps + c);
        const auto *id_1_36 = buffer.data(id_1 + 36 * ncomps + c);
        const auto *id_1_37 = buffer.data(id_1 + 37 * ncomps + c);
        const auto *id_1_38 = buffer.data(id_1 + 38 * ncomps + c);
        const auto *id_1_39 = buffer.data(id_1 + 39 * ncomps + c);
        const auto *id_1_40 = buffer.data(id_1 + 40 * ncomps + c);
        const auto *id_1_41 = buffer.data(id_1 + 41 * ncomps + c);
        const auto *id_1_42 = buffer.data(id_1 + 42 * ncomps + c);
        const auto *id_1_43 = buffer.data(id_1 + 43 * ncomps + c);
        const auto *id_1_44 = buffer.data(id_1 + 44 * ncomps + c);
        const auto *id_1_45 = buffer.data(id_1 + 45 * ncomps + c);
        const auto *id_1_46 = buffer.data(id_1 + 46 * ncomps + c);
        const auto *id_1_47 = buffer.data(id_1 + 47 * ncomps + c);
        const auto *id_1_48 = buffer.data(id_1 + 48 * ncomps + c);
        const auto *id_1_49 = buffer.data(id_1 + 49 * ncomps + c);
        const auto *id_1_50 = buffer.data(id_1 + 50 * ncomps + c);
        const auto *id_1_51 = buffer.data(id_1 + 51 * ncomps + c);
        const auto *id_1_52 = buffer.data(id_1 + 52 * ncomps + c);
        const auto *id_1_53 = buffer.data(id_1 + 53 * ncomps + c);
        const auto *id_1_54 = buffer.data(id_1 + 54 * ncomps + c);
        const auto *id_1_55 = buffer.data(id_1 + 55 * ncomps + c);
        const auto *id_1_56 = buffer.data(id_1 + 56 * ncomps + c);
        const auto *id_1_57 = buffer.data(id_1 + 57 * ncomps + c);
        const auto *id_1_58 = buffer.data(id_1 + 58 * ncomps + c);
        const auto *id_1_59 = buffer.data(id_1 + 59 * ncomps + c);
        const auto *id_1_60 = buffer.data(id_1 + 60 * ncomps + c);
        const auto *id_1_61 = buffer.data(id_1 + 61 * ncomps + c);
        const auto *id_1_62 = buffer.data(id_1 + 62 * ncomps + c);
        const auto *id_1_63 = buffer.data(id_1 + 63 * ncomps + c);
        const auto *id_1_64 = buffer.data(id_1 + 64 * ncomps + c);
        const auto *id_1_65 = buffer.data(id_1 + 65 * ncomps + c);
        const auto *id_1_66 = buffer.data(id_1 + 66 * ncomps + c);
        const auto *id_1_67 = buffer.data(id_1 + 67 * ncomps + c);
        const auto *id_1_68 = buffer.data(id_1 + 68 * ncomps + c);
        const auto *id_1_69 = buffer.data(id_1 + 69 * ncomps + c);
        const auto *id_1_70 = buffer.data(id_1 + 70 * ncomps + c);
        const auto *id_1_71 = buffer.data(id_1 + 71 * ncomps + c);
        const auto *id_1_75 = buffer.data(id_1 + 75 * ncomps + c);
        const auto *id_1_76 = buffer.data(id_1 + 76 * ncomps + c);
        const auto *id_1_77 = buffer.data(id_1 + 77 * ncomps + c);
        const auto *id_1_81 = buffer.data(id_1 + 81 * ncomps + c);
        const auto *id_1_82 = buffer.data(id_1 + 82 * ncomps + c);
        const auto *id_1_83 = buffer.data(id_1 + 83 * ncomps + c);
        const auto *id_1_89 = buffer.data(id_1 + 89 * ncomps + c);
        const auto *id_1_93 = buffer.data(id_1 + 93 * ncomps + c);
        const auto *id_1_94 = buffer.data(id_1 + 94 * ncomps + c);
        const auto *id_1_95 = buffer.data(id_1 + 95 * ncomps + c);
        const auto *id_1_101 = buffer.data(id_1 + 101 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, ab_x, hd_1_0, hd_1_1, hd_1_2, hd_0_0, hd_0_1, hd_0_2, \
                         id_1_0, id_1_1, id_1_2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * hd_1_0[k]
                     + hd_0_0[k]
                     + id_1_0[k];

            t_1[k] = ab_x[k] * hd_1_1[k]
                     + hd_0_1[k]
                     + id_1_1[k];

            t_2[k] = ab_x[k] * hd_1_2[k]
                     + hd_0_2[k]
                     + id_1_2[k];
        }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, ab_x, ab_y, hd_1_3, hd_1_4, hd_1_5, hd_0_3, \
                         hd_0_4, hd_0_5, id_1_3, id_1_4, id_1_5, \
                         id_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_3[k] = ab_x[k] * hd_1_3[k]
                     + hd_0_3[k]
                     + id_1_3[k];

            t_4[k] = ab_x[k] * hd_1_4[k]
                     + hd_0_4[k]
                     + id_1_4[k];

            t_5[k] = ab_x[k] * hd_1_5[k]
                     + hd_0_5[k]
                     + id_1_5[k];

            t_6[k] = ab_y[k] * hd_1_3[k]
                     + id_1_9[k];
        }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, hd_1_4, hd_1_5, hd_1_6, \
                         hd_0_6, id_1_6, id_1_10, id_1_11, id_1_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_7[k] = ab_y[k] * hd_1_4[k]
                     + id_1_10[k];

            t_8[k] = ab_y[k] * hd_1_5[k]
                     + id_1_11[k];

            t_9[k] = ab_z[k] * hd_1_5[k]
                     + id_1_17[k];

            t_10[k] = ab_x[k] * hd_1_6[k]
                      + hd_0_6[k]
                      + id_1_6[k];
        }

#pragma omp simd aligned(t_11, t_12, t_13, ab_x, hd_1_7, hd_1_8, hd_1_9, hd_0_7, hd_0_8, \
                         hd_0_9, id_1_7, id_1_8, id_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_11[k] = ab_x[k] * hd_1_7[k]
                      + hd_0_7[k]
                      + id_1_7[k];

            t_12[k] = ab_x[k] * hd_1_8[k]
                      + hd_0_8[k]
                      + id_1_8[k];

            t_13[k] = ab_x[k] * hd_1_9[k]
                      + hd_0_9[k]
                      + id_1_9[k];
        }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, ab_x, ab_y, hd_1_9, hd_1_10, hd_1_11, \
                         hd_0_10, hd_0_11, id_1_10, id_1_11, id_1_21, \
                         id_1_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_14[k] = ab_x[k] * hd_1_10[k]
                      + hd_0_10[k]
                      + id_1_10[k];

            t_15[k] = ab_x[k] * hd_1_11[k]
                      + hd_0_11[k]
                      + id_1_11[k];

            t_16[k] = ab_y[k] * hd_1_9[k]
                      + id_1_21[k];

            t_17[k] = ab_y[k] * hd_1_10[k]
                      + id_1_22[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, hd_1_11, hd_1_12, hd_1_13, \
                         hd_0_12, hd_0_13, id_1_12, id_1_13, id_1_23, \
                         id_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = ab_y[k] * hd_1_11[k]
                      + id_1_23[k];

            t_19[k] = ab_z[k] * hd_1_11[k]
                      + id_1_29[k];

            t_20[k] = ab_x[k] * hd_1_12[k]
                      + hd_0_12[k]
                      + id_1_12[k];

            t_21[k] = ab_x[k] * hd_1_13[k]
                      + hd_0_13[k]
                      + id_1_13[k];
        }

#pragma omp simd aligned(t_22, t_23, t_24, ab_x, hd_1_14, hd_1_15, hd_1_16, hd_0_14, hd_0_15, \
                         hd_0_16, id_1_14, id_1_15, id_1_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_22[k] = ab_x[k] * hd_1_14[k]
                      + hd_0_14[k]
                      + id_1_14[k];

            t_23[k] = ab_x[k] * hd_1_15[k]
                      + hd_0_15[k]
                      + id_1_15[k];

            t_24[k] = ab_x[k] * hd_1_16[k]
                      + hd_0_16[k]
                      + id_1_16[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, ab_x, ab_y, hd_1_15, hd_1_16, hd_1_17, \
                         hd_0_17, id_1_17, id_1_27, id_1_28, id_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_x[k] * hd_1_17[k]
                      + hd_0_17[k]
                      + id_1_17[k];

            t_26[k] = ab_y[k] * hd_1_15[k]
                      + id_1_27[k];

            t_27[k] = ab_y[k] * hd_1_16[k]
                      + id_1_28[k];

            t_28[k] = ab_y[k] * hd_1_17[k]
                      + id_1_29[k];
        }

#pragma omp simd aligned(t_29, t_30, t_31, ab_x, ab_z, hd_1_17, hd_1_18, hd_1_19, hd_0_18, \
                         hd_0_19, id_1_18, id_1_19, id_1_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_29[k] = ab_z[k] * hd_1_17[k]
                      + id_1_35[k];

            t_30[k] = ab_x[k] * hd_1_18[k]
                      + hd_0_18[k]
                      + id_1_18[k];

            t_31[k] = ab_x[k] * hd_1_19[k]
                      + hd_0_19[k]
                      + id_1_19[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, ab_x, hd_1_20, hd_1_21, hd_1_22, hd_0_20, hd_0_21, \
                         hd_0_22, id_1_20, id_1_21, id_1_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_x[k] * hd_1_20[k]
                      + hd_0_20[k]
                      + id_1_20[k];

            t_33[k] = ab_x[k] * hd_1_21[k]
                      + hd_0_21[k]
                      + id_1_21[k];

            t_34[k] = ab_x[k] * hd_1_22[k]
                      + hd_0_22[k]
                      + id_1_22[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, ab_x, ab_y, hd_1_21, hd_1_22, hd_1_23, \
                         hd_0_23, id_1_23, id_1_39, id_1_40, id_1_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * hd_1_23[k]
                      + hd_0_23[k]
                      + id_1_23[k];

            t_36[k] = ab_y[k] * hd_1_21[k]
                      + id_1_39[k];

            t_37[k] = ab_y[k] * hd_1_22[k]
                      + id_1_40[k];

            t_38[k] = ab_y[k] * hd_1_23[k]
                      + id_1_41[k];
        }

#pragma omp simd aligned(t_39, t_40, t_41, ab_x, ab_z, hd_1_23, hd_1_24, hd_1_25, hd_0_24, \
                         hd_0_25, id_1_24, id_1_25, id_1_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_39[k] = ab_z[k] * hd_1_23[k]
                      + id_1_47[k];

            t_40[k] = ab_x[k] * hd_1_24[k]
                      + hd_0_24[k]
                      + id_1_24[k];

            t_41[k] = ab_x[k] * hd_1_25[k]
                      + hd_0_25[k]
                      + id_1_25[k];
        }

#pragma omp simd aligned(t_42, t_43, t_44, ab_x, hd_1_26, hd_1_27, hd_1_28, hd_0_26, hd_0_27, \
                         hd_0_28, id_1_26, id_1_27, id_1_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_42[k] = ab_x[k] * hd_1_26[k]
                      + hd_0_26[k]
                      + id_1_26[k];

            t_43[k] = ab_x[k] * hd_1_27[k]
                      + hd_0_27[k]
                      + id_1_27[k];

            t_44[k] = ab_x[k] * hd_1_28[k]
                      + hd_0_28[k]
                      + id_1_28[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, hd_1_27, hd_1_28, hd_1_29, \
                         hd_0_29, id_1_29, id_1_45, id_1_46, id_1_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * hd_1_29[k]
                      + hd_0_29[k]
                      + id_1_29[k];

            t_46[k] = ab_y[k] * hd_1_27[k]
                      + id_1_45[k];

            t_47[k] = ab_y[k] * hd_1_28[k]
                      + id_1_46[k];

            t_48[k] = ab_y[k] * hd_1_29[k]
                      + id_1_47[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, ab_x, ab_z, hd_1_29, hd_1_30, hd_1_31, hd_0_30, \
                         hd_0_31, id_1_30, id_1_31, id_1_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_z[k] * hd_1_29[k]
                      + id_1_53[k];

            t_50[k] = ab_x[k] * hd_1_30[k]
                      + hd_0_30[k]
                      + id_1_30[k];

            t_51[k] = ab_x[k] * hd_1_31[k]
                      + hd_0_31[k]
                      + id_1_31[k];
        }

#pragma omp simd aligned(t_52, t_53, t_54, ab_x, hd_1_32, hd_1_33, hd_1_34, hd_0_32, hd_0_33, \
                         hd_0_34, id_1_32, id_1_33, id_1_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_52[k] = ab_x[k] * hd_1_32[k]
                      + hd_0_32[k]
                      + id_1_32[k];

            t_53[k] = ab_x[k] * hd_1_33[k]
                      + hd_0_33[k]
                      + id_1_33[k];

            t_54[k] = ab_x[k] * hd_1_34[k]
                      + hd_0_34[k]
                      + id_1_34[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, ab_x, ab_y, hd_1_33, hd_1_34, hd_1_35, \
                         hd_0_35, id_1_35, id_1_51, id_1_52, id_1_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = ab_x[k] * hd_1_35[k]
                      + hd_0_35[k]
                      + id_1_35[k];

            t_56[k] = ab_y[k] * hd_1_33[k]
                      + id_1_51[k];

            t_57[k] = ab_y[k] * hd_1_34[k]
                      + id_1_52[k];

            t_58[k] = ab_y[k] * hd_1_35[k]
                      + id_1_53[k];
        }

#pragma omp simd aligned(t_59, t_60, t_61, ab_x, ab_z, hd_1_35, hd_1_36, hd_1_37, hd_0_36, \
                         hd_0_37, id_1_36, id_1_37, id_1_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_59[k] = ab_z[k] * hd_1_35[k]
                      + id_1_59[k];

            t_60[k] = ab_x[k] * hd_1_36[k]
                      + hd_0_36[k]
                      + id_1_36[k];

            t_61[k] = ab_x[k] * hd_1_37[k]
                      + hd_0_37[k]
                      + id_1_37[k];
        }

#pragma omp simd aligned(t_62, t_63, t_64, ab_x, hd_1_38, hd_1_39, hd_1_40, hd_0_38, hd_0_39, \
                         hd_0_40, id_1_38, id_1_39, id_1_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_62[k] = ab_x[k] * hd_1_38[k]
                      + hd_0_38[k]
                      + id_1_38[k];

            t_63[k] = ab_x[k] * hd_1_39[k]
                      + hd_0_39[k]
                      + id_1_39[k];

            t_64[k] = ab_x[k] * hd_1_40[k]
                      + hd_0_40[k]
                      + id_1_40[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, ab_x, ab_y, hd_1_39, hd_1_40, hd_1_41, \
                         hd_0_41, id_1_41, id_1_63, id_1_64, id_1_65 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * hd_1_41[k]
                      + hd_0_41[k]
                      + id_1_41[k];

            t_66[k] = ab_y[k] * hd_1_39[k]
                      + id_1_63[k];

            t_67[k] = ab_y[k] * hd_1_40[k]
                      + id_1_64[k];

            t_68[k] = ab_y[k] * hd_1_41[k]
                      + id_1_65[k];
        }

#pragma omp simd aligned(t_69, t_70, t_71, ab_x, ab_z, hd_1_41, hd_1_42, hd_1_43, hd_0_42, \
                         hd_0_43, id_1_42, id_1_43, id_1_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_69[k] = ab_z[k] * hd_1_41[k]
                      + id_1_71[k];

            t_70[k] = ab_x[k] * hd_1_42[k]
                      + hd_0_42[k]
                      + id_1_42[k];

            t_71[k] = ab_x[k] * hd_1_43[k]
                      + hd_0_43[k]
                      + id_1_43[k];
        }

#pragma omp simd aligned(t_72, t_73, t_74, ab_x, hd_1_44, hd_1_45, hd_1_46, hd_0_44, hd_0_45, \
                         hd_0_46, id_1_44, id_1_45, id_1_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_72[k] = ab_x[k] * hd_1_44[k]
                      + hd_0_44[k]
                      + id_1_44[k];

            t_73[k] = ab_x[k] * hd_1_45[k]
                      + hd_0_45[k]
                      + id_1_45[k];

            t_74[k] = ab_x[k] * hd_1_46[k]
                      + hd_0_46[k]
                      + id_1_46[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, ab_x, ab_y, hd_1_45, hd_1_46, hd_1_47, \
                         hd_0_47, id_1_47, id_1_69, id_1_70, id_1_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = ab_x[k] * hd_1_47[k]
                      + hd_0_47[k]
                      + id_1_47[k];

            t_76[k] = ab_y[k] * hd_1_45[k]
                      + id_1_69[k];

            t_77[k] = ab_y[k] * hd_1_46[k]
                      + id_1_70[k];

            t_78[k] = ab_y[k] * hd_1_47[k]
                      + id_1_71[k];
        }

#pragma omp simd aligned(t_79, t_80, t_81, ab_x, ab_z, hd_1_47, hd_1_48, hd_1_49, hd_0_48, \
                         hd_0_49, id_1_48, id_1_49, id_1_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_79[k] = ab_z[k] * hd_1_47[k]
                      + id_1_77[k];

            t_80[k] = ab_x[k] * hd_1_48[k]
                      + hd_0_48[k]
                      + id_1_48[k];

            t_81[k] = ab_x[k] * hd_1_49[k]
                      + hd_0_49[k]
                      + id_1_49[k];
        }

#pragma omp simd aligned(t_82, t_83, t_84, ab_x, hd_1_50, hd_1_51, hd_1_52, hd_0_50, hd_0_51, \
                         hd_0_52, id_1_50, id_1_51, id_1_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_82[k] = ab_x[k] * hd_1_50[k]
                      + hd_0_50[k]
                      + id_1_50[k];

            t_83[k] = ab_x[k] * hd_1_51[k]
                      + hd_0_51[k]
                      + id_1_51[k];

            t_84[k] = ab_x[k] * hd_1_52[k]
                      + hd_0_52[k]
                      + id_1_52[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, ab_x, ab_y, hd_1_51, hd_1_52, hd_1_53, \
                         hd_0_53, id_1_53, id_1_75, id_1_76, id_1_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_x[k] * hd_1_53[k]
                      + hd_0_53[k]
                      + id_1_53[k];

            t_86[k] = ab_y[k] * hd_1_51[k]
                      + id_1_75[k];

            t_87[k] = ab_y[k] * hd_1_52[k]
                      + id_1_76[k];

            t_88[k] = ab_y[k] * hd_1_53[k]
                      + id_1_77[k];
        }

#pragma omp simd aligned(t_89, t_90, t_91, ab_x, ab_z, hd_1_53, hd_1_54, hd_1_55, hd_0_54, \
                         hd_0_55, id_1_54, id_1_55, id_1_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_89[k] = ab_z[k] * hd_1_53[k]
                      + id_1_83[k];

            t_90[k] = ab_x[k] * hd_1_54[k]
                      + hd_0_54[k]
                      + id_1_54[k];

            t_91[k] = ab_x[k] * hd_1_55[k]
                      + hd_0_55[k]
                      + id_1_55[k];
        }

#pragma omp simd aligned(t_92, t_93, t_94, ab_x, hd_1_56, hd_1_57, hd_1_58, hd_0_56, hd_0_57, \
                         hd_0_58, id_1_56, id_1_57, id_1_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_92[k] = ab_x[k] * hd_1_56[k]
                      + hd_0_56[k]
                      + id_1_56[k];

            t_93[k] = ab_x[k] * hd_1_57[k]
                      + hd_0_57[k]
                      + id_1_57[k];

            t_94[k] = ab_x[k] * hd_1_58[k]
                      + hd_0_58[k]
                      + id_1_58[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, ab_x, ab_y, hd_1_57, hd_1_58, hd_1_59, \
                         hd_0_59, id_1_59, id_1_81, id_1_82, id_1_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_x[k] * hd_1_59[k]
                      + hd_0_59[k]
                      + id_1_59[k];

            t_96[k] = ab_y[k] * hd_1_57[k]
                      + id_1_81[k];

            t_97[k] = ab_y[k] * hd_1_58[k]
                      + id_1_82[k];

            t_98[k] = ab_y[k] * hd_1_59[k]
                      + id_1_83[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, ab_x, ab_z, hd_1_59, hd_1_60, hd_1_61, hd_0_60, \
                         hd_0_61, id_1_60, id_1_61, id_1_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_z[k] * hd_1_59[k]
                      + id_1_89[k];

            t_100[k] = ab_x[k] * hd_1_60[k]
                       + hd_0_60[k]
                       + id_1_60[k];

            t_101[k] = ab_x[k] * hd_1_61[k]
                       + hd_0_61[k]
                       + id_1_61[k];
        }

#pragma omp simd aligned(t_102, t_103, t_104, ab_x, hd_1_62, hd_1_63, hd_1_64, hd_0_62, \
                         hd_0_63, hd_0_64, id_1_62, id_1_63, id_1_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_102[k] = ab_x[k] * hd_1_62[k]
                       + hd_0_62[k]
                       + id_1_62[k];

            t_103[k] = ab_x[k] * hd_1_63[k]
                       + hd_0_63[k]
                       + id_1_63[k];

            t_104[k] = ab_x[k] * hd_1_64[k]
                       + hd_0_64[k]
                       + id_1_64[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, ab_x, ab_y, hd_1_63, hd_1_64, hd_1_65, \
                         hd_0_65, id_1_65, id_1_93, id_1_94, id_1_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * hd_1_65[k]
                       + hd_0_65[k]
                       + id_1_65[k];

            t_106[k] = ab_y[k] * hd_1_63[k]
                       + id_1_93[k];

            t_107[k] = ab_y[k] * hd_1_64[k]
                       + id_1_94[k];

            t_108[k] = ab_y[k] * hd_1_65[k]
                       + id_1_95[k];
        }

#pragma omp simd aligned(t_109, t_110, t_111, ab_x, ab_z, hd_1_65, hd_1_66, hd_1_67, hd_0_66, \
                         hd_0_67, id_1_66, id_1_67, id_1_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_109[k] = ab_z[k] * hd_1_65[k]
                       + id_1_101[k];

            t_110[k] = ab_x[k] * hd_1_66[k]
                       + hd_0_66[k]
                       + id_1_66[k];

            t_111[k] = ab_x[k] * hd_1_67[k]
                       + hd_0_67[k]
                       + id_1_67[k];
        }

#pragma omp simd aligned(t_112, t_113, t_114, ab_x, hd_1_68, hd_1_69, hd_1_70, hd_0_68, \
                         hd_0_69, hd_0_70, id_1_68, id_1_69, id_1_70 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_112[k] = ab_x[k] * hd_1_68[k]
                       + hd_0_68[k]
                       + id_1_68[k];

            t_113[k] = ab_x[k] * hd_1_69[k]
                       + hd_0_69[k]
                       + id_1_69[k];

            t_114[k] = ab_x[k] * hd_1_70[k]
                       + hd_0_70[k]
                       + id_1_70[k];
        }
    }
}

static auto
compute_hrr_geom_100x_hf_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                             const size_t target, const size_t hd_1,
                                             const size_t hd_0, const size_t id_1,
                                             const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
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
        auto *t_150 = buffer.data(target + 150 * ncomps + c);
        auto *t_151 = buffer.data(target + 151 * ncomps + c);
        auto *t_152 = buffer.data(target + 152 * ncomps + c);
        auto *t_153 = buffer.data(target + 153 * ncomps + c);
        auto *t_154 = buffer.data(target + 154 * ncomps + c);
        auto *t_155 = buffer.data(target + 155 * ncomps + c);
        auto *t_156 = buffer.data(target + 156 * ncomps + c);
        auto *t_157 = buffer.data(target + 157 * ncomps + c);
        auto *t_158 = buffer.data(target + 158 * ncomps + c);
        auto *t_159 = buffer.data(target + 159 * ncomps + c);
        auto *t_160 = buffer.data(target + 160 * ncomps + c);
        auto *t_161 = buffer.data(target + 161 * ncomps + c);
        auto *t_162 = buffer.data(target + 162 * ncomps + c);
        auto *t_163 = buffer.data(target + 163 * ncomps + c);
        auto *t_164 = buffer.data(target + 164 * ncomps + c);
        auto *t_165 = buffer.data(target + 165 * ncomps + c);
        auto *t_166 = buffer.data(target + 166 * ncomps + c);
        auto *t_167 = buffer.data(target + 167 * ncomps + c);
        auto *t_168 = buffer.data(target + 168 * ncomps + c);
        auto *t_169 = buffer.data(target + 169 * ncomps + c);
        auto *t_170 = buffer.data(target + 170 * ncomps + c);
        auto *t_171 = buffer.data(target + 171 * ncomps + c);
        auto *t_172 = buffer.data(target + 172 * ncomps + c);
        auto *t_173 = buffer.data(target + 173 * ncomps + c);
        auto *t_174 = buffer.data(target + 174 * ncomps + c);
        auto *t_175 = buffer.data(target + 175 * ncomps + c);
        auto *t_176 = buffer.data(target + 176 * ncomps + c);
        auto *t_177 = buffer.data(target + 177 * ncomps + c);
        auto *t_178 = buffer.data(target + 178 * ncomps + c);
        auto *t_179 = buffer.data(target + 179 * ncomps + c);
        auto *t_180 = buffer.data(target + 180 * ncomps + c);
        auto *t_181 = buffer.data(target + 181 * ncomps + c);
        auto *t_182 = buffer.data(target + 182 * ncomps + c);
        auto *t_183 = buffer.data(target + 183 * ncomps + c);
        auto *t_184 = buffer.data(target + 184 * ncomps + c);
        auto *t_185 = buffer.data(target + 185 * ncomps + c);
        auto *t_186 = buffer.data(target + 186 * ncomps + c);
        auto *t_187 = buffer.data(target + 187 * ncomps + c);
        auto *t_188 = buffer.data(target + 188 * ncomps + c);
        auto *t_189 = buffer.data(target + 189 * ncomps + c);
        auto *t_190 = buffer.data(target + 190 * ncomps + c);
        auto *t_191 = buffer.data(target + 191 * ncomps + c);
        auto *t_192 = buffer.data(target + 192 * ncomps + c);
        auto *t_193 = buffer.data(target + 193 * ncomps + c);
        auto *t_194 = buffer.data(target + 194 * ncomps + c);
        auto *t_195 = buffer.data(target + 195 * ncomps + c);
        auto *t_196 = buffer.data(target + 196 * ncomps + c);
        auto *t_197 = buffer.data(target + 197 * ncomps + c);
        auto *t_198 = buffer.data(target + 198 * ncomps + c);
        auto *t_199 = buffer.data(target + 199 * ncomps + c);
        auto *t_200 = buffer.data(target + 200 * ncomps + c);
        auto *t_201 = buffer.data(target + 201 * ncomps + c);
        auto *t_202 = buffer.data(target + 202 * ncomps + c);
        auto *t_203 = buffer.data(target + 203 * ncomps + c);
        auto *t_204 = buffer.data(target + 204 * ncomps + c);
        auto *t_205 = buffer.data(target + 205 * ncomps + c);
        auto *t_206 = buffer.data(target + 206 * ncomps + c);
        auto *t_207 = buffer.data(target + 207 * ncomps + c);
        auto *t_208 = buffer.data(target + 208 * ncomps + c);
        auto *t_209 = buffer.data(target + 209 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *hd_1_69 = buffer.data(hd_1 + 69 * ncomps + c);
        const auto *hd_1_70 = buffer.data(hd_1 + 70 * ncomps + c);
        const auto *hd_1_71 = buffer.data(hd_1 + 71 * ncomps + c);
        const auto *hd_1_72 = buffer.data(hd_1 + 72 * ncomps + c);
        const auto *hd_1_73 = buffer.data(hd_1 + 73 * ncomps + c);
        const auto *hd_1_74 = buffer.data(hd_1 + 74 * ncomps + c);
        const auto *hd_1_75 = buffer.data(hd_1 + 75 * ncomps + c);
        const auto *hd_1_76 = buffer.data(hd_1 + 76 * ncomps + c);
        const auto *hd_1_77 = buffer.data(hd_1 + 77 * ncomps + c);
        const auto *hd_1_78 = buffer.data(hd_1 + 78 * ncomps + c);
        const auto *hd_1_79 = buffer.data(hd_1 + 79 * ncomps + c);
        const auto *hd_1_80 = buffer.data(hd_1 + 80 * ncomps + c);
        const auto *hd_1_81 = buffer.data(hd_1 + 81 * ncomps + c);
        const auto *hd_1_82 = buffer.data(hd_1 + 82 * ncomps + c);
        const auto *hd_1_83 = buffer.data(hd_1 + 83 * ncomps + c);
        const auto *hd_1_84 = buffer.data(hd_1 + 84 * ncomps + c);
        const auto *hd_1_85 = buffer.data(hd_1 + 85 * ncomps + c);
        const auto *hd_1_86 = buffer.data(hd_1 + 86 * ncomps + c);
        const auto *hd_1_87 = buffer.data(hd_1 + 87 * ncomps + c);
        const auto *hd_1_88 = buffer.data(hd_1 + 88 * ncomps + c);
        const auto *hd_1_89 = buffer.data(hd_1 + 89 * ncomps + c);
        const auto *hd_1_90 = buffer.data(hd_1 + 90 * ncomps + c);
        const auto *hd_1_91 = buffer.data(hd_1 + 91 * ncomps + c);
        const auto *hd_1_92 = buffer.data(hd_1 + 92 * ncomps + c);
        const auto *hd_1_93 = buffer.data(hd_1 + 93 * ncomps + c);
        const auto *hd_1_94 = buffer.data(hd_1 + 94 * ncomps + c);
        const auto *hd_1_95 = buffer.data(hd_1 + 95 * ncomps + c);
        const auto *hd_1_96 = buffer.data(hd_1 + 96 * ncomps + c);
        const auto *hd_1_97 = buffer.data(hd_1 + 97 * ncomps + c);
        const auto *hd_1_98 = buffer.data(hd_1 + 98 * ncomps + c);
        const auto *hd_1_99 = buffer.data(hd_1 + 99 * ncomps + c);
        const auto *hd_1_100 = buffer.data(hd_1 + 100 * ncomps + c);
        const auto *hd_1_101 = buffer.data(hd_1 + 101 * ncomps + c);
        const auto *hd_1_102 = buffer.data(hd_1 + 102 * ncomps + c);
        const auto *hd_1_103 = buffer.data(hd_1 + 103 * ncomps + c);
        const auto *hd_1_104 = buffer.data(hd_1 + 104 * ncomps + c);
        const auto *hd_1_105 = buffer.data(hd_1 + 105 * ncomps + c);
        const auto *hd_1_106 = buffer.data(hd_1 + 106 * ncomps + c);
        const auto *hd_1_107 = buffer.data(hd_1 + 107 * ncomps + c);
        const auto *hd_1_108 = buffer.data(hd_1 + 108 * ncomps + c);
        const auto *hd_1_109 = buffer.data(hd_1 + 109 * ncomps + c);
        const auto *hd_1_110 = buffer.data(hd_1 + 110 * ncomps + c);
        const auto *hd_1_111 = buffer.data(hd_1 + 111 * ncomps + c);
        const auto *hd_1_112 = buffer.data(hd_1 + 112 * ncomps + c);
        const auto *hd_1_113 = buffer.data(hd_1 + 113 * ncomps + c);
        const auto *hd_1_114 = buffer.data(hd_1 + 114 * ncomps + c);
        const auto *hd_1_115 = buffer.data(hd_1 + 115 * ncomps + c);
        const auto *hd_1_116 = buffer.data(hd_1 + 116 * ncomps + c);
        const auto *hd_1_117 = buffer.data(hd_1 + 117 * ncomps + c);
        const auto *hd_1_118 = buffer.data(hd_1 + 118 * ncomps + c);
        const auto *hd_1_119 = buffer.data(hd_1 + 119 * ncomps + c);
        const auto *hd_1_120 = buffer.data(hd_1 + 120 * ncomps + c);
        const auto *hd_1_121 = buffer.data(hd_1 + 121 * ncomps + c);
        const auto *hd_1_122 = buffer.data(hd_1 + 122 * ncomps + c);
        const auto *hd_1_123 = buffer.data(hd_1 + 123 * ncomps + c);
        const auto *hd_1_124 = buffer.data(hd_1 + 124 * ncomps + c);
        const auto *hd_1_125 = buffer.data(hd_1 + 125 * ncomps + c);

        const auto *hd_0_71 = buffer.data(hd_0 + 71 * ncomps + c);
        const auto *hd_0_72 = buffer.data(hd_0 + 72 * ncomps + c);
        const auto *hd_0_73 = buffer.data(hd_0 + 73 * ncomps + c);
        const auto *hd_0_74 = buffer.data(hd_0 + 74 * ncomps + c);
        const auto *hd_0_75 = buffer.data(hd_0 + 75 * ncomps + c);
        const auto *hd_0_76 = buffer.data(hd_0 + 76 * ncomps + c);
        const auto *hd_0_77 = buffer.data(hd_0 + 77 * ncomps + c);
        const auto *hd_0_78 = buffer.data(hd_0 + 78 * ncomps + c);
        const auto *hd_0_79 = buffer.data(hd_0 + 79 * ncomps + c);
        const auto *hd_0_80 = buffer.data(hd_0 + 80 * ncomps + c);
        const auto *hd_0_81 = buffer.data(hd_0 + 81 * ncomps + c);
        const auto *hd_0_82 = buffer.data(hd_0 + 82 * ncomps + c);
        const auto *hd_0_83 = buffer.data(hd_0 + 83 * ncomps + c);
        const auto *hd_0_84 = buffer.data(hd_0 + 84 * ncomps + c);
        const auto *hd_0_85 = buffer.data(hd_0 + 85 * ncomps + c);
        const auto *hd_0_86 = buffer.data(hd_0 + 86 * ncomps + c);
        const auto *hd_0_87 = buffer.data(hd_0 + 87 * ncomps + c);
        const auto *hd_0_88 = buffer.data(hd_0 + 88 * ncomps + c);
        const auto *hd_0_89 = buffer.data(hd_0 + 89 * ncomps + c);
        const auto *hd_0_90 = buffer.data(hd_0 + 90 * ncomps + c);
        const auto *hd_0_91 = buffer.data(hd_0 + 91 * ncomps + c);
        const auto *hd_0_92 = buffer.data(hd_0 + 92 * ncomps + c);
        const auto *hd_0_93 = buffer.data(hd_0 + 93 * ncomps + c);
        const auto *hd_0_94 = buffer.data(hd_0 + 94 * ncomps + c);
        const auto *hd_0_95 = buffer.data(hd_0 + 95 * ncomps + c);
        const auto *hd_0_96 = buffer.data(hd_0 + 96 * ncomps + c);
        const auto *hd_0_97 = buffer.data(hd_0 + 97 * ncomps + c);
        const auto *hd_0_98 = buffer.data(hd_0 + 98 * ncomps + c);
        const auto *hd_0_99 = buffer.data(hd_0 + 99 * ncomps + c);
        const auto *hd_0_100 = buffer.data(hd_0 + 100 * ncomps + c);
        const auto *hd_0_101 = buffer.data(hd_0 + 101 * ncomps + c);
        const auto *hd_0_102 = buffer.data(hd_0 + 102 * ncomps + c);
        const auto *hd_0_103 = buffer.data(hd_0 + 103 * ncomps + c);
        const auto *hd_0_104 = buffer.data(hd_0 + 104 * ncomps + c);
        const auto *hd_0_105 = buffer.data(hd_0 + 105 * ncomps + c);
        const auto *hd_0_106 = buffer.data(hd_0 + 106 * ncomps + c);
        const auto *hd_0_107 = buffer.data(hd_0 + 107 * ncomps + c);
        const auto *hd_0_108 = buffer.data(hd_0 + 108 * ncomps + c);
        const auto *hd_0_109 = buffer.data(hd_0 + 109 * ncomps + c);
        const auto *hd_0_110 = buffer.data(hd_0 + 110 * ncomps + c);
        const auto *hd_0_111 = buffer.data(hd_0 + 111 * ncomps + c);
        const auto *hd_0_112 = buffer.data(hd_0 + 112 * ncomps + c);
        const auto *hd_0_113 = buffer.data(hd_0 + 113 * ncomps + c);
        const auto *hd_0_114 = buffer.data(hd_0 + 114 * ncomps + c);
        const auto *hd_0_115 = buffer.data(hd_0 + 115 * ncomps + c);
        const auto *hd_0_116 = buffer.data(hd_0 + 116 * ncomps + c);
        const auto *hd_0_117 = buffer.data(hd_0 + 117 * ncomps + c);
        const auto *hd_0_118 = buffer.data(hd_0 + 118 * ncomps + c);
        const auto *hd_0_119 = buffer.data(hd_0 + 119 * ncomps + c);
        const auto *hd_0_120 = buffer.data(hd_0 + 120 * ncomps + c);
        const auto *hd_0_121 = buffer.data(hd_0 + 121 * ncomps + c);
        const auto *hd_0_122 = buffer.data(hd_0 + 122 * ncomps + c);
        const auto *hd_0_123 = buffer.data(hd_0 + 123 * ncomps + c);
        const auto *hd_0_124 = buffer.data(hd_0 + 124 * ncomps + c);
        const auto *hd_0_125 = buffer.data(hd_0 + 125 * ncomps + c);

        const auto *id_1_71 = buffer.data(id_1 + 71 * ncomps + c);
        const auto *id_1_72 = buffer.data(id_1 + 72 * ncomps + c);
        const auto *id_1_73 = buffer.data(id_1 + 73 * ncomps + c);
        const auto *id_1_74 = buffer.data(id_1 + 74 * ncomps + c);
        const auto *id_1_75 = buffer.data(id_1 + 75 * ncomps + c);
        const auto *id_1_76 = buffer.data(id_1 + 76 * ncomps + c);
        const auto *id_1_77 = buffer.data(id_1 + 77 * ncomps + c);
        const auto *id_1_78 = buffer.data(id_1 + 78 * ncomps + c);
        const auto *id_1_79 = buffer.data(id_1 + 79 * ncomps + c);
        const auto *id_1_80 = buffer.data(id_1 + 80 * ncomps + c);
        const auto *id_1_81 = buffer.data(id_1 + 81 * ncomps + c);
        const auto *id_1_82 = buffer.data(id_1 + 82 * ncomps + c);
        const auto *id_1_83 = buffer.data(id_1 + 83 * ncomps + c);
        const auto *id_1_84 = buffer.data(id_1 + 84 * ncomps + c);
        const auto *id_1_85 = buffer.data(id_1 + 85 * ncomps + c);
        const auto *id_1_86 = buffer.data(id_1 + 86 * ncomps + c);
        const auto *id_1_87 = buffer.data(id_1 + 87 * ncomps + c);
        const auto *id_1_88 = buffer.data(id_1 + 88 * ncomps + c);
        const auto *id_1_89 = buffer.data(id_1 + 89 * ncomps + c);
        const auto *id_1_90 = buffer.data(id_1 + 90 * ncomps + c);
        const auto *id_1_91 = buffer.data(id_1 + 91 * ncomps + c);
        const auto *id_1_92 = buffer.data(id_1 + 92 * ncomps + c);
        const auto *id_1_93 = buffer.data(id_1 + 93 * ncomps + c);
        const auto *id_1_94 = buffer.data(id_1 + 94 * ncomps + c);
        const auto *id_1_95 = buffer.data(id_1 + 95 * ncomps + c);
        const auto *id_1_96 = buffer.data(id_1 + 96 * ncomps + c);
        const auto *id_1_97 = buffer.data(id_1 + 97 * ncomps + c);
        const auto *id_1_98 = buffer.data(id_1 + 98 * ncomps + c);
        const auto *id_1_99 = buffer.data(id_1 + 99 * ncomps + c);
        const auto *id_1_100 = buffer.data(id_1 + 100 * ncomps + c);
        const auto *id_1_101 = buffer.data(id_1 + 101 * ncomps + c);
        const auto *id_1_102 = buffer.data(id_1 + 102 * ncomps + c);
        const auto *id_1_103 = buffer.data(id_1 + 103 * ncomps + c);
        const auto *id_1_104 = buffer.data(id_1 + 104 * ncomps + c);
        const auto *id_1_105 = buffer.data(id_1 + 105 * ncomps + c);
        const auto *id_1_106 = buffer.data(id_1 + 106 * ncomps + c);
        const auto *id_1_107 = buffer.data(id_1 + 107 * ncomps + c);
        const auto *id_1_108 = buffer.data(id_1 + 108 * ncomps + c);
        const auto *id_1_109 = buffer.data(id_1 + 109 * ncomps + c);
        const auto *id_1_110 = buffer.data(id_1 + 110 * ncomps + c);
        const auto *id_1_111 = buffer.data(id_1 + 111 * ncomps + c);
        const auto *id_1_112 = buffer.data(id_1 + 112 * ncomps + c);
        const auto *id_1_113 = buffer.data(id_1 + 113 * ncomps + c);
        const auto *id_1_114 = buffer.data(id_1 + 114 * ncomps + c);
        const auto *id_1_115 = buffer.data(id_1 + 115 * ncomps + c);
        const auto *id_1_116 = buffer.data(id_1 + 116 * ncomps + c);
        const auto *id_1_117 = buffer.data(id_1 + 117 * ncomps + c);
        const auto *id_1_118 = buffer.data(id_1 + 118 * ncomps + c);
        const auto *id_1_119 = buffer.data(id_1 + 119 * ncomps + c);
        const auto *id_1_120 = buffer.data(id_1 + 120 * ncomps + c);
        const auto *id_1_121 = buffer.data(id_1 + 121 * ncomps + c);
        const auto *id_1_122 = buffer.data(id_1 + 122 * ncomps + c);
        const auto *id_1_123 = buffer.data(id_1 + 123 * ncomps + c);
        const auto *id_1_124 = buffer.data(id_1 + 124 * ncomps + c);
        const auto *id_1_125 = buffer.data(id_1 + 125 * ncomps + c);
        const auto *id_1_129 = buffer.data(id_1 + 129 * ncomps + c);
        const auto *id_1_130 = buffer.data(id_1 + 130 * ncomps + c);
        const auto *id_1_131 = buffer.data(id_1 + 131 * ncomps + c);
        const auto *id_1_135 = buffer.data(id_1 + 135 * ncomps + c);
        const auto *id_1_136 = buffer.data(id_1 + 136 * ncomps + c);
        const auto *id_1_137 = buffer.data(id_1 + 137 * ncomps + c);
        const auto *id_1_141 = buffer.data(id_1 + 141 * ncomps + c);
        const auto *id_1_142 = buffer.data(id_1 + 142 * ncomps + c);
        const auto *id_1_143 = buffer.data(id_1 + 143 * ncomps + c);
        const auto *id_1_147 = buffer.data(id_1 + 147 * ncomps + c);
        const auto *id_1_148 = buffer.data(id_1 + 148 * ncomps + c);
        const auto *id_1_149 = buffer.data(id_1 + 149 * ncomps + c);
        const auto *id_1_153 = buffer.data(id_1 + 153 * ncomps + c);
        const auto *id_1_154 = buffer.data(id_1 + 154 * ncomps + c);
        const auto *id_1_155 = buffer.data(id_1 + 155 * ncomps + c);
        const auto *id_1_159 = buffer.data(id_1 + 159 * ncomps + c);
        const auto *id_1_160 = buffer.data(id_1 + 160 * ncomps + c);
        const auto *id_1_161 = buffer.data(id_1 + 161 * ncomps + c);
        const auto *id_1_167 = buffer.data(id_1 + 167 * ncomps + c);

#pragma omp simd aligned(t_115, t_116, t_117, t_118, ab_x, ab_y, hd_1_69, hd_1_70, hd_1_71, \
                         hd_0_71, id_1_71, id_1_99, id_1_100, \
                         id_1_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = ab_x[k] * hd_1_71[k]
                       + hd_0_71[k]
                       + id_1_71[k];

            t_116[k] = ab_y[k] * hd_1_69[k]
                       + id_1_99[k];

            t_117[k] = ab_y[k] * hd_1_70[k]
                       + id_1_100[k];

            t_118[k] = ab_y[k] * hd_1_71[k]
                       + id_1_101[k];
        }

#pragma omp simd aligned(t_119, t_120, t_121, ab_x, ab_z, hd_1_71, hd_1_72, hd_1_73, hd_0_72, \
                         hd_0_73, id_1_72, id_1_73, id_1_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_119[k] = ab_z[k] * hd_1_71[k]
                       + id_1_107[k];

            t_120[k] = ab_x[k] * hd_1_72[k]
                       + hd_0_72[k]
                       + id_1_72[k];

            t_121[k] = ab_x[k] * hd_1_73[k]
                       + hd_0_73[k]
                       + id_1_73[k];
        }

#pragma omp simd aligned(t_122, t_123, t_124, ab_x, hd_1_74, hd_1_75, hd_1_76, hd_0_74, \
                         hd_0_75, hd_0_76, id_1_74, id_1_75, id_1_76 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_122[k] = ab_x[k] * hd_1_74[k]
                       + hd_0_74[k]
                       + id_1_74[k];

            t_123[k] = ab_x[k] * hd_1_75[k]
                       + hd_0_75[k]
                       + id_1_75[k];

            t_124[k] = ab_x[k] * hd_1_76[k]
                       + hd_0_76[k]
                       + id_1_76[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, ab_x, ab_y, hd_1_75, hd_1_76, hd_1_77, \
                         hd_0_77, id_1_77, id_1_105, id_1_106, \
                         id_1_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = ab_x[k] * hd_1_77[k]
                       + hd_0_77[k]
                       + id_1_77[k];

            t_126[k] = ab_y[k] * hd_1_75[k]
                       + id_1_105[k];

            t_127[k] = ab_y[k] * hd_1_76[k]
                       + id_1_106[k];

            t_128[k] = ab_y[k] * hd_1_77[k]
                       + id_1_107[k];
        }

#pragma omp simd aligned(t_129, t_130, t_131, ab_x, ab_z, hd_1_77, hd_1_78, hd_1_79, hd_0_78, \
                         hd_0_79, id_1_78, id_1_79, id_1_113 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_129[k] = ab_z[k] * hd_1_77[k]
                       + id_1_113[k];

            t_130[k] = ab_x[k] * hd_1_78[k]
                       + hd_0_78[k]
                       + id_1_78[k];

            t_131[k] = ab_x[k] * hd_1_79[k]
                       + hd_0_79[k]
                       + id_1_79[k];
        }

#pragma omp simd aligned(t_132, t_133, t_134, ab_x, hd_1_80, hd_1_81, hd_1_82, hd_0_80, \
                         hd_0_81, hd_0_82, id_1_80, id_1_81, id_1_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_132[k] = ab_x[k] * hd_1_80[k]
                       + hd_0_80[k]
                       + id_1_80[k];

            t_133[k] = ab_x[k] * hd_1_81[k]
                       + hd_0_81[k]
                       + id_1_81[k];

            t_134[k] = ab_x[k] * hd_1_82[k]
                       + hd_0_82[k]
                       + id_1_82[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, ab_x, ab_y, hd_1_81, hd_1_82, hd_1_83, \
                         hd_0_83, id_1_83, id_1_111, id_1_112, \
                         id_1_113 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = ab_x[k] * hd_1_83[k]
                       + hd_0_83[k]
                       + id_1_83[k];

            t_136[k] = ab_y[k] * hd_1_81[k]
                       + id_1_111[k];

            t_137[k] = ab_y[k] * hd_1_82[k]
                       + id_1_112[k];

            t_138[k] = ab_y[k] * hd_1_83[k]
                       + id_1_113[k];
        }

#pragma omp simd aligned(t_139, t_140, t_141, ab_x, ab_z, hd_1_83, hd_1_84, hd_1_85, hd_0_84, \
                         hd_0_85, id_1_84, id_1_85, id_1_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_139[k] = ab_z[k] * hd_1_83[k]
                       + id_1_119[k];

            t_140[k] = ab_x[k] * hd_1_84[k]
                       + hd_0_84[k]
                       + id_1_84[k];

            t_141[k] = ab_x[k] * hd_1_85[k]
                       + hd_0_85[k]
                       + id_1_85[k];
        }

#pragma omp simd aligned(t_142, t_143, t_144, ab_x, hd_1_86, hd_1_87, hd_1_88, hd_0_86, \
                         hd_0_87, hd_0_88, id_1_86, id_1_87, id_1_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_142[k] = ab_x[k] * hd_1_86[k]
                       + hd_0_86[k]
                       + id_1_86[k];

            t_143[k] = ab_x[k] * hd_1_87[k]
                       + hd_0_87[k]
                       + id_1_87[k];

            t_144[k] = ab_x[k] * hd_1_88[k]
                       + hd_0_88[k]
                       + id_1_88[k];
        }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, ab_x, ab_y, hd_1_87, hd_1_88, hd_1_89, \
                         hd_0_89, id_1_89, id_1_117, id_1_118, \
                         id_1_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = ab_x[k] * hd_1_89[k]
                       + hd_0_89[k]
                       + id_1_89[k];

            t_146[k] = ab_y[k] * hd_1_87[k]
                       + id_1_117[k];

            t_147[k] = ab_y[k] * hd_1_88[k]
                       + id_1_118[k];

            t_148[k] = ab_y[k] * hd_1_89[k]
                       + id_1_119[k];
        }

#pragma omp simd aligned(t_149, t_150, t_151, ab_x, ab_z, hd_1_89, hd_1_90, hd_1_91, hd_0_90, \
                         hd_0_91, id_1_90, id_1_91, id_1_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_149[k] = ab_z[k] * hd_1_89[k]
                       + id_1_125[k];

            t_150[k] = ab_x[k] * hd_1_90[k]
                       + hd_0_90[k]
                       + id_1_90[k];

            t_151[k] = ab_x[k] * hd_1_91[k]
                       + hd_0_91[k]
                       + id_1_91[k];
        }

#pragma omp simd aligned(t_152, t_153, t_154, ab_x, hd_1_92, hd_1_93, hd_1_94, hd_0_92, \
                         hd_0_93, hd_0_94, id_1_92, id_1_93, id_1_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_152[k] = ab_x[k] * hd_1_92[k]
                       + hd_0_92[k]
                       + id_1_92[k];

            t_153[k] = ab_x[k] * hd_1_93[k]
                       + hd_0_93[k]
                       + id_1_93[k];

            t_154[k] = ab_x[k] * hd_1_94[k]
                       + hd_0_94[k]
                       + id_1_94[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, ab_x, ab_y, hd_1_93, hd_1_94, hd_1_95, \
                         hd_0_95, id_1_95, id_1_129, id_1_130, \
                         id_1_131 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * hd_1_95[k]
                       + hd_0_95[k]
                       + id_1_95[k];

            t_156[k] = ab_y[k] * hd_1_93[k]
                       + id_1_129[k];

            t_157[k] = ab_y[k] * hd_1_94[k]
                       + id_1_130[k];

            t_158[k] = ab_y[k] * hd_1_95[k]
                       + id_1_131[k];
        }

#pragma omp simd aligned(t_159, t_160, t_161, ab_x, ab_z, hd_1_95, hd_1_96, hd_1_97, hd_0_96, \
                         hd_0_97, id_1_96, id_1_97, id_1_137 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_159[k] = ab_z[k] * hd_1_95[k]
                       + id_1_137[k];

            t_160[k] = ab_x[k] * hd_1_96[k]
                       + hd_0_96[k]
                       + id_1_96[k];

            t_161[k] = ab_x[k] * hd_1_97[k]
                       + hd_0_97[k]
                       + id_1_97[k];
        }

#pragma omp simd aligned(t_162, t_163, t_164, ab_x, hd_1_98, hd_1_99, hd_1_100, hd_0_98, \
                         hd_0_99, hd_0_100, id_1_98, id_1_99, \
                         id_1_100 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_162[k] = ab_x[k] * hd_1_98[k]
                       + hd_0_98[k]
                       + id_1_98[k];

            t_163[k] = ab_x[k] * hd_1_99[k]
                       + hd_0_99[k]
                       + id_1_99[k];

            t_164[k] = ab_x[k] * hd_1_100[k]
                       + hd_0_100[k]
                       + id_1_100[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, ab_x, ab_y, hd_1_99, hd_1_100, hd_1_101, \
                         hd_0_101, id_1_101, id_1_135, id_1_136, \
                         id_1_137 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = ab_x[k] * hd_1_101[k]
                       + hd_0_101[k]
                       + id_1_101[k];

            t_166[k] = ab_y[k] * hd_1_99[k]
                       + id_1_135[k];

            t_167[k] = ab_y[k] * hd_1_100[k]
                       + id_1_136[k];

            t_168[k] = ab_y[k] * hd_1_101[k]
                       + id_1_137[k];
        }

#pragma omp simd aligned(t_169, t_170, t_171, ab_x, ab_z, hd_1_101, hd_1_102, hd_1_103, \
                         hd_0_102, hd_0_103, id_1_102, id_1_103, \
                         id_1_143 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_169[k] = ab_z[k] * hd_1_101[k]
                       + id_1_143[k];

            t_170[k] = ab_x[k] * hd_1_102[k]
                       + hd_0_102[k]
                       + id_1_102[k];

            t_171[k] = ab_x[k] * hd_1_103[k]
                       + hd_0_103[k]
                       + id_1_103[k];
        }

#pragma omp simd aligned(t_172, t_173, t_174, ab_x, hd_1_104, hd_1_105, hd_1_106, hd_0_104, \
                         hd_0_105, hd_0_106, id_1_104, id_1_105, \
                         id_1_106 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_172[k] = ab_x[k] * hd_1_104[k]
                       + hd_0_104[k]
                       + id_1_104[k];

            t_173[k] = ab_x[k] * hd_1_105[k]
                       + hd_0_105[k]
                       + id_1_105[k];

            t_174[k] = ab_x[k] * hd_1_106[k]
                       + hd_0_106[k]
                       + id_1_106[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, ab_x, ab_y, hd_1_105, hd_1_106, hd_1_107, \
                         hd_0_107, id_1_107, id_1_141, id_1_142, \
                         id_1_143 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = ab_x[k] * hd_1_107[k]
                       + hd_0_107[k]
                       + id_1_107[k];

            t_176[k] = ab_y[k] * hd_1_105[k]
                       + id_1_141[k];

            t_177[k] = ab_y[k] * hd_1_106[k]
                       + id_1_142[k];

            t_178[k] = ab_y[k] * hd_1_107[k]
                       + id_1_143[k];
        }

#pragma omp simd aligned(t_179, t_180, t_181, ab_x, ab_z, hd_1_107, hd_1_108, hd_1_109, \
                         hd_0_108, hd_0_109, id_1_108, id_1_109, \
                         id_1_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_179[k] = ab_z[k] * hd_1_107[k]
                       + id_1_149[k];

            t_180[k] = ab_x[k] * hd_1_108[k]
                       + hd_0_108[k]
                       + id_1_108[k];

            t_181[k] = ab_x[k] * hd_1_109[k]
                       + hd_0_109[k]
                       + id_1_109[k];
        }

#pragma omp simd aligned(t_182, t_183, t_184, ab_x, hd_1_110, hd_1_111, hd_1_112, hd_0_110, \
                         hd_0_111, hd_0_112, id_1_110, id_1_111, \
                         id_1_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_182[k] = ab_x[k] * hd_1_110[k]
                       + hd_0_110[k]
                       + id_1_110[k];

            t_183[k] = ab_x[k] * hd_1_111[k]
                       + hd_0_111[k]
                       + id_1_111[k];

            t_184[k] = ab_x[k] * hd_1_112[k]
                       + hd_0_112[k]
                       + id_1_112[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, ab_x, ab_y, hd_1_111, hd_1_112, hd_1_113, \
                         hd_0_113, id_1_113, id_1_147, id_1_148, \
                         id_1_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_x[k] * hd_1_113[k]
                       + hd_0_113[k]
                       + id_1_113[k];

            t_186[k] = ab_y[k] * hd_1_111[k]
                       + id_1_147[k];

            t_187[k] = ab_y[k] * hd_1_112[k]
                       + id_1_148[k];

            t_188[k] = ab_y[k] * hd_1_113[k]
                       + id_1_149[k];
        }

#pragma omp simd aligned(t_189, t_190, t_191, ab_x, ab_z, hd_1_113, hd_1_114, hd_1_115, \
                         hd_0_114, hd_0_115, id_1_114, id_1_115, \
                         id_1_155 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_189[k] = ab_z[k] * hd_1_113[k]
                       + id_1_155[k];

            t_190[k] = ab_x[k] * hd_1_114[k]
                       + hd_0_114[k]
                       + id_1_114[k];

            t_191[k] = ab_x[k] * hd_1_115[k]
                       + hd_0_115[k]
                       + id_1_115[k];
        }

#pragma omp simd aligned(t_192, t_193, t_194, ab_x, hd_1_116, hd_1_117, hd_1_118, hd_0_116, \
                         hd_0_117, hd_0_118, id_1_116, id_1_117, \
                         id_1_118 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_192[k] = ab_x[k] * hd_1_116[k]
                       + hd_0_116[k]
                       + id_1_116[k];

            t_193[k] = ab_x[k] * hd_1_117[k]
                       + hd_0_117[k]
                       + id_1_117[k];

            t_194[k] = ab_x[k] * hd_1_118[k]
                       + hd_0_118[k]
                       + id_1_118[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, ab_x, ab_y, hd_1_117, hd_1_118, hd_1_119, \
                         hd_0_119, id_1_119, id_1_153, id_1_154, \
                         id_1_155 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = ab_x[k] * hd_1_119[k]
                       + hd_0_119[k]
                       + id_1_119[k];

            t_196[k] = ab_y[k] * hd_1_117[k]
                       + id_1_153[k];

            t_197[k] = ab_y[k] * hd_1_118[k]
                       + id_1_154[k];

            t_198[k] = ab_y[k] * hd_1_119[k]
                       + id_1_155[k];
        }

#pragma omp simd aligned(t_199, t_200, t_201, ab_x, ab_z, hd_1_119, hd_1_120, hd_1_121, \
                         hd_0_120, hd_0_121, id_1_120, id_1_121, \
                         id_1_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_199[k] = ab_z[k] * hd_1_119[k]
                       + id_1_161[k];

            t_200[k] = ab_x[k] * hd_1_120[k]
                       + hd_0_120[k]
                       + id_1_120[k];

            t_201[k] = ab_x[k] * hd_1_121[k]
                       + hd_0_121[k]
                       + id_1_121[k];
        }

#pragma omp simd aligned(t_202, t_203, t_204, ab_x, hd_1_122, hd_1_123, hd_1_124, hd_0_122, \
                         hd_0_123, hd_0_124, id_1_122, id_1_123, \
                         id_1_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_202[k] = ab_x[k] * hd_1_122[k]
                       + hd_0_122[k]
                       + id_1_122[k];

            t_203[k] = ab_x[k] * hd_1_123[k]
                       + hd_0_123[k]
                       + id_1_123[k];

            t_204[k] = ab_x[k] * hd_1_124[k]
                       + hd_0_124[k]
                       + id_1_124[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, ab_x, ab_y, hd_1_123, hd_1_124, hd_1_125, \
                         hd_0_125, id_1_125, id_1_159, id_1_160, \
                         id_1_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = ab_x[k] * hd_1_125[k]
                       + hd_0_125[k]
                       + id_1_125[k];

            t_206[k] = ab_y[k] * hd_1_123[k]
                       + id_1_159[k];

            t_207[k] = ab_y[k] * hd_1_124[k]
                       + id_1_160[k];

            t_208[k] = ab_y[k] * hd_1_125[k]
                       + id_1_161[k];
        }

#pragma omp simd aligned(t_209, ab_z, hd_1_125, id_1_167 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_209[k] = ab_z[k] * hd_1_125[k]
                       + id_1_167[k];
        }
    }
}

auto
compute_hrr_geom_100x_hf_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                      const size_t target, const size_t hd_1, const size_t hd_0,
                                      const size_t id_1, const size_t ncomps,
                                      const size_t nmax) -> void
{
    compute_hrr_geom_100x_hf_out_of_first_piece0(buffer, coordinates, target, hd_1, hd_0, id_1,
                                                 ncomps, nmax);

    compute_hrr_geom_100x_hf_out_of_first_piece1(buffer, coordinates, target, hd_1, hd_0, id_1,
                                                 ncomps, nmax);
}

}  // namespace simdtrf
