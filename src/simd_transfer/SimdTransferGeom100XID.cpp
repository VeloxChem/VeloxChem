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


#include "SimdTransferGeom100XID.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_geom_100x_id_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                             const size_t target, const size_t ip_1,
                                             const size_t ip_0, const size_t kp_1,
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
        auto *t_115 = buffer.data(target + 115 * ncomps + c);
        auto *t_116 = buffer.data(target + 116 * ncomps + c);
        auto *t_117 = buffer.data(target + 117 * ncomps + c);
        auto *t_118 = buffer.data(target + 118 * ncomps + c);
        auto *t_119 = buffer.data(target + 119 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ip_1_0 = buffer.data(ip_1 + 0 * ncomps + c);
        const auto *ip_1_1 = buffer.data(ip_1 + 1 * ncomps + c);
        const auto *ip_1_2 = buffer.data(ip_1 + 2 * ncomps + c);
        const auto *ip_1_3 = buffer.data(ip_1 + 3 * ncomps + c);
        const auto *ip_1_4 = buffer.data(ip_1 + 4 * ncomps + c);
        const auto *ip_1_5 = buffer.data(ip_1 + 5 * ncomps + c);
        const auto *ip_1_6 = buffer.data(ip_1 + 6 * ncomps + c);
        const auto *ip_1_7 = buffer.data(ip_1 + 7 * ncomps + c);
        const auto *ip_1_8 = buffer.data(ip_1 + 8 * ncomps + c);
        const auto *ip_1_9 = buffer.data(ip_1 + 9 * ncomps + c);
        const auto *ip_1_10 = buffer.data(ip_1 + 10 * ncomps + c);
        const auto *ip_1_11 = buffer.data(ip_1 + 11 * ncomps + c);
        const auto *ip_1_12 = buffer.data(ip_1 + 12 * ncomps + c);
        const auto *ip_1_13 = buffer.data(ip_1 + 13 * ncomps + c);
        const auto *ip_1_14 = buffer.data(ip_1 + 14 * ncomps + c);
        const auto *ip_1_15 = buffer.data(ip_1 + 15 * ncomps + c);
        const auto *ip_1_16 = buffer.data(ip_1 + 16 * ncomps + c);
        const auto *ip_1_17 = buffer.data(ip_1 + 17 * ncomps + c);
        const auto *ip_1_18 = buffer.data(ip_1 + 18 * ncomps + c);
        const auto *ip_1_19 = buffer.data(ip_1 + 19 * ncomps + c);
        const auto *ip_1_20 = buffer.data(ip_1 + 20 * ncomps + c);
        const auto *ip_1_21 = buffer.data(ip_1 + 21 * ncomps + c);
        const auto *ip_1_22 = buffer.data(ip_1 + 22 * ncomps + c);
        const auto *ip_1_23 = buffer.data(ip_1 + 23 * ncomps + c);
        const auto *ip_1_24 = buffer.data(ip_1 + 24 * ncomps + c);
        const auto *ip_1_25 = buffer.data(ip_1 + 25 * ncomps + c);
        const auto *ip_1_26 = buffer.data(ip_1 + 26 * ncomps + c);
        const auto *ip_1_27 = buffer.data(ip_1 + 27 * ncomps + c);
        const auto *ip_1_28 = buffer.data(ip_1 + 28 * ncomps + c);
        const auto *ip_1_29 = buffer.data(ip_1 + 29 * ncomps + c);
        const auto *ip_1_30 = buffer.data(ip_1 + 30 * ncomps + c);
        const auto *ip_1_31 = buffer.data(ip_1 + 31 * ncomps + c);
        const auto *ip_1_32 = buffer.data(ip_1 + 32 * ncomps + c);
        const auto *ip_1_33 = buffer.data(ip_1 + 33 * ncomps + c);
        const auto *ip_1_34 = buffer.data(ip_1 + 34 * ncomps + c);
        const auto *ip_1_35 = buffer.data(ip_1 + 35 * ncomps + c);
        const auto *ip_1_36 = buffer.data(ip_1 + 36 * ncomps + c);
        const auto *ip_1_37 = buffer.data(ip_1 + 37 * ncomps + c);
        const auto *ip_1_38 = buffer.data(ip_1 + 38 * ncomps + c);
        const auto *ip_1_39 = buffer.data(ip_1 + 39 * ncomps + c);
        const auto *ip_1_40 = buffer.data(ip_1 + 40 * ncomps + c);
        const auto *ip_1_41 = buffer.data(ip_1 + 41 * ncomps + c);
        const auto *ip_1_42 = buffer.data(ip_1 + 42 * ncomps + c);
        const auto *ip_1_43 = buffer.data(ip_1 + 43 * ncomps + c);
        const auto *ip_1_44 = buffer.data(ip_1 + 44 * ncomps + c);
        const auto *ip_1_45 = buffer.data(ip_1 + 45 * ncomps + c);
        const auto *ip_1_46 = buffer.data(ip_1 + 46 * ncomps + c);
        const auto *ip_1_47 = buffer.data(ip_1 + 47 * ncomps + c);
        const auto *ip_1_48 = buffer.data(ip_1 + 48 * ncomps + c);
        const auto *ip_1_49 = buffer.data(ip_1 + 49 * ncomps + c);
        const auto *ip_1_50 = buffer.data(ip_1 + 50 * ncomps + c);
        const auto *ip_1_51 = buffer.data(ip_1 + 51 * ncomps + c);
        const auto *ip_1_52 = buffer.data(ip_1 + 52 * ncomps + c);
        const auto *ip_1_53 = buffer.data(ip_1 + 53 * ncomps + c);
        const auto *ip_1_54 = buffer.data(ip_1 + 54 * ncomps + c);
        const auto *ip_1_55 = buffer.data(ip_1 + 55 * ncomps + c);
        const auto *ip_1_56 = buffer.data(ip_1 + 56 * ncomps + c);
        const auto *ip_1_57 = buffer.data(ip_1 + 57 * ncomps + c);
        const auto *ip_1_58 = buffer.data(ip_1 + 58 * ncomps + c);
        const auto *ip_1_59 = buffer.data(ip_1 + 59 * ncomps + c);

        const auto *ip_0_0 = buffer.data(ip_0 + 0 * ncomps + c);
        const auto *ip_0_1 = buffer.data(ip_0 + 1 * ncomps + c);
        const auto *ip_0_2 = buffer.data(ip_0 + 2 * ncomps + c);
        const auto *ip_0_3 = buffer.data(ip_0 + 3 * ncomps + c);
        const auto *ip_0_4 = buffer.data(ip_0 + 4 * ncomps + c);
        const auto *ip_0_5 = buffer.data(ip_0 + 5 * ncomps + c);
        const auto *ip_0_6 = buffer.data(ip_0 + 6 * ncomps + c);
        const auto *ip_0_7 = buffer.data(ip_0 + 7 * ncomps + c);
        const auto *ip_0_8 = buffer.data(ip_0 + 8 * ncomps + c);
        const auto *ip_0_9 = buffer.data(ip_0 + 9 * ncomps + c);
        const auto *ip_0_10 = buffer.data(ip_0 + 10 * ncomps + c);
        const auto *ip_0_11 = buffer.data(ip_0 + 11 * ncomps + c);
        const auto *ip_0_12 = buffer.data(ip_0 + 12 * ncomps + c);
        const auto *ip_0_13 = buffer.data(ip_0 + 13 * ncomps + c);
        const auto *ip_0_14 = buffer.data(ip_0 + 14 * ncomps + c);
        const auto *ip_0_15 = buffer.data(ip_0 + 15 * ncomps + c);
        const auto *ip_0_16 = buffer.data(ip_0 + 16 * ncomps + c);
        const auto *ip_0_17 = buffer.data(ip_0 + 17 * ncomps + c);
        const auto *ip_0_18 = buffer.data(ip_0 + 18 * ncomps + c);
        const auto *ip_0_19 = buffer.data(ip_0 + 19 * ncomps + c);
        const auto *ip_0_20 = buffer.data(ip_0 + 20 * ncomps + c);
        const auto *ip_0_21 = buffer.data(ip_0 + 21 * ncomps + c);
        const auto *ip_0_22 = buffer.data(ip_0 + 22 * ncomps + c);
        const auto *ip_0_23 = buffer.data(ip_0 + 23 * ncomps + c);
        const auto *ip_0_24 = buffer.data(ip_0 + 24 * ncomps + c);
        const auto *ip_0_25 = buffer.data(ip_0 + 25 * ncomps + c);
        const auto *ip_0_26 = buffer.data(ip_0 + 26 * ncomps + c);
        const auto *ip_0_27 = buffer.data(ip_0 + 27 * ncomps + c);
        const auto *ip_0_28 = buffer.data(ip_0 + 28 * ncomps + c);
        const auto *ip_0_29 = buffer.data(ip_0 + 29 * ncomps + c);
        const auto *ip_0_30 = buffer.data(ip_0 + 30 * ncomps + c);
        const auto *ip_0_31 = buffer.data(ip_0 + 31 * ncomps + c);
        const auto *ip_0_32 = buffer.data(ip_0 + 32 * ncomps + c);
        const auto *ip_0_33 = buffer.data(ip_0 + 33 * ncomps + c);
        const auto *ip_0_34 = buffer.data(ip_0 + 34 * ncomps + c);
        const auto *ip_0_35 = buffer.data(ip_0 + 35 * ncomps + c);
        const auto *ip_0_36 = buffer.data(ip_0 + 36 * ncomps + c);
        const auto *ip_0_37 = buffer.data(ip_0 + 37 * ncomps + c);
        const auto *ip_0_38 = buffer.data(ip_0 + 38 * ncomps + c);
        const auto *ip_0_39 = buffer.data(ip_0 + 39 * ncomps + c);
        const auto *ip_0_40 = buffer.data(ip_0 + 40 * ncomps + c);
        const auto *ip_0_41 = buffer.data(ip_0 + 41 * ncomps + c);
        const auto *ip_0_42 = buffer.data(ip_0 + 42 * ncomps + c);
        const auto *ip_0_43 = buffer.data(ip_0 + 43 * ncomps + c);
        const auto *ip_0_44 = buffer.data(ip_0 + 44 * ncomps + c);
        const auto *ip_0_45 = buffer.data(ip_0 + 45 * ncomps + c);
        const auto *ip_0_46 = buffer.data(ip_0 + 46 * ncomps + c);
        const auto *ip_0_47 = buffer.data(ip_0 + 47 * ncomps + c);
        const auto *ip_0_48 = buffer.data(ip_0 + 48 * ncomps + c);
        const auto *ip_0_49 = buffer.data(ip_0 + 49 * ncomps + c);
        const auto *ip_0_50 = buffer.data(ip_0 + 50 * ncomps + c);
        const auto *ip_0_51 = buffer.data(ip_0 + 51 * ncomps + c);
        const auto *ip_0_52 = buffer.data(ip_0 + 52 * ncomps + c);
        const auto *ip_0_53 = buffer.data(ip_0 + 53 * ncomps + c);
        const auto *ip_0_54 = buffer.data(ip_0 + 54 * ncomps + c);
        const auto *ip_0_55 = buffer.data(ip_0 + 55 * ncomps + c);
        const auto *ip_0_56 = buffer.data(ip_0 + 56 * ncomps + c);
        const auto *ip_0_57 = buffer.data(ip_0 + 57 * ncomps + c);
        const auto *ip_0_58 = buffer.data(ip_0 + 58 * ncomps + c);
        const auto *ip_0_59 = buffer.data(ip_0 + 59 * ncomps + c);

        const auto *kp_1_0 = buffer.data(kp_1 + 0 * ncomps + c);
        const auto *kp_1_1 = buffer.data(kp_1 + 1 * ncomps + c);
        const auto *kp_1_2 = buffer.data(kp_1 + 2 * ncomps + c);
        const auto *kp_1_3 = buffer.data(kp_1 + 3 * ncomps + c);
        const auto *kp_1_4 = buffer.data(kp_1 + 4 * ncomps + c);
        const auto *kp_1_5 = buffer.data(kp_1 + 5 * ncomps + c);
        const auto *kp_1_6 = buffer.data(kp_1 + 6 * ncomps + c);
        const auto *kp_1_7 = buffer.data(kp_1 + 7 * ncomps + c);
        const auto *kp_1_8 = buffer.data(kp_1 + 8 * ncomps + c);
        const auto *kp_1_9 = buffer.data(kp_1 + 9 * ncomps + c);
        const auto *kp_1_10 = buffer.data(kp_1 + 10 * ncomps + c);
        const auto *kp_1_11 = buffer.data(kp_1 + 11 * ncomps + c);
        const auto *kp_1_12 = buffer.data(kp_1 + 12 * ncomps + c);
        const auto *kp_1_13 = buffer.data(kp_1 + 13 * ncomps + c);
        const auto *kp_1_14 = buffer.data(kp_1 + 14 * ncomps + c);
        const auto *kp_1_15 = buffer.data(kp_1 + 15 * ncomps + c);
        const auto *kp_1_16 = buffer.data(kp_1 + 16 * ncomps + c);
        const auto *kp_1_17 = buffer.data(kp_1 + 17 * ncomps + c);
        const auto *kp_1_18 = buffer.data(kp_1 + 18 * ncomps + c);
        const auto *kp_1_19 = buffer.data(kp_1 + 19 * ncomps + c);
        const auto *kp_1_20 = buffer.data(kp_1 + 20 * ncomps + c);
        const auto *kp_1_21 = buffer.data(kp_1 + 21 * ncomps + c);
        const auto *kp_1_22 = buffer.data(kp_1 + 22 * ncomps + c);
        const auto *kp_1_23 = buffer.data(kp_1 + 23 * ncomps + c);
        const auto *kp_1_24 = buffer.data(kp_1 + 24 * ncomps + c);
        const auto *kp_1_25 = buffer.data(kp_1 + 25 * ncomps + c);
        const auto *kp_1_26 = buffer.data(kp_1 + 26 * ncomps + c);
        const auto *kp_1_27 = buffer.data(kp_1 + 27 * ncomps + c);
        const auto *kp_1_28 = buffer.data(kp_1 + 28 * ncomps + c);
        const auto *kp_1_29 = buffer.data(kp_1 + 29 * ncomps + c);
        const auto *kp_1_30 = buffer.data(kp_1 + 30 * ncomps + c);
        const auto *kp_1_31 = buffer.data(kp_1 + 31 * ncomps + c);
        const auto *kp_1_32 = buffer.data(kp_1 + 32 * ncomps + c);
        const auto *kp_1_33 = buffer.data(kp_1 + 33 * ncomps + c);
        const auto *kp_1_34 = buffer.data(kp_1 + 34 * ncomps + c);
        const auto *kp_1_35 = buffer.data(kp_1 + 35 * ncomps + c);
        const auto *kp_1_36 = buffer.data(kp_1 + 36 * ncomps + c);
        const auto *kp_1_37 = buffer.data(kp_1 + 37 * ncomps + c);
        const auto *kp_1_38 = buffer.data(kp_1 + 38 * ncomps + c);
        const auto *kp_1_39 = buffer.data(kp_1 + 39 * ncomps + c);
        const auto *kp_1_40 = buffer.data(kp_1 + 40 * ncomps + c);
        const auto *kp_1_41 = buffer.data(kp_1 + 41 * ncomps + c);
        const auto *kp_1_42 = buffer.data(kp_1 + 42 * ncomps + c);
        const auto *kp_1_43 = buffer.data(kp_1 + 43 * ncomps + c);
        const auto *kp_1_44 = buffer.data(kp_1 + 44 * ncomps + c);
        const auto *kp_1_45 = buffer.data(kp_1 + 45 * ncomps + c);
        const auto *kp_1_46 = buffer.data(kp_1 + 46 * ncomps + c);
        const auto *kp_1_47 = buffer.data(kp_1 + 47 * ncomps + c);
        const auto *kp_1_48 = buffer.data(kp_1 + 48 * ncomps + c);
        const auto *kp_1_49 = buffer.data(kp_1 + 49 * ncomps + c);
        const auto *kp_1_50 = buffer.data(kp_1 + 50 * ncomps + c);
        const auto *kp_1_51 = buffer.data(kp_1 + 51 * ncomps + c);
        const auto *kp_1_52 = buffer.data(kp_1 + 52 * ncomps + c);
        const auto *kp_1_53 = buffer.data(kp_1 + 53 * ncomps + c);
        const auto *kp_1_54 = buffer.data(kp_1 + 54 * ncomps + c);
        const auto *kp_1_55 = buffer.data(kp_1 + 55 * ncomps + c);
        const auto *kp_1_56 = buffer.data(kp_1 + 56 * ncomps + c);
        const auto *kp_1_57 = buffer.data(kp_1 + 57 * ncomps + c);
        const auto *kp_1_58 = buffer.data(kp_1 + 58 * ncomps + c);
        const auto *kp_1_59 = buffer.data(kp_1 + 59 * ncomps + c);
        const auto *kp_1_62 = buffer.data(kp_1 + 62 * ncomps + c);
        const auto *kp_1_64 = buffer.data(kp_1 + 64 * ncomps + c);
        const auto *kp_1_65 = buffer.data(kp_1 + 65 * ncomps + c);
        const auto *kp_1_67 = buffer.data(kp_1 + 67 * ncomps + c);
        const auto *kp_1_68 = buffer.data(kp_1 + 68 * ncomps + c);
        const auto *kp_1_70 = buffer.data(kp_1 + 70 * ncomps + c);
        const auto *kp_1_71 = buffer.data(kp_1 + 71 * ncomps + c);
        const auto *kp_1_73 = buffer.data(kp_1 + 73 * ncomps + c);
        const auto *kp_1_74 = buffer.data(kp_1 + 74 * ncomps + c);
        const auto *kp_1_76 = buffer.data(kp_1 + 76 * ncomps + c);
        const auto *kp_1_77 = buffer.data(kp_1 + 77 * ncomps + c);
        const auto *kp_1_80 = buffer.data(kp_1 + 80 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, ab_x, ab_y, ip_1_0, ip_1_1, ip_1_2, ip_0_0, \
                         ip_0_1, ip_0_2, kp_1_0, kp_1_1, kp_1_2, \
                         kp_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * ip_1_0[k]
                     + ip_0_0[k]
                     + kp_1_0[k];

            t_1[k] = ab_x[k] * ip_1_1[k]
                     + ip_0_1[k]
                     + kp_1_1[k];

            t_2[k] = ab_x[k] * ip_1_2[k]
                     + ip_0_2[k]
                     + kp_1_2[k];

            t_3[k] = ab_y[k] * ip_1_1[k]
                     + kp_1_4[k];
        }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, ab_x, ab_y, ab_z, ip_1_2, ip_1_3, ip_1_4, ip_0_3, \
                         ip_0_4, kp_1_3, kp_1_4, kp_1_5, kp_1_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_4[k] = ab_y[k] * ip_1_2[k]
                     + kp_1_5[k];

            t_5[k] = ab_z[k] * ip_1_2[k]
                     + kp_1_8[k];

            t_6[k] = ab_x[k] * ip_1_3[k]
                     + ip_0_3[k]
                     + kp_1_3[k];

            t_7[k] = ab_x[k] * ip_1_4[k]
                     + ip_0_4[k]
                     + kp_1_4[k];
        }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, ab_x, ab_y, ab_z, ip_1_4, ip_1_5, ip_0_5, \
                         kp_1_5, kp_1_10, kp_1_11, kp_1_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_8[k] = ab_x[k] * ip_1_5[k]
                     + ip_0_5[k]
                     + kp_1_5[k];

            t_9[k] = ab_y[k] * ip_1_4[k]
                     + kp_1_10[k];

            t_10[k] = ab_y[k] * ip_1_5[k]
                      + kp_1_11[k];

            t_11[k] = ab_z[k] * ip_1_5[k]
                      + kp_1_14[k];
        }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, ab_x, ab_y, ip_1_6, ip_1_7, ip_1_8, ip_0_6, \
                         ip_0_7, ip_0_8, kp_1_6, kp_1_7, kp_1_8, \
                         kp_1_13 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_12[k] = ab_x[k] * ip_1_6[k]
                      + ip_0_6[k]
                      + kp_1_6[k];

            t_13[k] = ab_x[k] * ip_1_7[k]
                      + ip_0_7[k]
                      + kp_1_7[k];

            t_14[k] = ab_x[k] * ip_1_8[k]
                      + ip_0_8[k]
                      + kp_1_8[k];

            t_15[k] = ab_y[k] * ip_1_7[k]
                      + kp_1_13[k];
        }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, ab_x, ab_y, ab_z, ip_1_8, ip_1_9, ip_1_10, \
                         ip_0_9, ip_0_10, kp_1_9, kp_1_10, kp_1_14, \
                         kp_1_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_16[k] = ab_y[k] * ip_1_8[k]
                      + kp_1_14[k];

            t_17[k] = ab_z[k] * ip_1_8[k]
                      + kp_1_17[k];

            t_18[k] = ab_x[k] * ip_1_9[k]
                      + ip_0_9[k]
                      + kp_1_9[k];

            t_19[k] = ab_x[k] * ip_1_10[k]
                      + ip_0_10[k]
                      + kp_1_10[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, ab_x, ab_y, ab_z, ip_1_10, ip_1_11, ip_0_11, \
                         kp_1_11, kp_1_19, kp_1_20, kp_1_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * ip_1_11[k]
                      + ip_0_11[k]
                      + kp_1_11[k];

            t_21[k] = ab_y[k] * ip_1_10[k]
                      + kp_1_19[k];

            t_22[k] = ab_y[k] * ip_1_11[k]
                      + kp_1_20[k];

            t_23[k] = ab_z[k] * ip_1_11[k]
                      + kp_1_23[k];
        }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, ab_x, ab_y, ip_1_12, ip_1_13, ip_1_14, \
                         ip_0_12, ip_0_13, ip_0_14, kp_1_12, kp_1_13, kp_1_14, \
                         kp_1_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_24[k] = ab_x[k] * ip_1_12[k]
                      + ip_0_12[k]
                      + kp_1_12[k];

            t_25[k] = ab_x[k] * ip_1_13[k]
                      + ip_0_13[k]
                      + kp_1_13[k];

            t_26[k] = ab_x[k] * ip_1_14[k]
                      + ip_0_14[k]
                      + kp_1_14[k];

            t_27[k] = ab_y[k] * ip_1_13[k]
                      + kp_1_22[k];
        }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, ip_1_14, ip_1_15, ip_1_16, \
                         ip_0_15, ip_0_16, kp_1_15, kp_1_16, kp_1_23, \
                         kp_1_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_28[k] = ab_y[k] * ip_1_14[k]
                      + kp_1_23[k];

            t_29[k] = ab_z[k] * ip_1_14[k]
                      + kp_1_26[k];

            t_30[k] = ab_x[k] * ip_1_15[k]
                      + ip_0_15[k]
                      + kp_1_15[k];

            t_31[k] = ab_x[k] * ip_1_16[k]
                      + ip_0_16[k]
                      + kp_1_16[k];
        }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, ip_1_16, ip_1_17, ip_0_17, \
                         kp_1_17, kp_1_25, kp_1_26, kp_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_32[k] = ab_x[k] * ip_1_17[k]
                      + ip_0_17[k]
                      + kp_1_17[k];

            t_33[k] = ab_y[k] * ip_1_16[k]
                      + kp_1_25[k];

            t_34[k] = ab_y[k] * ip_1_17[k]
                      + kp_1_26[k];

            t_35[k] = ab_z[k] * ip_1_17[k]
                      + kp_1_29[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, ab_x, ab_y, ip_1_18, ip_1_19, ip_1_20, \
                         ip_0_18, ip_0_19, ip_0_20, kp_1_18, kp_1_19, kp_1_20, \
                         kp_1_31 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = ab_x[k] * ip_1_18[k]
                      + ip_0_18[k]
                      + kp_1_18[k];

            t_37[k] = ab_x[k] * ip_1_19[k]
                      + ip_0_19[k]
                      + kp_1_19[k];

            t_38[k] = ab_x[k] * ip_1_20[k]
                      + ip_0_20[k]
                      + kp_1_20[k];

            t_39[k] = ab_y[k] * ip_1_19[k]
                      + kp_1_31[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, ip_1_20, ip_1_21, ip_1_22, \
                         ip_0_21, ip_0_22, kp_1_21, kp_1_22, kp_1_32, \
                         kp_1_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_y[k] * ip_1_20[k]
                      + kp_1_32[k];

            t_41[k] = ab_z[k] * ip_1_20[k]
                      + kp_1_35[k];

            t_42[k] = ab_x[k] * ip_1_21[k]
                      + ip_0_21[k]
                      + kp_1_21[k];

            t_43[k] = ab_x[k] * ip_1_22[k]
                      + ip_0_22[k]
                      + kp_1_22[k];
        }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, ab_x, ab_y, ab_z, ip_1_22, ip_1_23, ip_0_23, \
                         kp_1_23, kp_1_34, kp_1_35, kp_1_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_44[k] = ab_x[k] * ip_1_23[k]
                      + ip_0_23[k]
                      + kp_1_23[k];

            t_45[k] = ab_y[k] * ip_1_22[k]
                      + kp_1_34[k];

            t_46[k] = ab_y[k] * ip_1_23[k]
                      + kp_1_35[k];

            t_47[k] = ab_z[k] * ip_1_23[k]
                      + kp_1_38[k];
        }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, ab_x, ab_y, ip_1_24, ip_1_25, ip_1_26, \
                         ip_0_24, ip_0_25, ip_0_26, kp_1_24, kp_1_25, kp_1_26, \
                         kp_1_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_48[k] = ab_x[k] * ip_1_24[k]
                      + ip_0_24[k]
                      + kp_1_24[k];

            t_49[k] = ab_x[k] * ip_1_25[k]
                      + ip_0_25[k]
                      + kp_1_25[k];

            t_50[k] = ab_x[k] * ip_1_26[k]
                      + ip_0_26[k]
                      + kp_1_26[k];

            t_51[k] = ab_y[k] * ip_1_25[k]
                      + kp_1_37[k];
        }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, ab_x, ab_y, ab_z, ip_1_26, ip_1_27, ip_1_28, \
                         ip_0_27, ip_0_28, kp_1_27, kp_1_28, kp_1_38, \
                         kp_1_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_52[k] = ab_y[k] * ip_1_26[k]
                      + kp_1_38[k];

            t_53[k] = ab_z[k] * ip_1_26[k]
                      + kp_1_41[k];

            t_54[k] = ab_x[k] * ip_1_27[k]
                      + ip_0_27[k]
                      + kp_1_27[k];

            t_55[k] = ab_x[k] * ip_1_28[k]
                      + ip_0_28[k]
                      + kp_1_28[k];
        }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, ip_1_28, ip_1_29, ip_0_29, \
                         kp_1_29, kp_1_40, kp_1_41, kp_1_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_56[k] = ab_x[k] * ip_1_29[k]
                      + ip_0_29[k]
                      + kp_1_29[k];

            t_57[k] = ab_y[k] * ip_1_28[k]
                      + kp_1_40[k];

            t_58[k] = ab_y[k] * ip_1_29[k]
                      + kp_1_41[k];

            t_59[k] = ab_z[k] * ip_1_29[k]
                      + kp_1_44[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, ab_x, ab_y, ip_1_30, ip_1_31, ip_1_32, \
                         ip_0_30, ip_0_31, ip_0_32, kp_1_30, kp_1_31, kp_1_32, \
                         kp_1_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * ip_1_30[k]
                      + ip_0_30[k]
                      + kp_1_30[k];

            t_61[k] = ab_x[k] * ip_1_31[k]
                      + ip_0_31[k]
                      + kp_1_31[k];

            t_62[k] = ab_x[k] * ip_1_32[k]
                      + ip_0_32[k]
                      + kp_1_32[k];

            t_63[k] = ab_y[k] * ip_1_31[k]
                      + kp_1_46[k];
        }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, ab_x, ab_y, ab_z, ip_1_32, ip_1_33, ip_1_34, \
                         ip_0_33, ip_0_34, kp_1_33, kp_1_34, kp_1_47, \
                         kp_1_50 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_64[k] = ab_y[k] * ip_1_32[k]
                      + kp_1_47[k];

            t_65[k] = ab_z[k] * ip_1_32[k]
                      + kp_1_50[k];

            t_66[k] = ab_x[k] * ip_1_33[k]
                      + ip_0_33[k]
                      + kp_1_33[k];

            t_67[k] = ab_x[k] * ip_1_34[k]
                      + ip_0_34[k]
                      + kp_1_34[k];
        }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, ab_x, ab_y, ab_z, ip_1_34, ip_1_35, ip_0_35, \
                         kp_1_35, kp_1_49, kp_1_50, kp_1_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_68[k] = ab_x[k] * ip_1_35[k]
                      + ip_0_35[k]
                      + kp_1_35[k];

            t_69[k] = ab_y[k] * ip_1_34[k]
                      + kp_1_49[k];

            t_70[k] = ab_y[k] * ip_1_35[k]
                      + kp_1_50[k];

            t_71[k] = ab_z[k] * ip_1_35[k]
                      + kp_1_53[k];
        }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, ab_x, ab_y, ip_1_36, ip_1_37, ip_1_38, \
                         ip_0_36, ip_0_37, ip_0_38, kp_1_36, kp_1_37, kp_1_38, \
                         kp_1_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_72[k] = ab_x[k] * ip_1_36[k]
                      + ip_0_36[k]
                      + kp_1_36[k];

            t_73[k] = ab_x[k] * ip_1_37[k]
                      + ip_0_37[k]
                      + kp_1_37[k];

            t_74[k] = ab_x[k] * ip_1_38[k]
                      + ip_0_38[k]
                      + kp_1_38[k];

            t_75[k] = ab_y[k] * ip_1_37[k]
                      + kp_1_52[k];
        }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, ab_x, ab_y, ab_z, ip_1_38, ip_1_39, ip_1_40, \
                         ip_0_39, ip_0_40, kp_1_39, kp_1_40, kp_1_53, \
                         kp_1_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_76[k] = ab_y[k] * ip_1_38[k]
                      + kp_1_53[k];

            t_77[k] = ab_z[k] * ip_1_38[k]
                      + kp_1_56[k];

            t_78[k] = ab_x[k] * ip_1_39[k]
                      + ip_0_39[k]
                      + kp_1_39[k];

            t_79[k] = ab_x[k] * ip_1_40[k]
                      + ip_0_40[k]
                      + kp_1_40[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, ab_x, ab_y, ab_z, ip_1_40, ip_1_41, ip_0_41, \
                         kp_1_41, kp_1_55, kp_1_56, kp_1_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_x[k] * ip_1_41[k]
                      + ip_0_41[k]
                      + kp_1_41[k];

            t_81[k] = ab_y[k] * ip_1_40[k]
                      + kp_1_55[k];

            t_82[k] = ab_y[k] * ip_1_41[k]
                      + kp_1_56[k];

            t_83[k] = ab_z[k] * ip_1_41[k]
                      + kp_1_59[k];
        }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, ab_x, ab_y, ip_1_42, ip_1_43, ip_1_44, \
                         ip_0_42, ip_0_43, ip_0_44, kp_1_42, kp_1_43, kp_1_44, \
                         kp_1_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_84[k] = ab_x[k] * ip_1_42[k]
                      + ip_0_42[k]
                      + kp_1_42[k];

            t_85[k] = ab_x[k] * ip_1_43[k]
                      + ip_0_43[k]
                      + kp_1_43[k];

            t_86[k] = ab_x[k] * ip_1_44[k]
                      + ip_0_44[k]
                      + kp_1_44[k];

            t_87[k] = ab_y[k] * ip_1_43[k]
                      + kp_1_58[k];
        }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, ab_x, ab_y, ab_z, ip_1_44, ip_1_45, ip_1_46, \
                         ip_0_45, ip_0_46, kp_1_45, kp_1_46, kp_1_59, \
                         kp_1_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_88[k] = ab_y[k] * ip_1_44[k]
                      + kp_1_59[k];

            t_89[k] = ab_z[k] * ip_1_44[k]
                      + kp_1_62[k];

            t_90[k] = ab_x[k] * ip_1_45[k]
                      + ip_0_45[k]
                      + kp_1_45[k];

            t_91[k] = ab_x[k] * ip_1_46[k]
                      + ip_0_46[k]
                      + kp_1_46[k];
        }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, ab_x, ab_y, ab_z, ip_1_46, ip_1_47, ip_0_47, \
                         kp_1_47, kp_1_64, kp_1_65, kp_1_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_92[k] = ab_x[k] * ip_1_47[k]
                      + ip_0_47[k]
                      + kp_1_47[k];

            t_93[k] = ab_y[k] * ip_1_46[k]
                      + kp_1_64[k];

            t_94[k] = ab_y[k] * ip_1_47[k]
                      + kp_1_65[k];

            t_95[k] = ab_z[k] * ip_1_47[k]
                      + kp_1_68[k];
        }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, ab_x, ab_y, ip_1_48, ip_1_49, ip_1_50, \
                         ip_0_48, ip_0_49, ip_0_50, kp_1_48, kp_1_49, kp_1_50, \
                         kp_1_67 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_96[k] = ab_x[k] * ip_1_48[k]
                      + ip_0_48[k]
                      + kp_1_48[k];

            t_97[k] = ab_x[k] * ip_1_49[k]
                      + ip_0_49[k]
                      + kp_1_49[k];

            t_98[k] = ab_x[k] * ip_1_50[k]
                      + ip_0_50[k]
                      + kp_1_50[k];

            t_99[k] = ab_y[k] * ip_1_49[k]
                      + kp_1_67[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, ab_x, ab_y, ab_z, ip_1_50, ip_1_51, \
                         ip_1_52, ip_0_51, ip_0_52, kp_1_51, kp_1_52, kp_1_68, \
                         kp_1_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = ab_y[k] * ip_1_50[k]
                       + kp_1_68[k];

            t_101[k] = ab_z[k] * ip_1_50[k]
                       + kp_1_71[k];

            t_102[k] = ab_x[k] * ip_1_51[k]
                       + ip_0_51[k]
                       + kp_1_51[k];

            t_103[k] = ab_x[k] * ip_1_52[k]
                       + ip_0_52[k]
                       + kp_1_52[k];
        }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, ab_x, ab_y, ab_z, ip_1_52, ip_1_53, \
                         ip_0_53, kp_1_53, kp_1_70, kp_1_71, kp_1_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_104[k] = ab_x[k] * ip_1_53[k]
                       + ip_0_53[k]
                       + kp_1_53[k];

            t_105[k] = ab_y[k] * ip_1_52[k]
                       + kp_1_70[k];

            t_106[k] = ab_y[k] * ip_1_53[k]
                       + kp_1_71[k];

            t_107[k] = ab_z[k] * ip_1_53[k]
                       + kp_1_74[k];
        }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, ab_x, ab_y, ip_1_54, ip_1_55, ip_1_56, \
                         ip_0_54, ip_0_55, ip_0_56, kp_1_54, kp_1_55, kp_1_56, \
                         kp_1_73 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_108[k] = ab_x[k] * ip_1_54[k]
                       + ip_0_54[k]
                       + kp_1_54[k];

            t_109[k] = ab_x[k] * ip_1_55[k]
                       + ip_0_55[k]
                       + kp_1_55[k];

            t_110[k] = ab_x[k] * ip_1_56[k]
                       + ip_0_56[k]
                       + kp_1_56[k];

            t_111[k] = ab_y[k] * ip_1_55[k]
                       + kp_1_73[k];
        }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, ab_x, ab_y, ab_z, ip_1_56, ip_1_57, \
                         ip_1_58, ip_0_57, ip_0_58, kp_1_57, kp_1_58, kp_1_74, \
                         kp_1_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_112[k] = ab_y[k] * ip_1_56[k]
                       + kp_1_74[k];

            t_113[k] = ab_z[k] * ip_1_56[k]
                       + kp_1_77[k];

            t_114[k] = ab_x[k] * ip_1_57[k]
                       + ip_0_57[k]
                       + kp_1_57[k];

            t_115[k] = ab_x[k] * ip_1_58[k]
                       + ip_0_58[k]
                       + kp_1_58[k];
        }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, ab_x, ab_y, ab_z, ip_1_58, ip_1_59, \
                         ip_0_59, kp_1_59, kp_1_76, kp_1_77, kp_1_80 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_116[k] = ab_x[k] * ip_1_59[k]
                       + ip_0_59[k]
                       + kp_1_59[k];

            t_117[k] = ab_y[k] * ip_1_58[k]
                       + kp_1_76[k];

            t_118[k] = ab_y[k] * ip_1_59[k]
                       + kp_1_77[k];

            t_119[k] = ab_z[k] * ip_1_59[k]
                       + kp_1_80[k];
        }
    }
}

static auto
compute_hrr_geom_100x_id_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                             const size_t target, const size_t ip_1,
                                             const size_t ip_0, const size_t kp_1,
                                             const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ip_1_60 = buffer.data(ip_1 + 60 * ncomps + c);
        const auto *ip_1_61 = buffer.data(ip_1 + 61 * ncomps + c);
        const auto *ip_1_62 = buffer.data(ip_1 + 62 * ncomps + c);
        const auto *ip_1_63 = buffer.data(ip_1 + 63 * ncomps + c);
        const auto *ip_1_64 = buffer.data(ip_1 + 64 * ncomps + c);
        const auto *ip_1_65 = buffer.data(ip_1 + 65 * ncomps + c);
        const auto *ip_1_66 = buffer.data(ip_1 + 66 * ncomps + c);
        const auto *ip_1_67 = buffer.data(ip_1 + 67 * ncomps + c);
        const auto *ip_1_68 = buffer.data(ip_1 + 68 * ncomps + c);
        const auto *ip_1_69 = buffer.data(ip_1 + 69 * ncomps + c);
        const auto *ip_1_70 = buffer.data(ip_1 + 70 * ncomps + c);
        const auto *ip_1_71 = buffer.data(ip_1 + 71 * ncomps + c);
        const auto *ip_1_72 = buffer.data(ip_1 + 72 * ncomps + c);
        const auto *ip_1_73 = buffer.data(ip_1 + 73 * ncomps + c);
        const auto *ip_1_74 = buffer.data(ip_1 + 74 * ncomps + c);
        const auto *ip_1_75 = buffer.data(ip_1 + 75 * ncomps + c);
        const auto *ip_1_76 = buffer.data(ip_1 + 76 * ncomps + c);
        const auto *ip_1_77 = buffer.data(ip_1 + 77 * ncomps + c);
        const auto *ip_1_78 = buffer.data(ip_1 + 78 * ncomps + c);
        const auto *ip_1_79 = buffer.data(ip_1 + 79 * ncomps + c);
        const auto *ip_1_80 = buffer.data(ip_1 + 80 * ncomps + c);
        const auto *ip_1_81 = buffer.data(ip_1 + 81 * ncomps + c);
        const auto *ip_1_82 = buffer.data(ip_1 + 82 * ncomps + c);
        const auto *ip_1_83 = buffer.data(ip_1 + 83 * ncomps + c);

        const auto *ip_0_60 = buffer.data(ip_0 + 60 * ncomps + c);
        const auto *ip_0_61 = buffer.data(ip_0 + 61 * ncomps + c);
        const auto *ip_0_62 = buffer.data(ip_0 + 62 * ncomps + c);
        const auto *ip_0_63 = buffer.data(ip_0 + 63 * ncomps + c);
        const auto *ip_0_64 = buffer.data(ip_0 + 64 * ncomps + c);
        const auto *ip_0_65 = buffer.data(ip_0 + 65 * ncomps + c);
        const auto *ip_0_66 = buffer.data(ip_0 + 66 * ncomps + c);
        const auto *ip_0_67 = buffer.data(ip_0 + 67 * ncomps + c);
        const auto *ip_0_68 = buffer.data(ip_0 + 68 * ncomps + c);
        const auto *ip_0_69 = buffer.data(ip_0 + 69 * ncomps + c);
        const auto *ip_0_70 = buffer.data(ip_0 + 70 * ncomps + c);
        const auto *ip_0_71 = buffer.data(ip_0 + 71 * ncomps + c);
        const auto *ip_0_72 = buffer.data(ip_0 + 72 * ncomps + c);
        const auto *ip_0_73 = buffer.data(ip_0 + 73 * ncomps + c);
        const auto *ip_0_74 = buffer.data(ip_0 + 74 * ncomps + c);
        const auto *ip_0_75 = buffer.data(ip_0 + 75 * ncomps + c);
        const auto *ip_0_76 = buffer.data(ip_0 + 76 * ncomps + c);
        const auto *ip_0_77 = buffer.data(ip_0 + 77 * ncomps + c);
        const auto *ip_0_78 = buffer.data(ip_0 + 78 * ncomps + c);
        const auto *ip_0_79 = buffer.data(ip_0 + 79 * ncomps + c);
        const auto *ip_0_80 = buffer.data(ip_0 + 80 * ncomps + c);
        const auto *ip_0_81 = buffer.data(ip_0 + 81 * ncomps + c);
        const auto *ip_0_82 = buffer.data(ip_0 + 82 * ncomps + c);
        const auto *ip_0_83 = buffer.data(ip_0 + 83 * ncomps + c);

        const auto *kp_1_60 = buffer.data(kp_1 + 60 * ncomps + c);
        const auto *kp_1_61 = buffer.data(kp_1 + 61 * ncomps + c);
        const auto *kp_1_62 = buffer.data(kp_1 + 62 * ncomps + c);
        const auto *kp_1_63 = buffer.data(kp_1 + 63 * ncomps + c);
        const auto *kp_1_64 = buffer.data(kp_1 + 64 * ncomps + c);
        const auto *kp_1_65 = buffer.data(kp_1 + 65 * ncomps + c);
        const auto *kp_1_66 = buffer.data(kp_1 + 66 * ncomps + c);
        const auto *kp_1_67 = buffer.data(kp_1 + 67 * ncomps + c);
        const auto *kp_1_68 = buffer.data(kp_1 + 68 * ncomps + c);
        const auto *kp_1_69 = buffer.data(kp_1 + 69 * ncomps + c);
        const auto *kp_1_70 = buffer.data(kp_1 + 70 * ncomps + c);
        const auto *kp_1_71 = buffer.data(kp_1 + 71 * ncomps + c);
        const auto *kp_1_72 = buffer.data(kp_1 + 72 * ncomps + c);
        const auto *kp_1_73 = buffer.data(kp_1 + 73 * ncomps + c);
        const auto *kp_1_74 = buffer.data(kp_1 + 74 * ncomps + c);
        const auto *kp_1_75 = buffer.data(kp_1 + 75 * ncomps + c);
        const auto *kp_1_76 = buffer.data(kp_1 + 76 * ncomps + c);
        const auto *kp_1_77 = buffer.data(kp_1 + 77 * ncomps + c);
        const auto *kp_1_78 = buffer.data(kp_1 + 78 * ncomps + c);
        const auto *kp_1_79 = buffer.data(kp_1 + 79 * ncomps + c);
        const auto *kp_1_80 = buffer.data(kp_1 + 80 * ncomps + c);
        const auto *kp_1_81 = buffer.data(kp_1 + 81 * ncomps + c);
        const auto *kp_1_82 = buffer.data(kp_1 + 82 * ncomps + c);
        const auto *kp_1_83 = buffer.data(kp_1 + 83 * ncomps + c);
        const auto *kp_1_85 = buffer.data(kp_1 + 85 * ncomps + c);
        const auto *kp_1_86 = buffer.data(kp_1 + 86 * ncomps + c);
        const auto *kp_1_88 = buffer.data(kp_1 + 88 * ncomps + c);
        const auto *kp_1_89 = buffer.data(kp_1 + 89 * ncomps + c);
        const auto *kp_1_91 = buffer.data(kp_1 + 91 * ncomps + c);
        const auto *kp_1_92 = buffer.data(kp_1 + 92 * ncomps + c);
        const auto *kp_1_94 = buffer.data(kp_1 + 94 * ncomps + c);
        const auto *kp_1_95 = buffer.data(kp_1 + 95 * ncomps + c);
        const auto *kp_1_97 = buffer.data(kp_1 + 97 * ncomps + c);
        const auto *kp_1_98 = buffer.data(kp_1 + 98 * ncomps + c);
        const auto *kp_1_100 = buffer.data(kp_1 + 100 * ncomps + c);
        const auto *kp_1_101 = buffer.data(kp_1 + 101 * ncomps + c);
        const auto *kp_1_103 = buffer.data(kp_1 + 103 * ncomps + c);
        const auto *kp_1_104 = buffer.data(kp_1 + 104 * ncomps + c);
        const auto *kp_1_107 = buffer.data(kp_1 + 107 * ncomps + c);

#pragma omp simd aligned(t_120, t_121, t_122, t_123, ab_x, ab_y, ip_1_60, ip_1_61, ip_1_62, \
                         ip_0_60, ip_0_61, ip_0_62, kp_1_60, kp_1_61, kp_1_62, \
                         kp_1_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = ab_x[k] * ip_1_60[k]
                       + ip_0_60[k]
                       + kp_1_60[k];

            t_121[k] = ab_x[k] * ip_1_61[k]
                       + ip_0_61[k]
                       + kp_1_61[k];

            t_122[k] = ab_x[k] * ip_1_62[k]
                       + ip_0_62[k]
                       + kp_1_62[k];

            t_123[k] = ab_y[k] * ip_1_61[k]
                       + kp_1_79[k];
        }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, ab_x, ab_y, ab_z, ip_1_62, ip_1_63, \
                         ip_1_64, ip_0_63, ip_0_64, kp_1_63, kp_1_64, kp_1_80, \
                         kp_1_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_124[k] = ab_y[k] * ip_1_62[k]
                       + kp_1_80[k];

            t_125[k] = ab_z[k] * ip_1_62[k]
                       + kp_1_83[k];

            t_126[k] = ab_x[k] * ip_1_63[k]
                       + ip_0_63[k]
                       + kp_1_63[k];

            t_127[k] = ab_x[k] * ip_1_64[k]
                       + ip_0_64[k]
                       + kp_1_64[k];
        }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, ab_x, ab_y, ab_z, ip_1_64, ip_1_65, \
                         ip_0_65, kp_1_65, kp_1_85, kp_1_86, kp_1_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_128[k] = ab_x[k] * ip_1_65[k]
                       + ip_0_65[k]
                       + kp_1_65[k];

            t_129[k] = ab_y[k] * ip_1_64[k]
                       + kp_1_85[k];

            t_130[k] = ab_y[k] * ip_1_65[k]
                       + kp_1_86[k];

            t_131[k] = ab_z[k] * ip_1_65[k]
                       + kp_1_89[k];
        }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, ab_x, ab_y, ip_1_66, ip_1_67, ip_1_68, \
                         ip_0_66, ip_0_67, ip_0_68, kp_1_66, kp_1_67, kp_1_68, \
                         kp_1_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_132[k] = ab_x[k] * ip_1_66[k]
                       + ip_0_66[k]
                       + kp_1_66[k];

            t_133[k] = ab_x[k] * ip_1_67[k]
                       + ip_0_67[k]
                       + kp_1_67[k];

            t_134[k] = ab_x[k] * ip_1_68[k]
                       + ip_0_68[k]
                       + kp_1_68[k];

            t_135[k] = ab_y[k] * ip_1_67[k]
                       + kp_1_88[k];
        }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, ab_x, ab_y, ab_z, ip_1_68, ip_1_69, \
                         ip_1_70, ip_0_69, ip_0_70, kp_1_69, kp_1_70, kp_1_89, \
                         kp_1_92 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_136[k] = ab_y[k] * ip_1_68[k]
                       + kp_1_89[k];

            t_137[k] = ab_z[k] * ip_1_68[k]
                       + kp_1_92[k];

            t_138[k] = ab_x[k] * ip_1_69[k]
                       + ip_0_69[k]
                       + kp_1_69[k];

            t_139[k] = ab_x[k] * ip_1_70[k]
                       + ip_0_70[k]
                       + kp_1_70[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, ab_x, ab_y, ab_z, ip_1_70, ip_1_71, \
                         ip_0_71, kp_1_71, kp_1_91, kp_1_92, kp_1_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = ab_x[k] * ip_1_71[k]
                       + ip_0_71[k]
                       + kp_1_71[k];

            t_141[k] = ab_y[k] * ip_1_70[k]
                       + kp_1_91[k];

            t_142[k] = ab_y[k] * ip_1_71[k]
                       + kp_1_92[k];

            t_143[k] = ab_z[k] * ip_1_71[k]
                       + kp_1_95[k];
        }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, ab_x, ab_y, ip_1_72, ip_1_73, ip_1_74, \
                         ip_0_72, ip_0_73, ip_0_74, kp_1_72, kp_1_73, kp_1_74, \
                         kp_1_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_144[k] = ab_x[k] * ip_1_72[k]
                       + ip_0_72[k]
                       + kp_1_72[k];

            t_145[k] = ab_x[k] * ip_1_73[k]
                       + ip_0_73[k]
                       + kp_1_73[k];

            t_146[k] = ab_x[k] * ip_1_74[k]
                       + ip_0_74[k]
                       + kp_1_74[k];

            t_147[k] = ab_y[k] * ip_1_73[k]
                       + kp_1_94[k];
        }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, ab_x, ab_y, ab_z, ip_1_74, ip_1_75, \
                         ip_1_76, ip_0_75, ip_0_76, kp_1_75, kp_1_76, kp_1_95, \
                         kp_1_98 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_148[k] = ab_y[k] * ip_1_74[k]
                       + kp_1_95[k];

            t_149[k] = ab_z[k] * ip_1_74[k]
                       + kp_1_98[k];

            t_150[k] = ab_x[k] * ip_1_75[k]
                       + ip_0_75[k]
                       + kp_1_75[k];

            t_151[k] = ab_x[k] * ip_1_76[k]
                       + ip_0_76[k]
                       + kp_1_76[k];
        }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, ab_x, ab_y, ab_z, ip_1_76, ip_1_77, \
                         ip_0_77, kp_1_77, kp_1_97, kp_1_98, kp_1_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_152[k] = ab_x[k] * ip_1_77[k]
                       + ip_0_77[k]
                       + kp_1_77[k];

            t_153[k] = ab_y[k] * ip_1_76[k]
                       + kp_1_97[k];

            t_154[k] = ab_y[k] * ip_1_77[k]
                       + kp_1_98[k];

            t_155[k] = ab_z[k] * ip_1_77[k]
                       + kp_1_101[k];
        }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, ab_x, ab_y, ip_1_78, ip_1_79, ip_1_80, \
                         ip_0_78, ip_0_79, ip_0_80, kp_1_78, kp_1_79, kp_1_80, \
                         kp_1_100 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_156[k] = ab_x[k] * ip_1_78[k]
                       + ip_0_78[k]
                       + kp_1_78[k];

            t_157[k] = ab_x[k] * ip_1_79[k]
                       + ip_0_79[k]
                       + kp_1_79[k];

            t_158[k] = ab_x[k] * ip_1_80[k]
                       + ip_0_80[k]
                       + kp_1_80[k];

            t_159[k] = ab_y[k] * ip_1_79[k]
                       + kp_1_100[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, ab_x, ab_y, ab_z, ip_1_80, ip_1_81, \
                         ip_1_82, ip_0_81, ip_0_82, kp_1_81, kp_1_82, kp_1_101, \
                         kp_1_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_y[k] * ip_1_80[k]
                       + kp_1_101[k];

            t_161[k] = ab_z[k] * ip_1_80[k]
                       + kp_1_104[k];

            t_162[k] = ab_x[k] * ip_1_81[k]
                       + ip_0_81[k]
                       + kp_1_81[k];

            t_163[k] = ab_x[k] * ip_1_82[k]
                       + ip_0_82[k]
                       + kp_1_82[k];
        }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, ab_x, ab_y, ab_z, ip_1_82, ip_1_83, \
                         ip_0_83, kp_1_83, kp_1_103, kp_1_104, \
                         kp_1_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_164[k] = ab_x[k] * ip_1_83[k]
                       + ip_0_83[k]
                       + kp_1_83[k];

            t_165[k] = ab_y[k] * ip_1_82[k]
                       + kp_1_103[k];

            t_166[k] = ab_y[k] * ip_1_83[k]
                       + kp_1_104[k];

            t_167[k] = ab_z[k] * ip_1_83[k]
                       + kp_1_107[k];
        }
    }
}

auto
compute_hrr_geom_100x_id_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                      const size_t target, const size_t ip_1, const size_t ip_0,
                                      const size_t kp_1, const size_t ncomps,
                                      const size_t nmax) -> void
{
    compute_hrr_geom_100x_id_out_of_first_piece0(buffer, coordinates, target, ip_1, ip_0, kp_1,
                                                 ncomps, nmax);

    compute_hrr_geom_100x_id_out_of_first_piece1(buffer, coordinates, target, ip_1, ip_0, kp_1,
                                                 ncomps, nmax);
}

}  // namespace simdtrf
