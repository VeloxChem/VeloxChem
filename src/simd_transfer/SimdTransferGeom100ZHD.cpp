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


#include "SimdTransferGeom100ZHD.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_100z_hd_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                      const size_t target, const size_t hp_1, const size_t hp_0,
                                      const size_t ip_1, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

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
        const auto *hp_1_45 = buffer.data(hp_1 + 45 * ncomps + c);
        const auto *hp_1_46 = buffer.data(hp_1 + 46 * ncomps + c);
        const auto *hp_1_47 = buffer.data(hp_1 + 47 * ncomps + c);
        const auto *hp_1_48 = buffer.data(hp_1 + 48 * ncomps + c);
        const auto *hp_1_49 = buffer.data(hp_1 + 49 * ncomps + c);
        const auto *hp_1_50 = buffer.data(hp_1 + 50 * ncomps + c);
        const auto *hp_1_51 = buffer.data(hp_1 + 51 * ncomps + c);
        const auto *hp_1_52 = buffer.data(hp_1 + 52 * ncomps + c);
        const auto *hp_1_53 = buffer.data(hp_1 + 53 * ncomps + c);
        const auto *hp_1_54 = buffer.data(hp_1 + 54 * ncomps + c);
        const auto *hp_1_55 = buffer.data(hp_1 + 55 * ncomps + c);
        const auto *hp_1_56 = buffer.data(hp_1 + 56 * ncomps + c);
        const auto *hp_1_57 = buffer.data(hp_1 + 57 * ncomps + c);
        const auto *hp_1_58 = buffer.data(hp_1 + 58 * ncomps + c);
        const auto *hp_1_59 = buffer.data(hp_1 + 59 * ncomps + c);
        const auto *hp_1_60 = buffer.data(hp_1 + 60 * ncomps + c);
        const auto *hp_1_61 = buffer.data(hp_1 + 61 * ncomps + c);
        const auto *hp_1_62 = buffer.data(hp_1 + 62 * ncomps + c);

        const auto *hp_0_2 = buffer.data(hp_0 + 2 * ncomps + c);
        const auto *hp_0_5 = buffer.data(hp_0 + 5 * ncomps + c);
        const auto *hp_0_8 = buffer.data(hp_0 + 8 * ncomps + c);
        const auto *hp_0_11 = buffer.data(hp_0 + 11 * ncomps + c);
        const auto *hp_0_14 = buffer.data(hp_0 + 14 * ncomps + c);
        const auto *hp_0_17 = buffer.data(hp_0 + 17 * ncomps + c);
        const auto *hp_0_20 = buffer.data(hp_0 + 20 * ncomps + c);
        const auto *hp_0_23 = buffer.data(hp_0 + 23 * ncomps + c);
        const auto *hp_0_26 = buffer.data(hp_0 + 26 * ncomps + c);
        const auto *hp_0_29 = buffer.data(hp_0 + 29 * ncomps + c);
        const auto *hp_0_32 = buffer.data(hp_0 + 32 * ncomps + c);
        const auto *hp_0_35 = buffer.data(hp_0 + 35 * ncomps + c);
        const auto *hp_0_38 = buffer.data(hp_0 + 38 * ncomps + c);
        const auto *hp_0_41 = buffer.data(hp_0 + 41 * ncomps + c);
        const auto *hp_0_44 = buffer.data(hp_0 + 44 * ncomps + c);
        const auto *hp_0_47 = buffer.data(hp_0 + 47 * ncomps + c);
        const auto *hp_0_50 = buffer.data(hp_0 + 50 * ncomps + c);
        const auto *hp_0_53 = buffer.data(hp_0 + 53 * ncomps + c);
        const auto *hp_0_56 = buffer.data(hp_0 + 56 * ncomps + c);
        const auto *hp_0_59 = buffer.data(hp_0 + 59 * ncomps + c);
        const auto *hp_0_62 = buffer.data(hp_0 + 62 * ncomps + c);

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
        const auto *ip_1_60 = buffer.data(ip_1 + 60 * ncomps + c);
        const auto *ip_1_61 = buffer.data(ip_1 + 61 * ncomps + c);
        const auto *ip_1_62 = buffer.data(ip_1 + 62 * ncomps + c);
        const auto *ip_1_64 = buffer.data(ip_1 + 64 * ncomps + c);
        const auto *ip_1_65 = buffer.data(ip_1 + 65 * ncomps + c);
        const auto *ip_1_67 = buffer.data(ip_1 + 67 * ncomps + c);
        const auto *ip_1_68 = buffer.data(ip_1 + 68 * ncomps + c);
        const auto *ip_1_70 = buffer.data(ip_1 + 70 * ncomps + c);
        const auto *ip_1_71 = buffer.data(ip_1 + 71 * ncomps + c);
        const auto *ip_1_73 = buffer.data(ip_1 + 73 * ncomps + c);
        const auto *ip_1_74 = buffer.data(ip_1 + 74 * ncomps + c);
        const auto *ip_1_76 = buffer.data(ip_1 + 76 * ncomps + c);
        const auto *ip_1_77 = buffer.data(ip_1 + 77 * ncomps + c);
        const auto *ip_1_79 = buffer.data(ip_1 + 79 * ncomps + c);
        const auto *ip_1_80 = buffer.data(ip_1 + 80 * ncomps + c);
        const auto *ip_1_83 = buffer.data(ip_1 + 83 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, hp_1_0, hp_1_1, hp_1_2, ip_1_0, \
                         ip_1_1, ip_1_2, ip_1_4, ip_1_5 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * hp_1_0[k]
                     + ip_1_0[k];

            t_1[k] = ab_x[k] * hp_1_1[k]
                     + ip_1_1[k];

            t_2[k] = ab_x[k] * hp_1_2[k]
                     + ip_1_2[k];

            t_3[k] = ab_y[k] * hp_1_1[k]
                     + ip_1_4[k];

            t_4[k] = ab_y[k] * hp_1_2[k]
                     + ip_1_5[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_z, hp_1_2, hp_1_3, hp_1_4, hp_1_5, \
                         hp_0_2, ip_1_3, ip_1_4, ip_1_5, ip_1_8 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_z[k] * hp_1_2[k]
                     + hp_0_2[k]
                     + ip_1_8[k];

            t_6[k] = ab_x[k] * hp_1_3[k]
                     + ip_1_3[k];

            t_7[k] = ab_x[k] * hp_1_4[k]
                     + ip_1_4[k];

            t_8[k] = ab_x[k] * hp_1_5[k]
                     + ip_1_5[k];
        }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_y, ab_z, hp_1_4, hp_1_5, hp_1_6, \
                         hp_0_5, ip_1_6, ip_1_10, ip_1_11, ip_1_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_9[k] = ab_y[k] * hp_1_4[k]
                     + ip_1_10[k];

            t_10[k] = ab_y[k] * hp_1_5[k]
                      + ip_1_11[k];

            t_11[k] = ab_z[k] * hp_1_5[k]
                      + hp_0_5[k]
                      + ip_1_14[k];

            t_12[k] = ab_x[k] * hp_1_6[k]
                      + ip_1_6[k];
        }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, hp_1_7, hp_1_8, \
                         hp_0_8, ip_1_7, ip_1_8, ip_1_13, ip_1_14, \
                         ip_1_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_13[k] = ab_x[k] * hp_1_7[k]
                      + ip_1_7[k];

            t_14[k] = ab_x[k] * hp_1_8[k]
                      + ip_1_8[k];

            t_15[k] = ab_y[k] * hp_1_7[k]
                      + ip_1_13[k];

            t_16[k] = ab_y[k] * hp_1_8[k]
                      + ip_1_14[k];

            t_17[k] = ab_z[k] * hp_1_8[k]
                      + hp_0_8[k]
                      + ip_1_17[k];
        }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, hp_1_9, hp_1_10, hp_1_11, \
                         ip_1_9, ip_1_10, ip_1_11, ip_1_19, ip_1_20 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_18[k] = ab_x[k] * hp_1_9[k]
                      + ip_1_9[k];

            t_19[k] = ab_x[k] * hp_1_10[k]
                      + ip_1_10[k];

            t_20[k] = ab_x[k] * hp_1_11[k]
                      + ip_1_11[k];

            t_21[k] = ab_y[k] * hp_1_10[k]
                      + ip_1_19[k];

            t_22[k] = ab_y[k] * hp_1_11[k]
                      + ip_1_20[k];
        }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_z, hp_1_11, hp_1_12, hp_1_13, \
                         hp_1_14, hp_0_11, ip_1_12, ip_1_13, ip_1_14, \
                         ip_1_23 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_23[k] = ab_z[k] * hp_1_11[k]
                      + hp_0_11[k]
                      + ip_1_23[k];

            t_24[k] = ab_x[k] * hp_1_12[k]
                      + ip_1_12[k];

            t_25[k] = ab_x[k] * hp_1_13[k]
                      + ip_1_13[k];

            t_26[k] = ab_x[k] * hp_1_14[k]
                      + ip_1_14[k];
        }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, hp_1_13, hp_1_14, hp_1_15, \
                         hp_0_14, ip_1_15, ip_1_22, ip_1_23, ip_1_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_27[k] = ab_y[k] * hp_1_13[k]
                      + ip_1_22[k];

            t_28[k] = ab_y[k] * hp_1_14[k]
                      + ip_1_23[k];

            t_29[k] = ab_z[k] * hp_1_14[k]
                      + hp_0_14[k]
                      + ip_1_26[k];

            t_30[k] = ab_x[k] * hp_1_15[k]
                      + ip_1_15[k];
        }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, hp_1_16, hp_1_17, \
                         hp_0_17, ip_1_16, ip_1_17, ip_1_25, ip_1_26, \
                         ip_1_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_31[k] = ab_x[k] * hp_1_16[k]
                      + ip_1_16[k];

            t_32[k] = ab_x[k] * hp_1_17[k]
                      + ip_1_17[k];

            t_33[k] = ab_y[k] * hp_1_16[k]
                      + ip_1_25[k];

            t_34[k] = ab_y[k] * hp_1_17[k]
                      + ip_1_26[k];

            t_35[k] = ab_z[k] * hp_1_17[k]
                      + hp_0_17[k]
                      + ip_1_29[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_x, ab_y, hp_1_18, hp_1_19, hp_1_20, \
                         ip_1_18, ip_1_19, ip_1_20, ip_1_31, ip_1_32 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = ab_x[k] * hp_1_18[k]
                      + ip_1_18[k];

            t_37[k] = ab_x[k] * hp_1_19[k]
                      + ip_1_19[k];

            t_38[k] = ab_x[k] * hp_1_20[k]
                      + ip_1_20[k];

            t_39[k] = ab_y[k] * hp_1_19[k]
                      + ip_1_31[k];

            t_40[k] = ab_y[k] * hp_1_20[k]
                      + ip_1_32[k];
        }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_x, ab_z, hp_1_20, hp_1_21, hp_1_22, \
                         hp_1_23, hp_0_20, ip_1_21, ip_1_22, ip_1_23, \
                         ip_1_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_41[k] = ab_z[k] * hp_1_20[k]
                      + hp_0_20[k]
                      + ip_1_35[k];

            t_42[k] = ab_x[k] * hp_1_21[k]
                      + ip_1_21[k];

            t_43[k] = ab_x[k] * hp_1_22[k]
                      + ip_1_22[k];

            t_44[k] = ab_x[k] * hp_1_23[k]
                      + ip_1_23[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, hp_1_22, hp_1_23, hp_1_24, \
                         hp_0_23, ip_1_24, ip_1_34, ip_1_35, ip_1_38 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_y[k] * hp_1_22[k]
                      + ip_1_34[k];

            t_46[k] = ab_y[k] * hp_1_23[k]
                      + ip_1_35[k];

            t_47[k] = ab_z[k] * hp_1_23[k]
                      + hp_0_23[k]
                      + ip_1_38[k];

            t_48[k] = ab_x[k] * hp_1_24[k]
                      + ip_1_24[k];
        }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, hp_1_25, hp_1_26, \
                         hp_0_26, ip_1_25, ip_1_26, ip_1_37, ip_1_38, \
                         ip_1_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_49[k] = ab_x[k] * hp_1_25[k]
                      + ip_1_25[k];

            t_50[k] = ab_x[k] * hp_1_26[k]
                      + ip_1_26[k];

            t_51[k] = ab_y[k] * hp_1_25[k]
                      + ip_1_37[k];

            t_52[k] = ab_y[k] * hp_1_26[k]
                      + ip_1_38[k];

            t_53[k] = ab_z[k] * hp_1_26[k]
                      + hp_0_26[k]
                      + ip_1_41[k];
        }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_x, ab_y, hp_1_27, hp_1_28, hp_1_29, \
                         ip_1_27, ip_1_28, ip_1_29, ip_1_40, ip_1_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_54[k] = ab_x[k] * hp_1_27[k]
                      + ip_1_27[k];

            t_55[k] = ab_x[k] * hp_1_28[k]
                      + ip_1_28[k];

            t_56[k] = ab_x[k] * hp_1_29[k]
                      + ip_1_29[k];

            t_57[k] = ab_y[k] * hp_1_28[k]
                      + ip_1_40[k];

            t_58[k] = ab_y[k] * hp_1_29[k]
                      + ip_1_41[k];
        }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, ab_x, ab_z, hp_1_29, hp_1_30, hp_1_31, \
                         hp_1_32, hp_0_29, ip_1_30, ip_1_31, ip_1_32, \
                         ip_1_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_59[k] = ab_z[k] * hp_1_29[k]
                      + hp_0_29[k]
                      + ip_1_44[k];

            t_60[k] = ab_x[k] * hp_1_30[k]
                      + ip_1_30[k];

            t_61[k] = ab_x[k] * hp_1_31[k]
                      + ip_1_31[k];

            t_62[k] = ab_x[k] * hp_1_32[k]
                      + ip_1_32[k];
        }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, ab_x, ab_y, ab_z, hp_1_31, hp_1_32, hp_1_33, \
                         hp_0_32, ip_1_33, ip_1_46, ip_1_47, ip_1_50 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_63[k] = ab_y[k] * hp_1_31[k]
                      + ip_1_46[k];

            t_64[k] = ab_y[k] * hp_1_32[k]
                      + ip_1_47[k];

            t_65[k] = ab_z[k] * hp_1_32[k]
                      + hp_0_32[k]
                      + ip_1_50[k];

            t_66[k] = ab_x[k] * hp_1_33[k]
                      + ip_1_33[k];
        }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, ab_x, ab_y, ab_z, hp_1_34, hp_1_35, \
                         hp_0_35, ip_1_34, ip_1_35, ip_1_49, ip_1_50, \
                         ip_1_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_67[k] = ab_x[k] * hp_1_34[k]
                      + ip_1_34[k];

            t_68[k] = ab_x[k] * hp_1_35[k]
                      + ip_1_35[k];

            t_69[k] = ab_y[k] * hp_1_34[k]
                      + ip_1_49[k];

            t_70[k] = ab_y[k] * hp_1_35[k]
                      + ip_1_50[k];

            t_71[k] = ab_z[k] * hp_1_35[k]
                      + hp_0_35[k]
                      + ip_1_53[k];
        }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, hp_1_36, hp_1_37, hp_1_38, \
                         ip_1_36, ip_1_37, ip_1_38, ip_1_52, ip_1_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_72[k] = ab_x[k] * hp_1_36[k]
                      + ip_1_36[k];

            t_73[k] = ab_x[k] * hp_1_37[k]
                      + ip_1_37[k];

            t_74[k] = ab_x[k] * hp_1_38[k]
                      + ip_1_38[k];

            t_75[k] = ab_y[k] * hp_1_37[k]
                      + ip_1_52[k];

            t_76[k] = ab_y[k] * hp_1_38[k]
                      + ip_1_53[k];
        }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, ab_x, ab_z, hp_1_38, hp_1_39, hp_1_40, \
                         hp_1_41, hp_0_38, ip_1_39, ip_1_40, ip_1_41, \
                         ip_1_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_77[k] = ab_z[k] * hp_1_38[k]
                      + hp_0_38[k]
                      + ip_1_56[k];

            t_78[k] = ab_x[k] * hp_1_39[k]
                      + ip_1_39[k];

            t_79[k] = ab_x[k] * hp_1_40[k]
                      + ip_1_40[k];

            t_80[k] = ab_x[k] * hp_1_41[k]
                      + ip_1_41[k];
        }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, ab_x, ab_y, ab_z, hp_1_40, hp_1_41, hp_1_42, \
                         hp_0_41, ip_1_42, ip_1_55, ip_1_56, ip_1_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_81[k] = ab_y[k] * hp_1_40[k]
                      + ip_1_55[k];

            t_82[k] = ab_y[k] * hp_1_41[k]
                      + ip_1_56[k];

            t_83[k] = ab_z[k] * hp_1_41[k]
                      + hp_0_41[k]
                      + ip_1_59[k];

            t_84[k] = ab_x[k] * hp_1_42[k]
                      + ip_1_42[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, hp_1_43, hp_1_44, \
                         hp_0_44, ip_1_43, ip_1_44, ip_1_58, ip_1_59, \
                         ip_1_62 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_x[k] * hp_1_43[k]
                      + ip_1_43[k];

            t_86[k] = ab_x[k] * hp_1_44[k]
                      + ip_1_44[k];

            t_87[k] = ab_y[k] * hp_1_43[k]
                      + ip_1_58[k];

            t_88[k] = ab_y[k] * hp_1_44[k]
                      + ip_1_59[k];

            t_89[k] = ab_z[k] * hp_1_44[k]
                      + hp_0_44[k]
                      + ip_1_62[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, ab_y, hp_1_45, hp_1_46, hp_1_47, \
                         ip_1_45, ip_1_46, ip_1_47, ip_1_64, ip_1_65 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * hp_1_45[k]
                      + ip_1_45[k];

            t_91[k] = ab_x[k] * hp_1_46[k]
                      + ip_1_46[k];

            t_92[k] = ab_x[k] * hp_1_47[k]
                      + ip_1_47[k];

            t_93[k] = ab_y[k] * hp_1_46[k]
                      + ip_1_64[k];

            t_94[k] = ab_y[k] * hp_1_47[k]
                      + ip_1_65[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, ab_x, ab_z, hp_1_47, hp_1_48, hp_1_49, \
                         hp_1_50, hp_0_47, ip_1_48, ip_1_49, ip_1_50, \
                         ip_1_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_z[k] * hp_1_47[k]
                      + hp_0_47[k]
                      + ip_1_68[k];

            t_96[k] = ab_x[k] * hp_1_48[k]
                      + ip_1_48[k];

            t_97[k] = ab_x[k] * hp_1_49[k]
                      + ip_1_49[k];

            t_98[k] = ab_x[k] * hp_1_50[k]
                      + ip_1_50[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, ab_x, ab_y, ab_z, hp_1_49, hp_1_50, \
                         hp_1_51, hp_0_50, ip_1_51, ip_1_67, ip_1_68, \
                         ip_1_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_y[k] * hp_1_49[k]
                      + ip_1_67[k];

            t_100[k] = ab_y[k] * hp_1_50[k]
                       + ip_1_68[k];

            t_101[k] = ab_z[k] * hp_1_50[k]
                       + hp_0_50[k]
                       + ip_1_71[k];

            t_102[k] = ab_x[k] * hp_1_51[k]
                       + ip_1_51[k];
        }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, ab_x, ab_y, ab_z, hp_1_52, \
                         hp_1_53, hp_0_53, ip_1_52, ip_1_53, ip_1_70, ip_1_71, \
                         ip_1_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_103[k] = ab_x[k] * hp_1_52[k]
                       + ip_1_52[k];

            t_104[k] = ab_x[k] * hp_1_53[k]
                       + ip_1_53[k];

            t_105[k] = ab_y[k] * hp_1_52[k]
                       + ip_1_70[k];

            t_106[k] = ab_y[k] * hp_1_53[k]
                       + ip_1_71[k];

            t_107[k] = ab_z[k] * hp_1_53[k]
                       + hp_0_53[k]
                       + ip_1_74[k];
        }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, ab_x, ab_y, hp_1_54, hp_1_55, \
                         hp_1_56, ip_1_54, ip_1_55, ip_1_56, ip_1_73, \
                         ip_1_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_108[k] = ab_x[k] * hp_1_54[k]
                       + ip_1_54[k];

            t_109[k] = ab_x[k] * hp_1_55[k]
                       + ip_1_55[k];

            t_110[k] = ab_x[k] * hp_1_56[k]
                       + ip_1_56[k];

            t_111[k] = ab_y[k] * hp_1_55[k]
                       + ip_1_73[k];

            t_112[k] = ab_y[k] * hp_1_56[k]
                       + ip_1_74[k];
        }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, ab_x, ab_z, hp_1_56, hp_1_57, hp_1_58, \
                         hp_1_59, hp_0_56, ip_1_57, ip_1_58, ip_1_59, \
                         ip_1_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_113[k] = ab_z[k] * hp_1_56[k]
                       + hp_0_56[k]
                       + ip_1_77[k];

            t_114[k] = ab_x[k] * hp_1_57[k]
                       + ip_1_57[k];

            t_115[k] = ab_x[k] * hp_1_58[k]
                       + ip_1_58[k];

            t_116[k] = ab_x[k] * hp_1_59[k]
                       + ip_1_59[k];
        }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, ab_x, ab_y, ab_z, hp_1_58, hp_1_59, \
                         hp_1_60, hp_0_59, ip_1_60, ip_1_76, ip_1_77, \
                         ip_1_80 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_117[k] = ab_y[k] * hp_1_58[k]
                       + ip_1_76[k];

            t_118[k] = ab_y[k] * hp_1_59[k]
                       + ip_1_77[k];

            t_119[k] = ab_z[k] * hp_1_59[k]
                       + hp_0_59[k]
                       + ip_1_80[k];

            t_120[k] = ab_x[k] * hp_1_60[k]
                       + ip_1_60[k];
        }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, ab_x, ab_y, ab_z, hp_1_61, \
                         hp_1_62, hp_0_62, ip_1_61, ip_1_62, ip_1_79, ip_1_80, \
                         ip_1_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_121[k] = ab_x[k] * hp_1_61[k]
                       + ip_1_61[k];

            t_122[k] = ab_x[k] * hp_1_62[k]
                       + ip_1_62[k];

            t_123[k] = ab_y[k] * hp_1_61[k]
                       + ip_1_79[k];

            t_124[k] = ab_y[k] * hp_1_62[k]
                       + ip_1_80[k];

            t_125[k] = ab_z[k] * hp_1_62[k]
                       + hp_0_62[k]
                       + ip_1_83[k];
        }
    }
}

}  // namespace simdtrf
