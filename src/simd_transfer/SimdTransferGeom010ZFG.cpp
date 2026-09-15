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


#include "SimdTransferGeom010ZFG.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_geom_010z_fg_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                const size_t target, const size_t dg_1, const size_t dg_0,
                                const size_t dh_1, const size_t ncomps,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

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
        const auto *dg_1_10 = buffer.data(dg_1 + 10 * ncomps + c);
        const auto *dg_1_11 = buffer.data(dg_1 + 11 * ncomps + c);
        const auto *dg_1_12 = buffer.data(dg_1 + 12 * ncomps + c);
        const auto *dg_1_13 = buffer.data(dg_1 + 13 * ncomps + c);
        const auto *dg_1_14 = buffer.data(dg_1 + 14 * ncomps + c);
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
        const auto *dg_1_25 = buffer.data(dg_1 + 25 * ncomps + c);
        const auto *dg_1_26 = buffer.data(dg_1 + 26 * ncomps + c);
        const auto *dg_1_27 = buffer.data(dg_1 + 27 * ncomps + c);
        const auto *dg_1_28 = buffer.data(dg_1 + 28 * ncomps + c);
        const auto *dg_1_29 = buffer.data(dg_1 + 29 * ncomps + c);
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
        const auto *dg_1_40 = buffer.data(dg_1 + 40 * ncomps + c);
        const auto *dg_1_41 = buffer.data(dg_1 + 41 * ncomps + c);
        const auto *dg_1_42 = buffer.data(dg_1 + 42 * ncomps + c);
        const auto *dg_1_43 = buffer.data(dg_1 + 43 * ncomps + c);
        const auto *dg_1_44 = buffer.data(dg_1 + 44 * ncomps + c);
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
        const auto *dg_1_59 = buffer.data(dg_1 + 59 * ncomps + c);
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
        const auto *dg_1_74 = buffer.data(dg_1 + 74 * ncomps + c);
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

        const auto *dg_0_75 = buffer.data(dg_0 + 75 * ncomps + c);
        const auto *dg_0_76 = buffer.data(dg_0 + 76 * ncomps + c);
        const auto *dg_0_77 = buffer.data(dg_0 + 77 * ncomps + c);
        const auto *dg_0_78 = buffer.data(dg_0 + 78 * ncomps + c);
        const auto *dg_0_79 = buffer.data(dg_0 + 79 * ncomps + c);
        const auto *dg_0_80 = buffer.data(dg_0 + 80 * ncomps + c);

        const auto *dh_1_0 = buffer.data(dh_1 + 0 * ncomps + c);
        const auto *dh_1_1 = buffer.data(dh_1 + 1 * ncomps + c);
        const auto *dh_1_2 = buffer.data(dh_1 + 2 * ncomps + c);
        const auto *dh_1_3 = buffer.data(dh_1 + 3 * ncomps + c);
        const auto *dh_1_4 = buffer.data(dh_1 + 4 * ncomps + c);
        const auto *dh_1_5 = buffer.data(dh_1 + 5 * ncomps + c);
        const auto *dh_1_6 = buffer.data(dh_1 + 6 * ncomps + c);
        const auto *dh_1_7 = buffer.data(dh_1 + 7 * ncomps + c);
        const auto *dh_1_8 = buffer.data(dh_1 + 8 * ncomps + c);
        const auto *dh_1_9 = buffer.data(dh_1 + 9 * ncomps + c);
        const auto *dh_1_10 = buffer.data(dh_1 + 10 * ncomps + c);
        const auto *dh_1_11 = buffer.data(dh_1 + 11 * ncomps + c);
        const auto *dh_1_12 = buffer.data(dh_1 + 12 * ncomps + c);
        const auto *dh_1_13 = buffer.data(dh_1 + 13 * ncomps + c);
        const auto *dh_1_14 = buffer.data(dh_1 + 14 * ncomps + c);
        const auto *dh_1_21 = buffer.data(dh_1 + 21 * ncomps + c);
        const auto *dh_1_22 = buffer.data(dh_1 + 22 * ncomps + c);
        const auto *dh_1_23 = buffer.data(dh_1 + 23 * ncomps + c);
        const auto *dh_1_24 = buffer.data(dh_1 + 24 * ncomps + c);
        const auto *dh_1_25 = buffer.data(dh_1 + 25 * ncomps + c);
        const auto *dh_1_26 = buffer.data(dh_1 + 26 * ncomps + c);
        const auto *dh_1_27 = buffer.data(dh_1 + 27 * ncomps + c);
        const auto *dh_1_28 = buffer.data(dh_1 + 28 * ncomps + c);
        const auto *dh_1_29 = buffer.data(dh_1 + 29 * ncomps + c);
        const auto *dh_1_30 = buffer.data(dh_1 + 30 * ncomps + c);
        const auto *dh_1_31 = buffer.data(dh_1 + 31 * ncomps + c);
        const auto *dh_1_32 = buffer.data(dh_1 + 32 * ncomps + c);
        const auto *dh_1_33 = buffer.data(dh_1 + 33 * ncomps + c);
        const auto *dh_1_34 = buffer.data(dh_1 + 34 * ncomps + c);
        const auto *dh_1_35 = buffer.data(dh_1 + 35 * ncomps + c);
        const auto *dh_1_42 = buffer.data(dh_1 + 42 * ncomps + c);
        const auto *dh_1_43 = buffer.data(dh_1 + 43 * ncomps + c);
        const auto *dh_1_44 = buffer.data(dh_1 + 44 * ncomps + c);
        const auto *dh_1_45 = buffer.data(dh_1 + 45 * ncomps + c);
        const auto *dh_1_46 = buffer.data(dh_1 + 46 * ncomps + c);
        const auto *dh_1_47 = buffer.data(dh_1 + 47 * ncomps + c);
        const auto *dh_1_48 = buffer.data(dh_1 + 48 * ncomps + c);
        const auto *dh_1_49 = buffer.data(dh_1 + 49 * ncomps + c);
        const auto *dh_1_50 = buffer.data(dh_1 + 50 * ncomps + c);
        const auto *dh_1_51 = buffer.data(dh_1 + 51 * ncomps + c);
        const auto *dh_1_52 = buffer.data(dh_1 + 52 * ncomps + c);
        const auto *dh_1_53 = buffer.data(dh_1 + 53 * ncomps + c);
        const auto *dh_1_54 = buffer.data(dh_1 + 54 * ncomps + c);
        const auto *dh_1_55 = buffer.data(dh_1 + 55 * ncomps + c);
        const auto *dh_1_56 = buffer.data(dh_1 + 56 * ncomps + c);
        const auto *dh_1_63 = buffer.data(dh_1 + 63 * ncomps + c);
        const auto *dh_1_64 = buffer.data(dh_1 + 64 * ncomps + c);
        const auto *dh_1_65 = buffer.data(dh_1 + 65 * ncomps + c);
        const auto *dh_1_66 = buffer.data(dh_1 + 66 * ncomps + c);
        const auto *dh_1_67 = buffer.data(dh_1 + 67 * ncomps + c);
        const auto *dh_1_68 = buffer.data(dh_1 + 68 * ncomps + c);
        const auto *dh_1_69 = buffer.data(dh_1 + 69 * ncomps + c);
        const auto *dh_1_70 = buffer.data(dh_1 + 70 * ncomps + c);
        const auto *dh_1_71 = buffer.data(dh_1 + 71 * ncomps + c);
        const auto *dh_1_72 = buffer.data(dh_1 + 72 * ncomps + c);
        const auto *dh_1_73 = buffer.data(dh_1 + 73 * ncomps + c);
        const auto *dh_1_74 = buffer.data(dh_1 + 74 * ncomps + c);
        const auto *dh_1_75 = buffer.data(dh_1 + 75 * ncomps + c);
        const auto *dh_1_76 = buffer.data(dh_1 + 76 * ncomps + c);
        const auto *dh_1_77 = buffer.data(dh_1 + 77 * ncomps + c);
        const auto *dh_1_78 = buffer.data(dh_1 + 78 * ncomps + c);
        const auto *dh_1_79 = buffer.data(dh_1 + 79 * ncomps + c);
        const auto *dh_1_80 = buffer.data(dh_1 + 80 * ncomps + c);
        const auto *dh_1_81 = buffer.data(dh_1 + 81 * ncomps + c);
        const auto *dh_1_82 = buffer.data(dh_1 + 82 * ncomps + c);
        const auto *dh_1_84 = buffer.data(dh_1 + 84 * ncomps + c);
        const auto *dh_1_85 = buffer.data(dh_1 + 85 * ncomps + c);
        const auto *dh_1_86 = buffer.data(dh_1 + 86 * ncomps + c);
        const auto *dh_1_87 = buffer.data(dh_1 + 87 * ncomps + c);
        const auto *dh_1_88 = buffer.data(dh_1 + 88 * ncomps + c);
        const auto *dh_1_89 = buffer.data(dh_1 + 89 * ncomps + c);
        const auto *dh_1_90 = buffer.data(dh_1 + 90 * ncomps + c);
        const auto *dh_1_91 = buffer.data(dh_1 + 91 * ncomps + c);
        const auto *dh_1_92 = buffer.data(dh_1 + 92 * ncomps + c);
        const auto *dh_1_93 = buffer.data(dh_1 + 93 * ncomps + c);
        const auto *dh_1_94 = buffer.data(dh_1 + 94 * ncomps + c);
        const auto *dh_1_95 = buffer.data(dh_1 + 95 * ncomps + c);
        const auto *dh_1_96 = buffer.data(dh_1 + 96 * ncomps + c);
        const auto *dh_1_97 = buffer.data(dh_1 + 97 * ncomps + c);
        const auto *dh_1_98 = buffer.data(dh_1 + 98 * ncomps + c);
        const auto *dh_1_99 = buffer.data(dh_1 + 99 * ncomps + c);
        const auto *dh_1_100 = buffer.data(dh_1 + 100 * ncomps + c);
        const auto *dh_1_101 = buffer.data(dh_1 + 101 * ncomps + c);
        const auto *dh_1_102 = buffer.data(dh_1 + 102 * ncomps + c);
        const auto *dh_1_103 = buffer.data(dh_1 + 103 * ncomps + c);
        const auto *dh_1_105 = buffer.data(dh_1 + 105 * ncomps + c);
        const auto *dh_1_106 = buffer.data(dh_1 + 106 * ncomps + c);
        const auto *dh_1_107 = buffer.data(dh_1 + 107 * ncomps + c);
        const auto *dh_1_108 = buffer.data(dh_1 + 108 * ncomps + c);
        const auto *dh_1_109 = buffer.data(dh_1 + 109 * ncomps + c);
        const auto *dh_1_110 = buffer.data(dh_1 + 110 * ncomps + c);
        const auto *dh_1_111 = buffer.data(dh_1 + 111 * ncomps + c);
        const auto *dh_1_112 = buffer.data(dh_1 + 112 * ncomps + c);
        const auto *dh_1_113 = buffer.data(dh_1 + 113 * ncomps + c);
        const auto *dh_1_114 = buffer.data(dh_1 + 114 * ncomps + c);
        const auto *dh_1_115 = buffer.data(dh_1 + 115 * ncomps + c);
        const auto *dh_1_116 = buffer.data(dh_1 + 116 * ncomps + c);
        const auto *dh_1_117 = buffer.data(dh_1 + 117 * ncomps + c);
        const auto *dh_1_118 = buffer.data(dh_1 + 118 * ncomps + c);
        const auto *dh_1_119 = buffer.data(dh_1 + 119 * ncomps + c);
        const auto *dh_1_120 = buffer.data(dh_1 + 120 * ncomps + c);
        const auto *dh_1_121 = buffer.data(dh_1 + 121 * ncomps + c);
        const auto *dh_1_122 = buffer.data(dh_1 + 122 * ncomps + c);
        const auto *dh_1_123 = buffer.data(dh_1 + 123 * ncomps + c);
        const auto *dh_1_124 = buffer.data(dh_1 + 124 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, dg_1_0, dg_1_1, dg_1_2, dg_1_3, \
                         dg_1_4, dh_1_0, dh_1_1, dh_1_2, dh_1_3, \
                         dh_1_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * dg_1_0[k]
                     + dh_1_0[k];

            t_1[k] = -ab_x[k] * dg_1_1[k]
                     + dh_1_1[k];

            t_2[k] = -ab_x[k] * dg_1_2[k]
                     + dh_1_2[k];

            t_3[k] = -ab_x[k] * dg_1_3[k]
                     + dh_1_3[k];

            t_4[k] = -ab_x[k] * dg_1_4[k]
                     + dh_1_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, dg_1_5, dg_1_6, dg_1_7, dg_1_8, \
                         dg_1_9, dh_1_5, dh_1_6, dh_1_7, dh_1_8, \
                         dh_1_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * dg_1_5[k]
                     + dh_1_5[k];

            t_6[k] = -ab_x[k] * dg_1_6[k]
                     + dh_1_6[k];

            t_7[k] = -ab_x[k] * dg_1_7[k]
                     + dh_1_7[k];

            t_8[k] = -ab_x[k] * dg_1_8[k]
                     + dh_1_8[k];

            t_9[k] = -ab_x[k] * dg_1_9[k]
                     + dh_1_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, dg_1_10, dg_1_11, dg_1_12, \
                         dg_1_13, dg_1_14, dh_1_10, dh_1_11, dh_1_12, dh_1_13, \
                         dh_1_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * dg_1_10[k]
                      + dh_1_10[k];

            t_11[k] = -ab_x[k] * dg_1_11[k]
                      + dh_1_11[k];

            t_12[k] = -ab_x[k] * dg_1_12[k]
                      + dh_1_12[k];

            t_13[k] = -ab_x[k] * dg_1_13[k]
                      + dh_1_13[k];

            t_14[k] = -ab_x[k] * dg_1_14[k]
                      + dh_1_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, dg_1_15, dg_1_16, dg_1_17, \
                         dg_1_18, dg_1_19, dh_1_21, dh_1_22, dh_1_23, dh_1_24, \
                         dh_1_25 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * dg_1_15[k]
                      + dh_1_21[k];

            t_16[k] = -ab_x[k] * dg_1_16[k]
                      + dh_1_22[k];

            t_17[k] = -ab_x[k] * dg_1_17[k]
                      + dh_1_23[k];

            t_18[k] = -ab_x[k] * dg_1_18[k]
                      + dh_1_24[k];

            t_19[k] = -ab_x[k] * dg_1_19[k]
                      + dh_1_25[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, dg_1_20, dg_1_21, dg_1_22, \
                         dg_1_23, dg_1_24, dh_1_26, dh_1_27, dh_1_28, dh_1_29, \
                         dh_1_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * dg_1_20[k]
                      + dh_1_26[k];

            t_21[k] = -ab_x[k] * dg_1_21[k]
                      + dh_1_27[k];

            t_22[k] = -ab_x[k] * dg_1_22[k]
                      + dh_1_28[k];

            t_23[k] = -ab_x[k] * dg_1_23[k]
                      + dh_1_29[k];

            t_24[k] = -ab_x[k] * dg_1_24[k]
                      + dh_1_30[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, dg_1_25, dg_1_26, dg_1_27, \
                         dg_1_28, dg_1_29, dh_1_31, dh_1_32, dh_1_33, dh_1_34, \
                         dh_1_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * dg_1_25[k]
                      + dh_1_31[k];

            t_26[k] = -ab_x[k] * dg_1_26[k]
                      + dh_1_32[k];

            t_27[k] = -ab_x[k] * dg_1_27[k]
                      + dh_1_33[k];

            t_28[k] = -ab_x[k] * dg_1_28[k]
                      + dh_1_34[k];

            t_29[k] = -ab_x[k] * dg_1_29[k]
                      + dh_1_35[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, dg_1_30, dg_1_31, dg_1_32, \
                         dg_1_33, dg_1_34, dh_1_42, dh_1_43, dh_1_44, dh_1_45, \
                         dh_1_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * dg_1_30[k]
                      + dh_1_42[k];

            t_31[k] = -ab_x[k] * dg_1_31[k]
                      + dh_1_43[k];

            t_32[k] = -ab_x[k] * dg_1_32[k]
                      + dh_1_44[k];

            t_33[k] = -ab_x[k] * dg_1_33[k]
                      + dh_1_45[k];

            t_34[k] = -ab_x[k] * dg_1_34[k]
                      + dh_1_46[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, dg_1_35, dg_1_36, dg_1_37, \
                         dg_1_38, dg_1_39, dh_1_47, dh_1_48, dh_1_49, dh_1_50, \
                         dh_1_51 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * dg_1_35[k]
                      + dh_1_47[k];

            t_36[k] = -ab_x[k] * dg_1_36[k]
                      + dh_1_48[k];

            t_37[k] = -ab_x[k] * dg_1_37[k]
                      + dh_1_49[k];

            t_38[k] = -ab_x[k] * dg_1_38[k]
                      + dh_1_50[k];

            t_39[k] = -ab_x[k] * dg_1_39[k]
                      + dh_1_51[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, dg_1_40, dg_1_41, dg_1_42, \
                         dg_1_43, dg_1_44, dh_1_52, dh_1_53, dh_1_54, dh_1_55, \
                         dh_1_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * dg_1_40[k]
                      + dh_1_52[k];

            t_41[k] = -ab_x[k] * dg_1_41[k]
                      + dh_1_53[k];

            t_42[k] = -ab_x[k] * dg_1_42[k]
                      + dh_1_54[k];

            t_43[k] = -ab_x[k] * dg_1_43[k]
                      + dh_1_55[k];

            t_44[k] = -ab_x[k] * dg_1_44[k]
                      + dh_1_56[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, dg_1_45, dg_1_46, dg_1_47, \
                         dg_1_48, dg_1_49, dh_1_63, dh_1_64, dh_1_65, dh_1_66, \
                         dh_1_67 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * dg_1_45[k]
                      + dh_1_63[k];

            t_46[k] = -ab_x[k] * dg_1_46[k]
                      + dh_1_64[k];

            t_47[k] = -ab_x[k] * dg_1_47[k]
                      + dh_1_65[k];

            t_48[k] = -ab_x[k] * dg_1_48[k]
                      + dh_1_66[k];

            t_49[k] = -ab_x[k] * dg_1_49[k]
                      + dh_1_67[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, dg_1_50, dg_1_51, dg_1_52, \
                         dg_1_53, dg_1_54, dh_1_68, dh_1_69, dh_1_70, dh_1_71, \
                         dh_1_72 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * dg_1_50[k]
                      + dh_1_68[k];

            t_51[k] = -ab_x[k] * dg_1_51[k]
                      + dh_1_69[k];

            t_52[k] = -ab_x[k] * dg_1_52[k]
                      + dh_1_70[k];

            t_53[k] = -ab_x[k] * dg_1_53[k]
                      + dh_1_71[k];

            t_54[k] = -ab_x[k] * dg_1_54[k]
                      + dh_1_72[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, dg_1_55, dg_1_56, dg_1_57, \
                         dg_1_58, dg_1_59, dh_1_73, dh_1_74, dh_1_75, dh_1_76, \
                         dh_1_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * dg_1_55[k]
                      + dh_1_73[k];

            t_56[k] = -ab_x[k] * dg_1_56[k]
                      + dh_1_74[k];

            t_57[k] = -ab_x[k] * dg_1_57[k]
                      + dh_1_75[k];

            t_58[k] = -ab_x[k] * dg_1_58[k]
                      + dh_1_76[k];

            t_59[k] = -ab_x[k] * dg_1_59[k]
                      + dh_1_77[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, dg_1_60, dg_1_61, dg_1_62, \
                         dg_1_63, dg_1_64, dh_1_84, dh_1_85, dh_1_86, dh_1_87, \
                         dh_1_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * dg_1_60[k]
                      + dh_1_84[k];

            t_61[k] = -ab_x[k] * dg_1_61[k]
                      + dh_1_85[k];

            t_62[k] = -ab_x[k] * dg_1_62[k]
                      + dh_1_86[k];

            t_63[k] = -ab_x[k] * dg_1_63[k]
                      + dh_1_87[k];

            t_64[k] = -ab_x[k] * dg_1_64[k]
                      + dh_1_88[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, dg_1_65, dg_1_66, dg_1_67, \
                         dg_1_68, dg_1_69, dh_1_89, dh_1_90, dh_1_91, dh_1_92, \
                         dh_1_93 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * dg_1_65[k]
                      + dh_1_89[k];

            t_66[k] = -ab_x[k] * dg_1_66[k]
                      + dh_1_90[k];

            t_67[k] = -ab_x[k] * dg_1_67[k]
                      + dh_1_91[k];

            t_68[k] = -ab_x[k] * dg_1_68[k]
                      + dh_1_92[k];

            t_69[k] = -ab_x[k] * dg_1_69[k]
                      + dh_1_93[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, dg_1_70, dg_1_71, dg_1_72, \
                         dg_1_73, dg_1_74, dh_1_94, dh_1_95, dh_1_96, dh_1_97, \
                         dh_1_98 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * dg_1_70[k]
                      + dh_1_94[k];

            t_71[k] = -ab_x[k] * dg_1_71[k]
                      + dh_1_95[k];

            t_72[k] = -ab_x[k] * dg_1_72[k]
                      + dh_1_96[k];

            t_73[k] = -ab_x[k] * dg_1_73[k]
                      + dh_1_97[k];

            t_74[k] = -ab_x[k] * dg_1_74[k]
                      + dh_1_98[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, dg_1_75, dg_1_76, dg_1_77, \
                         dg_1_78, dg_1_79, dh_1_105, dh_1_106, dh_1_107, dh_1_108, \
                         dh_1_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * dg_1_75[k]
                      + dh_1_105[k];

            t_76[k] = -ab_x[k] * dg_1_76[k]
                      + dh_1_106[k];

            t_77[k] = -ab_x[k] * dg_1_77[k]
                      + dh_1_107[k];

            t_78[k] = -ab_x[k] * dg_1_78[k]
                      + dh_1_108[k];

            t_79[k] = -ab_x[k] * dg_1_79[k]
                      + dh_1_109[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, dg_1_80, dg_1_81, dg_1_82, \
                         dg_1_83, dg_1_84, dh_1_110, dh_1_111, dh_1_112, dh_1_113, \
                         dh_1_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * dg_1_80[k]
                      + dh_1_110[k];

            t_81[k] = -ab_x[k] * dg_1_81[k]
                      + dh_1_111[k];

            t_82[k] = -ab_x[k] * dg_1_82[k]
                      + dh_1_112[k];

            t_83[k] = -ab_x[k] * dg_1_83[k]
                      + dh_1_113[k];

            t_84[k] = -ab_x[k] * dg_1_84[k]
                      + dh_1_114[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, dg_1_85, dg_1_86, dg_1_87, \
                         dg_1_88, dg_1_89, dh_1_115, dh_1_116, dh_1_117, dh_1_118, \
                         dh_1_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * dg_1_85[k]
                      + dh_1_115[k];

            t_86[k] = -ab_x[k] * dg_1_86[k]
                      + dh_1_116[k];

            t_87[k] = -ab_x[k] * dg_1_87[k]
                      + dh_1_117[k];

            t_88[k] = -ab_x[k] * dg_1_88[k]
                      + dh_1_118[k];

            t_89[k] = -ab_x[k] * dg_1_89[k]
                      + dh_1_119[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_y, dg_1_45, dg_1_46, dg_1_47, \
                         dg_1_48, dg_1_49, dh_1_64, dh_1_66, dh_1_67, dh_1_69, \
                         dh_1_70 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_y[k] * dg_1_45[k]
                      + dh_1_64[k];

            t_91[k] = -ab_y[k] * dg_1_46[k]
                      + dh_1_66[k];

            t_92[k] = -ab_y[k] * dg_1_47[k]
                      + dh_1_67[k];

            t_93[k] = -ab_y[k] * dg_1_48[k]
                      + dh_1_69[k];

            t_94[k] = -ab_y[k] * dg_1_49[k]
                      + dh_1_70[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_y, dg_1_50, dg_1_51, dg_1_52, \
                         dg_1_53, dg_1_54, dh_1_71, dh_1_73, dh_1_74, dh_1_75, \
                         dh_1_76 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_y[k] * dg_1_50[k]
                      + dh_1_71[k];

            t_96[k] = -ab_y[k] * dg_1_51[k]
                      + dh_1_73[k];

            t_97[k] = -ab_y[k] * dg_1_52[k]
                      + dh_1_74[k];

            t_98[k] = -ab_y[k] * dg_1_53[k]
                      + dh_1_75[k];

            t_99[k] = -ab_y[k] * dg_1_54[k]
                      + dh_1_76[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_y, dg_1_55, dg_1_56, dg_1_57, \
                         dg_1_58, dg_1_59, dh_1_78, dh_1_79, dh_1_80, dh_1_81, \
                         dh_1_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_y[k] * dg_1_55[k]
                       + dh_1_78[k];

            t_101[k] = -ab_y[k] * dg_1_56[k]
                       + dh_1_79[k];

            t_102[k] = -ab_y[k] * dg_1_57[k]
                       + dh_1_80[k];

            t_103[k] = -ab_y[k] * dg_1_58[k]
                       + dh_1_81[k];

            t_104[k] = -ab_y[k] * dg_1_59[k]
                       + dh_1_82[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_y, dg_1_60, dg_1_61, dg_1_62, \
                         dg_1_63, dg_1_64, dh_1_85, dh_1_87, dh_1_88, dh_1_90, \
                         dh_1_91 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_y[k] * dg_1_60[k]
                       + dh_1_85[k];

            t_106[k] = -ab_y[k] * dg_1_61[k]
                       + dh_1_87[k];

            t_107[k] = -ab_y[k] * dg_1_62[k]
                       + dh_1_88[k];

            t_108[k] = -ab_y[k] * dg_1_63[k]
                       + dh_1_90[k];

            t_109[k] = -ab_y[k] * dg_1_64[k]
                       + dh_1_91[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_y, dg_1_65, dg_1_66, dg_1_67, \
                         dg_1_68, dg_1_69, dh_1_92, dh_1_94, dh_1_95, dh_1_96, \
                         dh_1_97 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = -ab_y[k] * dg_1_65[k]
                       + dh_1_92[k];

            t_111[k] = -ab_y[k] * dg_1_66[k]
                       + dh_1_94[k];

            t_112[k] = -ab_y[k] * dg_1_67[k]
                       + dh_1_95[k];

            t_113[k] = -ab_y[k] * dg_1_68[k]
                       + dh_1_96[k];

            t_114[k] = -ab_y[k] * dg_1_69[k]
                       + dh_1_97[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_y, dg_1_70, dg_1_71, dg_1_72, \
                         dg_1_73, dg_1_74, dh_1_99, dh_1_100, dh_1_101, dh_1_102, \
                         dh_1_103 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = -ab_y[k] * dg_1_70[k]
                       + dh_1_99[k];

            t_116[k] = -ab_y[k] * dg_1_71[k]
                       + dh_1_100[k];

            t_117[k] = -ab_y[k] * dg_1_72[k]
                       + dh_1_101[k];

            t_118[k] = -ab_y[k] * dg_1_73[k]
                       + dh_1_102[k];

            t_119[k] = -ab_y[k] * dg_1_74[k]
                       + dh_1_103[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_y, dg_1_75, dg_1_76, dg_1_77, \
                         dg_1_78, dg_1_79, dh_1_106, dh_1_108, dh_1_109, dh_1_111, \
                         dh_1_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_y[k] * dg_1_75[k]
                       + dh_1_106[k];

            t_121[k] = -ab_y[k] * dg_1_76[k]
                       + dh_1_108[k];

            t_122[k] = -ab_y[k] * dg_1_77[k]
                       + dh_1_109[k];

            t_123[k] = -ab_y[k] * dg_1_78[k]
                       + dh_1_111[k];

            t_124[k] = -ab_y[k] * dg_1_79[k]
                       + dh_1_112[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_y, dg_1_80, dg_1_81, dg_1_82, \
                         dg_1_83, dg_1_84, dh_1_113, dh_1_115, dh_1_116, dh_1_117, \
                         dh_1_118 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = -ab_y[k] * dg_1_80[k]
                       + dh_1_113[k];

            t_126[k] = -ab_y[k] * dg_1_81[k]
                       + dh_1_115[k];

            t_127[k] = -ab_y[k] * dg_1_82[k]
                       + dh_1_116[k];

            t_128[k] = -ab_y[k] * dg_1_83[k]
                       + dh_1_117[k];

            t_129[k] = -ab_y[k] * dg_1_84[k]
                       + dh_1_118[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_y, dg_1_85, dg_1_86, dg_1_87, \
                         dg_1_88, dg_1_89, dh_1_120, dh_1_121, dh_1_122, dh_1_123, \
                         dh_1_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = -ab_y[k] * dg_1_85[k]
                       + dh_1_120[k];

            t_131[k] = -ab_y[k] * dg_1_86[k]
                       + dh_1_121[k];

            t_132[k] = -ab_y[k] * dg_1_87[k]
                       + dh_1_122[k];

            t_133[k] = -ab_y[k] * dg_1_88[k]
                       + dh_1_123[k];

            t_134[k] = -ab_y[k] * dg_1_89[k]
                       + dh_1_124[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, ab_z, dg_1_75, dg_1_76, dg_1_77, dg_0_75, \
                         dg_0_76, dg_0_77, dh_1_107, dh_1_109, \
                         dh_1_110 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = -ab_z[k] * dg_1_75[k]
                       + dg_0_75[k]
                       + dh_1_107[k];

            t_136[k] = -ab_z[k] * dg_1_76[k]
                       + dg_0_76[k]
                       + dh_1_109[k];

            t_137[k] = -ab_z[k] * dg_1_77[k]
                       + dg_0_77[k]
                       + dh_1_110[k];
        }

#pragma omp simd aligned(t_138, t_139, t_140, ab_z, dg_1_78, dg_1_79, dg_1_80, dg_0_78, \
                         dg_0_79, dg_0_80, dh_1_112, dh_1_113, \
                         dh_1_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_138[k] = -ab_z[k] * dg_1_78[k]
                       + dg_0_78[k]
                       + dh_1_112[k];

            t_139[k] = -ab_z[k] * dg_1_79[k]
                       + dg_0_79[k]
                       + dh_1_113[k];

            t_140[k] = -ab_z[k] * dg_1_80[k]
                       + dg_0_80[k]
                       + dh_1_114[k];
        }
    }
}

static auto
compute_hrr_geom_010z_fg_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                const size_t target, const size_t dg_1, const size_t dg_0,
                                const size_t dh_1, const size_t ncomps,
                                const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_141 = buffer.data(target + 141 * ncomps + c);
        auto *t_142 = buffer.data(target + 142 * ncomps + c);
        auto *t_143 = buffer.data(target + 143 * ncomps + c);
        auto *t_144 = buffer.data(target + 144 * ncomps + c);
        auto *t_145 = buffer.data(target + 145 * ncomps + c);
        auto *t_146 = buffer.data(target + 146 * ncomps + c);
        auto *t_147 = buffer.data(target + 147 * ncomps + c);
        auto *t_148 = buffer.data(target + 148 * ncomps + c);
        auto *t_149 = buffer.data(target + 149 * ncomps + c);

        const auto *ab_z = coordinates.data(8);

        const auto *dg_1_81 = buffer.data(dg_1 + 81 * ncomps + c);
        const auto *dg_1_82 = buffer.data(dg_1 + 82 * ncomps + c);
        const auto *dg_1_83 = buffer.data(dg_1 + 83 * ncomps + c);
        const auto *dg_1_84 = buffer.data(dg_1 + 84 * ncomps + c);
        const auto *dg_1_85 = buffer.data(dg_1 + 85 * ncomps + c);
        const auto *dg_1_86 = buffer.data(dg_1 + 86 * ncomps + c);
        const auto *dg_1_87 = buffer.data(dg_1 + 87 * ncomps + c);
        const auto *dg_1_88 = buffer.data(dg_1 + 88 * ncomps + c);
        const auto *dg_1_89 = buffer.data(dg_1 + 89 * ncomps + c);

        const auto *dg_0_81 = buffer.data(dg_0 + 81 * ncomps + c);
        const auto *dg_0_82 = buffer.data(dg_0 + 82 * ncomps + c);
        const auto *dg_0_83 = buffer.data(dg_0 + 83 * ncomps + c);
        const auto *dg_0_84 = buffer.data(dg_0 + 84 * ncomps + c);
        const auto *dg_0_85 = buffer.data(dg_0 + 85 * ncomps + c);
        const auto *dg_0_86 = buffer.data(dg_0 + 86 * ncomps + c);
        const auto *dg_0_87 = buffer.data(dg_0 + 87 * ncomps + c);
        const auto *dg_0_88 = buffer.data(dg_0 + 88 * ncomps + c);
        const auto *dg_0_89 = buffer.data(dg_0 + 89 * ncomps + c);

        const auto *dh_1_116 = buffer.data(dh_1 + 116 * ncomps + c);
        const auto *dh_1_117 = buffer.data(dh_1 + 117 * ncomps + c);
        const auto *dh_1_118 = buffer.data(dh_1 + 118 * ncomps + c);
        const auto *dh_1_119 = buffer.data(dh_1 + 119 * ncomps + c);
        const auto *dh_1_121 = buffer.data(dh_1 + 121 * ncomps + c);
        const auto *dh_1_122 = buffer.data(dh_1 + 122 * ncomps + c);
        const auto *dh_1_123 = buffer.data(dh_1 + 123 * ncomps + c);
        const auto *dh_1_124 = buffer.data(dh_1 + 124 * ncomps + c);
        const auto *dh_1_125 = buffer.data(dh_1 + 125 * ncomps + c);

#pragma omp simd aligned(t_141, t_142, t_143, ab_z, dg_1_81, dg_1_82, dg_1_83, dg_0_81, \
                         dg_0_82, dg_0_83, dh_1_116, dh_1_117, \
                         dh_1_118 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_141[k] = -ab_z[k] * dg_1_81[k]
                       + dg_0_81[k]
                       + dh_1_116[k];

            t_142[k] = -ab_z[k] * dg_1_82[k]
                       + dg_0_82[k]
                       + dh_1_117[k];

            t_143[k] = -ab_z[k] * dg_1_83[k]
                       + dg_0_83[k]
                       + dh_1_118[k];
        }

#pragma omp simd aligned(t_144, t_145, t_146, ab_z, dg_1_84, dg_1_85, dg_1_86, dg_0_84, \
                         dg_0_85, dg_0_86, dh_1_119, dh_1_121, \
                         dh_1_122 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_144[k] = -ab_z[k] * dg_1_84[k]
                       + dg_0_84[k]
                       + dh_1_119[k];

            t_145[k] = -ab_z[k] * dg_1_85[k]
                       + dg_0_85[k]
                       + dh_1_121[k];

            t_146[k] = -ab_z[k] * dg_1_86[k]
                       + dg_0_86[k]
                       + dh_1_122[k];
        }

#pragma omp simd aligned(t_147, t_148, t_149, ab_z, dg_1_87, dg_1_88, dg_1_89, dg_0_87, \
                         dg_0_88, dg_0_89, dh_1_123, dh_1_124, \
                         dh_1_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_147[k] = -ab_z[k] * dg_1_87[k]
                       + dg_0_87[k]
                       + dh_1_123[k];

            t_148[k] = -ab_z[k] * dg_1_88[k]
                       + dg_0_88[k]
                       + dh_1_124[k];

            t_149[k] = -ab_z[k] * dg_1_89[k]
                       + dg_0_89[k]
                       + dh_1_125[k];
        }
    }
}

auto
compute_hrr_geom_010z_fg(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                         const size_t target, const size_t dg_1, const size_t dg_0,
                         const size_t dh_1, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_geom_010z_fg_piece0(buffer, coordinates, target, dg_1, dg_0, dh_1, ncomps,
                                    nmax);

    compute_hrr_geom_010z_fg_piece1(buffer, coordinates, target, dg_1, dg_0, dh_1, ncomps,
                                    nmax);
}

}  // namespace simdtrf
