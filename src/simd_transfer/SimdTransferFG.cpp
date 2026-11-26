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


#include "SimdTransferFG.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_fg_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t dg, const size_t dh, const size_t ncomps,
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
        auto *t_141 = buffer.data(target + 141 * ncomps + c);
        auto *t_142 = buffer.data(target + 142 * ncomps + c);
        auto *t_143 = buffer.data(target + 143 * ncomps + c);
        auto *t_144 = buffer.data(target + 144 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *dg_0 = buffer.data(dg + 0 * ncomps + c);
        const auto *dg_1 = buffer.data(dg + 1 * ncomps + c);
        const auto *dg_2 = buffer.data(dg + 2 * ncomps + c);
        const auto *dg_3 = buffer.data(dg + 3 * ncomps + c);
        const auto *dg_4 = buffer.data(dg + 4 * ncomps + c);
        const auto *dg_5 = buffer.data(dg + 5 * ncomps + c);
        const auto *dg_6 = buffer.data(dg + 6 * ncomps + c);
        const auto *dg_7 = buffer.data(dg + 7 * ncomps + c);
        const auto *dg_8 = buffer.data(dg + 8 * ncomps + c);
        const auto *dg_9 = buffer.data(dg + 9 * ncomps + c);
        const auto *dg_10 = buffer.data(dg + 10 * ncomps + c);
        const auto *dg_11 = buffer.data(dg + 11 * ncomps + c);
        const auto *dg_12 = buffer.data(dg + 12 * ncomps + c);
        const auto *dg_13 = buffer.data(dg + 13 * ncomps + c);
        const auto *dg_14 = buffer.data(dg + 14 * ncomps + c);
        const auto *dg_15 = buffer.data(dg + 15 * ncomps + c);
        const auto *dg_16 = buffer.data(dg + 16 * ncomps + c);
        const auto *dg_17 = buffer.data(dg + 17 * ncomps + c);
        const auto *dg_18 = buffer.data(dg + 18 * ncomps + c);
        const auto *dg_19 = buffer.data(dg + 19 * ncomps + c);
        const auto *dg_20 = buffer.data(dg + 20 * ncomps + c);
        const auto *dg_21 = buffer.data(dg + 21 * ncomps + c);
        const auto *dg_22 = buffer.data(dg + 22 * ncomps + c);
        const auto *dg_23 = buffer.data(dg + 23 * ncomps + c);
        const auto *dg_24 = buffer.data(dg + 24 * ncomps + c);
        const auto *dg_25 = buffer.data(dg + 25 * ncomps + c);
        const auto *dg_26 = buffer.data(dg + 26 * ncomps + c);
        const auto *dg_27 = buffer.data(dg + 27 * ncomps + c);
        const auto *dg_28 = buffer.data(dg + 28 * ncomps + c);
        const auto *dg_29 = buffer.data(dg + 29 * ncomps + c);
        const auto *dg_30 = buffer.data(dg + 30 * ncomps + c);
        const auto *dg_31 = buffer.data(dg + 31 * ncomps + c);
        const auto *dg_32 = buffer.data(dg + 32 * ncomps + c);
        const auto *dg_33 = buffer.data(dg + 33 * ncomps + c);
        const auto *dg_34 = buffer.data(dg + 34 * ncomps + c);
        const auto *dg_35 = buffer.data(dg + 35 * ncomps + c);
        const auto *dg_36 = buffer.data(dg + 36 * ncomps + c);
        const auto *dg_37 = buffer.data(dg + 37 * ncomps + c);
        const auto *dg_38 = buffer.data(dg + 38 * ncomps + c);
        const auto *dg_39 = buffer.data(dg + 39 * ncomps + c);
        const auto *dg_40 = buffer.data(dg + 40 * ncomps + c);
        const auto *dg_41 = buffer.data(dg + 41 * ncomps + c);
        const auto *dg_42 = buffer.data(dg + 42 * ncomps + c);
        const auto *dg_43 = buffer.data(dg + 43 * ncomps + c);
        const auto *dg_44 = buffer.data(dg + 44 * ncomps + c);
        const auto *dg_45 = buffer.data(dg + 45 * ncomps + c);
        const auto *dg_46 = buffer.data(dg + 46 * ncomps + c);
        const auto *dg_47 = buffer.data(dg + 47 * ncomps + c);
        const auto *dg_48 = buffer.data(dg + 48 * ncomps + c);
        const auto *dg_49 = buffer.data(dg + 49 * ncomps + c);
        const auto *dg_50 = buffer.data(dg + 50 * ncomps + c);
        const auto *dg_51 = buffer.data(dg + 51 * ncomps + c);
        const auto *dg_52 = buffer.data(dg + 52 * ncomps + c);
        const auto *dg_53 = buffer.data(dg + 53 * ncomps + c);
        const auto *dg_54 = buffer.data(dg + 54 * ncomps + c);
        const auto *dg_55 = buffer.data(dg + 55 * ncomps + c);
        const auto *dg_56 = buffer.data(dg + 56 * ncomps + c);
        const auto *dg_57 = buffer.data(dg + 57 * ncomps + c);
        const auto *dg_58 = buffer.data(dg + 58 * ncomps + c);
        const auto *dg_59 = buffer.data(dg + 59 * ncomps + c);
        const auto *dg_60 = buffer.data(dg + 60 * ncomps + c);
        const auto *dg_61 = buffer.data(dg + 61 * ncomps + c);
        const auto *dg_62 = buffer.data(dg + 62 * ncomps + c);
        const auto *dg_63 = buffer.data(dg + 63 * ncomps + c);
        const auto *dg_64 = buffer.data(dg + 64 * ncomps + c);
        const auto *dg_65 = buffer.data(dg + 65 * ncomps + c);
        const auto *dg_66 = buffer.data(dg + 66 * ncomps + c);
        const auto *dg_67 = buffer.data(dg + 67 * ncomps + c);
        const auto *dg_68 = buffer.data(dg + 68 * ncomps + c);
        const auto *dg_69 = buffer.data(dg + 69 * ncomps + c);
        const auto *dg_70 = buffer.data(dg + 70 * ncomps + c);
        const auto *dg_71 = buffer.data(dg + 71 * ncomps + c);
        const auto *dg_72 = buffer.data(dg + 72 * ncomps + c);
        const auto *dg_73 = buffer.data(dg + 73 * ncomps + c);
        const auto *dg_74 = buffer.data(dg + 74 * ncomps + c);
        const auto *dg_75 = buffer.data(dg + 75 * ncomps + c);
        const auto *dg_76 = buffer.data(dg + 76 * ncomps + c);
        const auto *dg_77 = buffer.data(dg + 77 * ncomps + c);
        const auto *dg_78 = buffer.data(dg + 78 * ncomps + c);
        const auto *dg_79 = buffer.data(dg + 79 * ncomps + c);
        const auto *dg_80 = buffer.data(dg + 80 * ncomps + c);
        const auto *dg_81 = buffer.data(dg + 81 * ncomps + c);
        const auto *dg_82 = buffer.data(dg + 82 * ncomps + c);
        const auto *dg_83 = buffer.data(dg + 83 * ncomps + c);
        const auto *dg_84 = buffer.data(dg + 84 * ncomps + c);
        const auto *dg_85 = buffer.data(dg + 85 * ncomps + c);
        const auto *dg_86 = buffer.data(dg + 86 * ncomps + c);
        const auto *dg_87 = buffer.data(dg + 87 * ncomps + c);
        const auto *dg_88 = buffer.data(dg + 88 * ncomps + c);
        const auto *dg_89 = buffer.data(dg + 89 * ncomps + c);

        const auto *dh_0 = buffer.data(dh + 0 * ncomps + c);
        const auto *dh_1 = buffer.data(dh + 1 * ncomps + c);
        const auto *dh_2 = buffer.data(dh + 2 * ncomps + c);
        const auto *dh_3 = buffer.data(dh + 3 * ncomps + c);
        const auto *dh_4 = buffer.data(dh + 4 * ncomps + c);
        const auto *dh_5 = buffer.data(dh + 5 * ncomps + c);
        const auto *dh_6 = buffer.data(dh + 6 * ncomps + c);
        const auto *dh_7 = buffer.data(dh + 7 * ncomps + c);
        const auto *dh_8 = buffer.data(dh + 8 * ncomps + c);
        const auto *dh_9 = buffer.data(dh + 9 * ncomps + c);
        const auto *dh_10 = buffer.data(dh + 10 * ncomps + c);
        const auto *dh_11 = buffer.data(dh + 11 * ncomps + c);
        const auto *dh_12 = buffer.data(dh + 12 * ncomps + c);
        const auto *dh_13 = buffer.data(dh + 13 * ncomps + c);
        const auto *dh_14 = buffer.data(dh + 14 * ncomps + c);
        const auto *dh_21 = buffer.data(dh + 21 * ncomps + c);
        const auto *dh_22 = buffer.data(dh + 22 * ncomps + c);
        const auto *dh_23 = buffer.data(dh + 23 * ncomps + c);
        const auto *dh_24 = buffer.data(dh + 24 * ncomps + c);
        const auto *dh_25 = buffer.data(dh + 25 * ncomps + c);
        const auto *dh_26 = buffer.data(dh + 26 * ncomps + c);
        const auto *dh_27 = buffer.data(dh + 27 * ncomps + c);
        const auto *dh_28 = buffer.data(dh + 28 * ncomps + c);
        const auto *dh_29 = buffer.data(dh + 29 * ncomps + c);
        const auto *dh_30 = buffer.data(dh + 30 * ncomps + c);
        const auto *dh_31 = buffer.data(dh + 31 * ncomps + c);
        const auto *dh_32 = buffer.data(dh + 32 * ncomps + c);
        const auto *dh_33 = buffer.data(dh + 33 * ncomps + c);
        const auto *dh_34 = buffer.data(dh + 34 * ncomps + c);
        const auto *dh_35 = buffer.data(dh + 35 * ncomps + c);
        const auto *dh_42 = buffer.data(dh + 42 * ncomps + c);
        const auto *dh_43 = buffer.data(dh + 43 * ncomps + c);
        const auto *dh_44 = buffer.data(dh + 44 * ncomps + c);
        const auto *dh_45 = buffer.data(dh + 45 * ncomps + c);
        const auto *dh_46 = buffer.data(dh + 46 * ncomps + c);
        const auto *dh_47 = buffer.data(dh + 47 * ncomps + c);
        const auto *dh_48 = buffer.data(dh + 48 * ncomps + c);
        const auto *dh_49 = buffer.data(dh + 49 * ncomps + c);
        const auto *dh_50 = buffer.data(dh + 50 * ncomps + c);
        const auto *dh_51 = buffer.data(dh + 51 * ncomps + c);
        const auto *dh_52 = buffer.data(dh + 52 * ncomps + c);
        const auto *dh_53 = buffer.data(dh + 53 * ncomps + c);
        const auto *dh_54 = buffer.data(dh + 54 * ncomps + c);
        const auto *dh_55 = buffer.data(dh + 55 * ncomps + c);
        const auto *dh_56 = buffer.data(dh + 56 * ncomps + c);
        const auto *dh_63 = buffer.data(dh + 63 * ncomps + c);
        const auto *dh_64 = buffer.data(dh + 64 * ncomps + c);
        const auto *dh_65 = buffer.data(dh + 65 * ncomps + c);
        const auto *dh_66 = buffer.data(dh + 66 * ncomps + c);
        const auto *dh_67 = buffer.data(dh + 67 * ncomps + c);
        const auto *dh_68 = buffer.data(dh + 68 * ncomps + c);
        const auto *dh_69 = buffer.data(dh + 69 * ncomps + c);
        const auto *dh_70 = buffer.data(dh + 70 * ncomps + c);
        const auto *dh_71 = buffer.data(dh + 71 * ncomps + c);
        const auto *dh_72 = buffer.data(dh + 72 * ncomps + c);
        const auto *dh_73 = buffer.data(dh + 73 * ncomps + c);
        const auto *dh_74 = buffer.data(dh + 74 * ncomps + c);
        const auto *dh_75 = buffer.data(dh + 75 * ncomps + c);
        const auto *dh_76 = buffer.data(dh + 76 * ncomps + c);
        const auto *dh_77 = buffer.data(dh + 77 * ncomps + c);
        const auto *dh_78 = buffer.data(dh + 78 * ncomps + c);
        const auto *dh_79 = buffer.data(dh + 79 * ncomps + c);
        const auto *dh_80 = buffer.data(dh + 80 * ncomps + c);
        const auto *dh_81 = buffer.data(dh + 81 * ncomps + c);
        const auto *dh_82 = buffer.data(dh + 82 * ncomps + c);
        const auto *dh_84 = buffer.data(dh + 84 * ncomps + c);
        const auto *dh_85 = buffer.data(dh + 85 * ncomps + c);
        const auto *dh_86 = buffer.data(dh + 86 * ncomps + c);
        const auto *dh_87 = buffer.data(dh + 87 * ncomps + c);
        const auto *dh_88 = buffer.data(dh + 88 * ncomps + c);
        const auto *dh_89 = buffer.data(dh + 89 * ncomps + c);
        const auto *dh_90 = buffer.data(dh + 90 * ncomps + c);
        const auto *dh_91 = buffer.data(dh + 91 * ncomps + c);
        const auto *dh_92 = buffer.data(dh + 92 * ncomps + c);
        const auto *dh_93 = buffer.data(dh + 93 * ncomps + c);
        const auto *dh_94 = buffer.data(dh + 94 * ncomps + c);
        const auto *dh_95 = buffer.data(dh + 95 * ncomps + c);
        const auto *dh_96 = buffer.data(dh + 96 * ncomps + c);
        const auto *dh_97 = buffer.data(dh + 97 * ncomps + c);
        const auto *dh_98 = buffer.data(dh + 98 * ncomps + c);
        const auto *dh_99 = buffer.data(dh + 99 * ncomps + c);
        const auto *dh_100 = buffer.data(dh + 100 * ncomps + c);
        const auto *dh_101 = buffer.data(dh + 101 * ncomps + c);
        const auto *dh_102 = buffer.data(dh + 102 * ncomps + c);
        const auto *dh_103 = buffer.data(dh + 103 * ncomps + c);
        const auto *dh_105 = buffer.data(dh + 105 * ncomps + c);
        const auto *dh_106 = buffer.data(dh + 106 * ncomps + c);
        const auto *dh_107 = buffer.data(dh + 107 * ncomps + c);
        const auto *dh_108 = buffer.data(dh + 108 * ncomps + c);
        const auto *dh_109 = buffer.data(dh + 109 * ncomps + c);
        const auto *dh_110 = buffer.data(dh + 110 * ncomps + c);
        const auto *dh_111 = buffer.data(dh + 111 * ncomps + c);
        const auto *dh_112 = buffer.data(dh + 112 * ncomps + c);
        const auto *dh_113 = buffer.data(dh + 113 * ncomps + c);
        const auto *dh_114 = buffer.data(dh + 114 * ncomps + c);
        const auto *dh_115 = buffer.data(dh + 115 * ncomps + c);
        const auto *dh_116 = buffer.data(dh + 116 * ncomps + c);
        const auto *dh_117 = buffer.data(dh + 117 * ncomps + c);
        const auto *dh_118 = buffer.data(dh + 118 * ncomps + c);
        const auto *dh_119 = buffer.data(dh + 119 * ncomps + c);
        const auto *dh_120 = buffer.data(dh + 120 * ncomps + c);
        const auto *dh_121 = buffer.data(dh + 121 * ncomps + c);
        const auto *dh_122 = buffer.data(dh + 122 * ncomps + c);
        const auto *dh_123 = buffer.data(dh + 123 * ncomps + c);
        const auto *dh_124 = buffer.data(dh + 124 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, dg_0, dg_1, dg_2, dg_3, dg_4, dh_0, \
                         dh_1, dh_2, dh_3, dh_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * dg_0[k]
                     + dh_0[k];

            t_1[k] = -ab_x[k] * dg_1[k]
                     + dh_1[k];

            t_2[k] = -ab_x[k] * dg_2[k]
                     + dh_2[k];

            t_3[k] = -ab_x[k] * dg_3[k]
                     + dh_3[k];

            t_4[k] = -ab_x[k] * dg_4[k]
                     + dh_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, dg_5, dg_6, dg_7, dg_8, dg_9, dh_5, \
                         dh_6, dh_7, dh_8, dh_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * dg_5[k]
                     + dh_5[k];

            t_6[k] = -ab_x[k] * dg_6[k]
                     + dh_6[k];

            t_7[k] = -ab_x[k] * dg_7[k]
                     + dh_7[k];

            t_8[k] = -ab_x[k] * dg_8[k]
                     + dh_8[k];

            t_9[k] = -ab_x[k] * dg_9[k]
                     + dh_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, dg_10, dg_11, dg_12, dg_13, \
                         dg_14, dh_10, dh_11, dh_12, dh_13, dh_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * dg_10[k]
                      + dh_10[k];

            t_11[k] = -ab_x[k] * dg_11[k]
                      + dh_11[k];

            t_12[k] = -ab_x[k] * dg_12[k]
                      + dh_12[k];

            t_13[k] = -ab_x[k] * dg_13[k]
                      + dh_13[k];

            t_14[k] = -ab_x[k] * dg_14[k]
                      + dh_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, dg_15, dg_16, dg_17, dg_18, \
                         dg_19, dh_21, dh_22, dh_23, dh_24, dh_25 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * dg_15[k]
                      + dh_21[k];

            t_16[k] = -ab_x[k] * dg_16[k]
                      + dh_22[k];

            t_17[k] = -ab_x[k] * dg_17[k]
                      + dh_23[k];

            t_18[k] = -ab_x[k] * dg_18[k]
                      + dh_24[k];

            t_19[k] = -ab_x[k] * dg_19[k]
                      + dh_25[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, dg_20, dg_21, dg_22, dg_23, \
                         dg_24, dh_26, dh_27, dh_28, dh_29, dh_30 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * dg_20[k]
                      + dh_26[k];

            t_21[k] = -ab_x[k] * dg_21[k]
                      + dh_27[k];

            t_22[k] = -ab_x[k] * dg_22[k]
                      + dh_28[k];

            t_23[k] = -ab_x[k] * dg_23[k]
                      + dh_29[k];

            t_24[k] = -ab_x[k] * dg_24[k]
                      + dh_30[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, dg_25, dg_26, dg_27, dg_28, \
                         dg_29, dh_31, dh_32, dh_33, dh_34, dh_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * dg_25[k]
                      + dh_31[k];

            t_26[k] = -ab_x[k] * dg_26[k]
                      + dh_32[k];

            t_27[k] = -ab_x[k] * dg_27[k]
                      + dh_33[k];

            t_28[k] = -ab_x[k] * dg_28[k]
                      + dh_34[k];

            t_29[k] = -ab_x[k] * dg_29[k]
                      + dh_35[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, dg_30, dg_31, dg_32, dg_33, \
                         dg_34, dh_42, dh_43, dh_44, dh_45, dh_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * dg_30[k]
                      + dh_42[k];

            t_31[k] = -ab_x[k] * dg_31[k]
                      + dh_43[k];

            t_32[k] = -ab_x[k] * dg_32[k]
                      + dh_44[k];

            t_33[k] = -ab_x[k] * dg_33[k]
                      + dh_45[k];

            t_34[k] = -ab_x[k] * dg_34[k]
                      + dh_46[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, dg_35, dg_36, dg_37, dg_38, \
                         dg_39, dh_47, dh_48, dh_49, dh_50, dh_51 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * dg_35[k]
                      + dh_47[k];

            t_36[k] = -ab_x[k] * dg_36[k]
                      + dh_48[k];

            t_37[k] = -ab_x[k] * dg_37[k]
                      + dh_49[k];

            t_38[k] = -ab_x[k] * dg_38[k]
                      + dh_50[k];

            t_39[k] = -ab_x[k] * dg_39[k]
                      + dh_51[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, dg_40, dg_41, dg_42, dg_43, \
                         dg_44, dh_52, dh_53, dh_54, dh_55, dh_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * dg_40[k]
                      + dh_52[k];

            t_41[k] = -ab_x[k] * dg_41[k]
                      + dh_53[k];

            t_42[k] = -ab_x[k] * dg_42[k]
                      + dh_54[k];

            t_43[k] = -ab_x[k] * dg_43[k]
                      + dh_55[k];

            t_44[k] = -ab_x[k] * dg_44[k]
                      + dh_56[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, dg_45, dg_46, dg_47, dg_48, \
                         dg_49, dh_63, dh_64, dh_65, dh_66, dh_67 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * dg_45[k]
                      + dh_63[k];

            t_46[k] = -ab_x[k] * dg_46[k]
                      + dh_64[k];

            t_47[k] = -ab_x[k] * dg_47[k]
                      + dh_65[k];

            t_48[k] = -ab_x[k] * dg_48[k]
                      + dh_66[k];

            t_49[k] = -ab_x[k] * dg_49[k]
                      + dh_67[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, dg_50, dg_51, dg_52, dg_53, \
                         dg_54, dh_68, dh_69, dh_70, dh_71, dh_72 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * dg_50[k]
                      + dh_68[k];

            t_51[k] = -ab_x[k] * dg_51[k]
                      + dh_69[k];

            t_52[k] = -ab_x[k] * dg_52[k]
                      + dh_70[k];

            t_53[k] = -ab_x[k] * dg_53[k]
                      + dh_71[k];

            t_54[k] = -ab_x[k] * dg_54[k]
                      + dh_72[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, dg_55, dg_56, dg_57, dg_58, \
                         dg_59, dh_73, dh_74, dh_75, dh_76, dh_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * dg_55[k]
                      + dh_73[k];

            t_56[k] = -ab_x[k] * dg_56[k]
                      + dh_74[k];

            t_57[k] = -ab_x[k] * dg_57[k]
                      + dh_75[k];

            t_58[k] = -ab_x[k] * dg_58[k]
                      + dh_76[k];

            t_59[k] = -ab_x[k] * dg_59[k]
                      + dh_77[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, dg_60, dg_61, dg_62, dg_63, \
                         dg_64, dh_84, dh_85, dh_86, dh_87, dh_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * dg_60[k]
                      + dh_84[k];

            t_61[k] = -ab_x[k] * dg_61[k]
                      + dh_85[k];

            t_62[k] = -ab_x[k] * dg_62[k]
                      + dh_86[k];

            t_63[k] = -ab_x[k] * dg_63[k]
                      + dh_87[k];

            t_64[k] = -ab_x[k] * dg_64[k]
                      + dh_88[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, dg_65, dg_66, dg_67, dg_68, \
                         dg_69, dh_89, dh_90, dh_91, dh_92, dh_93 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * dg_65[k]
                      + dh_89[k];

            t_66[k] = -ab_x[k] * dg_66[k]
                      + dh_90[k];

            t_67[k] = -ab_x[k] * dg_67[k]
                      + dh_91[k];

            t_68[k] = -ab_x[k] * dg_68[k]
                      + dh_92[k];

            t_69[k] = -ab_x[k] * dg_69[k]
                      + dh_93[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, dg_70, dg_71, dg_72, dg_73, \
                         dg_74, dh_94, dh_95, dh_96, dh_97, dh_98 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * dg_70[k]
                      + dh_94[k];

            t_71[k] = -ab_x[k] * dg_71[k]
                      + dh_95[k];

            t_72[k] = -ab_x[k] * dg_72[k]
                      + dh_96[k];

            t_73[k] = -ab_x[k] * dg_73[k]
                      + dh_97[k];

            t_74[k] = -ab_x[k] * dg_74[k]
                      + dh_98[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, dg_75, dg_76, dg_77, dg_78, \
                         dg_79, dh_105, dh_106, dh_107, dh_108, \
                         dh_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * dg_75[k]
                      + dh_105[k];

            t_76[k] = -ab_x[k] * dg_76[k]
                      + dh_106[k];

            t_77[k] = -ab_x[k] * dg_77[k]
                      + dh_107[k];

            t_78[k] = -ab_x[k] * dg_78[k]
                      + dh_108[k];

            t_79[k] = -ab_x[k] * dg_79[k]
                      + dh_109[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, dg_80, dg_81, dg_82, dg_83, \
                         dg_84, dh_110, dh_111, dh_112, dh_113, \
                         dh_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * dg_80[k]
                      + dh_110[k];

            t_81[k] = -ab_x[k] * dg_81[k]
                      + dh_111[k];

            t_82[k] = -ab_x[k] * dg_82[k]
                      + dh_112[k];

            t_83[k] = -ab_x[k] * dg_83[k]
                      + dh_113[k];

            t_84[k] = -ab_x[k] * dg_84[k]
                      + dh_114[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, dg_85, dg_86, dg_87, dg_88, \
                         dg_89, dh_115, dh_116, dh_117, dh_118, \
                         dh_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * dg_85[k]
                      + dh_115[k];

            t_86[k] = -ab_x[k] * dg_86[k]
                      + dh_116[k];

            t_87[k] = -ab_x[k] * dg_87[k]
                      + dh_117[k];

            t_88[k] = -ab_x[k] * dg_88[k]
                      + dh_118[k];

            t_89[k] = -ab_x[k] * dg_89[k]
                      + dh_119[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_y, dg_45, dg_46, dg_47, dg_48, \
                         dg_49, dh_64, dh_66, dh_67, dh_69, dh_70 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_y[k] * dg_45[k]
                      + dh_64[k];

            t_91[k] = -ab_y[k] * dg_46[k]
                      + dh_66[k];

            t_92[k] = -ab_y[k] * dg_47[k]
                      + dh_67[k];

            t_93[k] = -ab_y[k] * dg_48[k]
                      + dh_69[k];

            t_94[k] = -ab_y[k] * dg_49[k]
                      + dh_70[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_y, dg_50, dg_51, dg_52, dg_53, \
                         dg_54, dh_71, dh_73, dh_74, dh_75, dh_76 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_y[k] * dg_50[k]
                      + dh_71[k];

            t_96[k] = -ab_y[k] * dg_51[k]
                      + dh_73[k];

            t_97[k] = -ab_y[k] * dg_52[k]
                      + dh_74[k];

            t_98[k] = -ab_y[k] * dg_53[k]
                      + dh_75[k];

            t_99[k] = -ab_y[k] * dg_54[k]
                      + dh_76[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_y, dg_55, dg_56, dg_57, dg_58, \
                         dg_59, dh_78, dh_79, dh_80, dh_81, dh_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_y[k] * dg_55[k]
                       + dh_78[k];

            t_101[k] = -ab_y[k] * dg_56[k]
                       + dh_79[k];

            t_102[k] = -ab_y[k] * dg_57[k]
                       + dh_80[k];

            t_103[k] = -ab_y[k] * dg_58[k]
                       + dh_81[k];

            t_104[k] = -ab_y[k] * dg_59[k]
                       + dh_82[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_y, dg_60, dg_61, dg_62, dg_63, \
                         dg_64, dh_85, dh_87, dh_88, dh_90, dh_91 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_y[k] * dg_60[k]
                       + dh_85[k];

            t_106[k] = -ab_y[k] * dg_61[k]
                       + dh_87[k];

            t_107[k] = -ab_y[k] * dg_62[k]
                       + dh_88[k];

            t_108[k] = -ab_y[k] * dg_63[k]
                       + dh_90[k];

            t_109[k] = -ab_y[k] * dg_64[k]
                       + dh_91[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_y, dg_65, dg_66, dg_67, dg_68, \
                         dg_69, dh_92, dh_94, dh_95, dh_96, dh_97 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = -ab_y[k] * dg_65[k]
                       + dh_92[k];

            t_111[k] = -ab_y[k] * dg_66[k]
                       + dh_94[k];

            t_112[k] = -ab_y[k] * dg_67[k]
                       + dh_95[k];

            t_113[k] = -ab_y[k] * dg_68[k]
                       + dh_96[k];

            t_114[k] = -ab_y[k] * dg_69[k]
                       + dh_97[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_y, dg_70, dg_71, dg_72, dg_73, \
                         dg_74, dh_99, dh_100, dh_101, dh_102, dh_103 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = -ab_y[k] * dg_70[k]
                       + dh_99[k];

            t_116[k] = -ab_y[k] * dg_71[k]
                       + dh_100[k];

            t_117[k] = -ab_y[k] * dg_72[k]
                       + dh_101[k];

            t_118[k] = -ab_y[k] * dg_73[k]
                       + dh_102[k];

            t_119[k] = -ab_y[k] * dg_74[k]
                       + dh_103[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_y, dg_75, dg_76, dg_77, dg_78, \
                         dg_79, dh_106, dh_108, dh_109, dh_111, \
                         dh_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_y[k] * dg_75[k]
                       + dh_106[k];

            t_121[k] = -ab_y[k] * dg_76[k]
                       + dh_108[k];

            t_122[k] = -ab_y[k] * dg_77[k]
                       + dh_109[k];

            t_123[k] = -ab_y[k] * dg_78[k]
                       + dh_111[k];

            t_124[k] = -ab_y[k] * dg_79[k]
                       + dh_112[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_y, dg_80, dg_81, dg_82, dg_83, \
                         dg_84, dh_113, dh_115, dh_116, dh_117, \
                         dh_118 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = -ab_y[k] * dg_80[k]
                       + dh_113[k];

            t_126[k] = -ab_y[k] * dg_81[k]
                       + dh_115[k];

            t_127[k] = -ab_y[k] * dg_82[k]
                       + dh_116[k];

            t_128[k] = -ab_y[k] * dg_83[k]
                       + dh_117[k];

            t_129[k] = -ab_y[k] * dg_84[k]
                       + dh_118[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_y, dg_85, dg_86, dg_87, dg_88, \
                         dg_89, dh_120, dh_121, dh_122, dh_123, \
                         dh_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = -ab_y[k] * dg_85[k]
                       + dh_120[k];

            t_131[k] = -ab_y[k] * dg_86[k]
                       + dh_121[k];

            t_132[k] = -ab_y[k] * dg_87[k]
                       + dh_122[k];

            t_133[k] = -ab_y[k] * dg_88[k]
                       + dh_123[k];

            t_134[k] = -ab_y[k] * dg_89[k]
                       + dh_124[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_z, dg_75, dg_76, dg_77, dg_78, \
                         dg_79, dh_107, dh_109, dh_110, dh_112, \
                         dh_113 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = -ab_z[k] * dg_75[k]
                       + dh_107[k];

            t_136[k] = -ab_z[k] * dg_76[k]
                       + dh_109[k];

            t_137[k] = -ab_z[k] * dg_77[k]
                       + dh_110[k];

            t_138[k] = -ab_z[k] * dg_78[k]
                       + dh_112[k];

            t_139[k] = -ab_z[k] * dg_79[k]
                       + dh_113[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_z, dg_80, dg_81, dg_82, dg_83, \
                         dg_84, dh_114, dh_116, dh_117, dh_118, \
                         dh_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = -ab_z[k] * dg_80[k]
                       + dh_114[k];

            t_141[k] = -ab_z[k] * dg_81[k]
                       + dh_116[k];

            t_142[k] = -ab_z[k] * dg_82[k]
                       + dh_117[k];

            t_143[k] = -ab_z[k] * dg_83[k]
                       + dh_118[k];

            t_144[k] = -ab_z[k] * dg_84[k]
                       + dh_119[k];
        }
    }
}

static auto
compute_hrr_fg_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t dg, const size_t dh, const size_t ncomps,
                      const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_145 = buffer.data(target + 145 * ncomps + c);
        auto *t_146 = buffer.data(target + 146 * ncomps + c);
        auto *t_147 = buffer.data(target + 147 * ncomps + c);
        auto *t_148 = buffer.data(target + 148 * ncomps + c);
        auto *t_149 = buffer.data(target + 149 * ncomps + c);

        const auto *ab_z = coordinates.data(8);

        const auto *dg_85 = buffer.data(dg + 85 * ncomps + c);
        const auto *dg_86 = buffer.data(dg + 86 * ncomps + c);
        const auto *dg_87 = buffer.data(dg + 87 * ncomps + c);
        const auto *dg_88 = buffer.data(dg + 88 * ncomps + c);
        const auto *dg_89 = buffer.data(dg + 89 * ncomps + c);

        const auto *dh_121 = buffer.data(dh + 121 * ncomps + c);
        const auto *dh_122 = buffer.data(dh + 122 * ncomps + c);
        const auto *dh_123 = buffer.data(dh + 123 * ncomps + c);
        const auto *dh_124 = buffer.data(dh + 124 * ncomps + c);
        const auto *dh_125 = buffer.data(dh + 125 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_z, dg_85, dg_86, dg_87, dg_88, \
                         dg_89, dh_121, dh_122, dh_123, dh_124, \
                         dh_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = -ab_z[k] * dg_85[k]
                       + dh_121[k];

            t_146[k] = -ab_z[k] * dg_86[k]
                       + dh_122[k];

            t_147[k] = -ab_z[k] * dg_87[k]
                       + dh_123[k];

            t_148[k] = -ab_z[k] * dg_88[k]
                       + dh_124[k];

            t_149[k] = -ab_z[k] * dg_89[k]
                       + dh_125[k];
        }
    }
}

auto
compute_hrr_fg(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t dg, const size_t dh, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_fg_piece0(buffer, coordinates, target, dg, dh, ncomps, nmax);

    compute_hrr_fg_piece1(buffer, coordinates, target, dg, dh, ncomps, nmax);
}

}  // namespace simdtrf
