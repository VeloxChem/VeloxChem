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


#include "SimdTransferMF.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_mf_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t md, const size_t nd,
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

        const auto *md_0 = buffer.data(md + 0 * ncomps + c);
        const auto *md_1 = buffer.data(md + 1 * ncomps + c);
        const auto *md_2 = buffer.data(md + 2 * ncomps + c);
        const auto *md_3 = buffer.data(md + 3 * ncomps + c);
        const auto *md_4 = buffer.data(md + 4 * ncomps + c);
        const auto *md_5 = buffer.data(md + 5 * ncomps + c);
        const auto *md_6 = buffer.data(md + 6 * ncomps + c);
        const auto *md_7 = buffer.data(md + 7 * ncomps + c);
        const auto *md_8 = buffer.data(md + 8 * ncomps + c);
        const auto *md_9 = buffer.data(md + 9 * ncomps + c);
        const auto *md_10 = buffer.data(md + 10 * ncomps + c);
        const auto *md_11 = buffer.data(md + 11 * ncomps + c);
        const auto *md_12 = buffer.data(md + 12 * ncomps + c);
        const auto *md_13 = buffer.data(md + 13 * ncomps + c);
        const auto *md_14 = buffer.data(md + 14 * ncomps + c);
        const auto *md_15 = buffer.data(md + 15 * ncomps + c);
        const auto *md_16 = buffer.data(md + 16 * ncomps + c);
        const auto *md_17 = buffer.data(md + 17 * ncomps + c);
        const auto *md_18 = buffer.data(md + 18 * ncomps + c);
        const auto *md_19 = buffer.data(md + 19 * ncomps + c);
        const auto *md_20 = buffer.data(md + 20 * ncomps + c);
        const auto *md_21 = buffer.data(md + 21 * ncomps + c);
        const auto *md_22 = buffer.data(md + 22 * ncomps + c);
        const auto *md_23 = buffer.data(md + 23 * ncomps + c);
        const auto *md_24 = buffer.data(md + 24 * ncomps + c);
        const auto *md_25 = buffer.data(md + 25 * ncomps + c);
        const auto *md_26 = buffer.data(md + 26 * ncomps + c);
        const auto *md_27 = buffer.data(md + 27 * ncomps + c);
        const auto *md_28 = buffer.data(md + 28 * ncomps + c);
        const auto *md_29 = buffer.data(md + 29 * ncomps + c);
        const auto *md_30 = buffer.data(md + 30 * ncomps + c);
        const auto *md_31 = buffer.data(md + 31 * ncomps + c);
        const auto *md_32 = buffer.data(md + 32 * ncomps + c);
        const auto *md_33 = buffer.data(md + 33 * ncomps + c);
        const auto *md_34 = buffer.data(md + 34 * ncomps + c);
        const auto *md_35 = buffer.data(md + 35 * ncomps + c);
        const auto *md_36 = buffer.data(md + 36 * ncomps + c);
        const auto *md_37 = buffer.data(md + 37 * ncomps + c);
        const auto *md_38 = buffer.data(md + 38 * ncomps + c);
        const auto *md_39 = buffer.data(md + 39 * ncomps + c);
        const auto *md_40 = buffer.data(md + 40 * ncomps + c);
        const auto *md_41 = buffer.data(md + 41 * ncomps + c);
        const auto *md_42 = buffer.data(md + 42 * ncomps + c);
        const auto *md_43 = buffer.data(md + 43 * ncomps + c);
        const auto *md_44 = buffer.data(md + 44 * ncomps + c);
        const auto *md_45 = buffer.data(md + 45 * ncomps + c);
        const auto *md_46 = buffer.data(md + 46 * ncomps + c);
        const auto *md_47 = buffer.data(md + 47 * ncomps + c);
        const auto *md_48 = buffer.data(md + 48 * ncomps + c);
        const auto *md_49 = buffer.data(md + 49 * ncomps + c);
        const auto *md_50 = buffer.data(md + 50 * ncomps + c);
        const auto *md_51 = buffer.data(md + 51 * ncomps + c);
        const auto *md_52 = buffer.data(md + 52 * ncomps + c);
        const auto *md_53 = buffer.data(md + 53 * ncomps + c);
        const auto *md_54 = buffer.data(md + 54 * ncomps + c);
        const auto *md_55 = buffer.data(md + 55 * ncomps + c);
        const auto *md_56 = buffer.data(md + 56 * ncomps + c);
        const auto *md_57 = buffer.data(md + 57 * ncomps + c);
        const auto *md_58 = buffer.data(md + 58 * ncomps + c);
        const auto *md_59 = buffer.data(md + 59 * ncomps + c);
        const auto *md_60 = buffer.data(md + 60 * ncomps + c);
        const auto *md_61 = buffer.data(md + 61 * ncomps + c);
        const auto *md_62 = buffer.data(md + 62 * ncomps + c);
        const auto *md_63 = buffer.data(md + 63 * ncomps + c);
        const auto *md_64 = buffer.data(md + 64 * ncomps + c);
        const auto *md_65 = buffer.data(md + 65 * ncomps + c);
        const auto *md_66 = buffer.data(md + 66 * ncomps + c);
        const auto *md_67 = buffer.data(md + 67 * ncomps + c);
        const auto *md_68 = buffer.data(md + 68 * ncomps + c);
        const auto *md_69 = buffer.data(md + 69 * ncomps + c);
        const auto *md_70 = buffer.data(md + 70 * ncomps + c);
        const auto *md_71 = buffer.data(md + 71 * ncomps + c);
        const auto *md_72 = buffer.data(md + 72 * ncomps + c);
        const auto *md_73 = buffer.data(md + 73 * ncomps + c);
        const auto *md_74 = buffer.data(md + 74 * ncomps + c);
        const auto *md_75 = buffer.data(md + 75 * ncomps + c);
        const auto *md_76 = buffer.data(md + 76 * ncomps + c);
        const auto *md_77 = buffer.data(md + 77 * ncomps + c);
        const auto *md_78 = buffer.data(md + 78 * ncomps + c);
        const auto *md_79 = buffer.data(md + 79 * ncomps + c);
        const auto *md_80 = buffer.data(md + 80 * ncomps + c);
        const auto *md_81 = buffer.data(md + 81 * ncomps + c);
        const auto *md_82 = buffer.data(md + 82 * ncomps + c);
        const auto *md_83 = buffer.data(md + 83 * ncomps + c);
        const auto *md_84 = buffer.data(md + 84 * ncomps + c);
        const auto *md_85 = buffer.data(md + 85 * ncomps + c);
        const auto *md_86 = buffer.data(md + 86 * ncomps + c);
        const auto *md_87 = buffer.data(md + 87 * ncomps + c);
        const auto *md_88 = buffer.data(md + 88 * ncomps + c);

        const auto *nd_0 = buffer.data(nd + 0 * ncomps + c);
        const auto *nd_1 = buffer.data(nd + 1 * ncomps + c);
        const auto *nd_2 = buffer.data(nd + 2 * ncomps + c);
        const auto *nd_3 = buffer.data(nd + 3 * ncomps + c);
        const auto *nd_4 = buffer.data(nd + 4 * ncomps + c);
        const auto *nd_5 = buffer.data(nd + 5 * ncomps + c);
        const auto *nd_6 = buffer.data(nd + 6 * ncomps + c);
        const auto *nd_7 = buffer.data(nd + 7 * ncomps + c);
        const auto *nd_8 = buffer.data(nd + 8 * ncomps + c);
        const auto *nd_9 = buffer.data(nd + 9 * ncomps + c);
        const auto *nd_10 = buffer.data(nd + 10 * ncomps + c);
        const auto *nd_11 = buffer.data(nd + 11 * ncomps + c);
        const auto *nd_12 = buffer.data(nd + 12 * ncomps + c);
        const auto *nd_13 = buffer.data(nd + 13 * ncomps + c);
        const auto *nd_14 = buffer.data(nd + 14 * ncomps + c);
        const auto *nd_15 = buffer.data(nd + 15 * ncomps + c);
        const auto *nd_16 = buffer.data(nd + 16 * ncomps + c);
        const auto *nd_17 = buffer.data(nd + 17 * ncomps + c);
        const auto *nd_18 = buffer.data(nd + 18 * ncomps + c);
        const auto *nd_19 = buffer.data(nd + 19 * ncomps + c);
        const auto *nd_20 = buffer.data(nd + 20 * ncomps + c);
        const auto *nd_21 = buffer.data(nd + 21 * ncomps + c);
        const auto *nd_22 = buffer.data(nd + 22 * ncomps + c);
        const auto *nd_23 = buffer.data(nd + 23 * ncomps + c);
        const auto *nd_24 = buffer.data(nd + 24 * ncomps + c);
        const auto *nd_25 = buffer.data(nd + 25 * ncomps + c);
        const auto *nd_26 = buffer.data(nd + 26 * ncomps + c);
        const auto *nd_27 = buffer.data(nd + 27 * ncomps + c);
        const auto *nd_28 = buffer.data(nd + 28 * ncomps + c);
        const auto *nd_29 = buffer.data(nd + 29 * ncomps + c);
        const auto *nd_30 = buffer.data(nd + 30 * ncomps + c);
        const auto *nd_31 = buffer.data(nd + 31 * ncomps + c);
        const auto *nd_32 = buffer.data(nd + 32 * ncomps + c);
        const auto *nd_33 = buffer.data(nd + 33 * ncomps + c);
        const auto *nd_34 = buffer.data(nd + 34 * ncomps + c);
        const auto *nd_35 = buffer.data(nd + 35 * ncomps + c);
        const auto *nd_36 = buffer.data(nd + 36 * ncomps + c);
        const auto *nd_37 = buffer.data(nd + 37 * ncomps + c);
        const auto *nd_38 = buffer.data(nd + 38 * ncomps + c);
        const auto *nd_39 = buffer.data(nd + 39 * ncomps + c);
        const auto *nd_40 = buffer.data(nd + 40 * ncomps + c);
        const auto *nd_41 = buffer.data(nd + 41 * ncomps + c);
        const auto *nd_42 = buffer.data(nd + 42 * ncomps + c);
        const auto *nd_43 = buffer.data(nd + 43 * ncomps + c);
        const auto *nd_44 = buffer.data(nd + 44 * ncomps + c);
        const auto *nd_45 = buffer.data(nd + 45 * ncomps + c);
        const auto *nd_46 = buffer.data(nd + 46 * ncomps + c);
        const auto *nd_47 = buffer.data(nd + 47 * ncomps + c);
        const auto *nd_48 = buffer.data(nd + 48 * ncomps + c);
        const auto *nd_49 = buffer.data(nd + 49 * ncomps + c);
        const auto *nd_50 = buffer.data(nd + 50 * ncomps + c);
        const auto *nd_51 = buffer.data(nd + 51 * ncomps + c);
        const auto *nd_52 = buffer.data(nd + 52 * ncomps + c);
        const auto *nd_53 = buffer.data(nd + 53 * ncomps + c);
        const auto *nd_54 = buffer.data(nd + 54 * ncomps + c);
        const auto *nd_55 = buffer.data(nd + 55 * ncomps + c);
        const auto *nd_56 = buffer.data(nd + 56 * ncomps + c);
        const auto *nd_57 = buffer.data(nd + 57 * ncomps + c);
        const auto *nd_58 = buffer.data(nd + 58 * ncomps + c);
        const auto *nd_59 = buffer.data(nd + 59 * ncomps + c);
        const auto *nd_60 = buffer.data(nd + 60 * ncomps + c);
        const auto *nd_61 = buffer.data(nd + 61 * ncomps + c);
        const auto *nd_62 = buffer.data(nd + 62 * ncomps + c);
        const auto *nd_63 = buffer.data(nd + 63 * ncomps + c);
        const auto *nd_64 = buffer.data(nd + 64 * ncomps + c);
        const auto *nd_65 = buffer.data(nd + 65 * ncomps + c);
        const auto *nd_66 = buffer.data(nd + 66 * ncomps + c);
        const auto *nd_67 = buffer.data(nd + 67 * ncomps + c);
        const auto *nd_68 = buffer.data(nd + 68 * ncomps + c);
        const auto *nd_69 = buffer.data(nd + 69 * ncomps + c);
        const auto *nd_70 = buffer.data(nd + 70 * ncomps + c);
        const auto *nd_71 = buffer.data(nd + 71 * ncomps + c);
        const auto *nd_72 = buffer.data(nd + 72 * ncomps + c);
        const auto *nd_73 = buffer.data(nd + 73 * ncomps + c);
        const auto *nd_74 = buffer.data(nd + 74 * ncomps + c);
        const auto *nd_75 = buffer.data(nd + 75 * ncomps + c);
        const auto *nd_76 = buffer.data(nd + 76 * ncomps + c);
        const auto *nd_77 = buffer.data(nd + 77 * ncomps + c);
        const auto *nd_78 = buffer.data(nd + 78 * ncomps + c);
        const auto *nd_79 = buffer.data(nd + 79 * ncomps + c);
        const auto *nd_80 = buffer.data(nd + 80 * ncomps + c);
        const auto *nd_81 = buffer.data(nd + 81 * ncomps + c);
        const auto *nd_82 = buffer.data(nd + 82 * ncomps + c);
        const auto *nd_83 = buffer.data(nd + 83 * ncomps + c);
        const auto *nd_84 = buffer.data(nd + 84 * ncomps + c);
        const auto *nd_85 = buffer.data(nd + 85 * ncomps + c);
        const auto *nd_86 = buffer.data(nd + 86 * ncomps + c);
        const auto *nd_87 = buffer.data(nd + 87 * ncomps + c);
        const auto *nd_88 = buffer.data(nd + 88 * ncomps + c);
        const auto *nd_89 = buffer.data(nd + 89 * ncomps + c);
        const auto *nd_93 = buffer.data(nd + 93 * ncomps + c);
        const auto *nd_94 = buffer.data(nd + 94 * ncomps + c);
        const auto *nd_95 = buffer.data(nd + 95 * ncomps + c);
        const auto *nd_99 = buffer.data(nd + 99 * ncomps + c);
        const auto *nd_100 = buffer.data(nd + 100 * ncomps + c);
        const auto *nd_101 = buffer.data(nd + 101 * ncomps + c);
        const auto *nd_105 = buffer.data(nd + 105 * ncomps + c);
        const auto *nd_106 = buffer.data(nd + 106 * ncomps + c);
        const auto *nd_107 = buffer.data(nd + 107 * ncomps + c);
        const auto *nd_111 = buffer.data(nd + 111 * ncomps + c);
        const auto *nd_112 = buffer.data(nd + 112 * ncomps + c);
        const auto *nd_113 = buffer.data(nd + 113 * ncomps + c);
        const auto *nd_119 = buffer.data(nd + 119 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, md_0, md_1, md_2, md_3, md_4, nd_0, \
                         nd_1, nd_2, nd_3, nd_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * md_0[k]
                     + nd_0[k];

            t_1[k] = ab_x[k] * md_1[k]
                     + nd_1[k];

            t_2[k] = ab_x[k] * md_2[k]
                     + nd_2[k];

            t_3[k] = ab_x[k] * md_3[k]
                     + nd_3[k];

            t_4[k] = ab_x[k] * md_4[k]
                     + nd_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ab_y, ab_z, md_3, md_4, md_5, nd_5, \
                         nd_9, nd_10, nd_11, nd_17 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * md_5[k]
                     + nd_5[k];

            t_6[k] = ab_y[k] * md_3[k]
                     + nd_9[k];

            t_7[k] = ab_y[k] * md_4[k]
                     + nd_10[k];

            t_8[k] = ab_y[k] * md_5[k]
                     + nd_11[k];

            t_9[k] = ab_z[k] * md_5[k]
                     + nd_17[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, md_6, md_7, md_8, md_9, md_10, \
                         nd_6, nd_7, nd_8, nd_9, nd_10 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_x[k] * md_6[k]
                      + nd_6[k];

            t_11[k] = ab_x[k] * md_7[k]
                      + nd_7[k];

            t_12[k] = ab_x[k] * md_8[k]
                      + nd_8[k];

            t_13[k] = ab_x[k] * md_9[k]
                      + nd_9[k];

            t_14[k] = ab_x[k] * md_10[k]
                      + nd_10[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, ab_y, ab_z, md_9, md_10, md_11, \
                         nd_11, nd_21, nd_22, nd_23, nd_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_x[k] * md_11[k]
                      + nd_11[k];

            t_16[k] = ab_y[k] * md_9[k]
                      + nd_21[k];

            t_17[k] = ab_y[k] * md_10[k]
                      + nd_22[k];

            t_18[k] = ab_y[k] * md_11[k]
                      + nd_23[k];

            t_19[k] = ab_z[k] * md_11[k]
                      + nd_29[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, md_12, md_13, md_14, md_15, \
                         md_16, nd_12, nd_13, nd_14, nd_15, nd_16 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * md_12[k]
                      + nd_12[k];

            t_21[k] = ab_x[k] * md_13[k]
                      + nd_13[k];

            t_22[k] = ab_x[k] * md_14[k]
                      + nd_14[k];

            t_23[k] = ab_x[k] * md_15[k]
                      + nd_15[k];

            t_24[k] = ab_x[k] * md_16[k]
                      + nd_16[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, ab_y, ab_z, md_15, md_16, md_17, \
                         nd_17, nd_27, nd_28, nd_29, nd_35 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_x[k] * md_17[k]
                      + nd_17[k];

            t_26[k] = ab_y[k] * md_15[k]
                      + nd_27[k];

            t_27[k] = ab_y[k] * md_16[k]
                      + nd_28[k];

            t_28[k] = ab_y[k] * md_17[k]
                      + nd_29[k];

            t_29[k] = ab_z[k] * md_17[k]
                      + nd_35[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, md_18, md_19, md_20, md_21, \
                         md_22, nd_18, nd_19, nd_20, nd_21, nd_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = ab_x[k] * md_18[k]
                      + nd_18[k];

            t_31[k] = ab_x[k] * md_19[k]
                      + nd_19[k];

            t_32[k] = ab_x[k] * md_20[k]
                      + nd_20[k];

            t_33[k] = ab_x[k] * md_21[k]
                      + nd_21[k];

            t_34[k] = ab_x[k] * md_22[k]
                      + nd_22[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, ab_y, ab_z, md_21, md_22, md_23, \
                         nd_23, nd_39, nd_40, nd_41, nd_47 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * md_23[k]
                      + nd_23[k];

            t_36[k] = ab_y[k] * md_21[k]
                      + nd_39[k];

            t_37[k] = ab_y[k] * md_22[k]
                      + nd_40[k];

            t_38[k] = ab_y[k] * md_23[k]
                      + nd_41[k];

            t_39[k] = ab_z[k] * md_23[k]
                      + nd_47[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, md_24, md_25, md_26, md_27, \
                         md_28, nd_24, nd_25, nd_26, nd_27, nd_28 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_x[k] * md_24[k]
                      + nd_24[k];

            t_41[k] = ab_x[k] * md_25[k]
                      + nd_25[k];

            t_42[k] = ab_x[k] * md_26[k]
                      + nd_26[k];

            t_43[k] = ab_x[k] * md_27[k]
                      + nd_27[k];

            t_44[k] = ab_x[k] * md_28[k]
                      + nd_28[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, ab_y, ab_z, md_27, md_28, md_29, \
                         nd_29, nd_45, nd_46, nd_47, nd_53 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * md_29[k]
                      + nd_29[k];

            t_46[k] = ab_y[k] * md_27[k]
                      + nd_45[k];

            t_47[k] = ab_y[k] * md_28[k]
                      + nd_46[k];

            t_48[k] = ab_y[k] * md_29[k]
                      + nd_47[k];

            t_49[k] = ab_z[k] * md_29[k]
                      + nd_53[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, md_30, md_31, md_32, md_33, \
                         md_34, nd_30, nd_31, nd_32, nd_33, nd_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = ab_x[k] * md_30[k]
                      + nd_30[k];

            t_51[k] = ab_x[k] * md_31[k]
                      + nd_31[k];

            t_52[k] = ab_x[k] * md_32[k]
                      + nd_32[k];

            t_53[k] = ab_x[k] * md_33[k]
                      + nd_33[k];

            t_54[k] = ab_x[k] * md_34[k]
                      + nd_34[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, md_33, md_34, md_35, \
                         nd_35, nd_51, nd_52, nd_53, nd_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = ab_x[k] * md_35[k]
                      + nd_35[k];

            t_56[k] = ab_y[k] * md_33[k]
                      + nd_51[k];

            t_57[k] = ab_y[k] * md_34[k]
                      + nd_52[k];

            t_58[k] = ab_y[k] * md_35[k]
                      + nd_53[k];

            t_59[k] = ab_z[k] * md_35[k]
                      + nd_59[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, md_36, md_37, md_38, md_39, \
                         md_40, nd_36, nd_37, nd_38, nd_39, nd_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * md_36[k]
                      + nd_36[k];

            t_61[k] = ab_x[k] * md_37[k]
                      + nd_37[k];

            t_62[k] = ab_x[k] * md_38[k]
                      + nd_38[k];

            t_63[k] = ab_x[k] * md_39[k]
                      + nd_39[k];

            t_64[k] = ab_x[k] * md_40[k]
                      + nd_40[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, ab_y, ab_z, md_39, md_40, md_41, \
                         nd_41, nd_63, nd_64, nd_65, nd_71 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * md_41[k]
                      + nd_41[k];

            t_66[k] = ab_y[k] * md_39[k]
                      + nd_63[k];

            t_67[k] = ab_y[k] * md_40[k]
                      + nd_64[k];

            t_68[k] = ab_y[k] * md_41[k]
                      + nd_65[k];

            t_69[k] = ab_z[k] * md_41[k]
                      + nd_71[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, md_42, md_43, md_44, md_45, \
                         md_46, nd_42, nd_43, nd_44, nd_45, nd_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = ab_x[k] * md_42[k]
                      + nd_42[k];

            t_71[k] = ab_x[k] * md_43[k]
                      + nd_43[k];

            t_72[k] = ab_x[k] * md_44[k]
                      + nd_44[k];

            t_73[k] = ab_x[k] * md_45[k]
                      + nd_45[k];

            t_74[k] = ab_x[k] * md_46[k]
                      + nd_46[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, ab_y, ab_z, md_45, md_46, md_47, \
                         nd_47, nd_69, nd_70, nd_71, nd_77 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = ab_x[k] * md_47[k]
                      + nd_47[k];

            t_76[k] = ab_y[k] * md_45[k]
                      + nd_69[k];

            t_77[k] = ab_y[k] * md_46[k]
                      + nd_70[k];

            t_78[k] = ab_y[k] * md_47[k]
                      + nd_71[k];

            t_79[k] = ab_z[k] * md_47[k]
                      + nd_77[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, md_48, md_49, md_50, md_51, \
                         md_52, nd_48, nd_49, nd_50, nd_51, nd_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_x[k] * md_48[k]
                      + nd_48[k];

            t_81[k] = ab_x[k] * md_49[k]
                      + nd_49[k];

            t_82[k] = ab_x[k] * md_50[k]
                      + nd_50[k];

            t_83[k] = ab_x[k] * md_51[k]
                      + nd_51[k];

            t_84[k] = ab_x[k] * md_52[k]
                      + nd_52[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, md_51, md_52, md_53, \
                         nd_53, nd_75, nd_76, nd_77, nd_83 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_x[k] * md_53[k]
                      + nd_53[k];

            t_86[k] = ab_y[k] * md_51[k]
                      + nd_75[k];

            t_87[k] = ab_y[k] * md_52[k]
                      + nd_76[k];

            t_88[k] = ab_y[k] * md_53[k]
                      + nd_77[k];

            t_89[k] = ab_z[k] * md_53[k]
                      + nd_83[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, md_54, md_55, md_56, md_57, \
                         md_58, nd_54, nd_55, nd_56, nd_57, nd_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * md_54[k]
                      + nd_54[k];

            t_91[k] = ab_x[k] * md_55[k]
                      + nd_55[k];

            t_92[k] = ab_x[k] * md_56[k]
                      + nd_56[k];

            t_93[k] = ab_x[k] * md_57[k]
                      + nd_57[k];

            t_94[k] = ab_x[k] * md_58[k]
                      + nd_58[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, ab_y, ab_z, md_57, md_58, md_59, \
                         nd_59, nd_81, nd_82, nd_83, nd_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_x[k] * md_59[k]
                      + nd_59[k];

            t_96[k] = ab_y[k] * md_57[k]
                      + nd_81[k];

            t_97[k] = ab_y[k] * md_58[k]
                      + nd_82[k];

            t_98[k] = ab_y[k] * md_59[k]
                      + nd_83[k];

            t_99[k] = ab_z[k] * md_59[k]
                      + nd_89[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, md_60, md_61, md_62, md_63, \
                         md_64, nd_60, nd_61, nd_62, nd_63, nd_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = ab_x[k] * md_60[k]
                       + nd_60[k];

            t_101[k] = ab_x[k] * md_61[k]
                       + nd_61[k];

            t_102[k] = ab_x[k] * md_62[k]
                       + nd_62[k];

            t_103[k] = ab_x[k] * md_63[k]
                       + nd_63[k];

            t_104[k] = ab_x[k] * md_64[k]
                       + nd_64[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, ab_y, ab_z, md_63, md_64, \
                         md_65, nd_65, nd_93, nd_94, nd_95, nd_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * md_65[k]
                       + nd_65[k];

            t_106[k] = ab_y[k] * md_63[k]
                       + nd_93[k];

            t_107[k] = ab_y[k] * md_64[k]
                       + nd_94[k];

            t_108[k] = ab_y[k] * md_65[k]
                       + nd_95[k];

            t_109[k] = ab_z[k] * md_65[k]
                       + nd_101[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, md_66, md_67, md_68, md_69, \
                         md_70, nd_66, nd_67, nd_68, nd_69, nd_70 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_x[k] * md_66[k]
                       + nd_66[k];

            t_111[k] = ab_x[k] * md_67[k]
                       + nd_67[k];

            t_112[k] = ab_x[k] * md_68[k]
                       + nd_68[k];

            t_113[k] = ab_x[k] * md_69[k]
                       + nd_69[k];

            t_114[k] = ab_x[k] * md_70[k]
                       + nd_70[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, ab_y, ab_z, md_69, md_70, \
                         md_71, nd_71, nd_99, nd_100, nd_101, nd_107 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = ab_x[k] * md_71[k]
                       + nd_71[k];

            t_116[k] = ab_y[k] * md_69[k]
                       + nd_99[k];

            t_117[k] = ab_y[k] * md_70[k]
                       + nd_100[k];

            t_118[k] = ab_y[k] * md_71[k]
                       + nd_101[k];

            t_119[k] = ab_z[k] * md_71[k]
                       + nd_107[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, md_72, md_73, md_74, md_75, \
                         md_76, nd_72, nd_73, nd_74, nd_75, nd_76 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = ab_x[k] * md_72[k]
                       + nd_72[k];

            t_121[k] = ab_x[k] * md_73[k]
                       + nd_73[k];

            t_122[k] = ab_x[k] * md_74[k]
                       + nd_74[k];

            t_123[k] = ab_x[k] * md_75[k]
                       + nd_75[k];

            t_124[k] = ab_x[k] * md_76[k]
                       + nd_76[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, ab_y, ab_z, md_75, md_76, \
                         md_77, nd_77, nd_105, nd_106, nd_107, nd_113 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = ab_x[k] * md_77[k]
                       + nd_77[k];

            t_126[k] = ab_y[k] * md_75[k]
                       + nd_105[k];

            t_127[k] = ab_y[k] * md_76[k]
                       + nd_106[k];

            t_128[k] = ab_y[k] * md_77[k]
                       + nd_107[k];

            t_129[k] = ab_z[k] * md_77[k]
                       + nd_113[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, md_78, md_79, md_80, md_81, \
                         md_82, nd_78, nd_79, nd_80, nd_81, nd_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = ab_x[k] * md_78[k]
                       + nd_78[k];

            t_131[k] = ab_x[k] * md_79[k]
                       + nd_79[k];

            t_132[k] = ab_x[k] * md_80[k]
                       + nd_80[k];

            t_133[k] = ab_x[k] * md_81[k]
                       + nd_81[k];

            t_134[k] = ab_x[k] * md_82[k]
                       + nd_82[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, ab_y, ab_z, md_81, md_82, \
                         md_83, nd_83, nd_111, nd_112, nd_113, nd_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = ab_x[k] * md_83[k]
                       + nd_83[k];

            t_136[k] = ab_y[k] * md_81[k]
                       + nd_111[k];

            t_137[k] = ab_y[k] * md_82[k]
                       + nd_112[k];

            t_138[k] = ab_y[k] * md_83[k]
                       + nd_113[k];

            t_139[k] = ab_z[k] * md_83[k]
                       + nd_119[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, md_84, md_85, md_86, md_87, \
                         md_88, nd_84, nd_85, nd_86, nd_87, nd_88 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = ab_x[k] * md_84[k]
                       + nd_84[k];

            t_141[k] = ab_x[k] * md_85[k]
                       + nd_85[k];

            t_142[k] = ab_x[k] * md_86[k]
                       + nd_86[k];

            t_143[k] = ab_x[k] * md_87[k]
                       + nd_87[k];

            t_144[k] = ab_x[k] * md_88[k]
                       + nd_88[k];
        }
    }
}

static auto
compute_hrr_mf_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t md, const size_t nd,
                                   const size_t ncomps, const size_t nmax) -> void
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
        auto *t_210 = buffer.data(target + 210 * ncomps + c);
        auto *t_211 = buffer.data(target + 211 * ncomps + c);
        auto *t_212 = buffer.data(target + 212 * ncomps + c);
        auto *t_213 = buffer.data(target + 213 * ncomps + c);
        auto *t_214 = buffer.data(target + 214 * ncomps + c);
        auto *t_215 = buffer.data(target + 215 * ncomps + c);
        auto *t_216 = buffer.data(target + 216 * ncomps + c);
        auto *t_217 = buffer.data(target + 217 * ncomps + c);
        auto *t_218 = buffer.data(target + 218 * ncomps + c);
        auto *t_219 = buffer.data(target + 219 * ncomps + c);
        auto *t_220 = buffer.data(target + 220 * ncomps + c);
        auto *t_221 = buffer.data(target + 221 * ncomps + c);
        auto *t_222 = buffer.data(target + 222 * ncomps + c);
        auto *t_223 = buffer.data(target + 223 * ncomps + c);
        auto *t_224 = buffer.data(target + 224 * ncomps + c);
        auto *t_225 = buffer.data(target + 225 * ncomps + c);
        auto *t_226 = buffer.data(target + 226 * ncomps + c);
        auto *t_227 = buffer.data(target + 227 * ncomps + c);
        auto *t_228 = buffer.data(target + 228 * ncomps + c);
        auto *t_229 = buffer.data(target + 229 * ncomps + c);
        auto *t_230 = buffer.data(target + 230 * ncomps + c);
        auto *t_231 = buffer.data(target + 231 * ncomps + c);
        auto *t_232 = buffer.data(target + 232 * ncomps + c);
        auto *t_233 = buffer.data(target + 233 * ncomps + c);
        auto *t_234 = buffer.data(target + 234 * ncomps + c);
        auto *t_235 = buffer.data(target + 235 * ncomps + c);
        auto *t_236 = buffer.data(target + 236 * ncomps + c);
        auto *t_237 = buffer.data(target + 237 * ncomps + c);
        auto *t_238 = buffer.data(target + 238 * ncomps + c);
        auto *t_239 = buffer.data(target + 239 * ncomps + c);
        auto *t_240 = buffer.data(target + 240 * ncomps + c);
        auto *t_241 = buffer.data(target + 241 * ncomps + c);
        auto *t_242 = buffer.data(target + 242 * ncomps + c);
        auto *t_243 = buffer.data(target + 243 * ncomps + c);
        auto *t_244 = buffer.data(target + 244 * ncomps + c);
        auto *t_245 = buffer.data(target + 245 * ncomps + c);
        auto *t_246 = buffer.data(target + 246 * ncomps + c);
        auto *t_247 = buffer.data(target + 247 * ncomps + c);
        auto *t_248 = buffer.data(target + 248 * ncomps + c);
        auto *t_249 = buffer.data(target + 249 * ncomps + c);
        auto *t_250 = buffer.data(target + 250 * ncomps + c);
        auto *t_251 = buffer.data(target + 251 * ncomps + c);
        auto *t_252 = buffer.data(target + 252 * ncomps + c);
        auto *t_253 = buffer.data(target + 253 * ncomps + c);
        auto *t_254 = buffer.data(target + 254 * ncomps + c);
        auto *t_255 = buffer.data(target + 255 * ncomps + c);
        auto *t_256 = buffer.data(target + 256 * ncomps + c);
        auto *t_257 = buffer.data(target + 257 * ncomps + c);
        auto *t_258 = buffer.data(target + 258 * ncomps + c);
        auto *t_259 = buffer.data(target + 259 * ncomps + c);
        auto *t_260 = buffer.data(target + 260 * ncomps + c);
        auto *t_261 = buffer.data(target + 261 * ncomps + c);
        auto *t_262 = buffer.data(target + 262 * ncomps + c);
        auto *t_263 = buffer.data(target + 263 * ncomps + c);
        auto *t_264 = buffer.data(target + 264 * ncomps + c);
        auto *t_265 = buffer.data(target + 265 * ncomps + c);
        auto *t_266 = buffer.data(target + 266 * ncomps + c);
        auto *t_267 = buffer.data(target + 267 * ncomps + c);
        auto *t_268 = buffer.data(target + 268 * ncomps + c);
        auto *t_269 = buffer.data(target + 269 * ncomps + c);
        auto *t_270 = buffer.data(target + 270 * ncomps + c);
        auto *t_271 = buffer.data(target + 271 * ncomps + c);
        auto *t_272 = buffer.data(target + 272 * ncomps + c);
        auto *t_273 = buffer.data(target + 273 * ncomps + c);
        auto *t_274 = buffer.data(target + 274 * ncomps + c);
        auto *t_275 = buffer.data(target + 275 * ncomps + c);
        auto *t_276 = buffer.data(target + 276 * ncomps + c);
        auto *t_277 = buffer.data(target + 277 * ncomps + c);
        auto *t_278 = buffer.data(target + 278 * ncomps + c);
        auto *t_279 = buffer.data(target + 279 * ncomps + c);
        auto *t_280 = buffer.data(target + 280 * ncomps + c);
        auto *t_281 = buffer.data(target + 281 * ncomps + c);
        auto *t_282 = buffer.data(target + 282 * ncomps + c);
        auto *t_283 = buffer.data(target + 283 * ncomps + c);
        auto *t_284 = buffer.data(target + 284 * ncomps + c);
        auto *t_285 = buffer.data(target + 285 * ncomps + c);
        auto *t_286 = buffer.data(target + 286 * ncomps + c);
        auto *t_287 = buffer.data(target + 287 * ncomps + c);
        auto *t_288 = buffer.data(target + 288 * ncomps + c);
        auto *t_289 = buffer.data(target + 289 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *md_87 = buffer.data(md + 87 * ncomps + c);
        const auto *md_88 = buffer.data(md + 88 * ncomps + c);
        const auto *md_89 = buffer.data(md + 89 * ncomps + c);
        const auto *md_90 = buffer.data(md + 90 * ncomps + c);
        const auto *md_91 = buffer.data(md + 91 * ncomps + c);
        const auto *md_92 = buffer.data(md + 92 * ncomps + c);
        const auto *md_93 = buffer.data(md + 93 * ncomps + c);
        const auto *md_94 = buffer.data(md + 94 * ncomps + c);
        const auto *md_95 = buffer.data(md + 95 * ncomps + c);
        const auto *md_96 = buffer.data(md + 96 * ncomps + c);
        const auto *md_97 = buffer.data(md + 97 * ncomps + c);
        const auto *md_98 = buffer.data(md + 98 * ncomps + c);
        const auto *md_99 = buffer.data(md + 99 * ncomps + c);
        const auto *md_100 = buffer.data(md + 100 * ncomps + c);
        const auto *md_101 = buffer.data(md + 101 * ncomps + c);
        const auto *md_102 = buffer.data(md + 102 * ncomps + c);
        const auto *md_103 = buffer.data(md + 103 * ncomps + c);
        const auto *md_104 = buffer.data(md + 104 * ncomps + c);
        const auto *md_105 = buffer.data(md + 105 * ncomps + c);
        const auto *md_106 = buffer.data(md + 106 * ncomps + c);
        const auto *md_107 = buffer.data(md + 107 * ncomps + c);
        const auto *md_108 = buffer.data(md + 108 * ncomps + c);
        const auto *md_109 = buffer.data(md + 109 * ncomps + c);
        const auto *md_110 = buffer.data(md + 110 * ncomps + c);
        const auto *md_111 = buffer.data(md + 111 * ncomps + c);
        const auto *md_112 = buffer.data(md + 112 * ncomps + c);
        const auto *md_113 = buffer.data(md + 113 * ncomps + c);
        const auto *md_114 = buffer.data(md + 114 * ncomps + c);
        const auto *md_115 = buffer.data(md + 115 * ncomps + c);
        const auto *md_116 = buffer.data(md + 116 * ncomps + c);
        const auto *md_117 = buffer.data(md + 117 * ncomps + c);
        const auto *md_118 = buffer.data(md + 118 * ncomps + c);
        const auto *md_119 = buffer.data(md + 119 * ncomps + c);
        const auto *md_120 = buffer.data(md + 120 * ncomps + c);
        const auto *md_121 = buffer.data(md + 121 * ncomps + c);
        const auto *md_122 = buffer.data(md + 122 * ncomps + c);
        const auto *md_123 = buffer.data(md + 123 * ncomps + c);
        const auto *md_124 = buffer.data(md + 124 * ncomps + c);
        const auto *md_125 = buffer.data(md + 125 * ncomps + c);
        const auto *md_126 = buffer.data(md + 126 * ncomps + c);
        const auto *md_127 = buffer.data(md + 127 * ncomps + c);
        const auto *md_128 = buffer.data(md + 128 * ncomps + c);
        const auto *md_129 = buffer.data(md + 129 * ncomps + c);
        const auto *md_130 = buffer.data(md + 130 * ncomps + c);
        const auto *md_131 = buffer.data(md + 131 * ncomps + c);
        const auto *md_132 = buffer.data(md + 132 * ncomps + c);
        const auto *md_133 = buffer.data(md + 133 * ncomps + c);
        const auto *md_134 = buffer.data(md + 134 * ncomps + c);
        const auto *md_135 = buffer.data(md + 135 * ncomps + c);
        const auto *md_136 = buffer.data(md + 136 * ncomps + c);
        const auto *md_137 = buffer.data(md + 137 * ncomps + c);
        const auto *md_138 = buffer.data(md + 138 * ncomps + c);
        const auto *md_139 = buffer.data(md + 139 * ncomps + c);
        const auto *md_140 = buffer.data(md + 140 * ncomps + c);
        const auto *md_141 = buffer.data(md + 141 * ncomps + c);
        const auto *md_142 = buffer.data(md + 142 * ncomps + c);
        const auto *md_143 = buffer.data(md + 143 * ncomps + c);
        const auto *md_144 = buffer.data(md + 144 * ncomps + c);
        const auto *md_145 = buffer.data(md + 145 * ncomps + c);
        const auto *md_146 = buffer.data(md + 146 * ncomps + c);
        const auto *md_147 = buffer.data(md + 147 * ncomps + c);
        const auto *md_148 = buffer.data(md + 148 * ncomps + c);
        const auto *md_149 = buffer.data(md + 149 * ncomps + c);
        const auto *md_150 = buffer.data(md + 150 * ncomps + c);
        const auto *md_151 = buffer.data(md + 151 * ncomps + c);
        const auto *md_152 = buffer.data(md + 152 * ncomps + c);
        const auto *md_153 = buffer.data(md + 153 * ncomps + c);
        const auto *md_154 = buffer.data(md + 154 * ncomps + c);
        const auto *md_155 = buffer.data(md + 155 * ncomps + c);
        const auto *md_156 = buffer.data(md + 156 * ncomps + c);
        const auto *md_157 = buffer.data(md + 157 * ncomps + c);
        const auto *md_158 = buffer.data(md + 158 * ncomps + c);
        const auto *md_159 = buffer.data(md + 159 * ncomps + c);
        const auto *md_160 = buffer.data(md + 160 * ncomps + c);
        const auto *md_161 = buffer.data(md + 161 * ncomps + c);
        const auto *md_162 = buffer.data(md + 162 * ncomps + c);
        const auto *md_163 = buffer.data(md + 163 * ncomps + c);
        const auto *md_164 = buffer.data(md + 164 * ncomps + c);
        const auto *md_165 = buffer.data(md + 165 * ncomps + c);
        const auto *md_166 = buffer.data(md + 166 * ncomps + c);
        const auto *md_167 = buffer.data(md + 167 * ncomps + c);
        const auto *md_168 = buffer.data(md + 168 * ncomps + c);
        const auto *md_169 = buffer.data(md + 169 * ncomps + c);
        const auto *md_170 = buffer.data(md + 170 * ncomps + c);
        const auto *md_171 = buffer.data(md + 171 * ncomps + c);
        const auto *md_172 = buffer.data(md + 172 * ncomps + c);
        const auto *md_173 = buffer.data(md + 173 * ncomps + c);

        const auto *nd_89 = buffer.data(nd + 89 * ncomps + c);
        const auto *nd_90 = buffer.data(nd + 90 * ncomps + c);
        const auto *nd_91 = buffer.data(nd + 91 * ncomps + c);
        const auto *nd_92 = buffer.data(nd + 92 * ncomps + c);
        const auto *nd_93 = buffer.data(nd + 93 * ncomps + c);
        const auto *nd_94 = buffer.data(nd + 94 * ncomps + c);
        const auto *nd_95 = buffer.data(nd + 95 * ncomps + c);
        const auto *nd_96 = buffer.data(nd + 96 * ncomps + c);
        const auto *nd_97 = buffer.data(nd + 97 * ncomps + c);
        const auto *nd_98 = buffer.data(nd + 98 * ncomps + c);
        const auto *nd_99 = buffer.data(nd + 99 * ncomps + c);
        const auto *nd_100 = buffer.data(nd + 100 * ncomps + c);
        const auto *nd_101 = buffer.data(nd + 101 * ncomps + c);
        const auto *nd_102 = buffer.data(nd + 102 * ncomps + c);
        const auto *nd_103 = buffer.data(nd + 103 * ncomps + c);
        const auto *nd_104 = buffer.data(nd + 104 * ncomps + c);
        const auto *nd_105 = buffer.data(nd + 105 * ncomps + c);
        const auto *nd_106 = buffer.data(nd + 106 * ncomps + c);
        const auto *nd_107 = buffer.data(nd + 107 * ncomps + c);
        const auto *nd_108 = buffer.data(nd + 108 * ncomps + c);
        const auto *nd_109 = buffer.data(nd + 109 * ncomps + c);
        const auto *nd_110 = buffer.data(nd + 110 * ncomps + c);
        const auto *nd_111 = buffer.data(nd + 111 * ncomps + c);
        const auto *nd_112 = buffer.data(nd + 112 * ncomps + c);
        const auto *nd_113 = buffer.data(nd + 113 * ncomps + c);
        const auto *nd_114 = buffer.data(nd + 114 * ncomps + c);
        const auto *nd_115 = buffer.data(nd + 115 * ncomps + c);
        const auto *nd_116 = buffer.data(nd + 116 * ncomps + c);
        const auto *nd_117 = buffer.data(nd + 117 * ncomps + c);
        const auto *nd_118 = buffer.data(nd + 118 * ncomps + c);
        const auto *nd_119 = buffer.data(nd + 119 * ncomps + c);
        const auto *nd_120 = buffer.data(nd + 120 * ncomps + c);
        const auto *nd_121 = buffer.data(nd + 121 * ncomps + c);
        const auto *nd_122 = buffer.data(nd + 122 * ncomps + c);
        const auto *nd_123 = buffer.data(nd + 123 * ncomps + c);
        const auto *nd_124 = buffer.data(nd + 124 * ncomps + c);
        const auto *nd_125 = buffer.data(nd + 125 * ncomps + c);
        const auto *nd_126 = buffer.data(nd + 126 * ncomps + c);
        const auto *nd_127 = buffer.data(nd + 127 * ncomps + c);
        const auto *nd_128 = buffer.data(nd + 128 * ncomps + c);
        const auto *nd_129 = buffer.data(nd + 129 * ncomps + c);
        const auto *nd_130 = buffer.data(nd + 130 * ncomps + c);
        const auto *nd_131 = buffer.data(nd + 131 * ncomps + c);
        const auto *nd_132 = buffer.data(nd + 132 * ncomps + c);
        const auto *nd_133 = buffer.data(nd + 133 * ncomps + c);
        const auto *nd_134 = buffer.data(nd + 134 * ncomps + c);
        const auto *nd_135 = buffer.data(nd + 135 * ncomps + c);
        const auto *nd_136 = buffer.data(nd + 136 * ncomps + c);
        const auto *nd_137 = buffer.data(nd + 137 * ncomps + c);
        const auto *nd_138 = buffer.data(nd + 138 * ncomps + c);
        const auto *nd_139 = buffer.data(nd + 139 * ncomps + c);
        const auto *nd_140 = buffer.data(nd + 140 * ncomps + c);
        const auto *nd_141 = buffer.data(nd + 141 * ncomps + c);
        const auto *nd_142 = buffer.data(nd + 142 * ncomps + c);
        const auto *nd_143 = buffer.data(nd + 143 * ncomps + c);
        const auto *nd_144 = buffer.data(nd + 144 * ncomps + c);
        const auto *nd_145 = buffer.data(nd + 145 * ncomps + c);
        const auto *nd_146 = buffer.data(nd + 146 * ncomps + c);
        const auto *nd_147 = buffer.data(nd + 147 * ncomps + c);
        const auto *nd_148 = buffer.data(nd + 148 * ncomps + c);
        const auto *nd_149 = buffer.data(nd + 149 * ncomps + c);
        const auto *nd_150 = buffer.data(nd + 150 * ncomps + c);
        const auto *nd_151 = buffer.data(nd + 151 * ncomps + c);
        const auto *nd_152 = buffer.data(nd + 152 * ncomps + c);
        const auto *nd_153 = buffer.data(nd + 153 * ncomps + c);
        const auto *nd_154 = buffer.data(nd + 154 * ncomps + c);
        const auto *nd_155 = buffer.data(nd + 155 * ncomps + c);
        const auto *nd_156 = buffer.data(nd + 156 * ncomps + c);
        const auto *nd_157 = buffer.data(nd + 157 * ncomps + c);
        const auto *nd_158 = buffer.data(nd + 158 * ncomps + c);
        const auto *nd_159 = buffer.data(nd + 159 * ncomps + c);
        const auto *nd_160 = buffer.data(nd + 160 * ncomps + c);
        const auto *nd_161 = buffer.data(nd + 161 * ncomps + c);
        const auto *nd_162 = buffer.data(nd + 162 * ncomps + c);
        const auto *nd_163 = buffer.data(nd + 163 * ncomps + c);
        const auto *nd_164 = buffer.data(nd + 164 * ncomps + c);
        const auto *nd_165 = buffer.data(nd + 165 * ncomps + c);
        const auto *nd_166 = buffer.data(nd + 166 * ncomps + c);
        const auto *nd_167 = buffer.data(nd + 167 * ncomps + c);
        const auto *nd_168 = buffer.data(nd + 168 * ncomps + c);
        const auto *nd_169 = buffer.data(nd + 169 * ncomps + c);
        const auto *nd_170 = buffer.data(nd + 170 * ncomps + c);
        const auto *nd_171 = buffer.data(nd + 171 * ncomps + c);
        const auto *nd_172 = buffer.data(nd + 172 * ncomps + c);
        const auto *nd_173 = buffer.data(nd + 173 * ncomps + c);
        const auto *nd_177 = buffer.data(nd + 177 * ncomps + c);
        const auto *nd_178 = buffer.data(nd + 178 * ncomps + c);
        const auto *nd_179 = buffer.data(nd + 179 * ncomps + c);
        const auto *nd_183 = buffer.data(nd + 183 * ncomps + c);
        const auto *nd_184 = buffer.data(nd + 184 * ncomps + c);
        const auto *nd_185 = buffer.data(nd + 185 * ncomps + c);
        const auto *nd_189 = buffer.data(nd + 189 * ncomps + c);
        const auto *nd_190 = buffer.data(nd + 190 * ncomps + c);
        const auto *nd_191 = buffer.data(nd + 191 * ncomps + c);
        const auto *nd_195 = buffer.data(nd + 195 * ncomps + c);
        const auto *nd_196 = buffer.data(nd + 196 * ncomps + c);
        const auto *nd_197 = buffer.data(nd + 197 * ncomps + c);
        const auto *nd_201 = buffer.data(nd + 201 * ncomps + c);
        const auto *nd_202 = buffer.data(nd + 202 * ncomps + c);
        const auto *nd_203 = buffer.data(nd + 203 * ncomps + c);
        const auto *nd_207 = buffer.data(nd + 207 * ncomps + c);
        const auto *nd_208 = buffer.data(nd + 208 * ncomps + c);
        const auto *nd_209 = buffer.data(nd + 209 * ncomps + c);
        const auto *nd_215 = buffer.data(nd + 215 * ncomps + c);
        const auto *nd_219 = buffer.data(nd + 219 * ncomps + c);
        const auto *nd_220 = buffer.data(nd + 220 * ncomps + c);
        const auto *nd_221 = buffer.data(nd + 221 * ncomps + c);
        const auto *nd_227 = buffer.data(nd + 227 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, ab_y, ab_z, md_87, md_88, \
                         md_89, nd_89, nd_117, nd_118, nd_119, nd_125 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = ab_x[k] * md_89[k]
                       + nd_89[k];

            t_146[k] = ab_y[k] * md_87[k]
                       + nd_117[k];

            t_147[k] = ab_y[k] * md_88[k]
                       + nd_118[k];

            t_148[k] = ab_y[k] * md_89[k]
                       + nd_119[k];

            t_149[k] = ab_z[k] * md_89[k]
                       + nd_125[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, md_90, md_91, md_92, md_93, \
                         md_94, nd_90, nd_91, nd_92, nd_93, nd_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * md_90[k]
                       + nd_90[k];

            t_151[k] = ab_x[k] * md_91[k]
                       + nd_91[k];

            t_152[k] = ab_x[k] * md_92[k]
                       + nd_92[k];

            t_153[k] = ab_x[k] * md_93[k]
                       + nd_93[k];

            t_154[k] = ab_x[k] * md_94[k]
                       + nd_94[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, ab_y, ab_z, md_93, md_94, \
                         md_95, nd_95, nd_129, nd_130, nd_131, nd_137 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * md_95[k]
                       + nd_95[k];

            t_156[k] = ab_y[k] * md_93[k]
                       + nd_129[k];

            t_157[k] = ab_y[k] * md_94[k]
                       + nd_130[k];

            t_158[k] = ab_y[k] * md_95[k]
                       + nd_131[k];

            t_159[k] = ab_z[k] * md_95[k]
                       + nd_137[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, md_96, md_97, md_98, md_99, \
                         md_100, nd_96, nd_97, nd_98, nd_99, nd_100 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_x[k] * md_96[k]
                       + nd_96[k];

            t_161[k] = ab_x[k] * md_97[k]
                       + nd_97[k];

            t_162[k] = ab_x[k] * md_98[k]
                       + nd_98[k];

            t_163[k] = ab_x[k] * md_99[k]
                       + nd_99[k];

            t_164[k] = ab_x[k] * md_100[k]
                       + nd_100[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, ab_y, ab_z, md_99, md_100, \
                         md_101, nd_101, nd_135, nd_136, nd_137, \
                         nd_143 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = ab_x[k] * md_101[k]
                       + nd_101[k];

            t_166[k] = ab_y[k] * md_99[k]
                       + nd_135[k];

            t_167[k] = ab_y[k] * md_100[k]
                       + nd_136[k];

            t_168[k] = ab_y[k] * md_101[k]
                       + nd_137[k];

            t_169[k] = ab_z[k] * md_101[k]
                       + nd_143[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, md_102, md_103, md_104, \
                         md_105, md_106, nd_102, nd_103, nd_104, nd_105, \
                         nd_106 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = ab_x[k] * md_102[k]
                       + nd_102[k];

            t_171[k] = ab_x[k] * md_103[k]
                       + nd_103[k];

            t_172[k] = ab_x[k] * md_104[k]
                       + nd_104[k];

            t_173[k] = ab_x[k] * md_105[k]
                       + nd_105[k];

            t_174[k] = ab_x[k] * md_106[k]
                       + nd_106[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, ab_y, ab_z, md_105, md_106, \
                         md_107, nd_107, nd_141, nd_142, nd_143, \
                         nd_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = ab_x[k] * md_107[k]
                       + nd_107[k];

            t_176[k] = ab_y[k] * md_105[k]
                       + nd_141[k];

            t_177[k] = ab_y[k] * md_106[k]
                       + nd_142[k];

            t_178[k] = ab_y[k] * md_107[k]
                       + nd_143[k];

            t_179[k] = ab_z[k] * md_107[k]
                       + nd_149[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, md_108, md_109, md_110, \
                         md_111, md_112, nd_108, nd_109, nd_110, nd_111, \
                         nd_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = ab_x[k] * md_108[k]
                       + nd_108[k];

            t_181[k] = ab_x[k] * md_109[k]
                       + nd_109[k];

            t_182[k] = ab_x[k] * md_110[k]
                       + nd_110[k];

            t_183[k] = ab_x[k] * md_111[k]
                       + nd_111[k];

            t_184[k] = ab_x[k] * md_112[k]
                       + nd_112[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, ab_y, ab_z, md_111, md_112, \
                         md_113, nd_113, nd_147, nd_148, nd_149, \
                         nd_155 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_x[k] * md_113[k]
                       + nd_113[k];

            t_186[k] = ab_y[k] * md_111[k]
                       + nd_147[k];

            t_187[k] = ab_y[k] * md_112[k]
                       + nd_148[k];

            t_188[k] = ab_y[k] * md_113[k]
                       + nd_149[k];

            t_189[k] = ab_z[k] * md_113[k]
                       + nd_155[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, md_114, md_115, md_116, \
                         md_117, md_118, nd_114, nd_115, nd_116, nd_117, \
                         nd_118 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = ab_x[k] * md_114[k]
                       + nd_114[k];

            t_191[k] = ab_x[k] * md_115[k]
                       + nd_115[k];

            t_192[k] = ab_x[k] * md_116[k]
                       + nd_116[k];

            t_193[k] = ab_x[k] * md_117[k]
                       + nd_117[k];

            t_194[k] = ab_x[k] * md_118[k]
                       + nd_118[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, ab_y, ab_z, md_117, md_118, \
                         md_119, nd_119, nd_153, nd_154, nd_155, \
                         nd_161 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = ab_x[k] * md_119[k]
                       + nd_119[k];

            t_196[k] = ab_y[k] * md_117[k]
                       + nd_153[k];

            t_197[k] = ab_y[k] * md_118[k]
                       + nd_154[k];

            t_198[k] = ab_y[k] * md_119[k]
                       + nd_155[k];

            t_199[k] = ab_z[k] * md_119[k]
                       + nd_161[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, md_120, md_121, md_122, \
                         md_123, md_124, nd_120, nd_121, nd_122, nd_123, \
                         nd_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = ab_x[k] * md_120[k]
                       + nd_120[k];

            t_201[k] = ab_x[k] * md_121[k]
                       + nd_121[k];

            t_202[k] = ab_x[k] * md_122[k]
                       + nd_122[k];

            t_203[k] = ab_x[k] * md_123[k]
                       + nd_123[k];

            t_204[k] = ab_x[k] * md_124[k]
                       + nd_124[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, ab_y, ab_z, md_123, md_124, \
                         md_125, nd_125, nd_159, nd_160, nd_161, \
                         nd_167 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = ab_x[k] * md_125[k]
                       + nd_125[k];

            t_206[k] = ab_y[k] * md_123[k]
                       + nd_159[k];

            t_207[k] = ab_y[k] * md_124[k]
                       + nd_160[k];

            t_208[k] = ab_y[k] * md_125[k]
                       + nd_161[k];

            t_209[k] = ab_z[k] * md_125[k]
                       + nd_167[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, md_126, md_127, md_128, \
                         md_129, md_130, nd_126, nd_127, nd_128, nd_129, \
                         nd_130 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = ab_x[k] * md_126[k]
                       + nd_126[k];

            t_211[k] = ab_x[k] * md_127[k]
                       + nd_127[k];

            t_212[k] = ab_x[k] * md_128[k]
                       + nd_128[k];

            t_213[k] = ab_x[k] * md_129[k]
                       + nd_129[k];

            t_214[k] = ab_x[k] * md_130[k]
                       + nd_130[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, ab_y, ab_z, md_129, md_130, \
                         md_131, nd_131, nd_171, nd_172, nd_173, \
                         nd_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = ab_x[k] * md_131[k]
                       + nd_131[k];

            t_216[k] = ab_y[k] * md_129[k]
                       + nd_171[k];

            t_217[k] = ab_y[k] * md_130[k]
                       + nd_172[k];

            t_218[k] = ab_y[k] * md_131[k]
                       + nd_173[k];

            t_219[k] = ab_z[k] * md_131[k]
                       + nd_179[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, md_132, md_133, md_134, \
                         md_135, md_136, nd_132, nd_133, nd_134, nd_135, \
                         nd_136 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = ab_x[k] * md_132[k]
                       + nd_132[k];

            t_221[k] = ab_x[k] * md_133[k]
                       + nd_133[k];

            t_222[k] = ab_x[k] * md_134[k]
                       + nd_134[k];

            t_223[k] = ab_x[k] * md_135[k]
                       + nd_135[k];

            t_224[k] = ab_x[k] * md_136[k]
                       + nd_136[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, ab_y, ab_z, md_135, md_136, \
                         md_137, nd_137, nd_177, nd_178, nd_179, \
                         nd_185 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = ab_x[k] * md_137[k]
                       + nd_137[k];

            t_226[k] = ab_y[k] * md_135[k]
                       + nd_177[k];

            t_227[k] = ab_y[k] * md_136[k]
                       + nd_178[k];

            t_228[k] = ab_y[k] * md_137[k]
                       + nd_179[k];

            t_229[k] = ab_z[k] * md_137[k]
                       + nd_185[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, md_138, md_139, md_140, \
                         md_141, md_142, nd_138, nd_139, nd_140, nd_141, \
                         nd_142 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = ab_x[k] * md_138[k]
                       + nd_138[k];

            t_231[k] = ab_x[k] * md_139[k]
                       + nd_139[k];

            t_232[k] = ab_x[k] * md_140[k]
                       + nd_140[k];

            t_233[k] = ab_x[k] * md_141[k]
                       + nd_141[k];

            t_234[k] = ab_x[k] * md_142[k]
                       + nd_142[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_x, ab_y, ab_z, md_141, md_142, \
                         md_143, nd_143, nd_183, nd_184, nd_185, \
                         nd_191 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = ab_x[k] * md_143[k]
                       + nd_143[k];

            t_236[k] = ab_y[k] * md_141[k]
                       + nd_183[k];

            t_237[k] = ab_y[k] * md_142[k]
                       + nd_184[k];

            t_238[k] = ab_y[k] * md_143[k]
                       + nd_185[k];

            t_239[k] = ab_z[k] * md_143[k]
                       + nd_191[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, md_144, md_145, md_146, \
                         md_147, md_148, nd_144, nd_145, nd_146, nd_147, \
                         nd_148 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = ab_x[k] * md_144[k]
                       + nd_144[k];

            t_241[k] = ab_x[k] * md_145[k]
                       + nd_145[k];

            t_242[k] = ab_x[k] * md_146[k]
                       + nd_146[k];

            t_243[k] = ab_x[k] * md_147[k]
                       + nd_147[k];

            t_244[k] = ab_x[k] * md_148[k]
                       + nd_148[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, ab_y, ab_z, md_147, md_148, \
                         md_149, nd_149, nd_189, nd_190, nd_191, \
                         nd_197 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = ab_x[k] * md_149[k]
                       + nd_149[k];

            t_246[k] = ab_y[k] * md_147[k]
                       + nd_189[k];

            t_247[k] = ab_y[k] * md_148[k]
                       + nd_190[k];

            t_248[k] = ab_y[k] * md_149[k]
                       + nd_191[k];

            t_249[k] = ab_z[k] * md_149[k]
                       + nd_197[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_x, md_150, md_151, md_152, \
                         md_153, md_154, nd_150, nd_151, nd_152, nd_153, \
                         nd_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = ab_x[k] * md_150[k]
                       + nd_150[k];

            t_251[k] = ab_x[k] * md_151[k]
                       + nd_151[k];

            t_252[k] = ab_x[k] * md_152[k]
                       + nd_152[k];

            t_253[k] = ab_x[k] * md_153[k]
                       + nd_153[k];

            t_254[k] = ab_x[k] * md_154[k]
                       + nd_154[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, ab_y, ab_z, md_153, md_154, \
                         md_155, nd_155, nd_195, nd_196, nd_197, \
                         nd_203 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = ab_x[k] * md_155[k]
                       + nd_155[k];

            t_256[k] = ab_y[k] * md_153[k]
                       + nd_195[k];

            t_257[k] = ab_y[k] * md_154[k]
                       + nd_196[k];

            t_258[k] = ab_y[k] * md_155[k]
                       + nd_197[k];

            t_259[k] = ab_z[k] * md_155[k]
                       + nd_203[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, md_156, md_157, md_158, \
                         md_159, md_160, nd_156, nd_157, nd_158, nd_159, \
                         nd_160 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = ab_x[k] * md_156[k]
                       + nd_156[k];

            t_261[k] = ab_x[k] * md_157[k]
                       + nd_157[k];

            t_262[k] = ab_x[k] * md_158[k]
                       + nd_158[k];

            t_263[k] = ab_x[k] * md_159[k]
                       + nd_159[k];

            t_264[k] = ab_x[k] * md_160[k]
                       + nd_160[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, ab_y, ab_z, md_159, md_160, \
                         md_161, nd_161, nd_201, nd_202, nd_203, \
                         nd_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = ab_x[k] * md_161[k]
                       + nd_161[k];

            t_266[k] = ab_y[k] * md_159[k]
                       + nd_201[k];

            t_267[k] = ab_y[k] * md_160[k]
                       + nd_202[k];

            t_268[k] = ab_y[k] * md_161[k]
                       + nd_203[k];

            t_269[k] = ab_z[k] * md_161[k]
                       + nd_209[k];
        }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, md_162, md_163, md_164, \
                         md_165, md_166, nd_162, nd_163, nd_164, nd_165, \
                         nd_166 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_270[k] = ab_x[k] * md_162[k]
                       + nd_162[k];

            t_271[k] = ab_x[k] * md_163[k]
                       + nd_163[k];

            t_272[k] = ab_x[k] * md_164[k]
                       + nd_164[k];

            t_273[k] = ab_x[k] * md_165[k]
                       + nd_165[k];

            t_274[k] = ab_x[k] * md_166[k]
                       + nd_166[k];
        }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, ab_y, ab_z, md_165, md_166, \
                         md_167, nd_167, nd_207, nd_208, nd_209, \
                         nd_215 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_275[k] = ab_x[k] * md_167[k]
                       + nd_167[k];

            t_276[k] = ab_y[k] * md_165[k]
                       + nd_207[k];

            t_277[k] = ab_y[k] * md_166[k]
                       + nd_208[k];

            t_278[k] = ab_y[k] * md_167[k]
                       + nd_209[k];

            t_279[k] = ab_z[k] * md_167[k]
                       + nd_215[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_x, md_168, md_169, md_170, \
                         md_171, md_172, nd_168, nd_169, nd_170, nd_171, \
                         nd_172 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = ab_x[k] * md_168[k]
                       + nd_168[k];

            t_281[k] = ab_x[k] * md_169[k]
                       + nd_169[k];

            t_282[k] = ab_x[k] * md_170[k]
                       + nd_170[k];

            t_283[k] = ab_x[k] * md_171[k]
                       + nd_171[k];

            t_284[k] = ab_x[k] * md_172[k]
                       + nd_172[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, ab_y, ab_z, md_171, md_172, \
                         md_173, nd_173, nd_219, nd_220, nd_221, \
                         nd_227 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = ab_x[k] * md_173[k]
                       + nd_173[k];

            t_286[k] = ab_y[k] * md_171[k]
                       + nd_219[k];

            t_287[k] = ab_y[k] * md_172[k]
                       + nd_220[k];

            t_288[k] = ab_y[k] * md_173[k]
                       + nd_221[k];

            t_289[k] = ab_z[k] * md_173[k]
                       + nd_227[k];
        }
    }
}

static auto
compute_hrr_mf_out_of_first_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t md, const size_t nd,
                                   const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_290 = buffer.data(target + 290 * ncomps + c);
        auto *t_291 = buffer.data(target + 291 * ncomps + c);
        auto *t_292 = buffer.data(target + 292 * ncomps + c);
        auto *t_293 = buffer.data(target + 293 * ncomps + c);
        auto *t_294 = buffer.data(target + 294 * ncomps + c);
        auto *t_295 = buffer.data(target + 295 * ncomps + c);
        auto *t_296 = buffer.data(target + 296 * ncomps + c);
        auto *t_297 = buffer.data(target + 297 * ncomps + c);
        auto *t_298 = buffer.data(target + 298 * ncomps + c);
        auto *t_299 = buffer.data(target + 299 * ncomps + c);
        auto *t_300 = buffer.data(target + 300 * ncomps + c);
        auto *t_301 = buffer.data(target + 301 * ncomps + c);
        auto *t_302 = buffer.data(target + 302 * ncomps + c);
        auto *t_303 = buffer.data(target + 303 * ncomps + c);
        auto *t_304 = buffer.data(target + 304 * ncomps + c);
        auto *t_305 = buffer.data(target + 305 * ncomps + c);
        auto *t_306 = buffer.data(target + 306 * ncomps + c);
        auto *t_307 = buffer.data(target + 307 * ncomps + c);
        auto *t_308 = buffer.data(target + 308 * ncomps + c);
        auto *t_309 = buffer.data(target + 309 * ncomps + c);
        auto *t_310 = buffer.data(target + 310 * ncomps + c);
        auto *t_311 = buffer.data(target + 311 * ncomps + c);
        auto *t_312 = buffer.data(target + 312 * ncomps + c);
        auto *t_313 = buffer.data(target + 313 * ncomps + c);
        auto *t_314 = buffer.data(target + 314 * ncomps + c);
        auto *t_315 = buffer.data(target + 315 * ncomps + c);
        auto *t_316 = buffer.data(target + 316 * ncomps + c);
        auto *t_317 = buffer.data(target + 317 * ncomps + c);
        auto *t_318 = buffer.data(target + 318 * ncomps + c);
        auto *t_319 = buffer.data(target + 319 * ncomps + c);
        auto *t_320 = buffer.data(target + 320 * ncomps + c);
        auto *t_321 = buffer.data(target + 321 * ncomps + c);
        auto *t_322 = buffer.data(target + 322 * ncomps + c);
        auto *t_323 = buffer.data(target + 323 * ncomps + c);
        auto *t_324 = buffer.data(target + 324 * ncomps + c);
        auto *t_325 = buffer.data(target + 325 * ncomps + c);
        auto *t_326 = buffer.data(target + 326 * ncomps + c);
        auto *t_327 = buffer.data(target + 327 * ncomps + c);
        auto *t_328 = buffer.data(target + 328 * ncomps + c);
        auto *t_329 = buffer.data(target + 329 * ncomps + c);
        auto *t_330 = buffer.data(target + 330 * ncomps + c);
        auto *t_331 = buffer.data(target + 331 * ncomps + c);
        auto *t_332 = buffer.data(target + 332 * ncomps + c);
        auto *t_333 = buffer.data(target + 333 * ncomps + c);
        auto *t_334 = buffer.data(target + 334 * ncomps + c);
        auto *t_335 = buffer.data(target + 335 * ncomps + c);
        auto *t_336 = buffer.data(target + 336 * ncomps + c);
        auto *t_337 = buffer.data(target + 337 * ncomps + c);
        auto *t_338 = buffer.data(target + 338 * ncomps + c);
        auto *t_339 = buffer.data(target + 339 * ncomps + c);
        auto *t_340 = buffer.data(target + 340 * ncomps + c);
        auto *t_341 = buffer.data(target + 341 * ncomps + c);
        auto *t_342 = buffer.data(target + 342 * ncomps + c);
        auto *t_343 = buffer.data(target + 343 * ncomps + c);
        auto *t_344 = buffer.data(target + 344 * ncomps + c);
        auto *t_345 = buffer.data(target + 345 * ncomps + c);
        auto *t_346 = buffer.data(target + 346 * ncomps + c);
        auto *t_347 = buffer.data(target + 347 * ncomps + c);
        auto *t_348 = buffer.data(target + 348 * ncomps + c);
        auto *t_349 = buffer.data(target + 349 * ncomps + c);
        auto *t_350 = buffer.data(target + 350 * ncomps + c);
        auto *t_351 = buffer.data(target + 351 * ncomps + c);
        auto *t_352 = buffer.data(target + 352 * ncomps + c);
        auto *t_353 = buffer.data(target + 353 * ncomps + c);
        auto *t_354 = buffer.data(target + 354 * ncomps + c);
        auto *t_355 = buffer.data(target + 355 * ncomps + c);
        auto *t_356 = buffer.data(target + 356 * ncomps + c);
        auto *t_357 = buffer.data(target + 357 * ncomps + c);
        auto *t_358 = buffer.data(target + 358 * ncomps + c);
        auto *t_359 = buffer.data(target + 359 * ncomps + c);
        auto *t_360 = buffer.data(target + 360 * ncomps + c);
        auto *t_361 = buffer.data(target + 361 * ncomps + c);
        auto *t_362 = buffer.data(target + 362 * ncomps + c);
        auto *t_363 = buffer.data(target + 363 * ncomps + c);
        auto *t_364 = buffer.data(target + 364 * ncomps + c);
        auto *t_365 = buffer.data(target + 365 * ncomps + c);
        auto *t_366 = buffer.data(target + 366 * ncomps + c);
        auto *t_367 = buffer.data(target + 367 * ncomps + c);
        auto *t_368 = buffer.data(target + 368 * ncomps + c);
        auto *t_369 = buffer.data(target + 369 * ncomps + c);
        auto *t_370 = buffer.data(target + 370 * ncomps + c);
        auto *t_371 = buffer.data(target + 371 * ncomps + c);
        auto *t_372 = buffer.data(target + 372 * ncomps + c);
        auto *t_373 = buffer.data(target + 373 * ncomps + c);
        auto *t_374 = buffer.data(target + 374 * ncomps + c);
        auto *t_375 = buffer.data(target + 375 * ncomps + c);
        auto *t_376 = buffer.data(target + 376 * ncomps + c);
        auto *t_377 = buffer.data(target + 377 * ncomps + c);
        auto *t_378 = buffer.data(target + 378 * ncomps + c);
        auto *t_379 = buffer.data(target + 379 * ncomps + c);
        auto *t_380 = buffer.data(target + 380 * ncomps + c);
        auto *t_381 = buffer.data(target + 381 * ncomps + c);
        auto *t_382 = buffer.data(target + 382 * ncomps + c);
        auto *t_383 = buffer.data(target + 383 * ncomps + c);
        auto *t_384 = buffer.data(target + 384 * ncomps + c);
        auto *t_385 = buffer.data(target + 385 * ncomps + c);
        auto *t_386 = buffer.data(target + 386 * ncomps + c);
        auto *t_387 = buffer.data(target + 387 * ncomps + c);
        auto *t_388 = buffer.data(target + 388 * ncomps + c);
        auto *t_389 = buffer.data(target + 389 * ncomps + c);
        auto *t_390 = buffer.data(target + 390 * ncomps + c);
        auto *t_391 = buffer.data(target + 391 * ncomps + c);
        auto *t_392 = buffer.data(target + 392 * ncomps + c);
        auto *t_393 = buffer.data(target + 393 * ncomps + c);
        auto *t_394 = buffer.data(target + 394 * ncomps + c);
        auto *t_395 = buffer.data(target + 395 * ncomps + c);
        auto *t_396 = buffer.data(target + 396 * ncomps + c);
        auto *t_397 = buffer.data(target + 397 * ncomps + c);
        auto *t_398 = buffer.data(target + 398 * ncomps + c);
        auto *t_399 = buffer.data(target + 399 * ncomps + c);
        auto *t_400 = buffer.data(target + 400 * ncomps + c);
        auto *t_401 = buffer.data(target + 401 * ncomps + c);
        auto *t_402 = buffer.data(target + 402 * ncomps + c);
        auto *t_403 = buffer.data(target + 403 * ncomps + c);
        auto *t_404 = buffer.data(target + 404 * ncomps + c);
        auto *t_405 = buffer.data(target + 405 * ncomps + c);
        auto *t_406 = buffer.data(target + 406 * ncomps + c);
        auto *t_407 = buffer.data(target + 407 * ncomps + c);
        auto *t_408 = buffer.data(target + 408 * ncomps + c);
        auto *t_409 = buffer.data(target + 409 * ncomps + c);
        auto *t_410 = buffer.data(target + 410 * ncomps + c);
        auto *t_411 = buffer.data(target + 411 * ncomps + c);
        auto *t_412 = buffer.data(target + 412 * ncomps + c);
        auto *t_413 = buffer.data(target + 413 * ncomps + c);
        auto *t_414 = buffer.data(target + 414 * ncomps + c);
        auto *t_415 = buffer.data(target + 415 * ncomps + c);
        auto *t_416 = buffer.data(target + 416 * ncomps + c);
        auto *t_417 = buffer.data(target + 417 * ncomps + c);
        auto *t_418 = buffer.data(target + 418 * ncomps + c);
        auto *t_419 = buffer.data(target + 419 * ncomps + c);
        auto *t_420 = buffer.data(target + 420 * ncomps + c);
        auto *t_421 = buffer.data(target + 421 * ncomps + c);
        auto *t_422 = buffer.data(target + 422 * ncomps + c);
        auto *t_423 = buffer.data(target + 423 * ncomps + c);
        auto *t_424 = buffer.data(target + 424 * ncomps + c);
        auto *t_425 = buffer.data(target + 425 * ncomps + c);
        auto *t_426 = buffer.data(target + 426 * ncomps + c);
        auto *t_427 = buffer.data(target + 427 * ncomps + c);
        auto *t_428 = buffer.data(target + 428 * ncomps + c);
        auto *t_429 = buffer.data(target + 429 * ncomps + c);
        auto *t_430 = buffer.data(target + 430 * ncomps + c);
        auto *t_431 = buffer.data(target + 431 * ncomps + c);
        auto *t_432 = buffer.data(target + 432 * ncomps + c);
        auto *t_433 = buffer.data(target + 433 * ncomps + c);
        auto *t_434 = buffer.data(target + 434 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *md_174 = buffer.data(md + 174 * ncomps + c);
        const auto *md_175 = buffer.data(md + 175 * ncomps + c);
        const auto *md_176 = buffer.data(md + 176 * ncomps + c);
        const auto *md_177 = buffer.data(md + 177 * ncomps + c);
        const auto *md_178 = buffer.data(md + 178 * ncomps + c);
        const auto *md_179 = buffer.data(md + 179 * ncomps + c);
        const auto *md_180 = buffer.data(md + 180 * ncomps + c);
        const auto *md_181 = buffer.data(md + 181 * ncomps + c);
        const auto *md_182 = buffer.data(md + 182 * ncomps + c);
        const auto *md_183 = buffer.data(md + 183 * ncomps + c);
        const auto *md_184 = buffer.data(md + 184 * ncomps + c);
        const auto *md_185 = buffer.data(md + 185 * ncomps + c);
        const auto *md_186 = buffer.data(md + 186 * ncomps + c);
        const auto *md_187 = buffer.data(md + 187 * ncomps + c);
        const auto *md_188 = buffer.data(md + 188 * ncomps + c);
        const auto *md_189 = buffer.data(md + 189 * ncomps + c);
        const auto *md_190 = buffer.data(md + 190 * ncomps + c);
        const auto *md_191 = buffer.data(md + 191 * ncomps + c);
        const auto *md_192 = buffer.data(md + 192 * ncomps + c);
        const auto *md_193 = buffer.data(md + 193 * ncomps + c);
        const auto *md_194 = buffer.data(md + 194 * ncomps + c);
        const auto *md_195 = buffer.data(md + 195 * ncomps + c);
        const auto *md_196 = buffer.data(md + 196 * ncomps + c);
        const auto *md_197 = buffer.data(md + 197 * ncomps + c);
        const auto *md_198 = buffer.data(md + 198 * ncomps + c);
        const auto *md_199 = buffer.data(md + 199 * ncomps + c);
        const auto *md_200 = buffer.data(md + 200 * ncomps + c);
        const auto *md_201 = buffer.data(md + 201 * ncomps + c);
        const auto *md_202 = buffer.data(md + 202 * ncomps + c);
        const auto *md_203 = buffer.data(md + 203 * ncomps + c);
        const auto *md_204 = buffer.data(md + 204 * ncomps + c);
        const auto *md_205 = buffer.data(md + 205 * ncomps + c);
        const auto *md_206 = buffer.data(md + 206 * ncomps + c);
        const auto *md_207 = buffer.data(md + 207 * ncomps + c);
        const auto *md_208 = buffer.data(md + 208 * ncomps + c);
        const auto *md_209 = buffer.data(md + 209 * ncomps + c);
        const auto *md_210 = buffer.data(md + 210 * ncomps + c);
        const auto *md_211 = buffer.data(md + 211 * ncomps + c);
        const auto *md_212 = buffer.data(md + 212 * ncomps + c);
        const auto *md_213 = buffer.data(md + 213 * ncomps + c);
        const auto *md_214 = buffer.data(md + 214 * ncomps + c);
        const auto *md_215 = buffer.data(md + 215 * ncomps + c);
        const auto *md_216 = buffer.data(md + 216 * ncomps + c);
        const auto *md_217 = buffer.data(md + 217 * ncomps + c);
        const auto *md_218 = buffer.data(md + 218 * ncomps + c);
        const auto *md_219 = buffer.data(md + 219 * ncomps + c);
        const auto *md_220 = buffer.data(md + 220 * ncomps + c);
        const auto *md_221 = buffer.data(md + 221 * ncomps + c);
        const auto *md_222 = buffer.data(md + 222 * ncomps + c);
        const auto *md_223 = buffer.data(md + 223 * ncomps + c);
        const auto *md_224 = buffer.data(md + 224 * ncomps + c);
        const auto *md_225 = buffer.data(md + 225 * ncomps + c);
        const auto *md_226 = buffer.data(md + 226 * ncomps + c);
        const auto *md_227 = buffer.data(md + 227 * ncomps + c);
        const auto *md_228 = buffer.data(md + 228 * ncomps + c);
        const auto *md_229 = buffer.data(md + 229 * ncomps + c);
        const auto *md_230 = buffer.data(md + 230 * ncomps + c);
        const auto *md_231 = buffer.data(md + 231 * ncomps + c);
        const auto *md_232 = buffer.data(md + 232 * ncomps + c);
        const auto *md_233 = buffer.data(md + 233 * ncomps + c);
        const auto *md_234 = buffer.data(md + 234 * ncomps + c);
        const auto *md_235 = buffer.data(md + 235 * ncomps + c);
        const auto *md_236 = buffer.data(md + 236 * ncomps + c);
        const auto *md_237 = buffer.data(md + 237 * ncomps + c);
        const auto *md_238 = buffer.data(md + 238 * ncomps + c);
        const auto *md_239 = buffer.data(md + 239 * ncomps + c);
        const auto *md_240 = buffer.data(md + 240 * ncomps + c);
        const auto *md_241 = buffer.data(md + 241 * ncomps + c);
        const auto *md_242 = buffer.data(md + 242 * ncomps + c);
        const auto *md_243 = buffer.data(md + 243 * ncomps + c);
        const auto *md_244 = buffer.data(md + 244 * ncomps + c);
        const auto *md_245 = buffer.data(md + 245 * ncomps + c);
        const auto *md_246 = buffer.data(md + 246 * ncomps + c);
        const auto *md_247 = buffer.data(md + 247 * ncomps + c);
        const auto *md_248 = buffer.data(md + 248 * ncomps + c);
        const auto *md_249 = buffer.data(md + 249 * ncomps + c);
        const auto *md_250 = buffer.data(md + 250 * ncomps + c);
        const auto *md_251 = buffer.data(md + 251 * ncomps + c);
        const auto *md_252 = buffer.data(md + 252 * ncomps + c);
        const auto *md_253 = buffer.data(md + 253 * ncomps + c);
        const auto *md_254 = buffer.data(md + 254 * ncomps + c);
        const auto *md_255 = buffer.data(md + 255 * ncomps + c);
        const auto *md_256 = buffer.data(md + 256 * ncomps + c);
        const auto *md_257 = buffer.data(md + 257 * ncomps + c);
        const auto *md_258 = buffer.data(md + 258 * ncomps + c);
        const auto *md_259 = buffer.data(md + 259 * ncomps + c);
        const auto *md_260 = buffer.data(md + 260 * ncomps + c);
        const auto *md_261 = buffer.data(md + 261 * ncomps + c);
        const auto *md_262 = buffer.data(md + 262 * ncomps + c);

        const auto *nd_174 = buffer.data(nd + 174 * ncomps + c);
        const auto *nd_175 = buffer.data(nd + 175 * ncomps + c);
        const auto *nd_176 = buffer.data(nd + 176 * ncomps + c);
        const auto *nd_177 = buffer.data(nd + 177 * ncomps + c);
        const auto *nd_178 = buffer.data(nd + 178 * ncomps + c);
        const auto *nd_179 = buffer.data(nd + 179 * ncomps + c);
        const auto *nd_180 = buffer.data(nd + 180 * ncomps + c);
        const auto *nd_181 = buffer.data(nd + 181 * ncomps + c);
        const auto *nd_182 = buffer.data(nd + 182 * ncomps + c);
        const auto *nd_183 = buffer.data(nd + 183 * ncomps + c);
        const auto *nd_184 = buffer.data(nd + 184 * ncomps + c);
        const auto *nd_185 = buffer.data(nd + 185 * ncomps + c);
        const auto *nd_186 = buffer.data(nd + 186 * ncomps + c);
        const auto *nd_187 = buffer.data(nd + 187 * ncomps + c);
        const auto *nd_188 = buffer.data(nd + 188 * ncomps + c);
        const auto *nd_189 = buffer.data(nd + 189 * ncomps + c);
        const auto *nd_190 = buffer.data(nd + 190 * ncomps + c);
        const auto *nd_191 = buffer.data(nd + 191 * ncomps + c);
        const auto *nd_192 = buffer.data(nd + 192 * ncomps + c);
        const auto *nd_193 = buffer.data(nd + 193 * ncomps + c);
        const auto *nd_194 = buffer.data(nd + 194 * ncomps + c);
        const auto *nd_195 = buffer.data(nd + 195 * ncomps + c);
        const auto *nd_196 = buffer.data(nd + 196 * ncomps + c);
        const auto *nd_197 = buffer.data(nd + 197 * ncomps + c);
        const auto *nd_198 = buffer.data(nd + 198 * ncomps + c);
        const auto *nd_199 = buffer.data(nd + 199 * ncomps + c);
        const auto *nd_200 = buffer.data(nd + 200 * ncomps + c);
        const auto *nd_201 = buffer.data(nd + 201 * ncomps + c);
        const auto *nd_202 = buffer.data(nd + 202 * ncomps + c);
        const auto *nd_203 = buffer.data(nd + 203 * ncomps + c);
        const auto *nd_204 = buffer.data(nd + 204 * ncomps + c);
        const auto *nd_205 = buffer.data(nd + 205 * ncomps + c);
        const auto *nd_206 = buffer.data(nd + 206 * ncomps + c);
        const auto *nd_207 = buffer.data(nd + 207 * ncomps + c);
        const auto *nd_208 = buffer.data(nd + 208 * ncomps + c);
        const auto *nd_209 = buffer.data(nd + 209 * ncomps + c);
        const auto *nd_210 = buffer.data(nd + 210 * ncomps + c);
        const auto *nd_211 = buffer.data(nd + 211 * ncomps + c);
        const auto *nd_212 = buffer.data(nd + 212 * ncomps + c);
        const auto *nd_213 = buffer.data(nd + 213 * ncomps + c);
        const auto *nd_214 = buffer.data(nd + 214 * ncomps + c);
        const auto *nd_215 = buffer.data(nd + 215 * ncomps + c);
        const auto *nd_216 = buffer.data(nd + 216 * ncomps + c);
        const auto *nd_217 = buffer.data(nd + 217 * ncomps + c);
        const auto *nd_218 = buffer.data(nd + 218 * ncomps + c);
        const auto *nd_219 = buffer.data(nd + 219 * ncomps + c);
        const auto *nd_220 = buffer.data(nd + 220 * ncomps + c);
        const auto *nd_221 = buffer.data(nd + 221 * ncomps + c);
        const auto *nd_222 = buffer.data(nd + 222 * ncomps + c);
        const auto *nd_223 = buffer.data(nd + 223 * ncomps + c);
        const auto *nd_224 = buffer.data(nd + 224 * ncomps + c);
        const auto *nd_225 = buffer.data(nd + 225 * ncomps + c);
        const auto *nd_226 = buffer.data(nd + 226 * ncomps + c);
        const auto *nd_227 = buffer.data(nd + 227 * ncomps + c);
        const auto *nd_228 = buffer.data(nd + 228 * ncomps + c);
        const auto *nd_229 = buffer.data(nd + 229 * ncomps + c);
        const auto *nd_230 = buffer.data(nd + 230 * ncomps + c);
        const auto *nd_231 = buffer.data(nd + 231 * ncomps + c);
        const auto *nd_232 = buffer.data(nd + 232 * ncomps + c);
        const auto *nd_233 = buffer.data(nd + 233 * ncomps + c);
        const auto *nd_234 = buffer.data(nd + 234 * ncomps + c);
        const auto *nd_235 = buffer.data(nd + 235 * ncomps + c);
        const auto *nd_236 = buffer.data(nd + 236 * ncomps + c);
        const auto *nd_237 = buffer.data(nd + 237 * ncomps + c);
        const auto *nd_238 = buffer.data(nd + 238 * ncomps + c);
        const auto *nd_239 = buffer.data(nd + 239 * ncomps + c);
        const auto *nd_240 = buffer.data(nd + 240 * ncomps + c);
        const auto *nd_241 = buffer.data(nd + 241 * ncomps + c);
        const auto *nd_242 = buffer.data(nd + 242 * ncomps + c);
        const auto *nd_243 = buffer.data(nd + 243 * ncomps + c);
        const auto *nd_244 = buffer.data(nd + 244 * ncomps + c);
        const auto *nd_245 = buffer.data(nd + 245 * ncomps + c);
        const auto *nd_246 = buffer.data(nd + 246 * ncomps + c);
        const auto *nd_247 = buffer.data(nd + 247 * ncomps + c);
        const auto *nd_248 = buffer.data(nd + 248 * ncomps + c);
        const auto *nd_249 = buffer.data(nd + 249 * ncomps + c);
        const auto *nd_250 = buffer.data(nd + 250 * ncomps + c);
        const auto *nd_251 = buffer.data(nd + 251 * ncomps + c);
        const auto *nd_252 = buffer.data(nd + 252 * ncomps + c);
        const auto *nd_253 = buffer.data(nd + 253 * ncomps + c);
        const auto *nd_254 = buffer.data(nd + 254 * ncomps + c);
        const auto *nd_255 = buffer.data(nd + 255 * ncomps + c);
        const auto *nd_256 = buffer.data(nd + 256 * ncomps + c);
        const auto *nd_257 = buffer.data(nd + 257 * ncomps + c);
        const auto *nd_258 = buffer.data(nd + 258 * ncomps + c);
        const auto *nd_259 = buffer.data(nd + 259 * ncomps + c);
        const auto *nd_260 = buffer.data(nd + 260 * ncomps + c);
        const auto *nd_261 = buffer.data(nd + 261 * ncomps + c);
        const auto *nd_262 = buffer.data(nd + 262 * ncomps + c);
        const auto *nd_263 = buffer.data(nd + 263 * ncomps + c);
        const auto *nd_269 = buffer.data(nd + 269 * ncomps + c);
        const auto *nd_273 = buffer.data(nd + 273 * ncomps + c);
        const auto *nd_274 = buffer.data(nd + 274 * ncomps + c);
        const auto *nd_275 = buffer.data(nd + 275 * ncomps + c);
        const auto *nd_279 = buffer.data(nd + 279 * ncomps + c);
        const auto *nd_280 = buffer.data(nd + 280 * ncomps + c);
        const auto *nd_281 = buffer.data(nd + 281 * ncomps + c);
        const auto *nd_285 = buffer.data(nd + 285 * ncomps + c);
        const auto *nd_286 = buffer.data(nd + 286 * ncomps + c);
        const auto *nd_287 = buffer.data(nd + 287 * ncomps + c);
        const auto *nd_291 = buffer.data(nd + 291 * ncomps + c);
        const auto *nd_292 = buffer.data(nd + 292 * ncomps + c);
        const auto *nd_293 = buffer.data(nd + 293 * ncomps + c);
        const auto *nd_297 = buffer.data(nd + 297 * ncomps + c);
        const auto *nd_298 = buffer.data(nd + 298 * ncomps + c);
        const auto *nd_299 = buffer.data(nd + 299 * ncomps + c);
        const auto *nd_303 = buffer.data(nd + 303 * ncomps + c);
        const auto *nd_304 = buffer.data(nd + 304 * ncomps + c);
        const auto *nd_305 = buffer.data(nd + 305 * ncomps + c);
        const auto *nd_309 = buffer.data(nd + 309 * ncomps + c);
        const auto *nd_310 = buffer.data(nd + 310 * ncomps + c);
        const auto *nd_311 = buffer.data(nd + 311 * ncomps + c);
        const auto *nd_317 = buffer.data(nd + 317 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, md_174, md_175, md_176, \
                         md_177, md_178, nd_174, nd_175, nd_176, nd_177, \
                         nd_178 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = ab_x[k] * md_174[k]
                       + nd_174[k];

            t_291[k] = ab_x[k] * md_175[k]
                       + nd_175[k];

            t_292[k] = ab_x[k] * md_176[k]
                       + nd_176[k];

            t_293[k] = ab_x[k] * md_177[k]
                       + nd_177[k];

            t_294[k] = ab_x[k] * md_178[k]
                       + nd_178[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_x, ab_y, ab_z, md_177, md_178, \
                         md_179, nd_179, nd_225, nd_226, nd_227, \
                         nd_233 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = ab_x[k] * md_179[k]
                       + nd_179[k];

            t_296[k] = ab_y[k] * md_177[k]
                       + nd_225[k];

            t_297[k] = ab_y[k] * md_178[k]
                       + nd_226[k];

            t_298[k] = ab_y[k] * md_179[k]
                       + nd_227[k];

            t_299[k] = ab_z[k] * md_179[k]
                       + nd_233[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, md_180, md_181, md_182, \
                         md_183, md_184, nd_180, nd_181, nd_182, nd_183, \
                         nd_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = ab_x[k] * md_180[k]
                       + nd_180[k];

            t_301[k] = ab_x[k] * md_181[k]
                       + nd_181[k];

            t_302[k] = ab_x[k] * md_182[k]
                       + nd_182[k];

            t_303[k] = ab_x[k] * md_183[k]
                       + nd_183[k];

            t_304[k] = ab_x[k] * md_184[k]
                       + nd_184[k];
        }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, ab_y, ab_z, md_183, md_184, \
                         md_185, nd_185, nd_231, nd_232, nd_233, \
                         nd_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_305[k] = ab_x[k] * md_185[k]
                       + nd_185[k];

            t_306[k] = ab_y[k] * md_183[k]
                       + nd_231[k];

            t_307[k] = ab_y[k] * md_184[k]
                       + nd_232[k];

            t_308[k] = ab_y[k] * md_185[k]
                       + nd_233[k];

            t_309[k] = ab_z[k] * md_185[k]
                       + nd_239[k];
        }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_x, md_186, md_187, md_188, \
                         md_189, md_190, nd_186, nd_187, nd_188, nd_189, \
                         nd_190 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_310[k] = ab_x[k] * md_186[k]
                       + nd_186[k];

            t_311[k] = ab_x[k] * md_187[k]
                       + nd_187[k];

            t_312[k] = ab_x[k] * md_188[k]
                       + nd_188[k];

            t_313[k] = ab_x[k] * md_189[k]
                       + nd_189[k];

            t_314[k] = ab_x[k] * md_190[k]
                       + nd_190[k];
        }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_x, ab_y, ab_z, md_189, md_190, \
                         md_191, nd_191, nd_237, nd_238, nd_239, \
                         nd_245 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_315[k] = ab_x[k] * md_191[k]
                       + nd_191[k];

            t_316[k] = ab_y[k] * md_189[k]
                       + nd_237[k];

            t_317[k] = ab_y[k] * md_190[k]
                       + nd_238[k];

            t_318[k] = ab_y[k] * md_191[k]
                       + nd_239[k];

            t_319[k] = ab_z[k] * md_191[k]
                       + nd_245[k];
        }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_x, md_192, md_193, md_194, \
                         md_195, md_196, nd_192, nd_193, nd_194, nd_195, \
                         nd_196 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_320[k] = ab_x[k] * md_192[k]
                       + nd_192[k];

            t_321[k] = ab_x[k] * md_193[k]
                       + nd_193[k];

            t_322[k] = ab_x[k] * md_194[k]
                       + nd_194[k];

            t_323[k] = ab_x[k] * md_195[k]
                       + nd_195[k];

            t_324[k] = ab_x[k] * md_196[k]
                       + nd_196[k];
        }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_x, ab_y, ab_z, md_195, md_196, \
                         md_197, nd_197, nd_243, nd_244, nd_245, \
                         nd_251 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_325[k] = ab_x[k] * md_197[k]
                       + nd_197[k];

            t_326[k] = ab_y[k] * md_195[k]
                       + nd_243[k];

            t_327[k] = ab_y[k] * md_196[k]
                       + nd_244[k];

            t_328[k] = ab_y[k] * md_197[k]
                       + nd_245[k];

            t_329[k] = ab_z[k] * md_197[k]
                       + nd_251[k];
        }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_x, md_198, md_199, md_200, \
                         md_201, md_202, nd_198, nd_199, nd_200, nd_201, \
                         nd_202 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_330[k] = ab_x[k] * md_198[k]
                       + nd_198[k];

            t_331[k] = ab_x[k] * md_199[k]
                       + nd_199[k];

            t_332[k] = ab_x[k] * md_200[k]
                       + nd_200[k];

            t_333[k] = ab_x[k] * md_201[k]
                       + nd_201[k];

            t_334[k] = ab_x[k] * md_202[k]
                       + nd_202[k];
        }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_x, ab_y, ab_z, md_201, md_202, \
                         md_203, nd_203, nd_249, nd_250, nd_251, \
                         nd_257 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_335[k] = ab_x[k] * md_203[k]
                       + nd_203[k];

            t_336[k] = ab_y[k] * md_201[k]
                       + nd_249[k];

            t_337[k] = ab_y[k] * md_202[k]
                       + nd_250[k];

            t_338[k] = ab_y[k] * md_203[k]
                       + nd_251[k];

            t_339[k] = ab_z[k] * md_203[k]
                       + nd_257[k];
        }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_x, md_204, md_205, md_206, \
                         md_207, md_208, nd_204, nd_205, nd_206, nd_207, \
                         nd_208 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_340[k] = ab_x[k] * md_204[k]
                       + nd_204[k];

            t_341[k] = ab_x[k] * md_205[k]
                       + nd_205[k];

            t_342[k] = ab_x[k] * md_206[k]
                       + nd_206[k];

            t_343[k] = ab_x[k] * md_207[k]
                       + nd_207[k];

            t_344[k] = ab_x[k] * md_208[k]
                       + nd_208[k];
        }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_x, ab_y, ab_z, md_207, md_208, \
                         md_209, nd_209, nd_255, nd_256, nd_257, \
                         nd_263 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_345[k] = ab_x[k] * md_209[k]
                       + nd_209[k];

            t_346[k] = ab_y[k] * md_207[k]
                       + nd_255[k];

            t_347[k] = ab_y[k] * md_208[k]
                       + nd_256[k];

            t_348[k] = ab_y[k] * md_209[k]
                       + nd_257[k];

            t_349[k] = ab_z[k] * md_209[k]
                       + nd_263[k];
        }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_x, md_210, md_211, md_212, \
                         md_213, md_214, nd_210, nd_211, nd_212, nd_213, \
                         nd_214 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_350[k] = ab_x[k] * md_210[k]
                       + nd_210[k];

            t_351[k] = ab_x[k] * md_211[k]
                       + nd_211[k];

            t_352[k] = ab_x[k] * md_212[k]
                       + nd_212[k];

            t_353[k] = ab_x[k] * md_213[k]
                       + nd_213[k];

            t_354[k] = ab_x[k] * md_214[k]
                       + nd_214[k];
        }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_x, ab_y, ab_z, md_213, md_214, \
                         md_215, nd_215, nd_261, nd_262, nd_263, \
                         nd_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_355[k] = ab_x[k] * md_215[k]
                       + nd_215[k];

            t_356[k] = ab_y[k] * md_213[k]
                       + nd_261[k];

            t_357[k] = ab_y[k] * md_214[k]
                       + nd_262[k];

            t_358[k] = ab_y[k] * md_215[k]
                       + nd_263[k];

            t_359[k] = ab_z[k] * md_215[k]
                       + nd_269[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, md_216, md_217, md_218, \
                         md_219, md_220, nd_216, nd_217, nd_218, nd_219, \
                         nd_220 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = ab_x[k] * md_216[k]
                       + nd_216[k];

            t_361[k] = ab_x[k] * md_217[k]
                       + nd_217[k];

            t_362[k] = ab_x[k] * md_218[k]
                       + nd_218[k];

            t_363[k] = ab_x[k] * md_219[k]
                       + nd_219[k];

            t_364[k] = ab_x[k] * md_220[k]
                       + nd_220[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, ab_y, ab_z, md_219, md_220, \
                         md_221, nd_221, nd_273, nd_274, nd_275, \
                         nd_281 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = ab_x[k] * md_221[k]
                       + nd_221[k];

            t_366[k] = ab_y[k] * md_219[k]
                       + nd_273[k];

            t_367[k] = ab_y[k] * md_220[k]
                       + nd_274[k];

            t_368[k] = ab_y[k] * md_221[k]
                       + nd_275[k];

            t_369[k] = ab_z[k] * md_221[k]
                       + nd_281[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_x, md_222, md_223, md_224, \
                         md_225, md_226, nd_222, nd_223, nd_224, nd_225, \
                         nd_226 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = ab_x[k] * md_222[k]
                       + nd_222[k];

            t_371[k] = ab_x[k] * md_223[k]
                       + nd_223[k];

            t_372[k] = ab_x[k] * md_224[k]
                       + nd_224[k];

            t_373[k] = ab_x[k] * md_225[k]
                       + nd_225[k];

            t_374[k] = ab_x[k] * md_226[k]
                       + nd_226[k];
        }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_x, ab_y, ab_z, md_225, md_226, \
                         md_227, nd_227, nd_279, nd_280, nd_281, \
                         nd_287 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_375[k] = ab_x[k] * md_227[k]
                       + nd_227[k];

            t_376[k] = ab_y[k] * md_225[k]
                       + nd_279[k];

            t_377[k] = ab_y[k] * md_226[k]
                       + nd_280[k];

            t_378[k] = ab_y[k] * md_227[k]
                       + nd_281[k];

            t_379[k] = ab_z[k] * md_227[k]
                       + nd_287[k];
        }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_x, md_228, md_229, md_230, \
                         md_231, md_232, nd_228, nd_229, nd_230, nd_231, \
                         nd_232 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_380[k] = ab_x[k] * md_228[k]
                       + nd_228[k];

            t_381[k] = ab_x[k] * md_229[k]
                       + nd_229[k];

            t_382[k] = ab_x[k] * md_230[k]
                       + nd_230[k];

            t_383[k] = ab_x[k] * md_231[k]
                       + nd_231[k];

            t_384[k] = ab_x[k] * md_232[k]
                       + nd_232[k];
        }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_x, ab_y, ab_z, md_231, md_232, \
                         md_233, nd_233, nd_285, nd_286, nd_287, \
                         nd_293 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_385[k] = ab_x[k] * md_233[k]
                       + nd_233[k];

            t_386[k] = ab_y[k] * md_231[k]
                       + nd_285[k];

            t_387[k] = ab_y[k] * md_232[k]
                       + nd_286[k];

            t_388[k] = ab_y[k] * md_233[k]
                       + nd_287[k];

            t_389[k] = ab_z[k] * md_233[k]
                       + nd_293[k];
        }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_x, md_234, md_235, md_236, \
                         md_237, md_238, nd_234, nd_235, nd_236, nd_237, \
                         nd_238 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_390[k] = ab_x[k] * md_234[k]
                       + nd_234[k];

            t_391[k] = ab_x[k] * md_235[k]
                       + nd_235[k];

            t_392[k] = ab_x[k] * md_236[k]
                       + nd_236[k];

            t_393[k] = ab_x[k] * md_237[k]
                       + nd_237[k];

            t_394[k] = ab_x[k] * md_238[k]
                       + nd_238[k];
        }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_x, ab_y, ab_z, md_237, md_238, \
                         md_239, nd_239, nd_291, nd_292, nd_293, \
                         nd_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_395[k] = ab_x[k] * md_239[k]
                       + nd_239[k];

            t_396[k] = ab_y[k] * md_237[k]
                       + nd_291[k];

            t_397[k] = ab_y[k] * md_238[k]
                       + nd_292[k];

            t_398[k] = ab_y[k] * md_239[k]
                       + nd_293[k];

            t_399[k] = ab_z[k] * md_239[k]
                       + nd_299[k];
        }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_x, md_240, md_241, md_242, \
                         md_243, md_244, nd_240, nd_241, nd_242, nd_243, \
                         nd_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_400[k] = ab_x[k] * md_240[k]
                       + nd_240[k];

            t_401[k] = ab_x[k] * md_241[k]
                       + nd_241[k];

            t_402[k] = ab_x[k] * md_242[k]
                       + nd_242[k];

            t_403[k] = ab_x[k] * md_243[k]
                       + nd_243[k];

            t_404[k] = ab_x[k] * md_244[k]
                       + nd_244[k];
        }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_x, ab_y, ab_z, md_243, md_244, \
                         md_245, nd_245, nd_297, nd_298, nd_299, \
                         nd_305 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_405[k] = ab_x[k] * md_245[k]
                       + nd_245[k];

            t_406[k] = ab_y[k] * md_243[k]
                       + nd_297[k];

            t_407[k] = ab_y[k] * md_244[k]
                       + nd_298[k];

            t_408[k] = ab_y[k] * md_245[k]
                       + nd_299[k];

            t_409[k] = ab_z[k] * md_245[k]
                       + nd_305[k];
        }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_x, md_246, md_247, md_248, \
                         md_249, md_250, nd_246, nd_247, nd_248, nd_249, \
                         nd_250 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_410[k] = ab_x[k] * md_246[k]
                       + nd_246[k];

            t_411[k] = ab_x[k] * md_247[k]
                       + nd_247[k];

            t_412[k] = ab_x[k] * md_248[k]
                       + nd_248[k];

            t_413[k] = ab_x[k] * md_249[k]
                       + nd_249[k];

            t_414[k] = ab_x[k] * md_250[k]
                       + nd_250[k];
        }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_x, ab_y, ab_z, md_249, md_250, \
                         md_251, nd_251, nd_303, nd_304, nd_305, \
                         nd_311 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_415[k] = ab_x[k] * md_251[k]
                       + nd_251[k];

            t_416[k] = ab_y[k] * md_249[k]
                       + nd_303[k];

            t_417[k] = ab_y[k] * md_250[k]
                       + nd_304[k];

            t_418[k] = ab_y[k] * md_251[k]
                       + nd_305[k];

            t_419[k] = ab_z[k] * md_251[k]
                       + nd_311[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_x, md_252, md_253, md_254, \
                         md_255, md_256, nd_252, nd_253, nd_254, nd_255, \
                         nd_256 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = ab_x[k] * md_252[k]
                       + nd_252[k];

            t_421[k] = ab_x[k] * md_253[k]
                       + nd_253[k];

            t_422[k] = ab_x[k] * md_254[k]
                       + nd_254[k];

            t_423[k] = ab_x[k] * md_255[k]
                       + nd_255[k];

            t_424[k] = ab_x[k] * md_256[k]
                       + nd_256[k];
        }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_x, ab_y, ab_z, md_255, md_256, \
                         md_257, nd_257, nd_309, nd_310, nd_311, \
                         nd_317 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = ab_x[k] * md_257[k]
                       + nd_257[k];

            t_426[k] = ab_y[k] * md_255[k]
                       + nd_309[k];

            t_427[k] = ab_y[k] * md_256[k]
                       + nd_310[k];

            t_428[k] = ab_y[k] * md_257[k]
                       + nd_311[k];

            t_429[k] = ab_z[k] * md_257[k]
                       + nd_317[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_x, md_258, md_259, md_260, \
                         md_261, md_262, nd_258, nd_259, nd_260, nd_261, \
                         nd_262 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = ab_x[k] * md_258[k]
                       + nd_258[k];

            t_431[k] = ab_x[k] * md_259[k]
                       + nd_259[k];

            t_432[k] = ab_x[k] * md_260[k]
                       + nd_260[k];

            t_433[k] = ab_x[k] * md_261[k]
                       + nd_261[k];

            t_434[k] = ab_x[k] * md_262[k]
                       + nd_262[k];
        }
    }
}

static auto
compute_hrr_mf_out_of_first_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t md, const size_t nd,
                                   const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_435 = buffer.data(target + 435 * ncomps + c);
        auto *t_436 = buffer.data(target + 436 * ncomps + c);
        auto *t_437 = buffer.data(target + 437 * ncomps + c);
        auto *t_438 = buffer.data(target + 438 * ncomps + c);
        auto *t_439 = buffer.data(target + 439 * ncomps + c);
        auto *t_440 = buffer.data(target + 440 * ncomps + c);
        auto *t_441 = buffer.data(target + 441 * ncomps + c);
        auto *t_442 = buffer.data(target + 442 * ncomps + c);
        auto *t_443 = buffer.data(target + 443 * ncomps + c);
        auto *t_444 = buffer.data(target + 444 * ncomps + c);
        auto *t_445 = buffer.data(target + 445 * ncomps + c);
        auto *t_446 = buffer.data(target + 446 * ncomps + c);
        auto *t_447 = buffer.data(target + 447 * ncomps + c);
        auto *t_448 = buffer.data(target + 448 * ncomps + c);
        auto *t_449 = buffer.data(target + 449 * ncomps + c);
        auto *t_450 = buffer.data(target + 450 * ncomps + c);
        auto *t_451 = buffer.data(target + 451 * ncomps + c);
        auto *t_452 = buffer.data(target + 452 * ncomps + c);
        auto *t_453 = buffer.data(target + 453 * ncomps + c);
        auto *t_454 = buffer.data(target + 454 * ncomps + c);
        auto *t_455 = buffer.data(target + 455 * ncomps + c);
        auto *t_456 = buffer.data(target + 456 * ncomps + c);
        auto *t_457 = buffer.data(target + 457 * ncomps + c);
        auto *t_458 = buffer.data(target + 458 * ncomps + c);
        auto *t_459 = buffer.data(target + 459 * ncomps + c);
        auto *t_460 = buffer.data(target + 460 * ncomps + c);
        auto *t_461 = buffer.data(target + 461 * ncomps + c);
        auto *t_462 = buffer.data(target + 462 * ncomps + c);
        auto *t_463 = buffer.data(target + 463 * ncomps + c);
        auto *t_464 = buffer.data(target + 464 * ncomps + c);
        auto *t_465 = buffer.data(target + 465 * ncomps + c);
        auto *t_466 = buffer.data(target + 466 * ncomps + c);
        auto *t_467 = buffer.data(target + 467 * ncomps + c);
        auto *t_468 = buffer.data(target + 468 * ncomps + c);
        auto *t_469 = buffer.data(target + 469 * ncomps + c);
        auto *t_470 = buffer.data(target + 470 * ncomps + c);
        auto *t_471 = buffer.data(target + 471 * ncomps + c);
        auto *t_472 = buffer.data(target + 472 * ncomps + c);
        auto *t_473 = buffer.data(target + 473 * ncomps + c);
        auto *t_474 = buffer.data(target + 474 * ncomps + c);
        auto *t_475 = buffer.data(target + 475 * ncomps + c);
        auto *t_476 = buffer.data(target + 476 * ncomps + c);
        auto *t_477 = buffer.data(target + 477 * ncomps + c);
        auto *t_478 = buffer.data(target + 478 * ncomps + c);
        auto *t_479 = buffer.data(target + 479 * ncomps + c);
        auto *t_480 = buffer.data(target + 480 * ncomps + c);
        auto *t_481 = buffer.data(target + 481 * ncomps + c);
        auto *t_482 = buffer.data(target + 482 * ncomps + c);
        auto *t_483 = buffer.data(target + 483 * ncomps + c);
        auto *t_484 = buffer.data(target + 484 * ncomps + c);
        auto *t_485 = buffer.data(target + 485 * ncomps + c);
        auto *t_486 = buffer.data(target + 486 * ncomps + c);
        auto *t_487 = buffer.data(target + 487 * ncomps + c);
        auto *t_488 = buffer.data(target + 488 * ncomps + c);
        auto *t_489 = buffer.data(target + 489 * ncomps + c);
        auto *t_490 = buffer.data(target + 490 * ncomps + c);
        auto *t_491 = buffer.data(target + 491 * ncomps + c);
        auto *t_492 = buffer.data(target + 492 * ncomps + c);
        auto *t_493 = buffer.data(target + 493 * ncomps + c);
        auto *t_494 = buffer.data(target + 494 * ncomps + c);
        auto *t_495 = buffer.data(target + 495 * ncomps + c);
        auto *t_496 = buffer.data(target + 496 * ncomps + c);
        auto *t_497 = buffer.data(target + 497 * ncomps + c);
        auto *t_498 = buffer.data(target + 498 * ncomps + c);
        auto *t_499 = buffer.data(target + 499 * ncomps + c);
        auto *t_500 = buffer.data(target + 500 * ncomps + c);
        auto *t_501 = buffer.data(target + 501 * ncomps + c);
        auto *t_502 = buffer.data(target + 502 * ncomps + c);
        auto *t_503 = buffer.data(target + 503 * ncomps + c);
        auto *t_504 = buffer.data(target + 504 * ncomps + c);
        auto *t_505 = buffer.data(target + 505 * ncomps + c);
        auto *t_506 = buffer.data(target + 506 * ncomps + c);
        auto *t_507 = buffer.data(target + 507 * ncomps + c);
        auto *t_508 = buffer.data(target + 508 * ncomps + c);
        auto *t_509 = buffer.data(target + 509 * ncomps + c);
        auto *t_510 = buffer.data(target + 510 * ncomps + c);
        auto *t_511 = buffer.data(target + 511 * ncomps + c);
        auto *t_512 = buffer.data(target + 512 * ncomps + c);
        auto *t_513 = buffer.data(target + 513 * ncomps + c);
        auto *t_514 = buffer.data(target + 514 * ncomps + c);
        auto *t_515 = buffer.data(target + 515 * ncomps + c);
        auto *t_516 = buffer.data(target + 516 * ncomps + c);
        auto *t_517 = buffer.data(target + 517 * ncomps + c);
        auto *t_518 = buffer.data(target + 518 * ncomps + c);
        auto *t_519 = buffer.data(target + 519 * ncomps + c);
        auto *t_520 = buffer.data(target + 520 * ncomps + c);
        auto *t_521 = buffer.data(target + 521 * ncomps + c);
        auto *t_522 = buffer.data(target + 522 * ncomps + c);
        auto *t_523 = buffer.data(target + 523 * ncomps + c);
        auto *t_524 = buffer.data(target + 524 * ncomps + c);
        auto *t_525 = buffer.data(target + 525 * ncomps + c);
        auto *t_526 = buffer.data(target + 526 * ncomps + c);
        auto *t_527 = buffer.data(target + 527 * ncomps + c);
        auto *t_528 = buffer.data(target + 528 * ncomps + c);
        auto *t_529 = buffer.data(target + 529 * ncomps + c);
        auto *t_530 = buffer.data(target + 530 * ncomps + c);
        auto *t_531 = buffer.data(target + 531 * ncomps + c);
        auto *t_532 = buffer.data(target + 532 * ncomps + c);
        auto *t_533 = buffer.data(target + 533 * ncomps + c);
        auto *t_534 = buffer.data(target + 534 * ncomps + c);
        auto *t_535 = buffer.data(target + 535 * ncomps + c);
        auto *t_536 = buffer.data(target + 536 * ncomps + c);
        auto *t_537 = buffer.data(target + 537 * ncomps + c);
        auto *t_538 = buffer.data(target + 538 * ncomps + c);
        auto *t_539 = buffer.data(target + 539 * ncomps + c);
        auto *t_540 = buffer.data(target + 540 * ncomps + c);
        auto *t_541 = buffer.data(target + 541 * ncomps + c);
        auto *t_542 = buffer.data(target + 542 * ncomps + c);
        auto *t_543 = buffer.data(target + 543 * ncomps + c);
        auto *t_544 = buffer.data(target + 544 * ncomps + c);
        auto *t_545 = buffer.data(target + 545 * ncomps + c);
        auto *t_546 = buffer.data(target + 546 * ncomps + c);
        auto *t_547 = buffer.data(target + 547 * ncomps + c);
        auto *t_548 = buffer.data(target + 548 * ncomps + c);
        auto *t_549 = buffer.data(target + 549 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *md_261 = buffer.data(md + 261 * ncomps + c);
        const auto *md_262 = buffer.data(md + 262 * ncomps + c);
        const auto *md_263 = buffer.data(md + 263 * ncomps + c);
        const auto *md_264 = buffer.data(md + 264 * ncomps + c);
        const auto *md_265 = buffer.data(md + 265 * ncomps + c);
        const auto *md_266 = buffer.data(md + 266 * ncomps + c);
        const auto *md_267 = buffer.data(md + 267 * ncomps + c);
        const auto *md_268 = buffer.data(md + 268 * ncomps + c);
        const auto *md_269 = buffer.data(md + 269 * ncomps + c);
        const auto *md_270 = buffer.data(md + 270 * ncomps + c);
        const auto *md_271 = buffer.data(md + 271 * ncomps + c);
        const auto *md_272 = buffer.data(md + 272 * ncomps + c);
        const auto *md_273 = buffer.data(md + 273 * ncomps + c);
        const auto *md_274 = buffer.data(md + 274 * ncomps + c);
        const auto *md_275 = buffer.data(md + 275 * ncomps + c);
        const auto *md_276 = buffer.data(md + 276 * ncomps + c);
        const auto *md_277 = buffer.data(md + 277 * ncomps + c);
        const auto *md_278 = buffer.data(md + 278 * ncomps + c);
        const auto *md_279 = buffer.data(md + 279 * ncomps + c);
        const auto *md_280 = buffer.data(md + 280 * ncomps + c);
        const auto *md_281 = buffer.data(md + 281 * ncomps + c);
        const auto *md_282 = buffer.data(md + 282 * ncomps + c);
        const auto *md_283 = buffer.data(md + 283 * ncomps + c);
        const auto *md_284 = buffer.data(md + 284 * ncomps + c);
        const auto *md_285 = buffer.data(md + 285 * ncomps + c);
        const auto *md_286 = buffer.data(md + 286 * ncomps + c);
        const auto *md_287 = buffer.data(md + 287 * ncomps + c);
        const auto *md_288 = buffer.data(md + 288 * ncomps + c);
        const auto *md_289 = buffer.data(md + 289 * ncomps + c);
        const auto *md_290 = buffer.data(md + 290 * ncomps + c);
        const auto *md_291 = buffer.data(md + 291 * ncomps + c);
        const auto *md_292 = buffer.data(md + 292 * ncomps + c);
        const auto *md_293 = buffer.data(md + 293 * ncomps + c);
        const auto *md_294 = buffer.data(md + 294 * ncomps + c);
        const auto *md_295 = buffer.data(md + 295 * ncomps + c);
        const auto *md_296 = buffer.data(md + 296 * ncomps + c);
        const auto *md_297 = buffer.data(md + 297 * ncomps + c);
        const auto *md_298 = buffer.data(md + 298 * ncomps + c);
        const auto *md_299 = buffer.data(md + 299 * ncomps + c);
        const auto *md_300 = buffer.data(md + 300 * ncomps + c);
        const auto *md_301 = buffer.data(md + 301 * ncomps + c);
        const auto *md_302 = buffer.data(md + 302 * ncomps + c);
        const auto *md_303 = buffer.data(md + 303 * ncomps + c);
        const auto *md_304 = buffer.data(md + 304 * ncomps + c);
        const auto *md_305 = buffer.data(md + 305 * ncomps + c);
        const auto *md_306 = buffer.data(md + 306 * ncomps + c);
        const auto *md_307 = buffer.data(md + 307 * ncomps + c);
        const auto *md_308 = buffer.data(md + 308 * ncomps + c);
        const auto *md_309 = buffer.data(md + 309 * ncomps + c);
        const auto *md_310 = buffer.data(md + 310 * ncomps + c);
        const auto *md_311 = buffer.data(md + 311 * ncomps + c);
        const auto *md_312 = buffer.data(md + 312 * ncomps + c);
        const auto *md_313 = buffer.data(md + 313 * ncomps + c);
        const auto *md_314 = buffer.data(md + 314 * ncomps + c);
        const auto *md_315 = buffer.data(md + 315 * ncomps + c);
        const auto *md_316 = buffer.data(md + 316 * ncomps + c);
        const auto *md_317 = buffer.data(md + 317 * ncomps + c);
        const auto *md_318 = buffer.data(md + 318 * ncomps + c);
        const auto *md_319 = buffer.data(md + 319 * ncomps + c);
        const auto *md_320 = buffer.data(md + 320 * ncomps + c);
        const auto *md_321 = buffer.data(md + 321 * ncomps + c);
        const auto *md_322 = buffer.data(md + 322 * ncomps + c);
        const auto *md_323 = buffer.data(md + 323 * ncomps + c);
        const auto *md_324 = buffer.data(md + 324 * ncomps + c);
        const auto *md_325 = buffer.data(md + 325 * ncomps + c);
        const auto *md_326 = buffer.data(md + 326 * ncomps + c);
        const auto *md_327 = buffer.data(md + 327 * ncomps + c);
        const auto *md_328 = buffer.data(md + 328 * ncomps + c);
        const auto *md_329 = buffer.data(md + 329 * ncomps + c);

        const auto *nd_263 = buffer.data(nd + 263 * ncomps + c);
        const auto *nd_264 = buffer.data(nd + 264 * ncomps + c);
        const auto *nd_265 = buffer.data(nd + 265 * ncomps + c);
        const auto *nd_266 = buffer.data(nd + 266 * ncomps + c);
        const auto *nd_267 = buffer.data(nd + 267 * ncomps + c);
        const auto *nd_268 = buffer.data(nd + 268 * ncomps + c);
        const auto *nd_269 = buffer.data(nd + 269 * ncomps + c);
        const auto *nd_270 = buffer.data(nd + 270 * ncomps + c);
        const auto *nd_271 = buffer.data(nd + 271 * ncomps + c);
        const auto *nd_272 = buffer.data(nd + 272 * ncomps + c);
        const auto *nd_273 = buffer.data(nd + 273 * ncomps + c);
        const auto *nd_274 = buffer.data(nd + 274 * ncomps + c);
        const auto *nd_275 = buffer.data(nd + 275 * ncomps + c);
        const auto *nd_276 = buffer.data(nd + 276 * ncomps + c);
        const auto *nd_277 = buffer.data(nd + 277 * ncomps + c);
        const auto *nd_278 = buffer.data(nd + 278 * ncomps + c);
        const auto *nd_279 = buffer.data(nd + 279 * ncomps + c);
        const auto *nd_280 = buffer.data(nd + 280 * ncomps + c);
        const auto *nd_281 = buffer.data(nd + 281 * ncomps + c);
        const auto *nd_282 = buffer.data(nd + 282 * ncomps + c);
        const auto *nd_283 = buffer.data(nd + 283 * ncomps + c);
        const auto *nd_284 = buffer.data(nd + 284 * ncomps + c);
        const auto *nd_285 = buffer.data(nd + 285 * ncomps + c);
        const auto *nd_286 = buffer.data(nd + 286 * ncomps + c);
        const auto *nd_287 = buffer.data(nd + 287 * ncomps + c);
        const auto *nd_288 = buffer.data(nd + 288 * ncomps + c);
        const auto *nd_289 = buffer.data(nd + 289 * ncomps + c);
        const auto *nd_290 = buffer.data(nd + 290 * ncomps + c);
        const auto *nd_291 = buffer.data(nd + 291 * ncomps + c);
        const auto *nd_292 = buffer.data(nd + 292 * ncomps + c);
        const auto *nd_293 = buffer.data(nd + 293 * ncomps + c);
        const auto *nd_294 = buffer.data(nd + 294 * ncomps + c);
        const auto *nd_295 = buffer.data(nd + 295 * ncomps + c);
        const auto *nd_296 = buffer.data(nd + 296 * ncomps + c);
        const auto *nd_297 = buffer.data(nd + 297 * ncomps + c);
        const auto *nd_298 = buffer.data(nd + 298 * ncomps + c);
        const auto *nd_299 = buffer.data(nd + 299 * ncomps + c);
        const auto *nd_300 = buffer.data(nd + 300 * ncomps + c);
        const auto *nd_301 = buffer.data(nd + 301 * ncomps + c);
        const auto *nd_302 = buffer.data(nd + 302 * ncomps + c);
        const auto *nd_303 = buffer.data(nd + 303 * ncomps + c);
        const auto *nd_304 = buffer.data(nd + 304 * ncomps + c);
        const auto *nd_305 = buffer.data(nd + 305 * ncomps + c);
        const auto *nd_306 = buffer.data(nd + 306 * ncomps + c);
        const auto *nd_307 = buffer.data(nd + 307 * ncomps + c);
        const auto *nd_308 = buffer.data(nd + 308 * ncomps + c);
        const auto *nd_309 = buffer.data(nd + 309 * ncomps + c);
        const auto *nd_310 = buffer.data(nd + 310 * ncomps + c);
        const auto *nd_311 = buffer.data(nd + 311 * ncomps + c);
        const auto *nd_312 = buffer.data(nd + 312 * ncomps + c);
        const auto *nd_313 = buffer.data(nd + 313 * ncomps + c);
        const auto *nd_314 = buffer.data(nd + 314 * ncomps + c);
        const auto *nd_315 = buffer.data(nd + 315 * ncomps + c);
        const auto *nd_316 = buffer.data(nd + 316 * ncomps + c);
        const auto *nd_317 = buffer.data(nd + 317 * ncomps + c);
        const auto *nd_318 = buffer.data(nd + 318 * ncomps + c);
        const auto *nd_319 = buffer.data(nd + 319 * ncomps + c);
        const auto *nd_320 = buffer.data(nd + 320 * ncomps + c);
        const auto *nd_321 = buffer.data(nd + 321 * ncomps + c);
        const auto *nd_322 = buffer.data(nd + 322 * ncomps + c);
        const auto *nd_323 = buffer.data(nd + 323 * ncomps + c);
        const auto *nd_324 = buffer.data(nd + 324 * ncomps + c);
        const auto *nd_325 = buffer.data(nd + 325 * ncomps + c);
        const auto *nd_326 = buffer.data(nd + 326 * ncomps + c);
        const auto *nd_327 = buffer.data(nd + 327 * ncomps + c);
        const auto *nd_328 = buffer.data(nd + 328 * ncomps + c);
        const auto *nd_329 = buffer.data(nd + 329 * ncomps + c);
        const auto *nd_333 = buffer.data(nd + 333 * ncomps + c);
        const auto *nd_334 = buffer.data(nd + 334 * ncomps + c);
        const auto *nd_335 = buffer.data(nd + 335 * ncomps + c);
        const auto *nd_339 = buffer.data(nd + 339 * ncomps + c);
        const auto *nd_340 = buffer.data(nd + 340 * ncomps + c);
        const auto *nd_341 = buffer.data(nd + 341 * ncomps + c);
        const auto *nd_345 = buffer.data(nd + 345 * ncomps + c);
        const auto *nd_346 = buffer.data(nd + 346 * ncomps + c);
        const auto *nd_347 = buffer.data(nd + 347 * ncomps + c);
        const auto *nd_351 = buffer.data(nd + 351 * ncomps + c);
        const auto *nd_352 = buffer.data(nd + 352 * ncomps + c);
        const auto *nd_353 = buffer.data(nd + 353 * ncomps + c);
        const auto *nd_357 = buffer.data(nd + 357 * ncomps + c);
        const auto *nd_358 = buffer.data(nd + 358 * ncomps + c);
        const auto *nd_359 = buffer.data(nd + 359 * ncomps + c);
        const auto *nd_363 = buffer.data(nd + 363 * ncomps + c);
        const auto *nd_364 = buffer.data(nd + 364 * ncomps + c);
        const auto *nd_365 = buffer.data(nd + 365 * ncomps + c);
        const auto *nd_369 = buffer.data(nd + 369 * ncomps + c);
        const auto *nd_370 = buffer.data(nd + 370 * ncomps + c);
        const auto *nd_371 = buffer.data(nd + 371 * ncomps + c);
        const auto *nd_375 = buffer.data(nd + 375 * ncomps + c);
        const auto *nd_376 = buffer.data(nd + 376 * ncomps + c);
        const auto *nd_377 = buffer.data(nd + 377 * ncomps + c);
        const auto *nd_381 = buffer.data(nd + 381 * ncomps + c);
        const auto *nd_382 = buffer.data(nd + 382 * ncomps + c);
        const auto *nd_383 = buffer.data(nd + 383 * ncomps + c);
        const auto *nd_387 = buffer.data(nd + 387 * ncomps + c);
        const auto *nd_388 = buffer.data(nd + 388 * ncomps + c);
        const auto *nd_389 = buffer.data(nd + 389 * ncomps + c);
        const auto *nd_395 = buffer.data(nd + 395 * ncomps + c);

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_x, ab_y, ab_z, md_261, md_262, \
                         md_263, nd_263, nd_315, nd_316, nd_317, \
                         nd_323 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = ab_x[k] * md_263[k]
                       + nd_263[k];

            t_436[k] = ab_y[k] * md_261[k]
                       + nd_315[k];

            t_437[k] = ab_y[k] * md_262[k]
                       + nd_316[k];

            t_438[k] = ab_y[k] * md_263[k]
                       + nd_317[k];

            t_439[k] = ab_z[k] * md_263[k]
                       + nd_323[k];
        }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ab_x, md_264, md_265, md_266, \
                         md_267, md_268, nd_264, nd_265, nd_266, nd_267, \
                         nd_268 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = ab_x[k] * md_264[k]
                       + nd_264[k];

            t_441[k] = ab_x[k] * md_265[k]
                       + nd_265[k];

            t_442[k] = ab_x[k] * md_266[k]
                       + nd_266[k];

            t_443[k] = ab_x[k] * md_267[k]
                       + nd_267[k];

            t_444[k] = ab_x[k] * md_268[k]
                       + nd_268[k];
        }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ab_x, ab_y, ab_z, md_267, md_268, \
                         md_269, nd_269, nd_321, nd_322, nd_323, \
                         nd_329 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_445[k] = ab_x[k] * md_269[k]
                       + nd_269[k];

            t_446[k] = ab_y[k] * md_267[k]
                       + nd_321[k];

            t_447[k] = ab_y[k] * md_268[k]
                       + nd_322[k];

            t_448[k] = ab_y[k] * md_269[k]
                       + nd_323[k];

            t_449[k] = ab_z[k] * md_269[k]
                       + nd_329[k];
        }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ab_x, md_270, md_271, md_272, \
                         md_273, md_274, nd_270, nd_271, nd_272, nd_273, \
                         nd_274 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_450[k] = ab_x[k] * md_270[k]
                       + nd_270[k];

            t_451[k] = ab_x[k] * md_271[k]
                       + nd_271[k];

            t_452[k] = ab_x[k] * md_272[k]
                       + nd_272[k];

            t_453[k] = ab_x[k] * md_273[k]
                       + nd_273[k];

            t_454[k] = ab_x[k] * md_274[k]
                       + nd_274[k];
        }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ab_x, ab_y, ab_z, md_273, md_274, \
                         md_275, nd_275, nd_333, nd_334, nd_335, \
                         nd_341 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_455[k] = ab_x[k] * md_275[k]
                       + nd_275[k];

            t_456[k] = ab_y[k] * md_273[k]
                       + nd_333[k];

            t_457[k] = ab_y[k] * md_274[k]
                       + nd_334[k];

            t_458[k] = ab_y[k] * md_275[k]
                       + nd_335[k];

            t_459[k] = ab_z[k] * md_275[k]
                       + nd_341[k];
        }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ab_x, md_276, md_277, md_278, \
                         md_279, md_280, nd_276, nd_277, nd_278, nd_279, \
                         nd_280 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_460[k] = ab_x[k] * md_276[k]
                       + nd_276[k];

            t_461[k] = ab_x[k] * md_277[k]
                       + nd_277[k];

            t_462[k] = ab_x[k] * md_278[k]
                       + nd_278[k];

            t_463[k] = ab_x[k] * md_279[k]
                       + nd_279[k];

            t_464[k] = ab_x[k] * md_280[k]
                       + nd_280[k];
        }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ab_x, ab_y, ab_z, md_279, md_280, \
                         md_281, nd_281, nd_339, nd_340, nd_341, \
                         nd_347 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_465[k] = ab_x[k] * md_281[k]
                       + nd_281[k];

            t_466[k] = ab_y[k] * md_279[k]
                       + nd_339[k];

            t_467[k] = ab_y[k] * md_280[k]
                       + nd_340[k];

            t_468[k] = ab_y[k] * md_281[k]
                       + nd_341[k];

            t_469[k] = ab_z[k] * md_281[k]
                       + nd_347[k];
        }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ab_x, md_282, md_283, md_284, \
                         md_285, md_286, nd_282, nd_283, nd_284, nd_285, \
                         nd_286 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_470[k] = ab_x[k] * md_282[k]
                       + nd_282[k];

            t_471[k] = ab_x[k] * md_283[k]
                       + nd_283[k];

            t_472[k] = ab_x[k] * md_284[k]
                       + nd_284[k];

            t_473[k] = ab_x[k] * md_285[k]
                       + nd_285[k];

            t_474[k] = ab_x[k] * md_286[k]
                       + nd_286[k];
        }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ab_x, ab_y, ab_z, md_285, md_286, \
                         md_287, nd_287, nd_345, nd_346, nd_347, \
                         nd_353 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_475[k] = ab_x[k] * md_287[k]
                       + nd_287[k];

            t_476[k] = ab_y[k] * md_285[k]
                       + nd_345[k];

            t_477[k] = ab_y[k] * md_286[k]
                       + nd_346[k];

            t_478[k] = ab_y[k] * md_287[k]
                       + nd_347[k];

            t_479[k] = ab_z[k] * md_287[k]
                       + nd_353[k];
        }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ab_x, md_288, md_289, md_290, \
                         md_291, md_292, nd_288, nd_289, nd_290, nd_291, \
                         nd_292 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_480[k] = ab_x[k] * md_288[k]
                       + nd_288[k];

            t_481[k] = ab_x[k] * md_289[k]
                       + nd_289[k];

            t_482[k] = ab_x[k] * md_290[k]
                       + nd_290[k];

            t_483[k] = ab_x[k] * md_291[k]
                       + nd_291[k];

            t_484[k] = ab_x[k] * md_292[k]
                       + nd_292[k];
        }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ab_x, ab_y, ab_z, md_291, md_292, \
                         md_293, nd_293, nd_351, nd_352, nd_353, \
                         nd_359 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_485[k] = ab_x[k] * md_293[k]
                       + nd_293[k];

            t_486[k] = ab_y[k] * md_291[k]
                       + nd_351[k];

            t_487[k] = ab_y[k] * md_292[k]
                       + nd_352[k];

            t_488[k] = ab_y[k] * md_293[k]
                       + nd_353[k];

            t_489[k] = ab_z[k] * md_293[k]
                       + nd_359[k];
        }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ab_x, md_294, md_295, md_296, \
                         md_297, md_298, nd_294, nd_295, nd_296, nd_297, \
                         nd_298 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_490[k] = ab_x[k] * md_294[k]
                       + nd_294[k];

            t_491[k] = ab_x[k] * md_295[k]
                       + nd_295[k];

            t_492[k] = ab_x[k] * md_296[k]
                       + nd_296[k];

            t_493[k] = ab_x[k] * md_297[k]
                       + nd_297[k];

            t_494[k] = ab_x[k] * md_298[k]
                       + nd_298[k];
        }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ab_x, ab_y, ab_z, md_297, md_298, \
                         md_299, nd_299, nd_357, nd_358, nd_359, \
                         nd_365 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_495[k] = ab_x[k] * md_299[k]
                       + nd_299[k];

            t_496[k] = ab_y[k] * md_297[k]
                       + nd_357[k];

            t_497[k] = ab_y[k] * md_298[k]
                       + nd_358[k];

            t_498[k] = ab_y[k] * md_299[k]
                       + nd_359[k];

            t_499[k] = ab_z[k] * md_299[k]
                       + nd_365[k];
        }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ab_x, md_300, md_301, md_302, \
                         md_303, md_304, nd_300, nd_301, nd_302, nd_303, \
                         nd_304 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_500[k] = ab_x[k] * md_300[k]
                       + nd_300[k];

            t_501[k] = ab_x[k] * md_301[k]
                       + nd_301[k];

            t_502[k] = ab_x[k] * md_302[k]
                       + nd_302[k];

            t_503[k] = ab_x[k] * md_303[k]
                       + nd_303[k];

            t_504[k] = ab_x[k] * md_304[k]
                       + nd_304[k];
        }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ab_x, ab_y, ab_z, md_303, md_304, \
                         md_305, nd_305, nd_363, nd_364, nd_365, \
                         nd_371 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_505[k] = ab_x[k] * md_305[k]
                       + nd_305[k];

            t_506[k] = ab_y[k] * md_303[k]
                       + nd_363[k];

            t_507[k] = ab_y[k] * md_304[k]
                       + nd_364[k];

            t_508[k] = ab_y[k] * md_305[k]
                       + nd_365[k];

            t_509[k] = ab_z[k] * md_305[k]
                       + nd_371[k];
        }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ab_x, md_306, md_307, md_308, \
                         md_309, md_310, nd_306, nd_307, nd_308, nd_309, \
                         nd_310 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_510[k] = ab_x[k] * md_306[k]
                       + nd_306[k];

            t_511[k] = ab_x[k] * md_307[k]
                       + nd_307[k];

            t_512[k] = ab_x[k] * md_308[k]
                       + nd_308[k];

            t_513[k] = ab_x[k] * md_309[k]
                       + nd_309[k];

            t_514[k] = ab_x[k] * md_310[k]
                       + nd_310[k];
        }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ab_x, ab_y, ab_z, md_309, md_310, \
                         md_311, nd_311, nd_369, nd_370, nd_371, \
                         nd_377 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_515[k] = ab_x[k] * md_311[k]
                       + nd_311[k];

            t_516[k] = ab_y[k] * md_309[k]
                       + nd_369[k];

            t_517[k] = ab_y[k] * md_310[k]
                       + nd_370[k];

            t_518[k] = ab_y[k] * md_311[k]
                       + nd_371[k];

            t_519[k] = ab_z[k] * md_311[k]
                       + nd_377[k];
        }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ab_x, md_312, md_313, md_314, \
                         md_315, md_316, nd_312, nd_313, nd_314, nd_315, \
                         nd_316 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_520[k] = ab_x[k] * md_312[k]
                       + nd_312[k];

            t_521[k] = ab_x[k] * md_313[k]
                       + nd_313[k];

            t_522[k] = ab_x[k] * md_314[k]
                       + nd_314[k];

            t_523[k] = ab_x[k] * md_315[k]
                       + nd_315[k];

            t_524[k] = ab_x[k] * md_316[k]
                       + nd_316[k];
        }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, ab_x, ab_y, ab_z, md_315, md_316, \
                         md_317, nd_317, nd_375, nd_376, nd_377, \
                         nd_383 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_525[k] = ab_x[k] * md_317[k]
                       + nd_317[k];

            t_526[k] = ab_y[k] * md_315[k]
                       + nd_375[k];

            t_527[k] = ab_y[k] * md_316[k]
                       + nd_376[k];

            t_528[k] = ab_y[k] * md_317[k]
                       + nd_377[k];

            t_529[k] = ab_z[k] * md_317[k]
                       + nd_383[k];
        }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, ab_x, md_318, md_319, md_320, \
                         md_321, md_322, nd_318, nd_319, nd_320, nd_321, \
                         nd_322 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_530[k] = ab_x[k] * md_318[k]
                       + nd_318[k];

            t_531[k] = ab_x[k] * md_319[k]
                       + nd_319[k];

            t_532[k] = ab_x[k] * md_320[k]
                       + nd_320[k];

            t_533[k] = ab_x[k] * md_321[k]
                       + nd_321[k];

            t_534[k] = ab_x[k] * md_322[k]
                       + nd_322[k];
        }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, ab_x, ab_y, ab_z, md_321, md_322, \
                         md_323, nd_323, nd_381, nd_382, nd_383, \
                         nd_389 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_535[k] = ab_x[k] * md_323[k]
                       + nd_323[k];

            t_536[k] = ab_y[k] * md_321[k]
                       + nd_381[k];

            t_537[k] = ab_y[k] * md_322[k]
                       + nd_382[k];

            t_538[k] = ab_y[k] * md_323[k]
                       + nd_383[k];

            t_539[k] = ab_z[k] * md_323[k]
                       + nd_389[k];
        }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, ab_x, md_324, md_325, md_326, \
                         md_327, md_328, nd_324, nd_325, nd_326, nd_327, \
                         nd_328 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_540[k] = ab_x[k] * md_324[k]
                       + nd_324[k];

            t_541[k] = ab_x[k] * md_325[k]
                       + nd_325[k];

            t_542[k] = ab_x[k] * md_326[k]
                       + nd_326[k];

            t_543[k] = ab_x[k] * md_327[k]
                       + nd_327[k];

            t_544[k] = ab_x[k] * md_328[k]
                       + nd_328[k];
        }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, ab_x, ab_y, ab_z, md_327, md_328, \
                         md_329, nd_329, nd_387, nd_388, nd_389, \
                         nd_395 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_545[k] = ab_x[k] * md_329[k]
                       + nd_329[k];

            t_546[k] = ab_y[k] * md_327[k]
                       + nd_387[k];

            t_547[k] = ab_y[k] * md_328[k]
                       + nd_388[k];

            t_548[k] = ab_y[k] * md_329[k]
                       + nd_389[k];

            t_549[k] = ab_z[k] * md_329[k]
                       + nd_395[k];
        }
    }
}

auto
compute_hrr_mf_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t md, const size_t nd,
                            const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_mf_out_of_first_piece0(buffer, coordinates, target, md, nd, ncomps, nmax);

    compute_hrr_mf_out_of_first_piece1(buffer, coordinates, target, md, nd, ncomps, nmax);

    compute_hrr_mf_out_of_first_piece2(buffer, coordinates, target, md, nd, ncomps, nmax);

    compute_hrr_mf_out_of_first_piece3(buffer, coordinates, target, md, nd, ncomps, nmax);
}

}  // namespace simdtrf
