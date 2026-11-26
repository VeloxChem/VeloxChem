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


#include "SimdTransferLG.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_lg_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t lf, const size_t mf,
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

        const auto *lf_0 = buffer.data(lf + 0 * ncomps + c);
        const auto *lf_1 = buffer.data(lf + 1 * ncomps + c);
        const auto *lf_2 = buffer.data(lf + 2 * ncomps + c);
        const auto *lf_3 = buffer.data(lf + 3 * ncomps + c);
        const auto *lf_4 = buffer.data(lf + 4 * ncomps + c);
        const auto *lf_5 = buffer.data(lf + 5 * ncomps + c);
        const auto *lf_6 = buffer.data(lf + 6 * ncomps + c);
        const auto *lf_7 = buffer.data(lf + 7 * ncomps + c);
        const auto *lf_8 = buffer.data(lf + 8 * ncomps + c);
        const auto *lf_9 = buffer.data(lf + 9 * ncomps + c);
        const auto *lf_10 = buffer.data(lf + 10 * ncomps + c);
        const auto *lf_11 = buffer.data(lf + 11 * ncomps + c);
        const auto *lf_12 = buffer.data(lf + 12 * ncomps + c);
        const auto *lf_13 = buffer.data(lf + 13 * ncomps + c);
        const auto *lf_14 = buffer.data(lf + 14 * ncomps + c);
        const auto *lf_15 = buffer.data(lf + 15 * ncomps + c);
        const auto *lf_16 = buffer.data(lf + 16 * ncomps + c);
        const auto *lf_17 = buffer.data(lf + 17 * ncomps + c);
        const auto *lf_18 = buffer.data(lf + 18 * ncomps + c);
        const auto *lf_19 = buffer.data(lf + 19 * ncomps + c);
        const auto *lf_20 = buffer.data(lf + 20 * ncomps + c);
        const auto *lf_21 = buffer.data(lf + 21 * ncomps + c);
        const auto *lf_22 = buffer.data(lf + 22 * ncomps + c);
        const auto *lf_23 = buffer.data(lf + 23 * ncomps + c);
        const auto *lf_24 = buffer.data(lf + 24 * ncomps + c);
        const auto *lf_25 = buffer.data(lf + 25 * ncomps + c);
        const auto *lf_26 = buffer.data(lf + 26 * ncomps + c);
        const auto *lf_27 = buffer.data(lf + 27 * ncomps + c);
        const auto *lf_28 = buffer.data(lf + 28 * ncomps + c);
        const auto *lf_29 = buffer.data(lf + 29 * ncomps + c);
        const auto *lf_30 = buffer.data(lf + 30 * ncomps + c);
        const auto *lf_31 = buffer.data(lf + 31 * ncomps + c);
        const auto *lf_32 = buffer.data(lf + 32 * ncomps + c);
        const auto *lf_33 = buffer.data(lf + 33 * ncomps + c);
        const auto *lf_34 = buffer.data(lf + 34 * ncomps + c);
        const auto *lf_35 = buffer.data(lf + 35 * ncomps + c);
        const auto *lf_36 = buffer.data(lf + 36 * ncomps + c);
        const auto *lf_37 = buffer.data(lf + 37 * ncomps + c);
        const auto *lf_38 = buffer.data(lf + 38 * ncomps + c);
        const auto *lf_39 = buffer.data(lf + 39 * ncomps + c);
        const auto *lf_40 = buffer.data(lf + 40 * ncomps + c);
        const auto *lf_41 = buffer.data(lf + 41 * ncomps + c);
        const auto *lf_42 = buffer.data(lf + 42 * ncomps + c);
        const auto *lf_43 = buffer.data(lf + 43 * ncomps + c);
        const auto *lf_44 = buffer.data(lf + 44 * ncomps + c);
        const auto *lf_45 = buffer.data(lf + 45 * ncomps + c);
        const auto *lf_46 = buffer.data(lf + 46 * ncomps + c);
        const auto *lf_47 = buffer.data(lf + 47 * ncomps + c);
        const auto *lf_48 = buffer.data(lf + 48 * ncomps + c);
        const auto *lf_49 = buffer.data(lf + 49 * ncomps + c);
        const auto *lf_50 = buffer.data(lf + 50 * ncomps + c);
        const auto *lf_51 = buffer.data(lf + 51 * ncomps + c);
        const auto *lf_52 = buffer.data(lf + 52 * ncomps + c);
        const auto *lf_53 = buffer.data(lf + 53 * ncomps + c);
        const auto *lf_54 = buffer.data(lf + 54 * ncomps + c);
        const auto *lf_55 = buffer.data(lf + 55 * ncomps + c);
        const auto *lf_56 = buffer.data(lf + 56 * ncomps + c);
        const auto *lf_57 = buffer.data(lf + 57 * ncomps + c);
        const auto *lf_58 = buffer.data(lf + 58 * ncomps + c);
        const auto *lf_59 = buffer.data(lf + 59 * ncomps + c);
        const auto *lf_60 = buffer.data(lf + 60 * ncomps + c);
        const auto *lf_61 = buffer.data(lf + 61 * ncomps + c);
        const auto *lf_62 = buffer.data(lf + 62 * ncomps + c);
        const auto *lf_63 = buffer.data(lf + 63 * ncomps + c);
        const auto *lf_64 = buffer.data(lf + 64 * ncomps + c);
        const auto *lf_65 = buffer.data(lf + 65 * ncomps + c);
        const auto *lf_66 = buffer.data(lf + 66 * ncomps + c);
        const auto *lf_67 = buffer.data(lf + 67 * ncomps + c);
        const auto *lf_68 = buffer.data(lf + 68 * ncomps + c);
        const auto *lf_69 = buffer.data(lf + 69 * ncomps + c);
        const auto *lf_70 = buffer.data(lf + 70 * ncomps + c);
        const auto *lf_71 = buffer.data(lf + 71 * ncomps + c);
        const auto *lf_72 = buffer.data(lf + 72 * ncomps + c);
        const auto *lf_73 = buffer.data(lf + 73 * ncomps + c);
        const auto *lf_74 = buffer.data(lf + 74 * ncomps + c);
        const auto *lf_75 = buffer.data(lf + 75 * ncomps + c);
        const auto *lf_76 = buffer.data(lf + 76 * ncomps + c);
        const auto *lf_77 = buffer.data(lf + 77 * ncomps + c);
        const auto *lf_78 = buffer.data(lf + 78 * ncomps + c);
        const auto *lf_79 = buffer.data(lf + 79 * ncomps + c);
        const auto *lf_80 = buffer.data(lf + 80 * ncomps + c);
        const auto *lf_81 = buffer.data(lf + 81 * ncomps + c);
        const auto *lf_82 = buffer.data(lf + 82 * ncomps + c);
        const auto *lf_83 = buffer.data(lf + 83 * ncomps + c);
        const auto *lf_84 = buffer.data(lf + 84 * ncomps + c);
        const auto *lf_85 = buffer.data(lf + 85 * ncomps + c);
        const auto *lf_86 = buffer.data(lf + 86 * ncomps + c);
        const auto *lf_87 = buffer.data(lf + 87 * ncomps + c);
        const auto *lf_88 = buffer.data(lf + 88 * ncomps + c);
        const auto *lf_89 = buffer.data(lf + 89 * ncomps + c);
        const auto *lf_90 = buffer.data(lf + 90 * ncomps + c);
        const auto *lf_91 = buffer.data(lf + 91 * ncomps + c);
        const auto *lf_92 = buffer.data(lf + 92 * ncomps + c);
        const auto *lf_93 = buffer.data(lf + 93 * ncomps + c);
        const auto *lf_94 = buffer.data(lf + 94 * ncomps + c);
        const auto *lf_95 = buffer.data(lf + 95 * ncomps + c);
        const auto *lf_96 = buffer.data(lf + 96 * ncomps + c);
        const auto *lf_97 = buffer.data(lf + 97 * ncomps + c);
        const auto *lf_98 = buffer.data(lf + 98 * ncomps + c);
        const auto *lf_99 = buffer.data(lf + 99 * ncomps + c);

        const auto *mf_0 = buffer.data(mf + 0 * ncomps + c);
        const auto *mf_1 = buffer.data(mf + 1 * ncomps + c);
        const auto *mf_2 = buffer.data(mf + 2 * ncomps + c);
        const auto *mf_3 = buffer.data(mf + 3 * ncomps + c);
        const auto *mf_4 = buffer.data(mf + 4 * ncomps + c);
        const auto *mf_5 = buffer.data(mf + 5 * ncomps + c);
        const auto *mf_6 = buffer.data(mf + 6 * ncomps + c);
        const auto *mf_7 = buffer.data(mf + 7 * ncomps + c);
        const auto *mf_8 = buffer.data(mf + 8 * ncomps + c);
        const auto *mf_9 = buffer.data(mf + 9 * ncomps + c);
        const auto *mf_10 = buffer.data(mf + 10 * ncomps + c);
        const auto *mf_11 = buffer.data(mf + 11 * ncomps + c);
        const auto *mf_12 = buffer.data(mf + 12 * ncomps + c);
        const auto *mf_13 = buffer.data(mf + 13 * ncomps + c);
        const auto *mf_14 = buffer.data(mf + 14 * ncomps + c);
        const auto *mf_15 = buffer.data(mf + 15 * ncomps + c);
        const auto *mf_16 = buffer.data(mf + 16 * ncomps + c);
        const auto *mf_17 = buffer.data(mf + 17 * ncomps + c);
        const auto *mf_18 = buffer.data(mf + 18 * ncomps + c);
        const auto *mf_19 = buffer.data(mf + 19 * ncomps + c);
        const auto *mf_20 = buffer.data(mf + 20 * ncomps + c);
        const auto *mf_21 = buffer.data(mf + 21 * ncomps + c);
        const auto *mf_22 = buffer.data(mf + 22 * ncomps + c);
        const auto *mf_23 = buffer.data(mf + 23 * ncomps + c);
        const auto *mf_24 = buffer.data(mf + 24 * ncomps + c);
        const auto *mf_25 = buffer.data(mf + 25 * ncomps + c);
        const auto *mf_26 = buffer.data(mf + 26 * ncomps + c);
        const auto *mf_27 = buffer.data(mf + 27 * ncomps + c);
        const auto *mf_28 = buffer.data(mf + 28 * ncomps + c);
        const auto *mf_29 = buffer.data(mf + 29 * ncomps + c);
        const auto *mf_30 = buffer.data(mf + 30 * ncomps + c);
        const auto *mf_31 = buffer.data(mf + 31 * ncomps + c);
        const auto *mf_32 = buffer.data(mf + 32 * ncomps + c);
        const auto *mf_33 = buffer.data(mf + 33 * ncomps + c);
        const auto *mf_34 = buffer.data(mf + 34 * ncomps + c);
        const auto *mf_35 = buffer.data(mf + 35 * ncomps + c);
        const auto *mf_36 = buffer.data(mf + 36 * ncomps + c);
        const auto *mf_37 = buffer.data(mf + 37 * ncomps + c);
        const auto *mf_38 = buffer.data(mf + 38 * ncomps + c);
        const auto *mf_39 = buffer.data(mf + 39 * ncomps + c);
        const auto *mf_40 = buffer.data(mf + 40 * ncomps + c);
        const auto *mf_41 = buffer.data(mf + 41 * ncomps + c);
        const auto *mf_42 = buffer.data(mf + 42 * ncomps + c);
        const auto *mf_43 = buffer.data(mf + 43 * ncomps + c);
        const auto *mf_44 = buffer.data(mf + 44 * ncomps + c);
        const auto *mf_45 = buffer.data(mf + 45 * ncomps + c);
        const auto *mf_46 = buffer.data(mf + 46 * ncomps + c);
        const auto *mf_47 = buffer.data(mf + 47 * ncomps + c);
        const auto *mf_48 = buffer.data(mf + 48 * ncomps + c);
        const auto *mf_49 = buffer.data(mf + 49 * ncomps + c);
        const auto *mf_50 = buffer.data(mf + 50 * ncomps + c);
        const auto *mf_51 = buffer.data(mf + 51 * ncomps + c);
        const auto *mf_52 = buffer.data(mf + 52 * ncomps + c);
        const auto *mf_53 = buffer.data(mf + 53 * ncomps + c);
        const auto *mf_54 = buffer.data(mf + 54 * ncomps + c);
        const auto *mf_55 = buffer.data(mf + 55 * ncomps + c);
        const auto *mf_56 = buffer.data(mf + 56 * ncomps + c);
        const auto *mf_57 = buffer.data(mf + 57 * ncomps + c);
        const auto *mf_58 = buffer.data(mf + 58 * ncomps + c);
        const auto *mf_59 = buffer.data(mf + 59 * ncomps + c);
        const auto *mf_60 = buffer.data(mf + 60 * ncomps + c);
        const auto *mf_61 = buffer.data(mf + 61 * ncomps + c);
        const auto *mf_62 = buffer.data(mf + 62 * ncomps + c);
        const auto *mf_63 = buffer.data(mf + 63 * ncomps + c);
        const auto *mf_64 = buffer.data(mf + 64 * ncomps + c);
        const auto *mf_65 = buffer.data(mf + 65 * ncomps + c);
        const auto *mf_66 = buffer.data(mf + 66 * ncomps + c);
        const auto *mf_67 = buffer.data(mf + 67 * ncomps + c);
        const auto *mf_68 = buffer.data(mf + 68 * ncomps + c);
        const auto *mf_69 = buffer.data(mf + 69 * ncomps + c);
        const auto *mf_70 = buffer.data(mf + 70 * ncomps + c);
        const auto *mf_71 = buffer.data(mf + 71 * ncomps + c);
        const auto *mf_72 = buffer.data(mf + 72 * ncomps + c);
        const auto *mf_73 = buffer.data(mf + 73 * ncomps + c);
        const auto *mf_74 = buffer.data(mf + 74 * ncomps + c);
        const auto *mf_75 = buffer.data(mf + 75 * ncomps + c);
        const auto *mf_76 = buffer.data(mf + 76 * ncomps + c);
        const auto *mf_77 = buffer.data(mf + 77 * ncomps + c);
        const auto *mf_78 = buffer.data(mf + 78 * ncomps + c);
        const auto *mf_79 = buffer.data(mf + 79 * ncomps + c);
        const auto *mf_80 = buffer.data(mf + 80 * ncomps + c);
        const auto *mf_81 = buffer.data(mf + 81 * ncomps + c);
        const auto *mf_82 = buffer.data(mf + 82 * ncomps + c);
        const auto *mf_83 = buffer.data(mf + 83 * ncomps + c);
        const auto *mf_84 = buffer.data(mf + 84 * ncomps + c);
        const auto *mf_85 = buffer.data(mf + 85 * ncomps + c);
        const auto *mf_86 = buffer.data(mf + 86 * ncomps + c);
        const auto *mf_87 = buffer.data(mf + 87 * ncomps + c);
        const auto *mf_88 = buffer.data(mf + 88 * ncomps + c);
        const auto *mf_89 = buffer.data(mf + 89 * ncomps + c);
        const auto *mf_90 = buffer.data(mf + 90 * ncomps + c);
        const auto *mf_91 = buffer.data(mf + 91 * ncomps + c);
        const auto *mf_92 = buffer.data(mf + 92 * ncomps + c);
        const auto *mf_93 = buffer.data(mf + 93 * ncomps + c);
        const auto *mf_94 = buffer.data(mf + 94 * ncomps + c);
        const auto *mf_95 = buffer.data(mf + 95 * ncomps + c);
        const auto *mf_96 = buffer.data(mf + 96 * ncomps + c);
        const auto *mf_97 = buffer.data(mf + 97 * ncomps + c);
        const auto *mf_98 = buffer.data(mf + 98 * ncomps + c);
        const auto *mf_99 = buffer.data(mf + 99 * ncomps + c);
        const auto *mf_106 = buffer.data(mf + 106 * ncomps + c);
        const auto *mf_107 = buffer.data(mf + 107 * ncomps + c);
        const auto *mf_108 = buffer.data(mf + 108 * ncomps + c);
        const auto *mf_109 = buffer.data(mf + 109 * ncomps + c);
        const auto *mf_116 = buffer.data(mf + 116 * ncomps + c);
        const auto *mf_117 = buffer.data(mf + 117 * ncomps + c);
        const auto *mf_118 = buffer.data(mf + 118 * ncomps + c);
        const auto *mf_119 = buffer.data(mf + 119 * ncomps + c);
        const auto *mf_126 = buffer.data(mf + 126 * ncomps + c);
        const auto *mf_127 = buffer.data(mf + 127 * ncomps + c);
        const auto *mf_128 = buffer.data(mf + 128 * ncomps + c);
        const auto *mf_129 = buffer.data(mf + 129 * ncomps + c);
        const auto *mf_139 = buffer.data(mf + 139 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, lf_0, lf_1, lf_2, lf_3, lf_4, mf_0, \
                         mf_1, mf_2, mf_3, mf_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * lf_0[k]
                     + mf_0[k];

            t_1[k] = ab_x[k] * lf_1[k]
                     + mf_1[k];

            t_2[k] = ab_x[k] * lf_2[k]
                     + mf_2[k];

            t_3[k] = ab_x[k] * lf_3[k]
                     + mf_3[k];

            t_4[k] = ab_x[k] * lf_4[k]
                     + mf_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, lf_5, lf_6, lf_7, lf_8, lf_9, mf_5, \
                         mf_6, mf_7, mf_8, mf_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * lf_5[k]
                     + mf_5[k];

            t_6[k] = ab_x[k] * lf_6[k]
                     + mf_6[k];

            t_7[k] = ab_x[k] * lf_7[k]
                     + mf_7[k];

            t_8[k] = ab_x[k] * lf_8[k]
                     + mf_8[k];

            t_9[k] = ab_x[k] * lf_9[k]
                     + mf_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_y, ab_z, lf_6, lf_7, lf_8, lf_9, \
                         mf_16, mf_17, mf_18, mf_19, mf_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_y[k] * lf_6[k]
                      + mf_16[k];

            t_11[k] = ab_y[k] * lf_7[k]
                      + mf_17[k];

            t_12[k] = ab_y[k] * lf_8[k]
                      + mf_18[k];

            t_13[k] = ab_y[k] * lf_9[k]
                      + mf_19[k];

            t_14[k] = ab_z[k] * lf_9[k]
                      + mf_29[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, lf_10, lf_11, lf_12, lf_13, \
                         lf_14, mf_10, mf_11, mf_12, mf_13, mf_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_x[k] * lf_10[k]
                      + mf_10[k];

            t_16[k] = ab_x[k] * lf_11[k]
                      + mf_11[k];

            t_17[k] = ab_x[k] * lf_12[k]
                      + mf_12[k];

            t_18[k] = ab_x[k] * lf_13[k]
                      + mf_13[k];

            t_19[k] = ab_x[k] * lf_14[k]
                      + mf_14[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, lf_15, lf_16, lf_17, lf_18, \
                         lf_19, mf_15, mf_16, mf_17, mf_18, mf_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_x[k] * lf_15[k]
                      + mf_15[k];

            t_21[k] = ab_x[k] * lf_16[k]
                      + mf_16[k];

            t_22[k] = ab_x[k] * lf_17[k]
                      + mf_17[k];

            t_23[k] = ab_x[k] * lf_18[k]
                      + mf_18[k];

            t_24[k] = ab_x[k] * lf_19[k]
                      + mf_19[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_y, ab_z, lf_16, lf_17, lf_18, lf_19, \
                         mf_36, mf_37, mf_38, mf_39, mf_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = ab_y[k] * lf_16[k]
                      + mf_36[k];

            t_26[k] = ab_y[k] * lf_17[k]
                      + mf_37[k];

            t_27[k] = ab_y[k] * lf_18[k]
                      + mf_38[k];

            t_28[k] = ab_y[k] * lf_19[k]
                      + mf_39[k];

            t_29[k] = ab_z[k] * lf_19[k]
                      + mf_49[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, lf_20, lf_21, lf_22, lf_23, \
                         lf_24, mf_20, mf_21, mf_22, mf_23, mf_24 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = ab_x[k] * lf_20[k]
                      + mf_20[k];

            t_31[k] = ab_x[k] * lf_21[k]
                      + mf_21[k];

            t_32[k] = ab_x[k] * lf_22[k]
                      + mf_22[k];

            t_33[k] = ab_x[k] * lf_23[k]
                      + mf_23[k];

            t_34[k] = ab_x[k] * lf_24[k]
                      + mf_24[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, lf_25, lf_26, lf_27, lf_28, \
                         lf_29, mf_25, mf_26, mf_27, mf_28, mf_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = ab_x[k] * lf_25[k]
                      + mf_25[k];

            t_36[k] = ab_x[k] * lf_26[k]
                      + mf_26[k];

            t_37[k] = ab_x[k] * lf_27[k]
                      + mf_27[k];

            t_38[k] = ab_x[k] * lf_28[k]
                      + mf_28[k];

            t_39[k] = ab_x[k] * lf_29[k]
                      + mf_29[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_y, ab_z, lf_26, lf_27, lf_28, lf_29, \
                         mf_46, mf_47, mf_48, mf_49, mf_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = ab_y[k] * lf_26[k]
                      + mf_46[k];

            t_41[k] = ab_y[k] * lf_27[k]
                      + mf_47[k];

            t_42[k] = ab_y[k] * lf_28[k]
                      + mf_48[k];

            t_43[k] = ab_y[k] * lf_29[k]
                      + mf_49[k];

            t_44[k] = ab_z[k] * lf_29[k]
                      + mf_59[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, lf_30, lf_31, lf_32, lf_33, \
                         lf_34, mf_30, mf_31, mf_32, mf_33, mf_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = ab_x[k] * lf_30[k]
                      + mf_30[k];

            t_46[k] = ab_x[k] * lf_31[k]
                      + mf_31[k];

            t_47[k] = ab_x[k] * lf_32[k]
                      + mf_32[k];

            t_48[k] = ab_x[k] * lf_33[k]
                      + mf_33[k];

            t_49[k] = ab_x[k] * lf_34[k]
                      + mf_34[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, lf_35, lf_36, lf_37, lf_38, \
                         lf_39, mf_35, mf_36, mf_37, mf_38, mf_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = ab_x[k] * lf_35[k]
                      + mf_35[k];

            t_51[k] = ab_x[k] * lf_36[k]
                      + mf_36[k];

            t_52[k] = ab_x[k] * lf_37[k]
                      + mf_37[k];

            t_53[k] = ab_x[k] * lf_38[k]
                      + mf_38[k];

            t_54[k] = ab_x[k] * lf_39[k]
                      + mf_39[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_y, ab_z, lf_36, lf_37, lf_38, lf_39, \
                         mf_66, mf_67, mf_68, mf_69, mf_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = ab_y[k] * lf_36[k]
                      + mf_66[k];

            t_56[k] = ab_y[k] * lf_37[k]
                      + mf_67[k];

            t_57[k] = ab_y[k] * lf_38[k]
                      + mf_68[k];

            t_58[k] = ab_y[k] * lf_39[k]
                      + mf_69[k];

            t_59[k] = ab_z[k] * lf_39[k]
                      + mf_79[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, lf_40, lf_41, lf_42, lf_43, \
                         lf_44, mf_40, mf_41, mf_42, mf_43, mf_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = ab_x[k] * lf_40[k]
                      + mf_40[k];

            t_61[k] = ab_x[k] * lf_41[k]
                      + mf_41[k];

            t_62[k] = ab_x[k] * lf_42[k]
                      + mf_42[k];

            t_63[k] = ab_x[k] * lf_43[k]
                      + mf_43[k];

            t_64[k] = ab_x[k] * lf_44[k]
                      + mf_44[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, lf_45, lf_46, lf_47, lf_48, \
                         lf_49, mf_45, mf_46, mf_47, mf_48, mf_49 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = ab_x[k] * lf_45[k]
                      + mf_45[k];

            t_66[k] = ab_x[k] * lf_46[k]
                      + mf_46[k];

            t_67[k] = ab_x[k] * lf_47[k]
                      + mf_47[k];

            t_68[k] = ab_x[k] * lf_48[k]
                      + mf_48[k];

            t_69[k] = ab_x[k] * lf_49[k]
                      + mf_49[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_y, ab_z, lf_46, lf_47, lf_48, lf_49, \
                         mf_76, mf_77, mf_78, mf_79, mf_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = ab_y[k] * lf_46[k]
                      + mf_76[k];

            t_71[k] = ab_y[k] * lf_47[k]
                      + mf_77[k];

            t_72[k] = ab_y[k] * lf_48[k]
                      + mf_78[k];

            t_73[k] = ab_y[k] * lf_49[k]
                      + mf_79[k];

            t_74[k] = ab_z[k] * lf_49[k]
                      + mf_89[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, lf_50, lf_51, lf_52, lf_53, \
                         lf_54, mf_50, mf_51, mf_52, mf_53, mf_54 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = ab_x[k] * lf_50[k]
                      + mf_50[k];

            t_76[k] = ab_x[k] * lf_51[k]
                      + mf_51[k];

            t_77[k] = ab_x[k] * lf_52[k]
                      + mf_52[k];

            t_78[k] = ab_x[k] * lf_53[k]
                      + mf_53[k];

            t_79[k] = ab_x[k] * lf_54[k]
                      + mf_54[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, lf_55, lf_56, lf_57, lf_58, \
                         lf_59, mf_55, mf_56, mf_57, mf_58, mf_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_x[k] * lf_55[k]
                      + mf_55[k];

            t_81[k] = ab_x[k] * lf_56[k]
                      + mf_56[k];

            t_82[k] = ab_x[k] * lf_57[k]
                      + mf_57[k];

            t_83[k] = ab_x[k] * lf_58[k]
                      + mf_58[k];

            t_84[k] = ab_x[k] * lf_59[k]
                      + mf_59[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_y, ab_z, lf_56, lf_57, lf_58, lf_59, \
                         mf_86, mf_87, mf_88, mf_89, mf_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = ab_y[k] * lf_56[k]
                      + mf_86[k];

            t_86[k] = ab_y[k] * lf_57[k]
                      + mf_87[k];

            t_87[k] = ab_y[k] * lf_58[k]
                      + mf_88[k];

            t_88[k] = ab_y[k] * lf_59[k]
                      + mf_89[k];

            t_89[k] = ab_z[k] * lf_59[k]
                      + mf_99[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, lf_60, lf_61, lf_62, lf_63, \
                         lf_64, mf_60, mf_61, mf_62, mf_63, mf_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = ab_x[k] * lf_60[k]
                      + mf_60[k];

            t_91[k] = ab_x[k] * lf_61[k]
                      + mf_61[k];

            t_92[k] = ab_x[k] * lf_62[k]
                      + mf_62[k];

            t_93[k] = ab_x[k] * lf_63[k]
                      + mf_63[k];

            t_94[k] = ab_x[k] * lf_64[k]
                      + mf_64[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, lf_65, lf_66, lf_67, lf_68, \
                         lf_69, mf_65, mf_66, mf_67, mf_68, mf_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = ab_x[k] * lf_65[k]
                      + mf_65[k];

            t_96[k] = ab_x[k] * lf_66[k]
                      + mf_66[k];

            t_97[k] = ab_x[k] * lf_67[k]
                      + mf_67[k];

            t_98[k] = ab_x[k] * lf_68[k]
                      + mf_68[k];

            t_99[k] = ab_x[k] * lf_69[k]
                      + mf_69[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_y, ab_z, lf_66, lf_67, lf_68, \
                         lf_69, mf_106, mf_107, mf_108, mf_109, \
                         mf_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = ab_y[k] * lf_66[k]
                       + mf_106[k];

            t_101[k] = ab_y[k] * lf_67[k]
                       + mf_107[k];

            t_102[k] = ab_y[k] * lf_68[k]
                       + mf_108[k];

            t_103[k] = ab_y[k] * lf_69[k]
                       + mf_109[k];

            t_104[k] = ab_z[k] * lf_69[k]
                       + mf_119[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, lf_70, lf_71, lf_72, lf_73, \
                         lf_74, mf_70, mf_71, mf_72, mf_73, mf_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = ab_x[k] * lf_70[k]
                       + mf_70[k];

            t_106[k] = ab_x[k] * lf_71[k]
                       + mf_71[k];

            t_107[k] = ab_x[k] * lf_72[k]
                       + mf_72[k];

            t_108[k] = ab_x[k] * lf_73[k]
                       + mf_73[k];

            t_109[k] = ab_x[k] * lf_74[k]
                       + mf_74[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, lf_75, lf_76, lf_77, lf_78, \
                         lf_79, mf_75, mf_76, mf_77, mf_78, mf_79 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = ab_x[k] * lf_75[k]
                       + mf_75[k];

            t_111[k] = ab_x[k] * lf_76[k]
                       + mf_76[k];

            t_112[k] = ab_x[k] * lf_77[k]
                       + mf_77[k];

            t_113[k] = ab_x[k] * lf_78[k]
                       + mf_78[k];

            t_114[k] = ab_x[k] * lf_79[k]
                       + mf_79[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_y, ab_z, lf_76, lf_77, lf_78, \
                         lf_79, mf_116, mf_117, mf_118, mf_119, \
                         mf_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = ab_y[k] * lf_76[k]
                       + mf_116[k];

            t_116[k] = ab_y[k] * lf_77[k]
                       + mf_117[k];

            t_117[k] = ab_y[k] * lf_78[k]
                       + mf_118[k];

            t_118[k] = ab_y[k] * lf_79[k]
                       + mf_119[k];

            t_119[k] = ab_z[k] * lf_79[k]
                       + mf_129[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, lf_80, lf_81, lf_82, lf_83, \
                         lf_84, mf_80, mf_81, mf_82, mf_83, mf_84 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = ab_x[k] * lf_80[k]
                       + mf_80[k];

            t_121[k] = ab_x[k] * lf_81[k]
                       + mf_81[k];

            t_122[k] = ab_x[k] * lf_82[k]
                       + mf_82[k];

            t_123[k] = ab_x[k] * lf_83[k]
                       + mf_83[k];

            t_124[k] = ab_x[k] * lf_84[k]
                       + mf_84[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, lf_85, lf_86, lf_87, lf_88, \
                         lf_89, mf_85, mf_86, mf_87, mf_88, mf_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = ab_x[k] * lf_85[k]
                       + mf_85[k];

            t_126[k] = ab_x[k] * lf_86[k]
                       + mf_86[k];

            t_127[k] = ab_x[k] * lf_87[k]
                       + mf_87[k];

            t_128[k] = ab_x[k] * lf_88[k]
                       + mf_88[k];

            t_129[k] = ab_x[k] * lf_89[k]
                       + mf_89[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_y, ab_z, lf_86, lf_87, lf_88, \
                         lf_89, mf_126, mf_127, mf_128, mf_129, \
                         mf_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = ab_y[k] * lf_86[k]
                       + mf_126[k];

            t_131[k] = ab_y[k] * lf_87[k]
                       + mf_127[k];

            t_132[k] = ab_y[k] * lf_88[k]
                       + mf_128[k];

            t_133[k] = ab_y[k] * lf_89[k]
                       + mf_129[k];

            t_134[k] = ab_z[k] * lf_89[k]
                       + mf_139[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, lf_90, lf_91, lf_92, lf_93, \
                         lf_94, mf_90, mf_91, mf_92, mf_93, mf_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = ab_x[k] * lf_90[k]
                       + mf_90[k];

            t_136[k] = ab_x[k] * lf_91[k]
                       + mf_91[k];

            t_137[k] = ab_x[k] * lf_92[k]
                       + mf_92[k];

            t_138[k] = ab_x[k] * lf_93[k]
                       + mf_93[k];

            t_139[k] = ab_x[k] * lf_94[k]
                       + mf_94[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, lf_95, lf_96, lf_97, lf_98, \
                         lf_99, mf_95, mf_96, mf_97, mf_98, mf_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = ab_x[k] * lf_95[k]
                       + mf_95[k];

            t_141[k] = ab_x[k] * lf_96[k]
                       + mf_96[k];

            t_142[k] = ab_x[k] * lf_97[k]
                       + mf_97[k];

            t_143[k] = ab_x[k] * lf_98[k]
                       + mf_98[k];

            t_144[k] = ab_x[k] * lf_99[k]
                       + mf_99[k];
        }
    }
}

static auto
compute_hrr_lg_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t lf, const size_t mf,
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

        const auto *lf_96 = buffer.data(lf + 96 * ncomps + c);
        const auto *lf_97 = buffer.data(lf + 97 * ncomps + c);
        const auto *lf_98 = buffer.data(lf + 98 * ncomps + c);
        const auto *lf_99 = buffer.data(lf + 99 * ncomps + c);
        const auto *lf_100 = buffer.data(lf + 100 * ncomps + c);
        const auto *lf_101 = buffer.data(lf + 101 * ncomps + c);
        const auto *lf_102 = buffer.data(lf + 102 * ncomps + c);
        const auto *lf_103 = buffer.data(lf + 103 * ncomps + c);
        const auto *lf_104 = buffer.data(lf + 104 * ncomps + c);
        const auto *lf_105 = buffer.data(lf + 105 * ncomps + c);
        const auto *lf_106 = buffer.data(lf + 106 * ncomps + c);
        const auto *lf_107 = buffer.data(lf + 107 * ncomps + c);
        const auto *lf_108 = buffer.data(lf + 108 * ncomps + c);
        const auto *lf_109 = buffer.data(lf + 109 * ncomps + c);
        const auto *lf_110 = buffer.data(lf + 110 * ncomps + c);
        const auto *lf_111 = buffer.data(lf + 111 * ncomps + c);
        const auto *lf_112 = buffer.data(lf + 112 * ncomps + c);
        const auto *lf_113 = buffer.data(lf + 113 * ncomps + c);
        const auto *lf_114 = buffer.data(lf + 114 * ncomps + c);
        const auto *lf_115 = buffer.data(lf + 115 * ncomps + c);
        const auto *lf_116 = buffer.data(lf + 116 * ncomps + c);
        const auto *lf_117 = buffer.data(lf + 117 * ncomps + c);
        const auto *lf_118 = buffer.data(lf + 118 * ncomps + c);
        const auto *lf_119 = buffer.data(lf + 119 * ncomps + c);
        const auto *lf_120 = buffer.data(lf + 120 * ncomps + c);
        const auto *lf_121 = buffer.data(lf + 121 * ncomps + c);
        const auto *lf_122 = buffer.data(lf + 122 * ncomps + c);
        const auto *lf_123 = buffer.data(lf + 123 * ncomps + c);
        const auto *lf_124 = buffer.data(lf + 124 * ncomps + c);
        const auto *lf_125 = buffer.data(lf + 125 * ncomps + c);
        const auto *lf_126 = buffer.data(lf + 126 * ncomps + c);
        const auto *lf_127 = buffer.data(lf + 127 * ncomps + c);
        const auto *lf_128 = buffer.data(lf + 128 * ncomps + c);
        const auto *lf_129 = buffer.data(lf + 129 * ncomps + c);
        const auto *lf_130 = buffer.data(lf + 130 * ncomps + c);
        const auto *lf_131 = buffer.data(lf + 131 * ncomps + c);
        const auto *lf_132 = buffer.data(lf + 132 * ncomps + c);
        const auto *lf_133 = buffer.data(lf + 133 * ncomps + c);
        const auto *lf_134 = buffer.data(lf + 134 * ncomps + c);
        const auto *lf_135 = buffer.data(lf + 135 * ncomps + c);
        const auto *lf_136 = buffer.data(lf + 136 * ncomps + c);
        const auto *lf_137 = buffer.data(lf + 137 * ncomps + c);
        const auto *lf_138 = buffer.data(lf + 138 * ncomps + c);
        const auto *lf_139 = buffer.data(lf + 139 * ncomps + c);
        const auto *lf_140 = buffer.data(lf + 140 * ncomps + c);
        const auto *lf_141 = buffer.data(lf + 141 * ncomps + c);
        const auto *lf_142 = buffer.data(lf + 142 * ncomps + c);
        const auto *lf_143 = buffer.data(lf + 143 * ncomps + c);
        const auto *lf_144 = buffer.data(lf + 144 * ncomps + c);
        const auto *lf_145 = buffer.data(lf + 145 * ncomps + c);
        const auto *lf_146 = buffer.data(lf + 146 * ncomps + c);
        const auto *lf_147 = buffer.data(lf + 147 * ncomps + c);
        const auto *lf_148 = buffer.data(lf + 148 * ncomps + c);
        const auto *lf_149 = buffer.data(lf + 149 * ncomps + c);
        const auto *lf_150 = buffer.data(lf + 150 * ncomps + c);
        const auto *lf_151 = buffer.data(lf + 151 * ncomps + c);
        const auto *lf_152 = buffer.data(lf + 152 * ncomps + c);
        const auto *lf_153 = buffer.data(lf + 153 * ncomps + c);
        const auto *lf_154 = buffer.data(lf + 154 * ncomps + c);
        const auto *lf_155 = buffer.data(lf + 155 * ncomps + c);
        const auto *lf_156 = buffer.data(lf + 156 * ncomps + c);
        const auto *lf_157 = buffer.data(lf + 157 * ncomps + c);
        const auto *lf_158 = buffer.data(lf + 158 * ncomps + c);
        const auto *lf_159 = buffer.data(lf + 159 * ncomps + c);
        const auto *lf_160 = buffer.data(lf + 160 * ncomps + c);
        const auto *lf_161 = buffer.data(lf + 161 * ncomps + c);
        const auto *lf_162 = buffer.data(lf + 162 * ncomps + c);
        const auto *lf_163 = buffer.data(lf + 163 * ncomps + c);
        const auto *lf_164 = buffer.data(lf + 164 * ncomps + c);
        const auto *lf_165 = buffer.data(lf + 165 * ncomps + c);
        const auto *lf_166 = buffer.data(lf + 166 * ncomps + c);
        const auto *lf_167 = buffer.data(lf + 167 * ncomps + c);
        const auto *lf_168 = buffer.data(lf + 168 * ncomps + c);
        const auto *lf_169 = buffer.data(lf + 169 * ncomps + c);
        const auto *lf_170 = buffer.data(lf + 170 * ncomps + c);
        const auto *lf_171 = buffer.data(lf + 171 * ncomps + c);
        const auto *lf_172 = buffer.data(lf + 172 * ncomps + c);
        const auto *lf_173 = buffer.data(lf + 173 * ncomps + c);
        const auto *lf_174 = buffer.data(lf + 174 * ncomps + c);
        const auto *lf_175 = buffer.data(lf + 175 * ncomps + c);
        const auto *lf_176 = buffer.data(lf + 176 * ncomps + c);
        const auto *lf_177 = buffer.data(lf + 177 * ncomps + c);
        const auto *lf_178 = buffer.data(lf + 178 * ncomps + c);
        const auto *lf_179 = buffer.data(lf + 179 * ncomps + c);
        const auto *lf_180 = buffer.data(lf + 180 * ncomps + c);
        const auto *lf_181 = buffer.data(lf + 181 * ncomps + c);
        const auto *lf_182 = buffer.data(lf + 182 * ncomps + c);
        const auto *lf_183 = buffer.data(lf + 183 * ncomps + c);
        const auto *lf_184 = buffer.data(lf + 184 * ncomps + c);
        const auto *lf_185 = buffer.data(lf + 185 * ncomps + c);
        const auto *lf_186 = buffer.data(lf + 186 * ncomps + c);
        const auto *lf_187 = buffer.data(lf + 187 * ncomps + c);
        const auto *lf_188 = buffer.data(lf + 188 * ncomps + c);
        const auto *lf_189 = buffer.data(lf + 189 * ncomps + c);
        const auto *lf_190 = buffer.data(lf + 190 * ncomps + c);
        const auto *lf_191 = buffer.data(lf + 191 * ncomps + c);
        const auto *lf_192 = buffer.data(lf + 192 * ncomps + c);
        const auto *lf_193 = buffer.data(lf + 193 * ncomps + c);
        const auto *lf_194 = buffer.data(lf + 194 * ncomps + c);

        const auto *mf_100 = buffer.data(mf + 100 * ncomps + c);
        const auto *mf_101 = buffer.data(mf + 101 * ncomps + c);
        const auto *mf_102 = buffer.data(mf + 102 * ncomps + c);
        const auto *mf_103 = buffer.data(mf + 103 * ncomps + c);
        const auto *mf_104 = buffer.data(mf + 104 * ncomps + c);
        const auto *mf_105 = buffer.data(mf + 105 * ncomps + c);
        const auto *mf_106 = buffer.data(mf + 106 * ncomps + c);
        const auto *mf_107 = buffer.data(mf + 107 * ncomps + c);
        const auto *mf_108 = buffer.data(mf + 108 * ncomps + c);
        const auto *mf_109 = buffer.data(mf + 109 * ncomps + c);
        const auto *mf_110 = buffer.data(mf + 110 * ncomps + c);
        const auto *mf_111 = buffer.data(mf + 111 * ncomps + c);
        const auto *mf_112 = buffer.data(mf + 112 * ncomps + c);
        const auto *mf_113 = buffer.data(mf + 113 * ncomps + c);
        const auto *mf_114 = buffer.data(mf + 114 * ncomps + c);
        const auto *mf_115 = buffer.data(mf + 115 * ncomps + c);
        const auto *mf_116 = buffer.data(mf + 116 * ncomps + c);
        const auto *mf_117 = buffer.data(mf + 117 * ncomps + c);
        const auto *mf_118 = buffer.data(mf + 118 * ncomps + c);
        const auto *mf_119 = buffer.data(mf + 119 * ncomps + c);
        const auto *mf_120 = buffer.data(mf + 120 * ncomps + c);
        const auto *mf_121 = buffer.data(mf + 121 * ncomps + c);
        const auto *mf_122 = buffer.data(mf + 122 * ncomps + c);
        const auto *mf_123 = buffer.data(mf + 123 * ncomps + c);
        const auto *mf_124 = buffer.data(mf + 124 * ncomps + c);
        const auto *mf_125 = buffer.data(mf + 125 * ncomps + c);
        const auto *mf_126 = buffer.data(mf + 126 * ncomps + c);
        const auto *mf_127 = buffer.data(mf + 127 * ncomps + c);
        const auto *mf_128 = buffer.data(mf + 128 * ncomps + c);
        const auto *mf_129 = buffer.data(mf + 129 * ncomps + c);
        const auto *mf_130 = buffer.data(mf + 130 * ncomps + c);
        const auto *mf_131 = buffer.data(mf + 131 * ncomps + c);
        const auto *mf_132 = buffer.data(mf + 132 * ncomps + c);
        const auto *mf_133 = buffer.data(mf + 133 * ncomps + c);
        const auto *mf_134 = buffer.data(mf + 134 * ncomps + c);
        const auto *mf_135 = buffer.data(mf + 135 * ncomps + c);
        const auto *mf_136 = buffer.data(mf + 136 * ncomps + c);
        const auto *mf_137 = buffer.data(mf + 137 * ncomps + c);
        const auto *mf_138 = buffer.data(mf + 138 * ncomps + c);
        const auto *mf_139 = buffer.data(mf + 139 * ncomps + c);
        const auto *mf_140 = buffer.data(mf + 140 * ncomps + c);
        const auto *mf_141 = buffer.data(mf + 141 * ncomps + c);
        const auto *mf_142 = buffer.data(mf + 142 * ncomps + c);
        const auto *mf_143 = buffer.data(mf + 143 * ncomps + c);
        const auto *mf_144 = buffer.data(mf + 144 * ncomps + c);
        const auto *mf_145 = buffer.data(mf + 145 * ncomps + c);
        const auto *mf_146 = buffer.data(mf + 146 * ncomps + c);
        const auto *mf_147 = buffer.data(mf + 147 * ncomps + c);
        const auto *mf_148 = buffer.data(mf + 148 * ncomps + c);
        const auto *mf_149 = buffer.data(mf + 149 * ncomps + c);
        const auto *mf_150 = buffer.data(mf + 150 * ncomps + c);
        const auto *mf_151 = buffer.data(mf + 151 * ncomps + c);
        const auto *mf_152 = buffer.data(mf + 152 * ncomps + c);
        const auto *mf_153 = buffer.data(mf + 153 * ncomps + c);
        const auto *mf_154 = buffer.data(mf + 154 * ncomps + c);
        const auto *mf_155 = buffer.data(mf + 155 * ncomps + c);
        const auto *mf_156 = buffer.data(mf + 156 * ncomps + c);
        const auto *mf_157 = buffer.data(mf + 157 * ncomps + c);
        const auto *mf_158 = buffer.data(mf + 158 * ncomps + c);
        const auto *mf_159 = buffer.data(mf + 159 * ncomps + c);
        const auto *mf_160 = buffer.data(mf + 160 * ncomps + c);
        const auto *mf_161 = buffer.data(mf + 161 * ncomps + c);
        const auto *mf_162 = buffer.data(mf + 162 * ncomps + c);
        const auto *mf_163 = buffer.data(mf + 163 * ncomps + c);
        const auto *mf_164 = buffer.data(mf + 164 * ncomps + c);
        const auto *mf_165 = buffer.data(mf + 165 * ncomps + c);
        const auto *mf_166 = buffer.data(mf + 166 * ncomps + c);
        const auto *mf_167 = buffer.data(mf + 167 * ncomps + c);
        const auto *mf_168 = buffer.data(mf + 168 * ncomps + c);
        const auto *mf_169 = buffer.data(mf + 169 * ncomps + c);
        const auto *mf_170 = buffer.data(mf + 170 * ncomps + c);
        const auto *mf_171 = buffer.data(mf + 171 * ncomps + c);
        const auto *mf_172 = buffer.data(mf + 172 * ncomps + c);
        const auto *mf_173 = buffer.data(mf + 173 * ncomps + c);
        const auto *mf_174 = buffer.data(mf + 174 * ncomps + c);
        const auto *mf_175 = buffer.data(mf + 175 * ncomps + c);
        const auto *mf_176 = buffer.data(mf + 176 * ncomps + c);
        const auto *mf_177 = buffer.data(mf + 177 * ncomps + c);
        const auto *mf_178 = buffer.data(mf + 178 * ncomps + c);
        const auto *mf_179 = buffer.data(mf + 179 * ncomps + c);
        const auto *mf_180 = buffer.data(mf + 180 * ncomps + c);
        const auto *mf_181 = buffer.data(mf + 181 * ncomps + c);
        const auto *mf_182 = buffer.data(mf + 182 * ncomps + c);
        const auto *mf_183 = buffer.data(mf + 183 * ncomps + c);
        const auto *mf_184 = buffer.data(mf + 184 * ncomps + c);
        const auto *mf_185 = buffer.data(mf + 185 * ncomps + c);
        const auto *mf_186 = buffer.data(mf + 186 * ncomps + c);
        const auto *mf_187 = buffer.data(mf + 187 * ncomps + c);
        const auto *mf_188 = buffer.data(mf + 188 * ncomps + c);
        const auto *mf_189 = buffer.data(mf + 189 * ncomps + c);
        const auto *mf_190 = buffer.data(mf + 190 * ncomps + c);
        const auto *mf_191 = buffer.data(mf + 191 * ncomps + c);
        const auto *mf_192 = buffer.data(mf + 192 * ncomps + c);
        const auto *mf_193 = buffer.data(mf + 193 * ncomps + c);
        const auto *mf_194 = buffer.data(mf + 194 * ncomps + c);
        const auto *mf_196 = buffer.data(mf + 196 * ncomps + c);
        const auto *mf_197 = buffer.data(mf + 197 * ncomps + c);
        const auto *mf_198 = buffer.data(mf + 198 * ncomps + c);
        const auto *mf_199 = buffer.data(mf + 199 * ncomps + c);
        const auto *mf_209 = buffer.data(mf + 209 * ncomps + c);
        const auto *mf_216 = buffer.data(mf + 216 * ncomps + c);
        const auto *mf_217 = buffer.data(mf + 217 * ncomps + c);
        const auto *mf_218 = buffer.data(mf + 218 * ncomps + c);
        const auto *mf_219 = buffer.data(mf + 219 * ncomps + c);
        const auto *mf_226 = buffer.data(mf + 226 * ncomps + c);
        const auto *mf_227 = buffer.data(mf + 227 * ncomps + c);
        const auto *mf_228 = buffer.data(mf + 228 * ncomps + c);
        const auto *mf_229 = buffer.data(mf + 229 * ncomps + c);
        const auto *mf_236 = buffer.data(mf + 236 * ncomps + c);
        const auto *mf_237 = buffer.data(mf + 237 * ncomps + c);
        const auto *mf_238 = buffer.data(mf + 238 * ncomps + c);
        const auto *mf_239 = buffer.data(mf + 239 * ncomps + c);
        const auto *mf_246 = buffer.data(mf + 246 * ncomps + c);
        const auto *mf_247 = buffer.data(mf + 247 * ncomps + c);
        const auto *mf_248 = buffer.data(mf + 248 * ncomps + c);
        const auto *mf_249 = buffer.data(mf + 249 * ncomps + c);
        const auto *mf_259 = buffer.data(mf + 259 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_y, ab_z, lf_96, lf_97, lf_98, \
                         lf_99, mf_136, mf_137, mf_138, mf_139, \
                         mf_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = ab_y[k] * lf_96[k]
                       + mf_136[k];

            t_146[k] = ab_y[k] * lf_97[k]
                       + mf_137[k];

            t_147[k] = ab_y[k] * lf_98[k]
                       + mf_138[k];

            t_148[k] = ab_y[k] * lf_99[k]
                       + mf_139[k];

            t_149[k] = ab_z[k] * lf_99[k]
                       + mf_149[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, lf_100, lf_101, lf_102, \
                         lf_103, lf_104, mf_100, mf_101, mf_102, mf_103, \
                         mf_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * lf_100[k]
                       + mf_100[k];

            t_151[k] = ab_x[k] * lf_101[k]
                       + mf_101[k];

            t_152[k] = ab_x[k] * lf_102[k]
                       + mf_102[k];

            t_153[k] = ab_x[k] * lf_103[k]
                       + mf_103[k];

            t_154[k] = ab_x[k] * lf_104[k]
                       + mf_104[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, lf_105, lf_106, lf_107, \
                         lf_108, lf_109, mf_105, mf_106, mf_107, mf_108, \
                         mf_109 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * lf_105[k]
                       + mf_105[k];

            t_156[k] = ab_x[k] * lf_106[k]
                       + mf_106[k];

            t_157[k] = ab_x[k] * lf_107[k]
                       + mf_107[k];

            t_158[k] = ab_x[k] * lf_108[k]
                       + mf_108[k];

            t_159[k] = ab_x[k] * lf_109[k]
                       + mf_109[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_y, ab_z, lf_106, lf_107, \
                         lf_108, lf_109, mf_156, mf_157, mf_158, mf_159, \
                         mf_169 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_y[k] * lf_106[k]
                       + mf_156[k];

            t_161[k] = ab_y[k] * lf_107[k]
                       + mf_157[k];

            t_162[k] = ab_y[k] * lf_108[k]
                       + mf_158[k];

            t_163[k] = ab_y[k] * lf_109[k]
                       + mf_159[k];

            t_164[k] = ab_z[k] * lf_109[k]
                       + mf_169[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, lf_110, lf_111, lf_112, \
                         lf_113, lf_114, mf_110, mf_111, mf_112, mf_113, \
                         mf_114 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = ab_x[k] * lf_110[k]
                       + mf_110[k];

            t_166[k] = ab_x[k] * lf_111[k]
                       + mf_111[k];

            t_167[k] = ab_x[k] * lf_112[k]
                       + mf_112[k];

            t_168[k] = ab_x[k] * lf_113[k]
                       + mf_113[k];

            t_169[k] = ab_x[k] * lf_114[k]
                       + mf_114[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, lf_115, lf_116, lf_117, \
                         lf_118, lf_119, mf_115, mf_116, mf_117, mf_118, \
                         mf_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = ab_x[k] * lf_115[k]
                       + mf_115[k];

            t_171[k] = ab_x[k] * lf_116[k]
                       + mf_116[k];

            t_172[k] = ab_x[k] * lf_117[k]
                       + mf_117[k];

            t_173[k] = ab_x[k] * lf_118[k]
                       + mf_118[k];

            t_174[k] = ab_x[k] * lf_119[k]
                       + mf_119[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_y, ab_z, lf_116, lf_117, \
                         lf_118, lf_119, mf_166, mf_167, mf_168, mf_169, \
                         mf_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = ab_y[k] * lf_116[k]
                       + mf_166[k];

            t_176[k] = ab_y[k] * lf_117[k]
                       + mf_167[k];

            t_177[k] = ab_y[k] * lf_118[k]
                       + mf_168[k];

            t_178[k] = ab_y[k] * lf_119[k]
                       + mf_169[k];

            t_179[k] = ab_z[k] * lf_119[k]
                       + mf_179[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, lf_120, lf_121, lf_122, \
                         lf_123, lf_124, mf_120, mf_121, mf_122, mf_123, \
                         mf_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = ab_x[k] * lf_120[k]
                       + mf_120[k];

            t_181[k] = ab_x[k] * lf_121[k]
                       + mf_121[k];

            t_182[k] = ab_x[k] * lf_122[k]
                       + mf_122[k];

            t_183[k] = ab_x[k] * lf_123[k]
                       + mf_123[k];

            t_184[k] = ab_x[k] * lf_124[k]
                       + mf_124[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, lf_125, lf_126, lf_127, \
                         lf_128, lf_129, mf_125, mf_126, mf_127, mf_128, \
                         mf_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = ab_x[k] * lf_125[k]
                       + mf_125[k];

            t_186[k] = ab_x[k] * lf_126[k]
                       + mf_126[k];

            t_187[k] = ab_x[k] * lf_127[k]
                       + mf_127[k];

            t_188[k] = ab_x[k] * lf_128[k]
                       + mf_128[k];

            t_189[k] = ab_x[k] * lf_129[k]
                       + mf_129[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_y, ab_z, lf_126, lf_127, \
                         lf_128, lf_129, mf_176, mf_177, mf_178, mf_179, \
                         mf_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = ab_y[k] * lf_126[k]
                       + mf_176[k];

            t_191[k] = ab_y[k] * lf_127[k]
                       + mf_177[k];

            t_192[k] = ab_y[k] * lf_128[k]
                       + mf_178[k];

            t_193[k] = ab_y[k] * lf_129[k]
                       + mf_179[k];

            t_194[k] = ab_z[k] * lf_129[k]
                       + mf_189[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, lf_130, lf_131, lf_132, \
                         lf_133, lf_134, mf_130, mf_131, mf_132, mf_133, \
                         mf_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = ab_x[k] * lf_130[k]
                       + mf_130[k];

            t_196[k] = ab_x[k] * lf_131[k]
                       + mf_131[k];

            t_197[k] = ab_x[k] * lf_132[k]
                       + mf_132[k];

            t_198[k] = ab_x[k] * lf_133[k]
                       + mf_133[k];

            t_199[k] = ab_x[k] * lf_134[k]
                       + mf_134[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, lf_135, lf_136, lf_137, \
                         lf_138, lf_139, mf_135, mf_136, mf_137, mf_138, \
                         mf_139 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = ab_x[k] * lf_135[k]
                       + mf_135[k];

            t_201[k] = ab_x[k] * lf_136[k]
                       + mf_136[k];

            t_202[k] = ab_x[k] * lf_137[k]
                       + mf_137[k];

            t_203[k] = ab_x[k] * lf_138[k]
                       + mf_138[k];

            t_204[k] = ab_x[k] * lf_139[k]
                       + mf_139[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_y, ab_z, lf_136, lf_137, \
                         lf_138, lf_139, mf_186, mf_187, mf_188, mf_189, \
                         mf_199 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = ab_y[k] * lf_136[k]
                       + mf_186[k];

            t_206[k] = ab_y[k] * lf_137[k]
                       + mf_187[k];

            t_207[k] = ab_y[k] * lf_138[k]
                       + mf_188[k];

            t_208[k] = ab_y[k] * lf_139[k]
                       + mf_189[k];

            t_209[k] = ab_z[k] * lf_139[k]
                       + mf_199[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, lf_140, lf_141, lf_142, \
                         lf_143, lf_144, mf_140, mf_141, mf_142, mf_143, \
                         mf_144 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = ab_x[k] * lf_140[k]
                       + mf_140[k];

            t_211[k] = ab_x[k] * lf_141[k]
                       + mf_141[k];

            t_212[k] = ab_x[k] * lf_142[k]
                       + mf_142[k];

            t_213[k] = ab_x[k] * lf_143[k]
                       + mf_143[k];

            t_214[k] = ab_x[k] * lf_144[k]
                       + mf_144[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, lf_145, lf_146, lf_147, \
                         lf_148, lf_149, mf_145, mf_146, mf_147, mf_148, \
                         mf_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = ab_x[k] * lf_145[k]
                       + mf_145[k];

            t_216[k] = ab_x[k] * lf_146[k]
                       + mf_146[k];

            t_217[k] = ab_x[k] * lf_147[k]
                       + mf_147[k];

            t_218[k] = ab_x[k] * lf_148[k]
                       + mf_148[k];

            t_219[k] = ab_x[k] * lf_149[k]
                       + mf_149[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_y, ab_z, lf_146, lf_147, \
                         lf_148, lf_149, mf_196, mf_197, mf_198, mf_199, \
                         mf_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = ab_y[k] * lf_146[k]
                       + mf_196[k];

            t_221[k] = ab_y[k] * lf_147[k]
                       + mf_197[k];

            t_222[k] = ab_y[k] * lf_148[k]
                       + mf_198[k];

            t_223[k] = ab_y[k] * lf_149[k]
                       + mf_199[k];

            t_224[k] = ab_z[k] * lf_149[k]
                       + mf_209[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, lf_150, lf_151, lf_152, \
                         lf_153, lf_154, mf_150, mf_151, mf_152, mf_153, \
                         mf_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = ab_x[k] * lf_150[k]
                       + mf_150[k];

            t_226[k] = ab_x[k] * lf_151[k]
                       + mf_151[k];

            t_227[k] = ab_x[k] * lf_152[k]
                       + mf_152[k];

            t_228[k] = ab_x[k] * lf_153[k]
                       + mf_153[k];

            t_229[k] = ab_x[k] * lf_154[k]
                       + mf_154[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, lf_155, lf_156, lf_157, \
                         lf_158, lf_159, mf_155, mf_156, mf_157, mf_158, \
                         mf_159 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = ab_x[k] * lf_155[k]
                       + mf_155[k];

            t_231[k] = ab_x[k] * lf_156[k]
                       + mf_156[k];

            t_232[k] = ab_x[k] * lf_157[k]
                       + mf_157[k];

            t_233[k] = ab_x[k] * lf_158[k]
                       + mf_158[k];

            t_234[k] = ab_x[k] * lf_159[k]
                       + mf_159[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_y, ab_z, lf_156, lf_157, \
                         lf_158, lf_159, mf_216, mf_217, mf_218, mf_219, \
                         mf_229 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = ab_y[k] * lf_156[k]
                       + mf_216[k];

            t_236[k] = ab_y[k] * lf_157[k]
                       + mf_217[k];

            t_237[k] = ab_y[k] * lf_158[k]
                       + mf_218[k];

            t_238[k] = ab_y[k] * lf_159[k]
                       + mf_219[k];

            t_239[k] = ab_z[k] * lf_159[k]
                       + mf_229[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, lf_160, lf_161, lf_162, \
                         lf_163, lf_164, mf_160, mf_161, mf_162, mf_163, \
                         mf_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = ab_x[k] * lf_160[k]
                       + mf_160[k];

            t_241[k] = ab_x[k] * lf_161[k]
                       + mf_161[k];

            t_242[k] = ab_x[k] * lf_162[k]
                       + mf_162[k];

            t_243[k] = ab_x[k] * lf_163[k]
                       + mf_163[k];

            t_244[k] = ab_x[k] * lf_164[k]
                       + mf_164[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, lf_165, lf_166, lf_167, \
                         lf_168, lf_169, mf_165, mf_166, mf_167, mf_168, \
                         mf_169 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = ab_x[k] * lf_165[k]
                       + mf_165[k];

            t_246[k] = ab_x[k] * lf_166[k]
                       + mf_166[k];

            t_247[k] = ab_x[k] * lf_167[k]
                       + mf_167[k];

            t_248[k] = ab_x[k] * lf_168[k]
                       + mf_168[k];

            t_249[k] = ab_x[k] * lf_169[k]
                       + mf_169[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_y, ab_z, lf_166, lf_167, \
                         lf_168, lf_169, mf_226, mf_227, mf_228, mf_229, \
                         mf_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = ab_y[k] * lf_166[k]
                       + mf_226[k];

            t_251[k] = ab_y[k] * lf_167[k]
                       + mf_227[k];

            t_252[k] = ab_y[k] * lf_168[k]
                       + mf_228[k];

            t_253[k] = ab_y[k] * lf_169[k]
                       + mf_229[k];

            t_254[k] = ab_z[k] * lf_169[k]
                       + mf_239[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, lf_170, lf_171, lf_172, \
                         lf_173, lf_174, mf_170, mf_171, mf_172, mf_173, \
                         mf_174 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = ab_x[k] * lf_170[k]
                       + mf_170[k];

            t_256[k] = ab_x[k] * lf_171[k]
                       + mf_171[k];

            t_257[k] = ab_x[k] * lf_172[k]
                       + mf_172[k];

            t_258[k] = ab_x[k] * lf_173[k]
                       + mf_173[k];

            t_259[k] = ab_x[k] * lf_174[k]
                       + mf_174[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, lf_175, lf_176, lf_177, \
                         lf_178, lf_179, mf_175, mf_176, mf_177, mf_178, \
                         mf_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = ab_x[k] * lf_175[k]
                       + mf_175[k];

            t_261[k] = ab_x[k] * lf_176[k]
                       + mf_176[k];

            t_262[k] = ab_x[k] * lf_177[k]
                       + mf_177[k];

            t_263[k] = ab_x[k] * lf_178[k]
                       + mf_178[k];

            t_264[k] = ab_x[k] * lf_179[k]
                       + mf_179[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_y, ab_z, lf_176, lf_177, \
                         lf_178, lf_179, mf_236, mf_237, mf_238, mf_239, \
                         mf_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = ab_y[k] * lf_176[k]
                       + mf_236[k];

            t_266[k] = ab_y[k] * lf_177[k]
                       + mf_237[k];

            t_267[k] = ab_y[k] * lf_178[k]
                       + mf_238[k];

            t_268[k] = ab_y[k] * lf_179[k]
                       + mf_239[k];

            t_269[k] = ab_z[k] * lf_179[k]
                       + mf_249[k];
        }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, lf_180, lf_181, lf_182, \
                         lf_183, lf_184, mf_180, mf_181, mf_182, mf_183, \
                         mf_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_270[k] = ab_x[k] * lf_180[k]
                       + mf_180[k];

            t_271[k] = ab_x[k] * lf_181[k]
                       + mf_181[k];

            t_272[k] = ab_x[k] * lf_182[k]
                       + mf_182[k];

            t_273[k] = ab_x[k] * lf_183[k]
                       + mf_183[k];

            t_274[k] = ab_x[k] * lf_184[k]
                       + mf_184[k];
        }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, lf_185, lf_186, lf_187, \
                         lf_188, lf_189, mf_185, mf_186, mf_187, mf_188, \
                         mf_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_275[k] = ab_x[k] * lf_185[k]
                       + mf_185[k];

            t_276[k] = ab_x[k] * lf_186[k]
                       + mf_186[k];

            t_277[k] = ab_x[k] * lf_187[k]
                       + mf_187[k];

            t_278[k] = ab_x[k] * lf_188[k]
                       + mf_188[k];

            t_279[k] = ab_x[k] * lf_189[k]
                       + mf_189[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_y, ab_z, lf_186, lf_187, \
                         lf_188, lf_189, mf_246, mf_247, mf_248, mf_249, \
                         mf_259 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = ab_y[k] * lf_186[k]
                       + mf_246[k];

            t_281[k] = ab_y[k] * lf_187[k]
                       + mf_247[k];

            t_282[k] = ab_y[k] * lf_188[k]
                       + mf_248[k];

            t_283[k] = ab_y[k] * lf_189[k]
                       + mf_249[k];

            t_284[k] = ab_z[k] * lf_189[k]
                       + mf_259[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, lf_190, lf_191, lf_192, \
                         lf_193, lf_194, mf_190, mf_191, mf_192, mf_193, \
                         mf_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = ab_x[k] * lf_190[k]
                       + mf_190[k];

            t_286[k] = ab_x[k] * lf_191[k]
                       + mf_191[k];

            t_287[k] = ab_x[k] * lf_192[k]
                       + mf_192[k];

            t_288[k] = ab_x[k] * lf_193[k]
                       + mf_193[k];

            t_289[k] = ab_x[k] * lf_194[k]
                       + mf_194[k];
        }
    }
}

static auto
compute_hrr_lg_out_of_first_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t lf, const size_t mf,
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

        const auto *lf_195 = buffer.data(lf + 195 * ncomps + c);
        const auto *lf_196 = buffer.data(lf + 196 * ncomps + c);
        const auto *lf_197 = buffer.data(lf + 197 * ncomps + c);
        const auto *lf_198 = buffer.data(lf + 198 * ncomps + c);
        const auto *lf_199 = buffer.data(lf + 199 * ncomps + c);
        const auto *lf_200 = buffer.data(lf + 200 * ncomps + c);
        const auto *lf_201 = buffer.data(lf + 201 * ncomps + c);
        const auto *lf_202 = buffer.data(lf + 202 * ncomps + c);
        const auto *lf_203 = buffer.data(lf + 203 * ncomps + c);
        const auto *lf_204 = buffer.data(lf + 204 * ncomps + c);
        const auto *lf_205 = buffer.data(lf + 205 * ncomps + c);
        const auto *lf_206 = buffer.data(lf + 206 * ncomps + c);
        const auto *lf_207 = buffer.data(lf + 207 * ncomps + c);
        const auto *lf_208 = buffer.data(lf + 208 * ncomps + c);
        const auto *lf_209 = buffer.data(lf + 209 * ncomps + c);
        const auto *lf_210 = buffer.data(lf + 210 * ncomps + c);
        const auto *lf_211 = buffer.data(lf + 211 * ncomps + c);
        const auto *lf_212 = buffer.data(lf + 212 * ncomps + c);
        const auto *lf_213 = buffer.data(lf + 213 * ncomps + c);
        const auto *lf_214 = buffer.data(lf + 214 * ncomps + c);
        const auto *lf_215 = buffer.data(lf + 215 * ncomps + c);
        const auto *lf_216 = buffer.data(lf + 216 * ncomps + c);
        const auto *lf_217 = buffer.data(lf + 217 * ncomps + c);
        const auto *lf_218 = buffer.data(lf + 218 * ncomps + c);
        const auto *lf_219 = buffer.data(lf + 219 * ncomps + c);
        const auto *lf_220 = buffer.data(lf + 220 * ncomps + c);
        const auto *lf_221 = buffer.data(lf + 221 * ncomps + c);
        const auto *lf_222 = buffer.data(lf + 222 * ncomps + c);
        const auto *lf_223 = buffer.data(lf + 223 * ncomps + c);
        const auto *lf_224 = buffer.data(lf + 224 * ncomps + c);
        const auto *lf_225 = buffer.data(lf + 225 * ncomps + c);
        const auto *lf_226 = buffer.data(lf + 226 * ncomps + c);
        const auto *lf_227 = buffer.data(lf + 227 * ncomps + c);
        const auto *lf_228 = buffer.data(lf + 228 * ncomps + c);
        const auto *lf_229 = buffer.data(lf + 229 * ncomps + c);
        const auto *lf_230 = buffer.data(lf + 230 * ncomps + c);
        const auto *lf_231 = buffer.data(lf + 231 * ncomps + c);
        const auto *lf_232 = buffer.data(lf + 232 * ncomps + c);
        const auto *lf_233 = buffer.data(lf + 233 * ncomps + c);
        const auto *lf_234 = buffer.data(lf + 234 * ncomps + c);
        const auto *lf_235 = buffer.data(lf + 235 * ncomps + c);
        const auto *lf_236 = buffer.data(lf + 236 * ncomps + c);
        const auto *lf_237 = buffer.data(lf + 237 * ncomps + c);
        const auto *lf_238 = buffer.data(lf + 238 * ncomps + c);
        const auto *lf_239 = buffer.data(lf + 239 * ncomps + c);
        const auto *lf_240 = buffer.data(lf + 240 * ncomps + c);
        const auto *lf_241 = buffer.data(lf + 241 * ncomps + c);
        const auto *lf_242 = buffer.data(lf + 242 * ncomps + c);
        const auto *lf_243 = buffer.data(lf + 243 * ncomps + c);
        const auto *lf_244 = buffer.data(lf + 244 * ncomps + c);
        const auto *lf_245 = buffer.data(lf + 245 * ncomps + c);
        const auto *lf_246 = buffer.data(lf + 246 * ncomps + c);
        const auto *lf_247 = buffer.data(lf + 247 * ncomps + c);
        const auto *lf_248 = buffer.data(lf + 248 * ncomps + c);
        const auto *lf_249 = buffer.data(lf + 249 * ncomps + c);
        const auto *lf_250 = buffer.data(lf + 250 * ncomps + c);
        const auto *lf_251 = buffer.data(lf + 251 * ncomps + c);
        const auto *lf_252 = buffer.data(lf + 252 * ncomps + c);
        const auto *lf_253 = buffer.data(lf + 253 * ncomps + c);
        const auto *lf_254 = buffer.data(lf + 254 * ncomps + c);
        const auto *lf_255 = buffer.data(lf + 255 * ncomps + c);
        const auto *lf_256 = buffer.data(lf + 256 * ncomps + c);
        const auto *lf_257 = buffer.data(lf + 257 * ncomps + c);
        const auto *lf_258 = buffer.data(lf + 258 * ncomps + c);
        const auto *lf_259 = buffer.data(lf + 259 * ncomps + c);
        const auto *lf_260 = buffer.data(lf + 260 * ncomps + c);
        const auto *lf_261 = buffer.data(lf + 261 * ncomps + c);
        const auto *lf_262 = buffer.data(lf + 262 * ncomps + c);
        const auto *lf_263 = buffer.data(lf + 263 * ncomps + c);
        const auto *lf_264 = buffer.data(lf + 264 * ncomps + c);
        const auto *lf_265 = buffer.data(lf + 265 * ncomps + c);
        const auto *lf_266 = buffer.data(lf + 266 * ncomps + c);
        const auto *lf_267 = buffer.data(lf + 267 * ncomps + c);
        const auto *lf_268 = buffer.data(lf + 268 * ncomps + c);
        const auto *lf_269 = buffer.data(lf + 269 * ncomps + c);
        const auto *lf_270 = buffer.data(lf + 270 * ncomps + c);
        const auto *lf_271 = buffer.data(lf + 271 * ncomps + c);
        const auto *lf_272 = buffer.data(lf + 272 * ncomps + c);
        const auto *lf_273 = buffer.data(lf + 273 * ncomps + c);
        const auto *lf_274 = buffer.data(lf + 274 * ncomps + c);
        const auto *lf_275 = buffer.data(lf + 275 * ncomps + c);
        const auto *lf_276 = buffer.data(lf + 276 * ncomps + c);
        const auto *lf_277 = buffer.data(lf + 277 * ncomps + c);
        const auto *lf_278 = buffer.data(lf + 278 * ncomps + c);
        const auto *lf_279 = buffer.data(lf + 279 * ncomps + c);
        const auto *lf_280 = buffer.data(lf + 280 * ncomps + c);
        const auto *lf_281 = buffer.data(lf + 281 * ncomps + c);
        const auto *lf_282 = buffer.data(lf + 282 * ncomps + c);
        const auto *lf_283 = buffer.data(lf + 283 * ncomps + c);
        const auto *lf_284 = buffer.data(lf + 284 * ncomps + c);
        const auto *lf_285 = buffer.data(lf + 285 * ncomps + c);
        const auto *lf_286 = buffer.data(lf + 286 * ncomps + c);
        const auto *lf_287 = buffer.data(lf + 287 * ncomps + c);
        const auto *lf_288 = buffer.data(lf + 288 * ncomps + c);
        const auto *lf_289 = buffer.data(lf + 289 * ncomps + c);

        const auto *mf_195 = buffer.data(mf + 195 * ncomps + c);
        const auto *mf_196 = buffer.data(mf + 196 * ncomps + c);
        const auto *mf_197 = buffer.data(mf + 197 * ncomps + c);
        const auto *mf_198 = buffer.data(mf + 198 * ncomps + c);
        const auto *mf_199 = buffer.data(mf + 199 * ncomps + c);
        const auto *mf_200 = buffer.data(mf + 200 * ncomps + c);
        const auto *mf_201 = buffer.data(mf + 201 * ncomps + c);
        const auto *mf_202 = buffer.data(mf + 202 * ncomps + c);
        const auto *mf_203 = buffer.data(mf + 203 * ncomps + c);
        const auto *mf_204 = buffer.data(mf + 204 * ncomps + c);
        const auto *mf_205 = buffer.data(mf + 205 * ncomps + c);
        const auto *mf_206 = buffer.data(mf + 206 * ncomps + c);
        const auto *mf_207 = buffer.data(mf + 207 * ncomps + c);
        const auto *mf_208 = buffer.data(mf + 208 * ncomps + c);
        const auto *mf_209 = buffer.data(mf + 209 * ncomps + c);
        const auto *mf_210 = buffer.data(mf + 210 * ncomps + c);
        const auto *mf_211 = buffer.data(mf + 211 * ncomps + c);
        const auto *mf_212 = buffer.data(mf + 212 * ncomps + c);
        const auto *mf_213 = buffer.data(mf + 213 * ncomps + c);
        const auto *mf_214 = buffer.data(mf + 214 * ncomps + c);
        const auto *mf_215 = buffer.data(mf + 215 * ncomps + c);
        const auto *mf_216 = buffer.data(mf + 216 * ncomps + c);
        const auto *mf_217 = buffer.data(mf + 217 * ncomps + c);
        const auto *mf_218 = buffer.data(mf + 218 * ncomps + c);
        const auto *mf_219 = buffer.data(mf + 219 * ncomps + c);
        const auto *mf_220 = buffer.data(mf + 220 * ncomps + c);
        const auto *mf_221 = buffer.data(mf + 221 * ncomps + c);
        const auto *mf_222 = buffer.data(mf + 222 * ncomps + c);
        const auto *mf_223 = buffer.data(mf + 223 * ncomps + c);
        const auto *mf_224 = buffer.data(mf + 224 * ncomps + c);
        const auto *mf_225 = buffer.data(mf + 225 * ncomps + c);
        const auto *mf_226 = buffer.data(mf + 226 * ncomps + c);
        const auto *mf_227 = buffer.data(mf + 227 * ncomps + c);
        const auto *mf_228 = buffer.data(mf + 228 * ncomps + c);
        const auto *mf_229 = buffer.data(mf + 229 * ncomps + c);
        const auto *mf_230 = buffer.data(mf + 230 * ncomps + c);
        const auto *mf_231 = buffer.data(mf + 231 * ncomps + c);
        const auto *mf_232 = buffer.data(mf + 232 * ncomps + c);
        const auto *mf_233 = buffer.data(mf + 233 * ncomps + c);
        const auto *mf_234 = buffer.data(mf + 234 * ncomps + c);
        const auto *mf_235 = buffer.data(mf + 235 * ncomps + c);
        const auto *mf_236 = buffer.data(mf + 236 * ncomps + c);
        const auto *mf_237 = buffer.data(mf + 237 * ncomps + c);
        const auto *mf_238 = buffer.data(mf + 238 * ncomps + c);
        const auto *mf_239 = buffer.data(mf + 239 * ncomps + c);
        const auto *mf_240 = buffer.data(mf + 240 * ncomps + c);
        const auto *mf_241 = buffer.data(mf + 241 * ncomps + c);
        const auto *mf_242 = buffer.data(mf + 242 * ncomps + c);
        const auto *mf_243 = buffer.data(mf + 243 * ncomps + c);
        const auto *mf_244 = buffer.data(mf + 244 * ncomps + c);
        const auto *mf_245 = buffer.data(mf + 245 * ncomps + c);
        const auto *mf_246 = buffer.data(mf + 246 * ncomps + c);
        const auto *mf_247 = buffer.data(mf + 247 * ncomps + c);
        const auto *mf_248 = buffer.data(mf + 248 * ncomps + c);
        const auto *mf_249 = buffer.data(mf + 249 * ncomps + c);
        const auto *mf_250 = buffer.data(mf + 250 * ncomps + c);
        const auto *mf_251 = buffer.data(mf + 251 * ncomps + c);
        const auto *mf_252 = buffer.data(mf + 252 * ncomps + c);
        const auto *mf_253 = buffer.data(mf + 253 * ncomps + c);
        const auto *mf_254 = buffer.data(mf + 254 * ncomps + c);
        const auto *mf_255 = buffer.data(mf + 255 * ncomps + c);
        const auto *mf_256 = buffer.data(mf + 256 * ncomps + c);
        const auto *mf_257 = buffer.data(mf + 257 * ncomps + c);
        const auto *mf_258 = buffer.data(mf + 258 * ncomps + c);
        const auto *mf_259 = buffer.data(mf + 259 * ncomps + c);
        const auto *mf_260 = buffer.data(mf + 260 * ncomps + c);
        const auto *mf_261 = buffer.data(mf + 261 * ncomps + c);
        const auto *mf_262 = buffer.data(mf + 262 * ncomps + c);
        const auto *mf_263 = buffer.data(mf + 263 * ncomps + c);
        const auto *mf_264 = buffer.data(mf + 264 * ncomps + c);
        const auto *mf_265 = buffer.data(mf + 265 * ncomps + c);
        const auto *mf_266 = buffer.data(mf + 266 * ncomps + c);
        const auto *mf_267 = buffer.data(mf + 267 * ncomps + c);
        const auto *mf_268 = buffer.data(mf + 268 * ncomps + c);
        const auto *mf_269 = buffer.data(mf + 269 * ncomps + c);
        const auto *mf_270 = buffer.data(mf + 270 * ncomps + c);
        const auto *mf_271 = buffer.data(mf + 271 * ncomps + c);
        const auto *mf_272 = buffer.data(mf + 272 * ncomps + c);
        const auto *mf_273 = buffer.data(mf + 273 * ncomps + c);
        const auto *mf_274 = buffer.data(mf + 274 * ncomps + c);
        const auto *mf_275 = buffer.data(mf + 275 * ncomps + c);
        const auto *mf_276 = buffer.data(mf + 276 * ncomps + c);
        const auto *mf_277 = buffer.data(mf + 277 * ncomps + c);
        const auto *mf_278 = buffer.data(mf + 278 * ncomps + c);
        const auto *mf_279 = buffer.data(mf + 279 * ncomps + c);
        const auto *mf_280 = buffer.data(mf + 280 * ncomps + c);
        const auto *mf_281 = buffer.data(mf + 281 * ncomps + c);
        const auto *mf_282 = buffer.data(mf + 282 * ncomps + c);
        const auto *mf_283 = buffer.data(mf + 283 * ncomps + c);
        const auto *mf_284 = buffer.data(mf + 284 * ncomps + c);
        const auto *mf_285 = buffer.data(mf + 285 * ncomps + c);
        const auto *mf_286 = buffer.data(mf + 286 * ncomps + c);
        const auto *mf_287 = buffer.data(mf + 287 * ncomps + c);
        const auto *mf_288 = buffer.data(mf + 288 * ncomps + c);
        const auto *mf_289 = buffer.data(mf + 289 * ncomps + c);
        const auto *mf_296 = buffer.data(mf + 296 * ncomps + c);
        const auto *mf_297 = buffer.data(mf + 297 * ncomps + c);
        const auto *mf_298 = buffer.data(mf + 298 * ncomps + c);
        const auto *mf_299 = buffer.data(mf + 299 * ncomps + c);
        const auto *mf_306 = buffer.data(mf + 306 * ncomps + c);
        const auto *mf_307 = buffer.data(mf + 307 * ncomps + c);
        const auto *mf_308 = buffer.data(mf + 308 * ncomps + c);
        const auto *mf_309 = buffer.data(mf + 309 * ncomps + c);
        const auto *mf_316 = buffer.data(mf + 316 * ncomps + c);
        const auto *mf_317 = buffer.data(mf + 317 * ncomps + c);
        const auto *mf_318 = buffer.data(mf + 318 * ncomps + c);
        const auto *mf_319 = buffer.data(mf + 319 * ncomps + c);
        const auto *mf_326 = buffer.data(mf + 326 * ncomps + c);
        const auto *mf_327 = buffer.data(mf + 327 * ncomps + c);
        const auto *mf_328 = buffer.data(mf + 328 * ncomps + c);
        const auto *mf_329 = buffer.data(mf + 329 * ncomps + c);
        const auto *mf_336 = buffer.data(mf + 336 * ncomps + c);
        const auto *mf_337 = buffer.data(mf + 337 * ncomps + c);
        const auto *mf_338 = buffer.data(mf + 338 * ncomps + c);
        const auto *mf_339 = buffer.data(mf + 339 * ncomps + c);
        const auto *mf_346 = buffer.data(mf + 346 * ncomps + c);
        const auto *mf_347 = buffer.data(mf + 347 * ncomps + c);
        const auto *mf_348 = buffer.data(mf + 348 * ncomps + c);
        const auto *mf_349 = buffer.data(mf + 349 * ncomps + c);
        const auto *mf_359 = buffer.data(mf + 359 * ncomps + c);
        const auto *mf_366 = buffer.data(mf + 366 * ncomps + c);
        const auto *mf_367 = buffer.data(mf + 367 * ncomps + c);
        const auto *mf_368 = buffer.data(mf + 368 * ncomps + c);
        const auto *mf_369 = buffer.data(mf + 369 * ncomps + c);
        const auto *mf_379 = buffer.data(mf + 379 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, lf_195, lf_196, lf_197, \
                         lf_198, lf_199, mf_195, mf_196, mf_197, mf_198, \
                         mf_199 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = ab_x[k] * lf_195[k]
                       + mf_195[k];

            t_291[k] = ab_x[k] * lf_196[k]
                       + mf_196[k];

            t_292[k] = ab_x[k] * lf_197[k]
                       + mf_197[k];

            t_293[k] = ab_x[k] * lf_198[k]
                       + mf_198[k];

            t_294[k] = ab_x[k] * lf_199[k]
                       + mf_199[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_y, ab_z, lf_196, lf_197, \
                         lf_198, lf_199, mf_256, mf_257, mf_258, mf_259, \
                         mf_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = ab_y[k] * lf_196[k]
                       + mf_256[k];

            t_296[k] = ab_y[k] * lf_197[k]
                       + mf_257[k];

            t_297[k] = ab_y[k] * lf_198[k]
                       + mf_258[k];

            t_298[k] = ab_y[k] * lf_199[k]
                       + mf_259[k];

            t_299[k] = ab_z[k] * lf_199[k]
                       + mf_269[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, lf_200, lf_201, lf_202, \
                         lf_203, lf_204, mf_200, mf_201, mf_202, mf_203, \
                         mf_204 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = ab_x[k] * lf_200[k]
                       + mf_200[k];

            t_301[k] = ab_x[k] * lf_201[k]
                       + mf_201[k];

            t_302[k] = ab_x[k] * lf_202[k]
                       + mf_202[k];

            t_303[k] = ab_x[k] * lf_203[k]
                       + mf_203[k];

            t_304[k] = ab_x[k] * lf_204[k]
                       + mf_204[k];
        }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, lf_205, lf_206, lf_207, \
                         lf_208, lf_209, mf_205, mf_206, mf_207, mf_208, \
                         mf_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_305[k] = ab_x[k] * lf_205[k]
                       + mf_205[k];

            t_306[k] = ab_x[k] * lf_206[k]
                       + mf_206[k];

            t_307[k] = ab_x[k] * lf_207[k]
                       + mf_207[k];

            t_308[k] = ab_x[k] * lf_208[k]
                       + mf_208[k];

            t_309[k] = ab_x[k] * lf_209[k]
                       + mf_209[k];
        }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_y, ab_z, lf_206, lf_207, \
                         lf_208, lf_209, mf_266, mf_267, mf_268, mf_269, \
                         mf_279 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_310[k] = ab_y[k] * lf_206[k]
                       + mf_266[k];

            t_311[k] = ab_y[k] * lf_207[k]
                       + mf_267[k];

            t_312[k] = ab_y[k] * lf_208[k]
                       + mf_268[k];

            t_313[k] = ab_y[k] * lf_209[k]
                       + mf_269[k];

            t_314[k] = ab_z[k] * lf_209[k]
                       + mf_279[k];
        }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_x, lf_210, lf_211, lf_212, \
                         lf_213, lf_214, mf_210, mf_211, mf_212, mf_213, \
                         mf_214 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_315[k] = ab_x[k] * lf_210[k]
                       + mf_210[k];

            t_316[k] = ab_x[k] * lf_211[k]
                       + mf_211[k];

            t_317[k] = ab_x[k] * lf_212[k]
                       + mf_212[k];

            t_318[k] = ab_x[k] * lf_213[k]
                       + mf_213[k];

            t_319[k] = ab_x[k] * lf_214[k]
                       + mf_214[k];
        }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_x, lf_215, lf_216, lf_217, \
                         lf_218, lf_219, mf_215, mf_216, mf_217, mf_218, \
                         mf_219 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_320[k] = ab_x[k] * lf_215[k]
                       + mf_215[k];

            t_321[k] = ab_x[k] * lf_216[k]
                       + mf_216[k];

            t_322[k] = ab_x[k] * lf_217[k]
                       + mf_217[k];

            t_323[k] = ab_x[k] * lf_218[k]
                       + mf_218[k];

            t_324[k] = ab_x[k] * lf_219[k]
                       + mf_219[k];
        }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_y, ab_z, lf_216, lf_217, \
                         lf_218, lf_219, mf_286, mf_287, mf_288, mf_289, \
                         mf_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_325[k] = ab_y[k] * lf_216[k]
                       + mf_286[k];

            t_326[k] = ab_y[k] * lf_217[k]
                       + mf_287[k];

            t_327[k] = ab_y[k] * lf_218[k]
                       + mf_288[k];

            t_328[k] = ab_y[k] * lf_219[k]
                       + mf_289[k];

            t_329[k] = ab_z[k] * lf_219[k]
                       + mf_299[k];
        }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_x, lf_220, lf_221, lf_222, \
                         lf_223, lf_224, mf_220, mf_221, mf_222, mf_223, \
                         mf_224 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_330[k] = ab_x[k] * lf_220[k]
                       + mf_220[k];

            t_331[k] = ab_x[k] * lf_221[k]
                       + mf_221[k];

            t_332[k] = ab_x[k] * lf_222[k]
                       + mf_222[k];

            t_333[k] = ab_x[k] * lf_223[k]
                       + mf_223[k];

            t_334[k] = ab_x[k] * lf_224[k]
                       + mf_224[k];
        }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_x, lf_225, lf_226, lf_227, \
                         lf_228, lf_229, mf_225, mf_226, mf_227, mf_228, \
                         mf_229 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_335[k] = ab_x[k] * lf_225[k]
                       + mf_225[k];

            t_336[k] = ab_x[k] * lf_226[k]
                       + mf_226[k];

            t_337[k] = ab_x[k] * lf_227[k]
                       + mf_227[k];

            t_338[k] = ab_x[k] * lf_228[k]
                       + mf_228[k];

            t_339[k] = ab_x[k] * lf_229[k]
                       + mf_229[k];
        }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_y, ab_z, lf_226, lf_227, \
                         lf_228, lf_229, mf_296, mf_297, mf_298, mf_299, \
                         mf_309 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_340[k] = ab_y[k] * lf_226[k]
                       + mf_296[k];

            t_341[k] = ab_y[k] * lf_227[k]
                       + mf_297[k];

            t_342[k] = ab_y[k] * lf_228[k]
                       + mf_298[k];

            t_343[k] = ab_y[k] * lf_229[k]
                       + mf_299[k];

            t_344[k] = ab_z[k] * lf_229[k]
                       + mf_309[k];
        }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_x, lf_230, lf_231, lf_232, \
                         lf_233, lf_234, mf_230, mf_231, mf_232, mf_233, \
                         mf_234 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_345[k] = ab_x[k] * lf_230[k]
                       + mf_230[k];

            t_346[k] = ab_x[k] * lf_231[k]
                       + mf_231[k];

            t_347[k] = ab_x[k] * lf_232[k]
                       + mf_232[k];

            t_348[k] = ab_x[k] * lf_233[k]
                       + mf_233[k];

            t_349[k] = ab_x[k] * lf_234[k]
                       + mf_234[k];
        }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_x, lf_235, lf_236, lf_237, \
                         lf_238, lf_239, mf_235, mf_236, mf_237, mf_238, \
                         mf_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_350[k] = ab_x[k] * lf_235[k]
                       + mf_235[k];

            t_351[k] = ab_x[k] * lf_236[k]
                       + mf_236[k];

            t_352[k] = ab_x[k] * lf_237[k]
                       + mf_237[k];

            t_353[k] = ab_x[k] * lf_238[k]
                       + mf_238[k];

            t_354[k] = ab_x[k] * lf_239[k]
                       + mf_239[k];
        }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_y, ab_z, lf_236, lf_237, \
                         lf_238, lf_239, mf_306, mf_307, mf_308, mf_309, \
                         mf_319 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_355[k] = ab_y[k] * lf_236[k]
                       + mf_306[k];

            t_356[k] = ab_y[k] * lf_237[k]
                       + mf_307[k];

            t_357[k] = ab_y[k] * lf_238[k]
                       + mf_308[k];

            t_358[k] = ab_y[k] * lf_239[k]
                       + mf_309[k];

            t_359[k] = ab_z[k] * lf_239[k]
                       + mf_319[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, lf_240, lf_241, lf_242, \
                         lf_243, lf_244, mf_240, mf_241, mf_242, mf_243, \
                         mf_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = ab_x[k] * lf_240[k]
                       + mf_240[k];

            t_361[k] = ab_x[k] * lf_241[k]
                       + mf_241[k];

            t_362[k] = ab_x[k] * lf_242[k]
                       + mf_242[k];

            t_363[k] = ab_x[k] * lf_243[k]
                       + mf_243[k];

            t_364[k] = ab_x[k] * lf_244[k]
                       + mf_244[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, lf_245, lf_246, lf_247, \
                         lf_248, lf_249, mf_245, mf_246, mf_247, mf_248, \
                         mf_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = ab_x[k] * lf_245[k]
                       + mf_245[k];

            t_366[k] = ab_x[k] * lf_246[k]
                       + mf_246[k];

            t_367[k] = ab_x[k] * lf_247[k]
                       + mf_247[k];

            t_368[k] = ab_x[k] * lf_248[k]
                       + mf_248[k];

            t_369[k] = ab_x[k] * lf_249[k]
                       + mf_249[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_y, ab_z, lf_246, lf_247, \
                         lf_248, lf_249, mf_316, mf_317, mf_318, mf_319, \
                         mf_329 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = ab_y[k] * lf_246[k]
                       + mf_316[k];

            t_371[k] = ab_y[k] * lf_247[k]
                       + mf_317[k];

            t_372[k] = ab_y[k] * lf_248[k]
                       + mf_318[k];

            t_373[k] = ab_y[k] * lf_249[k]
                       + mf_319[k];

            t_374[k] = ab_z[k] * lf_249[k]
                       + mf_329[k];
        }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_x, lf_250, lf_251, lf_252, \
                         lf_253, lf_254, mf_250, mf_251, mf_252, mf_253, \
                         mf_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_375[k] = ab_x[k] * lf_250[k]
                       + mf_250[k];

            t_376[k] = ab_x[k] * lf_251[k]
                       + mf_251[k];

            t_377[k] = ab_x[k] * lf_252[k]
                       + mf_252[k];

            t_378[k] = ab_x[k] * lf_253[k]
                       + mf_253[k];

            t_379[k] = ab_x[k] * lf_254[k]
                       + mf_254[k];
        }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_x, lf_255, lf_256, lf_257, \
                         lf_258, lf_259, mf_255, mf_256, mf_257, mf_258, \
                         mf_259 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_380[k] = ab_x[k] * lf_255[k]
                       + mf_255[k];

            t_381[k] = ab_x[k] * lf_256[k]
                       + mf_256[k];

            t_382[k] = ab_x[k] * lf_257[k]
                       + mf_257[k];

            t_383[k] = ab_x[k] * lf_258[k]
                       + mf_258[k];

            t_384[k] = ab_x[k] * lf_259[k]
                       + mf_259[k];
        }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_y, ab_z, lf_256, lf_257, \
                         lf_258, lf_259, mf_326, mf_327, mf_328, mf_329, \
                         mf_339 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_385[k] = ab_y[k] * lf_256[k]
                       + mf_326[k];

            t_386[k] = ab_y[k] * lf_257[k]
                       + mf_327[k];

            t_387[k] = ab_y[k] * lf_258[k]
                       + mf_328[k];

            t_388[k] = ab_y[k] * lf_259[k]
                       + mf_329[k];

            t_389[k] = ab_z[k] * lf_259[k]
                       + mf_339[k];
        }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_x, lf_260, lf_261, lf_262, \
                         lf_263, lf_264, mf_260, mf_261, mf_262, mf_263, \
                         mf_264 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_390[k] = ab_x[k] * lf_260[k]
                       + mf_260[k];

            t_391[k] = ab_x[k] * lf_261[k]
                       + mf_261[k];

            t_392[k] = ab_x[k] * lf_262[k]
                       + mf_262[k];

            t_393[k] = ab_x[k] * lf_263[k]
                       + mf_263[k];

            t_394[k] = ab_x[k] * lf_264[k]
                       + mf_264[k];
        }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_x, lf_265, lf_266, lf_267, \
                         lf_268, lf_269, mf_265, mf_266, mf_267, mf_268, \
                         mf_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_395[k] = ab_x[k] * lf_265[k]
                       + mf_265[k];

            t_396[k] = ab_x[k] * lf_266[k]
                       + mf_266[k];

            t_397[k] = ab_x[k] * lf_267[k]
                       + mf_267[k];

            t_398[k] = ab_x[k] * lf_268[k]
                       + mf_268[k];

            t_399[k] = ab_x[k] * lf_269[k]
                       + mf_269[k];
        }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_y, ab_z, lf_266, lf_267, \
                         lf_268, lf_269, mf_336, mf_337, mf_338, mf_339, \
                         mf_349 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_400[k] = ab_y[k] * lf_266[k]
                       + mf_336[k];

            t_401[k] = ab_y[k] * lf_267[k]
                       + mf_337[k];

            t_402[k] = ab_y[k] * lf_268[k]
                       + mf_338[k];

            t_403[k] = ab_y[k] * lf_269[k]
                       + mf_339[k];

            t_404[k] = ab_z[k] * lf_269[k]
                       + mf_349[k];
        }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_x, lf_270, lf_271, lf_272, \
                         lf_273, lf_274, mf_270, mf_271, mf_272, mf_273, \
                         mf_274 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_405[k] = ab_x[k] * lf_270[k]
                       + mf_270[k];

            t_406[k] = ab_x[k] * lf_271[k]
                       + mf_271[k];

            t_407[k] = ab_x[k] * lf_272[k]
                       + mf_272[k];

            t_408[k] = ab_x[k] * lf_273[k]
                       + mf_273[k];

            t_409[k] = ab_x[k] * lf_274[k]
                       + mf_274[k];
        }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_x, lf_275, lf_276, lf_277, \
                         lf_278, lf_279, mf_275, mf_276, mf_277, mf_278, \
                         mf_279 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_410[k] = ab_x[k] * lf_275[k]
                       + mf_275[k];

            t_411[k] = ab_x[k] * lf_276[k]
                       + mf_276[k];

            t_412[k] = ab_x[k] * lf_277[k]
                       + mf_277[k];

            t_413[k] = ab_x[k] * lf_278[k]
                       + mf_278[k];

            t_414[k] = ab_x[k] * lf_279[k]
                       + mf_279[k];
        }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_y, ab_z, lf_276, lf_277, \
                         lf_278, lf_279, mf_346, mf_347, mf_348, mf_349, \
                         mf_359 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_415[k] = ab_y[k] * lf_276[k]
                       + mf_346[k];

            t_416[k] = ab_y[k] * lf_277[k]
                       + mf_347[k];

            t_417[k] = ab_y[k] * lf_278[k]
                       + mf_348[k];

            t_418[k] = ab_y[k] * lf_279[k]
                       + mf_349[k];

            t_419[k] = ab_z[k] * lf_279[k]
                       + mf_359[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_x, lf_280, lf_281, lf_282, \
                         lf_283, lf_284, mf_280, mf_281, mf_282, mf_283, \
                         mf_284 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = ab_x[k] * lf_280[k]
                       + mf_280[k];

            t_421[k] = ab_x[k] * lf_281[k]
                       + mf_281[k];

            t_422[k] = ab_x[k] * lf_282[k]
                       + mf_282[k];

            t_423[k] = ab_x[k] * lf_283[k]
                       + mf_283[k];

            t_424[k] = ab_x[k] * lf_284[k]
                       + mf_284[k];
        }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_x, lf_285, lf_286, lf_287, \
                         lf_288, lf_289, mf_285, mf_286, mf_287, mf_288, \
                         mf_289 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = ab_x[k] * lf_285[k]
                       + mf_285[k];

            t_426[k] = ab_x[k] * lf_286[k]
                       + mf_286[k];

            t_427[k] = ab_x[k] * lf_287[k]
                       + mf_287[k];

            t_428[k] = ab_x[k] * lf_288[k]
                       + mf_288[k];

            t_429[k] = ab_x[k] * lf_289[k]
                       + mf_289[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_y, ab_z, lf_286, lf_287, \
                         lf_288, lf_289, mf_366, mf_367, mf_368, mf_369, \
                         mf_379 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = ab_y[k] * lf_286[k]
                       + mf_366[k];

            t_431[k] = ab_y[k] * lf_287[k]
                       + mf_367[k];

            t_432[k] = ab_y[k] * lf_288[k]
                       + mf_368[k];

            t_433[k] = ab_y[k] * lf_289[k]
                       + mf_369[k];

            t_434[k] = ab_z[k] * lf_289[k]
                       + mf_379[k];
        }
    }
}

static auto
compute_hrr_lg_out_of_first_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t lf, const size_t mf,
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
        auto *t_550 = buffer.data(target + 550 * ncomps + c);
        auto *t_551 = buffer.data(target + 551 * ncomps + c);
        auto *t_552 = buffer.data(target + 552 * ncomps + c);
        auto *t_553 = buffer.data(target + 553 * ncomps + c);
        auto *t_554 = buffer.data(target + 554 * ncomps + c);
        auto *t_555 = buffer.data(target + 555 * ncomps + c);
        auto *t_556 = buffer.data(target + 556 * ncomps + c);
        auto *t_557 = buffer.data(target + 557 * ncomps + c);
        auto *t_558 = buffer.data(target + 558 * ncomps + c);
        auto *t_559 = buffer.data(target + 559 * ncomps + c);
        auto *t_560 = buffer.data(target + 560 * ncomps + c);
        auto *t_561 = buffer.data(target + 561 * ncomps + c);
        auto *t_562 = buffer.data(target + 562 * ncomps + c);
        auto *t_563 = buffer.data(target + 563 * ncomps + c);
        auto *t_564 = buffer.data(target + 564 * ncomps + c);
        auto *t_565 = buffer.data(target + 565 * ncomps + c);
        auto *t_566 = buffer.data(target + 566 * ncomps + c);
        auto *t_567 = buffer.data(target + 567 * ncomps + c);
        auto *t_568 = buffer.data(target + 568 * ncomps + c);
        auto *t_569 = buffer.data(target + 569 * ncomps + c);
        auto *t_570 = buffer.data(target + 570 * ncomps + c);
        auto *t_571 = buffer.data(target + 571 * ncomps + c);
        auto *t_572 = buffer.data(target + 572 * ncomps + c);
        auto *t_573 = buffer.data(target + 573 * ncomps + c);
        auto *t_574 = buffer.data(target + 574 * ncomps + c);
        auto *t_575 = buffer.data(target + 575 * ncomps + c);
        auto *t_576 = buffer.data(target + 576 * ncomps + c);
        auto *t_577 = buffer.data(target + 577 * ncomps + c);
        auto *t_578 = buffer.data(target + 578 * ncomps + c);
        auto *t_579 = buffer.data(target + 579 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *lf_290 = buffer.data(lf + 290 * ncomps + c);
        const auto *lf_291 = buffer.data(lf + 291 * ncomps + c);
        const auto *lf_292 = buffer.data(lf + 292 * ncomps + c);
        const auto *lf_293 = buffer.data(lf + 293 * ncomps + c);
        const auto *lf_294 = buffer.data(lf + 294 * ncomps + c);
        const auto *lf_295 = buffer.data(lf + 295 * ncomps + c);
        const auto *lf_296 = buffer.data(lf + 296 * ncomps + c);
        const auto *lf_297 = buffer.data(lf + 297 * ncomps + c);
        const auto *lf_298 = buffer.data(lf + 298 * ncomps + c);
        const auto *lf_299 = buffer.data(lf + 299 * ncomps + c);
        const auto *lf_300 = buffer.data(lf + 300 * ncomps + c);
        const auto *lf_301 = buffer.data(lf + 301 * ncomps + c);
        const auto *lf_302 = buffer.data(lf + 302 * ncomps + c);
        const auto *lf_303 = buffer.data(lf + 303 * ncomps + c);
        const auto *lf_304 = buffer.data(lf + 304 * ncomps + c);
        const auto *lf_305 = buffer.data(lf + 305 * ncomps + c);
        const auto *lf_306 = buffer.data(lf + 306 * ncomps + c);
        const auto *lf_307 = buffer.data(lf + 307 * ncomps + c);
        const auto *lf_308 = buffer.data(lf + 308 * ncomps + c);
        const auto *lf_309 = buffer.data(lf + 309 * ncomps + c);
        const auto *lf_310 = buffer.data(lf + 310 * ncomps + c);
        const auto *lf_311 = buffer.data(lf + 311 * ncomps + c);
        const auto *lf_312 = buffer.data(lf + 312 * ncomps + c);
        const auto *lf_313 = buffer.data(lf + 313 * ncomps + c);
        const auto *lf_314 = buffer.data(lf + 314 * ncomps + c);
        const auto *lf_315 = buffer.data(lf + 315 * ncomps + c);
        const auto *lf_316 = buffer.data(lf + 316 * ncomps + c);
        const auto *lf_317 = buffer.data(lf + 317 * ncomps + c);
        const auto *lf_318 = buffer.data(lf + 318 * ncomps + c);
        const auto *lf_319 = buffer.data(lf + 319 * ncomps + c);
        const auto *lf_320 = buffer.data(lf + 320 * ncomps + c);
        const auto *lf_321 = buffer.data(lf + 321 * ncomps + c);
        const auto *lf_322 = buffer.data(lf + 322 * ncomps + c);
        const auto *lf_323 = buffer.data(lf + 323 * ncomps + c);
        const auto *lf_324 = buffer.data(lf + 324 * ncomps + c);
        const auto *lf_325 = buffer.data(lf + 325 * ncomps + c);
        const auto *lf_326 = buffer.data(lf + 326 * ncomps + c);
        const auto *lf_327 = buffer.data(lf + 327 * ncomps + c);
        const auto *lf_328 = buffer.data(lf + 328 * ncomps + c);
        const auto *lf_329 = buffer.data(lf + 329 * ncomps + c);
        const auto *lf_330 = buffer.data(lf + 330 * ncomps + c);
        const auto *lf_331 = buffer.data(lf + 331 * ncomps + c);
        const auto *lf_332 = buffer.data(lf + 332 * ncomps + c);
        const auto *lf_333 = buffer.data(lf + 333 * ncomps + c);
        const auto *lf_334 = buffer.data(lf + 334 * ncomps + c);
        const auto *lf_335 = buffer.data(lf + 335 * ncomps + c);
        const auto *lf_336 = buffer.data(lf + 336 * ncomps + c);
        const auto *lf_337 = buffer.data(lf + 337 * ncomps + c);
        const auto *lf_338 = buffer.data(lf + 338 * ncomps + c);
        const auto *lf_339 = buffer.data(lf + 339 * ncomps + c);
        const auto *lf_340 = buffer.data(lf + 340 * ncomps + c);
        const auto *lf_341 = buffer.data(lf + 341 * ncomps + c);
        const auto *lf_342 = buffer.data(lf + 342 * ncomps + c);
        const auto *lf_343 = buffer.data(lf + 343 * ncomps + c);
        const auto *lf_344 = buffer.data(lf + 344 * ncomps + c);
        const auto *lf_345 = buffer.data(lf + 345 * ncomps + c);
        const auto *lf_346 = buffer.data(lf + 346 * ncomps + c);
        const auto *lf_347 = buffer.data(lf + 347 * ncomps + c);
        const auto *lf_348 = buffer.data(lf + 348 * ncomps + c);
        const auto *lf_349 = buffer.data(lf + 349 * ncomps + c);
        const auto *lf_350 = buffer.data(lf + 350 * ncomps + c);
        const auto *lf_351 = buffer.data(lf + 351 * ncomps + c);
        const auto *lf_352 = buffer.data(lf + 352 * ncomps + c);
        const auto *lf_353 = buffer.data(lf + 353 * ncomps + c);
        const auto *lf_354 = buffer.data(lf + 354 * ncomps + c);
        const auto *lf_355 = buffer.data(lf + 355 * ncomps + c);
        const auto *lf_356 = buffer.data(lf + 356 * ncomps + c);
        const auto *lf_357 = buffer.data(lf + 357 * ncomps + c);
        const auto *lf_358 = buffer.data(lf + 358 * ncomps + c);
        const auto *lf_359 = buffer.data(lf + 359 * ncomps + c);
        const auto *lf_360 = buffer.data(lf + 360 * ncomps + c);
        const auto *lf_361 = buffer.data(lf + 361 * ncomps + c);
        const auto *lf_362 = buffer.data(lf + 362 * ncomps + c);
        const auto *lf_363 = buffer.data(lf + 363 * ncomps + c);
        const auto *lf_364 = buffer.data(lf + 364 * ncomps + c);
        const auto *lf_365 = buffer.data(lf + 365 * ncomps + c);
        const auto *lf_366 = buffer.data(lf + 366 * ncomps + c);
        const auto *lf_367 = buffer.data(lf + 367 * ncomps + c);
        const auto *lf_368 = buffer.data(lf + 368 * ncomps + c);
        const auto *lf_369 = buffer.data(lf + 369 * ncomps + c);
        const auto *lf_370 = buffer.data(lf + 370 * ncomps + c);
        const auto *lf_371 = buffer.data(lf + 371 * ncomps + c);
        const auto *lf_372 = buffer.data(lf + 372 * ncomps + c);
        const auto *lf_373 = buffer.data(lf + 373 * ncomps + c);
        const auto *lf_374 = buffer.data(lf + 374 * ncomps + c);
        const auto *lf_375 = buffer.data(lf + 375 * ncomps + c);
        const auto *lf_376 = buffer.data(lf + 376 * ncomps + c);
        const auto *lf_377 = buffer.data(lf + 377 * ncomps + c);
        const auto *lf_378 = buffer.data(lf + 378 * ncomps + c);
        const auto *lf_379 = buffer.data(lf + 379 * ncomps + c);
        const auto *lf_380 = buffer.data(lf + 380 * ncomps + c);
        const auto *lf_381 = buffer.data(lf + 381 * ncomps + c);
        const auto *lf_382 = buffer.data(lf + 382 * ncomps + c);
        const auto *lf_383 = buffer.data(lf + 383 * ncomps + c);
        const auto *lf_384 = buffer.data(lf + 384 * ncomps + c);
        const auto *lf_385 = buffer.data(lf + 385 * ncomps + c);
        const auto *lf_386 = buffer.data(lf + 386 * ncomps + c);
        const auto *lf_387 = buffer.data(lf + 387 * ncomps + c);
        const auto *lf_388 = buffer.data(lf + 388 * ncomps + c);
        const auto *lf_389 = buffer.data(lf + 389 * ncomps + c);

        const auto *mf_290 = buffer.data(mf + 290 * ncomps + c);
        const auto *mf_291 = buffer.data(mf + 291 * ncomps + c);
        const auto *mf_292 = buffer.data(mf + 292 * ncomps + c);
        const auto *mf_293 = buffer.data(mf + 293 * ncomps + c);
        const auto *mf_294 = buffer.data(mf + 294 * ncomps + c);
        const auto *mf_295 = buffer.data(mf + 295 * ncomps + c);
        const auto *mf_296 = buffer.data(mf + 296 * ncomps + c);
        const auto *mf_297 = buffer.data(mf + 297 * ncomps + c);
        const auto *mf_298 = buffer.data(mf + 298 * ncomps + c);
        const auto *mf_299 = buffer.data(mf + 299 * ncomps + c);
        const auto *mf_300 = buffer.data(mf + 300 * ncomps + c);
        const auto *mf_301 = buffer.data(mf + 301 * ncomps + c);
        const auto *mf_302 = buffer.data(mf + 302 * ncomps + c);
        const auto *mf_303 = buffer.data(mf + 303 * ncomps + c);
        const auto *mf_304 = buffer.data(mf + 304 * ncomps + c);
        const auto *mf_305 = buffer.data(mf + 305 * ncomps + c);
        const auto *mf_306 = buffer.data(mf + 306 * ncomps + c);
        const auto *mf_307 = buffer.data(mf + 307 * ncomps + c);
        const auto *mf_308 = buffer.data(mf + 308 * ncomps + c);
        const auto *mf_309 = buffer.data(mf + 309 * ncomps + c);
        const auto *mf_310 = buffer.data(mf + 310 * ncomps + c);
        const auto *mf_311 = buffer.data(mf + 311 * ncomps + c);
        const auto *mf_312 = buffer.data(mf + 312 * ncomps + c);
        const auto *mf_313 = buffer.data(mf + 313 * ncomps + c);
        const auto *mf_314 = buffer.data(mf + 314 * ncomps + c);
        const auto *mf_315 = buffer.data(mf + 315 * ncomps + c);
        const auto *mf_316 = buffer.data(mf + 316 * ncomps + c);
        const auto *mf_317 = buffer.data(mf + 317 * ncomps + c);
        const auto *mf_318 = buffer.data(mf + 318 * ncomps + c);
        const auto *mf_319 = buffer.data(mf + 319 * ncomps + c);
        const auto *mf_320 = buffer.data(mf + 320 * ncomps + c);
        const auto *mf_321 = buffer.data(mf + 321 * ncomps + c);
        const auto *mf_322 = buffer.data(mf + 322 * ncomps + c);
        const auto *mf_323 = buffer.data(mf + 323 * ncomps + c);
        const auto *mf_324 = buffer.data(mf + 324 * ncomps + c);
        const auto *mf_325 = buffer.data(mf + 325 * ncomps + c);
        const auto *mf_326 = buffer.data(mf + 326 * ncomps + c);
        const auto *mf_327 = buffer.data(mf + 327 * ncomps + c);
        const auto *mf_328 = buffer.data(mf + 328 * ncomps + c);
        const auto *mf_329 = buffer.data(mf + 329 * ncomps + c);
        const auto *mf_330 = buffer.data(mf + 330 * ncomps + c);
        const auto *mf_331 = buffer.data(mf + 331 * ncomps + c);
        const auto *mf_332 = buffer.data(mf + 332 * ncomps + c);
        const auto *mf_333 = buffer.data(mf + 333 * ncomps + c);
        const auto *mf_334 = buffer.data(mf + 334 * ncomps + c);
        const auto *mf_335 = buffer.data(mf + 335 * ncomps + c);
        const auto *mf_336 = buffer.data(mf + 336 * ncomps + c);
        const auto *mf_337 = buffer.data(mf + 337 * ncomps + c);
        const auto *mf_338 = buffer.data(mf + 338 * ncomps + c);
        const auto *mf_339 = buffer.data(mf + 339 * ncomps + c);
        const auto *mf_340 = buffer.data(mf + 340 * ncomps + c);
        const auto *mf_341 = buffer.data(mf + 341 * ncomps + c);
        const auto *mf_342 = buffer.data(mf + 342 * ncomps + c);
        const auto *mf_343 = buffer.data(mf + 343 * ncomps + c);
        const auto *mf_344 = buffer.data(mf + 344 * ncomps + c);
        const auto *mf_345 = buffer.data(mf + 345 * ncomps + c);
        const auto *mf_346 = buffer.data(mf + 346 * ncomps + c);
        const auto *mf_347 = buffer.data(mf + 347 * ncomps + c);
        const auto *mf_348 = buffer.data(mf + 348 * ncomps + c);
        const auto *mf_349 = buffer.data(mf + 349 * ncomps + c);
        const auto *mf_350 = buffer.data(mf + 350 * ncomps + c);
        const auto *mf_351 = buffer.data(mf + 351 * ncomps + c);
        const auto *mf_352 = buffer.data(mf + 352 * ncomps + c);
        const auto *mf_353 = buffer.data(mf + 353 * ncomps + c);
        const auto *mf_354 = buffer.data(mf + 354 * ncomps + c);
        const auto *mf_355 = buffer.data(mf + 355 * ncomps + c);
        const auto *mf_356 = buffer.data(mf + 356 * ncomps + c);
        const auto *mf_357 = buffer.data(mf + 357 * ncomps + c);
        const auto *mf_358 = buffer.data(mf + 358 * ncomps + c);
        const auto *mf_359 = buffer.data(mf + 359 * ncomps + c);
        const auto *mf_360 = buffer.data(mf + 360 * ncomps + c);
        const auto *mf_361 = buffer.data(mf + 361 * ncomps + c);
        const auto *mf_362 = buffer.data(mf + 362 * ncomps + c);
        const auto *mf_363 = buffer.data(mf + 363 * ncomps + c);
        const auto *mf_364 = buffer.data(mf + 364 * ncomps + c);
        const auto *mf_365 = buffer.data(mf + 365 * ncomps + c);
        const auto *mf_366 = buffer.data(mf + 366 * ncomps + c);
        const auto *mf_367 = buffer.data(mf + 367 * ncomps + c);
        const auto *mf_368 = buffer.data(mf + 368 * ncomps + c);
        const auto *mf_369 = buffer.data(mf + 369 * ncomps + c);
        const auto *mf_370 = buffer.data(mf + 370 * ncomps + c);
        const auto *mf_371 = buffer.data(mf + 371 * ncomps + c);
        const auto *mf_372 = buffer.data(mf + 372 * ncomps + c);
        const auto *mf_373 = buffer.data(mf + 373 * ncomps + c);
        const auto *mf_374 = buffer.data(mf + 374 * ncomps + c);
        const auto *mf_375 = buffer.data(mf + 375 * ncomps + c);
        const auto *mf_376 = buffer.data(mf + 376 * ncomps + c);
        const auto *mf_377 = buffer.data(mf + 377 * ncomps + c);
        const auto *mf_378 = buffer.data(mf + 378 * ncomps + c);
        const auto *mf_379 = buffer.data(mf + 379 * ncomps + c);
        const auto *mf_380 = buffer.data(mf + 380 * ncomps + c);
        const auto *mf_381 = buffer.data(mf + 381 * ncomps + c);
        const auto *mf_382 = buffer.data(mf + 382 * ncomps + c);
        const auto *mf_383 = buffer.data(mf + 383 * ncomps + c);
        const auto *mf_384 = buffer.data(mf + 384 * ncomps + c);
        const auto *mf_385 = buffer.data(mf + 385 * ncomps + c);
        const auto *mf_386 = buffer.data(mf + 386 * ncomps + c);
        const auto *mf_387 = buffer.data(mf + 387 * ncomps + c);
        const auto *mf_388 = buffer.data(mf + 388 * ncomps + c);
        const auto *mf_389 = buffer.data(mf + 389 * ncomps + c);
        const auto *mf_396 = buffer.data(mf + 396 * ncomps + c);
        const auto *mf_397 = buffer.data(mf + 397 * ncomps + c);
        const auto *mf_398 = buffer.data(mf + 398 * ncomps + c);
        const auto *mf_399 = buffer.data(mf + 399 * ncomps + c);
        const auto *mf_406 = buffer.data(mf + 406 * ncomps + c);
        const auto *mf_407 = buffer.data(mf + 407 * ncomps + c);
        const auto *mf_408 = buffer.data(mf + 408 * ncomps + c);
        const auto *mf_409 = buffer.data(mf + 409 * ncomps + c);
        const auto *mf_416 = buffer.data(mf + 416 * ncomps + c);
        const auto *mf_417 = buffer.data(mf + 417 * ncomps + c);
        const auto *mf_418 = buffer.data(mf + 418 * ncomps + c);
        const auto *mf_419 = buffer.data(mf + 419 * ncomps + c);
        const auto *mf_426 = buffer.data(mf + 426 * ncomps + c);
        const auto *mf_427 = buffer.data(mf + 427 * ncomps + c);
        const auto *mf_428 = buffer.data(mf + 428 * ncomps + c);
        const auto *mf_429 = buffer.data(mf + 429 * ncomps + c);
        const auto *mf_436 = buffer.data(mf + 436 * ncomps + c);
        const auto *mf_437 = buffer.data(mf + 437 * ncomps + c);
        const auto *mf_438 = buffer.data(mf + 438 * ncomps + c);
        const auto *mf_439 = buffer.data(mf + 439 * ncomps + c);
        const auto *mf_449 = buffer.data(mf + 449 * ncomps + c);
        const auto *mf_456 = buffer.data(mf + 456 * ncomps + c);
        const auto *mf_457 = buffer.data(mf + 457 * ncomps + c);
        const auto *mf_458 = buffer.data(mf + 458 * ncomps + c);
        const auto *mf_459 = buffer.data(mf + 459 * ncomps + c);
        const auto *mf_466 = buffer.data(mf + 466 * ncomps + c);
        const auto *mf_467 = buffer.data(mf + 467 * ncomps + c);
        const auto *mf_468 = buffer.data(mf + 468 * ncomps + c);
        const auto *mf_469 = buffer.data(mf + 469 * ncomps + c);
        const auto *mf_479 = buffer.data(mf + 479 * ncomps + c);

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_x, lf_290, lf_291, lf_292, \
                         lf_293, lf_294, mf_290, mf_291, mf_292, mf_293, \
                         mf_294 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = ab_x[k] * lf_290[k]
                       + mf_290[k];

            t_436[k] = ab_x[k] * lf_291[k]
                       + mf_291[k];

            t_437[k] = ab_x[k] * lf_292[k]
                       + mf_292[k];

            t_438[k] = ab_x[k] * lf_293[k]
                       + mf_293[k];

            t_439[k] = ab_x[k] * lf_294[k]
                       + mf_294[k];
        }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ab_x, lf_295, lf_296, lf_297, \
                         lf_298, lf_299, mf_295, mf_296, mf_297, mf_298, \
                         mf_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = ab_x[k] * lf_295[k]
                       + mf_295[k];

            t_441[k] = ab_x[k] * lf_296[k]
                       + mf_296[k];

            t_442[k] = ab_x[k] * lf_297[k]
                       + mf_297[k];

            t_443[k] = ab_x[k] * lf_298[k]
                       + mf_298[k];

            t_444[k] = ab_x[k] * lf_299[k]
                       + mf_299[k];
        }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ab_y, ab_z, lf_296, lf_297, \
                         lf_298, lf_299, mf_376, mf_377, mf_378, mf_379, \
                         mf_389 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_445[k] = ab_y[k] * lf_296[k]
                       + mf_376[k];

            t_446[k] = ab_y[k] * lf_297[k]
                       + mf_377[k];

            t_447[k] = ab_y[k] * lf_298[k]
                       + mf_378[k];

            t_448[k] = ab_y[k] * lf_299[k]
                       + mf_379[k];

            t_449[k] = ab_z[k] * lf_299[k]
                       + mf_389[k];
        }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ab_x, lf_300, lf_301, lf_302, \
                         lf_303, lf_304, mf_300, mf_301, mf_302, mf_303, \
                         mf_304 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_450[k] = ab_x[k] * lf_300[k]
                       + mf_300[k];

            t_451[k] = ab_x[k] * lf_301[k]
                       + mf_301[k];

            t_452[k] = ab_x[k] * lf_302[k]
                       + mf_302[k];

            t_453[k] = ab_x[k] * lf_303[k]
                       + mf_303[k];

            t_454[k] = ab_x[k] * lf_304[k]
                       + mf_304[k];
        }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ab_x, lf_305, lf_306, lf_307, \
                         lf_308, lf_309, mf_305, mf_306, mf_307, mf_308, \
                         mf_309 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_455[k] = ab_x[k] * lf_305[k]
                       + mf_305[k];

            t_456[k] = ab_x[k] * lf_306[k]
                       + mf_306[k];

            t_457[k] = ab_x[k] * lf_307[k]
                       + mf_307[k];

            t_458[k] = ab_x[k] * lf_308[k]
                       + mf_308[k];

            t_459[k] = ab_x[k] * lf_309[k]
                       + mf_309[k];
        }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ab_y, ab_z, lf_306, lf_307, \
                         lf_308, lf_309, mf_386, mf_387, mf_388, mf_389, \
                         mf_399 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_460[k] = ab_y[k] * lf_306[k]
                       + mf_386[k];

            t_461[k] = ab_y[k] * lf_307[k]
                       + mf_387[k];

            t_462[k] = ab_y[k] * lf_308[k]
                       + mf_388[k];

            t_463[k] = ab_y[k] * lf_309[k]
                       + mf_389[k];

            t_464[k] = ab_z[k] * lf_309[k]
                       + mf_399[k];
        }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ab_x, lf_310, lf_311, lf_312, \
                         lf_313, lf_314, mf_310, mf_311, mf_312, mf_313, \
                         mf_314 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_465[k] = ab_x[k] * lf_310[k]
                       + mf_310[k];

            t_466[k] = ab_x[k] * lf_311[k]
                       + mf_311[k];

            t_467[k] = ab_x[k] * lf_312[k]
                       + mf_312[k];

            t_468[k] = ab_x[k] * lf_313[k]
                       + mf_313[k];

            t_469[k] = ab_x[k] * lf_314[k]
                       + mf_314[k];
        }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ab_x, lf_315, lf_316, lf_317, \
                         lf_318, lf_319, mf_315, mf_316, mf_317, mf_318, \
                         mf_319 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_470[k] = ab_x[k] * lf_315[k]
                       + mf_315[k];

            t_471[k] = ab_x[k] * lf_316[k]
                       + mf_316[k];

            t_472[k] = ab_x[k] * lf_317[k]
                       + mf_317[k];

            t_473[k] = ab_x[k] * lf_318[k]
                       + mf_318[k];

            t_474[k] = ab_x[k] * lf_319[k]
                       + mf_319[k];
        }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ab_y, ab_z, lf_316, lf_317, \
                         lf_318, lf_319, mf_396, mf_397, mf_398, mf_399, \
                         mf_409 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_475[k] = ab_y[k] * lf_316[k]
                       + mf_396[k];

            t_476[k] = ab_y[k] * lf_317[k]
                       + mf_397[k];

            t_477[k] = ab_y[k] * lf_318[k]
                       + mf_398[k];

            t_478[k] = ab_y[k] * lf_319[k]
                       + mf_399[k];

            t_479[k] = ab_z[k] * lf_319[k]
                       + mf_409[k];
        }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ab_x, lf_320, lf_321, lf_322, \
                         lf_323, lf_324, mf_320, mf_321, mf_322, mf_323, \
                         mf_324 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_480[k] = ab_x[k] * lf_320[k]
                       + mf_320[k];

            t_481[k] = ab_x[k] * lf_321[k]
                       + mf_321[k];

            t_482[k] = ab_x[k] * lf_322[k]
                       + mf_322[k];

            t_483[k] = ab_x[k] * lf_323[k]
                       + mf_323[k];

            t_484[k] = ab_x[k] * lf_324[k]
                       + mf_324[k];
        }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ab_x, lf_325, lf_326, lf_327, \
                         lf_328, lf_329, mf_325, mf_326, mf_327, mf_328, \
                         mf_329 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_485[k] = ab_x[k] * lf_325[k]
                       + mf_325[k];

            t_486[k] = ab_x[k] * lf_326[k]
                       + mf_326[k];

            t_487[k] = ab_x[k] * lf_327[k]
                       + mf_327[k];

            t_488[k] = ab_x[k] * lf_328[k]
                       + mf_328[k];

            t_489[k] = ab_x[k] * lf_329[k]
                       + mf_329[k];
        }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ab_y, ab_z, lf_326, lf_327, \
                         lf_328, lf_329, mf_406, mf_407, mf_408, mf_409, \
                         mf_419 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_490[k] = ab_y[k] * lf_326[k]
                       + mf_406[k];

            t_491[k] = ab_y[k] * lf_327[k]
                       + mf_407[k];

            t_492[k] = ab_y[k] * lf_328[k]
                       + mf_408[k];

            t_493[k] = ab_y[k] * lf_329[k]
                       + mf_409[k];

            t_494[k] = ab_z[k] * lf_329[k]
                       + mf_419[k];
        }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ab_x, lf_330, lf_331, lf_332, \
                         lf_333, lf_334, mf_330, mf_331, mf_332, mf_333, \
                         mf_334 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_495[k] = ab_x[k] * lf_330[k]
                       + mf_330[k];

            t_496[k] = ab_x[k] * lf_331[k]
                       + mf_331[k];

            t_497[k] = ab_x[k] * lf_332[k]
                       + mf_332[k];

            t_498[k] = ab_x[k] * lf_333[k]
                       + mf_333[k];

            t_499[k] = ab_x[k] * lf_334[k]
                       + mf_334[k];
        }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ab_x, lf_335, lf_336, lf_337, \
                         lf_338, lf_339, mf_335, mf_336, mf_337, mf_338, \
                         mf_339 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_500[k] = ab_x[k] * lf_335[k]
                       + mf_335[k];

            t_501[k] = ab_x[k] * lf_336[k]
                       + mf_336[k];

            t_502[k] = ab_x[k] * lf_337[k]
                       + mf_337[k];

            t_503[k] = ab_x[k] * lf_338[k]
                       + mf_338[k];

            t_504[k] = ab_x[k] * lf_339[k]
                       + mf_339[k];
        }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ab_y, ab_z, lf_336, lf_337, \
                         lf_338, lf_339, mf_416, mf_417, mf_418, mf_419, \
                         mf_429 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_505[k] = ab_y[k] * lf_336[k]
                       + mf_416[k];

            t_506[k] = ab_y[k] * lf_337[k]
                       + mf_417[k];

            t_507[k] = ab_y[k] * lf_338[k]
                       + mf_418[k];

            t_508[k] = ab_y[k] * lf_339[k]
                       + mf_419[k];

            t_509[k] = ab_z[k] * lf_339[k]
                       + mf_429[k];
        }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ab_x, lf_340, lf_341, lf_342, \
                         lf_343, lf_344, mf_340, mf_341, mf_342, mf_343, \
                         mf_344 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_510[k] = ab_x[k] * lf_340[k]
                       + mf_340[k];

            t_511[k] = ab_x[k] * lf_341[k]
                       + mf_341[k];

            t_512[k] = ab_x[k] * lf_342[k]
                       + mf_342[k];

            t_513[k] = ab_x[k] * lf_343[k]
                       + mf_343[k];

            t_514[k] = ab_x[k] * lf_344[k]
                       + mf_344[k];
        }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ab_x, lf_345, lf_346, lf_347, \
                         lf_348, lf_349, mf_345, mf_346, mf_347, mf_348, \
                         mf_349 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_515[k] = ab_x[k] * lf_345[k]
                       + mf_345[k];

            t_516[k] = ab_x[k] * lf_346[k]
                       + mf_346[k];

            t_517[k] = ab_x[k] * lf_347[k]
                       + mf_347[k];

            t_518[k] = ab_x[k] * lf_348[k]
                       + mf_348[k];

            t_519[k] = ab_x[k] * lf_349[k]
                       + mf_349[k];
        }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ab_y, ab_z, lf_346, lf_347, \
                         lf_348, lf_349, mf_426, mf_427, mf_428, mf_429, \
                         mf_439 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_520[k] = ab_y[k] * lf_346[k]
                       + mf_426[k];

            t_521[k] = ab_y[k] * lf_347[k]
                       + mf_427[k];

            t_522[k] = ab_y[k] * lf_348[k]
                       + mf_428[k];

            t_523[k] = ab_y[k] * lf_349[k]
                       + mf_429[k];

            t_524[k] = ab_z[k] * lf_349[k]
                       + mf_439[k];
        }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, ab_x, lf_350, lf_351, lf_352, \
                         lf_353, lf_354, mf_350, mf_351, mf_352, mf_353, \
                         mf_354 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_525[k] = ab_x[k] * lf_350[k]
                       + mf_350[k];

            t_526[k] = ab_x[k] * lf_351[k]
                       + mf_351[k];

            t_527[k] = ab_x[k] * lf_352[k]
                       + mf_352[k];

            t_528[k] = ab_x[k] * lf_353[k]
                       + mf_353[k];

            t_529[k] = ab_x[k] * lf_354[k]
                       + mf_354[k];
        }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, ab_x, lf_355, lf_356, lf_357, \
                         lf_358, lf_359, mf_355, mf_356, mf_357, mf_358, \
                         mf_359 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_530[k] = ab_x[k] * lf_355[k]
                       + mf_355[k];

            t_531[k] = ab_x[k] * lf_356[k]
                       + mf_356[k];

            t_532[k] = ab_x[k] * lf_357[k]
                       + mf_357[k];

            t_533[k] = ab_x[k] * lf_358[k]
                       + mf_358[k];

            t_534[k] = ab_x[k] * lf_359[k]
                       + mf_359[k];
        }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, ab_y, ab_z, lf_356, lf_357, \
                         lf_358, lf_359, mf_436, mf_437, mf_438, mf_439, \
                         mf_449 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_535[k] = ab_y[k] * lf_356[k]
                       + mf_436[k];

            t_536[k] = ab_y[k] * lf_357[k]
                       + mf_437[k];

            t_537[k] = ab_y[k] * lf_358[k]
                       + mf_438[k];

            t_538[k] = ab_y[k] * lf_359[k]
                       + mf_439[k];

            t_539[k] = ab_z[k] * lf_359[k]
                       + mf_449[k];
        }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, ab_x, lf_360, lf_361, lf_362, \
                         lf_363, lf_364, mf_360, mf_361, mf_362, mf_363, \
                         mf_364 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_540[k] = ab_x[k] * lf_360[k]
                       + mf_360[k];

            t_541[k] = ab_x[k] * lf_361[k]
                       + mf_361[k];

            t_542[k] = ab_x[k] * lf_362[k]
                       + mf_362[k];

            t_543[k] = ab_x[k] * lf_363[k]
                       + mf_363[k];

            t_544[k] = ab_x[k] * lf_364[k]
                       + mf_364[k];
        }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, ab_x, lf_365, lf_366, lf_367, \
                         lf_368, lf_369, mf_365, mf_366, mf_367, mf_368, \
                         mf_369 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_545[k] = ab_x[k] * lf_365[k]
                       + mf_365[k];

            t_546[k] = ab_x[k] * lf_366[k]
                       + mf_366[k];

            t_547[k] = ab_x[k] * lf_367[k]
                       + mf_367[k];

            t_548[k] = ab_x[k] * lf_368[k]
                       + mf_368[k];

            t_549[k] = ab_x[k] * lf_369[k]
                       + mf_369[k];
        }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, ab_y, ab_z, lf_366, lf_367, \
                         lf_368, lf_369, mf_456, mf_457, mf_458, mf_459, \
                         mf_469 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_550[k] = ab_y[k] * lf_366[k]
                       + mf_456[k];

            t_551[k] = ab_y[k] * lf_367[k]
                       + mf_457[k];

            t_552[k] = ab_y[k] * lf_368[k]
                       + mf_458[k];

            t_553[k] = ab_y[k] * lf_369[k]
                       + mf_459[k];

            t_554[k] = ab_z[k] * lf_369[k]
                       + mf_469[k];
        }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, ab_x, lf_370, lf_371, lf_372, \
                         lf_373, lf_374, mf_370, mf_371, mf_372, mf_373, \
                         mf_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_555[k] = ab_x[k] * lf_370[k]
                       + mf_370[k];

            t_556[k] = ab_x[k] * lf_371[k]
                       + mf_371[k];

            t_557[k] = ab_x[k] * lf_372[k]
                       + mf_372[k];

            t_558[k] = ab_x[k] * lf_373[k]
                       + mf_373[k];

            t_559[k] = ab_x[k] * lf_374[k]
                       + mf_374[k];
        }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, ab_x, lf_375, lf_376, lf_377, \
                         lf_378, lf_379, mf_375, mf_376, mf_377, mf_378, \
                         mf_379 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_560[k] = ab_x[k] * lf_375[k]
                       + mf_375[k];

            t_561[k] = ab_x[k] * lf_376[k]
                       + mf_376[k];

            t_562[k] = ab_x[k] * lf_377[k]
                       + mf_377[k];

            t_563[k] = ab_x[k] * lf_378[k]
                       + mf_378[k];

            t_564[k] = ab_x[k] * lf_379[k]
                       + mf_379[k];
        }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, ab_y, ab_z, lf_376, lf_377, \
                         lf_378, lf_379, mf_466, mf_467, mf_468, mf_469, \
                         mf_479 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_565[k] = ab_y[k] * lf_376[k]
                       + mf_466[k];

            t_566[k] = ab_y[k] * lf_377[k]
                       + mf_467[k];

            t_567[k] = ab_y[k] * lf_378[k]
                       + mf_468[k];

            t_568[k] = ab_y[k] * lf_379[k]
                       + mf_469[k];

            t_569[k] = ab_z[k] * lf_379[k]
                       + mf_479[k];
        }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, ab_x, lf_380, lf_381, lf_382, \
                         lf_383, lf_384, mf_380, mf_381, mf_382, mf_383, \
                         mf_384 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_570[k] = ab_x[k] * lf_380[k]
                       + mf_380[k];

            t_571[k] = ab_x[k] * lf_381[k]
                       + mf_381[k];

            t_572[k] = ab_x[k] * lf_382[k]
                       + mf_382[k];

            t_573[k] = ab_x[k] * lf_383[k]
                       + mf_383[k];

            t_574[k] = ab_x[k] * lf_384[k]
                       + mf_384[k];
        }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, ab_x, lf_385, lf_386, lf_387, \
                         lf_388, lf_389, mf_385, mf_386, mf_387, mf_388, \
                         mf_389 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_575[k] = ab_x[k] * lf_385[k]
                       + mf_385[k];

            t_576[k] = ab_x[k] * lf_386[k]
                       + mf_386[k];

            t_577[k] = ab_x[k] * lf_387[k]
                       + mf_387[k];

            t_578[k] = ab_x[k] * lf_388[k]
                       + mf_388[k];

            t_579[k] = ab_x[k] * lf_389[k]
                       + mf_389[k];
        }
    }
}

static auto
compute_hrr_lg_out_of_first_piece4(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t lf, const size_t mf,
                                   const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_580 = buffer.data(target + 580 * ncomps + c);
        auto *t_581 = buffer.data(target + 581 * ncomps + c);
        auto *t_582 = buffer.data(target + 582 * ncomps + c);
        auto *t_583 = buffer.data(target + 583 * ncomps + c);
        auto *t_584 = buffer.data(target + 584 * ncomps + c);
        auto *t_585 = buffer.data(target + 585 * ncomps + c);
        auto *t_586 = buffer.data(target + 586 * ncomps + c);
        auto *t_587 = buffer.data(target + 587 * ncomps + c);
        auto *t_588 = buffer.data(target + 588 * ncomps + c);
        auto *t_589 = buffer.data(target + 589 * ncomps + c);
        auto *t_590 = buffer.data(target + 590 * ncomps + c);
        auto *t_591 = buffer.data(target + 591 * ncomps + c);
        auto *t_592 = buffer.data(target + 592 * ncomps + c);
        auto *t_593 = buffer.data(target + 593 * ncomps + c);
        auto *t_594 = buffer.data(target + 594 * ncomps + c);
        auto *t_595 = buffer.data(target + 595 * ncomps + c);
        auto *t_596 = buffer.data(target + 596 * ncomps + c);
        auto *t_597 = buffer.data(target + 597 * ncomps + c);
        auto *t_598 = buffer.data(target + 598 * ncomps + c);
        auto *t_599 = buffer.data(target + 599 * ncomps + c);
        auto *t_600 = buffer.data(target + 600 * ncomps + c);
        auto *t_601 = buffer.data(target + 601 * ncomps + c);
        auto *t_602 = buffer.data(target + 602 * ncomps + c);
        auto *t_603 = buffer.data(target + 603 * ncomps + c);
        auto *t_604 = buffer.data(target + 604 * ncomps + c);
        auto *t_605 = buffer.data(target + 605 * ncomps + c);
        auto *t_606 = buffer.data(target + 606 * ncomps + c);
        auto *t_607 = buffer.data(target + 607 * ncomps + c);
        auto *t_608 = buffer.data(target + 608 * ncomps + c);
        auto *t_609 = buffer.data(target + 609 * ncomps + c);
        auto *t_610 = buffer.data(target + 610 * ncomps + c);
        auto *t_611 = buffer.data(target + 611 * ncomps + c);
        auto *t_612 = buffer.data(target + 612 * ncomps + c);
        auto *t_613 = buffer.data(target + 613 * ncomps + c);
        auto *t_614 = buffer.data(target + 614 * ncomps + c);
        auto *t_615 = buffer.data(target + 615 * ncomps + c);
        auto *t_616 = buffer.data(target + 616 * ncomps + c);
        auto *t_617 = buffer.data(target + 617 * ncomps + c);
        auto *t_618 = buffer.data(target + 618 * ncomps + c);
        auto *t_619 = buffer.data(target + 619 * ncomps + c);
        auto *t_620 = buffer.data(target + 620 * ncomps + c);
        auto *t_621 = buffer.data(target + 621 * ncomps + c);
        auto *t_622 = buffer.data(target + 622 * ncomps + c);
        auto *t_623 = buffer.data(target + 623 * ncomps + c);
        auto *t_624 = buffer.data(target + 624 * ncomps + c);
        auto *t_625 = buffer.data(target + 625 * ncomps + c);
        auto *t_626 = buffer.data(target + 626 * ncomps + c);
        auto *t_627 = buffer.data(target + 627 * ncomps + c);
        auto *t_628 = buffer.data(target + 628 * ncomps + c);
        auto *t_629 = buffer.data(target + 629 * ncomps + c);
        auto *t_630 = buffer.data(target + 630 * ncomps + c);
        auto *t_631 = buffer.data(target + 631 * ncomps + c);
        auto *t_632 = buffer.data(target + 632 * ncomps + c);
        auto *t_633 = buffer.data(target + 633 * ncomps + c);
        auto *t_634 = buffer.data(target + 634 * ncomps + c);
        auto *t_635 = buffer.data(target + 635 * ncomps + c);
        auto *t_636 = buffer.data(target + 636 * ncomps + c);
        auto *t_637 = buffer.data(target + 637 * ncomps + c);
        auto *t_638 = buffer.data(target + 638 * ncomps + c);
        auto *t_639 = buffer.data(target + 639 * ncomps + c);
        auto *t_640 = buffer.data(target + 640 * ncomps + c);
        auto *t_641 = buffer.data(target + 641 * ncomps + c);
        auto *t_642 = buffer.data(target + 642 * ncomps + c);
        auto *t_643 = buffer.data(target + 643 * ncomps + c);
        auto *t_644 = buffer.data(target + 644 * ncomps + c);
        auto *t_645 = buffer.data(target + 645 * ncomps + c);
        auto *t_646 = buffer.data(target + 646 * ncomps + c);
        auto *t_647 = buffer.data(target + 647 * ncomps + c);
        auto *t_648 = buffer.data(target + 648 * ncomps + c);
        auto *t_649 = buffer.data(target + 649 * ncomps + c);
        auto *t_650 = buffer.data(target + 650 * ncomps + c);
        auto *t_651 = buffer.data(target + 651 * ncomps + c);
        auto *t_652 = buffer.data(target + 652 * ncomps + c);
        auto *t_653 = buffer.data(target + 653 * ncomps + c);
        auto *t_654 = buffer.data(target + 654 * ncomps + c);
        auto *t_655 = buffer.data(target + 655 * ncomps + c);
        auto *t_656 = buffer.data(target + 656 * ncomps + c);
        auto *t_657 = buffer.data(target + 657 * ncomps + c);
        auto *t_658 = buffer.data(target + 658 * ncomps + c);
        auto *t_659 = buffer.data(target + 659 * ncomps + c);
        auto *t_660 = buffer.data(target + 660 * ncomps + c);
        auto *t_661 = buffer.data(target + 661 * ncomps + c);
        auto *t_662 = buffer.data(target + 662 * ncomps + c);
        auto *t_663 = buffer.data(target + 663 * ncomps + c);
        auto *t_664 = buffer.data(target + 664 * ncomps + c);
        auto *t_665 = buffer.data(target + 665 * ncomps + c);
        auto *t_666 = buffer.data(target + 666 * ncomps + c);
        auto *t_667 = buffer.data(target + 667 * ncomps + c);
        auto *t_668 = buffer.data(target + 668 * ncomps + c);
        auto *t_669 = buffer.data(target + 669 * ncomps + c);
        auto *t_670 = buffer.data(target + 670 * ncomps + c);
        auto *t_671 = buffer.data(target + 671 * ncomps + c);
        auto *t_672 = buffer.data(target + 672 * ncomps + c);
        auto *t_673 = buffer.data(target + 673 * ncomps + c);
        auto *t_674 = buffer.data(target + 674 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *lf_386 = buffer.data(lf + 386 * ncomps + c);
        const auto *lf_387 = buffer.data(lf + 387 * ncomps + c);
        const auto *lf_388 = buffer.data(lf + 388 * ncomps + c);
        const auto *lf_389 = buffer.data(lf + 389 * ncomps + c);
        const auto *lf_390 = buffer.data(lf + 390 * ncomps + c);
        const auto *lf_391 = buffer.data(lf + 391 * ncomps + c);
        const auto *lf_392 = buffer.data(lf + 392 * ncomps + c);
        const auto *lf_393 = buffer.data(lf + 393 * ncomps + c);
        const auto *lf_394 = buffer.data(lf + 394 * ncomps + c);
        const auto *lf_395 = buffer.data(lf + 395 * ncomps + c);
        const auto *lf_396 = buffer.data(lf + 396 * ncomps + c);
        const auto *lf_397 = buffer.data(lf + 397 * ncomps + c);
        const auto *lf_398 = buffer.data(lf + 398 * ncomps + c);
        const auto *lf_399 = buffer.data(lf + 399 * ncomps + c);
        const auto *lf_400 = buffer.data(lf + 400 * ncomps + c);
        const auto *lf_401 = buffer.data(lf + 401 * ncomps + c);
        const auto *lf_402 = buffer.data(lf + 402 * ncomps + c);
        const auto *lf_403 = buffer.data(lf + 403 * ncomps + c);
        const auto *lf_404 = buffer.data(lf + 404 * ncomps + c);
        const auto *lf_405 = buffer.data(lf + 405 * ncomps + c);
        const auto *lf_406 = buffer.data(lf + 406 * ncomps + c);
        const auto *lf_407 = buffer.data(lf + 407 * ncomps + c);
        const auto *lf_408 = buffer.data(lf + 408 * ncomps + c);
        const auto *lf_409 = buffer.data(lf + 409 * ncomps + c);
        const auto *lf_410 = buffer.data(lf + 410 * ncomps + c);
        const auto *lf_411 = buffer.data(lf + 411 * ncomps + c);
        const auto *lf_412 = buffer.data(lf + 412 * ncomps + c);
        const auto *lf_413 = buffer.data(lf + 413 * ncomps + c);
        const auto *lf_414 = buffer.data(lf + 414 * ncomps + c);
        const auto *lf_415 = buffer.data(lf + 415 * ncomps + c);
        const auto *lf_416 = buffer.data(lf + 416 * ncomps + c);
        const auto *lf_417 = buffer.data(lf + 417 * ncomps + c);
        const auto *lf_418 = buffer.data(lf + 418 * ncomps + c);
        const auto *lf_419 = buffer.data(lf + 419 * ncomps + c);
        const auto *lf_420 = buffer.data(lf + 420 * ncomps + c);
        const auto *lf_421 = buffer.data(lf + 421 * ncomps + c);
        const auto *lf_422 = buffer.data(lf + 422 * ncomps + c);
        const auto *lf_423 = buffer.data(lf + 423 * ncomps + c);
        const auto *lf_424 = buffer.data(lf + 424 * ncomps + c);
        const auto *lf_425 = buffer.data(lf + 425 * ncomps + c);
        const auto *lf_426 = buffer.data(lf + 426 * ncomps + c);
        const auto *lf_427 = buffer.data(lf + 427 * ncomps + c);
        const auto *lf_428 = buffer.data(lf + 428 * ncomps + c);
        const auto *lf_429 = buffer.data(lf + 429 * ncomps + c);
        const auto *lf_430 = buffer.data(lf + 430 * ncomps + c);
        const auto *lf_431 = buffer.data(lf + 431 * ncomps + c);
        const auto *lf_432 = buffer.data(lf + 432 * ncomps + c);
        const auto *lf_433 = buffer.data(lf + 433 * ncomps + c);
        const auto *lf_434 = buffer.data(lf + 434 * ncomps + c);
        const auto *lf_435 = buffer.data(lf + 435 * ncomps + c);
        const auto *lf_436 = buffer.data(lf + 436 * ncomps + c);
        const auto *lf_437 = buffer.data(lf + 437 * ncomps + c);
        const auto *lf_438 = buffer.data(lf + 438 * ncomps + c);
        const auto *lf_439 = buffer.data(lf + 439 * ncomps + c);
        const auto *lf_440 = buffer.data(lf + 440 * ncomps + c);
        const auto *lf_441 = buffer.data(lf + 441 * ncomps + c);
        const auto *lf_442 = buffer.data(lf + 442 * ncomps + c);
        const auto *lf_443 = buffer.data(lf + 443 * ncomps + c);
        const auto *lf_444 = buffer.data(lf + 444 * ncomps + c);
        const auto *lf_445 = buffer.data(lf + 445 * ncomps + c);
        const auto *lf_446 = buffer.data(lf + 446 * ncomps + c);
        const auto *lf_447 = buffer.data(lf + 447 * ncomps + c);
        const auto *lf_448 = buffer.data(lf + 448 * ncomps + c);
        const auto *lf_449 = buffer.data(lf + 449 * ncomps + c);

        const auto *mf_390 = buffer.data(mf + 390 * ncomps + c);
        const auto *mf_391 = buffer.data(mf + 391 * ncomps + c);
        const auto *mf_392 = buffer.data(mf + 392 * ncomps + c);
        const auto *mf_393 = buffer.data(mf + 393 * ncomps + c);
        const auto *mf_394 = buffer.data(mf + 394 * ncomps + c);
        const auto *mf_395 = buffer.data(mf + 395 * ncomps + c);
        const auto *mf_396 = buffer.data(mf + 396 * ncomps + c);
        const auto *mf_397 = buffer.data(mf + 397 * ncomps + c);
        const auto *mf_398 = buffer.data(mf + 398 * ncomps + c);
        const auto *mf_399 = buffer.data(mf + 399 * ncomps + c);
        const auto *mf_400 = buffer.data(mf + 400 * ncomps + c);
        const auto *mf_401 = buffer.data(mf + 401 * ncomps + c);
        const auto *mf_402 = buffer.data(mf + 402 * ncomps + c);
        const auto *mf_403 = buffer.data(mf + 403 * ncomps + c);
        const auto *mf_404 = buffer.data(mf + 404 * ncomps + c);
        const auto *mf_405 = buffer.data(mf + 405 * ncomps + c);
        const auto *mf_406 = buffer.data(mf + 406 * ncomps + c);
        const auto *mf_407 = buffer.data(mf + 407 * ncomps + c);
        const auto *mf_408 = buffer.data(mf + 408 * ncomps + c);
        const auto *mf_409 = buffer.data(mf + 409 * ncomps + c);
        const auto *mf_410 = buffer.data(mf + 410 * ncomps + c);
        const auto *mf_411 = buffer.data(mf + 411 * ncomps + c);
        const auto *mf_412 = buffer.data(mf + 412 * ncomps + c);
        const auto *mf_413 = buffer.data(mf + 413 * ncomps + c);
        const auto *mf_414 = buffer.data(mf + 414 * ncomps + c);
        const auto *mf_415 = buffer.data(mf + 415 * ncomps + c);
        const auto *mf_416 = buffer.data(mf + 416 * ncomps + c);
        const auto *mf_417 = buffer.data(mf + 417 * ncomps + c);
        const auto *mf_418 = buffer.data(mf + 418 * ncomps + c);
        const auto *mf_419 = buffer.data(mf + 419 * ncomps + c);
        const auto *mf_420 = buffer.data(mf + 420 * ncomps + c);
        const auto *mf_421 = buffer.data(mf + 421 * ncomps + c);
        const auto *mf_422 = buffer.data(mf + 422 * ncomps + c);
        const auto *mf_423 = buffer.data(mf + 423 * ncomps + c);
        const auto *mf_424 = buffer.data(mf + 424 * ncomps + c);
        const auto *mf_425 = buffer.data(mf + 425 * ncomps + c);
        const auto *mf_426 = buffer.data(mf + 426 * ncomps + c);
        const auto *mf_427 = buffer.data(mf + 427 * ncomps + c);
        const auto *mf_428 = buffer.data(mf + 428 * ncomps + c);
        const auto *mf_429 = buffer.data(mf + 429 * ncomps + c);
        const auto *mf_430 = buffer.data(mf + 430 * ncomps + c);
        const auto *mf_431 = buffer.data(mf + 431 * ncomps + c);
        const auto *mf_432 = buffer.data(mf + 432 * ncomps + c);
        const auto *mf_433 = buffer.data(mf + 433 * ncomps + c);
        const auto *mf_434 = buffer.data(mf + 434 * ncomps + c);
        const auto *mf_435 = buffer.data(mf + 435 * ncomps + c);
        const auto *mf_436 = buffer.data(mf + 436 * ncomps + c);
        const auto *mf_437 = buffer.data(mf + 437 * ncomps + c);
        const auto *mf_438 = buffer.data(mf + 438 * ncomps + c);
        const auto *mf_439 = buffer.data(mf + 439 * ncomps + c);
        const auto *mf_440 = buffer.data(mf + 440 * ncomps + c);
        const auto *mf_441 = buffer.data(mf + 441 * ncomps + c);
        const auto *mf_442 = buffer.data(mf + 442 * ncomps + c);
        const auto *mf_443 = buffer.data(mf + 443 * ncomps + c);
        const auto *mf_444 = buffer.data(mf + 444 * ncomps + c);
        const auto *mf_445 = buffer.data(mf + 445 * ncomps + c);
        const auto *mf_446 = buffer.data(mf + 446 * ncomps + c);
        const auto *mf_447 = buffer.data(mf + 447 * ncomps + c);
        const auto *mf_448 = buffer.data(mf + 448 * ncomps + c);
        const auto *mf_449 = buffer.data(mf + 449 * ncomps + c);
        const auto *mf_476 = buffer.data(mf + 476 * ncomps + c);
        const auto *mf_477 = buffer.data(mf + 477 * ncomps + c);
        const auto *mf_478 = buffer.data(mf + 478 * ncomps + c);
        const auto *mf_479 = buffer.data(mf + 479 * ncomps + c);
        const auto *mf_486 = buffer.data(mf + 486 * ncomps + c);
        const auto *mf_487 = buffer.data(mf + 487 * ncomps + c);
        const auto *mf_488 = buffer.data(mf + 488 * ncomps + c);
        const auto *mf_489 = buffer.data(mf + 489 * ncomps + c);
        const auto *mf_496 = buffer.data(mf + 496 * ncomps + c);
        const auto *mf_497 = buffer.data(mf + 497 * ncomps + c);
        const auto *mf_498 = buffer.data(mf + 498 * ncomps + c);
        const auto *mf_499 = buffer.data(mf + 499 * ncomps + c);
        const auto *mf_506 = buffer.data(mf + 506 * ncomps + c);
        const auto *mf_507 = buffer.data(mf + 507 * ncomps + c);
        const auto *mf_508 = buffer.data(mf + 508 * ncomps + c);
        const auto *mf_509 = buffer.data(mf + 509 * ncomps + c);
        const auto *mf_516 = buffer.data(mf + 516 * ncomps + c);
        const auto *mf_517 = buffer.data(mf + 517 * ncomps + c);
        const auto *mf_518 = buffer.data(mf + 518 * ncomps + c);
        const auto *mf_519 = buffer.data(mf + 519 * ncomps + c);
        const auto *mf_526 = buffer.data(mf + 526 * ncomps + c);
        const auto *mf_527 = buffer.data(mf + 527 * ncomps + c);
        const auto *mf_528 = buffer.data(mf + 528 * ncomps + c);
        const auto *mf_529 = buffer.data(mf + 529 * ncomps + c);
        const auto *mf_536 = buffer.data(mf + 536 * ncomps + c);
        const auto *mf_537 = buffer.data(mf + 537 * ncomps + c);
        const auto *mf_538 = buffer.data(mf + 538 * ncomps + c);
        const auto *mf_539 = buffer.data(mf + 539 * ncomps + c);
        const auto *mf_549 = buffer.data(mf + 549 * ncomps + c);

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, ab_y, ab_z, lf_386, lf_387, \
                         lf_388, lf_389, mf_476, mf_477, mf_478, mf_479, \
                         mf_489 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_580[k] = ab_y[k] * lf_386[k]
                       + mf_476[k];

            t_581[k] = ab_y[k] * lf_387[k]
                       + mf_477[k];

            t_582[k] = ab_y[k] * lf_388[k]
                       + mf_478[k];

            t_583[k] = ab_y[k] * lf_389[k]
                       + mf_479[k];

            t_584[k] = ab_z[k] * lf_389[k]
                       + mf_489[k];
        }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, ab_x, lf_390, lf_391, lf_392, \
                         lf_393, lf_394, mf_390, mf_391, mf_392, mf_393, \
                         mf_394 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_585[k] = ab_x[k] * lf_390[k]
                       + mf_390[k];

            t_586[k] = ab_x[k] * lf_391[k]
                       + mf_391[k];

            t_587[k] = ab_x[k] * lf_392[k]
                       + mf_392[k];

            t_588[k] = ab_x[k] * lf_393[k]
                       + mf_393[k];

            t_589[k] = ab_x[k] * lf_394[k]
                       + mf_394[k];
        }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, ab_x, lf_395, lf_396, lf_397, \
                         lf_398, lf_399, mf_395, mf_396, mf_397, mf_398, \
                         mf_399 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_590[k] = ab_x[k] * lf_395[k]
                       + mf_395[k];

            t_591[k] = ab_x[k] * lf_396[k]
                       + mf_396[k];

            t_592[k] = ab_x[k] * lf_397[k]
                       + mf_397[k];

            t_593[k] = ab_x[k] * lf_398[k]
                       + mf_398[k];

            t_594[k] = ab_x[k] * lf_399[k]
                       + mf_399[k];
        }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, ab_y, ab_z, lf_396, lf_397, \
                         lf_398, lf_399, mf_486, mf_487, mf_488, mf_489, \
                         mf_499 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_595[k] = ab_y[k] * lf_396[k]
                       + mf_486[k];

            t_596[k] = ab_y[k] * lf_397[k]
                       + mf_487[k];

            t_597[k] = ab_y[k] * lf_398[k]
                       + mf_488[k];

            t_598[k] = ab_y[k] * lf_399[k]
                       + mf_489[k];

            t_599[k] = ab_z[k] * lf_399[k]
                       + mf_499[k];
        }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, ab_x, lf_400, lf_401, lf_402, \
                         lf_403, lf_404, mf_400, mf_401, mf_402, mf_403, \
                         mf_404 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_600[k] = ab_x[k] * lf_400[k]
                       + mf_400[k];

            t_601[k] = ab_x[k] * lf_401[k]
                       + mf_401[k];

            t_602[k] = ab_x[k] * lf_402[k]
                       + mf_402[k];

            t_603[k] = ab_x[k] * lf_403[k]
                       + mf_403[k];

            t_604[k] = ab_x[k] * lf_404[k]
                       + mf_404[k];
        }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, ab_x, lf_405, lf_406, lf_407, \
                         lf_408, lf_409, mf_405, mf_406, mf_407, mf_408, \
                         mf_409 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_605[k] = ab_x[k] * lf_405[k]
                       + mf_405[k];

            t_606[k] = ab_x[k] * lf_406[k]
                       + mf_406[k];

            t_607[k] = ab_x[k] * lf_407[k]
                       + mf_407[k];

            t_608[k] = ab_x[k] * lf_408[k]
                       + mf_408[k];

            t_609[k] = ab_x[k] * lf_409[k]
                       + mf_409[k];
        }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, ab_y, ab_z, lf_406, lf_407, \
                         lf_408, lf_409, mf_496, mf_497, mf_498, mf_499, \
                         mf_509 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_610[k] = ab_y[k] * lf_406[k]
                       + mf_496[k];

            t_611[k] = ab_y[k] * lf_407[k]
                       + mf_497[k];

            t_612[k] = ab_y[k] * lf_408[k]
                       + mf_498[k];

            t_613[k] = ab_y[k] * lf_409[k]
                       + mf_499[k];

            t_614[k] = ab_z[k] * lf_409[k]
                       + mf_509[k];
        }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, ab_x, lf_410, lf_411, lf_412, \
                         lf_413, lf_414, mf_410, mf_411, mf_412, mf_413, \
                         mf_414 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_615[k] = ab_x[k] * lf_410[k]
                       + mf_410[k];

            t_616[k] = ab_x[k] * lf_411[k]
                       + mf_411[k];

            t_617[k] = ab_x[k] * lf_412[k]
                       + mf_412[k];

            t_618[k] = ab_x[k] * lf_413[k]
                       + mf_413[k];

            t_619[k] = ab_x[k] * lf_414[k]
                       + mf_414[k];
        }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, ab_x, lf_415, lf_416, lf_417, \
                         lf_418, lf_419, mf_415, mf_416, mf_417, mf_418, \
                         mf_419 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_620[k] = ab_x[k] * lf_415[k]
                       + mf_415[k];

            t_621[k] = ab_x[k] * lf_416[k]
                       + mf_416[k];

            t_622[k] = ab_x[k] * lf_417[k]
                       + mf_417[k];

            t_623[k] = ab_x[k] * lf_418[k]
                       + mf_418[k];

            t_624[k] = ab_x[k] * lf_419[k]
                       + mf_419[k];
        }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, ab_y, ab_z, lf_416, lf_417, \
                         lf_418, lf_419, mf_506, mf_507, mf_508, mf_509, \
                         mf_519 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_625[k] = ab_y[k] * lf_416[k]
                       + mf_506[k];

            t_626[k] = ab_y[k] * lf_417[k]
                       + mf_507[k];

            t_627[k] = ab_y[k] * lf_418[k]
                       + mf_508[k];

            t_628[k] = ab_y[k] * lf_419[k]
                       + mf_509[k];

            t_629[k] = ab_z[k] * lf_419[k]
                       + mf_519[k];
        }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, ab_x, lf_420, lf_421, lf_422, \
                         lf_423, lf_424, mf_420, mf_421, mf_422, mf_423, \
                         mf_424 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_630[k] = ab_x[k] * lf_420[k]
                       + mf_420[k];

            t_631[k] = ab_x[k] * lf_421[k]
                       + mf_421[k];

            t_632[k] = ab_x[k] * lf_422[k]
                       + mf_422[k];

            t_633[k] = ab_x[k] * lf_423[k]
                       + mf_423[k];

            t_634[k] = ab_x[k] * lf_424[k]
                       + mf_424[k];
        }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, ab_x, lf_425, lf_426, lf_427, \
                         lf_428, lf_429, mf_425, mf_426, mf_427, mf_428, \
                         mf_429 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_635[k] = ab_x[k] * lf_425[k]
                       + mf_425[k];

            t_636[k] = ab_x[k] * lf_426[k]
                       + mf_426[k];

            t_637[k] = ab_x[k] * lf_427[k]
                       + mf_427[k];

            t_638[k] = ab_x[k] * lf_428[k]
                       + mf_428[k];

            t_639[k] = ab_x[k] * lf_429[k]
                       + mf_429[k];
        }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, ab_y, ab_z, lf_426, lf_427, \
                         lf_428, lf_429, mf_516, mf_517, mf_518, mf_519, \
                         mf_529 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_640[k] = ab_y[k] * lf_426[k]
                       + mf_516[k];

            t_641[k] = ab_y[k] * lf_427[k]
                       + mf_517[k];

            t_642[k] = ab_y[k] * lf_428[k]
                       + mf_518[k];

            t_643[k] = ab_y[k] * lf_429[k]
                       + mf_519[k];

            t_644[k] = ab_z[k] * lf_429[k]
                       + mf_529[k];
        }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, ab_x, lf_430, lf_431, lf_432, \
                         lf_433, lf_434, mf_430, mf_431, mf_432, mf_433, \
                         mf_434 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_645[k] = ab_x[k] * lf_430[k]
                       + mf_430[k];

            t_646[k] = ab_x[k] * lf_431[k]
                       + mf_431[k];

            t_647[k] = ab_x[k] * lf_432[k]
                       + mf_432[k];

            t_648[k] = ab_x[k] * lf_433[k]
                       + mf_433[k];

            t_649[k] = ab_x[k] * lf_434[k]
                       + mf_434[k];
        }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, ab_x, lf_435, lf_436, lf_437, \
                         lf_438, lf_439, mf_435, mf_436, mf_437, mf_438, \
                         mf_439 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_650[k] = ab_x[k] * lf_435[k]
                       + mf_435[k];

            t_651[k] = ab_x[k] * lf_436[k]
                       + mf_436[k];

            t_652[k] = ab_x[k] * lf_437[k]
                       + mf_437[k];

            t_653[k] = ab_x[k] * lf_438[k]
                       + mf_438[k];

            t_654[k] = ab_x[k] * lf_439[k]
                       + mf_439[k];
        }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, ab_y, ab_z, lf_436, lf_437, \
                         lf_438, lf_439, mf_526, mf_527, mf_528, mf_529, \
                         mf_539 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_655[k] = ab_y[k] * lf_436[k]
                       + mf_526[k];

            t_656[k] = ab_y[k] * lf_437[k]
                       + mf_527[k];

            t_657[k] = ab_y[k] * lf_438[k]
                       + mf_528[k];

            t_658[k] = ab_y[k] * lf_439[k]
                       + mf_529[k];

            t_659[k] = ab_z[k] * lf_439[k]
                       + mf_539[k];
        }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, ab_x, lf_440, lf_441, lf_442, \
                         lf_443, lf_444, mf_440, mf_441, mf_442, mf_443, \
                         mf_444 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_660[k] = ab_x[k] * lf_440[k]
                       + mf_440[k];

            t_661[k] = ab_x[k] * lf_441[k]
                       + mf_441[k];

            t_662[k] = ab_x[k] * lf_442[k]
                       + mf_442[k];

            t_663[k] = ab_x[k] * lf_443[k]
                       + mf_443[k];

            t_664[k] = ab_x[k] * lf_444[k]
                       + mf_444[k];
        }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, ab_x, lf_445, lf_446, lf_447, \
                         lf_448, lf_449, mf_445, mf_446, mf_447, mf_448, \
                         mf_449 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_665[k] = ab_x[k] * lf_445[k]
                       + mf_445[k];

            t_666[k] = ab_x[k] * lf_446[k]
                       + mf_446[k];

            t_667[k] = ab_x[k] * lf_447[k]
                       + mf_447[k];

            t_668[k] = ab_x[k] * lf_448[k]
                       + mf_448[k];

            t_669[k] = ab_x[k] * lf_449[k]
                       + mf_449[k];
        }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, ab_y, ab_z, lf_446, lf_447, \
                         lf_448, lf_449, mf_536, mf_537, mf_538, mf_539, \
                         mf_549 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_670[k] = ab_y[k] * lf_446[k]
                       + mf_536[k];

            t_671[k] = ab_y[k] * lf_447[k]
                       + mf_537[k];

            t_672[k] = ab_y[k] * lf_448[k]
                       + mf_538[k];

            t_673[k] = ab_y[k] * lf_449[k]
                       + mf_539[k];

            t_674[k] = ab_z[k] * lf_449[k]
                       + mf_549[k];
        }
    }
}

auto
compute_hrr_lg_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t lf, const size_t mf,
                            const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_lg_out_of_first_piece0(buffer, coordinates, target, lf, mf, ncomps, nmax);

    compute_hrr_lg_out_of_first_piece1(buffer, coordinates, target, lf, mf, ncomps, nmax);

    compute_hrr_lg_out_of_first_piece2(buffer, coordinates, target, lf, mf, ncomps, nmax);

    compute_hrr_lg_out_of_first_piece3(buffer, coordinates, target, lf, mf, ncomps, nmax);

    compute_hrr_lg_out_of_first_piece4(buffer, coordinates, target, lf, mf, ncomps, nmax);
}

}  // namespace simdtrf
