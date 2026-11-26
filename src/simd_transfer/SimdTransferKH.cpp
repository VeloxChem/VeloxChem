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


#include "SimdTransferKH.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_kh_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t kg, const size_t lg,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *kg_0 = buffer.data(kg + 0 * ncomps + c);
        const auto *kg_1 = buffer.data(kg + 1 * ncomps + c);
        const auto *kg_2 = buffer.data(kg + 2 * ncomps + c);
        const auto *kg_3 = buffer.data(kg + 3 * ncomps + c);
        const auto *kg_4 = buffer.data(kg + 4 * ncomps + c);
        const auto *kg_5 = buffer.data(kg + 5 * ncomps + c);
        const auto *kg_6 = buffer.data(kg + 6 * ncomps + c);
        const auto *kg_7 = buffer.data(kg + 7 * ncomps + c);
        const auto *kg_8 = buffer.data(kg + 8 * ncomps + c);
        const auto *kg_9 = buffer.data(kg + 9 * ncomps + c);
        const auto *kg_10 = buffer.data(kg + 10 * ncomps + c);
        const auto *kg_11 = buffer.data(kg + 11 * ncomps + c);
        const auto *kg_12 = buffer.data(kg + 12 * ncomps + c);
        const auto *kg_13 = buffer.data(kg + 13 * ncomps + c);
        const auto *kg_14 = buffer.data(kg + 14 * ncomps + c);
        const auto *kg_15 = buffer.data(kg + 15 * ncomps + c);
        const auto *kg_16 = buffer.data(kg + 16 * ncomps + c);
        const auto *kg_17 = buffer.data(kg + 17 * ncomps + c);
        const auto *kg_18 = buffer.data(kg + 18 * ncomps + c);
        const auto *kg_19 = buffer.data(kg + 19 * ncomps + c);
        const auto *kg_20 = buffer.data(kg + 20 * ncomps + c);
        const auto *kg_21 = buffer.data(kg + 21 * ncomps + c);
        const auto *kg_22 = buffer.data(kg + 22 * ncomps + c);
        const auto *kg_23 = buffer.data(kg + 23 * ncomps + c);
        const auto *kg_24 = buffer.data(kg + 24 * ncomps + c);
        const auto *kg_25 = buffer.data(kg + 25 * ncomps + c);
        const auto *kg_26 = buffer.data(kg + 26 * ncomps + c);
        const auto *kg_27 = buffer.data(kg + 27 * ncomps + c);
        const auto *kg_28 = buffer.data(kg + 28 * ncomps + c);
        const auto *kg_29 = buffer.data(kg + 29 * ncomps + c);
        const auto *kg_30 = buffer.data(kg + 30 * ncomps + c);
        const auto *kg_31 = buffer.data(kg + 31 * ncomps + c);
        const auto *kg_32 = buffer.data(kg + 32 * ncomps + c);
        const auto *kg_33 = buffer.data(kg + 33 * ncomps + c);
        const auto *kg_34 = buffer.data(kg + 34 * ncomps + c);
        const auto *kg_35 = buffer.data(kg + 35 * ncomps + c);
        const auto *kg_36 = buffer.data(kg + 36 * ncomps + c);
        const auto *kg_37 = buffer.data(kg + 37 * ncomps + c);
        const auto *kg_38 = buffer.data(kg + 38 * ncomps + c);
        const auto *kg_39 = buffer.data(kg + 39 * ncomps + c);
        const auto *kg_40 = buffer.data(kg + 40 * ncomps + c);
        const auto *kg_41 = buffer.data(kg + 41 * ncomps + c);
        const auto *kg_42 = buffer.data(kg + 42 * ncomps + c);
        const auto *kg_43 = buffer.data(kg + 43 * ncomps + c);
        const auto *kg_44 = buffer.data(kg + 44 * ncomps + c);
        const auto *kg_45 = buffer.data(kg + 45 * ncomps + c);
        const auto *kg_46 = buffer.data(kg + 46 * ncomps + c);
        const auto *kg_47 = buffer.data(kg + 47 * ncomps + c);
        const auto *kg_48 = buffer.data(kg + 48 * ncomps + c);
        const auto *kg_49 = buffer.data(kg + 49 * ncomps + c);
        const auto *kg_50 = buffer.data(kg + 50 * ncomps + c);
        const auto *kg_51 = buffer.data(kg + 51 * ncomps + c);
        const auto *kg_52 = buffer.data(kg + 52 * ncomps + c);
        const auto *kg_53 = buffer.data(kg + 53 * ncomps + c);
        const auto *kg_54 = buffer.data(kg + 54 * ncomps + c);
        const auto *kg_55 = buffer.data(kg + 55 * ncomps + c);
        const auto *kg_56 = buffer.data(kg + 56 * ncomps + c);
        const auto *kg_57 = buffer.data(kg + 57 * ncomps + c);
        const auto *kg_58 = buffer.data(kg + 58 * ncomps + c);
        const auto *kg_59 = buffer.data(kg + 59 * ncomps + c);
        const auto *kg_60 = buffer.data(kg + 60 * ncomps + c);
        const auto *kg_61 = buffer.data(kg + 61 * ncomps + c);
        const auto *kg_62 = buffer.data(kg + 62 * ncomps + c);
        const auto *kg_63 = buffer.data(kg + 63 * ncomps + c);
        const auto *kg_64 = buffer.data(kg + 64 * ncomps + c);
        const auto *kg_65 = buffer.data(kg + 65 * ncomps + c);
        const auto *kg_66 = buffer.data(kg + 66 * ncomps + c);
        const auto *kg_67 = buffer.data(kg + 67 * ncomps + c);
        const auto *kg_68 = buffer.data(kg + 68 * ncomps + c);
        const auto *kg_69 = buffer.data(kg + 69 * ncomps + c);
        const auto *kg_70 = buffer.data(kg + 70 * ncomps + c);
        const auto *kg_71 = buffer.data(kg + 71 * ncomps + c);
        const auto *kg_72 = buffer.data(kg + 72 * ncomps + c);
        const auto *kg_73 = buffer.data(kg + 73 * ncomps + c);
        const auto *kg_74 = buffer.data(kg + 74 * ncomps + c);
        const auto *kg_75 = buffer.data(kg + 75 * ncomps + c);
        const auto *kg_76 = buffer.data(kg + 76 * ncomps + c);
        const auto *kg_77 = buffer.data(kg + 77 * ncomps + c);
        const auto *kg_78 = buffer.data(kg + 78 * ncomps + c);
        const auto *kg_79 = buffer.data(kg + 79 * ncomps + c);
        const auto *kg_80 = buffer.data(kg + 80 * ncomps + c);
        const auto *kg_81 = buffer.data(kg + 81 * ncomps + c);
        const auto *kg_82 = buffer.data(kg + 82 * ncomps + c);
        const auto *kg_83 = buffer.data(kg + 83 * ncomps + c);
        const auto *kg_84 = buffer.data(kg + 84 * ncomps + c);
        const auto *kg_85 = buffer.data(kg + 85 * ncomps + c);
        const auto *kg_86 = buffer.data(kg + 86 * ncomps + c);
        const auto *kg_87 = buffer.data(kg + 87 * ncomps + c);
        const auto *kg_88 = buffer.data(kg + 88 * ncomps + c);
        const auto *kg_89 = buffer.data(kg + 89 * ncomps + c);
        const auto *kg_90 = buffer.data(kg + 90 * ncomps + c);
        const auto *kg_91 = buffer.data(kg + 91 * ncomps + c);
        const auto *kg_92 = buffer.data(kg + 92 * ncomps + c);
        const auto *kg_93 = buffer.data(kg + 93 * ncomps + c);
        const auto *kg_94 = buffer.data(kg + 94 * ncomps + c);
        const auto *kg_95 = buffer.data(kg + 95 * ncomps + c);
        const auto *kg_96 = buffer.data(kg + 96 * ncomps + c);
        const auto *kg_97 = buffer.data(kg + 97 * ncomps + c);
        const auto *kg_98 = buffer.data(kg + 98 * ncomps + c);
        const auto *kg_99 = buffer.data(kg + 99 * ncomps + c);
        const auto *kg_100 = buffer.data(kg + 100 * ncomps + c);
        const auto *kg_101 = buffer.data(kg + 101 * ncomps + c);
        const auto *kg_102 = buffer.data(kg + 102 * ncomps + c);
        const auto *kg_103 = buffer.data(kg + 103 * ncomps + c);
        const auto *kg_104 = buffer.data(kg + 104 * ncomps + c);

        const auto *lg_0 = buffer.data(lg + 0 * ncomps + c);
        const auto *lg_1 = buffer.data(lg + 1 * ncomps + c);
        const auto *lg_2 = buffer.data(lg + 2 * ncomps + c);
        const auto *lg_3 = buffer.data(lg + 3 * ncomps + c);
        const auto *lg_4 = buffer.data(lg + 4 * ncomps + c);
        const auto *lg_5 = buffer.data(lg + 5 * ncomps + c);
        const auto *lg_6 = buffer.data(lg + 6 * ncomps + c);
        const auto *lg_7 = buffer.data(lg + 7 * ncomps + c);
        const auto *lg_8 = buffer.data(lg + 8 * ncomps + c);
        const auto *lg_9 = buffer.data(lg + 9 * ncomps + c);
        const auto *lg_10 = buffer.data(lg + 10 * ncomps + c);
        const auto *lg_11 = buffer.data(lg + 11 * ncomps + c);
        const auto *lg_12 = buffer.data(lg + 12 * ncomps + c);
        const auto *lg_13 = buffer.data(lg + 13 * ncomps + c);
        const auto *lg_14 = buffer.data(lg + 14 * ncomps + c);
        const auto *lg_15 = buffer.data(lg + 15 * ncomps + c);
        const auto *lg_16 = buffer.data(lg + 16 * ncomps + c);
        const auto *lg_17 = buffer.data(lg + 17 * ncomps + c);
        const auto *lg_18 = buffer.data(lg + 18 * ncomps + c);
        const auto *lg_19 = buffer.data(lg + 19 * ncomps + c);
        const auto *lg_20 = buffer.data(lg + 20 * ncomps + c);
        const auto *lg_21 = buffer.data(lg + 21 * ncomps + c);
        const auto *lg_22 = buffer.data(lg + 22 * ncomps + c);
        const auto *lg_23 = buffer.data(lg + 23 * ncomps + c);
        const auto *lg_24 = buffer.data(lg + 24 * ncomps + c);
        const auto *lg_25 = buffer.data(lg + 25 * ncomps + c);
        const auto *lg_26 = buffer.data(lg + 26 * ncomps + c);
        const auto *lg_27 = buffer.data(lg + 27 * ncomps + c);
        const auto *lg_28 = buffer.data(lg + 28 * ncomps + c);
        const auto *lg_29 = buffer.data(lg + 29 * ncomps + c);
        const auto *lg_30 = buffer.data(lg + 30 * ncomps + c);
        const auto *lg_31 = buffer.data(lg + 31 * ncomps + c);
        const auto *lg_32 = buffer.data(lg + 32 * ncomps + c);
        const auto *lg_33 = buffer.data(lg + 33 * ncomps + c);
        const auto *lg_34 = buffer.data(lg + 34 * ncomps + c);
        const auto *lg_35 = buffer.data(lg + 35 * ncomps + c);
        const auto *lg_36 = buffer.data(lg + 36 * ncomps + c);
        const auto *lg_37 = buffer.data(lg + 37 * ncomps + c);
        const auto *lg_38 = buffer.data(lg + 38 * ncomps + c);
        const auto *lg_39 = buffer.data(lg + 39 * ncomps + c);
        const auto *lg_40 = buffer.data(lg + 40 * ncomps + c);
        const auto *lg_41 = buffer.data(lg + 41 * ncomps + c);
        const auto *lg_42 = buffer.data(lg + 42 * ncomps + c);
        const auto *lg_43 = buffer.data(lg + 43 * ncomps + c);
        const auto *lg_44 = buffer.data(lg + 44 * ncomps + c);
        const auto *lg_45 = buffer.data(lg + 45 * ncomps + c);
        const auto *lg_46 = buffer.data(lg + 46 * ncomps + c);
        const auto *lg_47 = buffer.data(lg + 47 * ncomps + c);
        const auto *lg_48 = buffer.data(lg + 48 * ncomps + c);
        const auto *lg_49 = buffer.data(lg + 49 * ncomps + c);
        const auto *lg_50 = buffer.data(lg + 50 * ncomps + c);
        const auto *lg_51 = buffer.data(lg + 51 * ncomps + c);
        const auto *lg_52 = buffer.data(lg + 52 * ncomps + c);
        const auto *lg_53 = buffer.data(lg + 53 * ncomps + c);
        const auto *lg_54 = buffer.data(lg + 54 * ncomps + c);
        const auto *lg_55 = buffer.data(lg + 55 * ncomps + c);
        const auto *lg_56 = buffer.data(lg + 56 * ncomps + c);
        const auto *lg_57 = buffer.data(lg + 57 * ncomps + c);
        const auto *lg_58 = buffer.data(lg + 58 * ncomps + c);
        const auto *lg_59 = buffer.data(lg + 59 * ncomps + c);
        const auto *lg_60 = buffer.data(lg + 60 * ncomps + c);
        const auto *lg_61 = buffer.data(lg + 61 * ncomps + c);
        const auto *lg_62 = buffer.data(lg + 62 * ncomps + c);
        const auto *lg_63 = buffer.data(lg + 63 * ncomps + c);
        const auto *lg_64 = buffer.data(lg + 64 * ncomps + c);
        const auto *lg_65 = buffer.data(lg + 65 * ncomps + c);
        const auto *lg_66 = buffer.data(lg + 66 * ncomps + c);
        const auto *lg_67 = buffer.data(lg + 67 * ncomps + c);
        const auto *lg_68 = buffer.data(lg + 68 * ncomps + c);
        const auto *lg_69 = buffer.data(lg + 69 * ncomps + c);
        const auto *lg_70 = buffer.data(lg + 70 * ncomps + c);
        const auto *lg_71 = buffer.data(lg + 71 * ncomps + c);
        const auto *lg_72 = buffer.data(lg + 72 * ncomps + c);
        const auto *lg_73 = buffer.data(lg + 73 * ncomps + c);
        const auto *lg_74 = buffer.data(lg + 74 * ncomps + c);
        const auto *lg_75 = buffer.data(lg + 75 * ncomps + c);
        const auto *lg_76 = buffer.data(lg + 76 * ncomps + c);
        const auto *lg_77 = buffer.data(lg + 77 * ncomps + c);
        const auto *lg_78 = buffer.data(lg + 78 * ncomps + c);
        const auto *lg_79 = buffer.data(lg + 79 * ncomps + c);
        const auto *lg_80 = buffer.data(lg + 80 * ncomps + c);
        const auto *lg_81 = buffer.data(lg + 81 * ncomps + c);
        const auto *lg_82 = buffer.data(lg + 82 * ncomps + c);
        const auto *lg_83 = buffer.data(lg + 83 * ncomps + c);
        const auto *lg_84 = buffer.data(lg + 84 * ncomps + c);
        const auto *lg_85 = buffer.data(lg + 85 * ncomps + c);
        const auto *lg_86 = buffer.data(lg + 86 * ncomps + c);
        const auto *lg_87 = buffer.data(lg + 87 * ncomps + c);
        const auto *lg_88 = buffer.data(lg + 88 * ncomps + c);
        const auto *lg_89 = buffer.data(lg + 89 * ncomps + c);
        const auto *lg_90 = buffer.data(lg + 90 * ncomps + c);
        const auto *lg_91 = buffer.data(lg + 91 * ncomps + c);
        const auto *lg_92 = buffer.data(lg + 92 * ncomps + c);
        const auto *lg_93 = buffer.data(lg + 93 * ncomps + c);
        const auto *lg_94 = buffer.data(lg + 94 * ncomps + c);
        const auto *lg_95 = buffer.data(lg + 95 * ncomps + c);
        const auto *lg_96 = buffer.data(lg + 96 * ncomps + c);
        const auto *lg_97 = buffer.data(lg + 97 * ncomps + c);
        const auto *lg_98 = buffer.data(lg + 98 * ncomps + c);
        const auto *lg_99 = buffer.data(lg + 99 * ncomps + c);
        const auto *lg_100 = buffer.data(lg + 100 * ncomps + c);
        const auto *lg_101 = buffer.data(lg + 101 * ncomps + c);
        const auto *lg_102 = buffer.data(lg + 102 * ncomps + c);
        const auto *lg_103 = buffer.data(lg + 103 * ncomps + c);
        const auto *lg_104 = buffer.data(lg + 104 * ncomps + c);
        const auto *lg_115 = buffer.data(lg + 115 * ncomps + c);
        const auto *lg_116 = buffer.data(lg + 116 * ncomps + c);
        const auto *lg_117 = buffer.data(lg + 117 * ncomps + c);
        const auto *lg_118 = buffer.data(lg + 118 * ncomps + c);
        const auto *lg_119 = buffer.data(lg + 119 * ncomps + c);
        const auto *lg_130 = buffer.data(lg + 130 * ncomps + c);
        const auto *lg_131 = buffer.data(lg + 131 * ncomps + c);
        const auto *lg_132 = buffer.data(lg + 132 * ncomps + c);
        const auto *lg_133 = buffer.data(lg + 133 * ncomps + c);
        const auto *lg_134 = buffer.data(lg + 134 * ncomps + c);
        const auto *lg_149 = buffer.data(lg + 149 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, kg_0, kg_1, kg_2, kg_3, kg_4, lg_0, \
                         lg_1, lg_2, lg_3, lg_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * kg_0[k]
                     + lg_0[k];

            t_1[k] = ab_x[k] * kg_1[k]
                     + lg_1[k];

            t_2[k] = ab_x[k] * kg_2[k]
                     + lg_2[k];

            t_3[k] = ab_x[k] * kg_3[k]
                     + lg_3[k];

            t_4[k] = ab_x[k] * kg_4[k]
                     + lg_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, kg_5, kg_6, kg_7, kg_8, kg_9, lg_5, \
                         lg_6, lg_7, lg_8, lg_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * kg_5[k]
                     + lg_5[k];

            t_6[k] = ab_x[k] * kg_6[k]
                     + lg_6[k];

            t_7[k] = ab_x[k] * kg_7[k]
                     + lg_7[k];

            t_8[k] = ab_x[k] * kg_8[k]
                     + lg_8[k];

            t_9[k] = ab_x[k] * kg_9[k]
                     + lg_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, kg_10, kg_11, kg_12, kg_13, \
                         kg_14, lg_10, lg_11, lg_12, lg_13, lg_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_x[k] * kg_10[k]
                      + lg_10[k];

            t_11[k] = ab_x[k] * kg_11[k]
                      + lg_11[k];

            t_12[k] = ab_x[k] * kg_12[k]
                      + lg_12[k];

            t_13[k] = ab_x[k] * kg_13[k]
                      + lg_13[k];

            t_14[k] = ab_x[k] * kg_14[k]
                      + lg_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_y, kg_10, kg_11, kg_12, kg_13, \
                         kg_14, lg_25, lg_26, lg_27, lg_28, lg_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_y[k] * kg_10[k]
                      + lg_25[k];

            t_16[k] = ab_y[k] * kg_11[k]
                      + lg_26[k];

            t_17[k] = ab_y[k] * kg_12[k]
                      + lg_27[k];

            t_18[k] = ab_y[k] * kg_13[k]
                      + lg_28[k];

            t_19[k] = ab_y[k] * kg_14[k]
                      + lg_29[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, ab_x, ab_z, kg_14, kg_15, kg_16, kg_17, \
                         lg_15, lg_16, lg_17, lg_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_z[k] * kg_14[k]
                      + lg_44[k];

            t_21[k] = ab_x[k] * kg_15[k]
                      + lg_15[k];

            t_22[k] = ab_x[k] * kg_16[k]
                      + lg_16[k];

            t_23[k] = ab_x[k] * kg_17[k]
                      + lg_17[k];
        }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, ab_x, kg_18, kg_19, kg_20, kg_21, \
                         kg_22, lg_18, lg_19, lg_20, lg_21, lg_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_24[k] = ab_x[k] * kg_18[k]
                      + lg_18[k];

            t_25[k] = ab_x[k] * kg_19[k]
                      + lg_19[k];

            t_26[k] = ab_x[k] * kg_20[k]
                      + lg_20[k];

            t_27[k] = ab_x[k] * kg_21[k]
                      + lg_21[k];

            t_28[k] = ab_x[k] * kg_22[k]
                      + lg_22[k];
        }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, ab_x, kg_23, kg_24, kg_25, kg_26, \
                         kg_27, lg_23, lg_24, lg_25, lg_26, lg_27 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_29[k] = ab_x[k] * kg_23[k]
                      + lg_23[k];

            t_30[k] = ab_x[k] * kg_24[k]
                      + lg_24[k];

            t_31[k] = ab_x[k] * kg_25[k]
                      + lg_25[k];

            t_32[k] = ab_x[k] * kg_26[k]
                      + lg_26[k];

            t_33[k] = ab_x[k] * kg_27[k]
                      + lg_27[k];
        }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, ab_x, ab_y, kg_25, kg_26, kg_28, kg_29, \
                         lg_28, lg_29, lg_55, lg_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_34[k] = ab_x[k] * kg_28[k]
                      + lg_28[k];

            t_35[k] = ab_x[k] * kg_29[k]
                      + lg_29[k];

            t_36[k] = ab_y[k] * kg_25[k]
                      + lg_55[k];

            t_37[k] = ab_y[k] * kg_26[k]
                      + lg_56[k];
        }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, ab_y, ab_z, kg_27, kg_28, kg_29, lg_57, \
                         lg_58, lg_59, lg_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_38[k] = ab_y[k] * kg_27[k]
                      + lg_57[k];

            t_39[k] = ab_y[k] * kg_28[k]
                      + lg_58[k];

            t_40[k] = ab_y[k] * kg_29[k]
                      + lg_59[k];

            t_41[k] = ab_z[k] * kg_29[k]
                      + lg_74[k];
        }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, ab_x, kg_30, kg_31, kg_32, kg_33, \
                         kg_34, lg_30, lg_31, lg_32, lg_33, lg_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_42[k] = ab_x[k] * kg_30[k]
                      + lg_30[k];

            t_43[k] = ab_x[k] * kg_31[k]
                      + lg_31[k];

            t_44[k] = ab_x[k] * kg_32[k]
                      + lg_32[k];

            t_45[k] = ab_x[k] * kg_33[k]
                      + lg_33[k];

            t_46[k] = ab_x[k] * kg_34[k]
                      + lg_34[k];
        }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, ab_x, kg_35, kg_36, kg_37, kg_38, \
                         kg_39, lg_35, lg_36, lg_37, lg_38, lg_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_47[k] = ab_x[k] * kg_35[k]
                      + lg_35[k];

            t_48[k] = ab_x[k] * kg_36[k]
                      + lg_36[k];

            t_49[k] = ab_x[k] * kg_37[k]
                      + lg_37[k];

            t_50[k] = ab_x[k] * kg_38[k]
                      + lg_38[k];

            t_51[k] = ab_x[k] * kg_39[k]
                      + lg_39[k];
        }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, ab_x, kg_40, kg_41, kg_42, kg_43, \
                         kg_44, lg_40, lg_41, lg_42, lg_43, lg_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_52[k] = ab_x[k] * kg_40[k]
                      + lg_40[k];

            t_53[k] = ab_x[k] * kg_41[k]
                      + lg_41[k];

            t_54[k] = ab_x[k] * kg_42[k]
                      + lg_42[k];

            t_55[k] = ab_x[k] * kg_43[k]
                      + lg_43[k];

            t_56[k] = ab_x[k] * kg_44[k]
                      + lg_44[k];
        }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, ab_y, kg_40, kg_41, kg_42, kg_43, \
                         kg_44, lg_70, lg_71, lg_72, lg_73, lg_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_57[k] = ab_y[k] * kg_40[k]
                      + lg_70[k];

            t_58[k] = ab_y[k] * kg_41[k]
                      + lg_71[k];

            t_59[k] = ab_y[k] * kg_42[k]
                      + lg_72[k];

            t_60[k] = ab_y[k] * kg_43[k]
                      + lg_73[k];

            t_61[k] = ab_y[k] * kg_44[k]
                      + lg_74[k];
        }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, ab_x, ab_z, kg_44, kg_45, kg_46, kg_47, \
                         lg_45, lg_46, lg_47, lg_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_62[k] = ab_z[k] * kg_44[k]
                      + lg_89[k];

            t_63[k] = ab_x[k] * kg_45[k]
                      + lg_45[k];

            t_64[k] = ab_x[k] * kg_46[k]
                      + lg_46[k];

            t_65[k] = ab_x[k] * kg_47[k]
                      + lg_47[k];
        }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, ab_x, kg_48, kg_49, kg_50, kg_51, \
                         kg_52, lg_48, lg_49, lg_50, lg_51, lg_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_66[k] = ab_x[k] * kg_48[k]
                      + lg_48[k];

            t_67[k] = ab_x[k] * kg_49[k]
                      + lg_49[k];

            t_68[k] = ab_x[k] * kg_50[k]
                      + lg_50[k];

            t_69[k] = ab_x[k] * kg_51[k]
                      + lg_51[k];

            t_70[k] = ab_x[k] * kg_52[k]
                      + lg_52[k];
        }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, ab_x, kg_53, kg_54, kg_55, kg_56, \
                         kg_57, lg_53, lg_54, lg_55, lg_56, lg_57 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_71[k] = ab_x[k] * kg_53[k]
                      + lg_53[k];

            t_72[k] = ab_x[k] * kg_54[k]
                      + lg_54[k];

            t_73[k] = ab_x[k] * kg_55[k]
                      + lg_55[k];

            t_74[k] = ab_x[k] * kg_56[k]
                      + lg_56[k];

            t_75[k] = ab_x[k] * kg_57[k]
                      + lg_57[k];
        }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, ab_x, ab_y, kg_55, kg_56, kg_58, kg_59, \
                         lg_58, lg_59, lg_100, lg_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_76[k] = ab_x[k] * kg_58[k]
                      + lg_58[k];

            t_77[k] = ab_x[k] * kg_59[k]
                      + lg_59[k];

            t_78[k] = ab_y[k] * kg_55[k]
                      + lg_100[k];

            t_79[k] = ab_y[k] * kg_56[k]
                      + lg_101[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, ab_y, ab_z, kg_57, kg_58, kg_59, lg_102, \
                         lg_103, lg_104, lg_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_y[k] * kg_57[k]
                      + lg_102[k];

            t_81[k] = ab_y[k] * kg_58[k]
                      + lg_103[k];

            t_82[k] = ab_y[k] * kg_59[k]
                      + lg_104[k];

            t_83[k] = ab_z[k] * kg_59[k]
                      + lg_119[k];
        }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, ab_x, kg_60, kg_61, kg_62, kg_63, \
                         kg_64, lg_60, lg_61, lg_62, lg_63, lg_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_84[k] = ab_x[k] * kg_60[k]
                      + lg_60[k];

            t_85[k] = ab_x[k] * kg_61[k]
                      + lg_61[k];

            t_86[k] = ab_x[k] * kg_62[k]
                      + lg_62[k];

            t_87[k] = ab_x[k] * kg_63[k]
                      + lg_63[k];

            t_88[k] = ab_x[k] * kg_64[k]
                      + lg_64[k];
        }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, ab_x, kg_65, kg_66, kg_67, kg_68, \
                         kg_69, lg_65, lg_66, lg_67, lg_68, lg_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_89[k] = ab_x[k] * kg_65[k]
                      + lg_65[k];

            t_90[k] = ab_x[k] * kg_66[k]
                      + lg_66[k];

            t_91[k] = ab_x[k] * kg_67[k]
                      + lg_67[k];

            t_92[k] = ab_x[k] * kg_68[k]
                      + lg_68[k];

            t_93[k] = ab_x[k] * kg_69[k]
                      + lg_69[k];
        }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, ab_x, kg_70, kg_71, kg_72, kg_73, \
                         kg_74, lg_70, lg_71, lg_72, lg_73, lg_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_94[k] = ab_x[k] * kg_70[k]
                      + lg_70[k];

            t_95[k] = ab_x[k] * kg_71[k]
                      + lg_71[k];

            t_96[k] = ab_x[k] * kg_72[k]
                      + lg_72[k];

            t_97[k] = ab_x[k] * kg_73[k]
                      + lg_73[k];

            t_98[k] = ab_x[k] * kg_74[k]
                      + lg_74[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, ab_y, kg_70, kg_71, kg_72, kg_73, \
                         kg_74, lg_115, lg_116, lg_117, lg_118, \
                         lg_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_y[k] * kg_70[k]
                      + lg_115[k];

            t_100[k] = ab_y[k] * kg_71[k]
                       + lg_116[k];

            t_101[k] = ab_y[k] * kg_72[k]
                       + lg_117[k];

            t_102[k] = ab_y[k] * kg_73[k]
                       + lg_118[k];

            t_103[k] = ab_y[k] * kg_74[k]
                       + lg_119[k];
        }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, ab_x, ab_z, kg_74, kg_75, kg_76, kg_77, \
                         lg_75, lg_76, lg_77, lg_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_104[k] = ab_z[k] * kg_74[k]
                       + lg_134[k];

            t_105[k] = ab_x[k] * kg_75[k]
                       + lg_75[k];

            t_106[k] = ab_x[k] * kg_76[k]
                       + lg_76[k];

            t_107[k] = ab_x[k] * kg_77[k]
                       + lg_77[k];
        }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, ab_x, kg_78, kg_79, kg_80, kg_81, \
                         kg_82, lg_78, lg_79, lg_80, lg_81, lg_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_108[k] = ab_x[k] * kg_78[k]
                       + lg_78[k];

            t_109[k] = ab_x[k] * kg_79[k]
                       + lg_79[k];

            t_110[k] = ab_x[k] * kg_80[k]
                       + lg_80[k];

            t_111[k] = ab_x[k] * kg_81[k]
                       + lg_81[k];

            t_112[k] = ab_x[k] * kg_82[k]
                       + lg_82[k];
        }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, ab_x, kg_83, kg_84, kg_85, kg_86, \
                         kg_87, lg_83, lg_84, lg_85, lg_86, lg_87 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_113[k] = ab_x[k] * kg_83[k]
                       + lg_83[k];

            t_114[k] = ab_x[k] * kg_84[k]
                       + lg_84[k];

            t_115[k] = ab_x[k] * kg_85[k]
                       + lg_85[k];

            t_116[k] = ab_x[k] * kg_86[k]
                       + lg_86[k];

            t_117[k] = ab_x[k] * kg_87[k]
                       + lg_87[k];
        }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, ab_x, ab_y, kg_85, kg_86, kg_88, kg_89, \
                         lg_88, lg_89, lg_130, lg_131 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_118[k] = ab_x[k] * kg_88[k]
                       + lg_88[k];

            t_119[k] = ab_x[k] * kg_89[k]
                       + lg_89[k];

            t_120[k] = ab_y[k] * kg_85[k]
                       + lg_130[k];

            t_121[k] = ab_y[k] * kg_86[k]
                       + lg_131[k];
        }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, ab_y, ab_z, kg_87, kg_88, kg_89, lg_132, \
                         lg_133, lg_134, lg_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_122[k] = ab_y[k] * kg_87[k]
                       + lg_132[k];

            t_123[k] = ab_y[k] * kg_88[k]
                       + lg_133[k];

            t_124[k] = ab_y[k] * kg_89[k]
                       + lg_134[k];

            t_125[k] = ab_z[k] * kg_89[k]
                       + lg_149[k];
        }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, ab_x, kg_90, kg_91, kg_92, kg_93, \
                         kg_94, lg_90, lg_91, lg_92, lg_93, lg_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_126[k] = ab_x[k] * kg_90[k]
                       + lg_90[k];

            t_127[k] = ab_x[k] * kg_91[k]
                       + lg_91[k];

            t_128[k] = ab_x[k] * kg_92[k]
                       + lg_92[k];

            t_129[k] = ab_x[k] * kg_93[k]
                       + lg_93[k];

            t_130[k] = ab_x[k] * kg_94[k]
                       + lg_94[k];
        }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, ab_x, kg_95, kg_96, kg_97, kg_98, \
                         kg_99, lg_95, lg_96, lg_97, lg_98, lg_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_131[k] = ab_x[k] * kg_95[k]
                       + lg_95[k];

            t_132[k] = ab_x[k] * kg_96[k]
                       + lg_96[k];

            t_133[k] = ab_x[k] * kg_97[k]
                       + lg_97[k];

            t_134[k] = ab_x[k] * kg_98[k]
                       + lg_98[k];

            t_135[k] = ab_x[k] * kg_99[k]
                       + lg_99[k];
        }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, ab_x, kg_100, kg_101, kg_102, \
                         kg_103, kg_104, lg_100, lg_101, lg_102, lg_103, \
                         lg_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_136[k] = ab_x[k] * kg_100[k]
                       + lg_100[k];

            t_137[k] = ab_x[k] * kg_101[k]
                       + lg_101[k];

            t_138[k] = ab_x[k] * kg_102[k]
                       + lg_102[k];

            t_139[k] = ab_x[k] * kg_103[k]
                       + lg_103[k];

            t_140[k] = ab_x[k] * kg_104[k]
                       + lg_104[k];
        }
    }
}

static auto
compute_hrr_kh_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t kg, const size_t lg,
                                   const size_t ncomps, const size_t nmax) -> void
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *kg_100 = buffer.data(kg + 100 * ncomps + c);
        const auto *kg_101 = buffer.data(kg + 101 * ncomps + c);
        const auto *kg_102 = buffer.data(kg + 102 * ncomps + c);
        const auto *kg_103 = buffer.data(kg + 103 * ncomps + c);
        const auto *kg_104 = buffer.data(kg + 104 * ncomps + c);
        const auto *kg_105 = buffer.data(kg + 105 * ncomps + c);
        const auto *kg_106 = buffer.data(kg + 106 * ncomps + c);
        const auto *kg_107 = buffer.data(kg + 107 * ncomps + c);
        const auto *kg_108 = buffer.data(kg + 108 * ncomps + c);
        const auto *kg_109 = buffer.data(kg + 109 * ncomps + c);
        const auto *kg_110 = buffer.data(kg + 110 * ncomps + c);
        const auto *kg_111 = buffer.data(kg + 111 * ncomps + c);
        const auto *kg_112 = buffer.data(kg + 112 * ncomps + c);
        const auto *kg_113 = buffer.data(kg + 113 * ncomps + c);
        const auto *kg_114 = buffer.data(kg + 114 * ncomps + c);
        const auto *kg_115 = buffer.data(kg + 115 * ncomps + c);
        const auto *kg_116 = buffer.data(kg + 116 * ncomps + c);
        const auto *kg_117 = buffer.data(kg + 117 * ncomps + c);
        const auto *kg_118 = buffer.data(kg + 118 * ncomps + c);
        const auto *kg_119 = buffer.data(kg + 119 * ncomps + c);
        const auto *kg_120 = buffer.data(kg + 120 * ncomps + c);
        const auto *kg_121 = buffer.data(kg + 121 * ncomps + c);
        const auto *kg_122 = buffer.data(kg + 122 * ncomps + c);
        const auto *kg_123 = buffer.data(kg + 123 * ncomps + c);
        const auto *kg_124 = buffer.data(kg + 124 * ncomps + c);
        const auto *kg_125 = buffer.data(kg + 125 * ncomps + c);
        const auto *kg_126 = buffer.data(kg + 126 * ncomps + c);
        const auto *kg_127 = buffer.data(kg + 127 * ncomps + c);
        const auto *kg_128 = buffer.data(kg + 128 * ncomps + c);
        const auto *kg_129 = buffer.data(kg + 129 * ncomps + c);
        const auto *kg_130 = buffer.data(kg + 130 * ncomps + c);
        const auto *kg_131 = buffer.data(kg + 131 * ncomps + c);
        const auto *kg_132 = buffer.data(kg + 132 * ncomps + c);
        const auto *kg_133 = buffer.data(kg + 133 * ncomps + c);
        const auto *kg_134 = buffer.data(kg + 134 * ncomps + c);
        const auto *kg_135 = buffer.data(kg + 135 * ncomps + c);
        const auto *kg_136 = buffer.data(kg + 136 * ncomps + c);
        const auto *kg_137 = buffer.data(kg + 137 * ncomps + c);
        const auto *kg_138 = buffer.data(kg + 138 * ncomps + c);
        const auto *kg_139 = buffer.data(kg + 139 * ncomps + c);
        const auto *kg_140 = buffer.data(kg + 140 * ncomps + c);
        const auto *kg_141 = buffer.data(kg + 141 * ncomps + c);
        const auto *kg_142 = buffer.data(kg + 142 * ncomps + c);
        const auto *kg_143 = buffer.data(kg + 143 * ncomps + c);
        const auto *kg_144 = buffer.data(kg + 144 * ncomps + c);
        const auto *kg_145 = buffer.data(kg + 145 * ncomps + c);
        const auto *kg_146 = buffer.data(kg + 146 * ncomps + c);
        const auto *kg_147 = buffer.data(kg + 147 * ncomps + c);
        const auto *kg_148 = buffer.data(kg + 148 * ncomps + c);
        const auto *kg_149 = buffer.data(kg + 149 * ncomps + c);
        const auto *kg_150 = buffer.data(kg + 150 * ncomps + c);
        const auto *kg_151 = buffer.data(kg + 151 * ncomps + c);
        const auto *kg_152 = buffer.data(kg + 152 * ncomps + c);
        const auto *kg_153 = buffer.data(kg + 153 * ncomps + c);
        const auto *kg_154 = buffer.data(kg + 154 * ncomps + c);
        const auto *kg_155 = buffer.data(kg + 155 * ncomps + c);
        const auto *kg_156 = buffer.data(kg + 156 * ncomps + c);
        const auto *kg_157 = buffer.data(kg + 157 * ncomps + c);
        const auto *kg_158 = buffer.data(kg + 158 * ncomps + c);
        const auto *kg_159 = buffer.data(kg + 159 * ncomps + c);
        const auto *kg_160 = buffer.data(kg + 160 * ncomps + c);
        const auto *kg_161 = buffer.data(kg + 161 * ncomps + c);
        const auto *kg_162 = buffer.data(kg + 162 * ncomps + c);
        const auto *kg_163 = buffer.data(kg + 163 * ncomps + c);
        const auto *kg_164 = buffer.data(kg + 164 * ncomps + c);
        const auto *kg_165 = buffer.data(kg + 165 * ncomps + c);
        const auto *kg_166 = buffer.data(kg + 166 * ncomps + c);
        const auto *kg_167 = buffer.data(kg + 167 * ncomps + c);
        const auto *kg_168 = buffer.data(kg + 168 * ncomps + c);
        const auto *kg_169 = buffer.data(kg + 169 * ncomps + c);
        const auto *kg_170 = buffer.data(kg + 170 * ncomps + c);
        const auto *kg_171 = buffer.data(kg + 171 * ncomps + c);
        const auto *kg_172 = buffer.data(kg + 172 * ncomps + c);
        const auto *kg_173 = buffer.data(kg + 173 * ncomps + c);
        const auto *kg_174 = buffer.data(kg + 174 * ncomps + c);
        const auto *kg_175 = buffer.data(kg + 175 * ncomps + c);
        const auto *kg_176 = buffer.data(kg + 176 * ncomps + c);
        const auto *kg_177 = buffer.data(kg + 177 * ncomps + c);
        const auto *kg_178 = buffer.data(kg + 178 * ncomps + c);
        const auto *kg_179 = buffer.data(kg + 179 * ncomps + c);
        const auto *kg_180 = buffer.data(kg + 180 * ncomps + c);
        const auto *kg_181 = buffer.data(kg + 181 * ncomps + c);
        const auto *kg_182 = buffer.data(kg + 182 * ncomps + c);
        const auto *kg_183 = buffer.data(kg + 183 * ncomps + c);
        const auto *kg_184 = buffer.data(kg + 184 * ncomps + c);
        const auto *kg_185 = buffer.data(kg + 185 * ncomps + c);
        const auto *kg_186 = buffer.data(kg + 186 * ncomps + c);
        const auto *kg_187 = buffer.data(kg + 187 * ncomps + c);
        const auto *kg_188 = buffer.data(kg + 188 * ncomps + c);
        const auto *kg_189 = buffer.data(kg + 189 * ncomps + c);
        const auto *kg_190 = buffer.data(kg + 190 * ncomps + c);
        const auto *kg_191 = buffer.data(kg + 191 * ncomps + c);
        const auto *kg_192 = buffer.data(kg + 192 * ncomps + c);
        const auto *kg_193 = buffer.data(kg + 193 * ncomps + c);
        const auto *kg_194 = buffer.data(kg + 194 * ncomps + c);
        const auto *kg_195 = buffer.data(kg + 195 * ncomps + c);
        const auto *kg_196 = buffer.data(kg + 196 * ncomps + c);
        const auto *kg_197 = buffer.data(kg + 197 * ncomps + c);
        const auto *kg_198 = buffer.data(kg + 198 * ncomps + c);
        const auto *kg_199 = buffer.data(kg + 199 * ncomps + c);
        const auto *kg_200 = buffer.data(kg + 200 * ncomps + c);
        const auto *kg_201 = buffer.data(kg + 201 * ncomps + c);
        const auto *kg_202 = buffer.data(kg + 202 * ncomps + c);

        const auto *lg_105 = buffer.data(lg + 105 * ncomps + c);
        const auto *lg_106 = buffer.data(lg + 106 * ncomps + c);
        const auto *lg_107 = buffer.data(lg + 107 * ncomps + c);
        const auto *lg_108 = buffer.data(lg + 108 * ncomps + c);
        const auto *lg_109 = buffer.data(lg + 109 * ncomps + c);
        const auto *lg_110 = buffer.data(lg + 110 * ncomps + c);
        const auto *lg_111 = buffer.data(lg + 111 * ncomps + c);
        const auto *lg_112 = buffer.data(lg + 112 * ncomps + c);
        const auto *lg_113 = buffer.data(lg + 113 * ncomps + c);
        const auto *lg_114 = buffer.data(lg + 114 * ncomps + c);
        const auto *lg_115 = buffer.data(lg + 115 * ncomps + c);
        const auto *lg_116 = buffer.data(lg + 116 * ncomps + c);
        const auto *lg_117 = buffer.data(lg + 117 * ncomps + c);
        const auto *lg_118 = buffer.data(lg + 118 * ncomps + c);
        const auto *lg_119 = buffer.data(lg + 119 * ncomps + c);
        const auto *lg_120 = buffer.data(lg + 120 * ncomps + c);
        const auto *lg_121 = buffer.data(lg + 121 * ncomps + c);
        const auto *lg_122 = buffer.data(lg + 122 * ncomps + c);
        const auto *lg_123 = buffer.data(lg + 123 * ncomps + c);
        const auto *lg_124 = buffer.data(lg + 124 * ncomps + c);
        const auto *lg_125 = buffer.data(lg + 125 * ncomps + c);
        const auto *lg_126 = buffer.data(lg + 126 * ncomps + c);
        const auto *lg_127 = buffer.data(lg + 127 * ncomps + c);
        const auto *lg_128 = buffer.data(lg + 128 * ncomps + c);
        const auto *lg_129 = buffer.data(lg + 129 * ncomps + c);
        const auto *lg_130 = buffer.data(lg + 130 * ncomps + c);
        const auto *lg_131 = buffer.data(lg + 131 * ncomps + c);
        const auto *lg_132 = buffer.data(lg + 132 * ncomps + c);
        const auto *lg_133 = buffer.data(lg + 133 * ncomps + c);
        const auto *lg_134 = buffer.data(lg + 134 * ncomps + c);
        const auto *lg_135 = buffer.data(lg + 135 * ncomps + c);
        const auto *lg_136 = buffer.data(lg + 136 * ncomps + c);
        const auto *lg_137 = buffer.data(lg + 137 * ncomps + c);
        const auto *lg_138 = buffer.data(lg + 138 * ncomps + c);
        const auto *lg_139 = buffer.data(lg + 139 * ncomps + c);
        const auto *lg_140 = buffer.data(lg + 140 * ncomps + c);
        const auto *lg_141 = buffer.data(lg + 141 * ncomps + c);
        const auto *lg_142 = buffer.data(lg + 142 * ncomps + c);
        const auto *lg_143 = buffer.data(lg + 143 * ncomps + c);
        const auto *lg_144 = buffer.data(lg + 144 * ncomps + c);
        const auto *lg_145 = buffer.data(lg + 145 * ncomps + c);
        const auto *lg_146 = buffer.data(lg + 146 * ncomps + c);
        const auto *lg_147 = buffer.data(lg + 147 * ncomps + c);
        const auto *lg_148 = buffer.data(lg + 148 * ncomps + c);
        const auto *lg_149 = buffer.data(lg + 149 * ncomps + c);
        const auto *lg_150 = buffer.data(lg + 150 * ncomps + c);
        const auto *lg_151 = buffer.data(lg + 151 * ncomps + c);
        const auto *lg_152 = buffer.data(lg + 152 * ncomps + c);
        const auto *lg_153 = buffer.data(lg + 153 * ncomps + c);
        const auto *lg_154 = buffer.data(lg + 154 * ncomps + c);
        const auto *lg_155 = buffer.data(lg + 155 * ncomps + c);
        const auto *lg_156 = buffer.data(lg + 156 * ncomps + c);
        const auto *lg_157 = buffer.data(lg + 157 * ncomps + c);
        const auto *lg_158 = buffer.data(lg + 158 * ncomps + c);
        const auto *lg_159 = buffer.data(lg + 159 * ncomps + c);
        const auto *lg_160 = buffer.data(lg + 160 * ncomps + c);
        const auto *lg_161 = buffer.data(lg + 161 * ncomps + c);
        const auto *lg_162 = buffer.data(lg + 162 * ncomps + c);
        const auto *lg_163 = buffer.data(lg + 163 * ncomps + c);
        const auto *lg_164 = buffer.data(lg + 164 * ncomps + c);
        const auto *lg_165 = buffer.data(lg + 165 * ncomps + c);
        const auto *lg_166 = buffer.data(lg + 166 * ncomps + c);
        const auto *lg_167 = buffer.data(lg + 167 * ncomps + c);
        const auto *lg_168 = buffer.data(lg + 168 * ncomps + c);
        const auto *lg_169 = buffer.data(lg + 169 * ncomps + c);
        const auto *lg_170 = buffer.data(lg + 170 * ncomps + c);
        const auto *lg_171 = buffer.data(lg + 171 * ncomps + c);
        const auto *lg_172 = buffer.data(lg + 172 * ncomps + c);
        const auto *lg_173 = buffer.data(lg + 173 * ncomps + c);
        const auto *lg_174 = buffer.data(lg + 174 * ncomps + c);
        const auto *lg_175 = buffer.data(lg + 175 * ncomps + c);
        const auto *lg_176 = buffer.data(lg + 176 * ncomps + c);
        const auto *lg_177 = buffer.data(lg + 177 * ncomps + c);
        const auto *lg_178 = buffer.data(lg + 178 * ncomps + c);
        const auto *lg_179 = buffer.data(lg + 179 * ncomps + c);
        const auto *lg_180 = buffer.data(lg + 180 * ncomps + c);
        const auto *lg_181 = buffer.data(lg + 181 * ncomps + c);
        const auto *lg_182 = buffer.data(lg + 182 * ncomps + c);
        const auto *lg_183 = buffer.data(lg + 183 * ncomps + c);
        const auto *lg_184 = buffer.data(lg + 184 * ncomps + c);
        const auto *lg_185 = buffer.data(lg + 185 * ncomps + c);
        const auto *lg_186 = buffer.data(lg + 186 * ncomps + c);
        const auto *lg_187 = buffer.data(lg + 187 * ncomps + c);
        const auto *lg_188 = buffer.data(lg + 188 * ncomps + c);
        const auto *lg_189 = buffer.data(lg + 189 * ncomps + c);
        const auto *lg_190 = buffer.data(lg + 190 * ncomps + c);
        const auto *lg_191 = buffer.data(lg + 191 * ncomps + c);
        const auto *lg_192 = buffer.data(lg + 192 * ncomps + c);
        const auto *lg_193 = buffer.data(lg + 193 * ncomps + c);
        const auto *lg_194 = buffer.data(lg + 194 * ncomps + c);
        const auto *lg_195 = buffer.data(lg + 195 * ncomps + c);
        const auto *lg_196 = buffer.data(lg + 196 * ncomps + c);
        const auto *lg_197 = buffer.data(lg + 197 * ncomps + c);
        const auto *lg_198 = buffer.data(lg + 198 * ncomps + c);
        const auto *lg_199 = buffer.data(lg + 199 * ncomps + c);
        const auto *lg_200 = buffer.data(lg + 200 * ncomps + c);
        const auto *lg_201 = buffer.data(lg + 201 * ncomps + c);
        const auto *lg_202 = buffer.data(lg + 202 * ncomps + c);
        const auto *lg_205 = buffer.data(lg + 205 * ncomps + c);
        const auto *lg_206 = buffer.data(lg + 206 * ncomps + c);
        const auto *lg_207 = buffer.data(lg + 207 * ncomps + c);
        const auto *lg_208 = buffer.data(lg + 208 * ncomps + c);
        const auto *lg_209 = buffer.data(lg + 209 * ncomps + c);
        const auto *lg_224 = buffer.data(lg + 224 * ncomps + c);
        const auto *lg_235 = buffer.data(lg + 235 * ncomps + c);
        const auto *lg_236 = buffer.data(lg + 236 * ncomps + c);
        const auto *lg_237 = buffer.data(lg + 237 * ncomps + c);
        const auto *lg_238 = buffer.data(lg + 238 * ncomps + c);
        const auto *lg_239 = buffer.data(lg + 239 * ncomps + c);
        const auto *lg_250 = buffer.data(lg + 250 * ncomps + c);
        const auto *lg_251 = buffer.data(lg + 251 * ncomps + c);
        const auto *lg_252 = buffer.data(lg + 252 * ncomps + c);
        const auto *lg_253 = buffer.data(lg + 253 * ncomps + c);
        const auto *lg_254 = buffer.data(lg + 254 * ncomps + c);
        const auto *lg_265 = buffer.data(lg + 265 * ncomps + c);
        const auto *lg_266 = buffer.data(lg + 266 * ncomps + c);
        const auto *lg_267 = buffer.data(lg + 267 * ncomps + c);
        const auto *lg_268 = buffer.data(lg + 268 * ncomps + c);
        const auto *lg_269 = buffer.data(lg + 269 * ncomps + c);
        const auto *lg_284 = buffer.data(lg + 284 * ncomps + c);

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, ab_y, kg_100, kg_101, kg_102, \
                         kg_103, kg_104, lg_160, lg_161, lg_162, lg_163, \
                         lg_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_141[k] = ab_y[k] * kg_100[k]
                       + lg_160[k];

            t_142[k] = ab_y[k] * kg_101[k]
                       + lg_161[k];

            t_143[k] = ab_y[k] * kg_102[k]
                       + lg_162[k];

            t_144[k] = ab_y[k] * kg_103[k]
                       + lg_163[k];

            t_145[k] = ab_y[k] * kg_104[k]
                       + lg_164[k];
        }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, ab_x, ab_z, kg_104, kg_105, kg_106, \
                         kg_107, lg_105, lg_106, lg_107, lg_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_146[k] = ab_z[k] * kg_104[k]
                       + lg_179[k];

            t_147[k] = ab_x[k] * kg_105[k]
                       + lg_105[k];

            t_148[k] = ab_x[k] * kg_106[k]
                       + lg_106[k];

            t_149[k] = ab_x[k] * kg_107[k]
                       + lg_107[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, kg_108, kg_109, kg_110, \
                         kg_111, kg_112, lg_108, lg_109, lg_110, lg_111, \
                         lg_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * kg_108[k]
                       + lg_108[k];

            t_151[k] = ab_x[k] * kg_109[k]
                       + lg_109[k];

            t_152[k] = ab_x[k] * kg_110[k]
                       + lg_110[k];

            t_153[k] = ab_x[k] * kg_111[k]
                       + lg_111[k];

            t_154[k] = ab_x[k] * kg_112[k]
                       + lg_112[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, kg_113, kg_114, kg_115, \
                         kg_116, kg_117, lg_113, lg_114, lg_115, lg_116, \
                         lg_117 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * kg_113[k]
                       + lg_113[k];

            t_156[k] = ab_x[k] * kg_114[k]
                       + lg_114[k];

            t_157[k] = ab_x[k] * kg_115[k]
                       + lg_115[k];

            t_158[k] = ab_x[k] * kg_116[k]
                       + lg_116[k];

            t_159[k] = ab_x[k] * kg_117[k]
                       + lg_117[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, ab_x, ab_y, kg_115, kg_116, kg_118, \
                         kg_119, lg_118, lg_119, lg_175, lg_176 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_x[k] * kg_118[k]
                       + lg_118[k];

            t_161[k] = ab_x[k] * kg_119[k]
                       + lg_119[k];

            t_162[k] = ab_y[k] * kg_115[k]
                       + lg_175[k];

            t_163[k] = ab_y[k] * kg_116[k]
                       + lg_176[k];
        }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, ab_y, ab_z, kg_117, kg_118, kg_119, \
                         lg_177, lg_178, lg_179, lg_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_164[k] = ab_y[k] * kg_117[k]
                       + lg_177[k];

            t_165[k] = ab_y[k] * kg_118[k]
                       + lg_178[k];

            t_166[k] = ab_y[k] * kg_119[k]
                       + lg_179[k];

            t_167[k] = ab_z[k] * kg_119[k]
                       + lg_194[k];
        }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, ab_x, kg_120, kg_121, kg_122, \
                         kg_123, kg_124, lg_120, lg_121, lg_122, lg_123, \
                         lg_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_168[k] = ab_x[k] * kg_120[k]
                       + lg_120[k];

            t_169[k] = ab_x[k] * kg_121[k]
                       + lg_121[k];

            t_170[k] = ab_x[k] * kg_122[k]
                       + lg_122[k];

            t_171[k] = ab_x[k] * kg_123[k]
                       + lg_123[k];

            t_172[k] = ab_x[k] * kg_124[k]
                       + lg_124[k];
        }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, ab_x, kg_125, kg_126, kg_127, \
                         kg_128, kg_129, lg_125, lg_126, lg_127, lg_128, \
                         lg_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_173[k] = ab_x[k] * kg_125[k]
                       + lg_125[k];

            t_174[k] = ab_x[k] * kg_126[k]
                       + lg_126[k];

            t_175[k] = ab_x[k] * kg_127[k]
                       + lg_127[k];

            t_176[k] = ab_x[k] * kg_128[k]
                       + lg_128[k];

            t_177[k] = ab_x[k] * kg_129[k]
                       + lg_129[k];
        }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, ab_x, kg_130, kg_131, kg_132, \
                         kg_133, kg_134, lg_130, lg_131, lg_132, lg_133, \
                         lg_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_178[k] = ab_x[k] * kg_130[k]
                       + lg_130[k];

            t_179[k] = ab_x[k] * kg_131[k]
                       + lg_131[k];

            t_180[k] = ab_x[k] * kg_132[k]
                       + lg_132[k];

            t_181[k] = ab_x[k] * kg_133[k]
                       + lg_133[k];

            t_182[k] = ab_x[k] * kg_134[k]
                       + lg_134[k];
        }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, ab_y, kg_130, kg_131, kg_132, \
                         kg_133, kg_134, lg_190, lg_191, lg_192, lg_193, \
                         lg_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_183[k] = ab_y[k] * kg_130[k]
                       + lg_190[k];

            t_184[k] = ab_y[k] * kg_131[k]
                       + lg_191[k];

            t_185[k] = ab_y[k] * kg_132[k]
                       + lg_192[k];

            t_186[k] = ab_y[k] * kg_133[k]
                       + lg_193[k];

            t_187[k] = ab_y[k] * kg_134[k]
                       + lg_194[k];
        }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, ab_x, ab_z, kg_134, kg_135, kg_136, \
                         kg_137, lg_135, lg_136, lg_137, lg_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_188[k] = ab_z[k] * kg_134[k]
                       + lg_209[k];

            t_189[k] = ab_x[k] * kg_135[k]
                       + lg_135[k];

            t_190[k] = ab_x[k] * kg_136[k]
                       + lg_136[k];

            t_191[k] = ab_x[k] * kg_137[k]
                       + lg_137[k];
        }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, ab_x, kg_138, kg_139, kg_140, \
                         kg_141, kg_142, lg_138, lg_139, lg_140, lg_141, \
                         lg_142 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_192[k] = ab_x[k] * kg_138[k]
                       + lg_138[k];

            t_193[k] = ab_x[k] * kg_139[k]
                       + lg_139[k];

            t_194[k] = ab_x[k] * kg_140[k]
                       + lg_140[k];

            t_195[k] = ab_x[k] * kg_141[k]
                       + lg_141[k];

            t_196[k] = ab_x[k] * kg_142[k]
                       + lg_142[k];
        }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, ab_x, kg_143, kg_144, kg_145, \
                         kg_146, kg_147, lg_143, lg_144, lg_145, lg_146, \
                         lg_147 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_197[k] = ab_x[k] * kg_143[k]
                       + lg_143[k];

            t_198[k] = ab_x[k] * kg_144[k]
                       + lg_144[k];

            t_199[k] = ab_x[k] * kg_145[k]
                       + lg_145[k];

            t_200[k] = ab_x[k] * kg_146[k]
                       + lg_146[k];

            t_201[k] = ab_x[k] * kg_147[k]
                       + lg_147[k];
        }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, ab_x, ab_y, kg_145, kg_146, kg_148, \
                         kg_149, lg_148, lg_149, lg_205, lg_206 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_202[k] = ab_x[k] * kg_148[k]
                       + lg_148[k];

            t_203[k] = ab_x[k] * kg_149[k]
                       + lg_149[k];

            t_204[k] = ab_y[k] * kg_145[k]
                       + lg_205[k];

            t_205[k] = ab_y[k] * kg_146[k]
                       + lg_206[k];
        }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, ab_y, ab_z, kg_147, kg_148, kg_149, \
                         lg_207, lg_208, lg_209, lg_224 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_206[k] = ab_y[k] * kg_147[k]
                       + lg_207[k];

            t_207[k] = ab_y[k] * kg_148[k]
                       + lg_208[k];

            t_208[k] = ab_y[k] * kg_149[k]
                       + lg_209[k];

            t_209[k] = ab_z[k] * kg_149[k]
                       + lg_224[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, kg_150, kg_151, kg_152, \
                         kg_153, kg_154, lg_150, lg_151, lg_152, lg_153, \
                         lg_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = ab_x[k] * kg_150[k]
                       + lg_150[k];

            t_211[k] = ab_x[k] * kg_151[k]
                       + lg_151[k];

            t_212[k] = ab_x[k] * kg_152[k]
                       + lg_152[k];

            t_213[k] = ab_x[k] * kg_153[k]
                       + lg_153[k];

            t_214[k] = ab_x[k] * kg_154[k]
                       + lg_154[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, kg_155, kg_156, kg_157, \
                         kg_158, kg_159, lg_155, lg_156, lg_157, lg_158, \
                         lg_159 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = ab_x[k] * kg_155[k]
                       + lg_155[k];

            t_216[k] = ab_x[k] * kg_156[k]
                       + lg_156[k];

            t_217[k] = ab_x[k] * kg_157[k]
                       + lg_157[k];

            t_218[k] = ab_x[k] * kg_158[k]
                       + lg_158[k];

            t_219[k] = ab_x[k] * kg_159[k]
                       + lg_159[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, kg_160, kg_161, kg_162, \
                         kg_163, kg_164, lg_160, lg_161, lg_162, lg_163, \
                         lg_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = ab_x[k] * kg_160[k]
                       + lg_160[k];

            t_221[k] = ab_x[k] * kg_161[k]
                       + lg_161[k];

            t_222[k] = ab_x[k] * kg_162[k]
                       + lg_162[k];

            t_223[k] = ab_x[k] * kg_163[k]
                       + lg_163[k];

            t_224[k] = ab_x[k] * kg_164[k]
                       + lg_164[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_y, kg_160, kg_161, kg_162, \
                         kg_163, kg_164, lg_235, lg_236, lg_237, lg_238, \
                         lg_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = ab_y[k] * kg_160[k]
                       + lg_235[k];

            t_226[k] = ab_y[k] * kg_161[k]
                       + lg_236[k];

            t_227[k] = ab_y[k] * kg_162[k]
                       + lg_237[k];

            t_228[k] = ab_y[k] * kg_163[k]
                       + lg_238[k];

            t_229[k] = ab_y[k] * kg_164[k]
                       + lg_239[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, ab_x, ab_z, kg_164, kg_165, kg_166, \
                         kg_167, lg_165, lg_166, lg_167, lg_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = ab_z[k] * kg_164[k]
                       + lg_254[k];

            t_231[k] = ab_x[k] * kg_165[k]
                       + lg_165[k];

            t_232[k] = ab_x[k] * kg_166[k]
                       + lg_166[k];

            t_233[k] = ab_x[k] * kg_167[k]
                       + lg_167[k];
        }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, ab_x, kg_168, kg_169, kg_170, \
                         kg_171, kg_172, lg_168, lg_169, lg_170, lg_171, \
                         lg_172 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_234[k] = ab_x[k] * kg_168[k]
                       + lg_168[k];

            t_235[k] = ab_x[k] * kg_169[k]
                       + lg_169[k];

            t_236[k] = ab_x[k] * kg_170[k]
                       + lg_170[k];

            t_237[k] = ab_x[k] * kg_171[k]
                       + lg_171[k];

            t_238[k] = ab_x[k] * kg_172[k]
                       + lg_172[k];
        }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, ab_x, kg_173, kg_174, kg_175, \
                         kg_176, kg_177, lg_173, lg_174, lg_175, lg_176, \
                         lg_177 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_239[k] = ab_x[k] * kg_173[k]
                       + lg_173[k];

            t_240[k] = ab_x[k] * kg_174[k]
                       + lg_174[k];

            t_241[k] = ab_x[k] * kg_175[k]
                       + lg_175[k];

            t_242[k] = ab_x[k] * kg_176[k]
                       + lg_176[k];

            t_243[k] = ab_x[k] * kg_177[k]
                       + lg_177[k];
        }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, ab_x, ab_y, kg_175, kg_176, kg_178, \
                         kg_179, lg_178, lg_179, lg_250, lg_251 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_244[k] = ab_x[k] * kg_178[k]
                       + lg_178[k];

            t_245[k] = ab_x[k] * kg_179[k]
                       + lg_179[k];

            t_246[k] = ab_y[k] * kg_175[k]
                       + lg_250[k];

            t_247[k] = ab_y[k] * kg_176[k]
                       + lg_251[k];
        }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, ab_y, ab_z, kg_177, kg_178, kg_179, \
                         lg_252, lg_253, lg_254, lg_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_248[k] = ab_y[k] * kg_177[k]
                       + lg_252[k];

            t_249[k] = ab_y[k] * kg_178[k]
                       + lg_253[k];

            t_250[k] = ab_y[k] * kg_179[k]
                       + lg_254[k];

            t_251[k] = ab_z[k] * kg_179[k]
                       + lg_269[k];
        }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, ab_x, kg_180, kg_181, kg_182, \
                         kg_183, kg_184, lg_180, lg_181, lg_182, lg_183, \
                         lg_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_252[k] = ab_x[k] * kg_180[k]
                       + lg_180[k];

            t_253[k] = ab_x[k] * kg_181[k]
                       + lg_181[k];

            t_254[k] = ab_x[k] * kg_182[k]
                       + lg_182[k];

            t_255[k] = ab_x[k] * kg_183[k]
                       + lg_183[k];

            t_256[k] = ab_x[k] * kg_184[k]
                       + lg_184[k];
        }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, ab_x, kg_185, kg_186, kg_187, \
                         kg_188, kg_189, lg_185, lg_186, lg_187, lg_188, \
                         lg_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_257[k] = ab_x[k] * kg_185[k]
                       + lg_185[k];

            t_258[k] = ab_x[k] * kg_186[k]
                       + lg_186[k];

            t_259[k] = ab_x[k] * kg_187[k]
                       + lg_187[k];

            t_260[k] = ab_x[k] * kg_188[k]
                       + lg_188[k];

            t_261[k] = ab_x[k] * kg_189[k]
                       + lg_189[k];
        }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, ab_x, kg_190, kg_191, kg_192, \
                         kg_193, kg_194, lg_190, lg_191, lg_192, lg_193, \
                         lg_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_262[k] = ab_x[k] * kg_190[k]
                       + lg_190[k];

            t_263[k] = ab_x[k] * kg_191[k]
                       + lg_191[k];

            t_264[k] = ab_x[k] * kg_192[k]
                       + lg_192[k];

            t_265[k] = ab_x[k] * kg_193[k]
                       + lg_193[k];

            t_266[k] = ab_x[k] * kg_194[k]
                       + lg_194[k];
        }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, ab_y, kg_190, kg_191, kg_192, \
                         kg_193, kg_194, lg_265, lg_266, lg_267, lg_268, \
                         lg_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_267[k] = ab_y[k] * kg_190[k]
                       + lg_265[k];

            t_268[k] = ab_y[k] * kg_191[k]
                       + lg_266[k];

            t_269[k] = ab_y[k] * kg_192[k]
                       + lg_267[k];

            t_270[k] = ab_y[k] * kg_193[k]
                       + lg_268[k];

            t_271[k] = ab_y[k] * kg_194[k]
                       + lg_269[k];
        }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, ab_x, ab_z, kg_194, kg_195, kg_196, \
                         kg_197, lg_195, lg_196, lg_197, lg_284 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_272[k] = ab_z[k] * kg_194[k]
                       + lg_284[k];

            t_273[k] = ab_x[k] * kg_195[k]
                       + lg_195[k];

            t_274[k] = ab_x[k] * kg_196[k]
                       + lg_196[k];

            t_275[k] = ab_x[k] * kg_197[k]
                       + lg_197[k];
        }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, ab_x, kg_198, kg_199, kg_200, \
                         kg_201, kg_202, lg_198, lg_199, lg_200, lg_201, \
                         lg_202 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_276[k] = ab_x[k] * kg_198[k]
                       + lg_198[k];

            t_277[k] = ab_x[k] * kg_199[k]
                       + lg_199[k];

            t_278[k] = ab_x[k] * kg_200[k]
                       + lg_200[k];

            t_279[k] = ab_x[k] * kg_201[k]
                       + lg_201[k];

            t_280[k] = ab_x[k] * kg_202[k]
                       + lg_202[k];
        }
    }
}

static auto
compute_hrr_kh_out_of_first_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t kg, const size_t lg,
                                   const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_281 = buffer.data(target + 281 * ncomps + c);
        auto *t_282 = buffer.data(target + 282 * ncomps + c);
        auto *t_283 = buffer.data(target + 283 * ncomps + c);
        auto *t_284 = buffer.data(target + 284 * ncomps + c);
        auto *t_285 = buffer.data(target + 285 * ncomps + c);
        auto *t_286 = buffer.data(target + 286 * ncomps + c);
        auto *t_287 = buffer.data(target + 287 * ncomps + c);
        auto *t_288 = buffer.data(target + 288 * ncomps + c);
        auto *t_289 = buffer.data(target + 289 * ncomps + c);
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *kg_203 = buffer.data(kg + 203 * ncomps + c);
        const auto *kg_204 = buffer.data(kg + 204 * ncomps + c);
        const auto *kg_205 = buffer.data(kg + 205 * ncomps + c);
        const auto *kg_206 = buffer.data(kg + 206 * ncomps + c);
        const auto *kg_207 = buffer.data(kg + 207 * ncomps + c);
        const auto *kg_208 = buffer.data(kg + 208 * ncomps + c);
        const auto *kg_209 = buffer.data(kg + 209 * ncomps + c);
        const auto *kg_210 = buffer.data(kg + 210 * ncomps + c);
        const auto *kg_211 = buffer.data(kg + 211 * ncomps + c);
        const auto *kg_212 = buffer.data(kg + 212 * ncomps + c);
        const auto *kg_213 = buffer.data(kg + 213 * ncomps + c);
        const auto *kg_214 = buffer.data(kg + 214 * ncomps + c);
        const auto *kg_215 = buffer.data(kg + 215 * ncomps + c);
        const auto *kg_216 = buffer.data(kg + 216 * ncomps + c);
        const auto *kg_217 = buffer.data(kg + 217 * ncomps + c);
        const auto *kg_218 = buffer.data(kg + 218 * ncomps + c);
        const auto *kg_219 = buffer.data(kg + 219 * ncomps + c);
        const auto *kg_220 = buffer.data(kg + 220 * ncomps + c);
        const auto *kg_221 = buffer.data(kg + 221 * ncomps + c);
        const auto *kg_222 = buffer.data(kg + 222 * ncomps + c);
        const auto *kg_223 = buffer.data(kg + 223 * ncomps + c);
        const auto *kg_224 = buffer.data(kg + 224 * ncomps + c);
        const auto *kg_225 = buffer.data(kg + 225 * ncomps + c);
        const auto *kg_226 = buffer.data(kg + 226 * ncomps + c);
        const auto *kg_227 = buffer.data(kg + 227 * ncomps + c);
        const auto *kg_228 = buffer.data(kg + 228 * ncomps + c);
        const auto *kg_229 = buffer.data(kg + 229 * ncomps + c);
        const auto *kg_230 = buffer.data(kg + 230 * ncomps + c);
        const auto *kg_231 = buffer.data(kg + 231 * ncomps + c);
        const auto *kg_232 = buffer.data(kg + 232 * ncomps + c);
        const auto *kg_233 = buffer.data(kg + 233 * ncomps + c);
        const auto *kg_234 = buffer.data(kg + 234 * ncomps + c);
        const auto *kg_235 = buffer.data(kg + 235 * ncomps + c);
        const auto *kg_236 = buffer.data(kg + 236 * ncomps + c);
        const auto *kg_237 = buffer.data(kg + 237 * ncomps + c);
        const auto *kg_238 = buffer.data(kg + 238 * ncomps + c);
        const auto *kg_239 = buffer.data(kg + 239 * ncomps + c);
        const auto *kg_240 = buffer.data(kg + 240 * ncomps + c);
        const auto *kg_241 = buffer.data(kg + 241 * ncomps + c);
        const auto *kg_242 = buffer.data(kg + 242 * ncomps + c);
        const auto *kg_243 = buffer.data(kg + 243 * ncomps + c);
        const auto *kg_244 = buffer.data(kg + 244 * ncomps + c);
        const auto *kg_245 = buffer.data(kg + 245 * ncomps + c);
        const auto *kg_246 = buffer.data(kg + 246 * ncomps + c);
        const auto *kg_247 = buffer.data(kg + 247 * ncomps + c);
        const auto *kg_248 = buffer.data(kg + 248 * ncomps + c);
        const auto *kg_249 = buffer.data(kg + 249 * ncomps + c);
        const auto *kg_250 = buffer.data(kg + 250 * ncomps + c);
        const auto *kg_251 = buffer.data(kg + 251 * ncomps + c);
        const auto *kg_252 = buffer.data(kg + 252 * ncomps + c);
        const auto *kg_253 = buffer.data(kg + 253 * ncomps + c);
        const auto *kg_254 = buffer.data(kg + 254 * ncomps + c);
        const auto *kg_255 = buffer.data(kg + 255 * ncomps + c);
        const auto *kg_256 = buffer.data(kg + 256 * ncomps + c);
        const auto *kg_257 = buffer.data(kg + 257 * ncomps + c);
        const auto *kg_258 = buffer.data(kg + 258 * ncomps + c);
        const auto *kg_259 = buffer.data(kg + 259 * ncomps + c);
        const auto *kg_260 = buffer.data(kg + 260 * ncomps + c);
        const auto *kg_261 = buffer.data(kg + 261 * ncomps + c);
        const auto *kg_262 = buffer.data(kg + 262 * ncomps + c);
        const auto *kg_263 = buffer.data(kg + 263 * ncomps + c);
        const auto *kg_264 = buffer.data(kg + 264 * ncomps + c);
        const auto *kg_265 = buffer.data(kg + 265 * ncomps + c);
        const auto *kg_266 = buffer.data(kg + 266 * ncomps + c);
        const auto *kg_267 = buffer.data(kg + 267 * ncomps + c);
        const auto *kg_268 = buffer.data(kg + 268 * ncomps + c);
        const auto *kg_269 = buffer.data(kg + 269 * ncomps + c);
        const auto *kg_270 = buffer.data(kg + 270 * ncomps + c);
        const auto *kg_271 = buffer.data(kg + 271 * ncomps + c);
        const auto *kg_272 = buffer.data(kg + 272 * ncomps + c);
        const auto *kg_273 = buffer.data(kg + 273 * ncomps + c);
        const auto *kg_274 = buffer.data(kg + 274 * ncomps + c);
        const auto *kg_275 = buffer.data(kg + 275 * ncomps + c);
        const auto *kg_276 = buffer.data(kg + 276 * ncomps + c);
        const auto *kg_277 = buffer.data(kg + 277 * ncomps + c);
        const auto *kg_278 = buffer.data(kg + 278 * ncomps + c);
        const auto *kg_279 = buffer.data(kg + 279 * ncomps + c);
        const auto *kg_280 = buffer.data(kg + 280 * ncomps + c);
        const auto *kg_281 = buffer.data(kg + 281 * ncomps + c);
        const auto *kg_282 = buffer.data(kg + 282 * ncomps + c);
        const auto *kg_283 = buffer.data(kg + 283 * ncomps + c);
        const auto *kg_284 = buffer.data(kg + 284 * ncomps + c);
        const auto *kg_285 = buffer.data(kg + 285 * ncomps + c);
        const auto *kg_286 = buffer.data(kg + 286 * ncomps + c);
        const auto *kg_287 = buffer.data(kg + 287 * ncomps + c);
        const auto *kg_288 = buffer.data(kg + 288 * ncomps + c);
        const auto *kg_289 = buffer.data(kg + 289 * ncomps + c);
        const auto *kg_290 = buffer.data(kg + 290 * ncomps + c);
        const auto *kg_291 = buffer.data(kg + 291 * ncomps + c);
        const auto *kg_292 = buffer.data(kg + 292 * ncomps + c);
        const auto *kg_293 = buffer.data(kg + 293 * ncomps + c);
        const auto *kg_294 = buffer.data(kg + 294 * ncomps + c);
        const auto *kg_295 = buffer.data(kg + 295 * ncomps + c);
        const auto *kg_296 = buffer.data(kg + 296 * ncomps + c);
        const auto *kg_297 = buffer.data(kg + 297 * ncomps + c);
        const auto *kg_298 = buffer.data(kg + 298 * ncomps + c);
        const auto *kg_299 = buffer.data(kg + 299 * ncomps + c);
        const auto *kg_300 = buffer.data(kg + 300 * ncomps + c);
        const auto *kg_301 = buffer.data(kg + 301 * ncomps + c);
        const auto *kg_302 = buffer.data(kg + 302 * ncomps + c);
        const auto *kg_303 = buffer.data(kg + 303 * ncomps + c);
        const auto *kg_304 = buffer.data(kg + 304 * ncomps + c);

        const auto *lg_203 = buffer.data(lg + 203 * ncomps + c);
        const auto *lg_204 = buffer.data(lg + 204 * ncomps + c);
        const auto *lg_205 = buffer.data(lg + 205 * ncomps + c);
        const auto *lg_206 = buffer.data(lg + 206 * ncomps + c);
        const auto *lg_207 = buffer.data(lg + 207 * ncomps + c);
        const auto *lg_208 = buffer.data(lg + 208 * ncomps + c);
        const auto *lg_209 = buffer.data(lg + 209 * ncomps + c);
        const auto *lg_210 = buffer.data(lg + 210 * ncomps + c);
        const auto *lg_211 = buffer.data(lg + 211 * ncomps + c);
        const auto *lg_212 = buffer.data(lg + 212 * ncomps + c);
        const auto *lg_213 = buffer.data(lg + 213 * ncomps + c);
        const auto *lg_214 = buffer.data(lg + 214 * ncomps + c);
        const auto *lg_215 = buffer.data(lg + 215 * ncomps + c);
        const auto *lg_216 = buffer.data(lg + 216 * ncomps + c);
        const auto *lg_217 = buffer.data(lg + 217 * ncomps + c);
        const auto *lg_218 = buffer.data(lg + 218 * ncomps + c);
        const auto *lg_219 = buffer.data(lg + 219 * ncomps + c);
        const auto *lg_220 = buffer.data(lg + 220 * ncomps + c);
        const auto *lg_221 = buffer.data(lg + 221 * ncomps + c);
        const auto *lg_222 = buffer.data(lg + 222 * ncomps + c);
        const auto *lg_223 = buffer.data(lg + 223 * ncomps + c);
        const auto *lg_224 = buffer.data(lg + 224 * ncomps + c);
        const auto *lg_225 = buffer.data(lg + 225 * ncomps + c);
        const auto *lg_226 = buffer.data(lg + 226 * ncomps + c);
        const auto *lg_227 = buffer.data(lg + 227 * ncomps + c);
        const auto *lg_228 = buffer.data(lg + 228 * ncomps + c);
        const auto *lg_229 = buffer.data(lg + 229 * ncomps + c);
        const auto *lg_230 = buffer.data(lg + 230 * ncomps + c);
        const auto *lg_231 = buffer.data(lg + 231 * ncomps + c);
        const auto *lg_232 = buffer.data(lg + 232 * ncomps + c);
        const auto *lg_233 = buffer.data(lg + 233 * ncomps + c);
        const auto *lg_234 = buffer.data(lg + 234 * ncomps + c);
        const auto *lg_235 = buffer.data(lg + 235 * ncomps + c);
        const auto *lg_236 = buffer.data(lg + 236 * ncomps + c);
        const auto *lg_237 = buffer.data(lg + 237 * ncomps + c);
        const auto *lg_238 = buffer.data(lg + 238 * ncomps + c);
        const auto *lg_239 = buffer.data(lg + 239 * ncomps + c);
        const auto *lg_240 = buffer.data(lg + 240 * ncomps + c);
        const auto *lg_241 = buffer.data(lg + 241 * ncomps + c);
        const auto *lg_242 = buffer.data(lg + 242 * ncomps + c);
        const auto *lg_243 = buffer.data(lg + 243 * ncomps + c);
        const auto *lg_244 = buffer.data(lg + 244 * ncomps + c);
        const auto *lg_245 = buffer.data(lg + 245 * ncomps + c);
        const auto *lg_246 = buffer.data(lg + 246 * ncomps + c);
        const auto *lg_247 = buffer.data(lg + 247 * ncomps + c);
        const auto *lg_248 = buffer.data(lg + 248 * ncomps + c);
        const auto *lg_249 = buffer.data(lg + 249 * ncomps + c);
        const auto *lg_250 = buffer.data(lg + 250 * ncomps + c);
        const auto *lg_251 = buffer.data(lg + 251 * ncomps + c);
        const auto *lg_252 = buffer.data(lg + 252 * ncomps + c);
        const auto *lg_253 = buffer.data(lg + 253 * ncomps + c);
        const auto *lg_254 = buffer.data(lg + 254 * ncomps + c);
        const auto *lg_255 = buffer.data(lg + 255 * ncomps + c);
        const auto *lg_256 = buffer.data(lg + 256 * ncomps + c);
        const auto *lg_257 = buffer.data(lg + 257 * ncomps + c);
        const auto *lg_258 = buffer.data(lg + 258 * ncomps + c);
        const auto *lg_259 = buffer.data(lg + 259 * ncomps + c);
        const auto *lg_260 = buffer.data(lg + 260 * ncomps + c);
        const auto *lg_261 = buffer.data(lg + 261 * ncomps + c);
        const auto *lg_262 = buffer.data(lg + 262 * ncomps + c);
        const auto *lg_263 = buffer.data(lg + 263 * ncomps + c);
        const auto *lg_264 = buffer.data(lg + 264 * ncomps + c);
        const auto *lg_265 = buffer.data(lg + 265 * ncomps + c);
        const auto *lg_266 = buffer.data(lg + 266 * ncomps + c);
        const auto *lg_267 = buffer.data(lg + 267 * ncomps + c);
        const auto *lg_268 = buffer.data(lg + 268 * ncomps + c);
        const auto *lg_269 = buffer.data(lg + 269 * ncomps + c);
        const auto *lg_270 = buffer.data(lg + 270 * ncomps + c);
        const auto *lg_271 = buffer.data(lg + 271 * ncomps + c);
        const auto *lg_272 = buffer.data(lg + 272 * ncomps + c);
        const auto *lg_273 = buffer.data(lg + 273 * ncomps + c);
        const auto *lg_274 = buffer.data(lg + 274 * ncomps + c);
        const auto *lg_275 = buffer.data(lg + 275 * ncomps + c);
        const auto *lg_276 = buffer.data(lg + 276 * ncomps + c);
        const auto *lg_277 = buffer.data(lg + 277 * ncomps + c);
        const auto *lg_278 = buffer.data(lg + 278 * ncomps + c);
        const auto *lg_279 = buffer.data(lg + 279 * ncomps + c);
        const auto *lg_280 = buffer.data(lg + 280 * ncomps + c);
        const auto *lg_281 = buffer.data(lg + 281 * ncomps + c);
        const auto *lg_282 = buffer.data(lg + 282 * ncomps + c);
        const auto *lg_283 = buffer.data(lg + 283 * ncomps + c);
        const auto *lg_284 = buffer.data(lg + 284 * ncomps + c);
        const auto *lg_285 = buffer.data(lg + 285 * ncomps + c);
        const auto *lg_286 = buffer.data(lg + 286 * ncomps + c);
        const auto *lg_287 = buffer.data(lg + 287 * ncomps + c);
        const auto *lg_288 = buffer.data(lg + 288 * ncomps + c);
        const auto *lg_289 = buffer.data(lg + 289 * ncomps + c);
        const auto *lg_290 = buffer.data(lg + 290 * ncomps + c);
        const auto *lg_291 = buffer.data(lg + 291 * ncomps + c);
        const auto *lg_292 = buffer.data(lg + 292 * ncomps + c);
        const auto *lg_293 = buffer.data(lg + 293 * ncomps + c);
        const auto *lg_294 = buffer.data(lg + 294 * ncomps + c);
        const auto *lg_295 = buffer.data(lg + 295 * ncomps + c);
        const auto *lg_296 = buffer.data(lg + 296 * ncomps + c);
        const auto *lg_297 = buffer.data(lg + 297 * ncomps + c);
        const auto *lg_298 = buffer.data(lg + 298 * ncomps + c);
        const auto *lg_299 = buffer.data(lg + 299 * ncomps + c);
        const auto *lg_300 = buffer.data(lg + 300 * ncomps + c);
        const auto *lg_301 = buffer.data(lg + 301 * ncomps + c);
        const auto *lg_302 = buffer.data(lg + 302 * ncomps + c);
        const auto *lg_303 = buffer.data(lg + 303 * ncomps + c);
        const auto *lg_304 = buffer.data(lg + 304 * ncomps + c);
        const auto *lg_314 = buffer.data(lg + 314 * ncomps + c);
        const auto *lg_325 = buffer.data(lg + 325 * ncomps + c);
        const auto *lg_326 = buffer.data(lg + 326 * ncomps + c);
        const auto *lg_327 = buffer.data(lg + 327 * ncomps + c);
        const auto *lg_328 = buffer.data(lg + 328 * ncomps + c);
        const auto *lg_329 = buffer.data(lg + 329 * ncomps + c);
        const auto *lg_340 = buffer.data(lg + 340 * ncomps + c);
        const auto *lg_341 = buffer.data(lg + 341 * ncomps + c);
        const auto *lg_342 = buffer.data(lg + 342 * ncomps + c);
        const auto *lg_343 = buffer.data(lg + 343 * ncomps + c);
        const auto *lg_344 = buffer.data(lg + 344 * ncomps + c);
        const auto *lg_355 = buffer.data(lg + 355 * ncomps + c);
        const auto *lg_356 = buffer.data(lg + 356 * ncomps + c);
        const auto *lg_357 = buffer.data(lg + 357 * ncomps + c);
        const auto *lg_358 = buffer.data(lg + 358 * ncomps + c);
        const auto *lg_359 = buffer.data(lg + 359 * ncomps + c);
        const auto *lg_370 = buffer.data(lg + 370 * ncomps + c);
        const auto *lg_371 = buffer.data(lg + 371 * ncomps + c);
        const auto *lg_372 = buffer.data(lg + 372 * ncomps + c);
        const auto *lg_373 = buffer.data(lg + 373 * ncomps + c);
        const auto *lg_374 = buffer.data(lg + 374 * ncomps + c);
        const auto *lg_385 = buffer.data(lg + 385 * ncomps + c);
        const auto *lg_386 = buffer.data(lg + 386 * ncomps + c);
        const auto *lg_387 = buffer.data(lg + 387 * ncomps + c);
        const auto *lg_388 = buffer.data(lg + 388 * ncomps + c);
        const auto *lg_389 = buffer.data(lg + 389 * ncomps + c);
        const auto *lg_404 = buffer.data(lg + 404 * ncomps + c);

#pragma omp simd aligned(t_281, t_282, t_283, t_284, t_285, ab_x, kg_203, kg_204, kg_205, \
                         kg_206, kg_207, lg_203, lg_204, lg_205, lg_206, \
                         lg_207 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_281[k] = ab_x[k] * kg_203[k]
                       + lg_203[k];

            t_282[k] = ab_x[k] * kg_204[k]
                       + lg_204[k];

            t_283[k] = ab_x[k] * kg_205[k]
                       + lg_205[k];

            t_284[k] = ab_x[k] * kg_206[k]
                       + lg_206[k];

            t_285[k] = ab_x[k] * kg_207[k]
                       + lg_207[k];
        }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, ab_x, ab_y, kg_205, kg_206, kg_208, \
                         kg_209, lg_208, lg_209, lg_280, lg_281 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_286[k] = ab_x[k] * kg_208[k]
                       + lg_208[k];

            t_287[k] = ab_x[k] * kg_209[k]
                       + lg_209[k];

            t_288[k] = ab_y[k] * kg_205[k]
                       + lg_280[k];

            t_289[k] = ab_y[k] * kg_206[k]
                       + lg_281[k];
        }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, ab_y, ab_z, kg_207, kg_208, kg_209, \
                         lg_282, lg_283, lg_284, lg_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = ab_y[k] * kg_207[k]
                       + lg_282[k];

            t_291[k] = ab_y[k] * kg_208[k]
                       + lg_283[k];

            t_292[k] = ab_y[k] * kg_209[k]
                       + lg_284[k];

            t_293[k] = ab_z[k] * kg_209[k]
                       + lg_299[k];
        }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, ab_x, kg_210, kg_211, kg_212, \
                         kg_213, kg_214, lg_210, lg_211, lg_212, lg_213, \
                         lg_214 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_294[k] = ab_x[k] * kg_210[k]
                       + lg_210[k];

            t_295[k] = ab_x[k] * kg_211[k]
                       + lg_211[k];

            t_296[k] = ab_x[k] * kg_212[k]
                       + lg_212[k];

            t_297[k] = ab_x[k] * kg_213[k]
                       + lg_213[k];

            t_298[k] = ab_x[k] * kg_214[k]
                       + lg_214[k];
        }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, ab_x, kg_215, kg_216, kg_217, \
                         kg_218, kg_219, lg_215, lg_216, lg_217, lg_218, \
                         lg_219 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_299[k] = ab_x[k] * kg_215[k]
                       + lg_215[k];

            t_300[k] = ab_x[k] * kg_216[k]
                       + lg_216[k];

            t_301[k] = ab_x[k] * kg_217[k]
                       + lg_217[k];

            t_302[k] = ab_x[k] * kg_218[k]
                       + lg_218[k];

            t_303[k] = ab_x[k] * kg_219[k]
                       + lg_219[k];
        }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, ab_x, kg_220, kg_221, kg_222, \
                         kg_223, kg_224, lg_220, lg_221, lg_222, lg_223, \
                         lg_224 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_304[k] = ab_x[k] * kg_220[k]
                       + lg_220[k];

            t_305[k] = ab_x[k] * kg_221[k]
                       + lg_221[k];

            t_306[k] = ab_x[k] * kg_222[k]
                       + lg_222[k];

            t_307[k] = ab_x[k] * kg_223[k]
                       + lg_223[k];

            t_308[k] = ab_x[k] * kg_224[k]
                       + lg_224[k];
        }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, ab_y, kg_220, kg_221, kg_222, \
                         kg_223, kg_224, lg_295, lg_296, lg_297, lg_298, \
                         lg_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_309[k] = ab_y[k] * kg_220[k]
                       + lg_295[k];

            t_310[k] = ab_y[k] * kg_221[k]
                       + lg_296[k];

            t_311[k] = ab_y[k] * kg_222[k]
                       + lg_297[k];

            t_312[k] = ab_y[k] * kg_223[k]
                       + lg_298[k];

            t_313[k] = ab_y[k] * kg_224[k]
                       + lg_299[k];
        }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, ab_x, ab_z, kg_224, kg_225, kg_226, \
                         kg_227, lg_225, lg_226, lg_227, lg_314 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_314[k] = ab_z[k] * kg_224[k]
                       + lg_314[k];

            t_315[k] = ab_x[k] * kg_225[k]
                       + lg_225[k];

            t_316[k] = ab_x[k] * kg_226[k]
                       + lg_226[k];

            t_317[k] = ab_x[k] * kg_227[k]
                       + lg_227[k];
        }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, ab_x, kg_228, kg_229, kg_230, \
                         kg_231, kg_232, lg_228, lg_229, lg_230, lg_231, \
                         lg_232 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_318[k] = ab_x[k] * kg_228[k]
                       + lg_228[k];

            t_319[k] = ab_x[k] * kg_229[k]
                       + lg_229[k];

            t_320[k] = ab_x[k] * kg_230[k]
                       + lg_230[k];

            t_321[k] = ab_x[k] * kg_231[k]
                       + lg_231[k];

            t_322[k] = ab_x[k] * kg_232[k]
                       + lg_232[k];
        }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, ab_x, kg_233, kg_234, kg_235, \
                         kg_236, kg_237, lg_233, lg_234, lg_235, lg_236, \
                         lg_237 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_323[k] = ab_x[k] * kg_233[k]
                       + lg_233[k];

            t_324[k] = ab_x[k] * kg_234[k]
                       + lg_234[k];

            t_325[k] = ab_x[k] * kg_235[k]
                       + lg_235[k];

            t_326[k] = ab_x[k] * kg_236[k]
                       + lg_236[k];

            t_327[k] = ab_x[k] * kg_237[k]
                       + lg_237[k];
        }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, ab_x, ab_y, kg_235, kg_236, kg_238, \
                         kg_239, lg_238, lg_239, lg_325, lg_326 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_328[k] = ab_x[k] * kg_238[k]
                       + lg_238[k];

            t_329[k] = ab_x[k] * kg_239[k]
                       + lg_239[k];

            t_330[k] = ab_y[k] * kg_235[k]
                       + lg_325[k];

            t_331[k] = ab_y[k] * kg_236[k]
                       + lg_326[k];
        }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, ab_y, ab_z, kg_237, kg_238, kg_239, \
                         lg_327, lg_328, lg_329, lg_344 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_332[k] = ab_y[k] * kg_237[k]
                       + lg_327[k];

            t_333[k] = ab_y[k] * kg_238[k]
                       + lg_328[k];

            t_334[k] = ab_y[k] * kg_239[k]
                       + lg_329[k];

            t_335[k] = ab_z[k] * kg_239[k]
                       + lg_344[k];
        }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, t_340, ab_x, kg_240, kg_241, kg_242, \
                         kg_243, kg_244, lg_240, lg_241, lg_242, lg_243, \
                         lg_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_336[k] = ab_x[k] * kg_240[k]
                       + lg_240[k];

            t_337[k] = ab_x[k] * kg_241[k]
                       + lg_241[k];

            t_338[k] = ab_x[k] * kg_242[k]
                       + lg_242[k];

            t_339[k] = ab_x[k] * kg_243[k]
                       + lg_243[k];

            t_340[k] = ab_x[k] * kg_244[k]
                       + lg_244[k];
        }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, t_345, ab_x, kg_245, kg_246, kg_247, \
                         kg_248, kg_249, lg_245, lg_246, lg_247, lg_248, \
                         lg_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_341[k] = ab_x[k] * kg_245[k]
                       + lg_245[k];

            t_342[k] = ab_x[k] * kg_246[k]
                       + lg_246[k];

            t_343[k] = ab_x[k] * kg_247[k]
                       + lg_247[k];

            t_344[k] = ab_x[k] * kg_248[k]
                       + lg_248[k];

            t_345[k] = ab_x[k] * kg_249[k]
                       + lg_249[k];
        }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, ab_x, kg_250, kg_251, kg_252, \
                         kg_253, kg_254, lg_250, lg_251, lg_252, lg_253, \
                         lg_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_346[k] = ab_x[k] * kg_250[k]
                       + lg_250[k];

            t_347[k] = ab_x[k] * kg_251[k]
                       + lg_251[k];

            t_348[k] = ab_x[k] * kg_252[k]
                       + lg_252[k];

            t_349[k] = ab_x[k] * kg_253[k]
                       + lg_253[k];

            t_350[k] = ab_x[k] * kg_254[k]
                       + lg_254[k];
        }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, ab_y, kg_250, kg_251, kg_252, \
                         kg_253, kg_254, lg_340, lg_341, lg_342, lg_343, \
                         lg_344 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_351[k] = ab_y[k] * kg_250[k]
                       + lg_340[k];

            t_352[k] = ab_y[k] * kg_251[k]
                       + lg_341[k];

            t_353[k] = ab_y[k] * kg_252[k]
                       + lg_342[k];

            t_354[k] = ab_y[k] * kg_253[k]
                       + lg_343[k];

            t_355[k] = ab_y[k] * kg_254[k]
                       + lg_344[k];
        }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, ab_x, ab_z, kg_254, kg_255, kg_256, \
                         kg_257, lg_255, lg_256, lg_257, lg_359 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_356[k] = ab_z[k] * kg_254[k]
                       + lg_359[k];

            t_357[k] = ab_x[k] * kg_255[k]
                       + lg_255[k];

            t_358[k] = ab_x[k] * kg_256[k]
                       + lg_256[k];

            t_359[k] = ab_x[k] * kg_257[k]
                       + lg_257[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, kg_258, kg_259, kg_260, \
                         kg_261, kg_262, lg_258, lg_259, lg_260, lg_261, \
                         lg_262 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = ab_x[k] * kg_258[k]
                       + lg_258[k];

            t_361[k] = ab_x[k] * kg_259[k]
                       + lg_259[k];

            t_362[k] = ab_x[k] * kg_260[k]
                       + lg_260[k];

            t_363[k] = ab_x[k] * kg_261[k]
                       + lg_261[k];

            t_364[k] = ab_x[k] * kg_262[k]
                       + lg_262[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, kg_263, kg_264, kg_265, \
                         kg_266, kg_267, lg_263, lg_264, lg_265, lg_266, \
                         lg_267 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = ab_x[k] * kg_263[k]
                       + lg_263[k];

            t_366[k] = ab_x[k] * kg_264[k]
                       + lg_264[k];

            t_367[k] = ab_x[k] * kg_265[k]
                       + lg_265[k];

            t_368[k] = ab_x[k] * kg_266[k]
                       + lg_266[k];

            t_369[k] = ab_x[k] * kg_267[k]
                       + lg_267[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, ab_x, ab_y, kg_265, kg_266, kg_268, \
                         kg_269, lg_268, lg_269, lg_355, lg_356 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = ab_x[k] * kg_268[k]
                       + lg_268[k];

            t_371[k] = ab_x[k] * kg_269[k]
                       + lg_269[k];

            t_372[k] = ab_y[k] * kg_265[k]
                       + lg_355[k];

            t_373[k] = ab_y[k] * kg_266[k]
                       + lg_356[k];
        }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, ab_y, ab_z, kg_267, kg_268, kg_269, \
                         lg_357, lg_358, lg_359, lg_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_374[k] = ab_y[k] * kg_267[k]
                       + lg_357[k];

            t_375[k] = ab_y[k] * kg_268[k]
                       + lg_358[k];

            t_376[k] = ab_y[k] * kg_269[k]
                       + lg_359[k];

            t_377[k] = ab_z[k] * kg_269[k]
                       + lg_374[k];
        }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, t_382, ab_x, kg_270, kg_271, kg_272, \
                         kg_273, kg_274, lg_270, lg_271, lg_272, lg_273, \
                         lg_274 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_378[k] = ab_x[k] * kg_270[k]
                       + lg_270[k];

            t_379[k] = ab_x[k] * kg_271[k]
                       + lg_271[k];

            t_380[k] = ab_x[k] * kg_272[k]
                       + lg_272[k];

            t_381[k] = ab_x[k] * kg_273[k]
                       + lg_273[k];

            t_382[k] = ab_x[k] * kg_274[k]
                       + lg_274[k];
        }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, ab_x, kg_275, kg_276, kg_277, \
                         kg_278, kg_279, lg_275, lg_276, lg_277, lg_278, \
                         lg_279 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_383[k] = ab_x[k] * kg_275[k]
                       + lg_275[k];

            t_384[k] = ab_x[k] * kg_276[k]
                       + lg_276[k];

            t_385[k] = ab_x[k] * kg_277[k]
                       + lg_277[k];

            t_386[k] = ab_x[k] * kg_278[k]
                       + lg_278[k];

            t_387[k] = ab_x[k] * kg_279[k]
                       + lg_279[k];
        }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, t_392, ab_x, kg_280, kg_281, kg_282, \
                         kg_283, kg_284, lg_280, lg_281, lg_282, lg_283, \
                         lg_284 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_388[k] = ab_x[k] * kg_280[k]
                       + lg_280[k];

            t_389[k] = ab_x[k] * kg_281[k]
                       + lg_281[k];

            t_390[k] = ab_x[k] * kg_282[k]
                       + lg_282[k];

            t_391[k] = ab_x[k] * kg_283[k]
                       + lg_283[k];

            t_392[k] = ab_x[k] * kg_284[k]
                       + lg_284[k];
        }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, t_397, ab_y, kg_280, kg_281, kg_282, \
                         kg_283, kg_284, lg_370, lg_371, lg_372, lg_373, \
                         lg_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_393[k] = ab_y[k] * kg_280[k]
                       + lg_370[k];

            t_394[k] = ab_y[k] * kg_281[k]
                       + lg_371[k];

            t_395[k] = ab_y[k] * kg_282[k]
                       + lg_372[k];

            t_396[k] = ab_y[k] * kg_283[k]
                       + lg_373[k];

            t_397[k] = ab_y[k] * kg_284[k]
                       + lg_374[k];
        }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, ab_x, ab_z, kg_284, kg_285, kg_286, \
                         kg_287, lg_285, lg_286, lg_287, lg_389 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_398[k] = ab_z[k] * kg_284[k]
                       + lg_389[k];

            t_399[k] = ab_x[k] * kg_285[k]
                       + lg_285[k];

            t_400[k] = ab_x[k] * kg_286[k]
                       + lg_286[k];

            t_401[k] = ab_x[k] * kg_287[k]
                       + lg_287[k];
        }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, ab_x, kg_288, kg_289, kg_290, \
                         kg_291, kg_292, lg_288, lg_289, lg_290, lg_291, \
                         lg_292 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_402[k] = ab_x[k] * kg_288[k]
                       + lg_288[k];

            t_403[k] = ab_x[k] * kg_289[k]
                       + lg_289[k];

            t_404[k] = ab_x[k] * kg_290[k]
                       + lg_290[k];

            t_405[k] = ab_x[k] * kg_291[k]
                       + lg_291[k];

            t_406[k] = ab_x[k] * kg_292[k]
                       + lg_292[k];
        }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, ab_x, kg_293, kg_294, kg_295, \
                         kg_296, kg_297, lg_293, lg_294, lg_295, lg_296, \
                         lg_297 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_407[k] = ab_x[k] * kg_293[k]
                       + lg_293[k];

            t_408[k] = ab_x[k] * kg_294[k]
                       + lg_294[k];

            t_409[k] = ab_x[k] * kg_295[k]
                       + lg_295[k];

            t_410[k] = ab_x[k] * kg_296[k]
                       + lg_296[k];

            t_411[k] = ab_x[k] * kg_297[k]
                       + lg_297[k];
        }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, ab_x, ab_y, kg_295, kg_296, kg_298, \
                         kg_299, lg_298, lg_299, lg_385, lg_386 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_412[k] = ab_x[k] * kg_298[k]
                       + lg_298[k];

            t_413[k] = ab_x[k] * kg_299[k]
                       + lg_299[k];

            t_414[k] = ab_y[k] * kg_295[k]
                       + lg_385[k];

            t_415[k] = ab_y[k] * kg_296[k]
                       + lg_386[k];
        }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, ab_y, ab_z, kg_297, kg_298, kg_299, \
                         lg_387, lg_388, lg_389, lg_404 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_416[k] = ab_y[k] * kg_297[k]
                       + lg_387[k];

            t_417[k] = ab_y[k] * kg_298[k]
                       + lg_388[k];

            t_418[k] = ab_y[k] * kg_299[k]
                       + lg_389[k];

            t_419[k] = ab_z[k] * kg_299[k]
                       + lg_404[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_x, kg_300, kg_301, kg_302, \
                         kg_303, kg_304, lg_300, lg_301, lg_302, lg_303, \
                         lg_304 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = ab_x[k] * kg_300[k]
                       + lg_300[k];

            t_421[k] = ab_x[k] * kg_301[k]
                       + lg_301[k];

            t_422[k] = ab_x[k] * kg_302[k]
                       + lg_302[k];

            t_423[k] = ab_x[k] * kg_303[k]
                       + lg_303[k];

            t_424[k] = ab_x[k] * kg_304[k]
                       + lg_304[k];
        }
    }
}

static auto
compute_hrr_kh_out_of_first_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t kg, const size_t lg,
                                   const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *kg_305 = buffer.data(kg + 305 * ncomps + c);
        const auto *kg_306 = buffer.data(kg + 306 * ncomps + c);
        const auto *kg_307 = buffer.data(kg + 307 * ncomps + c);
        const auto *kg_308 = buffer.data(kg + 308 * ncomps + c);
        const auto *kg_309 = buffer.data(kg + 309 * ncomps + c);
        const auto *kg_310 = buffer.data(kg + 310 * ncomps + c);
        const auto *kg_311 = buffer.data(kg + 311 * ncomps + c);
        const auto *kg_312 = buffer.data(kg + 312 * ncomps + c);
        const auto *kg_313 = buffer.data(kg + 313 * ncomps + c);
        const auto *kg_314 = buffer.data(kg + 314 * ncomps + c);
        const auto *kg_315 = buffer.data(kg + 315 * ncomps + c);
        const auto *kg_316 = buffer.data(kg + 316 * ncomps + c);
        const auto *kg_317 = buffer.data(kg + 317 * ncomps + c);
        const auto *kg_318 = buffer.data(kg + 318 * ncomps + c);
        const auto *kg_319 = buffer.data(kg + 319 * ncomps + c);
        const auto *kg_320 = buffer.data(kg + 320 * ncomps + c);
        const auto *kg_321 = buffer.data(kg + 321 * ncomps + c);
        const auto *kg_322 = buffer.data(kg + 322 * ncomps + c);
        const auto *kg_323 = buffer.data(kg + 323 * ncomps + c);
        const auto *kg_324 = buffer.data(kg + 324 * ncomps + c);
        const auto *kg_325 = buffer.data(kg + 325 * ncomps + c);
        const auto *kg_326 = buffer.data(kg + 326 * ncomps + c);
        const auto *kg_327 = buffer.data(kg + 327 * ncomps + c);
        const auto *kg_328 = buffer.data(kg + 328 * ncomps + c);
        const auto *kg_329 = buffer.data(kg + 329 * ncomps + c);
        const auto *kg_330 = buffer.data(kg + 330 * ncomps + c);
        const auto *kg_331 = buffer.data(kg + 331 * ncomps + c);
        const auto *kg_332 = buffer.data(kg + 332 * ncomps + c);
        const auto *kg_333 = buffer.data(kg + 333 * ncomps + c);
        const auto *kg_334 = buffer.data(kg + 334 * ncomps + c);
        const auto *kg_335 = buffer.data(kg + 335 * ncomps + c);
        const auto *kg_336 = buffer.data(kg + 336 * ncomps + c);
        const auto *kg_337 = buffer.data(kg + 337 * ncomps + c);
        const auto *kg_338 = buffer.data(kg + 338 * ncomps + c);
        const auto *kg_339 = buffer.data(kg + 339 * ncomps + c);
        const auto *kg_340 = buffer.data(kg + 340 * ncomps + c);
        const auto *kg_341 = buffer.data(kg + 341 * ncomps + c);
        const auto *kg_342 = buffer.data(kg + 342 * ncomps + c);
        const auto *kg_343 = buffer.data(kg + 343 * ncomps + c);
        const auto *kg_344 = buffer.data(kg + 344 * ncomps + c);
        const auto *kg_345 = buffer.data(kg + 345 * ncomps + c);
        const auto *kg_346 = buffer.data(kg + 346 * ncomps + c);
        const auto *kg_347 = buffer.data(kg + 347 * ncomps + c);
        const auto *kg_348 = buffer.data(kg + 348 * ncomps + c);
        const auto *kg_349 = buffer.data(kg + 349 * ncomps + c);
        const auto *kg_350 = buffer.data(kg + 350 * ncomps + c);
        const auto *kg_351 = buffer.data(kg + 351 * ncomps + c);
        const auto *kg_352 = buffer.data(kg + 352 * ncomps + c);
        const auto *kg_353 = buffer.data(kg + 353 * ncomps + c);
        const auto *kg_354 = buffer.data(kg + 354 * ncomps + c);
        const auto *kg_355 = buffer.data(kg + 355 * ncomps + c);
        const auto *kg_356 = buffer.data(kg + 356 * ncomps + c);
        const auto *kg_357 = buffer.data(kg + 357 * ncomps + c);
        const auto *kg_358 = buffer.data(kg + 358 * ncomps + c);
        const auto *kg_359 = buffer.data(kg + 359 * ncomps + c);
        const auto *kg_360 = buffer.data(kg + 360 * ncomps + c);
        const auto *kg_361 = buffer.data(kg + 361 * ncomps + c);
        const auto *kg_362 = buffer.data(kg + 362 * ncomps + c);
        const auto *kg_363 = buffer.data(kg + 363 * ncomps + c);
        const auto *kg_364 = buffer.data(kg + 364 * ncomps + c);
        const auto *kg_365 = buffer.data(kg + 365 * ncomps + c);
        const auto *kg_366 = buffer.data(kg + 366 * ncomps + c);
        const auto *kg_367 = buffer.data(kg + 367 * ncomps + c);
        const auto *kg_368 = buffer.data(kg + 368 * ncomps + c);
        const auto *kg_369 = buffer.data(kg + 369 * ncomps + c);
        const auto *kg_370 = buffer.data(kg + 370 * ncomps + c);
        const auto *kg_371 = buffer.data(kg + 371 * ncomps + c);
        const auto *kg_372 = buffer.data(kg + 372 * ncomps + c);
        const auto *kg_373 = buffer.data(kg + 373 * ncomps + c);
        const auto *kg_374 = buffer.data(kg + 374 * ncomps + c);
        const auto *kg_375 = buffer.data(kg + 375 * ncomps + c);
        const auto *kg_376 = buffer.data(kg + 376 * ncomps + c);
        const auto *kg_377 = buffer.data(kg + 377 * ncomps + c);
        const auto *kg_378 = buffer.data(kg + 378 * ncomps + c);
        const auto *kg_379 = buffer.data(kg + 379 * ncomps + c);
        const auto *kg_380 = buffer.data(kg + 380 * ncomps + c);
        const auto *kg_381 = buffer.data(kg + 381 * ncomps + c);
        const auto *kg_382 = buffer.data(kg + 382 * ncomps + c);
        const auto *kg_383 = buffer.data(kg + 383 * ncomps + c);
        const auto *kg_384 = buffer.data(kg + 384 * ncomps + c);
        const auto *kg_385 = buffer.data(kg + 385 * ncomps + c);
        const auto *kg_386 = buffer.data(kg + 386 * ncomps + c);
        const auto *kg_387 = buffer.data(kg + 387 * ncomps + c);
        const auto *kg_388 = buffer.data(kg + 388 * ncomps + c);
        const auto *kg_389 = buffer.data(kg + 389 * ncomps + c);
        const auto *kg_390 = buffer.data(kg + 390 * ncomps + c);
        const auto *kg_391 = buffer.data(kg + 391 * ncomps + c);
        const auto *kg_392 = buffer.data(kg + 392 * ncomps + c);
        const auto *kg_393 = buffer.data(kg + 393 * ncomps + c);
        const auto *kg_394 = buffer.data(kg + 394 * ncomps + c);
        const auto *kg_395 = buffer.data(kg + 395 * ncomps + c);
        const auto *kg_396 = buffer.data(kg + 396 * ncomps + c);
        const auto *kg_397 = buffer.data(kg + 397 * ncomps + c);
        const auto *kg_398 = buffer.data(kg + 398 * ncomps + c);
        const auto *kg_399 = buffer.data(kg + 399 * ncomps + c);
        const auto *kg_400 = buffer.data(kg + 400 * ncomps + c);
        const auto *kg_401 = buffer.data(kg + 401 * ncomps + c);
        const auto *kg_402 = buffer.data(kg + 402 * ncomps + c);
        const auto *kg_403 = buffer.data(kg + 403 * ncomps + c);
        const auto *kg_404 = buffer.data(kg + 404 * ncomps + c);

        const auto *lg_305 = buffer.data(lg + 305 * ncomps + c);
        const auto *lg_306 = buffer.data(lg + 306 * ncomps + c);
        const auto *lg_307 = buffer.data(lg + 307 * ncomps + c);
        const auto *lg_308 = buffer.data(lg + 308 * ncomps + c);
        const auto *lg_309 = buffer.data(lg + 309 * ncomps + c);
        const auto *lg_310 = buffer.data(lg + 310 * ncomps + c);
        const auto *lg_311 = buffer.data(lg + 311 * ncomps + c);
        const auto *lg_312 = buffer.data(lg + 312 * ncomps + c);
        const auto *lg_313 = buffer.data(lg + 313 * ncomps + c);
        const auto *lg_314 = buffer.data(lg + 314 * ncomps + c);
        const auto *lg_315 = buffer.data(lg + 315 * ncomps + c);
        const auto *lg_316 = buffer.data(lg + 316 * ncomps + c);
        const auto *lg_317 = buffer.data(lg + 317 * ncomps + c);
        const auto *lg_318 = buffer.data(lg + 318 * ncomps + c);
        const auto *lg_319 = buffer.data(lg + 319 * ncomps + c);
        const auto *lg_320 = buffer.data(lg + 320 * ncomps + c);
        const auto *lg_321 = buffer.data(lg + 321 * ncomps + c);
        const auto *lg_322 = buffer.data(lg + 322 * ncomps + c);
        const auto *lg_323 = buffer.data(lg + 323 * ncomps + c);
        const auto *lg_324 = buffer.data(lg + 324 * ncomps + c);
        const auto *lg_325 = buffer.data(lg + 325 * ncomps + c);
        const auto *lg_326 = buffer.data(lg + 326 * ncomps + c);
        const auto *lg_327 = buffer.data(lg + 327 * ncomps + c);
        const auto *lg_328 = buffer.data(lg + 328 * ncomps + c);
        const auto *lg_329 = buffer.data(lg + 329 * ncomps + c);
        const auto *lg_330 = buffer.data(lg + 330 * ncomps + c);
        const auto *lg_331 = buffer.data(lg + 331 * ncomps + c);
        const auto *lg_332 = buffer.data(lg + 332 * ncomps + c);
        const auto *lg_333 = buffer.data(lg + 333 * ncomps + c);
        const auto *lg_334 = buffer.data(lg + 334 * ncomps + c);
        const auto *lg_335 = buffer.data(lg + 335 * ncomps + c);
        const auto *lg_336 = buffer.data(lg + 336 * ncomps + c);
        const auto *lg_337 = buffer.data(lg + 337 * ncomps + c);
        const auto *lg_338 = buffer.data(lg + 338 * ncomps + c);
        const auto *lg_339 = buffer.data(lg + 339 * ncomps + c);
        const auto *lg_340 = buffer.data(lg + 340 * ncomps + c);
        const auto *lg_341 = buffer.data(lg + 341 * ncomps + c);
        const auto *lg_342 = buffer.data(lg + 342 * ncomps + c);
        const auto *lg_343 = buffer.data(lg + 343 * ncomps + c);
        const auto *lg_344 = buffer.data(lg + 344 * ncomps + c);
        const auto *lg_345 = buffer.data(lg + 345 * ncomps + c);
        const auto *lg_346 = buffer.data(lg + 346 * ncomps + c);
        const auto *lg_347 = buffer.data(lg + 347 * ncomps + c);
        const auto *lg_348 = buffer.data(lg + 348 * ncomps + c);
        const auto *lg_349 = buffer.data(lg + 349 * ncomps + c);
        const auto *lg_350 = buffer.data(lg + 350 * ncomps + c);
        const auto *lg_351 = buffer.data(lg + 351 * ncomps + c);
        const auto *lg_352 = buffer.data(lg + 352 * ncomps + c);
        const auto *lg_353 = buffer.data(lg + 353 * ncomps + c);
        const auto *lg_354 = buffer.data(lg + 354 * ncomps + c);
        const auto *lg_355 = buffer.data(lg + 355 * ncomps + c);
        const auto *lg_356 = buffer.data(lg + 356 * ncomps + c);
        const auto *lg_357 = buffer.data(lg + 357 * ncomps + c);
        const auto *lg_358 = buffer.data(lg + 358 * ncomps + c);
        const auto *lg_359 = buffer.data(lg + 359 * ncomps + c);
        const auto *lg_360 = buffer.data(lg + 360 * ncomps + c);
        const auto *lg_361 = buffer.data(lg + 361 * ncomps + c);
        const auto *lg_362 = buffer.data(lg + 362 * ncomps + c);
        const auto *lg_363 = buffer.data(lg + 363 * ncomps + c);
        const auto *lg_364 = buffer.data(lg + 364 * ncomps + c);
        const auto *lg_365 = buffer.data(lg + 365 * ncomps + c);
        const auto *lg_366 = buffer.data(lg + 366 * ncomps + c);
        const auto *lg_367 = buffer.data(lg + 367 * ncomps + c);
        const auto *lg_368 = buffer.data(lg + 368 * ncomps + c);
        const auto *lg_369 = buffer.data(lg + 369 * ncomps + c);
        const auto *lg_370 = buffer.data(lg + 370 * ncomps + c);
        const auto *lg_371 = buffer.data(lg + 371 * ncomps + c);
        const auto *lg_372 = buffer.data(lg + 372 * ncomps + c);
        const auto *lg_373 = buffer.data(lg + 373 * ncomps + c);
        const auto *lg_374 = buffer.data(lg + 374 * ncomps + c);
        const auto *lg_375 = buffer.data(lg + 375 * ncomps + c);
        const auto *lg_376 = buffer.data(lg + 376 * ncomps + c);
        const auto *lg_377 = buffer.data(lg + 377 * ncomps + c);
        const auto *lg_378 = buffer.data(lg + 378 * ncomps + c);
        const auto *lg_379 = buffer.data(lg + 379 * ncomps + c);
        const auto *lg_380 = buffer.data(lg + 380 * ncomps + c);
        const auto *lg_381 = buffer.data(lg + 381 * ncomps + c);
        const auto *lg_382 = buffer.data(lg + 382 * ncomps + c);
        const auto *lg_383 = buffer.data(lg + 383 * ncomps + c);
        const auto *lg_384 = buffer.data(lg + 384 * ncomps + c);
        const auto *lg_385 = buffer.data(lg + 385 * ncomps + c);
        const auto *lg_386 = buffer.data(lg + 386 * ncomps + c);
        const auto *lg_387 = buffer.data(lg + 387 * ncomps + c);
        const auto *lg_388 = buffer.data(lg + 388 * ncomps + c);
        const auto *lg_389 = buffer.data(lg + 389 * ncomps + c);
        const auto *lg_390 = buffer.data(lg + 390 * ncomps + c);
        const auto *lg_391 = buffer.data(lg + 391 * ncomps + c);
        const auto *lg_392 = buffer.data(lg + 392 * ncomps + c);
        const auto *lg_393 = buffer.data(lg + 393 * ncomps + c);
        const auto *lg_394 = buffer.data(lg + 394 * ncomps + c);
        const auto *lg_395 = buffer.data(lg + 395 * ncomps + c);
        const auto *lg_396 = buffer.data(lg + 396 * ncomps + c);
        const auto *lg_397 = buffer.data(lg + 397 * ncomps + c);
        const auto *lg_398 = buffer.data(lg + 398 * ncomps + c);
        const auto *lg_399 = buffer.data(lg + 399 * ncomps + c);
        const auto *lg_400 = buffer.data(lg + 400 * ncomps + c);
        const auto *lg_401 = buffer.data(lg + 401 * ncomps + c);
        const auto *lg_402 = buffer.data(lg + 402 * ncomps + c);
        const auto *lg_403 = buffer.data(lg + 403 * ncomps + c);
        const auto *lg_404 = buffer.data(lg + 404 * ncomps + c);
        const auto *lg_419 = buffer.data(lg + 419 * ncomps + c);
        const auto *lg_430 = buffer.data(lg + 430 * ncomps + c);
        const auto *lg_431 = buffer.data(lg + 431 * ncomps + c);
        const auto *lg_432 = buffer.data(lg + 432 * ncomps + c);
        const auto *lg_433 = buffer.data(lg + 433 * ncomps + c);
        const auto *lg_434 = buffer.data(lg + 434 * ncomps + c);
        const auto *lg_445 = buffer.data(lg + 445 * ncomps + c);
        const auto *lg_446 = buffer.data(lg + 446 * ncomps + c);
        const auto *lg_447 = buffer.data(lg + 447 * ncomps + c);
        const auto *lg_448 = buffer.data(lg + 448 * ncomps + c);
        const auto *lg_449 = buffer.data(lg + 449 * ncomps + c);
        const auto *lg_460 = buffer.data(lg + 460 * ncomps + c);
        const auto *lg_461 = buffer.data(lg + 461 * ncomps + c);
        const auto *lg_462 = buffer.data(lg + 462 * ncomps + c);
        const auto *lg_463 = buffer.data(lg + 463 * ncomps + c);
        const auto *lg_464 = buffer.data(lg + 464 * ncomps + c);
        const auto *lg_475 = buffer.data(lg + 475 * ncomps + c);
        const auto *lg_476 = buffer.data(lg + 476 * ncomps + c);
        const auto *lg_477 = buffer.data(lg + 477 * ncomps + c);
        const auto *lg_478 = buffer.data(lg + 478 * ncomps + c);
        const auto *lg_479 = buffer.data(lg + 479 * ncomps + c);
        const auto *lg_490 = buffer.data(lg + 490 * ncomps + c);
        const auto *lg_491 = buffer.data(lg + 491 * ncomps + c);
        const auto *lg_492 = buffer.data(lg + 492 * ncomps + c);
        const auto *lg_493 = buffer.data(lg + 493 * ncomps + c);
        const auto *lg_494 = buffer.data(lg + 494 * ncomps + c);
        const auto *lg_505 = buffer.data(lg + 505 * ncomps + c);
        const auto *lg_506 = buffer.data(lg + 506 * ncomps + c);
        const auto *lg_507 = buffer.data(lg + 507 * ncomps + c);
        const auto *lg_508 = buffer.data(lg + 508 * ncomps + c);
        const auto *lg_509 = buffer.data(lg + 509 * ncomps + c);

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_x, kg_305, kg_306, kg_307, \
                         kg_308, kg_309, lg_305, lg_306, lg_307, lg_308, \
                         lg_309 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = ab_x[k] * kg_305[k]
                       + lg_305[k];

            t_426[k] = ab_x[k] * kg_306[k]
                       + lg_306[k];

            t_427[k] = ab_x[k] * kg_307[k]
                       + lg_307[k];

            t_428[k] = ab_x[k] * kg_308[k]
                       + lg_308[k];

            t_429[k] = ab_x[k] * kg_309[k]
                       + lg_309[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_x, kg_310, kg_311, kg_312, \
                         kg_313, kg_314, lg_310, lg_311, lg_312, lg_313, \
                         lg_314 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = ab_x[k] * kg_310[k]
                       + lg_310[k];

            t_431[k] = ab_x[k] * kg_311[k]
                       + lg_311[k];

            t_432[k] = ab_x[k] * kg_312[k]
                       + lg_312[k];

            t_433[k] = ab_x[k] * kg_313[k]
                       + lg_313[k];

            t_434[k] = ab_x[k] * kg_314[k]
                       + lg_314[k];
        }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_y, kg_310, kg_311, kg_312, \
                         kg_313, kg_314, lg_400, lg_401, lg_402, lg_403, \
                         lg_404 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = ab_y[k] * kg_310[k]
                       + lg_400[k];

            t_436[k] = ab_y[k] * kg_311[k]
                       + lg_401[k];

            t_437[k] = ab_y[k] * kg_312[k]
                       + lg_402[k];

            t_438[k] = ab_y[k] * kg_313[k]
                       + lg_403[k];

            t_439[k] = ab_y[k] * kg_314[k]
                       + lg_404[k];
        }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, ab_x, ab_z, kg_314, kg_315, kg_316, \
                         kg_317, lg_315, lg_316, lg_317, lg_419 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = ab_z[k] * kg_314[k]
                       + lg_419[k];

            t_441[k] = ab_x[k] * kg_315[k]
                       + lg_315[k];

            t_442[k] = ab_x[k] * kg_316[k]
                       + lg_316[k];

            t_443[k] = ab_x[k] * kg_317[k]
                       + lg_317[k];
        }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, t_448, ab_x, kg_318, kg_319, kg_320, \
                         kg_321, kg_322, lg_318, lg_319, lg_320, lg_321, \
                         lg_322 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_444[k] = ab_x[k] * kg_318[k]
                       + lg_318[k];

            t_445[k] = ab_x[k] * kg_319[k]
                       + lg_319[k];

            t_446[k] = ab_x[k] * kg_320[k]
                       + lg_320[k];

            t_447[k] = ab_x[k] * kg_321[k]
                       + lg_321[k];

            t_448[k] = ab_x[k] * kg_322[k]
                       + lg_322[k];
        }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, t_453, ab_x, kg_323, kg_324, kg_325, \
                         kg_326, kg_327, lg_323, lg_324, lg_325, lg_326, \
                         lg_327 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_449[k] = ab_x[k] * kg_323[k]
                       + lg_323[k];

            t_450[k] = ab_x[k] * kg_324[k]
                       + lg_324[k];

            t_451[k] = ab_x[k] * kg_325[k]
                       + lg_325[k];

            t_452[k] = ab_x[k] * kg_326[k]
                       + lg_326[k];

            t_453[k] = ab_x[k] * kg_327[k]
                       + lg_327[k];
        }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, ab_x, ab_y, kg_325, kg_326, kg_328, \
                         kg_329, lg_328, lg_329, lg_430, lg_431 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_454[k] = ab_x[k] * kg_328[k]
                       + lg_328[k];

            t_455[k] = ab_x[k] * kg_329[k]
                       + lg_329[k];

            t_456[k] = ab_y[k] * kg_325[k]
                       + lg_430[k];

            t_457[k] = ab_y[k] * kg_326[k]
                       + lg_431[k];
        }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, ab_y, ab_z, kg_327, kg_328, kg_329, \
                         lg_432, lg_433, lg_434, lg_449 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_458[k] = ab_y[k] * kg_327[k]
                       + lg_432[k];

            t_459[k] = ab_y[k] * kg_328[k]
                       + lg_433[k];

            t_460[k] = ab_y[k] * kg_329[k]
                       + lg_434[k];

            t_461[k] = ab_z[k] * kg_329[k]
                       + lg_449[k];
        }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, t_466, ab_x, kg_330, kg_331, kg_332, \
                         kg_333, kg_334, lg_330, lg_331, lg_332, lg_333, \
                         lg_334 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_462[k] = ab_x[k] * kg_330[k]
                       + lg_330[k];

            t_463[k] = ab_x[k] * kg_331[k]
                       + lg_331[k];

            t_464[k] = ab_x[k] * kg_332[k]
                       + lg_332[k];

            t_465[k] = ab_x[k] * kg_333[k]
                       + lg_333[k];

            t_466[k] = ab_x[k] * kg_334[k]
                       + lg_334[k];
        }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, t_471, ab_x, kg_335, kg_336, kg_337, \
                         kg_338, kg_339, lg_335, lg_336, lg_337, lg_338, \
                         lg_339 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_467[k] = ab_x[k] * kg_335[k]
                       + lg_335[k];

            t_468[k] = ab_x[k] * kg_336[k]
                       + lg_336[k];

            t_469[k] = ab_x[k] * kg_337[k]
                       + lg_337[k];

            t_470[k] = ab_x[k] * kg_338[k]
                       + lg_338[k];

            t_471[k] = ab_x[k] * kg_339[k]
                       + lg_339[k];
        }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, t_476, ab_x, kg_340, kg_341, kg_342, \
                         kg_343, kg_344, lg_340, lg_341, lg_342, lg_343, \
                         lg_344 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_472[k] = ab_x[k] * kg_340[k]
                       + lg_340[k];

            t_473[k] = ab_x[k] * kg_341[k]
                       + lg_341[k];

            t_474[k] = ab_x[k] * kg_342[k]
                       + lg_342[k];

            t_475[k] = ab_x[k] * kg_343[k]
                       + lg_343[k];

            t_476[k] = ab_x[k] * kg_344[k]
                       + lg_344[k];
        }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, ab_y, kg_340, kg_341, kg_342, \
                         kg_343, kg_344, lg_445, lg_446, lg_447, lg_448, \
                         lg_449 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_477[k] = ab_y[k] * kg_340[k]
                       + lg_445[k];

            t_478[k] = ab_y[k] * kg_341[k]
                       + lg_446[k];

            t_479[k] = ab_y[k] * kg_342[k]
                       + lg_447[k];

            t_480[k] = ab_y[k] * kg_343[k]
                       + lg_448[k];

            t_481[k] = ab_y[k] * kg_344[k]
                       + lg_449[k];
        }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, ab_x, ab_z, kg_344, kg_345, kg_346, \
                         kg_347, lg_345, lg_346, lg_347, lg_464 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_482[k] = ab_z[k] * kg_344[k]
                       + lg_464[k];

            t_483[k] = ab_x[k] * kg_345[k]
                       + lg_345[k];

            t_484[k] = ab_x[k] * kg_346[k]
                       + lg_346[k];

            t_485[k] = ab_x[k] * kg_347[k]
                       + lg_347[k];
        }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, t_490, ab_x, kg_348, kg_349, kg_350, \
                         kg_351, kg_352, lg_348, lg_349, lg_350, lg_351, \
                         lg_352 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_486[k] = ab_x[k] * kg_348[k]
                       + lg_348[k];

            t_487[k] = ab_x[k] * kg_349[k]
                       + lg_349[k];

            t_488[k] = ab_x[k] * kg_350[k]
                       + lg_350[k];

            t_489[k] = ab_x[k] * kg_351[k]
                       + lg_351[k];

            t_490[k] = ab_x[k] * kg_352[k]
                       + lg_352[k];
        }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, ab_x, kg_353, kg_354, kg_355, \
                         kg_356, kg_357, lg_353, lg_354, lg_355, lg_356, \
                         lg_357 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_491[k] = ab_x[k] * kg_353[k]
                       + lg_353[k];

            t_492[k] = ab_x[k] * kg_354[k]
                       + lg_354[k];

            t_493[k] = ab_x[k] * kg_355[k]
                       + lg_355[k];

            t_494[k] = ab_x[k] * kg_356[k]
                       + lg_356[k];

            t_495[k] = ab_x[k] * kg_357[k]
                       + lg_357[k];
        }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, ab_x, ab_y, kg_355, kg_356, kg_358, \
                         kg_359, lg_358, lg_359, lg_460, lg_461 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_496[k] = ab_x[k] * kg_358[k]
                       + lg_358[k];

            t_497[k] = ab_x[k] * kg_359[k]
                       + lg_359[k];

            t_498[k] = ab_y[k] * kg_355[k]
                       + lg_460[k];

            t_499[k] = ab_y[k] * kg_356[k]
                       + lg_461[k];
        }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, ab_y, ab_z, kg_357, kg_358, kg_359, \
                         lg_462, lg_463, lg_464, lg_479 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_500[k] = ab_y[k] * kg_357[k]
                       + lg_462[k];

            t_501[k] = ab_y[k] * kg_358[k]
                       + lg_463[k];

            t_502[k] = ab_y[k] * kg_359[k]
                       + lg_464[k];

            t_503[k] = ab_z[k] * kg_359[k]
                       + lg_479[k];
        }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, ab_x, kg_360, kg_361, kg_362, \
                         kg_363, kg_364, lg_360, lg_361, lg_362, lg_363, \
                         lg_364 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_504[k] = ab_x[k] * kg_360[k]
                       + lg_360[k];

            t_505[k] = ab_x[k] * kg_361[k]
                       + lg_361[k];

            t_506[k] = ab_x[k] * kg_362[k]
                       + lg_362[k];

            t_507[k] = ab_x[k] * kg_363[k]
                       + lg_363[k];

            t_508[k] = ab_x[k] * kg_364[k]
                       + lg_364[k];
        }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, t_513, ab_x, kg_365, kg_366, kg_367, \
                         kg_368, kg_369, lg_365, lg_366, lg_367, lg_368, \
                         lg_369 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_509[k] = ab_x[k] * kg_365[k]
                       + lg_365[k];

            t_510[k] = ab_x[k] * kg_366[k]
                       + lg_366[k];

            t_511[k] = ab_x[k] * kg_367[k]
                       + lg_367[k];

            t_512[k] = ab_x[k] * kg_368[k]
                       + lg_368[k];

            t_513[k] = ab_x[k] * kg_369[k]
                       + lg_369[k];
        }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, t_518, ab_x, kg_370, kg_371, kg_372, \
                         kg_373, kg_374, lg_370, lg_371, lg_372, lg_373, \
                         lg_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_514[k] = ab_x[k] * kg_370[k]
                       + lg_370[k];

            t_515[k] = ab_x[k] * kg_371[k]
                       + lg_371[k];

            t_516[k] = ab_x[k] * kg_372[k]
                       + lg_372[k];

            t_517[k] = ab_x[k] * kg_373[k]
                       + lg_373[k];

            t_518[k] = ab_x[k] * kg_374[k]
                       + lg_374[k];
        }

#pragma omp simd aligned(t_519, t_520, t_521, t_522, t_523, ab_y, kg_370, kg_371, kg_372, \
                         kg_373, kg_374, lg_475, lg_476, lg_477, lg_478, \
                         lg_479 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_519[k] = ab_y[k] * kg_370[k]
                       + lg_475[k];

            t_520[k] = ab_y[k] * kg_371[k]
                       + lg_476[k];

            t_521[k] = ab_y[k] * kg_372[k]
                       + lg_477[k];

            t_522[k] = ab_y[k] * kg_373[k]
                       + lg_478[k];

            t_523[k] = ab_y[k] * kg_374[k]
                       + lg_479[k];
        }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, ab_x, ab_z, kg_374, kg_375, kg_376, \
                         kg_377, lg_375, lg_376, lg_377, lg_494 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_524[k] = ab_z[k] * kg_374[k]
                       + lg_494[k];

            t_525[k] = ab_x[k] * kg_375[k]
                       + lg_375[k];

            t_526[k] = ab_x[k] * kg_376[k]
                       + lg_376[k];

            t_527[k] = ab_x[k] * kg_377[k]
                       + lg_377[k];
        }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, t_532, ab_x, kg_378, kg_379, kg_380, \
                         kg_381, kg_382, lg_378, lg_379, lg_380, lg_381, \
                         lg_382 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_528[k] = ab_x[k] * kg_378[k]
                       + lg_378[k];

            t_529[k] = ab_x[k] * kg_379[k]
                       + lg_379[k];

            t_530[k] = ab_x[k] * kg_380[k]
                       + lg_380[k];

            t_531[k] = ab_x[k] * kg_381[k]
                       + lg_381[k];

            t_532[k] = ab_x[k] * kg_382[k]
                       + lg_382[k];
        }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, t_537, ab_x, kg_383, kg_384, kg_385, \
                         kg_386, kg_387, lg_383, lg_384, lg_385, lg_386, \
                         lg_387 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_533[k] = ab_x[k] * kg_383[k]
                       + lg_383[k];

            t_534[k] = ab_x[k] * kg_384[k]
                       + lg_384[k];

            t_535[k] = ab_x[k] * kg_385[k]
                       + lg_385[k];

            t_536[k] = ab_x[k] * kg_386[k]
                       + lg_386[k];

            t_537[k] = ab_x[k] * kg_387[k]
                       + lg_387[k];
        }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, ab_x, ab_y, kg_385, kg_386, kg_388, \
                         kg_389, lg_388, lg_389, lg_490, lg_491 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_538[k] = ab_x[k] * kg_388[k]
                       + lg_388[k];

            t_539[k] = ab_x[k] * kg_389[k]
                       + lg_389[k];

            t_540[k] = ab_y[k] * kg_385[k]
                       + lg_490[k];

            t_541[k] = ab_y[k] * kg_386[k]
                       + lg_491[k];
        }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, ab_y, ab_z, kg_387, kg_388, kg_389, \
                         lg_492, lg_493, lg_494, lg_509 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_542[k] = ab_y[k] * kg_387[k]
                       + lg_492[k];

            t_543[k] = ab_y[k] * kg_388[k]
                       + lg_493[k];

            t_544[k] = ab_y[k] * kg_389[k]
                       + lg_494[k];

            t_545[k] = ab_z[k] * kg_389[k]
                       + lg_509[k];
        }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, t_550, ab_x, kg_390, kg_391, kg_392, \
                         kg_393, kg_394, lg_390, lg_391, lg_392, lg_393, \
                         lg_394 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_546[k] = ab_x[k] * kg_390[k]
                       + lg_390[k];

            t_547[k] = ab_x[k] * kg_391[k]
                       + lg_391[k];

            t_548[k] = ab_x[k] * kg_392[k]
                       + lg_392[k];

            t_549[k] = ab_x[k] * kg_393[k]
                       + lg_393[k];

            t_550[k] = ab_x[k] * kg_394[k]
                       + lg_394[k];
        }

#pragma omp simd aligned(t_551, t_552, t_553, t_554, t_555, ab_x, kg_395, kg_396, kg_397, \
                         kg_398, kg_399, lg_395, lg_396, lg_397, lg_398, \
                         lg_399 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_551[k] = ab_x[k] * kg_395[k]
                       + lg_395[k];

            t_552[k] = ab_x[k] * kg_396[k]
                       + lg_396[k];

            t_553[k] = ab_x[k] * kg_397[k]
                       + lg_397[k];

            t_554[k] = ab_x[k] * kg_398[k]
                       + lg_398[k];

            t_555[k] = ab_x[k] * kg_399[k]
                       + lg_399[k];
        }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, ab_x, kg_400, kg_401, kg_402, \
                         kg_403, kg_404, lg_400, lg_401, lg_402, lg_403, \
                         lg_404 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_556[k] = ab_x[k] * kg_400[k]
                       + lg_400[k];

            t_557[k] = ab_x[k] * kg_401[k]
                       + lg_401[k];

            t_558[k] = ab_x[k] * kg_402[k]
                       + lg_402[k];

            t_559[k] = ab_x[k] * kg_403[k]
                       + lg_403[k];

            t_560[k] = ab_x[k] * kg_404[k]
                       + lg_404[k];
        }

#pragma omp simd aligned(t_561, t_562, t_563, t_564, t_565, ab_y, kg_400, kg_401, kg_402, \
                         kg_403, kg_404, lg_505, lg_506, lg_507, lg_508, \
                         lg_509 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_561[k] = ab_y[k] * kg_400[k]
                       + lg_505[k];

            t_562[k] = ab_y[k] * kg_401[k]
                       + lg_506[k];

            t_563[k] = ab_y[k] * kg_402[k]
                       + lg_507[k];

            t_564[k] = ab_y[k] * kg_403[k]
                       + lg_508[k];

            t_565[k] = ab_y[k] * kg_404[k]
                       + lg_509[k];
        }
    }
}

static auto
compute_hrr_kh_out_of_first_piece4(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t kg, const size_t lg,
                                   const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
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
        auto *t_675 = buffer.data(target + 675 * ncomps + c);
        auto *t_676 = buffer.data(target + 676 * ncomps + c);
        auto *t_677 = buffer.data(target + 677 * ncomps + c);
        auto *t_678 = buffer.data(target + 678 * ncomps + c);
        auto *t_679 = buffer.data(target + 679 * ncomps + c);
        auto *t_680 = buffer.data(target + 680 * ncomps + c);
        auto *t_681 = buffer.data(target + 681 * ncomps + c);
        auto *t_682 = buffer.data(target + 682 * ncomps + c);
        auto *t_683 = buffer.data(target + 683 * ncomps + c);
        auto *t_684 = buffer.data(target + 684 * ncomps + c);
        auto *t_685 = buffer.data(target + 685 * ncomps + c);
        auto *t_686 = buffer.data(target + 686 * ncomps + c);
        auto *t_687 = buffer.data(target + 687 * ncomps + c);
        auto *t_688 = buffer.data(target + 688 * ncomps + c);
        auto *t_689 = buffer.data(target + 689 * ncomps + c);
        auto *t_690 = buffer.data(target + 690 * ncomps + c);
        auto *t_691 = buffer.data(target + 691 * ncomps + c);
        auto *t_692 = buffer.data(target + 692 * ncomps + c);
        auto *t_693 = buffer.data(target + 693 * ncomps + c);
        auto *t_694 = buffer.data(target + 694 * ncomps + c);
        auto *t_695 = buffer.data(target + 695 * ncomps + c);
        auto *t_696 = buffer.data(target + 696 * ncomps + c);
        auto *t_697 = buffer.data(target + 697 * ncomps + c);
        auto *t_698 = buffer.data(target + 698 * ncomps + c);
        auto *t_699 = buffer.data(target + 699 * ncomps + c);
        auto *t_700 = buffer.data(target + 700 * ncomps + c);
        auto *t_701 = buffer.data(target + 701 * ncomps + c);
        auto *t_702 = buffer.data(target + 702 * ncomps + c);
        auto *t_703 = buffer.data(target + 703 * ncomps + c);
        auto *t_704 = buffer.data(target + 704 * ncomps + c);
        auto *t_705 = buffer.data(target + 705 * ncomps + c);
        auto *t_706 = buffer.data(target + 706 * ncomps + c);
        auto *t_707 = buffer.data(target + 707 * ncomps + c);
        auto *t_708 = buffer.data(target + 708 * ncomps + c);
        auto *t_709 = buffer.data(target + 709 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *kg_404 = buffer.data(kg + 404 * ncomps + c);
        const auto *kg_405 = buffer.data(kg + 405 * ncomps + c);
        const auto *kg_406 = buffer.data(kg + 406 * ncomps + c);
        const auto *kg_407 = buffer.data(kg + 407 * ncomps + c);
        const auto *kg_408 = buffer.data(kg + 408 * ncomps + c);
        const auto *kg_409 = buffer.data(kg + 409 * ncomps + c);
        const auto *kg_410 = buffer.data(kg + 410 * ncomps + c);
        const auto *kg_411 = buffer.data(kg + 411 * ncomps + c);
        const auto *kg_412 = buffer.data(kg + 412 * ncomps + c);
        const auto *kg_413 = buffer.data(kg + 413 * ncomps + c);
        const auto *kg_414 = buffer.data(kg + 414 * ncomps + c);
        const auto *kg_415 = buffer.data(kg + 415 * ncomps + c);
        const auto *kg_416 = buffer.data(kg + 416 * ncomps + c);
        const auto *kg_417 = buffer.data(kg + 417 * ncomps + c);
        const auto *kg_418 = buffer.data(kg + 418 * ncomps + c);
        const auto *kg_419 = buffer.data(kg + 419 * ncomps + c);
        const auto *kg_420 = buffer.data(kg + 420 * ncomps + c);
        const auto *kg_421 = buffer.data(kg + 421 * ncomps + c);
        const auto *kg_422 = buffer.data(kg + 422 * ncomps + c);
        const auto *kg_423 = buffer.data(kg + 423 * ncomps + c);
        const auto *kg_424 = buffer.data(kg + 424 * ncomps + c);
        const auto *kg_425 = buffer.data(kg + 425 * ncomps + c);
        const auto *kg_426 = buffer.data(kg + 426 * ncomps + c);
        const auto *kg_427 = buffer.data(kg + 427 * ncomps + c);
        const auto *kg_428 = buffer.data(kg + 428 * ncomps + c);
        const auto *kg_429 = buffer.data(kg + 429 * ncomps + c);
        const auto *kg_430 = buffer.data(kg + 430 * ncomps + c);
        const auto *kg_431 = buffer.data(kg + 431 * ncomps + c);
        const auto *kg_432 = buffer.data(kg + 432 * ncomps + c);
        const auto *kg_433 = buffer.data(kg + 433 * ncomps + c);
        const auto *kg_434 = buffer.data(kg + 434 * ncomps + c);
        const auto *kg_435 = buffer.data(kg + 435 * ncomps + c);
        const auto *kg_436 = buffer.data(kg + 436 * ncomps + c);
        const auto *kg_437 = buffer.data(kg + 437 * ncomps + c);
        const auto *kg_438 = buffer.data(kg + 438 * ncomps + c);
        const auto *kg_439 = buffer.data(kg + 439 * ncomps + c);
        const auto *kg_440 = buffer.data(kg + 440 * ncomps + c);
        const auto *kg_441 = buffer.data(kg + 441 * ncomps + c);
        const auto *kg_442 = buffer.data(kg + 442 * ncomps + c);
        const auto *kg_443 = buffer.data(kg + 443 * ncomps + c);
        const auto *kg_444 = buffer.data(kg + 444 * ncomps + c);
        const auto *kg_445 = buffer.data(kg + 445 * ncomps + c);
        const auto *kg_446 = buffer.data(kg + 446 * ncomps + c);
        const auto *kg_447 = buffer.data(kg + 447 * ncomps + c);
        const auto *kg_448 = buffer.data(kg + 448 * ncomps + c);
        const auto *kg_449 = buffer.data(kg + 449 * ncomps + c);
        const auto *kg_450 = buffer.data(kg + 450 * ncomps + c);
        const auto *kg_451 = buffer.data(kg + 451 * ncomps + c);
        const auto *kg_452 = buffer.data(kg + 452 * ncomps + c);
        const auto *kg_453 = buffer.data(kg + 453 * ncomps + c);
        const auto *kg_454 = buffer.data(kg + 454 * ncomps + c);
        const auto *kg_455 = buffer.data(kg + 455 * ncomps + c);
        const auto *kg_456 = buffer.data(kg + 456 * ncomps + c);
        const auto *kg_457 = buffer.data(kg + 457 * ncomps + c);
        const auto *kg_458 = buffer.data(kg + 458 * ncomps + c);
        const auto *kg_459 = buffer.data(kg + 459 * ncomps + c);
        const auto *kg_460 = buffer.data(kg + 460 * ncomps + c);
        const auto *kg_461 = buffer.data(kg + 461 * ncomps + c);
        const auto *kg_462 = buffer.data(kg + 462 * ncomps + c);
        const auto *kg_463 = buffer.data(kg + 463 * ncomps + c);
        const auto *kg_464 = buffer.data(kg + 464 * ncomps + c);
        const auto *kg_465 = buffer.data(kg + 465 * ncomps + c);
        const auto *kg_466 = buffer.data(kg + 466 * ncomps + c);
        const auto *kg_467 = buffer.data(kg + 467 * ncomps + c);
        const auto *kg_468 = buffer.data(kg + 468 * ncomps + c);
        const auto *kg_469 = buffer.data(kg + 469 * ncomps + c);
        const auto *kg_470 = buffer.data(kg + 470 * ncomps + c);
        const auto *kg_471 = buffer.data(kg + 471 * ncomps + c);
        const auto *kg_472 = buffer.data(kg + 472 * ncomps + c);
        const auto *kg_473 = buffer.data(kg + 473 * ncomps + c);
        const auto *kg_474 = buffer.data(kg + 474 * ncomps + c);
        const auto *kg_475 = buffer.data(kg + 475 * ncomps + c);
        const auto *kg_476 = buffer.data(kg + 476 * ncomps + c);
        const auto *kg_477 = buffer.data(kg + 477 * ncomps + c);
        const auto *kg_478 = buffer.data(kg + 478 * ncomps + c);
        const auto *kg_479 = buffer.data(kg + 479 * ncomps + c);
        const auto *kg_480 = buffer.data(kg + 480 * ncomps + c);
        const auto *kg_481 = buffer.data(kg + 481 * ncomps + c);
        const auto *kg_482 = buffer.data(kg + 482 * ncomps + c);
        const auto *kg_483 = buffer.data(kg + 483 * ncomps + c);
        const auto *kg_484 = buffer.data(kg + 484 * ncomps + c);
        const auto *kg_485 = buffer.data(kg + 485 * ncomps + c);
        const auto *kg_486 = buffer.data(kg + 486 * ncomps + c);
        const auto *kg_487 = buffer.data(kg + 487 * ncomps + c);
        const auto *kg_488 = buffer.data(kg + 488 * ncomps + c);
        const auto *kg_489 = buffer.data(kg + 489 * ncomps + c);
        const auto *kg_490 = buffer.data(kg + 490 * ncomps + c);
        const auto *kg_491 = buffer.data(kg + 491 * ncomps + c);
        const auto *kg_492 = buffer.data(kg + 492 * ncomps + c);
        const auto *kg_493 = buffer.data(kg + 493 * ncomps + c);
        const auto *kg_494 = buffer.data(kg + 494 * ncomps + c);
        const auto *kg_495 = buffer.data(kg + 495 * ncomps + c);
        const auto *kg_496 = buffer.data(kg + 496 * ncomps + c);
        const auto *kg_497 = buffer.data(kg + 497 * ncomps + c);
        const auto *kg_498 = buffer.data(kg + 498 * ncomps + c);
        const auto *kg_499 = buffer.data(kg + 499 * ncomps + c);
        const auto *kg_500 = buffer.data(kg + 500 * ncomps + c);
        const auto *kg_501 = buffer.data(kg + 501 * ncomps + c);
        const auto *kg_502 = buffer.data(kg + 502 * ncomps + c);
        const auto *kg_503 = buffer.data(kg + 503 * ncomps + c);
        const auto *kg_504 = buffer.data(kg + 504 * ncomps + c);
        const auto *kg_505 = buffer.data(kg + 505 * ncomps + c);
        const auto *kg_506 = buffer.data(kg + 506 * ncomps + c);
        const auto *kg_507 = buffer.data(kg + 507 * ncomps + c);
        const auto *kg_508 = buffer.data(kg + 508 * ncomps + c);
        const auto *kg_509 = buffer.data(kg + 509 * ncomps + c);

        const auto *lg_405 = buffer.data(lg + 405 * ncomps + c);
        const auto *lg_406 = buffer.data(lg + 406 * ncomps + c);
        const auto *lg_407 = buffer.data(lg + 407 * ncomps + c);
        const auto *lg_408 = buffer.data(lg + 408 * ncomps + c);
        const auto *lg_409 = buffer.data(lg + 409 * ncomps + c);
        const auto *lg_410 = buffer.data(lg + 410 * ncomps + c);
        const auto *lg_411 = buffer.data(lg + 411 * ncomps + c);
        const auto *lg_412 = buffer.data(lg + 412 * ncomps + c);
        const auto *lg_413 = buffer.data(lg + 413 * ncomps + c);
        const auto *lg_414 = buffer.data(lg + 414 * ncomps + c);
        const auto *lg_415 = buffer.data(lg + 415 * ncomps + c);
        const auto *lg_416 = buffer.data(lg + 416 * ncomps + c);
        const auto *lg_417 = buffer.data(lg + 417 * ncomps + c);
        const auto *lg_418 = buffer.data(lg + 418 * ncomps + c);
        const auto *lg_419 = buffer.data(lg + 419 * ncomps + c);
        const auto *lg_420 = buffer.data(lg + 420 * ncomps + c);
        const auto *lg_421 = buffer.data(lg + 421 * ncomps + c);
        const auto *lg_422 = buffer.data(lg + 422 * ncomps + c);
        const auto *lg_423 = buffer.data(lg + 423 * ncomps + c);
        const auto *lg_424 = buffer.data(lg + 424 * ncomps + c);
        const auto *lg_425 = buffer.data(lg + 425 * ncomps + c);
        const auto *lg_426 = buffer.data(lg + 426 * ncomps + c);
        const auto *lg_427 = buffer.data(lg + 427 * ncomps + c);
        const auto *lg_428 = buffer.data(lg + 428 * ncomps + c);
        const auto *lg_429 = buffer.data(lg + 429 * ncomps + c);
        const auto *lg_430 = buffer.data(lg + 430 * ncomps + c);
        const auto *lg_431 = buffer.data(lg + 431 * ncomps + c);
        const auto *lg_432 = buffer.data(lg + 432 * ncomps + c);
        const auto *lg_433 = buffer.data(lg + 433 * ncomps + c);
        const auto *lg_434 = buffer.data(lg + 434 * ncomps + c);
        const auto *lg_435 = buffer.data(lg + 435 * ncomps + c);
        const auto *lg_436 = buffer.data(lg + 436 * ncomps + c);
        const auto *lg_437 = buffer.data(lg + 437 * ncomps + c);
        const auto *lg_438 = buffer.data(lg + 438 * ncomps + c);
        const auto *lg_439 = buffer.data(lg + 439 * ncomps + c);
        const auto *lg_440 = buffer.data(lg + 440 * ncomps + c);
        const auto *lg_441 = buffer.data(lg + 441 * ncomps + c);
        const auto *lg_442 = buffer.data(lg + 442 * ncomps + c);
        const auto *lg_443 = buffer.data(lg + 443 * ncomps + c);
        const auto *lg_444 = buffer.data(lg + 444 * ncomps + c);
        const auto *lg_445 = buffer.data(lg + 445 * ncomps + c);
        const auto *lg_446 = buffer.data(lg + 446 * ncomps + c);
        const auto *lg_447 = buffer.data(lg + 447 * ncomps + c);
        const auto *lg_448 = buffer.data(lg + 448 * ncomps + c);
        const auto *lg_449 = buffer.data(lg + 449 * ncomps + c);
        const auto *lg_450 = buffer.data(lg + 450 * ncomps + c);
        const auto *lg_451 = buffer.data(lg + 451 * ncomps + c);
        const auto *lg_452 = buffer.data(lg + 452 * ncomps + c);
        const auto *lg_453 = buffer.data(lg + 453 * ncomps + c);
        const auto *lg_454 = buffer.data(lg + 454 * ncomps + c);
        const auto *lg_455 = buffer.data(lg + 455 * ncomps + c);
        const auto *lg_456 = buffer.data(lg + 456 * ncomps + c);
        const auto *lg_457 = buffer.data(lg + 457 * ncomps + c);
        const auto *lg_458 = buffer.data(lg + 458 * ncomps + c);
        const auto *lg_459 = buffer.data(lg + 459 * ncomps + c);
        const auto *lg_460 = buffer.data(lg + 460 * ncomps + c);
        const auto *lg_461 = buffer.data(lg + 461 * ncomps + c);
        const auto *lg_462 = buffer.data(lg + 462 * ncomps + c);
        const auto *lg_463 = buffer.data(lg + 463 * ncomps + c);
        const auto *lg_464 = buffer.data(lg + 464 * ncomps + c);
        const auto *lg_465 = buffer.data(lg + 465 * ncomps + c);
        const auto *lg_466 = buffer.data(lg + 466 * ncomps + c);
        const auto *lg_467 = buffer.data(lg + 467 * ncomps + c);
        const auto *lg_468 = buffer.data(lg + 468 * ncomps + c);
        const auto *lg_469 = buffer.data(lg + 469 * ncomps + c);
        const auto *lg_470 = buffer.data(lg + 470 * ncomps + c);
        const auto *lg_471 = buffer.data(lg + 471 * ncomps + c);
        const auto *lg_472 = buffer.data(lg + 472 * ncomps + c);
        const auto *lg_473 = buffer.data(lg + 473 * ncomps + c);
        const auto *lg_474 = buffer.data(lg + 474 * ncomps + c);
        const auto *lg_475 = buffer.data(lg + 475 * ncomps + c);
        const auto *lg_476 = buffer.data(lg + 476 * ncomps + c);
        const auto *lg_477 = buffer.data(lg + 477 * ncomps + c);
        const auto *lg_478 = buffer.data(lg + 478 * ncomps + c);
        const auto *lg_479 = buffer.data(lg + 479 * ncomps + c);
        const auto *lg_480 = buffer.data(lg + 480 * ncomps + c);
        const auto *lg_481 = buffer.data(lg + 481 * ncomps + c);
        const auto *lg_482 = buffer.data(lg + 482 * ncomps + c);
        const auto *lg_483 = buffer.data(lg + 483 * ncomps + c);
        const auto *lg_484 = buffer.data(lg + 484 * ncomps + c);
        const auto *lg_485 = buffer.data(lg + 485 * ncomps + c);
        const auto *lg_486 = buffer.data(lg + 486 * ncomps + c);
        const auto *lg_487 = buffer.data(lg + 487 * ncomps + c);
        const auto *lg_488 = buffer.data(lg + 488 * ncomps + c);
        const auto *lg_489 = buffer.data(lg + 489 * ncomps + c);
        const auto *lg_490 = buffer.data(lg + 490 * ncomps + c);
        const auto *lg_491 = buffer.data(lg + 491 * ncomps + c);
        const auto *lg_492 = buffer.data(lg + 492 * ncomps + c);
        const auto *lg_493 = buffer.data(lg + 493 * ncomps + c);
        const auto *lg_494 = buffer.data(lg + 494 * ncomps + c);
        const auto *lg_495 = buffer.data(lg + 495 * ncomps + c);
        const auto *lg_496 = buffer.data(lg + 496 * ncomps + c);
        const auto *lg_497 = buffer.data(lg + 497 * ncomps + c);
        const auto *lg_498 = buffer.data(lg + 498 * ncomps + c);
        const auto *lg_499 = buffer.data(lg + 499 * ncomps + c);
        const auto *lg_500 = buffer.data(lg + 500 * ncomps + c);
        const auto *lg_501 = buffer.data(lg + 501 * ncomps + c);
        const auto *lg_502 = buffer.data(lg + 502 * ncomps + c);
        const auto *lg_503 = buffer.data(lg + 503 * ncomps + c);
        const auto *lg_504 = buffer.data(lg + 504 * ncomps + c);
        const auto *lg_505 = buffer.data(lg + 505 * ncomps + c);
        const auto *lg_506 = buffer.data(lg + 506 * ncomps + c);
        const auto *lg_507 = buffer.data(lg + 507 * ncomps + c);
        const auto *lg_508 = buffer.data(lg + 508 * ncomps + c);
        const auto *lg_509 = buffer.data(lg + 509 * ncomps + c);
        const auto *lg_520 = buffer.data(lg + 520 * ncomps + c);
        const auto *lg_521 = buffer.data(lg + 521 * ncomps + c);
        const auto *lg_522 = buffer.data(lg + 522 * ncomps + c);
        const auto *lg_523 = buffer.data(lg + 523 * ncomps + c);
        const auto *lg_524 = buffer.data(lg + 524 * ncomps + c);
        const auto *lg_539 = buffer.data(lg + 539 * ncomps + c);
        const auto *lg_550 = buffer.data(lg + 550 * ncomps + c);
        const auto *lg_551 = buffer.data(lg + 551 * ncomps + c);
        const auto *lg_552 = buffer.data(lg + 552 * ncomps + c);
        const auto *lg_553 = buffer.data(lg + 553 * ncomps + c);
        const auto *lg_554 = buffer.data(lg + 554 * ncomps + c);
        const auto *lg_565 = buffer.data(lg + 565 * ncomps + c);
        const auto *lg_566 = buffer.data(lg + 566 * ncomps + c);
        const auto *lg_567 = buffer.data(lg + 567 * ncomps + c);
        const auto *lg_568 = buffer.data(lg + 568 * ncomps + c);
        const auto *lg_569 = buffer.data(lg + 569 * ncomps + c);
        const auto *lg_580 = buffer.data(lg + 580 * ncomps + c);
        const auto *lg_581 = buffer.data(lg + 581 * ncomps + c);
        const auto *lg_582 = buffer.data(lg + 582 * ncomps + c);
        const auto *lg_583 = buffer.data(lg + 583 * ncomps + c);
        const auto *lg_584 = buffer.data(lg + 584 * ncomps + c);
        const auto *lg_595 = buffer.data(lg + 595 * ncomps + c);
        const auto *lg_596 = buffer.data(lg + 596 * ncomps + c);
        const auto *lg_597 = buffer.data(lg + 597 * ncomps + c);
        const auto *lg_598 = buffer.data(lg + 598 * ncomps + c);
        const auto *lg_599 = buffer.data(lg + 599 * ncomps + c);
        const auto *lg_610 = buffer.data(lg + 610 * ncomps + c);
        const auto *lg_611 = buffer.data(lg + 611 * ncomps + c);
        const auto *lg_612 = buffer.data(lg + 612 * ncomps + c);
        const auto *lg_613 = buffer.data(lg + 613 * ncomps + c);
        const auto *lg_614 = buffer.data(lg + 614 * ncomps + c);
        const auto *lg_625 = buffer.data(lg + 625 * ncomps + c);
        const auto *lg_626 = buffer.data(lg + 626 * ncomps + c);
        const auto *lg_629 = buffer.data(lg + 629 * ncomps + c);

#pragma omp simd aligned(t_566, t_567, t_568, t_569, ab_x, ab_z, kg_404, kg_405, kg_406, \
                         kg_407, lg_405, lg_406, lg_407, lg_524 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_566[k] = ab_z[k] * kg_404[k]
                       + lg_524[k];

            t_567[k] = ab_x[k] * kg_405[k]
                       + lg_405[k];

            t_568[k] = ab_x[k] * kg_406[k]
                       + lg_406[k];

            t_569[k] = ab_x[k] * kg_407[k]
                       + lg_407[k];
        }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, ab_x, kg_408, kg_409, kg_410, \
                         kg_411, kg_412, lg_408, lg_409, lg_410, lg_411, \
                         lg_412 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_570[k] = ab_x[k] * kg_408[k]
                       + lg_408[k];

            t_571[k] = ab_x[k] * kg_409[k]
                       + lg_409[k];

            t_572[k] = ab_x[k] * kg_410[k]
                       + lg_410[k];

            t_573[k] = ab_x[k] * kg_411[k]
                       + lg_411[k];

            t_574[k] = ab_x[k] * kg_412[k]
                       + lg_412[k];
        }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, ab_x, kg_413, kg_414, kg_415, \
                         kg_416, kg_417, lg_413, lg_414, lg_415, lg_416, \
                         lg_417 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_575[k] = ab_x[k] * kg_413[k]
                       + lg_413[k];

            t_576[k] = ab_x[k] * kg_414[k]
                       + lg_414[k];

            t_577[k] = ab_x[k] * kg_415[k]
                       + lg_415[k];

            t_578[k] = ab_x[k] * kg_416[k]
                       + lg_416[k];

            t_579[k] = ab_x[k] * kg_417[k]
                       + lg_417[k];
        }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, ab_x, ab_y, kg_415, kg_416, kg_418, \
                         kg_419, lg_418, lg_419, lg_520, lg_521 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_580[k] = ab_x[k] * kg_418[k]
                       + lg_418[k];

            t_581[k] = ab_x[k] * kg_419[k]
                       + lg_419[k];

            t_582[k] = ab_y[k] * kg_415[k]
                       + lg_520[k];

            t_583[k] = ab_y[k] * kg_416[k]
                       + lg_521[k];
        }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, ab_y, ab_z, kg_417, kg_418, kg_419, \
                         lg_522, lg_523, lg_524, lg_539 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_584[k] = ab_y[k] * kg_417[k]
                       + lg_522[k];

            t_585[k] = ab_y[k] * kg_418[k]
                       + lg_523[k];

            t_586[k] = ab_y[k] * kg_419[k]
                       + lg_524[k];

            t_587[k] = ab_z[k] * kg_419[k]
                       + lg_539[k];
        }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, t_592, ab_x, kg_420, kg_421, kg_422, \
                         kg_423, kg_424, lg_420, lg_421, lg_422, lg_423, \
                         lg_424 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_588[k] = ab_x[k] * kg_420[k]
                       + lg_420[k];

            t_589[k] = ab_x[k] * kg_421[k]
                       + lg_421[k];

            t_590[k] = ab_x[k] * kg_422[k]
                       + lg_422[k];

            t_591[k] = ab_x[k] * kg_423[k]
                       + lg_423[k];

            t_592[k] = ab_x[k] * kg_424[k]
                       + lg_424[k];
        }

#pragma omp simd aligned(t_593, t_594, t_595, t_596, t_597, ab_x, kg_425, kg_426, kg_427, \
                         kg_428, kg_429, lg_425, lg_426, lg_427, lg_428, \
                         lg_429 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_593[k] = ab_x[k] * kg_425[k]
                       + lg_425[k];

            t_594[k] = ab_x[k] * kg_426[k]
                       + lg_426[k];

            t_595[k] = ab_x[k] * kg_427[k]
                       + lg_427[k];

            t_596[k] = ab_x[k] * kg_428[k]
                       + lg_428[k];

            t_597[k] = ab_x[k] * kg_429[k]
                       + lg_429[k];
        }

#pragma omp simd aligned(t_598, t_599, t_600, t_601, t_602, ab_x, kg_430, kg_431, kg_432, \
                         kg_433, kg_434, lg_430, lg_431, lg_432, lg_433, \
                         lg_434 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_598[k] = ab_x[k] * kg_430[k]
                       + lg_430[k];

            t_599[k] = ab_x[k] * kg_431[k]
                       + lg_431[k];

            t_600[k] = ab_x[k] * kg_432[k]
                       + lg_432[k];

            t_601[k] = ab_x[k] * kg_433[k]
                       + lg_433[k];

            t_602[k] = ab_x[k] * kg_434[k]
                       + lg_434[k];
        }

#pragma omp simd aligned(t_603, t_604, t_605, t_606, t_607, ab_y, kg_430, kg_431, kg_432, \
                         kg_433, kg_434, lg_550, lg_551, lg_552, lg_553, \
                         lg_554 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_603[k] = ab_y[k] * kg_430[k]
                       + lg_550[k];

            t_604[k] = ab_y[k] * kg_431[k]
                       + lg_551[k];

            t_605[k] = ab_y[k] * kg_432[k]
                       + lg_552[k];

            t_606[k] = ab_y[k] * kg_433[k]
                       + lg_553[k];

            t_607[k] = ab_y[k] * kg_434[k]
                       + lg_554[k];
        }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, ab_x, ab_z, kg_434, kg_435, kg_436, \
                         kg_437, lg_435, lg_436, lg_437, lg_569 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_608[k] = ab_z[k] * kg_434[k]
                       + lg_569[k];

            t_609[k] = ab_x[k] * kg_435[k]
                       + lg_435[k];

            t_610[k] = ab_x[k] * kg_436[k]
                       + lg_436[k];

            t_611[k] = ab_x[k] * kg_437[k]
                       + lg_437[k];
        }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, t_616, ab_x, kg_438, kg_439, kg_440, \
                         kg_441, kg_442, lg_438, lg_439, lg_440, lg_441, \
                         lg_442 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_612[k] = ab_x[k] * kg_438[k]
                       + lg_438[k];

            t_613[k] = ab_x[k] * kg_439[k]
                       + lg_439[k];

            t_614[k] = ab_x[k] * kg_440[k]
                       + lg_440[k];

            t_615[k] = ab_x[k] * kg_441[k]
                       + lg_441[k];

            t_616[k] = ab_x[k] * kg_442[k]
                       + lg_442[k];
        }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, t_621, ab_x, kg_443, kg_444, kg_445, \
                         kg_446, kg_447, lg_443, lg_444, lg_445, lg_446, \
                         lg_447 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_617[k] = ab_x[k] * kg_443[k]
                       + lg_443[k];

            t_618[k] = ab_x[k] * kg_444[k]
                       + lg_444[k];

            t_619[k] = ab_x[k] * kg_445[k]
                       + lg_445[k];

            t_620[k] = ab_x[k] * kg_446[k]
                       + lg_446[k];

            t_621[k] = ab_x[k] * kg_447[k]
                       + lg_447[k];
        }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, ab_x, ab_y, kg_445, kg_446, kg_448, \
                         kg_449, lg_448, lg_449, lg_565, lg_566 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_622[k] = ab_x[k] * kg_448[k]
                       + lg_448[k];

            t_623[k] = ab_x[k] * kg_449[k]
                       + lg_449[k];

            t_624[k] = ab_y[k] * kg_445[k]
                       + lg_565[k];

            t_625[k] = ab_y[k] * kg_446[k]
                       + lg_566[k];
        }

#pragma omp simd aligned(t_626, t_627, t_628, t_629, ab_y, ab_z, kg_447, kg_448, kg_449, \
                         lg_567, lg_568, lg_569, lg_584 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_626[k] = ab_y[k] * kg_447[k]
                       + lg_567[k];

            t_627[k] = ab_y[k] * kg_448[k]
                       + lg_568[k];

            t_628[k] = ab_y[k] * kg_449[k]
                       + lg_569[k];

            t_629[k] = ab_z[k] * kg_449[k]
                       + lg_584[k];
        }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, ab_x, kg_450, kg_451, kg_452, \
                         kg_453, kg_454, lg_450, lg_451, lg_452, lg_453, \
                         lg_454 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_630[k] = ab_x[k] * kg_450[k]
                       + lg_450[k];

            t_631[k] = ab_x[k] * kg_451[k]
                       + lg_451[k];

            t_632[k] = ab_x[k] * kg_452[k]
                       + lg_452[k];

            t_633[k] = ab_x[k] * kg_453[k]
                       + lg_453[k];

            t_634[k] = ab_x[k] * kg_454[k]
                       + lg_454[k];
        }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, ab_x, kg_455, kg_456, kg_457, \
                         kg_458, kg_459, lg_455, lg_456, lg_457, lg_458, \
                         lg_459 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_635[k] = ab_x[k] * kg_455[k]
                       + lg_455[k];

            t_636[k] = ab_x[k] * kg_456[k]
                       + lg_456[k];

            t_637[k] = ab_x[k] * kg_457[k]
                       + lg_457[k];

            t_638[k] = ab_x[k] * kg_458[k]
                       + lg_458[k];

            t_639[k] = ab_x[k] * kg_459[k]
                       + lg_459[k];
        }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, ab_x, kg_460, kg_461, kg_462, \
                         kg_463, kg_464, lg_460, lg_461, lg_462, lg_463, \
                         lg_464 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_640[k] = ab_x[k] * kg_460[k]
                       + lg_460[k];

            t_641[k] = ab_x[k] * kg_461[k]
                       + lg_461[k];

            t_642[k] = ab_x[k] * kg_462[k]
                       + lg_462[k];

            t_643[k] = ab_x[k] * kg_463[k]
                       + lg_463[k];

            t_644[k] = ab_x[k] * kg_464[k]
                       + lg_464[k];
        }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, ab_y, kg_460, kg_461, kg_462, \
                         kg_463, kg_464, lg_580, lg_581, lg_582, lg_583, \
                         lg_584 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_645[k] = ab_y[k] * kg_460[k]
                       + lg_580[k];

            t_646[k] = ab_y[k] * kg_461[k]
                       + lg_581[k];

            t_647[k] = ab_y[k] * kg_462[k]
                       + lg_582[k];

            t_648[k] = ab_y[k] * kg_463[k]
                       + lg_583[k];

            t_649[k] = ab_y[k] * kg_464[k]
                       + lg_584[k];
        }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, ab_x, ab_z, kg_464, kg_465, kg_466, \
                         kg_467, lg_465, lg_466, lg_467, lg_599 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_650[k] = ab_z[k] * kg_464[k]
                       + lg_599[k];

            t_651[k] = ab_x[k] * kg_465[k]
                       + lg_465[k];

            t_652[k] = ab_x[k] * kg_466[k]
                       + lg_466[k];

            t_653[k] = ab_x[k] * kg_467[k]
                       + lg_467[k];
        }

#pragma omp simd aligned(t_654, t_655, t_656, t_657, t_658, ab_x, kg_468, kg_469, kg_470, \
                         kg_471, kg_472, lg_468, lg_469, lg_470, lg_471, \
                         lg_472 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_654[k] = ab_x[k] * kg_468[k]
                       + lg_468[k];

            t_655[k] = ab_x[k] * kg_469[k]
                       + lg_469[k];

            t_656[k] = ab_x[k] * kg_470[k]
                       + lg_470[k];

            t_657[k] = ab_x[k] * kg_471[k]
                       + lg_471[k];

            t_658[k] = ab_x[k] * kg_472[k]
                       + lg_472[k];
        }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, ab_x, kg_473, kg_474, kg_475, \
                         kg_476, kg_477, lg_473, lg_474, lg_475, lg_476, \
                         lg_477 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_659[k] = ab_x[k] * kg_473[k]
                       + lg_473[k];

            t_660[k] = ab_x[k] * kg_474[k]
                       + lg_474[k];

            t_661[k] = ab_x[k] * kg_475[k]
                       + lg_475[k];

            t_662[k] = ab_x[k] * kg_476[k]
                       + lg_476[k];

            t_663[k] = ab_x[k] * kg_477[k]
                       + lg_477[k];
        }

#pragma omp simd aligned(t_664, t_665, t_666, t_667, ab_x, ab_y, kg_475, kg_476, kg_478, \
                         kg_479, lg_478, lg_479, lg_595, lg_596 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_664[k] = ab_x[k] * kg_478[k]
                       + lg_478[k];

            t_665[k] = ab_x[k] * kg_479[k]
                       + lg_479[k];

            t_666[k] = ab_y[k] * kg_475[k]
                       + lg_595[k];

            t_667[k] = ab_y[k] * kg_476[k]
                       + lg_596[k];
        }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, ab_y, ab_z, kg_477, kg_478, kg_479, \
                         lg_597, lg_598, lg_599, lg_614 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_668[k] = ab_y[k] * kg_477[k]
                       + lg_597[k];

            t_669[k] = ab_y[k] * kg_478[k]
                       + lg_598[k];

            t_670[k] = ab_y[k] * kg_479[k]
                       + lg_599[k];

            t_671[k] = ab_z[k] * kg_479[k]
                       + lg_614[k];
        }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, t_676, ab_x, kg_480, kg_481, kg_482, \
                         kg_483, kg_484, lg_480, lg_481, lg_482, lg_483, \
                         lg_484 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_672[k] = ab_x[k] * kg_480[k]
                       + lg_480[k];

            t_673[k] = ab_x[k] * kg_481[k]
                       + lg_481[k];

            t_674[k] = ab_x[k] * kg_482[k]
                       + lg_482[k];

            t_675[k] = ab_x[k] * kg_483[k]
                       + lg_483[k];

            t_676[k] = ab_x[k] * kg_484[k]
                       + lg_484[k];
        }

#pragma omp simd aligned(t_677, t_678, t_679, t_680, t_681, ab_x, kg_485, kg_486, kg_487, \
                         kg_488, kg_489, lg_485, lg_486, lg_487, lg_488, \
                         lg_489 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_677[k] = ab_x[k] * kg_485[k]
                       + lg_485[k];

            t_678[k] = ab_x[k] * kg_486[k]
                       + lg_486[k];

            t_679[k] = ab_x[k] * kg_487[k]
                       + lg_487[k];

            t_680[k] = ab_x[k] * kg_488[k]
                       + lg_488[k];

            t_681[k] = ab_x[k] * kg_489[k]
                       + lg_489[k];
        }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, t_686, ab_x, kg_490, kg_491, kg_492, \
                         kg_493, kg_494, lg_490, lg_491, lg_492, lg_493, \
                         lg_494 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_682[k] = ab_x[k] * kg_490[k]
                       + lg_490[k];

            t_683[k] = ab_x[k] * kg_491[k]
                       + lg_491[k];

            t_684[k] = ab_x[k] * kg_492[k]
                       + lg_492[k];

            t_685[k] = ab_x[k] * kg_493[k]
                       + lg_493[k];

            t_686[k] = ab_x[k] * kg_494[k]
                       + lg_494[k];
        }

#pragma omp simd aligned(t_687, t_688, t_689, t_690, t_691, ab_y, kg_490, kg_491, kg_492, \
                         kg_493, kg_494, lg_610, lg_611, lg_612, lg_613, \
                         lg_614 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_687[k] = ab_y[k] * kg_490[k]
                       + lg_610[k];

            t_688[k] = ab_y[k] * kg_491[k]
                       + lg_611[k];

            t_689[k] = ab_y[k] * kg_492[k]
                       + lg_612[k];

            t_690[k] = ab_y[k] * kg_493[k]
                       + lg_613[k];

            t_691[k] = ab_y[k] * kg_494[k]
                       + lg_614[k];
        }

#pragma omp simd aligned(t_692, t_693, t_694, t_695, ab_x, ab_z, kg_494, kg_495, kg_496, \
                         kg_497, lg_495, lg_496, lg_497, lg_629 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_692[k] = ab_z[k] * kg_494[k]
                       + lg_629[k];

            t_693[k] = ab_x[k] * kg_495[k]
                       + lg_495[k];

            t_694[k] = ab_x[k] * kg_496[k]
                       + lg_496[k];

            t_695[k] = ab_x[k] * kg_497[k]
                       + lg_497[k];
        }

#pragma omp simd aligned(t_696, t_697, t_698, t_699, t_700, ab_x, kg_498, kg_499, kg_500, \
                         kg_501, kg_502, lg_498, lg_499, lg_500, lg_501, \
                         lg_502 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_696[k] = ab_x[k] * kg_498[k]
                       + lg_498[k];

            t_697[k] = ab_x[k] * kg_499[k]
                       + lg_499[k];

            t_698[k] = ab_x[k] * kg_500[k]
                       + lg_500[k];

            t_699[k] = ab_x[k] * kg_501[k]
                       + lg_501[k];

            t_700[k] = ab_x[k] * kg_502[k]
                       + lg_502[k];
        }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, t_705, ab_x, kg_503, kg_504, kg_505, \
                         kg_506, kg_507, lg_503, lg_504, lg_505, lg_506, \
                         lg_507 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_701[k] = ab_x[k] * kg_503[k]
                       + lg_503[k];

            t_702[k] = ab_x[k] * kg_504[k]
                       + lg_504[k];

            t_703[k] = ab_x[k] * kg_505[k]
                       + lg_505[k];

            t_704[k] = ab_x[k] * kg_506[k]
                       + lg_506[k];

            t_705[k] = ab_x[k] * kg_507[k]
                       + lg_507[k];
        }

#pragma omp simd aligned(t_706, t_707, t_708, t_709, ab_x, ab_y, kg_505, kg_506, kg_508, \
                         kg_509, lg_508, lg_509, lg_625, lg_626 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_706[k] = ab_x[k] * kg_508[k]
                       + lg_508[k];

            t_707[k] = ab_x[k] * kg_509[k]
                       + lg_509[k];

            t_708[k] = ab_y[k] * kg_505[k]
                       + lg_625[k];

            t_709[k] = ab_y[k] * kg_506[k]
                       + lg_626[k];
        }
    }
}

static auto
compute_hrr_kh_out_of_first_piece5(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t kg, const size_t lg,
                                   const size_t ncomps, const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_710 = buffer.data(target + 710 * ncomps + c);
        auto *t_711 = buffer.data(target + 711 * ncomps + c);
        auto *t_712 = buffer.data(target + 712 * ncomps + c);
        auto *t_713 = buffer.data(target + 713 * ncomps + c);
        auto *t_714 = buffer.data(target + 714 * ncomps + c);
        auto *t_715 = buffer.data(target + 715 * ncomps + c);
        auto *t_716 = buffer.data(target + 716 * ncomps + c);
        auto *t_717 = buffer.data(target + 717 * ncomps + c);
        auto *t_718 = buffer.data(target + 718 * ncomps + c);
        auto *t_719 = buffer.data(target + 719 * ncomps + c);
        auto *t_720 = buffer.data(target + 720 * ncomps + c);
        auto *t_721 = buffer.data(target + 721 * ncomps + c);
        auto *t_722 = buffer.data(target + 722 * ncomps + c);
        auto *t_723 = buffer.data(target + 723 * ncomps + c);
        auto *t_724 = buffer.data(target + 724 * ncomps + c);
        auto *t_725 = buffer.data(target + 725 * ncomps + c);
        auto *t_726 = buffer.data(target + 726 * ncomps + c);
        auto *t_727 = buffer.data(target + 727 * ncomps + c);
        auto *t_728 = buffer.data(target + 728 * ncomps + c);
        auto *t_729 = buffer.data(target + 729 * ncomps + c);
        auto *t_730 = buffer.data(target + 730 * ncomps + c);
        auto *t_731 = buffer.data(target + 731 * ncomps + c);
        auto *t_732 = buffer.data(target + 732 * ncomps + c);
        auto *t_733 = buffer.data(target + 733 * ncomps + c);
        auto *t_734 = buffer.data(target + 734 * ncomps + c);
        auto *t_735 = buffer.data(target + 735 * ncomps + c);
        auto *t_736 = buffer.data(target + 736 * ncomps + c);
        auto *t_737 = buffer.data(target + 737 * ncomps + c);
        auto *t_738 = buffer.data(target + 738 * ncomps + c);
        auto *t_739 = buffer.data(target + 739 * ncomps + c);
        auto *t_740 = buffer.data(target + 740 * ncomps + c);
        auto *t_741 = buffer.data(target + 741 * ncomps + c);
        auto *t_742 = buffer.data(target + 742 * ncomps + c);
        auto *t_743 = buffer.data(target + 743 * ncomps + c);
        auto *t_744 = buffer.data(target + 744 * ncomps + c);
        auto *t_745 = buffer.data(target + 745 * ncomps + c);
        auto *t_746 = buffer.data(target + 746 * ncomps + c);
        auto *t_747 = buffer.data(target + 747 * ncomps + c);
        auto *t_748 = buffer.data(target + 748 * ncomps + c);
        auto *t_749 = buffer.data(target + 749 * ncomps + c);
        auto *t_750 = buffer.data(target + 750 * ncomps + c);
        auto *t_751 = buffer.data(target + 751 * ncomps + c);
        auto *t_752 = buffer.data(target + 752 * ncomps + c);
        auto *t_753 = buffer.data(target + 753 * ncomps + c);
        auto *t_754 = buffer.data(target + 754 * ncomps + c);
        auto *t_755 = buffer.data(target + 755 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *kg_507 = buffer.data(kg + 507 * ncomps + c);
        const auto *kg_508 = buffer.data(kg + 508 * ncomps + c);
        const auto *kg_509 = buffer.data(kg + 509 * ncomps + c);
        const auto *kg_510 = buffer.data(kg + 510 * ncomps + c);
        const auto *kg_511 = buffer.data(kg + 511 * ncomps + c);
        const auto *kg_512 = buffer.data(kg + 512 * ncomps + c);
        const auto *kg_513 = buffer.data(kg + 513 * ncomps + c);
        const auto *kg_514 = buffer.data(kg + 514 * ncomps + c);
        const auto *kg_515 = buffer.data(kg + 515 * ncomps + c);
        const auto *kg_516 = buffer.data(kg + 516 * ncomps + c);
        const auto *kg_517 = buffer.data(kg + 517 * ncomps + c);
        const auto *kg_518 = buffer.data(kg + 518 * ncomps + c);
        const auto *kg_519 = buffer.data(kg + 519 * ncomps + c);
        const auto *kg_520 = buffer.data(kg + 520 * ncomps + c);
        const auto *kg_521 = buffer.data(kg + 521 * ncomps + c);
        const auto *kg_522 = buffer.data(kg + 522 * ncomps + c);
        const auto *kg_523 = buffer.data(kg + 523 * ncomps + c);
        const auto *kg_524 = buffer.data(kg + 524 * ncomps + c);
        const auto *kg_525 = buffer.data(kg + 525 * ncomps + c);
        const auto *kg_526 = buffer.data(kg + 526 * ncomps + c);
        const auto *kg_527 = buffer.data(kg + 527 * ncomps + c);
        const auto *kg_528 = buffer.data(kg + 528 * ncomps + c);
        const auto *kg_529 = buffer.data(kg + 529 * ncomps + c);
        const auto *kg_530 = buffer.data(kg + 530 * ncomps + c);
        const auto *kg_531 = buffer.data(kg + 531 * ncomps + c);
        const auto *kg_532 = buffer.data(kg + 532 * ncomps + c);
        const auto *kg_533 = buffer.data(kg + 533 * ncomps + c);
        const auto *kg_534 = buffer.data(kg + 534 * ncomps + c);
        const auto *kg_535 = buffer.data(kg + 535 * ncomps + c);
        const auto *kg_536 = buffer.data(kg + 536 * ncomps + c);
        const auto *kg_537 = buffer.data(kg + 537 * ncomps + c);
        const auto *kg_538 = buffer.data(kg + 538 * ncomps + c);
        const auto *kg_539 = buffer.data(kg + 539 * ncomps + c);

        const auto *lg_510 = buffer.data(lg + 510 * ncomps + c);
        const auto *lg_511 = buffer.data(lg + 511 * ncomps + c);
        const auto *lg_512 = buffer.data(lg + 512 * ncomps + c);
        const auto *lg_513 = buffer.data(lg + 513 * ncomps + c);
        const auto *lg_514 = buffer.data(lg + 514 * ncomps + c);
        const auto *lg_515 = buffer.data(lg + 515 * ncomps + c);
        const auto *lg_516 = buffer.data(lg + 516 * ncomps + c);
        const auto *lg_517 = buffer.data(lg + 517 * ncomps + c);
        const auto *lg_518 = buffer.data(lg + 518 * ncomps + c);
        const auto *lg_519 = buffer.data(lg + 519 * ncomps + c);
        const auto *lg_520 = buffer.data(lg + 520 * ncomps + c);
        const auto *lg_521 = buffer.data(lg + 521 * ncomps + c);
        const auto *lg_522 = buffer.data(lg + 522 * ncomps + c);
        const auto *lg_523 = buffer.data(lg + 523 * ncomps + c);
        const auto *lg_524 = buffer.data(lg + 524 * ncomps + c);
        const auto *lg_525 = buffer.data(lg + 525 * ncomps + c);
        const auto *lg_526 = buffer.data(lg + 526 * ncomps + c);
        const auto *lg_527 = buffer.data(lg + 527 * ncomps + c);
        const auto *lg_528 = buffer.data(lg + 528 * ncomps + c);
        const auto *lg_529 = buffer.data(lg + 529 * ncomps + c);
        const auto *lg_530 = buffer.data(lg + 530 * ncomps + c);
        const auto *lg_531 = buffer.data(lg + 531 * ncomps + c);
        const auto *lg_532 = buffer.data(lg + 532 * ncomps + c);
        const auto *lg_533 = buffer.data(lg + 533 * ncomps + c);
        const auto *lg_534 = buffer.data(lg + 534 * ncomps + c);
        const auto *lg_535 = buffer.data(lg + 535 * ncomps + c);
        const auto *lg_536 = buffer.data(lg + 536 * ncomps + c);
        const auto *lg_537 = buffer.data(lg + 537 * ncomps + c);
        const auto *lg_538 = buffer.data(lg + 538 * ncomps + c);
        const auto *lg_539 = buffer.data(lg + 539 * ncomps + c);
        const auto *lg_627 = buffer.data(lg + 627 * ncomps + c);
        const auto *lg_628 = buffer.data(lg + 628 * ncomps + c);
        const auto *lg_629 = buffer.data(lg + 629 * ncomps + c);
        const auto *lg_640 = buffer.data(lg + 640 * ncomps + c);
        const auto *lg_641 = buffer.data(lg + 641 * ncomps + c);
        const auto *lg_642 = buffer.data(lg + 642 * ncomps + c);
        const auto *lg_643 = buffer.data(lg + 643 * ncomps + c);
        const auto *lg_644 = buffer.data(lg + 644 * ncomps + c);
        const auto *lg_655 = buffer.data(lg + 655 * ncomps + c);
        const auto *lg_656 = buffer.data(lg + 656 * ncomps + c);
        const auto *lg_657 = buffer.data(lg + 657 * ncomps + c);
        const auto *lg_658 = buffer.data(lg + 658 * ncomps + c);
        const auto *lg_659 = buffer.data(lg + 659 * ncomps + c);
        const auto *lg_674 = buffer.data(lg + 674 * ncomps + c);

#pragma omp simd aligned(t_710, t_711, t_712, t_713, ab_y, ab_z, kg_507, kg_508, kg_509, \
                         lg_627, lg_628, lg_629, lg_644 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_710[k] = ab_y[k] * kg_507[k]
                       + lg_627[k];

            t_711[k] = ab_y[k] * kg_508[k]
                       + lg_628[k];

            t_712[k] = ab_y[k] * kg_509[k]
                       + lg_629[k];

            t_713[k] = ab_z[k] * kg_509[k]
                       + lg_644[k];
        }

#pragma omp simd aligned(t_714, t_715, t_716, t_717, t_718, ab_x, kg_510, kg_511, kg_512, \
                         kg_513, kg_514, lg_510, lg_511, lg_512, lg_513, \
                         lg_514 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_714[k] = ab_x[k] * kg_510[k]
                       + lg_510[k];

            t_715[k] = ab_x[k] * kg_511[k]
                       + lg_511[k];

            t_716[k] = ab_x[k] * kg_512[k]
                       + lg_512[k];

            t_717[k] = ab_x[k] * kg_513[k]
                       + lg_513[k];

            t_718[k] = ab_x[k] * kg_514[k]
                       + lg_514[k];
        }

#pragma omp simd aligned(t_719, t_720, t_721, t_722, t_723, ab_x, kg_515, kg_516, kg_517, \
                         kg_518, kg_519, lg_515, lg_516, lg_517, lg_518, \
                         lg_519 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_719[k] = ab_x[k] * kg_515[k]
                       + lg_515[k];

            t_720[k] = ab_x[k] * kg_516[k]
                       + lg_516[k];

            t_721[k] = ab_x[k] * kg_517[k]
                       + lg_517[k];

            t_722[k] = ab_x[k] * kg_518[k]
                       + lg_518[k];

            t_723[k] = ab_x[k] * kg_519[k]
                       + lg_519[k];
        }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, ab_x, kg_520, kg_521, kg_522, \
                         kg_523, kg_524, lg_520, lg_521, lg_522, lg_523, \
                         lg_524 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_724[k] = ab_x[k] * kg_520[k]
                       + lg_520[k];

            t_725[k] = ab_x[k] * kg_521[k]
                       + lg_521[k];

            t_726[k] = ab_x[k] * kg_522[k]
                       + lg_522[k];

            t_727[k] = ab_x[k] * kg_523[k]
                       + lg_523[k];

            t_728[k] = ab_x[k] * kg_524[k]
                       + lg_524[k];
        }

#pragma omp simd aligned(t_729, t_730, t_731, t_732, t_733, ab_y, kg_520, kg_521, kg_522, \
                         kg_523, kg_524, lg_640, lg_641, lg_642, lg_643, \
                         lg_644 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_729[k] = ab_y[k] * kg_520[k]
                       + lg_640[k];

            t_730[k] = ab_y[k] * kg_521[k]
                       + lg_641[k];

            t_731[k] = ab_y[k] * kg_522[k]
                       + lg_642[k];

            t_732[k] = ab_y[k] * kg_523[k]
                       + lg_643[k];

            t_733[k] = ab_y[k] * kg_524[k]
                       + lg_644[k];
        }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, ab_x, ab_z, kg_524, kg_525, kg_526, \
                         kg_527, lg_525, lg_526, lg_527, lg_659 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_734[k] = ab_z[k] * kg_524[k]
                       + lg_659[k];

            t_735[k] = ab_x[k] * kg_525[k]
                       + lg_525[k];

            t_736[k] = ab_x[k] * kg_526[k]
                       + lg_526[k];

            t_737[k] = ab_x[k] * kg_527[k]
                       + lg_527[k];
        }

#pragma omp simd aligned(t_738, t_739, t_740, t_741, t_742, ab_x, kg_528, kg_529, kg_530, \
                         kg_531, kg_532, lg_528, lg_529, lg_530, lg_531, \
                         lg_532 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_738[k] = ab_x[k] * kg_528[k]
                       + lg_528[k];

            t_739[k] = ab_x[k] * kg_529[k]
                       + lg_529[k];

            t_740[k] = ab_x[k] * kg_530[k]
                       + lg_530[k];

            t_741[k] = ab_x[k] * kg_531[k]
                       + lg_531[k];

            t_742[k] = ab_x[k] * kg_532[k]
                       + lg_532[k];
        }

#pragma omp simd aligned(t_743, t_744, t_745, t_746, t_747, ab_x, kg_533, kg_534, kg_535, \
                         kg_536, kg_537, lg_533, lg_534, lg_535, lg_536, \
                         lg_537 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_743[k] = ab_x[k] * kg_533[k]
                       + lg_533[k];

            t_744[k] = ab_x[k] * kg_534[k]
                       + lg_534[k];

            t_745[k] = ab_x[k] * kg_535[k]
                       + lg_535[k];

            t_746[k] = ab_x[k] * kg_536[k]
                       + lg_536[k];

            t_747[k] = ab_x[k] * kg_537[k]
                       + lg_537[k];
        }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, ab_x, ab_y, kg_535, kg_536, kg_538, \
                         kg_539, lg_538, lg_539, lg_655, lg_656 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_748[k] = ab_x[k] * kg_538[k]
                       + lg_538[k];

            t_749[k] = ab_x[k] * kg_539[k]
                       + lg_539[k];

            t_750[k] = ab_y[k] * kg_535[k]
                       + lg_655[k];

            t_751[k] = ab_y[k] * kg_536[k]
                       + lg_656[k];
        }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, ab_y, ab_z, kg_537, kg_538, kg_539, \
                         lg_657, lg_658, lg_659, lg_674 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_752[k] = ab_y[k] * kg_537[k]
                       + lg_657[k];

            t_753[k] = ab_y[k] * kg_538[k]
                       + lg_658[k];

            t_754[k] = ab_y[k] * kg_539[k]
                       + lg_659[k];

            t_755[k] = ab_z[k] * kg_539[k]
                       + lg_674[k];
        }
    }
}

auto
compute_hrr_kh_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t kg, const size_t lg,
                            const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_kh_out_of_first_piece0(buffer, coordinates, target, kg, lg, ncomps, nmax);

    compute_hrr_kh_out_of_first_piece1(buffer, coordinates, target, kg, lg, ncomps, nmax);

    compute_hrr_kh_out_of_first_piece2(buffer, coordinates, target, kg, lg, ncomps, nmax);

    compute_hrr_kh_out_of_first_piece3(buffer, coordinates, target, kg, lg, ncomps, nmax);

    compute_hrr_kh_out_of_first_piece4(buffer, coordinates, target, kg, lg, ncomps, nmax);

    compute_hrr_kh_out_of_first_piece5(buffer, coordinates, target, kg, lg, ncomps, nmax);
}

}  // namespace simdtrf
