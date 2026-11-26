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


#include "SimdTransferIH.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_ih_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ig, const size_t kg,
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

        const auto *ig_0 = buffer.data(ig + 0 * ncomps + c);
        const auto *ig_1 = buffer.data(ig + 1 * ncomps + c);
        const auto *ig_2 = buffer.data(ig + 2 * ncomps + c);
        const auto *ig_3 = buffer.data(ig + 3 * ncomps + c);
        const auto *ig_4 = buffer.data(ig + 4 * ncomps + c);
        const auto *ig_5 = buffer.data(ig + 5 * ncomps + c);
        const auto *ig_6 = buffer.data(ig + 6 * ncomps + c);
        const auto *ig_7 = buffer.data(ig + 7 * ncomps + c);
        const auto *ig_8 = buffer.data(ig + 8 * ncomps + c);
        const auto *ig_9 = buffer.data(ig + 9 * ncomps + c);
        const auto *ig_10 = buffer.data(ig + 10 * ncomps + c);
        const auto *ig_11 = buffer.data(ig + 11 * ncomps + c);
        const auto *ig_12 = buffer.data(ig + 12 * ncomps + c);
        const auto *ig_13 = buffer.data(ig + 13 * ncomps + c);
        const auto *ig_14 = buffer.data(ig + 14 * ncomps + c);
        const auto *ig_15 = buffer.data(ig + 15 * ncomps + c);
        const auto *ig_16 = buffer.data(ig + 16 * ncomps + c);
        const auto *ig_17 = buffer.data(ig + 17 * ncomps + c);
        const auto *ig_18 = buffer.data(ig + 18 * ncomps + c);
        const auto *ig_19 = buffer.data(ig + 19 * ncomps + c);
        const auto *ig_20 = buffer.data(ig + 20 * ncomps + c);
        const auto *ig_21 = buffer.data(ig + 21 * ncomps + c);
        const auto *ig_22 = buffer.data(ig + 22 * ncomps + c);
        const auto *ig_23 = buffer.data(ig + 23 * ncomps + c);
        const auto *ig_24 = buffer.data(ig + 24 * ncomps + c);
        const auto *ig_25 = buffer.data(ig + 25 * ncomps + c);
        const auto *ig_26 = buffer.data(ig + 26 * ncomps + c);
        const auto *ig_27 = buffer.data(ig + 27 * ncomps + c);
        const auto *ig_28 = buffer.data(ig + 28 * ncomps + c);
        const auto *ig_29 = buffer.data(ig + 29 * ncomps + c);
        const auto *ig_30 = buffer.data(ig + 30 * ncomps + c);
        const auto *ig_31 = buffer.data(ig + 31 * ncomps + c);
        const auto *ig_32 = buffer.data(ig + 32 * ncomps + c);
        const auto *ig_33 = buffer.data(ig + 33 * ncomps + c);
        const auto *ig_34 = buffer.data(ig + 34 * ncomps + c);
        const auto *ig_35 = buffer.data(ig + 35 * ncomps + c);
        const auto *ig_36 = buffer.data(ig + 36 * ncomps + c);
        const auto *ig_37 = buffer.data(ig + 37 * ncomps + c);
        const auto *ig_38 = buffer.data(ig + 38 * ncomps + c);
        const auto *ig_39 = buffer.data(ig + 39 * ncomps + c);
        const auto *ig_40 = buffer.data(ig + 40 * ncomps + c);
        const auto *ig_41 = buffer.data(ig + 41 * ncomps + c);
        const auto *ig_42 = buffer.data(ig + 42 * ncomps + c);
        const auto *ig_43 = buffer.data(ig + 43 * ncomps + c);
        const auto *ig_44 = buffer.data(ig + 44 * ncomps + c);
        const auto *ig_45 = buffer.data(ig + 45 * ncomps + c);
        const auto *ig_46 = buffer.data(ig + 46 * ncomps + c);
        const auto *ig_47 = buffer.data(ig + 47 * ncomps + c);
        const auto *ig_48 = buffer.data(ig + 48 * ncomps + c);
        const auto *ig_49 = buffer.data(ig + 49 * ncomps + c);
        const auto *ig_50 = buffer.data(ig + 50 * ncomps + c);
        const auto *ig_51 = buffer.data(ig + 51 * ncomps + c);
        const auto *ig_52 = buffer.data(ig + 52 * ncomps + c);
        const auto *ig_53 = buffer.data(ig + 53 * ncomps + c);
        const auto *ig_54 = buffer.data(ig + 54 * ncomps + c);
        const auto *ig_55 = buffer.data(ig + 55 * ncomps + c);
        const auto *ig_56 = buffer.data(ig + 56 * ncomps + c);
        const auto *ig_57 = buffer.data(ig + 57 * ncomps + c);
        const auto *ig_58 = buffer.data(ig + 58 * ncomps + c);
        const auto *ig_59 = buffer.data(ig + 59 * ncomps + c);
        const auto *ig_60 = buffer.data(ig + 60 * ncomps + c);
        const auto *ig_61 = buffer.data(ig + 61 * ncomps + c);
        const auto *ig_62 = buffer.data(ig + 62 * ncomps + c);
        const auto *ig_63 = buffer.data(ig + 63 * ncomps + c);
        const auto *ig_64 = buffer.data(ig + 64 * ncomps + c);
        const auto *ig_65 = buffer.data(ig + 65 * ncomps + c);
        const auto *ig_66 = buffer.data(ig + 66 * ncomps + c);
        const auto *ig_67 = buffer.data(ig + 67 * ncomps + c);
        const auto *ig_68 = buffer.data(ig + 68 * ncomps + c);
        const auto *ig_69 = buffer.data(ig + 69 * ncomps + c);
        const auto *ig_70 = buffer.data(ig + 70 * ncomps + c);
        const auto *ig_71 = buffer.data(ig + 71 * ncomps + c);
        const auto *ig_72 = buffer.data(ig + 72 * ncomps + c);
        const auto *ig_73 = buffer.data(ig + 73 * ncomps + c);
        const auto *ig_74 = buffer.data(ig + 74 * ncomps + c);
        const auto *ig_75 = buffer.data(ig + 75 * ncomps + c);
        const auto *ig_76 = buffer.data(ig + 76 * ncomps + c);
        const auto *ig_77 = buffer.data(ig + 77 * ncomps + c);
        const auto *ig_78 = buffer.data(ig + 78 * ncomps + c);
        const auto *ig_79 = buffer.data(ig + 79 * ncomps + c);
        const auto *ig_80 = buffer.data(ig + 80 * ncomps + c);
        const auto *ig_81 = buffer.data(ig + 81 * ncomps + c);
        const auto *ig_82 = buffer.data(ig + 82 * ncomps + c);
        const auto *ig_83 = buffer.data(ig + 83 * ncomps + c);
        const auto *ig_84 = buffer.data(ig + 84 * ncomps + c);
        const auto *ig_85 = buffer.data(ig + 85 * ncomps + c);
        const auto *ig_86 = buffer.data(ig + 86 * ncomps + c);
        const auto *ig_87 = buffer.data(ig + 87 * ncomps + c);
        const auto *ig_88 = buffer.data(ig + 88 * ncomps + c);
        const auto *ig_89 = buffer.data(ig + 89 * ncomps + c);
        const auto *ig_90 = buffer.data(ig + 90 * ncomps + c);
        const auto *ig_91 = buffer.data(ig + 91 * ncomps + c);
        const auto *ig_92 = buffer.data(ig + 92 * ncomps + c);
        const auto *ig_93 = buffer.data(ig + 93 * ncomps + c);
        const auto *ig_94 = buffer.data(ig + 94 * ncomps + c);
        const auto *ig_95 = buffer.data(ig + 95 * ncomps + c);
        const auto *ig_96 = buffer.data(ig + 96 * ncomps + c);
        const auto *ig_97 = buffer.data(ig + 97 * ncomps + c);
        const auto *ig_98 = buffer.data(ig + 98 * ncomps + c);
        const auto *ig_99 = buffer.data(ig + 99 * ncomps + c);
        const auto *ig_100 = buffer.data(ig + 100 * ncomps + c);
        const auto *ig_101 = buffer.data(ig + 101 * ncomps + c);
        const auto *ig_102 = buffer.data(ig + 102 * ncomps + c);
        const auto *ig_103 = buffer.data(ig + 103 * ncomps + c);
        const auto *ig_104 = buffer.data(ig + 104 * ncomps + c);

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
        const auto *kg_115 = buffer.data(kg + 115 * ncomps + c);
        const auto *kg_116 = buffer.data(kg + 116 * ncomps + c);
        const auto *kg_117 = buffer.data(kg + 117 * ncomps + c);
        const auto *kg_118 = buffer.data(kg + 118 * ncomps + c);
        const auto *kg_119 = buffer.data(kg + 119 * ncomps + c);
        const auto *kg_130 = buffer.data(kg + 130 * ncomps + c);
        const auto *kg_131 = buffer.data(kg + 131 * ncomps + c);
        const auto *kg_132 = buffer.data(kg + 132 * ncomps + c);
        const auto *kg_133 = buffer.data(kg + 133 * ncomps + c);
        const auto *kg_134 = buffer.data(kg + 134 * ncomps + c);
        const auto *kg_149 = buffer.data(kg + 149 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ig_0, ig_1, ig_2, ig_3, ig_4, kg_0, \
                         kg_1, kg_2, kg_3, kg_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * ig_0[k]
                     + kg_0[k];

            t_1[k] = ab_x[k] * ig_1[k]
                     + kg_1[k];

            t_2[k] = ab_x[k] * ig_2[k]
                     + kg_2[k];

            t_3[k] = ab_x[k] * ig_3[k]
                     + kg_3[k];

            t_4[k] = ab_x[k] * ig_4[k]
                     + kg_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ig_5, ig_6, ig_7, ig_8, ig_9, kg_5, \
                         kg_6, kg_7, kg_8, kg_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * ig_5[k]
                     + kg_5[k];

            t_6[k] = ab_x[k] * ig_6[k]
                     + kg_6[k];

            t_7[k] = ab_x[k] * ig_7[k]
                     + kg_7[k];

            t_8[k] = ab_x[k] * ig_8[k]
                     + kg_8[k];

            t_9[k] = ab_x[k] * ig_9[k]
                     + kg_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, ig_10, ig_11, ig_12, ig_13, \
                         ig_14, kg_10, kg_11, kg_12, kg_13, kg_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_x[k] * ig_10[k]
                      + kg_10[k];

            t_11[k] = ab_x[k] * ig_11[k]
                      + kg_11[k];

            t_12[k] = ab_x[k] * ig_12[k]
                      + kg_12[k];

            t_13[k] = ab_x[k] * ig_13[k]
                      + kg_13[k];

            t_14[k] = ab_x[k] * ig_14[k]
                      + kg_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_y, ig_10, ig_11, ig_12, ig_13, \
                         ig_14, kg_25, kg_26, kg_27, kg_28, kg_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_y[k] * ig_10[k]
                      + kg_25[k];

            t_16[k] = ab_y[k] * ig_11[k]
                      + kg_26[k];

            t_17[k] = ab_y[k] * ig_12[k]
                      + kg_27[k];

            t_18[k] = ab_y[k] * ig_13[k]
                      + kg_28[k];

            t_19[k] = ab_y[k] * ig_14[k]
                      + kg_29[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, ab_x, ab_z, ig_14, ig_15, ig_16, ig_17, \
                         kg_15, kg_16, kg_17, kg_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_z[k] * ig_14[k]
                      + kg_44[k];

            t_21[k] = ab_x[k] * ig_15[k]
                      + kg_15[k];

            t_22[k] = ab_x[k] * ig_16[k]
                      + kg_16[k];

            t_23[k] = ab_x[k] * ig_17[k]
                      + kg_17[k];
        }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, ab_x, ig_18, ig_19, ig_20, ig_21, \
                         ig_22, kg_18, kg_19, kg_20, kg_21, kg_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_24[k] = ab_x[k] * ig_18[k]
                      + kg_18[k];

            t_25[k] = ab_x[k] * ig_19[k]
                      + kg_19[k];

            t_26[k] = ab_x[k] * ig_20[k]
                      + kg_20[k];

            t_27[k] = ab_x[k] * ig_21[k]
                      + kg_21[k];

            t_28[k] = ab_x[k] * ig_22[k]
                      + kg_22[k];
        }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, ab_x, ig_23, ig_24, ig_25, ig_26, \
                         ig_27, kg_23, kg_24, kg_25, kg_26, kg_27 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_29[k] = ab_x[k] * ig_23[k]
                      + kg_23[k];

            t_30[k] = ab_x[k] * ig_24[k]
                      + kg_24[k];

            t_31[k] = ab_x[k] * ig_25[k]
                      + kg_25[k];

            t_32[k] = ab_x[k] * ig_26[k]
                      + kg_26[k];

            t_33[k] = ab_x[k] * ig_27[k]
                      + kg_27[k];
        }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, ab_x, ab_y, ig_25, ig_26, ig_28, ig_29, \
                         kg_28, kg_29, kg_55, kg_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_34[k] = ab_x[k] * ig_28[k]
                      + kg_28[k];

            t_35[k] = ab_x[k] * ig_29[k]
                      + kg_29[k];

            t_36[k] = ab_y[k] * ig_25[k]
                      + kg_55[k];

            t_37[k] = ab_y[k] * ig_26[k]
                      + kg_56[k];
        }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, ab_y, ab_z, ig_27, ig_28, ig_29, kg_57, \
                         kg_58, kg_59, kg_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_38[k] = ab_y[k] * ig_27[k]
                      + kg_57[k];

            t_39[k] = ab_y[k] * ig_28[k]
                      + kg_58[k];

            t_40[k] = ab_y[k] * ig_29[k]
                      + kg_59[k];

            t_41[k] = ab_z[k] * ig_29[k]
                      + kg_74[k];
        }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, ab_x, ig_30, ig_31, ig_32, ig_33, \
                         ig_34, kg_30, kg_31, kg_32, kg_33, kg_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_42[k] = ab_x[k] * ig_30[k]
                      + kg_30[k];

            t_43[k] = ab_x[k] * ig_31[k]
                      + kg_31[k];

            t_44[k] = ab_x[k] * ig_32[k]
                      + kg_32[k];

            t_45[k] = ab_x[k] * ig_33[k]
                      + kg_33[k];

            t_46[k] = ab_x[k] * ig_34[k]
                      + kg_34[k];
        }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, ab_x, ig_35, ig_36, ig_37, ig_38, \
                         ig_39, kg_35, kg_36, kg_37, kg_38, kg_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_47[k] = ab_x[k] * ig_35[k]
                      + kg_35[k];

            t_48[k] = ab_x[k] * ig_36[k]
                      + kg_36[k];

            t_49[k] = ab_x[k] * ig_37[k]
                      + kg_37[k];

            t_50[k] = ab_x[k] * ig_38[k]
                      + kg_38[k];

            t_51[k] = ab_x[k] * ig_39[k]
                      + kg_39[k];
        }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, ab_x, ig_40, ig_41, ig_42, ig_43, \
                         ig_44, kg_40, kg_41, kg_42, kg_43, kg_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_52[k] = ab_x[k] * ig_40[k]
                      + kg_40[k];

            t_53[k] = ab_x[k] * ig_41[k]
                      + kg_41[k];

            t_54[k] = ab_x[k] * ig_42[k]
                      + kg_42[k];

            t_55[k] = ab_x[k] * ig_43[k]
                      + kg_43[k];

            t_56[k] = ab_x[k] * ig_44[k]
                      + kg_44[k];
        }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, ab_y, ig_40, ig_41, ig_42, ig_43, \
                         ig_44, kg_70, kg_71, kg_72, kg_73, kg_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_57[k] = ab_y[k] * ig_40[k]
                      + kg_70[k];

            t_58[k] = ab_y[k] * ig_41[k]
                      + kg_71[k];

            t_59[k] = ab_y[k] * ig_42[k]
                      + kg_72[k];

            t_60[k] = ab_y[k] * ig_43[k]
                      + kg_73[k];

            t_61[k] = ab_y[k] * ig_44[k]
                      + kg_74[k];
        }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, ab_x, ab_z, ig_44, ig_45, ig_46, ig_47, \
                         kg_45, kg_46, kg_47, kg_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_62[k] = ab_z[k] * ig_44[k]
                      + kg_89[k];

            t_63[k] = ab_x[k] * ig_45[k]
                      + kg_45[k];

            t_64[k] = ab_x[k] * ig_46[k]
                      + kg_46[k];

            t_65[k] = ab_x[k] * ig_47[k]
                      + kg_47[k];
        }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, ab_x, ig_48, ig_49, ig_50, ig_51, \
                         ig_52, kg_48, kg_49, kg_50, kg_51, kg_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_66[k] = ab_x[k] * ig_48[k]
                      + kg_48[k];

            t_67[k] = ab_x[k] * ig_49[k]
                      + kg_49[k];

            t_68[k] = ab_x[k] * ig_50[k]
                      + kg_50[k];

            t_69[k] = ab_x[k] * ig_51[k]
                      + kg_51[k];

            t_70[k] = ab_x[k] * ig_52[k]
                      + kg_52[k];
        }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, ab_x, ig_53, ig_54, ig_55, ig_56, \
                         ig_57, kg_53, kg_54, kg_55, kg_56, kg_57 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_71[k] = ab_x[k] * ig_53[k]
                      + kg_53[k];

            t_72[k] = ab_x[k] * ig_54[k]
                      + kg_54[k];

            t_73[k] = ab_x[k] * ig_55[k]
                      + kg_55[k];

            t_74[k] = ab_x[k] * ig_56[k]
                      + kg_56[k];

            t_75[k] = ab_x[k] * ig_57[k]
                      + kg_57[k];
        }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, ab_x, ab_y, ig_55, ig_56, ig_58, ig_59, \
                         kg_58, kg_59, kg_100, kg_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_76[k] = ab_x[k] * ig_58[k]
                      + kg_58[k];

            t_77[k] = ab_x[k] * ig_59[k]
                      + kg_59[k];

            t_78[k] = ab_y[k] * ig_55[k]
                      + kg_100[k];

            t_79[k] = ab_y[k] * ig_56[k]
                      + kg_101[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, ab_y, ab_z, ig_57, ig_58, ig_59, kg_102, \
                         kg_103, kg_104, kg_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_y[k] * ig_57[k]
                      + kg_102[k];

            t_81[k] = ab_y[k] * ig_58[k]
                      + kg_103[k];

            t_82[k] = ab_y[k] * ig_59[k]
                      + kg_104[k];

            t_83[k] = ab_z[k] * ig_59[k]
                      + kg_119[k];
        }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, ab_x, ig_60, ig_61, ig_62, ig_63, \
                         ig_64, kg_60, kg_61, kg_62, kg_63, kg_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_84[k] = ab_x[k] * ig_60[k]
                      + kg_60[k];

            t_85[k] = ab_x[k] * ig_61[k]
                      + kg_61[k];

            t_86[k] = ab_x[k] * ig_62[k]
                      + kg_62[k];

            t_87[k] = ab_x[k] * ig_63[k]
                      + kg_63[k];

            t_88[k] = ab_x[k] * ig_64[k]
                      + kg_64[k];
        }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, ab_x, ig_65, ig_66, ig_67, ig_68, \
                         ig_69, kg_65, kg_66, kg_67, kg_68, kg_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_89[k] = ab_x[k] * ig_65[k]
                      + kg_65[k];

            t_90[k] = ab_x[k] * ig_66[k]
                      + kg_66[k];

            t_91[k] = ab_x[k] * ig_67[k]
                      + kg_67[k];

            t_92[k] = ab_x[k] * ig_68[k]
                      + kg_68[k];

            t_93[k] = ab_x[k] * ig_69[k]
                      + kg_69[k];
        }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, ab_x, ig_70, ig_71, ig_72, ig_73, \
                         ig_74, kg_70, kg_71, kg_72, kg_73, kg_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_94[k] = ab_x[k] * ig_70[k]
                      + kg_70[k];

            t_95[k] = ab_x[k] * ig_71[k]
                      + kg_71[k];

            t_96[k] = ab_x[k] * ig_72[k]
                      + kg_72[k];

            t_97[k] = ab_x[k] * ig_73[k]
                      + kg_73[k];

            t_98[k] = ab_x[k] * ig_74[k]
                      + kg_74[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, ab_y, ig_70, ig_71, ig_72, ig_73, \
                         ig_74, kg_115, kg_116, kg_117, kg_118, \
                         kg_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_y[k] * ig_70[k]
                      + kg_115[k];

            t_100[k] = ab_y[k] * ig_71[k]
                       + kg_116[k];

            t_101[k] = ab_y[k] * ig_72[k]
                       + kg_117[k];

            t_102[k] = ab_y[k] * ig_73[k]
                       + kg_118[k];

            t_103[k] = ab_y[k] * ig_74[k]
                       + kg_119[k];
        }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, ab_x, ab_z, ig_74, ig_75, ig_76, ig_77, \
                         kg_75, kg_76, kg_77, kg_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_104[k] = ab_z[k] * ig_74[k]
                       + kg_134[k];

            t_105[k] = ab_x[k] * ig_75[k]
                       + kg_75[k];

            t_106[k] = ab_x[k] * ig_76[k]
                       + kg_76[k];

            t_107[k] = ab_x[k] * ig_77[k]
                       + kg_77[k];
        }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, ab_x, ig_78, ig_79, ig_80, ig_81, \
                         ig_82, kg_78, kg_79, kg_80, kg_81, kg_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_108[k] = ab_x[k] * ig_78[k]
                       + kg_78[k];

            t_109[k] = ab_x[k] * ig_79[k]
                       + kg_79[k];

            t_110[k] = ab_x[k] * ig_80[k]
                       + kg_80[k];

            t_111[k] = ab_x[k] * ig_81[k]
                       + kg_81[k];

            t_112[k] = ab_x[k] * ig_82[k]
                       + kg_82[k];
        }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, ab_x, ig_83, ig_84, ig_85, ig_86, \
                         ig_87, kg_83, kg_84, kg_85, kg_86, kg_87 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_113[k] = ab_x[k] * ig_83[k]
                       + kg_83[k];

            t_114[k] = ab_x[k] * ig_84[k]
                       + kg_84[k];

            t_115[k] = ab_x[k] * ig_85[k]
                       + kg_85[k];

            t_116[k] = ab_x[k] * ig_86[k]
                       + kg_86[k];

            t_117[k] = ab_x[k] * ig_87[k]
                       + kg_87[k];
        }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, ab_x, ab_y, ig_85, ig_86, ig_88, ig_89, \
                         kg_88, kg_89, kg_130, kg_131 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_118[k] = ab_x[k] * ig_88[k]
                       + kg_88[k];

            t_119[k] = ab_x[k] * ig_89[k]
                       + kg_89[k];

            t_120[k] = ab_y[k] * ig_85[k]
                       + kg_130[k];

            t_121[k] = ab_y[k] * ig_86[k]
                       + kg_131[k];
        }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, ab_y, ab_z, ig_87, ig_88, ig_89, kg_132, \
                         kg_133, kg_134, kg_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_122[k] = ab_y[k] * ig_87[k]
                       + kg_132[k];

            t_123[k] = ab_y[k] * ig_88[k]
                       + kg_133[k];

            t_124[k] = ab_y[k] * ig_89[k]
                       + kg_134[k];

            t_125[k] = ab_z[k] * ig_89[k]
                       + kg_149[k];
        }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, ab_x, ig_90, ig_91, ig_92, ig_93, \
                         ig_94, kg_90, kg_91, kg_92, kg_93, kg_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_126[k] = ab_x[k] * ig_90[k]
                       + kg_90[k];

            t_127[k] = ab_x[k] * ig_91[k]
                       + kg_91[k];

            t_128[k] = ab_x[k] * ig_92[k]
                       + kg_92[k];

            t_129[k] = ab_x[k] * ig_93[k]
                       + kg_93[k];

            t_130[k] = ab_x[k] * ig_94[k]
                       + kg_94[k];
        }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, ab_x, ig_95, ig_96, ig_97, ig_98, \
                         ig_99, kg_95, kg_96, kg_97, kg_98, kg_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_131[k] = ab_x[k] * ig_95[k]
                       + kg_95[k];

            t_132[k] = ab_x[k] * ig_96[k]
                       + kg_96[k];

            t_133[k] = ab_x[k] * ig_97[k]
                       + kg_97[k];

            t_134[k] = ab_x[k] * ig_98[k]
                       + kg_98[k];

            t_135[k] = ab_x[k] * ig_99[k]
                       + kg_99[k];
        }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, ab_x, ig_100, ig_101, ig_102, \
                         ig_103, ig_104, kg_100, kg_101, kg_102, kg_103, \
                         kg_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_136[k] = ab_x[k] * ig_100[k]
                       + kg_100[k];

            t_137[k] = ab_x[k] * ig_101[k]
                       + kg_101[k];

            t_138[k] = ab_x[k] * ig_102[k]
                       + kg_102[k];

            t_139[k] = ab_x[k] * ig_103[k]
                       + kg_103[k];

            t_140[k] = ab_x[k] * ig_104[k]
                       + kg_104[k];
        }
    }
}

static auto
compute_hrr_ih_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ig, const size_t kg,
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

        const auto *ig_100 = buffer.data(ig + 100 * ncomps + c);
        const auto *ig_101 = buffer.data(ig + 101 * ncomps + c);
        const auto *ig_102 = buffer.data(ig + 102 * ncomps + c);
        const auto *ig_103 = buffer.data(ig + 103 * ncomps + c);
        const auto *ig_104 = buffer.data(ig + 104 * ncomps + c);
        const auto *ig_105 = buffer.data(ig + 105 * ncomps + c);
        const auto *ig_106 = buffer.data(ig + 106 * ncomps + c);
        const auto *ig_107 = buffer.data(ig + 107 * ncomps + c);
        const auto *ig_108 = buffer.data(ig + 108 * ncomps + c);
        const auto *ig_109 = buffer.data(ig + 109 * ncomps + c);
        const auto *ig_110 = buffer.data(ig + 110 * ncomps + c);
        const auto *ig_111 = buffer.data(ig + 111 * ncomps + c);
        const auto *ig_112 = buffer.data(ig + 112 * ncomps + c);
        const auto *ig_113 = buffer.data(ig + 113 * ncomps + c);
        const auto *ig_114 = buffer.data(ig + 114 * ncomps + c);
        const auto *ig_115 = buffer.data(ig + 115 * ncomps + c);
        const auto *ig_116 = buffer.data(ig + 116 * ncomps + c);
        const auto *ig_117 = buffer.data(ig + 117 * ncomps + c);
        const auto *ig_118 = buffer.data(ig + 118 * ncomps + c);
        const auto *ig_119 = buffer.data(ig + 119 * ncomps + c);
        const auto *ig_120 = buffer.data(ig + 120 * ncomps + c);
        const auto *ig_121 = buffer.data(ig + 121 * ncomps + c);
        const auto *ig_122 = buffer.data(ig + 122 * ncomps + c);
        const auto *ig_123 = buffer.data(ig + 123 * ncomps + c);
        const auto *ig_124 = buffer.data(ig + 124 * ncomps + c);
        const auto *ig_125 = buffer.data(ig + 125 * ncomps + c);
        const auto *ig_126 = buffer.data(ig + 126 * ncomps + c);
        const auto *ig_127 = buffer.data(ig + 127 * ncomps + c);
        const auto *ig_128 = buffer.data(ig + 128 * ncomps + c);
        const auto *ig_129 = buffer.data(ig + 129 * ncomps + c);
        const auto *ig_130 = buffer.data(ig + 130 * ncomps + c);
        const auto *ig_131 = buffer.data(ig + 131 * ncomps + c);
        const auto *ig_132 = buffer.data(ig + 132 * ncomps + c);
        const auto *ig_133 = buffer.data(ig + 133 * ncomps + c);
        const auto *ig_134 = buffer.data(ig + 134 * ncomps + c);
        const auto *ig_135 = buffer.data(ig + 135 * ncomps + c);
        const auto *ig_136 = buffer.data(ig + 136 * ncomps + c);
        const auto *ig_137 = buffer.data(ig + 137 * ncomps + c);
        const auto *ig_138 = buffer.data(ig + 138 * ncomps + c);
        const auto *ig_139 = buffer.data(ig + 139 * ncomps + c);
        const auto *ig_140 = buffer.data(ig + 140 * ncomps + c);
        const auto *ig_141 = buffer.data(ig + 141 * ncomps + c);
        const auto *ig_142 = buffer.data(ig + 142 * ncomps + c);
        const auto *ig_143 = buffer.data(ig + 143 * ncomps + c);
        const auto *ig_144 = buffer.data(ig + 144 * ncomps + c);
        const auto *ig_145 = buffer.data(ig + 145 * ncomps + c);
        const auto *ig_146 = buffer.data(ig + 146 * ncomps + c);
        const auto *ig_147 = buffer.data(ig + 147 * ncomps + c);
        const auto *ig_148 = buffer.data(ig + 148 * ncomps + c);
        const auto *ig_149 = buffer.data(ig + 149 * ncomps + c);
        const auto *ig_150 = buffer.data(ig + 150 * ncomps + c);
        const auto *ig_151 = buffer.data(ig + 151 * ncomps + c);
        const auto *ig_152 = buffer.data(ig + 152 * ncomps + c);
        const auto *ig_153 = buffer.data(ig + 153 * ncomps + c);
        const auto *ig_154 = buffer.data(ig + 154 * ncomps + c);
        const auto *ig_155 = buffer.data(ig + 155 * ncomps + c);
        const auto *ig_156 = buffer.data(ig + 156 * ncomps + c);
        const auto *ig_157 = buffer.data(ig + 157 * ncomps + c);
        const auto *ig_158 = buffer.data(ig + 158 * ncomps + c);
        const auto *ig_159 = buffer.data(ig + 159 * ncomps + c);
        const auto *ig_160 = buffer.data(ig + 160 * ncomps + c);
        const auto *ig_161 = buffer.data(ig + 161 * ncomps + c);
        const auto *ig_162 = buffer.data(ig + 162 * ncomps + c);
        const auto *ig_163 = buffer.data(ig + 163 * ncomps + c);
        const auto *ig_164 = buffer.data(ig + 164 * ncomps + c);
        const auto *ig_165 = buffer.data(ig + 165 * ncomps + c);
        const auto *ig_166 = buffer.data(ig + 166 * ncomps + c);
        const auto *ig_167 = buffer.data(ig + 167 * ncomps + c);
        const auto *ig_168 = buffer.data(ig + 168 * ncomps + c);
        const auto *ig_169 = buffer.data(ig + 169 * ncomps + c);
        const auto *ig_170 = buffer.data(ig + 170 * ncomps + c);
        const auto *ig_171 = buffer.data(ig + 171 * ncomps + c);
        const auto *ig_172 = buffer.data(ig + 172 * ncomps + c);
        const auto *ig_173 = buffer.data(ig + 173 * ncomps + c);
        const auto *ig_174 = buffer.data(ig + 174 * ncomps + c);
        const auto *ig_175 = buffer.data(ig + 175 * ncomps + c);
        const auto *ig_176 = buffer.data(ig + 176 * ncomps + c);
        const auto *ig_177 = buffer.data(ig + 177 * ncomps + c);
        const auto *ig_178 = buffer.data(ig + 178 * ncomps + c);
        const auto *ig_179 = buffer.data(ig + 179 * ncomps + c);
        const auto *ig_180 = buffer.data(ig + 180 * ncomps + c);
        const auto *ig_181 = buffer.data(ig + 181 * ncomps + c);
        const auto *ig_182 = buffer.data(ig + 182 * ncomps + c);
        const auto *ig_183 = buffer.data(ig + 183 * ncomps + c);
        const auto *ig_184 = buffer.data(ig + 184 * ncomps + c);
        const auto *ig_185 = buffer.data(ig + 185 * ncomps + c);
        const auto *ig_186 = buffer.data(ig + 186 * ncomps + c);
        const auto *ig_187 = buffer.data(ig + 187 * ncomps + c);
        const auto *ig_188 = buffer.data(ig + 188 * ncomps + c);
        const auto *ig_189 = buffer.data(ig + 189 * ncomps + c);
        const auto *ig_190 = buffer.data(ig + 190 * ncomps + c);
        const auto *ig_191 = buffer.data(ig + 191 * ncomps + c);
        const auto *ig_192 = buffer.data(ig + 192 * ncomps + c);
        const auto *ig_193 = buffer.data(ig + 193 * ncomps + c);
        const auto *ig_194 = buffer.data(ig + 194 * ncomps + c);
        const auto *ig_195 = buffer.data(ig + 195 * ncomps + c);
        const auto *ig_196 = buffer.data(ig + 196 * ncomps + c);
        const auto *ig_197 = buffer.data(ig + 197 * ncomps + c);
        const auto *ig_198 = buffer.data(ig + 198 * ncomps + c);
        const auto *ig_199 = buffer.data(ig + 199 * ncomps + c);
        const auto *ig_200 = buffer.data(ig + 200 * ncomps + c);
        const auto *ig_201 = buffer.data(ig + 201 * ncomps + c);
        const auto *ig_202 = buffer.data(ig + 202 * ncomps + c);

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
        const auto *kg_205 = buffer.data(kg + 205 * ncomps + c);
        const auto *kg_206 = buffer.data(kg + 206 * ncomps + c);
        const auto *kg_207 = buffer.data(kg + 207 * ncomps + c);
        const auto *kg_208 = buffer.data(kg + 208 * ncomps + c);
        const auto *kg_209 = buffer.data(kg + 209 * ncomps + c);
        const auto *kg_224 = buffer.data(kg + 224 * ncomps + c);
        const auto *kg_235 = buffer.data(kg + 235 * ncomps + c);
        const auto *kg_236 = buffer.data(kg + 236 * ncomps + c);
        const auto *kg_237 = buffer.data(kg + 237 * ncomps + c);
        const auto *kg_238 = buffer.data(kg + 238 * ncomps + c);
        const auto *kg_239 = buffer.data(kg + 239 * ncomps + c);
        const auto *kg_250 = buffer.data(kg + 250 * ncomps + c);
        const auto *kg_251 = buffer.data(kg + 251 * ncomps + c);
        const auto *kg_252 = buffer.data(kg + 252 * ncomps + c);
        const auto *kg_253 = buffer.data(kg + 253 * ncomps + c);
        const auto *kg_254 = buffer.data(kg + 254 * ncomps + c);
        const auto *kg_265 = buffer.data(kg + 265 * ncomps + c);
        const auto *kg_266 = buffer.data(kg + 266 * ncomps + c);
        const auto *kg_267 = buffer.data(kg + 267 * ncomps + c);
        const auto *kg_268 = buffer.data(kg + 268 * ncomps + c);
        const auto *kg_269 = buffer.data(kg + 269 * ncomps + c);
        const auto *kg_284 = buffer.data(kg + 284 * ncomps + c);

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, ab_y, ig_100, ig_101, ig_102, \
                         ig_103, ig_104, kg_160, kg_161, kg_162, kg_163, \
                         kg_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_141[k] = ab_y[k] * ig_100[k]
                       + kg_160[k];

            t_142[k] = ab_y[k] * ig_101[k]
                       + kg_161[k];

            t_143[k] = ab_y[k] * ig_102[k]
                       + kg_162[k];

            t_144[k] = ab_y[k] * ig_103[k]
                       + kg_163[k];

            t_145[k] = ab_y[k] * ig_104[k]
                       + kg_164[k];
        }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, ab_x, ab_z, ig_104, ig_105, ig_106, \
                         ig_107, kg_105, kg_106, kg_107, kg_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_146[k] = ab_z[k] * ig_104[k]
                       + kg_179[k];

            t_147[k] = ab_x[k] * ig_105[k]
                       + kg_105[k];

            t_148[k] = ab_x[k] * ig_106[k]
                       + kg_106[k];

            t_149[k] = ab_x[k] * ig_107[k]
                       + kg_107[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, ig_108, ig_109, ig_110, \
                         ig_111, ig_112, kg_108, kg_109, kg_110, kg_111, \
                         kg_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * ig_108[k]
                       + kg_108[k];

            t_151[k] = ab_x[k] * ig_109[k]
                       + kg_109[k];

            t_152[k] = ab_x[k] * ig_110[k]
                       + kg_110[k];

            t_153[k] = ab_x[k] * ig_111[k]
                       + kg_111[k];

            t_154[k] = ab_x[k] * ig_112[k]
                       + kg_112[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, ig_113, ig_114, ig_115, \
                         ig_116, ig_117, kg_113, kg_114, kg_115, kg_116, \
                         kg_117 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * ig_113[k]
                       + kg_113[k];

            t_156[k] = ab_x[k] * ig_114[k]
                       + kg_114[k];

            t_157[k] = ab_x[k] * ig_115[k]
                       + kg_115[k];

            t_158[k] = ab_x[k] * ig_116[k]
                       + kg_116[k];

            t_159[k] = ab_x[k] * ig_117[k]
                       + kg_117[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, ab_x, ab_y, ig_115, ig_116, ig_118, \
                         ig_119, kg_118, kg_119, kg_175, kg_176 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_x[k] * ig_118[k]
                       + kg_118[k];

            t_161[k] = ab_x[k] * ig_119[k]
                       + kg_119[k];

            t_162[k] = ab_y[k] * ig_115[k]
                       + kg_175[k];

            t_163[k] = ab_y[k] * ig_116[k]
                       + kg_176[k];
        }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, ab_y, ab_z, ig_117, ig_118, ig_119, \
                         kg_177, kg_178, kg_179, kg_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_164[k] = ab_y[k] * ig_117[k]
                       + kg_177[k];

            t_165[k] = ab_y[k] * ig_118[k]
                       + kg_178[k];

            t_166[k] = ab_y[k] * ig_119[k]
                       + kg_179[k];

            t_167[k] = ab_z[k] * ig_119[k]
                       + kg_194[k];
        }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, ab_x, ig_120, ig_121, ig_122, \
                         ig_123, ig_124, kg_120, kg_121, kg_122, kg_123, \
                         kg_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_168[k] = ab_x[k] * ig_120[k]
                       + kg_120[k];

            t_169[k] = ab_x[k] * ig_121[k]
                       + kg_121[k];

            t_170[k] = ab_x[k] * ig_122[k]
                       + kg_122[k];

            t_171[k] = ab_x[k] * ig_123[k]
                       + kg_123[k];

            t_172[k] = ab_x[k] * ig_124[k]
                       + kg_124[k];
        }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, ab_x, ig_125, ig_126, ig_127, \
                         ig_128, ig_129, kg_125, kg_126, kg_127, kg_128, \
                         kg_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_173[k] = ab_x[k] * ig_125[k]
                       + kg_125[k];

            t_174[k] = ab_x[k] * ig_126[k]
                       + kg_126[k];

            t_175[k] = ab_x[k] * ig_127[k]
                       + kg_127[k];

            t_176[k] = ab_x[k] * ig_128[k]
                       + kg_128[k];

            t_177[k] = ab_x[k] * ig_129[k]
                       + kg_129[k];
        }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, ab_x, ig_130, ig_131, ig_132, \
                         ig_133, ig_134, kg_130, kg_131, kg_132, kg_133, \
                         kg_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_178[k] = ab_x[k] * ig_130[k]
                       + kg_130[k];

            t_179[k] = ab_x[k] * ig_131[k]
                       + kg_131[k];

            t_180[k] = ab_x[k] * ig_132[k]
                       + kg_132[k];

            t_181[k] = ab_x[k] * ig_133[k]
                       + kg_133[k];

            t_182[k] = ab_x[k] * ig_134[k]
                       + kg_134[k];
        }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, ab_y, ig_130, ig_131, ig_132, \
                         ig_133, ig_134, kg_190, kg_191, kg_192, kg_193, \
                         kg_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_183[k] = ab_y[k] * ig_130[k]
                       + kg_190[k];

            t_184[k] = ab_y[k] * ig_131[k]
                       + kg_191[k];

            t_185[k] = ab_y[k] * ig_132[k]
                       + kg_192[k];

            t_186[k] = ab_y[k] * ig_133[k]
                       + kg_193[k];

            t_187[k] = ab_y[k] * ig_134[k]
                       + kg_194[k];
        }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, ab_x, ab_z, ig_134, ig_135, ig_136, \
                         ig_137, kg_135, kg_136, kg_137, kg_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_188[k] = ab_z[k] * ig_134[k]
                       + kg_209[k];

            t_189[k] = ab_x[k] * ig_135[k]
                       + kg_135[k];

            t_190[k] = ab_x[k] * ig_136[k]
                       + kg_136[k];

            t_191[k] = ab_x[k] * ig_137[k]
                       + kg_137[k];
        }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, ab_x, ig_138, ig_139, ig_140, \
                         ig_141, ig_142, kg_138, kg_139, kg_140, kg_141, \
                         kg_142 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_192[k] = ab_x[k] * ig_138[k]
                       + kg_138[k];

            t_193[k] = ab_x[k] * ig_139[k]
                       + kg_139[k];

            t_194[k] = ab_x[k] * ig_140[k]
                       + kg_140[k];

            t_195[k] = ab_x[k] * ig_141[k]
                       + kg_141[k];

            t_196[k] = ab_x[k] * ig_142[k]
                       + kg_142[k];
        }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, ab_x, ig_143, ig_144, ig_145, \
                         ig_146, ig_147, kg_143, kg_144, kg_145, kg_146, \
                         kg_147 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_197[k] = ab_x[k] * ig_143[k]
                       + kg_143[k];

            t_198[k] = ab_x[k] * ig_144[k]
                       + kg_144[k];

            t_199[k] = ab_x[k] * ig_145[k]
                       + kg_145[k];

            t_200[k] = ab_x[k] * ig_146[k]
                       + kg_146[k];

            t_201[k] = ab_x[k] * ig_147[k]
                       + kg_147[k];
        }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, ab_x, ab_y, ig_145, ig_146, ig_148, \
                         ig_149, kg_148, kg_149, kg_205, kg_206 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_202[k] = ab_x[k] * ig_148[k]
                       + kg_148[k];

            t_203[k] = ab_x[k] * ig_149[k]
                       + kg_149[k];

            t_204[k] = ab_y[k] * ig_145[k]
                       + kg_205[k];

            t_205[k] = ab_y[k] * ig_146[k]
                       + kg_206[k];
        }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, ab_y, ab_z, ig_147, ig_148, ig_149, \
                         kg_207, kg_208, kg_209, kg_224 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_206[k] = ab_y[k] * ig_147[k]
                       + kg_207[k];

            t_207[k] = ab_y[k] * ig_148[k]
                       + kg_208[k];

            t_208[k] = ab_y[k] * ig_149[k]
                       + kg_209[k];

            t_209[k] = ab_z[k] * ig_149[k]
                       + kg_224[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, ig_150, ig_151, ig_152, \
                         ig_153, ig_154, kg_150, kg_151, kg_152, kg_153, \
                         kg_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = ab_x[k] * ig_150[k]
                       + kg_150[k];

            t_211[k] = ab_x[k] * ig_151[k]
                       + kg_151[k];

            t_212[k] = ab_x[k] * ig_152[k]
                       + kg_152[k];

            t_213[k] = ab_x[k] * ig_153[k]
                       + kg_153[k];

            t_214[k] = ab_x[k] * ig_154[k]
                       + kg_154[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, ig_155, ig_156, ig_157, \
                         ig_158, ig_159, kg_155, kg_156, kg_157, kg_158, \
                         kg_159 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = ab_x[k] * ig_155[k]
                       + kg_155[k];

            t_216[k] = ab_x[k] * ig_156[k]
                       + kg_156[k];

            t_217[k] = ab_x[k] * ig_157[k]
                       + kg_157[k];

            t_218[k] = ab_x[k] * ig_158[k]
                       + kg_158[k];

            t_219[k] = ab_x[k] * ig_159[k]
                       + kg_159[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, ig_160, ig_161, ig_162, \
                         ig_163, ig_164, kg_160, kg_161, kg_162, kg_163, \
                         kg_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = ab_x[k] * ig_160[k]
                       + kg_160[k];

            t_221[k] = ab_x[k] * ig_161[k]
                       + kg_161[k];

            t_222[k] = ab_x[k] * ig_162[k]
                       + kg_162[k];

            t_223[k] = ab_x[k] * ig_163[k]
                       + kg_163[k];

            t_224[k] = ab_x[k] * ig_164[k]
                       + kg_164[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_y, ig_160, ig_161, ig_162, \
                         ig_163, ig_164, kg_235, kg_236, kg_237, kg_238, \
                         kg_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = ab_y[k] * ig_160[k]
                       + kg_235[k];

            t_226[k] = ab_y[k] * ig_161[k]
                       + kg_236[k];

            t_227[k] = ab_y[k] * ig_162[k]
                       + kg_237[k];

            t_228[k] = ab_y[k] * ig_163[k]
                       + kg_238[k];

            t_229[k] = ab_y[k] * ig_164[k]
                       + kg_239[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, ab_x, ab_z, ig_164, ig_165, ig_166, \
                         ig_167, kg_165, kg_166, kg_167, kg_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = ab_z[k] * ig_164[k]
                       + kg_254[k];

            t_231[k] = ab_x[k] * ig_165[k]
                       + kg_165[k];

            t_232[k] = ab_x[k] * ig_166[k]
                       + kg_166[k];

            t_233[k] = ab_x[k] * ig_167[k]
                       + kg_167[k];
        }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, ab_x, ig_168, ig_169, ig_170, \
                         ig_171, ig_172, kg_168, kg_169, kg_170, kg_171, \
                         kg_172 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_234[k] = ab_x[k] * ig_168[k]
                       + kg_168[k];

            t_235[k] = ab_x[k] * ig_169[k]
                       + kg_169[k];

            t_236[k] = ab_x[k] * ig_170[k]
                       + kg_170[k];

            t_237[k] = ab_x[k] * ig_171[k]
                       + kg_171[k];

            t_238[k] = ab_x[k] * ig_172[k]
                       + kg_172[k];
        }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, ab_x, ig_173, ig_174, ig_175, \
                         ig_176, ig_177, kg_173, kg_174, kg_175, kg_176, \
                         kg_177 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_239[k] = ab_x[k] * ig_173[k]
                       + kg_173[k];

            t_240[k] = ab_x[k] * ig_174[k]
                       + kg_174[k];

            t_241[k] = ab_x[k] * ig_175[k]
                       + kg_175[k];

            t_242[k] = ab_x[k] * ig_176[k]
                       + kg_176[k];

            t_243[k] = ab_x[k] * ig_177[k]
                       + kg_177[k];
        }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, ab_x, ab_y, ig_175, ig_176, ig_178, \
                         ig_179, kg_178, kg_179, kg_250, kg_251 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_244[k] = ab_x[k] * ig_178[k]
                       + kg_178[k];

            t_245[k] = ab_x[k] * ig_179[k]
                       + kg_179[k];

            t_246[k] = ab_y[k] * ig_175[k]
                       + kg_250[k];

            t_247[k] = ab_y[k] * ig_176[k]
                       + kg_251[k];
        }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, ab_y, ab_z, ig_177, ig_178, ig_179, \
                         kg_252, kg_253, kg_254, kg_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_248[k] = ab_y[k] * ig_177[k]
                       + kg_252[k];

            t_249[k] = ab_y[k] * ig_178[k]
                       + kg_253[k];

            t_250[k] = ab_y[k] * ig_179[k]
                       + kg_254[k];

            t_251[k] = ab_z[k] * ig_179[k]
                       + kg_269[k];
        }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, ab_x, ig_180, ig_181, ig_182, \
                         ig_183, ig_184, kg_180, kg_181, kg_182, kg_183, \
                         kg_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_252[k] = ab_x[k] * ig_180[k]
                       + kg_180[k];

            t_253[k] = ab_x[k] * ig_181[k]
                       + kg_181[k];

            t_254[k] = ab_x[k] * ig_182[k]
                       + kg_182[k];

            t_255[k] = ab_x[k] * ig_183[k]
                       + kg_183[k];

            t_256[k] = ab_x[k] * ig_184[k]
                       + kg_184[k];
        }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, ab_x, ig_185, ig_186, ig_187, \
                         ig_188, ig_189, kg_185, kg_186, kg_187, kg_188, \
                         kg_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_257[k] = ab_x[k] * ig_185[k]
                       + kg_185[k];

            t_258[k] = ab_x[k] * ig_186[k]
                       + kg_186[k];

            t_259[k] = ab_x[k] * ig_187[k]
                       + kg_187[k];

            t_260[k] = ab_x[k] * ig_188[k]
                       + kg_188[k];

            t_261[k] = ab_x[k] * ig_189[k]
                       + kg_189[k];
        }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, ab_x, ig_190, ig_191, ig_192, \
                         ig_193, ig_194, kg_190, kg_191, kg_192, kg_193, \
                         kg_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_262[k] = ab_x[k] * ig_190[k]
                       + kg_190[k];

            t_263[k] = ab_x[k] * ig_191[k]
                       + kg_191[k];

            t_264[k] = ab_x[k] * ig_192[k]
                       + kg_192[k];

            t_265[k] = ab_x[k] * ig_193[k]
                       + kg_193[k];

            t_266[k] = ab_x[k] * ig_194[k]
                       + kg_194[k];
        }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, ab_y, ig_190, ig_191, ig_192, \
                         ig_193, ig_194, kg_265, kg_266, kg_267, kg_268, \
                         kg_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_267[k] = ab_y[k] * ig_190[k]
                       + kg_265[k];

            t_268[k] = ab_y[k] * ig_191[k]
                       + kg_266[k];

            t_269[k] = ab_y[k] * ig_192[k]
                       + kg_267[k];

            t_270[k] = ab_y[k] * ig_193[k]
                       + kg_268[k];

            t_271[k] = ab_y[k] * ig_194[k]
                       + kg_269[k];
        }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, ab_x, ab_z, ig_194, ig_195, ig_196, \
                         ig_197, kg_195, kg_196, kg_197, kg_284 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_272[k] = ab_z[k] * ig_194[k]
                       + kg_284[k];

            t_273[k] = ab_x[k] * ig_195[k]
                       + kg_195[k];

            t_274[k] = ab_x[k] * ig_196[k]
                       + kg_196[k];

            t_275[k] = ab_x[k] * ig_197[k]
                       + kg_197[k];
        }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, ab_x, ig_198, ig_199, ig_200, \
                         ig_201, ig_202, kg_198, kg_199, kg_200, kg_201, \
                         kg_202 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_276[k] = ab_x[k] * ig_198[k]
                       + kg_198[k];

            t_277[k] = ab_x[k] * ig_199[k]
                       + kg_199[k];

            t_278[k] = ab_x[k] * ig_200[k]
                       + kg_200[k];

            t_279[k] = ab_x[k] * ig_201[k]
                       + kg_201[k];

            t_280[k] = ab_x[k] * ig_202[k]
                       + kg_202[k];
        }
    }
}

static auto
compute_hrr_ih_out_of_first_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ig, const size_t kg,
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

        const auto *ig_203 = buffer.data(ig + 203 * ncomps + c);
        const auto *ig_204 = buffer.data(ig + 204 * ncomps + c);
        const auto *ig_205 = buffer.data(ig + 205 * ncomps + c);
        const auto *ig_206 = buffer.data(ig + 206 * ncomps + c);
        const auto *ig_207 = buffer.data(ig + 207 * ncomps + c);
        const auto *ig_208 = buffer.data(ig + 208 * ncomps + c);
        const auto *ig_209 = buffer.data(ig + 209 * ncomps + c);
        const auto *ig_210 = buffer.data(ig + 210 * ncomps + c);
        const auto *ig_211 = buffer.data(ig + 211 * ncomps + c);
        const auto *ig_212 = buffer.data(ig + 212 * ncomps + c);
        const auto *ig_213 = buffer.data(ig + 213 * ncomps + c);
        const auto *ig_214 = buffer.data(ig + 214 * ncomps + c);
        const auto *ig_215 = buffer.data(ig + 215 * ncomps + c);
        const auto *ig_216 = buffer.data(ig + 216 * ncomps + c);
        const auto *ig_217 = buffer.data(ig + 217 * ncomps + c);
        const auto *ig_218 = buffer.data(ig + 218 * ncomps + c);
        const auto *ig_219 = buffer.data(ig + 219 * ncomps + c);
        const auto *ig_220 = buffer.data(ig + 220 * ncomps + c);
        const auto *ig_221 = buffer.data(ig + 221 * ncomps + c);
        const auto *ig_222 = buffer.data(ig + 222 * ncomps + c);
        const auto *ig_223 = buffer.data(ig + 223 * ncomps + c);
        const auto *ig_224 = buffer.data(ig + 224 * ncomps + c);
        const auto *ig_225 = buffer.data(ig + 225 * ncomps + c);
        const auto *ig_226 = buffer.data(ig + 226 * ncomps + c);
        const auto *ig_227 = buffer.data(ig + 227 * ncomps + c);
        const auto *ig_228 = buffer.data(ig + 228 * ncomps + c);
        const auto *ig_229 = buffer.data(ig + 229 * ncomps + c);
        const auto *ig_230 = buffer.data(ig + 230 * ncomps + c);
        const auto *ig_231 = buffer.data(ig + 231 * ncomps + c);
        const auto *ig_232 = buffer.data(ig + 232 * ncomps + c);
        const auto *ig_233 = buffer.data(ig + 233 * ncomps + c);
        const auto *ig_234 = buffer.data(ig + 234 * ncomps + c);
        const auto *ig_235 = buffer.data(ig + 235 * ncomps + c);
        const auto *ig_236 = buffer.data(ig + 236 * ncomps + c);
        const auto *ig_237 = buffer.data(ig + 237 * ncomps + c);
        const auto *ig_238 = buffer.data(ig + 238 * ncomps + c);
        const auto *ig_239 = buffer.data(ig + 239 * ncomps + c);
        const auto *ig_240 = buffer.data(ig + 240 * ncomps + c);
        const auto *ig_241 = buffer.data(ig + 241 * ncomps + c);
        const auto *ig_242 = buffer.data(ig + 242 * ncomps + c);
        const auto *ig_243 = buffer.data(ig + 243 * ncomps + c);
        const auto *ig_244 = buffer.data(ig + 244 * ncomps + c);
        const auto *ig_245 = buffer.data(ig + 245 * ncomps + c);
        const auto *ig_246 = buffer.data(ig + 246 * ncomps + c);
        const auto *ig_247 = buffer.data(ig + 247 * ncomps + c);
        const auto *ig_248 = buffer.data(ig + 248 * ncomps + c);
        const auto *ig_249 = buffer.data(ig + 249 * ncomps + c);
        const auto *ig_250 = buffer.data(ig + 250 * ncomps + c);
        const auto *ig_251 = buffer.data(ig + 251 * ncomps + c);
        const auto *ig_252 = buffer.data(ig + 252 * ncomps + c);
        const auto *ig_253 = buffer.data(ig + 253 * ncomps + c);
        const auto *ig_254 = buffer.data(ig + 254 * ncomps + c);
        const auto *ig_255 = buffer.data(ig + 255 * ncomps + c);
        const auto *ig_256 = buffer.data(ig + 256 * ncomps + c);
        const auto *ig_257 = buffer.data(ig + 257 * ncomps + c);
        const auto *ig_258 = buffer.data(ig + 258 * ncomps + c);
        const auto *ig_259 = buffer.data(ig + 259 * ncomps + c);
        const auto *ig_260 = buffer.data(ig + 260 * ncomps + c);
        const auto *ig_261 = buffer.data(ig + 261 * ncomps + c);
        const auto *ig_262 = buffer.data(ig + 262 * ncomps + c);
        const auto *ig_263 = buffer.data(ig + 263 * ncomps + c);
        const auto *ig_264 = buffer.data(ig + 264 * ncomps + c);
        const auto *ig_265 = buffer.data(ig + 265 * ncomps + c);
        const auto *ig_266 = buffer.data(ig + 266 * ncomps + c);
        const auto *ig_267 = buffer.data(ig + 267 * ncomps + c);
        const auto *ig_268 = buffer.data(ig + 268 * ncomps + c);
        const auto *ig_269 = buffer.data(ig + 269 * ncomps + c);
        const auto *ig_270 = buffer.data(ig + 270 * ncomps + c);
        const auto *ig_271 = buffer.data(ig + 271 * ncomps + c);
        const auto *ig_272 = buffer.data(ig + 272 * ncomps + c);
        const auto *ig_273 = buffer.data(ig + 273 * ncomps + c);
        const auto *ig_274 = buffer.data(ig + 274 * ncomps + c);
        const auto *ig_275 = buffer.data(ig + 275 * ncomps + c);
        const auto *ig_276 = buffer.data(ig + 276 * ncomps + c);
        const auto *ig_277 = buffer.data(ig + 277 * ncomps + c);
        const auto *ig_278 = buffer.data(ig + 278 * ncomps + c);
        const auto *ig_279 = buffer.data(ig + 279 * ncomps + c);
        const auto *ig_280 = buffer.data(ig + 280 * ncomps + c);
        const auto *ig_281 = buffer.data(ig + 281 * ncomps + c);
        const auto *ig_282 = buffer.data(ig + 282 * ncomps + c);
        const auto *ig_283 = buffer.data(ig + 283 * ncomps + c);
        const auto *ig_284 = buffer.data(ig + 284 * ncomps + c);
        const auto *ig_285 = buffer.data(ig + 285 * ncomps + c);
        const auto *ig_286 = buffer.data(ig + 286 * ncomps + c);
        const auto *ig_287 = buffer.data(ig + 287 * ncomps + c);
        const auto *ig_288 = buffer.data(ig + 288 * ncomps + c);
        const auto *ig_289 = buffer.data(ig + 289 * ncomps + c);
        const auto *ig_290 = buffer.data(ig + 290 * ncomps + c);
        const auto *ig_291 = buffer.data(ig + 291 * ncomps + c);
        const auto *ig_292 = buffer.data(ig + 292 * ncomps + c);
        const auto *ig_293 = buffer.data(ig + 293 * ncomps + c);
        const auto *ig_294 = buffer.data(ig + 294 * ncomps + c);
        const auto *ig_295 = buffer.data(ig + 295 * ncomps + c);
        const auto *ig_296 = buffer.data(ig + 296 * ncomps + c);
        const auto *ig_297 = buffer.data(ig + 297 * ncomps + c);
        const auto *ig_298 = buffer.data(ig + 298 * ncomps + c);
        const auto *ig_299 = buffer.data(ig + 299 * ncomps + c);
        const auto *ig_300 = buffer.data(ig + 300 * ncomps + c);
        const auto *ig_301 = buffer.data(ig + 301 * ncomps + c);
        const auto *ig_302 = buffer.data(ig + 302 * ncomps + c);
        const auto *ig_303 = buffer.data(ig + 303 * ncomps + c);
        const auto *ig_304 = buffer.data(ig + 304 * ncomps + c);

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
        const auto *kg_314 = buffer.data(kg + 314 * ncomps + c);
        const auto *kg_325 = buffer.data(kg + 325 * ncomps + c);
        const auto *kg_326 = buffer.data(kg + 326 * ncomps + c);
        const auto *kg_327 = buffer.data(kg + 327 * ncomps + c);
        const auto *kg_328 = buffer.data(kg + 328 * ncomps + c);
        const auto *kg_329 = buffer.data(kg + 329 * ncomps + c);
        const auto *kg_340 = buffer.data(kg + 340 * ncomps + c);
        const auto *kg_341 = buffer.data(kg + 341 * ncomps + c);
        const auto *kg_342 = buffer.data(kg + 342 * ncomps + c);
        const auto *kg_343 = buffer.data(kg + 343 * ncomps + c);
        const auto *kg_344 = buffer.data(kg + 344 * ncomps + c);
        const auto *kg_355 = buffer.data(kg + 355 * ncomps + c);
        const auto *kg_356 = buffer.data(kg + 356 * ncomps + c);
        const auto *kg_357 = buffer.data(kg + 357 * ncomps + c);
        const auto *kg_358 = buffer.data(kg + 358 * ncomps + c);
        const auto *kg_359 = buffer.data(kg + 359 * ncomps + c);
        const auto *kg_370 = buffer.data(kg + 370 * ncomps + c);
        const auto *kg_371 = buffer.data(kg + 371 * ncomps + c);
        const auto *kg_372 = buffer.data(kg + 372 * ncomps + c);
        const auto *kg_373 = buffer.data(kg + 373 * ncomps + c);
        const auto *kg_374 = buffer.data(kg + 374 * ncomps + c);
        const auto *kg_385 = buffer.data(kg + 385 * ncomps + c);
        const auto *kg_386 = buffer.data(kg + 386 * ncomps + c);
        const auto *kg_387 = buffer.data(kg + 387 * ncomps + c);
        const auto *kg_388 = buffer.data(kg + 388 * ncomps + c);
        const auto *kg_389 = buffer.data(kg + 389 * ncomps + c);
        const auto *kg_404 = buffer.data(kg + 404 * ncomps + c);

#pragma omp simd aligned(t_281, t_282, t_283, t_284, t_285, ab_x, ig_203, ig_204, ig_205, \
                         ig_206, ig_207, kg_203, kg_204, kg_205, kg_206, \
                         kg_207 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_281[k] = ab_x[k] * ig_203[k]
                       + kg_203[k];

            t_282[k] = ab_x[k] * ig_204[k]
                       + kg_204[k];

            t_283[k] = ab_x[k] * ig_205[k]
                       + kg_205[k];

            t_284[k] = ab_x[k] * ig_206[k]
                       + kg_206[k];

            t_285[k] = ab_x[k] * ig_207[k]
                       + kg_207[k];
        }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, ab_x, ab_y, ig_205, ig_206, ig_208, \
                         ig_209, kg_208, kg_209, kg_280, kg_281 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_286[k] = ab_x[k] * ig_208[k]
                       + kg_208[k];

            t_287[k] = ab_x[k] * ig_209[k]
                       + kg_209[k];

            t_288[k] = ab_y[k] * ig_205[k]
                       + kg_280[k];

            t_289[k] = ab_y[k] * ig_206[k]
                       + kg_281[k];
        }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, ab_y, ab_z, ig_207, ig_208, ig_209, \
                         kg_282, kg_283, kg_284, kg_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = ab_y[k] * ig_207[k]
                       + kg_282[k];

            t_291[k] = ab_y[k] * ig_208[k]
                       + kg_283[k];

            t_292[k] = ab_y[k] * ig_209[k]
                       + kg_284[k];

            t_293[k] = ab_z[k] * ig_209[k]
                       + kg_299[k];
        }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, ab_x, ig_210, ig_211, ig_212, \
                         ig_213, ig_214, kg_210, kg_211, kg_212, kg_213, \
                         kg_214 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_294[k] = ab_x[k] * ig_210[k]
                       + kg_210[k];

            t_295[k] = ab_x[k] * ig_211[k]
                       + kg_211[k];

            t_296[k] = ab_x[k] * ig_212[k]
                       + kg_212[k];

            t_297[k] = ab_x[k] * ig_213[k]
                       + kg_213[k];

            t_298[k] = ab_x[k] * ig_214[k]
                       + kg_214[k];
        }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, ab_x, ig_215, ig_216, ig_217, \
                         ig_218, ig_219, kg_215, kg_216, kg_217, kg_218, \
                         kg_219 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_299[k] = ab_x[k] * ig_215[k]
                       + kg_215[k];

            t_300[k] = ab_x[k] * ig_216[k]
                       + kg_216[k];

            t_301[k] = ab_x[k] * ig_217[k]
                       + kg_217[k];

            t_302[k] = ab_x[k] * ig_218[k]
                       + kg_218[k];

            t_303[k] = ab_x[k] * ig_219[k]
                       + kg_219[k];
        }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, ab_x, ig_220, ig_221, ig_222, \
                         ig_223, ig_224, kg_220, kg_221, kg_222, kg_223, \
                         kg_224 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_304[k] = ab_x[k] * ig_220[k]
                       + kg_220[k];

            t_305[k] = ab_x[k] * ig_221[k]
                       + kg_221[k];

            t_306[k] = ab_x[k] * ig_222[k]
                       + kg_222[k];

            t_307[k] = ab_x[k] * ig_223[k]
                       + kg_223[k];

            t_308[k] = ab_x[k] * ig_224[k]
                       + kg_224[k];
        }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, ab_y, ig_220, ig_221, ig_222, \
                         ig_223, ig_224, kg_295, kg_296, kg_297, kg_298, \
                         kg_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_309[k] = ab_y[k] * ig_220[k]
                       + kg_295[k];

            t_310[k] = ab_y[k] * ig_221[k]
                       + kg_296[k];

            t_311[k] = ab_y[k] * ig_222[k]
                       + kg_297[k];

            t_312[k] = ab_y[k] * ig_223[k]
                       + kg_298[k];

            t_313[k] = ab_y[k] * ig_224[k]
                       + kg_299[k];
        }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, ab_x, ab_z, ig_224, ig_225, ig_226, \
                         ig_227, kg_225, kg_226, kg_227, kg_314 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_314[k] = ab_z[k] * ig_224[k]
                       + kg_314[k];

            t_315[k] = ab_x[k] * ig_225[k]
                       + kg_225[k];

            t_316[k] = ab_x[k] * ig_226[k]
                       + kg_226[k];

            t_317[k] = ab_x[k] * ig_227[k]
                       + kg_227[k];
        }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, ab_x, ig_228, ig_229, ig_230, \
                         ig_231, ig_232, kg_228, kg_229, kg_230, kg_231, \
                         kg_232 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_318[k] = ab_x[k] * ig_228[k]
                       + kg_228[k];

            t_319[k] = ab_x[k] * ig_229[k]
                       + kg_229[k];

            t_320[k] = ab_x[k] * ig_230[k]
                       + kg_230[k];

            t_321[k] = ab_x[k] * ig_231[k]
                       + kg_231[k];

            t_322[k] = ab_x[k] * ig_232[k]
                       + kg_232[k];
        }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, ab_x, ig_233, ig_234, ig_235, \
                         ig_236, ig_237, kg_233, kg_234, kg_235, kg_236, \
                         kg_237 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_323[k] = ab_x[k] * ig_233[k]
                       + kg_233[k];

            t_324[k] = ab_x[k] * ig_234[k]
                       + kg_234[k];

            t_325[k] = ab_x[k] * ig_235[k]
                       + kg_235[k];

            t_326[k] = ab_x[k] * ig_236[k]
                       + kg_236[k];

            t_327[k] = ab_x[k] * ig_237[k]
                       + kg_237[k];
        }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, ab_x, ab_y, ig_235, ig_236, ig_238, \
                         ig_239, kg_238, kg_239, kg_325, kg_326 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_328[k] = ab_x[k] * ig_238[k]
                       + kg_238[k];

            t_329[k] = ab_x[k] * ig_239[k]
                       + kg_239[k];

            t_330[k] = ab_y[k] * ig_235[k]
                       + kg_325[k];

            t_331[k] = ab_y[k] * ig_236[k]
                       + kg_326[k];
        }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, ab_y, ab_z, ig_237, ig_238, ig_239, \
                         kg_327, kg_328, kg_329, kg_344 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_332[k] = ab_y[k] * ig_237[k]
                       + kg_327[k];

            t_333[k] = ab_y[k] * ig_238[k]
                       + kg_328[k];

            t_334[k] = ab_y[k] * ig_239[k]
                       + kg_329[k];

            t_335[k] = ab_z[k] * ig_239[k]
                       + kg_344[k];
        }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, t_340, ab_x, ig_240, ig_241, ig_242, \
                         ig_243, ig_244, kg_240, kg_241, kg_242, kg_243, \
                         kg_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_336[k] = ab_x[k] * ig_240[k]
                       + kg_240[k];

            t_337[k] = ab_x[k] * ig_241[k]
                       + kg_241[k];

            t_338[k] = ab_x[k] * ig_242[k]
                       + kg_242[k];

            t_339[k] = ab_x[k] * ig_243[k]
                       + kg_243[k];

            t_340[k] = ab_x[k] * ig_244[k]
                       + kg_244[k];
        }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, t_345, ab_x, ig_245, ig_246, ig_247, \
                         ig_248, ig_249, kg_245, kg_246, kg_247, kg_248, \
                         kg_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_341[k] = ab_x[k] * ig_245[k]
                       + kg_245[k];

            t_342[k] = ab_x[k] * ig_246[k]
                       + kg_246[k];

            t_343[k] = ab_x[k] * ig_247[k]
                       + kg_247[k];

            t_344[k] = ab_x[k] * ig_248[k]
                       + kg_248[k];

            t_345[k] = ab_x[k] * ig_249[k]
                       + kg_249[k];
        }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, ab_x, ig_250, ig_251, ig_252, \
                         ig_253, ig_254, kg_250, kg_251, kg_252, kg_253, \
                         kg_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_346[k] = ab_x[k] * ig_250[k]
                       + kg_250[k];

            t_347[k] = ab_x[k] * ig_251[k]
                       + kg_251[k];

            t_348[k] = ab_x[k] * ig_252[k]
                       + kg_252[k];

            t_349[k] = ab_x[k] * ig_253[k]
                       + kg_253[k];

            t_350[k] = ab_x[k] * ig_254[k]
                       + kg_254[k];
        }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, ab_y, ig_250, ig_251, ig_252, \
                         ig_253, ig_254, kg_340, kg_341, kg_342, kg_343, \
                         kg_344 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_351[k] = ab_y[k] * ig_250[k]
                       + kg_340[k];

            t_352[k] = ab_y[k] * ig_251[k]
                       + kg_341[k];

            t_353[k] = ab_y[k] * ig_252[k]
                       + kg_342[k];

            t_354[k] = ab_y[k] * ig_253[k]
                       + kg_343[k];

            t_355[k] = ab_y[k] * ig_254[k]
                       + kg_344[k];
        }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, ab_x, ab_z, ig_254, ig_255, ig_256, \
                         ig_257, kg_255, kg_256, kg_257, kg_359 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_356[k] = ab_z[k] * ig_254[k]
                       + kg_359[k];

            t_357[k] = ab_x[k] * ig_255[k]
                       + kg_255[k];

            t_358[k] = ab_x[k] * ig_256[k]
                       + kg_256[k];

            t_359[k] = ab_x[k] * ig_257[k]
                       + kg_257[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, ig_258, ig_259, ig_260, \
                         ig_261, ig_262, kg_258, kg_259, kg_260, kg_261, \
                         kg_262 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = ab_x[k] * ig_258[k]
                       + kg_258[k];

            t_361[k] = ab_x[k] * ig_259[k]
                       + kg_259[k];

            t_362[k] = ab_x[k] * ig_260[k]
                       + kg_260[k];

            t_363[k] = ab_x[k] * ig_261[k]
                       + kg_261[k];

            t_364[k] = ab_x[k] * ig_262[k]
                       + kg_262[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, ig_263, ig_264, ig_265, \
                         ig_266, ig_267, kg_263, kg_264, kg_265, kg_266, \
                         kg_267 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = ab_x[k] * ig_263[k]
                       + kg_263[k];

            t_366[k] = ab_x[k] * ig_264[k]
                       + kg_264[k];

            t_367[k] = ab_x[k] * ig_265[k]
                       + kg_265[k];

            t_368[k] = ab_x[k] * ig_266[k]
                       + kg_266[k];

            t_369[k] = ab_x[k] * ig_267[k]
                       + kg_267[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, ab_x, ab_y, ig_265, ig_266, ig_268, \
                         ig_269, kg_268, kg_269, kg_355, kg_356 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = ab_x[k] * ig_268[k]
                       + kg_268[k];

            t_371[k] = ab_x[k] * ig_269[k]
                       + kg_269[k];

            t_372[k] = ab_y[k] * ig_265[k]
                       + kg_355[k];

            t_373[k] = ab_y[k] * ig_266[k]
                       + kg_356[k];
        }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, ab_y, ab_z, ig_267, ig_268, ig_269, \
                         kg_357, kg_358, kg_359, kg_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_374[k] = ab_y[k] * ig_267[k]
                       + kg_357[k];

            t_375[k] = ab_y[k] * ig_268[k]
                       + kg_358[k];

            t_376[k] = ab_y[k] * ig_269[k]
                       + kg_359[k];

            t_377[k] = ab_z[k] * ig_269[k]
                       + kg_374[k];
        }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, t_382, ab_x, ig_270, ig_271, ig_272, \
                         ig_273, ig_274, kg_270, kg_271, kg_272, kg_273, \
                         kg_274 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_378[k] = ab_x[k] * ig_270[k]
                       + kg_270[k];

            t_379[k] = ab_x[k] * ig_271[k]
                       + kg_271[k];

            t_380[k] = ab_x[k] * ig_272[k]
                       + kg_272[k];

            t_381[k] = ab_x[k] * ig_273[k]
                       + kg_273[k];

            t_382[k] = ab_x[k] * ig_274[k]
                       + kg_274[k];
        }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, ab_x, ig_275, ig_276, ig_277, \
                         ig_278, ig_279, kg_275, kg_276, kg_277, kg_278, \
                         kg_279 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_383[k] = ab_x[k] * ig_275[k]
                       + kg_275[k];

            t_384[k] = ab_x[k] * ig_276[k]
                       + kg_276[k];

            t_385[k] = ab_x[k] * ig_277[k]
                       + kg_277[k];

            t_386[k] = ab_x[k] * ig_278[k]
                       + kg_278[k];

            t_387[k] = ab_x[k] * ig_279[k]
                       + kg_279[k];
        }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, t_392, ab_x, ig_280, ig_281, ig_282, \
                         ig_283, ig_284, kg_280, kg_281, kg_282, kg_283, \
                         kg_284 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_388[k] = ab_x[k] * ig_280[k]
                       + kg_280[k];

            t_389[k] = ab_x[k] * ig_281[k]
                       + kg_281[k];

            t_390[k] = ab_x[k] * ig_282[k]
                       + kg_282[k];

            t_391[k] = ab_x[k] * ig_283[k]
                       + kg_283[k];

            t_392[k] = ab_x[k] * ig_284[k]
                       + kg_284[k];
        }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, t_397, ab_y, ig_280, ig_281, ig_282, \
                         ig_283, ig_284, kg_370, kg_371, kg_372, kg_373, \
                         kg_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_393[k] = ab_y[k] * ig_280[k]
                       + kg_370[k];

            t_394[k] = ab_y[k] * ig_281[k]
                       + kg_371[k];

            t_395[k] = ab_y[k] * ig_282[k]
                       + kg_372[k];

            t_396[k] = ab_y[k] * ig_283[k]
                       + kg_373[k];

            t_397[k] = ab_y[k] * ig_284[k]
                       + kg_374[k];
        }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, ab_x, ab_z, ig_284, ig_285, ig_286, \
                         ig_287, kg_285, kg_286, kg_287, kg_389 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_398[k] = ab_z[k] * ig_284[k]
                       + kg_389[k];

            t_399[k] = ab_x[k] * ig_285[k]
                       + kg_285[k];

            t_400[k] = ab_x[k] * ig_286[k]
                       + kg_286[k];

            t_401[k] = ab_x[k] * ig_287[k]
                       + kg_287[k];
        }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, ab_x, ig_288, ig_289, ig_290, \
                         ig_291, ig_292, kg_288, kg_289, kg_290, kg_291, \
                         kg_292 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_402[k] = ab_x[k] * ig_288[k]
                       + kg_288[k];

            t_403[k] = ab_x[k] * ig_289[k]
                       + kg_289[k];

            t_404[k] = ab_x[k] * ig_290[k]
                       + kg_290[k];

            t_405[k] = ab_x[k] * ig_291[k]
                       + kg_291[k];

            t_406[k] = ab_x[k] * ig_292[k]
                       + kg_292[k];
        }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, ab_x, ig_293, ig_294, ig_295, \
                         ig_296, ig_297, kg_293, kg_294, kg_295, kg_296, \
                         kg_297 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_407[k] = ab_x[k] * ig_293[k]
                       + kg_293[k];

            t_408[k] = ab_x[k] * ig_294[k]
                       + kg_294[k];

            t_409[k] = ab_x[k] * ig_295[k]
                       + kg_295[k];

            t_410[k] = ab_x[k] * ig_296[k]
                       + kg_296[k];

            t_411[k] = ab_x[k] * ig_297[k]
                       + kg_297[k];
        }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, ab_x, ab_y, ig_295, ig_296, ig_298, \
                         ig_299, kg_298, kg_299, kg_385, kg_386 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_412[k] = ab_x[k] * ig_298[k]
                       + kg_298[k];

            t_413[k] = ab_x[k] * ig_299[k]
                       + kg_299[k];

            t_414[k] = ab_y[k] * ig_295[k]
                       + kg_385[k];

            t_415[k] = ab_y[k] * ig_296[k]
                       + kg_386[k];
        }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, ab_y, ab_z, ig_297, ig_298, ig_299, \
                         kg_387, kg_388, kg_389, kg_404 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_416[k] = ab_y[k] * ig_297[k]
                       + kg_387[k];

            t_417[k] = ab_y[k] * ig_298[k]
                       + kg_388[k];

            t_418[k] = ab_y[k] * ig_299[k]
                       + kg_389[k];

            t_419[k] = ab_z[k] * ig_299[k]
                       + kg_404[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_x, ig_300, ig_301, ig_302, \
                         ig_303, ig_304, kg_300, kg_301, kg_302, kg_303, \
                         kg_304 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = ab_x[k] * ig_300[k]
                       + kg_300[k];

            t_421[k] = ab_x[k] * ig_301[k]
                       + kg_301[k];

            t_422[k] = ab_x[k] * ig_302[k]
                       + kg_302[k];

            t_423[k] = ab_x[k] * ig_303[k]
                       + kg_303[k];

            t_424[k] = ab_x[k] * ig_304[k]
                       + kg_304[k];
        }
    }
}

static auto
compute_hrr_ih_out_of_first_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ig, const size_t kg,
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

        const auto *ig_305 = buffer.data(ig + 305 * ncomps + c);
        const auto *ig_306 = buffer.data(ig + 306 * ncomps + c);
        const auto *ig_307 = buffer.data(ig + 307 * ncomps + c);
        const auto *ig_308 = buffer.data(ig + 308 * ncomps + c);
        const auto *ig_309 = buffer.data(ig + 309 * ncomps + c);
        const auto *ig_310 = buffer.data(ig + 310 * ncomps + c);
        const auto *ig_311 = buffer.data(ig + 311 * ncomps + c);
        const auto *ig_312 = buffer.data(ig + 312 * ncomps + c);
        const auto *ig_313 = buffer.data(ig + 313 * ncomps + c);
        const auto *ig_314 = buffer.data(ig + 314 * ncomps + c);
        const auto *ig_315 = buffer.data(ig + 315 * ncomps + c);
        const auto *ig_316 = buffer.data(ig + 316 * ncomps + c);
        const auto *ig_317 = buffer.data(ig + 317 * ncomps + c);
        const auto *ig_318 = buffer.data(ig + 318 * ncomps + c);
        const auto *ig_319 = buffer.data(ig + 319 * ncomps + c);
        const auto *ig_320 = buffer.data(ig + 320 * ncomps + c);
        const auto *ig_321 = buffer.data(ig + 321 * ncomps + c);
        const auto *ig_322 = buffer.data(ig + 322 * ncomps + c);
        const auto *ig_323 = buffer.data(ig + 323 * ncomps + c);
        const auto *ig_324 = buffer.data(ig + 324 * ncomps + c);
        const auto *ig_325 = buffer.data(ig + 325 * ncomps + c);
        const auto *ig_326 = buffer.data(ig + 326 * ncomps + c);
        const auto *ig_327 = buffer.data(ig + 327 * ncomps + c);
        const auto *ig_328 = buffer.data(ig + 328 * ncomps + c);
        const auto *ig_329 = buffer.data(ig + 329 * ncomps + c);
        const auto *ig_330 = buffer.data(ig + 330 * ncomps + c);
        const auto *ig_331 = buffer.data(ig + 331 * ncomps + c);
        const auto *ig_332 = buffer.data(ig + 332 * ncomps + c);
        const auto *ig_333 = buffer.data(ig + 333 * ncomps + c);
        const auto *ig_334 = buffer.data(ig + 334 * ncomps + c);
        const auto *ig_335 = buffer.data(ig + 335 * ncomps + c);
        const auto *ig_336 = buffer.data(ig + 336 * ncomps + c);
        const auto *ig_337 = buffer.data(ig + 337 * ncomps + c);
        const auto *ig_338 = buffer.data(ig + 338 * ncomps + c);
        const auto *ig_339 = buffer.data(ig + 339 * ncomps + c);
        const auto *ig_340 = buffer.data(ig + 340 * ncomps + c);
        const auto *ig_341 = buffer.data(ig + 341 * ncomps + c);
        const auto *ig_342 = buffer.data(ig + 342 * ncomps + c);
        const auto *ig_343 = buffer.data(ig + 343 * ncomps + c);
        const auto *ig_344 = buffer.data(ig + 344 * ncomps + c);
        const auto *ig_345 = buffer.data(ig + 345 * ncomps + c);
        const auto *ig_346 = buffer.data(ig + 346 * ncomps + c);
        const auto *ig_347 = buffer.data(ig + 347 * ncomps + c);
        const auto *ig_348 = buffer.data(ig + 348 * ncomps + c);
        const auto *ig_349 = buffer.data(ig + 349 * ncomps + c);
        const auto *ig_350 = buffer.data(ig + 350 * ncomps + c);
        const auto *ig_351 = buffer.data(ig + 351 * ncomps + c);
        const auto *ig_352 = buffer.data(ig + 352 * ncomps + c);
        const auto *ig_353 = buffer.data(ig + 353 * ncomps + c);
        const auto *ig_354 = buffer.data(ig + 354 * ncomps + c);
        const auto *ig_355 = buffer.data(ig + 355 * ncomps + c);
        const auto *ig_356 = buffer.data(ig + 356 * ncomps + c);
        const auto *ig_357 = buffer.data(ig + 357 * ncomps + c);
        const auto *ig_358 = buffer.data(ig + 358 * ncomps + c);
        const auto *ig_359 = buffer.data(ig + 359 * ncomps + c);
        const auto *ig_360 = buffer.data(ig + 360 * ncomps + c);
        const auto *ig_361 = buffer.data(ig + 361 * ncomps + c);
        const auto *ig_362 = buffer.data(ig + 362 * ncomps + c);
        const auto *ig_363 = buffer.data(ig + 363 * ncomps + c);
        const auto *ig_364 = buffer.data(ig + 364 * ncomps + c);
        const auto *ig_365 = buffer.data(ig + 365 * ncomps + c);
        const auto *ig_366 = buffer.data(ig + 366 * ncomps + c);
        const auto *ig_367 = buffer.data(ig + 367 * ncomps + c);
        const auto *ig_368 = buffer.data(ig + 368 * ncomps + c);
        const auto *ig_369 = buffer.data(ig + 369 * ncomps + c);
        const auto *ig_370 = buffer.data(ig + 370 * ncomps + c);
        const auto *ig_371 = buffer.data(ig + 371 * ncomps + c);
        const auto *ig_372 = buffer.data(ig + 372 * ncomps + c);
        const auto *ig_373 = buffer.data(ig + 373 * ncomps + c);
        const auto *ig_374 = buffer.data(ig + 374 * ncomps + c);
        const auto *ig_375 = buffer.data(ig + 375 * ncomps + c);
        const auto *ig_376 = buffer.data(ig + 376 * ncomps + c);
        const auto *ig_377 = buffer.data(ig + 377 * ncomps + c);
        const auto *ig_378 = buffer.data(ig + 378 * ncomps + c);
        const auto *ig_379 = buffer.data(ig + 379 * ncomps + c);
        const auto *ig_380 = buffer.data(ig + 380 * ncomps + c);
        const auto *ig_381 = buffer.data(ig + 381 * ncomps + c);
        const auto *ig_382 = buffer.data(ig + 382 * ncomps + c);
        const auto *ig_383 = buffer.data(ig + 383 * ncomps + c);
        const auto *ig_384 = buffer.data(ig + 384 * ncomps + c);
        const auto *ig_385 = buffer.data(ig + 385 * ncomps + c);
        const auto *ig_386 = buffer.data(ig + 386 * ncomps + c);
        const auto *ig_387 = buffer.data(ig + 387 * ncomps + c);
        const auto *ig_388 = buffer.data(ig + 388 * ncomps + c);
        const auto *ig_389 = buffer.data(ig + 389 * ncomps + c);
        const auto *ig_390 = buffer.data(ig + 390 * ncomps + c);
        const auto *ig_391 = buffer.data(ig + 391 * ncomps + c);
        const auto *ig_392 = buffer.data(ig + 392 * ncomps + c);
        const auto *ig_393 = buffer.data(ig + 393 * ncomps + c);
        const auto *ig_394 = buffer.data(ig + 394 * ncomps + c);
        const auto *ig_395 = buffer.data(ig + 395 * ncomps + c);
        const auto *ig_396 = buffer.data(ig + 396 * ncomps + c);
        const auto *ig_397 = buffer.data(ig + 397 * ncomps + c);
        const auto *ig_398 = buffer.data(ig + 398 * ncomps + c);
        const auto *ig_399 = buffer.data(ig + 399 * ncomps + c);
        const auto *ig_400 = buffer.data(ig + 400 * ncomps + c);
        const auto *ig_401 = buffer.data(ig + 401 * ncomps + c);
        const auto *ig_402 = buffer.data(ig + 402 * ncomps + c);
        const auto *ig_403 = buffer.data(ig + 403 * ncomps + c);
        const auto *ig_404 = buffer.data(ig + 404 * ncomps + c);

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
        const auto *kg_419 = buffer.data(kg + 419 * ncomps + c);
        const auto *kg_430 = buffer.data(kg + 430 * ncomps + c);
        const auto *kg_431 = buffer.data(kg + 431 * ncomps + c);
        const auto *kg_432 = buffer.data(kg + 432 * ncomps + c);
        const auto *kg_433 = buffer.data(kg + 433 * ncomps + c);
        const auto *kg_434 = buffer.data(kg + 434 * ncomps + c);
        const auto *kg_445 = buffer.data(kg + 445 * ncomps + c);
        const auto *kg_446 = buffer.data(kg + 446 * ncomps + c);
        const auto *kg_447 = buffer.data(kg + 447 * ncomps + c);
        const auto *kg_448 = buffer.data(kg + 448 * ncomps + c);
        const auto *kg_449 = buffer.data(kg + 449 * ncomps + c);
        const auto *kg_460 = buffer.data(kg + 460 * ncomps + c);
        const auto *kg_461 = buffer.data(kg + 461 * ncomps + c);
        const auto *kg_462 = buffer.data(kg + 462 * ncomps + c);
        const auto *kg_463 = buffer.data(kg + 463 * ncomps + c);
        const auto *kg_464 = buffer.data(kg + 464 * ncomps + c);
        const auto *kg_475 = buffer.data(kg + 475 * ncomps + c);
        const auto *kg_476 = buffer.data(kg + 476 * ncomps + c);
        const auto *kg_477 = buffer.data(kg + 477 * ncomps + c);
        const auto *kg_478 = buffer.data(kg + 478 * ncomps + c);
        const auto *kg_479 = buffer.data(kg + 479 * ncomps + c);
        const auto *kg_490 = buffer.data(kg + 490 * ncomps + c);
        const auto *kg_491 = buffer.data(kg + 491 * ncomps + c);
        const auto *kg_492 = buffer.data(kg + 492 * ncomps + c);
        const auto *kg_493 = buffer.data(kg + 493 * ncomps + c);
        const auto *kg_494 = buffer.data(kg + 494 * ncomps + c);
        const auto *kg_505 = buffer.data(kg + 505 * ncomps + c);
        const auto *kg_506 = buffer.data(kg + 506 * ncomps + c);
        const auto *kg_507 = buffer.data(kg + 507 * ncomps + c);
        const auto *kg_508 = buffer.data(kg + 508 * ncomps + c);
        const auto *kg_509 = buffer.data(kg + 509 * ncomps + c);

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_x, ig_305, ig_306, ig_307, \
                         ig_308, ig_309, kg_305, kg_306, kg_307, kg_308, \
                         kg_309 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = ab_x[k] * ig_305[k]
                       + kg_305[k];

            t_426[k] = ab_x[k] * ig_306[k]
                       + kg_306[k];

            t_427[k] = ab_x[k] * ig_307[k]
                       + kg_307[k];

            t_428[k] = ab_x[k] * ig_308[k]
                       + kg_308[k];

            t_429[k] = ab_x[k] * ig_309[k]
                       + kg_309[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_x, ig_310, ig_311, ig_312, \
                         ig_313, ig_314, kg_310, kg_311, kg_312, kg_313, \
                         kg_314 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = ab_x[k] * ig_310[k]
                       + kg_310[k];

            t_431[k] = ab_x[k] * ig_311[k]
                       + kg_311[k];

            t_432[k] = ab_x[k] * ig_312[k]
                       + kg_312[k];

            t_433[k] = ab_x[k] * ig_313[k]
                       + kg_313[k];

            t_434[k] = ab_x[k] * ig_314[k]
                       + kg_314[k];
        }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_y, ig_310, ig_311, ig_312, \
                         ig_313, ig_314, kg_400, kg_401, kg_402, kg_403, \
                         kg_404 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = ab_y[k] * ig_310[k]
                       + kg_400[k];

            t_436[k] = ab_y[k] * ig_311[k]
                       + kg_401[k];

            t_437[k] = ab_y[k] * ig_312[k]
                       + kg_402[k];

            t_438[k] = ab_y[k] * ig_313[k]
                       + kg_403[k];

            t_439[k] = ab_y[k] * ig_314[k]
                       + kg_404[k];
        }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, ab_x, ab_z, ig_314, ig_315, ig_316, \
                         ig_317, kg_315, kg_316, kg_317, kg_419 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = ab_z[k] * ig_314[k]
                       + kg_419[k];

            t_441[k] = ab_x[k] * ig_315[k]
                       + kg_315[k];

            t_442[k] = ab_x[k] * ig_316[k]
                       + kg_316[k];

            t_443[k] = ab_x[k] * ig_317[k]
                       + kg_317[k];
        }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, t_448, ab_x, ig_318, ig_319, ig_320, \
                         ig_321, ig_322, kg_318, kg_319, kg_320, kg_321, \
                         kg_322 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_444[k] = ab_x[k] * ig_318[k]
                       + kg_318[k];

            t_445[k] = ab_x[k] * ig_319[k]
                       + kg_319[k];

            t_446[k] = ab_x[k] * ig_320[k]
                       + kg_320[k];

            t_447[k] = ab_x[k] * ig_321[k]
                       + kg_321[k];

            t_448[k] = ab_x[k] * ig_322[k]
                       + kg_322[k];
        }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, t_453, ab_x, ig_323, ig_324, ig_325, \
                         ig_326, ig_327, kg_323, kg_324, kg_325, kg_326, \
                         kg_327 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_449[k] = ab_x[k] * ig_323[k]
                       + kg_323[k];

            t_450[k] = ab_x[k] * ig_324[k]
                       + kg_324[k];

            t_451[k] = ab_x[k] * ig_325[k]
                       + kg_325[k];

            t_452[k] = ab_x[k] * ig_326[k]
                       + kg_326[k];

            t_453[k] = ab_x[k] * ig_327[k]
                       + kg_327[k];
        }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, ab_x, ab_y, ig_325, ig_326, ig_328, \
                         ig_329, kg_328, kg_329, kg_430, kg_431 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_454[k] = ab_x[k] * ig_328[k]
                       + kg_328[k];

            t_455[k] = ab_x[k] * ig_329[k]
                       + kg_329[k];

            t_456[k] = ab_y[k] * ig_325[k]
                       + kg_430[k];

            t_457[k] = ab_y[k] * ig_326[k]
                       + kg_431[k];
        }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, ab_y, ab_z, ig_327, ig_328, ig_329, \
                         kg_432, kg_433, kg_434, kg_449 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_458[k] = ab_y[k] * ig_327[k]
                       + kg_432[k];

            t_459[k] = ab_y[k] * ig_328[k]
                       + kg_433[k];

            t_460[k] = ab_y[k] * ig_329[k]
                       + kg_434[k];

            t_461[k] = ab_z[k] * ig_329[k]
                       + kg_449[k];
        }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, t_466, ab_x, ig_330, ig_331, ig_332, \
                         ig_333, ig_334, kg_330, kg_331, kg_332, kg_333, \
                         kg_334 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_462[k] = ab_x[k] * ig_330[k]
                       + kg_330[k];

            t_463[k] = ab_x[k] * ig_331[k]
                       + kg_331[k];

            t_464[k] = ab_x[k] * ig_332[k]
                       + kg_332[k];

            t_465[k] = ab_x[k] * ig_333[k]
                       + kg_333[k];

            t_466[k] = ab_x[k] * ig_334[k]
                       + kg_334[k];
        }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, t_471, ab_x, ig_335, ig_336, ig_337, \
                         ig_338, ig_339, kg_335, kg_336, kg_337, kg_338, \
                         kg_339 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_467[k] = ab_x[k] * ig_335[k]
                       + kg_335[k];

            t_468[k] = ab_x[k] * ig_336[k]
                       + kg_336[k];

            t_469[k] = ab_x[k] * ig_337[k]
                       + kg_337[k];

            t_470[k] = ab_x[k] * ig_338[k]
                       + kg_338[k];

            t_471[k] = ab_x[k] * ig_339[k]
                       + kg_339[k];
        }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, t_476, ab_x, ig_340, ig_341, ig_342, \
                         ig_343, ig_344, kg_340, kg_341, kg_342, kg_343, \
                         kg_344 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_472[k] = ab_x[k] * ig_340[k]
                       + kg_340[k];

            t_473[k] = ab_x[k] * ig_341[k]
                       + kg_341[k];

            t_474[k] = ab_x[k] * ig_342[k]
                       + kg_342[k];

            t_475[k] = ab_x[k] * ig_343[k]
                       + kg_343[k];

            t_476[k] = ab_x[k] * ig_344[k]
                       + kg_344[k];
        }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, ab_y, ig_340, ig_341, ig_342, \
                         ig_343, ig_344, kg_445, kg_446, kg_447, kg_448, \
                         kg_449 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_477[k] = ab_y[k] * ig_340[k]
                       + kg_445[k];

            t_478[k] = ab_y[k] * ig_341[k]
                       + kg_446[k];

            t_479[k] = ab_y[k] * ig_342[k]
                       + kg_447[k];

            t_480[k] = ab_y[k] * ig_343[k]
                       + kg_448[k];

            t_481[k] = ab_y[k] * ig_344[k]
                       + kg_449[k];
        }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, ab_x, ab_z, ig_344, ig_345, ig_346, \
                         ig_347, kg_345, kg_346, kg_347, kg_464 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_482[k] = ab_z[k] * ig_344[k]
                       + kg_464[k];

            t_483[k] = ab_x[k] * ig_345[k]
                       + kg_345[k];

            t_484[k] = ab_x[k] * ig_346[k]
                       + kg_346[k];

            t_485[k] = ab_x[k] * ig_347[k]
                       + kg_347[k];
        }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, t_490, ab_x, ig_348, ig_349, ig_350, \
                         ig_351, ig_352, kg_348, kg_349, kg_350, kg_351, \
                         kg_352 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_486[k] = ab_x[k] * ig_348[k]
                       + kg_348[k];

            t_487[k] = ab_x[k] * ig_349[k]
                       + kg_349[k];

            t_488[k] = ab_x[k] * ig_350[k]
                       + kg_350[k];

            t_489[k] = ab_x[k] * ig_351[k]
                       + kg_351[k];

            t_490[k] = ab_x[k] * ig_352[k]
                       + kg_352[k];
        }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, ab_x, ig_353, ig_354, ig_355, \
                         ig_356, ig_357, kg_353, kg_354, kg_355, kg_356, \
                         kg_357 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_491[k] = ab_x[k] * ig_353[k]
                       + kg_353[k];

            t_492[k] = ab_x[k] * ig_354[k]
                       + kg_354[k];

            t_493[k] = ab_x[k] * ig_355[k]
                       + kg_355[k];

            t_494[k] = ab_x[k] * ig_356[k]
                       + kg_356[k];

            t_495[k] = ab_x[k] * ig_357[k]
                       + kg_357[k];
        }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, ab_x, ab_y, ig_355, ig_356, ig_358, \
                         ig_359, kg_358, kg_359, kg_460, kg_461 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_496[k] = ab_x[k] * ig_358[k]
                       + kg_358[k];

            t_497[k] = ab_x[k] * ig_359[k]
                       + kg_359[k];

            t_498[k] = ab_y[k] * ig_355[k]
                       + kg_460[k];

            t_499[k] = ab_y[k] * ig_356[k]
                       + kg_461[k];
        }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, ab_y, ab_z, ig_357, ig_358, ig_359, \
                         kg_462, kg_463, kg_464, kg_479 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_500[k] = ab_y[k] * ig_357[k]
                       + kg_462[k];

            t_501[k] = ab_y[k] * ig_358[k]
                       + kg_463[k];

            t_502[k] = ab_y[k] * ig_359[k]
                       + kg_464[k];

            t_503[k] = ab_z[k] * ig_359[k]
                       + kg_479[k];
        }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, ab_x, ig_360, ig_361, ig_362, \
                         ig_363, ig_364, kg_360, kg_361, kg_362, kg_363, \
                         kg_364 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_504[k] = ab_x[k] * ig_360[k]
                       + kg_360[k];

            t_505[k] = ab_x[k] * ig_361[k]
                       + kg_361[k];

            t_506[k] = ab_x[k] * ig_362[k]
                       + kg_362[k];

            t_507[k] = ab_x[k] * ig_363[k]
                       + kg_363[k];

            t_508[k] = ab_x[k] * ig_364[k]
                       + kg_364[k];
        }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, t_513, ab_x, ig_365, ig_366, ig_367, \
                         ig_368, ig_369, kg_365, kg_366, kg_367, kg_368, \
                         kg_369 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_509[k] = ab_x[k] * ig_365[k]
                       + kg_365[k];

            t_510[k] = ab_x[k] * ig_366[k]
                       + kg_366[k];

            t_511[k] = ab_x[k] * ig_367[k]
                       + kg_367[k];

            t_512[k] = ab_x[k] * ig_368[k]
                       + kg_368[k];

            t_513[k] = ab_x[k] * ig_369[k]
                       + kg_369[k];
        }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, t_518, ab_x, ig_370, ig_371, ig_372, \
                         ig_373, ig_374, kg_370, kg_371, kg_372, kg_373, \
                         kg_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_514[k] = ab_x[k] * ig_370[k]
                       + kg_370[k];

            t_515[k] = ab_x[k] * ig_371[k]
                       + kg_371[k];

            t_516[k] = ab_x[k] * ig_372[k]
                       + kg_372[k];

            t_517[k] = ab_x[k] * ig_373[k]
                       + kg_373[k];

            t_518[k] = ab_x[k] * ig_374[k]
                       + kg_374[k];
        }

#pragma omp simd aligned(t_519, t_520, t_521, t_522, t_523, ab_y, ig_370, ig_371, ig_372, \
                         ig_373, ig_374, kg_475, kg_476, kg_477, kg_478, \
                         kg_479 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_519[k] = ab_y[k] * ig_370[k]
                       + kg_475[k];

            t_520[k] = ab_y[k] * ig_371[k]
                       + kg_476[k];

            t_521[k] = ab_y[k] * ig_372[k]
                       + kg_477[k];

            t_522[k] = ab_y[k] * ig_373[k]
                       + kg_478[k];

            t_523[k] = ab_y[k] * ig_374[k]
                       + kg_479[k];
        }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, ab_x, ab_z, ig_374, ig_375, ig_376, \
                         ig_377, kg_375, kg_376, kg_377, kg_494 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_524[k] = ab_z[k] * ig_374[k]
                       + kg_494[k];

            t_525[k] = ab_x[k] * ig_375[k]
                       + kg_375[k];

            t_526[k] = ab_x[k] * ig_376[k]
                       + kg_376[k];

            t_527[k] = ab_x[k] * ig_377[k]
                       + kg_377[k];
        }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, t_532, ab_x, ig_378, ig_379, ig_380, \
                         ig_381, ig_382, kg_378, kg_379, kg_380, kg_381, \
                         kg_382 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_528[k] = ab_x[k] * ig_378[k]
                       + kg_378[k];

            t_529[k] = ab_x[k] * ig_379[k]
                       + kg_379[k];

            t_530[k] = ab_x[k] * ig_380[k]
                       + kg_380[k];

            t_531[k] = ab_x[k] * ig_381[k]
                       + kg_381[k];

            t_532[k] = ab_x[k] * ig_382[k]
                       + kg_382[k];
        }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, t_537, ab_x, ig_383, ig_384, ig_385, \
                         ig_386, ig_387, kg_383, kg_384, kg_385, kg_386, \
                         kg_387 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_533[k] = ab_x[k] * ig_383[k]
                       + kg_383[k];

            t_534[k] = ab_x[k] * ig_384[k]
                       + kg_384[k];

            t_535[k] = ab_x[k] * ig_385[k]
                       + kg_385[k];

            t_536[k] = ab_x[k] * ig_386[k]
                       + kg_386[k];

            t_537[k] = ab_x[k] * ig_387[k]
                       + kg_387[k];
        }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, ab_x, ab_y, ig_385, ig_386, ig_388, \
                         ig_389, kg_388, kg_389, kg_490, kg_491 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_538[k] = ab_x[k] * ig_388[k]
                       + kg_388[k];

            t_539[k] = ab_x[k] * ig_389[k]
                       + kg_389[k];

            t_540[k] = ab_y[k] * ig_385[k]
                       + kg_490[k];

            t_541[k] = ab_y[k] * ig_386[k]
                       + kg_491[k];
        }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, ab_y, ab_z, ig_387, ig_388, ig_389, \
                         kg_492, kg_493, kg_494, kg_509 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_542[k] = ab_y[k] * ig_387[k]
                       + kg_492[k];

            t_543[k] = ab_y[k] * ig_388[k]
                       + kg_493[k];

            t_544[k] = ab_y[k] * ig_389[k]
                       + kg_494[k];

            t_545[k] = ab_z[k] * ig_389[k]
                       + kg_509[k];
        }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, t_550, ab_x, ig_390, ig_391, ig_392, \
                         ig_393, ig_394, kg_390, kg_391, kg_392, kg_393, \
                         kg_394 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_546[k] = ab_x[k] * ig_390[k]
                       + kg_390[k];

            t_547[k] = ab_x[k] * ig_391[k]
                       + kg_391[k];

            t_548[k] = ab_x[k] * ig_392[k]
                       + kg_392[k];

            t_549[k] = ab_x[k] * ig_393[k]
                       + kg_393[k];

            t_550[k] = ab_x[k] * ig_394[k]
                       + kg_394[k];
        }

#pragma omp simd aligned(t_551, t_552, t_553, t_554, t_555, ab_x, ig_395, ig_396, ig_397, \
                         ig_398, ig_399, kg_395, kg_396, kg_397, kg_398, \
                         kg_399 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_551[k] = ab_x[k] * ig_395[k]
                       + kg_395[k];

            t_552[k] = ab_x[k] * ig_396[k]
                       + kg_396[k];

            t_553[k] = ab_x[k] * ig_397[k]
                       + kg_397[k];

            t_554[k] = ab_x[k] * ig_398[k]
                       + kg_398[k];

            t_555[k] = ab_x[k] * ig_399[k]
                       + kg_399[k];
        }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, ab_x, ig_400, ig_401, ig_402, \
                         ig_403, ig_404, kg_400, kg_401, kg_402, kg_403, \
                         kg_404 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_556[k] = ab_x[k] * ig_400[k]
                       + kg_400[k];

            t_557[k] = ab_x[k] * ig_401[k]
                       + kg_401[k];

            t_558[k] = ab_x[k] * ig_402[k]
                       + kg_402[k];

            t_559[k] = ab_x[k] * ig_403[k]
                       + kg_403[k];

            t_560[k] = ab_x[k] * ig_404[k]
                       + kg_404[k];
        }

#pragma omp simd aligned(t_561, t_562, t_563, t_564, t_565, ab_y, ig_400, ig_401, ig_402, \
                         ig_403, ig_404, kg_505, kg_506, kg_507, kg_508, \
                         kg_509 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_561[k] = ab_y[k] * ig_400[k]
                       + kg_505[k];

            t_562[k] = ab_y[k] * ig_401[k]
                       + kg_506[k];

            t_563[k] = ab_y[k] * ig_402[k]
                       + kg_507[k];

            t_564[k] = ab_y[k] * ig_403[k]
                       + kg_508[k];

            t_565[k] = ab_y[k] * ig_404[k]
                       + kg_509[k];
        }
    }
}

static auto
compute_hrr_ih_out_of_first_piece4(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t ig, const size_t kg,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ig_404 = buffer.data(ig + 404 * ncomps + c);
        const auto *ig_405 = buffer.data(ig + 405 * ncomps + c);
        const auto *ig_406 = buffer.data(ig + 406 * ncomps + c);
        const auto *ig_407 = buffer.data(ig + 407 * ncomps + c);
        const auto *ig_408 = buffer.data(ig + 408 * ncomps + c);
        const auto *ig_409 = buffer.data(ig + 409 * ncomps + c);
        const auto *ig_410 = buffer.data(ig + 410 * ncomps + c);
        const auto *ig_411 = buffer.data(ig + 411 * ncomps + c);
        const auto *ig_412 = buffer.data(ig + 412 * ncomps + c);
        const auto *ig_413 = buffer.data(ig + 413 * ncomps + c);
        const auto *ig_414 = buffer.data(ig + 414 * ncomps + c);
        const auto *ig_415 = buffer.data(ig + 415 * ncomps + c);
        const auto *ig_416 = buffer.data(ig + 416 * ncomps + c);
        const auto *ig_417 = buffer.data(ig + 417 * ncomps + c);
        const auto *ig_418 = buffer.data(ig + 418 * ncomps + c);
        const auto *ig_419 = buffer.data(ig + 419 * ncomps + c);

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
        const auto *kg_520 = buffer.data(kg + 520 * ncomps + c);
        const auto *kg_521 = buffer.data(kg + 521 * ncomps + c);
        const auto *kg_522 = buffer.data(kg + 522 * ncomps + c);
        const auto *kg_523 = buffer.data(kg + 523 * ncomps + c);
        const auto *kg_524 = buffer.data(kg + 524 * ncomps + c);
        const auto *kg_539 = buffer.data(kg + 539 * ncomps + c);

#pragma omp simd aligned(t_566, t_567, t_568, t_569, ab_x, ab_z, ig_404, ig_405, ig_406, \
                         ig_407, kg_405, kg_406, kg_407, kg_524 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_566[k] = ab_z[k] * ig_404[k]
                       + kg_524[k];

            t_567[k] = ab_x[k] * ig_405[k]
                       + kg_405[k];

            t_568[k] = ab_x[k] * ig_406[k]
                       + kg_406[k];

            t_569[k] = ab_x[k] * ig_407[k]
                       + kg_407[k];
        }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, ab_x, ig_408, ig_409, ig_410, \
                         ig_411, ig_412, kg_408, kg_409, kg_410, kg_411, \
                         kg_412 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_570[k] = ab_x[k] * ig_408[k]
                       + kg_408[k];

            t_571[k] = ab_x[k] * ig_409[k]
                       + kg_409[k];

            t_572[k] = ab_x[k] * ig_410[k]
                       + kg_410[k];

            t_573[k] = ab_x[k] * ig_411[k]
                       + kg_411[k];

            t_574[k] = ab_x[k] * ig_412[k]
                       + kg_412[k];
        }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, ab_x, ig_413, ig_414, ig_415, \
                         ig_416, ig_417, kg_413, kg_414, kg_415, kg_416, \
                         kg_417 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_575[k] = ab_x[k] * ig_413[k]
                       + kg_413[k];

            t_576[k] = ab_x[k] * ig_414[k]
                       + kg_414[k];

            t_577[k] = ab_x[k] * ig_415[k]
                       + kg_415[k];

            t_578[k] = ab_x[k] * ig_416[k]
                       + kg_416[k];

            t_579[k] = ab_x[k] * ig_417[k]
                       + kg_417[k];
        }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, ab_x, ab_y, ig_415, ig_416, ig_418, \
                         ig_419, kg_418, kg_419, kg_520, kg_521 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_580[k] = ab_x[k] * ig_418[k]
                       + kg_418[k];

            t_581[k] = ab_x[k] * ig_419[k]
                       + kg_419[k];

            t_582[k] = ab_y[k] * ig_415[k]
                       + kg_520[k];

            t_583[k] = ab_y[k] * ig_416[k]
                       + kg_521[k];
        }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, ab_y, ab_z, ig_417, ig_418, ig_419, \
                         kg_522, kg_523, kg_524, kg_539 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_584[k] = ab_y[k] * ig_417[k]
                       + kg_522[k];

            t_585[k] = ab_y[k] * ig_418[k]
                       + kg_523[k];

            t_586[k] = ab_y[k] * ig_419[k]
                       + kg_524[k];

            t_587[k] = ab_z[k] * ig_419[k]
                       + kg_539[k];
        }
    }
}

auto
compute_hrr_ih_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t ig, const size_t kg,
                            const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_ih_out_of_first_piece0(buffer, coordinates, target, ig, kg, ncomps, nmax);

    compute_hrr_ih_out_of_first_piece1(buffer, coordinates, target, ig, kg, ncomps, nmax);

    compute_hrr_ih_out_of_first_piece2(buffer, coordinates, target, ig, kg, ncomps, nmax);

    compute_hrr_ih_out_of_first_piece3(buffer, coordinates, target, ig, kg, ncomps, nmax);

    compute_hrr_ih_out_of_first_piece4(buffer, coordinates, target, ig, kg, ncomps, nmax);
}

static auto
compute_hrr_ih_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t ig, const size_t kg, const size_t ncomps,
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

        const auto *ig_0 = buffer.data(ig + 0 * ncomps + c);
        const auto *ig_1 = buffer.data(ig + 1 * ncomps + c);
        const auto *ig_2 = buffer.data(ig + 2 * ncomps + c);
        const auto *ig_3 = buffer.data(ig + 3 * ncomps + c);
        const auto *ig_4 = buffer.data(ig + 4 * ncomps + c);
        const auto *ig_5 = buffer.data(ig + 5 * ncomps + c);
        const auto *ig_6 = buffer.data(ig + 6 * ncomps + c);
        const auto *ig_7 = buffer.data(ig + 7 * ncomps + c);
        const auto *ig_8 = buffer.data(ig + 8 * ncomps + c);
        const auto *ig_9 = buffer.data(ig + 9 * ncomps + c);
        const auto *ig_10 = buffer.data(ig + 10 * ncomps + c);
        const auto *ig_11 = buffer.data(ig + 11 * ncomps + c);
        const auto *ig_12 = buffer.data(ig + 12 * ncomps + c);
        const auto *ig_13 = buffer.data(ig + 13 * ncomps + c);
        const auto *ig_14 = buffer.data(ig + 14 * ncomps + c);
        const auto *ig_15 = buffer.data(ig + 15 * ncomps + c);
        const auto *ig_16 = buffer.data(ig + 16 * ncomps + c);
        const auto *ig_17 = buffer.data(ig + 17 * ncomps + c);
        const auto *ig_18 = buffer.data(ig + 18 * ncomps + c);
        const auto *ig_19 = buffer.data(ig + 19 * ncomps + c);
        const auto *ig_20 = buffer.data(ig + 20 * ncomps + c);
        const auto *ig_21 = buffer.data(ig + 21 * ncomps + c);
        const auto *ig_22 = buffer.data(ig + 22 * ncomps + c);
        const auto *ig_23 = buffer.data(ig + 23 * ncomps + c);
        const auto *ig_24 = buffer.data(ig + 24 * ncomps + c);
        const auto *ig_25 = buffer.data(ig + 25 * ncomps + c);
        const auto *ig_26 = buffer.data(ig + 26 * ncomps + c);
        const auto *ig_27 = buffer.data(ig + 27 * ncomps + c);
        const auto *ig_28 = buffer.data(ig + 28 * ncomps + c);
        const auto *ig_29 = buffer.data(ig + 29 * ncomps + c);
        const auto *ig_30 = buffer.data(ig + 30 * ncomps + c);
        const auto *ig_31 = buffer.data(ig + 31 * ncomps + c);
        const auto *ig_32 = buffer.data(ig + 32 * ncomps + c);
        const auto *ig_33 = buffer.data(ig + 33 * ncomps + c);
        const auto *ig_34 = buffer.data(ig + 34 * ncomps + c);
        const auto *ig_35 = buffer.data(ig + 35 * ncomps + c);
        const auto *ig_36 = buffer.data(ig + 36 * ncomps + c);
        const auto *ig_37 = buffer.data(ig + 37 * ncomps + c);
        const auto *ig_38 = buffer.data(ig + 38 * ncomps + c);
        const auto *ig_39 = buffer.data(ig + 39 * ncomps + c);
        const auto *ig_40 = buffer.data(ig + 40 * ncomps + c);
        const auto *ig_41 = buffer.data(ig + 41 * ncomps + c);
        const auto *ig_42 = buffer.data(ig + 42 * ncomps + c);
        const auto *ig_43 = buffer.data(ig + 43 * ncomps + c);
        const auto *ig_44 = buffer.data(ig + 44 * ncomps + c);
        const auto *ig_45 = buffer.data(ig + 45 * ncomps + c);
        const auto *ig_46 = buffer.data(ig + 46 * ncomps + c);
        const auto *ig_47 = buffer.data(ig + 47 * ncomps + c);
        const auto *ig_48 = buffer.data(ig + 48 * ncomps + c);
        const auto *ig_49 = buffer.data(ig + 49 * ncomps + c);
        const auto *ig_50 = buffer.data(ig + 50 * ncomps + c);
        const auto *ig_51 = buffer.data(ig + 51 * ncomps + c);
        const auto *ig_52 = buffer.data(ig + 52 * ncomps + c);
        const auto *ig_53 = buffer.data(ig + 53 * ncomps + c);
        const auto *ig_54 = buffer.data(ig + 54 * ncomps + c);
        const auto *ig_55 = buffer.data(ig + 55 * ncomps + c);
        const auto *ig_56 = buffer.data(ig + 56 * ncomps + c);
        const auto *ig_57 = buffer.data(ig + 57 * ncomps + c);
        const auto *ig_58 = buffer.data(ig + 58 * ncomps + c);
        const auto *ig_59 = buffer.data(ig + 59 * ncomps + c);
        const auto *ig_60 = buffer.data(ig + 60 * ncomps + c);
        const auto *ig_61 = buffer.data(ig + 61 * ncomps + c);
        const auto *ig_62 = buffer.data(ig + 62 * ncomps + c);
        const auto *ig_63 = buffer.data(ig + 63 * ncomps + c);
        const auto *ig_64 = buffer.data(ig + 64 * ncomps + c);
        const auto *ig_65 = buffer.data(ig + 65 * ncomps + c);
        const auto *ig_66 = buffer.data(ig + 66 * ncomps + c);
        const auto *ig_67 = buffer.data(ig + 67 * ncomps + c);
        const auto *ig_68 = buffer.data(ig + 68 * ncomps + c);
        const auto *ig_69 = buffer.data(ig + 69 * ncomps + c);
        const auto *ig_70 = buffer.data(ig + 70 * ncomps + c);
        const auto *ig_71 = buffer.data(ig + 71 * ncomps + c);
        const auto *ig_72 = buffer.data(ig + 72 * ncomps + c);
        const auto *ig_73 = buffer.data(ig + 73 * ncomps + c);
        const auto *ig_74 = buffer.data(ig + 74 * ncomps + c);
        const auto *ig_75 = buffer.data(ig + 75 * ncomps + c);
        const auto *ig_76 = buffer.data(ig + 76 * ncomps + c);
        const auto *ig_77 = buffer.data(ig + 77 * ncomps + c);
        const auto *ig_78 = buffer.data(ig + 78 * ncomps + c);
        const auto *ig_79 = buffer.data(ig + 79 * ncomps + c);
        const auto *ig_80 = buffer.data(ig + 80 * ncomps + c);
        const auto *ig_81 = buffer.data(ig + 81 * ncomps + c);
        const auto *ig_82 = buffer.data(ig + 82 * ncomps + c);
        const auto *ig_83 = buffer.data(ig + 83 * ncomps + c);
        const auto *ig_84 = buffer.data(ig + 84 * ncomps + c);
        const auto *ig_85 = buffer.data(ig + 85 * ncomps + c);
        const auto *ig_86 = buffer.data(ig + 86 * ncomps + c);
        const auto *ig_87 = buffer.data(ig + 87 * ncomps + c);
        const auto *ig_88 = buffer.data(ig + 88 * ncomps + c);
        const auto *ig_89 = buffer.data(ig + 89 * ncomps + c);
        const auto *ig_90 = buffer.data(ig + 90 * ncomps + c);
        const auto *ig_91 = buffer.data(ig + 91 * ncomps + c);
        const auto *ig_92 = buffer.data(ig + 92 * ncomps + c);
        const auto *ig_93 = buffer.data(ig + 93 * ncomps + c);
        const auto *ig_94 = buffer.data(ig + 94 * ncomps + c);
        const auto *ig_95 = buffer.data(ig + 95 * ncomps + c);
        const auto *ig_96 = buffer.data(ig + 96 * ncomps + c);
        const auto *ig_97 = buffer.data(ig + 97 * ncomps + c);
        const auto *ig_98 = buffer.data(ig + 98 * ncomps + c);
        const auto *ig_99 = buffer.data(ig + 99 * ncomps + c);
        const auto *ig_100 = buffer.data(ig + 100 * ncomps + c);
        const auto *ig_101 = buffer.data(ig + 101 * ncomps + c);
        const auto *ig_102 = buffer.data(ig + 102 * ncomps + c);
        const auto *ig_103 = buffer.data(ig + 103 * ncomps + c);
        const auto *ig_104 = buffer.data(ig + 104 * ncomps + c);

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
        const auto *kg_115 = buffer.data(kg + 115 * ncomps + c);
        const auto *kg_116 = buffer.data(kg + 116 * ncomps + c);
        const auto *kg_117 = buffer.data(kg + 117 * ncomps + c);
        const auto *kg_118 = buffer.data(kg + 118 * ncomps + c);
        const auto *kg_119 = buffer.data(kg + 119 * ncomps + c);
        const auto *kg_130 = buffer.data(kg + 130 * ncomps + c);
        const auto *kg_131 = buffer.data(kg + 131 * ncomps + c);
        const auto *kg_132 = buffer.data(kg + 132 * ncomps + c);
        const auto *kg_133 = buffer.data(kg + 133 * ncomps + c);
        const auto *kg_134 = buffer.data(kg + 134 * ncomps + c);
        const auto *kg_149 = buffer.data(kg + 149 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ig_0, ig_1, ig_2, ig_3, ig_4, kg_0, \
                         kg_1, kg_2, kg_3, kg_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * ig_0[k]
                     + kg_0[k];

            t_1[k] = ab_x[k] * ig_1[k]
                     + kg_1[k];

            t_2[k] = ab_x[k] * ig_2[k]
                     + kg_2[k];

            t_3[k] = ab_x[k] * ig_3[k]
                     + kg_3[k];

            t_4[k] = ab_x[k] * ig_4[k]
                     + kg_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ig_5, ig_6, ig_7, ig_8, ig_9, kg_5, \
                         kg_6, kg_7, kg_8, kg_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * ig_5[k]
                     + kg_5[k];

            t_6[k] = ab_x[k] * ig_6[k]
                     + kg_6[k];

            t_7[k] = ab_x[k] * ig_7[k]
                     + kg_7[k];

            t_8[k] = ab_x[k] * ig_8[k]
                     + kg_8[k];

            t_9[k] = ab_x[k] * ig_9[k]
                     + kg_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, ig_10, ig_11, ig_12, ig_13, \
                         ig_14, kg_10, kg_11, kg_12, kg_13, kg_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_x[k] * ig_10[k]
                      + kg_10[k];

            t_11[k] = ab_x[k] * ig_11[k]
                      + kg_11[k];

            t_12[k] = ab_x[k] * ig_12[k]
                      + kg_12[k];

            t_13[k] = ab_x[k] * ig_13[k]
                      + kg_13[k];

            t_14[k] = ab_x[k] * ig_14[k]
                      + kg_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_y, ig_10, ig_11, ig_12, ig_13, \
                         ig_14, kg_25, kg_26, kg_27, kg_28, kg_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_y[k] * ig_10[k]
                      + kg_25[k];

            t_16[k] = ab_y[k] * ig_11[k]
                      + kg_26[k];

            t_17[k] = ab_y[k] * ig_12[k]
                      + kg_27[k];

            t_18[k] = ab_y[k] * ig_13[k]
                      + kg_28[k];

            t_19[k] = ab_y[k] * ig_14[k]
                      + kg_29[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, ab_x, ab_z, ig_14, ig_15, ig_16, ig_17, \
                         kg_15, kg_16, kg_17, kg_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_z[k] * ig_14[k]
                      + kg_44[k];

            t_21[k] = ab_x[k] * ig_15[k]
                      + kg_15[k];

            t_22[k] = ab_x[k] * ig_16[k]
                      + kg_16[k];

            t_23[k] = ab_x[k] * ig_17[k]
                      + kg_17[k];
        }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, ab_x, ig_18, ig_19, ig_20, ig_21, \
                         ig_22, kg_18, kg_19, kg_20, kg_21, kg_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_24[k] = ab_x[k] * ig_18[k]
                      + kg_18[k];

            t_25[k] = ab_x[k] * ig_19[k]
                      + kg_19[k];

            t_26[k] = ab_x[k] * ig_20[k]
                      + kg_20[k];

            t_27[k] = ab_x[k] * ig_21[k]
                      + kg_21[k];

            t_28[k] = ab_x[k] * ig_22[k]
                      + kg_22[k];
        }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, ab_x, ig_23, ig_24, ig_25, ig_26, \
                         ig_27, kg_23, kg_24, kg_25, kg_26, kg_27 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_29[k] = ab_x[k] * ig_23[k]
                      + kg_23[k];

            t_30[k] = ab_x[k] * ig_24[k]
                      + kg_24[k];

            t_31[k] = ab_x[k] * ig_25[k]
                      + kg_25[k];

            t_32[k] = ab_x[k] * ig_26[k]
                      + kg_26[k];

            t_33[k] = ab_x[k] * ig_27[k]
                      + kg_27[k];
        }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, ab_x, ab_y, ig_25, ig_26, ig_28, ig_29, \
                         kg_28, kg_29, kg_55, kg_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_34[k] = ab_x[k] * ig_28[k]
                      + kg_28[k];

            t_35[k] = ab_x[k] * ig_29[k]
                      + kg_29[k];

            t_36[k] = ab_y[k] * ig_25[k]
                      + kg_55[k];

            t_37[k] = ab_y[k] * ig_26[k]
                      + kg_56[k];
        }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, ab_y, ab_z, ig_27, ig_28, ig_29, kg_57, \
                         kg_58, kg_59, kg_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_38[k] = ab_y[k] * ig_27[k]
                      + kg_57[k];

            t_39[k] = ab_y[k] * ig_28[k]
                      + kg_58[k];

            t_40[k] = ab_y[k] * ig_29[k]
                      + kg_59[k];

            t_41[k] = ab_z[k] * ig_29[k]
                      + kg_74[k];
        }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, ab_x, ig_30, ig_31, ig_32, ig_33, \
                         ig_34, kg_30, kg_31, kg_32, kg_33, kg_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_42[k] = ab_x[k] * ig_30[k]
                      + kg_30[k];

            t_43[k] = ab_x[k] * ig_31[k]
                      + kg_31[k];

            t_44[k] = ab_x[k] * ig_32[k]
                      + kg_32[k];

            t_45[k] = ab_x[k] * ig_33[k]
                      + kg_33[k];

            t_46[k] = ab_x[k] * ig_34[k]
                      + kg_34[k];
        }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, ab_x, ig_35, ig_36, ig_37, ig_38, \
                         ig_39, kg_35, kg_36, kg_37, kg_38, kg_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_47[k] = ab_x[k] * ig_35[k]
                      + kg_35[k];

            t_48[k] = ab_x[k] * ig_36[k]
                      + kg_36[k];

            t_49[k] = ab_x[k] * ig_37[k]
                      + kg_37[k];

            t_50[k] = ab_x[k] * ig_38[k]
                      + kg_38[k];

            t_51[k] = ab_x[k] * ig_39[k]
                      + kg_39[k];
        }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, ab_x, ig_40, ig_41, ig_42, ig_43, \
                         ig_44, kg_40, kg_41, kg_42, kg_43, kg_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_52[k] = ab_x[k] * ig_40[k]
                      + kg_40[k];

            t_53[k] = ab_x[k] * ig_41[k]
                      + kg_41[k];

            t_54[k] = ab_x[k] * ig_42[k]
                      + kg_42[k];

            t_55[k] = ab_x[k] * ig_43[k]
                      + kg_43[k];

            t_56[k] = ab_x[k] * ig_44[k]
                      + kg_44[k];
        }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, ab_y, ig_40, ig_41, ig_42, ig_43, \
                         ig_44, kg_70, kg_71, kg_72, kg_73, kg_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_57[k] = ab_y[k] * ig_40[k]
                      + kg_70[k];

            t_58[k] = ab_y[k] * ig_41[k]
                      + kg_71[k];

            t_59[k] = ab_y[k] * ig_42[k]
                      + kg_72[k];

            t_60[k] = ab_y[k] * ig_43[k]
                      + kg_73[k];

            t_61[k] = ab_y[k] * ig_44[k]
                      + kg_74[k];
        }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, ab_x, ab_z, ig_44, ig_45, ig_46, ig_47, \
                         kg_45, kg_46, kg_47, kg_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_62[k] = ab_z[k] * ig_44[k]
                      + kg_89[k];

            t_63[k] = ab_x[k] * ig_45[k]
                      + kg_45[k];

            t_64[k] = ab_x[k] * ig_46[k]
                      + kg_46[k];

            t_65[k] = ab_x[k] * ig_47[k]
                      + kg_47[k];
        }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, ab_x, ig_48, ig_49, ig_50, ig_51, \
                         ig_52, kg_48, kg_49, kg_50, kg_51, kg_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_66[k] = ab_x[k] * ig_48[k]
                      + kg_48[k];

            t_67[k] = ab_x[k] * ig_49[k]
                      + kg_49[k];

            t_68[k] = ab_x[k] * ig_50[k]
                      + kg_50[k];

            t_69[k] = ab_x[k] * ig_51[k]
                      + kg_51[k];

            t_70[k] = ab_x[k] * ig_52[k]
                      + kg_52[k];
        }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, ab_x, ig_53, ig_54, ig_55, ig_56, \
                         ig_57, kg_53, kg_54, kg_55, kg_56, kg_57 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_71[k] = ab_x[k] * ig_53[k]
                      + kg_53[k];

            t_72[k] = ab_x[k] * ig_54[k]
                      + kg_54[k];

            t_73[k] = ab_x[k] * ig_55[k]
                      + kg_55[k];

            t_74[k] = ab_x[k] * ig_56[k]
                      + kg_56[k];

            t_75[k] = ab_x[k] * ig_57[k]
                      + kg_57[k];
        }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, ab_x, ab_y, ig_55, ig_56, ig_58, ig_59, \
                         kg_58, kg_59, kg_100, kg_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_76[k] = ab_x[k] * ig_58[k]
                      + kg_58[k];

            t_77[k] = ab_x[k] * ig_59[k]
                      + kg_59[k];

            t_78[k] = ab_y[k] * ig_55[k]
                      + kg_100[k];

            t_79[k] = ab_y[k] * ig_56[k]
                      + kg_101[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, ab_y, ab_z, ig_57, ig_58, ig_59, kg_102, \
                         kg_103, kg_104, kg_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_y[k] * ig_57[k]
                      + kg_102[k];

            t_81[k] = ab_y[k] * ig_58[k]
                      + kg_103[k];

            t_82[k] = ab_y[k] * ig_59[k]
                      + kg_104[k];

            t_83[k] = ab_z[k] * ig_59[k]
                      + kg_119[k];
        }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, ab_x, ig_60, ig_61, ig_62, ig_63, \
                         ig_64, kg_60, kg_61, kg_62, kg_63, kg_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_84[k] = ab_x[k] * ig_60[k]
                      + kg_60[k];

            t_85[k] = ab_x[k] * ig_61[k]
                      + kg_61[k];

            t_86[k] = ab_x[k] * ig_62[k]
                      + kg_62[k];

            t_87[k] = ab_x[k] * ig_63[k]
                      + kg_63[k];

            t_88[k] = ab_x[k] * ig_64[k]
                      + kg_64[k];
        }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, ab_x, ig_65, ig_66, ig_67, ig_68, \
                         ig_69, kg_65, kg_66, kg_67, kg_68, kg_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_89[k] = ab_x[k] * ig_65[k]
                      + kg_65[k];

            t_90[k] = ab_x[k] * ig_66[k]
                      + kg_66[k];

            t_91[k] = ab_x[k] * ig_67[k]
                      + kg_67[k];

            t_92[k] = ab_x[k] * ig_68[k]
                      + kg_68[k];

            t_93[k] = ab_x[k] * ig_69[k]
                      + kg_69[k];
        }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, ab_x, ig_70, ig_71, ig_72, ig_73, \
                         ig_74, kg_70, kg_71, kg_72, kg_73, kg_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_94[k] = ab_x[k] * ig_70[k]
                      + kg_70[k];

            t_95[k] = ab_x[k] * ig_71[k]
                      + kg_71[k];

            t_96[k] = ab_x[k] * ig_72[k]
                      + kg_72[k];

            t_97[k] = ab_x[k] * ig_73[k]
                      + kg_73[k];

            t_98[k] = ab_x[k] * ig_74[k]
                      + kg_74[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, ab_y, ig_70, ig_71, ig_72, ig_73, \
                         ig_74, kg_115, kg_116, kg_117, kg_118, \
                         kg_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_y[k] * ig_70[k]
                      + kg_115[k];

            t_100[k] = ab_y[k] * ig_71[k]
                       + kg_116[k];

            t_101[k] = ab_y[k] * ig_72[k]
                       + kg_117[k];

            t_102[k] = ab_y[k] * ig_73[k]
                       + kg_118[k];

            t_103[k] = ab_y[k] * ig_74[k]
                       + kg_119[k];
        }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, ab_x, ab_z, ig_74, ig_75, ig_76, ig_77, \
                         kg_75, kg_76, kg_77, kg_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_104[k] = ab_z[k] * ig_74[k]
                       + kg_134[k];

            t_105[k] = ab_x[k] * ig_75[k]
                       + kg_75[k];

            t_106[k] = ab_x[k] * ig_76[k]
                       + kg_76[k];

            t_107[k] = ab_x[k] * ig_77[k]
                       + kg_77[k];
        }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, ab_x, ig_78, ig_79, ig_80, ig_81, \
                         ig_82, kg_78, kg_79, kg_80, kg_81, kg_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_108[k] = ab_x[k] * ig_78[k]
                       + kg_78[k];

            t_109[k] = ab_x[k] * ig_79[k]
                       + kg_79[k];

            t_110[k] = ab_x[k] * ig_80[k]
                       + kg_80[k];

            t_111[k] = ab_x[k] * ig_81[k]
                       + kg_81[k];

            t_112[k] = ab_x[k] * ig_82[k]
                       + kg_82[k];
        }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, ab_x, ig_83, ig_84, ig_85, ig_86, \
                         ig_87, kg_83, kg_84, kg_85, kg_86, kg_87 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_113[k] = ab_x[k] * ig_83[k]
                       + kg_83[k];

            t_114[k] = ab_x[k] * ig_84[k]
                       + kg_84[k];

            t_115[k] = ab_x[k] * ig_85[k]
                       + kg_85[k];

            t_116[k] = ab_x[k] * ig_86[k]
                       + kg_86[k];

            t_117[k] = ab_x[k] * ig_87[k]
                       + kg_87[k];
        }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, ab_x, ab_y, ig_85, ig_86, ig_88, ig_89, \
                         kg_88, kg_89, kg_130, kg_131 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_118[k] = ab_x[k] * ig_88[k]
                       + kg_88[k];

            t_119[k] = ab_x[k] * ig_89[k]
                       + kg_89[k];

            t_120[k] = ab_y[k] * ig_85[k]
                       + kg_130[k];

            t_121[k] = ab_y[k] * ig_86[k]
                       + kg_131[k];
        }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, ab_y, ab_z, ig_87, ig_88, ig_89, kg_132, \
                         kg_133, kg_134, kg_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_122[k] = ab_y[k] * ig_87[k]
                       + kg_132[k];

            t_123[k] = ab_y[k] * ig_88[k]
                       + kg_133[k];

            t_124[k] = ab_y[k] * ig_89[k]
                       + kg_134[k];

            t_125[k] = ab_z[k] * ig_89[k]
                       + kg_149[k];
        }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, ab_x, ig_90, ig_91, ig_92, ig_93, \
                         ig_94, kg_90, kg_91, kg_92, kg_93, kg_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_126[k] = ab_x[k] * ig_90[k]
                       + kg_90[k];

            t_127[k] = ab_x[k] * ig_91[k]
                       + kg_91[k];

            t_128[k] = ab_x[k] * ig_92[k]
                       + kg_92[k];

            t_129[k] = ab_x[k] * ig_93[k]
                       + kg_93[k];

            t_130[k] = ab_x[k] * ig_94[k]
                       + kg_94[k];
        }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, ab_x, ig_95, ig_96, ig_97, ig_98, \
                         ig_99, kg_95, kg_96, kg_97, kg_98, kg_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_131[k] = ab_x[k] * ig_95[k]
                       + kg_95[k];

            t_132[k] = ab_x[k] * ig_96[k]
                       + kg_96[k];

            t_133[k] = ab_x[k] * ig_97[k]
                       + kg_97[k];

            t_134[k] = ab_x[k] * ig_98[k]
                       + kg_98[k];

            t_135[k] = ab_x[k] * ig_99[k]
                       + kg_99[k];
        }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, ab_x, ig_100, ig_101, ig_102, \
                         ig_103, ig_104, kg_100, kg_101, kg_102, kg_103, \
                         kg_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_136[k] = ab_x[k] * ig_100[k]
                       + kg_100[k];

            t_137[k] = ab_x[k] * ig_101[k]
                       + kg_101[k];

            t_138[k] = ab_x[k] * ig_102[k]
                       + kg_102[k];

            t_139[k] = ab_x[k] * ig_103[k]
                       + kg_103[k];

            t_140[k] = ab_x[k] * ig_104[k]
                       + kg_104[k];
        }
    }
}

static auto
compute_hrr_ih_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t ig, const size_t kg, const size_t ncomps,
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

        const auto *ig_100 = buffer.data(ig + 100 * ncomps + c);
        const auto *ig_101 = buffer.data(ig + 101 * ncomps + c);
        const auto *ig_102 = buffer.data(ig + 102 * ncomps + c);
        const auto *ig_103 = buffer.data(ig + 103 * ncomps + c);
        const auto *ig_104 = buffer.data(ig + 104 * ncomps + c);
        const auto *ig_105 = buffer.data(ig + 105 * ncomps + c);
        const auto *ig_106 = buffer.data(ig + 106 * ncomps + c);
        const auto *ig_107 = buffer.data(ig + 107 * ncomps + c);
        const auto *ig_108 = buffer.data(ig + 108 * ncomps + c);
        const auto *ig_109 = buffer.data(ig + 109 * ncomps + c);
        const auto *ig_110 = buffer.data(ig + 110 * ncomps + c);
        const auto *ig_111 = buffer.data(ig + 111 * ncomps + c);
        const auto *ig_112 = buffer.data(ig + 112 * ncomps + c);
        const auto *ig_113 = buffer.data(ig + 113 * ncomps + c);
        const auto *ig_114 = buffer.data(ig + 114 * ncomps + c);
        const auto *ig_115 = buffer.data(ig + 115 * ncomps + c);
        const auto *ig_116 = buffer.data(ig + 116 * ncomps + c);
        const auto *ig_117 = buffer.data(ig + 117 * ncomps + c);
        const auto *ig_118 = buffer.data(ig + 118 * ncomps + c);
        const auto *ig_119 = buffer.data(ig + 119 * ncomps + c);
        const auto *ig_120 = buffer.data(ig + 120 * ncomps + c);
        const auto *ig_121 = buffer.data(ig + 121 * ncomps + c);
        const auto *ig_122 = buffer.data(ig + 122 * ncomps + c);
        const auto *ig_123 = buffer.data(ig + 123 * ncomps + c);
        const auto *ig_124 = buffer.data(ig + 124 * ncomps + c);
        const auto *ig_125 = buffer.data(ig + 125 * ncomps + c);
        const auto *ig_126 = buffer.data(ig + 126 * ncomps + c);
        const auto *ig_127 = buffer.data(ig + 127 * ncomps + c);
        const auto *ig_128 = buffer.data(ig + 128 * ncomps + c);
        const auto *ig_129 = buffer.data(ig + 129 * ncomps + c);
        const auto *ig_130 = buffer.data(ig + 130 * ncomps + c);
        const auto *ig_131 = buffer.data(ig + 131 * ncomps + c);
        const auto *ig_132 = buffer.data(ig + 132 * ncomps + c);
        const auto *ig_133 = buffer.data(ig + 133 * ncomps + c);
        const auto *ig_134 = buffer.data(ig + 134 * ncomps + c);
        const auto *ig_135 = buffer.data(ig + 135 * ncomps + c);
        const auto *ig_136 = buffer.data(ig + 136 * ncomps + c);
        const auto *ig_137 = buffer.data(ig + 137 * ncomps + c);
        const auto *ig_138 = buffer.data(ig + 138 * ncomps + c);
        const auto *ig_139 = buffer.data(ig + 139 * ncomps + c);
        const auto *ig_140 = buffer.data(ig + 140 * ncomps + c);
        const auto *ig_141 = buffer.data(ig + 141 * ncomps + c);
        const auto *ig_142 = buffer.data(ig + 142 * ncomps + c);
        const auto *ig_143 = buffer.data(ig + 143 * ncomps + c);
        const auto *ig_144 = buffer.data(ig + 144 * ncomps + c);
        const auto *ig_145 = buffer.data(ig + 145 * ncomps + c);
        const auto *ig_146 = buffer.data(ig + 146 * ncomps + c);
        const auto *ig_147 = buffer.data(ig + 147 * ncomps + c);
        const auto *ig_148 = buffer.data(ig + 148 * ncomps + c);
        const auto *ig_149 = buffer.data(ig + 149 * ncomps + c);
        const auto *ig_150 = buffer.data(ig + 150 * ncomps + c);
        const auto *ig_151 = buffer.data(ig + 151 * ncomps + c);
        const auto *ig_152 = buffer.data(ig + 152 * ncomps + c);
        const auto *ig_153 = buffer.data(ig + 153 * ncomps + c);
        const auto *ig_154 = buffer.data(ig + 154 * ncomps + c);
        const auto *ig_155 = buffer.data(ig + 155 * ncomps + c);
        const auto *ig_156 = buffer.data(ig + 156 * ncomps + c);
        const auto *ig_157 = buffer.data(ig + 157 * ncomps + c);
        const auto *ig_158 = buffer.data(ig + 158 * ncomps + c);
        const auto *ig_159 = buffer.data(ig + 159 * ncomps + c);
        const auto *ig_160 = buffer.data(ig + 160 * ncomps + c);
        const auto *ig_161 = buffer.data(ig + 161 * ncomps + c);
        const auto *ig_162 = buffer.data(ig + 162 * ncomps + c);
        const auto *ig_163 = buffer.data(ig + 163 * ncomps + c);
        const auto *ig_164 = buffer.data(ig + 164 * ncomps + c);
        const auto *ig_165 = buffer.data(ig + 165 * ncomps + c);
        const auto *ig_166 = buffer.data(ig + 166 * ncomps + c);
        const auto *ig_167 = buffer.data(ig + 167 * ncomps + c);
        const auto *ig_168 = buffer.data(ig + 168 * ncomps + c);
        const auto *ig_169 = buffer.data(ig + 169 * ncomps + c);
        const auto *ig_170 = buffer.data(ig + 170 * ncomps + c);
        const auto *ig_171 = buffer.data(ig + 171 * ncomps + c);
        const auto *ig_172 = buffer.data(ig + 172 * ncomps + c);
        const auto *ig_173 = buffer.data(ig + 173 * ncomps + c);
        const auto *ig_174 = buffer.data(ig + 174 * ncomps + c);
        const auto *ig_175 = buffer.data(ig + 175 * ncomps + c);
        const auto *ig_176 = buffer.data(ig + 176 * ncomps + c);
        const auto *ig_177 = buffer.data(ig + 177 * ncomps + c);
        const auto *ig_178 = buffer.data(ig + 178 * ncomps + c);
        const auto *ig_179 = buffer.data(ig + 179 * ncomps + c);
        const auto *ig_180 = buffer.data(ig + 180 * ncomps + c);
        const auto *ig_181 = buffer.data(ig + 181 * ncomps + c);
        const auto *ig_182 = buffer.data(ig + 182 * ncomps + c);
        const auto *ig_183 = buffer.data(ig + 183 * ncomps + c);
        const auto *ig_184 = buffer.data(ig + 184 * ncomps + c);
        const auto *ig_185 = buffer.data(ig + 185 * ncomps + c);
        const auto *ig_186 = buffer.data(ig + 186 * ncomps + c);
        const auto *ig_187 = buffer.data(ig + 187 * ncomps + c);
        const auto *ig_188 = buffer.data(ig + 188 * ncomps + c);
        const auto *ig_189 = buffer.data(ig + 189 * ncomps + c);
        const auto *ig_190 = buffer.data(ig + 190 * ncomps + c);
        const auto *ig_191 = buffer.data(ig + 191 * ncomps + c);
        const auto *ig_192 = buffer.data(ig + 192 * ncomps + c);
        const auto *ig_193 = buffer.data(ig + 193 * ncomps + c);
        const auto *ig_194 = buffer.data(ig + 194 * ncomps + c);
        const auto *ig_195 = buffer.data(ig + 195 * ncomps + c);
        const auto *ig_196 = buffer.data(ig + 196 * ncomps + c);
        const auto *ig_197 = buffer.data(ig + 197 * ncomps + c);
        const auto *ig_198 = buffer.data(ig + 198 * ncomps + c);
        const auto *ig_199 = buffer.data(ig + 199 * ncomps + c);
        const auto *ig_200 = buffer.data(ig + 200 * ncomps + c);
        const auto *ig_201 = buffer.data(ig + 201 * ncomps + c);
        const auto *ig_202 = buffer.data(ig + 202 * ncomps + c);

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
        const auto *kg_205 = buffer.data(kg + 205 * ncomps + c);
        const auto *kg_206 = buffer.data(kg + 206 * ncomps + c);
        const auto *kg_207 = buffer.data(kg + 207 * ncomps + c);
        const auto *kg_208 = buffer.data(kg + 208 * ncomps + c);
        const auto *kg_209 = buffer.data(kg + 209 * ncomps + c);
        const auto *kg_224 = buffer.data(kg + 224 * ncomps + c);
        const auto *kg_235 = buffer.data(kg + 235 * ncomps + c);
        const auto *kg_236 = buffer.data(kg + 236 * ncomps + c);
        const auto *kg_237 = buffer.data(kg + 237 * ncomps + c);
        const auto *kg_238 = buffer.data(kg + 238 * ncomps + c);
        const auto *kg_239 = buffer.data(kg + 239 * ncomps + c);
        const auto *kg_250 = buffer.data(kg + 250 * ncomps + c);
        const auto *kg_251 = buffer.data(kg + 251 * ncomps + c);
        const auto *kg_252 = buffer.data(kg + 252 * ncomps + c);
        const auto *kg_253 = buffer.data(kg + 253 * ncomps + c);
        const auto *kg_254 = buffer.data(kg + 254 * ncomps + c);
        const auto *kg_265 = buffer.data(kg + 265 * ncomps + c);
        const auto *kg_266 = buffer.data(kg + 266 * ncomps + c);
        const auto *kg_267 = buffer.data(kg + 267 * ncomps + c);
        const auto *kg_268 = buffer.data(kg + 268 * ncomps + c);
        const auto *kg_269 = buffer.data(kg + 269 * ncomps + c);
        const auto *kg_284 = buffer.data(kg + 284 * ncomps + c);

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, ab_y, ig_100, ig_101, ig_102, \
                         ig_103, ig_104, kg_160, kg_161, kg_162, kg_163, \
                         kg_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_141[k] = ab_y[k] * ig_100[k]
                       + kg_160[k];

            t_142[k] = ab_y[k] * ig_101[k]
                       + kg_161[k];

            t_143[k] = ab_y[k] * ig_102[k]
                       + kg_162[k];

            t_144[k] = ab_y[k] * ig_103[k]
                       + kg_163[k];

            t_145[k] = ab_y[k] * ig_104[k]
                       + kg_164[k];
        }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, ab_x, ab_z, ig_104, ig_105, ig_106, \
                         ig_107, kg_105, kg_106, kg_107, kg_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_146[k] = ab_z[k] * ig_104[k]
                       + kg_179[k];

            t_147[k] = ab_x[k] * ig_105[k]
                       + kg_105[k];

            t_148[k] = ab_x[k] * ig_106[k]
                       + kg_106[k];

            t_149[k] = ab_x[k] * ig_107[k]
                       + kg_107[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, ig_108, ig_109, ig_110, \
                         ig_111, ig_112, kg_108, kg_109, kg_110, kg_111, \
                         kg_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * ig_108[k]
                       + kg_108[k];

            t_151[k] = ab_x[k] * ig_109[k]
                       + kg_109[k];

            t_152[k] = ab_x[k] * ig_110[k]
                       + kg_110[k];

            t_153[k] = ab_x[k] * ig_111[k]
                       + kg_111[k];

            t_154[k] = ab_x[k] * ig_112[k]
                       + kg_112[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, ig_113, ig_114, ig_115, \
                         ig_116, ig_117, kg_113, kg_114, kg_115, kg_116, \
                         kg_117 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * ig_113[k]
                       + kg_113[k];

            t_156[k] = ab_x[k] * ig_114[k]
                       + kg_114[k];

            t_157[k] = ab_x[k] * ig_115[k]
                       + kg_115[k];

            t_158[k] = ab_x[k] * ig_116[k]
                       + kg_116[k];

            t_159[k] = ab_x[k] * ig_117[k]
                       + kg_117[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, ab_x, ab_y, ig_115, ig_116, ig_118, \
                         ig_119, kg_118, kg_119, kg_175, kg_176 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_x[k] * ig_118[k]
                       + kg_118[k];

            t_161[k] = ab_x[k] * ig_119[k]
                       + kg_119[k];

            t_162[k] = ab_y[k] * ig_115[k]
                       + kg_175[k];

            t_163[k] = ab_y[k] * ig_116[k]
                       + kg_176[k];
        }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, ab_y, ab_z, ig_117, ig_118, ig_119, \
                         kg_177, kg_178, kg_179, kg_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_164[k] = ab_y[k] * ig_117[k]
                       + kg_177[k];

            t_165[k] = ab_y[k] * ig_118[k]
                       + kg_178[k];

            t_166[k] = ab_y[k] * ig_119[k]
                       + kg_179[k];

            t_167[k] = ab_z[k] * ig_119[k]
                       + kg_194[k];
        }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, ab_x, ig_120, ig_121, ig_122, \
                         ig_123, ig_124, kg_120, kg_121, kg_122, kg_123, \
                         kg_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_168[k] = ab_x[k] * ig_120[k]
                       + kg_120[k];

            t_169[k] = ab_x[k] * ig_121[k]
                       + kg_121[k];

            t_170[k] = ab_x[k] * ig_122[k]
                       + kg_122[k];

            t_171[k] = ab_x[k] * ig_123[k]
                       + kg_123[k];

            t_172[k] = ab_x[k] * ig_124[k]
                       + kg_124[k];
        }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, ab_x, ig_125, ig_126, ig_127, \
                         ig_128, ig_129, kg_125, kg_126, kg_127, kg_128, \
                         kg_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_173[k] = ab_x[k] * ig_125[k]
                       + kg_125[k];

            t_174[k] = ab_x[k] * ig_126[k]
                       + kg_126[k];

            t_175[k] = ab_x[k] * ig_127[k]
                       + kg_127[k];

            t_176[k] = ab_x[k] * ig_128[k]
                       + kg_128[k];

            t_177[k] = ab_x[k] * ig_129[k]
                       + kg_129[k];
        }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, ab_x, ig_130, ig_131, ig_132, \
                         ig_133, ig_134, kg_130, kg_131, kg_132, kg_133, \
                         kg_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_178[k] = ab_x[k] * ig_130[k]
                       + kg_130[k];

            t_179[k] = ab_x[k] * ig_131[k]
                       + kg_131[k];

            t_180[k] = ab_x[k] * ig_132[k]
                       + kg_132[k];

            t_181[k] = ab_x[k] * ig_133[k]
                       + kg_133[k];

            t_182[k] = ab_x[k] * ig_134[k]
                       + kg_134[k];
        }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, ab_y, ig_130, ig_131, ig_132, \
                         ig_133, ig_134, kg_190, kg_191, kg_192, kg_193, \
                         kg_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_183[k] = ab_y[k] * ig_130[k]
                       + kg_190[k];

            t_184[k] = ab_y[k] * ig_131[k]
                       + kg_191[k];

            t_185[k] = ab_y[k] * ig_132[k]
                       + kg_192[k];

            t_186[k] = ab_y[k] * ig_133[k]
                       + kg_193[k];

            t_187[k] = ab_y[k] * ig_134[k]
                       + kg_194[k];
        }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, ab_x, ab_z, ig_134, ig_135, ig_136, \
                         ig_137, kg_135, kg_136, kg_137, kg_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_188[k] = ab_z[k] * ig_134[k]
                       + kg_209[k];

            t_189[k] = ab_x[k] * ig_135[k]
                       + kg_135[k];

            t_190[k] = ab_x[k] * ig_136[k]
                       + kg_136[k];

            t_191[k] = ab_x[k] * ig_137[k]
                       + kg_137[k];
        }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, ab_x, ig_138, ig_139, ig_140, \
                         ig_141, ig_142, kg_138, kg_139, kg_140, kg_141, \
                         kg_142 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_192[k] = ab_x[k] * ig_138[k]
                       + kg_138[k];

            t_193[k] = ab_x[k] * ig_139[k]
                       + kg_139[k];

            t_194[k] = ab_x[k] * ig_140[k]
                       + kg_140[k];

            t_195[k] = ab_x[k] * ig_141[k]
                       + kg_141[k];

            t_196[k] = ab_x[k] * ig_142[k]
                       + kg_142[k];
        }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, ab_x, ig_143, ig_144, ig_145, \
                         ig_146, ig_147, kg_143, kg_144, kg_145, kg_146, \
                         kg_147 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_197[k] = ab_x[k] * ig_143[k]
                       + kg_143[k];

            t_198[k] = ab_x[k] * ig_144[k]
                       + kg_144[k];

            t_199[k] = ab_x[k] * ig_145[k]
                       + kg_145[k];

            t_200[k] = ab_x[k] * ig_146[k]
                       + kg_146[k];

            t_201[k] = ab_x[k] * ig_147[k]
                       + kg_147[k];
        }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, ab_x, ab_y, ig_145, ig_146, ig_148, \
                         ig_149, kg_148, kg_149, kg_205, kg_206 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_202[k] = ab_x[k] * ig_148[k]
                       + kg_148[k];

            t_203[k] = ab_x[k] * ig_149[k]
                       + kg_149[k];

            t_204[k] = ab_y[k] * ig_145[k]
                       + kg_205[k];

            t_205[k] = ab_y[k] * ig_146[k]
                       + kg_206[k];
        }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, ab_y, ab_z, ig_147, ig_148, ig_149, \
                         kg_207, kg_208, kg_209, kg_224 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_206[k] = ab_y[k] * ig_147[k]
                       + kg_207[k];

            t_207[k] = ab_y[k] * ig_148[k]
                       + kg_208[k];

            t_208[k] = ab_y[k] * ig_149[k]
                       + kg_209[k];

            t_209[k] = ab_z[k] * ig_149[k]
                       + kg_224[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, ig_150, ig_151, ig_152, \
                         ig_153, ig_154, kg_150, kg_151, kg_152, kg_153, \
                         kg_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = ab_x[k] * ig_150[k]
                       + kg_150[k];

            t_211[k] = ab_x[k] * ig_151[k]
                       + kg_151[k];

            t_212[k] = ab_x[k] * ig_152[k]
                       + kg_152[k];

            t_213[k] = ab_x[k] * ig_153[k]
                       + kg_153[k];

            t_214[k] = ab_x[k] * ig_154[k]
                       + kg_154[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, ig_155, ig_156, ig_157, \
                         ig_158, ig_159, kg_155, kg_156, kg_157, kg_158, \
                         kg_159 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = ab_x[k] * ig_155[k]
                       + kg_155[k];

            t_216[k] = ab_x[k] * ig_156[k]
                       + kg_156[k];

            t_217[k] = ab_x[k] * ig_157[k]
                       + kg_157[k];

            t_218[k] = ab_x[k] * ig_158[k]
                       + kg_158[k];

            t_219[k] = ab_x[k] * ig_159[k]
                       + kg_159[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, ig_160, ig_161, ig_162, \
                         ig_163, ig_164, kg_160, kg_161, kg_162, kg_163, \
                         kg_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = ab_x[k] * ig_160[k]
                       + kg_160[k];

            t_221[k] = ab_x[k] * ig_161[k]
                       + kg_161[k];

            t_222[k] = ab_x[k] * ig_162[k]
                       + kg_162[k];

            t_223[k] = ab_x[k] * ig_163[k]
                       + kg_163[k];

            t_224[k] = ab_x[k] * ig_164[k]
                       + kg_164[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_y, ig_160, ig_161, ig_162, \
                         ig_163, ig_164, kg_235, kg_236, kg_237, kg_238, \
                         kg_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = ab_y[k] * ig_160[k]
                       + kg_235[k];

            t_226[k] = ab_y[k] * ig_161[k]
                       + kg_236[k];

            t_227[k] = ab_y[k] * ig_162[k]
                       + kg_237[k];

            t_228[k] = ab_y[k] * ig_163[k]
                       + kg_238[k];

            t_229[k] = ab_y[k] * ig_164[k]
                       + kg_239[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, ab_x, ab_z, ig_164, ig_165, ig_166, \
                         ig_167, kg_165, kg_166, kg_167, kg_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = ab_z[k] * ig_164[k]
                       + kg_254[k];

            t_231[k] = ab_x[k] * ig_165[k]
                       + kg_165[k];

            t_232[k] = ab_x[k] * ig_166[k]
                       + kg_166[k];

            t_233[k] = ab_x[k] * ig_167[k]
                       + kg_167[k];
        }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, ab_x, ig_168, ig_169, ig_170, \
                         ig_171, ig_172, kg_168, kg_169, kg_170, kg_171, \
                         kg_172 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_234[k] = ab_x[k] * ig_168[k]
                       + kg_168[k];

            t_235[k] = ab_x[k] * ig_169[k]
                       + kg_169[k];

            t_236[k] = ab_x[k] * ig_170[k]
                       + kg_170[k];

            t_237[k] = ab_x[k] * ig_171[k]
                       + kg_171[k];

            t_238[k] = ab_x[k] * ig_172[k]
                       + kg_172[k];
        }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, ab_x, ig_173, ig_174, ig_175, \
                         ig_176, ig_177, kg_173, kg_174, kg_175, kg_176, \
                         kg_177 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_239[k] = ab_x[k] * ig_173[k]
                       + kg_173[k];

            t_240[k] = ab_x[k] * ig_174[k]
                       + kg_174[k];

            t_241[k] = ab_x[k] * ig_175[k]
                       + kg_175[k];

            t_242[k] = ab_x[k] * ig_176[k]
                       + kg_176[k];

            t_243[k] = ab_x[k] * ig_177[k]
                       + kg_177[k];
        }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, ab_x, ab_y, ig_175, ig_176, ig_178, \
                         ig_179, kg_178, kg_179, kg_250, kg_251 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_244[k] = ab_x[k] * ig_178[k]
                       + kg_178[k];

            t_245[k] = ab_x[k] * ig_179[k]
                       + kg_179[k];

            t_246[k] = ab_y[k] * ig_175[k]
                       + kg_250[k];

            t_247[k] = ab_y[k] * ig_176[k]
                       + kg_251[k];
        }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, ab_y, ab_z, ig_177, ig_178, ig_179, \
                         kg_252, kg_253, kg_254, kg_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_248[k] = ab_y[k] * ig_177[k]
                       + kg_252[k];

            t_249[k] = ab_y[k] * ig_178[k]
                       + kg_253[k];

            t_250[k] = ab_y[k] * ig_179[k]
                       + kg_254[k];

            t_251[k] = ab_z[k] * ig_179[k]
                       + kg_269[k];
        }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, ab_x, ig_180, ig_181, ig_182, \
                         ig_183, ig_184, kg_180, kg_181, kg_182, kg_183, \
                         kg_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_252[k] = ab_x[k] * ig_180[k]
                       + kg_180[k];

            t_253[k] = ab_x[k] * ig_181[k]
                       + kg_181[k];

            t_254[k] = ab_x[k] * ig_182[k]
                       + kg_182[k];

            t_255[k] = ab_x[k] * ig_183[k]
                       + kg_183[k];

            t_256[k] = ab_x[k] * ig_184[k]
                       + kg_184[k];
        }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, ab_x, ig_185, ig_186, ig_187, \
                         ig_188, ig_189, kg_185, kg_186, kg_187, kg_188, \
                         kg_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_257[k] = ab_x[k] * ig_185[k]
                       + kg_185[k];

            t_258[k] = ab_x[k] * ig_186[k]
                       + kg_186[k];

            t_259[k] = ab_x[k] * ig_187[k]
                       + kg_187[k];

            t_260[k] = ab_x[k] * ig_188[k]
                       + kg_188[k];

            t_261[k] = ab_x[k] * ig_189[k]
                       + kg_189[k];
        }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, ab_x, ig_190, ig_191, ig_192, \
                         ig_193, ig_194, kg_190, kg_191, kg_192, kg_193, \
                         kg_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_262[k] = ab_x[k] * ig_190[k]
                       + kg_190[k];

            t_263[k] = ab_x[k] * ig_191[k]
                       + kg_191[k];

            t_264[k] = ab_x[k] * ig_192[k]
                       + kg_192[k];

            t_265[k] = ab_x[k] * ig_193[k]
                       + kg_193[k];

            t_266[k] = ab_x[k] * ig_194[k]
                       + kg_194[k];
        }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, ab_y, ig_190, ig_191, ig_192, \
                         ig_193, ig_194, kg_265, kg_266, kg_267, kg_268, \
                         kg_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_267[k] = ab_y[k] * ig_190[k]
                       + kg_265[k];

            t_268[k] = ab_y[k] * ig_191[k]
                       + kg_266[k];

            t_269[k] = ab_y[k] * ig_192[k]
                       + kg_267[k];

            t_270[k] = ab_y[k] * ig_193[k]
                       + kg_268[k];

            t_271[k] = ab_y[k] * ig_194[k]
                       + kg_269[k];
        }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, ab_x, ab_z, ig_194, ig_195, ig_196, \
                         ig_197, kg_195, kg_196, kg_197, kg_284 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_272[k] = ab_z[k] * ig_194[k]
                       + kg_284[k];

            t_273[k] = ab_x[k] * ig_195[k]
                       + kg_195[k];

            t_274[k] = ab_x[k] * ig_196[k]
                       + kg_196[k];

            t_275[k] = ab_x[k] * ig_197[k]
                       + kg_197[k];
        }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, ab_x, ig_198, ig_199, ig_200, \
                         ig_201, ig_202, kg_198, kg_199, kg_200, kg_201, \
                         kg_202 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_276[k] = ab_x[k] * ig_198[k]
                       + kg_198[k];

            t_277[k] = ab_x[k] * ig_199[k]
                       + kg_199[k];

            t_278[k] = ab_x[k] * ig_200[k]
                       + kg_200[k];

            t_279[k] = ab_x[k] * ig_201[k]
                       + kg_201[k];

            t_280[k] = ab_x[k] * ig_202[k]
                       + kg_202[k];
        }
    }
}

static auto
compute_hrr_ih_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t ig, const size_t kg, const size_t ncomps,
                      const size_t nmax) -> void
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

        const auto *ig_203 = buffer.data(ig + 203 * ncomps + c);
        const auto *ig_204 = buffer.data(ig + 204 * ncomps + c);
        const auto *ig_205 = buffer.data(ig + 205 * ncomps + c);
        const auto *ig_206 = buffer.data(ig + 206 * ncomps + c);
        const auto *ig_207 = buffer.data(ig + 207 * ncomps + c);
        const auto *ig_208 = buffer.data(ig + 208 * ncomps + c);
        const auto *ig_209 = buffer.data(ig + 209 * ncomps + c);
        const auto *ig_210 = buffer.data(ig + 210 * ncomps + c);
        const auto *ig_211 = buffer.data(ig + 211 * ncomps + c);
        const auto *ig_212 = buffer.data(ig + 212 * ncomps + c);
        const auto *ig_213 = buffer.data(ig + 213 * ncomps + c);
        const auto *ig_214 = buffer.data(ig + 214 * ncomps + c);
        const auto *ig_215 = buffer.data(ig + 215 * ncomps + c);
        const auto *ig_216 = buffer.data(ig + 216 * ncomps + c);
        const auto *ig_217 = buffer.data(ig + 217 * ncomps + c);
        const auto *ig_218 = buffer.data(ig + 218 * ncomps + c);
        const auto *ig_219 = buffer.data(ig + 219 * ncomps + c);
        const auto *ig_220 = buffer.data(ig + 220 * ncomps + c);
        const auto *ig_221 = buffer.data(ig + 221 * ncomps + c);
        const auto *ig_222 = buffer.data(ig + 222 * ncomps + c);
        const auto *ig_223 = buffer.data(ig + 223 * ncomps + c);
        const auto *ig_224 = buffer.data(ig + 224 * ncomps + c);
        const auto *ig_225 = buffer.data(ig + 225 * ncomps + c);
        const auto *ig_226 = buffer.data(ig + 226 * ncomps + c);
        const auto *ig_227 = buffer.data(ig + 227 * ncomps + c);
        const auto *ig_228 = buffer.data(ig + 228 * ncomps + c);
        const auto *ig_229 = buffer.data(ig + 229 * ncomps + c);
        const auto *ig_230 = buffer.data(ig + 230 * ncomps + c);
        const auto *ig_231 = buffer.data(ig + 231 * ncomps + c);
        const auto *ig_232 = buffer.data(ig + 232 * ncomps + c);
        const auto *ig_233 = buffer.data(ig + 233 * ncomps + c);
        const auto *ig_234 = buffer.data(ig + 234 * ncomps + c);
        const auto *ig_235 = buffer.data(ig + 235 * ncomps + c);
        const auto *ig_236 = buffer.data(ig + 236 * ncomps + c);
        const auto *ig_237 = buffer.data(ig + 237 * ncomps + c);
        const auto *ig_238 = buffer.data(ig + 238 * ncomps + c);
        const auto *ig_239 = buffer.data(ig + 239 * ncomps + c);
        const auto *ig_240 = buffer.data(ig + 240 * ncomps + c);
        const auto *ig_241 = buffer.data(ig + 241 * ncomps + c);
        const auto *ig_242 = buffer.data(ig + 242 * ncomps + c);
        const auto *ig_243 = buffer.data(ig + 243 * ncomps + c);
        const auto *ig_244 = buffer.data(ig + 244 * ncomps + c);
        const auto *ig_245 = buffer.data(ig + 245 * ncomps + c);
        const auto *ig_246 = buffer.data(ig + 246 * ncomps + c);
        const auto *ig_247 = buffer.data(ig + 247 * ncomps + c);
        const auto *ig_248 = buffer.data(ig + 248 * ncomps + c);
        const auto *ig_249 = buffer.data(ig + 249 * ncomps + c);
        const auto *ig_250 = buffer.data(ig + 250 * ncomps + c);
        const auto *ig_251 = buffer.data(ig + 251 * ncomps + c);
        const auto *ig_252 = buffer.data(ig + 252 * ncomps + c);
        const auto *ig_253 = buffer.data(ig + 253 * ncomps + c);
        const auto *ig_254 = buffer.data(ig + 254 * ncomps + c);
        const auto *ig_255 = buffer.data(ig + 255 * ncomps + c);
        const auto *ig_256 = buffer.data(ig + 256 * ncomps + c);
        const auto *ig_257 = buffer.data(ig + 257 * ncomps + c);
        const auto *ig_258 = buffer.data(ig + 258 * ncomps + c);
        const auto *ig_259 = buffer.data(ig + 259 * ncomps + c);
        const auto *ig_260 = buffer.data(ig + 260 * ncomps + c);
        const auto *ig_261 = buffer.data(ig + 261 * ncomps + c);
        const auto *ig_262 = buffer.data(ig + 262 * ncomps + c);
        const auto *ig_263 = buffer.data(ig + 263 * ncomps + c);
        const auto *ig_264 = buffer.data(ig + 264 * ncomps + c);
        const auto *ig_265 = buffer.data(ig + 265 * ncomps + c);
        const auto *ig_266 = buffer.data(ig + 266 * ncomps + c);
        const auto *ig_267 = buffer.data(ig + 267 * ncomps + c);
        const auto *ig_268 = buffer.data(ig + 268 * ncomps + c);
        const auto *ig_269 = buffer.data(ig + 269 * ncomps + c);
        const auto *ig_270 = buffer.data(ig + 270 * ncomps + c);
        const auto *ig_271 = buffer.data(ig + 271 * ncomps + c);
        const auto *ig_272 = buffer.data(ig + 272 * ncomps + c);
        const auto *ig_273 = buffer.data(ig + 273 * ncomps + c);
        const auto *ig_274 = buffer.data(ig + 274 * ncomps + c);
        const auto *ig_275 = buffer.data(ig + 275 * ncomps + c);
        const auto *ig_276 = buffer.data(ig + 276 * ncomps + c);
        const auto *ig_277 = buffer.data(ig + 277 * ncomps + c);
        const auto *ig_278 = buffer.data(ig + 278 * ncomps + c);
        const auto *ig_279 = buffer.data(ig + 279 * ncomps + c);
        const auto *ig_280 = buffer.data(ig + 280 * ncomps + c);
        const auto *ig_281 = buffer.data(ig + 281 * ncomps + c);
        const auto *ig_282 = buffer.data(ig + 282 * ncomps + c);
        const auto *ig_283 = buffer.data(ig + 283 * ncomps + c);
        const auto *ig_284 = buffer.data(ig + 284 * ncomps + c);
        const auto *ig_285 = buffer.data(ig + 285 * ncomps + c);
        const auto *ig_286 = buffer.data(ig + 286 * ncomps + c);
        const auto *ig_287 = buffer.data(ig + 287 * ncomps + c);
        const auto *ig_288 = buffer.data(ig + 288 * ncomps + c);
        const auto *ig_289 = buffer.data(ig + 289 * ncomps + c);
        const auto *ig_290 = buffer.data(ig + 290 * ncomps + c);
        const auto *ig_291 = buffer.data(ig + 291 * ncomps + c);
        const auto *ig_292 = buffer.data(ig + 292 * ncomps + c);
        const auto *ig_293 = buffer.data(ig + 293 * ncomps + c);
        const auto *ig_294 = buffer.data(ig + 294 * ncomps + c);
        const auto *ig_295 = buffer.data(ig + 295 * ncomps + c);
        const auto *ig_296 = buffer.data(ig + 296 * ncomps + c);
        const auto *ig_297 = buffer.data(ig + 297 * ncomps + c);
        const auto *ig_298 = buffer.data(ig + 298 * ncomps + c);
        const auto *ig_299 = buffer.data(ig + 299 * ncomps + c);
        const auto *ig_300 = buffer.data(ig + 300 * ncomps + c);
        const auto *ig_301 = buffer.data(ig + 301 * ncomps + c);
        const auto *ig_302 = buffer.data(ig + 302 * ncomps + c);
        const auto *ig_303 = buffer.data(ig + 303 * ncomps + c);
        const auto *ig_304 = buffer.data(ig + 304 * ncomps + c);

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
        const auto *kg_314 = buffer.data(kg + 314 * ncomps + c);
        const auto *kg_325 = buffer.data(kg + 325 * ncomps + c);
        const auto *kg_326 = buffer.data(kg + 326 * ncomps + c);
        const auto *kg_327 = buffer.data(kg + 327 * ncomps + c);
        const auto *kg_328 = buffer.data(kg + 328 * ncomps + c);
        const auto *kg_329 = buffer.data(kg + 329 * ncomps + c);
        const auto *kg_340 = buffer.data(kg + 340 * ncomps + c);
        const auto *kg_341 = buffer.data(kg + 341 * ncomps + c);
        const auto *kg_342 = buffer.data(kg + 342 * ncomps + c);
        const auto *kg_343 = buffer.data(kg + 343 * ncomps + c);
        const auto *kg_344 = buffer.data(kg + 344 * ncomps + c);
        const auto *kg_355 = buffer.data(kg + 355 * ncomps + c);
        const auto *kg_356 = buffer.data(kg + 356 * ncomps + c);
        const auto *kg_357 = buffer.data(kg + 357 * ncomps + c);
        const auto *kg_358 = buffer.data(kg + 358 * ncomps + c);
        const auto *kg_359 = buffer.data(kg + 359 * ncomps + c);
        const auto *kg_370 = buffer.data(kg + 370 * ncomps + c);
        const auto *kg_371 = buffer.data(kg + 371 * ncomps + c);
        const auto *kg_372 = buffer.data(kg + 372 * ncomps + c);
        const auto *kg_373 = buffer.data(kg + 373 * ncomps + c);
        const auto *kg_374 = buffer.data(kg + 374 * ncomps + c);
        const auto *kg_385 = buffer.data(kg + 385 * ncomps + c);
        const auto *kg_386 = buffer.data(kg + 386 * ncomps + c);
        const auto *kg_387 = buffer.data(kg + 387 * ncomps + c);
        const auto *kg_388 = buffer.data(kg + 388 * ncomps + c);
        const auto *kg_389 = buffer.data(kg + 389 * ncomps + c);
        const auto *kg_404 = buffer.data(kg + 404 * ncomps + c);

#pragma omp simd aligned(t_281, t_282, t_283, t_284, t_285, ab_x, ig_203, ig_204, ig_205, \
                         ig_206, ig_207, kg_203, kg_204, kg_205, kg_206, \
                         kg_207 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_281[k] = ab_x[k] * ig_203[k]
                       + kg_203[k];

            t_282[k] = ab_x[k] * ig_204[k]
                       + kg_204[k];

            t_283[k] = ab_x[k] * ig_205[k]
                       + kg_205[k];

            t_284[k] = ab_x[k] * ig_206[k]
                       + kg_206[k];

            t_285[k] = ab_x[k] * ig_207[k]
                       + kg_207[k];
        }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, ab_x, ab_y, ig_205, ig_206, ig_208, \
                         ig_209, kg_208, kg_209, kg_280, kg_281 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_286[k] = ab_x[k] * ig_208[k]
                       + kg_208[k];

            t_287[k] = ab_x[k] * ig_209[k]
                       + kg_209[k];

            t_288[k] = ab_y[k] * ig_205[k]
                       + kg_280[k];

            t_289[k] = ab_y[k] * ig_206[k]
                       + kg_281[k];
        }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, ab_y, ab_z, ig_207, ig_208, ig_209, \
                         kg_282, kg_283, kg_284, kg_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = ab_y[k] * ig_207[k]
                       + kg_282[k];

            t_291[k] = ab_y[k] * ig_208[k]
                       + kg_283[k];

            t_292[k] = ab_y[k] * ig_209[k]
                       + kg_284[k];

            t_293[k] = ab_z[k] * ig_209[k]
                       + kg_299[k];
        }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, ab_x, ig_210, ig_211, ig_212, \
                         ig_213, ig_214, kg_210, kg_211, kg_212, kg_213, \
                         kg_214 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_294[k] = ab_x[k] * ig_210[k]
                       + kg_210[k];

            t_295[k] = ab_x[k] * ig_211[k]
                       + kg_211[k];

            t_296[k] = ab_x[k] * ig_212[k]
                       + kg_212[k];

            t_297[k] = ab_x[k] * ig_213[k]
                       + kg_213[k];

            t_298[k] = ab_x[k] * ig_214[k]
                       + kg_214[k];
        }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, ab_x, ig_215, ig_216, ig_217, \
                         ig_218, ig_219, kg_215, kg_216, kg_217, kg_218, \
                         kg_219 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_299[k] = ab_x[k] * ig_215[k]
                       + kg_215[k];

            t_300[k] = ab_x[k] * ig_216[k]
                       + kg_216[k];

            t_301[k] = ab_x[k] * ig_217[k]
                       + kg_217[k];

            t_302[k] = ab_x[k] * ig_218[k]
                       + kg_218[k];

            t_303[k] = ab_x[k] * ig_219[k]
                       + kg_219[k];
        }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, ab_x, ig_220, ig_221, ig_222, \
                         ig_223, ig_224, kg_220, kg_221, kg_222, kg_223, \
                         kg_224 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_304[k] = ab_x[k] * ig_220[k]
                       + kg_220[k];

            t_305[k] = ab_x[k] * ig_221[k]
                       + kg_221[k];

            t_306[k] = ab_x[k] * ig_222[k]
                       + kg_222[k];

            t_307[k] = ab_x[k] * ig_223[k]
                       + kg_223[k];

            t_308[k] = ab_x[k] * ig_224[k]
                       + kg_224[k];
        }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, ab_y, ig_220, ig_221, ig_222, \
                         ig_223, ig_224, kg_295, kg_296, kg_297, kg_298, \
                         kg_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_309[k] = ab_y[k] * ig_220[k]
                       + kg_295[k];

            t_310[k] = ab_y[k] * ig_221[k]
                       + kg_296[k];

            t_311[k] = ab_y[k] * ig_222[k]
                       + kg_297[k];

            t_312[k] = ab_y[k] * ig_223[k]
                       + kg_298[k];

            t_313[k] = ab_y[k] * ig_224[k]
                       + kg_299[k];
        }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, ab_x, ab_z, ig_224, ig_225, ig_226, \
                         ig_227, kg_225, kg_226, kg_227, kg_314 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_314[k] = ab_z[k] * ig_224[k]
                       + kg_314[k];

            t_315[k] = ab_x[k] * ig_225[k]
                       + kg_225[k];

            t_316[k] = ab_x[k] * ig_226[k]
                       + kg_226[k];

            t_317[k] = ab_x[k] * ig_227[k]
                       + kg_227[k];
        }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, ab_x, ig_228, ig_229, ig_230, \
                         ig_231, ig_232, kg_228, kg_229, kg_230, kg_231, \
                         kg_232 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_318[k] = ab_x[k] * ig_228[k]
                       + kg_228[k];

            t_319[k] = ab_x[k] * ig_229[k]
                       + kg_229[k];

            t_320[k] = ab_x[k] * ig_230[k]
                       + kg_230[k];

            t_321[k] = ab_x[k] * ig_231[k]
                       + kg_231[k];

            t_322[k] = ab_x[k] * ig_232[k]
                       + kg_232[k];
        }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, ab_x, ig_233, ig_234, ig_235, \
                         ig_236, ig_237, kg_233, kg_234, kg_235, kg_236, \
                         kg_237 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_323[k] = ab_x[k] * ig_233[k]
                       + kg_233[k];

            t_324[k] = ab_x[k] * ig_234[k]
                       + kg_234[k];

            t_325[k] = ab_x[k] * ig_235[k]
                       + kg_235[k];

            t_326[k] = ab_x[k] * ig_236[k]
                       + kg_236[k];

            t_327[k] = ab_x[k] * ig_237[k]
                       + kg_237[k];
        }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, ab_x, ab_y, ig_235, ig_236, ig_238, \
                         ig_239, kg_238, kg_239, kg_325, kg_326 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_328[k] = ab_x[k] * ig_238[k]
                       + kg_238[k];

            t_329[k] = ab_x[k] * ig_239[k]
                       + kg_239[k];

            t_330[k] = ab_y[k] * ig_235[k]
                       + kg_325[k];

            t_331[k] = ab_y[k] * ig_236[k]
                       + kg_326[k];
        }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, ab_y, ab_z, ig_237, ig_238, ig_239, \
                         kg_327, kg_328, kg_329, kg_344 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_332[k] = ab_y[k] * ig_237[k]
                       + kg_327[k];

            t_333[k] = ab_y[k] * ig_238[k]
                       + kg_328[k];

            t_334[k] = ab_y[k] * ig_239[k]
                       + kg_329[k];

            t_335[k] = ab_z[k] * ig_239[k]
                       + kg_344[k];
        }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, t_340, ab_x, ig_240, ig_241, ig_242, \
                         ig_243, ig_244, kg_240, kg_241, kg_242, kg_243, \
                         kg_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_336[k] = ab_x[k] * ig_240[k]
                       + kg_240[k];

            t_337[k] = ab_x[k] * ig_241[k]
                       + kg_241[k];

            t_338[k] = ab_x[k] * ig_242[k]
                       + kg_242[k];

            t_339[k] = ab_x[k] * ig_243[k]
                       + kg_243[k];

            t_340[k] = ab_x[k] * ig_244[k]
                       + kg_244[k];
        }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, t_345, ab_x, ig_245, ig_246, ig_247, \
                         ig_248, ig_249, kg_245, kg_246, kg_247, kg_248, \
                         kg_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_341[k] = ab_x[k] * ig_245[k]
                       + kg_245[k];

            t_342[k] = ab_x[k] * ig_246[k]
                       + kg_246[k];

            t_343[k] = ab_x[k] * ig_247[k]
                       + kg_247[k];

            t_344[k] = ab_x[k] * ig_248[k]
                       + kg_248[k];

            t_345[k] = ab_x[k] * ig_249[k]
                       + kg_249[k];
        }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, ab_x, ig_250, ig_251, ig_252, \
                         ig_253, ig_254, kg_250, kg_251, kg_252, kg_253, \
                         kg_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_346[k] = ab_x[k] * ig_250[k]
                       + kg_250[k];

            t_347[k] = ab_x[k] * ig_251[k]
                       + kg_251[k];

            t_348[k] = ab_x[k] * ig_252[k]
                       + kg_252[k];

            t_349[k] = ab_x[k] * ig_253[k]
                       + kg_253[k];

            t_350[k] = ab_x[k] * ig_254[k]
                       + kg_254[k];
        }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, ab_y, ig_250, ig_251, ig_252, \
                         ig_253, ig_254, kg_340, kg_341, kg_342, kg_343, \
                         kg_344 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_351[k] = ab_y[k] * ig_250[k]
                       + kg_340[k];

            t_352[k] = ab_y[k] * ig_251[k]
                       + kg_341[k];

            t_353[k] = ab_y[k] * ig_252[k]
                       + kg_342[k];

            t_354[k] = ab_y[k] * ig_253[k]
                       + kg_343[k];

            t_355[k] = ab_y[k] * ig_254[k]
                       + kg_344[k];
        }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, ab_x, ab_z, ig_254, ig_255, ig_256, \
                         ig_257, kg_255, kg_256, kg_257, kg_359 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_356[k] = ab_z[k] * ig_254[k]
                       + kg_359[k];

            t_357[k] = ab_x[k] * ig_255[k]
                       + kg_255[k];

            t_358[k] = ab_x[k] * ig_256[k]
                       + kg_256[k];

            t_359[k] = ab_x[k] * ig_257[k]
                       + kg_257[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, ig_258, ig_259, ig_260, \
                         ig_261, ig_262, kg_258, kg_259, kg_260, kg_261, \
                         kg_262 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = ab_x[k] * ig_258[k]
                       + kg_258[k];

            t_361[k] = ab_x[k] * ig_259[k]
                       + kg_259[k];

            t_362[k] = ab_x[k] * ig_260[k]
                       + kg_260[k];

            t_363[k] = ab_x[k] * ig_261[k]
                       + kg_261[k];

            t_364[k] = ab_x[k] * ig_262[k]
                       + kg_262[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, ig_263, ig_264, ig_265, \
                         ig_266, ig_267, kg_263, kg_264, kg_265, kg_266, \
                         kg_267 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = ab_x[k] * ig_263[k]
                       + kg_263[k];

            t_366[k] = ab_x[k] * ig_264[k]
                       + kg_264[k];

            t_367[k] = ab_x[k] * ig_265[k]
                       + kg_265[k];

            t_368[k] = ab_x[k] * ig_266[k]
                       + kg_266[k];

            t_369[k] = ab_x[k] * ig_267[k]
                       + kg_267[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, ab_x, ab_y, ig_265, ig_266, ig_268, \
                         ig_269, kg_268, kg_269, kg_355, kg_356 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = ab_x[k] * ig_268[k]
                       + kg_268[k];

            t_371[k] = ab_x[k] * ig_269[k]
                       + kg_269[k];

            t_372[k] = ab_y[k] * ig_265[k]
                       + kg_355[k];

            t_373[k] = ab_y[k] * ig_266[k]
                       + kg_356[k];
        }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, ab_y, ab_z, ig_267, ig_268, ig_269, \
                         kg_357, kg_358, kg_359, kg_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_374[k] = ab_y[k] * ig_267[k]
                       + kg_357[k];

            t_375[k] = ab_y[k] * ig_268[k]
                       + kg_358[k];

            t_376[k] = ab_y[k] * ig_269[k]
                       + kg_359[k];

            t_377[k] = ab_z[k] * ig_269[k]
                       + kg_374[k];
        }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, t_382, ab_x, ig_270, ig_271, ig_272, \
                         ig_273, ig_274, kg_270, kg_271, kg_272, kg_273, \
                         kg_274 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_378[k] = ab_x[k] * ig_270[k]
                       + kg_270[k];

            t_379[k] = ab_x[k] * ig_271[k]
                       + kg_271[k];

            t_380[k] = ab_x[k] * ig_272[k]
                       + kg_272[k];

            t_381[k] = ab_x[k] * ig_273[k]
                       + kg_273[k];

            t_382[k] = ab_x[k] * ig_274[k]
                       + kg_274[k];
        }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, ab_x, ig_275, ig_276, ig_277, \
                         ig_278, ig_279, kg_275, kg_276, kg_277, kg_278, \
                         kg_279 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_383[k] = ab_x[k] * ig_275[k]
                       + kg_275[k];

            t_384[k] = ab_x[k] * ig_276[k]
                       + kg_276[k];

            t_385[k] = ab_x[k] * ig_277[k]
                       + kg_277[k];

            t_386[k] = ab_x[k] * ig_278[k]
                       + kg_278[k];

            t_387[k] = ab_x[k] * ig_279[k]
                       + kg_279[k];
        }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, t_392, ab_x, ig_280, ig_281, ig_282, \
                         ig_283, ig_284, kg_280, kg_281, kg_282, kg_283, \
                         kg_284 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_388[k] = ab_x[k] * ig_280[k]
                       + kg_280[k];

            t_389[k] = ab_x[k] * ig_281[k]
                       + kg_281[k];

            t_390[k] = ab_x[k] * ig_282[k]
                       + kg_282[k];

            t_391[k] = ab_x[k] * ig_283[k]
                       + kg_283[k];

            t_392[k] = ab_x[k] * ig_284[k]
                       + kg_284[k];
        }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, t_397, ab_y, ig_280, ig_281, ig_282, \
                         ig_283, ig_284, kg_370, kg_371, kg_372, kg_373, \
                         kg_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_393[k] = ab_y[k] * ig_280[k]
                       + kg_370[k];

            t_394[k] = ab_y[k] * ig_281[k]
                       + kg_371[k];

            t_395[k] = ab_y[k] * ig_282[k]
                       + kg_372[k];

            t_396[k] = ab_y[k] * ig_283[k]
                       + kg_373[k];

            t_397[k] = ab_y[k] * ig_284[k]
                       + kg_374[k];
        }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, ab_x, ab_z, ig_284, ig_285, ig_286, \
                         ig_287, kg_285, kg_286, kg_287, kg_389 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_398[k] = ab_z[k] * ig_284[k]
                       + kg_389[k];

            t_399[k] = ab_x[k] * ig_285[k]
                       + kg_285[k];

            t_400[k] = ab_x[k] * ig_286[k]
                       + kg_286[k];

            t_401[k] = ab_x[k] * ig_287[k]
                       + kg_287[k];
        }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, ab_x, ig_288, ig_289, ig_290, \
                         ig_291, ig_292, kg_288, kg_289, kg_290, kg_291, \
                         kg_292 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_402[k] = ab_x[k] * ig_288[k]
                       + kg_288[k];

            t_403[k] = ab_x[k] * ig_289[k]
                       + kg_289[k];

            t_404[k] = ab_x[k] * ig_290[k]
                       + kg_290[k];

            t_405[k] = ab_x[k] * ig_291[k]
                       + kg_291[k];

            t_406[k] = ab_x[k] * ig_292[k]
                       + kg_292[k];
        }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, ab_x, ig_293, ig_294, ig_295, \
                         ig_296, ig_297, kg_293, kg_294, kg_295, kg_296, \
                         kg_297 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_407[k] = ab_x[k] * ig_293[k]
                       + kg_293[k];

            t_408[k] = ab_x[k] * ig_294[k]
                       + kg_294[k];

            t_409[k] = ab_x[k] * ig_295[k]
                       + kg_295[k];

            t_410[k] = ab_x[k] * ig_296[k]
                       + kg_296[k];

            t_411[k] = ab_x[k] * ig_297[k]
                       + kg_297[k];
        }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, ab_x, ab_y, ig_295, ig_296, ig_298, \
                         ig_299, kg_298, kg_299, kg_385, kg_386 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_412[k] = ab_x[k] * ig_298[k]
                       + kg_298[k];

            t_413[k] = ab_x[k] * ig_299[k]
                       + kg_299[k];

            t_414[k] = ab_y[k] * ig_295[k]
                       + kg_385[k];

            t_415[k] = ab_y[k] * ig_296[k]
                       + kg_386[k];
        }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, ab_y, ab_z, ig_297, ig_298, ig_299, \
                         kg_387, kg_388, kg_389, kg_404 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_416[k] = ab_y[k] * ig_297[k]
                       + kg_387[k];

            t_417[k] = ab_y[k] * ig_298[k]
                       + kg_388[k];

            t_418[k] = ab_y[k] * ig_299[k]
                       + kg_389[k];

            t_419[k] = ab_z[k] * ig_299[k]
                       + kg_404[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_x, ig_300, ig_301, ig_302, \
                         ig_303, ig_304, kg_300, kg_301, kg_302, kg_303, \
                         kg_304 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = ab_x[k] * ig_300[k]
                       + kg_300[k];

            t_421[k] = ab_x[k] * ig_301[k]
                       + kg_301[k];

            t_422[k] = ab_x[k] * ig_302[k]
                       + kg_302[k];

            t_423[k] = ab_x[k] * ig_303[k]
                       + kg_303[k];

            t_424[k] = ab_x[k] * ig_304[k]
                       + kg_304[k];
        }
    }
}

static auto
compute_hrr_ih_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t ig, const size_t kg, const size_t ncomps,
                      const size_t nmax) -> void
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

        const auto *ig_305 = buffer.data(ig + 305 * ncomps + c);
        const auto *ig_306 = buffer.data(ig + 306 * ncomps + c);
        const auto *ig_307 = buffer.data(ig + 307 * ncomps + c);
        const auto *ig_308 = buffer.data(ig + 308 * ncomps + c);
        const auto *ig_309 = buffer.data(ig + 309 * ncomps + c);
        const auto *ig_310 = buffer.data(ig + 310 * ncomps + c);
        const auto *ig_311 = buffer.data(ig + 311 * ncomps + c);
        const auto *ig_312 = buffer.data(ig + 312 * ncomps + c);
        const auto *ig_313 = buffer.data(ig + 313 * ncomps + c);
        const auto *ig_314 = buffer.data(ig + 314 * ncomps + c);
        const auto *ig_315 = buffer.data(ig + 315 * ncomps + c);
        const auto *ig_316 = buffer.data(ig + 316 * ncomps + c);
        const auto *ig_317 = buffer.data(ig + 317 * ncomps + c);
        const auto *ig_318 = buffer.data(ig + 318 * ncomps + c);
        const auto *ig_319 = buffer.data(ig + 319 * ncomps + c);
        const auto *ig_320 = buffer.data(ig + 320 * ncomps + c);
        const auto *ig_321 = buffer.data(ig + 321 * ncomps + c);
        const auto *ig_322 = buffer.data(ig + 322 * ncomps + c);
        const auto *ig_323 = buffer.data(ig + 323 * ncomps + c);
        const auto *ig_324 = buffer.data(ig + 324 * ncomps + c);
        const auto *ig_325 = buffer.data(ig + 325 * ncomps + c);
        const auto *ig_326 = buffer.data(ig + 326 * ncomps + c);
        const auto *ig_327 = buffer.data(ig + 327 * ncomps + c);
        const auto *ig_328 = buffer.data(ig + 328 * ncomps + c);
        const auto *ig_329 = buffer.data(ig + 329 * ncomps + c);
        const auto *ig_330 = buffer.data(ig + 330 * ncomps + c);
        const auto *ig_331 = buffer.data(ig + 331 * ncomps + c);
        const auto *ig_332 = buffer.data(ig + 332 * ncomps + c);
        const auto *ig_333 = buffer.data(ig + 333 * ncomps + c);
        const auto *ig_334 = buffer.data(ig + 334 * ncomps + c);
        const auto *ig_335 = buffer.data(ig + 335 * ncomps + c);
        const auto *ig_336 = buffer.data(ig + 336 * ncomps + c);
        const auto *ig_337 = buffer.data(ig + 337 * ncomps + c);
        const auto *ig_338 = buffer.data(ig + 338 * ncomps + c);
        const auto *ig_339 = buffer.data(ig + 339 * ncomps + c);
        const auto *ig_340 = buffer.data(ig + 340 * ncomps + c);
        const auto *ig_341 = buffer.data(ig + 341 * ncomps + c);
        const auto *ig_342 = buffer.data(ig + 342 * ncomps + c);
        const auto *ig_343 = buffer.data(ig + 343 * ncomps + c);
        const auto *ig_344 = buffer.data(ig + 344 * ncomps + c);
        const auto *ig_345 = buffer.data(ig + 345 * ncomps + c);
        const auto *ig_346 = buffer.data(ig + 346 * ncomps + c);
        const auto *ig_347 = buffer.data(ig + 347 * ncomps + c);
        const auto *ig_348 = buffer.data(ig + 348 * ncomps + c);
        const auto *ig_349 = buffer.data(ig + 349 * ncomps + c);
        const auto *ig_350 = buffer.data(ig + 350 * ncomps + c);
        const auto *ig_351 = buffer.data(ig + 351 * ncomps + c);
        const auto *ig_352 = buffer.data(ig + 352 * ncomps + c);
        const auto *ig_353 = buffer.data(ig + 353 * ncomps + c);
        const auto *ig_354 = buffer.data(ig + 354 * ncomps + c);
        const auto *ig_355 = buffer.data(ig + 355 * ncomps + c);
        const auto *ig_356 = buffer.data(ig + 356 * ncomps + c);
        const auto *ig_357 = buffer.data(ig + 357 * ncomps + c);
        const auto *ig_358 = buffer.data(ig + 358 * ncomps + c);
        const auto *ig_359 = buffer.data(ig + 359 * ncomps + c);
        const auto *ig_360 = buffer.data(ig + 360 * ncomps + c);
        const auto *ig_361 = buffer.data(ig + 361 * ncomps + c);
        const auto *ig_362 = buffer.data(ig + 362 * ncomps + c);
        const auto *ig_363 = buffer.data(ig + 363 * ncomps + c);
        const auto *ig_364 = buffer.data(ig + 364 * ncomps + c);
        const auto *ig_365 = buffer.data(ig + 365 * ncomps + c);
        const auto *ig_366 = buffer.data(ig + 366 * ncomps + c);
        const auto *ig_367 = buffer.data(ig + 367 * ncomps + c);
        const auto *ig_368 = buffer.data(ig + 368 * ncomps + c);
        const auto *ig_369 = buffer.data(ig + 369 * ncomps + c);
        const auto *ig_370 = buffer.data(ig + 370 * ncomps + c);
        const auto *ig_371 = buffer.data(ig + 371 * ncomps + c);
        const auto *ig_372 = buffer.data(ig + 372 * ncomps + c);
        const auto *ig_373 = buffer.data(ig + 373 * ncomps + c);
        const auto *ig_374 = buffer.data(ig + 374 * ncomps + c);
        const auto *ig_375 = buffer.data(ig + 375 * ncomps + c);
        const auto *ig_376 = buffer.data(ig + 376 * ncomps + c);
        const auto *ig_377 = buffer.data(ig + 377 * ncomps + c);
        const auto *ig_378 = buffer.data(ig + 378 * ncomps + c);
        const auto *ig_379 = buffer.data(ig + 379 * ncomps + c);
        const auto *ig_380 = buffer.data(ig + 380 * ncomps + c);
        const auto *ig_381 = buffer.data(ig + 381 * ncomps + c);
        const auto *ig_382 = buffer.data(ig + 382 * ncomps + c);
        const auto *ig_383 = buffer.data(ig + 383 * ncomps + c);
        const auto *ig_384 = buffer.data(ig + 384 * ncomps + c);
        const auto *ig_385 = buffer.data(ig + 385 * ncomps + c);
        const auto *ig_386 = buffer.data(ig + 386 * ncomps + c);
        const auto *ig_387 = buffer.data(ig + 387 * ncomps + c);
        const auto *ig_388 = buffer.data(ig + 388 * ncomps + c);
        const auto *ig_389 = buffer.data(ig + 389 * ncomps + c);
        const auto *ig_390 = buffer.data(ig + 390 * ncomps + c);
        const auto *ig_391 = buffer.data(ig + 391 * ncomps + c);
        const auto *ig_392 = buffer.data(ig + 392 * ncomps + c);
        const auto *ig_393 = buffer.data(ig + 393 * ncomps + c);
        const auto *ig_394 = buffer.data(ig + 394 * ncomps + c);
        const auto *ig_395 = buffer.data(ig + 395 * ncomps + c);
        const auto *ig_396 = buffer.data(ig + 396 * ncomps + c);
        const auto *ig_397 = buffer.data(ig + 397 * ncomps + c);
        const auto *ig_398 = buffer.data(ig + 398 * ncomps + c);
        const auto *ig_399 = buffer.data(ig + 399 * ncomps + c);
        const auto *ig_400 = buffer.data(ig + 400 * ncomps + c);
        const auto *ig_401 = buffer.data(ig + 401 * ncomps + c);
        const auto *ig_402 = buffer.data(ig + 402 * ncomps + c);
        const auto *ig_403 = buffer.data(ig + 403 * ncomps + c);
        const auto *ig_404 = buffer.data(ig + 404 * ncomps + c);

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
        const auto *kg_419 = buffer.data(kg + 419 * ncomps + c);
        const auto *kg_430 = buffer.data(kg + 430 * ncomps + c);
        const auto *kg_431 = buffer.data(kg + 431 * ncomps + c);
        const auto *kg_432 = buffer.data(kg + 432 * ncomps + c);
        const auto *kg_433 = buffer.data(kg + 433 * ncomps + c);
        const auto *kg_434 = buffer.data(kg + 434 * ncomps + c);
        const auto *kg_445 = buffer.data(kg + 445 * ncomps + c);
        const auto *kg_446 = buffer.data(kg + 446 * ncomps + c);
        const auto *kg_447 = buffer.data(kg + 447 * ncomps + c);
        const auto *kg_448 = buffer.data(kg + 448 * ncomps + c);
        const auto *kg_449 = buffer.data(kg + 449 * ncomps + c);
        const auto *kg_460 = buffer.data(kg + 460 * ncomps + c);
        const auto *kg_461 = buffer.data(kg + 461 * ncomps + c);
        const auto *kg_462 = buffer.data(kg + 462 * ncomps + c);
        const auto *kg_463 = buffer.data(kg + 463 * ncomps + c);
        const auto *kg_464 = buffer.data(kg + 464 * ncomps + c);
        const auto *kg_475 = buffer.data(kg + 475 * ncomps + c);
        const auto *kg_476 = buffer.data(kg + 476 * ncomps + c);
        const auto *kg_477 = buffer.data(kg + 477 * ncomps + c);
        const auto *kg_478 = buffer.data(kg + 478 * ncomps + c);
        const auto *kg_479 = buffer.data(kg + 479 * ncomps + c);
        const auto *kg_490 = buffer.data(kg + 490 * ncomps + c);
        const auto *kg_491 = buffer.data(kg + 491 * ncomps + c);
        const auto *kg_492 = buffer.data(kg + 492 * ncomps + c);
        const auto *kg_493 = buffer.data(kg + 493 * ncomps + c);
        const auto *kg_494 = buffer.data(kg + 494 * ncomps + c);
        const auto *kg_505 = buffer.data(kg + 505 * ncomps + c);
        const auto *kg_506 = buffer.data(kg + 506 * ncomps + c);
        const auto *kg_507 = buffer.data(kg + 507 * ncomps + c);
        const auto *kg_508 = buffer.data(kg + 508 * ncomps + c);
        const auto *kg_509 = buffer.data(kg + 509 * ncomps + c);

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_x, ig_305, ig_306, ig_307, \
                         ig_308, ig_309, kg_305, kg_306, kg_307, kg_308, \
                         kg_309 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = ab_x[k] * ig_305[k]
                       + kg_305[k];

            t_426[k] = ab_x[k] * ig_306[k]
                       + kg_306[k];

            t_427[k] = ab_x[k] * ig_307[k]
                       + kg_307[k];

            t_428[k] = ab_x[k] * ig_308[k]
                       + kg_308[k];

            t_429[k] = ab_x[k] * ig_309[k]
                       + kg_309[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_x, ig_310, ig_311, ig_312, \
                         ig_313, ig_314, kg_310, kg_311, kg_312, kg_313, \
                         kg_314 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = ab_x[k] * ig_310[k]
                       + kg_310[k];

            t_431[k] = ab_x[k] * ig_311[k]
                       + kg_311[k];

            t_432[k] = ab_x[k] * ig_312[k]
                       + kg_312[k];

            t_433[k] = ab_x[k] * ig_313[k]
                       + kg_313[k];

            t_434[k] = ab_x[k] * ig_314[k]
                       + kg_314[k];
        }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_y, ig_310, ig_311, ig_312, \
                         ig_313, ig_314, kg_400, kg_401, kg_402, kg_403, \
                         kg_404 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = ab_y[k] * ig_310[k]
                       + kg_400[k];

            t_436[k] = ab_y[k] * ig_311[k]
                       + kg_401[k];

            t_437[k] = ab_y[k] * ig_312[k]
                       + kg_402[k];

            t_438[k] = ab_y[k] * ig_313[k]
                       + kg_403[k];

            t_439[k] = ab_y[k] * ig_314[k]
                       + kg_404[k];
        }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, ab_x, ab_z, ig_314, ig_315, ig_316, \
                         ig_317, kg_315, kg_316, kg_317, kg_419 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = ab_z[k] * ig_314[k]
                       + kg_419[k];

            t_441[k] = ab_x[k] * ig_315[k]
                       + kg_315[k];

            t_442[k] = ab_x[k] * ig_316[k]
                       + kg_316[k];

            t_443[k] = ab_x[k] * ig_317[k]
                       + kg_317[k];
        }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, t_448, ab_x, ig_318, ig_319, ig_320, \
                         ig_321, ig_322, kg_318, kg_319, kg_320, kg_321, \
                         kg_322 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_444[k] = ab_x[k] * ig_318[k]
                       + kg_318[k];

            t_445[k] = ab_x[k] * ig_319[k]
                       + kg_319[k];

            t_446[k] = ab_x[k] * ig_320[k]
                       + kg_320[k];

            t_447[k] = ab_x[k] * ig_321[k]
                       + kg_321[k];

            t_448[k] = ab_x[k] * ig_322[k]
                       + kg_322[k];
        }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, t_453, ab_x, ig_323, ig_324, ig_325, \
                         ig_326, ig_327, kg_323, kg_324, kg_325, kg_326, \
                         kg_327 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_449[k] = ab_x[k] * ig_323[k]
                       + kg_323[k];

            t_450[k] = ab_x[k] * ig_324[k]
                       + kg_324[k];

            t_451[k] = ab_x[k] * ig_325[k]
                       + kg_325[k];

            t_452[k] = ab_x[k] * ig_326[k]
                       + kg_326[k];

            t_453[k] = ab_x[k] * ig_327[k]
                       + kg_327[k];
        }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, ab_x, ab_y, ig_325, ig_326, ig_328, \
                         ig_329, kg_328, kg_329, kg_430, kg_431 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_454[k] = ab_x[k] * ig_328[k]
                       + kg_328[k];

            t_455[k] = ab_x[k] * ig_329[k]
                       + kg_329[k];

            t_456[k] = ab_y[k] * ig_325[k]
                       + kg_430[k];

            t_457[k] = ab_y[k] * ig_326[k]
                       + kg_431[k];
        }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, ab_y, ab_z, ig_327, ig_328, ig_329, \
                         kg_432, kg_433, kg_434, kg_449 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_458[k] = ab_y[k] * ig_327[k]
                       + kg_432[k];

            t_459[k] = ab_y[k] * ig_328[k]
                       + kg_433[k];

            t_460[k] = ab_y[k] * ig_329[k]
                       + kg_434[k];

            t_461[k] = ab_z[k] * ig_329[k]
                       + kg_449[k];
        }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, t_466, ab_x, ig_330, ig_331, ig_332, \
                         ig_333, ig_334, kg_330, kg_331, kg_332, kg_333, \
                         kg_334 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_462[k] = ab_x[k] * ig_330[k]
                       + kg_330[k];

            t_463[k] = ab_x[k] * ig_331[k]
                       + kg_331[k];

            t_464[k] = ab_x[k] * ig_332[k]
                       + kg_332[k];

            t_465[k] = ab_x[k] * ig_333[k]
                       + kg_333[k];

            t_466[k] = ab_x[k] * ig_334[k]
                       + kg_334[k];
        }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, t_471, ab_x, ig_335, ig_336, ig_337, \
                         ig_338, ig_339, kg_335, kg_336, kg_337, kg_338, \
                         kg_339 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_467[k] = ab_x[k] * ig_335[k]
                       + kg_335[k];

            t_468[k] = ab_x[k] * ig_336[k]
                       + kg_336[k];

            t_469[k] = ab_x[k] * ig_337[k]
                       + kg_337[k];

            t_470[k] = ab_x[k] * ig_338[k]
                       + kg_338[k];

            t_471[k] = ab_x[k] * ig_339[k]
                       + kg_339[k];
        }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, t_476, ab_x, ig_340, ig_341, ig_342, \
                         ig_343, ig_344, kg_340, kg_341, kg_342, kg_343, \
                         kg_344 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_472[k] = ab_x[k] * ig_340[k]
                       + kg_340[k];

            t_473[k] = ab_x[k] * ig_341[k]
                       + kg_341[k];

            t_474[k] = ab_x[k] * ig_342[k]
                       + kg_342[k];

            t_475[k] = ab_x[k] * ig_343[k]
                       + kg_343[k];

            t_476[k] = ab_x[k] * ig_344[k]
                       + kg_344[k];
        }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, ab_y, ig_340, ig_341, ig_342, \
                         ig_343, ig_344, kg_445, kg_446, kg_447, kg_448, \
                         kg_449 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_477[k] = ab_y[k] * ig_340[k]
                       + kg_445[k];

            t_478[k] = ab_y[k] * ig_341[k]
                       + kg_446[k];

            t_479[k] = ab_y[k] * ig_342[k]
                       + kg_447[k];

            t_480[k] = ab_y[k] * ig_343[k]
                       + kg_448[k];

            t_481[k] = ab_y[k] * ig_344[k]
                       + kg_449[k];
        }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, ab_x, ab_z, ig_344, ig_345, ig_346, \
                         ig_347, kg_345, kg_346, kg_347, kg_464 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_482[k] = ab_z[k] * ig_344[k]
                       + kg_464[k];

            t_483[k] = ab_x[k] * ig_345[k]
                       + kg_345[k];

            t_484[k] = ab_x[k] * ig_346[k]
                       + kg_346[k];

            t_485[k] = ab_x[k] * ig_347[k]
                       + kg_347[k];
        }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, t_490, ab_x, ig_348, ig_349, ig_350, \
                         ig_351, ig_352, kg_348, kg_349, kg_350, kg_351, \
                         kg_352 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_486[k] = ab_x[k] * ig_348[k]
                       + kg_348[k];

            t_487[k] = ab_x[k] * ig_349[k]
                       + kg_349[k];

            t_488[k] = ab_x[k] * ig_350[k]
                       + kg_350[k];

            t_489[k] = ab_x[k] * ig_351[k]
                       + kg_351[k];

            t_490[k] = ab_x[k] * ig_352[k]
                       + kg_352[k];
        }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, ab_x, ig_353, ig_354, ig_355, \
                         ig_356, ig_357, kg_353, kg_354, kg_355, kg_356, \
                         kg_357 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_491[k] = ab_x[k] * ig_353[k]
                       + kg_353[k];

            t_492[k] = ab_x[k] * ig_354[k]
                       + kg_354[k];

            t_493[k] = ab_x[k] * ig_355[k]
                       + kg_355[k];

            t_494[k] = ab_x[k] * ig_356[k]
                       + kg_356[k];

            t_495[k] = ab_x[k] * ig_357[k]
                       + kg_357[k];
        }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, ab_x, ab_y, ig_355, ig_356, ig_358, \
                         ig_359, kg_358, kg_359, kg_460, kg_461 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_496[k] = ab_x[k] * ig_358[k]
                       + kg_358[k];

            t_497[k] = ab_x[k] * ig_359[k]
                       + kg_359[k];

            t_498[k] = ab_y[k] * ig_355[k]
                       + kg_460[k];

            t_499[k] = ab_y[k] * ig_356[k]
                       + kg_461[k];
        }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, ab_y, ab_z, ig_357, ig_358, ig_359, \
                         kg_462, kg_463, kg_464, kg_479 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_500[k] = ab_y[k] * ig_357[k]
                       + kg_462[k];

            t_501[k] = ab_y[k] * ig_358[k]
                       + kg_463[k];

            t_502[k] = ab_y[k] * ig_359[k]
                       + kg_464[k];

            t_503[k] = ab_z[k] * ig_359[k]
                       + kg_479[k];
        }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, ab_x, ig_360, ig_361, ig_362, \
                         ig_363, ig_364, kg_360, kg_361, kg_362, kg_363, \
                         kg_364 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_504[k] = ab_x[k] * ig_360[k]
                       + kg_360[k];

            t_505[k] = ab_x[k] * ig_361[k]
                       + kg_361[k];

            t_506[k] = ab_x[k] * ig_362[k]
                       + kg_362[k];

            t_507[k] = ab_x[k] * ig_363[k]
                       + kg_363[k];

            t_508[k] = ab_x[k] * ig_364[k]
                       + kg_364[k];
        }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, t_513, ab_x, ig_365, ig_366, ig_367, \
                         ig_368, ig_369, kg_365, kg_366, kg_367, kg_368, \
                         kg_369 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_509[k] = ab_x[k] * ig_365[k]
                       + kg_365[k];

            t_510[k] = ab_x[k] * ig_366[k]
                       + kg_366[k];

            t_511[k] = ab_x[k] * ig_367[k]
                       + kg_367[k];

            t_512[k] = ab_x[k] * ig_368[k]
                       + kg_368[k];

            t_513[k] = ab_x[k] * ig_369[k]
                       + kg_369[k];
        }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, t_518, ab_x, ig_370, ig_371, ig_372, \
                         ig_373, ig_374, kg_370, kg_371, kg_372, kg_373, \
                         kg_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_514[k] = ab_x[k] * ig_370[k]
                       + kg_370[k];

            t_515[k] = ab_x[k] * ig_371[k]
                       + kg_371[k];

            t_516[k] = ab_x[k] * ig_372[k]
                       + kg_372[k];

            t_517[k] = ab_x[k] * ig_373[k]
                       + kg_373[k];

            t_518[k] = ab_x[k] * ig_374[k]
                       + kg_374[k];
        }

#pragma omp simd aligned(t_519, t_520, t_521, t_522, t_523, ab_y, ig_370, ig_371, ig_372, \
                         ig_373, ig_374, kg_475, kg_476, kg_477, kg_478, \
                         kg_479 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_519[k] = ab_y[k] * ig_370[k]
                       + kg_475[k];

            t_520[k] = ab_y[k] * ig_371[k]
                       + kg_476[k];

            t_521[k] = ab_y[k] * ig_372[k]
                       + kg_477[k];

            t_522[k] = ab_y[k] * ig_373[k]
                       + kg_478[k];

            t_523[k] = ab_y[k] * ig_374[k]
                       + kg_479[k];
        }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, ab_x, ab_z, ig_374, ig_375, ig_376, \
                         ig_377, kg_375, kg_376, kg_377, kg_494 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_524[k] = ab_z[k] * ig_374[k]
                       + kg_494[k];

            t_525[k] = ab_x[k] * ig_375[k]
                       + kg_375[k];

            t_526[k] = ab_x[k] * ig_376[k]
                       + kg_376[k];

            t_527[k] = ab_x[k] * ig_377[k]
                       + kg_377[k];
        }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, t_532, ab_x, ig_378, ig_379, ig_380, \
                         ig_381, ig_382, kg_378, kg_379, kg_380, kg_381, \
                         kg_382 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_528[k] = ab_x[k] * ig_378[k]
                       + kg_378[k];

            t_529[k] = ab_x[k] * ig_379[k]
                       + kg_379[k];

            t_530[k] = ab_x[k] * ig_380[k]
                       + kg_380[k];

            t_531[k] = ab_x[k] * ig_381[k]
                       + kg_381[k];

            t_532[k] = ab_x[k] * ig_382[k]
                       + kg_382[k];
        }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, t_537, ab_x, ig_383, ig_384, ig_385, \
                         ig_386, ig_387, kg_383, kg_384, kg_385, kg_386, \
                         kg_387 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_533[k] = ab_x[k] * ig_383[k]
                       + kg_383[k];

            t_534[k] = ab_x[k] * ig_384[k]
                       + kg_384[k];

            t_535[k] = ab_x[k] * ig_385[k]
                       + kg_385[k];

            t_536[k] = ab_x[k] * ig_386[k]
                       + kg_386[k];

            t_537[k] = ab_x[k] * ig_387[k]
                       + kg_387[k];
        }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, ab_x, ab_y, ig_385, ig_386, ig_388, \
                         ig_389, kg_388, kg_389, kg_490, kg_491 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_538[k] = ab_x[k] * ig_388[k]
                       + kg_388[k];

            t_539[k] = ab_x[k] * ig_389[k]
                       + kg_389[k];

            t_540[k] = ab_y[k] * ig_385[k]
                       + kg_490[k];

            t_541[k] = ab_y[k] * ig_386[k]
                       + kg_491[k];
        }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, ab_y, ab_z, ig_387, ig_388, ig_389, \
                         kg_492, kg_493, kg_494, kg_509 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_542[k] = ab_y[k] * ig_387[k]
                       + kg_492[k];

            t_543[k] = ab_y[k] * ig_388[k]
                       + kg_493[k];

            t_544[k] = ab_y[k] * ig_389[k]
                       + kg_494[k];

            t_545[k] = ab_z[k] * ig_389[k]
                       + kg_509[k];
        }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, t_550, ab_x, ig_390, ig_391, ig_392, \
                         ig_393, ig_394, kg_390, kg_391, kg_392, kg_393, \
                         kg_394 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_546[k] = ab_x[k] * ig_390[k]
                       + kg_390[k];

            t_547[k] = ab_x[k] * ig_391[k]
                       + kg_391[k];

            t_548[k] = ab_x[k] * ig_392[k]
                       + kg_392[k];

            t_549[k] = ab_x[k] * ig_393[k]
                       + kg_393[k];

            t_550[k] = ab_x[k] * ig_394[k]
                       + kg_394[k];
        }

#pragma omp simd aligned(t_551, t_552, t_553, t_554, t_555, ab_x, ig_395, ig_396, ig_397, \
                         ig_398, ig_399, kg_395, kg_396, kg_397, kg_398, \
                         kg_399 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_551[k] = ab_x[k] * ig_395[k]
                       + kg_395[k];

            t_552[k] = ab_x[k] * ig_396[k]
                       + kg_396[k];

            t_553[k] = ab_x[k] * ig_397[k]
                       + kg_397[k];

            t_554[k] = ab_x[k] * ig_398[k]
                       + kg_398[k];

            t_555[k] = ab_x[k] * ig_399[k]
                       + kg_399[k];
        }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, ab_x, ig_400, ig_401, ig_402, \
                         ig_403, ig_404, kg_400, kg_401, kg_402, kg_403, \
                         kg_404 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_556[k] = ab_x[k] * ig_400[k]
                       + kg_400[k];

            t_557[k] = ab_x[k] * ig_401[k]
                       + kg_401[k];

            t_558[k] = ab_x[k] * ig_402[k]
                       + kg_402[k];

            t_559[k] = ab_x[k] * ig_403[k]
                       + kg_403[k];

            t_560[k] = ab_x[k] * ig_404[k]
                       + kg_404[k];
        }

#pragma omp simd aligned(t_561, t_562, t_563, t_564, t_565, ab_y, ig_400, ig_401, ig_402, \
                         ig_403, ig_404, kg_505, kg_506, kg_507, kg_508, \
                         kg_509 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_561[k] = ab_y[k] * ig_400[k]
                       + kg_505[k];

            t_562[k] = ab_y[k] * ig_401[k]
                       + kg_506[k];

            t_563[k] = ab_y[k] * ig_402[k]
                       + kg_507[k];

            t_564[k] = ab_y[k] * ig_403[k]
                       + kg_508[k];

            t_565[k] = ab_y[k] * ig_404[k]
                       + kg_509[k];
        }
    }
}

static auto
compute_hrr_ih_piece4(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t ig, const size_t kg, const size_t ncomps,
                      const size_t nmax) -> void
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *ig_404 = buffer.data(ig + 404 * ncomps + c);
        const auto *ig_405 = buffer.data(ig + 405 * ncomps + c);
        const auto *ig_406 = buffer.data(ig + 406 * ncomps + c);
        const auto *ig_407 = buffer.data(ig + 407 * ncomps + c);
        const auto *ig_408 = buffer.data(ig + 408 * ncomps + c);
        const auto *ig_409 = buffer.data(ig + 409 * ncomps + c);
        const auto *ig_410 = buffer.data(ig + 410 * ncomps + c);
        const auto *ig_411 = buffer.data(ig + 411 * ncomps + c);
        const auto *ig_412 = buffer.data(ig + 412 * ncomps + c);
        const auto *ig_413 = buffer.data(ig + 413 * ncomps + c);
        const auto *ig_414 = buffer.data(ig + 414 * ncomps + c);
        const auto *ig_415 = buffer.data(ig + 415 * ncomps + c);
        const auto *ig_416 = buffer.data(ig + 416 * ncomps + c);
        const auto *ig_417 = buffer.data(ig + 417 * ncomps + c);
        const auto *ig_418 = buffer.data(ig + 418 * ncomps + c);
        const auto *ig_419 = buffer.data(ig + 419 * ncomps + c);

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
        const auto *kg_520 = buffer.data(kg + 520 * ncomps + c);
        const auto *kg_521 = buffer.data(kg + 521 * ncomps + c);
        const auto *kg_522 = buffer.data(kg + 522 * ncomps + c);
        const auto *kg_523 = buffer.data(kg + 523 * ncomps + c);
        const auto *kg_524 = buffer.data(kg + 524 * ncomps + c);
        const auto *kg_539 = buffer.data(kg + 539 * ncomps + c);

#pragma omp simd aligned(t_566, t_567, t_568, t_569, ab_x, ab_z, ig_404, ig_405, ig_406, \
                         ig_407, kg_405, kg_406, kg_407, kg_524 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_566[k] = ab_z[k] * ig_404[k]
                       + kg_524[k];

            t_567[k] = ab_x[k] * ig_405[k]
                       + kg_405[k];

            t_568[k] = ab_x[k] * ig_406[k]
                       + kg_406[k];

            t_569[k] = ab_x[k] * ig_407[k]
                       + kg_407[k];
        }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, ab_x, ig_408, ig_409, ig_410, \
                         ig_411, ig_412, kg_408, kg_409, kg_410, kg_411, \
                         kg_412 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_570[k] = ab_x[k] * ig_408[k]
                       + kg_408[k];

            t_571[k] = ab_x[k] * ig_409[k]
                       + kg_409[k];

            t_572[k] = ab_x[k] * ig_410[k]
                       + kg_410[k];

            t_573[k] = ab_x[k] * ig_411[k]
                       + kg_411[k];

            t_574[k] = ab_x[k] * ig_412[k]
                       + kg_412[k];
        }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, ab_x, ig_413, ig_414, ig_415, \
                         ig_416, ig_417, kg_413, kg_414, kg_415, kg_416, \
                         kg_417 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_575[k] = ab_x[k] * ig_413[k]
                       + kg_413[k];

            t_576[k] = ab_x[k] * ig_414[k]
                       + kg_414[k];

            t_577[k] = ab_x[k] * ig_415[k]
                       + kg_415[k];

            t_578[k] = ab_x[k] * ig_416[k]
                       + kg_416[k];

            t_579[k] = ab_x[k] * ig_417[k]
                       + kg_417[k];
        }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, ab_x, ab_y, ig_415, ig_416, ig_418, \
                         ig_419, kg_418, kg_419, kg_520, kg_521 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_580[k] = ab_x[k] * ig_418[k]
                       + kg_418[k];

            t_581[k] = ab_x[k] * ig_419[k]
                       + kg_419[k];

            t_582[k] = ab_y[k] * ig_415[k]
                       + kg_520[k];

            t_583[k] = ab_y[k] * ig_416[k]
                       + kg_521[k];
        }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, ab_y, ab_z, ig_417, ig_418, ig_419, \
                         kg_522, kg_523, kg_524, kg_539 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_584[k] = ab_y[k] * ig_417[k]
                       + kg_522[k];

            t_585[k] = ab_y[k] * ig_418[k]
                       + kg_523[k];

            t_586[k] = ab_y[k] * ig_419[k]
                       + kg_524[k];

            t_587[k] = ab_z[k] * ig_419[k]
                       + kg_539[k];
        }
    }
}

auto
compute_hrr_ih(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t ig, const size_t kg, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_ih_piece0(buffer, coordinates, target, ig, kg, ncomps, nmax);

    compute_hrr_ih_piece1(buffer, coordinates, target, ig, kg, ncomps, nmax);

    compute_hrr_ih_piece2(buffer, coordinates, target, ig, kg, ncomps, nmax);

    compute_hrr_ih_piece3(buffer, coordinates, target, ig, kg, ncomps, nmax);

    compute_hrr_ih_piece4(buffer, coordinates, target, ig, kg, ncomps, nmax);
}

}  // namespace simdtrf
