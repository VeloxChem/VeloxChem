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


#include "SimdTransferHH.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

static auto
compute_hrr_hh_out_of_first_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t hg, const size_t ig,
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

        const auto *hg_0 = buffer.data(hg + 0 * ncomps + c);
        const auto *hg_1 = buffer.data(hg + 1 * ncomps + c);
        const auto *hg_2 = buffer.data(hg + 2 * ncomps + c);
        const auto *hg_3 = buffer.data(hg + 3 * ncomps + c);
        const auto *hg_4 = buffer.data(hg + 4 * ncomps + c);
        const auto *hg_5 = buffer.data(hg + 5 * ncomps + c);
        const auto *hg_6 = buffer.data(hg + 6 * ncomps + c);
        const auto *hg_7 = buffer.data(hg + 7 * ncomps + c);
        const auto *hg_8 = buffer.data(hg + 8 * ncomps + c);
        const auto *hg_9 = buffer.data(hg + 9 * ncomps + c);
        const auto *hg_10 = buffer.data(hg + 10 * ncomps + c);
        const auto *hg_11 = buffer.data(hg + 11 * ncomps + c);
        const auto *hg_12 = buffer.data(hg + 12 * ncomps + c);
        const auto *hg_13 = buffer.data(hg + 13 * ncomps + c);
        const auto *hg_14 = buffer.data(hg + 14 * ncomps + c);
        const auto *hg_15 = buffer.data(hg + 15 * ncomps + c);
        const auto *hg_16 = buffer.data(hg + 16 * ncomps + c);
        const auto *hg_17 = buffer.data(hg + 17 * ncomps + c);
        const auto *hg_18 = buffer.data(hg + 18 * ncomps + c);
        const auto *hg_19 = buffer.data(hg + 19 * ncomps + c);
        const auto *hg_20 = buffer.data(hg + 20 * ncomps + c);
        const auto *hg_21 = buffer.data(hg + 21 * ncomps + c);
        const auto *hg_22 = buffer.data(hg + 22 * ncomps + c);
        const auto *hg_23 = buffer.data(hg + 23 * ncomps + c);
        const auto *hg_24 = buffer.data(hg + 24 * ncomps + c);
        const auto *hg_25 = buffer.data(hg + 25 * ncomps + c);
        const auto *hg_26 = buffer.data(hg + 26 * ncomps + c);
        const auto *hg_27 = buffer.data(hg + 27 * ncomps + c);
        const auto *hg_28 = buffer.data(hg + 28 * ncomps + c);
        const auto *hg_29 = buffer.data(hg + 29 * ncomps + c);
        const auto *hg_30 = buffer.data(hg + 30 * ncomps + c);
        const auto *hg_31 = buffer.data(hg + 31 * ncomps + c);
        const auto *hg_32 = buffer.data(hg + 32 * ncomps + c);
        const auto *hg_33 = buffer.data(hg + 33 * ncomps + c);
        const auto *hg_34 = buffer.data(hg + 34 * ncomps + c);
        const auto *hg_35 = buffer.data(hg + 35 * ncomps + c);
        const auto *hg_36 = buffer.data(hg + 36 * ncomps + c);
        const auto *hg_37 = buffer.data(hg + 37 * ncomps + c);
        const auto *hg_38 = buffer.data(hg + 38 * ncomps + c);
        const auto *hg_39 = buffer.data(hg + 39 * ncomps + c);
        const auto *hg_40 = buffer.data(hg + 40 * ncomps + c);
        const auto *hg_41 = buffer.data(hg + 41 * ncomps + c);
        const auto *hg_42 = buffer.data(hg + 42 * ncomps + c);
        const auto *hg_43 = buffer.data(hg + 43 * ncomps + c);
        const auto *hg_44 = buffer.data(hg + 44 * ncomps + c);
        const auto *hg_45 = buffer.data(hg + 45 * ncomps + c);
        const auto *hg_46 = buffer.data(hg + 46 * ncomps + c);
        const auto *hg_47 = buffer.data(hg + 47 * ncomps + c);
        const auto *hg_48 = buffer.data(hg + 48 * ncomps + c);
        const auto *hg_49 = buffer.data(hg + 49 * ncomps + c);
        const auto *hg_50 = buffer.data(hg + 50 * ncomps + c);
        const auto *hg_51 = buffer.data(hg + 51 * ncomps + c);
        const auto *hg_52 = buffer.data(hg + 52 * ncomps + c);
        const auto *hg_53 = buffer.data(hg + 53 * ncomps + c);
        const auto *hg_54 = buffer.data(hg + 54 * ncomps + c);
        const auto *hg_55 = buffer.data(hg + 55 * ncomps + c);
        const auto *hg_56 = buffer.data(hg + 56 * ncomps + c);
        const auto *hg_57 = buffer.data(hg + 57 * ncomps + c);
        const auto *hg_58 = buffer.data(hg + 58 * ncomps + c);
        const auto *hg_59 = buffer.data(hg + 59 * ncomps + c);
        const auto *hg_60 = buffer.data(hg + 60 * ncomps + c);
        const auto *hg_61 = buffer.data(hg + 61 * ncomps + c);
        const auto *hg_62 = buffer.data(hg + 62 * ncomps + c);
        const auto *hg_63 = buffer.data(hg + 63 * ncomps + c);
        const auto *hg_64 = buffer.data(hg + 64 * ncomps + c);
        const auto *hg_65 = buffer.data(hg + 65 * ncomps + c);
        const auto *hg_66 = buffer.data(hg + 66 * ncomps + c);
        const auto *hg_67 = buffer.data(hg + 67 * ncomps + c);
        const auto *hg_68 = buffer.data(hg + 68 * ncomps + c);
        const auto *hg_69 = buffer.data(hg + 69 * ncomps + c);
        const auto *hg_70 = buffer.data(hg + 70 * ncomps + c);
        const auto *hg_71 = buffer.data(hg + 71 * ncomps + c);
        const auto *hg_72 = buffer.data(hg + 72 * ncomps + c);
        const auto *hg_73 = buffer.data(hg + 73 * ncomps + c);
        const auto *hg_74 = buffer.data(hg + 74 * ncomps + c);
        const auto *hg_75 = buffer.data(hg + 75 * ncomps + c);
        const auto *hg_76 = buffer.data(hg + 76 * ncomps + c);
        const auto *hg_77 = buffer.data(hg + 77 * ncomps + c);
        const auto *hg_78 = buffer.data(hg + 78 * ncomps + c);
        const auto *hg_79 = buffer.data(hg + 79 * ncomps + c);
        const auto *hg_80 = buffer.data(hg + 80 * ncomps + c);
        const auto *hg_81 = buffer.data(hg + 81 * ncomps + c);
        const auto *hg_82 = buffer.data(hg + 82 * ncomps + c);
        const auto *hg_83 = buffer.data(hg + 83 * ncomps + c);
        const auto *hg_84 = buffer.data(hg + 84 * ncomps + c);
        const auto *hg_85 = buffer.data(hg + 85 * ncomps + c);
        const auto *hg_86 = buffer.data(hg + 86 * ncomps + c);
        const auto *hg_87 = buffer.data(hg + 87 * ncomps + c);
        const auto *hg_88 = buffer.data(hg + 88 * ncomps + c);
        const auto *hg_89 = buffer.data(hg + 89 * ncomps + c);
        const auto *hg_90 = buffer.data(hg + 90 * ncomps + c);
        const auto *hg_91 = buffer.data(hg + 91 * ncomps + c);
        const auto *hg_92 = buffer.data(hg + 92 * ncomps + c);
        const auto *hg_93 = buffer.data(hg + 93 * ncomps + c);
        const auto *hg_94 = buffer.data(hg + 94 * ncomps + c);
        const auto *hg_95 = buffer.data(hg + 95 * ncomps + c);
        const auto *hg_96 = buffer.data(hg + 96 * ncomps + c);
        const auto *hg_97 = buffer.data(hg + 97 * ncomps + c);
        const auto *hg_98 = buffer.data(hg + 98 * ncomps + c);
        const auto *hg_99 = buffer.data(hg + 99 * ncomps + c);
        const auto *hg_100 = buffer.data(hg + 100 * ncomps + c);
        const auto *hg_101 = buffer.data(hg + 101 * ncomps + c);
        const auto *hg_102 = buffer.data(hg + 102 * ncomps + c);
        const auto *hg_103 = buffer.data(hg + 103 * ncomps + c);
        const auto *hg_104 = buffer.data(hg + 104 * ncomps + c);

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
        const auto *ig_115 = buffer.data(ig + 115 * ncomps + c);
        const auto *ig_116 = buffer.data(ig + 116 * ncomps + c);
        const auto *ig_117 = buffer.data(ig + 117 * ncomps + c);
        const auto *ig_118 = buffer.data(ig + 118 * ncomps + c);
        const auto *ig_119 = buffer.data(ig + 119 * ncomps + c);
        const auto *ig_130 = buffer.data(ig + 130 * ncomps + c);
        const auto *ig_131 = buffer.data(ig + 131 * ncomps + c);
        const auto *ig_132 = buffer.data(ig + 132 * ncomps + c);
        const auto *ig_133 = buffer.data(ig + 133 * ncomps + c);
        const auto *ig_134 = buffer.data(ig + 134 * ncomps + c);
        const auto *ig_149 = buffer.data(ig + 149 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, hg_0, hg_1, hg_2, hg_3, hg_4, ig_0, \
                         ig_1, ig_2, ig_3, ig_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = ab_x[k] * hg_0[k]
                     + ig_0[k];

            t_1[k] = ab_x[k] * hg_1[k]
                     + ig_1[k];

            t_2[k] = ab_x[k] * hg_2[k]
                     + ig_2[k];

            t_3[k] = ab_x[k] * hg_3[k]
                     + ig_3[k];

            t_4[k] = ab_x[k] * hg_4[k]
                     + ig_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, hg_5, hg_6, hg_7, hg_8, hg_9, ig_5, \
                         ig_6, ig_7, ig_8, ig_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = ab_x[k] * hg_5[k]
                     + ig_5[k];

            t_6[k] = ab_x[k] * hg_6[k]
                     + ig_6[k];

            t_7[k] = ab_x[k] * hg_7[k]
                     + ig_7[k];

            t_8[k] = ab_x[k] * hg_8[k]
                     + ig_8[k];

            t_9[k] = ab_x[k] * hg_9[k]
                     + ig_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, hg_10, hg_11, hg_12, hg_13, \
                         hg_14, ig_10, ig_11, ig_12, ig_13, ig_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = ab_x[k] * hg_10[k]
                      + ig_10[k];

            t_11[k] = ab_x[k] * hg_11[k]
                      + ig_11[k];

            t_12[k] = ab_x[k] * hg_12[k]
                      + ig_12[k];

            t_13[k] = ab_x[k] * hg_13[k]
                      + ig_13[k];

            t_14[k] = ab_x[k] * hg_14[k]
                      + ig_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_y, hg_10, hg_11, hg_12, hg_13, \
                         hg_14, ig_25, ig_26, ig_27, ig_28, ig_29 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = ab_y[k] * hg_10[k]
                      + ig_25[k];

            t_16[k] = ab_y[k] * hg_11[k]
                      + ig_26[k];

            t_17[k] = ab_y[k] * hg_12[k]
                      + ig_27[k];

            t_18[k] = ab_y[k] * hg_13[k]
                      + ig_28[k];

            t_19[k] = ab_y[k] * hg_14[k]
                      + ig_29[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, ab_x, ab_z, hg_14, hg_15, hg_16, hg_17, \
                         ig_15, ig_16, ig_17, ig_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = ab_z[k] * hg_14[k]
                      + ig_44[k];

            t_21[k] = ab_x[k] * hg_15[k]
                      + ig_15[k];

            t_22[k] = ab_x[k] * hg_16[k]
                      + ig_16[k];

            t_23[k] = ab_x[k] * hg_17[k]
                      + ig_17[k];
        }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, ab_x, hg_18, hg_19, hg_20, hg_21, \
                         hg_22, ig_18, ig_19, ig_20, ig_21, ig_22 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_24[k] = ab_x[k] * hg_18[k]
                      + ig_18[k];

            t_25[k] = ab_x[k] * hg_19[k]
                      + ig_19[k];

            t_26[k] = ab_x[k] * hg_20[k]
                      + ig_20[k];

            t_27[k] = ab_x[k] * hg_21[k]
                      + ig_21[k];

            t_28[k] = ab_x[k] * hg_22[k]
                      + ig_22[k];
        }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, ab_x, hg_23, hg_24, hg_25, hg_26, \
                         hg_27, ig_23, ig_24, ig_25, ig_26, ig_27 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_29[k] = ab_x[k] * hg_23[k]
                      + ig_23[k];

            t_30[k] = ab_x[k] * hg_24[k]
                      + ig_24[k];

            t_31[k] = ab_x[k] * hg_25[k]
                      + ig_25[k];

            t_32[k] = ab_x[k] * hg_26[k]
                      + ig_26[k];

            t_33[k] = ab_x[k] * hg_27[k]
                      + ig_27[k];
        }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, ab_x, ab_y, hg_25, hg_26, hg_28, hg_29, \
                         ig_28, ig_29, ig_55, ig_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_34[k] = ab_x[k] * hg_28[k]
                      + ig_28[k];

            t_35[k] = ab_x[k] * hg_29[k]
                      + ig_29[k];

            t_36[k] = ab_y[k] * hg_25[k]
                      + ig_55[k];

            t_37[k] = ab_y[k] * hg_26[k]
                      + ig_56[k];
        }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, ab_y, ab_z, hg_27, hg_28, hg_29, ig_57, \
                         ig_58, ig_59, ig_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_38[k] = ab_y[k] * hg_27[k]
                      + ig_57[k];

            t_39[k] = ab_y[k] * hg_28[k]
                      + ig_58[k];

            t_40[k] = ab_y[k] * hg_29[k]
                      + ig_59[k];

            t_41[k] = ab_z[k] * hg_29[k]
                      + ig_74[k];
        }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, ab_x, hg_30, hg_31, hg_32, hg_33, \
                         hg_34, ig_30, ig_31, ig_32, ig_33, ig_34 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_42[k] = ab_x[k] * hg_30[k]
                      + ig_30[k];

            t_43[k] = ab_x[k] * hg_31[k]
                      + ig_31[k];

            t_44[k] = ab_x[k] * hg_32[k]
                      + ig_32[k];

            t_45[k] = ab_x[k] * hg_33[k]
                      + ig_33[k];

            t_46[k] = ab_x[k] * hg_34[k]
                      + ig_34[k];
        }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, ab_x, hg_35, hg_36, hg_37, hg_38, \
                         hg_39, ig_35, ig_36, ig_37, ig_38, ig_39 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_47[k] = ab_x[k] * hg_35[k]
                      + ig_35[k];

            t_48[k] = ab_x[k] * hg_36[k]
                      + ig_36[k];

            t_49[k] = ab_x[k] * hg_37[k]
                      + ig_37[k];

            t_50[k] = ab_x[k] * hg_38[k]
                      + ig_38[k];

            t_51[k] = ab_x[k] * hg_39[k]
                      + ig_39[k];
        }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, ab_x, hg_40, hg_41, hg_42, hg_43, \
                         hg_44, ig_40, ig_41, ig_42, ig_43, ig_44 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_52[k] = ab_x[k] * hg_40[k]
                      + ig_40[k];

            t_53[k] = ab_x[k] * hg_41[k]
                      + ig_41[k];

            t_54[k] = ab_x[k] * hg_42[k]
                      + ig_42[k];

            t_55[k] = ab_x[k] * hg_43[k]
                      + ig_43[k];

            t_56[k] = ab_x[k] * hg_44[k]
                      + ig_44[k];
        }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, ab_y, hg_40, hg_41, hg_42, hg_43, \
                         hg_44, ig_70, ig_71, ig_72, ig_73, ig_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_57[k] = ab_y[k] * hg_40[k]
                      + ig_70[k];

            t_58[k] = ab_y[k] * hg_41[k]
                      + ig_71[k];

            t_59[k] = ab_y[k] * hg_42[k]
                      + ig_72[k];

            t_60[k] = ab_y[k] * hg_43[k]
                      + ig_73[k];

            t_61[k] = ab_y[k] * hg_44[k]
                      + ig_74[k];
        }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, ab_x, ab_z, hg_44, hg_45, hg_46, hg_47, \
                         ig_45, ig_46, ig_47, ig_89 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_62[k] = ab_z[k] * hg_44[k]
                      + ig_89[k];

            t_63[k] = ab_x[k] * hg_45[k]
                      + ig_45[k];

            t_64[k] = ab_x[k] * hg_46[k]
                      + ig_46[k];

            t_65[k] = ab_x[k] * hg_47[k]
                      + ig_47[k];
        }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, ab_x, hg_48, hg_49, hg_50, hg_51, \
                         hg_52, ig_48, ig_49, ig_50, ig_51, ig_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_66[k] = ab_x[k] * hg_48[k]
                      + ig_48[k];

            t_67[k] = ab_x[k] * hg_49[k]
                      + ig_49[k];

            t_68[k] = ab_x[k] * hg_50[k]
                      + ig_50[k];

            t_69[k] = ab_x[k] * hg_51[k]
                      + ig_51[k];

            t_70[k] = ab_x[k] * hg_52[k]
                      + ig_52[k];
        }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, ab_x, hg_53, hg_54, hg_55, hg_56, \
                         hg_57, ig_53, ig_54, ig_55, ig_56, ig_57 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_71[k] = ab_x[k] * hg_53[k]
                      + ig_53[k];

            t_72[k] = ab_x[k] * hg_54[k]
                      + ig_54[k];

            t_73[k] = ab_x[k] * hg_55[k]
                      + ig_55[k];

            t_74[k] = ab_x[k] * hg_56[k]
                      + ig_56[k];

            t_75[k] = ab_x[k] * hg_57[k]
                      + ig_57[k];
        }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, ab_x, ab_y, hg_55, hg_56, hg_58, hg_59, \
                         ig_58, ig_59, ig_100, ig_101 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_76[k] = ab_x[k] * hg_58[k]
                      + ig_58[k];

            t_77[k] = ab_x[k] * hg_59[k]
                      + ig_59[k];

            t_78[k] = ab_y[k] * hg_55[k]
                      + ig_100[k];

            t_79[k] = ab_y[k] * hg_56[k]
                      + ig_101[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, ab_y, ab_z, hg_57, hg_58, hg_59, ig_102, \
                         ig_103, ig_104, ig_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = ab_y[k] * hg_57[k]
                      + ig_102[k];

            t_81[k] = ab_y[k] * hg_58[k]
                      + ig_103[k];

            t_82[k] = ab_y[k] * hg_59[k]
                      + ig_104[k];

            t_83[k] = ab_z[k] * hg_59[k]
                      + ig_119[k];
        }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, ab_x, hg_60, hg_61, hg_62, hg_63, \
                         hg_64, ig_60, ig_61, ig_62, ig_63, ig_64 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_84[k] = ab_x[k] * hg_60[k]
                      + ig_60[k];

            t_85[k] = ab_x[k] * hg_61[k]
                      + ig_61[k];

            t_86[k] = ab_x[k] * hg_62[k]
                      + ig_62[k];

            t_87[k] = ab_x[k] * hg_63[k]
                      + ig_63[k];

            t_88[k] = ab_x[k] * hg_64[k]
                      + ig_64[k];
        }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, ab_x, hg_65, hg_66, hg_67, hg_68, \
                         hg_69, ig_65, ig_66, ig_67, ig_68, ig_69 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_89[k] = ab_x[k] * hg_65[k]
                      + ig_65[k];

            t_90[k] = ab_x[k] * hg_66[k]
                      + ig_66[k];

            t_91[k] = ab_x[k] * hg_67[k]
                      + ig_67[k];

            t_92[k] = ab_x[k] * hg_68[k]
                      + ig_68[k];

            t_93[k] = ab_x[k] * hg_69[k]
                      + ig_69[k];
        }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, ab_x, hg_70, hg_71, hg_72, hg_73, \
                         hg_74, ig_70, ig_71, ig_72, ig_73, ig_74 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_94[k] = ab_x[k] * hg_70[k]
                      + ig_70[k];

            t_95[k] = ab_x[k] * hg_71[k]
                      + ig_71[k];

            t_96[k] = ab_x[k] * hg_72[k]
                      + ig_72[k];

            t_97[k] = ab_x[k] * hg_73[k]
                      + ig_73[k];

            t_98[k] = ab_x[k] * hg_74[k]
                      + ig_74[k];
        }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, ab_y, hg_70, hg_71, hg_72, hg_73, \
                         hg_74, ig_115, ig_116, ig_117, ig_118, \
                         ig_119 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_99[k] = ab_y[k] * hg_70[k]
                      + ig_115[k];

            t_100[k] = ab_y[k] * hg_71[k]
                       + ig_116[k];

            t_101[k] = ab_y[k] * hg_72[k]
                       + ig_117[k];

            t_102[k] = ab_y[k] * hg_73[k]
                       + ig_118[k];

            t_103[k] = ab_y[k] * hg_74[k]
                       + ig_119[k];
        }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, ab_x, ab_z, hg_74, hg_75, hg_76, hg_77, \
                         ig_75, ig_76, ig_77, ig_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_104[k] = ab_z[k] * hg_74[k]
                       + ig_134[k];

            t_105[k] = ab_x[k] * hg_75[k]
                       + ig_75[k];

            t_106[k] = ab_x[k] * hg_76[k]
                       + ig_76[k];

            t_107[k] = ab_x[k] * hg_77[k]
                       + ig_77[k];
        }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, ab_x, hg_78, hg_79, hg_80, hg_81, \
                         hg_82, ig_78, ig_79, ig_80, ig_81, ig_82 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_108[k] = ab_x[k] * hg_78[k]
                       + ig_78[k];

            t_109[k] = ab_x[k] * hg_79[k]
                       + ig_79[k];

            t_110[k] = ab_x[k] * hg_80[k]
                       + ig_80[k];

            t_111[k] = ab_x[k] * hg_81[k]
                       + ig_81[k];

            t_112[k] = ab_x[k] * hg_82[k]
                       + ig_82[k];
        }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, ab_x, hg_83, hg_84, hg_85, hg_86, \
                         hg_87, ig_83, ig_84, ig_85, ig_86, ig_87 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_113[k] = ab_x[k] * hg_83[k]
                       + ig_83[k];

            t_114[k] = ab_x[k] * hg_84[k]
                       + ig_84[k];

            t_115[k] = ab_x[k] * hg_85[k]
                       + ig_85[k];

            t_116[k] = ab_x[k] * hg_86[k]
                       + ig_86[k];

            t_117[k] = ab_x[k] * hg_87[k]
                       + ig_87[k];
        }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, ab_x, ab_y, hg_85, hg_86, hg_88, hg_89, \
                         ig_88, ig_89, ig_130, ig_131 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_118[k] = ab_x[k] * hg_88[k]
                       + ig_88[k];

            t_119[k] = ab_x[k] * hg_89[k]
                       + ig_89[k];

            t_120[k] = ab_y[k] * hg_85[k]
                       + ig_130[k];

            t_121[k] = ab_y[k] * hg_86[k]
                       + ig_131[k];
        }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, ab_y, ab_z, hg_87, hg_88, hg_89, ig_132, \
                         ig_133, ig_134, ig_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_122[k] = ab_y[k] * hg_87[k]
                       + ig_132[k];

            t_123[k] = ab_y[k] * hg_88[k]
                       + ig_133[k];

            t_124[k] = ab_y[k] * hg_89[k]
                       + ig_134[k];

            t_125[k] = ab_z[k] * hg_89[k]
                       + ig_149[k];
        }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, ab_x, hg_90, hg_91, hg_92, hg_93, \
                         hg_94, ig_90, ig_91, ig_92, ig_93, ig_94 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_126[k] = ab_x[k] * hg_90[k]
                       + ig_90[k];

            t_127[k] = ab_x[k] * hg_91[k]
                       + ig_91[k];

            t_128[k] = ab_x[k] * hg_92[k]
                       + ig_92[k];

            t_129[k] = ab_x[k] * hg_93[k]
                       + ig_93[k];

            t_130[k] = ab_x[k] * hg_94[k]
                       + ig_94[k];
        }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, ab_x, hg_95, hg_96, hg_97, hg_98, \
                         hg_99, ig_95, ig_96, ig_97, ig_98, ig_99 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_131[k] = ab_x[k] * hg_95[k]
                       + ig_95[k];

            t_132[k] = ab_x[k] * hg_96[k]
                       + ig_96[k];

            t_133[k] = ab_x[k] * hg_97[k]
                       + ig_97[k];

            t_134[k] = ab_x[k] * hg_98[k]
                       + ig_98[k];

            t_135[k] = ab_x[k] * hg_99[k]
                       + ig_99[k];
        }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, ab_x, hg_100, hg_101, hg_102, \
                         hg_103, hg_104, ig_100, ig_101, ig_102, ig_103, \
                         ig_104 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_136[k] = ab_x[k] * hg_100[k]
                       + ig_100[k];

            t_137[k] = ab_x[k] * hg_101[k]
                       + ig_101[k];

            t_138[k] = ab_x[k] * hg_102[k]
                       + ig_102[k];

            t_139[k] = ab_x[k] * hg_103[k]
                       + ig_103[k];

            t_140[k] = ab_x[k] * hg_104[k]
                       + ig_104[k];
        }
    }
}

static auto
compute_hrr_hh_out_of_first_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t hg, const size_t ig,
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

        const auto *hg_100 = buffer.data(hg + 100 * ncomps + c);
        const auto *hg_101 = buffer.data(hg + 101 * ncomps + c);
        const auto *hg_102 = buffer.data(hg + 102 * ncomps + c);
        const auto *hg_103 = buffer.data(hg + 103 * ncomps + c);
        const auto *hg_104 = buffer.data(hg + 104 * ncomps + c);
        const auto *hg_105 = buffer.data(hg + 105 * ncomps + c);
        const auto *hg_106 = buffer.data(hg + 106 * ncomps + c);
        const auto *hg_107 = buffer.data(hg + 107 * ncomps + c);
        const auto *hg_108 = buffer.data(hg + 108 * ncomps + c);
        const auto *hg_109 = buffer.data(hg + 109 * ncomps + c);
        const auto *hg_110 = buffer.data(hg + 110 * ncomps + c);
        const auto *hg_111 = buffer.data(hg + 111 * ncomps + c);
        const auto *hg_112 = buffer.data(hg + 112 * ncomps + c);
        const auto *hg_113 = buffer.data(hg + 113 * ncomps + c);
        const auto *hg_114 = buffer.data(hg + 114 * ncomps + c);
        const auto *hg_115 = buffer.data(hg + 115 * ncomps + c);
        const auto *hg_116 = buffer.data(hg + 116 * ncomps + c);
        const auto *hg_117 = buffer.data(hg + 117 * ncomps + c);
        const auto *hg_118 = buffer.data(hg + 118 * ncomps + c);
        const auto *hg_119 = buffer.data(hg + 119 * ncomps + c);
        const auto *hg_120 = buffer.data(hg + 120 * ncomps + c);
        const auto *hg_121 = buffer.data(hg + 121 * ncomps + c);
        const auto *hg_122 = buffer.data(hg + 122 * ncomps + c);
        const auto *hg_123 = buffer.data(hg + 123 * ncomps + c);
        const auto *hg_124 = buffer.data(hg + 124 * ncomps + c);
        const auto *hg_125 = buffer.data(hg + 125 * ncomps + c);
        const auto *hg_126 = buffer.data(hg + 126 * ncomps + c);
        const auto *hg_127 = buffer.data(hg + 127 * ncomps + c);
        const auto *hg_128 = buffer.data(hg + 128 * ncomps + c);
        const auto *hg_129 = buffer.data(hg + 129 * ncomps + c);
        const auto *hg_130 = buffer.data(hg + 130 * ncomps + c);
        const auto *hg_131 = buffer.data(hg + 131 * ncomps + c);
        const auto *hg_132 = buffer.data(hg + 132 * ncomps + c);
        const auto *hg_133 = buffer.data(hg + 133 * ncomps + c);
        const auto *hg_134 = buffer.data(hg + 134 * ncomps + c);
        const auto *hg_135 = buffer.data(hg + 135 * ncomps + c);
        const auto *hg_136 = buffer.data(hg + 136 * ncomps + c);
        const auto *hg_137 = buffer.data(hg + 137 * ncomps + c);
        const auto *hg_138 = buffer.data(hg + 138 * ncomps + c);
        const auto *hg_139 = buffer.data(hg + 139 * ncomps + c);
        const auto *hg_140 = buffer.data(hg + 140 * ncomps + c);
        const auto *hg_141 = buffer.data(hg + 141 * ncomps + c);
        const auto *hg_142 = buffer.data(hg + 142 * ncomps + c);
        const auto *hg_143 = buffer.data(hg + 143 * ncomps + c);
        const auto *hg_144 = buffer.data(hg + 144 * ncomps + c);
        const auto *hg_145 = buffer.data(hg + 145 * ncomps + c);
        const auto *hg_146 = buffer.data(hg + 146 * ncomps + c);
        const auto *hg_147 = buffer.data(hg + 147 * ncomps + c);
        const auto *hg_148 = buffer.data(hg + 148 * ncomps + c);
        const auto *hg_149 = buffer.data(hg + 149 * ncomps + c);
        const auto *hg_150 = buffer.data(hg + 150 * ncomps + c);
        const auto *hg_151 = buffer.data(hg + 151 * ncomps + c);
        const auto *hg_152 = buffer.data(hg + 152 * ncomps + c);
        const auto *hg_153 = buffer.data(hg + 153 * ncomps + c);
        const auto *hg_154 = buffer.data(hg + 154 * ncomps + c);
        const auto *hg_155 = buffer.data(hg + 155 * ncomps + c);
        const auto *hg_156 = buffer.data(hg + 156 * ncomps + c);
        const auto *hg_157 = buffer.data(hg + 157 * ncomps + c);
        const auto *hg_158 = buffer.data(hg + 158 * ncomps + c);
        const auto *hg_159 = buffer.data(hg + 159 * ncomps + c);
        const auto *hg_160 = buffer.data(hg + 160 * ncomps + c);
        const auto *hg_161 = buffer.data(hg + 161 * ncomps + c);
        const auto *hg_162 = buffer.data(hg + 162 * ncomps + c);
        const auto *hg_163 = buffer.data(hg + 163 * ncomps + c);
        const auto *hg_164 = buffer.data(hg + 164 * ncomps + c);
        const auto *hg_165 = buffer.data(hg + 165 * ncomps + c);
        const auto *hg_166 = buffer.data(hg + 166 * ncomps + c);
        const auto *hg_167 = buffer.data(hg + 167 * ncomps + c);
        const auto *hg_168 = buffer.data(hg + 168 * ncomps + c);
        const auto *hg_169 = buffer.data(hg + 169 * ncomps + c);
        const auto *hg_170 = buffer.data(hg + 170 * ncomps + c);
        const auto *hg_171 = buffer.data(hg + 171 * ncomps + c);
        const auto *hg_172 = buffer.data(hg + 172 * ncomps + c);
        const auto *hg_173 = buffer.data(hg + 173 * ncomps + c);
        const auto *hg_174 = buffer.data(hg + 174 * ncomps + c);
        const auto *hg_175 = buffer.data(hg + 175 * ncomps + c);
        const auto *hg_176 = buffer.data(hg + 176 * ncomps + c);
        const auto *hg_177 = buffer.data(hg + 177 * ncomps + c);
        const auto *hg_178 = buffer.data(hg + 178 * ncomps + c);
        const auto *hg_179 = buffer.data(hg + 179 * ncomps + c);
        const auto *hg_180 = buffer.data(hg + 180 * ncomps + c);
        const auto *hg_181 = buffer.data(hg + 181 * ncomps + c);
        const auto *hg_182 = buffer.data(hg + 182 * ncomps + c);
        const auto *hg_183 = buffer.data(hg + 183 * ncomps + c);
        const auto *hg_184 = buffer.data(hg + 184 * ncomps + c);
        const auto *hg_185 = buffer.data(hg + 185 * ncomps + c);
        const auto *hg_186 = buffer.data(hg + 186 * ncomps + c);
        const auto *hg_187 = buffer.data(hg + 187 * ncomps + c);
        const auto *hg_188 = buffer.data(hg + 188 * ncomps + c);
        const auto *hg_189 = buffer.data(hg + 189 * ncomps + c);
        const auto *hg_190 = buffer.data(hg + 190 * ncomps + c);
        const auto *hg_191 = buffer.data(hg + 191 * ncomps + c);
        const auto *hg_192 = buffer.data(hg + 192 * ncomps + c);
        const auto *hg_193 = buffer.data(hg + 193 * ncomps + c);
        const auto *hg_194 = buffer.data(hg + 194 * ncomps + c);
        const auto *hg_195 = buffer.data(hg + 195 * ncomps + c);
        const auto *hg_196 = buffer.data(hg + 196 * ncomps + c);
        const auto *hg_197 = buffer.data(hg + 197 * ncomps + c);
        const auto *hg_198 = buffer.data(hg + 198 * ncomps + c);
        const auto *hg_199 = buffer.data(hg + 199 * ncomps + c);
        const auto *hg_200 = buffer.data(hg + 200 * ncomps + c);
        const auto *hg_201 = buffer.data(hg + 201 * ncomps + c);
        const auto *hg_202 = buffer.data(hg + 202 * ncomps + c);

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
        const auto *ig_205 = buffer.data(ig + 205 * ncomps + c);
        const auto *ig_206 = buffer.data(ig + 206 * ncomps + c);
        const auto *ig_207 = buffer.data(ig + 207 * ncomps + c);
        const auto *ig_208 = buffer.data(ig + 208 * ncomps + c);
        const auto *ig_209 = buffer.data(ig + 209 * ncomps + c);
        const auto *ig_224 = buffer.data(ig + 224 * ncomps + c);
        const auto *ig_235 = buffer.data(ig + 235 * ncomps + c);
        const auto *ig_236 = buffer.data(ig + 236 * ncomps + c);
        const auto *ig_237 = buffer.data(ig + 237 * ncomps + c);
        const auto *ig_238 = buffer.data(ig + 238 * ncomps + c);
        const auto *ig_239 = buffer.data(ig + 239 * ncomps + c);
        const auto *ig_250 = buffer.data(ig + 250 * ncomps + c);
        const auto *ig_251 = buffer.data(ig + 251 * ncomps + c);
        const auto *ig_252 = buffer.data(ig + 252 * ncomps + c);
        const auto *ig_253 = buffer.data(ig + 253 * ncomps + c);
        const auto *ig_254 = buffer.data(ig + 254 * ncomps + c);
        const auto *ig_265 = buffer.data(ig + 265 * ncomps + c);
        const auto *ig_266 = buffer.data(ig + 266 * ncomps + c);
        const auto *ig_267 = buffer.data(ig + 267 * ncomps + c);
        const auto *ig_268 = buffer.data(ig + 268 * ncomps + c);
        const auto *ig_269 = buffer.data(ig + 269 * ncomps + c);
        const auto *ig_284 = buffer.data(ig + 284 * ncomps + c);

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, ab_y, hg_100, hg_101, hg_102, \
                         hg_103, hg_104, ig_160, ig_161, ig_162, ig_163, \
                         ig_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_141[k] = ab_y[k] * hg_100[k]
                       + ig_160[k];

            t_142[k] = ab_y[k] * hg_101[k]
                       + ig_161[k];

            t_143[k] = ab_y[k] * hg_102[k]
                       + ig_162[k];

            t_144[k] = ab_y[k] * hg_103[k]
                       + ig_163[k];

            t_145[k] = ab_y[k] * hg_104[k]
                       + ig_164[k];
        }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, ab_x, ab_z, hg_104, hg_105, hg_106, \
                         hg_107, ig_105, ig_106, ig_107, ig_179 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_146[k] = ab_z[k] * hg_104[k]
                       + ig_179[k];

            t_147[k] = ab_x[k] * hg_105[k]
                       + ig_105[k];

            t_148[k] = ab_x[k] * hg_106[k]
                       + ig_106[k];

            t_149[k] = ab_x[k] * hg_107[k]
                       + ig_107[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, hg_108, hg_109, hg_110, \
                         hg_111, hg_112, ig_108, ig_109, ig_110, ig_111, \
                         ig_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = ab_x[k] * hg_108[k]
                       + ig_108[k];

            t_151[k] = ab_x[k] * hg_109[k]
                       + ig_109[k];

            t_152[k] = ab_x[k] * hg_110[k]
                       + ig_110[k];

            t_153[k] = ab_x[k] * hg_111[k]
                       + ig_111[k];

            t_154[k] = ab_x[k] * hg_112[k]
                       + ig_112[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, hg_113, hg_114, hg_115, \
                         hg_116, hg_117, ig_113, ig_114, ig_115, ig_116, \
                         ig_117 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = ab_x[k] * hg_113[k]
                       + ig_113[k];

            t_156[k] = ab_x[k] * hg_114[k]
                       + ig_114[k];

            t_157[k] = ab_x[k] * hg_115[k]
                       + ig_115[k];

            t_158[k] = ab_x[k] * hg_116[k]
                       + ig_116[k];

            t_159[k] = ab_x[k] * hg_117[k]
                       + ig_117[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, ab_x, ab_y, hg_115, hg_116, hg_118, \
                         hg_119, ig_118, ig_119, ig_175, ig_176 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = ab_x[k] * hg_118[k]
                       + ig_118[k];

            t_161[k] = ab_x[k] * hg_119[k]
                       + ig_119[k];

            t_162[k] = ab_y[k] * hg_115[k]
                       + ig_175[k];

            t_163[k] = ab_y[k] * hg_116[k]
                       + ig_176[k];
        }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, ab_y, ab_z, hg_117, hg_118, hg_119, \
                         ig_177, ig_178, ig_179, ig_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_164[k] = ab_y[k] * hg_117[k]
                       + ig_177[k];

            t_165[k] = ab_y[k] * hg_118[k]
                       + ig_178[k];

            t_166[k] = ab_y[k] * hg_119[k]
                       + ig_179[k];

            t_167[k] = ab_z[k] * hg_119[k]
                       + ig_194[k];
        }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, ab_x, hg_120, hg_121, hg_122, \
                         hg_123, hg_124, ig_120, ig_121, ig_122, ig_123, \
                         ig_124 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_168[k] = ab_x[k] * hg_120[k]
                       + ig_120[k];

            t_169[k] = ab_x[k] * hg_121[k]
                       + ig_121[k];

            t_170[k] = ab_x[k] * hg_122[k]
                       + ig_122[k];

            t_171[k] = ab_x[k] * hg_123[k]
                       + ig_123[k];

            t_172[k] = ab_x[k] * hg_124[k]
                       + ig_124[k];
        }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, ab_x, hg_125, hg_126, hg_127, \
                         hg_128, hg_129, ig_125, ig_126, ig_127, ig_128, \
                         ig_129 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_173[k] = ab_x[k] * hg_125[k]
                       + ig_125[k];

            t_174[k] = ab_x[k] * hg_126[k]
                       + ig_126[k];

            t_175[k] = ab_x[k] * hg_127[k]
                       + ig_127[k];

            t_176[k] = ab_x[k] * hg_128[k]
                       + ig_128[k];

            t_177[k] = ab_x[k] * hg_129[k]
                       + ig_129[k];
        }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, ab_x, hg_130, hg_131, hg_132, \
                         hg_133, hg_134, ig_130, ig_131, ig_132, ig_133, \
                         ig_134 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_178[k] = ab_x[k] * hg_130[k]
                       + ig_130[k];

            t_179[k] = ab_x[k] * hg_131[k]
                       + ig_131[k];

            t_180[k] = ab_x[k] * hg_132[k]
                       + ig_132[k];

            t_181[k] = ab_x[k] * hg_133[k]
                       + ig_133[k];

            t_182[k] = ab_x[k] * hg_134[k]
                       + ig_134[k];
        }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, ab_y, hg_130, hg_131, hg_132, \
                         hg_133, hg_134, ig_190, ig_191, ig_192, ig_193, \
                         ig_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_183[k] = ab_y[k] * hg_130[k]
                       + ig_190[k];

            t_184[k] = ab_y[k] * hg_131[k]
                       + ig_191[k];

            t_185[k] = ab_y[k] * hg_132[k]
                       + ig_192[k];

            t_186[k] = ab_y[k] * hg_133[k]
                       + ig_193[k];

            t_187[k] = ab_y[k] * hg_134[k]
                       + ig_194[k];
        }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, ab_x, ab_z, hg_134, hg_135, hg_136, \
                         hg_137, ig_135, ig_136, ig_137, ig_209 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_188[k] = ab_z[k] * hg_134[k]
                       + ig_209[k];

            t_189[k] = ab_x[k] * hg_135[k]
                       + ig_135[k];

            t_190[k] = ab_x[k] * hg_136[k]
                       + ig_136[k];

            t_191[k] = ab_x[k] * hg_137[k]
                       + ig_137[k];
        }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, ab_x, hg_138, hg_139, hg_140, \
                         hg_141, hg_142, ig_138, ig_139, ig_140, ig_141, \
                         ig_142 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_192[k] = ab_x[k] * hg_138[k]
                       + ig_138[k];

            t_193[k] = ab_x[k] * hg_139[k]
                       + ig_139[k];

            t_194[k] = ab_x[k] * hg_140[k]
                       + ig_140[k];

            t_195[k] = ab_x[k] * hg_141[k]
                       + ig_141[k];

            t_196[k] = ab_x[k] * hg_142[k]
                       + ig_142[k];
        }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, ab_x, hg_143, hg_144, hg_145, \
                         hg_146, hg_147, ig_143, ig_144, ig_145, ig_146, \
                         ig_147 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_197[k] = ab_x[k] * hg_143[k]
                       + ig_143[k];

            t_198[k] = ab_x[k] * hg_144[k]
                       + ig_144[k];

            t_199[k] = ab_x[k] * hg_145[k]
                       + ig_145[k];

            t_200[k] = ab_x[k] * hg_146[k]
                       + ig_146[k];

            t_201[k] = ab_x[k] * hg_147[k]
                       + ig_147[k];
        }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, ab_x, ab_y, hg_145, hg_146, hg_148, \
                         hg_149, ig_148, ig_149, ig_205, ig_206 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_202[k] = ab_x[k] * hg_148[k]
                       + ig_148[k];

            t_203[k] = ab_x[k] * hg_149[k]
                       + ig_149[k];

            t_204[k] = ab_y[k] * hg_145[k]
                       + ig_205[k];

            t_205[k] = ab_y[k] * hg_146[k]
                       + ig_206[k];
        }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, ab_y, ab_z, hg_147, hg_148, hg_149, \
                         ig_207, ig_208, ig_209, ig_224 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_206[k] = ab_y[k] * hg_147[k]
                       + ig_207[k];

            t_207[k] = ab_y[k] * hg_148[k]
                       + ig_208[k];

            t_208[k] = ab_y[k] * hg_149[k]
                       + ig_209[k];

            t_209[k] = ab_z[k] * hg_149[k]
                       + ig_224[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, hg_150, hg_151, hg_152, \
                         hg_153, hg_154, ig_150, ig_151, ig_152, ig_153, \
                         ig_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = ab_x[k] * hg_150[k]
                       + ig_150[k];

            t_211[k] = ab_x[k] * hg_151[k]
                       + ig_151[k];

            t_212[k] = ab_x[k] * hg_152[k]
                       + ig_152[k];

            t_213[k] = ab_x[k] * hg_153[k]
                       + ig_153[k];

            t_214[k] = ab_x[k] * hg_154[k]
                       + ig_154[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, hg_155, hg_156, hg_157, \
                         hg_158, hg_159, ig_155, ig_156, ig_157, ig_158, \
                         ig_159 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = ab_x[k] * hg_155[k]
                       + ig_155[k];

            t_216[k] = ab_x[k] * hg_156[k]
                       + ig_156[k];

            t_217[k] = ab_x[k] * hg_157[k]
                       + ig_157[k];

            t_218[k] = ab_x[k] * hg_158[k]
                       + ig_158[k];

            t_219[k] = ab_x[k] * hg_159[k]
                       + ig_159[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, hg_160, hg_161, hg_162, \
                         hg_163, hg_164, ig_160, ig_161, ig_162, ig_163, \
                         ig_164 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = ab_x[k] * hg_160[k]
                       + ig_160[k];

            t_221[k] = ab_x[k] * hg_161[k]
                       + ig_161[k];

            t_222[k] = ab_x[k] * hg_162[k]
                       + ig_162[k];

            t_223[k] = ab_x[k] * hg_163[k]
                       + ig_163[k];

            t_224[k] = ab_x[k] * hg_164[k]
                       + ig_164[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_y, hg_160, hg_161, hg_162, \
                         hg_163, hg_164, ig_235, ig_236, ig_237, ig_238, \
                         ig_239 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = ab_y[k] * hg_160[k]
                       + ig_235[k];

            t_226[k] = ab_y[k] * hg_161[k]
                       + ig_236[k];

            t_227[k] = ab_y[k] * hg_162[k]
                       + ig_237[k];

            t_228[k] = ab_y[k] * hg_163[k]
                       + ig_238[k];

            t_229[k] = ab_y[k] * hg_164[k]
                       + ig_239[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, ab_x, ab_z, hg_164, hg_165, hg_166, \
                         hg_167, ig_165, ig_166, ig_167, ig_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = ab_z[k] * hg_164[k]
                       + ig_254[k];

            t_231[k] = ab_x[k] * hg_165[k]
                       + ig_165[k];

            t_232[k] = ab_x[k] * hg_166[k]
                       + ig_166[k];

            t_233[k] = ab_x[k] * hg_167[k]
                       + ig_167[k];
        }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, ab_x, hg_168, hg_169, hg_170, \
                         hg_171, hg_172, ig_168, ig_169, ig_170, ig_171, \
                         ig_172 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_234[k] = ab_x[k] * hg_168[k]
                       + ig_168[k];

            t_235[k] = ab_x[k] * hg_169[k]
                       + ig_169[k];

            t_236[k] = ab_x[k] * hg_170[k]
                       + ig_170[k];

            t_237[k] = ab_x[k] * hg_171[k]
                       + ig_171[k];

            t_238[k] = ab_x[k] * hg_172[k]
                       + ig_172[k];
        }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, ab_x, hg_173, hg_174, hg_175, \
                         hg_176, hg_177, ig_173, ig_174, ig_175, ig_176, \
                         ig_177 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_239[k] = ab_x[k] * hg_173[k]
                       + ig_173[k];

            t_240[k] = ab_x[k] * hg_174[k]
                       + ig_174[k];

            t_241[k] = ab_x[k] * hg_175[k]
                       + ig_175[k];

            t_242[k] = ab_x[k] * hg_176[k]
                       + ig_176[k];

            t_243[k] = ab_x[k] * hg_177[k]
                       + ig_177[k];
        }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, ab_x, ab_y, hg_175, hg_176, hg_178, \
                         hg_179, ig_178, ig_179, ig_250, ig_251 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_244[k] = ab_x[k] * hg_178[k]
                       + ig_178[k];

            t_245[k] = ab_x[k] * hg_179[k]
                       + ig_179[k];

            t_246[k] = ab_y[k] * hg_175[k]
                       + ig_250[k];

            t_247[k] = ab_y[k] * hg_176[k]
                       + ig_251[k];
        }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, ab_y, ab_z, hg_177, hg_178, hg_179, \
                         ig_252, ig_253, ig_254, ig_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_248[k] = ab_y[k] * hg_177[k]
                       + ig_252[k];

            t_249[k] = ab_y[k] * hg_178[k]
                       + ig_253[k];

            t_250[k] = ab_y[k] * hg_179[k]
                       + ig_254[k];

            t_251[k] = ab_z[k] * hg_179[k]
                       + ig_269[k];
        }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, ab_x, hg_180, hg_181, hg_182, \
                         hg_183, hg_184, ig_180, ig_181, ig_182, ig_183, \
                         ig_184 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_252[k] = ab_x[k] * hg_180[k]
                       + ig_180[k];

            t_253[k] = ab_x[k] * hg_181[k]
                       + ig_181[k];

            t_254[k] = ab_x[k] * hg_182[k]
                       + ig_182[k];

            t_255[k] = ab_x[k] * hg_183[k]
                       + ig_183[k];

            t_256[k] = ab_x[k] * hg_184[k]
                       + ig_184[k];
        }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, ab_x, hg_185, hg_186, hg_187, \
                         hg_188, hg_189, ig_185, ig_186, ig_187, ig_188, \
                         ig_189 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_257[k] = ab_x[k] * hg_185[k]
                       + ig_185[k];

            t_258[k] = ab_x[k] * hg_186[k]
                       + ig_186[k];

            t_259[k] = ab_x[k] * hg_187[k]
                       + ig_187[k];

            t_260[k] = ab_x[k] * hg_188[k]
                       + ig_188[k];

            t_261[k] = ab_x[k] * hg_189[k]
                       + ig_189[k];
        }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, ab_x, hg_190, hg_191, hg_192, \
                         hg_193, hg_194, ig_190, ig_191, ig_192, ig_193, \
                         ig_194 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_262[k] = ab_x[k] * hg_190[k]
                       + ig_190[k];

            t_263[k] = ab_x[k] * hg_191[k]
                       + ig_191[k];

            t_264[k] = ab_x[k] * hg_192[k]
                       + ig_192[k];

            t_265[k] = ab_x[k] * hg_193[k]
                       + ig_193[k];

            t_266[k] = ab_x[k] * hg_194[k]
                       + ig_194[k];
        }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, ab_y, hg_190, hg_191, hg_192, \
                         hg_193, hg_194, ig_265, ig_266, ig_267, ig_268, \
                         ig_269 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_267[k] = ab_y[k] * hg_190[k]
                       + ig_265[k];

            t_268[k] = ab_y[k] * hg_191[k]
                       + ig_266[k];

            t_269[k] = ab_y[k] * hg_192[k]
                       + ig_267[k];

            t_270[k] = ab_y[k] * hg_193[k]
                       + ig_268[k];

            t_271[k] = ab_y[k] * hg_194[k]
                       + ig_269[k];
        }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, ab_x, ab_z, hg_194, hg_195, hg_196, \
                         hg_197, ig_195, ig_196, ig_197, ig_284 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_272[k] = ab_z[k] * hg_194[k]
                       + ig_284[k];

            t_273[k] = ab_x[k] * hg_195[k]
                       + ig_195[k];

            t_274[k] = ab_x[k] * hg_196[k]
                       + ig_196[k];

            t_275[k] = ab_x[k] * hg_197[k]
                       + ig_197[k];
        }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, ab_x, hg_198, hg_199, hg_200, \
                         hg_201, hg_202, ig_198, ig_199, ig_200, ig_201, \
                         ig_202 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_276[k] = ab_x[k] * hg_198[k]
                       + ig_198[k];

            t_277[k] = ab_x[k] * hg_199[k]
                       + ig_199[k];

            t_278[k] = ab_x[k] * hg_200[k]
                       + ig_200[k];

            t_279[k] = ab_x[k] * hg_201[k]
                       + ig_201[k];

            t_280[k] = ab_x[k] * hg_202[k]
                       + ig_202[k];
        }
    }
}

static auto
compute_hrr_hh_out_of_first_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t hg, const size_t ig,
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

        const auto *hg_203 = buffer.data(hg + 203 * ncomps + c);
        const auto *hg_204 = buffer.data(hg + 204 * ncomps + c);
        const auto *hg_205 = buffer.data(hg + 205 * ncomps + c);
        const auto *hg_206 = buffer.data(hg + 206 * ncomps + c);
        const auto *hg_207 = buffer.data(hg + 207 * ncomps + c);
        const auto *hg_208 = buffer.data(hg + 208 * ncomps + c);
        const auto *hg_209 = buffer.data(hg + 209 * ncomps + c);
        const auto *hg_210 = buffer.data(hg + 210 * ncomps + c);
        const auto *hg_211 = buffer.data(hg + 211 * ncomps + c);
        const auto *hg_212 = buffer.data(hg + 212 * ncomps + c);
        const auto *hg_213 = buffer.data(hg + 213 * ncomps + c);
        const auto *hg_214 = buffer.data(hg + 214 * ncomps + c);
        const auto *hg_215 = buffer.data(hg + 215 * ncomps + c);
        const auto *hg_216 = buffer.data(hg + 216 * ncomps + c);
        const auto *hg_217 = buffer.data(hg + 217 * ncomps + c);
        const auto *hg_218 = buffer.data(hg + 218 * ncomps + c);
        const auto *hg_219 = buffer.data(hg + 219 * ncomps + c);
        const auto *hg_220 = buffer.data(hg + 220 * ncomps + c);
        const auto *hg_221 = buffer.data(hg + 221 * ncomps + c);
        const auto *hg_222 = buffer.data(hg + 222 * ncomps + c);
        const auto *hg_223 = buffer.data(hg + 223 * ncomps + c);
        const auto *hg_224 = buffer.data(hg + 224 * ncomps + c);
        const auto *hg_225 = buffer.data(hg + 225 * ncomps + c);
        const auto *hg_226 = buffer.data(hg + 226 * ncomps + c);
        const auto *hg_227 = buffer.data(hg + 227 * ncomps + c);
        const auto *hg_228 = buffer.data(hg + 228 * ncomps + c);
        const auto *hg_229 = buffer.data(hg + 229 * ncomps + c);
        const auto *hg_230 = buffer.data(hg + 230 * ncomps + c);
        const auto *hg_231 = buffer.data(hg + 231 * ncomps + c);
        const auto *hg_232 = buffer.data(hg + 232 * ncomps + c);
        const auto *hg_233 = buffer.data(hg + 233 * ncomps + c);
        const auto *hg_234 = buffer.data(hg + 234 * ncomps + c);
        const auto *hg_235 = buffer.data(hg + 235 * ncomps + c);
        const auto *hg_236 = buffer.data(hg + 236 * ncomps + c);
        const auto *hg_237 = buffer.data(hg + 237 * ncomps + c);
        const auto *hg_238 = buffer.data(hg + 238 * ncomps + c);
        const auto *hg_239 = buffer.data(hg + 239 * ncomps + c);
        const auto *hg_240 = buffer.data(hg + 240 * ncomps + c);
        const auto *hg_241 = buffer.data(hg + 241 * ncomps + c);
        const auto *hg_242 = buffer.data(hg + 242 * ncomps + c);
        const auto *hg_243 = buffer.data(hg + 243 * ncomps + c);
        const auto *hg_244 = buffer.data(hg + 244 * ncomps + c);
        const auto *hg_245 = buffer.data(hg + 245 * ncomps + c);
        const auto *hg_246 = buffer.data(hg + 246 * ncomps + c);
        const auto *hg_247 = buffer.data(hg + 247 * ncomps + c);
        const auto *hg_248 = buffer.data(hg + 248 * ncomps + c);
        const auto *hg_249 = buffer.data(hg + 249 * ncomps + c);
        const auto *hg_250 = buffer.data(hg + 250 * ncomps + c);
        const auto *hg_251 = buffer.data(hg + 251 * ncomps + c);
        const auto *hg_252 = buffer.data(hg + 252 * ncomps + c);
        const auto *hg_253 = buffer.data(hg + 253 * ncomps + c);
        const auto *hg_254 = buffer.data(hg + 254 * ncomps + c);
        const auto *hg_255 = buffer.data(hg + 255 * ncomps + c);
        const auto *hg_256 = buffer.data(hg + 256 * ncomps + c);
        const auto *hg_257 = buffer.data(hg + 257 * ncomps + c);
        const auto *hg_258 = buffer.data(hg + 258 * ncomps + c);
        const auto *hg_259 = buffer.data(hg + 259 * ncomps + c);
        const auto *hg_260 = buffer.data(hg + 260 * ncomps + c);
        const auto *hg_261 = buffer.data(hg + 261 * ncomps + c);
        const auto *hg_262 = buffer.data(hg + 262 * ncomps + c);
        const auto *hg_263 = buffer.data(hg + 263 * ncomps + c);
        const auto *hg_264 = buffer.data(hg + 264 * ncomps + c);
        const auto *hg_265 = buffer.data(hg + 265 * ncomps + c);
        const auto *hg_266 = buffer.data(hg + 266 * ncomps + c);
        const auto *hg_267 = buffer.data(hg + 267 * ncomps + c);
        const auto *hg_268 = buffer.data(hg + 268 * ncomps + c);
        const auto *hg_269 = buffer.data(hg + 269 * ncomps + c);
        const auto *hg_270 = buffer.data(hg + 270 * ncomps + c);
        const auto *hg_271 = buffer.data(hg + 271 * ncomps + c);
        const auto *hg_272 = buffer.data(hg + 272 * ncomps + c);
        const auto *hg_273 = buffer.data(hg + 273 * ncomps + c);
        const auto *hg_274 = buffer.data(hg + 274 * ncomps + c);
        const auto *hg_275 = buffer.data(hg + 275 * ncomps + c);
        const auto *hg_276 = buffer.data(hg + 276 * ncomps + c);
        const auto *hg_277 = buffer.data(hg + 277 * ncomps + c);
        const auto *hg_278 = buffer.data(hg + 278 * ncomps + c);
        const auto *hg_279 = buffer.data(hg + 279 * ncomps + c);
        const auto *hg_280 = buffer.data(hg + 280 * ncomps + c);
        const auto *hg_281 = buffer.data(hg + 281 * ncomps + c);
        const auto *hg_282 = buffer.data(hg + 282 * ncomps + c);
        const auto *hg_283 = buffer.data(hg + 283 * ncomps + c);
        const auto *hg_284 = buffer.data(hg + 284 * ncomps + c);
        const auto *hg_285 = buffer.data(hg + 285 * ncomps + c);
        const auto *hg_286 = buffer.data(hg + 286 * ncomps + c);
        const auto *hg_287 = buffer.data(hg + 287 * ncomps + c);
        const auto *hg_288 = buffer.data(hg + 288 * ncomps + c);
        const auto *hg_289 = buffer.data(hg + 289 * ncomps + c);
        const auto *hg_290 = buffer.data(hg + 290 * ncomps + c);
        const auto *hg_291 = buffer.data(hg + 291 * ncomps + c);
        const auto *hg_292 = buffer.data(hg + 292 * ncomps + c);
        const auto *hg_293 = buffer.data(hg + 293 * ncomps + c);
        const auto *hg_294 = buffer.data(hg + 294 * ncomps + c);
        const auto *hg_295 = buffer.data(hg + 295 * ncomps + c);
        const auto *hg_296 = buffer.data(hg + 296 * ncomps + c);
        const auto *hg_297 = buffer.data(hg + 297 * ncomps + c);
        const auto *hg_298 = buffer.data(hg + 298 * ncomps + c);
        const auto *hg_299 = buffer.data(hg + 299 * ncomps + c);
        const auto *hg_300 = buffer.data(hg + 300 * ncomps + c);
        const auto *hg_301 = buffer.data(hg + 301 * ncomps + c);
        const auto *hg_302 = buffer.data(hg + 302 * ncomps + c);
        const auto *hg_303 = buffer.data(hg + 303 * ncomps + c);
        const auto *hg_304 = buffer.data(hg + 304 * ncomps + c);

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
        const auto *ig_314 = buffer.data(ig + 314 * ncomps + c);
        const auto *ig_325 = buffer.data(ig + 325 * ncomps + c);
        const auto *ig_326 = buffer.data(ig + 326 * ncomps + c);
        const auto *ig_327 = buffer.data(ig + 327 * ncomps + c);
        const auto *ig_328 = buffer.data(ig + 328 * ncomps + c);
        const auto *ig_329 = buffer.data(ig + 329 * ncomps + c);
        const auto *ig_340 = buffer.data(ig + 340 * ncomps + c);
        const auto *ig_341 = buffer.data(ig + 341 * ncomps + c);
        const auto *ig_342 = buffer.data(ig + 342 * ncomps + c);
        const auto *ig_343 = buffer.data(ig + 343 * ncomps + c);
        const auto *ig_344 = buffer.data(ig + 344 * ncomps + c);
        const auto *ig_355 = buffer.data(ig + 355 * ncomps + c);
        const auto *ig_356 = buffer.data(ig + 356 * ncomps + c);
        const auto *ig_357 = buffer.data(ig + 357 * ncomps + c);
        const auto *ig_358 = buffer.data(ig + 358 * ncomps + c);
        const auto *ig_359 = buffer.data(ig + 359 * ncomps + c);
        const auto *ig_370 = buffer.data(ig + 370 * ncomps + c);
        const auto *ig_371 = buffer.data(ig + 371 * ncomps + c);
        const auto *ig_372 = buffer.data(ig + 372 * ncomps + c);
        const auto *ig_373 = buffer.data(ig + 373 * ncomps + c);
        const auto *ig_374 = buffer.data(ig + 374 * ncomps + c);
        const auto *ig_385 = buffer.data(ig + 385 * ncomps + c);
        const auto *ig_386 = buffer.data(ig + 386 * ncomps + c);
        const auto *ig_387 = buffer.data(ig + 387 * ncomps + c);
        const auto *ig_388 = buffer.data(ig + 388 * ncomps + c);
        const auto *ig_389 = buffer.data(ig + 389 * ncomps + c);
        const auto *ig_404 = buffer.data(ig + 404 * ncomps + c);

#pragma omp simd aligned(t_281, t_282, t_283, t_284, t_285, ab_x, hg_203, hg_204, hg_205, \
                         hg_206, hg_207, ig_203, ig_204, ig_205, ig_206, \
                         ig_207 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_281[k] = ab_x[k] * hg_203[k]
                       + ig_203[k];

            t_282[k] = ab_x[k] * hg_204[k]
                       + ig_204[k];

            t_283[k] = ab_x[k] * hg_205[k]
                       + ig_205[k];

            t_284[k] = ab_x[k] * hg_206[k]
                       + ig_206[k];

            t_285[k] = ab_x[k] * hg_207[k]
                       + ig_207[k];
        }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, ab_x, ab_y, hg_205, hg_206, hg_208, \
                         hg_209, ig_208, ig_209, ig_280, ig_281 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_286[k] = ab_x[k] * hg_208[k]
                       + ig_208[k];

            t_287[k] = ab_x[k] * hg_209[k]
                       + ig_209[k];

            t_288[k] = ab_y[k] * hg_205[k]
                       + ig_280[k];

            t_289[k] = ab_y[k] * hg_206[k]
                       + ig_281[k];
        }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, ab_y, ab_z, hg_207, hg_208, hg_209, \
                         ig_282, ig_283, ig_284, ig_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = ab_y[k] * hg_207[k]
                       + ig_282[k];

            t_291[k] = ab_y[k] * hg_208[k]
                       + ig_283[k];

            t_292[k] = ab_y[k] * hg_209[k]
                       + ig_284[k];

            t_293[k] = ab_z[k] * hg_209[k]
                       + ig_299[k];
        }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, ab_x, hg_210, hg_211, hg_212, \
                         hg_213, hg_214, ig_210, ig_211, ig_212, ig_213, \
                         ig_214 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_294[k] = ab_x[k] * hg_210[k]
                       + ig_210[k];

            t_295[k] = ab_x[k] * hg_211[k]
                       + ig_211[k];

            t_296[k] = ab_x[k] * hg_212[k]
                       + ig_212[k];

            t_297[k] = ab_x[k] * hg_213[k]
                       + ig_213[k];

            t_298[k] = ab_x[k] * hg_214[k]
                       + ig_214[k];
        }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, ab_x, hg_215, hg_216, hg_217, \
                         hg_218, hg_219, ig_215, ig_216, ig_217, ig_218, \
                         ig_219 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_299[k] = ab_x[k] * hg_215[k]
                       + ig_215[k];

            t_300[k] = ab_x[k] * hg_216[k]
                       + ig_216[k];

            t_301[k] = ab_x[k] * hg_217[k]
                       + ig_217[k];

            t_302[k] = ab_x[k] * hg_218[k]
                       + ig_218[k];

            t_303[k] = ab_x[k] * hg_219[k]
                       + ig_219[k];
        }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, ab_x, hg_220, hg_221, hg_222, \
                         hg_223, hg_224, ig_220, ig_221, ig_222, ig_223, \
                         ig_224 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_304[k] = ab_x[k] * hg_220[k]
                       + ig_220[k];

            t_305[k] = ab_x[k] * hg_221[k]
                       + ig_221[k];

            t_306[k] = ab_x[k] * hg_222[k]
                       + ig_222[k];

            t_307[k] = ab_x[k] * hg_223[k]
                       + ig_223[k];

            t_308[k] = ab_x[k] * hg_224[k]
                       + ig_224[k];
        }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, ab_y, hg_220, hg_221, hg_222, \
                         hg_223, hg_224, ig_295, ig_296, ig_297, ig_298, \
                         ig_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_309[k] = ab_y[k] * hg_220[k]
                       + ig_295[k];

            t_310[k] = ab_y[k] * hg_221[k]
                       + ig_296[k];

            t_311[k] = ab_y[k] * hg_222[k]
                       + ig_297[k];

            t_312[k] = ab_y[k] * hg_223[k]
                       + ig_298[k];

            t_313[k] = ab_y[k] * hg_224[k]
                       + ig_299[k];
        }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, ab_x, ab_z, hg_224, hg_225, hg_226, \
                         hg_227, ig_225, ig_226, ig_227, ig_314 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_314[k] = ab_z[k] * hg_224[k]
                       + ig_314[k];

            t_315[k] = ab_x[k] * hg_225[k]
                       + ig_225[k];

            t_316[k] = ab_x[k] * hg_226[k]
                       + ig_226[k];

            t_317[k] = ab_x[k] * hg_227[k]
                       + ig_227[k];
        }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, ab_x, hg_228, hg_229, hg_230, \
                         hg_231, hg_232, ig_228, ig_229, ig_230, ig_231, \
                         ig_232 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_318[k] = ab_x[k] * hg_228[k]
                       + ig_228[k];

            t_319[k] = ab_x[k] * hg_229[k]
                       + ig_229[k];

            t_320[k] = ab_x[k] * hg_230[k]
                       + ig_230[k];

            t_321[k] = ab_x[k] * hg_231[k]
                       + ig_231[k];

            t_322[k] = ab_x[k] * hg_232[k]
                       + ig_232[k];
        }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, ab_x, hg_233, hg_234, hg_235, \
                         hg_236, hg_237, ig_233, ig_234, ig_235, ig_236, \
                         ig_237 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_323[k] = ab_x[k] * hg_233[k]
                       + ig_233[k];

            t_324[k] = ab_x[k] * hg_234[k]
                       + ig_234[k];

            t_325[k] = ab_x[k] * hg_235[k]
                       + ig_235[k];

            t_326[k] = ab_x[k] * hg_236[k]
                       + ig_236[k];

            t_327[k] = ab_x[k] * hg_237[k]
                       + ig_237[k];
        }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, ab_x, ab_y, hg_235, hg_236, hg_238, \
                         hg_239, ig_238, ig_239, ig_325, ig_326 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_328[k] = ab_x[k] * hg_238[k]
                       + ig_238[k];

            t_329[k] = ab_x[k] * hg_239[k]
                       + ig_239[k];

            t_330[k] = ab_y[k] * hg_235[k]
                       + ig_325[k];

            t_331[k] = ab_y[k] * hg_236[k]
                       + ig_326[k];
        }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, ab_y, ab_z, hg_237, hg_238, hg_239, \
                         ig_327, ig_328, ig_329, ig_344 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_332[k] = ab_y[k] * hg_237[k]
                       + ig_327[k];

            t_333[k] = ab_y[k] * hg_238[k]
                       + ig_328[k];

            t_334[k] = ab_y[k] * hg_239[k]
                       + ig_329[k];

            t_335[k] = ab_z[k] * hg_239[k]
                       + ig_344[k];
        }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, t_340, ab_x, hg_240, hg_241, hg_242, \
                         hg_243, hg_244, ig_240, ig_241, ig_242, ig_243, \
                         ig_244 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_336[k] = ab_x[k] * hg_240[k]
                       + ig_240[k];

            t_337[k] = ab_x[k] * hg_241[k]
                       + ig_241[k];

            t_338[k] = ab_x[k] * hg_242[k]
                       + ig_242[k];

            t_339[k] = ab_x[k] * hg_243[k]
                       + ig_243[k];

            t_340[k] = ab_x[k] * hg_244[k]
                       + ig_244[k];
        }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, t_345, ab_x, hg_245, hg_246, hg_247, \
                         hg_248, hg_249, ig_245, ig_246, ig_247, ig_248, \
                         ig_249 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_341[k] = ab_x[k] * hg_245[k]
                       + ig_245[k];

            t_342[k] = ab_x[k] * hg_246[k]
                       + ig_246[k];

            t_343[k] = ab_x[k] * hg_247[k]
                       + ig_247[k];

            t_344[k] = ab_x[k] * hg_248[k]
                       + ig_248[k];

            t_345[k] = ab_x[k] * hg_249[k]
                       + ig_249[k];
        }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, ab_x, hg_250, hg_251, hg_252, \
                         hg_253, hg_254, ig_250, ig_251, ig_252, ig_253, \
                         ig_254 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_346[k] = ab_x[k] * hg_250[k]
                       + ig_250[k];

            t_347[k] = ab_x[k] * hg_251[k]
                       + ig_251[k];

            t_348[k] = ab_x[k] * hg_252[k]
                       + ig_252[k];

            t_349[k] = ab_x[k] * hg_253[k]
                       + ig_253[k];

            t_350[k] = ab_x[k] * hg_254[k]
                       + ig_254[k];
        }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, ab_y, hg_250, hg_251, hg_252, \
                         hg_253, hg_254, ig_340, ig_341, ig_342, ig_343, \
                         ig_344 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_351[k] = ab_y[k] * hg_250[k]
                       + ig_340[k];

            t_352[k] = ab_y[k] * hg_251[k]
                       + ig_341[k];

            t_353[k] = ab_y[k] * hg_252[k]
                       + ig_342[k];

            t_354[k] = ab_y[k] * hg_253[k]
                       + ig_343[k];

            t_355[k] = ab_y[k] * hg_254[k]
                       + ig_344[k];
        }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, ab_x, ab_z, hg_254, hg_255, hg_256, \
                         hg_257, ig_255, ig_256, ig_257, ig_359 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_356[k] = ab_z[k] * hg_254[k]
                       + ig_359[k];

            t_357[k] = ab_x[k] * hg_255[k]
                       + ig_255[k];

            t_358[k] = ab_x[k] * hg_256[k]
                       + ig_256[k];

            t_359[k] = ab_x[k] * hg_257[k]
                       + ig_257[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, hg_258, hg_259, hg_260, \
                         hg_261, hg_262, ig_258, ig_259, ig_260, ig_261, \
                         ig_262 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = ab_x[k] * hg_258[k]
                       + ig_258[k];

            t_361[k] = ab_x[k] * hg_259[k]
                       + ig_259[k];

            t_362[k] = ab_x[k] * hg_260[k]
                       + ig_260[k];

            t_363[k] = ab_x[k] * hg_261[k]
                       + ig_261[k];

            t_364[k] = ab_x[k] * hg_262[k]
                       + ig_262[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, hg_263, hg_264, hg_265, \
                         hg_266, hg_267, ig_263, ig_264, ig_265, ig_266, \
                         ig_267 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = ab_x[k] * hg_263[k]
                       + ig_263[k];

            t_366[k] = ab_x[k] * hg_264[k]
                       + ig_264[k];

            t_367[k] = ab_x[k] * hg_265[k]
                       + ig_265[k];

            t_368[k] = ab_x[k] * hg_266[k]
                       + ig_266[k];

            t_369[k] = ab_x[k] * hg_267[k]
                       + ig_267[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, ab_x, ab_y, hg_265, hg_266, hg_268, \
                         hg_269, ig_268, ig_269, ig_355, ig_356 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = ab_x[k] * hg_268[k]
                       + ig_268[k];

            t_371[k] = ab_x[k] * hg_269[k]
                       + ig_269[k];

            t_372[k] = ab_y[k] * hg_265[k]
                       + ig_355[k];

            t_373[k] = ab_y[k] * hg_266[k]
                       + ig_356[k];
        }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, ab_y, ab_z, hg_267, hg_268, hg_269, \
                         ig_357, ig_358, ig_359, ig_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_374[k] = ab_y[k] * hg_267[k]
                       + ig_357[k];

            t_375[k] = ab_y[k] * hg_268[k]
                       + ig_358[k];

            t_376[k] = ab_y[k] * hg_269[k]
                       + ig_359[k];

            t_377[k] = ab_z[k] * hg_269[k]
                       + ig_374[k];
        }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, t_382, ab_x, hg_270, hg_271, hg_272, \
                         hg_273, hg_274, ig_270, ig_271, ig_272, ig_273, \
                         ig_274 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_378[k] = ab_x[k] * hg_270[k]
                       + ig_270[k];

            t_379[k] = ab_x[k] * hg_271[k]
                       + ig_271[k];

            t_380[k] = ab_x[k] * hg_272[k]
                       + ig_272[k];

            t_381[k] = ab_x[k] * hg_273[k]
                       + ig_273[k];

            t_382[k] = ab_x[k] * hg_274[k]
                       + ig_274[k];
        }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, ab_x, hg_275, hg_276, hg_277, \
                         hg_278, hg_279, ig_275, ig_276, ig_277, ig_278, \
                         ig_279 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_383[k] = ab_x[k] * hg_275[k]
                       + ig_275[k];

            t_384[k] = ab_x[k] * hg_276[k]
                       + ig_276[k];

            t_385[k] = ab_x[k] * hg_277[k]
                       + ig_277[k];

            t_386[k] = ab_x[k] * hg_278[k]
                       + ig_278[k];

            t_387[k] = ab_x[k] * hg_279[k]
                       + ig_279[k];
        }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, t_392, ab_x, hg_280, hg_281, hg_282, \
                         hg_283, hg_284, ig_280, ig_281, ig_282, ig_283, \
                         ig_284 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_388[k] = ab_x[k] * hg_280[k]
                       + ig_280[k];

            t_389[k] = ab_x[k] * hg_281[k]
                       + ig_281[k];

            t_390[k] = ab_x[k] * hg_282[k]
                       + ig_282[k];

            t_391[k] = ab_x[k] * hg_283[k]
                       + ig_283[k];

            t_392[k] = ab_x[k] * hg_284[k]
                       + ig_284[k];
        }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, t_397, ab_y, hg_280, hg_281, hg_282, \
                         hg_283, hg_284, ig_370, ig_371, ig_372, ig_373, \
                         ig_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_393[k] = ab_y[k] * hg_280[k]
                       + ig_370[k];

            t_394[k] = ab_y[k] * hg_281[k]
                       + ig_371[k];

            t_395[k] = ab_y[k] * hg_282[k]
                       + ig_372[k];

            t_396[k] = ab_y[k] * hg_283[k]
                       + ig_373[k];

            t_397[k] = ab_y[k] * hg_284[k]
                       + ig_374[k];
        }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, ab_x, ab_z, hg_284, hg_285, hg_286, \
                         hg_287, ig_285, ig_286, ig_287, ig_389 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_398[k] = ab_z[k] * hg_284[k]
                       + ig_389[k];

            t_399[k] = ab_x[k] * hg_285[k]
                       + ig_285[k];

            t_400[k] = ab_x[k] * hg_286[k]
                       + ig_286[k];

            t_401[k] = ab_x[k] * hg_287[k]
                       + ig_287[k];
        }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, ab_x, hg_288, hg_289, hg_290, \
                         hg_291, hg_292, ig_288, ig_289, ig_290, ig_291, \
                         ig_292 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_402[k] = ab_x[k] * hg_288[k]
                       + ig_288[k];

            t_403[k] = ab_x[k] * hg_289[k]
                       + ig_289[k];

            t_404[k] = ab_x[k] * hg_290[k]
                       + ig_290[k];

            t_405[k] = ab_x[k] * hg_291[k]
                       + ig_291[k];

            t_406[k] = ab_x[k] * hg_292[k]
                       + ig_292[k];
        }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, ab_x, hg_293, hg_294, hg_295, \
                         hg_296, hg_297, ig_293, ig_294, ig_295, ig_296, \
                         ig_297 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_407[k] = ab_x[k] * hg_293[k]
                       + ig_293[k];

            t_408[k] = ab_x[k] * hg_294[k]
                       + ig_294[k];

            t_409[k] = ab_x[k] * hg_295[k]
                       + ig_295[k];

            t_410[k] = ab_x[k] * hg_296[k]
                       + ig_296[k];

            t_411[k] = ab_x[k] * hg_297[k]
                       + ig_297[k];
        }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, ab_x, ab_y, hg_295, hg_296, hg_298, \
                         hg_299, ig_298, ig_299, ig_385, ig_386 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_412[k] = ab_x[k] * hg_298[k]
                       + ig_298[k];

            t_413[k] = ab_x[k] * hg_299[k]
                       + ig_299[k];

            t_414[k] = ab_y[k] * hg_295[k]
                       + ig_385[k];

            t_415[k] = ab_y[k] * hg_296[k]
                       + ig_386[k];
        }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, ab_y, ab_z, hg_297, hg_298, hg_299, \
                         ig_387, ig_388, ig_389, ig_404 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_416[k] = ab_y[k] * hg_297[k]
                       + ig_387[k];

            t_417[k] = ab_y[k] * hg_298[k]
                       + ig_388[k];

            t_418[k] = ab_y[k] * hg_299[k]
                       + ig_389[k];

            t_419[k] = ab_z[k] * hg_299[k]
                       + ig_404[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_x, hg_300, hg_301, hg_302, \
                         hg_303, hg_304, ig_300, ig_301, ig_302, ig_303, \
                         ig_304 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = ab_x[k] * hg_300[k]
                       + ig_300[k];

            t_421[k] = ab_x[k] * hg_301[k]
                       + ig_301[k];

            t_422[k] = ab_x[k] * hg_302[k]
                       + ig_302[k];

            t_423[k] = ab_x[k] * hg_303[k]
                       + ig_303[k];

            t_424[k] = ab_x[k] * hg_304[k]
                       + ig_304[k];
        }
    }
}

static auto
compute_hrr_hh_out_of_first_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                   const size_t target, const size_t hg, const size_t ig,
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

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *hg_305 = buffer.data(hg + 305 * ncomps + c);
        const auto *hg_306 = buffer.data(hg + 306 * ncomps + c);
        const auto *hg_307 = buffer.data(hg + 307 * ncomps + c);
        const auto *hg_308 = buffer.data(hg + 308 * ncomps + c);
        const auto *hg_309 = buffer.data(hg + 309 * ncomps + c);
        const auto *hg_310 = buffer.data(hg + 310 * ncomps + c);
        const auto *hg_311 = buffer.data(hg + 311 * ncomps + c);
        const auto *hg_312 = buffer.data(hg + 312 * ncomps + c);
        const auto *hg_313 = buffer.data(hg + 313 * ncomps + c);
        const auto *hg_314 = buffer.data(hg + 314 * ncomps + c);

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
        const auto *ig_400 = buffer.data(ig + 400 * ncomps + c);
        const auto *ig_401 = buffer.data(ig + 401 * ncomps + c);
        const auto *ig_402 = buffer.data(ig + 402 * ncomps + c);
        const auto *ig_403 = buffer.data(ig + 403 * ncomps + c);
        const auto *ig_404 = buffer.data(ig + 404 * ncomps + c);
        const auto *ig_419 = buffer.data(ig + 419 * ncomps + c);

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_x, hg_305, hg_306, hg_307, \
                         hg_308, hg_309, ig_305, ig_306, ig_307, ig_308, \
                         ig_309 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = ab_x[k] * hg_305[k]
                       + ig_305[k];

            t_426[k] = ab_x[k] * hg_306[k]
                       + ig_306[k];

            t_427[k] = ab_x[k] * hg_307[k]
                       + ig_307[k];

            t_428[k] = ab_x[k] * hg_308[k]
                       + ig_308[k];

            t_429[k] = ab_x[k] * hg_309[k]
                       + ig_309[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_x, hg_310, hg_311, hg_312, \
                         hg_313, hg_314, ig_310, ig_311, ig_312, ig_313, \
                         ig_314 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = ab_x[k] * hg_310[k]
                       + ig_310[k];

            t_431[k] = ab_x[k] * hg_311[k]
                       + ig_311[k];

            t_432[k] = ab_x[k] * hg_312[k]
                       + ig_312[k];

            t_433[k] = ab_x[k] * hg_313[k]
                       + ig_313[k];

            t_434[k] = ab_x[k] * hg_314[k]
                       + ig_314[k];
        }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_y, hg_310, hg_311, hg_312, \
                         hg_313, hg_314, ig_400, ig_401, ig_402, ig_403, \
                         ig_404 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = ab_y[k] * hg_310[k]
                       + ig_400[k];

            t_436[k] = ab_y[k] * hg_311[k]
                       + ig_401[k];

            t_437[k] = ab_y[k] * hg_312[k]
                       + ig_402[k];

            t_438[k] = ab_y[k] * hg_313[k]
                       + ig_403[k];

            t_439[k] = ab_y[k] * hg_314[k]
                       + ig_404[k];
        }

#pragma omp simd aligned(t_440, ab_z, hg_314, ig_419 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = ab_z[k] * hg_314[k]
                       + ig_419[k];
        }
    }
}

auto
compute_hrr_hh_out_of_first(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                            const size_t target, const size_t hg, const size_t ig,
                            const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_hh_out_of_first_piece0(buffer, coordinates, target, hg, ig, ncomps, nmax);

    compute_hrr_hh_out_of_first_piece1(buffer, coordinates, target, hg, ig, ncomps, nmax);

    compute_hrr_hh_out_of_first_piece2(buffer, coordinates, target, hg, ig, ncomps, nmax);

    compute_hrr_hh_out_of_first_piece3(buffer, coordinates, target, hg, ig, ncomps, nmax);
}

static auto
compute_hrr_hh_piece0(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t gh, const size_t gi, const size_t ncomps,
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

        const auto *gh_0 = buffer.data(gh + 0 * ncomps + c);
        const auto *gh_1 = buffer.data(gh + 1 * ncomps + c);
        const auto *gh_2 = buffer.data(gh + 2 * ncomps + c);
        const auto *gh_3 = buffer.data(gh + 3 * ncomps + c);
        const auto *gh_4 = buffer.data(gh + 4 * ncomps + c);
        const auto *gh_5 = buffer.data(gh + 5 * ncomps + c);
        const auto *gh_6 = buffer.data(gh + 6 * ncomps + c);
        const auto *gh_7 = buffer.data(gh + 7 * ncomps + c);
        const auto *gh_8 = buffer.data(gh + 8 * ncomps + c);
        const auto *gh_9 = buffer.data(gh + 9 * ncomps + c);
        const auto *gh_10 = buffer.data(gh + 10 * ncomps + c);
        const auto *gh_11 = buffer.data(gh + 11 * ncomps + c);
        const auto *gh_12 = buffer.data(gh + 12 * ncomps + c);
        const auto *gh_13 = buffer.data(gh + 13 * ncomps + c);
        const auto *gh_14 = buffer.data(gh + 14 * ncomps + c);
        const auto *gh_15 = buffer.data(gh + 15 * ncomps + c);
        const auto *gh_16 = buffer.data(gh + 16 * ncomps + c);
        const auto *gh_17 = buffer.data(gh + 17 * ncomps + c);
        const auto *gh_18 = buffer.data(gh + 18 * ncomps + c);
        const auto *gh_19 = buffer.data(gh + 19 * ncomps + c);
        const auto *gh_20 = buffer.data(gh + 20 * ncomps + c);
        const auto *gh_21 = buffer.data(gh + 21 * ncomps + c);
        const auto *gh_22 = buffer.data(gh + 22 * ncomps + c);
        const auto *gh_23 = buffer.data(gh + 23 * ncomps + c);
        const auto *gh_24 = buffer.data(gh + 24 * ncomps + c);
        const auto *gh_25 = buffer.data(gh + 25 * ncomps + c);
        const auto *gh_26 = buffer.data(gh + 26 * ncomps + c);
        const auto *gh_27 = buffer.data(gh + 27 * ncomps + c);
        const auto *gh_28 = buffer.data(gh + 28 * ncomps + c);
        const auto *gh_29 = buffer.data(gh + 29 * ncomps + c);
        const auto *gh_30 = buffer.data(gh + 30 * ncomps + c);
        const auto *gh_31 = buffer.data(gh + 31 * ncomps + c);
        const auto *gh_32 = buffer.data(gh + 32 * ncomps + c);
        const auto *gh_33 = buffer.data(gh + 33 * ncomps + c);
        const auto *gh_34 = buffer.data(gh + 34 * ncomps + c);
        const auto *gh_35 = buffer.data(gh + 35 * ncomps + c);
        const auto *gh_36 = buffer.data(gh + 36 * ncomps + c);
        const auto *gh_37 = buffer.data(gh + 37 * ncomps + c);
        const auto *gh_38 = buffer.data(gh + 38 * ncomps + c);
        const auto *gh_39 = buffer.data(gh + 39 * ncomps + c);
        const auto *gh_40 = buffer.data(gh + 40 * ncomps + c);
        const auto *gh_41 = buffer.data(gh + 41 * ncomps + c);
        const auto *gh_42 = buffer.data(gh + 42 * ncomps + c);
        const auto *gh_43 = buffer.data(gh + 43 * ncomps + c);
        const auto *gh_44 = buffer.data(gh + 44 * ncomps + c);
        const auto *gh_45 = buffer.data(gh + 45 * ncomps + c);
        const auto *gh_46 = buffer.data(gh + 46 * ncomps + c);
        const auto *gh_47 = buffer.data(gh + 47 * ncomps + c);
        const auto *gh_48 = buffer.data(gh + 48 * ncomps + c);
        const auto *gh_49 = buffer.data(gh + 49 * ncomps + c);
        const auto *gh_50 = buffer.data(gh + 50 * ncomps + c);
        const auto *gh_51 = buffer.data(gh + 51 * ncomps + c);
        const auto *gh_52 = buffer.data(gh + 52 * ncomps + c);
        const auto *gh_53 = buffer.data(gh + 53 * ncomps + c);
        const auto *gh_54 = buffer.data(gh + 54 * ncomps + c);
        const auto *gh_55 = buffer.data(gh + 55 * ncomps + c);
        const auto *gh_56 = buffer.data(gh + 56 * ncomps + c);
        const auto *gh_57 = buffer.data(gh + 57 * ncomps + c);
        const auto *gh_58 = buffer.data(gh + 58 * ncomps + c);
        const auto *gh_59 = buffer.data(gh + 59 * ncomps + c);
        const auto *gh_60 = buffer.data(gh + 60 * ncomps + c);
        const auto *gh_61 = buffer.data(gh + 61 * ncomps + c);
        const auto *gh_62 = buffer.data(gh + 62 * ncomps + c);
        const auto *gh_63 = buffer.data(gh + 63 * ncomps + c);
        const auto *gh_64 = buffer.data(gh + 64 * ncomps + c);
        const auto *gh_65 = buffer.data(gh + 65 * ncomps + c);
        const auto *gh_66 = buffer.data(gh + 66 * ncomps + c);
        const auto *gh_67 = buffer.data(gh + 67 * ncomps + c);
        const auto *gh_68 = buffer.data(gh + 68 * ncomps + c);
        const auto *gh_69 = buffer.data(gh + 69 * ncomps + c);
        const auto *gh_70 = buffer.data(gh + 70 * ncomps + c);
        const auto *gh_71 = buffer.data(gh + 71 * ncomps + c);
        const auto *gh_72 = buffer.data(gh + 72 * ncomps + c);
        const auto *gh_73 = buffer.data(gh + 73 * ncomps + c);
        const auto *gh_74 = buffer.data(gh + 74 * ncomps + c);
        const auto *gh_75 = buffer.data(gh + 75 * ncomps + c);
        const auto *gh_76 = buffer.data(gh + 76 * ncomps + c);
        const auto *gh_77 = buffer.data(gh + 77 * ncomps + c);
        const auto *gh_78 = buffer.data(gh + 78 * ncomps + c);
        const auto *gh_79 = buffer.data(gh + 79 * ncomps + c);
        const auto *gh_80 = buffer.data(gh + 80 * ncomps + c);
        const auto *gh_81 = buffer.data(gh + 81 * ncomps + c);
        const auto *gh_82 = buffer.data(gh + 82 * ncomps + c);
        const auto *gh_83 = buffer.data(gh + 83 * ncomps + c);
        const auto *gh_84 = buffer.data(gh + 84 * ncomps + c);
        const auto *gh_85 = buffer.data(gh + 85 * ncomps + c);
        const auto *gh_86 = buffer.data(gh + 86 * ncomps + c);
        const auto *gh_87 = buffer.data(gh + 87 * ncomps + c);
        const auto *gh_88 = buffer.data(gh + 88 * ncomps + c);
        const auto *gh_89 = buffer.data(gh + 89 * ncomps + c);
        const auto *gh_90 = buffer.data(gh + 90 * ncomps + c);
        const auto *gh_91 = buffer.data(gh + 91 * ncomps + c);
        const auto *gh_92 = buffer.data(gh + 92 * ncomps + c);
        const auto *gh_93 = buffer.data(gh + 93 * ncomps + c);
        const auto *gh_94 = buffer.data(gh + 94 * ncomps + c);
        const auto *gh_95 = buffer.data(gh + 95 * ncomps + c);
        const auto *gh_96 = buffer.data(gh + 96 * ncomps + c);
        const auto *gh_97 = buffer.data(gh + 97 * ncomps + c);
        const auto *gh_98 = buffer.data(gh + 98 * ncomps + c);
        const auto *gh_99 = buffer.data(gh + 99 * ncomps + c);
        const auto *gh_100 = buffer.data(gh + 100 * ncomps + c);
        const auto *gh_101 = buffer.data(gh + 101 * ncomps + c);
        const auto *gh_102 = buffer.data(gh + 102 * ncomps + c);
        const auto *gh_103 = buffer.data(gh + 103 * ncomps + c);
        const auto *gh_104 = buffer.data(gh + 104 * ncomps + c);
        const auto *gh_105 = buffer.data(gh + 105 * ncomps + c);
        const auto *gh_106 = buffer.data(gh + 106 * ncomps + c);
        const auto *gh_107 = buffer.data(gh + 107 * ncomps + c);
        const auto *gh_108 = buffer.data(gh + 108 * ncomps + c);
        const auto *gh_109 = buffer.data(gh + 109 * ncomps + c);
        const auto *gh_110 = buffer.data(gh + 110 * ncomps + c);
        const auto *gh_111 = buffer.data(gh + 111 * ncomps + c);
        const auto *gh_112 = buffer.data(gh + 112 * ncomps + c);
        const auto *gh_113 = buffer.data(gh + 113 * ncomps + c);
        const auto *gh_114 = buffer.data(gh + 114 * ncomps + c);
        const auto *gh_115 = buffer.data(gh + 115 * ncomps + c);
        const auto *gh_116 = buffer.data(gh + 116 * ncomps + c);
        const auto *gh_117 = buffer.data(gh + 117 * ncomps + c);
        const auto *gh_118 = buffer.data(gh + 118 * ncomps + c);
        const auto *gh_119 = buffer.data(gh + 119 * ncomps + c);
        const auto *gh_120 = buffer.data(gh + 120 * ncomps + c);
        const auto *gh_121 = buffer.data(gh + 121 * ncomps + c);
        const auto *gh_122 = buffer.data(gh + 122 * ncomps + c);
        const auto *gh_123 = buffer.data(gh + 123 * ncomps + c);
        const auto *gh_124 = buffer.data(gh + 124 * ncomps + c);
        const auto *gh_125 = buffer.data(gh + 125 * ncomps + c);
        const auto *gh_126 = buffer.data(gh + 126 * ncomps + c);
        const auto *gh_127 = buffer.data(gh + 127 * ncomps + c);
        const auto *gh_128 = buffer.data(gh + 128 * ncomps + c);
        const auto *gh_129 = buffer.data(gh + 129 * ncomps + c);
        const auto *gh_130 = buffer.data(gh + 130 * ncomps + c);
        const auto *gh_131 = buffer.data(gh + 131 * ncomps + c);
        const auto *gh_132 = buffer.data(gh + 132 * ncomps + c);
        const auto *gh_133 = buffer.data(gh + 133 * ncomps + c);
        const auto *gh_134 = buffer.data(gh + 134 * ncomps + c);
        const auto *gh_135 = buffer.data(gh + 135 * ncomps + c);
        const auto *gh_136 = buffer.data(gh + 136 * ncomps + c);
        const auto *gh_137 = buffer.data(gh + 137 * ncomps + c);
        const auto *gh_138 = buffer.data(gh + 138 * ncomps + c);
        const auto *gh_139 = buffer.data(gh + 139 * ncomps + c);
        const auto *gh_140 = buffer.data(gh + 140 * ncomps + c);
        const auto *gh_141 = buffer.data(gh + 141 * ncomps + c);
        const auto *gh_142 = buffer.data(gh + 142 * ncomps + c);
        const auto *gh_143 = buffer.data(gh + 143 * ncomps + c);
        const auto *gh_144 = buffer.data(gh + 144 * ncomps + c);

        const auto *gi_0 = buffer.data(gi + 0 * ncomps + c);
        const auto *gi_1 = buffer.data(gi + 1 * ncomps + c);
        const auto *gi_2 = buffer.data(gi + 2 * ncomps + c);
        const auto *gi_3 = buffer.data(gi + 3 * ncomps + c);
        const auto *gi_4 = buffer.data(gi + 4 * ncomps + c);
        const auto *gi_5 = buffer.data(gi + 5 * ncomps + c);
        const auto *gi_6 = buffer.data(gi + 6 * ncomps + c);
        const auto *gi_7 = buffer.data(gi + 7 * ncomps + c);
        const auto *gi_8 = buffer.data(gi + 8 * ncomps + c);
        const auto *gi_9 = buffer.data(gi + 9 * ncomps + c);
        const auto *gi_10 = buffer.data(gi + 10 * ncomps + c);
        const auto *gi_11 = buffer.data(gi + 11 * ncomps + c);
        const auto *gi_12 = buffer.data(gi + 12 * ncomps + c);
        const auto *gi_13 = buffer.data(gi + 13 * ncomps + c);
        const auto *gi_14 = buffer.data(gi + 14 * ncomps + c);
        const auto *gi_15 = buffer.data(gi + 15 * ncomps + c);
        const auto *gi_16 = buffer.data(gi + 16 * ncomps + c);
        const auto *gi_17 = buffer.data(gi + 17 * ncomps + c);
        const auto *gi_18 = buffer.data(gi + 18 * ncomps + c);
        const auto *gi_19 = buffer.data(gi + 19 * ncomps + c);
        const auto *gi_20 = buffer.data(gi + 20 * ncomps + c);
        const auto *gi_28 = buffer.data(gi + 28 * ncomps + c);
        const auto *gi_29 = buffer.data(gi + 29 * ncomps + c);
        const auto *gi_30 = buffer.data(gi + 30 * ncomps + c);
        const auto *gi_31 = buffer.data(gi + 31 * ncomps + c);
        const auto *gi_32 = buffer.data(gi + 32 * ncomps + c);
        const auto *gi_33 = buffer.data(gi + 33 * ncomps + c);
        const auto *gi_34 = buffer.data(gi + 34 * ncomps + c);
        const auto *gi_35 = buffer.data(gi + 35 * ncomps + c);
        const auto *gi_36 = buffer.data(gi + 36 * ncomps + c);
        const auto *gi_37 = buffer.data(gi + 37 * ncomps + c);
        const auto *gi_38 = buffer.data(gi + 38 * ncomps + c);
        const auto *gi_39 = buffer.data(gi + 39 * ncomps + c);
        const auto *gi_40 = buffer.data(gi + 40 * ncomps + c);
        const auto *gi_41 = buffer.data(gi + 41 * ncomps + c);
        const auto *gi_42 = buffer.data(gi + 42 * ncomps + c);
        const auto *gi_43 = buffer.data(gi + 43 * ncomps + c);
        const auto *gi_44 = buffer.data(gi + 44 * ncomps + c);
        const auto *gi_45 = buffer.data(gi + 45 * ncomps + c);
        const auto *gi_46 = buffer.data(gi + 46 * ncomps + c);
        const auto *gi_47 = buffer.data(gi + 47 * ncomps + c);
        const auto *gi_48 = buffer.data(gi + 48 * ncomps + c);
        const auto *gi_56 = buffer.data(gi + 56 * ncomps + c);
        const auto *gi_57 = buffer.data(gi + 57 * ncomps + c);
        const auto *gi_58 = buffer.data(gi + 58 * ncomps + c);
        const auto *gi_59 = buffer.data(gi + 59 * ncomps + c);
        const auto *gi_60 = buffer.data(gi + 60 * ncomps + c);
        const auto *gi_61 = buffer.data(gi + 61 * ncomps + c);
        const auto *gi_62 = buffer.data(gi + 62 * ncomps + c);
        const auto *gi_63 = buffer.data(gi + 63 * ncomps + c);
        const auto *gi_64 = buffer.data(gi + 64 * ncomps + c);
        const auto *gi_65 = buffer.data(gi + 65 * ncomps + c);
        const auto *gi_66 = buffer.data(gi + 66 * ncomps + c);
        const auto *gi_67 = buffer.data(gi + 67 * ncomps + c);
        const auto *gi_68 = buffer.data(gi + 68 * ncomps + c);
        const auto *gi_69 = buffer.data(gi + 69 * ncomps + c);
        const auto *gi_70 = buffer.data(gi + 70 * ncomps + c);
        const auto *gi_71 = buffer.data(gi + 71 * ncomps + c);
        const auto *gi_72 = buffer.data(gi + 72 * ncomps + c);
        const auto *gi_73 = buffer.data(gi + 73 * ncomps + c);
        const auto *gi_74 = buffer.data(gi + 74 * ncomps + c);
        const auto *gi_75 = buffer.data(gi + 75 * ncomps + c);
        const auto *gi_76 = buffer.data(gi + 76 * ncomps + c);
        const auto *gi_84 = buffer.data(gi + 84 * ncomps + c);
        const auto *gi_85 = buffer.data(gi + 85 * ncomps + c);
        const auto *gi_86 = buffer.data(gi + 86 * ncomps + c);
        const auto *gi_87 = buffer.data(gi + 87 * ncomps + c);
        const auto *gi_88 = buffer.data(gi + 88 * ncomps + c);
        const auto *gi_89 = buffer.data(gi + 89 * ncomps + c);
        const auto *gi_90 = buffer.data(gi + 90 * ncomps + c);
        const auto *gi_91 = buffer.data(gi + 91 * ncomps + c);
        const auto *gi_92 = buffer.data(gi + 92 * ncomps + c);
        const auto *gi_93 = buffer.data(gi + 93 * ncomps + c);
        const auto *gi_94 = buffer.data(gi + 94 * ncomps + c);
        const auto *gi_95 = buffer.data(gi + 95 * ncomps + c);
        const auto *gi_96 = buffer.data(gi + 96 * ncomps + c);
        const auto *gi_97 = buffer.data(gi + 97 * ncomps + c);
        const auto *gi_98 = buffer.data(gi + 98 * ncomps + c);
        const auto *gi_99 = buffer.data(gi + 99 * ncomps + c);
        const auto *gi_100 = buffer.data(gi + 100 * ncomps + c);
        const auto *gi_101 = buffer.data(gi + 101 * ncomps + c);
        const auto *gi_102 = buffer.data(gi + 102 * ncomps + c);
        const auto *gi_103 = buffer.data(gi + 103 * ncomps + c);
        const auto *gi_104 = buffer.data(gi + 104 * ncomps + c);
        const auto *gi_112 = buffer.data(gi + 112 * ncomps + c);
        const auto *gi_113 = buffer.data(gi + 113 * ncomps + c);
        const auto *gi_114 = buffer.data(gi + 114 * ncomps + c);
        const auto *gi_115 = buffer.data(gi + 115 * ncomps + c);
        const auto *gi_116 = buffer.data(gi + 116 * ncomps + c);
        const auto *gi_117 = buffer.data(gi + 117 * ncomps + c);
        const auto *gi_118 = buffer.data(gi + 118 * ncomps + c);
        const auto *gi_119 = buffer.data(gi + 119 * ncomps + c);
        const auto *gi_120 = buffer.data(gi + 120 * ncomps + c);
        const auto *gi_121 = buffer.data(gi + 121 * ncomps + c);
        const auto *gi_122 = buffer.data(gi + 122 * ncomps + c);
        const auto *gi_123 = buffer.data(gi + 123 * ncomps + c);
        const auto *gi_124 = buffer.data(gi + 124 * ncomps + c);
        const auto *gi_125 = buffer.data(gi + 125 * ncomps + c);
        const auto *gi_126 = buffer.data(gi + 126 * ncomps + c);
        const auto *gi_127 = buffer.data(gi + 127 * ncomps + c);
        const auto *gi_128 = buffer.data(gi + 128 * ncomps + c);
        const auto *gi_129 = buffer.data(gi + 129 * ncomps + c);
        const auto *gi_130 = buffer.data(gi + 130 * ncomps + c);
        const auto *gi_131 = buffer.data(gi + 131 * ncomps + c);
        const auto *gi_132 = buffer.data(gi + 132 * ncomps + c);
        const auto *gi_140 = buffer.data(gi + 140 * ncomps + c);
        const auto *gi_141 = buffer.data(gi + 141 * ncomps + c);
        const auto *gi_142 = buffer.data(gi + 142 * ncomps + c);
        const auto *gi_143 = buffer.data(gi + 143 * ncomps + c);
        const auto *gi_144 = buffer.data(gi + 144 * ncomps + c);
        const auto *gi_145 = buffer.data(gi + 145 * ncomps + c);
        const auto *gi_146 = buffer.data(gi + 146 * ncomps + c);
        const auto *gi_147 = buffer.data(gi + 147 * ncomps + c);
        const auto *gi_148 = buffer.data(gi + 148 * ncomps + c);
        const auto *gi_149 = buffer.data(gi + 149 * ncomps + c);
        const auto *gi_150 = buffer.data(gi + 150 * ncomps + c);
        const auto *gi_151 = buffer.data(gi + 151 * ncomps + c);
        const auto *gi_152 = buffer.data(gi + 152 * ncomps + c);
        const auto *gi_153 = buffer.data(gi + 153 * ncomps + c);
        const auto *gi_154 = buffer.data(gi + 154 * ncomps + c);
        const auto *gi_155 = buffer.data(gi + 155 * ncomps + c);
        const auto *gi_156 = buffer.data(gi + 156 * ncomps + c);
        const auto *gi_157 = buffer.data(gi + 157 * ncomps + c);
        const auto *gi_158 = buffer.data(gi + 158 * ncomps + c);
        const auto *gi_159 = buffer.data(gi + 159 * ncomps + c);
        const auto *gi_160 = buffer.data(gi + 160 * ncomps + c);
        const auto *gi_168 = buffer.data(gi + 168 * ncomps + c);
        const auto *gi_169 = buffer.data(gi + 169 * ncomps + c);
        const auto *gi_170 = buffer.data(gi + 170 * ncomps + c);
        const auto *gi_171 = buffer.data(gi + 171 * ncomps + c);
        const auto *gi_172 = buffer.data(gi + 172 * ncomps + c);
        const auto *gi_173 = buffer.data(gi + 173 * ncomps + c);
        const auto *gi_174 = buffer.data(gi + 174 * ncomps + c);
        const auto *gi_175 = buffer.data(gi + 175 * ncomps + c);
        const auto *gi_176 = buffer.data(gi + 176 * ncomps + c);
        const auto *gi_177 = buffer.data(gi + 177 * ncomps + c);
        const auto *gi_178 = buffer.data(gi + 178 * ncomps + c);
        const auto *gi_179 = buffer.data(gi + 179 * ncomps + c);
        const auto *gi_180 = buffer.data(gi + 180 * ncomps + c);
        const auto *gi_181 = buffer.data(gi + 181 * ncomps + c);
        const auto *gi_182 = buffer.data(gi + 182 * ncomps + c);
        const auto *gi_183 = buffer.data(gi + 183 * ncomps + c);
        const auto *gi_184 = buffer.data(gi + 184 * ncomps + c);
        const auto *gi_185 = buffer.data(gi + 185 * ncomps + c);
        const auto *gi_186 = buffer.data(gi + 186 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, gh_0, gh_1, gh_2, gh_3, gh_4, gi_0, \
                         gi_1, gi_2, gi_3, gi_4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * gh_0[k]
                     + gi_0[k];

            t_1[k] = -ab_x[k] * gh_1[k]
                     + gi_1[k];

            t_2[k] = -ab_x[k] * gh_2[k]
                     + gi_2[k];

            t_3[k] = -ab_x[k] * gh_3[k]
                     + gi_3[k];

            t_4[k] = -ab_x[k] * gh_4[k]
                     + gi_4[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, gh_5, gh_6, gh_7, gh_8, gh_9, gi_5, \
                         gi_6, gi_7, gi_8, gi_9 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * gh_5[k]
                     + gi_5[k];

            t_6[k] = -ab_x[k] * gh_6[k]
                     + gi_6[k];

            t_7[k] = -ab_x[k] * gh_7[k]
                     + gi_7[k];

            t_8[k] = -ab_x[k] * gh_8[k]
                     + gi_8[k];

            t_9[k] = -ab_x[k] * gh_9[k]
                     + gi_9[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, gh_10, gh_11, gh_12, gh_13, \
                         gh_14, gi_10, gi_11, gi_12, gi_13, gi_14 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * gh_10[k]
                      + gi_10[k];

            t_11[k] = -ab_x[k] * gh_11[k]
                      + gi_11[k];

            t_12[k] = -ab_x[k] * gh_12[k]
                      + gi_12[k];

            t_13[k] = -ab_x[k] * gh_13[k]
                      + gi_13[k];

            t_14[k] = -ab_x[k] * gh_14[k]
                      + gi_14[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, gh_15, gh_16, gh_17, gh_18, \
                         gh_19, gi_15, gi_16, gi_17, gi_18, gi_19 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * gh_15[k]
                      + gi_15[k];

            t_16[k] = -ab_x[k] * gh_16[k]
                      + gi_16[k];

            t_17[k] = -ab_x[k] * gh_17[k]
                      + gi_17[k];

            t_18[k] = -ab_x[k] * gh_18[k]
                      + gi_18[k];

            t_19[k] = -ab_x[k] * gh_19[k]
                      + gi_19[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, gh_20, gh_21, gh_22, gh_23, \
                         gh_24, gi_20, gi_28, gi_29, gi_30, gi_31 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * gh_20[k]
                      + gi_20[k];

            t_21[k] = -ab_x[k] * gh_21[k]
                      + gi_28[k];

            t_22[k] = -ab_x[k] * gh_22[k]
                      + gi_29[k];

            t_23[k] = -ab_x[k] * gh_23[k]
                      + gi_30[k];

            t_24[k] = -ab_x[k] * gh_24[k]
                      + gi_31[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, gh_25, gh_26, gh_27, gh_28, \
                         gh_29, gi_32, gi_33, gi_34, gi_35, gi_36 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * gh_25[k]
                      + gi_32[k];

            t_26[k] = -ab_x[k] * gh_26[k]
                      + gi_33[k];

            t_27[k] = -ab_x[k] * gh_27[k]
                      + gi_34[k];

            t_28[k] = -ab_x[k] * gh_28[k]
                      + gi_35[k];

            t_29[k] = -ab_x[k] * gh_29[k]
                      + gi_36[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, gh_30, gh_31, gh_32, gh_33, \
                         gh_34, gi_37, gi_38, gi_39, gi_40, gi_41 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_x[k] * gh_30[k]
                      + gi_37[k];

            t_31[k] = -ab_x[k] * gh_31[k]
                      + gi_38[k];

            t_32[k] = -ab_x[k] * gh_32[k]
                      + gi_39[k];

            t_33[k] = -ab_x[k] * gh_33[k]
                      + gi_40[k];

            t_34[k] = -ab_x[k] * gh_34[k]
                      + gi_41[k];
        }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, gh_35, gh_36, gh_37, gh_38, \
                         gh_39, gi_42, gi_43, gi_44, gi_45, gi_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_35[k] = -ab_x[k] * gh_35[k]
                      + gi_42[k];

            t_36[k] = -ab_x[k] * gh_36[k]
                      + gi_43[k];

            t_37[k] = -ab_x[k] * gh_37[k]
                      + gi_44[k];

            t_38[k] = -ab_x[k] * gh_38[k]
                      + gi_45[k];

            t_39[k] = -ab_x[k] * gh_39[k]
                      + gi_46[k];
        }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, gh_40, gh_41, gh_42, gh_43, \
                         gh_44, gi_47, gi_48, gi_56, gi_57, gi_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_40[k] = -ab_x[k] * gh_40[k]
                      + gi_47[k];

            t_41[k] = -ab_x[k] * gh_41[k]
                      + gi_48[k];

            t_42[k] = -ab_x[k] * gh_42[k]
                      + gi_56[k];

            t_43[k] = -ab_x[k] * gh_43[k]
                      + gi_57[k];

            t_44[k] = -ab_x[k] * gh_44[k]
                      + gi_58[k];
        }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, gh_45, gh_46, gh_47, gh_48, \
                         gh_49, gi_59, gi_60, gi_61, gi_62, gi_63 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_45[k] = -ab_x[k] * gh_45[k]
                      + gi_59[k];

            t_46[k] = -ab_x[k] * gh_46[k]
                      + gi_60[k];

            t_47[k] = -ab_x[k] * gh_47[k]
                      + gi_61[k];

            t_48[k] = -ab_x[k] * gh_48[k]
                      + gi_62[k];

            t_49[k] = -ab_x[k] * gh_49[k]
                      + gi_63[k];
        }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, gh_50, gh_51, gh_52, gh_53, \
                         gh_54, gi_64, gi_65, gi_66, gi_67, gi_68 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_50[k] = -ab_x[k] * gh_50[k]
                      + gi_64[k];

            t_51[k] = -ab_x[k] * gh_51[k]
                      + gi_65[k];

            t_52[k] = -ab_x[k] * gh_52[k]
                      + gi_66[k];

            t_53[k] = -ab_x[k] * gh_53[k]
                      + gi_67[k];

            t_54[k] = -ab_x[k] * gh_54[k]
                      + gi_68[k];
        }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, gh_55, gh_56, gh_57, gh_58, \
                         gh_59, gi_69, gi_70, gi_71, gi_72, gi_73 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_55[k] = -ab_x[k] * gh_55[k]
                      + gi_69[k];

            t_56[k] = -ab_x[k] * gh_56[k]
                      + gi_70[k];

            t_57[k] = -ab_x[k] * gh_57[k]
                      + gi_71[k];

            t_58[k] = -ab_x[k] * gh_58[k]
                      + gi_72[k];

            t_59[k] = -ab_x[k] * gh_59[k]
                      + gi_73[k];
        }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, gh_60, gh_61, gh_62, gh_63, \
                         gh_64, gi_74, gi_75, gi_76, gi_84, gi_85 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_60[k] = -ab_x[k] * gh_60[k]
                      + gi_74[k];

            t_61[k] = -ab_x[k] * gh_61[k]
                      + gi_75[k];

            t_62[k] = -ab_x[k] * gh_62[k]
                      + gi_76[k];

            t_63[k] = -ab_x[k] * gh_63[k]
                      + gi_84[k];

            t_64[k] = -ab_x[k] * gh_64[k]
                      + gi_85[k];
        }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, gh_65, gh_66, gh_67, gh_68, \
                         gh_69, gi_86, gi_87, gi_88, gi_89, gi_90 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_65[k] = -ab_x[k] * gh_65[k]
                      + gi_86[k];

            t_66[k] = -ab_x[k] * gh_66[k]
                      + gi_87[k];

            t_67[k] = -ab_x[k] * gh_67[k]
                      + gi_88[k];

            t_68[k] = -ab_x[k] * gh_68[k]
                      + gi_89[k];

            t_69[k] = -ab_x[k] * gh_69[k]
                      + gi_90[k];
        }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, gh_70, gh_71, gh_72, gh_73, \
                         gh_74, gi_91, gi_92, gi_93, gi_94, gi_95 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_70[k] = -ab_x[k] * gh_70[k]
                      + gi_91[k];

            t_71[k] = -ab_x[k] * gh_71[k]
                      + gi_92[k];

            t_72[k] = -ab_x[k] * gh_72[k]
                      + gi_93[k];

            t_73[k] = -ab_x[k] * gh_73[k]
                      + gi_94[k];

            t_74[k] = -ab_x[k] * gh_74[k]
                      + gi_95[k];
        }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, gh_75, gh_76, gh_77, gh_78, \
                         gh_79, gi_96, gi_97, gi_98, gi_99, gi_100 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_75[k] = -ab_x[k] * gh_75[k]
                      + gi_96[k];

            t_76[k] = -ab_x[k] * gh_76[k]
                      + gi_97[k];

            t_77[k] = -ab_x[k] * gh_77[k]
                      + gi_98[k];

            t_78[k] = -ab_x[k] * gh_78[k]
                      + gi_99[k];

            t_79[k] = -ab_x[k] * gh_79[k]
                      + gi_100[k];
        }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, gh_80, gh_81, gh_82, gh_83, \
                         gh_84, gi_101, gi_102, gi_103, gi_104, \
                         gi_112 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_80[k] = -ab_x[k] * gh_80[k]
                      + gi_101[k];

            t_81[k] = -ab_x[k] * gh_81[k]
                      + gi_102[k];

            t_82[k] = -ab_x[k] * gh_82[k]
                      + gi_103[k];

            t_83[k] = -ab_x[k] * gh_83[k]
                      + gi_104[k];

            t_84[k] = -ab_x[k] * gh_84[k]
                      + gi_112[k];
        }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, gh_85, gh_86, gh_87, gh_88, \
                         gh_89, gi_113, gi_114, gi_115, gi_116, \
                         gi_117 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_85[k] = -ab_x[k] * gh_85[k]
                      + gi_113[k];

            t_86[k] = -ab_x[k] * gh_86[k]
                      + gi_114[k];

            t_87[k] = -ab_x[k] * gh_87[k]
                      + gi_115[k];

            t_88[k] = -ab_x[k] * gh_88[k]
                      + gi_116[k];

            t_89[k] = -ab_x[k] * gh_89[k]
                      + gi_117[k];
        }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, gh_90, gh_91, gh_92, gh_93, \
                         gh_94, gi_118, gi_119, gi_120, gi_121, \
                         gi_122 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_90[k] = -ab_x[k] * gh_90[k]
                      + gi_118[k];

            t_91[k] = -ab_x[k] * gh_91[k]
                      + gi_119[k];

            t_92[k] = -ab_x[k] * gh_92[k]
                      + gi_120[k];

            t_93[k] = -ab_x[k] * gh_93[k]
                      + gi_121[k];

            t_94[k] = -ab_x[k] * gh_94[k]
                      + gi_122[k];
        }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, gh_95, gh_96, gh_97, gh_98, \
                         gh_99, gi_123, gi_124, gi_125, gi_126, \
                         gi_127 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_95[k] = -ab_x[k] * gh_95[k]
                      + gi_123[k];

            t_96[k] = -ab_x[k] * gh_96[k]
                      + gi_124[k];

            t_97[k] = -ab_x[k] * gh_97[k]
                      + gi_125[k];

            t_98[k] = -ab_x[k] * gh_98[k]
                      + gi_126[k];

            t_99[k] = -ab_x[k] * gh_99[k]
                      + gi_127[k];
        }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, gh_100, gh_101, gh_102, \
                         gh_103, gh_104, gi_128, gi_129, gi_130, gi_131, \
                         gi_132 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_100[k] = -ab_x[k] * gh_100[k]
                       + gi_128[k];

            t_101[k] = -ab_x[k] * gh_101[k]
                       + gi_129[k];

            t_102[k] = -ab_x[k] * gh_102[k]
                       + gi_130[k];

            t_103[k] = -ab_x[k] * gh_103[k]
                       + gi_131[k];

            t_104[k] = -ab_x[k] * gh_104[k]
                       + gi_132[k];
        }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, gh_105, gh_106, gh_107, \
                         gh_108, gh_109, gi_140, gi_141, gi_142, gi_143, \
                         gi_144 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_105[k] = -ab_x[k] * gh_105[k]
                       + gi_140[k];

            t_106[k] = -ab_x[k] * gh_106[k]
                       + gi_141[k];

            t_107[k] = -ab_x[k] * gh_107[k]
                       + gi_142[k];

            t_108[k] = -ab_x[k] * gh_108[k]
                       + gi_143[k];

            t_109[k] = -ab_x[k] * gh_109[k]
                       + gi_144[k];
        }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, gh_110, gh_111, gh_112, \
                         gh_113, gh_114, gi_145, gi_146, gi_147, gi_148, \
                         gi_149 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_110[k] = -ab_x[k] * gh_110[k]
                       + gi_145[k];

            t_111[k] = -ab_x[k] * gh_111[k]
                       + gi_146[k];

            t_112[k] = -ab_x[k] * gh_112[k]
                       + gi_147[k];

            t_113[k] = -ab_x[k] * gh_113[k]
                       + gi_148[k];

            t_114[k] = -ab_x[k] * gh_114[k]
                       + gi_149[k];
        }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, gh_115, gh_116, gh_117, \
                         gh_118, gh_119, gi_150, gi_151, gi_152, gi_153, \
                         gi_154 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_115[k] = -ab_x[k] * gh_115[k]
                       + gi_150[k];

            t_116[k] = -ab_x[k] * gh_116[k]
                       + gi_151[k];

            t_117[k] = -ab_x[k] * gh_117[k]
                       + gi_152[k];

            t_118[k] = -ab_x[k] * gh_118[k]
                       + gi_153[k];

            t_119[k] = -ab_x[k] * gh_119[k]
                       + gi_154[k];
        }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, gh_120, gh_121, gh_122, \
                         gh_123, gh_124, gi_155, gi_156, gi_157, gi_158, \
                         gi_159 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_120[k] = -ab_x[k] * gh_120[k]
                       + gi_155[k];

            t_121[k] = -ab_x[k] * gh_121[k]
                       + gi_156[k];

            t_122[k] = -ab_x[k] * gh_122[k]
                       + gi_157[k];

            t_123[k] = -ab_x[k] * gh_123[k]
                       + gi_158[k];

            t_124[k] = -ab_x[k] * gh_124[k]
                       + gi_159[k];
        }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, gh_125, gh_126, gh_127, \
                         gh_128, gh_129, gi_160, gi_168, gi_169, gi_170, \
                         gi_171 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_125[k] = -ab_x[k] * gh_125[k]
                       + gi_160[k];

            t_126[k] = -ab_x[k] * gh_126[k]
                       + gi_168[k];

            t_127[k] = -ab_x[k] * gh_127[k]
                       + gi_169[k];

            t_128[k] = -ab_x[k] * gh_128[k]
                       + gi_170[k];

            t_129[k] = -ab_x[k] * gh_129[k]
                       + gi_171[k];
        }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, gh_130, gh_131, gh_132, \
                         gh_133, gh_134, gi_172, gi_173, gi_174, gi_175, \
                         gi_176 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_130[k] = -ab_x[k] * gh_130[k]
                       + gi_172[k];

            t_131[k] = -ab_x[k] * gh_131[k]
                       + gi_173[k];

            t_132[k] = -ab_x[k] * gh_132[k]
                       + gi_174[k];

            t_133[k] = -ab_x[k] * gh_133[k]
                       + gi_175[k];

            t_134[k] = -ab_x[k] * gh_134[k]
                       + gi_176[k];
        }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, gh_135, gh_136, gh_137, \
                         gh_138, gh_139, gi_177, gi_178, gi_179, gi_180, \
                         gi_181 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_135[k] = -ab_x[k] * gh_135[k]
                       + gi_177[k];

            t_136[k] = -ab_x[k] * gh_136[k]
                       + gi_178[k];

            t_137[k] = -ab_x[k] * gh_137[k]
                       + gi_179[k];

            t_138[k] = -ab_x[k] * gh_138[k]
                       + gi_180[k];

            t_139[k] = -ab_x[k] * gh_139[k]
                       + gi_181[k];
        }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, gh_140, gh_141, gh_142, \
                         gh_143, gh_144, gi_182, gi_183, gi_184, gi_185, \
                         gi_186 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_140[k] = -ab_x[k] * gh_140[k]
                       + gi_182[k];

            t_141[k] = -ab_x[k] * gh_141[k]
                       + gi_183[k];

            t_142[k] = -ab_x[k] * gh_142[k]
                       + gi_184[k];

            t_143[k] = -ab_x[k] * gh_143[k]
                       + gi_185[k];

            t_144[k] = -ab_x[k] * gh_144[k]
                       + gi_186[k];
        }
    }
}

static auto
compute_hrr_hh_piece1(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t gh, const size_t gi, const size_t ncomps,
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

        const auto *gh_145 = buffer.data(gh + 145 * ncomps + c);
        const auto *gh_146 = buffer.data(gh + 146 * ncomps + c);
        const auto *gh_147 = buffer.data(gh + 147 * ncomps + c);
        const auto *gh_148 = buffer.data(gh + 148 * ncomps + c);
        const auto *gh_149 = buffer.data(gh + 149 * ncomps + c);
        const auto *gh_150 = buffer.data(gh + 150 * ncomps + c);
        const auto *gh_151 = buffer.data(gh + 151 * ncomps + c);
        const auto *gh_152 = buffer.data(gh + 152 * ncomps + c);
        const auto *gh_153 = buffer.data(gh + 153 * ncomps + c);
        const auto *gh_154 = buffer.data(gh + 154 * ncomps + c);
        const auto *gh_155 = buffer.data(gh + 155 * ncomps + c);
        const auto *gh_156 = buffer.data(gh + 156 * ncomps + c);
        const auto *gh_157 = buffer.data(gh + 157 * ncomps + c);
        const auto *gh_158 = buffer.data(gh + 158 * ncomps + c);
        const auto *gh_159 = buffer.data(gh + 159 * ncomps + c);
        const auto *gh_160 = buffer.data(gh + 160 * ncomps + c);
        const auto *gh_161 = buffer.data(gh + 161 * ncomps + c);
        const auto *gh_162 = buffer.data(gh + 162 * ncomps + c);
        const auto *gh_163 = buffer.data(gh + 163 * ncomps + c);
        const auto *gh_164 = buffer.data(gh + 164 * ncomps + c);
        const auto *gh_165 = buffer.data(gh + 165 * ncomps + c);
        const auto *gh_166 = buffer.data(gh + 166 * ncomps + c);
        const auto *gh_167 = buffer.data(gh + 167 * ncomps + c);
        const auto *gh_168 = buffer.data(gh + 168 * ncomps + c);
        const auto *gh_169 = buffer.data(gh + 169 * ncomps + c);
        const auto *gh_170 = buffer.data(gh + 170 * ncomps + c);
        const auto *gh_171 = buffer.data(gh + 171 * ncomps + c);
        const auto *gh_172 = buffer.data(gh + 172 * ncomps + c);
        const auto *gh_173 = buffer.data(gh + 173 * ncomps + c);
        const auto *gh_174 = buffer.data(gh + 174 * ncomps + c);
        const auto *gh_175 = buffer.data(gh + 175 * ncomps + c);
        const auto *gh_176 = buffer.data(gh + 176 * ncomps + c);
        const auto *gh_177 = buffer.data(gh + 177 * ncomps + c);
        const auto *gh_178 = buffer.data(gh + 178 * ncomps + c);
        const auto *gh_179 = buffer.data(gh + 179 * ncomps + c);
        const auto *gh_180 = buffer.data(gh + 180 * ncomps + c);
        const auto *gh_181 = buffer.data(gh + 181 * ncomps + c);
        const auto *gh_182 = buffer.data(gh + 182 * ncomps + c);
        const auto *gh_183 = buffer.data(gh + 183 * ncomps + c);
        const auto *gh_184 = buffer.data(gh + 184 * ncomps + c);
        const auto *gh_185 = buffer.data(gh + 185 * ncomps + c);
        const auto *gh_186 = buffer.data(gh + 186 * ncomps + c);
        const auto *gh_187 = buffer.data(gh + 187 * ncomps + c);
        const auto *gh_188 = buffer.data(gh + 188 * ncomps + c);
        const auto *gh_189 = buffer.data(gh + 189 * ncomps + c);
        const auto *gh_190 = buffer.data(gh + 190 * ncomps + c);
        const auto *gh_191 = buffer.data(gh + 191 * ncomps + c);
        const auto *gh_192 = buffer.data(gh + 192 * ncomps + c);
        const auto *gh_193 = buffer.data(gh + 193 * ncomps + c);
        const auto *gh_194 = buffer.data(gh + 194 * ncomps + c);
        const auto *gh_195 = buffer.data(gh + 195 * ncomps + c);
        const auto *gh_196 = buffer.data(gh + 196 * ncomps + c);
        const auto *gh_197 = buffer.data(gh + 197 * ncomps + c);
        const auto *gh_198 = buffer.data(gh + 198 * ncomps + c);
        const auto *gh_199 = buffer.data(gh + 199 * ncomps + c);
        const auto *gh_200 = buffer.data(gh + 200 * ncomps + c);
        const auto *gh_201 = buffer.data(gh + 201 * ncomps + c);
        const auto *gh_202 = buffer.data(gh + 202 * ncomps + c);
        const auto *gh_203 = buffer.data(gh + 203 * ncomps + c);
        const auto *gh_204 = buffer.data(gh + 204 * ncomps + c);
        const auto *gh_205 = buffer.data(gh + 205 * ncomps + c);
        const auto *gh_206 = buffer.data(gh + 206 * ncomps + c);
        const auto *gh_207 = buffer.data(gh + 207 * ncomps + c);
        const auto *gh_208 = buffer.data(gh + 208 * ncomps + c);
        const auto *gh_209 = buffer.data(gh + 209 * ncomps + c);
        const auto *gh_210 = buffer.data(gh + 210 * ncomps + c);
        const auto *gh_211 = buffer.data(gh + 211 * ncomps + c);
        const auto *gh_212 = buffer.data(gh + 212 * ncomps + c);
        const auto *gh_213 = buffer.data(gh + 213 * ncomps + c);
        const auto *gh_214 = buffer.data(gh + 214 * ncomps + c);
        const auto *gh_215 = buffer.data(gh + 215 * ncomps + c);
        const auto *gh_216 = buffer.data(gh + 216 * ncomps + c);
        const auto *gh_217 = buffer.data(gh + 217 * ncomps + c);
        const auto *gh_218 = buffer.data(gh + 218 * ncomps + c);
        const auto *gh_219 = buffer.data(gh + 219 * ncomps + c);
        const auto *gh_220 = buffer.data(gh + 220 * ncomps + c);
        const auto *gh_221 = buffer.data(gh + 221 * ncomps + c);
        const auto *gh_222 = buffer.data(gh + 222 * ncomps + c);
        const auto *gh_223 = buffer.data(gh + 223 * ncomps + c);
        const auto *gh_224 = buffer.data(gh + 224 * ncomps + c);
        const auto *gh_225 = buffer.data(gh + 225 * ncomps + c);
        const auto *gh_226 = buffer.data(gh + 226 * ncomps + c);
        const auto *gh_227 = buffer.data(gh + 227 * ncomps + c);
        const auto *gh_228 = buffer.data(gh + 228 * ncomps + c);
        const auto *gh_229 = buffer.data(gh + 229 * ncomps + c);
        const auto *gh_230 = buffer.data(gh + 230 * ncomps + c);
        const auto *gh_231 = buffer.data(gh + 231 * ncomps + c);
        const auto *gh_232 = buffer.data(gh + 232 * ncomps + c);
        const auto *gh_233 = buffer.data(gh + 233 * ncomps + c);
        const auto *gh_234 = buffer.data(gh + 234 * ncomps + c);
        const auto *gh_235 = buffer.data(gh + 235 * ncomps + c);
        const auto *gh_236 = buffer.data(gh + 236 * ncomps + c);
        const auto *gh_237 = buffer.data(gh + 237 * ncomps + c);
        const auto *gh_238 = buffer.data(gh + 238 * ncomps + c);
        const auto *gh_239 = buffer.data(gh + 239 * ncomps + c);
        const auto *gh_240 = buffer.data(gh + 240 * ncomps + c);
        const auto *gh_241 = buffer.data(gh + 241 * ncomps + c);
        const auto *gh_242 = buffer.data(gh + 242 * ncomps + c);
        const auto *gh_243 = buffer.data(gh + 243 * ncomps + c);
        const auto *gh_244 = buffer.data(gh + 244 * ncomps + c);
        const auto *gh_245 = buffer.data(gh + 245 * ncomps + c);
        const auto *gh_246 = buffer.data(gh + 246 * ncomps + c);
        const auto *gh_247 = buffer.data(gh + 247 * ncomps + c);
        const auto *gh_248 = buffer.data(gh + 248 * ncomps + c);
        const auto *gh_249 = buffer.data(gh + 249 * ncomps + c);
        const auto *gh_250 = buffer.data(gh + 250 * ncomps + c);
        const auto *gh_251 = buffer.data(gh + 251 * ncomps + c);
        const auto *gh_252 = buffer.data(gh + 252 * ncomps + c);
        const auto *gh_253 = buffer.data(gh + 253 * ncomps + c);
        const auto *gh_254 = buffer.data(gh + 254 * ncomps + c);
        const auto *gh_255 = buffer.data(gh + 255 * ncomps + c);
        const auto *gh_256 = buffer.data(gh + 256 * ncomps + c);
        const auto *gh_257 = buffer.data(gh + 257 * ncomps + c);
        const auto *gh_258 = buffer.data(gh + 258 * ncomps + c);
        const auto *gh_259 = buffer.data(gh + 259 * ncomps + c);
        const auto *gh_260 = buffer.data(gh + 260 * ncomps + c);
        const auto *gh_261 = buffer.data(gh + 261 * ncomps + c);
        const auto *gh_262 = buffer.data(gh + 262 * ncomps + c);
        const auto *gh_263 = buffer.data(gh + 263 * ncomps + c);
        const auto *gh_264 = buffer.data(gh + 264 * ncomps + c);
        const auto *gh_265 = buffer.data(gh + 265 * ncomps + c);
        const auto *gh_266 = buffer.data(gh + 266 * ncomps + c);
        const auto *gh_267 = buffer.data(gh + 267 * ncomps + c);
        const auto *gh_268 = buffer.data(gh + 268 * ncomps + c);
        const auto *gh_269 = buffer.data(gh + 269 * ncomps + c);
        const auto *gh_270 = buffer.data(gh + 270 * ncomps + c);
        const auto *gh_271 = buffer.data(gh + 271 * ncomps + c);
        const auto *gh_272 = buffer.data(gh + 272 * ncomps + c);
        const auto *gh_273 = buffer.data(gh + 273 * ncomps + c);
        const auto *gh_274 = buffer.data(gh + 274 * ncomps + c);
        const auto *gh_275 = buffer.data(gh + 275 * ncomps + c);
        const auto *gh_276 = buffer.data(gh + 276 * ncomps + c);
        const auto *gh_277 = buffer.data(gh + 277 * ncomps + c);
        const auto *gh_278 = buffer.data(gh + 278 * ncomps + c);
        const auto *gh_279 = buffer.data(gh + 279 * ncomps + c);
        const auto *gh_280 = buffer.data(gh + 280 * ncomps + c);
        const auto *gh_281 = buffer.data(gh + 281 * ncomps + c);
        const auto *gh_282 = buffer.data(gh + 282 * ncomps + c);
        const auto *gh_283 = buffer.data(gh + 283 * ncomps + c);
        const auto *gh_284 = buffer.data(gh + 284 * ncomps + c);
        const auto *gh_285 = buffer.data(gh + 285 * ncomps + c);
        const auto *gh_286 = buffer.data(gh + 286 * ncomps + c);
        const auto *gh_287 = buffer.data(gh + 287 * ncomps + c);
        const auto *gh_288 = buffer.data(gh + 288 * ncomps + c);
        const auto *gh_289 = buffer.data(gh + 289 * ncomps + c);

        const auto *gi_187 = buffer.data(gi + 187 * ncomps + c);
        const auto *gi_188 = buffer.data(gi + 188 * ncomps + c);
        const auto *gi_196 = buffer.data(gi + 196 * ncomps + c);
        const auto *gi_197 = buffer.data(gi + 197 * ncomps + c);
        const auto *gi_198 = buffer.data(gi + 198 * ncomps + c);
        const auto *gi_199 = buffer.data(gi + 199 * ncomps + c);
        const auto *gi_200 = buffer.data(gi + 200 * ncomps + c);
        const auto *gi_201 = buffer.data(gi + 201 * ncomps + c);
        const auto *gi_202 = buffer.data(gi + 202 * ncomps + c);
        const auto *gi_203 = buffer.data(gi + 203 * ncomps + c);
        const auto *gi_204 = buffer.data(gi + 204 * ncomps + c);
        const auto *gi_205 = buffer.data(gi + 205 * ncomps + c);
        const auto *gi_206 = buffer.data(gi + 206 * ncomps + c);
        const auto *gi_207 = buffer.data(gi + 207 * ncomps + c);
        const auto *gi_208 = buffer.data(gi + 208 * ncomps + c);
        const auto *gi_209 = buffer.data(gi + 209 * ncomps + c);
        const auto *gi_210 = buffer.data(gi + 210 * ncomps + c);
        const auto *gi_211 = buffer.data(gi + 211 * ncomps + c);
        const auto *gi_212 = buffer.data(gi + 212 * ncomps + c);
        const auto *gi_213 = buffer.data(gi + 213 * ncomps + c);
        const auto *gi_214 = buffer.data(gi + 214 * ncomps + c);
        const auto *gi_215 = buffer.data(gi + 215 * ncomps + c);
        const auto *gi_216 = buffer.data(gi + 216 * ncomps + c);
        const auto *gi_224 = buffer.data(gi + 224 * ncomps + c);
        const auto *gi_225 = buffer.data(gi + 225 * ncomps + c);
        const auto *gi_226 = buffer.data(gi + 226 * ncomps + c);
        const auto *gi_227 = buffer.data(gi + 227 * ncomps + c);
        const auto *gi_228 = buffer.data(gi + 228 * ncomps + c);
        const auto *gi_229 = buffer.data(gi + 229 * ncomps + c);
        const auto *gi_230 = buffer.data(gi + 230 * ncomps + c);
        const auto *gi_231 = buffer.data(gi + 231 * ncomps + c);
        const auto *gi_232 = buffer.data(gi + 232 * ncomps + c);
        const auto *gi_233 = buffer.data(gi + 233 * ncomps + c);
        const auto *gi_234 = buffer.data(gi + 234 * ncomps + c);
        const auto *gi_235 = buffer.data(gi + 235 * ncomps + c);
        const auto *gi_236 = buffer.data(gi + 236 * ncomps + c);
        const auto *gi_237 = buffer.data(gi + 237 * ncomps + c);
        const auto *gi_238 = buffer.data(gi + 238 * ncomps + c);
        const auto *gi_239 = buffer.data(gi + 239 * ncomps + c);
        const auto *gi_240 = buffer.data(gi + 240 * ncomps + c);
        const auto *gi_241 = buffer.data(gi + 241 * ncomps + c);
        const auto *gi_242 = buffer.data(gi + 242 * ncomps + c);
        const auto *gi_243 = buffer.data(gi + 243 * ncomps + c);
        const auto *gi_244 = buffer.data(gi + 244 * ncomps + c);
        const auto *gi_252 = buffer.data(gi + 252 * ncomps + c);
        const auto *gi_253 = buffer.data(gi + 253 * ncomps + c);
        const auto *gi_254 = buffer.data(gi + 254 * ncomps + c);
        const auto *gi_255 = buffer.data(gi + 255 * ncomps + c);
        const auto *gi_256 = buffer.data(gi + 256 * ncomps + c);
        const auto *gi_257 = buffer.data(gi + 257 * ncomps + c);
        const auto *gi_258 = buffer.data(gi + 258 * ncomps + c);
        const auto *gi_259 = buffer.data(gi + 259 * ncomps + c);
        const auto *gi_260 = buffer.data(gi + 260 * ncomps + c);
        const auto *gi_261 = buffer.data(gi + 261 * ncomps + c);
        const auto *gi_262 = buffer.data(gi + 262 * ncomps + c);
        const auto *gi_263 = buffer.data(gi + 263 * ncomps + c);
        const auto *gi_264 = buffer.data(gi + 264 * ncomps + c);
        const auto *gi_265 = buffer.data(gi + 265 * ncomps + c);
        const auto *gi_266 = buffer.data(gi + 266 * ncomps + c);
        const auto *gi_267 = buffer.data(gi + 267 * ncomps + c);
        const auto *gi_268 = buffer.data(gi + 268 * ncomps + c);
        const auto *gi_269 = buffer.data(gi + 269 * ncomps + c);
        const auto *gi_270 = buffer.data(gi + 270 * ncomps + c);
        const auto *gi_271 = buffer.data(gi + 271 * ncomps + c);
        const auto *gi_272 = buffer.data(gi + 272 * ncomps + c);
        const auto *gi_280 = buffer.data(gi + 280 * ncomps + c);
        const auto *gi_281 = buffer.data(gi + 281 * ncomps + c);
        const auto *gi_282 = buffer.data(gi + 282 * ncomps + c);
        const auto *gi_283 = buffer.data(gi + 283 * ncomps + c);
        const auto *gi_284 = buffer.data(gi + 284 * ncomps + c);
        const auto *gi_285 = buffer.data(gi + 285 * ncomps + c);
        const auto *gi_286 = buffer.data(gi + 286 * ncomps + c);
        const auto *gi_287 = buffer.data(gi + 287 * ncomps + c);
        const auto *gi_288 = buffer.data(gi + 288 * ncomps + c);
        const auto *gi_289 = buffer.data(gi + 289 * ncomps + c);
        const auto *gi_290 = buffer.data(gi + 290 * ncomps + c);
        const auto *gi_291 = buffer.data(gi + 291 * ncomps + c);
        const auto *gi_292 = buffer.data(gi + 292 * ncomps + c);
        const auto *gi_293 = buffer.data(gi + 293 * ncomps + c);
        const auto *gi_294 = buffer.data(gi + 294 * ncomps + c);
        const auto *gi_295 = buffer.data(gi + 295 * ncomps + c);
        const auto *gi_296 = buffer.data(gi + 296 * ncomps + c);
        const auto *gi_297 = buffer.data(gi + 297 * ncomps + c);
        const auto *gi_298 = buffer.data(gi + 298 * ncomps + c);
        const auto *gi_299 = buffer.data(gi + 299 * ncomps + c);
        const auto *gi_300 = buffer.data(gi + 300 * ncomps + c);
        const auto *gi_308 = buffer.data(gi + 308 * ncomps + c);
        const auto *gi_309 = buffer.data(gi + 309 * ncomps + c);
        const auto *gi_310 = buffer.data(gi + 310 * ncomps + c);
        const auto *gi_311 = buffer.data(gi + 311 * ncomps + c);
        const auto *gi_312 = buffer.data(gi + 312 * ncomps + c);
        const auto *gi_313 = buffer.data(gi + 313 * ncomps + c);
        const auto *gi_314 = buffer.data(gi + 314 * ncomps + c);
        const auto *gi_315 = buffer.data(gi + 315 * ncomps + c);
        const auto *gi_316 = buffer.data(gi + 316 * ncomps + c);
        const auto *gi_317 = buffer.data(gi + 317 * ncomps + c);
        const auto *gi_318 = buffer.data(gi + 318 * ncomps + c);
        const auto *gi_319 = buffer.data(gi + 319 * ncomps + c);
        const auto *gi_320 = buffer.data(gi + 320 * ncomps + c);
        const auto *gi_321 = buffer.data(gi + 321 * ncomps + c);
        const auto *gi_322 = buffer.data(gi + 322 * ncomps + c);
        const auto *gi_323 = buffer.data(gi + 323 * ncomps + c);
        const auto *gi_324 = buffer.data(gi + 324 * ncomps + c);
        const auto *gi_325 = buffer.data(gi + 325 * ncomps + c);
        const auto *gi_326 = buffer.data(gi + 326 * ncomps + c);
        const auto *gi_327 = buffer.data(gi + 327 * ncomps + c);
        const auto *gi_328 = buffer.data(gi + 328 * ncomps + c);
        const auto *gi_336 = buffer.data(gi + 336 * ncomps + c);
        const auto *gi_337 = buffer.data(gi + 337 * ncomps + c);
        const auto *gi_338 = buffer.data(gi + 338 * ncomps + c);
        const auto *gi_339 = buffer.data(gi + 339 * ncomps + c);
        const auto *gi_340 = buffer.data(gi + 340 * ncomps + c);
        const auto *gi_341 = buffer.data(gi + 341 * ncomps + c);
        const auto *gi_342 = buffer.data(gi + 342 * ncomps + c);
        const auto *gi_343 = buffer.data(gi + 343 * ncomps + c);
        const auto *gi_344 = buffer.data(gi + 344 * ncomps + c);
        const auto *gi_345 = buffer.data(gi + 345 * ncomps + c);
        const auto *gi_346 = buffer.data(gi + 346 * ncomps + c);
        const auto *gi_347 = buffer.data(gi + 347 * ncomps + c);
        const auto *gi_348 = buffer.data(gi + 348 * ncomps + c);
        const auto *gi_349 = buffer.data(gi + 349 * ncomps + c);
        const auto *gi_350 = buffer.data(gi + 350 * ncomps + c);
        const auto *gi_351 = buffer.data(gi + 351 * ncomps + c);
        const auto *gi_352 = buffer.data(gi + 352 * ncomps + c);
        const auto *gi_353 = buffer.data(gi + 353 * ncomps + c);
        const auto *gi_354 = buffer.data(gi + 354 * ncomps + c);
        const auto *gi_355 = buffer.data(gi + 355 * ncomps + c);
        const auto *gi_356 = buffer.data(gi + 356 * ncomps + c);
        const auto *gi_364 = buffer.data(gi + 364 * ncomps + c);
        const auto *gi_365 = buffer.data(gi + 365 * ncomps + c);
        const auto *gi_366 = buffer.data(gi + 366 * ncomps + c);
        const auto *gi_367 = buffer.data(gi + 367 * ncomps + c);
        const auto *gi_368 = buffer.data(gi + 368 * ncomps + c);
        const auto *gi_369 = buffer.data(gi + 369 * ncomps + c);
        const auto *gi_370 = buffer.data(gi + 370 * ncomps + c);
        const auto *gi_371 = buffer.data(gi + 371 * ncomps + c);
        const auto *gi_372 = buffer.data(gi + 372 * ncomps + c);
        const auto *gi_373 = buffer.data(gi + 373 * ncomps + c);
        const auto *gi_374 = buffer.data(gi + 374 * ncomps + c);
        const auto *gi_375 = buffer.data(gi + 375 * ncomps + c);
        const auto *gi_376 = buffer.data(gi + 376 * ncomps + c);
        const auto *gi_377 = buffer.data(gi + 377 * ncomps + c);
        const auto *gi_378 = buffer.data(gi + 378 * ncomps + c);
        const auto *gi_379 = buffer.data(gi + 379 * ncomps + c);
        const auto *gi_380 = buffer.data(gi + 380 * ncomps + c);

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, gh_145, gh_146, gh_147, \
                         gh_148, gh_149, gi_187, gi_188, gi_196, gi_197, \
                         gi_198 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_145[k] = -ab_x[k] * gh_145[k]
                       + gi_187[k];

            t_146[k] = -ab_x[k] * gh_146[k]
                       + gi_188[k];

            t_147[k] = -ab_x[k] * gh_147[k]
                       + gi_196[k];

            t_148[k] = -ab_x[k] * gh_148[k]
                       + gi_197[k];

            t_149[k] = -ab_x[k] * gh_149[k]
                       + gi_198[k];
        }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, gh_150, gh_151, gh_152, \
                         gh_153, gh_154, gi_199, gi_200, gi_201, gi_202, \
                         gi_203 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_150[k] = -ab_x[k] * gh_150[k]
                       + gi_199[k];

            t_151[k] = -ab_x[k] * gh_151[k]
                       + gi_200[k];

            t_152[k] = -ab_x[k] * gh_152[k]
                       + gi_201[k];

            t_153[k] = -ab_x[k] * gh_153[k]
                       + gi_202[k];

            t_154[k] = -ab_x[k] * gh_154[k]
                       + gi_203[k];
        }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, gh_155, gh_156, gh_157, \
                         gh_158, gh_159, gi_204, gi_205, gi_206, gi_207, \
                         gi_208 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_155[k] = -ab_x[k] * gh_155[k]
                       + gi_204[k];

            t_156[k] = -ab_x[k] * gh_156[k]
                       + gi_205[k];

            t_157[k] = -ab_x[k] * gh_157[k]
                       + gi_206[k];

            t_158[k] = -ab_x[k] * gh_158[k]
                       + gi_207[k];

            t_159[k] = -ab_x[k] * gh_159[k]
                       + gi_208[k];
        }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, gh_160, gh_161, gh_162, \
                         gh_163, gh_164, gi_209, gi_210, gi_211, gi_212, \
                         gi_213 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_160[k] = -ab_x[k] * gh_160[k]
                       + gi_209[k];

            t_161[k] = -ab_x[k] * gh_161[k]
                       + gi_210[k];

            t_162[k] = -ab_x[k] * gh_162[k]
                       + gi_211[k];

            t_163[k] = -ab_x[k] * gh_163[k]
                       + gi_212[k];

            t_164[k] = -ab_x[k] * gh_164[k]
                       + gi_213[k];
        }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, gh_165, gh_166, gh_167, \
                         gh_168, gh_169, gi_214, gi_215, gi_216, gi_224, \
                         gi_225 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_165[k] = -ab_x[k] * gh_165[k]
                       + gi_214[k];

            t_166[k] = -ab_x[k] * gh_166[k]
                       + gi_215[k];

            t_167[k] = -ab_x[k] * gh_167[k]
                       + gi_216[k];

            t_168[k] = -ab_x[k] * gh_168[k]
                       + gi_224[k];

            t_169[k] = -ab_x[k] * gh_169[k]
                       + gi_225[k];
        }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, gh_170, gh_171, gh_172, \
                         gh_173, gh_174, gi_226, gi_227, gi_228, gi_229, \
                         gi_230 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_170[k] = -ab_x[k] * gh_170[k]
                       + gi_226[k];

            t_171[k] = -ab_x[k] * gh_171[k]
                       + gi_227[k];

            t_172[k] = -ab_x[k] * gh_172[k]
                       + gi_228[k];

            t_173[k] = -ab_x[k] * gh_173[k]
                       + gi_229[k];

            t_174[k] = -ab_x[k] * gh_174[k]
                       + gi_230[k];
        }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, gh_175, gh_176, gh_177, \
                         gh_178, gh_179, gi_231, gi_232, gi_233, gi_234, \
                         gi_235 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_175[k] = -ab_x[k] * gh_175[k]
                       + gi_231[k];

            t_176[k] = -ab_x[k] * gh_176[k]
                       + gi_232[k];

            t_177[k] = -ab_x[k] * gh_177[k]
                       + gi_233[k];

            t_178[k] = -ab_x[k] * gh_178[k]
                       + gi_234[k];

            t_179[k] = -ab_x[k] * gh_179[k]
                       + gi_235[k];
        }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, gh_180, gh_181, gh_182, \
                         gh_183, gh_184, gi_236, gi_237, gi_238, gi_239, \
                         gi_240 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_180[k] = -ab_x[k] * gh_180[k]
                       + gi_236[k];

            t_181[k] = -ab_x[k] * gh_181[k]
                       + gi_237[k];

            t_182[k] = -ab_x[k] * gh_182[k]
                       + gi_238[k];

            t_183[k] = -ab_x[k] * gh_183[k]
                       + gi_239[k];

            t_184[k] = -ab_x[k] * gh_184[k]
                       + gi_240[k];
        }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, gh_185, gh_186, gh_187, \
                         gh_188, gh_189, gi_241, gi_242, gi_243, gi_244, \
                         gi_252 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_185[k] = -ab_x[k] * gh_185[k]
                       + gi_241[k];

            t_186[k] = -ab_x[k] * gh_186[k]
                       + gi_242[k];

            t_187[k] = -ab_x[k] * gh_187[k]
                       + gi_243[k];

            t_188[k] = -ab_x[k] * gh_188[k]
                       + gi_244[k];

            t_189[k] = -ab_x[k] * gh_189[k]
                       + gi_252[k];
        }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, gh_190, gh_191, gh_192, \
                         gh_193, gh_194, gi_253, gi_254, gi_255, gi_256, \
                         gi_257 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_190[k] = -ab_x[k] * gh_190[k]
                       + gi_253[k];

            t_191[k] = -ab_x[k] * gh_191[k]
                       + gi_254[k];

            t_192[k] = -ab_x[k] * gh_192[k]
                       + gi_255[k];

            t_193[k] = -ab_x[k] * gh_193[k]
                       + gi_256[k];

            t_194[k] = -ab_x[k] * gh_194[k]
                       + gi_257[k];
        }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, gh_195, gh_196, gh_197, \
                         gh_198, gh_199, gi_258, gi_259, gi_260, gi_261, \
                         gi_262 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_195[k] = -ab_x[k] * gh_195[k]
                       + gi_258[k];

            t_196[k] = -ab_x[k] * gh_196[k]
                       + gi_259[k];

            t_197[k] = -ab_x[k] * gh_197[k]
                       + gi_260[k];

            t_198[k] = -ab_x[k] * gh_198[k]
                       + gi_261[k];

            t_199[k] = -ab_x[k] * gh_199[k]
                       + gi_262[k];
        }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, gh_200, gh_201, gh_202, \
                         gh_203, gh_204, gi_263, gi_264, gi_265, gi_266, \
                         gi_267 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_200[k] = -ab_x[k] * gh_200[k]
                       + gi_263[k];

            t_201[k] = -ab_x[k] * gh_201[k]
                       + gi_264[k];

            t_202[k] = -ab_x[k] * gh_202[k]
                       + gi_265[k];

            t_203[k] = -ab_x[k] * gh_203[k]
                       + gi_266[k];

            t_204[k] = -ab_x[k] * gh_204[k]
                       + gi_267[k];
        }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, gh_205, gh_206, gh_207, \
                         gh_208, gh_209, gi_268, gi_269, gi_270, gi_271, \
                         gi_272 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_205[k] = -ab_x[k] * gh_205[k]
                       + gi_268[k];

            t_206[k] = -ab_x[k] * gh_206[k]
                       + gi_269[k];

            t_207[k] = -ab_x[k] * gh_207[k]
                       + gi_270[k];

            t_208[k] = -ab_x[k] * gh_208[k]
                       + gi_271[k];

            t_209[k] = -ab_x[k] * gh_209[k]
                       + gi_272[k];
        }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, gh_210, gh_211, gh_212, \
                         gh_213, gh_214, gi_280, gi_281, gi_282, gi_283, \
                         gi_284 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_210[k] = -ab_x[k] * gh_210[k]
                       + gi_280[k];

            t_211[k] = -ab_x[k] * gh_211[k]
                       + gi_281[k];

            t_212[k] = -ab_x[k] * gh_212[k]
                       + gi_282[k];

            t_213[k] = -ab_x[k] * gh_213[k]
                       + gi_283[k];

            t_214[k] = -ab_x[k] * gh_214[k]
                       + gi_284[k];
        }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, gh_215, gh_216, gh_217, \
                         gh_218, gh_219, gi_285, gi_286, gi_287, gi_288, \
                         gi_289 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_215[k] = -ab_x[k] * gh_215[k]
                       + gi_285[k];

            t_216[k] = -ab_x[k] * gh_216[k]
                       + gi_286[k];

            t_217[k] = -ab_x[k] * gh_217[k]
                       + gi_287[k];

            t_218[k] = -ab_x[k] * gh_218[k]
                       + gi_288[k];

            t_219[k] = -ab_x[k] * gh_219[k]
                       + gi_289[k];
        }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, gh_220, gh_221, gh_222, \
                         gh_223, gh_224, gi_290, gi_291, gi_292, gi_293, \
                         gi_294 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_220[k] = -ab_x[k] * gh_220[k]
                       + gi_290[k];

            t_221[k] = -ab_x[k] * gh_221[k]
                       + gi_291[k];

            t_222[k] = -ab_x[k] * gh_222[k]
                       + gi_292[k];

            t_223[k] = -ab_x[k] * gh_223[k]
                       + gi_293[k];

            t_224[k] = -ab_x[k] * gh_224[k]
                       + gi_294[k];
        }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, gh_225, gh_226, gh_227, \
                         gh_228, gh_229, gi_295, gi_296, gi_297, gi_298, \
                         gi_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_225[k] = -ab_x[k] * gh_225[k]
                       + gi_295[k];

            t_226[k] = -ab_x[k] * gh_226[k]
                       + gi_296[k];

            t_227[k] = -ab_x[k] * gh_227[k]
                       + gi_297[k];

            t_228[k] = -ab_x[k] * gh_228[k]
                       + gi_298[k];

            t_229[k] = -ab_x[k] * gh_229[k]
                       + gi_299[k];
        }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, gh_230, gh_231, gh_232, \
                         gh_233, gh_234, gi_300, gi_308, gi_309, gi_310, \
                         gi_311 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_230[k] = -ab_x[k] * gh_230[k]
                       + gi_300[k];

            t_231[k] = -ab_x[k] * gh_231[k]
                       + gi_308[k];

            t_232[k] = -ab_x[k] * gh_232[k]
                       + gi_309[k];

            t_233[k] = -ab_x[k] * gh_233[k]
                       + gi_310[k];

            t_234[k] = -ab_x[k] * gh_234[k]
                       + gi_311[k];
        }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_x, gh_235, gh_236, gh_237, \
                         gh_238, gh_239, gi_312, gi_313, gi_314, gi_315, \
                         gi_316 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_235[k] = -ab_x[k] * gh_235[k]
                       + gi_312[k];

            t_236[k] = -ab_x[k] * gh_236[k]
                       + gi_313[k];

            t_237[k] = -ab_x[k] * gh_237[k]
                       + gi_314[k];

            t_238[k] = -ab_x[k] * gh_238[k]
                       + gi_315[k];

            t_239[k] = -ab_x[k] * gh_239[k]
                       + gi_316[k];
        }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, gh_240, gh_241, gh_242, \
                         gh_243, gh_244, gi_317, gi_318, gi_319, gi_320, \
                         gi_321 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_240[k] = -ab_x[k] * gh_240[k]
                       + gi_317[k];

            t_241[k] = -ab_x[k] * gh_241[k]
                       + gi_318[k];

            t_242[k] = -ab_x[k] * gh_242[k]
                       + gi_319[k];

            t_243[k] = -ab_x[k] * gh_243[k]
                       + gi_320[k];

            t_244[k] = -ab_x[k] * gh_244[k]
                       + gi_321[k];
        }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, gh_245, gh_246, gh_247, \
                         gh_248, gh_249, gi_322, gi_323, gi_324, gi_325, \
                         gi_326 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_245[k] = -ab_x[k] * gh_245[k]
                       + gi_322[k];

            t_246[k] = -ab_x[k] * gh_246[k]
                       + gi_323[k];

            t_247[k] = -ab_x[k] * gh_247[k]
                       + gi_324[k];

            t_248[k] = -ab_x[k] * gh_248[k]
                       + gi_325[k];

            t_249[k] = -ab_x[k] * gh_249[k]
                       + gi_326[k];
        }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_x, gh_250, gh_251, gh_252, \
                         gh_253, gh_254, gi_327, gi_328, gi_336, gi_337, \
                         gi_338 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_250[k] = -ab_x[k] * gh_250[k]
                       + gi_327[k];

            t_251[k] = -ab_x[k] * gh_251[k]
                       + gi_328[k];

            t_252[k] = -ab_x[k] * gh_252[k]
                       + gi_336[k];

            t_253[k] = -ab_x[k] * gh_253[k]
                       + gi_337[k];

            t_254[k] = -ab_x[k] * gh_254[k]
                       + gi_338[k];
        }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, gh_255, gh_256, gh_257, \
                         gh_258, gh_259, gi_339, gi_340, gi_341, gi_342, \
                         gi_343 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_255[k] = -ab_x[k] * gh_255[k]
                       + gi_339[k];

            t_256[k] = -ab_x[k] * gh_256[k]
                       + gi_340[k];

            t_257[k] = -ab_x[k] * gh_257[k]
                       + gi_341[k];

            t_258[k] = -ab_x[k] * gh_258[k]
                       + gi_342[k];

            t_259[k] = -ab_x[k] * gh_259[k]
                       + gi_343[k];
        }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, gh_260, gh_261, gh_262, \
                         gh_263, gh_264, gi_344, gi_345, gi_346, gi_347, \
                         gi_348 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_260[k] = -ab_x[k] * gh_260[k]
                       + gi_344[k];

            t_261[k] = -ab_x[k] * gh_261[k]
                       + gi_345[k];

            t_262[k] = -ab_x[k] * gh_262[k]
                       + gi_346[k];

            t_263[k] = -ab_x[k] * gh_263[k]
                       + gi_347[k];

            t_264[k] = -ab_x[k] * gh_264[k]
                       + gi_348[k];
        }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, gh_265, gh_266, gh_267, \
                         gh_268, gh_269, gi_349, gi_350, gi_351, gi_352, \
                         gi_353 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_265[k] = -ab_x[k] * gh_265[k]
                       + gi_349[k];

            t_266[k] = -ab_x[k] * gh_266[k]
                       + gi_350[k];

            t_267[k] = -ab_x[k] * gh_267[k]
                       + gi_351[k];

            t_268[k] = -ab_x[k] * gh_268[k]
                       + gi_352[k];

            t_269[k] = -ab_x[k] * gh_269[k]
                       + gi_353[k];
        }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, gh_270, gh_271, gh_272, \
                         gh_273, gh_274, gi_354, gi_355, gi_356, gi_364, \
                         gi_365 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_270[k] = -ab_x[k] * gh_270[k]
                       + gi_354[k];

            t_271[k] = -ab_x[k] * gh_271[k]
                       + gi_355[k];

            t_272[k] = -ab_x[k] * gh_272[k]
                       + gi_356[k];

            t_273[k] = -ab_x[k] * gh_273[k]
                       + gi_364[k];

            t_274[k] = -ab_x[k] * gh_274[k]
                       + gi_365[k];
        }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, gh_275, gh_276, gh_277, \
                         gh_278, gh_279, gi_366, gi_367, gi_368, gi_369, \
                         gi_370 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_275[k] = -ab_x[k] * gh_275[k]
                       + gi_366[k];

            t_276[k] = -ab_x[k] * gh_276[k]
                       + gi_367[k];

            t_277[k] = -ab_x[k] * gh_277[k]
                       + gi_368[k];

            t_278[k] = -ab_x[k] * gh_278[k]
                       + gi_369[k];

            t_279[k] = -ab_x[k] * gh_279[k]
                       + gi_370[k];
        }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_x, gh_280, gh_281, gh_282, \
                         gh_283, gh_284, gi_371, gi_372, gi_373, gi_374, \
                         gi_375 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_280[k] = -ab_x[k] * gh_280[k]
                       + gi_371[k];

            t_281[k] = -ab_x[k] * gh_281[k]
                       + gi_372[k];

            t_282[k] = -ab_x[k] * gh_282[k]
                       + gi_373[k];

            t_283[k] = -ab_x[k] * gh_283[k]
                       + gi_374[k];

            t_284[k] = -ab_x[k] * gh_284[k]
                       + gi_375[k];
        }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, gh_285, gh_286, gh_287, \
                         gh_288, gh_289, gi_376, gi_377, gi_378, gi_379, \
                         gi_380 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_285[k] = -ab_x[k] * gh_285[k]
                       + gi_376[k];

            t_286[k] = -ab_x[k] * gh_286[k]
                       + gi_377[k];

            t_287[k] = -ab_x[k] * gh_287[k]
                       + gi_378[k];

            t_288[k] = -ab_x[k] * gh_288[k]
                       + gi_379[k];

            t_289[k] = -ab_x[k] * gh_289[k]
                       + gi_380[k];
        }
    }
}

static auto
compute_hrr_hh_piece2(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t gh, const size_t gi, const size_t ncomps,
                      const size_t nmax) -> void
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

        const auto *gh_210 = buffer.data(gh + 210 * ncomps + c);
        const auto *gh_211 = buffer.data(gh + 211 * ncomps + c);
        const auto *gh_212 = buffer.data(gh + 212 * ncomps + c);
        const auto *gh_213 = buffer.data(gh + 213 * ncomps + c);
        const auto *gh_214 = buffer.data(gh + 214 * ncomps + c);
        const auto *gh_215 = buffer.data(gh + 215 * ncomps + c);
        const auto *gh_216 = buffer.data(gh + 216 * ncomps + c);
        const auto *gh_217 = buffer.data(gh + 217 * ncomps + c);
        const auto *gh_218 = buffer.data(gh + 218 * ncomps + c);
        const auto *gh_219 = buffer.data(gh + 219 * ncomps + c);
        const auto *gh_220 = buffer.data(gh + 220 * ncomps + c);
        const auto *gh_221 = buffer.data(gh + 221 * ncomps + c);
        const auto *gh_222 = buffer.data(gh + 222 * ncomps + c);
        const auto *gh_223 = buffer.data(gh + 223 * ncomps + c);
        const auto *gh_224 = buffer.data(gh + 224 * ncomps + c);
        const auto *gh_225 = buffer.data(gh + 225 * ncomps + c);
        const auto *gh_226 = buffer.data(gh + 226 * ncomps + c);
        const auto *gh_227 = buffer.data(gh + 227 * ncomps + c);
        const auto *gh_228 = buffer.data(gh + 228 * ncomps + c);
        const auto *gh_229 = buffer.data(gh + 229 * ncomps + c);
        const auto *gh_230 = buffer.data(gh + 230 * ncomps + c);
        const auto *gh_231 = buffer.data(gh + 231 * ncomps + c);
        const auto *gh_232 = buffer.data(gh + 232 * ncomps + c);
        const auto *gh_233 = buffer.data(gh + 233 * ncomps + c);
        const auto *gh_234 = buffer.data(gh + 234 * ncomps + c);
        const auto *gh_235 = buffer.data(gh + 235 * ncomps + c);
        const auto *gh_236 = buffer.data(gh + 236 * ncomps + c);
        const auto *gh_237 = buffer.data(gh + 237 * ncomps + c);
        const auto *gh_238 = buffer.data(gh + 238 * ncomps + c);
        const auto *gh_239 = buffer.data(gh + 239 * ncomps + c);
        const auto *gh_240 = buffer.data(gh + 240 * ncomps + c);
        const auto *gh_241 = buffer.data(gh + 241 * ncomps + c);
        const auto *gh_242 = buffer.data(gh + 242 * ncomps + c);
        const auto *gh_243 = buffer.data(gh + 243 * ncomps + c);
        const auto *gh_244 = buffer.data(gh + 244 * ncomps + c);
        const auto *gh_245 = buffer.data(gh + 245 * ncomps + c);
        const auto *gh_246 = buffer.data(gh + 246 * ncomps + c);
        const auto *gh_247 = buffer.data(gh + 247 * ncomps + c);
        const auto *gh_248 = buffer.data(gh + 248 * ncomps + c);
        const auto *gh_249 = buffer.data(gh + 249 * ncomps + c);
        const auto *gh_250 = buffer.data(gh + 250 * ncomps + c);
        const auto *gh_251 = buffer.data(gh + 251 * ncomps + c);
        const auto *gh_252 = buffer.data(gh + 252 * ncomps + c);
        const auto *gh_253 = buffer.data(gh + 253 * ncomps + c);
        const auto *gh_254 = buffer.data(gh + 254 * ncomps + c);
        const auto *gh_255 = buffer.data(gh + 255 * ncomps + c);
        const auto *gh_256 = buffer.data(gh + 256 * ncomps + c);
        const auto *gh_257 = buffer.data(gh + 257 * ncomps + c);
        const auto *gh_258 = buffer.data(gh + 258 * ncomps + c);
        const auto *gh_259 = buffer.data(gh + 259 * ncomps + c);
        const auto *gh_260 = buffer.data(gh + 260 * ncomps + c);
        const auto *gh_261 = buffer.data(gh + 261 * ncomps + c);
        const auto *gh_262 = buffer.data(gh + 262 * ncomps + c);
        const auto *gh_263 = buffer.data(gh + 263 * ncomps + c);
        const auto *gh_264 = buffer.data(gh + 264 * ncomps + c);
        const auto *gh_265 = buffer.data(gh + 265 * ncomps + c);
        const auto *gh_266 = buffer.data(gh + 266 * ncomps + c);
        const auto *gh_267 = buffer.data(gh + 267 * ncomps + c);
        const auto *gh_268 = buffer.data(gh + 268 * ncomps + c);
        const auto *gh_269 = buffer.data(gh + 269 * ncomps + c);
        const auto *gh_270 = buffer.data(gh + 270 * ncomps + c);
        const auto *gh_271 = buffer.data(gh + 271 * ncomps + c);
        const auto *gh_272 = buffer.data(gh + 272 * ncomps + c);
        const auto *gh_273 = buffer.data(gh + 273 * ncomps + c);
        const auto *gh_274 = buffer.data(gh + 274 * ncomps + c);
        const auto *gh_275 = buffer.data(gh + 275 * ncomps + c);
        const auto *gh_276 = buffer.data(gh + 276 * ncomps + c);
        const auto *gh_277 = buffer.data(gh + 277 * ncomps + c);
        const auto *gh_278 = buffer.data(gh + 278 * ncomps + c);
        const auto *gh_279 = buffer.data(gh + 279 * ncomps + c);
        const auto *gh_280 = buffer.data(gh + 280 * ncomps + c);
        const auto *gh_281 = buffer.data(gh + 281 * ncomps + c);
        const auto *gh_282 = buffer.data(gh + 282 * ncomps + c);
        const auto *gh_283 = buffer.data(gh + 283 * ncomps + c);
        const auto *gh_284 = buffer.data(gh + 284 * ncomps + c);
        const auto *gh_285 = buffer.data(gh + 285 * ncomps + c);
        const auto *gh_286 = buffer.data(gh + 286 * ncomps + c);
        const auto *gh_287 = buffer.data(gh + 287 * ncomps + c);
        const auto *gh_288 = buffer.data(gh + 288 * ncomps + c);
        const auto *gh_289 = buffer.data(gh + 289 * ncomps + c);
        const auto *gh_290 = buffer.data(gh + 290 * ncomps + c);
        const auto *gh_291 = buffer.data(gh + 291 * ncomps + c);
        const auto *gh_292 = buffer.data(gh + 292 * ncomps + c);
        const auto *gh_293 = buffer.data(gh + 293 * ncomps + c);
        const auto *gh_294 = buffer.data(gh + 294 * ncomps + c);
        const auto *gh_295 = buffer.data(gh + 295 * ncomps + c);
        const auto *gh_296 = buffer.data(gh + 296 * ncomps + c);
        const auto *gh_297 = buffer.data(gh + 297 * ncomps + c);
        const auto *gh_298 = buffer.data(gh + 298 * ncomps + c);
        const auto *gh_299 = buffer.data(gh + 299 * ncomps + c);
        const auto *gh_300 = buffer.data(gh + 300 * ncomps + c);
        const auto *gh_301 = buffer.data(gh + 301 * ncomps + c);
        const auto *gh_302 = buffer.data(gh + 302 * ncomps + c);
        const auto *gh_303 = buffer.data(gh + 303 * ncomps + c);
        const auto *gh_304 = buffer.data(gh + 304 * ncomps + c);
        const auto *gh_305 = buffer.data(gh + 305 * ncomps + c);
        const auto *gh_306 = buffer.data(gh + 306 * ncomps + c);
        const auto *gh_307 = buffer.data(gh + 307 * ncomps + c);
        const auto *gh_308 = buffer.data(gh + 308 * ncomps + c);
        const auto *gh_309 = buffer.data(gh + 309 * ncomps + c);
        const auto *gh_310 = buffer.data(gh + 310 * ncomps + c);
        const auto *gh_311 = buffer.data(gh + 311 * ncomps + c);
        const auto *gh_312 = buffer.data(gh + 312 * ncomps + c);
        const auto *gh_313 = buffer.data(gh + 313 * ncomps + c);
        const auto *gh_314 = buffer.data(gh + 314 * ncomps + c);

        const auto *gi_281 = buffer.data(gi + 281 * ncomps + c);
        const auto *gi_283 = buffer.data(gi + 283 * ncomps + c);
        const auto *gi_284 = buffer.data(gi + 284 * ncomps + c);
        const auto *gi_286 = buffer.data(gi + 286 * ncomps + c);
        const auto *gi_287 = buffer.data(gi + 287 * ncomps + c);
        const auto *gi_288 = buffer.data(gi + 288 * ncomps + c);
        const auto *gi_290 = buffer.data(gi + 290 * ncomps + c);
        const auto *gi_291 = buffer.data(gi + 291 * ncomps + c);
        const auto *gi_292 = buffer.data(gi + 292 * ncomps + c);
        const auto *gi_293 = buffer.data(gi + 293 * ncomps + c);
        const auto *gi_295 = buffer.data(gi + 295 * ncomps + c);
        const auto *gi_296 = buffer.data(gi + 296 * ncomps + c);
        const auto *gi_297 = buffer.data(gi + 297 * ncomps + c);
        const auto *gi_298 = buffer.data(gi + 298 * ncomps + c);
        const auto *gi_299 = buffer.data(gi + 299 * ncomps + c);
        const auto *gi_301 = buffer.data(gi + 301 * ncomps + c);
        const auto *gi_302 = buffer.data(gi + 302 * ncomps + c);
        const auto *gi_303 = buffer.data(gi + 303 * ncomps + c);
        const auto *gi_304 = buffer.data(gi + 304 * ncomps + c);
        const auto *gi_305 = buffer.data(gi + 305 * ncomps + c);
        const auto *gi_306 = buffer.data(gi + 306 * ncomps + c);
        const auto *gi_309 = buffer.data(gi + 309 * ncomps + c);
        const auto *gi_311 = buffer.data(gi + 311 * ncomps + c);
        const auto *gi_312 = buffer.data(gi + 312 * ncomps + c);
        const auto *gi_314 = buffer.data(gi + 314 * ncomps + c);
        const auto *gi_315 = buffer.data(gi + 315 * ncomps + c);
        const auto *gi_316 = buffer.data(gi + 316 * ncomps + c);
        const auto *gi_318 = buffer.data(gi + 318 * ncomps + c);
        const auto *gi_319 = buffer.data(gi + 319 * ncomps + c);
        const auto *gi_320 = buffer.data(gi + 320 * ncomps + c);
        const auto *gi_321 = buffer.data(gi + 321 * ncomps + c);
        const auto *gi_323 = buffer.data(gi + 323 * ncomps + c);
        const auto *gi_324 = buffer.data(gi + 324 * ncomps + c);
        const auto *gi_325 = buffer.data(gi + 325 * ncomps + c);
        const auto *gi_326 = buffer.data(gi + 326 * ncomps + c);
        const auto *gi_327 = buffer.data(gi + 327 * ncomps + c);
        const auto *gi_329 = buffer.data(gi + 329 * ncomps + c);
        const auto *gi_330 = buffer.data(gi + 330 * ncomps + c);
        const auto *gi_331 = buffer.data(gi + 331 * ncomps + c);
        const auto *gi_332 = buffer.data(gi + 332 * ncomps + c);
        const auto *gi_333 = buffer.data(gi + 333 * ncomps + c);
        const auto *gi_334 = buffer.data(gi + 334 * ncomps + c);
        const auto *gi_337 = buffer.data(gi + 337 * ncomps + c);
        const auto *gi_339 = buffer.data(gi + 339 * ncomps + c);
        const auto *gi_340 = buffer.data(gi + 340 * ncomps + c);
        const auto *gi_342 = buffer.data(gi + 342 * ncomps + c);
        const auto *gi_343 = buffer.data(gi + 343 * ncomps + c);
        const auto *gi_344 = buffer.data(gi + 344 * ncomps + c);
        const auto *gi_346 = buffer.data(gi + 346 * ncomps + c);
        const auto *gi_347 = buffer.data(gi + 347 * ncomps + c);
        const auto *gi_348 = buffer.data(gi + 348 * ncomps + c);
        const auto *gi_349 = buffer.data(gi + 349 * ncomps + c);
        const auto *gi_351 = buffer.data(gi + 351 * ncomps + c);
        const auto *gi_352 = buffer.data(gi + 352 * ncomps + c);
        const auto *gi_353 = buffer.data(gi + 353 * ncomps + c);
        const auto *gi_354 = buffer.data(gi + 354 * ncomps + c);
        const auto *gi_355 = buffer.data(gi + 355 * ncomps + c);
        const auto *gi_357 = buffer.data(gi + 357 * ncomps + c);
        const auto *gi_358 = buffer.data(gi + 358 * ncomps + c);
        const auto *gi_359 = buffer.data(gi + 359 * ncomps + c);
        const auto *gi_360 = buffer.data(gi + 360 * ncomps + c);
        const auto *gi_361 = buffer.data(gi + 361 * ncomps + c);
        const auto *gi_362 = buffer.data(gi + 362 * ncomps + c);
        const auto *gi_365 = buffer.data(gi + 365 * ncomps + c);
        const auto *gi_367 = buffer.data(gi + 367 * ncomps + c);
        const auto *gi_368 = buffer.data(gi + 368 * ncomps + c);
        const auto *gi_370 = buffer.data(gi + 370 * ncomps + c);
        const auto *gi_371 = buffer.data(gi + 371 * ncomps + c);
        const auto *gi_372 = buffer.data(gi + 372 * ncomps + c);
        const auto *gi_374 = buffer.data(gi + 374 * ncomps + c);
        const auto *gi_375 = buffer.data(gi + 375 * ncomps + c);
        const auto *gi_376 = buffer.data(gi + 376 * ncomps + c);
        const auto *gi_377 = buffer.data(gi + 377 * ncomps + c);
        const auto *gi_379 = buffer.data(gi + 379 * ncomps + c);
        const auto *gi_380 = buffer.data(gi + 380 * ncomps + c);
        const auto *gi_381 = buffer.data(gi + 381 * ncomps + c);
        const auto *gi_382 = buffer.data(gi + 382 * ncomps + c);
        const auto *gi_383 = buffer.data(gi + 383 * ncomps + c);
        const auto *gi_384 = buffer.data(gi + 384 * ncomps + c);
        const auto *gi_385 = buffer.data(gi + 385 * ncomps + c);
        const auto *gi_386 = buffer.data(gi + 386 * ncomps + c);
        const auto *gi_387 = buffer.data(gi + 387 * ncomps + c);
        const auto *gi_388 = buffer.data(gi + 388 * ncomps + c);
        const auto *gi_389 = buffer.data(gi + 389 * ncomps + c);
        const auto *gi_390 = buffer.data(gi + 390 * ncomps + c);
        const auto *gi_392 = buffer.data(gi + 392 * ncomps + c);
        const auto *gi_393 = buffer.data(gi + 393 * ncomps + c);
        const auto *gi_394 = buffer.data(gi + 394 * ncomps + c);
        const auto *gi_395 = buffer.data(gi + 395 * ncomps + c);
        const auto *gi_396 = buffer.data(gi + 396 * ncomps + c);
        const auto *gi_397 = buffer.data(gi + 397 * ncomps + c);
        const auto *gi_398 = buffer.data(gi + 398 * ncomps + c);
        const auto *gi_399 = buffer.data(gi + 399 * ncomps + c);
        const auto *gi_400 = buffer.data(gi + 400 * ncomps + c);
        const auto *gi_401 = buffer.data(gi + 401 * ncomps + c);
        const auto *gi_402 = buffer.data(gi + 402 * ncomps + c);
        const auto *gi_403 = buffer.data(gi + 403 * ncomps + c);
        const auto *gi_404 = buffer.data(gi + 404 * ncomps + c);
        const auto *gi_405 = buffer.data(gi + 405 * ncomps + c);
        const auto *gi_406 = buffer.data(gi + 406 * ncomps + c);
        const auto *gi_407 = buffer.data(gi + 407 * ncomps + c);
        const auto *gi_408 = buffer.data(gi + 408 * ncomps + c);
        const auto *gi_409 = buffer.data(gi + 409 * ncomps + c);
        const auto *gi_410 = buffer.data(gi + 410 * ncomps + c);
        const auto *gi_411 = buffer.data(gi + 411 * ncomps + c);
        const auto *gi_412 = buffer.data(gi + 412 * ncomps + c);
        const auto *gi_413 = buffer.data(gi + 413 * ncomps + c);
        const auto *gi_414 = buffer.data(gi + 414 * ncomps + c);
        const auto *gi_415 = buffer.data(gi + 415 * ncomps + c);
        const auto *gi_416 = buffer.data(gi + 416 * ncomps + c);
        const auto *gi_417 = buffer.data(gi + 417 * ncomps + c);
        const auto *gi_418 = buffer.data(gi + 418 * ncomps + c);

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, gh_290, gh_291, gh_292, \
                         gh_293, gh_294, gi_381, gi_382, gi_383, gi_384, \
                         gi_392 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_290[k] = -ab_x[k] * gh_290[k]
                       + gi_381[k];

            t_291[k] = -ab_x[k] * gh_291[k]
                       + gi_382[k];

            t_292[k] = -ab_x[k] * gh_292[k]
                       + gi_383[k];

            t_293[k] = -ab_x[k] * gh_293[k]
                       + gi_384[k];

            t_294[k] = -ab_x[k] * gh_294[k]
                       + gi_392[k];
        }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_x, gh_295, gh_296, gh_297, \
                         gh_298, gh_299, gi_393, gi_394, gi_395, gi_396, \
                         gi_397 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_295[k] = -ab_x[k] * gh_295[k]
                       + gi_393[k];

            t_296[k] = -ab_x[k] * gh_296[k]
                       + gi_394[k];

            t_297[k] = -ab_x[k] * gh_297[k]
                       + gi_395[k];

            t_298[k] = -ab_x[k] * gh_298[k]
                       + gi_396[k];

            t_299[k] = -ab_x[k] * gh_299[k]
                       + gi_397[k];
        }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, gh_300, gh_301, gh_302, \
                         gh_303, gh_304, gi_398, gi_399, gi_400, gi_401, \
                         gi_402 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_300[k] = -ab_x[k] * gh_300[k]
                       + gi_398[k];

            t_301[k] = -ab_x[k] * gh_301[k]
                       + gi_399[k];

            t_302[k] = -ab_x[k] * gh_302[k]
                       + gi_400[k];

            t_303[k] = -ab_x[k] * gh_303[k]
                       + gi_401[k];

            t_304[k] = -ab_x[k] * gh_304[k]
                       + gi_402[k];
        }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, gh_305, gh_306, gh_307, \
                         gh_308, gh_309, gi_403, gi_404, gi_405, gi_406, \
                         gi_407 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_305[k] = -ab_x[k] * gh_305[k]
                       + gi_403[k];

            t_306[k] = -ab_x[k] * gh_306[k]
                       + gi_404[k];

            t_307[k] = -ab_x[k] * gh_307[k]
                       + gi_405[k];

            t_308[k] = -ab_x[k] * gh_308[k]
                       + gi_406[k];

            t_309[k] = -ab_x[k] * gh_309[k]
                       + gi_407[k];
        }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_x, gh_310, gh_311, gh_312, \
                         gh_313, gh_314, gi_408, gi_409, gi_410, gi_411, \
                         gi_412 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_310[k] = -ab_x[k] * gh_310[k]
                       + gi_408[k];

            t_311[k] = -ab_x[k] * gh_311[k]
                       + gi_409[k];

            t_312[k] = -ab_x[k] * gh_312[k]
                       + gi_410[k];

            t_313[k] = -ab_x[k] * gh_313[k]
                       + gi_411[k];

            t_314[k] = -ab_x[k] * gh_314[k]
                       + gi_412[k];
        }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_y, gh_210, gh_211, gh_212, \
                         gh_213, gh_214, gi_281, gi_283, gi_284, gi_286, \
                         gi_287 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_315[k] = -ab_y[k] * gh_210[k]
                       + gi_281[k];

            t_316[k] = -ab_y[k] * gh_211[k]
                       + gi_283[k];

            t_317[k] = -ab_y[k] * gh_212[k]
                       + gi_284[k];

            t_318[k] = -ab_y[k] * gh_213[k]
                       + gi_286[k];

            t_319[k] = -ab_y[k] * gh_214[k]
                       + gi_287[k];
        }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_y, gh_215, gh_216, gh_217, \
                         gh_218, gh_219, gi_288, gi_290, gi_291, gi_292, \
                         gi_293 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_320[k] = -ab_y[k] * gh_215[k]
                       + gi_288[k];

            t_321[k] = -ab_y[k] * gh_216[k]
                       + gi_290[k];

            t_322[k] = -ab_y[k] * gh_217[k]
                       + gi_291[k];

            t_323[k] = -ab_y[k] * gh_218[k]
                       + gi_292[k];

            t_324[k] = -ab_y[k] * gh_219[k]
                       + gi_293[k];
        }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_y, gh_220, gh_221, gh_222, \
                         gh_223, gh_224, gi_295, gi_296, gi_297, gi_298, \
                         gi_299 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_325[k] = -ab_y[k] * gh_220[k]
                       + gi_295[k];

            t_326[k] = -ab_y[k] * gh_221[k]
                       + gi_296[k];

            t_327[k] = -ab_y[k] * gh_222[k]
                       + gi_297[k];

            t_328[k] = -ab_y[k] * gh_223[k]
                       + gi_298[k];

            t_329[k] = -ab_y[k] * gh_224[k]
                       + gi_299[k];
        }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_y, gh_225, gh_226, gh_227, \
                         gh_228, gh_229, gi_301, gi_302, gi_303, gi_304, \
                         gi_305 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_330[k] = -ab_y[k] * gh_225[k]
                       + gi_301[k];

            t_331[k] = -ab_y[k] * gh_226[k]
                       + gi_302[k];

            t_332[k] = -ab_y[k] * gh_227[k]
                       + gi_303[k];

            t_333[k] = -ab_y[k] * gh_228[k]
                       + gi_304[k];

            t_334[k] = -ab_y[k] * gh_229[k]
                       + gi_305[k];
        }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_y, gh_230, gh_231, gh_232, \
                         gh_233, gh_234, gi_306, gi_309, gi_311, gi_312, \
                         gi_314 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_335[k] = -ab_y[k] * gh_230[k]
                       + gi_306[k];

            t_336[k] = -ab_y[k] * gh_231[k]
                       + gi_309[k];

            t_337[k] = -ab_y[k] * gh_232[k]
                       + gi_311[k];

            t_338[k] = -ab_y[k] * gh_233[k]
                       + gi_312[k];

            t_339[k] = -ab_y[k] * gh_234[k]
                       + gi_314[k];
        }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_y, gh_235, gh_236, gh_237, \
                         gh_238, gh_239, gi_315, gi_316, gi_318, gi_319, \
                         gi_320 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_340[k] = -ab_y[k] * gh_235[k]
                       + gi_315[k];

            t_341[k] = -ab_y[k] * gh_236[k]
                       + gi_316[k];

            t_342[k] = -ab_y[k] * gh_237[k]
                       + gi_318[k];

            t_343[k] = -ab_y[k] * gh_238[k]
                       + gi_319[k];

            t_344[k] = -ab_y[k] * gh_239[k]
                       + gi_320[k];
        }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_y, gh_240, gh_241, gh_242, \
                         gh_243, gh_244, gi_321, gi_323, gi_324, gi_325, \
                         gi_326 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_345[k] = -ab_y[k] * gh_240[k]
                       + gi_321[k];

            t_346[k] = -ab_y[k] * gh_241[k]
                       + gi_323[k];

            t_347[k] = -ab_y[k] * gh_242[k]
                       + gi_324[k];

            t_348[k] = -ab_y[k] * gh_243[k]
                       + gi_325[k];

            t_349[k] = -ab_y[k] * gh_244[k]
                       + gi_326[k];
        }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_y, gh_245, gh_246, gh_247, \
                         gh_248, gh_249, gi_327, gi_329, gi_330, gi_331, \
                         gi_332 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_350[k] = -ab_y[k] * gh_245[k]
                       + gi_327[k];

            t_351[k] = -ab_y[k] * gh_246[k]
                       + gi_329[k];

            t_352[k] = -ab_y[k] * gh_247[k]
                       + gi_330[k];

            t_353[k] = -ab_y[k] * gh_248[k]
                       + gi_331[k];

            t_354[k] = -ab_y[k] * gh_249[k]
                       + gi_332[k];
        }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_y, gh_250, gh_251, gh_252, \
                         gh_253, gh_254, gi_333, gi_334, gi_337, gi_339, \
                         gi_340 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_355[k] = -ab_y[k] * gh_250[k]
                       + gi_333[k];

            t_356[k] = -ab_y[k] * gh_251[k]
                       + gi_334[k];

            t_357[k] = -ab_y[k] * gh_252[k]
                       + gi_337[k];

            t_358[k] = -ab_y[k] * gh_253[k]
                       + gi_339[k];

            t_359[k] = -ab_y[k] * gh_254[k]
                       + gi_340[k];
        }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_y, gh_255, gh_256, gh_257, \
                         gh_258, gh_259, gi_342, gi_343, gi_344, gi_346, \
                         gi_347 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_360[k] = -ab_y[k] * gh_255[k]
                       + gi_342[k];

            t_361[k] = -ab_y[k] * gh_256[k]
                       + gi_343[k];

            t_362[k] = -ab_y[k] * gh_257[k]
                       + gi_344[k];

            t_363[k] = -ab_y[k] * gh_258[k]
                       + gi_346[k];

            t_364[k] = -ab_y[k] * gh_259[k]
                       + gi_347[k];
        }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_y, gh_260, gh_261, gh_262, \
                         gh_263, gh_264, gi_348, gi_349, gi_351, gi_352, \
                         gi_353 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_365[k] = -ab_y[k] * gh_260[k]
                       + gi_348[k];

            t_366[k] = -ab_y[k] * gh_261[k]
                       + gi_349[k];

            t_367[k] = -ab_y[k] * gh_262[k]
                       + gi_351[k];

            t_368[k] = -ab_y[k] * gh_263[k]
                       + gi_352[k];

            t_369[k] = -ab_y[k] * gh_264[k]
                       + gi_353[k];
        }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_y, gh_265, gh_266, gh_267, \
                         gh_268, gh_269, gi_354, gi_355, gi_357, gi_358, \
                         gi_359 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_370[k] = -ab_y[k] * gh_265[k]
                       + gi_354[k];

            t_371[k] = -ab_y[k] * gh_266[k]
                       + gi_355[k];

            t_372[k] = -ab_y[k] * gh_267[k]
                       + gi_357[k];

            t_373[k] = -ab_y[k] * gh_268[k]
                       + gi_358[k];

            t_374[k] = -ab_y[k] * gh_269[k]
                       + gi_359[k];
        }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_y, gh_270, gh_271, gh_272, \
                         gh_273, gh_274, gi_360, gi_361, gi_362, gi_365, \
                         gi_367 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_375[k] = -ab_y[k] * gh_270[k]
                       + gi_360[k];

            t_376[k] = -ab_y[k] * gh_271[k]
                       + gi_361[k];

            t_377[k] = -ab_y[k] * gh_272[k]
                       + gi_362[k];

            t_378[k] = -ab_y[k] * gh_273[k]
                       + gi_365[k];

            t_379[k] = -ab_y[k] * gh_274[k]
                       + gi_367[k];
        }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_y, gh_275, gh_276, gh_277, \
                         gh_278, gh_279, gi_368, gi_370, gi_371, gi_372, \
                         gi_374 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_380[k] = -ab_y[k] * gh_275[k]
                       + gi_368[k];

            t_381[k] = -ab_y[k] * gh_276[k]
                       + gi_370[k];

            t_382[k] = -ab_y[k] * gh_277[k]
                       + gi_371[k];

            t_383[k] = -ab_y[k] * gh_278[k]
                       + gi_372[k];

            t_384[k] = -ab_y[k] * gh_279[k]
                       + gi_374[k];
        }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_y, gh_280, gh_281, gh_282, \
                         gh_283, gh_284, gi_375, gi_376, gi_377, gi_379, \
                         gi_380 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_385[k] = -ab_y[k] * gh_280[k]
                       + gi_375[k];

            t_386[k] = -ab_y[k] * gh_281[k]
                       + gi_376[k];

            t_387[k] = -ab_y[k] * gh_282[k]
                       + gi_377[k];

            t_388[k] = -ab_y[k] * gh_283[k]
                       + gi_379[k];

            t_389[k] = -ab_y[k] * gh_284[k]
                       + gi_380[k];
        }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_y, gh_285, gh_286, gh_287, \
                         gh_288, gh_289, gi_381, gi_382, gi_383, gi_385, \
                         gi_386 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_390[k] = -ab_y[k] * gh_285[k]
                       + gi_381[k];

            t_391[k] = -ab_y[k] * gh_286[k]
                       + gi_382[k];

            t_392[k] = -ab_y[k] * gh_287[k]
                       + gi_383[k];

            t_393[k] = -ab_y[k] * gh_288[k]
                       + gi_385[k];

            t_394[k] = -ab_y[k] * gh_289[k]
                       + gi_386[k];
        }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_y, gh_290, gh_291, gh_292, \
                         gh_293, gh_294, gi_387, gi_388, gi_389, gi_390, \
                         gi_393 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_395[k] = -ab_y[k] * gh_290[k]
                       + gi_387[k];

            t_396[k] = -ab_y[k] * gh_291[k]
                       + gi_388[k];

            t_397[k] = -ab_y[k] * gh_292[k]
                       + gi_389[k];

            t_398[k] = -ab_y[k] * gh_293[k]
                       + gi_390[k];

            t_399[k] = -ab_y[k] * gh_294[k]
                       + gi_393[k];
        }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_y, gh_295, gh_296, gh_297, \
                         gh_298, gh_299, gi_395, gi_396, gi_398, gi_399, \
                         gi_400 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_400[k] = -ab_y[k] * gh_295[k]
                       + gi_395[k];

            t_401[k] = -ab_y[k] * gh_296[k]
                       + gi_396[k];

            t_402[k] = -ab_y[k] * gh_297[k]
                       + gi_398[k];

            t_403[k] = -ab_y[k] * gh_298[k]
                       + gi_399[k];

            t_404[k] = -ab_y[k] * gh_299[k]
                       + gi_400[k];
        }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_y, gh_300, gh_301, gh_302, \
                         gh_303, gh_304, gi_402, gi_403, gi_404, gi_405, \
                         gi_407 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_405[k] = -ab_y[k] * gh_300[k]
                       + gi_402[k];

            t_406[k] = -ab_y[k] * gh_301[k]
                       + gi_403[k];

            t_407[k] = -ab_y[k] * gh_302[k]
                       + gi_404[k];

            t_408[k] = -ab_y[k] * gh_303[k]
                       + gi_405[k];

            t_409[k] = -ab_y[k] * gh_304[k]
                       + gi_407[k];
        }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_y, gh_305, gh_306, gh_307, \
                         gh_308, gh_309, gi_408, gi_409, gi_410, gi_411, \
                         gi_413 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_410[k] = -ab_y[k] * gh_305[k]
                       + gi_408[k];

            t_411[k] = -ab_y[k] * gh_306[k]
                       + gi_409[k];

            t_412[k] = -ab_y[k] * gh_307[k]
                       + gi_410[k];

            t_413[k] = -ab_y[k] * gh_308[k]
                       + gi_411[k];

            t_414[k] = -ab_y[k] * gh_309[k]
                       + gi_413[k];
        }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_y, gh_310, gh_311, gh_312, \
                         gh_313, gh_314, gi_414, gi_415, gi_416, gi_417, \
                         gi_418 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_415[k] = -ab_y[k] * gh_310[k]
                       + gi_414[k];

            t_416[k] = -ab_y[k] * gh_311[k]
                       + gi_415[k];

            t_417[k] = -ab_y[k] * gh_312[k]
                       + gi_416[k];

            t_418[k] = -ab_y[k] * gh_313[k]
                       + gi_417[k];

            t_419[k] = -ab_y[k] * gh_314[k]
                       + gi_418[k];
        }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_z, gh_294, gh_295, gh_296, \
                         gh_297, gh_298, gi_394, gi_396, gi_397, gi_399, \
                         gi_400 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_420[k] = -ab_z[k] * gh_294[k]
                       + gi_394[k];

            t_421[k] = -ab_z[k] * gh_295[k]
                       + gi_396[k];

            t_422[k] = -ab_z[k] * gh_296[k]
                       + gi_397[k];

            t_423[k] = -ab_z[k] * gh_297[k]
                       + gi_399[k];

            t_424[k] = -ab_z[k] * gh_298[k]
                       + gi_400[k];
        }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_z, gh_299, gh_300, gh_301, \
                         gh_302, gh_303, gi_401, gi_403, gi_404, gi_405, \
                         gi_406 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_425[k] = -ab_z[k] * gh_299[k]
                       + gi_401[k];

            t_426[k] = -ab_z[k] * gh_300[k]
                       + gi_403[k];

            t_427[k] = -ab_z[k] * gh_301[k]
                       + gi_404[k];

            t_428[k] = -ab_z[k] * gh_302[k]
                       + gi_405[k];

            t_429[k] = -ab_z[k] * gh_303[k]
                       + gi_406[k];
        }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_z, gh_304, gh_305, gh_306, \
                         gh_307, gh_308, gi_408, gi_409, gi_410, gi_411, \
                         gi_412 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_430[k] = -ab_z[k] * gh_304[k]
                       + gi_408[k];

            t_431[k] = -ab_z[k] * gh_305[k]
                       + gi_409[k];

            t_432[k] = -ab_z[k] * gh_306[k]
                       + gi_410[k];

            t_433[k] = -ab_z[k] * gh_307[k]
                       + gi_411[k];

            t_434[k] = -ab_z[k] * gh_308[k]
                       + gi_412[k];
        }
    }
}

static auto
compute_hrr_hh_piece3(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
                      const size_t gh, const size_t gi, const size_t ncomps,
                      const size_t nmax) -> void
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

        const auto *ab_z = coordinates.data(8);

        const auto *gh_309 = buffer.data(gh + 309 * ncomps + c);
        const auto *gh_310 = buffer.data(gh + 310 * ncomps + c);
        const auto *gh_311 = buffer.data(gh + 311 * ncomps + c);
        const auto *gh_312 = buffer.data(gh + 312 * ncomps + c);
        const auto *gh_313 = buffer.data(gh + 313 * ncomps + c);
        const auto *gh_314 = buffer.data(gh + 314 * ncomps + c);

        const auto *gi_414 = buffer.data(gi + 414 * ncomps + c);
        const auto *gi_415 = buffer.data(gi + 415 * ncomps + c);
        const auto *gi_416 = buffer.data(gi + 416 * ncomps + c);
        const auto *gi_417 = buffer.data(gi + 417 * ncomps + c);
        const auto *gi_418 = buffer.data(gi + 418 * ncomps + c);
        const auto *gi_419 = buffer.data(gi + 419 * ncomps + c);

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_z, gh_309, gh_310, gh_311, \
                         gh_312, gh_313, gi_414, gi_415, gi_416, gi_417, \
                         gi_418 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_435[k] = -ab_z[k] * gh_309[k]
                       + gi_414[k];

            t_436[k] = -ab_z[k] * gh_310[k]
                       + gi_415[k];

            t_437[k] = -ab_z[k] * gh_311[k]
                       + gi_416[k];

            t_438[k] = -ab_z[k] * gh_312[k]
                       + gi_417[k];

            t_439[k] = -ab_z[k] * gh_313[k]
                       + gi_418[k];
        }

#pragma omp simd aligned(t_440, ab_z, gh_314, gi_419 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_440[k] = -ab_z[k] * gh_314[k]
                       + gi_419[k];
        }
    }
}

auto
compute_hrr_hh(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t gh, const size_t gi, const size_t ncomps, const size_t nmax) -> void
{
    compute_hrr_hh_piece0(buffer, coordinates, target, gh, gi, ncomps, nmax);

    compute_hrr_hh_piece1(buffer, coordinates, target, gh, gi, ncomps, nmax);

    compute_hrr_hh_piece2(buffer, coordinates, target, gh, gi, ncomps, nmax);

    compute_hrr_hh_piece3(buffer, coordinates, target, gh, gi, ncomps, nmax);
}

}  // namespace simdtrf
